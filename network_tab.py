import glob
import itertools

import pandas as pd
import numpy as np

from datetime import datetime, timedelta
from functools import partial

from bokeh.io import show, curdoc
from bokeh.plotting import figure

from bokeh.models import CategoricalColorMapper, HoverTool, ColumnDataSource, Panel
from bokeh.models.widgets import MultiSelect,RadioButtonGroup, Button, CheckboxButtonGroup, CheckboxGroup, Slider, RangeSlider, Tabs

from bokeh.layouts import gridplot, column, row , widgetbox
from bokeh.palettes import Category20_16

from obspy.signal import PPSD
from obspy.signal.spectral_estimation import get_nlnm, get_nhnm 

import fnmatch

maxnstation = 29*3

def read(pastdays_select=None,
         bufferdays_select=None, 
         #whitelists = ['CH.*..HG*','XE.*.HN*'],
         whitelists = ['SV.*.00.HN*','XE.*.00.HN*'],
         slinktooldir = '/home/sysop/slinktool/latencies/',
         psddir = '/home/sysop/msnoise/PSD/NPZ/',
         src=None,
         flights=None):

    dfs = []
    ppsds = []
    buffersizes = {}

    files = []
    for w in whitelists:
        parts = (slinktooldir,w)
        filenames = '%s/*/*/*/*.D/%s*.????.???'%parts       
        files += glob.glob(filenames)     
    yeardays = list(set(['.'.join(f.split('.')[-2:]) for f in files]))
    mseedids = list(set(['.'.join(f.split('/')[-1].split('.')[:4]) for f in files]))
    #print(files)#yeardays)
    #print(mseedids)

    now = datetime.now()
    for n in range(int(pastdays_select.value[0]),int(1+pastdays_select.value[1])):
        if len(list(buffersizes.keys()))>=maxnstation:
            break
        d = timedelta(days = n)
        year = (now-d).year
        day = (now-d).timetuple().tm_yday
        if '%s.%s'%(year,day) not in yeardays: 
            continue
        
        parts = (slinktooldir,year,year,day)
        dayfilenames = '%s/%d/*/*/*.D/*.%d.%d'%parts      
        dayfilenames = dayfilenames.replace('//','/')
        #print(files)
        #print(dayfilenames)
        dayfilenames = [f for f in files if fnmatch.fnmatch(f, dayfilenames)]
        if len(dayfilenames)==0:
            continue
        #print(dayfilenames)

        for whitelist in whitelists:            
            if len(list(buffersizes.keys()))>=maxnstation:
                break
            # Read in data from slinktool
            parts = (slinktooldir,year,whitelist,year,day)
            filenames = '%s/%d/*/*/*.D/%s.%d.%d'%parts       
            filenames = filenames.replace('//','/')
            #print(filenames)
            filenames = [f for f in dayfilenames if fnmatch.fnmatch(f, filenames)]
            #print(filenames)

            for filename in filenames[::-1]:
                if len(list(buffersizes.keys()))>=maxnstation:
                    break
                mseedid = '.'.join(filename.split('/')[-1].split('.')[:4])
                if mseedid not in buffersizes:
                    buffersizes[mseedid]=0
                if buffersizes[mseedid]+1 > np.ceil(bufferdays_select.value):
                    continue
                print('loading %s (%s lines)'%(filename,int(bufferdays_select.value*24*60*60)))
                buffersizes[mseedid] += 1
                df = pd.read_csv(filename, 
                                 nrows=int(bufferdays_select.value*24*60*60),
                                 sep=' ',
                                 #XE_CH01_00_HNZ, 100 samples, 100 Hz, 2020,141,19:01:44.040000 (latency ~1.3802 sec)#
                                 names= ['name', 'reclen', 'unitlen', 'freq', 'unitfreq', 'time', '(latency', 'arr_delay', 'unitdelay'],
                                 index_col=False)[['arr_delay', 'time', 'name']]
                
                # Trim
                df['time'] = pd.to_datetime(df['time'], format='%Y,%j,%H:%M:%S.%f')
                earliest = max(df['time'])-timedelta(days=bufferdays_select.value)
                df = df.loc[df['time']>=earliest]
                
                # clean
                df['name'] = df['name'].map(lambda x: x.rstrip(','))
                df['name'] = df['name'].str.replace("_",".")

                # Convert
                df['Delays'] = df['arr_delay'].map(lambda x: x.lstrip('~')).astype(float)
                
                # Compute reception time and latency
                df['arr_time'] = pd.to_datetime(df['time']) + pd.to_timedelta(df['Delays'], unit='s')
                df['Latencies'] = df['arr_time'].diff() / np.timedelta64(1, 's') 
                
                # Concatenate data
                dfs.append(df[['name','arr_time','Latencies', 'Delays']])
                
                # get psd
                print(filename.replace(slinktooldir,psddir)+'.npz')
                filename = glob.glob(filename.replace(slinktooldir,psddir)+'.npz')
                print(filename)
                if len(filename):
                    print('loading',filename[0])
                    ppsd = PPSD.load_npz(filename[0])
                    rows = []
                    times = []
                    indexes = np.argsort(ppsd.times_processed)
                    for indexpsd in indexes:
                        time=ppsd.times_processed[indexpsd]
                        if time<min(df['time']):
                            continue
                        if time>max(df['time']):
                            break
                        times += [time.datetime]
                        rows += [[mseedid, ppsd.period_bin_centers, ppsd.psd_values[indexpsd]]]
                    if len(times):
                        ppsds.append(pd.DataFrame(np.array(rows),
                                              index=np.array(times),
                                              columns=['name', 'PSD_periods', 'PSD']))
                #print('OK')
    # Concatenate all data into one DataFrame
    flights = pd.concat(dfs, ignore_index=True)
    ppsds = pd.concat(ppsds, ignore_index=True)

    # Available carrier list
    available_carriers = list(flights['name'].unique())

    # Sort the list in-place (alphabetical order)
    available_carriers.sort()

    return [flights, ppsds], available_carriers
        
# Dataset based on selected carriers, selected start and end of delays,
# and selected width of bin
def make_dataset(flights, 
                 carrier_list, 
                 range_start, 
                 range_end, 
                 bin_width, 
                 output='Latencies',
                 plottype='Cum'):

    data = {'proportion':[], 
               'center':[],  
               'name':[],
               'fullname':[],
               'color':[]}
    names = []
    #range_extent = int((range_end - range_start) / bin_width)
    range_extent = np.logspace(np.log10(max([0.01,range_start])), np.log10(range_end), bin_width, endpoint=True)
    if range_start<-0.009:
        range_extent = np.append(np.sort(-1*np.logspace(np.log10(0.009), # 0.01 creates a point at 0, do not do that
                                                        np.log10(range_start*-1), 
                                                        bin_width, 
                                                        endpoint=True)),
                                 range_extent)
    if output == 'PSD':
        
        names+=['N(H,L)NM']
        data['fullname'] += ['NLNM']
        data['fullname'] += ['NHNM']
        
        for d in [get_nlnm(),get_nhnm()]:
            data['proportion'] += [d[1][(d[0] >= range_start) & (d[0] <= range_end)]]
            data['center'] += [d[0][(d[0] >= range_start) & (d[0] <= range_end)]]
            data['name']+=['N(H,L)M']
            data['color'] += ['black']
        
        for i, carrier_name in enumerate(carrier_list):
            subset = flights[1][flights[1]['name'] == carrier_name ]
            d=[]
            for i, row in subset.iterrows():
                psd = row['PSD'][(row['PSD_periods'] >= range_start) & (row['PSD_periods'] <= range_end)]
                psd_periods = row['PSD_periods'][(row['PSD_periods'] >= range_start) & (row['PSD_periods'] <= range_end)]
                d+=[list(psd)]
            if len(d):
                data['proportion'] += [np.percentile(d, 83, axis=0)] 
                data['center'] += [psd_periods]            
                # Assign the carrier for labels
                data['fullname'] += [carrier_name]
                seedid = carrier_name.split('.')
                seedid[1] = '*'
                seedid[2] = '*'
                data['name'] += ['.'.join(seedid)]
                # Color each carrier differently
                if data['name'][-1] not in names:
                    names+=[data['name'][-1]]
                data['color'] += [Category20_16[names.index(data['name'][-1])]]
            
        return ColumnDataSource(data = data)

    
    data = {'proportion':[], 
            'steproportion':[],
               'left':[], 
               'right':[],
               'center':[],  
               'leftright':[], 
               'f_interval':[],
               'step_interval':[],
               'f_proportion':[], 
               'step_proportion':[], 
               'name':[],
               'fullname':[],
               'color':[]}
    
    # Iterate through all the carriers
    for i, carrier_name in enumerate(carrier_list):

        # Subset to the carrier
        subset = flights[0][flights[0]['name'] == carrier_name ]
        subset = subset[output]
        subset = subset[subset >= range_start]
        subset = subset[subset <= range_end]

        # Create a histogram with 5 minute bins
        arr_hist, edges = np.histogram(subset, 
                                       bins = range_extent, 
                                       range = [range_start, range_end])

        # Divide the counts by the total to get a proportion
        if not np.nansum(arr_hist)> 0 :
            continue
        data['proportion'] += [arr_hist / np.nansum(arr_hist) * np.sign(edges[:-1])] 
        data['left'] += [edges[:-1]] 
        data['right'] += [edges[1:]] 
        data['center'] += [np.abs((edges[:-1]+edges[1:])/2) ]

        if 'Cum' in plottype:
            data['proportion'][-1] = np.nancumsum(data['proportion'][-1])

        # Format the proportion 
        data['f_proportion'] += [['%0.5f' % proportion for proportion in data['proportion'][-1]]]

        # Format the interval
        data['f_interval'] += [['%.3f to %.3f sec' % (left, right) for left, right in zip(data['left'][-1], data['right'][-1])]]
        
        # Step data
        data['steproportion'] += [list(itertools.chain.from_iterable(zip(data['proportion'][-1],data['proportion'][-1])))]
        data['leftright'] += [list(itertools.chain.from_iterable(zip(data['left'][-1],data['right'][-1])))]
        data['step_proportion'] += [list(itertools.chain.from_iterable(zip(data['f_proportion'][-1],data['f_proportion'][-1])))]
        data['step_interval'] += [list(itertools.chain.from_iterable(zip(data['f_interval'][-1],data['f_interval'][-1])))]
        

        # Assign the carrier for labels
        data['fullname'] += [carrier_name]
        seedid = carrier_name.split('.')
        seedid[1] = '*'
        seedid[2] = '*'
        data['name'] += ['.'.join(seedid)]

        # Color each carrier differently
        if data['name'][-1] not in names:
            names+=[data['name'][-1]]
        data['color'] += [Category20_16[names.index(data['name'][-1])]]
        

    # Overall dataframe
    #by_carrier = by_carrier.sort_values(['name', 'left'])
    
    return ColumnDataSource(data = data)

# Styling for a plot
def style(p):
    # Title
    #p.title.align = 'center'
    #p.title.text_font_size = '20pt'
    ##p.title.text_font = 'serif'

    # Axis titles
    #p.xaxis.axis_label_text_font_size = '14pt'
    p.xaxis.axis_label_text_font_style = 'bold'
    #p.yaxis.axis_label_text_font_size = '14pt'
    p.yaxis.axis_label_text_font_style = 'bold'

    # Tick labels
    #p.xaxis.major_label_text_font_size = '12pt'
    #p.yaxis.major_label_text_font_size = '12pt'

    return p

# Function to make the plot
def make_plot(src, output='Latencies'):
    # Blank plot with correct labels
    if output is 'PSD':
        p = figure(sizing_mode="stretch_both",#plot_width = 700, plot_height = 300, 
                   x_axis_type="log",
                  title = 'High noise %s by Channel'%output,
                  x_axis_label = 'Period (sec)', y_axis_label = '%s (dB)'%output)
        
    else:
        p = figure(sizing_mode="stretch_both",#plot_width = 700, plot_height = 300, 
                   x_axis_type="log",
                  title = 'Histogram of Packet %s by Channel'%output,
                  x_axis_label = '%s (sec)'%output, y_axis_label = 'Proportion')

    if True:
        p.multi_line('center',#leftright',
                     'proportion',#steproportion', 
                     source = src, 
                     line_color = 'color', 
                     alpha=0.5,
                     legend_group = 'name',
                     hover_line_color='black',
                     hover_line_alpha=1.0,
                     hover_alpha = 1, 
                     hover_line_width = 3,
                     #line_width=0.5
                     )
        p.legend.location = 'bottom_left'
        p.legend.orientation = "horizontal"
    else:
        # Quad glyphs to create a histogram
        p.quad(source = src, 
               bottom = 0, 
               top = 'proportion', 
               left = 'left', 
               right = 'right',
               line_color = 'color', 
               fill_alpha = 0, 
               hover_fill_color = 'pink', 
               hover_fill_alpha = 0, 
               #legend_group = 'name',
               #color = 'color', 
               line_width=0.5)
    
    # Hover tool with vline mode
    hover = HoverTool(tooltips=[('SeedId', '@fullname'), 
                                #('%s'%output, '$step_interval'),
                                #('Proportion', '$step_proportion')
                                ],
                      #mode='vline'
                      )

    p.add_tools(hover)

    p.legend.click_policy = 'hide'

    # Styling
    p = style(p)

    return p

def needsreadupdate(attr,old,new,
                    refresh_button=None):
    refresh_button.disabled = False
    refresh_button.button_type = "warning"
    refresh_button.label = "Refresh time window and buffer"

def readupdate(new,#attr, old, new,
               src=None,
               selections=None,
               flights=None,
               type_switch=None,
               range_select=None,
               binwidth_select=None,
               pastdays_select=None,
               bufferdays_select=None,
               refresh_button=None):

    refresh_button.update(disabled = True,
            button_type = "danger",
            label = 'Reloading, please wait...')

    flights, available_carriers = read(pastdays_select,bufferdays_select)
    update(None,None,None,#attr, old, new, 
           output='Delays',
           flights=flights,
           src=src,
            type_switch=type_switch,
           selections=selections,
           range_select=range_select,
           binwidth_select=binwidth_select) 
    update(None,None,None,#attr, old, new,
           output='Latencies',
           flights=flights,
            type_switch=type_switch,
           src=src,
           selections=selections,
           range_select=range_select,
           binwidth_select=binwidth_select) 
    
    refresh_button.disabled = True
    refresh_button.button_type = "success"
    refresh_button.label = 'Time window and buffer are up to date'



# Update the plot based on selections
def update(attr, old, new, 
           output='Delays',
           type_switch=None,
           src=None,
           selections=None,
           flights=None,
           range_select=None,
           binwidth_select=None):
    
    to_plot = [s.labels[i] for s in selections for i in s.active]
    plottype=type_switch.labels[type_switch.active]
    
    new_src = {output: make_dataset(flights,
                                    to_plot,
                                    plottype = plottype,
                                    range_start = range_select.value[0],
                                    range_end = range_select.value[1],
                                    bin_width = binwidth_select.value,
                                    output=output)}

    src[output].data.update(new_src[output].data)


# Initiate the tab
def make_tab():
    # Slider to select width of bin
    pastdays_select = RangeSlider(start = 0, end = 999, 
                                  value = (0,999), step = 1, 
                                  title = 'Past Days',
                                  sizing_mode="stretch_both")
    # Slider to select buffer size
    bufferdays_select = Slider(start = .01, end = 9, 
                             value = 0.01, step = .1, 
                             title = 'Buffer Size (days)',
                             sizing_mode="stretch_both")
    # Re-read
    refresh_button = Button(label="Time window and buffer are up to date", 
                            button_type="success",
                            sizing_mode="stretch_both")
    refresh_button.disabled = True

    # read data
    flights, available_carriers = read(pastdays_select,bufferdays_select)

    # CheckboxGroup to select carrier to display
    locationcodes = np.unique(['.'.join(seedid.split('.')[:3]) for seedid in available_carriers])
    selections = []
    for locationcode in locationcodes[:maxnstation]:
        matching = [s for s in available_carriers if locationcode in s]
        active = [i for i, m in enumerate(matching)]# if "Z" == m[-1]]
        selections += [CheckboxButtonGroup(labels=matching, 
                                          active = active,
                                          sizing_mode="stretch_both")]
        #selections += [MultiSelect(#title="Option:", 
        #                           value=matching,
        #                           active=active)

        
    # Find the initially selected carrieres
    initial_carriers = [s.labels[i] for s in selections for i in s.active]

    # Slider to select width of bin
    binwidth_select = Slider(start = 16, end = 160, 
                         step = 16, value = 80,
                         title = 'Bin number',
                         sizing_mode="stretch_both")
    
    # RangeSlider control to select start and end of plotted delays
    range_select = RangeSlider(start = -1, end = 999, value = (-.2, 99),
                               step = .1, title = 'Range (sec)',
                               sizing_mode="stretch_both")
 
    # Switch from lines to hists
    type_switch = RadioButtonGroup(labels=["Histogram", "Cumulated dist."], active=0,sizing_mode="stretch_both")

    # Find the initially selected carrieres
    plottype = type_switch.labels[type_switch.active]

    src={}
    for output in ['Latencies','Delays', 'PSD']:
        src[output] = make_dataset(flights,
                                   initial_carriers,
                          range_start = range_select.value[0],
                          range_end = range_select.value[1],
                          bin_width = binwidth_select.value,
                          output=output,
                          plottype=plottype)

    callback = partial(update,
            output='Delays',
             type_switch=type_switch,
            flights=flights,
            src=src,
            selections=selections,
            range_select=range_select,
            binwidth_select=binwidth_select)

    callbacklat = partial(update,
            output='Latencies',
            type_switch=type_switch,
            flights=flights,
            src=src,
            range_select=range_select,
            selections=selections,
            binwidth_select=binwidth_select)
    
    callbackpsd = partial(update,
            output='PSD',
            type_switch=type_switch,
            flights=flights,
            src=src,
            range_select=range_select,
            selections=selections,
            binwidth_select=binwidth_select)
    
    callbackneedsread = partial(needsreadupdate,
                                refresh_button=refresh_button)
    
    callbackread = partial(readupdate,
                           src=src,
                           selections=selections,
                           range_select=range_select,
                           type_switch=type_switch,
                           binwidth_select=binwidth_select,
                           pastdays_select=pastdays_select,
                           bufferdays_select=bufferdays_select,
                           refresh_button=refresh_button)

    [s.on_change('active', callback, callbacklat, callbackpsd) for s in selections]
    type_switch.on_change('active', callback, callbacklat, callbackpsd)
    binwidth_select.on_change('value',  callback, callbacklat, callbackpsd)
    range_select.on_change('value',  callback, callbacklat, callbackpsd)
    pastdays_select.on_change('value', callbackneedsread)
    bufferdays_select.on_change('value', callbackneedsread)
    refresh_button.on_click(callbackread)

    p={}
    for output in ['PSD','Latencies','Delays']:
        p[output] = make_plot(src[output],
                              output=output)

    # Create a row layout
    graphs = [p[k] for k in p]
    controls = [type_switch,
                binwidth_select,
                range_select,
                refresh_button, 
                pastdays_select,
                bufferdays_select,
                *selections[:maxnstation]
                ]
    graphslayout = column(children=graphs,
                          sizing_mode="stretch_both"
                          )
    controlslayout = column(children=controls,
                               sizing_mode='fixed',#stretch_width',
                               width=400,
                               )
    layout = row(children=[graphslayout,
                           controlslayout],
                 sizing_mode='stretch_both'
                 )

    # Make a tab with the layout 
    return Panel(child = layout, 
                 title = 'Channels')


#tabs = Tabs(tabs=[tab])

## Add it to the current document (displays plot)
#curdoc().add_root(tabs)
