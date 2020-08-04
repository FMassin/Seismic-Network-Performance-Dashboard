from bokeh.io import curdoc
from bokeh.models.widgets import Tabs

from network_tab import make_tab as network_tab

tabs = Tabs(tabs=[network_tab()])

# Add it to the current document (displays plot)
curdoc().add_root(tabs)
curdoc().title = "Seismic Network Performances"

#from bokeh.plotting import output_file,save
#output_file("/home/sysop/dashboard.html")
#save(tabs)

