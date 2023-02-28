#%%
import datetime
from argparse import Namespace
import copy

import main

#%%
default_args = Namespace(
    download=False,
    plot=False,
    table=False,
    map=False,
    output_dir="latest",
)

#%%
# Download data
# download_args = copy.deepcopy(default_args)
# download_args.download = True
# main.main(download_args)

#%%
# Update plot until previous calendar week
plot_args = copy.deepcopy(default_args)
plot_args.plot = True
plot_args.start_date = (datetime.date.today() - datetime.timedelta(days=12*7)).strftime('%Y-%m-%d')
plot_args.end_date = (datetime.date.today() - datetime.timedelta(days=7)).strftime('%Y-%m-%d')
main.main(plot_args)

#%%
# Update map for previous 4 weeks
map_args = copy.deepcopy(default_args)
map_args.map = True
map_args.start_date = (datetime.date.today() - datetime.timedelta(days=5*7)).strftime('%Y-%m-%d')
map_args.end_date = (datetime.date.today()).strftime('%Y-%m-%d')
main.main(map_args)
