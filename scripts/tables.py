# %%
import datetime
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import os
# from .variant_to_pango import variant_to_lineage
from variant_to_pango import variant_to_lineage


def load_data():
    # df = pd.read_csv('data/cases_seq_by_cw_region.csv', parse_dates=['date'])
    df = pd.read_csv('../data/cases_seq_by_cw_region.csv',
                     parse_dates=['date'])
    df.set_index('date', inplace=True)
    return df


def restrict_dates(df, start, end):
    df = df.loc[(df.index >= start) & (df.index <= end)]
    return df


def save_fig(fig, output_dir, filename):
    results_dir = f'plots/{output_dir}'
    if not os.path.isdir(results_dir):
        os.makedirs(results_dir)
    fig.savefig(results_dir + '/' + filename)


# %%
df = load_data()
# %%
df = restrict_dates(df, datetime.datetime(2021, 4, 1),
                    datetime.datetime(2021, 5, 1))
# %%
df = df.reset_index()
df = df.set_index(['date', 'region'])
df.sum(level='region')
# %%
