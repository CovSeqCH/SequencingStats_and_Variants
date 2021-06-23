import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import os
from .variant_to_pango import variant_to_lineage


def load_data():
    df = pd.read_csv('data/cases_seq_by_cw_region.csv', parse_dates=['date'])
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


def generate_plots(df, output_dir):
    data = df[df.region == 0]
    fig, ax = plt.subplots()
    ax.plot(data.index, data.sequences/data.cases)
    ax.axhline(0.1, ls=':')
    ax.set_ylim(0,)
    ax.set_title("Proportion of sequenced cases by week in Switzerland")
    ax.set_ylabel("Proportion of sequenced cases")
    ax.set_xlabel("First day of calendar week")
    ax.set_xticks(data.index)
    fig.autofmt_xdate()
    save_fig(fig, output_dir, 'sequence_share_CH.pdf')

    data = df[df.region == 0]
    fig, ax = plt.subplots()
    ax.plot(data.index, data.sequences)
    ax.axhline(1000, ls=':')
    ax.set_ylim(0,)
    ax.set_title("Sequenced cases by week in Switzerland")
    ax.set_ylabel("Sequenced cases")
    ax.set_xlabel("First day of calendar week")
    ax.set_xticks(data.index)
    fig.autofmt_xdate()
    save_fig(fig, output_dir, 'sequences_CH.pdf')

    fig, ax = plt.subplots()
    for i in range(1, 7):
        data = df[df.region == i]
        ax.plot(data.index, data.sequences/data.cases, label=i)
    ax.axhline(0.1, ls=':')
    ax.set_title(
        "Proportion of sequenced cases by week and surveillance region")
    ax.set_ylabel("Proportion of sequenced cases")
    ax.set_xlabel("First day of calendar week")
    ax.set_ylim(0,)
    ax.legend()
    ax.set_xticks(data.index)
    fig.autofmt_xdate()
    save_fig(fig, output_dir, 'sequence_share_regions.pdf')

    data = df[df.region == 0]
    fig, ax = plt.subplots()
    for key in variant_to_lineage.keys():
        ax.plot(data.index, data[key]/data.sequences, label=key)
    ax.set_title("Share of variants by week in Switzerland")
    ax.set_ylabel("Share of all sequences")
    ax.set_xlabel("First day of calendar week")
    ax.set_ylim(0,)
    ax.legend()
    ax.set_xticks(data.index)
    fig.autofmt_xdate()
    save_fig(fig, output_dir, 'variant_share_CH.pdf')
