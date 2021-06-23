# %%
from datetime import datetime
import pandas as pd
import matplotlib.pyplot as plt
import os
# %%


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
    gen = df[df.region == 0]
    fig, ax = plt.subplots()
    ax.plot(gen.index, gen.sequences/gen.cases)
    ax.axhline(0.1, ls=':')
    ax.set_ylim(0,)
    ax.set_title("Proportion of sequenced cases by week in Switzerland")
    ax.set_ylabel("Proportion of sequenced cases")
    fig.autofmt_xdate()
    save_fig(fig, output_dir, 'sequence_share_CH.pdf')
    # %%
    gen = df[df.region == 0]
    fig, ax = plt.subplots()
    ax.plot(gen.index, gen.sequences)
    ax.axhline(1000, ls=':')
    ax.set_ylim(0,)
    ax.set_title("Sequenced cases by week in Switzerland")
    ax.set_ylabel("Sequenced cases")
    fig.autofmt_xdate()
    save_fig(fig, output_dir, 'sequences_CH.pdf')
    # %%
    fig, ax = plt.subplots()
    for i in range(1, 7):
        gen = df[df.region == i]
        ax.plot(gen.index, gen.sequences/gen.cases, label=i)
    ax.axhline(0.1, ls=':')
    ax.set_title(
        "Proportion of sequenced cases by week and surveillance region")
    ax.set_ylabel("Proportion of sequenced cases")
    ax.set_ylim(0,)
    ax.legend()
    fig.autofmt_xdate()
    save_fig(fig, output_dir, 'sequence_share_regions.pdf')
    # %%
    gen = df[df.region == 0]
    fig, ax = plt.subplots()
    ax.bar(gen.index, gen.delta/gen.sequences)
    ax.set_title("Share of Delta of sequenced cases by week in Switzerland")
    ax.set_ylabel("Share of Delta")
    fig.autofmt_xdate()
    save_fig(fig, output_dir, 'share_CH.pdf')
