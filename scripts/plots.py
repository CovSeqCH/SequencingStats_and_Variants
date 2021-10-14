import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import os
import datetime as dt
from .variant_to_pango import variant_to_lineage
# from variant_to_pango import variant_to_lineage


def load_data():
    df = pd.read_csv('data/cases_seq_by_cw_region.csv', parse_dates=['date'])
    # df = pd.read_csv('../data/cases_seq_by_cw_region.csv',
    #  parse_dates=['date'])
    df.set_index('date', inplace=True)
    return df


def restrict_dates(df, start, end):
    df = df.loc[(df.index >= start) & (df.index < end)]
    return df


def save_fig(fig, output_dir, filename):
    file_types = ['pdf', 'png']
    for file_type in file_types:
        results_dir = f'plots/{output_dir}/{file_type}'
        if not os.path.isdir(results_dir):
            os.makedirs(results_dir)
        fig.savefig(f'{results_dir}/{filename}.{file_type}', dpi=400)


def generate_plots(df, output_dir):
    data = df[df.region == 0]

    today = dt.date.today()
    monday_4_weeks_ago = today + dt.timedelta(days=-today.weekday(), weeks=-4)

    fig, ax1 = plt.subplots()
    ax1.plot(data.index, data.sequences/data.cases, label='Fraction sequenced')
    ax1.set_ylim(0,)
    ax1.set_title("Swiss sequences by calendar week of sample")
    ax1.set_ylabel("Fraction sequenced")
    ax1.set_xlabel(f"Sampling date (Calendar week) - Generated: {today}")
    ax1.axvspan(monday_4_weeks_ago, today + dt.timedelta(days=10), color='grey', alpha=0.3, label='Incomplete data')
    ax1.legend(loc=3)
    ax2 = ax1.twinx()
    ax2.plot(data.index, data.sequences,
             label='Total number of sequences', color='g', ls='-.')
    ax2.plot(data.index, data.in_program,
             label='Sequences by consortium', color='r', ls=':')
    ax2.set_ylabel("Number of sequences")
    ax2.set_ylim(0,)
    ax2.legend(loc=4)
    ax1.set_xlim(df.index.min(), df.index.max())
    ax1.set_xticks(data.index)
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%-W'))
    save_fig(fig, output_dir, 'sequence_share_CH')

    fig, ax = plt.subplots()
    for i in range(1, 7):
        data = df[df.region == i]
        ax.plot(data.index, data.sequences/data.cases, label=f'Region {i}')
    ax.axhline(0.1, ls=':')
    ax.axvspan(monday_4_weeks_ago, today + dt.timedelta(days=10), color='grey', alpha=0.3, label='Incomplete data')
    ax.set_title(
        "Case fraction sequenced by sample week and surveillance region")
    ax.set_ylabel("Proportion of cases sequenced")
    ax.set_xlabel(f"Sampling date (Calendar week) - Generated: {today}")
    ax.set_ylim(0,)
    ax.set_xlim(df.index.min(), df.index.max())
    ax.legend()
    ax.set_xticks(data.index)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%-W'))
    save_fig(fig, output_dir, 'sequence_share_regions')

    fig, ax = plt.subplots()
    for i in range(1, 7):
        data = df[df.region == i]
        ax.plot(data.index, data.sequences, label=f'Region {i}')
    ax.axvspan(monday_4_weeks_ago, today + dt.timedelta(days=10), color='grey', alpha=0.3, label='Incomplete data')
    ax.set_title("Sequenced cases by sample week and surveillance region")
    ax.set_ylabel("Total number of sequences")
    ax.set_xlabel(f"Sampling date (Calendar week) - Generated: {today}")
    ax.set_ylim(0,)
    ax.set_xlim(df.index.min(), df.index.max())
    ax.legend()
    ax.set_xticks(data.index)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%-W'))
    save_fig(fig, output_dir, 'sequences_regions')

    data = df[df.region == 0]
    fig, ax = plt.subplots()
    for key in variant_to_lineage.keys():
        ax.plot(data.index, data[key]/data.sequences,
                label=key)
    ax.axvspan(monday_4_weeks_ago, today + dt.timedelta(days=10), color='grey', alpha=0.3, label='Incomplete data')
    ax.set_title("Fraction of variants by sample week in Switzerland")
    ax.set_ylabel("Fraction of sequences")
    ax.set_xlabel(f"Sampling date (Calendar week) - Generated: {today}")
    ax.set_ylim(0,)
    ax.set_xlim(df.index.min(), df.index.max())
    ax.legend()
    ax.set_xticks(data.index)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%-W'))
    save_fig(fig, output_dir, 'variant_share_CH')

    data = df[df.region == 0]
    fig, ax = plt.subplots()
    for key in variant_to_lineage.keys():
        ax.plot(data.index, data[key]/data.sequences*data.cases,
                label=key)
    ax.axvspan(monday_4_weeks_ago, today + dt.timedelta(days=10), color='grey', alpha=0.3, label='Incomplete data')
    ax.set_title(
        "Estimated number of cases per week per variant in Switzerland")
    ax.set_ylabel("Estimated number of cases (fraction x cases)")
    ax.set_xlabel(f"Sampling date (Calendar week) - Generated: {today}")
    ax.set_ylim(0,)
    ax.set_xlim(df.index.min(), df.index.max())
    ax.legend()
    ax.set_xticks(data.index)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%-W'))
    save_fig(fig, output_dir, 'variant_estimate_CH')
