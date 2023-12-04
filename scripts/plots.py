# %%
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.ticker as mtick
import os
import datetime as dt
from variant_to_pango import variant_to_lineage
from matplotlib.colors import to_rgb


def load_data():
    df = pd.read_csv("data/cases_seq_by_cw_region.csv", parse_dates=["date"])
    # df = pd.read_csv('../data/cases_seq_by_cw_region.csv',
    #  parse_dates=['date'])
    df.set_index("date", inplace=True)
    return df


def restrict_dates(df, start, end):
    df = df.loc[(df.index >= start) & (df.index < end)]
    return df


def save_fig(fig, output_dir, filename):
    file_types = ["pdf", "png"]
    for file_type in file_types:
        results_dir = f"plots/{output_dir}/{file_type}"
        if not os.path.isdir(results_dir):
            os.makedirs(results_dir)
        fig.savefig(f"{results_dir}/{filename}.{file_type}", dpi=400)


# %%
# %%
def generate_plots(df, output_dir):
    data = df[df.region == 0]

    today = dt.date.today()
    monday_4_weeks_ago = today + dt.timedelta(days=-today.weekday(), weeks=-4)

    fig, ax1 = plt.subplots(figsize=(7, 4))
    ax1.set_title("Swiss sequences by calendar week of sample collection date")
    ax1.set_xlabel(f"Sampling date - Generated: {today}")
    ax1.plot(
        data.index,
        data.sequences,
        label="Total number of sequences",
        color="g",
        ls="-.",
    )
    ax1.plot(
        data.index, data.in_program, label="Sequences by consortium", color="r", ls=":"
    )
    ax1.set_ylabel("Number of sequences")
    ax1.set_ylim(
        0,
    )
    ax1.set_xlim(df.index.min(), df.index.max())
    for label in ax1.get_xticklabels(which="major"):
        label.set(rotation=30, horizontalalignment="right")
    ax1.legend()
    plt.tight_layout()
    save_fig(fig, output_dir, "sequence_share_CH")

    fig, ax = plt.subplots(figsize=(7, 4))
    for i in range(1, 7):
        data = df[df.region == i]
        ax.plot(data.index, data.sequences, label=f"Region {i}")
    ax.set_title("Sequenced cases by sample week and surveillance region")
    ax.set_ylabel("Total number of sequences")
    ax.set_xlabel(f"Sampling date - Generated: {today}")
    ax.set_ylim(
        0,
    )
    ax.set_xlim(df.index.min(), df.index.max())
    for label in ax.get_xticklabels(which="major"):
        label.set(rotation=30, horizontalalignment="right")
    ax.legend()
    plt.tight_layout()
    save_fig(fig, output_dir, "sequences_regions")

    data = df[df.region == 0]
    fig, ax = plt.subplots(figsize=(7, 4))
    print(data)
    variants = list(variant_to_lineage.keys())
    print(variants)
    variants.append("others")
    for key in variants:
        if key in data.columns:
            if data[key].sum() > 10:
                try:
                    ax.plot(data.index, data[key] / data.sequences, label=key)
                except KeyError:
                    pass
    ax.set_title("Fraction of variants by sample week in Switzerland")
    ax.set_ylabel("Fraction of sequences")
    ax.yaxis.set_major_formatter(mtick.PercentFormatter(1.0))
    ax.set_xlabel(f"Sampling date (Calendar week) - Generated: {today}")
    ax.set_ylim(
        0,
    )
    ax.set_xlim(df.index.min(), df.index.max())
    for label in ax.get_xticklabels(which="major"):
        label.set(rotation=30, horizontalalignment="right")
    ax.legend(bbox_to_anchor=(1.04, 1), loc="upper left")
    # ax.legend()
    plt.tight_layout()
    save_fig(fig, output_dir, "variant_share_CH")

    data = df[df.region == 0]
    fig, ax = plt.subplots(figsize=(7, 4))
    variants = list(variant_to_lineage.keys())
    variants.append("others")
    for key in variants:
        if key in data.columns:
            if data[key].sum() > 10:
                try:
                    ax.plot(
                        data.index, data[key] / data.sequences * data.cases, label=key
                    )
                except KeyError:
                    pass
    ax.axvspan("2023-01-01", "2023-12-31", alpha=0.2, color="grey", label="No case data")
    ax.set_title("Estimated number of cases per week per variant in Switzerland")
    ax.set_ylabel("Estimated number of cases (fraction x cases)")
    ax.set_xlabel(f"Sampling date (Calendar week) - Generated: {today}")
    ax.set_ylim(
        0,
    )
    ax.set_xlim(df.index.min(), df.index.max())
    for label in ax.get_xticklabels(which="major"):
        label.set(rotation=30, horizontalalignment="right")
    ax.legend()
    plt.tight_layout()
    save_fig(fig, output_dir, "variant_estimate_CH")


# %%
# Phase 1
df = load_data()
df
df = restrict_dates(df, "2021-03-01", "2022-03-31")
output_dir = "summary/phase1"
generate_plots(df, output_dir)
# %%

# Phase 2
df = load_data()
df = restrict_dates(df, "2022-04-01", "2022-12-31")
output_dir = "summary/phase2"
generate_plots(df, output_dir)
# %%
# Phase 3
df = load_data()
df = restrict_dates(df, "2023-01-01", "2023-11-10")
output_dir = "summary/phase3"
generate_plots(df, output_dir)
# %%

# Entire period
df = load_data()
df = restrict_dates(df, "2021-03-01", "2023-11-10")
output_dir = "summary/entire_period"
generate_plots(df, output_dir)

# %%
