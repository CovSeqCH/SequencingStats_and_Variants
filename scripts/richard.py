import pandas as pd
import datetime
import matplotlib.pyplot as plt
from urllib.request import urlretrieve
import numpy as np
from CH_cases import parse, aggregates
​
​
tmp = {
    9: 817,
    10: 1179,
    11: 1230,
    12: 989,
    13: 1378,
    14: 1330,
    15: 1405,
    16: 1554,
    17: 1437
}
sequences_surveillance_program = pd.DataFrame(tmp, index=['sequences']).T
​
​
regions_to_analyze = list(aggregates.keys())
​


def CW_to_date(cw):
    return datetime.datetime.strptime(f"2021-W{cw}-1", '%G-W%V-%u')


​


def date_to_CW(d):
    return d.isocalendar()[1]


​


def datestr_to_CW(d):
    return date_to_CW(datetime.datetime.strptime(d, '%Y-%m-%d'))


​
# aws s3 cp  s3://nextstrain-ncov-private/metadata.tsv.gz data/metadata.tsv.gz
meta = pd.read_csv('data/metadata.tsv.gz', sep='\t', index_col=0)
​
# sequence data
CH_meta = meta.loc[(meta.country == 'Switzerland')
                   & (meta.date >= '2021-01-04')]
CH_meta.loc[CH_meta.division == 'Basel-Land', 'division'] = 'Basel-Landschaft'
CH_meta = CH_meta.loc[CH_meta.date.apply(lambda x:len(x) == 10)]
CH_meta.loc[:, "CW"] = CH_meta.date.apply(
    lambda x: date_to_CW(datetime.datetime.strptime(x, '%Y-%m-%d')))
CH_meta_total = CH_meta.groupby("CW").count()["virus"]
​
​
# case data
cases = parse()
for region in regions_to_analyze:
    cases[region] = cases[region].loc[cases[region].date > '2021-01-01']
    cases[region]['CW'] = cases[region].date.apply(datestr_to_CW)
    cases[region]['new_cases'] = [0] + list(np.diff(cases[region]['cases']))
​
cases_by_week = {}
sequences_by_week = {}
for region in regions_to_analyze:
    cases_by_week[region] = cases[region][[
        "CW", "new_cases"]].groupby('CW').sum()
    if region in aggregates:
        sequences_by_week[region] = CH_meta.loc[CH_meta.division.apply(
            lambda x: x in aggregates[region]), :].groupby("CW").count()["virus"]
    else:
        sequences_by_week[region] = CH_meta.loc[CH_meta.division == region, :].groupby(
            "CW").count()["virus"]
​
fig, ax1 = plt.subplots()
min_date = "2021-03-01"
max_date = "2021-05-02"
total_seqs = ((CH_meta.date >= min_date) & (CH_meta.date <= max_date)).sum()
plt.title(f"total {min_date}--{max_date}: {total_seqs}")
ax1.plot(CH_meta_total, label='Total CH')
ax1.plot(sequences_surveillance_program, label='Surveillance Program')
ax1.set_xlabel("calendar week")
ax1.set_xticks(np.arange(min(CH_meta_total.index),
               max(CH_meta_total.index), 2))
ax1.set_ylabel("sequences from Switzerland")
ax1.legend(loc=3)
ax2 = ax1.twinx()
ax2.plot(CH_meta_total.index, (CH_meta_total/cases_by_week['CH']['new_cases'])[
         CH_meta_total.index], c='C2', label='fraction sequenced')
ax2.set_ylabel("fraction of cases sequenced")
ax2.legend(loc=4)
ax2.set_ylim(0, 0.4)
ax1.set_xlim(1, 17)
plt.savefig('sequences_from_CH.png')
​
fig, ax1 = plt.subplots()
for ri, region in enumerate(regions_to_analyze):
    ax1.plot(sequences_by_week[region].index, (sequences_by_week[region]/cases_by_week[region]['new_cases'])[
             sequences_by_week[region].index], label=region, lw=3 if region == 'CH' else 2, c="k" if region == 'CH' else f"C{ri}")
    ax1.set_ylabel("fraction of cases sequenced")
ax1.set_xticks(np.arange(min(CH_meta_total.index),
               max(CH_meta_total.index), 2))
ax1.set_xlabel("calendar week")
ax1.set_xlim(1, 17)
plt.legend()
plt.savefig("fraction_sequenced_by_region.png")
​
​
fig, ax1 = plt.subplots()
for ri, region in enumerate(regions_to_analyze):
    ax1.plot(sequences_by_week[region].index, sequences_by_week[region], label=region,
             lw=3 if region == 'CH' else 2, c="k" if region == 'CH' else f"C{ri}")
    ax1.set_ylabel("sequences by region")
​
sum_of_regions = pd.DataFrame(
    [sequences_by_week[k] for k in sequences_by_week if k != 'CH']).fillna(0).sum(axis=0)
ax1.set_xticks(np.arange(min(CH_meta_total.index),
               max(CH_meta_total.index), 2))
ax1.set_xlabel("calendar week")
ax1.set_xlim(1, 17)
plt.legend()
plt.savefig("sequences_by_region.png")
​
​
fig = plt.figure()
# ,  "B.1.619",  "B.1.620", "B.1.1.318"]
variants = ["B.1.351", "P.1", "B.1.617.2", "B.1.617.1"]
all_VoCs = []
for v in variants:
    CH_meta_variant = CH_meta.loc[CH_meta.pango_lineage == v, :].groupby(["CW"]).count()[
        "virus"]
    CH_meta_variant.name = v
    all_VoCs.append(CH_meta_variant)
    plt.plot(CH_meta_variant, label=v)
​
plt.xlabel("calendar week")
plt.ylabel("sequences from Switzerland")
plt.xticks(np.arange(min(CH_meta_total.index), max(CH_meta_total.index), 2))
plt.xlim(1, 17)
plt.legend(loc=2, ncol=2)
plt.savefig('VoC_in_CH.png')
​
​
plt.figure()
all_VoCs = pd.DataFrame(all_VoCs).fillna(0)
plt.plot(all_VoCs.sum(axis=0)/CH_meta_total*100)
plt.xlabel("calendar week")
plt.ylabel("% of VoCs in Switzerland (excl B.1.1.7)")
plt.xticks(np.arange(min(CH_meta_total.index), max(CH_meta_total.index), 2))
plt.xlim(1, 17)
plt.savefig("fraction_VoC.png")
​
min_date = "2021-03-01"
max_date = "2021-03-31"
recent = CH_meta.loc[(CH_meta.pango_lineage.apply(lambda x:x in variants)) & (CH_meta.date >= min_date) & (
    CH_meta.date <= max_date), ["division", "pango_lineage", "virus"]].groupby(["pango_lineage", "division"]).count()
d = pd.concat([recent.loc[v] if v in recent.index else pd.DataFrame(
    [0]*len(recent.columns)) for v in variants], axis=1).fillna(0).astype(int)
d.columns = variants
d.loc[:, "total_VoC"] = d.sum(axis=1)
d.loc[:, "total_sequences"] = CH_meta.loc[(CH_meta.date >= min_date) & (
    CH_meta.date <= max_date), ["division", "virus"]].groupby(["division"]).count()["virus"]
d.loc["sum"] = d.sum(axis=0)
d.to_csv(f'VoC_last_{min_date}--{max_date}.tsv')
​
min_date = "2021-04-01"
max_date = "2021-04-30"
recent = CH_meta.loc[(CH_meta.pango_lineage.apply(lambda x:x in variants)) & (CH_meta.date >= min_date) & (
    CH_meta.date <= max_date), ["division", "pango_lineage", "virus"]].groupby(["pango_lineage", "division"]).count()
d = pd.concat([recent.loc[v] if v in recent.index else pd.DataFrame(
    [0]*len(recent.columns)) for v in variants], axis=1).fillna(0).astype(int)
d.columns = variants
d.loc[:, "total_VoC"] = d.sum(axis=1)
d.loc[:, "total_sequences"] = CH_meta.loc[(CH_meta.date >= min_date) & (
    CH_meta.date <= max_date), ["division", "virus"]].groupby(["division"]).count()["virus"]
d.loc["sum"] = d.sum(axis=0)
d.to_csv(f'VoC_last_{min_date}--{max_date}.tsv')
