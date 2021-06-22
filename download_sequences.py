# %%
import requests
import pandas as pd
import numpy as np
from surveillance_region_map import name_to_region as region_mapping
from variant_to_pango import variant_to_lineage
from download_cases import cases_by_cw
import matplotlib.pyplot as plt
# %%
url = 'https://cov-spectrum.ethz.ch/api/resource/sample2'
params = {'country': 'Switzerland',
          'fields': 'date,division,pangolinLineage'}
r = requests.get(url, params)
# %%
df = pd.read_json(r.content)
# Unroll the count compression, one row per sample
df = df.loc[df.index.repeat(df['count'])]  # %%
df = df.set_index('date')
df = df.assign(reg=lambda x: x.division.map(
    region_mapping).fillna(0).astype('int64'))
# %%


def lineage_to_variant(lineage: str) -> str:
    for val, keys in variant_to_lineage.items():
        if lineage in keys:
            return val
    return 'others'


df = df.assign(variant=lambda x: x.pangolinLineage.map(lineage_to_variant))
# %%
# Variants by surveillance region
variants_by_reg = pd.pivot_table(df, values='pangolinLineage', columns='variant', aggfunc='count', fill_value=0, index=[
    'reg', pd.Grouper(level='date', freq='W-MON', label='left', closed='left')])
variants_by_reg
# %%
# Add total sequences
seq_by_reg = pd.concat([variants_by_reg.loc[1:], pd.concat(
    {0: variants_by_reg.sum(level=1)}, names=['reg'])]).sort_index()
seq_by_reg = seq_by_reg.assign(sequences=lambda x: x.sum(axis=1))
# %%
# Add cases from BAG dashboard
seq_and_cases = seq_by_reg.join(
    cases_by_cw(), how='outer').fillna(0).astype('int64')
# %%
# Export dataset
seq_and_cases.to_csv('cases_seq_by_cw_region.csv')
seq_and_cases.xs(1, level='reg')[-10:]
# %%
# Plots
gen = seq_and_cases.loc[1].iloc[-20:]
fig, ax = plt.subplots()
ax.plot(gen.index, gen.Sequences/gen.Cases)
ax.set_title("Proportion of sequenced cases by week in region 1")
ax.set_ylabel("Proportion of sequenced cases")
fig.autofmt_xdate()
# %%
gen = seq_and_cases.loc[0].iloc[-10:]
fig, ax = plt.subplots()
ax.scatter(gen.index, gen.Delta/gen.Sequences)
ax.set_title("Share of Delta of sequenced cases by week in Switzerland")
ax.set_ylabel("Share of Delta")
fig.autofmt_xdate()
# %%
fig = gen.Sequences.plot()
fig.set_title("Sequences from region 1 (Geneva etc.) by calendar week")
fig.set_ylabel("# of sequences")
