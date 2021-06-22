# %%
import requests
import pandas as pd
from scripts.surveillance_region_map import name_to_region as region_mapping
from scripts.variant_to_pango import variant_to_lineage
from scripts.download_cases import cases_by_cw
# %%
url = 'https://cov-spectrum.ethz.ch/api/resource/sample2'
params = {'country': 'Switzerland',
          'fields': 'date,division,pangolinLineage'}
r = requests.get(url, params)
# %%
df = pd.read_json(r.content)
# Unroll the count compression, one row per sample
# df = df.loc[df.index.repeat(df['count'])]  # %%
df = df.set_index('date')
df = df.assign(region=lambda x: x.division.map(
    region_mapping).fillna(0).astype('int64'))
# %%


def lineage_to_variant(lineage):
    for val, keys in variant_to_lineage.items():
        if lineage in keys:
            return val
    return 'others'


df = df.assign(variant=lambda x: x.pangolinLineage.map(lineage_to_variant))
# %%
# Variants by surveillance region
variants_by_reg = pd.pivot_table(df, values='count', columns='variant', aggfunc='sum', fill_value=0, index=[
    'region', pd.Grouper(level='date', freq='W-MON', label='left', closed='left')])
variants_by_reg
# %%
# Add total sequences
seq_by_reg = pd.concat([variants_by_reg.loc[1:], pd.concat(
    {0: variants_by_reg.sum(level=1)}, names=['region'])]).sort_index()
seq_by_reg = seq_by_reg.assign(sequences=lambda x: x.sum(axis=1))
seq_by_reg
# %%
# Add cases from BAG dashboard
seq_and_cases = seq_by_reg.join(
    cases_by_cw(), how='outer').fillna(0).astype('int64')
seq_and_cases
# %%
# Export dataset
seq_and_cases.to_csv('data/cases_seq_by_cw_region.csv')
seq_and_cases.xs(1, level='region')[-10:]
