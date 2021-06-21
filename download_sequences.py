# %%
import requests
import pandas as pd
# %%
url = 'https://cov-spectrum.ethz.ch/api/resource/sample2'
params = {'country': 'Switzerland',
          'fields': 'date,region,country,division,ageGroup,sex,hospitalized,deceased,pangolinLineage'}
r = requests.get(url, params)
r.status_code
# %%
df = pd.read_json(r.content)
# Unroll the count compression, one row per sample
df = df.loc[df.index.repeat(df['count'])]
# Todo: Surveillance vs Non-Surveillance
seq = df
seq = seq.set_index('date')
seq.head()
# %%
# Sequences per calendar week (week beginning and including Monday)
seq_count_by_week = pd.DataFrame(seq.resample(
    'W-MON', label='left', closed='left').country.count())
seq_count_by_week.rename(columns={'country': 'sequences'}, inplace=True)
seq_count_by_week.index.rename('start_of_cw', inplace=True)
seq_count_by_week = seq_count_by_week.assign(
    calendar_week=lambda x: x.index.weekofyear)
seq_count_by_week
# %%
# Sequences per calendar month (labelled by last day of month)
seq_count_by_week = pd.DataFrame(seq.resample('M').country.count())
seq_count_by_week.rename(columns={'country': 'sequences'}, inplace=True)
seq_count_by_week.index.rename('end_of_month', inplace=True)
seq_count_by_week = seq_count_by_week.assign(month=lambda x: x.index.month)
seq_count_by_week
# %%
# # %%
# # %%
# # Plots
# # 1. Table of sequences by calendar week

# # %%
# df.columns
# # %%
# df.head()

# # %%
# df.pangolinLineage.value_counts()[0:60]
# # %%
# df.dtypes
# # %%
# df.pangolinLineage.value_counts()['AY.1']
# # %%
