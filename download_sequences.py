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
len(r.content)
# %%
df = pd.read_json(r.content)
# %%
# Unroll the count compression, one row per sample
df = df.loc[df.index.repeat(df['count'])]
# %%
# Todo: Surveillance vs Non-Surveillance
df.plot('date',)
# %%
# Plots
# 1. Sequencing over time
# 2. Lineages over time
# 3.
# Need to create the views, resample by time or aggregate in some other way

df.describe()

# %%
df.columns
# %%
df.head()

# %%
df.pangolinLineage.value_counts()[0:60]
# %%
df.dtypes
# %%
df.pangolinLineage.value_counts()['AY.1']
# %%
