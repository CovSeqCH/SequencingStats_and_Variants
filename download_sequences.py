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
df = df.loc[df.index.repeat(df['count'])]  # %%
df = df.set_index('date')
# Region definitions
region_mapping = {
    'Basel-Stadt': 3,
    'Solothurn': 3,
    'Aargau': 3,
    'Vaud': 1,
    'Bern': 2,
    'Basel-Land': 3,
    'Zürich': 5,
    'Schwyz': 4,
    'Thurgau': 5,
    'Neuchâtel': 1,
    'Schaffhausen': 5,
    'Lucerne': 4,
    'Graubünden': 6,
    'Uri': 4,
    'Geneva': 1,
    'Ticino': 6,
    'Valais': 1,
    'Sankt Gallen': 5,
    'Switzerland': pd.NA,
    'Obwalden': 4,
    'Jura': 2,
    'Fribourg': 2,
    'Zug': 4,
    'Glarus': 5,
    'Nidwalden': 4,
    'Appenzell Ausserrhoden': 5,
    'VALAIS': 1,
    'Graub�nden': 6
}
df = df.assign(reg=lambda x: x.division.map(region_mapping))
# Todo: Surveillance vs Non-Surveillance
# %%


def select_region(df, region=0) -> pd.DataFrame:
    out = df.copy(deep=True)
    if region == 0:
        return out
    return out[out.reg == region]


seq = select_region(df)
seq.head()
# %%
# Sequences per calendar week (week beginning and including Monday)
# Select e.g. region 2: seq = select_region(df,2)
seq_count_by_week = pd.DataFrame(seq.resample(
    'W-MON', label='left', closed='left').country.count())
seq_count_by_week.rename(columns={'country': 'sequences'}, inplace=True)
seq_count_by_week.index.rename('start_of_cw', inplace=True)
seq_count_by_week = seq_count_by_week.assign(
    calendar_week=lambda x: x.index.weekofyear)
seq_count_by_week
# %%
fig = seq_count_by_week.sequences.plot()
fig.set_title("Sequences by calendar week")
fig.set_ylabel("# of sequences")
# %%
# Sequences per calendar month (labelled by last day of month)
seq_count_by_month = pd.DataFrame(seq.resample('M').country.count())
seq_count_by_month.rename(columns={'country': 'sequences'}, inplace=True)
seq_count_by_month.index.rename('end_of_month', inplace=True)
seq_count_by_month = seq_count_by_month.assign(month=lambda x: x.index.month)
seq_count_by_month
# %%
fig = seq_count_by_month.sequences.plot()
fig.set_title("Sequences by calendar month")
fig.set_ylabel("# of sequences")

# Now can run above plots for each region

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

# %%
