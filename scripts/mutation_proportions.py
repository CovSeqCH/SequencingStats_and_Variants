#%%
import json
from functools import reduce
from itertools import product
from pdb import set_trace
from typing import Dict

import pandas as pd
import requests

#%%
# List of mutations to watch

# sotrovimab mutations of interest (as provided by Pauline Vetter)
mois = ["", "S:377", "S:337H", "S:337L", "S:337R", "S:337T", "S:340", "S:340A", "S:340K", "S:340G", "S:356", "S:356T"]

regions = [{}, {"region": "Europe"}, {"country": "Switzerland"}]

combinations = list(product(mois, regions))
#%%
ACCESS_KEY = "9Cb3CqmrFnVjO3XCxQLO6gUnKPd"
def get_mutation_numbers(mutation: str = "", geo_filter: Dict = {}) -> pd.DataFrame:
    """
    Get mutation numbers for a given mutation and a geo filter
    """
    url = "https://cov-spectrum.ethz.ch/gisaid/api/v1/sample/aggregated"
    params = {"fields": "date", "aaMutations": mutation, "dateFrom": "2019-12-01","accessKey": ACCESS_KEY}
    params.update(geo_filter)
    r = requests.get(url, params)
    j = json.loads(r.text)
    data = json.dumps(j["data"])
    if data == "[]":
        data = '[{"date":"2020-12-21","count":0}]'
    df = pd.read_json(data)
    df["date"] = pd.to_datetime(df["date"], errors="coerce")
    df = df[df.date.notnull()]
    df = df.set_index("date")
    geo = (list(geo_filter.values()) or ["Global"])[0]
    mutation_name = mutation or "all"
    df.rename(columns={"count": f"{geo}_{mutation_name}"}, inplace=True)
    return df.groupby(
        [pd.Grouper(level="date", freq="W-MON", label="left", closed="left", sort=True)]
    ).sum()


#%%
data_frames = list(
    map(lambda x: get_mutation_numbers(mutation=x[0], geo_filter=x[1]), combinations)
)
#%%
df_merged = (
    reduce(lambda left, right: pd.merge(left, right, on="date", how="outer"), data_frames)
    .fillna(0)
    .astype(int)
)
df_merged.to_csv("data/mutation_counts.csv")
