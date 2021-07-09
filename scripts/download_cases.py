# %%
import requests
import pandas as pd
import json
from scripts.surveillance_region_map import iso_canton_to_region as code_mapping


def cases_by_cw() -> pd.DataFrame:
    url = 'https://www.covid19.admin.ch/api/data/context'
    r = requests.get(url)
    url = json.loads(r.content)[ 'sources']['individual']['json']['daily']['cases']
    df = pd.read_json(url, convert_dates=['datum'])
    df = df.assign(region=lambda x: x.geoRegion.map(code_mapping))
    df = df.set_index(['region', 'geoRegion', 'datum'])
    # Slice to exclude -1 which includes Liechtenstein and the like
    cases_by_cw = df.groupby([pd.Grouper(level='region'), pd.Grouper(
        level='datum', freq='W-MON', label='left', closed='left')]).sum().entries.loc[0:, slice(None)].rename('cases')
    cases_by_cw = pd.DataFrame(cases_by_cw).rename_axis(
        index=['region', 'date'])
    return cases_by_cw

def cases_by_day() -> pd.DataFrame:
    url = 'https://www.covid19.admin.ch/api/data/context'
    r = requests.get(url)
    url = json.loads(r.content)['sources']['individual']['json']['daily']['cases']
    df = pd.read_json(url, convert_dates=['datum'])
    df = df.assign(region=lambda x: x.geoRegion.map(code_mapping))
    df = df.set_index(['region', 'geoRegion', 'datum'])
    # Slice to exclude -1 which includes Liechtenstein and the like
    cases_by_day = df.groupby([pd.Grouper(level='region'), pd.Grouper(
        level='datum', freq='D', label='left', closed='left')]).sum().entries.loc[0:, slice(None)].rename('cases')
    cases_by_day = pd.DataFrame(cases_by_day).rename_axis(
        index=['region', 'date'])
    return cases_by_day
