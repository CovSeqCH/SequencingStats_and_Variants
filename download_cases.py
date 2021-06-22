# %%
import requests
import pandas as pd
import json
from surveillance_region_map import iso_canton_to_region as code_mapping


def cases_by_cw() -> pd.DataFrame:
    url = 'https://www.covid19.admin.ch/api/data/context'
    r = requests.get(url)
    r.status_code
    context = json.loads(r.content)
    url = context['sources']['individual']['json']['daily']['cases']
    df = pd.read_json(url, convert_dates=['datum'])
    df = df.assign(reg=lambda x: x.geoRegion.map(code_mapping))
    df = df.set_index(['reg', 'geoRegion', 'datum'])
    # cases = mi.sum(level=[0, 2]).entries.reset_index()
    # Slice to exclude -1 which includes Liechtenstein and the like
    cases_by_cw = df.groupby([pd.Grouper(level='reg'), pd.Grouper(
        level='datum', freq='W-MON', label='left', closed='left')]).sum().entries.loc[0:, slice(None)].rename('Cases')
    cases_by_cw = pd.DataFrame(cases_by_cw).rename_axis(index=['reg', 'date'])
    return cases_by_cw


cases_by_cw()
