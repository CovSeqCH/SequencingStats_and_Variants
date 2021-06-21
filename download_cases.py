# %%
import requests
import pandas as pd
import json
import pprint
# %%
url = 'https://www.covid19.admin.ch/api/data/context'
r = requests.get(url)
r.status_code
# %%
context = json.loads(r.content)
# %%
url = context['sources']['individual']['json']['daily']['cases']
# %%
df = pd.read_json(url, convert_dates=['datum'])
# %%
code_mapping = {
    'CH': 0,
    'CHFL': -1,
    'AG': 3,
    'AI': 5,
    'AR': 5,
    'BE': 2,
    'BL': 3,
    'BS': 3,
    'FL': -1,
    'FR': 2,
    'GE': 1,
    'GL': 5,
    'GR': 6,
    'JU': 2,
    'LU': 4,
    'NE': 1,
    'NW': 4,
    'OW': 4,
    'SG': 5,
    'SH': 5,
    'SO': 3,
    'SZ': 4,
    'TG': 5,
    'TI': 6,
    'UR': 4,
    'VD': 1,
    'VS': 1,
    'ZG': 4,
    'ZH': 5
}
# %%
df = df.assign(reg=lambda x: x.geoRegion.map(code_mapping))
mi = df.set_index(['reg', 'geoRegion', 'datum'])
cases = mi.sum(level=[0, 2]).entries.reset_index()
cases_by_cw = mi.groupby([pd.Grouper(level='reg'), pd.Grouper(
    level='datum', freq='W-MON', label='left', closed='left')]).sum().entries.loc[0:, slice(None)]
cases_by_cw.reset_index()
# %%
# dataset sequencing activity and regions
## cw, reg, seq, cases, (variants,...)
