#%%
import datetime as dt
import json

import pandas as pd
import requests
from pango_aliasor.aliasor import Aliasor

from scripts.download_cases import cases_by_cw
from scripts.participating_labs import excluded_labs
from scripts.surveillance_region_map import name_to_region as region_mapping
from scripts.variant_to_pango import variant_to_lineage

#%%
ACCESS_KEY = "9Cb3CqmrFnVjO3XCxQLO6gUnKPd"

def generate_csv():
    url = "https://lapis.cov-spectrum.org/gisaid/v1/sample/aggregated"
    params = {"country": "Switzerland", "fields": "date,division,nextcladePangoLineage,submittingLab", "accessKey": ACCESS_KEY}
    r = requests.get(url, params)
    j = json.loads(r.text)
    df = pd.read_json(json.dumps(j["data"]))
    df = df.set_index("date")
    df = df.assign(region=lambda x: x.division.map(region_mapping).fillna(0).astype("int64"))

    # %%
    aliasor = Aliasor()

    def lineage_to_variant(lineage):
        """
        Takes Pango lineage string and returns variant name
        Uses variant to lineage mapping from variant_to_pango.py
        Most specific lineages are at top
        Expand via aliasor and see if it matches
        """
        for val, keys in variant_to_lineage.items():
            for key in keys:
                if str(aliasor.uncompress(lineage)).startswith(key):
                    return val
        if debug:
            print(f'{lineage} assigned "others"')
        # if lineage is None or lineage == "Unassigned":
        #     return "None"
        return "others"

    debug = False
    recent = df[df.index > "20210830"]
    recent.assign(variant=lambda x: x.nextcladePangoLineage.map(lineage_to_variant))
    debug = False

    df = df.assign(variant=lambda x: x.nextcladePangoLineage.map(lineage_to_variant))

    #%%
    def lab_to_program(lab):
        return "in_program" if lab not in excluded_labs else "not_in_program"

    df = df.assign(program=lambda x: x.submittingLab.map(lab_to_program))
    # %%
    # Variants by submitting lab
    variants_by_lab = pd.pivot_table(
        df,
        values="count",
        columns="program",
        aggfunc="sum",
        fill_value=0,
        index=["region", pd.Grouper(level="date", freq="W-MON", label="left", closed="left")],
    )
    variants_by_lab
    # %%
    # Variants by surveillance region
    variants_by_reg = pd.pivot_table(
        df,
        values="count",
        columns="variant",
        aggfunc="sum",
        fill_value=0,
        index=["region", pd.Grouper(level="date", freq="W-MON", label="left", closed="left")],
    )
    # %%
    # Add total sequences
    variants = variants_by_reg.join(variants_by_lab, how="outer")
    seq_by_reg = pd.concat(
        [variants.loc[1:], pd.concat({0: variants.sum(level=1)}, names=["region"])]
    ).sort_index()
    seq_by_reg = seq_by_reg.assign(sequences=lambda x: seq_by_reg.iloc[:, :-2].sum(axis=1))
    # %%
    # Add cases from BAG dashboard
    seq_and_cases = seq_by_reg.join(cases_by_cw(), how="outer").fillna(0).astype("int64")
    # seq_and_cases
    # %%
    # Export dataset
    seq_and_cases.to_csv("data/cases_seq_by_cw_region.csv")
    # %%
    variants_by_week = pd.pivot_table(
        df,
        values="count",
        columns="nextcladePangoLineage",
        aggfunc="sum",
        fill_value=0,
        index=[pd.Grouper(level="date", freq="W-MON", label="left", closed="left")],
    )
    variants_by_week.to_csv("data/variants_by_week.csv")
    # variants_by_week.set_index('date')
    recent_common_variants = (
        variants_by_week[-10:].sum(axis=0)[1:].sort_values(ascending=False)[0:50]
    )
    variants_by_week["total"] = variants_by_week.sum(axis=1)
    variants_by_week["others"] = variants_by_week["total"] - variants_by_week[
        recent_common_variants.index
    ].sum(axis=1)
    # print(variants_by_week)
    list_of_common_variants = recent_common_variants.index.tolist()
    variants_by_week["total_common"] = variants_by_week[list_of_common_variants].sum(axis=1)
    variants_by_week["others"] = variants_by_week["total"] - variants_by_week["total_common"]
    list_of_common_variants.extend(["others", "total"])
    print(f"Most common lineages: {list_of_common_variants}")
    variants_by_week.to_csv(
        "data/recent_common_variants_by_week.csv", columns=list_of_common_variants
    )
    variants_by_week.div(variants_by_week["total"], axis=0).to_csv(
        "data/recent_common_variants_by_week_relative.csv",
        columns=recent_common_variants.index,
        float_format="%.4f",
    )
    # seq_and_cases.xs(1, level='region')[-10:]
