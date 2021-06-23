# data-ingest

Download raw data and turn it into useful dataset.

Private repo at first, but may be made public, or at least the dataset that are generated.

Analyses and graph generation will have their own repo.

# How to run yourself
Do fresh data download as follows:
```
./main.py --download
```

Sample usage to create plots for May in folder `/plots/may/`:
```
./main.py --plot --start-date 2021-04-01 --end-date 2021-06-07 --output-dir may
```

The output `.csv` is produced in the script /scripts/download_sequences.py

The plots are produced in the script /scripts/plots.py

Variant definitions are in /scripts/variant_to_pango.py

Region definitions in /scripts/surveillance_region_map.py

# CSV definitions

- **reg**: `int` for surveillance region, 0 stands for all of Switzerland
- **date**: Beginning=Monday of respective calendar week
- **Alpha,Beta,...**: Variants according to WHO definition
- **Others**: Anything that's not accounted for in named variant columns
- **Sequences**: Total of all sequences
- **Cases**: Case count according to BAG dashboard

# Todo

- [x] Get sequence metadata from covSpectrum
- [x] Get case data from BAG website
- [x] Turn pango lineages into WHO VOC designations
- [x] Add region 1-6 definitions
- [x] Define output dataset
