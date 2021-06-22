# data-ingest

Download raw data and turn it into useful dataset.

Private repo at first, but may be made public, or at least the dataset that are generated.

Analyses and graph generation will have their own repo.

# CSV definitions

- *reg*: `int` for surveillance region, 0 stands for all of Switzerland
- *date*: Beginning=Monday of respective calendar week
- *Alpha,Beta,...*: Variants according to WHO definition
- *Others*: Anything that's not accounted for in named variant columns
- *Sequences*: Total of all sequences
- *Cases*: Case count according to BAG dashboard

# Todo

- [ ] Get sequence metadata from covSpectrum
- [ ] Get case data from BAG website
- [ ] Turn pango lineages into WHO VOC designations
- [ ] Add region 1-6 definitions
- [ ] Define output dataset
