set -x
aws s3 cp  s3://nextstrain-ncov-private/metadata.tsv.gz - | \
gzcat - | \
tsv-filter -H --str-in-fld country:Switzerland | \
tsv-filter -H --regex pango_lineage:'(B.1.1.529|BA.)' | \
tsv-select -H -f strain,date,date_submitted,country,division,Nextstrain_clade,pango_lineage,originating_lab,submitting_lab  | \
keep-header - -- sort -k2 \
>omicron_swiss.tsv 
