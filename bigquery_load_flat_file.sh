bq mk variant_table.WI-20210121-strain-annotation-integer-impact-gene-pheno-v3;
  chmod 775 WI.20210121.strain.annotation-integer-impact-gene-pheno-v3.tsv;
  bq load --source_format=CSV  \
        --skip_leading_rows 1 --autodetect variant_table.WI-20210121-strain-annotation-integer-impact-gene-pheno-v3 \
        /Users/ryanmckeown/Documents/andersen_lab/annotation_concordance/cendr_incorporation/flat_file/WI.20210121.strain.annotation-integer-impact-gene-pheno-v3.tsv
