#!/bin/bash
cd /cmap/bin/sushi/scripts/

# read in flagged arguments
case "$1" in
  fastq2readcount)
    shift
    Rscript fastq2readcount.R "${@}"
    ;;
  filter_counts)
    shift
    Rscript filter_counts.R "${@}"
    ;;
  CBnormalize)
    shift
    Rscript CBnormalize.R "${@}"
    ;;
  compute_l2fc)
    shift
    Rscript compute_l2fc.R "${@}"
    ;;
  collapse_replicates)
    shift
    Rscript collapse_replicates.R "${@}"
    ;;
  replicate_QC)
    shift
    Rscript replicate_QC.R "${@}"
    ;;
  seq_to_mts)
    shift
    source activate prism
    python seq_to_mts.py "${@}"
    ;;
  generate_biomarkers)
    shift
    Rscript generate_biomarkers.R "${@}"
    ;;
  --help|-h)
    printf "Available commands:"
    echo *.R | sed 's/.R /\n/g'
    ;;
  *)
    printf "Unknown parameter: %s \n" "$1"
    printf "Available commands:\n\t"
    echo *.R | sed -r 's/.R\s?/\n\t/g'
    exit -1
    shift
    ;;
esac

exit_code=$?
#echo "$exit_code"
exit $exit_code
