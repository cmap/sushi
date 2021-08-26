#!/bin/bash
source /broad/software/scripts/useuse
reuse .bcl2fastq2-v2.20.0


OUT_DIR=/xchip/prism/bcl2fastq/
PSEQ_obelix=/cmap/obelix/pod/prismSeq/

#optional
if test $# -lt 1; then
  printf "Usage ./bcl2fastq.sh [options]\nOptions include:\n"
  printf -- "-s, --seq_code \t\t Sequencer run code E.g. JNV9V \n"
  printf -- "-p, --proj_code \t Project code to prep fastq/ directory in /cmap/obelix \n"
  printf -- "-o, --out_dir \t Path to output files/folders \n"
  printf -- "-h, --help \t\t Print this help text\n"
  exit 1
fi

while test $# -gt 0; do
  case "$1" in
    -h|--help)
      printf "Usage ./bcl2fastq.sh [options]\nOptions include:\n"
      printf -- "-s, --seq_code \t\t Sequencer run code E.g. JNV9V \n"
      printf -- "-p, --proj_code \t Project code to prep fastq/ directory in /cmap/obelix \n"
      printf -- "-o, --out_dir \t Path to output files/folders \n"
      printf -- "-h, --help \t\t Print this help text\n"
      exit 0
      ;;
    -s|--seq_code)
      shift
      SEQ_CODE=$1
      ;;
    -p|--proj_code)
      shift
      #echo $1
      PROJ_CODE=$1
      ;;
    -o|--out_dir)
      shift
      #echo $1
      OUT_DIR=$1
      ;;
    *)
      printf "Unknown parameter: %s \n" "$1"
      shift
      ;;
  esac
  shift
done

if [ ! -d $OUT_DIR ]
then
  mkdir $OUT_DIR
fi

RUNFOLDER_DIR=$(echo /xchip/prism/MiSeq\ Outputs/*-$SEQ_CODE)

echo $RUNFOLDER_DIR

bcl2fastq --runfolder-dir "$RUNFOLDER_DIR" --output-dir $OUT_DIR/$SEQ_CODE --minimum-trimmed-read-length 35 --mask-short-adapter-reads 22 --create-fastq-for-index-reads

if [ ! -d $PSEQ_obelix/$PROJ_CODE ]
then
  mkdir $PSEQ_obelix/$PROJ_CODE
  mkdir $PSEQ_obelix/$PROJ_CODE/fastq/
fi

echo Copying fastq files from $OUT_DIR/$SEQ_CODE/ to $PSEQ_obelix/$PROJ_CODE/fastq/
cp $OUT_DIR/$SEQ_CODE/*.fastq.gz $PSEQ_obelix/$PROJ_CODE/fastq/
