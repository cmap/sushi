#!/bin/bash

echo Starting collate_fastq_reads...

export API_KEY=$(cat /local/jenkins/.clue_api_key)
export API_URL="https://api.clue.io/api/"

if [ -z "$BUILD_DIR" ]
then
	echo BUILD_DIR not specified
    exit -1
fi

case "$SEQ_TYPE" in
  HiSeq|NovaSeq)
    INDEX_1="barcode_1"
    INDEX_2="barcode_2"
    BARCODE_SUFFIX="unmapped.1."
    ;;
  MiSeq)
    INDEX_1="_I1_"
    INDEX_2="_I2_"
    BARCODE_SUFFIX="_R1_"
    ;;
  DRAGEN)
    INDEX_1="_I1_"
    INDEX_2="_I2_"
    BARCODE_SUFFIX="_R1_"
    ;;
  *)
	echo "Unknown SEQ_TYPE"
    exit -1
    ;;
esac

if [[ "$SEQ_TYPE" == "NovaSeq" ]]
then
	export REVERSE_INDEX2=TRUE
    echo NovaSeq
fi

if [[ "$SEQ_TYPE" == "DRAGEN" ]]
then
	export REVERSE_INDEX2=TRUE
    echo DRAGEN
fi

if [ -z "$INDEX_1" ]
then
	echo INDEX_1 parameter empty
    exit -1
fi

if [ -z "$INDEX_2" ]
then
	echo INDEX_2 parameter empty
    exit -1
fi

if [ -z "$BARCODE_SUFFIX" ]
then
	echo BARCODE_SUFFIX parameter empty
    exit -1
fi

#Enforces abs paths
if [[ "$RAW_COUNTS_UNCOLLAPSED" = /* ]]
then
	RAW_COUNTS_UNCOLLAPSED=$(ls $RAW_COUNTS_UNCOLLAPSED)
else
	RAW_COUNTS_UNCOLLAPSED=$BUILD_DIR/$RAW_COUNTS_UNCOLLAPSED
fi

#Enforces abs paths
if [[ "$SAMPLE_META" = /* ]]
then
  SAMPLE_META=$(ls $SAMPLE_META)
else
  SAMPLE_META=$BUILD_DIR/$SAMPLE_META
fi

echo Build dir is: $BUILD_DIR

PROJECT_DIR=$(dirname "$BUILD_DIR")
PROJECT_CODE=$(basename "$PROJECT_DIR")

echo Project Code: $PROJECT_CODE
echo REVERSE_INDEX2 is: $REVERSE_INDEX2

args=(
--raw_counts_uncollapsed "$RAW_COUNTS_UNCOLLAPSED"
--sample_meta "$SAMPLE_META"
--out "$BUILD_DIR"
--sequencing_index_cols="$SEQUENCING_INDEX_COLS"
--id_cols "$ID_COLS" 
--reverse_index2 "$REVERSE_INDEX2"
)

echo Rscript collate_fastq_reads.R "${args[@]}"
Rscript collate_fastq_reads.R "${args[@]}"
