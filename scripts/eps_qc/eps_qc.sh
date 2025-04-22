#!/bin/bash

echo Starting EPS_QC...

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

# Define an array of parameters and their corresponding values
parameters=(
  "BUILD_DIR:$BUILD_DIR"
  "COUNT_THRESHOLD:$COUNT_THRESHOLD"
  "BUILD_NAME:$BUILD_NAME"
  "CONTROL_TYPES:$CTL_TYPES"
  "DAYS:$DAYS"
  "PSEUDOCOUNT:$PSEUDOCOUNT"
)

args=(
--out "$BUILD_DIR"
--build_dir "$BUILD_DIR"
--control_type "$CTL_TYPES"
--count_threshold "$COUNT_THRESHOLD"
--pseudocount "$PSEUDOCOUNT"
--days "$DAYS"
--name "$BUILD_NAME"
)

echo Rscript EPS_QC.R "${args[@]}"

Rscript EPS_QC.R "${args[@]}"
