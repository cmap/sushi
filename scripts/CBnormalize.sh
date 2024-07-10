#!/bin/bash

echo Starting fastq2readcounts...

if [ -z "$BUILD_DIR" ]
then
	echo BUILD_DIR not specified
    exit -1
fi

if [ -z "$READ_COUNTS" ]
then
	echo READ_COUNTS parameter empty
    exit -1
fi

if [ -z "$SAMPLE_META" ]
then
	echo SAMPLE_META parameter empty
    exit -1
fi



#Enforces abs paths
if [[ "$SAMPLE_META" = /* ]]
then
	SAMPLE_META=$(ls $SAMPLE_META)
else
	SAMPLE_META=$BUILD_DIR/$SAMPLE_META
fi

#Enforces abs paths
if [[ "$READ_COUNTS" = /* ]]
then
	READ_COUNTS=$(ls $READ_COUNTS)
else
	READ_COUNTS=$BUILD_DIR/$READ_COUNTS
fi


echo Build dir is: $BUILD_DIR
echo SAMPLE_META is: $SAMPLE_META


echo $RUN_NORM


if [[ "$RUN_NORM" == "true" ]]
then
	echo "Running module"

	echo Rscript CBnormalize.R -c $READ_COUNTS	\
    --CB_meta $CONTROL_BARCODE_META \
    --pseudocount $PSEUDOCOUNT \
    --out $BUILD_DIR

	Rscript CBnormalize.R -c $READ_COUNTS	\
	--CB_meta $CONTROL_BARCODE_META \
    --pseudocount $PSEUDOCOUNT \
	--out $BUILD_DIR

	COUNTS="normalized_counts.csv"
#    COUNT_COL_NAME="normalized_n"
else
	echo "Not running module"
    COUNTS=$READ_COUNTS
    COUNT_COL_NAME="n"
    echo $COUNTS
fi

#curl "http://suffix:8889/view/sushi/job/filteredCounts_QC/buildWithParameters?token=cIJQNoc7aL&BUILD_DIR=${BUILD_DIR}&CELL_SET_META=${CELL_SET_META}&CONTROL_BARCODE_META=${CONTROL_BARCODE_META}&COUNT_THRESHOLD=${COUNT_THRESHOLD}&COUNT_COL_NAME=${COUNT_COL_NAME}&SIG_COLS=${SIG_COLS}"
#curl "http://vercingetorix-r8:8889/view/sushi-podman/job/compute_LFC-podman/buildWithParameters?token=jjrDRMJ8Yh&BUILD_DIR=${BUILD_DIR}&BUILD_NAME=${BUILD_NAME}&CTL_TYPES=${CTL_TYPES}&COUNTS=${COUNTS}&COUNT_COL_NAME=${COUNT_COL_NAME}&SIG_COLS=${SIG_COLS}&SAMPLE_COLS=${SAMPLE_COLS}&CONTROL_COLS=${CONTROL_COLS}&COUNT_THRESHOLD=${COUNT_THRESHOLD}&DAYS=${DAYS}&CONVERT_SUSHI=${CONVERT_SUSHI}"

#Rscript -e "library(cmapR);packageVersion('cmapR')"

