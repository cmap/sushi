#!/bin/bash

echo Starting CBnormalize...

if [ -z "$BUILD_DIR" ]
then
	echo BUILD_DIR not specified
    exit -1
fi

if [ -z "$FILTERED_COUNTS" ]
then
	echo FILTERED_COUNTS parameter empty
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
if [[ "$FILTERED_COUNTS" = /* ]]
then
	FILTERED_COUNTS=$(ls $FILTERED_COUNTS)
else
	FILTERED_COUNTS=$BUILD_DIR/$FILTERED_COUNTS
fi


echo Build dir is: $BUILD_DIR
echo SAMPLE_META is: $SAMPLE_META


echo $RUN_NORM


if [[ "$RUN_NORM" == "true" ]]
then
	echo "Running module"

	echo Rscript CBnormalize.R -c $FILTERED_COUNTS	\
    --CB_meta $CONTROL_BARCODE_META \
    --pseudocount $PSEUDOCOUNT \
    --out $BUILD_DIR

	Rscript CBnormalize.R -c $FILTERED_COUNTS	\
	--CB_meta $CONTROL_BARCODE_META \
    --pseudocount $PSEUDOCOUNT \
	--out $BUILD_DIR

	COUNTS="normalized_counts.csv"

else
	echo "Not running module"
    COUNTS=$FILTERED_COUNTS
    COUNT_COL_NAME="n"
    echo $COUNTS
fi

#curl "http://suffix:8889/view/sushi/job/filteredCounts_QC/buildWithParameters?token=cIJQNoc7aL&BUILD_DIR=${BUILD_DIR}&CELL_SET_META=${CELL_SET_META}&CONTROL_BARCODE_META=${CONTROL_BARCODE_META}&COUNT_THRESHOLD=${COUNT_THRESHOLD}&COUNT_COL_NAME=${COUNT_COL_NAME}&SIG_COLS=${SIG_COLS}"
#curl "http://vercingetorix-r8:8889/view/sushi-podman/job/compute_LFC-podman/buildWithParameters?token=jjrDRMJ8Yh&BUILD_DIR=${BUILD_DIR}&BUILD_NAME=${BUILD_NAME}&CTL_TYPES=${CTL_TYPES}&COUNTS=${COUNTS}&COUNT_COL_NAME=${COUNT_COL_NAME}&SIG_COLS=${SIG_COLS}&SAMPLE_COLS=${SAMPLE_COLS}&CONTROL_COLS=${CONTROL_COLS}&COUNT_THRESHOLD=${COUNT_THRESHOLD}&DAYS=${DAYS}&CONVERT_SUSHI=${CONVERT_SUSHI}"

#Rscript -e "library(cmapR);packageVersion('cmapR')"
