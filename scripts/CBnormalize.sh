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

#Enforces abs paths
if [[ "$CONTROL_BARCODE_META" = /* ]]
then
	CONTROL_BARCODE_META=$(ls $CONTROL_BARCODE_META)
else
	CONTROL_BARCODE_META=/data/$CONTROL_BARCODE_META
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
