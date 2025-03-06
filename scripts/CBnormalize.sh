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
	CONTROL_BARCODE_META=$BUILD_DIR/$CONTROL_BARCODE_META
fi

#Enforces abs paths
if [[ "$UNKNOWN_BARCODE_COUNTS" = /* ]]
then
  UNKNOWN_BARCODE_COUNTS=$(ls $UNKNOWN_BARCODE_COUNTS)
else
  UNKNOWN_BARCODE_COUNTS=$BUILD_DIR/$UNKNOWN_BARCODE_COUNTS
fi

#Enforces abs paths
if [[ "$ANNOTATED_COUNTS" = /* ]]
then
  ANNOTATED_COUNTS=$(ls $ANNOTATED_COUNTS)
else
  ANNOTATED_COUNTS=$BUILD_DIR/$ANNOTATED_COUNTS
fi

#Enforces abs paths
if [[ "$CELL_SET_AND_POOL_META" = /* ]]
then
  CELL_SET_AND_POOL_META=$(ls $CELL_SET_AND_POOL_META)
else
  CELL_SET_AND_POOL_META=$BUILD_DIR/$CELL_SET_AND_POOL_META
fi


echo Build dir is: $BUILD_DIR
echo SAMPLE_META is: $SAMPLE_META


echo $RUN_NORM


echo "Running normalization module"

echo Rscript CBnormalize.R -c $FILTERED_COUNTS	\
--CB_meta $CONTROL_BARCODE_META \
--pseudocount $PSEUDOCOUNT \
--id_cols $ID_COLS \
--out $BUILD_DIR \
--filtered_counts $FILTERED_COUNTS \
--annotated_counts $ANNOTATED_COUNTS \
--unknown_barcode_counts $UNKNOWN_BARCODE_COUNTS \
--cell_set_meta $CELL_SET_AND_POOL_META

Rscript CBnormalize.R -c $FILTERED_COUNTS	\
--CB_meta $CONTROL_BARCODE_META \
--pseudocount $PSEUDOCOUNT \
--id_cols $ID_COLS \
--out $BUILD_DIR \
--filtered_counts $FILTERED_COUNTS \
--annotated_counts $ANNOTATED_COUNTS \
--unknown_barcode_counts $UNKNOWN_BARCODE_COUNTS \
--cell_set_meta $CELL_SET_AND_POOL_META

