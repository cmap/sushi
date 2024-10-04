#!/bin/bash

echo Starting filter_counts...

if [ -z "$BUILD_DIR" ]
then
	echo BUILD_DIR not specified
    exit -1
fi

if [ -z "$RAW_COUNTS" ]
then
	echo RAW_COUNTS parameter empty
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
if [[ "$PRISM_BARCODE_COUNTS" = /* ]]
then
	PRISM_BARCODE_COUNTS=$(ls $PRISM_BARCODE_COUNTS)
else
	PRISM_BARCODE_COUNTS=$BUILD_DIR/$PRISM_BARCODE_COUNTS
fi

echo $CELL_LINE_META

#Enforces abs paths
if [[ "$CELL_LINE_META" = /* ]]
then
	CELL_LINE_META=$(ls $CELL_LINE_META)
else
	CELL_LINE_META=$BUILD_DIR/$CELL_LINE_META
fi

echo $CELL_SET_META

#Enforces abs paths
if [[ "$CELL_SET_META" = /* ]]
then
	CELL_SET_META=$(ls $CELL_SET_META)
else
	CELL_SET_META=$BUILD_DIR/$CELL_SET_META
fi

#Enforces abs paths
if [[ "$CONTROL_BARCODE_META" = /* ]]
then
  CONTROL_BARCODE_META=$(ls $CONTROL_BARCODE_META)
else
	CONTROL_BARCODE_META=$BUILD_DIR/$CONTROL_BARCODE_META
fi

#Enforces abs paths
if [[ "$ASSAY_POOL_META" = /* ]]
then
	ASSAY_POOL_META=$(ls $ASSAY_POOL_META)
else
	ASSAY_POOL_META=$BUILD_DIR/$ASSAY_POOL_META
fi

echo Build dir is: $BUILD_DIR
echo SAMPLE_META is: $SAMPLE_META
echo PRISM_BARCODE_COUNTS is: $PRISM_BARCODE_COUNTS
echo CELL_LINE_META is: $CELL_LINE_META
echo CONTROL_BARCODE_META is: $CONTROL_BARCODE_META
echo CELL_SET_META is: $CELL_SET_META
echo ID_COLS is: $ID_COLS

args=(
--prism_barcode_counts "$PRISM_BARCODE_COUNTS"
--sample_meta "$SAMPLE_META"
--cell_line_meta "$CELL_LINE_META"
--CB_meta "$CONTROL_BARCODE_META"
--cell_set_meta "$CELL_SET_META"
--id_cols "$ID_COLS"
--out "$BUILD_DIR"
--pool_id "$PULL_POOL_ID"
--rm_data "$REMOVE_DATA"
--assay_pool_meta "$ASSAY_POOL_META"
)

echo Rscript filter_counts.R "${args[@]}"

Rscript filter_counts.R "${args[@]}"