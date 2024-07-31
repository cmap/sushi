#!/bin/bash

echo Starting filteredCounts_QC...

if [ -z "$BUILD_DIR" ]
then
	echo BUILD_DIR not specified
    exit -1
fi


if [ -z "$ANNOTATED_COUNTS" ]
then
	echo ANNOTATED_COUNTS parameter empty
    exit -1
fi

if [ -z "$NORMALIZED_COUNTS" ]
then
	echo NORMALIZED_COUNTS parameter empty
    exit -1
fi

if [ -z "$SAMPLE_META" ]
then
	echo SAMPLE_META parameter empty
    exit -1
fi

if [ -z "$CELL_SET_META" ]
then
	echo CELL_SET_META parameter empty
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
if [[ "$ANNOTATED_COUNTS" = /* ]]
then
	ANNOTATED_COUNTS=$(ls $ANNOTATED_COUNTS)
else
	ANNOTATED_COUNTS=$BUILD_DIR/$ANNOTATED_COUNTS
fi

#Enforces abs paths
if [[ "$CELL_SET_META" = /* ]]
then
	CELL_SET_META=$(ls $CELL_SET_META)
else
	CELL_SET_META=$BUILD_DIR/$CELL_SET_META
fi

#Enforces abs paths
if [[ "$NORMALIZED_COUNTS" = /* ]]
then
	COUNTS=$(ls $NORMALIZED_COUNTS)
else
	COUNTS=$BUILD_DIR/$NORMALIZED_COUNTS
fi

#Enforces abs paths
if [[ "$RAW_COUNTS" = /* ]]
then
	RAW_COUNTS=$(ls $RAW_COUNTS)
else
	RAW_COUNTS=$BUILD_DIR/$RAW_COUNTS
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
echo ANNOTATED_COUNTS is: $ANNOTATED_COUNTS
echo NORMALIZED_COUNTS is: $NORMALIZED_COUNTS
echo CELL_SET_META is: $CELL_SET_META
echo CONTROL_BARCODE_META is: $CONTROL_BARCODE_META
echo COUNT_THRESHOLD is: $COUNT_THRESHOLD
echo COUNT_COL_NAME is: $COUNT_COL_NAME
echo RAW_COUNTS is: $RAW_COUNTS

args=(
--sample_meta "$SAMPLE_META"
--annotated_counts "$ANNOTATED_COUNTS"
--normalized_counts "$NORMALIZED_COUNTS"
--sig_cols "$SIG_COLS"
--cell_set_meta "$CELL_SET_META"
--CB_meta "$CONTROL_BARCODE_META"
--out "$BUILD_DIR"
--count_threshold "$COUNT_THRESHOLD"
--count_col_name "$COUNT_COL_NAME"
--control_type "$CTL_TYPES"
--raw_counts "$RAW_COUNTS"
)

if [[ "$REVERSE_INDEX2" == "true" ]]
then
	args+=(--reverse_index2)
fi

echo Rscript filteredCounts_QC.R --sample_meta $SAMPLE_META \
--annotated_counts $ANNOTATED_COUNTS \
--normalized_counts $NORMALIZED_COUNTS \
--cell_set_meta $CELL_SET_META \
--CB_meta $CONTROL_BARCODE_META \
--sig_cols $SIG_COLS \
--out $BUILD_DIR \
--count_threshold $COUNT_THRESHOLD \
--count_col_name $COUNT_COL_NAME \
--reverse_index2 $REVERSE_INDEX2 \
--control_type $CTL_TYPES \
--raw_counts $RAW_COUNTS


Rscript -e "library(cmapR);packageVersion('cmapR')"

Rscript filteredCounts_QC.R "${args[@]}"