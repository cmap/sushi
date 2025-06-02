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
if [[ "$NORMALIZED_COUNTS" = /* ]]
then
	NORMALIZED_COUNTS=$(ls $NORMALIZED_COUNTS)
else
	NORMALIZED_COUNTS=$BUILD_DIR/$NORMALIZED_COUNTS
fi

#Enforces abs paths
if [[ "$PRISM_BARCODE_COUNTS" = /* ]]
then
	PRISM_BARCODE_COUNTS=$(ls $PRISM_BARCODE_COUNTS)
else
	PRISM_BARCODE_COUNTS=$BUILD_DIR/$PRISM_BARCODE_COUNTS
fi

#Enforces abs paths
if [[ "$UNKNOWN_BARCODE_COUNTS" = /* ]]
then
	UNKNOWN_BARCODE_COUNTS=$(ls $UNKNOWN_BARCODE_COUNTS)
else
	UNKNOWN_BARCODE_COUNTS=$BUILD_DIR/$UNKNOWN_BARCODE_COUNTS
fi

#Enforces abs paths
if [[ "$RAW_COUNTS_UNCOLLAPSED" = /* ]]
then
	RAW_COUNTS_UNCOLLAPSED=$(ls $RAW_COUNTS_UNCOLLAPSED)
else
	RAW_COUNTS_UNCOLLAPSED=$BUILD_DIR/$RAW_COUNTS_UNCOLLAPSED
fi

#Enforces abs paths
if [[ "$LFC" = /* ]]
then
	LFC=$(ls $LFC)
else
	LFC=$BUILD_DIR/$LFC
fi

#Enforces abs paths
if [[ "$CONTROL_BARCODE_META" = /* ]]
then
	CONTROL_BARCODE_META=$(ls $CONTROL_BARCODE_META)
else
	CONTROL_BARCODE_META=$BUILD_DIR/$CONTROL_BARCODE_META
fi

echo Build dir is: $BUILD_DIR
echo SAMPLE_META is: $SAMPLE_META
echo ANNOTATED_COUNTS is: $ANNOTATED_COUNTS
echo NORMALIZED_COUNTS is: $NORMALIZED_COUNTS
echo CONTROL_BARCODE_META is: $CONTROL_BARCODE_META
echo COUNT_THRESHOLD is: $COUNT_THRESHOLD
echo RAW_COUNTS_UNCOLLAPSED is: $RAW_COUNTS_UNCOLLAPSED
echo LFC is: $LFC
echo REVERSE_INDEX2 is: $REVERSE_INDEX2
echo PSEUDOCOUNT is: $PSEUDOCOUNT

args=(
--sample_meta "$SAMPLE_META"
--annotated_counts "$ANNOTATED_COUNTS"
--normalized_counts "$NORMALIZED_COUNTS"
--sig_cols "$SIG_COLS"
--out "$BUILD_DIR"
--count_threshold "$COUNT_THRESHOLD"
--control_type "$CTL_TYPES"
--raw_counts_uncollapsed "$RAW_COUNTS_UNCOLLAPSED"
--prism_barcode_counts "$PRISM_BARCODE_COUNTS"
--unknown_barcode_counts "$UNKNOWN_BARCODE_COUNTS"
--lfc "$LFC"
--id_cols "$ID_COLS"
--barcode_col "$BARCODE_COL"
--reverse_index2 "$REVERSE_INDEX2"
--chunk_size "$CHUNK_SIZE"
--cell_line_cols "$CELL_LINE_COLS"
--pseudocount "$PSEUDOCOUNT"
)

echo Rscript filter_counts_qc/filter_counts_qc.R "${args[@]}"

Rscript filter_counts_qc/filter_counts_qc.R "${args[@]}"