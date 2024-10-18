#!/bin/bash

echo Starting compute_l2fc...

if [ -z "$BUILD_DIR" ]
then
	echo BUILD_DIR not specified
    exit -1
fi

if [ -z "$NORMALIZED_COUNTS" ]
then
	echo NORMALIZED_COUNTS parameter empty
    exit -1
fi

#Enforces abs paths
if [[ "$NORMALIZED_COUNTS" = /* ]]
then
	NORMALIZED_COUNTS=$(ls $NORMALIZED_COUNTS)
else
	NORMALIZED_COUNTS=$BUILD_DIR/$NORMALIZED_COUNTS
fi

echo Build dir is: $BUILD_DIR
echo NORMALIZED_COUNTS is: $NORMALIZED_COUNTS

args=(
--normalized_counts "$NORMALIZED_COUNTS"
--out "$BUILD_DIR"
--control_type $CTL_TYPES 
--count_col_name $COUNT_COL_NAME 
--sig_cols $SIG_COLS 
--ctrl_cols $CONTROL_COLS 
--count_threshold $COUNT_THRESHOLD 
--cell_line_cols $CELL_LINE_COLS
)

echo Rscript compute_l2fc.R "${args[@]}"
Rscript compute_l2fc.R "${args[@]}"
