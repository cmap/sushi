#!/bin/bash

echo Starting compute_l2fc...

if [ -z "$BUILD_DIR" ]
then
	echo BUILD_DIR not specified
    exit -1
fi

if [ -z "$LFC" ]
then
	echo LFC parameter empty
    exit -1
fi

#Enforces abs paths
if [[ "$LFC" = /* ]]
then
	LFC=$(ls $LFC)
else
	LFC=$LFC/$LFC
fi

echo Build dir is: $BUILD_DIR
echo LFC is: $LFC

args=(
--l2fc "$LFC"
--build_dir "$BUILD_DIR"
--sig_cols "$SIG_COLS"
--cell_line_cols "$CELL_LINE_COLS"
--l2fc_column "$L2FC_COLUMN"
)

echo Rscript dose_response.R "${args[@]}"
Rscript dose_response.R "${args[@]}"
