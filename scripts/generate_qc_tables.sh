#!/bin/bash

echo Starting compute_l2fc...

if [ -z "$BUILD_DIR" ]
then
	echo BUILD_DIR not specified
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

if [[ "$ANNOTATED_COUNTS" = /* ]]
then
  ANNOTATED_COUNTS=$(ls $ANNOTATED_COUNTS)
else
  ANNOTATED_COUNTS=$BUILD_DIR/$ANNOTATED_COUNTS
fi

args=(
--normalized_counts "$NORMALIZED_COUNTS"
--annotated_counts "$ANNOTATED_COUNTS"
--negcon_type "$CTL_TYPES"
--poscon_type "$POSCON_TYPE"
--cell_set_and_pool_meta "$CELL_SET_AND_POOL_META"
--out "$BUILD_DIR"
)

echo Rscript compute_l2fc.R "${args[@]}"
Rscript compute_l2fc.R "${args[@]}"
