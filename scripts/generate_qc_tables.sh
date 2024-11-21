#!/bin/bash

echo Starting generate_qc_tables...

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
echo ANNOTATED_COUNTS is: $ANNOTATED_COUNTS

if [[ "$CELL_SET_AND_POOL_META" = /* ]]
then
  CELL_SET_AND_POOL_META=$(ls $CELL_SET_AND_POOL_META)
else
  CELL_SET_AND_POOL_META=$BUILD_DIR/$CELL_SET_AND_POOL_META
fi
echo CELL_SET_AND_POOL_META is: $CELL_SET_AND_POOL_META

args=(
--normalized_counts "$NORMALIZED_COUNTS"
--annotated_counts "$ANNOTATED_COUNTS"
--negcon_type "$CTL_TYPES"
--poscon_type "$POSCON_TYPE"
--cell_set_and_pool_meta "$CELL_SET_AND_POOL_META"
--out "$BUILD_DIR"
)

echo Rscript generate_qc_tables.R "${args[@]}"
Rscript qc.R "${args[@]}"