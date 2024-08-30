#!/bin/bash

echo Starting metadata join...

if [ -z "$BUILD_DIR" ]
then
	echo BUILD_DIR not specified
    exit -1
fi

if [ -z "$L2FC" ]
then
	echo LFC parameter empty
    exit -1

fi

if [ -z "$COLLAPSED_VALUES" ]
then
	echo Collapsed l2fc parameter empty
    exit -1
fi

if [ -z "$ASSAY_POOL_META" ]
then
	echo ASSAY_POOL_META parameter empty
    exit -1
fi

#Enforces abs paths
if [[ "$LFC" = /* ]]
then
	LFC=$(ls $LFC)
else
	LFC=$BUILD_DIR/$LFC
fi

#Enforces abs paths
if [[ "$COLLAPSED_VALUES" = /* ]]
then
	COLLAPSED_VALUES=$(ls $COLLAPSED_VALUES)
else
	COLLAPSED_VALUES=$BUILD_DIR/$COLLAPSED_VALUES
fi

#Enforces abs paths
if [[ "$ASSAY_POOL_META" = /* ]]
then
	ASSAY_POOL_META=$(ls $ASSAY_POOL_META)
else
	ASSAY_POOL_META=$BUILD_DIR/$ASSAY_POOL_META
fi

#Enforces abs paths
if [[ "$SAMPLE_META" = /* ]]
then
	SAMPLE_META=$(ls $SAMPLE_META)
else
	SAMPLE_META=$BUILD_DIR/$SAMPLE_META
fi

echo Build dir is: $BUILD_DIR
echo LFC is: $LFC
echo COLLAPSED_VALUES is: $COLLAPSED_VALUES
echo SAMPLE_META is: $SAMPLE_META

echo Rscript join_metadata.R -c $LFC	\
--collapsed_l2fc $COLLAPSED_VALUES \
--assay_pool_meta $ASSAY_POOL_META \
--out $BUILD_DIR \
--sig_cols $SIG_COLS \
--sample_meta $SAMPLE_META

Rscript join_metadata.R -c $LFC	\
--collapsed_l2fc $COLLAPSED_VALUES \
--assay_pool_meta $ASSAY_POOL_META \
--out $BUILD_DIR \
--sig_cols $SIG_COLS \
--sample_meta $SAMPLE_META