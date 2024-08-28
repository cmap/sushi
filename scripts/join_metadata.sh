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

if [ -z "$COLLAPSED_L2FC" ]
then
	echo Collapsed l2fc parameter empty
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
if [[ "$COLLAPSED_L2FC" = /* ]]
then
	COLLAPSED_L2FC=$(ls $COLLAPSED_L2FC)
else
	COLLAPSED_L2FC=$BUILD_DIR/$COLLAPSED_L2FC
fi

echo Build dir is: $BUILD_DIR
echo LFC is: $LFC
echo COLLAPSED_L2FC is: $COLLAPSED_L2FC

echo Rscript join_metadata.R -c $LFC	\
--collapsed_l2fc $COLLAPSED_L2FC \
--out $BUILD_DIR \
--sig_cols $SIG_COLS

Rscript join_metadata.R -c $LFC	\
--collapsed_l2fc $COLLAPSED_L2FC \
--out $BUILD_DIR \
--sig_cols $SIG_COLS