#!/bin/bash

echo Starting replicate collapse...

if [ -z "$BUILD_DIR" ]
then
	echo BUILD_DIR not specified
    exit -1
fi

if [ -z "$COLLAPSED_VALUES" ]
then
	echo COLLAPSED_VALUES parameter empty
    exit -1
fi


#Enforces abs paths
if [[ "$COLLAPSED_VALUES" = /* ]]
then
	COLLAPSED_VALUES=$(ls $COLLAPSED_VALUES)
else
	COLLAPSED_VALUES=$BUILD_DIR/$COLLAPSED_VALUES
fi


echo Build dir is: $BUILD_DIR
echo COLLAPSED_VALUES is: $COLLAPSED_VALUES

echo Rscript collapse_replicates.R -c $COLLAPSED_VALUES	\
--out $BUILD_DIR \
--sig_cols $SIG_COLS


Rscript collapse_replicates.R -c $COLLAPSED_VALUES	\
--out $BUILD_DIR \
--sig_cols $SIG_COLS
