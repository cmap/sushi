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
if [[ "$RAW_COUNTS" = /* ]]
then
	RAW_COUNTS=$(ls $RAW_COUNTS)
else
	RAW_COUNTS=$BUILD_DIR/$RAW_COUNTS
fi

echo $CELL_LINE_META

#Enforces abs paths
if [[ "$CELL_LINE_META" = /* ]]
then
	CELL_LINE_META=$(ls $CELL_LINE_META)
else
	CELL_LINE_META=/data/$CELL_LINE_META
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
	CELL_SET_META=$(ls $CONTROL_BARCODE_META)
else
	CELL_SET_META=/data/$CONTROL_BARCODE_META
fi


echo Build dir is: $BUILD_DIR
echo SAMPLE_META is: $SAMPLE_META
echo RAW_COUNTS is: $RAW_COUNTS
echo CELL_LINE_META is: $CELL_LINE_META
echo CONTROL_BARCODE_META is: $CONTROL_BARCODE_META
echo CELL_SET_META is: $CELL_SET_META
echo REVERSE_INDEX2 is: $REVERSE_INDEX2

args=(
-c "$RAW_COUNTS"
--sample_meta "$SAMPLE_META"
--id_cols "$ID_COLS"
--cell_line_meta "$CELL_LINE_META"
--CB_meta "$CONTROL_BARCODE_META"
--cell_set_meta "$CELL_SET_META"
--out "$BUILD_DIR"
--count_threshold "$COUNT_THRESHOLD"
--sequencing_index_cols "$SEQUENCING_INDEX_COLS"
)

if [[ "$REVERSE_INDEX2" == "true" ]]
then
	args+=(--reverse_index2)
fi

echo Checking whether pool ID flag was passed...

if [[ "$PULL_POOL_ID" == "true" ]]; then
    args+=(--pool_id)
    ASSAY_POOL_META="$BUILD_DIR/assay_pool_meta.csv"
    args+=(--assay_pool_meta "$ASSAY_POOL_META")
fi


if [[ "$REMOVE_DATA" == "true" ]]
then
	args+=(--rm_data)
fi

echo Rscript filter_counts.R "${args[@]}"

# Define an array of parameters and their corresponding values
parameters=(
  "BUILD_DIR:$BUILD_DIR"
  "SAMPLE_META:$SAMPLE_META"
  "SEQ_TYPE:$SEQ_TYPE"
  "ID_COLS:$ID_COLS"
  "SIG_COLS:$SIG_COLS"
  "CONTROL_COLS:$CONTROL_COLS"
  "COUNT_THRESHOLD:$COUNT_THRESHOLD"
  "PSEUDOCOUNT:$PSEUDOCOUNT"
  "COUNT_COL_NAME:$COUNT_COL_NAME"
  "CELL_SET_META:$CELL_SET_META"
  "CONTROL_BARCODE_META:$CONTROL_BARCODE_META"
  "CELL_LINE_META:$CELL_LINE_META"
  "REVERSE_INDEX2:$REVERSE_INDEX2"
  "REMOVE_DATA:$REMOVE_DATA"
  "PULL_POOL_ID:$PULL_POOL_ID"
  "BUILD_NAME:$BUILD_NAME"
  "SEQUENCING_INDEX_COLS:$SEQUENCING_INDEX_COLS"
)

# Overwrite the CSV file with parameters each time a new build is kicked off
> $BUILD_DIR/config_from_filtercounts.csv

for param in "${parameters[@]}"; do
  # Split the parameter and value using ':'
  IFS=':' read -ra parts <<< "$param"
  parameter="${parts[0]}"
  value="${parts[1]}"
  echo "$parameter: $value"
  echo "$parameter,$value" >> $BUILD_DIR/config_from_filtercounts.csv
done


Rscript -e "library(cmapR);packageVersion('cmapR')"


Rscript filter_counts.R "${args[@]}"