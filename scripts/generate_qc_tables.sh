#!/bin/bash

echo Starting generate_qc_tables...

if [ -z "$BUILD_DIR" ]
then
	echo BUILD_DIR not specified
    exit -1
fi
echo Build dir is: $BUILD_DIR

#Enforces abs paths

enforce_abs_path() {
  local var_name=$1
  local value=${!var_name}

  if [[ "$value" = /* ]]; then
    eval "$var_name=$(ls "$value")"
  else
    eval "$var_name=$BUILD_DIR/$value"
  fi
  echo "$var_name is: ${!var_name}"
}

enforce_abs_path NORMALIZED_COUNTS
enforce_abs_path ANNOTATED_COUNTS
enforce_abs_path CELL_SET_AND_POOL_META
enforce_abs_path FILTERED_COUNTS
enforce_abs_path CONTROL_BARCODE_META
enforce_abs_path UNKNOWN_BARCODE_COUNTS

args=(
--normalized_counts "$NORMALIZED_COUNTS"
--annotated_counts "$ANNOTATED_COUNTS"
--filtered_counts "$FILTERED_COUNTS"
--negcon_type "$CTL_TYPES"
--poscon_type "$POSCON_TYPE"
--cell_set_and_pool_meta "$CELL_SET_AND_POOL_META"
--out "$BUILD_DIR"
--id_cols "$ID_COLS"
--cell_line_cols "$CELL_LINE_COLS"
--count_threshold "$COUNT_THRESHOLD"
--control_barcode_meta "$CONTROL_BARCODE_META"
--unknown_barcode_counts "$UNKNOWN_BARCODE_COUNTS"
)

echo Rscript qc_tables.R "${args[@]}"
Rscript qc_tables.R "${args[@]}"
