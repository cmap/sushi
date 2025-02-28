#!/bin/bash

echo Starting dose response...

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

enforce_abs_path LFC

echo Build dir is: $BUILD_DIR
echo LFC is: $LFC

args=(
--l2fc "$LFC"
--build_dir "$BUILD_DIR"
--sig_cols "$SIG_COLS"
--cell_line_cols "$CELL_LINE_COLS"
--l2fc_column "$L2FC_COLUMN"
--cap_for_viability "$VIABILITY_CAP"
--out_dir "$BUILD_DIR/drc"
)

echo Rscript dose_response.R "${args[@]}"
Rscript dose_response.R "${args[@]}"
