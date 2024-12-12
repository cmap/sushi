#!/bin/bash

echo Starting biomarker analysis...

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

enforce_abs_path_biomarker() {
  local var_name=$1
  local value=${!var_name}

  if [[ "$value" = /* ]]; then
    eval "$var_name=$(ls "$value")"
  else
    eval "$var_name=/cmap_obelix/pod/prismSeq/data/biomarker/current/$value"
  fi
  echo "$var_name is: ${!var_name}"
}

enforce_abs_path COLLAPSED_LFC
enforce_abs_path_biomarker BIOMARKER_FILE

args=(
--collapsed_lfc "$COLLAPSED_LFC"
--out "$BUILD_DIR"
--sig_cols "$SIG_COLS"
--univariate_biomarker "$UNIVARIATE_BIOMARKER"
--multivariate_biomarker "$MULTIVARIATE_BIOMARKER"
--biomarker_file "$BIOMARKER_FILE"
)

echo Rscript biomarker.R "${args[@]}"
Rscript biomarker.R "${args[@]}"
