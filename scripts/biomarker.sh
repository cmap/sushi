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

enforce_abs_path COLLAPSED_LFC

args=(
--collapsed_lfc "$COLLAPSED_LFC"
--build_dir "$BUILD_DIR"
--sig_cols "$SIG_COLS"
--univariate_biomarker "$UNIVARIATE_BIOMARKER"
--multivariate_biomarker "$MULTIVARIATE_BIOMARKER"
--biomarker_file "$BIOMARKER_FILE"
--auc_column "$AUC_COLUMN"
)

echo Rscript biomarker.R "${args[@]}"
Rscript biomarker.R "${args[@]}"
