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
  elif [[ -d "$BUILD_DIR/$value" ]]; then
    eval "$var_name=$BUILD_DIR/$value"
  elif [[ -d "$BUILD_DIR/drc/$value" ]]; then
    eval "$var_name=$BUILD_DIR/drc/$value"
  elif [[ -d "$BUILD_DIR/biomarker/$value" ]]; then
    eval "$var_name=$BUILD_DIR/biomarker/$value"
  else
    echo "Error: Path $value does not exist."
    return 1
  fi

  echo "$var_name is: ${!var_name}"
}


enforce_abs_path COLLAPSED_LFC
enforce_abs_path DRC_PATH

args=(
--collapsed_lfc "$COLLAPSED_LFC"
--auc_path "$DRC_PATH"
--build_dir "$BUILD_DIR"
--sig_cols "$SIG_COLS"
--univariate_biomarker "$UNIVARIATE_BIOMARKER"
--multivariate_biomarker "$MULTIVARIATE_BIOMARKER"
--biomarker_file "$BIOMARKER_FILE"
--auc_column "$AUC_COLUMN"
--lfc_biomarker "$LFC_BIOMARKER"
--auc_biomarker "$AUC_BIOMARKER"
)

echo Rscript biomarker.R "${args[@]}"
Rscript biomarker.R "${args[@]}"
