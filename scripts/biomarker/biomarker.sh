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

  echo "DEBUG: Resolving path for $var_name (value: '$value')"

  # Check if value is already an absolute path
  if [[ "$value" = /* ]]; then
    echo "DEBUG: Absolute path detected: $value"
    if [[ -e "$value" ]]; then
      eval "$var_name=\"$value\""
      echo "DEBUG: $var_name resolved to absolute path: ${!var_name}"
      return 0
    else
      echo "ERROR: Absolute path '$value' does not exist."
      return 1
    fi
  fi

  # Try appending to possible base directories
  for base_path in "$BUILD_DIR" "$BUILD_DIR/drc" "$BUILD_DIR/biomarker"; do
    local full_path="$base_path/$value"
    echo "DEBUG: Checking relative path: $full_path"
    if [[ -e "$full_path" ]]; then
      eval "$var_name=\"$full_path\""
      echo "DEBUG: $var_name resolved to: ${!var_name}"
      return 0
    fi
  done

  # If no path is found
  echo "ERROR: Path '$value' could not be resolved under any base directory."
  return 1
}

enforce_abs_path COLLAPSED_LFC
enforce_abs_path DR_PATH

args=(
--collapsed_lfc "$COLLAPSED_LFC"
--drc_file "$DR_PATH"
--build_dir "$BUILD_DIR"
--sig_cols "$SIG_COLS"
--univariate_biomarker "$UNIVARIATE_BIOMARKER"
--multivariate_biomarker "$MULTIVARIATE_BIOMARKER"
--biomarker_file "$BIOMARKER_FILE"
--dr_column "$DR_COLUMN"
--collapsed_l2fc_column "$COLLAPSED_L2FC_COLUMN"
--lfc_biomarker "$LFC_BIOMARKER"
--auc_biomarker "$AUC_BIOMARKER"
--out_path "$BUILD_DIR/biomarker"
)

echo Rscript biomarker.R "${args[@]}"
Rscript biomarker.R "${args[@]}"
