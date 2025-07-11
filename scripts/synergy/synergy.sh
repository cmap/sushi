
echo Starting synergy module...

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
enforce_abs_path LFC

args=(
--normalized_counts "$NORMALIZED_COUNTS"
--l2fc "$LFC"
--cell_line_cols "$CELL_LINE_COLS"
--ctrl_cols "$CONTROL_COLS"
--sig_cols "$SIG_COLS"
--control_type "$CTL_TYPES"
--count_col_name "$COUNT_COL_NAME" 
--count_threshold "$COUNT_THRESHOLD"
--l2fc_col "$L2FC_COLUMN"
--viab_cap "$VIABILITY_CAP"
--out "$BUILD_DIR"
)

echo Rscript synergy.R "${args[@]}"
Rscript synergy.R "${args[@]}"