echo Starting growth correction module ....

if [ -z "$BUILD_DIR" ]
then
	echo BUILD_DIR not specified
    exit -1
fi
echo Build dir is: $BUILD_DIR

# Function to help build paths
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
enforce_abs_path GROWTH_ANNOTATIONS

args=(
--l2fc "$LFC"
--growth_annotations "$GROWTH_ANNOTATIONS"
--sig_cols "$SIG_COLS"
)

echo Rscript growth_correction/growth_correction.R "${args[@]}"
Rscript growth_correction/growth_correction.R "${args[@]}"