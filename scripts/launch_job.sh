#!/bin/bash -e

# Function to get parameter from JSON or fallback to Jenkins parameter
get_param() {
  local param_name="$1"
  local json_value=$(jq -r ".$param_name // empty" "$CONFIG_FILE")

  if [ -z "$json_value" ]; then
    echo "${!param_name}"
  else
    echo "$json_value"
  fi
}

# Check for the presence of config file
CONFIG_FILE="$BUILD_DIR/config.json"

if [ -f "$CONFIG_FILE" ]; then
  echo "Configuration file found. Loading parameters from $CONFIG_FILE."
else
  echo "Configuration file not found. Using Jenkins parameters."
fi

# List of parameters
PARAMS=(
  SEQ_TYPE API_URL BUILD_DIR INDEX_1 INDEX_2 BARCODE_SUFFIX REVERSE_INDEX2
  SAMPLE_META CONTROL_BARCODE_META CTL_TYPES ID_COLS SIG_COLS
  RUN_NORM CONTROL_COLS COUNT_THRESHOLD COUNT_COL_NAME BUILD_NAME
  CONVERT_SUSHI PULL_POOL_ID RUN_EPS_QC PSEUDOCOUNT REMOVE_DATA DAYS
  SEQUENCING_INDEX_COLS RAW_COUNTS CELL_SET_META CELL_LINE_META FILTERED_COUNTS
  LFC COUNTS ANNOTATED_COUNTS COLLAPSED_VALUES NORMALIZED_COUNTS ASSAY_POOL_META
)

# Load parameters
for param in "${PARAMS[@]}"; do
  declare "$param=$(get_param "$param")"
done

# API_KEY is an exception and loaded directly from file (for now)
API_KEY=$(cat /local/jenkins/.clue_api_key)

# Verify that the Git plugin has checked out the correct branch
echo "Current Git branch:"
git -C "$WORKSPACE" rev-parse --abbrev-ref HEAD

# Check if the image exists locally
if /usr/bin/podman images | grep -qE '^localhost/sushi-podman\s+latest\s'; then
  echo "The image 'localhost/sushi-podman:latest' exists locally."
else
  echo "Error: The image 'localhost/sushi-podman:latest' does not exist locally."
  exit 1
fi

# Get the script name from the argument
SCRIPT_NAME="$1"
if [ -z "$SCRIPT_NAME" ]; then
  echo "Error: No script name provided."
  exit 1
fi

# Debugging: Verify the script exists in the Jenkins workspace
echo "Verifying the script location in the Jenkins workspace:"
ls -l "$WORKSPACE/scripts/$SCRIPT_NAME"

# Ensure the script has execute permissions
chmod +x "$WORKSPACE/scripts/$SCRIPT_NAME"

# Run the Podman command using the image tag
echo "Running in container:"
/usr/bin/podman run --rm --user root \
  --entrypoint /bin/bash \
  -e SEQ_TYPE="$SEQ_TYPE" \
  -e API_URL="$API_URL" \
  -e API_KEY="$API_KEY" \
  -e BUILD_DIR="$BUILD_DIR" \
  -e INDEX_1="$INDEX_1" \
  -e INDEX_2="$INDEX_2" \
  -e BARCODE_SUFFIX="$BARCODE_SUFFIX" \
  -e REVERSE_INDEX2="$REVERSE_INDEX2" \
  -e SAMPLE_META="$SAMPLE_META" \
  -e CONTROL_BARCODE_META="$CONTROL_BARCODE_META" \
  -e CTL_TYPES="$CTL_TYPES" \
  -e ID_COLS="$ID_COLS" \
  -e SIG_COLS="$SIG_COLS" \
  -e RUN_NORM="$RUN_NORM" \
  -e CONTROL_COLS="$CONTROL_COLS" \
  -e COUNT_THRESHOLD="$COUNT_THRESHOLD" \
  -e COUNT_COL_NAME="$COUNT_COL_NAME" \
  -e BUILD_NAME="$BUILD_NAME" \
  -e CONVERT_SUSHI="$CONVERT_SUSHI" \
  -e PULL_POOL_ID="$PULL_POOL_ID" \
  -e RUN_EPS_QC="$RUN_EPS_QC" \
  -e PSEUDOCOUNT="$PSEUDOCOUNT" \
  -e REMOVE_DATA="$REMOVE_DATA" \
  -e DAYS="$DAYS" \
  -e SEQUENCING_INDEX_COLS="$SEQUENCING_INDEX_COLS" \
  -e RAW_COUNTS="$RAW_COUNTS" \
  -e CELL_SET_META="$CELL_SET_META" \
  -e CELL_LINE_META="$CELL_LINE_META" \
  -e FILTERED_COUNTS="$FILTERED_COUNTS" \
  -e LFC="$LFC" \
  -e COUNTS="$COUNTS" \
  -e ANNOTATED_COUNTS="$ANNOTATED_COUNTS" \
  -e COLLAPSED_VALUES="$COLLAPSED_VALUES" \
  -e NORMALIZED_COUNTS="$NORMALIZED_COUNTS" \
  -e ASSAY_POOL_META="$ASSAY_POOL_META" \
  -v "$WORKSPACE:/workspace" \
  -v /local/jenkins/.clue_api_key:/local/jenkins/.clue_api_key \
  -v /cmap/data/vdb/prismSeq:/data \
  -v "$BUILD_DIR:$BUILD_DIR" \
  -w /workspace/scripts \
  localhost/sushi-podman:latest \
  ./"$SCRIPT_NAME"
