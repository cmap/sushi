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
  SEQ_TYPE BUILD_DIR INDEX_1 INDEX_2 BARCODE_SUFFIX REVERSE_INDEX2 BUILD_NAME
  CONVERT_SUSHI PULL_POOL_ID RUN_EPS_QC REMOVE_DATA DAYS COUNTS API_URL
  # sushi input files
  RAW_COUNTS_UNCOLLAPSED SAMPLE_META CELL_SET_AND_POOL_META CELL_LINE_META CONTROL_BARCODE_META
  # susi files
  PRISM_BARCODE_COUNTS UNKNOWN_BARCODE_COUNTS ANNOTATED_COUNTS FILTERED_COUNTS NORMALIZED_COUNTS 
  LFC COLLAPSED_LFC
  # collate_fastq_reads parameters
  SEQUENCING_INDEX_COLS ID_COLS LOW_ABUNDANCE_THRESHOLD CHUNK_SIZE BARCODE_COL
  # normalize parameters
  PSEUDOCOUNT
  # compute_l2fc parameters
  CELL_LINE_COLS SIG_COLS CONTROL_COLS COUNT_COL_NAME CTL_TYPES COUNT_THRESHOLD
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
  -e CELL_SET_AND_POOL_META="$CELL_SET_AND_POOL_META" \
  -e CELL_LINE_META="$CELL_LINE_META" \
  -e CONTROL_BARCODE_META="$CONTROL_BARCODE_META" \
  -e ASSAY_POOL_META="$ASSAY_POOL_META" \
  -e RAW_COUNTS_UNCOLLAPSED="$RAW_COUNTS_UNCOLLAPSED"\
  -e PRISM_BARCODE_COUNTS="$PRISM_BARCODE_COUNTS"\
  -e UNKNOWN_BARCODE_COUNTS="$UNKNOWN_BARCODE_COUNTS"\
  -e ANNOTATED_COUNTS="$ANNOTATED_COUNTS" \
  -e FILTERED_COUNTS="$FILTERED_COUNTS" \
  -e NORMALIZED_COUNTS="$NORMALIZED_COUNTS" \
  -e LFC="$LFC" \
  -e COLLAPSED_LFC="$COLLAPSED_LFC" \
  -e SEQUENCING_INDEX_COLS="$SEQUENCING_INDEX_COLS" \
  -e ID_COLS="$ID_COLS" \
  -e CELL_LINE_COLS="$CELL_LINE_COLS" \
  -e SIG_COLS="$SIG_COLS" \
  -e CONTROL_COLS="$CONTROL_COLS" \
  -e BARCODE_COL="$BARCODE_COL" \
  -e PSEUDOCOUNT="$PSEUDOCOUNT" \
  -e COUNT_COL_NAME="$COUNT_COL_NAME" \
  -e CTL_TYPES="$CTL_TYPES" \
  -e CHUNK_SIZE="$CHUNK_SIZE" \
  -e COUNT_THRESHOLD="$COUNT_THRESHOLD" \
  -e LOW_ABUNDANCE_THRESHOLD="$LOW_ABUNDANCE_THRESHOLD" \
  -e BUILD_NAME="$BUILD_NAME" \
  -e CONVERT_SUSHI="$CONVERT_SUSHI" \
  -e PULL_POOL_ID="$PULL_POOL_ID" \
  -e RUN_EPS_QC="$RUN_EPS_QC" \
  -e REMOVE_DATA="$REMOVE_DATA" \
  -e DAYS="$DAYS" \
  -e COUNTS="$COUNTS" \
  -v "$WORKSPACE:/workspace" \
  -v /cmap/tools/analysis2clue/credentials:/root/.aws/credentials:ro \
  -v /local/jenkins/.clue_api_key:/local/jenkins/.clue_api_key:ro \
  -v /cmap/data/vdb/prismSeq:/data \
  -v "$BUILD_DIR:$BUILD_DIR" \
  -w /workspace/scripts \
  localhost/sushi-podman:latest \
  ./"$SCRIPT_NAME"
