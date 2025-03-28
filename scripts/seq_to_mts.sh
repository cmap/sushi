#!/bin/bash

echo Starting seq_to_mts...

if [ -z "$BUILD_DIR" ]; then
    echo "BUILD_DIR not specified"
    exit -1
fi

echo "Build dir is: $BUILD_DIR"
echo "Output dir is: $BUILD_DIR/sync_to_s3"

# Create the output directory if it does not exist
if [ ! -d "$BUILD_DIR/sync_to_s3" ]; then
    echo "Creating output directory: $BUILD_DIR/sync_to_s3"
    mkdir -p "$BUILD_DIR/sync_to_s3"
else
    echo "Output directory already exists: $BUILD_DIR/sync_to_s3"
fi

args=(
--build_path "$BUILD_DIR"
--out "$BUILD_DIR/sync_to_s3"
--build_name "$BUILD_NAME"
--days "$DAYS"
--config "config.json"
--api_key "$API_KEY"
)

echo python3 seq_to_mts.py "${args[@]}"

python3 seq_to_mts.py "${args[@]}"
