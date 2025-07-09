#!/bin/bash

echo Starting load_bq...

if [ -z "$BUILD_DIR" ]; then
    echo "BUILD_DIR not specified"
    exit -1
fi

echo "Build dir is: $BUILD_DIR"

args=(
--build_path "$BUILD_DIR"
)

export GOOGLE_APPLICATION_CREDENTIALS=/local/jenkins/.gcp_credentials.json

echo python3 load_bq/load_bq.py "${args[@]}"

python3 load_bq/load_bq.py "${args[@]}"
