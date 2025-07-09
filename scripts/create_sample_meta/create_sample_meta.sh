#!/bin/bash

echo Starting create_sample_meta...

export API_KEY=$(cat /local/jenkins/.clue_api_key)
export API_URL="https://api.clue.io/api/"

if [ -z "$BUILD_DIR" ]
then
	echo BUILD_DIR not specified
    exit -1
fi

echo Build dir is: $BUILD_DIR

echo $R_LIBS

args=(
--screen "$SCREEN"
--build_dir "$BUILD_DIR"
--api_key "$API_KEY"
--pert_plates "$PERT_PLATES"
--build_name "$BUILD_NAME"
)

echo python3 create_sample_meta/create_sample_meta.py "${args[@]}"

python3 create_sample_meta/create_sample_meta.py "${args[@]}"
