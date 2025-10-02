#!/bin/bash

echo Starting generate_qc_tables_2...

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

enforce_abs_path LFC

args=(
--out "$BUILD_DIR"
--l2fc "$LFC"
)

echo Rscript qc_tables_2/qc_tables_2.R "${args[@]}"
Rscript qc_tables_2/qc_tables_2.R "${args[@]}"
