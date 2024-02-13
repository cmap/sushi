#!/usr/bin/env bash

# read in flagged arguments
while getopts ":d:c:p:o:" arg; do
  case $arg in
    d) # specify input folder
      data_file=${OPTARG};;
    c) # specifcy output folder
      criteria_file=${OPTARG};;
    p) # specifcy output folder
      project_code=${OPTARG};;
    o) # specify build nam e
      output_file=${OPTARG};;
  esac
done

chmod +x /main.py
chmod +x /functions/remove_functions.py
export HDF5_USE_FILE_LOCKING=FALSE
echo "${data_file}" "${criteria_file}" "${project_code}" "${output_file}"

args=(
  -d "${data_file}"
  -c "${criteria_file}"
  -p "${project_code}"
  -o "${output_file}"
)

if [[ ! -z $exclude ]]
then
  args+=(-x "${exclude}")
fi

python /main.py "${args[@]}"

exit_code=$?

echo "$exit_code"
exit $exit_code
