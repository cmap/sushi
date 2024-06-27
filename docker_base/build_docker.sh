#!/usr/bin/env bash

# Ensure script stops on first error
set -e

# Print the current working directory
echo "Current working directory: $(pwd)"

# List files in the current directory for debugging
ls -la

# Specify the directory containing the Dockerfile
BUILD_CONTEXT_DIR="./docker_base"

# Ensure the build context directory exists
if [ ! -d "$BUILD_CONTEXT_DIR" ]; then
  echo "Error: Build context directory '$BUILD_CONTEXT_DIR' does not exist."
  exit 1
fi

# Change to the build context directory
cd ${BUILD_CONTEXT_DIR}

# Run the podman build command with the specified context directory
podman build -t sushi-podman .
