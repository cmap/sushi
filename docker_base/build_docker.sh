#!/usr/bin/env bash

# Ensure script stops on first error
set -e

# Change the version number for each new build
# Specify the directory containing the Dockerfile
BUILD_CONTEXT_DIR="./docker_base"

# Run the podman build command with the specified context directory
podman build --platform linux/amd64 -t podman_base_image:develop --rm=true -f ${BUILD_CONTEXT_DIR}/Dockerfile ${BUILD_CONTEXT_DIR}
