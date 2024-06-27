#!/usr/bin/env bash

# Ensure script stops on first error
set -e

# Specify the directory containing the Dockerfile
BUILD_CONTEXT_DIR="./docker_base"

# Run the podman build command with the specified context directory
podman build -t sushi-podman -f ${BUILD_CONTEXT_DIR}/Dockerfile ${BUILD_CONTEXT_DIR}
