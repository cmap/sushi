#!/usr/bin/env bash
set -e

echo "PWD=$(pwd)"
echo "Contents of repo root:"
ls -1

# Build, using the repo root (.) as the context,
# and pointing Docker at the Dockerfile under ./docker/
podman build \
  --pull \
  --no-cache \
  -f docker/Dockerfile \
  -t sushi-podman:latest \
  .
