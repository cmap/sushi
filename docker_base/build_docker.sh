#!/usr/bin/env bash

#change the version number for each new build
podman build --platform linux/amd64 -t prismcmap/base-clue-pseq:develop --rm=true .
