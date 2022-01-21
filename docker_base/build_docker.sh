#!/usr/bin/env bash

#change the version number for each new build
docker build -t prismcmap/base-clue-pseq:latest -t prismcmap/base-clue-pseq:v0.0.1 --rm=true .
