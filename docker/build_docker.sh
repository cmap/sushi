#!/usr/bin/env bash

cd ../

docker build -t prismcmap/sushi:develop --rm=true -f docker/Dockerfile . --no-cache
