#!/bin/bash

set -ex

docker build -t dragonfly/sizespectra .

docker run --rm -it -v $PWD/..:/work -w /work/ \
  dragonfly/sizespectra ./build.sh