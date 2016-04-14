#!/bin/bash
eval $(docker-machine env default)
set -ex

docker build -t dragonfly/sizespectra .

docker run --rm -i -v $PWD/..:/work -w /work/large-marine-sizespectrum \
  dragonfly/sizespectra ./build.sh
  