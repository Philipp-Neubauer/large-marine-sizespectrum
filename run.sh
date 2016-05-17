#!/bin/bash
eval $(docker-machine env default)
set -ex

#docker pull mmmmk/size-spectra-docker

docker run --rm -i -v $PWD/..:/work -w /work/large-marine-sizespectrum \
  mmmmk/size-spectra-docker ./build.sh
  