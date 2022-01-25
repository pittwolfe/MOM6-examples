#!/bin/bash

# uncomment if we want reproduction
#REPRO="REPRO=1"

if [ ! -z "$REPRO" ]
then
  DIR="repro"
else
  DIR="opt"
fi

mkdir -p shared/$DIR
cd shared/$DIR
rm -f path_names
../../../../src/mkmf/bin/list_paths ../../../../src/FMS2
../../../../src/mkmf/bin/mkmf -t ../../seawulf-intel.mk -p libfms.a -c "-Duse_libMPI -Duse_netCDF -DSPMD" path_names
make NETCDF=4 $REPRO libfms.a -j
cd ../..

