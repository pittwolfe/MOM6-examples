#!/bin/bash

# uncomment if we want reproduction
#REPRO="REPRO=1"

if [ ! -z "$REPRO" ]
then
  DIR="repro"
else
  DIR="opt"
fi

mkdir -p ocean_only/$DIR
cd ocean_only/$DIR

if [[ -e Makefile ]] ; then
  make clean
  rm Makefile
fi

rm -f path_names
../../../../src/mkmf/bin/list_paths -l ./ ../../../../src/MOM6/{config_src/infra/FMS2,config_src/memory/dynamic_symmetric,config_src/drivers/solo_driver,config_src/external,src/{*,*/*}}/
../../../../src/mkmf/bin/mkmf -t ../../seawulf-intel.mk -o "-I../../shared/$DIR" -p MOM6 -l "-L../../shared/$DIR -lfms" path_names
make NETCDF=4 $REPRO MOM6 -j
