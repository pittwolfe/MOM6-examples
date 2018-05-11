#!/bin/bash

mkdir -p build/intel/ocean_only/repro/
(cd build/intel/ocean_only/repro/; rm -f path_names; \
../../../../src/mkmf/bin/list_paths ./ ../../../../src/MOM6/{config_src/dynamic,config_src/solo_driver,src/{*,*/*}}/ ; \
../../../../src/mkmf/bin/mkmf -t ../../linux-intel.mk -o '-I../../shared/repro' -p MOM6 -l '-L../../shared/repro -lfms' -c '-Duse_libMPI -Duse_netCDF -DSPMD' path_names)

(cd build/intel/ocean_only/repro/; make NETCDF=4 REPRO=1 MOM6 -j)
