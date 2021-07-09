#!/bin/bash

mkdir -p build/intel/shared/repro/
(cd build/intel/shared/repro/; rm -f path_names; \
../../../../src/mkmf/bin/list_paths ../../../../src/FMS; \
../../../../src/mkmf/bin/mkmf -t ../../linux-intel.mk -p libfms.a -c "-Duse_libMPI -Duse_netCDF -DSPMD" path_names)
(cd build/intel/shared/repro/; make NETCDF=4 REPRO=1 libfms.a -j)
