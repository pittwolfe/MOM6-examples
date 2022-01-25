#!/bin/bash

mkdir -p build/gnu/shared/repro/
(cd build/gnu/shared/repro/; rm -f path_names; \
../../../../src/mkmf/bin/list_paths ../../../../src/FMS2; \
../../../../src/mkmf/bin/mkmf -t ../../../../macOS-gnu10-openmpi.mk -p libfms.a -c "-Duse_libMPI -Duse_netCDF -DSPMD" path_names)
(cd build/gnu/shared/repro/; make NETCDF=4 REPRO=1 libfms.a -j)
