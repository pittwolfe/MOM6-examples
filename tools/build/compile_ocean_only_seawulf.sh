#!/bin/bash

mkdir -p build/gnu/ocean_only/repro/
(cd build/gnu/ocean_only/repro/; 
rm -f path_names; \
../../../../src/mkmf/bin/list_paths -l ./ ../../../../src/MOM6/{config_src/infra/FMS1,config_src/memory/dynamic_symmetric,config_src/drivers/solo_driver,config_src/external,src/{*,*/*}}/ ; \
../../../../src/mkmf/bin/mkmf -t ../../../../macOS-gnu10-openmpi.mk -o '-I../../shared/repro' -p MOM6 -l '-L../../shared/repro -lfms' path_names)
(cd build/intel/ocean_only/repro/; make NETCDF=4 REPRO=1 MOM6 -j)
