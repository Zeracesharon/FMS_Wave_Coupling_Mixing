#!/bin/csh

set version=REPRO
set platform="ncrc6.intel23"
set rootdir = `dirname $0`
source ${rootdir}/envs.${platform}

foreach exec ( ww3_grid ww3_ounf ww3_multi ww3_prnc )

    mkdir -p build/$platform/$exec/$version

    (cd build/$platform/$exec/$version; rm -f path_names; \
    ../../../../src/mkmf/bin/list_paths -l ./ ../../../../src/WW3/model/$exec)
    (cd build/${platform}/$exec/$version; \
    ../../../../src/mkmf/bin/mkmf -t ../../../../src/mkmf/templates/ncrc-intel-WW3.mk -p $exec path_names)
    (cd build/${platform}/$exec/$version; make NETCDF=3 REPRO=1 $exec)

end
