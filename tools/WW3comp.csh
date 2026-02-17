#!/bin/csh
#

set version=REPRO
set platform="ncrc6.intel23"
set rootdir = `dirname $0`
source ${rootdir}/envs.${platform}

mkdir -p build/${platform}/ww3lib/$version

(cd build/${platform}/ww3lib/$version; rm -f path_names; \
../../../../src/mkmf/bin/list_paths -l ./ ../../../../src/WW3/model/ww3_multi)
(cd build/${platform}/ww3lib/$version; \
../../../../src/mkmf/bin/mkmf -t ../../../../src/mkmf/templates/ncrc-intel-WW3.mk -p libww3.a -c "-Duse_libMPI -Duse_netCDF -DSPMD " path_names)
(cd build/${platform}/ww3lib/$version; make NETCDF=4 REPRO=1 libww3.a -j)
