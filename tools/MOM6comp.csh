#!/bin/csh
#

set version=REPRO
set platform="ncrc6.intel23"
set rootdir = `dirname $0`
source ${rootdir}/envs.${platform}

mkdir -p build/${platform}/wave_ice_ocean/$version

(cd build/${platform}/wave_ice_ocean/$version/; rm -f path_names; \
../../../../src/mkmf/bin/list_paths -l ./ \
                                       ../../../../src/Icepack/columnphysics \
                                       ../../../../src/icebergs/src \
                                       ../../../../src/MOM6/config_src/memory/dynamic_symmetric \
				       ../../../../src/MOM6/config_src/drivers/FMS_cap \
				       ../../../../src/MOM6/config_src/external \
				       ../../../../src/MOM6/config_src/infra/FMS2 \
                                       ../../../../src/MOM6/src/{*,*/*}/ \
				       ../../../../src/{FMS/coupler,FMS/include} \
				       ../../../../src/coupler/{shared,full} \
				       ../../../../src/SIS2/config_src/dynamic_symmetric \
                                       ../../../../src/SIS2/src \
				       ../../../../src/{atmos_null,land_null,ice_param,WW3/model/CPL}/ \
)
(cd build/${platform}/wave_ice_ocean/$version/; \
../../../../src/mkmf/bin/mkmf -t ../../../../src/mkmf/templates/ncrc-intel.mk -o "-I../../fms/${version} -I../../ww3lib/${version}" -p MOM6 -l "-L../../fms/${version} -lfms -L../../ww3lib/${version} -lww3" -c '-Duse_AM3_physics -D_USE_LEGACY_LAND_ -DUSE_FMS2_IO' path_names )
(cd build/${platform}/wave_ice_ocean/$version/; make NETCDF=3 $version=1 MOM6 -j)
