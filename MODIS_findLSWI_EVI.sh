#!/bin/bash

gdalwarp -of netCDF -s_srs '+proj=sinu +R=6371007.181 +nadgrids=@null +wktext' -r near -t_srs '+proj=longlat +datum=WGS84 +no_defs' HDF4_EOS:EOS_GRID:"MOD09A1_2019-04-23_Down.hdf":MOD_Grid_500m_Surface_Reflectance:sur_refl_b01 band1.nc
gdalwarp -of netCDF -s_srs '+proj=sinu +R=6371007.181 +nadgrids=@null +wktext' -r near -t_srs '+proj=longlat +datum=WGS84 +no_defs' HDF4_EOS:EOS_GRID:"MOD09A1_2019-04-23_Down.hdf":MOD_Grid_500m_Surface_Reflectance:sur_refl_b02 band2.nc
gdalwarp -of netCDF -s_srs '+proj=sinu +R=6371007.181 +nadgrids=@null +wktext' -r near -t_srs '+proj=longlat +datum=WGS84 +no_defs' HDF4_EOS:EOS_GRID:"MOD09A1_2019-04-23_Down.hdf":MOD_Grid_500m_Surface_Reflectance:sur_refl_b03 band3.nc
gdalwarp -of netCDF -s_srs '+proj=sinu +R=6371007.181 +nadgrids=@null +wktext' -r near -t_srs '+proj=longlat +datum=WGS84 +no_defs' HDF4_EOS:EOS_GRID:"MOD09A1_2019-04-23_Down.hdf":MOD_Grid_500m_Surface_Reflectance:sur_refl_b06 band6.nc

ncap2 -s 'Band1=Band1/10000.0' band1.nc band1_v1.nc
ncap2 -s 'Band1=Band1/10000.0' band2.nc band2_v1.nc
ncap2 -s 'Band1=Band1/10000.0' band3.nc band3_v1.nc
ncap2 -s 'Band1=Band1/10000.0' band6.nc band6_v1.nc

rm -rf band1.nc
rm -rf band2.nc
rm -rf band3.nc
rm -rf band6.nc

ncrename -v Band1,band1 band1_v1.nc
ncrename -v Band1,band2 band2_v1.nc
ncrename -v Band1,band3 band3_v1.nc
ncrename -v Band1,band6 band6_v1.nc

ncks -A -v band2 band2_v1.nc band1_v1.nc
ncks -A -v band3 band3_v1.nc band1_v1.nc
ncks -A -v band6 band6_v1.nc band1_v1.nc

rm -rf band2_v1.nc
rm -rf band3_v1.nc
rm -rf band6_v1.nc

mv band1_v1.nc All_bands.nc


ncap2 -s 'CHISL=band2-band6' All_bands.nc All_bands_v2.nc
ncap2 -s 'ZN=band2+band6' All_bands_v2.nc All_bands_v3.nc
ncap2 -s 'LSWI=CHISL/ZN' All_bands_v3.nc All_bands_v4.nc
ncap2 -s 'CHISL=band2-band1' All_bands_v4.nc All_bands_v5.nc
ncap2 -s 'CHISL=CHISL*2.5' All_bands_v5.nc All_bands_v6.nc
ncap2 -s 'ZN=6.0*band1-7.5*band3' All_bands_v6.nc All_bands_v7.nc
ncap2 -s 'ZN=ZN+band2+1' All_bands_v7.nc All_bands_v8.nc
ncap2 -s 'EVI=CHISL/ZN' All_bands_v8.nc All_bands_2019-04-23_Down.nc
rm -rf All_bands_v*
