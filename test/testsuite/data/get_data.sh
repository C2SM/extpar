#!/bin/bash

# this scripts downloads test data for the Extpar testsuite

# echo on
set -x

# go to data directory in case script is invoked from top-level directory
test -d src
if [ $? -eq 0 ] ; then
  cd data
fi

# mch
test -d mch || exit 1
cd mch/c7_globe
wget --quiet 'ftp://iacftp.ethz.ch/pub_read/stelliom/external_parameter_mch_c7_PR417.nc'
cd -

test -d mch || exit 1
cd mch/c1_aster
wget --quiet 'ftp://iacftp.ethz.ch/pub_read/stelliom/external_parameter_mch_c1_PR417.nc'
cd -

# clm
test -d clm || exit 1
cd clm/12km_globe
wget --quiet 'ftp://iacftp.ethz.ch/pub_read/stelliom/external_parameter_12km_globe_PR417.nc'
cd -

cd clm/ecoclimap_sg
wget --quiet 'ftp://iacftp.ethz.ch/pub_read/stelliom/icon_grid_bolivia.nc'
wget --quiet 'ftp://iacftp.ethz.ch/pub_read/stelliom/external_parameter_icon_eco_PR417.nc'
cd -

# dwd
test -d dwd || exit 1

cd dwd/icon_d2
wget --quiet 'ftp://iacftp.ethz.ch/pub_read/stelliom/icon_grid_DOM01.nc'
wget --quiet 'ftp://iacftp.ethz.ch/pub_read/stelliom/external_parameter_icon_d2_PR417.nc'
cd -

cd dwd/icon_ecci
wget --quiet 'ftp://iacftp.ethz.ch/pub_read/stelliom/icon_grid_bolivia.nc'
wget --quiet 'ftp://iacftp.ethz.ch/pub_read/stelliom/clim_t2m_icon_ecci.nc'
wget --quiet 'ftp://iacftp.ethz.ch/pub_read/stelliom/clim_tsea_icon_ecci.nc'
wget --quiet 'ftp://iacftp.ethz.ch/pub_read/stelliom/external_parameter_icon_ecci_PR417.nc'
cd -

# mpim
test -d mpim || exit 1
cd mpim/icon_r2b4
wget --quiet 'http://icon-downloads.mpimet.mpg.de/grids/public/mpim/0013/icon_grid_0013_R02B04_G.nc'
wget --quiet 'ftp://iacftp.ethz.ch/pub_read/stelliom/clim_t2m_icon_r2b4.nc'
wget --quiet 'ftp://iacftp.ethz.ch/pub_read/stelliom/clim_tsea_icon_r2b4.nc'
wget --quiet 'ftp://iacftp.ethz.ch/pub_read/stelliom/external_parameter_icon_mpim_PR417.nc'
cd -

# ecmwf
test -d ecmwf || exit 1
cd ecmwf/corine_icon
wget --quiet 'ftp://iacftp.ethz.ch/pub_read/stelliom/corine/icon_grid_0099_R19B10.nc'
wget --quiet 'ftp://iacftp.ethz.ch/pub_read/stelliom/corine/external_parameter_icon_corine_PR417.nc'
wget --quiet 'ftp://iacftp.ethz.ch/pub_read/stelliom/corine/clim_tsea_corine.nc'
wget --quiet 'ftp://iacftp.ethz.ch/pub_read/stelliom/corine/clim_t2m_corine.nc'
cd -
# done
