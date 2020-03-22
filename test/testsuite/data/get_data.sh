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
wget --quiet 'ftp://iacftp.ethz.ch/pub_read/juckerj/external_parameter_mch_c7_PR159.nc'
cd -

test -d mch || exit 1
cd mch/c1_aster
wget --quiet 'ftp://iacftp.ethz.ch/pub_read/juckerj/external_parameter_mch_c1_PR159.nc'
cd -

# clm
test -d clm || exit 1
cd clm/12km_globe
wget --quiet 'ftp://iacftp.ethz.ch/pub_read/silvertk/external_parameter_12km_globe_PR143.nc'
cd -

# dwd
test -d dwd || exit 1
cd dwd/cde2_globe
wget --quiet 'ftp://iacftp.ethz.ch/pub_read/juckerj/external_parameter_cde2_dwd_PR159.nc'
cd -

cd dwd/icon_d2
wget --quiet 'ftp://iacftp.ethz.ch/pub_read/juckerj/icon_grid_DOM01.nc'
wget --quiet 'ftp://iacftp.ethz.ch/pub_read/juckerj/ei_2t_an1986-2015_domain2_DOM01_BUFFER.nc'
wget --quiet 'ftp://iacftp.ethz.ch/pub_read/juckerj/ei_an1986-2015_domain2_DOM01_BUFFER.nc'
wget --quiet 'ftp://iacftp.ethz.ch/pub_read/juckerj/external_parameter_icon_domain2_DOM01_tiles.nc'
cd -

# mpim
test -d mpim || exit 1
cd mpim/icon_r2b4
wget --quiet 'http://icon-downloads.mpimet.mpg.de/grids/public/mpim/0013/icon_grid_0013_R02B04_G.nc'
wget --quiet 'ftp://iacftp.ethz.ch/pub_read/silvertk/ei_t2m_an1986-2015_0013_R02B04_G_BUFFER.nc'
wget --quiet 'ftp://iacftp.ethz.ch/pub_read/silvertk/ei_sst_an1986-2015_0013_R02B04_G_BUFFER.nc'
wget --quiet 'ftp://iacftp.ethz.ch/pub_read/silvertk/external_parameter_icon_mpim_PR125.nc'
cd -

# intel references
cd intel

# mch
test -d mch || exit 1
cd mch/c7_globe
wget --quiet 'ftp://iacftp.ethz.ch/pub_read/juckerj/external_parameter_intel_mch_c7_PR159.nc'
cd -

test -d mch || exit 1
cd mch/c1_aster
wget --quiet 'ftp://iacftp.ethz.ch/pub_read/juckerj/external_parameter_intel_mch_c1_PR159.nc'
cd -

# clm
test -d clm || exit 1
cd clm/12km_globe
wget --quiet 'ftp://iacftp.ethz.ch/pub_read/juckerj/external_parameter_intel_12km_globe_PR143.nc'
cd -

# dwd
test -d dwd || exit 1
cd dwd/cde2_globe
wget --quiet 'ftp://iacftp.ethz.ch/pub_read/silvertk/external_parameter_intel_cde2_globe_PR159.nc'
cd -

# mpim
test -d mpim || exit 1
cd mpim/icon_r2b4
wget --quiet 'ftp://iacftp.ethz.ch/pub_read/silvertk/external_parameter_intel_icon_mpim_PR138.nc'
cd -
cd ..

# done
