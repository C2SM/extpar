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
wget --quiet 'ftp://iacftp.ethz.ch/pub_read/silvertk/external_parameter_mch_c7_globe_5.2.nc'
cd -

test -d mch || exit 1
cd mch/c7_aster
wget --quiet 'ftp://iacftp.ethz.ch/pub_read/silvertk/external_parameter_mch_c7_aster_5.1.nc'
cd -

test -d mch || exit 1
cd mch/c1_aster
wget --quiet 'ftp://iacftp.ethz.ch/pub_read/silvertk/external_parameter_mch_c1_5.1.nc'
cd -

# clm
test -d clm || exit 1
cd clm/12km_globe
wget --quiet 'ftp://iacftp.ethz.ch/pub_read/silvertk/external_parameter_12km_globe_5.1b.nc'
cd -

# dwd
test -d dwd || exit 1
cd dwd/cde2_globe
wget --quiet 'ftp://iacftp.ethz.ch/pub_read/silvertk/external_parameter_dwd_5.1.nc'
cd -

# mpim
test -d mpim || exit 1
cd mpim/icon_r2b4
wget --quiet 'http://icon-downloads.mpimet.mpg.de/grids/public/mpim/0013/icon_grid_0013_R02B04_G.nc'
wget --quiet 'ftp://iacftp.ethz.ch/pub_read/silvertk/ei_t2m_an1986-2015_0013_R02B04_G_BUFFER.nc'
wget --quiet 'ftp://iacftp.ethz.ch/pub_read/silvertk/ei_sst_an1986-2015_0013_R02B04_G_BUFFER.nc'
wget --quiet 'ftp://iacftp.ethz.ch/pub_read/silvertk/external_parameter_icon_mpim_5.2b.nc'
cd -

# done
