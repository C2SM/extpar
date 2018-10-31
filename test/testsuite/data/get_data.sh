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
wget --quiet 'ftp://iacftp.ethz.ch/pub_read/silvertk/external_parameter_mch_c7.globe.nc'
cd -

test -d mch || exit 1
cd mch/c7_aster
wget --quiet 'ftp://iacftp.ethz.ch/pub_read/silvertk/external_parameter_mch_c7.aster.nc'
cd -

test -d mch || exit 1
cd mch/c1_aster
wget --quiet 'ftp://iacftp.ethz.ch/pub_read/silvertk/external_parameter_mch_cosmo1.nc'
cd -

# clm
test -d clm || exit 1
cd clm/12km_globe
wget --quiet 'ftp://iacftp.ethz.ch/pub_read/silvertk/external_parameter_12km.globe.nc'
cd -

# dwd
test -d dwd || exit 1
cd dwd/cde2_globe
wget --quiet 'ftp://iacftp.ethz.ch/pub_read/silvertk/external_parameter_dwd.nc'
cd -
# done
