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
cd mch/test_1
wget --quiet 'ftp://iacftp.ethz.ch/pub_read/silvertk/external_parameter_mch_cosmo7.nc'
cd -

# dwd
test -d dwd || exit 1
cd dwd/test_1
wget --quiet 'ftp://iacftp.ethz.ch/pub_read/silvertk/external_parameter_dwd.nc'
cd -
# done
