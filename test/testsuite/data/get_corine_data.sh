#!/bin/bash

# this scripts downloads test data for the Extpar testsuite

# echo on
set -x

# go to data directory in case script is invoked from top-level directory
test -d src
if [ $? -eq 0 ] ; then
  cd data
fi

# ecmwf
test -d ecmwf || exit 1
cd ecmwf/corine_icon
wget --quiet 'ftp://iacftp.ethz.ch/pub_read/juckerj/corine/icon_grid_0099_R19B10.nc'
wget --quiet 'ftp://iacftp.ethz.ch/pub_read/juckerj/corine/external_parameter_icon_corine.nc'
wget --quiet 'ftp://iacftp.ethz.ch/pub_read/juckerj/corine/ei_an1986-2015_0099_R19B10_BUFFER.nc'
wget --quiet 'ftp://iacftp.ethz.ch/pub_read/juckerj/corine/ei_2t_an1986-2015_0099_R19B10_BUFFER.nc'

cd -

