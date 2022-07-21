#!/usr/bin/env python3 
import argparse
import numpy as np

# extpar modules from lib
import grid_def


parser = argparse.ArgumentParser(description='Process some integers.')

parser.add_argument(
    '--input_cosmo_grid',
    type=str,
    required=True,
    help=
    'Fortran Namelist "INPUT_COSMO_GRID"')
args = parser.parse_args()

tg = grid_def.CosmoGrid(args.input_cosmo_grid)

zlonmax = np.amax(tg.lons)
zlonmin = np.amin(tg.lons)
zlatmin = np.amin(tg.lats)
zlatmax = np.amax(tg.lats)

# safety check
if zlatmax > 60.0 or zlatmax < -60.0:
    raise ValueError('Domains using Aster cannot exceed 60 N or 60 S')

aster_tiles_lon = np.empty([12,20])
aster_tiles_lat = np.empty([12,20])

aster_lon=-180.0
aster_lat=60.0
for j in range(0,20):
    for i in range(0,12):
        aster_tiles_lon[i,j] = aster_lon + float(i * 30)
        aster_tiles_lat[i,j] = aster_lat - float(j * 6)

ilon_min = 0
ilon_max = 0
ilat_min = 0
ilat_max = 0

for j in range(0,20):
    for i in range(0,12):
        if aster_tiles_lon[i,j] < zlonmin:
            ilon_min = i
        if aster_tiles_lon[i,j] < zlonmax: 
            ilon_max = i
        if aster_tiles_lat[i,j] > zlatmin: 
            ilat_max = j
        if aster_tiles_lat[i,j] > zlatmax:
            ilat_min = j

ntiles_column = ilon_max - ilon_min + 1
ntiles_row = ilat_max - ilat_min + 1

aster_files = np.empty(240,int)

icount = 0
for j in range(ilat_min,ilat_max+1):
    for i in range(ilon_min,ilon_max+1):
        aster_files[icount] = int(1 + i + 12*j)
        icount += 1

print(f'NTILES_COLUMN={ntiles_column}')
print(f'NTILES_ROW={ntiles_row}')
for idx in range(1,icount + 1):
    print(f"raw_data_aster_T{idx:03}='ASTER_orig_T{aster_files[idx-1]:03}.nc'")
