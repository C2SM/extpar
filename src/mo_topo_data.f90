!   Fortran module with data files for GLOBE and ASTER data
!
!  History:
!  Version                  Date                 Name
!  ------------------------ -------------------- -----------------
!  V2_0                     2012/11/21           Martina Messmer
!   initial release
!   replaces the older module mo_globe_data.f90
!   introduce the ASTER topography as external parameters
!   contains a deallocation subroutine
!  V2_0_3                   2014/09/17           Burkhardt Rockel
!   Added use of directory information to access raw data files
!
!  Code Description:
!  Language: Fortran 90
!=================================================================
MODULE mo_topo_data

  ! Modules used:

  USE mo_kind,               ONLY: wp, i4

  USE mo_grid_structures,    ONLY:  reg_lonlat_grid

  USE mo_utilities_extpar,   ONLY:  abort_extpar

  USE mo_io_utilities,       ONLY:  check_netcdf

  USE mo_io_units,           ONLY:  filename_max

  USE mo_topo_tg_fields, ONLY: fr_land_topo,  &
       &                       hh_topo,       &
       &                       hh_topo_max,   &
       &                       hh_topo_min,   &
       &                       stdh_topo,     &
       &                       theta_topo,    &
       &                       aniso_topo,    &
       &                       slope_topo,    &
       &                       z0_topo,       &
       &                       slope_asp_topo,&
       &                       slope_ang_topo,&
       &                       horizon_topo,  &
       &                       skyview_topo

  USE netcdf,       ONLY :    &
       nf90_open,              &
       nf90_close,             &
       nf90_inquire_dimension, &
       nf90_inquire_variable,  &
       nf90_get_att,           &
       nf90_inquire_dimension, &
       nf90_inq_dimid,         &
       nf90_inq_varid,         &
       nf90_get_var,           &
       NF90_NOWRITE

  IMPLICIT NONE

  PRIVATE


  PUBLIC :: ntiles,                  &   ! number of tiles in GLOBE / ASTER
       max_tiles,               &   ! maximal possible number of tiles that can be read
       nc_tot,                  &   ! total number of columns in GLOBE / ASTER data set
       nr_tot,                  &   ! total number of rows in GLOBE / ASTER data set
       nc_tile ,                &   ! total number of columns in one GLOBE / ASTER tile
       tiles_lon_min,           &   ! starting longitude of every GLOBE / ASTER tile
       tiles_lon_max,           &   ! ending longitude of every GLOBE / ASTER tile
       tiles_lat_min,           &   ! starting latitude of every GLOBE / ASTER tile
       tiles_lat_max,           &   ! ending latitude of every GLOBE / ASTER tile
       tiles_ncolumns,          &   ! number of columns (lonitude increments) in each GLOBE / ASTER tile
       tiles_nrows,             &   ! number of rows (latitude increments) in each GLOBE /ASTER tile
       topo_tiles_grid,         &
       topo_grid,               &
       raw_topo_line,           &
       h_tile_row,              &
       raw_topo_block,          &
       allocate_raw_topo_fields,&
       fill_topo_data,           &  ! subroutine (intent(in) and intent(out))
       num_tiles,                &  ! integer function
       allocate_topo_data,       &  ! subroutine (only intent(in))
       get_fill_value,           &
       get_varname,              &
       undef_topo,               &
       varname,                  &
       itopo_type,               &
       topo_gl,                  &
       topo_aster,               &
       aster_lat_min,            &
       aster_lon_min,            &
       aster_lat_max,            &
       aster_lon_max,            &
       ntiles_row,               &
       ntiles_column,            &
       lradtopo,                 &
       nhori,                    &
       deallocate_topo_fields

  TYPE(reg_lonlat_grid), ALLOCATABLE  :: topo_tiles_grid(:)
  TYPE(reg_lonlat_grid)               :: topo_grid

  INTEGER(i4), ALLOCATABLE :: tiles_ncolumns(:)
  INTEGER(i4), ALLOCATABLE :: tiles_nrows(:)
  INTEGER(i4), ALLOCATABLE :: h_tile_row(:)

  INTEGER(i4) :: ntiles
  INTEGER(i4) :: nc_tot
  INTEGER(i4) :: nr_tot
  INTEGER(i4) :: nc_tile
  INTEGER(i4) :: undef_topo
  INTEGER(i4) :: itopo_type
  INTEGER(i4), PROTECTED :: ntiles_row
  INTEGER(i4), PROTECTED :: ntiles_column
  INTEGER(i4) :: nhori

  INTEGER, PARAMETER:: topo_gl = 1
  INTEGER, PARAMETER:: topo_aster = 2
  INTEGER, PARAMETER:: max_tiles = 1000
  !INTEGER, PARAMETER:: ntiles_topo_column = 4
  !INTEGER, PARAMETER:: ntiles_topo_row    = 4
  !INTEGER, PARAMETER:: ntiles_aster_column = 3      ! needs to be changed if ASTER tiles are added
  !INTEGER, PARAMETER:: ntiles_aster_row    = 10     ! needs to be changed if ASTER tiles are added

  REAL(wp), ALLOCATABLE    :: tiles_lon_min(:)
  REAL(wp), ALLOCATABLE    :: tiles_lon_max(:)
  REAL(wp), ALLOCATABLE    :: tiles_lat_min(:)
  REAL(wp), ALLOCATABLE    :: tiles_lat_max(:)
  REAL(wp), ALLOCATABLE    :: raw_topo_line(:)
  REAL(wp), ALLOCATABLE    :: raw_topo_block(:,:)


  REAL(wp):: aster_lat_min
  REAL(wp):: aster_lat_max
  REAL(wp):: aster_lon_min
  REAL(wp):: aster_lon_max

  LOGICAL      :: lradtopo


  CHARACTER(LEN=80) :: varname


CONTAINS


  SUBROUTINE num_tiles(itopo_type, columns, rows, ntiles) ! it gives the value of the number of tiles depending

    INTEGER, INTENT(in)  :: itopo_type
    INTEGER, INTENT(in)  :: columns
    INTEGER, INTENT(in)  :: rows
    INTEGER, INTENT(out) :: ntiles           ! if the user chooses GLOBE or ASTER

    ntiles_column = columns
    ntiles_row    = rows
    ntiles = ntiles_column * ntiles_row
    PRINT*, 'number of tiles is: ', ntiles

  END SUBROUTINE num_tiles

  SUBROUTINE allocate_topo_data(ntiles)
    ! As it is unknown so far whether GLOBE or ASTER is chosen all parameters must be allocated in a second step.

    IMPLICIT NONE
    INTEGER, INTENT (IN) :: ntiles       ! number of tiles: 36 for ASTER and 16 for GLOBE
    INTEGER :: errorcode

    ALLOCATE (tiles_lon_min(1:ntiles), STAT = errorcode)
    IF (errorcode /= 0) CALL abort_extpar('Cant allocate the vector tiles_lon_min')
    ALLOCATE (tiles_lon_max(1:ntiles), STAT = errorcode)
    IF (errorcode /= 0) CALL abort_extpar('Cant allocate the vector tiles_lon_max')
    ALLOCATE (tiles_lat_min(1:ntiles), STAT = errorcode)
    IF (errorcode /= 0) CALL abort_extpar('Cant allocate the vector tiles_lat_min')
    ALLOCATE (tiles_lat_max(1:ntiles), STAT = errorcode)
    IF (errorcode /= 0) CALL abort_extpar('Cant allocate the vector tiles_lat_max')

    ALLOCATE (tiles_ncolumns(1:ntiles), STAT = errorcode)
    IF (errorcode /= 0) CALL abort_extpar('Cant allocate the vector tiles_ncolumns')
    ALLOCATE (tiles_nrows(1:ntiles), STAT = errorcode)
    IF (errorcode /= 0) CALL abort_extpar('Cant allocate the vector tiles_nrows')

    ALLOCATE (topo_tiles_grid(1:ntiles), STAT = errorcode)
    IF (errorcode /= 0) CALL abort_extpar('Cant allocate the vector topo_tiles_grid')

    tiles_lon_min   = 0.0
    tiles_lon_max   = 0.0
    tiles_lat_min   = 0.0
    tiles_lat_max   = 0.0
    tiles_ncolumns  = 0
    tiles_nrows     = 0
    aster_lat_min   = 0.0
    aster_lat_max   = 0.0
    aster_lon_min   = 0.0
    aster_lon_max   = 0.0

  END SUBROUTINE allocate_topo_data




  SUBROUTINE fill_topo_data(raw_data_orography_path,topo_files,  &
       tiles_lon_min,  &
       tiles_lon_max,  &
                                ! the allocated vectors need to be filled with the respective value.
       tiles_lat_min,  &
       tiles_lat_max,  &
       nc_tot,         &
       nr_tot,         &
       nc_tile)
    IMPLICIT NONE
    SAVE
    CHARACTER (len=filename_max),INTENT(IN) :: raw_data_orography_path
    CHARACTER (len=filename_max),INTENT(IN) :: topo_files(1:ntiles)
    REAL(wp), INTENT(OUT)   :: tiles_lon_min(1:ntiles)
    REAL(wp), INTENT(OUT)   :: tiles_lon_max(1:ntiles)
    REAL(wp), INTENT(OUT)   :: tiles_lat_min(1:ntiles)
    REAL(wp), INTENT(OUT)   :: tiles_lat_max(1:ntiles)
    INTEGER(i4), INTENT(OUT):: nc_tot, nr_tot, nc_tile
    INTEGER(i4)    :: i, errorcode     ! i is a counter, errorcode is used to check if allocation was successful
    INTEGER(i4)    :: ncid
    INTEGER(i4)    :: dimID_lat, dimID_lon, varID_lat, varID_lon
    REAL(wp)       :: half_gridp       ! distance of half a grid point as the grid point is centered on a GLOBE / ASTER pixel


    SELECT CASE (itopo_type)                ! Also topo could additionally be used for SELECT CASE (must first be read in)
    CASE(topo_aster)                                         ! ASTER topography, as it has 36 tiles at the moment.
      PRINT*, 'ASTER is used as topography'
      half_gridp = 1./(3600.*2.)           ! the resolution of the ASTER data is 1./3600. degrees as it is half a grid point
      ! it is additionally divided by 2
    CASE (topo_gl)                                           ! GLOBE topography is composed of 16 tiles
      PRINT*, 'GLOBE is used as topography'
      half_gridp = 1./(120.*2.)                              ! GLOBE resolution is 1./120. degrees (30 arc-seconds)
    END SELECT

    DO i = 1,ntiles

      !_br 17.09.14       CALL check_netcdf(nf90_open(path =TRIM(topo_files(i)), mode = nf90_nowrite, ncid = ncid))
      ! ASTER/GLOBE file is opened (intent(out) is only the ncid)
      CALL check_netcdf(nf90_open(path =TRIM(raw_data_orography_path)//TRIM(topo_files(i)), mode = nf90_nowrite, ncid = ncid), __FILE__, __LINE__)
      CALL check_netcdf(nf90_inq_dimid(ncid,"lon", dimID_lon), __FILE__, __LINE__)
      CALL check_netcdf(nf90_inq_dimid(ncid,"lat", dimID_lat), __FILE__, __LINE__)
      CALL check_netcdf(nf90_inquire_dimension(ncid,dimID_lon, len = tiles_ncolumns(i)), __FILE__, __LINE__)
      CALL check_netcdf(nf90_inquire_dimension(ncid,dimID_lat, len = tiles_nrows(i)), __FILE__, __LINE__)
      CALL check_netcdf(nf90_inq_varid(ncid, "lon", varID_lon), __FILE__, __LINE__)
      CALL check_netcdf(nf90_inq_varid(ncid, "lat", varID_lat), __FILE__, __LINE__)
      CALL check_netcdf(nf90_get_var(ncid, varID_lon, tiles_lon_min(i), start = (/1/)), __FILE__, __LINE__)
      ! reads in the first longitude value of tile i
      CALL check_netcdf(nf90_get_var(ncid, varID_lon, tiles_lon_max(i), start = (/tiles_ncolumns(i)/)), __FILE__, __LINE__)
      ! reads in the last longitude value of tile i
      CALL check_netcdf(nf90_get_var(ncid, varID_lat, tiles_lat_max(i), start = (/1/)), __FILE__, __LINE__)
      ! reads in the first latitude value of tile i
      CALL check_netcdf(nf90_get_var(ncid, varID_lat, tiles_lat_min(i), start = (/tiles_nrows(i)/)), __FILE__, __LINE__)
      ! reads in the last latitude value of tile i
      CALL check_netcdf(nf90_close(ncid), __FILE__, __LINE__)
      ! the netcdf file is closed again
      tiles_lon_min(i) = REAL(NINT(tiles_lon_min(i) - half_gridp)) !< half of a grid point must be
      tiles_lon_max(i) = REAL(NINT(tiles_lon_max(i) + half_gridp)) !< added, as the ASTER/GLOBE data
      tiles_lat_min(i) = REAL(NINT(tiles_lat_min(i) + half_gridp)) !< is located at the pixel center
      tiles_lat_max(i) = REAL(NINT(tiles_lat_max(i) - half_gridp))
    END DO


    SELECT CASE(itopo_type)
    CASE(topo_aster)
      aster_lat_min = MINVAL(tiles_lat_min)
      aster_lat_max = MAXVAL(tiles_lat_max)
      aster_lon_min = MINVAL(tiles_lon_min)
      aster_lon_max = MAXVAL(tiles_lon_max)
    END SELECT

    nc_tot = 0
    nr_tot = 0

    DO i = 1,ntiles
      nc_tot = nc_tot + tiles_ncolumns(i)
      nr_tot = nr_tot + tiles_nrows(i)
    END DO
    nc_tot = nc_tot/ntiles_row
    nr_tot = nr_tot/ntiles_column
    nc_tile = tiles_ncolumns(1)


    SELECT CASE(itopo_type)
    CASE(topo_gl)
      WHERE (tiles_lon_max.GT.(-1*half_gridp).AND.tiles_lon_max.LT.half_gridp) tiles_lon_max = 0.00000
      WHERE (tiles_lat_max.GT.(-1*half_gridp).AND.tiles_lat_max.LT.half_gridp) tiles_lat_max = 0.00000
      ! There are probably some rounding problems, which are removed by this procedure.
    END SELECT


    ALLOCATE (raw_topo_line(1:nc_tot), STAT = errorcode)
    IF (errorcode /= 0) CALL abort_extpar('Cant allocate the vector raw_topo_line')
    ALLOCATE (h_tile_row(1:nc_tile), STAT = errorcode)
    IF (errorcode /= 0) CALL abort_extpar('Cant allocate the vector h_tile_row')

  END SUBROUTINE fill_topo_data




  SUBROUTINE allocate_raw_topo_fields(nrows,ncolumns)
    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: nrows
    INTEGER, INTENT(IN)  :: ncolumns

    INTEGER              :: errorcode

    ALLOCATE (raw_topo_block(1:ncolumns, 1:nrows), STAT = errorcode)
    IF (errorcode /= 0) CALL abort_extpar('Cant allocate the array raw_topo_block')
    raw_topo_block = 0.0

  END SUBROUTINE allocate_raw_topo_fields



  SUBROUTINE get_fill_value(topo_file_1, undef_topo)
    CHARACTER (len=filename_max), INTENT(in) :: topo_file_1
    INTEGER, INTENT(out)           :: undef_topo
    INTEGER(i4)               :: ncid

    SELECT CASE(itopo_type)

    CASE(topo_aster)
      CALL check_netcdf(nf90_open(path = topo_file_1, mode = nf90_nowrite, ncid = ncid))
      CALL check_netcdf(nf90_get_att(ncid, 3, "_FillValue", undef_topo), __FILE__, __LINE__)
      CALL check_netcdf(nf90_close(ncid))

    CASE(topo_gl)
      CALL check_netcdf(nf90_open(path = topo_file_1 , mode = nf90_nowrite, ncid = ncid))
      CALL check_netcdf(nf90_get_att(ncid, 1, "_FillValue", undef_topo), __FILE__, __LINE__)
      CALL check_netcdf(nf90_close(ncid))

    END SELECT

  END SUBROUTINE get_fill_value



  SUBROUTINE get_varname(topo_file_1,varname)
    IMPLICIT NONE
    SAVE
    CHARACTER(len=*), INTENT(IN) :: topo_file_1
    CHARACTER(len=*), INTENT(OUT)   :: varname
    INTEGER(i4)               :: ncid, type, ndims
    INTEGER(i4)               :: dimids(2)

#ifdef DEBUG
    logical :: lexist
    print*,"get_varname: ", trim(topo_file_1)
    inquire (file=trim(topo_file_1), exist=lexist)
    print*,"get_varname: file exists? ", lexist
#endif

    SELECT CASE(itopo_type)

    CASE(topo_aster)
      CALL check_netcdf(nf90_open(path = trim(topo_file_1), mode = nf90_nowrite, ncid = ncid))
      CALL check_netcdf(nf90_inquire_variable(ncid,3,varname,type,ndims,dimids))
      CALL check_netcdf(nf90_close(ncid))
    CASE(topo_gl)
      CALL check_netcdf(nf90_open(path = trim(topo_file_1), mode = nf90_nowrite, ncid = ncid), __FILE__, __LINE__)
      CALL check_netcdf(nf90_inquire_variable(ncid,1,varname,type,ndims,dimids), __FILE__, __LINE__)
      CALL check_netcdf(nf90_close(ncid), __FILE__, __LINE__)
    END SELECT

  END SUBROUTINE get_varname

  SUBROUTINE deallocate_topo_fields()

    IMPLICIT NONE

    INTEGER :: errorcode


    DEALLOCATE (topo_tiles_grid, STAT = errorcode)
    IF (errorcode /= 0) CALL abort_extpar('Cant deallocate the vector topo_tiles_grid')
    DEALLOCATE (tiles_lon_min, STAT = errorcode)
    IF (errorcode /= 0) CALL abort_extpar('Cant deallocate the vector tiles_lon_min')
    DEALLOCATE (tiles_lon_max, STAT = errorcode)
    IF (errorcode /= 0) CALL abort_extpar('Cant deallocate the vector tiles_lon_max')
    DEALLOCATE (tiles_lat_min, STAT = errorcode)
    IF (errorcode /= 0) CALL abort_extpar('Cant deallocate the vector tiles_lat_min')
    DEALLOCATE (tiles_lat_max, STAT = errorcode)
    IF (errorcode /= 0) CALL abort_extpar('Cant deallocate the vector tiles_lat_max')
    DEALLOCATE (hh_topo, STAT = errorcode)
    IF (errorcode /= 0) CALL abort_extpar('Cant deallocate the vector hh_topo')
    DEALLOCATE (hh_topo_max, STAT = errorcode)
    IF (errorcode /= 0) CALL abort_extpar('Cant deallocate the vector hh_topo_max')
    DEALLOCATE (hh_topo_min, STAT = errorcode)
    IF (errorcode /= 0) CALL abort_extpar('Cant deallocate the vector hh_topo_min')
    DEALLOCATE (stdh_topo, STAT = errorcode)
    IF (errorcode /= 0) CALL abort_extpar('Cant deallocate the vector stdh_topo')
    DEALLOCATE (theta_topo, STAT = errorcode)
    IF (errorcode /= 0) CALL abort_extpar('Cant deallocate the vector theta_topo')
    DEALLOCATE (aniso_topo, STAT = errorcode)
    IF (errorcode /= 0) CALL abort_extpar('Cant deallocate the vector aniso_topo')
    DEALLOCATE (slope_topo, STAT = errorcode)
    IF (errorcode /= 0) CALL abort_extpar('Cant deallocate the vector slope_topo')
    DEALLOCATE (fr_land_topo, STAT = errorcode)
    IF (errorcode /= 0) CALL abort_extpar('Cant deallocate the vector fr_land_topo')
    DEALLOCATE (z0_topo, STAT = errorcode)
    IF (errorcode /= 0) CALL abort_extpar('Cant deallocate the vector z0_topo')
    DEALLOCATE (slope_asp_topo, STAT = errorcode)
    IF (errorcode /= 0) CALL abort_extpar('Cant deallocate the vector slope_asp_topo')
    DEALLOCATE (slope_ang_topo, STAT = errorcode)
    IF (errorcode /= 0) CALL abort_extpar('Cant deallocate the vector slope_ang_topo')
    DEALLOCATE (horizon_topo, STAT = errorcode)
    IF (errorcode /= 0) CALL abort_extpar('Cant deallocate the vector horizon_topo')
    DEALLOCATE (skyview_topo, STAT = errorcode)
    IF (errorcode /= 0) CALL abort_extpar('Cant deallocate the vector skyview_topo')

  END SUBROUTINE deallocate_topo_fields

END MODULE mo_topo_data
