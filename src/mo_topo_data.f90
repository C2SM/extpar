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

 USE mo_logging
 USE mo_kind,                   ONLY: wp, i4, i2
                                
 USE mo_grid_structures,        ONLY: reg_lonlat_grid
                               
 USE mo_io_utilities,           ONLY: check_netcdf

 USE mo_topo_tg_fields, ONLY:         fr_land_topo,  &
      &                               hh_topo,       &
      &                               hh_topo_max,   &
      &                               hh_topo_min,   &
      &                               sgsl,          &
      &                               stdh_topo,     &
      &                               theta_topo,    &
      &                               aniso_topo,    &
      &                               slope_topo,    &
      &                               z0_topo,       &
      &                               slope_asp_topo,&
      &                               slope_ang_topo,&
      &                               horizon_topo,  &
      &                               skyview_topo

 USE netcdf,                    ONLY:    &
      &                               nf90_open,              &
      &                               nf90_close,             &
      &                               nf90_inquire_dimension, &
      &                               nf90_inquire_variable,  &
      &                               nf90_inquire_dimension, &
      &                               nf90_inq_dimid,         &
      &                               nf90_inq_varid,         &
      &                               nf90_get_var,           &
      &                               nf90_get_att,           &
      &                               nf90_enotvar,           &
      &                               nf90_nowrite

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: ntiles,                  &   ! number of tiles in GLOBE / ASTER
       &    max_tiles,               &   ! maximal possible number of tiles that can be read
       &    nc_tot,                  &   ! total number of columns in GLOBE / ASTER data set
       &    nr_tot,                  &   ! total number of rows in GLOBE / ASTER data set
       &    nc_tile ,                &   ! total number of columns in one GLOBE / ASTER tile
       &    tiles_lon_min,           &   ! starting longitude of every GLOBE / ASTER tile
       &    tiles_lon_max,           &   ! ending longitude of every GLOBE / ASTER tile
       &    tiles_lat_min,           &   ! starting latitude of every GLOBE / ASTER tile
       &    tiles_lat_max,           &   ! ending latitude of every GLOBE / ASTER tile
       &    tiles_ncolumns,          &   ! number of columns (lonitude increments) in each GLOBE / ASTER tile
       &    tiles_nrows,             &   ! number of rows (latitude increments) in each GLOBE /ASTER tile
       &    topo_tiles_grid,         &
       &    topo_grid,               &
       &    raw_topo_line,           &
       &    h_tile_row,              &
       &    raw_topo_block,          &
       &    raw_sgsl_block,          &
       &    allocate_raw_topo_fields,&
       &    fill_topo_data,           &  ! subroutine (intent(in) and intent(out))
       &    num_tiles,                &  ! integer function
       &    allocate_topo_data,       &  ! subroutine (only intent(in))
       &    get_fill_value,           &
       &    get_fill_value_sgsl,       &
       &    get_varname,              &
       &    get_varname_sgsl,         &
       &    undef_topo,               &
       &    varname,                  &
       &    itopo_type,               &
       &    topo_gl,                  &
       &    topo_aster,               &
       &    aster_lat_min,            &
       &    aster_lon_min,            &
       &    aster_lat_max,            &
       &    aster_lon_max,            &
       &    ntiles_row,               &
       &    ntiles_column,            &
       &    lradtopo,                 &
       &    nhori,                    &
       &    deallocate_topo_fields

  SAVE

  TYPE(reg_lonlat_grid), ALLOCATABLE  :: topo_tiles_grid(:)
  TYPE(reg_lonlat_grid)               :: topo_grid

  INTEGER(KIND=i4), ALLOCATABLE :: tiles_ncolumns(:), &
       &                           tiles_nrows(:), &
       &                           h_tile_row(:)

  INTEGER(KIND=i4)              :: ntiles, &
       &                           nc_tot, &
       &                           nr_tot, &
       &                           nc_tile, &
       &                           undef_topo, &
       &                           itopo_type, &
       &                           ntiles_row, &
       &                           ntiles_column, &
       &                           nhori

  INTEGER(KIND=i4), PARAMETER   :: topo_gl = 1, &
       &                           topo_aster = 2, &
       &                           topo_merit = 3, &
       &                           max_tiles = 1000

  REAL(KIND=wp), ALLOCATABLE    :: tiles_lon_min(:), &
       &                           tiles_lon_max(:), &
       &                           tiles_lat_min(:), &
       &                           tiles_lat_max(:), &
       &                           raw_topo_line(:), &
       &                           raw_topo_block(:,:), &
       &                           raw_sgsl_block(:,:)

  REAL(KIND=wp)::                  aster_lat_min, &
       &                           aster_lat_max, &
       &                           aster_lon_min, &
       &                           aster_lon_max

  REAL(KIND=wp)::                  merit_lat_min, &
       &                           merit_lat_max, &
       &                           merit_lon_min, &
       &                           merit_lon_max

  LOGICAL                       :: lradtopo

  CHARACTER(LEN=80)             :: varname


  CONTAINS

  SUBROUTINE num_tiles(columns, rows, ntiles) ! it gives the value of the number of tiles depending
    SAVE
    INTEGER(KIND=i4), INTENT(IN) :: columns, &
         &                          rows
    INTEGER, INTENT(OUT)         :: ntiles           ! if the user chooses GLOBE, ASTER, MERIT

    ntiles_column = columns
    ntiles_row    = rows
    ntiles = ntiles_column * ntiles_row

  END SUBROUTINE num_tiles

  SUBROUTINE allocate_topo_data(ntiles)
    ! As it is unknown so far whether GLOBE or ASTER is chosen all parameters must be allocated in a second step.

    IMPLICIT NONE

    INTEGER, INTENT (IN) :: ntiles       ! number of tiles: 36 for ASTER and 16 for GLOBE
    INTEGER              :: errorcode
    
    CALL logging%info('Enter routine: allocate_topo_data')
   
    ALLOCATE (tiles_lon_min(1:ntiles), STAT = errorcode)
         IF (errorcode.NE.0) CALL logging%error('Cant allocate the vector tiles_lon_min',__FILE__,__LINE__)
    ALLOCATE (tiles_lon_max(1:ntiles), STAT = errorcode)
         IF (errorcode.NE.0) CALL logging%error('Cant allocate the vector tiles_lon_max',__FILE__,__LINE__)
    ALLOCATE (tiles_lat_min(1:ntiles), STAT = errorcode)
         IF (errorcode.NE.0) CALL logging%error('Cant allocate the vector tiles_lat_min',__FILE__,__LINE__)
    ALLOCATE (tiles_lat_max(1:ntiles), STAT = errorcode)
         IF (errorcode.NE.0) CALL logging%error('Cant allocate the vector tiles_lat_max',__FILE__,__LINE__)
    
    ALLOCATE (tiles_ncolumns(1:ntiles), STAT = errorcode)
         IF (errorcode.NE.0) CALL logging%error('Cant allocate the vector tiles_ncolumns',__FILE__,__LINE__)
    ALLOCATE (tiles_nrows(1:ntiles), STAT = errorcode)
         IF (errorcode.NE.0) CALL logging%error('Cant allocate the vector tiles_nrows',__FILE__,__LINE__)
    
    ALLOCATE (topo_tiles_grid(1:ntiles), STAT = errorcode)
         IF (errorcode.NE.0) CALL logging%error('Cant allocate the vector topo_tiles_grid',__FILE__,__LINE__)
    
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
  
    merit_lat_min   = 0.0
    merit_lat_max   = 0.0
    merit_lon_min   = 0.0
    merit_lon_max   = 0.0


  END SUBROUTINE allocate_topo_data

  SUBROUTINE fill_topo_data(raw_data_orography_path,topo_files,  &
                                                  tiles_lon_min,  &
                                                  tiles_lon_max,  &
                                                  tiles_lat_min,  &
                                                  tiles_lat_max,  &
                                                  nc_tot,         &
                                                  nr_tot,         &
                                                  nc_tile)
   IMPLICIT NONE

   SAVE
   CHARACTER (len=*),INTENT(IN) :: raw_data_orography_path, &
        &                          topo_files(1:ntiles)

   REAL(KIND=wp), INTENT(OUT)   :: tiles_lon_min(1:ntiles), &
        &                          tiles_lon_max(1:ntiles), &
        &                          tiles_lat_min(1:ntiles), &
        &                          tiles_lat_max(1:ntiles)

   INTEGER(KIND=i4), INTENT(OUT):: nc_tot, nr_tot, nc_tile

   INTEGER(KIND=i4)             :: i, errorcode, ncid, &
        &                          dimID_lat, dimID_lon, varID_lat, varID_lon
                               
   REAL(KIND=wp)                :: half_gridp! distance of half a grid point as the grid point is centered on a GLOBE / ASTER pixel

   CALL logging%info('Enter routine: fill_topo_data')

   SELECT CASE (itopo_type)                ! Also topo could additionally be used for SELECT CASE (must first be read in)
     CASE(topo_aster)                                         ! ASTER topography, as it has 36 tiles at the moment.
       CALL logging%info('ASTER is used as topography')
       half_gridp = 1./(3600.*2.)           ! the resolution of the ASTER data is 1./3600. degrees as it is half a grid point
                                            ! it is additionally divided by 2
     CASE (topo_gl)                                           ! GLOBE topography is composed of 16 tiles
       CALL logging%info('GLOBE is used as topography')
       half_gridp = 1./(120.*2.)                              ! GLOBE resolution is 1./120. degrees (30 arc-seconds)

     CASE(topo_merit)                                         ! ASTER topography, as it has 36 tiles at the moment.
       CALL logging%info('MERIT is used as topography')
       half_gridp = 1./(1200.*2.)           ! the resolution of the MERIT data is 1./1200. degrees as it is half a grid point
                                            ! it is additionally divided by 2
   END SELECT

   DO i = 1,ntiles
     CALL check_netcdf(nf90_open(path =TRIM(raw_data_orography_path)//TRIM(topo_files(i)), mode = nf90_nowrite, ncid = ncid))
     CALL check_netcdf(nf90_inq_dimid(ncid,"lon", dimID_lon))
     CALL check_netcdf(nf90_inq_dimid(ncid,"lat", dimID_lat))
     CALL check_netcdf(nf90_inquire_dimension(ncid,dimID_lon, len = tiles_ncolumns(i)))
     CALL check_netcdf(nf90_inquire_dimension(ncid,dimID_lat, len = tiles_nrows(i)))
     CALL check_netcdf(nf90_inq_varid(ncid, "lon", varID_lon))
     CALL check_netcdf(nf90_inq_varid(ncid, "lat", varID_lat))
     CALL check_netcdf(nf90_get_var(ncid, varID_lon, tiles_lon_min(i), start = (/1/)))
     ! reads in the first longitude value of tile i
     CALL check_netcdf(nf90_get_var(ncid, varID_lon, tiles_lon_max(i), start = (/tiles_ncolumns(i)/)))
     ! reads in the last longitude value of tile i
     CALL check_netcdf(nf90_get_var(ncid, varID_lat, tiles_lat_max(i), start = (/1/)))
     ! reads in the first latitude value of tile i
     CALL check_netcdf(nf90_get_var(ncid, varID_lat, tiles_lat_min(i), start = (/tiles_nrows(i)/)))
     ! reads in the last latitude value of tile i
     CALL check_netcdf(nf90_close(ncid))
     ! the netcdf file is closed again
     tiles_lon_min(i) = REAL(NINT(tiles_lon_min(i) - half_gridp)) !< half of a grid point must be
     tiles_lon_max(i) = REAL(NINT(tiles_lon_max(i) + half_gridp)) !< added, as the ASTER/GLOBE/MERIT data
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

     CASE(topo_merit)
       merit_lat_min = MINVAL(tiles_lat_min)
       merit_lat_max = MAXVAL(tiles_lat_max)
       merit_lon_min = MINVAL(tiles_lon_min)
       merit_lon_max = MAXVAL(tiles_lon_max)
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
   END SELECT

   ALLOCATE (raw_topo_line(1:nc_tot), STAT = errorcode)
         IF (errorcode.NE.0) CALL logging%error('Cant allocate the vector raw_topo_line',__FILE__,__LINE__)
   ALLOCATE (h_tile_row(1:nc_tile), STAT = errorcode)
         IF (errorcode.NE.0) CALL logging%error('Cant allocate the vector h_tile_row',__FILE__,__LINE__)

   CALL logging%info('Exit routine: fill_topo_data')

  END SUBROUTINE fill_topo_data

  SUBROUTINE allocate_raw_topo_fields(nrows,ncolumns,lcompute_sgsl)

    INTEGER(KIND=i4), INTENT(IN)  :: nrows, ncolumns

    LOGICAL, INTENT(IN)           :: lcompute_sgsl

    INTEGER(KIND=i4)              :: errorcode


    ALLOCATE (raw_topo_block(1:ncolumns, 1:nrows), STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant allocate the array raw_topo_block',__FILE__,__LINE__)
    raw_topo_block = 0.0

    IF (lcompute_sgsl) THEN
      ALLOCATE (raw_sgsl_block(1:ncolumns, 1:nrows), STAT = errorcode)
      IF (errorcode.NE.0) CALL logging%error('Cant allocate the array raw_sgsl_block',__FILE__,__LINE__)
      raw_sgsl_block = 0.0
    ENDIF

  END SUBROUTINE allocate_raw_topo_fields

  SUBROUTINE get_fill_value(topo_file_1, undef_topo)

    CHARACTER (len=*), INTENT(in)  :: topo_file_1
    INTEGER(KIND=i4), INTENT(out)  :: undef_topo

    INTEGER(KIND=i4)               :: ncid, varid, status

    CALL check_netcdf(nf90_open(path = topo_file_1, mode = nf90_nowrite, ncid = ncid))
    status = nf90_inq_varid(ncid, "altitude", varid)
    IF (status == NF90_ENOTVAR) THEN
      status = nf90_inq_varid(ncid, "Z", varid)      
      IF (status == NF90_ENOTVAR) THEN
        CALL logging%error('Could not find "altitude" or "Z" in topography file '//TRIM(topo_file_1),__FILE__,__LINE__)
      ELSE
        CALL check_netcdf(status, __FILE__, __LINE__)      
      ENDIF
    ELSE
      CALL check_netcdf(status, __FILE__, __LINE__)      
    ENDIF
    CALL check_netcdf(nf90_get_att(ncid, varid, "_FillValue", undef_topo), __FILE__, __LINE__)
    CALL check_netcdf(nf90_close(ncid))

  END SUBROUTINE get_fill_value

  SUBROUTINE get_varname(topo_file_1,varname)

    CHARACTER(len=*), INTENT(IN) :: topo_file_1
    CHARACTER(len=*), INTENT(OUT):: varname
    INTEGER(KIND=i4)             :: ncid, type, ndims, &
        &                           dimids(2)

    SELECT CASE(itopo_type)
      CASE(topo_aster)
        CALL check_netcdf(nf90_open(path = trim(topo_file_1), mode = nf90_nowrite, ncid = ncid))
        CALL check_netcdf(nf90_inquire_variable(ncid,3,varname,type,ndims,dimids))
        CALL check_netcdf(nf90_close(ncid))
      CASE(topo_gl)
        CALL check_netcdf(nf90_open(path = trim(topo_file_1), mode = nf90_nowrite, ncid = ncid), __FILE__, __LINE__)
        CALL check_netcdf(nf90_inquire_variable(ncid,1,varname,type,ndims,dimids), __FILE__, __LINE__)
        CALL check_netcdf(nf90_close(ncid), __FILE__, __LINE__)
      CASE(topo_merit)
        CALL check_netcdf(nf90_open(path = trim(topo_file_1), mode = nf90_nowrite, ncid = ncid))
        CALL check_netcdf(nf90_inquire_variable(ncid,3,varname,type,ndims,dimids))
        CALL check_netcdf(nf90_close(ncid))
    END SELECT

  END SUBROUTINE get_varname

  SUBROUTINE deallocate_topo_fields(lcompute_sgsl)

    IMPLICIT NONE
    
    LOGICAL, INTENT(IN):: lcompute_sgsl

    INTEGER(KIND=i4)   :: errorcode

    CALL logging%info('Enter routine: deallocate_topo_fields')

    DEALLOCATE (topo_tiles_grid, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector topo_tiles_grid',__FILE__,__LINE__)
    DEALLOCATE (tiles_lon_min, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector tiles_lon_min',__FILE__,__LINE__)
    DEALLOCATE (tiles_lon_max, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector tiles_lon_max',__FILE__,__LINE__)
    DEALLOCATE (tiles_lat_min, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector tiles_lat_min',__FILE__,__LINE__)
    DEALLOCATE (tiles_lat_max, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector tiles_lat_max',__FILE__,__LINE__)
    DEALLOCATE (hh_topo, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector hh_topo',__FILE__,__LINE__)
    DEALLOCATE (hh_topo_max, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector hh_topo_max',__FILE__,__LINE__)
    DEALLOCATE (hh_topo_min, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector hh_topo_min',__FILE__,__LINE__)
    DEALLOCATE (stdh_topo, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector stdh_topo',__FILE__,__LINE__)
    DEALLOCATE (theta_topo, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector theta_topo',__FILE__,__LINE__)
    DEALLOCATE (aniso_topo, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector aniso_topo',__FILE__,__LINE__)
    DEALLOCATE (slope_topo, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector slope_topo',__FILE__,__LINE__)
    DEALLOCATE (fr_land_topo, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector fr_land_topo',__FILE__,__LINE__)
    DEALLOCATE (z0_topo, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector z0_topo',__FILE__,__LINE__)
    DEALLOCATE (slope_asp_topo, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector slope_asp_topo',__FILE__,__LINE__)
    DEALLOCATE (slope_ang_topo, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector slope_ang_topo',__FILE__,__LINE__)
    DEALLOCATE (horizon_topo, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector horizon_topo',__FILE__,__LINE__)
    DEALLOCATE (skyview_topo, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector skyview_topo',__FILE__,__LINE__)

    IF (lcompute_sgsl) THEN
      DEALLOCATE (sgsl, STAT = errorcode)
      IF (errorcode.NE.0) CALL logging%error('Cant deallocate the array sgsl',__FILE__,__LINE__)
    ENDIF

  END SUBROUTINE deallocate_topo_fields

  SUBROUTINE get_fill_value_sgsl(sgsl_file_1,undef_sgsl)

    IMPLICIT NONE

    SAVE

    CHARACTER (len=*), INTENT(IN) :: sgsl_file_1     
    REAL(KIND=wp), INTENT(OUT)    :: undef_sgsl

    INTEGER(KIND=i2)              :: fillval
    INTEGER(KIND=i4)              :: ncid
    REAL(KIND=wp)                 :: scale_factor

    SELECT CASE(itopo_type)
    
      CASE(topo_aster)
        CALL check_netcdf(nf90_open(path = sgsl_file_1, mode = nf90_nowrite, ncid = ncid))
        CALL check_netcdf(nf90_get_att(ncid, 3, "_FillValue", fillval))
        CALL check_netcdf(nf90_get_att(ncid, 3, "scale_factor", scale_factor))
        undef_sgsl = fillval * scale_factor
        CALL check_netcdf(nf90_close(ncid))

      CASE(topo_gl)
        CALL check_netcdf(nf90_open(path = sgsl_file_1 , mode = nf90_nowrite, ncid = ncid))
        CALL check_netcdf(nf90_get_att(ncid, 3, "_FillValue", fillval))
        CALL check_netcdf(nf90_get_att(ncid, 3, "scale_factor", scale_factor))
        undef_sgsl = fillval * scale_factor
        CALL check_netcdf(nf90_close(ncid))

      CASE(topo_merit)
        CALL check_netcdf(nf90_open(path = sgsl_file_1, mode = nf90_nowrite, ncid = ncid))
        CALL check_netcdf(nf90_get_att(ncid, 3, "_FillValue", fillval))
        CALL check_netcdf(nf90_get_att(ncid, 3, "scale_factor", scale_factor))
        undef_sgsl = fillval * scale_factor
        CALL check_netcdf(nf90_close(ncid))

    END SELECT

  END SUBROUTINE get_fill_value_sgsl

  SUBROUTINE get_varname_sgsl(sgsl_file_1,varname)
  
  IMPLICIT NONE

  SAVE
  CHARACTER (len=*), INTENT(IN) :: sgsl_file_1     
  CHARACTER(LEN=*),INTENT(OUT)   :: varname
  INTEGER(KIND=i4)               :: ncid, type, ndims
  INTEGER(KIND=i4)               :: dimids(2)

  SELECT CASE(itopo_type)
  
   CASE(topo_aster)
     CALL check_netcdf(nf90_open(path = sgsl_file_1, mode = nf90_nowrite, ncid = ncid))
     CALL check_netcdf(nf90_inquire_variable(ncid,3,varname,type,ndims,dimids))
     CALL check_netcdf(nf90_close(ncid))

   CASE(topo_gl)
   
     CALL check_netcdf(nf90_open(path = sgsl_file_1, mode = nf90_nowrite, ncid = ncid))
  
     CALL check_netcdf(nf90_inquire_variable(ncid,3,varname,type,ndims,dimids))
  
     CALL check_netcdf(nf90_close(ncid))
     varname = TRIM(varname)
  END SELECT

   CASE(topo_merit)
     CALL check_netcdf(nf90_open(path = sgsl_file_1, mode = nf90_nowrite, ncid = ncid))
     CALL check_netcdf(nf90_inquire_variable(ncid,3,varname,type,ndims,dimids))
     CALL check_netcdf(nf90_close(ncid))

  END SUBROUTINE get_varname_sgsl

END MODULE mo_topo_data
