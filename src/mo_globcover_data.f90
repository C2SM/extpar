!+ Fortran Module with data fields for the globcover data
!
! Current Code Owner: DWD, Hermann Asensio
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_3         2011/04/19 Hermann Asensio
!  Initial release
! V2_0_3       2014/09/17 Burkhardt Rockel
!  Added use of directory information to access raw data files
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!> Fortran Module with data fields for the Globcover data
!> \author Hermann Asensio
!!
MODULE mo_globcover_data
                                                        
  USE mo_logging
  USE mo_kind,                  ONLY: wp,i4


  USE mo_grid_structures,       ONLY: reg_lonlat_grid

  USE mo_io_utilities,          ONLY: check_netcdf, &
       &                              join_path
                               
  USE netcdf,                   ONLY:      &
       &                              nf90_open,              &
       &                              nf90_close,             &
       &                              nf90_inquire_dimension, &
       &                              nf90_inquire_dimension, &
       &                              nf90_inq_dimid,         &
       &                              nf90_inq_varid,         &
       &                              nf90_get_var,           &
       &                              nf90_nowrite

  USE mo_lu_tg_fields,          ONLY: fr_land_lu,       &
       &                              ice_lu,           &
       &                              z0_lu,            &
       &                              z0_tot,           &
       &                              root_lu,          &
       &                              plcov_mn_lu,      &
       &                              plcov_mx_lu,      &
       &                              lai_mn_lu,        &
       &                              lai_mx_lu,        &
       &                              rs_min_lu,        &
       &                              urban_lu,         &
       &                              for_d_lu,         &
       &                              for_e_lu,         &
       &                              skinc_lu,         &
       &                              emissivity_lu,    &
       &                              fr_ocean_lu,      &
       &                              lu_class_fraction,&
       &                              lu_class_npixel,  &
       &                              lu_tot_npixel


  IMPLICIT NONE

  PRIVATE

  PUBLIC :: globcover_grid,                &
            lon_globcover,                 &
            lat_globcover,                 &
            allocate_raw_globcover_fields, &
            ntiles_globcover,              &   ! number of tiles in GLOBCOVER
            ncolumn_tiles,                 &   ! number of columns in tile matrix
            nrow_tiles,                    &   ! number of rows in tile matrix
            len_lu_lon,                    &   ! Number of columns in land-use dataset
            len_lu_lat,                    &   ! Number of rows in land-use dataset
            max_tiles_lu,                  &   ! maximal possible number of tiles that can be read
            nc_tiles_lu,                   &   ! total number of columns in one GLOBCOVER tile
            lu_tiles_lon_min,              &   ! starting longitude of every GLOBCOVER tile
            lu_tiles_lon_max,              &   ! ending longitude of every GLOBCOVER tile
            lu_tiles_lat_min,              &   ! starting latitude of every GLOBCOVER tile
            lu_tiles_lat_max,              &   ! ending latitude of every GLOBCOVER tile
            lu_tiles_ncolumns,             &   ! number of columns (lonitude increments) in each GLOBCOVER tile
            lu_tiles_nrows,                &   ! number of rows (latitude increments) in each GLOBCOVER tile
            globcover_tiles_grid,          &
            fill_globcover_data,           &   ! subroutine (intent(in) and intent(out))
            allocate_globcover_data,       &
            deallocate_globcover_fields

  TYPE(reg_lonlat_grid) :: globcover_grid !< structure with defenition of the raw data grid for the whole globcover dataset

  REAL (KIND=wp),  ALLOCATABLE         :: lon_globcover(:), &  !< longitude of globcover raw data
       &                                 lat_globcover(:) !< latitude of globcover raw data

  SAVE

  TYPE(reg_lonlat_grid), ALLOCATABLE  :: globcover_tiles_grid(:)

  INTEGER(KIND=i4), ALLOCATABLE       :: lu_tiles_ncolumns(:), & 
       &                                 lu_tiles_nrows(:)
                                      
  INTEGER(KIND=i4)                    :: ntiles_globcover = 6, &   ! number of GLOBCOVER tiles
       &                                 ncolumn_tiles = 3, &   ! number of used tiles E-W
       &                                 nrow_tiles = 2, &      ! number of used tiles N-S
       &                                 len_lu_lon, & 
       &                                 len_lu_lat, & 
       &                                 nc_tiles_lu
                                      
  INTEGER(KIND=i4), PARAMETER         :: max_tiles_lu = 1000
                                      
  REAL(KIND=wp), ALLOCATABLE          :: lu_tiles_lon_min(:), & 
      &                                  lu_tiles_lon_max(:), & 
      &                                  lu_tiles_lat_min(:), & 
      &                                  lu_tiles_lat_max(:)

  CONTAINS

  !> allocate raw data fields
  SUBROUTINE allocate_raw_globcover_fields(nrows,ncolumns)
    
    IMPLICIT NONE
    
    INTEGER (KIND=i4), INTENT(IN) :: nrows, &  !< number of rows
         &                           ncolumns !< number of columns

    INTEGER(KIND=i4)              :: errorcode !< error status variable

    CALL logging%info('Enter routine: allocate_raw_globcover_fields')

    ALLOCATE (lon_globcover(1:ncolumns), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array lon_globcover',__FILE__,__LINE__)
    lon_globcover = 0.0

    ALLOCATE (lat_globcover(1:nrows), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array lat_globcover',__FILE__,__LINE__)
    lat_globcover = 0.0

  END  SUBROUTINE allocate_raw_globcover_fields

  SUBROUTINE allocate_globcover_data(ntiles)   

    IMPLICIT NONE

    INTEGER, INTENT (IN) :: ntiles       ! number of tiles: 6 for GLOBCOVER
    INTEGER(KIND=i4)     :: errorcode

    CALL logging%info('Enter routine: allocate_globcover_data')

    ALLOCATE (lu_tiles_lon_min(1:ntiles), STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant allocate the vector lu_tiles_lon_min',__FILE__,__LINE__)
    ALLOCATE (lu_tiles_lon_max(1:ntiles), STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant allocate the vector lu_tiles_lon_max',__FILE__,__LINE__)
    ALLOCATE (lu_tiles_lat_min(1:ntiles), STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant allocate the vector lu_tiles_lat_min',__FILE__,__LINE__)
    ALLOCATE (lu_tiles_lat_max(1:ntiles), STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant allocate the vector lu_tiles_lat_max',__FILE__,__LINE__)
    
    ALLOCATE (lu_tiles_ncolumns(1:ntiles), STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant allocate the vector lu_tiles_ncolumns',__FILE__,__LINE__)
    ALLOCATE (lu_tiles_nrows(1:ntiles), STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant allocate the vector lu_tiles_nrows',__FILE__,__LINE__)

    ALLOCATE (globcover_tiles_grid(1:ntiles), STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant allocate the vector globcover_tiles_grid',__FILE__,__LINE__)
   
    lu_tiles_lon_min   = 0.0
    lu_tiles_lon_max   = 0.0
    lu_tiles_lat_min   = 0.0
    lu_tiles_lat_max   = 0.0
    lu_tiles_ncolumns  = 0
    lu_tiles_nrows     = 0
   
  END SUBROUTINE allocate_globcover_data

  SUBROUTINE fill_globcover_data(raw_data_lu_path,       &
       &                          raw_data_lu_filename,   &
       &                          lu_tiles_lon_min,       &
       &                          lu_tiles_lon_max,       &  ! the allocated vectors need to be filled with the respective value
       &                          lu_tiles_lat_min,       &
       &                          lu_tiles_lat_max,       &
       &                          nc_tiles_lu)

    IMPLICIT NONE

    SAVE

    CHARACTER (len=*),INTENT(IN) :: raw_data_lu_path, & 
         &                                     raw_data_lu_filename(:)

    REAL(KIND=wp), INTENT(OUT)   :: lu_tiles_lon_min(1:ntiles_globcover), &  
         &                          lu_tiles_lon_max(1:ntiles_globcover), & 
         &                          lu_tiles_lat_min(1:ntiles_globcover), & 
         &                          lu_tiles_lat_max(1:ntiles_globcover)

    INTEGER(KIND=i4), INTENT(OUT):: nc_tiles_lu

    INTEGER(KIND=i4)             :: i, ncid, &        ! i is a counter
         &                          dimID_lat, dimID_lon, varID_lat, varID_lon                  

    REAL(KIND=wp)                :: half_gridp ! distance of half a grid point as the grid point is centered on a GLOBCOVER pixel
     
    half_gridp = 0.001388888889

    DO i = 1,ntiles_globcover
      ! open GLOBCOVER file
      CALL check_netcdf(nf90_open(path =join_path(raw_data_lu_path,raw_data_lu_filename(i)), mode = nf90_nowrite, ncid = ncid))
      CALL check_netcdf(nf90_inq_dimid(ncid,"lon", dimID_lon))
      CALL check_netcdf(nf90_inq_dimid(ncid,"lat", dimID_lat))
      CALL check_netcdf(nf90_inquire_dimension(ncid,dimID_lon, len = lu_tiles_ncolumns(i)))          
      CALL check_netcdf(nf90_inquire_dimension(ncid,dimID_lat, len = lu_tiles_nrows(i))) 
      CALL check_netcdf(nf90_inq_varid(ncid, "lon", varID_lon))
      CALL check_netcdf(nf90_inq_varid(ncid, "lat", varID_lat))
      CALL check_netcdf(nf90_get_var(ncid, varID_lon, lu_tiles_lon_min(i), start = (/1/)))            
      ! reads in the first longitude value of tile i
      CALL check_netcdf(nf90_get_var(ncid, varID_lon, lu_tiles_lon_max(i), start = (/lu_tiles_ncolumns(i)/))) 
      ! reads in the last longitude value of tile i
      CALL check_netcdf(nf90_get_var(ncid, varID_lat, lu_tiles_lat_max(i), start = (/1/)))            
      ! reads in the first latitude value of tile i
      CALL check_netcdf(nf90_get_var(ncid, varID_lat, lu_tiles_lat_min(i), start = (/lu_tiles_nrows(i)/))) 
      ! reads in the last latitude value of tile i
      CALL check_netcdf(nf90_close(ncid))  ! the netcdf file is closed again
      lu_tiles_lon_min(i) = lu_tiles_lon_min(i) - half_gridp !< half of a grid point must be
      lu_tiles_lon_max(i) = lu_tiles_lon_max(i) + half_gridp !< added, as the GLOBCOVER data
      lu_tiles_lat_min(i) = lu_tiles_lat_min(i) - half_gridp !< is located at the pixel center
      lu_tiles_lat_max(i) = lu_tiles_lat_max(i) + half_gridp
      
      len_lu_lon=lu_tiles_ncolumns(i)
      len_lu_lat=lu_tiles_nrows(i)
    ENDDO

    nc_tiles_lu = lu_tiles_ncolumns(1)

  END SUBROUTINE fill_globcover_data

  SUBROUTINE  deallocate_globcover_fields()

    IMPLICIT NONE     

    INTEGER(KIND=i4) :: errorcode

    CALL logging%info('Enter routine: deallocate_globcover_data')

    DEALLOCATE (lu_tiles_lon_min, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector lu_tiles_lon_min',__FILE__,__LINE__)
    DEALLOCATE (lu_tiles_lon_max, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector lu_tiles_lon_max',__FILE__,__LINE__)
    DEALLOCATE (lu_tiles_lat_min, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector lu_tiles_lat_min',__FILE__,__LINE__)
    DEALLOCATE (lu_tiles_lat_max, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector lu_tiles_lat_max',__FILE__,__LINE__)

    DEALLOCATE (lu_tiles_ncolumns, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector lu_tiles_ncolumns',__FILE__,__LINE__)
    DEALLOCATE (lu_tiles_nrows, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector lu_tiles_nrows',__FILE__,__LINE__)

    DEALLOCATE (globcover_tiles_grid, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector globcover_tiles_grid',__FILE__,__LINE__)

    DEALLOCATE (lat_globcover, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector lat_globcover',__FILE__,__LINE__)
    DEALLOCATE (lon_globcover, STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector lon_globcover',__FILE__,__LINE__)

  END SUBROUTINE deallocate_globcover_fields

END MODULE mo_globcover_data
