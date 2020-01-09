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

!> kind parameters are defined in MODULE data_parameters
USE mo_kind, ONLY: wp, &
                   i4, &
                   i8

!> abort_extpar defined in MODULE utilities_extpar
USE mo_utilities_extpar, ONLY: abort_extpar

USE mo_grid_structures, ONLY: reg_lonlat_grid

! >mes

USE mo_io_utilities,       ONLY:  check_netcdf

USE mo_io_units,           ONLY:  filename_max

USE netcdf,       ONLY :     &
     nf90_open,              &
     nf90_close,             &
     nf90_inquire,           &
     nf90_inquire_dimension, &
     nf90_inquire_variable,  &
     nf90_inq_attname,       &
     nf90_inquire_attribute, &
     nf90_get_att,           &
     nf90_inquire_dimension, &
     nf90_inq_dimid,         &
     nf90_inq_varid,         &
     nf90_get_var,           &
     nf90_noerr,             &
     nf90_strerror,          &
     nf90_create,            &
     nf90_def_dim,           &
     nf90_def_var,           &
     nf90_enddef,            &
     nf90_redef,             &
     nf90_put_att,           &
     nf90_put_var,           &
     NF90_CHAR,              &
     NF90_DOUBLE,            &
     NF90_FLOAT,             &
     NF90_INT,               &
     NF90_BYTE,              &
     NF90_SHORT,             &
     NF90_GLOBAL,            &
     NF90_UNLIMITED,         &
     NF90_CLOBBER,           &
     NF90_NOWRITE

! <mes

IMPLICIT NONE

PRIVATE

PUBLIC :: globcover_grid,                &
          lon_globcover,                 &
          lat_globcover,                 &
          allocate_raw_globcover_fields, &
! >mes
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
          deallocate_landuse_data

TYPE(reg_lonlat_grid) :: globcover_grid !< structure with defenition of the raw data grid for the whole globcover dataset

REAL (KIND=wp), ALLOCATABLE    :: lon_globcover(:) !< longitude of globcover raw data
REAL (KIND=wp), ALLOCATABLE    :: lat_globcover(:) !< latitude of globcover raw data

SAVE
TYPE(reg_lonlat_grid), ALLOCATABLE  :: globcover_tiles_grid(:)

INTEGER(KIND=i4), ALLOCATABLE :: lu_tiles_ncolumns(:)
INTEGER(KIND=i4), ALLOCATABLE :: lu_tiles_nrows(:)
INTEGER(KIND=i4)              :: ntiles_globcover = 6  ! number of GLOBCOVER tiles
INTEGER(KIND=i4)              :: ncolumn_tiles = 3  ! number of used tiles E-W
INTEGER(KIND=i4)              :: nrow_tiles = 2     ! number of used tiles N-S
INTEGER(KIND=i4)              :: len_lu_lon
INTEGER(KIND=i4)              :: len_lu_lat
INTEGER, PARAMETER            :: max_tiles_lu = 1000
INTEGER(KIND=i4)              :: nc_tiles_lu

REAL(KIND=wp), ALLOCATABLE    :: lu_tiles_lon_min(:)
REAL(KIND=wp), ALLOCATABLE    :: lu_tiles_lon_max(:)
REAL(KIND=wp), ALLOCATABLE    :: lu_tiles_lat_min(:)
REAL(KIND=wp), ALLOCATABLE    :: lu_tiles_lat_max(:)
REAL(KIND=wp), ALLOCATABLE    :: raw_lu_line(:)
REAL(KIND=wp), ALLOCATABLE    :: raw_lu_block(:,:)
! <mes

CONTAINS

  !> allocate raw data fields
  SUBROUTINE allocate_raw_globcover_fields(nrows,ncolumns)
  IMPLICIT NONE
  INTEGER (KIND=i8), INTENT(IN) :: nrows !< number of rows
  INTEGER (KIND=i8), INTENT(IN) :: ncolumns !< number of columns

  INTEGER :: errorcode !< error status variable

    ALLOCATE (lon_globcover(1:ncolumns), STAT=errorcode)
        IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array lon_globcover')
    lon_globcover = 0.0

     ALLOCATE (lat_globcover(1:nrows), STAT=errorcode)
        IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array lat_globcover')
    lat_globcover = 0.0

  END  SUBROUTINE allocate_raw_globcover_fields

  !----------------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------
! >mes

     SUBROUTINE allocate_globcover_data(ntiles)   
       IMPLICIT NONE
       INTEGER, INTENT (IN) :: ntiles       ! number of tiles: 6 for GLOBCOVER
       INTEGER :: errorcode

       ALLOCATE (lu_tiles_lon_min(1:ntiles), STAT = errorcode)
       IF (errorcode.NE.0) CALL abort_extpar('Cant allocate the vector lu_tiles_lon_min')
       ALLOCATE (lu_tiles_lon_max(1:ntiles), STAT = errorcode)
       IF (errorcode.NE.0) CALL abort_extpar('Cant allocate the vector lu_tiles_lon_max')
       ALLOCATE (lu_tiles_lat_min(1:ntiles), STAT = errorcode)
       IF (errorcode.NE.0) CALL abort_extpar('Cant allocate the vector lu_tiles_lat_min')
       ALLOCATE (lu_tiles_lat_max(1:ntiles), STAT = errorcode)
       IF (errorcode.NE.0) CALL abort_extpar('Cant allocate the vector lu_tiles_lat_max')
       
       ALLOCATE (lu_tiles_ncolumns(1:ntiles), STAT = errorcode)
       IF (errorcode.NE.0) CALL abort_extpar('Cant allocate the vector lu_tiles_ncolumns')
       ALLOCATE (lu_tiles_nrows(1:ntiles), STAT = errorcode)
       IF (errorcode.NE.0) CALL abort_extpar('Cant allocate the vector lu_tiles_nrows')

       ALLOCATE (globcover_tiles_grid(1:ntiles), STAT = errorcode)
       IF (errorcode.NE.0) CALL abort_extpar('Cant allocate the vector globcover_tiles_grid')
   
   lu_tiles_lon_min   = 0.0
   lu_tiles_lon_max   = 0.0
   lu_tiles_lat_min   = 0.0
   lu_tiles_lat_max   = 0.0
   lu_tiles_ncolumns  = 0
   lu_tiles_nrows     = 0
   
   END SUBROUTINE allocate_globcover_data

! <mes

  !----------------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------
! >mes

   SUBROUTINE fill_globcover_data(raw_data_lu_path,       &
                                  raw_data_lu_filename,   &
                                  lu_tiles_lon_min,       &
                                  lu_tiles_lon_max,       &  ! the allocated vectors need to be filled with the respective value
                                  lu_tiles_lat_min,       &
                                  lu_tiles_lat_max,       &
                                  nc_tiles_lu)
   IMPLICIT NONE
   SAVE
   CHARACTER (len=filename_max),INTENT(IN) :: raw_data_lu_path
   CHARACTER (len=filename_max),INTENT(IN) :: raw_data_lu_filename(:)
   REAL(KIND=wp), INTENT(OUT)   :: lu_tiles_lon_min(1:ntiles_globcover) 
   REAL(KIND=wp), INTENT(OUT)   :: lu_tiles_lon_max(1:ntiles_globcover)
   REAL(KIND=wp), INTENT(OUT)   :: lu_tiles_lat_min(1:ntiles_globcover)
   REAL(KIND=wp), INTENT(OUT)   :: lu_tiles_lat_max(1:ntiles_globcover)
   INTEGER(KIND=i4), INTENT(OUT):: nc_tiles_lu
   CHARACTER(len=2)    :: num
   CHARACTER(len=80)   :: path
   INTEGER(KIND=i4)    :: i, errorcode        ! i is a counter, errorcode is used to check if allocation was successful
   INTEGER(KIND=i4)    :: ncid
   INTEGER(KIND=i4)    :: dimID_lat, dimID_lon, varID_lat, varID_lon                  
   REAL(KIND=wp)       :: half_gridp          ! distance of half a grid point as the grid point is centered on a GLOBCOVER pixel
    
   half_gridp = 0.001388888889
!   print*, half_gridp

     DO i = 1,ntiles_globcover
       ! open GLOBCOVER file
       CALL check_netcdf(nf90_open(path =TRIM(raw_data_lu_path)//TRIM(raw_data_lu_filename(i)), mode = nf90_nowrite, ncid = ncid))
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
       print*, 'GLOBCOVER TILE ',ntiles_globcover,': NLON,NLAT ',len_lu_lon,len_lu_lat
     ENDDO

     nc_tiles_lu = lu_tiles_ncolumns(1)


   END SUBROUTINE fill_globcover_data

! <mes


   SUBROUTINE  deallocate_landuse_data()

     USE mo_lu_tg_fields, ONLY: fr_land_lu,       &
                       &        ice_lu,           &
                       &        z0_lu,            &
                       &        z0_tot,           &
                       &        root_lu,          &
                       &        plcov_mn_lu,      &
                       &        plcov_mx_lu,      &
                       &        lai_mn_lu,        &
                       &        lai_mx_lu,        &
                       &        rs_min_lu,        &
                       &        urban_lu,         &
                       &        for_d_lu,         &
                       &        for_e_lu,         &
                       &        skinc_lu,         &
                       &        emissivity_lu,    &
                       &        fr_ocean_lu,      &
                       &        lu_class_fraction,&
                       &        lu_class_npixel,  &
                       &        lu_tot_npixel
     
     IMPLICIT NONE     

     INTEGER :: errorcode

       DEALLOCATE (lu_tiles_lon_min, STAT = errorcode)
       IF (errorcode.NE.0) CALL abort_extpar('Cant deallocate the vector lu_tiles_lon_min')
       DEALLOCATE (lu_tiles_lon_max, STAT = errorcode)
       IF (errorcode.NE.0) CALL abort_extpar('Cant deallocate the vector lu_tiles_lon_max')
       DEALLOCATE (lu_tiles_lat_min, STAT = errorcode)
       IF (errorcode.NE.0) CALL abort_extpar('Cant deallocate the vector lu_tiles_lat_min')
       DEALLOCATE (lu_tiles_lat_max, STAT = errorcode)
       IF (errorcode.NE.0) CALL abort_extpar('Cant deallocate the vector lu_tiles_lat_max')
       DEALLOCATE (lu_tiles_ncolumns, STAT = errorcode)
       IF (errorcode.NE.0) CALL abort_extpar('Cant deallocate the vector lu_tiles_ncolumns')
       DEALLOCATE (lu_tiles_nrows, STAT = errorcode)
       IF (errorcode.NE.0) CALL abort_extpar('Cant deallocate the vector lu_tiles_nrows')
       DEALLOCATE (lat_globcover, STAT = errorcode)
       IF (errorcode.NE.0) CALL abort_extpar('Cant deallocate the vector lat_globcover')
       DEALLOCATE (lon_globcover, STAT = errorcode)
       IF (errorcode.NE.0) CALL abort_extpar('Cant deallocate the vector lon_globcover')
       DEALLOCATE (fr_land_lu, STAT = errorcode)
       IF (errorcode.NE.0) CALL abort_extpar('Cant deallocate the vector fr_land_lu')
       DEALLOCATE (ice_lu, STAT = errorcode)
       IF (errorcode.NE.0) CALL abort_extpar('Cant deallocate the vector ice_lu')
       DEALLOCATE (z0_lu, STAT = errorcode)
       IF (errorcode.NE.0) CALL abort_extpar('Cant deallocate the vector z0_lu')
       DEALLOCATE (z0_tot, STAT = errorcode)
       IF (errorcode.NE.0) CALL abort_extpar('Cant deallocate the vector z0_tot')
       DEALLOCATE (root_lu, STAT = errorcode)
       IF (errorcode.NE.0) CALL abort_extpar('Cant deallocate the vector root_lu')
       DEALLOCATE (plcov_mn_lu, STAT = errorcode)
       IF (errorcode.NE.0) CALL abort_extpar('Cant deallocate the vector plcov_mn_lu')
       DEALLOCATE (plcov_mx_lu, STAT = errorcode)
       IF (errorcode.NE.0) CALL abort_extpar('Cant deallocate the vector plcov_mx_lu')
       DEALLOCATE (lai_mn_lu, STAT = errorcode)
       IF (errorcode.NE.0) CALL abort_extpar('Cant deallocate the vector lai_mn_lu')
       DEALLOCATE (lai_mx_lu, STAT = errorcode)
       IF (errorcode.NE.0) CALL abort_extpar('Cant deallocate the vector lai_mx_lu')
       DEALLOCATE (rs_min_lu, STAT = errorcode)
       IF (errorcode.NE.0) CALL abort_extpar('Cant deallocate the vector rs_min_lu')
       DEALLOCATE (urban_lu, STAT = errorcode)
       IF (errorcode.NE.0) CALL abort_extpar('Cant deallocate the vector urban_lu')
       DEALLOCATE (for_d_lu, STAT = errorcode)
       IF (errorcode.NE.0) CALL abort_extpar('Cant deallocate the vector for_d_lu')
       DEALLOCATE (for_e_lu, STAT = errorcode)
       IF (errorcode.NE.0) CALL abort_extpar('Cant deallocate the vector for_e_lu')
       DEALLOCATE (skinc_lu, STAT = errorcode)
       IF (errorcode.NE.0) CALL abort_extpar('Cant deallocate the vector skinc_lu')
       DEALLOCATE (emissivity_lu, STAT = errorcode)
       IF (errorcode.NE.0) CALL abort_extpar('Cant deallocate the vector emissivity_lu')
       DEALLOCATE (fr_ocean_lu, STAT = errorcode)
       IF (errorcode.NE.0) CALL abort_extpar('Cant deallocate the vector fr_ocean_lu')
       DEALLOCATE (lu_class_fraction, STAT = errorcode)
       IF (errorcode.NE.0) CALL abort_extpar('Cant deallocate the vector lu_class_fraction')
       DEALLOCATE (lu_class_npixel, STAT = errorcode)
       IF (errorcode.NE.0) CALL abort_extpar('Cant deallocate the vector lu_class_npixel')
       DEALLOCATE (lu_tot_npixel, STAT = errorcode)
       IF (errorcode.NE.0) CALL abort_extpar('Cant deallocate the vector lu_tot_npixel')

     END SUBROUTINE deallocate_landuse_data


END MODULE mo_globcover_data

