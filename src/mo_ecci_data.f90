!+ Fortran Module with data fields for the ecci data
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
!> Fortran Module with data fields for the Ecci data
!> \author Hermann Asensio
!!
MODULE mo_ecci_data

  USE mo_logging
  USE netcdf
  USE mo_kind,                  ONLY: wp, i4

  USE mo_grid_structures,       ONLY: reg_lonlat_grid

  USE mo_io_utilities,          ONLY: check_netcdf, &
       &                              join_path

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

  PUBLIC :: ecci_grid,                &
       &    lon_ecci,                 &
       &    lat_ecci,                 &
       &    allocate_raw_ecci_fields, &
       &    ntiles_ecci,              &   ! number of tiles in ECCI
       &    ncolumn_tiles_ecci,                 &   ! number of columns in tile matrix
       &    nrow_tiles_ecci,                    &   ! number of rows in tile matrix
       &    len_lu_lon_ecci,                    &   ! Number of columns in land-use dataset
       &    len_lu_lat_ecci,                    &   ! Number of rows in land-use dataset
       &    max_tiles_lu_ecci,                  &   ! maximal possible number of tiles that can be read
       &    nc_tiles_lu_ecci,                   &   ! total number of columns in one ECCI tile
       &    lu_tiles_lon_min_ecci,              &   ! starting longitude of every ECCI tile
       &    lu_tiles_lon_max_ecci,              &   ! ending longitude of every ECCI tile
       &    lu_tiles_lat_min_ecci,              &   ! starting latitude of every ECCI tile
       &    lu_tiles_lat_max_ecci,              &   ! ending latitude of every ECCI tile
       &    lu_tiles_ncolumns_ecci,             &   ! number of columns (lonitude increments) in each ECCI tile
       &    lu_tiles_nrows_ecci,                &   ! number of rows (latitude increments) in each ECCI tile
       &    ecci_tiles_grid,          &
       &    fill_ecci_data,           &   ! subroutine (intent(in) and intent(out))
       &    allocate_ecci_data,       &
       &    deallocate_ecci_fields

  TYPE(reg_lonlat_grid)          :: ecci_grid !< structure with defenition of the raw data grid for the whole ecci dataset

  REAL (KIND=wp), ALLOCATABLE    :: lon_ecci(:), & !< longitude of ecci raw data
       &                            lat_ecci(:) !< latitude of ecci raw data

  SAVE
  TYPE(reg_lonlat_grid), ALLOCATABLE  :: ecci_tiles_grid(:)

  INTEGER(KIND=i4)              :: ntiles_ecci = 6, &  ! number of ECCI tiles
       &                           nc_tiles_lu_ecci, &
       &                           ncolumn_tiles_ecci = 3, &  ! number of used tiles E-W
       &                           nrow_tiles_ecci = 2, &     ! number of used tiles N-S
       &                           len_lu_lon_ecci, &
       &                           len_lu_lat_ecci

  INTEGER(KIND=i4), PARAMETER   :: max_tiles_lu_ecci = 1000

  INTEGER(KIND=i4), ALLOCATABLE :: lu_tiles_ncolumns_ecci(:), &
       &                           lu_tiles_nrows_ecci(:)

  REAL(KIND=wp), ALLOCATABLE    :: lu_tiles_lon_min_ecci(:), &
       &                           lu_tiles_lon_max_ecci(:), &
       &                           lu_tiles_lat_min_ecci(:), &
       &                           lu_tiles_lat_max_ecci(:)

  CONTAINS

  !> allocate raw data fields
  SUBROUTINE allocate_raw_ecci_fields(nrows,ncolumns)

    INTEGER (KIND=i4), INTENT(IN) :: nrows, & !< number of rows
         &                           ncolumns !< number of columns

    INTEGER(KIND=i4)              :: errorcode !< error status variable

    ALLOCATE (lon_ecci(1:ncolumns), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array lon_ecci',__FILE__,__LINE__)
    lon_ecci = 0.0

    ALLOCATE (lat_ecci(1:nrows), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array lat_ecci',__FILE__,__LINE__)
    lat_ecci = 0.0

  END  SUBROUTINE allocate_raw_ecci_fields

  SUBROUTINE allocate_ecci_data(ntiles_ecci)

    INTEGER(KIND=i4), INTENT (IN) :: ntiles_ecci       ! number of tiles: 6 for ECCI
    INTEGER(KIND=i4)              :: errorcode

    ALLOCATE (lu_tiles_lon_min_ecci(1:ntiles_ecci), STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant allocate the vector lu_tiles_lon_min_ecci',__FILE__,__LINE__)
    ALLOCATE (lu_tiles_lon_max_ecci(1:ntiles_ecci), STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant allocate the vector lu_tiles_lon_max_ecci',__FILE__,__LINE__)
    ALLOCATE (lu_tiles_lat_min_ecci(1:ntiles_ecci), STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant allocate the vector lu_tiles_lat_min_ecci',__FILE__,__LINE__)
    ALLOCATE (lu_tiles_lat_max_ecci(1:ntiles_ecci), STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant allocate the vector lu_tiles_lat_max_ecci',__FILE__,__LINE__)

    ALLOCATE (lu_tiles_ncolumns_ecci(1:ntiles_ecci), STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant allocate the vector lu_tiles_ncolumns_ecci',__FILE__,__LINE__)
    ALLOCATE (lu_tiles_nrows_ecci(1:ntiles_ecci), STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant allocate the vector lu_tiles_nrows_ecci',__FILE__,__LINE__)

    ALLOCATE (ecci_tiles_grid(1:ntiles_ecci), STAT = errorcode)
    IF (errorcode.NE.0) CALL logging%error('Cant allocate the vector ecci_tiles_grid',__FILE__,__LINE__)

    lu_tiles_lon_min_ecci   = 0.0
    lu_tiles_lon_max_ecci   = 0.0
    lu_tiles_lat_min_ecci   = 0.0
    lu_tiles_lat_max_ecci   = 0.0
    lu_tiles_ncolumns_ecci  = 0
    lu_tiles_nrows_ecci     = 0

  END SUBROUTINE allocate_ecci_data

  SUBROUTINE fill_ecci_data(raw_data_lu_path,       &
                             raw_data_lu_filename,   &
                             lu_tiles_lon_min_ecci,       &
                             lu_tiles_lon_max_ecci,       &  ! the allocated vectors need to be filled with the respective value
                             lu_tiles_lat_min_ecci,       &
                             lu_tiles_lat_max_ecci,       &
                             nc_tiles_lu)

    SAVE
    CHARACTER (len=*),INTENT(IN) :: raw_data_lu_path, &
         &                          raw_data_lu_filename(:)

    REAL(KIND=wp), INTENT(OUT)   :: lu_tiles_lon_min_ecci(1:ntiles_ecci), &
         &                          lu_tiles_lon_max_ecci(1:ntiles_ecci), &
         &                          lu_tiles_lat_min_ecci(1:ntiles_ecci), &
         &                          lu_tiles_lat_max_ecci(1:ntiles_ecci)

    INTEGER(KIND=i4), INTENT(OUT):: nc_tiles_lu

    ! local variables
    INTEGER(KIND=i4)             :: i, &
         &                          ncid, &
         &                          dimID_lat, dimID_lon, varID_lat, varID_lon

    REAL(KIND=wp)                :: half_gridp

    half_gridp = 0.001388888889

    DO i = 1,ntiles_ecci
      ! open ECCI file
      CALL check_netcdf(nf90_open(path =join_path(raw_data_lu_path,raw_data_lu_filename(i)), mode = nf90_nowrite, ncid = ncid))
      CALL check_netcdf(nf90_inq_dimid(ncid,"lon", dimID_lon))
      CALL check_netcdf(nf90_inq_dimid(ncid,"lat", dimID_lat))
      CALL check_netcdf(nf90_inquire_dimension(ncid,dimID_lon, len = lu_tiles_ncolumns_ecci(i)))
      CALL check_netcdf(nf90_inquire_dimension(ncid,dimID_lat, len = lu_tiles_nrows_ecci(i)))
      CALL check_netcdf(nf90_inq_varid(ncid, "lon", varID_lon))
      CALL check_netcdf(nf90_inq_varid(ncid, "lat", varID_lat))
      CALL check_netcdf(nf90_get_var(ncid, varID_lon, lu_tiles_lon_min_ecci(i), start = (/1/)))
      ! reads in the first longitude value of tile i
      CALL check_netcdf(nf90_get_var(ncid, varID_lon, lu_tiles_lon_max_ecci(i), start = (/lu_tiles_ncolumns_ecci(i)/)))
      ! reads in the last longitude value of tile i
      CALL check_netcdf(nf90_get_var(ncid, varID_lat, lu_tiles_lat_max_ecci(i), start = (/1/)))
      ! reads in the first latitude value of tile i
      CALL check_netcdf(nf90_get_var(ncid, varID_lat, lu_tiles_lat_min_ecci(i), start = (/lu_tiles_nrows_ecci(i)/)))
      ! reads in the last latitude value of tile i
      CALL check_netcdf(nf90_close(ncid))  ! the netcdf file is closed again
      lu_tiles_lon_min_ecci(i) = lu_tiles_lon_min_ecci(i) - half_gridp !< half of a grid point must be
      lu_tiles_lon_max_ecci(i) = lu_tiles_lon_max_ecci(i) + half_gridp !< added, as the ECCI data
      lu_tiles_lat_min_ecci(i) = lu_tiles_lat_min_ecci(i) - half_gridp !< is located at the pixel center
      lu_tiles_lat_max_ecci(i) = lu_tiles_lat_max_ecci(i) + half_gridp

      len_lu_lon_ecci=lu_tiles_ncolumns_ecci(i)
      len_lu_lat_ecci=lu_tiles_nrows_ecci(i)
    ENDDO

    nc_tiles_lu = lu_tiles_ncolumns_ecci(1)

  END SUBROUTINE fill_ecci_data

  SUBROUTINE  deallocate_ecci_fields()

     INTEGER(KIND=i4) :: errorcode

     DEALLOCATE (lu_tiles_lon_min_ecci, STAT = errorcode)
     IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector lu_tiles_lon_min_ecci',__FILE__,__LINE__)
     DEALLOCATE (lu_tiles_lon_max_ecci, STAT = errorcode)
     IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector lu_tiles_lon_max_ecci',__FILE__,__LINE__)
     DEALLOCATE (lu_tiles_lat_min_ecci, STAT = errorcode)
     IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector lu_tiles_lat_min_ecci',__FILE__,__LINE__)
     DEALLOCATE (lu_tiles_lat_max_ecci, STAT = errorcode)
     IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector lu_tiles_lat_max_ecci',__FILE__,__LINE__)

     DEALLOCATE (lu_tiles_ncolumns_ecci, STAT = errorcode)
     IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector lu_tiles_ncolumns_ecci',__FILE__,__LINE__)
     DEALLOCATE (lu_tiles_nrows_ecci, STAT = errorcode)
     IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector lu_tiles_nrows_ecci',__FILE__,__LINE__)

     DEALLOCATE (ecci_tiles_grid, STAT = errorcode)
     IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector ecci_tiles_grid',__FILE__,__LINE__)

     DEALLOCATE (lat_ecci, STAT = errorcode)
     IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector lat_ecci',__FILE__,__LINE__)
     DEALLOCATE (lon_ecci, STAT = errorcode)
     IF (errorcode.NE.0) CALL logging%error('Cant deallocate the vector lon_ecci',__FILE__,__LINE__)


  END SUBROUTINE deallocate_ecci_fields

END MODULE mo_ecci_data
