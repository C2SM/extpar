!+ Fortran Module with data fields for the isa data
!
!
! Description:
! Fortran Module with data fields for the isa data
!
! Current Code Owner: DWD, <Name of person responsible for this code>
!    <smail, phone, fax and email>
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V2_4         2015-05-21 Hendrik Wouters
!  Initial release
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!=======================================================================
!> Fortran Module with data fields for the Globcover data
!> \author Hermann Asensio
!!
MODULE mo_isa_data

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

PUBLIC :: isa_grid,                &
          lon_isa,                 &
          lat_isa,                 &
          allocate_raw_isa_fields, &
! >mes
          ntiles_isa,              &   ! number of tiles in ISA
          len_isa_lon,                    &   ! Number of columns in land-use dataset
          len_isa_lat,                    &   ! Number of rows in land-use dataset
          max_tiles_isa,                  &   ! maximal possible number of tiles that can be read
          nc_tiles_isa,                   &   ! total number of columns in one ISA tile
          isa_tiles_lon_min,              &   ! starting longitude of every ISA tile
          isa_tiles_lon_max,              &   ! ending longitude of every ISA tile
          isa_tiles_lat_min,              &   ! starting latitude of every ISA tile
          isa_tiles_lat_max,              &   ! ending latitude of every ISA tile
          isa_tiles_ncolumns,             &   ! number of columns (lonitude increments) in each ISA tile
          isa_tiles_nrows,                &   ! number of rows (latitude increments) in each ISA tile
          isa_tiles_grid,          &
          fill_isa_data,           &   ! subroutine (intent(in) and intent(out))
          allocate_isa_data,       &
          deallocate_isa_data

PUBLIC :: isa_type !_br 14.04.16


PUBLIC :: undef_isa, minimal_isa

TYPE(reg_lonlat_grid) :: isa_grid !< structure with defenition of the raw data grid for the whole isa dataset

REAL (KIND=wp), ALLOCATABLE    :: lon_isa(:) !< longitude of isa raw data
REAL (KIND=wp), ALLOCATABLE    :: lat_isa(:) !< latitude of isa raw data

SAVE
TYPE(reg_lonlat_grid), ALLOCATABLE  :: isa_tiles_grid(:)

INTEGER(KIND=i4), ALLOCATABLE :: isa_tiles_ncolumns(:)
INTEGER(KIND=i4), ALLOCATABLE :: isa_tiles_nrows(:)
INTEGER(KIND=i4)    :: ntiles_isa    !  now in INPUT_LU  number of ISA tiles = 6
INTEGER(KIND=i4)    :: len_isa_lon,len_isa_lat
INTEGER, PARAMETER            :: max_tiles_isa = 1000
INTEGER(KIND=i4)              :: nc_tiles_isa

REAL(KIND=wp), ALLOCATABLE    :: isa_tiles_lon_min(:)
REAL(KIND=wp), ALLOCATABLE    :: isa_tiles_lon_max(:)
REAL(KIND=wp), ALLOCATABLE    :: isa_tiles_lat_min(:)
REAL(KIND=wp), ALLOCATABLE    :: isa_tiles_lat_max(:)
! <mes

REAL (KIND=wp) :: undef_isa = 0.0  !< undefined value for ISA data
REAL (KIND=wp) :: minimal_isa = 0.0 !< minimal ISA value 

INTEGER (KIND=i4)            :: isa_type = 1 !_br 14.04.16

CONTAINS

  !> allocate raw data fields
  SUBROUTINE allocate_raw_isa_fields(nrows,ncolumns)
  IMPLICIT NONE
  INTEGER (KIND=i8), INTENT(IN) :: nrows !< number of rows
  INTEGER (KIND=i8), INTENT(IN) :: ncolumns !< number of columns

  INTEGER :: errorcode !< error status variable

    ALLOCATE (lon_isa(1:ncolumns), STAT=errorcode)
        IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array lon_isa')
    lon_isa = 0.0

     ALLOCATE (lat_isa(1:nrows), STAT=errorcode)
        IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array lat_isa')
    lat_isa = 0.0

  END  SUBROUTINE allocate_raw_isa_fields

  !----------------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------
! >mes

     SUBROUTINE allocate_isa_data(ntiles)   
       IMPLICIT NONE
       INTEGER, INTENT (IN) :: ntiles       ! number of tiles: 6 for ISA
       INTEGER :: errorcode

       ALLOCATE (isa_tiles_lon_min(1:ntiles), STAT = errorcode)
       IF (errorcode.NE.0) CALL abort_extpar('Cant allocate the vector isa_tiles_lon_min')
       ALLOCATE (isa_tiles_lon_max(1:ntiles), STAT = errorcode)
       IF (errorcode.NE.0) CALL abort_extpar('Cant allocate the vector isa_tiles_lon_max')
       ALLOCATE (isa_tiles_lat_min(1:ntiles), STAT = errorcode)
       IF (errorcode.NE.0) CALL abort_extpar('Cant allocate the vector isa_tiles_lat_min')
       ALLOCATE (isa_tiles_lat_max(1:ntiles), STAT = errorcode)
       IF (errorcode.NE.0) CALL abort_extpar('Cant allocate the vector isa_tiles_lat_max')
       
       ALLOCATE (isa_tiles_ncolumns(1:ntiles), STAT = errorcode)
       IF (errorcode.NE.0) CALL abort_extpar('Cant allocate the vector isa_tiles_ncolumns')
       ALLOCATE (isa_tiles_nrows(1:ntiles), STAT = errorcode)
       IF (errorcode.NE.0) CALL abort_extpar('Cant allocate the vector isa_tiles_nrows')

       ALLOCATE (isa_tiles_grid(1:ntiles), STAT = errorcode)
       IF (errorcode.NE.0) CALL abort_extpar('Cant allocate the vector isa_tiles_grid')
   
   isa_tiles_lon_min   = 0.0
   isa_tiles_lon_max   = 0.0
   isa_tiles_lat_min   = 0.0
   isa_tiles_lat_max   = 0.0
   isa_tiles_ncolumns  = 0
   isa_tiles_nrows     = 0
   
   END SUBROUTINE allocate_isa_data

! <mes

  !----------------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------
! >mes

   SUBROUTINE fill_isa_data(raw_data_isa_path,       &
                                  raw_data_isa_filename,   &
                                  isa_tiles_lon_min,       &
                                  isa_tiles_lon_max,       &    ! the allocated vectors need to be filled with the respective valu &
!& e
                                  isa_tiles_lat_min,       &
                                  isa_tiles_lat_max,       &
                                  nc_tiles_isa)
   IMPLICIT NONE
   SAVE
   CHARACTER (len=filename_max),INTENT(IN) :: raw_data_isa_path
   CHARACTER (len=filename_max),INTENT(IN) :: raw_data_isa_filename(:)
   REAL(KIND=wp), INTENT(OUT)   :: isa_tiles_lon_min(1:ntiles_isa) 
   REAL(KIND=wp), INTENT(OUT)   :: isa_tiles_lon_max(1:ntiles_isa)
   REAL(KIND=wp), INTENT(OUT)   :: isa_tiles_lat_min(1:ntiles_isa)
   REAL(KIND=wp), INTENT(OUT)   :: isa_tiles_lat_max(1:ntiles_isa)
   INTEGER(KIND=i4), INTENT(OUT):: nc_tiles_isa
   INTEGER(KIND=i4)    :: i    ! i is a counter, errorcode is used to check if allocation was succ &
!& essful
   INTEGER(KIND=i4)    :: ncid
   INTEGER(KIND=i4)    :: dimID_lat, dimID_lon, varID_lat, varID_lon                  
   REAL(KIND=wp)       :: half_gridp                          ! distance of half a grid point as the grid point is centered on a  &
!& ISA pixel
    
   half_gridp = 0.001388888889
!   print*, half_gridp

     DO i = 1,ntiles_isa
!_br 12.04.16 included the path raw_data_isa_path
       CALL check_netcdf(nf90_open(path =TRIM(raw_data_isa_path)//TRIM(raw_data_isa_filename(i)), &
                                   mode = nf90_nowrite, ncid = ncid))    ! ISA file is opened 
!_br 12.04.16 end
       CALL check_netcdf(nf90_inq_dimid(ncid,"lon", dimID_lon))
       CALL check_netcdf(nf90_inq_dimid(ncid,"lat", dimID_lat))
       CALL check_netcdf(nf90_inquire_dimension(ncid,dimID_lon, len = isa_tiles_ncolumns(i)))          
       CALL check_netcdf(nf90_inquire_dimension(ncid,dimID_lat, len = isa_tiles_nrows(i))) 
       CALL check_netcdf(nf90_inq_varid(ncid, "lon", varID_lon))
       CALL check_netcdf(nf90_inq_varid(ncid, "lat", varID_lat))
       CALL check_netcdf(nf90_get_var(ncid, varID_lon, isa_tiles_lon_min(i), start = (/1/)))            ! reads in the first longi &
!& tude value of tile i
       CALL check_netcdf(nf90_get_var(ncid, varID_lon, isa_tiles_lon_max(i), start = (/isa_tiles_ncolumns(i)/))) ! reads in the las &
!& t longitude value of tile i
       CALL check_netcdf(nf90_get_var(ncid, varID_lat, isa_tiles_lat_max(i), start = (/1/)))            ! reads in the first latit &
!& ude value of tile i
       CALL check_netcdf(nf90_get_var(ncid, varID_lat, isa_tiles_lat_min(i), start = (/isa_tiles_nrows(i)/))) ! reads in the last l &
!& atitude value of tile i
       CALL check_netcdf(nf90_close(ncid))                                                             ! the netcdf file is close &
!& d again
       isa_tiles_lon_min(i) = isa_tiles_lon_min(i) - half_gridp !< half of a grid point must be
       isa_tiles_lon_max(i) = isa_tiles_lon_max(i) + half_gridp !< added, as the ISA data
       isa_tiles_lat_min(i) = isa_tiles_lat_min(i) - half_gridp !< is located at the pixel center
       isa_tiles_lat_max(i) = isa_tiles_lat_max(i) + half_gridp
       
       len_isa_lon=isa_tiles_ncolumns(i)
       len_isa_lat=isa_tiles_nrows(i)

   print*, 'ISA TILE ',ntiles_isa,': NLON,NLAT ',len_isa_lon,len_isa_lat

     ENDDO

     nc_tiles_isa = isa_tiles_ncolumns(1)


   END SUBROUTINE fill_isa_data

! <mes


   SUBROUTINE  deallocate_isa_data()

     USE mo_isa_tg_fields, ONLY: &
                       &        isa_field,         &
                       &        isa_tot_npixel
     
     IMPLICIT NONE     

     INTEGER :: errorcode

       DEALLOCATE (isa_tiles_lon_min, STAT = errorcode)
       IF (errorcode.NE.0) CALL abort_extpar('Cant deallocate the vector isa_tiles_lon_min')
       DEALLOCATE (isa_tiles_lon_max, STAT = errorcode)
       IF (errorcode.NE.0) CALL abort_extpar('Cant deallocate the vector isa_tiles_lon_max')
       DEALLOCATE (isa_tiles_lat_min, STAT = errorcode)
       IF (errorcode.NE.0) CALL abort_extpar('Cant deallocate the vector isa_tiles_lat_min')
       DEALLOCATE (isa_tiles_lat_max, STAT = errorcode)
       IF (errorcode.NE.0) CALL abort_extpar('Cant deallocate the vector isa_tiles_lat_max')
       DEALLOCATE (isa_tiles_ncolumns, STAT = errorcode)
       IF (errorcode.NE.0) CALL abort_extpar('Cant deallocate the vector isa_tiles_ncolumns')
       DEALLOCATE (isa_tiles_nrows, STAT = errorcode)
       IF (errorcode.NE.0) CALL abort_extpar('Cant deallocate the vector isa_tiles_nrows')
       DEALLOCATE (lat_isa, STAT = errorcode)
       IF (errorcode.NE.0) CALL abort_extpar('Cant deallocate the vector lat_isa')
       DEALLOCATE (lon_isa, STAT = errorcode)
       IF (errorcode.NE.0) CALL abort_extpar('Cant deallocate the vector lon_isa')
       DEALLOCATE (isa_field, STAT = errorcode)
       IF (errorcode.NE.0) CALL abort_extpar('Cant deallocate the vector isa_field')
       DEALLOCATE (isa_tot_npixel, STAT = errorcode)
       IF (errorcode.NE.0) CALL abort_extpar('Cant deallocate the vector isa_tot_npixel')

     END SUBROUTINE deallocate_isa_data


END MODULE mo_isa_data

