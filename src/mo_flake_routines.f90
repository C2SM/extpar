!+ Fortran module with routines for flake data settings
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_0         2010/12/21 Hermann Asensio
!  Initial release
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!> Fortran module with routines for flake data settings
!> \author Hermann Asensio
!>
MODULE mo_flake_routines

!> kind parameters are defined in MODULE data_parameters
USE mo_kind, ONLY: wp, &
                   i8, &
                   i4

USE netcdf,      ONLY :   &
  nf90_open,              &
  nf90_close,             &
  nf90_inquire,           &
  nf90_inquire_dimension, &
  nf90_inquire_variable,  &
  nf90_inq_attname,       &
  nf90_inquire_attribute, &
  nf90_get_att,           &
  nf90_inquire_dimension, &
  nf90_inq_varid,          &
  nf90_get_var,            &
  nf90_noerr,              &
  nf90_strerror

USE netcdf,      ONLY:     &
  nf90_create,             &
  nf90_def_dim,            &
  nf90_def_var,            &
  nf90_enddef,             &
  nf90_redef,              &
  nf90_put_att,            &
  nf90_put_var

 
USE netcdf,      ONLY :   &
  NF90_CHAR,               &
  NF90_DOUBLE,             &
  NF90_FLOAT,              &
  NF90_INT,                &
  NF90_BYTE,               &
  NF90_SHORT


USE netcdf,      ONLY :   &
  NF90_GLOBAL,             &
  NF90_UNLIMITED,          &
  NF90_CLOBBER,            &
  NF90_NOWRITE


!> abort_extpar defined in MODULE utilities_extpar
USE mo_utilities_extpar, ONLY: abort_extpar
USE mo_io_utilities,           ONLY: check_netcdf
USE mo_io_units,          ONLY: filename_max

USE mo_grid_structures,        ONLY: reg_lonlat_grid
USE mo_base_geometry,          ONLY: geographical_coordinates

IMPLICIT NONE

PRIVATE

PUBLIC :: read_namelists_extpar_flake,  &
 &         get_dimension_flake_data,        &
 &         get_lonlat_flake_data


CONTAINS

!---------------------------------------------------------------------------
!> subroutine to read namelist for orography data settings for EXTPAR 
SUBROUTINE read_namelists_extpar_flake(namelist_file, &
                                         raw_data_flake_path, &
                                         raw_data_flake_filename, &
                                         flake_buffer_file, &
                                         flake_output_file)


  USE mo_utilities_extpar, ONLY: free_un ! function to get free unit number

  CHARACTER (len=filename_max), INTENT(IN) :: namelist_file !< filename with namelists for for EXTPAR settings
CHARACTER (len=filename_max), INTENT(OUT) :: raw_data_flake_path        !< path to raw data
CHARACTER (len=filename_max), INTENT(OUT) :: raw_data_flake_filename !< filename flake raw data

CHARACTER (len=filename_max), INTENT(OUT) :: flake_buffer_file !< name for flake buffer file
CHARACTER (len=filename_max), INTENT(OUT) :: flake_output_file !< name for flake output file


!> namelist with land use data input, flake
NAMELIST /flake_raw_data/ raw_data_flake_path, raw_data_flake_filename
!> namelist with filenames for land use data output, flake data
NAMELIST /flake_io_extpar/ flake_buffer_file, flake_output_file

   INTEGER           :: nuin !< unit number
   INTEGER (KIND=i4) :: ierr !< error flag


   nuin = free_un()  ! functioin free_un returns free Fortran unit number
   OPEN(nuin,FILE=TRIM(namelist_file), IOSTAT=ierr)

   READ(nuin, NML=flake_raw_data, IOSTAT=ierr)
   READ(nuin, NML=flake_io_extpar, IOSTAT=ierr)
   
   CLOSE(nuin)


END SUBROUTINE read_namelists_extpar_flake
!---------------------------------------------------------------------------

        !> inquire dimension information for flake raw data 
        SUBROUTINE get_dimension_flake_data(path_flake_file, &
                                          nlon_flake, &
                                          nlat_flake)


        CHARACTER (len=*), INTENT(in) :: path_flake_file         !< filename with path for flake raw data
        INTEGER (KIND=i8), INTENT(out) :: nlon_flake !< number of grid elements in zonal direction for flake data
        INTEGER (KIND=i8), INTENT(out) :: nlat_flake !< number of grid elements in meridional direction for flake data

        !local variables
        INTEGER :: ncid                             !< netcdf unit file number
        INTEGER :: ndimension                       !< number of dimensions in netcdf file
        INTEGER :: nVars                            !< number of variables in netcdf file
        INTEGER :: nGlobalAtts                      !< number of gloabal Attributes in netcdf file
        INTEGER :: unlimdimid                       !< id of unlimited dimension (e.g. time) in netcdf file

        INTEGER :: dimid                            !< id of dimension
        CHARACTER (len=80) :: dimname               !< name of dimensiona
        INTEGER :: length                           !< length of dimension


          ! open netcdf file 
        CALL check_netcdf( nf90_open(TRIM(path_flake_file),NF90_NOWRITE, ncid))

       ! look for numbers of dimensions, Variable, Attributes, and the dimid for the unlimited dimension (probably time)
       !; nf90_inquire input: ncid; nf90_inquire output: ndimension, nVars, nGlobalAtts,unlimdimid
       CALL check_netcdf (nf90_inquire(ncid,ndimension, nVars, nGlobalAtts,unlimdimid))
       !print *,'ncid,ndimension, nVars, nGlobalAtts,unlimdimid',ncid,ndimension, nVars, nGlobalAtts,unlimdimid


       !; the dimid in netcdf-files is counted from 1 to ndimension
       !; look for the name and length of the dimension with f90_inquire_dimension
       !; nf90_inquire_dimension input: ncid, dimid; nf90_inquire_dimension output: name, length
       DO dimid=1,ndimension
                    !print *,'dimension loop dimid ',dimid
         CALL check_netcdf( nf90_inquire_dimension(ncid,dimid, dimname, length) )
                     !print*, 'ncid,dimid, dimname, length',ncid,dimid, trim(dimname), length
         IF ( trim(dimname) == 'lon') nlon_flake=length          ! here I know that the name of zonal dimension is 'lon'
         IF ( trim(dimname) == 'lat') nlat_flake=length          ! here I know that the name of meridional dimension is 'lat'
       ENDDO


       ! close netcdf file 
       CALL check_netcdf( nf90_close( ncid))

       END SUBROUTINE get_dimension_flake_data

!----------------------------------------------------------------------------------------------------------------  
        !> get coordinates for flake raw data 
        SUBROUTINE get_lonlat_flake_data(path_flake_file, &
                                          nlon_flake, &
                                          nlat_flake, &
                                          lon_flake,  &
                                          lat_flake,  &
                                          flake_grid)

       USE mo_grid_structures, ONLY: reg_lonlat_grid

        CHARACTER (len=*), INTENT(in) :: path_flake_file         !< filename with path for flake raw data
        INTEGER (KIND=i8), INTENT(in) :: nlon_flake !< number of grid elements in zonal direction for flake data
        INTEGER (KIND=i8), INTENT(in) :: nlat_flake !< number of grid elements in meridional direction for flake data
        REAL (KIND=wp), INTENT(out)    :: lon_flake(1:nlon_flake) !< longitude of flake raw data
        REAL (KIND=wp), INTENT(out)    :: lat_flake(1:nlat_flake) !< latitude of flake raw data
        TYPE(reg_lonlat_grid), INTENT(OUT) :: flake_grid 
                                              !< structure with defenition of the raw data grid for the whole GLOBE dataset
        
        !local variables
        INTEGER :: ncid                             !< netcdf unit file number

        CHARACTER (LEN=80) :: varname  !< name of variable
        INTEGER :: varid               !< id of variable

        ! open netcdf file 
        CALL check_netcdf( nf90_open(TRIM(path_flake_file),NF90_NOWRITE, ncid))

        varname = 'lon' ! I know that the longitude coordinates for the GLC2000 data are stored in a variable called 'lon'

        CALL check_netcdf( nf90_inq_varid(ncid, TRIM(varname), varid))

        CALL check_netcdf(nf90_get_var(ncid, varid,  lon_flake))

        varname = 'lat' ! I know that the longitude coordinates for the GLC2000 data are stored in a variable called 'lon'

        CALL check_netcdf( nf90_inq_varid(ncid, TRIM(varname), varid))

        CALL check_netcdf(nf90_get_var(ncid, varid,  lat_flake))

        ! close netcdf file 
        CALL check_netcdf( nf90_close( ncid))

        ! define the values for the structure flake_grid
 
        flake_grid%start_lon_reg = lon_flake(1)
        flake_grid%end_lon_reg   = lon_flake(nlon_flake)
        flake_grid%start_lat_reg = lat_flake(1)
        flake_grid%end_lat_reg   = lat_flake(nlat_flake)
        flake_grid%dlon_reg      = (lon_flake(nlon_flake) - lon_flake(1)) / (nlon_flake - 1)
        flake_grid%dlat_reg      = (lat_flake(nlat_flake) - lat_flake(1)) / (nlat_flake - 1)
        flake_grid%nlon_reg      = nlon_flake
        flake_grid%nlat_reg      = nlat_flake

       END SUBROUTINE get_lonlat_flake_data

        !----------------------------------------------------------------------------------------------------------------

END MODULE mo_flake_routines
