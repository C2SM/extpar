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

  USE mo_logging
  USE mo_kind,                  ONLY: wp,i4

  USE netcdf,                   ONLY:   &
       &                              nf90_open,              &
       &                              nf90_close,             &
       &                              nf90_inquire,           &
       &                              nf90_inquire_dimension, &
       &                              nf90_inquire_dimension, &
       &                              nf90_inq_varid,          &
       &                              nf90_get_var,            &
       &                              nf90_nowrite

  USE mo_io_utilities,          ONLY: check_netcdf
  USE mo_io_units,              ONLY: filename_max
  USE mo_utilities_extpar,      ONLY: free_un
  USE mo_grid_structures,       ONLY: reg_lonlat_grid

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: read_namelists_extpar_flake, &
     &       get_dimension_flake_data,   &
     &       get_lonlat_flake_data


  CONTAINS

  !---------------------------------------------------------------------------
  !> subroutine to read namelist for orography data settings for EXTPAR 
  SUBROUTINE read_namelists_extpar_flake(namelist_file, &
                                           raw_data_flake_path, &
                                           raw_data_flake_filename, &
                                           flake_buffer_file, &
                                           flake_output_file)

    CHARACTER (len=filename_max), INTENT(IN)  :: namelist_file !< filename with namelists for for EXTPAR settings
    CHARACTER (len=filename_max), INTENT(OUT) :: raw_data_flake_path, &        !< path to raw data
         &                                       raw_data_flake_filename, & !< filename flake raw data

         &                                       flake_buffer_file, & !< name for flake buffer file
         &                                       flake_output_file !< name for flake output file
    INTEGER (KIND=i4)                          :: ierr, nuin
    
    !> namelist with land use data input, flake
    NAMELIST /flake_raw_data/ raw_data_flake_path, raw_data_flake_filename
    !> namelist with filenames for land use data output, flake data
    NAMELIST /flake_io_extpar/ flake_buffer_file, flake_output_file

    nuin = free_un()

    OPEN(nuin,FILE=TRIM(namelist_file), IOSTAT=ierr)
    IF (ierr /= 0) THEN
      WRITE(message_text,*)'Cannot open ', TRIM(namelist_file)
      CALL logging%error(message_text,__FILE__, __LINE__) 
    ENDIF

    READ(nuin, NML=flake_raw_data, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      CALL logging%error('Cannot read in namelist flake_raw_data',__FILE__, __LINE__) 
    ENDIF

    READ(nuin, NML=flake_io_extpar, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      CALL logging%error('Cannot read in namelist flake_io_extpar',__FILE__, __LINE__) 
    ENDIF
    
    CLOSE(nuin)

  END SUBROUTINE read_namelists_extpar_flake
!---------------------------------------------------------------------------

        !> inquire dimension information for flake raw data 
        SUBROUTINE get_dimension_flake_data(path_flake_file, &
                                          nlon_flake, &
                                          nlat_flake)


        CHARACTER (len=*), INTENT(in) :: path_flake_file         !< filename with path for flake raw data
        INTEGER (KIND=i4), INTENT(out) :: nlon_flake !< number of grid elements in zonal direction for flake data
        INTEGER (KIND=i4), INTENT(out) :: nlat_flake !< number of grid elements in meridional direction for flake data

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


    CHARACTER (len=*), INTENT(in)      :: path_flake_file         !< filename with path for flake raw data
    INTEGER (KIND=i4), INTENT(in)      :: nlon_flake, & !< number of grid elements in zonal direction for flake data
         &                                nlat_flake !< number of grid elements in meridional direction for flake data
    REAL (KIND=wp), INTENT(out)        :: lon_flake(1:nlon_flake), & !< longitude of flake raw data
         &                                lat_flake(1:nlat_flake) !< latitude of flake raw data

    TYPE(reg_lonlat_grid), INTENT(OUT) :: flake_grid 
                                          !< structure with defenition of the raw data grid for the whole GLOBE dataset
    
    !local variables
    INTEGER(KIND=i4)                   :: ncid, varid

    CHARACTER (LEN=80)                 :: varname  !< name of variable

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

END MODULE mo_flake_routines
