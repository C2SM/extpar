!+ Fortran module with routines for landuse data
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_0         2010/12/21 Hermann Asensio
!  Initial release
! V1_3         2011/04/19 Hermann Asensio
! introduce Globcover 2009 land use data set for external parameters
! V2_0         2013/06/04 Martina Messmer
!  adjusted to the Globcover tiles
!  introduction of a deallocation subroutine
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!> Fortran module with routines for landuse data 
!> \author Hermann Asensio
!>
MODULE mo_landuse_routines

!> kind parameters are defined in MODULE data_parameters
USE mo_kind, ONLY: wp, &
                   i8, &
                   i4, &
                   i2

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
  nf90_inq_varid,         &
  nf90_get_var,           &
  nf90_noerr,             &
  nf90_strerror

USE netcdf,      ONLY:     &
  nf90_create,             &
  nf90_def_dim,            &
  nf90_def_var,            &
  nf90_enddef,             &
  nf90_redef,              &
  nf90_put_att,            &
  nf90_put_var

 
USE netcdf,      ONLY :    &
  NF90_CHAR,               &
  NF90_DOUBLE,             &
  NF90_FLOAT,              &
  NF90_INT,                &
  NF90_BYTE,               &
  NF90_SHORT


USE netcdf,      ONLY :    &
  NF90_GLOBAL,             &
  NF90_UNLIMITED,          &
  NF90_CLOBBER,            &
  NF90_NOWRITE


!> abort_extpar defined in MODULE utilities_extpar
USE mo_utilities_extpar, ONLY: abort_extpar
USE mo_io_utilities,     ONLY: check_netcdf
USE mo_io_units,         ONLY: filename_max


USE mo_grid_structures,  ONLY: reg_lonlat_grid
USE mo_base_geometry,    ONLY: geographical_coordinates

IMPLICIT NONE

PRIVATE

PUBLIC :: read_namelists_extpar_land_use,  &
 &         get_dimension_glc2000_data,     &
 &         get_lonlat_glc2000_data

PUBLIC :: get_dimension_glcc_data, &
  &       get_lonlat_glcc_data

PUBLIC :: get_dimension_globcover_data, &
  &       get_lonlat_globcover_data,    &
  &       get_globcover_tiles_grid,     &
  &       det_band_globcover_data,      &
  &       get_globcover_data_block,     &
  &       get_globcover_tile_block_indices

PUBLIC :: get_dimension_ecoclimap_data, &
  &       get_lonlat_ecoclimap_data

CONTAINS

!---------------------------------------------------------------------------
!> subroutine to read namelist for orography data settings for EXTPAR 
SUBROUTINE read_namelists_extpar_land_use(namelist_file,            &
                                         i_landuse_data,            &
                                         raw_data_lu_path,          &
                                         raw_data_lu_filename,      &
                                         ntiles_globcover,          &
                                         ilookup_table_lu,          &
                                         lu_buffer_file,            &
                                         lu_output_file,            &
                                         raw_data_glcc_path_opt,    &
                                         raw_data_glcc_filename_opt,&
                                         ilookup_table_glcc_opt,    &
                                         glcc_buffer_file_opt,      &
                                         glcc_output_file_opt)


  USE mo_utilities_extpar, ONLY: free_un ! function to get free unit number
  
! >mes
  USE mo_globcover_data,   ONLY: max_tiles_lu, ncolumn_tiles, nrow_tiles
!  USE mo_globcover_data,   ONLY: ntiles_globcover
!<mes
  
  CHARACTER (len=filename_max), INTENT(IN) :: namelist_file !< filename with namelists for for EXTPAR settings
! land use
  INTEGER, INTENT(OUT) :: i_landuse_data !< integer switch to choose a land use raw data set
                                         !! 1 Globcover2009, 2 GLC2000, 3 GLCC
  CHARACTER (len=filename_max), INTENT(OUT) :: raw_data_lu_path        !< path to raw data
  CHARACTER (len=filename_max), INTENT(OUT) :: raw_data_lu_filename(1:max_tiles_lu) !< filename lu raw data
  INTEGER, INTENT(OUT) :: ilookup_table_lu !< integer switch to choose a lookup table

  CHARACTER (len=filename_max), INTENT(OUT) :: lu_buffer_file !< name for landuse buffer file
  CHARACTER (len=filename_max), INTENT(OUT) :: lu_output_file !< name for landuse output file
  INTEGER, INTENT(OUT) :: ntiles_globcover
!--
  CHARACTER (len=filename_max), INTENT(OUT), OPTIONAL :: raw_data_glcc_path_opt        !< path to raw data
  CHARACTER (len=filename_max), INTENT(OUT), OPTIONAL :: raw_data_glcc_filename_opt !< filename glc2000 raw data
  INTEGER, INTENT(OUT), OPTIONAL :: ilookup_table_glcc_opt  !< integer switch to choose a lookup table

  CHARACTER (len=filename_max), INTENT(OUT), OPTIONAL :: glcc_buffer_file_opt    !< name for glcc buffer file
  CHARACTER (len=filename_max), INTENT(OUT), OPTIONAL :: glcc_output_file_opt    !< name for glcc output file
!--
  CHARACTER (len=filename_max) :: raw_data_glcc_path        !< path to raw data
  CHARACTER (len=filename_max) :: raw_data_glcc_filename !< filename glc2000 raw data
  INTEGER :: ilookup_table_glcc  !< integer switch to choose a lookup table


  CHARACTER (len=filename_max) :: glcc_buffer_file    !< name for glcc buffer file
  CHARACTER (len=filename_max) :: glcc_output_file    !< name for glcc output file

  !> namelist with land use data input
  NAMELIST /lu_raw_data/ raw_data_lu_path, raw_data_lu_filename, i_landuse_data, ilookup_table_lu, &
                         ntiles_globcover, ncolumn_tiles
  !> namelist with filenames for land use data output
  NAMELIST /lu_io_extpar/ lu_buffer_file, lu_output_file

  !> namelist with land use data input, glcc
  NAMELIST /glcc_raw_data/ raw_data_glcc_path, raw_data_glcc_filename, ilookup_table_glcc
  
  !> namelist with filenames for land use data output. glcc data
  NAMELIST /glcc_io_extpar/  glcc_buffer_file, glcc_output_file


  INTEGER           :: nuin !< unit number
  INTEGER (KIND=i4) :: ierr !< error flag
  ntiles_globcover = 6

  nuin = free_un()  ! functioin free_un returns free Fortran unit number
  OPEN(nuin,FILE=TRIM(namelist_file), IOSTAT=ierr)

  READ(nuin, NML=lu_raw_data, IOSTAT=ierr)

  READ(nuin, NML=lu_io_extpar, IOSTAT=ierr)
   
  ! If optional argument is present for output, copy the value from the local variable to the output argument variable
  IF (PRESENT(raw_data_glcc_path_opt)) THEN
    READ(nuin, NML=glcc_raw_data, IOSTAT=ierr)
    READ(nuin, NML=glcc_io_extpar, IOSTAT=ierr)
  ENDIF

  CLOSE(nuin)

    ! If optional argument is present for output, copy the value from the local variable to the output argument variable
  IF (PRESENT(raw_data_glcc_path_opt)) raw_data_glcc_path_opt = TRIM(raw_data_glcc_path)


  ! If optional argument is present for output, copy the value from the local variable to the output argument variable
  IF (PRESENT(raw_data_glcc_filename_opt)) raw_data_glcc_filename_opt = TRIM(raw_data_glcc_filename)

  ! If optional argument is present for output, copy the value from the local variable to the output argument variable
  IF (PRESENT(ilookup_table_glcc_opt)) ilookup_table_glcc = ilookup_table_glcc_opt 


    ! If optional argument is present for output, copy the value from the local variable to the output argument variable
  IF (PRESENT(glcc_buffer_file_opt)) glcc_buffer_file_opt = TRIM(glcc_buffer_file)

   ! If optional argument is present for output, copy the value from the local variable to the output argument variable
  IF (PRESENT(glcc_output_file_opt)) glcc_output_file_opt = TRIM(glcc_output_file)

  IF (ntiles_globcover == 1) THEN
    nrow_tiles = 1
    ncolumn_tiles = 1
  ELSE IF (ntiles_globcover /= 6) THEN
    nrow_tiles = ntiles_globcover / ncolumn_tiles
    IF (ncolumn_tiles * nrow_tiles /= ntiles_globcover) THEN
      CALL abort_extpar('ncolumn_tiles specified does not fit to ntiles_globcover')
    ENDIF
  ENDIF


END SUBROUTINE read_namelists_extpar_land_use
!---------------------------------------------------------------------------
        !> inquire dimension information for glc2000 raw data 
        SUBROUTINE get_dimension_glc2000_data(path_glc2000_file, &
                                          nlon_glc2000,          &
                                          nlat_glc2000)


        CHARACTER (len=*), INTENT(in) :: path_glc2000_file(:)         !< filename with path for glc2000 raw data
        INTEGER (KIND=i8), INTENT(out) :: nlon_glc2000 !< number of grid elements in zonal direction for glc2000 data
        INTEGER (KIND=i8), INTENT(out) :: nlat_glc2000 !< number of grid elements in meridional direction for glc2000 data

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
        CALL check_netcdf( nf90_open(TRIM(path_glc2000_file(1)),NF90_NOWRITE, ncid))


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
         IF ( trim(dimname) == 'lon') nlon_glc2000=length          ! here I know that the name of zonal dimension is 'lon'
         IF ( trim(dimname) == 'lat') nlat_glc2000=length          ! here I know that the name of meridional dimension is 'lat'
       ENDDO


       ! close netcdf file 
       CALL check_netcdf( nf90_close( ncid))

       END SUBROUTINE get_dimension_glc2000_data

!----------------------------------------------------------------------------------------------------------------  
        !> get coordinates for glc2000 raw data 
        SUBROUTINE get_lonlat_glc2000_data(path_glc2000_file, &
                                          nlon_glc2000,       &
                                          nlat_glc2000,       &
                                          lon_glc2000,        &
                                          lat_glc2000,        &
                                          glc2000_grid)

       USE mo_grid_structures, ONLY: reg_lonlat_grid

        CHARACTER (len=*), INTENT(in) :: path_glc2000_file(:)         !< filename with path for glc2000 raw data
        INTEGER (KIND=i8), INTENT(in) :: nlon_glc2000 !< number of grid elements in zonal direction for glc2000 data
        INTEGER (KIND=i8), INTENT(in) :: nlat_glc2000 !< number of grid elements in meridional direction for glc2000 data
        REAL (KIND=wp), INTENT(out)    :: lon_glc2000(1:nlon_glc2000) !< longitude of glc2000 raw data
        REAL (KIND=wp), INTENT(out)    :: lat_glc2000(1:nlat_glc2000) !< latitude of glc2000 raw data
        TYPE(reg_lonlat_grid), INTENT(OUT) :: glc2000_grid !< structure with defenition of the raw data grid
                                                           !  for the whole GLC2000 dataset
        
        !local variables
        INTEGER :: ncid                             !< netcdf unit file number

        CHARACTER (LEN=80) :: varname  !< name of variable
        INTEGER :: varid               !< id of variable

        ! open netcdf file 
        CALL check_netcdf( nf90_open(TRIM(path_glc2000_file(1)),NF90_NOWRITE, ncid))

        varname = 'lon' ! I know that the longitude coordinates for the GLC2000 data are stored in a variable called 'lon'

        CALL check_netcdf( nf90_inq_varid(ncid, TRIM(varname), varid))

        CALL check_netcdf(nf90_get_var(ncid, varid,  lon_glc2000))

        varname = 'lat' ! I know that the longitude coordinates for the GLC2000 data are stored in a variable called 'lon'

        CALL check_netcdf( nf90_inq_varid(ncid, TRIM(varname), varid))

        CALL check_netcdf(nf90_get_var(ncid, varid,  lat_glc2000))

        ! close netcdf file 
        CALL check_netcdf( nf90_close( ncid))

        ! define the values for the structure glc2000_grid
 
        glc2000_grid%start_lon_reg = lon_glc2000(1)
        glc2000_grid%end_lon_reg   = lon_glc2000(nlon_glc2000)
        glc2000_grid%start_lat_reg = lat_glc2000(1)
        glc2000_grid%end_lat_reg   = lat_glc2000(nlat_glc2000)
        glc2000_grid%dlon_reg      = (lon_glc2000(nlon_glc2000) - lon_glc2000(1)) / (nlon_glc2000 - 1)
        glc2000_grid%dlat_reg      = (lat_glc2000(nlat_glc2000) - lat_glc2000(1)) / (nlat_glc2000 - 1)
        glc2000_grid%nlon_reg      = nlon_glc2000
        glc2000_grid%nlat_reg      = nlat_glc2000

       END SUBROUTINE get_lonlat_glc2000_data

        !----------------------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------------------  
        !> get one row of glc2000 raw data 
        SUBROUTINE get_row_glc2000_data(path_glc2000_file, &
                                          nlon_glc2000,    &
                                          data_row,        &
                                          glc2000_data_row)

       USE mo_grid_structures, ONLY: reg_lonlat_grid

        CHARACTER (len=*), INTENT(in) :: path_glc2000_file         !< filename with path for glc2000 raw data
        INTEGER , INTENT(in) :: nlon_glc2000 !< number of grid elements in zonal direction for glc2000 data
        INTEGER , INTENT(in) :: data_row !< number or row for data to read in

        INTEGER , INTENT(OUT) :: glc2000_data_row(1:nlon_glc2000)


        !local variables
        INTEGER :: ncid                             !< netcdf unit file number

        CHARACTER (LEN=80) :: varname  !< name of variable
        INTEGER :: varid               !< id of variable

       ! open netcdf file 
        CALL check_netcdf( nf90_open(TRIM(path_glc2000_file),NF90_NOWRITE, ncid))

        varname = 'glc2000byte' ! I know that the longitude coordinates for the GLC2000 data are stored
                                ! in a variable called 'glc2000byte'

         CALL check_netcdf( nf90_inq_varid(ncid, TRIM(varname), varid))

         CALL check_netcdf(nf90_get_var(ncid, varid,  glc2000_data_row,  &
                       start=(/1,data_row/),count=(/nlon_glc2000,1/)))

       ! close netcdf file 
       CALL check_netcdf( nf90_close( ncid))

       END SUBROUTINE get_row_glc2000_data

      !----------------------------------------------------------------------------------------------------------------
      !----------------------------------------------------------------------------------------------------------------
      !> open netcdf-file and get netcdf unit file number
      SUBROUTINE open_netcdf_glc2000(path_glc2000_file, &
                                        ncid)
        CHARACTER (len=*), INTENT(in) :: path_glc2000_file         !< filename with path to GLC2000 tile
        INTEGER, INTENT(out) :: ncid                             !< netcdf unit file number

        ! open netcdf file 
         CALL check_netcdf( nf90_open(TRIM(path_glc2000_file),NF90_NOWRITE, ncid))

      END SUBROUTINE open_netcdf_glc2000

      !> close netcdf-file 
      SUBROUTINE close_netcdf_glc2000(ncid)
         INTEGER, INTENT(in) :: ncid                             !< netcdf unit file number

        ! close netcdf file 
        CALL check_netcdf( nf90_close( ncid))

       END SUBROUTINE close_netcdf_glc2000

       !----------------------------------------------------------------------------------------------------------------
       !----------------------------------------------------------------------------------------------------------------

       !----------------------------------------------------------------------------------------------------------------
       !----------------------------------------------------------------------------------------------------------------


       !> inquire dimension information for glcc raw data 
       SUBROUTINE get_dimension_glcc_data(path_glcc_file, &
                                          nlon_glcc,      &
                                          nlat_glcc)


        CHARACTER (len=*), INTENT(in) :: path_glcc_file(1)         !< filename with path for glcc raw data
        INTEGER (KIND=i8), INTENT(out) :: nlon_glcc !< number of grid elements in zonal direction for glcc data
        INTEGER (KIND=i8), INTENT(out) :: nlat_glcc !< number of grid elements in meridional direction for glcc data

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
        CALL check_netcdf( nf90_open(TRIM(path_glcc_file(1)),NF90_NOWRITE, ncid))


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
           IF ( TRIM(dimname) == 'lon') nlon_glcc=length          ! here I know that the name of zonal dimension is 'lon'
           IF ( TRIM(dimname) == 'lat') nlat_glcc=length          ! here I know that the name of meridional dimension is 'lat'
         ENDDO


         ! close netcdf file 
         CALL check_netcdf( nf90_close( ncid))

       END SUBROUTINE get_dimension_glcc_data

      !> get coordinates for glcc raw data 
      SUBROUTINE get_lonlat_glcc_data(path_glcc_file,&
                                          nlon_glcc, &
                                          nlat_glcc, &
                                          lon_glcc,  &
                                          lat_glcc,  &
                                          glcc_grid)

      USE mo_grid_structures, ONLY: reg_lonlat_grid

        CHARACTER (len=*), INTENT(in) :: path_glcc_file(1)         !< filename with path for glcc raw data
        INTEGER (KIND=i8), INTENT(in) :: nlon_glcc !< number of grid elements in zonal direction for glcc data
        INTEGER (KIND=i8), INTENT(in) :: nlat_glcc !< number of grid elements in meridional direction for glcc data
        REAL (KIND=wp), INTENT(out)    :: lon_glcc(1:nlon_glcc) !< longitude of glcc raw data
        REAL (KIND=wp), INTENT(out)    :: lat_glcc(1:nlat_glcc) !< latitude of glcc raw data
        TYPE(reg_lonlat_grid), INTENT(OUT) :: glcc_grid !< structure with defenition of the raw data grid
                                                        !  for the whole GLCC dataset

        !local variables
        INTEGER :: ncid                             !< netcdf unit file number

        CHARACTER (LEN=80) :: varname  !< name of variable
        INTEGER :: varid               !< id of variable

        ! open netcdf file 
        CALL check_netcdf( nf90_open(TRIM(path_glcc_file(1)),NF90_NOWRITE, ncid))

        varname = 'lon' ! I know that the longitude coordinates for the GLC2000 data are stored in a variable called 'lon'

        CALL check_netcdf( nf90_inq_varid(ncid, TRIM(varname), varid))

        CALL check_netcdf(nf90_get_var(ncid, varid,  lon_glcc))

        varname = 'lat' ! I know that the longitude coordinates for the GLC2000 data are stored in a variable called 'lon'

        CALL check_netcdf( nf90_inq_varid(ncid, TRIM(varname), varid))

        CALL check_netcdf(nf90_get_var(ncid, varid,  lat_glcc))

        ! close netcdf file 
        CALL check_netcdf( nf90_close( ncid))

        ! define the values for the structure glcc_grid

        glcc_grid%start_lon_reg = lon_glcc(1)
        glcc_grid%end_lon_reg   = lon_glcc(nlon_glcc)
        glcc_grid%start_lat_reg = lat_glcc(1)
        glcc_grid%end_lat_reg   = lat_glcc(nlat_glcc)
        glcc_grid%dlon_reg      = (lon_glcc(nlon_glcc) - lon_glcc(1)) / (nlon_glcc - 1)
        glcc_grid%dlat_reg      = (lat_glcc(nlat_glcc) - lat_glcc(1)) / (nlat_glcc - 1)
        glcc_grid%nlon_reg      = nlon_glcc
        glcc_grid%nlat_reg      = nlat_glcc

      END SUBROUTINE get_lonlat_glcc_data

      !----------------------------------------------------------------------------------------------------------------
      !----------------------------------------------------------------------------------------------------------------  
      !> get one row of glcc raw data 
      SUBROUTINE get_row_glcc_data(path_glcc_file,   &
                                          nlon_glcc, &
                                          data_row,  &
                                          glcc_data_row)

      USE mo_grid_structures, ONLY: reg_lonlat_grid

        CHARACTER (len=*), INTENT(in) :: path_glcc_file         !< filename with path for glcc raw data
        INTEGER , INTENT(in) :: nlon_glcc !< number of grid elements in zonal direction for glcc data
        INTEGER , INTENT(in) :: data_row !< number or row for data to read in

        INTEGER , INTENT(OUT) :: glcc_data_row(1:nlon_glcc)

        !local variables
        INTEGER :: ncid                             !< netcdf unit file number

        CHARACTER (LEN=80) :: varname  !< name of variable
        INTEGER :: varid               !< id of variable

        ! open netcdf file 
        CALL check_netcdf( nf90_open(TRIM(path_glcc_file),NF90_NOWRITE, ncid))

        varname = 'glccbyte' ! I know that the longitude coordinates for the GLC2000 data
                             ! are stored in a variable called 'glccbyte'

        CALL check_netcdf( nf90_inq_varid(ncid, TRIM(varname), varid))

        CALL check_netcdf(nf90_get_var(ncid, varid,  glcc_data_row,  &
          &             start=(/1,data_row/),count=(/nlon_glcc,1/)))

        ! close netcdf file 
        CALL check_netcdf( nf90_close( ncid))

      END SUBROUTINE get_row_glcc_data

      !----------------------------------------------------------------------------------------------------------------
      !----------------------------------------------------------------------------------------------------------------
      !> open netcdf-file and get netcdf unit file number
      SUBROUTINE open_netcdf_glcc(path_glcc_file, &
                                        ncid)
        CHARACTER (len=*), INTENT(in) :: path_glcc_file         !< filename with path to GLCC tile
        INTEGER, INTENT(out) :: ncid                             !< netcdf unit file number

        ! open netcdf file 
        CALL check_netcdf( nf90_open(TRIM(path_glcc_file),NF90_NOWRITE, ncid))

       END SUBROUTINE open_netcdf_glcc

       !> close netcdf-file 
       SUBROUTINE close_netcdf_glcc(ncid)
         INTEGER, INTENT(in) :: ncid                             !< netcdf unit file number
         ! close netcdf file 
         CALL check_netcdf( nf90_close( ncid))
       END SUBROUTINE close_netcdf_glcc

               !> inquire dimension information for globcover raw data
        SUBROUTINE get_dimension_globcover_data(nlon_globcover, &
                                          nlat_globcover)

        USE mo_globcover_data,   ONLY: max_tiles_lu,     &
                                       ntiles_globcover, &
                                       ncolumn_tiles,    &
                                       nrow_tiles,       &
                                       len_lu_lon, len_lu_lat

        INTEGER (KIND=i8), INTENT(OUT) :: nlon_globcover !< number of grid elements in zonal direction for globcover data
        INTEGER (KIND=i8), INTENT(OUT) :: nlat_globcover !< number of grid elements in meridional direction for globcover data

        !local variables
        INTEGER, PARAMETER :: nx=43200
        INTEGER, PARAMETER :: ny=27900

        IF(ntiles_globcover == 1) THEN
          nlon_globcover = len_lu_lon
          nlat_globcover = len_lu_lat
        ELSE IF(ntiles_globcover == 6) THEN
          nlon_globcover = 3 * nx
          nlat_globcover = 2 * ny
        ELSE
          nlon_globcover = ncolumn_tiles * nx
          nlat_globcover = nrow_tiles * ny
        END IF

       END SUBROUTINE get_dimension_globcover_data

       !----------------------------------------------------------------------------------------------------------------
        SUBROUTINE get_dimension_ecoclimap_data(nlon_ecoclimap, &
                                          nlat_ecoclimap)

        INTEGER (KIND=i8), INTENT(OUT) :: nlon_ecoclimap !< number of grid elements in zonal direction for globcover data
        INTEGER (KIND=i8), INTENT(OUT) :: nlat_ecoclimap !< number of grid elements in meridional direction for globcover data

        !local variables
	INTEGER, PARAMETER :: nx=43200
        INTEGER, PARAMETER :: ny=21600

        nlon_ecoclimap = nx
        nlat_ecoclimap = ny

       END SUBROUTINE get_dimension_ecoclimap_data

       !----------------------------------------------------------------------------------------------------------------

        !> get coordinates for globcover raw data
        SUBROUTINE get_lonlat_globcover_data(nlon_globcover, &
           &                                 nlat_globcover, &
           &                                 lon_globcover,  &
           &                                 lat_globcover,  &
           &                                 globcover_grid)

       USE mo_grid_structures, ONLY: reg_lonlat_grid
       USE mo_globcover_data,  ONLY: lu_tiles_lat_min, &
                                     lu_tiles_lat_max, &
                                     lu_tiles_lon_min, &
                                     lu_tiles_lon_max

        INTEGER (KIND=i8), INTENT(IN) :: nlon_globcover !< number of grid elements in zonal direction for globcover data
        INTEGER (KIND=i8), INTENT(IN) :: nlat_globcover !< number of grid elements in meridional direction for globcover data
        REAL (KIND=wp), INTENT(OUT)    :: lon_globcover(1:nlon_globcover) !< longitude of globcover raw data
        REAL (KIND=wp), INTENT(OUT)    :: lat_globcover(1:nlat_globcover) !< latitude of globcover raw data
        TYPE(reg_lonlat_grid), INTENT(OUT) :: globcover_grid !< structure with defenition of the raw data grid for
                                                             !  the whole GLOBECOVER dataset

        !local variables
!        REAL, PARAMETER ::  xmin_glc = -180.001388888889 ! area of glcover data: western longitude
!        REAL, PARAMETER ::  xmax_glc  =  179.998611111111 ! area of glcover data: eastern longitude
!        REAL, PARAMETER ::  ymax_glc  =   90.001388888888! area of glcover data: northern latitude
!        REAL, PARAMETER ::  ymin_glc = -64.9986111111111! area of glcover data: southern latitude

!        REAL, PARAMETER :: dx_glc =    0.0027777777  ! grid element size of glcover data pixel in zonal direction
!        REAL, PARAMETER :: dy_glc =   -0.0027777777  ! grid element size of glcover data pixel in meridional directionon
! >mes
        REAL (KIND=wp) ::  xmin_glc ! area of glcover data: western longitude
        REAL (KIND=wp) ::  xmax_glc ! area of glcover data: eastern longitude
        REAL (KIND=wp) ::  ymax_glc ! area of glcover data: northern latitude
        REAL (KIND=wp) ::  ymin_glc ! area of glcover data: southern latitude

        REAL (KIND=wp):: dx_glc  ! grid element size of glcover data pixel in zonal direction
        REAL (KIND=wp):: dy_glc  ! grid element size of glcover data pixel in meridional directionon

        INTEGER (KIND=i8) :: jx,jy


        xmin_glc = MINVAL(lu_tiles_lon_min)
        xmax_glc = MAXVAL(lu_tiles_lon_max)
        ymax_glc = MAXVAL(lu_tiles_lat_max)
        ymin_glc = MINVAL(lu_tiles_lat_min)
        dx_glc   = (xmax_glc - xmin_glc)/ REAL(nlon_globcover) !_br 04.04.14
        dy_glc   = -1.0 * (ymax_glc - ymin_glc) / REAL(nlat_globcover) !_br 04.04.14

! <mes

           DO jx=1,nlon_globcover
              lon_globcover(jx)  = xmin_glc + 0.5*dx_glc + (jx-1)*dx_glc
           ENDDO
           DO jy=1,nlat_globcover
            lat_globcover(jy) = ymax_glc + 0.5*dy_glc + (jy-1)*dy_glc !note negative increment!
           ENDDO

        ! define the values for the structure globcover_grid

       globcover_grid%start_lon_reg = lon_globcover(1)
       globcover_grid%end_lon_reg   = lon_globcover(nlon_globcover)
       globcover_grid%start_lat_reg = lat_globcover(1)
       globcover_grid%end_lat_reg   = lat_globcover(nlat_globcover)
!        globcover_grid%start_lon_reg = xmin_glc
!        globcover_grid%end_lon_reg   = xmax_glc
!        globcover_grid%start_lat_reg = ymax_glc
!        globcover_grid%end_lat_reg   = ymin_glc
        globcover_grid%dlon_reg      = dx_glc ! (lon_globcover(nlon_globcover) - lon_globcover(1)) / (nlon_globcover - 1)
        globcover_grid%dlat_reg      = dy_glc ! (lat_globcover(nlat_globcover) - lat_globcover(1)) / (nlat_globcover - 1)
        globcover_grid%nlon_reg      = nlon_globcover
        globcover_grid%nlat_reg      = nlat_globcover

       END SUBROUTINE get_lonlat_globcover_data

        !----------------------------------------------------------------------------------------------------------------
        !> get coordinates for ecoclimap raw data
        SUBROUTINE get_lonlat_ecoclimap_data(nlon_ecoclimap, &
           &                               nlat_ecoclimap, &
           &                               lon_ecoclimap,  &
           &                               lat_ecoclimap,  &
           &                               ecoclimap_grid)

       USE mo_grid_structures, ONLY: reg_lonlat_grid

        INTEGER (KIND=i8), INTENT(IN) :: nlon_ecoclimap !< number of grid elements in zonal direction for ecoclimap data
        INTEGER (KIND=i8), INTENT(IN) :: nlat_ecoclimap !< number of grid elements in meridional direction for ecoclimap data
        REAL (KIND=wp), INTENT(OUT)    :: lon_ecoclimap(1:nlon_ecoclimap) !< longitude of ecoclimap raw data
        REAL (KIND=wp), INTENT(OUT)    :: lat_ecoclimap(1:nlat_ecoclimap) !< latitude of ecoclimap raw data
        TYPE(reg_lonlat_grid), INTENT(OUT) :: ecoclimap_grid !< structure with defenition of the raw data grid
                                                             !  for the whole ecoclimap dataset

        !local variables
        REAL, PARAMETER ::  xmin_glc = -179.995833335    ! area of glcover data: western longitude
        REAL, PARAMETER ::  xmax_glc  = 179.995689334969 ! area of glcover data: eastern longitude
        REAL, PARAMETER ::  ymax_glc  = 89.995761335     ! area of glcover data: northern latitude
        REAL, PARAMETER ::  ymin_glc = -89.995833334    ! area of glcover data: southern latitude

        REAL, PARAMETER :: dx_glc =    0.0083333333  ! grid element size of glcover data pixel in zonal direction
        REAL, PARAMETER :: dy_glc =   -0.0083333333  ! grid element size of glcover data pixel in meridional directionon

        INTEGER (KIND=i8) :: jx,jy
           DO jx=1,nlon_ecoclimap
              lon_ecoclimap(jx)  = xmin_glc + (jx-1)*dx_glc
           ENDDO
           DO jy=1,nlat_ecoclimap
            lat_ecoclimap(jy) = ymax_glc + (jy-1)*dy_glc ! note negative increment!
           ENDDO

        ! define the values for the structure ecoclimap_grid

        ecoclimap_grid%start_lon_reg = lon_ecoclimap(1)
        ecoclimap_grid%end_lon_reg   = lon_ecoclimap(nlon_ecoclimap)
        ecoclimap_grid%start_lat_reg = lat_ecoclimap(1)
        ecoclimap_grid%end_lat_reg   = lat_ecoclimap(nlat_ecoclimap)
        ecoclimap_grid%dlon_reg      = dx_glc ! (lon_ecoclimap(nlon_ecoclimap) - lon_ecoclimap(1)) / (nlon_ecoclimap - 1)
        ecoclimap_grid%dlat_reg      = dy_glc ! (lat_ecoclimap(nlat_ecoclimap) - lat_ecoclimap(1)) / (nlat_ecoclimap - 1)
        ecoclimap_grid%nlon_reg      = nlon_ecoclimap
        ecoclimap_grid%nlat_reg      = nlat_ecoclimap

       END SUBROUTINE get_lonlat_ecoclimap_data

        !----------------------------------------------------------------------------------------------------------------
       !----------------------------------------------------------------------------------------------------------------
! >mes
        SUBROUTINE get_globcover_tiles_grid(globcover_tiles_grid)
          USE mo_globcover_data,  ONLY: ntiles_globcover,   &
                                        lu_tiles_lon_min,   &
                                        lu_tiles_lon_max,   &
                                        lu_tiles_lat_min,   &
                                        lu_tiles_lat_max,   &
                                        lu_tiles_ncolumns,  &
                                        lu_tiles_nrows

          TYPE(reg_lonlat_grid), INTENT(OUT):: globcover_tiles_grid(1:ntiles_globcover)

          INTEGER::  k      !counter

          REAL(KIND=wp) :: dlon
          REAL(KIND=wp) :: dlat

          DO k = 1,ntiles_globcover

            dlon = (lu_tiles_lon_max(k) - lu_tiles_lon_min(k)) / REAL(lu_tiles_ncolumns(k)) !_br 04.04.14
            dlat = -1. * (lu_tiles_lat_max(k) - lu_tiles_lat_min(k)) / REAL(lu_tiles_nrows(k))   
            !< latitude from north to south, negative increment !_br 04.04.14

            globcover_tiles_grid(k)%start_lon_reg = lu_tiles_lon_min(k) + 0.5*dlon
            globcover_tiles_grid(k)%end_lon_reg = lu_tiles_lon_max(k) - 0.5*dlon 

            globcover_tiles_grid(k)%start_lat_reg = lu_tiles_lat_max(k) + 0.5*dlat
            globcover_tiles_grid(k)%end_lat_reg = lu_tiles_lat_min(k) - 0.5*dlat 
            globcover_tiles_grid(k)%dlon_reg = dlon
            globcover_tiles_grid(k)%dlat_reg = dlat
            globcover_tiles_grid(k)%nlon_reg = lu_tiles_ncolumns(k)
            globcover_tiles_grid(k)%nlat_reg = lu_tiles_nrows(k)

          END DO

        END SUBROUTINE get_globcover_tiles_grid
! <mes

        !----------------------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------------------
! > mes
        !> determine grid description of band for GLOBCOVER

        SUBROUTINE det_band_globcover_data(globcover_grid,start_globcover_row,ta_grid)
          TYPE(reg_lonlat_grid),INTENT(IN) :: globcover_grid ! sturcture with the definition of the global data grid
                                                             ! of the GLOBCOVER data
          INTEGER, INTENT(IN) :: start_globcover_row         ! number of the start row of band 
          TYPE(reg_lonlat_grid), INTENT(OUT):: ta_grid       ! structure with definition of the target area grid.

          INTEGER(KIND=i4):: nrows = 2500              ! number of rows, set to 2500 as default
          ! band from east to west for the whole globe, like the complete globcover grid

          ta_grid%dlon_reg = globcover_grid%dlon_reg
          ta_grid%dlat_reg = globcover_grid%dlat_reg

          ta_grid%start_lon_reg = globcover_grid%start_lon_reg
          ta_grid%end_lon_reg = globcover_grid%end_lon_reg
          ta_grid%nlon_reg = globcover_grid%nlon_reg

          !latitude from north to south, negative increment
          ta_grid%nlat_reg = nrows
          ta_grid%start_lat_reg = globcover_grid%start_lat_reg + ta_grid%dlat_reg * (start_globcover_row-1)   
          !< latitude from north to south, note the negative increment!
          ta_grid%end_lat_reg = ta_grid%start_lat_reg + ta_grid%dlat_reg * (nrows - 1)   
          !< latitude from north to south, note the negative increment!
          
         ! check for the southern bound of the globcover data
         IF (ta_grid%end_lat_reg < globcover_grid%end_lat_reg) THEN ! band is at the southern bound
           ta_grid%end_lat_reg = globcover_grid%end_lat_reg
           ta_grid%nlat_reg    = NINT(((ta_grid%end_lat_reg - ta_grid%start_lat_reg) / ta_grid%dlat_reg)) + 1
         END IF

       END SUBROUTINE det_band_globcover_data
! < mes

        !----------------------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------------------
! >mes
! get GLOBCOVER data block for a given target area from the tile block indices
       SUBROUTINE get_globcover_data_block(ta_grid,              &
                                           globcover_tiles_grid, &
                                           ncids_globcover,      &
                                           lu_block)

         USE mo_grid_structures,  ONLY: reg_lonlat_grid  ! Definition of DATA Typeto describe a regular lonlat grid
         USE mo_globcover_data,   ONLY: ntiles_globcover, &
                                        nc_tiles_lu
         TYPE(reg_lonlat_grid), INTENT(IN)  :: ta_grid !< structure with definition of the target area grid
                                                       !  (dlon must be the same as for the whole GLOBCOVER dataset)
         TYPE(reg_lonlat_grid), INTENT(IN) :: globcover_tiles_grid(1:ntiles_globcover)
         !< structure with defenition of the raw data grid for the 16 GLOBECOVER tiles
         INTEGER , INTENT(IN) :: ncids_globcover(1:ntiles_globcover)  
         !< ncid for the GLOBCOVER tiles, the netcdf files have to be opened previously
         INTEGER (KIND=i2), INTENT(OUT) :: lu_block(1:ta_grid%nlon_reg,1:ta_grid%nlat_reg) !< a block of GLOBCOVER data 

       !local variables
       INTEGER (KIND=i4) :: globcover_startrow(1:ntiles_globcover) !< startrow indices for each GLOBCOVER tile
       INTEGER (KIND=i4) :: globcover_endrow(1:ntiles_globcover) !< endrow indices for each GLOBCOVER tile
       INTEGER (KIND=i4) :: globcover_startcolumn(1:ntiles_globcover) !< starcolumn indices for each GLOBCOVER tile
       INTEGER (KIND=i4) :: globcover_endcolumn(1:ntiles_globcover) !< endcolumn indices for each GLOBCOVER tile

       INTEGER (KIND=i4) :: ta_start_ie(1:ntiles_globcover)
       !< indices of target area block for first column of each GLOBCOVER tile
       INTEGER (KIND=i4) :: ta_end_ie(1:ntiles_globcover)
       !< indices of target area block for last column of each GLOBCOVER tile
       INTEGER (KIND=i4) :: ta_start_je(1:ntiles_globcover)
       !< indices of target area block for first row of each GLOBCOVER tile
       INTEGER (KIND=i4) :: ta_end_je(1:ntiles_globcover)
       !< indices of target area block for last row of each GLOBCOVER tile


       INTEGER (KIND=i2), ALLOCATABLE :: raw_lu_block(:,:) !< a block with GLOBCOVER data
       INTEGER :: varid               !< id of variable
       CHARACTER (LEN=80) :: varname  !< name of variable

       INTEGER :: nrows !< number of rows ! dimensions for raw_lu_block
       INTEGER :: ncolumns !< number of columns ! dimensions for raw_lu_block



       INTEGER :: i,j,k     ! counter
       INTEGER :: errorcode !< error status variable

       varname = 'GLOBCOVER'   ! I know that in the GLOBCOVER netcdf files the LU data is stored in a variable "GLOBCOVER"

       CALL get_globcover_tile_block_indices(ta_grid,              &
            &                                globcover_tiles_grid, &  
            &                                globcover_startrow,   &
            &                                globcover_endrow,     & 
            &                                globcover_startcolumn,&
            &                                globcover_endcolumn,  &
            &                                ta_start_ie,          & 
            &                                ta_end_ie,            &
            &                                ta_start_je,          &
            &                                ta_end_je)


       DO k = 1, ntiles_globcover
           IF ((globcover_startrow(k)/=0).AND.(globcover_startcolumn(k)/=0)) THEN
             nrows = globcover_endrow(k) - globcover_startrow(k) + 1
             ncolumns = globcover_endcolumn(k) - globcover_startcolumn(k) + 1
 
           ALLOCATE (raw_lu_block(1:ncolumns,1:nrows), STAT=errorcode)
             IF(errorcode/=0) CALL abort_extpar('Cant allocate the array raw_lu_block')

             CALL check_netcdf(nf90_inq_varid(ncids_globcover(k),TRIM(varname),varid)) ! get the varid of the altitude variable
             ! get the data into the raw_lu_block
             CALL check_netcdf(nf90_get_var(ncids_globcover(k), varid,  raw_lu_block,     & 
             &     start=(/globcover_startcolumn(k),globcover_startrow(k)/),count=(/ncolumns,nrows/)))
            

            lu_block(ta_start_ie(k):ta_end_ie(k),ta_start_je(k):ta_end_je(k)) = raw_lu_block(1:ncolumns,1:nrows)
            
             DEALLOCATE (raw_lu_block, STAT=errorcode)
              IF(errorcode/=0) CALL abort_extpar('Cant deallocate the array raw_lu_block')
          
           ENDIF
         ENDDO

       END SUBROUTINE get_globcover_data_block
! < mes      
       !----------------------------------------------------------------------------------------------------------------
       !----------------------------------------------------------------------------------------------------------------
! >mes
       ! get startrow, endrow, startcolumn and endcolumn of each GLOBCOVER tile (raw data) for a 
       ! given target area (ta_grid) and get start_indices (lon, lat) and end_indices of the target
       ! area for each GLOBCOVER tile
       ! The GLOBCOVER raw data are split in 6 tiles, so the target area may overlap several tiles.
       ! This subroutine determines the necessary indices to read in the GLOBCOVER data into the
       ! target area.
       ! GLOBCOVER tiles which are outside the target block will get indices with the value '0'

SUBROUTINE get_globcover_tile_block_indices(ta_grid,              &
         &                                  globcover_tiles_grid, &
         &                                  globcover_startrow,   &
         &                                  globcover_endrow,     & 
         &                                  globcover_startcolumn,&
         &                                  globcover_endcolumn,  &
         &                                  ta_start_ie,          &
         &                                  ta_end_ie,            &
         &                                  ta_start_je,          &
         &                                  ta_end_je)

USE mo_globcover_data, ONLY : ntiles_globcover,  &          !< GLOBCOVER raw data has 6 tiles
                              lu_tiles_lon_min,  &
                              lu_tiles_lon_max,  &
                              lu_tiles_lat_min,  &
                              lu_tiles_lat_max,  &
                              lu_tiles_ncolumns, &
                              lu_tiles_nrows

       USE mo_grid_structures, ONLY: reg_lonlat_grid  !< Definition of data type to describe a regular (lonlat) grid
       TYPE(reg_lonlat_grid), INTENT(IN) :: ta_grid !< structure with definition of the target area grid
                                                    !  (dlon must be the same as for the whole GLOBCOVER dataset)
       TYPE(reg_lonlat_grid), INTENT(IN) :: globcover_tiles_grid(1:ntiles_globcover) 
       !< structure with defenition of the raw data grid for the 6 GLOBCOVER tiles

       INTEGER (KIND=i4), INTENT(OUT) :: globcover_startrow(1:ntiles_globcover)
                                         !< startrow indices for each GLOBCOVER tile
       INTEGER (KIND=i4), INTENT(OUT) :: globcover_endrow(1:ntiles_globcover)
                                         !< endrow indices for each GLOBCOVER tile

       INTEGER (KIND=i4), INTENT(OUT) :: globcover_startcolumn(1:ntiles_globcover)
                                         !< starcolumn indices for each GLOBCOVER tile
       INTEGER (KIND=i4), INTENT(OUT) :: globcover_endcolumn(1:ntiles_globcover)
                                         !< endcolumn indices for each GLOBCOVER tile

       INTEGER (KIND=i4), INTENT(OUT) :: ta_start_ie(1:ntiles_globcover)    
       !< indices of target area block for first column of each GLOBCOVER tile
       INTEGER (KIND=i4), INTENT(OUT) :: ta_end_ie(1:ntiles_globcover)      
       !< indices of target area block for last column of each GLOBCOVER tile
       INTEGER (KIND=i4), INTENT(OUT) :: ta_start_je(1:ntiles_globcover)  
       !< indices of target area block for first row of each GLOBCOVER tile
       INTEGER (KIND=i4), INTENT(OUT) :: ta_end_je(1:ntiles_globcover)   
       !< indices of target area block for last row of each GLOBCOVER tile

       
       INTEGER (KIND=i4) :: index_k !< index of GLOBCOVER tile which contains point_geo

       ! local variables

       INTEGER  :: i          ! index for tiles (i,j,m,n,o)
       INTEGER  :: j 
       INTEGER  :: m
       INTEGER  :: n
       INTEGER  :: o
       INTEGER  :: t_i_start 
       INTEGER  :: t_i_end
       INTEGER  :: t_j_start
       INTEGER  :: t_j_end
       INTEGER  :: undefined

       REAL (KIND=wp) :: point_lon_coor

       INTEGER (KIND=i4) :: startrow ! startrow for tile
       INTEGER (KIND=i4) :: endrow 
       INTEGER (KIND=i4) :: startcolumn
       INTEGER (KIND=i4) :: endcolumn

       REAL (KIND=wp) :: dlon
       REAL (KIND=wp) :: dlat

       INTEGER :: k

       undefined = 0
       globcover_startrow     = undefined
       globcover_endrow       = undefined
       globcover_startcolumn  = undefined
       globcover_endcolumn    = undefined
       ta_start_ie = undefined 
       ta_end_ie   = undefined
       ta_start_je = undefined
       ta_end_je   = undefined

       k=1                      ! determin dlon and dlat (are the same for all tiles)
       dlon = ta_grid%dlon_reg
       dlat = ta_grid%dlat_reg

        DO k = 1,ntiles_globcover   !loop over the tiles which overlap the target area

          startcolumn = NINT((ta_grid%start_lon_reg - globcover_tiles_grid(k)%start_lon_reg)/dlon) + 1
          ! here I want nearest index (NINT)

          IF (startcolumn < 1) THEN 
            globcover_startcolumn(k) = 1
            ! get the start index of the subtile for the target area block
            ta_start_ie(k) = NINT ((globcover_tiles_grid(k)%start_lon_reg - ta_grid%start_lon_reg)/dlon) + 1
            ! index of target area block

          ELSE IF (startcolumn > lu_tiles_ncolumns(k)) THEN
            globcover_startcolumn(k) = 0
            ta_start_ie(k) = 0
          ELSE
            globcover_startcolumn(k) = startcolumn
            ta_start_ie(k) = 1
          ENDIF

         ! get endcolumn for tile k
         endcolumn = NINT((ta_grid%end_lon_reg - globcover_tiles_grid(k)%start_lon_reg)/dlon) +1

         IF (endcolumn > lu_tiles_ncolumns(k)) THEN 
           globcover_endcolumn(k) = lu_tiles_ncolumns(k)
           ! get the end index of the subtile for the target area block
           ta_end_ie(k) = NINT ((globcover_tiles_grid(k)%end_lon_reg - ta_grid%start_lon_reg)/dlon) + 1
           ! index of target area block
         ELSE IF (endcolumn < 1) THEN
           globcover_endcolumn(k) = 0
           ta_end_ie(k) = 0
         ELSE
           globcover_endcolumn(k) = endcolumn
           ta_end_ie(k) = ta_grid%nlon_reg
         ENDIF

         ! get startrow for tile k
         startrow = NINT((ta_grid%start_lat_reg - globcover_tiles_grid(k)%start_lat_reg)/dlat) + 1
        
         IF (startrow < 1) THEN 
           globcover_startrow(k) = 1
           ! get the start index of the subtile for the target area block
           ta_start_je(k) = NINT ((globcover_tiles_grid(k)%start_lat_reg  - ta_grid%start_lat_reg)/dlat) + 1
           ! index of target area block
         ELSE IF (startrow > lu_tiles_nrows(k)) THEN
           globcover_startrow(k) = 0
           ta_start_je(k) = 0
         ELSE
           globcover_startrow(k) = startrow
           ta_start_je(k) = 1
         ENDIF

         ! get endrow for tile k
         endrow   = NINT(( ta_grid%end_lat_reg - globcover_tiles_grid(k)%start_lat_reg )/dlat)  + 1
        
         IF (endrow > lu_tiles_nrows(k)) THEN 
           globcover_endrow(k) = lu_tiles_nrows(k)
           ! get the start index of the subtile for the target area block
           ta_end_je(k) = NINT ((globcover_tiles_grid(k)%end_lat_reg -  ta_grid%start_lat_reg )/dlat) + 1
           ! index of target area block

         ELSE IF (endrow < 1) THEN
           globcover_endrow(k) = 0
           ta_end_je(k) = 0
         ELSE
           globcover_endrow(k) = endrow
           ta_end_je(k) =  ta_grid%nlat_reg
         ENDIF
 
       ENDDO  ! loop over the tiles 

     END SUBROUTINE get_globcover_tile_block_indices

! <mes

       !----------------------------------------------------------------------------------------------------------------
       !----------------------------------------------------------------------------------------------------------------

        !> get one row of globcover raw data
        SUBROUTINE get_row_globcover_data(path_globcover_file, &
                                          nlon_globcover,      &
                                          data_row,            &
                                          globcover_data_row)

       USE mo_grid_structures, ONLY: reg_lonlat_grid

        CHARACTER (LEN=filename_max), INTENT(IN) :: path_globcover_file         !< filename with path for globcover raw data
        INTEGER , INTENT(IN) :: nlon_globcover !< number of grid elements in zonal direction for globcover data
        INTEGER , INTENT(IN) :: data_row !< number or row for data to read in

        INTEGER , INTENT(OUT):: globcover_data_row(1:nlon_globcover)

        !local variables
        INTEGER :: ncid                             !< netcdf unit file number

        CHARACTER (LEN=80) :: varname  !< name of variable
        INTEGER :: varid               !< id of variable

       ! open netcdf file
        CALL check_netcdf( nf90_open(TRIM(path_globcover_file),NF90_NOWRITE, ncid))

        varname = 'GLOBCOVER' ! I know that the globcover data are stored in a variable called 'GLOBCOVER'

         CALL check_netcdf( nf90_inq_varid(ncid, TRIM(varname), varid))

         CALL check_netcdf(nf90_get_var(ncid, varid,  globcover_data_row,  &
                       start=(/1,data_row/),count=(/nlon_globcover,1/)))

       ! close netcdf file
       CALL check_netcdf( nf90_close( ncid))

       END SUBROUTINE get_row_globcover_data

      !----------------------------------------------------------------------------------------------------------------
      !> get one row of ecoclimap raw data
      SUBROUTINE get_row_ecoclimap_data(path_ecoclimap_file, &
                                          nlon_ecoclimap, &
                                          data_row, &
                                          ecoclimap_data_row)

        USE mo_grid_structures, ONLY: reg_lonlat_grid

        CHARACTER (LEN=filename_max), INTENT(IN) :: path_ecoclimap_file   !< filename with path for ecoclimap raw data
        INTEGER , INTENT(IN) :: nlon_ecoclimap !< number of grid elements in zonal direction for ecoclimap data
        INTEGER , INTENT(IN) :: data_row !< number or row for data to read in

        INTEGER , INTENT(OUT) :: ecoclimap_data_row(1:nlon_ecoclimap)

        !local variables
        INTEGER :: ncid                             !< netcdf unit file number

        CHARACTER (LEN=80) :: varname  !< name of variable
        INTEGER :: varid               !< id of variable

        ! open netcdf file
        CALL check_netcdf( nf90_open(TRIM(path_ecoclimap_file),NF90_NOWRITE, ncid))

        varname = 'landuse' ! I know that the ecoclimap data are stored in a variable called 'LANDUSE'

        CALL check_netcdf( nf90_inq_varid(ncid, TRIM(varname), varid))

        CALL check_netcdf(nf90_get_var(ncid, varid,  ecoclimap_data_row,  &
                       start=(/1,data_row/),count=(/nlon_ecoclimap,1/)))

        ! close netcdf file

        CALL check_netcdf( nf90_close( ncid))

      END SUBROUTINE get_row_ecoclimap_data

      !----------------------------------------------------------------------------------------------------------------


END MODULE mo_landuse_routines

