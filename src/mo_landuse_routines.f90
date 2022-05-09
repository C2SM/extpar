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

  USE mo_logging
  USE mo_kind,                  ONLY: wp, i4, i2

  USE netcdf,                   ONLY: &
                                      nf90_close,   &
                                      nf90_get_var, &
                                      nf90_inquire, &
                                      nf90_inquire_dimension, &
                                      nf90_inq_varid, &
                                      nf90_nowrite,   &
                                      nf90_open

  USE mo_io_utilities,          ONLY: check_netcdf

  USE mo_io_units,              ONLY: filename_max


  USE mo_grid_structures,       ONLY: reg_lonlat_grid

  USE mo_utilities_extpar,      ONLY: free_un

  USE mo_globcover_data,        ONLY: max_tiles_lu, ncolumn_tiles, nrow_tiles, &
       &                              ntiles_globcover, &
       &                              len_lu_lon, len_lu_lat, &
       &                              lu_tiles_lat_min, &
       &                              lu_tiles_lat_max, &
       &                              lu_tiles_lon_min, &
       &                              lu_tiles_lon_max, &
       &                              lu_tiles_ncolumns,  &
       &                              lu_tiles_nrows

  USE mo_ecci_data,             ONLY: ncolumn_tiles_ecci, &
       &                              nrow_tiles_ecci, ntiles_ecci, &
       &                              len_lu_lon_ecci, len_lu_lat_ecci, &
       &                              lu_tiles_lat_min_ecci, &
       &                              lu_tiles_lat_max_ecci, &
       &                              lu_tiles_lon_min_ecci, &
       &                              lu_tiles_lon_max_ecci, &
       &                              lu_tiles_ncolumns_ecci,  &
                                      lu_tiles_nrows_ecci
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

  PUBLIC :: get_dimension_ecci_data, &
    &       get_lonlat_ecci_data,    &
    &       get_ecci_tiles_grid,     &
    &       det_band_ecci_data,      &
    &       get_ecci_data_block,     &
    &       get_ecci_tile_block_indices

  PUBLIC :: get_dimension_ecoclimap_data, &
    &       get_lonlat_ecoclimap_data

  CONTAINS

  !> subroutine to read namelist for orography data settings for EXTPAR
  SUBROUTINE read_namelists_extpar_land_use(namelist_file,            &
                                           i_landuse_data,            &
                                           l_use_corine,              &
                                           raw_data_lu_path,          &
                                           raw_data_lu_filename,      &
                                           ilookup_table_lu,          &
                                           lu_buffer_file,            &
                                           raw_data_glcc_path_opt,    &
                                           raw_data_glcc_filename_opt,&
                                           ilookup_table_glcc_opt,    &
                                           glcc_buffer_file_opt)


     CHARACTER (len=*), INTENT(IN)                      :: namelist_file !< filename with namelists for for EXTPAR settings

    LOGICAL, INTENT(OUT)                                :: l_use_corine  !< flag to use corine datasete instead of globcover


    INTEGER(KIND=i4), INTENT(OUT)                       :: i_landuse_data, &  !< integer switch to choose a land use raw data set
         &                                                 ilookup_table_lu !< integer switch to choose a lookup table

    INTEGER, INTENT(OUT), OPTIONAL                      :: ilookup_table_glcc_opt  !< integer switch to choose a lookup table

    CHARACTER (len=filename_max), INTENT(OUT)           :: raw_data_lu_path, &         !< path to raw data
         &                                                 raw_data_lu_filename(1:max_tiles_lu), &  !< filename lu raw data
         &                                                 lu_buffer_file  !< name for landuse buffer file

    CHARACTER (len=filename_max), INTENT(OUT), OPTIONAL :: raw_data_glcc_path_opt, &         !< path to raw data
         &                                                 raw_data_glcc_filename_opt, &  !< filename glc2000 raw data
         &                                                 glcc_buffer_file_opt     !< name for glcc buffer file

    CHARACTER (len=filename_max)                       :: raw_data_glcc_path, &         !< path to raw data
         &                                                raw_data_glcc_filename, &  !< filename glc2000 raw data
         &                                                glcc_buffer_file     !< name for glcc buffer file

    INTEGER(KIND=i4)                                   :: ilookup_table_glcc, &   !< integer switch to choose a lookup table
         &                                                nuin, &  !< unit number
         &                                                ierr !< error flag

    !> namelist with land use data input
    NAMELIST /lu_raw_data/ raw_data_lu_path, raw_data_lu_filename, i_landuse_data, ilookup_table_lu, &
                           ntiles_globcover, ncolumn_tiles, l_use_corine
    !> namelist with filenames for land use data output
    NAMELIST /lu_io_extpar/ lu_buffer_file

    !> namelist with land use data input, glcc
    NAMELIST /glcc_raw_data/ raw_data_glcc_path, raw_data_glcc_filename, ilookup_table_glcc

    !> namelist with filenames for land use data output. glcc data
    NAMELIST /glcc_io_extpar/  glcc_buffer_file

    nuin = free_un()  ! functioin free_un returns free Fortran unit number
    OPEN(nuin,FILE=TRIM(namelist_file), IOSTAT=ierr)
    IF (ierr /= 0) THEN
      WRITE(message_text,*)'Cannot open ', TRIM(namelist_file)
      CALL logging%error(message_text,__FILE__, __LINE__)
    ENDIF

    READ(nuin, NML=lu_raw_data, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      CALL logging%error('Cannot read in namelist lu_raw_data',__FILE__, __LINE__)
    ENDIF

    READ(nuin, NML=lu_io_extpar, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      CALL logging%error('Cannot read in namelist lu_io_data',__FILE__, __LINE__)
    ENDIF

    ! If optional argument is present for output, copy the value from the local variable to the output argument variable
    IF (PRESENT(raw_data_glcc_path_opt)) THEN
      READ(nuin, NML=glcc_raw_data, IOSTAT=ierr)
      IF (ierr /= 0) THEN
        CALL logging%error('Cannot read in namelist glcc_path_opt',__FILE__, __LINE__)
      ENDIF
      READ(nuin, NML=glcc_io_extpar, IOSTAT=ierr)
      IF (ierr /= 0) THEN
        CALL logging%error('Cannot read in namelist glcc_io_extpar',__FILE__, __LINE__)
      ENDIF
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

    IF (ntiles_globcover == 1) THEN
      nrow_tiles = 1
      ncolumn_tiles = 1
    ELSE IF (ntiles_globcover /= 6) THEN
      nrow_tiles = ntiles_globcover / ncolumn_tiles
      IF (ncolumn_tiles * nrow_tiles /= ntiles_globcover) THEN
        CALL logging%error('ncolumn_tiles specified does not fit to ntiles_globcover',__FILE__,__LINE__)
      ENDIF
    ENDIF

    ! prohibit use of corine for other cases thaan globcover
    IF (l_use_corine .AND. i_landuse_data /= 1 ) THEN
      CALL logging%error('Corine dataset can only be use in combination with Globcover dataset')
    ENDIF

  END SUBROUTINE read_namelists_extpar_land_use

  !> inquire dimension information for glc2000 raw data
  SUBROUTINE get_dimension_glc2000_data(path_glc2000_file, &
       &                                nlon_glc2000,          &
       &                                nlat_glc2000)


    CHARACTER (len=*), INTENT(in)  :: path_glc2000_file(:)         !< filename with path for glc2000 raw data
    INTEGER (KIND=i4), INTENT(out) :: nlon_glc2000, &  !< number of grid elements in zonal direction for glc2000 data
         &                            nlat_glc2000 !< number of grid elements in meridional direction for glc2000 data

    !local variables
    INTEGER(KIND=i4)               :: ncid, &                              !< netcdf unit file number
         &                            ndimension, &                        !< number of dimensions in netcdf file
         &                            nVars, &                             !< number of variables in netcdf file
         &                            nGlobalAtts, &                       !< number of gloabal Attributes in netcdf file
         &                            unlimdimid, &                        !< id of unlimited dimension (e.g. time) in netcdf file
         &                            dimid, &                             !< id of dimension
         &                            length                           !< length of dimension

    CHARACTER (len=80)             :: dimname               !< name of dimensiona

    ! open netcdf file
    CALL check_netcdf( nf90_open(TRIM(path_glc2000_file(1)),NF90_NOWRITE, ncid))

    ! look for numbers of dimensions, Variable, Attributes, and the dimid for the unlimited dimension (probably time)
    !; nf90_inquire input: ncid; nf90_inquire output: ndimension, nVars, nGlobalAtts,unlimdimid
    CALL check_netcdf (nf90_inquire(ncid,ndimension, nVars, nGlobalAtts,unlimdimid))

    !; the dimid in netcdf-files is counted from 1 to ndimension
    !; look for the name and length of the dimension with f90_inquire_dimension
    !; nf90_inquire_dimension input: ncid, dimid; nf90_inquire_dimension output: name, length
    DO dimid=1,ndimension
      CALL check_netcdf( nf90_inquire_dimension(ncid,dimid, dimname, length) )
      IF ( trim(dimname) == 'lon') nlon_glc2000=length          ! here I know that the name of zonal dimension is 'lon'
      IF ( trim(dimname) == 'lat') nlat_glc2000=length          ! here I know that the name of meridional dimension is 'lat'
    ENDDO

    ! close netcdf file
    CALL check_netcdf( nf90_close( ncid))

  END SUBROUTINE get_dimension_glc2000_data

  !> get coordinates for glc2000 raw data
  SUBROUTINE get_lonlat_glc2000_data(path_glc2000_file, &
       &                             nlon_glc2000,       &
       &                             nlat_glc2000,       &
       &                             lon_glc2000,        &
       &                             lat_glc2000,        &
       &                             glc2000_grid)


    CHARACTER (len=*), INTENT(in)      :: path_glc2000_file(:)         !< filename with path for glc2000 raw data

    INTEGER (KIND=i4), INTENT(in)      :: nlon_glc2000, &  !< number of grid elements in zonal direction for glc2000 data
         &                                nlat_glc2000 !< number of grid elements in meridional direction for glc2000 data

    REAL (KIND=wp), INTENT(out)        :: lon_glc2000(1:nlon_glc2000), &  !< longitude of glc2000 raw data
         &                                lat_glc2000(1:nlat_glc2000) !< latitude of glc2000 raw data

    TYPE(reg_lonlat_grid), INTENT(OUT) :: glc2000_grid !< structure with defenition of the raw data grid
                                                       !  for the whole GLC2000 dataset
    !local variables
    INTEGER(KIND=i4)                   :: ncid, varid
    CHARACTER (LEN=80)                 :: varname  !< name of variable

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

  !> inquire dimension information for glcc raw data
  SUBROUTINE get_dimension_glcc_data(path_glcc_file, &
                                      nlon_glcc,      &
                                      nlat_glcc)


    CHARACTER (len=*), INTENT(in)  :: path_glcc_file(1)         !< filename with path for glcc raw data
    INTEGER (KIND=i4), INTENT(out) :: nlon_glcc, &  !< number of grid elements in zonal direction for glcc data
         &                            nlat_glcc !< number of grid elements in meridional direction for glcc data

    !local variables
    INTEGER(KIND=i4)               :: ncid, &                              !< netcdf unit file number
         &                            ndimension, &                        !< number of dimensions in netcdf file
         &                            nVars, &                             !< number of variables in netcdf file
         &                            nGlobalAtts, &                       !< number of gloabal Attributes in netcdf file
         &                            unlimdimid, &                        !< id of unlimited dimension (e.g. time) in netcdf file
         &                            dimid, &                             !< id of dimension
         &                            length                           !< length of dimension

    CHARACTER (len=80)            :: dimname               !< name of dimensiona

    ! open netcdf file
    CALL check_netcdf( nf90_open(TRIM(path_glcc_file(1)),NF90_NOWRITE, ncid))

    ! look for numbers of dimensions, Variable, Attributes, and the dimid for the unlimited dimension (probably time)
    !; nf90_inquire input: ncid; nf90_inquire output: ndimension, nVars, nGlobalAtts,unlimdimid
    CALL check_netcdf (nf90_inquire(ncid,ndimension, nVars, nGlobalAtts,unlimdimid))

    !; the dimid in netcdf-files is counted from 1 to ndimension
    !; look for the name and length of the dimension with f90_inquire_dimension
    !; nf90_inquire_dimension input: ncid, dimid; nf90_inquire_dimension output: name, length
    DO dimid=1,ndimension
      CALL check_netcdf( nf90_inquire_dimension(ncid,dimid, dimname, length) )
      IF ( TRIM(dimname) == 'lon') nlon_glcc=length          ! here I know that the name of zonal dimension is 'lon'
      IF ( TRIM(dimname) == 'lat') nlat_glcc=length          ! here I know that the name of meridional dimension is 'lat'
    ENDDO

    ! close netcdf file
    CALL check_netcdf( nf90_close( ncid))

  END SUBROUTINE get_dimension_glcc_data

      !> get coordinates for glcc raw data
  SUBROUTINE get_lonlat_glcc_data(path_glcc_file,&
       &                              nlon_glcc, &
       &                              nlat_glcc, &
       &                              lon_glcc,  &
       &                              lat_glcc,  &
       &                              glcc_grid)

    CHARACTER (len=*), INTENT(in)     :: path_glcc_file(1)         !< filename with path for glcc raw data
    INTEGER (KIND=i4), INTENT(in)     :: nlon_glcc, &  !< number of grid elements in zonal direction for glcc data
         &                               nlat_glcc !< number of grid elements in meridional direction for glcc data

    REAL (KIND=wp), INTENT(out)       :: lon_glcc(1:nlon_glcc), &  !< longitude of glcc raw data
         &                               lat_glcc(1:nlat_glcc) !< latitude of glcc raw data

    TYPE(reg_lonlat_grid), INTENT(OUT):: glcc_grid !< structure with defenition of the raw data grid
                                                    !  for the whole GLCC dataset
    !local variables
    INTEGER(KIND=i4)                  :: ncid,varid

    CHARACTER (LEN=80)                :: varname  !< name of variable

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

    WRITE(message_text,'(a,2f18.12)') 'GLCC longitude bounds', glcc_grid%start_lon_reg, glcc_grid%end_lon_reg
    CALL logging%info(message_text)
    WRITE(message_text,'(a,2f18.12)') 'GLCC latitude bounds', glcc_grid%start_lat_reg, glcc_grid%end_lat_reg
    CALL logging%info(message_text)

  END SUBROUTINE get_lonlat_glcc_data

  !> inquire dimension information for globcover raw data
  SUBROUTINE get_dimension_globcover_data(nlon_globcover, &
                                    nlat_globcover)


    INTEGER (KIND=i4), INTENT(OUT) :: nlon_globcover, &  !< number of grid elements in zonal direction for globcover data
         &                            nlat_globcover !< number of grid elements in meridional direction for globcover data

    !local variables
    INTEGER(KIND=i4), PARAMETER    :: nx=43200, &
         &                            ny=27900

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

  SUBROUTINE get_dimension_ecoclimap_data(nlon_ecoclimap, &
       &                            nlat_ecoclimap)

    INTEGER (KIND=i4), INTENT(OUT) :: nlon_ecoclimap, &  !< number of grid elements in zonal direction for globcover data
         &                            nlat_ecoclimap !< number of grid elements in meridional direction for globcover data

    !local variables
    INTEGER(KIND=i4), PARAMETER    :: nx=43200, &
         &                            ny=21600

    nlon_ecoclimap = nx
    nlat_ecoclimap = ny

  END SUBROUTINE get_dimension_ecoclimap_data

  !> inquire dimension information for ecci raw data
  SUBROUTINE get_dimension_ecci_data(nlon_ecci, &
                                    nlat_ecci)


  INTEGER (KIND=i4), INTENT(OUT) :: nlon_ecci, & !< number of grid elements in zonal direction for ecci data
       &                            nlat_ecci !< number of grid elements in meridional direction for ecci data

  !local variables
  INTEGER(KIND=i4), PARAMETER    :: nx=43200, &
       &                            ny=32400

  IF(ntiles_ecci == 1) THEN
    nlon_ecci = len_lu_lon_ecci
    nlat_ecci = len_lu_lat_ecci
  ELSE IF(ntiles_ecci == 6) THEN
    nlon_ecci = 3 * nx
    nlat_ecci = 2 * ny
  ELSE
    nlon_ecci = ncolumn_tiles_ecci * nx
    nlat_ecci = nrow_tiles_ecci * ny
  END IF

 END SUBROUTINE get_dimension_ecci_data

  !> get coordinates for globcover raw data
  SUBROUTINE get_lonlat_globcover_data(nlon_globcover, &
       &                               nlat_globcover, &
       &                               lon_globcover,  &
       &                               lat_globcover,  &
       &                               globcover_grid)


    INTEGER (KIND=i4), INTENT(IN)      :: nlon_globcover, &  !< number of grid elements in zonal direction for globcover data
         &                                nlat_globcover !< number of grid elements in meridional direction for globcover data

    REAL (KIND=wp), INTENT(OUT)        :: lon_globcover(1:nlon_globcover), &  !< longitude of globcover raw data
         &                                lat_globcover(1:nlat_globcover) !< latitude of globcover raw data

    TYPE(reg_lonlat_grid), INTENT(OUT) :: globcover_grid !< structure with defenition of the raw data grid for

    REAL (KIND=wp)                     :: xmin_glc, &  ! area of glcover data: western longitude
         &                                xmax_glc, &  ! area of glcover data: eastern longitude
         &                                ymax_glc, &  ! area of glcover data: northern latitude
         &                                ymin_glc, &  ! area of glcover data: southern latitude
         &                                dx_glc, &   ! grid element size of glcover data pixel in zonal direction
         &                                dy_glc  ! grid element size of glcover data pixel in meridional directionon

    INTEGER (KIND=i4)                  :: jx,jy


    xmin_glc = MINVAL(lu_tiles_lon_min)
    xmax_glc = MAXVAL(lu_tiles_lon_max)
    ymax_glc = MAXVAL(lu_tiles_lat_max)
    ymin_glc = MINVAL(lu_tiles_lat_min)
    dx_glc   = (xmax_glc - xmin_glc)/ REAL(nlon_globcover)
    dy_glc   = -1.0 * (ymax_glc - ymin_glc) / REAL(nlat_globcover)

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
    globcover_grid%dlon_reg      = dx_glc ! (lon_globcover(nlon_globcover) - lon_globcover(1)) / (nlon_globcover - 1)
    globcover_grid%dlat_reg      = dy_glc ! (lat_globcover(nlat_globcover) - lat_globcover(1)) / (nlat_globcover - 1)
    globcover_grid%nlon_reg      = nlon_globcover
    globcover_grid%nlat_reg      = nlat_globcover

  END SUBROUTINE get_lonlat_globcover_data

  !> get coordinates for ecoclimap raw data
  SUBROUTINE get_lonlat_ecoclimap_data(nlon_ecoclimap, &
       &                               nlat_ecoclimap, &
       &                               lon_ecoclimap,  &
       &                               lat_ecoclimap,  &
       &                               ecoclimap_grid)


    INTEGER (KIND=i4), INTENT(IN)      :: nlon_ecoclimap, &  !< number of grid elements in zonal direction for ecoclimap data
         &                                nlat_ecoclimap !< number of grid elements in meridional direction for ecoclimap data

    REAL (KIND=wp), INTENT(OUT)        :: lon_ecoclimap(1:nlon_ecoclimap), &  !< longitude of ecoclimap raw data
         &                                lat_ecoclimap(1:nlat_ecoclimap) !< latitude of ecoclimap raw data

    TYPE(reg_lonlat_grid), INTENT(OUT) :: ecoclimap_grid !< structure with defenition of the raw data grid

    !local variables
    REAL(KIND=i4), PARAMETER          ::  xmin_glc = -179.995833335, &     ! area of glcover data: western longitude
         &                                xmax_glc  = 179.995689334969, &  ! area of glcover data: eastern longitude
         &                                ymax_glc  = 89.995761335, &      ! area of glcover data: northern latitude
         &                                ymin_glc = -89.995833334, &     ! area of glcover data: southern latitude
         &                                dx_glc =    0.0083333333, &   ! grid element size of glcover data pixel in zonal direction
         &                                dy_glc =   -0.0083333333  ! grid element size of glcover data pixel in meridional directionon

    INTEGER (KIND=i4)                  :: jx,jy

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

  !> get coordinates for ecci raw data
  SUBROUTINE get_lonlat_ecci_data(nlon_ecci, &
     &                                 nlat_ecci, &
     &                                 lon_ecci,  &
     &                                 lat_ecci,  &
     &                                 ecci_grid)


    INTEGER (KIND=i4), INTENT(IN)      :: nlon_ecci, & !< number of grid elements in zonal direction for ecci data
         &                                nlat_ecci !< number of grid elements in meridional direction for ecci data
    REAL (KIND=wp), INTENT(OUT)        :: lon_ecci(1:nlon_ecci), & !< longitude of ecci raw data
         &                                lat_ecci(1:nlat_ecci) !< latitude of ecci raw data

    TYPE(reg_lonlat_grid), INTENT(OUT) :: ecci_grid !< structure with defenition of the raw data grid for
                                                         !  the whole GLOBECOVER dataset
    ! local variables
    REAL (KIND=wp)                     :: xmin_ecci, & ! area of glcover data: western longitude
         &                                xmax_ecci, & ! area of glcover data: eastern longitude
         &                                ymax_ecci, & ! area of glcover data: northern latitude
         &                                ymin_ecci, & ! area of glcover data: southern latitude
         &                                dx_ecci, &  ! grid element size of glcover data pixel in zonal direction
         &                                dy_ecci  ! grid element size of glcover data pixel in meridional directionon

    INTEGER (KIND=i4)                  :: jx,jy

    xmin_ecci = MINVAL(lu_tiles_lon_min_ecci)
    xmax_ecci = MAXVAL(lu_tiles_lon_max_ecci)
    ymax_ecci = MAXVAL(lu_tiles_lat_max_ecci)
    ymin_ecci = MINVAL(lu_tiles_lat_min_ecci)
    dx_ecci   = (xmax_ecci - xmin_ecci)/ REAL(nlon_ecci) !_br 04.04.14
    dy_ecci   = -1.0_wp * (ymax_ecci - ymin_ecci) / REAL(nlat_ecci) !_br 04.04.14

    DO jx=1,nlon_ecci
       lon_ecci(jx)  = xmin_ecci + 0.5_wp*dx_ecci + (jx-1)*dx_ecci
    ENDDO
    DO jy=1,nlat_ecci
     lat_ecci(jy) = ymax_ecci + 0.5_wp*dy_ecci + (jy-1)*dy_ecci !note negative increment!
    ENDDO

    ! define the values for the structure ecci_grid
    ecci_grid%start_lon_reg = lon_ecci(1)
    ecci_grid%end_lon_reg   = lon_ecci(nlon_ecci)
    ecci_grid%start_lat_reg = lat_ecci(1)
    ecci_grid%end_lat_reg   = lat_ecci(nlat_ecci)
    ecci_grid%dlon_reg      = dx_ecci ! (lon_ecci(nlon_ecci) - lon_ecci(1)) / (nlon_ecci - 1)
    ecci_grid%dlat_reg      = dy_ecci ! (lat_ecci(nlat_ecci) - lat_ecci(1)) / (nlat_ecci - 1)
    ecci_grid%nlon_reg      = nlon_ecci
    ecci_grid%nlat_reg      = nlat_ecci

  END SUBROUTINE get_lonlat_ecci_data

  SUBROUTINE get_globcover_tiles_grid(globcover_tiles_grid)

    TYPE(reg_lonlat_grid), INTENT(OUT):: globcover_tiles_grid(1:ntiles_globcover)

    INTEGER(KIND=i4)                  ::  k      !counter

    REAL(KIND=wp)                     :: dlon, dlat

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
 
  SUBROUTINE get_ecci_tiles_grid(ecci_tiles_grid)

    TYPE(reg_lonlat_grid), INTENT(OUT):: ecci_tiles_grid(1:ntiles_ecci)

    INTEGER::  k      !counter

    REAL(KIND=wp) :: dlon
    REAL(KIND=wp) :: dlat

    DO k = 1,ntiles_ecci

      dlon = (lu_tiles_lon_max_ecci(k) - lu_tiles_lon_min_ecci(k)) / REAL(lu_tiles_ncolumns_ecci(k)) !_br 04.04.14
      dlat = -1._wp * (lu_tiles_lat_max_ecci(k) - lu_tiles_lat_min_ecci(k)) / REAL(lu_tiles_nrows_ecci(k))   
      !< latitude from north to south, negative increment !_br 04.04.14

      ecci_tiles_grid(k)%start_lon_reg = lu_tiles_lon_min_ecci(k) + 0.5_wp*dlon
      ecci_tiles_grid(k)%end_lon_reg = lu_tiles_lon_max_ecci(k) - 0.5_wp*dlon 

      ecci_tiles_grid(k)%start_lat_reg = lu_tiles_lat_max_ecci(k) + 0.5_wp*dlat
      ecci_tiles_grid(k)%end_lat_reg = lu_tiles_lat_min_ecci(k) - 0.5_wp*dlat 
      ecci_tiles_grid(k)%dlon_reg = dlon
      ecci_tiles_grid(k)%dlat_reg = dlat
      ecci_tiles_grid(k)%nlon_reg = lu_tiles_ncolumns_ecci(k)
      ecci_tiles_grid(k)%nlat_reg = lu_tiles_nrows_ecci(k)

    END DO

  END SUBROUTINE get_ecci_tiles_grid

  !> determine grid description of band for GLOBCOVER
  SUBROUTINE det_band_globcover_data(globcover_grid,start_globcover_row,ta_grid)

    TYPE(reg_lonlat_grid),INTENT(IN)  :: globcover_grid ! sturcture with the definition of the global data grid
    INTEGER, INTENT(IN)               :: start_globcover_row         ! number of the start row of band
    TYPE(reg_lonlat_grid), INTENT(OUT):: ta_grid       ! structure with definition of the target area grid.

    INTEGER(KIND=i4)                  :: nrows = 2500              ! number of rows, set to 2500 as default
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

  !> determine grid description of band for ECCI
  SUBROUTINE det_band_ecci_data(ecci_grid,start_ecci_row,ta_grid)

    TYPE(reg_lonlat_grid),INTENT(IN)  :: ecci_grid ! sturcture with the definition of the global data grid
    INTEGER(KIND=i4), INTENT(IN)      :: start_ecci_row         ! number of the start row of band 
    TYPE(reg_lonlat_grid), INTENT(OUT):: ta_grid       ! structure with definition of the target area grid.

    ! local variables
    INTEGER(KIND=i4)                  :: nrows = 2500              ! number of rows, set to 2500 as default

    ! band from east to west for the whole globe, like the complete ecci grid
    ta_grid%dlon_reg = ecci_grid%dlon_reg
    ta_grid%dlat_reg = ecci_grid%dlat_reg

    ta_grid%start_lon_reg = ecci_grid%start_lon_reg
    ta_grid%end_lon_reg = ecci_grid%end_lon_reg
    ta_grid%nlon_reg = ecci_grid%nlon_reg

    !latitude from north to south, negative increment
    ta_grid%nlat_reg = nrows
    ta_grid%start_lat_reg = ecci_grid%start_lat_reg + ta_grid%dlat_reg * (start_ecci_row-1._wp)   
    !< latitude from north to south, note the negative increment!
    ta_grid%end_lat_reg = ta_grid%start_lat_reg + ta_grid%dlat_reg * (nrows - 1._wp)   
    !< latitude from north to south, note the negative increment!
    
    ! check for the southern bound of the ecci data
    IF (ta_grid%end_lat_reg < ecci_grid%end_lat_reg) THEN ! band is at the southern bound
      ta_grid%end_lat_reg = ecci_grid%end_lat_reg
      ta_grid%nlat_reg    = NINT(((ta_grid%end_lat_reg - ta_grid%start_lat_reg) / ta_grid%dlat_reg)) + 1
    END IF

  END SUBROUTINE det_band_ecci_data

  ! get GLOBCOVER data block for a given target area from the tile block indices
  SUBROUTINE get_globcover_data_block(ta_grid,              &
       &                              globcover_tiles_grid, &
       &                              ncids_globcover,      &
       &                              lu_block)

    TYPE(reg_lonlat_grid), INTENT(IN) :: ta_grid, &  !< structure with definition of the target area grid
         &                               globcover_tiles_grid(1:ntiles_globcover)

    !< structure with defenition of the raw data grid for the 16 GLOBECOVER tiles
    INTEGER(KIND=i4) , INTENT(IN)     :: ncids_globcover(1:ntiles_globcover)
    !< ncid for the GLOBCOVER tiles, the netcdf files have to be opened previously
    INTEGER (KIND=i2), INTENT(OUT)    :: lu_block(1:ta_grid%nlon_reg,1:ta_grid%nlat_reg) !< a block of GLOBCOVER data

    !local variables
    INTEGER (KIND=i4)                 :: globcover_startrow(1:ntiles_globcover), &  !< startrow indices for each GLOBCOVER tile
         &                               globcover_endrow(1:ntiles_globcover), &  !< endrow indices for each GLOBCOVER tile
         &                               globcover_startcolumn(1:ntiles_globcover), &  !< starcolumn indices for each GLOBCOVER tile
         &                               globcover_endcolumn(1:ntiles_globcover), &  !< endcolumn indices for each GLOBCOVER tile
         &                               ta_start_ie(1:ntiles_globcover), &
         &                               ta_end_ie(1:ntiles_globcover), &
         &                               ta_start_je(1:ntiles_globcover), &
         &                               ta_end_je(1:ntiles_globcover)


    INTEGER(KIND=i4)                  :: varid, nrows, ncolumns, k ,errorcode
    INTEGER (KIND=i2), ALLOCATABLE    :: raw_lu_block(:,:) !< a block with GLOBCOVER data
    CHARACTER (LEN=80)                :: varname  !< name of variable

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
        IF(errorcode/=0) CALL logging%error('Cant allocate the array raw_lu_block',__FILE__,__LINE__)

        CALL check_netcdf(nf90_inq_varid(ncids_globcover(k),TRIM(varname),varid)) ! get the varid of the altitude variable
        ! get the data into the raw_lu_block
        CALL check_netcdf(nf90_get_var(ncids_globcover(k), varid,  raw_lu_block,     &
        &     start=(/globcover_startcolumn(k),globcover_startrow(k)/),count=(/ncolumns,nrows/)))


        lu_block(ta_start_ie(k):ta_end_ie(k),ta_start_je(k):ta_end_je(k)) = raw_lu_block(1:ncolumns,1:nrows)

        DEALLOCATE (raw_lu_block, STAT=errorcode)
        IF(errorcode/=0) CALL logging%error('Cant deallocate the array raw_lu_block',__FILE__,__LINE__)

      ENDIF
    ENDDO

  END SUBROUTINE get_globcover_data_block
 
  ! get ECCI data block for a given target area from the tile block indices
  SUBROUTINE get_ecci_data_block(ta_grid,              &
                                 ecci_tiles_grid, &
                                 ncids_ecci,      &
                                 lu_block)

    TYPE(reg_lonlat_grid), INTENT(IN)  :: ta_grid !< structure with definition of the target area grid
                                                  !  (dlon must be the same as for the whole ECCI dataset)
    TYPE(reg_lonlat_grid), INTENT(IN)  :: ecci_tiles_grid(1:ntiles_ecci)
    !< structure with defenition of the raw data grid for the 16 GLOBECOVER tiles
    INTEGER(KIND=i4) , INTENT(IN)      :: ncids_ecci(1:ntiles_ecci)  
    !< ncid for the ECCI tiles, the netcdf files have to be opened previously
    INTEGER (KIND=i2), INTENT(OUT)     :: lu_block(1:ta_grid%nlon_reg,1:ta_grid%nlat_reg) !< a block of ECCI data 

    !local variables
    INTEGER (KIND=i4)                  :: ecci_startrow(1:ntiles_ecci), & !< startrow indices for each ECCI tile
         &                                ecci_endrow(1:ntiles_ecci), & !< endrow indices for each ECCI tile
         &                                ecci_startcolumn(1:ntiles_ecci), & !< starcolumn indices for each ECCI tile
         &                                ecci_endcolumn(1:ntiles_ecci), & !< endcolumn indices for each ECCI tile
         &                                ta_start_ie(1:ntiles_ecci), &
         &                                ta_end_ie(1:ntiles_ecci), &
         &                                ta_start_je(1:ntiles_ecci), &
         &                                ta_end_je(1:ntiles_ecci)


    INTEGER (KIND=i2), ALLOCATABLE     :: raw_lu_block(:,:) !< a block with ECCI data
    INTEGER(KIND=i4)                   :: varid, nrows, ncolumns, &
         &                                k, errorcode
    CHARACTER (LEN=80)                 :: varname  !< name of variable

    varname = 'lccs_class'   ! I know that in the ECCI netcdf files the LU data is stored in a variable "ECCI"

    CALL get_ecci_tile_block_indices(ta_grid,              &
         &                                ecci_tiles_grid, &  
         &                                ecci_startrow,   &
         &                                ecci_endrow,     & 
         &                                ecci_startcolumn,&
         &                                ecci_endcolumn,  &
         &                                ta_start_ie,          & 
         &                                ta_end_ie,            &
         &                                ta_start_je,          &
         &                                ta_end_je)

    DO k = 1, ntiles_ecci
      IF ((ecci_startrow(k)/=0).AND.(ecci_startcolumn(k)/=0)) THEN
        nrows = ecci_endrow(k) - ecci_startrow(k) + 1
        ncolumns = ecci_endcolumn(k) - ecci_startcolumn(k) + 1
      
        ALLOCATE (raw_lu_block(1:ncolumns,1:nrows), STAT=errorcode)
        IF(errorcode/=0) CALL logging%error('Cant allocate the array raw_lu_block',__FILE__,__LINE__)
      
        CALL check_netcdf(nf90_inq_varid(ncids_ecci(k),TRIM(varname),varid)) ! get the varid of the altitude variable
        ! get the data into the raw_lu_block
        CALL check_netcdf(nf90_get_var(ncids_ecci(k), varid,  raw_lu_block,     & 
        &     start=(/ecci_startcolumn(k),ecci_startrow(k)/),count=(/ncolumns,nrows/)))
      
        lu_block(ta_start_ie(k):ta_end_ie(k),ta_start_je(k):ta_end_je(k)) = raw_lu_block(1:ncolumns,1:nrows)
        
        DEALLOCATE (raw_lu_block, STAT=errorcode)
        IF(errorcode/=0) CALL logging%error('Cant deallocate the array raw_lu_block',__FILE__,__LINE__) 
      ENDIF
    ENDDO

  END SUBROUTINE get_ecci_data_block
 
  ! get startrow, endrow, startcolumn and endcolumn of each GLOBCOVER tile (raw data) for a 
  ! given target area (ta_grid) and get start_indices (lon, lat) and end_indices of the target
  ! area for each GLOBCOVER tile
  ! The GLOBCOVER raw data are split in 6 tiles, so the target area may overlap several tiles.
  ! This subroutine determines the necessary indices to read in the GLOBCOVER data into the
  ! target area.
  ! GLOBCOVER tiles which are outside the target block will get indices with the value '0'

  SUBROUTINE get_globcover_tile_block_indices(ta_grid,              &
       &                                      globcover_tiles_grid, &
       &                                      globcover_startrow,   &
       &                                      globcover_endrow,     &
       &                                      globcover_startcolumn,&
       &                                      globcover_endcolumn,  &
       &                                      ta_start_ie,          &
       &                                      ta_end_ie,            &
       &                                      ta_start_je,          &
       &                                      ta_end_je)

    TYPE(reg_lonlat_grid), INTENT(IN) :: ta_grid, &  !< structure with definition of the target area grid
         &                               globcover_tiles_grid(1:ntiles_globcover)

    INTEGER (KIND=i4), INTENT(OUT)    :: globcover_startrow(1:ntiles_globcover), &
         &                               globcover_endrow(1:ntiles_globcover), &
         &                               globcover_startcolumn(1:ntiles_globcover), &
         &                               globcover_endcolumn(1:ntiles_globcover), &
         &                               ta_start_ie(1:ntiles_globcover), &
         &                               ta_end_ie(1:ntiles_globcover), &
         &                               ta_start_je(1:ntiles_globcover), &
         &                               ta_end_je(1:ntiles_globcover)

    ! local variables

    INTEGER (KIND=i4)                :: startrow, &  ! startrow for tile
         &                              endrow, &
         &                              startcolumn, &
         &                              endcolumn, &
         &                              undefined, &
         &                              k

    REAL (KIND=wp)                  :: dlon, dlat

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

       ! get startrow, endrow, startcolumn and endcolumn of each ECCI tile (raw data) for a 
       ! given target area (ta_grid) and get start_indices (lon, lat) and end_indices of the target
       ! area for each ECCI tile
       ! The ECCI raw data are split in 6 tiles, so the target area may overlap several tiles.
       ! This subroutine determines the necessary indices to read in the ECCI data into the
       ! target area.
       ! ECCI tiles which are outside the target block will get indices with the value '0'

  SUBROUTINE get_ecci_tile_block_indices(ta_grid,              &
           &                                  ecci_tiles_grid, &
           &                                  ecci_startrow,   &
           &                                  ecci_endrow,     & 
           &                                  ecci_startcolumn,&
           &                                  ecci_endcolumn,  &
           &                                  ta_start_ie,          &
           &                                  ta_end_ie,            &
           &                                  ta_start_je,          &
           &                                  ta_end_je)


    TYPE(reg_lonlat_grid), INTENT(IN) :: ta_grid, & !< structure with definition of the target area grid
         &                               ecci_tiles_grid(1:ntiles_ecci) 

    INTEGER (KIND=i4), INTENT(OUT)    :: ecci_startrow(1:ntiles_ecci), &
         &                               ecci_endrow(1:ntiles_ecci), &
         &                               ecci_startcolumn(1:ntiles_ecci), &
         &                               ecci_endcolumn(1:ntiles_ecci), &
         &                               ta_start_ie(1:ntiles_ecci), &    
         &                               ta_end_ie(1:ntiles_ecci), &      
         &                               ta_start_je(1:ntiles_ecci), &  
         &                               ta_end_je(1:ntiles_ecci)   

    

    ! local variables
    INTEGER (KIND=i4)                 :: k, & 
         &                               undefined, &
         &                               startrow, & ! startrow for tile
         &                               endrow, & 
         &                               startcolumn, &
         &                               endcolumn

    REAL (KIND=wp)                    :: dlon, &
         &                               dlat


    undefined = 0
    ecci_startrow     = undefined
    ecci_endrow       = undefined
    ecci_startcolumn  = undefined
    ecci_endcolumn    = undefined
    ta_start_ie = undefined 
    ta_end_ie   = undefined
    ta_start_je = undefined
    ta_end_je   = undefined

    k=1                      ! determin dlon and dlat (are the same for all tiles)
    dlon = ta_grid%dlon_reg
    dlat = ta_grid%dlat_reg

    DO k = 1,ntiles_ecci   !loop over the tiles which overlap the target area

      startcolumn = NINT((ta_grid%start_lon_reg - ecci_tiles_grid(k)%start_lon_reg)/dlon) + 1
      ! here I want nearest index (NINT)

      IF (startcolumn < 1) THEN 
        ecci_startcolumn(k) = 1
        ! get the start index of the subtile for the target area block
        ta_start_ie(k) = NINT ((ecci_tiles_grid(k)%start_lon_reg - ta_grid%start_lon_reg)/dlon) + 1
        ! index of target area block

      ELSE IF (startcolumn > lu_tiles_ncolumns_ecci(k)) THEN
        ecci_startcolumn(k) = 0
        ta_start_ie(k) = 0
      ELSE
        ecci_startcolumn(k) = startcolumn
        ta_start_ie(k) = 1
      ENDIF

      ! get endcolumn for tile k
      endcolumn = NINT((ta_grid%end_lon_reg - ecci_tiles_grid(k)%start_lon_reg)/dlon) +1

      IF (endcolumn > lu_tiles_ncolumns_ecci(k)) THEN 
        ecci_endcolumn(k) = lu_tiles_ncolumns_ecci(k)
        ! get the end index of the subtile for the target area block
        ta_end_ie(k) = NINT ((ecci_tiles_grid(k)%end_lon_reg - ta_grid%start_lon_reg)/dlon) + 1
        ! index of target area block
      ELSE IF (endcolumn < 1) THEN
        ecci_endcolumn(k) = 0
        ta_end_ie(k) = 0
      ELSE
        ecci_endcolumn(k) = endcolumn
        ta_end_ie(k) = ta_grid%nlon_reg
      ENDIF

      ! get startrow for tile k
      startrow = NINT((ta_grid%start_lat_reg - ecci_tiles_grid(k)%start_lat_reg)/dlat) + 1
     
      IF (startrow < 1) THEN 
        ecci_startrow(k) = 1
        ! get the start index of the subtile for the target area block
        ta_start_je(k) = NINT ((ecci_tiles_grid(k)%start_lat_reg  - ta_grid%start_lat_reg)/dlat) + 1
        ! index of target area block
      ELSE IF (startrow > lu_tiles_nrows_ecci(k)) THEN
        ecci_startrow(k) = 0
        ta_start_je(k) = 0
      ELSE
        ecci_startrow(k) = startrow
        ta_start_je(k) = 1
      ENDIF
      
      ! get endrow for tile k
      endrow   = NINT(( ta_grid%end_lat_reg - ecci_tiles_grid(k)%start_lat_reg )/dlat)  + 1
     
      IF (endrow > lu_tiles_nrows_ecci(k)) THEN 
        ecci_endrow(k) = lu_tiles_nrows_ecci(k)
        ! get the start index of the subtile for the target area block
        ta_end_je(k) = NINT ((ecci_tiles_grid(k)%end_lat_reg -  ta_grid%start_lat_reg )/dlat) + 1
        ! index of target area block
      
      ELSE IF (endrow < 1) THEN
        ecci_endrow(k) = 0
        ta_end_je(k) = 0
      ELSE
        ecci_endrow(k) = endrow
        ta_end_je(k) =  ta_grid%nlat_reg
      ENDIF
 
    ENDDO  ! loop over the tiles 

  END SUBROUTINE get_ecci_tile_block_indices

END MODULE mo_landuse_routines
