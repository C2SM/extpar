!+ Fortran module with NDVI data handling routines
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_0         2010/12/21 Hermann Asensio
!  Initial release
! V2_0_3       2015-02-23 Juergen Helmert, Daniel Luethi
!  Increase working precision in grid size computation         
!  compute raw data lat/lon values (values on file are inaccurate)
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!> Fortran module with NDVI data handling routines
!> \author Hermann Asensio
!>
MODULE mo_ndvi_routines

  USE mo_logging
  USE mo_kind,                  ONLY: wp,i4

  USE mo_io_units,              ONLY: filename_max

  USE netcdf,                   ONLY :   &
                                      nf90_open,              &
                                      nf90_close,             &
                                      nf90_inquire,           &
                                      nf90_inquire_dimension, &
                                      nf90_inquire_variable,  &
                                      nf90_inq_attname,       &
                                      nf90_inquire_attribute, &
                                      nf90_get_att,           &
                                      nf90_inquire_dimension, &
                                      nf90_get_var,           &
                                      nf90_char,              &
                                      nf90_nowrite

  USE mo_io_utilities,           ONLY: check_netcdf

  USE mo_utilities_extpar,       ONLY: free_un ! function to get free unit number

  IMPLICIT NONE

  PRIVATE
  !
  PUBLIC :: open_netcdf_NDVI_data, &
       &    close_netcdf_NDVI_data, &
       &    read_namelists_extpar_ndvi, &
       &    read_ndvi_data_input_namelist, &
       &    get_dimension_NDVI_data, &
       &    get_NDVI_data_coordinates, &
       &    get_one_row_NDVI_data, &
       &    get_block_NDVI_data, &
       &    get_pixel_NDVI_data

  CONTAINS

  !> subroutine to read namelist for NDVI data settings for EXTPAR 
  SUBROUTINE read_namelists_extpar_ndvi(namelist_file, &
       &                                raw_data_ndvi_path, &
       &                                raw_data_ndvi_filename, &
       &                                ndvi_buffer_file, &
       &                                ndvi_output_file)



    CHARACTER (len=*), INTENT(IN) :: namelist_file !< filename with namelists for for EXTPAR settings

    CHARACTER (len=filename_max)             :: raw_data_ndvi_path, &        !< path to raw data
         &                                      raw_data_ndvi_filename, & !< filename NDVI raw data
         &                                      ndvi_buffer_file, & !< name for NDVI buffer file
         &                                      ndvi_output_file !< name for NDVI output file
       
    INTEGER (KIND=i4)                        :: ierr, nuin

    !> namelist with filenames for NDVI data input
    NAMELIST /ndvi_raw_data/ raw_data_ndvi_path, raw_data_ndvi_filename
    !> namelist with filenames for NDVI data output
    NAMELIST /ndvi_io_extpar/ ndvi_buffer_file, ndvi_output_file

    nuin = free_un()  ! functioin free_un returns free Fortran unit number
    OPEN(nuin,FILE=TRIM(namelist_file), IOSTAT=ierr)
    IF (ierr /= 0) THEN
      WRITE(message_text,*)'Cannot open ', TRIM(namelist_file)
      CALL logging%error(message_text,__FILE__, __LINE__) 
    ENDIF

    READ(nuin, NML=ndvi_raw_data, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      CALL logging%error('Cannot read in namelist ndvi_raw_data',__FILE__, __LINE__) 
    ENDIF

    READ(nuin, NML=ndvi_io_extpar, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      CALL logging%error('Cannot read in namelist ndvi_io_extpar',__FILE__, __LINE__) 
    ENDIF

    CLOSE(nuin)

  END SUBROUTINE read_namelists_extpar_ndvi

  !> open netcdf-file and get netcdf unit file number
  SUBROUTINE open_netcdf_NDVI_data(path_ndvi_file, ncid) 
            
    CHARACTER (len=*), INTENT(IN) :: path_ndvi_file         !< filename with path for NDVI raw data
    INTEGER, INTENT(OUT)          :: ncid                             !< netcdf unit file number

    !! open netcdf file 
    CALL check_netcdf( nf90_open(TRIM(path_ndvi_file),NF90_NOWRITE, ncid))

  END SUBROUTINE open_netcdf_NDVI_data

  !> close netcdf-file 
  SUBROUTINE close_netcdf_NDVI_data(ncid)

    INTEGER(KIND=i4), INTENT(in) :: ncid                             !< netcdf unit file number

    !! close netcdf file 
    CALL check_netcdf( nf90_close( ncid))

  END SUBROUTINE close_netcdf_NDVI_data

  !> read namelist with settings for COSMO target grid
  !> \author Hermann Asensio
  SUBROUTINE read_ndvi_data_input_namelist(input_namelist_file,        &
       &                                   raw_data_path,             &
       &                                   raw_data_ndvi_filename,    &
       &                                   outputgrid_ndvi_filename)

    CHARACTER (len=*), INTENT(in)             :: input_namelist_file !< file with input namelist 

    CHARACTER (len=filename_max), INTENT(out) :: raw_data_path, &        !< path to raw data
         &                                       raw_data_ndvi_filename, & !< filename ndvi raw data
         &                                       outputgrid_ndvi_filename !< output filename

    INTEGER (KIND=i4) :: ierr, nuin 

    !>Define the namelist group
    NAMELIST /ndvi_data_input/ raw_data_path, raw_data_ndvi_filename, outputgrid_ndvi_filename

    nuin = free_un()  ! functioin free_un returns free Fortran unit number

    OPEN(nuin,FILE=input_namelist_file, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      WRITE(message_text,*)'Cannot open ', TRIM(input_namelist_file)
      CALL logging%error(message_text,__FILE__, __LINE__) 
    ENDIF

    READ(nuin, NML=ndvi_data_input, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      CALL logging%error('Cannot read in namelist ndvi_data_input',__FILE__, __LINE__) 
    ENDIF

    CLOSE(nuin)

  END SUBROUTINE read_ndvi_data_input_namelist

  !> inquire dimension information for NDVI raw data 
  SUBROUTINE get_dimension_NDVI_data(ncid, &
       &                             nlon_ndvi, &
       &                             nlat_ndvi, &
       &                             ntime_ndvi)

    INTEGER(KIND=i4), INTENT(in)   :: ncid                             !< netcdf unit file number

    INTEGER (KIND=i4), INTENT(out) :: nlon_ndvi, & !< number of grid elements in zonal direction for NDVI data
         &                            nlat_ndvi, & !< number of grid elements in meridional direction for NDVI data
         &                            ntime_ndvi !< number of dates with NDVI data

    !local variables
    INTEGER(KIND=i4)               :: ndimension, &    !< number of dimensions in netcdf file
         &                            nVars, &         !< number of variables in netcdf file
         &                            nGlobalAtts, &   !< number of gloabal Attributes in netcdf file
         &                            unlimdimid, &    !< id of unlimited dimension (e.g. time) in netcdf file
         &                            dimid, &         !< id of dimension
         &                            length        !< length of dimension

    CHARACTER (len=80)             :: dimname               !< name of dimensiona


    !! look for numbers of dimensions, Variable, Attributes, and the dimid for the unlimited dimension (probably time)
    !; nf90_inquire input: ncid; nf90_inquire output: ndimension, nVars, nGlobalAtts,unlimdimid
    CALL check_netcdf (nf90_inquire(ncid,ndimension, nVars, nGlobalAtts,unlimdimid))

    !; the dimid in netcdf-files is counted from 1 to ndimension
    !; look for the name and length of the dimension with f90_inquire_dimension
    !; nf90_inquire_dimension input: ncid, dimid; nf90_inquire_dimension output: name, length
    DO dimid=1,ndimension
      CALL check_netcdf( nf90_inquire_dimension(ncid,dimid, dimname, length) )
      IF ( TRIM(dimname) == 'lon') nlon_ndvi=length          ! here I know that the name of zonal dimension is 'lon'
      IF ( TRIM(dimname) == 'lat') nlat_ndvi=length          ! here I know that the name of meridional dimension is 'lat'
      IF ( TRIM(dimname) == 'time') ntime_ndvi=length        ! here I know that the name of time dimension is 'time'
    ENDDO

  END SUBROUTINE get_dimension_NDVI_data

  !> read coordinates for NDVI raw data from netcdf file
  SUBROUTINE get_NDVI_data_coordinates(ncid,          &
       &                               nlon_ndvi,     &
       &                               nlat_ndvi,     &  
       &                               startlon_ndvi, &
       &                               startlat_ndvi, &
       &                               dlon_ndvi,     &
       &                               dlat_ndvi,     &
       &                               lon,           &
       &                               lat)

    INTEGER (KIND=i4), INTENT(IN):: nlon_ndvi, & !< number of grid elements in zonal direction for NDVI data
         &                          nlat_ndvi, & !< number of grid elements in meridional direction for NDVI data
         &                          ncid                             !< netcdf unit file number

    REAL (KIND=wp), INTENT(OUT)  :: startlon_ndvi, & !< longitude of lower left grid element for NDVI data 
         &                          startlat_ndvi, & !< latitude of lower left grid element for NDVI data
         &                          dlon_ndvi, & !< grid point distance in zonal direction (in degrees) for NDVI data
         &                          dlat_ndvi, & !< grid point distance in meridional direction (in degrees) for NDVI data
         &                          lon(1:nlon_ndvi), &      !< longitude of ndvi raw data in geographical system
         &                          lat(1:nlat_ndvi)      !< latitude of ndvi raw date in geographical system

    !local variables
    INTEGER(KIND=i4)             :: ndimension, &                       !< number of dimensions in netcdf file
         &                          nVars, &                            !< number of variables in netcdf file
         &                          nGlobalAtts, &                      !< number of gloabal Attributes in netcdf file
         &                          unlimdimid, &                       !< id of unlimited dimension (e.g. time) in netcdf file
         &                          length, &                           !< length of dimension
         &                          varid, &                            !< id of variable
         &                          xtype, &                            !< netcdf type of variable/attribute
         &                          ndim, &                             !< number of dimensions of variable
         &                          errorcode, &                        !< error status variable
         &                          attnum, &                           !< counter for attribute number
         &                          i, &                                !< loop index
         &                          nAtts                            !< number of attributes for a netcdf variable

    INTEGER(KIND=i4),ALLOCATABLE:: var_dimids(:)       !< id of variable dimensions, vector, maximal dimension ndimension

    CHARACTER (len=80)          :: dimname, &               !< name of dimension
         &                         varname, &               !< name of variable
         &                         attname, &               !< name of attribute
         &                         attributetext         !< attribute text

    !! look for numbers of dimensions, Variable, Attributes, and the dimid for unlimited dimension (probably time)
    !; nf90_inquire input: ncid; nf90_inquire output: ndimension, nVars, nGlobalAtts,unlimdimid
    CALL check_netcdf (nf90_inquire(ncid,ndimension, nVars, nGlobalAtts,unlimdimid))
    ALLOCATE (var_dimids(ndimension), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array var_dimids',__FILE__,__LINE__)
    var_dimids = 0

    !; the varid in netcdf files is counted from 1 to nVars
    !; look for variable names, number of dimension, var_dimids etc with nf90_inquire_variable
    !; nf90_inquire_variable input: ncid, varid; nf90_inquire_variable output name xtype, ndim, var_dimids, nAtts
    variables: DO varid=1,nVars
      CALL check_netcdf(nf90_inquire_variable(ncid,varid,varname,xtype, ndim, var_dimids, nAtts))

      IF (nAtts.gt.0) THEN
        ! get Attribute name
        DO attnum=1,nAtts
          ! input nf90_inq_attname: ncid, varid, attnum, output nf90_inq_attname name
          CALL check_netcdf ( nf90_inq_attname(ncid, varid, attnum, attname ))
          ! input nf90_inquire_attribute: ncid, varid, name; nf90_inquire_attribute output: xtype, length, attnum
          CALL check_netcdf(nf90_inquire_attribute(ncid, varid, attname, xtype, length))
          ! nf90_get_att input: ncid, varid, name; nf90_get_att ou put: attributetext
          ! note, attributetext should be the right type of variable
          getattribute: SELECT CASE (xtype)
            CASE (NF90_CHAR) ! for character attributes
              CALL check_netcdf(nf90_get_att(ncid, varid, attname, attributetext))
            END SELECT getattribute
        END DO ! done with attributes
      ENDIF

      IF (trim(varname) == 'lon') THEN      ! here I know that the variable with longitude coordinates is called 'lon'
        CALL check_netcdf( nf90_inquire_dimension(ncid,var_dimids(1), dimname, len=length) )
        IF (length /= nlon_ndvi)  CALL logging%error('nlon_ndvi is not equal data dimension',__FILE__,__LINE__)
        CALL check_netcdf(nf90_get_var(ncid,varid,lon,start=(/1/),count=(/length/))) ! read from netcdf file into lon(:)
      ENDIF

      IF (trim(varname) == 'lat') THEN       ! here I know that the variable with latitude coordinates is called 'lat'
        CALL check_netcdf( nf90_inquire_dimension(ncid,var_dimids(1), dimname, len=length) )
        IF (length /= nlat_ndvi)  CALL logging%error('nlat_ndvi is not equal data dimension',__FILE__,__LINE__)
        CALL check_netcdf(nf90_get_var(ncid,varid,lat,start=(/1/),count=(/length/))) ! read from netcdf file into lat(:)
      ENDIF
    ENDDO variables

    startlon_ndvi = lon(1) ! longitude of the upper left grid element
    startlat_ndvi = lat(1) ! latitude of the upper left grid element

    dlon_ndvi = 360._wp/real(nlon_ndvi) ! dlon_ndvi in degrees
    dlat_ndvi = -180._wp/real(nlat_ndvi) ! dlat_ndvi in degrees

    !! recalculate the values for raw data longitude and latitude, since values read from netcdf file are inaccurate
    DO i=1,nlon_ndvi-1
      lon(i+1) = lon(1) + REAL(i,wp)*dlon_ndvi
    ENDDO
    DO i=1,nlat_ndvi-1
      lat(i+1) = lat(1) + REAL(i,wp)*dlat_ndvi
    ENDDO

  END SUBROUTINE get_NDVI_data_coordinates


  !> get one row of NDVI raw data from netcdf file (along zonal direction)
  SUBROUTINE get_one_row_NDVI_data(ncid,                &
       &                           nlon_ndvi,           & 
       &                           nlat_ndvi,           &
       &                           ntime_ndvi,          &
       &                           row_index,           &
       &                           time_index,          &
       &                           ndvi_raw_data_lonrow)

    INTEGER(KIND=i4), INTENT(IN)  :: ncid, &                             !< netcdf unit file number
         &                           nlon_ndvi, &       !< number of grid elements in zonal direction for NDVI data
         &                           nlat_ndvi, &       !< number of grid elements in meridional direction for NDVI data
         &                           ntime_ndvi, &      !< number of dates with NDVI data
         &                           row_index, &       !< the index of the data row to read in
         &                           time_index      !< the index of the time (month) to read in
                                  
    REAL (KIND=wp), INTENT(out)   :: ndvi_raw_data_lonrow(1:nlon_ndvi)  !< longitude of ndvi raw data in geographical system


    !local variables
    INTEGER(KIND=i4)              :: ndimension, &                       !< number of dimensions in netcdf file
         &                           nVars, &                            !< number of variables in netcdf file
         &                           nGlobalAtts, &                      !< number of gloabal Attributes in netcdf file
         &                           unlimdimid, &                       !< id of unlimited dimension (e.g. time) in netcdf file
         &                           length, &                           !< length of dimension
         &                           dim_lon, &                           !< length of dimension lon
         &                           dim_lat, &                           !< length of dimension lat
         &                           dim_time, &                           !< length of dimension time
         &                           varid, &                            !< id of variable
         &                           xtype, &                            !< netcdf type of variable/attribute
         &                           ndim, &                             !< number of dimensions of variable
         &                           nAtts, &                            !< number of attributes for a netcdf variable
         &                           errorcode, &                        !< error status variable
         &                           attnum                           !< counter for attribute number

    INTEGER(KIND=i4), ALLOCATABLE :: var_dimids(:)       !< id of variable dimensions, vector, maximal dimension ndimension

    CHARACTER (len=80)            :: dimname, &               !< name of dimension
         &                           varname, &               !< name of variable
         &                           attname, &               !< name of attribute
         &                           attributetext         !< attribute text


    !! look for numbers of dimensions, Variable, Attributes, and the dimid for the unlimited dimension (probably time)
    !; nf90_inquire input: ncid; nf90_inquire output: ndimension, nVars, nGlobalAtts,unlimdimid
    CALL check_netcdf (nf90_inquire(ncid,ndimension, nVars, nGlobalAtts,unlimdimid))
    ALLOCATE (var_dimids(ndimension), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array var_dimids',__FILE__,__LINE__)
    var_dimids = 0

    !; the varid in netcdf files is counted from 1 to nVars
    !; look for variable names, number of dimension, var_dimids etc with nf90_inquire_variable
    !; nf90_inquire_variable input: ncid, varid; nf90_inquire_variable output name xtype, ndim, var_dimids, nAtts
    variables: DO varid=1,nVars
      CALL check_netcdf(nf90_inquire_variable(ncid,varid,varname,xtype, ndim, var_dimids, nAtts))

      IF (nAtts.GT.0) THEN
        ! get Attribute name
        DO attnum=1,nAtts
          ! input nf90_inq_attname: ncid, varid, attnum, output nf90_inq_attname name
          CALL check_netcdf ( nf90_inq_attname(ncid, varid, attnum, attname ))
          ! input nf90_inquire_attribute: ncid, varid, name; nf90_inquire_attribute output: xtype, length, attnum
          CALL check_netcdf(nf90_inquire_attribute(ncid, varid, attname, xtype, length))
          ! nf90_get_att input: ncid, varid, name; nf90_get_att ou put: attributetext
          ! note, attributetext should be the right type of variable
          getattribute: SELECT CASE (xtype)
            CASE (NF90_CHAR) ! for character attributes
              CALL check_netcdf(nf90_get_att(ncid, varid, attname, attributetext))
          END SELECT getattribute
        END DO ! done with attributes
      ENDIF

      IF (trim(varname) == 'ndvi') THEN     ! here I know that the variable with latitude coordinates is called 'ndvi'

        CALL check_netcdf( nf90_inquire_dimension(ncid,var_dimids(1), dimname, len=dim_lon) )
        IF (dim_lon /= nlon_ndvi)   CALL logging%error('nlon_ndvi is not equal data dimension',__FILE__,__LINE__)
        CALL check_netcdf( nf90_inquire_dimension(ncid,var_dimids(2), dimname, len=dim_lat) )
        IF (dim_lat /= nlat_ndvi)   CALL logging%error('nlat_ndvi is not equal data dimension',__FILE__,__LINE__)

        CALL check_netcdf( nf90_inquire_dimension(ncid,var_dimids(3), dimname, len=dim_time) )
        IF (dim_time /= ntime_ndvi)  CALL logging%error('nlon_ndvi is not equal data dimension',__FILE__,__LINE__)

        CALL check_netcdf(nf90_get_var(ncid,varid,ndvi_raw_data_lonrow,                          &
             start=(/1,row_index,time_index/),count=(/nlon_ndvi,1,1/))) ! read from netcdf file into ndvi_raw_data_lonrow(:)
      ENDIF
    ENDDO variables

  END SUBROUTINE get_one_row_NDVI_data

  !> get a block of NDVI raw data from netcdf file (given startrow, endrow, startcolumn, endcolumn)
  SUBROUTINE get_block_NDVI_data(ncid,                &
       &                         nlon_ndvi,           & 
       &                         nlat_ndvi,           &
       &                         ntime_ndvi,          &
       &                         startcolumn_index,   &
       &                         endcolumn_index,     &
       &                         startrow_index,      &
       &                         endrow_index,        &
       &                         ncolumns,            &
       &                         nrows,               &
       &                         time_index,          &
       &                         ndvi_data_block)

    INTEGER (KIND=i4), INTENT(IN) :: nlon_ndvi, &       !< number of grid elements in zonal direction for NDVI data
         &                           nlat_ndvi, &       !< number of grid elements in meridional direction for NDVI data
         &                           ntime_ndvi, &      !< number of dates with NDVI data
         &                           startcolumn_index, &    !< the index of the startcolumn of data to read in
         &                           endcolumn_index, &      !< the index of the endcolumn of data to read in
         &                           startrow_index, &       !< the index of the startrow of data to read in
         &                           endrow_index, &         !< the index of the endrow of data to read in
         &                           ncolumns, &           !< number of columns of data block
         &                           nrows, &              !< number of rows of data block
         &                           ncid, &                             !< netcdf unit file number
         &                           time_index      !< the index of the time (month) to read in

    REAL (KIND=wp), INTENT(OUT)   :: ndvi_data_block(1:ncolumns,1:nrows) !< longitude of ndvi raw data in geographical system


    !local variables
    INTEGER(KIND=i4)              :: ndimension, &                       !< number of dimensions in netcdf file
         &                           nVars, &                            !< number of variables in netcdf file
         &                           nGlobalAtts, &                      !< number of gloabal Attributes in netcdf file
         &                           unlimdimid, &                       !< id of unlimited dimension (e.g. time) in netcdf file
         &                           length, &                           !< length of dimension
         &                           dim_lon, &                           !< length of dimension lon
         &                           dim_lat, &                           !< length of dimension lat
         &                           dim_time, &                           !< length of dimension time
         &                           varid, &                            !< id of variable
         &                           xtype, &                            !< netcdf type of variable/attribute
         &                           ndim, &                             !< number of dimensions of variable
         &                           nAtts, &                            !< number of attributes for a netcdf variable
         &                           errorcode, &                        !< error status variable
         &                           attnum                           !< counter for attribute number

    INTEGER(KIND=i4), ALLOCATABLE :: var_dimids(:)       !< id of variable dimensions, vector, maximal dimension ndimension
                                  
    CHARACTER (LEN=80)            :: dimname, &               !< name of dimension
         &                           varname, &               !< name of variable
         &                           attname, &               !< name of attribute
         &                           attributetext         !< attribute text


    IF ( (startcolumn_index < 1) .or. (startcolumn_index > nlon_ndvi)) then
      CALL logging%error('startcolumn_index out of range',__FILE__,__LINE__)
    ENDIF
    IF ( (endcolumn_index < 1) .or. (endcolumn_index > nlon_ndvi)) then
      CALL logging%error('endcolumn_index out of range',__FILE__,__LINE__)
    ENDIF
    IF ( (startrow_index < 1) .or. (startrow_index > nlat_ndvi)) then
      CALL logging%error('startrow_index out of range',__FILE__,__LINE__)
    ENDIF

    IF ( (endrow_index < 1) .or. (endrow_index > nlat_ndvi)) then
      CALL logging%error('endrow_index out of range',__FILE__,__LINE__)
    ENDIF

    !! look for numbers of dimensions, Variable, Attributes, and the dimid for the unlimited dimension (probably time)
    !; nf90_inquire input: ncid; nf90_inquire output: ndimension, nVars, nGlobalAtts,unlimdimid
    CALL check_netcdf (nf90_inquire(ncid,ndimension, nVars, nGlobalAtts,unlimdimid))
    ALLOCATE (var_dimids(ndimension), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array var_dimids',__FILE__,__LINE__)
    var_dimids = 0

    !; the varid in netcdf files is counted from 1 to nVars
    !; look for variable names, number of dimension, var_dimids etc with nf90_inquire_variable
    !; nf90_inquire_variable input: ncid, varid; nf90_inquire_variable output name xtype, ndim, var_dimids, nAtts
    variables: DO varid=1,nVars
      CALL check_netcdf(nf90_inquire_variable(ncid,varid,varname,xtype, ndim, var_dimids, nAtts))
      IF (nAtts.gt.0) THEN
        ! get Attribute name
        DO attnum=1,nAtts
          ! input nf90_inq_attname: ncid, varid, attnum, output nf90_inq_attname name
          CALL check_netcdf ( nf90_inq_attname(ncid, varid, attnum, attname ))
          ! input nf90_inquire_attribute: ncid, varid, name; nf90_inquire_attribute output: xtype, length, attnum
          CALL check_netcdf(nf90_inquire_attribute(ncid, varid, attname, xtype, length))
          ! nf90_get_att input: ncid, varid, name; nf90_get_att ou put: attributetext
          ! note, attributetext should be the right type of variable
          getattribute: SELECT CASE (xtype)
            CASE (NF90_CHAR) ! for character attributes
              CALL check_netcdf(nf90_get_att(ncid, varid, attname, attributetext))
          END SELECT getattribute
        END DO ! done with attributes
      ENDIF

      IF (TRIM(varname) == 'ndvi') THEN    ! here I know that the variable with latitude coordinates is called 'ndvi'
        CALL check_netcdf( nf90_inquire_dimension(ncid,var_dimids(1), dimname, len=dim_lon) )
        IF (dim_lon /= nlon_ndvi)   CALL logging%error('nlon_ndvi is not equal data dimension',__FILE__,__LINE__)

        CALL check_netcdf( nf90_inquire_dimension(ncid,var_dimids(2), dimname, len=dim_lat) )
        IF (dim_lat /= nlat_ndvi)   CALL logging%error('nlat_ndvi is not equal data dimension',__FILE__,__LINE__)

        CALL check_netcdf( nf90_inquire_dimension(ncid,var_dimids(3), dimname, len=dim_time) )
        IF (dim_time /= ntime_ndvi)  CALL logging%error('nlon_ndvi is not equal data dimension',__FILE__,__LINE__)

        CALL check_netcdf(nf90_get_var(ncid,varid,ndvi_data_block,                          &
             start=(/startcolumn_index,startrow_index,time_index/),count=(/ncolumns,nrows,1/))) 
        !< read from netcdf file into ndvi_raw_data_lonrow(:)
      ENDIF
    ENDDO variables

  END SUBROUTINE get_block_NDVI_data

  !> get a pixel of NDVI raw data from netcdf file (given grid index)
  SUBROUTINE get_pixel_NDVI_data(ncid,                &
       &                         nlon_ndvi,           & 
       &                         nlat_ndvi,           &
       &                         ntime_ndvi,          &
       &                         column_index,        &
       &                         row_index,           &
       &                         time_index,          &
       &                         ndvi_pixel_data)

    INTEGER (KIND=i4), INTENT(in) :: nlon_ndvi, &       !< number of grid elements in zonal direction for NDVI data
         &                           nlat_ndvi, &       !< number of grid elements in meridional direction for NDVI data
         &                           ntime_ndvi, &      !< number of dates with NDVI data
         &                           column_index, &    !< the index of the column of data to read in
         &                           row_index, &       !< the index of the trow of data to read in
         &                           time_index, &      !< the index of the time (month) to read in
         &                           ncid                             !< netcdf unit file number

    REAL (KIND=wp), INTENT(OUT)   :: ndvi_pixel_data      !< value of ndvi raw data pixel

    !local variables
    INTEGER(KIND=i4)              :: ndimension, &                       !< number of dimensions in netcdf file
         &                           nVars, &                            !< number of variables in netcdf file
         &                           nGlobalAtts, &                      !< number of gloabal Attributes in netcdf file
         &                           unlimdimid, &                       !< id of unlimited dimension (e.g. time) in netcdf file
         &                           dim_lon, &                          !< length of dimension lon
         &                           dim_lat, &                          !< length of dimension lat
         &                           dim_time, &                         !< length of dimension time
         &                           varid, &                            !< id of variable
         &                           xtype, &                            !< netcdf type of variable/attribute
         &                           ndim, &                             !< number of dimensions of variable
         &                           nAtts, &                            !< number of attributes for a netcdf variable
         &                           errorcode                        !< error status variable

    CHARACTER (len=80)           :: dimname, &               !< name of dimension
         &                          varname               !< name of variable

    INTEGER(KIND=i4),ALLOCATABLE :: var_dimids(:)       !< id of variable dimensions, vector, maximal dimension ndimension

    IF ( (column_index < 1) .or. (column_index > nlon_ndvi)) then
      CALL logging%error('column_index out of range',__FILE__,__LINE__)
    ENDIF
    IF ( (row_index < 1) .or. (row_index > nlat_ndvi)) then
      CALL logging%error('row_index out of range',__FILE__,__LINE__)
    ENDIF

    !! look for numbers of dimensions, Variable, Attributes, and the dimid for the unlimited dimension (probably time)
    !; nf90_inquire input: ncid; nf90_inquire output: ndimension, nVars, nGlobalAtts,unlimdimid
    CALL check_netcdf (nf90_inquire(ncid,ndimension, nVars, nGlobalAtts,unlimdimid))
    ALLOCATE (var_dimids(ndimension), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array var_dimids',__FILE__,__LINE__)
    var_dimids = 0

    !; the varid in netcdf files is counted from 1 to nVars
    !; look for variable names, number of dimension, var_dimids etc with nf90_inquire_variable
    !; nf90_inquire_variable input: ncid, varid; nf90_inquire_variable output name xtype, ndim, var_dimids, nAtts
    variables: DO varid=1,nVars
      CALL check_netcdf(nf90_inquire_variable(ncid,varid,varname,xtype, ndim, var_dimids, nAtts))
      IF (trim(varname) == 'ndvi') THEN   ! here I know that the variable with latitude coordinates is called 'ndvi'
        CALL check_netcdf( nf90_inquire_dimension(ncid,var_dimids(1), dimname, len=dim_lon) )
        IF (dim_lon /= nlon_ndvi)   CALL logging%error('nlon_ndvi is not equal data dimension',__FILE__,__LINE__)

        CALL check_netcdf( nf90_inquire_dimension(ncid,var_dimids(2), dimname, len=dim_lat) )
        IF (dim_lat /= nlat_ndvi)   CALL logging%error('nlat_ndvi is not equal data dimension',__FILE__,__LINE__)

        CALL check_netcdf( nf90_inquire_dimension(ncid,var_dimids(3), dimname, len=dim_time) )
        IF (dim_time /= ntime_ndvi)  CALL logging%error('nlon_ndvi is not equal data dimension',__FILE__,__LINE__)

        CALL check_netcdf(nf90_get_var(ncid,varid,ndvi_pixel_data,                      &
             start=(/column_index,row_index,time_index/))) ! read from netcdf file into ndvi_pixel_data
      ENDIF
    ENDDO variables

  END SUBROUTINE get_pixel_NDVI_data

END MODULE mo_ndvi_routines
