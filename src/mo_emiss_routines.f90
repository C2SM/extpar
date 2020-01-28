!+ Fortran module with EMISS data handling routines
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
!> Fortran module with EMISS data handling routines
!> \author Hermann Asensio
!>
MODULE mo_emiss_routines

  USE mo_logging
  USE mo_kind,                  ONLY: wp, i4

  USE netcdf,                   ONLY:   &
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

  USE mo_utilities_extpar,      ONLY: free_un ! function to get free unit number

  USE mo_io_utilities,          ONLY: check_netcdf

  USE mo_io_units,              ONLY: filename_max


  IMPLICIT NONE

  PRIVATE
   
  PUBLIC :: open_netcdf_EMISS_data, &
       &    close_netcdf_EMISS_data, &
       &    read_namelists_extpar_emiss, &
       &    read_emiss_data_input_namelist, &
       &    get_dimension_EMISS_data, &
       &    get_EMISS_data_coordinates, &
       &    get_one_row_EMISS_data, &
       &    get_block_EMISS_data, &
       &    get_pixel_EMISS_data


  CONTAINS

  !---------------------------------------------------------------------------
  !> subroutine to read namelist for EMISS data settings for EXTPAR 
  SUBROUTINE read_namelists_extpar_emiss(namelist_file, &
       raw_data_emiss_path, &
       raw_data_emiss_filename, &
       emiss_buffer_file, &
       emiss_output_file)

    CHARACTER (len=filename_max), INTENT(IN) :: namelist_file !< filename with namelists for for EXTPAR settings

    CHARACTER (len=filename_max)             :: raw_data_emiss_path, &         !< path to raw data
         &                                      raw_data_emiss_filename, &  !< filename EMISS raw data
         &                                      emiss_buffer_file, &  !< name for EMISS buffer file
         &                                      emiss_output_file !< name for EMISS output file

    INTEGER(KIND=i4)                         :: nuin, ierr

    !> namelist with filenames for EMISS data input
    NAMELIST /emiss_raw_data/ raw_data_emiss_path, raw_data_emiss_filename
    !> namelist with filenames for EMISS data output
    NAMELIST /emiss_io_extpar/ emiss_buffer_file, emiss_output_file

    nuin = free_un()  ! functioin free_un returns free Fortran unit number
    OPEN(nuin,FILE=TRIM(namelist_file), IOSTAT=ierr)
    IF (ierr /= 0) THEN
      WRITE(message_text,*)'Cannot open ', TRIM(namelist_file)
      CALL logging%error(message_text,__FILE__, __LINE__) 
    ENDIF

    READ(nuin, NML=emiss_raw_data, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      CALL logging%error('Cannot read in namelist emiss_raw_data',__FILE__, __LINE__) 
    ENDIF

    READ(nuin, NML=emiss_io_extpar, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      CALL logging%error('Cannot read in namelist emiss_io_extpar',__FILE__, __LINE__) 
    ENDIF

    CLOSE(nuin)

  END SUBROUTINE read_namelists_extpar_emiss
  !---------------------------------------------------------------------------

  !> open netcdf-file and get netcdf unit file number
  SUBROUTINE open_netcdf_EMISS_data(path_emiss_file, ncid)
            
    CHARACTER (len=*), INTENT(in) :: path_emiss_file         !< filename with path for EMISS raw data
    INTEGER, INTENT(out)          :: ncid                             !< netcdf unit file number

    CALL check_netcdf( nf90_open(TRIM(path_emiss_file),NF90_NOWRITE, ncid))

  END SUBROUTINE open_netcdf_EMISS_data

  !> close netcdf-file 
  SUBROUTINE close_netcdf_EMISS_data(ncid)
    INTEGER, INTENT(in) :: ncid                             !< netcdf unit file number

    !! close netcdf file 
    call check_netcdf( nf90_close( ncid))

  END SUBROUTINE close_netcdf_EMISS_data


  !> read namelist with settings for COSMO target grid
  !> \author Hermann Asensio
  SUBROUTINE read_emiss_data_input_namelist(input_namelist_file,        &
       &                                    raw_data_path,             &
       &                                    raw_data_emiss_filename,    &
       &                                    outputgrid_emiss_filename)

    CHARACTER (len=*), INTENT(in)             :: input_namelist_file !< file with input namelist 
    CHARACTER (len=filename_max), INTENT(out) :: raw_data_path, &         !< path to raw data
         &                                       raw_data_emiss_filename, &  !< filename emiss raw data
         &                                       outputgrid_emiss_filename !< output filename

    INTEGER (KIND=i4) :: ierr, nuin

    !>Define the namelist group
    NAMELIST /emiss_data_input/ raw_data_path, raw_data_emiss_filename, outputgrid_emiss_filename

    nuin = free_un()  ! functioin free_un returns free Fortran unit number
    OPEN(nuin,FILE=input_namelist_file, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      WRITE(message_text,*)'Cannot open ', TRIM(input_namelist_file)
      CALL logging%error(message_text,__FILE__, __LINE__) 
    ENDIF

    READ(nuin, NML=emiss_data_input, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      CALL logging%error('Cannot read in namelist emmis_data_input',__FILE__, __LINE__) 
    ENDIF

    CLOSE(nuin)

  END SUBROUTINE read_emiss_data_input_namelist


  !> inquire dimension information for EMISS raw data 
  SUBROUTINE get_dimension_EMISS_data(ncid, &
       &                              nlon_emiss, &
       &                              nlat_emiss, &
       &                              ntime_emiss)

    INTEGER(KIND=i4), INTENT(in)   :: ncid                             !< netcdf unit file number

    INTEGER (KIND=i4), INTENT(out) :: nlon_emiss, &  !< number of grid elements in zonal direction for EMISS data
         &                            nlat_emiss, &  !< number of grid elements in meridional direction for EMISS data
         &                            ntime_emiss !< number of dates with EMISS data

    !local variables
    INTEGER(KIND=i4)               :: ndimension, &                        !< number of dimensions in netcdf file
         &                            nVars, &                             !< number of variables in netcdf file
         &                            nGlobalAtts, &                       !< number of gloabal Attributes in netcdf file
         &                            unlimdimid, &                        !< id of unlimited dimension (e.g. time) in netcdf file
         &                            dimid, &                             !< id of dimension
         &                            length                           !< length of dimension

    CHARACTER (len=80)             :: dimname               !< name of dimensiona



    CALL check_netcdf (nf90_inquire(ncid,ndimension, nVars, nGlobalAtts,unlimdimid))

    !; the dimid in netcdf-files is counted from 1 to ndimension
    !; look for the name and length of the dimension with f90_inquire_dimension
    !; nf90_inquire_dimension input: ncid, dimid; nf90_inquire_dimension output: name, length
    DO dimid=1,ndimension
      CALL check_netcdf( nf90_inquire_dimension(ncid,dimid, dimname, length) )
      IF ( TRIM(dimname) == 'lon') nlon_emiss=length          ! here I know that the name of zonal dimension is 'lon'
      IF ( TRIM(dimname) == 'lat') nlat_emiss=length          ! here I know that the name of meridional dimension is 'lat'
      IF ( TRIM(dimname) == 'time') ntime_emiss=length        ! here I know that the name of time dimension is 'time'
    ENDDO

  END SUBROUTINE get_dimension_EMISS_data

  !> read coordinates for EMISS raw data from netcdf file
  SUBROUTINE get_EMISS_data_coordinates(ncid,          &
       &                                nlon_emiss,     &
       &                                nlat_emiss,     &  
       &                                startlon_emiss, &
       &                                startlat_emiss, &
       &                                dlon_emiss,     &
       &                                dlat_emiss,     &
       &                                lon,           &
       &                                lat)

    INTEGER(KIND=i4), INTENT(in):: ncid, &                              !< netcdf unit file number
                                   nlon_emiss, &  !< number of grid elements in zonal direction for EMISS data
                                   nlat_emiss !< number of grid elements in meridional direction for EMISS data

    REAL (KIND=wp), INTENT(out) :: startlon_emiss, &  !< longitude of lower left grid element for EMISS data 
         &                         startlat_emiss, &  !< latitude of lower left grid element for EMISS data
         &                         dlon_emiss, &  !< grid point distance in zonal direction (in degrees) for EMISS data
         &                         dlat_emiss, &  !< grid point distance in meridional direction (in degrees) for EMISS data
         &                         lon(1:nlon_emiss), &       !< longitude of emiss raw data in geographical system
         &                         lat(1:nlat_emiss)      !< latitude of emiss raw date in geographical system

    !local variables
    INTEGER(KIND=i4)            :: ndimension, &                        !< number of dimensions in netcdf file
         &                         nVars, &                             !< number of variables in netcdf file
         &                         nGlobalAtts, &                       !< number of gloabal Attributes in netcdf file
         &                         unlimdimid, &                        !< id of unlimited dimension (e.g. time) in netcdf file
         &                         length, &                            !< length of dimension
         &                         varid, &                             !< id of variable
         &                         xtype, &                             !< netcdf type of variable/attribute
         &                         ndim, &                              !< number of dimensions of variable
         &                         errorcode, &                         !< error status variable
         &                         attnum, &                            !< counter for attribute number
         &                         nAtts, &                             !< number of attributes for a netcdf variable
         &                         i                                !< loop index

    INTEGER, ALLOCATABLE        :: var_dimids(:)       !< id of variable dimensions, vector, maximal dimension ndimension

    CHARACTER (len=80)          :: dimname, &                !< name of dimension
         &                         varname, &                !< name of variable
         &                         attname, &                !< name of attribute
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
        IF (length /= nlon_emiss)  CALL logging%error('nlon_emiss is not equal data dimension',__FILE__,__LINE__)
        CALL check_netcdf(nf90_get_var(ncid,varid,lon,start=(/1/),count=(/length/))) ! read from netcdf file into lon(:)
      ENDIF

      IF (trim(varname) == 'lat') THEN       ! here I know that the variable with latitude coordinates is called 'lat'
        CALL check_netcdf( nf90_inquire_dimension(ncid,var_dimids(1), dimname, len=length) )
        IF (length /= nlat_emiss)  CALL logging%error('nlat_emiss is not equal data dimension',__FILE__,__LINE__)
        CALL check_netcdf(nf90_get_var(ncid,varid,lat,start=(/1/),count=(/length/))) ! read from netcdf file into lat(:)
      ENDIF

    ENDDO variables

    startlon_emiss = lon(1) ! longitude of the upper left grid element
    startlat_emiss = lat(1) ! latitude of the upper left grid element

    dlon_emiss = 360._wp/real(nlon_emiss) ! dlon_emiss in degrees
    dlat_emiss = -180._wp/real(nlat_emiss) ! dlat_emiss in degrees

    !! recalculate the values for raw data longitude and latitude, since values read from netcdf file are inaccurate
    DO i=1,nlon_emiss-1
      lon(i+1) = lon(1) + REAL(i,wp)*dlon_emiss
    ENDDO
    DO i=1,nlat_emiss-1
      lat(i+1) = lat(1) + REAL(i,wp)*dlat_emiss
    ENDDO

  END SUBROUTINE get_EMISS_data_coordinates

  !> get one row of EMISS raw data from netcdf file (along zonal direction)
  SUBROUTINE get_one_row_EMISS_data(ncid,                &
       &                            nlon_emiss,           & 
       &                            nlat_emiss,           &
       &                            ntime_emiss,          &
       &                            row_index,           &
       &                            time_index,          &
       &                            emiss_raw_data_lonrow)

    INTEGER (KIND=i4), INTENT(in) :: nlon_emiss, &        !< number of grid elements in zonal direction for EMISS data
         &                           nlat_emiss, &        !< number of grid elements in meridional direction for EMISS data
         &                           ntime_emiss, &       !< number of dates with EMISS data
         &                           row_index, &        !< the index of the data row to read in
         &                           ncid, &
         &                           time_index      !< the index of the time (month) to read in

    REAL (KIND=wp), INTENT(out)   :: emiss_raw_data_lonrow(1:nlon_emiss)  !< longitude of emiss raw data in geographical system


    !local variables
    INTEGER(KIND=i4)              :: ndimension, &                        !< number of dimensions in netcdf file
         &                           nVars, &                             !< number of variables in netcdf file
         &                           nGlobalAtts, &                       !< number of gloabal Attributes in netcdf file
         &                           unlimdimid, &                        !< id of unlimited dimension (e.g. time) in netcdf file
         &                           length, &                            !< length of dimension
         &                           dim_lon, &                            !< length of dimension lon
         &                           dim_lat, &                            !< length of dimension lat
         &                           dim_time, &                            !< length of dimension time
         &                           varid, &                             !< id of variable
         &                           xtype, &                             !< netcdf type of variable/attribute
         &                           ndim, &                              !< number of dimensions of variable
         &                           nAtts, &                             !< number of attributes for a netcdf variable
         &                           errorcode, &                         !< error status variable
         &                           attnum                           !< counter for attribute number

    INTEGER, ALLOCATABLE          :: var_dimids(:)       !< id of variable dimensions, vector, maximal dimension ndimension

    CHARACTER (len=80)            :: varname, &                !< name of variable
         &                           dimname, &                !< name of dimension
         &                           attname, &                !< name of attribute
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

      IF (trim(varname) == 'bbemis_longwave') THEN     ! here I know that the variable with latitude coordinates is called 'emiss'

        CALL check_netcdf( nf90_inquire_dimension(ncid,var_dimids(1), dimname, len=dim_lon) )
        IF (dim_lon /= nlon_emiss)   CALL logging%error('nlon_emiss is not equal data dimension',__FILE__,__LINE__)
        CALL check_netcdf( nf90_inquire_dimension(ncid,var_dimids(2), dimname, len=dim_lat) )
        IF (dim_lat /= nlat_emiss)   CALL logging%error('nlat_emiss is not equal data dimension',__FILE__,__LINE__)

        CALL check_netcdf( nf90_inquire_dimension(ncid,var_dimids(3), dimname, len=dim_time) )
        IF (dim_time /= ntime_emiss)  CALL logging%error('nlon_emiss is not equal data dimension',__FILE__,__LINE__)

        CALL check_netcdf(nf90_get_var(ncid,varid,emiss_raw_data_lonrow,                          &
             start=(/1,row_index,time_index/),count=(/nlon_emiss,1,1/))) ! read from netcdf file into emiss_raw_data_lonrow(:)
      ENDIF
    ENDDO variables

  END SUBROUTINE get_one_row_EMISS_data

  !> get a block of EMISS raw data from netcdf file (given startrow, endrow, startcolumn, endcolumn)
  SUBROUTINE get_block_EMISS_data(ncid,                &
       nlon_emiss,           & 
       nlat_emiss,           &
       ntime_emiss,          &
       startcolumn_index,   &
       endcolumn_index,     &
       startrow_index,      &
       endrow_index,        &
       ncolumns,            &
       nrows,               &
       time_index,          &
       emiss_data_block)

    INTEGER, INTENT(in) :: ncid                             !< netcdf unit file number
    INTEGER (KIND=i4), INTENT(in) :: nlon_emiss       !< number of grid elements in zonal direction for EMISS data
    INTEGER (KIND=i4), INTENT(in) :: nlat_emiss       !< number of grid elements in meridional direction for EMISS data
    INTEGER (KIND=i4), INTENT(in) :: ntime_emiss      !< number of dates with EMISS data
    INTEGER (KIND=i4), INTENT(in) :: startcolumn_index    !< the index of the startcolumn of data to read in
    INTEGER (KIND=i4), INTENT(in) :: endcolumn_index      !< the index of the endcolumn of data to read in
    INTEGER (KIND=i4), INTENT(in) :: startrow_index       !< the index of the startrow of data to read in
    INTEGER (KIND=i4), INTENT(in) :: endrow_index         !< the index of the endrow of data to read in
    INTEGER (KIND=i4), INTENT(in) ::   ncolumns           !< number of columns of data block
    INTEGER (KIND=i4), INTENT(in) ::   nrows              !< number of rows of data block

    INTEGER (KIND=i4), INTENT(in) :: time_index      !< the index of the time (month) to read in

    REAL (KIND=wp), INTENT(out) :: emiss_data_block(1:ncolumns,1:nrows) !< longitude of emiss raw data in geographical system


    !local variables
    INTEGER :: ndimension                       !< number of dimensions in netcdf file
    INTEGER :: nVars                            !< number of variables in netcdf file
    INTEGER :: nGlobalAtts                      !< number of gloabal Attributes in netcdf file
    INTEGER :: unlimdimid                       !< id of unlimited dimension (e.g. time) in netcdf file

    CHARACTER (LEN=80) :: dimname               !< name of dimension
    INTEGER :: length                           !< length of dimension
    INTEGER :: dim_lon                           !< length of dimension lon
    INTEGER :: dim_lat                           !< length of dimension lat
    INTEGER :: dim_time                           !< length of dimension time



    INTEGER :: varid                            !< id of variable
    CHARACTER (LEN=80) :: varname               !< name of variable
    INTEGER :: xtype                            !< netcdf type of variable/attribute
    INTEGER :: ndim                             !< number of dimensions of variable
    INTEGER, ALLOCATABLE :: var_dimids(:)       !< id of variable dimensions, vector, maximal dimension ndimension
    INTEGER :: nAtts                            !< number of attributes for a netcdf variable

    INTEGER :: errorcode                        !< error status variable

    INTEGER :: attnum                           !< counter for attribute number
    CHARACTER (LEN=80) :: attname               !< name of attribute
    CHARACTER (LEN=80) :: attributetext         !< attribute text


    IF ( (startcolumn_index < 1) .or. (startcolumn_index > nlon_emiss)) then
      CALL logging%error('startcolumn_index out of range')
    ENDIF
    IF ( (endcolumn_index < 1) .or. (endcolumn_index > nlon_emiss)) then
      CALL logging%error('endcolumn_index out of range')
    ENDIF
    IF ( (startrow_index < 1) .or. (startrow_index > nlat_emiss)) then
      CALL logging%error('startrow_index out of range')
    ENDIF

    IF ( (endrow_index < 1) .or. (endrow_index > nlat_emiss)) then
      CALL logging%error('endrow_index out of range')
    ENDIF

    !ncolumns = endcolumn_index - startcolumn_index + 1
    !nrows    = endrow_index - startrow_index + 1

    !! look for numbers of dimensions, Variable, Attributes, and the dimid for the unlimited dimension (probably time)
    !; nf90_inquire input: ncid; nf90_inquire output: ndimension, nVars, nGlobalAtts,unlimdimid
    CALL check_netcdf (nf90_inquire(ncid,ndimension, nVars, nGlobalAtts,unlimdimid))
    !print *,'ncid,ndimension, nVars, nGlobalAtts,unlimdimid',ncid,ndimension, nVars, nGlobalAtts,unlimdimid
    ALLOCATE (var_dimids(ndimension), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array var_dimids')
    var_dimids = 0


    !; the varid in netcdf files is counted from 1 to nVars
    !; look for variable names, number of dimension, var_dimids etc with nf90_inquire_variable
    !; nf90_inquire_variable input: ncid, varid; nf90_inquire_variable output name xtype, ndim, var_dimids, nAtts
    !print *,'nVars ',nVars
    variables: DO varid=1,nVars
      !print *,'variable loop, varid ',varid
      CALL check_netcdf(nf90_inquire_variable(ncid,varid,varname,xtype, ndim, var_dimids, nAtts))
      !print *,'------------------------------------'
      !print *,'ncid,varid,vrname',ncid,varid,trim(varname)
      !print *, 'xtype',xtype
      !print *, 'ndim', ndim
      !print *, 'var_dimids', var_dimids
      !print *, 'nAtts', nAtts
      !print *,'------------------------------------'

      IF (nAtts.gt.0) THEN
        ! get Attribute name
        DO attnum=1,nAtts
          ! input nf90_inq_attname: ncid, varid, attnum, output nf90_inq_attname name
          CALL check_netcdf ( nf90_inq_attname(ncid, varid, attnum, attname ))
          ! print *,'ncid,varid,attnumm, attname',ncid,varid,attnum,trim(attname)
          ! input nf90_inquire_attribute: ncid, varid, name; nf90_inquire_attribute output: xtype, length, attnum
          CALL check_netcdf(nf90_inquire_attribute(ncid, varid, attname, xtype, length))
          !print *, 'ncid, varid, attname, xtype, length, attnum',ncid,varid,trim(attname),' ',xtype,length,attnum
          ! nf90_get_att input: ncid, varid, name; nf90_get_att ou put: attributetext
          ! note, attributetext should be the right type of variable
          getattribute: SELECT CASE (xtype)
          CASE (NF90_CHAR) ! for character attributes
            CALL check_netcdf(nf90_get_att(ncid, varid, attname, attributetext))
            ! print *,'get attribute: ncid, varid, attname, attributetext',ncid,varid,trim(attname),' ',trim(attributetext)
          END SELECT getattribute
        END DO ! done with attributes
      ENDIF


      IF (TRIM(varname) == 'bbemis_longwave') THEN    ! here I know that the variable with latitude coordinates is called 'emiss'

        !print *,'------------------------------------'
        !print *,'ncid,varid,vrname',ncid,varid,trim(varname)
        !print *, 'xtype',xtype
        !print *, 'ndim', ndim
        !print *, 'var_dimids', var_dimids
        !print *, 'nAtts', nAtts
        !print *,'------------------------------------'

        CALL check_netcdf( nf90_inquire_dimension(ncid,var_dimids(1), dimname, len=dim_lon) )
        ! print *, 'var_dimids(1), dimname, dim_lon: ', var_dimids(1), TRIM(dimname), dim_lon
        IF (dim_lon /= nlon_emiss)   CALL logging%error('nlon_emiss is not equal data dimension')

        CALL check_netcdf( nf90_inquire_dimension(ncid,var_dimids(2), dimname, len=dim_lat) )
        ! print *, 'var_dimids(2), dimname, dim_lat: ', var_dimids(2), TRIM(dimname), dim_lat
        IF (dim_lat /= nlat_emiss)   CALL logging%error('nlat_emiss is not equal data dimension')

        CALL check_netcdf( nf90_inquire_dimension(ncid,var_dimids(3), dimname, len=dim_time) )
        ! print *, 'var_dimids(3), dimname, dim_time: ', var_dimids(3), TRIM(dimname), dim_time
        IF (dim_time /= ntime_emiss)  CALL logging%error('nlon_emiss is not equal data dimension')


        !print *,'dimt ',dimt
        CALL check_netcdf(nf90_get_var(ncid,varid,emiss_data_block,                          &
             start=(/startcolumn_index,startrow_index,time_index/),count=(/ncolumns,nrows,1/))) 
        !< read from netcdf file into emiss_raw_data_lonrow(:)
      ENDIF

    ENDDO variables


  END SUBROUTINE get_block_EMISS_data

  !> get a pixel of EMISS raw data from netcdf file (given grid index)
  SUBROUTINE get_pixel_EMISS_data(ncid,                &
       &                          nlon_emiss,           & 
       &                          nlat_emiss,           &
       &                          ntime_emiss,          &
       &                          column_index,        &
       &                          row_index,           &
       &                          time_index,          &
       &                          emiss_pixel_data)

    INTEGER (KIND=i4), INTENT(in) :: nlon_emiss, &        !< number of grid elements in zonal direction for EMISS data
         &                           nlat_emiss, &        !< number of grid elements in meridional direction for EMISS data
         &                           ntime_emiss, &       !< number of dates with EMISS data
         &                           column_index, &     !< the index of the column of data to read in
         &                           row_index, &        !< the index of the trow of data to read in
         &                           time_index, &       !< the index of the time (month) to read in
         &                           ncid

    REAL (KIND=wp), INTENT(out)   :: emiss_pixel_data      !< value of emiss raw data pixel


    !local variables
    INTEGER(KIND=i4)              :: ndimension, &                        !< number of dimensions in netcdf file
         &                           nVars, &                             !< number of variables in netcdf file
         &                           nGlobalAtts, &                       !< number of gloabal Attributes in netcdf file
         &                           unlimdimid, &                        !< id of unlimited dimension (e.g. time) in netcdf file
         &                           dim_lon, &                           !< length of dimension lon
         &                           dim_lat, &                           !< length of dimension lat
         &                           dim_time, &                          !< length of dimension time
         &                           varid, &                             !< id of variable
         &                           xtype, &                             !< netcdf type of variable/attribute
         &                           ndim, &                              !< number of dimensions of variable
         &                           nAtts, &                             !< number of attributes for a netcdf variable
         &                           errorcode                        !< error status variable

    INTEGER(KIND=i4), ALLOCATABLE :: var_dimids(:)       !< id of variable dimensions, vector, maximal dimension ndimension

    CHARACTER (len=80)            :: dimname, &                !< name of dimension
         &                           varname               !< name of variable

    IF ( (column_index < 1) .or. (column_index > nlon_emiss)) THEN
      CALL logging%error('column_index out of range',__FILE__,__LINE__)
    ENDIF
    IF ( (row_index < 1) .or. (row_index > nlat_emiss)) THEN
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

      IF(trim(varname) == 'bbemis_longwavebbemis_longwave') THEN 
        CALL check_netcdf( nf90_inquire_dimension(ncid,var_dimids(1), dimname, len=dim_lon) )
        IF (dim_lon /= nlon_emiss)   CALL logging%error('nlon_emiss is not equal data dimension',__FILE__,__LINE__)

        CALL check_netcdf( nf90_inquire_dimension(ncid,var_dimids(2), dimname, len=dim_lat) )
        IF (dim_lat /= nlat_emiss)   CALL logging%error('nlat_emiss is not equal data dimension',__FILE__,__LINE__)

        CALL check_netcdf( nf90_inquire_dimension(ncid,var_dimids(3), dimname, len=dim_time) )
        IF (dim_time /= ntime_emiss)  CALL logging%error('nlon_emiss is not equal data dimension',__FILE__,__LINE__)

        CALL check_netcdf(nf90_get_var(ncid,varid,emiss_pixel_data,                      &
             start=(/column_index,row_index,time_index/))) ! read from netcdf file into emiss_pixel_data
      ENDIF
    ENDDO variables

  END SUBROUTINE get_pixel_EMISS_data

END MODULE mo_emiss_routines
