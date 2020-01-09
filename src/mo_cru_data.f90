!+ Fortran modules with data fields for CRU temperature data
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_0         2010/12/21 Hermann Asensio
!  Initial release
! V2_0         2013/06/04 Martina Messmer
!  introduction of a finer CRU temperature data set (CLM Community)
!  new parameter for the CRU temperature elevation
!  introduction of a deallocate subroutine
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!> Fortran modules with data fields for CRU temperature data
!> \author Hermann Asensio
!
MODULE mo_cru_data

  USE mo_kind,             ONLY: wp, i4, i8
  USE mo_logging,          
  USE mo_utilities_extpar, ONLY: abort_extpar
  USE mo_io_units,         ONLY: filename_max  
  USE mo_io_utilities,     ONLY: check_netcdf
  USE mo_grid_structures,  ONLY: reg_lonlat_grid
  
  USE netcdf,     ONLY :   &
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

  USE netcdf,     ONLY:     &
       nf90_create,             &
       nf90_def_dim,            &
       nf90_def_var,            &
       nf90_enddef,             &
       nf90_redef,              &
       nf90_put_att,            &
       nf90_put_var

  USE netcdf,     ONLY :   &
       NF90_CHAR,               &
       NF90_DOUBLE,             &
       NF90_FLOAT,              &
       NF90_INT,                &
       NF90_BYTE,               &
       NF90_SHORT
  
  USE netcdf,     ONLY :   &
       NF90_GLOBAL,             &
       NF90_UNLIMITED,          &
       NF90_CLOBBER,            &
       NF90_NOWRITE

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: cru_grid
  PUBLIC :: allocate_cru_data, &
       deallocate_cru_data, &
       read_cru_data_input_namelist, &
       read_namelists_extpar_t_clim, &
       get_dimension_cru_data, &
       get_cru_grid_and_data, &
       lon_cru, &
       lat_cru, &
       cru_raw_data, &
       cru_raw_elev

  TYPE(reg_lonlat_grid) :: cru_grid !< structure with defenition of the raw data grid for the AOT dataset

  REAL (wp), ALLOCATABLE :: lon_cru(:) !< longitude of aot grid
  REAL (wp), ALLOCATABLE :: lat_cru(:) !< latitude of aot grid

  REAL (wp), ALLOCATABLE :: cru_raw_data(:,:,:) !< aerosol optical thickness, aot(ncolumns,nrows,ntime) 
  REAL (wp), ALLOCATABLE :: cru_raw_elev(:,:,:) !< surface height in cru (ncolumns,nrows,ntime)

CONTAINS

  !---------------------------------------------------------------------------
  !> subroutine to read namelist for t_clim data settings for EXTPAR 
  SUBROUTINE read_namelists_extpar_t_clim(namelist_file,            &
       &                                  it_cl_type,               &
       &                                  raw_data_t_clim_path,     &
       &                                  raw_data_t_clim_filename, &
       &                                  t_clim_buffer_file,       &
       &                                  t_clim_output_file)

    CHARACTER (len=*), INTENT(IN)  :: namelist_file !< filename with namelists for for EXTPAR settings
    INTEGER (i8),      INTENT(OUT) :: it_cl_type    !< integer switch to choose a land use raw data set
    ! 1 CRU fine (new), 2 CRU coarse (old) temperature climatology
    CHARACTER (len=filename_max), INTENT(OUT) :: raw_data_t_clim_path        !< path to raw data
    CHARACTER (len=filename_max), INTENT(OUT) :: raw_data_t_clim_filename    !< filename temperature climatology raw data
    
    CHARACTER (len=filename_max), INTENT(OUT) :: t_clim_buffer_file !< name for temperature climatology buffer
    CHARACTER (len=filename_max), INTENT(OUT) :: t_clim_output_file !< name for temperature climatology output file
    
    !> namelist with filename for temperature climatlogy data output
    NAMELIST /t_clim_raw_data/ raw_data_t_clim_path, raw_data_t_clim_filename, it_cl_type
    
    !> namelist with filename for temperature climatlogy data output
    NAMELIST /t_clim_io_extpar/ t_clim_buffer_file, t_clim_output_file

    INTEGER :: nuin !< unit number
    INTEGER :: ierr !< error flag

    message_text = ''

    it_cl_type = -1
    
    raw_data_t_clim_path = ''
    raw_data_t_clim_filename = ''

    t_clim_buffer_file = ''
    t_clim_output_file = ''

    OPEN(NEWUNIT=nuin,FILE=TRIM(namelist_file), IOSTAT=ierr, IOMSG=message_text)
    IF (ierr /= 0) THEN
      CALL abort_extpar('CRU namelist open error: '//TRIM(message_text), __FILE__, __LINE__)
    ENDIF
    READ(nuin, NML=t_clim_raw_data, IOSTAT=ierr, IOMSG=message_text)
    IF (ierr /= 0) THEN
      WRITE(0,NML=t_clim_raw_data)
      CALL abort_extpar('CRU raw data namelist read error: '//TRIM(message_text), __FILE__, __LINE__)      
    ENDIF
    READ(nuin, NML=t_clim_io_extpar, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      WRITE(0,NML=t_clim_io_extpar)
      CALL abort_extpar('CRU io namelist read error: '//TRIM(message_text), __FILE__, __LINE__)      
    ENDIF
    CLOSE(nuin)
    
  END SUBROUTINE read_namelists_extpar_t_clim
  !---------------------------------------------------------------------------


  !> subroutine to allocate aot data fields
  SUBROUTINE allocate_cru_data(nrows,ncolumns,ntime)
    IMPLICIT NONE
    INTEGER (i8), INTENT(IN) :: nrows !< number of rows
    INTEGER (i8), INTENT(IN) :: ncolumns !< number of columns
    INTEGER (i8), INTENT(IN) :: ntime !< number of times

    INTEGER :: errorcode !< error status variable


    ALLOCATE (lon_cru(1:ncolumns), STAT=errorcode)
    IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array lon_cru')
    lon_cru = 0.0

    ALLOCATE (lat_cru(1:nrows), STAT=errorcode)
    IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array lat_cru')
    lat_cru = 0.0

    ALLOCATE (cru_raw_data(1:ncolumns,1:nrows,1:ntime),STAT=errorcode)
    IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array cru_raw_data')
    cru_raw_data = 0.0

    ALLOCATE (cru_raw_elev(1:ncolumns,1:nrows,1:ntime),STAT=errorcode)
    IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array cru_raw_elev')
    cru_raw_elev = 0.0

  END SUBROUTINE allocate_cru_data


  !> subroutine to read namelist
  !! \author Hermann Asensio
  SUBROUTINE read_cru_data_input_namelist(input_namelist_file, cru_filename)


    CHARACTER (LEN=*), INTENT(IN)  :: input_namelist_file !< file with input namelist 
    CHARACTER (LEN=filename_max), INTENT(OUT) :: cru_filename  !< filename aot raw data

    !>Define the namelist group
    NAMELIST /cru_file_info/ cru_filename

    INTEGER :: ierr !< error flag
    INTEGER :: nuin !< unit number

    open(NEWUNIT=nuin,FILE=TRIM(input_namelist_file), IOSTAT=ierr)
    read(nuin, NML=cru_file_info, IOSTAT=ierr)
    close(nuin)

  END SUBROUTINE read_cru_data_input_namelist

  !> get dimension information of aot data from netcdf file
  SUBROUTINE get_dimension_cru_data(cru_filename, &
       nrows,        &
       ncolumns,     &
       ntime)

    IMPLICIT NONE
    CHARACTER (LEN=*), INTENT(IN)  ::  cru_filename  !< filename aot raw data

    INTEGER (i8), INTENT(OUT) :: nrows !< number of rows
    INTEGER (i8), INTENT(OUT) :: ncolumns !< number of columns
    INTEGER (i8), INTENT(OUT) :: ntime !< number of times


    !local variables
    INTEGER :: ncid                             !< netcdf unit file number
    INTEGER :: ndimension                       !< number of dimensions in netcdf file
    INTEGER :: nVars                            !< number of variables in netcdf file
    INTEGER :: nGlobalAtts                      !< number of gloabal Attributes in netcdf file
    INTEGER :: unlimdimid                       !< id of unlimited dimension (e.g. time) in netcdf file

    INTEGER :: dimid                            !< id of dimension
    CHARACTER (len=80) :: dimname               !< name of dimensiona
    INTEGER :: length                           !< length of dimension

    ntime = 1

    ! open netcdf file 
    call check_netcdf( nf90_open(TRIM(cru_filename),NF90_NOWRITE, ncid))

    ! look for numbers of dimensions, Variable, Attributes, and the dimid for the unlimited dimension (probably time)
    call check_netcdf (nf90_inquire(ncid,ndimension, nVars, nGlobalAtts,unlimdimid))

    !; the dimid in netcdf-files is counted from 1 to ndimension
    !; look for the name and length of the dimension with f90_inquire_dimension
    !; nf90_inquire_dimension input: ncid, dimid; nf90_inquire_dimension output: name, length
    do dimid=1,ndimension
      !print *,'dimension loop dimid ',dimid
      call check_netcdf( nf90_inquire_dimension(ncid,dimid, dimname, length) )
      !print*, 'ncid,dimid, dimname, length',ncid,dimid, trim(dimname), length
      if ( trim(dimname) == 'lon') ncolumns=length          ! here I know that the name of zonal dimension is 'lon'
      if ( trim(dimname) == 'lat') nrows=length             ! here I know that the name of meridional dimension is 'lat'
      if ( trim(dimname) == 'time') ntime=length            ! here I know that the name of time dimension is "time"
    enddo

    ! close netcdf file 
    call check_netcdf( nf90_close( ncid))




  END SUBROUTINE get_dimension_cru_data


  !> get all aot data and coordinates and grid description
  SUBROUTINE get_cru_grid_and_data(cru_filename,   &
       raw_data_t_id,    &
       nrows,         &
       ncolumns,      &
       ntime)


    USE mo_cru_target_fields, ONLY: i_t_cru_fine, &
         &                              i_t_cru_coarse

    IMPLICIT NONE
    CHARACTER (LEN=*), INTENT(IN) :: cru_filename  !< filename aot raw data
    INTEGER (i8), INTENT(IN) :: raw_data_t_id    !< gives the data id (CRU fine (1) and CRU coarse (2))
    INTEGER (i8), INTENT(IN) :: nrows !< number of rows
    INTEGER (i8), INTENT(IN) :: ncolumns !< number of columns
    INTEGER (i8), INTENT(IN) :: ntime !< number of times

    ! the 'output' is via global variables, \TODO maybe change this i/o

    !local variables
    INTEGER :: ncid                             !< netcdf unit file number
    INTEGER :: ndimension                       !< number of dimensions in netcdf file
    INTEGER :: nVars                            !< number of variables in netcdf file
    INTEGER :: nGlobalAtts                      !< number of gloabal Attributes in netcdf file
    INTEGER :: unlimdimid                       !< id of unlimited dimension (e.g. time) in netcdf file

    CHARACTER (LEN=80) :: varname               !< name of variable
    !< new CRU temperature has 2 variables (temperature and elevation)
    !< old CRU temperatue has 1 variable (temperature)
    INTEGER :: varid                            !< id of variable
    INTEGER :: xtype                            !< netcdf type of variable/attribute
    INTEGER :: ndim                             !< number of dimensions of variable
    INTEGER, ALLOCATABLE :: var_dimids(:)       !< id of variable dimensions, vector, maximal dimension ndimension
    INTEGER :: nAtts                            !< number of attributes for a netcdf variable

    INTEGER :: errorcode                        !< error status variable

    CHARACTER (LEN=80) :: cooname(2) !< name of coordinates
    INTEGER :: coovarid(2)           !< varid of coordinats

    CHARACTER (LEN=80) :: attname !< name of attribute
    REAL :: scale_factor

    INTEGER :: n !< counter

    cooname(1) = 'lon'
    cooname(2) = 'lat'
    !        varname = 'tem' ! the name of the temperature variable in the netcdf file



    ! open netcdf file 
    call check_netcdf( nf90_open(TRIM(cru_filename),NF90_NOWRITE, ncid))

    DO n=1,2
      call check_netcdf( nf90_inq_varid(ncid, TRIM(cooname(n)), coovarid(n)))
    ENDDO

    CALL check_netcdf(nf90_get_var(ncid, coovarid(1),  lon_cru))
    CALL check_netcdf(nf90_get_var(ncid, coovarid(2),  lat_cru))

    CALL check_netcdf(nf90_inquire(ncid,ndimension, nVars, nGlobalAtts,unlimdimid))

    ALLOCATE (var_dimids(ndimension), STAT=errorcode)
    IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array var_dimids')
    var_dimids = 0

    variables: DO varid = 1,nVars
      CALL check_netcdf(nf90_inquire_variable(ncid,varid,varname,xtype, ndim, var_dimids, nAtts))
      IF (verbose >= idbg_high ) WRITE(logging%fileunit,*) 'varname: ', TRIM(varname)

      getvar: SELECT CASE(TRIM(varname))
      CASE ('tem')      ! the variable of the old CRU file is called tem
        CALL check_netcdf( nf90_get_var(ncid, varid, cru_raw_data))

        attname="scale_factor"
        CALL check_netcdf( nf90_get_att(ncid, varid, TRIM(attname), scale_factor))

      CASE ('T_CL')     ! the variable of the new CRU temperature is called T_CL
        CALL check_netcdf( nf90_get_var(ncid, varid, cru_raw_data))

      CASE ('HSURF')  ! the variable of the new CRU elevation is called HSURF
        CALL check_netcdf( nf90_get_var(ncid, varid, cru_raw_elev))

      END SELECT getvar

    ENDDO variables

    SELECT CASE(raw_data_t_id)
    CASE(i_t_cru_coarse)
      cru_raw_data = cru_raw_data * scale_factor
    END SELECT



    cru_grid%start_lon_reg = lon_cru(1)
    cru_grid%end_lon_reg = lon_cru(ncolumns)
    cru_grid%start_lat_reg = lat_cru(1)
    cru_grid%end_lat_reg = lat_cru(nrows)
    cru_grid%dlon_reg = (lon_cru(ncolumns) -  lon_cru(1) ) / (ncolumns - 1)
    cru_grid%dlat_reg = (lat_cru(nrows) - lat_cru(1) ) / (nrows -1) 
    cru_grid%nlon_reg = ncolumns
    cru_grid%nlat_reg = nrows




    ! close netcdf file 
    CALL check_netcdf( nf90_close(ncid))





  END SUBROUTINE get_cru_grid_and_data


  SUBROUTINE deallocate_cru_data()

    USE mo_cru_target_fields, ONLY: crutemp

    IMPLICIT NONE

    INTEGER :: errorcode !< error status variable

    DEALLOCATE (lon_cru, STAT=errorcode)
    IF(errorcode.NE.0) CALL abort_extpar('Cant deallocate the array lon_aot')
    DEALLOCATE (lat_cru, STAT=errorcode)
    IF(errorcode.NE.0) CALL abort_extpar('Cant deallocate the array lat_cru')
    DEALLOCATE (cru_raw_data, STAT=errorcode)
    IF(errorcode.NE.0) CALL abort_extpar('Cant deallocate the array cru_raw_data')
    DEALLOCATE (crutemp, STAT=errorcode)
    IF(errorcode.NE.0) CALL abort_extpar('Cant deallocate the array crutemp')
    DEALLOCATE (cru_raw_elev,STAT=errorcode)
    IF(errorcode.NE.0) CALL abort_extpar('Cant deallocate the array cru_raw_elev')
  END SUBROUTINE deallocate_cru_data



END MODULE mo_cru_data
