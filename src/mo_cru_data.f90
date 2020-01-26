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

  USE mo_logging          
  USE mo_kind,                  ONLY: wp, i4
  USE mo_io_units,              ONLY: filename_max  
  USE mo_io_utilities,          ONLY: check_netcdf
  USE mo_grid_structures,       ONLY: reg_lonlat_grid
  USE netcdf,                   ONLY: &
      &                               nf90_close,    &
      &                               nf90_get_att,  &
      &                               nf90_get_var,  &
      &                               nf90_inquire,  &
      &                               nf90_inquire_dimension, &
      &                               nf90_inquire_variable,  &
      &                               nf90_inq_varid,&
      &                               nf90_nowrite,  &
      &                               nf90_open
  
  USE mo_cru_target_fields,     ONLY: i_t_cru_coarse, &
      &                               crutemp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: cru_grid, &
       &    allocate_cru_data, &
       &    deallocate_cru_data, &
       &    read_cru_data_input_namelist, &
       &    read_namelists_extpar_t_clim, &
       &    get_dimension_cru_data, &
       &    get_cru_grid_and_data, &
       &    lon_cru, &
       &    lat_cru, &
       &    cru_raw_data, &
       &    cru_raw_elev

  TYPE(reg_lonlat_grid)       :: cru_grid !< structure with defenition of the raw data grid for the AOT dataset

  REAL (KIND=wp), ALLOCATABLE :: lon_cru(:), & !< longitude of aot grid
       &                         lat_cru(:), & !< latitude of aot grid
       &                         cru_raw_data(:,:,:), & !< aerosol optical thickness, aot(ncolumns,nrows,ntime) 
       &                         cru_raw_elev(:,:,:) !< surface height in cru (ncolumns,nrows,ntime)

  CONTAINS

  !---------------------------------------------------------------------------
  !> subroutine to read namelist for t_clim data settings for EXTPAR 
  SUBROUTINE read_namelists_extpar_t_clim(namelist_file,            &
       &                                  it_cl_type,               &
       &                                  raw_data_t_clim_path,     &
       &                                  raw_data_t_clim_filename, &
       &                                  t_clim_buffer_file,       &
       &                                  t_clim_output_file)

    CHARACTER (len=*), INTENT(IN)             :: namelist_file !< filename with namelists for for EXTPAR settings
    INTEGER (KIND=i4),      INTENT(OUT)       :: it_cl_type    !< integer switch to choose a land use raw data set
    ! 1 CRU fine (new), 2 CRU coarse (old) temperature climatology
    CHARACTER (len=filename_max), INTENT(OUT) :: raw_data_t_clim_path, &        !< path to raw data
         &                                       raw_data_t_clim_filename, &    !< filename temperature climatology raw data
         &                                       t_clim_buffer_file, & !< name for temperature climatology buffer
         &                                       t_clim_output_file !< name for temperature climatology output file
    

    INTEGER (KIND=i4)                         :: nuin, ierr

    !> namelist with filename for temperature climatlogy data output
    NAMELIST /t_clim_raw_data/ raw_data_t_clim_path, raw_data_t_clim_filename, it_cl_type

    !> namelist with filename for temperature climatlogy data output
    NAMELIST /t_clim_io_extpar/ t_clim_buffer_file, t_clim_output_file

    it_cl_type = -1
    
    raw_data_t_clim_path = ''
    raw_data_t_clim_filename = ''

    t_clim_buffer_file = ''
    t_clim_output_file = ''

    OPEN(NEWUNIT=nuin,FILE=TRIM(namelist_file), IOSTAT=ierr)
    IF (ierr /= 0) THEN
      CALL logging%error('CRU namelist open error ', __FILE__, __LINE__)
    ENDIF
    READ(nuin, NML=t_clim_raw_data, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      WRITE(0,NML=t_clim_raw_data)
      CALL logging%error('CRU raw data namelist read error ', __FILE__, __LINE__)      
    ENDIF
    READ(nuin, NML=t_clim_io_extpar, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      WRITE(0,NML=t_clim_io_extpar)
      CALL logging%error('CRU io namelist read error ', __FILE__, __LINE__)      
    ENDIF
    CLOSE(nuin)
    
  END SUBROUTINE read_namelists_extpar_t_clim
  !---------------------------------------------------------------------------

  !> subroutine to allocate aot data fields
  SUBROUTINE allocate_cru_data(nrows,ncolumns,ntime)

    IMPLICIT NONE
    INTEGER (KIND=i4), INTENT(IN) :: nrows, & !< number of rows
         &                           ncolumns, & !< number of columns
         &                           ntime !< number of times

    INTEGER (KIND=i4)             :: errorcode !< error status variable

    CALL logging%info('Enter routine: allocate_cru_data')

    ALLOCATE (lon_cru(1:ncolumns), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array lon_cru',__FILE__,__LINE__)
    lon_cru = 0.0

    ALLOCATE (lat_cru(1:nrows), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array lat_cru',__FILE__,__LINE__)
    lat_cru = 0.0

    ALLOCATE (cru_raw_data(1:ncolumns,1:nrows,1:ntime),STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array cru_raw_data',__FILE__,__LINE__)
    cru_raw_data = 0.0

    ALLOCATE (cru_raw_elev(1:ncolumns,1:nrows,1:ntime),STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array cru_raw_elev',__FILE__,__LINE__)
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

    INTEGER (KIND=i4), INTENT(OUT) :: nrows, & !< number of rows
         &                            ncolumns, & !< number of columns
         &                            ntime !< number of times

    !local variables
    INTEGER(KIND=i4)               :: ncid, &        !< netcdf unit file number
         &                            ndimension, &  !< number of dimensions in netcdf file
         &                            nVars, &       !< number of variables in netcdf file
         &                            nGlobalAtts, & !< number of gloabal Attributes in netcdf file
         &                            unlimdimid, &  !< id of unlimited dimension (e.g. time) in netcdf file
         &                            dimid, &       !< id of dimension
         &                            length         !< length of dimension

    CHARACTER (len=80)             :: dimname        !< name of dimensiona

    ntime = 1
    ! open netcdf file 
    CALL check_netcdf( nf90_open(TRIM(cru_filename),NF90_NOWRITE, ncid))

    ! look for numbers of dimensions, Variable, Attributes, and the dimid for the unlimited dimension (probably time)
    CALL check_netcdf (nf90_inquire(ncid,ndimension, nVars, nGlobalAtts,unlimdimid))

    DO dimid=1,ndimension
      CALL check_netcdf( nf90_inquire_dimension(ncid,dimid, dimname, length) )
      IF ( trim(dimname) == 'lon') ncolumns=length          ! here I know that the name of zonal dimension is 'lon'
      IF ( trim(dimname) == 'lat') nrows=length             ! here I know that the name of meridional dimension is 'lat'
      IF ( trim(dimname) == 'time') ntime=length            ! here I know that the name of time dimension is "time"
    ENDDO
    ! close netcdf file 
    CALL check_netcdf( nf90_close( ncid))

  END SUBROUTINE get_dimension_cru_data


  !> get all aot data and coordinates and grid description
  SUBROUTINE get_cru_grid_and_data(cru_filename,   &
       raw_data_t_id,    &
       nrows,         &
       ncolumns)

    IMPLICIT NONE

    CHARACTER (LEN=*), INTENT(IN) :: cru_filename  !< filename aot raw data
    INTEGER (KIND=i4), INTENT(IN) :: raw_data_t_id, &    !< gives the data id (CRU fine (1) and CRU coarse (2))
         &                           nrows, & !< number of rows
         &                           ncolumns !< number of columns

    !local variables
    INTEGER(KIND=i4)              :: ncid, &        !< netcdf unit file number
         &                           ndimension, &  !< number of dimensions in netcdf file
         &                           nVars, &       !< number of variables in netcdf file
         &                           nGlobalAtts, & !< number of gloabal Attributes in netcdf file
         &                           unlimdimid, &  !< id of unlimited dimension (e.g. time) in netcdf file
         &                           varid, &          !< id of variable
         &                           xtype, &          !< netcdf type of variable/attribute
         &                           ndim, &           !< number of dimensions of variable
         &                           errorcode, &      !< error status variable
         &                           coovarid(2), &    !< varid of coordinats
         &                           nAtts, &
         &                           n !< counter

    INTEGER, ALLOCATABLE          :: var_dimids(:)       !< id of variable dimensions, vector, maximal dimension ndimension

    CHARACTER (LEN=80)            :: varname, &    !< name of variable
        &                            cooname(2), & !< name of coordinates
        &                            attname !< name of attribute

    REAL                          :: scale_factor

    CALL logging%info('Enter routine: get_cru_grid_and_data')

    cooname(1) = 'lon'
    cooname(2) = 'lat'

    ! open netcdf file 
    call check_netcdf( nf90_open(TRIM(cru_filename),NF90_NOWRITE, ncid))

    DO n=1,2
      call check_netcdf( nf90_inq_varid(ncid, TRIM(cooname(n)), coovarid(n)))
    ENDDO

    CALL check_netcdf(nf90_get_var(ncid, coovarid(1),  lon_cru))
    CALL check_netcdf(nf90_get_var(ncid, coovarid(2),  lat_cru))

    CALL check_netcdf(nf90_inquire(ncid,ndimension, nVars, nGlobalAtts,unlimdimid))

    ALLOCATE (var_dimids(ndimension), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array var_dimids',__FILE__,__LINE__)
    var_dimids = 0

    variables: DO varid = 1,nVars
      CALL check_netcdf(nf90_inquire_variable(ncid,varid,varname,xtype, ndim, var_dimids, nAtts))

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

    IMPLICIT NONE

    INTEGER(KIND=i4) :: errorcode !< error status variable

    CALL logging%info('Enter routine: deallocate_cru_data')

    DEALLOCATE (lon_cru, STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant deallocate the array lon_aot',__FILE__,__LINE__)
    DEALLOCATE (lat_cru, STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant deallocate the array lat_cru',__FILE__,__LINE__)
    DEALLOCATE (cru_raw_data, STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant deallocate the array cru_raw_data',__FILE__,__LINE__)
    DEALLOCATE (crutemp, STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant deallocate the array crutemp',__FILE__,__LINE__)
    DEALLOCATE (cru_raw_elev,STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant deallocate the array cru_raw_elev',__FILE__,__LINE__)

  END SUBROUTINE deallocate_cru_data

END MODULE mo_cru_data
