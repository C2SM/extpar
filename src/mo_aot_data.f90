!+  Fortran module with data fields for Aerosol optical thickness data
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_0         2010/12/21 Hermann Asensio
!  Initial release
! V4_0         2016/08/17 authors from RHM and Daniel Lthi
!  Added support for MACv2 aerosol data fields
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!> Fortran module with data fields for Aerosol optical thickness data
!> \author Hermann Asensio
!
MODULE mo_aot_data

  !> kind parameters are defined in MODULE data_parameters
  USE mo_logging
  USE mo_kind,                  ONLY: wp,i4

  USE netcdf,                   ONLY: nf90_open,              &
                                      nf90_close,             &
                                      nf90_inquire,           &
                                      nf90_inquire_dimension, &
                                      nf90_inquire_dimension, &
                                      nf90_inq_varid,         &
                                      nf90_get_var,           &
                                      nf90_nowrite

  USE mo_utilities_extpar,      ONLY: free_un

  USE mo_io_utilities,          ONLY: check_netcdf

  USE mo_grid_structures,       ONLY: reg_lonlat_grid

  USE mo_io_units,              ONLY: filename_max

  USE mo_aot_target_fields,     ONLY: aot_tg, &
    &                                 MAC_aot_tg, &
    &                                 MAC_ssa_tg, &
    &                                 MAC_asy_tg
                             
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: aot_grid, ntime_aot, nspb_aot, iaot_type

  PUBLIC :: allocate_aot_data, &
            deallocate_aot_data, &
            read_namelists_extpar_aerosol, &
            read_aot_data_input_namelist, &
            get_dimension_aot_data, &
            get_aot_grid_and_data, &
            lon_aot, &
            lat_aot, &
            aot_data, &
            ntype_aot, &
            n_spectr,  &
            aot_varname, &
            aot_longname, &
            aot_shortname, &
            MAC_data    !------new kinne-----

  TYPE(reg_lonlat_grid)       :: aot_grid !< definition of the raw data grid for the AOT dataset
  REAL (KIND=wp), ALLOCATABLE :: lon_aot(:), & !< longitude of aot grid
    &                            lat_aot(:), & !< latitude of aot grid
    &                            aot_data(:,:,:,:), & !< aerosol optical thickness, aot(ntype,ncolumns,nrows,ntime) 
    &                            MAC_data(:,:,:,:,:) !< normalized aerosol optical properties, aot(ntype,ncolumns,nrows,ntime,itype) 

  INTEGER (KIND=i4), PARAMETER :: ntype_aot = 5, & !< number of types of aerosols
    &                             ntime_aot = 12, & !< 12 monthly mean data of aeorsol optical thickness
    &                             nspb_aot = 9 !< 9 spectral bands of aeorsol optical thickness
  INTEGER (KIND=i4)             :: n_spectr

  CHARACTER (len=32)            :: aot_varname(ntype_aot) = &    !< variable name for aerosolt type
    &                               (/ 'bc   ', 'dust ', 'org  ', 'so4  ', 'ssalt' /)
  CHARACTER (len=80)            :: aot_longname(ntype_aot) = &   !< long name for aereosol type
    &                               (/ 'aerosol optical thickness of black carbon  ', &
    &                               'aerosol optical thickness of dust          ', &
    &                               'aerosol optical thickness of organic matter', &
    &                               'aerosol optical thickness of sulfate       ', &
    &                               'aerosol optical thickness of sea salt      ' /)
  CHARACTER (len=20)            :: aot_shortname(ntype_aot)= &   !< short name for aereosol type
    &                               (/ 'AER_BC              ', &
    &                                'AER_DUST            ', &
    &                                'AER_ORG             ', &
    &                                'AER_SO4             ', &
    &                                'AER_SS              ' /)

  INTEGER (KIND=i4)             :: iaot_type = 1

  CONTAINS

  !---------------------------------------------------------------------------
  !> subroutine to read namelist for aerosol optical thickness data settings for EXTPAR 
  SUBROUTINE read_namelists_extpar_aerosol(namelist_file, &
    &                                      iaot_type,    &
    &                                    raw_data_aot_path, &
    &                                    raw_data_aot_filename, &
    &                                    aot_buffer_file)
  
    CHARACTER (LEN=*), INTENT(IN)            :: namelist_file !< filename with namelists for for EXTPAR settings

    ! aerosol optical thickness
    CHARACTER (LEN=filename_max)             :: raw_data_aot_path, &
      &                                         raw_data_aot_filename, &
      &                                         aot_buffer_file, &
      &                                         filename
    INTEGER (KIND=i4)                        :: iaot_type, nuin, ierr

!> namelist with filenames for aerosol optical thickness data input
    NAMELIST /aerosol_raw_data/ raw_data_aot_path, raw_data_aot_filename, iaot_type

!> namelist with filenames for aerosol optical thickness data output
    NAMELIST /aerosol_io_extpar/ aot_buffer_file

   nuin = free_un()  ! functioin free_un returns free Fortran unit number
   filename = TRIM(namelist_file)

   OPEN(nuin,FILE=filename, IOSTAT=ierr)

   READ(nuin, NML=aerosol_raw_data, IOSTAT=ierr)
   READ(nuin, NML=aerosol_io_extpar, IOSTAT=ierr)

   CLOSE(nuin)

    IF (ierr /= 0) THEN
      WRITE(message_text,*)'Cannot read ', filename
      CALL logging%error(message_text,__FILE__, __LINE__) 
    ENDIF

  END SUBROUTINE read_namelists_extpar_aerosol
!---------------------------------------------------------------------------

!> subroutine to allocate aot data fields
  SUBROUTINE allocate_aot_data(iaot_type,nrows,ncolumns,ntime,ntype,n_spectr)

  IMPLICIT NONE

  INTEGER (KIND=i4), INTENT(IN) :: ntype, & !< number of types of aerosols
    &                              nrows, & !< number of rows
    &                              ncolumns, & !< number of columns
    &                              ntime, & !< number of times
    &                              n_spectr, & !< number of times
    &                              iaot_type !< if =4 MACv2 new

  INTEGER                       :: errorcode !< error status variable

    CALL logging%info('Enter routine: allocate_aot_data')

    ALLOCATE (lon_aot(1:ncolumns+1), STAT=errorcode)
!DWD   ALLOCATE (lon_aot(1:ncolumns), STAT=errorcode)
        IF(errorcode.NE.0) CALL logging%error('Cant allocate the array lon_aot',__FILE__,__LINE__)
    lon_aot = 0.0

     ALLOCATE (lat_aot(1:nrows), STAT=errorcode)
        IF(errorcode.NE.0) CALL logging%error('Cant allocate the array lat_aot',__FILE__,__LINE__)
    lat_aot = 0.0

    IF (iaot_type.NE.4) THEN
      ALLOCATE (aot_data(1:ncolumns+1,1:nrows,1:ntime,1:ntype),STAT=errorcode)
      IF(errorcode.NE.0) CALL logging%error('Cant allocate the array aot_data',__FILE__,__LINE__)
      aot_data = 0.0
    ELSE
      ALLOCATE (MAC_data(1:ncolumns+1,1:nrows,1:n_spectr,1:ntime,1:ntype),STAT=errorcode)
      IF(errorcode.NE.0) CALL logging%error('Cant allocate the array MAC_data',__FILE__,__LINE__)
      MAC_data = 0.0
    ENDIF

  END SUBROUTINE allocate_aot_data

  !> subroutine to read namelist with settings for aerosol data
  !! \author Hermann Asensio
  SUBROUTINE read_aot_data_input_namelist(input_namelist_file, aot_filename)
  

    CHARACTER (LEN=*), INTENT(IN)             :: input_namelist_file !< file with input namelist 
    CHARACTER (LEN=filename_max), INTENT(OUT) :: aot_filename  !< filename aot raw data

    INTEGER (KIND=i4)                         :: ierr, nuin

    !>Define the namelist group
    NAMELIST /AOT_file_info/ aot_filename

    nuin = free_un()  ! functioin free_un returns free Fortran unit number
    OPEN(nuin,FILE=TRIM(input_namelist_file), IOSTAT=ierr)
    IF (ierr /= 0) THEN
      WRITE(message_text,*)'Cannot open ', TRIM(input_namelist_file)
      CALL logging%error(message_text,__FILE__, __LINE__) 
    ENDIF

    READ(nuin, NML=AOT_file_info, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      CALL logging%error('Cannot read in namelist AOT_file_info',__FILE__, __LINE__) 
    ENDIF
    CLOSE(nuin)


   END SUBROUTINE read_aot_data_input_namelist

   !> get dimension information of aot data from netcdf file
   SUBROUTINE get_dimension_aot_data(aot_filename, &
     &                               iaot_type,    &
     &                               nrows,        &
     &                               ncolumns,     &
     &                               ntime,        &
     &                               ntype,        &
     &                               n_spectr)

     IMPLICIT NONE

     CHARACTER (LEN=*), INTENT(IN)           ::  aot_filename  !< filename aot raw data
     INTEGER (KIND=i4), INTENT(IN)           :: iaot_type !< if =4 MACv2

     INTEGER (KIND=i4), INTENT(OUT)          :: ntype, & !< number of types of aerosols
       &                                        nrows, & !< number of rows
       &                                        ncolumns, & !< number of columns
       &                                        ntime, & !< number of times
       &                                        n_spectr !< number of times

       !local variables
     CHARACTER (LEN=filename_max)            :: filename
     CHARACTER (LEN=80)                      :: dimname              
     INTEGER                                 :: ncid, ndimension, nVars, &
       &                                        nGlobalAtts, unlimdimid, dimid, &
       &                                        length

     CALL logging%info('Enter routine: get_dimension_aot_data')

      ! open netcdf file
     filename = TRIM(aot_filename)
     CALL check_netcdf( nf90_open(filename,NF90_NOWRITE, ncid))
     
     ! look for numbers of dimensions, Variable, Attributes, and the dimid for the unlimited dimension (probably time)
     CALL check_netcdf (nf90_inquire(ncid,ndimension, nVars, nGlobalAtts,unlimdimid))

     DO dimid=1,ndimension
       CALL check_netcdf( nf90_inquire_dimension(ncid,dimid, dimname, length) )
       IF ( TRIM(dimname) == 'lon') ncolumns=length 
       IF ( TRIM(dimname) == 'lat') nrows=length    
       IF ( TRIM(dimname) == 'time') ntime=length   
     ENDDO

     IF (iaot_type == 4) THEN
       ntype=3
       n_spectr=9
     ELSE
       ntype=ntype_aot 
       n_spectr=1
     ENDIF

     ! close netcdf file 
     CALL check_netcdf( nf90_close( ncid))

   END SUBROUTINE get_dimension_aot_data 

   !> get all aot data and coordinates and grid description
  SUBROUTINE get_aot_grid_and_data(iaot_type, &
    &                                aot_filename, &
    &                                nrows,        &
    &                                ncolumns,     &
    &                                ntime,        &
    &                                ntype,        &
    &                                n_spectr,     &
    &                                aot_grid,     &
    &                                lon_aot,      &
    &                                lat_aot,      &
    &                                aot_data,     &
                                     MAC_data)
    IMPLICIT NONE

    CHARACTER (LEN=*), INTENT(IN)             ::  aot_filename  !< filename aot raw data
    INTEGER (KIND=i4), INTENT(IN)             :: iaot_type, & !< if =0 MACv2 new
      &                                          ntype, &  !< number of types of aerosols
      &                                          nrows, &  !< number of rows
      &                                          ncolumns, &  !< number of columns
      &                                          ntime, &  !< number of times
      &                                          n_spectr !< number of spectral bands
    
    TYPE(reg_lonlat_grid), INTENT(INOUT)      :: aot_grid
    
    REAL (KIND=wp), INTENT(INOUT)             :: lon_aot(1:ncolumns+1), & 
                                                 lat_aot(1:nrows), & 
                                                 aot_data(:,:,:,:), & 
                                                 MAC_data(:,:,:,:,:)
    !local variables
    REAL, ALLOCATABLE                         :: aot_data_stype(:,:,:), &
      &                                          MAC_data_stype(:,:,:,:)

    INTEGER                                   :: ncid,n, coovarid(2), &
      &                                          varid(ntype)

    CHARACTER (LEN=80)                        :: varname(ntype), &  !< name of variable
      &                                          cooname(2) !< name of coordinates
    
    CALL logging%info('Enter routine: get_aot_grid_and_data')

    cooname(1) = 'lon'
    cooname(2) = 'lat'

    ! I know the names of tha variables already
    IF (iaot_type == 4) THEN
      varname(1) = 'AOT'
      varname(2) = 'SSA'
      varname(3) = 'ASY'
    ELSE
      varname(1) = 'black_carbon'
      varname(2) = 'dust'
      varname(3) = 'organic'
      varname(4) = 'sulfate'
      varname(5) = 'sea_salt'
    ENDIF

    ! open netcdf file 
    CALL check_netcdf( nf90_open(TRIM(aot_filename),NF90_NOWRITE, ncid))

    DO n=1,2
      CALL check_netcdf( nf90_inq_varid(ncid, TRIM(cooname(n)), coovarid(n)))
    ENDDO

    CALL check_netcdf(nf90_get_var(ncid, coovarid(1),  lon_aot(1:ncolumns)))
    CALL check_netcdf(nf90_get_var(ncid, coovarid(2),  lat_aot))

    IF (iaot_type == 4) THEN
      ALLOCATE (MAC_data_stype(ncolumns,nrows,n_spectr,ntime))
      DO n=1,ntype
        CALL check_netcdf( nf90_inq_varid(ncid, TRIM(varname(n)), varid(n)))

        CALL check_netcdf(nf90_get_var(ncid, varid(n),  MAC_data_stype))

        MAC_data(1:ncolumns,:,:,:,n) = MAC_data_stype(1:ncolumns,:,:,:)
      ENDDO
      MAC_data(ncolumns+1,:,:,:,:) = MAC_data(1,:,:,:,:)
      DEALLOCATE (MAC_data_stype)
    ELSE
      ALLOCATE (aot_data_stype(ncolumns,nrows,ntime))
    DO n=1,ntype
        CALL check_netcdf( nf90_inq_varid(ncid, TRIM(varname(n)), varid(n)))

    CALL check_netcdf(nf90_get_var(ncid, varid(n),  aot_data_stype))

        aot_data(1:ncolumns,:,:,n) = aot_data_stype(1:ncolumns,:,:)
      ENDDO
      aot_data(ncolumns+1,:,:,:) = aot_data(1,:,:,:)
      DEALLOCATE (aot_data_stype)
    ENDIF

    CALL check_netcdf( nf90_close( ncid))
    ! close netcdf file 
    
    IF (iaot_type /= 4) THEN
     ! extend aot_data by 1 column so that the field covers the whole globe
     aot_data(ncolumns+1,:,:,:) = aot_data(1,:,:,:)
    ENDIF
     ! set aot_grid values
    aot_grid%start_lon_reg = lon_aot(1)
    aot_grid%start_lat_reg = lat_aot(1)
    aot_grid%dlon_reg = (lon_aot(ncolumns) -  lon_aot(1) ) / (ncolumns - 1)
    aot_grid%dlat_reg = (lat_aot(nrows) - lat_aot(1) ) / (nrows -1) 
    aot_grid%nlon_reg = ncolumns+1
    aot_grid%nlat_reg = nrows

    aot_grid%end_lon_reg = lon_aot(ncolumns) + aot_grid%dlon_reg
    aot_grid%end_lat_reg = lat_aot(nrows)
    lon_aot(ncolumns+1)=lon_aot(ncolumns) + aot_grid%dlon_reg

  END SUBROUTINE get_aot_grid_and_data

  SUBROUTINE deallocate_aot_data()

    IMPLICIT NONE
   
    INTEGER :: errorcode

    CALL logging%info('Enter routine: deallocate_aot_data')

    DEALLOCATE (lon_aot, STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant deallocate the array lon_aot',__FILE__,__LINE__)
    DEALLOCATE (lat_aot, STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant deallocate the array lat_aot',__FILE__,__LINE__)

    IF (iaot_type == 4) THEN
      DEALLOCATE (MAC_data, STAT=errorcode)
      IF(errorcode.NE.0) CALL logging%error('Cant deallocate the array MAC_data',__FILE__,__LINE__)
      DEALLOCATE (MAC_aot_tg, STAT=errorcode)
      IF(errorcode.NE.0) CALL logging%error('Cant deallocate the array MAC_aot_tg',__FILE__,__LINE__)
      DEALLOCATE (MAC_ssa_tg, STAT=errorcode)
      IF(errorcode.NE.0) CALL logging%error('Cant deallocate the array MAC_ssa_tg',__FILE__,__LINE__)
      DEALLOCATE (MAC_asy_tg, STAT=errorcode)
      IF(errorcode.NE.0) CALL logging%error('Cant deallocate the array MAC_asy_tg',__FILE__,__LINE__)
    ELSE
    DEALLOCATE (aot_data, STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant deallocate the array aot_data',__FILE__,__LINE__)
    DEALLOCATE (aot_tg, STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant deallocate the array aot_tg',__FILE__,__LINE__)
    ENDIF

  END SUBROUTINE deallocate_aot_data

END MODULE mo_aot_data
