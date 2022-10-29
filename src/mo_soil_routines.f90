!+ Fortran module with for soil data specification subroutines
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_0         2010/12/21 Hermann Asensio
!  Initial release
! V1_10        2013/03/26 Juergen Helmert
!  Adaption in DSMW grid increment computation         
! V2_0         1013/06/04 Martina Messmer
!  adaptations in a way that HWSD data set can be used (top- and subsoil)
!  Code received from Juergen Helmert
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!> Fortran module with for soil data specification subroutines
!> \author Hermann Asensio
!>
MODULE mo_soil_routines


  USE mo_logging
  USE mo_kind,                  ONLY: wp, i4

  USE netcdf,                   ONLY:   &
                                      nf90_open,              &
                                      nf90_close,             &
                                      nf90_inquire,           &
                                      nf90_inquire_dimension, &
                                      nf90_inquire_variable,  &
                                      nf90_inquire_dimension, &
                                      nf90_get_var,           &
                                      nf90_nowrite

  USE mo_io_utilities,          ONLY: check_netcdf
                               
  USE mo_io_units,              ONLY: filename_max
                               
  USE mo_utilities_extpar,      ONLY: free_un ! function to get free unit number

  USE mo_soil_data,             ONLY: lon_full, &       !< longitide coordinates of the regular grid in the geographical (lonlat)
       &                              lat_full, &       !< latitude coordinates of the regular grid in the geographical (lonlat)
       &                              lon_soil, &
       &                              lat_soil, &
       &                              soil_texslo, &
       &                              dsmw_soil_unit, &
       &                              n_unit, &
       &                              dsmw_grid

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: get_soil_data
  PUBLIC :: get_dimension_soil_data
  PUBLIC :: read_namelists_extpar_soil

  SAVE
  PUBLIC :: nlon_soil, nlat_soil, n_unit

  INTEGER  (KIND=i4) :: nlon_soil, nlat_soil

  CONTAINS


  !---------------------------------------------------------------------------
  !> subroutine to read namelist for soil data settings for EXTPAR 
  SUBROUTINE read_namelists_extpar_soil(namelist_file,                 &
                                         isoil_data,                 &
                                         raw_data_soil_path,         &
                                         raw_data_soil_filename,     &
                                         soil_buffer_file,           &
                                         soil_buffer_file_consistent,&
                                         soil_output_file_consistent)
  
    CHARACTER (len=*), INTENT(IN)  :: namelist_file !< filename with namelists for for EXTPAR settings
    
    INTEGER (KIND=i4), INTENT(OUT) :: isoil_data !< 1 for FAO data , 2 for HWSD
    
    CHARACTER (len=filename_max)   :: raw_data_soil_path, &        !< path to raw data
         &                            raw_data_soil_filename, & !< filename soil raw data
         &                            soil_buffer_file, &  !< name for soil buffer file
         &                            soil_buffer_file_consistent, & !< name for soil buffer file after consistency check
         &                            soil_output_file_consistent, & !< name for soil output file after consistency check
         &                            raw_data_deep_soil_filename

    INTEGER (KIND=i4)              :: ierr, nuin

    LOGICAL                        :: ldeep_soil

    !>Define the namelist group for soil raw data
    NAMELIST /soil_raw_data/ isoil_data, ldeep_soil, raw_data_soil_path, raw_data_soil_filename, raw_data_deep_soil_filename

    !> namelist with filenames for output of soil data
    NAMELIST /soil_io_extpar/ soil_buffer_file, soil_buffer_file_consistent, soil_output_file_consistent

    nuin = free_un()  ! functioin free_un returns free Fortran unit number
    OPEN(nuin,FILE=TRIM(namelist_file), IOSTAT=ierr)
    IF (ierr /= 0) THEN
      WRITE(message_text,*)'Cannot open ', TRIM(namelist_file)
      CALL logging%error(message_text,__FILE__, __LINE__) 
    ENDIF

    READ(nuin, NML=soil_raw_data, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      CALL logging%error('Cannot read in namelist soil_raw_data',__FILE__, __LINE__) 
    ENDIF

    READ(nuin, NML=soil_io_extpar, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      CALL logging%error('Cannot read in namelist soil_io_extpar',__FILE__, __LINE__) 
    ENDIF

    CLOSE(nuin)

    IF (ldeep_soil) THEN
      CALL logging%warning('Option ldeep_soil was removed with Release 5.11')
      CALL logging%error('ldeep_soil not supported anymore',__FILE__, __LINE__) 
    ENDIF
  
  END SUBROUTINE read_namelists_extpar_soil

    !> inquire dimension information for NDVI raw data 
  SUBROUTINE get_dimension_soil_data(path_soil_file, &
                                      nlon_full, &
                                      nlat_full, &
                                      n_unit)

    CHARACTER (len=*), INTENT(IN)  :: path_soil_file         !< filename with path for soil raw data

    INTEGER (KIND=i4), INTENT(OUT) :: nlon_full, & !< number of grid elements in zonal direction for soil data
         &                            nlat_full, & !< number of grid elements in meridional direction for soil data
         &                            n_unit    !< number of soil unit numbers (legend)

    !local variables
    INTEGER(KIND=i4)               :: ncid, &                             !< netcdf unit file number
         &                            ndimension, &                       !< number of dimensions in netcdf file
         &                            nVars, &                            !< number of variables in netcdf file
         &                            nGlobalAtts, &                      !< number of gloabal Attributes in netcdf file
         &                            unlimdimid, &                       !< id of unlimited dimension (e.g. time) in netcdf file
         &                            dimid, &                            !< id of dimension
         &                            length, &                           !< length of dimension
         &                            xtype, &                            !< netcdf type of variable/attribute
         &                            ndim, &                             !< number of dimensions of variable
         &                            varid, &                            !< id of variable
         &                            nAtts                            !< number of attributes for a netcdf variable

    INTEGER(KIND=i4), ALLOCATABLE  :: var_dimids(:)       !< id of variable dimensions, vector, maximal dimension ndimension

    CHARACTER (LEN=80)             :: varname, &               !< name of variable
         &                            dimname               !< name of dimensiona
    
    ! open netcdf file 
    CALL check_netcdf(nf90_open(TRIM(path_soil_file),NF90_NOWRITE, ncid))
    
    ! look for numbers of dimensions, Variable, Attributes, and the dimid for the one possible unlimited dimension 
    ! (probably time)
    !; nf90_inquire input: ncid; nf90_inquire output: ndimension, nVars, nGlobalAtts,unlimdimid
    CALL check_netcdf (nf90_inquire(ncid,ndimension, nVars, nGlobalAtts,unlimdimid))
    ALLOCATE (var_dimids(ndimension))

    !; the dimid in netcdf-files is counted from 1 to ndimension
    !; look for the name and length of the dimension with f90_inquire_dimension
    !; nf90_inquire_dimension input: ncid, dimid; nf90_inquire_dimension output: name, length
    DO dimid=1,ndimension
      CALL check_netcdf( nf90_inquire_dimension(ncid,dimid, dimname, length) )
      IF ( TRIM(dimname) == 'lon') nlon_full=length          ! here I know that the name of zonal dimension is 'lon'
      IF ( TRIM(dimname) == 'lat') nlat_full=length          ! here I know that the name of meridional dimension is 'lat'
      IF ( TRIM(dimname) == 'soil_unit') n_unit=length        ! here I know that the name of unit dimension is 'soil_unit'
    ENDDO
    
    ! Allocate lon_full, lat_full
    ALLOCATE (lon_full(nlon_full), lat_full(nlat_full))
    
    ! read lon_full and lat_full from raw data file
    variables: DO varid=1,nVars
      CALL check_netcdf(nf90_inquire_variable(ncid,varid,varname,xtype, ndim, var_dimids, nAtts))
    
      getvar: SELECT CASE(TRIM(varname))
    
          CASE('lon')     !  here I know that the variable with the longitude coordinates is called 'lon'
            CALL check_netcdf(nf90_get_var(ncid,varid,lon_full(:)) )
              
          CASE('lat')     !  here I know that the variable with the latitude coordinates is called 'lat'
            CALL check_netcdf(nf90_get_var(ncid,varid,lat_full(:)) )
          
      END SELECT getvar
    ENDDO variables
    
    dsmw_grid%nlon_reg = nlon_full
    dsmw_grid%nlat_reg = nlat_full
    
    dsmw_grid%start_lon_reg = lon_full(1)
    dsmw_grid%end_lon_reg   = lon_full(nlon_full)
    
    dsmw_grid%start_lat_reg = lat_full(1)
    dsmw_grid%end_lat_reg   = lat_full(nlat_full)
    
    IF (dsmw_grid%nlon_reg /= 0) THEN
      dsmw_grid%dlon_reg = (dsmw_grid%end_lon_reg - dsmw_grid%start_lon_reg)/(dsmw_grid%nlon_reg-1._wp)
    ENDIF 
    
    IF (dsmw_grid%nlat_reg /= 0) THEN
      dsmw_grid%dlat_reg = (dsmw_grid%end_lat_reg - dsmw_grid%start_lat_reg)/(dsmw_grid%nlat_reg-1._wp)
    ENDIF ! in case of latitude orientation from north to south dlat is negative!
    
    ! close netcdf file 
    CALL check_netcdf( nf90_close(ncid))

  END SUBROUTINE get_dimension_soil_data

  !> get coordintates, legend and data for soil raw data 
  SUBROUTINE get_soil_data(path_soil_file,start)

    CHARACTER (len=*), INTENT(in) :: path_soil_file                !< filename with path for soil raw data
    INTEGER , INTENT(IN)          :: start(2)

    !local variables
    INTEGER(KIND=i4)              :: ncid, &                             !< netcdf unit file number
         &                           ndimension, &                       !< number of dimensions in netcdf file
         &                           nVars, &                            !< number of variables in netcdf file
         &                           nGlobalAtts, &                      !< number of gloabal Attributes in netcdf file
         &                           unlimdimid, &                       !< id of unlimited dimension (e.g. time) in netcdf file
         &                           dimid, &                            !< id of dimension
         &                           length, &                           !< length of dimension
         &                           varid, &                            !< id of variable
         &                           xtype, &                            !< netcdf type of variable/attribute
         &                           ndim, &                             !< number of dimensions of variable
         &                           nAtts, &                            !< number of attributes for a netcdf variable
         &                           errorcode                        !< error status variable

    INTEGER, ALLOCATABLE          :: var_dimids(:)       !< id of variable dimensions, vector, maximal dimension ndimension
    
    CHARACTER (LEN=80)            :: dimname, &               !< name of dimension
         &                           varname               !< name of variable

    ! open netcdf file
    CALL check_netcdf( nf90_open(TRIM(path_soil_file),NF90_NOWRITE, ncid))
    ! look for numbers of dimensions, Variable, Attributes, and the dimid for the unlimited dimension (probably time)
    ! nf90_inquire input: ncid; nf90_inquire output: ndimension, nVars, nGlobalAtts,unlimdimid
    CALL check_netcdf (nf90_inquire(ncid,ndimension, nVars, nGlobalAtts,unlimdimid))

    ALLOCATE (var_dimids(ndimension), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array var_dimids',__FILE__,__LINE__)
    var_dimids = 0
    !; the dimid in netcdf-files is counted from 1 to ndimension
    !; look for the name and length of the dimension with f90_inquire_dimension
    !; nf90_inquire_dimension input: ncid, dimid; nf90_inquire_dimension output: name, length
    DO dimid=1,ndimension
      CALL check_netcdf( nf90_inquire_dimension(ncid,dimid, dimname, length) )
      IF ( TRIM(dimname) == 'soil_unit') n_unit=length        ! here I know that the name of unit dimension is 'soil_unit'
    ENDDO
    !; the varid in netcdf files is counted from 1 to nVars
    !; look for variable names, number of dimension, var_dimids etc with nf90_inquire_variable
    !; nf90_inquire_variable input: ncid, varid; nf90_inquire_variable output name xtype, ndim, var_dimids, nAtts

    variables: DO varid=1,nVars
      CALL check_netcdf(nf90_inquire_variable(ncid,varid,varname,xtype, ndim, var_dimids, nAtts))
      getvar: SELECT CASE(TRIM(varname))
             
          CASE('code')    !  here I know that the variable with the dsmw_code is called 'code'
            CALL check_netcdf(nf90_get_var(ncid,varid,soil_texslo(:)%dsmw_code) )
             
          CASE('tex_coarse')!  here I know that the variable with the  coarse texture part is called 'tex_coarse'
            CALL check_netcdf(nf90_get_var(ncid,varid,soil_texslo(:)%tex_coarse) )
             
          CASE('tex_medium')!  here I know that the variable with the medium texture part is called 'tex_medium'
            CALL check_netcdf(nf90_get_var(ncid,varid,soil_texslo(:)%tex_medium) )
             
          CASE('tex_fine')!  here I know that the variable with the fine texture part is called 'tex_fine'
            CALL check_netcdf(nf90_get_var(ncid,varid,soil_texslo(:)%tex_fine) )
             
          CASE('undefined')!  here I know that the variable with the undefined soil part is called 'undefined'
            CALL check_netcdf(nf90_get_var(ncid,varid,soil_texslo(:)%part_undefined) )
             
          CASE('flat')!  here I know that the variable with the flat area part is called 'flat'
            CALL check_netcdf(nf90_get_var(ncid,varid,soil_texslo(:)%flat) )
             
          CASE('hilly')!  here I know that the variable with the hilly area part is called 'hilly'
            CALL check_netcdf(nf90_get_var(ncid,varid,soil_texslo(:)%hilly) )
             
          CASE('steep')!  here I know that the variable with the steep area part is called 'steep'
            CALL check_netcdf(nf90_get_var(ncid,varid,soil_texslo(:)%steep) )                         
             
          CASE('DSMW')  !  here I know that the variable with the DSMW soil units is called 'DSMW'
            CALL check_netcdf(nf90_get_var(ncid,varid,dsmw_soil_unit,start=start,count=(/ nlon_soil, nlat_soil /)))
             
          CASE('Soil')  !  here I know that the variable with the DSMW soil units is called 'Soil'
            CALL check_netcdf(nf90_get_var(ncid,varid,dsmw_soil_unit,start=start,count=(/ nlon_soil, nlat_soil /)) )
             
      END SELECT getvar
    ENDDO variables
         
    ! close netcdf file 
    CALL check_netcdf( nf90_close( ncid))
         
    ! Fill the structure soil_raw_data_grid
    dsmw_grid%nlon_reg = nlon_soil
    dsmw_grid%nlat_reg = nlat_soil

    dsmw_grid%start_lon_reg = lon_soil(1)
    dsmw_grid%end_lon_reg   = lon_soil(nlon_soil)

    dsmw_grid%start_lat_reg = lat_soil(1)
    dsmw_grid%end_lat_reg   = lat_soil(nlat_soil)
       
  END SUBROUTINE get_soil_data
     !------------------------------------------------------------------------------------------------------------

END MODULE mo_soil_routines
