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


 USE mo_kind, ONLY: wp
 USE mo_kind, ONLY: i4
 USE mo_kind, ONLY: i8

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

USE netcdf,      ONLY:    &
  nf90_create,            &
  nf90_def_dim,           &
  nf90_def_var,           &
  nf90_enddef,            &
  nf90_redef,             &
  nf90_put_att,           &
  nf90_put_var

 
USE netcdf,      ONLY :   &
  NF90_CHAR,              &
  NF90_DOUBLE,            &
  NF90_FLOAT,             &
  NF90_INT,               &
  NF90_BYTE,              &
  NF90_SHORT


USE netcdf,      ONLY :   &
  NF90_GLOBAL,            &
  NF90_UNLIMITED,         &
  NF90_CLOBBER,           &
  NF90_NOWRITE


!> abort_extpar defined in MODULE utilities_extpar
USE mo_utilities_extpar, ONLY: abort_extpar


USE mo_io_utilities,     ONLY: check_netcdf

USE mo_io_units,         ONLY: filename_max

USE mo_grid_structures,  ONLY: reg_lonlat_grid

USE mo_soil_data,        ONLY: dsmw_legend, n_unit


IMPLICIT NONE

PRIVATE

PUBLIC :: get_soil_data
PUBLIC :: get_deep_soil_data
PUBLIC :: get_dimension_soil_data
PUBLIC :: read_namelists_extpar_soil

SAVE
PUBLIC :: nlon_soil, nlat_soil, n_unit

INTEGER  (KIND=i4) :: nlon_soil
INTEGER  (KIND=i4) :: nlat_soil

CONTAINS


!---------------------------------------------------------------------------
!> subroutine to read namelist for soil data settings for EXTPAR 
  SUBROUTINE read_namelists_extpar_soil(namelist_file,                 &
                                         isoil_data,                 &
                                         ldeep_soil,                 &
                                         raw_data_soil_path,         &
                                         raw_data_soil_filename,     &
                                         raw_data_deep_soil_filename,&
                                         soil_buffer_file,           &
                                         soil_output_file,           &
                                         soil_buffer_file_consistent,&
                                         soil_output_file_consistent)

    USE mo_utilities_extpar, ONLY: free_un ! function to get free unit number

  
    CHARACTER (len=filename_max), INTENT(IN) :: namelist_file !< filename with namelists for for EXTPAR settings

    INTEGER (KIND=i4), INTENT(OUT) :: isoil_data !< 1 for FAO data , 2 for HWSD
    LOGICAL,           INTENT(OUT) :: ldeep_soil !< can only be true if HWSD data is chosen

    ! soil
    CHARACTER (len=filename_max) :: raw_data_soil_path        !< path to raw data
    CHARACTER (len=filename_max) :: raw_data_soil_filename !< filename soil raw data
    CHARACTER (len=filename_max) :: raw_data_deep_soil_filename !< filename deep soil raw data
    CHARACTER (len=filename_max) :: soil_buffer_file  !< name for soil buffer file
    CHARACTER (len=filename_max) :: soil_output_file  !< name for soil output file
    CHARACTER (len=filename_max) :: soil_buffer_file_consistent !< name for soil buffer file after consistency check
    CHARACTER (len=filename_max) :: soil_output_file_consistent !< name for soil output file after consistency check


    !>Define the namelist group for soil raw data
    NAMELIST /soil_raw_data/ isoil_data, ldeep_soil, raw_data_soil_path, raw_data_soil_filename, raw_data_deep_soil_filename

    !> namelist with filenames for output of soil data
    NAMELIST /soil_io_extpar/ soil_buffer_file, soil_output_file, soil_buffer_file_consistent, soil_output_file_consistent

    INTEGER           :: nuin !< unit number
    INTEGER (KIND=i4) :: ierr !< error flag


    nuin = free_un()  ! functioin free_un returns free Fortran unit number
    OPEN(nuin,FILE=TRIM(namelist_file), IOSTAT=ierr)

    READ(nuin, NML=soil_raw_data, IOSTAT=ierr)
    READ(nuin, NML=soil_io_extpar, IOSTAT=ierr)

    CLOSE(nuin)
  
  END SUBROUTINE read_namelists_extpar_soil
!---------------------------------------------------------------------------

     !> inquire dimension information for NDVI raw data 
     SUBROUTINE get_dimension_soil_data(path_soil_file, &
                                          nlon_full, &
                                          nlat_full, &
                                          n_unit)

       !! here the coordintates, legend and data are read into global variables from the "Soil_data" Module
       USE mo_soil_data, ONLY: lon_full       !< longitide coordinates of the regular grid in the geographical (lonlat)
       USE mo_soil_data, ONLY: lat_full       !< latitude coordinates of the regular grid in the geographical (lonlat)

       USE mo_soil_data, ONLY: dsmw_grid      !< structure with the definition of the soil raw data grid

       CHARACTER (len=*), INTENT(IN) :: path_soil_file         !< filename with path for soil raw data
       INTEGER (KIND=i4), INTENT(OUT) :: nlon_full !< number of grid elements in zonal direction for soil data
       INTEGER (KIND=i4), INTENT(OUT) :: nlat_full !< number of grid elements in meridional direction for soil data
       INTEGER (KIND=i4), INTENT(OUT) :: n_unit    !< number of soil unit numbers (legend)

       !local variables
       INTEGER :: ncid                             !< netcdf unit file number
       INTEGER :: ndimension                       !< number of dimensions in netcdf file
       INTEGER :: nVars                            !< number of variables in netcdf file
       INTEGER :: nGlobalAtts                      !< number of gloabal Attributes in netcdf file
       INTEGER :: unlimdimid                       !< id of unlimited dimension (e.g. time) in netcdf file

       INTEGER :: dimid                            !< id of dimension
       CHARACTER (LEN=80) :: dimname               !< name of dimensiona
       INTEGER :: length                           !< length of dimension

       INTEGER :: varid                            !< id of variable

       CHARACTER (LEN=80) :: varname               !< name of variable
       INTEGER :: xtype                            !< netcdf type of variable/attribute
       INTEGER :: ndim                             !< number of dimensions of variable
       INTEGER, ALLOCATABLE :: var_dimids(:)       !< id of variable dimensions, vector, maximal dimension ndimension
       INTEGER :: nAtts                            !< number of attributes for a netcdf variable

       INTEGER :: errorcode                        !< error status variable

       ! open netcdf file 
       CALL check_netcdf(nf90_open(TRIM(path_soil_file),NF90_NOWRITE, ncid))

       ! look for numbers of dimensions, Variable, Attributes, and the dimid for the one possible unlimited dimension 
       ! (probably time)
       !; nf90_inquire input: ncid; nf90_inquire output: ndimension, nVars, nGlobalAtts,unlimdimid
       CALL check_netcdf (nf90_inquire(ncid,ndimension, nVars, nGlobalAtts,unlimdimid))
       !print *,'ncid,ndimension, nVars, nGlobalAtts,unlimdimid',ncid,ndimension, nVars, nGlobalAtts,unlimdimid
       ALLOCATE (var_dimids(ndimension), STAT=errorcode)


       !; the dimid in netcdf-files is counted from 1 to ndimension
       !; look for the name and length of the dimension with f90_inquire_dimension
       !; nf90_inquire_dimension input: ncid, dimid; nf90_inquire_dimension output: name, length
       DO dimid=1,ndimension
                    !print *,'dimension loop dimid ',dimid
         CALL check_netcdf( nf90_inquire_dimension(ncid,dimid, dimname, length) )
                     !print*, 'ncid,dimid, dimname, length',ncid,dimid, trim(dimname), length
         IF ( TRIM(dimname) == 'lon') nlon_full=length          ! here I know that the name of zonal dimension is 'lon'
         IF ( TRIM(dimname) == 'lat') nlat_full=length          ! here I know that the name of meridional dimension is 'lat'
         IF ( TRIM(dimname) == 'soil_unit') n_unit=length        ! here I know that the name of unit dimension is 'soil_unit'
       ENDDO

       ! Allocate lon_full, lat_full
       ALLOCATE (lon_full(nlon_full), lat_full(nlat_full))

       ! read lon_full and lat_full from raw data file
       variables: DO varid=1,nVars
         !print *,'variable loop, varid ',varid
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

       !HA debug
       PRINT *,'get_soil_data, dsmw_grid: ', dsmw_grid
       ! close netcdf file 
       CALL check_netcdf( nf90_close(ncid))

     END SUBROUTINE get_dimension_soil_data

!----------------------------------------------------------------------------------------------------------------       

     !> get coordintates, legend and data for soil raw data 
     SUBROUTINE get_soil_data(path_soil_file,start)
       !! here the coordintates, legend and data are read into global variables from the "Soil_data" Module
       USE mo_soil_data, ONLY: lon_soil       !< longitide coordinates of the regular grid in the geographical (lonlat)
       USE mo_soil_data, ONLY: lat_soil       !< latitude coordinates of the regular grid in the geographical (lonlat)

       USE mo_soil_data, ONLY: dsmw_grid      !< structure with the definition of the soil raw data grid
       USE mo_soil_data, ONLY: soil_texslo    !< legend for DSMW with texture and slope information
       USE mo_soil_data, ONLY: dsmw_soil_unit !< FAO Digital Soil Map of the World, the values represent the 
                                               !< soil unit number (see for legend in variable soil_texslo)
       !USE SOIL_data, ONLY: reg_lonlat_grid !< Definition of Data Type to describe a regular lonlat grid
       !USE SOIL_data, ONLY: dsmw_legend     !< Definition of Data Type to describe the legend for the FAO DSMW

       CHARACTER (len=*), INTENT(in) :: path_soil_file                !< filename with path for soil raw data
       INTEGER , INTENT(IN)          :: start(2)

       !local variables
       INTEGER :: ncid                             !< netcdf unit file number
       INTEGER :: ndimension                       !< number of dimensions in netcdf file
       INTEGER :: nVars                            !< number of variables in netcdf file
       INTEGER :: nGlobalAtts                      !< number of gloabal Attributes in netcdf file
       INTEGER :: unlimdimid                       !< id of unlimited dimension (e.g. time) in netcdf file

       INTEGER :: dimid                            !< id of dimension
       CHARACTER (LEN=80) :: dimname               !< name of dimension
       INTEGER :: length                           !< length of dimension

       INTEGER :: varid                            !< id of variable

       CHARACTER (LEN=80) :: varname               !< name of variable
       INTEGER :: xtype                            !< netcdf type of variable/attribute
       INTEGER :: ndim                             !< number of dimensions of variable
       INTEGER, ALLOCATABLE :: var_dimids(:)       !< id of variable dimensions, vector, maximal dimension ndimension
       INTEGER :: nAtts                            !< number of attributes for a netcdf variable
        
       INTEGER :: errorcode                        !< error status variable


       INTEGER (KIND=i4) :: n_unit    !< number of soil unit numbers (legend)

 

       ! open netcdf file
       CALL check_netcdf( nf90_open(TRIM(path_soil_file),NF90_NOWRITE, ncid))
       ! look for numbers of dimensions, Variable, Attributes, and the dimid for the unlimited dimension (probably time)
       ! nf90_inquire input: ncid; nf90_inquire output: ndimension, nVars, nGlobalAtts,unlimdimid
       CALL check_netcdf (nf90_inquire(ncid,ndimension, nVars, nGlobalAtts,unlimdimid))
       !print *,'ncid,ndimension, nVars, nGlobalAtts,unlimdimid',ncid,ndimension, nVars, nGlobalAtts,unlimdimid

       ALLOCATE (var_dimids(ndimension), STAT=errorcode)
       IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array var_dimids')
       var_dimids = 0
       !; the dimid in netcdf-files is counted from 1 to ndimension
       !; look for the name and length of the dimension with f90_inquire_dimension
       !; nf90_inquire_dimension input: ncid, dimid; nf90_inquire_dimension output: name, length
       DO dimid=1,ndimension
                    !print *,'dimension loop dimid ',dimid
         CALL check_netcdf( nf90_inquire_dimension(ncid,dimid, dimname, length) )
         !print*, 'ncid,dimid, dimname, length',ncid,dimid, trim(dimname), length
         IF ( TRIM(dimname) == 'soil_unit') n_unit=length        ! here I know that the name of unit dimension is 'soil_unit'
       ENDDO
       ! the deep soil has the same dimensions as the top soil

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
       
       !HA debug
       PRINT *,'get_soil_data, dsmw_grid: ', dsmw_grid



     END SUBROUTINE get_soil_data
     !------------------------------------------------------------------------------------------------------------

     !> get coordintates, legend and data for soil raw data 
     SUBROUTINE get_deep_soil_data(path_deep_soil_file,start)
     !! here the coordintates, legend and data are read into global variables from the "Soil_data" Module
     USE mo_soil_data, ONLY: dsmw_grid      !< structure with the definition of the soil raw data grid
     USE mo_soil_data, ONLY: soil_texslo_deep    !< legend for DSMW with texture and slope information
     USE mo_soil_data, ONLY: dsmw_deep_soil_unit !< FAO Digital Soil Map of the World, the values represent the 
                                                    !< soil unit number (see for legend in variable soil_texslo)
     !USE SOIL_data, ONLY: reg_lonlat_grid !< Definition of Data Type to describe a regular lonlat grid
     !USE SOIL_data, ONLY: dsmw_legend     !< Definition of Data Type to describe the legend for the FAO DSMW


     INTEGER , INTENT(IN)          :: start(2)
     CHARACTER (LEN=*), INTENT(IN) :: path_deep_soil_file !< filename with path for deep soil raw data

     !local variables
     INTEGER :: ncid_deep                        !< netcdf unit file number
     INTEGER :: ndimension_deep                       !< number of dimensions in netcdf file
     INTEGER :: nVars_deep                            !< number of variables in netcdf file
     INTEGER :: nGlobalAtts_deep                      !< number of gloabal Attributes in netcdf file
     INTEGER :: unlimdimid_deep                       !< id of unlimited dimension (e.g. time) in netcdf file

     INTEGER :: dimid_deep                            !< id of dimension
     CHARACTER (LEN=80) :: dimname               !< name of dimension
     INTEGER :: length                           !< length of dimension

     INTEGER :: varid_deep                            !< id of variable
     CHARACTER (LEN=80) :: varname_deep               !< name of variable
     INTEGER :: xtype_deep                            !< netcdf type of variable/attribute
     INTEGER :: ndim_deep                             !< number of dimensions of variable
     INTEGER, ALLOCATABLE :: var_dimids_deep(:)  !< id of variable dimensions, vector, maximal dimension ndimension

     INTEGER :: nAtts_deep                       !< number of attributes for a netcdf variable
        
     INTEGER :: errorcode                        !< error status variable

               
     CALL check_netcdf( nf90_open(TRIM(path_deep_soil_file),NF90_NOWRITE, ncid_deep))

     CALL check_netcdf (nf90_inquire(ncid_deep,ndimension_deep, nVars_deep, nGlobalAtts_deep,unlimdimid_deep))
        
     ALLOCATE (var_dimids_deep(ndimension_deep), STAT=errorcode)
     IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array var_dimids_deep')
     var_dimids_deep = 0
        
     variables_deep: DO varid_deep=1,nVars_deep
       !print *,'variable loop, varid ',varid
       CALL check_netcdf(nf90_inquire_variable(ncid_deep,varid_deep,varname_deep, & !_br 21.02.14 splitted too long line
                            xtype_deep, ndim_deep, var_dimids_deep, nAtts_deep))    !_br 21.02.14

       getvar_deep: SELECT CASE(TRIM(varname_deep))
            
       CASE('code')    !  here I know that the variable with the dsmw_code is called 'code'
         CALL check_netcdf(nf90_get_var(ncid_deep,varid_deep,soil_texslo_deep(:)%dsmw_code) )
     
       CASE('tex_coarse')!  here I know that the variable with the  coarse texture part is called 'tex_coarse'
         CALL check_netcdf(nf90_get_var(ncid_deep,varid_deep,soil_texslo_deep(:)%tex_coarse) )
            
       CASE('tex_medium')!  here I know that the variable with the medium texture part is called 'tex_medium'
         CALL check_netcdf(nf90_get_var(ncid_deep,varid_deep,soil_texslo_deep(:)%tex_medium) )
            
       CASE('tex_fine')!  here I know that the variable with the fine texture part is called 'tex_fine'
         CALL check_netcdf(nf90_get_var(ncid_deep,varid_deep,soil_texslo_deep(:)%tex_fine) )

       CASE('undefined')!  here I know that the variable with the undefined soil part is called 'undefined'
         CALL check_netcdf(nf90_get_var(ncid_deep,varid_deep,soil_texslo_deep(:)%part_undefined) )

       CASE('flat')!  here I know that the variable with the flat area part is called 'flat'
         CALL check_netcdf(nf90_get_var(ncid_deep,varid_deep,soil_texslo_deep(:)%flat) )

       CASE('hilly')!  here I know that the variable with the hilly area part is called 'hilly'
         CALL check_netcdf(nf90_get_var(ncid_deep,varid_deep,soil_texslo_deep(:)%hilly) )

       CASE('steep')!  here I know that the variable with the steep area part is called 'steep'
         CALL check_netcdf(nf90_get_var(ncid_deep,varid_deep,soil_texslo_deep(:)%steep) )
            
       CASE('Soil')  !  here I know that the variable with the DSMW soil units is called 'Soil'
         CALL check_netcdf(nf90_get_var(ncid_deep,varid_deep,dsmw_deep_soil_unit,start=start,count=(/ nlon_soil, nlat_soil /)) )
 
       END SELECT getvar_deep
          
     ENDDO variables_deep
        
     CALL check_netcdf( nf90_close( ncid_deep))

    END SUBROUTINE get_deep_soil_data


END MODULE mo_soil_routines
