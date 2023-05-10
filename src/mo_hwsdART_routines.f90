MODULE mo_hwsdART_routines

USE mo_logging
!> kind parameters are defined in MODULE data_parameters
USE mo_kind, ONLY: wp
USE mo_kind, ONLY: i4

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

USE mo_GRID_structures,  ONLY: reg_lonlat_grid


IMPLICIT NONE

PRIVATE

PUBLIC :: get_hwsdART_data
PUBLIC :: get_dimension_hwsdART_data
PUBLIC :: read_namelists_extpar_hwsdART

CONTAINS


!---------------------------------------------------------------------------
!> subroutine to read namelist for hwsdART data settings for EXTPAR 
SUBROUTINE read_namelists_extpar_hwsdART(namelist_file,                 &
                                         raw_data_hwsdART_path,         &
                                         raw_data_hwsdART_filename,     &
                                         hwsdART_buffer_file,           &
                                         hwsdART_output_file)

  USE mo_utilities_extpar, ONLY: free_un ! function to get free unit number

  
  CHARACTER (len=filename_max), INTENT(IN) :: namelist_file !< filename with namelists for for EXTPAR settings

 ! hwsdART
   CHARACTER (len=filename_max) :: raw_data_hwsdART_path          !< path to raw data
   CHARACTER (len=filename_max) :: raw_data_hwsdART_filename      !< filename hwsdART raw data
   CHARACTER (len=filename_max) :: hwsdART_buffer_file            !< name for hwsdART buffer file
   CHARACTER (len=filename_max) :: hwsdART_output_file            !< name for hwsdART output file


   !>Define the namelist group for hwsdART raw data
   NAMELIST /hwsdART_nml/ raw_data_hwsdART_path, raw_data_hwsdART_filename,hwsdART_buffer_file, hwsdART_output_file

   INTEGER           :: nuin !< unit number
   INTEGER (KIND=i4) :: ierr !< error flag


   nuin = free_un()  ! functioin free_un returns free Fortran unit number
   OPEN(nuin,FILE=TRIM(namelist_file), IOSTAT=ierr)

   READ(nuin, NML=hwsdART_nml, IOSTAT=ierr)

   CLOSE(nuin)
  

END SUBROUTINE read_namelists_extpar_hwsdART
!---------------------------------------------------------------------------

        !> inquire dimension information
        SUBROUTINE get_dimension_hwsdART_data(path_hwsdART_file, &
                                          nlon_hwsdART, &
                                          nlat_hwsdART)


        CHARACTER (len=*), INTENT(in) :: path_hwsdART_file         !< filename with path for hwsdART raw data
        INTEGER (KIND=i4), INTENT(out) :: nlon_hwsdART !< number of grid elements in zonal direction for hwsdART data
        INTEGER (KIND=i4), INTENT(out) :: nlat_hwsdART !< number of grid elements in meridional direction for hwsdART data

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
        call check_netcdf(nf90_open(TRIM(path_hwsdART_file),NF90_NOWRITE, ncid))

       ! look for numbers of dimensions, Variable, Attributes, and the dimid for the one possible unlimited dimension 
       ! (probably time)
       !; nf90_inquire input: ncid; nf90_inquire output: ndimension, nVars, nGlobalAtts,unlimdimid
       call check_netcdf (nf90_inquire(ncid,ndimension, nVars, nGlobalAtts,unlimdimid))



       !; the dimid in netcdf-files is counted from 1 to ndimension
       !; look for the name and length of the dimension with f90_inquire_dimension
       !; nf90_inquire_dimension input: ncid, dimid; nf90_inquire_dimension output: name, length
       do dimid=1,ndimension
         call check_netcdf( nf90_inquire_dimension(ncid,dimid, dimname, length) )

         if ( trim(dimname) == 'lon') nlon_hwsdART=length          ! here I know that the name of zonal dimension is 'lon'
         if ( trim(dimname) == 'lat') nlat_hwsdART=length          ! here I know that the name of meridional dimension is 'lat'
       enddo


       ! close netcdf file 
       call check_netcdf( nf90_close(ncid))

       END SUBROUTINE get_dimension_hwsdART_data

!----------------------------------------------------------------------------------------------------------------       


        !> get coordintates, legend and data for hwsdART raw data 
        SUBROUTINE get_hwsdART_data(path_hwsdART_file)
        !! here the coordintates, legend and data are read into global variables from the "hwsdART_data" Module
        USE mo_hwsdART_data, ONLY: lon_hwsdART       !< longitide coordinates of the regular grid in the geographical (lonlat)
        USE mo_hwsdART_data, ONLY: lat_hwsdART       !< latitude coordinates of the regular grid in the geographical (lonlat)

        USE mo_hwsdART_data, ONLY: hwsdART_grid      !< structure with the definition of the hwsdART raw data grid
        USE mo_hwsdART_data, ONLY: hwsdART_soil_unit !< The values represent the hwsdART unit number

        CHARACTER (len=*), INTENT(in) :: path_hwsdART_file                !< filename with path for hwsdART raw data





        !local variables
        INTEGER :: ncid                             !< netcdf unit file number
        INTEGER :: ndimension                       !< number of dimensions in netcdf file
        INTEGER :: nVars                            !< number of variables in netcdf file
        INTEGER :: nGlobalAtts                      !< number of gloabal Attributes in netcdf file
        INTEGER :: unlimdimid                       !< id of unlimited dimension (e.g. time) in netcdf file

        INTEGER :: dimid                            !< id of dimension
        CHARACTER (len=80) :: dimname               !< name of dimension
        INTEGER :: length                           !< length of dimension

        INTEGER :: varid                            !< id of variable

        CHARACTER (len=80) :: varname               !< name of variable
        INTEGER :: xtype                            !< netcdf type of variable/attribute
        INTEGER :: ndim                             !< number of dimensions of variable
        INTEGER, ALLOCATABLE :: var_dimids(:)       !< id of variable dimensions, vector, maximal dimension ndimension
        INTEGER :: nAtts                            !< number of attributes for a netcdf variable
        
        INTEGER :: errorcode                        !< error status variable


        INTEGER (KIND=i4) :: nlon_hwsdART !< number of grid elements in zonal direction for hwsdART data
        INTEGER (KIND=i4) :: nlat_hwsdART !< number of grid elements in meridional direction for hwsdART data

 

       ! open netcdf file
        call check_netcdf( nf90_open(TRIM(path_hwsdART_file),NF90_NOWRITE, ncid))
       ! look for numbers of dimensions, Variable, Attributes, and the dimid for the one possible unlimited dimension (probably time)
       ! nf90_inquire input: ncid; nf90_inquire output: ndimension, nVars, nGlobalAtts,unlimdimid
       call check_netcdf (nf90_inquire(ncid,ndimension, nVars, nGlobalAtts,unlimdimid))


        ALLOCATE (var_dimids(ndimension), STAT=errorcode)
          IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array var_dimids')
           var_dimids = 0
       !; the dimid in netcdf-files is counted from 1 to ndimension
       !; look for the name and length of the dimension with f90_inquire_dimension
       !; nf90_inquire_dimension input: ncid, dimid; nf90_inquire_dimension output: name, length
       do dimid=1,ndimension

         call check_netcdf( nf90_inquire_dimension(ncid,dimid, dimname, length) )

         if ( trim(dimname) == 'lon') nlon_hwsdART=length          ! here I know that the name of zonal dimension is 'lon'
         if ( trim(dimname) == 'lat') nlat_hwsdART=length          ! here I know that the name of meridional dimension is 'lat'
       enddo
       ! the deep hwsdART has the same dimensions as the top hwsdART



                !; the varid in netcdf files is counted from 1 to nVars
                !; look for variable names, number of dimension, var_dimids etc with nf90_inquire_variable
                !; nf90_inquire_variable input: ncid, varid; nf90_inquire_variable output name xtype, ndim, var_dimids, nAtts


         variables: DO varid=1,nVars
           CALL check_netcdf(nf90_inquire_variable(ncid,varid,varname,xtype, ndim, var_dimids, nAtts))
           
           getvar: SELECT CASE(TRIM(varname))
             
           CASE('lon')     !  here I know that the variable with the longitude coordinates is called 'lon'
             CALL check_netcdf(nf90_get_var(ncid,varid,lon_hwsdART(:)) )
             
           CASE('lat')     !  here I know that the variable with the latitude coordinates is called 'lat'
             CALL check_netcdf(nf90_get_var(ncid,varid,lat_hwsdART(:)) )
             
           CASE('LU')    !  here I know that the variable with the land use index is called 'LU'
             CALL check_netcdf(nf90_get_var(ncid,varid,hwsdART_soil_unit(:,:)) )
                
           END SELECT getvar
           
         ENDDO variables
         
         ! close netcdf file 
         call check_netcdf( nf90_close( ncid))
         
                 
       ! Fill the structure hwsdART_raw_data_grid

       hwsdART_grid%nlon_reg = nlon_hwsdART
       hwsdART_grid%nlat_reg = nlat_hwsdART

       hwsdART_grid%start_lon_reg = lon_hwsdART(1)
       hwsdART_grid%end_lon_reg   = lon_hwsdART(nlon_hwsdART)

       hwsdART_grid%start_lat_reg = lat_hwsdART(1)
       hwsdART_grid%end_lat_reg   = lat_hwsdART(nlat_hwsdART)
       
       IF (hwsdART_grid%nlon_reg /= 0) THEN
       hwsdART_grid%dlon_reg = (hwsdART_grid%end_lon_reg - hwsdART_grid%start_lon_reg)/(hwsdART_grid%nlon_reg-1._wp)
       ENDIF 

       IF (hwsdART_grid%nlat_reg /= 0) THEN
       hwsdART_grid%dlat_reg = (hwsdART_grid%end_lat_reg - hwsdART_grid%start_lat_reg)/(hwsdART_grid%nlat_reg-1._wp)
       ENDIF ! in case of latitude orientation from north to south dlat is negative!

     WRITE(message_text,*) 'get_hwsdART_data, hwsdART_grid: ', hwsdART_grid
     CALL logging%info(message_text)  

     END SUBROUTINE get_hwsdART_data
     !------------------------------------------------------------------------------------------------------------

END MODULE mo_hwsdART_routines
