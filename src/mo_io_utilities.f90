!+ Fortran module for input/ouptut utilities
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_0         2010/12/21 Hermann Asensio
!  Initial release
! V1_2         2011/03/25 Hermann Asensio
!  Update documentation
!  use capital letters in variable names for netcdf files
!  enhance keylen_max too 100 characters
!  add netcdf-attribute missing_value
! V1_3         2011/04/19 Hermann Asensio
!  set time dimension for netcdf output as "unlimited"
! V1_4         2011/04/21 Hermann Asensio
!  clean up
! V1_7         2013/01/25 Guenther Zaengl
!   Parallel threads for ICON and COSMO using Open-MP,
!   Several bug fixes and optimizations for ICON search algorithm,
!   particularly for the special case of non-contiguous domains;
!   simplified namelist control for ICON
! V2_0         2013/08/08 Daniel Luethi
!   added possibility to add name of data set to var_meta_info
! V2_0_3       2015-01-12 Juergen Helmert
!   added code for number_of_grid_used global attribute
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!> Fortran module for input/ouptut utilities
!! some useful functions and data structures for netcdf i/o
MODULE mo_io_utilities

  !> kind parameters are defined in MODULE data_parameters
  USE mo_kind, ONLY: wp
  USE mo_kind, ONLY: i8
  USE mo_kind, ONLY: i4

  USE netcdf

  USE mo_utilities_extpar, ONLY: abort_extpar

  IMPLICIT NONE

  PUBLIC :: dim_meta_info
  PUBLIC :: var_meta_info

  PUBLIC :: netcdf_attributes
  PUBLIC :: netcdf_char_attributes
  PUBLIC :: netcdf_real_attributes
  PUBLIC :: netcdf_int_attributes

  PUBLIC :: netcdf_grid_mapping

  PUBLIC :: check_netcdf

  PUBLIC :: netcdf_def_grid_mapping

  PUBLIC :: netcdf_put_var
  PUBLIC :: netcdf_put_int_i4_1d
  PUBLIC :: netcdf_put_int_i4_2d
  PUBLIC :: netcdf_put_int_i4_3d
  PUBLIC :: netcdf_put_int_i4_4d
  PUBLIC :: netcdf_put_int_i4_5d
  PUBLIC :: netcdf_put_int_i8_1d
  PUBLIC :: netcdf_put_int_i8_2d
  PUBLIC :: netcdf_put_int_i8_3d
  PUBLIC :: netcdf_put_int_i8_4d
  PUBLIC :: netcdf_put_int_i8_5d
  PUBLIC :: netcdf_put_real_1d
  PUBLIC :: netcdf_put_real_2d
  PUBLIC :: netcdf_put_real_3d
  PUBLIC :: netcdf_put_real_4d
  PUBLIC :: netcdf_put_real_5d
  PUBLIC :: open_new_netcdf_file
  PUBLIC :: close_netcdf_file

  PUBLIC :: netcdf_get_var
  PUBLIC :: netcdf_get_var_real_3d, netcdf_get_var_real_4d, netcdf_get_var_real_5d
  PUBLIC :: netcdf_get_var_int_3d_i4 ,netcdf_get_var_int_3d_i8, netcdf_get_var_int_4d

  PUBLIC :: vartype_int
  PUBLIC :: vartype_real
  PUBLIC :: vartype_char

  PUBLIC :: get_date_const_field
  PUBLIC :: set_date_mm_extpar_field

  INTEGER, PARAMETER :: keylen_max = 100 !< maximum length for the length of keys of type character

  !> structure to save dimension information
  TYPE dim_meta_info
    CHARACTER (len=12) :: dimname !< name of dimension
    INTEGER            :: dimsize !< size of dimension
    INTEGER            :: dimid   !< netcdf id for dimension
  END TYPE dim_meta_info

  !> structure to save meta information for target variables (name, units, etc)
  TYPE var_meta_info
    CHARACTER (len=20)           :: varname       !< name of variable
    INTEGER                      :: n_dim         !< number of dimensions
    TYPE(dim_meta_info), POINTER :: diminfo(:)    !< pointer to dimensions of variable
    INTEGER                      :: vartype       !< type of variable, 1 for INTEGER, 2 for REAL, 3 for CHARACTER
    !INTEGER                     :: varid         !< netcdf varid of variable
    CHARACTER (len=keylen_max)   :: standard_name !< netcdf attribute for standard name
    CHARACTER (len=keylen_max)   :: long_name     !< netcdf attribute for long name
    CHARACTER (len=20)           :: units         !< netcdf attribute for units
    CHARACTER (len=keylen_max)   :: grid_mapping  !< netcdf attribute grid mapping
    CHARACTER (len=keylen_max)   :: coordinates   !< netcdf attribute coordinates
    CHARACTER (len=keylen_max)   :: shortName     !< GRIB API shortName key
    CHARACTER (len=10)           :: stepType      !< GRIB API type of temporal information (edition independant)
    CHARACTER (len=keylen_max)   :: data_set      !< name of source data set
  END TYPE var_meta_info

  INTEGER, PARAMETER :: vartype_int = 1
  INTEGER, PARAMETER :: vartype_real = 2
  INTEGER, PARAMETER :: vartype_char = 3

  !> structure to save global netcdf attributes
  TYPE netcdf_attributes
    CHARACTER (len=80)  :: attname
    CHARACTER (len=255) :: attributetext
  END TYPE netcdf_attributes

  !> structure to store character type netcdf attributes
  TYPE netcdf_char_attributes
    CHARACTER (len=80)  :: attname
    character (len=255) :: attributetext
  END TYPE netcdf_char_attributes

  !> structure to store real type netcdf attributes
  TYPE netcdf_real_attributes
    CHARACTER (len=80)  :: attname
    REAL                :: att_value_r
  END TYPE netcdf_real_attributes

  !> structure to store integer type netcdf attributes
  TYPE netcdf_int_attributes
    CHARACTER (len=80)  :: attname
    INTEGER             :: att_value_int
  END TYPE netcdf_int_attributes

  !> structure to store netcdf grid mapping information according to cf conventions
  TYPE netcdf_grid_mapping
    CHARACTER (len=80)                        :: grid_mapping_varname !< name for variable in netcdf file with grid_mapping data
    TYPE(netcdf_char_attributes)              :: grid_mapping_name    !< netcdf attribute with grid mapping name according to cf,
                                                                      !< see e.g.
                                                                      !< http://cf-pcmdi.llnl.gov/documents/cf-conventions/1.4/apf.html
    INTEGER                                   :: n_r_att              !< number of real attributes
    TYPE(netcdf_real_attributes), ALLOCATABLE :: map_param(:)         !< store mapping parameters for the netcdf output
  END TYPE netcdf_grid_mapping

  !> put attributes to netcdf file
  INTERFACE netcdf_put_att
    MODULE PROCEDURE netcdf_put_att_real
    MODULE PROCEDURE netcdf_put_att_int
  END INTERFACE

  !> put variables to netcdf file
  INTERFACE netcdf_put_var
    MODULE PROCEDURE netcdf_put_int_scalar
    MODULE PROCEDURE netcdf_put_int_i4_1d
    MODULE PROCEDURE netcdf_put_int_i4_2d
    MODULE PROCEDURE netcdf_put_int_i4_3d
    MODULE PROCEDURE netcdf_put_int_i4_4d
    MODULE PROCEDURE netcdf_put_int_i4_5d
    MODULE PROCEDURE netcdf_put_int_i8_1d
    MODULE PROCEDURE netcdf_put_int_i8_2d
    MODULE PROCEDURE netcdf_put_int_i8_3d
    MODULE PROCEDURE netcdf_put_int_i8_4d
    MODULE PROCEDURE netcdf_put_int_i8_5d
    MODULE PROCEDURE netcdf_put_real_scalar
    MODULE PROCEDURE netcdf_put_real_1d
    MODULE PROCEDURE netcdf_put_real_2d
    MODULE PROCEDURE netcdf_put_real_3d
    MODULE PROCEDURE netcdf_put_real_4d
    MODULE PROCEDURE netcdf_put_real_5d
  END INTERFACE

  !> get netcdf variables
  INTERFACE netcdf_get_var
    MODULE PROCEDURE netcdf_get_var_int_3d_i8
    MODULE PROCEDURE netcdf_get_var_int_3d_i4
    MODULE PROCEDURE netcdf_get_var_int_4d
    MODULE PROCEDURE netcdf_get_var_real_2d
    MODULE PROCEDURE netcdf_get_var_real_3d
    MODULE PROCEDURE netcdf_get_var_real_4d
    MODULE PROCEDURE netcdf_get_var_real_5d
  END INTERFACE

CONTAINS

  !> subroutine to check error status of netcdf function error codes
  SUBROUTINE check_netcdf(status, filename, line_number)
    INTEGER,          INTENT(in)           :: status
    CHARACTER(len=*), INTENT(in), OPTIONAL :: filename
    INTEGER,          INTENT(in), OPTIONAL :: line_number

    CHARACTER(len=132) abort_message

    IF (status /= nf90_noerr) then
      IF (PRESENT(filename) .AND. PRESENT(line_number)) THEN
        WRITE(abort_message,'(a,i4,a)') trim(filename), line_number, ': '//TRIM(nf90_strerror(status))
      ELSE
        WRITE(abort_message,'(a)') TRIM(nf90_strerror(status))
      ENDIF
      CALL abort_extpar(abort_message)
    END IF

  END SUBROUTINE check_netcdf

  !-----------------------------------------------------------------------------

  !> specific subroutine to put some standard attributes to an integer type netcdf variable
  SUBROUTINE netcdf_put_att_int(ncid, varid, varinfo, fill_value_i)

    INTEGER, INTENT(IN)             :: ncid         !< id for netcdf file
    INTEGER, INTENT(IN)             :: varid        !< variable id for netcdf variable
    TYPE(var_meta_info), INTENT(IN) :: varinfo      !< structure with information on the meta data
    INTEGER, INTENT(IN)             :: fill_value_i !< undefined value

    ! local variables
    INTEGER :: fill_value
    INTEGER :: missing_value

    CALL check_netcdf( nf90_put_att(ncid,     &
         &             varid,                 &
         &             TRIM('standard_name'), &
         &             TRIM(varinfo%standard_name)), __FILE__, __LINE__ )

    CALL check_netcdf( nf90_put_att(ncid, &
         &             varid,             &
         &             TRIM('long_name'), &
         &             TRIM(varinfo%long_name)), __FILE__, __LINE__ )

    IF (TRIM(varinfo%units) /= "-") THEN
      CALL check_netcdf( nf90_put_att(ncid, &
           &             varid,             &
           &             TRIM('units'),     &
           &             TRIM(varinfo%units)), __FILE__, __LINE__ )
    ENDIF

    IF (TRIM(varinfo%grid_mapping) /= "-") THEN
      CALL check_netcdf( nf90_put_att(ncid,    &
           &             varid,                &
           &             TRIM('grid_mapping'), &
           &             TRIM(varinfo%grid_mapping)), __FILE__, __LINE__ )


    ENDIF

    IF (TRIM(varinfo%coordinates) /= "-") THEN
      CALL check_netcdf( nf90_put_att(ncid,   &
           &             varid,               &
           &             TRIM('coordinates'), &
           &             TRIM(varinfo%coordinates)), __FILE__, __LINE__ )
    ENDIF

    fill_value = fill_value_i
    CALL check_netcdf( nf90_put_att(ncid,  &
         &             varid,              &
         &             TRIM('_FillValue'), &
         &             fill_value), __FILE__, __LINE__ )

    missing_value  = fill_value_i !type conversion
    CALL check_netcdf( nf90_put_att(ncid,     &
         &             varid,                 &
         &             TRIM('missing_value'), &
         &             missing_value), __FILE__, __LINE__ )

    IF (TRIM(varinfo%data_set) /= "-") THEN
      CALL check_netcdf( nf90_put_att(ncid, &
           &             varid,             &
           &             TRIM('data_set'),  &
           &             TRIM(varinfo%data_set)), __FILE__, __LINE__ )
    ENDIF

  END SUBROUTINE netcdf_put_att_int

  !> specific subroutine to put some standard attributes to an real type netcdf variable
  SUBROUTINE netcdf_put_att_real(ncid, varid, varinfo, fill_value_r)
    INTEGER, INTENT(IN)             :: ncid         !< id for netcdf file
    INTEGER, INTENT(IN)             :: varid        !< variable id for netcdf variable
    TYPE(var_meta_info), INTENT(IN) :: varinfo      !< structure with information on the meta data
    REAL(KIND=wp), INTENT(IN)       :: fill_value_r !< undefined value

    ! local variable
    REAL :: fill_value
    REAL :: missing_value

    CALL check_netcdf( nf90_put_att(ncid,     &
         &             varid,                 &
         &             TRIM('standard_name'), &
         &             TRIM(varinfo%standard_name)), __FILE__, __LINE__ )

    CALL check_netcdf( nf90_put_att(ncid, &
         &             varid,             &
         &             TRIM('long_name'), &
         &             TRIM(varinfo%long_name)), __FILE__, __LINE__ )

    IF (TRIM(varinfo%units) /= "-") THEN
      CALL check_netcdf( nf90_put_att(ncid, &
           &             varid,             &
           &             TRIM('units'),     &
           &             TRIM(varinfo%units)), __FILE__, __LINE__ )
    ENDIF

    IF (TRIM(varinfo%grid_mapping) /= "-") THEN
      CALL check_netcdf( nf90_put_att(ncid,    &
           &             varid,                &
           &             TRIM('grid_mapping'), &
           &             TRIM(varinfo%grid_mapping)), __FILE__, __LINE__ )
    ENDIF

    IF (TRIM(varinfo%coordinates) /= "-") THEN
      CALL check_netcdf( nf90_put_att(ncid, &
           &             varid, &
           &             TRIM('coordinates'), &
           &             TRIM(varinfo%coordinates)), __FILE__, __LINE__ )
    ENDIF

    fill_value = fill_value_r !type conversion
    CALL check_netcdf( nf90_put_att(ncid,  &
         &             varid,              &
         &             TRIM('_FillValue'), &
         &             fill_value), __FILE__, __LINE__ )

    missing_value  = fill_value_r !type conversion
    CALL check_netcdf( nf90_put_att(ncid,     &
         &             varid,                 &
         &             TRIM('missing_value'), &
         &             missing_value), __FILE__, __LINE__ )

    IF (TRIM(varinfo%data_set) /= "-") THEN
      CALL check_netcdf( nf90_put_att(ncid, &
           &             varid,             &
           &             TRIM('data_set'),  &
           &             TRIM(varinfo%data_set)), __FILE__, __LINE__ )
    ENDIF


  END SUBROUTINE netcdf_put_att_real

  !-----------------------------------------------------------------------------

  !> specific subroutine to define scalar real variable for netcdf
  SUBROUTINE netcdf_put_real_scalar(ncid, var_real, meta, fill_value_r)

    INTEGER, INTENT(IN) :: ncid !< id for netcdf file
    REAL (KIND=wp), INTENT(IN) ::  var_real
    TYPE (var_meta_info), INTENT(IN) :: meta   !< addtional information for array
    REAL(KIND=wp), INTENT(IN) :: fill_value_r

    !local variables
    INTEGER :: varid  !< netcdf varid of variable
    CHARACTER (len=20) :: varname    !< name of variable
    ! redef netcdf file
    CALL check_netcdf(nf90_redef(ncid), __FILE__, __LINE__ )

    ! get varname
    varname = TRIM(meta%varname)
    ! define netcdf variable
    CALL check_netcdf( nf90_def_var(ncid,       &
         &                          varname,    &
         &                          NF90_FLOAT, &
         &                          varid), __FILE__, __LINE__ )
    ! put standard attributes to variable
    CALL netcdf_put_att(ncid, varid, meta, fill_value_r)

    ! end of definition
    CALL check_netcdf(nf90_enddef(ncid), __FILE__, __LINE__ )

    ! put variable to netcdf file

    CALL check_netcdf(nf90_put_var(ncid,varid,var_real), __FILE__, __LINE__ )

  END SUBROUTINE netcdf_put_real_scalar

  !-----------------------------------------------------------------------------

  !> specific subroutine to define 1d real variable for netcdf
  SUBROUTINE netcdf_put_real_1d(ncid, var_real_1d, meta_1d, fill_value_r)

    INTEGER, INTENT(IN)              :: ncid            !< id for netcdf file
    REAL (KIND=wp), INTENT(IN)       :: var_real_1d(:)
    TYPE (var_meta_info), INTENT(IN) :: meta_1d         !< addtional information for array
    REAL(KIND=wp), INTENT(IN)        :: fill_value_r

    !local variables
    INTEGER            :: dimid_1d(1) !< dimension ids for 1d variable
    INTEGER            :: varid       !< netcdf varid of variable
    INTEGER            :: errorcode
    CHARACTER (len=20) :: varname     !< name of variable

    ! redef netcdf file
    CALL check_netcdf(nf90_redef(ncid), __FILE__, __LINE__ )

    ! get dimid
    errorcode = nf90_inq_dimid(ncid, &
         &                     TRIM(meta_1d%diminfo(1)%dimname), &
         &                     dimid_1d(1))

    IF (errorcode /= nf90_noerr) THEN
      CALL check_netcdf( nf90_def_dim(ncid,                            &
           &                          TRIM(meta_1d%diminfo(1)%dimname),&
           &                          meta_1d%diminfo(1)%dimsize,      &
           &                          dimid_1d(1)))
    ENDIF

    ! get varname
    varname = TRIM(meta_1d%varname)

    ! define netcdf variable
    CALL check_netcdf( nf90_def_var(ncid,       &
         &                          varname,    &
         &                          NF90_FLOAT, &
         &                          dimid_1d,   &
         &                          varid), __FILE__, __LINE__ )

    ! put standard attributes to variable
    CALL netcdf_put_att(ncid, varid, meta_1d, fill_value_r)

    ! end of definition
    CALL check_netcdf(nf90_enddef(ncid), __FILE__, __LINE__ )

    ! put variable to netcdf file

    CALL check_netcdf(nf90_put_var(ncid,varid,var_real_1d), __FILE__, __LINE__ )

  END SUBROUTINE netcdf_put_real_1d

  !-----------------------------------------------------------------------------

  !> specific subroutine to define 2d real variable for netcdf
  SUBROUTINE netcdf_put_real_2d(ncid, var_real_2d, meta_2d, fill_value_r)

    INTEGER, INTENT(IN)              :: ncid              !< id for netcdf file
    REAL (KIND=wp), INTENT(IN)       ::  var_real_2d(:,:)
    TYPE (var_meta_info), INTENT(IN) :: meta_2d           !< addtional information for array
    REAL(KIND=wp), INTENT(IN)        :: fill_value_r

    !local variables
    INTEGER            :: dimid_2d(2) !< dimension ids for 2d variable
    INTEGER            :: varid       !< netcdf varid of variable
    INTEGER            :: errorcode
    CHARACTER (len=20) :: varname     !< name of variable

    ! redef netcdf file
    CALL check_netcdf(nf90_redef(ncid))

    ! get dimid
    errorcode = nf90_inq_dimid( ncid,                             &
         &                      TRIM(meta_2d%diminfo(1)%dimname), &
         &                      dimid_2d(1))

    IF (errorcode /= nf90_noerr) THEN
      CALL check_netcdf( nf90_def_dim( ncid,                             &
           &                           TRIM(meta_2d%diminfo(1)%dimname), &
           &                           meta_2d%diminfo(1)%dimsize,       &
           &                           dimid_2d(1)), __FILE__, __LINE__ )
    ENDIF

    errorcode = nf90_inq_dimid( ncid, &
         &                      TRIM(meta_2d%diminfo(2)%dimname), &
         &                      dimid_2d(2))

    IF (errorcode /= nf90_noerr) THEN
      CALL check_netcdf( nf90_def_dim( ncid,                             &
           &                           TRIM(meta_2d%diminfo(2)%dimname), &
           &                           meta_2d%diminfo(2)%dimsize,       &
           &                           dimid_2d(2)), __FILE__, __LINE__ )
    ENDIF

    ! get varname
    varname = TRIM(meta_2d%varname)
    ! define netcdf variable
    CALL check_netcdf( nf90_def_var(ncid,       &
         &                          varname,    &
         &                          NF90_FLOAT, &
         &                          dimid_2d,   &
                                    varid), __FILE__, __LINE__ )
    ! put standard attributes to variable
    CALL netcdf_put_att(ncid, varid, meta_2d, fill_value_r)

    ! end of definition
    CALL check_netcdf(nf90_enddef(ncid), __FILE__, __LINE__ )

    ! put variable to netcdf file

    CALL check_netcdf(nf90_put_var(ncid,varid,var_real_2d), __FILE__, __LINE__ )

  END SUBROUTINE netcdf_put_real_2d

  !-----------------------------------------------------------------------------

  !> specific subroutine to define 3d real variable for netcdf
  SUBROUTINE netcdf_put_real_3d(ncid, var_real_3d, meta_3d, fill_value_r)

    INTEGER, INTENT(IN)              :: ncid               !< id for netcdf file
    REAL (KIND=wp), INTENT(IN)       :: var_real_3d(:,:,:)
    TYPE (var_meta_info), INTENT(IN) :: meta_3d            !< addtional information for array
    REAL(KIND=wp), INTENT(IN)        :: fill_value_r

    !local variables
    INTEGER            :: dimid_3d(3)  !< dimension ids for 3d variable
    INTEGER            :: varid        !< netcdf varid of variable
    INTEGER            :: errorcode
    CHARACTER (len=20) :: varname      !< name of variable

    ! redef netcdf file
    CALL check_netcdf(nf90_redef(ncid))

    ! get dimid
    errorcode = nf90_inq_dimid( ncid,                             &
         &                      TRIM(meta_3d%diminfo(1)%dimname), &
         &                      dimid_3d(1) )

    IF (errorcode /= nf90_noerr) THEN
      CALL check_netcdf( nf90_def_dim( ncid,                             &
           &                           TRIM(meta_3d%diminfo(1)%dimname), &
           &                           meta_3d%diminfo(1)%dimsize,       &
           &                           dimid_3d(1)), __FILE__, __LINE__ )
    ENDIF

    errorcode = nf90_inq_dimid( ncid,                             &
         &                      TRIM(meta_3d%diminfo(2)%dimname), &
         &                      dimid_3d(2))
    IF (errorcode /= nf90_noerr) THEN
      CALL check_netcdf( nf90_def_dim( ncid,                             &
           &                           TRIM(meta_3d%diminfo(2)%dimname), &
           &                           meta_3d%diminfo(2)%dimsize,       &
           &                           dimid_3d(2)), __FILE__, __LINE__ )
    ENDIF

    errorcode = nf90_inq_dimid( ncid,                             &
         &                      TRIM(meta_3d%diminfo(3)%dimname), &
         &                      dimid_3d(3) )

    IF (errorcode /= nf90_noerr) THEN
      CALL check_netcdf( nf90_def_dim( ncid,                             &
           &                           TRIM(meta_3d%diminfo(3)%dimname), &
           &                           meta_3d%diminfo(3)%dimsize,       &
           &                           dimid_3d(3)), __FILE__, __LINE__ )
    ENDIF

    ! get varname
    varname = TRIM(meta_3d%varname)

    ! define netcdf variable
    CALL check_netcdf( nf90_def_var(ncid,       &
         &                          varname,    &
         &                          NF90_FLOAT, &
         &                          dimid_3d,   &
                                    varid), __FILE__, __LINE__ )

    ! put standard attributes to variable
    CALL netcdf_put_att(ncid, varid, meta_3d, fill_value_r)

    ! end of definition
    CALL check_netcdf(nf90_enddef(ncid), __FILE__, __LINE__ )

    ! put variable to netcdf file

    CALL check_netcdf(nf90_put_var(ncid,varid,var_real_3d), __FILE__, __LINE__ )

  END SUBROUTINE netcdf_put_real_3d

  !-----------------------------------------------------------------------------

  !> specific subroutine to define 4d real variable for netcdf
  SUBROUTINE netcdf_put_real_4d(ncid, var_real_4d, meta_4d, fill_value_r)

    INTEGER, INTENT(IN)              :: ncid                  !< id for netcdf file
    REAL (KIND=wp), INTENT(IN)       ::  var_real_4d(:,:,:,:)
    TYPE (var_meta_info), INTENT(IN) :: meta_4d               !< addtional information for array
    REAL(KIND=wp), INTENT(IN)        :: fill_value_r

    !local variables
    INTEGER            :: n
    INTEGER            :: dimid_4d(4) !< dimension ids for 4D variable
    INTEGER            :: varid       !< netcdf varid of variable
    INTEGER            :: errorcode
    CHARACTER (len=20) :: varname     !< name of variable

    ! redef netcdf file
    CALL check_netcdf(nf90_redef(ncid), __FILE__, __LINE__ )

    ! get dimid
    errorcode = nf90_inq_dimid( ncid,                             &
         &                      TRIM(meta_4d%diminfo(1)%dimname), &
         &                      dimid_4d(1))

    IF (errorcode /= nf90_noerr) THEN
      CALL check_netcdf( nf90_def_dim( ncid,                             &
           &                           TRIM(meta_4d%diminfo(1)%dimname), &
           &                           meta_4d%diminfo(1)%dimsize,       &
           &                           dimid_4d(1)), __FILE__, __LINE__ )
    ENDIF

    errorcode = nf90_inq_dimid(ncid, &
         &                        TRIM(meta_4d%diminfo(2)%dimname), &
         &                        dimid_4d(2))

    IF (errorcode /= nf90_noerr) THEN
      CALL check_netcdf( nf90_def_dim( ncid,                             &
           &                           TRIM(meta_4d%diminfo(2)%dimname), &
           &                           meta_4d%diminfo(2)%dimsize,       &
           &                           dimid_4d(2)), __FILE__, __LINE__ )
    ENDIF

    errorcode = nf90_inq_dimid( ncid,                             &
         &                      TRIM(meta_4d%diminfo(3)%dimname), &
         &                      dimid_4d(3))

    IF (errorcode /= nf90_noerr) THEN
      CALL check_netcdf( nf90_def_dim( ncid,                             &
           &                           TRIM(meta_4d%diminfo(3)%dimname), &
           &                           meta_4d%diminfo(3)%dimsize,       &
           &                           dimid_4d(3)), __FILE__, __LINE__ )
    ENDIF

    errorcode = nf90_inq_dimid( ncid,                             &
         &                      TRIM(meta_4d%diminfo(4)%dimname), &
         &                      dimid_4d(4))

    IF (errorcode /= nf90_noerr) THEN
      CALL check_netcdf( nf90_def_dim( ncid,                             &
           &                           TRIM(meta_4d%diminfo(4)%dimname), &
           &                           meta_4d%diminfo(4)%dimsize,       &
           &                           dimid_4d(4)), __FILE__, __LINE__ )
    ENDIF

    ! get varname
    varname = TRIM(meta_4d%varname)

    ! define netcdf variable
    CALL check_netcdf( nf90_def_var(ncid,       &
         &                          varname,    &
         &                          NF90_FLOAT, &
         &                          dimid_4d,   &
                                    varid), __FILE__, __LINE__ )

    ! put standard attributes to variable
    CALL netcdf_put_att(ncid, varid, meta_4d, fill_value_r)

    ! end of definition
    CALL check_netcdf(nf90_enddef(ncid), __FILE__, __LINE__ )

    ! put variable to netcdf file

    !    CALL check_netcdf(nf90_put_var(ncid,varid,var_real_4d))
    DO n=1,size(var_real_4d,4)
      CALL check_netcdf(nf90_put_var(ncid,varid,var_real_4d(:,:,:,n),start=(/1,1,1,n/)), &
         &              __FILE__, __LINE__ )
    ENDDO

  END SUBROUTINE netcdf_put_real_4d

  !-----------------------------------------------------------------------------

  !> specific subroutine to define 5d real variable for netcdf
  SUBROUTINE netcdf_put_real_5d(ncid, var_real_5d, meta_5d, fill_value_r)

    INTEGER, INTENT(IN)              :: ncid                   !< id for netcdf file
    REAL (KIND=wp), INTENT(IN)       :: var_real_5d(:,:,:,:,:)
    TYPE (var_meta_info), INTENT(IN) :: meta_5d                !< addtional information for array
    REAL(KIND=wp), INTENT(IN)        :: fill_value_r

    !local variables
    INTEGER            :: dimid_5d(5)  !< dimension ids for 5d variable
    INTEGER            :: varid        !< netcdf varid of variable
    INTEGER            :: errorcode
    CHARACTER (len=20) :: varname      !< name of variable

    ! redef netcdf file
    CALL check_netcdf(nf90_redef(ncid), __FILE__, __LINE__ )

    ! get dimid
    errorcode = nf90_inq_dimid( ncid,                             &
         &                      TRIM(meta_5d%diminfo(1)%dimname), &
         &                      dimid_5d(1))

    IF (errorcode /= nf90_noerr) THEN
      CALL check_netcdf( nf90_def_dim( ncid,                             &
           &                           TRIM(meta_5d%diminfo(1)%dimname), &
           &                           meta_5d%diminfo(1)%dimsize,       &
           &                           dimid_5d(1)), __FILE__, __LINE__ )
    ENDIF

    errorcode = nf90_inq_dimid( ncid,                             &
         &                      TRIM(meta_5d%diminfo(2)%dimname), &
         &                      dimid_5d(2))
    IF (errorcode /= nf90_noerr) THEN
      CALL check_netcdf( nf90_def_dim( ncid,                             &
           &                           TRIM(meta_5d%diminfo(2)%dimname), &
           &                           meta_5d%diminfo(2)%dimsize,       &
           &                           dimid_5d(2)), __FILE__, __LINE__ )
    ENDIF

    errorcode = nf90_inq_dimid( ncid,                             &
         &                      TRIM(meta_5d%diminfo(3)%dimname), &
         &                      dimid_5d(3))
    IF (errorcode /= nf90_noerr) THEN
      CALL check_netcdf( nf90_def_dim( ncid,                             &
           &                           TRIM(meta_5d%diminfo(3)%dimname), &
           &                           meta_5d%diminfo(3)%dimsize,       &
           &                           dimid_5d(3)), __FILE__, __LINE__ )
    ENDIF

    errorcode = nf90_inq_dimid( ncid,                             &
         &                      TRIM(meta_5d%diminfo(4)%dimname), &
         &                      dimid_5d(4))

    IF (errorcode /= nf90_noerr) THEN
      CALL check_netcdf( nf90_def_dim( ncid,                             &
           &                           TRIM(meta_5d%diminfo(4)%dimname), &
           &                           meta_5d%diminfo(4)%dimsize,       &
           &                           dimid_5d(4)), __FILE__, __LINE__ )
    ENDIF

    errorcode = nf90_inq_dimid( ncid,                             &
         &                      TRIM(meta_5d%diminfo(5)%dimname), &
         &                      dimid_5d(5))

    IF (errorcode /= nf90_noerr) THEN
      CALL check_netcdf( nf90_def_dim( ncid,                             &
           &                           TRIM(meta_5d%diminfo(5)%dimname), &
           &                           meta_5d%diminfo(5)%dimsize,       &
           &                           dimid_5d(5)), __FILE__, __LINE__ )
    ENDIF

    ! get varname
    varname = TRIM(meta_5d%varname)

    ! define netcdf variable
    PRINT *,'put_real_5d: ',varname
    CALL check_netcdf( nf90_def_var(ncid,       &
         &                          varname,    &
         &                          NF90_FLOAT, &
         &                          dimid_5d,   &
                                    varid), __FILE__, __LINE__ )

    ! put standard attributes to variable
    CALL netcdf_put_att(ncid, varid, meta_5d, fill_value_r)

    ! end of definition
    CALL check_netcdf(nf90_enddef(ncid), __FILE__, __LINE__ )

    ! put variable to netcdf file
    CALL check_netcdf(nf90_put_var(ncid,varid,var_real_5d), __FILE__, __LINE__ )

  END SUBROUTINE netcdf_put_real_5d

  !-----------------------------------------------------------------------------

  !> specific subroutine to define scalar int variable for netcdf
  SUBROUTINE netcdf_put_int_scalar(ncid, var_int, meta, fill_value_i)

    INTEGER, INTENT(IN)              :: ncid         !< id for netcdf file
    INTEGER (KIND=i8)                :: var_int
    TYPE (var_meta_info), INTENT(IN) :: meta         !< addtional information for array
    INTEGER, INTENT(IN)              :: fill_value_i

    !local variables
    INTEGER            :: varid    !< netcdf varid of variable
    CHARACTER (len=20) :: varname  !< name of variable

    ! redef netcdf file
    CALL check_netcdf(nf90_redef(ncid), __FILE__, __LINE__ )

    ! get varname
    varname = TRIM(meta%varname)

    ! define netcdf variable
    CALL check_netcdf( nf90_def_var(ncid,     &
         &                          varname,  &
         &                          NF90_INT, &
                                    varid), __FILE__, __LINE__ )
    ! put standard attributes to variable
    CALL netcdf_put_att(ncid, varid, meta, fill_value_i)

    ! end of definition
    CALL check_netcdf(nf90_enddef(ncid), __FILE__, __LINE__ )

    ! put variable to netcdf file

    CALL check_netcdf(nf90_put_var(ncid,varid,var_int), __FILE__, __LINE__ )

  END SUBROUTINE netcdf_put_int_scalar

  !-----------------------------------------------------------------------------

  !> specific subroutine to define 1d int variable for netcdf
  SUBROUTINE netcdf_put_int_i8_1d(ncid, var_int_1d, meta_1d, fill_value_i)

    INTEGER, INTENT(IN) :: ncid !< id for netcdf file
    INTEGER (KIND=i8), INTENT(IN) ::  var_int_1d(:)
    TYPE (var_meta_info), INTENT(IN) :: meta_1d   !< addtional information for array
    INTEGER (KIND=i8), INTENT(IN) :: fill_value_i

    !local variables
    INTEGER            :: dimid_1d(1) !< dimension ids for 1d variable
    INTEGER            :: varid       !< netcdf varid of variable
    INTEGER            :: errorcode
    INTEGER            :: fill_value
    CHARACTER (len=20) :: varname     !< name of variable

    ! redef netcdf file
    CALL check_netcdf(nf90_redef(ncid), __FILE__, __LINE__ )

    ! get dimid
    errorcode = nf90_inq_dimid( ncid,                             &
         &                      TRIM(meta_1d%diminfo(1)%dimname), &
         &                      dimid_1d(1))

    IF (errorcode /= nf90_noerr) THEN
      CALL check_netcdf( nf90_def_dim( ncid,                             &
           &                           TRIM(meta_1d%diminfo(1)%dimname), &
           &                           meta_1d%diminfo(1)%dimsize,       &
           &                           dimid_1d(1)), __FILE__, __LINE__ )
    ENDIF

    ! get varname
    varname = TRIM(meta_1d%varname)

    ! define netcdf variable
    CALL check_netcdf( nf90_def_var(ncid,     &
         &                          varname,  &
         &                          NF90_INT, &
         &                          dimid_1d, &
                                    varid), __FILE__, __LINE__ )

    ! put standard attributes to variable
    fill_value = fill_value_i
    CALL netcdf_put_att(ncid, varid, meta_1d, fill_value)

    ! end of definition
    CALL check_netcdf(nf90_enddef(ncid), __FILE__, __LINE__ )

    ! put variable to netcdf file
    CALL check_netcdf(nf90_put_var(ncid,varid,var_int_1d), __FILE__, __LINE__ )

  END SUBROUTINE netcdf_put_int_i8_1d

  !-----------------------------------------------------------------------------

  !> specific subroutine to define 2d int variable for netcdf
  SUBROUTINE netcdf_put_int_i8_2d(ncid, var_int_2d, meta_2d, fill_value_i)

    INTEGER, INTENT(IN)              :: ncid !< id for netcdf file
    INTEGER (KIND=i8), INTENT(IN)    ::  var_int_2d(:,:)
    TYPE (var_meta_info), INTENT(IN) :: meta_2d   !< addtional information for array
    INTEGER (KIND=i8), INTENT(IN)    :: fill_value_i

    !local variables
    INTEGER            :: dimid_2d(2) !< dimension ids for 2d variable
    INTEGER            :: varid       !< netcdf varid of variable
    INTEGER            :: errorcode
    INTEGER            :: fill_value
    CHARACTER (len=20) :: varname     !< name of variable

    ! redef netcdf file
    CALL check_netcdf(nf90_redef(ncid), __FILE__, __LINE__ )

    ! get dimid
    errorcode = nf90_inq_dimid( ncid,                             &
         &                      TRIM(meta_2d%diminfo(1)%dimname), &
         &                      dimid_2d(1))

    IF (errorcode /= nf90_noerr) THEN
      CALL check_netcdf( nf90_def_dim( ncid,                             &
           &                           TRIM(meta_2d%diminfo(1)%dimname), &
           &                           meta_2d%diminfo(1)%dimsize,       &
           &                           dimid_2d(1)), __FILE__, __LINE__ )
    ENDIF

    errorcode = nf90_inq_dimid( ncid,                             &
         &                      TRIM(meta_2d%diminfo(2)%dimname), &
         &                      dimid_2d(2))

    IF (errorcode /= nf90_noerr) THEN
      CALL check_netcdf( nf90_def_dim( ncid,                             &
           &                           TRIM(meta_2d%diminfo(2)%dimname), &
           &                           meta_2d%diminfo(2)%dimsize,       &
           &                           dimid_2d(2)), __FILE__, __LINE__ )
    ENDIF

    ! get varname
    varname = TRIM(meta_2d%varname)

    ! define netcdf variable
    CALL check_netcdf( nf90_def_var(ncid, &
         & varname,            &
         & NF90_INT,           &
         & dimid_2d,           &
         varid), __FILE__, __LINE__ )

    ! put standard attributes to variable
    fill_value = fill_value_i
    CALL netcdf_put_att(ncid, varid, meta_2d, fill_value)

    ! end of definition
    CALL check_netcdf(nf90_enddef(ncid), __FILE__, __LINE__ )

    ! put variable to netcdf file
    CALL check_netcdf(nf90_put_var(ncid,varid,var_int_2d), __FILE__, __LINE__ )

  END SUBROUTINE netcdf_put_int_i8_2d

  !-----------------------------------------------------------------------------

  !> specific subroutine to define 3d int variable for netcdf
  SUBROUTINE netcdf_put_int_i8_3d(ncid, var_int_3d, meta_3d, fill_value_i)

    INTEGER, INTENT(IN)              :: ncid              !< id for netcdf file
    INTEGER (KIND=i8), INTENT(IN)    :: var_int_3d(:,:,:)
    TYPE (var_meta_info), INTENT(IN) :: meta_3d           !< addtional information for array
    INTEGER (KIND=i8), INTENT(IN)    :: fill_value_i

    !local variables
    INTEGER            :: dimid_3d(3) !< dimension ids for 3d variable
    INTEGER            :: varid       !< netcdf varid of variable
    INTEGER            :: errorcode
    INTEGER            :: fill_value
    CHARACTER (len=20) :: varname     !< name of variable

    ! redef netcdf file
    CALL check_netcdf(nf90_redef(ncid), __FILE__, __LINE__ )

    ! get dimid
    errorcode = nf90_inq_dimid( ncid,                             &
         &                      TRIM(meta_3d%diminfo(1)%dimname), &
         &                      dimid_3d(1))

    IF (errorcode /= nf90_noerr) THEN
      CALL check_netcdf( nf90_def_dim( ncid,                             &
           &                           TRIM(meta_3d%diminfo(1)%dimname), &
           &                           meta_3d%diminfo(1)%dimsize,       &
           &                           dimid_3d(1)), __FILE__, __LINE__ )
    ENDIF

    errorcode = nf90_inq_dimid( ncid,                             &
         &                      TRIM(meta_3d%diminfo(2)%dimname), &
         &                      dimid_3d(2))
    IF (errorcode /= nf90_noerr) THEN
      CALL check_netcdf( nf90_def_dim( ncid,                             &
           &                           TRIM(meta_3d%diminfo(2)%dimname), &
           &                           meta_3d%diminfo(2)%dimsize,       &
           &                           dimid_3d(2)), __FILE__, __LINE__ )
    ENDIF

    errorcode = nf90_inq_dimid( ncid,                             &
         &                      TRIM(meta_3d%diminfo(3)%dimname), &
         &                      dimid_3d(3))

    IF (errorcode /= nf90_noerr) THEN
      CALL check_netcdf( nf90_def_dim( ncid,                             &
           &                           TRIM(meta_3d%diminfo(3)%dimname), &
           &                           meta_3d%diminfo(3)%dimsize,       &
           &                           dimid_3d(3)), __FILE__, __LINE__ )
    ENDIF

    ! get varname
    varname = TRIM(meta_3d%varname)

    ! define netcdf variable
    CALL check_netcdf( nf90_def_var(ncid,     &
         &                          varname,  &
         &                          NF90_INT, &
         &                          dimid_3d, &
                                    varid), __FILE__, __LINE__ )

    ! put standard attributes to variable
    fill_value = fill_value_i
    CALL netcdf_put_att(ncid, varid, meta_3d, fill_value)

    ! end of definition
    CALL check_netcdf(nf90_enddef(ncid), __FILE__, __LINE__ )

    ! put variable to netcdf file

    CALL check_netcdf(nf90_put_var(ncid,varid,var_int_3d), __FILE__, __LINE__ )

  END SUBROUTINE netcdf_put_int_i8_3d

  !-----------------------------------------------------------------------------

  !> specific subroutine to define 4d int variable for netcdf
  SUBROUTINE netcdf_put_int_i8_4d(ncid, var_int_4d, meta_4d, fill_value_i)

    INTEGER, INTENT(IN)              :: ncid !< id for netcdf file
    INTEGER (KIND=i8), INTENT(IN)    ::  var_int_4d(:,:,:,:)
    TYPE (var_meta_info), INTENT(IN) :: meta_4d   !< addtional information for array
    INTEGER (KIND=i8), INTENT(IN)    :: fill_value_i

    !local variables
    INTEGER            :: dimid_4d(4)  !< dimension ids for 4D variable
    INTEGER            :: varid        !< netcdf varid of variable
    INTEGER            :: errorcode
    INTEGER            :: fill_value
    CHARACTER (len=20) :: varname      !< name of variable

    ! redef netcdf file
    CALL check_netcdf(nf90_redef(ncid), __FILE__, __LINE__ )

    ! get dimid
    errorcode = nf90_inq_dimid( ncid,                             &
         &                      TRIM(meta_4d%diminfo(1)%dimname), &
         &                      dimid_4d(1))

    IF (errorcode /= nf90_noerr) THEN
      CALL check_netcdf( nf90_def_dim( ncid,                             &
           &                           TRIM(meta_4d%diminfo(1)%dimname), &
           &                           meta_4d%diminfo(1)%dimsize,       &
           &                           dimid_4d(1)), __FILE__, __LINE__ )
    ENDIF

    errorcode = nf90_inq_dimid( ncid,                             &
         &                      TRIM(meta_4d%diminfo(2)%dimname), &
         &                      dimid_4d(2))

    IF (errorcode /= nf90_noerr) THEN
      CALL check_netcdf( nf90_def_dim( ncid,                             &
           &                           TRIM(meta_4d%diminfo(2)%dimname), &
           &                           meta_4d%diminfo(2)%dimsize,       &
           &                           dimid_4d(2)), __FILE__, __LINE__ )
    ENDIF

    errorcode = nf90_inq_dimid( ncid,                             &
         &                      TRIM(meta_4d%diminfo(3)%dimname), &
         &                      dimid_4d(3))

    IF (errorcode /= nf90_noerr) THEN
      CALL check_netcdf( nf90_def_dim( ncid,                             &
           &                           TRIM(meta_4d%diminfo(3)%dimname), &
           &                           meta_4d%diminfo(3)%dimsize,       &
           &                           dimid_4d(3)), __FILE__, __LINE__ )
    ENDIF

    errorcode = nf90_inq_dimid( ncid,                             &
         &                      TRIM(meta_4d%diminfo(4)%dimname), &
         &                      dimid_4d(4))

    IF (errorcode /= nf90_noerr) THEN
      CALL check_netcdf( nf90_def_dim( ncid,                             &
           &                           TRIM(meta_4d%diminfo(4)%dimname), &
           &                           meta_4d%diminfo(4)%dimsize,       &
           &                           dimid_4d(4)), __FILE__, __LINE__ )
    ENDIF

    ! get varname
    varname = TRIM(meta_4d%varname)

    ! define netcdf variable
    CALL check_netcdf( nf90_def_var( ncid,     &
         &                           varname,  &
         &                           NF90_INT, &
         &                           dimid_4d, &
                                     varid), __FILE__, __LINE__ )

    ! put standard attributes to variable
    fill_value = fill_value_i
    CALL netcdf_put_att( ncid, varid, meta_4d, fill_value)

    ! end of definition
    CALL check_netcdf(nf90_enddef( ncid), __FILE__, __LINE__ )

    ! put variable to netcdf file

    CALL check_netcdf(nf90_put_var( ncid,varid,var_int_4d), __FILE__, __LINE__ )

  END SUBROUTINE netcdf_put_int_i8_4d

  !-----------------------------------------------------------------------------

  !> specific subroutine to define 5d int variable for netcdf
  SUBROUTINE netcdf_put_int_i8_5d(ncid, var_int_5d, meta_5d, fill_value_i)

    INTEGER, INTENT(IN)              :: ncid                   !< id for netcdf file
    INTEGER (KIND=i8), INTENT(IN)    ::  var_int_5d(:,:,:,:,:)
    TYPE (var_meta_info), INTENT(IN) :: meta_5d                !< addtional information for array
    INTEGER (KIND=i8), INTENT(IN)    :: fill_value_i

    !local variables
    INTEGER            :: dimid_5d(5)  !< dimension ids for 5d variable
    INTEGER            :: varid        !< netcdf varid of variable
    INTEGER            :: errorcode
    INTEGER            :: fill_value
    CHARACTER (len=20) :: varname      !< name of variable

    ! redef netcdf file
    CALL check_netcdf(nf90_redef(ncid))

    ! get dimid
    errorcode = nf90_inq_dimid( ncid,                             &
         &                      TRIM(meta_5d%diminfo(1)%dimname), &
         &                      dimid_5d(1))

    IF (errorcode /= nf90_noerr) THEN
      CALL check_netcdf( nf90_def_dim( ncid,                             &
           &                           TRIM(meta_5d%diminfo(1)%dimname), &
           &                           meta_5d%diminfo(1)%dimsize,       &
           &                           dimid_5d(1)), __FILE__, __LINE__ )
    ENDIF

    errorcode = nf90_inq_dimid( ncid,                            &
         &                     TRIM(meta_5d%diminfo(2)%dimname), &
         &                      dimid_5d(2))
    IF (errorcode /= nf90_noerr) THEN
      CALL check_netcdf( nf90_def_dim( ncid,                             &
           &                           TRIM(meta_5d%diminfo(2)%dimname), &
           &                           meta_5d%diminfo(2)%dimsize,       &
           &                           dimid_5d(2)), __FILE__, __LINE__ )
    ENDIF

    errorcode = nf90_inq_dimid( ncid,                             &
         &                      TRIM(meta_5d%diminfo(3)%dimname), &
         &                      dimid_5d(3))

    IF (errorcode /= nf90_noerr) THEN
      CALL check_netcdf( nf90_def_dim( ncid,                             &
           &                           TRIM(meta_5d%diminfo(3)%dimname), &
           &                           meta_5d%diminfo(3)%dimsize,       &
           &                           dimid_5d(3)), __FILE__, __LINE__ )
    ENDIF

    errorcode = nf90_inq_dimid( ncid,                             &
         &                      TRIM(meta_5d%diminfo(4)%dimname), &
         &                      dimid_5d(4))

    IF (errorcode /= nf90_noerr) THEN
      CALL check_netcdf( nf90_def_dim( ncid,                             &
           &                           TRIM(meta_5d%diminfo(4)%dimname), &
           &                           meta_5d%diminfo(4)%dimsize,       &
           &                           dimid_5d(4)), __FILE__, __LINE__ )
    ENDIF

    errorcode = nf90_inq_dimid( ncid,                             &
         &                      TRIM(meta_5d%diminfo(5)%dimname), &
         &                      dimid_5d(5))

    IF (errorcode /= nf90_noerr) THEN
      CALL check_netcdf( nf90_def_dim( ncid,                             &
           &                           TRIM(meta_5d%diminfo(5)%dimname), &
           &                           meta_5d%diminfo(5)%dimsize,       &
           &                           dimid_5d(5)), __FILE__, __LINE__ )
    ENDIF

    ! get varname
    varname = TRIM(meta_5d%varname)

    ! define netcdf variable
    CALL check_netcdf( nf90_def_var( ncid,     &
         &                           varname,  &
         &                           NF90_INT, &
         &                           dimid_5d, &
                                     varid), __FILE__, __LINE__ )

    ! put standard attributes to variable
    fill_value = fill_value_i
    CALL netcdf_put_att(ncid, varid, meta_5d, fill_value)

    ! end of definition
    CALL check_netcdf(nf90_enddef(ncid), __FILE__, __LINE__ )

    ! put variable to netcdf file

    CALL check_netcdf(nf90_put_var(ncid,varid,var_int_5d), __FILE__, __LINE__ )

  END SUBROUTINE netcdf_put_int_i8_5d

  !-----------------------------------------------------------------------------

  !> specific subroutine to define 1d int variable for netcdf
  SUBROUTINE netcdf_put_int_i4_1d(ncid, var_int_1d, meta_1d, fill_value_i)

    INTEGER, INTENT(IN)              :: ncid           !< id for netcdf file
    INTEGER (KIND=i4), INTENT(IN)    ::  var_int_1d(:)
    TYPE (var_meta_info), INTENT(IN) :: meta_1d        !< addtional information for array
    INTEGER, INTENT(IN)              :: fill_value_i

    !local variables
    INTEGER            :: dimid_1d(1)  !< dimension ids for 1d variable
    INTEGER            :: varid        !< netcdf varid of variable
    INTEGER            :: errorcode
    INTEGER            :: fill_value
    CHARACTER (len=20) :: varname      !< name of variable

    ! redef netcdf file
    CALL check_netcdf(nf90_redef(ncid), __FILE__, __LINE__ )

    ! get dimid
    errorcode = nf90_inq_dimid( ncid,                             &
         &                      TRIM(meta_1d%diminfo(1)%dimname), &
         &                      dimid_1d(1))

    IF (errorcode /= nf90_noerr) THEN
      CALL check_netcdf( nf90_def_dim( ncid,                             &
           &                           TRIM(meta_1d%diminfo(1)%dimname), &
           &                           meta_1d%diminfo(1)%dimsize,       &
           &                           dimid_1d(1)), __FILE__, __LINE__ )
    ENDIF

    ! get varname
    varname = TRIM(meta_1d%varname)

    ! define netcdf variable
    CALL check_netcdf( nf90_def_var( ncid,     &
         &                           varname,  &
         &                           NF90_INT, &
         &                           dimid_1d, &
                                     varid), __FILE__, __LINE__ )
    ! put standard attributes to variable
    fill_value = fill_value_i
    CALL netcdf_put_att(ncid, varid, meta_1d, fill_value)

    ! end of definition
    CALL check_netcdf(nf90_enddef(ncid), __FILE__, __LINE__ )

    ! put variable to netcdf file

    CALL check_netcdf(nf90_put_var(ncid,varid,var_int_1d), __FILE__, __LINE__ )

  END SUBROUTINE netcdf_put_int_i4_1d

  !-----------------------------------------------------------------------------

  !> specific subroutine to define 2d int variable for netcdf
  SUBROUTINE netcdf_put_int_i4_2d(ncid, var_int_2d, meta_2d, fill_value_i)

    INTEGER, INTENT(IN)              :: ncid             !< id for netcdf file
    INTEGER (KIND=i4), INTENT(IN)    ::  var_int_2d(:,:)
    TYPE (var_meta_info), INTENT(IN) :: meta_2d          !< addtional information for array
    INTEGER (KIND=i4), INTENT(IN)    :: fill_value_i

    !local variables
    INTEGER            :: dimid_2d(2) !< dimension ids for 2d variable
    INTEGER            :: varid       !< netcdf varid of variable
    INTEGER            :: errorcode
    INTEGER            :: fill_value
    CHARACTER (len=20) :: varname     !< name of variable

    ! redef netcdf file
    CALL check_netcdf(nf90_redef(ncid))

    ! get dimid
    errorcode = nf90_inq_dimid( ncid,                             &
         &                      TRIM(meta_2d%diminfo(1)%dimname), &
         &                      dimid_2d(1))
    IF (errorcode /= nf90_noerr) THEN
      CALL check_netcdf( nf90_def_dim( ncid,                             &
           &                           TRIM(meta_2d%diminfo(1)%dimname), &
           &                           meta_2d%diminfo(1)%dimsize,       &
           &                           dimid_2d(1)), __FILE__, __LINE__ )
    ENDIF

    errorcode = nf90_inq_dimid( ncid,                             &
         &                      TRIM(meta_2d%diminfo(2)%dimname), &
         &                      dimid_2d(2))
    IF (errorcode /= nf90_noerr) THEN
      CALL check_netcdf( nf90_def_dim( ncid,                             &
           &                           TRIM(meta_2d%diminfo(2)%dimname), &
           &                           meta_2d%diminfo(2)%dimsize,       &
           &                           dimid_2d(2)), __FILE__, __LINE__ )
    ENDIF

    ! get varname
    varname = TRIM(meta_2d%varname)

    ! define netcdf variable
    CALL check_netcdf( nf90_def_var( ncid,     &
         &                           varname,  &
         &                           NF90_INT, &
         &                           dimid_2d, &
                                     varid), __FILE__, __LINE__ )

    ! put standard attributes to variable
    fill_value = fill_value_i
    CALL netcdf_put_att(ncid, varid, meta_2d, fill_value)

    ! end of definition
    CALL check_netcdf(nf90_enddef(ncid), __FILE__, __LINE__ )

    ! put variable to netcdf file

    CALL check_netcdf(nf90_put_var(ncid,varid,var_int_2d), __FILE__, __LINE__ )

  END SUBROUTINE netcdf_put_int_i4_2d

  !-----------------------------------------------------------------------------

  !> specific subroutine to define 3d int variable for netcdf
  SUBROUTINE netcdf_put_int_i4_3d(ncid, var_int_3d, meta_3d, fill_value_i)

    INTEGER, INTENT(IN)              :: ncid              !< id for netcdf file
    INTEGER (KIND=i4), INTENT(IN)    :: var_int_3d(:,:,:)
    TYPE (var_meta_info), INTENT(IN) :: meta_3d           !< addtional information for array
    INTEGER (KIND=i4), INTENT(IN)    :: fill_value_i

    !local variables
    INTEGER            :: dimid_3d(3) !< dimension ids for 3d variable
    INTEGER            :: varid       !< netcdf varid of variable
    INTEGER            :: errorcode
    INTEGER            :: fill_value
    CHARACTER (len=20) :: varname     !< name of variable

    ! redef netcdf file
    CALL check_netcdf(nf90_redef(ncid), __FILE__, __LINE__ )

    ! get dimid
    errorcode = nf90_inq_dimid( ncid,                             &
         &                      TRIM(meta_3d%diminfo(1)%dimname), &
         &                      dimid_3d(1))
    IF (errorcode /= nf90_noerr) THEN
      CALL check_netcdf( nf90_def_dim( ncid,                             &
           &                           TRIM(meta_3d%diminfo(1)%dimname), &
           &                           meta_3d%diminfo(1)%dimsize,       &
           &                           dimid_3d(1)), __FILE__, __LINE__ )
    ENDIF

    errorcode = nf90_inq_dimid( ncid,                             &
         &                      TRIM(meta_3d%diminfo(2)%dimname), &
         &                      dimid_3d(2))

    IF (errorcode /= nf90_noerr) THEN
      CALL check_netcdf( nf90_def_dim( ncid,                             &
           &                           TRIM(meta_3d%diminfo(2)%dimname), &
           &                           meta_3d%diminfo(2)%dimsize,       &
           &                           dimid_3d(2)), __FILE__, __LINE__ )
    ENDIF

    errorcode = nf90_inq_dimid( ncid,                             &
         &                      TRIM(meta_3d%diminfo(3)%dimname), &
         &                      dimid_3d(3))

    IF (errorcode /= nf90_noerr) THEN
      CALL check_netcdf( nf90_def_dim( ncid,                             &
           &                           TRIM(meta_3d%diminfo(3)%dimname), &
           &                           meta_3d%diminfo(3)%dimsize,       &
           &                           dimid_3d(3)), __FILE__, __LINE__ )
    ENDIF

    ! get varname
    varname = TRIM(meta_3d%varname)

    ! define netcdf variable
    CALL check_netcdf( nf90_def_var( ncid,     &
         &                           varname,  &
         &                           NF90_INT, &
         &                           dimid_3d, &
                                     varid), __FILE__, __LINE__ )

    ! put standard attributes to variable
    fill_value = fill_value_i
    CALL netcdf_put_att(ncid, varid, meta_3d, fill_value)

    ! end of definition
    CALL check_netcdf(nf90_enddef(ncid), __FILE__, __LINE__ )

    ! put variable to netcdf file

    CALL check_netcdf(nf90_put_var(ncid,varid,var_int_3d), __FILE__, __LINE__ )

  END SUBROUTINE netcdf_put_int_i4_3d

  !-----------------------------------------------------------------------------

  !> specific subroutine to define 4d int variable for netcdf
  SUBROUTINE netcdf_put_int_i4_4d(ncid, var_int_4d, meta_4d, fill_value_i)

    INTEGER, INTENT(IN)              :: ncid                 !< id for netcdf file
    INTEGER (KIND=i4), INTENT(IN)    ::  var_int_4d(:,:,:,:)
    TYPE (var_meta_info), INTENT(IN) :: meta_4d              !< addtional information for array
    INTEGER (KIND=i4), INTENT(IN)    :: fill_value_i

    !local variables

    INTEGER            :: dimid_4d(4) !< dimension ids for 4D variable
    INTEGER            :: varid       !< netcdf varid of variable
    INTEGER            :: errorcode
    INTEGER            :: fill_value
    CHARACTER (len=20) :: varname     !< name of variable

    ! redef netcdf file
    CALL check_netcdf(nf90_redef(ncid), __FILE__, __LINE__ )

    ! get dimid
    errorcode = nf90_inq_dimid( ncid,                             &
         &                      TRIM(meta_4d%diminfo(1)%dimname), &
         &                      dimid_4d(1))

    IF (errorcode /= nf90_noerr) THEN
      CALL check_netcdf( nf90_def_dim( ncid,                             &
           &                           TRIM(meta_4d%diminfo(1)%dimname), &
           &                           meta_4d%diminfo(1)%dimsize,       &
           &                           dimid_4d(1)), __FILE__, __LINE__ )
    ENDIF

    errorcode = nf90_inq_dimid( ncid,                             &
         &                      TRIM(meta_4d%diminfo(2)%dimname), &
         &                      dimid_4d(2))

    IF (errorcode /= nf90_noerr) THEN
      CALL check_netcdf( nf90_def_dim( ncid,                             &
           &                           TRIM(meta_4d%diminfo(2)%dimname), &
           &                           meta_4d%diminfo(2)%dimsize,       &
           &                           dimid_4d(2)), __FILE__, __LINE__ )
    ENDIF

    errorcode = nf90_inq_dimid( ncid,                             &
         &                      TRIM(meta_4d%diminfo(3)%dimname), &
         &                      dimid_4d(3))

    IF (errorcode /= nf90_noerr) THEN
      CALL check_netcdf( nf90_def_dim( ncid,                             &
           &                           TRIM(meta_4d%diminfo(3)%dimname), &
           &                           meta_4d%diminfo(3)%dimsize,       &
           &                           dimid_4d(3)), __FILE__, __LINE__ )
    ENDIF

    errorcode = nf90_inq_dimid( ncid,                             &
         &                      TRIM(meta_4d%diminfo(4)%dimname), &
         &                      dimid_4d(4))

    IF (errorcode /= nf90_noerr) THEN
      CALL check_netcdf( nf90_def_dim( ncid,                             &
           &                           TRIM(meta_4d%diminfo(4)%dimname), &
           &                           meta_4d%diminfo(4)%dimsize,       &
           &                           dimid_4d(4)), __FILE__, __LINE__ )
    ENDIF

    ! get varname
    varname = TRIM(meta_4d%varname)

    ! define netcdf variable
    CALL check_netcdf( nf90_def_var( ncid,     &
         &                           varname,  &
         &                           NF90_INT, &
         &                           dimid_4d, &
                                     varid), __FILE__, __LINE__ )

    ! put standard attributes to variable
    fill_value = fill_value_i
    CALL netcdf_put_att(ncid, varid, meta_4d, fill_value)

    ! end of definition
    CALL check_netcdf(nf90_enddef(ncid), __FILE__, __LINE__ )

    ! put variable to netcdf file

    CALL check_netcdf(nf90_put_var(ncid,varid,var_int_4d), __FILE__, __LINE__ )

  END SUBROUTINE netcdf_put_int_i4_4d

  !-----------------------------------------------------------------------------

  !> specific subroutine to define 5d int variable for netcdf
  SUBROUTINE netcdf_put_int_i4_5d(ncid, var_int_5d, meta_5d, fill_value_i)

    INTEGER, INTENT(IN)              :: ncid !< id for netcdf file
    INTEGER (KIND=i4), INTENT(IN)    ::  var_int_5d(:,:,:,:,:)
    TYPE (var_meta_info), INTENT(IN) :: meta_5d   !< addtional information for array
    INTEGER (KIND=i4), INTENT(IN)    :: fill_value_i

    !local variables
    INTEGER            :: dimid_5d(5) !< dimension ids for 5d variable
    INTEGER            :: varid       !< netcdf varid of variable
    INTEGER            :: errorcode
    INTEGER            :: fill_value
    CHARACTER (len=20) :: varname     !< name of variable

    ! redef netcdf file
    CALL check_netcdf(nf90_redef(ncid), __FILE__, __LINE__ )

    ! get dimid
    errorcode = nf90_inq_dimid( ncid,                             &
         &                      TRIM(meta_5d%diminfo(1)%dimname), &
         &                      dimid_5d(1) )

    IF (errorcode /= nf90_noerr) THEN
      CALL check_netcdf( nf90_def_dim( ncid,                             &
           &                           TRIM(meta_5d%diminfo(1)%dimname), &
           &                           meta_5d%diminfo(1)%dimsize,       &
           &                           dimid_5d(1)), __FILE__, __LINE__ )
    ENDIF

    errorcode = nf90_inq_dimid( ncid,                             &
         &                      TRIM(meta_5d%diminfo(2)%dimname), &
         &                      dimid_5d(2))

    IF (errorcode /= nf90_noerr) THEN
      CALL check_netcdf( nf90_def_dim( ncid,                             &
           &                           TRIM(meta_5d%diminfo(2)%dimname), &
           &                           meta_5d%diminfo(2)%dimsize,       &
           &                           dimid_5d(2)), __FILE__, __LINE__ )
    ENDIF

    errorcode = nf90_inq_dimid( ncid,                             &
         &                      TRIM(meta_5d%diminfo(3)%dimname), &
         &                      dimid_5d(3))
    IF (errorcode /= nf90_noerr) THEN
      CALL check_netcdf( nf90_def_dim( ncid,                             &
           &                           TRIM(meta_5d%diminfo(3)%dimname), &
           &                           meta_5d%diminfo(3)%dimsize,       &
           &                           dimid_5d(3)), __FILE__, __LINE__ )
    ENDIF

    errorcode = nf90_inq_dimid( ncid,                             &
         &                      TRIM(meta_5d%diminfo(4)%dimname), &
         &                      dimid_5d(4))
    IF (errorcode /= nf90_noerr) THEN
      CALL check_netcdf( nf90_def_dim( ncid,                             &
           &                           TRIM(meta_5d%diminfo(4)%dimname), &
           &                           meta_5d%diminfo(4)%dimsize,       &
           &                           dimid_5d(4)), __FILE__, __LINE__ )
    ENDIF

    errorcode = nf90_inq_dimid( ncid,                             &
         &                      TRIM(meta_5d%diminfo(5)%dimname), &
         &                      dimid_5d(5))
    IF (errorcode /= nf90_noerr) THEN
      CALL check_netcdf( nf90_def_dim( ncid,                             &
           &                           TRIM(meta_5d%diminfo(5)%dimname), &
           &                           meta_5d%diminfo(5)%dimsize,       &
           &                           dimid_5d(5)), __FILE__, __LINE__ )
    ENDIF

    ! get varname
    varname = TRIM(meta_5d%varname)
    ! define netcdf variable
    CALL check_netcdf( nf90_def_var( ncid,     &
         &                           varname,  &
         &                           NF90_INT, &
         &                           dimid_5d, &
                                     varid), __FILE__, __LINE__ )
    ! put standard attributes to variable
    fill_value = fill_value_i
    CALL netcdf_put_att(ncid, varid, meta_5d, fill_value)

    ! end of definition
    CALL check_netcdf(nf90_enddef(ncid), __FILE__, __LINE__ )

    ! put variable to netcdf file

    CALL check_netcdf(nf90_put_var(ncid,varid,var_int_5d), __FILE__, __LINE__ )

  END SUBROUTINE netcdf_put_int_i4_5d

  !----------------------------------------------------------------------------------------

  !> define variable with grid defintion for netcdf cf standard
  !! following the cf conventions
  !! see  see e.g. http://cf-pcmdi.llnl.gov/documents/cf-conventions/1.4/

  SUBROUTINE netcdf_def_grid_mapping(ncid, nc_grid_def, varid)

    INTEGER, INTENT(IN)                   :: ncid        !< id for netcdf file
    TYPE(netcdf_grid_mapping), INTENT(IN) :: nc_grid_def !< mapping parameters for netcdf
    INTEGER, INTENT(OUT)                  :: varid       !< netcdf varid of variable

    ! local variables
    INTEGER :: n !< counter

    ! redef netcdf file
    CALL check_netcdf(nf90_redef(ncid), __FILE__, __LINE__ )

    ! define "dummy" variable with grid defintion attributes for netcdf, following the cf conventions
    ! see  see e.g. http://cf-pcmdi.llnl.gov/documents/cf-conventions/1.4/
    CALL check_netcdf( nf90_def_var( ncid,TRIM(nc_grid_def%grid_mapping_varname), &
         &                           NF90_CHAR, varid), __FILE__, __LINE__ )

    CALL check_netcdf( nf90_put_att( ncid,varid, TRIM(nc_grid_def%grid_mapping_name%attname), &
         &                           TRIM(nc_grid_def%grid_mapping_name%attributetext)), __FILE__, __LINE__ )

    DO n=1, nc_grid_def%n_r_att
      CALL check_netcdf( nf90_put_att( ncid,varid,                             &
           &                           TRIM(nc_grid_def%map_param(n)%attname), &
           &                           nc_grid_def%map_param(n)%att_value_r), __FILE__, __LINE__ )
    ENDDO

    ! end of definition
    CALL check_netcdf(nf90_enddef(ncid), __FILE__, __LINE__ )

  END SUBROUTINE netcdf_def_grid_mapping

  !-----------------------------------------------------------------------------

  !> specific subroutine to read 2D real variable from netcdf file
  SUBROUTINE netcdf_get_var_real_2d(path_netcdf_file,var_real_2d_meta,var_real_2d)

    CHARACTER (len=*), INTENT(IN)   :: path_netcdf_file
    TYPE(var_meta_info), INTENT(IN) :: var_real_2d_meta !< meta information for variable

    REAL (KIND=wp), INTENT(OUT)     :: var_real_2d(1:var_real_2d_meta%diminfo(1)%dimsize, &
        &                                          1:var_real_2d_meta%diminfo(2)%dimsize)


    !local variables
    INTEGER            :: ncid    !< id for netcdf file
    INTEGER            :: n       !< counter

    CHARACTER (len=20) :: varname !< name of variable
    INTEGER            :: varid   !< netcdf varid of variable

    CHARACTER (len=12) :: dimname !< name of dimension
    INTEGER            :: dimid   !< id of dimension

    INTEGER            :: ndim    !< number of dimensions of variable
    INTEGER            :: length  !< length of dimension

    ! open netcdf file
    CALL check_netcdf(nf90_open(TRIM(path_netcdf_file),NF90_NOWRITE, ncid), __FILE__, __LINE__ )

    ! first get information for variable
    varname = TRIM(var_real_2d_meta%varname)
    CALL check_netcdf(nf90_inq_varid(ncid, TRIM(varname), varid), __FILE__, __LINE__ )

    ! second  check for dimension size
    ndim = var_real_2d_meta%n_dim
    DO n=1,ndim
      dimname = TRIM(var_real_2d_meta%diminfo(n)%dimname)
      CALL check_netcdf(nf90_inq_dimid(ncid,TRIM(dimname), dimid), __FILE__, __LINE__ )
      CALL check_netcdf(nf90_inquire_dimension(ncid,dimid,len=length), __FILE__, __LINE__ )
      IF (length /= var_real_2d_meta%diminfo(n)%dimsize) THEN
        write(0,*) 'netcdf_get_var_real_2d',n,length,var_real_2d_meta%diminfo(n)%dimsize
        CALL abort_extpar('Dimension size of input file in variable does not match')
      ENDIF
    ENDDO

    ! third get variable
    CALL check_netcdf(nf90_get_var(ncid,varid,var_real_2d), __FILE__, __LINE__ )

    ! close netcdf file
    CALL check_netcdf(nf90_close(ncid), __FILE__, __LINE__ )

  END SUBROUTINE netcdf_get_var_real_2d

  !-----------------------------------------------------------------------------

  !> specific subroutine to read 3D real variable from netcdf file
  SUBROUTINE netcdf_get_var_real_3d(path_netcdf_file,var_real_3d_meta,var_real_3d)

    CHARACTER (len=*), INTENT(IN)   :: path_netcdf_file
    TYPE(var_meta_info), INTENT(IN) :: var_real_3d_meta !< meta information for variable

    REAL (KIND=wp), INTENT(OUT)     :: var_real_3d(1:var_real_3d_meta%diminfo(1)%dimsize, &
         &                                         1:var_real_3d_meta%diminfo(2)%dimsize, &
         &                                         1:var_real_3d_meta%diminfo(3)%dimsize) !< 3D real variable

    !local variables
    INTEGER            :: ncid    !< id for netcdf file
    INTEGER            :: n       !< counter

    CHARACTER (len=20) :: varname !< name of variable
    INTEGER            :: varid   !< netcdf varid of variable

    CHARACTER (len=12) :: dimname !< name of dimension
    INTEGER            :: dimid   !< id of dimension

    INTEGER            :: ndim    !< number of dimensions of variable
    INTEGER            :: length  !< length of dimension


    ! open netcdf file
    CALL check_netcdf(nf90_open(TRIM(path_netcdf_file),NF90_NOWRITE, ncid))

    ! first get information for variable
    varname = TRIM(var_real_3d_meta%varname)
    CALL check_netcdf(nf90_inq_varid(ncid, TRIM(varname), varid), __FILE__, __LINE__ )

    ! second  check for dimension size
    ndim = var_real_3d_meta%n_dim
    DO n=1,ndim
      dimname = TRIM(var_real_3d_meta%diminfo(n)%dimname)
      CALL check_netcdf(nf90_inq_dimid(ncid,TRIM(dimname), dimid), __FILE__, __LINE__ )
      CALL check_netcdf(nf90_inquire_dimension(ncid,dimid,len=length), __FILE__, __LINE__ )
      IF (length /= var_real_3d_meta%diminfo(n)%dimsize) THEN
        write(0,*) 'netcdf_get_var_real_3d',n,length,var_real_3d_meta%diminfo(n)%dimsize
        CALL abort_extpar('Dimension size of input file in variable does not match')
      ENDIF
    ENDDO

    ! third get variable
    CALL check_netcdf(nf90_get_var(ncid,varid,var_real_3d), __FILE__, __LINE__ )

    ! close netcdf file
    CALL check_netcdf(nf90_close(ncid), __FILE__, __LINE__ )

  END SUBROUTINE netcdf_get_var_real_3d

  !-----------------------------------------------------------------------------

  !> specific subroutine to read 4D real variable from netcdf file
  SUBROUTINE netcdf_get_var_real_4d(path_netcdf_file,var_real_4d_meta,var_real_4d)

    CHARACTER (len=*), INTENT(IN)   :: path_netcdf_file
    TYPE(var_meta_info), INTENT(IN) :: var_real_4d_meta !< meta information for variable

    REAL (KIND=wp), INTENT(OUT)     :: var_real_4d(1:var_real_4d_meta%diminfo(1)%dimsize, &
         &                                         1:var_real_4d_meta%diminfo(2)%dimsize, &
         &                                         1:var_real_4d_meta%diminfo(3)%dimsize, &
         &                                         1:var_real_4d_meta%diminfo(4)%dimsize) !< 4D real variable

    !local variables
    INTEGER            :: ncid    !< id for netcdf file
    INTEGER            :: n       !< counter

    CHARACTER (len=20) :: varname !< name of variable
    INTEGER            :: varid   !< netcdf varid of variable

    CHARACTER (len=12) :: dimname !< name of dimension
    INTEGER            :: dimid   !< id of dimension

    INTEGER            :: ndim    !< number of dimensions of variable
    INTEGER            :: length  !< length of dimension

    ! open netcdf file
    CALL check_netcdf(nf90_open(TRIM(path_netcdf_file),NF90_NOWRITE, ncid), __FILE__, __LINE__ )

    ! first get information for variable
    varname = TRIM(var_real_4d_meta%varname)

    CALL check_netcdf(nf90_inq_varid(ncid, TRIM(varname), varid), __FILE__, __LINE__ )

    ! second  check for dimension size
    ndim = var_real_4d_meta%n_dim
    DO n=1,ndim
      dimname = TRIM(var_real_4d_meta%diminfo(n)%dimname)
      CALL check_netcdf(nf90_inq_dimid(ncid,TRIM(dimname), dimid), __FILE__, __LINE__ )
      CALL check_netcdf(nf90_inquire_dimension(ncid,dimid,len=length), __FILE__, __LINE__ )
      IF (length /= var_real_4d_meta%diminfo(n)%dimsize) THEN
        write(0,*) 'netcdf_get_var_real_4d',n,length,var_real_4d_meta%diminfo(n)%dimsize
        CALL abort_extpar('Dimension size of input file in variable does not match')
      ENDIF
    ENDDO

    ! third get variable
    CALL check_netcdf(nf90_get_var(ncid,varid,var_real_4d), __FILE__, __LINE__ )

    ! close netcdf file
    CALL check_netcdf(nf90_close(ncid), __FILE__, __LINE__)

  END SUBROUTINE netcdf_get_var_real_4d

  !-----------------------------------------------------------------------------

  !> specific subroutine to read 4D real variable from netcdf file
  SUBROUTINE netcdf_get_var_real_5d(path_netcdf_file,var_real_5d_meta,var_real_5d)

    CHARACTER (len=*), INTENT(IN)   :: path_netcdf_file
    TYPE(var_meta_info), INTENT(IN) :: var_real_5d_meta !< meta information for variable

    REAL (KIND=wp), INTENT(OUT)     :: var_real_5d(1:var_real_5d_meta%diminfo(1)%dimsize, &
         &                                         1:var_real_5d_meta%diminfo(2)%dimsize, &
         &                                         1:var_real_5d_meta%diminfo(3)%dimsize, &
         &                                         1:var_real_5d_meta%diminfo(4)%dimsize, &
         &                                         1:var_real_5d_meta%diminfo(5)%dimsize) !< 4D real variable

    !local variables
    INTEGER            :: ncid    !< id for netcdf file
    INTEGER            :: n       !< counter

    CHARACTER (len=20) :: varname !< name of variable
    INTEGER            :: varid   !< netcdf varid of variable

    CHARACTER (len=12) :: dimname !< name of dimension
    INTEGER            :: dimid   !< id of dimension

    INTEGER            :: ndim    !< number of dimensions of variable
    INTEGER            :: length  !< length of dimension

    ! open netcdf file
    CALL check_netcdf(nf90_open(TRIM(path_netcdf_file),NF90_NOWRITE, ncid), __FILE__, __LINE__ )

    ! first get information for variable
    varname = TRIM(var_real_5d_meta%varname)

    CALL check_netcdf(nf90_inq_varid(ncid, TRIM(varname), varid), __FILE__, __LINE__ )

    ! second  check for dimension size
    ndim = var_real_5d_meta%n_dim
    DO n=1,ndim
      dimname = TRIM(var_real_5d_meta%diminfo(n)%dimname)
      CALL check_netcdf(nf90_inq_dimid(ncid,TRIM(dimname), dimid), __FILE__, __LINE__ )
      CALL check_netcdf(nf90_inquire_dimension(ncid,dimid,len=length), __FILE__, __LINE__ )
      IF (length /= var_real_5d_meta%diminfo(n)%dimsize) THEN
        write(0,*) 'netcdf_get_var_real_5d',n,length,var_real_5d_meta%diminfo(n)%dimsize
        CALL abort_extpar('Dimension size of input file in variable does not match')
      ENDIF
    ENDDO

    ! third get variable
    CALL check_netcdf(nf90_get_var(ncid,varid,var_real_5d), __FILE__, __LINE__ )

    ! close netcdf file
    CALL check_netcdf(nf90_close(ncid), __FILE__, __LINE__ )

  END SUBROUTINE netcdf_get_var_real_5d

  !-----------------------------------------------------------------------------

  !> specific subroutine to read 3D integer variable from netcdf file
  SUBROUTINE netcdf_get_var_int_3d_i8(path_netcdf_file,var_int_3d_meta,var_int_3d)

    CHARACTER (len=*), INTENT(IN)   :: path_netcdf_file
    TYPE(var_meta_info), INTENT(IN) :: var_int_3d_meta !< meta information for variable

    INTEGER (KIND=i8), INTENT(OUT)  :: var_int_3d(1:var_int_3d_meta%diminfo(1)%dimsize, &
         &                                        1:var_int_3d_meta%diminfo(2)%dimsize, &
         &                                        1:var_int_3d_meta%diminfo(3)%dimsize) !< 3D integer variable

    !local variables
    INTEGER            :: ncid    !< id for netcdf file
    INTEGER            :: n       !< counter

    CHARACTER (len=20) :: varname !< name of variable
    INTEGER            :: varid   !< netcdf varid of variable

    CHARACTER (len=12) :: dimname !< name of dimension
    INTEGER            :: dimid   !< id of dimension

    INTEGER            :: ndim    !< number of dimensions of variable
    INTEGER            :: length  !< length of dimension

    ! open netcdf file
    CALL check_netcdf(nf90_open(TRIM(path_netcdf_file),NF90_NOWRITE, ncid), __FILE__, __LINE__ )

    ! first get information for variable
    varname = TRIM(var_int_3d_meta%varname)
    CALL check_netcdf(nf90_inq_varid(ncid, TRIM(varname), varid), __FILE__, __LINE__ )

    ! second  check for dimension size
    ndim = var_int_3d_meta%n_dim
    DO n=1,ndim
      dimname = TRIM(var_int_3d_meta%diminfo(n)%dimname)
      CALL check_netcdf(nf90_inq_dimid(ncid,TRIM(dimname), dimid), __FILE__, __LINE__ )
      CALL check_netcdf(nf90_inquire_dimension(ncid,dimid,len=length), __FILE__, __LINE__ )
      IF (length /= var_int_3d_meta%diminfo(n)%dimsize) THEN
        write(0,*) 'netcdf_get_var_int_3di8',n,length,var_int_3d_meta%diminfo(n)%dimsize
        CALL abort_extpar('Dimension size of input file in variable does not match')
      ENDIF
    ENDDO

    ! third get variable
    CALL check_netcdf(nf90_get_var(ncid,varid,var_int_3d), __FILE__, __LINE__ )

    ! close netcdf file
    CALL check_netcdf(nf90_close(ncid), __FILE__, __LINE__ )

  END SUBROUTINE netcdf_get_var_int_3d_i8

  !-----------------------------------------------------------------------------

  !> specific subroutine to read 3D integer variable from netcdf file
  SUBROUTINE netcdf_get_var_int_3d_i4(path_netcdf_file,var_int_3d_meta,var_int_3d)

    CHARACTER (len=*), INTENT(IN)   :: path_netcdf_file
    TYPE(var_meta_info), INTENT(IN) :: var_int_3d_meta !< meta information for variable

    INTEGER (KIND=i4), INTENT(OUT)  :: var_int_3d(1:var_int_3d_meta%diminfo(1)%dimsize, &
         &                                        1:var_int_3d_meta%diminfo(2)%dimsize, &
         &                                        1:var_int_3d_meta%diminfo(3)%dimsize) !< 3D integer variable

    !local variables
    INTEGER            :: ncid    !< id for netcdf file
    INTEGER            :: n       !< counter

    CHARACTER (len=20) :: varname !< name of variable
    INTEGER            :: varid   !< netcdf varid of variable

    CHARACTER (len=12) :: dimname !< name of dimension
    INTEGER            :: dimid   !< id of dimension

    INTEGER            :: ndim    !< number of dimensions of variable
    INTEGER            :: length  !< length of dimension

    ! open netcdf file
    CALL check_netcdf(nf90_open(TRIM(path_netcdf_file),NF90_NOWRITE, ncid))

    ! first get information for variable
    varname = TRIM(var_int_3d_meta%varname)
    CALL check_netcdf(nf90_inq_varid(ncid, TRIM(varname), varid), __FILE__, __LINE__ )

    ! second  check for dimension size
    ndim = var_int_3d_meta%n_dim
    DO n=1,ndim
      dimname = TRIM(var_int_3d_meta%diminfo(n)%dimname)
      CALL check_netcdf(nf90_inq_dimid(ncid,TRIM(dimname), dimid), __FILE__, __LINE__ )
      CALL check_netcdf(nf90_inquire_dimension(ncid,dimid,len=length), __FILE__, __LINE__ )
      IF (length /= var_int_3d_meta%diminfo(n)%dimsize) THEN
        write(0,*) 'netcdf_get_var_int_3d',n,length,var_int_3d_meta%diminfo(n)%dimsize
        CALL abort_extpar('Dimension size of input file in variable does not match')
      ENDIF
    ENDDO

    ! third get variable
    CALL check_netcdf(nf90_get_var(ncid,varid,var_int_3d), __FILE__, __LINE__ )

    ! close netcdf file
    CALL check_netcdf(nf90_close(ncid), __FILE__, __LINE__ )

  END SUBROUTINE netcdf_get_var_int_3d_i4

  !-----------------------------------------------------------------------------

  !> specific subroutine to read 4D integer variable from netcdf file
  SUBROUTINE netcdf_get_var_int_4d(path_netcdf_file,var_int_4d_meta,var_int_4d)

    CHARACTER (len=*), INTENT(IN)   :: path_netcdf_file
    TYPE(var_meta_info), INTENT(IN) :: var_int_4d_meta !< meta information for variable

    INTEGER (KIND=i8) :: var_int_4d(1:var_int_4d_meta%diminfo(1)%dimsize, &
         &                          1:var_int_4d_meta%diminfo(2)%dimsize, &
         &                          1:var_int_4d_meta%diminfo(3)%dimsize, &
         &                          1:var_int_4d_meta%diminfo(4)%dimsize) !< 4D integer variable

    !local variables
    INTEGER :: ncid !< id for netcdf file
    INTEGER :: n !< counter

    CHARACTER (len=20) :: varname !< name of variable
    INTEGER :: varid !< netcdf varid of variable

    CHARACTER (len=12) :: dimname  !< name of dimension
    INTEGER :: dimid  !< id of dimension

    INTEGER :: ndim  !< number of dimensions of variable
    INTEGER :: length!< length of dimension

    ! open netcdf file
    CALL check_netcdf(nf90_open(TRIM(path_netcdf_file),NF90_NOWRITE, ncid), __FILE__, __LINE__ )

    ! first get information for variable
    varname = TRIM(var_int_4d_meta%varname)
    CALL check_netcdf(nf90_inq_varid(ncid, TRIM(varname), varid), __FILE__, __LINE__ )

    ! second  check for dimension size
    ndim = var_int_4d_meta%n_dim
    DO n=1,ndim
      dimname = TRIM(var_int_4d_meta%diminfo(n)%dimname)
      CALL check_netcdf(nf90_inq_dimid(ncid,TRIM(dimname), dimid), __FILE__, __LINE__ )
      CALL check_netcdf(nf90_inquire_dimension(ncid,dimid,len=length), __FILE__, __LINE__ )
      IF (length /= var_int_4d_meta%diminfo(n)%dimsize) THEN
        write(0,*) 'netcdf_get_var_int_4d',n,length,var_int_4d_meta%diminfo(n)%dimsize
        CALL abort_extpar('Dimension size of input file in variable does not match')
      ENDIF
    ENDDO

    ! third get variable
    CALL check_netcdf(nf90_get_var(ncid,varid,var_int_4d), __FILE__, __LINE__ )

    ! close netcdf file
    CALL check_netcdf(nf90_close(ncid), __FILE__, __LINE__ )

  END SUBROUTINE netcdf_get_var_int_4d

  !-----------------------------------------------------------------------------

  !> set date for an invariant field of the external parameters in a GRIB message
  !! the convention at DWD is to set the date for the invariant fields to
  !! year 1, january 1, 00:00 hour
  SUBROUTINE get_date_const_field(dataDate,dataTime)

    INTEGER (KIND=i8), INTENT(OUT)  :: dataDate
                                                 !< date, for edition independent use of GRIB_API dataDate
                                                 !< as Integer in the format ccyymmdd
    INTEGER (KIND=i8), INTENT(OUT)  :: dataTime  !< time, for edition independent use GRIB_API dataTime in the format hhmm

    !local variables

    INTEGER  :: cc
    INTEGER  :: yy
    INTEGER  :: mm
    INTEGER  :: dd
    INTEGER  :: hh
    INTEGER  :: minute

    cc     = 0
    yy     = 1
    mm     = 1
    dd     = 1
    hh     = 0
    minute = 0

    ! dataDate format ccyymmdd
    ! dataTime format hhmm (hour minute)
    dataDate = dd + (100*mm) + (10000*yy) + (1000000*cc)
    dataTime = minute + (100*hh)

  END  SUBROUTINE get_date_const_field

  !-----------------------------------------------------------------------------

  !> set date for an monthly climatology field of the external parameters in a GRIB message with GRIB_API
  !! the convention at DWD is to set the date for the invariant fields to
  !! year 1111, month mm, day 11, 00:00 hour
  !! the grib message should have been previously defined, pass the grib_id to this subroutine
  SUBROUTINE set_date_mm_extpar_field(mm,dataDate,dataTime)

    INTEGER, INTENT(IN)             :: mm        !< month
    INTEGER (KIND=i8), INTENT(OUT)  :: dataDate
                                                 !< date, for edition independent use of GRIB_API dataDate
                                                 !< as Integer in the format ccyymmdd
    INTEGER (KIND=i8), INTENT(OUT)  :: dataTime  !< time, for edition independent use GRIB_API dataTime in the format hhmm
    !local variables

    INTEGER  :: cc
    INTEGER  :: yy
    INTEGER  :: dd
    INTEGER  :: hh
    INTEGER  :: minute

    cc = 11
    yy = 11
    dd = 11
    hh = 0
    minute=0
    ! dataDate format ccyymmdd
    ! dataTime format hhmm (hour minute)
    dataDate = dd + (100*mm) + (10000*yy) + (1000000*cc)
    dataTime = minute + (100*hh)

  END  SUBROUTINE set_date_mm_extpar_field

  !-----------------------------------------------------------------------------

  !> open netcdf-file and get netcdf unit file number ncid
  SUBROUTINE open_new_netcdf_file(netcdf_filename, dim_list, global_attributes, time, ncid)

    CHARACTER (len=*), INTENT(IN)                 :: netcdf_filename       !< filename for the netcdf file
    TYPE(dim_meta_info), INTENT(IN)               :: dim_list(:)           !< dimensions for netcdf file
    TYPE(netcdf_attributes), INTENT(IN), OPTIONAL :: global_attributes(:)  !< structure with global attributes
    REAL (KIND=wp), INTENT(IN), OPTIONAL          :: time(:)               !< time variable

    INTEGER, INTENT(OUT)                          :: ncid                  !< netcdf unit file number

    ! local variables
    INTEGER                :: ndims      !< number of dimension
    INTEGER                :: ng_att     !< number of global attributes
    INTEGER, ALLOCATABLE   :: dimids(:)   !< list of netcdf dim ids
    INTEGER                :: call_mode
    INTEGER                :: dimsize
    INTEGER                :: errorcode
    INTEGER                :: n           !< counter

    CHARACTER (len=20)     :: varname     !< name of variable
    CHARACTER (len=12)     :: dimname     !< name of dimension

    INTEGER                :: varid_time  !< netcdf varid of variable
    INTEGER                :: dimid_time

    INTEGER                :: varid_mlev  !< netcdf varid of variable
    INTEGER                :: dimid_mlev
    INTEGER, ALLOCATABLE   :: mlev(:)
    INTEGER                :: number_of_grid_used_int
    CHARACTER (LEN=4)      :: number_of_grid_used_char

    ! Comment from Merge write an option for netcdf file type
    call_mode = NF90_CLOBBER + NF90_64BIT_OFFSET
    CALL check_netcdf( nf90_create(TRIM(netcdf_filename),call_mode,ncid), __FILE__, __LINE__ )
    !    CALL check_netcdf( nf90_create(TRIM(netcdf_filename),NF90_CLOBBER,ncid))
    CALL check_netcdf( nf90_create(TRIM(netcdf_filename),NF90_NETCDF4,ncid), __FILE__, __LINE__ )

    ndims = SIZE(dim_list)

    dimid_time = -1

    ALLOCATE(dimids(1:ndims), STAT=errorcode)
    IF(errorcode /= 0) CALL abort_extpar('Cant allocate dimids')

    ! define dimensions
    DO n=1, ndims
      dimsize = dim_list(n)%dimsize
      IF (TRIM(dim_list(n)%dimname)=='time') dimsize = NF90_UNLIMITED
      CALL check_netcdf( nf90_def_dim( ncid,                      &
           &                          TRIM(dim_list(n)%dimname), &
           &                          dimsize,                   &
           &                          dimids(n)), __FILE__, __LINE__ )
      IF (TRIM(dim_list(n)%dimname)=='time') dimid_time=dimids(n)
    ENDDO

    ! put global attributes
    IF (PRESENT(global_attributes)) THEN
      ng_att=SIZE(global_attributes)

      DO n=1, ng_att
        IF(TRIM(global_attributes(n)%attname)=='number_of_grid_used') THEN
          number_of_grid_used_char=TRIM(global_attributes(n)%attributetext)
          READ(number_of_grid_used_char,'(I4)') number_of_grid_used_int
          CALL check_netcdf( nf90_put_att( ncid, NF90_GLOBAL,                  &
               &                          TRIM(global_attributes(n)%attname), &
               &                          number_of_grid_used_int), __FILE__, __LINE__ )
        ELSE
          CALL check_netcdf( nf90_put_att( ncid,NF90_GLOBAL,                   &
               &                          TRIM(global_attributes(n)%attname), &
               &                          TRIM(global_attributes(n)%attributetext)), __FILE__, __LINE__ )
        END IF
      ENDDO
    ENDIF

    ! end of definition
    CALL check_netcdf(nf90_enddef(ncid), __FILE__, __LINE__ )

    IF (PRESENT(time)) THEN
      CALL check_netcdf(nf90_redef(ncid), __FILE__, __LINE__ )
      ! dimsize=SIZE(time)
      dimname='time'
      varname='time'
      IF (dimid_time == -1) THEN
        CALL check_netcdf(nf90_def_dim(ncid,dimname,NF90_UNLIMITED, dimid_time), __FILE__, __LINE__ )
      ENDIF
      ! define netcdf variable
      CALL check_netcdf( nf90_def_var(ncid,varname,NF90_FLOAT,dimid_time,varid_time), __FILE__, __LINE__ )

      dimname='mlev'
      varname='mlev'
      dimsize=1
      ALLOCATE(mlev(1:dimsize), STAT=errorcode)
      IF(errorcode /= 0) CALL abort_extpar('Cant allocate mlev')
      mlev=1
      CALL check_netcdf(nf90_def_dim(ncid,dimname,dimsize, dimid_mlev), __FILE__, __LINE__ )
      CALL check_netcdf( nf90_def_var(ncid,varname,NF90_FLOAT,dimid_mlev,varid_mlev), __FILE__, __LINE__ )

      CALL check_netcdf(nf90_enddef(ncid))

      ! put time variable to netcdf file
      CALL check_netcdf(nf90_put_var(ncid,varid_time,time), __FILE__, __LINE__ )
      CALL check_netcdf(nf90_put_var(ncid,varid_mlev,mlev), __FILE__, __LINE__ )

    ENDIF

  END SUBROUTINE open_new_netcdf_file

  !-----------------------------------------------------------------------------

  !> close netcdf-file with unit file number ncid
  SUBROUTINE close_netcdf_file(ncid)
    INTEGER, INTENT(IN) :: ncid                       !< netcdf unit file number
    !! close netcdf file
    CALL check_netcdf( nf90_close( ncid), __FILE__, __LINE__ )
  END SUBROUTINE close_netcdf_file

END MODULE mo_io_utilities
