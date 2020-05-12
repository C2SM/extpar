!+ Fortran module for netcdf output of AHF data (anthropogenic heat flux)
!
!
! Description:
! Fortran module for netcdf output of AHF data (anthropogenic heat flux)
!
! Current Code Owner: DWD, <Name of person responsible for this code>
!    <smail, phone, fax and email>
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V2_4         2015-05-21 Hendrik Wouters
!  Initial release 
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!=======================================================================
!> Fortran module for netcdf output of AHF data
!> \author Hermann Asensio
MODULE mo_ahf_output_nc

  USE mo_logging
  USE mo_kind,                  ONLY: wp, i4

  USE mo_grid_structures,       ONLY: target_grid_def

  USE mo_io_utilities,          ONLY: netcdf_attributes, &
       &                              dim_meta_info, &
       &                              netcdf_put_var, &
       &                              open_new_netcdf_file, &
       &                              close_netcdf_file, &
       &                              netcdf_get_var

  USE mo_var_meta_data,         ONLY: dim_3d_tg, &
       &                              def_dimension_info_buffer, &
       &                              def_com_target_fields_meta, &  
       &                              ahf_field_meta, &
       &                              def_ahf_meta, &
       &                              lon_geo_meta, &
       &                              lat_geo_meta, &
       &                              def_com_target_fields_meta

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: write_netcdf_buffer_ahf

  PUBLIC :: read_netcdf_buffer_ahf

  CONTAINS

  SUBROUTINE write_netcdf_buffer_ahf(netcdf_filename,  &
   &                                     tg,         &
   &                                     undefined, &
   &                                     lon_geo,     &
   &                                     lat_geo, &
   &                                     ahf_field)




    CHARACTER (len=*), INTENT(IN)      :: netcdf_filename !< filename for the netcdf file
    TYPE(target_grid_def), INTENT(IN)  :: tg !< structure with target grid description
    REAL(KIND=wp), INTENT(IN)          :: undefined, &       !< value to indicate undefined grid elements 
         &                                lon_geo(:,:,:), &  !< longitude coordinates of the target grid in the geographical system
         &                                lat_geo(:,:,:), &  !< latitude coordinates of the target grid in the geographical system
         &                                ahf_field(:,:,:) !< field for monthly mean ahf  data (12 months)

    ! local variables
    INTEGER(KIND=i4)                   :: ndims, ncid, errorcode  
    INTEGER(KIND=i4), PARAMETER        :: nglob_atts=6
                                      
    TYPE(dim_meta_info), ALLOCATABLE   :: dim_list(:) !< dimensions for netcdf file
    TYPE(netcdf_attributes)            :: global_attributes(nglob_atts)

    CALL logging%info('Enter routine: write_netcdf_buffer_ahf')

    !-------------------------------------------------------------
    ! define global attributes
    CALL set_global_att_ahf(global_attributes)

    !set up dimensions for buffer
    CALL  def_dimension_info_buffer(tg)
    
    ! define meta information for target field variables lon_geo, lat_geo 
    CALL def_com_target_fields_meta(dim_3d_tg)
    ! lon_geo_meta and lat_geo_meta
    !define meta information for various AHF data related variables for netcdf output
    CALL def_ahf_meta(dim_3d_tg)
    ! dim_ahf_tg, ahf_field_meta, ahf_field_mom_meta, ahf_ratio_mom_meta
    
    ! set up dimensions for netcdf output 
    ndims = 3
    ALLOCATE(dim_list(1:ndims),STAT=errorcode)
    IF (errorcode /= 0 ) CALL logging%error('Cant allocate array dim_list',__FILE__,__LINE__)
      
    dim_list(1) = dim_3d_tg(1) ! ie
    dim_list(2) = dim_3d_tg(2) ! je
    dim_list(3) = dim_3d_tg(3) ! ke

    !-----------------------------------------------------------------

    CALL open_new_netcdf_file(netcdf_filename=TRIM(netcdf_filename),   &
      &                       dim_list=dim_list,                  &
      &                       global_attributes=global_attributes, &
      &                       ncid=ncid)

    ! lon
    CALL netcdf_put_var(ncid,lon_geo,lon_geo_meta,undefined)

    ! lat
    CALL netcdf_put_var(ncid,lat_geo,lat_geo_meta,undefined)

    ! ahf_field
    CALL netcdf_put_var(ncid,ahf_field,ahf_field_meta,undefined)

    CALL close_netcdf_file(ncid)

    CALL logging%info('Exit routine: write_netcdf_buffer_ahf')

  END SUBROUTINE write_netcdf_buffer_ahf
   !-----------------------------------------------------------------
   !-----------------------------------------------------------------

  !> set global attributes for netcdf with AHF data
  SUBROUTINE set_global_att_ahf(global_attributes)
    TYPE(netcdf_attributes), INTENT(INOUT) :: global_attributes(1:6)

    !local variables
    CHARACTER(len=10) :: ydate
    CHARACTER(len=10) :: ytime
    CHARACTER(len=2)  :: cc
    CHARACTER(len=2)  :: yy
    CHARACTER(len=2)  :: mm
    CHARACTER(len=2)  :: dd
    CHARACTER(len=2)  :: hh
    CHARACTER(len=2)  :: minute

    ! define global attributes
    
    global_attributes(1)%attname = 'title'
    global_attributes(1)%attributetext='AHF data '
    global_attributes(2)%attname = 'institution'
    global_attributes(2)%attributetext='KU Leuven'

    global_attributes(3)%attname = 'source'
    global_attributes(3)%attributetext='NCAR Flanner2009'

    CALL DATE_AND_TIME(ydate,ytime)
    READ(ydate,'(4A2)') cc,yy,mm,dd
    READ(ytime,'(2A2)') hh, minute

    ydate=TRIM(cc)//TRIM(yy)//'-'//TRIM(mm)//'-'//TRIM(dd)
    ytime=TRIM(hh)//':'//TRIM(minute) 

    global_attributes(4)%attname = 'history'
    global_attributes(4)%attributetext=TRIM(ydate)//'T'//TRIM(ytime)//' ahf_to_buffer'

    global_attributes(5)%attname = 'references'
    global_attributes(5)%attributetext='http://www.cgd.ucar.edu/tss/ahf/'

    global_attributes(6)%attname = 'comment'
    global_attributes(6)%attributetext=''

  END SUBROUTINE set_global_att_ahf
  !-----------------------------------------------------------------------

  SUBROUTINE read_netcdf_buffer_ahf(netcdf_filename,  &
   &                                     tg,         &
   &                                     ahf_field)

    CHARACTER (len=*), INTENT(IN)      :: netcdf_filename !< filename for the netcdf file
    TYPE(target_grid_def), INTENT(IN)  :: tg !< structure with target grid description
    REAL (KIND=wp), INTENT(OUT)        :: ahf_field(:,:,:) !< field for ahf 

    ! local variables
    INTEGER(KIND=i4), PARAMETER        :: nglob_atts=6

    CALL logging%info('Enter routine: read_netcdf_buffer_ahf')

    !set up dimensions for buffer
    CALL  def_dimension_info_buffer(tg)

    ! define meta information for target field variables lon_geo, lat_geo 
    CALL def_com_target_fields_meta(dim_3d_tg)
    ! lon_geo_meta and lat_geo_meta
    !define meta information for various AHF data related variables for netcdf output
    CALL def_ahf_meta(dim_3d_tg)
    ! dim_ahf_tg, ahf_field_meta

    CALL netcdf_get_var(TRIM(netcdf_filename),ahf_field_meta,ahf_field)

    CALL logging%info('Exit routine: read_netcdf_buffer_ahf')

   END SUBROUTINE read_netcdf_buffer_ahf
                                                                     
END MODULE mo_ahf_output_nc
