!+  Fortran module for netcdf output of soil data
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_0         2010/12/21 Hermann Asensio
!  Initial release
! V1_4         2011/04/21 Hermann Asensio
!  clean up
! V2_0         1013/06/04 Martina Messmer
!  adaptations in a way that HWSD data set can be used (top- and subsoil)
!  Code received from Juergen Helmert
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!> Fortran module for netcdf output of soil data
!> \author Hermann Asensio
MODULE mo_soil_output_nc

  USE mo_logging
  USE mo_kind,                  ONLY: wp, i4

  USE mo_grid_structures,       ONLY: target_grid_def

  USE mo_io_utilities,         ONLY: netcdf_attributes, &
       &                             dim_meta_info, &
       &                             netcdf_put_var, &
       &                             open_new_netcdf_file, &
       &                             close_netcdf_file, &
       &                             netcdf_get_var, &
       &                             dim_meta_info

  USE mo_soil_data,            ONLY: FAO_data, HWSD_data

  USE mo_var_meta_data,        ONLY: def_dimension_info_buffer, &
       &                             lon_geo_meta, &
       &                             lat_geo_meta, &
       &                             def_com_target_fields_meta, &  
       &                             def_soil_meta, &
       &                             fr_land_soil_meta, &
       &                             soiltype_fao_meta, &
       &                             def_com_target_fields_meta, &
       &                             soiltype_FAO_deep_meta, &
       &                             dim_3d_tg,&
       &                             lon_geo_meta, &
       &                             lat_geo_meta, &
       &                             soiltype_hwsd_meta, &
       &                             soiltype_HWSD_deep_meta
  
  
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: write_netcdf_soil_buffer, &
       &    read_netcdf_soil_buffer

  CONTAINS

  !----------------------------------------------------------------------------
  !> set global attributes for netcdf with soiltype data
  SUBROUTINE set_global_att_soiltype(global_attributes,isoil_data)


    TYPE(netcdf_attributes), INTENT(INOUT) :: global_attributes(1:6)
    INTEGER (KIND=i4), INTENT(IN):: isoil_data

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
    global_attributes(1)%attributetext='Soil type'
    global_attributes(2)%attname = 'institution'
    global_attributes(2)%attributetext='DWD'

    global_attributes(3)%attname = 'source'
    IF (isoil_data == FAO_data) THEN
    global_attributes(3)%attributetext='FAO Digital Soil Map of the World'
    ELSE
      global_attributes(3)%attributetext='HWSD Digital Soil Map of the World'
    ENDIF

    CALL DATE_AND_TIME(ydate,ytime)
    READ(ydate,'(4A2)') cc,yy,mm,dd
    READ(ytime,'(2A2)') hh, minute

    ydate=TRIM(cc)//TRIM(yy)//'-'//TRIM(mm)//'-'//TRIM(dd)
    ytime=TRIM(hh)//':'//TRIM(minute) 

    global_attributes(4)%attname = 'history'
    global_attributes(4)%attributetext=TRIM(ydate)//'T'//TRIM(ytime)//' soil_to_buffer'

    global_attributes(5)%attname = 'references'
    IF (isoil_data == FAO_data) THEN
    global_attributes(5)%attributetext='FAO Digital Soil Map of the World'
    ELSE
      global_attributes(5)%attributetext='HWSD Digital Soil Map of the World'
    ENDIF

    global_attributes(6)%attname = 'comment'
    IF (isoil_data /= HWSD_data) THEN
    global_attributes(6)%attributetext='1 ice, 2 rock, 3 sand, 4 sandy loam, 5 loam,&
                                      & 6 loamy clay, 7 clay, 8 histosol(e.g. peat), 9 sea point'
    ELSE
      global_attributes(6)%attributetext='HWSD soil types'
    ENDIF

  END SUBROUTINE set_global_att_soiltype
  !----------------------------------------------------------------------------

  !> create a netcdf file for the FAO soil data in the buffer
  SUBROUTINE write_netcdf_soil_buffer(netcdf_filename,  &
       &                                 tg,           &
       &                                 isoil_data,   &
       &                                 ldeep_soil,   &
       &                                 undefined,    &
       &                                 undef_int,    &
       &                                 lon_geo,      &
       &                                 lat_geo,      &
       &                                 fr_land_soil, &
       &                                 soiltype_fao, &
       &                                 soiltype_hwsd, &
       &                                 soiltype_fao_deep, &
       &                                 soiltype_hwsd_deep  )

    CHARACTER (len=*), INTENT(IN)             :: netcdf_filename !< filename for the netcdf file
                                              
    TYPE(target_grid_def), INTENT(INOUT)      :: tg !< structure which contains various target grid fields
                                              
    INTEGER (KIND=i4), INTENT(IN)             :: isoil_data, &
         &                                       undef_int       !< value to indicate undefined grid elements
                                              
    REAL (KIND=wp), INTENT(IN)                :: lon_geo(:,:,:), & !< longitude coordinates of the target grid
         &                                       undefined, & !< value to indicate undefined grid elements 
         &                                       lat_geo(:,:,:)!< latitude coordinates of the target grid
                                              
    LOGICAL,       INTENT(IN)                 :: ldeep_soil
                                              
    REAL(KIND=wp), INTENT(INOUT)              :: fr_land_soil(:,:,:) !< fraction land due to FAO Digital Soil map of the World
                                              
    INTEGER(KIND=i4), INTENT(INOUT)           :: soiltype_fao(:,:,:), & !< soiltype due to FAO Digital Soil map of the World
         &                                       soiltype_hwsd(:,:,:) !< soiltype due to FAO Digital Soil map of the World

    INTEGER(KIND=i4), INTENT(INOUT), OPTIONAL :: soiltype_fao_deep(:,:,:), & !< soiltype due to FAO Digital Soil map of the World
                                                 soiltype_hwsd_deep(:,:,:) !< soiltype due to HWSD

    ! local variables
    INTEGER(KIND=i4)                          :: ndims, ncid, errorcode

    INTEGER(KIND=i4), PARAMETER               :: nglob_atts=6

    TYPE(dim_meta_info), ALLOCATABLE          :: dim_list(:) !< dimensions for netcdf file

    TYPE(netcdf_attributes)                   :: global_attributes(nglob_atts)

    CALL logging%info('Enter routine: write_netcdf_soil_buffer')
    !-------------------------------------------------------------
    ! define global attributes
    CALL set_global_att_soiltype(global_attributes,isoil_data)

    !set up dimensions for buffer
    CALL  def_dimension_info_buffer(tg)
    ! dim_3d_tg

    ! define meta information for target field variables lon_geo, lat_geo 
    CALL def_com_target_fields_meta(dim_3d_tg)
    ! lon_geo_meta and lat_geo_meta

    CALL def_soil_meta(dim_3d_tg,isoil_data)
    !  fr_land_soil_meta, soiltype_fao_meta

    !set up dimensions for buffer
    ndims = 3
    ALLOCATE(dim_list(1:ndims),STAT=errorcode)
    IF (errorcode /= 0 ) CALL logging%error('Cant allocate array dim_list',__FILE__,__LINE__)
    dim_list = dim_3d_tg

    !-----------------------------------------------------------------

    CALL open_new_netcdf_file(netcdf_filename=TRIM(netcdf_filename),   &
      &                       dim_list=dim_list,                  &
      &                       global_attributes=global_attributes, &
      &                       ncid=ncid)


    ! soiltype FAO
    CALL netcdf_put_var(ncid,soiltype_fao,soiltype_fao_meta,undef_int)
    ! soil Id HWSD 
    CALL netcdf_put_var(ncid,soiltype_hwsd,soiltype_hwsd_meta,undef_int)

    IF (ldeep_soil) THEN
      CALL netcdf_put_var(ncid,soiltype_fao_deep,soiltype_fao_deep_meta,undef_int)
      CALL netcdf_put_var(ncid,soiltype_hwsd_deep,soiltype_hwsd_deep_meta,undef_int)
    ENDIF

     ! lon
    CALL netcdf_put_var(ncid,lon_geo,lon_geo_meta,undefined)

    ! lat
    CALL netcdf_put_var(ncid,lat_geo,lat_geo_meta,undefined)

     ! fr_land_soil
    CALL netcdf_put_var(ncid,fr_land_soil,fr_land_soil_meta,undefined)

    CALL close_netcdf_file(ncid)

    CALL logging%info('Exit routine: write_netcdf_soil_buffer')

  END SUBROUTINE write_netcdf_soil_buffer
   !----------------------------------------------------------------------------
  !> read FAO soil data buffer fields
  SUBROUTINE read_netcdf_soil_buffer(netcdf_filename,   &
       &                                 tg,           &
       &                                 isoil_data,   &
       &                                 fr_land_soil, &
       &                                 soiltype_fao, &
       &                                 soiltype_hwsd, &
       &                                 soiltype_fao_deep, &
       &                                 soiltype_hwsd_deep  )


    CHARACTER (len=*), INTENT(IN)           :: netcdf_filename !< filename for the netcdf file
    TYPE(target_grid_def), INTENT(IN)       :: tg !< structure which contains various target grid fields
    INTEGER (KIND=i4), INTENT(IN)           :: isoil_data

    REAL(KIND=wp), INTENT(OUT)              :: fr_land_soil(:,:,:) !< fraction land due to FAO Digital Soil map of the World
                                            
    INTEGER(KIND=i4), INTENT(OUT)           :: soiltype_fao(:,:,:), & !< soiltype due to FAO Digital Soil map of the World
         &                                     soiltype_hwsd(:,:,:) !< soiltype due to FAO Digital Soil map of the World

    INTEGER(KIND=i4), INTENT(OUT), OPTIONAL :: soiltype_FAO_deep(:,:,:), & !< soiltype due to FAO Digital Soil map of the World
         &                                     soiltype_HWSD_deep(:,:,:) !< soiltype due to FAO Digital Soil map of the World

    CALL logging%info('Enter routine: read_netcdf_soil_buffer')

    !set up dimensions for buffer
    CALL  def_dimension_info_buffer(tg)
    ! dim_3d_tg

    ! define meta information for target field variables lon_geo, lat_geo 
    CALL def_com_target_fields_meta(dim_3d_tg)

    ! lon_geo_meta and lat_geo_meta
    CALL def_soil_meta(dim_3d_tg, isoil_data)
    !  fr_land_soil_meta, soiltype_fao_meta
    CALL netcdf_get_var(TRIM(netcdf_filename),fr_land_soil_meta,fr_land_soil)

    CALL netcdf_get_var(TRIM(netcdf_filename),soiltype_fao_meta,soiltype_fao)
    CALL netcdf_get_var(TRIM(netcdf_filename),soiltype_hwsd_meta,soiltype_hwsd)

    IF(PRESENT(soiltype_FAO_deep)) THEN
      CALL netcdf_get_var(TRIM(netcdf_filename),soiltype_FAO_deep_meta,soiltype_FAO_deep)
      CALL netcdf_get_var(TRIM(netcdf_filename),soiltype_HWSD_deep_meta,soiltype_HWSD_deep)
    ENDIF

    CALL logging%info('Exit routine: read_netcdf_soil_buffer')

  END SUBROUTINE read_netcdf_soil_buffer

END MODULE mo_soil_output_nc
