MODULE mo_python_output_nc

  USE mo_logging
  USE mo_kind,                  ONLY: wp, i4

  USE mo_grid_structures,       ONLY: target_grid_def

  USE mo_io_utilities,          ONLY: netcdf_get_var

  USE mo_var_meta_data,         ONLY: dim_3d_tg, &
       &                              def_dimension_info_buffer, & 
       &                              def_com_target_fields_meta, &   
  ! cru
       &                              crutemp_meta,               &
       &                              def_crutemp_meta,           &
       &                              cruelev_meta,               &
       &                              def_cruelev_meta, &
  ! albedo
       &                              alb_field_mom_meta, &
       &                              alnid_field_mom_meta, &
       &                              aluvd_field_mom_meta, &
       &                              alb_dry_meta, &
       &                              alb_sat_meta, &
       &                              def_alb_meta, &
  ! ndvi
       &                              ndvi_max_meta, &
       &                              ndvi_field_mom_meta, &
       &                              ndvi_ratio_mom_meta, &
       &                              def_ndvi_meta, &
  ! emiss
       &                              def_emiss_meta, & 
       &                              emiss_field_mom_meta, &
       &                              emiss_ratio_mom_meta, &
       &                              emiss_max_meta, &
  ! era
       &                              sst_field_meta, &
       &                              wsnow_field_meta,&
       &                              t2m_field_meta, & 
       &                              hsurf_field_meta,&
       &                              def_era_meta, &
  ! ahf
       &                              ahf_field_meta, &
       &                              def_ahf_meta, &
  ! isa
       &                              def_isa_fields_meta, &
       &                              isa_tot_npixel_meta,isa_field_meta, &
       &                              isa_field_meta


  IMPLICIT NONE

  PRIVATE

  PUBLIC :: &
  ! emiss
       &    read_netcdf_buffer_emiss, &
  ! albedo
       &    read_netcdf_buffer_alb, &
  ! ndvi
       &    read_netcdf_buffer_ndvi, &
  ! cru
       &    read_netcdf_buffer_cru, &
  ! era
       &    read_netcdf_buffer_era, &
  ! ahf
       &    read_netcdf_buffer_ahf, &
  ! isa
       &    read_netcdf_buffer_isa

  CONTAINS

  SUBROUTINE read_netcdf_buffer_emiss(netcdf_filename,  &
       &                              tg,         &
       &                              ntime, &
       &                              emiss_max,  &
       &                              emiss_field_mom,&
       &                              emiss_ratio_mom)


    CHARACTER (len=*), INTENT(IN)      :: netcdf_filename !< filename for the netcdf file
    TYPE(target_grid_def), INTENT(IN)  :: tg !< structure with target grid description
    INTEGER (KIND=i4), INTENT(INOUT)   :: ntime !< number of times of emiss data (12 monthly mean values)
    REAL (KIND=wp), INTENT(OUT)        :: emiss_max(:,:,:), &  !< field for emiss maximum
         &                                emiss_field_mom(:,:,:,:), &  !< field for monthly mean emiss data (12 months)
         &                                emiss_ratio_mom(:,:,:,:) !< field for monthly emiss ratio (12 months)

    CALL logging%info('Enter routine: read_netcdf_buffer_emiss')

    !set up dimensions for buffer
    CALL  def_dimension_info_buffer(tg)
    ! dim_3d_tg
    ! define meta information for target field variables lon_geo, lat_geo 
    CALL def_com_target_fields_meta(dim_3d_tg)
    ! lon_geo_meta and lat_geo_meta
    !define meta information for various EMISS data related variables for netcdf output
    CALL def_emiss_meta(ntime,dim_3d_tg)
    ! dim_emiss_tg, emiss_max_meta, emiss_field_mom_meta, emiss_ratio_mom_meta

    CALL netcdf_get_var(TRIM(netcdf_filename),emiss_max_meta,emiss_max)

    CALL netcdf_get_var(TRIM(netcdf_filename),emiss_field_mom_meta,emiss_field_mom)

    CALL netcdf_get_var(TRIM(netcdf_filename),emiss_ratio_mom_meta,emiss_ratio_mom)

    CALL logging%info('Exit routine: read_netcdf_buffer_emiss')

  END SUBROUTINE read_netcdf_buffer_emiss

  SUBROUTINE read_netcdf_buffer_alb(netcdf_filename,  &
       &                            tg,         &
       &                            ntime, &
       &                            alb_field_mom, &
       &                            alnid_field_mom, &
       &                            aluvd_field_mom, &
       &                            alb_dry, &
       &                            alb_sat)

    CHARACTER (len=*), INTENT(IN)         :: netcdf_filename !< filename for the netcdf file
    TYPE(target_grid_def), INTENT(IN)     :: tg !< structure with target grid description
    INTEGER (KIND=i4), INTENT(OUT)        :: ntime !< number of times of input data (12 monthly mean values)
    REAL (KIND=wp), INTENT(OUT), OPTIONAL :: alb_field_mom(:,:,:,:), & !< field for monthly mean albedo data (12 months)
      &                                      alnid_field_mom(:,:,:,:), &
      &                                      aluvd_field_mom(:,:,:,:), &
      &                                      alb_dry(:,:,:), &
      &                                      alb_sat(:,:,:)

    CALL logging%info('Enter routine: read_netcdf_buffer_alb')
    !set up dimensions for buffer
    CALL  def_dimension_info_buffer(tg)

    ! dim_3d_tg
    ! define meta information for target field variables lon_geo, lat_geo 
    CALL def_com_target_fields_meta(dim_3d_tg)
    ! lon_geo_meta and lat_geo_meta

    !define albedo meta information, related variables for netcdf output
    CALL def_alb_meta(ntime,dim_3d_tg)

    IF (PRESENT(alb_field_mom)) THEN
      CALL netcdf_get_var(TRIM(netcdf_filename),alb_field_mom_meta,alb_field_mom)
    ENDIF
    IF (PRESENT(alnid_field_mom)) THEN
      CALL netcdf_get_var(TRIM(netcdf_filename),alnid_field_mom_meta,alnid_field_mom)
    ENDIF
    IF (PRESENT(aluvd_field_mom)) THEN
      CALL netcdf_get_var(TRIM(netcdf_filename),aluvd_field_mom_meta,aluvd_field_mom)
    ENDIF

    IF (PRESENT(alb_dry)) THEN
      CALL netcdf_get_var(TRIM(netcdf_filename),alb_dry_meta,alb_dry)
    ENDIF
    IF (PRESENT(alb_sat)) THEN
      CALL netcdf_get_var(TRIM(netcdf_filename),alb_sat_meta,alb_sat)
    ENDIF

    CALL logging%info('Exit routine: read_netcdf_buffer_alb')

  END SUBROUTINE read_netcdf_buffer_alb

  SUBROUTINE read_netcdf_buffer_ndvi(netcdf_filename,  &
       &                             tg,         &
       &                             ntime, &
       &                             ndvi_max,  &
       &                             ndvi_field_mom,&
       &                             ndvi_ratio_mom)


    CHARACTER (len=*), INTENT(IN)      :: netcdf_filename !< filename for the netcdf file
    TYPE(target_grid_def), INTENT(IN)  :: tg !< structure with target grid description
    INTEGER (KIND=i4), INTENT(INOUT)   :: ntime !< number of times of ndvi data (12 monthly mean values)

    REAL (KIND=wp), INTENT(OUT)        :: ndvi_max(:,:,:), & !< field for ndvi maximum
         &                                ndvi_field_mom(:,:,:,:), & !< field for monthly mean ndvi data (12 months)
         &                                ndvi_ratio_mom(:,:,:,:) !< field for monthly ndvi ratio (12 months)

    CALL logging%info('Enter routine: read_netcdf_buffer_ndvi')

    !set up dimensions for buffer
    CALL  def_dimension_info_buffer(tg)
    ! dim_3d_tg
    ! define meta information for target field variables lon_geo, lat_geo 
    CALL def_com_target_fields_meta(dim_3d_tg)
    ! lon_geo_meta and lat_geo_meta
    !define meta information for various NDVI data related variables for netcdf output
    CALL def_ndvi_meta(ntime,dim_3d_tg)
    ! dim_ndvi_tg, ndvi_max_meta, ndvi_field_mom_meta, ndvi_ratio_mom_meta

    CALL netcdf_get_var(TRIM(netcdf_filename),ndvi_max_meta,ndvi_max)

    CALL netcdf_get_var(TRIM(netcdf_filename),ndvi_field_mom_meta,ndvi_field_mom)

    CALL netcdf_get_var(TRIM(netcdf_filename),ndvi_ratio_mom_meta,ndvi_ratio_mom)
    
    CALL logging%info('Exit routine: read_netcdf_buffer_ndvi')

  END SUBROUTINE read_netcdf_buffer_ndvi

  SUBROUTINE read_netcdf_buffer_cru(netcdf_filename,  &
       &                            tg,         &
       &                            crutemp,    &
       &                            cruelev)

    CHARACTER (len=*), INTENT(IN)      :: netcdf_filename !< filename for the netcdf file
    TYPE(target_grid_def), INTENT(IN)  :: tg !< structure with target grid description
    REAL(KIND=wp), INTENT(OUT)         :: crutemp(:,:,:)  !< cru climatological temperature , crutemp(ie,je,ke)
    REAL(KIND=wp), OPTIONAL,INTENT(OUT):: cruelev(:,:,:)  !< cru elevation , cruelev(ie,je,ke)

    CALL logging%info('Enter routine: read_netcdf_buffer_cru')

    !set up dimensions for buffer
    CALL  def_dimension_info_buffer(tg)

    ! define meta information for variable crutemp for netcdf output
    CALL def_crutemp_meta(dim_3d_tg)
    ! crutemp_meta

    CALL netcdf_get_var(TRIM(netcdf_filename),crutemp_meta,crutemp)

    IF(PRESENT(cruelev)) THEN
      CALL def_cruelev_meta(dim_3d_tg)
      CALL netcdf_get_var(TRIM(netcdf_filename),cruelev_meta,cruelev)
    ENDIF

    CALL logging%info('Exit routine: read_netcdf_buffer_cru')
  END SUBROUTINE read_netcdf_buffer_cru

  SUBROUTINE read_netcdf_buffer_era(era_buffer_file,  &
       &                            sst_file_legacy, &
       &                            t2m_file_legacy, &
       &                            tg,         &
       &                            ntime, &
       &                            l_use_unified_era, &
       &                            sst_field,&
       &                            wsnow_field, &
       &                            t2m_field, &
       &                            hsurf_field)



    CHARACTER (len=*), INTENT(IN)      :: era_buffer_file, & !< name of unified era-file
         &                                sst_file_legacy, & !< name of legacy sst-file
         &                                t2m_file_legacy !< name of legacy t2m-file

    TYPE(target_grid_def), INTENT(IN)  :: tg !< structure with target grid description
    LOGICAL, INTENT(IN)                :: l_use_unified_era !< use buffer file from extpar_era_to_buffer.py
    INTEGER (KIND=i4), INTENT(INOUT)   :: ntime !< number of times of sst data (12 monthly mean values)

    REAL (KIND=wp), INTENT(OUT)        :: sst_field(:,:,:,:), & !< field for monthly mean sst data (12 months)
         &                                wsnow_field(:,:,:,:), & !< field for monthly sst ratio (12 months)
         &                                t2m_field(:,:,:,:), & !< field for monthly mean t2m data (12 months)
         &                                hsurf_field(:,:,:) !< field for hsurf

    ! local variables
    INTEGER(KIND=i4), PARAMETER        :: nglob_atts=6
    
    CALL logging%info('Enter routine: read_netcdf_buffer_era')
    
    !set up dimensions for buffer
    CALL  def_dimension_info_buffer(tg)
    ! dim_3d_tg
    ! define meta information for target field variables lon_geo, lat_geo 
    CALL def_com_target_fields_meta(dim_3d_tg)
    ! lon_geo_meta and lat_geo_meta
    !define meta information for various SST data related variables for netcdf output
    CALL def_era_meta(ntime,dim_3d_tg)
    ! dim_sst_tg, sst_field_meta, wsnow_field_meta

    IF (l_use_unified_era) THEN
      CALL netcdf_get_var(TRIM(era_buffer_file),sst_field_meta,sst_field)

      CALL netcdf_get_var(TRIM(era_buffer_file),wsnow_field_meta,wsnow_field)

      CALL netcdf_get_var(TRIM(era_buffer_file),t2m_field_meta,t2m_field)

      CALL netcdf_get_var(TRIM(era_buffer_file),hsurf_field_meta,hsurf_field)

    ELSE
      CALL netcdf_get_var(TRIM(sst_file_legacy),sst_field_meta,sst_field)

      CALL netcdf_get_var(TRIM(sst_file_legacy),wsnow_field_meta,wsnow_field)

      CALL netcdf_get_var(TRIM(t2m_file_legacy),t2m_field_meta,t2m_field)

      CALL netcdf_get_var(TRIM(t2m_file_legacy),hsurf_field_meta,hsurf_field)
    ENDIF

    CALL logging%info('Exit routine: read_netcdf_buffer_era')

  END SUBROUTINE read_netcdf_buffer_era

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

  SUBROUTINE read_netcdf_buffer_isa(netcdf_filename,  &
       &                            tg,         &
       &                            isa_field,  &
       &                            isa_tot_npixel)


    CHARACTER (len=*), INTENT(IN)     :: netcdf_filename
    TYPE(target_grid_def), INTENT(IN) :: tg
    INTEGER (KIND=i4), INTENT(OUT)    :: isa_tot_npixel(:,:,:)
    REAL (KIND=wp), INTENT(OUT)       :: isa_field(:,:,:)   !< urban fraction due to isa data

    CALL logging%info('Enter routine: read_netcdf_buffer_isa')

    !set up dimensions for buffer
    CALL  def_dimension_info_buffer(tg)
    ! dim_3d_tg

    ! define meta information for various isa related variables for netcdf output
    CALL def_isa_fields_meta(dim_3d_tg)

    ! define meta information for target field variables lon_geo, lat_geo 
    CALL def_com_target_fields_meta(dim_3d_tg)
    ! lon_geo_meta and lat_geo_meta

    CALL netcdf_get_var(TRIM(netcdf_filename),isa_tot_npixel_meta,isa_tot_npixel)

    CALL netcdf_get_var(TRIM(netcdf_filename),isa_field_meta,isa_field)

    CALL logging%info('Exit routine: read_netcdf_buffer_isa')

  END SUBROUTINE read_netcdf_buffer_isa

END MODULE mo_python_output_nc
