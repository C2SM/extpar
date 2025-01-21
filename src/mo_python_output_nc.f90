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
  ! edgar
       &                              edgar_emi_bc_meta, &
       &                              edgar_emi_oc_meta, &
       &                              edgar_emi_so2_meta, &
       &                              def_edgar_meta, &
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
       &                              isa_field_meta, &
       &                              isa_field_meta, &
  ! hhs      
       &                              hhs_ksat_field_meta, &
       &                              def_hhs_ksat_meta, &
!
       &                              hhs_ormc_field_meta, &
       &                              def_hhs_ormc_meta, &      
!
       &                              hhs_alfa_field_meta, &
       &                              def_hhs_alfa_meta, &   
!
       &                              hhs_critw_field_meta, &
       &                              def_hhs_critw_meta, &   
!
       &                              hhs_fieldc_field_meta, &
       &                              def_hhs_fieldc_meta, &   
!
       &                              hhs_n_field_meta, &
       &                              def_hhs_n_meta, &   
!
       &                              hhs_satf_field_meta, &
       &                              def_hhs_satf_meta, &   
!
       &                              hhs_stc_field_meta, &
       &                              def_hhs_stc_meta, &   
!
       &                              hhs_wcav_field_meta, &
       &                              def_hhs_wcav_meta, &
!
       &                              hhs_wcpf2_field_meta, &
       &                              def_hhs_wcpf2_meta, &   
!
       &                              hhs_wcpf3_field_meta, &
       &                              def_hhs_wcpf3_meta, &
!
       &                              hhs_wcpf42_field_meta, &
       &                              def_hhs_wcpf42_meta, &   
!
       &                              hhs_wcres_field_meta, &
       &                              def_hhs_wcres_meta, &
!
       &                              hhs_wcsat_field_meta, &
       &                              def_hhs_wcsat_meta 
       
       
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: &
  ! emiss
       &    read_netcdf_buffer_emiss, &
  ! albedo
       &    read_netcdf_buffer_alb, &
  ! ndvi
       &    read_netcdf_buffer_ndvi, &
  ! edgar
       &    read_netcdf_buffer_edgar, &
  ! cru
       &    read_netcdf_buffer_cru, &
  ! era
       &    read_netcdf_buffer_era, &
  ! ahf
       &    read_netcdf_buffer_ahf, &
  ! isa
       &    read_netcdf_buffer_isa, &
  ! hhs     
       &    read_netcdf_buffer_hhs_ksat, &
       &    read_netcdf_buffer_hhs_ormc, &
       &    read_netcdf_buffer_hhs_alfa, &
       &    read_netcdf_buffer_hhs_critw, &       
       &    read_netcdf_buffer_hhs_fieldc, &  
       &    read_netcdf_buffer_hhs_n, &
       &    read_netcdf_buffer_hhs_satf, &
       &    read_netcdf_buffer_hhs_stc, &
       &    read_netcdf_buffer_hhs_wcav, &
       &    read_netcdf_buffer_hhs_wcpf2, &       
       &    read_netcdf_buffer_hhs_wcpf3, &  
       &    read_netcdf_buffer_hhs_wcpf42, &
       &    read_netcdf_buffer_hhs_wcres, &  
       &    read_netcdf_buffer_hhs_wcsat
       
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

  SUBROUTINE read_netcdf_buffer_edgar(netcdf_filename,  &
       &                             tg,                &
       &                             edgar_emi_bc,      &
       &                             edgar_emi_oc,      &
       &                             edgar_emi_so2)


    CHARACTER (len=*), INTENT(IN)      :: netcdf_filename !< filename for the netcdf file
    TYPE(target_grid_def), INTENT(IN)  :: tg !< structure with target grid description

    REAL (KIND=wp), INTENT(OUT)        :: edgar_emi_bc(:,:,:), & !< field for black carbon emission from edgar
         &                                edgar_emi_oc(:,:,:), & !< field for organic carbon emission from edgar
         &                                edgar_emi_so2(:,:,:)   !< field for sulfur dioxide emission from edgar

    CALL logging%info('Enter routine: read_netcdf_buffer_edgar')

    !set up dimensions for buffer
    CALL  def_dimension_info_buffer(tg)
    ! dim_3d_tg
    ! define meta information for target field variables lon_geo, lat_geo 
    CALL def_com_target_fields_meta(dim_3d_tg)
    ! lon_geo_meta and lat_geo_meta
    !define meta information for various NDVI data related variables for netcdf output
    CALL def_edgar_meta(dim_3d_tg)
    ! dim_ndvi_tg, ndvi_max_meta, ndvi_field_mom_meta, ndvi_ratio_mom_meta

    CALL netcdf_get_var(TRIM(netcdf_filename),edgar_emi_bc_meta, edgar_emi_bc)

    CALL netcdf_get_var(TRIM(netcdf_filename),edgar_emi_oc_meta, edgar_emi_oc)

    CALL netcdf_get_var(TRIM(netcdf_filename),edgar_emi_so2_meta,edgar_emi_so2)
    
    CALL logging%info('Exit routine: read_netcdf_buffer_edgar')

  END SUBROUTINE read_netcdf_buffer_edgar

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
       &                            isa_field)


    CHARACTER (len=*), INTENT(IN)     :: netcdf_filename
    TYPE(target_grid_def), INTENT(IN) :: tg
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

    CALL netcdf_get_var(TRIM(netcdf_filename),isa_field_meta,isa_field)

    CALL logging%info('Exit routine: read_netcdf_buffer_isa')

  END SUBROUTINE read_netcdf_buffer_isa

  SUBROUTINE read_netcdf_buffer_hhs_ksat(netcdf_filename,  &
   &                                     tg,         &
   &                                     hhs_ksat_field)

    CHARACTER (len=*), INTENT(IN)      :: netcdf_filename !< filename for the netcdf file
    TYPE(target_grid_def), INTENT(IN)  :: tg !< structure with target grid description
    REAL (KIND=wp), INTENT(OUT)        :: hhs_ksat_field(:,:,:) !< field for hhs_ksat 

    ! local variables
    INTEGER(KIND=i4), PARAMETER        :: nglob_atts=6

    CALL logging%info('Enter routine: read_netcdf_buffer_hhs_ksat')

    !set up dimensions for buffer
    CALL  def_dimension_info_buffer(tg)

    ! define meta information for target field variables lon_geo, lat_geo 
    CALL def_com_target_fields_meta(dim_3d_tg)
    ! lon_geo_meta and lat_geo_meta
    !define meta information for various HHS_KSAT data related variables for netcdf output
    CALL def_hhs_ksat_meta(dim_3d_tg)
    ! dim_hhs_ksat_tg, hhs_ksat_field_meta

    CALL netcdf_get_var(TRIM(netcdf_filename),hhs_ksat_field_meta,hhs_ksat_field)

    CALL logging%info('Exit routine: read_netcdf_buffer_hhs_ksat')

  END SUBROUTINE read_netcdf_buffer_hhs_ksat

  SUBROUTINE read_netcdf_buffer_hhs_ormc(netcdf_filename,  &
   &                                     tg,         &
   &                                     hhs_ormc_field)

    CHARACTER (len=*), INTENT(IN)      :: netcdf_filename !< filename for the netcdf file
    TYPE(target_grid_def), INTENT(IN)  :: tg !< structure with target grid description
    REAL (KIND=wp), INTENT(OUT)        :: hhs_ormc_field(:,:,:) !< field for hhs_ormc 

    ! local variables
    INTEGER(KIND=i4), PARAMETER        :: nglob_atts=6

    CALL logging%info('Enter routine: read_netcdf_buffer_hhs_ormc')

    !set up dimensions for buffer
    CALL  def_dimension_info_buffer(tg)

    ! define meta information for target field variables lon_geo, lat_geo 
    CALL def_com_target_fields_meta(dim_3d_tg)
    ! lon_geo_meta and lat_geo_meta
    !define meta information for various HHS_ORMC data related variables for netcdf output
    CALL def_hhs_ormc_meta(dim_3d_tg)
    ! dim_hhs_ormc_tg, hhs_ormc_field_meta

    CALL netcdf_get_var(TRIM(netcdf_filename),hhs_ormc_field_meta,hhs_ormc_field)

    CALL logging%info('Exit routine: read_netcdf_buffer_hhs_ormc')

  END SUBROUTINE read_netcdf_buffer_hhs_ormc
  
  SUBROUTINE read_netcdf_buffer_hhs_alfa(netcdf_filename,  &
   &                                     tg,         &
   &                                     hhs_alfa_field)

    CHARACTER (len=*), INTENT(IN)      :: netcdf_filename !< filename for the netcdf file
    TYPE(target_grid_def), INTENT(IN)  :: tg !< structure with target grid description
    REAL (KIND=wp), INTENT(OUT)        :: hhs_alfa_field(:,:,:) !< field for hhs_alfa 

    ! local variables
    INTEGER(KIND=i4), PARAMETER        :: nglob_atts=6

    CALL logging%info('Enter routine: read_netcdf_buffer_hhs_alfa')

    !set up dimensions for buffer
    CALL  def_dimension_info_buffer(tg)

    ! define meta information for target field variables lon_geo, lat_geo 
    CALL def_com_target_fields_meta(dim_3d_tg)
    ! lon_geo_meta and lat_geo_meta
    !define meta information for various HHS_ALFA data related variables for netcdf output
    CALL def_hhs_alfa_meta(dim_3d_tg)
    ! dim_hhs_alfa_tg, hhs_alfa_field_meta

    CALL netcdf_get_var(TRIM(netcdf_filename),hhs_alfa_field_meta,hhs_alfa_field)

    CALL logging%info('Exit routine: read_netcdf_buffer_hhs_alfa')

  END SUBROUTINE read_netcdf_buffer_hhs_alfa

    SUBROUTINE read_netcdf_buffer_hhs_critw(netcdf_filename,  &
   &                                     tg,         &
   &                                     hhs_critw_field)

    CHARACTER (len=*), INTENT(IN)      :: netcdf_filename !< filename for the netcdf file
    TYPE(target_grid_def), INTENT(IN)  :: tg !< structure with target grid description
    REAL (KIND=wp), INTENT(OUT)        :: hhs_critw_field(:,:,:) !< field for hhs_critw 

    ! local variables
    INTEGER(KIND=i4), PARAMETER        :: nglob_atts=6

    CALL logging%info('Enter routine: read_netcdf_buffer_hhs_critw')

    !set up dimensions for buffer
    CALL  def_dimension_info_buffer(tg)

    ! define meta information for target field variables lon_geo, lat_geo 
    CALL def_com_target_fields_meta(dim_3d_tg)
    ! lon_geo_meta and lat_geo_meta
    !define meta information for various HHS_CRITW data related variables for netcdf output
    CALL def_hhs_critw_meta(dim_3d_tg)
    ! dim_hhs_critw_tg, hhs_critw_field_meta

    CALL netcdf_get_var(TRIM(netcdf_filename),hhs_critw_field_meta,hhs_critw_field)

    CALL logging%info('Exit routine: read_netcdf_buffer_hhs_critw')

  END SUBROUTINE read_netcdf_buffer_hhs_critw

    SUBROUTINE read_netcdf_buffer_hhs_fieldc(netcdf_filename,  &
   &                                     tg,         &
   &                                     hhs_fieldc_field)

    CHARACTER (len=*), INTENT(IN)      :: netcdf_filename !< filename for the netcdf file
    TYPE(target_grid_def), INTENT(IN)  :: tg !< structure with target grid description
    REAL (KIND=wp), INTENT(OUT)        :: hhs_fieldc_field(:,:,:) !< field for hhs_fieldc 

    ! local variables
    INTEGER(KIND=i4), PARAMETER        :: nglob_atts=6

    CALL logging%info('Enter routine: read_netcdf_buffer_hhs_fieldc')

    !set up dimensions for buffer
    CALL  def_dimension_info_buffer(tg)

    ! define meta information for target field variables lon_geo, lat_geo 
    CALL def_com_target_fields_meta(dim_3d_tg)
    ! lon_geo_meta and lat_geo_meta
    !define meta information for various HHS_FIELDC data related variables for netcdf output
    CALL def_hhs_fieldc_meta(dim_3d_tg)
    ! dim_hhs_fieldc_tg, hhs_fieldc_field_meta

    CALL netcdf_get_var(TRIM(netcdf_filename),hhs_fieldc_field_meta,hhs_fieldc_field)

    CALL logging%info('Exit routine: read_netcdf_buffer_hhs_fieldc')

  END SUBROUTINE read_netcdf_buffer_hhs_fieldc

    SUBROUTINE read_netcdf_buffer_hhs_n(netcdf_filename,  &
   &                                     tg,         &
   &                                     hhs_n_field)

    CHARACTER (len=*), INTENT(IN)      :: netcdf_filename !< filename for the netcdf file
    TYPE(target_grid_def), INTENT(IN)  :: tg !< structure with target grid description
    REAL (KIND=wp), INTENT(OUT)        :: hhs_n_field(:,:,:) !< field for hhs_n 

    ! local variables
    INTEGER(KIND=i4), PARAMETER        :: nglob_atts=6

    CALL logging%info('Enter routine: read_netcdf_buffer_hhs_n')

    !set up dimensions for buffer
    CALL  def_dimension_info_buffer(tg)

    ! define meta information for target field variables lon_geo, lat_geo 
    CALL def_com_target_fields_meta(dim_3d_tg)
    ! lon_geo_meta and lat_geo_meta
    !define meta information for various HHS_N data related variables for netcdf output
    CALL def_hhs_n_meta(dim_3d_tg)
    ! dim_hhs_n_tg, hhs_n_field_meta

    CALL netcdf_get_var(TRIM(netcdf_filename),hhs_n_field_meta,hhs_n_field)

    CALL logging%info('Exit routine: read_netcdf_buffer_hhs_n')

  END SUBROUTINE read_netcdf_buffer_hhs_n

    SUBROUTINE read_netcdf_buffer_hhs_satf(netcdf_filename,  &
   &                                     tg,         &
   &                                     hhs_satf_field)

    CHARACTER (len=*), INTENT(IN)      :: netcdf_filename !< filename for the netcdf file
    TYPE(target_grid_def), INTENT(IN)  :: tg !< structure with target grid description
    REAL (KIND=wp), INTENT(OUT)        :: hhs_satf_field(:,:,:) !< field for hhs_satf 

    ! local variables
    INTEGER(KIND=i4), PARAMETER        :: nglob_atts=6

    CALL logging%info('Enter routine: read_netcdf_buffer_hhs_satf')

    !set up dimensions for buffer
    CALL  def_dimension_info_buffer(tg)

    ! define meta information for target field variables lon_geo, lat_geo 
    CALL def_com_target_fields_meta(dim_3d_tg)
    ! lon_geo_meta and lat_geo_meta
    !define meta information for various HHS_SATF data related variables for netcdf output
    CALL def_hhs_satf_meta(dim_3d_tg)
    ! dim_hhs_satf_tg, hhs_satf_field_meta

    CALL netcdf_get_var(TRIM(netcdf_filename),hhs_satf_field_meta,hhs_satf_field)

    CALL logging%info('Exit routine: read_netcdf_buffer_hhs_satf')

  END SUBROUTINE read_netcdf_buffer_hhs_satf

    SUBROUTINE read_netcdf_buffer_hhs_stc(netcdf_filename,  &
   &                                     tg,         &
   &                                     hhs_stc_field)

    CHARACTER (len=*), INTENT(IN)      :: netcdf_filename !< filename for the netcdf file
    TYPE(target_grid_def), INTENT(IN)  :: tg !< structure with target grid description
    REAL (KIND=wp), INTENT(OUT)        :: hhs_stc_field(:,:,:) !< field for hhs_stc 

    ! local variables
    INTEGER(KIND=i4), PARAMETER        :: nglob_atts=6

    CALL logging%info('Enter routine: read_netcdf_buffer_hhs_stc')

    !set up dimensions for buffer
    CALL  def_dimension_info_buffer(tg)

    ! define meta information for target field variables lon_geo, lat_geo 
    CALL def_com_target_fields_meta(dim_3d_tg)
    ! lon_geo_meta and lat_geo_meta
    !define meta information for various HHS_STC data related variables for netcdf output
    CALL def_hhs_stc_meta(dim_3d_tg)
    ! dim_hhs_stc_tg, hhs_stc_field_meta

    CALL netcdf_get_var(TRIM(netcdf_filename),hhs_stc_field_meta,hhs_stc_field)

    CALL logging%info('Exit routine: read_netcdf_buffer_hhs_stc')

  END SUBROUTINE read_netcdf_buffer_hhs_stc

    SUBROUTINE read_netcdf_buffer_hhs_wcav(netcdf_filename,  &
   &                                     tg,         &
   &                                     hhs_wcav_field)

    CHARACTER (len=*), INTENT(IN)      :: netcdf_filename !< filename for the netcdf file
    TYPE(target_grid_def), INTENT(IN)  :: tg !< structure with target grid description
    REAL (KIND=wp), INTENT(OUT)        :: hhs_wcav_field(:,:,:) !< field for hhs_wcav 

    ! local variables
    INTEGER(KIND=i4), PARAMETER        :: nglob_atts=6

    CALL logging%info('Enter routine: read_netcdf_buffer_hhs_wcav')

    !set up dimensions for buffer
    CALL  def_dimension_info_buffer(tg)

    ! define meta information for target field variables lon_geo, lat_geo 
    CALL def_com_target_fields_meta(dim_3d_tg)
    ! lon_geo_meta and lat_geo_meta
    !define meta information for various HHS_WCAV data related variables for netcdf output
    CALL def_hhs_wcav_meta(dim_3d_tg)
    ! dim_hhs_wcav_tg, hhs_wcav_field_meta

    CALL netcdf_get_var(TRIM(netcdf_filename),hhs_wcav_field_meta,hhs_wcav_field)

    CALL logging%info('Exit routine: read_netcdf_buffer_hhs_wcav')

  END SUBROUTINE read_netcdf_buffer_hhs_wcav

    SUBROUTINE read_netcdf_buffer_hhs_wcpf2(netcdf_filename,  &
   &                                     tg,         &
   &                                     hhs_wcpf2_field)

    CHARACTER (len=*), INTENT(IN)      :: netcdf_filename !< filename for the netcdf file
    TYPE(target_grid_def), INTENT(IN)  :: tg !< structure with target grid description
    REAL (KIND=wp), INTENT(OUT)        :: hhs_wcpf2_field(:,:,:) !< field for hhs_wcpf2 

    ! local variables
    INTEGER(KIND=i4), PARAMETER        :: nglob_atts=6

    CALL logging%info('Enter routine: read_netcdf_buffer_hhs_wcpf2')

    !set up dimensions for buffer
    CALL  def_dimension_info_buffer(tg)

    ! define meta information for target field variables lon_geo, lat_geo 
    CALL def_com_target_fields_meta(dim_3d_tg)
    ! lon_geo_meta and lat_geo_meta
    !define meta information for various HHS_WCPF2 data related variables for netcdf output
    CALL def_hhs_wcpf2_meta(dim_3d_tg)
    ! dim_hhs_wcpf2_tg, hhs_wcpf2_field_meta

    CALL netcdf_get_var(TRIM(netcdf_filename),hhs_wcpf2_field_meta,hhs_wcpf2_field)

    CALL logging%info('Exit routine: read_netcdf_buffer_hhs_wcpf2')

  END SUBROUTINE read_netcdf_buffer_hhs_wcpf2

      SUBROUTINE read_netcdf_buffer_hhs_wcpf3(netcdf_filename,  &
   &                                     tg,         &
   &                                     hhs_wcpf3_field)

    CHARACTER (len=*), INTENT(IN)      :: netcdf_filename !< filename for the netcdf file
    TYPE(target_grid_def), INTENT(IN)  :: tg !< structure with target grid description
    REAL (KIND=wp), INTENT(OUT)        :: hhs_wcpf3_field(:,:,:) !< field for hhs_wcpf3 

    ! local variables
    INTEGER(KIND=i4), PARAMETER        :: nglob_atts=6

    CALL logging%info('Enter routine: read_netcdf_buffer_hhs_wcpf3')

    !set up dimensions for buffer
    CALL  def_dimension_info_buffer(tg)

    ! define meta information for target field variables lon_geo, lat_geo 
    CALL def_com_target_fields_meta(dim_3d_tg)
    ! lon_geo_meta and lat_geo_meta
    !define meta information for various HHS_WCPF3 data related variables for netcdf output
    CALL def_hhs_wcpf3_meta(dim_3d_tg)
    ! dim_hhs_wcpf3_tg, hhs_wcpf3_field_meta

    CALL netcdf_get_var(TRIM(netcdf_filename),hhs_wcpf3_field_meta,hhs_wcpf3_field)

    CALL logging%info('Exit routine: read_netcdf_buffer_hhs_wcpf3')

  END SUBROUTINE read_netcdf_buffer_hhs_wcpf3

      SUBROUTINE read_netcdf_buffer_hhs_wcpf42(netcdf_filename,  &
   &                                     tg,         &
   &                                     hhs_wcpf42_field)

    CHARACTER (len=*), INTENT(IN)      :: netcdf_filename !< filename for the netcdf file
    TYPE(target_grid_def), INTENT(IN)  :: tg !< structure with target grid description
    REAL (KIND=wp), INTENT(OUT)        :: hhs_wcpf42_field(:,:,:) !< field for hhs_wcpf42 

    ! local variables
    INTEGER(KIND=i4), PARAMETER        :: nglob_atts=6

    CALL logging%info('Enter routine: read_netcdf_buffer_hhs_wcpf42')

    !set up dimensions for buffer
    CALL  def_dimension_info_buffer(tg)

    ! define meta information for target field variables lon_geo, lat_geo 
    CALL def_com_target_fields_meta(dim_3d_tg)
    ! lon_geo_meta and lat_geo_meta
    !define meta information for various HHS_WCPF42 data related variables for netcdf output
    CALL def_hhs_wcpf42_meta(dim_3d_tg)
    ! dim_hhs_wcpf42_tg, hhs_wcpf42_field_meta

    CALL netcdf_get_var(TRIM(netcdf_filename),hhs_wcpf42_field_meta,hhs_wcpf42_field)

    CALL logging%info('Exit routine: read_netcdf_buffer_hhs_wcpf42')

  END SUBROUTINE read_netcdf_buffer_hhs_wcpf42

      SUBROUTINE read_netcdf_buffer_hhs_wcres(netcdf_filename,  &
   &                                     tg,         &
   &                                     hhs_wcres_field)

    CHARACTER (len=*), INTENT(IN)      :: netcdf_filename !< filename for the netcdf file
    TYPE(target_grid_def), INTENT(IN)  :: tg !< structure with target grid description
    REAL (KIND=wp), INTENT(OUT)        :: hhs_wcres_field(:,:,:) !< field for hhs_wcres 

    ! local variables
    INTEGER(KIND=i4), PARAMETER        :: nglob_atts=6

    CALL logging%info('Enter routine: read_netcdf_buffer_hhs_wcres')

    !set up dimensions for buffer
    CALL  def_dimension_info_buffer(tg)

    ! define meta information for target field variables lon_geo, lat_geo 
    CALL def_com_target_fields_meta(dim_3d_tg)
    ! lon_geo_meta and lat_geo_meta
    !define meta information for various HHS_WCRES data related variables for netcdf output
    CALL def_hhs_wcres_meta(dim_3d_tg)
    ! dim_hhs_wcres_tg, hhs_wcres_field_meta

    CALL netcdf_get_var(TRIM(netcdf_filename),hhs_wcres_field_meta,hhs_wcres_field)

    CALL logging%info('Exit routine: read_netcdf_buffer_hhs_wcres')

  END SUBROUTINE read_netcdf_buffer_hhs_wcres

        SUBROUTINE read_netcdf_buffer_hhs_wcsat(netcdf_filename,  &
   &                                     tg,         &
   &                                     hhs_wcsat_field)

    CHARACTER (len=*), INTENT(IN)      :: netcdf_filename !< filename for the netcdf file
    TYPE(target_grid_def), INTENT(IN)  :: tg !< structure with target grid description
    REAL (KIND=wp), INTENT(OUT)        :: hhs_wcsat_field(:,:,:) !< field for hhs_wcsat 

    ! local variables
    INTEGER(KIND=i4), PARAMETER        :: nglob_atts=6

    CALL logging%info('Enter routine: read_netcdf_buffer_hhs_wcsat')

    !set up dimensions for buffer
    CALL  def_dimension_info_buffer(tg)

    ! define meta information for target field variables lon_geo, lat_geo 
    CALL def_com_target_fields_meta(dim_3d_tg)
    ! lon_geo_meta and lat_geo_meta
    !define meta information for various HHS_WCSAT data related variables for netcdf output
    CALL def_hhs_wcsat_meta(dim_3d_tg)
    ! dim_hhs_wcsat_tg, hhs_wcsat_field_meta

    CALL netcdf_get_var(TRIM(netcdf_filename),hhs_wcsat_field_meta,hhs_wcsat_field)

    CALL logging%info('Exit routine: read_netcdf_buffer_hhs_wcsat')

  END SUBROUTINE read_netcdf_buffer_hhs_wcsat
  
END MODULE mo_python_output_nc
