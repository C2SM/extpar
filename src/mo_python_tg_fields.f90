MODULE mo_python_tg_fields

  USE mo_logging
  USE mo_kind,                  ONLY: wp, i4
  USE mo_array_cache,           ONLY: allocate_cached
  USE mo_grid_structures,       ONLY: target_grid_def
  USE mo_io_utilities,          ONLY: var_meta_info

  IMPLICIT NONE

  PRIVATE

  PUBLIC ::  &
  ! emiss
    &        emiss_field, &
    &        emiss_max, &
    &        emiss_field_mom, &
    &        emiss_ratio_mom, &
    &        allocate_emiss_target_fields, &
  ! ndvi
    &        ndvi_field, &
    &        ndvi_max, &
    &        ndvi_field_mom, &
    &        ndvi_ratio_mom, &
    &        allocate_ndvi_target_fields, &
  ! edgar
    &        edgar_emi_bc, &
    &        edgar_emi_oc, &
    &        edgar_emi_so2, &
    &        edgar_emi_nox, &
    &        edgar_emi_nh3, &
    &        allocate_edgar_target_fields, &
  ! cdnc
    &        cdnc,                        &
    &        allocate_cdnc_target_fields, &
  ! cru
    &        allocate_cru_target_fields, &
    &        crutemp,crutemp2, cruelev,  &
    &        meta_crutemp, meta_cruelev, &
  ! albedo
    &        alb_dry, &
    &        alb_sat, &
    &        alb_field_mom, &
    &        alnid_field_mom, &
    &        aluvd_field_mom, &
    &        allocate_alb_target_fields, &
    &        alb_interpol, &
  ! era
             sst_field, &
    &        wsnow_field, &
    &        t2m_field, &
    &        hsurf_field, &
    &        allocate_era_target_fields, &
  ! ahf
    &        allocate_ahf_target_fields, &
    &        ahf_field, &
  ! isa      
    &        allocate_isa_target_fields, &
    &        isa_field


  REAL(KIND=wp), POINTER :: &
  ! emiss
       &                    emiss_field(:,:,:), & !< field for emiss data
       &                    emiss_max(:,:,:), & !< field for emiss maximum
       &                    emiss_field_mom(:,:,:,:), & !< field for monthly mean emiss data (12 months)
       &                    emiss_ratio_mom(:,:,:,:), & !< field for monthly emiss ratio (12 months)
  ! ndvi
       &                    ndvi_field(:,:,:), & !< field for ndvi data
       &                    ndvi_max(:,:,:), & !< field for ndvi maximum
       &                    ndvi_field_mom(:,:,:,:), & !< field for monthly mean ndvi data (12 months)
       &                    ndvi_ratio_mom(:,:,:,:), & !< field for monthly ndvi ratio (12 months)
  ! edgar
       &                    edgar_emi_bc(:,:,:), & !< field for black carbon emission from edgar
       &                    edgar_emi_oc(:,:,:), & !< field for organic carbon emission from edgar
       &                    edgar_emi_so2(:,:,:), & !< field for sulfur dioxide emission from edgar
       &                    edgar_emi_nox(:,:,:), & !< field for nitrogen oxides emission from edgar
       &                    edgar_emi_nh3(:,:,:), & !< field for ammonia emission from edgar
  ! cdnc
       &                    cdnc(:,:,:,:), & !< field for cloud droplet number (12 months)
  ! cru
       &                    crutemp(:,:,:), & !< cru climatological temperature , crutemp(ie,je,ke)
       &                    crutemp2(:,:,:), & !< cru climatological temperature , crutemp(ie,je,ke)
       &                    cruelev(:,:,:), & !< cru climatological temperature , cruelev(ie,je,ke)
  ! albedo
       &                    alb_field_mom(:,:,:,:), & !< field for monthly mean albedo data (12 months)
       &                    alnid_field_mom(:,:,:,:), &
       &                    aluvd_field_mom(:,:,:,:), &
       &                    alb_interpol(:,:,:,:), & !<  field for interpolated albedo
       &                    alb_dry(:,:,:), & !< field for dry soil albedo
       &                    alb_sat(:,:,:), & !< field for saturated soil albedo
  ! era
       &                    sst_field(:,:,:,:), & !< field for sst data (12 months)
       &                    wsnow_field(:,:,:,:), & !< field for wsnow data (12 months)
       &                    t2m_field(:,:,:,:), & !< field for wsnow data (12 months)
       &                    hsurf_field(:,:,:), & !< field for wsnow data (12 months)
  ! ahf
       &                    ahf_field(:,:,:), & !< fields for artifical heat flux (12 months)
  ! isa
       &                    isa_field(:,:,:) !< fraction land due to land use raw data

  TYPE(var_meta_info)    :: meta_crutemp, meta_cruelev

  CONTAINS

  !> allocate fields for GLOBE target data
  SUBROUTINE allocate_emiss_target_fields(tg, nt, l_use_array_cache)

    TYPE(target_grid_def), INTENT(IN) :: tg  !< structure with target grid description
    INTEGER (KIND=i4), INTENT(in)     :: nt !< number of timesteps (12 for monthly mean values)
    LOGICAL, INTENT(in)               :: l_use_array_cache
    
    INTEGER(KIND=i4)                  :: errorcode !< error status variable

    errorcode = 0
    
    IF (l_use_array_cache) THEN
       CALL allocate_cached('emiss_field', emiss_field, [tg%ie,tg%je,tg%ke])
    ELSE
       allocate(emiss_field(tg%ie,tg%je,tg%ke), stat=errorcode)
    ENDIF
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array emiss_field',__FILE__,__LINE__)
    emiss_field = 0.0

    IF (l_use_array_cache) THEN
       CALL allocate_cached('emiss_max', emiss_max, [tg%ie,tg%je,tg%ke])
    ELSE
       allocate(emiss_max(tg%ie,tg%je,tg%ke), stat=errorcode)
    ENDIF
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array emiss_max',__FILE__,__LINE__)
    emiss_max = 0.0

    IF (l_use_array_cache) THEN
       CALL allocate_cached('emiss_field_mom', emiss_field_mom, [tg%ie,tg%je,tg%ke,nt])
    ELSE
       allocate(emiss_field_mom(tg%ie,tg%je,tg%ke,nt), stat=errorcode)
    ENDIF
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array emiss_field_mom',__FILE__,__LINE__)
    emiss_field_mom = 0.0

    IF (l_use_array_cache) THEN
       CALL allocate_cached('emiss_ratio_mom', emiss_ratio_mom, [tg%ie,tg%je,tg%ke,nt])
    ELSE
       allocate(emiss_ratio_mom(tg%ie,tg%je,tg%ke,nt), stat=errorcode)
    ENDIF
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array emiss_ratio_mom',__FILE__,__LINE__)
    emiss_ratio_mom = 0.0

  END SUBROUTINE allocate_emiss_target_fields

  SUBROUTINE allocate_ndvi_target_fields(tg,nt, l_use_array_cache)

    TYPE(target_grid_def), INTENT(IN) :: tg  !< structure with target grid description
    INTEGER (KIND=i4), INTENT(IN)     :: nt !< number of timesteps (12 for monthly mean values)
    LOGICAL, INTENT(in)               :: l_use_array_cache
    
    INTEGER(KIND=i4)                  :: errorcode !< error status variable

    errorcode = 0
    
    CALL logging%info('Enter routine: allocate_ndvi_target_fields')

    IF (l_use_array_cache) THEN
       call allocate_cached('ndvi_field', ndvi_field, [tg%ie,tg%je,tg%ke])
    ELSE
       allocate(ndvi_field(tg%ie,tg%je,tg%ke), stat=errorcode)
    ENDIF
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array ndvi_field',__FILE__,__LINE__)
    ndvi_field = 0.0

    IF (l_use_array_cache) THEN
       call allocate_cached('ndvi_max', ndvi_max, [tg%ie,tg%je,tg%ke])
    ELSE
       allocate(ndvi_max(tg%ie,tg%je,tg%ke), stat=errorcode)
    ENDIF
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array ndvi_max',__FILE__,__LINE__)
    ndvi_max = 0.0

    IF (l_use_array_cache) THEN
       call allocate_cached('ndvi_field_mom', ndvi_field_mom, [tg%ie,tg%je,tg%ke,nt])
    ELSE
       allocate(ndvi_field_mom(tg%ie,tg%je,tg%ke,nt), stat=errorcode)
    ENDIF
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array ndvi_field_mom',__FILE__,__LINE__)
    ndvi_field_mom = 0.0

    IF (l_use_array_cache) THEN
       call allocate_cached('ndvi_ratio_mom', ndvi_ratio_mom, [tg%ie,tg%je,tg%ke,nt])
    ELSE
       allocate(ndvi_ratio_mom(tg%ie,tg%je,tg%ke,nt), stat=errorcode)
    ENDIF
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array ndvi_ratio_mom',__FILE__,__LINE__)
    ndvi_ratio_mom = 0.0

  END SUBROUTINE allocate_ndvi_target_fields

  SUBROUTINE allocate_edgar_target_fields(tg, l_use_array_cache)

    TYPE(target_grid_def), INTENT(IN) :: tg  !< structure with target grid description
    LOGICAL, INTENT(in)               :: l_use_array_cache
    
    INTEGER(KIND=i4)                  :: errorcode !< error status variable

    errorcode = 0

    CALL logging%info('Enter routine: allocate_edgar_target_fields')

    IF (l_use_array_cache) THEN
       call allocate_cached('emi_bc', edgar_emi_bc, [tg%ie,tg%je,tg%ke])
    ELSE
       allocate(edgar_emi_bc(tg%ie,tg%je,tg%ke), stat=errorcode)
    ENDIF
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array edgar_emi_bc',__FILE__,__LINE__)
    edgar_emi_bc = 0.0

    IF (l_use_array_cache) THEN
       call allocate_cached('emi_oc', edgar_emi_oc, [tg%ie,tg%je,tg%ke])
    ELSE
       allocate(edgar_emi_oc(tg%ie,tg%je,tg%ke), stat=errorcode)
    ENDIF
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array edgar_emi_oc',__FILE__,__LINE__)
    edgar_emi_oc = 0.0

    IF (l_use_array_cache) THEN
       call allocate_cached('emi_so2', edgar_emi_so2, [tg%ie,tg%je,tg%ke])
    ELSE
       allocate(edgar_emi_so2(tg%ie,tg%je,tg%ke), stat=errorcode)
    ENDIF
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array edgar_emi_so2',__FILE__,__LINE__)
    edgar_emi_so2 = 0.0

    IF (l_use_array_cache) THEN
       call allocate_cached('emi_nox', edgar_emi_nox, [tg%ie,tg%je,tg%ke])
    ELSE
       allocate(edgar_emi_nox(tg%ie,tg%je,tg%ke), stat=errorcode)
    ENDIF
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array edgar_emi_nox',__FILE__,__LINE__)
    edgar_emi_nox = 0.0

    IF (l_use_array_cache) THEN
       call allocate_cached('emi_nh3', edgar_emi_nh3, [tg%ie,tg%je,tg%ke])
    ELSE
       allocate(edgar_emi_nh3(tg%ie,tg%je,tg%ke), stat=errorcode)
    ENDIF
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array edgar_emi_nh3',__FILE__,__LINE__)
    edgar_emi_nh3 = 0.0

  END SUBROUTINE allocate_edgar_target_fields

  SUBROUTINE allocate_cdnc_target_fields(tg, nt, l_use_array_cache)

    TYPE(target_grid_def), INTENT(IN) :: tg  !< structure with target grid description
    INTEGER (KIND=i4), INTENT(IN)     :: nt  !< number of timesteps (12 for monthly mean values)
    LOGICAL, INTENT(in)               :: l_use_array_cache
    
    INTEGER(KIND=i4)                  :: errorcode !< error status variable

    errorcode = 0

    CALL logging%info('Enter routine: allocate_cdnc_target_fields')

    IF (l_use_array_cache) THEN
       call allocate_cached('cdnc', cdnc, [tg%ie,tg%je,tg%ke,nt])
    ELSE
       allocate(cdnc(tg%ie,tg%je,tg%ke,nt), stat=errorcode)
    ENDIF
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array cdnc',__FILE__,__LINE__)
    cdnc = 0.0

  END SUBROUTINE allocate_cdnc_target_fields

  SUBROUTINE allocate_cru_target_fields(tg, l_use_array_cache)

    TYPE(target_grid_def), INTENT(IN) :: tg  !< structure with target grid description
    LOGICAL, INTENT(in)               :: l_use_array_cache

    INTEGER                           :: errorcode !< error status variable

    errorcode = 0
    
    CALL logging%info('Enter routine: allocate_cru_target_fields')

    IF (l_use_array_cache) then
       CALL allocate_cached('crutemp', crutemp, [tg%ie,tg%je,tg%ke])
    ELSE
       ALLOCATE(crutemp(tg%ie,tg%je,tg%ke), stat=errorcode)
    ENDIF

    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array crutemp',__FILE__,__LINE__)
    crutemp = 0.0

    IF (l_use_array_cache) then
       CALL allocate_cached('crutemp2', crutemp2, [tg%ie,tg%je,tg%ke])
    ELSE
       ALLOCATE(crutemp2(tg%ie,tg%je,tg%ke), stat=errorcode)
    ENDIF

    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array crutemp',__FILE__,__LINE__)
    crutemp2 = 0.0

    IF (l_use_array_cache) then
       CALL allocate_cached('cruelev', cruelev, [tg%ie,tg%je,tg%ke])
    ELSE
       ALLOCATE(cruelev(tg%ie,tg%je,tg%ke), stat=errorcode)
    ENDIF

    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array cruelev',__FILE__,__LINE__)
    cruelev = 0.0

    meta_crutemp%varname = 'tem_clim'
    meta_crutemp%n_dim = 3

    ALLOCATE(meta_crutemp%diminfo(meta_crutemp%n_dim), stat=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array meta_crutemp%diminfo',__FILE__,__LINE__)

    meta_crutemp%diminfo(1)%dimname = 'ie'
    meta_crutemp%diminfo(1)%dimsize = tg%ie
    meta_crutemp%diminfo(2)%dimname = 'je'
    meta_crutemp%diminfo(2)%dimsize = tg%je
    meta_crutemp%diminfo(3)%dimname = 'ke'
    meta_crutemp%diminfo(3)%dimsize = tg%ke
    meta_crutemp%vartype = 2 ! REAL variable
    meta_crutemp%standard_name = 'CRU T'
    meta_crutemp%long_name = 'CRU near surface temperature climatology'
    meta_crutemp%units = ''

    meta_cruelev%varname = 'elev_clim'
    meta_cruelev%n_dim = 3

    ALLOCATE(meta_cruelev%diminfo(meta_cruelev%n_dim), stat=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array meta_cruelev%diminfo',__FILE__,__LINE__)

    meta_cruelev%diminfo(1)%dimname = 'ie'
    meta_cruelev%diminfo(1)%dimsize = tg%ie
    meta_cruelev%diminfo(2)%dimname = 'je'
    meta_cruelev%diminfo(2)%dimsize = tg%je
    meta_cruelev%diminfo(3)%dimname = 'ke'
    meta_cruelev%diminfo(3)%dimsize = tg%ke
    meta_cruelev%vartype = 2 ! REAL variable
    meta_cruelev%standard_name = 'HSURF'
    meta_cruelev%long_name = 'CRU grid cell elevation'
    meta_cruelev%units = 'm'

  END SUBROUTINE allocate_cru_target_fields

  SUBROUTINE allocate_alb_target_fields(tg,nt,raw_id, l_use_array_cache)

    TYPE(target_grid_def), INTENT(IN) :: tg  !< structure with target grid description
    INTEGER (KIND=i4), INTENT(IN)     :: nt, & !< number of timesteps (12 for monthly mean values)
      &                                  raw_id !< type of albedo treatment
    LOGICAL, INTENT(in)               :: l_use_array_cache
    
    INTEGER(KIND=i4)                  :: errorcode !< error status variable

    errorcode = 0
    
    CALL logging%info('Enter routine: allocate_alb_target_fields')

    IF (l_use_array_cache) THEN
       CALL allocate_cached('alb_field_mom', alb_field_mom, [tg%ie,tg%je,tg%ke,nt])
    ELSE
       allocate(alb_field_mom(tg%ie,tg%je,tg%ke,nt), stat=errorcode)
    ENDIF
    IF(errorcode.NE.0) CALL logging%error('Cant allocate array alb_field_mom',__FILE__,__LINE__)
    alb_field_mom = 0.0

  !> the following fields are always used in the interface write_netcdf_cosmo_grid_extpar
  !> and must be allocated even IF not used
    IF (raw_id == 2) THEN
      IF (l_use_array_cache) THEN
         CALL allocate_cached('alb_dry', alb_dry, [tg%ie,tg%je,tg%ke])
      ELSE
         allocate(alb_dry(tg%ie,tg%je,tg%ke), stat=errorcode)
      ENDIF
      IF(errorcode.NE.0) CALL logging%error('Cant allocate array alb_dry',__FILE__,__LINE__)
      alb_dry = 0.0

      IF (l_use_array_cache) THEN
         CALL allocate_cached('alb_sat', alb_sat, [tg%ie,tg%je,tg%ke])
      ELSE
         allocate(alb_sat(tg%ie,tg%je,tg%ke), stat=errorcode)
      ENDIF
      IF(errorcode.NE.0) CALL logging%error('Cant allocate array alb_sat',__FILE__,__LINE__)
      alb_sat = 0.0

      IF (l_use_array_cache) THEN
         CALL allocate_cached('alnid_field_mom', alnid_field_mom, [0,0,0,0])
      ELSE
         allocate(alnid_field_mom(0,0,0,0), stat=errorcode)
      ENDIF
        IF(errorcode.NE.0) CALL logging%error('Cant allocate array alnid_field_mom',__FILE__,__LINE__)
      alnid_field_mom = 0.0

      IF (l_use_array_cache) THEN
         CALL allocate_cached('aluvd_field_mom', aluvd_field_mom, [0,0,0,0])
      ELSE
         allocate(aluvd_field_mom(0,0,0,0), stat=errorcode)
      ENDIF
      IF(errorcode.NE.0) CALL logging%error('Cant allocate array aluvd_field_mom',__FILE__,__LINE__)
      aluvd_field_mom = 0.0

      IF (l_use_array_cache) THEN
         CALL allocate_cached('alb_interpol', alb_interpol, [0,0,0,0])
      ELSE
         allocate(alb_interpol(0,0,0,0), stat=errorcode)
      ENDIF
      IF(errorcode.NE.0) CALL logging%error('Cant allocate array alb_interpol',__FILE__,__LINE__)
      alb_interpol = 0.0

    ELSE ! raw_id == 2

      IF (l_use_array_cache) THEN
         CALL allocate_cached('alb_dry', alb_dry, [0,0,0])
      ELSE
         allocate(alb_dry(0,0,0), stat=errorcode)
      ENDIF
      IF(errorcode.NE.0) CALL logging%error('Cant allocate array alb_dry',__FILE__,__LINE__)
      alb_dry = 0.0

      IF (l_use_array_cache) THEN
         CALL allocate_cached('alb_sat', alb_sat, [0,0,0])
      ELSE
         allocate(alb_sat(0,0,0), stat=errorcode)
      ENDIF
      IF(errorcode.NE.0) CALL logging%error('Cant allocate array alb_sat',__FILE__,__LINE__)
      alb_sat = 0.0

      IF (l_use_array_cache) THEN
         CALL allocate_cached('alnid_field_mom', alnid_field_mom, [tg%ie,tg%je,tg%ke,nt])
      ELSE
         allocate(alnid_field_mom(tg%ie,tg%je,tg%ke,nt), stat=errorcode)
      ENDIF
      IF(errorcode.NE.0) CALL logging%error('Cant allocate array alnid_field_mom',__FILE__,__LINE__)
      alnid_field_mom = 0.0

      IF (l_use_array_cache) THEN
         CALL allocate_cached('aluvd_field_mom', aluvd_field_mom, [tg%ie,tg%je,tg%ke,nt])
      ELSE
         allocate(aluvd_field_mom(tg%ie,tg%je,tg%ke,nt), stat=errorcode)
      ENDIF
      IF(errorcode.NE.0) CALL logging%error('Cant allocate array aluvd_field_mom',__FILE__,__LINE__)
      aluvd_field_mom = 0.0

      IF (l_use_array_cache) THEN
         CALL allocate_cached('alb_interpol', alb_interpol, [tg%ie,tg%je,tg%ke,nt])
      ELSE
         allocate(alb_interpol(tg%ie,tg%je,tg%ke,nt), stat=errorcode)
      ENDIF
      IF(errorcode.NE.0) CALL logging%error('Cant allocate array alb_interpol',__FILE__,__LINE__)
      alb_interpol = 0.0
    ENDIF

  END SUBROUTINE allocate_alb_target_fields

  !> allocate fields for ERA target data
  SUBROUTINE allocate_era_target_fields(tg,nt, l_use_array_cache)

    TYPE(target_grid_def), INTENT(IN) :: tg  !< structure with target grid description
    INTEGER (KIND=i4), INTENT(in)     :: nt !< number of timesteps (12 for monthly mean values)
    LOGICAL, INTENT(in)               :: l_use_array_cache

    INTEGER(KIND=i4)                  :: errorcode !< error status variable

    CALL logging%info('Enter routine: allocate_era_target_fields')

    IF (l_use_array_cache) then
      call allocate_cached('sst_field', sst_field, [tg%ie,tg%je,tg%ke,nt])
    ELSE
      allocate(sst_field(tg%ie,tg%je,tg%ke,nt), stat=errorcode)
    ENDIF
      IF(errorcode.NE.0) CALL logging%error('Cant allocate the array sst_field',__FILE__,__LINE__)
      sst_field = 0.0

    IF (l_use_array_cache) then
      call allocate_cached('wsnow_field', wsnow_field, [tg%ie,tg%je,tg%ke,nt])
    ELSE
      allocate(wsnow_field(tg%ie,tg%je,tg%ke,nt), stat=errorcode)
    ENDIF
      IF(errorcode.NE.0) CALL logging%error('Cant allocate the array wsnow_field',__FILE__,__LINE__)
      wsnow_field = 0.0

    IF (l_use_array_cache) then
      call allocate_cached('t2m_field', t2m_field, [tg%ie,tg%je,tg%ke,nt])
    ELSE
      allocate(t2m_field(tg%ie,tg%je,tg%ke,nt), stat=errorcode)
    ENDIF
      IF(errorcode.NE.0) CALL logging%error('Cant allocate the array t2m_field',__FILE__,__LINE__)
      t2m_field = 0.0

    IF (l_use_array_cache) then
      call allocate_cached('hsurf_field', hsurf_field, [tg%ie,tg%je,tg%ke])
    ELSE
      allocate(hsurf_field(tg%ie,tg%je,tg%ke), stat=errorcode)
    ENDIF
      IF(errorcode.NE.0) CALL logging%error('Cant allocate the array hsurf_field',__FILE__,__LINE__)
      hsurf_field = 0.0

    CALL logging%info('Exit routine: allocate_era_target_fields')

  END SUBROUTINE allocate_era_target_fields

  SUBROUTINE allocate_ahf_target_fields(tg, l_use_array_cache)

    TYPE(target_grid_def), INTENT(IN) :: tg  !< structure with target grid description
    LOGICAL, INTENT(in)               :: l_use_array_cache
    INTEGER                           :: errorcode !< error status variable

    IF (l_use_array_cache) then
     call allocate_cached('ahf_field', ahf_field, [tg%ie,tg%je,tg%ke])
    ELSE
     allocate(ahf_field(tg%ie,tg%je,tg%ke), stat=errorcode)
    ENDIF
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array ahf_field',__FILE__,__LINE__)
    ahf_field = 0.0

  END SUBROUTINE allocate_ahf_target_fields

  SUBROUTINE allocate_isa_target_fields(tg, l_use_array_cache)

    TYPE(target_grid_def), INTENT(IN) :: tg  !< structure with target grid description
    LOGICAL, INTENT(in)               :: l_use_array_cache
    
    INTEGER(KIND=i4)                   :: errorcode !< error status variable

    errorcode = 0
    
    CALL logging%info('Enter routine: allocate_isa_target_fields')

    IF (l_use_array_cache) then
     CALL allocate_cached('isa_field', isa_field, [tg%ie,tg%je,tg%ke])
    ELSE
      ALLOCATE(isa_field(tg%ie,tg%je,tg%ke), stat=errorcode)
    ENDIF
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array isa_field',__FILE__,__LINE__)
    isa_field = 0.0

  END SUBROUTINE allocate_isa_target_fields

END MODULE mo_python_tg_fields
