!+ Fortran module with netcdf output routines for external parameters
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_0         2010/12/21 Hermann Asensio
!  Initial release
! V1_1         2011/01/20 Hermann Asensio
!  write optical thickness of aeorosls into seperate variables 
!  in the netcdf output for ICON and COSMO grids        
! V1_2         2011/03/25 Hermann Asensio
!  Adapt netcdf output to standards of COSMO-CLM
! V1_3         2011/04/19 Hermann Asensio
!  Bug fix in netcdf output for COSMO grids
! introduce Globcover 2009 land use data set for external parameters
! V1_4         2011/04/21 Anne Roche
!  implementation of orography smoothing
! V1_7         2013/01/25 Guenther Zaengl 
!   Parallel threads for ICON and COSMO using Open-MP, 
!   Several bug fixes and optimizations for ICON search algorithm, 
!   particularly for the special case of non-contiguous domains; 
!   simplified namelist control for ICON  
! V1_8         2013-03-12 Frank Brenner
!  introduced MODIS albedo dataset(s) as new external parameter(s)         
! V1_11        2013/04/16 Juergen Helmert
!  Small adaptions in output of global attributes
! V2_0         1013-06-04 Anne Roches
!  introduced topographical corrected radiaton parameters as new 
!  external parameters
! V2_0         2013-06-04 Martina Messmer/Daniel Luethi
!  introduced HWSD data set as new external parameters (Juergen Helmert)
!  SSO parameters changed to optional parameters        
!  allow for alternative AOT data sets via raw_data_aot_data
!  allow for prescribed soil albedo instead of albedo       
! V2_0_3       2015-01-12 Juergen Helmert
!  Bugfix correction covers CSCS SVN r5907-r6359
! V3_0         2015-05-21 Juergen Helmert
!  Adaptions for urban fields         
!  correction of double allocation in case of lrad=.FALSE. (B. Rockel)
! V4_0         2016/08/17 Daniel Lthi and authors from RHM
!  added support for subgrid-scale slope parameter
!  added support for MACv2 spectrally stratified monthly aerosol fields
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!> Fortran module with netcdf output routines for external parameters 
!> ouptut routines
!> \author Hermann Asensio
MODULE mo_extpar_output_nc

  USE mo_logging
  USE mo_var_meta_data
  USE mo_kind,                  ONLY: wp, i4
  USE info_extpar,              ONLY: INFO_RevisionHash, INFO_CodeIsModified, &
       &                              INFO_PackageName
                                
  USE mo_grid_structures,       ONLY: rotated_lonlat_grid, &
       &                              icosahedral_triangular_grid, &
       &                              target_grid_def

  USE mo_io_utilities,          ONLY: netcdf_attributes, &
       &                              dim_meta_info, &
       &                              netcdf_put_var, &
       &                              open_new_netcdf_file, &
       &                              close_netcdf_file, &
       &                              netcdf_def_grid_mapping, &
                                      set_date_mm_extpar_field
                                   
  USE mo_io_units,             ONLY: filename_max

  USE mo_albedo_data,          ONLY: ntime_alb, &
       &                             ialb_type, undef_alb_bs

  USE mo_ndvi_data,            ONLY: ntime_ndvi
  USE mo_emiss_data,           ONLY: ntime_emiss

  USE mo_aot_data,             ONLY: ntype_aot, ntime_aot,n_spectr, &
       &                             iaot_type
                               
  USE mo_soil_data,            ONLY: HWSD_data
  USE mo_topo_data,            ONLY: itopo_type, topo_aster, topo_gl

  USE mo_ecoclimap_data,       ONLY: ntime_ecoclimap

  USE mo_topo_data,            ONLY: topo_aster, topo_gl
  
  USE mo_cosmo_grid, ONLY: lon_rot, lat_rot

  USE mo_physical_constants,   ONLY: grav

  USE mo_icon_grid_data,       ONLY: icon_grid_region, &
       &                             clon, clat, &
       &                             clon_vertices, clat_vertices, &
       &                             allocate_icon_coor

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: write_netcdf_cosmo_grid_extpar
  PUBLIC :: write_netcdf_icon_grid_extpar

  CONTAINS


  !> netcdf output of external Parameters for COSMO
  SUBROUTINE write_netcdf_cosmo_grid_extpar(netcdf_filename,     &
       &                                    cosmo_grid,          &
       &                                    tg,                  &
       &                                    isoil_data,          &
       &                                    ldeep_soil,          &
       &                                    itopo_type,          &
       &                                    lsso,                &
       &                                    l_use_isa,           &
       &                                    l_use_ahf,           &
       &                                    l_use_sgsl,          &
       &                                    lscale_separation,   &
       &                                    y_orofilt,           &
       &                                    lrad,                &
       &                                    nhori,               &
       &                                    undefined,           &
       &                                    name_lookup_table_lu,&
       &                                    lu_dataset,          &
       &                                    i_landuse_data,      &
       &                                    nclass_lu,           &
       &                                    lon_geo,             &
       &                                    lat_geo,             &
       &                                    fr_land_lu,          &
       &                                    lu_class_fraction,   &
       &                                    ice_lu,              &
       &                                    z0_lu,               &
       &                                    z0_topo,             &
       &                                    z012_lu,             &
       &                                    root_lu,             &
       &                                    plcov_mn_lu,         &
       &                                    plcov_mx_lu,         &
       &                                    plcov12_lu,          &
       &                                    lai_mn_lu,           & 
       &                                    lai_mx_lu,           &
       &                                    lai12_lu,            &
       &                                    rs_min_lu,           &
       &                                    urban_lu,            &
       &                                    for_d_lu,            &
       &                                    for_e_lu,            &
       &                                    skinc_lu,            &
       &                                    emissivity_lu,       &
       &                                    lake_depth,          &
       &                                    fr_lake,             &
       &                                    soiltype_fao,        &
       &                                    ndvi_max,            &
       &                                    ndvi_field_mom,      &
       &                                    ndvi_ratio_mom,      &
       &                                    emiss_field_mom,     &
       &                                    hh_topo,             &
       &                                    stdh_topo,           &
       &                                    aot_tg,              &
       &                                    MAC_aot_tg,          &
       &                                    MAC_ssa_tg,          &
       &                                    MAC_asy_tg,          & 
       &                                    crutemp,             &
       &                                    alb_field_mom,       &
       &                                    alnid_field_mom,     &
       &                                    aluvd_field_mom,     &
       &                                    alb_dry,             &
       &                                    alb_sat,             &
       &                                    fr_sand,             &
       &                                    fr_silt,             &
       &                                    fr_clay,             &
       &                                    fr_oc,               &
       &                                    fr_bd,               &
       &                                    soiltype_deep,       &
       &                                    fr_sand_deep,        &
       &                                    fr_silt_deep,        &
       &                                    fr_clay_deep,        &
       &                                    fr_oc_deep,          &
       &                                    fr_bd_deep,          &
       &                                    theta_topo,          &
       &                                    aniso_topo,          &
       &                                    slope_topo,          & 
       &                                    slope_asp_topo,      &
       &                                    slope_ang_topo,      &
       &                                    horizon_topo,        &
       &                                    skyview_topo,        &
       &                                    isa_field,           &
       &                                    ahf_field,           &
       &                                    sgsl             )


    CHARACTER (len=*), INTENT(IN)         :: netcdf_filename, & !< filename for the netcdf file
         &                                   y_orofilt, &
         &                                   name_lookup_table_lu, & !< name of lookup table
         &                                   lu_dataset !< name of landuse data set

    TYPE(rotated_lonlat_grid), INTENT(IN) :: cosmo_grid !< structure which contains the definition of the COSMO grid
    TYPE(target_grid_def), INTENT(IN)     :: tg !< structure with target grid description

    INTEGER (KIND=i4),      INTENT(IN)    :: isoil_data, &
         &                                   itopo_type, &
         &                                   nhori, &
         &                                   i_landuse_data, &
         &                                   nclass_lu, &
         &                                   soiltype_fao(:,:,:) !< soiltype due to FAO Digital Soil map of the World

    INTEGER(KIND=i4),INTENT(IN), OPTIONAL :: soiltype_deep(:,:,:) !< soiltype due to FAO Digital Soil map of the World

    LOGICAL,               INTENT(IN)     :: ldeep_soil, &
         &                                   l_use_isa, &
         &                                   l_use_ahf, &
         &                                   l_use_sgsl, &
         &                                   lsso, &
         &                                   lscale_separation, &
         &                                   lrad

    
    REAL (KIND=wp), INTENT(IN)           :: lon_geo(:,:,:), &  !<longitude coordinates of the target grid in the geographical system
         &                                  lat_geo(:,:,:), &  !< latitude coordinates of the target grid in the geographical system
         &                                  lu_class_fraction(:,:,:,:), &  
         &                                  fr_land_lu(:,:,:), & !< fraction land due to lu raw data
         &                                  ice_lu(:,:,:), &     !< fraction of ice due to lu raw data
         &                                  z0_lu(:,:,:), &      !< roughness length
         &                                  z0_glc2000(:,:,:), & !< roughness length due to lu land use data
         &                                  z0_topo(:,:,:), &      !< roughness length due to lu land use data
         &                                  root_lu(:,:,:), &    !< root depth due to lu land use data
         &                                  plcov_mx_lu(:,:,:), &!< plant cover maximum due to lu land use data
         &                                  plcov_mn_lu(:,:,:), &!< plant cover minimum due to lu land use data
         &                                  lai_mx_lu(:,:,:), &  !< Leaf Area Index maximum due to lu land use data
         &                                  lai_mn_lu(:,:,:), &  !< Leaf Area Index minimum due to lu land use data
         &                                  rs_min_lu(:,:,:), &  !< minimal stomata resistance due to lu land use data
         &                                  urban_lu(:,:,:), &   !< urban fraction due to lu land use data
         &                                  for_d_lu(:,:,:), &   !< deciduous forest (fraction) due to lu land use data
         &                                  for_e_lu(:,:,:), &   !< evergreen forest (fraction) due to lu land use data
         &                                  skinc_lu(:,:,:), &   !< skin conductivity due to lu land use data
         &                                  emissivity_lu(:,:,:), & !< longwave emissivity due to lu land use data
         &                                  lake_depth(:,:,:), & !< lake depth
         &                                  fr_lake(:,:,:), &     !< fraction of fresh water (lakes)
         &                                  alb_field_mom(:,:,:,:), & !< field for monthly mean albedo data
         &                                  alnid_field_mom(:,:,:,:), &
         &                                  aluvd_field_mom(:,:,:,:), &
         &                                  ndvi_max(:,:,:), & !< field for ndvi maximum
         &                                  ndvi_field_mom(:,:,:,:), & !< field for monthly mean ndvi data (12 months)
         &                                  ndvi_ratio_mom(:,:,:,:), & !< field for monthly ndvi ratio (12 months)
         &                                  hh_topo(:,:,:), &  !< mean height 
         &                                  stdh_topo(:,:,:), & !< standard deviation of subgrid scale orographic height
         &                                  aot_tg(:,:,:,:,:), & !< aerosol optical thickness, aot_tg(ie,je,ke,ntype,ntime)
         &                                  MAC_aot_tg(:,:,:,:), &
         &                                  MAC_ssa_tg(:,:,:,:), &
         &                                  MAC_asy_tg(:,:,:,:), &
         &                                  crutemp(:,:,:), &  !< cru climatological temperature , crutemp(ie,je,ke) 
         &                                  theta_topo(:,:,:), & !< sso parameter, angle of principal axis
         &                                  aniso_topo(:,:,:), & !< sso parameter, anisotropie factor
         &                                  slope_topo(:,:,:), & !< sso parameter, mean slope
         &                                  plcov12_lu(:,:,:,:), & !<  plcov ecoclimap
         &                                  lai12_lu(:,:,:,:), & !<  lai ecoclimap
         &                                  z012_lu(:,:,:,:), & !<  lai ecoclimap
         &                                  undefined

    REAL (KIND=wp), INTENT(IN), OPTIONAL::  emiss_field_mom(:,:,:,:), & !< field for monthly mean emiss data (12 months)
         &                                  alb_dry(:,:,:), & !< field for soil albedo data
         &                                  alb_sat(:,:,:), &
         &                                  fr_sand(:,:,:), &   !< sand fraction due to HWSD
         &                                  fr_silt(:,:,:), &   !< silt fraction due to HWSD
         &                                  fr_clay(:,:,:), &   !< clay fraction due to HWSD
         &                                  fr_oc(:,:,:), &     !< oc fraction due to HWSD
         &                                  fr_bd(:,:,:), &     !< bulk density due to HWSD
         &                                  fr_sand_deep(:,:,:), &   !< sand fraction due to HWSD
         &                                  fr_silt_deep(:,:,:), &   !< silt fraction due to HWSD
         &                                  fr_clay_deep(:,:,:), &   !< clay fraction due to HWSD
         &                                  fr_oc_deep(:,:,:), &     !< oc fraction due to HWSD
         &                                  fr_bd_deep(:,:,:), &     !< bulk density due to HWSD
         &                                  slope_asp_topo(:,:,:), &   !< lradtopo parameter, slope_aspect
         &                                  slope_ang_topo(:,:,:), &   !< lradtopo parameter, slope_angle
         &                                  horizon_topo  (:,:,:,:), & !< lradtopo parameter, horizon
         &                                  skyview_topo  (:,:,:), &   !< lradtopo parameter, skyview
         &                                  isa_field(:,:,:), & !< impervious surface area
         &                                  ahf_field(:,:,:), & !< field for ahf
         &                                  sgsl(:,:,:) !< field for subgrid-scale slope

    ! local variables
    REAL(KIND=wp), ALLOCATABLE          :: z012tot(:,:,:), & !<  z0 ecoclimap plant+oro
         &                                 var_real_2d(:,:), &
         &                                 var_real_hor(:,:,:), &
         &                                 var_real_MAC(:,:,:,:), &
         &                                 time(:)   
                                        
    INTEGER (KIND=i4)                   :: dataDate, &  !< date, for edition independent use of GRIB_API dataDate as Integer in the format ccyymmdd
         &                                 dataTime, &  !< time, for edition independent use GRIB_API dataTime in the format hhmm
         &                                 ndims, &  
         &                                 ncid, &
         &                                 varid, &
         &                                 errorcode, n !< error status variable

    INTEGER(KIND=i4), PARAMETER         :: nglob_atts=11
                                        
    LOGICAL                             :: l_use_emiss=.FALSE. !< flag if additional CAMEL emissivity data are present
    TYPE(dim_meta_info), ALLOCATABLE    :: dim_list(:) !< dimensions for netcdf file
    TYPE(netcdf_attributes)             :: global_attributes(nglob_atts)

    CHARACTER (len=filename_max)        :: namelist_file !< filename with namelists for for EXTPAR settings for optional output

    CHARACTER (len=80)::                   grid_mapping, & !< netcdf attribute grid mapping
         &                                 coordinates  !< netcdf attribute coordinates

    CALL logging%info('Enter routine: write_netcdf_cosmo_grid_extpar')

    !-------------------------------------------------------------
    ! define global attributes
    CALL set_global_att_extpar(global_attributes,name_lookup_table_lu,lu_dataset,isoil_data,lscale_separation,y_orofilt)
    WRITE(message_text,*) '----------------   NetCDF global_attributes ----------------------'
    CALL logging%info(message_text)
    DO n=1,nglob_atts
      WRITE(message_text,*) global_attributes(n)
      CALL logging%info(message_text)
    END DO
    WRITE(message_text,*) '------------------------------------------------------------------'
    CALL logging%info(message_text)

    !set up dimensions for buffer
    CALL  def_dimension_info_buffer(tg,nhori=nhori)
    ! dim_3d_tg

    !set up dimensions for COSMO grid
    CALL def_dimension_info_cosmo(cosmo_grid,nhori=nhori)
    ! dim_rlon_cosmo, dim_rlat_cosmo, dim_2d_cosmo, rlon_meta, rlat_meta

    ! set mapping parameters for netcdf
    grid_mapping="rotated_pole"
    coordinates="lon lat"
    CALL set_nc_grid_def_cosmo(cosmo_grid,grid_mapping)
    ! nc_grid_def_cosmo
    ! define meta information for target field variables lon_geo, lat_geo 
    CALL def_com_target_fields_meta(dim_2d_cosmo,coordinates,grid_mapping)
    ! lon_geo_meta and lat_geo_meta

    ! define meta information for various land use related variables (GLC2000) for netcdf output
    CALL def_isa_fields_meta(dim_2d_cosmo,coordinates,grid_mapping)

    ! define meta information for various land use related variables for netcdf output
    IF (i_landuse_data .EQ. 4) THEN
      CALL  def_ecoclimap_fields_meta(ntime_ecoclimap,nclass_lu,dim_2d_cosmo,coordinates,grid_mapping) 
    ELSE
      CALL def_lu_fields_meta(nclass_lu,dim_2d_cosmo,lu_dataset,coordinates,grid_mapping)
    ENDIF

    CALL def_soil_meta(dim_2d_cosmo,isoil_data,coordinates,grid_mapping)
    !  fr_land_soil_meta, soiltype_fao_meta

    CALL def_alb_meta(ntime_alb,dim_2d_cosmo,coordinates,grid_mapping)

    !define meta information for various AHF data related variables for netcdf output
    CALL def_ahf_meta(dim_2d_cosmo,coordinates,grid_mapping)
    ! dim_ahf_tg, ahf_field_meta

    !define meta information for various NDVI data related variables for netcdf output
    CALL def_ndvi_meta(ntime_ndvi,dim_2d_cosmo,coordinates,grid_mapping)
    ! dim_ndvi_tg, ndvi_max_meta, ndvi_field_mom_meta, ndvi_ratio_mom_meta

    !define meta information for various EMISS data related variables for netcdf output
    CALL def_emiss_meta(ntime_emiss,dim_2d_cosmo,coordinates,grid_mapping)
    ! dim_emiss_tg, emiss_max_meta, emiss_field_mom_meta, emiss_ratio_mom_meta

    ! define meta information for various TOPO data related variables for netcdf output
    IF(lrad) THEN
      CALL def_topo_meta(dim_2d_cosmo,itopo_type,coordinates=coordinates,grid_mapping=grid_mapping,diminfohor=dim_3d_cosmo)
    ELSE
      CALL def_topo_meta(dim_2d_cosmo,itopo_type,coordinates=coordinates,grid_mapping=grid_mapping)
    ENDIF

    !define meta information for subgrid-scale slope data related variables for netcdf output
    IF (l_use_sgsl) THEN
      CALL def_sgsl_meta(dim_2d_cosmo,itopo_type,coordinates,grid_mapping)
    ENDIF

    ! define dimensions and meta information for variable aot_tg for netcdf output
    CALL def_aot_tg_meta(ntime_aot,ntype_aot,dim_2d_cosmo,coordinates,grid_mapping)
    ! dim_aot_tg and aot_tg_meta
    ! dim_aot_ty, aer_bc_meta, aer_dust_meta, aer_org_meta, aer_so4_meta, aer_ss_meta

    ! define meta information for variable crutemp for netcdf output
    CALL def_crutemp_meta(dim_2d_cosmo,coordinates,grid_mapping)
    ! crutemp_meta

    ! define meta information for various land use related variables (FLAKE) for netcdf output

    CALL def_flake_fields_meta(dim_2d_cosmo,coordinates,grid_mapping)
    ! lake_depth_meta, fr_lake_meta, &
    !  &       flake_tot_npixel_meta

    ALLOCATE(var_real_2d(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot), STAT=errorcode)
    IF (errorcode /= 0 ) CALL logging%error('Cant allocate var_real_2d',__FILE__,__LINE__)

    ! z0 tot veg + topo
    IF (i_landuse_data .EQ. 4) THEN
      ALLOCATE(z012tot(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1:12), STAT=errorcode)
      IF (errorcode /= 0 ) CALL logging%error('Cant allocate z012tot' ,__FILE__,__LINE__ )
      z012tot(:,:,:) =  z012_lu(:,:,1,:) 
      DO n=1, ntime_ecoclimap   
        z012tot(:,:,n) =  z012tot(:,:,n) + z0_topo(:,:,1)
      ENDDO
    ENDIF
    !>
    !set up dimensions for buffer netcdf output 

    ndims = 4
    IF(lrad) ndims = ndims + 1

    ALLOCATE(time(1:ntime_aot),STAT=errorcode)
    IF (errorcode /= 0 ) CALL logging%error('Cant allocate array time',__FILE__,__LINE__)
    DO n=1,ntime_aot
      CALL set_date_mm_extpar_field(n,dataDate,dataTime)
      time(n) = REAL(dataDate,wp) + REAL(dataTime,wp)/10000. ! units = "day as %Y%m%d.%f"
    ENDDO
    ALLOCATE(dim_list(1:ndims),STAT=errorcode)
    IF (errorcode /= 0 ) CALL logging%error('Cant allocate array dim_list',__FILE__,__LINE__)

    dim_list(1) = dim_rlat_cosmo(1) ! rlat
    dim_list(2) = dim_rlon_cosmo(1) ! rlon
    dim_list(3)%dimname = 'nclass_lu'  ! nclass_lu
    dim_list(3)%dimsize = nclass_lu 
    dim_list(4)%dimname = 'time'
    dim_list(4)%dimsize = ntime_aot

    IF (lrad) THEN
      dim_list(5)%dimname = 'nhori'
      dim_list(5)%dimsize = nhori
    ENDIF

    !-----------------------------------------------------------------
    CALL open_new_netcdf_file(netcdf_filename=TRIM(netcdf_filename),   &
         &                       dim_list=dim_list,                  &
         &                       global_attributes=global_attributes, &
         &                       time=time,                           &
         &                       ncid=ncid)
    !-----------------------------------------------------------------
    ! start with real 1d variables

    ! rlon
    CALL netcdf_put_var(ncid,lon_rot(1:cosmo_grid%nlon_rot),rlon_meta,undefined)

    CALL netcdf_put_var(ncid,lat_rot(1:cosmo_grid%nlat_rot),rlat_meta,undefined)

    ! fr_land_lu
    var_real_2d(:,:) = fr_land_lu(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1)
    CALL netcdf_put_var(ncid,var_real_2d,fr_land_lu_meta,undefined)

    ! ice_lu
    var_real_2d(:,:) = ice_lu(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1)
    CALL netcdf_put_var(ncid,var_real_2d,ice_lu_meta,undefined)

    IF (i_landuse_data .EQ. 4) THEN
      CALL netcdf_put_var(ncid,                                  &
           & plcov12_lu(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1,1:ntime_ecoclimap), &
           & plcov12_lu_meta,undefined)

      CALL netcdf_put_var(ncid,                                  &
           & lai12_lu(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1,1:ntime_ecoclimap), &
           & lai12_lu_meta,undefined)
      CALL netcdf_put_var(ncid,                                  &
           & z012_lu(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1,1:ntime_ecoclimap), &
           & z012_lu_meta,undefined)
      CALL netcdf_put_var(ncid,                                  &
           & z012tot(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1:ntime_ecoclimap), &
           & z012_tot_meta,undefined)

      ! plcov_mx_lu
      var_real_2d(:,:) = SUM(MAXVAL(plcov12_lu,DIM=4),DIM=3)
      CALL netcdf_put_var(ncid,var_real_2d,plcov_mx_lu_meta,undefined)

      ! plcov_mn_lu
      var_real_2d(:,:) = SUM(MINVAL(plcov12_lu,DIM=4),DIM=3)
      CALL netcdf_put_var(ncid,var_real_2d,plcov_mn_lu_meta,undefined)

      ! lai_mx_lu
      var_real_2d(:,:) = SUM(MAXVAL(lai12_lu,DIM=4),DIM=3)
      CALL netcdf_put_var(ncid,var_real_2d,lai_mx_lu_meta,undefined)

      ! lai_mn_lu
      var_real_2d(:,:) = SUM(MINVAL(lai12_lu,DIM=4),DIM=3)
      CALL netcdf_put_var(ncid,var_real_2d,lai_mn_lu_meta,undefined)

      ! z0_lu
      var_real_2d(:,:) = SUM(z012tot,DIM=3)/REAL(ntime_ecoclimap,wp)
      CALL netcdf_put_var(ncid,var_real_2d,z0_lu_meta,undefined)
    ELSE 
      ! plcov_mn_lu
      var_real_2d(:,:) = plcov_mn_lu(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1)
      CALL netcdf_put_var(ncid,var_real_2d,plcov_mn_lu_meta,undefined)

      ! plcov_mx_lu
      var_real_2d(:,:) = plcov_mx_lu(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1)
      CALL netcdf_put_var(ncid,var_real_2d,plcov_mx_lu_meta,undefined)

      ! lai_mn_lu
      var_real_2d(:,:) = lai_mn_lu(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1)
      CALL netcdf_put_var(ncid,var_real_2d,lai_mn_lu_meta,undefined)

      ! lai_mx_lu
      var_real_2d(:,:) = lai_mx_lu(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1)
      CALL netcdf_put_var(ncid,var_real_2d,lai_mx_lu_meta,undefined)

      ! z0_lu
      var_real_2d(:,:) = z0_lu(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1)
      CALL netcdf_put_var(ncid,var_real_2d,z0_lu_meta,undefined)
    ENDIF

    ! emissivity_lu
    var_real_2d(:,:) = emissivity_lu(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1)
    CALL netcdf_put_var(ncid,var_real_2d,emissivity_lu_meta,undefined)

    ! rs_min_lu
    var_real_2d(:,:) = rs_min_lu(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1)
    CALL netcdf_put_var(ncid,var_real_2d,rs_min_lu_meta,undefined)

    ! urban_lu
    var_real_2d(:,:) = urban_lu(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1)
    CALL netcdf_put_var(ncid,var_real_2d,urban_lu_meta,undefined)

    IF (l_use_isa) THEN
      ! isa_field
      var_real_2d(:,:) = isa_field(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1)
      CALL netcdf_put_var(ncid,var_real_2d,isa_field_meta,undefined)
    END IF

    IF (l_use_ahf) THEN
      ! hw-marker. maybe provide a switch to either include/exclude ahf
      ! ahf_field
      var_real_2d(:,:) = ahf_field(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1)
      CALL netcdf_put_var(ncid,var_real_2d,ahf_field_meta,undefined)
    END IF

    IF (l_use_sgsl) THEN
      ! hw-marker. maybe provide a switch to either include/exclude ahf
      ! ahf_field
      var_real_2d(:,:) = sgsl(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1)
      CALL netcdf_put_var(ncid,var_real_2d,sgsl_meta,undefined)
    END IF

    ! for_d_lu
    var_real_2d(:,:) = for_d_lu(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1)
    CALL netcdf_put_var(ncid,var_real_2d,for_d_lu_meta,undefined)

    ! for_e_lu
    var_real_2d(:,:) = for_e_lu(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1)
    CALL netcdf_put_var(ncid,var_real_2d,for_e_lu_meta,undefined)

    ! skinc_lu
    var_real_2d(:,:) = skinc_lu(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1)
    CALL netcdf_put_var(ncid,var_real_2d,skinc_lu_meta,undefined)

    ! root_lu
    var_real_2d(:,:) = root_lu(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1)
    CALL netcdf_put_var(ncid,var_real_2d,root_lu_meta,undefined)

    ! lon
    var_real_2d(:,:) = lon_geo(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1)
    CALL netcdf_put_var(ncid,var_real_2d,lon_geo_meta,undefined)

    ! lat
    var_real_2d(:,:) = lat_geo(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1)
    CALL netcdf_put_var(ncid,var_real_2d,lat_geo_meta,undefined)

    ! ndvi_max
    var_real_2d(:,:) = ndvi_max(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1)
    CALL netcdf_put_var(ncid,var_real_2d,ndvi_max_meta,undefined)

    ! hh_topo
    var_real_2d(:,:) = hh_topo(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1)
    CALL netcdf_put_var(ncid,var_real_2d,hh_topo_meta,undefined)
    !roa nc>
    ! hh_fis
    var_real_2d(:,:) = grav * hh_topo(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1)
    CALL netcdf_put_var(ncid,var_real_2d,hh_fis_meta,undefined)
    !roa nc<

    ! stdh_topo
    var_real_2d(:,:) = stdh_topo(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1)
    CALL netcdf_put_var(ncid,var_real_2d,stdh_topo_meta,undefined)

    ! theta_topo
    IF (lsso) THEN
      var_real_2d(:,:) = theta_topo(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1)
      CALL netcdf_put_var(ncid,var_real_2d,theta_topo_meta,undefined)
    ENDIF

    ! aniso_topo
    IF (lsso) THEN
      var_real_2d(:,:) = aniso_topo(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1)
      CALL netcdf_put_var(ncid,var_real_2d,aniso_topo_meta,undefined)
    ENDIF

    ! slope_topo
    IF (lsso) THEN
      var_real_2d(:,:) = slope_topo(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1)
      CALL netcdf_put_var(ncid,var_real_2d,slope_topo_meta,undefined)
    ENDIF

    ! slope_asp_topo
    IF (lrad) THEN
      var_real_2d(:,:) = slope_asp_topo(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1)
      CALL netcdf_put_var(ncid,var_real_2d, slope_asp_topo_meta,undefined)
    ENDIF

    ! slope_ang_topo
    IF (lrad) THEN
      var_real_2d(:,:) = slope_ang_topo(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1)
      CALL netcdf_put_var(ncid,var_real_2d, slope_ang_topo_meta,undefined)
    ENDIF

    ! horizon_topo
    IF (lrad) THEN
      ALLOCATE(var_real_hor(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1:nhori), STAT=errorcode)
      IF (errorcode /= 0 ) CALL logging%error('Cant allocate var_real_hor',__FILE__,__LINE__)
      var_real_hor(:,:,:) = horizon_topo(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1,1:nhori)
      CALL netcdf_put_var(ncid,var_real_hor, horizon_topo_meta,undefined)
      DEALLOCATE(var_real_hor)
    ENDIF

    ! skyview
    IF (lrad) THEN
      var_real_2d(:,:) = skyview_topo(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1)
      CALL netcdf_put_var(ncid,var_real_2d, skyview_topo_meta,undefined)
    ENDIF

    ! crutemp
    var_real_2d(:,:) = crutemp(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1)
    CALL netcdf_put_var(ncid,var_real_2d,crutemp_meta,undefined)

    ! fr_lake
    var_real_2d(:,:) = fr_lake(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1)
    CALL netcdf_put_var(ncid,var_real_2d,fr_lake_meta,undefined)

    ! lake_depth
    var_real_2d(:,:) = lake_depth(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1)
    CALL netcdf_put_var(ncid,var_real_2d,lake_depth_meta,undefined)

    !-----------------------------------------------------------------
    ! soiltype
    ! change HA: write soiltype as real variable as requested by CLM community
    !CALL netcdf_put_var(ncid,soiltype_fao(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1), &
    ! &                 soiltype_fao_meta,undef_int)

    var_real_2d(:,:) = soiltype_fao(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1)
    CALL netcdf_put_var(ncid,var_real_2d,soiltype_fao_meta,undefined)

    ! soiltype_deep
    IF (ldeep_soil) THEN
      var_real_2d(:,:) = soiltype_deep(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1)
      CALL netcdf_put_var(ncid,var_real_2d,soiltype_FAO_deep_meta,undefined)
    ENDIF

    ! fr_sand
    IF (isoil_data == HWSD_data) THEN
      var_real_2d(:,:) = fr_sand(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1)
      CALL netcdf_put_var(ncid,var_real_2d,HWSD_SAND_meta,undefined)
    ENDIF

    ! fr_silt
    IF (isoil_data == HWSD_data) THEN
      var_real_2d(:,:) = fr_silt(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1)
      CALL netcdf_put_var(ncid,var_real_2d,HWSD_SILT_meta,undefined)
    ENDIF

    ! fr_clay
    IF (isoil_data == HWSD_data) THEN
      var_real_2d(:,:) = fr_clay(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1)
      CALL netcdf_put_var(ncid,var_real_2d,HWSD_CLAY_meta,undefined)
    ENDIF

    ! fr_oc
    IF (isoil_data == HWSD_data) THEN
      var_real_2d(:,:) = fr_oc(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1)
      CALL netcdf_put_var(ncid,var_real_2d,HWSD_OC_meta,undefined)
    ENDIF

    ! fr_bd
    IF (isoil_data == HWSD_data) THEN
      var_real_2d(:,:) = fr_bd(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1)
      CALL netcdf_put_var(ncid,var_real_2d,HWSD_BD_meta,undefined)
    ENDIF

    ! fr_sand_deep
    IF (ldeep_soil) THEN
      var_real_2d(:,:) = fr_sand_deep(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1)
      CALL netcdf_put_var(ncid,var_real_2d,HWSD_SAND_deep_meta,undefined)
    ENDIF

    ! fr_silt_deep
    IF (ldeep_soil) THEN
      var_real_2d(:,:) = fr_silt_deep(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1)
      CALL netcdf_put_var(ncid,var_real_2d,HWSD_SILT_deep_meta,undefined)
    ENDIF

    ! fr_clay_deep
    IF (ldeep_soil) THEN
      var_real_2d(:,:) = fr_clay_deep(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1)
      CALL netcdf_put_var(ncid,var_real_2d,HWSD_CLAY_deep_meta,undefined)
    ENDIF

    ! fr_oc_deep
    IF (ldeep_soil) THEN
      var_real_2d(:,:) = fr_oc_deep(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1)
      CALL netcdf_put_var(ncid,var_real_2d,HWSD_OC_deep_meta,undefined)
    ENDIF

    ! fr_bd_deep
    IF (ldeep_soil) THEN
      var_real_2d(:,:) = fr_bd_deep(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1)
      CALL netcdf_put_var(ncid,var_real_2d,HWSD_BD_deep_meta,undefined)
    ENDIF

    !-----------------------------------------------------------------
    ! lu_class_fraction
    CALL netcdf_put_var(ncid,&
         & lu_class_fraction(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1,1:nclass_lu), &
         & lu_class_fraction_meta, &
         & undefined)


    IF (ialb_type == 2) THEN
      !-----------------------------------------------------------------
      ! alb_dry
      CALL netcdf_put_var(ncid,&
           & alb_dry(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1), &
           & alb_dry_meta, &
           & undef_alb_bs)

      ! alb_sat
      CALL netcdf_put_var(ncid,&
           & alb_sat(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1), &
           & alb_sat_meta, &
           & undef_alb_bs)

    ELSE IF (ialb_type == 1) THEN
      !-----------------------------------------------------------------
      ! alb_field_mom
      CALL alb_field_mom_meta%overwrite_units('1')
      CALL netcdf_put_var(ncid,&
           & alb_field_mom(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1,1:ntime_alb)/100., &
           & alb_field_mom_meta, &
           & undefined)

      !-----------------------------------------------------------------
      ! alnid_field_mom
      CALL alnid_field_mom_meta%overwrite_units('1')
      CALL netcdf_put_var(ncid,&
           & alnid_field_mom(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1,1:ntime_alb)/100., &
           & alnid_field_mom_meta, &
           & undefined)

      !-----------------------------------------------------------------
      ! aluvd_field_mom
      CALL aluvd_field_mom_meta%overwrite_units('1')
      CALL netcdf_put_var(ncid,&
           & aluvd_field_mom(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1,1:ntime_alb)/100., &
           & aluvd_field_mom_meta, &
           & undefined)

    ELSE IF (ialb_type == 3) THEN
      !-----------------------------------------------------------------
      ! alb_field_mom
      CALL alb_field_mom_meta%overwrite_units('1')
      CALL netcdf_put_var(ncid,&
           & alb_field_mom(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1,1:ntime_alb)/100., &
           & alb_field_mom_meta, &
           & undefined)
    ENDIF

    !-----------------------------------------------------------------
    ! ndvi_field_mom
    CALL netcdf_put_var(ncid,&
         & ndvi_field_mom(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1,1:ntime_ndvi), &
         & ndvi_field_mom_meta, &
         & undefined)

    ! ndvi_ratio_mom
    CALL netcdf_put_var(ncid,&
         & ndvi_ratio_mom(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1,1:ntime_ndvi), &
         & ndvi_ratio_mom_meta, &
         & undefined)

  namelist_file="INPUT_EMISS"
  INQUIRE(FILE=TRIM(namelist_file), EXIST=l_use_emiss)
  if (l_use_emiss) then
  ! emiss_field_mom
    CALL netcdf_put_var(ncid,&
         & emiss_field_mom(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1,1:ntime_emiss), &
         & emiss_field_mom_meta, &
         & undefined)
  end if
    !-----------------------------------------------------------------
    ! aot
    IF (iaot_type == 4) THEN
      ALLOCATE(var_real_MAC(cosmo_grid%nlon_rot, cosmo_grid%nlat_rot,n_spectr,ntime_aot))

      var_real_MAC(:,:,:,:)=MAC_aot_tg(:,:,:,:)
      CALL netcdf_put_var(ncid,var_real_MAC,aot_tg_MAC_meta,undefined)

      var_real_MAC(:,:,:,:)=MAC_ssa_tg(:,:,:,:)
      CALL netcdf_put_var(ncid,var_real_MAC,ssa_tg_MAC_meta,undefined)

      var_real_MAC(:,:,:,:)=MAC_asy_tg(:,:,:,:)
      CALL netcdf_put_var(ncid,var_real_MAC,asy_tg_MAC_meta,undefined)

    ELSE
      ALLOCATE(var_real_hor(cosmo_grid%nlon_rot,cosmo_grid%nlat_rot,ntime_aot))
      n=1 ! aot_bc
      var_real_hor(:,:,:)=aot_tg(:,:,1,1,:)
      !     CALL netcdf_put_var(ncid,aot_tg(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1,1,1:ntime_aot), &
      !       &                 aer_bc_meta, undefined)
      CALL netcdf_put_var(ncid,var_real_hor,aer_bc_meta, undefined)

      n=2 ! aot_dust
      var_real_hor(:,:,:)=aot_tg(:,:,1,2,:)
      !     CALL netcdf_put_var(ncid,aot_tg(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1,2,1:ntime_aot), &
      !       &                 aer_dust_meta, undefined)
      CALL netcdf_put_var(ncid,var_real_hor,aer_dust_meta, undefined)

      n=3 ! aot_org
      var_real_hor(:,:,:)=aot_tg(:,:,1,3,:)
      !     CALL netcdf_put_var(ncid,aot_tg(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1,3,1:ntime_aot), &
      !       &                 aer_org_meta, undefined)
      CALL netcdf_put_var(ncid,var_real_hor,aer_org_meta, undefined)

      n=4 ! aot_so4
      var_real_hor(:,:,:)=aot_tg(:,:,1,4,:)
      !     CALL netcdf_put_var(ncid,aot_tg(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1,4,1:ntime_aot), &
      !       &                 aer_so4_meta, undefined)
      CALL netcdf_put_var(ncid,var_real_hor,aer_so4_meta, undefined)

      n=5 ! aot_ss
      var_real_hor(:,:,:)=aot_tg(:,:,1,5,:)
      !     CALL netcdf_put_var(ncid,aot_tg(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1,5,1:ntime_aot), &
      !       &                 aer_ss_meta, undefined)
      CALL netcdf_put_var(ncid,var_real_hor,aer_ss_meta, undefined)
    ENDIF

    !-----------------------------------------------------------------
    CALL netcdf_def_grid_mapping(ncid, nc_grid_def_cosmo, varid)

    CALL close_netcdf_file(ncid)

    CALL logging%info('Exit routine: write_netcdf_cosmo_grid_extpar')

  END SUBROUTINE write_netcdf_cosmo_grid_extpar
  !-----------------------------------------------------------------------

  !> netcdf output of GLC2000 derived ICON fields
  SUBROUTINE write_netcdf_icon_grid_extpar(netcdf_filename,     &
       &                                   icon_grid,           &
       &                                   tg,                  &
       &                                   isoil_data,          &
       &                                   ldeep_soil,          &
       &                                   itopo_type,          &
       &                                   lsso,                &
       &                                   l_use_isa,           &
       &                                   l_use_ahf,           &
       &                                   undefined,           &
       &                                   undef_int,           &
       &                                   name_lookup_table_lu,&
       &                                   lu_dataset,          &
       &                                   nclass_lu,           &
       &                                   lon_geo,             &
       &                                   lat_geo,             &
       &                                   fr_land_lu,          &
       &                                   lu_class_fraction,   &
       &                                   ice_lu,              &
       &                                   z0_lu,               &
       &                                   root_lu,             &  
       &                                   plcov_mx_lu,         &
       &                                   lai_mx_lu,           &
       &                                   rs_min_lu,           &
       &                                   urban_lu,            &
       &                                   for_d_lu,            &
       &                                   for_e_lu,            &
       &                                   skinc_lu,            &
       &                                   emissivity_lu,       &
       &                                   lake_depth,          &
       &                                   fr_lake,             &
       &                                   soiltype_fao,        &
       &                                   ndvi_max,            &
       &                                   ndvi_field_mom,      &
       &                                   ndvi_ratio_mom,      &
       &                                   emiss_field_mom,      &
       &                                   hh_topo,             &
       &                                   hh_topo_max,         &
       &                                   hh_topo_min,         &       
       &                                   stdh_topo,           &
       &                                   theta_topo,          &
       &                                   aniso_topo,          & 
       &                                   slope_topo,          &
       &                                   aot_tg,              &
       &                                   crutemp,             &
       &                                   alb_field_mom,       &
       &                                   alnid_field_mom,     &
       &                                   aluvd_field_mom,     &      
       &                                   fr_sand,             &
       &                                   fr_silt,             &
       &                                   fr_clay,             &
       &                                   fr_oc,               &
       &                                   fr_bd,               &
       &                                   soiltype_deep,       &
       &                                   fr_sand_deep,        &
       &                                   fr_silt_deep,        &
       &                                   fr_clay_deep,        &
       &                                   fr_oc_deep,          &
       &                                   fr_bd_deep,          &
       &                                   isa_field,           &
       &                                   ahf_field,           &
       &                                   sst_field,           &
       &                                   wsnow_field,         &          
       &                                   t2m_field,           &          
       &                                   hsurf_field          )    



    CHARACTER (len=*), INTENT(IN)                 :: netcdf_filename, &      !< filename for the netcdf file
         &                                           name_lookup_table_lu, & !< name of lookup table
         &                                           lu_dataset !< name of landuse data set

    TYPE(icosahedral_triangular_grid), INTENT(IN) :: icon_grid !< definition of the ICON grid
    TYPE(target_grid_def), INTENT(IN)             :: tg                    !< target grid description
    LOGICAL,               INTENT(IN)             :: ldeep_soil, &
         &                                           l_use_isa, &
         &                                           l_use_ahf, &
         &                                           lsso

    INTEGER, INTENT(IN)                           :: undef_int, &       !< value to indicate undefined grid elements
         &                                           isoil_data, &
         &                                           itopo_type, &
         &                                           soiltype_fao(:,:,:), & !< soiltype due to FAO Digital Soil map of the World
         &                                           nclass_lu !< number of classes for the land use description

    INTEGER(KIND=i4), INTENT(IN), OPTIONAL       :: soiltype_deep(:,:,:) !< soiltype due to FAO Digital Soil map of the World

    REAL(KIND=wp), INTENT(IN)                    :: undefined, &       !< value to indicate undefined grid elements 
         &                                          lon_geo(:,:,:), &  !< longitude coordinates of the target grid in the geographical system
         &                                          lat_geo(:,:,:), &  !< latitude coordinates of the target grid in the geographical system
         &                                          lu_class_fraction(:,:,:,:), &  
         &                                          fr_land_lu(:,:,:), & !< fraction land due to lu raw data
         &                                          ice_lu(:,:,:), &     !< fraction of ice due to lu raw data
         &                                          z0_lu(:,:,:), &      !< roughness length due to lu land use data
         &                                          root_lu(:,:,:), &    !< root depth due to lu land use data
         &                                          plcov_mx_lu(:,:,:), &!< plant cover maximum due to lu land use data
         &                                          lai_mx_lu(:,:,:), &  !< Leaf Area Index maximum due to lu land use data
         &                                          rs_min_lu(:,:,:), &  !< minimal stomata resistance due to lu land use data
         &                                          urban_lu(:,:,:), &   !< urban fraction due to lu land use data
         &                                          for_d_lu(:,:,:), &   !< deciduous forest (fraction) due to lu land use data
         &                                          for_e_lu(:,:,:), &   !< evergreen forest (fraction) due to lu land use data
         &                                          emissivity_lu(:,:,:), & !< longwave emissivity due to lu land use data
         &                                          skinc_lu(:,:,:), &   !< skin conductivity due to lu land use data
         &                                          lake_depth(:,:,:), & !< lake depth
         &                                          fr_lake(:,:,:), &     !< fraction of fresh water (lakes)
         &                                          alb_field_mom(:,:,:,:), &!< field for monthly mean albedo data
         &                                          alnid_field_mom(:,:,:,:), &
         &                                          aluvd_field_mom(:,:,:,:), &
         &                                          ndvi_max(:,:,:), & !< field for ndvi maximum
         &                                          ndvi_field_mom(:,:,:,:), & !< field for monthly mean ndvi data (12 months)
         &                                          ndvi_ratio_mom(:,:,:,:), & !< field for monthly ndvi ratio (12 months)
         &                                          sst_field(:,:,:,:), & !< field for monthly mean sst data (12 months)
         &                                          wsnow_field(:,:,:,:), & !< field for monthly mean wsnow data (12 months)
         &                                          t2m_field(:,:,:,:), & !< field for monthly mean wsnow data (12 months)
         &                                          hsurf_field(:,:,:), & !< field for monthly mean wsnow data (12 months)
         &                                          hh_topo(:,:,:), &  !< mean height
         &                                          hh_topo_max(:,:,:), &  !< max height on a gridpoint
         &                                          hh_topo_min(:,:,:), &  !< min height on a gridpoint
         &                                          stdh_topo(:,:,:), & !< standard deviation of subgrid scale orographic height
         &                                          aot_tg(:,:,:,:,:), & !< aerosol optical thickness, aot_tg(ie,je,ke,ntype,ntime)
         &                                          crutemp(:,:,:)  !< cru climatological temperature , crutemp(ie,je,ke) 

    REAL (KIND=wp), INTENT(IN), OPTIONAL         :: emiss_field_mom(:,:,:,:), & !< field for monthly mean emiss data (12 months)
         &                                          fr_sand(:,:,:), &   !< sand fraction due to HWSD
         &                                          fr_silt(:,:,:), &   !< silt fraction due to HWSD
         &                                          fr_clay(:,:,:), &   !< clay fraction due to HWSD
         &                                          fr_oc(:,:,:), &     !< oc fraction due to HWSD
         &                                          fr_bd(:,:,:), &     !< bulk density due to HWSD
         &                                          fr_sand_deep(:,:,:), &   !< sand fraction due to HWSD
         &                                          fr_silt_deep(:,:,:), &   !< silt fraction due to HWSD
         &                                          fr_clay_deep(:,:,:), &   !< clay fraction due to HWSD
         &                                          fr_oc_deep(:,:,:), &     !< oc fraction due to HWSD
         &                                          fr_bd_deep(:,:,:), &     !< bulk density due to HWSD
         &                                          theta_topo(:,:,:), & !< sso parameter, angle of principal axis
         &                                          aniso_topo(:,:,:), & !< sso parameter, anisotropie factor
         &                                          slope_topo(:,:,:), & !< sso parameter, mean slope
         &                                          isa_field(:,:,:), & !< field for isa 
         &                                          ahf_field(:,:,:) !< field for ahf 

    ! local variables
    INTEGER(KIND=i4)                             :: ndims, ncid, &
         &                                          dataDate, &  !< date, for edition independent use of GRIB_API dataDate as Integer in the format ccyymmdd
         &                                          dataTime, &  !< time, for edition independent use GRIB_API dataTime in the format hhmm
         &                                          errorcode, &
         &                                          vert_id, &
         &                                          n,nc,nv

    TYPE(dim_meta_info), ALLOCATABLE             :: dim_list(:) !< dimensions for netcdf file
    TYPE(dim_meta_info), TARGET                  :: dim_1d_icon(1:1)
    REAL (KIND=wp), ALLOCATABLE                  :: time(:) !< time variable
                                                 
    LOGICAL                                      :: l_use_emiss=.FALSE. !< flag if additional CAMEL emissivity data are present
    INTEGER(KIND=i4), PARAMETER                  :: nglob_atts=9
    TYPE(netcdf_attributes)                      :: global_attributes(nglob_atts)
                                                 
    CHARACTER (len=80)                           :: grid_mapping !< netcdf attribute grid mapping
    CHARACTER (len=filename_max)                 :: namelist_file !< filename with namelists for for EXTPAR settings for optional output

    !-------------------------------------------------------------
    !set up dimensions for buffer netcdf output 
    ndims = 4
    ALLOCATE(dim_list(1:ndims),STAT=errorcode)
    IF (errorcode /= 0 ) CALL logging%error('Cant allocate array dim_list', __FILE__, __LINE__)
    ALLOCATE(time(1:ntime_aot),STAT=errorcode)
    IF (errorcode /= 0 ) CALL logging%error('Cant allocate array time', __FILE__, __LINE__)
    DO n = 1, ntime_aot
      CALL set_date_mm_extpar_field(n,dataDate,dataTime)
      time(n) = REAL(dataDate,wp) + REAL(dataTime,wp)/10000. ! units = "day as %Y%m%d.%f"
    ENDDO
    dim_list(1)%dimname = 'cell'
    dim_list(1)%dimsize = icon_grid%ncell 
    dim_list(2)%dimname = 'nv' ! icon_grid%nvertex_per_cell
    dim_list(2)%dimsize = icon_grid%nvertex_per_cell 
    dim_list(3)%dimname = 'vertex'
    dim_list(3)%dimsize = icon_grid%nvertex
    dim_list(4)%dimname = 'time'
    dim_list(4)%dimsize = ntime_aot

    dim_1d_icon   =  dim_list(1) ! cell

    ! set Icon coordinates for output
    CALL allocate_icon_coor(icon_grid%ncell, icon_grid%nvertex_per_cell)

    clon(:) = icon_grid_region%cells%center(:)%lon
    clat(:) = icon_grid_region%cells%center(:)%lat

    DO nc = 1, icon_grid%ncell
      DO nv = 1, icon_grid%nvertex_per_cell
        vert_id =  icon_grid_region%cells%vertex_index(nc,nv)
        clon_vertices(nv,nc) =  icon_grid_region%verts%vertex(vert_id)%lon
        clat_vertices(nv,nc) =  icon_grid_region%verts%vertex(vert_id)%lat
      ENDDO
    ENDDO

    ! define global attributes
    CALL set_global_att_icon(icon_grid,global_attributes,itopo_type,name_lookup_table_lu,lu_dataset,isoil_data)

! need to overload ...    
!    DO nc = 1, nglob_atts
!      CALL logging%info('NetCDF attributes: '//TRIM(global_attributes(nc)), __FILE__, __LINE__)
!    END DO

    !set up dimensions for buffer
    CALL  def_dimension_info_buffer(tg)
    ! dim_3d_tg

    !set up dimensions for ICON grid
    CALL def_dimension_info_icon(icon_grid)
    ! dim_icon

    ! define meta information for various land use related variables (GLC2000) for netcdf output
    CALL def_isa_fields_meta(dim_1d_icon)

    ! define meta information for various land use related variables for netcdf output
    CALL def_lu_fields_meta(nclass_lu,dim_1d_icon,lu_dataset=lu_dataset)
    ! dim_lu_tg
    ! fr_land_lu_meta, lu_tot_npixel_meta, &
    !  &       lu_class_fraction_meta, lu_class_npixel_meta, &
    !  &       ice_lu_meta, z0_lu_meta, &
    !  &       plcov_mx_lu_meta, plcov_mn_lu_meta, &
    !  &       lai_mx_lu_meta, lai_mn_lu_meta, &
    !  &       rs_min_lu_meta, urban_lu_meta, &
    !  &       for_d_lu_meta, for_e_lu_meta, &
    !  &       emissivity_lu_meta, root_lu_meta

    ! define meta information for target field variables lon_geo, lat_geo 
    CALL def_com_target_fields_meta(dim_1d_icon)
    ! lon_geo_meta and lat_geo_meta

    ! set mapping parameters for netcdf
    grid_mapping="lon_lat_on_sphere"

    CALL set_nc_grid_def_icon(grid_mapping)
    ! nc_grid_def_icon

    CALL def_soil_meta(dim_1d_icon, isoil_data)
    !  fr_land_soil_meta, soiltype_fao_meta

    CALL def_alb_meta(ntime_alb,dim_1d_icon)

    !define meta information for various NDVI data related variables for netcdf output
    CALL def_ahf_meta(dim_1d_icon)
    ! dim_ahf_tg, ahf_field_meta

    !define meta information for various NDVI data related variables for netcdf output
    CALL def_ndvi_meta(ntime_ndvi,dim_1d_icon)
    ! dim_ndvi_tg, ndvi_max_meta, ndvi_field_mom_meta, ndvi_ratio_mom_meta

    !define meta information for various EMISS data related variables for netcdf output
    CALL def_emiss_meta(ntime_emiss,dim_1d_icon)
    ! dim_emiss_tg, emiss_max_meta, emiss_field_mom_meta, emiss_ratio_mom_meta

    CALL def_era_meta(ntime_ndvi,dim_1d_icon)

    ! define meta information for various TOPO data related variables for netcdf output
    CALL def_topo_meta(dim_1d_icon,itopo_type)

    !  hh_topo_meta, fr_land_topo_meta, &
    !         stdh_topo_meta, theta_topo_meta, &
    !         aniso_topo_meta, slope_topo_meta, &
    !         hh_vert_meta, npixel_vert_meta
    !\TODO HA: this is a "quick fix" for ICON, find a better solution
    CALL hh_topo_meta%overwrite_varname('topography_c')

    CALL def_topo_vertex_meta(icon_grid%nvertex)
    ! dim_buffer_vertex
    !  hh_vert_meta, npixel_vert_meta

    ! define dimensions and meta information for variable aot_tg for netcdf output
    CALL def_aot_tg_meta(ntime_aot,ntype_aot,dim_1d_icon)
    ! dim_aot_tg and aot_tg_meta
    ! dim_aot_ty, aer_bc_meta, aer_dust_meta, aer_org_meta, aer_so4_meta, aer_ss_meta

    ! define meta information for variable crutemp for netcdf output
    CALL def_crutemp_meta(dim_1d_icon)
    ! crutemp_meta

    ! define meta information for various land use related variables (FLAKE) for netcdf output
    CALL def_flake_fields_meta(dim_1d_icon)
    ! lake_depth_meta, fr_lake_meta, &
    !  &       flake_tot_npixel_meta


    ! ** provisional fix for ICON netcdf variable names in order to ensure backward compatibility
    CALL ice_lu_meta%overwrite_varname('ICE')
    CALL alb_field_mom_meta%overwrite_varname('ALB')
    CALL alnid_field_mom_meta%overwrite_varname('ALNID')
    CALL aluvd_field_mom_meta%overwrite_varname('ALUVD')
    CALL aer_bc_meta%overwrite_varname('AER_BC')
    CALL aer_dust_meta%overwrite_varname('AER_DUST')
    CALL aer_org_meta%overwrite_varname('AER_ORG')
    CALL aer_so4_meta%overwrite_varname('AER_SO4')
    CALL aer_ss_meta%overwrite_varname('AER_SS')
    ! ** 


    !-----------------------------------------------------------------
   CALL open_new_netcdf_file(netcdf_filename=TRIM(netcdf_filename),   &
         &                       dim_list=dim_list,                  &
         &                       global_attributes=global_attributes, &
         &                       time=time,                           &
         &                       ncid=ncid)
    !-----------------------------------------------------------------
    ! soiltype
!    write(0,*) trim(soiltype_fao_meta%varname)
!    CALL netcdf_put_var(ncid,soiltype_fao(1:icon_grid%ncell,1,1),soiltype_fao_meta,undef_int)

    ! soiltype_deep
    IF (ldeep_soil) THEN
      CALL netcdf_put_var(ncid,soiltype_deep(1:icon_grid%ncell,1,1),soiltype_FAO_deep_meta,undef_int)
    ENDIF

    ! fr_sand
    IF (isoil_data == HWSD_data) THEN
      CALL netcdf_put_var(ncid,fr_sand(1:icon_grid%ncell,1,1),HWSD_SAND_meta,undefined)
    ENDIF

    ! fr_silt
    IF (isoil_data == HWSD_data) THEN
      CALL netcdf_put_var(ncid,fr_silt(1:icon_grid%ncell,1,1),HWSD_SILT_meta,undefined)
    ENDIF

    ! fr_clay
    IF (isoil_data == HWSD_data) THEN
      CALL netcdf_put_var(ncid,fr_clay(1:icon_grid%ncell,1,1),HWSD_CLAY_meta,undefined)
    ENDIF

    ! fr_oc
    IF (isoil_data == HWSD_data) THEN
      CALL netcdf_put_var(ncid,fr_oc(1:icon_grid%ncell,1,1),HWSD_OC_meta,undefined)
    ENDIF

    ! fr_bd
    IF (isoil_data == HWSD_data) THEN
      CALL netcdf_put_var(ncid,fr_bd(1:icon_grid%ncell,1,1),HWSD_BD_meta,undefined)
    ENDIF

    ! fr_sand_deep
    IF (ldeep_soil) THEN
      CALL netcdf_put_var(ncid,fr_sand_deep(1:icon_grid%ncell,1,1),HWSD_SAND_deep_meta,undefined)
    ENDIF

    ! fr_silt_deep
    IF (ldeep_soil) THEN
      CALL netcdf_put_var(ncid,fr_silt_deep(1:icon_grid%ncell,1,1),HWSD_SILT_deep_meta,undefined)
    ENDIF

    ! fr_clay_deep
    IF (ldeep_soil) THEN
      CALL netcdf_put_var(ncid,fr_clay_deep(1:icon_grid%ncell,1,1),HWSD_CLAY_deep_meta,undefined)
    ENDIF

    ! fr_oc_deep
    IF (ldeep_soil) THEN
      CALL netcdf_put_var(ncid,fr_oc_deep(1:icon_grid%ncell,1,1),HWSD_OC_deep_meta,undefined)
    ENDIF

    ! fr_bd_deep
    IF (ldeep_soil) THEN
      CALL netcdf_put_var(ncid,fr_bd_deep(1:icon_grid%ncell,1,1),HWSD_BD_deep_meta,undefined)
    ENDIF

    !-----------------------------------------------------------------
    ! soiltype  -> Integer Field!!
    CALL netcdf_put_var(ncid,soiltype_fao(1:icon_grid%ncell,1,1),soiltype_fao_meta,undef_int)

    n=1 ! fr_land_lu
    CALL netcdf_put_var(ncid,fr_land_lu(1:icon_grid%ncell,1,1),fr_land_lu_meta,undefined)

    n=2 ! ice_lu
    CALL netcdf_put_var(ncid,ice_lu(1:icon_grid%ncell,1,1),ice_lu_meta,undefined)

    n=3 ! plcov_mx_lu
    CALL netcdf_put_var(ncid,plcov_mx_lu(1:icon_grid%ncell,1,1),plcov_mx_lu_meta,undefined)

    n=4 ! lai_mx_lu
    CALL netcdf_put_var(ncid,lai_mx_lu(1:icon_grid%ncell,1,1),lai_mx_lu_meta,undefined)

    n=5 ! rs_min_lu
    CALL netcdf_put_var(ncid,rs_min_lu(1:icon_grid%ncell,1,1),rs_min_lu_meta,undefined)

    n=6 ! urban_lu
    CALL netcdf_put_var(ncid,urban_lu(1:icon_grid%ncell,1,1),urban_lu_meta,undefined)

    n=7 ! for_d_lu
    CALL netcdf_put_var(ncid,for_d_lu(1:icon_grid%ncell,1,1),for_d_lu_meta,undefined)

    n=8 ! for_e_lu
    CALL netcdf_put_var(ncid,for_e_lu(1:icon_grid%ncell,1,1),for_e_lu_meta,undefined)

    n=9 ! skinc_lu
    CALL netcdf_put_var(ncid, skinc_lu(1:icon_grid%ncell,1,1),skinc_lu_meta,undefined)

    n=10 ! emissivity_lu
    CALL netcdf_put_var(ncid, emissivity_lu(1:icon_grid%ncell,1,1),emissivity_lu_meta,undefined)

    n=11 ! root_lu
    CALL netcdf_put_var(ncid,root_lu(1:icon_grid%ncell,1,1),root_lu_meta,undefined)

    n=12 ! z0_lu
    CALL netcdf_put_var(ncid,z0_lu(1:icon_grid%ncell,1,1),z0_lu_meta,undefined)

    n=13 ! lon
    CALL netcdf_put_var(ncid,lon_geo(1:icon_grid%ncell,1,1),lon_geo_meta,undefined)

    n=14 ! lat
    CALL netcdf_put_var(ncid,lat_geo(1:icon_grid%ncell,1,1),lat_geo_meta,undefined)

    n=15 ! ndvi_max
    CALL netcdf_put_var(ncid,ndvi_max(1:icon_grid%ncell,1,1),ndvi_max_meta,undefined)

    n=16 ! hh_topo
    CALL netcdf_put_var(ncid,hh_topo(1:icon_grid%ncell,1,1),hh_topo_meta,undefined)

    n=17 ! hh_topo
    CALL netcdf_put_var(ncid,hh_topo_max(1:icon_grid%ncell,1,1),hh_topo_max_meta,undefined)

    n=18 ! hh_topo
    CALL netcdf_put_var(ncid,hh_topo_min(1:icon_grid%ncell,1,1),hh_topo_min_meta,undefined)
    
    n=19 ! stdh_topo
    CALL netcdf_put_var(ncid,stdh_topo(1:icon_grid%ncell,1,1),stdh_topo_meta,undefined)

    IF (lsso) THEN
      n=20 ! theta_topo
      CALL netcdf_put_var(ncid,theta_topo(1:icon_grid%ncell,1,1),theta_topo_meta,undefined)
    ENDIF

    IF (lsso) THEN
      n=21 ! aniso_topo
      CALL netcdf_put_var(ncid,aniso_topo(1:icon_grid%ncell,1,1),aniso_topo_meta,undefined)
    ENDIF

    IF (lsso) THEN
      n=22 ! slope_topo
      CALL netcdf_put_var(ncid,slope_topo(1:icon_grid%ncell,1,1),slope_topo_meta,undefined)
    ENDIF

    n=23 ! crutemp
    CALL netcdf_put_var(ncid,crutemp(1:icon_grid%ncell,1,1),crutemp_meta,undefined)

    n=24 ! fr_lake
    CALL netcdf_put_var(ncid,fr_lake(1:icon_grid%ncell,1,1),fr_lake_meta,undefined)

    n=25 ! lake_depth
    CALL netcdf_put_var(ncid,lake_depth(1:icon_grid%ncell,1,1),lake_depth_meta,undefined)

    IF (l_use_ahf) THEN
      n=26 ! ahf
      CALL netcdf_put_var(ncid,ahf_field(1:icon_grid%ncell,1,1),ahf_field_meta,undefined)
    END IF


    IF (l_use_isa) THEN
      n=27 ! isa
      CALL netcdf_put_var(ncid,isa_field(1:icon_grid%ncell,1,1),isa_field_meta,undefined)
    END IF


    ! hh_vert not demanded for output 

    n=1 ! lu_class_fraction

    CALL netcdf_put_var(ncid,lu_class_fraction(1:icon_grid%ncell,1,1,1:nclass_lu),&
         &                 lu_class_fraction_meta,undefined)

    n=2 ! ndvi_field_mom
    CALL netcdf_put_var(ncid,ndvi_field_mom(1:icon_grid%ncell,1,1,1:ntime_ndvi), &
         &                 ndvi_field_mom_meta, undefined)

    n=3 ! ndvi_ratio_mom
    CALL netcdf_put_var(ncid,ndvi_ratio_mom(1:icon_grid%ncell,1,1,1:ntime_ndvi), &
         &                 ndvi_ratio_mom_meta, undefined)

  namelist_file="INPUT_EMISS"
  INQUIRE(FILE=TRIM(namelist_file), EXIST=l_use_emiss)
  if (l_use_emiss) then
    n=4 ! emiss_field_mom
    CALL netcdf_put_var(ncid,emiss_field_mom(1:icon_grid%ncell,1,1,1:ntime_emiss), &
         &                 emiss_field_mom_meta, undefined)
  end if

    n=1 ! aot_bc
    CALL netcdf_put_var(ncid,aot_tg(1:icon_grid%ncell,1,1,1,1:ntime_aot), &
         &                 aer_bc_meta, undefined)
    n=2 ! aot_dust
    CALL netcdf_put_var(ncid,aot_tg(1:icon_grid%ncell,1,1,2,1:ntime_aot), &
         &                 aer_dust_meta, undefined)

    n=3 ! aot_org
    CALL netcdf_put_var(ncid,aot_tg(1:icon_grid%ncell,1,1,3,1:ntime_aot), &
         &                 aer_org_meta, undefined)

    n=4 ! aot_so4
    CALL netcdf_put_var(ncid,aot_tg(1:icon_grid%ncell,1,1,4,1:ntime_aot), &
         &                 aer_so4_meta, undefined)

    n=5 ! aot_ss
    CALL netcdf_put_var(ncid,aot_tg(1:icon_grid%ncell,1,1,5,1:ntime_aot), &
         &                 aer_ss_meta, undefined)

    n=6 ! alb_field_mom
    CALL netcdf_put_var(ncid,alb_field_mom(1:icon_grid%ncell,1,1,1:ntime_alb), &
         &                 alb_field_mom_meta, undefined)

    n=7 ! alnid_field_mom
    CALL netcdf_put_var(ncid,alnid_field_mom(1:icon_grid%ncell,1,1,1:ntime_alb), &
         &                 alnid_field_mom_meta, undefined)

    n=8 ! aluvd_field_mom
    CALL netcdf_put_var(ncid,aluvd_field_mom(1:icon_grid%ncell,1,1,1:ntime_alb), &
         &                 aluvd_field_mom_meta, undefined)


    n=9 ! sst_field
    CALL netcdf_put_var(ncid,sst_field(1:icon_grid%ncell,1,1,1:ntime_ndvi), &
         &                 sst_field_meta, undefined)

    n=10 ! wsnow_field
    CALL netcdf_put_var(ncid,wsnow_field(1:icon_grid%ncell,1,1,1:ntime_ndvi), &
         &                 wsnow_field_meta, undefined)

    n=11 ! t2m_field
    CALL netcdf_put_var(ncid,t2m_field(1:icon_grid%ncell,1,1,1:ntime_ndvi), &
         &                 t2m_field_meta, undefined)

    n=12 ! hsurf_field
    CALL netcdf_put_var(ncid,hsurf_field(1:icon_grid%ncell,1,1), &
         &                 hsurf_field_meta, undefined)


    ! write out ICON grid cell coordinates (in radians) to netcdf file
    CALL  netcdf_put_var(ncid,clon,clon_meta,undefined)
    CALL  netcdf_put_var(ncid,clat,clat_meta,undefined)
!!$     CALL  netcdf_put_var(ncid,clon_vertices,clon_vertices_meta,undefined)
!!$     CALL  netcdf_put_var(ncid,clat_vertices,clat_vertices_meta,undefined)


    !-----------------------------------------------------------------

    CALL close_netcdf_file(ncid)
    
  END SUBROUTINE write_netcdf_icon_grid_extpar
  !-----------------------------------------------------------------------

  !----------------------------------------------------------------------- 
  !-----------------------------------------------------------------------
  !> set global attributes for netcdf with lu data
  SUBROUTINE set_global_att_icon(icon_grid,global_attributes,itopo_type,name_lookup_table_lu,lu_dataset,isoil_data)



    TYPE(netcdf_attributes), INTENT(INOUT) :: global_attributes(1:9)
    TYPE(icosahedral_triangular_grid), INTENT(IN) :: icon_grid !< structure which contains the definition of the ICON grid
    INTEGER (KIND=i4),     INTENT(IN) :: itopo_type,isoil_data
    INTEGER env_len, status
    CHARACTER (LEN=*),INTENT(IN) :: name_lookup_table_lu
    CHARACTER (LEN=*),INTENT(IN) :: lu_dataset
    CHARACTER (LEN=filename_max) :: env_str

    !local variables
    CHARACTER(len=4)  :: number_Of_Grid_Used_string
    CHARACTER(len=10) :: ydate
    CHARACTER(len=10) :: ytime
    CHARACTER(len=2)  :: cc
    CHARACTER(len=2)  :: yy
    CHARACTER(len=2)  :: mm
    CHARACTER(len=2)  :: dd
    CHARACTER(len=2)  :: hh
    CHARACTER(len=2)  :: minute
    CHARACTER(len=1 ) :: uuid(16)    !   UUID of unstructured grids 

    ! define global attributes

    global_attributes(1)%attname = 'title'
    global_attributes(1)%attributetext='external parameter'

    global_attributes(2)%attname = 'institution'
    global_attributes(2)%attributetext='Deutscher Wetterdienst'

    global_attributes(3)%attname = 'rawdata'
    SELECT CASE(itopo_type)
    CASE(1) ! topo_gl
      global_attributes(3)%attributetext=TRIM(lu_dataset)//', FAO DSMW, GLOBE, Lake Database'
      IF(isoil_data == HWSD_data) global_attributes(3)%attributetext=TRIM(lu_dataset)//', HWSD, GLOBE, Lake Database'
    CASE(2) ! topo_aster
      global_attributes(3)%attributetext=TRIM(lu_dataset)//', FAO DSMW, ASTER, Lake Database'
      IF(isoil_data == HWSD_data) global_attributes(3)%attributetext=TRIM(lu_dataset)//', HWSD, ASTER, Lake Database'   
    END SELECT

    global_attributes(4)%attname = 'note'
    global_attributes(4)%attributetext='Landuse data look-up table: '//TRIM(name_lookup_table_lu)

    CALL DATE_AND_TIME(ydate,ytime)
    READ(ydate,'(4A2)') cc,yy,mm,dd
    READ(ytime,'(2A2)') hh, minute

    ydate=TRIM(cc)//TRIM(yy)//'-'//TRIM(mm)//'-'//TRIM(dd)
    ytime=TRIM(hh)//':'//TRIM(minute) 

    global_attributes(5)%attname = 'history'
    global_attributes(5)%attributetext=TRIM(ydate)//'T'//TRIM(ytime)//' extpar_consistency_check'

    global_attributes(6)%attname = 'version'
    global_attributes(6)%attributetext = TRIM(INFO_PackageName)

    global_attributes(7)%attname = 'comment'
    CALL get_environment_VARIABLE( "progdir", env_str, env_len, status)
    global_attributes(7)%attributetext='binaries in '//TRIM(env_str)

    WRITE(number_Of_Grid_Used_string,'(I4)')  icon_grid%number_Of_Grid_Used
    global_attributes(8)%attname = 'number_of_grid_used'
    global_attributes(8)%attributetext=number_Of_Grid_Used_string

    CALL decode_uuid (icon_grid%uuidOfHGrid, uuid) 

    global_attributes(9)%attname = 'uuidOfHGrid'
    global_attributes(9)%attributetext=icon_grid%uuidOfHGrid

    CLOSE(1)

  END SUBROUTINE set_global_att_icon
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !> set global attributes for netcdf with lu data
  SUBROUTINE set_global_att_extpar(global_attributes,name_lookup_table_lu,lu_dataset,isoil_data,lscale_separation,y_orofilt)

    TYPE(netcdf_attributes), INTENT(INOUT) :: global_attributes(1:11)
    CHARACTER (LEN=*),INTENT(IN) :: name_lookup_table_lu
    CHARACTER (LEN=*),INTENT(IN) :: lu_dataset
    INTEGER,          INTENT(IN) :: isoil_data
    LOGICAL,          INTENT(IN) :: lscale_separation
    CHARACTER (LEN=*),INTENT(IN) :: y_orofilt
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
    global_attributes(1)%attributetext='external parameter'

    global_attributes(2)%attname = 'institution'
    global_attributes(2)%attributetext='COSMO Consortium and CLM Community'
    !    global_attributes(2)%attributetext='Deutscher Wetterdienst'

    global_attributes(3)%attname = 'rawdata'
    SELECT CASE(itopo_type)
    CASE(topo_aster)
      IF (isoil_data >= HWSD_data) THEN
        global_attributes(3)%attributetext=TRIM(lu_dataset)//', HWSD, ASTER, Lake Database'
      ELSE
        global_attributes(3)%attributetext=TRIM(lu_dataset)//', FAO DSMW, ASTER, Lake Database'
      ENDIF
    CASE(topo_gl)
      IF (isoil_data >= HWSD_data) THEN
        global_attributes(3)%attributetext=TRIM(lu_dataset)//', HWSD, GLOBE, Lake Database'
      ELSE
        global_attributes(3)%attributetext=TRIM(lu_dataset)//', FAO DSMW, GLOBE, Lake Database'
      ENDIF
    END SELECT
    global_attributes(4)%attname = 'note'
    global_attributes(4)%attributetext='Landuse data look-up table: '//TRIM(name_lookup_table_lu)

    global_attributes(5)%attname = 'scale_separation'
    IF (lscale_separation) THEN
      global_attributes(5)%attributetext = 'scale separation: ON'
    ELSE 
      global_attributes(5)%attributetext = 'scale separation: OFF'
    ENDIF

    global_attributes(6)%attname = 'filter_options'
    global_attributes(6)%attributetext = TRIM(y_orofilt)

    CALL DATE_AND_TIME(ydate,ytime)
    READ(ydate,'(4A2)') cc,yy,mm,dd
    READ(ytime,'(2A2)') hh, minute

    ydate=TRIM(cc)//TRIM(yy)//'-'//TRIM(mm)//'-'//TRIM(dd)
    ytime=TRIM(hh)//':'//TRIM(minute) 

    global_attributes(7)%attname = 'history'
    global_attributes(7)%attributetext=TRIM(ydate)//'T'//TRIM(ytime)//' extpar_consistency_check'

    global_attributes(8)%attname = 'comment'
    global_attributes(8)%attributetext='Generation of external parameters '// &
#ifdef CLM
         &                             'for COSMO-CLM through WebPEP'
#else
         &                             'for numerical atmospheric models COSMO and ICON.'
#endif

    global_attributes(9)%attname = 'version'
    global_attributes(9)%attributetext = TRIM(INFO_PackageName)

    global_attributes(10)%attname = 'Revision Hash'
    global_attributes(10)%attributetext = TRIM(INFO_RevisionHash)//" ("//TRIM(INFO_CodeIsModified)//")"

    global_attributes(11)%attname = 'Conventions'
    global_attributes(11)%attributetext = 'CF-1.5'


  END SUBROUTINE set_global_att_extpar
  !-----------------------------------------------------------------------


  SUBROUTINE decode_uuid (uuid_str, uuid)
    CHARACTER(LEN=*), INTENT(IN)  :: uuid_str   ! uuid encoded as string
    CHARACTER(LEN=1), INTENT(OUT) :: uuid(:)    ! decoded uuid

    INTEGER          :: i, j, l, n, b
    CHARACTER(LEN=2) :: buf

    uuid(:) = ACHAR (0)
    l = VERIFY (uuid_str, "0123456789ABCDEFabcdef-")
    IF (l > 0) THEN
      WRITE (message_text,*) "Warning: invalid character in uuid: '", uuid_str(l:l),"'"
      CALL logging%error(message_text, __FILE__, __LINE__)
    END IF
    n = LEN  (uuid_str)
    i = 1
    j = 0
    DO WHILE (i < n)
      buf = uuid_str(i:i+1)
      IF (buf(1:1) == "-") THEN
        i = i + 1                     ! Skip over dashes
        CYCLE
      END IF
      i = i + 2
      READ (buf,'(Z2)') b
      j = j + 1
      IF (j > SIZE (uuid)) CALL logging%error("uuid input too long!", __FILE__, __LINE__)
      uuid(j) = ACHAR (b)
    END DO
    IF (i == n) CALL logging%error("uuid bad length", __FILE__, __LINE__)
  END SUBROUTINE decode_uuid

END MODULE mo_extpar_output_nc
