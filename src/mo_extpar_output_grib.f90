!+ Fortran module with grib output routines for external parameters
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_0         2010/12/21 Hermann Asensio
!  Initial release
! V1_1         2011/01/20 Hermann Asensio
! small bug fixes accroding to Fortran compiler warnings 
! V1_3         2011/04/19 Hermann Asensio
! add support for GRIB1 and GRIB2
! V1_8         2013/03/12 Frank Brenner
!  introduced MODIS albedo dataset(s) as new external parameter(s) 
! V2_0_3       2015-01-12 Juergen Helmert
!  Bugfix correction covers CSCS SVN r5907-r6359
! V3_0         2015-05-18 Juergen Helmert 
!  Remove tableVersion from GRIB2-Output use Template default  
! V4_0         2016/08/23 Daniel Lüthi
!  added support for MACv2 spectrally stratified monthly aerosol fields       
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!> Fortran module with grib output routines for external parameters 
!> ouptut routines
!> \author Hermann Asensio
MODULE mo_extpar_output_grib


  !> kind parameters are defined in MODULE data_parameters
  USE mo_kind, ONLY: wp
  USE mo_kind, ONLY: i8
  USE mo_kind, ONLY: i4

  !> data type structures form module GRID_structures
  USE mo_grid_structures, ONLY: reg_lonlat_grid
  USE mo_grid_structures, ONLY: rotated_lonlat_grid
  USE mo_grid_structures, ONLY: icosahedral_triangular_grid
  USE mo_grid_structures, ONLY: gme_triangular_grid
  USE mo_grid_structures, ONLY: target_grid_def

  USE mo_io_utilities, ONLY: var_meta_info

  USE mo_io_utilities, ONLY: dim_meta_info

  USE mo_io_utilities, ONLY: vartype_int 
  USE mo_io_utilities, ONLY: vartype_real
  USE mo_io_utilities, ONLY: vartype_char

  USE mo_io_utilities, ONLY: get_date_const_field, set_date_mm_extpar_field

  USE mo_exception,         ONLY: message_text, message, finish

  USE mo_io_units,          ONLY: filename_max


  !> abort_extpar defined in MODULE utilities_extpar
  USE mo_utilities_extpar, ONLY: abort_extpar
  USE mo_utilities_extpar, ONLY: get_rot_spol_coor


  USE mo_glc2000_lookup_tables, ONLY: nclass_glc2000
  USE mo_ndvi_data, ONLY: ntime_ndvi
  USE mo_aot_data, ONLY: ntype_aot, ntime_aot, nspb_aot, iaot_type
  USE mo_albedo_data, ONLY: ntime_alb

  USE mo_soil_data,   ONLY: HWSD_data

  USE grib_api 
  USE mo_io_grib_api
  USE mo_var_meta_data, ONLY: aer_bc_meta
  USE mo_var_meta_data, ONLY: aer_dust_meta
  USE mo_var_meta_data, ONLY: aer_org_meta
  USE mo_var_meta_data, ONLY: aer_so4_meta
  USE mo_var_meta_data, ONLY: aer_ss_meta



  IMPLICIT NONE

  PRIVATE:: decode_uuid

  PUBLIC :: write_cosmo_grid_extpar_grib
  PUBLIC :: write_gme_grid_extpar_grib



  CONTAINS


  !> grib output of external parameters for COSMO
  SUBROUTINE write_cosmo_grid_extpar_grib(grib_filename, &
    &                                     grib_sample,   &
    &                                     cosmo_grid,    &
    &                                     tg,            &
    &                                     isoil_data,    &
    &                                     ldeep_soil,    &
    &                                     itopo_type,    &
    &                                     lsso,          &
    &                                     lrad,          &
    &                                     nhori,         &
    &                                     undefined,     &
    &                                     undef_int,     &
    &                                     lon_geo,       &
    &                                     lat_geo,       &   
    &                                     fr_land_lu,    & ! &                                     ice_lu, &
    &                                     z0_lu,         &
    &                                     root_lu,       &
    &                                     plcov_mn_lu,   &
    &                                     plcov_mx_lu,   &
    &                                     lai_mn_lu,     &
    &                                     lai_mx_lu,     &
    &                                     rs_min_lu,     &
    &                                     urban_lu,      &
    &                                     for_d_lu,      &
    &                                     for_e_lu,      &
    &                                     emissivity_lu, &
    &                                     lake_depth,    &
    &                                     fr_lake,       &
    &                                     soiltype_fao,  &
    &                                     ndvi_max,      &
    &                                     ndvi_field_mom,&
    &                                     ndvi_ratio_mom,&
    &                                     hh_topo,       &
    &                                     stdh_topo,     &
    &                                     aot_tg,        &
    &                                     MAC_aot_tg,    &
    &                                     MAC_ssa_tg,    &
    &                                     MAC_asy_tg,    & 
    &                                     crutemp,       &
    &                                     alb_field_mom, &
    &                                     fr_sand,       & 
    &                                     fr_silt,       &
    &                                     fr_clay,       &
    &                                     fr_oc,         &
    &                                     fr_bd,         &
    &                                     theta_topo,    &
    &                                     aniso_topo,    &
    &                                     slope_topo,    &
    &                                     slope_asp_topo,&
    &                                     slope_ang_topo,&
    &                                     horizon_topo,  &
    &                                     skyview_topo)
 
  USE mo_var_meta_data, ONLY: lon_geo_meta, &
    &                         lat_geo_meta, &
    &                         no_raw_data_pixel_meta, &
    &                         def_com_target_fields_meta  
   
 
  USE mo_cosmo_grid, ONLY: lon_rot, lat_rot

      
  USE mo_var_meta_data, ONLY: dim_rlon_cosmo, &
    &                         dim_rlat_cosmo, &
    &                         dim_2d_cosmo,   &
    &                         rlon_meta,      &
    &                         rlat_meta,      &
    &                         def_dimension_info_cosmo


  USE mo_var_meta_data, ONLY: dim_3d_tg, &
    &                         dim_4d_tg, &
    &                         def_dimension_info_buffer


  USE mo_var_meta_data, ONLY: def_lu_fields_meta

  USE mo_var_meta_data, ONLY: dim_lu_tg

  USE mo_var_meta_data, ONLY: fr_land_lu_meta, lu_tot_npixel_meta, &
    &       lu_class_fraction_meta, lu_class_npixel_meta, &
    &       ice_lu_meta, z0_lu_meta, &
    &       plcov_mx_lu_meta, plcov_mn_lu_meta, &
    &       lai_mx_lu_meta, lai_mn_lu_meta, &
    &       rs_min_lu_meta, urban_lu_meta, &
    &       for_d_lu_meta, for_e_lu_meta, &
    &       emissivity_lu_meta, root_lu_meta

  USE mo_var_meta_data, ONLY: def_flake_fields_meta
  USE mo_var_meta_data, ONLY: lake_depth_meta, fr_lake_meta

  USE mo_var_meta_data, ONLY: def_soil_meta
  USE mo_var_meta_data, ONLY: fr_land_soil_meta, soiltype_fao_meta,&
      &                       HWSD_SAND_META,HWSD_SILT_META,       &
      &                       HWSD_CLAY_META,HWSD_OC_META,         &
      &                       HWSD_BD_META,HWSD_DM_META

  USE mo_var_meta_data, ONLY: dim_alb_tg
  USE mo_var_meta_data, ONLY: alb_field_mom_meta, &
      &                       def_alb_meta,       &
      &                       alb_interpol_meta
  
  USE mo_var_meta_data, ONLY: dim_ndvi_tg
  USE mo_var_meta_data, ONLY: ndvi_max_meta,       &
      &                       ndvi_field_mom_meta, &
      &                       ndvi_ratio_mom_meta, &
      &                       def_ndvi_meta

 USE mo_var_meta_data, ONLY: def_topo_meta, def_topo_vertex_meta
 USE mo_var_meta_data, ONLY: dim_buffer_cell, dim_buffer_vertex

 USE mo_var_meta_data, ONLY: hh_topo_meta, fr_land_topo_meta,        &
   &                         stdh_topo_meta, theta_topo_meta,        &
   &                         aniso_topo_meta, slope_topo_meta,       &
   &                         hh_vert_meta, npixel_vert_meta,           &
   &                         hh_fis_meta, npixel_vert_meta,            &
   &                         hh_fis_meta, slope_asp_topo_meta,        & 
   &                         slope_ang_topo_meta, horizon_topo_meta, &
   &                         skyview_topo_meta

 USE mo_var_meta_data, ONLY: dim_aot_tg,  &
   &                         aot_tg_meta, &
   &                         def_aot_tg_meta
 USE mo_var_meta_data, ONLY: aot_type_shortname

 USE mo_var_meta_data, ONLY: aot_tg_MAC_meta,&
   &                         ssa_tg_MAC_meta,&
   &                         asy_tg_MAC_meta

 USE mo_var_meta_data, ONLY: crutemp_meta, &
    &                        def_crutemp_meta

 USE mo_physical_constants, ONLY: grav


  CHARACTER (len=*), INTENT(IN) :: grib_filename !< filename for the grib file
  CHARACTER (len=*), INTENT(IN) :: grib_sample  !< name for grib sample  (sample to be found in $GRIB_SAMPLES_PATH)

  TYPE(rotated_lonlat_grid), INTENT(IN) :: cosmo_grid !< structure which contains the definition of the COSMO grid
  TYPE(target_grid_def), INTENT(IN) :: tg !< structure with target grid description
  INTEGER (KIND=i4), INTENT(IN)     :: isoil_data
  LOGICAL,         INTENT(IN)       :: ldeep_soil
  INTEGER (KIND=i4), INTENT(IN)     :: itopo_type
  LOGICAL,         INTENT(IN)       :: lsso  
  LOGICAL,         INTENT(IN)       :: lrad  
  INTEGER(KIND=i4),INTENT(IN)       :: nhori
  REAL(KIND=wp), INTENT(IN)         :: undefined       !< value to indicate undefined grid elements 
  INTEGER, INTENT(IN)               :: undef_int       !< value to indicate undefined grid elements
  REAL (KIND=wp), INTENT(IN) :: lon_geo(:,:,:)  !< longitude coordinates of the target grid in the geographical system
  REAL (KIND=wp), INTENT(IN) :: lat_geo(:,:,:)  !< latitude coordinates of the target grid in the geographical system
  REAL (KIND=wp), INTENT(IN)  :: fr_land_lu(:,:,:) !< fraction land due to lu raw data
  !REAL (KIND=wp), INTENT(IN)  :: ice_lu(:,:,:)     !< fraction of ice due to lu raw data
  REAL (KIND=wp), INTENT(IN)  :: z0_lu(:,:,:)      !< roughness length due to lu land use data
  REAL (KIND=wp), INTENT(IN)  :: root_lu(:,:,:)    !< root depth due to lu land use data
  REAL (KIND=wp), INTENT(IN)  :: plcov_mn_lu(:,:,:)!< plant cover maximum due to lu land use data
  REAL (KIND=wp), INTENT(IN)  :: plcov_mx_lu(:,:,:)!< plant cover minimum due to lu land use data
  REAL (KIND=wp), INTENT(IN)  :: lai_mn_lu(:,:,:)  !< Leaf Area Index maximum due to lu land use data
  REAL (KIND=wp), INTENT(IN)  :: lai_mx_lu(:,:,:)  !< Leaf Area Index minimum due to lu land use data
  REAL (KIND=wp), INTENT(IN)  :: rs_min_lu(:,:,:)  !< minimal stomata resistance due to lu land use data
  REAL (KIND=wp), INTENT(IN)  :: urban_lu(:,:,:)   !< urban fraction due to lu land use data
  REAL (KIND=wp), INTENT(IN)  :: for_d_lu(:,:,:)   !< deciduous forest (fraction) due to lu land use data
  REAL (KIND=wp), INTENT(IN)  :: for_e_lu(:,:,:)   !< evergreen forest (fraction) due to lu land use data
  REAL (KIND=wp), INTENT(IN)  :: emissivity_lu(:,:,:) !< longwave emissivity due to lu land use data

  REAL (KIND=wp), INTENT(IN)  :: lake_depth(:,:,:) !< lake depth
  REAL (KIND=wp), INTENT(IN)  :: fr_lake(:,:,:)     !< fraction of fresh water (lakes)

  INTEGER(KIND=i4), INTENT(IN) :: soiltype_fao(:,:,:) !< soiltype due to FAO Digital Soil map of the World

  REAL (KIND=wp), INTENT(IN) :: alb_field_mom(:,:,:,:) !< field for monthly mean albedo data 

  REAL (KIND=wp), INTENT(IN) :: ndvi_max(:,:,:) !< field for ndvi maximum
  REAL (KIND=wp), INTENT(IN) :: ndvi_field_mom(:,:,:,:) !< field for monthly mean ndvi data (12 months)
  REAL (KIND=wp), INTENT(IN) :: ndvi_ratio_mom(:,:,:,:) !< field for monthly ndvi ratio (12 months)

  REAL(KIND=wp), INTENT(IN)  :: hh_topo(:,:,:)  !< mean height 
  REAL(KIND=wp), INTENT(IN)  :: stdh_topo(:,:,:) !< standard deviation of subgrid scale orographic height
  REAL (KIND=wp), INTENT(IN)  :: aot_tg(:,:,:,:,:) !< aerosol optical thickness, aot_tg(ie,je,ke,ntype,ntime)
  REAL (KIND=wp), INTENT(IN)  :: MAC_aot_tg(:,:,:,:)
  REAL (KIND=wp), INTENT(IN)  :: MAC_ssa_tg(:,:,:,:)
  REAL (KIND=wp), INTENT(IN)  :: MAC_asy_tg(:,:,:,:)
  REAL(KIND=wp), INTENT(IN)  :: crutemp(:,:,:)  !< cru climatological temperature , crutemp(ie,je,ke) 
  REAL(KIND=wp), INTENT(IN), OPTIONAL :: fr_sand(:,:,:)   !< sand fraction due to HWSD
  REAL(KIND=wp), INTENT(IN), OPTIONAL :: fr_silt(:,:,:)   !< silt fraction due to HWSD
  REAL(KIND=wp), INTENT(IN), OPTIONAL :: fr_clay(:,:,:)   !< clay fraction due to HWSD
  REAL(KIND=wp), INTENT(IN), OPTIONAL :: fr_oc(:,:,:)   !< oc fraction due to HWSD
  REAL(KIND=wp), INTENT(IN), OPTIONAL :: fr_bd(:,:,:)   !< bulk density due to HWSD

  REAL(KIND=wp), INTENT(IN), OPTIONAL  :: theta_topo(:,:,:) !< sso parameter, angle of principal axis
  REAL(KIND=wp), INTENT(IN), OPTIONAL  :: aniso_topo(:,:,:) !< sso parameter, anisotropie factor
  REAL(KIND=wp), INTENT(IN), OPTIONAL  :: slope_topo(:,:,:) !< sso parameter, mean slope
  REAL(KIND=wp), INTENT(IN), OPTIONAL  :: slope_asp_topo(:,:,:)   !< lradtopo parameter, slope_aspect
  REAL(KIND=wp), INTENT(IN), OPTIONAL  :: slope_ang_topo(:,:,:)   !< lradtopo parameter, slope_angle
  REAL(KIND=wp), INTENT(IN), OPTIONAL  :: horizon_topo  (:,:,:,:) !< lradtopo parameter, horizon
  REAL(KIND=wp), INTENT(IN), OPTIONAL  :: skyview_topo  (:,:,:)   !< lradtopo parameter, skyview

  ! local variables

  INTEGER :: outfile_id
  INTEGER :: errorcode
  REAL (KIND=wp)  :: extpar_cosmo_buffer(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1:1) !< field to write out to GRIB file
  INTEGER :: gribid_dest !< id of grib message (GRIB_API)
  INTEGER (KIND=i8)  :: dataDate
  INTEGER (KIND=i8)  :: dataTime
  INTEGER :: mm ! month
  INTEGER :: ntype ! type of aerosol
  INTEGER :: ll ! level
  TYPE(var_meta_info) :: field_meta



  ! prepar meta data for output
  !set up dimensions for buffer
  CALL  def_dimension_info_buffer(tg,nhori=nhori)
  ! dim_3d_tg
  
  ! define meta information for various land use related variables (GLC2000) for netcdf output
  CALL def_lu_fields_meta(tg,nclass_glc2000,dim_3d_tg)
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
  CALL def_com_target_fields_meta(dim_3d_tg)
  ! lon_geo_meta and lat_geo_meta

  CALL def_soil_meta(dim_3d_tg,isoil_data)
  !  fr_land_soil_meta, soiltype_fao_meta

  CALL def_alb_meta(tg,ntime_alb,dim_3d_tg)

  !define meta information for various NDVI data related variables for netcdf output
  CALL def_ndvi_meta(tg,ntime_ndvi,dim_3d_tg)
  ! dim_ndvi_tg, ndvi_max_meta, ndvi_field_mom_meta, ndvi_ratio_mom_meta

  ! define meta information for various TOPO data related variables for netcdf output
  IF (lrad) THEN
  CALL def_topo_meta(dim_3d_tg,itopo_type,diminfohor=dim_4d_tg)
  !  hh_topo_meta, fr_land_topo_meta, &
  !         stdh_topo_meta, theta_topo_meta, &
  !         aniso_topo_meta, slope_topo_meta, &
  !         hh_vert_meta, npixel_vert_meta
  ELSE
    CALL def_topo_meta(dim_3d_tg,itopo_type)
    !  hh_topo_meta, fr_land_topo_meta, &
    !         stdh_topo_meta, theta_topo_meta, &
    !         aniso_topo_meta, slope_topo_meta, &
    !         hh_vert_meta, npixel_vert_meta
  ENDIF
  ! define dimensions and meta information for variable aot_tg for netcdf output
  CALL def_aot_tg_meta(tg,ntime_aot,ntype_aot,dim_3d_tg)
  ! dim_aot_tg and aot_tg_meta

  ! define meta information for variable crutemp for netcdf output
  CALL def_crutemp_meta(dim_3d_tg)
  ! crutemp_meta

   CALL def_flake_fields_meta(dim_3d_tg)
   ! lake_depth_meta, fr_lake_meta



  ! open a new grib file
  CALL grib_open_file(outfile_id,TRIM(grib_filename),'w')


  ! get dateDate according to DWD convention for invariant external parameter fields
  CALL get_date_const_field(dataDate,dataTime)

  ! write out fields to grib file
  PRINT *,'fr_land_lu_meta%shortName: ',fr_land_lu_meta%shortName

  CALL write_extpar_cosmo_field_grib(outfile_id,TRIM(grib_sample),&
  & cosmo_grid,fr_land_lu,fr_land_lu_meta,dataDate,dataTime)

  PRINT *,'hh_topo_meta%shortName: ',hh_topo_meta%shortName

  CALL write_extpar_cosmo_field_grib(outfile_id,TRIM(grib_sample),&
  & cosmo_grid,hh_topo,hh_topo_meta,dataDate,dataTime)
  
  ! output FIS as well
  PRINT *,'hh_fis_meta%shortName: ',hh_fis_meta%shortName
  CALL write_extpar_cosmo_field_grib(outfile_id,TRIM(grib_sample),&
  & cosmo_grid,hh_topo*grav,hh_fis_meta,dataDate,dataTime)

  CALL write_extpar_cosmo_field_grib(outfile_id,TRIM(grib_sample),&
  & cosmo_grid,stdh_topo,stdh_topo_meta,dataDate,dataTime)
  
  IF (lsso) THEN
  CALL write_extpar_cosmo_field_grib(outfile_id,TRIM(grib_sample),&
  & cosmo_grid,theta_topo,theta_topo_meta,dataDate,dataTime)
  ENDIF

  IF (lsso) THEN
  CALL write_extpar_cosmo_field_grib(outfile_id,TRIM(grib_sample),&
  & cosmo_grid,aniso_topo,aniso_topo_meta,dataDate,dataTime)
  ENDIF

  IF (lsso) THEN
  CALL write_extpar_cosmo_field_grib(outfile_id,TRIM(grib_sample),&
  & cosmo_grid,slope_topo,slope_topo_meta,dataDate,dataTime)
  ENDIF

  IF (lrad) THEN
    CALL write_extpar_cosmo_field_grib(outfile_id,TRIM(grib_sample),&
    & cosmo_grid,slope_asp_topo,slope_asp_topo_meta,dataDate,dataTime)
  ENDIF
  IF (lrad) THEN
    CALL write_extpar_cosmo_field_grib(outfile_id,TRIM(grib_sample),&
    & cosmo_grid,slope_ang_topo,slope_ang_topo_meta,dataDate,dataTime)
  ENDIF
!ROATODO: check 4D
  IF (lrad) THEN
    DO mm=1,nhori
      extpar_cosmo_buffer(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1:1) = &
      & horizon_topo(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1:1,mm) 
      CALL write_extpar_cosmo_real_1lev_grib(outfile_id,TRIM(grib_sample),&
      & cosmo_grid,extpar_cosmo_buffer,horizon_topo_meta,dataDate,dataTime,mm)
!      CALL write_extpar_cosmo_field_grib(outfile_id,TRIM(grib_sample),&
!      & cosmo_grid,extpar_cosmo_buffer,horizon_topo_meta%shortName,dataDate,dataTime)
    ENDDO
  ENDIF
  IF (lrad) THEN
    CALL write_extpar_cosmo_field_grib(outfile_id,TRIM(grib_sample),&
    & cosmo_grid,skyview_topo,skyview_topo_meta,dataDate,dataTime)
  ENDIF
  
  CALL write_extpar_cosmo_field_grib(outfile_id,TRIM(grib_sample),&
  & cosmo_grid,plcov_mn_lu,plcov_mn_lu_meta,dataDate,dataTime)

  CALL write_extpar_cosmo_field_grib(outfile_id,TRIM(grib_sample),&
  & cosmo_grid,plcov_mx_lu,plcov_mx_lu_meta,dataDate,dataTime)
  
  CALL write_extpar_cosmo_field_grib(outfile_id,TRIM(grib_sample), &
  & cosmo_grid,lai_mn_lu,lai_mn_lu_meta,dataDate,dataTime)

  CALL write_extpar_cosmo_field_grib(outfile_id,TRIM(grib_sample),&
  & cosmo_grid,lai_mx_lu,lai_mx_lu_meta,dataDate,dataTime)

  CALL write_extpar_cosmo_field_grib(outfile_id,TRIM(grib_sample), &
  & cosmo_grid,rs_min_lu,rs_min_lu_meta,dataDate,dataTime)
  
  CALL write_extpar_cosmo_field_grib(outfile_id,TRIM(grib_sample), &
  & cosmo_grid,for_e_lu,for_e_lu_meta,dataDate,dataTime)
  
  CALL write_extpar_cosmo_field_grib(outfile_id,TRIM(grib_sample), &
  & cosmo_grid,for_d_lu,for_d_lu_meta,dataDate,dataTime)
  
  CALL write_extpar_cosmo_field_grib(outfile_id,TRIM(grib_sample), &
  & cosmo_grid,emissivity_lu,emissivity_lu_meta,dataDate,dataTime)

  CALL write_extpar_cosmo_field_grib(outfile_id,TRIM(grib_sample),&
  & cosmo_grid,root_lu,root_lu_meta,dataDate,dataTime)

  CALL write_extpar_cosmo_field_grib(outfile_id,TRIM(grib_sample),&
  & cosmo_grid,z0_lu,z0_lu_meta,dataDate,dataTime)

  CALL write_extpar_cosmo_field_grib(outfile_id,TRIM(grib_sample),&
  & cosmo_grid,lake_depth,lake_depth_meta,dataDate,dataTime)
  
  CALL write_extpar_cosmo_field_grib(outfile_id,TRIM(grib_sample),&
  & cosmo_grid,fr_lake,fr_lake_meta,dataDate,dataTime)

  CALL write_extpar_cosmo_field_grib(outfile_id,TRIM(grib_sample),&
  & cosmo_grid,soiltype_fao,soiltype_fao_meta,dataDate,dataTime)

  IF(isoil_data == HWSD_data) THEN
  CALL write_extpar_cosmo_field_grib(outfile_id,TRIM(grib_sample),&
  & cosmo_grid,fr_sand,HWSD_SAND_meta,dataDate,dataTime)
  ENDIF

  IF(isoil_data == HWSD_data) THEN
  CALL write_extpar_cosmo_field_grib(outfile_id,TRIM(grib_sample),&
  & cosmo_grid,fr_silt,HWSD_SILT_meta,dataDate,dataTime)
  ENDIF

  IF(isoil_data == HWSD_data) THEN
  CALL write_extpar_cosmo_field_grib(outfile_id,TRIM(grib_sample),&
  & cosmo_grid,fr_clay,HWSD_CLAY_meta,dataDate,dataTime)
  ENDIF

  IF(isoil_data == HWSD_data) THEN
  CALL write_extpar_cosmo_field_grib(outfile_id,TRIM(grib_sample),&
  & cosmo_grid,fr_oc,HWSD_OC_meta,dataDate,dataTime)
  ENDIF

  IF(isoil_data == HWSD_data) THEN
  CALL write_extpar_cosmo_field_grib(outfile_id,TRIM(grib_sample),&
  & cosmo_grid,fr_bd,HWSD_BD_meta,dataDate,dataTime)
  ENDIF

  CALL write_extpar_cosmo_field_grib(outfile_id,TRIM(grib_sample),&
  & cosmo_grid,crutemp,crutemp_meta,dataDate,dataTime)
  
  CALL write_extpar_cosmo_field_grib(outfile_id,TRIM(grib_sample),&
  & cosmo_grid,lon_geo,lon_geo_meta,dataDate,dataTime)
  
  CALL write_extpar_cosmo_field_grib(outfile_id,TRIM(grib_sample),&
  & cosmo_grid,lat_geo,lat_geo_meta,dataDate,dataTime)

 ! urban_lu
 ! ndvi_max
 CALL write_extpar_cosmo_field_grib(outfile_id,TRIM(grib_sample),&
 & cosmo_grid,ndvi_max,ndvi_max_meta,dataDate,dataTime)

 ! ndvi_field_mom
 DO mm=1,12
    ! get dataDate and dataTime according to DWD convention for external parameters of climatological monthly mean fields
    CALL set_date_mm_extpar_field(mm,dataDate,dataTime)
    extpar_cosmo_buffer(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1:1) = &
  & ndvi_field_mom(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1:1,mm) 
    CALL write_extpar_cosmo_field_grib(outfile_id,TRIM(grib_sample),&
    & cosmo_grid,extpar_cosmo_buffer,ndvi_field_mom_meta,dataDate,dataTime)
 ENDDO

 ! ndvi_ratio_mom
 DO mm=1,12
    ! get dataDate and dataTime according to DWD convention for external parameters of climatological monthly mean fields
    CALL set_date_mm_extpar_field(mm,dataDate,dataTime)
    extpar_cosmo_buffer(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1:1) = &
  & ndvi_ratio_mom(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1:1,mm) 
    CALL write_extpar_cosmo_field_grib(outfile_id,TRIM(grib_sample),&
    & cosmo_grid,extpar_cosmo_buffer,ndvi_ratio_mom_meta,dataDate,dataTime)
 ENDDO

 DO mm=1,12
    CALL set_date_mm_extpar_field(mm,dataDate,dataTime)
    extpar_cosmo_buffer(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1:1) = &
  & alb_field_mom(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1:1,mm)
    CALL write_extpar_cosmo_field_grib(outfile_id,TRIM(grib_sample),&
    & cosmo_grid,extpar_cosmo_buffer,alb_field_mom_meta,dataDate,dataTime)
 ENDDO 

 IF (iaot_type == 4) THEN
   DO mm=1,12
     CALL set_date_mm_extpar_field(mm,dataDate,dataTime)
     DO ll=1,nspb_aot
       extpar_cosmo_buffer(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1:1) = &
       &  MAC_aot_tg(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,ll:ll,mm)
       CALL write_extpar_cosmo_real_1lev_grib(outfile_id,TRIM(grib_sample), &
       & cosmo_grid,extpar_cosmo_buffer,aot_tg_MAC_meta,dataDate,dataTime,ll)
       extpar_cosmo_buffer(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1:1) = &
       &  MAC_ssa_tg(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,ll:ll,mm)
       CALL write_extpar_cosmo_real_1lev_grib(outfile_id,TRIM(grib_sample), &
       & cosmo_grid,extpar_cosmo_buffer,ssa_tg_MAC_meta,dataDate,dataTime,ll)
       extpar_cosmo_buffer(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1:1) = &
       &  MAC_asy_tg(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,ll:ll,mm)
       CALL write_extpar_cosmo_real_1lev_grib(outfile_id,TRIM(grib_sample), &
       & cosmo_grid,extpar_cosmo_buffer,asy_tg_MAC_meta,dataDate,dataTime,ll)
     ENDDO
   ENDDO
 ELSE
   ! aot_tg
   DO ntype=1,ntype_aot
     SELECT CASE ( ntype )
     CASE ( 1 )
       field_meta = aer_bc_meta
     CASE ( 2 )
       field_meta = aer_dust_meta
     CASE ( 3 )
       field_meta = aer_org_meta
     CASE ( 4 )
       field_meta = aer_so4_meta
     CASE ( 5 )
       field_meta = aer_ss_meta
     END SELECT
     DO mm=1,12
   ! get dataDate and dataTime according to DWD convention for external parameters of climatological monthly mean fields
      CALL set_date_mm_extpar_field(mm,dataDate,dataTime)
      extpar_cosmo_buffer(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1:1) = &
    & aot_tg(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1:1,ntype,mm)
      CALL write_extpar_cosmo_field_grib(outfile_id,TRIM(grib_sample), &
      & cosmo_grid,extpar_cosmo_buffer,field_meta,dataDate,dataTime)
     ENDDO
   ENDDO
 ENDIF

 CALL grib_close_file(outfile_id)



  END SUBROUTINE write_cosmo_grid_extpar_grib
  !-----------------------------------------------------------------------

  !> grib output of external parameters for GME
  SUBROUTINE write_gme_grid_extpar_grib(grib_filename,   &
    &                                     grib_sample,   &
    &                                     gme_grid,      &
    &                                     tg,            &
    &                                     isoil_data,    &
    &                                     ldeep_soil,    &
    &                                     itopo_type,    &
    &                                     lsso,          &
    &                                     undefined,     &
    &                                     undef_int,     &
    &                                     lon_geo,       &
    &                                     lat_geo,       & 
    &                                     fr_land_lu,    & ! &                                     ice_lu, &
    &                                     z0_lu,         &
    &                                     root_lu,       &
    &                                     plcov_mn_lu,   &
    &                                     plcov_mx_lu,   &
    &                                     lai_mn_lu,     &
    &                                     lai_mx_lu,     &
    &                                     rs_min_lu,     &
    &                                     urban_lu,      &
    &                                     for_d_lu,      &
    &                                     for_e_lu,      &
    &                                     emissivity_lu, &
    &                                     lake_depth,    &
    &                                     fr_lake,       &
    &                                     soiltype_fao,  &
    &                                     ndvi_max,      &
    &                                     ndvi_field_mom,&
    &                                     ndvi_ratio_mom,&
    &                                     hh_topo,      &
    &                                     stdh_topo,    &
    &                                     aot_tg,        &
    &                                     crutemp,       &
    &                                     alb_field_mom, &
    &                                     theta_topo,   &
    &                                     aniso_topo,   &
    &                                     slope_topo)
 
  USE mo_var_meta_data, ONLY: lon_geo_meta,           &
    &                         lat_geo_meta,           &
    &                         no_raw_data_pixel_meta, &
    &                         def_com_target_fields_meta  
   
 
  USE mo_gme_grid, ONLY: gme_real_field, gme_int_field, gme_i4_field
  USE mo_gme_grid, ONLY: cp_buf2gme
      
  USE mo_var_meta_data, ONLY: dim_3d_tg, &
    &                         def_dimension_info_buffer

  USE mo_var_meta_data, ONLY: def_lu_fields_meta

  USE mo_var_meta_data, ONLY: dim_lu_tg

  USE mo_var_meta_data, ONLY: fr_land_lu_meta, lu_tot_npixel_meta, &
    &       lu_class_fraction_meta, lu_class_npixel_meta, &
    &       ice_lu_meta, z0_lu_meta, &
    &       plcov_mx_lu_meta, plcov_mn_lu_meta, &
    &       lai_mx_lu_meta, lai_mn_lu_meta, &
    &       rs_min_lu_meta, urban_lu_meta, &
    &       for_d_lu_meta, for_e_lu_meta, &
    &       emissivity_lu_meta, root_lu_meta

  USE mo_var_meta_data, ONLY: def_flake_fields_meta
  USE mo_var_meta_data, ONLY: lake_depth_meta, fr_lake_meta

  USE mo_var_meta_data, ONLY: def_soil_meta
  USE mo_var_meta_data, ONLY: fr_land_soil_meta, soiltype_fao_meta

  USE mo_var_meta_data, ONLY: dim_alb_tg
  USE mo_var_meta_data, ONLY: alb_field_mom_meta, &
      &                       def_alb_meta
  
  USE mo_var_meta_data, ONLY: dim_ndvi_tg
  USE mo_var_meta_data, ONLY: ndvi_max_meta,         &
      &                         ndvi_field_mom_meta, &
      &                         ndvi_ratio_mom_meta, &
      &                         def_ndvi_meta

 USE mo_var_meta_data, ONLY: def_topo_meta, def_topo_vertex_meta
 USE mo_var_meta_data, ONLY: dim_buffer_cell, dim_buffer_vertex

 USE mo_var_meta_data, ONLY: hh_topo_meta, fr_land_topo_meta, &
   &       stdh_topo_meta, theta_topo_meta,                   &
   &       aniso_topo_meta, slope_topo_meta,                  &
   &       hh_vert_meta, npixel_vert_meta,                      &
   &       hh_fis_meta



  USE mo_var_meta_data, ONLY: dim_aot_tg,  &
    &                         aot_tg_meta, &
    &                         def_aot_tg_meta
  USE mo_var_meta_data, ONLY: aot_type_shortname

  USE mo_var_meta_data, ONLY: crutemp_meta, &
    &                         def_crutemp_meta

  USE mo_physical_constants, ONLY: grav


  CHARACTER (len=filename_max), INTENT(IN) :: grib_filename !< filename for the grib file
  CHARACTER (len=filename_max), INTENT(IN) :: grib_sample  !< name for grib sample  (sample to be found in $GRIB_SAMPLES_PATH)

  TYPE(gme_triangular_grid), INTENT(IN) :: gme_grid !< structure which contains the definition of the GME grid
  TYPE(target_grid_def), INTENT(IN) :: tg !< structure with target grid description
  INTEGER (KIND=i4),     INTENT(IN) :: isoil_data
  LOGICAL ,              INTENT(IN) :: ldeep_soil
  INTEGER (KIND=i4),     INTENT(IN) :: itopo_type
  LOGICAL,               INTENT(IN) :: lsso

  REAL(KIND=wp), INTENT(IN)          :: undefined       !< value to indicate undefined grid elements 
  INTEGER, INTENT(IN)                :: undef_int       !< value to indicate undefined grid elements
  REAL (KIND=wp), INTENT(IN) :: lon_geo(:,:,:)  !< longitude coordinates of the target grid in the geographical system
  REAL (KIND=wp), INTENT(IN) :: lat_geo(:,:,:)  !< latitude coordinates of the target grid in the geographical system
  REAL (KIND=wp), INTENT(IN)  :: fr_land_lu(:,:,:) !< fraction land due to lu raw data
  !REAL (KIND=wp), INTENT(IN)  :: ice_lu(:,:,:)     !< fraction of ice due to lu raw data
  REAL (KIND=wp), INTENT(IN)  :: z0_lu(:,:,:)      !< roughness length due to lu land use data
  REAL (KIND=wp), INTENT(IN)  :: root_lu(:,:,:)    !< root depth due to lu land use data
  REAL (KIND=wp), INTENT(IN)  :: plcov_mn_lu(:,:,:)!< plant cover maximum due to lu land use data
  REAL (KIND=wp), INTENT(IN)  :: plcov_mx_lu(:,:,:)!< plant cover minimum due to lu land use data
  REAL (KIND=wp), INTENT(IN)  :: lai_mn_lu(:,:,:)  !< Leaf Area Index maximum due to lu land use data
  REAL (KIND=wp), INTENT(IN)  :: lai_mx_lu(:,:,:)  !< Leaf Area Index minimum due to lu land use data
  REAL (KIND=wp), INTENT(IN)  :: rs_min_lu(:,:,:)  !< minimal stomata resistance due to lu land use data
  REAL (KIND=wp), INTENT(IN)  :: urban_lu(:,:,:)   !< urban fraction due to lu land use data
  REAL (KIND=wp), INTENT(IN)  :: for_d_lu(:,:,:)   !< deciduous forest (fraction) due to lu land use data
  REAL (KIND=wp), INTENT(IN)  :: for_e_lu(:,:,:)   !< evergreen forest (fraction) due to lu land use data
  REAL (KIND=wp), INTENT(IN)  :: emissivity_lu(:,:,:) !< longwave emissivity due to lu land use data

  REAL (KIND=wp), INTENT(IN)  :: lake_depth(:,:,:) !< lake depth
  REAL (KIND=wp), INTENT(IN)  :: fr_lake(:,:,:)     !< fraction of fresh water (lakes)

  INTEGER(KIND=i4), INTENT(IN) :: soiltype_fao(:,:,:) !< soiltype due to FAO Digital Soil map of the World

  REAL (KIND=wp), INTENT(IN) :: alb_field_mom(:,:,:,:) !< field for monthly albedo data

  REAL (KIND=wp), INTENT(IN) :: ndvi_max(:,:,:) !< field for ndvi maximum
  REAL (KIND=wp), INTENT(IN) :: ndvi_field_mom(:,:,:,:) !< field for monthly mean ndvi data (12 months)
  REAL (KIND=wp), INTENT(IN) :: ndvi_ratio_mom(:,:,:,:) !< field for monthly ndvi ratio (12 months)

  REAL(KIND=wp), INTENT(IN)  :: hh_topo(:,:,:)  !< mean height 
  REAL(KIND=wp), INTENT(IN)  :: stdh_topo(:,:,:) !< standard deviation of subgrid scale orographic height
  REAL (KIND=wp), INTENT(IN)  :: aot_tg(:,:,:,:,:) !< aerosol optical thickness, aot_tg(ie,je,ke,,ntype,ntime)
  REAL(KIND=wp), INTENT(IN)  :: crutemp(:,:,:)  !< cru climatological temperature , crutemp(ie,je,ke) 
  REAL(KIND=wp), INTENT(IN), OPTIONAL  :: theta_topo(:,:,:) !< sso parameter, angle of principal axis
  REAL(KIND=wp), INTENT(IN), OPTIONAL  :: aniso_topo(:,:,:) !< sso parameter, anisotropie factor
  REAL(KIND=wp), INTENT(IN), OPTIONAL  :: slope_topo(:,:,:) !< sso parameter, mean slope


  ! local variables

  INTEGER :: outfile_id
  INTEGER :: errorcode

  INTEGER (KIND=i8)  :: dataDate
  INTEGER (KIND=i8)  :: dataTime
  INTEGER :: mm ! month
  INTEGER :: ntype ! type of aerosol
  TYPE(var_meta_info) :: field_meta

  REAL (KIND=wp) :: real_buffer(1:gme_grid%nip1,1:gme_grid%nip1,1:gme_grid%nd)

  INTEGER :: i,j,k,t



  ! prepar meta data for output
  !set up dimensions for buffer
  CALL  def_dimension_info_buffer(tg)
  ! dim_3d_tg
  
  ! define meta information for various land use related variables (GLC2000) for netcdf output
  CALL def_lu_fields_meta(tg,nclass_glc2000,dim_3d_tg)
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
  CALL def_com_target_fields_meta(dim_3d_tg)
  ! lon_geo_meta and lat_geo_meta

  CALL def_soil_meta(dim_3d_tg,isoil_data)
  !  fr_land_soil_meta, soiltype_fao_meta

  ! define meta information for various albedo variables for nc output
  CALL def_alb_meta(tg,ntime_alb,dim_3d_tg)

  !define meta information for various NDVI data related variables for netcdf output
  CALL def_ndvi_meta(tg,ntime_ndvi,dim_3d_tg)
  ! dim_ndvi_tg, ndvi_max_meta, ndvi_field_mom_meta, ndvi_ratio_mom_meta

  ! define meta information for various TOPO data related variables for netcdf output
  CALL def_topo_meta(dim_3d_tg,itopo_type)
  !  hh_topo_meta, fr_land_topo_meta, &
  !         stdh_topo_meta, theta_topo_meta, &
  !         aniso_topo_meta, slope_topo_meta, &
  !         hh_vert_meta, npixel_vert_meta

  ! define dimensions and meta information for variable aot_tg for netcdf output
  CALL def_aot_tg_meta(tg,ntime_aot,ntype_aot,dim_3d_tg)
  ! dim_aot_tg and aot_tg_meta

  ! define meta information for variable crutemp for netcdf output
  CALL def_crutemp_meta(dim_3d_tg)
  ! crutemp_meta

   CALL def_flake_fields_meta(dim_3d_tg)
   ! lake_depth_meta, fr_lake_meta



  ! open a new grib file
  CALL grib_open_file(outfile_id,TRIM(grib_filename),'w')


  ! get dateDate according to DWD convention for invariant external parameter fields
  CALL get_date_const_field(dataDate,dataTime)

 
  ! write out fields to grib file
  PRINT *,'fr_land_lu_meta%shortName: ',fr_land_lu_meta%shortName

  CALL write_extpar_gme_field_grib(outfile_id,TRIM(grib_sample),&
  & gme_grid,fr_land_lu,fr_land_lu_meta,dataDate,dataTime)

  PRINT *,'hh_topo_meta%shortName: ',hh_topo_meta%shortName

  CALL write_extpar_gme_field_grib(outfile_id,TRIM(grib_sample),&
  & gme_grid,hh_topo,hh_topo_meta,dataDate,dataTime)
  ! output FIS as well
  PRINT *,'hh_fis_meta%shortName: ',hh_fis_meta%shortName
  CALL write_extpar_gme_field_grib(outfile_id,TRIM(grib_sample),&
  & gme_grid,hh_topo*grav,hh_fis_meta,dataDate,dataTime)

  CALL write_extpar_gme_field_grib(outfile_id,TRIM(grib_sample),&
  & gme_grid,stdh_topo,stdh_topo_meta,dataDate,dataTime)
  
  IF (lsso) THEN
  CALL write_extpar_gme_field_grib(outfile_id,TRIM(grib_sample),&
  & gme_grid,theta_topo,theta_topo_meta,dataDate,dataTime)
  ENDIF

  IF (lsso) THEN
  CALL write_extpar_gme_field_grib(outfile_id,TRIM(grib_sample),&
  & gme_grid,aniso_topo,aniso_topo_meta,dataDate,dataTime)
  ENDIF

  IF (lsso) THEN
  CALL write_extpar_gme_field_grib(outfile_id,TRIM(grib_sample), &
  & gme_grid,slope_topo,slope_topo_meta,dataDate,dataTime)
  ENDIF

  CALL write_extpar_gme_field_grib(outfile_id,TRIM(grib_sample),&
  & gme_grid,plcov_mn_lu,plcov_mn_lu_meta,dataDate,dataTime)
  
  CALL write_extpar_gme_field_grib(outfile_id,TRIM(grib_sample),&
  & gme_grid,plcov_mx_lu,plcov_mx_lu_meta,dataDate,dataTime)
  
  CALL write_extpar_gme_field_grib(outfile_id,TRIM(grib_sample), &
  & gme_grid,lai_mn_lu,lai_mn_lu_meta,dataDate,dataTime)
  
  CALL write_extpar_gme_field_grib(outfile_id,TRIM(grib_sample),&
  & gme_grid,lai_mx_lu,lai_mx_lu_meta,dataDate,dataTime)

  CALL write_extpar_gme_field_grib(outfile_id,TRIM(grib_sample), &
  & gme_grid,rs_min_lu,rs_min_lu_meta,dataDate,dataTime)
  
  CALL write_extpar_gme_field_grib(outfile_id,TRIM(grib_sample), &
  & gme_grid,for_e_lu,for_e_lu_meta,dataDate,dataTime)
  
  CALL write_extpar_gme_field_grib(outfile_id,TRIM(grib_sample), &
  & gme_grid,for_d_lu,for_d_lu_meta,dataDate,dataTime)
  
  CALL write_extpar_gme_field_grib(outfile_id,TRIM(grib_sample), &
  & gme_grid,emissivity_lu,emissivity_lu_meta,dataDate,dataTime)
  
  CALL write_extpar_gme_field_grib(outfile_id,TRIM(grib_sample), &
  & gme_grid,root_lu,root_lu_meta,dataDate,dataTime)

  CALL write_extpar_gme_field_grib(outfile_id,TRIM(grib_sample), &
  & gme_grid,z0_lu,z0_lu_meta,dataDate,dataTime)

  CALL write_extpar_gme_field_grib(outfile_id,TRIM(grib_sample),&
  & gme_grid,lake_depth,lake_depth_meta,dataDate,dataTime)
  
  CALL write_extpar_gme_field_grib(outfile_id,TRIM(grib_sample),&
  & gme_grid,fr_lake,fr_lake_meta,dataDate,dataTime)

  CALL write_extpar_gme_field_grib(outfile_id,TRIM(grib_sample),&
  & gme_grid,soiltype_fao,soiltype_fao_meta,dataDate,dataTime)

  CALL write_extpar_gme_field_grib(outfile_id,TRIM(grib_sample),&
  & gme_grid,crutemp,crutemp_meta,dataDate,dataTime)
  
  CALL write_extpar_gme_field_grib(outfile_id,TRIM(grib_sample), &
  & gme_grid,lon_geo,lon_geo_meta,dataDate,dataTime)
  
  CALL write_extpar_gme_field_grib(outfile_id,TRIM(grib_sample),&
  & gme_grid,lat_geo,lat_geo_meta,dataDate,dataTime)


 ! urban_lu
 ! ndvi_max
 CALL write_extpar_gme_field_grib(outfile_id,TRIM(grib_sample), &
 & gme_grid,ndvi_max,ndvi_max_meta,dataDate,dataTime)
 ! ndvi_field_mom
 
 DO mm=1,12
    ! get dataDate and dataTime according to DWD convention for external parameters of climatological monthly mean fields
    CALL set_date_mm_extpar_field(mm,dataDate,dataTime)
    real_buffer(1:gme_grid%nip1,1:gme_grid%nip1,1:gme_grid%nd) = &
  & ndvi_field_mom(1:gme_grid%nip1,1:gme_grid%nip1,1:gme_grid%nd,mm)
    CALL write_extpar_gme_field_grib(outfile_id,TRIM(grib_sample),&
    & gme_grid,real_buffer,ndvi_field_mom_meta,dataDate,dataTime)
 ENDDO
 ! ndvi_ratio_mom
 DO mm=1,12
    ! get dataDate and dataTime according to DWD convention for external parameters of climatological monthly mean fields
    CALL set_date_mm_extpar_field(mm,dataDate,dataTime)
    real_buffer(1:gme_grid%nip1,1:gme_grid%nip1,1:gme_grid%nd) = &
  & ndvi_ratio_mom(1:gme_grid%nip1,1:gme_grid%nip1,1:gme_grid%nd,mm) 
    CALL write_extpar_gme_field_grib(outfile_id,TRIM(grib_sample),&
    & gme_grid,real_buffer,ndvi_ratio_mom_meta,dataDate,dataTime)
 ENDDO

 ! alb_field_mom
 DO mm=1,12
    CALL set_date_mm_extpar_field(mm,dataDate,dataTime)
    real_buffer(1:gme_grid%nip1,1:gme_grid%nip1,1:gme_grid%nd) = &
  & alb_field_mom(1:gme_grid%nip1,1:gme_grid%nip1,1:gme_grid%nd,mm)
    CALL write_extpar_gme_field_grib(outfile_id,TRIM(grib_sample),&
    & gme_grid,real_buffer,alb_field_mom_meta,dataDate,dataTime)
 ENDDO

 ! aot_tg

 DO ntype=1,ntype_aot
   SELECT CASE ( ntype )
   CASE ( 1 )
     field_meta = aer_bc_meta
   CASE ( 2 )
     field_meta = aer_dust_meta
   CASE ( 3 )
     field_meta = aer_org_meta
   CASE ( 4 )
     field_meta = aer_so4_meta
   CASE ( 5 )
     field_meta = aer_ss_meta
   END SELECT
   DO mm=1,12
   ! get dataDate and dataTime according to DWD convention for external parameters of climatological monthly mean fields
    CALL set_date_mm_extpar_field(mm,dataDate,dataTime)
    real_buffer(1:gme_grid%nip1,1:gme_grid%nip1,1:gme_grid%nd) = &
  & aot_tg(1:gme_grid%nip1,1:gme_grid%nip1,1:gme_grid%nd,ntype,mm)
    CALL write_extpar_gme_field_grib(outfile_id,TRIM(grib_sample), &
    & gme_grid,real_buffer,field_meta,dataDate,dataTime)
   ENDDO
 ENDDO
 ! ice_lu


  CALL grib_close_file(outfile_id)



  END SUBROUTINE write_gme_grid_extpar_grib
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !> grib output of external parameters for GME
  SUBROUTINE write_ICON_grid_extpar_grib(grib_filename,  &
    &                                     grib_sample,       &
    &                                     ICON_grid,       &
    &                                     tg,         &
    &                                     isoil_data,    &
    &                                     ldeep_soil,    &
    &                                     itopo_type,    &
    &                                     lsso,          &
    &                                     undefined, &
    &                                     undef_int,   &
    &                                     lon_geo,     &
    &                                     lat_geo, & 
    &                                     fr_land_lu, & 
    &                                     lu_class_fraction,    &
    &                                     ice_lu, &
    &                                     z0_lu, &
    &                                     root_lu, &
    &                                     plcov_mn_lu, &
    &                                     plcov_mx_lu, &
    &                                     lai_mn_lu, &
    &                                     lai_mx_lu, &
    &                                     rs_min_lu, &
    &                                     urban_lu,  &
    &                                     for_d_lu,  &
    &                                     for_e_lu, &
    &                                     emissivity_lu, &
    &                                     lake_depth, &
    &                                     fr_lake,    &
    &                                     soiltype_fao, &
    &                                     ndvi_max,  &
    &                                     ndvi_field_mom,&
    &                                     ndvi_ratio_mom, &
    &                                     hh_topo,            &
    &                                     stdh_topo,          &
    &                                     theta_topo,         &
    &                                     aniso_topo,         &
    &                                     slope_topo,   &
    &                                     vertex_param,        &
    &                                     aot_tg, &
    &                                     crutemp, &
    &                                     alb_field_mom, &
    &                                     alnid_field_mom, &
    &                                     aluvd_field_mom)
 
  USE mo_var_meta_data, ONLY: clon_meta, &
    &                         clat_meta, &
    &                         no_raw_data_pixel_meta, &
    &                         def_com_target_fields_meta  

  USE mo_var_meta_data, ONLY:  dim_icon, &
    &                          def_dimension_info_icon
   
   USE mo_var_meta_data, ONLY: nc_grid_def_icon, &
    &                         set_nc_grid_def_icon

     
  USE mo_var_meta_data, ONLY: dim_3d_tg, &
    &                         def_dimension_info_buffer

  USE mo_var_meta_data, ONLY: def_lu_fields_meta

  USE mo_var_meta_data, ONLY: dim_lu_tg

  USE mo_var_meta_data, ONLY: fr_land_lu_meta, lu_tot_npixel_meta, &
    &       lu_class_fraction_meta, lu_class_npixel_meta, &
    &       ice_lu_meta, z0_lu_meta, &
    &       plcov_mx_lu_meta, plcov_mn_lu_meta, &
    &       lai_mx_lu_meta, lai_mn_lu_meta, &
    &       rs_min_lu_meta, urban_lu_meta, &
    &       for_d_lu_meta, for_e_lu_meta, &
    &       emissivity_lu_meta, root_lu_meta

  USE mo_var_meta_data, ONLY: def_flake_fields_meta
  USE mo_var_meta_data, ONLY: lake_depth_meta, fr_lake_meta

  USE mo_var_meta_data, ONLY: def_soil_meta
  USE mo_var_meta_data, ONLY: fr_land_soil_meta, soiltype_fao_meta

  USE mo_var_meta_data, ONLY: dim_alb_tg
  USE mo_var_meta_data, ONLY: alb_field_mom_meta, alnid_field_mom_meta, aluvd_field_mom_meta, &
      &                       def_alb_meta


  
  USE mo_var_meta_data, ONLY: dim_ndvi_tg
  USE mo_var_meta_data, ONLY: ndvi_max_meta, &
      &                         ndvi_field_mom_meta, &
      &                         ndvi_ratio_mom_meta,&
      &                         def_ndvi_meta

  USE mo_var_meta_data, ONLY: def_topo_meta, def_topo_vertex_meta
  USE mo_var_meta_data, ONLY: dim_buffer_cell, dim_buffer_vertex

  USE mo_var_meta_data, ONLY: hh_topo_meta, fr_land_topo_meta, &
   &       stdh_topo_meta, theta_topo_meta, &
   &       aniso_topo_meta, slope_topo_meta, &
   &       hh_vert_meta, npixel_vert_meta,    &
   &       hh_fis_meta

  USE mo_var_meta_data, ONLY: def_topo_meta, def_topo_vertex_meta
  USE mo_var_meta_data, ONLY: dim_buffer_cell, dim_buffer_vertex


  USE mo_var_meta_data, ONLY: dim_aot_tg, &
    &                         aot_tg_meta, &
    &                         def_aot_tg_meta
  USE mo_var_meta_data, ONLY: aot_type_shortname,aer_bc_meta,aer_dust_meta,&
                                   aer_org_meta,aer_so4_meta,aer_ss_meta

  USE mo_var_meta_data, ONLY: crutemp_meta, &
    &                         def_crutemp_meta

  USE mo_physical_constants, ONLY: grav
 
  USE mo_icon_grid_data, ONLY: icon_grid_region

  USE mo_icon_grid_data, ONLY: clon, clat
  USE mo_icon_grid_data, ONLY: clon_vertices, clat_vertices
  USE mo_icon_grid_data, ONLY: allocate_icon_coor

  USE mo_topo_tg_fields, ONLY: add_parameters_domain

  CHARACTER (len=*), INTENT(IN) :: grib_filename !< filename for the grib file
  CHARACTER (len=*), INTENT(IN) :: grib_sample  !< name for grib sample  (sample to be found in $GRIB_SAMPLES_PATH)
  
  CHARACTER datum*8,zeit*10


  TYPE(icosahedral_triangular_grid), INTENT(IN) :: icon_grid !< structure which contains the definition of the ICON grid
  TYPE(target_grid_def), INTENT(IN) :: tg !< structure with target grid description

  INTEGER (KIND=i4),     INTENT(IN) :: isoil_data
  LOGICAL ,              INTENT(IN) :: ldeep_soil
  INTEGER (KIND=i4),     INTENT(IN) :: itopo_type
  LOGICAL,               INTENT(IN) :: lsso

  REAL(KIND=wp), INTENT(IN)          :: undefined       !< value to indicate undefined grid elements 
  INTEGER, INTENT(IN)                :: undef_int       !< value to indicate undefined grid elements
  REAL (KIND=wp), INTENT(IN) :: lon_geo(:,:,:)  !< longitude coordinates of the target grid in the geographical system
  REAL (KIND=wp), INTENT(IN) :: lat_geo(:,:,:)  !< latitude coordinates of the target grid in the geographical system
  REAL (KIND=wp), INTENT(IN)  :: fr_land_lu(:,:,:) !< fraction land due to lu raw data
  REAL (KIND=wp), INTENT(IN)  :: lu_class_fraction(:,:,:,:)  !< fraction land use class according to look-up table
  REAL (KIND=wp), INTENT(IN)  :: ice_lu(:,:,:)     !< fraction of ice due to lu raw data
  REAL (KIND=wp), INTENT(IN)  :: z0_lu(:,:,:)      !< roughness length due to lu land use data
  REAL (KIND=wp), INTENT(IN)  :: root_lu(:,:,:)    !< root depth due to lu land use data
  REAL (KIND=wp), INTENT(IN)  :: plcov_mn_lu(:,:,:)!< plant cover maximum due to lu land use data
  REAL (KIND=wp), INTENT(IN)  :: plcov_mx_lu(:,:,:)!< plant cover minimum due to lu land use data
  REAL (KIND=wp), INTENT(IN)  :: lai_mn_lu(:,:,:)  !< Leaf Area Index maximum due to lu land use data
  REAL (KIND=wp), INTENT(IN)  :: lai_mx_lu(:,:,:)  !< Leaf Area Index minimum due to lu land use data
  REAL (KIND=wp), INTENT(IN)  :: rs_min_lu(:,:,:)  !< minimal stomata resistance due to lu land use data
  REAL (KIND=wp), INTENT(IN)  :: urban_lu(:,:,:)   !< urban fraction due to lu land use data
  REAL (KIND=wp), INTENT(IN)  :: for_d_lu(:,:,:)   !< deciduous forest (fraction) due to lu land use data
  REAL (KIND=wp), INTENT(IN)  :: for_e_lu(:,:,:)   !< evergreen forest (fraction) due to lu land use data
  REAL (KIND=wp), INTENT(IN)  :: emissivity_lu(:,:,:) !< longwave emissivity due to lu land use data

  REAL (KIND=wp), INTENT(IN)  :: lake_depth(:,:,:) !< lake depth
  REAL (KIND=wp), INTENT(IN)  :: fr_lake(:,:,:)     !< fraction of fresh water (lakes)

  INTEGER(KIND=i4), INTENT(IN) :: soiltype_fao(:,:,:) !< soiltype due to FAO Digital Soil map of the World

  REAL (KIND=wp), INTENT(IN) :: alb_field_mom(:,:,:,:) !< field for monthly albedo data
  REAL (KIND=wp), INTENT(IN) :: alnid_field_mom(:,:,:,:) !< field for monthly albedo data
  REAL (KIND=wp), INTENT(IN) :: aluvd_field_mom(:,:,:,:) !< field for monthly albedo data

  REAL (KIND=wp), INTENT(IN) :: ndvi_max(:,:,:) !< field for ndvi maximum
  REAL (KIND=wp), INTENT(IN) :: ndvi_field_mom(:,:,:,:) !< field for monthly mean ndvi data (12 months)
  REAL (KIND=wp), INTENT(IN) :: ndvi_ratio_mom(:,:,:,:) !< field for monthly ndvi ratio (12 months)

  REAL(KIND=wp), INTENT(IN)  :: hh_topo(:,:,:)  !< mean height 
  REAL(KIND=wp), INTENT(IN)  :: stdh_topo(:,:,:) !< standard deviation of subgrid scale orographic height
  REAL(KIND=wp), INTENT(IN)  :: theta_topo(:,:,:) !< sso parameter, angle of principal axis
  REAL(KIND=wp), INTENT(IN)  :: aniso_topo(:,:,:) !< sso parameter, anisotropie factor
  REAL(KIND=wp), INTENT(IN)  :: slope_topo(:,:,:) !< sso parameter, mean slope
  TYPE(add_parameters_domain), INTENT(IN) :: vertex_param !< additional external parameters for ICON domain
  REAL(KIND=wp), INTENT(IN)  :: aot_tg(:,:,:,:,:) !< aerosol optical thickness, aot_tg(ie,je,ke,,ntype,ntime)
  REAL(KIND=wp), INTENT(IN)  :: crutemp(:,:,:)  !< cru climatological temperature , crutemp(ie,je,ke) 

  ! local variables

  INTEGER :: outfile_id
  INTEGER :: gribid_in !< id of grib message (GRIB_API)
  INTEGER :: gribid_dest !< id of grib message (GRIB_API)
  INTEGER :: igribclone,ogrib,igrib,infile, outfile
  INTEGER :: errorcode

  INTEGER (KIND=i8)  :: dataDate,dataDate_default
  INTEGER (KIND=i8)  :: dataTime,dataTime_default
  INTEGER :: mm ! month
  INTEGER :: ntype ! type of aerosol
  INTEGER,PARAMETER :: ntype_aot=5 ! type of aerosol

  REAL (KIND=wp) :: real_buffer(1:icon_grid%ncell,1,1)

  INTEGER :: i,j,k,t
  INTEGER :: d_run,d_p,d_e,d_li
  INTEGER :: lcY,lcM,lcD,lcH,lcMM,lcS
  INTEGER :: list(23) ! List of land-use class code numbers
  CHARACTER(LEN=1 ) :: uuid(16)    !   UUID of unstructured grids 

  CALL def_topo_meta(dim_3d_tg,itopo_type)
    !  hh_topo_meta, fr_land_topo_meta, &
    !         stdh_topo_meta, theta_topo_meta, &
    !         aniso_topo_meta, slope_topo_meta, &
    !         hh_vert_meta, npixel_vert_meta


  ! open a new grib file
  CALL grib_open_file(outfile_id,TRIM(grib_filename),'w')
  CALL grib_new_from_samples(gribid_dest,"GRIB2")

! set defaults from grib2 template
    !----------
    ! Section 1
    !----------
    ! grib1/0.table


  WRITE(0,*) "GRIB2-File ",TRIM(grib_filename)," opend. Setting first default values"

  dataDate = 10101
  dataTime = 0

  CALL grib_set (gribid_dest,'centre'                     ,78) ! DWD
  CALL grib_set (gribid_dest,'subCentre'                  ,255)! DWD

!Only needed if a wrong value in the Template exist: was the case for early grib_api-versions
!   call grib_set (gribid_dest,'tablesVersion',11)
!  CALL grib_set (gribid_dest,'tablesVersion',5)


  ! grib2/tables/4/1.2.table :: 0=ana,1=fc_start,2=veri_time
  CALL grib_set (gribid_dest,'significanceOfReferenceTime', 0)
  CALL grib_set (gribid_dest,'dataDate',dataDate)
  CALL grib_set (gribid_dest,'dataTime',dataTime)

  ! grib2/tables/5/1.3.table   :: 0=oper,1=test,2=exp
  CALL grib_set (gribid_dest,'productionStatusOfProcessedData', 0)
  ! grib2/tables/5/1.4.table   :: 0=ana,1=fc,...
  CALL grib_set (gribid_dest,'typeOfProcessedData', 0)
  CALL grib_set (gribid_dest,'stepType','instant')
  CALL grib_set (gribid_dest,'startStep',0)
  CALL grib_set (gribid_dest,'endStep',0)
  CALL grib_set (gribid_dest,'stepUnits','M')

  !---------------------
  ! Section 2: local use
  !---------------------

  CALL date_and_time(datum,zeit)

  WRITE (*,*) 'GRIB2 local creation date and time is ',datum,zeit

  READ(datum(1:4),'(I4)')  lcY
  READ(datum(5:6),'(I2)')  lcM
  READ(datum(7:8),'(I2)')  lcD

  READ(zeit(1:2),'(I2)')  lcH
  READ(zeit(3:4),'(I2)')  lcMM
  READ(zeit(5:6),'(I2)')  lcS



  CALL grib_set (gribid_dest,'grib2LocalSectionPresent', 1)
  CALL grib_set (gribid_dest,'localDefinitionNumber',254) !254=deterministic
!  CALL grib_set_missing (gribid_dest,'localHostIdentifier')
  CALL grib_set (gribid_dest,'localHostIdentifier'    ,255) ! 1=oper.,2=backup,???
  CALL grib_set (gribid_dest,'localCreationDateYear'  ,lcY)
  CALL grib_set (gribid_dest,'localCreationDateMonth' ,lcM)
  CALL grib_set (gribid_dest,'localCreationDateDay'   ,lcD)
  CALL grib_set (gribid_dest,'localCreationDateHour'  ,lcH)
  CALL grib_set (gribid_dest,'localCreationDateMinute',lcMM)
  CALL grib_set (gribid_dest,'localCreationDateSecond',lcS)

  CALL grib_set (gribid_dest,'localNumberOfExperiment',255)
  CALL grib_set (gribid_dest,'localInformationNumber' ,255)  ! ZEN

  !------------------------------
  ! Section 3 (GDS)
  ! template.3.*: grid definition
  !------------------------------
!  CALL grib_set (gribid_dest,'numberOfDataPoints',0)   ! will be set later
  ! grib2/tables/4/3.11.table
  CALL grib_set (gribid_dest,'interpretationOfNumberOfPoints',0)
  ! grib2/tables/4/3.1.table
!  CALL grib_set (gribid_dest,'gridType','unstructured_grid')    ! not yet implemented
  CALL grib_set (gribid_dest,'gridDefinitionTemplateNumber',101)
  ! grib2/tables/4/3.2.table :: R = 6,371,229.0 m
  CALL grib_set (gribid_dest,'shapeOfTheEarth'             ,6)
  CALL grib_set (gribid_dest,'numberOfGridUsed'            ,icon_grid%number_Of_Grid_Used) ! Temporary value
  ! grib2/tables/4/3.8.table
  CALL grib_set (gribid_dest,'numberOfGridInReference'     ,1)  ! Triangle centers
!!$    call grib_set (gribid_dest,'numberOfGridInReference'     ,2)  ! Triangle edge
!!$    call grib_set (gribid_dest,'numberOfGridInReference'     ,3)  ! Triangle corner
!!! This requires GRIB_API version >= 1.10.4
  CALL decode_uuid (icon_grid%uuidOfHGrid, uuid) 
  CALL grib_set (gribid_dest,'uuidOfHGrid'                 ,uuid)
  CALL grib_set (gribid_dest,'numberOfDataPoints'          ,icon_grid%ncell)

  !----------
  ! Section 4 (mostly)
  !----------
  ! template.4
  ! grib2/tables/4/4.0.table :: 0=Deterministic, 1=Ensemble, 7=fg/anaerr...
  CALL grib_set (gribid_dest,'productDefinitionTemplateNumber', 0)
!  CALL grib_set (gribid_dest,'productDefinitionTemplateNumber', 1)

  ! grib2/tables/4/0.0.table
  CALL grib_set (gribid_dest,'discipline'       ,0)
  ! grib2/tables/4/4.1.*.table
  CALL grib_set (gribid_dest,'parameterCategory',0)
  ! grib2/tables/4/4.2.*.*.table
  CALL grib_set (gribid_dest,'parameterNumber'  ,0)

  ! grib2/tables/4/4.3.table :: 0=ana,1=ini,2=fc,7=anaerr, 9=clim.,...
  ! Initialize using known DWD-defined generating process identifiers:
  CALL grib_set (gribid_dest,'typeOfGeneratingProcess', 196) ! for all invar data
 

  ! grib2/backgroundProcess.table :: 0=main,1=pre,2=ana,3=test
  CALL grib_set (gribid_dest,'backgroundProcess', 0)
    ! call grib_set (gribid_dest,'backgroundGeneratingProcessIdentifier', 0)
  ! grib2/generatingProcessIdentifier.table
  CALL grib_set (gribid_dest,'generatingProcessIdentifier', 1)

  ! grib2/tables/4/4.4.table
  CALL grib_set (gribid_dest,'indicatorOfUnitOfTimeRange',1)
  CALL grib_set (gribid_dest,'stepUnits',1)

  CALL grib_set (gribid_dest,'numberOfValues', 0)
  !---------- 
  ! Section 5
  !----------
  CALL grib_set (gribid_dest,'packingType' ,'grid_simple')
  !----------
  ! Section 7
  !----------
  CALL grib_set (gribid_dest,'bitsPerValue', 16 )!grib% isec4% bits)
! end defaults
!!$
  CALL grib_set (gribid_dest,'typeOfFirstFixedSurface',1) ![Ground or water surface  (grib2/tables/5/4.5.table) ]
  CALL grib_set (gribid_dest,'typeOfSecondFixedSurface',255) ![Missing (grib2/tables/5/4.5.table) ]
!!$    CALL grib_set (gribid_dest,'scaleFactorOfSecondFixedSurface','MISSING')
!!$    CALL grib_set (gribid_dest,'scaledValueOfSecondFixedSurface','MISSING')

  WRITE (0,*) 'Now writing EXTPAR fields into ',TRIM(grib_filename)



  WRITE (0,*) fr_land_lu_meta%shortName," -> ",TRIM(grib_filename)
  CALL write_extpar_ICON_field_grib(outfile_id,gribid_dest,ICON_grid,&
         fr_land_lu,fr_land_lu_meta%shortName,dataDate,dataTime) 

  CALL grib_set (gribid_dest,'typeOfSecondFixedSurface',101) 

  WRITE (0,*) hh_topo_meta%shortName," -> ",TRIM(grib_filename)
  CALL write_extpar_ICON_field_grib(outfile_id,gribid_dest,ICON_grid,&
         hh_topo,hh_topo_meta%shortName,dataDate,dataTime)


  CALL grib_set (gribid_dest,'typeOfSecondFixedSurface',255) ![Missing (grib2/tables/5/4.5.table) ]
! FIS not needed
!!$ write(0,*) hh_fis_meta%shortName," -> ",TRIM(grib_filename)
!!$  CALL write_extpar_ICON_field_grib(outfile_id,gribid_dest,ICON_grid,&
!!$         hh_topo*grav,hh_fis_meta%shortName,dataDate,dataTime)

  WRITE (0,*) stdh_topo_meta%shortName," -> ",TRIM(grib_filename)
  CALL write_extpar_ICON_field_grib(outfile_id,gribid_dest,ICON_grid,&
         stdh_topo,stdh_topo_meta%shortName,dataDate,dataTime)
 
  WRITE (0,*) theta_topo_meta%shortName," -> ",TRIM(grib_filename)
  CALL write_extpar_ICON_field_grib(outfile_id,gribid_dest,ICON_grid,&
         theta_topo,theta_topo_meta%shortName,dataDate,dataTime)

  WRITE (0,*) aniso_topo_meta%shortName," -> ",TRIM(grib_filename)
  CALL write_extpar_ICON_field_grib(outfile_id,gribid_dest,ICON_grid,&
         aniso_topo,aniso_topo_meta%shortName,dataDate,dataTime)

  WRITE (0,*) slope_topo_meta%shortName," -> ",TRIM(grib_filename)
  CALL write_extpar_ICON_field_grib(outfile_id,gribid_dest,ICON_grid,&
         slope_topo,slope_topo_meta%shortName,dataDate,dataTime)

  WRITE (0,*) for_e_lu_meta%shortName," -> ",TRIM(grib_filename)
  CALL write_extpar_ICON_field_grib(outfile_id,gribid_dest,ICON_grid,&
         for_e_lu,for_e_lu_meta%shortName,dataDate,dataTime)

  WRITE (0,*) for_d_lu_meta%shortName," -> ",TRIM(grib_filename)
  CALL write_extpar_ICON_field_grib(outfile_id,gribid_dest,ICON_grid,&
           for_d_lu,for_d_lu_meta%shortName,dataDate,dataTime)

  WRITE (0,*) emissivity_lu_meta%shortName," -> ",TRIM(grib_filename)
  CALL write_extpar_ICON_field_grib(outfile_id,gribid_dest,ICON_grid,&
           emissivity_lu,emissivity_lu_meta%shortName,dataDate,dataTime)

  WRITE (0,*) root_lu_meta%shortName," -> ",TRIM(grib_filename)
  CALL write_extpar_ICON_field_grib(outfile_id,gribid_dest,ICON_grid,&
           root_lu,root_lu_meta%shortName,dataDate,dataTime)

  WRITE (0,*) z0_lu_meta%shortName," -> ",TRIM(grib_filename)
  CALL write_extpar_ICON_field_grib(outfile_id,gribid_dest,ICON_grid,&
           z0_lu,z0_lu_meta%shortName,dataDate,dataTime)

  WRITE (0,*) lake_depth_meta%shortName," -> ",TRIM(grib_filename)
  CALL write_extpar_ICON_field_grib(outfile_id,gribid_dest,ICON_grid,&
          lake_depth,lake_depth_meta%shortName,dataDate,dataTime)

  WRITE (0,*) fr_lake_meta%shortName," -> ",TRIM(grib_filename)
  CALL write_extpar_ICON_field_grib(outfile_id,gribid_dest,ICON_grid,&
          fr_lake,fr_lake_meta%shortName,dataDate,dataTime)
 
  WRITE (0,*) soiltype_fao_meta%shortName," -> ",TRIM(grib_filename)
  CALL write_extpar_ICON_field_grib(outfile_id,gribid_dest,ICON_grid,&
           soiltype_fao,soiltype_fao_meta%shortName,dataDate,dataTime)

  WRITE (0,*) clon_meta%shortName," -> ",TRIM(grib_filename)
  CALL write_extpar_ICON_field_grib(outfile_id,gribid_dest,ICON_grid,&
           lon_geo,clon_meta%shortName,dataDate,dataTime)
  WRITE (0,*) clat_meta%shortName," -> ",TRIM(grib_filename)
  CALL write_extpar_ICON_field_grib(outfile_id,gribid_dest,ICON_grid,&
            lat_geo,clat_meta%shortName,dataDate,dataTime)

! FR_ICE now stored in FR_LUC with Antarctica
!!$ write(0,*) ice_lu_meta%shortName," -> ",TRIM(grib_filename)
!!$  CALL write_extpar_ICON_field_grib(outfile_id,gribid_dest,ICON_grid,&
!!$            ice_lu,ice_lu_meta%shortName,dataDate,dataTime)

  CALL grib_set (gribid_dest,'stepType','max')
!    call grib_set (gribid_dest,'startStep',0)
!    call grib_set (gribid_dest,'endStep',12)
  WRITE (0,*) ndvi_max_meta%shortName," -> ",TRIM(grib_filename)
  CALL write_extpar_ICON_field_grib(outfile_id,gribid_dest,ICON_grid,&
             ndvi_max,ndvi_max_meta%shortName,dataDate,dataTime)

  WRITE (0,*) plcov_mx_lu_meta%shortName," -> ",TRIM(grib_filename)
  CALL write_extpar_ICON_field_grib(outfile_id,gribid_dest,ICON_grid,&
         plcov_mx_lu,plcov_mx_lu_meta%shortName,dataDate,dataTime)
 
  WRITE (0,*) lai_mx_lu_meta%shortName," -> ",TRIM(grib_filename)
  CALL write_extpar_ICON_field_grib(outfile_id,gribid_dest,ICON_grid,&
         lai_mx_lu,lai_mx_lu_meta%shortName,dataDate,dataTime)

  CALL grib_set (gribid_dest,'stepType','min')
!    call grib_set (gribid_dest,'startStep',0)
!    call grib_set (gribid_dest,'endStep',12)
  WRITE (0,*) rs_min_lu_meta%shortName," -> ",TRIM(grib_filename)
  CALL write_extpar_ICON_field_grib(outfile_id,gribid_dest,ICON_grid,&
         rs_min_lu,rs_min_lu_meta%shortName,dataDate,dataTime)

  CALL grib_set (gribid_dest,'stepType','avg')
!    call grib_set (gribid_dest,'startStep',0)
!    call grib_set (gribid_dest,'endStep',12)
  WRITE (0,*) crutemp_meta%shortName," -> ",TRIM(grib_filename)
  CALL write_extpar_ICON_field_grib(outfile_id,gribid_dest,ICON_grid,&
           crutemp,crutemp_meta%shortName,dataDate,dataTime)


!--> Climatology Loop
!
  CALL grib_set (gribid_dest,'typeOfGeneratingProcess', 9) ! for climate data only

  dataDate_default=11110111


  DO mm=1,12 
  
    dataDate=dataDate_default +(mm-1)*100
    dataTime=1100
    CALL grib_set (gribid_dest,'dataDate',dataDate) 
    CALL grib_set (gribid_dest,'dataTime',dataTime)
    CALL grib_set (gribid_dest,'endStep',0)
    CALL grib_set (gribid_dest,'stepUnits',1)
    WRITE (0,*) ndvi_field_mom_meta%shortName,' of ',mm," -> ",TRIM(grib_filename)
    real_buffer(:,:,:)=ndvi_field_mom(:,:,:,mm)
    CALL write_extpar_ICON_field_grib(outfile_id,gribid_dest,ICON_grid,&
             real_buffer,ndvi_field_mom_meta%shortName,dataDate,dataTime)
  END DO

  DO mm=1,12 
  
    dataDate=dataDate_default +(mm-1)*100
    dataTime=1100
    CALL grib_set (gribid_dest,'dataDate',dataDate) 
    CALL grib_set (gribid_dest,'dataTime',dataTime)
    CALL grib_set (gribid_dest,'endStep',0)
    CALL grib_set (gribid_dest,'stepUnits',1)
    WRITE (0,*) ndvi_ratio_mom_meta%shortName,' of ',mm," -> ",TRIM(grib_filename)
    real_buffer(:,:,:)=ndvi_ratio_mom(:,:,:,mm)
    CALL write_extpar_ICON_field_grib(outfile_id,gribid_dest,ICON_grid,&
             real_buffer,ndvi_ratio_mom_meta%shortName,dataDate,dataTime)
  END DO

  DO mm=1,12 

    dataDate=dataDate_default +(mm-1)*100
    dataTime=1100
    CALL grib_set (gribid_dest,'dataDate',dataDate) 
    CALL grib_set (gribid_dest,'dataTime',dataTime)
    CALL grib_set (gribid_dest,'endStep',0)
    CALL grib_set (gribid_dest,'stepUnits',1)
    WRITE (0,*) alb_field_mom_meta%shortName,' of ',mm," -> ",TRIM(grib_filename)
    real_buffer(:,:,:)=alb_field_mom(:,:,:,mm)
    CALL write_extpar_ICON_field_grib(outfile_id,gribid_dest,ICON_grid,&
             real_buffer,alb_field_mom_meta%shortName,dataDate,dataTime)
  END DO

  DO mm=1,12 

    dataDate=dataDate_default +(mm-1)*100
    dataTime=1100
    CALL grib_set (gribid_dest,'dataDate',dataDate) 
    CALL grib_set (gribid_dest,'dataTime',dataTime)
    CALL grib_set (gribid_dest,'endStep',0)
    CALL grib_set (gribid_dest,'stepUnits',1)
    WRITE (0,*) alnid_field_mom_meta%shortName,' of ',mm," -> ",TRIM(grib_filename)
    real_buffer(:,:,:)=alnid_field_mom(:,:,:,mm)
    CALL write_extpar_ICON_field_grib(outfile_id,gribid_dest,ICON_grid,&
            real_buffer,alnid_field_mom_meta%shortName,dataDate,dataTime)
  END DO

  DO mm=1,12 

    dataDate=dataDate_default +(mm-1)*100
    dataTime=1100
    CALL grib_set (gribid_dest,'dataDate',dataDate) 
    CALL grib_set (gribid_dest,'dataTime',dataTime)
    CALL grib_set (gribid_dest,'endStep',0)
    CALL grib_set (gribid_dest,'stepUnits',1)
    WRITE (0,*) aluvd_field_mom_meta%shortName,' of ',mm," -> ",TRIM(grib_filename)
    real_buffer(:,:,:)=aluvd_field_mom(:,:,:,mm)
    CALL write_extpar_ICON_field_grib(outfile_id,gribid_dest,ICON_grid,&
            real_buffer,aluvd_field_mom_meta%shortName,dataDate,dataTime)
  END DO

!Aerosol
    CALL grib_set (gribid_dest,'productDefinitionTemplateNumber', 42)
  DO ntype=1,ntype_aot
    DO mm=1,12
      dataDate=dataDate_default +(mm-1)*100
      dataTime=1100
      CALL grib_set (gribid_dest,'dataDate',dataDate) 
      CALL grib_set (gribid_dest,'dataTime',dataTime)
      CALL grib_set (gribid_dest,'endStep',0)
      CALL grib_set (gribid_dest,'stepUnits',1)
      WRITE (0,*) aot_type_shortname(ntype),' of ',mm," -> ",TRIM(grib_filename)
      real_buffer(:,:,:)=aot_tg(:,:,:,ntype,mm)
      CALL write_extpar_ICON_field_grib(outfile_id,gribid_dest,ICON_grid,&
             real_buffer,aot_type_shortname(ntype),dataDate,dataTime)
    ENDDO
  END DO
 
  CALL grib_set (gribid_dest,'typeOfGeneratingProcess', 196) ! reset typeOfGeneratingProcess to general value

! Warning: Should the last write due to different grib settings
! Land-Use for TILES - now encoded as FR_LAND for different periods
!GLC2000
!!  01   'evergreen broadleaf tree     ' , &    !        evergreen broadleaf forest
!!  02   'deciduous broadleaf tree clos' , &    !        deciduous broadleaf forest
!!  03   'deciduous broadleaf tree open' , &    !        woody savannas & savanna
!!  04   'evergreen needleleaf tree    ' , &    !        evergreen needleleaf forest
!!  05   'deciduous needleleaf tree    ' , &    !        deciduous needleleaf forest
!!  06   'mixed leaf tree              ' , &    !        mixed forest
!!  07   'fresh water flooded tree     ' , &    !        evergreen broadleaf
!!  08   'saline water flooded tree    ' , &    !        evergreen broadleaf
!!  09   'mosaic tree / other nat.veg. ' , &    !        shrubland ?
!!  10   'burnt tree cover             ' , &    ! 
!!  11   'evergreen shrubs closed-open ' , &    !        shrubland open-closed
!!  12   'deciduous shrubs closed-open ' , &    !        savannas
!!  13   'herbaceous cover closed-open ' , &    !        grassland
!!  14   'sparse herbaceous or grass   ' , &    !
!!  15   'flooded shrub or herbaceous  ' , &    !        wetlands
!!  16   'cultivated & managed areas   ' , &    !        croplands
!!  17   'mosaic crop/tree/natural veg.' , &    !        cropland/other vegetation mosaic
!!  18   'mosaic crop/shrub or grass   ' , &    !        cropland/other vegetation mosaic
!!  19   'bare areas                   ' , &    !        barren or sparsely vegetated
!!  20   'water bodies                 ' , &    !        water
!!  21   'snow & ice                   ' , &    !        snow & ice
!!  22   'artificial surfaces          ' , &    !        urban and built-up areas
!!  23   'undefined                    ' /      ! 


!!$  DO mm=1,22
!!$     list(mm) = mm
!!$  ENDDO
!!$     list(23)=39

!GCV2009
!! The globcover dataset contains the following land use classification scheme
! class value         description
! 01  11   'irrigated croplands
! 02  14   'rainfed croplands                             '
! 03  20   'mosaic cropland (50-70%) - vegetation (20-50%)'
! 04  30   'mosaic vegetation (50-70%) - cropland (20-50%)'
! 05  40   'closed broadleaved evergreen forest           '
! 06  50   'closed broadleaved deciduous forest           '
! 07  60   'open broadleaved deciduous forest             '
! 08  70   'closed needleleaved evergreen forest          '
! 09  90   'open needleleaved decid. or evergr. forest    '
! 10  100  'mixed broadleaved and needleleaved forest     '
! 11  110  'mosaic shrubland (50-70%) - grassland (20-50%)'
! 12  120  'mosaic grassland (50-70%) - shrubland (20-50%)'
! 13  130  'closed to open shrubland                      '
! 14  140  'closed to open herbaceous vegetation          '
! 15  150  'sparse vegetation                             '
! 16  160  'closed to open forest regulary flooded        '
! 17  170  'closed forest or shrubland permanently flooded'
! 18  180  'closed to open grassland regularly flooded    '
! 19  190  'artificial surfaces                           '
! 20  200  'bare areas                                    '
! 21  210  'water bodies                                  '
! 22  220  'permanent snow and ice                        '
! 23  230  'undefined                                    
  DO mm=1,5
     list (mm) = mm + 23
  ENDDO
     list (6) = 2
     list (7) = 3
  DO mm=8,13
     list (mm) = mm + 21
  ENDDO
     list(14) = 13
  DO mm=15,18
     list (mm) = mm + 20
  ENDDO
     list(19) = 22
  DO mm=20,22
     list (mm) = mm - 1
  ENDDO
     list(23) = 39

  PRINT *, 'Liste fuer GCV2009'
  PRINT *, '------------------'
  DO mm=1,9
     PRINT *,' ####  ',mm,' #### ',list(mm)
  ENDDO
  DO mm=10,23
     PRINT *,' #### ',mm,' #### ',list(mm)
  ENDDO

!!$  call grib_open_file(outfile,&
!!$       'part.grb2','w')
!!$  CALL grib_new_from_samples(ogrib,"GRIB2")
!!$
!!$  !     a new grib message is loaded from file
!!$  !     igrib is the grib id to be used in subsequent calls

  CALL grib_clone(gribid_dest,ogrib)
  dataDate = 10101
  dataTime = 0
  CALL grib_set (ogrib,'dataDate',dataDate) 
  CALL grib_set (ogrib,'dataTime',dataTime)
  CALL grib_set (ogrib,'stepType','instant')
  CALL grib_set (ogrib,'startStep',0)
  CALL grib_set (ogrib,'endStep',0)
  CALL grib_set (ogrib,'stepUnits',1)

!  CALL grib_set (ogrib,'tablesVersion',11)
  CALL grib_set (ogrib,'productDefinitionTemplateNumber', 53)
  CALL grib_set (ogrib,'discipline',2)
  CALL grib_set (ogrib,'parameterCategory',0)
  CALL grib_set (ogrib,'parameterNumber',36)
  CALL grib_set (ogrib,'NV',0)
  CALL grib_set (ogrib,'PVPresent',0)
  CALL grib_set (ogrib,'typeOfLevel','surface')
  CALL grib_set (ogrib,'partitionTable',243) ! GlobCover, 242-for GLC2000
  CALL grib_set (ogrib,'numberOfPartitions',23)
  CALL grib_set (ogrib,'partitionItems',list)
  CALL grib_set (ogrib,'localInformationNumber',1) !1=ESA GLOBCOVER GCV2009, 2=GLC2000

  DO mm=1,23
    CALL grib_set (ogrib,'partitionNumber',list(mm))
    WRITE(0,*) lu_class_fraction_meta%shortName,' of ',mm," -> ",TRIM(grib_filename)
    real_buffer(:,:,:)=lu_class_fraction(:,:,:,mm)
    CALL write_extpar_ICON_field_grib(outfile_id,ogrib,ICON_grid,&
             real_buffer,lu_class_fraction_meta%shortName,dataDate,dataTime)
  END DO
  CALL grib_release(ogrib)

!!$ dataDate = 10101
!!$ dataTime = 0
!!$ call grib_set (gribid_dest,'dataDate',dataDate) 
!!$ call grib_set (gribid_dest,'dataTime',dataTime)
!!$   call grib_set (gribid_dest,'endStep',0)
!!$   call grib_set (gribid_dest,'stepUnits',1)
!!$    call grib_set (gribid_dest,'stepType','instant')
!!$ call grib_set (gribid_dest,'numberOfGridInReference'     ,2)  ! Triangle corner ! was "=3"
!!$    call grib_set (gribid_dest,'productDefinitionTemplateNumber', 0) ! Reset PDT 
!!$ write(0,*) hh_topo_meta%shortName," -> ",TRIM(grib_filename)
!!$  CALL write_extpar_ICON_field_grib(outfile_id,gribid_dest,ICON_grid,&
!!$         vertex_param%hh_vert,hh_vert_meta%shortName,dataDate,dataTime)


!Aerosol
!!$ DO ntype=1,ntype_aot
!!$   DO mm=1,12
!!$!    call grib_set (gribid_dest,'stepType','avg')
!!$    call grib_set (gribid_dest,'startStep',mm-1)
!!$    call grib_set (gribid_dest,'endStep',mm)
!!$!    call grib_set (gribid_dest,'stepRange',1)
!!$!    call grib_set (gribid_dest,'stepUnits','M')
!!$ write(0,*) aot_type_shortname(ntype),' of ',mm," -> ",TRIM(grib_filename)
!!$  real_buffer(:,:,:)=aot_tg(:,:,:,ntype,mm)
!!$  CALL write_extpar_ICON_field_grib(outfile_id,gribid_dest,ICON_grid,&
!!$             real_buffer,aot_type_shortname(ntype),dataDate,dataTime)
!!$   ENDDO
!!$ ENDDO
!!$! 
!-->




!!$  ! output FIS as well
!!$
!!$
!!$  
!!$  CALL write_extpar_ICON_field_grib(gribid_dest,TRIM(grib_sample), &
!!$  & ICON_grid,lon_geo,lon_geo_meta%shortName,dataDate,dataTime)
!!$  
!!$  CALL write_extpar_ICON_field_grib(gribid_dest,TRIM(grib_sample),&
!!$  & ICON_grid,lat_geo,lat_geo_meta%shortName,dataDate,dataTime)
!---------


!!$    CALL grib_clone(gribid_dest,igribclone) ! clone sample before modifying it
!!$    CALL grib_set(igribclone,'shortName',TRIM(hh_topo_meta%shortName),errorcode) 
!!$
!!$    real_buffer=0._wp
!!$ 
!!$    DO k=1,1
!!$    DO j=1,1
!!$    DO i=1,icon_grid%ncell
!!$       real_buffer(i) =hh_topo(i,j,k)   ! put data to 1D array before putting it to the GRIB record
!!$    ENDDO
!!$    ENDDO
!!$    ENDDO
!!$
!!$    CALL grib_set(igribclone,'values',real_buffer,errorcode) ! put data to GRIB message
!!$    CALL grib_write(igribclone,outfile_id,errorcode) ! write out GRIB message to file
!!$    CALL grib_release(igribclone) ! free memory of grib message

!---------

!!$    CALL grib_clone(gribid_dest,igribclone) ! clone sample before modifying it
!!$    CALL grib_set(igribclone,'shortName',TRIM(fr_land_lu_meta%shortName),errorcode) 
!!$ 
!!$    real_buffer=0._wp
!!$
!!$    DO k=1,1
!!$    DO j=1,1
!!$    DO i=1,icon_grid%ncell
!!$       real_buffer(i) =fr_land_lu(i,j,k)   ! put data to 1D array before putting it to the GRIB record
!!$    ENDDO
!!$    ENDDO
!!$    ENDDO
!!$
!!$    CALL grib_set(igribclone,'values',real_buffer,errorcode) ! put data to GRIB message
!!$    CALL grib_write(igribclone,outfile_id,errorcode) ! write out GRIB message to file
!!$    CALL grib_release(igribclone) ! free memory of grib message
!------------

    CALL grib_close_file(outfile_id)


  END SUBROUTINE write_ICON_grid_extpar_grib
  !-----------------------------------------------------------------------
    

  SUBROUTINE decode_uuid (uuid_str, uuid)
    CHARACTER(LEN=*), INTENT(IN)  :: uuid_str   ! uuid encoded as string
    CHARACTER(LEN=1), INTENT(OUT) :: uuid(:)    ! decoded uuid

    INTEGER          :: i, j, l, n, b
    CHARACTER(LEN=2) :: buf

    uuid(:) = ACHAR (0)
    l = VERIFY (uuid_str, "0123456789ABCDEFabcdef-")
    IF (l > 0) THEN
      WRITE (0,*) "Warning: invalid character in uuid: '", uuid_str(l:l),"'"
      RETURN
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
      IF (j > SIZE (uuid)) CALL finish ("decode_uuid", "uuid input too long!")
      uuid(j) = ACHAR (b)
    END DO
    IF (i == n) CALL finish ("decode_uuid", "uuid bad length")
  END SUBROUTINE decode_uuid
 
END Module mo_extpar_output_grib


