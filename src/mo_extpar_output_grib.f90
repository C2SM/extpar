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

  USE mo_io_units,          ONLY: filename_max


  !> abort_extpar defined in MODULE utilities_extpar
  USE mo_utilities_extpar, ONLY: abort_extpar
  USE mo_utilities_extpar, ONLY: get_rot_spol_coor


  USE mo_glc2000_lookup_tables, ONLY: nclass_glc2000
  USE mo_ndvi_data, ONLY: ntime_ndvi
  USE mo_aot_data, ONLY: ntype_aot, ntime_aot
  USE mo_albedo_data, ONLY: ntime_alb

  USE mo_soil_data,   ONLY: soil_data, FAO_data, HWSD_data

  USE grib_api 
  USE mo_io_grib_api



  IMPLICIT NONE

  PRIVATE

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
    &                                     crutemp,       &
    &                                     alb_field_mom, &
    &                                     fr_sand,       & 
    &                                     fr_silt,       &
    &                                     fr_clay,       &
    &                                     fr_oc,         &
    &                                     fr_bd,         &
    &                                     fr_dm,         &
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
    &                        aot_tg_meta, &
    &                        def_aot_tg_meta
 USE mo_var_meta_data, ONLY: aot_type_shortname

 USE mo_var_meta_data, ONLY: crutemp_meta, &
    &                        def_crutemp_meta

  USE mo_physical_constants, ONLY: grav


  CHARACTER (len=*), INTENT(IN) :: grib_filename !< filename for the grib file
  CHARACTER (len=*), INTENT(IN) :: grib_sample  !< name for grib sample  (sample to be found in $GRIB_SAMPLES_PATH)

  TYPE(rotated_lonlat_grid), INTENT(IN) :: cosmo_grid !< structure which contains the definition of the COSMO grid
  TYPE(target_grid_def), INTENT(IN) :: tg !< structure with target grid description
  INTEGER (KIND=i4), INTENT(IN)     :: isoil_data
  LOGICAL,         INTENT(IN)       :: ldeep_soil 
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
  REAL(KIND=wp), INTENT(IN)  :: crutemp(:,:,:)  !< cru climatological temperature , crutemp(ie,je,ke) 
  REAL(KIND=wp), INTENT(IN), OPTIONAL :: fr_sand(:,:,:)   !< sand fraction due to HWSD
  REAL(KIND=wp), INTENT(IN), OPTIONAL :: fr_silt(:,:,:)   !< silt fraction due to HWSD
  REAL(KIND=wp), INTENT(IN), OPTIONAL :: fr_clay(:,:,:)   !< clay fraction due to HWSD
  REAL(KIND=wp), INTENT(IN), OPTIONAL :: fr_oc(:,:,:)   !< oc fraction due to HWSD
  REAL(KIND=wp), INTENT(IN), OPTIONAL :: fr_bd(:,:,:)   !< bulk density due to HWSD
  REAL(KIND=wp), INTENT(IN), OPTIONAL :: fr_dm(:,:,:)   !< bulk density due to HWSD

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

  INTEGER (KIND=i8)  :: dataDate
  INTEGER (KIND=i8)  :: dataTime
  INTEGER :: mm ! month
  INTEGER :: ntype ! type of aerosol



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
  CALL def_topo_meta(dim_3d_tg,diminfohor=dim_4d_tg)
  !  hh_topo_meta, fr_land_topo_meta, &
  !         stdh_topo_meta, theta_topo_meta, &
  !         aniso_topo_meta, slope_topo_meta, &
  !         hh_vert_meta, npixel_vert_meta
  ELSE
    CALL def_topo_meta(dim_3d_tg)
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
  & cosmo_grid,fr_land_lu,fr_land_lu_meta%shortName,dataDate,dataTime)

  PRINT *,'hh_topo_meta%shortName: ',hh_topo_meta%shortName

  CALL write_extpar_cosmo_field_grib(outfile_id,TRIM(grib_sample),&
  & cosmo_grid,hh_topo,hh_topo_meta%shortName,dataDate,dataTime)
  
  ! output FIS as well
  PRINT *,'hh_fis_meta%shortName: ',hh_fis_meta%shortName
  CALL write_extpar_cosmo_field_grib(outfile_id,TRIM(grib_sample),&
  & cosmo_grid,hh_topo*grav,hh_fis_meta%shortName,dataDate,dataTime)

  CALL write_extpar_cosmo_field_grib(outfile_id,TRIM(grib_sample),&
  & cosmo_grid,stdh_topo,stdh_topo_meta%shortName,dataDate,dataTime)
  
  IF (lsso) THEN
  CALL write_extpar_cosmo_field_grib(outfile_id,TRIM(grib_sample),&
  & cosmo_grid,theta_topo,theta_topo_meta%shortName,dataDate,dataTime)
  ENDIF

  IF (lsso) THEN
  CALL write_extpar_cosmo_field_grib(outfile_id,TRIM(grib_sample),&
  & cosmo_grid,aniso_topo,aniso_topo_meta%shortName,dataDate,dataTime)
  ENDIF

  IF (lsso) THEN
  CALL write_extpar_cosmo_field_grib(outfile_id,TRIM(grib_sample),&
  & cosmo_grid,slope_topo,slope_topo_meta%shortName,dataDate,dataTime)
  ENDIF

  IF (lrad) THEN
    CALL write_extpar_cosmo_field_grib(outfile_id,TRIM(grib_sample),&
    & cosmo_grid,slope_asp_topo,slope_asp_topo_meta%shortName,dataDate,dataTime)
  ENDIF
  IF (lrad) THEN
    CALL write_extpar_cosmo_field_grib(outfile_id,TRIM(grib_sample),&
    & cosmo_grid,slope_ang_topo,slope_ang_topo_meta%shortName,dataDate,dataTime)
  ENDIF
!ROATODO: check 4D
  IF (lrad) THEN
    DO mm=1,nhori
      extpar_cosmo_buffer(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1:1) = &
      & horizon_topo(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1:1,mm) 
      CALL write_extpar_cosmo_field_grib(outfile_id,TRIM(grib_sample),&
      & cosmo_grid,extpar_cosmo_buffer,horizon_topo_meta%shortName,dataDate,dataTime)
    ENDDO
  ENDIF
  IF (lrad) THEN
    CALL write_extpar_cosmo_field_grib(outfile_id,TRIM(grib_sample),&
    & cosmo_grid,skyview_topo,skyview_topo_meta%shortName,dataDate,dataTime)
  ENDIF
  
  CALL write_extpar_cosmo_field_grib(outfile_id,TRIM(grib_sample),&
  & cosmo_grid,plcov_mn_lu,plcov_mn_lu_meta%shortName,dataDate,dataTime)

  CALL write_extpar_cosmo_field_grib(outfile_id,TRIM(grib_sample),&
  & cosmo_grid,plcov_mx_lu,plcov_mx_lu_meta%shortName,dataDate,dataTime)
  
  CALL write_extpar_cosmo_field_grib(outfile_id,TRIM(grib_sample), &
  & cosmo_grid,lai_mn_lu,lai_mn_lu_meta%shortName,dataDate,dataTime)

  CALL write_extpar_cosmo_field_grib(outfile_id,TRIM(grib_sample),&
  & cosmo_grid,lai_mx_lu,lai_mx_lu_meta%shortName,dataDate,dataTime)

  CALL write_extpar_cosmo_field_grib(outfile_id,TRIM(grib_sample), &
  & cosmo_grid,rs_min_lu,rs_min_lu_meta%shortName,dataDate,dataTime)
  
  CALL write_extpar_cosmo_field_grib(outfile_id,TRIM(grib_sample), &
  & cosmo_grid,for_e_lu,for_e_lu_meta%shortName,dataDate,dataTime)
  
  CALL write_extpar_cosmo_field_grib(outfile_id,TRIM(grib_sample), &
  & cosmo_grid,for_d_lu,for_d_lu_meta%shortName,dataDate,dataTime)
  
  CALL write_extpar_cosmo_field_grib(outfile_id,TRIM(grib_sample), &
  & cosmo_grid,emissivity_lu,emissivity_lu_meta%shortName,dataDate,dataTime)

  CALL write_extpar_cosmo_field_grib(outfile_id,TRIM(grib_sample),&
  & cosmo_grid,root_lu,root_lu_meta%shortName,dataDate,dataTime)

  CALL write_extpar_cosmo_field_grib(outfile_id,TRIM(grib_sample),&
  & cosmo_grid,z0_lu,z0_lu_meta%shortName,dataDate,dataTime)

  CALL write_extpar_cosmo_field_grib(outfile_id,TRIM(grib_sample),&
  & cosmo_grid,lake_depth,lake_depth_meta%shortName,dataDate,dataTime)
  
  CALL write_extpar_cosmo_field_grib(outfile_id,TRIM(grib_sample),&
  & cosmo_grid,fr_lake,fr_lake_meta%shortName,dataDate,dataTime)

  CALL write_extpar_cosmo_field_grib(outfile_id,TRIM(grib_sample),&
  & cosmo_grid,soiltype_fao,soiltype_fao_meta%shortName,dataDate,dataTime)

  IF(soil_data == HWSD_data) THEN
  CALL write_extpar_cosmo_field_grib(outfile_id,TRIM(grib_sample),&
  & cosmo_grid,fr_sand,HWSD_SAND_meta%shortName,dataDate,dataTime)
  ENDIF

  IF(soil_data == HWSD_data) THEN
  CALL write_extpar_cosmo_field_grib(outfile_id,TRIM(grib_sample),&
  & cosmo_grid,fr_silt,HWSD_SILT_meta%shortName,dataDate,dataTime)
  ENDIF

  IF(soil_data == HWSD_data) THEN
  CALL write_extpar_cosmo_field_grib(outfile_id,TRIM(grib_sample),&
  & cosmo_grid,fr_clay,HWSD_CLAY_meta%shortName,dataDate,dataTime)
  ENDIF

  IF(soil_data == HWSD_data) THEN
  CALL write_extpar_cosmo_field_grib(outfile_id,TRIM(grib_sample),&
  & cosmo_grid,fr_oc,HWSD_OC_meta%shortName,dataDate,dataTime)
  ENDIF

  IF(soil_data == HWSD_data) THEN
  CALL write_extpar_cosmo_field_grib(outfile_id,TRIM(grib_sample),&
  & cosmo_grid,fr_bd,HWSD_BD_meta%shortName,dataDate,dataTime)
  ENDIF

  IF(soil_data == HWSD_data) THEN
  CALL write_extpar_cosmo_field_grib(outfile_id,TRIM(grib_sample),&
  & cosmo_grid,fr_dm,HWSD_DM_meta%shortName,dataDate,dataTime)
  ENDIF


  CALL write_extpar_cosmo_field_grib(outfile_id,TRIM(grib_sample),&
  & cosmo_grid,crutemp,crutemp_meta%shortName,dataDate,dataTime)
  
  CALL write_extpar_cosmo_field_grib(outfile_id,TRIM(grib_sample),&
  & cosmo_grid,lon_geo,lon_geo_meta%shortName,dataDate,dataTime)
  
  CALL write_extpar_cosmo_field_grib(outfile_id,TRIM(grib_sample),&
  & cosmo_grid,lat_geo,lat_geo_meta%shortName,dataDate,dataTime)

 ! urban_lu
 ! ndvi_max
 CALL write_extpar_cosmo_field_grib(outfile_id,TRIM(grib_sample),&
 & cosmo_grid,ndvi_max,ndvi_max_meta%shortName,dataDate,dataTime)

 ! ndvi_field_mom
 DO mm=1,12
    ! get dataDate and dataTime according to DWD convention for external parameters of climatological monthly mean fields
    CALL set_date_mm_extpar_field(mm,dataDate,dataTime)
    extpar_cosmo_buffer(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1:1) = &
  & ndvi_field_mom(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1:1,mm) 
    CALL write_extpar_cosmo_field_grib(outfile_id,TRIM(grib_sample),&
    & cosmo_grid,extpar_cosmo_buffer,ndvi_field_mom_meta%shortName,dataDate,dataTime)
 ENDDO

 ! ndvi_ratio_mom
 DO mm=1,12
    ! get dataDate and dataTime according to DWD convention for external parameters of climatological monthly mean fields
    CALL set_date_mm_extpar_field(mm,dataDate,dataTime)
    extpar_cosmo_buffer(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1:1) = &
  & ndvi_ratio_mom(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1:1,mm) 
    CALL write_extpar_cosmo_field_grib(outfile_id,TRIM(grib_sample),&
    & cosmo_grid,extpar_cosmo_buffer,ndvi_ratio_mom_meta%shortName,dataDate,dataTime)
 ENDDO

 DO mm=1,12
    CALL set_date_mm_extpar_field(mm,dataDate,dataTime)
    extpar_cosmo_buffer(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1:1) = &
  & alb_field_mom(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1:1,mm)
    CALL write_extpar_cosmo_field_grib(outfile_id,TRIM(grib_sample),&
    & cosmo_grid,extpar_cosmo_buffer,alb_field_mom_meta%shortName,dataDate,dataTime)
 ENDDO 

 ! aot_tg

 DO ntype=1,ntype_aot
   DO mm=1,12
   ! get dataDate and dataTime according to DWD convention for external parameters of climatological monthly mean fields
    CALL set_date_mm_extpar_field(mm,dataDate,dataTime)
    extpar_cosmo_buffer(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1:1) = &
  & aot_tg(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1:1,ntype,mm)
    CALL write_extpar_cosmo_field_grib(outfile_id,TRIM(grib_sample), &
    & cosmo_grid,extpar_cosmo_buffer,aot_type_shortname(ntype),dataDate,dataTime)
   ENDDO
 ENDDO
 ! ice_lu


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
  CALL def_topo_meta(dim_3d_tg)
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
  & gme_grid,fr_land_lu,TRIM(fr_land_lu_meta%shortName),dataDate,dataTime)

  PRINT *,'hh_topo_meta%shortName: ',hh_topo_meta%shortName

  CALL write_extpar_gme_field_grib(outfile_id,TRIM(grib_sample),&
  & gme_grid,hh_topo,hh_topo_meta%shortName,dataDate,dataTime)
  ! output FIS as well
  PRINT *,'hh_fis_meta%shortName: ',hh_fis_meta%shortName
  CALL write_extpar_gme_field_grib(outfile_id,TRIM(grib_sample),&
  & gme_grid,hh_topo*grav,hh_fis_meta%shortName,dataDate,dataTime)

  CALL write_extpar_gme_field_grib(outfile_id,TRIM(grib_sample),&
  & gme_grid,stdh_topo,stdh_topo_meta%shortName,dataDate,dataTime)
  
  IF (lsso) THEN
  CALL write_extpar_gme_field_grib(outfile_id,TRIM(grib_sample),&
  & gme_grid,theta_topo,theta_topo_meta%shortName,dataDate,dataTime)
  ENDIF

  IF (lsso) THEN
  CALL write_extpar_gme_field_grib(outfile_id,TRIM(grib_sample),&
  & gme_grid,aniso_topo,aniso_topo_meta%shortName,dataDate,dataTime)
  ENDIF

  IF (lsso) THEN
  CALL write_extpar_gme_field_grib(outfile_id,TRIM(grib_sample), &
  & gme_grid,slope_topo,slope_topo_meta%shortName,dataDate,dataTime)
  ENDIF

  CALL write_extpar_gme_field_grib(outfile_id,TRIM(grib_sample),&
  & gme_grid,plcov_mn_lu,plcov_mn_lu_meta%shortName,dataDate,dataTime)
  
  CALL write_extpar_gme_field_grib(outfile_id,TRIM(grib_sample),&
  & gme_grid,plcov_mx_lu,plcov_mx_lu_meta%shortName,dataDate,dataTime)
  
  CALL write_extpar_gme_field_grib(outfile_id,TRIM(grib_sample), &
  & gme_grid,lai_mn_lu,lai_mn_lu_meta%shortName,dataDate,dataTime)
  
  CALL write_extpar_gme_field_grib(outfile_id,TRIM(grib_sample),&
  & gme_grid,lai_mx_lu,lai_mx_lu_meta%shortName,dataDate,dataTime)

  CALL write_extpar_gme_field_grib(outfile_id,TRIM(grib_sample), &
  & gme_grid,rs_min_lu,rs_min_lu_meta%shortName,dataDate,dataTime)
  
  CALL write_extpar_gme_field_grib(outfile_id,TRIM(grib_sample), &
  & gme_grid,for_e_lu,for_e_lu_meta%shortName,dataDate,dataTime)
  
  CALL write_extpar_gme_field_grib(outfile_id,TRIM(grib_sample), &
  & gme_grid,for_d_lu,for_d_lu_meta%shortName,dataDate,dataTime)
  
  CALL write_extpar_gme_field_grib(outfile_id,TRIM(grib_sample), &
  & gme_grid,emissivity_lu,emissivity_lu_meta%shortName,dataDate,dataTime)
  
  CALL write_extpar_gme_field_grib(outfile_id,TRIM(grib_sample), &
  & gme_grid,root_lu,root_lu_meta%shortName,dataDate,dataTime)

  CALL write_extpar_gme_field_grib(outfile_id,TRIM(grib_sample), &
  & gme_grid,z0_lu,z0_lu_meta%shortName,dataDate,dataTime)

  CALL write_extpar_gme_field_grib(outfile_id,TRIM(grib_sample),&
  & gme_grid,lake_depth,lake_depth_meta%shortName,dataDate,dataTime)
  
  CALL write_extpar_gme_field_grib(outfile_id,TRIM(grib_sample),&
  & gme_grid,fr_lake,fr_lake_meta%shortName,dataDate,dataTime)

  CALL write_extpar_gme_field_grib(outfile_id,TRIM(grib_sample),&
  & gme_grid,soiltype_fao,soiltype_fao_meta%shortName,dataDate,dataTime)

  CALL write_extpar_gme_field_grib(outfile_id,TRIM(grib_sample),&
  & gme_grid,crutemp,crutemp_meta%shortName,dataDate,dataTime)
  
  CALL write_extpar_gme_field_grib(outfile_id,TRIM(grib_sample), &
  & gme_grid,lon_geo,lon_geo_meta%shortName,dataDate,dataTime)
  
  CALL write_extpar_gme_field_grib(outfile_id,TRIM(grib_sample),&
  & gme_grid,lat_geo,lat_geo_meta%shortName,dataDate,dataTime)


 ! urban_lu
 ! ndvi_max
 CALL write_extpar_gme_field_grib(outfile_id,TRIM(grib_sample), &
 & gme_grid,ndvi_max,ndvi_max_meta%shortName,dataDate,dataTime)
 ! ndvi_field_mom
 
 DO mm=1,12
    ! get dataDate and dataTime according to DWD convention for external parameters of climatological monthly mean fields
    CALL set_date_mm_extpar_field(mm,dataDate,dataTime)
    real_buffer(1:gme_grid%nip1,1:gme_grid%nip1,1:gme_grid%nd) = &
  & ndvi_field_mom(1:gme_grid%nip1,1:gme_grid%nip1,1:gme_grid%nd,mm)
    CALL write_extpar_gme_field_grib(outfile_id,TRIM(grib_sample),&
    & gme_grid,real_buffer,ndvi_field_mom_meta%shortName,dataDate,dataTime)
 ENDDO
 ! ndvi_ratio_mom
 DO mm=1,12
    ! get dataDate and dataTime according to DWD convention for external parameters of climatological monthly mean fields
    CALL set_date_mm_extpar_field(mm,dataDate,dataTime)
    real_buffer(1:gme_grid%nip1,1:gme_grid%nip1,1:gme_grid%nd) = &
  & ndvi_ratio_mom(1:gme_grid%nip1,1:gme_grid%nip1,1:gme_grid%nd,mm) 
    CALL write_extpar_gme_field_grib(outfile_id,TRIM(grib_sample),&
    & gme_grid,real_buffer,ndvi_ratio_mom_meta%shortName,dataDate,dataTime)
 ENDDO

 ! alb_field_mom
 DO mm=1,12
    CALL set_date_mm_extpar_field(mm,dataDate,dataTime)
    real_buffer(1:gme_grid%nip1,1:gme_grid%nip1,1:gme_grid%nd) = &
  & alb_field_mom(1:gme_grid%nip1,1:gme_grid%nip1,1:gme_grid%nd,mm)
    CALL write_extpar_gme_field_grib(outfile_id,TRIM(grib_sample),&
    & gme_grid,real_buffer,alb_field_mom_meta%shortName,dataDate,dataTime)
 ENDDO

 ! aot_tg

 DO ntype=1,ntype_aot
   DO mm=1,12
   ! get dataDate and dataTime according to DWD convention for external parameters of climatological monthly mean fields
    CALL set_date_mm_extpar_field(mm,dataDate,dataTime)
    real_buffer(1:gme_grid%nip1,1:gme_grid%nip1,1:gme_grid%nd) = &
  & aot_tg(1:gme_grid%nip1,1:gme_grid%nip1,1:gme_grid%nd,ntype,mm)
    CALL write_extpar_gme_field_grib(outfile_id,TRIM(grib_sample), &
    & gme_grid,real_buffer,aot_type_shortname(ntype),dataDate,dataTime)
   ENDDO
 ENDDO
 ! ice_lu


  CALL grib_close_file(outfile_id)



  END SUBROUTINE write_gme_grid_extpar_grib
  !-----------------------------------------------------------------------



    
 
END Module mo_extpar_output_grib


