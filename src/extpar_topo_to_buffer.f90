 !+ Fortran main program to read in GLOBE orography data and aggregate to target grid
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_0         2010/12/21 Hermann Asensio
!  Initial release
! V1_1         2011/01/20 Hermann Asensio
!  small bug fixes accroding to Fortran compiler warnings        
! V1_2         2011/03/25 Hermann Asensio
!  update to support ICON refinement grids
! V1_4         2011/04/21 Anne Roches     
!  implementation of orography smoothing  
! V1_7         2013/01/25 Guenther Zaengl 
!  Parallel threads for ICON and COSMO using Open-MP, 
!  Several bug fixes and optimizations for ICON search algorithm, 
!  particularly for the special case of non-contiguous domains; 
!  simplified namelist control for ICON  
! V2_0         2013/06/04 Martina Messmer, Anne Roches
!  introduction of the ASTER topography raw data set for external parameters
!  switch to choose if SSO parameters are desired or not
! V2_0         2013/06/04 Anne Roches
!  Implementation of the topographical corrected radiation parameters
! V1_14        2014-07-18 Juergen Helmert
!  Combined COSMO Release
! V2_1         2015-01-12 Juergen Helmert 
!  Bugfix correction covers CSCS SVN r5907-r6359
! V2_6         2016-10-07 Juergen Helmert
!  Add namelist switch lfilter_topo         
! V2_10        2018-02-19 Juergen Helmert 
!  lsubtract_mean_slope, ERA-I surface temp for land points          
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!> Fortran main program to read in GLOBE orography data and aggregate to target grid
!>  
!! @par extpar_topo_to_buffer 
!!
!! This program reads in the GLOBE/ASTER orography data set and aggregates the orographic height to the target grid 
!! and computes the subgrid-scale orography parameters (SSO) required by the SSO-parameterization.
!!
!> Purpose: read in GLOBE/ASTER orography data and aggregate to COSMO grid
!> \author Hermann Asensio
PROGRAM extpar_topo_to_buffer

  ! Load the library information data:
  USE info_extpar, ONLY: info_define, info_print


  !> kind parameters are defined in MODULE data_parameters
  USE mo_kind, ONLY: wp
  USE mo_kind, ONLY: i8
  USE mo_kind, ONLY: i4


  USE mo_target_grid_data, ONLY: lon_geo, &
    &                            lat_geo, &
    &                            no_raw_data_pixel

  USE mo_target_grid_data, ONLY: tg  !< structure with target grid description

  USE mo_target_grid_routines, ONLY: init_target_grid

  USE mo_grid_structures, ONLY: target_grid_def,  &
    &                            reg_lonlat_grid,  &
    &                            rotated_lonlat_grid
  USE mo_grid_structures, ONLY: igrid_icon
  USE mo_grid_structures, ONLY: igrid_cosmo
  USE mo_grid_structures, ONLY: igrid_gme

  USE mo_cosmo_grid, ONLY: COSMO_grid, &
    &                       lon_rot, &
    &                       lat_rot, &
    &                       allocate_cosmo_rc, &
    &                       get_cosmo_grid_info, &
    &                       calculate_cosmo_domain_coordinates

  USE mo_icon_grid_data, ONLY: ICON_grid    !< structure which contains the definition of the ICON grid
    
  USE mo_icon_grid_data, ONLY: icon_grid_region

  USE mo_base_geometry,  ONLY: geographical_coordinates, &
    &                          cartesian_coordinates

  USE mo_additional_geometry,  ONLY: cc2gc,            &
    &                                gc2cc,            &
    &                                arc_length,       &
    &                                cos_arc_length,   &
    &                                inter_section,    &
    &                                vector_product,   &
    &                                point_in_polygon_sp

  USE mo_icon_domain,          ONLY: icon_domain,          &
    &                                grid_cells,           &
    &                                grid_vertices,        &
    &                                construct_icon_domain,&
    &                                destruct_icon_domain

  USE mo_io_units,          ONLY: filename_max

  USE mo_exception,         ONLY: message_text, message, finish

  USE mo_utilities_extpar, ONLY: abort_extpar

  USE mo_math_constants,  ONLY: pi, pi_2, dbl_eps,rad2deg

  USE mo_topo_routines, ONLY: read_namelists_extpar_orography
  USE mo_topo_routines, ONLY: read_namelists_extpar_scale_sep

  USE mo_topo_tg_fields, ONLY:  fr_land_topo,    &
    &                         hh_topo,            &
    &                         stdh_topo,          &
    &                         theta_topo,         &
    &                         aniso_topo,         &
    &                         slope_topo,         &
    &                         z0_topo,             &
    &                         allocate_topo_target_fields,&
    &                         slope_asp_topo,     &
    &                         slope_ang_topo,     &
    &                         horizon_topo,       &
    &                         skyview_topo

  USE mo_topo_tg_fields, ONLY: add_parameters_domain, &
    &                           vertex_param,          &
    &                           allocate_additional_hh_param

! mes > -------------------------------------------------------------
  USE mo_topo_data,      ONLY:  topo_aster,        &
    &                           topo_gl,           &
    &                           itopo_type,        &    
    &                           topo_tiles_grid,   &
    &                           topo_grid,         &
    &                           ntiles,            &
    &                           max_tiles,         &
    &                           nc_tot,            &
    &                           nr_tot,            &
    &                           nc_tile,           &
    &                           tiles_lon_min,     &
    &                           tiles_lon_max,     &
    &                           tiles_lat_min,     &
    &                           tiles_lat_max,     &
    &                           aster_lat_min,     &
    &                           aster_lat_max,     &
    &                           aster_lon_min,     &
    &                           aster_lon_max,     &
    &                           num_tiles,         &
    &                           allocate_topo_data,&
    &                           fill_topo_data,    &
    &                           lradtopo,          &
    &                           nhori,             &
    &                           deallocate_topo_fields

! mes < -------------------------------------------------------------

  USE mo_topo_routines, ONLY:  read_topo_data_input_namelist, &
    &                           det_topo_tiles_grid,           &
    &                           det_topo_grid,                  &
    &                           get_topo_tile_nr,              &
    &                           get_topo_tile_block_indices

  USE mo_agg_topo_icon,  ONLY: agg_topo_data_to_target_grid_icon
  USE mo_agg_topo_cosmo, ONLY: agg_topo_data_to_target_grid_cosmo

  USE mo_topo_output_nc, ONLY: write_netcdf_buffer_topo
  USE mo_topo_output_nc, ONLY: write_netcdf_icon_grid_topo
  USE mo_topo_output_nc, ONLY: write_netcdf_cosmo_grid_topo
!roa >
  USE mo_oro_filter, ONLY: read_namelists_extpar_orosmooth
  USE mo_lradtopo,   ONLY: read_namelists_extpar_lradtopo, &
    &                      compute_lradtopo
!roa <

  IMPLICIT NONE

  CHARACTER(len=filename_max) :: filename
  CHARACTER(len=filename_max) :: netcdf_filename
  CHARACTER(len=filename_max) :: namelist_grid_def

  CHARACTER(len=filename_max) :: input_namelist_file
  CHARACTER (len=filename_max) :: input_namelist_cosmo_grid !< file with input namelist with COSMO grid definition
  CHARACTER (len=filename_max) :: namelist_topo_data_input !< file with input namelist with GLOBE data information
  CHARACTER (len=filename_max) :: namelist_scale_sep_data_input!< file with input namelist with scale separated data information
!roa >
  CHARACTER (len=filename_max) :: namelist_oro_smooth       !< file with orography smoothing information (switches)
  CHARACTER (len=filename_max) :: namelist_lrad             !< file with opo information (switches)
!roa <  

    
  CHARACTER (len=filename_max) :: raw_data_path        !< path to raw data
  CHARACTER (LEN=filename_max) :: topo_files(1:max_tiles)  !< filenames globe raw data

  CHARACTER (len=filename_max) :: orography_buffer_file !< name for orography buffer file
  CHARACTER (len=filename_max) :: orography_output_file !< name for orography output file

  CHARACTER (len=filename_max) :: raw_data_orography_path        !< path to raw data
  CHARACTER (len=filename_max) :: raw_data_scale_sep_orography_path !< path to raw data
  CHARACTER (LEN=filename_max) :: scale_sep_files(1:max_tiles)  !< filenames globe raw data
      
  CHARACTER (len=filename_max) :: netcdf_out_filename      !< filename for netcdf file with GLOBE data on COSMO grid

  REAL (KIND=wp) :: point_lon_geo !< longitude of a point in geographical system
  REAL (KIND=wp) :: point_lat_geo !< latitude of a point in geographical system

  REAL(KIND=wp) :: undefined !< value to indicate undefined grid elements 
  INTEGER (KIND=i4) :: undef_int   !< value for undefined integer
  INTEGER (KIND=i4) :: default_value !< default value
  INTEGER (KIND=i4) :: index_tile   !< index for dummy test
  TYPE(geographical_coordinates) :: DWD_location !< geographical coordinates of DWD for dummy test
  
!_br 21.02.14 begin (some compilers do not like variable dimensions)
!  INTEGER (KIND=i4) :: topo_startrow(1:ntiles)    !< startrow indeces for each GLOBE tile
!  INTEGER (KIND=i4) :: topo_endrow(1:ntiles)      !< endrow indeces for each GLOBE tile
!  INTEGER (KIND=i4) :: topo_startcolumn(1:ntiles)  !< starcolumn indeces for each GLOBE tile

!  INTEGER (KIND=i4) :: topo_endcolumn(1:ntiles)   !< endcolumn indeces for each GLOBE tile

  INTEGER (KIND=i4), ALLOCATABLE :: topo_startrow(:)    !< startrow indeces for each GLOBE tile
  INTEGER (KIND=i4), ALLOCATABLE :: topo_endrow(:)      !< endrow indeces for each GLOBE tile
  INTEGER (KIND=i4), ALLOCATABLE :: topo_startcolumn(:)  !< starcolumn indeces for each GLOBE tile
  INTEGER (KIND=i4), ALLOCATABLE :: topo_endcolumn(:)   !< endcolumn indeces for each GLOBE tile
!_br 21.02.14 end

  TYPE(geographical_coordinates) :: ur   !< upper right point for test block
  TYPE(geographical_coordinates) :: ll   !< lower left point for test block

  INTEGER :: k !< counter
  INTEGER :: ie !< counter
  INTEGER :: je !< counter
  INTEGER :: ke !< counter

  INTEGER (KIND=i4) :: startrow_index = 0      !< the index of the GLOBE data row of the first data block row
  INTEGER (KIND=i4) :: endrow_index  = 0       !< the index of the GLOBE data row of the last data block row
  INTEGER (KIND=i4) :: startcolumn_index    !< the index of the startcolumn of data to read in
  INTEGER (KIND=i4) :: endcolumn_index      !< the index of the endcolumn of data to read in
  INTEGER  (KIND=i4) :: point_lon_index !< longitude index of point for regular lon-lat grid
  INTEGER  (KIND=i4) :: point_lat_index !< latitude index of point for regular lon-lat grid

  REAL (KIND=wp) :: topo_target_value  !< interpolated altitude from GLOBE data

  INTEGER (KIND=i4) :: igrid_type  !< target grid type, 1 for ICON, 2 for COSMO, 3 for GME grid

  ! variables for the ICON grid 
  INTEGER :: nvertex  !< total number of vertices
  INTEGER :: nv ! counter

  REAL :: timestart
  REAL :: timeend
  REAL :: timediff

 !mes > -------------------------
 INTEGER (KIND=i4) :: ntiles_column        !< number of tile columns in total domain
 INTEGER (KIND=i4) :: ntiles_row           !< number of tile rows in total domain
 LOGICAL           ::lsso_param
 LOGICAL           :: lscale_separation
 LOGICAL           ::lfilter_topo
 LOGICAL           ::lsubtract_mean_slope
 !mes <

  !roa >
  LOGICAL           :: &
    lfilter_oro,     &
    lxso_first
  
  INTEGER(KIND=i4)  :: &
    ilow_pass_oro,   & 
    numfilt_oro,     &
    ifill_valley,    &
    ilow_pass_xso,   &
    numfilt_xso

  REAL(KIND=wp)     :: &
    eps_filter,      &
    rfill_valley,    &   
    rxso_mask     
!roa <


  ALLOCATE (topo_startrow(1:ntiles), topo_endrow(1:ntiles),topo_startcolumn(1:ntiles),topo_endcolumn(1:ntiles))  !_br 21.02.14
!_br 21.02.14 for clean programming this should be deallocated somewhere

 ! Print the default information to stdout:
  CALL info_define ('topo_to_buffer')      ! Pre-define the program name as binary name
  CALL info_print ()                     ! Print the information to stdout
  !--------------------------------------------------------------------------------------------------------
 
  namelist_grid_def = 'INPUT_grid_org'

  namelist_lrad     = 'INPUT_RADTOPO'

  CALL read_namelists_extpar_lradtopo(namelist_lrad,lradtopo,nhori)

  CALL init_target_grid(namelist_grid_def,lrad=lradtopo)

  PRINT *,' target grid tg: ',tg%ie, tg%je, tg%ke, tg%minlon, tg%maxlon, tg%minlat, tg%maxlat
  igrid_type = tg%igrid_type
  ! get information on target grid
  ! get GLOBE raw data information

  ! Checks
   IF (igrid_type /= igrid_cosmo) THEN 
      lradtopo = .FALSE.
      lfilter_oro = .FALSE.
   print*,"The model grid is not COSMO: Assume "
   print*,"lradtopo = .FALSE. "
   print*,"lfilter_oro = .FALSE. "
   END IF
 

  !--------------------------------------------------------------------------------------------------------

  ! read namelist with globe data information

  namelist_topo_data_input = 'INPUT_ORO'
  CALL read_namelists_extpar_orography(namelist_topo_data_input, &
    &                                  raw_data_orography_path,   &
    &                                  topo_files,                &  !mes>
    &                                  ntiles_column,             &
    &                                  ntiles_row,                &
    &                                  itopo_type,                &
    &                                  lsso_param,                &
    &                                  lfilter_topo,               &
    &                                  lsubtract_mean_slope,     &
    &                                  orography_buffer_file,     &
    &                                  orography_output_file)

!< *mes
  namelist_scale_sep_data_input = 'INPUT_SCALE_SEP'
  CALL read_namelists_extpar_scale_sep(namelist_scale_sep_data_input,     &
    &                                  raw_data_scale_sep_orography_path, &
    &                                  scale_sep_files,                   &
    &                                  lscale_separation) 

  IF (lscale_separation .AND. itopo_type.eq.2) THEN   !_br 21.02.14 replaced eq by eqv
    lscale_separation = .FALSE.
    PRINT*, '*** Scale separation can only be used with GLOBE as raw topography ***'
  ENDIF
!> *mes>

  CALL num_tiles(itopo_type, ntiles_column, ntiles_row, ntiles)        
 ! gives back the number of tiles that are available 16 for GLOBE or 36 for ASTER
  
!mes <

! mes >
  CALL allocate_topo_data(ntiles)                  ! allocates the data using ntiles

  CALL fill_topo_data(raw_data_orography_path,topo_files, &! the allocated vectors need to be filled with the respective value.
                                           tiles_lon_min, &
                                           tiles_lon_max, &    
                                           tiles_lat_min, &
                                           tiles_lat_max, &
                                           nc_tot,        &
                                           nr_tot,        &
                                           nc_tile)

  SELECT CASE(itopo_type)
    CASE(topo_aster)
      PRINT*,'edges of ASTER domain: ', aster_lon_min,' ', aster_lon_max,' ', aster_lat_min,' ',aster_lat_max
    END SELECT
! mes <

!mes >
  PRINT*,' itopo_type :', itopo_type

  SELECT CASE (itopo_type)
  CASE (topo_aster)
   IF (lon_geo (tg%ie,tg%je,tg%ke) > aster_lon_max .OR. lon_geo(1,1,1) < aster_lon_min) THEN
   PRINT*, 'ASTER min lon is: ', aster_lon_min, ' and ASTER max lon is: ', aster_lon_max
   CALL abort_extpar('The chosen longitude edges are not within the ASTER domain.')
   END IF
   IF (lat_geo(tg%ie,tg%je,tg%ke) > aster_lat_max .OR. lat_geo(1,1,1) < aster_lat_min) THEN
   PRINT*, 'ASTER min lat is: ', aster_lat_min, ' and ASTER max lat is: ', aster_lat_max
   CALL abort_extpar('The chosen latitude edges are not within the ASTER domain.')
   END IF
  END SELECT


!roa >
  namelist_oro_smooth = 'INPUT_OROSMOOTH'

  IF(lfilter_oro) CALL read_namelists_extpar_orosmooth(namelist_oro_smooth,      &
                                           lfilter_oro,          &
                                           ilow_pass_oro,        &
                                           numfilt_oro,          &
                                           eps_filter,           &
                                           ifill_valley,         &
                                           rfill_valley,         &
                                           ilow_pass_xso,        &
                                           numfilt_xso,          &
                                           lxso_first,           &
                                           rxso_mask)
!roa <

  IF (lradtopo .AND. (.NOT. lfilter_oro)) THEN
    PRINT *,' Warning *** lradtopo should not be used without orography filtering *** '
    PRINT *,'                            (consistency problem)                        ' 
  ENDIF

  CALL det_topo_tiles_grid(topo_tiles_grid)
  !HA debug
  DO k=1,ntiles
    print *,'topo files: ', TRIM(topo_files(k))
  ENDDO
    print *,'topo_tiles_grid(1): ', topo_tiles_grid(1)

  CALL det_topo_grid(topo_grid)
  !HA debug
  print *,'topo_grid: ', topo_grid


  ! allocate globe fields for target grid
  !--------------------------------------------------------------------------------------------------------
  CALL allocate_topo_target_fields(tg,nhori)

  ! allocate additional fields for icon grid
  !--------------------------------------------------------------------------------------------------------

  SELECT CASE(igrid_type)
    CASE(igrid_icon) ! ICON GRID
    ! allocate addtional target fields
    nvertex = icon_grid_region%nverts
    CALL  allocate_additional_hh_param(nvertex)
  END SELECT

  ! call the aggregation routine
  !--------------------------------------------------------------------------------------------------------
 IF (igrid_type == igrid_icon) THEN ! ICON GRID
   IF (lsso_param) THEN
     PRINT *,'CALL agg_topo_data_to_target_grid with SSO'
     CALL agg_topo_data_to_target_grid_icon(topo_tiles_grid,  &
       &                                topo_grid,        &
       &                                tg,               &
       &                                topo_files,       &
       &                                lsso_param,       &
       &                                lsubtract_mean_slope, &
       &                                lfilter_oro,      &
       &                                lfilter_topo,      &
       &                                ilow_pass_oro,    &
       &                                numfilt_oro,      &
       &                                eps_filter,       &
       &                                ifill_valley,     &
       &                                rfill_valley,     &
       &                                ilow_pass_xso,    &
       &                                numfilt_xso,      &
       &                                lxso_first,       &
       &                                rxso_mask,        & 
!roa<
       &                                hh_topo,         &
       &                                stdh_topo,       &
       &                                fr_land_topo,    &
       &                                z0_topo,          &
       &                                no_raw_data_pixel,&
       &                                theta_topo,&
       &                                aniso_topo,&
       &                                slope_topo        )
   ELSE
     PRINT *,'CALL agg_topo_data_to_target_grid without SSO'
     CALL agg_topo_data_to_target_grid_icon(topo_tiles_grid,  &
       &                                topo_grid,        &
       &                                tg,               &
       &                                topo_files,       &
       &                                lsso_param,       &
!roa>
       &                                lsubtract_mean_slope, &
       &                                lfilter_oro,      &
       &                                lfilter_topo,      &
       &                                ilow_pass_oro,    &
       &                                numfilt_oro,      &
       &                                eps_filter,       &
       &                                ifill_valley,     &
       &                                rfill_valley,     &
       &                                ilow_pass_xso,    &
       &                                numfilt_xso,      &
       &                                lxso_first,       &
       &                                rxso_mask,        & 
!roa<
       &                                hh_topo,         &
       &                                stdh_topo,       &
       &                                fr_land_topo,    &
       &                                z0_topo,          &
       &                                no_raw_data_pixel )
 
       ENDIF

    ELSE  ! COSMO/GME GRID
   IF (lsso_param) THEN
     IF (lscale_separation) THEN
       PRINT *,'CALL agg_topo_data_to_target_grid_cosmo with SSO'
       CALL agg_topo_data_to_target_grid_cosmo(topo_tiles_grid, &
       &                                topo_grid,        &
       &                                tg,               &
       &                                topo_files,       &
       &                                lsso_param,       &
!< *mes                                                                                                                                                                                                                                   
       &                                lscale_separation,&
!> *mes                                                                                                                                                                                                                                   
!roa>                                                                                                                                                                                                                                     
       &                                lfilter_oro,      &
       &                                ilow_pass_oro,    &
       &                                numfilt_oro,      &
       &                                eps_filter,       &
       &                                ifill_valley,     &
       &                                rfill_valley,     &
       &                                ilow_pass_xso,    &
       &                                numfilt_xso,      &
       &                                lxso_first,       &
       &                                rxso_mask,        &
!roa<                                                                                                                                                                                                                                     
       &                                hh_topo,          &
       &                                stdh_topo,        &
       &                                fr_land_topo,     &
       &                                z0_topo,          &
       &                                no_raw_data_pixel,&
       &                                theta_topo,       &
       &                                aniso_topo,       &
       &                                slope_topo,       &
!< *mes                                                                                                                                                                                                                                   
       &                                raw_data_orography_path=raw_data_orography_path,& !_br 17.09.14                                                                                                                                   
       &                                raw_data_scale_sep_orography_path=raw_data_scale_sep_orography_path,& !_br 17.09.14                                                                                                               
       &                                scale_sep_files = scale_sep_files)
!> *mes                                                                                                                                                                                                                                   
     ELSE
       CALL agg_topo_data_to_target_grid_cosmo(topo_tiles_grid, &
       &                                topo_grid,        &
       &                                tg,               &
       &                                topo_files,       &
       &                                lsso_param,       &
       &                                lscale_separation,&
!roa>                                                                                                                                                                                                                                     
       &                                lfilter_oro,      &
       &                                ilow_pass_oro,    &
       &                                numfilt_oro,      &
       &                                eps_filter,       &
       &                                ifill_valley,     &
       &                                rfill_valley,     &
       &                                ilow_pass_xso,    &
       &                                numfilt_xso,      &
       &                                lxso_first,       &
       &                                rxso_mask,        &
!roa<                                                                                                                                                                                                                                     
       &                                hh_topo,          &
       &                                stdh_topo,        &
       &                                fr_land_topo,     &
       &                                z0_topo,          &
       &                                no_raw_data_pixel,&
       &                                theta_topo,       &
       &                                aniso_topo,       &
       &                                slope_topo,       &
       &                                raw_data_orography_path=raw_data_orography_path) !_br 17.09.14)                                                                                                                                   
     ENDIF
   ELSE
     IF (lscale_separation) THEN
       CALL agg_topo_data_to_target_grid_cosmo(topo_tiles_grid, &
       &                                topo_grid,        &
       &                                tg,               &
       &                                topo_files,       &
       &                                lsso_param,       &
!< *mes                                                                                                                                                                                                                                   
       &                                lscale_separation,&
!> *mes                                                                                                                                                                                                                                   
!roa>                                                                                                                                                                                                                                     
       &                                lfilter_oro,      &
       &                                ilow_pass_oro,    &
       &                                numfilt_oro,      &
       &                                eps_filter,       &
       &                                ifill_valley,     &
       &                                rfill_valley,     &
       &                                ilow_pass_xso,    &
       &                                numfilt_xso,      &
       &                                lxso_first,       &
       &                                rxso_mask,        &
!roa<                                                                                                                                                                                                                                     
       &                                hh_topo,          &
       &                                stdh_topo,        &
       &                                fr_land_topo,     &
       &                                z0_topo,          &
       &                                no_raw_data_pixel,&
       &                                raw_data_orography_path=raw_data_orography_path,& !_br 17.09.14                                                                                                                                   
       &                                raw_data_scale_sep_orography_path=raw_data_scale_sep_orography_path,& !_br 17.09.14                                                                                                               
       &                                scale_sep_files = scale_sep_files)
!                                                                                                                                                                                                                                         
     ELSE
       PRINT *,'CALL agg_topo_data_to_target_grid_cosmo without SSO'
       CALL agg_topo_data_to_target_grid_cosmo(topo_tiles_grid, &
       &                                topo_grid,        &
       &                                tg,               &
       &                                topo_files,       &
       &                                lsso_param,       &
       &                                lscale_separation,&
!roa>                                                                                                                                                                                                                                     
       &                                lfilter_oro,      &
       &                                ilow_pass_oro,    &
       &                                numfilt_oro,      &
       &                                eps_filter,       &
       &                                ifill_valley,     &
       &                                rfill_valley,     &
       &                                ilow_pass_xso,    &
       &                                numfilt_xso,      &
       &                                lxso_first,       &
       &                                rxso_mask,        &
!roa<                                                                                                                                                                                                                                     
       &                                hh_topo,         &
       &                                stdh_topo,       &
       &                                fr_land_topo,    &
       &                                z0_topo,          &
       &                                no_raw_data_pixel, &
       &                                raw_data_orography_path=raw_data_orography_path) !_br 17.09.14)                                                                                                                                   
     ENDIF

  ENDIF

END IF !igrid_type
   
   ! if the target domain has a higher resolution of than the GLOBE data set (30'') some grid elements might not
   ! be set by the routine agg_topo_data_to_target_grid, (no_raw_data_pixel(ie,je,ke) == 0 in this case
   ! loop overa all grid elements to check and perform a bilinear interplation if necessary
   k = 0
   undefined = -999.9


       PRINT*,'lsubtract_mean_slope is set to ',lsubtract_mean_slope


       PRINT *,'Maximum number of TOPO raw data pixel in a target grid element: '
       PRINT *,'MAXVAL(no_raw_data_pixel): ', MAXVAL(no_raw_data_pixel)

       PRINT *,'Minimal number of TOPO raw data pixel in a target grid element: '
       PRINT *,'MINVAL(no_raw_data_pixel): ', MINVAL(no_raw_data_pixel)

        PRINT *,'Maximum stdh_topo: '
       PRINT *,'MAXVAL(stdh_topo): ', MAXVAL(stdh_topo)

       PRINT *,'Minimal stdh_topo: '
       PRINT *,'MINVAL(stdh_topo): ', MINVAL(stdh_topo)

  ! consistency for small grid sizes, do not use estimates of variance for small sample size
!   IF ( (MAXVAL(no_raw_data_pixel)< 10).OR. (MINVAL(no_raw_data_pixel)==0)) THEN
   IF (MAXVAL(no_raw_data_pixel)< 10) THEN !_dl 22.8.14
     IF (lsso_param) THEN
       stdh_topo  = 0.  
       theta_topo = 0.
       aniso_topo = 0.
       slope_topo = 0.
     ENDIF
     z0_topo     = 0.
   ENDIF

   DO ke=1,tg%ke
   DO je=1,tg%je
   DO ie=1,tg%ie
      IF (no_raw_data_pixel(ie,je,ke) < 1 ) THEN
        IF (lsso_param) THEN
          stdh_topo(ie,je,ke)  = 0.  
          theta_topo(ie,je,ke) = 0.
          aniso_topo(ie,je,ke) = 0.
          slope_topo(ie,je,ke) = 0.
        ENDIF
         z0_topo(ie,je,ke)     = 0.
       ENDIF
   ENDDO
   ENDDO
   ENDDO


!   SELECT CASE(igrid_type)
!     CASE(igrid_icon) ! ICON GRID
!     nvertex = icon_domain_grid%nverts
!     je=1
!     ke=1
!     DO nv=1, nvertex
!       k=0
!       IF (vertex_param%npixel_vert(nv,je,ke) == 0) THEN ! interpolate from raw data in this case
!         k = k + 1 ! count number of grid element for which a bilinear interpolation is done
!          point_lon_geo =  rad2deg * icon_domain_grid%verts%vertex(nv)%lon 
!          point_lat_geo =  rad2deg * icon_domain_grid%verts%vertex(nv)%lat
!          CALL bilinear_interpol_topo_to_target_point(topo_files, topo_grid, &
!            &                                 topo_files,      &
!            &                                 undefined,    &
!            &                                 point_lon_geo,      &
!            &                                 point_lat_geo,      &
!            &                                 nrows,              &
!            &                                 h_band,    &
!            &                                 startrow_index,     &
!            &                                 endrow_index,       &
!            &                                 topo_target_value)
!
!          vertex_param%hh_vert(nv,je,ke) = topo_target_value
!       ENDIF  
!     ENDDO
!     print *,'bilinear interpolation of GLOBE data used for no of grid vertices: ',k
!   END SELECT

   ! compute the lradtopo parameters if needed
   IF ( lradtopo ) THEN
     CALL CPU_TIME(timestart)
     CALL compute_lradtopo(nhori,tg,hh_topo,slope_asp_topo,slope_ang_topo,horizon_topo,skyview_topo)
     CALL CPU_TIME(timeend)
      timediff = timeend - timestart
      PRINT *,'lradtopo calculations were done in: ', timediff
   ENDIF


   !HA debug
   PRINT *,'agg_topo_data_to_target_grid finished'

   ! output to netcdf file
   undefined = -999.9
   undef_int = -999

   netcdf_filename = TRIM(orography_buffer_file)
   print *, 'filename: ',TRIM(netcdf_filename)

   SELECT CASE(igrid_type)
     CASE(igrid_icon) ! ICON GRID
       IF (lsso_param) THEN
         CALL write_netcdf_buffer_topo(netcdf_filename,  &
          &                                tg,            &
          &                                undefined,     &
          &                                undef_int,     &
          &                                igrid_type,    &
          &                                lon_geo,       &
          &                                lat_geo,       &
          &                                fr_land_topo, &
          &                                hh_topo,      &
          &                                stdh_topo,    &
          &                                z0_topo,       &
          &                                lradtopo,      &
          &                                nhori,         &
          &                                theta_topo=theta_topo,&
          &                                aniso_topo=aniso_topo,&
          &                                slope_topo=slope_topo,&
          &                                vertex_param=vertex_param)
       ELSE
         CALL write_netcdf_buffer_topo(netcdf_filename,  &
          &                                tg,            &
          &                                undefined,     &
          &                                undef_int,     &
          &                                igrid_type,    &
          &                                lon_geo,       &
          &                                lat_geo,       &
          &                                fr_land_topo, &
          &                                hh_topo,      &
          &                                stdh_topo,    &
          &                                z0_topo,       &
          &                                lradtopo,      &
          &                                nhori,         &
          &                                vertex_param=vertex_param)
       ENDIF

     CASE DEFAULT

     IF (lradtopo) THEN
       IF (lsso_param) THEN
         CALL write_netcdf_buffer_topo(netcdf_filename,      &
         &                                     tg,            &
         &                                     undefined,     &
         &                                     undef_int,     &
         &                                     igrid_type,    &
         &                                     lon_geo,       &
         &                                     lat_geo,       &
         &                                     fr_land_topo, &
         &                                     hh_topo,      &
         &                                     stdh_topo,    &
         &                                     z0_topo,       &
         &                                     lradtopo,      &
         &                                     nhori,         &
         &                                     theta_topo=theta_topo,          &
         &                                     aniso_topo=aniso_topo,          &
         &                                     slope_topo=slope_topo,          &
         &                                     slope_asp_topo=slope_asp_topo, &
         &                                     slope_ang_topo= slope_ang_topo,&
         &                                     horizon_topo=horizon_topo,     &
         &                                     skyview_topo=skyview_topo)
       ELSE
         CALL write_netcdf_buffer_topo(netcdf_filename,      &
         &                                     tg,            &
         &                                     undefined,     &
         &                                     undef_int,     &
         &                                     igrid_type,    &
         &                                     lon_geo,       &
         &                                     lat_geo,       &
         &                                     fr_land_topo, &
         &                                     hh_topo,      &
         &                                     stdh_topo,    &
         &                                     z0_topo,       &
         &                                     lradtopo,      &
         &                                     nhori,         &
         &                                     slope_asp_topo=slope_asp_topo, &
         &                                     slope_ang_topo= slope_ang_topo,&
         &                                     horizon_topo=horizon_topo,     &
         &                                     skyview_topo=skyview_topo)
       ENDIF

     ELSE

       IF (lsso_param) THEN

         CALL write_netcdf_buffer_topo(netcdf_filename,        &
         &                                     tg,              &
         &                                     undefined,       &
         &                                     undef_int,       &
         &                                     igrid_type,      &
         &                                     lon_geo,         &
         &                                     lat_geo,         &
         &                                     fr_land_topo,   &
         &                                     hh_topo,        &
         &                                     stdh_topo,      &
         &                                     z0_topo,         &
         &                                     lradtopo,        &
         &                                     nhori,           &
         &                                     theta_topo=theta_topo,&
         &                                     aniso_topo=aniso_topo,&
         &                                     slope_topo=slope_topo)

       ELSE

         CALL write_netcdf_buffer_topo(netcdf_filename,        &
         &                                     tg,              &
         &                                     undefined,       &
         &                                     undef_int,       &
         &                                     igrid_type,      &
         &                                     lon_geo,         &
         &                                     lat_geo,         &
         &                                     fr_land_topo,   &
         &                                     hh_topo,        &
         &                                     stdh_topo,      &
         &                                     z0_topo,         &
         &                                     lradtopo,        &
         &                                     nhori)
       ENDIF
     ENDIF


   END SELECT


   SELECT CASE(igrid_type)
     CASE(igrid_icon) ! ICON GRID
        
       netcdf_filename = TRIM(orography_output_file)
       PRINT *,'write out ', TRIM(netcdf_filename)

       IF (lsso_param) THEN
         CALL write_netcdf_icon_grid_topo(netcdf_filename,&
          &                                icon_grid,      &
          &                                tg,             &
          &                                undefined,      &
          &                                undef_int,      &
          &                                lon_geo,        &
          &                                lat_geo,        &
          &                                fr_land_topo,  &
          &                                hh_topo,       &
          &                                stdh_topo,     &
          &                                z0_topo,        &
          &                                vertex_param,   &
          &                                theta_topo=theta_topo,&
          &                                aniso_topo=aniso_topo,&
          &                                slope_topo=slope_topo)
       ELSE
          CALL write_netcdf_icon_grid_topo(netcdf_filename,&
          &                                icon_grid,       &
          &                                tg,              &
          &                                undefined,       &
          &                                undef_int,       &
          &                                lon_geo,         &
          &                                lat_geo,         &
          &                                fr_land_topo,   &
          &                                hh_topo,        &
          &                                stdh_topo,      &
          &                                z0_topo,         &
          &                                vertex_param)
        ENDIF

     CASE(igrid_cosmo) ! COSMO grid

       netcdf_filename = TRIM(orography_output_file)
       PRINT *,'write out ', TRIM(netcdf_filename)

       IF(lradtopo) THEN
         IF (lsso_param) THEN
           CALL write_netcdf_cosmo_grid_topo(netcdf_filename,   &
             &                                     cosmo_grid,   &
             &                                     tg,           &
             &                                     undefined,    &
             &                                     undef_int,    &
             &                                     lon_geo,      &
             &                                     lat_geo,      &
             &                                     fr_land_topo,&
             &                                     hh_topo,     &
             &                                     stdh_topo,   &
             &                                     z0_topo,      &
             &                                     lradtopo,     &
             &                                     nhori,        &
             &                                     theta_topo=theta_topo,           &
             &                                     aniso_topo=aniso_topo,           &
             &                                     slope_topo=slope_topo,           &
             &                                     slope_asp_topo=slope_asp_topo,   &
             &                                     slope_ang_topo=slope_ang_topo,   &
             &                                     horizon_topo=horizon_topo,       &
             &                                     skyview_topo=skyview_topo)
         ELSE
           CALL write_netcdf_cosmo_grid_topo(netcdf_filename,   &
             &                                     cosmo_grid,   &
             &                                     tg,           &
             &                                     undefined,    &
             &                                     undef_int,    &
             &                                     lon_geo,      &
             &                                     lat_geo,      &
             &                                     fr_land_topo,&
             &                                     hh_topo,     &
             &                                     stdh_topo,   &
             &                                     z0_topo,      &
             &                                     lradtopo,     &
             &                                     nhori,        &
             &                                     slope_ang_topo=slope_ang_topo,   &
             &                                     horizon_topo=horizon_topo,       &
             &                                     skyview_topo=skyview_topo)
         ENDIF

       ELSE
         IF (lsso_param) THEN
           CALL write_netcdf_cosmo_grid_topo(netcdf_filename,&
             &                                cosmo_grid,     &
             &                                tg,             &
             &                                undefined,      &
             &                                undef_int,      &
             &                                lon_geo,        &
             &                                lat_geo,        &
             &                                fr_land_topo,  &
             &                                hh_topo,       &
             &                                stdh_topo,     &
             &                                z0_topo,        &
             &                                lradtopo,       &
             &                                nhori,          &
             &                                theta_topo=theta_topo,&
             &                                aniso_topo=aniso_topo,&
             &                                slope_topo=slope_topo)
         ELSE
           CALL write_netcdf_cosmo_grid_topo(netcdf_filename,&
             &                                cosmo_grid,     &
             &                                tg,             &
             &                                undefined,      &
             &                                undef_int,      &
             &                                lon_geo,        &
             &                                lat_geo,        &
             &                                fr_land_topo,  &
             &                                hh_topo,       &
             &                                stdh_topo,     &
             &                                z0_topo,        &
             &                                lradtopo,       &
             &                                nhori)
         ENDIF
       ENDIF


     CASE(igrid_gme) ! GME grid   
   END SELECT

   CALL deallocate_topo_fields()

  PRINT *,'============= topo_to_buffer done ==============='


END PROGRAM extpar_topo_to_buffer

