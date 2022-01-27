!+ Fortran module to aggregate GLOBE orogrphy data to the target grid
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_0         2010/12/21 Hermann Asensio
!  Initial releases
! V1_1         2011/01/20 Hermann Asensio
!  correct calculation of standard deviation of height
! V1_2         2011/03/25 Hermann Asensio
!  update to support ICON refinement grids
!  bug fix in interpolation routine and handling of undefined GLOBE values
! V1_4         2011/04/21 Anne Roches
!  implementation of orography smoothing
! V1_7         2013/01/25 Guenther Zaengl
!   Parallel threads for ICON and COSMO using Open-MP,
!   Several bug fixes and optimizations for ICON search algorithm,
!   particularly for the special case of non-contiguous domains;
!   simplified namelist control for ICON
!   Potential bugfix in treatment of GLOBE data for COSMO -
!   no impact on results detected so far
! V2_0        2013/04/09 Martina Messmer
!  introduce ASTER topography for external parameters
!  Change all 'globe' to topo in globe_files, remove all 'globe' in
!  change mo_GLOBE_data to mo_topo_data, globe_tiles_grid to
!  topo_tiles_grid, globe_files to topo_files, globe_grid to
!  topo_grid and change ntiles_gl to ntiles to obtain a more
!  dynamical code.
! V1_14        2014-07-18 Juergen Helmert
!  Combined COSMO Release
! V2_1         2015-01-12 Juergen Helmert
!  Bugfix correction covers CSCS SVN r5907-r6359
! V2_2         2015-02-23 Juergen Helmert
!  Small change for Intel compiler
! V2_5         2015-08-04 Juergen Helmert
!  Bugfix for index triple  j_n, j_c, j_s if mlat > 1, ICON only (Th. Raddatz)
! V2_6         2016-10-07 Juergen Helmert
!  Bugfix for orography block loop
!  Filtering of orography and different treatment
!  of SSO parameters (D. Reinert, G. Zaengl, J. Helmert)
! V2_8         2017-08-30 G. Zaengl
!  Revision of SSO parameter determination for ICON
! V2_9         2017-11-13 J. Helmert
!  Introduce Namelist switch lsubtract_mean_slope (= .FALSE.) for SSO
! V2_10        2018-02-19 Juergen Helmert
!  lsubtract_mean_slope, ERA-I surface temp for land points
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!> Fortran module to aggregate GLOBE orogrphy data to the target grid
!> \author Hermann Asensio
MODULE mo_agg_topo_icon

  USE mo_logging
  USE mo_kind,                  ONLY: wp,i4

  USE mo_io_units,              ONLY: filename_max

  USE mo_grid_structures,       ONLY: reg_lonlat_grid, &
       &                              rotated_lonlat_grid, &
       &                              icosahedral_triangular_grid, &
       &                              target_grid_def, &
       &                              igrid_icon

  USE mo_search_ll_grid,        ONLY: find_reg_lonlat_grid_element_index, &
       &                              find_rotated_lonlat_grid_element_index

  USE mo_icon_domain,           ONLY: icon_domain, grid_cells

  USE mo_topo_data,             ONLY: ntiles,         & !< there are 16/240 GLOBE/ASTER/MERIT tiles
       &                              max_tiles,      &
       &                              nc_tot,         & !< number of total GLOBE/ASTER columns un a latitude circle
       &                              nr_tot,         & !< total number of rows in GLOBE/ASTER/MERIT data
       &                              get_fill_value, &
       &                              itopo_type,     &
       &                              topo_gl,        &
       &                              topo_aster,     &
       &                              topo_merit,     &
       &                              get_varname

  USE mo_topo_sso,              ONLY: auxiliary_sso_parameter_icon,&
       &                              calculate_sso

  USE mo_topo_routines,         ONLY: open_netcdf_TOPO_tile, &
       &                              close_netcdf_TOPO_tile, &
       &                              get_topo_data_block, &
       &                              det_band_gd

  USE mo_target_grid_data,      ONLY: lon_geo, & !< longitude coordinates of the grid in the geographical system
       &                              lat_geo, &    !< latitude coordinates of the grid in the geographical system
       &                              search_res ! resolution of ICON grid search index list

  USE mo_icon_grid_data,        ONLY: icon_grid, & !< structure which contains the definition of the ICON grid
       &                              icon_grid_region

  USE mo_topo_tg_fields,        ONLY: vertex_param          !< this structure contains the fields
  USE mo_search_icongrid,       ONLY: walk_to_nc, find_nearest_vert

  USE mo_base_geometry,         ONLY: geographical_coordinates, &
       &                              cartesian_coordinates

  USE mo_additional_geometry,   ONLY: gc2cc

  USE mo_math_constants,        ONLY: rad2deg, deg2rad
  USE mo_physical_constants,    ONLY: re !< av. radius of the earth [m]

  USE mo_oro_filter,            ONLY: do_orosmooth

  USE mo_bilinterpol,           ONLY: get_4_surrounding_raw_data_indices, &
       &                              calc_weight_bilinear_interpol, &
       &                              calc_value_bilinear_interpol

 USE mo_io_utilities,           ONLY: join_path

  IMPLICIT NONE

  PUBLIC :: agg_topo_data_to_target_grid_icon

CONTAINS

  SUBROUTINE agg_topo_data_to_target_grid_icon(topo_tiles_grid,  &
       &                                       topo_grid,        &
       &                                       tg,               &
       &                                       topo_files,       &
       &                                       lsso_param,       &
       &                                       lsubtract_mean_slope, &
       &                                       lfilter_oro,      &
       &                                       ilow_pass_oro,    &
       &                                       numfilt_oro,      &
       &                                       eps_filter,       &
       &                                       ifill_valley,     &
       &                                       rfill_valley,     &
       &                                       ilow_pass_xso,    &
       &                                       numfilt_xso,      &
       &                                       lxso_first,       &
       &                                       rxso_mask,        &
       &                                       hh_target,        &
       &                                       hh_target_max,    &
       &                                       hh_target_min,    &
       &                                       stdh_target,      &
       &                                       fr_land_topo,     &
       &                                       z0_topo,          &
       &                                       no_raw_data_pixel,&
       &                                       theta_target,     &
       &                                       aniso_target,     &
       &                                       slope_target,     &
       &                                       raw_data_orography_path)

    TYPE(target_grid_def), INTENT(IN)        :: tg              !< !< structure with target grid description

    CHARACTER (LEN=*), INTENT(IN)            :: topo_files(1:max_tiles), &  !< filenames globe/aster raw data
         &                                      raw_data_orography_path

    LOGICAL, INTENT(IN)                      :: lsso_param, &
         &                                      lxso_first

    REAL(KIND=wp),    INTENT(IN)             :: eps_filter, &               !< smoothing param 
         &                                      rfill_valley, &             !< mask for valley filling (threshold value)
         &                                      rxso_mask                !< mask for eXtra SmOothing (threshold value)
    INTEGER(KIND=i4), INTENT(IN)             :: ifill_valley, &             !< fill valleys before or after oro smoothing
         &                                      ilow_pass_oro, &            !< type of oro smoothing and
         &                                      numfilt_oro, &              !< number of applications of the filter
         &                                      ilow_pass_xso, &            !< type of oro eXtra SmOothing for steep
         &                                      numfilt_xso              !< number of applications of the eXtra filter

    LOGICAL, INTENT(INOUT)                   :: lfilter_oro, &  !< oro smoothing to be performed? (TRUE/FALSE)
         &                                      lsubtract_mean_slope
  
    REAL(KIND=wp), INTENT(OUT)               :: hh_target(1:tg%ie,1:tg%je,1:tg%ke), &
         &                                      hh_target_max(1:tg%ie,1:tg%je,1:tg%ke), &
         &                                      hh_target_min(1:tg%ie,1:tg%je,1:tg%ke), &
         &                                      stdh_target(1:tg%ie,1:tg%je,1:tg%ke), &
         &                                      z0_topo(1:tg%ie,1:tg%je,1:tg%ke), &
         &                                      fr_land_topo(1:tg%ie,1:tg%je,1:tg%ke), &
         &                                      theta_target(1:tg%ie,1:tg%je,1:tg%ke), & !< sso parameter, angle of principal axis
         &                                      aniso_target(1:tg%ie,1:tg%je,1:tg%ke), & !< sso parameter, anisotropie factor
         &                                      slope_target(1:tg%ie,1:tg%je,1:tg%ke) !< sso parameter, mean slope


    INTEGER (KIND=i4), INTENT(OUT)           :: no_raw_data_pixel(1:tg%ie,1:tg%je,1:tg%ke)

    ! local variables
    TYPE(reg_lonlat_grid)                    :: topo_tiles_grid(1:ntiles), &!< structure w/ def of the raw data grid
         &                                      ta_grid, &
         &                                      topo_grid !< structure with defenition of the raw data grid for the whole GLOBE/ASTER dataset

    TYPE(geographical_coordinates)           :: target_geo_co  !< structure for geographical coordinates of raw data pixel

    TYPE(cartesian_coordinates)              :: target_cc_co

    REAL (KIND=wp)                           :: lon_topo(1:nc_tot), &   !< longitude coordinates of the GLOBE grid
         &                                      lon_red(nc_tot), lon_diff, lontopo, &
         &                                      lat_topo(1:nr_tot), &   !< latititude coordinates of the GLOBE grid
         &                                      hh_red(0:nc_tot+1,1:3), &  !< topographic height on reduced grid
         &                                      dhdxdx(1:nc_tot),dhdx(1:nc_tot),dhdy(1:nc_tot), &  !< x-gradient square for one latitude row
         &                                      dhdydy(1:nc_tot), &  !< y-gradient square for one latitude row
         &                                      dhdxdy(1:nc_tot), &  !< dxdy for one latitude row
         &                                      hh1_target(1:tg%ie,1:tg%je,1:tg%ke), &  !< mean height of grid element
         &                                      hh2_target(1:tg%ie,1:tg%je,1:tg%ke), &  !< square mean height of grid element
         &                                      hsmooth(1:tg%ie,1:tg%je,1:tg%ke), &  !< mean smoothed height of grid element
         &                                      h11(1:tg%ie,1:tg%je,1:tg%ke), & !< help variables
         &                                      h12(1:tg%ie,1:tg%je,1:tg%ke), & !< help variables
         &                                      h22(1:tg%ie,1:tg%je,1:tg%ke), & !< help variables
         &                                      hx(1:tg%ie,1:tg%je,1:tg%ke),hy(1:tg%ie,1:tg%je,1:tg%ke), &
         &                                      dx, dy, dx0, &    !  grid distance for gradient calculation (in [m])
         &                                      d2x, d2y, &       ! 2 times grid distance for gradient calculation (in [m])
         &                                      row_lat(1:3), &    ! latitude of the row for the topographic height array hh
         &                                      znorm, znfi2sum, zarg, & ! help variables for the estiamtion of the variance
         &                                      point_lon_geo, &       !< longitude coordinate in geographical system of input point
         &                                      point_lat_geo, &       !< latitude coordinate in geographical system of input point
         &                                      point_lon, point_lat, &
         &                                      topo_target_value, &  !< interpolated altitude from GLOBE data
         &                                      fr_land_pixel, &  !< interpolated fr_land from GLOBE data
         &                                      dnorm, &  !< scale factor
         &                                      zlnorm = 2250. , &    !< scale factor [m]
         &                                      alpha  = 1.E-05, & !< scale factor [1/m]
         &                                      factor, & !< factor
         &                                      zhp = 10.0_wp, &    !< height of Prandtl-layer [m]
         &                                      z0_topography, &   !< rougness length according to Erdmann Heise Formula
         &                                      dxrat, &                                       ! ratio of dy to dx when given in m
         &                                      dlon0, coslat, &
         &                                      wgt, wgtsum, &
         &                                      a2, a3, b   ! temporary coeficients
       
    REAL(KIND=wp)                            :: topo_rawdata(3)
    REAL(KIND=wp), ALLOCATABLE               :: sum_topo(:,:,:,:), sum_topo_sq(:,:,:,:), sum_topo_x(:,:,:,:)

    INTEGER(KIND=i4)                         :: nc_tot_p1, nc_red, ijlist(nc_tot), &
         &                                      ncids_topo(1:ntiles), &
         &                                      h_parallel(1:nc_tot), &  !< one line with GLOBE/ASTER data
         &                                      hh(0:nc_tot+1,1:3), & !< topographic height for gradient calculations
         &                                      ndata(1:tg%ie,1:tg%je,1:tg%ke), &  !< number of raw data pixel with land point
         &                                      undef_topo, &
         &                                      default_topo, &
         &                                      i, j, & ! counters
         &                                      ie, je, ke, &  ! indices for grid elements
         &                                      i_vert, j_vert, k_vert, & ! indeces for ICON grid vertices
         &                                      i1, i2, &
         &                                      nv, & ! counter
         &                                      j_n, j_c, j_s, & ! counter for northern, central and southern row
         &                                      j_new, & ! counter for swapping indices j_n, j_c, j_s
         &                                      block_row_start, &
         &                                      block_row, &
         &                                      errorcode, &
         &                                      start_cell_id, &
         &                                      nt, &
         &                                      ij, np, istart, iend, &
         &                                      mlat, & ! row number for GLOBE data
         &                                      nr_tot_fraction, &
         &                                      workload_fraction_index(10), &
         &                                      load_idx

    INTEGER(KIND=i4), ALLOCATABLE            :: ie_vec(:), iev_vec(:), &  ! indices for target grid elements
         &                                      h_block(:,:) !< a block of GLOBE/ASTER altitude data

    ! Some stuff for OpenMP parallelization
    INTEGER(KIND=i4)                         :: num_blocks, ib, blk_len, nlon_sub, &
         &                                      istartlon, iendlon, ishift, il
    !$ INTEGER :: omp_get_max_threads, omp_get_thread_num, thread_id
    !$ INTEGER (i4), ALLOCATABLE :: start_cell_arr(:)

    ! global data flag
    LOGICAL                                  :: lskip
                                             
    CHARACTER (LEN=80)                       :: varname_topo  !< name of variable for topo data
    CHARACTER(len=filename_max)              :: topo_file_1

    CALL logging%info('Enter routine: agg_topo_data_to_target_grid_icon')

    WRITE(message_text,*) 'lsso_param: ', lsso_param
    CALL logging%info(message_text)

    nc_red = 0
    nc_tot_p1 = nc_tot + 1
    topo_file_1 = topo_files(1)

    nc_tot_p1 = nc_tot + 1
    topo_file_1 = join_path(raw_data_orography_path,topo_files(1))

    ke = 1
    j_n = 1 ! index for northern row
    j_c = 2 ! index for central row
    j_s = 3 ! index for southern row

    CALL get_fill_value(topo_file_1,undef_topo)
    default_topo = 0

    SELECT CASE(itopo_type)
      CASE(topo_aster)
        hh = default_topo
        hh_red = default_topo
      CASE(topo_gl)
        hh = undef_topo
        hh_red = undef_topo
      CASE(topo_merit)
        hh = default_topo
        hh_red = default_topo
    END SELECT

    ! initialize some variables
    no_raw_data_pixel     = 0
    ndata      = 0
    z0_topo    = 0.0_wp
    hh_target   = 0.0_wp
    hh_target_max = -1.0e36_wp
    hh_target_min =  1.0e36_wp
    hh1_target  = 0.0_wp
    hh2_target  = 0.0_wp
    hx = 0.0_wp
    hy = 0.0_wp
    stdh_target = 0.0_wp
    IF (lsso_param) THEN
      theta_target = 0.0_wp
      aniso_target = 0.0_wp
      slope_target = 0.0_wp
    ENDIF

    IF (lsubtract_mean_slope) THEN
      ALLOCATE(sum_topo(3,tg%ie,tg%je,tg%ke), sum_topo_sq(3,tg%ie,tg%je,tg%ke), sum_topo_x(3,tg%ie,tg%je,tg%ke))
      sum_topo = 0.0_wp
      sum_topo_sq = 0.0_wp
      sum_topo_x = 0.0_wp
      topo_rawdata = 0.0_wp
    ENDIF

    h11         = 0.0_wp
    h12         = 0.0_wp
    h22         = 0.0_wp

    hsmooth     = 0.0_wp

    vertex_param%hh_vert = 0.0_wp
    vertex_param%npixel_vert = 0

    ! calculate the longitude coordinate of the GLOBE columns
    DO i =1, nc_tot
      lon_topo(i) = topo_grid%start_lon_reg + (i-1) * topo_grid%dlon_reg
    ENDDO

    ! calculate the latitiude coordinate of the GLOBE columns
    DO j = 1, nr_tot
      lat_topo(j) = topo_grid%start_lat_reg + (j-1) * topo_grid%dlat_reg
    ENDDO

    ALLOCATE(ie_vec(nc_tot),iev_vec(nc_tot))
    ie_vec(:) = 0
    iev_vec(:) = 0
    start_cell_id = 1

    nt = 1
    dx0 =  topo_tiles_grid(nt)%dlon_reg * deg2rad * re ! longitudinal distance between to topo grid elemtens at equator
    dy = topo_tiles_grid(nt)%dlat_reg * deg2rad * re
    d2y = 2.0_wp * dy

    ! first open the GLOBE netcdf files
    DO nt=1,ntiles
      CALL open_netcdf_TOPO_tile(join_path(raw_data_orography_path,topo_files(nt)), ncids_topo(nt))
    ENDDO
    mlat = 1
    block_row_start = mlat

    CALL det_band_gd(topo_grid,block_row_start, ta_grid)

    CALL get_varname(topo_file_1,varname_topo)

    IF(ALLOCATED(h_block)) THEN
      DEALLOCATE(h_block, stat=errorcode)
      IF(errorcode/=0) CALL logging%error('cant deallocate the h_block',__FILE__,__LINE__)
    ENDIF
    ALLOCATE (h_block(1:ta_grid%nlon_reg,1:ta_grid%nlat_reg), stat=errorcode)
    IF(errorcode/=0) CALL logging%error('cant allocate h_block',__FILE__,__LINE__)

    CALL get_topo_data_block(topo_file_1,     &
         &                   ta_grid,         &
         &                   topo_tiles_grid, &
         &                   ncids_topo,      &
         &                   h_block,         &
         &                   .FALSE.)


    block_row = 1

    ! determine start and end longitude of search
    istartlon = 1_i4
    iendlon = nc_tot
    DO i = 1, nc_tot
      point_lon = lon_topo(i)
      IF (point_lon < tg%minlon) istartlon = i + 1_i4
      IF (point_lon > tg%maxlon) THEN
        iendlon = i - 1_i4
        EXIT
      ENDIF
    ENDDO

    nlon_sub = iendlon - istartlon + 1

    num_blocks = 1
    !$ num_blocks = omp_get_max_threads()
    IF (MOD(nlon_sub,num_blocks)== 0) THEN
      blk_len = nlon_sub/num_blocks
    ELSE
      blk_len = nlon_sub/num_blocks + 1
    ENDIF
    !$ allocate(start_cell_arr(num_blocks))
    !$ start_cell_arr(:) = 1
    WRITE(message_text,*) 'nlon_sub: ',nlon_sub,' num_blocks: ',num_blocks, ' blk_len: ',blk_len
    CALL logging%info(message_text)

    nr_tot_fraction = nr_tot / 10
    load_idx = 1
    DO i = 1, 10
      workload_fraction_index(i) = nr_tot_fraction * i
    ENDDO

    CALL logging%info('Start loop over topo rows...')
    !-----------------------------------------------------------------------------
    topo_rows: DO mlat=1,nr_tot    !mes ><
      !   topo_rows: do mlat=1,20
      !-----------------------------------------------------------------------------
      !-----------------------------------------------------------------------------

      IF (mlat == workload_fraction_index(load_idx) ) THEN
        WRITE(message_text,'(I4,A)') load_idx * 10, '%'
        CALL logging%info(message_text)
        load_idx = load_idx + 1
      ENDIF

      block_row= block_row + 1
      IF(block_row > ta_grid%nlat_reg) THEN ! read in new block
        block_row_start = mlat +1
        block_row = 1
        CALL det_band_gd(topo_grid,block_row_start, ta_grid)
        lskip = .FALSE.
        IF (mlat > 1) THEN
          IF (ta_grid%end_lat_reg > tg%maxlat .OR. ta_grid%start_lat_reg < tg%minlat) lskip = .TRUE.
        ENDIF
        IF (.NOT. lskip) THEN
          IF(ALLOCATED(h_block)) THEN
            DEALLOCATE(h_block, stat=errorcode)
            IF(errorcode/=0) CALL logging%error('cant deallocate the h_block',__FILE__,__LINE__)
          ENDIF
          ALLOCATE (h_block(1:ta_grid%nlon_reg,1:ta_grid%nlat_reg), stat=errorcode)
          IF(errorcode/=0) CALL logging%error('cant allocate h_block',__FILE__,__LINE__)

          CALL get_topo_data_block(varname_topo,     &            !mes ><
               &                   ta_grid,         &
               &                   topo_tiles_grid, &
               &                   ncids_topo,     &
               &                   h_block,        &
               &                   .TRUE.)
        ENDIF
      ENDIF

      IF (mlat==1) THEN  !first row of topo data

        h_parallel(1:nc_tot) = h_block(1:nc_tot,block_row-1)
        row_lat(j_c) = topo_grid%start_lat_reg + (mlat-1) * topo_grid%dlat_reg

        hh(1:nc_tot,j_c) = h_parallel(1:nc_tot)  ! put data to "central row"
        hh(0,j_c)        = h_parallel(nc_tot) ! western wrap at -180/180 degree longitude
        hh(nc_tot_p1,j_c) = h_parallel(1)      ! eastern wrap at -180/180 degree longitude

        ! set undefined values to 0 altitude (default)
        WHERE (hh(:,j_c) == undef_topo)
          hh(:,j_c) = default_topo
        END WHERE
      ENDIF

      IF (mlat==2) j_s = 1 ! compensate the loss of index 1 in j_n,j_c,j_s by the handling of the case mlat==1
      row_lat(j_s) = topo_grid%start_lat_reg + mlat * topo_grid%dlat_reg  !  ((mlat+1)-1)

      lskip = .FALSE.

      IF (row_lat(j_s) > tg%maxlat .OR. row_lat(j_s) < tg%minlat) lskip = .TRUE.

      IF (lskip) THEN
        ! swap indices of the hh array for next data row before skipping the loop
        j_new = j_n ! the new data will be written in the former "northern" array
        j_n = j_c   ! the "center" row will become "northern" row
        j_c = j_s   ! the "southern" row will become "center" row
        j_s = j_new ! the new data will be written in the "southern" row
        CYCLE topo_rows
      ENDIF

      IF(mlat /= nr_tot) THEN !  read raw data south of "central" row except when you are at the most southern raw data line

        h_parallel(1:nc_tot) = h_block(1:nc_tot,block_row)
        hh(1:nc_tot,j_s) = h_parallel(1:nc_tot) ! put data to "southern row"
        hh(0,j_s)        = h_parallel(nc_tot) ! western wrap at -180/180 degree longitude
        hh(nc_tot_p1,j_s) = h_parallel(1)      ! eastern wrap at -180/180 degree longitude

        !dr test
        ! set undefined values to 0 altitude (default)
        WHERE (hh(:,j_s) == undef_topo)
          hh(:,j_s) = default_topo
        END WHERE

        ! compute hh_red
        dxrat = 1.0_wp/(COS(row_lat(j_c)*deg2rad))
        nc_red = NINT(REAL(nc_tot,wp)/dxrat)
        dxrat = REAL(nc_tot,wp)/REAL(nc_red,wp)
        np = INT((dxrat-1)/2.0_wp) +1
        dlon0 = ABS(topo_grid%dlat_reg)*dxrat
        ijlist(:) = 0
        DO i = 1, nc_red
          lon_red(i) = lon_topo(1)+(lon_topo(i)-lon_topo(1))*dxrat
          wgtsum = 0.0_wp
          hh_red(i,:) = 0.0_wp
          ij = INT(1+(i-1)*dxrat)
          ijlist(ij) = i
          istart = ij-np
          iend   = ij+np
          DO j = istart-1_i4,iend+1_i4
            ij = j
            IF (ij > nc_tot) THEN
              ij = ij - nc_tot
              lontopo = lon_topo(ij) + 360.0_wp
            ELSE IF (ij < 1) THEN
              ij = ij + nc_tot
              lontopo = lon_topo(ij) - 360.0_wp
            ELSE
              lontopo = lon_topo(ij)
            ENDIF
            lon_diff = ABS(lon_red(i)-lontopo)
            IF (lon_diff < 0.5_wp*dlon0) THEN
              wgt = MIN(5.0_wp,2.0_wp*dxrat,0.5_wp*dlon0/MAX(1.e-6_wp,lon_diff))
            ELSE IF (lon_diff-ABS(topo_grid%dlon_reg) < 0.5_wp*dlon0) THEN
              wgt = ABS(lon_diff-0.5_wp*dlon0)/ABS(topo_grid%dlon_reg)
            ELSE
              wgt = 0.0_wp
            ENDIF
            hh_red(i,1:3) = hh_red(i,1:3) + wgt*hh(ij,1:3)
            wgtsum = wgtsum + wgt
          ENDDO
          hh_red(i,1:3) = hh_red(i,1:3)/wgtsum
        ENDDO

        hh_red(0,1:3)        = hh_red(nc_red,1:3) ! western wrap at -180/180 degree longitude
        hh_red(nc_red+1,1:3) = hh_red(1, 1:3)      ! eastern wrap at -180/180 degree longitude

      ENDIF

      dx      = dx0 * COS(row_lat(j_c) * deg2rad)  ! longitudinal distance between to globe grid elemtens
      d2x = 2.0_wp * dx
      d2y = 2.0_wp * dy
      IF (mlat==1) THEN ! most northern row of raw data
        j_n = j_c  ! put the index of "northern row" to the same index as "central row"
        d2y = dy   ! adjust d2y in this case too
      ELSEIF (mlat==nr_tot) THEN ! most southern row of raw data
        j_s = j_c  ! put the index of "southern row" to the same index as "central row"
        d2y = dy   ! adjust d2y in this case too
      ENDIF

      IF(lsso_param) THEN
        CALL auxiliary_sso_parameter_icon(d2x,d2y,j_n,j_c,j_s,hh_red,INT(nc_red),dxrat,dhdx,dhdy,dhdxdx,dhdydy,dhdxdy)
      ENDIF
      ie_vec(istartlon:iendlon) = 0
      iev_vec(istartlon:iendlon) = 0

      point_lat = row_lat(j_c)

!$omp parallel do private(ib,il,ij,i,i1,i2,ishift,point_lon,thread_id,start_cell_id,target_geo_co,target_cc_co)
      DO ib = 1, num_blocks

        !$   thread_id = omp_get_thread_num()+1
        !$   start_cell_id = start_cell_arr(thread_id)
        ishift = NINT((istartlon-1)/dxrat)+(ib-1)*NINT(blk_len/dxrat)
        ij = NINT(blk_len/dxrat)
        IF (ib==num_blocks) THEN
          IF (tg%maxlon > 179.5_wp) THEN
            ij = nc_red
          ELSE
            ij = MIN(nc_red,NINT(blk_len/dxrat,i4))
          ENDIF
        ENDIF

        ! loop over one latitude circle of the raw data
        columns1: DO il = 1_i4, ij
          i = ishift+il
          IF (i >= nc_red) THEN
            CYCLE columns1
          ENDIF
          ! find the corresponding target grid indices
          point_lon = lon_red(i)

          ! reset start cell when entering a new row or when the previous data point was outside
          ! the model domain
          IF (il == 1 .OR. start_cell_id == 0) THEN
            i1 = NINT(point_lon*search_res)
            i2 = NINT(point_lat*search_res)
            start_cell_id = tg%search_index(i1,i2)
            IF (start_cell_id == 0) EXIT ! in this case, the whole row is empty; may happen with merged (non-contiguous) domains
          ENDIF

          target_geo_co%lon = point_lon * deg2rad ! note that the icon coordinates do not have the unit degree but radians
          target_geo_co%lat = point_lat * deg2rad
          target_cc_co = gc2cc(target_geo_co)
          CALL walk_to_nc(icon_grid_region,   &
               target_cc_co,     &
               start_cell_id,    &
               icon_grid%nvertex_per_cell, &
               icon_grid%nedges_per_vertex, &
               ie_vec(i))

          ! additional get the nearest vertex index for accumulating height values there
          IF (ie_vec(i) /= 0_i4) THEN
            CALL  find_nearest_vert(icon_grid_region, &
                 target_cc_co,                  &
                 ie_vec(i),        &
                 icon_grid%nvertex_per_cell,    &
                 iev_vec(i))
          ENDIF
        ENDDO columns1
        !$   start_cell_arr(thread_id) = start_cell_id
      ENDDO !num_blocks
!$omp end parallel do

      DO i=istartlon,iendlon

        ! call here the attribution of raw data pixel to target grid for different grid types
        IF (ijlist(i)/=0) THEN
          ie = ie_vec(ijlist(i))
          je = 1
          ke = 1
        ELSE
          CYCLE
        ENDIF

        ! get the nearest vertex index for accumulating height values there
        ! aggregate the vertex parameter here
        i_vert = iev_vec(ijlist(i))
        j_vert = 1
        k_vert = 1
        IF ((i_vert /=0)) THEN ! raw data pixel within target grid
          vertex_param%npixel_vert(i_vert,j_vert,k_vert) =  &
               vertex_param%npixel_vert(i_vert,j_vert,k_vert) + 1

          vertex_param%hh_vert(i_vert,j_vert,k_vert) =  &
               vertex_param%hh_vert(i_vert,j_vert,k_vert) +  hh_red(ijlist(i),j_c)
          !dr note that the following was equivalent to adding hh(i,j_s) except for mlat=1
          !dr is that correct. actually this part could be removed since it is no longer required
          !dr by icon anyway.
          !dr         vertex_param%hh_vert(i_vert,j_vert,k_vert) +  h_parallel(i)
        ENDIF

        IF ((ie /= 0).AND.(je/=0).AND.(ke/=0))THEN
          ! raw data pixel within target grid, see output of routine find_rotated_lonlat_grid_element_index
          no_raw_data_pixel(ie,je,ke) = no_raw_data_pixel(ie,je,ke) + 1

          !  summation of variables

          IF (hh_red(ijlist(i),j_c) /= default_topo) THEN
            ndata(ie,je,ke)      = ndata(ie,je,ke) + 1
            hh_target(ie,je,ke)  = hh_target(ie,je,ke) + hh_red(ijlist(i),j_c)

            IF (lsubtract_mean_slope) THEN
              topo_rawdata(1) = hh_red(ijlist(i),j_c)
              topo_rawdata(2) = lon_red(ijlist(i))
              IF (rad2deg*icon_grid_region%cells%center(ie)%lon - lon_red(ijlist(i)) > 180.0_wp) THEN
                topo_rawdata(2) = topo_rawdata(2) + 360.0_wp
              ELSE IF (rad2deg*icon_grid_region%cells%center(ie)%lon - lon_red(ijlist(i)) < -180.0_wp) THEN
                topo_rawdata(2) = topo_rawdata(2) - 360.0_wp
              ENDIF
              topo_rawdata(3) = row_lat(j_c)
              ! Compute intermediate sum, sum of squares and sum of cross products
              sum_topo(:,ie,je,ke) = sum_topo(:,ie,je,ke) + topo_rawdata(:)
              sum_topo_sq(:,ie,je,ke) = sum_topo_sq(:,ie,je,ke) + topo_rawdata(:) * topo_rawdata(:)
              sum_topo_x(1,ie,je,ke) = sum_topo_x(1,ie,je,ke) + topo_rawdata(2) * topo_rawdata(3)
              sum_topo_x(2,ie,je,ke) = sum_topo_x(2,ie,je,ke) + topo_rawdata(3) * topo_rawdata(1)
              sum_topo_x(3,ie,je,ke) = sum_topo_x(3,ie,je,ke) + topo_rawdata(1) * topo_rawdata(2)
            ELSE
              hh2_target(ie,je,ke) = hh2_target(ie,je,ke) + (hh_red(ijlist(i),j_c) * hh_red(ijlist(i),j_c))
            ENDIF

            hh_target_min(ie,je,ke) = MIN(hh_target_min(ie,je,ke), hh_red(ijlist(i),j_c))
            hh_target_max(ie,je,ke) = MAX(hh_target_max(ie,je,ke), hh_red(ijlist(i),j_c))

            IF(lsso_param) THEN
              h11(ie,je,ke)        = h11(ie,je,ke) + dhdxdx(ijlist(i))
              h12(ie,je,ke)        = h12(ie,je,ke) + dhdxdy(ijlist(i))
              h22(ie,je,ke)        = h22(ie,je,ke) + dhdydy(ijlist(i))
              hx(ie,je,ke)         = hx(ie,je,ke)  + dhdx(ijlist(i))
              hy(ie,je,ke)         = hy(ie,je,ke)  + dhdy(ijlist(i))
            ENDIF
          ENDIF
        ENDIF

      ENDDO ! loop over one latitude circle of the raw data


      ! swap indices of the hh array for next data row
      j_new = j_n ! the new data will be written in the former "northern" array
      j_n = j_c   ! the "center" row will become "northern" row
      j_c = j_s   ! the "southern" row will become "center" row
      j_s = j_new ! the new data will be written in the "southern" row

      !-----------------------------------------------------------------------------
      !-----------------------------------------------------------------------------
    ENDDO topo_rows
    !-----------------------------------------------------------------------------

    DEALLOCATE(ie_vec,iev_vec)
    !$     DEALLOCATE(start_cell_arr)

    CALL logging%info('...done')

    WRITE(message_text,*)'Maximum number of TOPO raw data pixel in a target grid element: '
    CALL logging%info(message_text)
    WRITE(message_text,*)'MAXVAL(no_raw_data_pixel): ', MAXVAL(no_raw_data_pixel)
    CALL logging%info(message_text)
    WRITE(message_text,*)'Index of target grid element: ', MAXLOC(no_raw_data_pixel)
    CALL logging%info(message_text)
    WRITE(message_text,*)'Maximum number of TOPO land pixel in a target grid element: '
    CALL logging%info(message_text)
    WRITE(message_text,*)'MAXVAL(ndata): ', MAXVAL(ndata)
    CALL logging%info(message_text)
    WRITE(message_text,*)   'Index of target grid element: ', MAXLOC(ndata)
    CALL logging%info(message_text)
    WRITE(message_text,*)   'Minimal number of TOPO raw data pixel in a target grid element: '
    CALL logging%info(message_text)
    WRITE(message_text,*) 'MINVAL(no_raw_data_pixel): ', MINVAL(no_raw_data_pixel)
    CALL logging%info(message_text)
    WRITE(message_text,*)   'Index of target grid element: ', MINLOC(no_raw_data_pixel)
    CALL logging%info(message_text)
    WRITE(message_text,*) 'Minimal number of TOPO land pixel in a target grid element: '
    CALL logging%info(message_text)
    WRITE(message_text,*) 'MINVAL(ndata): ', MINVAL(ndata)
    CALL logging%info(message_text)
    WRITE(message_text,*) 'Index of target grid element: ', MINLOC(ndata)
    CALL logging%info(message_text)

    hh1_target = hh_target ! save values of hh_target for computations of standard deviation

    ! Average height
    DO ke=1, tg%ke
      DO je=1, tg%je
        DO ie=1, tg%ie
          IF (no_raw_data_pixel(ie,je,ke) /= 0 .AND. ndata(ie,je,ke) /= 0)  THEN
            ! avoid division by zero for small target grids
            ! average height, oceans point counted as 0 height
            hh_target(ie,je,ke) = hh_target(ie,je,ke)/no_raw_data_pixel(ie,je,ke)
            hx(ie,je,ke) = hx(ie,je,ke)/ndata(ie,je,ke)
            hy(ie,je,ke) = hy(ie,je,ke)/ndata(ie,je,ke)
            ! fraction land
            fr_land_topo(ie,je,ke) = REAL(ndata(ie,je,ke),wp)/REAL(no_raw_data_pixel(ie,je,ke),wp) 
          ELSE
            hh_target(ie,je,ke) = REAL(default_topo)
            fr_land_topo(ie,je,ke) = 0.0_wp
          ENDIF

          ! set still not covered points to 0
          IF (hh_target_max(ie,je,ke) < -1.0e+35_wp ) hh_target_max(ie,je,ke) = 0.0_wp
          IF (hh_target_min(ie,je,ke) > 1.0e+35_wp) hh_target_min(ie,je,ke) = 0.0_wp

        ENDDO
      ENDDO
    ENDDO

    hsmooth = hh_target

    ! oro filt here
    lfilter_oro = .FALSE.
    IF (lfilter_oro)  CALL do_orosmooth(tg,                                 &
         &                                      hh_target,        &
         &                                      fr_land_topo,    &
         &                                      ilow_pass_oro,    &
         &                                      numfilt_oro,      &
         &                                      eps_filter,       &
         &                                      ifill_valley,     &
         &                                      rfill_valley,     &
         &                                      ilow_pass_xso,    &
         &                                      numfilt_xso,      &
         &                                      lxso_first,       &
         &                                      rxso_mask,        &
         &                                      hsmooth           )


    ! Average height for vertices
    DO ke=1, 1
      DO je=1, 1
        DO ie=1, icon_grid_region%nverts

          IF (vertex_param%npixel_vert(ie,je,ke) /= 0) THEN ! avoid division by zero for small target grids
            vertex_param%hh_vert(ie,je,ke) =  &
                 vertex_param%hh_vert(ie,je,ke)/vertex_param%npixel_vert(ie,je,ke) ! average height
          ELSE
            vertex_param%hh_vert(ie,je,ke) = REAL(default_topo)
          ENDIF
        ENDDO
      ENDDO
    ENDDO

    !     Standard deviation of height.
    DO ke=1, tg%ke
      DO je=1, tg%je
        DO ie=1, tg%ie
          ! estimation of variance
          IF (no_raw_data_pixel(ie,je,ke) > 1) THEN
            znorm = 1.0_wp/(REAL(no_raw_data_pixel(ie,je,ke),wp) * (REAL(no_raw_data_pixel(ie,je,ke),wp)-1))
          ELSE
            znorm = 0.0_wp
          ENDIF

          IF (lsubtract_mean_slope) THEN
            coslat = COS(icon_grid_region%cells%center(ie)%lat)
            a2 = 111111.0_wp * coslat * hx(ie,je,ke)
            a3 = 111111.0_wp * hy(ie,je,ke)
            b = 111111.0_wp * rad2deg * (  coslat * hx(ie,je,ke) * icon_grid_region%cells%center(ie)%lon   &
              &                          +          hy(ie,je,ke) * icon_grid_region%cells%center(ie)%lat )   &
              & - hh_target(ie,je,ke)
            
            hh2_target(ie,je,ke) =  sum_topo_sq(1,ie,je,ke) + ndata(ie,je,ke) * b*b  &
              &                    + a2*a2* sum_topo_sq(2,ie,je,ke) + a3*a3* sum_topo_sq(3,ie,je,ke)  &
              &                    + 2.0_wp*(a2*a3*sum_topo_x(1,ie,je,ke) - a3*sum_topo_x(2,ie,je,ke)  &
              &                    -  a2*sum_topo_x(3,ie,je,ke)   &
              &                    + b * (sum_topo(1,ie,je,ke) - a2*sum_topo(2,ie,je,ke) - a3*sum_topo(3,ie,je,ke)))
            
            znfi2sum = no_raw_data_pixel(ie,je,ke) * hh2_target(ie,je,ke)
            zarg     = znfi2sum * znorm
          ELSE
            znfi2sum = no_raw_data_pixel(ie,je,ke) * hh2_target(ie,je,ke)
            zarg     = ( znfi2sum - (hh1_target(ie,je,ke)*hh1_target(ie,je,ke))) * znorm
          ENDIF
          zarg = MAX(zarg,0.0_wp) ! truncation errors may cause zarg < 0.0
          stdh_target(ie,je,ke) = SQRT(zarg)
        ENDDO
      ENDDO
    ENDDO

    IF (lsubtract_mean_slope) THEN
      DEALLOCATE(sum_topo, sum_topo_sq, sum_topo_x)
    ENDIF

    IF (lsso_param) THEN

      CALL calculate_sso(tg,no_raw_data_pixel,    &
           &                   h11,h12,h22,stdh_target,&
           &                   theta_target,           &
           &                   aniso_target,           &
           &                   slope_target)

    ENDIF
    !----------------------------------------------------------------------------------
    ! calculate roughness length
    ! first zo_topo with "Erdmann Heise formula"
    !----------------------------------------------------------------------------------


    dnorm = 60000.0_wp         ! dummy value for normation of Erdmann Heise formula

    !---------------------------------------------------------------------------------
    ! Erdman Heise Formel
    !---------------------------------------------------------------------------------
    factor= alpha*ATAN(dnorm/zlnorm) !  alpha  = 1.E-05 [1/m] ,  zlnorm = 2250 [m]
    DO ke=1, tg%ke
      DO je=1, tg%je
        DO ie=1, tg%ie
          z0_topography = factor*stdh_target(ie,je,ke)**2.0_wp
          z0_topography = MIN(z0_topography,zhp-1.0_wp)
          z0_topo(ie,je,ke) = z0_topography
        ENDDO
      ENDDO
    ENDDO

    !roa>
    ! set the orography variable hh_target to the smoothed orography variable
    ! hsmooth in case of orogrpahy smoothing in extpar
    IF (lfilter_oro) THEN
      hh_target (:,:,:) = hsmooth (:,:,:)
    ENDIF


    ! bilinear interpolation of the orography in case of target grid points having
    ! no corresponding points in GLOBE

    CALL get_fill_value(topo_file_1,undef_topo)
    CALL get_varname(topo_file_1,varname_topo)

    DO ke=1, tg%ke
      DO je=1, tg%je
        DO ie=1, tg%ie
          IF (no_raw_data_pixel(ie,je,ke) == 0) THEN  ! bilinear interpolation to target grid

            point_lon_geo = lon_geo(ie,je,ke)
            point_lat_geo = lat_geo(ie,je,ke)

            CALL bilinear_interpol_topo_to_target_point(topo_grid,       &
                 &                                      topo_tiles_grid, &
                 &                                      ncids_topo,     &
                 &                                      lon_topo,       &
                 &                                      lat_topo,       &
                 &                                      point_lon_geo,   &
                 &                                      point_lat_geo,   &
                 &                                      fr_land_pixel,   &
                 &                                      topo_target_value, undef_topo, varname_topo)

            fr_land_topo(ie,je,ke) = fr_land_pixel
            hh_target(ie,je,ke) = topo_target_value

            IF (lsso_param) THEN
              theta_target(ie,je,ke) = 0.0_wp
              aniso_target(ie,je,ke) = 0.0_wp
              slope_target(ie,je,ke) = 0.0_wp
            ENDIF
            stdh_target(ie,je,ke)  = 0.0_wp
            z0_topo(ie,je,ke)      = 0.0_wp
          ENDIF
        ENDDO
      ENDDO
    ENDDO

    je=1
    ke=1

    DO nv=1, icon_grid_region%nverts
      IF (vertex_param%npixel_vert(nv,je,ke) == 0) THEN ! interpolate from raw data in this case
        point_lon_geo =  rad2deg * icon_grid_region%verts%vertex(nv)%lon
        point_lat_geo =  rad2deg * icon_grid_region%verts%vertex(nv)%lat

        CALL bilinear_interpol_topo_to_target_point(topo_grid,       &
             &                                      topo_tiles_grid, &
             &                                      ncids_topo,     &
             &                                      lon_topo,       &
             &                                      lat_topo,       &
             &                                      point_lon_geo,   &
             &                                      point_lat_geo,   &
             &                                      fr_land_pixel,   &
             &                                      topo_target_value, undef_topo, varname_topo)

        vertex_param%hh_vert(nv,je,ke) = topo_target_value
      ENDIF
    ENDDO

    ! close the GLOBE netcdf files
    DO nt=1,ntiles
      CALL close_netcdf_TOPO_tile(ncids_topo(nt))
    ENDDO

  CALL logging%info('Exit routine: agg_topo_data_to_target_grid_icon')

  END SUBROUTINE agg_topo_data_to_target_grid_icon

  !----------------------------------------------------------------------------------------------------------------

  !> subroutine for bilenar interpolation from GLOBE data (regular lonlat grid) to a single target point
  !!
  !! the GLOBE data are passed to the subroutine in the topo_data_block 2D-Array, which is re-read
  !! from the raw data file if the target point is out of the range of the data block.
  !! (If the data block is not too small, repeated I/O to the hard disk is avoided, reading from memory is much faster.)
  !!
  !! the definition of the regular lon-lat grid requires
  !! - the coordinates of the north-western point of the domain ("upper left") startlon_reg_lonlat and startlat_reg_lonlat
  !! - the increment dlon_reg_lonlat and dlat_reg_lonlat(implict assuming that the grid definiton goes
  !!   from the west to the east and from the north to the south)
  !! - the number of grid elements nlon_reg_lonlat and nlat_reg_lonlat for both directions
  SUBROUTINE bilinear_interpol_topo_to_target_point( topo_grid,      &
       topo_tiles_grid,&
       ncids_topo,    &
       lon_topo,      &
       lat_topo,      &
       point_lon_geo,  &
       point_lat_geo,  &
       fr_land_pixel,  &
       topo_target_value, undef_topo, varname_topo)


    TYPE(reg_lonlat_grid), INTENT(IN) :: topo_grid, & !< structure with defenition of the raw data grid for whole GLOBE/ASTER dataset
         &                               topo_tiles_grid(1:ntiles)!< structure w/ def of the raw data for 16/36 GLOBE/ASTER tiles

    INTEGER (KIND=i4), INTENT(IN)     :: ncids_topo(1:ntiles), & !< ncid for topo tiles
         &                               undef_topo

    REAL (KIND=wp), INTENT(IN)        :: lon_topo(1:nc_tot), &    !< longitude coordinates of the GLOBE grid
         &                               lat_topo(1:nr_tot), &    !< latititude coordinates of the GLOBE grid
         &                               point_lon_geo, &        !< longitude coordinate in geographical system of input point
         &                               point_lat_geo       !< latitude coordinate in geographical system of input point

    REAL (KIND=wp), INTENT(OUT)       :: fr_land_pixel, &   !< interpolated fr_land from GLOBE data
         &                               topo_target_value  !< interpolated altitude from GLOBE data

    CHARACTER (LEN=80),    INTENT(IN) :: varname_topo  !< name of variable

    ! local variables
    INTEGER (i4), ALLOCATABLE         :: h_block(:,:) !< a block of GLOBE altitude data

    TYPE(reg_lonlat_grid)             :: ta_grid

    INTEGER (KIND=i4)                 :: western_column, &      !< the index of the western_column of data to read in
         &                               eastern_column, &      !< the index of the eastern_column of data to read in
         &                               northern_row, &        !< the index of the northern_row of data to read in
         &                               southern_row, &        !< the index of the southern_row of data to read in
         &                               default_topo, & 
         &                               errorcode

    REAL (KIND=wp)                    :: bwlon, &   !< weight for bilinear interpolation
         &                               bwlat, &   !< weight for bilinear interpolation
         &                               topo_point_sw, &        !< value of the GLOBE raw data pixel south west
         &                               topo_point_se, &        !< value of the GLOBE raw data pixel south east
         &                               topo_point_ne, &        !< value of the GLOBE raw data pixel north east
         &                               topo_point_nw       !< value of the GLOBE raw data pixel north west

    LOGICAL                          :: gldata=.TRUE. ! GLOBE data are global

    default_topo = 0

    ! get four surrounding raw data indices
    CALL  get_4_surrounding_raw_data_indices(topo_grid,     &
         &                                      lon_topo,     &
         &                                      lat_topo,     &
         &                                      gldata,        &
         &                                      point_lon_geo, &
         &                                      point_lat_geo, &
         &                                      western_column,&
         &                                      eastern_column,&
         &                                      northern_row,  &
         &                                      southern_row)

    ta_grid%dlon_reg = topo_grid%dlon_reg
    ta_grid%dlat_reg = topo_grid%dlat_reg
    ta_grid%nlon_reg = eastern_column - western_column + 1
    ta_grid%nlat_reg = southern_row - northern_row + 1
    ta_grid%start_lon_reg = lon_topo(western_column)
    ta_grid%end_lon_reg = lon_topo(eastern_column)
    ta_grid%start_lat_reg = lat_topo(northern_row)
    ta_grid%end_lat_reg  = lat_topo(southern_row)

    ! calculate weight for bilinear interpolation
    CALL calc_weight_bilinear_interpol(point_lon_geo,             &
         &                                point_lat_geo,             &
         &                                lon_topo(western_column), &
         &                                lon_topo(eastern_column), &
         &                                lat_topo(northern_row),   &
         &                                lat_topo(southern_row),   &
         &                                bwlon,                     &
         &                                bwlat)

    ALLOCATE (h_block(western_column:eastern_column,northern_row:southern_row), STAT=errorcode)
    IF(errorcode/=0) CALL logging%error('Cant allocate h_block',__FILE__,__LINE__)
    CALL get_topo_data_block(varname_topo,     &   !mes ><
         &                       ta_grid,         &
         &                       topo_tiles_grid, &
         &                       ncids_topo,     &
         &                       h_block,        &
         &                       .TRUE.)
    ! check for undefined GLOBE data, which indicate ocean grid element

    IF( h_block(western_column,southern_row) == undef_topo) THEN
      topo_point_sw = 0.0_wp
      h_block(western_column,southern_row) = default_topo
    ELSE
      topo_point_sw = 1.0_wp
    ENDIF

    IF( h_block(eastern_column,southern_row) == undef_topo) THEN
      topo_point_se = 0.0_wp
      h_block(eastern_column,southern_row) = default_topo
    ELSE
      topo_point_se = 1.0_wp
    ENDIF

    IF( h_block(eastern_column,northern_row) == undef_topo) THEN
      topo_point_ne = 0.0_wp
      h_block(eastern_column,northern_row) = default_topo
    ELSE
      topo_point_ne = 1.0_wp
    ENDIF

    IF( h_block(western_column,northern_row) == undef_topo) THEN
      topo_point_nw = 0.0_wp
      h_block(western_column,northern_row) = default_topo
    ELSE
      topo_point_nw = 1.0_wp
    ENDIF


    ! perform the interpolation, first for fraction land
    fr_land_pixel = calc_value_bilinear_interpol(bwlon,          &
         &  bwlat,          &
         &  topo_point_sw, &
         &  topo_point_se, &
         &  topo_point_ne, &
         &  topo_point_nw)

    topo_point_sw = h_block(western_column,southern_row)
    topo_point_se = h_block(eastern_column,southern_row)
    topo_point_ne = h_block(eastern_column,northern_row)
    topo_point_nw = h_block(western_column,northern_row)

    ! perform the interpolation for height
    topo_target_value = calc_value_bilinear_interpol(bwlon,          &
         &       bwlat,          &
         &       topo_point_sw, &
         &       topo_point_se, &
         &       topo_point_ne, &
         &       topo_point_nw)

  END SUBROUTINE bilinear_interpol_topo_to_target_point

END MODULE mo_agg_topo_icon
