!+ Fortran module to aggregate globcover land use data to a target grid
!
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_3         2011/04/19 Hermann Asensio
!  Initial release
! V1_4         2011/04/21 Hermann Asensio
!  clean up
! V1_7         2013/01/25 Guenther Zaengl
!   Parallel threads for ICON and COSMO using Open-MP,
!   Several bug fixes and optimizations for ICON search algorithm,
!   particularly for the special case of non-contiguous domains;
!   simplified namelist control for ICON
! V2_0         2013/06/04 Martina Messmer
!   adaptations to the new Globecover 2009 tiles
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!> Fortran module to aggregate globcover land use data to a target grid
!!
!> \author Hermann Asensio
MODULE mo_agg_globcover

  USE mo_logging
  USE mo_kind,                  ONLY: wp, i4, i2
                                
  USE mo_grid_structures,       ONLY: reg_lonlat_grid,     &
       &                              igrid_icon,          &
       &                              igrid_cosmo,         &
       &                              target_grid_def   !< type definition of structure for tg       
                                
  USE mo_search_ll_grid,        ONLY: find_reg_lonlat_grid_element_index
                                
  USE mo_io_utilities,          ONLY: check_netcdf

  USE mo_search_target_grid,    ONLY: find_nearest_target_grid_element

  USE mo_globcover_data,        ONLY: globcover_grid, &
       &                              lon_globcover,  &
       &                              lat_globcover,  &
       &                              ntiles_globcover

  USE mo_globcover_lookup_tables, ONLY: name_lookup_table_globcover, &
       &                                init_globcover_lookup_tables, &
       &                                get_name_globcover_lookup_tables, get_globcover_idx, &
       &                                z0_lt_globcover, lnz0_lt_globcover, plc_mn_lt_globcover, &
       &                                plc_mx_lt_globcover, lai_mn_lt_globcover, &
       &                                lai_mx_lt_globcover, rd_lt_globcover, skinc_lt_globcover, &
       &                                emiss_lt_globcover, rs_min_lt_globcover, globcover_look_up, &
       &                                get_corinecover_idx, corinecover_look_up

  USE mo_landuse_routines,      ONLY: det_band_globcover_data, &
       &                              get_globcover_data_block

  USE mo_math_constants,        ONLY: deg2rad

  USE mo_physical_constants,    ONLY: re

  USE mo_target_grid_data,      ONLY: lon_geo, & !< longitude coordinates of the COSMO grid in the geographical system
       &                              lat_geo, & !< latitude coordinates of the COSMO grid in the geographical system
       &                              search_res !< resolution of ICON grid search index list

  USE netcdf,                   ONLY: nf90_open,      &
       &                              nf90_close,     &
       &                              nf90_inq_varid, &
       &                              nf90_get_var,   &
       &                              nf90_nowrite

#ifdef _OPENMP
  USE omp_lib
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: agg_globcover_data_to_target_grid

  REAL(wp), PARAMETER :: rs_min_undef=999. !< undefined value for minimal stomata resistance


  CONTAINS

  !> Subroutine to aggregate globcover data to the target grid
  SUBROUTINE agg_globcover_data_to_target_grid(globcover_file,          &
       &                                          ilookup_table_globcover, &
       &                                          l_use_corine,            &
       &                                          undefined,               &
       &                                          globcover_tiles_grid,    &
       &                                          tg,                      &
       &                                          nclass_globcover,        &
       &                                          globcover_class_fraction,&
       &                                          globcover_class_npixel,  &
       &                                          globcover_tot_npixel,    &
       &                                          fr_land_globcover ,      &
       &                                          ice_globcover,           &
       &                                          z0_globcover,            &
       &                                          root_globcover,          &
       &                                          plcov_mn_globcover,      &
       &                                          plcov_mx_globcover,      &
       &                                          lai_mn_globcover,        &
       &                                          lai_mx_globcover,        &
       &                                          rs_min_globcover,        &
       &                                          urban_globcover,         &
       &                                          for_d_globcover,         &
       &                                          for_e_globcover,         &
       &                                          skinc_globcover,         &
       &                                          emissivity_globcover)

    !-------------------------------------------------------------------------------------

    CHARACTER (LEN=*), INTENT(IN)            :: globcover_file(:)  !< filename globcover raw data
    REAL (KIND=wp), INTENT(IN)               :: undefined            !< undef value
    TYPE(reg_lonlat_grid), INTENT(IN)        :: globcover_tiles_grid(:)  ! grid structure of globcover tiles
    TYPE(target_grid_def), INTENT(IN)        :: tg  !< structure with target grid description

    INTEGER(KIND=i4), INTENT(IN)            :: nclass_globcover, &  !< globcover has 23 classes for the land use description
         &                                     ilookup_table_globcover

    LOGICAL, INTENT(IN)                     :: l_use_corine

    INTEGER (KIND=i4), INTENT(OUT)          :: globcover_class_npixel(:,:,:,:), & 
         &                                     globcover_tot_npixel(:,:,:)

    REAL (KIND=wp), INTENT(OUT)             :: globcover_class_fraction(:,:,:,:), & 
         &                                     fr_land_globcover(:,:,:), &  !< fraction land due to globcover raw data
         &                                     ice_globcover(:,:,:), &      !< fraction of ice due to globcover raw data
         &                                     z0_globcover(:,:,:), &       !< roughness length due to globcover land use data
         &                                     root_globcover(:,:,:), &     !< root depth due to globcover land use data
         &                                     plcov_mx_globcover(:,:,:), & !< plant cover maximum due to globcover land use data
         &                                     plcov_mn_globcover(:,:,:), & !< plant cover minimum due to globcover land use data
         &                                     lai_mx_globcover(:,:,:), &   !< Leaf Area Index maximum due to globcover land use data
         &                                     lai_mn_globcover(:,:,:), &   !< Leaf Area Index minimum due to globcover land use data
         &                                     rs_min_globcover(:,:,:), &   !< minimal stomata resistance due to globcover land use data
         &                                     urban_globcover(:,:,:), &    !< urban fraction due to globcover land use data
         &                                     for_d_globcover(:,:,:), &    !< deciduous forest (fraction) due to globcover land use data
         &                                     for_e_globcover(:,:,:), &    !< evergreen forest (fraction) due to globcover land use data
         &                                     skinc_globcover(:,:,:), &    !< skin conductivity due to globcover land use data
         &                                     emissivity_globcover(:,:,:) !< longwave emissivity due to globcover land use da

    !local variables
    ! structure with definition of the target area grid (dlon must be the same for the whole GLOBCOVER dataset)
    TYPE(reg_lonlat_grid)                   :: ta_grid

    INTEGER (KIND=i4)                       :: undefined_integer, &  ! undef value
         &                                     l, &       ! counters
         &                                     nt, &            ! counter
         &                                     i_col, j_row, &  ! counter
         &                                     i_lu, j_lu, & 
         &                                     ie, je, ke, &   ! indices for target grid elements
         &                                     lu, &   ! land use class
         &                                     start_cell_id, &  !< ID of starting cell for ICON search
         &                                     ii1, ii2, & 
         &                                     nclass, &  ! index in array of globcover tables
         &                                     ncid_globcover(1:ntiles_globcover), &             !< netcdf unit file number
         &                                     varid_globcover, &                !< id of variable
         &                                     varid_gc(1:ntiles_globcover), & 
         &                                     nlon, iendlon, & 
         &                                     block_row_start, & 
         &                                     block_row, & 
         &                                     mlat, & 
         &                                     tile, & 
         &                                      k_error, errorcode   ! error return code

    INTEGER (KIND=i4), ALLOCATABLE         :: ie_vec(:), je_vec(:), ke_vec(:)  ! indices for target grid elements

    INTEGER (KIND=i2)                      :: globcover_data_row(globcover_grid%nlon_reg), &
         &                                    globcover_data_pixel(1:1,1:1)

    INTEGER (KIND=i2), ALLOCATABLE         :: lu_block(:,:)  ! a block of GLOBCOVER land use data

    REAL (KIND=wp)                          :: default_real, & 
         &                                     a_class(1:tg%ie,1:tg%je,1:tg%ke,1:nclass_globcover), & 
         &                                     a_weight(1:tg%ie,1:tg%je,1:tg%ke), & 
         &                                     apix, &       !< area of a raw data pixel
         &                                     apix_e, &       !< area of a raw data pixel at equator
         &                                     point_lon, point_lat, & 
         &                                     pland, &           !< land cover                      (-)
         &                                     pice, &            !< ice fraction                    (-)
         &                                     plnz0, &           !< logarithm of roughness length   (m)
         &                                     proot, &           !< root depth                      (m)
         &                                     pmn, &             !< minimal plant cover             (-)
         &                                     pmx, &             !< maximum plant cover             (-)
         &                                     plaimn, &          !< minimal leaf area index         (m**2/m**2)
         &                                     plaimx, &          !< maximum leaf area index         (m**2/m**2)
         &                                     purb, &            !< urbanisation                    (-)
         &                                     pfor_d, &          !< deciduous forest                (-)
         &                                     pfor_e, &          !< evergreen forest                (-)
         &                                     pskinc, &          !< skin conductivity               (W m-2 K-1)
         &                                     pemissivity, &     !< surface thermal emissivity      (-)
         &                                     prs_min, &         !< minimum stomata resistance      (s/m)
         &                                     lnhp, & 
         &                                     pwz0, &  ! weighted summand for z0
         &                                     area_tot, &    ! total area
         &                                     area_land, &   ! area with land
         &                                     area_plcov, &  ! area covered with plants
         &                                     bound_north_cosmo, &  !< northern boundary for COSMO target domain
         &                                     bound_south_cosmo, &  !< southern boundary for COSMO target domain
         &                                     bound_west_cosmo, &   !< western  boundary for COSMO target domain
         &                                     bound_east_cosmo, &   !< eastern  boundary for COSMO target domain
         &                                     hp ! height of Prandtl-layer


    LOGICAL                                 :: l_opn_gc_file(1:ntiles_globcover)
    CHARACTER (LEN=80)                      :: varname  !< name of variable

    ! Some stuff for OpenMP parallelization
#ifdef _OPENMP    
    REAL(KIND=wp)                           :: region_start, region_end, region_wallclock, loop_start, loop_end, loop_wallclock
#endif
    INTEGER(KIND=i4)                        :: num_blocks, ib, il, blk_len, istartlon, ishift
    INTEGER(KIND=i4)                        :: nblocks1, nblocks2, blk_len1, blk_len2, nlon_sub1, nlon_sub2, istartlon2, iendlon2

    !$   INTEGER :: omp_get_max_threads, omp_get_thread_num, thread_id
    !$   INTEGER (i4), ALLOCATABLE :: start_cell_arr(:)

    CALL logging%info('Enter routine: agg_globcover_data_to_target_grid')

    IF (l_use_corine) THEN
      CALL logging%info('Running with corine')
    ELSE
      CALL logging%info('Running with globcover')
    ENDIF

    fr_land_globcover(:,:,:)    = 0.0_wp
    ice_globcover(:,:,:)        = 0.0_wp
    z0_globcover(:,:,:)         = 0.0_wp
    root_globcover(:,:,:)       = 0.0_wp
    plcov_mx_globcover(:,:,:)   = 0.0_wp
    plcov_mn_globcover(:,:,:)   = 0.0_wp
    lai_mx_globcover(:,:,:)     = 0.0_wp
    lai_mn_globcover(:,:,:)     = 0.0_wp
    rs_min_globcover(:,:,:)     = 0.0_wp
    urban_globcover(:,:,:)      = 0.0_wp
    for_d_globcover(:,:,:)      = 0.0_wp
    for_e_globcover(:,:,:)      = 0.0_wp
    emissivity_globcover(:,:,:) = 0.0_wp

    apix_e  = re * re * deg2rad* ABS(globcover_grid%dlon_reg) * deg2rad * ABS(globcover_grid%dlat_reg)

    hp   = 30.0      ! height of Prandtl-layer
    lnhp = LOG(hp)

    default_real = 0.0
    undefined_integer= NINT(undefined)

    globcover_class_fraction = default_real
    globcover_class_npixel   = undefined_integer
    globcover_tot_npixel = undefined_integer

    a_weight = default_real
    a_class  = default_real

    SELECT CASE(tg%igrid_type)
      CASE(igrid_icon)  ! ICON GRID
        ke = 1
      CASE(igrid_cosmo)  ! COSMO GRID
        ke = 1
        bound_north_cosmo = MAXVAL(lat_geo) + 0.05_wp  ! add some "buffer"
        bound_north_cosmo = MIN(bound_north_cosmo,90.0_wp)
        bound_south_cosmo = MINVAL(lat_geo) - 0.05_wp  ! add some "buffer"
        bound_south_cosmo = MAX(bound_south_cosmo,-90.0_wp)

        bound_east_cosmo = MAXVAL(lon_geo) + 0.25_wp  ! add some "buffer"
        bound_east_cosmo = MIN(bound_east_cosmo,180.0_wp)
        bound_west_cosmo = MINVAL(lon_geo) - 0.25_wp  ! add some "buffer"
        bound_west_cosmo = MAX(bound_west_cosmo,-180.0_wp)

    END SELECT

    ! init lookup tables
    CALL init_globcover_lookup_tables(nclass_globcover, &
         &      ilookup_table_globcover, &
         &      z0_lt_globcover,            &
         &      lnz0_lt_globcover,          &
         &      plc_mn_lt_globcover,        &
         &      plc_mx_lt_globcover,        &
         &      lai_mn_lt_globcover,        &
         &      lai_mx_lt_globcover,        &
         &      rd_lt_globcover,            &
         &      skinc_lt_globcover,         &
         &      emiss_lt_globcover,         &
         &      rs_min_lt_globcover)

    CALL get_name_globcover_lookup_tables(ilookup_table_globcover, name_lookup_table_globcover)

    nlon = globcover_grid%nlon_reg
    ALLOCATE(ie_vec(nlon),je_vec(nlon),ke_vec(nlon))
    ie_vec(:) = 0
    je_vec(:) = 0
    ke_vec(:) = 0
    start_cell_id = 1

    ! open netcdf file
    DO nt = 1,ntiles_globcover
      CALL check_netcdf( nf90_open(TRIM(globcover_file(nt)),NF90_NOWRITE, ncid_globcover(nt)))
    END DO

    varname = 'GLOBCOVER' ! I know that the globcover data are stored in a variable called 'GLOBCOVER'

    CALL check_netcdf(nf90_inq_varid(ncid_globcover(1), TRIM(varname), varid_globcover))

    mlat = 1
    block_row_start = mlat

    CALL det_band_globcover_data(globcover_grid,block_row_start,ta_grid)

    IF (ALLOCATED(lu_block)) THEN
      DEALLOCATE(lu_block, STAT = errorcode)
      IF(errorcode /= 0) CALL logging%error('logging%error deallocate the lu_block',__FILE__,__LINE__)
    END IF
    ALLOCATE (lu_block(1:ta_grid%nlon_reg, 1:ta_grid%nlat_reg), STAT = errorcode)
    IF (errorcode /= 0) CALL logging%error('logging%error allocate the lu_block',__FILE__,__LINE__)

    CALL get_globcover_data_block(ta_grid,               &
         globcover_tiles_grid,  &
         ncid_globcover,       &
         lu_block)

    block_row = 0

    ! Determine start and end longitude of search
    istartlon = 1
    iendlon = globcover_grid%nlon_reg
    IF (tg%igrid_type == igrid_icon) THEN
      DO i_col = 1, globcover_grid%nlon_reg
        point_lon = lon_globcover(i_col)
        IF (point_lon < tg%minlon) istartlon = i_col + 1
        IF (point_lon > tg%maxlon) THEN
          iendlon = i_col - 1
          EXIT
        ENDIF
      ENDDO
      ! second search for shifted longitudes to detect regional domains crossing the dateline
      ! (needed to optimize the 'domain decomposition' for this case)
      IF (tg%maxlon - tg%minlon > 360._wp .AND. tg%maxlon_s - tg%minlon_s < 360._wp) THEN
        DO i_col = 1, globcover_grid%nlon_reg
          point_lon = lon_globcover(i_col)
          IF (tg%maxlon_s > 180._wp .AND. point_lon + 360._wp < tg%maxlon_s .OR. &
              tg%maxlon_s < 180._wp .AND. point_lon < tg%maxlon_s) iendlon2 = i_col + 1
          IF (tg%minlon_s > 180._wp .AND. point_lon + 360._wp > tg%minlon_s .OR. &
              tg%minlon_s < 180._wp .AND. point_lon > tg%minlon_s) THEN
            istartlon2 = i_col - 1
            EXIT
          ENDIF
        ENDDO
        WRITE(message_text,*) 'Limited-area domain crossing the dateline detected'
        CALL logging%info(message_text)
        WRITE(message_text,*) 'End and start index of partial domains', iendlon2, istartlon2
        CALL logging%info(message_text)
      ELSE
        iendlon2 = 0
        istartlon2 = 1
      ENDIF
    ELSE IF (tg%igrid_type == igrid_cosmo) THEN
      DO i_col = 1, globcover_grid%nlon_reg
        point_lon = lon_globcover(i_col)
        IF (point_lon < bound_west_cosmo) istartlon = i_col + 1
        IF (point_lon > bound_east_cosmo) THEN
          iendlon = i_col - 1
          EXIT
        ENDIF
      ENDDO
      iendlon2 = 0
      istartlon2 = 1
    ENDIF

    IF (iendlon2 == 0) THEN
       nlon_sub1 = iendlon - istartlon + 1
       nlon_sub2 = 0
     ELSE
       nlon_sub1 = iendlon2 - istartlon + 1
       nlon_sub2 = iendlon - istartlon2 + 1
     ENDIF

     num_blocks = 1
     blk_len2   = 0
     !$ num_blocks = omp_get_max_threads()
     IF (num_blocks > 1 .AND. nlon_sub2 > 0) THEN
       nblocks1 = NINT(REAL(num_blocks*nlon_sub1,wp)/REAL(nlon_sub1+nlon_sub2,wp))
       nblocks2 = num_blocks - nblocks1
     ELSE
       nblocks1 = num_blocks
       nblocks2 = 0
     ENDIF
     IF (MOD(nlon_sub1,nblocks1)== 0) THEN
       blk_len1 = nlon_sub1/nblocks1
     ELSE
       blk_len1 = nlon_sub1/nblocks1 + 1
     ENDIF
     IF (nblocks2 > 0) THEN
       IF (MOD(nlon_sub2,num_blocks)== 0) THEN
         blk_len2 = nlon_sub2/nblocks2
       ELSE
         blk_len2 = nlon_sub2/nblocks2 + 1
       ENDIF
     ELSE
       blk_len2 = 0
     ENDIF
     !$ allocate(start_cell_arr(num_blocks))
     !$ start_cell_arr(:) = 1

     WRITE(message_text,*) 'nlon_sub1/2, nblocks1/2, blk_len1/2: ',nlon_sub1, nlon_sub2, nblocks1, nblocks2, blk_len1, blk_len2
     CALL logging%info(message_text)

#ifdef _OPENMP
    region_wallclock = 0.0_wp
#endif

    CALL logging%info('Start loop over globcover rows...')

    globcover_rows: DO mlat = 1,globcover_grid%nlat_reg

#ifdef _OPENMP
      loop_start = omp_get_wtime()
#endif
      block_row= block_row + 1
      IF(block_row > ta_grid%nlat_reg) THEN ! read in new block
        block_row_start = mlat
        block_row = 1
        CALL det_band_globcover_data(globcover_grid,block_row_start,ta_grid)
        IF(ALLOCATED(lu_block)) THEN
          DEALLOCATE(lu_block, STAT=errorcode)
          IF(errorcode/=0) CALL logging%error('logging%error deallocate the lu_block',__FILE__,__LINE__)
        ENDIF
        ALLOCATE (lu_block(1:ta_grid%nlon_reg,1:ta_grid%nlat_reg), STAT=errorcode)
        IF(errorcode/=0) CALL logging%error('logging%error allocate lu_block',__FILE__,__LINE__)
        CALL get_globcover_data_block(ta_grid,              &
             &                        globcover_tiles_grid, &
             &                        ncid_globcover,       &
             &                        lu_block)
      ENDIF

      point_lat = lat_globcover(mlat)

      IF (tg%igrid_type == igrid_cosmo) THEN ! CASE COSMO grid, omit some I/O from hard disk when out of target domain
        IF ((point_lat > bound_north_cosmo).OR.(point_lat < bound_south_cosmo) ) THEN ! raw data out of target grid
          CYCLE globcover_rows
        ENDIF
      ELSE IF (tg%igrid_type == igrid_icon) THEN
        IF (point_lat > tg%maxlat .OR. point_lat < tg%minlat) THEN
          CYCLE globcover_rows
        ENDIF
      ENDIF ! grid type

      globcover_data_row(1:nlon) = lu_block(1:nlon,block_row)
      apix = apix_e * COS(point_lat * deg2rad) ! area of raw data pixel (in [m**2])

      ie_vec(istartlon:iendlon) = 0
      IF (tg%igrid_type /= igrid_icon) THEN
        je_vec(:) = 0
        ke_vec(:) = 0
      ENDIF
#ifdef _OPENMP
      region_start = omp_get_wtime()
!$OMP PARALLEL DO PRIVATE(ib,il,i_col,ii1,ii2,ishift,blk_len,point_lon,thread_id,start_cell_id)
#endif
      DO ib = 1, num_blocks

        !$     thread_id = omp_get_thread_num()+1
        !$     start_cell_id = start_cell_arr(thread_id)

        IF (ib <= nblocks1) THEN
          ishift = istartlon-1+(ib-1)*blk_len1
          blk_len = blk_len1
        ELSE
          ishift = istartlon2-1+(ib-(nblocks1+1))*blk_len2
          blk_len = blk_len2
        ENDIF
        ! Prevent truncation errors near the dateline
        IF (ib == num_blocks .AND. tg%maxlon > 179._wp) blk_len = nlon-ishift

        columns1: DO il = 1, blk_len
          i_col = ishift+il
          IF (i_col > iendlon) CYCLE columns1

          ! find the corresponding target grid indices
          point_lon = lon_globcover(i_col)

          ! Reset start cell when entering a new row or when the previous data point was outside
          ! the model domain
          IF (tg%igrid_type == igrid_icon .AND. (il == 1 .OR. start_cell_id == 0)) THEN
            ii1 = NINT(point_lon*search_res)
            ii2 = NINT(point_lat*search_res)
            start_cell_id = tg%search_index(ii1,ii2)
            IF (start_cell_id == 0) EXIT ! in this case, the whole row is empty; may happen with merged (non-contiguous) domains
          ENDIF

          CALL  find_nearest_target_grid_element( point_lon,     &
               &      point_lat,     &
               &      tg,            &
               &      start_cell_id, &
               &      ie_vec(i_col), &
               &      je_vec(i_col), &
               &      ke_vec(i_col)  )

        ENDDO columns1
        !$     start_cell_arr(thread_id) = start_cell_id
      ENDDO
#ifdef _OPENMP
!$OMP END PARALLEL DO
      region_end = omp_get_wtime()      
      region_wallclock = region_end - region_start
#endif
      
      columns2: DO i_col = istartlon, iendlon
        ! find the corresponding target grid indices

        ie = ie_vec(i_col)
        je = je_vec(i_col)
        ke = ke_vec(i_col)

        IF ((ie /= 0).AND.(je /= 0).AND.(ke /= 0)) THEN
          ! raw data pixel within target grid, see output of routine find_rotated_lonlat_grid_element_index

          lu = globcover_data_row(i_col)                        ! land use class
          IF (l_use_corine) THEN
            CALL corinecover_look_up(lu,    &
                &      nclass_globcover,   &
                &      lnz0_lt_globcover,  &
                &      plc_mn_lt_globcover,&
                &      plc_mx_lt_globcover,&
                &      lai_mn_lt_globcover,&
                &      lai_mx_lt_globcover,&
                &      rd_lt_globcover,    &
                &      skinc_lt_globcover, &
                &      emiss_lt_globcover, &
                &      rs_min_lt_globcover,&
                &      pland,              &
                &      pice,               &
                &      plnz0,              &
                &      proot,              &
                &      pmn,                &
                &      pmx,                &
                &      plaimn,             &
                &      plaimx,             &
                &      purb,               &
                &      pfor_d,             &
                &      pfor_e,             &
                &      pskinc,             &
                &      pemissivity,        &
                &      prs_min,            &
                &      k_error)
          ELSE
            CALL globcover_look_up(lu,    &
                &      nclass_globcover,   &
                &      lnz0_lt_globcover,  &
                &      plc_mn_lt_globcover,&
                &      plc_mx_lt_globcover,&
                &      lai_mn_lt_globcover,&
                &      lai_mx_lt_globcover,&
                &      rd_lt_globcover,    &
                &      skinc_lt_globcover, &
                &      emiss_lt_globcover, &
                &      rs_min_lt_globcover,&
                &      pland,              &
                &      pice,               &
                &      plnz0,              &
                &      proot,              &
                &      pmn,                &
                &      pmx,                &
                &      plaimn,             &
                &      plaimx,             &
                &      purb,               &
                &      pfor_d,             &
                &      pfor_e,             &
                &      pskinc,             &
                &      pemissivity,        &
                &      prs_min,            &
                &      k_error)
          ENDIF

          IF (k_error == 0) THEN ! valid land use class

            globcover_tot_npixel(ie,je,ke) = globcover_tot_npixel(ie,je,ke) + 1
            a_weight(ie,je,ke) = a_weight(ie,je,ke) + apix  ! sum up area for weight
            IF(l_use_corine ) THEN
              CALL get_corinecover_idx(lu,nclass)
            ELSE
              CALL get_globcover_idx(lu,nclass)
            ENDIF
            globcover_class_npixel(ie,je,ke,nclass) = globcover_class_npixel(ie,je,ke,nclass) + 1
            a_class(ie,je,ke,nclass) = a_class(ie,je,ke,nclass) + apix   ! sum area of valid land use pixels
            ! (use as weight later)
            emissivity_globcover(ie,je,ke) =  emissivity_globcover(ie,je,ke) + apix * pemissivity
            IF (pland >  0.0) THEN ! only for land pixel

              ! weighted with whole area
              fr_land_globcover(ie,je,ke) = fr_land_globcover(ie,je,ke) + apix * pland
              ice_globcover(ie,je,ke) = ice_globcover(ie,je,ke) + apix * pice
              urban_globcover(ie,je,ke) = urban_globcover(ie,je,ke) + apix * purb
              ! z0 is averaged logarithmic
              IF ( lnhp /= plnz0) THEN ! z0 is averaged logarithmic
                pwz0 = 1./(lnhp - plnz0)
              ELSE
                pwz0 = 0.
              ENDIF
              z0_globcover(ie,je,ke)      = z0_globcover(ie,je,ke) + apix * pwz0
              plcov_mn_globcover(ie,je,ke) = plcov_mn_globcover(ie,je,ke) + apix * pmn
              plcov_mx_globcover(ie,je,ke) = plcov_mx_globcover(ie,je,ke) + apix * pmx
              skinc_globcover(ie,je,ke)    = skinc_globcover(ie,je,ke)    + apix * pskinc
              ! the following fields are weighted with the plant cover
              root_globcover(ie,je,ke) = root_globcover(ie,je,ke) + apix * pmx * proot
              lai_mn_globcover(ie,je,ke) = lai_mn_globcover(ie,je,ke) + apix * pmx * plaimn
              lai_mx_globcover(ie,je,ke) = lai_mx_globcover(ie,je,ke) + apix * pmx * plaimx
              rs_min_globcover(ie,je,ke) = rs_min_globcover(ie,je,ke) + apix * pmx * prs_min
              for_d_globcover(ie,je,ke) = for_d_globcover(ie,je,ke) + apix * pmx * pfor_d
              for_e_globcover(ie,je,ke) = for_e_globcover(ie,je,ke) + apix * pmx * pfor_e

            END IF
          ENDIF
        ENDIF
        ! end loops
      ENDDO columns2
#ifdef _OPENMP
      IF (mlat == 1.OR.(MOD(mlat,3600) == 0)) THEN
        loop_end = omp_get_wtime()
        loop_wallclock = loop_end-loop_start
        WRITE(message_text, '(a,i6,a,f7.2,a,f18.12,a,f18.12,a)') &
             & 'GLOBCOVER row:', mlat, ' latitude: ', lat_globcover(mlat), &
             & ' time: ', loop_wallclock, ' s openmp: ', region_wallclock, ' s'      
        CALL logging%info(message_text)
      ENDIF
#else
      IF (MOD(mlat,200) == 0) THEN
        WRITE(message_text,'(a,i6,a,f7.2)') 'GLOBCOVER row:', mlat, ' latitude: ', lat_globcover(mlat)
        CALL logging%info(message_text)
      ENDIF
#endif
    ENDDO globcover_rows

    DEALLOCATE(ie_vec,je_vec,ke_vec)
    !$   DEALLOCATE(start_cell_arr)
    CALL logging%info('...done')

    ! calculate globcover_class_fraction (globcover_class_fraction/globcover_class_npixel)
    DO ke=1, tg%ke
      DO je=1, tg%je
        DO ie=1, tg%ie
          area_tot = a_weight(ie,je,ke)
          area_land = fr_land_globcover(ie,je,ke)

          ! weight by total area
          IF (area_tot > 0.0) THEN
            emissivity_globcover(ie,je,ke) = emissivity_globcover(ie,je,ke) / area_tot
            DO l=1,nclass_globcover
              globcover_class_fraction(ie,je,ke,l) =  a_class(ie,je,ke,l) /  area_tot
              ! area fraction of each land use class
            ENDDO
          ENDIF

          IF (area_land > 0.0 ) THEN
            area_plcov = plcov_mx_globcover(ie,je,ke) ! area covered with plants

            z0_globcover(ie,je,ke) = z0_globcover(ie,je,ke) / area_land
            z0_globcover(ie,je,ke) = hp * EXP(-1./z0_globcover(ie,je,ke))

            ! weight by total area
            fr_land_globcover(ie,je,ke) = fr_land_globcover(ie,je,ke) / area_tot
            ice_globcover(ie,je,ke)   = ice_globcover(ie,je,ke)   / area_tot


            ! weight by land area
            urban_globcover(ie,je,ke) = urban_globcover(ie,je,ke) / area_land
            for_d_globcover(ie,je,ke) = for_d_globcover(ie,je,ke) / area_land
            for_e_globcover(ie,je,ke) = for_e_globcover(ie,je,ke) / area_land
            plcov_mn_globcover(ie,je,ke) = plcov_mn_globcover(ie,je,ke) / area_land

            plcov_mx_globcover(ie,je,ke) = plcov_mx_globcover(ie,je,ke) / area_land
            skinc_globcover(ie,je,ke) = skinc_globcover(ie,je,ke) / area_land
            ! weight by area covered with plants
            IF (area_plcov > 0.0) THEN
              root_globcover(ie,je,ke) = root_globcover(ie,je,ke)     / area_plcov
              lai_mn_globcover(ie,je,ke) = lai_mn_globcover(ie,je,ke) / area_plcov
              lai_mx_globcover(ie,je,ke) = lai_mx_globcover(ie,je,ke) / area_plcov
              rs_min_globcover(ie,je,ke) = rs_min_globcover(ie,je,ke) / area_plcov
            ELSE ! may occur for ice grid elements
              root_globcover(ie,je,ke) = undefined
              lai_mn_globcover(ie,je,ke) = undefined
              lai_mx_globcover(ie,je,ke) = undefined
              rs_min_globcover(ie,je,ke) = rs_min_undef
            ENDIF

          ELSE IF ( area_tot > 0.0) THEN ! only sea pixels were found
            fr_land_globcover(ie,je,ke) = undefined
            z0_globcover(ie,je,ke)      = undefined
            ice_globcover(ie,je,ke)     = undefined
            urban_globcover(ie,je,ke)   = undefined
            for_d_globcover(ie,je,ke)   = undefined
            for_e_globcover(ie,je,ke)   = undefined
            plcov_mx_globcover(ie,je,ke)= undefined
            plcov_mn_globcover(ie,je,ke)= undefined

            root_globcover(ie,je,ke)    = undefined
            lai_mx_globcover(ie,je,ke)  = undefined
            lai_mn_globcover(ie,je,ke)  = undefined
            skinc_globcover(ie,je,ke)   = undefined
            rs_min_globcover(ie,je,ke)  = undefined

          ENDIF
        ENDDO
      ENDDO
    ENDDO

    l_opn_gc_file = .FALSE.

    DO ke=1, tg%ke
      DO je=1, tg%je
        DO ie=1, tg%ie
          IF (globcover_tot_npixel(ie,je,ke) == 0) THEN ! nearest neighbor search

            point_lon = lon_geo(ie,je,ke)
            point_lat = lat_geo(ie,je,ke)

            CALL find_reg_lonlat_grid_element_index(point_lon,      &
                 &                                     point_lat,      &
                 &                                     globcover_grid, &
                 &                                     i_lu,           &
                 &                                     j_lu,           &
                 &                                     ntiles_globcover,&
                 &                                     tile = tile,    &
                 &                                     regular_tiles_grid_info=globcover_tiles_grid)

            IF ((i_lu /= 0).AND.(j_lu /= 0))THEN

              i_col = i_lu
              j_row = j_lu

              IF (.NOT. l_opn_gc_file(tile)) THEN
                CALL check_netcdf(nf90_open(TRIM(globcover_file(tile)),NF90_NOWRITE, ncid_globcover(tile)), &
                     &            __FILE__, &
                     &            __LINE__)
                CALL check_netcdf(nf90_inq_varid(ncid_globcover(tile),TRIM(varname), varid_gc(tile)))
                l_opn_gc_file(tile) = .TRUE.
              ENDIF

              CALL check_netcdf(nf90_get_var(ncid_globcover(tile), varid_gc(tile),  &
                   &               globcover_data_pixel,  &
                   &               start=(/ i_col,j_row /),count=(/ 1,1 /)), __FILE__, __LINE__)

              lu = globcover_data_pixel(1,1)

               IF(l_use_corine ) THEN
                 CALL corinecover_look_up(lu, &
                       &      nclass_globcover, &
                       &      lnz0_lt_globcover,          &
                       &      plc_mn_lt_globcover,        &
                       &      plc_mx_lt_globcover,        &
                       &      lai_mn_lt_globcover,        &
                       &      lai_mx_lt_globcover,        &
                       &      rd_lt_globcover,            &
                       &      skinc_lt_globcover,         &
                       &      emiss_lt_globcover,         &
                       &      rs_min_lt_globcover,        &
                       &      pland,          &
                       &      pice,           &
                       &      plnz0,          &
                       &      proot,          &
                       &      pmn,            &
                       &      pmx,            &
                       &      plaimn,         &
                       &      plaimx,         &
                       &      purb,           &
                       &      pfor_d,         &
                       &      pfor_e,         &
                       &      pskinc,         &
                       &      pemissivity,    &
                       &      prs_min,        &
                       &      k_error)
               ELSE
                  CALL globcover_look_up(lu, &
                       &      nclass_globcover, &
                       &      lnz0_lt_globcover,          &
                       &      plc_mn_lt_globcover,        &
                       &      plc_mx_lt_globcover,        &
                       &      lai_mn_lt_globcover,        &
                       &      lai_mx_lt_globcover,        &
                       &      rd_lt_globcover,            &
                       &      skinc_lt_globcover,         &
                       &      emiss_lt_globcover,         &
                       &      rs_min_lt_globcover,        &
                       &      pland,          &
                       &      pice,           &
                       &      plnz0,          &
                       &      proot,          &
                       &      pmn,            &
                       &      pmx,            &
                       &      plaimn,         &
                       &      plaimx,         &
                       &      purb,           &
                       &      pfor_d,         &
                       &      pfor_e,         &
                       &      pskinc,         &
                       &      pemissivity,    &
                       &      prs_min,        &
                       &      k_error)
               ENDIF

            ELSE
              lu = 0
              k_error = 1
            ENDIF

            IF (k_error == 0) THEN ! valid land use class
              apix = apix_e * COS(point_lat * deg2rad) ! area of raw data pixel (in [m**2])
              IF(l_use_corine) THEN
                CALL get_corinecover_idx(lu,nclass)
              ELSE
                CALL get_globcover_idx(lu,nclass)
              ENDIF
              globcover_class_npixel(ie,je,ke,nclass) = globcover_class_npixel(ie,je,ke,nclass) + 1
              a_class(ie,je,ke,nclass) = a_class(ie,je,ke,nclass) + apix   ! sum area of valid land use pixels
              emissivity_globcover(ie,je,ke) =  pemissivity
              IF (pland >  0.0) THEN ! only for land pixel
                fr_land_globcover(ie,je,ke) = pland
                ice_globcover(ie,je,ke) =  pice
                urban_globcover(ie,je,ke) = purb
                IF ( lnhp /= plnz0) THEN ! log z0
                  pwz0 = 1./(lnhp - plnz0)
                ELSE
                  pwz0 = 0.
                ENDIF
                z0_globcover(ie,je,ke)      =  pwz0
                z0_globcover(ie,je,ke) = hp * EXP(-1./z0_globcover(ie,je,ke))
                plcov_mn_globcover(ie,je,ke) = pmn
                plcov_mx_globcover(ie,je,ke) = pmx
                skinc_globcover(ie,je,ke)    = pskinc 
                root_globcover(ie,je,ke) = proot
                lai_mn_globcover(ie,je,ke) = plaimn
                lai_mx_globcover(ie,je,ke) = plaimx
                rs_min_globcover(ie,je,ke) = prs_min
                for_d_globcover(ie,je,ke) = pfor_d * pmx
                for_e_globcover(ie,je,ke) = pfor_e * pmx
                IF (pice == 1.0) THEN  ! ice land use class
                  root_globcover(ie,je,ke) = undefined
                  lai_mn_globcover(ie,je,ke) = undefined
                  lai_mx_globcover(ie,je,ke) = undefined
                  rs_min_globcover(ie,je,ke) = rs_min_undef
                ENDIF
              ENDIF
            ELSE ! not a valid land use class
              fr_land_globcover(ie,je,ke) = undefined
              z0_globcover(ie,je,ke)      = undefined
              ice_globcover(ie,je,ke)     = undefined
              urban_globcover(ie,je,ke)   = undefined
              for_d_globcover(ie,je,ke)   = undefined
              for_e_globcover(ie,je,ke)   = undefined
              plcov_mx_globcover(ie,je,ke)= undefined
              plcov_mn_globcover(ie,je,ke)= undefined
              root_globcover(ie,je,ke)    = undefined
              lai_mx_globcover(ie,je,ke)  = undefined
              lai_mn_globcover(ie,je,ke)  = undefined
              skinc_globcover(ie,je,ke)   = undefined
              rs_min_globcover(ie,je,ke)  = undefined
            ENDIF
          ENDIF ! nearest neighbour search
        ENDDO
      ENDDO
    ENDDO

    DO tile=1,ntiles_globcover
      IF (l_opn_gc_file(tile)) THEN
        CALL check_netcdf(nf90_close(ncid_globcover(tile)))
      ENDIF
    ENDDO

    CALL logging%info('Exit routine: agg_globcover_data_to_target_grid')

  END SUBROUTINE agg_globcover_data_to_target_grid

END MODULE mo_agg_globcover
