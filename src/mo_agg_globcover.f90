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

  USE mo_kind, ONLY: dp, wp, i8, i2

  USE mo_utilities_extpar, ONLY: abort_extpar

  USE mo_grid_structures, ONLY: reg_lonlat_grid,     &
       &                        igrid_icon,          &
       &                        igrid_cosmo,         &
       &                        target_grid_def   !< type definition of structure for tg       

  USE mo_search_ll_grid, ONLY: find_reg_lonlat_grid_element_index

  USE mo_io_units,       ONLY: filename_max

  USE mo_io_utilities,   ONLY: check_netcdf

  USE mo_search_target_grid, ONLY: find_nearest_target_grid_element

  USE mo_globcover_data, ONLY: globcover_grid, &
       &                       lon_globcover,  &
       &                       lat_globcover,  &
       &                       ntiles_globcover

  USE mo_globcover_lookup_tables, ONLY: name_lookup_table_globcover, &
       &                                init_globcover_lookup_tables, &
       &                                get_name_globcover_lookup_tables, get_globcover_idx, &
       &                                z0_lt_globcover, lnz0_lt_globcover, plc_mn_lt_globcover, &
       &                                plc_mx_lt_globcover, lai_mn_lt_globcover, &
       &                                lai_mx_lt_globcover, rd_lt_globcover, emiss_lt_globcover, &
       &                                rs_min_lt_globcover, globcover_look_up

  USE mo_landuse_routines, ONLY: det_band_globcover_data, &
       &                         get_globcover_data_block

  USE mo_math_constants, ONLY: deg2rad

  USE mo_physical_constants, ONLY: re

  USE mo_target_grid_data, ONLY: lon_geo, & !< longitude coordinates of the COSMO grid in the geographical system
       &                         lat_geo, & !< latitude coordinates of the COSMO grid in the geographical system
       &                         search_res !< resolution of ICON grid search index list


  USE netcdf, ONLY: nf90_open,      &
       &            nf90_close,     &
       &            nf90_inq_varid, &
       &            nf90_get_var,   &
       &            nf90_nowrite

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
       &                                          emissivity_globcover)

    !-------------------------------------------------------------------------------------

    CHARACTER (LEN=filename_max), INTENT(IN) :: globcover_file(:)  !< filename globcover raw data
    INTEGER, INTENT(IN) :: ilookup_table_globcover
    REAL (wp), INTENT(IN) :: undefined            !< undef value
    TYPE(reg_lonlat_grid), INTENT(IN) :: globcover_tiles_grid(:)  ! grid structure of globcover tiles
    TYPE(target_grid_def), INTENT(IN) :: tg  !< structure with target grid description
    INTEGER, INTENT(IN) :: nclass_globcover !< globcover has 23 classes for the land use description
    REAL (wp), INTENT(OUT)  :: globcover_class_fraction(:,:,:,:)
    !< fraction for each globcover class on target grid (dimension (ie,je,ke,nclass_globcover))

    INTEGER (i8), INTENT(OUT) :: globcover_class_npixel(:,:,:,:)
    !< number of raw data pixels for each globcover class on target grid (dimension (ie,je,ke,nclass_globcover))


    INTEGER (i8), INTENT(OUT) :: globcover_tot_npixel(:,:,:)
    !< total number of globcover raw data pixels on target grid (dimension (ie,je,ke))


    REAL (wp), INTENT(OUT)  :: fr_land_globcover(:,:,:) !< fraction land due to globcover raw data
    REAL (wp), INTENT(OUT)  :: ice_globcover(:,:,:)     !< fraction of ice due to globcover raw data
    REAL (wp), INTENT(OUT)  :: z0_globcover(:,:,:)      !< roughness length due to globcover land use data
    REAL (wp), INTENT(OUT)  :: root_globcover(:,:,:)    !< root depth due to globcover land use data
    REAL (wp), INTENT(OUT)  :: plcov_mx_globcover(:,:,:)!< plant cover maximum due to globcover land use data
    REAL (wp), INTENT(OUT)  :: plcov_mn_globcover(:,:,:)!< plant cover minimum due to globcover land use data
    REAL (wp), INTENT(OUT)  :: lai_mx_globcover(:,:,:)  !< Leaf Area Index maximum due to globcover land use data
    REAL (wp), INTENT(OUT)  :: lai_mn_globcover(:,:,:)  !< Leaf Area Index minimum due to globcover land use data
    REAL (wp), INTENT(OUT)  :: rs_min_globcover(:,:,:)  !< minimal stomata resistance due to globcover land use data
    REAL (wp), INTENT(OUT)  :: urban_globcover(:,:,:)   !< urban fraction due to globcover land use data
    REAL (wp), INTENT(OUT)  :: for_d_globcover(:,:,:)   !< deciduous forest (fraction) due to globcover land use data
    REAL (wp), INTENT(OUT)  :: for_e_globcover(:,:,:)   !< evergreen forest (fraction) due to globcover land use data
    REAL (wp), INTENT(OUT)  :: emissivity_globcover(:,:,:) !< longwave emissivity due to globcover land use da

    ! structure with definition of the target area grid (dlon must be the same for the whole GLOBCOVER dataset)
    TYPE(reg_lonlat_grid):: ta_grid

    INTEGER (i8) :: undefined_integer ! undef value
    REAL (wp)    :: default_real


    INTEGER :: l      ! counters
    INTEGER :: nt           ! counter
    INTEGER :: i_col, j_row ! counter
    INTEGER (i8) :: i_lu, j_lu
    INTEGER (i8) :: ie, je, ke  ! indices for target grid elements
    INTEGER (i8), ALLOCATABLE :: ie_vec(:), je_vec(:), ke_vec(:)  ! indices for target grid elements
    INTEGER (i8) :: start_cell_id !< ID of starting cell for ICON search
    INTEGER (i8) :: ii1, ii2

    REAL (wp)    :: a_weight(1:tg%ie,1:tg%je,1:tg%ke)
    !< area weight of all raw data pixels in target grid
    REAL (wp)    :: a_class(1:tg%ie,1:tg%je,1:tg%ke,1:nclass_globcover)
    !< area for each land use class grid  in target grid element (for a area weight)
    INTEGER (i2), ALLOCATABLE:: lu_block(:,:)  ! a block of GLOBCOVER land use data

    REAL (wp)    :: apix      !< area of a raw data pixel
    REAL (wp)    :: apix_e      !< area of a raw data pixel at equator

    INTEGER (i2) :: globcover_data_row(globcover_grid%nlon_reg)
    INTEGER (i2) :: globcover_data_pixel(1:1,1:1)
    INTEGER :: lu  ! land use class
    INTEGER :: nclass ! index in array of globcover tables
    INTEGER :: ncid_globcover(1:ntiles_globcover)            !< netcdf unit file number
    CHARACTER (LEN=80) :: varname  !< name of variable
    INTEGER :: varid_globcover               !< id of variable
    INTEGER :: varid_gc(1:ntiles_globcover)
    LOGICAL :: l_opn_gc_file(1:ntiles_globcover)
    INTEGER :: nlon
    INTEGER :: block_row_start
    INTEGER :: block_row
    INTEGER :: mlat
    INTEGER :: tile

    REAL(wp)   :: point_lon, point_lat

    REAL (wp) :: pland          !< land cover                      (-)
    REAL (wp) :: pice           !< ice fraction                    (-)
    REAL (wp) :: plnz0          !< logarithm of roughness length   (m)
    REAL (wp) :: proot          !< root depth                      (m)
    REAL (wp) :: pmn            !< minimal plant cover             (-)
    REAL (wp) :: pmx            !< maximum plant cover             (-)
    REAL (wp) :: plaimn         !< minimal leaf area index         (m**2/m**2)
    REAL (wp) :: plaimx         !< maximum leaf area index         (m**2/m**2)
    REAL (wp) :: purb           !< urbanisation                    (-)
    REAL (wp) :: pfor_d         !< deciduous forest                (-)
    REAL (wp) :: pfor_e         !< evergreen forest                (-)
    REAL (wp) :: pemissivity    !< surface thermal emissivity      (-)
    REAL (wp) :: prs_min        !< minimum stomata resistance      (s/m)

    INTEGER        :: k_error, errorcode   ! error return code

    REAL           :: hp ! height of Prandtl-layer
    REAL (wp) :: lnhp
    REAL (wp) :: pwz0 ! weighted summand for z0

    REAL (wp) :: area_tot   ! total area
    REAL (wp) :: area_land  ! area with land
    REAL (wp) :: area_plcov ! area covered with plants

    REAL (wp) :: bound_north_cosmo !< northern boundary for COSMO target domain
    REAL (wp) :: bound_south_cosmo !< southern boundary for COSMO target domain
    REAL (wp) :: bound_west_cosmo  !< western  boundary for COSMO target domain
    REAL (wp) :: bound_east_cosmo  !< eastern  boundary for COSMO target domain

    ! Some stuff for OpenMP parallelization
#ifdef _OPENMP    
    REAL(dp) :: region_start, region_end, region_wallclock, loop_start, loop_end, loop_wallclock
#endif
    INTEGER :: num_blocks, ib, il, blk_len, istartlon, iendlon, nlon_sub, ishift
    !$   INTEGER :: omp_get_max_threads, omp_get_thread_num, thread_id
    !$   INTEGER (i8), ALLOCATABLE :: start_cell_arr(:)

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
    ! area of globcover raw data pixel at equator
    PRINT *,'area pixel at equator: ',apix_e

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
    PRINT *,'open GLOBCOVER files'
    ! >mes
    DO nt = 1,ntiles_globcover
      CALL check_netcdf( nf90_open(TRIM(globcover_file(nt)),NF90_NOWRITE, ncid_globcover(nt)))
    END DO
    ! <mes

    varname = 'GLOBCOVER' ! I know that the globcover data are stored in a variable called 'GLOBCOVER'

    CALL check_netcdf(nf90_inq_varid(ncid_globcover(1), TRIM(varname), varid_globcover))

    ! >mes
    mlat = 1
    block_row_start = mlat

    CALL det_band_globcover_data(globcover_grid,block_row_start,ta_grid)

    print*, 'ta_grid: ', ta_grid

    IF (ALLOCATED(lu_block)) THEN
      DEALLOCATE(lu_block, STAT = errorcode)
      IF(errorcode /= 0) CALL abort_extpar('Cant deallocate the lu_block')
    END IF
    ALLOCATE (lu_block(1:ta_grid%nlon_reg, 1:ta_grid%nlat_reg), STAT = errorcode)
    IF (errorcode /= 0) CALL abort_extpar('Cant allocate the lu_block')

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
    ELSE IF (tg%igrid_type == igrid_cosmo) THEN
      DO i_col = 1, globcover_grid%nlon_reg
        point_lon = lon_globcover(i_col)
        IF (point_lon < bound_west_cosmo) istartlon = i_col + 1
        IF (point_lon > bound_east_cosmo) THEN
          iendlon = i_col - 1
          EXIT
        ENDIF
      ENDDO
    ENDIF
    nlon_sub = iendlon - istartlon + 1

    num_blocks = 1
    !$   num_blocks = omp_get_max_threads()
    IF (MOD(nlon_sub,num_blocks)== 0) THEN
      blk_len = nlon_sub/num_blocks
    ELSE
      blk_len = nlon_sub/num_blocks + 1
    ENDIF
    !$   ALLOCATE(start_cell_arr(num_blocks))
    !$   start_cell_arr(:) = 1
    PRINT*, 'nlon_sub, num_blocks, blk_len: ',nlon_sub, num_blocks, blk_len

#ifdef _OPENMP
    region_wallclock = 0.0_dp
#endif
    print*, 'Start loop over GLOBCOVER rows'
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
          IF(errorcode/=0) CALL abort_extpar('Cant deallocate the lu_block')
        ENDIF
        ALLOCATE (lu_block(1:ta_grid%nlon_reg,1:ta_grid%nlat_reg), STAT=errorcode)
        IF(errorcode/=0) CALL abort_extpar('Cant allocate lu_block')
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
!$OMP PARALLEL DO PRIVATE(ib,il,i_col,ii1,ii2,ishift,point_lon,thread_id,start_cell_id)
#endif
      DO ib = 1, num_blocks

        !$     thread_id = omp_get_thread_num()+1
        !$     start_cell_id = start_cell_arr(thread_id)
        ishift = istartlon-1+(ib-1)*blk_len

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

          CALL globcover_look_up(lu,    &
               &      nclass_globcover,   &
               &      lnz0_lt_globcover,  &
               &      plc_mn_lt_globcover,&
               &      plc_mx_lt_globcover,&
               &      lai_mn_lt_globcover,&
               &      lai_mx_lt_globcover,&
               &      rd_lt_globcover,    &
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
               &      pemissivity,        &
               &      prs_min,            &
               &      k_error)


          IF (k_error == 0) THEN ! valid land use class

            globcover_tot_npixel(ie,je,ke) = globcover_tot_npixel(ie,je,ke) + 1
            a_weight(ie,je,ke) = a_weight(ie,je,ke) + apix  ! sum up area for weight

            CALL get_globcover_idx(lu,nclass)
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
        PRINT '(a,i6,a,f7.2,a,f18.12,a,f18.12,a)', &
             & 'GLOBCOVER row:', mlat, ' latitude: ', lat_globcover(mlat), ' time: ', loop_wallclock, ' s openmp: ', region_wallclock, ' s'      
      ENDIF
#else
      IF (MOD(mlat,200) == 0) PRINT '(a,i6,a,f7.2)', 'GLOBCOVER row:', mlat, ' latitude: ', lat_globcover(mlat)
#endif
    ENDDO globcover_rows

    DEALLOCATE(ie_vec,je_vec,ke_vec)
    !$   DEALLOCATE(start_cell_arr)

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
                 &                                     tile = tile,    &
                 &                                     regular_tiles_grid_info=globcover_tiles_grid)

            IF ((i_lu /= 0).AND.(j_lu /= 0))THEN

              i_col = i_lu
              j_row = j_lu

              IF (.NOT. l_opn_gc_file(tile)) THEN
                CALL check_netcdf(nf90_open(TRIM(globcover_file(tile)),NF90_NOWRITE, ncid_globcover(tile)), __FILE__, __LINE__)
              ENDIF

              CALL check_netcdf(nf90_get_var(ncid_globcover(tile), varid_gc(tile),  &
                   &               globcover_data_pixel,  &
                   &               start=(/ i_col,j_row /),count=(/ 1,1 /)), __FILE__, __LINE__)

              lu = globcover_data_pixel(1,1)


              CALL globcover_look_up(lu, &
                   &      nclass_globcover, &
                   &      lnz0_lt_globcover,          &
                   &      plc_mn_lt_globcover,        &
                   &      plc_mx_lt_globcover,        &
                   &      lai_mn_lt_globcover,        &
                   &      lai_mx_lt_globcover,        &
                   &      rd_lt_globcover,            &
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
                   &      pemissivity,    &
                   &      prs_min,        &
                   &      k_error)
            ELSE
              lu = 0
              k_error = 1
            ENDIF

            IF (k_error == 0) THEN ! valid land use class
              apix = apix_e * COS(point_lat * deg2rad) ! area of raw data pixel (in [m**2])
              CALL get_globcover_idx(lu,nclass)
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

PRINT*,' MAX ICE_GLOBCOVER: ', MAXVAL(ice_globcover)

  END SUBROUTINE agg_globcover_data_to_target_grid

END MODULE mo_agg_globcover
