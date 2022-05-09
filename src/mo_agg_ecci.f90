!+ Fortran module to aggregate ecci land use data to a target grid
!
!
MODULE mo_agg_ecci

  USE mo_logging
  USE mo_kind,                  ONLY:  wp, i4, i2

  USE mo_grid_structures,       ONLY: reg_lonlat_grid,     &
       &                              igrid_icon,          &
       &                              igrid_cosmo,         &
       &                              target_grid_def   !< type definition of structure for tg       

  USE mo_search_ll_grid,        ONLY: find_reg_lonlat_grid_element_index
                                
  USE mo_io_utilities,          ONLY: check_netcdf

  USE mo_search_target_grid,    ONLY: find_nearest_target_grid_element

  USE mo_ecci_data,             ONLY: ecci_grid, &
       &                              lon_ecci,  &
       &                              lat_ecci,  &
       &                              ntiles_ecci

  USE mo_ecci_lookup_tables,    ONLY: name_lookup_table_ecci, &
       &                              init_ecci_lookup_tables, &
       &                              get_name_ecci_lookup_tables, get_ecci_idx, &
       &                              z0_lt_ecci, lnz0_lt_ecci, plc_mn_lt_ecci, &
       &                              plc_mx_lt_ecci, lai_mn_lt_ecci, &
       &                              lai_mx_lt_ecci, rd_lt_ecci, skinc_lt_ecci, &
       &                              emiss_lt_ecci, rs_min_lt_ecci, ecci_look_up

  USE mo_landuse_routines,      ONLY: det_band_ecci_data, &
       &                              get_ecci_data_block

  USE mo_math_constants,        ONLY: deg2rad

  USE mo_physical_constants,    ONLY: re

  USE mo_target_grid_data,      ONLY: lon_geo, & !< longitude coordinates of the COSMO grid in the geographical system
       &                              lat_geo, & !< latitude coordinates of the COSMO grid in the geographical system
       &                              search_res !< resolution of ICON grid search index list


  USE netcdf,                  ONLY: nf90_open,      &
       &                             nf90_close,     &
       &                             nf90_inq_varid, &
       &                             nf90_get_var,   &
       &                             nf90_nowrite


  IMPLICIT NONE

  PRIVATE

  PUBLIC :: agg_ecci_data_to_target_grid

  REAL(KIND=wp), PARAMETER :: rs_min_undef=999. !< undefined value for minimal stomata resistance


CONTAINS

  !> Subroutine to aggregate ecci data to the target grid
  SUBROUTINE agg_ecci_data_to_target_grid(ecci_file,          &
       &                                  ilookup_table_ecci, &
       &                                  undefined,               &
       &                                  ecci_tiles_grid,    &
       &                                  tg,                      &
       &                                  nclass_ecci,        &
       &                                  ecci_class_fraction,&
       &                                  ecci_class_npixel,  &
       &                                  ecci_tot_npixel,    &
       &                                  fr_land_ecci ,      &
       &                                  ice_ecci,           &
       &                                  z0_ecci,            &
       &                                  root_ecci,          &
       &                                  plcov_mn_ecci,      &
       &                                  plcov_mx_ecci,      &
       &                                  lai_mn_ecci,        &
       &                                  lai_mx_ecci,        &
       &                                  rs_min_ecci,        &
       &                                  urban_ecci,         &
       &                                  for_d_ecci,         &
       &                                  for_e_ecci,         &
       &                                  skinc_ecci,         &
       &                                  emissivity_ecci)

    CHARACTER (LEN=*), INTENT(IN)     :: ecci_file(:)  !< filename ecci raw data
    REAL (wp), INTENT(IN)             :: undefined            !< undef value
    TYPE(reg_lonlat_grid), INTENT(IN) :: ecci_tiles_grid(:)  ! grid structure of ecci tiles
    TYPE(target_grid_def), INTENT(IN) :: tg  !< structure with target grid description

    INTEGER (KIND=i4), INTENT(IN)     :: nclass_ecci, & !< ecci has 23 classes for the land use description
         &                              ilookup_table_ecci
    !< fraction for each ecci class on target grid (dimension (ie,je,ke,nclass_ecci))

    INTEGER (KIND=i4), INTENT(OUT)    :: ecci_class_npixel(:,:,:,:), &
         &                               ecci_tot_npixel(:,:,:)

    REAL (KIND=wp), INTENT(OUT)      :: fr_land_ecci(:,:,:), & !< fraction land due to ecci raw data
         &                              ice_ecci(:,:,:), &     !< fraction of ice due to ecci raw data
         &                              z0_ecci(:,:,:), &      !< roughness length due to ecci land use data
         &                              root_ecci(:,:,:), &    !< root depth due to ecci land use data
         &                              plcov_mx_ecci(:,:,:), &!< plant cover maximum due to ecci land use data
         &                              plcov_mn_ecci(:,:,:), &!< plant cover minimum due to ecci land use data
         &                              lai_mx_ecci(:,:,:), &  !< Leaf Area Index maximum due to ecci land use data
         &                              lai_mn_ecci(:,:,:), &  !< Leaf Area Index minimum due to ecci land use data
         &                              rs_min_ecci(:,:,:), &  !< minimal stomata resistance due to ecci land use data
         &                              urban_ecci(:,:,:), &   !< urban fraction due to ecci land use data
         &                              for_d_ecci(:,:,:), &   !< deciduous forest (fraction) due to ecci land use data
         &                              for_e_ecci(:,:,:), &   !< evergreen forest (fraction) due to ecci land use data
         &                              skinc_ecci(:,:,:), &   !< skin conductivity due to ecci land use data
         &                              emissivity_ecci(:,:,:), & !< longwave emissivity due to ecci land use da
         &                              ecci_class_fraction(:,:,:,:)

    ! local variables

    ! structure with definition of the target area grid (dlon must be the same for the whole ECCI dataset)
    TYPE(reg_lonlat_grid)            :: ta_grid

    REAL (KIND=wp)                   :: default_real, &
         &                              a_weight(1:tg%ie,1:tg%je,1:tg%ke), &
         &                              a_class(1:tg%ie,1:tg%je,1:tg%ke,1:nclass_ecci), &
         &                              point_lon, point_lat, &
         &                              pland, &          !< land cover                      (-)
         &                              pice, &           !< ice fraction                    (-)
         &                              plnz0, &          !< logarithm of roughness length   (m)
         &                              proot, &          !< root depth                      (m)
         &                              pmn, &            !< minimal plant cover             (-)
         &                              pmx, &            !< maximum plant cover             (-)
         &                              plaimn, &         !< minimal leaf area index         (m**2/m**2)
         &                              plaimx, &         !< maximum leaf area index         (m**2/m**2)
         &                              purb, &           !< urbanisation                    (-)
         &                              pfor_d, &         !< deciduous forest                (-)
         &                              pfor_e, &         !< evergreen forest                (-)
         &                              pskinc, &         !< skin conductivity               (W m-2 K-1)
         &                              pemissivity, &    !< surface thermal emissivity      (-)
         &                              prs_min, &        !< minimum stomata resistance      (s/m)
         &                              hp, & ! height of Prandtl-layer
         &                              lnhp, &
         &                              pwz0, & ! weighted summand for z0
         &                              area_tot, &   ! total area
         &                              area_land, &  ! area with land
         &                              area_plcov, & ! area covered with plants
         &                              bound_north_cosmo, & !< northern boundary for COSMO target domain
         &                              bound_south_cosmo, & !< southern boundary for COSMO target domain
         &                              bound_west_cosmo, &  !< western  boundary for COSMO target domain
         &                              bound_east_cosmo, &  !< eastern  boundary for COSMO target domain
         &                              apix, &      !< area of a raw data pixel
         &                              apix_e      !< area of a raw data pixel at equator

    INTEGER (KIND=i4)               :: undefined_integer, & ! undef value
         &                             l, &      ! counters
         &                             nt, &           ! counter
         &                             i_col, j_row, & ! counter
         &                             i_lu, j_lu, &
         &                             ie, je, ke, &  ! indices for target grid elements
         &                             start_cell_id, & !< ID of starting cell for ICON search
         &                             ii1, ii2, &
         &                             lu, &  ! land use class
         &                             nclass, & ! index in array of ecci tables
         &                             ncid_ecci(1:ntiles_ecci), & !< netcdf unit file number
         &                             varid_ecci, &  !< id of variable
         &                             varid_gc(1:ntiles_ecci), &
         &                             nlon, &
         &                             block_row_start, &
         &                             block_row, &
         &                             mlat, &
         &                             tile, &
         &                             k_error, errorcode   ! error return code

    INTEGER (KIND=i4), ALLOCATABLE  :: ie_vec(:), je_vec(:), ke_vec(:)  ! indices for target grid elements

    INTEGER (KIND=i2)               :: ecci_data_row(ecci_grid%nlon_reg), &
         &                             ecci_data_pixel(1:1,1:1)

    INTEGER (KIND=i2), ALLOCATABLE  :: lu_block(:,:)  ! a block of ECCI land use data

    CHARACTER (LEN=80)              :: varname  !< name of variable

    LOGICAL                         :: l_opn_gc_file(1:ntiles_ecci)

    ! Some stuff for OpenMP parallelization
    INTEGER :: num_blocks, ib, il, blk_len, istartlon, iendlon, nlon_sub, ishift

    fr_land_ecci(:,:,:)    = 0.0_wp
    ice_ecci(:,:,:)        = 0.0_wp
    z0_ecci(:,:,:)         = 0.0_wp
    root_ecci(:,:,:)       = 0.0_wp
    plcov_mx_ecci(:,:,:)   = 0.0_wp
    plcov_mn_ecci(:,:,:)   = 0.0_wp
    lai_mx_ecci(:,:,:)     = 0.0_wp
    lai_mn_ecci(:,:,:)     = 0.0_wp
    rs_min_ecci(:,:,:)     = 0.0_wp
    urban_ecci(:,:,:)      = 0.0_wp
    for_d_ecci(:,:,:)      = 0.0_wp
    for_e_ecci(:,:,:)      = 0.0_wp
    emissivity_ecci(:,:,:) = 0.0_wp

    apix_e  = re * re * deg2rad* ABS(ecci_grid%dlon_reg) * deg2rad * ABS(ecci_grid%dlat_reg)
    ! area of ecci raw data pixel at equator

    hp   = 30.0      ! height of Prandtl-layer
    lnhp = LOG(hp)

    default_real = 0.0
    undefined_integer= NINT(undefined)

    ecci_class_fraction = default_real
    ecci_class_npixel   = undefined_integer
    ecci_tot_npixel = undefined_integer

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
    CALL init_ecci_lookup_tables(nclass_ecci, &
         &                       ilookup_table_ecci, &
         &                       z0_lt_ecci,            &
         &                       lnz0_lt_ecci,          &
         &                       plc_mn_lt_ecci,        &
         &                       plc_mx_lt_ecci,        &
         &                       lai_mn_lt_ecci,        &
         &                       lai_mx_lt_ecci,        &
         &                       rd_lt_ecci,            &
         &                       skinc_lt_ecci,         &
         &                       emiss_lt_ecci,         &
         &                       rs_min_lt_ecci)

    CALL get_name_ecci_lookup_tables(ilookup_table_ecci, name_lookup_table_ecci)

    nlon = ecci_grid%nlon_reg
    ALLOCATE(ie_vec(nlon),je_vec(nlon),ke_vec(nlon))
    ie_vec(:) = 0
    je_vec(:) = 0
    ke_vec(:) = 0
    start_cell_id = 1

    ! open netcdf file
    DO nt = 1,ntiles_ecci
      CALL check_netcdf( nf90_open(TRIM(ecci_file(nt)),NF90_NOWRITE, ncid_ecci(nt)))
    END DO

    varname = 'lccs_class' ! I know that the ecci data are stored in a variable called 'lccs_class'

    CALL check_netcdf(nf90_inq_varid(ncid_ecci(1), TRIM(varname), varid_ecci))

    mlat = 1
    block_row_start = mlat

    CALL det_band_ecci_data(ecci_grid,block_row_start,ta_grid)

    IF (ALLOCATED(lu_block)) THEN
      DEALLOCATE(lu_block, STAT = errorcode)
      IF(errorcode /= 0) CALL logging%error('Cant deallocate the lu_block',__FILE__,__LINE__)
    END IF
    ALLOCATE (lu_block(1:ta_grid%nlon_reg, 1:ta_grid%nlat_reg), STAT = errorcode)
    IF (errorcode /= 0) CALL logging%error('Cant allocate the lu_block',__FILE__,__LINE__)

    CALL get_ecci_data_block(ta_grid,               &
         ecci_tiles_grid,  &
         ncid_ecci,       &
         lu_block)

    block_row = 0

    ! Determine start and end longitude of search
    istartlon = 1
    iendlon = ecci_grid%nlon_reg
    IF (tg%igrid_type == igrid_icon) THEN
      DO i_col = 1, ecci_grid%nlon_reg
        point_lon = lon_ecci(i_col)
        IF (point_lon < tg%minlon) istartlon = i_col + 1
        IF (point_lon > tg%maxlon) THEN
          iendlon = i_col - 1
          EXIT
        ENDIF
      ENDDO
    ELSE IF (tg%igrid_type == igrid_cosmo) THEN
      DO i_col = 1, ecci_grid%nlon_reg
        point_lon = lon_ecci(i_col)
        IF (point_lon < bound_west_cosmo) istartlon = i_col + 1
        IF (point_lon > bound_east_cosmo) THEN
          iendlon = i_col - 1
          EXIT
        ENDIF
      ENDDO
    ENDIF
    nlon_sub = iendlon - istartlon + 1

    num_blocks = 1
    IF (MOD(nlon_sub,num_blocks)== 0) THEN
      blk_len = nlon_sub/num_blocks
    ELSE
      blk_len = nlon_sub/num_blocks + 1
    ENDIF

    ecci_rows: DO mlat = 1,ecci_grid%nlat_reg

      block_row= block_row + 1
      IF(block_row > ta_grid%nlat_reg) THEN ! read in new block
        block_row_start = mlat
        block_row = 1
        CALL det_band_ecci_data(ecci_grid,block_row_start,ta_grid)
        IF(ALLOCATED(lu_block)) THEN
          DEALLOCATE(lu_block, STAT=errorcode)
          IF(errorcode/=0) CALL logging%error('Cant deallocate the lu_block',__FILE__,__LINE__)
        ENDIF
        ALLOCATE (lu_block(1:ta_grid%nlon_reg,1:ta_grid%nlat_reg), STAT=errorcode)
        IF(errorcode/=0) CALL logging%error('Cant allocate lu_block',__FILE__,__LINE__)
        CALL get_ecci_data_block(ta_grid,              &
             &                        ecci_tiles_grid, &
             &                        ncid_ecci,       &
             &                        lu_block)
      ENDIF

      point_lat = lat_ecci(mlat)

      IF (tg%igrid_type == igrid_cosmo) THEN ! CASE COSMO grid, omit some I/O from hard disk when out of target domain
        IF ((point_lat > bound_north_cosmo).OR.(point_lat < bound_south_cosmo) ) THEN ! raw data out of target grid
          CYCLE ecci_rows
        ENDIF
      ELSE IF (tg%igrid_type == igrid_icon) THEN
        IF (point_lat > tg%maxlat .OR. point_lat < tg%minlat) THEN
          CYCLE ecci_rows
        ENDIF
      ENDIF ! grid type

      ecci_data_row(1:nlon) = lu_block(1:nlon,block_row)
      apix = apix_e * COS(point_lat * deg2rad) ! area of raw data pixel (in [m**2])

      ie_vec(istartlon:iendlon) = 0
      IF (tg%igrid_type /= igrid_icon) THEN
        je_vec(:) = 0
        ke_vec(:) = 0
      ENDIF
      DO ib = 1, num_blocks

        ishift = istartlon-1+(ib-1)*blk_len

        columns1: DO il = 1, blk_len
          i_col = ishift+il
          IF (i_col > iendlon) CYCLE columns1

          ! find the corresponding target grid indices
          point_lon = lon_ecci(i_col)

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
      ENDDO
      
      columns2: DO i_col = istartlon, iendlon

        ! find the corresponding target grid indices
        ie = ie_vec(i_col)
        je = je_vec(i_col)
        ke = ke_vec(i_col)

        IF ((ie /= 0).AND.(je /= 0).AND.(ke /= 0)) THEN
          ! raw data pixel within target grid, see output of routine find_rotated_lonlat_grid_element_index

          lu = ecci_data_row(i_col)                        ! land use class

          CALL ecci_look_up(lu,    &
               &      nclass_ecci,   &
               &      lnz0_lt_ecci,  &
               &      plc_mn_lt_ecci,&
               &      plc_mx_lt_ecci,&
               &      lai_mn_lt_ecci,&
               &      lai_mx_lt_ecci,&
               &      rd_lt_ecci,    &
               &      skinc_lt_ecci, &
               &      emiss_lt_ecci, &
               &      rs_min_lt_ecci,&
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


          IF (k_error == 0) THEN ! valid land use class

            ecci_tot_npixel(ie,je,ke) = ecci_tot_npixel(ie,je,ke) + 1
            a_weight(ie,je,ke) = a_weight(ie,je,ke) + apix  ! sum up area for weight

            CALL get_ecci_idx(lu,nclass)
            ecci_class_npixel(ie,je,ke,nclass) = ecci_class_npixel(ie,je,ke,nclass) + 1
            a_class(ie,je,ke,nclass) = a_class(ie,je,ke,nclass) + apix   ! sum area of valid land use pixels
            ! (use as weight later)
            emissivity_ecci(ie,je,ke) =  emissivity_ecci(ie,je,ke) + apix * pemissivity
            IF (pland >  0.0) THEN ! only for land pixel

              ! weighted with whole area
              fr_land_ecci(ie,je,ke) = fr_land_ecci(ie,je,ke) + apix * pland
              ice_ecci(ie,je,ke) = ice_ecci(ie,je,ke) + apix * pice
              urban_ecci(ie,je,ke) = urban_ecci(ie,je,ke) + apix * purb
              ! z0 is averaged logarithmic
              IF ( lnhp /= plnz0) THEN ! z0 is averaged logarithmic
                pwz0 = 1./(lnhp - plnz0)
              ELSE
                pwz0 = 0.
              ENDIF
              z0_ecci(ie,je,ke)      = z0_ecci(ie,je,ke) + apix * pwz0
              plcov_mn_ecci(ie,je,ke) = plcov_mn_ecci(ie,je,ke) + apix * pmn
              plcov_mx_ecci(ie,je,ke) = plcov_mx_ecci(ie,je,ke) + apix * pmx
              skinc_ecci(ie,je,ke)    = skinc_ecci(ie,je,ke)    + apix * pskinc
              ! the following fields are weighted with the plant cover
              root_ecci(ie,je,ke) = root_ecci(ie,je,ke) + apix * pmx * proot
              lai_mn_ecci(ie,je,ke) = lai_mn_ecci(ie,je,ke) + apix * pmx * plaimn
              lai_mx_ecci(ie,je,ke) = lai_mx_ecci(ie,je,ke) + apix * pmx * plaimx
              rs_min_ecci(ie,je,ke) = rs_min_ecci(ie,je,ke) + apix * pmx * prs_min
              for_d_ecci(ie,je,ke) = for_d_ecci(ie,je,ke) + apix * pmx * pfor_d
              for_e_ecci(ie,je,ke) = for_e_ecci(ie,je,ke) + apix * pmx * pfor_e

            END IF

          ENDIF
        ENDIF

        ! end loops
      ENDDO columns2
    ENDDO ecci_rows

    DEALLOCATE(ie_vec,je_vec,ke_vec)

    ! calculate ecci_class_fraction (ecci_class_fraction/ecci_class_npixel)
    DO ke=1, tg%ke
      DO je=1, tg%je
        DO ie=1, tg%ie
          area_tot = a_weight(ie,je,ke)
          area_land = fr_land_ecci(ie,je,ke)

          ! weight by total area
          IF (area_tot > 0.0) THEN
            emissivity_ecci(ie,je,ke) = emissivity_ecci(ie,je,ke) / area_tot
            DO l=1,nclass_ecci
              ecci_class_fraction(ie,je,ke,l) =  a_class(ie,je,ke,l) /  area_tot
              ! area fraction of each land use class
            ENDDO
          ENDIF

          IF (area_land > 0.0 ) THEN
            area_plcov = plcov_mx_ecci(ie,je,ke) ! area covered with plants

            z0_ecci(ie,je,ke) = z0_ecci(ie,je,ke) / area_land
            z0_ecci(ie,je,ke) = hp * EXP(-1./z0_ecci(ie,je,ke))

            ! weight by total area
            fr_land_ecci(ie,je,ke) = fr_land_ecci(ie,je,ke) / area_tot
            ice_ecci(ie,je,ke)   = ice_ecci(ie,je,ke)   / area_tot


            ! weight by land area
            urban_ecci(ie,je,ke) = urban_ecci(ie,je,ke) / area_land
            for_d_ecci(ie,je,ke) = for_d_ecci(ie,je,ke) / area_land
            for_e_ecci(ie,je,ke) = for_e_ecci(ie,je,ke) / area_land
            plcov_mn_ecci(ie,je,ke) = plcov_mn_ecci(ie,je,ke) / area_land

            plcov_mx_ecci(ie,je,ke) = plcov_mx_ecci(ie,je,ke) / area_land
            skinc_ecci(ie,je,ke) = skinc_ecci(ie,je,ke) / area_land
            ! weight by area covered with plants
            IF (area_plcov > 0.0) THEN
              root_ecci(ie,je,ke) = root_ecci(ie,je,ke)     / area_plcov
              lai_mn_ecci(ie,je,ke) = lai_mn_ecci(ie,je,ke) / area_plcov
              lai_mx_ecci(ie,je,ke) = lai_mx_ecci(ie,je,ke) / area_plcov
              rs_min_ecci(ie,je,ke) = rs_min_ecci(ie,je,ke) / area_plcov
            ELSE ! may occur for ice grid elements
              root_ecci(ie,je,ke) = undefined
              lai_mn_ecci(ie,je,ke) = undefined
              lai_mx_ecci(ie,je,ke) = undefined
              rs_min_ecci(ie,je,ke) = rs_min_undef
            ENDIF

          ELSE IF ( area_tot > 0.0) THEN ! only sea pixels were found
            fr_land_ecci(ie,je,ke) = undefined
            z0_ecci(ie,je,ke)      = undefined
            ice_ecci(ie,je,ke)     = undefined
            urban_ecci(ie,je,ke)   = undefined
            for_d_ecci(ie,je,ke)   = undefined
            for_e_ecci(ie,je,ke)   = undefined
            plcov_mx_ecci(ie,je,ke)= undefined
            plcov_mn_ecci(ie,je,ke)= undefined

            root_ecci(ie,je,ke)    = undefined
            lai_mx_ecci(ie,je,ke)  = undefined
            lai_mn_ecci(ie,je,ke)  = undefined
            skinc_ecci(ie,je,ke)   = undefined
            rs_min_ecci(ie,je,ke)  = undefined

          ENDIF
        ENDDO
      ENDDO
    ENDDO

    l_opn_gc_file = .FALSE.

    DO ke=1, tg%ke
      DO je=1, tg%je
        DO ie=1, tg%ie
          IF (ecci_tot_npixel(ie,je,ke) == 0) THEN ! nearest neighbor search

            point_lon = lon_geo(ie,je,ke)
            point_lat = lat_geo(ie,je,ke)

            CALL find_reg_lonlat_grid_element_index(point_lon,      &
                 &                                     point_lat,      &
                 &                                     ecci_grid, &
                 &                                     i_lu,           &
                 &                                     j_lu,           &
                 &                                     ntiles_ecci,&
                 &                                     tile = tile,    &
                 &                                     regular_tiles_grid_info=ecci_tiles_grid)

            IF ((i_lu /= 0).AND.(j_lu /= 0))THEN

              i_col = i_lu
              j_row = j_lu

              IF (.NOT. l_opn_gc_file(tile)) THEN
                CALL check_netcdf(nf90_open(TRIM(ecci_file(tile)),NF90_NOWRITE, ncid_ecci(tile)), __FILE__, __LINE__)
                CALL check_netcdf(nf90_inq_varid(ncid_ecci(tile),TRIM(varname), varid_gc(tile)))
                l_opn_gc_file(tile) = .TRUE.
              ENDIF

              CALL check_netcdf(nf90_get_var(ncid_ecci(tile), varid_gc(tile),  &
                   &               ecci_data_pixel,  &
                   &               start=(/ i_col,j_row /),count=(/ 1,1 /)), __FILE__, __LINE__)

              lu = ecci_data_pixel(1,1)


              CALL ecci_look_up(lu, &
                   &      nclass_ecci, &
                   &      lnz0_lt_ecci,          &
                   &      plc_mn_lt_ecci,        &
                   &      plc_mx_lt_ecci,        &
                   &      lai_mn_lt_ecci,        &
                   &      lai_mx_lt_ecci,        &
                   &      rd_lt_ecci,            &
                   &      skinc_lt_ecci,         &
                   &      emiss_lt_ecci,         &
                   &      rs_min_lt_ecci,        &
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
              lu = 0
              k_error = 1
            ENDIF

            IF (k_error == 0) THEN ! valid land use class
              apix = apix_e * COS(point_lat * deg2rad) ! area of raw data pixel (in [m**2])
              CALL get_ecci_idx(lu,nclass)
              ecci_class_npixel(ie,je,ke,nclass) = ecci_class_npixel(ie,je,ke,nclass) + 1
              a_class(ie,je,ke,nclass) = a_class(ie,je,ke,nclass) + apix   ! sum area of valid land use pixels
              emissivity_ecci(ie,je,ke) =  pemissivity
              IF (pland >  0.0) THEN ! only for land pixel
                fr_land_ecci(ie,je,ke) = pland
                ice_ecci(ie,je,ke) =  pice
                urban_ecci(ie,je,ke) = purb
                IF ( lnhp /= plnz0) THEN ! log z0
                  pwz0 = 1./(lnhp - plnz0)
                ELSE
                  pwz0 = 0.
                ENDIF
                z0_ecci(ie,je,ke)      =  pwz0
                z0_ecci(ie,je,ke) = hp * EXP(-1./z0_ecci(ie,je,ke))
                plcov_mn_ecci(ie,je,ke) = pmn
                plcov_mx_ecci(ie,je,ke) = pmx
                skinc_ecci(ie,je,ke)    = pskinc 
                root_ecci(ie,je,ke) = proot
                lai_mn_ecci(ie,je,ke) = plaimn
                lai_mx_ecci(ie,je,ke) = plaimx
                rs_min_ecci(ie,je,ke) = prs_min
                for_d_ecci(ie,je,ke) = pfor_d * pmx
                for_e_ecci(ie,je,ke) = pfor_e * pmx
                IF (pice == 1.0) THEN  ! ice land use class
                  root_ecci(ie,je,ke) = undefined
                  lai_mn_ecci(ie,je,ke) = undefined
                  lai_mx_ecci(ie,je,ke) = undefined
                  rs_min_ecci(ie,je,ke) = rs_min_undef
                ENDIF
              ENDIF
            ELSE ! not a valid land use class
              fr_land_ecci(ie,je,ke) = undefined
              z0_ecci(ie,je,ke)      = undefined
              ice_ecci(ie,je,ke)     = undefined
              urban_ecci(ie,je,ke)   = undefined
              for_d_ecci(ie,je,ke)   = undefined
              for_e_ecci(ie,je,ke)   = undefined
              plcov_mx_ecci(ie,je,ke)= undefined
              plcov_mn_ecci(ie,je,ke)= undefined
              root_ecci(ie,je,ke)    = undefined
              lai_mx_ecci(ie,je,ke)  = undefined
              lai_mn_ecci(ie,je,ke)  = undefined
              skinc_ecci(ie,je,ke)   = undefined
              rs_min_ecci(ie,je,ke)  = undefined
            ENDIF
          ENDIF ! nearest neighbour search
        ENDDO
      ENDDO
    ENDDO

    DO tile=1,ntiles_ecci
      IF (l_opn_gc_file(tile)) THEN
        CALL check_netcdf(nf90_close(ncid_ecci(tile)))
      ENDIF
    ENDDO

  END SUBROUTINE agg_ecci_data_to_target_grid

END MODULE mo_agg_ecci
