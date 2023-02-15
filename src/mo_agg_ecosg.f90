!=======================================================================
!+ Fortran module to aggregate ECOSG land use data to a target grid 
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_0         2022/12/21 Andrzej Wyszogrodzki
!
!   Parallel threads for ICON and COSMO using Open-MP, 
!   Several bug fixes and optimizations for ICON search algorithm, 
!   particularly for the special case of non-contiguous domains; 
!   simplified namelist control for ICON  
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!> Fortran module to aggregate ECOSG land use data to a target grid 
!!
!> \author Andrzej Wyszogrodzki
!=======================================================================
MODULE mo_agg_ecosg

  USE mo_logging
  USE mo_kind,                  ONLY: wp, i4

  USE mo_io_units,              ONLY: filename_max

  USE mo_grid_structures,       ONLY: igrid_icon, &
       &                              igrid_cosmo, & 
       &                              target_grid_def
                                
  USE mo_search_ll_grid,        ONLY: find_reg_lonlat_grid_element_index
                                
  USE mo_io_utilities,          ONLY: check_netcdf

  USE mo_search_target_grid,    ONLY: find_nearest_target_grid_element

  USE netcdf,                   ONLY :   &
       &                              nf90_open,              &
       &                              nf90_close,             &
       &                              nf90_inq_varid,         &
       &                              nf90_get_var,           &
       &                              NF90_NOWRITE

  USE mo_ecosg_data,            ONLY: ecosg_grid, &
       &                              lon_ecosg,  &
       &                              lat_ecosg

  USE mo_ecosg_lookup_tables,    ONLY: name_lookup_table_ecosg, & 
       &                              init_ecosg_lookup_tables, &
       &                              get_name_ecosg_lookup_tables, & 
       &                              z0_lt_ecosg, &
       &                              lnz0_lt_ecosg, &
       &                              plc_mn_lt_ecosg, &
       &                              plc_mx_lt_ecosg, & 
       &                              lai_mn_lt_ecosg, &
       &                              lai_mx_lt_ecosg, &
       &                              rd_lt_ecosg, &
       &                              skinc_lt_ecosg, &
       &                              emiss_lt_ecosg, &
       &                              rs_min_lt_ecosg   

  USE mo_ecosg_lookup_tables,    ONLY: ecosg_look_up


  USE mo_math_constants,        ONLY: deg2rad
  USE mo_physical_constants,    ONLY: re

  USE mo_target_grid_data,      ONLY: lon_geo, & 
       &                              lat_geo, & 
       &                              search_res

  USE mo_utilities_extpar,      ONLY: check_input_file
  
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: agg_ecosg_data_to_target_grid

  REAL(KIND=wp), PARAMETER :: rs_min_undef=999. !< undefined value for minimal stomata resistance

  CONTAINS

  !> Subroutine to aggregate ecosg data to the target grid
  SUBROUTINE agg_ecosg_data_to_target_grid(ecosg_file,ilookup_table_ecosg, &
       &                                        undefined,            &
       &                                        tg,                   &
       &                                        nclass_ecosg,         &
       &                                        ecosg_class_fraction, & 
       &                                        ecosg_class_npixel,   &
       &                                        ecosg_tot_npixel,     &
       &                                        fr_land_ecosg ,       &  
       &                                        ice_ecosg,            &
       &                                        z0_ecosg,             &
       &                                        root_ecosg,           &
       &                                        plcov_mn_ecosg,       &
       &                                        plcov_mx_ecosg,       &
       &                                        lai_mn_ecosg,         &
       &                                        lai_mx_ecosg,         &
       &                                        rs_min_ecosg,         &
       &                                        urban_ecosg,          &
       &                                        for_d_ecosg,          &
       &                                        for_e_ecosg,          &
       &                                        skinc_ecosg,          &
       &                                        emissivity_ecosg    )    

    TYPE(target_grid_def), INTENT(IN)        :: tg  !< structure with target grid description

    CHARACTER (LEN=filename_max), INTENT(IN) :: ecosg_file  !< filename ecosg raw data

    INTEGER(KIND=i4), INTENT(IN)             :: ilookup_table_ecosg, & 
         &                                      nclass_ecosg

    REAL (KIND=wp), INTENT(IN)               :: undefined            !< undef value


    INTEGER (KIND=i4), INTENT(OUT)           :: ecosg_class_npixel(:,:,:,:),& 
         &                                      ecosg_tot_npixel(:,:,:)  

    REAL (KIND=wp), INTENT(OUT)              :: ecosg_class_fraction(:,:,:,:), & 
         &                                      fr_land_ecosg(:,:,:), &  !< fraction land due to ecosg raw data
         &                                      ice_ecosg(:,:,:), &      !< fraction of ice due to ecosg raw data
         &                                      z0_ecosg(:,:,:), &       !< roughness length due to ecosg land use data
         &                                      root_ecosg(:,:,:), &     !< root depth due to ecosg land use data
         &                                      plcov_mx_ecosg(:,:,:), & !< plant cover maximum due to ecosg land use data
         &                                      plcov_mn_ecosg(:,:,:), & !< plant cover minimum due to ecosg land use data
         &                                      lai_mx_ecosg(:,:,:), &   !< Leaf Area Index maximum due to ecosg land use data
         &                                      lai_mn_ecosg(:,:,:), &   !< Leaf Area Index minimum due to ecosg land use data
         &                                      rs_min_ecosg(:,:,:), &   !< minimal stomata resistance due to ecosg land use data
         &                                      urban_ecosg(:,:,:), &    !< urban fraction due to ecosg land use data
         &                                      for_d_ecosg(:,:,:), &    !< deciduous forest (fraction) due to ecosg land use data
         &                                      for_e_ecosg(:,:,:), &    !< evergreen forest (fraction) due to ecosg land use data
         &                                      skinc_ecosg(:,:,:), &    !< skin conductivity due to ecosg land use data
         &                                      emissivity_ecosg(:,:,:)  !< longwave emissivity due to ecosg land use da

    !local variables
    REAL (KIND=wp)                          :: default_real, & 
         &                                     a_weight(1:tg%ie,1:tg%je,1:tg%ke), &  !< area weight of all raw data pixels in target grid
         &                                     a_class(1:tg%ie,1:tg%je,1:tg%ke,1:nclass_ecosg), &  
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
         &                                     hp, &              ! height of Prandtl-layer
         &                                     lnhp, & 
         &                                     pwz0, &  ! weighted summand for z0
         &                                     area_tot, &    ! total area
         &                                     area_land, &   ! area with land
         &                                     area_plcov, &  ! area covered with plants
         &                                     bound_north_cosmo, &  !< northern boundary for COSMO target domain
         &                                     bound_south_cosmo, &  !< southern boundary for COSMO target domain
         &                                     bound_west_cosmo, &   !< western  boundary for COSMO target domain
         &                                     bound_east_cosmo  !< eastern  boundary for COSMO target domain

    INTEGER (KIND=i4)                       :: undefined_integer, &  ! undef value
         &                                     l, &  ! counters
         &                                     i_col, j_row, &  ! counter
         &                                     i_lu, j_lu, & 
         &                                     ie, je, ke, &   ! indices for target grid elements
         &                                     start_cell_id, &  !< ID of starting cell for ICON search
         &                                     i1, i2, & 
         &                                     ecosg_data_row(ecosg_grid%nlon_reg), & 
         &                                     ecosg_data_pixel(1:1,1:1), & 
         &                                     lu, &   ! land use class
         &                                     ncid_ecosg, &                              !< netcdf unit file number
         &                                     varid_ecosg, &                !< id of variable
         &                                     nlon, & 
         &                                     k_error      ! error return code
         
    INTEGER (KIND=i4), ALLOCATABLE          :: ie_vec(:), je_vec(:), ke_vec(:)  ! indices for target grid elements


    CHARACTER (LEN=80)                      :: varname  !< name of variable
    ! Some stuff for OpenMP parallelization
    INTEGER(KIND=i4)                        :: num_blocks, ib, il, blk_len, istartlon, iendlon, nlon_sub, ishift
    !$   INTEGER :: omp_get_max_threads, omp_get_thread_num, thread_id
    !$   INTEGER (i4), ALLOCATABLE :: start_cell_arr(:)

    CALL logging%info('Enter routine: agg_ecosg_data_to_target_grid')

    fr_land_ecosg(:,:,:)    = 0.0_wp
    ice_ecosg(:,:,:)        = 0.0_wp
    z0_ecosg(:,:,:)         = 0.0_wp
    root_ecosg(:,:,:)       = 0.0_wp
    plcov_mx_ecosg(:,:,:)   = 0.0_wp
    plcov_mn_ecosg(:,:,:)   = 0.0_wp
    lai_mx_ecosg(:,:,:)     = 0.0_wp
    lai_mn_ecosg(:,:,:)     = 0.0_wp
    rs_min_ecosg(:,:,:)     = 0.0_wp
    urban_ecosg(:,:,:)      = 0.0_wp
    for_d_ecosg(:,:,:)      = 0.0_wp
    for_e_ecosg(:,:,:)      = 0.0_wp
    skinc_ecosg(:,:,:)      = 0.0_wp
    emissivity_ecosg(:,:,:) = 0.0_wp

    apix_e  = re * re * deg2rad* ABS(ecosg_grid%dlon_reg) * deg2rad * ABS(ecosg_grid%dlat_reg) 
    ! area of ECOSG raw data pixel at equator

    hp   = 30.0      ! height of Prandtl-layer
    lnhp = LOG(hp)

    default_real = 0.0
    undefined_integer= NINT(undefined)

    ecosg_class_fraction = default_real
    ecosg_class_npixel   = undefined_integer
    ecosg_tot_npixel = undefined_integer

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

    emissivity_ecosg(:,:,:) = 0.0_wp
    
    ! init lookup tables
    CALL init_ecosg_lookup_tables(nclass_ecosg, &
         &      ilookup_table_ecosg, &
         &      z0_lt_ecosg,            &
         &      lnz0_lt_ecosg,          &
         &      plc_mn_lt_ecosg,        &
         &      plc_mx_lt_ecosg,        &
         &      lai_mn_lt_ecosg,        &
         &      lai_mx_lt_ecosg,        &
         &      rd_lt_ecosg,            &
         &      skinc_lt_ecosg,         &
         &      emiss_lt_ecosg,         &
         &      rs_min_lt_ecosg)
 
    CALL get_name_ecosg_lookup_tables(ilookup_table_ecosg, name_lookup_table_ecosg)  

    ! open netcdf file
    CALL check_input_file(TRIM(ecosg_file),__FILE__,__LINE__)
    CALL check_netcdf( nf90_open(TRIM(ecosg_file),NF90_NOWRITE, ncid_ecosg),__FILE__,__LINE__)

    varname = 'lc_class'

    CALL check_netcdf( nf90_inq_varid(ncid_ecosg, TRIM(varname), varid_ecosg),__FILE__,__LINE__)
    nlon = ecosg_grid%nlon_reg
    ALLOCATE(ie_vec(nlon),je_vec(nlon),ke_vec(nlon))
    ie_vec(:) = 0
    je_vec(:) = 0
    ke_vec(:) = 0
    start_cell_id = 1

    ! Determine start and end longitude of search
    istartlon = 1
    iendlon = ecosg_grid%nlon_reg
    IF (tg%igrid_type == igrid_icon) THEN
      DO i_col = 1, ecosg_grid%nlon_reg
        point_lon = lon_ecosg(i_col)
        IF (point_lon < tg%minlon) istartlon = i_col + 1
        IF (point_lon > tg%maxlon) THEN
          iendlon = i_col - 1
          EXIT
        ENDIF
      ENDDO
    ELSE IF (tg%igrid_type == igrid_cosmo) THEN
      DO i_col = 1, ecosg_grid%nlon_reg
        point_lon = lon_ecosg(i_col)
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
    WRITE(message_text,*) 'nlon_sub: ',nlon_sub,' num_blocks: ',num_blocks, ' blk_len: ',blk_len
    CALL logging%info(message_text)

    CALL logging%info('Start loop over ecosg rows...')

    ! loop over rows of ECOSG dataset
    rows: DO j_row=1,ecosg_grid%nlat_reg
      point_lat = lat_ecosg(j_row)

      IF (tg%igrid_type == igrid_cosmo) THEN ! CASE COSMO grid, omit I/O from hard disk when outside target domain
        IF ((point_lat > bound_north_cosmo).OR.(point_lat < bound_south_cosmo) ) THEN ! raw data out of target grid
          CYCLE rows
        ENDIF
      ELSE IF (tg%igrid_type == igrid_icon) THEN 
        IF (point_lat > tg%maxlat .OR. point_lat < tg%minlat) THEN
          CYCLE rows
        ENDIF
      ENDIF ! grid type

      ! read in pixels
      CALL check_netcdf(nf90_get_var(ncid_ecosg, varid_ecosg,  ecosg_data_row,  &
           &               start=(/1,j_row/),count=(/nlon,1/)),__FILE__,__LINE__)
      apix = apix_e * COS(point_lat * deg2rad) ! area of raw data pixel (in [m**2])
      ie_vec(istartlon:iendlon) = 0
      IF (tg%igrid_type /= igrid_icon) THEN
        je_vec(:) = 0
        ke_vec(:) = 0
      ENDIF

!$OMP PARALLEL DO PRIVATE(ib,il,i_col,i1,i2,ishift,point_lon,thread_id,start_cell_id)
      DO ib = 1, num_blocks

        !$     thread_id = omp_get_thread_num()+1
        !$     start_cell_id = start_cell_arr(thread_id)
        ishift = istartlon-1+(ib-1)*blk_len

        columns1: DO il = 1,blk_len
          i_col = ishift+il
          IF (i_col > iendlon) CYCLE columns1

          ! find the corresponding target grid indices
          point_lon = lon_ecosg(i_col)

          ! Reset start cell when entering a new row or when the previous data point was outside
          ! the model domain
          IF (tg%igrid_type == igrid_icon .AND. (il == 1 .OR. start_cell_id == 0)) THEN
            i1 = NINT(point_lon*search_res)
            i2 = NINT(point_lat*search_res)
            start_cell_id = tg%search_index(i1,i2)
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
!$OMP END PARALLEL DO

      columns2: DO i_col=istartlon,iendlon

        ie = ie_vec(i_col)
        je = je_vec(i_col)
        ke = ke_vec(i_col)

        IF ((ie /= 0).AND.(je/=0).AND.(ke/=0))THEN 
          lu = ecosg_data_row(i_col)                        ! land use class

          CALL ecosg_look_up(lu, &
               &      nclass_ecosg, &
               &      lnz0_lt_ecosg,          &
               &      plc_mn_lt_ecosg,        &
               &      plc_mx_lt_ecosg,        &
               &      lai_mn_lt_ecosg,        &
               &      lai_mx_lt_ecosg,        &
               &      rd_lt_ecosg,            &
               &      skinc_lt_ecosg,         &
               &      emiss_lt_ecosg,         &
               &      rs_min_lt_ecosg,        &
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

          IF (k_error == 0) THEN ! valid land use class

            ! raw data pixel within target grid, see output of routine find_rotated_lonlat_grid_element_index
            !- summation of variables
            ecosg_tot_npixel(ie,je,ke) = ecosg_tot_npixel(ie,je,ke) + 1
            a_weight(ie,je,ke) = a_weight(ie,je,ke) + apix  ! sum up area for weight

            ecosg_class_npixel(ie,je,ke,lu) = ecosg_class_npixel(ie,je,ke,lu) + 1
            a_class(ie,je,ke,lu) = a_class(ie,je,ke,lu) + apix   ! sum area of valid land use pixels 
            !(use as weight later)
            emissivity_ecosg(ie,je,ke) =  emissivity_ecosg(ie,je,ke) + apix * pemissivity
            IF (pland >  0.0) THEN ! only for land pixel

              ! weighted with whole area
              fr_land_ecosg(ie,je,ke) = fr_land_ecosg(ie,je,ke) + apix * pland
              ice_ecosg(ie,je,ke) = ice_ecosg(ie,je,ke) + apix * pice
              urban_ecosg(ie,je,ke) = urban_ecosg(ie,je,ke) + apix * purb
              ! z0 is averaged logarithmic
              IF ( lnhp /= plnz0) THEN ! z0 is averaged logarithmic
                pwz0 = 1./(lnhp - plnz0)
              ELSE
                pwz0 = 0.
              ENDIF
              z0_ecosg(ie,je,ke)      = z0_ecosg(ie,je,ke) + apix * pwz0
              plcov_mn_ecosg(ie,je,ke) = plcov_mn_ecosg(ie,je,ke) + apix * pmn
              plcov_mx_ecosg(ie,je,ke) = plcov_mx_ecosg(ie,je,ke) + apix * pmx
              skinc_ecosg(ie,je,ke)    = skinc_ecosg(ie,je,ke)    + apix * pskinc
              ! the following fields are weighted with the plant cover
              root_ecosg(ie,je,ke) = root_ecosg(ie,je,ke) + apix * pmx * proot
              lai_mn_ecosg(ie,je,ke) = lai_mn_ecosg(ie,je,ke) + apix * pmx * plaimn
              lai_mx_ecosg(ie,je,ke) = lai_mx_ecosg(ie,je,ke) + apix * pmx * plaimx
              rs_min_ecosg(ie,je,ke) = rs_min_ecosg(ie,je,ke) + apix * pmx * prs_min
              for_d_ecosg(ie,je,ke) = for_d_ecosg(ie,je,ke) + apix * pmx * pfor_d
              for_e_ecosg(ie,je,ke) = for_e_ecosg(ie,je,ke) + apix * pmx * pfor_e

            END IF
          ENDIF
        ENDIF

        ! end loops
      ENDDO columns2
    ENDDO rows

    DEALLOCATE(ie_vec,je_vec,ke_vec)
    !$   DEALLOCATE(start_cell_arr)
    CALL logging%info('...done')

    ! calculate ecosg_class_fraction (ecosg_class_fraction/ecosg_class_npixel)
    DO ke=1, tg%ke
      DO je=1, tg%je
        DO ie=1, tg%ie
          area_tot = a_weight(ie,je,ke)
          area_land = fr_land_ecosg(ie,je,ke)

          ! weight by total area
          IF (area_tot > 0.0) THEN
            emissivity_ecosg(ie,je,ke) = emissivity_ecosg(ie,je,ke) / area_tot
            DO l=1,nclass_ecosg
              ecosg_class_fraction(ie,je,ke,l) =  a_class(ie,je,ke,l) /  area_tot 
              ! area fraction of each land use class
            ENDDO
          ENDIF

          IF (area_land > 0.0 ) THEN
            area_plcov = plcov_mx_ecosg(ie,je,ke)
            !area_plcov = MAX(area_plcov, eps)      ! area covered with plants

            z0_ecosg(ie,je,ke) = z0_ecosg(ie,je,ke) / area_land
            z0_ecosg(ie,je,ke) = hp * EXP(-1./z0_ecosg(ie,je,ke))

            ! weight by total area
            fr_land_ecosg(ie,je,ke) = fr_land_ecosg(ie,je,ke) / area_tot
            ice_ecosg(ie,je,ke)   = ice_ecosg(ie,je,ke)   / area_tot


            ! weight by land area
            urban_ecosg(ie,je,ke) = urban_ecosg(ie,je,ke) / area_land
            for_d_ecosg(ie,je,ke) = for_d_ecosg(ie,je,ke) / area_land
            for_e_ecosg(ie,je,ke) = for_e_ecosg(ie,je,ke) / area_land
            plcov_mn_ecosg(ie,je,ke) = plcov_mn_ecosg(ie,je,ke) / area_land

            plcov_mx_ecosg(ie,je,ke) = plcov_mx_ecosg(ie,je,ke) / area_land
            skinc_ecosg(ie,je,ke) = skinc_ecosg(ie,je,ke) / area_land
            ! weight by area covered with plants

            IF (area_plcov > 0.0) THEN 
              root_ecosg(ie,je,ke) = root_ecosg(ie,je,ke)     / area_plcov
              lai_mn_ecosg(ie,je,ke) = lai_mn_ecosg(ie,je,ke) / area_plcov
              lai_mx_ecosg(ie,je,ke) = lai_mx_ecosg(ie,je,ke) / area_plcov
              rs_min_ecosg(ie,je,ke) = rs_min_ecosg(ie,je,ke) / area_plcov
            ELSE ! may occur for ice grid elements
              root_ecosg(ie,je,ke) = undefined
              lai_mn_ecosg(ie,je,ke) = undefined
              lai_mx_ecosg(ie,je,ke) = undefined
              rs_min_ecosg(ie,je,ke) = rs_min_undef
            ENDIF


          ELSE IF ( area_tot > 0.0) THEN ! only sea pixels were found
            fr_land_ecosg(ie,je,ke) = undefined
            z0_ecosg(ie,je,ke)      = undefined
            ice_ecosg(ie,je,ke)     = undefined
            urban_ecosg(ie,je,ke)   = undefined
            for_d_ecosg(ie,je,ke)   = undefined
            for_e_ecosg(ie,je,ke)   = undefined
            plcov_mx_ecosg(ie,je,ke)= undefined
            plcov_mn_ecosg(ie,je,ke)= undefined

            root_ecosg(ie,je,ke)    = undefined
            lai_mx_ecosg(ie,je,ke)  = undefined
            lai_mn_ecosg(ie,je,ke)  = undefined
            skinc_ecosg(ie,je,ke)   = undefined
            rs_min_ecosg(ie,je,ke)  = undefined

          ENDIF
        ENDDO
      ENDDO
    ENDDO

    DO ke=1, tg%ke
      DO je=1, tg%je
        DO ie=1, tg%ie
          IF (ecosg_tot_npixel(ie,je,ke) == 0) THEN ! nearest neighbor search

            point_lon = lon_geo(ie,je,ke)
            point_lat = lat_geo(ie,je,ke)
            CALL find_reg_lonlat_grid_element_index(point_lon,      &
                 &                                     point_lat,      &
                 &                                     ecosg_grid,  &
                 &                                     i_lu,    &
                 &                                     j_lu)
            ! get data
            i_col = i_lu
            j_row = j_lu
            CALL check_netcdf(nf90_get_var(ncid_ecosg, varid_ecosg,  ecosg_data_pixel,  &
                 &               start=(/ i_col,j_row /),count=(/ 1,1 /)),__FILE__,__LINE__)

            lu = ecosg_data_pixel(1,1)


            CALL ecosg_look_up(lu, &
                 &      nclass_ecosg, &
                 &      lnz0_lt_ecosg,          &
                 &      plc_mn_lt_ecosg,        &
                 &      plc_mx_lt_ecosg,        &
                 &      lai_mn_lt_ecosg,        &
                 &      lai_mx_lt_ecosg,        &
                 &      rd_lt_ecosg,            &
                 &      skinc_lt_ecosg,         &
                 &      emiss_lt_ecosg,         &
                 &      rs_min_lt_ecosg,        &
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


            IF (k_error == 0) THEN ! valid land use class
              apix = apix_e * COS(point_lat * deg2rad) ! area of raw data pixel (in [m**2])
              ecosg_class_npixel(ie,je,ke,lu) = ecosg_class_npixel(ie,je,ke,lu) + 1
              a_class(ie,je,ke,lu) = a_class(ie,je,ke,lu) + apix   ! sum area of valid land use pixels 
              emissivity_ecosg(ie,je,ke) =  pemissivity
              IF (pland >  0.0) THEN ! only for land pixel
                fr_land_ecosg(ie,je,ke) = pland
                ice_ecosg(ie,je,ke) =  pice
                urban_ecosg(ie,je,ke) = purb
                IF ( lnhp /= plnz0) THEN ! log z0 
                  pwz0 = 1./(lnhp - plnz0)
                ELSE
                  pwz0 = 0.
                ENDIF
                z0_ecosg(ie,je,ke)      =  pwz0
                z0_ecosg(ie,je,ke) = hp * EXP(-1./z0_ecosg(ie,je,ke))
                plcov_mn_ecosg(ie,je,ke) = pmn
                plcov_mx_ecosg(ie,je,ke) = pmx
                skinc_ecosg(ie,je,ke)    = pskinc
                root_ecosg(ie,je,ke) = proot
                lai_mn_ecosg(ie,je,ke) = plaimn
                lai_mx_ecosg(ie,je,ke) = plaimx
                rs_min_ecosg(ie,je,ke) = prs_min
                for_d_ecosg(ie,je,ke) = pfor_d * pmx
                for_e_ecosg(ie,je,ke) = pfor_e * pmx
                IF (pice == 1.0) THEN  ! ice land use class
                  root_ecosg(ie,je,ke) = undefined
                  lai_mn_ecosg(ie,je,ke) = undefined
                  lai_mx_ecosg(ie,je,ke) = undefined
                  rs_min_ecosg(ie,je,ke) = rs_min_undef
                ENDIF
              ENDIF
            ELSE ! not a valid land use class
              fr_land_ecosg(ie,je,ke) = undefined
              z0_ecosg(ie,je,ke)      = undefined
              ice_ecosg(ie,je,ke)     = undefined
              urban_ecosg(ie,je,ke)   = undefined
              for_d_ecosg(ie,je,ke)   = undefined
              for_e_ecosg(ie,je,ke)   = undefined
              plcov_mx_ecosg(ie,je,ke)= undefined
              plcov_mn_ecosg(ie,je,ke)= undefined
              root_ecosg(ie,je,ke)    = undefined
              lai_mx_ecosg(ie,je,ke)  = undefined
              lai_mn_ecosg(ie,je,ke)  = undefined
              skinc_ecosg(ie,je,ke)   = undefined
              rs_min_ecosg(ie,je,ke)  = undefined
            ENDIF
          ENDIF ! nearest neighbour search
        ENDDO
      ENDDO
    ENDDO

    ! close netcdf file 
    CALL check_netcdf( nf90_close(ncid_ecosg),__FILE__,__LINE__)

    CALL logging%info('Exit routine: agg_ecosg_data_to_target_grid')

  END SUBROUTINE agg_ecosg_data_to_target_grid

END MODULE mo_agg_ecosg
