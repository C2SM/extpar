!+ Fortran module to aggregate GLC2000 land use data to a target grid 
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
! V1_4         2011/04/21 Hermann Asensio
!  clean up
! V1_7         2013/01/25 Guenther Zaengl 
!   Parallel threads for ICON and COSMO using Open-MP, 
!   Several bug fixes and optimizations for ICON search algorithm, 
!   particularly for the special case of non-contiguous domains; 
!   simplified namelist control for ICON  
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!> Fortran module to aggregate GLC2000 land use data to a target grid 
!!
!> \author Hermann Asensio
MODULE mo_agg_glc2000

  USE mo_kind, ONLY: wp, i4
  USE mo_logging

  !> data type structures form module GRID_structures
  USE mo_grid_structures, ONLY: igrid_icon
  USE mo_grid_structures, ONLY: igrid_cosmo

  USE mo_search_ll_grid, ONLY: find_reg_lonlat_grid_element_index
  USE mo_io_units,       ONLY: filename_max
  USE mo_io_utilities,   ONLY: check_netcdf

  USE mo_search_target_grid, ONLY: find_nearest_target_grid_element


  USE netcdf,      ONLY :   &
    & nf90_open,              &
    & nf90_close,             &
    & nf90_inq_varid,         &
    & nf90_get_var,           &
    & NF90_NOWRITE

  

  IMPLICIT NONE

  PRIVATE
  

  PUBLIC :: agg_glc2000_data_to_target_grid

  REAL(KIND=wp), PARAMETER :: rs_min_undef=999. !< undefined value for minimal stomata resistance


  CONTAINS

  !> Subroutine to aggregate glc2000 data to the target grid
  SUBROUTINE agg_glc2000_data_to_target_grid(glc2000_file,ilookup_table_glc2000,undefined,       &
    &                                        tg,                                         &
    &                                        nclass_glc2000,                             &
    &                                        glc2000_class_fraction, & 
    &                                        glc2000_class_npixel, &
    &                                        glc2000_tot_npixel,   &
    &                                        fr_land_glc2000 ,     &  
    &                                        ice_glc2000,          &
    &                                        z0_glc2000, &
    &                                        root_glc2000, &
    &                                        plcov_mn_glc2000, &
    &                                        plcov_mx_glc2000, &
    &                                        lai_mn_glc2000,   &
    &                                        lai_mx_glc2000, &
    &                                        rs_min_glc2000, &
    &                                        urban_glc2000,  &
    &                                        for_d_glc2000,  &
    &                                        for_e_glc2000, &
    &                                        emissivity_glc2000    )    
    
    
    
    
    
  !-------------------------------------------------------------------------------------
  ! list of modules which are used as "input"
  USE mo_grid_structures, ONLY: target_grid_def   !< type definition of structure for tg

  ! raw data grid description, here for GLC2000 data
  USE mo_glc2000_data, ONLY: glc2000_grid, &
    &                          lon_glc2000,  &
    &                          lat_glc2000

  USE mo_glc2000_lookup_tables, ONLY: name_lookup_table_glc2000
  USE mo_glc2000_lookup_tables, ONLY: init_glc2000_lookup_tables, &
    &                                 get_name_glc2000_lookup_tables
  USE mo_glc2000_lookup_tables, ONLY:   z0_lt_glc2000, lnz0_lt_glc2000, plc_mn_lt_glc2000, plc_mx_lt_glc2000, & 
    &               lai_mn_lt_glc2000, lai_mx_lt_glc2000, rd_lt_glc2000, emiss_lt_glc2000, rs_min_lt_glc2000   


  USE mo_glc2000_lookup_tables, ONLY: glc2000_look_up

    USE mo_math_constants, ONLY: deg2rad
    USE mo_physical_constants, ONLY: re
    ! USE global data fields (coordinates)
    USE mo_target_grid_data, ONLY: lon_geo, & !< longitude coordinates of the COSMO grid in the geographical system 
      &                            lat_geo !< latitude coordinates of the COSMO grid in the geographical system
    USE mo_target_grid_data, ONLY: search_res !< resolution of ICON grid search index list

     CHARACTER (LEN=filename_max), INTENT(IN) :: glc2000_file(:)  !< filename glc2000 raw data
     INTEGER, INTENT(IN) :: ilookup_table_glc2000
     REAL (KIND=wp), INTENT(IN) :: undefined            !< undef value

     TYPE(target_grid_def), INTENT(IN) :: tg  !< structure with target grid description
     INTEGER, INTENT(IN) :: nclass_glc2000 !< GLC2000 has 23 classes for the land use description
     REAL (KIND=wp), INTENT(OUT)  :: glc2000_class_fraction(:,:,:,:)  !< fraction for each glc2000 class on target grid 

    INTEGER (KIND=i4), INTENT(OUT) :: glc2000_class_npixel(:,:,:,:) !< number of raw data pixels for each class on target grid


    INTEGER (KIND=i4), INTENT(OUT) :: glc2000_tot_npixel(:,:,:)  !< total number of glc2000 raw data pixels on target grid 


    REAL (KIND=wp), INTENT(OUT)  :: fr_land_glc2000(:,:,:) !< fraction land due to glc2000 raw data
    REAL (KIND=wp), INTENT(OUT)  :: ice_glc2000(:,:,:)     !< fraction of ice due to glc2000 raw data
    REAL (KIND=wp), INTENT(OUT)  :: z0_glc2000(:,:,:)      !< roughness length due to glc2000 land use data
    REAL (KIND=wp), INTENT(OUT)  :: root_glc2000(:,:,:)    !< root depth due to glc2000 land use data
    REAL (KIND=wp), INTENT(OUT)  :: plcov_mx_glc2000(:,:,:)!< plant cover maximum due to glc2000 land use data
    REAL (KIND=wp), INTENT(OUT)  :: plcov_mn_glc2000(:,:,:)!< plant cover minimum due to glc2000 land use data
    REAL (KIND=wp), INTENT(OUT)  :: lai_mx_glc2000(:,:,:)  !< Leaf Area Index maximum due to glc2000 land use data
    REAL (KIND=wp), INTENT(OUT)  :: lai_mn_glc2000(:,:,:)  !< Leaf Area Index minimum due to glc2000 land use data
    REAL (KIND=wp), INTENT(OUT)  :: rs_min_glc2000(:,:,:)  !< minimal stomata resistance due to glc2000 land use data
    REAL (KIND=wp), INTENT(OUT)  :: urban_glc2000(:,:,:)   !< urban fraction due to glc2000 land use data
    REAL (KIND=wp), INTENT(OUT)  :: for_d_glc2000(:,:,:)   !< deciduous forest (fraction) due to glc2000 land use data
    REAL (KIND=wp), INTENT(OUT)  :: for_e_glc2000(:,:,:)   !< evergreen forest (fraction) due to glc2000 land use data
    REAL (KIND=wp), INTENT(OUT)  :: emissivity_glc2000(:,:,:) !< longwave emissivity due to glc2000 land use da
     
     INTEGER (KIND=i4) :: undefined_integer ! undef value
     REAL (KIND=wp)    :: default_real


     INTEGER :: l ! counters
     INTEGER :: i_col, j_row ! counter
     INTEGER (KIND=i4) :: i_lu, j_lu
     INTEGER (KIND=i4) :: ie, je, ke  ! indices for target grid elements
     INTEGER (KIND=i4), ALLOCATABLE :: ie_vec(:), je_vec(:), ke_vec(:)  ! indices for target grid elements
     INTEGER (KIND=i4) :: start_cell_id !< ID of starting cell for ICON search
     INTEGER (KIND=i4) :: i1, i2

     REAL (KIND=wp)    :: a_weight(1:tg%ie,1:tg%je,1:tg%ke) !< area weight of all raw data pixels in target grid
     REAL (KIND=wp)    :: a_class(1:tg%ie,1:tg%je,1:tg%ke,1:nclass_glc2000) !< area for each land use class grid  
     
     REAL (KIND=wp)    :: apix      !< area of a raw data pixel
     REAL (KIND=wp)    :: apix_e      !< area of a raw data pixel at equator

     INTEGER :: glc2000_data_row(glc2000_grid%nlon_reg)
     INTEGER :: glc2000_data_pixel(1:1,1:1)
     INTEGER :: lu  ! land use class
     INTEGER :: ncid_glc2000                             !< netcdf unit file number
     CHARACTER (LEN=80) :: varname  !< name of variable
     INTEGER :: varid_glc2000               !< id of variable
     INTEGER :: nlon

     REAL(KIND=wp)   :: point_lon, point_lat

      REAL (KIND=wp) :: pland          !< land cover                      (-)
      REAL (KIND=wp) :: pice           !< ice fraction                    (-)
      REAL (KIND=wp) :: plnz0          !< logarithm of roughness length   (m)
      REAL (KIND=wp) :: proot          !< root depth                      (m)
      REAL (KIND=wp) :: pmn            !< minimal plant cover             (-)
      REAL (KIND=wp) :: pmx            !< maximum plant cover             (-)
      REAL (KIND=wp) :: plaimn         !< minimal leaf area index         (m**2/m**2)
      REAL (KIND=wp) :: plaimx         !< maximum leaf area index         (m**2/m**2)
      REAL (KIND=wp) :: purb           !< urbanisation                    (-)
      REAL (KIND=wp) :: pfor_d         !< deciduous forest                (-)
      REAL (KIND=wp) :: pfor_e         !< evergreen forest                (-)
      REAL (KIND=wp) :: pemissivity    !< surface thermal emissivity      (-)
      REAL (KIND=wp) :: prs_min        !< minimum stomata resistance      (s/m)

     INTEGER        :: k_error     ! error return code

     REAL           :: hp ! height of Prandtl-layer
     REAL (KIND=wp) :: lnhp
     REAL (KIND=wp) :: pwz0 ! weighted summand for z0

     REAL (KIND=wp) :: area_tot   ! total area
     REAL (KIND=wp) :: area_land  ! area with land
     REAL (KIND=wp) :: area_plcov ! area covered with plants

     REAL (KIND=wp) :: bound_north_cosmo !< northern boundary for COSMO target domain
     REAL (KIND=wp) :: bound_south_cosmo !< southern boundary for COSMO target domain
     REAL (KIND=wp) :: bound_west_cosmo  !< western  boundary for COSMO target domain
     REAL (KIND=wp) :: bound_east_cosmo  !< eastern  boundary for COSMO target domain

     ! Some stuff for OpenMP parallelization
     INTEGER :: num_blocks, ib, il, blk_len, istartlon, iendlon, nlon_sub, ishift
!$   INTEGER :: omp_get_max_threads, omp_get_thread_num, thread_id
!$   INTEGER (KIND=i4), ALLOCATABLE :: start_cell_arr(:)

     apix_e  = re * re * deg2rad* ABS(glc2000_grid%dlon_reg) * deg2rad * ABS(glc2000_grid%dlat_reg) 
! area of GLC2000 raw data pixel at equator
     IF (verbose >= idbg_low ) WRITE(logging%fileunit,*)'area pixel at equator: ',apix_e

     hp   = 30.0      ! height of Prandtl-layer
     lnhp = ALOG(hp)

     default_real = 0.0
     undefined_integer= NINT(undefined)

     glc2000_class_fraction = default_real
     glc2000_class_npixel   = undefined_integer
     glc2000_tot_npixel = undefined_integer
     emissivity_glc2000 = 0.0
     fr_land_glc2000 = 0.0
     ice_glc2000 = 0.0
     urban_glc2000 = 0.0
     z0_glc2000 = 0.0
     plcov_mn_glc2000 = 0.0
     plcov_mx_glc2000 = 0.0
     root_glc2000 = 0.0
     lai_mn_glc2000 = 0.0
     lai_mx_glc2000 = 0.0 
     rs_min_glc2000 = 0.0
     for_d_glc2000 = 0.0
     for_e_glc2000 = 0.0
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
     CALL init_glc2000_lookup_tables(nclass_glc2000, &
        &      ilookup_table_glc2000, &
        &      z0_lt_glc2000,            &
        &      lnz0_lt_glc2000,          &
        &      plc_mn_lt_glc2000,        &
        &      plc_mx_lt_glc2000,        &
        &      lai_mn_lt_glc2000,        &
        &      lai_mx_lt_glc2000,        &
        &      rd_lt_glc2000,            &
        &      emiss_lt_glc2000,         &
        &      rs_min_lt_glc2000)

      CALL get_name_glc2000_lookup_tables(ilookup_table_glc2000, name_lookup_table_glc2000)  
     ! open netcdf file 
     IF (verbose >= idbg_low ) WRITE(logging%fileunit,*)TRIM(glc2000_file(1))
     CALL check_netcdf( nf90_open(TRIM(glc2000_file(1)),NF90_NOWRITE, ncid_glc2000))

     varname = 'glc2000byte' 
! I know that the longitude coordinates for the GLC2000 data are stored in a variable called 'glc2000byte'

     CALL check_netcdf( nf90_inq_varid(ncid_glc2000, TRIM(varname), varid_glc2000))
     nlon = glc2000_grid%nlon_reg
     ALLOCATE(ie_vec(nlon),je_vec(nlon),ke_vec(nlon))
     ie_vec(:) = 0
     je_vec(:) = 0
     ke_vec(:) = 0
     start_cell_id = 1

     ! Determine start and end longitude of search
     istartlon = 1
     iendlon = glc2000_grid%nlon_reg
     IF (tg%igrid_type == igrid_icon) THEN
       DO i_col = 1, glc2000_grid%nlon_reg
         point_lon = lon_glc2000(i_col)
         IF (point_lon < tg%minlon) istartlon = i_col + 1
         IF (point_lon > tg%maxlon) THEN
           iendlon = i_col - 1
           EXIT
         ENDIF
       ENDDO
     ELSE IF (tg%igrid_type == igrid_cosmo) THEN
       DO i_col = 1, glc2000_grid%nlon_reg
         point_lon = lon_glc2000(i_col)
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
     IF (verbose >= idbg_low ) WRITE(logging%fileunit,*) 'nlon_sub, num_blocks, blk_len: ',nlon_sub, num_blocks, blk_len

     IF (verbose >= idbg_low ) WRITE(logging%fileunit,*)'Start loop over glc2000 dataset '
     ! loop over rows of GLC2000 dataset
     rows: DO j_row=1,glc2000_grid%nlat_reg
       point_lat = lat_glc2000(j_row)
        
       IF (tg%igrid_type == igrid_cosmo) THEN 
       ! CASE COSMO grid, save some I/O from hard disk if you are out of the target domain
         IF ((point_lat > bound_north_cosmo).OR.(point_lat < bound_south_cosmo) ) THEN ! raw data out of target grid
           CYCLE rows
         ENDIF
       ELSE IF (tg%igrid_type == igrid_icon) THEN 
         IF (point_lat > tg%maxlat .OR. point_lat < tg%minlat) THEN
           CYCLE rows
         ENDIF
       ENDIF ! grid type
        
      ! control output on progress
      IF (verbose >= idbg_high ) THEN
        IF (MOD(j_row,500) == 0  ) THEN
          WRITE(logging%fileunit,*)'nlon: ', nlon
          WRITE(logging%fileunit,*)'j_row: ', j_row
        ENDIF
      ENDIF

       CALL check_netcdf(nf90_get_var(ncid_glc2000, varid_glc2000,  glc2000_data_row,  &
         &               start=(/1,j_row/),count=(/nlon,1/)))
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
       point_lon = lon_glc2000(i_col)

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
         ! raw data pixel within target grid, see output of routine find_rotated_lonlat_grid_element_index
         !- summation of variables
         glc2000_tot_npixel(ie,je,ke) = glc2000_tot_npixel(ie,je,ke) + 1
         a_weight(ie,je,ke) = a_weight(ie,je,ke) + apix  ! sum up area for weight
         lu = glc2000_data_row(i_col)                        ! land use class

          CALL glc2000_look_up(lu, &
            &      nclass_glc2000, &
            &      lnz0_lt_glc2000,          &
            &      plc_mn_lt_glc2000,        &
            &      plc_mx_lt_glc2000,        &
            &      lai_mn_lt_glc2000,        &
            &      lai_mx_lt_glc2000,        &
            &      rd_lt_glc2000,            &
            &      emiss_lt_glc2000,         &
            &      rs_min_lt_glc2000,        &
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

          IF (k_error == 0) THEN ! valid land use class

            glc2000_class_npixel(ie,je,ke,lu) = glc2000_class_npixel(ie,je,ke,lu) + 1
            a_class(ie,je,ke,lu) = a_class(ie,je,ke,lu) + apix   ! sum area of valid land use pixels 
                                                                    !(use as weight later)
            emissivity_glc2000(ie,je,ke) =  emissivity_glc2000(ie,je,ke) + apix * pemissivity
            IF (pland >  0.0) THEN ! only for land pixel

            ! weighted with whole area
            fr_land_glc2000(ie,je,ke) = fr_land_glc2000(ie,je,ke) + apix * pland
            ice_glc2000(ie,je,ke) = ice_glc2000(ie,je,ke) + apix * pice
            urban_glc2000(ie,je,ke) = urban_glc2000(ie,je,ke) + apix * purb
            ! z0 is averaged logarithmic
            IF ( lnhp /= plnz0) THEN ! z0 is averaged logarithmic
              pwz0 = 1./(lnhp - plnz0)
            ELSE
              pwz0 = 0.
            ENDIF
            z0_glc2000(ie,je,ke)      = z0_glc2000(ie,je,ke) + apix * pwz0
            plcov_mn_glc2000(ie,je,ke) = plcov_mn_glc2000(ie,je,ke) + apix * pmn
            plcov_mx_glc2000(ie,je,ke) = plcov_mx_glc2000(ie,je,ke) + apix * pmx

            ! the following fields are weighted with the plant cover
            root_glc2000(ie,je,ke) = root_glc2000(ie,je,ke) + apix * pmx * proot
            lai_mn_glc2000(ie,je,ke) = lai_mn_glc2000(ie,je,ke) + apix * pmx * plaimn
            lai_mx_glc2000(ie,je,ke) = lai_mx_glc2000(ie,je,ke) + apix * pmx * plaimx
            rs_min_glc2000(ie,je,ke) = rs_min_glc2000(ie,je,ke) + apix * pmx * prs_min
            for_d_glc2000(ie,je,ke) = for_d_glc2000(ie,je,ke) + apix * pmx * pfor_d
            for_e_glc2000(ie,je,ke) = for_e_glc2000(ie,je,ke) + apix * pmx * pfor_e
                     
          END IF

        ENDIF
      ENDIF

      ! end loops
    ENDDO columns2
    ENDDO rows
   
     DEALLOCATE(ie_vec,je_vec,ke_vec)
!$   DEALLOCATE(start_cell_arr)

     ! calculate glc2000_class_fraction (glc2000_class_fraction/glc2000_class_npixel)
    DO ke=1, tg%ke
    DO je=1, tg%je
    DO ie=1, tg%ie
      area_tot = a_weight(ie,je,ke)
      area_land = fr_land_glc2000(ie,je,ke)

      ! weight by total area
      IF (area_tot > 0.0) THEN
        emissivity_glc2000(ie,je,ke) = emissivity_glc2000(ie,je,ke) / area_tot
        DO l=1,nclass_glc2000
          glc2000_class_fraction(ie,je,ke,l) =  a_class(ie,je,ke,l) /  area_tot 
          ! area fraction of each land use class
        ENDDO
      ENDIF

      IF (area_land > 0.0 ) THEN
        area_plcov = plcov_mx_glc2000(ie,je,ke) ! area covered with plants
        !area_plcov = MAX(area_plcov, eps)      ! area covered with plants

        z0_glc2000(ie,je,ke) = z0_glc2000(ie,je,ke) / area_land
        z0_glc2000(ie,je,ke) = hp * EXP(-1./z0_glc2000(ie,je,ke))

        ! weight by total area
        fr_land_glc2000(ie,je,ke) = fr_land_glc2000(ie,je,ke) / area_tot
        ice_glc2000(ie,je,ke)   = ice_glc2000(ie,je,ke)   / area_tot


        ! weight by land area
        urban_glc2000(ie,je,ke) = urban_glc2000(ie,je,ke) / area_land
        for_d_glc2000(ie,je,ke) = for_d_glc2000(ie,je,ke) / area_land
        for_e_glc2000(ie,je,ke) = for_e_glc2000(ie,je,ke) / area_land
        plcov_mn_glc2000(ie,je,ke) = plcov_mn_glc2000(ie,je,ke) / area_land

        plcov_mx_glc2000(ie,je,ke) = plcov_mx_glc2000(ie,je,ke) / area_land

        ! weight by area covered with plants
        IF (area_plcov > 0.0) THEN
          root_glc2000(ie,je,ke) = root_glc2000(ie,je,ke)     / area_plcov
          lai_mn_glc2000(ie,je,ke) = lai_mn_glc2000(ie,je,ke) / area_plcov
          lai_mx_glc2000(ie,je,ke) = lai_mx_glc2000(ie,je,ke) / area_plcov
          rs_min_glc2000(ie,je,ke) = rs_min_glc2000(ie,je,ke) / area_plcov
        ELSE ! may occur for ice grid elements
          root_glc2000(ie,je,ke) = undefined
          lai_mn_glc2000(ie,je,ke) = undefined
          lai_mx_glc2000(ie,je,ke) = undefined
          rs_min_glc2000(ie,je,ke) = rs_min_undef
        ENDIF



      ELSE IF ( area_tot > 0.0) THEN ! only sea pixels were found
        fr_land_glc2000(ie,je,ke) = undefined
        z0_glc2000(ie,je,ke)      = undefined
        ice_glc2000(ie,je,ke)     = undefined
        urban_glc2000(ie,je,ke)   = undefined
        for_d_glc2000(ie,je,ke)   = undefined
        for_e_glc2000(ie,je,ke)   = undefined
        plcov_mx_glc2000(ie,je,ke)= undefined
        plcov_mn_glc2000(ie,je,ke)= undefined

        root_glc2000(ie,je,ke)    = undefined
        lai_mx_glc2000(ie,je,ke)  = undefined
        lai_mn_glc2000(ie,je,ke)  = undefined
        rs_min_glc2000(ie,je,ke)  = undefined

      ENDIF
    ENDDO
    ENDDO
    ENDDO

    DO ke=1, tg%ke
    DO je=1, tg%je
    DO ie=1, tg%ie
      IF (glc2000_tot_npixel(ie,je,ke) == 0) THEN ! nearest neighbor search

      point_lon = lon_geo(ie,je,ke)
      point_lat = lat_geo(ie,je,ke)
       CALL find_reg_lonlat_grid_element_index(point_lon,      &
          &                                     point_lat,      &
          &                                     glc2000_grid,  &
          &                                     i_lu,    &
          &                                     j_lu)
        ! get data
        IF ((i_lu /= 0).AND.(j_lu/=0))THEN 

        i_col = i_lu
        j_row = j_lu
        CALL check_netcdf(nf90_get_var(ncid_glc2000, varid_glc2000,  glc2000_data_pixel,  &
          &               start=(/ i_col,j_row /),count=(/ 1,1 /)))

        lu = glc2000_data_pixel(1,1)

        
              CALL glc2000_look_up(lu, &
            &      nclass_glc2000, &
            &      lnz0_lt_glc2000,          &
            &      plc_mn_lt_glc2000,        &
            &      plc_mx_lt_glc2000,        &
            &      lai_mn_lt_glc2000,        &
            &      lai_mx_lt_glc2000,        &
            &      rd_lt_glc2000,            &
            &      emiss_lt_glc2000,         &
            &      rs_min_lt_glc2000,        &
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
            glc2000_class_npixel(ie,je,ke,lu) = glc2000_class_npixel(ie,je,ke,lu) + 1
            a_class(ie,je,ke,lu) = a_class(ie,je,ke,lu) + apix   ! sum area of valid land use pixels 
            emissivity_glc2000(ie,je,ke) =  pemissivity
            IF (pland >  0.0) THEN ! only for land pixel
              fr_land_glc2000(ie,je,ke) = pland
              ice_glc2000(ie,je,ke) =  pice
              urban_glc2000(ie,je,ke) = purb
              IF ( lnhp /= plnz0) THEN ! log z0 
                pwz0 = 1./(lnhp - plnz0)
              ELSE
                pwz0 = 0.
              ENDIF
              z0_glc2000(ie,je,ke)      =  pwz0
              z0_glc2000(ie,je,ke) = hp * EXP(-1./z0_glc2000(ie,je,ke))
              plcov_mn_glc2000(ie,je,ke) = pmn
              plcov_mx_glc2000(ie,je,ke) = pmx
              root_glc2000(ie,je,ke) = proot
              lai_mn_glc2000(ie,je,ke) = plaimn
              lai_mx_glc2000(ie,je,ke) = plaimx
              rs_min_glc2000(ie,je,ke) = prs_min
              for_d_glc2000(ie,je,ke) = pfor_d * pmx
              for_e_glc2000(ie,je,ke) = pfor_e * pmx
              IF (pice == 1.0) THEN  ! ice land use class
                root_glc2000(ie,je,ke) = undefined
                lai_mn_glc2000(ie,je,ke) = undefined
                lai_mx_glc2000(ie,je,ke) = undefined
                rs_min_glc2000(ie,je,ke) = rs_min_undef
              ENDIF
            ENDIF
          ELSE ! not a valid land use class
            fr_land_glc2000(ie,je,ke) = undefined
            z0_glc2000(ie,je,ke)      = undefined
            ice_glc2000(ie,je,ke)     = undefined
            urban_glc2000(ie,je,ke)   = undefined
            for_d_glc2000(ie,je,ke)   = undefined
            for_e_glc2000(ie,je,ke)   = undefined
            plcov_mx_glc2000(ie,je,ke)= undefined
            plcov_mn_glc2000(ie,je,ke)= undefined
            root_glc2000(ie,je,ke)    = undefined
            lai_mx_glc2000(ie,je,ke)  = undefined
            lai_mn_glc2000(ie,je,ke)  = undefined
            rs_min_glc2000(ie,je,ke)  = undefined
          ENDIF
      ENDIF ! nearest neighbour search
   ENDDO
   ENDDO
   ENDDO

  
    ! close netcdf file 
    call check_netcdf( nf90_close(ncid_glc2000))

  END SUBROUTINE agg_glc2000_data_to_target_grid

END MODULE mo_agg_glc2000
