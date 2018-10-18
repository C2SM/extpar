!+ Fortran module to aggregate lake depth data to a target grid 
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_0         2010/12/21 Hermann Asensio
!  Initial release
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
!> Fortran module to aggregate lake depth data to a target grid 
!!
!> \author Hermann Asensio
MODULE mo_agg_flake

  !> kind parameters are defined in MODULE data_parameters
  USE mo_kind, ONLY: wp
  USE mo_kind, ONLY: i8
  USE mo_kind, ONLY: i4

  !> abort_extpar defined in MODULE utilities_extpar
  USE mo_utilities_extpar, ONLY: abort_extpar


  !> data type structures form module GRID_structures
  USE mo_grid_structures, ONLY: reg_lonlat_grid, &
    &                           rotated_lonlat_grid
  
  USE mo_grid_structures, ONLY: igrid_icon
  USE mo_grid_structures, ONLY: igrid_cosmo

  USE mo_search_ll_grid, ONLY: find_reg_lonlat_grid_element_index, &
    &                          find_rotated_lonlat_grid_element_index
  USE mo_io_units,          ONLY: filename_max
  USE mo_io_utilities, ONLY: check_netcdf

  USE mo_search_target_grid, ONLY: find_nearest_target_grid_element


  USE netcdf,      ONLY :   &
    & nf90_open,              &
    & nf90_close,             &
    & nf90_inq_varid,         &
    & nf90_get_var,           &
    & nf90_get_att,           &
    & NF90_NOWRITE,           &
    & nf90_noerr,             &
    & nf90_strerror

  

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: agg_flake_data_to_target_grid


  CONTAINS

  !> Subroutine to aggregate flake data to the target grid
  SUBROUTINE agg_flake_data_to_target_grid(flake_file, &
    &                                      undefined,  &
    &                                      tg,         &
    &                                      lake_depth, &
    &                                      fr_lake,    &
    &                                      flake_tot_npixel)
    
    
  !-------------------------------------------------------------------------------------
  ! list of modules which are used as "input"
  USE mo_grid_structures, ONLY: target_grid_def   !< type definition of structure for tg

  ! raw data grid description, here for GLCC data
  USE mo_flake_data, ONLY: flake_grid, &
    &                          lon_flake,  &
    &                          lat_flake


    ! USE structure which contains the definition of the ICON grid
    USE  mo_icon_grid_data, ONLY: ICON_grid !< structure which contains the definition of the ICON grid

    ! USE structure which contains the definition of the COSMO grid
    USE  mo_cosmo_grid, ONLY: COSMO_grid !< structure which contains the definition of the COSMO grid

    USE mo_math_constants, ONLY: pi, rad2deg, deg2rad, eps
    USE mo_physical_constants, ONLY: re
       
    USE mo_flake_data, ONLY: flake_depth_undef !< default value for undefined lake depth
    USE mo_flake_data, ONLY: flake_depth_default  !< default value for default lake depth, 10 [m]


    ! USE global data fields (coordinates)
    USE mo_target_grid_data, ONLY: lon_geo, & !< longitude coordinates of the COSMO grid in the geographical system 
      &                            lat_geo !< latitude coordinates of the COSMO grid in the geographical system
    USE mo_target_grid_data, ONLY: search_res !< resolution of ICON grid search index list


     CHARACTER (LEN=filename_max), INTENT(IN) :: flake_file  !< filename flake raw data
     REAL (KIND=wp), INTENT(IN) :: undefined            !< undef value

     TYPE(target_grid_def), INTENT(IN) :: tg  !< structure with target grid description

     REAL (KIND=wp), INTENT(OUT)  :: lake_depth(:,:,:) !< lake depth
     REAL (KIND=wp), INTENT(OUT)  :: fr_lake(:,:,:)     !< fraction of fresh water (lakes)

     INTEGER (KIND=i8), INTENT(OUT) :: flake_tot_npixel(:,:,:)  !< total number of flake raw data pixels on target grid 
     INTEGER (KIND=i8) :: undefined_integer ! undef value
     REAL (KIND=wp)    :: default_real
     INTEGER :: i,j,k,l ! counters
     INTEGER :: i_col, j_row ! counter
     INTEGER (KIND=i8) :: ie, je, ke  ! indices for target grid elements
     INTEGER (KIND=i8), ALLOCATABLE :: ie_vec(:), je_vec(:), ke_vec(:)  ! indices for target grid elements
     INTEGER (KIND=i8) :: i1, i2
     INTEGER (KIND=i8) :: n_flake_pixel(1:tg%ie,1:tg%je,1:tg%ke)  !< number of raw data pixel with lakes
     REAL (KIND=wp)    :: a_weight(1:tg%ie,1:tg%je,1:tg%ke) !< area weight of all raw data pixels in target grid
     REAL (KIND=wp)    :: latw      !< latitude weight (for area weighted mean)
     REAL (KIND=wp)    :: apix      !< area of a raw data pixel
     REAL (KIND=wp)    :: apix_e      !< area of a raw data pixel at equator
     INTEGER :: flake_data_row(flake_grid%nlon_reg)
     INTEGER :: flake_data_pixel
     INTEGER :: flake_data_pixels(1:1,1:1)

     REAL :: scale_factor
     INTEGER :: ncid_flake                             !< netcdf unit file number
     CHARACTER (LEN=80) :: varname  !< name of variable
     INTEGER :: varid_flake               !< id of variable
     INTEGER :: nlon

     REAL(KIND=wp)   :: point_lon, point_lat

     INTEGER        :: k_error     ! error return code

     REAL (KIND=wp) :: area_tot   ! total area
     REAL (KIND=wp) :: area_land  ! area with land

     REAL (KIND=wp) :: bound_north_cosmo !< northern boundary for COSMO target domain
     REAL (KIND=wp) :: bound_south_cosmo !< southern boundary for COSMO target domain
     REAL (KIND=wp) :: bound_west_cosmo  !< western  boundary for COSMO target domain
     REAL (KIND=wp) :: bound_east_cosmo  !< eastern  boundary for COSMO target domain

     
     REAL (KIND=wp) :: lon_target ! longitude coordinate of target grid element
     REAL (KIND=wp) :: lat_target ! latitude coordinate of target grid element

     
     INTEGER (KIND=i8) :: flake_ir ! index of raw data pixel (lon axis)
     INTEGER (KIND=i8) :: flake_jr ! index of raw data pixel (lat axis)

     INTEGER (KIND=i8) :: start_cell_id ! start cell ID for ICON search

     ! Some stuff for OpenMP parallelization
     INTEGER :: num_blocks, ib, il, blk_len, istartlon, iendlon, nlon_sub, ishift
!$   INTEGER :: omp_get_max_threads, omp_get_thread_num, thread_id
!$   INTEGER (KIND=i8), ALLOCATABLE :: start_cell_arr(:)


     apix_e  = re * re * deg2rad* ABS(flake_grid%dlon_reg) * deg2rad * ABS(flake_grid%dlat_reg) 
! area of GLCC raw data pixel at equator
     PRINT *,'area pixel at equator: ',apix_e

     default_real = 0.0
     undefined_integer= NINT(undefined)

     lake_depth = default_real ! set to 0.
     fr_lake = default_real ! set to 0.
     flake_tot_npixel   = undefined_integer

     n_flake_pixel = 0

     a_weight = default_real
     
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

     ! open netcdf file 
     PRINT *,'open ',TRIM(flake_file)
     CALL check_netcdf( nf90_open(TRIM(flake_file),NF90_NOWRITE, ncid_flake))

     varname = 'DEPTH_LK' ! I know that the lake depth data are stored in a variable called 'DEPTH_LK'

     CALL check_netcdf( nf90_inq_varid(ncid_flake, TRIM(varname), varid_flake)) ! get varid for variable DEPTH_LK
     CALL check_netcdf( nf90_get_att(ncid_flake,varid_flake,"scale_factor",scale_factor)) ! get scale factor for lake depth
     nlon = flake_grid%nlon_reg
     ALLOCATE(ie_vec(nlon),je_vec(nlon),ke_vec(nlon))
     ie_vec(:) = 0
     je_vec(:) = 0
     ke_vec(:) = 0
     start_cell_id = 1

     ! Determine start and end longitude of search
     istartlon = 1
     iendlon = flake_grid%nlon_reg
     IF (tg%igrid_type == igrid_icon) THEN
       DO i_col = 1, flake_grid%nlon_reg
         point_lon = lon_flake(i_col)
         IF (point_lon < tg%minlon) istartlon = i_col + 1
         IF (point_lon > tg%maxlon) THEN
           iendlon = i_col - 1
           EXIT
         ENDIF
       ENDDO
     ELSE IF (tg%igrid_type == igrid_cosmo) THEN
       DO i_col = 1, flake_grid%nlon_reg
         point_lon = lon_flake(i_col)
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

     PRINT *,'Start loop over flake dataset ',flake_grid%nlat_reg
     ! loop over rows of GLCC dataset
     rows: DO j_row=1,flake_grid%nlat_reg

      if (MOD(j_row, 100) == 0) PRINT '(a,i6,f7.2)', ' FLAKE processing row: ', j_row, lat_flake(j_row)

       point_lat = lat_flake(j_row)
        
       IF (tg%igrid_type == igrid_cosmo) THEN ! CASE COSMO grid, save some I/O from hard disk if you are out or the target domain
         IF ((point_lat > bound_north_cosmo).OR.(point_lat < bound_south_cosmo) ) THEN ! raw data out of target grid
           CYCLE rows
         ENDIF
       ELSE IF (tg%igrid_type == igrid_icon) THEN 
         IF (point_lat > tg%maxlat .OR. point_lat < tg%minlat) THEN
           CYCLE rows
         ENDIF
       ENDIF ! grid type
        
       ! read in pixels
       CALL check_netcdf(nf90_get_var(ncid_flake, varid_flake,  flake_data_row,  &
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
       point_lon = lon_flake(i_col)

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
         flake_tot_npixel(ie,je,ke) = flake_tot_npixel(ie,je,ke) + 1
         a_weight(ie,je,ke) = a_weight(ie,je,ke) + apix  ! sum up area for weight

         IF (flake_data_row(i_col) > 0) THEN   ! lake pixel
           lake_depth(ie,je,ke) =  lake_depth(ie,je,ke) + flake_data_row(i_col) * scale_factor
           n_flake_pixel(ie,je,ke) = n_flake_pixel(ie,je,ke) + 1
           fr_lake(ie,je,ke)    = fr_lake(ie,je,ke) + apix
         ENDIF

       ENDIF

      ! end loops
    ENDDO columns2
    ENDDO rows

     DEALLOCATE(ie_vec,je_vec,ke_vec)
!$   DEALLOCATE(start_cell_arr)
 
    ! calculate flake_class_fraction (flake_class_fraction/flake_class_npixel)
    DO ke=1, tg%ke
    DO je=1, tg%je
    DO ie=1, tg%ie
      IF (a_weight(ie,je,ke) > 0.0) THEN
        fr_lake(ie,je,ke) = fr_lake(ie,je,ke) / a_weight(ie,je,ke)
        IF (n_flake_pixel(ie,je,ke)> 0) THEN ! flake data found
           lake_depth(ie,je,ke) =  lake_depth(ie,je,ke)/n_flake_pixel(ie,je,ke)
        ELSE
          lake_depth(ie,je,ke) = flake_depth_undef ! set depth to -1
        ENDIF
        
      ENDIF
    ENDDO
    ENDDO
    ENDDO

    DO ke=1, tg%ke
    DO je=1, tg%je
    DO ie=1, tg%ie
      IF (flake_tot_npixel(ie,je,ke) == 0) THEN ! nearest neighbour search
        lon_target = MIN(lon_geo(ie,je,ke),flake_grid%end_lon_reg)
        lat_target = lat_geo(ie,je,ke)
        ! nearest neighbour search
        CALL find_reg_lonlat_grid_element_index(lon_target,      &
          &                                     lat_target,      &
          &                                     flake_grid,  &
          &                                     flake_ir,    &
          &                                     flake_jr)
        ! get data
        i_col = flake_ir
        j_row = flake_jr
!        PRINT *,lon_target,lat_target,i_col,j_row
        CALL check_netcdf(nf90_get_var(ncid_flake, varid_flake,  flake_data_pixels,  &
          &               start=(/ i_col,j_row /),count=(/ 1,1 /)))

        IF (flake_data_pixels(1,1) > 0) THEN   ! lake pixel
          lake_depth(ie,je,ke) =   flake_data_pixels(1,1) * scale_factor
          fr_lake(ie,je,ke)    = 1.0

        ELSE
          lake_depth(ie,je,ke) = flake_depth_undef
          fr_lake(ie,je,ke)    = 0.0
        ENDIF
      ENDIF
    ENDDO
    ENDDO
    ENDDO


   ! close netcdf file 
    CALL check_netcdf( nf90_close(ncid_flake))


  END SUBROUTINE agg_flake_data_to_target_grid


END MODULE mo_agg_flake
