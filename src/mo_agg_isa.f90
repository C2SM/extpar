!+ Fortran module to aggregate isa land use data to a target grid
!
!
! Description:
! Fortran module to aggregate isa land use data to a target grid
!
! Current Code Owner: DWD, <Name of person responsible for this code>
!    <smail, phone, fax and email>
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V2_4         2015-05-21 Hendrik Wouters 
!  Initial release based on mo_agg_globcover.f90 V1_14
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!=======================================================================
!> Fortran module to aggregate isa land use data to a target grid
!!
!> \author Hendrik Wouters
MODULE mo_agg_isa

  USE mo_logging
  USE mo_kind,                  ONLY: wp, i4
  USE mo_io_units,              ONLY: filename_max
  USE mo_io_utilities,          ONLY: check_netcdf

  USE mo_grid_structures,       ONLY: igrid_icon, &
       &                              reg_lonlat_grid,&
       &                              target_grid_def, &
       &                              igrid_cosmo

  USE mo_search_ll_grid,        ONLY: find_reg_lonlat_grid_element_index

  USE mo_search_target_grid,    ONLY: find_nearest_target_grid_element


  USE netcdf,                   ONLY: &
       &                              nf90_open,              &
       &                              nf90_close,             &
       &                              nf90_inq_varid,         &
       &                              nf90_get_var,           &
       &                              nf90_nowrite

  USE mo_isa_data,              ONLY: isa_grid, &
       &                              lon_isa,  &
       &                              lat_isa,  &
       &                              ntiles_isa, &
       &                              max_tiles_isa

  USE mo_isa_routines,          ONLY: det_band_isa_data, &
         &                            get_isa_data_block
        
  USE mo_math_constants,        ONLY: deg2rad
  USE mo_physical_constants,    ONLY: re

  USE mo_target_grid_data,      ONLY: lon_geo, & !< longitude coordinates of the COSMO grid in the geographical system
    &                                 lat_geo, & !< latitude coordinates of the COSMO grid in the geographical system
    &                                 search_res

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: agg_isa_data_to_target_grid

  REAL(KIND=wp), PARAMETER :: rs_min_undef=999. !< undefined value for minimal stomata resistance

  CONTAINS

  !> Subroutine to aggregate isa data to the target grid
  SUBROUTINE agg_isa_data_to_target_grid(isa_file,          &
       &                                 undefined,               &
       &                                 isa_tiles_grid,    &
       &                                 tg,                      &
       &                                 isa_tot_npixel, &
       &                                 isa_field)

    CHARACTER (LEN=filename_max), INTENT(IN) :: isa_file(1:max_tiles_isa)  !< filename isa raw data
    TYPE(reg_lonlat_grid), INTENT(IN)        :: isa_tiles_grid(:)  ! grid structure of isa tiles
    TYPE(target_grid_def), INTENT(IN)        :: tg  !< structure with target grid description
    REAL (KIND=wp), INTENT(IN)               :: undefined            !< undef value

    INTEGER (KIND=i4), INTENT(OUT)           :: isa_tot_npixel(:,:,:)  
    REAL (KIND=wp), INTENT(OUT)              :: isa_field(:,:,:)   !< urban fraction due to isa land use data

    !local variables
    TYPE(reg_lonlat_grid)                    :: ta_grid

    INTEGER(KIND=i4)                         :: nt, &           ! counter
         &                                      i_col, j_row, & ! counter
         &                                      isa_data_row(isa_grid%nlon_reg), &
         &                                      isa_data_pixel(1:1,1:1), &
         &                                      isa, &   
         &                                      ncid_isa(1:ntiles_isa), &            !< netcdf unit file number
         &                                      varid_isa, &               !< id of variable
         &                                      nlon, &
         &                                      block_row_start, &
         &                                      block_row, &
         &                                      mlat, &
         &                                      tile, &
         &                                      errorcode, &   ! error return code
         &                                      i_isa, j_isa, &
         &                                      undefined_integer, & ! undef value
         &                                      ie, je, ke, &  ! indices for target grid elements
         &                                      start_cell_id, & !< ID of starting cell for ICON search
         &                                      i1, i2
                                             
    REAL(KIND=wp)                            :: a_weight(1:tg%ie,1:tg%je,1:tg%ke), &
          &                                     default_real, &
          &                                     apix, &      !< area of a raw data pixel
          &                                     apix_e, &      !< area of a raw data pixel at equator
          &                                     point_lon, point_lat, &
          &                                     area_tot, &   ! total area
          &                                     bound_north_cosmo, & !< northern boundary for COSMO target domain
          &                                     bound_south_cosmo, & !< southern boundary for COSMO target domain
          &                                     bound_west_cosmo, &  !< western  boundary for COSMO target domain
          &                                     bound_east_cosmo  !< eastern  boundary for COSMO target domain

    INTEGER (KIND=i4), ALLOCATABLE           :: ie_vec(:), je_vec(:), ke_vec(:)  ! indices for target grid elements
    REAL (KIND=wp), ALLOCATABLE              :: isa_block(:,:)  ! a block of ISA land use data
                                             
    CHARACTER (LEN=80)                       :: varname  !< name of variable

    ! Some stuff for OpenMP parallelization
    INTEGER(KIND=i4)                         :: num_blocks, ib, il, blk_len,&
          &                                     istartlon, iendlon, nlon_sub, ishift

    CALL logging%info('Enter routine: agg_isa_data_to_target_grid') 

    apix_e  = re * re * deg2rad* ABS(isa_grid%dlon_reg) * deg2rad * ABS(isa_grid%dlat_reg) 

    default_real = 0.0
    undefined_integer= NINT(undefined)

    isa_tot_npixel = undefined_integer

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

    nlon = isa_grid%nlon_reg
    ALLOCATE(ie_vec(nlon),je_vec(nlon),ke_vec(nlon))
    ie_vec(:) = 0
    je_vec(:) = 0
    ke_vec(:) = 0
    start_cell_id = 1
    
    ! open netcdf file
    DO nt = 1,ntiles_isa
      CALL check_netcdf( nf90_open(TRIM(isa_file(nt)),NF90_NOWRITE, ncid_isa(nt)))
    END DO

    varname = 'ISA' ! I know that the isa data are stored in a variable called 'ISA'

    CALL check_netcdf(nf90_inq_varid(ncid_isa(1), TRIM(varname), varid_isa))

    mlat = 1
    block_row_start = mlat

    CALL det_band_isa_data(isa_grid,block_row_start,ta_grid)

    IF (ALLOCATED(isa_block)) THEN
      DEALLOCATE(isa_block, STAT = errorcode)
      IF(errorcode /= 0) CALL logging%error('Cant deallocate the isa_block',__FILE__,__LINE__)
    END IF
    ALLOCATE (isa_block(1:ta_grid%nlon_reg, 1:ta_grid%nlat_reg), STAT = errorcode)
    IF (errorcode /= 0) CALL logging%error('Cant allocate the isa_block',__FILE__,__LINE__)

    CALL get_isa_data_block(ta_grid,               &
                                  isa_tiles_grid,  &
                                  ncid_isa,       &
                                  isa_block)

    block_row = 0

    ! Determine start and end longitude of search
    istartlon = 1
    iendlon = isa_grid%nlon_reg
    IF (tg%igrid_type == igrid_icon) THEN
      DO i_col = 1, isa_grid%nlon_reg
        point_lon = lon_isa(i_col)
        IF (point_lon < tg%minlon) istartlon = i_col + 1
        IF (point_lon > tg%maxlon) THEN
          iendlon = i_col - 1
          EXIT
        ENDIF
      ENDDO
    ELSE IF (tg%igrid_type == igrid_cosmo) THEN
      DO i_col = 1, isa_grid%nlon_reg
        point_lon = lon_isa(i_col)
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

    WRITE(message_text,*) 'nlon_sub: ',nlon_sub,' num_blocks: ',num_blocks, ' blk_len: ',blk_len
    CALL logging%info(message_text)

    CALL logging%info('Start loop over isa rows...')

    isa_rows: DO mlat = 1, isa_grid%nlat_reg
      block_row= block_row + 1
      IF(block_row > ta_grid%nlat_reg) THEN ! read in new block
        block_row_start = mlat
        block_row = 1
        CALL det_band_isa_data(isa_grid,block_row_start,ta_grid)
        IF(ALLOCATED(isa_block)) THEN
          DEALLOCATE(isa_block, STAT=errorcode)
          IF(errorcode/=0) CALL logging%error('Cant deallocate the isa_block')
        ENDIF
        ALLOCATE (isa_block(1:ta_grid%nlon_reg,1:ta_grid%nlat_reg), STAT=errorcode)
        IF(errorcode/=0) CALL logging%error('Cant allocate isa_block',__FILE__,__LINE__)
        CALL get_isa_data_block(ta_grid,              &
             &                  isa_tiles_grid, &
             &                  ncid_isa,       &
             &                  isa_block)
      ENDIF

      point_lat = lat_isa(mlat)
      IF (tg%igrid_type == igrid_cosmo) THEN ! CASE COSMO grid, save some I/O from hard disk if you are out or the target doma &
        IF ((point_lat > bound_north_cosmo).OR.(point_lat < bound_south_cosmo) ) THEN ! raw data out of target grid
          CYCLE isa_rows
        ENDIF
      ELSE IF (tg%igrid_type == igrid_icon) THEN
        IF (point_lat > tg%maxlat .OR. point_lat < tg%minlat) THEN
          CYCLE isa_rows
        ENDIF
      ENDIF ! grid type

      isa_data_row(1:nlon) = INT(isa_block(1:nlon,block_row)) !type conversion
      apix = apix_e * COS(point_lat * deg2rad) ! area of raw data pixel (in [m**2])

      ie_vec(istartlon:iendlon) = 0
      IF (tg%igrid_type /= igrid_icon) THEN
        je_vec(:) = 0
        ke_vec(:) = 0
      ENDIF

!OMP PARALLEL DO PRIVATE(ib,il,i_col,i1,i2,ishift,point_lon,thread_id,start_cell_id)
      DO ib = 1, num_blocks

  !     thread_id = omp_get_thread_num()+1
  !     start_cell_id = start_cell_arr(thread_id)
        ishift = istartlon-1+(ib-1)*blk_len

        columns1: DO il = 1,blk_len
          i_col = ishift+il
          IF (i_col > iendlon) CYCLE columns1
          ! find the corresponding target grid indices
          point_lon = lon_isa(i_col)
          
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
!     start_cell_arr(thread_id) = start_cell_id
      ENDDO !num_blocks
!OMP END PARALLEL DO

      columns2: DO i_col=istartlon,iendlon
        ! find the corresponding target grid indices
        ie = ie_vec(i_col)
        je = je_vec(i_col)
        ke = ke_vec(i_col)

        IF ((ie /= 0).AND.(je/=0).AND.(ke/=0))THEN
          ! raw data pixel within target grid, see output of routine find_rotated_lonlat_grid_element_index
          
          ! input data should be smaller than 100, otherwise it represents an invalid value (255 in our case)
          IF (isa_data_row(i_col) <= 100) THEN
            isa = isa_data_row(i_col)
            isa_tot_npixel(ie,je,ke) = isa_tot_npixel(ie,je,ke) + 1
            a_weight(ie,je,ke) = a_weight(ie,je,ke) + apix  ! sum up area for weight
              
            ! weighted with whole area
            isa_field(ie,je,ke) = isa_field(ie,je,ke) + 0.01_wp * apix * isa
          END IF
        ENDIF
        ! end loops
      ENDDO columns2
    ENDDO isa_rows

    CALL logging%info('...done')

    DEALLOCATE(ie_vec,je_vec,ke_vec)

    DO ke=1, tg%ke
      DO je=1, tg%je
        DO ie=1, tg%ie
          area_tot = a_weight(ie,je,ke)
          IF (area_tot > 0.0 ) THEN
            ! weight by land area
            isa_field(ie,je,ke) = isa_field(ie,je,ke) / area_tot
          ELSE
            isa_field(ie,je,ke)   = undefined
          ENDIF
        ENDDO
      ENDDO
    ENDDO

    DO ke=1, tg%ke
      DO je=1, tg%je
        DO ie=1, tg%ie
          IF (isa_tot_npixel(ie,je,ke) == 0) THEN ! nearest neighbor search
            point_lon = lon_geo(ie,je,ke)
            point_lat = lat_geo(ie,je,ke)

            CALL find_reg_lonlat_grid_element_index(point_lon,      &
              &                                     point_lat,      &
              &                                     isa_grid, &
              &                                     i_isa,           &
              &                                     j_isa,           &
              &                                     ntiles_isa,     &
              &                                     tile = tile,    &
              &                                     regular_tiles_grid_info=isa_tiles_grid)

            IF ((i_isa /= 0).AND.(j_isa /= 0))THEN
           
              i_col = i_isa
              j_row = j_isa
             
              IF (ntiles_isa > 1) THEN
                CALL check_netcdf(nf90_open(TRIM(isa_file(tile)),NF90_NOWRITE, ncid_isa(tile)))
                CALL check_netcdf(nf90_inq_varid(ncid_isa(tile), TRIM(varname), varid_isa))        
                CALL check_netcdf(nf90_get_var(ncid_isa(tile), varid_isa,  isa_data_pixel,  &
                  &               start=(/ i_col,j_row /),count=(/ 1,1 /)))
                CALL check_netcdf(nf90_close(ncid_isa(tile)))
              ELSE
                CALL check_netcdf(nf90_get_var(ncid_isa(tile), varid_isa,  isa_data_pixel,  &
                    &               start=(/ i_col,j_row /),count=(/ 1,1 /)))
              END IF
             
              ! values above 100 are invalid
              IF (isa_data_pixel(1,1) <= 100) THEN
                isa = isa_data_pixel(1,1)
                apix = apix_e * COS(point_lat * deg2rad) ! area of raw data pixel (in [m**2])
                isa_field(ie,je,ke) = 0.01_wp * isa
              ELSE
                isa = 0
                isa_field(ie,je,ke) = undefined
              ENDIF
             
            ELSE
              isa = 0
              isa_field(ie,je,ke) = undefined
            ENDIF
          ENDIF ! nearest neighbour search
        ENDDO
      ENDDO
    ENDDO

    CALL logging%info('Exit routine: agg_isa_data_to_target_grid')

  END SUBROUTINE agg_isa_data_to_target_grid

END MODULE mo_agg_isa
