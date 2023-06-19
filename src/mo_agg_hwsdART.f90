MODULE mo_agg_hwsdART
  USE mo_logging
  USE mo_kind, ONLY: wp,i4


  USE mo_hwsdART_data,    ONLY: type_clay_heavy, &
  & type_silty_clay, &
  & type_clay_light, &
  & type_silty_clay_loam, &
  & type_clay_loam, &
  & type_silt, &
  & type_silt_loam, &
  & type_sandy_clay,&
  & type_loam, &
  & type_sandy_clay_loam, &
  & type_sandy_loam, &
  & type_loamy_sand, &
  & type_sand, &
  & undef_hwsdARTtype, &
  & no_data
  
  USE mo_hwsdART_tg_fields, ONLY: fr_heavy_clay, &
  & fr_silty_clay, &
  & fr_light_clay, &
  & fr_silty_clay_loam, &
  & fr_clay_loam, &
  & fr_silt, &
  & fr_silt_loam, &
  & fr_sandy_clay, &
  & fr_loam, &
  & fr_sandy_clay_loam, &
  & fr_sandy_loam, &
  & fr_loamy_sand, &
  & fr_sand, &
  & fr_undef

  USE mo_grid_structures, ONLY: reg_lonlat_grid, &
   & target_grid_def, &
   & igrid_icon, &
   & igrid_cosmo 

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: agg_hwsdART_data_to_target_grid

    CONTAINS

       !> Subroutine to aggregate hwsdART data to target grid
       SUBROUTINE agg_hwsdART_data_to_target_grid(tg,          &
                  &                   hwsdART_soil_unit,    &
                  &                   hwsdART_grid,         &
                  &                   lon_hwsdART,          &
                  &                   lat_hwsdART)

       USE mo_target_grid_data, ONLY: no_raw_data_pixel, &
        & lon_geo, &
        & lat_geo, &
        & search_res !< resolution of ICON grid search index list

       USE mo_search_target_grid, ONLY: find_nearest_target_grid_element

       TYPE(target_grid_def), INTENT(IN) :: tg  !< structure with target grid description

       INTEGER (KIND=i4), INTENT(IN) :: hwsdART_soil_unit(:,:) 
       TYPE(reg_lonlat_grid), INTENT(IN) :: hwsdART_grid ! grid of hwsdART data
       
       REAL (KIND=wp), INTENT(IN)  :: lon_hwsdART(:), &          
       & lat_hwsdART(:)          

       ! local variables
       INTEGER (KIND=i4) :: I_clay_heavy(tg%ie,tg%je,tg%ke)      , & !< number of heavy clay pixels
        &  I_silty_clay(tg%ie,tg%je,tg%ke)      , & !< number of silty clay pixels
        &  I_clay_light(tg%ie,tg%je,tg%ke)      , & !< number of light clay pixels
        &  I_silty_clay_loam(tg%ie,tg%je,tg%ke) , & !< number of silty clay loam pixels
        &  I_clay_loam(tg%ie,tg%je,tg%ke)       , & !< number of clay loam pixels
        &  I_silt(tg%ie,tg%je,tg%ke)            , & !< number of silt pixels
        &  I_silt_loam(tg%ie,tg%je,tg%ke)       , & !< number of silt loam pixels
        &  I_sandy_clay(tg%ie,tg%je,tg%ke)      , & !< number of sandy clay pixels
        &  I_loam(tg%ie,tg%je,tg%ke)            , & !< number of loam pixels
        &  I_sandy_clay_loam(tg%ie,tg%je,tg%ke) , & !< number of sandy clay loam pixels
        &  I_sandy_loam(tg%ie,tg%je,tg%ke)      , & !< number of sandy loam pixels
        &  I_loamy_sand(tg%ie,tg%je,tg%ke)      , & !< number of loamy sand pixels
        &  I_sand(tg%ie,tg%je,tg%ke)            , & !< number of sand pixels
        &  I_undef(tg%ie,tg%je,tg%ke)           , &
        &   soil_unit, &
        &   undefined_integer, &
        &   start_cell_id, & !< ID of starting cell for ICON search
        &   ie, &  ! counter for grid element index
        &   je, &  ! counter for grid element index
        &   ke, & ! counter for grid element index
        &   i1, &
        &   i2

       REAL (KIND=wp) :: lon_pixel, & ! longitude coordinate of raw data pixel
       & lat_pixel, & ! latitude coordinate of raw data pixel
       & bound_north_cosmo, & !< northern boundary for COSMO target domain
       & bound_south_cosmo    !< southern boundary for COSMO target domain

       INTEGER :: ir, & ! counter
       & jr ! counter

       undefined_integer= 0

       no_raw_data_pixel = undefined_integer

       I_clay_heavy      = undefined_integer
       I_silty_clay      = undefined_integer
       I_clay_light      = undefined_integer
       I_silty_clay_loam = undefined_integer
       I_clay_loam       = undefined_integer
       I_silt            = undefined_integer
       I_silt_loam       = undefined_integer
       I_loam            = undefined_integer
       I_sandy_clay      = undefined_integer
       I_sandy_clay_loam = undefined_integer
       I_loamy_sand      = undefined_integer
       I_sandy_loam      = undefined_integer
       I_sand            = undefined_integer
       I_undef           = undefined_integer

       start_cell_id = 1


      SELECT CASE(tg%igrid_type)
        CASE(igrid_icon)  ! ICON GRID
          ke = 1
        CASE(igrid_cosmo)  ! COSMO GRID
          ke = 1
          bound_north_cosmo = MAXVAL(lat_geo) + 0.05_wp  ! add some "buffer"
          bound_north_cosmo = MIN(bound_north_cosmo,90.0_wp)
          bound_south_cosmo = MINVAL(lat_geo) - 0.05_wp  ! add some "buffer"
          bound_south_cosmo = MAX(bound_south_cosmo,-90.0_wp)
      END SELECT

       ! loop over raw data grid
      lat_loop: DO jr=1,hwsdART_grid%nlat_reg
      WRITE(message_text,*) 'jr = ',jr
      CALL logging%info(message_text)  
      raw_loop: DO ir=1,hwsdART_grid%nlon_reg
      ! find target data grid element index which is nearest to the raw data grid
      lon_pixel = lon_hwsdART(ir)
      lat_pixel = lat_hwsdART(jr)
      ! Reset start cell when entering a new row or when the previous data point was outside
      ! the model domain
      IF (tg%igrid_type == igrid_icon .AND. (ir == 1 .OR. start_cell_id == 0)) THEN
        i1 = NINT(lon_pixel*search_res)
        i2 = NINT(lat_pixel*search_res)
        start_cell_id = tg%search_index(i1,i2)
        !drieg: debug
!        IF (start_cell_id == 0) EXIT raw_loop ! in this case, the whole row is empty
      ENDIF
      
      CALL  find_nearest_target_grid_element(lon_pixel,     &
                                           & lat_pixel,     &
                                           & tg,            &
                                           & start_cell_id, &
                                           & ie,            &
                                           & je,            &
                                           & ke)

      IF ((ie == 0).OR.(je == 0)) CYCLE raw_loop ! Raw data pixel out of range of target domain

      ! count number of raw data pixel in the target grid element
      no_raw_data_pixel(ie,je,ke) = no_raw_data_pixel(ie,je,ke) + 1

      soil_unit = hwsdART_soil_unit(ir,jr)

      IF (soil_unit == type_clay_heavy) then      ! heavy clay
        I_clay_heavy(ie,je,ke) = I_clay_heavy(ie,je,ke) + 1
      ENDIF
      IF (soil_unit == type_silty_clay) then      ! silty_clay
        I_silty_clay(ie,je,ke) = I_silty_clay(ie,je,ke) + 1
      ENDIF
      IF (soil_unit == type_clay_light) then      ! clay light
        I_clay_light(ie,je,ke) = I_clay_light(ie,je,ke) + 1
      ENDIF
      IF (soil_unit == type_silty_clay_loam) then ! silty clay loam
        I_silty_clay_loam(ie,je,ke) = I_silty_clay_loam(ie,je,ke) + 1
      ENDIF
      IF (soil_unit == type_clay_loam) then       ! clay loam
        I_clay_loam(ie,je,ke) = I_clay_loam(ie,je,ke) + 1
      ENDIF
      IF (soil_unit == type_silt) then            ! silt
        I_silt(ie,je,ke) = I_silt(ie,je,ke) + 1
      ENDIF
      IF (soil_unit == type_silt_loam) then       ! silt loam
        I_silt_loam(ie,je,ke) = I_silt_loam(ie,je,ke) + 1
      ENDIF
      IF (soil_unit == type_sandy_clay) then      ! sandy clay
        I_sandy_clay(ie,je,ke) = I_sandy_clay(ie,je,ke) + 1
      ENDIF
      IF (soil_unit == type_loam) then            ! loam
        I_loam(ie,je,ke) = I_loam(ie,je,ke) + 1
      ENDIF
      IF (soil_unit == type_sandy_clay_loam) then ! sandy clay loam
        I_sandy_clay_loam(ie,je,ke) = I_sandy_clay_loam(ie,je,ke) + 1
      ENDIF
      IF (soil_unit == type_sandy_loam) then      ! sandy loam
        I_sandy_loam(ie,je,ke) = I_sandy_loam(ie,je,ke) + 1
      ENDIF
      IF (soil_unit == type_loamy_sand) then      ! loamy sand
        I_loamy_sand(ie,je,ke) = I_loamy_sand(ie,je,ke) + 1
      ENDIF
      IF (soil_unit == type_sand) then            ! sand
        I_sand(ie,je,ke) = I_sand(ie,je,ke) + 1
      ENDIF
      IF (soil_unit == undef_hwsdARTtype) then    ! ELSE
        I_undef(ie,je,ke) = I_undef(ie,je,ke) + 1
      ENDIF
      IF (soil_unit == no_data) then              ! ELSE
        I_undef(ie,je,ke) = I_undef(ie,je,ke) + 1
      ENDIF

       ENDDO raw_loop
     ENDDO lat_loop

     DO ke=1, tg%ke
       DO je=1, tg%je
       WRITE(message_text,*) 'je = ',je
       CALL logging%info(message_text)  
       target_grid: DO ie=1, tg%ie

         IF (no_raw_data_pixel(ie,je,ke) /= 0) THEN ! data for target grid element found

           fr_heavy_clay(ie,je,ke)     = REAL(I_clay_heavy(ie,je,ke))     &
                                       & / REAL(no_raw_data_pixel(ie,je,ke))
           fr_silty_clay(ie,je,ke)     = REAL(I_silty_clay(ie,je,ke))     &
                                       & / REAL(no_raw_data_pixel(ie,je,ke))
           fr_light_clay(ie,je,ke)     = REAL(I_clay_light(ie,je,ke))     &
                                       & / REAL(no_raw_data_pixel(ie,je,ke))
           fr_silty_clay_loam(ie,je,ke)= REAL(I_silty_clay_loam(ie,je,ke))&
                                       & / REAL(no_raw_data_pixel(ie,je,ke))
           fr_clay_loam(ie,je,ke)      = REAL(I_clay_loam(ie,je,ke))      &
                                       & / REAL(no_raw_data_pixel(ie,je,ke))
           fr_silt(ie,je,ke)           = REAL(I_silt(ie,je,ke))           &
                                       & / REAL(no_raw_data_pixel(ie,je,ke))
           fr_silt_loam(ie,je,ke)      = REAL(I_silt_loam(ie,je,ke))      &
                                       & / REAL(no_raw_data_pixel(ie,je,ke))
           fr_sandy_clay(ie,je,ke)     = REAL(I_sandy_clay(ie,je,ke))     &
                                       & / REAL(no_raw_data_pixel(ie,je,ke))
           fr_loam(ie,je,ke)           = REAL(I_loam(ie,je,ke))           &
                                       & / REAL(no_raw_data_pixel(ie,je,ke))
           fr_sandy_clay_loam(ie,je,ke)= REAL(I_sandy_clay_loam(ie,je,ke))&
                                       & / REAL(no_raw_data_pixel(ie,je,ke))
           fr_sandy_loam(ie,je,ke)     = REAL(I_sandy_loam(ie,je,ke))     &
                                       & / REAL(no_raw_data_pixel(ie,je,ke))
           fr_loamy_sand(ie,je,ke)     = REAL(I_loamy_sand(ie,je,ke))     &
                                       & / REAL(no_raw_data_pixel(ie,je,ke))
           fr_sand(ie,je,ke)           = REAL(I_sand(ie,je,ke))           &
                                       & / REAL(no_raw_data_pixel(ie,je,ke))
           fr_undef(ie,je,ke)          = REAL(I_undef(ie,je,ke))          &
                                       & / REAL(no_raw_data_pixel(ie,je,ke))
         ELSE
           fr_undef(ie,je,ke)          = 1.0
         ENDIF
       ENDDO target_grid
     ENDDO
   ENDDO
   END SUBROUTINE agg_hwsdART_data_to_target_grid
END MODULE mo_agg_hwsdART
