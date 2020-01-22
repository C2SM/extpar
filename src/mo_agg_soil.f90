!+ Fortran module to aggregate the FAO DSMW to the target grid
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_0         2010/12/21 Hermann Asensio
!  Initial release
! V1_1         2011/01/20 Hermann Asensio
!   small bug fixes accroding to Fortran compiler warnings
! V1_7         2013/01/25 Guenther Zaengl 
!   Parallel threads for ICON and COSMO using Open-MP, 
!   Several bug fixes and optimizations for ICON search algorithm, 
!   particularly for the special case of non-contiguous domains; 
!   simplified namelist control for ICON  
! V2_0         2013/06/04 Martina Messmer
!   adaptations such that additionally the HWSD data set 
!   can be aggregated (topsoil and subsoil)
!   Code received from Juergen Helmert.
! V1_14        2014-07-18 Juergen Helmert
!  Combined COSMO Release
! V2_1         2015-01-12 Juergen Helmert
!  Bugfix correction covers CSCS SVN r5907-r6359
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!> Fortran module to aggregate the FAO DSMW to the target grid
!> \author Hermann Asensio
MODULE mo_agg_soil

  USE mo_kind, ONLY: wp, i4, i4

  USE mo_soil_data,       ONLY: default_soiltype
  USE mo_soil_data,       ONLY: FAO_data, HWSD_data, HWSD_map, soil_data
  USE mo_grid_structures, ONLY: reg_lonlat_grid, &
       &                        target_grid_def, &
       &                        igrid_icon,      &
       &                        igrid_cosmo
  USE mo_logging

  USE mo_search_ll_grid, ONLY: find_reg_lonlat_grid_element_index

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: agg_soil_data_to_target_grid
  PUBLIC :: nearest_soil_data_to_target_grid

CONTAINS

  !> Subroutine to aggregate soil data to target grid
  SUBROUTINE agg_soil_data_to_target_grid(tg,          &
       &                   undefined,            &
       &                   soil_texslo,          &
       &                   dsmw_soil_unit,       &
       &                   dsmw_grid,            &
       &                   lon_soil,             &
       &                   lat_soil,             &
       &                   soiltype_fao,         &
       &                   soiltype_hwsd,         &
       &                   fr_land_soil)

    USE mo_soil_data,        ONLY: dsmw_legend

    USE mo_target_grid_data, ONLY: no_raw_data_pixel
    USE mo_target_grid_data, ONLY: lat_geo
    USE mo_target_grid_data, ONLY: search_res !< resolution of ICON grid search index list

    USE mo_search_target_grid, ONLY: find_nearest_target_grid_element

    TYPE(target_grid_def), INTENT(IN) :: tg  !< structure with target grid description

    REAL (wp), INTENT(IN) :: undefined  !< undefined value for soil type
    TYPE(dsmw_legend), INTENT(IN) :: soil_texslo(:)  !< legend for DSMW with texture and slope information, (1:n_unit)
    INTEGER (i4), INTENT(IN) :: dsmw_soil_unit(:,:) 
    !< FAO Digital Soil Map of the World, the values represent the soil unit number (see for legend in variable soil_texslo)
    TYPE(reg_lonlat_grid), INTENT(IN) :: dsmw_grid 
    !< structure with defenition of the raw data grid for the FAO Digital Soil Map of the World

    REAL (wp), INTENT(IN)  :: lon_soil(:)          
    !< longitide coordinates of the soil grid in the geographical (lonlat) system, dimension (nlon_reg)
    REAL (wp), INTENT(IN)  :: lat_soil(:)          
    !< latitude coordinates of the soil grid in the geographical (lonlat) system, dimension (nlat_reg)
    INTEGER(i4), INTENT(OUT) :: soiltype_fao(:,:,:) !< soiltype due to FAO Digital Soil map of the World
    INTEGER(i4), INTENT(OUT) :: soiltype_hwsd(:,:,:) !< store HWSD soiltype IDs
    REAL(wp), INTENT(OUT), OPTIONAL :: fr_land_soil(:,:,:) 
    !< fraction land due to FAO Digital Soil map of the World

    ! local variables

    ! varibles with results of aggregation
    REAL (wp)       :: texture(tg%ie,tg%je,tg%ke)  !< mean texture on COSMO grid

    ! utility variables and intermediate results
    REAL (wp)       :: Z_texture(tg%ie,tg%je,tg%ke)!< sum of texture values
    REAL (wp)       :: Z_slope(tg%ie,tg%je,tg%ke)  !< sum of slope values

    INTEGER (i4) :: I_sea(tg%ie,tg%je,tg%ke)    !< number of sea pixels
    INTEGER (i4) :: I_land(tg%ie,tg%je,tg%ke)    !< number of  land pixels

    INTEGER (i4) :: I_lake(tg%ie,tg%je,tg%ke)   !< number of lake pixels
    INTEGER (i4) :: I_ice(tg%ie,tg%je,tg%ke)    !< number of ice pixels
    INTEGER (i4) :: I_rock(tg%ie,tg%je,tg%ke)   !< number of rock pixels
    INTEGER (i4) :: I_salt(tg%ie,tg%je,tg%ke)   !< number of salt pixels
    INTEGER (i4) :: I_hist(tg%ie,tg%je,tg%ke)   !< number of histosol pixels
    INTEGER (i4) :: I_nodata(tg%ie,tg%je,tg%ke) !< number of nodata pixels
    INTEGER (i4) :: I_slope(tg%ie,tg%je,tg%ke)  !< number of pixels with defined slope
    INTEGER (i4) :: I_texture(tg%ie,tg%je,tg%ke)!< number of pixels with defined texture
    INTEGER (i4) :: I_undef_s(tg%ie,tg%je,tg%ke)!< number of undefined pixels with regard to slope
    INTEGER (i4) :: I_undef_t(tg%ie,tg%je,tg%ke)!< number of undefined pixels with regard to texture

    INTEGER (i4) :: undefined_integer

    INTEGER (i4) :: ir ! counter
    INTEGER (i4) :: jr ! counter
    INTEGER (i4) :: ie  ! counter for grid element index
    INTEGER (i4) :: je  ! counter for grid element index
    INTEGER (i4) :: ke ! counter for grid element index
    INTEGER (i4) :: i1, i2

    REAL (wp) :: lon_pixel ! longitude coordinate of raw data pixel
    REAL (wp) :: lat_pixel ! latitude coordinate of raw data pixel

    INTEGER (i4) :: soil_unit      ! soil unit number
    INTEGER (i4) :: soil_code      ! soil code number

    REAL (wp) :: zcoarse ! help variables
    REAL (wp) :: zmedium
    REAL (wp) :: zfine
    REAL (wp) :: zundef
    REAL (wp) :: zsum_tex
    REAL (wp) :: zflat
    REAL (wp) :: zhilly
    REAL (wp) :: zsteep

    REAL (wp) :: zsum_slope
    REAL (wp) :: zmix

    REAL (wp) :: zsoil ! help variable
    INTEGER (i4) :: itex ! help variable
    INTEGER (i4) :: isoil ! help variable

    INTEGER (i4)  :: dominant_part !< dominant part (undefined, sea or textured soil) for target grid element)

    INTEGER (i4) :: ocean          ! < soil code for ocean
    INTEGER (i4) :: inland_water   ! < soil code for inland water
    INTEGER (i4) :: glacier_ice    ! < soil code for glacier and ice
    INTEGER (i4) :: rock           ! < soil code for rock
    INTEGER (i4) :: salt           ! < soil code for salt
    INTEGER (i4) :: histosols      ! < soil code for histosols
    INTEGER (i4) :: no_data_flag   ! < soil code for no data flag
    INTEGER (i4) :: dunes          ! < soil code for dunes

    REAL (wp) :: bound_north_cosmo !< northern boundary for COSMO target domain
    REAL (wp) :: bound_south_cosmo !< southern boundary for COSMO target domain
    INTEGER (i4) :: start_cell_id !< ID of starting cell for ICON search

    !undefined_integer= NINT(undefined)

    undefined_integer= 0

    no_raw_data_pixel = undefined_integer
    texture = undefined

    !zw = undefined
    Z_texture = undefined
    Z_slope = undefined

    fr_land_soil = undefined

    I_sea = undefined_integer
    I_land = undefined_integer
    I_lake = undefined_integer
    I_ice = undefined_integer
    I_rock = undefined_integer
    I_salt = undefined_integer
    I_hist = undefined_integer
    I_nodata = undefined_integer
    I_slope = undefined_integer
    I_texture = undefined_integer
    I_undef_s = undefined_integer
    I_undef_t = undefined_integer 

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
    IF (verbose >= idbg_low ) WRITE(logging%fileunit,*) 'soil_data: ', soil_data

    lat_loop: DO jr=1,dsmw_grid%nlat_reg
      raw_loop: DO ir=1,dsmw_grid%nlon_reg
        ! find target data grid element index which is nearest to the raw data grid
        lon_pixel = lon_soil(ir)
        lat_pixel = lat_soil(jr)
        SELECT CASE(tg%igrid_type)
        CASE(igrid_cosmo)  ! COSMO GRID
          !>JH probably not working with ICON!!
          IF (lat_pixel>bound_north_cosmo .or. lat_pixel<bound_south_cosmo) THEN
            CYCLE lat_loop
          ENDIF
        END SELECT
        !<JH
        ! Reset start cell when entering a new row or when the previous data point was outside
        ! the model domain
        IF (tg%igrid_type == igrid_icon .AND. (ir == 1 .OR. start_cell_id == 0)) THEN
          i1 = NINT(lon_pixel*search_res)
          i2 = NINT(lat_pixel*search_res)
          start_cell_id = tg%search_index(i1,i2)
          IF (start_cell_id == 0) EXIT raw_loop ! in this case, the whole row is empty
        ENDIF

        CALL  find_nearest_target_grid_element(lon_pixel,    &
             & lat_pixel,     &
             & tg,            &
             & start_cell_id, &
             & ie,            &
             & je,            &
             & ke)

        IF ((ie == 0).OR.(je == 0)) CYCLE raw_loop ! Raw data pixel out of range of target domain

        ! count number of raw data pixel in the target grid element
        no_raw_data_pixel(ie,je,ke) = no_raw_data_pixel(ie,je,ke) + 1
        SELECT CASE(soil_data)
        CASE(FAO_data)
          soil_unit = dsmw_soil_unit(ir,jr)
        CASE(HWSD_data)
          soil_unit = dsmw_soil_unit(ir,jr)
          soiltype_hwsd(ie,je,ke) =    soil_unit
        CASE(HWSD_map)
          soil_unit = MAX(0,dsmw_soil_unit(ir,jr))         
        END SELECT

        zcoarse = 0.0 
        zmedium = 0.0
        zfine   = 0.0
        zundef  = 0.0
        zundef  = 0.0
        zflat   = 0.0
        zhilly  = 0.0
        zsteep  = 0.0

        SELECT CASE (soil_data)
        CASE(HWSD_data)
          ocean = 1 ! is 1 in raw data!
          no_data_flag = undefined_integer
          inland_water = undefined_integer
          glacier_ice = undefined_integer
          rock = undefined_integer
          salt = undefined_integer
          histosols = undefined_integer
          dunes = undefined_integer
        CASE(HWSD_map)
          ocean = 1
          inland_water = 9
          glacier_ice = 1
          rock = 2
          salt = 10
          histosols = 8
          no_data_flag = 255
          dunes = 11
        CASE(FAO_data)
          ocean = 0
          inland_water = 9000
          glacier_ice = 9001
          rock = 9002
          salt = 9003
          histosols = 9004
          no_data_flag = 9009
          dunes = 9005
        END SELECT


        ! decode soil unit to soil texture and slope information with the legend of DSMW (stored in soil_texslo datastructure)
        IF (soil_unit == ocean) then ! ocean
          I_sea(ie,je,ke) = I_sea(ie,je,ke) + 1
        ELSE
          I_land(ie,je,ke) =  I_land(ie,je,ke) + 1
        SELECT CASE (soil_data)
        CASE(FAO_data, HWSD_map)
          soil_code = soil_texslo(soil_unit)%dsmw_code ! the legend has some special cases for the "soil_code"
        CASE(HWSD_data)
          soil_code = MAX(0,dsmw_soil_unit(ir,jr))
END SELECT

          dsmwcode:  IF(soil_code == inland_water) THEN ! inland water
            I_lake(ie,je,ke)    = I_lake(ie,je,ke) + 1
          ELSEIF(soil_code == glacier_ice) THEN ! glacier or ice
            I_ice(ie,je,ke)     = I_ice(ie,je,ke) + 1
            zflat  = 0.0
            zhilly = 1.0
            zsteep = 0.0
            zundef = 0.0
          ELSEIF(soil_code == rock) THEN ! rock
            I_rock(ie,je,ke)    = I_rock(ie,je,ke) + 1
            zflat  = 0.0
            zhilly = 1.0
            zsteep = 0.0
            zundef = 0.0
          ELSEIF(soil_code == salt) THEN ! salt
            I_salt(ie,je,ke)    = I_salt(ie,je,ke) + 1  
            zflat  = 1.0
            zhilly = 0.0
            zsteep = 0.0
            zundef = 0.0
          ELSEIF(soil_code == histosols) THEN ! histosols (peat)
            I_hist(ie,je,ke)    = I_hist(ie,je,ke) + 1
            zundef  = soil_texslo(soil_unit)%part_undefined
            zflat   = soil_texslo(soil_unit)%flat
            zhilly  = soil_texslo(soil_unit)%hilly
            zsteep  = soil_texslo(soil_unit)%steep

          ELSEIF(soil_code == no_data_flag) THEN ! no data flag
            I_nodata(ie,je,ke)  = I_nodata(ie,je,ke) + 1
          ELSEIF(soil_code == dunes) THEN ! dunes or shifting sand
            ! set values for dunes and sand to coarse particles
            zundef = 0. ! 
            zcoarse = 1.0
            zmedium = 0.0
            zfine   = 0.0

            zflat   = soil_texslo(soil_unit)%flat
            zhilly  = soil_texslo(soil_unit)%hilly
            zsteep  = soil_texslo(soil_unit)%steep

            zsum_tex = zcoarse + zmedium + zfine     ! area of soil unit with defined texture
            ! put the textured part of the soil unit in a virtual bucket, 
            !steer and get a single value to describe the texture
            zmix = zfine * 0.0 + zmedium * 0.5 + zcoarse * 1.0 ! calculate a weighted mean
            zmix = zmix / zsum_tex
            Z_texture(ie,je,ke) = Z_texture(ie,je,ke) + zmix
            ! Z_texture(ie,je,ke) = dunes
            I_texture(ie,je,ke) = I_texture(ie,je,ke) + 1

          ELSE
            zcoarse = soil_texslo(soil_unit)%tex_coarse ! get texture from legend
            zmedium = soil_texslo(soil_unit)%tex_medium
            zfine   = soil_texslo(soil_unit)%tex_fine
            zundef  = soil_texslo(soil_unit)%part_undefined
            zflat   = soil_texslo(soil_unit)%flat      ! get slope information from legend
            zhilly  = soil_texslo(soil_unit)%hilly
            zsteep  = soil_texslo(soil_unit)%steep



            zsum_tex = zcoarse + zmedium + zfine     ! area of soil unit with defined texture
            IF (zsum_tex >= zundef) THEN ! well defined texture
              ! put the textured part of the soil unit in a virtual bucket, 
              !steer and get a single value to describe the texture
              zmix = zfine * 0.0 + zmedium * 0.5 + zcoarse * 1.0 ! calculate a weighted mean
              IF (zsum_tex /= 0.0) THEN
                zmix = zmix / zsum_tex
              ELSE
                zmix = 0.0
              ENDIF
            SELECT CASE(soil_data)
            CASE(FAO_data, HWSD_map)
              Z_texture(ie,je,ke) = Z_texture(ie,je,ke) + zmix
            CASE(HWSD_data)
              Z_texture(ie,je,ke) = dsmw_soil_unit(ir,jr)
            END SELECT

              I_texture(ie,je,ke) = I_texture(ie,je,ke) + 1
            ELSE
              I_undef_t(ie,je,ke) = I_undef_t(ie,je,ke) + 1
            ENDIF


          ENDIF dsmwcode



        ENDIF ! ocean

        !----------------------------------------------------------------------------------------------
        SELECT CASE (soil_data)
        CASE(FAO_data)
          zsum_slope = zflat   + zhilly  + zsteep  ! area of soil unit with defined slope information
          IF (zsum_slope > zundef) THEN ! well defined slope
            zmix = zflat * 0.0 + zhilly * 0.5 + zsteep * 1.0 ! generate a weighted mean value
            zmix = zmix / zsum_slope
            Z_slope(ie,je,ke) = Z_slope(ie,je,ke) + zmix
            I_slope(ie,je,ke) = I_slope(ie,je,ke) + 1
          ELSE
            I_undef_s(ie,je,ke) = I_undef_s(ie,je,ke) + 1
          ENDIF
        END SELECT
        !----------------------------------------------------------------------------------------------

      ENDDO raw_loop
    ENDDO lat_loop





    !----------------------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------------------

    ! loop through target grid to determine texture and slope for the target grid element


    dominant_part = 0

  SELECT CASE (soil_data)
  CASE(FAO_data, HWSD_map)
    texture = -99.           ! undefined flag
    fr_land_soil = -99.      ! undefined flag
    soiltype_fao = -99       ! undefined flag
  END SELECT
    DO ke=1, tg%ke
      DO je=1, tg%je
        target_grid: DO ie=1, tg%ie

          IF (no_raw_data_pixel(ie,je,ke) /= 0) THEN ! data for target grid element found

            ! fr_land_soil(ie,je,ke) =  I_land(ie,je,ke) / no_raw_data_pixel(ie,je,ke)  
            ! fr_land_soil has value "1" for land point and "0" for ocean points (lake cosidere as land here)

            IF (PRESENT(fr_land_soil)) THEN
              fr_land_soil(ie,je,ke) =  (I_land(ie,je,ke)-I_lake(ie,je,ke) ) / &
                   &   no_raw_data_pixel(ie,je,ke)  ! fr_land as water-land mask
            ENDIF

          SELECT CASE (soil_data)
          CASE(FAO_data, HWSD_map)
            texture(ie,je,ke) = -9. ! default for ocean, texture(ie,je,ke) is overwritten for other soil types

            ! set I_nodata as dominant part at start
            dominant_part = I_nodata(ie,je,ke) !  CASE(no_data_flag)
            IF (I_nodata(ie,je,ke) > I_sea(ie,je,ke) ) THEN ! set nodata flag 
              texture(ie,je,ke) = -9009.
            ENDIF


            ! texture information for target grid element
            !----------------------------------------------------------------------------------------------
            ! set I_nodata as dominant part at start
            dominant_part = I_nodata(ie,je,ke) !  CASE(no_data_flag)

            IF (I_undef_t(ie,je,ke) > dominant_part) then ! undefined texture part
              dominant_part = I_undef_t(ie,je,ke)
              texture(ie,je,ke) = -9012.
            ENDIF

            IF (I_lake(ie,je,ke) > dominant_part) then ! water pixel 
              texture(ie,je,ke) = -9000. 
            ENDIF

            !IF (I_lake(ie,je,ke) > dominant_part) then ! water pixel 
!dominant_part = I_lake(ie,je,ke)           ! avoid using lake pixels for partial "land grid elements"
            ! even if only very few raw data pixel define a valid soil type
            !texture(ie,je,ke) = -5.
            !ENDIF

            IF (I_ice(ie,je,ke) > dominant_part) then ! glacier or ice
              dominant_part = I_ice(ie,je,ke)
              texture(ie,je,ke) = -9001.
            ENDIF

            IF (I_rock(ie,je,ke) > dominant_part) then ! rock
              dominant_part = I_rock(ie,je,ke)
              texture(ie,je,ke) = -9002.
            ENDIF

            IF (I_salt(ie,je,ke) > dominant_part) then ! salt
              dominant_part = I_salt(ie,je,ke)
              texture(ie,je,ke) = -9003.
            ENDIF

            IF (I_hist(ie,je,ke) > dominant_part) then ! histosol (peat etc)
              dominant_part = I_hist(ie,je,ke)
              texture(ie,je,ke) = -9004.
            ENDIF

            !    IF (I_dunes(ie,je,ke) > dominant_part) then ! dunes/ shifting sands
            !      dominant_part = I_dunes(ie,je,ke)
            !      texture(ie,je,ke) = -9005. ! set to sand 
            !    ENDIF


            IF (I_texture(ie,je,ke) > dominant_part) then ! textured soil is dominant part in grid element
              dominant_part = I_texture(ie,je,ke)
              texture(ie,je,ke) = Z_texture(ie,je,ke) / real(I_texture(ie,je,ke))
            ENDIF
            !----------------------------------------------------------------------------------------------

            !----------------------------------------------------------------------------------------------
            ! convert from texture information to soil type 
            ! (1 ice, 2 rock, 3 sand, 4 sandy loam, 5 loam, 6 loamy clay, 7 clay, 8 histosol(peat),9 sea point)
            zsoil = texture(ie,je,ke)
            isoil = -99
            itex = NINT(100 * texture(ie,je,ke)) ! texture in percent as Integer

            IF (itex == -900100)   isoil =  1 ! ice, glacier (soil type 1) 
            IF (itex == -900200)   isoil =  2 ! rock, lithosols (soil type 2)
            IF (itex == -900300)   isoil =  3 ! salt, set soiltype to sand (soil type 3)

            IF (itex == -900400)   isoil =  8 ! histosol, e.g. peat (soil type 8)
            IF (itex == -900)      isoil = 9 ! undefined (ocean)

            IF (itex == -900500)   isoil = 3 ! shifting sands or dunes, set soiltype to sand (soil type 3)
            IF (itex == -900000)   isoil = 9 ! undefined (inland lake)
            IF (itex == -900900)   isoil = 5 !  default_soiltype ! undefined (nodata), set soiltype to loam (soil type )
            IF (itex == -901200)   isoil = 5 !  default_soiltype 
            ! undefined (dominant part undefined), set soiltype to loam (soil type 5)

            IF (itex >= 0)    isoil = 7 ! fine textured, clay (soil type 7)
            IF (itex >= 20)    isoil = 6 ! medium to fine textured, loamy clay (soil type 6)
            IF (itex >= 40)    isoil = 5 ! medium textured, loam (soil type 5)
            IF (itex >= 60)    isoil = 4 ! coarse to medium textured, sandy loam (soil type 4)
            IF (itex >= 80)    isoil = 3 ! coarse textured, sand (soil type 3)


            soiltype_fao(ie,je,ke) = isoil
          END SELECT
            if (soiltype_fao(ie,je,ke) < 1) then
              print*,'Aggregation Problem!!! - Soiltype < 1!',isoil, zsoil,itex,default_soiltype
            end if


            !----------------------------------------------------------------------------------------------

          ENDIF ! data for target grid element found ! data for target grid element found

        ENDDO target_grid
      ENDDO
    ENDDO


  END SUBROUTINE agg_soil_data_to_target_grid

  !> Subroutine to aggregate soil data to COSMO grid
  SUBROUTINE nearest_soil_data_to_target_grid(tg,              &
       &                             undefined,          &
       &                             soil_texslo,        &
       &                             dsmw_soil_unit,     &
       &                             dsmw_grid,          &
       &                             soiltype_fao,       &
       &                             soiltype_hwsd,       &
       &                             fr_land_soil)

    USE mo_soil_data,       ONLY: dsmw_legend

    USE mo_target_grid_data, ONLY: no_raw_data_pixel
    USE mo_target_grid_data, ONLY: lon_geo
    USE mo_target_grid_data, ONLY: lat_geo

    TYPE(target_grid_def), INTENT(IN) :: tg  !< structure with target grid description
    REAL (wp), INTENT(IN) :: undefined  !< undefined value for soil type

    TYPE(dsmw_legend), INTENT(IN) :: soil_texslo(:)  !< legend for DSMW with texture and slope information, (1:n_unit)
    INTEGER (i4), INTENT(IN) :: dsmw_soil_unit(:,:) 
    !< FAO Digital Soil Map of the World, the values represent the soil unit number (see for legend in variable soil_texslo)
    TYPE(reg_lonlat_grid), INTENT(IN) :: dsmw_grid 
    !< structure with defenition of the raw data grid for the FAO Digital Soil Map of the World

    !< latitude coordinates of the soil grid in the geographical (lonlat) system, dimension (nlat_reg)
    INTEGER(i4), INTENT(INOUT) :: soiltype_fao(:,:,:) !< soiltype due to FAO Digital Soil map of the World
    INTEGER(i4), INTENT(INOUT) :: soiltype_hwsd(:,:,:) !< store HWSD soil IDs
    REAL(wp), INTENT(INOUT), OPTIONAL:: fr_land_soil(:,:,:) 
    !< fraction land due to FAO Digital Soil map of the World

    ! local variables

    ! varibles with results of aggregation
    REAL (wp)       :: texture(tg%ie,tg%je,tg%ke)!< texture values

    INTEGER (i4) :: undefined_integer

    INTEGER (i4) :: ie  ! counter for grid element index
    INTEGER (i4) :: je  ! counter for grid element index
    INTEGER (i4) :: ke ! counter for grid element index

    REAL (wp) :: lon_target ! longitude coordinate of target grid element
    REAL (wp) :: lat_target ! latitude coordinate of target grid element


    INTEGER (i4) :: soil_ir ! index of raw data pixel (lon axis)
    INTEGER (i4) :: soil_jr ! index of raw data pixel (lat axis)

    INTEGER (i4) :: soil_unit      ! soil unit number
    INTEGER (i4) :: soil_code      ! soil code number

    REAL (wp) :: zcoarse ! help variables
    REAL (wp) :: zmedium
    REAL (wp) :: zfine
    REAL (wp) :: zundef
    REAL (wp) :: zsum_tex
    REAL (wp) :: zflat
    REAL (wp) :: zhilly
    REAL (wp) :: zsteep

    REAL (wp) :: zsum_slope
    REAL (wp) :: zmix

    REAL (wp) :: zsoil ! help variable
    INTEGER (i4) :: itex ! help variable
    INTEGER (i4) :: isoil ! help variable

    INTEGER (i4) :: ocean          ! < soil code for ocean
    INTEGER (i4) :: inland_water   ! < soil code for inland water
    INTEGER (i4) :: glacier_ice    ! < soil code for glacier and ice
    INTEGER (i4) :: rock           ! < soil code for rock
    INTEGER (i4) :: salt           ! < soil code for salt
    INTEGER (i4) :: histosols      ! < soil code for histosols
    INTEGER (i4) :: no_data_flag   ! < soil code for no data flag
    INTEGER (i4) :: dunes          ! < soil code for dunes

    undefined_integer= NINT(undefined)
    SELECT CASE(soil_data)
    CASE(FAO_data,HWSD_map)
      texture = -99.           ! undefined flag
    CASE(HWSD_data)
      texture = 0 ! undefined flag
    END SELECT

    zsum_tex = 0.0
    zmix = 0.0

    zcoarse = 0.0 
    zmedium = 0.0
    zfine   = 0.0

    zundef  = 0.0
    zsum_slope  = 0.0

    zflat   = 0.0
    zhilly  = 0.0
    zsteep  = 0.0

    SELECT CASE (soil_data)
    CASE(HWSD_data, HWSD_map)
      ocean = undefined_integer ! is 1 in raw data!
      no_data_flag = undefined_integer
      inland_water = undefined_integer
      glacier_ice =  undefined_integer
      rock = undefined_integer
      salt = undefined_integer
      histosols = undefined_integer
      dunes =  undefined_integer
!!$   ocean = 9
!!$    inland_water = 9
!!$    glacier_ice = 1
!!$    rock = 2
!!$    salt = 10
!!$    histosols = 8
!!$    no_data_flag = 255
!!$    dunes = 11
    CASE(FAO_data)
      ocean = -9
      inland_water = -9000
      glacier_ice = -9001
      rock = -9002
      salt = -9003
      histosols = -9004
      no_data_flag = -9009
      dunes = -9005
    END SELECT

    ! loop over target grid

    !HA debug
    !PRINT *,'start loop '
    ! loop through target grid,
    ! find nearest neighbour in raw data grid for no_raw_data_pixel == 0
    ! for latitudes < -60 degrees North put all non-ocean grid elements to ice (no data for antarctica) -> not here? 
    !this is for consistency check when land use information are available for glcc data (glc2000 do not include antarctica)
    ! convert from texture information to soiltyp (1 ice, 2 rock, 3 sand, 4 sandy loam, 5 loam, 6 loamy clay, 7 clay, 
    !8 histosol(peat),9 sea point)
    DO ke=1, tg%ke
      DO je=1, tg%je
        target_loop: DO ie=1, tg%ie 
          ! test if there is information of raw data grid in the target grid element, 
          ! if not, find nearest neighbour in target grid (soil data)

          IF (no_raw_data_pixel(ie,je,ke) == 0) THEN ! no data yet for target grid element

            lon_target = lon_geo(ie,je,ke)
            lat_target = lat_geo(ie,je,ke)
            !HA debug:
            !PRINT *,'CALL find_reg_lonlat_grid_element_index'
            !PRINT *,'lon_target: ', lon_target
            !PRINT *,'lat_target: ', lat_target
            !PRINT *,'dsmw_grid: ', dsmw_grid
            CALL find_reg_lonlat_grid_element_index(lon_target,&
                 lat_target,  &
                 dsmw_grid,   &
                 soil_ir,     &
                 soil_jr)
            !PRINT *,'find_reg_lonlat_grid_element_index done'
            IF ((soil_ir == 0).or.(soil_jr == 0) ) then ! problem, target grid element outside raw data grid
              texture(ie,je,ke) = -9. ! ocean
            ELSE

              SELECT CASE(soil_data)
              CASE(HWSD_data)
                soil_unit = dsmw_soil_unit(soil_ir,soil_jr)
                soiltype_hwsd(ie,je,ke) = dsmw_soil_unit(soil_ir,soil_jr)
              CASE(HWSD_map)
                soil_unit = MAX(0,dsmw_soil_unit(soil_ir,soil_jr))
              CASE(FAO_data)
                soil_unit = dsmw_soil_unit(soil_ir,soil_jr)
              END SELECT

              zcoarse = 0.0 
              zmedium = 0.0
              zfine   = 0.0
              zundef  = 0.0
              zflat   = 0.0
              zhilly  = 0.0
              zsteep  = 0.0
              ! decode soil unit to soil texture and slope information with the legend of DSMW 
              !(stored in soil_texslo datastructure)


              IF (soil_unit == 0 ) then ! ocean 
                texture(ie,je,ke) = ocean
                IF (PRESENT(fr_land_soil)) THEN
                  fr_land_soil(ie,je,ke) = 0. ! ocean point
                ENDIF
              ELSE

                IF (PRESENT(fr_land_soil)) THEN
                  fr_land_soil(ie,je,ke) = 1.0 ! land point
                ENDIF
               SELECT CASE(soil_data)
               CASE(FAO_data)
                 soil_code = soil_texslo(soil_unit)%dsmw_code ! the legend has some special cases for the "soil_code"
               CASE(HWSD_map, HWSD_data)
                 soil_code =  MAX(0,dsmw_soil_unit(soil_ir,soil_jr))
               END SELECT

                dsmwcode: IF(soil_code == ABS(INT(inland_water))) THEN ! inland water
                  texture(ie,je,ke) = inland_water
                  IF (PRESENT(fr_land_soil)) THEN
                    fr_land_soil(ie,je,ke) = 0. ! water point
                  ENDIF

                ELSEIF(soil_code == ABS(INT(glacier_ice))) THEN ! glacier or ice
                  texture(ie,je,ke) = glacier_ice
                  zflat  = 0.0
                  zhilly = 1.0
                  zsteep = 0.0
                  zundef = 0.0
                ELSEIF(soil_code == ABS(INT(rock))) THEN ! rock
                  texture(ie,je,ke) = rock
                  zflat  = 0.0
                  zhilly = 1.0
                  zsteep = 0.0
                  zundef = 0.0
                ELSEIF(soil_code == ABS(INT(salt))) THEN ! salt
                  texture(ie,je,ke) = salt
                  zflat  = 1.0
                  zhilly = 0.0
                  zsteep = 0.0
                  zundef = 0.0
                ELSEIF(soil_code == ABS(INT(histosols))) THEN ! histosols (peat)
                  texture(ie,je,ke) = histosols
                  zundef  = soil_texslo(soil_unit)%part_undefined
                  zflat   = soil_texslo(soil_unit)%flat
                  zhilly  = soil_texslo(soil_unit)%hilly
                  zsteep  = soil_texslo(soil_unit)%steep
                ELSEIF(soil_code == ABS(INT(dunes))) THEN ! dunes or shifting sand
                  ! set values for dunes and sand to coarse particles
                  zundef = 0. ! 
                  zcoarse = 1.0
                  zmedium = 0.0
                  zfine   = 0.0

                  zflat   = soil_texslo(soil_unit)%flat
                  zhilly  = soil_texslo(soil_unit)%hilly
                  zsteep  = soil_texslo(soil_unit)%steep

                  zsum_tex = zcoarse + zmedium + zfine     ! area of soil unit with defined texture
                  !IF (zsum_tex > zundef) THEN ! well defined texture
                  ! put the textured part of the soil unit in a virtual bucket, 
                  !steer and get a single value to describe the texture
                  zmix = zfine * 0.0 + zmedium * 0.5 + zcoarse * 1.0 ! calculate a weighted mean
                  zmix = zmix / zsum_tex
                  texture(ie,je,ke) = zmix ! texture(ie,je,ke) = zmix/ 1 ! only one raw data pixel

                  ! texture(ie,je,ke) = -9005.
                  ! zundef  = soil_texslo(soil_unit)%part_undefined
                  ! zflat   = soil_texslo(soil_unit)%flat
                  ! zhilly  = soil_texslo(soil_unit)%hilly
                  ! zsteep  = soil_texslo(soil_unit)%steep
                ELSEIF(soil_code == 9009) THEN ! no data flag FAO
                  texture(ie,je,ke) = no_data_flag
                ELSEIF(soil_code == -9999) THEN ! no data flag HWSD
                  texture(ie,je,ke) = no_data_flag
                  !---------------------------------------------------------------------------------------------- 
                ELSE
                  zcoarse = soil_texslo(soil_unit)%tex_coarse ! get texture from legend
                  zmedium = soil_texslo(soil_unit)%tex_medium
                  zfine   = soil_texslo(soil_unit)%tex_fine
                  zundef  = soil_texslo(soil_unit)%part_undefined
                  zflat   = soil_texslo(soil_unit)%flat      ! get slope information from legend
                  zhilly  = soil_texslo(soil_unit)%hilly
                  zsteep  = soil_texslo(soil_unit)%steep

                  zsum_tex = zcoarse + zmedium + zfine     ! area of soil unit with defined texture
                  IF (zsum_tex >= zundef) THEN ! well defined texture
                    ! put the textured part of the soil unit in a virtual bucket,
                    ! steer and get a single value to describe the texture
                    zmix = zfine * 0.0 + zmedium * 0.5 + zcoarse * 1.0 ! calculate a weighted mean
                    IF (zsum_tex /= 0.0) THEN
                      zmix = zmix / zsum_tex
                    ELSE
                      zmix = 0.0
                    ENDIF
                    texture(ie,je,ke) = zmix ! texture(ie,je,ke) = zmix/ 1 ! only one raw data pixel
                  ELSE
                    texture(ie,je,ke) = -9012.
!!$
!!$                         SELECT CASE(soil_data)
!!$                         CASE(FAO_data)
!!$                           texture(ie,je,ke) = zmix ! texture(ie,je,ke) = zmix/ 1 ! only one raw data pixel
!!$                         CASE(HWSD_data)
!!$                           texture(ie,je,ke) = dsmw_soil_unit(soil_ir,soil_jr)
!!$                         END SELECT
!!$                      ELSE
!!$                         SELECT CASE(soil_data)
!!$                         CASE(FAO_data)
!!$                           texture(ie,je,ke) = -9012.
!!$                         CASE(HWSD_data)
!!$                           texture(ie,je,ke) = 136.
!!$                         END SELECT

                  ENDIF
                  !----------------------------------------------------------------------------------------------
                ENDIF dsmwcode


                !----------------------------------------------------------------------------------------------
                zsum_slope = zflat   + zhilly  + zsteep  ! area of soil unit with defined slope information
                IF (zsum_slope > zundef) THEN ! well defined slope
                  zmix = zflat * 0.0 + zhilly * 0.5 + zsteep * 1.0 ! generate a weighted mean value
                  IF(zsum_slope /= 0.0) THEN
                    zmix = zmix / zsum_slope
                  ELSE
                    zmix = 0.0
                  ENDIF
                ENDIF

              ENDIF ! ocean

              ! texture(ie,je,ke) = -9.
              ! !HA debug


            ENDIF ! target grid check

            !----------------------------------------------------------------------------------------------
            ! convert from texture information to soil type 
            ! (1 ice, 2 rock, 3 sand, 4 sandy loam, 5 loam, 6 loamy clay, 7 clay, 8 histosol(peat),9 sea point)

            zsoil = texture(ie,je,ke)
            isoil = -99
            itex = NINT(100 * texture(ie,je,ke)) ! texture in percent as Integer

          SELECT CASE (soil_data)
          CASE(FAO_data, HWSD_map)
            IF (itex == -900100)   isoil =  1 ! ice, glacier (soil type 1) 
            IF (itex == -900200)   isoil =  2 ! rock, lithosols (soil type 2)
            IF (itex == -900300)   isoil =  3 ! salt, set soiltype to sand (soil type 3)
            IF (itex == -900400)   isoil =  8 ! histosol, e.g. peat (soil type 8)
            IF (itex == -900)      isoil = 9 ! undefined (ocean)

            IF (itex == -900500)   isoil = 3 ! shifting sands or dunes, set soiltype to sand (soil type 3)
            IF (itex == -900000)   isoil = 9 ! undefined (inland lake)


            IF (itex == -900500)   isoil = 3 ! shifting sands or dunes, set soiltype to sand (soil type 3)
            IF (itex == -900000)   isoil = 9 ! undefined (inland lake)

            IF (itex == -900900)   isoil = 5 !  default_soiltype ! undefined (nodata), set soiltype to loam (soil type )
            IF (itex == -901200)   isoil = 5 !  default_soiltype ! undefined (dominant part undefined), 
            ! set soiltype to loam (soil type 5)




            IF (itex >= 0)    isoil = 7 ! fine textured, clay (soil type 7)
            IF (itex >= 20)    isoil = 6 ! medium to fine textured, loamy clay (soil type 6)
            IF (itex >= 40)    isoil = 5 ! medium textured, loam (soil type 5)
            IF (itex >= 60)    isoil = 4 ! coarse to medium textured, sandy loam (soil type 4)
            IF (itex >= 80)    isoil = 3 ! coarse textured, sand (soil type 3)


            soiltype_fao(ie,je,ke) = isoil
          END SELECT
            if (soiltype_fao(ie,je,ke) < 1) then
              WRITE(logging%fileunit,*)'WARNING: ***Nearest neighbor check: Soiltype < 1!  ',isoil, zsoil,itex,default_soiltype, '***'
            end if

            !----------------------------------------------------------------------------------------------

          ENDIF ! no data yet for target grid element

        ENDDO target_loop
      ENDDO
    ENDDO


  END SUBROUTINE nearest_soil_data_to_target_grid

END MODULE mo_agg_soil
