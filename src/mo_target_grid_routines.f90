!+ Fortran module with routines to initialize and calculate target grid fields and coordinates
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
! V1_7         2013/01/25 Guenther Zaengl 
!    Parallel threads for ICON and COSMO using Open-MP, 
!    Several bug fixes and optimizations for ICON search algorithm, 
!    particularly for the special case of non-contiguous domains; 
!    simplified namelist control for ICON  
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!> Fortran module with routines to initialize and calculate target grid fields and coordinates
!> \author Hermann Asensio
!>
MODULE mo_target_grid_routines

  USE mo_kind,                  ONLY: wp, i4, i4
  USE mo_io_units,              ONLY: filename_max
  USE mo_utilities_extpar,      ONLY: abort_extpar
  USE mo_logging

  USE mo_read_extpar_namelists, ONLY: read_namelists_extpar_grid_def

  USE mo_math_constants,        ONLY: rad2deg, deg2rad
  USE mo_grid_structures,       ONLY: igrid_icon, &
    &                                 igrid_cosmo

  USE mo_icon_grid_data,        ONLY: icon_grid, & !< structure which contains the definition of the ICON grid
    &                                 icon_grid_region, &
    &                                 icon_dom_def

  USE mo_icon_grid_routines,    ONLY: get_icon_grid_info, &
    &                                 get_icon_domain_info, &
    &                                 init_icon_grid

  USE mo_cosmo_grid,            ONLY: cosmo_grid, &
    &                                 get_cosmo_grid_info, &
    &                                 lon_rot, lat_rot, &
    &                                 allocate_cosmo_rc, &
    &                                 calculate_cosmo_target_coordinates

  USE mo_target_grid_data,      ONLY: tg, lon_geo, lat_geo, &
    &                                 allocate_com_target_fields, &
    &                                 search_res

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: init_target_grid

  CONTAINS

  !> allocate target grid and determine coordinates for each grid element
  SUBROUTINE init_target_grid(namelist_grid_def,lrad)

    ! arguments
    CHARACTER(len=*), INTENT(IN)  :: namelist_grid_def
    LOGICAL, INTENT(IN), OPTIONAL :: lrad

    ! local variables
    CHARACTER (len=filename_max)  :: icon_coor_file, & !< filname of the ICON grid file with the coordinates
      &                             domain_def_namelist !< namelist file with domain definition

    INTEGER (KIND=i4)             :: igrid_type, & !< target grid type, 1 for ICON, 2 for COSMO
      &                              i1,i2,ii,i1s,i1e,i2s,i2e
    INTEGER                       :: i,j,k !< counters
    LOGICAL                       :: lzrad, lonadj

    REAL(wp)                      :: d1,d2
    REAL(wp), ALLOCATABLE         :: auxlon(:,:)


    CALL logging%info('Enter routine: init_target_grid')
    !----------------------------------------------------------------------------------------------
    lzrad = .FALSE.
    IF (PRESENT(lrad)) lzrad = lrad

    ! get information on target grid
    !-----------------------------------------------------------------------------------------------
     
    CALL read_namelists_extpar_grid_def(namelist_grid_def, &
                                           igrid_type, &
                                           domain_def_namelist)

    ! Checks
    IF ((igrid_type /= igrid_cosmo) .AND. lzrad) THEN
      CALL logging%error('lradtopo only implemented for the COSMO grid', __FILE__, __LINE__)
    ENDIF
    !--------------------------------------------------------------------------------------
    ! read in target grid information

    SELECT CASE(igrid_type)
      !-----------------------------------------------------------------
      CASE(igrid_icon) ! ICON GRID

        CALL get_icon_grid_info(domain_def_namelist,tg,icon_grid,icon_coor_file)
        CALL get_icon_domain_info(icon_coor_file,icon_dom_def)
        CALL init_icon_grid(icon_dom_def)

      CASE(igrid_cosmo) ! COSMO grid
        CALL get_cosmo_grid_info(domain_def_namelist,tg,COSMO_grid,lzrad)
        
        ! allocate structure cosmo_rot_coor
        CALL allocate_cosmo_rc(tg%ie,tg%je)
    END SELECT

    !-----------------------------------------------------------------
    !-----------------------------------------------------------------
    ! allocate target_fields
    CALL allocate_com_target_fields(tg)
    !-----------------------------------------------------------------

    SELECT CASE(igrid_type)
         !-----------------------------------------------------------------
      CASE(igrid_icon) ! ICON GRID
        k=1
        j=1
        DO i=1,tg%ie
          lon_geo(i,j,k) = rad2deg * icon_grid_region%cells%center(i)%lon ! convert von radians to degrees
          lat_geo(i,j,k) = rad2deg * icon_grid_region%cells%center(i)%lat ! convert von radians to degrees
        ENDDO

        tg%minlon = 999._wp
        tg%minlat = 999._wp
        tg%maxlon = -999._wp
        tg%maxlat = -999._wp

        DO i=1,icon_grid_region%nverts
          tg%minlon = MIN(tg%minlon,icon_grid_region%verts%vertex(i)%lon)
          tg%maxlon = MAX(tg%maxlon,icon_grid_region%verts%vertex(i)%lon)
          tg%minlat = MIN(tg%minlat,icon_grid_region%verts%vertex(i)%lat)
          tg%maxlat = MAX(tg%maxlat,icon_grid_region%verts%vertex(i)%lat)
        ENDDO

        tg%minlon = rad2deg * tg%minlon - 0.25_wp
        tg%minlat = rad2deg * tg%minlat - 0.05_wp
        tg%maxlon = rad2deg * tg%maxlon + 0.25_wp
        tg%maxlat = rad2deg * tg%maxlat + 0.05_wp

        ! Compute list for search start index; dimensions(lon,lat)
        i1s = -180*search_res
        i1e =  180*search_res
        i2s =  -90*search_res
        i2e =   90*search_res
        ALLOCATE(tg%search_index(i1s:i1e,i2s:i2e))
        tg%search_index(:,:) = 0_i4
        DO i=1,tg%ie
          i1 = NINT(lon_geo(i,j,k)*search_res)
          i2 = NINT(lat_geo(i,j,k)*search_res)
          IF (tg%search_index(i1,i2) == 0_i4) THEN
            tg%search_index(i1,i2) = i
          ELSE ! determine which point is closer to the target point
            ii = tg%search_index(i1,i2)
            d1 = (lat_geo(i,j,k) - REAL(i2,wp)/REAL(search_res,wp))**2 + &
                ( (lon_geo(i,j,k) - REAL(i1,wp)/REAL(search_res,wp))*COS(deg2rad*REAL(i2,wp)/REAL(search_res,wp)) )**2
            d2 = (lat_geo(ii,j,k) - REAL(i2,wp)/REAL(search_res,wp))**2 + &
                ( (lon_geo(ii,j,k) - REAL(i1,wp)/REAL(search_res,wp))*COS(deg2rad*REAL(i2,wp)/REAL(search_res,wp)) )**2
            IF (d1 < d2) tg%search_index(i1,i2) = i
          ENDIF
        ENDDO

        ! Fill gaps in index list
        i1s = MAX(-180*search_res,NINT(tg%minlon*search_res)-1)
        i1e = MIN( 180*search_res,NINT(tg%maxlon*search_res)+1)
        i2s = MAX( -90*search_res,NINT(tg%minlat*search_res)-1)
        i2e = MIN(  90*search_res,NINT(tg%maxlat*search_res)+1)

        DO i2 = i2s, i2e
          DO i1 = i1e,i1s+1,-1
            IF (tg%search_index(i1-1,i2) == 0 .AND. tg%search_index(i1,i2) /= 0) &
              tg%search_index(i1-1,i2) = tg%search_index(i1,i2)
          ENDDO
          IF (i1e == 180*search_res .AND. tg%search_index(i1e,i2) == 0) &
            tg%search_index(i1e,i2) = tg%search_index(i1s,i2)
          DO i1 = i1s,i1e-1
            IF (tg%search_index(i1+1,i2) == 0 .AND. tg%search_index(i1,i2) /= 0) &
              tg%search_index(i1+1,i2) = tg%search_index(i1,i2)
          ENDDO
        ENDDO
        ! Fill also empty rows adjacent to valid rows
        DO i2 = i2s+1,i2e
          IF (ALL(tg%search_index(i1s:i1e,i2-1) == 0) .AND. ANY(tg%search_index(i1s:i1e,i2) /= 0)) &
            tg%search_index(i1s:i1e,i2-1) = tg%search_index(i1s:i1e,i2)
        ENDDO
        DO i2 = i2e-1,i2s,-1
          IF (ALL(tg%search_index(i1s:i1e,i2+1) == 0) .AND. ANY(tg%search_index(i1s:i1e,i2) /= 0)) &
            tg%search_index(i1s:i1e,i2+1) = tg%search_index(i1s:i1e,i2)
        ENDDO

        !-----------------------------------------------------------------
      CASE(igrid_cosmo) ! COSMO grid  
              
        CALL calculate_cosmo_target_coordinates(tg,cosmo_grid,lon_geo,lat_geo,lon_rot,lat_rot)
        lonadj = .FALSE.
        lonadjLoop: DO i2=1,tg%je  !_br 21.02.14 changed lonadj to lonadjLoop because of name conflict
          DO i1=1,tg%ie - 1 
            IF (lon_geo(i1,i2,1) > lon_geo(i1+1,i2,1)) THEN
              lonadj=.TRUE.
              EXIT lonadjLoop  !_br 21.02.14
            ENDIF
          ENDDO
         ENDDO lonadjLoop  !_br 21.02.14
         IF (lonadj) THEN
           ALLOCATE ( auxlon(tg%ie,tg%je) )
           auxlon(:,:) = lon_geo(:,:,1)
           WHERE (auxlon < 0.) auxlon = auxlon + 360.
           tg%minlon = MINVAL(auxlon)
           tg%maxlon = MAXVAL(auxlon)
           DEALLOCATE ( auxlon )
         ELSE
           tg%minlon = MINVAL(lon_geo)
           tg%maxlon = MAXVAL(lon_geo)
         ENDIF
         tg%minlat = MINVAL(lat_geo)
         tg%maxlat = MAXVAL(lat_geo)

    END SELECT
    !-----------------------------------------------------------------

    CALL logging%info('Exit routine: init_target_grid')
  END SUBROUTINE init_target_grid

END MODULE mo_target_grid_routines
