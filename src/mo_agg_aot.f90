!+ Fortran module to aggregate aerosol optical thickness raw data to target grid
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_0         2010/12/21 Hermann Asensio
!  Initial release
! V1_3         2011/04/19 Hermann Asensio
!  suppurt unlimited time dimension for netcdf
! V4_0         2016/08/17 authors from RHM and Daniel Lüthi
!  Add support for MACv2 aerosol fields
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!> Fortran module to aggregate aerosol optical thickness raw data to target grid
!!
!! This module interpolates monthly aerosol optical thicknesses for five different type of aerosols 
!! <ul>
!!   <li> black carbon </li>
!!   <li> dust </li>
!!   <li> organic </li>
!!   <li> SO4 </li>
!!   <li> sea salt </li>
!! </ul>
!! from a global climatology from Ina Tegen (Tegen et al. 1997) to a target grid (COSMO/ICON). 
!! The raw data and the describing paper are available at NASA/GISS at the Global Aerosol Climatology Project 
!! (GACP http://gacp.giss.nasa.gov/data_sets/transport/). 
!!
!!
!! Tegen, I., P. Hollrigl, M. Chin, I. Fung, D. Jacob, and J. Penner 1997.
!!  <a href="http://pubs.giss.nasa.gov/abstracts/1997/Tegen_etal.html">
!!  Contribution of different aerosol species to the global aerosol extinction optical thickness:&
!!& Estimates from model results</a>.
!! J. Geophys. Res., <b>102</b>, 23895-23915.
!> \author Hermann Asensio
MODULE mo_agg_aot

  !> kind parameters are defined in MODULE data_parameters
  USE mo_logging
  USE mo_kind,                  ONLY: wp, i4

  USE mo_bilinterpol,           ONLY: get_4_surrounding_raw_data_indices, &
    &                                 calc_weight_bilinear_interpol, &
    &                                 calc_value_bilinear_interpol

  USE mo_target_grid_data,      ONLY: lon_geo, & 
    &                                 lat_geo, &
    &                                 tg

  USE mo_aot_data,              ONLY :lon_aot, &
    &                                 lat_aot, &
    &                                 aot_data, &
    &                                 aot_grid

  USE mo_aot_target_fields,     ONLY: aot_tg

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: agg_aot_data_to_target_grid

  CONTAINS

  !> Subroutine to aggregate aerosol optical thickness data to the target grid
  SUBROUTINE agg_aot_data_to_target_grid(ntime,ntype)
  !-------------------------------------------------------------------------------------

    INTEGER (KIND=i4), INTENT(IN) :: ntype, & !< number of types of aerosols
      &                              ntime!< number of times

    REAL (KIND=wp)                :: point_lon_geo, &       !< longitude coordinate in geographical system of input point 
      &                              point_lat_geo, &       !< latitude coordinate in geographical system of input point
      &                              data_array_sw(ntime,ntype), & !< data array values at south-western point
      &                              data_array_se(ntime,ntype), & !< data array values at south-eastern point
      &                              data_array_ne(ntime,ntype), & !< data array values at north-eastern point
      &                              data_array_nw(ntime,ntype), & !< data array values at north-western point
      &                              target_array_value(ntime,ntype), & !< interpolated values                        
      &                              bwlon, &
      &                              bwlat, &
      &                              bwlon2d(ntime,ntype),&
      &                              bwlat2d(ntime,ntype)


    INTEGER (KIND=i4)             :: western_column, &     !< the index of the western_column of raw data 
      &                              eastern_column, &     !< the index of the eastern_column of raw data 
      &                              northern_row, &       !< the index of the northern_row of raw data 
      &                              southern_row, &       !< the index of the southern_row of raw data 
      &                              i,j,k,l ! counters

    ! global data flag
    LOGICAL                       :: gldata=.TRUE. ! AOT data are global
    CALL logging%info('Enter routine: agg_aot_data_to_target_grid')

    target_array_value = -999.

    ! loop through all target grid elements
    DO i=1,tg%ie
      DO j=1,tg%je
        DO k=1,tg%ke
          point_lon_geo = lon_geo(i,j,1) 
          point_lat_geo = lat_geo(i,j,1)

          CALL  get_4_surrounding_raw_data_indices( aot_grid,      &
                                                    lon_aot,       &
                                                    lat_aot,       &
                                                    gldata,        &
                                                    point_lon_geo, &
                                                    point_lat_geo, &
                                                    western_column,&
                                                    eastern_column,&
                                                    northern_row,  &
                                                    southern_row)

          ! calculate weight for bilinear interpolation
          target_array_value = -999.
          IF ((western_column /= 0) ) THEN  ! point is not out of data grid range

            CALL calc_weight_bilinear_interpol(point_lon_geo, &
                                            point_lat_geo, &
                                            lon_aot(western_column),      &
                                            lon_aot(eastern_column),      &
                                            lat_aot(northern_row),     &
                                            lat_aot(southern_row),     &
                                            bwlon,         &
                                            bwlat)
            ! the weights are bwlon and bwlat
            ! put all relevant data to an array, dimension (ntime, ntype) for each grid point  
            ! perform the interpolation
            bwlon2d = bwlon
            bwlat2d = bwlat

            data_array_sw(1:ntime,1:ntype) = aot_data(western_column,southern_row,1:ntime,1:ntype) 
            data_array_se(1:ntime,1:ntype) = aot_data(eastern_column,southern_row,1:ntime,1:ntype)
            data_array_ne(1:ntime,1:ntype) = aot_data(eastern_column,northern_row,1:ntime,1:ntype)
            data_array_nw(1:ntime,1:ntype) = aot_data(western_column,northern_row,1:ntime,1:ntype)
            target_array_value = calc_value_bilinear_interpol(bwlon2d,bwlat2d,&
                    data_array_sw, data_array_se, data_array_ne, data_array_nw)
            aot_tg(i,j,k,1:ntype,1:ntime) = TRANSPOSE(target_array_value(1:ntime,1:ntype))
          ELSE
            aot_tg(i,j,k,1:ntype,1:ntime) = TRANSPOSE(target_array_value(1:ntime,1:ntype))
          ENDIF

        ENDDO
      ENDDO
    ENDDO ! loop through all target grid elements

    CALL logging%info('Exit routine: agg_aot_data_to_target_grid')

  END SUBROUTINE agg_aot_data_to_target_grid

END MODULE mo_agg_aot
