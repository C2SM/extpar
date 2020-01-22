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
  USE mo_kind, ONLY: wp
  USE mo_kind, ONLY: i4
  USE mo_kind, ONLY: i4

  USE mo_bilinterpol, ONLY:  get_4_surrounding_raw_data_indices, &
    &                        calc_weight_bilinear_interpol, &
    &                        calc_value_bilinear_interpol
  USE mo_logging


IMPLICIT NONE

PRIVATE

PUBLIC :: agg_aot_data_to_target_grid

  CONTAINS

  !> Subroutine to aggregate aerosol optical thickness data to the target grid
  SUBROUTINE agg_aot_data_to_target_grid(iaot_type,ntime,ntype,n_spectr)
  !-------------------------------------------------------------------------------------
  ! list of modules which are used as "input"

  ! USE global data fields (coordinates)
  USE mo_target_grid_data, ONLY: lon_geo, & !< longitude coordinates of the COSMO grid in the geographical system 
    &                            lat_geo !< latitude coordinates of the COSMO grid in the geographical system
  USE mo_target_grid_data, ONLY: tg !< structure with target grid description

  ! "input" fields
  USE mo_aot_data, ONLY : lon_aot, &
    &                     lat_aot, &
    &                     aot_data, &
    &                     aot_grid, &
    &                     MAC_data !------new MACv2-----
  !-------------------------------------------------------------------------------------
  ! list of modules and fields for "output"
  USE mo_aot_target_fields, ONLY: aot_tg,&
       &                          MAC_aot_tg,&
       &                          MAC_ssa_tg,&
       &                          MAC_asy_tg


  INTEGER (KIND=i4), INTENT(IN) :: iaot_type !< type of AOT source data
  INTEGER (KIND=i4), INTENT(IN) :: ntype !< number of types of aerosols
  INTEGER (KIND=i4), INTENT(IN) :: ntime !< number of times
  INTEGER (KIND=i4), INTENT(IN) :: n_spectr !< number of spectr new
  INTEGER :: i,j,k,l ! counters

  REAL (KIND=wp) :: point_lon_geo       !< longitude coordinate in geographical system of input point 
  REAL (KIND=wp) :: point_lat_geo       !< latitude coordinate in geographical system of input point
      
       
  INTEGER (KIND=i4) :: western_column     !< the index of the western_column of raw data 
  INTEGER (KIND=i4) :: eastern_column     !< the index of the eastern_column of raw data 
  INTEGER (KIND=i4) :: northern_row       !< the index of the northern_row of raw data 
  INTEGER (KIND=i4) :: southern_row       !< the index of the southern_row of raw data 


  ! global data flag
  LOGICAL :: gldata=.TRUE. ! AOT data are global


  REAL (KIND=wp) :: bwlon !< weight for bilinear interpolation
  REAL (KIND=wp) :: bwlat !< weight for bilinear interpolation
  REAL (KIND=wp) :: bwlon2d(ntime,ntype)
  REAL (KIND=wp) :: bwlat2d(ntime,ntype)

  REAL (KIND=wp) :: data_array_sw(ntime,ntype) !< data array values at south-western point
  REAL (KIND=wp) :: data_array_se(ntime,ntype) !< data array values at south-eastern point
  REAL (KIND=wp) :: data_array_ne(ntime,ntype) !< data array values at north-eastern point
  REAL (KIND=wp) :: data_array_nw(ntime,ntype) !< data array values at north-western point
  REAL (KIND=wp) :: target_array_value(ntime,ntype) !< interpolated values

  target_array_value = -999.

  IF (verbose >= idbg_low ) THEN
    WRITE(logging%fileunit,*)'MINVAL(lon_geo)', MINVAL(lon_geo)
    WRITE(logging%fileunit,*)'MAXVAL(lon_geo)', MAXVAL(lon_geo)
    WRITE(logging%fileunit,*)'MINVAL(lat_geo)', MINVAL(lat_geo)
    WRITE(logging%fileunit,*)'MAXVAL(lat_geo)', MAXVAL(lat_geo)
  ENDIF

  ! loop through all target grid elements
  DO i=1,tg%ie
    DO j=1,tg%je
      DO k=1,tg%ke

        point_lon_geo = lon_geo(i,j,1) 
        point_lat_geo = lat_geo(i,j,1)

        ! get four surrounding raw data indices
        !PRINT *,'CALL  get_4_surrounding_raw_data_indices'

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
        ! for the bilinear interpolation, use the four data pixels
         ! data(western_column,northern_row)
         ! data(western_column,southern_row)
         ! data(eastern_column,southern_row)
         ! data(eastern_column,northern_row)

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

          IF (iaot_type ==4) THEN
            DO l=1,n_spectr
              data_array_sw(1:ntime,1:ntype)= MAC_data(western_column,southern_row,l,1:ntime,1:ntype)
              data_array_se(1:ntime,1:ntype)= MAC_data(eastern_column,southern_row,l,1:ntime,1:ntype)
              data_array_ne(1:ntime,1:ntype)= MAC_data(eastern_column,northern_row,l,1:ntime,1:ntype)
              data_array_nw(1:ntime,1:ntype)= MAC_data(western_column,northern_row,l,1:ntime,1:ntype)
              target_array_value = calc_value_bilinear_interpol(bwlon2d,bwlat2d,&
                  data_array_sw, data_array_se, data_array_ne, data_array_nw)
              MAC_aot_tg(i,j,l,1:ntime)=target_array_value(1:ntime,1)
              MAC_ssa_tg(i,j,l,1:ntime)=target_array_value(1:ntime,2)
              MAC_asy_tg(i,j,l,1:ntime)=target_array_value(1:ntime,3)
            ENDDO
          ELSE
            data_array_sw(1:ntime,1:ntype) = aot_data(western_column,southern_row,1:ntime,1:ntype) 
            data_array_se(1:ntime,1:ntype) = aot_data(eastern_column,southern_row,1:ntime,1:ntype)
            data_array_ne(1:ntime,1:ntype) = aot_data(eastern_column,northern_row,1:ntime,1:ntype)
            data_array_nw(1:ntime,1:ntype) = aot_data(western_column,northern_row,1:ntime,1:ntype)
            target_array_value = calc_value_bilinear_interpol(bwlon2d,bwlat2d,&
                   data_array_sw, data_array_se, data_array_ne, data_array_nw)
            aot_tg(i,j,k,1:ntype,1:ntime) = TRANSPOSE(target_array_value(1:ntime,1:ntype))
          ENDIF
        ELSE
          IF (iaot_type == 4) THEN
            DO l = 1,n_spectr
              MAC_aot_tg(i,j,l,1:ntime)=target_array_value(1:ntime,1)
              MAC_ssa_tg(i,j,l,1:ntime)=target_array_value(1:ntime,2)
              MAC_asy_tg(i,j,l,1:ntime)=target_array_value(1:ntime,3)
            ENDDO
          ELSE
            aot_tg(i,j,k,1:ntype,1:ntime) = TRANSPOSE(target_array_value(1:ntime,1:ntype))
          ENDIF
        ENDIF

      ENDDO
    ENDDO
  ENDDO ! loop through all target grid elements

  END SUBROUTINE agg_aot_data_to_target_grid

END MODULE mo_agg_aot
