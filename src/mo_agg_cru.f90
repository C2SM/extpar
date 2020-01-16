!+  Fortran module to aggregate the CRU near surface climatology to target grid
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_0         2010/12/21 Hermann Asensio
!  Initial release
! V1_7         2013/01/25 Guenther Zaengl 
!   Parallel threads for ICON and COSMO using Open-MP, 
!   Several bug fixes and optimizations for ICON search algorithm, 
!   particularly for the special case of non-contiguous domains; 
!   simplified namelist control for ICON  
! V2_0         2013/06/04 Martina Messmer
!  adaptations to read also the higher resolved CRU data set (CLM Community)
!  new parameter for the CRU temperature elevation is introduced
!  aggregation of new parameter is performed in the same way as temperature
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!> Fortran module to aggregate the CRU near surface climatology to target grid
!!
!> \author Hermann Asensio
MODULE mo_agg_cru

  !> kind parameters are defined in MODULE data_parameters
  USE mo_kind, ONLY: wp
  USE mo_kind, ONLY: i8

  USE mo_bilinterpol, ONLY:  get_4_surrounding_raw_data_indices, &
    &                        calc_weight_bilinear_interpol, &
    &                        calc_value_bilinear_interpol

  USE mo_cru_target_fields, ONLY: i_t_cru_fine, i_t_cru_coarse

IMPLICIT NONE

PRIVATE

PUBLIC :: agg_cru_data_to_target_grid


CONTAINS

  !> Subroutine to aggregate CRU temperature data to the target grid
  SUBROUTINE agg_cru_data_to_target_grid(nrows,ncolumns,ntime,raw_data_t_id)
  !-------------------------------------------------------------------------------------
     ! USE global data fields (coordinates)
     USE mo_target_grid_data, ONLY: lon_geo, & !< longitude coordinates of the COSMO grid in the geographical system 
                                     lat_geo !< latitude coordinates of the COSMO grid in the geographical system
     USE mo_target_grid_data, ONLY: tg  !< structure with target grid description


     ! "input" fields
     USE mo_cru_data, ONLY : lon_cru,      &
       &                     lat_cru,      &
       &                     cru_raw_data, &
       &                     cru_raw_elev, &
       &                     cru_grid

     !-------------------------------------------------------------------------------------
     ! list of modules and fields for "output"

     USE mo_cru_target_fields, ONLY: crutemp, cruelev

                                  
     INTEGER (KIND=i8), INTENT(IN) :: nrows !< number of rows
     INTEGER (KIND=i8), INTENT(IN) :: ncolumns !< number of columns
     INTEGER (KIND=i8), INTENT(IN) :: ntime !< number of times
     INTEGER (KIND=i8), INTENT(IN) :: raw_data_t_id !< integer switch to decide which data set must be used. (CRU fine, CRU coarse)


     INTEGER :: i,j,k, t ! counters

      REAL (KIND=wp) :: point_lon_geo       !< longitude coordinate in geographical system of input point 
      REAL (KIND=wp) :: point_lat_geo       !< latitude coordinate in geographical system of input point
      
       
      INTEGER (KIND=i8) :: western_column     !< the index of the western_column of raw data 
      INTEGER (KIND=i8) :: eastern_column     !< the index of the eastern_column of raw data 
      INTEGER (KIND=i8) :: northern_row       !< the index of the northern_row of raw data 
      INTEGER (KIND=i8) :: southern_row       !< the index of the southern_row of raw data 

      REAL (KIND=wp) :: data_sw
      REAL (KIND=wp) :: data_se
      REAL (KIND=wp) :: data_ne
      REAL (KIND=wp) :: data_nw

      REAL (KIND=wp) :: target_value

      ! global data flag
      LOGICAL :: gldata=.TRUE. ! CRU data are global

      REAL (KIND=wp) :: bwlon !< weight for bilinear interpolation
      REAL (KIND=wp) :: bwlat !< weight for bilinear interpolation

      REAL (KIND=wp) :: tem_clim_raw(ncolumns,nrows)
      REAL (KIND=wp) :: elev_clim_raw(ncolumns,nrows)

      INTEGER, PARAMETER :: dpm(12) = [ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

       tem_clim_raw = 0.0
       DO j=1, nrows
       DO i=1, ncolumns

       SELECT CASE(raw_data_t_id)
       CASE(i_t_cru_coarse)
        DO t=1, ntime
          tem_clim_raw(i,j) = tem_clim_raw(i,j) +  dpm(t) * cru_raw_data(i,j,t)
        ENDDO

         tem_clim_raw(i,j) = 273.15 +  tem_clim_raw(i,j) / 365 ! unit in K instead degC, and yearly mean instead of monthly means

 

     CASE(i_t_cru_fine)
         DO t=1, ntime
         tem_clim_raw(i,j)  =  cru_raw_data(i,j,t)
         elev_clim_raw(i,j) =  cru_raw_elev(i,j,t)
       ENDDO
     END SELECT

   ENDDO
 ENDDO
      


 


      ! loop through all target grid elements
     DO i=1,tg%ie
     DO j=1,tg%je
     DO k=1,tg%ke

!     DO i=200,202
!     DO j=300,301
!     DO k=1,tg%ke


     point_lon_geo = lon_geo(i,j,k) 
     point_lat_geo = lat_geo(i,j,k)

     ! get four surrounding raw data indices

       CALL  get_4_surrounding_raw_data_indices(   cru_grid, &
                                                   lon_cru,           &
                                                   lat_cru,           &
                                                   gldata,            &
                                                   point_lon_geo,      &
                                                   point_lat_geo,      &
                                                   western_column,     &
                                                   eastern_column,     &
                                                   northern_row,       &
                                                   southern_row)
     ! for the bilinear interpolation, use the four data pixels
     ! data(western_column,northern_row)
     ! data(western_column,southern_row)
     ! data(eastern_column,southern_row)
     ! data(eastern_column,northern_row)

   ! PRINT *,'western_column: ', western_column
   ! PRINT *,'eastern_column: ', eastern_column
   ! PRINT *,'northern_row: ', northern_row
   ! PRINT *,'southern_row: ', southern_row



     ! calculate weight for bilinear interpolation
    target_value = -999.
    if ((western_column /= 0) ) then  ! point is not out of data grid range

     CALL calc_weight_bilinear_interpol(point_lon_geo, &
                                        point_lat_geo, &
                                        lon_cru(western_column),      &
                                        lon_cru(eastern_column),      &
                                        lat_cru(northern_row),     &
                                        lat_cru(southern_row),     &
                                        bwlon,         &
                                        bwlat)
     ! the weights are bwlon and bwlat
   !  PRINT *,'bwlon, bwlat', bwlon,bwlat
     ! put all relevent data to an array



     data_sw = tem_clim_raw(western_column,southern_row) 
     data_se = tem_clim_raw(eastern_column,southern_row)
     data_ne = tem_clim_raw(eastern_column,northern_row)
     data_nw = tem_clim_raw(western_column,northern_row)

     !  missing over water is 0

     ! perform the interpolation
     SELECT CASE(raw_data_t_id)
     CASE(i_t_cru_fine)
       if (data_sw .gt. 0 .and. data_se .gt. 0 .and.                 &
           data_ne .gt. 0 .and. data_nw .gt. 0) then 

           target_value = calc_value_bilinear_interpol(bwlon, bwlat, &
                                        data_sw, data_se, data_ne, data_nw)
       else
          !interpolation not possible as missing
!          print *, data_sw, data_se, data_ne, data_nw
       endif
     CASE(i_t_cru_coarse) 
       target_value = calc_value_bilinear_interpol(bwlon, bwlat, &
                                       data_sw, data_se, data_ne, data_nw)
     END SELECT

     ELSE
    !  PRINT *,'point_lon_geo: ', point_lon_geo
    !  PRINT *,'point_lat_geo: ', point_lat_geo
    !  PRINT *,'HA debug: western_column: ', western_column
    !  PRINT *,'HA debug: eastern_column: ', eastern_column
    !  PRINT *,'northern_row: ', northern_row
    !  PRINT *,'southern_row: ', southern_row
!
!
!      PRINT *,'Grid element i,j,k: ', i,j,k
!      PRINT *,' lon_geo(i,j,k) : ', lon_geo(i,j,k)
!      PRINT *,' lat_geo(i,j,k) : ', lat_geo(i,j,k)

     ENDIF
     crutemp(i,j,k) = target_value

     ENDDO
     ENDDO
     ENDDO ! loop through all target grid elements

     ! ELEVATION
     SELECT CASE(raw_data_t_id)
     CASE(i_t_cru_fine)
     ! loop through all target grid elements
       DO i=1,tg%ie
         DO j=1,tg%je
           DO k=1,tg%ke

             point_lon_geo = lon_geo(i,j,k)
             point_lat_geo = lat_geo(i,j,k)

             CALL  get_4_surrounding_raw_data_indices(   cru_grid, &
                                                   lon_cru,           &
                                                   lat_cru,           &
                                                   gldata,            &
                                                   point_lon_geo,      &
                                                   point_lat_geo,      &
                                                   western_column,     &
                                                   eastern_column,     &
                                                   northern_row,       &
                                                   southern_row)

             ! calculate weight for bilinear interpolation
             target_value = -999.
             if ((western_column /= 0) ) then  ! point is not out of data grid range
               
               
               CALL calc_weight_bilinear_interpol(point_lon_geo, &
                                        point_lat_geo, &
                                        lon_cru(western_column),      &
                                        lon_cru(eastern_column),      &
                                        lat_cru(northern_row),     &
                                        lat_cru(southern_row),     &
                                        bwlon,         &
                                        bwlat)

               data_sw = elev_clim_raw(western_column,southern_row)
               data_se = elev_clim_raw(eastern_column,southern_row)
               data_ne = elev_clim_raw(eastern_column,northern_row)
               data_nw = elev_clim_raw(western_column,northern_row)

               ! perform the interpolation
               
               target_value = calc_value_bilinear_interpol(bwlon, bwlat, &
                              data_sw, data_se, data_ne, data_nw)
             ELSE

             ENDIF

             cruelev(i,j,k) = target_value
             
           ENDDO
         ENDDO
       ENDDO ! loop through all target grid elements

    END SELECT


    END SUBROUTINE agg_cru_data_to_target_grid

END MODULE mo_agg_cru
