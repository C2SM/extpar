!+ Fortran module with netcdf output routines for FLAKE data on the target grid
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_0         2010/12/21 Hermann Asensio
!  Initial release
! V1_4         2011/04/21 Hermann Asensio
!  clean up
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!> Fortran module with netcdf output routines for FLAKE data on the target grid
!> ouptut routines
!> \author Hermann Asensio
MODULE mo_flake_output_nc

  USE mo_logging
  USE mo_kind,                  ONLY: wp, i4

  USE mo_grid_structures,       ONLY: target_grid_def

  USE mo_io_utilities,          ONLY: netcdf_attributes, &
       &                              dim_meta_info, &
       &                              netcdf_put_var, &
       &                              open_new_netcdf_file, &
       &                              netcdf_get_var, &
       &                              close_netcdf_file

  USE mo_var_meta_data,         ONLY: dim_3d_tg, &
      &                               def_dimension_info_buffer, &
      &                               lon_geo_meta, &
      &                               lat_geo_meta, &
      &                               def_flake_fields_meta, &
      &                               def_com_target_fields_meta, &
      &                               lake_depth_meta, fr_lake_meta, &
      &                               flake_tot_npixel_meta

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: write_netcdf_buffer_flake, &
       &    read_netcdf_buffer_flake


  CONTAINS

    !> netcdf output of FLAKE derived buffer fields
  SUBROUTINE write_netcdf_buffer_flake(netcdf_filename,  &
  &                                     tg,         &
  &                                     undefined, &
  &                                     undef_int,   &
  &                                     lon_geo,     &
  &                                     lat_geo, &
  &                                     lake_depth, &
  &                                     fr_lake,    &
  &                                     flake_tot_npixel)


  CHARACTER (len=*), INTENT(IN)      :: netcdf_filename !< filename for the netcdf file
  TYPE(target_grid_def), INTENT(IN)  :: tg !< structure with target grid description
  INTEGER(KIND=i4), INTENT(IN)       :: undef_int, &       !< value to indicate undefined grid elements
       &                                flake_tot_npixel(:,:,:)  
  REAL (KIND=wp), INTENT(IN)         :: lon_geo(:,:,:), &  !< longitude coordinates of the target grid in the geographical system
       &                                undefined, &       !< value to indicate undefined grid elements 
       &                                lat_geo(:,:,:), &  !< latitude coordinates of the target grid in the geographical system
       &                                lake_depth(:,:,:), & !< lake depth
                                        fr_lake(:,:,:)     !< fraction of fresh water (lakes)

  ! local variables
  INTEGER (KIND=i4)                  :: undefined_i, ndims, ncid, errorcode
  INTEGER, PARAMETER                 :: nglob_atts=6

  TYPE(dim_meta_info), ALLOCATABLE   :: dim_list(:) !< dimensions for netcdf file
  TYPE(netcdf_attributes)            :: global_attributes(nglob_atts)

  CALL logging%info('Enter routine: write_netcdf_buffer_flake')

  !-------------------------------------------------------------
  ! define global attributes
  CALL set_global_att_flake(global_attributes)

  !set up dimensions for buffer
  CALL  def_dimension_info_buffer(tg)
  ! dim_3d_tg

  ! define meta information for various land use related variables (FLAKE) for netcdf output
  CALL def_flake_fields_meta(dim_3d_tg)

  ! lake_depth_meta, fr_lake_meta, &
  !  &       flake_tot_npixel_meta

  ! define meta information for target field variables lon_geo, lat_geo 
  CALL def_com_target_fields_meta(dim_3d_tg)
  ! lon_geo_meta and lat_geo_meta
  !set up dimensions for buffer netcdf output 
  ndims = 3
  undefined_i = undef_int
  ALLOCATE(dim_list(1:ndims),STAT=errorcode)
  IF (errorcode /= 0 ) CALL logging%error('Cant allocate array dim_list',__FILE__,__LINE__)
  dim_list = dim_3d_tg

  !-----------------------------------------------------------------

  CALL open_new_netcdf_file(netcdf_filename=TRIM(netcdf_filename),   &
      &                       dim_list=dim_list,                  &
      &                       global_attributes=global_attributes, &
      &                       ncid=ncid)

  !-----------------------------------------------------------------
  ! 3D variables flake_tot_npixel
  ! lon
  CALL netcdf_put_var(ncid,lon_geo,lon_geo_meta,undefined)

  ! lat
  CALL netcdf_put_var(ncid,lat_geo,lat_geo_meta,undefined)

  ! fr_lake
  CALL netcdf_put_var(ncid,fr_lake,fr_lake_meta,undefined)

  ! lake_depth
  CALL netcdf_put_var(ncid,lake_depth,lake_depth_meta,undefined)

  ! flake_tot_npixel
  CALL netcdf_put_var(ncid,flake_tot_npixel,flake_tot_npixel_meta,undefined_i)


  CALL close_netcdf_file(ncid)

  CALL logging%info('Exit routine: write_netcdf_buffer_flake')

  END SUBROUTINE write_netcdf_buffer_flake

  !-----------------------------------------------------------------------
  !> set global attributes for netcdf with flake data
  SUBROUTINE set_global_att_flake(global_attributes)
    TYPE(netcdf_attributes), INTENT(INOUT) :: global_attributes(1:6)

    !local variables
    CHARACTER(len=10) :: ydate
    CHARACTER(len=10) :: ytime
    CHARACTER(len=2)  :: cc
    CHARACTER(len=2)  :: yy
    CHARACTER(len=2)  :: mm
    CHARACTER(len=2)  :: dd
    CHARACTER(len=2)  :: hh
    CHARACTER(len=2)  :: minute

    ! define global attributes
    global_attributes(1)%attname = 'title'
    global_attributes(1)%attributetext='Land Use data'
    global_attributes(2)%attname = 'institution'
    global_attributes(2)%attributetext='Deutscher Wetterdienst'

    global_attributes(3)%attname = 'source'
    global_attributes(3)%attributetext='FLAKE data'

    CALL DATE_AND_TIME(ydate,ytime)
    READ(ydate,'(4A2)') cc,yy,mm,dd
    READ(ytime,'(2A2)') hh, minute

    ydate=TRIM(cc)//TRIM(yy)//'-'//TRIM(mm)//'-'//TRIM(dd)
    ytime=TRIM(hh)//':'//TRIM(minute) 

    global_attributes(4)%attname = 'history'
    global_attributes(4)%attributetext=TRIM(ydate)//'T'//TRIM(ytime)//' flake_to_buffer'

    global_attributes(5)%attname = 'references'
    global_attributes(5)%attributetext=''

    global_attributes(6)%attname = 'comment'
    global_attributes(6)%attributetext=''

  END SUBROUTINE set_global_att_flake
  !----------------------------------------------------------------------------

  !> read FLAKE derived buffer fields
  SUBROUTINE read_netcdf_buffer_flake(netcdf_filename,  &
    &                                     lake_depth, &
    &                                     fr_lake,    &
    &                                     flake_tot_npixel)

    CHARACTER (len=*), INTENT(IN)      :: netcdf_filename !< filename for the netcdf file
    REAL (KIND=wp), INTENT(OUT)        :: lake_depth(:,:,:), & !< lake depth
         &                                fr_lake(:,:,:)     !< fraction of fresh water (lakes)

    INTEGER (KIND=i4), INTENT(OUT)     :: flake_tot_npixel(:,:,:)

    CALL logging%info('Enter routine: read_netcdf_buffer_flake')

    ! define meta information for various land use related variables (FLAKE) for netcdf output
    CALL def_flake_fields_meta(dim_3d_tg)

    ! define meta information for target field variables lon_geo, lat_geo 
    CALL def_com_target_fields_meta(dim_3d_tg)

    CALL netcdf_get_var(TRIM(netcdf_filename),lake_depth_meta,lake_depth)

    CALL netcdf_get_var(TRIM(netcdf_filename),flake_tot_npixel_meta,flake_tot_npixel)

    CALL netcdf_get_var(TRIM(netcdf_filename),fr_lake_meta,fr_lake)

    CALL logging%info('Exit routine: read_netcdf_buffer_flake')

  END SUBROUTINE read_netcdf_buffer_flake
 
END Module mo_flake_output_nc
