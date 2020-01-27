!+ Fortran module for netcdf output of DEM subgrid-scale slope data
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! Initial release
! V4_0         2016/07/28 Daniel Luethi
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!> Fortran module for netcdf output of DEM subgrid-scale slope data
!> \author Hermann Asensio
MODULE mo_sgsl_output_nc
  
  USE mo_logging
  USE mo_kind,                  ONLY: wp, i4

  USE mo_grid_structures,       ONLY: target_grid_def, &
       &                              igrid_icon, &
       &                              igrid_cosmo

  USE mo_cosmo_grid,            ONLY: cosmo_grid
  USE mo_icon_grid_data,        ONLY: ICON_grid

  USE mo_sgsl_data,             ONLY: idem_type, &
       &                              dem_gl, &
       &                              dem_aster

  USE mo_io_utilities,          ONLY: netcdf_attributes, &
       &                              dim_meta_info, &
       &                              netcdf_put_var, &
       &                              open_new_netcdf_file, &
       &                              netcdf_get_var, &
       &                              close_netcdf_file

  USE mo_var_meta_data,         ONLY: dim_3d_tg, dim_4d_tg,    &
       &                              def_dimension_info_buffer, &
       &                              lon_geo_meta, &
       &                              lat_geo_meta, &
       &                              def_sgsl_meta, &
       &                              sgsl_meta, &
       &                              def_sgsl_meta, &
       &                              def_com_target_fields_meta  
    

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: write_netcdf_buffer_sgsl

  PUBLIC :: read_netcdf_buffer_sgsl

  CONTAINS

    !> create a netcdf file for the fields derived from GLOBE data to the buffer 
  SUBROUTINE write_netcdf_buffer_sgsl(netcdf_filename,&
     &                                  tg,            &
     &                                  undefined,     &
     &                                  undef_int,     &
     &                                  igrid_type,    &
     &                                  lon_geo,       &
     &                                  lat_geo,       &
     &                                  sgsl)

   
     TYPE(target_grid_def)              :: tg !< structure with target grid description
     CHARACTER (len=*), INTENT(IN)      :: netcdf_filename !< filename for the netcdf file

     REAL(KIND=wp), INTENT(IN)          :: undefined, &       !< value to indicate undefined grid elements 
          &                                lon_geo(:,:,:), &  !< longitude coordinates of the target grid in the geographical system
          &                                lat_geo(:,:,:), &  !< latitude coordinates of the target grid in the geographical system
          &                                sgsl(:,:,:)  !< subgrid-scale slope parameter 

     INTEGER(KIND=i4), INTENT(IN)       :: undef_int, &       !< value to indicate undefined grid elements
          &                                igrid_type

     ! local variables
     INTEGER (KIND=i4)                  :: istart, iend, jstart, jend, &
         &                                ndims, & 
         &                                ncid, &
         &                                errorcode !< error status variable

     INTEGER(KIND=i4), PARAMETER        :: nglob_atts=6

     TYPE(dim_meta_info), ALLOCATABLE   :: dim_list(:) !< dimensions for netcdf file
     TYPE(netcdf_attributes)            :: global_attributes(nglob_atts)

    CALL logging%info('Enter routine: write_netcdf_buffer_sgsl')

    !-------------------------------------------------------------
    ! define global attributes
    CALL set_global_att_sgsl(global_attributes)

     !set up dimensions for buffer
    CALL  def_dimension_info_buffer(tg)
    ! dim_3d_tg, dim_4d_tg
    
    ! define meta information for target field variables lon_geo, lat_geo 
    CALL def_com_target_fields_meta(dim_3d_tg)
    ! lon_geo_meta and lat_geo_meta
    ! define meta information for various GLOBE data related variables for netcdf output
    CALL def_sgsl_meta(dim_3d_tg,idem_type)
    !set up dimensions for buffer netcdf output 
    ndims = SIZE(dim_3d_tg)
    ALLOCATE(dim_list(1:ndims),STAT=errorcode)
    IF (errorcode /= 0 ) CALL logging%error('Cant allocate array dim_list',__FILE__,__LINE__)
    dim_list = dim_3d_tg

    !-----------------------------------------------------------------
    CALL open_new_netcdf_file(netcdf_filename=TRIM(netcdf_filename), &
        &                       dim_list=dim_list,                   &
        &                       global_attributes=global_attributes, &
        &                       ncid=ncid)

    ! correct start and stop indices if needed
    SELECT CASE (igrid_type)
      CASE(igrid_cosmo)
        istart = 1
        jstart = 1
        iend   = cosmo_grid%nlon_rot
        jend   = cosmo_grid%nlat_rot      
      CASE(igrid_icon)
        istart = 1
        jstart = 1
        iend   = icon_grid%ncell
        jend   = 1
    END SELECT

      ! lon
    CALL netcdf_put_var(ncid,lon_geo(istart:iend,jstart:jend,:),lon_geo_meta,undefined)

      ! lat
    CALL netcdf_put_var(ncid,lat_geo(istart:iend,jstart:jend,:),lat_geo_meta,undefined)

    ! sgsl
    CALL netcdf_put_var(ncid,sgsl(istart:iend,jstart:jend,:),sgsl_meta,undefined)

    CALL close_netcdf_file(ncid)
     
    CALL logging%info('Exit routine: write_netcdf_buffer_sgsl')

  END SUBROUTINE write_netcdf_buffer_sgsl

  SUBROUTINE set_global_att_sgsl(global_attributes)

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
    SELECT CASE(idem_type)
      CASE(dem_aster)
        global_attributes(1)%attributetext='ASTER data '
      CASE(dem_gl)
        global_attributes(1)%attributetext='GLOBE data '
      END SELECT
    global_attributes(2)%attname = 'institution'
    global_attributes(2)%attributetext='Deutscher Wetterdienst'

    global_attributes(3)%attname = 'source'
    SELECT CASE(idem_type)
      CASE(dem_aster)
        global_attributes(3)%attributetext='ASTER,The Advanced Spaceborne Thermal Emission '// & !_br 21.02.14 splitted too long line
        & 'and Reflection Radiometer, 1 arc-second digital elevation model' !_br 21.02.14
      CASE(dem_gl)
        global_attributes(3)%attributetext='GLOBE, Global Land One-km Base Elevation'
      END SELECT
    

    CALL DATE_AND_TIME(ydate,ytime)
    READ(ydate,'(4A2)') cc,yy,mm,dd
    READ(ytime,'(2A2)') hh, minute

    ydate=TRIM(cc)//TRIM(yy)//'-'//TRIM(mm)//'-'//TRIM(dd)
    ytime=TRIM(hh)//':'//TRIM(minute) 

    global_attributes(4)%attname = 'history'
    global_attributes(4)%attributetext=TRIM(ydate)//'T'//TRIM(ytime)//' sgsl_to_buffer'

    global_attributes(5)%attname = 'references'
    SELECT CASE(idem_type)
      CASE(dem_aster)
        global_attributes(5)%attributetext='http://www.jspacesystems.or.jp/ersdac/GDEM/E/4.html'
      CASE(dem_gl)
        global_attributes(5)%attributetext='http://www.ngdc.noaa.gov/mgg/topo/globe.html'
      END SELECT

    global_attributes(6)%attname = 'comment'
    global_attributes(6)%attributetext=''

  END SUBROUTINE set_global_att_sgsl
  !-----------------------------------------------------------------------

  !> read netcdf file for the fields derived from GLOBE data from the buffer 
   SUBROUTINE read_netcdf_buffer_sgsl(netcdf_filename,&
     &                                 tg,             &
     &                                 sgsl)
   
     CHARACTER (len=*), INTENT(IN)         :: netcdf_filename !< filename for the netcdf file
     TYPE(target_grid_def), INTENT(IN)     :: tg !< structure with target grid description

     REAL(KIND=wp), INTENT(OUT)            :: sgsl(:,:,:)  !< mean height 

     CALL logging%info('Enter routine: read_netcdf_buffer')
     
     !set up dimensions for buffer
     CALL  def_dimension_info_buffer(tg)
     ! define meta information for target field variables lon_geo, lat_geo 
     CALL def_com_target_fields_meta(dim_3d_tg)
     ! lon_geo_meta and lat_geo_meta
     ! define meta information for various GLOBE data related variables for netcdf output
     CALL def_sgsl_meta(dim_3d_tg,idem_type)
     !set up dimensions for buffer netcdf output 

     CALL netcdf_get_var(TRIM(netcdf_filename),sgsl_meta,sgsl)

  END SUBROUTINE read_netcdf_buffer_sgsl

END MODULE mo_sgsl_output_nc

