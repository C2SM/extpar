!+ Fortran module with netcdf output routines for GLC2000 data on the target grid
!
!
! Description:
! Fortran module with netcdf output routines for GLC2000 data on the target grid
!
! Current Code Owner: DWD, <Name of person responsible for this code>
!    <smail, phone, fax and email>
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V2_4         2015-05-21 Hendrik Wouters
!  Initial release based on mo_landuse_output_nc.f90 V1_14
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!=======================================================================
!> Fortran module with netcdf output routines for ISA data on the target grid
!> ouptut routines
!> \author Hendrik Wouters



MODULE mo_isa_output_nc

  USE mo_logging
  USE mo_kind,                  ONLY: wp, i4

  USE mo_grid_structures,       ONLY: target_grid_def

  USE mo_io_utilities,          ONLY: netcdf_attributes, &
       &                              netcdf_put_var,&
       &                              open_new_netcdf_file, &
       &                              close_netcdf_file, &
       &                              netcdf_get_var, &
       &                              dim_meta_info

  USE mo_var_meta_data,        ONLY: dim_3d_tg, &
      &                              def_com_target_fields_meta,&
      &                              def_isa_fields_meta, &
      &                              isa_tot_npixel_meta,isa_field_meta, &
      &                              def_dimension_info_buffer, &
      &                              lon_geo_meta, &
      &                              lat_geo_meta, &
      &                              isa_field_meta, &
      &                              def_com_target_fields_meta  
  

  IMPLICIT NONE

  PRIVATE

  
  PUBLIC :: write_netcdf_buffer_isa, &
       &    read_netcdf_buffer_isa

  CONTAINS

  !> netcdf output of isa buffer fields
  SUBROUTINE write_netcdf_buffer_isa(netcdf_filename,  &
       &                             tg,         &
       &                             undefined, &
       &                             undef_int,   &
       &                             lon_geo,     &
       &                             lat_geo, &
       &                             isa_tot_npixel, &
       &                             isa_field )


    CHARACTER (len=*), INTENT(IN)      :: netcdf_filename !< filename for the netcdf file
    TYPE(target_grid_def), INTENT(IN)  :: tg !< structure with target grid description

    REAL(KIND=wp), INTENT(IN)          :: undefined, &       !< value to indicate undefined grid elements 
         &                                lon_geo(:,:,:), &  !< longitude coordinates of the target grid in the geographical system
         &                                lat_geo(:,:,:), &  !< latitude coordinates of the target grid in the geographical system
         &                                isa_field(:,:,:)   !< urban fraction due to isa data

    INTEGER (KIND=i4), INTENT(IN)      :: isa_tot_npixel(:,:,:), &  !< total number of isa raw data pixels on target grid
         &                                undef_int       !< value to indicate undefined grid elements

    ! local variables
    INTEGER(KIND=i4)                   :: ndims, ncid, undefined_i, errorcode

    INTEGER(KIND=i4), PARAMETER        :: nglob_atts=6

    TYPE(dim_meta_info), ALLOCATABLE   :: dim_list(:) !< dimensions for netcdf file

    TYPE(netcdf_attributes)            :: global_attributes(nglob_atts)

    CALL logging%info('Enter routine: write_netcdf_buffer_isa')

    !-------------------------------------------------------------
    ! define global attributes
    !CALL set_global_att_lu(i_landuse_data,ilookup_table_lu,global_attributes)
    CALL set_global_att_isa(global_attributes)

    !set up dimensions for buffer
    CALL  def_dimension_info_buffer(tg)
    ! dim_3d_tg

    ! define meta information for various land use related variables for netcdf output
    CALL def_isa_fields_meta(dim_3d_tg)

    ! define meta information for target field variables lon_geo, lat_geo 
    CALL def_com_target_fields_meta(dim_3d_tg)
    ! lon_geo_meta and lat_geo_meta


    !set up dimensions for buffer netcdf output 
    ndims = 3
    ALLOCATE(dim_list(1:ndims),STAT=errorcode)
    IF (errorcode /= 0 ) CALL logging%error('Cant allocate array dim_list',__FILE__,__LINE__)

    dim_list(1) = dim_3d_tg(1) ! ie
    dim_list(2) = dim_3d_tg(2) ! je
    dim_list(3) = dim_3d_tg(3) ! ke

    undefined_i = undef_int

    CALL open_new_netcdf_file(netcdf_filename=TRIM(netcdf_filename),   &
        &                       dim_list=dim_list,                  &
        &                       global_attributes=global_attributes, &
        &                       ncid=ncid)

    ! lon
    CALL netcdf_put_var(ncid,lon_geo,lon_geo_meta,undefined)

    ! lat
    CALL netcdf_put_var(ncid,lat_geo,lat_geo_meta,undefined)

    ! isa_field
    CALL netcdf_put_var(ncid,isa_field,isa_field_meta,undefined)

    ! isa_tot_npixel
    CALL netcdf_put_var(ncid,isa_tot_npixel,isa_tot_npixel_meta,undefined_i)

    CALL close_netcdf_file(ncid)

    CALL logging%info('Exit routine: write_netcdf_buffer_isa')

  END SUBROUTINE write_netcdf_buffer_isa
                                                                          
  !> set global attributes for netcdf with is data
  SUBROUTINE set_global_att_isa(global_attributes)

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
    global_attributes(1)%attributetext='Impervious Surface Area'
    global_attributes(2)%attname = 'institution'
    global_attributes(2)%attributetext='KU Leuven'

    global_attributes(3)%attname = 'source'
    global_attributes(3)%attributetext='NOAA'

    CALL DATE_AND_TIME(ydate,ytime)
    READ(ydate,'(4A2)') cc,yy,mm,dd
    READ(ytime,'(2A2)') hh, minute

    ydate=TRIM(cc)//TRIM(yy)//'-'//TRIM(mm)//'-'//TRIM(dd)
    ytime=TRIM(hh)//':'//TRIM(minute) 

    global_attributes(4)%attname = 'history'
    global_attributes(4)%attributetext=TRIM(ydate)//'T'//TRIM(ytime)//' extpar_isa_to_buffer'

    global_attributes(5)%attname = 'references'
    global_attributes(5)%attributetext='http://ngdc.noaa.gov/eog/dmsp/download_global_isa.html'

    global_attributes(6)%attname = 'comment'
    global_attributes(6)%attributetext='year 2O00 - 2001'

  END SUBROUTINE set_global_att_isa

  SUBROUTINE read_netcdf_buffer_isa(netcdf_filename,  &
       &                            tg,         &
       &                            isa_field,  &
       &                            isa_tot_npixel)


  CHARACTER (len=*), INTENT(IN)     :: netcdf_filename
  TYPE(target_grid_def), INTENT(IN) :: tg
  INTEGER (KIND=i4), INTENT(OUT)    :: isa_tot_npixel(:,:,:)
  REAL (KIND=wp), INTENT(OUT)       :: isa_field(:,:,:)   !< urban fraction due to isa data

  CALL logging%info('Enter routine: read_netcdf_buffer_isa')

  !set up dimensions for buffer
  CALL  def_dimension_info_buffer(tg)
  ! dim_3d_tg

  ! define meta information for various isa related variables for netcdf output
  CALL def_isa_fields_meta(dim_3d_tg)

  ! define meta information for target field variables lon_geo, lat_geo 
  CALL def_com_target_fields_meta(dim_3d_tg)
  ! lon_geo_meta and lat_geo_meta

  CALL netcdf_get_var(TRIM(netcdf_filename),isa_tot_npixel_meta,isa_tot_npixel)

  CALL netcdf_get_var(TRIM(netcdf_filename),isa_field_meta,isa_field)

  CALL logging%info('Exit routine: read_netcdf_buffer_isa')

  END SUBROUTINE read_netcdf_buffer_isa
 
END Module mo_isa_output_nc
