!+ Fortran module with netcdf output routines for GLC2000 data on the target grid
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_0         2010/12/21 Hermann Asensio
!  Initial release
! V1_3         2011/04/19 Hermann Asensio
! introduce Globcover 2009 land use data set for external parameters
! V1_4         2011/04/21 Hermann Asensio
!  clean up
! V2_0         2013/06/07 Martina Messmer
!  Adaptations in a way that the Globcover 2009 data set (6 tiles) can be read
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!> Fortran module with netcdf output routines for GLC2000 data on the target grid
!> ouptut routines
!> \author Hermann Asensio
MODULE mo_landuse_output_nc


  !> kind parameters are defined in MODULE data_parameters
  USE mo_kind, ONLY: wp, i4, i8

  !> data type structures form module GRID_structures
  USE mo_grid_structures, ONLY: igrid_icon,                  &
                                reg_lonlat_grid,             &
                                rotated_lonlat_grid,         &
                                icosahedral_triangular_grid, &
                                target_grid_def

  USE mo_io_utilities, ONLY: var_meta_info
  USE mo_io_utilities, ONLY: netcdf_attributes

  USE mo_io_utilities, ONLY: dim_meta_info

  USE mo_io_utilities, ONLY: netcdf_put_var
  USE mo_io_utilities, ONLY: open_new_netcdf_file
  USE mo_io_utilities, ONLY: close_netcdf_file
  USE mo_io_utilities, ONLY: netcdf_def_grid_mapping

  USE mo_io_utilities, ONLY: get_date_const_field
  USE mo_io_utilities, ONLY: set_date_mm_extpar_field

  USE mo_io_utilities, ONLY: vartype_int 
  USE mo_io_utilities, ONLY: vartype_real
  USE mo_io_utilities, ONLY: vartype_char

  !> abort_extpar defined in MODULE utilities_extpar
  USE mo_utilities_extpar, ONLY: abort_extpar

  USE mo_glc2000_lookup_tables, ONLY: nclass_glc2000

  USE mo_glcc_lookup_tables, ONLY: nclass_glcc


  IMPLICIT NONE

  PRIVATE

  PUBLIC :: write_netcdf_buffer_glc2000
  PUBLIC :: write_netcdf_cosmo_grid_glc2000
  PUBLIC :: write_netcdf_icon_grid_glc2000
  PUBLIC :: read_netcdf_buffer_glc2000

  
  PUBLIC :: write_netcdf_buffer_glcc
  PUBLIC :: write_netcdf_cosmo_grid_glcc
  PUBLIC :: write_netcdf_icon_grid_glcc
  PUBLIC :: read_netcdf_buffer_glcc

  PUBLIC :: write_netcdf_buffer_lu
  PUBLIC :: read_netcdf_buffer_lu

  PUBLIC :: write_netcdf_buffer_ecoclimap
  PUBLIC :: read_netcdf_buffer_ecoclimap

  CONTAINS

  !> netcdf output of landuse buffer fields
  SUBROUTINE write_netcdf_buffer_lu(netcdf_filename,  &
    &                                     lu_dataset, &
    &                                     tg,         &
    &                                     i_landuse_data, &
    &                                     ilookup_table_lu, &
    &                                     nclass_lu, &
    &                                     undefined, &
    &                                     undef_int,   &
    &                                     lon_geo,     &
    &                                     lat_geo, &
    &                                     fr_land_lu, &
    &                                     lu_class_fraction,    &
    &                                     lu_class_npixel, &
    &                                     lu_tot_npixel, &
    &                                     ice_lu, &
    &                                     z0_lu, &
    &                                     root_lu, &
    &                                     plcov_mn_lu, &
    &                                     plcov_mx_lu, &
    &                                     lai_mn_lu, &
    &                                     lai_mx_lu, &
    &                                     rs_min_lu, &
    &                                     urban_lu,  &
    &                                     for_d_lu,  &
    &                                     for_e_lu, &
    &                                     skinc_lu, &
    &                                     emissivity_lu)

  USE mo_var_meta_data, ONLY: dim_3d_tg, &
    &                         def_dimension_info_buffer


  USE mo_var_meta_data, ONLY: lon_geo_meta, &
    &                         lat_geo_meta, &
    &                         no_raw_data_pixel_meta, &
    &                         def_com_target_fields_meta  
  
  USE mo_var_meta_data, ONLY: def_lu_fields_meta

  USE mo_var_meta_data, ONLY: dim_lu_tg

  USE mo_var_meta_data, ONLY: fr_land_lu_meta, lu_tot_npixel_meta, &
    &       lu_class_fraction_meta, lu_class_npixel_meta, &
    &       ice_lu_meta, z0_lu_meta, &
    &       plcov_mx_lu_meta, plcov_mn_lu_meta, &
    &       lai_mx_lu_meta, lai_mn_lu_meta, &
    &       rs_min_lu_meta, urban_lu_meta, &
    &       for_d_lu_meta, for_e_lu_meta, &
    &       skinc_lu_meta, &
    &       emissivity_lu_meta, root_lu_meta

  CHARACTER (len=*), INTENT(IN)      :: netcdf_filename !< filename for the netcdf file
  CHARACTER (len=*), INTENT(IN)      :: lu_dataset !< name of landuse data set
  TYPE(target_grid_def), INTENT(IN) :: tg !< structure with target grid description
  INTEGER, INTENT(IN) :: i_landuse_data !<integer switch to choose a land use raw data set
  INTEGER, INTENT(IN) :: ilookup_table_lu !< integer switch to choose a lookup table
  INTEGER, INTENT(IN) :: nclass_lu !< number of land use classes 

  REAL(KIND=wp), INTENT(IN)          :: undefined       !< value to indicate undefined grid elements 
  INTEGER, INTENT(IN)                :: undef_int       !< value to indicate undefined grid elements
  REAL (KIND=wp), INTENT(IN) :: lon_geo(:,:,:)  !< longitude coordinates of the target grid in the geographical system
  REAL (KIND=wp), INTENT(IN) :: lat_geo(:,:,:)  !< latitude coordinates of the target grid in the geographical system
  REAL (KIND=wp), INTENT(IN)  :: lu_class_fraction(:,:,:,:)  
                                 !< fraction for each lu class on target grid (dim (ie,je,ke,nclass_lu))
  INTEGER (KIND=i8), INTENT(IN) :: lu_class_npixel(:,:,:,:) 
                                   !< number of raw data pixels for each lu class on target grid (dim (ie,je,ke,nclass_lu))
  INTEGER (KIND=i8), INTENT(IN) :: lu_tot_npixel(:,:,:)  
                                   !< total number of lu raw data pixels on target grid (dimension (ie,je,ke))
  REAL (KIND=wp), INTENT(IN)  :: fr_land_lu(:,:,:) !< fraction land due to lu raw data
  REAL (KIND=wp), INTENT(IN)  :: ice_lu(:,:,:)     !< fraction of ice due to lu raw data
  REAL (KIND=wp), INTENT(IN)  :: z0_lu(:,:,:)      !< roughness length due to lu land use data
  REAL (KIND=wp), INTENT(IN)  :: root_lu(:,:,:)    !< root depth due to lu land use data
  REAL (KIND=wp), INTENT(IN)  :: plcov_mx_lu(:,:,:)!< plant cover maximum due to lu land use data
  REAL (KIND=wp), INTENT(IN)  :: plcov_mn_lu(:,:,:)!< plant cover minimum due to lu land use data
  REAL (KIND=wp), INTENT(IN)  :: lai_mx_lu(:,:,:)  !< Leaf Area Index maximum due to lu land use data
  REAL (KIND=wp), INTENT(IN)  :: lai_mn_lu(:,:,:)  !< Leaf Area Index minimum due to lu land use data
  REAL (KIND=wp), INTENT(IN)  :: rs_min_lu(:,:,:)  !< minimal stomata resistance due to lu land use data
  REAL (KIND=wp), INTENT(IN)  :: urban_lu(:,:,:)   !< urban fraction due to lu land use data
  REAL (KIND=wp), INTENT(IN)  :: for_d_lu(:,:,:)   !< deciduous forest (fraction) due to lu land use data
  REAL (KIND=wp), INTENT(IN)  :: for_e_lu(:,:,:)   !< evergreen forest (fraction) due to lu land use data
  REAL (KIND=wp), INTENT(IN)  :: skinc_lu(:,:,:)   !< skin conductivity due to lu land use data
  REAL (KIND=wp), INTENT(IN)  :: emissivity_lu(:,:,:) !< longwave emissivity due to lu land use data


  ! local variables
  INTEGER :: ndims  
  INTEGER :: ncid
  INTEGER (KIND=i8) :: undefined_i

  TYPE(dim_meta_info), ALLOCATABLE :: dim_list(:) !< dimensions for netcdf file

  INTEGER, PARAMETER :: nglob_atts=6
  TYPE(netcdf_attributes) :: global_attributes(nglob_atts)

  INTEGER :: errorcode !< error status variable

  INTEGER :: n !< counter

  PRINT *,'ENTER write_netcdf_buffer_lu'

  !-------------------------------------------------------------
  ! define global attributes
  CALL set_global_att_lu(i_landuse_data,ilookup_table_lu,global_attributes)



  !set up dimensions for buffer
  CALL  def_dimension_info_buffer(tg)
  ! dim_3d_tg



  ! define meta information for various land use related variables for netcdf output
  CALL def_lu_fields_meta(tg,nclass_lu,dim_3d_tg,lu_dataset=lu_dataset)
  ! dim_lu_tg
  ! fr_land_lu_meta, lu_tot_npixel_meta, &
  !  &       lu_class_fraction_meta, lu_class_npixel_meta, &
  !  &       ice_lu_meta, z0_lu_meta, &
  !  &       plcov_mx_lu_meta, plcov_mn_lu_meta, &
  !  &       lai_mx_lu_meta, lai_mn_lu_meta, &
  !  &       rs_min_lu_meta, urban_lu_meta, &
  !  &       for_d_lu_meta, for_e_lu_meta, &
  !  &       emissivity_lu_meta, root_lu_meta



  ! define meta information for target field variables lon_geo, lat_geo 
  CALL def_com_target_fields_meta(dim_3d_tg)
  ! lon_geo_meta and lat_geo_meta


  !set up dimensions for buffer netcdf output 
  ndims = 4
  ALLOCATE(dim_list(1:ndims),STAT=errorcode)
  IF (errorcode /= 0 ) CALL abort_extpar('Cant allocate array dim_list')
  dim_list = dim_lu_tg

  !  dim_3d_buffer(:) = dim_list(:)

  undefined_i = undef_int

  !-----------------------------------------------------------------


  CALL open_new_netcdf_file(netcdf_filename=TRIM(netcdf_filename),   &
      &                       dim_list=dim_list,                  &
      &                       global_attributes=global_attributes, &
      &                       ncid=ncid)

  ! lon
  CALL netcdf_put_var(ncid,lon_geo,lon_geo_meta,undefined)

  ! lat
  CALL netcdf_put_var(ncid,lat_geo,lat_geo_meta,undefined)

  ! fr_land_lu
  CALL netcdf_put_var(ncid,fr_land_lu,fr_land_lu_meta,undefined)

  ! ice_lu
  PRINT *,'MAX ICE_LU BUFFER: ', MAXVAL(ICE_LU)
  CALL netcdf_put_var(ncid,ice_lu,ice_lu_meta,undefined)

  ! plcov_mx_lu
  CALL netcdf_put_var(ncid,plcov_mx_lu,plcov_mx_lu_meta,undefined)

  ! lai_mx_lu
  CALL netcdf_put_var(ncid,lai_mx_lu,lai_mx_lu_meta,undefined)

  ! rs_min_lu
  CALL netcdf_put_var(ncid,rs_min_lu,rs_min_lu_meta,undefined)

  ! urban_lu
  CALL netcdf_put_var(ncid,urban_lu,urban_lu_meta,undefined)

  ! for_d_lu
  CALL netcdf_put_var(ncid,for_d_lu,for_d_lu_meta,undefined)

  ! for_e_lu
  CALL netcdf_put_var(ncid,for_e_lu,for_e_lu_meta,undefined)

! skinc_lu
  CALL netcdf_put_var(ncid,skinc_lu,skinc_lu_meta,undefined)

  ! emissivity_lu
  CALL netcdf_put_var(ncid,emissivity_lu,emissivity_lu_meta,undefined)

  ! root_lu
  CALL netcdf_put_var(ncid,root_lu,root_lu_meta,undefined)

  ! z0_lu
  CALL netcdf_put_var(ncid,z0_lu,z0_lu_meta,undefined)

  ! lai_mn_lu
  CALL netcdf_put_var(ncid,lai_mn_lu,lai_mn_lu_meta,undefined)
    
  ! plcov_mn_lu
  CALL netcdf_put_var(ncid,plcov_mn_lu,plcov_mn_lu_meta,undefined)

  ! lu_tot_npixel
  CALL netcdf_put_var(ncid,lu_tot_npixel,lu_tot_npixel_meta,undefined_i)

  ! lu_class_fraction
  CALL netcdf_put_var(ncid,lu_class_fraction,lu_class_fraction_meta,undefined)
  !-----------------------------------------------------------------

  ! lu_class_npixel
  CALL netcdf_put_var(ncid,lu_class_npixel,lu_class_npixel_meta,undefined_i)
  !-----------------------------------------------------------------

  CALL close_netcdf_file(ncid)

END SUBROUTINE write_netcdf_buffer_lu

!----------------------------------------------------------------------------------
  !> netcdf output of landuse ecoclimap buffer fields
  SUBROUTINE write_netcdf_buffer_ecoclimap(netcdf_filename,  &
    &                                     tg,         &
    &                                     i_landuse_data, &
    &                                     ilookup_table_lu, &
    &                                     nclass_lu, &
    &                                     undefined, &
    &                                     undef_int,   &
    &                                     lon_geo,     &
    &                                     lat_geo, &
    &                                     fr_land_lu, &
    &                                     lu_class_fraction,    &
    &                                     lu_class_npixel, &
    &                                     lu_tot_npixel, &
    &                                     ice_lu, &
    &                                     z012_lu, &
    &                                     root_lu, &
    &                                     plcov12_lu, &
    &                                     lai12_lu, &
    &                                     rs_min_lu, &
    &                                     urban_lu,  &
    &                                     for_d_lu,  &
    &                                     for_e_lu, &
    &                                     emissivity_lu)

  USE mo_var_meta_data, ONLY: dim_3d_tg, &
    &                         def_dimension_info_buffer


  USE mo_var_meta_data, ONLY: lon_geo_meta, &
    &                         lat_geo_meta, &
    &                         no_raw_data_pixel_meta, &
    &                         def_com_target_fields_meta  
  
  USE mo_var_meta_data, ONLY: def_ecoclimap_fields_meta

  USE mo_var_meta_data, ONLY: dim_ecoclimap_tg,   dim_ecoclimap_tg2 
    

  USE mo_var_meta_data, ONLY: fr_land_lu_meta, lu_tot_npixel_meta, &
    &       lu_class_fraction_meta, lu_class_npixel_meta, &
    &       ice_lu_meta, z012_lu_meta,  &
    &       plcov12_lu_meta, &
    &       lai12_lu_meta, &
    &       rs_min_lu_meta, urban_lu_meta, &
    &       for_d_lu_meta, for_e_lu_meta, &
    &       emissivity_lu_meta, root_lu_meta

  CHARACTER (len=*), INTENT(IN)      :: netcdf_filename !< filename for the netcdf file
  TYPE(target_grid_def), INTENT(IN) :: tg !< structure with target grid description
  INTEGER, INTENT(IN) :: i_landuse_data !<integer switch to choose a land use raw data set
  INTEGER (KIND=i4) :: ntime !< number of times
  INTEGER, INTENT(IN) :: ilookup_table_lu !< integer switch to choose a lookup table
  INTEGER, INTENT(IN) :: nclass_lu !< number of land use classes 


  REAL(KIND=wp), INTENT(IN)          :: undefined       !< value to indicate undefined grid elements 
  INTEGER, INTENT(IN)                :: undef_int       !< value to indicate undefined grid elements
  REAL (KIND=wp), INTENT(IN) :: lon_geo(:,:,:)  !< longitude coordinates of the target grid in the geographical system
  REAL (KIND=wp), INTENT(IN) :: lat_geo(:,:,:)  !< latitude coordinates of the target grid in the geographical system
  REAL (KIND=wp), INTENT(IN)  :: lu_class_fraction(:,:,:,:)  
                                 !< fraction for each lu class on target grid (dim (ie,je,ke,nclass_lu))
  INTEGER (KIND=i8), INTENT(IN) :: lu_class_npixel(:,:,:,:) 
                                 !< number of raw data pixels for each lu class on target grid (dim (ie,je,ke,nclass_lu))
  INTEGER (KIND=i8), INTENT(IN) :: lu_tot_npixel(:,:,:)  
                                   !< total number of lu raw data pixels on target grid (dimension (ie,je,ke))
  REAL (KIND=wp), INTENT(IN)  :: fr_land_lu(:,:,:) !< fraction land due to lu raw data
  REAL (KIND=wp), INTENT(IN)  :: ice_lu(:,:,:)     !< fraction of ice due to lu raw data
  REAL (KIND=wp), INTENT(IN)  :: z012_lu(:,:,:,:)      !< roughness length due to lu land use data
  REAL (KIND=wp), INTENT(IN)  :: root_lu(:,:,:)    !< root depth due to lu land use data
  REAL (KIND=wp), INTENT(IN)  :: plcov12_lu(:,:,:,:)!< plant cover minimum due to lu land use data
  REAL (KIND=wp), INTENT(IN)  :: lai12_lu(:,:,:,:)  !< Leaf Area Index minimum due to lu land use data
  REAL (KIND=wp), INTENT(IN)  :: rs_min_lu(:,:,:)  !< minimal stomata resistance due to lu land use data
  REAL (KIND=wp), INTENT(IN)  :: urban_lu(:,:,:)   !< urban fraction due to lu land use data
  REAL (KIND=wp), INTENT(IN)  :: for_d_lu(:,:,:)   !< deciduous forest (fraction) due to lu land use data
  REAL (KIND=wp), INTENT(IN)  :: for_e_lu(:,:,:)   !< evergreen forest (fraction) due to lu land use data
  REAL (KIND=wp), INTENT(IN)  :: emissivity_lu(:,:,:) !< longwave emissivity due to lu land use data


  ! local variables
  INTEGER :: ndims  
  INTEGER :: ncid
  INTEGER (KIND=i8) :: undefined_i
!gs_24.04.12  
  REAL (KIND=wp),ALLOCATABLE :: time(:) !< time variable
  INTEGER (KIND=i8) :: dataDate  !< date, for edition independent use of GRIB_API dataDate as Integer in the format ccyymmdd
   INTEGER (KIND=i8) :: dataTime  !< time, for edition independent use GRIB_API dataTime in the format hhmm
!>
  TYPE(dim_meta_info), ALLOCATABLE :: dim_list(:) !< dimensions for netcdf file

  INTEGER, PARAMETER :: nglob_atts=6
  TYPE(netcdf_attributes) :: global_attributes(nglob_atts)

  INTEGER :: errorcode !< error status variable

  INTEGER :: n !< counter

  ntime = 12

  PRINT *,'ECOCLIMAP: ENTER write_netcdf_buffer_lu'

  !-------------------------------------------------------------
  ! define global attributes
  CALL set_global_att_lu(i_landuse_data,ilookup_table_lu,global_attributes)

  !set up dimensions for buffer
  CALL  def_dimension_info_buffer(tg)
  ! dim_3d_tg
  ! define meta information for various land use related variables for netcdf output
  CALL def_ecoclimap_fields_meta(tg,ntime,nclass_lu,dim_3d_tg)
  
  ! define meta information for target field variables lon_geo, lat_geo 
  CALL def_com_target_fields_meta(dim_3d_tg)
  ! lon_geo_meta and lat_geo_meta

!gs_24.04.12
    ALLOCATE(time(1:ntime),STAT=errorcode)
    IF (errorcode /= 0 ) CALL abort_extpar('Cant allocate array time')
    DO n=1,ntime
      CALL set_date_mm_extpar_field(n,dataDate,dataTime)
!      time(n) = REAL(dataDate,wp) + REAL(dataTime,wp)/10000. ! units = "day as %Y%m%d.%f"
       time(n) = REAL(n,wp) ! months !_gs 20.07.12
    ENDDO
!>

  !set up dimensions for buffer netcdf output 
  ndims = 5 
  ALLOCATE(dim_list(1:ndims),STAT=errorcode)
  IF (errorcode /= 0 ) CALL abort_extpar('Cant allocate array dim_list')
  dim_list = dim_ecoclimap_tg

  !  dim_3d_buffer(:) = dim_list(:)

  undefined_i = undef_int
  !-----------------------------------------------------------------

  PRINT *,'ECOCLIMAP: write ', TRIM(netcdf_filename) 

  CALL open_new_netcdf_file(netcdf_filename=TRIM(netcdf_filename),   &
      &                       dim_list=dim_list,                     &
      &                       global_attributes=global_attributes,   &
      &                       time=time,                             &
      &                       ncid=ncid)

  ! lon
  CALL netcdf_put_var(ncid,lon_geo,lon_geo_meta,undefined)

  ! lat
  CALL netcdf_put_var(ncid,lat_geo,lat_geo_meta,undefined)

  ! fr_land_lu
  CALL netcdf_put_var(ncid,fr_land_lu,fr_land_lu_meta,undefined)

  ! ice_lu
  CALL netcdf_put_var(ncid,ice_lu,ice_lu_meta,undefined)

  ! rs_min_lu
  CALL netcdf_put_var(ncid,rs_min_lu,rs_min_lu_meta,undefined)

  ! urban_lu
  CALL netcdf_put_var(ncid,urban_lu,urban_lu_meta,undefined)

  ! for_d_lu
  CALL netcdf_put_var(ncid,for_d_lu,for_d_lu_meta,undefined)

  ! for_e_lu
  CALL netcdf_put_var(ncid,for_e_lu,for_e_lu_meta,undefined)

  ! root_lu
  CALL netcdf_put_var(ncid,root_lu,root_lu_meta,undefined)

  ! z012_lu
  CALL netcdf_put_var(ncid,z012_lu,z012_lu_meta,undefined)

  ! plcov12_lu
  CALL netcdf_put_var(ncid,plcov12_lu,plcov12_lu_meta,undefined)

  ! lai12_lu
  CALL netcdf_put_var(ncid,lai12_lu,lai12_lu_meta,undefined)
    
  ! lu_tot_npixel
  CALL netcdf_put_var(ncid,lu_tot_npixel,lu_tot_npixel_meta,undefined_i)

  ! lu_class_fraction
  CALL netcdf_put_var(ncid,lu_class_fraction,lu_class_fraction_meta,undefined)
  !-----------------------------------------------------------------

  ! emissivity_lu
  CALL netcdf_put_var(ncid,emissivity_lu,emissivity_lu_meta,undefined)

  ! lu_class_npixel
  CALL netcdf_put_var(ncid,lu_class_npixel,lu_class_npixel_meta,undefined_i)
  !-----------------------------------------------------------------

  CALL close_netcdf_file(ncid)

  PRINT *,'ECOCLIMAP: wrote all data fields'

END SUBROUTINE write_netcdf_buffer_ecoclimap

 !-----------------------------------------------------------------------
  !> set global attributes for netcdf with lu data
  SUBROUTINE set_global_att_lu(i_landuse_data,ilookup_table_lu,global_attributes)
    USE mo_io_units, ONLY: filename_max
    USE mo_lu_tg_fields, ONLY :  i_lu_globcover, i_lu_glc2000, i_lu_glcc
    USE mo_lu_tg_fields, ONLY :  i_lu_ecoclimap
    USE mo_globcover_lookup_tables, ONLY: get_name_globcover_lookup_tables
    USE mo_ecoclimap_lookup_tables, ONLY: get_name_ecoclimap_lookup_tables
    USE mo_glc2000_lookup_tables, ONLY: get_name_glc2000_lookup_tables
    USE mo_glcc_lookup_tables, ONLY: get_name_glcc_lookup_tables
    INTEGER, INTENT(IN) :: i_landuse_data !<integer switch to choose a land use raw data set
    INTEGER, INTENT(IN) :: ilookup_table_lu !< integer switch to choose a lookup table
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

    CHARACTER(len=filename_max) :: name_lookup_table_lu



    ! define global attributes
    
    global_attributes(1)%attname = 'title'
    global_attributes(1)%attributetext='Land Use data'
    global_attributes(2)%attname = 'institution'
    global_attributes(2)%attributetext='Deutscher Wetterdienst'

    global_attributes(3)%attname = 'source'
    SELECT CASE (i_landuse_data)
      CASE (i_lu_ecoclimap)
         global_attributes(2)%attributetext='KIT/IMK-IFU'
         global_attributes(3)%attributetext='ECOCLIMAP_V2  data'
         CALL get_name_ecoclimap_lookup_tables(ilookup_table_lu, name_lookup_table_lu)
      CASE (i_lu_globcover)
        global_attributes(3)%attributetext='Globcover 2009 data'
         CALL get_name_globcover_lookup_tables(ilookup_table_lu, name_lookup_table_lu)
      CASE (i_lu_glc2000)
        global_attributes(3)%attributetext='GLC2000 data'
         CALL get_name_glc2000_lookup_tables(ilookup_table_lu, name_lookup_table_lu)
      CASE(i_lu_glcc)
        global_attributes(3)%attributetext='GLCC data'
         CALL get_name_glcc_lookup_tables(ilookup_table_lu, name_lookup_table_lu)
    END SELECT

    

    CALL DATE_AND_TIME(ydate,ytime)
    READ(ydate,'(4A2)') cc,yy,mm,dd
    READ(ytime,'(2A2)') hh, minute

    ydate=TRIM(cc)//TRIM(yy)//'-'//TRIM(mm)//'-'//TRIM(dd)
    ytime=TRIM(hh)//':'//TRIM(minute) 

    global_attributes(4)%attname = 'history'
    global_attributes(4)%attributetext=TRIM(ydate)//'T'//TRIM(ytime)//' extpar_landuse_to_buffer'

    global_attributes(5)%attname = 'references'
    global_attributes(5)%attributetext=''

    global_attributes(6)%attname = 'comment'
    global_attributes(6)%attributetext='Landuse data look-up table: '//TRIM(name_lookup_table_lu)

  END SUBROUTINE set_global_att_lu
  !-----------------------------------------------------------------------


  !> read land use derived buffer fields
  SUBROUTINE read_netcdf_buffer_lu(netcdf_filename,  &
    &                                     tg,         &
    &                                     nclass_lu, &
    &                                     undefined, &
    &                                     undef_int,   &
    &                                     fr_land_lu, &
    &                                     lu_class_fraction,    &
    &                                     lu_class_npixel, &
    &                                     lu_tot_npixel, &
    &                                     ice_lu, &
    &                                     z0_lu, &
    &                                     root_lu, &
    &                                     plcov_mn_lu, &
    &                                     plcov_mx_lu, &
    &                                     lai_mn_lu, &
    &                                     lai_mx_lu, &
    &                                     rs_min_lu, &
    &                                     urban_lu,  &
    &                                     for_d_lu,  &
    &                                     for_e_lu, &
    &                                     skinc_lu, &
    &                                     emissivity_lu)


  USE mo_var_meta_data, ONLY: dim_3d_tg, &
    &                         def_dimension_info_buffer


  USE mo_var_meta_data, ONLY: lon_geo_meta, &
    &                         lat_geo_meta, &
    &                         no_raw_data_pixel_meta, &
    &                         def_com_target_fields_meta  
  
  USE mo_var_meta_data, ONLY: def_lu_fields_meta

  USE mo_var_meta_data, ONLY: dim_lu_tg

  USE mo_var_meta_data, ONLY: fr_land_lu_meta, lu_tot_npixel_meta, &
    &       lu_class_fraction_meta, lu_class_npixel_meta, &
    &       ice_lu_meta, z0_lu_meta, &
    &       plcov_mx_lu_meta, plcov_mn_lu_meta, &
    &       lai_mx_lu_meta, lai_mn_lu_meta, &
    &       rs_min_lu_meta, urban_lu_meta, &
    &       for_d_lu_meta, for_e_lu_meta, &
    &       skinc_lu_meta,                &
    &       emissivity_lu_meta, root_lu_meta

  USE mo_io_utilities, ONLY: netcdf_get_var

  CHARACTER (len=*), INTENT(IN)      :: netcdf_filename !< filename for the netcdf file
  TYPE(target_grid_def), INTENT(IN) :: tg !< structure with target grid description
  INTEGER, INTENT(IN) :: nclass_lu !< number of land use classes 
  REAL(KIND=wp), INTENT(OUT)          :: undefined       !< value to indicate undefined grid elements 
  INTEGER, INTENT(OUT)                :: undef_int       !< value to indicate undefined grid elements
  REAL (KIND=wp), INTENT(OUT)  :: lu_class_fraction(:,:,:,:)  
                                  !< fraction for each lu class on target grid (dim (ie,je,ke,nclass_lu))
  INTEGER (KIND=i8), INTENT(OUT) :: lu_class_npixel(:,:,:,:) 
                                    !< number of raw data pixels for each lu class on target grid (dim (ie,je,ke,nclass_lu))
  INTEGER (KIND=i8), INTENT(OUT) :: lu_tot_npixel(:,:,:)  
                                    !< total number of lu raw data pixels on target grid (dimension (ie,je,ke))
  REAL (KIND=wp), INTENT(OUT)  :: fr_land_lu(:,:,:) !< fraction land due to lu raw data
  REAL (KIND=wp), INTENT(OUT)  :: ice_lu(:,:,:)     !< fraction of ice due to lu raw data
  REAL (KIND=wp), INTENT(OUT)  :: z0_lu(:,:,:)      !< roughness length due to lu land use data
  REAL (KIND=wp), INTENT(OUT)  :: root_lu(:,:,:)    !< root depth due to lu land use data
  REAL (KIND=wp), INTENT(OUT)  :: plcov_mn_lu(:,:,:)!< plant cover maximum due to lu land use data
  REAL (KIND=wp), INTENT(OUT)  :: plcov_mx_lu(:,:,:)!< plant cover minimum due to lu land use data
  REAL (KIND=wp), INTENT(OUT)  :: lai_mn_lu(:,:,:)  !< Leaf Area Index maximum due to lu land use data
  REAL (KIND=wp), INTENT(OUT)  :: lai_mx_lu(:,:,:)  !< Leaf Area Index minimum due to lu land use data
  REAL (KIND=wp), INTENT(OUT)  :: rs_min_lu(:,:,:)  !< minimal stomata resistance due to lu land use data
  REAL (KIND=wp), INTENT(OUT)  :: urban_lu(:,:,:)   !< urban fraction due to lu land use data
  REAL (KIND=wp), INTENT(OUT)  :: for_d_lu(:,:,:)   !< deciduous forest (fraction) due to lu land use data
  REAL (KIND=wp), INTENT(OUT)  :: for_e_lu(:,:,:)   !< evergreen forest (fraction) due to lu land use data
  REAL (KIND=wp), INTENT(OUT)  :: skinc_lu(:,:,:)   !< skin conductivity due to lu land use data
  REAL (KIND=wp), INTENT(OUT)  :: emissivity_lu(:,:,:) !< longwave emissivity due to lu land use data


  ! local variables
  INTEGER :: errorcode !< error status variable
  INTEGER :: n !< counter

  PRINT *,'ENTER read_netcdf_buffer_lu'

  !set up dimensions for buffer
  CALL  def_dimension_info_buffer(tg)
  ! dim_3d_tg

  ! define meta information for various land use related variables for netcdf output
  CALL def_lu_fields_meta(tg,nclass_lu,dim_3d_tg)
  ! dim_lu_tg
  ! fr_land_lu_meta, lu_tot_npixel_meta, &
  !  &       lu_class_fraction_meta, lu_class_npixel_meta, &
  !  &       ice_lu_meta, z0_lu_meta, &
  !  &       plcov_mx_lu_meta, plcov_mn_lu_meta, &
  !  &       lai_mx_lu_meta, lai_mn_lu_meta, &
  !  &       rs_min_lu_meta, urban_lu_meta, &
  !  &       for_d_lu_meta, for_e_lu_meta, &
  !  &       emissivity_lu_meta, root_lu_meta

  ! define meta information for target field variables lon_geo, lat_geo 
  CALL def_com_target_fields_meta(dim_3d_tg)
  ! lon_geo_meta and lat_geo_meta

  PRINT *,'CALL read netcdf data Land Use'

  CALL netcdf_get_var(TRIM(netcdf_filename),fr_land_lu_meta,fr_land_lu)
  PRINT *,'fr_land_lu read'

  CALL netcdf_get_var(TRIM(netcdf_filename),lu_tot_npixel_meta,lu_tot_npixel)
  PRINT *,'lu_tot_npixel read'

  CALL netcdf_get_var(TRIM(netcdf_filename),lu_class_fraction_meta,lu_class_fraction)
  PRINT *,'lu_class_fraction read'

  CALL netcdf_get_var(TRIM(netcdf_filename),lu_class_npixel_meta,lu_class_npixel)
  PRINT *,'lu_class_npixel read'

!  IF (tg%igrid_type /= igrid_icon) THEN
    CALL netcdf_get_var(TRIM(netcdf_filename),ice_lu_meta,ice_lu)
    PRINT *,'ice_lu read - MAX ICE_LU BUFFER: ', MAXVAL(ICE_LU)
!  ENDIF
  
  CALL netcdf_get_var(TRIM(netcdf_filename),z0_lu_meta,z0_lu)
  PRINT *,'z0_lu read'

  CALL netcdf_get_var(TRIM(netcdf_filename),plcov_mx_lu_meta,plcov_mx_lu)
  PRINT *,'plcov_mx_lu read'

  CALL netcdf_get_var(TRIM(netcdf_filename),plcov_mn_lu_meta,plcov_mn_lu)
  PRINT *,'plcov_mn_lu read'

  CALL netcdf_get_var(TRIM(netcdf_filename),lai_mx_lu_meta,lai_mx_lu)
  PRINT *,'lai_mx_lu read'

  CALL netcdf_get_var(TRIM(netcdf_filename),lai_mn_lu_meta,lai_mn_lu)
  PRINT *,'lai_mn_lu read'

  CALL netcdf_get_var(TRIM(netcdf_filename),rs_min_lu_meta,rs_min_lu)
  PRINT *,'rs_min_lu read'

  CALL netcdf_get_var(TRIM(netcdf_filename),urban_lu_meta,urban_lu)
  PRINT *,'urban_lu read'

  CALL netcdf_get_var(TRIM(netcdf_filename),for_d_lu_meta,for_d_lu)
  PRINT *,'for_d_lu read'

  CALL netcdf_get_var(TRIM(netcdf_filename),for_e_lu_meta,for_e_lu)
  PRINT *,'for_e_lu read'

  CALL netcdf_get_var(TRIM(netcdf_filename),skinc_lu_meta,skinc_lu)
  PRINT *,'skinc_lu read'

  CALL netcdf_get_var(TRIM(netcdf_filename),emissivity_lu_meta,emissivity_lu)
  PRINT *,'emissivity_lu read'

  CALL netcdf_get_var(TRIM(netcdf_filename),root_lu_meta,root_lu)
  PRINT *,'root_lu read'


END SUBROUTINE read_netcdf_buffer_lu
!-----------------------------------------------------------------------


!--------------------------------------------------------------------------
! read ECOCLIMAP buffer netcdf file
!--------------------------------------------------------------------------
  SUBROUTINE read_netcdf_buffer_ecoclimap(netcdf_filename,  &
    &                                     tg,         &
    &                                     nclass_lu, &
    &                                     undefined, &
    &                                     undef_int,   &
    &                                     fr_land_lu, &
    &                                     ecoclimap_class_fraction,    &
    &                                     lu_class_npixel, &
    &                                     lu_tot_npixel, &
    &                                     ice_lu, &
    &                                     z012_lu, &
    &                                     root_lu, &
    &                                     plcov12_lu, &
    &                                     lai12_lu, &
    &                                     rs_min_lu, &
    &                                     urban_lu,  &
    &                                     for_d_lu,  &
    &                                     for_e_lu, &
    &                                     emissivity_lu)


  USE mo_var_meta_data, ONLY: dim_3d_tg, &
    &                         def_dimension_info_buffer


  USE mo_var_meta_data, ONLY: lon_geo_meta, &
    &                         lat_geo_meta, &
    &                         no_raw_data_pixel_meta, &
    &                         def_com_target_fields_meta  
  
  !gs_23.04.12
  USE mo_var_meta_data, ONLY: def_ecoclimap_fields_meta
  USE mo_var_meta_data, ONLY: dim_ecoclimap_tg,   dim_ecoclimap_tg2
  !>

!  USE mo_var_meta_data, ONLY: dim_ecoclimap_tg

  USE mo_var_meta_data, ONLY: fr_land_lu_meta, lu_tot_npixel_meta, &
    &      lu_class_fraction_meta, lu_class_npixel_meta, &
    &       ice_lu_meta, z012_lu_meta, &
    &       plcov12_lu_meta,         &
    &       lai12_lu_meta,          &
    &       rs_min_lu_meta, urban_lu_meta, &
    &       for_d_lu_meta, for_e_lu_meta, &
    &       emissivity_lu_meta, root_lu_meta

  USE mo_io_utilities, ONLY: netcdf_get_var

  CHARACTER (len=*), INTENT(IN)      :: netcdf_filename !< filename for the netcdf file
  TYPE(target_grid_def), INTENT(IN) :: tg !< structure with target grid description
  INTEGER, INTENT(IN) :: nclass_lu !< number of land use classes 
  REAL(KIND=wp), INTENT(OUT)          :: undefined       !< value to indicate undefined grid elements 
  INTEGER, INTENT(OUT)                :: undef_int       !< value to indicate undefined grid elements
  REAL (KIND=wp), INTENT(OUT)  :: ecoclimap_class_fraction(:,:,:,:)  
                                  !< fraction for each lu class on target grid (dim (ie,je,ke,nclass_lu))
  INTEGER (KIND=i8), INTENT(OUT) :: lu_class_npixel(:,:,:,:) 
                                    !< number of raw data pixels for each lu class on target grid (dim (ie,je,ke,nclass_lu))
  INTEGER (KIND=i8), INTENT(OUT) :: lu_tot_npixel(:,:,:)  
                                    !< total number of lu raw data pixels on target grid (dimension (ie,je,ke))
  REAL (KIND=wp), INTENT(OUT)  :: fr_land_lu(:,:,:) !< fraction land due to lu raw data
  REAL (KIND=wp), INTENT(OUT)  :: ice_lu(:,:,:)     !< fraction of ice due to lu raw data
  REAL (KIND=wp), INTENT(OUT)  :: z012_lu(:,:,:,:)      !< roughness length due to lu land use data
  REAL (KIND=wp), INTENT(OUT)  :: root_lu(:,:,:)    !< root depth due to lu land use data
  REAL (KIND=wp), INTENT(OUT)  :: plcov12_lu(:,:,:,:)!< plant cover maximum due to lu land use data
  REAL (KIND=wp), INTENT(OUT)  :: lai12_lu(:,:,:,:)  !< Leaf Area Index maximum due to lu land use data
  REAL (KIND=wp), INTENT(OUT)  :: rs_min_lu(:,:,:)  !< minimal stomata resistance due to lu land use data
  REAL (KIND=wp), INTENT(OUT)  :: urban_lu(:,:,:)   !< urban fraction due to lu land use data
  REAL (KIND=wp), INTENT(OUT)  :: for_d_lu(:,:,:)   !< deciduous forest (fraction) due to lu land use data
  REAL (KIND=wp), INTENT(OUT)  :: for_e_lu(:,:,:)   !< evergreen forest (fraction) due to lu land use data
  REAL (KIND=wp), INTENT(OUT)  :: emissivity_lu(:,:,:) !< longwave emissivity due to lu land use data


  ! local variables
  INTEGER :: errorcode !< error status variable
  INTEGER :: n !< counter
  INTEGER (KIND=i4) :: ntime !< number of times
  PRINT *,'ECOCLIMAP read_netcdf_buffer_lu: ', TRIM(netcdf_filename)
  ntime=12
  !set up dimensions for buffer
  CALL  def_dimension_info_buffer(tg)
  ! dim_3d_tg

  ! define meta information for various land use related variables  for netcdf output
  CALL def_ecoclimap_fields_meta(tg,ntime,nclass_lu,dim_3d_tg)
  ! define meta information for target field variables lon_geo, lat_geo 
  CALL def_com_target_fields_meta(dim_3d_tg)
  ! lon_geo_meta and lat_geo_meta

  PRINT *,'ECOCLIMAP read netcdf data'

  CALL netcdf_get_var(TRIM(netcdf_filename),fr_land_lu_meta,fr_land_lu)
  PRINT *,'fr_land_lu read'

  CALL netcdf_get_var(TRIM(netcdf_filename),lu_tot_npixel_meta,lu_tot_npixel)
  PRINT *,'lu_tot_npixel read'

  CALL netcdf_get_var(TRIM(netcdf_filename),lu_class_fraction_meta,ecoclimap_class_fraction)
  PRINT *,'lu_class_fraction read'

  CALL netcdf_get_var(TRIM(netcdf_filename),lu_class_npixel_meta,lu_class_npixel)
  PRINT *,'lu_class_npixel read'

  IF (tg%igrid_type /= igrid_icon) THEN
    CALL netcdf_get_var(TRIM(netcdf_filename),ice_lu_meta,ice_lu)
    PRINT *,'ice_lu read'
  ENDIF

  CALL netcdf_get_var(TRIM(netcdf_filename),z012_lu_meta,z012_lu)
  PRINT *,'z0_lu read'

  CALL netcdf_get_var(TRIM(netcdf_filename),plcov12_lu_meta,plcov12_lu)
  PRINT *,'plcov12_lu read'

  CALL netcdf_get_var(TRIM(netcdf_filename),lai12_lu_meta,lai12_lu)
  PRINT *,'lai12_lu read'

  CALL netcdf_get_var(TRIM(netcdf_filename),rs_min_lu_meta,rs_min_lu)
  PRINT *,'rs_min_lu read'

  CALL netcdf_get_var(TRIM(netcdf_filename),urban_lu_meta,urban_lu)
  PRINT *,'urban_lu read'

  CALL netcdf_get_var(TRIM(netcdf_filename),for_d_lu_meta,for_d_lu)
  PRINT *,'for_d_lu read'

  CALL netcdf_get_var(TRIM(netcdf_filename),for_e_lu_meta,for_e_lu)
  PRINT *,'for_e_lu read'

  CALL netcdf_get_var(TRIM(netcdf_filename),emissivity_lu_meta,emissivity_lu)
  PRINT *,'emissivity_lu read'

  CALL netcdf_get_var(TRIM(netcdf_filename),root_lu_meta,root_lu)
  PRINT *,'root_lu read'


 END SUBROUTINE read_netcdf_buffer_ecoclimap


!-----------------------------------------------------------------------


    !> netcdf output of GLC2000 derived buffer fields
    SUBROUTINE write_netcdf_buffer_glc2000(netcdf_filename,  &
    &                                     tg,         &
    &                                     undefined, &
    &                                     undef_int,   &
    &                                     lon_geo,     &
    &                                     lat_geo, &
    &                                     fr_land_glc2000, &
    &                                     glc2000_class_fraction,    &
    &                                     glc2000_class_npixel, &
    &                                     glc2000_tot_npixel, &
    &                                     ice_glc2000, &
    &                                     z0_glc2000, &
    &                                     root_glc2000, &
    &                                     plcov_mn_glc2000, &
    &                                     plcov_mx_glc2000, &
    &                                     lai_mn_glc2000, &
    &                                     lai_mx_glc2000, &
    &                                     rs_min_glc2000, &
    &                                     urban_glc2000,  &
    &                                     for_d_glc2000,  &
    &                                     for_e_glc2000, &
    &                                     emissivity_glc2000)

  USE mo_var_meta_data, ONLY: dim_3d_tg, &
    &                         def_dimension_info_buffer


  USE mo_var_meta_data, ONLY: lon_geo_meta, &
    &                         lat_geo_meta, &
    &                         no_raw_data_pixel_meta, &
    &                         def_com_target_fields_meta  
  
  USE mo_var_meta_data, ONLY: def_glc2000_fields_meta

  USE mo_var_meta_data, ONLY: dim_glc2000_tg

  USE mo_var_meta_data, ONLY: fr_land_glc2000_meta, glc2000_tot_npixel_meta, &
    &       glc2000_class_fraction_meta, glc2000_class_npixel_meta, &
    &       ice_glc2000_meta, z0_glc2000_meta, &
    &       plcov_mx_glc2000_meta, plcov_mn_glc2000_meta, &
    &       lai_mx_glc2000_meta, lai_mn_glc2000_meta, &
    &       rs_min_glc2000_meta, urban_glc2000_meta, &
    &       for_d_glc2000_meta, for_e_glc2000_meta, &
    &       emissivity_glc2000_meta, root_glc2000_meta

  CHARACTER (len=*), INTENT(IN)      :: netcdf_filename !< filename for the netcdf file
  TYPE(target_grid_def), INTENT(IN) :: tg !< structure with target grid description
  REAL(KIND=wp), INTENT(IN)          :: undefined       !< value to indicate undefined grid elements 
  INTEGER, INTENT(IN)                :: undef_int       !< value to indicate undefined grid elements
  REAL (KIND=wp), INTENT(IN) :: lon_geo(:,:,:)  !< longitude coordinates of the target grid in the geographical system
  REAL (KIND=wp), INTENT(IN) :: lat_geo(:,:,:)  !< latitude coordinates of the target grid in the geographical system
  REAL (KIND=wp), INTENT(IN)  :: glc2000_class_fraction(:,:,:,:)  
                                 !< fraction for each glc2000 class on target grid (dim (ie,je,ke,nclass_glc2000))
  INTEGER (KIND=i8), INTENT(IN) :: glc2000_class_npixel(:,:,:,:) 
                       !< number of raw data pixels for each glc2000 class on target grid (dim (ie,je,ke,nclass_glc2000))
  INTEGER (KIND=i8), INTENT(IN) :: glc2000_tot_npixel(:,:,:)  
                                   !< total number of glc2000 raw data pixels on target grid (dimension (ie,je,ke))
  REAL (KIND=wp), INTENT(IN)  :: fr_land_glc2000(:,:,:) !< fraction land due to glc2000 raw data
  REAL (KIND=wp), INTENT(IN)  :: ice_glc2000(:,:,:)     !< fraction of ice due to glc2000 raw data
  REAL (KIND=wp), INTENT(IN)  :: z0_glc2000(:,:,:)      !< roughness length due to glc2000 land use data
  REAL (KIND=wp), INTENT(IN)  :: root_glc2000(:,:,:)    !< root depth due to glc2000 land use data
  REAL (KIND=wp), INTENT(IN)  :: plcov_mx_glc2000(:,:,:)!< plant cover maximum due to glc2000 land use data
  REAL (KIND=wp), INTENT(IN)  :: plcov_mn_glc2000(:,:,:)!< plant cover minimum due to glc2000 land use data
  REAL (KIND=wp), INTENT(IN)  :: lai_mx_glc2000(:,:,:)  !< Leaf Area Index maximum due to glc2000 land use data
  REAL (KIND=wp), INTENT(IN)  :: lai_mn_glc2000(:,:,:)  !< Leaf Area Index minimum due to glc2000 land use data
  REAL (KIND=wp), INTENT(IN)  :: rs_min_glc2000(:,:,:)  !< minimal stomata resistance due to glc2000 land use data
  REAL (KIND=wp), INTENT(IN)  :: urban_glc2000(:,:,:)   !< urban fraction due to glc2000 land use data
  REAL (KIND=wp), INTENT(IN)  :: for_d_glc2000(:,:,:)   !< deciduous forest (fraction) due to glc2000 land use data
  REAL (KIND=wp), INTENT(IN)  :: for_e_glc2000(:,:,:)   !< evergreen forest (fraction) due to glc2000 land use data
  REAL (KIND=wp), INTENT(IN)  :: emissivity_glc2000(:,:,:) !< longwave emissivity due to glc2000 land use data


  ! local variables
  INTEGER :: ndims  
  INTEGER :: ncid
  INTEGER (KIND=i8) :: undefined_i

  TYPE(dim_meta_info), ALLOCATABLE :: dim_list(:) !< dimensions for netcdf file

  INTEGER, PARAMETER :: nglob_atts=6
  TYPE(netcdf_attributes) :: global_attributes(nglob_atts)

  INTEGER :: errorcode !< error status variable

  INTEGER :: n !< counter

  PRINT *,'ENTER write_netcdf_buffer_glc2000'

  !-------------------------------------------------------------
  ! define global attributes
  CALL set_global_att_glc2000(global_attributes)

  !set up dimensions for buffer
  CALL  def_dimension_info_buffer(tg)
  ! dim_3d_tg

  ! define meta information for various land use related variables (GLC2000) for netcdf output
  CALL def_glc2000_fields_meta(tg,nclass_glc2000,dim_3d_tg)
  ! dim_glc2000_tg
  ! fr_land_glc2000_meta, glc2000_tot_npixel_meta, &
  !  &       glc2000_class_fraction_meta, glc2000_class_npixel_meta, &
  !  &       ice_glc2000_meta, z0_glc2000_meta, &
  !  &       plcov_mx_glc2000_meta, plcov_mn_glc2000_meta, &
  !  &       lai_mx_glc2000_meta, lai_mn_glc2000_meta, &
  !  &       rs_min_glc2000_meta, urban_glc2000_meta, &
  !  &       for_d_glc2000_meta, for_e_glc2000_meta, &
  !  &       emissivity_glc2000_meta, root_glc2000_meta

  ! define meta information for target field variables lon_geo, lat_geo 
  CALL def_com_target_fields_meta(dim_3d_tg)
  ! lon_geo_meta and lat_geo_meta

  !set up dimensions for buffer netcdf output 
  ndims = 4
  ALLOCATE(dim_list(1:ndims),STAT=errorcode)
  IF (errorcode /= 0 ) CALL abort_extpar('Cant allocate array dim_list')
  dim_list = dim_glc2000_tg

  !  dim_3d_buffer(:) = dim_list(:)

  undefined_i = undef_int

  !-----------------------------------------------------------------

  CALL open_new_netcdf_file(netcdf_filename=TRIM(netcdf_filename),   &
      &                       dim_list=dim_list,                  &
      &                       global_attributes=global_attributes, &
      &                       ncid=ncid)

  ! lon
  CALL netcdf_put_var(ncid,lon_geo,lon_geo_meta,undefined)

  ! lat
  CALL netcdf_put_var(ncid,lat_geo,lat_geo_meta,undefined)

  ! fr_land_glc2000
  CALL netcdf_put_var(ncid,fr_land_glc2000,fr_land_glc2000_meta,undefined)

  ! ice_glc2000
  CALL netcdf_put_var(ncid,ice_glc2000,ice_glc2000_meta,undefined)

  ! plcov_mx_glc2000
  CALL netcdf_put_var(ncid,plcov_mx_glc2000,plcov_mx_glc2000_meta,undefined)

  ! lai_mx_glc2000
  CALL netcdf_put_var(ncid,lai_mx_glc2000,lai_mx_glc2000_meta,undefined)

  ! rs_min_glc2000
  CALL netcdf_put_var(ncid,rs_min_glc2000,rs_min_glc2000_meta,undefined)

  ! urban_glc2000
  CALL netcdf_put_var(ncid,urban_glc2000,urban_glc2000_meta,undefined)

  ! for_d_glc2000
  CALL netcdf_put_var(ncid,for_d_glc2000,for_d_glc2000_meta,undefined)

  ! for_e_glc2000
  CALL netcdf_put_var(ncid,for_e_glc2000,for_e_glc2000_meta,undefined)

  ! emissivity_glc2000
  CALL netcdf_put_var(ncid,emissivity_glc2000,emissivity_glc2000_meta,undefined)

  ! root_glc2000
  CALL netcdf_put_var(ncid,root_glc2000,root_glc2000_meta,undefined)

  ! z0_glc2000
  CALL netcdf_put_var(ncid,z0_glc2000,z0_glc2000_meta,undefined)

  ! lai_mn_glc2000
  CALL netcdf_put_var(ncid,lai_mn_glc2000,lai_mn_glc2000_meta,undefined)
    
  ! plcov_mn_glc2000
  CALL netcdf_put_var(ncid,plcov_mn_glc2000,plcov_mn_glc2000_meta,undefined)

  ! glc2000_tot_npixel
  CALL netcdf_put_var(ncid,glc2000_tot_npixel,glc2000_tot_npixel_meta,undefined_i)

  ! glc2000_class_fraction
  CALL netcdf_put_var(ncid,glc2000_class_fraction,glc2000_class_fraction_meta,undefined)
  !-----------------------------------------------------------------

  ! glc2000_class_npixel
  CALL netcdf_put_var(ncid,glc2000_class_npixel,glc2000_class_npixel_meta,undefined_i)
  !-----------------------------------------------------------------


  CALL close_netcdf_file(ncid)

  END SUBROUTINE write_netcdf_buffer_glc2000
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
   !> netcdf output of GLC2000 derived COSMO fields
  SUBROUTINE write_netcdf_cosmo_grid_glc2000(netcdf_filename,  &
    &                                     cosmo_grid,       &
    &                                     tg,         &
    &                                     undefined, &
    &                                     undef_int,   &
    &                                     lon_geo,     &
    &                                     lat_geo, &
    &                                     fr_land_glc2000, &
    &                                     glc2000_class_fraction,    &
    &                                     glc2000_class_npixel, &
    &                                     glc2000_tot_npixel, &
    &                                     ice_glc2000, &
    &                                     z0_glc2000, &
    &                                     root_glc2000, &
    &                                     plcov_mn_glc2000, &
    &                                     plcov_mx_glc2000, &
    &                                     lai_mn_glc2000, &
    &                                     lai_mx_glc2000, &
    &                                     rs_min_glc2000, &
    &                                     urban_glc2000,  &
    &                                     for_d_glc2000,  &
    &                                     for_e_glc2000, &
    &                                     emissivity_glc2000)
  

  USE mo_var_meta_data, ONLY: dim_3d_tg, &
    &                         def_dimension_info_buffer


  USE mo_var_meta_data, ONLY: lon_geo_meta, &
    &                         lat_geo_meta, &
    &                         no_raw_data_pixel_meta, &
    &                         def_com_target_fields_meta  
   
  USE mo_var_meta_data, ONLY: nc_grid_def_cosmo, &
    &                         set_nc_grid_def_cosmo
    
  USE mo_var_meta_data, ONLY: dim_rlon_cosmo, &
    &                         dim_rlat_cosmo, &
    &                         dim_2d_cosmo,   &
    &                         rlon_meta,      &
    &                         rlat_meta,      &
    &                         def_dimension_info_cosmo

  USE mo_cosmo_grid, ONLY: lon_rot, lat_rot

  USE mo_var_meta_data, ONLY: def_glc2000_fields_meta

  USE mo_var_meta_data, ONLY: dim_glc2000_tg

  USE mo_var_meta_data, ONLY: fr_land_glc2000_meta, glc2000_tot_npixel_meta, &
    &       glc2000_class_fraction_meta, glc2000_class_npixel_meta, &
    &       ice_glc2000_meta, z0_glc2000_meta, &
    &       plcov_mx_glc2000_meta, plcov_mn_glc2000_meta, &
    &       lai_mx_glc2000_meta, lai_mn_glc2000_meta, &
    &       rs_min_glc2000_meta, urban_glc2000_meta, &
    &       for_d_glc2000_meta, for_e_glc2000_meta, &
    &       emissivity_glc2000_meta, root_glc2000_meta

  CHARACTER (len=*), INTENT(IN)      :: netcdf_filename !< filename for the netcdf file
  TYPE(rotated_lonlat_grid), INTENT(IN) :: cosmo_grid !< structure which contains the definition of the COSMO grid
  TYPE(target_grid_def), INTENT(IN) :: tg !< structure with target grid description
  REAL(KIND=wp), INTENT(IN)          :: undefined       !< value to indicate undefined grid elements 
  INTEGER, INTENT(IN)                :: undef_int       !< value to indicate undefined grid elements
  REAL (KIND=wp), INTENT(IN) :: lon_geo(:,:,:)  !< longitude coordinates of the target grid in the geographical system
  REAL (KIND=wp), INTENT(IN) :: lat_geo(:,:,:)  !< latitude coordinates of the target grid in the geographical system
  REAL (KIND=wp), INTENT(IN)  :: glc2000_class_fraction(:,:,:,:)  
                                 !< fraction for each glc2000 class on target grid (dim (ie,je,ke,nclass_glc2000))
  INTEGER (KIND=i8), INTENT(IN) :: glc2000_class_npixel(:,:,:,:) 
                         !< number of raw data pixels for each glc2000 class on target grid (dim (ie,je,ke,nclass_glc2000))
  INTEGER (KIND=i8), INTENT(IN) :: glc2000_tot_npixel(:,:,:)  
                                   !< total number of glc2000 raw data pixels on target grid (dimension (ie,je,ke))
  REAL (KIND=wp), INTENT(IN)  :: fr_land_glc2000(:,:,:) !< fraction land due to glc2000 raw data
  REAL (KIND=wp), INTENT(IN)  :: ice_glc2000(:,:,:)     !< fraction of ice due to glc2000 raw data
  REAL (KIND=wp), INTENT(IN)  :: z0_glc2000(:,:,:)      !< roughness length due to glc2000 land use data
  REAL (KIND=wp), INTENT(IN)  :: root_glc2000(:,:,:)    !< root depth due to glc2000 land use data
  REAL (KIND=wp), INTENT(IN)  :: plcov_mx_glc2000(:,:,:)!< plant cover maximum due to glc2000 land use data
  REAL (KIND=wp), INTENT(IN)  :: plcov_mn_glc2000(:,:,:)!< plant cover minimum due to glc2000 land use data
  REAL (KIND=wp), INTENT(IN)  :: lai_mx_glc2000(:,:,:)  !< Leaf Area Index maximum due to glc2000 land use data
  REAL (KIND=wp), INTENT(IN)  :: lai_mn_glc2000(:,:,:)  !< Leaf Area Index minimum due to glc2000 land use data
  REAL (KIND=wp), INTENT(IN)  :: rs_min_glc2000(:,:,:)  !< minimal stomata resistance due to glc2000 land use data
  REAL (KIND=wp), INTENT(IN)  :: urban_glc2000(:,:,:)   !< urban fraction due to glc2000 land use data
  REAL (KIND=wp), INTENT(IN)  :: for_d_glc2000(:,:,:)   !< deciduous forest (fraction) due to glc2000 land use data
  REAL (KIND=wp), INTENT(IN)  :: for_e_glc2000(:,:,:)   !< evergreen forest (fraction) due to glc2000 land use data
  REAL (KIND=wp), INTENT(IN)  :: emissivity_glc2000(:,:,:) !< longwave emissivity due to glc2000 land use data


  ! local variables

  INTEGER :: ndims  

  INTEGER :: ncid
  INTEGER :: varid
  INTEGER (KIND=i8) :: undefined_i

  TYPE(dim_meta_info), ALLOCATABLE :: dim_list(:) !< dimensions for netcdf file
  TYPE(dim_meta_info), TARGET :: dim_nclass(1:3)

  INTEGER, PARAMETER :: nglob_atts=6
  TYPE(netcdf_attributes) :: global_attributes(nglob_atts)

  INTEGER :: errorcode !< error status variable

  CHARACTER (len=80):: grid_mapping !< netcdf attribute grid mapping
  CHARACTER (len=80):: coordinates  !< netcdf attribute coordinates

  INTEGER :: n !< counter

  PRINT *,'Enter write_netcdf_cosmo_grid_glc2000'

  !-------------------------------------------------------------
  ! define global attributes
  CALL set_global_att_glc2000(global_attributes)

  !set up dimensions for buffer
  CALL  def_dimension_info_buffer(tg)
  ! dim_3d_tg

   !set up dimensions for COSMO grid
  CALL def_dimension_info_cosmo(cosmo_grid)
  ! dim_rlon_cosmo, dim_rlat_cosmo, dim_2d_cosmo, rlon_meta, rlat_meta

  ! set mapping parameters for netcdf
  grid_mapping="rotated_pole"
  coordinates="lon lat"
  CALL set_nc_grid_def_cosmo(cosmo_grid,grid_mapping)
  ! nc_grid_def_cosmo

  ! define meta information for various land use related variables (GLC2000) for netcdf output
  CALL def_glc2000_fields_meta(tg,nclass_glc2000,dim_2d_cosmo,coordinates,grid_mapping)
  ! dim_glc2000_tg
  ! fr_land_glc2000_meta, glc2000_tot_npixel_meta, &
  !  &       glc2000_class_fraction_meta, glc2000_class_npixel_meta, &
  !  &       ice_glc2000_meta, z0_glc2000_meta, &
  !  &       plcov_mx_glc2000_meta, plcov_mn_glc2000_meta, &
  !  &       lai_mx_glc2000_meta, lai_mn_glc2000_meta, &
  !  &       rs_min_glc2000_meta, urban_glc2000_meta, &
  !  &       for_d_glc2000_meta, for_e_glc2000_meta, &
  !  &       emissivity_glc2000_meta, root_glc2000_meta

  ! define meta information for target field variables lon_geo, lat_geo 
  CALL def_com_target_fields_meta(dim_2d_cosmo,coordinates,grid_mapping)
  ! lon_geo_meta and lat_geo_meta

  !set up dimensions for buffer netcdf output 
  ndims = 3
  ALLOCATE(dim_list(1:ndims),STAT=errorcode)
  IF (errorcode /= 0 ) CALL abort_extpar('Cant allocate array dim_list')

  dim_list(1) = dim_rlon_cosmo(1) ! rlon
  dim_list(2) = dim_rlat_cosmo(1) ! rlat
  dim_list(3)%dimname = 'nclass'
  dim_list(3)%dimsize = nclass_glc2000

  dim_nclass = dim_list

  undefined_i = undef_int

  !-----------------------------------------------------------------
   CALL open_new_netcdf_file(netcdf_filename=TRIM(netcdf_filename),   &
      &                       dim_list=dim_list,                  &
      &                       global_attributes=global_attributes, &
      &                       ncid=ncid)
    !-----------------------------------------------------------------

    ! rlon
    CALL netcdf_put_var(ncid,lon_rot(1:cosmo_grid%nlon_rot),rlon_meta,undefined)

    CALL netcdf_put_var(ncid,lat_rot(1:cosmo_grid%nlat_rot),rlat_meta,undefined)

    ! lon
  CALL netcdf_put_var(ncid,lon_geo(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1), &
    &                 lon_geo_meta,undefined)

  ! lat
  CALL netcdf_put_var(ncid,lat_geo(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1), &
    &                 lat_geo_meta,undefined)

  ! fr_land_glc2000
  CALL netcdf_put_var(ncid,fr_land_glc2000(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1), &
    &                 fr_land_glc2000_meta,undefined)

  ! ice_glc2000
  CALL netcdf_put_var(ncid,ice_glc2000(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1), &
    &                 ice_glc2000_meta,undefined)

  ! plcov_mx_glc2000
  CALL netcdf_put_var(ncid,plcov_mx_glc2000(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1), &
    &                 plcov_mx_glc2000_meta,undefined)

  ! lai_mx_glc2000
  CALL netcdf_put_var(ncid,lai_mx_glc2000(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1), &
    &                 lai_mx_glc2000_meta,undefined)

  ! rs_min_glc2000
  CALL netcdf_put_var(ncid,rs_min_glc2000(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1), &
    &                 rs_min_glc2000_meta,undefined)

  ! urban_glc2000
  CALL netcdf_put_var(ncid,urban_glc2000(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1), &
    &                 urban_glc2000_meta,undefined)

  ! for_d_glc2000
  CALL netcdf_put_var(ncid,for_d_glc2000(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1), &
    &                 for_d_glc2000_meta,undefined)

  ! for_e_glc2000
  CALL netcdf_put_var(ncid,for_e_glc2000(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1), &
    &                 for_e_glc2000_meta,undefined)

  ! emissivity_glc2000
  CALL netcdf_put_var(ncid,emissivity_glc2000(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1), &
    &                 emissivity_glc2000_meta,undefined)

  ! root_glc2000
  CALL netcdf_put_var(ncid,root_glc2000(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1), &
    &                 root_glc2000_meta,undefined)

  ! z0_glc2000
  CALL netcdf_put_var(ncid,z0_glc2000(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1), &
    &                 z0_glc2000_meta,undefined)

  ! lai_mn_glc2000
  CALL netcdf_put_var(ncid,lai_mn_glc2000(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1), &
    &                 lai_mn_glc2000_meta,undefined)
    
  ! plcov_mn_glc2000
  CALL netcdf_put_var(ncid,plcov_mn_glc2000(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1), &
    &                 plcov_mn_glc2000_meta,undefined)

  ! glc2000_tot_npixel
  CALL netcdf_put_var(ncid,glc2000_tot_npixel(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1), &
    &                 glc2000_tot_npixel_meta,undefined_i)

  ! glc2000_class_fraction
  CALL netcdf_put_var(ncid,glc2000_class_fraction(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1,1:nclass_glc2000), &
    &                 glc2000_class_fraction_meta,undefined)
  !-----------------------------------------------------------------

  ! glc2000_class_npixel
  CALL netcdf_put_var(ncid,glc2000_class_npixel(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1,1:nclass_glc2000), &
    &                 glc2000_class_npixel_meta,undefined_i)
  !-----------------------------------------------------------------
    CALL netcdf_def_grid_mapping(ncid, nc_grid_def_cosmo, varid)

    CALL close_netcdf_file(ncid)


  END SUBROUTINE write_netcdf_cosmo_grid_glc2000
  !-----------------------------------------------------------------------


  !> netcdf output of GLC2000 derived ICON fields
  SUBROUTINE write_netcdf_icon_grid_glc2000(netcdf_filename,  &
    &                                     icon_grid,       &
    &                                     tg,         &
    &                                     undefined, &
    &                                     undef_int,   &
    &                                     lon_geo,     &
    &                                     lat_geo, &
    &                                     fr_land_glc2000, &
    &                                     glc2000_class_fraction,    &
    &                                     glc2000_class_npixel, &
    &                                     glc2000_tot_npixel, &
    &                                     ice_glc2000, &
    &                                     z0_glc2000, &
    &                                     root_glc2000, &
    &                                     plcov_mn_glc2000, &
    &                                     plcov_mx_glc2000, &
    &                                     lai_mn_glc2000, &
    &                                     lai_mx_glc2000, &
    &                                     rs_min_glc2000, &
    &                                     urban_glc2000,  &
    &                                     for_d_glc2000,  &
    &                                     for_e_glc2000, &
    &                                     emissivity_glc2000)

  USE mo_var_meta_data, ONLY: dim_3d_tg, &
    &                         def_dimension_info_buffer


  USE mo_var_meta_data, ONLY: lon_geo_meta, &
    &                         lat_geo_meta, &
    &                         no_raw_data_pixel_meta, &
    &                         def_com_target_fields_meta  
   

  USE mo_var_meta_data, ONLY:  dim_icon, &
    &                          def_dimension_info_icon

  USE mo_var_meta_data, ONLY: nc_grid_def_icon, &
    &                         set_nc_grid_def_icon

  USE mo_var_meta_data, ONLY: def_glc2000_fields_meta

  USE mo_var_meta_data, ONLY: dim_glc2000_tg

  USE mo_var_meta_data, ONLY: fr_land_glc2000_meta, glc2000_tot_npixel_meta, &
    &       glc2000_class_fraction_meta, glc2000_class_npixel_meta, &
    &       ice_glc2000_meta, z0_glc2000_meta, &
    &       plcov_mx_glc2000_meta, plcov_mn_glc2000_meta, &
    &       lai_mx_glc2000_meta, lai_mn_glc2000_meta, &
    &       rs_min_glc2000_meta, urban_glc2000_meta, &
    &       for_d_glc2000_meta, for_e_glc2000_meta, &
    &       emissivity_glc2000_meta, root_glc2000_meta

  CHARACTER (len=*), INTENT(IN)      :: netcdf_filename !< filename for the netcdf file
  TYPE(icosahedral_triangular_grid), INTENT(IN) :: icon_grid !< structure which contains the definition of the ICON grid
  TYPE(target_grid_def), INTENT(IN) :: tg !< structure with target grid description
  REAL(KIND=wp), INTENT(IN)          :: undefined       !< value to indicate undefined grid elements 
  INTEGER, INTENT(IN)                :: undef_int       !< value to indicate undefined grid elements
  REAL (KIND=wp), INTENT(IN) :: lon_geo(:,:,:)  !< longitude coordinates of the target grid in the geographical system
  REAL (KIND=wp), INTENT(IN) :: lat_geo(:,:,:)  !< latitude coordinates of the target grid in the geographical system
  REAL (KIND=wp), INTENT(IN)  :: glc2000_class_fraction(:,:,:,:)  
                                 !< fraction for each glc2000 class on target grid (dim (ie,je,ke,nclass_glc2000))
  INTEGER (KIND=i8), INTENT(IN) :: glc2000_class_npixel(:,:,:,:) 
                           !< number of raw data pixels for each glc2000 class on target grid (dim (ie,je,ke,nclass_glc2000))
  INTEGER (KIND=i8), INTENT(IN) :: glc2000_tot_npixel(:,:,:)  
                                   !< total number of glc2000 raw data pixels on target grid (dim (ie,je,ke))
  REAL (KIND=wp), INTENT(IN)  :: fr_land_glc2000(:,:,:) !< fraction land due to glc2000 raw data
  REAL (KIND=wp), INTENT(IN)  :: ice_glc2000(:,:,:)     !< fraction of ice due to glc2000 raw data
  REAL (KIND=wp), INTENT(IN)  :: z0_glc2000(:,:,:)      !< roughness length due to glc2000 land use data
  REAL (KIND=wp), INTENT(IN)  :: root_glc2000(:,:,:)    !< root depth due to glc2000 land use data
  REAL (KIND=wp), INTENT(IN)  :: plcov_mx_glc2000(:,:,:)!< plant cover maximum due to glc2000 land use data
  REAL (KIND=wp), INTENT(IN)  :: plcov_mn_glc2000(:,:,:)!< plant cover minimum due to glc2000 land use data
  REAL (KIND=wp), INTENT(IN)  :: lai_mx_glc2000(:,:,:)  !< Leaf Area Index maximum due to glc2000 land use data
  REAL (KIND=wp), INTENT(IN)  :: lai_mn_glc2000(:,:,:)  !< Leaf Area Index minimum due to glc2000 land use data
  REAL (KIND=wp), INTENT(IN)  :: rs_min_glc2000(:,:,:)  !< minimal stomata resistance due to glc2000 land use data
  REAL (KIND=wp), INTENT(IN)  :: urban_glc2000(:,:,:)   !< urban fraction due to glc2000 land use data
  REAL (KIND=wp), INTENT(IN)  :: for_d_glc2000(:,:,:)   !< deciduous forest (fraction) due to glc2000 land use data
  REAL (KIND=wp), INTENT(IN)  :: for_e_glc2000(:,:,:)   !< evergreen forest (fraction) due to glc2000 land use data
  REAL (KIND=wp), INTENT(IN)  :: emissivity_glc2000(:,:,:) !< longwave emissivity due to glc2000 land use data


  ! local variables


  INTEGER :: ndims 
  INTEGER :: ncid
  INTEGER :: varid
  INTEGER (KIND=i8) :: undefined_i

  TYPE(dim_meta_info), ALLOCATABLE :: dim_list(:) !< dimensions for netcdf file
  TYPE(dim_meta_info), TARGET :: dim_2d_icon(1:2)
  TYPE(dim_meta_info), TARGET :: dim_1d_icon(1:1)

  INTEGER, PARAMETER :: nglob_atts=6
  TYPE(netcdf_attributes) :: global_attributes(nglob_atts)

  INTEGER :: errorcode !< error status variable

  CHARACTER (len=80):: grid_mapping !< netcdf attribute grid mapping
  CHARACTER (len=80):: coordinates  !< netcdf attribute coordinates

  INTEGER :: n !< counter

  !-------------------------------------------------------------
  ! define global attributes
  CALL set_global_att_glc2000(global_attributes)

  !set up dimensions for buffer
  CALL  def_dimension_info_buffer(tg)
  ! dim_3d_tg

  !set up dimensions for ICON grid
  CALL def_dimension_info_icon(icon_grid)
  ! dim_icon

  ! define meta information for various land use related variables (GLC2000) for netcdf output
  CALL def_glc2000_fields_meta(tg,nclass_glc2000,dim_icon)
  ! dim_glc2000_tg
  ! fr_land_glc2000_meta, glc2000_tot_npixel_meta, &
  !  &       glc2000_class_fraction_meta, glc2000_class_npixel_meta, &
  !  &       ice_glc2000_meta, z0_glc2000_meta, &
  !  &       plcov_mx_glc2000_meta, plcov_mn_glc2000_meta, &
  !  &       lai_mx_glc2000_meta, lai_mn_glc2000_meta, &
  !  &       rs_min_glc2000_meta, urban_glc2000_meta, &
  !  &       for_d_glc2000_meta, for_e_glc2000_meta, &
  !  &       emissivity_glc2000_meta, root_glc2000_meta

  ! define meta information for target field variables lon_geo, lat_geo 
  CALL def_com_target_fields_meta(dim_icon)
  ! lon_geo_meta and lat_geo_meta

  ! set mapping parameters for netcdf
  grid_mapping="lon_lat_on_sphere"
  coordinates="lon lat"

  CALL set_nc_grid_def_icon(grid_mapping)
  ! nc_grid_def_icon


  !set up dimensions for buffer netcdf output 
  ndims = 2
  ALLOCATE(dim_list(1:ndims),STAT=errorcode)
  IF (errorcode /= 0 ) CALL abort_extpar('Cant allocate array dim_list')

  dim_list(1) = dim_icon(1) ! cell
  dim_list(2)%dimname = 'nclass'
  dim_list(2)%dimsize = nclass_glc2000

  dim_1d_icon =  dim_icon(1) ! cell
  dim_2d_icon = dim_list

  undefined_i = undef_int

  !-----------------------------------------------------------------
   CALL open_new_netcdf_file(netcdf_filename=TRIM(netcdf_filename),   &
      &                       dim_list=dim_list,                  &
      &                       global_attributes=global_attributes, &
      &                       ncid=ncid)
  !-----------------------------------------------------------------

   ! lon
    CALL netcdf_put_var(ncid,lon_geo(1:icon_grid%ncell,1,1),lon_geo_meta,undefined)

    ! lat
    CALL netcdf_put_var(ncid,lat_geo(1:icon_grid%ncell,1,1),lat_geo_meta,undefined)

    !-----------------------------------------------------------------
    n=1 ! fr_land_glc2000
    CALL netcdf_put_var(ncid,fr_land_glc2000(1:icon_grid%ncell,1,1),fr_land_glc2000_meta,undefined)

    n=2 ! ice_glc2000
    CALL netcdf_put_var(ncid,ice_glc2000(1:icon_grid%ncell,1,1),ice_glc2000_meta,undefined)

    n=3 ! plcov_mx_glc2000
    CALL netcdf_put_var(ncid,plcov_mx_glc2000(1:icon_grid%ncell,1,1),plcov_mx_glc2000_meta,undefined)

    n=4 ! lai_mx_glc2000
    CALL netcdf_put_var(ncid,lai_mx_glc2000(1:icon_grid%ncell,1,1),lai_mx_glc2000_meta,undefined)

    n=5 ! rs_min_glc2000
    CALL netcdf_put_var(ncid,rs_min_glc2000(1:icon_grid%ncell,1,1),rs_min_glc2000_meta,undefined)

    n=6 ! urban_glc2000
    CALL netcdf_put_var(ncid,urban_glc2000(1:icon_grid%ncell,1,1),urban_glc2000_meta,undefined)

    n=7 ! for_d_glc2000
    CALL netcdf_put_var(ncid,for_d_glc2000(1:icon_grid%ncell,1,1),for_d_glc2000_meta,undefined)

    n=8 ! for_e_glc2000
    CALL netcdf_put_var(ncid,for_e_glc2000(1:icon_grid%ncell,1,1),for_e_glc2000_meta,undefined)

    n=9 ! emissivity_glc2000
    CALL netcdf_put_var(ncid, emissivity_glc2000(1:icon_grid%ncell,1,1),emissivity_glc2000_meta,undefined)

    n=10 ! root_glc2000
    CALL netcdf_put_var(ncid,root_glc2000(1:icon_grid%ncell,1,1),root_glc2000_meta,undefined)

    n=11 ! z0_glc2000
    CALL netcdf_put_var(ncid,z0_glc2000(1:icon_grid%ncell,1,1),z0_glc2000_meta,undefined)

  !-----------------------------------------------------------------
  CALL close_netcdf_file(ncid)


  END SUBROUTINE write_netcdf_icon_grid_glc2000
  !-----------------------------------------------------------------------

  !----------------------------------------------------------------------- 
  !-----------------------------------------------------------------------
  !> set global attributes for netcdf with glc2000 data
  SUBROUTINE set_global_att_glc2000(global_attributes)
    USE mo_glc2000_lookup_tables, ONLY: name_lookup_table_glc2000
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
    global_attributes(3)%attributetext='GLC2000 data'

    CALL DATE_AND_TIME(ydate,ytime)
    READ(ydate,'(4A2)') cc,yy,mm,dd
    READ(ytime,'(2A2)') hh, minute

    ydate=TRIM(cc)//TRIM(yy)//'-'//TRIM(mm)//'-'//TRIM(dd)
    ytime=TRIM(hh)//':'//TRIM(minute) 

    global_attributes(4)%attname = 'history'
    global_attributes(4)%attributetext=TRIM(ydate)//'T'//TRIM(ytime)//' glc2000_to_buffer'

    global_attributes(5)%attname = 'references'
    global_attributes(5)%attributetext=''

    global_attributes(6)%attname = 'comment'
    global_attributes(6)%attributetext='Landuse data look-up table: '//TRIM(name_lookup_table_glc2000)

  END SUBROUTINE set_global_att_glc2000
  !-----------------------------------------------------------------------
    !> netcdf output of GLCC derived buffer fields
    SUBROUTINE write_netcdf_buffer_glcc(netcdf_filename,  &
    &                                     tg,         &
    &                                     undefined, &
    &                                     undef_int,   &
    &                                     lon_geo,     &
    &                                     lat_geo, &
    &                                     fr_land_glcc, &
    &                                     glcc_class_fraction,    &
    &                                     glcc_class_npixel, &
    &                                     glcc_tot_npixel, &
    &                                     ice_glcc, &
    &                                     z0_glcc, &
    &                                     root_glcc, &
    &                                     plcov_mn_glcc, &
    &                                     plcov_mx_glcc, &
    &                                     lai_mn_glcc, &
    &                                     lai_mx_glcc, &
    &                                     rs_min_glcc, &
    &                                     urban_glcc,  &
    &                                     for_d_glcc,  &
    &                                     for_e_glcc, &
    &                                     emissivity_glcc)

  USE mo_var_meta_data, ONLY: dim_3d_tg, &
    &                         def_dimension_info_buffer


  USE mo_var_meta_data, ONLY: lon_geo_meta, &
    &                         lat_geo_meta, &
    &                         no_raw_data_pixel_meta, &
    &                         def_com_target_fields_meta  
  
  USE mo_var_meta_data, ONLY: def_glcc_fields_meta

  USE mo_var_meta_data, ONLY: dim_glcc_tg

  USE mo_var_meta_data, ONLY: fr_land_glcc_meta, glcc_tot_npixel_meta, &
    &       glcc_class_fraction_meta, glcc_class_npixel_meta, &
    &       ice_glcc_meta, z0_glcc_meta, &
    &       plcov_mx_glcc_meta, plcov_mn_glcc_meta, &
    &       lai_mx_glcc_meta, lai_mn_glcc_meta, &
    &       rs_min_glcc_meta, urban_glcc_meta, &
    &       for_d_glcc_meta, for_e_glcc_meta, &
    &       emissivity_glcc_meta, root_glcc_meta

  CHARACTER (len=*), INTENT(IN)      :: netcdf_filename !< filename for the netcdf file
  TYPE(target_grid_def), INTENT(IN) :: tg !< structure with target grid description
  REAL(KIND=wp), INTENT(IN)          :: undefined       !< value to indicate undefined grid elements 
  INTEGER, INTENT(IN)                :: undef_int       !< value to indicate undefined grid elements
  REAL (KIND=wp), INTENT(IN) :: lon_geo(:,:,:)  !< longitude coordinates of the target grid in the geographical system
  REAL (KIND=wp), INTENT(IN) :: lat_geo(:,:,:)  !< latitude coordinates of the target grid in the geographical system
  !< fraction for each glcc class on target grid (dimension (ie,je,ke,nclass_glcc))
  REAL (KIND=wp), INTENT(IN)  :: glcc_class_fraction(:,:,:,:)  
  !< number of raw data pixels for each glcc class on target grid (dimension (ie,je,ke,nclass_glcc))
  INTEGER (KIND=i8), INTENT(IN) :: glcc_class_npixel(:,:,:,:) 
  !< total number of glcc raw data pixels on target grid (dimension (ie,je,ke))
  INTEGER (KIND=i8), INTENT(IN) :: glcc_tot_npixel(:,:,:)  
  REAL (KIND=wp), INTENT(IN)  :: fr_land_glcc(:,:,:) !< fraction land due to glcc raw data
  REAL (KIND=wp), INTENT(IN)  :: ice_glcc(:,:,:)     !< fraction of ice due to glcc raw data
  REAL (KIND=wp), INTENT(IN)  :: z0_glcc(:,:,:)      !< roughness length due to glcc land use data
  REAL (KIND=wp), INTENT(IN)  :: root_glcc(:,:,:)    !< root depth due to glcc land use data
  REAL (KIND=wp), INTENT(IN)  :: plcov_mx_glcc(:,:,:)!< plant cover maximum due to glcc land use data
  REAL (KIND=wp), INTENT(IN)  :: plcov_mn_glcc(:,:,:)!< plant cover minimum due to glcc land use data
  REAL (KIND=wp), INTENT(IN)  :: lai_mx_glcc(:,:,:)  !< Leaf Area Index maximum due to glcc land use data
  REAL (KIND=wp), INTENT(IN)  :: lai_mn_glcc(:,:,:)  !< Leaf Area Index minimum due to glcc land use data
  REAL (KIND=wp), INTENT(IN)  :: rs_min_glcc(:,:,:)  !< minimal stomata resistance due to glcc land use data
  REAL (KIND=wp), INTENT(IN)  :: urban_glcc(:,:,:)   !< urban fraction due to glcc land use data
  REAL (KIND=wp), INTENT(IN)  :: for_d_glcc(:,:,:)   !< deciduous forest (fraction) due to glcc land use data
  REAL (KIND=wp), INTENT(IN)  :: for_e_glcc(:,:,:)   !< evergreen forest (fraction) due to glcc land use data
  REAL (KIND=wp), INTENT(IN)  :: emissivity_glcc(:,:,:) !< longwave emissivity due to glcc land use data


  ! local variables
  INTEGER :: ndims  
  INTEGER :: ncid
  INTEGER (KIND=i8) :: undefined_i

  TYPE(dim_meta_info), ALLOCATABLE :: dim_list(:) !< dimensions for netcdf file

  INTEGER, PARAMETER :: nglob_atts=6
  TYPE(netcdf_attributes) :: global_attributes(nglob_atts)

  INTEGER :: errorcode !< error status variable

  INTEGER :: n !< counter

  PRINT *,'ENTER write_netcdf_buffer_glcc'

  !-------------------------------------------------------------
  ! define global attributes
  CALL set_global_att_glcc(global_attributes)

  !set up dimensions for buffer
  CALL  def_dimension_info_buffer(tg)
  ! dim_3d_tg

  ! define meta information for various land use related variables (GLCC) for netcdf output
  CALL def_glcc_fields_meta(tg,nclass_glcc,dim_3d_tg)
  ! dim_glcc_tg
  ! fr_land_glcc_meta, glcc_tot_npixel_meta, &
  !  &       glcc_class_fraction_meta, glcc_class_npixel_meta, &
  !  &       ice_glcc_meta, z0_glcc_meta, &
  !  &       plcov_mx_glcc_meta, plcov_mn_glcc_meta, &
  !  &       lai_mx_glcc_meta, lai_mn_glcc_meta, &
  !  &       rs_min_glcc_meta, urban_glcc_meta, &
  !  &       for_d_glcc_meta, for_e_glcc_meta, &
  !  &       emissivity_glcc_meta, root_glcc_meta

  ! define meta information for target field variables lon_geo, lat_geo 
  CALL def_com_target_fields_meta(dim_3d_tg)
  ! lon_geo_meta and lat_geo_meta

  !set up dimensions for buffer netcdf output 
  ndims = 4
  ALLOCATE(dim_list(1:ndims),STAT=errorcode)
  IF (errorcode /= 0 ) CALL abort_extpar('Cant allocate array dim_list')
  dim_list = dim_glcc_tg

  !  dim_3d_buffer(:) = dim_list(:)

  undefined_i = undef_int

  !-----------------------------------------------------------------

  CALL open_new_netcdf_file(netcdf_filename=TRIM(netcdf_filename),   &
      &                       dim_list=dim_list,                  &
      &                       global_attributes=global_attributes, &
      &                       ncid=ncid)

  ! lon
  CALL netcdf_put_var(ncid,lon_geo,lon_geo_meta,undefined)

  ! lat
  CALL netcdf_put_var(ncid,lat_geo,lat_geo_meta,undefined)

  ! fr_land_glcc
  CALL netcdf_put_var(ncid,fr_land_glcc,fr_land_glcc_meta,undefined)

  ! ice_glcc
  CALL netcdf_put_var(ncid,ice_glcc,ice_glcc_meta,undefined)

  ! plcov_mx_glcc
  CALL netcdf_put_var(ncid,plcov_mx_glcc,plcov_mx_glcc_meta,undefined)

  ! lai_mx_glcc
  CALL netcdf_put_var(ncid,lai_mx_glcc,lai_mx_glcc_meta,undefined)

  ! rs_min_glcc
  CALL netcdf_put_var(ncid,rs_min_glcc,rs_min_glcc_meta,undefined)

  ! urban_glcc
  CALL netcdf_put_var(ncid,urban_glcc,urban_glcc_meta,undefined)

  ! for_d_glcc
  CALL netcdf_put_var(ncid,for_d_glcc,for_d_glcc_meta,undefined)

  ! for_e_glcc
  CALL netcdf_put_var(ncid,for_e_glcc,for_e_glcc_meta,undefined)

  ! emissivity_glcc
  CALL netcdf_put_var(ncid,emissivity_glcc,emissivity_glcc_meta,undefined)

  ! root_glcc
  CALL netcdf_put_var(ncid,root_glcc,root_glcc_meta,undefined)

  ! z0_glcc
  CALL netcdf_put_var(ncid,z0_glcc,z0_glcc_meta,undefined)

  ! lai_mn_glcc
  CALL netcdf_put_var(ncid,lai_mn_glcc,lai_mn_glcc_meta,undefined)
    
  ! plcov_mn_glcc
  CALL netcdf_put_var(ncid,plcov_mn_glcc,plcov_mn_glcc_meta,undefined)

  ! glcc_tot_npixel
  CALL netcdf_put_var(ncid,glcc_tot_npixel,glcc_tot_npixel_meta,undefined_i)

  ! glcc_class_fraction
  CALL netcdf_put_var(ncid,glcc_class_fraction,glcc_class_fraction_meta,undefined)
  !-----------------------------------------------------------------

  ! glcc_class_npixel
  CALL netcdf_put_var(ncid,glcc_class_npixel,glcc_class_npixel_meta,undefined_i)
  !-----------------------------------------------------------------


  CALL close_netcdf_file(ncid)

  END SUBROUTINE write_netcdf_buffer_glcc
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------

  !> netcdf output of GLCC derived COSMO fields
  SUBROUTINE write_netcdf_cosmo_grid_glcc(netcdf_filename,  &
    &                                     cosmo_grid,       &
    &                                     tg,         &
    &                                     undefined, &
    &                                     undef_int,   &
    &                                     lon_geo,     &
    &                                     lat_geo, &
    &                                     fr_land_glcc, &
    &                                     glcc_class_fraction,    &
    &                                     glcc_class_npixel, &
    &                                     glcc_tot_npixel, &
    &                                     ice_glcc, &
    &                                     z0_glcc, &
    &                                     root_glcc, &
    &                                     plcov_mn_glcc, &
    &                                     plcov_mx_glcc, &
    &                                     lai_mn_glcc, &
    &                                     lai_mx_glcc, &
    &                                     rs_min_glcc, &
    &                                     urban_glcc,  &
    &                                     for_d_glcc,  &
    &                                     for_e_glcc, &
    &                                     emissivity_glcc)
  

  USE mo_var_meta_data, ONLY: dim_3d_tg, &
    &                         def_dimension_info_buffer


  USE mo_var_meta_data, ONLY: lon_geo_meta, &
    &                         lat_geo_meta, &
    &                         no_raw_data_pixel_meta, &
    &                         def_com_target_fields_meta  
   
  USE mo_var_meta_data, ONLY: nc_grid_def_cosmo, &
    &                         set_nc_grid_def_cosmo
    
  USE mo_var_meta_data, ONLY: dim_rlon_cosmo, &
    &                         dim_rlat_cosmo, &
    &                         dim_2d_cosmo,   &
    &                         rlon_meta,      &
    &                         rlat_meta,      &
    &                         def_dimension_info_cosmo

  USE mo_cosmo_grid, ONLY: lon_rot, lat_rot

  USE mo_var_meta_data, ONLY: def_glcc_fields_meta

  USE mo_var_meta_data, ONLY: dim_glcc_tg

  USE mo_var_meta_data, ONLY: fr_land_glcc_meta, glcc_tot_npixel_meta, &
    &       glcc_class_fraction_meta, glcc_class_npixel_meta, &
    &       ice_glcc_meta, z0_glcc_meta, &
    &       plcov_mx_glcc_meta, plcov_mn_glcc_meta, &
    &       lai_mx_glcc_meta, lai_mn_glcc_meta, &
    &       rs_min_glcc_meta, urban_glcc_meta, &
    &       for_d_glcc_meta, for_e_glcc_meta, &
    &       emissivity_glcc_meta, root_glcc_meta

  CHARACTER (len=*), INTENT(IN)      :: netcdf_filename !< filename for the netcdf file
  TYPE(rotated_lonlat_grid), INTENT(IN) :: cosmo_grid !< structure which contains the definition of the COSMO grid
  TYPE(target_grid_def), INTENT(IN) :: tg !< structure with target grid description
  REAL(KIND=wp), INTENT(IN)          :: undefined       !< value to indicate undefined grid elements 
  INTEGER, INTENT(IN)                :: undef_int       !< value to indicate undefined grid elements
  REAL (KIND=wp), INTENT(IN) :: lon_geo(:,:,:)  !< longitude coordinates of the target grid in the geographical system
  REAL (KIND=wp), INTENT(IN) :: lat_geo(:,:,:)  !< latitude coordinates of the target grid in the geographical system
  REAL (KIND=wp), INTENT(IN)  :: glcc_class_fraction(:,:,:,:)  
                                 !< fraction for each glcc class on target grid (dim (ie,je,ke,nclass_glcc))
  INTEGER (KIND=i8), INTENT(IN) :: glcc_class_npixel(:,:,:,:) 
                             !< number of raw data pixels for each glcc class on target grid (dim (ie,je,ke,nclass_glcc))
  INTEGER (KIND=i8), INTENT(IN) :: glcc_tot_npixel(:,:,:)  
                                   !< total number of glcc raw data pixels on target grid (dim (ie,je,ke))
  REAL (KIND=wp), INTENT(IN)  :: fr_land_glcc(:,:,:) !< fraction land due to glcc raw data
  REAL (KIND=wp), INTENT(IN)  :: ice_glcc(:,:,:)     !< fraction of ice due to glcc raw data
  REAL (KIND=wp), INTENT(IN)  :: z0_glcc(:,:,:)      !< roughness length due to glcc land use data
  REAL (KIND=wp), INTENT(IN)  :: root_glcc(:,:,:)    !< root depth due to glcc land use data
  REAL (KIND=wp), INTENT(IN)  :: plcov_mx_glcc(:,:,:)!< plant cover maximum due to glcc land use data
  REAL (KIND=wp), INTENT(IN)  :: plcov_mn_glcc(:,:,:)!< plant cover minimum due to glcc land use data
  REAL (KIND=wp), INTENT(IN)  :: lai_mx_glcc(:,:,:)  !< Leaf Area Index maximum due to glcc land use data
  REAL (KIND=wp), INTENT(IN)  :: lai_mn_glcc(:,:,:)  !< Leaf Area Index minimum due to glcc land use data
  REAL (KIND=wp), INTENT(IN)  :: rs_min_glcc(:,:,:)  !< minimal stomata resistance due to glcc land use data
  REAL (KIND=wp), INTENT(IN)  :: urban_glcc(:,:,:)   !< urban fraction due to glcc land use data
  REAL (KIND=wp), INTENT(IN)  :: for_d_glcc(:,:,:)   !< deciduous forest (fraction) due to glcc land use data
  REAL (KIND=wp), INTENT(IN)  :: for_e_glcc(:,:,:)   !< evergreen forest (fraction) due to glcc land use data
  REAL (KIND=wp), INTENT(IN)  :: emissivity_glcc(:,:,:) !< longwave emissivity due to glcc land use data


  ! local variables

  INTEGER :: ndims  

  INTEGER :: ncid
  INTEGER :: varid
  INTEGER (KIND=i8) :: undefined_i

  TYPE(dim_meta_info), ALLOCATABLE :: dim_list(:) !< dimensions for netcdf file
  TYPE(dim_meta_info), TARGET :: dim_nclass(1:3)

  INTEGER, PARAMETER :: nglob_atts=6
  TYPE(netcdf_attributes) :: global_attributes(nglob_atts)

  INTEGER :: errorcode !< error status variable

  CHARACTER (len=80):: grid_mapping !< netcdf attribute grid mapping
  CHARACTER (len=80):: coordinates  !< netcdf attribute coordinates

  INTEGER :: n !< counter

  PRINT *,'Enter write_netcdf_cosmo_grid_glcc'

  !-------------------------------------------------------------
  ! define global attributes
  CALL set_global_att_glcc(global_attributes)

  !set up dimensions for buffer
  CALL  def_dimension_info_buffer(tg)
  ! dim_3d_tg

   !set up dimensions for COSMO grid
  CALL def_dimension_info_cosmo(cosmo_grid)
  ! dim_rlon_cosmo, dim_rlat_cosmo, dim_2d_cosmo, rlon_meta, rlat_meta

  ! set mapping parameters for netcdf
  grid_mapping="rotated_pole"
  coordinates="lon lat"
  CALL set_nc_grid_def_cosmo(cosmo_grid,grid_mapping)
  ! nc_grid_def_cosmo

  ! define meta information for various land use related variables (GLCC) for netcdf output
  CALL def_glcc_fields_meta(tg,nclass_glcc,dim_2d_cosmo,coordinates,grid_mapping)
  ! dim_glcc_tg
  ! fr_land_glcc_meta, glcc_tot_npixel_meta, &
  !  &       glcc_class_fraction_meta, glcc_class_npixel_meta, &
  !  &       ice_glcc_meta, z0_glcc_meta, &
  !  &       plcov_mx_glcc_meta, plcov_mn_glcc_meta, &
  !  &       lai_mx_glcc_meta, lai_mn_glcc_meta, &
  !  &       rs_min_glcc_meta, urban_glcc_meta, &
  !  &       for_d_glcc_meta, for_e_glcc_meta, &
  !  &       emissivity_glcc_meta, root_glcc_meta

  ! define meta information for target field variables lon_geo, lat_geo 
  CALL def_com_target_fields_meta(dim_2d_cosmo,coordinates,grid_mapping)
  ! lon_geo_meta and lat_geo_meta

  !set up dimensions for buffer netcdf output 
  ndims = 3
  ALLOCATE(dim_list(1:ndims),STAT=errorcode)
  IF (errorcode /= 0 ) CALL abort_extpar('Cant allocate array dim_list')

  dim_list(1) = dim_rlon_cosmo(1) ! rlon
  dim_list(2) = dim_rlat_cosmo(1) ! rlat
  dim_list(3)%dimname = 'nclass'
  dim_list(3)%dimsize = nclass_glcc

  dim_nclass = dim_list

  undefined_i = undef_int

  !-----------------------------------------------------------------
   CALL open_new_netcdf_file(netcdf_filename=TRIM(netcdf_filename),   &
      &                       dim_list=dim_list,                  &
      &                       global_attributes=global_attributes, &
      &                       ncid=ncid)
    !-----------------------------------------------------------------

    ! rlon
    CALL netcdf_put_var(ncid,lon_rot(1:cosmo_grid%nlon_rot),rlon_meta,undefined)

    CALL netcdf_put_var(ncid,lat_rot(1:cosmo_grid%nlat_rot),rlat_meta,undefined)

    ! lon
  CALL netcdf_put_var(ncid,lon_geo(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1), &
    &                 lon_geo_meta,undefined)

  ! lat
  CALL netcdf_put_var(ncid,lat_geo(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1), &
    &                 lat_geo_meta,undefined)

  ! fr_land_glcc
  CALL netcdf_put_var(ncid,fr_land_glcc(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1), &
    &                 fr_land_glcc_meta,undefined)

  ! ice_glcc
  CALL netcdf_put_var(ncid,ice_glcc(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1), &
    &                 ice_glcc_meta,undefined)

  ! plcov_mx_glcc
  CALL netcdf_put_var(ncid,plcov_mx_glcc(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1), &
    &                 plcov_mx_glcc_meta,undefined)

  ! lai_mx_glcc
  CALL netcdf_put_var(ncid,lai_mx_glcc(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1), &
    &                 lai_mx_glcc_meta,undefined)

  ! rs_min_glcc
  CALL netcdf_put_var(ncid,rs_min_glcc(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1), &
    &                 rs_min_glcc_meta,undefined)

  ! urban_glcc
  CALL netcdf_put_var(ncid,urban_glcc(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1), &
    &                 urban_glcc_meta,undefined)

  ! for_d_glcc
  CALL netcdf_put_var(ncid,for_d_glcc(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1), &
    &                 for_d_glcc_meta,undefined)

  ! for_e_glcc
  CALL netcdf_put_var(ncid,for_e_glcc(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1), &
    &                 for_e_glcc_meta,undefined)

  ! emissivity_glcc
  CALL netcdf_put_var(ncid,emissivity_glcc(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1), &
    &                 emissivity_glcc_meta,undefined)

  ! root_glcc
  CALL netcdf_put_var(ncid,root_glcc(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1), &
    &                 root_glcc_meta,undefined)

  ! z0_glcc
  CALL netcdf_put_var(ncid,z0_glcc(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1), &
    &                 z0_glcc_meta,undefined)

  ! lai_mn_glcc
  CALL netcdf_put_var(ncid,lai_mn_glcc(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1), &
    &                 lai_mn_glcc_meta,undefined)
    
  ! plcov_mn_glcc
  CALL netcdf_put_var(ncid,plcov_mn_glcc(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1), &
    &                 plcov_mn_glcc_meta,undefined)

  ! glcc_tot_npixel
  CALL netcdf_put_var(ncid,glcc_tot_npixel(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1), &
    &                 glcc_tot_npixel_meta,undefined_i)

  ! glcc_class_fraction
  CALL netcdf_put_var(ncid,glcc_class_fraction(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1,1:nclass_glcc), &
    &                 glcc_class_fraction_meta,undefined)
  !-----------------------------------------------------------------

  ! glcc_class_npixel
  CALL netcdf_put_var(ncid,glcc_class_npixel(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1,1:nclass_glcc), &
    &                 glcc_class_npixel_meta,undefined_i)
  !-----------------------------------------------------------------
    CALL netcdf_def_grid_mapping(ncid, nc_grid_def_cosmo, varid)

    CALL close_netcdf_file(ncid)


  END SUBROUTINE write_netcdf_cosmo_grid_glcc
  !-----------------------------------------------------------------------

  !> netcdf output of GLCC derived ICON fields
  SUBROUTINE write_netcdf_icon_grid_glcc(netcdf_filename,  &
    &                                     icon_grid,       &
    &                                     tg,         &
    &                                     undefined, &
    &                                     undef_int,   &
    &                                     lon_geo,     &
    &                                     lat_geo, &
    &                                     fr_land_glcc, &
    &                                     glcc_class_fraction,    &
    &                                     glcc_class_npixel, &
    &                                     glcc_tot_npixel, &
    &                                     ice_glcc, &
    &                                     z0_glcc, &
    &                                     root_glcc, &
    &                                     plcov_mn_glcc, &
    &                                     plcov_mx_glcc, &
    &                                     lai_mn_glcc, &
    &                                     lai_mx_glcc, &
    &                                     rs_min_glcc, &
    &                                     urban_glcc,  &
    &                                     for_d_glcc,  &
    &                                     for_e_glcc, &
    &                                     emissivity_glcc)


  USE mo_var_meta_data, ONLY: dim_3d_tg, &
    &                         def_dimension_info_buffer


  USE mo_var_meta_data, ONLY: lon_geo_meta, &
    &                         lat_geo_meta, &
    &                         no_raw_data_pixel_meta, &
    &                         def_com_target_fields_meta  
   

  USE mo_var_meta_data, ONLY:  dim_icon, &
    &                          def_dimension_info_icon

  USE mo_var_meta_data, ONLY: nc_grid_def_icon, &
    &                         set_nc_grid_def_icon

  USE mo_var_meta_data, ONLY: def_glcc_fields_meta

  USE mo_var_meta_data, ONLY: dim_glcc_tg

  USE mo_var_meta_data, ONLY: fr_land_glcc_meta, glcc_tot_npixel_meta, &
    &       glcc_class_fraction_meta, glcc_class_npixel_meta, &
    &       ice_glcc_meta, z0_glcc_meta, &
    &       plcov_mx_glcc_meta, plcov_mn_glcc_meta, &
    &       lai_mx_glcc_meta, lai_mn_glcc_meta, &
    &       rs_min_glcc_meta, urban_glcc_meta, &
    &       for_d_glcc_meta, for_e_glcc_meta, &
    &       emissivity_glcc_meta, root_glcc_meta

  CHARACTER (len=*), INTENT(IN)      :: netcdf_filename !< filename for the netcdf file
  TYPE(icosahedral_triangular_grid), INTENT(IN) :: icon_grid !< structure which contains the definition of the ICON grid
  TYPE(target_grid_def), INTENT(IN) :: tg !< structure with target grid description
  REAL(KIND=wp), INTENT(IN)          :: undefined       !< value to indicate undefined grid elements 
  INTEGER, INTENT(IN)                :: undef_int       !< value to indicate undefined grid elements
  REAL (KIND=wp), INTENT(IN) :: lon_geo(:,:,:)  !< longitude coordinates of the target grid in the geographical system
  REAL (KIND=wp), INTENT(IN) :: lat_geo(:,:,:)  !< latitude coordinates of the target grid in the geographical system
  REAL (KIND=wp), INTENT(IN)  :: glcc_class_fraction(:,:,:,:)  
                                 !< fraction for each glcc class on target grid (dim (ie,je,ke,nclass_glcc))
  INTEGER (KIND=i8), INTENT(IN) :: glcc_class_npixel(:,:,:,:) 
                             !< number of raw data pixels for each glcc class on target grid (dim (ie,je,ke,nclass_glcc))
  INTEGER (KIND=i8), INTENT(IN) :: glcc_tot_npixel(:,:,:)  
                                   !< total number of glcc raw data pixels on target grid (dim (ie,je,ke))
  REAL (KIND=wp), INTENT(IN)  :: fr_land_glcc(:,:,:) !< fraction land due to glcc raw data
  REAL (KIND=wp), INTENT(IN)  :: ice_glcc(:,:,:)     !< fraction of ice due to glcc raw data
  REAL (KIND=wp), INTENT(IN)  :: z0_glcc(:,:,:)      !< roughness length due to glcc land use data
  REAL (KIND=wp), INTENT(IN)  :: root_glcc(:,:,:)    !< root depth due to glcc land use data
  REAL (KIND=wp), INTENT(IN)  :: plcov_mx_glcc(:,:,:)!< plant cover maximum due to glcc land use data
  REAL (KIND=wp), INTENT(IN)  :: plcov_mn_glcc(:,:,:)!< plant cover minimum due to glcc land use data
  REAL (KIND=wp), INTENT(IN)  :: lai_mx_glcc(:,:,:)  !< Leaf Area Index maximum due to glcc land use data
  REAL (KIND=wp), INTENT(IN)  :: lai_mn_glcc(:,:,:)  !< Leaf Area Index minimum due to glcc land use data
  REAL (KIND=wp), INTENT(IN)  :: rs_min_glcc(:,:,:)  !< minimal stomata resistance due to glcc land use data
  REAL (KIND=wp), INTENT(IN)  :: urban_glcc(:,:,:)   !< urban fraction due to glcc land use data
  REAL (KIND=wp), INTENT(IN)  :: for_d_glcc(:,:,:)   !< deciduous forest (fraction) due to glcc land use data
  REAL (KIND=wp), INTENT(IN)  :: for_e_glcc(:,:,:)   !< evergreen forest (fraction) due to glcc land use data
  REAL (KIND=wp), INTENT(IN)  :: emissivity_glcc(:,:,:) !< longwave emissivity due to glcc land use data


  ! local variables


  INTEGER :: ndims 
  INTEGER :: ncid
  INTEGER :: varid
  INTEGER (KIND=i8) :: undefined_i

  TYPE(dim_meta_info), ALLOCATABLE :: dim_list(:) !< dimensions for netcdf file
  TYPE(dim_meta_info), TARGET :: dim_2d_icon(1:2)
  TYPE(dim_meta_info), TARGET :: dim_1d_icon(1:1)

  INTEGER, PARAMETER :: nglob_atts=6
  TYPE(netcdf_attributes) :: global_attributes(nglob_atts)

  INTEGER :: errorcode !< error status variable

  CHARACTER (len=80):: grid_mapping !< netcdf attribute grid mapping
  CHARACTER (len=80):: coordinates  !< netcdf attribute coordinates

  INTEGER :: n !< counter

  !-------------------------------------------------------------
  ! define global attributes
  CALL set_global_att_glcc(global_attributes)

  !set up dimensions for buffer
  CALL  def_dimension_info_buffer(tg)
  ! dim_3d_tg

  !set up dimensions for ICON grid
  CALL def_dimension_info_icon(icon_grid)
  ! dim_icon

  ! define meta information for various land use related variables (GLCC) for netcdf output
  CALL def_glcc_fields_meta(tg,nclass_glcc,dim_icon)
  ! dim_glcc_tg
  ! fr_land_glcc_meta, glcc_tot_npixel_meta, &
  !  &       glcc_class_fraction_meta, glcc_class_npixel_meta, &
  !  &       ice_glcc_meta, z0_glcc_meta, &
  !  &       plcov_mx_glcc_meta, plcov_mn_glcc_meta, &
  !  &       lai_mx_glcc_meta, lai_mn_glcc_meta, &
  !  &       rs_min_glcc_meta, urban_glcc_meta, &
  !  &       for_d_glcc_meta, for_e_glcc_meta, &
  !  &       emissivity_glcc_meta, root_glcc_meta

  ! define meta information for target field variables lon_geo, lat_geo 
  CALL def_com_target_fields_meta(dim_icon)
  ! lon_geo_meta and lat_geo_meta

  ! set mapping parameters for netcdf
  grid_mapping="lon_lat_on_sphere"
  coordinates="lon lat"

  CALL set_nc_grid_def_icon(grid_mapping)
  ! nc_grid_def_icon


  !set up dimensions for buffer netcdf output 
  ndims = 2
  ALLOCATE(dim_list(1:ndims),STAT=errorcode)
  IF (errorcode /= 0 ) CALL abort_extpar('Cant allocate array dim_list')

  dim_list(1) = dim_icon(1) ! cell
  dim_list(2)%dimname = 'nclass'
  dim_list(2)%dimsize = nclass_glcc

  dim_1d_icon =  dim_icon(1) ! cell
  dim_2d_icon = dim_list

  undefined_i = undef_int

  !-----------------------------------------------------------------
   CALL open_new_netcdf_file(netcdf_filename=TRIM(netcdf_filename),   &
      &                       dim_list=dim_list,                  &
      &                       global_attributes=global_attributes, &
      &                       ncid=ncid)
  !-----------------------------------------------------------------

   ! lon
    CALL netcdf_put_var(ncid,lon_geo(1:icon_grid%ncell,1,1),lon_geo_meta,undefined)

    ! lat
    CALL netcdf_put_var(ncid,lat_geo(1:icon_grid%ncell,1,1),lat_geo_meta,undefined)

    !-----------------------------------------------------------------
    n=1 ! fr_land_glcc
    CALL netcdf_put_var(ncid,fr_land_glcc(1:icon_grid%ncell,1,1),fr_land_glcc_meta,undefined)

    n=2 ! ice_glcc
    CALL netcdf_put_var(ncid,ice_glcc(1:icon_grid%ncell,1,1),ice_glcc_meta,undefined)

    n=3 ! plcov_mx_glcc
    CALL netcdf_put_var(ncid,plcov_mx_glcc(1:icon_grid%ncell,1,1),plcov_mx_glcc_meta,undefined)

    n=4 ! lai_mx_glcc
    CALL netcdf_put_var(ncid,lai_mx_glcc(1:icon_grid%ncell,1,1),lai_mx_glcc_meta,undefined)

    n=5 ! rs_min_glcc
    CALL netcdf_put_var(ncid,rs_min_glcc(1:icon_grid%ncell,1,1),rs_min_glcc_meta,undefined)

    n=6 ! urban_glcc
    CALL netcdf_put_var(ncid,urban_glcc(1:icon_grid%ncell,1,1),urban_glcc_meta,undefined)

    n=7 ! for_d_glcc
    CALL netcdf_put_var(ncid,for_d_glcc(1:icon_grid%ncell,1,1),for_d_glcc_meta,undefined)

    n=8 ! for_e_glcc
    CALL netcdf_put_var(ncid,for_e_glcc(1:icon_grid%ncell,1,1),for_e_glcc_meta,undefined)

    n=9 ! emissivity_glcc
    CALL netcdf_put_var(ncid, emissivity_glcc(1:icon_grid%ncell,1,1),emissivity_glcc_meta,undefined)

    n=10 ! root_glcc
    CALL netcdf_put_var(ncid,root_glcc(1:icon_grid%ncell,1,1),root_glcc_meta,undefined)

    n=11 ! z0_glcc
    CALL netcdf_put_var(ncid,z0_glcc(1:icon_grid%ncell,1,1),z0_glcc_meta,undefined)

  !-----------------------------------------------------------------
  CALL close_netcdf_file(ncid)


  END SUBROUTINE write_netcdf_icon_grid_glcc
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !> set global attributes for netcdf with glcc data
  SUBROUTINE set_global_att_glcc(global_attributes)
    USE mo_glcc_lookup_tables, ONLY: name_lookup_table_glcc
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
    global_attributes(3)%attributetext='GLCC data'

    CALL DATE_AND_TIME(ydate,ytime)
    READ(ydate,'(4A2)') cc,yy,mm,dd
    READ(ytime,'(2A2)') hh, minute

    ydate=TRIM(cc)//TRIM(yy)//'-'//TRIM(mm)//'-'//TRIM(dd)
    ytime=TRIM(hh)//':'//TRIM(minute) 

    global_attributes(4)%attname = 'history'
    global_attributes(4)%attributetext=TRIM(ydate)//'T'//TRIM(ytime)//' glcc_to_buffer'

    global_attributes(5)%attname = 'references'
    global_attributes(5)%attributetext=''

    global_attributes(6)%attname = 'comment'
    global_attributes(6)%attributetext='Landuse data look-up table: '//TRIM(name_lookup_table_glcc)


  END SUBROUTINE set_global_att_glcc
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------



  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------

  !> read GLC2000 derived buffer fields
    SUBROUTINE read_netcdf_buffer_glc2000(netcdf_filename,  &
    &                                     tg,         &
    &                                     undefined, &
    &                                     undef_int,   &
    &                                     fr_land_glc2000, &
    &                                     glc2000_class_fraction,    &
    &                                     glc2000_class_npixel, &
    &                                     glc2000_tot_npixel, &
    &                                     ice_glc2000, &
    &                                     z0_glc2000, &
    &                                     root_glc2000, &
    &                                     plcov_mn_glc2000, &
    &                                     plcov_mx_glc2000, &
    &                                     lai_mn_glc2000, &
    &                                     lai_mx_glc2000, &
    &                                     rs_min_glc2000, &
    &                                     urban_glc2000,  &
    &                                     for_d_glc2000,  &
    &                                     for_e_glc2000, &
    &                                     emissivity_glc2000)


  USE mo_var_meta_data, ONLY: dim_3d_tg, &
    &                         def_dimension_info_buffer


  USE mo_var_meta_data, ONLY: lon_geo_meta, &
    &                         lat_geo_meta, &
    &                         no_raw_data_pixel_meta, &
    &                         def_com_target_fields_meta  
  
  USE mo_var_meta_data, ONLY: def_glc2000_fields_meta

  USE mo_var_meta_data, ONLY: dim_glc2000_tg

  USE mo_var_meta_data, ONLY: fr_land_glc2000_meta, glc2000_tot_npixel_meta, &
    &       glc2000_class_fraction_meta, glc2000_class_npixel_meta, &
    &       ice_glc2000_meta, z0_glc2000_meta, &
    &       plcov_mx_glc2000_meta, plcov_mn_glc2000_meta, &
    &       lai_mx_glc2000_meta, lai_mn_glc2000_meta, &
    &       rs_min_glc2000_meta, urban_glc2000_meta, &
    &       for_d_glc2000_meta, for_e_glc2000_meta, &
    &       emissivity_glc2000_meta, root_glc2000_meta

  !USE mo_io_utilities, ONLY: netcdf_get_var_real_3d, netcdf_get_var_real_4d
  !USE mo_io_utilities, ONLY: netcdf_get_var_int_3d, netcdf_get_var_int_4d

  USE mo_io_utilities, ONLY: netcdf_get_var

  CHARACTER (len=*), INTENT(IN)      :: netcdf_filename !< filename for the netcdf file
  TYPE(target_grid_def), INTENT(IN) :: tg !< structure with target grid description
  REAL(KIND=wp), INTENT(OUT)          :: undefined       !< value to indicate undefined grid elements 
  INTEGER, INTENT(OUT)                :: undef_int       !< value to indicate undefined grid elements
  REAL (KIND=wp), INTENT(OUT)  :: glc2000_class_fraction(:,:,:,:)  
                                  !< fraction for each glc2000 class on target grid (dim (ie,je,ke,nclass_glc2000))
  INTEGER (KIND=i8), INTENT(OUT) :: glc2000_class_npixel(:,:,:,:) 
                            !< number of raw data pixels for each glc2000 class on target grid (dim (ie,je,ke,nclass_glc2000))
  INTEGER (KIND=i8), INTENT(OUT) :: glc2000_tot_npixel(:,:,:)
                                    !< total number of glc2000 raw data pixels on target grid (dim (ie,je,ke))
  REAL (KIND=wp), INTENT(OUT)  :: fr_land_glc2000(:,:,:) !< fraction land due to glc2000 raw data
  REAL (KIND=wp), INTENT(OUT)  :: ice_glc2000(:,:,:)     !< fraction of ice due to glc2000 raw data
  REAL (KIND=wp), INTENT(OUT)  :: z0_glc2000(:,:,:)      !< roughness length due to glc2000 land use data
  REAL (KIND=wp), INTENT(OUT)  :: root_glc2000(:,:,:)    !< root depth due to glc2000 land use data
  REAL (KIND=wp), INTENT(OUT)  :: plcov_mn_glc2000(:,:,:)!< plant cover maximum due to glc2000 land use data
  REAL (KIND=wp), INTENT(OUT)  :: plcov_mx_glc2000(:,:,:)!< plant cover minimum due to glc2000 land use data
  REAL (KIND=wp), INTENT(OUT)  :: lai_mn_glc2000(:,:,:)  !< Leaf Area Index maximum due to glc2000 land use data
  REAL (KIND=wp), INTENT(OUT)  :: lai_mx_glc2000(:,:,:)  !< Leaf Area Index minimum due to glc2000 land use data
  REAL (KIND=wp), INTENT(OUT)  :: rs_min_glc2000(:,:,:)  !< minimal stomata resistance due to glc2000 land use data
  REAL (KIND=wp), INTENT(OUT)  :: urban_glc2000(:,:,:)   !< urban fraction due to glc2000 land use data
  REAL (KIND=wp), INTENT(OUT)  :: for_d_glc2000(:,:,:)   !< deciduous forest (fraction) due to glc2000 land use data
  REAL (KIND=wp), INTENT(OUT)  :: for_e_glc2000(:,:,:)   !< evergreen forest (fraction) due to glc2000 land use data
  REAL (KIND=wp), INTENT(OUT)  :: emissivity_glc2000(:,:,:) !< longwave emissivity due to glc2000 land use data


  ! local variables
  INTEGER :: errorcode !< error status variable
  INTEGER :: n !< counter

  PRINT *,'ENTER read_netcdf_buffer_glc2000'


  !set up dimensions for buffer
  CALL  def_dimension_info_buffer(tg)
  ! dim_3d_tg

  ! define meta information for various land use related variables (GLC2000) for netcdf output
  CALL def_glc2000_fields_meta(tg,nclass_glc2000,dim_3d_tg)
  ! dim_glc2000_tg
  ! fr_land_glc2000_meta, glc2000_tot_npixel_meta, &
  !  &       glc2000_class_fraction_meta, glc2000_class_npixel_meta, &
  !  &       ice_glc2000_meta, z0_glc2000_meta, &
  !  &       plcov_mx_glc2000_meta, plcov_mn_glc2000_meta, &
  !  &       lai_mx_glc2000_meta, lai_mn_glc2000_meta, &
  !  &       rs_min_glc2000_meta, urban_glc2000_meta, &
  !  &       for_d_glc2000_meta, for_e_glc2000_meta, &
  !  &       emissivity_glc2000_meta, root_glc2000_meta

  ! define meta information for target field variables lon_geo, lat_geo 
  CALL def_com_target_fields_meta(dim_3d_tg)
  ! lon_geo_meta and lat_geo_meta
  

  PRINT *,'CALL read netcdf data Land Use'

  CALL netcdf_get_var(TRIM(netcdf_filename),fr_land_glc2000_meta,fr_land_glc2000)
  PRINT *,'fr_land_glc2000 read'

  CALL netcdf_get_var(TRIM(netcdf_filename),glc2000_tot_npixel_meta,glc2000_tot_npixel)
  PRINT *,'glc2000_tot_npixel read'

  CALL netcdf_get_var(TRIM(netcdf_filename),glc2000_class_fraction_meta,glc2000_class_fraction)
  PRINT *,'glc2000_class_fraction read'

  CALL netcdf_get_var(TRIM(netcdf_filename),glc2000_class_npixel_meta,glc2000_class_npixel)
  PRINT *,'glc2000_class_npixel read'

  CALL netcdf_get_var(TRIM(netcdf_filename),ice_glc2000_meta,ice_glc2000)
  PRINT *,'ice_glc2000 read'

  CALL netcdf_get_var(TRIM(netcdf_filename),z0_glc2000_meta,z0_glc2000)
  PRINT *,'z0_glc2000 read'

  CALL netcdf_get_var(TRIM(netcdf_filename),plcov_mx_glc2000_meta,plcov_mx_glc2000)
  PRINT *,'plcov_mx_glc2000 read'

  CALL netcdf_get_var(TRIM(netcdf_filename),plcov_mn_glc2000_meta,plcov_mn_glc2000)
  PRINT *,'plcov_mn_glc2000 read'

  CALL netcdf_get_var(TRIM(netcdf_filename),lai_mx_glc2000_meta,lai_mx_glc2000)
  PRINT *,'lai_mx_glc2000 read'

  CALL netcdf_get_var(TRIM(netcdf_filename),lai_mn_glc2000_meta,lai_mn_glc2000)
  PRINT *,'lai_mn_glc2000 read'

  CALL netcdf_get_var(TRIM(netcdf_filename),rs_min_glc2000_meta,rs_min_glc2000)
  PRINT *,'rs_min_glc2000 read'

  CALL netcdf_get_var(TRIM(netcdf_filename),urban_glc2000_meta,urban_glc2000)
  PRINT *,'urban_glc2000 read'

  CALL netcdf_get_var(TRIM(netcdf_filename),for_d_glc2000_meta,for_d_glc2000)
  PRINT *,'for_d_glc2000 read'

  CALL netcdf_get_var(TRIM(netcdf_filename),for_e_glc2000_meta,for_e_glc2000)
  PRINT *,'for_e_glc2000 read'

  CALL netcdf_get_var(TRIM(netcdf_filename),emissivity_glc2000_meta,emissivity_glc2000)
  PRINT *,'emissivity_glc2000 read'

  CALL netcdf_get_var(TRIM(netcdf_filename),root_glc2000_meta,root_glc2000)
  PRINT *,'root_glc2000 read'


  END SUBROUTINE read_netcdf_buffer_glc2000
  !-----------------------------------------------------------------------

  
  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------

  !> read GLCC derived buffer fields
    SUBROUTINE read_netcdf_buffer_glcc(netcdf_filename,  &
    &                                     tg,         &
    &                                     undefined, &
    &                                     undef_int,   &
    &                                     fr_land_glcc, &
    &                                     glcc_class_fraction,    &
    &                                     glcc_class_npixel, &
    &                                     glcc_tot_npixel, &
    &                                     ice_glcc, &
    &                                     z0_glcc, &
    &                                     root_glcc, &
    &                                     plcov_mn_glcc, &
    &                                     plcov_mx_glcc, &
    &                                     lai_mn_glcc, &
    &                                     lai_mx_glcc, &
    &                                     rs_min_glcc, &
    &                                     urban_glcc,  &
    &                                     for_d_glcc,  &
    &                                     for_e_glcc, &
    &                                     emissivity_glcc)


  USE mo_var_meta_data, ONLY: dim_3d_tg, &
    &                         def_dimension_info_buffer


  USE mo_var_meta_data, ONLY: lon_geo_meta, &
    &                         lat_geo_meta, &
    &                         no_raw_data_pixel_meta, &
    &                         def_com_target_fields_meta  
  
  USE mo_var_meta_data, ONLY: def_glcc_fields_meta

  USE mo_var_meta_data, ONLY: dim_glcc_tg

  USE mo_var_meta_data, ONLY: fr_land_glcc_meta, glcc_tot_npixel_meta, &
    &       glcc_class_fraction_meta, glcc_class_npixel_meta, &
    &       ice_glcc_meta, z0_glcc_meta, &
    &       plcov_mx_glcc_meta, plcov_mn_glcc_meta, &
    &       lai_mx_glcc_meta, lai_mn_glcc_meta, &
    &       rs_min_glcc_meta, urban_glcc_meta, &
    &       for_d_glcc_meta, for_e_glcc_meta, &
    &       emissivity_glcc_meta, root_glcc_meta

  !USE mo_io_utilities, ONLY: netcdf_get_var_real_3d, netcdf_get_var_real_4d
  !USE mo_io_utilities, ONLY: netcdf_get_var_int_3d, netcdf_get_var_int_4d

  USE mo_io_utilities, ONLY: netcdf_get_var

  CHARACTER (len=*), INTENT(IN)      :: netcdf_filename !< filename for the netcdf file
  TYPE(target_grid_def), INTENT(IN) :: tg !< structure with target grid description
  REAL(KIND=wp), INTENT(OUT)          :: undefined       !< value to indicate undefined grid elements 
  INTEGER, INTENT(OUT)                :: undef_int       !< value to indicate undefined grid elements
  REAL (KIND=wp), INTENT(OUT)  :: glcc_class_fraction(:,:,:,:)  
                                  !< fraction for each glcc class on target grid (dim (ie,je,ke,nclass_glcc))
  INTEGER (KIND=i8), INTENT(OUT) :: glcc_class_npixel(:,:,:,:) 
                           !< number of raw data pixels for each glcc class on target grid (dim (ie,je,ke,nclass_glcc))
  INTEGER (KIND=i8), INTENT(OUT) :: glcc_tot_npixel(:,:,:)  
                                    !< total number of glcc raw data pixels on target grid (dim (ie,je,ke))
  REAL (KIND=wp), INTENT(OUT)  :: fr_land_glcc(:,:,:) !< fraction land due to glcc raw data
  REAL (KIND=wp), INTENT(OUT)  :: ice_glcc(:,:,:)     !< fraction of ice due to glcc raw data
  REAL (KIND=wp), INTENT(OUT)  :: z0_glcc(:,:,:)      !< roughness length due to glcc land use data
  REAL (KIND=wp), INTENT(OUT)  :: root_glcc(:,:,:)    !< root depth due to glcc land use data
  REAL (KIND=wp), INTENT(OUT)  :: plcov_mn_glcc(:,:,:)!< plant cover maximum due to glcc land use data
  REAL (KIND=wp), INTENT(OUT)  :: plcov_mx_glcc(:,:,:)!< plant cover minimum due to glcc land use data
  REAL (KIND=wp), INTENT(OUT)  :: lai_mn_glcc(:,:,:)  !< Leaf Area Index maximum due to glcc land use data
  REAL (KIND=wp), INTENT(OUT)  :: lai_mx_glcc(:,:,:)  !< Leaf Area Index minimum due to glcc land use data
  REAL (KIND=wp), INTENT(OUT)  :: rs_min_glcc(:,:,:)  !< minimal stomata resistance due to glcc land use data
  REAL (KIND=wp), INTENT(OUT)  :: urban_glcc(:,:,:)   !< urban fraction due to glcc land use data
  REAL (KIND=wp), INTENT(OUT)  :: for_d_glcc(:,:,:)   !< deciduous forest (fraction) due to glcc land use data
  REAL (KIND=wp), INTENT(OUT)  :: for_e_glcc(:,:,:)   !< evergreen forest (fraction) due to glcc land use data
  REAL (KIND=wp), INTENT(OUT)  :: emissivity_glcc(:,:,:) !< longwave emissivity due to glcc land use data


  ! local variables
  INTEGER :: errorcode !< error status variable
  INTEGER :: n !< counter

  PRINT *,'ENTER read_netcdf_buffer_glcc'


  !set up dimensions for buffer
  CALL  def_dimension_info_buffer(tg)
  ! dim_3d_tg

  ! define meta information for various land use related variables (GLCC) for netcdf output
  CALL def_glcc_fields_meta(tg,nclass_glcc,dim_3d_tg)
  ! dim_glcc_tg
  ! fr_land_glcc_meta, glcc_tot_npixel_meta, &
  !  &       glcc_class_fraction_meta, glcc_class_npixel_meta, &
  !  &       ice_glcc_meta, z0_glcc_meta, &
  !  &       plcov_mx_glcc_meta, plcov_mn_glcc_meta, &
  !  &       lai_mx_glcc_meta, lai_mn_glcc_meta, &
  !  &       rs_min_glcc_meta, urban_glcc_meta, &
  !  &       for_d_glcc_meta, for_e_glcc_meta, &
  !  &       emissivity_glcc_meta, root_glcc_meta

  ! define meta information for target field variables lon_geo, lat_geo 
  CALL def_com_target_fields_meta(dim_3d_tg)
  ! lon_geo_meta and lat_geo_meta
  

  PRINT *,'CALL read netcdf data Land Use'

  CALL netcdf_get_var(TRIM(netcdf_filename),fr_land_glcc_meta,fr_land_glcc)
  PRINT *,'fr_land_glcc read'

  CALL netcdf_get_var(TRIM(netcdf_filename),glcc_tot_npixel_meta,glcc_tot_npixel)
  PRINT *,'glcc_tot_npixel read'

  CALL netcdf_get_var(TRIM(netcdf_filename),glcc_class_fraction_meta,glcc_class_fraction)
  PRINT *,'glcc_class_fraction read'

  CALL netcdf_get_var(TRIM(netcdf_filename),glcc_class_npixel_meta,glcc_class_npixel)
  PRINT *,'glcc_class_npixel read'

  CALL netcdf_get_var(TRIM(netcdf_filename),ice_glcc_meta,ice_glcc)
  PRINT *,'ice_glcc read'

  CALL netcdf_get_var(TRIM(netcdf_filename),z0_glcc_meta,z0_glcc)
  PRINT *,'z0_glcc read'

  CALL netcdf_get_var(TRIM(netcdf_filename),plcov_mx_glcc_meta,plcov_mx_glcc)
  PRINT *,'plcov_mx_glcc read'

  CALL netcdf_get_var(TRIM(netcdf_filename),plcov_mn_glcc_meta,plcov_mn_glcc)
  PRINT *,'plcov_mn_glcc read'

  CALL netcdf_get_var(TRIM(netcdf_filename),lai_mx_glcc_meta,lai_mx_glcc)
  PRINT *,'lai_mx_glcc read'

  CALL netcdf_get_var(TRIM(netcdf_filename),lai_mn_glcc_meta,lai_mn_glcc)
  PRINT *,'lai_mn_glcc read'

  CALL netcdf_get_var(TRIM(netcdf_filename),rs_min_glcc_meta,rs_min_glcc)
  PRINT *,'rs_min_glcc read'

  CALL netcdf_get_var(TRIM(netcdf_filename),urban_glcc_meta,urban_glcc)
  PRINT *,'urban_glcc read'

  CALL netcdf_get_var(TRIM(netcdf_filename),for_d_glcc_meta,for_d_glcc)
  PRINT *,'for_d_glcc read'

  CALL netcdf_get_var(TRIM(netcdf_filename),for_e_glcc_meta,for_e_glcc)
  PRINT *,'for_e_glcc read'

  CALL netcdf_get_var(TRIM(netcdf_filename),emissivity_glcc_meta,emissivity_glcc)
  PRINT *,'emissivity_glcc read'

  CALL netcdf_get_var(TRIM(netcdf_filename),root_glcc_meta,root_glcc)
  PRINT *,'root_glcc read'


  END SUBROUTINE read_netcdf_buffer_glcc
  !-----------------------------------------------------------------------


    
 
END Module mo_landuse_output_nc

