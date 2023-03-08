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

  USE mo_logging
  USE mo_kind,                  ONLY: wp, i4

  USE mo_io_units,              ONLY: filename_max

  USE mo_grid_structures,       ONLY: igrid_icon, target_grid_def

  USE mo_io_utilities,          ONLY: netcdf_attributes, &
       &                              dim_meta_info, &
       &                              netcdf_put_var, &
       &                              open_new_netcdf_file, &
       &                              close_netcdf_file, &
       &                              netcdf_get_var, &
       &                              set_date_mm_extpar_field, &
       &                              my_get_var_int_4d

  USE mo_glc2000_lookup_tables, ONLY: nclass_glc2000, &
       &                              get_name_glc2000_lookup_tables, &
       &                              name_lookup_table_glc2000

  USE mo_glcc_lookup_tables,    ONLY: nclass_glcc, &
       &                              name_lookup_table_glcc, &
       &                              get_name_glcc_lookup_tables

  USE mo_ecosg_lookup_tables,   ONLY: nclass_ecosg, &
       &                              name_lookup_table_ecosg, &
       &                              get_name_ecosg_lookup_tables

  USE mo_var_meta_data,         ONLY: dim_3d_tg, &
       &                              def_dimension_info_buffer, &
       &                              lon_geo_meta, &
       &                              lat_geo_meta, &
       &                              def_com_target_fields_meta, &
       &                              def_lu_fields_meta, &
       &                              dim_lu_tg, &
       &                              fr_land_lu_meta, lu_tot_npixel_meta, &
       &                              lu_class_fraction_meta, lu_class_npixel_meta, &
       &                              ice_lu_meta, z0_lu_meta, &
       &                              plcov_mx_lu_meta, plcov_mn_lu_meta, &
       &                              lai_mx_lu_meta, lai_mn_lu_meta, &
       &                              rs_min_lu_meta, urban_lu_meta, &
       &                              for_d_lu_meta, for_e_lu_meta, &
       &                              skinc_lu_meta, &
       &                              emissivity_lu_meta, root_lu_meta, &
       &                              def_glc2000_fields_meta, &
       &                              dim_glc2000_tg, &
       &                              fr_land_glc2000_meta, glc2000_tot_npixel_meta, &
       &                              glc2000_class_fraction_meta, glc2000_class_npixel_meta, &
       &                              ice_glc2000_meta, z0_glc2000_meta, &
       &                              plcov_mx_glc2000_meta, plcov_mn_glc2000_meta, &
       &                              lai_mx_glc2000_meta, lai_mn_glc2000_meta, &
       &                              rs_min_glc2000_meta, urban_glc2000_meta, &
       &                              for_d_glc2000_meta, for_e_glc2000_meta, &
       &                              emissivity_glc2000_meta, root_glc2000_meta, &
       &                              def_glcc_fields_meta, &
       &                              fr_land_glcc_meta, glcc_tot_npixel_meta, &
       &                              glcc_class_fraction_meta, glcc_class_npixel_meta, &
       &                              ice_glcc_meta, z0_glcc_meta, &
       &                              plcov_mx_glcc_meta, plcov_mn_glcc_meta, &
       &                              lai_mx_glcc_meta, lai_mn_glcc_meta, &
       &                              rs_min_glcc_meta, urban_glcc_meta, &
       &                              for_d_glcc_meta, for_e_glcc_meta, &
       &                              emissivity_glcc_meta, root_glcc_meta, &
       &                              dim_glcc_tg

  USE mo_lu_tg_fields,        ONLY :  i_lu_globcover, i_lu_glc2000, i_lu_glcc, &
       &                              i_lu_ecosg

  USE mo_globcover_lookup_tables, ONLY: get_name_globcover_lookup_tables

  USE mo_terra_urb,               ONLY: l_terra_urb, &
       &                                terra_urb_write_netcdf, &
       &                                terra_urb_read_netcdf

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: write_netcdf_buffer_glc2000, &
       &    write_netcdf_buffer_glcc, &
       &    read_netcdf_buffer_glcc, &
       &    write_netcdf_buffer_lu, &
       &    read_netcdf_buffer_lu

  CONTAINS

  !> netcdf output of landuse buffer fields
  SUBROUTINE write_netcdf_buffer_lu(netcdf_filename,  &
       &                            lu_dataset, &
       &                            tg,         &
       &                            i_landuse_data, &
       &                            ilookup_table_lu, &
       &                            nclass_lu, &
       &                            undefined, &
       &                            undef_int,   &
       &                            lon_geo,     &
       &                            lat_geo, &
       &                            fr_land_lu, &
       &                            lu_class_fraction,    &
       &                            lu_class_npixel, &
       &                            lu_tot_npixel, &
       &                            ice_lu, &
       &                            z0_lu, &
       &                            root_lu, &
       &                            plcov_mn_lu, &
       &                            plcov_mx_lu, &
       &                            lai_mn_lu, &
       &                            lai_mx_lu, &
       &                            rs_min_lu, &
       &                            urban_lu,  &
       &                            for_d_lu,  &
       &                            for_e_lu, &
       &                            skinc_lu, &
       &                            emissivity_lu)


    CHARACTER (len=*), INTENT(IN)     :: netcdf_filename, & !< filename for the netcdf file
         &                               lu_dataset !< name of landuse data set

    TYPE(target_grid_def), INTENT(IN) :: tg !< structure with target grid description

    INTEGER (KIND=i4), INTENT(IN)     :: lu_class_npixel(:,:,:,:), &
         &                               lu_tot_npixel(:,:,:), &
         &                               i_landuse_data, & !<integer switch to choose a land use raw data set
         &                               ilookup_table_lu, & !< integer switch to choose a lookup table
         &                               nclass_lu, & !< number of land use classes
         &                               undef_int       !< value to indicate undefined grid elements

    REAL(KIND=wp), INTENT(IN)         :: undefined, &       !< value to indicate undefined grid elements
         &                               lon_geo(:,:,:), &  !< longitude coordinates of the target grid in the geographical system
         &                               lat_geo(:,:,:), &  !< latitude coordinates of the target grid in the geographical system
         &                               lu_class_fraction(:,:,:,:), &
         &                               fr_land_lu(:,:,:), & !< fraction land due to lu raw data
         &                               ice_lu(:,:,:), &     !< fraction of ice due to lu raw data
         &                               z0_lu(:,:,:), &      !< roughness length due to lu land use data
         &                               root_lu(:,:,:), &    !< root depth due to lu land use data
         &                               plcov_mx_lu(:,:,:), &!< plant cover maximum due to lu land use data
         &                               plcov_mn_lu(:,:,:), &!< plant cover minimum due to lu land use data
         &                               lai_mx_lu(:,:,:), &  !< Leaf Area Index maximum due to lu land use data
         &                               lai_mn_lu(:,:,:), &  !< Leaf Area Index minimum due to lu land use data
         &                               rs_min_lu(:,:,:), &  !< minimal stomata resistance due to lu land use data
         &                               urban_lu(:,:,:), &   !< urban fraction due to lu land use data
         &                               for_d_lu(:,:,:), &   !< deciduous forest (fraction) due to lu land use data
         &                               for_e_lu(:,:,:), &   !< evergreen forest (fraction) due to lu land use data
         &                               skinc_lu(:,:,:), &   !< skin conductivity due to lu land use data
         &                               emissivity_lu(:,:,:) !< longwave emissivity due to lu land use data


    ! local variables
    INTEGER (KIND=i4)                 :: undefined_i, ndims, ncid, errorcode

    TYPE(dim_meta_info), ALLOCATABLE  :: dim_list(:) !< dimensions for netcdf file

    INTEGER(KIND=i4), PARAMETER       :: nglob_atts=6

    TYPE(netcdf_attributes)           :: global_attributes(nglob_atts)

    CALL logging%info('Enter routine: write_netcdf_buffer_lu')

    !-------------------------------------------------------------
    ! define global attributes
    CALL set_global_att_lu(i_landuse_data,ilookup_table_lu,global_attributes)

    !set up dimensions for buffer
    CALL  def_dimension_info_buffer(tg)
    ! dim_3d_tg

    ! define meta information for various land use related variables for netcdf output
    CALL def_lu_fields_meta(nclass_lu,dim_3d_tg,lu_dataset=lu_dataset)

    ! define meta information for target field variables lon_geo, lat_geo
    CALL def_com_target_fields_meta(dim_3d_tg)
    ! lon_geo_meta and lat_geo_meta

    !set up dimensions for buffer netcdf output
    ndims = 4
    ALLOCATE(dim_list(1:ndims),STAT=errorcode)
    IF (errorcode /= 0 ) CALL logging%error('Cant allocate array dim_list',__FILE__,__LINE__)
    dim_list = dim_lu_tg

    undefined_i = undef_int

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

    IF (l_terra_urb) CALL terra_urb_write_netcdf(ncid, undefined)

    CALL close_netcdf_file(ncid)

    CALL logging%info('Exit routine: write_netcdf_buffer_lu')

  END SUBROUTINE write_netcdf_buffer_lu

  !-----------------------------------------------------------------------
  !> set global attributes for netcdf with lu data
  SUBROUTINE set_global_att_lu(i_landuse_data,ilookup_table_lu,global_attributes)

    INTEGER(KIND=i4), INTENT(IN)           :: i_landuse_data, & !<integer switch to choose a land use raw data set
         &                                    ilookup_table_lu !< integer switch to choose a lookup table

    TYPE(netcdf_attributes), INTENT(INOUT) :: global_attributes(1:6)

    !local variables
    CHARACTER(len=10)                      :: ydate, &
         &                                    ytime

    CHARACTER(len=2)                       :: cc, &
         &                                    yy, &
         &                                    mm, &
         &                                    dd, &
         &                                    hh, &
         &                                    minute

    CHARACTER(len=filename_max)           :: name_lookup_table_lu

    ! define global attributes
    global_attributes(1)%attname = 'title'
    global_attributes(1)%attributetext='Land Use data'
    global_attributes(2)%attname = 'institution'
    global_attributes(2)%attributetext='Deutscher Wetterdienst'

    global_attributes(3)%attname = 'source'
    SELECT CASE (i_landuse_data)
      CASE (i_lu_globcover)
        global_attributes(3)%attributetext='Globcover 2009 data'
         CALL get_name_globcover_lookup_tables(ilookup_table_lu, name_lookup_table_lu)
      CASE (i_lu_glc2000)
        global_attributes(3)%attributetext='GLC2000 data'
         CALL get_name_glc2000_lookup_tables(ilookup_table_lu, name_lookup_table_lu)
      CASE(i_lu_glcc)
        global_attributes(3)%attributetext='GLCC data'
         CALL get_name_glcc_lookup_tables(ilookup_table_lu, name_lookup_table_lu)
      CASE(i_lu_ecosg)
        global_attributes(3)%attributetext='ECOCLIMAP SG'
         CALL get_name_ecosg_lookup_tables(ilookup_table_lu, name_lookup_table_lu)
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

  !> read land use derived buffer fields
  SUBROUTINE read_netcdf_buffer_lu(netcdf_filename,  &
      &                            tg,         &
      &                            nclass_lu, &
      &                            fr_land_lu, &
      &                            lu_class_fraction,    &
      &                            lu_class_npixel, &
      &                            lu_tot_npixel, &
      &                            ice_lu, &
      &                            z0_lu, &
      &                            root_lu, &
      &                            plcov_mn_lu, &
      &                            plcov_mx_lu, &
      &                            lai_mn_lu, &
      &                            lai_mx_lu, &
      &                            rs_min_lu, &
      &                            urban_lu,  &
      &                            for_d_lu,  &
      &                            for_e_lu, &
      &                            skinc_lu, &
      &                            emissivity_lu)


    CHARACTER (len=*), INTENT(IN)      :: netcdf_filename !< filename for the netcdf file

    TYPE(target_grid_def), INTENT(IN)  :: tg !< structure with target grid description

    INTEGER(KIND=i4), INTENT(IN)       :: nclass_lu !< number of land use classes


    INTEGER (KIND=i4), INTENT(OUT)     :: lu_class_npixel(:,:,:,:), &
                                          !< number of raw data pixels for each lu class on target grid (dim (ie,je,ke,nclass_lu))
         &                                lu_tot_npixel(:,:,:)
                                          !< total number of lu raw data pixels on target grid (dimension (ie,je,ke))

    REAL (KIND=wp), INTENT(OUT)        :: fr_land_lu(:,:,:), & !< fraction land due to lu raw data
         &                                ice_lu(:,:,:), &     !< fraction of ice due to lu raw data
         &                                z0_lu(:,:,:), &      !< roughness length due to lu land use data
         &                                root_lu(:,:,:), &    !< root depth due to lu land use data
         &                                plcov_mn_lu(:,:,:), &!< plant cover maximum due to lu land use data
         &                                plcov_mx_lu(:,:,:), &!< plant cover minimum due to lu land use data
         &                                lai_mn_lu(:,:,:), &  !< Leaf Area Index maximum due to lu land use data
         &                                lai_mx_lu(:,:,:), &  !< Leaf Area Index minimum due to lu land use data
         &                                rs_min_lu(:,:,:), &  !< minimal stomata resistance due to lu land use data
         &                                urban_lu(:,:,:), &   !< urban fraction due to lu land use data
         &                                for_d_lu(:,:,:), &   !< deciduous forest (fraction) due to lu land use data
         &                                for_e_lu(:,:,:), &   !< evergreen forest (fraction) due to lu land use data
         &                                skinc_lu(:,:,:), &   !< skin conductivity due to lu land use data
         &                                emissivity_lu(:,:,:), & !< longwave emissivity due to lu land use data
         &                                lu_class_fraction(:,:,:,:)!< fraction for each lu class on target grid


    CALL logging%info('Enter routine: read_netcdf_buffer_lu')

    !set up dimensions for buffer
    CALL  def_dimension_info_buffer(tg)
    ! dim_3d_tg

    ! define meta information for various land use related variables for netcdf output
    CALL def_lu_fields_meta(nclass_lu,dim_3d_tg)

    ! define meta information for target field variables lon_geo, lat_geo
    CALL def_com_target_fields_meta(dim_3d_tg)
    ! lon_geo_meta and lat_geo_meta

    CALL netcdf_get_var(TRIM(netcdf_filename),fr_land_lu_meta,fr_land_lu)

    CALL netcdf_get_var(TRIM(netcdf_filename),lu_tot_npixel_meta,lu_tot_npixel)

    CALL netcdf_get_var(TRIM(netcdf_filename),lu_class_fraction_meta,lu_class_fraction)

    ! CALL netcdf_get_var(TRIM(netcdf_filename),lu_class_npixel_meta,lu_class_npixel)
    CALL my_get_var_int_4d(TRIM(netcdf_filename),lu_class_npixel_meta,lu_class_npixel)

    CALL netcdf_get_var(TRIM(netcdf_filename),ice_lu_meta,ice_lu)

    CALL netcdf_get_var(TRIM(netcdf_filename),z0_lu_meta,z0_lu)

    CALL netcdf_get_var(TRIM(netcdf_filename),plcov_mx_lu_meta,plcov_mx_lu)

    CALL netcdf_get_var(TRIM(netcdf_filename),plcov_mn_lu_meta,plcov_mn_lu)

    CALL netcdf_get_var(TRIM(netcdf_filename),lai_mx_lu_meta,lai_mx_lu)

    CALL netcdf_get_var(TRIM(netcdf_filename),lai_mn_lu_meta,lai_mn_lu)

    CALL netcdf_get_var(TRIM(netcdf_filename),rs_min_lu_meta,rs_min_lu)

    CALL netcdf_get_var(TRIM(netcdf_filename),urban_lu_meta,urban_lu)

    CALL netcdf_get_var(TRIM(netcdf_filename),for_d_lu_meta,for_d_lu)

    CALL netcdf_get_var(TRIM(netcdf_filename),for_e_lu_meta,for_e_lu)

    CALL netcdf_get_var(TRIM(netcdf_filename),skinc_lu_meta,skinc_lu)

    CALL netcdf_get_var(TRIM(netcdf_filename),emissivity_lu_meta,emissivity_lu)

    CALL netcdf_get_var(TRIM(netcdf_filename),root_lu_meta,root_lu)

    IF (l_terra_urb) CALL terra_urb_read_netcdf(netcdf_filename)

    CALL logging%info('Exit routine: read_netcdf_buffer_lu')

  END SUBROUTINE read_netcdf_buffer_lu

  !> netcdf output of GLC2000 derived buffer fields
  SUBROUTINE write_netcdf_buffer_glc2000(netcdf_filename,  &
       &                                 tg,         &
       &                                 undefined, &
       &                                 undef_int,   &
       &                                 lon_geo,     &
       &                                 lat_geo, &
       &                                 fr_land_glc2000, &
       &                                 glc2000_class_fraction,    &
       &                                 glc2000_class_npixel, &
       &                                 glc2000_tot_npixel, &
       &                                 ice_glc2000, &
       &                                 z0_glc2000, &
       &                                 root_glc2000, &
       &                                 plcov_mn_glc2000, &
       &                                 plcov_mx_glc2000, &
       &                                 lai_mn_glc2000, &
       &                                 lai_mx_glc2000, &
       &                                 rs_min_glc2000, &
       &                                 urban_glc2000,  &
       &                                 for_d_glc2000,  &
       &                                 for_e_glc2000, &
       &                                 emissivity_glc2000)


    CHARACTER (len=*), INTENT(IN)      :: netcdf_filename !< filename for the netcdf file

    TYPE(target_grid_def), INTENT(IN)  :: tg !< structure with target grid description

    REAL(KIND=wp), INTENT(IN)          :: undefined, &       !< value to indicate undefined grid elements
         &                                lon_geo(:,:,:), &  !< longitude coordinates of the target grid in the geographical system
         &                                lat_geo(:,:,:), &  !< latitude coordinates of the target grid in the geographical system
         &                                glc2000_class_fraction(:,:,:,:)

    INTEGER (KIND=i4), INTENT(IN)      :: glc2000_class_npixel(:,:,:,:), &
         &                                glc2000_tot_npixel(:,:,:), &
         &                                undef_int

    REAL (KIND=wp), INTENT(IN)         :: fr_land_glc2000(:,:,:), & !< fraction land due to glc2000 raw data
         &                                ice_glc2000(:,:,:), &     !< fraction of ice due to glc2000 raw data
         &                                z0_glc2000(:,:,:), &      !< roughness length due to glc2000 land use data
         &                                root_glc2000(:,:,:), &    !< root depth due to glc2000 land use data
         &                                plcov_mx_glc2000(:,:,:), &!< plant cover maximum due to glc2000 land use data
         &                                plcov_mn_glc2000(:,:,:), &!< plant cover minimum due to glc2000 land use data
         &                                lai_mx_glc2000(:,:,:), &  !< Leaf Area Index maximum due to glc2000 land use data
         &                                lai_mn_glc2000(:,:,:), &  !< Leaf Area Index minimum due to glc2000 land use data
         &                                rs_min_glc2000(:,:,:), &  !< minimal stomata resistance due to glc2000 land use data
         &                                urban_glc2000(:,:,:), &   !< urban fraction due to glc2000 land use data
         &                                for_d_glc2000(:,:,:), &   !< deciduous forest (fraction) due to glc2000 land use data
         &                                for_e_glc2000(:,:,:), &   !< evergreen forest (fraction) due to glc2000 land use data
         &                                emissivity_glc2000(:,:,:) !< longwave emissivity due to glc2000 land use data


    ! local variables
    INTEGER(KIND=i4)                   :: ndims, ncid, undefined_i, errorcode

    INTEGER(KIND=i4), PARAMETER        :: nglob_atts=6

    TYPE(dim_meta_info), ALLOCATABLE   :: dim_list(:) !< dimensions for netcdf file
    TYPE(netcdf_attributes)            :: global_attributes(nglob_atts)

    CALL logging%info('Enter routine: write_netcdf_buffer_glc2000')

    ! define global attributes
    CALL set_global_att_glc2000(global_attributes)

    !set up dimensions for buffer
    CALL  def_dimension_info_buffer(tg)
    ! dim_3d_tg

    ! define meta information for various land use related variables (GLC2000) for netcdf output
    CALL def_glc2000_fields_meta(nclass_glc2000,dim_3d_tg)

    ! define meta information for target field variables lon_geo, lat_geo
    CALL def_com_target_fields_meta(dim_3d_tg)
    ! lon_geo_meta and lat_geo_meta

    !set up dimensions for buffer netcdf output
    ndims = 4
    ALLOCATE(dim_list(1:ndims),STAT=errorcode)
    IF (errorcode /= 0 ) CALL logging%error('Cant allocate array dim_list',__FILE__,__LINE__)
    dim_list = dim_glc2000_tg

    undefined_i = undef_int

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

    CALL logging%info('Exit routine: write_netcdf_buffer_glc2000')

  END SUBROUTINE write_netcdf_buffer_glc2000
  !-----------------------------------------------------------------------
  !> set global attributes for netcdf with glc2000 data
  SUBROUTINE set_global_att_glc2000(global_attributes)

    TYPE(netcdf_attributes), INTENT(INOUT) :: global_attributes(1:6)

    !local variables
    CHARACTER(len=10)                      :: ydate, &
         &                                    ytime

    CHARACTER(len=2)                       :: cc, &
         &                                    yy, &
         &                                    mm, &
         &                                    dd, &
         &                                    hh, &
         &                                    minute

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

  !> netcdf output of GLCC derived buffer fields
  SUBROUTINE write_netcdf_buffer_glcc(netcdf_filename,  &
       &                              tg,         &
       &                              undefined, &
       &                              undef_int,   &
       &                              lon_geo,     &
       &                              lat_geo, &
       &                              fr_land_glcc, &
       &                              glcc_class_fraction,    &
       &                              glcc_class_npixel, &
       &                              glcc_tot_npixel, &
       &                              ice_glcc, &
       &                              z0_glcc, &
       &                              root_glcc, &
       &                              plcov_mn_glcc, &
       &                              plcov_mx_glcc, &
       &                              lai_mn_glcc, &
       &                              lai_mx_glcc, &
       &                              rs_min_glcc, &
       &                              urban_glcc,  &
       &                              for_d_glcc,  &
       &                              for_e_glcc, &
       &                              emissivity_glcc)

    CHARACTER (len=*), INTENT(IN)     :: netcdf_filename !< filename for the netcdf file

    TYPE(target_grid_def), INTENT(IN) :: tg !< structure with target grid description

    INTEGER(KIND=i4), INTENT(IN)      :: undef_int, &       !< value to indicate undefined grid elements
         &                               glcc_class_npixel(:,:,:,:), &
         &                               glcc_tot_npixel(:,:,:)

    REAL(KIND=wp), INTENT(IN)         :: undefined, &       !< value to indicate undefined grid elements
         &                               glcc_class_fraction(:,:,:,:), &
         &                               lon_geo(:,:,:), &  !< longitude coordinates of the target grid in the geographical system
         &                               lat_geo(:,:,:), &  !< latitude coordinates of the target grid in the geographical system
         &                               fr_land_glcc(:,:,:), & !< fraction land due to glcc raw data
         &                               ice_glcc(:,:,:), &     !< fraction of ice due to glcc raw data
         &                               z0_glcc(:,:,:), &      !< roughness length due to glcc land use data
         &                               root_glcc(:,:,:), &    !< root depth due to glcc land use data
         &                               plcov_mx_glcc(:,:,:), &!< plant cover maximum due to glcc land use data
         &                               plcov_mn_glcc(:,:,:), &!< plant cover minimum due to glcc land use data
         &                               lai_mx_glcc(:,:,:), &  !< Leaf Area Index maximum due to glcc land use data
         &                               lai_mn_glcc(:,:,:), &  !< Leaf Area Index minimum due to glcc land use data
         &                               rs_min_glcc(:,:,:), &  !< minimal stomata resistance due to glcc land use data
         &                               urban_glcc(:,:,:), &   !< urban fraction due to glcc land use data
         &                               for_d_glcc(:,:,:), &   !< deciduous forest (fraction) due to glcc land use data
         &                               for_e_glcc(:,:,:), &   !< evergreen forest (fraction) due to glcc land use data
         &                               emissivity_glcc(:,:,:) !< longwave emissivity due to glcc land use data


    ! local variables
    INTEGER (KIND=i4)                 :: undefined_i, ndims, ncid, errorcode
    INTEGER(KIND=i4), PARAMETER       :: nglob_atts=6

    TYPE(dim_meta_info), ALLOCATABLE  :: dim_list(:) !< dimensions for netcdf file

    TYPE(netcdf_attributes)           :: global_attributes(nglob_atts)

    CALL logging%info('Enter routine: write_netcdf_buffer_glcc')

    ! define global attributes
    CALL set_global_att_glcc(global_attributes)

    !set up dimensions for buffer
    CALL  def_dimension_info_buffer(tg)
    ! dim_3d_tg

    ! define meta information for various land use related variables (GLCC) for netcdf output
    CALL def_glcc_fields_meta(nclass_glcc,dim_3d_tg)

    ! define meta information for target field variables lon_geo, lat_geo
    CALL def_com_target_fields_meta(dim_3d_tg)
    ! lon_geo_meta and lat_geo_meta

    !set up dimensions for buffer netcdf output
    ndims = 4
    ALLOCATE(dim_list(1:ndims),STAT=errorcode)
    IF (errorcode /= 0 ) CALL logging%error('Cant allocate array dim_list',__FILE__,__LINE__)
    dim_list = dim_glcc_tg

    !  dim_3d_buffer(:) = dim_list(:)

    undefined_i = undef_int

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

    CALL close_netcdf_file(ncid)

    CALL logging%info('Exit routine: write_netcdf_buffer_glcc')

  END SUBROUTINE write_netcdf_buffer_glcc

  !-----------------------------------------------------------------------
  !> set global attributes for netcdf with glcc data
  SUBROUTINE set_global_att_glcc(global_attributes)

    TYPE(netcdf_attributes), INTENT(INOUT) :: global_attributes(1:6)

    !local variables
    CHARACTER(len=10)                      :: ydate, &
         &                                    ytime

    CHARACTER(len=2)                       :: cc, &
         &                                    yy, &
         &                                    mm, &
         &                                    dd, &
         &                                    hh, &
         &                                    minute

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

  !> read GLCC derived buffer fields
  SUBROUTINE read_netcdf_buffer_glcc(netcdf_filename,  &
       &                             tg,         &
       &                             fr_land_glcc, &
       &                             glcc_class_fraction,    &
       &                             glcc_class_npixel, &
       &                             glcc_tot_npixel, &
       &                             ice_glcc, &
       &                             z0_glcc, &
       &                             root_glcc, &
       &                             plcov_mn_glcc, &
       &                             plcov_mx_glcc, &
       &                             lai_mn_glcc, &
       &                             lai_mx_glcc, &
       &                             rs_min_glcc, &
       &                             urban_glcc,  &
       &                             for_d_glcc,  &
       &                             for_e_glcc, &
       &                             emissivity_glcc)

    CHARACTER (len=*), INTENT(IN)      :: netcdf_filename !< filename for the netcdf file

    TYPE(target_grid_def), INTENT(IN)  :: tg !< structure with target grid description

    INTEGER (KIND=i4), INTENT(OUT)     :: glcc_class_npixel(:,:,:,:), &
         &                                glcc_tot_npixel(:,:,:)

    REAL (KIND=wp), INTENT(OUT)        :: glcc_class_fraction(:,:,:,:), &
         &                                fr_land_glcc(:,:,:), & !< fraction land due to glcc raw data
         &                                ice_glcc(:,:,:), &     !< fraction of ice due to glcc raw data
         &                                z0_glcc(:,:,:), &      !< roughness length due to glcc land use data
         &                                root_glcc(:,:,:), &    !< root depth due to glcc land use data
         &                                plcov_mn_glcc(:,:,:), &!< plant cover maximum due to glcc land use data
         &                                plcov_mx_glcc(:,:,:), &!< plant cover minimum due to glcc land use data
         &                                lai_mn_glcc(:,:,:), &  !< Leaf Area Index maximum due to glcc land use data
         &                                lai_mx_glcc(:,:,:), &  !< Leaf Area Index minimum due to glcc land use data
         &                                rs_min_glcc(:,:,:), &  !< minimal stomata resistance due to glcc land use data
         &                                urban_glcc(:,:,:), &   !< urban fraction due to glcc land use data
         &                                for_d_glcc(:,:,:), &   !< deciduous forest (fraction) due to glcc land use data
         &                                for_e_glcc(:,:,:), &   !< evergreen forest (fraction) due to glcc land use data
         &                                emissivity_glcc(:,:,:) !< longwave emissivity due to glcc land use data

    CALL logging%info('Enter routine: read_netcdf_buffer_glcc')

    !set up dimensions for buffer
    CALL  def_dimension_info_buffer(tg)
    ! dim_3d_tg

    ! define meta information for various land use related variables (GLCC) for netcdf output
    CALL def_glcc_fields_meta(nclass_glcc,dim_3d_tg)

    ! define meta information for target field variables lon_geo, lat_geo
    CALL def_com_target_fields_meta(dim_3d_tg)
    ! lon_geo_meta and lat_geo_meta

    CALL netcdf_get_var(TRIM(netcdf_filename),fr_land_glcc_meta,fr_land_glcc)

    CALL netcdf_get_var(TRIM(netcdf_filename),glcc_tot_npixel_meta,glcc_tot_npixel)

    CALL netcdf_get_var(TRIM(netcdf_filename),glcc_class_fraction_meta,glcc_class_fraction)

    ! CALL netcdf_get_var(TRIM(netcdf_filename),glcc_class_npixel_meta,glcc_class_npixel)
    CALL my_get_var_int_4d(TRIM(netcdf_filename),glcc_class_npixel_meta,glcc_class_npixel)

    CALL netcdf_get_var(TRIM(netcdf_filename),ice_glcc_meta,ice_glcc)

    CALL netcdf_get_var(TRIM(netcdf_filename),z0_glcc_meta,z0_glcc)

    CALL netcdf_get_var(TRIM(netcdf_filename),plcov_mx_glcc_meta,plcov_mx_glcc)

    CALL netcdf_get_var(TRIM(netcdf_filename),plcov_mn_glcc_meta,plcov_mn_glcc)

    CALL netcdf_get_var(TRIM(netcdf_filename),lai_mx_glcc_meta,lai_mx_glcc)

    CALL netcdf_get_var(TRIM(netcdf_filename),lai_mn_glcc_meta,lai_mn_glcc)

    CALL netcdf_get_var(TRIM(netcdf_filename),rs_min_glcc_meta,rs_min_glcc)

    CALL netcdf_get_var(TRIM(netcdf_filename),urban_glcc_meta,urban_glcc)

    CALL netcdf_get_var(TRIM(netcdf_filename),for_d_glcc_meta,for_d_glcc)

    CALL netcdf_get_var(TRIM(netcdf_filename),for_e_glcc_meta,for_e_glcc)

    CALL netcdf_get_var(TRIM(netcdf_filename),emissivity_glcc_meta,emissivity_glcc)

    CALL netcdf_get_var(TRIM(netcdf_filename),root_glcc_meta,root_glcc)

    CALL logging%info('Exit routine: read_netcdf_buffer_glcc')

  END SUBROUTINE read_netcdf_buffer_glcc

END Module mo_landuse_output_nc
