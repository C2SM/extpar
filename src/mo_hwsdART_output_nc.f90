MODULE mo_hwsdART_output_nc

  USE mo_logging
  !> kind parameters are defined in MODULE data_parameters
  USE mo_kind, ONLY: wp, &
        &          i4

  USE mo_grid_structures, ONLY: reg_lonlat_grid, &
   & rotated_lonlat_grid, &
   & icosahedral_triangular_grid, &
   & target_grid_def

  USE mo_io_utilities, ONLY: var_meta_info, &
   &  netcdf_attributes, &
   & dim_meta_info, &
   & netcdf_put_var, &
   & open_new_netcdf_file, &
   & close_netcdf_file, &
   & netcdf_def_grid_mapping, &
   & vartype_int , &
   & vartype_real, &
   & vartype_char, &
   & check_netcdf, &
   & dim_meta_info
  

  USE mo_icon_grid_data, ONLY: icon_grid_region, &
  & clon, &
  & clat, &
  & clon_vertices, &
  & clat_vertices, &
  & allocate_icon_coor

  USE mo_hwsdART_tg_fields, ONLY: fr_heavy_clay, &
  & fr_silty_clay, &
  & fr_light_clay, &
  & fr_silty_clay_loam,  &
  & fr_clay_loam, &
  & fr_silt, &
  & fr_silt_loam, &
  & fr_sandy_clay, &
  & fr_loam, &
  & fr_sandy_clay_loam, &
  & fr_sandy_loam,        &
  & fr_loamy_sand, &
  & fr_sand,fr_undef

  USE mo_var_meta_data, ONLY: dim_3d_tg, &
    & def_dimension_info_buffer, &
    & lon_geo_meta, &
    & lat_geo_meta, &
    & no_raw_data_pixel_meta, &
    & def_com_target_fields_meta , & 
    & def_dimension_info_icon, &
    & nc_grid_def_icon, &
    & set_nc_grid_def_icon, &
    & nc_grid_def_cosmo, &
    & set_nc_grid_def_cosmo, &
    & dim_rlon_cosmo, &
    & dim_rlat_cosmo, &
    & dim_2d_cosmo,   &
    & rlon_meta,      &
    & rlat_meta,      &
    & def_dimension_info_cosmo, &
    & clon_meta, &
    & clat_meta, &
    & clon_vertices_meta, &
    & clat_vertices_meta

  USE mo_cosmo_grid, ONLY: lon_rot, lat_rot


  IMPLICIT NONE

  PRIVATE


  PUBLIC :: write_netcdf_hwsdART_icon_grid, &
  & write_netcdf_hwsdART_cosmo_grid

  TYPE(var_meta_info) :: fr_heavy_clay_meta, &
  &  fr_silty_clay_meta, &
  &  fr_light_clay_meta, &
  &  fr_silty_clay_loam_meta, &
  &  fr_clay_loam_meta, &
  &  fr_silt_meta, &
  &  fr_silt_loam_meta, &
  &  fr_sandy_clay_meta, &
  &  fr_loam_meta, &
  &  fr_sandy_clay_loam_meta, &
  &  fr_sandy_loam_meta, &
  &  fr_loamy_sand_meta, &
  &  fr_sand_meta, &
  &  fr_undef_meta
  
  character (len=1), parameter :: c_undef = "-" 

  CONTAINS


  !----------------------------------------------------------------------------
  !> set global attributes for netcdf with hwsdARTtype data
  SUBROUTINE set_global_att_hwsdARTtype(global_attributes,icon_grid)
    TYPE(icosahedral_triangular_grid), INTENT(IN),OPTIONAL :: icon_grid !< structure which contains the definition of the ICON grid
    TYPE(netcdf_attributes), INTENT(INOUT) :: global_attributes(1:8)

    !local variables
    CHARACTER(len=10) :: ydate, &
    &  ytime
    CHARACTER(len=2)  :: number_Of_Grid_Used_string, &
    & cc, &
    & yy, &
    & mm, &
    & dd, &
    & hh, &
    & minute
    
    CHARACTER(len=1 ) :: uuid(16)    !   UUID of unstructured grids 
    CHARACTER(len=16 ) :: uuidtxt    !   UUID of unstructured grids 

    INTEGER :: i

    ! define global attributes
      
    global_attributes(1)%attname = 'title'
    global_attributes(1)%attributetext='Percentage of HWSD soil type'
    global_attributes(2)%attname = 'institution'
    global_attributes(2)%attributetext='KIT'

    global_attributes(3)%attname = 'source'
    global_attributes(3)%attributetext='HWSD Harmonized World Soil Database'

    CALL DATE_AND_TIME(ydate,ytime)
    READ(ydate,'(4A2)') cc,yy,mm,dd
    READ(ytime,'(2A2)') hh, minute

    ydate=TRIM(cc)//TRIM(yy)//'-'//TRIM(mm)//'-'//TRIM(dd)
    ytime=TRIM(hh)//':'//TRIM(minute) 

    global_attributes(4)%attname = 'history'
    global_attributes(4)%attributetext=TRIM(ydate)//'T'//TRIM(ytime)//' hwsdART_to_buffer'

    global_attributes(5)%attname = 'references'
    global_attributes(5)%attributetext='HWSD Harmonized World Soil Database'

    global_attributes(6)%attname = 'comment'
    global_attributes(6)%attributetext='For generating mineral dust emissions in ICON-ART and COSMO-ART'
    
    IF(PRESENT(icon_grid)) then
    
      write(number_Of_Grid_Used_string,'(I2)')  icon_grid%number_Of_Grid_Used
      global_attributes(7)%attname = 'number_of_grid_used'
      global_attributes(7)%attributetext=number_Of_Grid_Used_string

      call decode_uuid (icon_grid%uuidOfHGrid, uuid) 
      do i=1,len(uuid)
         uuidtxt(i:i)=uuid(i)
      end do

      global_attributes(8)%attname = 'uuidOfHGrid'
      global_attributes(8)%attributetext=icon_grid%uuidOfHGrid
    ENDIF
  END SUBROUTINE set_global_att_hwsdARTtype
  !----------------------------------------------------------------------------

 !> create a netcdf file for data to the ICON grid
 SUBROUTINE write_netcdf_hwsdART_icon_grid(netcdf_filename, &
   &                                     icon_grid,      &
   &                                     tg,             &
   &                                     undefined       )







    CHARACTER (len=*), INTENT(IN)      :: netcdf_filename !< filename for the netcdf file
    TYPE(icosahedral_triangular_grid), INTENT(IN)  :: icon_grid      !< structure which contains the definition of the ICON grid
    TYPE(target_grid_def), INTENT(IN)  :: tg !< structure with target grid description
    REAL (KIND=wp), INTENT(IN) ::   undefined       !< value to indicate undefined grid elements

    ! local variables
    INTEGER, PARAMETER :: nglob_atts=8
    TYPE(netcdf_attributes) :: global_attributes(nglob_atts)
    INTEGER :: ncid, &                             !< netcdf unit file number
    &  ndims, &
    &  errorcode, & !< error status variable
    &  vert_id, &
    &  nc, &
    &  nv !< counters

    TYPE(dim_meta_info), ALLOCATABLE :: dim_list(:) !< dimensions for netcdf file
    TYPE(dim_meta_info), TARGET :: dim_1d_icon(1:1), &
    & dim_1d_icon_v(1:1)

    CHARACTER (len=80):: grid_mapping, & !< netcdf attribute grid mapping
    & coordinates  !< netcdf attribute coordinates

   


    WRITE(message_text,*) 'ENTER write_netcdf_hwsdART_icon_grid'
    CALL logging%info(message_text)  
    !-------------------------------------------------------------
    ! organize output

    !-------------------------------------------------------------
    !set up dimensions for buffer netcdf output 
    ndims = 3
    ALLOCATE(dim_list(1:ndims),STAT=errorcode)
    IF (errorcode /= 0 ) call logging%error ('Cant allocate array dim_list',__FILE__,__LINE__) 
    dim_list(1)%dimname = 'cell'
    dim_list(1)%dimsize = icon_grid%ncell 
    dim_list(2)%dimname = 'nv' ! icon_grid%nvertex_per_cell
    dim_list(2)%dimsize = icon_grid%nvertex_per_cell 
    dim_list(3)%dimname = 'vertex'
    dim_list(3)%dimsize = icon_grid%nvertex

    dim_1d_icon =  dim_list(1) ! cell
    dim_1d_icon_v =  dim_list(3) ! nvertex
    
    ! set Icon coordinates for output
    CALL allocate_icon_coor(icon_grid%ncell, icon_grid%nvertex_per_cell)

    clon(:) = icon_grid_region%cells%center(:)%lon
    clat(:) = icon_grid_region%cells%center(:)%lat

    DO nc=1,icon_grid%ncell
      DO nv=1,icon_grid%nvertex_per_cell
      vert_id =  icon_grid_region%cells%vertex_index(nc,nv)
      clon_vertices(nv,nc) =  icon_grid_region%verts%vertex(vert_id)%lon
      clat_vertices(nv,nc) =  icon_grid_region%verts%vertex(vert_id)%lat
      ENDDO
    ENDDO  
    
    
    
   !-------------------------------------------------------------
   ! define global attributes
   CALL set_global_att_hwsdARTtype(global_attributes,icon_grid)
   
   !set up dimensions for buffer
   CALL  def_dimension_info_buffer(tg)
   ! dim_3d_tg
   WRITE(message_text,*) 'def_com_target_fields_meta'
   CALL logging%info(message_text)  
   !set up dimensions for ICON grid
   CALL def_dimension_info_icon(icon_grid)
   ! dim_icon


  
  ! define meta information for target field variables lon_geo, lat_geo 
  CALL def_com_target_fields_meta(dim_1d_icon)
  ! lon_geo_meta and lat_geo_meta
  ! set mapping parameters for netcdf
  grid_mapping="lon_lat_on_sphere"
  coordinates="clon clat"
  CALL set_nc_grid_def_icon(grid_mapping)
    
       
    !-----------------------------------------------------------------
    WRITE(message_text,*) ' CALL open_new_netcdf_file'
    CALL logging%info(message_text)  
    CALL open_new_netcdf_file(netcdf_filename=TRIM(netcdf_filename),   &
      &                       dim_list=dim_list,                  &
      &                       global_attributes=global_attributes, &
      &                       ncid=ncid)
    !-----------------------------------------------------------------
    ! FILL METADATA
    CALL create_hwsdART_meta(fr_heavy_clay_meta,     dim_1d_icon,'fr_hcla', &
                          &  'Fraction of Heavy Clay')
    CALL create_hwsdART_meta(fr_silty_clay_meta,     dim_1d_icon,'fr_silc', &
                          &  'Fraction of Silty Clay')
    CALL create_hwsdART_meta(fr_light_clay_meta,     dim_1d_icon,'fr_lcla', &
                          &  'Fraction of Light Clay')
    CALL create_hwsdART_meta(fr_silty_clay_loam_meta,dim_1d_icon,'fr_sicl', &
                          &  'Fraction of Silty Clay Loam')
    CALL create_hwsdART_meta(fr_clay_loam_meta,      dim_1d_icon,'fr_cloa', &
                          &  'Fraction of Clay Loam')
    CALL create_hwsdART_meta(fr_silt_meta,           dim_1d_icon,'fr_silt', &
                          &  'Fraction of Silt')
    CALL create_hwsdART_meta(fr_silt_loam_meta,      dim_1d_icon,'fr_silo', &
                          &  'Fraction of Silty Loam')
    CALL create_hwsdART_meta(fr_sandy_clay_meta,     dim_1d_icon,'fr_scla', &
                          &  'Fraction of Sandy Clay')
    CALL create_hwsdART_meta(fr_loam_meta,           dim_1d_icon,'fr_loam', &
                          &  'Fraction of Loam')
    CALL create_hwsdART_meta(fr_sandy_clay_loam_meta,dim_1d_icon,'fr_sclo', &
                          &  'Fraction of Sandy Clay Loam')
    CALL create_hwsdART_meta(fr_sandy_loam_meta,     dim_1d_icon,'fr_sloa', &
                          &  'Fraction of Sandy Loam')
    CALL create_hwsdART_meta(fr_loamy_sand_meta,     dim_1d_icon,'fr_lsan', &
                          &  'Fraction of Loamy Sand')
    CALL create_hwsdART_meta(fr_sand_meta,           dim_1d_icon,'fr_sand', &
                          &  'Fraction of Sand')
    CALL create_hwsdART_meta(fr_undef_meta,          dim_1d_icon,'fr_udef', &
                          &  'Fraction of Undefined or Water')
                          
    !-----------------------------------------------------------------

    CALL netcdf_put_var(ncid,fr_heavy_clay(1:icon_grid%ncell,1,1),     fr_heavy_clay_meta,     undefined)
    CALL netcdf_put_var(ncid,fr_silty_clay(1:icon_grid%ncell,1,1),     fr_silty_clay_meta,     undefined)
    CALL netcdf_put_var(ncid,fr_light_clay(1:icon_grid%ncell,1,1),     fr_light_clay_meta,     undefined)
    CALL netcdf_put_var(ncid,fr_silty_clay_loam(1:icon_grid%ncell,1,1),fr_silty_clay_loam_meta,undefined)
    CALL netcdf_put_var(ncid,fr_clay_loam(1:icon_grid%ncell,1,1),      fr_clay_loam_meta,      undefined)
    CALL netcdf_put_var(ncid,fr_silt(1:icon_grid%ncell,1,1),           fr_silt_meta,           undefined)
    CALL netcdf_put_var(ncid,fr_silt_loam(1:icon_grid%ncell,1,1),      fr_silt_loam_meta,      undefined)
    CALL netcdf_put_var(ncid,fr_sandy_clay(1:icon_grid%ncell,1,1),     fr_sandy_clay_meta,     undefined)
    CALL netcdf_put_var(ncid,fr_loam(1:icon_grid%ncell,1,1),           fr_loam_meta,           undefined)
    CALL netcdf_put_var(ncid,fr_sandy_clay_loam(1:icon_grid%ncell,1,1),fr_sandy_clay_loam_meta,undefined)
    CALL netcdf_put_var(ncid,fr_sandy_loam(1:icon_grid%ncell,1,1),     fr_sandy_loam_meta,     undefined)
    CALL netcdf_put_var(ncid,fr_loamy_sand(1:icon_grid%ncell,1,1),     fr_loamy_sand_meta,     undefined)
    CALL netcdf_put_var(ncid,fr_sand(1:icon_grid%ncell,1,1),           fr_sand_meta,           undefined)
    CALL netcdf_put_var(ncid,fr_undef(1:icon_grid%ncell,1,1),          fr_undef_meta,          undefined)
    

  !<<RR>>
  WRITE(message_text,*) 'Putting var: clon'
  CALL logging%info(message_text)  
  CALL  netcdf_put_var(ncid,clon,clon_meta,undefined)
  WRITE(message_text,*) 'Putting var: clat'
  CALL logging%info(message_text)  
  CALL  netcdf_put_var(ncid,clat,clat_meta,undefined)
  WRITE(message_text,*) 'Putting var: clon_vertices'
  CALL logging%info(message_text)  
  CALL  netcdf_put_var(ncid,clon_vertices,clon_vertices_meta,undefined)
  WRITE(message_text,*) 'Putting var: clat_vertices'
  CALL logging%info(message_text)  
  CALL  netcdf_put_var(ncid,clat_vertices,clat_vertices_meta,undefined)
  !<<RR>>
  WRITE(message_text,*) 'Closing the file'
  CALL logging%info(message_text)  
  CALL close_netcdf_file(ncid)

END SUBROUTINE write_netcdf_hwsdART_icon_grid

  !> create a netcdf file for data to the COSMO grid
 SUBROUTINE write_netcdf_hwsdART_cosmo_grid(netcdf_filename,&
   &                                     cosmo_grid,     &
   &                                     tg,             &
   &                                     undefined       )



  CHARACTER (len=*), INTENT(IN)      :: netcdf_filename !< filename for the netcdf file
  TYPE(rotated_lonlat_grid), INTENT(IN)  :: cosmo_grid      !< structure which contains the definition of the COSMO grid
  TYPE(target_grid_def), INTENT(IN)  :: tg !< structure with target grid description
  REAL(KIND=wp), INTENT(IN)          :: undefined       !< value to indicate undefined grid elements
  ! local variables

  INTEGER, PARAMETER :: nglob_atts=8
  TYPE(netcdf_attributes) :: global_attributes(nglob_atts)
  TYPE(dim_meta_info), ALLOCATABLE :: dim_list(:) !< dimensions for netcdf file

  INTEGER :: ndims, &  
   &  ncid, &
   &  varid, &
   &  errorcode        !< error status variable


  CHARACTER (len=80):: grid_mapping, & !< netcdf attribute grid mapping
   &  coordinates  !< netcdf attribute coordinates


    
  !-------------------------------------------------------------
  ! define global attributes
  CALL set_global_att_hwsdARTtype(global_attributes)

  !set up dimensions for buffer
  CALL  def_dimension_info_buffer(tg)
  ! dim_2d_tg

  !set up dimensions for COSMO grid
  CALL def_dimension_info_cosmo(cosmo_grid)
  ! dim_rlon_cosmo, dim_rlat_cosmo, dim_2d_cosmo, rlon_meta, rlat_meta

  ! set mapping parameters for netcdf
   grid_mapping="rotated_pole"
   coordinates="lon lat"
   CALL set_nc_grid_def_cosmo(cosmo_grid,grid_mapping)
   ! nc_grid_def_cosmo

   ! define meta information for target field variables lon_geo, lat_geo 
  CALL def_com_target_fields_meta(dim_2d_cosmo,coordinates,grid_mapping)
  ! lon_geo_meta and lat_geo_meta

    !-----------------------------------------------------------------
    ! FILL METADATA
    CALL create_hwsdART_meta(fr_heavy_clay_meta,     dim_2d_cosmo,'fr_hcla', &
                          &  'Fraction of Heavy Clay')
    CALL create_hwsdART_meta(fr_silty_clay_meta,     dim_2d_cosmo,'fr_silc', &
                          &  'Fraction of Silty Clay')
    CALL create_hwsdART_meta(fr_light_clay_meta,     dim_2d_cosmo,'fr_lcla', &
                          &  'Fraction of Light Clay')
    CALL create_hwsdART_meta(fr_silty_clay_loam_meta,dim_2d_cosmo,'fr_sicl', &
                          &  'Fraction of Silty Clay Loam')
    CALL create_hwsdART_meta(fr_clay_loam_meta,      dim_2d_cosmo,'fr_cloa', &
                          &  'Fraction of Clay Loam')
    CALL create_hwsdART_meta(fr_silt_meta,           dim_2d_cosmo,'fr_silt', &
                          &  'Fraction of Silt')
    CALL create_hwsdART_meta(fr_silt_loam_meta,      dim_2d_cosmo,'fr_silo', &
                          &  'Fraction of Silt Loam')
    CALL create_hwsdART_meta(fr_sandy_clay_meta,     dim_2d_cosmo,'fr_scla', &
                          &  'Fraction of Sandy Clay')
    CALL create_hwsdART_meta(fr_loam_meta,           dim_2d_cosmo,'fr_loam', &
                          &  'Fraction of Loam')
    CALL create_hwsdART_meta(fr_sandy_clay_loam_meta,dim_2d_cosmo,'fr_sclo', &
                          &  'Fraction of Sandy Clay Loam')
    CALL create_hwsdART_meta(fr_sandy_loam_meta,     dim_2d_cosmo,'fr_sloa', &
                          &  'Fraction of Sandy Loam')
    CALL create_hwsdART_meta(fr_loamy_sand_meta,     dim_2d_cosmo,'fr_lsan', &
                          &  'Fraction of Loamy Sand')
    CALL create_hwsdART_meta(fr_sand_meta,           dim_2d_cosmo,'fr_sand', &
                          &  'Fraction of Sand')
    CALL create_hwsdART_meta(fr_undef_meta,          dim_2d_cosmo,'fr_udef', &
                          &  'Fraction of Undefined or Water')

   !-------------------------------------------------------------
   ! organize output

   ! set up dimensions for COSMO grid
    ndims = 2

    ALLOCATE(dim_list(1:ndims),STAT=errorcode)
    IF (errorcode /= 0 ) call logging%error ('Cant allocate array dim_list',__FILE__,__LINE__) 

    dim_list(1) = dim_rlon_cosmo(1) ! rlon
    dim_list(2) = dim_rlat_cosmo(1) ! rlat

    dim_rlon_cosmo(1) = dim_list(1)
    dim_rlat_cosmo(1) = dim_list(2)

    
   !-----------------------------------------------------------------
    WRITE(MESSAGE_TEXT,*) ' CALL open_new_netcdf_file'
    CALL logging%info(message_text)   
    CALL open_new_netcdf_file(netcdf_filename=TRIM(netcdf_filename),   &
      &                       dim_list=dim_list,                  &
      &                       global_attributes=global_attributes, &
      &                       ncid=ncid)
    !-----------------------------------------------------------------

    ! rlon
    !HA debug
    CALL netcdf_put_var(ncid,lon_rot(1:cosmo_grid%nlon_rot),rlon_meta,undefined)

    ! rlat
    CALL netcdf_put_var(ncid,lat_rot(1:cosmo_grid%nlat_rot),rlat_meta,undefined)

    CALL netcdf_put_var(ncid,fr_heavy_clay(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1),     fr_heavy_clay_meta,     undefined)
    CALL netcdf_put_var(ncid,fr_silty_clay(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1),     fr_silty_clay_meta,     undefined)
    CALL netcdf_put_var(ncid,fr_light_clay(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1),     fr_light_clay_meta,     undefined)
    CALL netcdf_put_var(ncid,fr_silty_clay_loam(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1),fr_silty_clay_loam_meta,undefined)
    CALL netcdf_put_var(ncid,fr_clay_loam(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1),      fr_clay_loam_meta,      undefined)
    CALL netcdf_put_var(ncid,fr_silt(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1),           fr_silt_meta,           undefined)
    CALL netcdf_put_var(ncid,fr_silt_loam(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1),      fr_silt_loam_meta,      undefined)
    CALL netcdf_put_var(ncid,fr_sandy_clay(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1),     fr_sandy_clay_meta,     undefined)
    CALL netcdf_put_var(ncid,fr_loam(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1),           fr_loam_meta,           undefined)
    CALL netcdf_put_var(ncid,fr_sandy_clay_loam(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1),fr_sandy_clay_loam_meta,undefined)
    CALL netcdf_put_var(ncid,fr_sandy_loam(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1),     fr_sandy_loam_meta,     undefined)
    CALL netcdf_put_var(ncid,fr_loamy_sand(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1),     fr_loamy_sand_meta,     undefined)
    CALL netcdf_put_var(ncid,fr_sand(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1),           fr_sand_meta,           undefined)
    CALL netcdf_put_var(ncid,fr_undef(1:cosmo_grid%nlon_rot,1:cosmo_grid%nlat_rot,1),          fr_undef_meta,          undefined)
    

    !-----------------------------------------------------------------
    CALL netcdf_def_grid_mapping(ncid, nc_grid_def_cosmo, varid)

    CALL close_netcdf_file(ncid)


 END SUBROUTINE write_netcdf_hwsdART_cosmo_grid
!-------------------------------------------------------------
SUBROUTINE create_hwsdART_meta(p_meta,diminfo,vname,lname)
    TYPE(var_meta_info),INTENT(out) :: p_meta
    TYPE(dim_meta_info), TARGET     :: diminfo(:)
    CHARACTER(LEN=*),INTENT(in)     :: vname, &
    & lname
    
  ! local variables
  character (len=80) :: gridmp, &
  & coord

  gridmp = c_undef
  coord = c_undef
    
    p_meta%varname       =  TRIM(vname)   !< name of variable
    p_meta%n_dim         =  SIZE(diminfo) !< number of dimensions
    p_meta%diminfo       => diminfo       !< pointer to dimensions of variable
    p_meta%vartype       =  vartype_real  !< type of variable, 1 for INTEGER, 2 for REAL, 3 for CHARACTER
    p_meta%standard_name =  TRIM(vname)//'.st'   !< netcdf attribute for standard name
    p_meta%long_name     =  TRIM(lname)   !< netcdf attribute for long name
    p_meta%units         =  '1'           !< netcdf attribute for units
    p_meta%grid_mapping  =  gridmp        !< netcdf attribute grid mapping
    p_meta%coordinates   =  coord         !< netcdf attribute coordinates
    p_meta%shortName     =  TRIM(vname)//'.sh'   !< GRIB API shortName key 
END SUBROUTINE create_hwsdART_meta
!-------------------------------------------------------------
!-------------------------------------------------------------
  subroutine decode_uuid (uuid_str, uuid)
    character(len=*), intent(in)  :: uuid_str   ! uuid encoded as string
    character(len=1), intent(out) :: uuid(:)    ! decoded uuid
    character(len=2) :: buf

    integer          :: i, &
    & j, &
    & l, &
    & n, &
      b


    uuid(:) = achar (0)
    l = verify (uuid_str, "0123456789ABCDEFabcdef-")
    if (l > 0) then
       write (0,*) "Warning: invalid character in uuid: '", uuid_str(l:l),"'"
       return
    end if
    n = len  (uuid_str)
    i = 1
    j = 0
    do while (i < n)
       buf = uuid_str(i:i+1)
       if (buf(1:1) == "-") then
          i = i + 1                     ! Skip over dashes
          cycle
       end if
       i = i + 2
       read (buf,'(Z2)') b
       j = j + 1
       if (j > size (uuid)) call logging%error ("decode_uuid: uuid input too long!",__FILE__,__LINE__)
       uuid(j) = achar (b)
    end do
    if (i == n) call logging%error ("decode_uuid: uuid bad length",__FILE__,__LINE__)
  end subroutine decode_uuid
!-------------------------------------------------------------
!-------------------------------------------------------------
END Module mo_hwsdART_output_nc

