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
MODULE mo_isa_routines

  USE mo_logging
  USE mo_kind,                  ONLY: wp, i4

  USE netcdf,                   ONLY: nf90_inq_varid, &
       &                              nf90_get_var

  USE mo_io_utilities,          ONLY: check_netcdf

  USE mo_io_units,              ONLY: filename_max

  USE mo_grid_structures,       ONLY: reg_lonlat_grid


  USE mo_utilities_extpar,      ONLY: free_un ! function to get free unit number
  
  USE mo_isa_data,              ONLY: max_tiles_isa, &
       &                              ntiles_isa, &
       &                              len_isa_lon, len_isa_lat, &
       &                              isa_tiles_lat_min, &
       &                              isa_tiles_lat_max, &
       &                              isa_tiles_lon_min, &
       &                              isa_tiles_lon_max, &
       &                              isa_tiles_ncolumns,  &
       &                              isa_tiles_nrows

  USE mo_grid_structures,       ONLY: reg_lonlat_grid  ! Definition of DATA Typeto describe a regular lonlat grid

  IMPLICIT NONE

  PRIVATE


  PUBLIC :: get_dimension_isa_data, &
       &    read_namelists_extpar_isa, &
       &    get_lonlat_isa_data,    &
       &    get_isa_tiles_grid,     &
       &    det_band_isa_data,      &
       &    get_isa_data_block,     &
       &    get_isa_tile_block_indices

  CONTAINS

  !> subroutine to read namelist for orography data settings for EXTPAR 
  SUBROUTINE read_namelists_extpar_isa(namelist_file,  &
       &                               isa_type,    &
       &                               raw_data_isa_path, &
       &                               raw_data_isa_filename, &
       &                               ntiles_isa, &
       &                               isa_buffer_file      )

  
    CHARACTER (len=*), INTENT(IN) :: namelist_file !< filename with namelists for for EXTPAR settings

    CHARACTER (len=filename_max), INTENT(OUT) :: raw_data_isa_path, &        !< path to raw data
         &                                       raw_data_isa_filename(1:max_tiles_isa), & !< filename isa raw data
         &                                       isa_buffer_file !< name for isa buffer file

    INTEGER, INTENT(OUT)                      :: ntiles_isa

    INTEGER (KIND=i4)                         :: isa_type, &  !< ID of dataset used !_br 14.04.16
         &                                       nuin, & !< unit number
         &                                       ierr !< error flag

    ! NAMELIST /isa_raw_data/ raw_data_isa_path, raw_data_isa_filename, i_isa_data, ilookup_table_isa, ntiles_isa
    NAMELIST /isa_raw_data/ raw_data_isa_path, raw_data_isa_filename, ntiles_isa, isa_type !_br 14.04.16
    !> namelist with filenames for isa data output
    NAMELIST /isa_io_extpar/ isa_buffer_file

    ntiles_isa=1
    nuin = free_un()  ! functioin free_un returns free Fortran unit number
    
    OPEN(nuin,FILE=TRIM(namelist_file), IOSTAT=ierr)
    IF (ierr /= 0) THEN
      WRITE(message_text,*)'Cannot open ', TRIM(namelist_file)
      CALL logging%error(message_text,__FILE__, __LINE__) 
    ENDIF

    READ(nuin, NML=isa_raw_data, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      CALL logging%error('Cannot read in namelist isa_raw_data',__FILE__, __LINE__) 
    ENDIF
    
    READ(nuin, NML=isa_io_extpar, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      CALL logging%error('Cannot read in namelist isa_io_extpar',__FILE__, __LINE__) 
    ENDIF
     
    CLOSE(nuin)

  END SUBROUTINE read_namelists_extpar_isa

  !> inquire dimension information for isa raw data
  SUBROUTINE get_dimension_isa_data(nlon_isa, nlat_isa)

    INTEGER (KIND=i4), INTENT(OUT) :: nlon_isa, & !< number of grid elements in zonal direction for isa data
         &                            nlat_isa !< number of grid elements in meridional direction for isa data

    !local variables
    INTEGER(KIND=i4), PARAMETER    :: nx=129600, &
         &                            ny=55800

    IF(ntiles_isa == 1) THEN
      nlon_isa = len_isa_lon
      nlat_isa = len_isa_lat
    ELSE
      nlon_isa = nx
      nlat_isa = ny
    END IF

  END SUBROUTINE get_dimension_isa_data

  !> get coordinates for isa raw data
  SUBROUTINE get_lonlat_isa_data(nlon_isa, &
       &                         nlat_isa, &
       &                         lon_isa,  &
       &                         lat_isa,  &
       &                         isa_grid)


    INTEGER (KIND=i4), INTENT(IN)     :: nlon_isa, & !< number of grid elements in zonal direction for isa data
         &                               nlat_isa !< number of grid elements in meridional direction for isa data
    REAL (KIND=wp), INTENT(OUT)       :: lon_isa(1:nlon_isa), & !< longitude of isa raw data
         &                               lat_isa(1:nlat_isa) !< latitude of isa raw data

    TYPE(reg_lonlat_grid), INTENT(OUT):: isa_grid

    REAL (KIND=wp)                    :: xmin_glc, & ! area of glcover data: western longitude
         &                               xmax_glc, & ! area of glcover data: eastern longitude
         &                               ymax_glc, & ! area of glcover data: northern latitude
         &                               ymin_glc, & ! area of glcover data: southern latitude
         &                               dx_glc, &  ! grid element size of glcover data pixel in zonal direction
         &                               dy_glc  ! grid element size of glcover data pixel in meridional directionon

    INTEGER (KIND=i4)                 :: jx,jy

    xmin_glc = MINVAL(isa_tiles_lon_min)
    xmax_glc = MAXVAL(isa_tiles_lon_max)
    ymax_glc = MAXVAL(isa_tiles_lat_max)
    ymin_glc = MINVAL(isa_tiles_lat_min)
    dx_glc   = (xmax_glc - xmin_glc)/ REAL(nlon_isa,wp)
    dy_glc   = -1.0 * (ymax_glc - ymin_glc) / REAL(nlat_isa,wp)

    DO jx=1,nlon_isa
      lon_isa(jx)  = xmin_glc + 0.5*dx_glc + (jx-1)*dx_glc
    ENDDO
    DO jy=1,nlat_isa
      lat_isa(jy) = ymax_glc + 0.5*dy_glc + (jy-1)*dy_glc !note negative increment!
    ENDDO

    ! define the values for the structure isa_grid
    isa_grid%start_lon_reg = lon_isa(1)
    isa_grid%end_lon_reg   = lon_isa(nlon_isa)
    isa_grid%start_lat_reg = lat_isa(1)
    isa_grid%end_lat_reg   = lat_isa(nlat_isa)
    isa_grid%dlon_reg      = dx_glc ! (lon_isa(nlon_isa) - lon_isa(1)) / (nlon_isa - 1)
    isa_grid%dlat_reg      = dy_glc ! (lat_isa(nlat_isa) - lat_isa(1)) / (nlat_isa - 1)
    isa_grid%nlon_reg      = nlon_isa
    isa_grid%nlat_reg      = nlat_isa

  END SUBROUTINE get_lonlat_isa_data

  SUBROUTINE get_isa_tiles_grid(isa_tiles_grid)

    TYPE(reg_lonlat_grid), INTENT(OUT):: isa_tiles_grid(1:ntiles_isa)

    INTEGER(KIND=i4)                  ::  k      !counter

    REAL(KIND=wp)                     :: dlon, dlat

    DO k = 1,ntiles_isa
      dlon = (isa_tiles_lon_max(k) - isa_tiles_lon_min(k)) / REAL(isa_tiles_ncolumns(k),wp)
      ! latitude from north to south, negative increment
      dlat = -1. * (isa_tiles_lat_max(k) - isa_tiles_lat_min(k)) / REAL(isa_tiles_nrows(k),wp)   

      isa_tiles_grid(k)%start_lon_reg = isa_tiles_lon_min(k) + 0.5*dlon
      isa_tiles_grid(k)%end_lon_reg = isa_tiles_lon_max(k) - 0.5*dlon 

      isa_tiles_grid(k)%start_lat_reg = isa_tiles_lat_max(k) + 0.5*dlat
      isa_tiles_grid(k)%end_lat_reg = isa_tiles_lat_min(k) - 0.5*dlat 
      isa_tiles_grid(k)%dlon_reg = dlon
      isa_tiles_grid(k)%dlat_reg = dlat
      isa_tiles_grid(k)%nlon_reg = isa_tiles_ncolumns(k)
      isa_tiles_grid(k)%nlat_reg = isa_tiles_nrows(k)
    END DO

  END SUBROUTINE get_isa_tiles_grid

        !> determine grid description of band for ISA
  SUBROUTINE det_band_isa_data(isa_grid,start_isa_row,ta_grid)

    TYPE(reg_lonlat_grid),INTENT(IN)  :: isa_grid ! sturcture with the definition of the global data grid of the GLOBC &
                                      
    INTEGER (KIND=i4),INTENT(IN)      :: start_isa_row         ! number of the start row of band 

    TYPE(reg_lonlat_grid), INTENT(OUT):: ta_grid       ! structure with definition of the target area grid.

    INTEGER(KIND=i4)                  :: nrows = 1000              ! number of rows, set to 1000 as default

    ta_grid%dlon_reg = isa_grid%dlon_reg
    ta_grid%dlat_reg = isa_grid%dlat_reg

    ta_grid%start_lon_reg = isa_grid%start_lon_reg
    ta_grid%end_lon_reg = isa_grid%end_lon_reg
    ta_grid%nlon_reg = isa_grid%nlon_reg

    !latitude from north to south, negative increment
    ta_grid%nlat_reg = nrows
    ta_grid%start_lat_reg = isa_grid%start_lat_reg + ta_grid%dlat_reg * (start_isa_row-1)   ! latitude from nor &
!& touth, note the negative increment!
    ta_grid%end_lat_reg = ta_grid%start_lat_reg + ta_grid%dlat_reg * (nrows - 1)   ! latitude from north to south, note the &
!&  ve increment!
     
    ! check for the southern bound of the isa data
    IF (ta_grid%end_lat_reg < isa_grid%end_lat_reg) THEN ! band is at the southern bound
      ta_grid%end_lat_reg = isa_grid%end_lat_reg
      ta_grid%nlat_reg    = NINT(((ta_grid%end_lat_reg - ta_grid%start_lat_reg) / ta_grid%dlat_reg)) + 1
    END IF

  END SUBROUTINE det_band_isa_data
      
  ! get ISA data block for a given target area from the tile block indices
  SUBROUTINE get_isa_data_block(ta_grid,              &
       &                        isa_tiles_grid, &
       &                        ncids_isa,      &
       &                        isa_block)

    TYPE(reg_lonlat_grid), INTENT(IN)  :: ta_grid, &
         &                                isa_tiles_grid(1:ntiles_isa)

    INTEGER(KIND=i4) , INTENT(IN)      :: ncids_isa(1:ntiles_isa)

    REAL (KIND=wp), INTENT(OUT)        :: isa_block(1:ta_grid%nlon_reg,1:ta_grid%nlat_reg)

    !local variables
    INTEGER (KIND=i4)                  :: isa_startrow(1:ntiles_isa), & !< startrow indices for each ISA tile
         &                                isa_endrow(1:ntiles_isa), & !< endrow indices for each ISA tile
         &                                isa_startcolumn(1:ntiles_isa), & !< starcolumn indices for each ISA tile
         &                                isa_endcolumn(1:ntiles_isa), & !< endcolumn indices for each ISA tile
         &                                ta_start_ie(1:ntiles_isa), &    !< indices of target area block for first column of each ISA &
         &                                ta_end_ie(1:ntiles_isa), &      !< indices of target area block for last column of each ISA  &
         &                                ta_start_je(1:ntiles_isa), &  !< indices of target area block for first row of each ISA tile &
         &                                ta_end_je(1:ntiles_isa), &   !< indices of target area block for last row of each ISA tile
         &                                varid, &               !< id of variable
         &                                nrows, & !< number of rows ! dimensions for raw_isa_block
         &                                ncolumns, & !< number of columns ! dimensions for raw_isa_block
         &                                k, &     ! counter
         &                                errorcode !< error status variable

    INTEGER (KIND=i4), ALLOCATABLE     :: raw_isa_block(:,:) !< a block with ISA data

    CHARACTER (LEN=80)                 :: varname  !< name of variable

    varname = 'ISA'   ! I know that in the ISA netcdf files the isa data is stored in a variable "ISA"

    CALL get_isa_tile_block_indices(ta_grid,              &
         &                                isa_tiles_grid, &  
         &                                isa_startrow,   &
         &                                isa_endrow,     & 
         &                                isa_startcolumn,&
         &                                isa_endcolumn,  &
         &                                ta_start_ie,          & 
         &                                ta_end_ie,            &
         &                                ta_start_je,          &
         &                                ta_end_je)


    DO k = 1, ntiles_isa
      IF ((isa_startrow(k)/=0).AND.(isa_startcolumn(k)/=0)) THEN
        nrows = isa_endrow(k) - isa_startrow(k) + 1
        ncolumns = isa_endcolumn(k) - isa_startcolumn(k) + 1
 
        ALLOCATE (raw_isa_block(1:ncolumns,1:nrows), STAT=errorcode)
        IF(errorcode/=0) CALL logging%error('Cant allocate the array raw_isa_block',__FILE__,__LINE__)

        CALL check_netcdf(nf90_inq_varid(ncids_isa(k),TRIM(varname),varid)) ! get the varid of the altitude variable
        ! get the data into the raw_isa_block
        CALL check_netcdf(nf90_get_var(ncids_isa(k), varid,  raw_isa_block,     & 
             &  start=(/isa_startcolumn(k),isa_startrow(k)/),count=(/ncolumns,nrows/)))

         isa_block(ta_start_ie(k):ta_end_ie(k),ta_start_je(k):ta_end_je(k)) = raw_isa_block(1:ncolumns,1:nrows)
         
         DEALLOCATE (raw_isa_block, STAT=errorcode)
         IF(errorcode/=0) CALL logging%error('Cant deallocate the array raw_isa_block',__FILE__,__LINE__)
       
      ENDIF
    ENDDO

  END SUBROUTINE get_isa_data_block

       ! get startrow, endrow, startcolumn and endcolumn of each ISA tile (raw data) for a 
       ! given target area (ta_grid) and get start_indices (lon, lat) and end_indices of the target
       ! area for each ISA tile
       ! The ISA raw data are split in 6 tiles, so the target area may overlap several tiles.
       ! This subroutine determines the necessary indices to read in the ISA data into the
       ! target area.
       ! ISA tiles which are outside the target block will get indices with the value '0'
  SUBROUTINE get_isa_tile_block_indices(ta_grid,              &
       &                                isa_tiles_grid, &
       &                                isa_startrow,   &
       &                                isa_endrow,     & 
       &                                isa_startcolumn,&
       &                                isa_endcolumn,  &
       &                                ta_start_ie,          &
       &                                ta_end_ie,            &
       &                                ta_start_je,          &
       &                                ta_end_je)

    TYPE(reg_lonlat_grid), INTENT(IN) :: ta_grid, & !< structure with definition of the target area grid (dlon must be the same a &
         &                               isa_tiles_grid(1:ntiles_isa) !< structure with defenition of the raw data &

    INTEGER (KIND=i4), INTENT(OUT)    :: isa_startrow(1:ntiles_isa), &    !< startrow indices for each ISA tile
         &                               isa_endrow(1:ntiles_isa), &      !< endrow indices for each ISA tile
         &                               isa_startcolumn(1:ntiles_isa), &  !< starcolumn indices for each ISA tile &
         &                               isa_endcolumn(1:ntiles_isa), &   !< endcolumn indices for each ISA tile
         &                               ta_start_ie(1:ntiles_isa), &    !< indices of target area block for first column of e &
         &                               ta_end_ie(1:ntiles_isa), &      !< indices of target area block for last column of ea &
         &                               ta_start_je(1:ntiles_isa), &  !< indices of target area block for first row of each G &
         &                               ta_end_je(1:ntiles_isa)   !< indices of target area block for last row of each GLO &

    ! local variables
    INTEGER (KIND=i4)                 :: startrow, & ! startrow for tile
                                         endrow, & 
                                         startcolumn, &
                                         endcolumn, k, undefined

    REAL (KIND=wp)                    :: dlon, dlat

    undefined = 0
    isa_startrow     = undefined
    isa_endrow       = undefined
    isa_startcolumn  = undefined
    isa_endcolumn    = undefined
    ta_start_ie = undefined 
    ta_end_ie   = undefined
    ta_start_je = undefined
    ta_end_je   = undefined

    k=1                      ! determin dlon and dlat (are the same for all tiles)
    dlon = ta_grid%dlon_reg
    dlat = ta_grid%dlat_reg

    DO k = 1,ntiles_isa   !loop over the tiles which overlap the target area
      startcolumn = NINT((ta_grid%start_lon_reg - isa_tiles_grid(k)%start_lon_reg)/dlon) +1
      IF (startcolumn < 1) THEN 
        isa_startcolumn(k) = 1
        ! get the start index of the subtile for the target area block
        ta_start_ie(k) = NINT ((isa_tiles_grid(k)%start_lon_reg - ta_grid%start_lon_reg)/dlon) + 1
      ELSE IF (startcolumn > isa_tiles_ncolumns(k)) THEN
        isa_startcolumn(k) = 0
        ta_start_ie(k) = 0
      ELSE
        isa_startcolumn(k) = startcolumn
        ta_start_ie(k) = 1
      ENDIF

      ! get endcolumn for tile k
      endcolumn = NINT((ta_grid%end_lon_reg - isa_tiles_grid(k)%start_lon_reg)/dlon) +1

      IF (endcolumn > isa_tiles_ncolumns(k)) THEN 
        isa_endcolumn(k) = isa_tiles_ncolumns(k)
        ! get the end index of the subtile for the target area block
        ta_end_ie(k) = NINT ((isa_tiles_grid(k)%end_lon_reg - ta_grid%start_lon_reg)/dlon) + 1
      ELSE IF (endcolumn < 1) THEN
        isa_endcolumn(k) = 0
        ta_end_ie(k) = 0
      ELSE
        isa_endcolumn(k) = endcolumn
        ta_end_ie(k) = ta_grid%nlon_reg
      ENDIF

      ! get startrow for tile k
      startrow = NINT((ta_grid%start_lat_reg - isa_tiles_grid(k)%start_lat_reg)/dlat) + 1
     
      IF (startrow < 1) THEN 
        isa_startrow(k) = 1
        ! get the start index of the subtile for the target area block
        ta_start_je(k) = NINT ((isa_tiles_grid(k)%start_lat_reg  - ta_grid%start_lat_reg)/dlat) + 1

      ELSE IF (startrow > isa_tiles_nrows(k)) THEN
        isa_startrow(k) = 0
        ta_start_je(k) = 0
      ELSE
        isa_startrow(k) = startrow
        ta_start_je(k) = 1
      ENDIF

      ! get endrow for tile k
      endrow   = NINT(( ta_grid%end_lat_reg - isa_tiles_grid(k)%start_lat_reg )/dlat)  + 1
     
      IF (endrow > isa_tiles_nrows(k)) THEN 
        isa_endrow(k) = isa_tiles_nrows(k)
        ! get the start index of the subtile for the target area block
        ta_end_je(k) = NINT ((isa_tiles_grid(k)%end_lat_reg -  ta_grid%start_lat_reg )/dlat) + 1

      ELSE IF (endrow < 1) THEN
        isa_endrow(k) = 0
        ta_end_je(k) = 0
      ELSE
        isa_endrow(k) = endrow
        ta_end_je(k) =  ta_grid%nlat_reg
      ENDIF
 
    ENDDO  ! loop over the tiles 

  END SUBROUTINE get_isa_tile_block_indices

END MODULE mo_isa_routines
