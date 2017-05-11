!+ Fortran module with routines and settings for DEM subgrid-scale slope data
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V4_0         2016/07/28 Daniel Luethi
!  Initial release
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!> Fortran module with routines and settings for DEM subgrid-scale slope data
!> \author Daniel Luethi
!>
MODULE mo_sgsl_routines

!> kind parameters are defined in MODULE data_parameters
USE mo_kind, ONLY: wp, &
                   i8, &
                   i4, &
                   ishort

USE netcdf,      ONLY :   &
  nf90_open,              &
  nf90_close,             &
  nf90_inquire,           &
  nf90_inquire_dimension, &
  nf90_inquire_variable,  &
  nf90_inq_attname,       &
  nf90_inquire_attribute, &
  nf90_get_att,           &
  nf90_inquire_dimension, &
  nf90_inq_varid,         &
  nf90_get_var,           &
  nf90_noerr,             &
  nf90_strerror

USE netcdf,      ONLY:     &
  nf90_create,             &
  nf90_def_dim,            &
  nf90_def_var,            &
  nf90_enddef,             &
  nf90_redef,              &
  nf90_put_att,            &
  nf90_put_var

 
USE netcdf,      ONLY :    &
  NF90_CHAR,               &
  NF90_DOUBLE,             &
  NF90_FLOAT,              &
  NF90_INT,                &
  NF90_BYTE,               &
  NF90_SHORT


USE netcdf,      ONLY :    &
  NF90_GLOBAL,             &
  NF90_UNLIMITED,          &
  NF90_CLOBBER,            &
  NF90_NOWRITE


!> abort_extpar defined in MODULE utilities_extpar
USE mo_utilities_extpar, ONLY: abort_extpar

USE mo_grid_structures,  ONLY: reg_lonlat_grid
USE mo_base_geometry,    ONLY: geographical_coordinates

USE mo_io_units,         ONLY: filename_max


USE mo_io_utilities,     ONLY: check_netcdf

IMPLICIT NONE

PRIVATE

PUBLIC :: read_namelists_extpar_sg_slope,&
          det_sgsl_tiles_grid,            &
          det_sgsl_grid,                  &
          get_sgsl_tile_block_indices,    &
          open_netcdf_sgsl_tile,          &
          close_netcdf_sgsl_tile

PUBLIC :: det_band_gd, get_sgsl_data_block

CONTAINS

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!> subroutine to read namelist for orography data settings for EXTPAR 
SUBROUTINE read_namelists_extpar_sg_slope (namelist_file,          &
                                           raw_data_sgsl_path,&
                                           sgsl_files,             &  !mes>
                                           ntiles_column,          &
                                           ntiles_row,             &
                                           idem_type,             &  
                                           sgsl_buffer_file)

  USE mo_utilities_extpar, ONLY: free_un ! function to get free unit number
  USE mo_sgsl_data,        ONLY: max_tiles  

  
  CHARACTER (LEN=filename_max), INTENT(IN) :: namelist_file !< filename with namelists for for EXTPAR settings
  ! orography

  
  CHARACTER (LEN=filename_max), INTENT(OUT) :: raw_data_sgsl_path        !< path to raw data
  CHARACTER (LEN=filename_max), INTENT(OUT) :: sgsl_files(1:max_tiles)                     !< filenames globe raw data
  INTEGER (KIND=i4), INTENT(OUT)  :: ntiles_column     !< number of tile columns
  INTEGER (KIND=i4), INTENT(OUT)  :: ntiles_row        !< number of tile rows
  INTEGER (KIND=i4), INTENT(OUT)  :: idem_type
  CHARACTER (len=filename_max), INTENT(OUT) :: sgsl_buffer_file !< name for subgrid slope buffer file

   INTEGER           :: nuin !< unit number
   INTEGER (KIND=i4) :: ierr !< error flag
   INTEGER :: i, nzylen

  !> namelist with filenames for orography data output
  NAMELIST /sgsl_io_extpar/ sgsl_buffer_file
  !> namelist with information on orography data input
  NAMELIST /sgsl_raw_data/ idem_type, raw_data_sgsl_path, ntiles_column, ntiles_row, sgsl_files
    

   nuin = free_un()  ! function free_un returns free Fortran unit number
   OPEN(nuin,FILE=TRIM(namelist_file), IOSTAT=ierr)
   READ(nuin, NML=sgsl_io_extpar, IOSTAT=ierr)
   READ(nuin, NML=sgsl_raw_data, IOSTAT=ierr)

   CLOSE(nuin, IOSTAT=ierr)

   nzylen=LEN_TRIM(raw_data_sgsl_path)
   IF( nzylen > 0 ) THEN
     IF( raw_data_sgsl_path(nzylen:nzylen) /= '/') THEN
       IF( nzylen < LEN(raw_data_sgsl_path) ) THEN
         raw_data_sgsl_path = raw_data_sgsl_path (1:nzylen)//'/'
       ENDIF
     ENDIF
   ENDIF



END SUBROUTINE read_namelists_extpar_sg_slope

!> determine GLOBE raw data grid
!> \author Hermann Asensio
SUBROUTINE det_sgsl_tiles_grid(sgsl_tiles_grid)
USE mo_sgsl_data, ONLY : ntiles , &    !< GLOBE raw data has 16 tiles and ASTER has 13
                         idem_type,    &
                         dem_aster,    &
                         dem_gl,       &
                         tiles_lon_min, &
                         tiles_lon_max, &
                         tiles_lat_min, &
                         tiles_lat_max, &
                         tiles_ncolumns, &
                         tiles_nrows

  TYPE(reg_lonlat_grid), INTENT(OUT) :: sgsl_tiles_grid(1:ntiles) 
    !< structure with definition of the raw data grid for the input tiles

  INTEGER :: k ! counter

  REAL (KIND=wp) :: lon0
  REAL (KIND=wp) :: lat0

  REAL (KIND=wp) :: dlon
  REAL (KIND=wp) :: dlat

  DO k=1,  ntiles ! determine the globe_tile_grid information from the namelist information

    dlon = (tiles_lon_max(k) - tiles_lon_min(k)) / FLOAT(tiles_ncolumns(k))
           
    dlat = -1. * (tiles_lat_max(k) - tiles_lat_min(k)) / FLOAT(tiles_nrows(k))

    ! latitude from north to south, negative increment
            

    sgsl_tiles_grid(k)%start_lon_reg  = tiles_lon_min(k) + 0.5 * dlon
    sgsl_tiles_grid(k)%end_lon_reg    = tiles_lon_max(k) - 0.5 * dlon
           
    sgsl_tiles_grid(k)%start_lat_reg  = tiles_lat_max(k) + 0.5 * dlat 
    ! latitude from north to south, note the negative increment!
    sgsl_tiles_grid(k)%end_lat_reg    = tiles_lat_min(k) - 0.5 * dlat 
    ! latitude from north to south, note the negative increment!

    sgsl_tiles_grid(k)%dlon_reg = dlon
    sgsl_tiles_grid(k)%dlat_reg = dlat

    sgsl_tiles_grid(k)%nlon_reg = tiles_ncolumns(k)
    sgsl_tiles_grid(k)%nlat_reg = tiles_nrows(k)

  ENDDO
            

END SUBROUTINE det_sgsl_tiles_grid



!> determine complete(global) GLOBE raw data grid 
!> \author Hermann Asensio
SUBROUTINE det_sgsl_grid(sgsl_grid)
  USE mo_sgsl_data, ONLY :  nc_tot,        &      
   &                        nr_tot,        &
   &                        idem_type,     &
   &                        dem_aster,     &
   &                        dem_gl,        &
   &                        demraw_lat_min,&
   &                        demraw_lat_max,&
   &                        demraw_lon_min,&
   &                        demraw_lon_max

  TYPE(reg_lonlat_grid), INTENT(OUT) :: sgsl_grid !< structure with definition of the global data grid of the GLOBE data 
  INTEGER :: k ! counter

  REAL (KIND=wp) :: lon0
  REAL (KIND=wp) :: lat0
  REAL (KIND=wp) :: dlon
  REAL (KIND=wp) :: dlat

  !mes > as ASTER does not cover the whole globe until now different procedures must be chosen for ASTER and GLOBE
  SELECT CASE(idem_type)
    CASE(dem_aster, dem_gl)

      dlon = (demraw_lon_max - demraw_lon_min) / FLOAT(nc_tot)

      dlat = -1. * (demraw_lat_max - demraw_lat_min) / FLOAT(nr_tot)
      ! latitude from north to south, negative increment

      sgsl_grid%start_lon_reg  =  demraw_lon_min + 0.5 * dlon
      sgsl_grid%end_lon_reg    =  demraw_lon_max - 0.5 * dlon

            
      sgsl_grid%start_lat_reg = demraw_lat_max + 0.5 * dlat 
      ! latitude from north to south, note the negative increment!
      sgsl_grid%end_lat_reg  =  demraw_lat_min - 0.5 * dlat
      ! latitude from north to south, note the negative increment!


  END SELECT

  sgsl_grid%dlon_reg = dlon
  sgsl_grid%dlat_reg = dlat

  sgsl_grid%nlon_reg = nc_tot
  sgsl_grid%nlat_reg = nr_tot


END SUBROUTINE det_sgsl_grid

!> determine grid description of band for GLOBE I/O 
!> \author Hermann Asensio
SUBROUTINE det_band_gd(sgsl_grid,start_sgsl_row, ta_grid)

  TYPE(reg_lonlat_grid), INTENT(IN) :: sgsl_grid !< structure with definition of the global data grid of the GLOBE data 
  INTEGER, INTENT(IN) :: start_sgsl_row !< number of the start row of band of sgsl_grid (global domain)
  TYPE(reg_lonlat_grid), INTENT(OUT) :: ta_grid !< structure with defenition of the target area grid

  INTEGER  :: nrows = 500 !< number of rows, set to 500 as default
  ! band from east to west for the whole globe, like the complete sgsl_grid

  ta_grid%dlon_reg = sgsl_grid%dlon_reg
  ta_grid%dlat_reg = sgsl_grid%dlat_reg

  ta_grid%start_lon_reg = sgsl_grid%start_lon_reg
  ta_grid%end_lon_reg =  sgsl_grid%end_lon_reg
  ta_grid%nlon_reg = sgsl_grid%nlon_reg
          
  ! latitude from north to south, negative increment
  ta_grid%nlat_reg = nrows
  ta_grid%start_lat_reg = sgsl_grid%start_lat_reg + ta_grid%dlat_reg * (start_sgsl_row - 1)  
  ! latitude from north to south, note the negative increment!
  ta_grid%end_lat_reg  =  ta_grid%start_lat_reg + ta_grid%dlat_reg * (nrows - 1)
  ! latitude from north to south, note the negative increment!
  ! check for south pole
  IF (ta_grid%end_lat_reg < sgsl_grid%end_lat_reg) THEN ! band is at south pole
    ta_grid%end_lat_reg =  sgsl_grid%end_lat_reg
    ta_grid%nlat_reg =  NINT(((ta_grid%end_lat_reg - ta_grid%start_lat_reg) / ta_grid%dlat_reg)) + 1
  ENDIF

END SUBROUTINE det_band_gd

!----------------------------------------------------------------------------------------------------------------

       !> get startrow, endrow, startcolumn and endcolumn of each GLOBE tile (raw data) for a 
       !! given target area (ta_grid) and
       !! get start_indices (lon, lat) and end_indices of the target area for each GLOBE tile
       !! The GLOBE raw data are split in 16 tiles (ASTER in 36), so the target area may overlap several tiles.
       !! This subroutine determines the necesarry indices to read in the GLOBE/ASTER data into the
       !! target area.
       !! GLOBE/ASTER tiles which are outside the target block will get indices with the value '0'

       SUBROUTINE get_sgsl_tile_block_indices(ta_grid,         &
         &                                     sgsl_tiles_grid, &
         &                                     sgsl_startrow,  &
         &                                     sgsl_endrow,    & 
         &                                     sgsl_startcolumn,&
         &                                     sgsl_endcolumn, &
         &                                     ta_start_ie, &
         &                                     ta_end_ie,   &
         &                                     ta_start_je, &
         &                                     ta_end_je)

        USE mo_sgsl_data, ONLY : ntiles ,     &    !< GLOBE raw data has 16 tiles, ASTER has 36
                                idem_type,    &
                                dem_aster,    &
                                dem_gl,       &
                                tiles_lon_min, &
                                tiles_lon_max, &
                                tiles_lat_min, &
                                tiles_lat_max, &
                                tiles_ncolumns,&
                                tiles_nrows

       USE mo_grid_structures, ONLY: reg_lonlat_grid  !< Definition of Data Type to describe a regular (lonlat) grid

       TYPE(reg_lonlat_grid), INTENT(IN)  :: ta_grid !< structure with definition of the target area grid (dlon must be the same as for the whole GLOBE dataset)


       TYPE(reg_lonlat_grid), INTENT(IN) :: sgsl_tiles_grid(1:ntiles) 
                                            !< structure with defenition of the raw data grid for the 16 GLOBE tiles

       INTEGER (KIND=i4), INTENT(OUT) :: sgsl_startrow(1:ntiles)    !< startrow indices for each GLOBE tile
       INTEGER (KIND=i4), INTENT(OUT) :: sgsl_endrow(1:ntiles)      !< endrow indices for each GLOBE tile

       INTEGER (KIND=i4), INTENT(OUT) :: sgsl_startcolumn(1:ntiles)  !< starcolumn indices for each GLOBE tile
       INTEGER (KIND=i4), INTENT(OUT) :: sgsl_endcolumn(1:ntiles)   !< endcolumn indices for each GLOBE tile

       INTEGER (KIND=i4), INTENT(OUT) :: ta_start_ie(1:ntiles)    
                                         !< indices of target area block for first column of each GLOBE tile
       INTEGER (KIND=i4), INTENT(OUT) :: ta_end_ie(1:ntiles)      
                                         !< indices of target area block for last column of each GLOBE tile
       INTEGER (KIND=i4), INTENT(OUT) :: ta_start_je(1:ntiles)  
                                         !< indices of target area block for first row of each GLOBE tile
       INTEGER (KIND=i4), INTENT(OUT) :: ta_end_je(1:ntiles)   
                                         !< indices of target area block for last row of each GLOBE tile

       
       INTEGER (KIND=i4) :: index_k !< index of GLOBE tile which contains point_geo

       ! local variables

       INTEGER  :: i          ! index for tiles (i,j,m,n,o)
       INTEGER  :: j 
       INTEGER  :: m
       INTEGER  :: n
       INTEGER  :: o
       INTEGER  :: t_i_start 
       INTEGER  :: t_i_end
       INTEGER  :: t_j_start
       INTEGER  :: t_j_end

       REAL  :: lon0_t ! startlon for dummy grid
       REAL  :: lat0_t ! startlat for dummy grid
       REAL  :: dlon_t ! dlon for dummy grid
       REAL  :: dlat_t ! dlat for dummy grid

       INTEGER  :: undefined

       REAL (KIND=wp) :: point_lon_coor

       REAL (KIND=wp) :: tb_ll_lon ! longitude coordinate for lower left corner of target block
       REAL (KIND=wp) :: tb_ll_lat ! longitude coordinate for lower left corner of target block

       REAL (KIND=wp) :: tb_ur_lon ! longitude coordinate for upper right corner of target block
       REAL (KIND=wp) :: tb_ur_lat ! longitude coordinate for upper right corner of target block


       INTEGER (KIND=i4) :: startrow ! startrow for tile
       INTEGER (KIND=i4) :: endrow 
       INTEGER (KIND=i4) :: startcolumn
       INTEGER (KIND=i4) :: endcolumn

       REAL (KIND=wp) :: dlon
       REAL (KIND=wp) :: dlat

       REAL (KIND=wp) :: stile_ll_lon ! longitude coordinate for lower left corner of subtile
       REAL (KIND=wp) :: stile_ll_lat ! latitued coordinate for lower left corner of subtile

       REAL (KIND=wp) :: stile_ur_lon ! longitude coordinate for upper right corner of subtile
       REAL (KIND=wp) :: stile_ur_lat ! latitude coordinate for upper right corner of subtile



       INTEGER :: k

       undefined = 0
       sgsl_startrow     = undefined
       sgsl_endrow       = undefined
       sgsl_startcolumn   = undefined
       sgsl_endcolumn    = undefined
       ta_start_ie = undefined 
       ta_end_ie = undefined
       ta_start_je = undefined
       ta_end_je = undefined

       k=1 ! determin dlon and dlat (are the same for all tiles)
       dlon = ta_grid%dlon_reg
       dlat = ta_grid%dlat_reg
       !dlon = (tiles_lon_max(k) - tiles_lon_min(k)) / FLOAT(tiles_ncolumns(k))
       !dlat =(tiles_lat_max(k) - tiles_lat_min(k)) / FLOAT(tiles_nrows(k))

       ! the GLOBE data are diveded in 16 tiles, 
       ! this defines a "dummy grid" to determine the tile index with a function
       ! lon from -180 to 180 with dlon 90 degrees
       ! lat from 100 to -100 with dlat 50 degrees
       lon0_t = -180. 
       lat0_t = 100.
       dlon_t = 90.
       dlat_t = 50.

!       SELECT CASE(idem_type)
!       CASE(dem_aster, dem_gl)
!        m = 1
!        n = 1
!        o = ntiles
!       END SELECT
       

        DO k = 1,ntiles
       
         ! get startcolumn for tile k
         startcolumn = NINT((ta_grid%start_lon_reg - sgsl_tiles_grid(k)%start_lon_reg)/dlon) + 1 
         !< here I want nearest index (NINT)
         IF (startcolumn < 1) THEN 
           sgsl_startcolumn(k) = 1
           ! get the start index of the subtile for the target area block
           ta_start_ie(k) = NINT ((sgsl_tiles_grid(k)%start_lon_reg - ta_grid%start_lon_reg)/dlon) + 1 
           !< index of target area block

         ELSE IF (startcolumn > tiles_ncolumns(k)) THEN
           sgsl_startcolumn(k) = 0
           ta_start_ie(k) = 0
         ELSE
           sgsl_startcolumn(k) = startcolumn
           ta_start_ie(k) = 1
         ENDIF

         ! get endcolumn for tile k
         endcolumn = NINT((ta_grid%end_lon_reg - sgsl_tiles_grid(k)%start_lon_reg)/dlon) + 1
         IF (endcolumn > tiles_ncolumns(k)) THEN 
           sgsl_endcolumn(k) = tiles_ncolumns(k)
           ! get the end index of the subtile for the target area block
           stile_ur_lon =  sgsl_tiles_grid(k)%end_lon_reg ! coordinates [degrees]
           ta_end_ie(k) = NINT ((sgsl_tiles_grid(k)%end_lon_reg - ta_grid%start_lon_reg)/dlon) + 1 
           !< index of target area block
         ELSE IF (endcolumn < 1) THEN
           sgsl_endcolumn(k) = 0
           ta_end_ie(k) = 0
         ELSE
           sgsl_endcolumn(k) = endcolumn
           ta_end_ie(k) = ta_grid%nlon_reg
         ENDIF


         ! get startrow for tile k
         startrow = NINT((ta_grid%start_lat_reg - sgsl_tiles_grid(k)%start_lat_reg)/dlat) + 1
        
         IF (startrow < 1) THEN 
           sgsl_startrow(k) = 1
           ! get the start index of the subtile for the target area block
           ta_start_je(k) = NINT ((sgsl_tiles_grid(k)%start_lat_reg  - ta_grid%start_lat_reg)/dlat) + 1 
           !< index of target area block
            
          ELSE IF (startrow > tiles_nrows(k)) THEN
           sgsl_startrow(k) = 0
           ta_start_je(k) = 0
         ELSE
           sgsl_startrow(k) = startrow
           ta_start_je(k) = 1
         ENDIF


         ! get endrow for tile k
         endrow   = NINT(( ta_grid%end_lat_reg - sgsl_tiles_grid(k)%start_lat_reg )/dlat)  + 1
        
         IF (endrow > tiles_nrows(k)) THEN 
           sgsl_endrow(k) = tiles_nrows(k)
           ! get the start index of the subtile for the target area block
           ta_end_je(k) = NINT ((sgsl_tiles_grid(k)%end_lat_reg -  ta_grid%start_lat_reg )/dlat) + 1 
           !< index of target area block

         ELSE IF (endrow < 1) THEN
           sgsl_endrow(k) = 0
           ta_end_je(k) = 0
         ELSE
           sgsl_endrow(k) = endrow
           ta_end_je(k) =  ta_grid%nlat_reg
         ENDIF
 

       ENDDO


       END SUBROUTINE get_sgsl_tile_block_indices


!----------------------------------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------------------------------
    !> open netcdf-file and get netcdf unit file number
       SUBROUTINE open_netcdf_sgsl_tile(path_sgsl_tile, &
                                        ncid)
         CHARACTER (len=*), INTENT(in) :: path_sgsl_tile         !< filename with path to GLOBE tile
         INTEGER, INTENT(out) :: ncid                             !< netcdf unit file number

        !! open netcdf file 
            call check_netcdf( nf90_open(TRIM(path_sgsl_tile),NF90_NOWRITE, ncid))

       END SUBROUTINE open_netcdf_sgsl_tile

        !> close netcdf-file 
       SUBROUTINE close_netcdf_sgsl_tile(ncid)
         INTEGER, INTENT(in) :: ncid                             !< netcdf unit file number

        !! close netcdf file 
            call check_netcdf( nf90_close( ncid))

       END SUBROUTINE close_netcdf_sgsl_tile

!-------------------------------------------------------------------------------   !> get GLOBE data block for a given target area from the tile block indices
       SUBROUTINE get_sgsl_data_block(sgsl_file_1,     &   !mes ><
         &                             ta_grid,         &
         &                             sgsl_tiles_grid, &
         &                             ncids_sgsl, &
         &                             sgsl_block)

       USE mo_grid_structures, ONLY: reg_lonlat_grid  !< Definition of Data Type to describe a regular (lonlat) grid

       USE mo_sgsl_data, ONLY : ntiles  !< there are 16 GLOBE tiles 
       USE mo_sgsl_data, ONLY : nc_tot     !< total number of columns in GLOBE data: 43200
       USE mo_sgsl_data, ONLY : nc_tile    !< number of columns in a GLOBE tile
       ! mes >
       USE mo_sgsl_data, ONLY : get_varname   ! gets the variable name of the elevation 
       USE mo_sgsl_data, ONLY : idem_type
       USE mo_sgsl_data, ONLY : dem_aster
       USE mo_sgsl_data, ONLY : dem_gl

       CHARACTER (len=filename_max), INTENT(IN)     :: sgsl_file_1
       ! mes <

       TYPE(reg_lonlat_grid), INTENT(IN)  :: ta_grid 
       !< structure with definition of the target area grid (dlon must be the same as for the whole GLOBE dataset)
       TYPE(reg_lonlat_grid), INTENT(IN) :: sgsl_tiles_grid(1:ntiles) 
       !< structure with defenition of the raw data grid for the 16 GLOBE tiles

       INTEGER , INTENT(IN) :: ncids_sgsl(1:ntiles)  
       !< ncid for the GLOBE tiles, the netcdf files have to be opened by a previous call of open_netcdf_GLOBE_tile

       REAL (KIND=wp), INTENT(OUT) :: sgsl_block(1:ta_grid%nlon_reg,1:ta_grid%nlat_reg) !< a block of GLOBE altitude data 
       !local variables

       INTEGER (KIND=i4) :: sgsl_startrow(1:ntiles) !< startrow indices for each GLOBE tile
       INTEGER (KIND=i4) :: sgsl_endrow(1:ntiles) !< endrow indices for each GLOBE tile
       INTEGER (KIND=i4) :: sgsl_startcolumn(1:ntiles) !< starcolumn indices for each GLOBE tile
       INTEGER (KIND=i4) :: sgsl_endcolumn(1:ntiles) !< endcolumn indices for each GLOBE tile

       INTEGER (KIND=i4) :: ta_start_ie(1:ntiles)    !< indices of target area block for first column of each GLOBE tile
       INTEGER (KIND=i4) :: ta_end_ie(1:ntiles)      !< indices of target area block for last column of each GLOBE tile
       INTEGER (KIND=i4) :: ta_start_je(1:ntiles)  !< indices of target area block for first row of each GLOBE tile
       INTEGER (KIND=i4) :: ta_end_je(1:ntiles)   !< indices of target area block for last row of each GLOBE tile


       INTEGER (KIND=ishort), ALLOCATABLE :: raw_sgsl_block(:,:) !< a block with GLOBE data
       REAL (KIND=wp) :: scale_factor

       INTEGER :: varid               !< id of variable
       CHARACTER (LEN=80) :: varname  !< name of variable

       INTEGER :: nrows !< number of rows ! dimensions for raw_sgsl_block
       INTEGER :: ncolumns !< number of columns ! dimensions for raw_sgsl_block



       INTEGER :: k ! counter
       INTEGER :: errorcode !< error status variable
       
! mes >
       
       CALL get_varname(sgsl_file_1,varname)
       !print*, TRIM(varname)

!       varname = 'S_ORO'  
        ! I know that in the S_ORO is used for the subgrid-scale slope data
! mes <
       
       
       CALL get_sgsl_tile_block_indices( ta_grid,         &
         &                                sgsl_tiles_grid,& 
         &                                sgsl_startrow,  &
         &                                sgsl_endrow,    & 
         &                                sgsl_startcolumn,&
         &                                sgsl_endcolumn,  &
         &                                ta_start_ie, &
         &                                ta_end_ie,   &
         &                                ta_start_je, &
         &                                ta_end_je)

        
          DO k=1,ntiles
           IF ((sgsl_startrow(k)/=0).AND.(sgsl_startcolumn(k)/=0)) THEN
             nrows = sgsl_endrow(k) - sgsl_startrow(k) + 1
             ncolumns = sgsl_endcolumn(k) - sgsl_startcolumn(k) + 1
 !            print*, k      
 !            print*, sgsl_startrow(k)
 !            print*, sgsl_endrow(k) 
 !            print*, k      
 !            print*, sgsl_startcolumn(k)
 !            print*, sgsl_endcolumn(k)
           
             ALLOCATE (raw_sgsl_block(1:ncolumns,1:nrows), STAT=errorcode)
             IF(errorcode/=0) CALL abort_extpar('Cant allocate the array raw_sgsl_block')
             CALL check_netcdf(nf90_inq_varid(ncids_sgsl(k),TRIM(varname),varid)) ! get the varid of the altitude variable
             ! get the data into the raw_sgsl_block
             CALL check_netcdf(nf90_get_var(ncids_sgsl(k),varid,raw_sgsl_block,     & 
                  &     start=(/sgsl_startcolumn(k),sgsl_startrow(k)/),count=(/ncolumns,nrows/)))
            

                sgsl_block(ta_start_ie(k):ta_end_ie(k),ta_start_je(k):ta_end_je(k)) = raw_sgsl_block(1:ncolumns,1:nrows)
             CALL check_netcdf(nf90_get_att(ncids_sgsl(k),varid, & 
                                            'scale_factor', scale_factor))

                sgsl_block(ta_start_ie(k):ta_end_ie(k),ta_start_je(k):ta_end_je(k)) = &
                                      raw_sgsl_block(1:ncolumns,1:nrows) * scale_factor


             !           Print*, h_block
            
             DEALLOCATE (raw_sgsl_block, STAT=errorcode)
             IF(errorcode/=0) CALL abort_extpar('Cant deallocate the array raw_sgsl_block')
          
           ENDIF
         ENDDO
      



       END SUBROUTINE get_sgsl_data_block

!----------------------------------------------------------------------------------------------------------------


END MODULE mo_sgsl_routines

