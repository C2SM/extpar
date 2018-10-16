!   Fortran module with data files for GLOBE and ASTER data
!
!  History:
!  Version                  Date                 Name
!  ------------------------ -------------------- -----------------
!  V4_0                     2016/07/28           Daniel Luethi
! 
!  Code Description:
!  Language: Fortran 90
!=================================================================
MODULE mo_sgsl_data

! Modules used: 

 USE mo_kind,               ONLY: wp,     &
                                  i4,     &
                                  i8,     &
                                  i2

 USE mo_grid_structures,    ONLY:  reg_lonlat_grid

 USE mo_utilities_extpar,   ONLY:  abort_extpar

 USE mo_io_utilities,       ONLY:  check_netcdf

 USE mo_io_units,           ONLY:  filename_max

 USE mo_sgsl_tg_fields,     ONLY:  sgsl

 USE netcdf,       ONLY :    &
     nf90_open,              &
     nf90_close,             &
     nf90_inquire,           &
     nf90_inquire_dimension, &
     nf90_inquire_variable,  &
     nf90_inq_attname,       &
     nf90_inquire_attribute, &
     nf90_get_att,           &
     nf90_inquire_dimension, &
     nf90_inq_dimid,         &
     nf90_inq_varid,         &
     nf90_get_var,           &
     nf90_noerr,             &
     nf90_strerror,          &
     nf90_create,            &
     nf90_def_dim,           &
     nf90_def_var,           &
     nf90_enddef,            &
     nf90_redef,             &
     nf90_put_att,           &
     nf90_put_var,           &
     NF90_CHAR,              &
     NF90_DOUBLE,            &
     NF90_FLOAT,             &
     NF90_INT,               &
     NF90_BYTE,              &
     NF90_SHORT,             &
     NF90_GLOBAL,            &
     NF90_UNLIMITED,         &
     NF90_CLOBBER,           &
     NF90_NOWRITE

IMPLICIT NONE

PRIVATE

 
PUBLIC :: ntiles,                  &   ! number of tiles in GLOBE / ASTER
          max_tiles,               &   ! maximal possible number of tiles that can be read
          nc_tot,                  &   ! total number of columns in GLOBE / ASTER data set
          nr_tot,                  &   ! total number of rows in GLOBE / ASTER data set 
          nc_tile ,                &   ! total number of columns in one GLOBE / ASTER tile
          tiles_lon_min,           &   ! starting longitude of every GLOBE / ASTER tile
          tiles_lon_max,           &   ! ending longitude of every GLOBE / ASTER tile
          tiles_lat_min,           &   ! starting latitude of every GLOBE / ASTER tile
          tiles_lat_max,           &   ! ending latitude of every GLOBE / ASTER tile
          tiles_ncolumns,          &   ! number of columns (lonitude increments) in each GLOBE / ASTER tile
          tiles_nrows,             &   ! number of rows (latitude increments) in each GLOBE /ASTER tile
          sgsl_tiles_grid,         &
          sgsl_grid,               &
          raw_sgsl_line,           &
          sgsl_tile_row,           &
          raw_sgsl_block,          &
          allocate_raw_sgsl_fields,&
          fill_sgsl_data,           &  ! subroutine (intent(in) and intent(out))
          num_tiles,                &  ! integer function
          allocate_sgsl_data,       &  ! subroutine (only intent(in))
          get_fill_value,           &
          get_varname,              &
          undef_sgsl,               &
          varname,                  &
          idem_type,                &
          dem_gl,                   &
          dem_aster,                &
          demraw_lat_min,           &
          demraw_lon_min,           &
          demraw_lat_max,           &
          demraw_lon_max,           &
          ntiles_row,               &
          ntiles_column,            &
          deallocate_sgsl_fields

SAVE
 
TYPE(reg_lonlat_grid), ALLOCATABLE  :: sgsl_tiles_grid(:)
TYPE(reg_lonlat_grid)               :: sgsl_grid  
      
INTEGER(KIND=i4), ALLOCATABLE :: tiles_ncolumns(:)
INTEGER(KIND=i4), ALLOCATABLE :: tiles_nrows(:)
INTEGER(KIND=i4), ALLOCATABLE :: sgsl_tile_row(:)

INTEGER(KIND=i4) :: ntiles
INTEGER(KIND=i4) :: nc_tot     
INTEGER(KIND=i4) :: nr_tot  
INTEGER(KIND=i4) :: nc_tile 
INTEGER(KIND=i4) :: idem_type
INTEGER(KIND=i4) :: ntiles_row
INTEGER(KIND=i4) :: ntiles_column

INTEGER, PARAMETER:: dem_gl = 1
INTEGER, PARAMETER:: dem_aster = 2
INTEGER, PARAMETER:: max_tiles = 1000

REAL(KIND=wp), ALLOCATABLE    :: tiles_lon_min(:)
REAL(KIND=wp), ALLOCATABLE    :: tiles_lon_max(:)
REAL(KIND=wp), ALLOCATABLE    :: tiles_lat_min(:)
REAL(KIND=wp), ALLOCATABLE    :: tiles_lat_max(:)
REAL(KIND=wp), ALLOCATABLE    :: raw_sgsl_line(:)
REAL(KIND=wp), ALLOCATABLE    :: raw_sgsl_block(:,:)


REAL(KIND=wp):: undef_sgsl
REAL(KIND=wp):: demraw_lat_min  
REAL(KIND=wp):: demraw_lat_max   
REAL(KIND=wp):: demraw_lon_min
REAL(KIND=wp):: demraw_lon_max

CHARACTER(LEN=80) :: varname

 
  CONTAINS


   SUBROUTINE num_tiles(columns,rows,ntiles) ! it gives the value of the number of tiles depending 
   IMPLICIT NONE
   SAVE
   INTEGER, INTENT(IN) :: columns
   INTEGER, INTENT(IN) :: rows
   INTEGER, INTENT(OUT):: ntiles           ! if the user chooses GLOBE or ASTER 

   ntiles_column = columns
   ntiles_row    = rows
   ntiles = ntiles_column * ntiles_row
   PRINT*, 'number of tiles is: ', ntiles

   END SUBROUTINE num_tiles



   SUBROUTINE allocate_sgsl_data(ntiles)   
! As it is unknown so far whether GLOBE or ASTER is chosen all parameters must be allocated in a second step.

   IMPLICIT NONE
   INTEGER, INTENT (IN) :: ntiles       ! number of tiles: 36 for ASTER and 16 for GLOBE
   INTEGER :: errorcode

   ALLOCATE (tiles_lon_min(1:ntiles), STAT = errorcode)
        IF (errorcode.NE.0) CALL abort_extpar('Cant allocate the vector tiles_lon_min')
   ALLOCATE (tiles_lon_max(1:ntiles), STAT = errorcode)
        IF (errorcode.NE.0) CALL abort_extpar('Cant allocate the vector tiles_lon_max')
   ALLOCATE (tiles_lat_min(1:ntiles), STAT = errorcode)
        IF (errorcode.NE.0) CALL abort_extpar('Cant allocate the vector tiles_lat_min')
   ALLOCATE (tiles_lat_max(1:ntiles), STAT = errorcode)
        IF (errorcode.NE.0) CALL abort_extpar('Cant allocate the vector tiles_lat_max')

   ALLOCATE (tiles_ncolumns(1:ntiles), STAT = errorcode)
        IF (errorcode.NE.0) CALL abort_extpar('Cant allocate the vector tiles_ncolumns')
   ALLOCATE (tiles_nrows(1:ntiles), STAT = errorcode)
        IF (errorcode.NE.0) CALL abort_extpar('Cant allocate the vector tiles_nrows')

  ALLOCATE (sgsl_tiles_grid(1:ntiles), STAT = errorcode)
        IF (errorcode.NE.0) CALL abort_extpar('Cant allocate the vector sgsl_tiles_grid')
   
   tiles_lon_min   = 0.0
   tiles_lon_max   = 0.0
   tiles_lat_min   = 0.0
   tiles_lat_max   = 0.0
   tiles_ncolumns  = 0
   tiles_nrows     = 0
   demraw_lat_min   = 0.0 
   demraw_lat_max   = 0.0 
   demraw_lon_min   = 0.0  
   demraw_lon_max   = 0.0  
   
   END SUBROUTINE allocate_sgsl_data




   SUBROUTINE fill_sgsl_data(raw_data_sgsl_path,sgsl_files,  &
                                                  tiles_lon_min,  &
                                                  tiles_lon_max,  &    
                                                  ! the allocated vectors need to be filled with the respective value.
                                                  tiles_lat_min,  &
                                                  tiles_lat_max,  &
                                                  nc_tot,         &
                                                  nr_tot,         &
                                                  nc_tile)
   IMPLICIT NONE
   SAVE
   CHARACTER (len=filename_max),INTENT(IN) :: raw_data_sgsl_path
   CHARACTER (len=filename_max),INTENT(IN) :: sgsl_files(1:ntiles)
   REAL(KIND=wp), INTENT(OUT)   :: tiles_lon_min(1:ntiles) 
   REAL(KIND=wp), INTENT(OUT)   :: tiles_lon_max(1:ntiles)
   REAL(KIND=wp), INTENT(OUT)   :: tiles_lat_min(1:ntiles)
   REAL(KIND=wp), INTENT(OUT)   :: tiles_lat_max(1:ntiles)
   INTEGER(KIND=i4), INTENT(OUT):: nc_tot, nr_tot, nc_tile
   CHARACTER(len=2)    :: num
   CHARACTER(len=80)   :: path
   INTEGER(KIND=i4)    :: i, errorcode     ! i is a counter, errorcode is used to check if allocation was successful
   INTEGER(KIND=i4)    :: ncid
   INTEGER(KIND=i4)    :: dimID_lat, dimID_lon, varID_lat, varID_lon                  
   REAL(KIND=wp)       :: half_gridp       ! distance of half a grid point as the grid point is centered on a GLOBE / ASTER pixel
  

    SELECT CASE (idem_type)                ! Also  could additionally be used for SELECT CASE (must first be read in)
     CASE(dem_aster)                       ! ASTER DEM, as it has 36 tiles at the moment.
       PRINT*, 'ASTER was used as DEM'
       half_gridp = 1./(3600.*2.)           ! the resolution of the ASTER data is 1./3600. degrees as it is half a grid point
                                            ! it is additionally divided by 2 
     CASE (dem_gl)                                           ! GLOBE DEM is composed of 16 tiles
       PRINT*, 'GLOBE was used as DEM'
       half_gridp = 1./(120.*2.)                              ! GLOBE resolution is 1./120. degrees (30 arc-seconds) 
     END SELECT

     DO i = 1,ntiles
     
       CALL check_netcdf(nf90_open(path =TRIM(raw_data_sgsl_path)//TRIM(sgsl_files(i)), mode = nf90_nowrite, ncid = ncid))    
       CALL check_netcdf(nf90_inq_dimid(ncid,"lon", dimID_lon))
       CALL check_netcdf(nf90_inq_dimid(ncid,"lat", dimID_lat))
       CALL check_netcdf(nf90_inquire_dimension(ncid,dimID_lon, len = tiles_ncolumns(i)))          
       CALL check_netcdf(nf90_inquire_dimension(ncid,dimID_lat, len = tiles_nrows(i))) 
       CALL check_netcdf(nf90_inq_varid(ncid, "lon", varID_lon))
       CALL check_netcdf(nf90_inq_varid(ncid, "lat", varID_lat))
       CALL check_netcdf(nf90_get_var(ncid, varID_lon, tiles_lon_min(i), start = (/1/)))            
       ! reads in the first longitude value of tile i
       CALL check_netcdf(nf90_get_var(ncid, varID_lon, tiles_lon_max(i), start = (/tiles_ncolumns(i)/))) 
       ! reads in the last longitude value of tile i
       CALL check_netcdf(nf90_get_var(ncid, varID_lat, tiles_lat_max(i), start = (/1/)))            
       ! reads in the first latitude value of tile i
       CALL check_netcdf(nf90_get_var(ncid, varID_lat, tiles_lat_min(i), start = (/tiles_nrows(i)/))) 
       ! reads in the last latitude value of tile i
       CALL check_netcdf(nf90_close(ncid))                                                         
       ! the netcdf file is closed again
       tiles_lon_min(i) = REAL(NINT(tiles_lon_min(i) - half_gridp)) !< half of a grid point must be
       tiles_lon_max(i) = REAL(NINT(tiles_lon_max(i) + half_gridp)) !< added, as the ASTER/GLOBE data
       tiles_lat_min(i) = REAL(NINT(tiles_lat_min(i) + half_gridp)) !< is located at the pixel center
       tiles_lat_max(i) = REAL(NINT(tiles_lat_max(i) - half_gridp))
     END DO


    SELECT CASE(idem_type)
     CASE(dem_gl)       
       WHERE (tiles_lon_max.GT.(-1*half_gridp).AND.tiles_lon_max.LT.half_gridp) tiles_lon_max = 0.00000   
       WHERE (tiles_lat_max.GT.(-1*half_gridp).AND.tiles_lat_max.LT.half_gridp) tiles_lat_max = 0.00000   
! There are probably some rounding problems, which are removed by this procedure.
     END SELECT


    SELECT CASE(idem_type)
     CASE(dem_aster, dem_gl) 
       demraw_lat_min = MINVAL(tiles_lat_min)
       demraw_lat_max = MAXVAL(tiles_lat_max)
       demraw_lon_min = MINVAL(tiles_lon_min)
       demraw_lon_max = MAXVAL(tiles_lon_max)
     END SELECT

     nc_tot = 0
     nr_tot = 0

     DO i = 1,ntiles
       nc_tot = nc_tot + tiles_ncolumns(i)                  
       nr_tot = nr_tot + tiles_nrows(i)
     END DO
     nc_tot = nc_tot/ntiles_row
     nr_tot = nr_tot/ntiles_column
     nc_tile = tiles_ncolumns(1)
    
  
   ALLOCATE (raw_sgsl_line(1:nc_tot), STAT = errorcode)
         IF (errorcode.NE.0) CALL abort_extpar('Cant allocate the vector raw_sgsl_line')
   ALLOCATE (sgsl_tile_row(1:nc_tile), STAT = errorcode)
         IF (errorcode.NE.0) CALL abort_extpar('Cant allocate the vector sgsl_tile_row')

   END SUBROUTINE fill_sgsl_data




   SUBROUTINE allocate_raw_sgsl_fields(nrows,ncolumns)
   IMPLICIT NONE
   INTEGER, INTENT(IN)  :: nrows
   INTEGER, INTENT(IN)  :: ncolumns

   INTEGER              :: errorcode

      ALLOCATE (raw_sgsl_block(1:ncolumns, 1:nrows), STAT = errorcode)
           IF (errorcode.NE.0) CALL abort_extpar('Cant allocate the array raw_sgsl_block')
              raw_sgsl_block = 0.0

   END SUBROUTINE allocate_raw_sgsl_fields



  SUBROUTINE get_fill_value(sgsl_file_1,undef_sgsl)
  IMPLICIT NONE
  SAVE
  CHARACTER (len=filename_max), INTENT(IN) :: sgsl_file_1     
  REAL(KIND=wp), INTENT(OUT)           :: undef_sgsl
  INTEGER(KIND=i2)           :: fillval
  REAL(KIND=wp)                  :: scale_factor
  INTEGER(KIND=i4)               :: ncid

  SELECT CASE(idem_type)
  
   CASE(dem_aster)
   CALL check_netcdf(nf90_open(path = sgsl_file_1, mode = nf90_nowrite, ncid = ncid))
   CALL check_netcdf(nf90_get_att(ncid, 3, "_FillValue", fillval))
   CALL check_netcdf(nf90_get_att(ncid, 3, "scale_factor", scale_factor))
   undef_sgsl = fillval * scale_factor
   CALL check_netcdf(nf90_close(ncid))

   CASE(dem_gl)
   CALL check_netcdf(nf90_open(path = sgsl_file_1 , mode = nf90_nowrite, ncid = ncid))
print *,'before get_att fillval'
   CALL check_netcdf(nf90_get_att(ncid, 3, "_FillValue", fillval))
print *,'before get_att scale_factor'
   CALL check_netcdf(nf90_get_att(ncid, 3, "scale_factor", scale_factor))
   undef_sgsl = fillval * scale_factor
   CALL check_netcdf(nf90_close(ncid))

  END SELECT

  END SUBROUTINE get_fill_value



  SUBROUTINE get_varname(sgsl_file_1,varname)
  IMPLICIT NONE
  SAVE
  CHARACTER (len=filename_max), INTENT(IN) :: sgsl_file_1     
  CHARACTER(LEN=*),INTENT(OUT)   :: varname
  INTEGER(KIND=i4)               :: ncid, type, ndims
  INTEGER(KIND=i4)               :: dimids(2)

  SELECT CASE(idem_type)
  
   CASE(dem_aster)
     CALL check_netcdf(nf90_open(path = sgsl_file_1, mode = nf90_nowrite, ncid = ncid))
     CALL check_netcdf(nf90_inquire_variable(ncid,3,varname,type,ndims,dimids))
     CALL check_netcdf(nf90_close(ncid))

   CASE(dem_gl)
   
     CALL check_netcdf(nf90_open(path = sgsl_file_1, mode = nf90_nowrite, ncid = ncid))
  
     CALL check_netcdf(nf90_inquire_variable(ncid,3,varname,type,ndims,dimids))
  
     CALL check_netcdf(nf90_close(ncid))
     varname = TRIM(varname)
  END SELECT

  END SUBROUTINE get_varname

   SUBROUTINE deallocate_sgsl_fields()

   IMPLICIT NONE

   INTEGER :: errorcode

   
   DEALLOCATE (sgsl_tiles_grid, STAT = errorcode)
   IF (errorcode.NE.0) CALL abort_extpar('Cant deallocate the vector sgsl_tiles_grid')
   DEALLOCATE (tiles_lon_min, STAT = errorcode)
   IF (errorcode.NE.0) CALL abort_extpar('Cant deallocate the vector tiles_lon_min')
   DEALLOCATE (tiles_lon_max, STAT = errorcode)
   IF (errorcode.NE.0) CALL abort_extpar('Cant deallocate the vector tiles_lon_max')
   DEALLOCATE (tiles_lat_min, STAT = errorcode)
   IF (errorcode.NE.0) CALL abort_extpar('Cant deallocate the vector tiles_lat_min')
   DEALLOCATE (tiles_lat_max, STAT = errorcode)
   IF (errorcode.NE.0) CALL abort_extpar('Cant deallocate the vector tiles_lat_max')
   DEALLOCATE (sgsl, STAT = errorcode)
   IF (errorcode.NE.0) CALL abort_extpar('Cant deallocate the array sgsl')
   
 END SUBROUTINE deallocate_sgsl_fields

END MODULE mo_sgsl_data
