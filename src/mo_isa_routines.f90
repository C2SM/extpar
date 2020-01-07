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

!> kind parameters are defined in MODULE data_parameters
USE mo_kind, ONLY: wp, &
                   i8, &
                   i4

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
USE mo_io_utilities,     ONLY: check_netcdf
USE mo_io_units,         ONLY: filename_max


USE mo_grid_structures,  ONLY: reg_lonlat_grid
USE mo_base_geometry,    ONLY: geographical_coordinates


IMPLICIT NONE

PRIVATE

PUBLIC :: read_namelists_extpar_isa


PUBLIC :: get_dimension_isa_data, &
  &       get_lonlat_isa_data,    &
  &       get_isa_tiles_grid,     &
  &       det_band_isa_data,      &
  &       get_isa_data_block,     &
  &       get_isa_tile_block_indices

CONTAINS

!---------------------------------------------------------------------------
!> subroutine to read namelist for orography data settings for EXTPAR 
SUBROUTINE read_namelists_extpar_isa(namelist_file,            &
    &                                  isa_type,    & !_br 14.04.16
                                       raw_data_isa_path,          &
                                         raw_data_isa_filename,      &
                                         ntiles_isa,          &
                                         isa_buffer_file      )


  USE mo_utilities_extpar, ONLY: free_un ! function to get free unit number
  
! >mes
  USE mo_isa_data,   ONLY: max_tiles_isa
!<mes
  
  CHARACTER (len=filename_max), INTENT(IN) :: namelist_file !< filename with namelists for for EXTPAR settings
! isa
  !INTEGER, INTENT(OUT) :: i_isa_data !< integer switch to choose a isa raw data set
                                         !! 1 Globcover2009, 2 GLC2000, 3 GLCC
  CHARACTER (len=filename_max), INTENT(OUT) :: raw_data_isa_path        !< path to raw data
  CHARACTER (len=filename_max), INTENT(OUT) :: raw_data_isa_filename(1:max_tiles_isa) !< filename isa raw data

  CHARACTER (len=filename_max), INTENT(OUT) :: isa_buffer_file !< name for isa buffer file
  INTEGER, INTENT(OUT) :: ntiles_isa
  INTEGER (KIND=i4)            :: isa_type  !< ID of dataset used !_br 14.04.16
!--
  !> namelist with isa data input
  !! HW: options from land-use left away
  ! NAMELIST /isa_raw_data/ raw_data_isa_path, raw_data_isa_filename, i_isa_data, ilookup_table_isa, ntiles_isa
  NAMELIST /isa_raw_data/ raw_data_isa_path, raw_data_isa_filename, ntiles_isa, isa_type !_br 14.04.16
  !> namelist with filenames for isa data output
  NAMELIST /isa_io_extpar/ isa_buffer_file

  INTEGER           :: nuin !< unit number
  INTEGER (KIND=i4) :: ierr !< error flag

  ntiles_isa=1
  nuin = free_un()  ! functioin free_un returns free Fortran unit number
  OPEN(nuin,FILE=TRIM(namelist_file), IOSTAT=ierr)

  READ(nuin, NML=isa_raw_data, IOSTAT=ierr)
  READ(nuin, NML=isa_io_extpar, IOSTAT=ierr)
   
  CLOSE(nuin)

!       print*,   namelist_file,        raw_data_isa_path,       &
!     &                                 raw_data_isa_filename,   &
!     &                                 ntiles_isa,       &
!     &                                 isa_buffer_file,         &
!     &                                 isa_output_file


END SUBROUTINE read_namelists_extpar_isa

               !> inquire dimension information for isa raw data
        SUBROUTINE get_dimension_isa_data(nlon_isa, &
                                          nlat_isa)

        USE mo_isa_data,   ONLY: max_tiles_isa,     &
                                       ntiles_isa, &
                                       len_isa_lon, len_isa_lat

        INTEGER (KIND=i8), INTENT(OUT) :: nlon_isa !< number of grid elements in zonal direction for isa data
        INTEGER (KIND=i8), INTENT(OUT) :: nlat_isa !< number of grid elements in meridional direction for isa data

        !local variables
        INTEGER, PARAMETER :: nx=129600
        INTEGER, PARAMETER :: ny=55800

        IF(ntiles_isa == 1) THEN
          nlon_isa = len_isa_lon
          nlat_isa = len_isa_lat
        ELSE
        nlon_isa = nx
        nlat_isa = ny
        END IF

       END SUBROUTINE get_dimension_isa_data

       !----------------------------------------------------------------------------------------------------------------
       !----------------------------------------------------------------------------------------------------------------

        !> get coordinates for isa raw data
        SUBROUTINE get_lonlat_isa_data(nlon_isa, &
           &                                 nlat_isa, &
           &                                 lon_isa,  &
           &                                 lat_isa,  &
           &                                 isa_grid)

       USE mo_grid_structures, ONLY: reg_lonlat_grid
       USE mo_isa_data,  ONLY: isa_tiles_lat_min, &
                                     isa_tiles_lat_max, &
                                     isa_tiles_lon_min, &
                                     isa_tiles_lon_max

        INTEGER (KIND=i8), INTENT(IN) :: nlon_isa !< number of grid elements in zonal direction for isa data
        INTEGER (KIND=i8), INTENT(IN) :: nlat_isa !< number of grid elements in meridional direction for isa data
        REAL (KIND=wp), INTENT(OUT)    :: lon_isa(1:nlon_isa) !< longitude of isa raw data
        REAL (KIND=wp), INTENT(OUT)    :: lat_isa(1:nlat_isa) !< latitude of isa raw data
        TYPE(reg_lonlat_grid), INTENT(OUT) :: isa_grid !< structure with defenition of the raw data grid for the whole GLOB &
!& ECOVER dataset

        !local variables
!        REAL, PARAMETER ::  xmin_glc = -180.001388888889 ! area of glcover data: western longitude
!        REAL, PARAMETER ::  xmax_glc  =  179.998611111111 ! area of glcover data: eastern longitude
!        REAL, PARAMETER ::  ymax_glc  =   90.001388888888! area of glcover data: northern latitude
!        REAL, PARAMETER ::  ymin_glc = -64.9986111111111! area of glcover data: southern latitude

!        REAL, PARAMETER :: dx_glc =    0.0027777777  ! grid element size of glcover data pixel in zonal direction
!        REAL, PARAMETER :: dy_glc =   -0.0027777777  ! grid element size of glcover data pixel in meridional directionon
! >mes
        REAL (KIND=wp) ::  xmin_glc ! area of glcover data: western longitude
        REAL (KIND=wp) ::  xmax_glc ! area of glcover data: eastern longitude
        REAL (KIND=wp) ::  ymax_glc ! area of glcover data: northern latitude
        REAL (KIND=wp) ::  ymin_glc ! area of glcover data: southern latitude

        REAL (KIND=wp):: dx_glc  ! grid element size of glcover data pixel in zonal direction
        REAL (KIND=wp):: dy_glc  ! grid element size of glcover data pixel in meridional directionon

        INTEGER (KIND=i8) :: jx,jy


        xmin_glc = MINVAL(isa_tiles_lon_min)
        xmax_glc = MAXVAL(isa_tiles_lon_max)
        ymax_glc = MAXVAL(isa_tiles_lat_max)
        ymin_glc = MINVAL(isa_tiles_lat_min)
        dx_glc   = (xmax_glc - xmin_glc)/ REAL(nlon_isa,wp)
        dy_glc   = -1.0 * (ymax_glc - ymin_glc) / REAL(nlat_isa,wp)

! <mes

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
!        isa_grid%start_lon_reg = xmin_glc
!        isa_grid%end_lon_reg   = xmax_glc
!        isa_grid%start_lat_reg = ymax_glc
!        isa_grid%end_lat_reg   = ymin_glc
        isa_grid%dlon_reg      = dx_glc ! (lon_isa(nlon_isa) - lon_isa(1)) / (nlon_isa - 1)
        isa_grid%dlat_reg      = dy_glc ! (lat_isa(nlat_isa) - lat_isa(1)) / (nlat_isa - 1)
        isa_grid%nlon_reg      = nlon_isa
        isa_grid%nlat_reg      = nlat_isa

       END SUBROUTINE get_lonlat_isa_data

        !----------------------------------------------------------------------------------------------------------------
       !----------------------------------------------------------------------------------------------------------------
! >mes
        SUBROUTINE get_isa_tiles_grid(isa_tiles_grid)
          USE mo_isa_data,  ONLY: ntiles_isa,   &
                                        isa_tiles_lon_min,   &
                                        isa_tiles_lon_max,   &
                                        isa_tiles_lat_min,   &
                                        isa_tiles_lat_max,   &
                                        isa_tiles_ncolumns,  &
                                        isa_tiles_nrows

          TYPE(reg_lonlat_grid), INTENT(OUT):: isa_tiles_grid(1:ntiles_isa)

          INTEGER::  k      !counter

          REAL(KIND=wp) :: dlon
          REAL(KIND=wp) :: dlat

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
! <mes

        !----------------------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------------------
! > mes
        !> determine grid description of band for ISA

        SUBROUTINE det_band_isa_data(isa_grid,start_isa_row,ta_grid)
          TYPE(reg_lonlat_grid),INTENT(IN) :: isa_grid ! sturcture with the definition of the global data grid of the GLOBC &
!& OVER data
          INTEGER, INTENT(IN) :: start_isa_row         ! number of the start row of band 
          TYPE(reg_lonlat_grid), INTENT(OUT):: ta_grid       ! structure with definition of the target area grid.

          INTEGER(KIND=i4):: nrows = 1000              ! number of rows, set to 1000 as default
          ! band from east to west for the whole globe, like the complete isa grid

          ta_grid%dlon_reg = isa_grid%dlon_reg
          ta_grid%dlat_reg = isa_grid%dlat_reg

          ta_grid%start_lon_reg = isa_grid%start_lon_reg
          ta_grid%end_lon_reg = isa_grid%end_lon_reg
          ta_grid%nlon_reg = isa_grid%nlon_reg

          !latitude from north to south, negative increment
          ta_grid%nlat_reg = nrows
          ta_grid%start_lat_reg = isa_grid%start_lat_reg + ta_grid%dlat_reg * (start_isa_row-1)   ! latitude from nor &
!& th to south, note the negative increment!
          ta_grid%end_lat_reg = ta_grid%start_lat_reg + ta_grid%dlat_reg * (nrows - 1)   ! latitude from north to south, note the &
!&  negative increment!
          
         ! check for the southern bound of the isa data
         IF (ta_grid%end_lat_reg < isa_grid%end_lat_reg) THEN ! band is at the southern bound
           ta_grid%end_lat_reg = isa_grid%end_lat_reg
           ta_grid%nlat_reg    = NINT(((ta_grid%end_lat_reg - ta_grid%start_lat_reg) / ta_grid%dlat_reg)) + 1
         END IF

       END SUBROUTINE det_band_isa_data
! < mes

        !----------------------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------------------
! >mes
! get ISA data block for a given target area from the tile block indices
       SUBROUTINE get_isa_data_block(ta_grid,              &
                                           isa_tiles_grid, &
                                           ncids_isa,      &
                                           isa_block)

         USE mo_grid_structures,  ONLY: reg_lonlat_grid  ! Definition of DATA Typeto describe a regular lonlat grid
         USE mo_isa_data,   ONLY: ntiles_isa
         TYPE(reg_lonlat_grid), INTENT(IN)  :: ta_grid !< structure with definition of the target area grid (dlon must be the sam &
!& e as for the whole ISA dataset)
         TYPE(reg_lonlat_grid), INTENT(IN) :: isa_tiles_grid(1:ntiles_isa) !< structure with defenition of the raw da &
!& ta grid for the 16 GLOBECOVER tiles
         INTEGER , INTENT(IN) :: ncids_isa(1:ntiles_isa)  !< ncid for the ISA tiles, the netcdf files have to b &
!& e opened previously
         REAL (KIND=wp), INTENT(OUT) :: isa_block(1:ta_grid%nlon_reg,1:ta_grid%nlat_reg) !< a block of ISA data 

       !local variables
       INTEGER (KIND=i4) :: isa_startrow(1:ntiles_isa) !< startrow indices for each ISA tile
       INTEGER (KIND=i4) :: isa_endrow(1:ntiles_isa) !< endrow indices for each ISA tile
       INTEGER (KIND=i4) :: isa_startcolumn(1:ntiles_isa) !< starcolumn indices for each ISA tile
       INTEGER (KIND=i4) :: isa_endcolumn(1:ntiles_isa) !< endcolumn indices for each ISA tile

       INTEGER (KIND=i4) :: ta_start_ie(1:ntiles_isa)    !< indices of target area block for first column of each ISA &
!&  tile
       INTEGER (KIND=i4) :: ta_end_ie(1:ntiles_isa)      !< indices of target area block for last column of each ISA  &
!& tile
       INTEGER (KIND=i4) :: ta_start_je(1:ntiles_isa)  !< indices of target area block for first row of each ISA tile &
!& 
       INTEGER (KIND=i4) :: ta_end_je(1:ntiles_isa)   !< indices of target area block for last row of each ISA tile


       INTEGER (KIND=i4), ALLOCATABLE :: raw_isa_block(:,:) !< a block with ISA data
       INTEGER :: varid               !< id of variable
       CHARACTER (LEN=80) :: varname  !< name of variable

       INTEGER :: nrows !< number of rows ! dimensions for raw_isa_block
       INTEGER :: ncolumns !< number of columns ! dimensions for raw_isa_block



       INTEGER :: k     ! counter
       INTEGER :: errorcode !< error status variable

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
 
!          print*,"(raw_isa_block: ncolumns ",ncolumns, "  nrows ",nrows
 
           ALLOCATE (raw_isa_block(1:ncolumns,1:nrows), STAT=errorcode)
             IF(errorcode/=0) CALL abort_extpar('Cant allocate the array raw_isa_block')

             CALL check_netcdf(nf90_inq_varid(ncids_isa(k),TRIM(varname),varid)) ! get the varid of the altitude variable
             ! get the data into the raw_isa_block
             CALL check_netcdf(nf90_get_var(ncids_isa(k), varid,  raw_isa_block,     & 
             &     start=(/isa_startcolumn(k),isa_startrow(k)/),count=(/ncolumns,nrows/)))
            

            isa_block(ta_start_ie(k):ta_end_ie(k),ta_start_je(k):ta_end_je(k)) = raw_isa_block(1:ncolumns,1:nrows)
            
             DEALLOCATE (raw_isa_block, STAT=errorcode)
              IF(errorcode/=0) CALL abort_extpar('Cant deallocate the array raw_isa_block')
          
           ENDIF
         ENDDO

       END SUBROUTINE get_isa_data_block
! < mes      
       !----------------------------------------------------------------------------------------------------------------
       !----------------------------------------------------------------------------------------------------------------
! >mes
       ! get startrow, endrow, startcolumn and endcolumn of each ISA tile (raw data) for a 
       ! given target area (ta_grid) and get start_indices (lon, lat) and end_indices of the target
       ! area for each ISA tile
       ! The ISA raw data are split in 6 tiles, so the target area may overlap several tiles.
       ! This subroutine determines the necessary indices to read in the ISA data into the
       ! target area.
       ! ISA tiles which are outside the target block will get indices with the value '0'

SUBROUTINE get_isa_tile_block_indices(ta_grid,              &
         &                                  isa_tiles_grid, &
         &                                  isa_startrow,   &
         &                                  isa_endrow,     & 
         &                                  isa_startcolumn,&
         &                                  isa_endcolumn,  &
         &                                  ta_start_ie,          &
         &                                  ta_end_ie,            &
         &                                  ta_start_je,          &
         &                                  ta_end_je)

USE mo_isa_data, ONLY : ntiles_isa,  &          !< ISA raw data has 6 tiles
                              isa_tiles_ncolumns, &
                              isa_tiles_nrows

       USE mo_grid_structures, ONLY: reg_lonlat_grid  !< Definition of data type to describe a regular (lonlat) grid
       TYPE(reg_lonlat_grid), INTENT(IN) :: ta_grid !< structure with definition of the target area grid (dlon must be the same a &
!& s for the whole ISA dataset)
       TYPE(reg_lonlat_grid), INTENT(IN) :: isa_tiles_grid(1:ntiles_isa) !< structure with defenition of the raw data &
!&  grid for the 6 ISA tiles

       INTEGER (KIND=i4), INTENT(OUT) :: isa_startrow(1:ntiles_isa)    !< startrow indices for each ISA tile
       INTEGER (KIND=i4), INTENT(OUT) :: isa_endrow(1:ntiles_isa)      !< endrow indices for each ISA tile

       INTEGER (KIND=i4), INTENT(OUT) :: isa_startcolumn(1:ntiles_isa)  !< starcolumn indices for each ISA tile &
!& 
       INTEGER (KIND=i4), INTENT(OUT) :: isa_endcolumn(1:ntiles_isa)   !< endcolumn indices for each ISA tile

       INTEGER (KIND=i4), INTENT(OUT) :: ta_start_ie(1:ntiles_isa)    !< indices of target area block for first column of e &
!& ach ISA tile
       INTEGER (KIND=i4), INTENT(OUT) :: ta_end_ie(1:ntiles_isa)      !< indices of target area block for last column of ea &
!& ch ISA tile
       INTEGER (KIND=i4), INTENT(OUT) :: ta_start_je(1:ntiles_isa)  !< indices of target area block for first row of each G &
!& LOBCOVER tile
       INTEGER (KIND=i4), INTENT(OUT) :: ta_end_je(1:ntiles_isa)   !< indices of target area block for last row of each GLO &
!& BCOVER tile

       ! local variables

       INTEGER  :: undefined

       INTEGER (KIND=i4) :: startrow ! startrow for tile
       INTEGER (KIND=i4) :: endrow 
       INTEGER (KIND=i4) :: startcolumn
       INTEGER (KIND=i4) :: endcolumn

       REAL (KIND=wp) :: dlon
       REAL (KIND=wp) :: dlat

       INTEGER :: k

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

          startcolumn = NINT((ta_grid%start_lon_reg - isa_tiles_grid(k)%start_lon_reg)/dlon) +1 ! here I want nearest index &
!&  (NINT)

          IF (startcolumn < 1) THEN 
            isa_startcolumn(k) = 1
            ! get the start index of the subtile for the target area block
            ta_start_ie(k) = NINT ((isa_tiles_grid(k)%start_lon_reg - ta_grid%start_lon_reg)/dlon) + 1 ! index of target ar &
!& ea block

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
           ta_end_ie(k) = NINT ((isa_tiles_grid(k)%end_lon_reg - ta_grid%start_lon_reg)/dlon) + 1                  ! index  &
!& of target area block
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
           ta_start_je(k) = NINT ((isa_tiles_grid(k)%start_lat_reg  - ta_grid%start_lat_reg)/dlat) + 1 ! index of target ar &
!& ea block
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
           ta_end_je(k) = NINT ((isa_tiles_grid(k)%end_lat_reg -  ta_grid%start_lat_reg )/dlat) + 1 ! index of target area  &
!& block

         ELSE IF (endrow < 1) THEN
           isa_endrow(k) = 0
           ta_end_je(k) = 0
         ELSE
           isa_endrow(k) = endrow
           ta_end_je(k) =  ta_grid%nlat_reg
         ENDIF
 
       ENDDO  ! loop over the tiles 

     END SUBROUTINE get_isa_tile_block_indices

! <mes

       !----------------------------------------------------------------------------------------------------------------
       !----------------------------------------------------------------------------------------------------------------
END MODULE mo_isa_routines




