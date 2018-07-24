!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     This module creates grid description files for input to the SCRIP code
!
!-----------------------------------------------------------------------

MODULE scripgrid_mod

  USE kinds_mod
  USE constants
  USE iounits
  USE netcdf
  USE netcdf_mod
  
  IMPLICIT NONE

  !-----------------------------------------------------------------------
  !     module variables that describe the grid

  INTEGER (kind=int_kind), parameter :: &
       grid_rank = 2,         &
       grid_corners = 4
  INTEGER (kind=int_kind) :: nx, ny, grid_size
  INTEGER (kind=int_kind), dimension(2) :: &
       grid_dims,    &  ! size of x, y dimensions
       grid_dim_ids     ! ids of the x, y dimensions
  INTEGER (kind=int_kind), ALLOCATABLE, DIMENSION(:) :: &
       grid_imask          ! land-sea mask
  REAL (kind=int_kind), ALLOCATABLE, DIMENSION(:) :: &
       grid_center_lat, &  ! lat/lon coordinates for
       grid_center_lon     ! each grid center in degrees
  REAL (kind=dbl_kind), ALLOCATABLE, DIMENSION(:,:) :: &
       grid_corner_lat, &  ! lat/lon coordinates for
       grid_corner_lon     ! each grid corner in degrees
  REAL (kind=dbl_kind), ALLOCATABLE, DIMENSION(:,:,:) :: &
       corner_lon, &
       corner_lat
  REAL (kind=dbl_kind), PARAMETER :: circle = 360.0

  !-----------------------------------------------------------------------
  !     module variables that describe the netcdf file

  INTEGER (kind=int_kind) :: &
       ncstat,           &   ! general netCDF status variable
       ncid_in

CONTAINS

  ! ==============================================================================
  
  SUBROUTINE convert(nm_in)
  
    ! -----------------------------------------------------------------------------
    ! - input variables
  
    CHARACTER(char_len), INTENT(in) ::  &
      nm_in
  
    ! -----------------------------------------------------------------------------
    ! - local variables
  
    CHARACTER(char_len) ::  &
      nemo_file, input_file, method, input_lon, input_lat, datagrid_file, &
      nemogrid_file, nemo_lon, nemo_lat, corn_lon, corn_lat, nemo_mask, input_mask
    INTEGER (kind=int_kind), dimension(2) :: &
      offset
    INTEGER (kind=int_kind) :: &
      iunit, nemo_mask_value, input_mask_value
  
    namelist /grid_inputs/ nemo_file, input_file, datagrid_file, nemogrid_file,  &
                           method, input_lon, input_lat, nemo_lon, nemo_lat,     &
                           nemo_mask, nemo_mask_value, input_mask, input_mask_value
  
    !-----------------------------------------------------------------------
    ! - namelist describing the processing
    !   note that mask_value is the minimum good value,
    !   so that where the mask is less than the value is masked

    nemo_file = "coordinates.nc"
    nemo_lon = "glamt"
    nemo_lat = "gphit"
    input_lon = "lon"
    input_lat = "lat"
    input_mask = "none"
    input_mask_value = 0
    datagrid_file = 'remap_data_grid.nc'
    nemogrid_file = 'remap_nemo_grid.nc'

    call get_unit(iunit)
    open(iunit, file=nm_in, status='old', form='formatted')
    read(iunit, nml=grid_inputs) 
    call release_unit(iunit)

    if (nemo_lon(1:4) .ne. 'glam' .or. nemo_lat(1:4) .ne. 'gphi') then
      write(6,*) 'lon name does not start with "glam" or lat name does not start with "gphi"'
      stop
    endif

    ! set up the names of the corner variables for a given input
    ! the offset represents what needs to be added to (i,j) to get to the correct 
    ! element in the corner arrays to correspond to the point northeast of the center
    if (nemo_lon(5:5) == "t") then
      corn_lon = "glamf"
      corn_lat = "gphif"
      offset = (/ 0,0 /)
    else if (nemo_lon(5:5) == "u") then
      corn_lon = "glamv"
      corn_lat = "gphiv"
      offset = (/ 1,0 /)
    else if (nemo_lon(5:5) == "v") then
      corn_lon = "glamu"
      corn_lat = "gphiu"
      offset = (/ 0,1 /)
    else
      write(6,*) 'unknown nemo_lon name'
      stop
    endif

    write(6,*) "processing " // trim(nemo_file)
    call convertNEMO(nemo_file, nemo_lon, nemo_lat, corn_lon, corn_lat, &
                     offset, nemogrid_file)

    write(6,*) "processing regular grid"
    call convertFLUX(input_file, input_lon, input_lat, &
                     input_mask, input_mask_value, datagrid_file)

  END SUBROUTINE convert
  
  ! ==============================================================================
  
  SUBROUTINE convertNEMO(grid_file_in, cent_lon, cent_lat, corn_lon, corn_lat, &
                         off, grid_file_out)
  
    !-----------------------------------------------------------------------
    !
    !     This routine converts a NEMO coordinates.nc file to a remapping grid file.
    !
    
    CHARACTER(char_len), INTENT(in) :: cent_lon, cent_lat, corn_lon, corn_lat
    INTEGER (kind=int_kind), INTENT(in), DIMENSION(2) :: off
    CHARACTER(char_len), INTENT(in) :: grid_file_out
    CHARACTER(char_len), INTENT(in) :: grid_file_in

    !-----------------------------------------------------------------------
    !     module variables that describe the grid
  
    CHARACTER(char_len), parameter ::  &
         grid_name = 'Remapped NEMO grid for SCRIP'
  
    !-----------------------------------------------------------------------
    !     grid coordinates and masks
  
    REAL (kind=dbl_kind), ALLOCATABLE, DIMENSION(:,:) :: &
         clon, clat, &             ! expanded corner arrays
         glam, &                   ! center longitude
         gphi, &                   ! center latitude
         glamc, &                  ! corner longitude
         gphic                     ! corner latitude
  
    !-----------------------------------------------------------------------
    !     other local variables
  
    INTEGER (kind=int_kind) :: i, j, n, iunit, im1, jm1, imid, isame, ic, jc
    INTEGER (kind=int_kind) :: varid_lam, varid_phi, varid_lamc, varid_phic
    INTEGER (kind=int_kind) :: jdim
    INTEGER (kind=int_kind), dimension(4) :: grid_dimids  ! input fields have 4 dims
    REAL (kind=dbl_kind) :: tmplon, dxt, dyt
  
    !-----------------------------------------------------------------------
    !     read in grid info
    !
    !     For NEMO input grids, assume that variable names are glam, glamc etc.
    !     Assume that 1st 2 dimensions of these variables are x and y directions.
    !     These assumptions are made by NEMO, so should be valid for coordinates.nc.
    !
    ! write in nf90 calls (without error handling) and then think about 
    ! making more readable by taking chunks into ncutil
    !
  
    ncstat = nf90_open( grid_file_in, NF90_NOWRITE, ncid_in )
    call netcdf_error_handler(ncstat)
  
    ! find dimids for 'glam'
    ! use dimids to get dimlengths
    ! allocate glam array
    ! get glam from file 
  
    ncstat = nf90_inq_varid( ncid_in, cent_lon, varid_lam )
    call netcdf_error_handler(ncstat)
    ncstat = nf90_inq_varid( ncid_in, corn_lon, varid_lamc )
    call netcdf_error_handler(ncstat)
    ncstat = nf90_inq_varid( ncid_in, cent_lat, varid_phi )
    call netcdf_error_handler(ncstat)
    ncstat = nf90_inq_varid( ncid_in, corn_lat, varid_phic )
    call netcdf_error_handler(ncstat)
  
    ncstat = nf90_inquire_variable( ncid_in, varid_lam, dimids=grid_dimids(:) )
    call netcdf_error_handler(ncstat)
    DO jdim = 1, SIZE(grid_dims)
       ncstat = nf90_inquire_dimension( ncid_in, grid_dimids(jdim),  &
                                        len=grid_dims(jdim) )
       call netcdf_error_handler(ncstat)
    END DO
    nx = grid_dims(1)
    ny = grid_dims(2)
    grid_size = nx * ny
    WRITE(*,FMT='("Input grid dimensions are:",2i6)') nx, ny
    
    ! assume that dimensions are all the same as glam
    ALLOCATE( glam(nx,ny), glamc(nx,ny), gphi(nx,ny), gphic(nx,ny) )
    ncstat = nf90_get_var( ncid_in, varid_lam, glam(:,:) )
    call netcdf_error_handler(ncstat)
    ncstat = nf90_get_var( ncid_in, varid_lamc, glamc(:,:) )
    call netcdf_error_handler(ncstat)
    ncstat = nf90_get_var( ncid_in, varid_phi, gphi(:,:) )
    call netcdf_error_handler(ncstat)
    ncstat = nf90_get_var( ncid_in, varid_phic, gphic(:,:) )
    call netcdf_error_handler(ncstat)
    
    !-----------------------------------------------------------------------
    ! - Mask is all ocean for now
  
    ALLOCATE( grid_imask(grid_size) )
    grid_imask(:) = 1
  
    !-----------------------------------------------------------------------
    ! corners are arranged as follows:    4 3
    !                                     1 2
    !
    ! Assume that cyclic grids have 2 wrap columns in coordinates.nc
    !  (this is the case for ORCA grids)
    !
  
    ! -----------------------------------------------------------------------------
    ! create a single pair of arrays for the corners where clon(1,1) corresponds
    ! to the south west corner of a box containing glam(1,1)
    ! various special cases then apply
    !   bottom row:  assume clon(:,j) = clon(:,j+1)
  
    ALLOCATE ( clon(nx+1,ny+1), clat(nx+1,ny+1) )

    ! first the easy internal points
    DO j = 2,ny
      DO i = 2,nx
        ic = i + off(1) - 1
        jc = j + off(2) - 1
        clon(i,j) = glamc(ic,jc)
        clat(i,j) = gphic(ic,jc)
      ENDDO
    ENDDO

    ! then the tricky boundary points
    imid = (nx-1)/2 + 1
    DO j = 1,ny+1,ny
      DO i = 1,nx+1,nx
        ic = i + off(1) - 1
        jc = j + off(2) - 1
        if (ic == 0 .and. jc == 0) then
          clon(i,j) = glamc(nx,1)
          clat(i,j) = gphic(nx,1) - (gphic(nx,2)-gphic(nx,1))
        else if (ic == nx+1 .and. jc == 0) then
          clon(i,j) = glamc(1,1)
          clat(i,j) = gphic(1,1) - (gphic(1,2)-gphic(1,1))
        else if (ic == 0 .and. jc == ny+1) then
          isame = 2*imid - nx + 1
          clon(i,j) = glamc(isame,jc-1)
          clat(i,j) = gphic(isame,jc-1)
        else if (ic == nx+1 .and. jc == ny+1) then
          isame = 2*imid
          clon(i,j) = glamc(isame,jc-1)
          clat(i,j) = gphic(isame,jc-1)
        else if (ic == 0) then
          clon(i,j) = glamc(nx,jc)
          clat(i,j) = gphic(nx,jc)
        else if (jc == 0) then
          clon(i,j) = glamc(ic,1)
          clat(i,j) = gphic(ic,1) - (gphic(ic,2)-gphic(ic,1))
        else if (ic == nx+1) then
          clon(i,j) = glamc(1,jc)
          clat(i,j) = gphic(1,jc)
        else if (jc == ny+1) then
          isame = 2*imid - ic + 1
          clon(i,j) = glamc(isame,jc-1)
          clat(i,j) = gphic(isame,jc-1)
        endif
      ENDDO
    ENDDO

    ALLOCATE ( corner_lon(4,nx,ny), corner_lat(4,nx,ny) )
  
    ! top-right corner
    corner_lon(3,:,:) = clon(2:nx+1,2:ny+1)
    corner_lat(3,:,:) = clat(2:nx+1,2:ny+1)
  
    ! top-left corner
    corner_lon(4,:,:) = clon(1:nx,2:ny+1)
    corner_lat(4,:,:) = clat(1:nx,2:ny+1)
  
    ! bottom-right corner 
    corner_lon(2,:,:) = clon(2:nx+1,1:ny)
    corner_lat(2,:,:) = clat(2:nx+1,1:ny)
  
    ! bottom-left corner
    corner_lon(1,:,:) = clon(1:nx,1:ny)
    corner_lat(1,:,:) = clat(1:nx,1:ny)
  
  ! For [N, E, W]-ward extrapolation near the poles, should we use stereographic (or
  ! similar) projection?  This issue will come for V,F interpolation, and for all
  ! grids with non-cyclic grids.
  
    ! -----------------------------------------------------------------------------
    !     correct for 0,2pi longitude crossings
    !      (In practice this means putting all corners into 0,2pi range
    !       and ensuring that no box corners are miles from each other.
    !       3pi/2 is used as threshold - I think this is quite arbitrary.)
  
    corner_lon(:,:,:) = MODULO( corner_lon(:,:,:), circle )
    DO n = 2, grid_corners
       WHERE    ( corner_lon(n,:,:) - corner_lon(n-1,:,:) < -three*circle*0.25 )
          corner_lon(n,:,:) = corner_lon(n,:,:) + circle
       ELSEWHERE( corner_lon(n,:,:) - corner_lon(n-1,:,:) >  three*circle*0.25 )
          corner_lon(n,:,:) = corner_lon(n,:,:) - circle
       END WHERE
    END DO
  
    ! -----------------------------------------------------------------------------
    ! - put longitudes on smooth grid 
  
 !  call mouldlon(glam,nx,ny)
 !  call mouldlon(corner_lon(1,:,:),nx,ny)
 !  call mouldlon(corner_lon(2,:,:),nx,ny)
 !  call mouldlon(corner_lon(3,:,:),nx,ny)
 !  call mouldlon(corner_lon(4,:,:),nx,ny)

    ! -----------------------------------------------------------------------------
    ! - reshape for SCRIP input format
  
    ALLOCATE( grid_center_lon(grid_size), grid_center_lat(grid_size) )
    
    grid_center_lon(:) = RESHAPE( glam(:,:), (/ grid_size /) )
    grid_center_lat(:) = RESHAPE( gphi(:,:), (/ grid_size /) )
  
    DEALLOCATE( glam, gphi, glamc, gphic )
  
    ALLOCATE( grid_corner_lon(4, grid_size), grid_corner_lat(4, grid_size) )
  
    grid_corner_lon(:,:) = RESHAPE( corner_lon(:,:,:), (/ 4, grid_size /) )
    grid_corner_lat(:,:) = RESHAPE( corner_lat(:,:,:), (/ 4, grid_size /) )
  
    DEALLOCATE( corner_lon, corner_lat )
  
    CALL createSCRIPgrid(grid_file_out, grid_name)
  
  END SUBROUTINE convertNEMO
  
  ! ==============================================================================
  
  SUBROUTINE convertFLUX(grid_file_in, name_lon, name_lat, &
                         name_mask, value_mask, grid_file_out)
  
    !-----------------------------------------------------------------------
    !
    !     This routine creates a remapping grid file from an input grid.
    !
    !-----------------------------------------------------------------------
    
    CHARACTER(char_len), INTENT(in) ::  &
         grid_file_in, name_lon, name_lat, name_mask, grid_file_out
    INTEGER (kind=int_kind) :: value_mask
  
    !-----------------------------------------------------------------------
    !     variables that describe the grid
  
    CHARACTER(char_len), parameter ::  &
         grid_name = 'Remapped regular grid for SCRIP'
  
    !-----------------------------------------------------------------------
    !     grid coordinates (note that a flux file just has lon and lat)
  
    REAL (kind=dbl_kind), ALLOCATABLE, DIMENSION(:) :: &
         lam, phi
    REAL (kind=dbl_kind), ALLOCATABLE, DIMENSION(:,:) :: &
         glam, &                   ! longitude
         gphi, &                   ! latitude
         glamc, &
         gphic
    REAL (kind=dbl_kind), ALLOCATABLE, DIMENSION(:,:) :: mask
  
    !-----------------------------------------------------------------------
    !     other local variables
  
    INTEGER (kind=int_kind) :: i, j, n, iunit, im1, jm1
    INTEGER (kind=int_kind) :: varid_lam, varid_phi, varid_mask
    INTEGER (kind=int_kind) :: jdim, nspace
    INTEGER (kind=int_kind), dimension(4) :: grid_dimids  ! input fields have 4 dims
    REAL (kind=dbl_kind) :: tmplon, dxt, dyt
  
    !-----------------------------------------------------------------------
    !     read in grid info
    !
    !     For NEMO input grids, assume that variable names are glam, glamc etc.
    !     Assume that 1st 2 dimensions of these variables are x and y directions.
    !     These assumptions are made by NEMO, so should be valid for coordinates.nc.
    !
    ! write in nf90 calls (without error handling) and then think about 
    ! making more readable by taking chunks into ncutil
  
    ncstat = nf90_open( grid_file_in, NF90_NOWRITE, ncid_in )
    call netcdf_error_handler(ncstat)
  
    ! find dimids for 'glamt'
    ! use dimids to get dimlengths
    ! allocate glam array
    ! get glam from file 
  
    ncstat = nf90_inq_varid( ncid_in, name_lat, varid_phi )
    call netcdf_error_handler(ncstat)
    ncstat = nf90_inq_varid( ncid_in, name_lon, varid_lam )
    call netcdf_error_handler(ncstat)
  
    ncstat = nf90_inquire_variable( ncid_in, varid_lam, ndims=nspace )
    call netcdf_error_handler(ncstat)

    if (nspace == 1) then
      ncstat = nf90_inquire_variable( ncid_in, varid_lam, dimids=grid_dimids(:1) )
      call netcdf_error_handler(ncstat)
      ncstat = nf90_inquire_variable( ncid_in, varid_phi, dimids=grid_dimids(2:) )
      call netcdf_error_handler(ncstat)
      ncstat = nf90_inquire_dimension( ncid_in, grid_dimids(1), len=grid_dims(1) )
      call netcdf_error_handler(ncstat)
      ncstat = nf90_inquire_dimension( ncid_in, grid_dimids(2), len=grid_dims(2) )
      call netcdf_error_handler(ncstat)
      nx = grid_dims(1)
      ny = grid_dims(2)
      grid_size = nx * ny
      WRITE(*,FMT='("Input grid dimensions are:",2i6)') nx, ny
    
      ALLOCATE( lam(nx), phi(ny) )
      write(6,*) 'double'
      ncstat = nf90_get_var( ncid_in, varid_lam, lam )
      call netcdf_error_handler(ncstat)
      ncstat = nf90_get_var( ncid_in, varid_phi, phi )
      call netcdf_error_handler(ncstat)
    
      ALLOCATE( glam(nx,ny), gphi(nx,ny))
      write(6,*) shape(lam),shape(phi)
      glam(:,:) = SPREAD(lam,2,ny)
      gphi(:,:) = SPREAD(phi,1,nx)
    else

      ncstat = nf90_inquire_variable( ncid_in, varid_lam, dimids=grid_dimids(:2) )
      call netcdf_error_handler(ncstat)
      ncstat = nf90_inquire_dimension( ncid_in, grid_dimids(1), len=grid_dims(1) )
      call netcdf_error_handler(ncstat)
      ncstat = nf90_inquire_dimension( ncid_in, grid_dimids(2), len=grid_dims(2) )
      call netcdf_error_handler(ncstat)
      nx = grid_dims(1)
      ny = grid_dims(2)
      grid_size = nx * ny
      WRITE(*,FMT='("Input grid dimensions are:",2i6)') nx, ny
    
      ALLOCATE( glam(nx,ny), gphi(nx,ny))
      ncstat = nf90_get_var( ncid_in, varid_lam, glam )
      call netcdf_error_handler(ncstat)
      ncstat = nf90_get_var( ncid_in, varid_phi, gphi )
      call netcdf_error_handler(ncstat)

    endif
    write(6,*) grid_size,nx,ny
  
    ALLOCATE(glamc(0:nx,0:ny), gphic(0:nx,0:ny) )

    ! - for now a simple average to get top right box corners
    ! - glamc(i,j), gphic(i,j) are top right coordinates of box containing 
    ! - glam(i,j),gphi(i,j)
    write(6,*) 'averaging'
    write(6,*) size(gphic),size(gphi)
    gphic(1:nx,1:ny-1) = 0.5*(gphi(:,1:ny-1)+gphi(:,2:ny))
    write(6,*) size(glamc),size(glam)
    glamc(1:nx-1,1:ny) = 0.5*(glam(1:nx-1,:)+glam(2:nx,:))
  
    ! - left and right column of longitudes 
    write(6,*) 'columns'
    glamc(nx,1:ny) = 1.5*glam(nx,:)-0.5*glam(nx-1,:)
    glamc( 0,1:ny) = 1.5*glam(1,:)-0.5*glam(2,:)
    glamc(nx, 0) = glamc(nx,1)
    glamc( 0, 0) = glamc( 0,1)

    ! - top and bottom row of latitudes by extrapolation
    write(6,*) 'rows'
    gphic(1:nx,ny) = 1.5*gphi(:,ny)-0.5*gphi(:,ny-1)
    gphic(1:nx, 0) = 1.5*gphi(:,1)-0.5*gphi(:,2)
    gphic( 0,ny) = gphic(1,ny)
    gphic( 0, 0) = gphic(1, 0)

    !-----------------------------------------------------------------------
  
    write(6,*) 'allocating'
    ALLOCATE( grid_imask(grid_size) )
    grid_imask(:) = 1
    write(6,*) name_mask
    if (trim(name_mask) /= "none") then
      write(6,*) 'masking'
      ncstat = nf90_inq_varid( ncid_in, name_mask, varid_mask )
      call netcdf_error_handler(ncstat)
      ALLOCATE( mask(nx,ny) )
      write(6,*) 'reading mask'
      ncstat = nf90_get_var( ncid_in, varid_mask, mask )
      call netcdf_error_handler(ncstat)
      write(6,*) 'setting mask'
      WHERE ( RESHAPE(mask(:,:),(/ grid_size /)) < value_mask)
        grid_imask = 0
      END WHERE
      write(6,*) 'masked'
    END IF
  
    !-----------------------------------------------------------------------
    ! corners are arranged as follows:    4 3
    !                                     1 2
  
    ALLOCATE ( corner_lon(4,nx,ny), corner_lat(4,nx,ny) )
  
    ! - bottom-left corner
    corner_lon(1,:,:) = glamc(0:nx-1, 0:ny-1 )
    corner_lat(1,:,:) = gphic(0:nx-1, 0:ny-1 )
  
    ! - bottom-right corner 
    corner_lon(2,:,:) = glamc(1:nx, 0:ny-1 )
    corner_lat(2,:,:) = gphic(1:nx, 0:ny-1 )
  
    ! - top-right corner
    corner_lon(3,:,:) = glamc(1:nx,1:ny)
    corner_lat(3,:,:) = gphic(1:nx,1:ny)
    write(6,*) corner_lat(3,nx-2:nx,ny)
  
    ! - top-left corner
    corner_lon(4,:,:) = glamc(0:nx-1, 1:ny )
    corner_lat(4,:,:) = gphic(0:nx-1, 1:ny )
  
  ! For [N, E, W]-ward extrapolation near the poles, should we use stereographic (or
  ! similar) projection?  This issue will come for V,F interpolation, and for all
  ! grids with non-cyclic grids.
  
    ! -----------------------------------------------------------------------------
    !     correct for 0,2pi longitude crossings
    !      (In practice this means putting all corners into 0,2pi range
    !       and ensuring that no box corners are miles from each other.
    !       3pi/2 is used as threshold - I think this is quite arbitrary.)
  
  ! corner_lon(:,:,:) = MODULO( corner_lon(:,:,:), circle )
  ! DO n = 2, grid_corners
  !    WHERE    ( corner_lon(n,:,:) - corner_lon(n-1,:,:) < -three*circle*0.25 )
  !       corner_lon(n,:,:) = corner_lon(n,:,:) + circle
  !    ELSEWHERE( corner_lon(n,:,:) - corner_lon(n-1,:,:) >  three*circle*0.25 )
  !       corner_lon(n,:,:) = corner_lon(n,:,:) - circle
  !    END WHERE
  ! END DO
  
    ! -----------------------------------------------------------------------------
    ! - reshape for SCRIP input format
  
    ALLOCATE( grid_center_lon(grid_size), grid_center_lat(grid_size) )
    
    grid_center_lon(:) = RESHAPE( glam(:,:), (/ grid_size /) )
    grid_center_lat(:) = RESHAPE( gphi(:,:), (/ grid_size /) )
  
    DEALLOCATE( glam, gphi, glamc, gphic )
  
    ALLOCATE( grid_corner_lon(4, grid_size), grid_corner_lat(4, grid_size) )
  
    grid_corner_lon(:,:) = RESHAPE( corner_lon(:,:,:), (/ 4, grid_size /) )
    grid_corner_lat(:,:) = RESHAPE( corner_lat(:,:,:), (/ 4, grid_size /) )
  
    DEALLOCATE( corner_lon, corner_lat )
  
    CALL createSCRIPgrid(grid_file_out, grid_name)
  
  END SUBROUTINE convertFLUX
  
  ! ==============================================================================
  
  SUBROUTINE mouldlon(lon_grid, nx, ny)

    ! -----------------------------------------------------------------------------
    ! - input variables

    INTEGER, INTENT(in) :: nx, ny
    REAL (kind=dbl_kind), INTENT(inout), DIMENSION(nx,ny) ::  &
      lon_grid
  
    ! -----------------------------------------------------------------------------
    ! - local variables

    INTEGER :: ix, iy
    REAL (kind=dbl_kind), DIMENSION(:,:), ALLOCATABLE ::  &
      dlon
    REAL :: step

    ! -----------------------------------------------------------------------------
    ! - try to eliminate any 360 degree steps in a grid of longitudes

    ALLOCATE(dlon(nx,ny))

    step = 0.75*circle
    dlon(:,:) = 0
    dlon(2:,:) = lon_grid(2:,:) - lon_grid(:nx-1,:)
    WHERE (dlon > -step .AND. dlon < step)
      dlon = 0.0
    ELSEWHERE
      dlon = -SIGN(circle,dlon)
    END WHERE

    ! - close your eyes this is nasty
    DO ix = 2,nx
      dlon(ix,:) = dlon(ix,:) + dlon(ix-1,:)
    END DO
    lon_grid = lon_grid + dlon

  END SUBROUTINE mouldlon
  
  ! ==============================================================================
  
  SUBROUTINE createSCRIPgrid(grid_file_out, grid_name)
  
    ! -----------------------------------------------------------------------------
    ! - input variables
  
    CHARACTER(char_len), INTENT(in) ::  &
         grid_name, grid_file_out
  
    ! -----------------------------------------------------------------------------
    ! - local variables that describe the netcdf file
  
    INTEGER (kind=int_kind) :: &
         nc_grid_id,       &   ! netCDF grid dataset id
         nc_gridsize_id,   &   ! netCDF grid size dim id
         nc_gridcorn_id,   &   ! netCDF grid corner dim id
         nc_gridrank_id,   &   ! netCDF grid rank dim id
         nc_griddims_id,   &   ! netCDF grid dimensions id
         nc_grdcntrlat_id, &   ! netCDF grid center lat id
         nc_grdcntrlon_id, &   ! netCDF grid center lon id
         nc_grdimask_id,   &   ! netCDF grid mask id
         nc_gridarea_id,   &   ! netCDF grid area id
         nc_grdcrnrlat_id, &   ! netCDF grid corner lat id
         nc_grdcrnrlon_id      ! netCDF grid corner lon id
  
    ! -----------------------------------------------------------------------------
    ! - create netCDF dataset for this grid
    ! - rewrite in nf90 
    ! - (bring out functional blocks into ncclear for readability)
  
    ncstat = nf90_create (grid_file_out, or(NF90_CLOBBER,NF90_NETCDF4), nc_grid_id)
    call netcdf_error_handler(ncstat)
    ncstat = nf90_put_att (nc_grid_id, NF90_GLOBAL, 'title', grid_name)
    call netcdf_error_handler(ncstat)
  
    ! -----------------------------------------------------------------------------
    ! - define grid size dimension
  
    ncstat = nf90_def_dim (nc_grid_id, 'grid_size', grid_size, nc_gridsize_id)
    call netcdf_error_handler(ncstat)
  
    ! -----------------------------------------------------------------------------
    ! - define grid rank dimension
  
    ncstat = nf90_def_dim (nc_grid_id, 'grid_rank', grid_rank, nc_gridrank_id)
    call netcdf_error_handler(ncstat)
  
    ! -----------------------------------------------------------------------------
    ! - define grid corner dimension
  
    ncstat = nf90_def_dim (nc_grid_id, 'grid_corners', grid_corners, nc_gridcorn_id)
    call netcdf_error_handler(ncstat)
  
    ! -----------------------------------------------------------------------------
    ! - define grid dim size array
  
    ncstat = nf90_def_var(nc_grid_id, 'grid_dims', NF90_INT, nc_gridrank_id, nc_griddims_id)
    call netcdf_error_handler(ncstat)
  
    ! -----------------------------------------------------------------------------
    ! - define grid mask
  
    ncstat = nf90_def_var(nc_grid_id, 'grid_imask', NF90_INT, &
                          nc_gridsize_id, nc_grdimask_id)
    call netcdf_error_handler(ncstat)
    ncstat = nf90_put_att(nc_grid_id, nc_grdimask_id, 'units', 'unitless')
    call netcdf_error_handler(ncstat)
  
    ! -----------------------------------------------------------------------------
    ! - define grid center latitude array
  
    ncstat = nf90_def_var(nc_grid_id, 'grid_center_lat', NF90_DOUBLE, &
                          nc_gridsize_id, nc_grdcntrlat_id)
    call netcdf_error_handler(ncstat)
    ncstat = nf90_put_att(nc_grid_id, nc_grdcntrlat_id, 'units', 'degrees')
    call netcdf_error_handler(ncstat)
  
    ! -----------------------------------------------------------------------------
    ! - define grid center longitude array
  
    ncstat = nf90_def_var(nc_grid_id, 'grid_center_lon', NF90_DOUBLE, &
                          nc_gridsize_id, nc_grdcntrlon_id)
    call netcdf_error_handler(ncstat)
    ncstat = nf90_put_att(nc_grid_id, nc_grdcntrlon_id, 'units', 'degrees')
    call netcdf_error_handler(ncstat)
  
    ! -----------------------------------------------------------------------------
    ! - define grid corner latitude array
  
    grid_dim_ids = (/ nc_gridcorn_id, nc_gridsize_id /)
    ncstat = nf90_def_var(nc_grid_id, 'grid_corner_lat', NF90_DOUBLE, &
                          grid_dim_ids, nc_grdcrnrlat_id)
    call netcdf_error_handler(ncstat)
    ncstat = nf90_put_att(nc_grid_id, nc_grdcrnrlat_id, 'units', 'degrees')
    call netcdf_error_handler(ncstat)
  
    ! -----------------------------------------------------------------------------
    ! - define grid corner longitude array
  
    ncstat = nf90_def_var(nc_grid_id, 'grid_corner_lon', NF90_DOUBLE, &
                          grid_dim_ids, nc_grdcrnrlon_id)
    call netcdf_error_handler(ncstat)
    ncstat = nf90_put_att(nc_grid_id, nc_grdcrnrlon_id, 'units', 'degrees')
    call netcdf_error_handler(ncstat)
  
    ! -----------------------------------------------------------------------------
    ! end definition stage
  
    ncstat = nf90_enddef(nc_grid_id)
    call netcdf_error_handler(ncstat)
  
    ! -----------------------------------------------------------------------------
    ! write grid data
  
    ncstat = nf90_put_var(nc_grid_id, nc_griddims_id, grid_dims)
    call netcdf_error_handler(ncstat)
    ncstat = nf90_put_var(nc_grid_id, nc_grdimask_id, grid_imask)
    call netcdf_error_handler(ncstat)
    ncstat = nf90_put_var(nc_grid_id, nc_grdcntrlat_id, grid_center_lat)
    call netcdf_error_handler(ncstat)
    ncstat = nf90_put_var(nc_grid_id, nc_grdcntrlon_id, grid_center_lon)
    call netcdf_error_handler(ncstat)
    ncstat = nf90_put_var(nc_grid_id, nc_grdcrnrlat_id, grid_corner_lat)
    call netcdf_error_handler(ncstat)
    ncstat = nf90_put_var(nc_grid_id, nc_grdcrnrlon_id, grid_corner_lon)
    call netcdf_error_handler(ncstat)
  
    ncstat = nf90_close(nc_grid_id)
    call netcdf_error_handler(ncstat)
  
    DEALLOCATE( grid_imask, grid_center_lon, grid_center_lat, &
                grid_corner_lon, grid_corner_lat )
  
  
  END SUBROUTINE createSCRIPgrid
  
END MODULE scripgrid_mod

