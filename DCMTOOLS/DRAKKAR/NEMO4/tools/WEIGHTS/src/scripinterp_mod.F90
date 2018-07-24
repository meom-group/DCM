! **************************************************************************

  module scripinterp_mod

! ==========================================================================

  use kinds_mod    ! defines common data types
  use constants    ! defines common constants
  use iounits      ! I/O unit manager
  use netcdf
  use netcdf_mod   ! netcdf I/O stuff
  use grids        ! module containing grid info
  use remap_vars   ! module containing remapping info
  use remap_mod    ! module containing remapping routines
  use remap_read   ! routines for reading remap files

  implicit none

  real(kind=dbl_kind), dimension(:), allocatable ::   &
    grid_out

  integer(kind=int_kind) ::                           &  ! netCDF ids for files and arrays
          ncstat, nc_outfile_id, nc_infile_id
  integer (kind=int_kind), dimension(4) ::    &
      input_dims, input_dimids, input_count
  real (kind=dbl_kind), dimension(:), allocatable ::   &
      scale
  integer (kind=int_kind), dimension(:), allocatable ::        &
      nc_xysize_id, nc_gridsize_id, nc_gridsize,           &
      nc_variable_id
  integer :: nc_lon_id, nc_lat_id, nc_array_id

  character (char_len) ::  &
          map_name      ! name for mapping from grid1 to grid2
  character (1), parameter ::  &
          separator = '|'
  
  ! - input namelist variables
  
  character (char_len) ::  &
    input_file,            &  ! filename containing input fields
    interp_file,           &  ! filename containing remap data (map1)
    input_name                ! name of variable to grid
  integer (kind=int_kind), dimension(4) :: &
    input_stride,  &          ! how much of input array to process
    input_start,   &          ! where to start
    input_stop                ! and where to stop
  character (char_len), dimension(4) ::  &
    input_vars                ! input variables to be copied
  
  ! - output namelist variables
  
  character (char_len) ::  &
    output_file,           &  ! filename for test results
    output_mode,           &  ! 'create' or 'append'
    output_name,           &  ! name of new grid variable
    output_lat,            &  ! as it says
    output_lon,            &  ! as it says
    output_ydir               ! choose to invert output arrays in y direction
  character (char_len), dimension(4) ::  &
    output_dims,           &    ! name of new grid variable
    output_vars                 ! variables to copy
  character (char_len), dimension(10) ::  &
    output_attributes,     &  ! attributes of stuff in output file
    output_scaling             ! scaling factor to apply to input data to get
                               ! correct units in output file
contains

! ==========================================================================

  subroutine process_grid(nm_in)
  
    !-----------------------------------------------------------------------
    ! - dummy variables

    character (char_len) ::  &
      nm_in                    ! name of input namelist file
  
    !-----------------------------------------------------------------------
    ! - local variables

    integer (kind=int_kind), dimension(4) ::    &
        astart, acount, plus_one
    integer (kind=int_kind), dimension(3) ::    &
        write_dims
    integer (kind=int_kind) ::                  &
      i1, i2, jdim, n, nx, ny, nloop,           &
      nc_input_id, nc_input_rank,               &
      vstart, vstride, numv

    real (kind=dbl_kind), dimension(:), allocatable  ::    &
      grid1_array
    real (kind=dbl_kind), dimension(:,:), allocatable  ::    &
      var_out

    plus_one(:) = 1

    !-----------------------------------------------------------------------

    call read_mappings(nm_in)

    !-----------------------------------------------------------------------
    ! - read input grid
    ! - WARNING - lots of assumptions here at the moment
  
    ncstat = nf90_open( input_file, NF90_NOWRITE, nc_infile_id )
    call netcdf_error_handler(ncstat,"open")
  
    ncstat = nf90_inq_varid( nc_infile_id, input_name, nc_input_id )
    call netcdf_error_handler(ncstat,"inq_varid")
  
    input_dims(:) = 0
    ncstat = nf90_inquire_variable( nc_infile_id, nc_input_id,   &
                  ndims=nc_input_rank, dimids=input_dimids(:) )
    call netcdf_error_handler(ncstat,"inquire_variable")

    do jdim = 1,nc_input_rank
      ncstat = nf90_inquire_dimension(nc_infile_id,              &
                    input_dimids(jdim), len=input_dims(jdim) )
      call netcdf_error_handler(ncstat,"inquire_dimension")
    enddo
  
    ! - dimids seem to be returned in storage order so the outer dimension of 
    ! - the array as described by the netcdf file becomes the first dimension
    ! - as far as f90 is concerned (confused? you will be!)

    do jdim = 1,nc_input_rank
      if (input_stop(jdim) == 0) then
        input_stop(jdim) = input_dims(jdim)
      endif
      input_count(jdim) = input_stop(jdim) - input_start(jdim) + 1
    enddo

    ! - rashly we assume x followed by y
    nx = input_dims(1)
    ny = input_dims(2)
    write(*,fmt='("Input grid dimensions are:",2i6)') nx, ny
    if (nx*ny /= grid1_size) then
      write(6,*) "mismatch between input grid and remap data"
      stop
    endif

    ! - calculate number of horizontal slices to process
    ! - at the moment this is not very general and will only work with 3 dimensions

    acount(1:nc_input_rank) = &
             (input_stop(1:nc_input_rank)-input_start(1:nc_input_rank)+1) / &
             input_stride(1:nc_input_rank)
    nloop = 1
    do jdim = 1,nc_input_rank
      nloop = nloop*acount(jdim)
    enddo
    nloop = nloop/grid1_size
    write(6,*) "total slices requested: ",nloop

    vstart = input_start(nc_input_rank)            ! ie extra var has outer dimension
    vstride = input_stride(nc_input_rank)

    ! - in general we cant read in the whole array so do it slice by slice
    ! - slow but sure

    write(6,*) "allocating input and output grids"
    allocate( grid1_array(grid1_size))
    allocate( grid_out(grid2_size) )

    numv = 0
    do n = 1,4
      if (trim(input_vars(n)) /= '-' .and. &
          trim(output_vars(n)) /= '-') numv = numv + 1
    enddo

    write_dims(1) = grid2_dims(1)
    write_dims(2) = grid2_dims(2)
    write_dims(3) = nloop
    call define_grid(write_dims(1:3) , 2+numv)

    astart(:) = input_start(:)
    astart(3) = astart(3) - input_stride(3)
    acount(:) = 1
    acount(1) = nx
    acount(2) = ny

    do n = 1,nloop

      write(6,*) "processing slice: ",n
      astart(3) = astart(3) + input_stride(3)
      ncstat = nf90_get_var(nc_infile_id, nc_input_id, grid1_array,  &
                            start=astart(1:nc_input_rank),              &
                            count=acount(1:nc_input_rank))
      call netcdf_error_handler(ncstat,"get_var")

      call calculate_grid(grid1_array, grid_out)

      call write_grid(grid_out, n, write_dims(1:2) , 2)

    enddo

    ! ---------------------------------------------------------------------
    ! - now for any extra variables to copy

    if (numv > 0) then

      write(6,*) "reading ",numv," extra variables"
      allocate( var_out(nloop,numv) )
  
      do n = 1,numv
        write(6,*) "looking for variable: ",trim(input_vars(n))
        ncstat = nf90_inq_varid( nc_infile_id, trim(input_vars(n)), nc_input_id )
        call netcdf_error_handler(ncstat,"inq_varid")
      
        input_dims(:) = 0
        ncstat = nf90_inquire_variable( nc_infile_id, nc_input_id,   &
                      ndims=nc_input_rank, dimids=input_dimids(:) )
        call netcdf_error_handler(ncstat,"inquire_variable")
  
        if (nc_input_rank /= 1) then
          write(6,*) 'sorry, only rank 1 variables can be copied'
          cycle
        endif
        ncstat = nf90_inquire_dimension(nc_infile_id,              &
                      input_dimids(1), len=input_dims(1) )
        call netcdf_error_handler(ncstat,"inquire_dimension")
  
        ncstat = nf90_get_var(nc_infile_id, nc_input_id, var_out(1:nloop,n),  &
                              start=(/ vstart /), stride=(/ vstride /))
        call netcdf_error_handler(ncstat,"get_var")
      enddo
  
      call write_extra(var_out, numv+2)
      deallocate(var_out)

    endif

    ncstat = nf90_close(nc_outfile_id)
    call netcdf_error_handler(ncstat,"out close")
    ncstat = nf90_close(nc_infile_id)
    call netcdf_error_handler(ncstat,"in close")

    ! ---------------------------------------------------------------------

    deallocate( grid1_array, grid_out)

    ! ---------------------------------------------------------------------

  end subroutine process_grid

  ! ==========================================================================

  subroutine define_grid(thedims, therank)

    !-----------------------------------------------------------------------
    ! - dummy variables
    
    integer (kind=int_kind) :: &
      therank
    integer (kind=int_kind), dimension(therank) :: &
      thedims

    !-----------------------------------------------------------------------
    ! - local variables

    integer ::               &
      k, n, ilon, ilat, icolon, i1, i2, natt, nvar, id, jd, kd, nd
    character (char_len) ::  &
      aname, vname, att
    real (kind=dbl_kind) ::   s

    ! - netcdf variables

    integer :: xtype

    !-----------------------------------------------------------------------
    ! -  define grid size dimensions
  
    allocate(nc_xysize_id(grid2_rank))
    allocate(nc_gridsize_id(therank))
    allocate(nc_gridsize(therank))
    allocate(nc_variable_id(therank-2))
    
    !-----------------------------------------------------------------------
    ! - setup a NetCDF file for output
  
    xtype = NF90_FLOAT

    write(6,*) 'creating output file'
    ncstat = nf90_create (output_file, or(NF90_CLOBBER,NF90_NETCDF4), nc_outfile_id)
    call netcdf_error_handler(ncstat,"create")
  
    write(6,*) 'setting global attributes'
    ncstat = nf90_put_att(nc_outfile_id, NF90_GLOBAL, 'title', map_name)
    call netcdf_error_handler(ncstat,"put_att")
  
    write(6,*) 'setting dimensions'
    do n=1,therank
      if (n .eq. therank .and. therank .gt. 2) then
        write(6,*) '    unlimited dim ',trim(output_dims(n)),' size: ',thedims(n)
        ncstat = nf90_def_dim (nc_outfile_id, output_dims(n), NF90_UNLIMITED, &
                               nc_gridsize_id(n))
      else
        write(6,*) '    dim ',trim(output_dims(n)),' size: ',thedims(n)
        ncstat = nf90_def_dim (nc_outfile_id, output_dims(n), thedims(n), &
                               nc_gridsize_id(n))
      endif
      call netcdf_error_handler(ncstat,"def_dim")
    end do
    nc_gridsize(:) = thedims(1:therank)

    ! - at the moment there is an assumption here that the ordering is (lon,lat)

    ilon = 1
    ilat = 2
    nc_xysize_id(1) = nc_gridsize_id(ilon)
    nc_xysize_id(2) = nc_gridsize_id(ilat)

    ! ----------------------------------------------------------------
    ! -  define grid center longitude array
  
    write(6,*) 'defining longitude variable'
    ncstat = nf90_def_var (nc_outfile_id, output_lon,   &
                         xtype, nc_xysize_id,   &
                         nc_lon_id)
    call netcdf_error_handler(ncstat,"def_var")
  
    ncstat = nf90_put_att (nc_outfile_id, nc_lon_id, 'units', 'degrees')
    call netcdf_error_handler(ncstat,"put_att")
  
    ! ----------------------------------------------------------------
    ! -  define grid center latitude array
  
    write(6,*) 'defining latitude variable'
    ncstat = nf90_def_var (nc_outfile_id, output_lat,     &
                         xtype, nc_xysize_id,   &
                         nc_lat_id)
    call netcdf_error_handler(ncstat,"def_var")
  
    ncstat = nf90_put_att (nc_outfile_id, nc_lat_id, 'units', 'degrees')
    call netcdf_error_handler(ncstat,"put_att")
  
    ! ----------------------------------------------------------------
    ! -  define copy variables array
  
    write(6,*) 'defining copy variables'
    do n = 3,therank
      ncstat = nf90_def_var (nc_outfile_id, output_vars(n-2),   &
                           xtype, nc_gridsize_id(n),    &
                           nc_variable_id(n-2))
      call netcdf_error_handler(ncstat,"def_var")
    enddo
  
    ! ----------------------------------------------------------------
    ! -  define output array
  
    write(6,*) 'defining grid variable'
    ncstat = nf90_def_var (nc_outfile_id, output_name,               &
                         xtype, nc_gridsize_id,    &
                         nc_array_id)
    call netcdf_error_handler(ncstat,"def_var")
  
    ! ----------------------------------------------------------------
    ! - output attributes has to come after all other definitions
    ! - this code currently a bit murky, needs a rewrite
  
    ncstat = nf90_inquire (nc_outfile_id, nVariables=nvar)
    call netcdf_error_handler(ncstat,"inquire")
    do n = 1,10
      att = trim(output_attributes(n))
      natt = len(att)
      if (att /= '-') then
        i1 = index(att,separator)
        aname = att(1:i1-1)
        do k = 1,nvar
          ncstat = nf90_inquire_variable(nc_outfile_id, k, vname)
          call netcdf_error_handler(ncstat,"inquire_variable")
          if (vname == aname) then
            i2 = index(att,separator,.true.)
            ncstat = nf90_put_att (nc_outfile_id, k, &
                                   att(i1+1:i2-1), att(i2+1:natt))
            call netcdf_error_handler(ncstat,"put_att")
            exit      ! from innermost do
          endif
        enddo
      endif
    enddo
  
    ! output scaling

    allocate (scale(nvar))
    scale(:) = 1.0

    do n = 1,10
      att = trim(output_scaling(n))
      natt = len(att)
      if (att /= '-') then
        i1 = index(att,separator)
        aname = att(1:i1-1)
        do k = 1,nvar
          ncstat = nf90_inquire_variable(nc_outfile_id, k, vname)
          call netcdf_error_handler(ncstat,"inquire_variable")
          if (vname == aname) then
            i2 = index(att,separator,.true.)
            read(att(i2+1:natt),*) scale(k)
            call netcdf_error_handler(ncstat,"put_att")
            exit      ! from innermost do
          endif
        enddo
      endif
    enddo
  
    ! ----------------------------------------------------------------
    ! -  end definition stage
  
    ncstat = nf90_enddef(nc_outfile_id)
    call netcdf_error_handler(ncstat,"enddef")
 
  end subroutine define_grid

  ! ==========================================================================

  subroutine write_grid(thegrid, thelevel, thedims, therank)

    !-----------------------------------------------------------------------
    ! - dummy variables
    
    integer (kind=int_kind), intent(in) :: &
      therank, thelevel
    real (kind=dbl_kind), dimension(:), intent(in) ::   &
      thegrid
    integer (kind=int_kind), dimension(therank) :: &
      thedims

    !-----------------------------------------------------------------------
    ! - local variables

    integer ::               &
      k, n, ilon, ilat, icolon, j1, j2, dj, natt, nvar, id, jd, kd, nd
    character (char_len) ::  &
      aname, vname, att
    real (kind=dbl_kind), dimension(:,:,:), allocatable ::   &
      data
    real (kind=dbl_kind) ::   s
    real (kind=dbl_kind), parameter :: todeg = 57.295779513082323
    integer (kind=int_kind), dimension(3) ::   &
      start

    ! - netcdf variables

    integer :: xtype

    !-----------------------------------------------------------------------
    ! - write results to NetCDF file
  
    allocate (data(thedims(1),thedims(2),1))
    if (output_ydir .eq. 'invert') then
      j1 = thedims(2)
      j2 = 1
      dj = -1
    else
      j1 = 1
      j2 = thedims(2)
      dj = 1
    endif

    if (thelevel .eq. 1) then

      ! -  grid center latitude array
  
      write(6,*) 'writing latitude variable'
      s = scale(nc_lat_id)
      nd = 0
      do jd = j1,j2,dj
        do id =1,thedims(1)
          nd = nd + 1
          data(id,jd,1) = s*todeg*grid2_center_lat(nd)
        enddo
      enddo
      ncstat = nf90_put_var(nc_outfile_id, nc_lat_id, data(:,:,1))
      call netcdf_error_handler(ncstat,"put_var")
  
      ! -  grid center longitude array
  
      write(6,*) 'writing longitude variable'
      s = scale(nc_lon_id)
      nd = 0
      do jd = j1,j2,dj
        do id =1,thedims(1)
          nd = nd + 1
          data(id,jd,1) = s*todeg*grid2_center_lon(nd)
        enddo
      enddo
      ncstat = nf90_put_var(nc_outfile_id, nc_lon_id, data(:,:,1))
      call netcdf_error_handler(ncstat,"put_var")

    endif

    !-----------------------------------------------------------------------
    ! -  new grid

    write(6,*) 'writing grid variable'
    n = therank
    s = scale(nc_array_id)
    nd = 0
    do jd = j1,j2,dj
      do id =1,thedims(1)
        nd = nd + 1
        data(id,jd,1) = thegrid(nd)
      enddo
    enddo
    write(6,*) 'scaling data '
    data(:,:,1) = s*data(:,:,1)
    start(:) = (/ 1, 1, thelevel /)
    ncstat = nf90_put_var(nc_outfile_id, nc_array_id, data, start)
    call netcdf_error_handler(ncstat,"put_var")
    deallocate(data)

  end subroutine write_grid

  ! ==========================================================================

  subroutine write_extra(thevars, therank)

    real (kind=dbl_kind), dimension(:,:), intent(in) ::   &
      thevars
    real (kind=dbl_kind), dimension(:), allocatable ::   &
      thedata
    integer (kind=int_kind), intent(in) :: &
      therank
    real (kind=dbl_kind) ::   s
    integer :: n

    allocate( thedata(size(thevars,1)) )

    ! -  copy variable arrays

    write(6,*) 'writing copy variables'
    do n = 3,therank
      s = scale(nc_variable_id(n-2))
      thedata(:) = s*thevars(:,n-2)
      ncstat = nf90_put_var(nc_outfile_id, nc_variable_id(n-2), thedata)
      call netcdf_error_handler(ncstat,"put_var")
    enddo

    deallocate( thedata )

  end subroutine write_extra

  ! ==========================================================================

  subroutine close_grid()

  !     close netCDF file
  
    write(6,*) 'closing file'
    ncstat = nf90_close(nc_outfile_id)
    call netcdf_error_handler(ncstat,"close")
  
  end subroutine close_grid

  ! ==========================================================================

  subroutine read_mappings(nm_in)

    !-----------------------------------------------------------------------
    ! - dummy variables

    character (char_len) ::  &
      nm_in                    ! name of input namelist file
  
    !-----------------------------------------------------------------------
    ! - local variables
  
    character (char_len) ::   &
              dim_name    ! netCDF dimension name
  
    integer (kind=int_kind) ::     &
        iunit                                                ! unit number for namelist file
  
    !-----------------------------------------------------------------------
    ! - namelist block

    namelist /interp_inputs/ input_file, interp_file, input_name, &
                             input_stride, input_start, input_stop, &
                             input_vars
    namelist /interp_outputs/ output_dims, output_file, output_mode, output_name, &
                              output_lat, output_lon, output_ydir,  &
                              output_scaling, output_vars, output_attributes
  
    !-----------------------------------------------------------------------
    ! - read namelist for file and mapping info
  
    input_stride(:) = 1
    input_start(:) = 1
    input_stop(:) = 0
    output_scaling(:) = '-'
    input_vars(:) = '-'
    output_lon = '-'
    output_lat = '-'
    output_vars(:) = '-'
    output_ydir = 'none'
    output_attributes(:) = '-'

    call get_unit(iunit)
    open(iunit, file=nm_in, status='old', form='formatted')
    read(iunit, nml=interp_inputs)
    read(iunit, nml=interp_outputs)
    call release_unit(iunit)
    write(*,nml=interp_inputs)
    write(*,nml=interp_outputs)
    if (trim(output_mode) == "create") then
      if (trim(output_lon) == '-' .or. trim(output_lat) == '-') then
        write(6,*) 'if creating, need to supply lon and lat names'
        stop
      endif
    endif

    !-----------------------------------------------------------------------
    ! - read remapping data
    ! - via the scrip package this sets variables:
    !    grid1_size, grid2_size:     sizes of input and output grids
    !    grid1_mask, grid2_mask:     masks
    !    grid1_rank, grid2_rank:     ranks
  
    call read_remap(map_name, interp_file)
  
  end subroutine read_mappings

  ! ==========================================================================

  subroutine calculate_grid(grid1_array, grid2_array)

    !-----------------------------------------------------------------------
    ! - dummy variables

    real (kind=dbl_kind), intent(in), dimension(:) ::    &
      grid1_array
    real (kind=dbl_kind), intent(out), dimension(:) ::    &
      grid2_array

    !-----------------------------------------------------------------------
    ! - local variables

    integer (kind=int_kind), dimension(:), allocatable :: &
        grid1_imask, grid2_imask, grid2_count

    real (kind=dbl_kind), dimension(:), allocatable :: &
        grid1_tmp,           &
        grad1_lat,           &
        grad1_lon,           &
        grad1_latlon,        &
        grad1_lat_zero,      &
        grad1_lon_zero,      &
        grid2_tmp1,           &
        grid2_tmp2

    real (kind=dbl_kind) ::  &
        delew, delns            ! variables for computing bicub gradients

    integer (kind=int_kind) ::                          &
        i,j,n,imin,imax,idiff,                          &
        ip1,im1,jp1,jm1,nx,ny,                          &  ! for computing bicub gradients
        in,is,ie,iw,ine,inw,ise,isw

    logical, parameter :: lat_gradient = .false.
  
    write(6,*) 'starting'

    !-----------------------------------------------------------------------
    ! - allocate arrays
  
    allocate (grid1_tmp      (grid1_size),   &
              grad1_lat      (grid1_size),    &
              grad1_lon      (grid1_size),    &
              grad1_lat_zero (grid1_size),    &
              grad1_lon_zero (grid1_size),    &
              grid1_imask    (grid1_size),   &
              grid2_tmp1     (grid2_size),   &
              grid2_tmp2     (grid2_size),   &
              grid2_imask    (grid2_size),   &
              grid2_count    (grid2_size))
  
    write(6,*) 'allocated'
    write(6,*) grid1_size,grid2_size

    grid1_imask(:) = 1
    grid2_imask(:) = 1
    where (grid1_mask)
      grid1_imask = 1
    elsewhere
      grid1_imask = 0
    endwhere
    where (grid2_mask)
      grid2_imask = 1
    elsewhere
      grid2_imask = 0
    endwhere

    write(6,*) 'masked'

    grad1_lat_zero = zero
    grad1_lon_zero = zero
    nx = input_dims(1)
    ny = input_dims(2)
    write(6,*) nx,ny
  
    !-----------------------------------------------------------------------
    ! - if bicubic, we need 3 gradients in logical space
  
    if (map_type == map_type_bicubic) then
  
      write(6,*) 'bicubic'
      write(6,*) grid1_size

      allocate (grad1_latlon (grid1_size)) 
  
      do n=1,grid1_size
  
        grad1_lat(n) = zero
        grad1_lon(n) = zero
        grad1_latlon(n) = zero
  
!       if (n.ge.8000) write(6,*) 0,grid1_mask(n),nx
        if (grid1_mask(n)) then
  
          delew = half
          delns = half
  
          j = (n-1)/nx + 1
          i = n - (j-1)*nx
  
          ip1 = i+1
          im1 = i-1
          jp1 = j+1
          jm1 = j-1
  
          if (ip1 > nx) ip1 = ip1 - nx
          if (im1 < 1 ) im1 = nx
          if (jp1 > ny) then
            jp1 = j
            delns = one
          endif
          if (jm1 < 1 ) then
            jm1 = j
            delns = one
          endif
  
          in  = (jp1-1)*nx + i
          is  = (jm1-1)*nx + i
          ie  = (j  -1)*nx + ip1
          iw  = (j  -1)*nx + im1
  
          ine = (jp1-1)*nx + ip1
          inw = (jp1-1)*nx + im1
          ise = (jm1-1)*nx + ip1
          isw = (jm1-1)*nx + im1
  
          ! -  compute i-gradient
  
          if (.not. grid1_mask(ie)) then
            ie = n
            delew = one
          endif
          if (.not. grid1_mask(iw)) then
            iw = n
            delew = one
          endif
   
          grad1_lat(n) = delew*(grid1_array(ie) - grid1_array(iw))
!         if (n.ge.8000) write(6,*) 1,grad1_lat(n)
  
          ! -  compute j-gradient
  
          if (.not. grid1_mask(in)) then
            in = n
            delns = one
          endif
          if (.not. grid1_mask(is)) then
            is = n
            delns = one
          endif
   
          grad1_lon(n) = delns*(grid1_array(in) - grid1_array(is))
!         if (n.ge.8000) write(6,*) 2,grad1_lon(n)
  
          ! -  compute ij-gradient
  
          delew = half
          if (jp1 == j .or. jm1 == j) then
            delns = one
          else 
            delns = half
          endif
  
          if (.not. grid1_mask(ine)) then
            if (in /= n) then
              ine = in
              delew = one
            else if (ie /= n) then
              ine = ie
              inw = iw
              if (inw == n) delew = one
              delns = one
            else
              ine = n
              inw = iw
              delew = one
              delns = one
            endif
          endif
  
          if (.not. grid1_mask(inw)) then
            if (in /= n) then
              inw = in
              delew = one
            else if (iw /= n) then
              inw = iw
              ine = ie
              if (ie == n) delew = one
              delns = one
            else
              inw = n
              ine = ie
              delew = one
              delns = one
            endif
          endif
  
          grad1_lat_zero(n) = delew*(grid1_array(ine) - grid1_array(inw))
!         if (n.ge.8000) write(6,*) 3,grad1_lat_zero(n)
  
          if (.not. grid1_mask(ise)) then
            if (is /= n) then
              ise = is
              delew = one
            else if (ie /= n) then
              ise = ie
              isw = iw
              if (isw == n) delew = one
              delns = one
            else
              ise = n
              isw = iw
              delew = one
              delns = one
            endif
          endif
  
          if (.not. grid1_mask(isw)) then
            if (is /= n) then
              isw = is
              delew = one
            else if (iw /= n) then
              isw = iw
              ise = ie
              if (ie == n) delew = one
              delns = one
            else
              isw = n
              ise = ie
              delew = one
              delns = one
            endif
          endif
  
          grad1_lon_zero(n) = delew*(grid1_array(ise) - grid1_array(isw))
          grad1_latlon(n) = delns*(grad1_lat_zero(n) - grad1_lon_zero(n))
!         if (n.ge.8000) write(6,*) 4,grad1_lon_zero(n),grad1_latlon(n)
  
        endif
      enddo
  
      write(6,*) 'remapping'
      call remap(grid2_array, wts_map1, grid2_add_map1, grid1_add_map1, &
                 grid1_array, src_grad1=grad1_lat,                    &
                 src_grad2=grad1_lon, src_grad3=grad1_latlon)
  
      print *,'Third order mapping from grid1 to grid2:'
      print *,'----------------------------------------'
      print *,'Grid1 min,max: ',minval(grid1_array),maxval(grid1_array)
      print *,'Grid2 min,max: ',minval(grid2_array  ),maxval(grid2_array  )
    
      ! -  Conservation Test
    
      print *,'Conservation:'
      print *,'Grid1 Integral = ',sum(grid1_array*grid1_area*grid1_frac)
      print *,'Grid2 Integral = ',sum(grid2_array  *grid2_area*grid2_frac)
    
    !-----------------------------------------------------------------------
    ! - a first-order map from grid1 to grid2
  
    else if (map_type /= map_type_conserv .AND.map_type /= map_type_bicubic) then

      write(6,*) 'bilinear or conservative'
  
      call remap(grid2_array, wts_map1, grid2_add_map1, grid1_add_map1,grid1_array)
  
      print *,'First order mapping from grid1 to grid2:'
      print *,'----------------------------------------'
      print *,'Grid1 min,max: ',minval(grid1_array),maxval(grid1_array)
      print *,'Grid2 min,max: ',minval(grid2_array  ),maxval(grid2_array  )
    
      ! -  Conservation Test
    
      print *,'Conservation:'
      print *,'Grid1 Integral = ',sum(grid1_array*grid1_area*grid1_frac)
      print *,'Grid2 Integral = ',sum(grid2_array  *grid2_area*grid2_frac)
    
    !-----------------------------------------------------------------------
    ! - conservative mappings:
    ! - a second-order map from grid1 to grid2 with only lat grads
  
    else if (map_type == map_type_conserv .AND. lat_gradient) then
  
      call remap(grid2_array, wts_map1, grid2_add_map1, grid1_add_map1,  &
                 grid1_array, src_grad1=grad1_lat,src_grad2=grad1_lon_zero) 
  
      select case (norm_opt)
      case (norm_opt_none)
        grid2_tmp2 = grid2_frac*grid2_area
        where (grid2_tmp2 /= zero)
          grid2_array = grid2_array/grid2_tmp2
        elsewhere
          grid2_array = zero
        end where
      case (norm_opt_frcarea)
      case (norm_opt_dstarea)
        where (grid2_frac /= zero)
          grid2_array = grid2_array/grid2_frac
        elsewhere
          grid2_array = zero
        end where
      end select
  
      print *,'Second order mapping from grid1 to grid2 (lat only):'
      print *,'----------------------------------------'
      print *,'Grid1 min,max: ',minval(grid1_array),maxval(grid1_array)
      print *,'Grid2 min,max: ',minval(grid2_array  ),maxval(grid2_array  )
  
      ! -  Conservation Test
  
      print *,'Conservation:'
      print *,'Grid1 Integral = ',sum(grid1_array*grid1_area*grid1_frac)
      print *,'Grid2 Integral = ',sum(grid2_array*grid2_area*grid2_frac)
  
    !-----------------------------------------------------------------------
    ! - conservative mappings:
    ! - a second-order map from grid1 to grid2 both gradients
  
    else if (map_type == map_type_conserv .AND..NOT. lat_gradient) then
  
      call remap(grid2_array,wts_map1,grid2_add_map1,grid1_add_map1,  &
                 grid1_array, src_grad1=grad1_lat,src_grad2=grad1_lon) 
  
      select case (norm_opt)
      case (norm_opt_none)
        grid2_tmp2 = grid2_frac*grid2_area
        where (grid2_tmp2 /= zero)
          grid2_array = grid2_array/grid2_tmp2
        elsewhere
          grid2_array = zero
        end where
      case (norm_opt_frcarea)
      case (norm_opt_dstarea)
        where (grid2_frac /= zero)
          grid2_array = grid2_array/grid2_frac
        elsewhere
          grid2_array = zero
        end where
      end select
  
      print *,'Second order mapping from grid1 to grid2:'
      print *,'-----------------------------------------'
      print *,'Grid1 min,max: ',minval(grid1_array),maxval(grid1_array)
      print *,'Grid2 min,max: ',minval(grid2_array  ),maxval(grid2_array  )
  
      ! -  Conservation Test
  
      print *,'Conservation:'
      print *,'Grid1 Integral = ',sum(grid1_array*grid1_area*grid1_frac)
      print *,'Grid2 Integral = ',sum(grid2_array*grid2_area*grid2_frac)
  
    endif
  
    !-----------------------------------------------------------------------
    !     calculate some statistics
  
    grid2_count = zero
    grid2_tmp1  = zero
    grid2_tmp2  = zero
  
    print *,'number of sparse matrix entries ',num_links_map1
    do n=1,num_links_map1
      grid2_count(grid2_add_map1(n)) = grid2_count(grid2_add_map1(n)) + 1
      if (wts_map1(1,n) > one .or. wts_map1(1,n) < zero) then
        grid2_tmp1(grid2_add_map1(n)) = grid2_tmp1(grid2_add_map1(n)) + 1
        grid2_tmp2(grid2_add_map1(n)) = max(abs(wts_map1(1,n)),grid2_tmp2(grid2_add_map1(n)) )
      endif
    end do
  
    do n=1,grid2_size
      if (grid2_tmp1(n) > zero) print *,n,grid2_tmp2(n)
    end do
  
    imin = minval(grid2_count, mask=(grid2_count > 0))
    imax = maxval(grid2_count)
    idiff =  (imax - imin)/10 + 1
    print *,'total number of dest cells ',grid2_size
    print *,'number of cells participating in remap ',count(grid2_count > zero)
    print *,'min no of entries/row = ',imin
    print *,'max no of entries/row = ',imax
  
    imax = imin + idiff
    do n=1,10
      print *,'num of rows with entries between ',imin,' - ',imax-1,   &
               count(grid2_count >= imin .and. grid2_count < imax)
      imin = imin + idiff
      imax = imax + idiff
    end do
  
    !-----------------------------------------------------------------------
    ! - deallocate arrays
  
    deallocate (grid1_tmp, grad1_lat, grad1_lon,    &
                grad1_lat_zero, grad1_lon_zero, grid1_imask,     &
                grid2_tmp1, grid2_tmp2,   &
                grid2_imask, grid2_count)
  
  end subroutine calculate_grid

  ! ==========================================================================

end module scripinterp_mod

