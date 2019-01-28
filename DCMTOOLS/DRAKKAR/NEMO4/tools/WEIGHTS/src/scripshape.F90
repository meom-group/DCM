# define ARGC
      PROGRAM scripshape
!
! program to take output from the SCRIP weights generator
! and rearrange the data into a series of 2D fields suitable
! for reading with iom_get in NEMO configurations using the
! interpolation on the fly option
!
      USE netcdf
      IMPLICIT none
      INTEGER    :: ncId, VarId, status
      INTEGER    :: start(4), count(4)
      CHARACTER(LEN=1) :: y
      INTEGER    :: nd, ns, nl, nw, sx, sy, dx, dy
      INTEGER    :: i, j, k, m, n, smax
!
!  ifort -O2 -o scripshape scripshape.f90 \
!        -I/nerc/packages/netcdfifort/v3.6.0-pl1/include \
!        -L/nerc/packages/netcdfifort/v3.6.0-pl1/lib -lnetcdf
!
!
      INTEGER(KIND=4), ALLOCATABLE :: src(:)
      INTEGER(KIND=4), ALLOCATABLE :: dst(:)
      REAL(KIND=8), ALLOCATABLE :: wgt(:,:)
      REAL(KIND=8), ALLOCATABLE :: src1(:,:),dst1(:,:),wgt1(:,:)
      LOGICAL  :: around, verbose
#if defined ARGC
      INTEGER(KIND=4) :: iargc
      EXTERNAL :: iargc
#endif

      CHARACTER(LEN=256) :: interp_file, output_file, name_file
      INTEGER            :: ew_wrap
      NAMELIST /shape_inputs/ interp_file, output_file, ew_wrap

! scripshape requires 1 arguments; the name of the file containing
! the input namelist.
! This namelist contains:
!     the name of the input file containing the weights ! produced by SCRIP in its format;
!     the name of the new output file which ! is to contain the reorganized fields ready for input to NEMO.
!     the east-west wrapping of the input grid (-1, 0, 1 and 2 are accepted values)
!
! E.g.
!     interp_file = 'data_nemo_bilin.nc'
!     output_file = 'weights_bilin.nc'
!     ew_wrap     = 2
!
#if defined ARGC
      IF (iargc() == 1) THEN
        CALL getarg(1, name_file)
      ELSE
        WRITE(*,*) 'Usage: scripshape namelist_file'
        STOP
      ENDIF
#else
      WRITE(6,*) 'enter name of namelist file'
      READ(5,*) name_file
#endif
      interp_file = 'none'
      output_file = 'none'
      ew_wrap = 0
      OPEN(12, FILE=name_file, STATUS='OLD', FORM='FORMATTED')
      READ(12, NML=shape_inputs)
      CLOSE(12)
!
      INQUIRE(FILE = TRIM(interp_file), EXIST=around)
      IF (.not.around) THEN
       WRITE(*,*) 'Input file: '//TRIM(interp_file)//' not found'
       STOP
      ENDIF
! 
      INQUIRE(FILE = TRIM(output_file), EXIST=around)
      IF (around) THEN
       WRITE(*,*) 'Output file: '//TRIM(output_file)//' exists'
       WRITE(*,*) 'Ok to overwrite (y/n)?'
       READ(5,'(a)') y
       IF ( y .ne. 'y' .AND. y .ne. 'Y' ) STOP
      ENDIF
!
      verbose = .true.
!
! Obtain grid size information from interp_file
!
      CALL ncgetsize
!
! Allocate array spaces
!
      ALLOCATE(src(nl), STAT=status)
      IF(status /= 0 ) CALL alloc_err('src')
      ALLOCATE(dst(nl), STAT=status)
      IF(status /= 0 ) CALL alloc_err('dst')
      ALLOCATE(wgt(nw,nl), STAT=status)
      IF(status /= 0 ) CALL alloc_err('wgt')
      ALLOCATE(src1(dx,dy), STAT=status)
      IF(status /= 0 ) CALL alloc_err('src1')
      ALLOCATE(dst1(dx,dy), STAT=status)
      IF(status /= 0 ) CALL alloc_err('dst1')
      ALLOCATE(wgt1(dx,dy), STAT=status)
      IF(status /= 0 ) CALL alloc_err('wgt1')
!
! Read all required data from interp_file
!
      CALL ncgetfields
!
! Check that dst is monotonically increasing
!
      DO k = 1,nl-1
       IF(dst(k+1).lt.dst(k)) THEN
        WRITE(*,*) 'non-monotonic at ',k
        WRITE(*,*) dst(k-4:k+16)
        STOP
       ENDIF
      END DO
!
! Remove references to the top row of src
!
      IF(verbose) WRITE(*,*) &
        'Removing references to the top row of the source grid'
      smax = (sy-1)*sx
      n = 0
      DO k = 1,nl
       IF(src(k).gt.smax-1) THEN
        src(k) = src(k)-sx
        n = n + 1
       ENDIF
      END DO
      IF(verbose) WRITE(*,*) n,' values changed (',100.*n/nl,'%)'
!
! Loop through weights for each corner in turn and
! rearrange weight fields into separate 2D fields for
! reading with iom_get in NEMO
!
      DO k = 1,nw
       DO n = 1,4

        i = 0
        j = 1
        DO m = n,nl,4
         i = i+1
         IF(i.gt.dx) THEN
          i = 1
          j = j + 1
         ENDIF
         src1(i,j) = src(m)
         dst1(i,j) = dst(m)
         wgt1(i,j) = wgt(k,m)
        END DO
!
! Write out this set which will be labelled with
! a 2 digit number equal to n+4*(k-1)
!
        CALL wrwgts
!   
       END DO
      END DO
      STOP
      CONTAINS
!
!----------------------------------------------------------------------*
      SUBROUTINE ncgetsize
!
! Access grid size information in interp_file and set the
! following integers:
!
!    nd = dst_grid_size
!    ns = src_grid_size
!    nl = num_links    
!    nw = num_wgts     
! sx,sy = src_grid_dims     
! dx,dy = dst_grid_dims     
!
      INTEGER idims(2)
!
       status = nf90_open(interp_file, nf90_NoWrite, ncid)
       IF(status /= nf90_NoErr) CALL handle_err(status)
!
       status = nf90_inq_dimid(ncid, 'dst_grid_size', VarId)
       IF(status /= nf90_NoErr) CALL handle_err(status)
       status = nf90_inquire_dimension(ncid, VarId, LEN = nd)
       IF(status /= nf90_NoErr) CALL handle_err(status)
!
       status = nf90_inq_dimid(ncid, 'src_grid_size', VarId)
       IF(status /= nf90_NoErr) CALL handle_err(status)
       status = nf90_inquire_dimension(ncid, VarId, LEN = ns)
       IF(status /= nf90_NoErr) CALL handle_err(status)
!
       status = nf90_inq_dimid(ncid, 'num_links', VarId)
       IF(status /= nf90_NoErr) CALL handle_err(status)
       status = nf90_inquire_dimension(ncid, VarId, LEN = nl)
       IF(status /= nf90_NoErr) CALL handle_err(status)
!
       status = nf90_inq_dimid(ncid, 'num_wgts', VarId)
       IF(status /= nf90_NoErr) CALL handle_err(status)
       status = nf90_inquire_dimension(ncid, VarId, LEN = nw)
       IF(status /= nf90_NoErr) CALL handle_err(status)
!
       start = 1
       count = 2
       status = nf90_inq_varid(ncid, 'src_grid_dims', VarId)
       IF(status /= nf90_NoErr) CALL handle_err(status)
       status = nf90_get_var(ncid, VarId, idims, start, count)
       IF(status /= nf90_NoErr) CALL handle_err(status)
       sx = idims(1) ; sy = idims(2)
!
       status = nf90_inq_varid(ncid, 'dst_grid_dims', VarId)
       IF(status /= nf90_NoErr) CALL handle_err(status)
       status = nf90_get_var(ncid, VarId, idims, start, count)
       IF(status /= nf90_NoErr) CALL handle_err(status)
       dx = idims(1) ; dy = idims(2)
!
       status = nf90_close(ncid)
       IF (status /= nf90_noerr) CALL handle_err(status)
!
       IF(verbose) THEN
         WRITE(*,*) 'Detected sizes: '
         WRITE(*,*) 'dst_grid_size: ', nd
         WRITE(*,*) 'src_grid_size: ', ns
         WRITE(*,*) 'num_links    : ', nl
         WRITE(*,*) 'num_wgts     : ', nw
         WRITE(*,*) 'src_grid_dims: ', sx, ' x ', sy
         WRITE(*,*) 'dst_grid_dims: ', dx, ' x ', dy
       ENDIF
!
      END SUBROUTINE ncgetsize

!----------------------------------------------------------------------*
      SUBROUTINE ncgetfields
!
! Read all required data from interp_file. The data read are:
!
! netcdf variable    size   internal array
!-----------------+-------+--------------
! src_address        nl     src
! dst_address        nl     dst
! remap_matrix     (nw,nl)  wgt
!
       status = nf90_open(interp_file, nf90_NoWrite, ncid)
       IF(status /= nf90_NoErr) CALL handle_err(status)
!
       status = nf90_inq_varid(ncid, 'src_address', VarId)
       IF(status /= nf90_NoErr) CALL handle_err(status)
!
! Read the values for src
       status = nf90_get_var(ncid, VarId, src, &
                         start = (/ 1 /),      &
                         count = (/ nl /))
       IF(status /= nf90_NoErr) CALL handle_err(status)
!
       status = nf90_inq_varid(ncid, 'dst_address', VarId)
       IF(status /= nf90_NoErr) CALL handle_err(status)
!
! Read the values for dst
       status = nf90_get_var(ncid, VarId, dst, &
                         start = (/ 1 /),      &
                         count = (/ nl /))
       IF(status /= nf90_NoErr) CALL handle_err(status)
!
       status = nf90_inq_varid(ncid, 'remap_matrix', VarId)
       IF(status /= nf90_NoErr) CALL handle_err(status)
!
! Read the values for wgt
       status = nf90_get_var(ncid, VarId, wgt, &
                         start = (/ 1, 1 /),   &
                         count = (/ nw, nl /))
       IF(status /= nf90_NoErr) CALL handle_err(status)
!
       status = nf90_close(ncid)
       IF (status /= nf90_noerr) CALL handle_err(status)
!
      END SUBROUTINE ncgetfields

!----------------------------------------------------------------------*
      SUBROUTINE handle_err(status)
!
! Simple netcdf error checking routine
!
      INTEGER, intent ( in) :: status
!
      IF(status /= nf90_noerr) THEN
       IF(trim(nf90_strerror(status)) .eq. 'Attribute not found') THEN
! ignore
       ELSE
        WRITE(*,*) trim(nf90_strerror(status))
        STOP "Stopped"
       END IF
      END IF
      END SUBROUTINE handle_err

!----------------------------------------------------------------------*
      SUBROUTINE alloc_err(arname)
!
! Simple allocation error checking routine
!
      CHARACTER(LEN=*) :: arname
!
      WRITE(*,*) 'Allocation error attempting to ALLOCATE '//arname
      STOP "Stopped"
      END SUBROUTINE alloc_err

!                                                                       
!----------------------------------------------------------------------*
      SUBROUTINE wrwgts
!
! Write out each set of 2D fields to output_file.
! Each call will write out a set of srcXX, dstXX and wgtXX fields
! where XX is a two digit number equal to n + 4*(k-1). The first
! and last calls to this routine initialise and close the output
! file respectively. The first call is detected when k*n=1 and the
! last call is detected when k*n=4*nw. The outfile file remains
! open between the first and last calls. 
!
      INTEGER :: status, ncid, ncin
      INTEGER :: Lontdid, Lattdid
      INTEGER :: tvid, tvid2, tvid3
      INTEGER :: ioldfill
      CHARACTER(LEN=2) :: cs
      SAVE ncid, Lontdid, Lattdid
!
      IF(k*n.eq.1) THEN
!
! Create output_file and set the dimensions
!
         status = nf90_create(output_file, or(nf90_Clobber,nf90_netcdf4), ncid)
         IF(status /= nf90_NoErr) CALL handle_err(status)
         status = nf90_set_fill(ncid, nf90_NoFill, ioldfill)
         IF(status /= nf90_NoErr) CALL handle_err(status)
!
         status = nf90_def_dim(ncid, "lon", dx, Lontdid)
         IF(status /= nf90_NoErr) CALL handle_err(status)
         status = nf90_def_dim(ncid, "lat", dy, Lattdid)
         IF(status /= nf90_NoErr) CALL handle_err(status)
!
         status = nf90_put_att(ncid, nf90_global, 'ew_wrap', ew_wrap)
         IF(status /= nf90_NoErr) CALL handle_err(status)
      ELSE
!
! Reenter define mode
!
         status = nf90_redef(ncid)
         IF(status /= nf90_NoErr) CALL handle_err(status)
      ENDIF
!
      WRITE(cs,'(i2.2)') n + 4*(k-1)
!
! Define new variables
!
      status = nf90_def_var(ncid, "src"//cs, nf90_double, &
                         (/ Lontdid, Lattdid /), tvid)
      IF(status /= nf90_NoErr) CALL handle_err(status)
!
      status = nf90_def_var(ncid, "dst"//cs, nf90_double, &
                         (/ Lontdid, Lattdid /), tvid2)
      IF(status /= nf90_NoErr) CALL handle_err(status)
!
      status = nf90_def_var(ncid, "wgt"//cs, nf90_double, &
                         (/ Lontdid, Lattdid /), tvid3)
      IF(status /= nf90_NoErr) CALL handle_err(status)
!
! Leave define mode
!
      status = nf90_enddef(ncid)
      IF(status /= nf90_NoErr) CALL handle_err(status)
!
! Write the data
!
      status = nf90_put_var(ncid, tvid, src1,    &
                            start = (/ 1, 1 /),  &
                            count = (/ dx, dy /) )
      IF(status /= nf90_NoErr) CALL handle_err(status)
!
      status = nf90_put_var(ncid, tvid2, dst1,   &
                            start = (/ 1, 1 /),  &
                            count = (/ dx, dy /) )
      IF(status /= nf90_NoErr) CALL handle_err(status)
!
      status = nf90_put_var(ncid, tvid3, wgt1,   &
                            start = (/ 1, 1  /), &
                            count = (/ dx, dy /) )
      IF(status /= nf90_NoErr) CALL handle_err(status)
!
      IF(k*n.eq.4*nw) THEN
!
! --     Reenter define mode
!
         status = nf90_redef(ncid)
         IF(status /= nf90_NoErr) CALL handle_err(status)
!
! --     Reopen interp_file and transfer some global attributes
!
         status = nf90_open(interp_file, nf90_NoWrite, ncin)
         IF(status /= nf90_NoErr) CALL handle_err(status)
!
         status = nf90_copy_att(ncin,NF90_GLOBAL,'title',        ncid,NF90_GLOBAL)
         IF(status /= nf90_NoErr) CALL handle_err(status)
!
         status = nf90_copy_att(ncin,NF90_GLOBAL,'normalization',ncid,NF90_GLOBAL)
         IF(status /= nf90_NoErr) CALL handle_err(status)
!
         status = nf90_copy_att(ncin,NF90_GLOBAL,'map_method',   ncid,NF90_GLOBAL)
         IF(status /= nf90_NoErr) CALL handle_err(status)
!
         status = nf90_copy_att(ncin,NF90_GLOBAL,'conventions',  ncid,NF90_GLOBAL)
         IF(status /= nf90_NoErr) CALL handle_err(status)
!
         status = nf90_copy_att(ncin,NF90_GLOBAL,'history',      ncid,NF90_GLOBAL)
         IF(status /= nf90_NoErr) CALL handle_err(status)
!
! --     Close interp_file
!
         status = nf90_close(ncin)
         IF(status /= nf90_NoErr) CALL handle_err(status)
!
! --     Close output_file
!
         status = nf90_close(ncid)
         IF(status /= nf90_NoErr) CALL handle_err(status)
      ENDIF

      END SUBROUTINE wrwgts
      END PROGRAM scripshape
