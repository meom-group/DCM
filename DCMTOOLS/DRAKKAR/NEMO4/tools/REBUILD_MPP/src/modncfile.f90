MODULE modncfile
  IMPLICIT NONE
  PUBLIC
  TYPE  :: ncfile                         ! logical structure reflecting file structure
     INTEGER(KIND=4)                              :: ncid    ! file ncid
     INTEGER(KIND=4)                              :: ndims   ! number of dims
     INTEGER(KIND=4)                              :: nvars   ! number of vars
     INTEGER(KIND=4)                              :: natts   ! number of global attributes
     INTEGER(KIND=4)                              :: iunlim  ! ID of unlimited dimension
     INTEGER(KIND=4)                              :: kdimid  ! ID of vertical dimension
     INTEGER(KIND=4)                              :: npi     ! i-size of file
     INTEGER(KIND=4)                              :: npj     ! j-size of file
     INTEGER(KIND=4)                              :: npk     ! k-size of file
     INTEGER(KIND=4)                              :: npt     ! t-size of file
     INTEGER(KIND=4)                              :: npb     ! time_bound size
     INTEGER(KIND=4)                              :: idx     ! x dimid
     INTEGER(KIND=4)                              :: idy     ! y dimid
     INTEGER(KIND=4)                              :: idz     ! z dimid
     INTEGER(KIND=4)                              :: idt     ! t dimid
     INTEGER(KIND=4)                              :: idb     ! time bounds dimid
     INTEGER(KIND=4), DIMENSION(:),   ALLOCATABLE :: ideflat ! deflate level (nvar)
     INTEGER(KIND=4), DIMENSION(:),   ALLOCATABLE :: nvatt   ! number of att of each variable (var)
     INTEGER(KIND=4), DIMENSION(:),   ALLOCATABLE :: nvid    ! varid of each variable (var)
     INTEGER(KIND=4), DIMENSION(:),   ALLOCATABLE :: nvdim   ! dimension of each variable (var)
     INTEGER(KIND=4), DIMENSION(:),   ALLOCATABLE :: itype   ! type of each variable (var)
     INTEGER(KIND=4), DIMENSION(:),   ALLOCATABLE :: nlen    ! len of each dimension ( ndims)
     INTEGER(KIND=4), DIMENSION(:,:), ALLOCATABLE :: ichunk  ! size of chunk (nvar, ndims)
     INTEGER(KIND=4), DIMENSION(:,:), ALLOCATABLE :: idimids ! dimids of each variable (nvar, ndims) 
     CHARACTER(LEN=80)                            :: c_fnam  ! name of working file
     CHARACTER(LEN=80), DIMENSION(:), ALLOCATABLE :: c_vnam  ! name of each variable (var)
     CHARACTER(LEN=80), DIMENSION(:), ALLOCATABLE :: c_dnam  ! name of each dimension (ndims)
     LOGICAL,           DIMENSION(:), ALLOCATABLE :: lconti  ! contiguous flag (nvar)
     !   extra information for global attribute 
     INTEGER(KIND=4)                :: number_total          ! DOMAIN_number_total
     INTEGER(KIND=4)                :: number                ! DOMAIN_number
     INTEGER(KIND=4), DIMENSION(2)  :: idimensions_ids       ! DOMAIN_dimensions_ids
     INTEGER(KIND=4), DIMENSION(2)  :: isize_global          ! DOMAIN_size_global
     INTEGER(KIND=4), DIMENSION(2)  :: isize_local           ! DOMAIN_size_local
     INTEGER(KIND=4), DIMENSION(2)  :: iposition_first       ! DOMAIN_position_first
     INTEGER(KIND=4), DIMENSION(2)  :: iposition_last        ! DOMAIN_position_last 
     INTEGER(KIND=4), DIMENSION(2)  :: ihalo_size_start      ! DOMAIN_halo_size_start
     INTEGER(KIND=4), DIMENSION(2)  :: ihalo_size_end        ! DOMAIN_halo_size_end
     CHARACTER(LEN=80)              :: c_type                ! DOMAIN_type
  END TYPE ncfile

  TYPE(ncfile)                      :: sf_in  ! current working in structure
  TYPE(ncfile)                      :: sf_out ! array of output structure

  CHARACTER(LEN=80) :: cf_in                              ! input file name
  CHARACTER(LEN=80) :: cf_root                            ! root input file name (for merge)
  CHARACTER(LEN=80) :: cf_coor                            ! coordinate file name if used (for merge)
  CHARACTER(LEN=255) :: c_dirout='./'                     ! output directory

  CHARACTER(LEN=80), DIMENSION(:), ALLOCATABLE :: cf_list ! list of input <file*>_0000.nc

  LOGICAL           :: lg_coord      = .FALSE. ! flag for use of coordinate file (option)
  LOGICAL           :: lg_coord_read = .FALSE. ! flag set when the coordinate file is already read
  LOGICAL           :: lg_rename     = .FALSE. ! flag for renaming the output file to drakkar standards
  LOGICAL           :: lg_agrif      = .FALSE. ! flag for telling the prog that files are produced with agrif
  LOGICAL           :: lg_nc3        = .FALSE. ! flag for creating netcdf3 files instead of netcdf4
  LOGICAL           :: lg_verbose    = .FALSE. ! flag for increased verbosity
  LOGICAL           :: lg_win        = .FALSE. ! flag for screening output.
  LOGICAL           :: lg_kmax       = .FALSE. ! flag for limiting the vertical levels on output

  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dglam, dgphi   ! longitude, latitude from coordinate file
  INTEGER(KIND=4) :: mmpirank=0, nndone =1    ! dummy variables to mimic MPP behaviour



END MODULE modncfile
