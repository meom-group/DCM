MODULE modutil
  USE netcdf
  USE modncfile
  IMPLICIT NONE
  !
  INTEGER(KIND=4), PUBLIC :: mini,minj,maxi,maxj   ! coordinates for the zoom
  INTEGER(KIND=4), PUBLIC :: maxk=32767            ! limitation for the level output (high default)
  INTEGER(KIND=4), PUBLIC :: nbloc=1               ! number of blocks files
  INTEGER(KIND=4), PUBLIC :: ndigit=4              ! number of digit for coding rank

  CHARACTER(LEN=20), PUBLIC :: c_pattern='_0000.nc'

CONTAINS
  SUBROUTINE RebuildNc ( kindex ) 
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE RebuildNc  ***
    !!
    !! ** Purpose :  Perfom rebuild on the current kindex file 
    !!
    !! ** Method  :   
    !!
    !! References :  
    !!----------------------------------------------------------------------
    INTEGER(KIND=4), INTENT(in) :: kindex 

    INTEGER(KIND=4)    :: ipos, ierr
    CHARACTER(LEN=255) :: cl_tmp

    !!----------------------------------------------------------------------
    ! rebuild merged file
    cf_in = cf_list (kindex )
    ipos    = INDEX(cf_in,TRIM(c_pattern), .true. )
    IF ( ipos == 0 ) THEN
       PRINT *,' ERROR : pattern ',TRIM(c_pattern),' not found in file name.'
       STOP 99
    ENDIF
    cf_root = cf_in(1:ipos-1)
    sf_in   = GetNcFile ( cf_in )

    IF ( lg_coord_read ) THEN   ! Read coordinates to fill holes in nav_lon, nav_lat
       ierr = GetCoord ()
       lg_coord_read = .FALSE.  ! read once for all
    ENDIF

    IF ( lg_verbose ) THEN
       PRINT 9000,  mmpirank, nndone, 'FILE_npk :'              , sf_in%npk
       PRINT 9000 , mmpirank, nndone, 'DOMAIN_number_total: '   , sf_in%number_total
       PRINT 9000 , mmpirank, nndone, 'DOMAIN_number: '         , sf_in%number
       PRINT 9001 , mmpirank, nndone, 'DOMAIN_dimensions_ids: ' , sf_in%idimensions_ids(:)
       PRINT 9001 , mmpirank, nndone, 'DOMAIN_size_global: '    , sf_in%isize_global(:)
       PRINT 9001 , mmpirank, nndone, 'DOMAIN_size_local:  '    , sf_in%isize_local(:)
       PRINT 9001 , mmpirank, nndone, 'DOMAIN_position_first: ' , sf_in%iposition_first(:)
       PRINT 9001 , mmpirank, nndone, 'DOMAIN_position_last: '  , sf_in%iposition_last(:)
       PRINT 9001 , mmpirank, nndone, 'DOMAIN_halo_size_start: ', sf_in%ihalo_size_start(:)
       PRINT 9001 , mmpirank, nndone, 'DOMAIN_halo_size_end: '  , sf_in%ihalo_size_end(:)
       PRINT 9002 , mmpirank, nndone, 'DOMAIN_type: '           , TRIM(sf_in%c_type)
9000   FORMAT( I5.5,1x,i5, 1x, a22,i5)
9001   FORMAT( I5.5,1x,i5, 1x, a22,i5,1x,i5)
9002   FORMAT( I5.5,1x,i5, 1x, a22,a)
    ENDIF

    ALLOCATE ( sf_out%nvdim   (sf_in%nvars) )
    ALLOCATE ( sf_out%nvid    (sf_in%nvars) )
    ALLOCATE ( sf_out%c_vnam  (sf_in%nvars) )
    ALLOCATE ( sf_out%nvatt   (sf_in%nvars) )
    ALLOCATE ( sf_out%itype   (sf_in%nvars) )
    ALLOCATE ( sf_out%c_dnam  (sf_in%ndims) )
    ALLOCATE ( sf_out%nlen    (sf_in%ndims) )
    ALLOCATE ( sf_out%idimids (sf_in%nvars, sf_in%ndims) )

    sf_out                    = sf_in    ! copy whole stucture at once then adjust
    sf_out%number_total       = 1        ! one single domain at the end
    sf_out%number             = -1       ! flag value
    sf_out%iposition_first(:) = (/1                    , 1                    /)
    sf_out%iposition_last(:)  = (/sf_in%isize_global(1), sf_in%isize_global(2)/)
    sf_out%isize_local(:)     = sf_out%isize_global
    sf_out%npi                = sf_in%isize_global(1)
    sf_out%npj                = sf_in%isize_global(2)

    IF ( lg_rename ) THEN
       ! transform <CONFIG-CASE>_<freq>_<typ>_yyyymmdd-YYYYMMDD to  <CONFIG-CASE>_yYYYYmMMdDD.<freq>_<typ>.nc
       sf_out%c_fnam = TRIM(RenameFile( cf_root ) )   !  c_dirout already included in renamed file
    ELSE
       sf_out%c_fnam = TRIM(c_dirout)//'/'//TRIM(cf_root)//'_merg.nc'
    ENDIF

    sf_out%nlen(sf_in%idimensions_ids(1))     = sf_in%isize_global(1)
    sf_out%nlen(sf_in%idimensions_ids(2))     = sf_in%isize_global(2)
    IF ( lg_kmax ) THEN
    ! limit vertical dimension to maxk
       IF ( sf_in%ndims >= 3 ) THEN  ! likely the 3 rd dim is vertical ( or ice cat ? )
         sf_out%nlen(sf_out%kdimid) = MIN(maxk, sf_out%nlen(sf_out%kdimid) )
       ENDIF
    ENDIF
print *, 'JMMMM : ', sf_out%nlen(sf_out%kdimid)
    sf_out%ihalo_size_start(:) = 0
    sf_out%ihalo_size_end(:)   = 0

    CALL CreateMergedFile ( sf_out                     )
    CALL FillMergedFile   ( sf_out, sf_in%number_total )
    DEALLOCATE ( sf_out%nvdim, sf_out%nvid, sf_out%c_vnam, sf_out%nvatt, &
         &         sf_out%itype, sf_out%c_dnam, sf_out%nlen, sf_out%idimids )

  END SUBROUTINE RebuildNc

  FUNCTION GetNcFile (cd_file)
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION GetNcFile  ***
    !!
    !! ** Purpose :  fills in the ncfile  structure corresponding to the file
    !!               given in argument 
    !!
    !! ** Method  :  Use NF90 function to get the ad-hoc information 
    !!
    !!----------------------------------------------------------------------
    CHARACTER(LEN=*), INTENT(in) :: cd_file
    TYPE(ncfile)                 :: GetNcFile

    INTEGER(KIND=4) :: jvar, jdim                      ! loop index
    INTEGER(KIND=4) :: ierr, idx, idy, idz, idt, idb   ! error status and dimids
    INTEGER(KIND=4) :: ilkdim=0                        ! dim index of k-dimension
    !!----------------------------------------------------------------------
    GetNcFile%c_fnam = cd_file
    ! Open nc file
    ierr = NF90_OPEN(cd_file, NF90_NOWRITE, GetNcFile%ncid )
    IF ( ierr /= NF90_NOERR ) print *, TRIM( NF90_STRERROR (ierr)), 'OPEN'
    ! inquire information about dataset
    ierr = NF90_INQUIRE(GetNcFile%ncid, nDimensions    = GetNcFile%ndims,  &
         &                              nVariables     = GetNcFile%nvars,  &
         &                              nAttributes    = GetNcFile%natts,  &
         &                              unlimitedDimId = GetNcFile%iunlim  )
    ! Allocate space
    ALLOCATE (GetNcFile%nvdim   (GetNcFile%nvars) )
    ALLOCATE (GetNcFile%nvid    (GetNcFile%nvars) )
    ALLOCATE (GetNcFile%c_vnam  (GetNcFile%nvars) )
    ALLOCATE (GetNcFile%nvatt   (GetNcFile%nvars) )
    ALLOCATE (GetNcFile%itype   (GetNcFile%nvars) )
    ALLOCATE (GetNcFile%lconti  (GetNcFile%nvars) )
    ALLOCATE (GetNcFile%ideflat (GetNcFile%nvars) )
    ALLOCATE (GetNcFile%c_dnam(GetNcFile%ndims) )
    ALLOCATE (GetNcFile%nlen    (GetNcFile%ndims) )
    ALLOCATE (GetNcFile%idimids (GetNcFile%nvars,GetNcFile%ndims) )
    ALLOCATE (GetNcFile%ichunk  (GetNcFile%nvars,GetNcFile%ndims) )

    ! Look for dimensions
    DO jdim = 1, GetNcFile%ndims
       ierr = NF90_INQUIRE_DIMENSION(GetNcFile%ncid,jdim,                 &
            &                         name   = GetNcFile%c_dnam(jdim),    &
            &                         len    = GetNcFile%nlen  (jdim) )
    ENDDO
    ! Look for variables
    DO jvar = 1, GetNcFile%nvars
       ierr = NF90_INQUIRE_VARIABLE (GetNcFile%ncid, jvar,                &
            &                         name   = GetNcFile%c_vnam(jvar),    &
            &                         xtype  = GetNcFile%itype(jvar),     &
            &                         ndims  = GetNcFile%nvdim(jvar),     &
            &                         dimids = GetNcFile%idimids(jvar,:), &
            &                         nAtts  = GetNcFile%nvatt(jvar)) ! ,     &
       !           &                         contiguous = GetNcFile%lconti(jvar),   &
       !           &                         chunksizes = GetNcFile%ichunk(jvar,:), &
       !           &                         deflate_level = GetNcFile%ideflat(jvar)      )
       IF ( GetNcFile%nvdim(jvar) == 4 ) THEN
          ! likely, the 3rd dimension of this variable is vertical
          ilkdim=GetNcFile%idimids(jvar,3)
          GetNcFile%npk = GetNcFile%nlen  (ilkdim)
          GetNcFile%kdimid = ilkdim
       ENDIF
    END DO
    ! Look for attributes
    ierr = NF90_GET_ATT (GetNcFile%ncid, NF90_GLOBAL, 'DOMAIN_number_total'   , GetNcFile%number_total       )
    ierr = NF90_GET_ATT (GetNcFile%ncid, NF90_GLOBAL, 'DOMAIN_number      '   , GetNcFile%number             )
    ierr = NF90_GET_ATT (GetNcFile%ncid, NF90_GLOBAL, 'DOMAIN_dimensions_ids' , GetNcFile%idimensions_ids(:) )
    ierr = NF90_GET_ATT (GetNcFile%ncid, NF90_GLOBAL, 'DOMAIN_size_global'    , GetNcFile%isize_global(:)    )
    ierr = NF90_GET_ATT (GetNcFile%ncid, NF90_GLOBAL, 'DOMAIN_size_local'     , GetNcFile%isize_local(:)     )
    ierr = NF90_GET_ATT (GetNcFile%ncid, NF90_GLOBAL, 'DOMAIN_position_first' , GetNcFile%iposition_first(:) )
    ierr = NF90_GET_ATT (GetNcFile%ncid, NF90_GLOBAL, 'DOMAIN_position_last'  , GetNcFile%iposition_last(:)  )
    ierr = NF90_GET_ATT (GetNcFile%ncid, NF90_GLOBAL, 'DOMAIN_halo_size_start', GetNcFile%ihalo_size_start(:))
    ierr = NF90_GET_ATT (GetNcFile%ncid, NF90_GLOBAL, 'DOMAIN_halo_size_end'  , GetNcFile%ihalo_size_end(:)  )
    ierr = NF90_GET_ATT (GetNcFile%ncid, NF90_GLOBAL, 'DOMAIN_type'           , GetNcFile%c_type             )

    ! DOMAIN_dimensions_ids gives ids for x, y 
    idx = GetNcFile%idimensions_ids(1)
    idy = GetNcFile%idimensions_ids(2)
    ! time is unlimited dim
    idt = GetNcFile%iunlim

    ! try to infer size of the domain assuming some basis:
    ! (1) 2D var are (x,y)
    ! (2) time dim is unlimited
    ! (3) allowed shape of var : x,y ; x,y,t ; x,y,z,t  ; x,y,z   ; t   ; z 

    ! Look for   bound dim id and other 3rd dimension ( z, or nicecat or icbcla, icbsect ... ) 
    idz=-1; idb=-1 
    ierr = NF90_INQ_DIMID(GetNcFile%ncid,'axis_nbounds',idb)   ! hard coded dimension name

    ! remaining dimensions are either z, nicecat or icbsect, icbcla
    !  note that there might be more than 1 3rd dimension ( case of ICB for instance)

    DO jvar = 1, GetNcFile%nvars
       GetNcFile%nvid(jvar) = jvar
    ENDDO

    GetNcFile%idx= GetNcFile%idimensions_ids(1)
    GetNcFile%idy= GetNcFile%idimensions_ids(2)
    GetNcFile%idt= GetNcFile%iunlim
    GetNcFile%idb= idb

    IF ( idx == -1 .OR. idy == -1 ) THEN 
       PRINT *, ' ERROR : no x, y dimensions found'
       STOP
    ENDIF

    ! fill in DOMAIN attributes 
    ! move above ...
    ! CORRECT for Halo :  ???
    GetNcfile%isize_local(:)     = GetNcfile%isize_local(:)     - GetNcFile%ihalo_size_start(:) - GetNcFile%ihalo_size_end(:)
    GetNcFile%iposition_first(:) = GetNcFile%iposition_first(:) + GetNcFile%ihalo_size_start(:)
    GetNcFile%iposition_last(:)  = GetNcFile%iposition_last(:)  - GetNcFile%ihalo_size_end(:)

    GetNcFile%npi       = GetNcfile%nlen(idx)
    GetNcFile%npj       = GetNcfile%nlen(idy)

  END FUNCTION GetNcFile

  FUNCTION RenameFile (cd_froot )
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION RenamFile  ***
    !!
    !! ** Purpose :  Rename NEMO output file from XIOS to DRAKKAR standards 
    !!
    !! ** Method  :  transform <CONFIG-CASE>_<freq>_<typ>_yyyymmdd-YYYYMMDD 
    !!                    to   <CONFIG-CASE>_yYYYYmMMdDD.<freq>_<typ>
    !!               update globale variable c_dirout to <freq>_OUTPUT
    !!
    !!----------------------------------------------------------------------
    CHARACTER(LEN=*), INTENT(in ) :: cd_froot
    CHARACTER(LEN=255)            :: RenameFile

    INTEGER(KIND=4)    :: ipos         ! for searching _ or - in root name
    CHARACTER(LEN=80 ) :: cl_confcase  ! hold the config-case part of the name
    CHARACTER(LEN=80 ) :: cl_freq      ! hold the freq  part of the name
    CHARACTER(LEN=80 ) :: cl_typ       ! hold the typ  part of the name
    CHARACTER(LEN=80 ) :: cl_tagi      ! hold the tag  part of the name in c_froot
    CHARACTER(LEN=80 ) :: cl_tago      ! hold the tag  part of the name in renamed file
    CHARACTER(LEN=255) :: cl_dum       ! working copy of cd_froot
    !!----------------------------------------------------------------------
    cl_dum = cd_froot

    ! look for config-case
    ipos        = INDEX( cl_dum,'_')
    IF ( lg_agrif) THEN  ! skip 1_ or 2_ 
      cl_dum      = cl_dum(ipos+1: )
      ipos        = INDEX( cl_dum,'_')
    ENDIF
     
    cl_confcase = cl_dum(1:ipos-1)
    cl_dum      = cl_dum(ipos+1: )

    ! look for freq
    ipos        = INDEX( cl_dum,'_')
    cl_freq     = cl_dum(1:ipos-1)
    cl_dum      = cl_dum(ipos+1: )
    c_dirout    = TRIM(cl_freq)//'_OUTPUT'
    CALL SYSTEM ("mkdir -p "//TRIM(c_dirout) )

    ! look for typ
    ipos        = INDEX( cl_dum,'_')
    cl_typ      = cl_dum(1:ipos-1)
    cl_dum      = cl_dum(ipos+1: )

    ! look for second date in the remainder of cl_dum
    ipos       = INDEX( cl_dum,'-')
    cl_tagi    = cl_dum(ipos+1 : )

    IF ( TRIM(cl_tagi(5:6)) /= '' ) THEN
      IF ( TRIM(cl_tagi(7:8)) /= '' ) THEN
         cl_tago    = 'y'//TRIM(cl_tagi(1:4))//'m'//TRIM(cl_tagi(5:6))//'d'//TRIM(cl_tagi(7:8))
      ELSE
         cl_tago    = 'y'//TRIM(cl_tagi(1:4))//'m'//TRIM(cl_tagi(5:6))
      ENDIF
    ELSE 
      cl_tago   = 'y'//TRIM(cl_tagi(1:4))
    ENDIF

    ! build the final name
    RenameFile = TRIM(c_dirout)//'/'//TRIM(cl_confcase)//'_'//TRIM(cl_tago)//'.'//TRIM(cl_freq)//'_'//TRIM(cl_typ)//'.nc'

    IF ( lg_verbose) PRINT '(I4.4,I4,1x,a)', mmpirank, nndone, TRIM(cf_root)//' renamed to '//TRIM(RenameFile)

  END FUNCTION RenameFile

  SUBROUTINE CreateMergedFile ( sd_fout )
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE CreateMergedFile  ***
    !!
    !! ** Purpose :  Create the structure of the output file according to
    !!               the structure passed in argument 
    !!
    !! ** Method  :  Just use NF90 functions, in the right order
    !!
    !!----------------------------------------------------------------------
    TYPE(ncfile),    INTENT(inout) :: sd_fout   ! output file structure

    INTEGER(KIND=4)                :: jdim, jvar, jatt  ! loop index
    INTEGER(KIND=4)                :: ierr              ! ncdf error status
    INTEGER(KIND=4)                :: id, iold          ! var/dim id and working integer
    INTEGER(KIND=4)                :: ipos, ildim       ! 
    INTEGER(KIND=4)                :: idimv
    INTEGER(KIND=4)                :: ideflate
    INTEGER(KIND=4), DIMENSION(4)  :: ichunk
    CHARACTER(LEN=256)             :: cl_att
    !!----------------------------------------------------------------------
    ! create file
    IF ( lg_nc3 ) THEN
       ierr = NF90_CREATE  (sd_fout%c_fnam, cmode=or(NF90_CLOBBER,NF90_64BIT_OFFSET), ncid=sd_fout%ncid )
    ELSE
       ierr = NF90_CREATE  (sd_fout%c_fnam, cmode=or(NF90_CLOBBER,NF90_NETCDF4), ncid=sd_fout%ncid )
    ENDIF
    ierr = NF90_SET_FILL(sd_fout%ncid, NF90_NOFILL, iold)

    ! create dimensions
    DO jdim = 1, sd_fout%ndims
       ildim = sd_fout%nlen(jdim)
       IF ( jdim == sd_fout%iunlim )  ildim = NF90_UNLIMITED     ! trick to set UNLIMITED time_counter dimension  
       ierr = NF90_DEF_DIM(sd_fout%ncid, sd_fout%c_dnam(jdim), ildim, id)
    ENDDO

    ! create variables
    DO jvar = 1, sd_fout%nvars
       idimv=sd_fout%nvdim(jvar)
       SELECT CASE ( idimv )
       CASE ( 1 )
          IF ( sd_fout%idimids(jvar,1) == sd_fout%iunlim ) THEN
             ! time counter variable
             ichunk(1) = 1
             ideflate = 0
          ELSE
             ! depth like variable : related to 3rd dimension (can be depth but also nicecat, icbsect, icbcla ...)
             ichunk(1) = sd_fout%nlen(sd_fout%idimids(jvar,1))
             ideflate = 0
          ENDIF
       CASE ( 2 )
          IF ( sd_fout%idimids(jvar,1) == sd_fout%idb ) THEN
             ! axis_bound, time variable or depth
             ichunk(1:2)=(/sd_fout%nlen(sd_fout%idb),1/)
             ideflate = 0
          ELSE IF (sd_fout%idimids(jvar,2) == sd_fout%iunlim ) THEN  ! z,t  variables ? or x,t or y t ? 
             ichunk(1:2)=(/sd_fout%nlen(sd_fout%idimids(jvar,1)),1/)
             ideflate = 0
          ELSE
             ! x, y variables ( nav_lon, nav_lat )
             ichunk(1:2)=sf_in%isize_local(:)
             ideflate = 0
          ENDIF
       CASE ( 3 )
          ! x, y, t 
          ichunk(1:3)=(/sf_in%isize_local(:),1 /)
          ideflate = 1
       CASE ( 4 )
          ! x, y, z, t
          ichunk=(/sf_in%isize_local(:), 1, 1 /)
          ideflate = 1
       CASE DEFAULT 
          ! kt, rdt
          ichunk=1
          ideflate = 0
       END SELECT

       IF ( lg_nc3 ) THEN
          ierr = NF90_DEF_VAR( sd_fout%ncid,                                &
               &             sd_fout%c_vnam(jvar),                          &
               &             sd_fout%itype(jvar),                           &
               &             sd_fout%idimids(jvar,1:sd_fout%nvdim(jvar) ) , &
               &             sd_fout%nvid(jvar)                             )
       ELSE
          !        print *,'vnam : ', TRIM(sd_fout%c_vnam(jvar))
          !        print *,' typ : ', sd_fout%itype(jvar)
          !        print *,'idim : ', sd_fout%idimids(jvar,1:sd_fout%nvdim(jvar) )
          !        print *,'chk  : ', ichunk(1:idimv)
          !        print *,'defl : ', ideflate
          !        print *,'idimv: ', idimv
          !        print *, 'SIZE:', size(  sd_fout%idimids(jvar,1:sd_fout%nvdim(jvar) ))
          IF ( idimv /= 0 ) THEN
             ierr = NF90_DEF_VAR( sd_fout%ncid,                                &
                  &             sd_fout%c_vnam(jvar),                          &
                  &             sd_fout%itype(jvar),                           &
                  &             sd_fout%idimids(jvar,1:sd_fout%nvdim(jvar) ) , &
                  &             sd_fout%nvid(jvar),                            &
                  &             chunksizes=ichunk(1:idimv),                    &
                  &             deflate_level=ideflate                         )
          ELSE
             ierr = NF90_DEF_VAR( sd_fout%ncid,                                &
                  &             sd_fout%c_vnam(jvar),                          &
                  &             sd_fout%itype(jvar),                           &
                  &             sd_fout%nvid(jvar))
          ENDIF
       ENDIF

       ! copy attribute
       DO jatt = 1, sd_fout%nvatt(jvar) 
          ierr = NF90_INQ_ATTNAME( sf_in%ncid, sf_in%nvid(jvar),jatt, cl_att) 
          ierr = NF90_COPY_ATT( sf_in%ncid, sf_in%nvid(jvar), cl_att, sd_fout%ncid, sd_fout%nvid(jvar))
       ENDDO
    ENDDO

    ! copy Global attributes
    DO jatt = 1, sf_in%natts
       ierr = NF90_INQ_ATTNAME( sf_in%ncid, NF90_GLOBAL, jatt, cl_att) 
       IF (   cl_att == 'name'      .OR.   &
            & cl_att == 'timeStamp' .OR.   &
            & cl_att == 'ibegin'    .OR.   &
            & cl_att == 'jbegin'    .OR.   &
            & cl_att == 'ni'        .OR.   &
            & cl_att == 'nj'               ) CYCLE
       IF ( INDEX(cl_att,'DOMAIN') == 0 ) THEN
          ierr = NF90_COPY_ATT( sf_in%ncid, NF90_GLOBAL, cl_att, sd_fout%ncid, NF90_GLOBAL)
       ENDIF
    ENDDO

    ! End define mode
    ierr = NF90_ENDDEF ( sd_fout%ncid )
    print *,'ENDEFF : ', NF90_STRERROR(ierr)
  END SUBROUTINE CreateMergedFile

  SUBROUTINE FillMergedFile (sd_fout, knrank ) 
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE FillMergedFile  ***
    !!
    !! ** Purpose :  Fill the meged file with data from individual rank
    !!
    !! ** Method  :  Inner loop is on rank (tested to be the most efficient)
    !!
    !!----------------------------------------------------------------------
    TYPE(ncfile),    INTENT(inout) :: sd_fout   ! output file structure
    INTEGER(KIND=4), INTENT(in   ) :: knrank    ! Total number of rank 

    INTEGER(KIND=4) :: j3, j4, jvar, jrank      ! loop index
    INTEGER(KIND=4) :: ii1, ii2, ij1, ij2       ! local working integer
    INTEGER(KIND=4) :: ierr                     ! ncdf error status
    INTEGER(KIND=4) :: ilen1, ilen2             ! working integer
    INTEGER(KIND=4) :: inpi, inpj
    INTEGER(KIND=4), DIMENSION(2) :: ihalos, ihaloe

    REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: dl_wrk_1d ! use dble precision always, reduction
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dl_wrk_2d ! to ad-hoc precision performed at write time
    TYPE(ncfile), DIMENSION(:),   ALLOCATABLE :: sl_fin    ! local subdomain file structure

    INTEGER(KIND=4) :: ibloc_sz, ikr1, ikr2,iknrank, ijrank
    LOGICAL, DIMENSION(:), ALLOCATABLE :: ll_good
    CHARACTER(80) :: cl_format
    !!----------------------------------------------------------------------
    !      ierr = PrintNcFile(sd_fout)
    !      ierr = PrintNcFile(sf_in)
    ibloc_sz=knrank/nbloc
    ikr1=1
    ikr2=0
    DO WHILE (ikr2 /=  knrank) 
       ikr2=ikr1+ibloc_sz
       IF ( ikr2 > knrank ) ikr2=knrank 
       iknrank=ikr2-ikr1+1

       ALLOCATE ( sl_fin(iknrank) , ll_good(iknrank))
       ll_good(:) = .TRUE.
       IF (lg_verbose) PRINT '(i5.5,I5, 1x, a,1x, i5,1x,a)', mmpirank,nndone, ' OPEN THE ',iknrank,' SUB-DOMAIN FILE'
       DO jrank =ikr1,ikr2
          ijrank=jrank - ikr1 + 1
          WRITE(cl_format,'(a,i1,a,i1,a)') '(a,I',ndigit,'.',ndigit,',a)'
          WRITE(cf_in,cl_format) TRIM(cf_root)//'_',jrank-1,'.nc'
          sl_fin(ijrank) = GetNcFile ( cf_in )
          ! screen the procs
          IF ( lg_win) THEN
             IF ( sl_fin(ijrank)%iposition_last(1)  < mini )  ll_good(ijrank)=.FALSE.
             IF ( sl_fin(ijrank)%iposition_first(1) > maxi )  ll_good(ijrank)=.FALSE.
             IF ( sl_fin(ijrank)%iposition_last(2)  < minj )  ll_good(ijrank)=.FALSE.
             IF ( sl_fin(ijrank)%iposition_first(2) > maxj )  ll_good(ijrank)=.FALSE.
          ENDIF
          print *, TRIM(cf_in),"/",iknrank, ll_good(ijrank)


       ENDDO

       DO jvar = 1 , sd_fout%nvars 
          PRINT '(i5.5,i5,1x,a,i5,1x, a,a,i3)', mmpirank,nndone, 'merge variable ', jvar, TRIM(sd_fout%c_vnam(jvar)), &
               &  ' dimension : ', sd_fout%nvdim(jvar)
          SELECT CASE ( sd_fout%nvdim(jvar) )
          CASE ( 1 )  ! depth or time
             ! variables not depending on x,y , just copy from any spliited file,
             !  current in sf_in  so far
             ALLOCATE ( dl_wrk_1d( sf_in%nlen( sf_in%idimids(jvar,1) ) ) )
             dl_wrk_1d=0.d0
             ierr= NF90_GET_VAR(sf_in%ncid,   sf_in%nvid(jvar),   dl_wrk_1d )
             ierr= NF90_PUT_VAR(sd_fout%ncid, sd_fout%nvid(jvar), dl_wrk_1d )
             DEALLOCATE ( dl_wrk_1d )
          CASE ( 2 )   ! in general can be nav_lon, nav_lat or  variable (:, axis_nbounds) 
                       ! in domain_cfg file can also be a z,t file
             IF ( sd_fout%idimids(jvar,1) == sd_fout%idb .OR. sd_fout%idimids(jvar,2) == sd_fout%iunlim ) THEN ! bounds variable,  same for all ranks
                ilen1=sf_in%nlen(sd_fout%idimids(jvar,1))
                ilen2=sf_in%nlen(sd_fout%idimids(jvar,2))
                ALLOCATE (dl_wrk_2d (ilen1,  ilen2 ) )
                ierr= NF90_GET_VAR(sf_in%ncid,   sf_in%nvid(jvar)  ,   dl_wrk_2d )
                ierr= NF90_PUT_VAR(sd_fout%ncid, sd_fout%nvid(jvar),   dl_wrk_2d )
                DEALLOCATE (dl_wrk_2d )
             ELSE ! nav_lon nav_lat
                DO jrank = ikr1,ikr2
                   ijrank=jrank - ikr1 + 1
                   IF ( ll_good(ijrank) ) THEN
                      ii1=sl_fin(ijrank)%iposition_first(1)  ; ii2=ii1+sl_fin(ijrank)%npi-1
                      ij1=sl_fin(ijrank)%iposition_first(2)  ; ij2=ij1+sl_fin(ijrank)%npj-1
                      ihalos(:) = sl_fin(ijrank)%ihalo_size_start(:)
                      ihaloe(:) = sl_fin(ijrank)%ihalo_size_end(:)
                      inpi=sl_fin(ijrank)%npi
                      inpj=sl_fin(ijrank)%npj
                      ALLOCATE (dl_wrk_2d (inpi,inpj ) )
                      dl_wrk_2d=0.d0
                      ierr= NF90_GET_VAR(sl_fin(ijrank)%ncid,  sl_fin(ijrank)%nvid(jvar),   dl_wrk_2d )
                      IF ( lg_coord ) THEN
                         SELECT CASE ( sl_fin(ijrank)%c_vnam(jvar) )
                         CASE ('nav_lon') 
                            WHERE ( dl_wrk_2d == 0 .OR. dl_wrk_2d == -1. .OR. dl_wrk_2d > 200. )  dl_wrk_2d = dglam(ii1:ii2, ij1:ij2 )
                         CASE ('nav_lat') 
                            WHERE ( dl_wrk_2d == 0 .OR. dl_wrk_2d == -1. .OR. dl_wrk_2d > 200. )  dl_wrk_2d = dgphi(ii1:ii2, ij1:ij2 )
                         END SELECT
                      ENDIF

                      ierr= NF90_PUT_VAR( sd_fout%ncid,  sd_fout%nvid(jvar), dl_wrk_2d ,                   &
                           & start=(/ii1,ij1/), count=(/inpi,inpj/) )
                      DEALLOCATE ( dl_wrk_2d )
                   ENDIF
                END DO
             ENDIF
          CASE ( 3 )   ! likely to be x, y, t fields 
             ! loop on dim 3
             DO j3 = 1, sd_fout%nlen(sd_fout%idimids(jvar,3))
                DO jrank = ikr1,ikr2
                   ijrank=jrank - ikr1 + 1
                   IF ( ll_good(ijrank) ) THEN
                      ihalos(:) = sl_fin(ijrank)%ihalo_size_start(:)
                      ihaloe(:) = sl_fin(ijrank)%ihalo_size_end(:)
                      ii1=sl_fin(ijrank)%iposition_first(1)  ; ii2=ii1+sl_fin(ijrank)%npi-1
                      ij1=sl_fin(ijrank)%iposition_first(2)  ; ij2=ij1+sl_fin(ijrank)%npj-1
                      inpi=sl_fin(ijrank)%npi-ihaloe(1) - ihalos(1)
                      inpj=sl_fin(ijrank)%npj-ihaloe(2) - ihalos(2)
                      ALLOCATE (dl_wrk_2d (inpi,inpj ) )
                      dl_wrk_2d=0.d0
                      ierr = NF90_GET_VAR( sl_fin(ijrank)%ncid,  sl_fin(ijrank)%nvid(jvar),  dl_wrk_2d ,    &
                           &  start=(/ 1 + ihalos(:), j3 /), count=(/inpi,inpj,1/) )
                      ierr = NF90_PUT_VAR( sd_fout%ncid, sd_fout%nvid(jvar), dl_wrk_2d ,                    &
                           & start=(/ ii1,ij1, j3/),  count=(/inpi,inpj,1/) )
                      DEALLOCATE ( dl_wrk_2d )
                   ENDIF
                ENDDO
             END DO
          CASE ( 4 )  ! x,y,z,t fields
             ! loop on dim 4 and 3
             DO j4 = 1, sd_fout%nlen(sd_fout%idimids(jvar,4))
                DO j3 = 1, sd_fout%nlen(sd_fout%idimids(jvar,3))
                   DO jrank = ikr1,ikr2
                      ijrank=jrank - ikr1 + 1
                      IF ( ll_good(ijrank) ) THEN
                         ihalos(:) = sl_fin(ijrank)%ihalo_size_start(:)
                         ihaloe(:) = sl_fin(ijrank)%ihalo_size_end(:)
                         ii1=sl_fin(ijrank)%iposition_first(1)  ; ii2=ii1+sl_fin(ijrank)%npi-1
                         ij1=sl_fin(ijrank)%iposition_first(2)  ; ij2=ij1+sl_fin(ijrank)%npj-1
                         inpi=sl_fin(ijrank)%npi-ihaloe(1) - ihalos(1)
                         inpj=sl_fin(ijrank)%npj-ihaloe(2) - ihalos(2)
                         ALLOCATE (dl_wrk_2d (inpi,inpj ) )
                         dl_wrk_2d=0.d0
                         ierr = NF90_GET_VAR( sl_fin(ijrank)%ncid,  sl_fin(ijrank)%nvid(jvar),  dl_wrk_2d ,  &
                              &  start=(/ 1 + ihalos(:), j3 ,j4/), count=(/inpi,inpj,1/) )
                         ierr = NF90_PUT_VAR( sd_fout%ncid, sd_fout%nvid(jvar), dl_wrk_2d ,                  &
                              & start=(/ ii1,ij1, j3,j4/),  count=(/inpi,inpj,1,1/) )
                         DEALLOCATE ( dl_wrk_2d )
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO
          END SELECT
       ENDDO

       DO jrank = ikr1,ikr2
          ijrank=jrank - ikr1 + 1
          ierr = NF90_CLOSE(sl_fin(ijrank)%ncid )
       ENDDO
       DEALLOCATE ( sl_fin, ll_good ) 

       ikr1=ikr2+1
       ierr = NF90_SYNC(sd_fout%ncid)
    ENDDO
    ierr = NF90_CLOSE(sd_fout%ncid )

  END SUBROUTINE FillMergedFile

  FUNCTION GetCoord( )
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE GetCoord ***
    !!
    !! ** Purpose :  read dglam, dgphi from coordinate file, according to grid point
    !!
    !! ** Method  :  netcdf read size, check size, allocate, read ,close
    !!
    !!----------------------------------------------------------------------
    INTEGER(KIND=4)   :: GetCoord

    INTEGER(KIND=4)   :: ierr, icid, id, ix, iy
    CHARACTER(LEN=30) :: cl_glam
    CHARACTER(LEN=30) :: cl_gphi
    LOGICAL           :: lchk = .FALSE.
    !!----------------------------------------------------------------------
    cl_glam='glamt'  ! default longitude name (T point)
    cl_gphi='gphit'  ! defaut latitude name   (T point)
    IF ( INDEX(cf_root, 'gridU'  ) /= 0 ) THEN ; cl_glam='glamu' ; cl_gphi='gphiu' ; 
    ENDIF
    IF ( INDEX(cf_root, 'grid_U' ) /= 0 ) THEN ; cl_glam='glamu' ; cl_gphi='gphiu' ; 
    ENDIF
    IF ( INDEX(cf_root, 'gridV'  ) /= 0 ) THEN ; cl_glam='glamv' ; cl_gphi='gphiv' ; 
    ENDIF
    IF ( INDEX(cf_root, 'grid_V' ) /= 0 ) THEN ; cl_glam='glamv' ; cl_gphi='gphiv' ; 
    ENDIF

    ierr = NF90_OPEN(cf_coor, NF90_NOWRITE, icid )
    ierr = NF90_INQ_DIMID(icid,'x', id ) ; ierr = NF90_INQUIRE_DIMENSION( icid, id, len=ix)
    ierr = NF90_INQ_DIMID(icid,'y', id ) ; ierr = NF90_INQUIRE_DIMENSION( icid, id, len=iy)
    IF ( ix /= sf_in%isize_global(1) ) lchk = .TRUE.
    IF ( iy /= sf_in%isize_global(2) ) lchk = .TRUE.
    IF ( lchk ) THEN 
       PRINT *, 'DIMENSIONS OF ',TRIM(cf_coor), ' did not match global dims'
       ierr = NF90_CLOSE(icid )
       GetCoord = 1
       STOP ' error in coordinates dimensions'
    ENDIF

    ALLOCATE ( dglam(ix,iy) , dgphi(ix,iy) )
    ierr = NF90_INQ_VARID(icid,cl_glam, id ) ; ierr = NF90_GET_VAR(icid, id, dglam )
    ierr = NF90_INQ_VARID(icid,cl_gphi, id ) ; ierr = NF90_GET_VAR(icid, id, dgphi )

    ierr = NF90_CLOSE(icid )
    GetCoord = 0

  END FUNCTION GetCoord

  FUNCTION PrintNcFile (sd_nc)
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION PrintNcFile  ***
    !!
    !! ** Purpose :  Print in clear language, the content of the ncfile structure
    !!               passed as argument
    !!
    !! ** Method  :  Print structure to std output
    !!
    !!----------------------------------------------------------------------
    TYPE(ncfile), INTENT(in) :: sd_nc
    INTEGER(KIND=4)          :: PrintNcFile

    INTEGER(KIND=4)          :: jvar, jdim, jatt   ! loop index
    !!----------------------------------------------------------------------

    PRINT *, '** DATASET **  '
    PRINT *, 'NAME      : ', TRIM(sd_nc%c_fnam)
    PRINT *,' NCID      : ', sd_nc%ncid
    PRINT *
    PRINT *, '** DIMENSIONS **'
    PRINT *,' NDIMS     : ', sd_nc%ndims
    PRINT *,' UNLIMITED : ', sd_nc%iunlim
    PRINT *,' NPI       : ', sd_nc%npi
    PRINT *,' NPJ       : ', sd_nc%npj
    PRINT *,' NPK **    : ', sd_nc%npk
    PRINT *,' NPT       : ', sd_nc%nlen(sd_nc%idt)
    PRINT *,' IDX       : ', sd_nc%idx
    PRINT *,' IDY       : ', sd_nc%idy
    PRINT *,' IDZ **    : ', sd_nc%idz
    PRINT *,' IDT       : ', sd_nc%idt
    PRINT *,' NLEN      : ', (/ (sd_nc%nlen(jdim) , jdim = 1, sd_nc%ndims) /)
    PRINT *,' NAME      : ', (/ (TRIM(sd_nc%c_dnam(jdim))//' ' , jdim = 1, sd_nc%ndims) /)
    PRINT *
    PRINT *,' ** VARIABLES **'
    PRINT *,' NVARS     : ', sd_nc%nvars
    DO jvar = 1, sd_nc%nvars
       PRINT *, '     Name  : ', TRIM(sd_nc%c_vnam(jvar))
       PRINT *, '     Varid : ', sd_nc%nvid(jvar)
       PRINT *, '     Typ   : ', sd_nc%itype(jvar)
       PRINT *, '     Dims  : ', sd_nc%nvdim(jvar)
       PRINT *, '       DimId : ', (/ (sd_nc%idimids(jvar,jdim), jdim =1,sd_nc%nvdim(jvar) ) /)
       PRINT *, '     Attr  : ', sd_nc%nvatt(jvar)
    ENDDO
    PRINT * 
    PRINT *,' ** ATTRIBUTES **'
    PRINT *,' NATTS     : ', sd_nc%natts
    PRINT *
    PRINT *,' ** DOMAIN attribute **'
    PRINT *,' DOMAIN_number_total    : ', sd_nc%number_total
    PRINT *,' DOMAIN_number          : ', sd_nc%number
    PRINT *,' DOMAIN_dimensions_ids  : ', sd_nc%idimensions_ids(:)
    PRINT *,' DOMAIN_size_global     : ', sd_nc%isize_global(:)
    PRINT *,' DOMAIN_size_local      : ', sd_nc%isize_local(:)
    PRINT *,' DOMAIN_position_first  : ', sd_nc%iposition_first(:)
    PRINT *,' DOMAIN_position_last   : ', sd_nc%iposition_last(:)
    PRINT *,' DOMAIN_halo_size_start : ', sd_nc%ihalo_size_start(:)
    PRINT *,' DOMAIN_halo_size_end   : ', sd_nc%ihalo_size_end(:)

    PrintNcFile = 0

  END FUNCTION  PrintNcFile
END MODULE modutil
