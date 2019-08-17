PROGRAM bdy_coord_create
  !!======================================================================
  !!                     ***  PROGRAM  bdy_coord_create  ***
  !!=====================================================================
  !!  ** Purpose : Create bdy_coordinates file from the nambdy_index in 
  !!               case of straight boundaries. Also read nambdy to know 
  !!               the number of bdy as well as the rimwidth of each bdy.
  !!
  !!  ** Method  : Read the namelist and setup the nbi, npj, nbr for each
  !!               grid type : T U V. 
  !!
  !! History :  1.0  : 06/2019  : J.M. Molines : 
  !!----------------------------------------------------------------------
  USE netcdf
  !!----------------------------------------------------------------------
  IMPLICIT NONE
  
  INTEGER(KIND=4), PARAMETER :: jpbdy=20
  INTEGER(KIND=4) :: jb
  INTEGER(KIND=4) :: narg, iargc, ijarg

  INTEGER(KIND=4) :: numnam=10
  INTEGER(KIND=4) :: numnamd=11
  INTEGER(KIND=4) :: nb_bdy, n0=0, ir
  INTEGER(KIND=4) :: nbdyind,nbdybeg,nbdyend
  INTEGER(KIND=4), DIMENSION(jpbdy) :: nn_rimwidth

  CHARACTER(LEN=80) :: cf_namlist
  CHARACTER(LEN=80) :: cldum
  CHARACTER(LEN=80) :: ctypebdy
  CHARACTER(LEN=80), DIMENSION(jpbdy) :: cn_coords_file

  LOGICAL, DIMENSION(jpbdy)  :: ln_coords_file
  LOGICAL                    :: ll_rim0=.FALSE.

  NAMELIST/nambdy/ nb_bdy, ln_coords_file, cn_coords_file, nn_rimwidth
  NAMELIST/nambdy_index/ctypebdy,nbdyind,nbdybeg,nbdyend 
  !!----------------------------------------------------------------------
  narg = iargc()

  IF ( narg == 0 ) THEN
     PRINT *,' usage : bdy_coord_create -n namelist [-p] [-rim0]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Create bdy_coordinates file from nambdy_index information.' 
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -n namelist : namelist is a namelist file holding the NEMO nambdy_index'
     PRINT *,'            and a reduce version of nambdy with nb_bdy, ln_coords_file,'
     PRINT *,'            cn_coords_file and nn_rimwidth.' 
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       -rim0 : create a nbr=0 for bdy+second order numerical scheme'
     PRINT *,'       -p : print a template namelist namlist.tmpl to be edited to fit user'
     PRINT *,'            requirement.' 
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       none '
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file :  defined in cn_coords_file'
     PRINT *,'         variables : nbit, nbiu, nbiv,   nbjt, nbju, nbjv, nbrt, nbru, nbrv'
     PRINT *,'      '
     PRINT *,'      '
     STOP
  ENDIF

  ijarg = 1 
  DO WHILE ( ijarg <= narg )
     CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1
     SELECT CASE ( cldum )
     CASE ( '-n'   ) ; CALL getarg(ijarg, cf_namlist ) ; ijarg=ijarg+1
        ! option
     CASE ( '-p'   ) ; CALL PrintNamelist ; STOP
     CASE ( '-rim0') ; ll_rim0=.TRUE. ; n0=1
     CASE DEFAULT    ; PRINT *, ' ERROR : ', TRIM(cldum),' : unknown option.'; STOP 1
     END SELECT
  ENDDO

!  (0) dummy allocation :
!  ALLOCATE ( ln_coords_file(1), cn_coords_file(1), nn_rimwidth(1) )
!  (1) first read of namelist
   OPEN(numnam, FILE=cf_namlist)
   READ(numnam,nambdy) 
   PRINT *, nb_bdy
   DO jb=1, nb_bdy
     PRINT *, 'BDY ',jb,' : ', ln_coords_file(jb), TRIM(cn_coords_file(jb) ), nn_rimwidth(jb)
     IF ( ll_rim0 ) THEN
       PRINT *,'        starts at rim=0'
     ELSE
       PRINT *,'        starts at rim=1'
     ENDIF
   ENDDO
   REWIND(numnam)  ! only once 
   DO jb =1, nb_bdy
     READ(numnam,nambdy_index)
     IF ( ln_coords_file(jb)  ) THEN
        PRINT *, 'Working with bdy set ', jb
        PRINT *, '   Creating file :', TRIM(cn_coords_file(jb) )
        PRINT *, '   Rim Width     :', nn_rimwidth(jb)
        PRINT *, '    Type         :', TRIM(ctypebdy)
        SELECT CASE (ctypebdy)
        CASE ('N','S') ; PRINT *,' J-bdy : ', nbdyind ;
                       ; PRINT *,'   from I =', nbdybeg, ', to I =', nbdyend
        CASE ('E','W') ; PRINT *,' I-bdy : ', nbdyind ;
                       ; PRINT *,'   from J =', nbdybeg, ', to J =', nbdyend
        CASE DEFAULT   ; PRINT *,' type ', TRIM(ctypebdy),' unsupported' ; STOP 99
        END SELECT
        CALL CreateBdyCoords( ctypebdy, nbdyind,nbdybeg,nbdyend, cn_coords_file(jb),nn_rimwidth(jb) )
     ENDIF
   ENDDO

CONTAINS
  
  SUBROUTINE CreateBdyCoords( cd_typ, kbdyind, kbdybeg, kbdyend, cd_nam, krim )
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE CreateBdyCoords  ***
    !!
    !! ** Purpose :   Create BDY coordinate file for specific boundary
    !!
    !! ** Method  :    Use cd_type to determine the relative indices T U V
    !!
    !!----------------------------------------------------------------------
    CHARACTER(LEN=*), INTENT(in) :: cd_typ        ! type of BDY ( N S E W )
    INTEGER(KIND=4),  INTENT(in) :: kbdyind       ! Index of the bdy N S : J index
    !                                             !                  E W : I index
    INTEGER(KIND=4),  INTENT(in) :: kbdybeg       !  Along bdy limits (either I or J )
    INTEGER(KIND=4),  INTENT(in) :: kbdyend 
    CHARACTER(LEN=*), INTENT(in) :: cd_nam        ! Name of the output file
    INTEGER(KIND=4),  INTENT(in) :: krim          ! Rim width
    !
    INTEGER(KIND=4)                            :: jl, jr
    INTEGER(KIND=4)                            :: ii
    INTEGER(KIND=4)                            :: ixsiz, ilen
    INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: nbit, nbjt, nbrt
    INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: nbiu, nbju, nbru
    INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: nbiv, nbjv, nbrv
    INTEGER(KIND=4)                            :: ierr, icid
    INTEGER(KIND=4)                            :: idit, idiu, idiv
    INTEGER(KIND=4)                            :: idjt, idju, idjv
    INTEGER(KIND=4)                            :: idrt, idru, idrv
    INTEGER(KIND=4)                            :: idy, idxt, idxu, idxv
    !!----------------------------------------------------------------------
    ilen  = (kbdyend - kbdybeg + 1)
    ixsiz = ilen * (krim + n0 )

    ALLOCATE ( nbit(ixsiz), nbjt(ixsiz), nbrt(ixsiz) )
    ALLOCATE ( nbiu(ixsiz), nbju(ixsiz), nbru(ixsiz) )
    ALLOCATE ( nbiv(ixsiz), nbjv(ixsiz), nbrv(ixsiz) )

    ! All files have y dimension = 1 and xbt, xbu, xbv dimensions = ixsiz
    ierr = NF90_CREATE(cd_nam, NF90_NETCDF4, icid)
    ! define dims
    ierr = NF90_DEF_DIM(icid, 'yb',1, idy)
    ierr = NF90_DEF_DIM(icid, 'xbt',ixsiz, idxt)
    ierr = NF90_DEF_DIM(icid, 'xbu',ixsiz, idxu)
    ierr = NF90_DEF_DIM(icid, 'xbv',ixsiz, idxv)

    ! define variables
    ierr = NF90_DEF_VAR(icid,'nbit',NF90_INT,  (/idxt,idy/), idit )
    ierr = NF90_DEF_VAR(icid,'nbiu',NF90_INT,  (/idxu,idy/), idiu )
    ierr = NF90_DEF_VAR(icid,'nbiv',NF90_INT,  (/idxv,idy/), idiv )

    ierr = NF90_DEF_VAR(icid,'nbjt',NF90_INT,  (/idxt,idy/), idjt )
    ierr = NF90_DEF_VAR(icid,'nbju',NF90_INT,  (/idxu,idy/), idju )
    ierr = NF90_DEF_VAR(icid,'nbjv',NF90_INT,  (/idxv,idy/), idjv )

    ierr = NF90_DEF_VAR(icid,'nbrt',NF90_SHORT,(/idxt,idy/), idrt )
    ierr = NF90_DEF_VAR(icid,'nbru',NF90_SHORT,(/idxu,idy/), idru )
    ierr = NF90_DEF_VAR(icid,'nbrv',NF90_SHORT,(/idxv,idy/), idrv )

    ierr = NF90_ENDDEF(icid)
    
    SELECT CASE ( cd_typ )
    CASE ('S')
        ! Boundary is set at V , J point, T U at J
        ii=0
        IF ( ll_rim0 ) THEN
           ir=0
           DO jl=1,ilen
              ii= ii + 1
              nbjv(ii) = kbdyind + ir -1 
              nbjt(ii) = nbjv(ii)
              nbju(ii) = nbjv(ii)

              nbiv(ii) = kbdybeg + jl -1
              nbit(ii) = nbiv(ii)
              nbiu(ii) = nbiv(ii)

              nbrv(ii) = ir
              nbrt(ii) = ir
              nbru(ii) = ir
           ENDDO
        ENDIF
              
        DO jr=1,krim
           DO jl=1,ilen
              ii = ii + 1
              nbjv(ii) = kbdyind + jr -1 
              nbjt(ii) = nbjv(ii)
              nbju(ii) = nbjv(ii)

              nbiv(ii) = kbdybeg + jl -1
              nbit(ii) = nbiv(ii)
              nbiu(ii) = nbiv(ii)

              nbrv(ii) = jr
              nbrt(ii) = jr
              nbru(ii) = jr
           ENDDO
        ENDDO
    CASE ('W')
        ! Boundary is set at U , I point, T U at I
        ii=0
        IF ( ll_rim0 ) THEN
           ir=0
           DO jl=1,ilen
              ii= ii + 1
              nbiu(ii) = kbdyind + ir -1
              nbit(ii) = nbiu(ii)
              nbiv(ii) = nbiu(ii)

              nbju(ii) = kbdybeg + jl -1
              nbjt(ii) = nbju(ii)
              nbjv(ii) = nbju(ii)

              nbrv(ii) = ir
              nbrt(ii) = ir
              nbru(ii) = ir
           ENDDO
        ENDIF
        DO jr=1,krim
           DO jl=1,ilen
              ii = ii + 1
              nbiu(ii) = kbdyind + jr -1
              nbit(ii) = nbiu(ii)
              nbiv(ii) = nbiu(ii)

              nbju(ii) = kbdybeg + jl -1
              nbjt(ii) = nbju(ii)
              nbjv(ii) = nbju(ii)

              nbrv(ii) = jr
              nbrt(ii) = jr
              nbru(ii) = jr
           ENDDO
        ENDDO

    CASE ('N')
        ii=0
        IF ( ll_rim0 ) THEN
           ir=0
           DO jl=1,ilen
              ii= ii + 1
              nbjv(ii) = kbdyind - ir + 1
              nbjt(ii) = nbjv(ii) + 1
              nbju(ii) = nbjv(ii) + 1

              nbiv(ii) = kbdybeg + jl - 1
              nbit(ii) = nbiv(ii)
              nbiu(ii) = nbiv(ii)

              nbrv(ii) = ir
              nbrt(ii) = ir
              nbru(ii) = ir
           ENDDO
        ENDIF
        DO jr=1,krim
           DO jl=1,ilen
              ii = ii + 1
              nbjv(ii) = kbdyind - jr + 1 
              nbjt(ii) = nbjv(ii) + 1 
              nbju(ii) = nbjv(ii) + 1
 
              nbiv(ii) = kbdybeg + jl - 1
              nbit(ii) = nbiv(ii)
              nbiu(ii) = nbiv(ii)

              nbrv(ii) = jr
              nbrt(ii) = jr
              nbru(ii) = jr
           ENDDO
        ENDDO
    CASE ('E')
        ii=0
        IF ( ll_rim0 ) THEN
           ir=0
           DO jl=1,ilen
              ii= ii + 1
              nbiu(ii) = kbdyind - ir + 1
              nbit(ii) = nbiu(ii) + 1
              nbiv(ii) = nbiu(ii) + 1

              nbju(ii) = kbdybeg + jl - 1
              nbjt(ii) = nbju(ii)
              nbjv(ii) = nbju(ii)

              nbrv(ii) = ir
              nbrt(ii) = ir
              nbru(ii) = ir
           ENDDO
        ENDIF
        DO jr=1,krim
           DO jl=1,ilen
              ii = ii + 1
              nbiu(ii) = kbdyind - jr + 1
              nbit(ii) = nbiu(ii) + 1
              nbiv(ii) = nbiu(ii) + 1

              nbju(ii) = kbdybeg + jl - 1
              nbjt(ii) = nbju(ii)
              nbjv(ii) = nbju(ii)

              nbrv(ii) = jr
              nbrt(ii) = jr
              nbru(ii) = jr
           ENDDO
        ENDDO
    END SELECT

    ierr = NF90_PUT_VAR(icid, idit, nbit)
    ierr = NF90_PUT_VAR(icid, idiu, nbiu)
    ierr = NF90_PUT_VAR(icid, idiv, nbiv)

    ierr = NF90_PUT_VAR(icid, idjt, nbjt)
    ierr = NF90_PUT_VAR(icid, idju, nbju)
    ierr = NF90_PUT_VAR(icid, idjv, nbjv)

    ierr = NF90_PUT_VAR(icid, idrt, nbrt)
    ierr = NF90_PUT_VAR(icid, idru, nbru)
    ierr = NF90_PUT_VAR(icid, idrv, nbrv)
    
    ierr = NF90_CLOSE(icid) 

    DEALLOCATE ( nbit, nbjt, nbrt )
    DEALLOCATE ( nbiu, nbju, nbru )
    DEALLOCATE ( nbiv, nbjv, nbrv )

  END SUBROUTINE CreateBdyCoords
    
  SUBROUTINE PrintNamelist()
     INTEGER(KIND=4), PARAMETER :: jpline=25
     INTEGER(KIND=4) :: jl
     CHARACTER(LEN=256),DIMENSION(jpline) :: cline
!-------------------------------------------------
     OPEN(numnamd,FILE='namelist.tmpl')

   
    cline(:)=(/'!-----------------  tmpl namelist for bdy_coord_create ------------------------------------',     &
             & '&nambdy ',                                                                                        &
             & '   nb_bdy         = 3         !  number of open boundary sets   ',                                &
             & '   ln_coords_file = .true.,.true.,.true.   !  =T : read bdy coordinates from file ',              &
             & '   cn_coords_file = ''sbdy_coordinate_file.nc'',''nbdy_coordinate_file.nc'',''wbdy_coordinate_file.nc'' ' , &
             & '   nn_rimwidth  = 9,9,9           !  width of the relaxation zone',                               &
             & '/', &
             & '&nambdy_index',                                                                                   &
             & '    ctypebdy = ''S''',                                                                            &
             & '    nbdyind  =  3',                                                                               &
             & '    nbdybeg  = 2',                                                                                &
             & '    nbdyend  = 933',                                                                              &
             & '/', &
             & '&nambdy_index',                                                                                   &
             & '    ctypebdy = ''N'' ',                                                                            &
             & '    nbdyind  =  399',                                                                               &
             & '    nbdybeg  = 2',                                                                                &
             & '    nbdyend  = 933',                                                                              &
             & '/', &
             & '&nambdy_index',                                                                                   &
             & '    ctypebdy = ''W''' ,                                                                          &
             & '    nbdyind  =  76',                                                                               &
             & '    nbdybeg  = 301',                                                                                &
             & '    nbdyend  = 358',                                                                              &
             & '/'/)

   DO jl=1, jpline
      WRITE(numnamd,'(a)') TRIM(cline(jl) )
   ENDDO
    CLOSE (numnamd)

  END SUBROUTINE PrintNamelist

END PROGRAM bdy_coord_create
