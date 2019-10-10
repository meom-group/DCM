PROGRAM bdy_mk_coordinates_mpp
  !!======================================================================
  !!                     ***  PROGRAM  bdy_mk_coordinates  ***
  !!=====================================================================
  !!  ** Purpose : Read a rim file (build with BMGTOOLS) where rims are
  !!               given the values 100 200 etc for rim 0 1 ... and produce
  !!               a BDY coordinates files ( nbit, nbjt, nbrt)
  !!
  !!  ** Method  : scan the input file in a clever way ...
  !!
  !! History :  1.0  : 08/2019  : J.M. Molines : 
  !!----------------------------------------------------------------------
  !!----------------------------------------------------------------------
  !!   routines      : description
  !!----------------------------------------------------------------------
  USE netcdf
  USE mpi
  !!----------------------------------------------------------------------
  !! $Id$
  !! Copyright (c) 2019, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !!----------------------------------------------------------------------
  IMPLICIT NONE
  INTEGER(KIND=4)               :: ji,jj,jr
  INTEGER(KIND=4)               :: narg, iargc, ijarg
  INTEGER(KIND=4)               :: npiglo, npjglo
  INTEGER(KIND=4)               :: iioff, ijoff
  INTEGER(KIND=4)               :: nrim,  irim
  INTEGER(KIND=4)               :: nxt, ii, nvar
  INTEGER(KIND=4)               :: ncid, id, ierr
  INTEGER(KIND=4)               :: idxt,idxu,idxv, idy
  INTEGER(KIND=4)               :: idit,idjt,idrt
  INTEGER(KIND=4)               :: idiu,idju,idru
  INTEGER(KIND=4)               :: idiv,idjv,idrv
  INTEGER(KIND=4)               :: mpisize, mpirank,mpierr

  INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: nbit, nbjt, nbrt
  INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: nbiu, nbju, nbru
  INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: nbiv, nbjv, nbrv
  INTEGER(KIND=4), DIMENSION(2)              :: idims

  REAL(KIND=4), DIMENSION(:,:),  ALLOCATABLE :: v2d,zwk

  CHARACTER(LEN=255) :: cf_rim 
  CHARACTER(LEN=255) :: cv_rim ="Bathymetry"  ! BMGTOOLS does it !
  CHARACTER(LEN=255) :: cf_coo 
  CHARACTER(LEN=255) :: c_dimx ="x"         
  CHARACTER(LEN=255) :: c_dimy ="y"        
  CHARACTER(LEN=255) :: c_dimt ="time_counter"           
  CHARACTER(LEN=255) :: c_dimz ="z"          
  CHARACTER(LEN=80)  :: cldum
  CHARACTER(LEN=80)  :: cbdy_typ='none'
  CHARACTER(LEN=255) :: cf_out='bdy_coordinates.nc'

  CHARACTER(LEN=20) :: cn_votemper='votemper', cn_vosaline='vosaline', cn_sossheig='sossheig'
  CHARACTER(LEN=20) :: cn_vozocrtx='vozocrtx', cn_vomecrty='vomecrty'
  CHARACTER(LEN=20) :: cn_ileadfra='ileadfra', cn_iicethic='iicethic', cn_isnowthi='isnowthi'

  LOGICAL            :: ll_data=.FALSE.
  LOGICAL            :: ll_rim =.FALSE.
  LOGICAL            :: ll_coor=.FALSE.
  LOGICAL            :: ll_votemper=.TRUE., ll_vosaline=.TRUE., ll_sossheig=.TRUE.
  LOGICAL            :: ll_vozocrtx=.TRUE., ll_vomecrty=.TRUE.
  LOGICAL            :: ll_ileadfra=.TRUE., ll_iicethic=.TRUE., ll_isnowthi=.TRUE.

  !!---------------------------------------------------------------------
  cf_coo = "None"
  narg=iargc()
  iioff=1
  ijoff=1
  IF ( narg == 0 ) THEN
     PRINT *,'   '
     PRINT *,'  usage : bdy_mk_coordinates_from_file_mpp -r RIM-file | -c BDY-coord  -t BDY-type'
     PRINT *,'          [-o BDY-coordinates-file] [-offset I-offset J-offset]'
     PRINT *,'          [-data] [-v LST-vars]'
     PRINT *,'   '
     PRINT *,'  PURPOSE :'
     PRINT *,'     Build bdy coordinates file from rim file. At this stage this program'
     PRINT *,'     assumes that boundary are recti-linear (however, oblique). Rim file is'
     PRINT *,'     a hand-prepared files (BMGTOOLS) using tmaskutil as a basis and where'
     PRINT *,'     rim number is indicated for T points, starting with 100 (rim-0) to '
     PRINT *,'     100*(rimwidth) + 100.' 
     PRINT *,'   '
     PRINT *,'  ARGUMENTS :'
     PRINT *,'      Use only one of -r or -c options.'
     PRINT *,'    -r RIM-file : gives the name of the original rim file'
     PRINT *,'    -c BDY-COORD : gives the name of a bdy_coordinate file holding the nbi,nbj'
     PRINT *,'    -t BDY-type : One of N S E W. Although the boundary geometry is not'
     PRINT *,'        strictly N S E or W, as far as we are dealing with strait lines'
     PRINT *,'        we can always figure out one of those type for the BDY.'
     PRINT *,'   '
     PRINT *,'  OPTIONS :'
     PRINT *,'    -o BDY-coordinates-file : name of output file instead of default '
     PRINT *,'        name : ',TRIM(cf_out)
     PRINT *,'    -offset I-offset J-offset : If specified, the resulting nbix, nbjx '
     PRINT *,'        from this program will be offset with the given (I,J) pait. This is '
     PRINT *,'        usefull for building the final configuration bdy coordinates, as in'
     PRINT *,'        general work is done locally on a smaller domain around the bdy.'
     PRINT *,'        Iglobal = Ilocal + I-offset -1'
     PRINT *,'        Jglobal = Jlocal + J-offset -1'
     PRINT *,'    -data  : This option forces the program to build bdy_data files from'
     PRINT *,'        pre processed input files for all required variables'
     PRINT *,'    -v VAR-lst : give a comma separated list of the variables you want to '
     PRINT *,'        process. Default is to process all variables.'
     PRINT *,'        Possible variables are votemper,vosaline,sossheig,vozocrtx,vomecrty,'
     PRINT *,'        iicethic,isnowthi and ileadfra.'
     PRINT *,'   '
     PRINT *,'  REQUIRED FILES :'
     PRINT *,'    none'
     PRINT *,'   '
     PRINT *,'  OUTPUT : '
     PRINT *,'    netcdf file : ',TRIM(cf_out)
     PRINT *,'      variables :  nbit,nbjt,nbr and equivalent for u and '
     PRINT *,'   '
     PRINT *,'  SEE ALSO :'
     PRINT *,'     bdy_mk_coordinate'
     PRINT *,'      '
     STOP
  ENDIF

  ijarg = 1 
  DO WHILE ( ijarg <= narg )
     CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1
     SELECT CASE ( cldum )
     CASE ( '-r'     ) ; CALL getarg(ijarg, cf_rim ) ; ijarg=ijarg+1 ; ll_rim =.TRUE.
     CASE ( '-c'     ) ; CALL getarg(ijarg, cf_coo ) ; ijarg=ijarg+1 ; ll_coor=.TRUE.
     CASE ( '-t'     ) ; CALL getarg(ijarg,cbdy_typ) ; ijarg=ijarg+1
        ! options
     CASE ( '-o'     ) ; CALL getarg(ijarg,cf_out  ) ; ijarg=ijarg+1
     CASE ( '-offset') ; CALL getarg(ijarg,cldum   ) ; ijarg=ijarg+1 ; READ(cldum,*) iioff
        ;              ; CALL getarg(ijarg,cldum   ) ; ijarg=ijarg+1 ; READ(cldum,*) ijoff
     CASE ( '-data  ') ; ll_data=.TRUE.
     CASE ( '-v'     ) ; CALL getarg(ijarg, cldum )  ; ijarg=ijarg+1
        ;              ; CALL ParseVars ( cldum   )
     CASE DEFAULT      ; PRINT *, ' ERROR : ', TRIM(cldum),' : unknown option.'; STOP 1

     END SELECT
  ENDDO
  IF ( cbdy_typ == 'none' ) THEN
     PRINT *,'  +++ E R R O R: You must specify a bdy-type with option -t '
     STOP 999
  ENDIF
  IF ( ll_rim .AND. ll_coor ) THEN
     PRINT *,'  +++ E R R O R: You must choose onmy one of -r or -c option'
     STOP 999
  ENDIF
  

  CALL MPI_INIT(mpierr)
  IF ( mpierr /= 0 ) PRINT *, ' ERROR in mpi_init' 
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, mpisize, mpierr)
  IF ( mpierr /= 0 ) PRINT *, ' ERROR in mpi_comm_size' 
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, mpirank, mpierr)
  IF ( mpierr /= 0 ) PRINT *, ' ERROR in mpi_comm_rank' 

PRINT *,' RANK : ', mpirank,mpisize

  IF ( ll_rim ) THEN
  ! read rim file
  ierr = NF90_OPEN(cf_rim,NF90_NOWRITE,ncid)
  ierr = NF90_INQ_DIMID(ncid,c_dimx,id)  ; ierr = NF90_INQUIRE_DIMENSION(ncid,id,len=npiglo)
  ierr = NF90_INQ_DIMID(ncid,c_dimy,id)  ; ierr = NF90_INQUIRE_DIMENSION(ncid,id,len=npjglo)

  ALLOCATE (v2d(npiglo,npjglo), zwk(npiglo,npjglo) )
  ierr = NF90_INQ_VARID(ncid,cv_rim,id)  ! BMGTOOLS impose Bathymetry ...
  ierr = NF90_GET_VAR(ncid,id,v2d )
  ierr = NF90_CLOSE(ncid)

  zwk=0.
  WHERE( v2d >= 100 ) zwk = 1.
  nrim=( MAXVAL(v2d) - 100 )/100 + 1

  nxt = SUM(zwk)
  ALLOCATE (nbit(nxt), nbjt(nxt), nbrt(nxt) )
  ALLOCATE (nbiu(nxt), nbju(nxt), nbru(nxt) )
  ALLOCATE (nbiv(nxt), nbjv(nxt), nbrv(nxt) )
  nbit=-9999
  nbjt=-9999
  nbrt=-9999


  PRINT *, ' Number of rims : ', nrim
  PRINT *, '           nXt  : ', nxt

  ii=0
  DO jr=1,nrim
     DO ji=1,npiglo
        DO jj=1,npjglo
           IF ( v2d(ji,jj) < 100 ) CYCLE
           irim = (v2d(ji,jj) - 100 )/100
           IF ( irim == jr-1  ) THEN
              ii = ii + 1
              nbit(ii) = ji
              nbjt(ii) = jj
              nbrt(ii) = irim
           ENDIF
        ENDDO
     ENDDO
  ENDDO

  ! now deduce the u, v values (
  SELECT CASE ( cbdy_typ )
  CASE ('N')
     nbiu=nbit
     nbju=nbjt
     nbru=nbrt
     nbiv=nbit
     nbjv=nbjt-1
     nbrv=nbrt
  CASE ('S','W')
     nbiu=nbit
     nbju=nbjt
     nbru=nbrt
     nbiv=nbit
     nbjv=nbjt
     nbrv=nbrt
  CASE ('E')
     nbiu=nbit-1
     nbju=nbjt
     nbru=nbrt
     nbiv=nbit
     nbjv=nbjt
     nbrv=nbrt
  END SELECT

  DEALLOCATE( v2d, zwk)
  ENDIF

  IF ( ll_coor) THEN
     ! read data in bdy_coord
     ierr = NF90_OPEN(cf_coo,NF90_NOWRITE,ncid)
     ! find out the name of the dimensions
     ierr = NF90_INQ_VARID(ncid,'nbit',id)
     ierr = NF90_INQUIRE_VARIABLE(ncid, id, dimids=idims)
     ierr = NF90_INQUIRE_DIMENSION(ncid,idims(1),len=nxt) 

     ALLOCATE (nbit(nxt), nbjt(nxt) )
     ALLOCATE (nbiu(nxt), nbju(nxt) )
     ALLOCATE (nbiv(nxt), nbjv(nxt) )
     ierr = NF90_INQ_VARID(ncid,'nbit',id) ; ierr=NF90_GET_VAR(ncid,id,nbit,start=(/1,1/), count=(/nxt,1/) )
     ierr = NF90_INQ_VARID(ncid,'nbjt',id) ; ierr=NF90_GET_VAR(ncid,id,nbjt,start=(/1,1/), count=(/nxt,1/) )
     ierr = NF90_INQ_VARID(ncid,'nbiu',id) ; ierr=NF90_GET_VAR(ncid,id,nbiu,start=(/1,1/), count=(/nxt,1/) )
     ierr = NF90_INQ_VARID(ncid,'nbju',id) ; ierr=NF90_GET_VAR(ncid,id,nbju,start=(/1,1/), count=(/nxt,1/) )
     ierr = NF90_INQ_VARID(ncid,'nbiv',id) ; ierr=NF90_GET_VAR(ncid,id,nbiv,start=(/1,1/), count=(/nxt,1/) )
     ierr = NF90_INQ_VARID(ncid,'nbjv',id) ; ierr=NF90_GET_VAR(ncid,id,nbjv,start=(/1,1/), count=(/nxt,1/) )
     ierr=NF90_CLOSE(ncid)
     ! reverse offset

     nbit=nbit- iioff+1 ; nbiu=nbiu-iioff+1 ;  nbiv=nbiv-iioff+1
     nbjt=nbjt- ijoff+1 ; nbju=nbju-ijoff+1 ;  nbjv=nbjv-ijoff+1
  ENDIF


  IF (ll_data) THEN
     IF ( ll_votemper ) CALL CreateData(cn_votemper,'T')
     IF ( ll_vosaline ) CALL CreateData(cn_vosaline,'T')
     IF ( ll_sossheig ) CALL CreateData(cn_sossheig,'T')
     IF ( ll_vozocrtx ) CALL CreateData(cn_vozocrtx,'U')
     IF ( ll_vomecrty ) CALL CreateData(cn_vomecrty,'V')
     IF ( ll_ileadfra ) CALL CreateData(cn_ileadfra,'T')
     IF ( ll_iicethic ) CALL CreateData(cn_iicethic,'T')
     IF ( ll_isnowthi ) CALL CreateData(cn_isnowthi,'T')
  ENDIF

  IF ( mpirank == 0 .AND. ll_rim ) THEN
  ! write output file (bdy_coordinates file)
  ierr = NF90_CREATE( cf_out, NF90_NETCDF4,ncid)
  ierr = NF90_DEF_DIM(ncid,'yb',1, idy)
  ierr = NF90_DEF_DIM(ncid,'xbT',nxt, idxt)
  ierr = NF90_DEF_DIM(ncid,'xbU',nxt, idxu)
  ierr = NF90_DEF_DIM(ncid,'xbV',nxt, idxv)

  ierr = NF90_DEF_VAR(ncid,'nbit',NF90_INT,(/idxt,idy/),idit)
  ierr = NF90_DEF_VAR(ncid,'nbjt',NF90_INT,(/idxt,idy/),idjt)
  ierr = NF90_DEF_VAR(ncid,'nbrt',NF90_INT,(/idxt,idy/),idrt)

  ierr = NF90_DEF_VAR(ncid,'nbiu',NF90_INT,(/idxu,idy/),idiu)
  ierr = NF90_DEF_VAR(ncid,'nbju',NF90_INT,(/idxu,idy/),idju)
  ierr = NF90_DEF_VAR(ncid,'nbru',NF90_INT,(/idxu,idy/),idru)

  ierr = NF90_DEF_VAR(ncid,'nbiv',NF90_INT,(/idxu,idy/),idiv)
  ierr = NF90_DEF_VAR(ncid,'nbjv',NF90_INT,(/idxu,idy/),idjv)
  ierr = NF90_DEF_VAR(ncid,'nbrv',NF90_INT,(/idxu,idy/),idrv)

  ierr = NF90_PUT_ATT(ncid,NF90_GLOBAL,'rimwidth',nrim-1)

  ierr = NF90_ENDDEF(ncid)
  ! apply offset
  nbit=nbit+ iioff-1 ; nbiu=nbiu+iioff-1 ;  nbiv=nbiv+iioff-1
  nbjt=nbjt+ ijoff-1 ; nbju=nbju+ijoff-1 ;  nbjv=nbjv+ijoff-1

  ierr = NF90_PUT_VAR( ncid,idit,nbit)
  ierr = NF90_PUT_VAR( ncid,idjt,nbjt)
  ierr = NF90_PUT_VAR( ncid,idrt,nbrt)

  ierr = NF90_PUT_VAR( ncid,idiu,nbiu)
  ierr = NF90_PUT_VAR( ncid,idju,nbju)
  ierr = NF90_PUT_VAR( ncid,idru,nbru)

  ierr = NF90_PUT_VAR( ncid,idiv,nbiv)
  ierr = NF90_PUT_VAR( ncid,idjv,nbjv)
  ierr = NF90_PUT_VAR( ncid,idrv,nbrv)

  ierr = NF90_CLOSE(ncid)
  ENDIF
  CALL MPI_FINALIZE(mpierr)

CONTAINS

  SUBROUTINE CreateData ( cd_var, cd_typ)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE CreateData  ***
    !!
    !! ** Purpose :  Create the bdydata netcdffile corresponding to cd_var
    !!               according to its position on the C-grid 
    !!
    !! ** Method  :  Read or use pre computed nbix, nbjx bdy points position
    !!               and extract the relevant variables from the data files,
    !!               assuming they already exists on a bubset of the model grid. 
    !!
    !!----------------------------------------------------------------------
    ! Arguments
    CHARACTER(LEN=*), INTENT(in) :: cd_var
    CHARACTER(LEN=1), INTENT(in) :: cd_typ
    ! local
    INTEGER(KIND=4) :: jf, jt, jx, jk  ! dummy loop index
    INTEGER(KIND=4) :: infile,  ierr, incid, id, idim
    INTEGER(KIND=4) :: incout, idout
    INTEGER(KIND=4) :: ipiglo, ipjglo, ipt, ipkglo, ii, ij
    INTEGER(KIND=4) :: iblock ! loop over processes 
    INTEGER(KIND=4), DIMENSION(nxt) :: nbi, nbj

    REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: zdata
    REAL(KIND=4), DIMENSION(:)    , ALLOCATABLE :: zbdydata
    REAL(KIND=4)                                :: zspval

    CHARACTER(LEN=255), DIMENSION(:), ALLOCATABLE :: cl_list
    CHARACTER(LEN=255)                            :: clf_in
    
    LOGICAL :: l_go
    !!----------------------------------------------------------------------
    CALL GetList( cd_var, cl_list, infile)    
    iblock = infile/mpisize
    PRINT *, ' ... Processing ',iblock, ' of ',mpisize, ' files from ', infile
    DO jf= 1, iblock+1
       clf_in=TRIM(cl_list((mpirank+1)+mpisize*(jf-1)))
       IF ( (mpirank+1)+mpisize*(jf-1) > infile ) EXIT
       PRINT *,mpirank,  clf_in
       ierr=NF90_OPEN(clf_in,NF90_NOWRITE,incid)
       ierr=NF90_INQ_VARID(incid,cd_var,id)
       ierr=NF90_INQUIRE_VARIABLE(incid,id,ndims=idim)

       ierr=NF90_INQ_DIMID(incid,c_dimx,id) ; ierr=NF90_INQUIRE_DIMENSION(incid,id,len=ipiglo)
       ierr=NF90_INQ_DIMID(incid,c_dimy,id) ; ierr=NF90_INQUIRE_DIMENSION(incid,id,len=ipjglo)
       ierr=NF90_INQ_DIMID(incid,c_dimt,id) ; ierr=NF90_INQUIRE_DIMENSION(incid,id,len=ipt)
       ipkglo=1
       IF ( idim ==4 )  THEN
          ierr=NF90_INQ_DIMID(incid,c_dimz,id) ; ierr=NF90_INQUIRE_DIMENSION(incid,id,len=ipkglo)
       ENDIF
       ALLOCATE(zbdydata(nxt) )
       IF (ll_rim ) THEN
       IF ( ipiglo /= npiglo .OR. ipjglo /= npjglo) THEN
          PRINT *, ' inconsistent domain size between rim file and data file '
          STOP 'ERROR 1'
       ENDIF
       ELSE 
         npiglo=ipiglo
         npjglo=ipjglo
       ENDIF
       ALLOCATE (zdata(npiglo,npjglo) )
       ierr = NF90_INQ_VARID(incid,cd_var,id)
       ierr = NF90_GET_ATT(incid,id,'_FillValue',zspval )
       IF ( ierr /= NF90_NOERR)      PRINT *,' INQ_GET_ATT :', TRIM(NF90_STRERROR(ierr))
       CALL CreateOutput(clf_in,cd_var,ipkglo, zspval, incout, idout)
       DO jt=1,ipt
          PRINT *,TRIM(cd_var), ' JT : ', jt,' / ', ipt
          l_go =.TRUE.
          DO jk=1, ipkglo
             zbdydata(:) = zspval
!           PRINT *,'JK : ', jk,' / ',ipkglo
             IF ( l_go ) THEN      
             IF ( ipkglo == 1 ) THEN
               ierr = NF90_GET_VAR(incid,id,zdata,start=(/1,1,jt/), count=(/npiglo,npjglo,1/) )
             ELSE
               ierr = NF90_GET_VAR(incid,id,zdata,start=(/1,1,jk,jt/), count=(/npiglo,npjglo,1,1/) )
             ENDIF
             IF ( ierr /= NF90_NOERR)      PRINT *,' GET_VAR:', TRIM(NF90_STRERROR(ierr))
             SELECT CASE ( cd_typ)
             CASE ('T') ; nbi=nbit ; nbj=nbjt
             CASE ('U') ; nbi=nbiu ; nbj=nbju
             CASE ('V') ; nbi=nbiv ; nbj=nbjv
             END SELECT
             DO jx=1,nxt
                ii=nbi(jx)
                ij=nbj(jx)
                zbdydata(jx) = zdata(ii,ij)
             ENDDO
             IF ( count ( (zbdydata == zspval)) == nxt ) l_go=.FALSE.
             ENDIF
             IF (ipkglo == 1 ) THEN
               ierr = NF90_PUT_VAR(incout,idout,zbdydata(:), start=(/1,1,jt/), count=(/nxt,1,1/) )
             ELSE
               ierr = NF90_PUT_VAR(incout,idout,zbdydata(:), start=(/1,1,jk,jt/), count=(/nxt,1,1,1/) )
             ENDIF
          ENDDO
       ENDDO
       ierr = NF90_CLOSE(incout)

       DEALLOCATE ( zdata, zbdydata )
       ierr = NF90_CLOSE(incid)
    ENDDO
  END SUBROUTINE CreateData

  SUBROUTINE CreateOutput( cd_file, cd_var, kpkglo, pspval,kncout, kidout)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE CreateOutput  ***
    !!
    !! ** Purpose :   Create bdy_data netcdf output file
    !!
    !! ** Method  :   Infer name from data file, create variable, according to
    !!                its 2D/3D character and return ncout and idout
    !!
    !!----------------------------------------------------------------------
    ! Arguments
    CHARACTER(LEN=*), INTENT(in) :: cd_file  ! name of the data files used for xtrac
    CHARACTER(LEN=*), INTENT(in) :: cd_var   ! name of the  variable to process

    INTEGER(KIND=4), INTENT(in )  :: kpkglo  ! number of dimensions for cd_var
    REAL(KIND=4),    INTENT(in )  :: pspval  ! Fill value
    INTEGER(KIND=4), INTENT(out)  :: kncout  ! netcdf id of the output file
    INTEGER(KIND=4), INTENT(out)  :: kidout  ! varid of the output variables
    ! Local 
    INTEGER(KIND=4)    :: ierr
    INTEGER(KIND=4)    :: idx, idy,idz,idt
    CHARACTER(LEN=255) :: cl_file
    !!----------------------------------------------------------------------
    cl_file='bdydta_'//TRIM(cd_file)
    ierr = NF90_CREATE(cl_file,NF90_NETCDF4,kncout)

    ierr = NF90_DEF_DIM(kncout,'yb',1  ,idy)
    ierr = NF90_DEF_DIM(kncout,'xt',nxt,idx)
    ierr = NF90_DEF_DIM(kncout,'time_counter',NF90_UNLIMITED,idt)
    IF ( kpkglo /= 1  ) THEN
       ierr = NF90_DEF_DIM(kncout,'z',kpkglo,idz)
       ierr = NF90_DEF_VAR(kncout,cd_var,NF90_FLOAT,(/idx,idy,idz,idt/), kidout )
    ELSE
       ierr = NF90_DEF_VAR(kncout,cd_var,NF90_FLOAT,(/idx,idy,idt/), kidout )
    ENDIF
       ierr = NF90_PUT_ATT(kncout,kidout,'_FillValue',pspval) 
    ierr = NF90_ENDDEF(kncout)

  END SUBROUTINE CreateOutput


  SUBROUTINE GetList ( cd_var, cd_list, kfile )
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE GetList  ***
    !!
    !! ** Purpose :  Allocate and fill in cd_list(:) with the names of all
    !!               files whose name starts with cd_var. Return also the 
    !!               number of files.
    !!
    !! ** Method  :  Use a system call that output to a hidden file and then
    !!               read the file to establish the list.
    !!
    !!----------------------------------------------------------------------
    ! Arguments
    CHARACTER(LEN=*),                               INTENT(in ) :: cd_var
    CHARACTER(LEN=255), DIMENSION(:), ALLOCATABLE , INTENT(out) :: cd_list
    INTEGER(KIND=4),                                INTENT(out) :: kfile
    ! Local
    INTEGER(KIND=4) :: jf
    INTEGER(KIND=4) :: inumlst=10
    CHARACTER(LEN=4) :: crank
    !!----------------------------------------------------------------------
    WRITE(crank,'(I4.4)') mpirank
    CALL SYSTEM( "find . -maxdepth 1 -name """ //  TRIM(cd_var) // "*"" -printf ""%f\n"" | sort -r  > .zbdyfile_list."//TRIM(crank))
    kfile=0
    OPEN(inumlst, FILE='.zbdyfile_list.'//TRIM(crank) )
    DO
       READ(inumlst,*,END=999)
       kfile=kfile+1
    ENDDO
999 PRINT *, 'Files to process : ', kfile
    ALLOCATE( cd_list(kfile) )
    REWIND(inumlst)
    DO jf=1,kfile
       READ(inumlst,"(a)")  cd_list(jf)
    ENDDO
    CLOSE(inumlst)

    CALL SYSTEM( "rm .zbdyfile_list."//TRIM(crank))

  END SUBROUTINE GetList
   SUBROUTINE ParseVars (cdum)
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE ParseVars  ***
      !!
      !! ** Purpose :  Decode variables names to be used
      !!
      !! ** Method  :  look for , in the argument string and set the number of
      !!         variable (nvaro), allocate cv_fix array and fill it with the
      !!         decoded  names.
      !!
      !!----------------------------------------------------------------------
      CHARACTER(LEN=*), INTENT(in) :: cdum

      CHARACTER(LEN=80), DIMENSION(100) :: cl_dum  ! 100 is arbitrary
      INTEGER  :: ji
      INTEGER  :: inchar,  i1=1
      !!----------------------------------------------------------------------

      inchar= LEN(TRIM(cdum))
      nvar = 1
      ! scan the input string and look for ',' as separator
      DO ji=1,inchar
         IF ( cdum(ji:ji) == ',' ) THEN
            cl_dum(nvar) = cdum(i1:ji-1)
            i1=ji+1
            nvar=nvar+1
         ENDIF
      ENDDO
      ! reset all flags  to false
      ll_votemper=.FALSE.; ll_vosaline=.FALSE.; ll_sossheig=.FALSE.
      ll_vozocrtx=.FALSE.; ll_vomecrty=.FALSE.
      ll_ileadfra=.FALSE.; ll_iicethic=.FALSE.; ll_isnowthi=.FALSE.

      ! last name of the list does not have a ','
      cl_dum(nvar) = cdum(i1:inchar)

      DO ji =1, nvar
        SELECT CASE ( cl_dum(ji) )
        CASE ( 'votemper' ) ; ll_votemper=.TRUE.
        CASE ( 'vosaline' ) ; ll_vosaline=.TRUE.
        CASE ( 'sossheig' ) ; ll_sossheig=.TRUE.
        CASE ( 'vozocrtx' ) ; ll_vozocrtx=.TRUE.
        CASE ( 'vomecrty' ) ; ll_vomecrty=.TRUE.
        CASE ( 'ileadfra' ) ; ll_ileadfra=.TRUE.
        CASE ( 'iicethic' ) ; ll_iicethic=.TRUE.
        CASE ( 'isnowthi' ) ; ll_isnowthi=.TRUE.
        CASE DEFAULT 
         PRINT *,' +++ E R R O R : incorrect specified variable :', TRIM( cl_dum(ji))
        END SELECT
      ENDDO
   END SUBROUTINE ParseVars

END PROGRAM bdy_mk_coordinates_mpp
