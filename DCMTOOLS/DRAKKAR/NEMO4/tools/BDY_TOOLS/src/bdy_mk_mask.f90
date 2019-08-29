PROGRAM bdy_mk_mask
  !!======================================================================
  !!                     ***  PROGRAM  bdy_mk_mask  ***
  !!=====================================================================
  !!  ** Purpose : build a bdy_mask_file from tmask, and masking points
  !!               points corresponding to rim=0 in the bdy_coordinates file
  !!
  !!  ** Method  : Read tmask and the list of bdy_coordinates file
  !!               use nbit, nbjt and nbrt for masking rim-0 points.
  !!
  !! History :  1.0  : 08/2019  : J.M. Molines : 
  !!----------------------------------------------------------------------
  !!----------------------------------------------------------------------
  !!   routines      : description
  !!----------------------------------------------------------------------
  USE netcdf
  !!----------------------------------------------------------------------
  !! $Id$
  !! Copyright (c) 2012, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !!----------------------------------------------------------------------
  INTEGER(KIND=4)    :: jb, jbp
  INTEGER(KIND=4)    :: iargc, ijarg, narg
  INTEGER(KIND=4)    :: npiglo, npjglo, nxt, nbdy
  INTEGER(KIND=4)    :: ierr, ncid, id
  INTEGER(KIND=4)    :: ncidb, idb
  INTEGER(KIND=4), DIMENSION(2) :: idimx
  INTEGER(KIND=4), DIMENSION(:),   ALLOCATABLE :: nbit,nbjt,nbrt
  INTEGER(KIND=4), DIMENSION(:,:), ALLOCATABLE :: ibdymsk

  CHARACTER(LEN=255) :: cf_mask
  CHARACTER(LEN=255) :: cf_bdymsk='bdy_mask.nc'
  CHARACTER(LEN=255) :: cldum, cmd
  CHARACTER(LEN=255), DIMENSION(:), ALLOCATABLE :: cf_lst
  !!----------------------------------------------------------------------
  narg=iargc()

  IF ( narg == 0 ) THEN
     PRINT *,' usage :  bdy_mk_mask -m MASK-file -b BDY1-coord-file,BDY2-coord-file,BDY3....'
     PRINT *,'          [ -o BDYMASK-file]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'          Create a BDY mask file (ie mask points outside the BDY) using tmask'
     PRINT *,'          as a basis. Add land points where rim number=0.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'        -m MASK-file : mesh-mask file with at tmask for the domain'
     PRINT *,'        -b BDY1-coord-file,BDY2-coord-file,BDY3.. : comma separated list of '
     PRINT *,'              the bdy coordinates files used in the config.'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'        -o BDYMASK-file: name of the output file. Default :', TRIM(cf_bdymsk) 
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       only those specified on the command line' 
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_bdymsk) 
     PRINT *,'         variables : ', bdy_mask
     PRINT *,'      '
     PRINT *,'     SEE ALSO :'
     PRINT *,'      ' 
     PRINT *,'      '
     STOP
  ENDIF

  ijarg = 1 
  DO WHILE ( ijarg <= narg )
     CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1
     SELECT CASE ( cldum )
     CASE ( '-m'   ) ; CALL getarg(ijarg, cf_mask ) ; ijarg=ijarg+1
     CASE ( '-b'   ) ; CALL getarg(ijarg, cldum   ) ; ijarg=ijarg+1
        ;            ; CALL GetList(cldum) 
        ! option
     CASE ( '-o'   ) ; CALL getarg(ijarg, cf_bdymsk) ; ijarg=ijarg+1
     CASE DEFAULT    ; PRINT *, ' ERROR : ', TRIM(cldum),' : unknown option.'; STOP 1
     END SELECT
  ENDDO


  ! copy tmask(1) on the outputfile:
  cmd="ncks -O -6 -d nav_lev,1 -v tmask "//TRIM(cf_mask)//" "//TRIM(cf_bdymsk)
  PRINT *,TRIM(cmd)
  CALL SYSTEM( cmd)
  cmd="ncrename -v tmask,bdy_msk "//TRIM(cf_bdymsk)
  PRINT *,TRIM(cmd)
  CALL SYSTEM (cmd)
  cmd="ncwa -O -a nav_lev "//TRIM(cf_bdymsk)//" ztmp.nc"
  PRINT *,TRIM(cmd)
  CALL SYSTEM (cmd)
  cmd="mv ztmp.nc "//TRIM(cf_bdymsk)
  PRINT *,TRIM(cmd)
  CALL SYSTEM (cmd)
  !
  ! Read tmask
  ierr = NF90_OPEN(cf_bdymsk, NF90_WRITE, ncid)
  ierr = NF90_INQ_DIMID(ncid,'x', id ) ; ierr = NF90_INQUIRE_DIMENSION(ncid,id,len=npiglo) 
  ierr = NF90_INQ_DIMID(ncid,'y', id ) ; ierr = NF90_INQUIRE_DIMENSION(ncid,id,len=npjglo) 
  ALLOCATE ( ibdymsk(npiglo,npjglo) )
  ierr = NF90_INQ_VARID( ncid,'bdy_mask',id)
  ierr = NF90_GET_VAR(ncid,id,ibdymsk,start=(/1,1,1/), count=(/npiglo,npjglo,1/) )
  ! keep file open for later write
  DO jb=1, nbdy
  ! read bdy_coord file
  print *,TRIM(cf_lst(jb) )
  ierr = NF90_OPEN(cf_lst(jb),NF90_NOWRITE,ncidb)
  ierr = NF90_INQ_VARID( ncidb,'nbit',idb)
  ierr = NF90_INQUIRE_VARIABLE(ncidb,idb,dimids=idimx)
  ierr = NF90_INQUIRE_DIMENSION(ncidb,idimx(1),len=nxt)
  PRINT *,' NXT=', nxt
  ALLOCATE (nbit(nxt), nbjt(nxt), nbrt(nxt) )
  ierr = NF90_GET_VAR(ncidb,idb,nbit,start=(/1,1/), count=(/nxt,1/) )
  ierr = NF90_INQ_VARID( ncidb,'nbjt',idb)
  ierr = NF90_GET_VAR(ncidb,idb,nbjt,start=(/1,1/), count=(/nxt,1/) )
  ierr = NF90_INQ_VARID( ncidb,'nbrt',idb)
  ierr = NF90_GET_VAR(ncidb,idb,nbrt,start=(/1,1/), count=(/nxt,1/) )
  ierr = NF90_CLOSE(ncidb)
  DO jbp=1,nxt
    IF ( nbrt(jbp) == 0 ) THEN
      ii=nbit(jbp)
      ij=nbjt(jbp)
      ibdymsk(ii,ij) = 0
    ENDIF
  ENDDO
  DEALLOCATE ( nbit, nbjt, nbrt)
  ENDDO
  ! write back ibdymsk
  ierr = NF90_PUT_VAR(ncid,id,ibdymsk,start=(/1,1,1/),count=(/npiglo,npjglo,1/) )
  ierr = NF90_CLOSE(ncid)

  CONTAINS
  SUBROUTINE GetList(cdum)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE GetList  ***
    !!
    !! ** Purpose :  Set up a file list given on the command line as 
    !!               comma separated list
    !!
    !! ** Method  :  
    !!----------------------------------------------------------------------
    CHARACTER(LEN=*) , INTENT(in) :: cdum
    INTEGER (KIND=4)  :: ji, i1
    INTEGER (KIND=4)  :: icur
    INTEGER (KIND=4)  :: inchar
  
    !!----------------------------------------------------------------------
    !!
    print *, TRIM(cdum)
    nbdy=0
    inchar = LEN(TRIM(cdum))
    print *, inchar
    DO ji=1,inchar
       IF ( cdum(ji:ji) == ',' ) THEN
          nbdy=nbdy+1
       ENDIF
    ENDDO
    nbdy=nbdy+1  ! one more than ',' !
    print *, "nbdy", nbdy
    ALLOCATE( cf_lst(nbdy) )
    nbdy=1
    i1=1
    DO ji=1,inchar
       IF ( cdum(ji:ji) == ',' ) THEN
          cf_lst(nbdy) = cdum(i1:ji-1)
          i1=ji+1
          nbdy=nbdy+1
       ENDIF
    ENDDO
    cf_lst(nbdy)=cdum(i1:)


  END SUBROUTINE GetList

END PROGRAM bdy_mk_mask
