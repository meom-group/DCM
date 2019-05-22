PROGRAM dcmtk_dom_doc
  !!======================================================================
  !!                     ***  PROGRAM  dcmtk_dom_doc  ***
  !!=====================================================================
  !!  ** Purpose : Add documentation to domain_cfg.nc files created for
  !!               NEMO4.
  !!
  !!  ** Method  : Define a namelist_cfg variable and add source_bathy, 
  !!               source_coord as global attributes.
  !!
  !! History :  1.0  : 05/2019  : J.M. Molines : 
  !!----------------------------------------------------------------------
  USE netcdf
  IMPLICIT NONE

  INTEGER :: narg, ijarg, iargc
  INTEGER :: inum=10, n, ji
  INTEGER :: ncid, id , ierr, idl, idlen, ilen
  INTEGER :: iformat

  CHARACTER(LEN=300), DIMENSION(:), ALLOCATABLE :: cv_namlist
  CHARACTER(LEN=300)                            :: cldum
  CHARACTER(LEN=80)                             :: cf_namlist
  CHARACTER(LEN=255)                            :: cf_bathy
  CHARACTER(LEN=255)                            :: cf_coord
  CHARACTER(LEN=255)                            :: cf_domcfg
  !!----------------------------------------------------------------------
  !! DCM4, MEOM 2019
  !! Copyright (c) 2019, J.-M. Molines
  !! Software governed by the CeCILL licence 
  !!----------------------------------------------------------------------
  narg = iargc()
  
  IF ( narg == 0 ) THEN
     PRINT *,' usage : dcmtk_dom_doc -b NAM-bathy -c NAM-coord -n NAMELIST-file '
     PRINT *,'                       -d DOMAIN_CFG-file'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'        Add information in the domain_cfg.nc file after its creation for'
     PRINT *,'        NEMO4. The additional information consists in a new netcdf variable'
     PRINT *,'        called namelist_cfg, holding the content of the used namelist_cfg,'
     PRINT *,'        and 2 additional global attributes : source_bathy and source_coord,'
     PRINT *,'        indicating the name of the used bathy_meter and coordinates files.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'        -b NAM-bathy : name of the bathymetric file.'
     PRINT *,'        -c NAM-coord : name of the coordinates file.'
     PRINT *,'        -n NAMELIST-file : name of the namelist_cfg.'
     PRINT *,'        -d DOMAIN_CFG-file : name of the domain_cfg file to document.'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'         NAMELIST-file and DOMAIN_CFG-file'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'         input DOMAIN_CFG-file is modified on output.'
     PRINT *,'      '
     PRINT *,'     SEE ALSO :'
     PRINT *,'      domain_cfg.exe '
     PRINT *,'      '
     STOP
  ENDIF

  ijarg = 1 
  DO WHILE ( ijarg <= narg )
     CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1
     SELECT CASE ( cldum )
     CASE ( '-b'   ) ; CALL getarg(ijarg, cf_bathy   ) ; ijarg=ijarg+1
     CASE ( '-c'   ) ; CALL getarg(ijarg, cf_coord   ) ; ijarg=ijarg+1
     CASE ( '-n'   ) ; CALL getarg(ijarg, cf_namlist ) ; ijarg=ijarg+1
     CASE ( '-d'   ) ; CALL getarg(ijarg, cf_domcfg  ) ; ijarg=ijarg+1
     CASE DEFAULT    ; PRINT *, ' ERROR : ', TRIM(cldum),' : unknown option.'; STOP 1
     END SELECT
  ENDDO
  
  ! Open namelist and count the number of lines before redifining namelist variable in domcfg file.
  OPEN(inum, file=cf_namlist)

  ! count the number of lines
  n=0 ; ilen=0
  DO
     READ(inum,'(a)',END=999) cldum   ! loop till the end of file
     ilen=MAX(ilen,LEN(TRIM(cldum)))
     n=n+1
  ENDDO

999 PRINT *,' Number of lines in ',TRIM(cf_namlist),' is : ', n

  ! allocate variable array
  ALLOCATE (cv_namlist(n) )

  REWIND (inum)
  DO ji = 1, n
     READ(inum,'(a)',END=999) cv_namlist(ji)
  ENDDO
  CLOSE(inum)
  
  ! Open domain_cfg file 
  ierr = NF90_OPEN(cf_domcfg,NF90_WRITE,ncid)
  ierr = NF90_REDEF(ncid)

  ierr=NF90_DEF_DIM(ncid,'nlines',n, idl)
  ierr=NF90_DEF_DIM(ncid,'nlen',ilen, idlen)

  ! note that namelist is a character 2D variables in the sense of netcdf 
  ! the program variable is thus an array of character with ilen char per line and n lines.
  ! If domain_cfg is netcdf4, then deflate the namelist
  ierr = NF90_INQUIRE(ncid,formatNum=iformat)
  IF ( iformat == NF90_FORMAT_NETCDF4 .OR. iformat == NF90_FORMAT_NETCDF4_CLASSIC ) THEN
     ierr=NF90_DEF_VAR(ncid,'namelist_cfg',NF90_CHAR,(/idlen,idl/), id, chunksizes=(/ilen,n/), deflate_level=9 )
  ELSE
     ierr=NF90_DEF_VAR(ncid,'namelist_cfg',NF90_CHAR,(/idlen,idl/), id)
  ENDIF
  ! add extra global attributes:
  ierr = NF90_PUT_ATT(ncid,NF90_GLOBAL,'source_bathy',cf_bathy)
  ierr = NF90_PUT_ATT(ncid,NF90_GLOBAL,'source_coord',cf_coord)

  ierr=NF90_ENDDEF(ncid)
  DO ji=1,n
     ierr=NF90_PUT_VAR(ncid,id,cv_namlist(ji)(1:ilen) ,start=(/1,ji/), count=(/ilen,1/) )
     PRINT *,NF90_STRERROR(ierr)
  ENDDO
  ierr=NF90_CLOSE(ncid)

END PROGRAM dcmtk_dom_doc
