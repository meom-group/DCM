PROGRAM icbrstsplit
  !!======================================================================
  !!                     ***  PROGRAM  icbrstsplit  ***
  !!=====================================================================
  !!  ** Purpose : Proceed with ICB restart file for changing the
  !!               domain decomposition
  !!
  !!  ** Method  : Read input restart files, look for files having icebergs,
  !!               In memory, recombine 2D fields and collates icebergs fields
  !!               Read target layout and select icebergs belonging to each target
  !!               subdomain, according to xi,yj. Write output files. For kount variable
  !!               set the value to imax_icb - nprocout + narea
  !!
  !! History :  1.0  : 07/2022  : J.M. Molines : 
  !!----------------------------------------------------------------------
  !!----------------------------------------------------------------------
  !!   routines      : description
  !!   CreateOutput  : Create and fill in output files 
  !!----------------------------------------------------------------------
  USE netcdf
  IMPLICIT NONE

  INTEGER(KIND=4)               :: jproc,  jb 
  INTEGER(KIND=4)               :: narg, ijarg, iargc
  INTEGER(KIND=4)               :: ncid, id, ierr, iuld,  nc, nk
  INTEGER(KIND=4)               :: nprocin, nprocout
  INTEGER(KIND=4)               :: npi,npj,ntot,nn, iloc
  INTEGER(KIND=4)               :: npil,npjl
  INTEGER(KIND=4)               :: imax_icb,imax_licb
  INTEGER(KIND=4)               :: nlci,nlcj,nldi,nldj,nlei,nlej,nimpp,njmpp
  INTEGER(KIND=4)               :: ii1,ii2,ij1,ij2, ii, ij
  INTEGER(KIND=4)               :: n1,n2
  INTEGER(KIND=4)               :: inum=10
  INTEGER(KIND=4), DIMENSION(2) :: idum
  INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: idum3  ! note that this 3 suppose nk=3 !!!

  ! Global array (full domain)
  REAL(KIND=8), DIMENSION(:,:),    ALLOCATABLE :: dcalving, dcalving_hflx
  REAL(KIND=8), DIMENSION(:,:,:),  ALLOCATABLE :: dstored_ice
  REAL(KIND=8), DIMENSION(:,:),    ALLOCATABLE :: dstored_heat

  REAL(KIND=8), DIMENSION(:),      ALLOCATABLE :: dlon, dlat,dxi,dyj,duvel,dvvel, dmass
  REAL(KIND=8), DIMENSION(:),      ALLOCATABLE :: dthickness, dwidth, dlength,day, dmass_scaling
  REAL(KIND=8), DIMENSION(:),      ALLOCATABLE :: dmass_of_bits, dheat_density

  INTEGER(KIND=4), DIMENSION(3)                :: nkount
  INTEGER(KIND=4), DIMENSION(:),   ALLOCATABLE :: nyear
  INTEGER(KIND=4), DIMENSION(:,:), ALLOCATABLE :: nnumber
  ! local domain
  REAL(KIND=8), DIMENSION(:,:),    ALLOCATABLE :: dlcalving, dlcalving_hflx
  REAL(KIND=8), DIMENSION(:,:,:),  ALLOCATABLE :: dlstored_ice
  REAL(KIND=8), DIMENSION(:,:),    ALLOCATABLE :: dlstored_heat

  REAL(KIND=8), DIMENSION(:),      ALLOCATABLE :: dllon, dllat,dlxi,dlyj,dluvel,dlvvel, dlmass
  REAL(KIND=8), DIMENSION(:),      ALLOCATABLE :: dlthickness, dlwidth, dllength,dlday, dlmass_scaling
  REAL(KIND=8), DIMENSION(:),      ALLOCATABLE :: dlmass_of_bits, dlheat_density
  INTEGER(KIND=4), DIMENSION(3)                :: ikount
  INTEGER(KIND=4), DIMENSION(:),   ALLOCATABLE :: iyear
  INTEGER(KIND=4), DIMENSION(:,:), ALLOCATABLE :: inumber

  CHARACTER(LEN=255) :: c_rootname
  CHARACTER(LEN=255) :: cf_rst
  CHARACTER(LEN=255) :: cf_layout
  CHARACTER(LEN=80)  :: cldum
  !!----------------------------------------------------------------------
  narg=iargc()

  IF ( narg == 0 ) THEN
     PRINT *,' usage : icbrstsplit.exe -r RST-root -l OUT-layout'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Change the domain decomposition for ICB restart file' 
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'      -r RST-root : give root name of input/output restart files'
     PRINT *,'                    (omiting _<RANK>.nc)' 
     PRINT *,'      -l OUT-layout : layout file (produced by NEM0) for the target'
     PRINT *,'                     domain decomposition.'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'        None so far'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'        none '
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'      Use <RST-root>_split for output file names'
     PRINT *,'     SEE ALSO :'
     PRINT *,'       '
     PRINT *,'      '
     STOP
  ENDIF

  ! Browse command line
  ijarg = 1 
  DO WHILE ( ijarg <= narg )
     CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1
     SELECT CASE ( cldum )
     CASE ( '-r'   ) ; CALL getarg(ijarg, c_rootname) ; ijarg=ijarg+1
     CASE ( '-l'   ) ; CALL getarg(ijarg, cf_layout ) ; ijarg=ijarg+1
        ! option
     CASE DEFAULT    ; PRINT *, ' ERROR : ', TRIM(cldum),' : unknown option.'; STOP 1
     END SELECT
  ENDDO

  PRINT *,'  RST-root : ', TRIM(c_rootname)
  PRINT *,'    LAYOUT : ', TRIM(cf_layout)
  ! Retrieve information for _0000 file:
  cf_rst=TRIM(c_rootname)//'_0000.nc'
  ierr = NF90_OPEN(cf_rst,NF90_NOWRITE,ncid)
  ierr = NF90_GET_ATT(ncid,NF90_GLOBAL,'DOMAIN_number_total',nprocin)
  ierr = NF90_GET_ATT(ncid,NF90_GLOBAL,'DOMAIN_size_global',idum) ; npi=idum(1) ; npj=idum(2)
  ierr = NF90_INQ_DIMID(ncid,'c',id)   ; ierr = NF90_INQUIRE_DIMENSION(ncid,id,len=nc)
  ierr = NF90_INQ_DIMID(ncid,'k',id)   ; ierr = NF90_INQUIRE_DIMENSION(ncid,id,len=nk)
  ALLOCATE( idum3(nk) )

  ierr = NF90_CLOSE(ncid)
  PRINT *, ' NUMBER OF ACTUAL DOMAINS : ', nprocin
  PRINT *, '                      NPI : ', npi
  PRINT *, '                      NPJ : ', npj
  ! Allocate 2D x,y fields
  ALLOCATE ( dcalving(npi,npj), dcalving_hflx(npi,npj), dstored_ice(npi,npj,nc), dstored_heat(npi,npj) )
  dcalving=0.d0 ; dcalving_hflx=0.d0 ; dstored_ice=0.d0 ; dstored_heat=0.d0
  nn=0
  ntot=0
  ! Look for total number of icebergs scanning once all input files
  DO jproc=1,nprocin
     WRITE(cf_rst,'(a,"_",i4.4,".nc" )') TRIM(c_rootname),jproc-1
     ierr = NF90_OPEN( cf_rst,NF90_NOWRITE,ncid)
     ! look for unlimited dimension
     ierr = NF90_INQUIRE(ncid,unlimitedDimId=iuld)
     IF (iuld /= -1 ) THEN
        ierr = NF90_INQUIRE_DIMENSION(ncid,iuld, len=iloc)
        nn = nn + iloc
     ENDIF
     ierr = NF90_CLOSE(ncid)
  ENDDO
  PRINT *, '  NN : ', nn
  ! Allocate global array with the total number of icebergs
  ALLOCATE ( dlon(nn), dlat(nn), dxi(nn),dyj(nn), duvel(nn),dvvel(nn), dmass(nn)  )
  ALLOCATE ( dthickness(nn), dwidth(nn), dlength(nn),day(nn), dmass_scaling(nn) ) 
  ALLOCATE ( dmass_of_bits(nn), dheat_density(nn) )
  ALLOCATE ( nyear(nn), nnumber(nk,nn))

  n1=1
  DO jproc=1,nprocin
     WRITE(cf_rst,'(a,"_",i4.4,".nc" )') TRIM(c_rootname),jproc-1
     !  print *, TRIM(cf_rst)
     ierr = NF90_OPEN( cf_rst,NF90_NOWRITE,ncid)
     ierr = NF90_INQ_DIMID(ncid,'x',id) ; ierr = NF90_INQUIRE_DIMENSION(ncid,id, len=npil) ; nlci=npil
     ierr = NF90_INQ_DIMID(ncid,'y',id) ; ierr = NF90_INQUIRE_DIMENSION(ncid,id, len=npjl) ; nlcj=npjl
     ierr = NF90_GET_ATT(ncid,NF90_GLOBAL,'DOMAIN_position_first' ,idum) ; nimpp = idum(1)      ; njmpp=idum(2) 
     ierr = NF90_GET_ATT(ncid,NF90_GLOBAL,'DOMAIN_halo_size_start',idum) ; nldi  = 1+idum(1)    ; nldj=1+idum(2)
     ierr = NF90_GET_ATT(ncid,NF90_GLOBAL,'DOMAIN_halo_size_end  ',idum) ; nlei  = nlci-idum(1) ; nlej=nlcj-idum(2)
     ALLOCATE ( dlcalving(npil,npjl), dlcalving_hflx(npil,npjl), dlstored_ice(npil,npjl,nc), dlstored_heat(npil,npjl) )
     ierr = NF90_INQ_VARID(ncid,'kount',id)
     ierr = NF90_GET_VAR(ncid,id,idum3)  ! not used for the moment
     ! Process 2D xy fields : calving, calving_hflx,stored_ice, stored_heat
     ierr = NF90_INQ_VARID(ncid,'calving',      id) ;  ierr = NF90_GET_VAR(ncid,id, dlcalving      ) 
     ierr = NF90_INQ_VARID(ncid,'calving_hflx', id) ;  ierr = NF90_GET_VAR(ncid,id, dlcalving_hflx )
     ierr = NF90_INQ_VARID(ncid,'stored_ice',   id) ;  ierr = NF90_GET_VAR(ncid,id, dlstored_ice   ) 
     ierr = NF90_INQ_VARID(ncid,'stored_heat',  id) ;  ierr = NF90_GET_VAR(ncid,id, dlstored_heat  ) 

     ii1 = nimpp+nldi-1 ; ii2= nimpp+nlei-1
     ij1 = njmpp+nldj-1 ; ij2= njmpp+nlej-1
     dcalving     (ii1:ii2,ij1:ij2)   = dlcalving     (nldi:nlei,nldj:nlej  )
     dcalving_hflx(ii1:ii2,ij1:ij2)   = dlcalving_hflx(nldi:nlei,nldj:nlej  )
     dstored_ice  (ii1:ii2,ij1:ij2,:) = dlstored_ice  (nldi:nlei,nldj:nlej,:)
     dstored_heat (ii1:ii2,ij1:ij2)   = dlstored_heat (nldi:nlei,nldj:nlej  )

     ! look for unlimited dimension
     ierr = NF90_INQUIRE(ncid,unlimitedDimId=iuld)
     IF (iuld /= -1 ) THEN
        ! Process files with icebergs
        ierr = NF90_INQUIRE_DIMENSION(ncid,iuld, len=iloc)
        ALLOCATE ( dllon(iloc), dllat(iloc), dlxi(iloc),dlyj(iloc), dluvel(iloc),dlvvel(iloc), dlmass(iloc)  )
        ALLOCATE ( dlthickness(iloc), dlwidth(iloc), dllength(iloc),dlday(iloc), dlmass_scaling(iloc) ) 
        ALLOCATE ( dlmass_of_bits(iloc), dlheat_density(iloc) )
        ALLOCATE ( iyear(iloc), inumber(nk,iloc) )

        ierr = NF90_INQ_VARID(ncid,'number',id)
        ierr = NF90_GET_VAR(ncid,id,inumber)
        imax_licb = MAXVAL(inumber(1,:))
        imax_icb  = MAX(imax_icb, imax_licb)
!       PRINT *, jproc, iloc, idum3(1), jproc-nprocin, MOD(idum3(1) - jproc, nprocin ), inumber(1,iloc) - nprocin + jproc
        ierr = NF90_INQ_VARID(ncid,'lon',id)          ; ierr = NF90_GET_VAR(ncid,id,dllon) 
        ierr = NF90_INQ_VARID(ncid,'lat',id)          ; ierr = NF90_GET_VAR(ncid,id,dllat) 
        ierr = NF90_INQ_VARID(ncid,'xi',id)           ; ierr = NF90_GET_VAR(ncid,id,dlxi) 
        ierr = NF90_INQ_VARID(ncid,'yj',id)           ; ierr = NF90_GET_VAR(ncid,id,dlyj) 
        ierr = NF90_INQ_VARID(ncid,'uvel',id)         ; ierr = NF90_GET_VAR(ncid,id,dluvel) 
        ierr = NF90_INQ_VARID(ncid,'vvel',id)         ; ierr = NF90_GET_VAR(ncid,id,dlvvel) 
        ierr = NF90_INQ_VARID(ncid,'mass',id)         ; ierr = NF90_GET_VAR(ncid,id,dlmass) 
        ierr = NF90_INQ_VARID(ncid,'thickness',id)    ; ierr = NF90_GET_VAR(ncid,id,dlthickness) 
        ierr = NF90_INQ_VARID(ncid,'width',id)        ; ierr = NF90_GET_VAR(ncid,id,dlwidth) 
        ierr = NF90_INQ_VARID(ncid,'length',id)       ; ierr = NF90_GET_VAR(ncid,id,dllength) 
        ierr = NF90_INQ_VARID(ncid,'year',id)         ; ierr = NF90_GET_VAR(ncid,id,iyear) 
        ierr = NF90_INQ_VARID(ncid,'day',id)          ; ierr = NF90_GET_VAR(ncid,id,dlday) 
        ierr = NF90_INQ_VARID(ncid,'mass_scaling',id) ; ierr = NF90_GET_VAR(ncid,id,dlmass_scaling) 
        ierr = NF90_INQ_VARID(ncid,'mass_of_bits',id) ; ierr = NF90_GET_VAR(ncid,id,dlmass_of_bits) 
        ierr = NF90_INQ_VARID(ncid,'heat_density',id) ; ierr = NF90_GET_VAR(ncid,id,dlheat_density) 

        ! now spread local value on global variable
        n2                   = n1+iloc-1
        dlon(n1:n2)          = dllon
        dlat(n1:n2)          = dllat
        ! xi, yj already in global coordinates
        dxi(n1:n2)           = dlxi
        dyj(n1:n2)           = dlyj
        duvel(n1:n2)         = dluvel
        dvvel(n1:n2)         = dlvvel
        dmass(n1:n2)         = dlmass
        dthickness(n1:n2)    = dlthickness
        dwidth(n1:n2)        = dlwidth
        dlength(n1:n2)       = dllength
        nyear(n1:n2)         = iyear
        day(n1:n2)           = dlday
        dmass_scaling(n1:n2) = dlmass_scaling
        dmass_of_bits(n1:n2) = dlmass_of_bits
        dheat_density(n1:n2) = dheat_density
        nnumber(:,n1:n2)     = inumber
        n1=n2+1
        DEALLOCATE ( dllon, dllat, dlxi,dlyj, dluvel,dlvvel, dlmass )
        DEALLOCATE ( dlthickness, dlwidth, dllength,dlday, dlmass_scaling ) 
        DEALLOCATE ( dlmass_of_bits, dlheat_density )
        DEALLOCATE ( inumber, iyear )
     ENDIF  ! files with  bergs
     ierr = NF90_CLOSE(ncid)
     DEALLOCATE ( dlcalving, dlcalving_hflx, dlstored_ice, dlstored_heat )
  ENDDO
  PRINT *, 'IMAX_ICB : ', imax_icb
  ! Now ready to write splitted file
  OPEN(inum, file=cf_layout)
  READ(inum,*)
  READ(inum,*) nprocout
  READ(inum,*)
  DO jproc=1,nprocout  ! Loop on target sub domains
     READ(inum,*) nlci,nlci,nlcj,nldi,nldj,nlei,nlej,nimpp,njmpp
     ALLOCATE( dlcalving(nlci,nlcj), dlcalving_hflx(nlci,nlcj), dlstored_ice(nlci,nlcj,nc), dlstored_heat(nlci,nlcj) )
     ii1=nimpp ; ii2=nimpp+nlci-1
     ij1=njmpp ; ij2=njmpp+nlcj-1
     dlcalving(:,:)      = dcalving     (ii1:ii2,ij1:ij2  )
     dlcalving_hflx(:,:) = dcalving_hflx(ii1:ii2,ij1:ij2  )
     dlstored_ice(:,:,:) = dstored_ice  (ii1:ii2,ij1:ij2,:)
     dlstored_heat(:,:)  = dstored_heat  (ii1:ii2,ij1:ij2 )
     ! check if there are icebergs in this file
     iloc=0
     DO jb=1, nn  ! scan the list of icebergs globally
        ii= INT(dxi(jb)+0.5)
        ij= INT(dyj(jb)+0.5)
        IF ( ii >= nldi+nimpp-1 .AND. ii <= nlei+nimpp-1 .AND. &
          &  ij >= nldj+njmpp-1 .AND. ij <= nlej+njmpp-1 ) THEN
           iloc=iloc + 1 
        ENDIF
     ENDDO
     ntot=ntot+iloc  ! for later checking that we do not loose icebergs in the process
     ! allocate local array
     ALLOCATE ( dllon(iloc), dllat(iloc), dlxi(iloc),dlyj(iloc), dluvel(iloc),dlvvel(iloc), dlmass(iloc)  )
     ALLOCATE ( dlthickness(iloc), dlwidth(iloc), dllength(iloc),dlday(iloc), dlmass_scaling(iloc) ) 
     ALLOCATE ( dlmass_of_bits(iloc), dlheat_density(iloc) )
     ALLOCATE ( iyear(iloc), inumber(nk,iloc) )
     iloc=0
     DO jb=1, nn
        ii= INT(dxi(jb)+0.5)
        ij= INT(dyj(jb)+0.5)
        IF ( ii .GE. nldi+nimpp-1 .AND. ii .LE. nlei+nimpp-1 .AND. &
             &     ij .GE. nldj+njmpp-1 .AND. ij .LE. nlej+njmpp-1 ) THEN
           iloc=iloc + 1 
           dllon(iloc)         = dlon(jb)          ; dllat(iloc)         = dlat(jb)
           dlxi(iloc)          = dxi(jb)           ; dlyj(iloc)          = dyj(jb)
           dluvel(iloc)        = duvel(jb)         ; dlvvel(iloc)        = dvvel(jb)
           dlmass(iloc)        = dmass(jb)         ; dlthickness(iloc)   = dthickness(jb)
           dlwidth(iloc)       = dwidth(jb)        ; dllength(iloc)      = dlength(jb)
           dlday(iloc)         = day(jb)           ; dlmass_scaling(iloc)= dmass_scaling(jb)
           dlmass_of_bits(iloc)= dmass_of_bits(jb) ; dlheat_density(iloc)= dheat_density(jb)
           iyear(iloc)         = nyear(jb)         ; inumber(:,iloc)     = nnumber(:,jb)
        ENDIF
     ENDDO
     ! write netcdf file
     WRITE(cf_rst,'(a,"_split_",i4.4,".nc" )') TRIM(c_rootname),jproc-1
     CALL CreateOutput(cf_rst)
     DEALLOCATE ( dllon, dllat, dlxi,dlyj, dluvel,dlvvel, dlmass )
     DEALLOCATE ( dlthickness, dlwidth, dllength,dlday, dlmass_scaling ) 
     DEALLOCATE ( dlmass_of_bits, dlheat_density )
     DEALLOCATE ( inumber, iyear )

     DEALLOCATE ( dlcalving, dlcalving_hflx, dlstored_ice, dlstored_heat )
  ENDDO
  PRINT *," NTOT (out) : ", ntot
  CLOSE(inum)

CONTAINS
  SUBROUTINE CreateOutput(cd_file)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE  CreateOutput ***
    !!
    !! ** Purpose :  Create and fill all target output icebergs restart file
    !!
    !! ** Method  :  This routine is called once for each output file.
    !!               It basically follow NEMO coding for easy comparison,
    !!               although not fully Dr Norm compliant...
    !!
    !!----------------------------------------------------------------------
    CHARACTER(LEN=*), INTENT(in) :: cd_file

    INTEGER(KIND=4) :: nret
    INTEGER(KIND=4) :: ix_dim, iy_dim, nc_dim, ik_dim, in_dim
    INTEGER(KIND=4) :: nkountid,ncalvid,ncalvhid,nsiceid,nsheatid
    INTEGER(KIND=4) :: nlonid,nlatid,nxid,nyid,nuvelid,nvvelid,nmassid,nthicknessid,nwidthid,nlengthid
    INTEGER(KIND=4) :: numberid,nyearid,ndayid,nscaling_id,nmass_of_bits_id,nheat_density_id
    !!----------------------------------------------------------------------

    nret = NF90_CREATE(TRIM(cf_rst),NF90_NETCDF4,ncid)
    nret = NF90_DEF_DIM(ncid, 'x', nlci, ix_dim)
    nret = NF90_DEF_DIM(ncid, 'y', nlcj, iy_dim)
    nret = NF90_DEF_DIM(ncid, 'c', nc, nc_dim)
    nret = NF90_DEF_DIM(ncid, 'k', nk , ik_dim)

    nret = NF90_PUT_ATT( ncid, NF90_GLOBAL, 'DOMAIN_number_total'   , nprocout           )
    nret = NF90_PUT_ATT( ncid, NF90_GLOBAL, 'DOMAIN_number'         , jproc-1            )
    nret = NF90_PUT_ATT( ncid, NF90_GLOBAL, 'DOMAIN_dimensions_ids' , (/1     , 2     /) )
    nret = NF90_PUT_ATT( ncid, NF90_GLOBAL, 'DOMAIN_size_global'    , (/npi, npj/)       )
    nret = NF90_PUT_ATT( ncid, NF90_GLOBAL, 'DOMAIN_size_local'     , (/nlci   , nlcj /) )
    nret = NF90_PUT_ATT( ncid, NF90_GLOBAL, 'DOMAIN_position_first' , (/nimpp , njmpp /) )
    nret = NF90_PUT_ATT( ncid, NF90_GLOBAL, 'DOMAIN_position_last'  , (/nimpp + nlci - 1 , njmpp + nlcj - 1  /) )
    nret = NF90_PUT_ATT( ncid, NF90_GLOBAL, 'DOMAIN_halo_size_start', (/nldi - 1         , nldj - 1          /) )
    nret = NF90_PUT_ATT( ncid, NF90_GLOBAL, 'DOMAIN_halo_size_end'  , (/nlci - nlei      , nlcj - nlej       /) )
    nret = NF90_PUT_ATT( ncid, NF90_GLOBAL, 'DOMAIN_type'           , 'BOX'              )
    IF (iloc /= 0 ) THEN  ! there are icebergs in this file
       nret = NF90_DEF_DIM(ncid, 'n', NF90_UNLIMITED, in_dim)
    ENDIF
    ! Variables always present:
    nret = NF90_DEF_VAR(ncid, 'kount'       , NF90_INT   , (/ ik_dim /), nkountid)
    nret = NF90_DEF_VAR(ncid, 'calving'     , NF90_DOUBLE, (/ ix_dim, iy_dim /), ncalvid)
    nret = NF90_DEF_VAR(ncid, 'calving_hflx', NF90_DOUBLE, (/ ix_dim, iy_dim /), ncalvhid)
    nret = NF90_DEF_VAR(ncid, 'stored_ice'  , NF90_DOUBLE, (/ ix_dim, iy_dim, nc_dim /), nsiceid)
    nret = NF90_DEF_VAR(ncid, 'stored_heat' , NF90_DOUBLE, (/ ix_dim, iy_dim /), nsheatid)

    ! Attributes
    nret = NF90_PUT_ATT(ncid, ncalvid , 'long_name', 'iceberg calving')
    nret = NF90_PUT_ATT(ncid, ncalvid , 'units', 'some')
    nret = NF90_PUT_ATT(ncid, ncalvhid, 'long_name', 'heat flux associated with iceberg calving')
    nret = NF90_PUT_ATT(ncid, ncalvhid, 'units', 'some')
    nret = NF90_PUT_ATT(ncid, nsiceid , 'long_name', 'stored ice used to calve icebergs')
    nret = NF90_PUT_ATT(ncid, nsiceid , 'units', 'kg/s')
    nret = NF90_PUT_ATT(ncid, nsheatid, 'long_name', 'heat in stored ice used to calve icebergs')
    nret = NF90_PUT_ATT(ncid, nsheatid, 'units', 'J/kg/s')

    IF (iloc /= 0 ) THEN
       ! Only add berg variables for this PE if we have anything to say

       ! Variables only in files containing bergs
       nret = NF90_DEF_VAR(ncid, 'lon'         , NF90_DOUBLE, in_dim, nlonid)
       nret = NF90_DEF_VAR(ncid, 'lat'         , NF90_DOUBLE, in_dim, nlatid)
       nret = NF90_DEF_VAR(ncid, 'xi'          , NF90_DOUBLE, in_dim, nxid)
       nret = NF90_DEF_VAR(ncid, 'yj'          , NF90_DOUBLE, in_dim, nyid)
       nret = NF90_DEF_VAR(ncid, 'uvel'        , NF90_DOUBLE, in_dim, nuvelid)
       nret = NF90_DEF_VAR(ncid, 'vvel'        , NF90_DOUBLE, in_dim, nvvelid)
       nret = NF90_DEF_VAR(ncid, 'mass'        , NF90_DOUBLE, in_dim, nmassid)
       nret = NF90_DEF_VAR(ncid, 'thickness'   , NF90_DOUBLE, in_dim, nthicknessid)
       nret = NF90_DEF_VAR(ncid, 'width'       , NF90_DOUBLE, in_dim, nwidthid)
       nret = NF90_DEF_VAR(ncid, 'length'      , NF90_DOUBLE, in_dim, nlengthid)
       nret = NF90_DEF_VAR(ncid, 'number'      , NF90_INT   , (/ik_dim,in_dim/), numberid)
       nret = NF90_DEF_VAR(ncid, 'year'        , NF90_INT   , in_dim, nyearid)
       nret = NF90_DEF_VAR(ncid, 'day'         , NF90_DOUBLE, in_dim, ndayid)
       nret = NF90_DEF_VAR(ncid, 'mass_scaling', NF90_DOUBLE, in_dim, nscaling_id)
       nret = NF90_DEF_VAR(ncid, 'mass_of_bits', NF90_DOUBLE, in_dim, nmass_of_bits_id)
       nret = NF90_DEF_VAR(ncid, 'heat_density', NF90_DOUBLE, in_dim, nheat_density_id)

       ! Attributes
       nret = NF90_PUT_ATT(ncid, nlonid, 'long_name', 'longitude')
       nret = NF90_PUT_ATT(ncid, nlonid, 'units', 'degrees_E')
       nret = NF90_PUT_ATT(ncid, nlatid, 'long_name', 'latitude')
       nret = NF90_PUT_ATT(ncid, nlatid, 'units', 'degrees_N')
       nret = NF90_PUT_ATT(ncid, nxid, 'long_name', 'x grid box position')
       nret = NF90_PUT_ATT(ncid, nxid, 'units', 'fractional')
       nret = NF90_PUT_ATT(ncid, nyid, 'long_name', 'y grid box position')
       nret = NF90_PUT_ATT(ncid, nyid, 'units', 'fractional')
       nret = NF90_PUT_ATT(ncid, nuvelid, 'long_name', 'zonal velocity')
       nret = NF90_PUT_ATT(ncid, nuvelid, 'units', 'm/s')
       nret = NF90_PUT_ATT(ncid, nvvelid, 'long_name', 'meridional velocity')
       nret = NF90_PUT_ATT(ncid, nvvelid, 'units', 'm/s')
       nret = NF90_PUT_ATT(ncid, nmassid, 'long_name', 'mass')
       nret = NF90_PUT_ATT(ncid, nmassid, 'units', 'kg')
       nret = NF90_PUT_ATT(ncid, nthicknessid, 'long_name', 'thickness')
       nret = NF90_PUT_ATT(ncid, nthicknessid, 'units', 'm')
       nret = NF90_PUT_ATT(ncid, nwidthid, 'long_name', 'width')
       nret = NF90_PUT_ATT(ncid, nwidthid, 'units', 'm')
       nret = NF90_PUT_ATT(ncid, nlengthid, 'long_name', 'length')
       nret = NF90_PUT_ATT(ncid, nlengthid, 'units', 'm')
       nret = NF90_PUT_ATT(ncid, numberid, 'long_name', 'iceberg number on this processor')
       nret = NF90_PUT_ATT(ncid, numberid, 'units', 'count')
       nret = NF90_PUT_ATT(ncid, nyearid, 'long_name', 'calendar year of calving event')
       nret = NF90_PUT_ATT(ncid, nyearid, 'units', 'years')
       nret = NF90_PUT_ATT(ncid, ndayid, 'long_name', 'year day of calving event')
       nret = NF90_PUT_ATT(ncid, ndayid, 'units', 'days')
       nret = NF90_PUT_ATT(ncid, nscaling_id, 'long_name', 'scaling factor for mass of calving berg')
       nret = NF90_PUT_ATT(ncid, nscaling_id, 'units', 'none')
       nret = NF90_PUT_ATT(ncid, nmass_of_bits_id, 'long_name', 'mass of bergy bits')
       nret = NF90_PUT_ATT(ncid, nmass_of_bits_id, 'units', 'kg')
       nret = NF90_PUT_ATT(ncid, nheat_density_id, 'long_name', 'heat density')
       nret = NF90_PUT_ATT(ncid, nheat_density_id, 'units', 'J/kg')

    ENDIF ! associated(first_berg)
    nret =  NF90_ENDDEF(ncid)
    ! kount is set so that there cannot be duplicate names (from NEMO icbrst code)
    idum3(1) =  imax_icb - nprocout + jproc
    idum3(2:nk) = 0
    nret = NF90_PUT_VAR( ncid, nkountid, idum3 )
    nret = NF90_PUT_VAR( ncid, nsheatid, dlstored_heat(:,:) )
    nret = NF90_PUT_VAR( ncid, ncalvid , dlcalving(:,:) )
    nret = NF90_PUT_VAR( ncid, ncalvhid, dlcalving_hflx(:,:) )
    nret = NF90_PUT_VAR( ncid, nsiceid, dlstored_ice(:,:,:) )
    IF ( iloc /= 0 ) THEN
       nret = NF90_PUT_VAR(ncid, numberid, inumber)
       nret = NF90_PUT_VAR(ncid, nscaling_id, dlmass_scaling)

       nret = NF90_PUT_VAR(ncid, nlonid, dllon )
       nret = NF90_PUT_VAR(ncid, nlatid, dllat )
       nret = NF90_PUT_VAR(ncid, nxid, dlxi )
       nret = NF90_PUT_VAR(ncid, nyid, dlyj )
       nret = NF90_PUT_VAR(ncid, nuvelid, dluvel )
       nret = NF90_PUT_VAR(ncid, nvvelid, dlvvel )
       nret = NF90_PUT_VAR(ncid, nmassid, dlmass)
       nret = NF90_PUT_VAR(ncid, nthicknessid, dlthickness)
       nret = NF90_PUT_VAR(ncid, nwidthid, dlwidth )
       nret = NF90_PUT_VAR(ncid, nlengthid, dllength )
       nret = NF90_PUT_VAR(ncid, nyearid, iyear )
       nret = NF90_PUT_VAR(ncid, ndayid, dlday )
       nret = NF90_PUT_VAR(ncid, nmass_of_bits_id, dlmass_of_bits)
       nret = NF90_PUT_VAR(ncid, nheat_density_id, dlheat_density)
    ENDIF

    nret = NF90_CLOSE(ncid)
  END SUBROUTINE CreateOutput

END PROGRAM icbrstsplit
