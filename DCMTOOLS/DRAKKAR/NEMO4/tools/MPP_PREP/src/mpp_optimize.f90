PROGRAM mpp_optimize_jm
   !!======================================================================
   !!                     ***  PROGRAM  mpp_optimize  ***
   !!=====================================================================
   !!  ** Purpose : Propose possible domain decompositions for a given 
   !!               bathymetric file, which is particularly intersting when
   !!               we want to eliminate land-only domain. 
   !!               All solution are proposed and written to output file.
   !!               The ratio between the effective number of computed 
   !!               point and the total number of points in the domain is 
   !!               given and is probably a major criteria for choosing a 
   !!               domain decomposition.
   !!
   !!  ** Method  : Use mpp_init like code for seting up the decomposition
   !!               and evaluate the efficiency of the decomposition.
   !! History
   !!       original  : 95-12 (Imbard M) for OPA8.1, CLIPPER
   !!       f90       : 03-06 (Molines JM), namelist as input
   !!                 : 05-05 (Molines JM), bathy in ncdf
   !!                 : 13-03 (Molines JM), Nemo-like coding and license.
   !!                 : 18-10 (Mathiot  P), upgrade the NEMO 4.0
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   routines      : description
   !!----------------------------------------------------------------------


   !!----------------------------------------------------------------------
   !! MPP-PREP, MEOM 2013
   !! $Id$
   !! Copyright (c) 2013, J.-M. Molines
   !! Software governed by the CeCILL licence (Licence/MPP-PREPCeCILL.txt)
   !!----------------------------------------------------------------------
   USE netcdf

   IMPLICIT NONE

   INTEGER, PARAMETER :: nn_hls=1             !: overlap between processors
   INTEGER            :: jperio               !: Periodic conditions read in domain_cfg file
   INTEGER            :: irm, ijpjmin         !: northen most row of domain is smaller

   ! Namelist declaration and definition
   ! -----------------------------------
   INTEGER ::  nn_procmax  =250    !: maximum number of proc. (Read from namelist)
   INTEGER ::  nn_procmin  = 1     !: maximum number of proc. (Read from namelist)
   LOGICAL ::  ln_memchk = .FALSE. ! add a memory constraint if true (obsolete)
   NAMELIST /namproc/ nn_procmax, nn_procmin, ln_memchk
   !
   INTEGER ::  nn_jpk = 46   !: vertical levels 
   INTEGER ::  nn_izoom = 1  !: I zoom indicator
   INTEGER ::  nn_jzoom = 1  !: J zoom indicator
   NAMELIST /namspace/ nn_jpk, nn_izoom, nn_jzoom
   !
   ! Following variables are used only if ln_memchk=.true.
   REAL(KIND=4) ::  required_memory, rppmpt !: not in namelist working array
   REAL(KIND=4) ::  rn_ppmcal = 225000000. !: maximum memory of one processor for a 
   !: given machine (in 8 byte words)
   REAL(KIND=4) ::  rn_ppmin  = 0.4        !: minimum ratio to fill the memory
   REAL(KIND=4) ::  rn_ppmax = 0.9         !: maximum ratio to fill the memory
   NAMELIST /namparam/ rn_ppmcal, rn_ppmin, rn_ppmax
   !
   CHARACTER(LEN=80) :: cn_var='none'   !: Variable name of the bathymetry
   CHARACTER(LEN=80) :: cn_x='x'        !: X dimension name
   CHARACTER(LEN=80) :: cn_y='y'        !: Y dimension name
   CHARACTER(LEN=80) :: cn_fbathy       !: File name of the netcdf bathymetry (namelist)
   NAMELIST /namfile/ cn_fbathy, cn_var, cn_x, cn_y
   !
   CHARACTER(LEN=80) :: cn_fovdta     !: root file name for keep output
   NAMELIST /namkeep/ cn_fovdta
   !
   INTEGER            :: numnam = 4       ! logical unit for namelist
   INTEGER            :: numout = 10        ! logical unit for output
   INTEGER            :: npiglo, npjglo   ! domain size
   INTEGER            :: npidta, npjdta   ! domain size

   INTEGER            :: ji, jj, jni, jnj ! dummy loop index
   INTEGER            :: jni2, jnj2, jjc  ! dummy loop index
   INTEGER            :: narg, ijarg      ! browsing command line

   ! Decomposition related arrays (using same meaning than in NEMO)
   INTEGER, DIMENSION(:,:), ALLOCATABLE :: nlcit, nlcjt ,nimppt, njmppt
   INTEGER, DIMENSION(:)  , ALLOCATABLE :: nlei_ocea, nldi_ocea
   INTEGER, DIMENSION(:)  , ALLOCATABLE :: nlej_ocea, nldj_ocea
   INTEGER, DIMENSION(:)  , ALLOCATABLE :: nlei_land, nldi_land
   INTEGER, DIMENSION(:)  , ALLOCATABLE :: nlej_land, nldj_land
   INTEGER                              :: nimpp, njmpp
   INTEGER                              :: ireci, irecj
   INTEGER                              :: nlci, nlcj
   INTEGER                              :: iresti, irestj
   !
   INTEGER :: ioce, isurf             !: number of ocean points cumulated, per_proc
   INTEGER :: ioce_opt                !: number of ocean points cumulated for optimal case
   INTEGER :: nland, nocea, nvalid    !: number of land, ocean, memory_valid  procs 
   INTEGER :: nland_opt               !: optimal number of land procs
   INTEGER :: ii1, ii2, ij1, ij2      !: limit of subdomain in global domain
   INTEGER :: ipi,     ipj            !: size of sub domain
   INTEGER :: ipi_opt, ipj_opt        !: size of sub domain for optimal case
   INTEGER :: inf10,     inf30,     inf50      !: 
   INTEGER :: inf10_opt, inf30_opt, inf50_opt  !:  in optimal case
   INTEGER :: npni_opt, npnj_opt      !: optimal domain decomposition

   INTEGER :: iminproci, imaxproci    !: limits of the processor loop
   INTEGER :: iminprocj, imaxprocj    !: can be reduded to  nkeepi, nkeepj

   ! Saving criteria
   REAL(KIND=4) :: ratio_min=99999.   !: keep only decomposition with ration less than ratio_min
   INTEGER      :: nocea_min = 1      !: minimum number of ocean procs for saving
   INTEGER      :: nmodulo = 1        !: Only keep solution multiple of nmodulo
   LOGICAL      :: ll_criteria=.TRUE. !:
   !
   REAL(KIND=4)                              ::  oce_cover
   REAL(KIND=4)                              ::  oce_cover_min,     oce_cover_max,     ratio
   REAL(KIND=4)                              ::  oce_cover_min_opt, oce_cover_max_opt, ratio_opt
   REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE ::  tmask     ! npiglo x npjglo
   REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE ::  bathy     ! npidta x npjdta

   ! CDF stuff
   INTEGER :: ncid, istatus, id
   LOGICAL ::  ll_good = .FALSE.

   CHARACTER(LEN=80) :: cf_namlist='namelist'
   CHARACTER(LEN=80) :: cf_out='processor.layout'
   CHARACTER(LEN=80) :: cdum                       ! dummy character variable

   ! Keep stuff      
   LOGICAL ::  ll_keep = .FALSE.
   INTEGER :: nkeepi, nkeepj          !: for option -keep : the retained decomposition
  !

   !!----------------------------------------------------------------------
   narg=iargc()
   ijarg=1
   IF ( narg == 0 ) THEN
      PRINT *,' try mpp_optimize -h for instructions !'
      STOP
   ENDIF
   !
   DO WHILE ( ijarg <= narg )
      CALL getarg(ijarg,cdum) ; ijarg=ijarg+1
      SELECT CASE ( cdum )
      CASE ('-h') 
         PRINT *,'  usage : mpp_optimize [ -h ]  [-keep jpni jpnj] [ -o file out ] '
         PRINT *,'               [ -modulo val ] [-r ratio] [-minocean procs] -n namelist'
         PRINT *,'      '
         PRINT *,'     PURPOSE :'
         PRINT *,'         This program is build to optimize the domain beakdown into'
         PRINT *,'         subdomain for mpp computing.'
         PRINT *,'         Once the grid size, and the land/sea mask is known, it looks'
         PRINT *,'         for all the possibilities within a range of setting parameters'
         PRINT *,'         and determine the optimal.'
         PRINT *,''
         PRINT *,'         Optimization is done with respect to the maximum number of'
         PRINT *,'         sea processors and to the maximum numbers of procs (nn_procmax)'
         PRINT *,'                '
         PRINT *,'         Optional optimization can be performed taking into account'
         PRINT *,'         the maximum available processor memory rn_ppmcal. This is'
         PRINT *,'         activated if ln_memchk is set true in the namelist'
         PRINT *,'      '
         PRINT *,'         Additional criteria can be given on the command line to reduce'
         PRINT *,'         the amount of possible choices.'
         PRINT *,'      '
         PRINT *,'     ARGUMENTS :'
         PRINT *,'         -n namelist : indicate the name of the namelist to use'
         PRINT *,'      '
         PRINT *,'     OPTIONS :'
         PRINT *,'         -h : print this help message'
         PRINT *,'         -keep jpni jpnj : print a file suitable for plotting,'
         PRINT *,'                 corresponding to the given decomposition'
         PRINT *,'         -o output file : give the name of the output file'
         PRINT *,'                 default is ',TRIM(cf_out)
         PRINT *,'         -modulo val : only retain decomposition whose total number'
         PRINT *,'                 of util processors (sea) are a multiple of val'
         PRINT *,'         -r ratio : only retain decomposition with a ratio computed/global'
         PRINT *,'                 less or equal to the given ratio'
         PRINT *,'         -minocean procs : only retain decomposition with a number of '
         PRINT *,'                 ocean procs greater of equal to procs'
         PRINT *,'      '
         PRINT *,'     REQUIRED FILES :'
         PRINT *,'       A bathymetric file and an ad-hoc namelist are required.'
         PRINT *,'       The file name of the bathymetry is specified in the namelist'
         PRINT *,'      '
         PRINT *,'     OUTPUT : '
         PRINT *,'       ',TRIM(cf_out),' : an ascii file with all found possibilities'
         PRINT *,'      '
         PRINT *,'     SEE ALSO :'
         PRINT *,'       script screen.ksh helps a lot in the analysis of the output file.'
         PRINT *,'      '
         STOP
      CASE ('-n' )
         CALL getarg(ijarg,cf_namlist) ; ijarg=ijarg+1
      CASE ('-o' )
         CALL getarg(ijarg,cf_out) ; ijarg=ijarg+1
      CASE ('-keep' )
         ll_keep=.TRUE.
         CALL getarg(ijarg,cdum) ; ijarg=ijarg+1 ; READ( cdum,*) nkeepi
         CALL getarg(ijarg,cdum) ; ijarg=ijarg+1 ; READ( cdum,*) nkeepj
      CASE ('-modulo' )
         CALL getarg(ijarg,cdum) ; ijarg=ijarg+1 ; READ( cdum,*) nmodulo
      CASE ('-r' )
         CALL getarg(ijarg,cdum) ; ijarg=ijarg+1 ; READ( cdum,*) ratio_min
      CASE ('-minocean' )
         CALL getarg(ijarg,cdum) ; ijarg=ijarg+1 ; READ( cdum,*) nocea_min
      END SELECT
   ENDDO

   ! Open and read the namelist
   OPEN(numnam,FILE=cf_namlist)
   REWIND(numnam)
   READ(numnam,namspace)

   REWIND(numnam)
   READ(numnam,namfile)

   REWIND(numnam)
   READ(numnam,namparam)

   REWIND(numnam)
   READ(numnam,namproc)

   REWIND(numnam)
   READ(numnam,namkeep)  ! only used for -keep option but still ...
   CLOSE(numnam)

   ! estimated code size expressed in number of 3D arrays (valid for OPA8.1) to be tuned for OPA9.0/Nemo
   rppmpt = 55.+73./nn_jpk

   ! Open bathy file an allocate required memory
   INQUIRE( FILE=cn_fbathy, EXIST=ll_good )
   IF( ll_good ) THEN
      istatus = NF90_OPEN(cn_fbathy, NF90_NOWRITE, ncid)
      istatus = NF90_INQ_DIMID(ncid, cn_x, id) ; istatus = NF90_INQUIRE_DIMENSION(ncid, id, len=npiglo)
      istatus = NF90_INQ_DIMID(ncid, cn_y, id) ; istatus = NF90_INQUIRE_DIMENSION(ncid, id, len=npjglo)
      npidta  = npiglo ; npjdta=npjglo
   ELSE
      PRINT *,' File missing : ', TRIM(cn_fbathy)
      STOP 42
   ENDIF

   ALLOCATE (tmask(npiglo,npjglo), bathy(npidta,npjdta) )
   ALLOCATE (nlcit(nn_procmax,nn_procmax), nlcjt(nn_procmax,nn_procmax) )
   ALLOCATE (nimppt(nn_procmax,nn_procmax), njmppt(nn_procmax,nn_procmax) )

   ! Open output file for results
   IF ( ll_keep ) THEN
      nn_procmax = nkeepi*nkeepj  ! reduce nn_procmax
      ! File will be open later
   ELSE
      OPEN(numout,FILE=cf_out)
      WRITE(numout,*)
      WRITE(numout,*) ' Domain decomposition optimization '
      WRITE(numout,*) ' ----------------------------------'
      WRITE(numout,*)
   ENDIF
   !
   ALLOCATE ( nlei_ocea(nn_procmax), nldi_ocea(nn_procmax), nlej_ocea(nn_procmax), nldj_ocea(nn_procmax) )
   ALLOCATE ( nlei_land(nn_procmax), nldi_land(nn_procmax), nlej_land(nn_procmax), nldj_land(nn_procmax) )
   !
   ! Read cdf bathy file
   IF ( cn_var == 'none' ) THEN  ! automatic detection of variable name according to partial step
       cn_var = 'bottom_level'
   ENDIF
   PRINT *,''
   PRINT *,' ocean/land file used is: ', TRIM(cn_fbathy)
   PRINT *,' variable used to find ocean domain is: ', TRIM(cn_var)
   PRINT *,' Dimensions (jpi x jpj) are: ',npiglo,'x',npjglo
   PRINT *,''

   istatus = NF90_INQ_VARID (ncid, cn_var, id)
   istatus = NF90_GET_VAR   (ncid, id,   bathy)
   istatus = NF90_INQ_VARID (ncid, 'jperio', id) 
   IF ( istatus /= NF90_NOERR) THEN
     PRINT *, " ERROR:  variable jperio not in ",TRIM(cn_fbathy)
     STOP 45
   ENDIF
   istatus = NF90_GET_VAR   (ncid, id,   jperio)
   istatus = NF90_CLOSE     (ncid)
   !
   ! Building the mask ( eventually on a smaller domain than the bathy)
   tmask(:,:) = bathy(nn_izoom:nn_izoom+npiglo -1,  nn_jzoom:nn_jzoom+npjglo -1)

   WHERE ( tmask > 0 ) 
      tmask = 1.
   ELSEWHERE
      tmask = 0.
   ENDWHERE

   !  Main loop on processors
   ! ------------------------
   ! initialization of working variables
   npni_opt=1       ; npnj_opt=1
   ipi_opt=npiglo ; ipj_opt=npjglo
   nland_opt=0   
   ioce_opt=0
   oce_cover_min_opt=0. ; oce_cover_max_opt=0.
   inf10_opt=0 ; inf30_opt=0 ; inf50_opt=0
   ratio_opt=1.

   nvalid=0       ! counter for valid case ( ln_memchk true )
   IF ( ll_keep ) THEN
      iminproci = nkeepi    ; imaxproci = iminproci
      iminprocj = nkeepj    ; imaxprocj = iminprocj
   ELSE
      iminproci = 1    ; imaxproci = MIN( nn_procmax, npiglo )
      iminprocj = 1    ; imaxprocj = MIN( nn_procmax, npjglo )
   ENDIF

   ! loop on all decomposition a priori
   PRINT *, 'Loop over all the decomposition (can take a while) ...'
   PRINT *, ''
   DO jni=iminproci, imaxproci
      DO jnj=iminprocj, imaxprocj
         ! Limitation of the maxumun number of PE's
         IF ( jni*jnj <=  nn_procmax .AND. jni*jnj >= nn_procmin )  THEN
            !
            !  1. Dimension arrays for subdomains
            ! -----------------------------------
            !
            ! Partition : size of sub-domain 
            ipi=(npiglo-2*nn_hls + (jni-1))/jni + 2*nn_hls
            ipj=(npjglo-2*nn_hls + (jnj-1))/jnj + 2*nn_hls
            !
            ! Memory optimization ?
            IF ( ln_memchk ) THEN
               required_memory=rppmpt*ipi*ipj*nn_jpk
               IF( required_memory > rn_ppmcal ) EXIT
               IF( required_memory > rn_ppmax*rn_ppmcal .OR. required_memory < rn_ppmin*rn_ppmcal) EXIT
            ENDIF
            nvalid=nvalid+1
            !
            ! Position of each sub domain   (jni x jni in total )
            nimpp  = 1                             ; njmpp  = 1
            ireci  = 2*nn_hls                      ; irecj  = 2*nn_hls
            iresti = 1 + MOD ( npiglo - ireci -1 , jni )  ; irestj = 1 + MOD ( npjglo - irecj -1 , jnj )
            !
            ! 
            nlcit(       1:iresti, 1:jnj) = ipi
            nlcit(iresti+1:jni   , 1:jnj) = ipi-1

            IF( jperio == 3 .OR. jperio == 4 .OR. jperio == 5 .OR. jperio == 6 ) THEN
             ! minimize the size of the last row to compensate for the north pole folding coast
               IF( jperio == 3 .OR. jperio == 4 )   ijpjmin = 5   ! V and F folding involves line jpj-3 that must not be south boundary
               IF( jperio == 5 .OR. jperio == 6 )   ijpjmin = 4   ! V and F folding involves line jpj-2 that must not be south boundary
               irm = jnj - irestj                                    ! total number of lines to be removed
               nlcjt(:,            jnj) = MAX( ijpjmin, ipj-irm )   ! we must have jpj >= ijpjmin in the last row
               irm = irm - ( ipj - nlcjt(1,jnj) )                   ! remaining number of lines to remove 
               irestj = jnj - 1 - irm
               nlcjt(1:jni,        1:irestj) = ipj
               nlcjt(1:jni, irestj+1:jnj-1 ) = ipj-1
            ELSE
               ijpjmin = 3
               nlcjt(1:jni,        1:irestj) = ipj
               nlcjt(1:jni, irestj+1:jnj   ) = ipj-1
            ENDIF

            nimppt(1:jni, 1:jnj) = nimpp
            njmppt(1:jni, 1:jnj) = njmpp

            IF( jni > 1 ) THEN
               DO jj=1,jnj
                  DO ji=2,jni
                     nimppt(ji,jj)=nimppt(ji-1,jj)+nlcit(ji-1,jj)-ireci
                  END DO
               END DO
            ENDIF

            IF( jnj > 1 ) THEN
               DO jj=2,jnj
                  DO ji=1,jni
                     njmppt(ji,jj)=njmppt(ji,jj-1)+nlcjt(ji,jj-1)-irecj
                  END DO
               END DO
            ENDIF
            !
            ! Loop on each subdomain to look for land proportion
            nland = 0
            nocea = 0
            ioce  = 0
            oce_cover_min = 1.e+20
            oce_cover_max = -1.e+20
            inf10=0
            inf30=0
            inf50=0
            !
            DO jnj2=1,jnj
              DO jni2=1,jni
                  nimpp = nimppt(jni2,jnj2)
                  njmpp = njmppt(jni2,jnj2)
                  nlci  = nlcit(jni2,jnj2)
                  nlcj  = nlcjt(jni2,jnj2)

                  isurf = 0
                  ! loop on inner point of sub-domain
                  !             DO jj=1+nn_hls, nlcjt(jni2,jnj2)-nn_hls
                  !                DO  ji=1+nn_hls, nlcit(jni2,jnj2)-nn_hls
                  !                   IF( tmask(ji+nimpp-1,jj+njmpp-1) == 1. ) isurf=isurf+1
                  !                END DO
                  !             END DO
!                 isurf = SUM( INT(tmask(ii1:ii2, ij1:ij2) ) )
                  isurf = SUM( INT(tmask(nimpp:nimpp+nlci-1, njmpp:njmpp+nlcj-1) ) )

                  ii1   = nimpp+nn_hls      ; ii2 = nimpp+nlci-1 -nn_hls
                  ij1   = njmpp+nn_hls      ; ij2 = njmpp+nlcj-1 -nn_hls

                  IF ( isurf == 0 ) THEN
                     nland = nland+1
                     nldi_land(nland) = ii1
                     nlei_land(nland) = ii2
                     nldj_land(nland) = ij1
                     nlej_land(nland) = ij2
                  ELSE
                     nocea = nocea+1
                     ioce  = ioce + isurf
                     nldi_ocea(nocea) = ii1
                     nlei_ocea(nocea) = ii2
                     nldj_ocea(nocea) = ij1
                     nlej_ocea(nocea) = ij2
                  ENDIF

                  ! ratio of wet points over total number of point per proc.
                  oce_cover = float(isurf)/float(ipi*ipj)

                  IF(oce_cover_min > oce_cover .AND. isurf /= 0) oce_cover_min=oce_cover
                  IF(oce_cover_max < oce_cover .AND. isurf /= 0) oce_cover_max=oce_cover
                  IF(oce_cover     < 0.1       .AND. isurf /= 0) inf10=inf10+1
                  IF(oce_cover     < 0.3       .AND. isurf /= 0) inf30=inf30+1
                  IF(oce_cover     < 0.5       .AND. isurf /= 0) inf50=inf50+1
                  !
               END DO  ! loop on processors
            END DO     ! loop on processors
            ! 
            ratio=float(nocea)*float(ipi*ipj)/float(npiglo*npjglo)

            ! criteria for printing results
            ll_criteria = ( ( MOD ( nocea, nmodulo ) == 0 ) .AND. &
                 &          ( ratio <= ratio_min          ) .AND. &
                 &          ( nocea >= nocea_min           )  )
            IF ( ll_keep ) THEN   ! the loop in done only once !
               WRITE(cdum,'(a,"-",i3.3,"x",i3.3,"_",i4.4)') TRIM(cn_fovdta), nkeepi, nkeepj, nocea
               OPEN(numout, file=cdum )
               WRITE(numout,'("# ocean ",i5)') nocea
               DO jjc=1, nocea
                  WRITE(numout,'("#",i5)') jjc
                  WRITE(numout,'(2i5)') nldi_ocea(jjc)-1+nn_izoom-1, nldj_ocea(jjc)-1+nn_jzoom -1
                  WRITE(numout,'(2i5)') nlei_ocea(jjc)+1+nn_izoom-1, nldj_ocea(jjc)-1+nn_jzoom -1
                  WRITE(numout,'(2i5)') nlei_ocea(jjc)+1+nn_izoom-1, nlej_ocea(jjc)+1+nn_jzoom -1
                  WRITE(numout,'(2i5)') nldi_ocea(jjc)-1+nn_izoom-1, nlej_ocea(jjc)+1+nn_jzoom -1
                  WRITE(numout,'(2i5)') nldi_ocea(jjc)-1+nn_izoom-1, nldj_ocea(jjc)-1+nn_jzoom -1
                  WRITE(numout,'(2i5)') 9999, 9999
               ENDDO
               !
               WRITE(numout,'("# land ",i5)') nland
               DO jjc=1, nland
                  WRITE(numout,'("# land ",i5)') jjc
                  WRITE(numout,'(2i5)') nldi_land(jjc)-1+nn_izoom-1, nldj_land(jjc)-1+nn_jzoom -1
                  WRITE(numout,'(2i5)') nlei_land(jjc)+1+nn_izoom-1, nldj_land(jjc)-1+nn_jzoom -1
                  WRITE(numout,'(2i5)') nlei_land(jjc)+1+nn_izoom-1, nlej_land(jjc)+1+nn_jzoom -1
                  WRITE(numout,'(2i5)') nldi_land(jjc)-1+nn_izoom-1, nlej_land(jjc)+1+nn_jzoom -1
                  WRITE(numout,'(2i5)') nldi_land(jjc)-1+nn_izoom-1, nldj_land(jjc)-1+nn_jzoom -1
                  WRITE(numout,'(2i5)') nlei_land(jjc)+1+nn_izoom-1, nlej_land(jjc)+1+nn_jzoom -1
                  WRITE(numout,'(2i5)') nldi_land(jjc)-1+nn_izoom-1, nlej_land(jjc)+1+nn_jzoom -1
                  WRITE(numout,'(2i5)') nlei_land(jjc)+1+nn_izoom-1, nldj_land(jjc)-1+nn_jzoom -1
                  WRITE(numout,'(2i5)') 9999, 9999
               ENDDO
               !
            ELSE
               IF ( ll_criteria ) THEN
                  WRITE(numout,*) ' iresti=',iresti,' irestj=',irestj
                  WRITE(numout,*) '--> Total number of domains ',jni*jnj
                  WRITE(numout,*) ' '
                  WRITE(numout,*) ' jpni=',jni ,' jpnj=',jnj
                  WRITE(numout,*) ' jpi= ',ipi ,' jpj= ',ipj
                  WRITE(numout,*) ' Number of ocean processors       ', nocea
                  WRITE(numout,*) ' Number of land processors        ', nland
                  WRITE(numout,*) ' Mean ocean coverage per domain   ', float(ioce)/float(nocea)/float(ipi*ipj)
                  WRITE(numout,*) ' Minimum ocean coverage           ', oce_cover_min
                  WRITE(numout,*) ' Maximum ocean coverage           ', oce_cover_max
                  WRITE(numout,*) ' nb of proc with coverage         < 10 % ', inf10
                  WRITE(numout,*) ' nb of proc with coverage 10 < nb < 30 % ', inf30 -inf10
                  WRITE(numout,*) ' nb of proc with coverage 30 < nb < 50 % ', inf50 -inf10 -inf30
                  WRITE(numout,*) ' Number of computed points        ', nocea*ipi*ipj
                  WRITE(numout,*) ' Overhead of computed points      ', nocea*ipi*ipj-npiglo*npjglo
                  WRITE(numout,*) ' % sup (computed / global)        ', ratio
                  WRITE(numout,*)
               ENDIF   ! note that indication of optimum does not take modulo into account (for information)
               ! 
               ! Look for optimum 
               IF( nland > nland_opt ) THEN
                  npni_opt          = jni
                  npnj_opt          = jnj
                  ipi_opt           = ipi
                  ipj_opt           = ipj
                  nland_opt         = nland
                  ioce_opt          = ioce
                  oce_cover_min_opt = oce_cover_min
                  oce_cover_max_opt = oce_cover_max
                  inf10_opt         = inf10
                  inf30_opt         = inf30
                  inf50_opt         = inf50
                  ratio_opt         = ratio
               ELSE IF( nland == nland_opt .AND. ratio_opt < ratio) THEN
                  npni_opt          = jni
                  npnj_opt          = jnj
                  ipi_opt           = ipi
                  ipj_opt           = ipj
                  ioce_opt          = ioce
                  oce_cover_min_opt = oce_cover_min
                  oce_cover_max_opt = oce_cover_max
                  inf10_opt         = inf10
                  inf30_opt         = inf30
                  inf50_opt         = inf50
                  ratio_opt         = ratio
               ENDIF
            ENDIF
         ENDIF
      END DO
   END DO
   !
   ! print optimal result
   IF ( .NOT. ll_keep ) THEN
      IF ( nvalid == 0 ) THEN
         WRITE(numout,*) ' no possible choice ...'
         WRITE(numout,*)
         WRITE(numout,*) 'insufficient number of processors for the available memory'
         STOP 
      ENDIF

      WRITE(numout,*) ' Optimal choice'
      WRITE(numout,*) ' =============='
      WRITE(numout,*) 
      WRITE(numout,*) '--> Total number of domains ',npni_opt*npnj_opt
      WRITE(numout,*) ' '
      WRITE(numout,*) ' jpni=',npni_opt ,' jpnj=',npnj_opt
      WRITE(numout,*) ' jpi= ',ipi_opt ,' jpj= ',ipj_opt
      WRITE(numout,*) 
      WRITE(numout,*) ' Number of ocean processors  ', npni_opt*npnj_opt-nland_opt
      WRITE(numout,*) ' Number of land processors   ', nland_opt
      WRITE(numout,*) ' Mean ocean coverage         ', float(ioce_opt)/float(npni_opt*npnj_opt-nland_opt)/float(ipi_opt*ipj_opt)
      WRITE(numout,*) ' Minimum ocean coverage      ', oce_cover_min_opt
      WRITE(numout,*) ' Maximum ocean coverage      ', oce_cover_max_opt
      WRITE(numout,*) ' nb de p recouvrement < 10 % ', inf10_opt
      WRITE(numout,*) ' nb de p      10 < nb < 30 % ', inf30_opt-inf10_opt
      WRITE(numout,*) ' nb de p      30 < nb < 50 % ', inf50_opt-inf10_opt -inf30_opt
      WRITE(numout,*) ' Number of computed points   ', (npni_opt*npnj_opt-nland_opt)*ipi_opt*ipj_opt
      WRITE(numout,*) ' Overhead of computed points ', (npni_opt*npnj_opt-nland_opt)*ipi_opt*ipj_opt-npiglo*npjglo
      WRITE(numout,*) ' % sup (computed / global)   ', ratio_opt
      WRITE(numout,*)
   ENDIF
   CLOSE(numout)
   !
   STOP
END PROGRAM mpp_optimize_jm
