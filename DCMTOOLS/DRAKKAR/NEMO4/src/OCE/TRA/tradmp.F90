MODULE tradmp
   !!======================================================================
   !!                       ***  MODULE  tradmp  ***
   !! Ocean physics: internal restoring trend on active tracers (T and S)
   !!======================================================================
   !! History :  OPA  ! 1991-03  (O. Marti, G. Madec)  Original code
   !!                 ! 1992-06  (M. Imbard)  doctor norme
   !!                 ! 1998-07  (M. Imbard, G. Madec) ORCA version
   !!            7.0  ! 2001-02  (M. Imbard)  add distance to coast, Original code
   !!            8.1  ! 2001-02  (G. Madec, E. Durand)  cleaning
   !!  NEMO      1.0  ! 2002-08  (G. Madec, E. Durand)  free form + modules
   !!            3.2  ! 2009-08  (G. Madec, C. Talandier)  DOCTOR norm for namelist parameter
   !!            3.3  ! 2010-06  (C. Ethe, G. Madec) merge TRA-TRC 
   !!            3.4  ! 2011-04  (G. Madec, C. Ethe) Merge of dtatem and dtasal + suppression of CPP keys
   !!            3.6  ! 2015-06  (T. Graham)  read restoring coefficient in a file
   !!            3.7  ! 2015-10  (G. Madec)  remove useless trends arrays
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   tra_dmp_alloc : allocate tradmp arrays
   !!   tra_dmp       : update the tracer trend with the internal damping
   !!   tra_dmp_init  : initialization, namlist read, parameters control
   !!----------------------------------------------------------------------
   USE oce            ! ocean: variables
   USE dom_oce        ! ocean: domain variables
   USE c1d            ! 1D vertical configuration
   USE trd_oce        ! trends: ocean variables
   USE trdtra         ! trends manager: tracers 
   USE zdf_oce        ! ocean: vertical physics
   USE phycst         ! physical constants
   USE dtatsd         ! data: temperature & salinity
   USE zdfmxl         ! vertical physics: mixed layer depth
   !
   USE in_out_manager ! I/O manager
   USE iom            ! XIOS
   USE lib_mpp        ! MPP library
   USE prtctl         ! Print control
   USE timing         ! Timing
#if defined key_drakkar
   USE fldread , ONLY : FLD_N      ! for setting dmpmask 
#endif

   IMPLICIT NONE
   PRIVATE

   PUBLIC   tra_dmp        ! called by step.F90
   PUBLIC   tra_dmp_init   ! called by nemogcm.F90

   !                                           !!* Namelist namtra_dmp : T & S newtonian damping *
   LOGICAL            , PUBLIC ::   ln_tradmp   !: internal damping flag
   INTEGER            , PUBLIC ::   nn_zdmp     !: = 0/1/2 flag for damping in the mixed layer
   CHARACTER(LEN=200) , PUBLIC ::   cn_resto    !: name of netcdf file containing restoration coefficient field
   !
#if defined key_drakkar
   INTEGER            ::   nn_hdmp     ! > 0 standard NEMO CODE
                                       ! = -2 = DRAKKAR customization
   INTEGER            ::   nn_file     ! = 1 create a damping.coeff NetCDF file 
   LOGICAL            ::   ln_dmpmask  !  flag for using mask_dmp file
   REAL(wp)           ::   rn_timsk    !  restoring time scale used with mask_dmp       [days] 
   CHARACTER(LEN=255) ::   cn_dir      !  directory for damping mask
   TYPE(FLD_N)        ::   sn_dmp      ! structure for setting the damping file
#endif
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   resto    !: restoring coeff. on T and S (s-1)

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: tradmp.F90 11536 2019-09-11 13:54:18Z smasson $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION tra_dmp_alloc()
      !!----------------------------------------------------------------------
      !!                ***  FUNCTION tra_dmp_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( resto(jpi,jpj,jpk), STAT= tra_dmp_alloc )
      !
      CALL mpp_sum ( 'tradmp', tra_dmp_alloc )
      IF( tra_dmp_alloc > 0 )   CALL ctl_warn('tra_dmp_alloc: allocation of arrays failed')
      !
   END FUNCTION tra_dmp_alloc


   SUBROUTINE tra_dmp( kt )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE tra_dmp  ***
      !!                  
      !! ** Purpose :   Compute the tracer trend due to a newtonian damping
      !!      of the tracer field towards given data field and add it to the
      !!      general tracer trends.
      !!
      !! ** Method  :   Newtonian damping towards t_dta and s_dta computed 
      !!      and add to the general tracer trends:
      !!                     ta = ta + resto * (t_dta - tb)
      !!                     sa = sa + resto * (s_dta - sb)
      !!         The trend is computed either throughout the water column
      !!      (nlmdmp=0) or in area of weak vertical mixing (nlmdmp=1) or
      !!      below the well mixed layer (nlmdmp=2)
      !!
      !! ** Action  : - tsa: tracer trends updated with the damping trend
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      !
      INTEGER ::   ji, jj, jk, jn   ! dummy loop indices
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpts)     ::  zts_dta
      REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::  ztrdts
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('tra_dmp')
      !
      IF( l_trdtra )   THEN                    !* Save ta and sa trends
         ALLOCATE( ztrdts(jpi,jpj,jpk,jpts) ) 
         ztrdts(:,:,:,:) = tsa(:,:,:,:) 
      ENDIF
      !                           !==  input T-S data at kt  ==!
      CALL dta_tsd( kt, zts_dta )            ! read and interpolates T-S data at kt
      !
      SELECT CASE ( nn_zdmp )     !==  type of damping  ==!
      !
      CASE( 0 )                        !*  newtonian damping throughout the water column  *!
         DO jn = 1, jpts
            DO jk = 1, jpkm1
               DO jj = 2, jpjm1
                  DO ji = fs_2, fs_jpim1   ! vector opt.
                     tsa(ji,jj,jk,jn) = tsa(ji,jj,jk,jn) + resto(ji,jj,jk) * ( zts_dta(ji,jj,jk,jn) - tsb(ji,jj,jk,jn) )
                  END DO
               END DO
            END DO
         END DO
         !
      CASE ( 1 )                       !*  no damping in the turbocline (avt > 5 cm2/s)  *!
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  IF( avt(ji,jj,jk) <= avt_c ) THEN
                     tsa(ji,jj,jk,jp_tem) = tsa(ji,jj,jk,jp_tem)   &
                        &                 + resto(ji,jj,jk) * ( zts_dta(ji,jj,jk,jp_tem) - tsb(ji,jj,jk,jp_tem) )
                     tsa(ji,jj,jk,jp_sal) = tsa(ji,jj,jk,jp_sal)   &
                        &                 + resto(ji,jj,jk) * ( zts_dta(ji,jj,jk,jp_sal) - tsb(ji,jj,jk,jp_sal) )
                  ENDIF
               END DO
            END DO
         END DO
         !
      CASE ( 2 )                       !*  no damping in the mixed layer   *!
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  IF( gdept_n(ji,jj,jk) >= hmlp (ji,jj) ) THEN
                     tsa(ji,jj,jk,jp_tem) = tsa(ji,jj,jk,jp_tem)   &
                        &                 + resto(ji,jj,jk) * ( zts_dta(ji,jj,jk,jp_tem) - tsb(ji,jj,jk,jp_tem) )
                     tsa(ji,jj,jk,jp_sal) = tsa(ji,jj,jk,jp_sal)   &
                        &                 + resto(ji,jj,jk) * ( zts_dta(ji,jj,jk,jp_sal) - tsb(ji,jj,jk,jp_sal) )
                  ENDIF
               END DO
            END DO
         END DO
         !
      END SELECT
      !
      IF( l_trdtra )   THEN       ! trend diagnostic
         ztrdts(:,:,:,:) = tsa(:,:,:,:) - ztrdts(:,:,:,:)
         CALL trd_tra( kt, 'TRA', jp_tem, jptra_dmp, ztrdts(:,:,:,jp_tem) )
         CALL trd_tra( kt, 'TRA', jp_sal, jptra_dmp, ztrdts(:,:,:,jp_sal) )
         DEALLOCATE( ztrdts ) 
      ENDIF
      !                           ! Control print
      IF(ln_ctl)   CALL prt_ctl( tab3d_1=tsa(:,:,:,jp_tem), clinfo1=' dmp  - Ta: ', mask1=tmask,   &
         &                       tab3d_2=tsa(:,:,:,jp_sal), clinfo2=       ' Sa: ', mask2=tmask, clinfo3='tra' )
      !
      IF( ln_timing )   CALL timing_stop('tra_dmp')
      !
   END SUBROUTINE tra_dmp


   SUBROUTINE tra_dmp_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_dmp_init  ***
      !! 
      !! ** Purpose :   Initialization for the newtonian damping 
      !!
      !! ** Method  :   read the namtra_dmp namelist and check the parameters
      !!----------------------------------------------------------------------
      INTEGER ::   ios, imask   ! local integers 
      !
      NAMELIST/namtra_dmp/ ln_tradmp, nn_zdmp, cn_resto
#if defined key_drakkar
      INTEGER :: jk ! dummy loop index
      NAMELIST/namtra_dmp_drk/ nn_hdmp , nn_file, ln_dmpmask, rn_timsk, cn_dir, sn_dmp
#endif
      !!----------------------------------------------------------------------
      !
      REWIND( numnam_ref )   ! Namelist namtra_dmp in reference namelist : T & S relaxation
      READ  ( numnam_ref, namtra_dmp, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namtra_dmp in reference namelist' )
      !
      REWIND( numnam_cfg )   ! Namelist namtra_dmp in configuration namelist : T & S relaxation
      READ  ( numnam_cfg, namtra_dmp, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namtra_dmp in configuration namelist' )
      IF(lwm) WRITE ( numond, namtra_dmp )
#if defined key_drakkar
      REWIND( numnam_ref )   ! Namelist namtra_dmp in reference namelist : T & S relaxation
      READ  ( numnam_ref, namtra_dmp_drk, IOSTAT = ios, ERR = 903)
903   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namtra_dmp_drk in reference namelist' )
      !
      REWIND( numnam_cfg )   ! Namelist namtra_dmp in configuration namelist : T & S relaxation
      READ  ( numnam_cfg, namtra_dmp_drk, IOSTAT = ios, ERR = 904 )
904   IF( ios >  0 )   CALL ctl_nam ( ios , 'namtra_dmp_drk in configuration namelist' )
      IF(lwm) WRITE ( numond, namtra_dmp_drk )
#endif
      !
      IF(lwp) THEN                  ! Namelist print
         WRITE(numout,*)
         WRITE(numout,*) 'tra_dmp_init : T and S newtonian relaxation'
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namtra_dmp : set relaxation parameters'
         WRITE(numout,*) '      Apply relaxation   or not       ln_tradmp   = ', ln_tradmp
         WRITE(numout,*) '         mixed layer damping option      nn_zdmp  = ', nn_zdmp
         WRITE(numout,*) '         Damping file name               cn_resto = ', cn_resto
#if defined key_drakkar
         WRITE(numout,*) '      T and S damping option          nn_hdmp   = ', nn_hdmp
         WRITE(numout,*) '      use a mask_dmp file (T/F)      ln_dmpmask = ', ln_dmpmask
         WRITE(numout,*) '      time scale (mask_dmp)            rn_timsk = ', rn_timsk
         WRITE(numout,*) '      create a damping.coeff file     nn_file   = ', nn_file
#endif
         WRITE(numout,*)
      ENDIF
      !
      IF( ln_tradmp ) THEN
         !                          ! Allocate arrays
         IF( tra_dmp_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'tra_dmp_init: unable to allocate arrays' )
         !
#if defined key_drakkar
         SELECT CASE ( nn_hdmp )
         CASE (  -2  )   ;   IF(lwp) WRITE(numout,*) '   Drakkar customization '
         CASE DEFAULT
                             IF(lwp) WRITE(numout,*) '   Standard Nemo relaxation '
         END SELECT
#endif
         SELECT CASE (nn_zdmp)      ! Check values of nn_zdmp
         CASE ( 0 )   ;   IF(lwp) WRITE(numout,*) '   tracer damping as specified by mask'
         CASE ( 1 )   ;   IF(lwp) WRITE(numout,*) '   no tracer damping in the mixing layer (kz > 5 cm2/s)'
         CASE ( 2 )   ;   IF(lwp) WRITE(numout,*) '   no tracer damping in the mixed  layer'
         CASE DEFAULT
            CALL ctl_stop('tra_dmp_init : wrong value of nn_zdmp')
         END SELECT
         !
         !!TG: Initialisation of dtatsd - Would it be better to have dmpdta routine
         !    so can damp to something other than intitial conditions files?
         !!gm: In principle yes. Nevertheless, we can't anticipate demands that have never been formulated.
         ! JMM : DRAKKAR version is doing so ! 
         IF( .NOT.ln_tsd_dmp ) THEN
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout, *)  '   read T-S data not initialized, we force ln_tsd_dmp=T'
            CALL dta_tsd_init( ld_tradmp=ln_tradmp )        ! forces the initialisation of T-S data
         ENDIF
#if defined key_drakkar
         IF ( nn_hdmp == -2 ) THEN
            CALL dtacof( nn_file, 'TRA', resto )
         ELSE
#endif
         !                          ! Read in mask from file
         CALL iom_open ( cn_resto, imask)
         CALL iom_get  ( imask, jpdom_autoglo, 'resto', resto )
         CALL iom_close( imask )
#if defined key_drakkar
         ENDIF
#endif
      ENDIF
      !
   END SUBROUTINE tra_dmp_init

#if defined key_drakkar
   SUBROUTINE dtacof(kn_file, cdtype , presto )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dtacof  ***
      !!
      !! ** Purpose :   Compute the damping coefficient
      !!
      !! ** Method  :   Arrays defining the damping are computed for each grid
      !!                point for temperature and salinity (resto)
      !!                Damping depends on distance to coast, depth and latitude
      !!
      !! ** Action  : - resto, the damping coeff. for T and S
      !!----------------------------------------------------------------------
      USE ioipsl
      !!
      INTEGER                         , INTENT(in   )  ::  kn_file    ! save the damping coef on a file or not
      CHARACTER(len=3)                , INTENT(in   )  ::  cdtype     ! =TRA, TRC or DYN (tracer/dynamics indicator)
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout)  ::  presto     ! restoring coeff. (s-1)
      !
      INTEGER  ::   ji, jj, jk                  ! dummy loop indices
      INTEGER  ::   ii0, ii1, ij0, ij1          ! local integers
      INTEGER  ::   imask        ! File handle 
      INTEGER  ::   jrelax                      ! width of buffer zone
      INTEGER  ::   inum                        ! Logical unit for reading dmp_mask
      REAL(wp) ::   ztrelax, ztvanish           ! restoring time scale
      REAL(wp) ::   zlon1, zlon2                ! Longitude min and max for patch
      REAL(wp) ::   zbw, zd1, zd2               ! Band width, depth limit
      REAL(wp) ::   zv1, zv2                    ! local scalars
      REAL(wp) ::   zinfl, zlon                 ! local scalars
      REAL(wp) ::   zlat, zlat0, zlat1, zlat2   !   -      -
      REAL(wp) ::   zsdmp, zbdmp                !   -      -
      CHARACTER(len=80)                   :: cfile, cvar
      REAL(wp),  DIMENSION(:,:,:), ALLOCATABLE :: zdct 
      !!----------------------------------------------------------------------
      !
      !
      ALLOCATE(zdct(jpi,jpj,jpk))

      presto(:,:,:) = 0._wp
      !
       !--------------------------------------------------------------------------------------------------------------
       ! 3D restoring in semi enclosed seas :
       ! Black Sea :
       ! Lonmin     ! Lonmax       ! Latmin        ! Latmax       ! Band      ! min relax time   !
       !  deg E     ! deg E        !  deg N        ! Deg N        !  deg      !   days           !
       zlon1 = 27.4 ; zlon2 = 42.0 ;  zlat1 = 41.0 ; zlat2 = 47.5 ; zbw = 0.  ; ztrelax = 180. 

       CALL resto_patch ( zlon1, zlon2, zlat1, zlat2, zbw, ztrelax , presto )

       ! Red Sea
       ! Lonmin     ! Lonmax       ! Latmin        ! Latmax       ! Band      ! min relax time   !
       !  deg E     ! deg E        !  deg N        ! Deg N        !  deg      !   days           !
       zlon1 = 29.4 ; zlon2 = 43.6 ;  zlat1 = 12.9 ; zlat2 = 30.3 ; zbw = 0.  ; ztrelax = 180. 

       CALL resto_patch ( zlon1, zlon2, zlat1, zlat2, zbw, ztrelax , presto )

       ! Persian Gulf
       ! Lonmin     ! Lonmax       ! Latmin        ! Latmax       ! Band      ! min relax time   !
       !  deg E     ! deg E        !  deg N        ! Deg N        !  deg      !   days           !
       zlon1 = 46.5 ; zlon2 = 57.0 ;  zlat1 = 23.0 ; zlat2 = 31.5 ; zbw = 1.  ; ztrelax = 180. 

       CALL resto_patch ( zlon1, zlon2, zlat1, zlat2, zbw, ztrelax , presto )
      
       ! ---------------------------------------------------------------------------------------------------------------
       ! 3D overflow restoring 
       ! Gibraltar (limited to 600-1300 m depth range )
       ! Lonmin     ! Lonmax       ! Latmin        ! Latmax       ! Radius    ! relax time   ! dep min    ! dep max    !
       !  deg E     ! deg E        !  deg N        ! Deg N        !  km       !   days       !  m         !   m        !
       zlon1 = -7.  ; zlon2 = -7.  ;  zlat1 = 36.0 ; zlat2 = 36.0 ; zbw = 80. ; ztrelax = 6. ; zd1 = 600. ; zd2 = 1300.

       CALL resto_patch ( zlon1, zlon2, zlat1, zlat2, zbw, ztrelax, presto , zd1, zd2 )

       ! Bab El Mandeb  ( all water column )
       ! Lonmin     ! Lonmax       ! Latmin        ! Latmax       ! Radius    ! relax time   ! dep min    ! dep max    !
       !  deg E     ! deg E        !  deg N        ! Deg N        !  km       !   days       !  m         !   m        !
       zlon1=44.75  ; zlon2 =44.75 ;  zlat1 = 11.5 ; zlat2 = 11.5 ; zbw =100. ; ztrelax = 6. ; zd1 = 0.   ; zd2 = 6000.!

       CALL resto_patch ( zlon1, zlon2, zlat1, zlat2, zbw, ztrelax, presto , zd1, zd2 )

       ! Ormuz Strait  ( all water column )
       ! Lonmin     ! Lonmax       ! Latmin        ! Latmax       ! Radius    ! relax time   ! dep min    ! dep max    !
       !  deg E     ! deg E        !  deg N        ! Deg N        !  km       !   days       !  m         !   m        !
       zlon1=57.75  ; zlon2 =57.75 ;  zlat1 = 25.0 ; zlat2 = 25.0 ; zbw =100. ; ztrelax = 6. ; zd1 = 0.   ; zd2 = 6000.!

       CALL resto_patch ( zlon1, zlon2, zlat1, zlat2, zbw, ztrelax, presto , zd1, zd2 )

       !---------------------------------------------------------------------------------------------------------------
       ! Use a mask_dmp.nc file 
       IF ( ln_dmpmask ) THEN 
         ! in this case, a real 0-1 mask is read into mask_dmp.nc file for 3D restoring following 
         ! particular geometry. This mask is to be build by preprocessing using its own criteria
         ! eg : used for restoring particular water mass. 
         ! The typical restoring time scale is introduced here.
          cfile=TRIM(cn_dir)//'/'//TRIM(sn_dmp%clname)
          cvar = TRIM(sn_dmp%clvar)
          CALL iom_open ( cfile, imask )
          CALL iom_get ( imask, jpdom_data, cvar, zdct )  ! use zdct as temporary array
          CALL iom_close (imask)
          WHERE ( zdct > 1 ) zdct = 0.  !  JMM : WHY ???
          ! it seems that this where is not well accepted on brodie => replaced by a loop
          presto(:,:,:) = presto(:,:,:) + zdct(:,:,:)/rn_timsk/86400.
           IF (lwp) WRITE(numout,*) 'dtacof : read dmp_mask.nc file '
           IF (lwp) WRITE(numout,*) '~~~~~'
       ENDIF
       ! Particular cases
       IF ( cn_cfg == 'natl' ) THEN  ! ALL NATL config have an eastern  Med sea truncated
         ! Eastern Med Sea truncated in NATL
         ! Lonmin     ! Lonmax       ! Latmin        ! Latmax       ! Band      ! relax time   !
         !  deg E     ! deg E        !  deg N        ! Deg N        !  deg      !   days       !
         zlon1= 20.0  ; zlon2 = 26.0 ;  zlat1 = 28.0 ; zlat2 = 45.0 ; zbw =2.   ; ztrelax = 3. 

         CALL resto_patch ( zlon1, zlon2, zlat1, zlat2, zbw, ztrelax , presto )
        
!         IF ( Agrif_Root() ) THEN
!         ! Additional modification in case of closed boundary (north and/or south )
!         IF ( .NOT. lk_obc ) THEN
!           jrelax = 28 ! number of point for the buffer zone
!                       ! use resto_patch with lon lat expressed in i,j points, but set ld_ij flag
!           ztrelax  = 3.
!           ztvanish= 100.
!           zv1 = ( ztvanish -  ztrelax ) / ( jrelax  - 1 )
!           DO jj=1, jrelax
!              zv2 = 1./ ( ztrelax + ( jj - 1 ) * zv1 ) / 86400.
!              ! South
!              presto (:, mj0(jj):mj1(jj)                , :) = zv2
!              ! North
!              presto (:, mj0(jpjglo-jj):mj1(jpjglo - jj), :) = zv2
!           ENDDO
!         ENDIF
!         ENDIF  ! agrif root
       ENDIF

      !                            !--------------------------------!
      IF( kn_file == 1 ) THEN      !  save damping coef. in a file  !
         !                         !--------------------------------!
         IF(lwp) WRITE(numout,*) '              create damping.coeff.nc file'
         IF( cdtype == 'TRA' ) cfile = 'damping.coeff'
         IF( cdtype == 'TRC' ) cfile = 'damping.coeff.trc'
         IF( cdtype == 'DYN' ) cfile = 'damping.coeff.dyn'
         cfile = TRIM( cfile )
         CALL iom_open  ( cfile, imask, ldwrt = .TRUE. )
         CALL iom_rstput( 0, 0, imask, 'Resto', presto )
         CALL iom_close ( imask )
      ENDIF
      !
      DEALLOCATE (zdct )
      !
      !
   END SUBROUTINE dtacof

   SUBROUTINE resto_patch ( plon1, plon2, plat1, plat2, pbw ,ptmax, presto, pz1, pz2 )
      !!------------------------------------------------------------------------
      !!                 ***  Routine resto_patch  ***
      !!
      !! ** Purpose :   modify resto array on a geographically defined zone.
      !!
      !! ** Method  :  Use glamt, gphit arrays. If the defined zone is outside 
      !!              the domain, resto is unchanged. If pz1 and pz2 are provided
      !!              then plon1, plat1 is taken as the position of a the center
      !!              of a circle with decay radius is pbw (in km) 
      !!
      !! ** Action  : IF not present pz1, pz2 : 
      !!              - plon1, plon2 : min and max longitude of the zone (Deg)
      !!              - plat1, plat2 : min and max latitude of the zone (Deg)
      !!              - pbw : band width of the linear decaying restoring (Deg)
      !!              - ptmax : restoring time scale for the inner zone (days)
      !!              IF present pz1 pz2
      !!              - plon1, plat1 : position of the center of the circle
      !!              - pbw = radius (km) of the restoring circle
      !!              - ptmax = time scale at maximum restoring
      !!              - pz1, pz2 : optional: if used, define the depth range (m)
      !!                          for restoring. If not all depths are considered
      !!------------------------------------------------------------------------
      REAL(wp),                   INTENT(in   ) :: plon1, plon2, plat1, plat2, pbw, ptmax
      REAL(wp), DIMENSION(:,:,:), INTENT(inout) :: presto 
      REAL(wp), OPTIONAL, INTENT(in) :: pz1, pz2
      !!
      INTEGER :: ji,jj, jk    ! dummy loop index
      INTEGER :: ik1, ik2     ! limiting vertical index corresponding to pz1,pz2
      INTEGER :: ij0, ij1, iiO, ii1
      INTEGER, DIMENSION(1)           :: iloc 

      REAL(wp) :: zv1, zv2, zv3, zv4, zcoef, ztmp, zdist, zradius2, zcoef2
      REAL(wp), DIMENSION(:,:), ALLOCATABLE :: zpatch
      REAL(wp), DIMENSION(:),   ALLOCATABLE :: zmask
      !!------------------------------------------------------------------------
      ALLOCATE(zpatch(jpi,jpj),  zmask(jpk))
 
      zpatch = 0._wp
      zcoef  = 1._wp/ptmax/86400._wp

      IF (PRESENT (pz1) ) THEN 
        ! horizontal extent
        zradius2 = pbw * pbw !  radius squared
        DO jj = 1, jpj
          DO ji = 1 , jpi
            zpatch(ji,jj) =  sin(gphit(ji,jj)*rad)* sin(plat1*rad)  &
       &                         + cos(gphit(ji,jj)*rad)* cos(plat1*rad)  &
       &                         * cos(rad*(plon1-glamt(ji,jj)))
          ENDDO
        ENDDO

        WHERE ( abs (zpatch ) > 1 ) zpatch = 1.
        DO jj = 1, jpj
          DO ji= 1, jpi 
             ztmp = zpatch(ji,jj)
             zdist = atan(sqrt( (1.-ztmp)/(1+ztmp)) )*2.*ra/1000.
             zpatch(ji,jj) = exp( - zdist*zdist/zradius2 )
          ENDDO
        ENDDO
        ! clean cut off
        WHERE (ABS(zpatch) < 0.01 ) zpatch = 0.
        ! Vertical limitation
        zmask = 1.
        WHERE ( gdept_1d < pz1 .OR. gdept_1d > pz2 ) zmask = 0.
        ! Look for first 1 from the top and keep corresponding ik1 index
        iloc=MAXLOC(zmask) ; ik1 = iloc(1)
        ! Reset the upper part of zmask to 1
        zmask(1:ik1) = 1.
        ! Look from the first 0 from the top and infer the ik2 index
        iloc=MINLOC(zmask) ; ik2 = iloc(1) - 1
        ! Now restore zmask to correct value using a 2 points ramp if possible
        zmask =0._wp
        IF (ik2 -ik1 > 4 ) THEN
          zmask(ik1       ) = 0.25_wp
          zmask(ik1+1     ) = 0.75_wp
          zmask(ik1+2:ik2-2) = 1.0_wp
          zmask(ik2-1     ) = 0.75_wp
          zmask(ik2       ) = 0.25_wp
        ELSE
          zmask(ik1:ik2) = 1._wp   ! all the water column is restored the same
        ENDIF

        DO jk=1, jpk
          presto(:,:,jk)= presto(:,:,jk) + zpatch * zcoef * zmask(jk) 
        ENDDO
        ! JMM : eventually add some checking to avoid locally large resto.

      ELSE
        ! horizontal extent
        zcoef2=1./(pbw +1.e-20 ) ! to avoid division by 0
        DO jj=1,jpj
          DO ji=1,jpi
             zv1=MAX(0., zcoef2*( glamt(ji,jj) - plon1)  )
             zv2=MAX(0., zcoef2*( plon2 - glamt(ji,jj))  )
             zv3=MAX(0., zcoef2*( gphit(ji,jj) - plat1)  )
             zv4=MAX(0., zcoef2*( plat2 - gphit(ji,jj))  )
             zpatch(ji,jj)= MIN( 1., MIN( 1., zv1,zv2,zv3,zv4 ) )
          ENDDO
        ENDDO
          ! resto all over the water column
          presto(:,:,1)= presto(:,:,1) + zpatch *zcoef
          WHERE (zpatch /= 0 ) presto(:,:,1) = MIN( presto(:,:,1), zcoef )
          DO jk=2, jpk
             presto(:,:,jk)=presto(:,:,1)
          ENDDO
      ENDIF
      !
      DEALLOCATE( zpatch, zmask)

      !
   END SUBROUTINE resto_patch
#endif  

   !!======================================================================
END MODULE tradmp
