MODULE sbcblk
   !!======================================================================
   !!                       ***  MODULE  sbcblk  ***
   !! Ocean forcing:  momentum, heat and freshwater flux formulation
   !!                         Aerodynamic Bulk Formulas
   !!                        SUCCESSOR OF "sbcblk_core"
   !!=====================================================================
   !! History :  1.0  !  2004-08  (U. Schweckendiek)  Original CORE code
   !!            2.0  !  2005-04  (L. Brodeau, A.M. Treguier)  improved CORE bulk and its user interface
   !!            3.0  !  2006-06  (G. Madec)  sbc rewritting
   !!             -   !  2006-12  (L. Brodeau)  Original code for turb_core
   !!            3.2  !  2009-04  (B. Lemaire)  Introduce iom_put
   !!            3.3  !  2010-10  (S. Masson)  add diurnal cycle
   !!            3.4  !  2011-11  (C. Harris)  Fill arrays required by CICE
   !!            3.7  !  2014-06  (L. Brodeau)  simplification and optimization of CORE bulk
   !!            4.0  !  2016-06  (L. Brodeau)  sbcblk_core becomes sbcblk and is not restricted to the CORE algorithm anymore
   !!                 !                        ==> based on AeroBulk (https://github.com/brodeau/aerobulk/)
   !!            4.0  !  2016-10  (G. Madec)  introduce a sbc_blk_init routine
   !!            4.0  !  2016-10  (M. Vancoppenolle)  Introduce conduction flux emulator (M. Vancoppenolle)
   !!            4.0  !  2019-03  (F. Lemarié & G. Samson)  add ABL compatibility (ln_abl=TRUE)
   !!            4.2  !  2020-12  (L. Brodeau) Introduction of various air-ice bulk parameterizations + improvements
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   sbc_blk_init  : initialisation of the chosen bulk formulation as ocean surface boundary condition
   !!   sbc_blk       : bulk formulation as ocean surface boundary condition
   !!   blk_oce_1     : computes pieces of momentum, heat and freshwater fluxes over ocean for ABL model  (ln_abl=TRUE)
   !!   blk_oce_2     : finalizes momentum, heat and freshwater fluxes computation over ocean after the ABL step  (ln_abl=TRUE)
   !!             sea-ice case only :
   !!   blk_ice_1   : provide the air-ice stress
   !!   blk_ice_2   : provide the heat and mass fluxes at air-ice interface
   !!   blk_ice_qcn   : provide ice surface temperature and snow/ice conduction flux (emulating conduction flux)
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers
   USE dom_oce        ! ocean space and time domain
   USE phycst         ! physical constants
   USE fldread        ! read input fields
   USE sbc_oce        ! Surface boundary condition: ocean fields
   USE trc_oce         ! share SMS/Ocean variables
   USE cyclone        ! Cyclone 10m wind form trac of cyclone centres
   USE sbcdcy         ! surface boundary condition: diurnal cycle
   USE sbcwave , ONLY :   cdn_wave ! wave module
   USE lib_fortran    ! to use key_nosignedzero and glob_sum
   !
#if defined key_si3
   USE sbc_ice        ! Surface boundary condition: ice fields #LB? ok to be in 'key_si3' ???
   USE ice     , ONLY :   u_ice, v_ice, jpl, a_i_b, at_i_b, t_su, rn_cnd_s, hfx_err_dif, nn_qtrice
   USE icevar         ! for CALL ice_var_snwblow
   USE sbcblk_algo_ice_an05
   USE sbcblk_algo_ice_lu12
   USE sbcblk_algo_ice_lg15
#endif
   USE sbcblk_algo_ncar     ! => turb_ncar     : NCAR - (formerly known as CORE, Large & Yeager, 2009)
   USE sbcblk_algo_coare3p0 ! => turb_coare3p0 : COAREv3.0 (Fairall et al. 2003)
   USE sbcblk_algo_coare3p6 ! => turb_coare3p6 : COAREv3.6 (Fairall et al. 2018 + Edson et al. 2013)
   USE sbcblk_algo_ecmwf    ! => turb_ecmwf    : ECMWF (IFS cycle 45r1)
   USE sbcblk_algo_andreas  ! => turb_andreas  : Andreas et al. 2015
   !
   USE iom            ! I/O manager library
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! distribued memory computing library
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE prtctl         ! Print control

   USE sbc_phy        ! Catalog of functions for physical/meteorological parameters in the marine boundary layer

   IMPLICIT NONE
   PRIVATE

   PUBLIC   sbc_blk_init  ! called in sbcmod
   PUBLIC   sbc_blk       ! called in sbcmod
   PUBLIC   blk_oce_1     ! called in sbcabl
   PUBLIC   blk_oce_2     ! called in sbcabl
#if defined key_si3
   PUBLIC   blk_ice_1     ! routine called in icesbc
   PUBLIC   blk_ice_2     ! routine called in icesbc
   PUBLIC   blk_ice_qcn   ! routine called in icesbc
#endif

   INTEGER , PUBLIC, PARAMETER ::   jp_wndi  =  1   ! index of 10m wind velocity (i-component) (m/s)    at T-point
   INTEGER , PUBLIC, PARAMETER ::   jp_wndj  =  2   ! index of 10m wind velocity (j-component) (m/s)    at T-point
   INTEGER , PUBLIC, PARAMETER ::   jp_tair  =  3   ! index of 10m air temperature             (Kelvin)
   INTEGER , PUBLIC, PARAMETER ::   jp_humi  =  4   ! index of specific humidity               (kg/kg)
   INTEGER , PUBLIC, PARAMETER ::   jp_qsr   =  5   ! index of solar heat                      (W/m2)
   INTEGER , PUBLIC, PARAMETER ::   jp_qlw   =  6   ! index of Long wave                       (W/m2)
   INTEGER , PUBLIC, PARAMETER ::   jp_prec  =  7   ! index of total precipitation (rain+snow) (Kg/m2/s)
   INTEGER , PUBLIC, PARAMETER ::   jp_snow  =  8   ! index of snow (solid prcipitation)       (kg/m2/s)
   INTEGER , PUBLIC, PARAMETER ::   jp_slp   =  9   ! index of sea level pressure              (Pa)
   INTEGER , PUBLIC, PARAMETER ::   jp_uoatm = 10   ! index of surface current (i-component)
   !                                                !          seen by the atmospheric forcing (m/s) at T-point
   INTEGER , PUBLIC, PARAMETER ::   jp_voatm = 11   ! index of surface current (j-component)
   !                                                !          seen by the atmospheric forcing (m/s) at T-point
   INTEGER , PUBLIC, PARAMETER ::   jp_cc    = 12   ! index of cloud cover                     (-)      range:0-1
   INTEGER , PUBLIC, PARAMETER ::   jp_hpgi  = 13   ! index of ABL geostrophic wind or hpg (i-component) (m/s) at T-point
   INTEGER , PUBLIC, PARAMETER ::   jp_hpgj  = 14   ! index of ABL geostrophic wind or hpg (j-component) (m/s) at T-point
#if defined key_drakkar
   LOGICAL                     ::   ln_clim_forcing = .false.   ! Logical flag for use of climatological forcing 
                                                  !     ( T= climatological, F=interannual)
   INTEGER , PUBLIC, PARAMETER ::   jp_wmod =15     ! index of wind module                     (m/s)
   INTEGER , PUBLIC, PARAMETER ::   jp_uw   =16     ! index of zonal pseudo stress             (m2/s2)
   INTEGER , PUBLIC, PARAMETER ::   jp_vw   =17     ! index of meridional pseudo stress        (m2/s2)
   INTEGER , PUBLIC, PARAMETER ::   jpfld    = 17   ! maximum number of files to read
#else
   INTEGER , PUBLIC, PARAMETER ::   jpfld    = 14   ! maximum number of files to read
#endif

   ! Warning: keep this structure allocatable for Agrif...
   TYPE(FLD), PUBLIC, ALLOCATABLE, DIMENSION(:) ::   sf   ! structure of input atmospheric fields (file informations, fields read)

   !                           !!* Namelist namsbc_blk : bulk parameters
   LOGICAL  ::   ln_NCAR        ! "NCAR"      algorithm   (Large and Yeager 2008)
   LOGICAL  ::   ln_COARE_3p0   ! "COARE 3.0" algorithm   (Fairall et al. 2003)
   LOGICAL  ::   ln_COARE_3p6   ! "COARE 3.6" algorithm   (Edson et al. 2013)
   LOGICAL  ::   ln_ECMWF       ! "ECMWF"     algorithm   (IFS cycle 45r1)
   LOGICAL  ::   ln_ANDREAS     ! "ANDREAS"   algorithm   (Andreas et al. 2015)
   !
   !#LB:
   LOGICAL  ::   ln_Cx_ice_cst             ! use constant air-ice bulk transfer coefficients (value given in namelist's rn_Cd_i, rn_Ce_i & rn_Ch_i)
   REAL(wp) ::   rn_Cd_i, rn_Ce_i, rn_Ch_i ! values for  "    "
   LOGICAL  ::   ln_Cx_ice_AN05            ! air-ice bulk transfer coefficients based on Andreas et al., 2005
   LOGICAL  ::   ln_Cx_ice_LU12            ! air-ice bulk transfer coefficients based on Lupkes et al., 2012
   LOGICAL  ::   ln_Cx_ice_LG15            ! air-ice bulk transfer coefficients based on Lupkes & Gryanik, 2015
   !#LB.
   !
   LOGICAL  ::   ln_crt_fbk     ! Add surface current feedback to the wind stress computation  (Renault et al. 2020)
   REAL(wp) ::   rn_stau_a      ! Alpha and Beta coefficients of Renault et al. 2020, eq. 10: Stau = Alpha * Wnd + Beta
   REAL(wp) ::   rn_stau_b      !
   !
   REAL(wp)         ::   rn_pfac   ! multiplication factor for precipitation
   REAL(wp), PUBLIC ::   rn_efac   ! multiplication factor for evaporation
   REAL(wp)         ::   rn_zqt    ! z(q,t) : height of humidity and temperature measurements
   REAL(wp)         ::   rn_zu     ! z(u)   : height of wind measurements
   !
   INTEGER          :: nn_iter_algo   !  Number of iterations in bulk param. algo ("stable ABL + weak wind" requires more)

   REAL(wp),         ALLOCATABLE, DIMENSION(:,:) ::   theta_zu, q_zu                 ! air temp. and spec. hum. at wind speed height (L15 bulk scheme)

#if defined key_si3
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: Cd_ice , Ch_ice , Ce_ice   !#LB transfert coefficients over ice
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: theta_zu_i, q_zu_i         !#LB fixme ! air temp. and spec. hum. over ice at wind speed height (L15 bulk scheme)
#endif


   LOGICAL  ::   ln_skin_cs     ! use the cool-skin (only available in ECMWF and COARE algorithms) !LB
   LOGICAL  ::   ln_skin_wl     ! use the warm-layer parameterization (only available in ECMWF and COARE algorithms) !LB
   LOGICAL  ::   ln_humi_sph    ! humidity read in files ("sn_humi") is specific humidity [kg/kg] if .true. !LB
   LOGICAL  ::   ln_humi_dpt    ! humidity read in files ("sn_humi") is dew-point temperature [K] if .true. !LB
   LOGICAL  ::   ln_humi_rlh    ! humidity read in files ("sn_humi") is relative humidity     [%] if .true. !LB
   LOGICAL  ::   ln_tair_pot    ! temperature read in files ("sn_tair") is already potential temperature (not absolute)
   !
   INTEGER  ::   nhumi          ! choice of the bulk algorithm
   !                            ! associated indices:
   INTEGER, PARAMETER :: np_humi_sph = 1
   INTEGER, PARAMETER :: np_humi_dpt = 2
   INTEGER, PARAMETER :: np_humi_rlh = 3

   INTEGER  ::   nblk           ! choice of the bulk algorithm
   !                            ! associated indices:
   INTEGER, PARAMETER ::   np_NCAR      = 1   ! "NCAR" algorithm        (Large and Yeager 2008)
   INTEGER, PARAMETER ::   np_COARE_3p0 = 2   ! "COARE 3.0" algorithm   (Fairall et al. 2003)
   INTEGER, PARAMETER ::   np_COARE_3p6 = 3   ! "COARE 3.6" algorithm   (Edson et al. 2013)
   INTEGER, PARAMETER ::   np_ECMWF     = 4   ! "ECMWF" algorithm       (IFS cycle 45r1)
   INTEGER, PARAMETER ::   np_ANDREAS   = 5   ! "ANDREAS" algorithm       (Andreas et al. 2015)

   !#LB:
#if defined key_si3
   ! Same, over sea-ice:
   INTEGER  ::   nblk_ice           ! choice of the bulk algorithm
   !                            ! associated indices:
   INTEGER, PARAMETER ::   np_ice_cst  = 1   ! constant transfer coefficients
   INTEGER, PARAMETER ::   np_ice_an05 = 2   ! Andreas et al., 2005
   INTEGER, PARAMETER ::   np_ice_lu12 = 3   ! Lupkes el al., 2012
   INTEGER, PARAMETER ::   np_ice_lg15 = 4   ! Lupkes & Gryanik, 2015
#endif
   !LB.



   !! * Substitutions
#  include "do_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: sbcblk.F90 15551 2021-11-28 20:19:36Z gsamson $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION sbc_blk_alloc()
      !!-------------------------------------------------------------------
      !!             ***  ROUTINE sbc_blk_alloc ***
      !!-------------------------------------------------------------------
      ALLOCATE( theta_zu(jpi,jpj), q_zu(jpi,jpj),  STAT=sbc_blk_alloc )
      CALL mpp_sum ( 'sbcblk', sbc_blk_alloc )
      IF( sbc_blk_alloc /= 0 )   CALL ctl_stop( 'STOP', 'sbc_blk_alloc: failed to allocate arrays' )
   END FUNCTION sbc_blk_alloc

#if defined key_si3
   INTEGER FUNCTION sbc_blk_ice_alloc()
      !!-------------------------------------------------------------------
      !!             ***  ROUTINE sbc_blk_ice_alloc ***
      !!-------------------------------------------------------------------
      ALLOCATE( Cd_ice (jpi,jpj), Ch_ice (jpi,jpj), Ce_ice (jpi,jpj),         &
         &      theta_zu_i(jpi,jpj), q_zu_i(jpi,jpj),  STAT=sbc_blk_ice_alloc )
      CALL mpp_sum ( 'sbcblk', sbc_blk_ice_alloc )
      IF( sbc_blk_ice_alloc /= 0 )   CALL ctl_stop( 'STOP', 'sbc_blk_ice_alloc: failed to allocate arrays' )
   END FUNCTION sbc_blk_ice_alloc
#endif


   SUBROUTINE sbc_blk_init
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE sbc_blk_init  ***
      !!
      !! ** Purpose :   choose and initialize a bulk formulae formulation
      !!
      !! ** Method  :
      !!
      !!----------------------------------------------------------------------
      INTEGER  ::   jfpr                  ! dummy loop indice and argument
      INTEGER  ::   ios, ierror, ioptio   ! Local integer
      !!
      CHARACTER(len=100)            ::   cn_dir                ! Root directory for location of atmospheric forcing files
      TYPE(FLD_N), DIMENSION(jpfld) ::   slf_i                 ! array of namelist informations on the fields to read
      TYPE(FLD_N) ::   sn_wndi, sn_wndj , sn_humi, sn_qsr      ! informations about the fields to be read
      TYPE(FLD_N) ::   sn_qlw , sn_tair , sn_prec, sn_snow     !       "                        "
      TYPE(FLD_N) ::   sn_slp , sn_uoatm, sn_voatm             !       "                        "
      TYPE(FLD_N) ::   sn_cc, sn_hpgi, sn_hpgj                 !       "                        "
      INTEGER     ::   ipka                                    ! number of levels in the atmospheric variable
      NAMELIST/namsbc_blk/ ln_NCAR, ln_COARE_3p0, ln_COARE_3p6, ln_ECMWF, ln_ANDREAS, &   ! bulk algorithm
         &                 rn_zqt, rn_zu, nn_iter_algo, ln_skin_cs, ln_skin_wl,       &
         &                 rn_pfac, rn_efac,                                          &
         &                 ln_crt_fbk, rn_stau_a, rn_stau_b,                          &   ! current feedback
         &                 ln_humi_sph, ln_humi_dpt, ln_humi_rlh, ln_tair_pot,        &
         &                 ln_Cx_ice_cst, rn_Cd_i, rn_Ce_i, rn_Ch_i,                  &
         &                 ln_Cx_ice_AN05, ln_Cx_ice_LU12, ln_Cx_ice_LG15,            &
         &                 cn_dir,                                                    &
         &                 sn_wndi, sn_wndj, sn_qsr, sn_qlw ,                         &   ! input fields
         &                 sn_tair, sn_humi, sn_prec, sn_snow, sn_slp,                &
         &                 sn_uoatm, sn_voatm, sn_cc, sn_hpgi, sn_hpgj
#if defined key_drakkar
      TYPE(FLD_N) :: sn_wmod, sn_uw, sn_vw
      NAMELIST/namsbc_blk_drk/ ln_clim_forcing,sn_wmod, sn_uw, sn_vw
#endif

      ! cool-skin / warm-layer !LB
      !!---------------------------------------------------------------------
      !
      !                                      ! allocate sbc_blk_core array
      IF( sbc_blk_alloc()     /= 0 )   CALL ctl_stop( 'STOP', 'sbc_blk : unable to allocate standard arrays' )
      !
#if defined key_si3
      IF( sbc_blk_ice_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'sbc_blk : unable to allocate standard ice arrays' )
#endif
      !
      !                             !** read bulk namelist
      READ  ( numnam_ref, namsbc_blk, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namsbc_blk in reference namelist' )
      !
      READ  ( numnam_cfg, namsbc_blk, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namsbc_blk in configuration namelist' )
      !
      IF(lwm) WRITE( numond, namsbc_blk )
#if defined key_drakkar
      !                             !** read bulk  drakkar namelist  
      READ  ( numnam_ref, namsbc_blk_drk, IOSTAT = ios, ERR = 903)
903   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namsbc_blk_drk in reference namelist' )
      !
      READ  ( numnam_cfg, namsbc_blk_drk, IOSTAT = ios, ERR = 904 )
904   IF( ios >  0 )   CALL ctl_nam ( ios , 'namsbc_blk_drk in configuration namelist' )
      !
      IF(lwm) WRITE( numond, namsbc_blk_drk )
#endif
      !
      !                             !** initialization of the chosen bulk formulae (+ check)
      !                                   !* select the bulk chosen in the namelist and check the choice
      ioptio = 0
      IF( ln_NCAR      ) THEN
         nblk =  np_NCAR        ;   ioptio = ioptio + 1
      ENDIF
      IF( ln_COARE_3p0 ) THEN
         nblk =  np_COARE_3p0   ;   ioptio = ioptio + 1
      ENDIF
      IF( ln_COARE_3p6 ) THEN
         nblk =  np_COARE_3p6   ;   ioptio = ioptio + 1
      ENDIF
      IF( ln_ECMWF     ) THEN
         nblk =  np_ECMWF       ;   ioptio = ioptio + 1
      ENDIF
      IF( ln_ANDREAS     ) THEN
         nblk =  np_ANDREAS       ;   ioptio = ioptio + 1
      ENDIF
      IF( ioptio /= 1 )   CALL ctl_stop( 'sbc_blk_init: Choose one and only one bulk algorithm' )

      !                             !** initialization of the cool-skin / warm-layer parametrization
      IF( ln_skin_cs .OR. ln_skin_wl ) THEN
         !! Some namelist sanity tests:
         IF( ln_NCAR )      &
            & CALL ctl_stop( 'sbc_blk_init: Cool-skin/warm-layer param. not compatible with NCAR algorithm' )
         IF( ln_ANDREAS )      &
            & CALL ctl_stop( 'sbc_blk_init: Cool-skin/warm-layer param. not compatible with ANDREAS algorithm' )
         !IF( nn_fsbc /= 1 ) &
         !   & CALL ctl_stop( 'sbc_blk_init: Please set "nn_fsbc" to 1 when using cool-skin/warm-layer param.')
      END IF

      IF( ln_skin_wl ) THEN
         !! Check if the frequency of downwelling solar flux input makes sense and if ln_dm2dc=T if it is daily!
         IF( (sn_qsr%freqh  < 0.).OR.(sn_qsr%freqh  > 24.) ) &
            & CALL ctl_stop( 'sbc_blk_init: Warm-layer param. (ln_skin_wl) not compatible with freq. of solar flux > daily' )
         IF( (sn_qsr%freqh == 24.).AND.(.NOT. ln_dm2dc) ) &
            & CALL ctl_stop( 'sbc_blk_init: Please set ln_dm2dc=T for warm-layer param. (ln_skin_wl) to work properly' )
      END IF

      ioptio = 0
      IF( ln_humi_sph ) THEN
         nhumi =  np_humi_sph    ;   ioptio = ioptio + 1
      ENDIF
      IF( ln_humi_dpt ) THEN
         nhumi =  np_humi_dpt    ;   ioptio = ioptio + 1
      ENDIF
      IF( ln_humi_rlh ) THEN
         nhumi =  np_humi_rlh    ;   ioptio = ioptio + 1
      ENDIF
      IF( ioptio /= 1 )   CALL ctl_stop( 'sbc_blk_init: Choose one and only one type of air humidity' )
      !
      IF( ln_dm2dc ) THEN                 !* check: diurnal cycle on Qsr
         IF( sn_qsr%freqh /= 24. )   CALL ctl_stop( 'sbc_blk_init: ln_dm2dc=T only with daily short-wave input' )
         IF( sn_qsr%ln_tint ) THEN
            CALL ctl_warn( 'sbc_blk_init: ln_dm2dc=T daily qsr time interpolation done by sbcdcy module',   &
               &           '              ==> We force time interpolation = .false. for qsr' )
            sn_qsr%ln_tint = .false.
         ENDIF
      ENDIF

#if defined key_si3
      ioptio = 0
      IF( ln_Cx_ice_cst ) THEN
         nblk_ice =  np_ice_cst     ;   ioptio = ioptio + 1
      ENDIF
      IF( ln_Cx_ice_AN05 ) THEN
         nblk_ice =  np_ice_an05   ;   ioptio = ioptio + 1
      ENDIF
      IF( ln_Cx_ice_LU12 ) THEN
         nblk_ice =  np_ice_lu12    ;   ioptio = ioptio + 1
      ENDIF
      IF( ln_Cx_ice_LG15 ) THEN
         nblk_ice =  np_ice_lg15   ;   ioptio = ioptio + 1
      ENDIF
      IF( ioptio /= 1 )   CALL ctl_stop( 'sbc_blk_init: Choose one and only one ice-atm bulk algorithm' )
#endif


      !                                   !* set the bulk structure
      !                                      !- store namelist information in an array
      !
      slf_i(jp_wndi ) = sn_wndi    ;   slf_i(jp_wndj ) = sn_wndj
      slf_i(jp_qsr  ) = sn_qsr     ;   slf_i(jp_qlw  ) = sn_qlw
      slf_i(jp_tair ) = sn_tair    ;   slf_i(jp_humi ) = sn_humi
      slf_i(jp_prec ) = sn_prec    ;   slf_i(jp_snow ) = sn_snow
      slf_i(jp_slp  ) = sn_slp     ;   slf_i(jp_cc   ) = sn_cc
      slf_i(jp_uoatm) = sn_uoatm   ;   slf_i(jp_voatm) = sn_voatm
      slf_i(jp_hpgi ) = sn_hpgi    ;   slf_i(jp_hpgj ) = sn_hpgj
      !
      IF( .NOT. ln_abl ) THEN   ! force to not use jp_hpgi and jp_hpgj, should already be done in namelist_* but we never know...
         slf_i(jp_hpgi)%clname = 'NOT USED'
         slf_i(jp_hpgj)%clname = 'NOT USED'
      ENDIF
#if defined key_drakkar
      slf_i(jp_wmod ) = sn_wmod    ;   slf_i(jp_uw ) = sn_uw  ; slf_i(jp_vw ) = sn_vw
      IF( .NOT. ln_clim_forcing) THEN 
         slf_i(jp_wmod)%clname = 'NOT USED'
         slf_i(jp_uw  )%clname = 'NOT USED'
         slf_i(jp_vw  )%clname = 'NOT USED'
      ELSE
         slf_i(jp_wndi)%clname = 'NOT USED'
         slf_i(jp_wndj)%clname = 'NOT USED'
      ENDIF
#endif
      !
      !                                      !- allocate the bulk structure
      ALLOCATE( sf(jpfld), STAT=ierror )
      IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'sbc_blk_init: unable to allocate sf structure' )
      !
      !                                      !- fill the bulk structure with namelist informations
      CALL fld_fill( sf, slf_i, cn_dir, 'sbc_blk_init', 'surface boundary condition -- bulk formulae', 'namsbc_blk' )
      sf(jp_wndi )%zsgn = -1._wp   ;   sf(jp_wndj )%zsgn = -1._wp   ! vector field at T point: overwrite default definition of zsgn
      sf(jp_uoatm)%zsgn = -1._wp   ;   sf(jp_voatm)%zsgn = -1._wp   ! vector field at T point: overwrite default definition of zsgn
      sf(jp_hpgi )%zsgn = -1._wp   ;   sf(jp_hpgj )%zsgn = -1._wp   ! vector field at T point: overwrite default definition of zsgn
      !
      DO jfpr= 1, jpfld
         !
         IF(   ln_abl    .AND.                                                      &
            &    ( jfpr == jp_wndi .OR. jfpr == jp_wndj .OR. jfpr == jp_humi .OR.   &
            &      jfpr == jp_hpgi .OR. jfpr == jp_hpgj .OR. jfpr == jp_tair     )  ) THEN
            ipka = jpka   ! ABL: some fields are 3D input
         ELSE
            ipka = 1
         ENDIF
         !
         ALLOCATE( sf(jfpr)%fnow(jpi,jpj,ipka) )
         !
         IF( TRIM(sf(jfpr)%clrootname) == 'NOT USED' ) THEN    !--  not used field  --!   (only now allocated and set to default)
            IF(     jfpr == jp_slp ) THEN
               sf(jfpr)%fnow(:,:,1:ipka) = 101325._wp   ! use standard pressure in Pa
            ELSEIF( jfpr == jp_prec .OR. jfpr == jp_snow .OR. jfpr == jp_uoatm .OR. jfpr == jp_voatm ) THEN
               sf(jfpr)%fnow(:,:,1:ipka) = 0._wp        ! no precip or no snow or no surface currents
            ELSEIF( jfpr == jp_wndi .OR. jfpr == jp_wndj ) THEN
               sf(jfpr)%fnow(:,:,1:ipka) = 0._wp
            ELSEIF( jfpr == jp_hpgi .OR. jfpr == jp_hpgj ) THEN
               IF( .NOT. ln_abl ) THEN
                  DEALLOCATE( sf(jfpr)%fnow )   ! deallocate as not used in this case
               ELSE
                  sf(jfpr)%fnow(:,:,1:ipka) = 0._wp
               ENDIF
            ELSEIF( jfpr == jp_cc  ) THEN
               sf(jp_cc)%fnow(:,:,1:ipka) = pp_cldf
#if defined key_drakkar
            ELSEIF( jfpr == jp_wmod .OR. jfpr == jp_uw .OR. jfpr == jp_vw ) THEN
               sf(jfpr)%fnow(:,:,1:ipka) = 0._wp
#endif
            ELSE
               WRITE(ctmp1,*) 'sbc_blk_init: no default value defined for field number', jfpr
               CALL ctl_stop( ctmp1 )
            ENDIF
         ELSE                                                  !-- used field  --!
            IF( sf(jfpr)%ln_tint )   ALLOCATE( sf(jfpr)%fdta(jpi,jpj,ipka,2) )   ! allocate array for temporal interpolation
            !
            IF( sf(jfpr)%freqh > 0. .AND. MOD( NINT(3600. * sf(jfpr)%freqh), nn_fsbc * NINT(rn_Dt) ) /= 0 )   &
         &  CALL ctl_warn( 'sbc_blk_init: sbcmod timestep rn_Dt*nn_fsbc is NOT a submultiple of atmospheric forcing frequency.',   &
         &                 '               This is not ideal. You should consider changing either rn_Dt or nn_fsbc value...' )
         ENDIF
      END DO
      !
      IF( ln_abl ) THEN       ! ABL: read 3D fields for wind, temperature, humidity and pressure gradient
         rn_zqt = ght_abl(2)          ! set the bulk altitude to ABL first level
         rn_zu  = ght_abl(2)
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) '   ABL formulation: overwrite rn_zqt & rn_zu with ABL first level altitude'
      ENDIF
      !
      !
      IF(lwp) THEN                     !** Control print
         !
         WRITE(numout,*)                  !* namelist
         WRITE(numout,*) '   Namelist namsbc_blk (other than data information):'
         WRITE(numout,*) '      "NCAR"      algorithm   (Large and Yeager 2008)      ln_NCAR      = ', ln_NCAR
         WRITE(numout,*) '      "COARE 3.0" algorithm   (Fairall et al. 2003)       ln_COARE_3p0 = ', ln_COARE_3p0
         WRITE(numout,*) '      "COARE 3.6" algorithm (Fairall 2018 + Edson al 2013) ln_COARE_3p6 = ', ln_COARE_3p6
         WRITE(numout,*) '      "ECMWF"     algorithm   (IFS cycle 45r1)             ln_ECMWF     = ', ln_ECMWF
         WRITE(numout,*) '      "ANDREAS"   algorithm   (Andreas et al. 2015)       ln_ANDREAS   = ', ln_ANDREAS
#if defined key_drakkar
         WRITE(numout,*) '      Use climatolofical forcing formulation           ln_clim_forcing = ', ln_clim_forcing
#endif
         WRITE(numout,*) '      Air temperature and humidity reference height (m)   rn_zqt       = ', rn_zqt
         WRITE(numout,*) '      Wind vector reference height (m)                    rn_zu        = ', rn_zu
         WRITE(numout,*) '      factor applied on precipitation (total & snow)      rn_pfac      = ', rn_pfac
         WRITE(numout,*) '      factor applied on evaporation                       rn_efac      = ', rn_efac
         WRITE(numout,*) '         (form absolute (=0) to relative winds(=1))'
         WRITE(numout,*) '      use surface current feedback on wind stress         ln_crt_fbk   = ', ln_crt_fbk
         IF(ln_crt_fbk) THEN
         WRITE(numout,*) '         Renault et al. 2020, eq. 10: Stau = Alpha * Wnd + Beta'
         WRITE(numout,*) '            Alpha                                         rn_stau_a    = ', rn_stau_a
         WRITE(numout,*) '            Beta                                          rn_stau_b    = ', rn_stau_b
         ENDIF
         !
         WRITE(numout,*)
         SELECT CASE( nblk )              !* Print the choice of bulk algorithm
         CASE( np_NCAR      )   ;   WRITE(numout,*) '   ==>>>   "NCAR" algorithm        (Large and Yeager 2008)'
         CASE( np_COARE_3p0 )   ;   WRITE(numout,*) '   ==>>>   "COARE 3.0" algorithm   (Fairall et al. 2003)'
         CASE( np_COARE_3p6 )   ;   WRITE(numout,*) '   ==>>>   "COARE 3.6" algorithm (Fairall 2018+Edson et al. 2013)'
         CASE( np_ECMWF     )   ;   WRITE(numout,*) '   ==>>>   "ECMWF" algorithm       (IFS cycle 45r1)'
         CASE( np_ANDREAS   )   ;   WRITE(numout,*) '   ==>>>   "ANDREAS" algorithm (Andreas et al. 2015)'
         END SELECT
         !
         WRITE(numout,*)
         WRITE(numout,*) '      use cool-skin  parameterization (SSST)  ln_skin_cs  = ', ln_skin_cs
         WRITE(numout,*) '      use warm-layer parameterization (SSST)  ln_skin_wl  = ', ln_skin_wl
         !
         WRITE(numout,*)
         SELECT CASE( nhumi )              !* Print the choice of air humidity
         CASE( np_humi_sph )   ;   WRITE(numout,*) '   ==>>>   air humidity is SPECIFIC HUMIDITY     [kg/kg]'
         CASE( np_humi_dpt )   ;   WRITE(numout,*) '   ==>>>   air humidity is DEW-POINT TEMPERATURE [K]'
         CASE( np_humi_rlh )   ;   WRITE(numout,*) '   ==>>>   air humidity is RELATIVE HUMIDITY     [%]'
         END SELECT
         !
         !#LB:
#if defined key_si3
         IF( nn_ice > 0 ) THEN
            WRITE(numout,*)
            WRITE(numout,*) '      use constant ice-atm bulk transfer coeff.           ln_Cx_ice_cst  = ', ln_Cx_ice_cst
            WRITE(numout,*) '      use ice-atm bulk coeff. from Andreas et al., 2005   ln_Cx_ice_AN05 = ', ln_Cx_ice_AN05
            WRITE(numout,*) '      use ice-atm bulk coeff. from Lupkes et al., 2012    ln_Cx_ice_LU12 = ', ln_Cx_ice_LU12
            WRITE(numout,*) '      use ice-atm bulk coeff. from Lupkes & Gryanik, 2015 ln_Cx_ice_LG15 = ', ln_Cx_ice_LG15
         ENDIF
         WRITE(numout,*)
         SELECT CASE( nblk_ice )              !* Print the choice of bulk algorithm
         CASE( np_ice_cst  )
            WRITE(numout,*) '   ==>>>   Constant bulk transfer coefficients over sea-ice:'
            WRITE(numout,*) '      => Cd_ice, Ce_ice, Ch_ice =', REAL(rn_Cd_i,4), REAL(rn_Ce_i,4), REAL(rn_Ch_i,4)
            IF( (rn_Cd_i<0._wp).OR.(rn_Cd_i>1.E-2_wp).OR.(rn_Ce_i<0._wp).OR.(rn_Ce_i>1.E-2_wp).OR.(rn_Ch_i<0._wp).OR.(rn_Ch_i>1.E-2_wp) ) &
               & CALL ctl_stop( 'Be realistic in your pick of Cd_ice, Ce_ice & Ch_ice ! (0 < Cx < 1.E-2)')
         CASE( np_ice_an05 )   ;   WRITE(numout,*) '   ==>>> bulk algo over ice: Andreas et al, 2005'
         CASE( np_ice_lu12 )   ;   WRITE(numout,*) '   ==>>> bulk algo over ice: Lupkes et al, 2012'
         CASE( np_ice_lg15 )   ;   WRITE(numout,*) '   ==>>> bulk algo over ice: Lupkes & Gryanik, 2015'
         END SELECT
#endif
         !#LB.
         !
      ENDIF
      !
   END SUBROUTINE sbc_blk_init


   SUBROUTINE sbc_blk( kt )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE sbc_blk  ***
      !!
      !! ** Purpose :   provide at each time step the surface ocean fluxes
      !!              (momentum, heat, freshwater and runoff)
      !!
      !! ** Method  :
      !!              (1) READ each fluxes in NetCDF files:
      !!      the wind velocity (i-component) at z=rn_zu  (m/s) at T-point
      !!      the wind velocity (j-component) at z=rn_zu  (m/s) at T-point
      !!      the specific humidity           at z=rn_zqt (kg/kg)
      !!      the air temperature             at z=rn_zqt (Kelvin)
      !!      the solar heat                              (W/m2)
      !!      the Long wave                               (W/m2)
      !!      the total precipitation (rain+snow)         (Kg/m2/s)
      !!      the snow (solid precipitation)              (kg/m2/s)
      !!      ABL dynamical forcing (i/j-components of either hpg or geostrophic winds)
      !!              (2) CALL blk_oce_1 and blk_oce_2
      !!
      !!      C A U T I O N : never mask the surface stress fields
      !!                      the stress is assumed to be in the (i,j) mesh referential
      !!
      !! ** Action  :   defined at each time-step at the air-sea interface
      !!              - utau, vtau  i- and j-component of the wind stress
      !!              - taum        wind stress module at T-point
      !!              - wndm        wind speed  module at T-point over free ocean or leads in presence of sea-ice
      !!              - qns, qsr    non-solar and solar heat fluxes
      !!              - emp         upward mass flux (evapo. - precip.)
      !!              - sfx         salt flux due to freezing/melting (non-zero only if ice is present)
      !!
      !! ** References :   Large & Yeager, 2004 / Large & Yeager, 2008
      !!                   Brodeau et al. Ocean Modelling 2010
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time step
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj) ::   zssq, zcd_du, zsen, zlat, zevp, zpre, ztheta
      REAL(wp) :: ztst
      LOGICAL  :: llerr
      !!----------------------------------------------------------------------
      !
      CALL fld_read( kt, nn_fsbc, sf )             ! input fields provided at the current time-step

      ! Sanity/consistence test on humidity at first time step to detect potential screw-up:
      IF( kt == nit000 ) THEN
         ! mean humidity over ocean on proc
         ztst = glob_sum( 'sbcblk', sf(jp_humi)%fnow(:,:,1) * e1e2t(:,:) * tmask(:,:,1) ) / glob_sum( 'sbcblk', e1e2t(:,:) * tmask(:,:,1) )
         llerr = .FALSE.
         SELECT CASE( nhumi )
         CASE( np_humi_sph ) ! specific humidity => expect: 0. <= something < 0.065 [kg/kg] (0.061 is saturation at 45degC !!!)
            IF( (ztst <   0._wp) .OR. (ztst > 0.065_wp) )   llerr = .TRUE.
         CASE( np_humi_dpt ) ! dew-point temperature => expect: 110. <= something < 320. [K]
            IF( (ztst < 110._wp) .OR. (ztst >  320._wp) )   llerr = .TRUE.
         CASE( np_humi_rlh ) ! relative humidity => expect: 0. <= something < 100. [%]
            IF( (ztst <   0._wp) .OR. (ztst >  100._wp) )   llerr = .TRUE.
         END SELECT
         IF(llerr) THEN
            WRITE(ctmp1,'("   Error on mean humidity value: ",f10.5)') ztst
            CALL ctl_stop( 'STOP', ctmp1, 'Something is wrong with air humidity!!!', &
               &   ' ==> check the unit in your input files'       , &
               &   ' ==> check consistence of namelist choice: specific? relative? dew-point?', &
               &   ' ==> ln_humi_sph -> [kg/kg] | ln_humi_rlh -> [%] | ln_humi_dpt -> [K] !!!' )
         ENDIF
         IF(lwp) THEN
            WRITE(numout,*) ''
            WRITE(numout,*) ' Global mean humidity at kt = nit000: ', ztst
            WRITE(numout,*) ' === Sanity/consistence test on air humidity sucessfuly passed! ==='
            WRITE(numout,*) ''
         ENDIF
      ENDIF   !IF( kt == nit000 )
      !                                            ! compute the surface ocean fluxes using bulk formulea
      IF( MOD( kt - 1, nn_fsbc ) == 0 ) THEN

         ! Specific humidity of air at z=rn_zqt
         SELECT CASE( nhumi )
         CASE( np_humi_sph )
            q_air_zt(:,:) = sf(jp_humi )%fnow(:,:,1)      ! what we read in file is already a spec. humidity!
         CASE( np_humi_dpt )
            IF((kt==nit000).AND.lwp) WRITE(numout,*) ' *** sbc_blk() => computing q_air out of dew-point and P !'
            q_air_zt(:,:) = q_sat( sf(jp_humi )%fnow(:,:,1), sf(jp_slp  )%fnow(:,:,1) )
         CASE( np_humi_rlh )
            IF((kt==nit000).AND.lwp) WRITE(numout,*) ' *** sbc_blk() => computing q_air out of RH, t_air and slp !' !LBrm
            q_air_zt(:,:) = q_air_rh( 0.01_wp*sf(jp_humi )%fnow(:,:,1), &
               &                      sf(jp_tair )%fnow(:,:,1), sf(jp_slp  )%fnow(:,:,1) ) !#LB: 0.01 => RH is % percent in file
         END SELECT

         ! Potential temperature of air at z=rn_zqt (most reanalysis products provide absolute temp., not potential temp.)
         IF( ln_tair_pot ) THEN
            ! temperature read into file is already potential temperature, do nothing...
            theta_air_zt(:,:) = sf(jp_tair )%fnow(:,:,1)
         ELSE
            ! temperature read into file is ABSOLUTE temperature (that's the case for ECMWF products for example...)
            IF((kt==nit000).AND.lwp) WRITE(numout,*) ' *** sbc_blk() => air temperature converted from ABSOLUTE to POTENTIAL!'
            zpre(:,:)         = pres_temp( q_air_zt(:,:), sf(jp_slp)%fnow(:,:,1), rn_zu, pta=sf(jp_tair)%fnow(:,:,1) )
            theta_air_zt(:,:) = theta_exner( sf(jp_tair)%fnow(:,:,1), zpre(:,:) )
         ENDIF
         !
         CALL blk_oce_1( kt, sf(jp_wndi )%fnow(:,:,1), sf(jp_wndj )%fnow(:,:,1),   &   !   <<= in
            &                theta_air_zt(:,:), q_air_zt(:,:),                     &   !   <<= in
            &                sf(jp_slp  )%fnow(:,:,1), sst_m, ssu_m, ssv_m,        &   !   <<= in
            &                sf(jp_uoatm)%fnow(:,:,1), sf(jp_voatm)%fnow(:,:,1),   &   !   <<= in
            &                sf(jp_qsr  )%fnow(:,:,1), sf(jp_qlw  )%fnow(:,:,1),   &   !   <<= in (wl/cs)
            &                tsk_m, zssq, zcd_du, zsen, zlat, zevp )                   !   =>> out

         CALL blk_oce_2(     theta_air_zt(:,:),                                    &   !   <<= in
            &                sf(jp_qlw  )%fnow(:,:,1), sf(jp_prec )%fnow(:,:,1),   &   !   <<= in
            &                sf(jp_snow )%fnow(:,:,1), tsk_m,                      &   !   <<= in
            &                zsen, zlat, zevp )                                        !   <=> in out
      ENDIF
      !
#if defined key_cice
      IF( MOD( kt - 1, nn_fsbc ) == 0 )   THEN
         qlw_ice(:,:,1)   = sf(jp_qlw )%fnow(:,:,1)
         IF( ln_dm2dc ) THEN
            qsr_ice(:,:,1) = sbc_dcy( sf(jp_qsr)%fnow(:,:,1) )
         ELSE
            qsr_ice(:,:,1) =          sf(jp_qsr)%fnow(:,:,1)
         ENDIF
         tatm_ice(:,:) = sf(jp_tair)%fnow(:,:,1)    !#LB: should it be POTENTIAL temperature (theta_air_zt) instead ????
         qatm_ice(:,:) = q_air_zt(:,:)
         tprecip(:,:)  = sf(jp_prec)%fnow(:,:,1) * rn_pfac
         sprecip(:,:)  = sf(jp_snow)%fnow(:,:,1) * rn_pfac
         wndi_ice(:,:) = sf(jp_wndi)%fnow(:,:,1)
         wndj_ice(:,:) = sf(jp_wndj)%fnow(:,:,1)
      ENDIF
#endif

#if defined key_top
      IF( ln_trcdc2dm )  THEN      !  diurnal cycle in TOP
         IF( MOD( kt - 1, nn_fsbc ) == 0 ) THEN
            IF( ln_dm2dc )  THEN
                qsr_mean(:,:) = ( 1. - albo )  * sf(jp_qsr)%fnow(:,:,1)  * tmask(:,:,1)
            ELSE
                ncpl_qsr_freq = sf(jp_qsr)%freqh * 3600 !   qsr_mean will be computed in TOP
            ENDIF
         ENDIF
      ENDIF
#endif
      !
   END SUBROUTINE sbc_blk


   SUBROUTINE blk_oce_1( kt, pwndi, pwndj, ptair, pqair,         &  ! inp
      &                      pslp , pst  , pu   , pv,            &  ! inp
      &                      puatm, pvatm, pdqsr , pdqlw ,       &  ! inp
      &                      ptsk , pssq , pcd_du, psen, plat, pevp ) ! out
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE blk_oce_1  ***
      !!
      !! ** Purpose :   if ln_blk=T, computes surface momentum, heat and freshwater fluxes
      !!                if ln_abl=T, computes Cd x |U|, Ch x |U|, Ce x |U| for ABL integration
      !!
      !! ** Method  :   bulk formulae using atmospheric fields from :
      !!                if ln_blk=T, atmospheric fields read in sbc_read
      !!                if ln_abl=T, the ABL model at previous time-step
      !!
      !! ** Outputs : - pssq    : surface humidity used to compute latent heat flux (kg/kg)
      !!              - pcd_du  : Cd x |dU| at T-points  (m/s)
      !!              - psen    : sensible heat flux (W/m^2)
      !!              - plat    : latent heat flux   (W/m^2)
      !!              - pevp    : evaporation        (mm/s) #lolo
      !!---------------------------------------------------------------------
      INTEGER , INTENT(in   )                 ::   kt     ! time step index
      REAL(wp), INTENT(in   ), DIMENSION(:,:) ::   pwndi  ! atmospheric wind at T-point              [m/s]
      REAL(wp), INTENT(in   ), DIMENSION(:,:) ::   pwndj  ! atmospheric wind at T-point              [m/s]
      REAL(wp), INTENT(in   ), DIMENSION(:,:) ::   pqair  ! specific humidity at T-points            [kg/kg]
      REAL(wp), INTENT(in   ), DIMENSION(:,:) ::   ptair  ! potential temperature at T-points        [Kelvin]
      REAL(wp), INTENT(in   ), DIMENSION(:,:) ::   pslp   ! sea-level pressure                       [Pa]
      REAL(wp), INTENT(in   ), DIMENSION(:,:) ::   pst    ! surface temperature                      [Celsius]
      REAL(wp), INTENT(in   ), DIMENSION(:,:) ::   pu     ! surface current at U-point (i-component) [m/s]
      REAL(wp), INTENT(in   ), DIMENSION(:,:) ::   pv     ! surface current at V-point (j-component) [m/s]
      REAL(wp), INTENT(in   ), DIMENSION(:,:) ::   puatm  ! surface current seen by the atm at T-point (i-component) [m/s]
      REAL(wp), INTENT(in   ), DIMENSION(:,:) ::   pvatm  ! surface current seen by the atm at T-point (j-component) [m/s]
      REAL(wp), INTENT(in   ), DIMENSION(:,:) ::   pdqsr  ! downwelling solar (shortwave) radiation at surface [W/m^2]
      REAL(wp), INTENT(in   ), DIMENSION(:,:) ::   pdqlw  ! downwelling longwave radiation at surface [W/m^2]
      REAL(wp), INTENT(  out), DIMENSION(:,:) ::   ptsk   ! skin temp. (or SST if CS & WL not used)  [Celsius]
      REAL(wp), INTENT(  out), DIMENSION(:,:) ::   pssq   ! specific humidity at pst                 [kg/kg]
      REAL(wp), INTENT(  out), DIMENSION(:,:) ::   pcd_du
      REAL(wp), INTENT(  out), DIMENSION(:,:) ::   psen
      REAL(wp), INTENT(  out), DIMENSION(:,:) ::   plat
      REAL(wp), INTENT(  out), DIMENSION(:,:) ::   pevp
      !
      INTEGER  ::   ji, jj               ! dummy loop indices
      REAL(wp) ::   zztmp                ! local variable
      REAL(wp) ::   zstmax, zstau
#if defined key_cyclone
      REAL(wp), DIMENSION(jpi,jpj) ::   zwnd_i, zwnd_j    ! wind speed components at T-point
#endif
      REAL(wp), DIMENSION(jpi,jpj) ::   ztau_i, ztau_j    ! wind stress components at T-point
      REAL(wp), DIMENSION(jpi,jpj) ::   zU_zu             ! bulk wind speed at height zu  [m/s]
      REAL(wp), DIMENSION(jpi,jpj) ::   zcd_oce           ! momentum transfert coefficient over ocean
      REAL(wp), DIMENSION(jpi,jpj) ::   zch_oce           ! sensible heat transfert coefficient over ocean
      REAL(wp), DIMENSION(jpi,jpj) ::   zce_oce           ! latent   heat transfert coefficient over ocean
      REAL(wp), DIMENSION(jpi,jpj) ::   zsspt             ! potential sea-surface temperature [K]
      REAL(wp), DIMENSION(jpi,jpj) ::   zpre, ztabs       ! air pressure [Pa] & absolute temperature [K]
      REAL(wp), DIMENSION(jpi,jpj) ::   zztmp1, zztmp2
#if defined key_drakkar
      REAL(wp), DIMENSION(jpi,jpj) ::   zwu, zwv
#endif
      !!---------------------------------------------------------------------
      !
      ! local scalars ( place there for vector optimisation purposes)
      !                           ! Temporary conversion from Celcius to Kelvin (and set minimum value far above 0 K)
      ptsk(:,:) = pst(:,:) + rt0  ! by default: skin temperature = "bulk SST" (will remain this way if NCAR algorithm used!)

      ! sea surface potential temperature [K]
      zsspt(:,:) = theta_exner( ptsk(:,:), pslp(:,:) )

      ! --- cloud cover --- !
      cloud_fra(:,:) = sf(jp_cc)%fnow(:,:,1)

      ! ----------------------------------------------------------------------------- !
      !      0   Wind components and module at T-point relative to the moving ocean   !
      ! ----------------------------------------------------------------------------- !

      ! ... components ( U10m - U_oce ) at T-point (unmasked)
#if defined key_cyclone
      zwnd_i(:,:) = 0._wp
      zwnd_j(:,:) = 0._wp
      CALL wnd_cyc( kt, zwnd_i, zwnd_j )    ! add analytical tropical cyclone (Vincent et al. JGR 2012)
      DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
         zwnd_i(ji,jj) = pwndi(ji,jj) + zwnd_i(ji,jj)
         zwnd_j(ji,jj) = pwndj(ji,jj) + zwnd_j(ji,jj)
         ! ... scalar wind at T-point (not masked)
         wndm(ji,jj) = SQRT( zwnd_i(ji,jj) * zwnd_i(ji,jj) + zwnd_j(ji,jj) * zwnd_j(ji,jj) )
      END_2D
#else
#if defined key_drakkar
    IF (ln_clim_forcing) THEN
       wndm(:,:) = sf(jp_wmod)%fnow(:,:,1)
    ELSE
#endif
      ! ... scalar wind module at T-point (not masked)
      DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
         wndm(ji,jj) = SQRT(  pwndi(ji,jj) * pwndi(ji,jj) + pwndj(ji,jj) * pwndj(ji,jj)  )
      END_2D
#if defined key_drakkar
    ENDIF
#endif
#endif
      ! ----------------------------------------------------------------------------- !
      !      I   Solar FLUX                                                           !
      ! ----------------------------------------------------------------------------- !

      ! ocean albedo assumed to be constant + modify now Qsr to include the diurnal cycle                    ! Short Wave
      zztmp = 1. - albo
      IF( ln_dm2dc ) THEN
         qsr(:,:) = zztmp * sbc_dcy( pdqsr(:,:) ) * tmask(:,:,1)
      ELSE
         qsr(:,:) = zztmp *          pdqsr(:,:)   * tmask(:,:,1)
      ENDIF


      ! ----------------------------------------------------------------------------- !
      !     II   Turbulent FLUXES                                                     !
      ! ----------------------------------------------------------------------------- !

      ! specific humidity at SST
      pssq(:,:) = rdct_qsat_salt * q_sat( ptsk(:,:), pslp(:,:) )

      IF( ln_skin_cs .OR. ln_skin_wl ) THEN
         !! Backup "bulk SST" and associated spec. hum.
         zztmp1(:,:) = zsspt(:,:)
         zztmp2(:,:) = pssq(:,:)
      ENDIF

      !! Time to call the user-selected bulk parameterization for
      !!  ==  transfer coefficients  ==!   Cd, Ch, Ce at T-point, and more...
      SELECT CASE( nblk )

      CASE( np_NCAR      )
         CALL turb_ncar    (     rn_zqt, rn_zu, zsspt, ptair, pssq, pqair, wndm, &
            &                zcd_oce, zch_oce, zce_oce, theta_zu, q_zu, zU_zu , &
            &                nb_iter=nn_iter_algo )
         !
      CASE( np_COARE_3p0 )
         CALL turb_coare3p0( kt, rn_zqt, rn_zu, zsspt, ptair, pssq, pqair, wndm, &
            &                ln_skin_cs, ln_skin_wl,                            &
            &                zcd_oce, zch_oce, zce_oce, theta_zu, q_zu, zU_zu,  &
            &                nb_iter=nn_iter_algo,                              &
            &                Qsw=qsr(:,:), rad_lw=pdqlw(:,:), slp=pslp(:,:) )
         !
      CASE( np_COARE_3p6 )
         CALL turb_coare3p6( kt, rn_zqt, rn_zu, zsspt, ptair, pssq, pqair, wndm, &
            &                ln_skin_cs, ln_skin_wl,                            &
            &                zcd_oce, zch_oce, zce_oce, theta_zu, q_zu, zU_zu,  &
            &                nb_iter=nn_iter_algo,                              &
            &                Qsw=qsr(:,:), rad_lw=pdqlw(:,:), slp=pslp(:,:) )
         !
      CASE( np_ECMWF     )
         CALL turb_ecmwf   ( kt, rn_zqt, rn_zu, zsspt, ptair, pssq, pqair, wndm, &
            &                ln_skin_cs, ln_skin_wl,                            &
            &                zcd_oce, zch_oce, zce_oce, theta_zu, q_zu, zU_zu,  &
            &                nb_iter=nn_iter_algo,                              &
            &                Qsw=qsr(:,:), rad_lw=pdqlw(:,:), slp=pslp(:,:) )
         !
      CASE( np_ANDREAS   )
         CALL turb_andreas (     rn_zqt, rn_zu, zsspt, ptair, pssq, pqair, wndm, &
            &                zcd_oce, zch_oce, zce_oce, theta_zu, q_zu, zU_zu , &
            &                nb_iter=nn_iter_algo   )
         !
      CASE DEFAULT
         CALL ctl_stop( 'STOP', 'sbc_oce: non-existing bulk parameterizaton selected' )
         !
      END SELECT

      IF( iom_use('Cd_oce') )   CALL iom_put("Cd_oce",   zcd_oce * tmask(:,:,1))
      IF( iom_use('Ce_oce') )   CALL iom_put("Ce_oce",   zce_oce * tmask(:,:,1))
      IF( iom_use('Ch_oce') )   CALL iom_put("Ch_oce",   zch_oce * tmask(:,:,1))
      !! LB: mainly here for debugging purpose:
      IF( iom_use('theta_zt') ) CALL iom_put("theta_zt", (ptair-rt0) * tmask(:,:,1)) ! potential temperature at z=zt
      IF( iom_use('q_zt') )     CALL iom_put("q_zt",     pqair       * tmask(:,:,1)) ! specific humidity       "
      IF( iom_use('theta_zu') ) CALL iom_put("theta_zu", (theta_zu -rt0) * tmask(:,:,1)) ! potential temperature at z=zu
      IF( iom_use('q_zu') )     CALL iom_put("q_zu",     q_zu        * tmask(:,:,1)) ! specific humidity       "
      IF( iom_use('ssq') )      CALL iom_put("ssq",      pssq        * tmask(:,:,1)) ! saturation specific humidity at z=0
      IF( iom_use('wspd_blk') ) CALL iom_put("wspd_blk", zU_zu       * tmask(:,:,1)) ! bulk wind speed at z=zu

      IF( ln_skin_cs .OR. ln_skin_wl ) THEN
         !! In the presence of sea-ice we forget about the cool-skin/warm-layer update of zsspt, pssq & ptsk:
         WHERE ( fr_i(:,:) > 0.001_wp )
            ! sea-ice present, we forget about the update, using what we backed up before call to turb_*()
            zsspt(:,:) = zztmp1(:,:)
            pssq(:,:)  = zztmp2(:,:)
         END WHERE
         ! apply potential temperature increment to abolute SST
         ptsk(:,:) = ptsk(:,:) + ( zsspt(:,:) - zztmp1(:,:) )
      END IF

      !  Turbulent fluxes over ocean  => BULK_FORMULA @ sbc_phy.F90
      ! -------------------------------------------------------------

      IF( ln_abl ) THEN         !==  ABL formulation  ==!   multiplication by rho_air and turbulent fluxes computation done in ablstp

         DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
            zztmp = zU_zu(ji,jj)
            wndm(ji,jj)   = zztmp                   ! Store zU_zu in wndm to compute ustar2 in ablmod
            pcd_du(ji,jj) = zztmp * zcd_oce(ji,jj)
            psen(ji,jj)   = zztmp * zch_oce(ji,jj)
            pevp(ji,jj)   = zztmp * zce_oce(ji,jj)
            zpre(ji,jj)   = pres_temp( pqair(ji,jj), pslp(ji,jj), rn_zu, ptpot=ptair(ji,jj), pta=ztabs(ji,jj) )
            rhoa(ji,jj)   = rho_air( ztabs(ji,jj), pqair(ji,jj), zpre(ji,jj) )
         END_2D

      ELSE                      !==  BLK formulation  ==!   turbulent fluxes computation

         DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
            zpre(ji,jj) = pres_temp( q_zu(ji,jj), pslp(ji,jj), rn_zu, ptpot=theta_zu(ji,jj), pta=ztabs(ji,jj) )
            rhoa(ji,jj) = rho_air( ztabs(ji,jj), q_zu(ji,jj), zpre(ji,jj) )
         END_2D

         CALL BULK_FORMULA( rn_zu, zsspt(:,:), pssq(:,:), theta_zu(:,:), q_zu(:,:), &
            &               zcd_oce(:,:), zch_oce(:,:), zce_oce(:,:),          &
            &               wndm(:,:), zU_zu(:,:), pslp(:,:), rhoa(:,:),       &
            &               taum(:,:), psen(:,:), plat(:,:),                   &
            &               pEvap=pevp(:,:), pfact_evap=rn_efac )

         psen(:,:) = psen(:,:) * tmask(:,:,1)
         plat(:,:) = plat(:,:) * tmask(:,:,1)
         taum(:,:) = taum(:,:) * tmask(:,:,1)
         pevp(:,:) = pevp(:,:) * tmask(:,:,1)

#if defined key_drakkar
        IF ( ln_clim_forcing ) THEN
            zwu(:,:) =  sf(jp_uw)%fnow(:,:,1)
            zwv(:,:) =  sf(jp_vw)%fnow(:,:,1)

            DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
               IF( wndm(ji,jj) > 0._wp ) THEN
                  zztmp = taum(ji,jj) / wndm(ji,jj) / wndm(ji,jj)
            ! note that key_cyclone is not supported with climatological forcing
                  ztau_i(ji,jj) = zztmp * zwu(ji,jj)
                  ztau_j(ji,jj) = zztmp * zwv(ji,jj)
               ELSE
                 ztau_i(ji,jj) = 0._wp
                 ztau_j(ji,jj) = 0._wp
               ENDIF
            END_2D
        ELSE
#endif

         DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
            IF( wndm(ji,jj) > 0._wp ) THEN
              zztmp = taum(ji,jj) / wndm(ji,jj)
#if defined key_cyclone
               ztau_i(ji,jj) = zztmp * zwnd_i(ji,jj)
               ztau_j(ji,jj) = zztmp * zwnd_j(ji,jj)
#else
               ztau_i(ji,jj) = zztmp * pwndi(ji,jj)
               ztau_j(ji,jj) = zztmp * pwndj(ji,jj)
#endif
            ELSE
               ztau_i(ji,jj) = 0._wp
               ztau_j(ji,jj) = 0._wp
            ENDIF
         END_2D
#if defined key_drakkar
        ENDIF
#endif

         IF( ln_crt_fbk ) THEN   ! aply eq. 10 and 11 of Renault et al. 2020 (doi: 10.1029/2019MS001715)
            zstmax = MIN( rn_stau_a * 3._wp + rn_stau_b, 0._wp )   ! set the max value of Stau corresponding to a wind of 3 m/s (<0)
            DO_2D( 0, 1, 0, 1 )   ! end at jpj and jpi, as ztau_j(ji,jj+1) ztau_i(ji+1,jj) used in the next loop
               zstau = MIN( rn_stau_a * wndm(ji,jj) + rn_stau_b, zstmax )   ! stau (<0) must be smaller than zstmax
               ztau_i(ji,jj) = ztau_i(ji,jj) + zstau * ( 0.5_wp * ( pu(ji-1,jj  ) + pu(ji,jj) ) - puatm(ji,jj) )
               ztau_j(ji,jj) = ztau_j(ji,jj) + zstau * ( 0.5_wp * ( pv(ji  ,jj-1) + pv(ji,jj) ) - pvatm(ji,jj) )
               taum(ji,jj) = SQRT( ztau_i(ji,jj) * ztau_i(ji,jj) + ztau_j(ji,jj) * ztau_j(ji,jj) )
            END_2D
         ENDIF

         ! ... utau, vtau at U- and V_points, resp.
         !     Note the use of 0.5*(2-umask) in order to unmask the stress along coastlines
         !     Note that coastal wind stress is not used in the code... so this extra care has no effect
         DO_2D( 0, 0, 0, 0 )              ! start loop at 2, in case ln_crt_fbk = T
            utau(ji,jj) = 0.5 * ( 2. - umask(ji,jj,1) ) * ( ztau_i(ji,jj) + ztau_i(ji+1,jj  ) ) &
               &              * MAX(tmask(ji,jj,1),tmask(ji+1,jj,1))
            vtau(ji,jj) = 0.5 * ( 2. - vmask(ji,jj,1) ) * ( ztau_j(ji,jj) + ztau_j(ji  ,jj+1) ) &
               &              * MAX(tmask(ji,jj,1),tmask(ji,jj+1,1))
         END_2D

         IF( ln_crt_fbk ) THEN
            CALL lbc_lnk( 'sbcblk', utau, 'U', -1._wp, vtau, 'V', -1._wp, taum, 'T', 1._wp )
         ELSE
            CALL lbc_lnk( 'sbcblk', utau, 'U', -1._wp, vtau, 'V', -1._wp )
         ENDIF

         ! Saving open-ocean wind-stress (module and components) on T-points:
         CALL iom_put( "taum_oce",   taum(:,:)*tmask(:,:,1) )   ! output wind stress module
         !#LB: These 2 lines below mostly here for 'STATION_ASF' test-case, otherwize "utau" (U-grid) and vtau" (V-grid) does the job in: [DYN/dynatf.F90])
         CALL iom_put( "utau_oce", ztau_i(:,:)*tmask(:,:,1) )  ! utau at T-points!
         CALL iom_put( "vtau_oce", ztau_j(:,:)*tmask(:,:,1) )  ! vtau at T-points!

         IF(sn_cfctl%l_prtctl) THEN
            CALL prt_ctl( tab2d_1=pssq   , clinfo1=' blk_oce_1: pssq   : ', mask1=tmask )
            CALL prt_ctl( tab2d_1=wndm   , clinfo1=' blk_oce_1: wndm   : ', mask1=tmask )
            CALL prt_ctl( tab2d_1=utau   , clinfo1=' blk_oce_1: utau   : ', mask1=umask,   &
               &          tab2d_2=vtau   , clinfo2='            vtau   : ', mask2=vmask )
            CALL prt_ctl( tab2d_1=zcd_oce, clinfo1=' blk_oce_1: Cd     : ', mask1=tmask )
         ENDIF
         !
      ENDIF ! ln_blk / ln_abl

      ptsk(:,:) = ( ptsk(:,:) - rt0 ) * tmask(:,:,1)  ! Back to Celsius

      IF( ln_skin_cs .OR. ln_skin_wl ) THEN
         CALL iom_put( "t_skin" ,  ptsk        )  ! T_skin in Celsius
         CALL iom_put( "dt_skin" , ptsk - pst  )  ! T_skin - SST temperature difference
      ENDIF
      !
   END SUBROUTINE blk_oce_1


   SUBROUTINE blk_oce_2( ptair, pdqlw, pprec, psnow, &   ! <<= in
      &                   ptsk, psen, plat, pevp     )   ! <<= in
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE blk_oce_2  ***
      !!
      !! ** Purpose :   finalize the momentum, heat and freshwater fluxes computation
      !!                at the ocean surface at each time step knowing Cd, Ch, Ce and
      !!                atmospheric variables (from ABL or external data)
      !!
#if defined key_drakkar
      !! ** Outputs : - emp     :  evaporation minus precipitation       (kg/m2/s)
      !!              - qns     : Non Solar heat flux over the ocean    (W/m2)
      !!                          (include heat content of the precip, evap
      !!                           latent heat flux from melting snow )
      !!              - qns_oce : Just the sensible, latent and LW (to be used in SI3)
      !!              - qsr_oce : = qsr   (to be used in SI3)
#else
      !! ** Outputs : - utau    : i-component of the stress at U-point  (N/m2)
      !!              - vtau    : j-component of the stress at V-point  (N/m2)
      !!              - taum    : Wind stress module at T-point         (N/m2)
      !!              - wndm    : Wind speed module at T-point          (m/s)
      !!              - qsr     : Solar heat flux over the ocean        (W/m2)
      !!              - qns     : Non Solar heat flux over the ocean    (W/m2)
      !!              - emp     : evaporation minus precipitation       (kg/m2/s)
      !!---------------------------------------------------------------------
#endif
      REAL(wp), INTENT(in), DIMENSION(:,:) ::   ptair   ! potential temperature of air #LB: confirm!
      REAL(wp), INTENT(in), DIMENSION(:,:) ::   pdqlw   ! downwelling longwave radiation at surface [W/m^2]
      REAL(wp), INTENT(in), DIMENSION(:,:) ::   pprec
      REAL(wp), INTENT(in), DIMENSION(:,:) ::   psnow
      REAL(wp), INTENT(in), DIMENSION(:,:) ::   ptsk   ! SKIN surface temperature   [Celsius]
      REAL(wp), INTENT(in), DIMENSION(:,:) ::   psen
      REAL(wp), INTENT(in), DIMENSION(:,:) ::   plat
      REAL(wp), INTENT(in), DIMENSION(:,:) ::   pevp
      !
      INTEGER  ::   ji, jj               ! dummy loop indices
      REAL(wp) ::   zztmp,zz1,zz2,zz3    ! local variable
      REAL(wp), DIMENSION(jpi,jpj) ::   zqlw              ! net long wave radiative heat flux
      REAL(wp), DIMENSION(jpi,jpj) ::   zcptrain, zcptsnw, zcptn ! Heat content per unit mass (J/kg)
      !!---------------------------------------------------------------------
      !
      ! Heat content per unit mass (J/kg)
      zcptrain(:,:) = (      ptair        - rt0 ) * rcp  * tmask(:,:,1)
      zcptsnw (:,:) = ( MIN( ptair, rt0 ) - rt0 ) * rcpi * tmask(:,:,1)
      zcptn   (:,:) =        ptsk                 * rcp  * tmask(:,:,1)
      !
      ! ----------------------------------------------------------------------------- !
      !     III    Net longwave radiative FLUX                                        !
      ! ----------------------------------------------------------------------------- !
      !! #LB: now moved after Turbulent fluxes because must use the skin temperature rather than bulk SST
      !! (ptsk is skin temperature if ln_skin_cs==.TRUE. .OR. ln_skin_wl==.TRUE.)
      zqlw(:,:) = qlw_net( pdqlw(:,:), ptsk(:,:)+rt0 )

      ! ----------------------------------------------------------------------------- !
      !     IV    Total FLUXES                                                       !
      ! ----------------------------------------------------------------------------- !
      !
      emp (:,:) = ( pevp(:,:) - pprec(:,:) * rn_pfac ) * tmask(:,:,1)      ! mass flux (evap. - precip.)
      !
      qns(:,:) = zqlw(:,:) + psen(:,:) + plat(:,:)                     &   ! Downward Non Solar
         &     - psnow(:,:) * rn_pfac * rLfus                          &   ! remove latent melting heat for solid precip
         &     - pevp(:,:) * zcptn(:,:)                                &   ! remove evap heat content at SST
         &     + ( pprec(:,:) - psnow(:,:) ) * rn_pfac * zcptrain(:,:) &   ! add liquid precip heat content at Tair
         &     + psnow(:,:) * rn_pfac * zcptsnw(:,:)                       ! add solid  precip heat content at min(Tair,Tsnow)
      qns(:,:) = qns(:,:) * tmask(:,:,1)
      !
#if defined key_si3
      IF ( nn_ice == 2 ) THEN
         qns_oce(:,:) = zqlw(:,:) + psen(:,:) + plat(:,:)                  ! non solar without emp (only needed by SI3)
         qsr_oce(:,:) = qsr(:,:)
      ENDIF
#endif
      !
      CALL iom_put( "rho_air"  , rhoa*tmask(:,:,1) )       ! output air density [kg/m^3]
      CALL iom_put( "evap_oce" , pevp )                    ! evaporation
      CALL iom_put( "qlw_oce"  , zqlw )                    ! output downward longwave heat over the ocean
      CALL iom_put( "qsb_oce"  , psen )                    ! output downward sensible heat over the ocean
      CALL iom_put( "qla_oce"  , plat )                    ! output downward latent   heat over the ocean
      tprecip(:,:) = pprec(:,:) * rn_pfac * tmask(:,:,1)   ! output total precipitation [kg/m2/s]
      sprecip(:,:) = psnow(:,:) * rn_pfac * tmask(:,:,1)   ! output solid precipitation [kg/m2/s]
      CALL iom_put( 'snowpre', sprecip )                   ! Snow
      CALL iom_put( 'precip' , tprecip )                   ! Total precipitation
      !
      IF ( nn_ice == 0 ) THEN
         CALL iom_put( "qemp_oce" , qns-zqlw-psen-plat )   ! output downward heat content of E-P over the ocean
         CALL iom_put( "qns_oce"  ,   qns  )               ! output downward non solar heat over the ocean
         CALL iom_put( "qsr_oce"  ,   qsr  )               ! output downward solar heat over the ocean
         CALL iom_put( "qt_oce"   ,   qns+qsr )            ! output total downward heat over the ocean
      ENDIF
      !
      IF(sn_cfctl%l_prtctl) THEN
         CALL prt_ctl(tab2d_1=zqlw , clinfo1=' blk_oce_2: zqlw  : ', mask1=tmask )
         CALL prt_ctl(tab2d_1=psen , clinfo1=' blk_oce_2: psen  : ', mask1=tmask )
         CALL prt_ctl(tab2d_1=plat , clinfo1=' blk_oce_2: plat  : ', mask1=tmask )
         CALL prt_ctl(tab2d_1=qns  , clinfo1=' blk_oce_2: qns   : ', mask1=tmask )
         CALL prt_ctl(tab2d_1=emp  , clinfo1=' blk_oce_2: emp   : ', mask1=tmask )
      ENDIF
      !
   END SUBROUTINE blk_oce_2


#if defined key_si3
   !!----------------------------------------------------------------------
   !!   'key_si3'                                       SI3 sea-ice model
   !!----------------------------------------------------------------------
   !!   blk_ice_1   : provide the air-ice stress
   !!   blk_ice_2   : provide the heat and mass fluxes at air-ice interface
   !!   blk_ice_qcn : provide ice surface temperature and snow/ice conduction flux (emulating conduction flux)
   !!----------------------------------------------------------------------

   SUBROUTINE blk_ice_1( pwndi, pwndj, ptair, pqair, pslp , puice, pvice, ptsui,  &   ! inputs
      &                  putaui, pvtaui, pseni, pevpi, pssqi, pcd_dui             )   ! optional outputs
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE blk_ice_1  ***
      !!
      !! ** Purpose :   provide the surface boundary condition over sea-ice
      !!
      !! ** Method  :   compute momentum using bulk formulation
      !!                formulea, ice variables and read atmospheric fields.
      !!                NB: ice drag coefficient is assumed to be a constant
      !!---------------------------------------------------------------------
      REAL(wp) , INTENT(in   ), DIMENSION(:,:  ) ::   pslp    ! sea-level pressure [Pa]
      REAL(wp) , INTENT(in   ), DIMENSION(:,:  ) ::   pwndi   ! atmospheric wind at T-point [m/s]
      REAL(wp) , INTENT(in   ), DIMENSION(:,:  ) ::   pwndj   ! atmospheric wind at T-point [m/s]
      REAL(wp) , INTENT(in   ), DIMENSION(:,:  ) ::   ptair   ! atmospheric potential temperature at T-point [K]
      REAL(wp) , INTENT(in   ), DIMENSION(:,:  ) ::   pqair   ! atmospheric specific humidity at T-point [kg/kg]
      REAL(wp) , INTENT(in   ), DIMENSION(:,:  ) ::   puice   ! sea-ice velocity on I or C grid [m/s]
      REAL(wp) , INTENT(in   ), DIMENSION(:,:  ) ::   pvice   ! "
      REAL(wp) , INTENT(in   ), DIMENSION(:,:  ) ::   ptsui   ! sea-ice surface temperature [K]
      REAL(wp) , INTENT(  out), DIMENSION(:,:  ), OPTIONAL ::   putaui  ! if ln_blk
      REAL(wp) , INTENT(  out), DIMENSION(:,:  ), OPTIONAL ::   pvtaui  ! if ln_blk
      REAL(wp) , INTENT(  out), DIMENSION(:,:  ), OPTIONAL ::   pseni   ! if ln_abl
      REAL(wp) , INTENT(  out), DIMENSION(:,:  ), OPTIONAL ::   pevpi   ! if ln_abl
      REAL(wp) , INTENT(  out), DIMENSION(:,:  ), OPTIONAL ::   pssqi   ! if ln_abl
      REAL(wp) , INTENT(  out), DIMENSION(:,:  ), OPTIONAL ::   pcd_dui ! if ln_abl
      !
      INTEGER  ::   ji, jj    ! dummy loop indices
      REAL(wp) ::   zootm_su                      ! sea-ice surface mean temperature
      REAL(wp) ::   zztmp1, zztmp2                ! temporary scalars
      REAL(wp), DIMENSION(jpi,jpj) :: ztmp, zsipt ! temporary array
#if defined key_drakkar
      REAL(wp), DIMENSION(jpi,jpj) :: zwu, zwv    ! temporary array
#endif
      !!---------------------------------------------------------------------
      !
      ! ------------------------------------------------------------ !
      !    Wind module relative to the moving ice ( U10m - U_ice )   !
      ! ------------------------------------------------------------ !
      ! C-grid ice dynamics :   U & V-points (same as ocean)
#if defined key_drakkar
      IF ( ln_clim_forcing) THEN
        wndm_ice(:,:) = sf(jp_wmod)%fnow(:,:,1)
      ELSE
#endif
      DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
         wndm_ice(ji,jj) = SQRT( pwndi(ji,jj) * pwndi(ji,jj) + pwndj(ji,jj) * pwndj(ji,jj) )
      END_2D
#if defined key_drakkar
      ENDIF
#endif
      !
      ! potential sea-ice surface temperature [K]
      zsipt(:,:) = theta_exner( ptsui(:,:), pslp(:,:) )

      ! sea-ice <-> atmosphere bulk transfer coefficients
      SELECT CASE( nblk_ice )

      CASE( np_ice_cst      )
         ! Constant bulk transfer coefficients over sea-ice:
         Cd_ice(:,:) = rn_Cd_i
         Ch_ice(:,:) = rn_Ch_i
         Ce_ice(:,:) = rn_Ce_i
         ! no height adjustment, keeping zt values:
         theta_zu_i(:,:) = ptair(:,:)
         q_zu_i(:,:)     = pqair(:,:)

      CASE( np_ice_an05 )  ! calculate new drag from Lupkes(2015) equations
         ztmp(:,:) = q_sat( ptsui(:,:), pslp(:,:), l_ice=.TRUE. ) ! temporary array for SSQ
         CALL turb_ice_an05( rn_zqt, rn_zu, zsipt, ptair, ztmp, pqair, wndm_ice,       &
            &                      Cd_ice, Ch_ice, Ce_ice, theta_zu_i, q_zu_i )
         !!
      CASE( np_ice_lu12 )
         ztmp(:,:) = q_sat( ptsui(:,:), pslp(:,:), l_ice=.TRUE. ) ! temporary array for SSQ
         CALL turb_ice_lu12( rn_zqt, rn_zu, zsipt, ptair, ztmp, pqair, wndm_ice, fr_i, &
            &                      Cd_ice, Ch_ice, Ce_ice, theta_zu_i, q_zu_i )
         !!
      CASE( np_ice_lg15 )  ! calculate new drag from Lupkes(2015) equations
         ztmp(:,:) = q_sat( ptsui(:,:), pslp(:,:), l_ice=.TRUE. ) ! temporary array for SSQ
         CALL turb_ice_lg15( rn_zqt, rn_zu, zsipt, ptair, ztmp, pqair, wndm_ice, fr_i, &
            &                      Cd_ice, Ch_ice, Ce_ice, theta_zu_i, q_zu_i )
         !!
      END SELECT

      IF( iom_use('Cd_ice').OR.iom_use('Ce_ice').OR.iom_use('Ch_ice').OR.iom_use('taum_ice').OR.iom_use('utau_ice').OR.iom_use('vtau_ice') ) &
         & ztmp(:,:) = ( 1._wp - MAX(0._wp, SIGN( 1._wp, 1.E-6_wp - fr_i )) )*tmask(:,:,1) ! mask for presence of ice !

      IF( iom_use('Cd_ice') ) CALL iom_put("Cd_ice", Cd_ice*ztmp)
      IF( iom_use('Ce_ice') ) CALL iom_put("Ce_ice", Ce_ice*ztmp)
      IF( iom_use('Ch_ice') ) CALL iom_put("Ch_ice", Ch_ice*ztmp)


      IF( ln_blk ) THEN
         ! ---------------------------------------------------- !
         !    Wind stress relative to nonmoving ice ( U10m )    !
         ! ---------------------------------------------------- !
         ! supress moving ice in wind stress computation as we don't know how to do it properly...
#if defined key_drakkar
      IF ( ln_clim_forcing ) THEN
         zwu(:,:) = sf(jp_uw)%fnow(:,:,1)
         zwv(:,:) = sf(jp_vw)%fnow(:,:,1)
         DO_2D( 0, 1, 0, 1 )    ! at T point
            zztmp1        = rhoa(ji,jj) * Cd_ice(ji,jj) 
            putaui(ji,jj) = zztmp1 * zwu (ji,jj)
            pvtaui(ji,jj) = zztmp1 * zwv (ji,jj)
         END_2D
      ELSE
#endif
         DO_2D( 0, 1, 0, 1 )    ! at T point
            zztmp1        = rhoa(ji,jj) * Cd_ice(ji,jj) * wndm_ice(ji,jj)
            putaui(ji,jj) = zztmp1 * pwndi(ji,jj)
            pvtaui(ji,jj) = zztmp1 * pwndj(ji,jj)
         END_2D
#if defined key_drakkar
      ENDIF
#endif

         !#LB: saving the module, and x-y components, of the ai wind-stress at T-points: NOT weighted by the ice concentration !!!
         IF(iom_use('taum_ice')) CALL iom_put('taum_ice', SQRT( putaui*putaui + pvtaui*pvtaui )*ztmp )
         !#LB: These 2 lines below mostly here for 'STATION_ASF' test-case, otherwize "utau_oi" (U-grid) and vtau_oi" (V-grid) does the job in: [ICE/icedyn_rhg_evp.F90])
         IF(iom_use('utau_ice')) CALL iom_put("utau_ice", putaui*ztmp)  ! utau at T-points!
         IF(iom_use('vtau_ice')) CALL iom_put("vtau_ice", pvtaui*ztmp)  ! vtau at T-points!

         !
         DO_2D( 0, 0, 0, 0 )    ! U & V-points (same as ocean).
            !#LB: QUESTION?? so SI3 expects wind stress vector to be provided at U & V points? Not at T-points ?
            ! take care of the land-sea mask to avoid "pollution" of coastal stress. p[uv]taui used in frazil and  rheology
            zztmp1 = 0.5_wp * ( 2. - umask(ji,jj,1) ) * MAX( tmask(ji,jj,1),tmask(ji+1,jj  ,1) )
            zztmp2 = 0.5_wp * ( 2. - vmask(ji,jj,1) ) * MAX( tmask(ji,jj,1),tmask(ji  ,jj+1,1) )
            putaui(ji,jj) = zztmp1 * ( putaui(ji,jj) + putaui(ji+1,jj  ) )
            pvtaui(ji,jj) = zztmp2 * ( pvtaui(ji,jj) + pvtaui(ji  ,jj+1) )
         END_2D
         CALL lbc_lnk( 'sbcblk', putaui, 'U', -1._wp, pvtaui, 'V', -1._wp )
         !
         IF(sn_cfctl%l_prtctl)  CALL prt_ctl( tab2d_1=putaui  , clinfo1=' blk_ice: putaui : ', mask1=umask   &
            &                               , tab2d_2=pvtaui  , clinfo2='          pvtaui : ', mask2=vmask )
      ELSE ! ln_abl

         DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
            pcd_dui(ji,jj) = wndm_ice(ji,jj) * Cd_ice(ji,jj)
            pseni  (ji,jj) = wndm_ice(ji,jj) * Ch_ice(ji,jj)
            pevpi  (ji,jj) = wndm_ice(ji,jj) * Ce_ice(ji,jj)
         END_2D
         pssqi(:,:) = q_sat( ptsui(:,:), pslp(:,:), l_ice=.TRUE. ) ; ! more accurate way to obtain ssq !

      ENDIF ! ln_blk  / ln_abl
      !
      IF(sn_cfctl%l_prtctl)  CALL prt_ctl(tab2d_1=wndm_ice  , clinfo1=' blk_ice: wndm_ice : ', mask1=tmask )
      !
   END SUBROUTINE blk_ice_1


   SUBROUTINE blk_ice_2( ptsu, phs, phi, palb, ptair, pqair, pslp, pdqlw, pprec, psnow  )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE blk_ice_2  ***
      !!
      !! ** Purpose :   provide the heat and mass fluxes at air-ice interface
      !!
      !! ** Method  :   compute heat and freshwater exchanged
      !!                between atmosphere and sea-ice using bulk formulation
      !!                formulea, ice variables and read atmmospheric fields.
      !!
      !! caution : the net upward water flux has with mm/day unit
      !!---------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:,:), INTENT(in)  ::   ptsu   ! sea ice surface temperature [K]
      REAL(wp), DIMENSION(:,:,:), INTENT(in)  ::   phs    ! snow thickness
      REAL(wp), DIMENSION(:,:,:), INTENT(in)  ::   phi    ! ice thickness
      REAL(wp), DIMENSION(:,:,:), INTENT(in)  ::   palb   ! ice albedo (all skies)
      REAL(wp), DIMENSION(:,:  ), INTENT(in)  ::   ptair  ! potential temperature of air #LB: okay ???
      REAL(wp), DIMENSION(:,:  ), INTENT(in)  ::   pqair  ! specific humidity of air
      REAL(wp), DIMENSION(:,:  ), INTENT(in)  ::   pslp
      REAL(wp), DIMENSION(:,:  ), INTENT(in)  ::   pdqlw
      REAL(wp), DIMENSION(:,:  ), INTENT(in)  ::   pprec
      REAL(wp), DIMENSION(:,:  ), INTENT(in)  ::   psnow
      !!
      INTEGER  ::   ji, jj, jl               ! dummy loop indices
      REAL(wp) ::   zst, zst3, zsq, zsipt    ! local variable
      REAL(wp) ::   zcoef_dqlw, zcoef_dqla   !   -      -
      REAL(wp) ::   zztmp, zzblk, zztmp1, z1_rLsub   !   -      -
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::   zmsk   ! temporary mask for prt_ctl
      REAL(wp), DIMENSION(jpi,jpj,jpl) ::   z_qlw         ! long wave heat flux over ice
      REAL(wp), DIMENSION(jpi,jpj,jpl) ::   z_qsb         ! sensible  heat flux over ice
      REAL(wp), DIMENSION(jpi,jpj,jpl) ::   z_dqlw        ! long wave heat sensitivity over ice
      REAL(wp), DIMENSION(jpi,jpj,jpl) ::   z_dqsb        ! sensible  heat sensitivity over ice
      REAL(wp), DIMENSION(jpi,jpj)     ::   zevap, zsnw   ! evaporation and snw distribution after wind blowing (SI3)
      REAL(wp), DIMENSION(jpi,jpj)     ::   ztmp, ztmp2
      REAL(wp), DIMENSION(jpi,jpj)     ::   ztri
      REAL(wp), DIMENSION(jpi,jpj)     ::   zcptrain, zcptsnw, zcptn ! Heat content per unit mass (J/kg)
      !!---------------------------------------------------------------------
      !
      zcoef_dqlw = 4._wp * emiss_i * stefan             ! local scalars
      zztmp = 1. / ( 1. - albo )
      dqla_ice(:,:,:) = 0._wp

      ! Heat content per unit mass (J/kg)
      zcptrain(:,:) = (      ptair        - rt0 ) * rcp  * tmask(:,:,1)
      zcptsnw (:,:) = ( MIN( ptair, rt0 ) - rt0 ) * rcpi * tmask(:,:,1)
      zcptn   (:,:) =        sst_m                * rcp  * tmask(:,:,1)
      !
      !                                     ! ========================== !
      DO jl = 1, jpl                        !  Loop over ice categories  !
         !                                  ! ========================== !
         DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )

            zst   = ptsu(ji,jj,jl)                                ! surface temperature of sea-ice [K]
            zsq   = q_sat( zst, pslp(ji,jj), l_ice=.TRUE. )       ! surface saturation specific humidity when ice present
            zsipt = theta_exner( zst, pslp(ji,jj) )               ! potential sea-ice surface temperature [K]  

            ! ----------------------------!
            !      I   Radiative FLUXES   !
            ! ----------------------------!
            ! Short Wave (sw)
            qsr_ice(ji,jj,jl) = zztmp * ( 1. - palb(ji,jj,jl) ) * qsr(ji,jj)

            ! Long  Wave (lw)
            zst3 = zst * zst * zst
            z_qlw(ji,jj,jl)   = emiss_i * ( pdqlw(ji,jj) - stefan * zst * zst3 ) * tmask(ji,jj,1)
            ! lw sensitivity
            z_dqlw(ji,jj,jl)  = zcoef_dqlw * zst3

            ! ----------------------------!
            !     II    Turbulent FLUXES  !
            ! ----------------------------!

            ! ... turbulent heat fluxes with Ch_ice recalculated in blk_ice_1

            ! Common term in bulk F. equations...
            zzblk = rhoa(ji,jj) * wndm_ice(ji,jj)

            ! Sensible Heat
            zztmp1 = zzblk * rCp_air * Ch_ice(ji,jj)
            z_qsb (ji,jj,jl) = zztmp1 * (zsipt - theta_zu_i(ji,jj))
            z_dqsb(ji,jj,jl) = zztmp1                        ! ==> Qsens sensitivity (Dqsb_ice/Dtn_ice)

            ! Latent Heat
            zztmp1 = zzblk * rLsub * Ce_ice(ji,jj)
            qla_ice(ji,jj,jl) = MAX( zztmp1 * (zsq - q_zu_i(ji,jj)) , 0._wp )   ! #LB: only sublimation (and not condensation) ???
            IF(qla_ice(ji,jj,jl)>0._wp) dqla_ice(ji,jj,jl) = zztmp1*dq_sat_dt_ice(zst, pslp(ji,jj)) ! ==> Qlat sensitivity  (dQlat/dT)
            !                                                                                       !#LB: dq_sat_dt_ice() in "sbc_phy.F90"
            !#LB: without this unjustified "condensation sensure":
            !qla_ice( ji,jj,jl) = zztmp1 * (zsq - q_zu_i(ji,jj))
            !dqla_ice(ji,jj,jl) = zztmp1 * dq_sat_dt_ice(zst, pslp(ji,jj)) ! ==> Qlat sensitivity  (dQlat/dT)


            ! ----------------------------!
            !     III    Total FLUXES     !
            ! ----------------------------!
            ! Downward Non Solar flux
            qns_ice (ji,jj,jl) =     z_qlw (ji,jj,jl) - z_qsb (ji,jj,jl) - qla_ice (ji,jj,jl)
            ! Total non solar heat flux sensitivity for ice
            dqns_ice(ji,jj,jl) = - ( z_dqlw(ji,jj,jl) + z_dqsb(ji,jj,jl) + dqla_ice(ji,jj,jl) ) !#LB: correct signs ????

         END_2D
         !
      END DO
      !
      tprecip(:,:) = pprec(:,:) * rn_pfac * tmask(:,:,1)  ! total precipitation [kg/m2/s]
      sprecip(:,:) = psnow(:,:) * rn_pfac * tmask(:,:,1)  ! solid precipitation [kg/m2/s]
      CALL iom_put( 'snowpre', sprecip )                  ! Snow precipitation
      CALL iom_put( 'precip' , tprecip )                  ! Total precipitation

      ! --- evaporation --- !
      z1_rLsub = 1._wp / rLsub
      evap_ice (:,:,:) = rn_efac * qla_ice (:,:,:) * z1_rLsub    ! sublimation
      devap_ice(:,:,:) = rn_efac * dqla_ice(:,:,:) * z1_rLsub    ! d(sublimation)/dT
      zevap    (:,:)   = emp(:,:) + tprecip(:,:)   ! evaporation over ocean  !LB: removed rn_efac here, correct???

      ! --- evaporation minus precipitation --- !
      zsnw(:,:) = 0._wp
      CALL ice_var_snwblow( (1.-at_i_b(:,:)), zsnw )  ! snow distribution over ice after wind blowing
      emp_oce(:,:) = ( 1._wp - at_i_b(:,:) ) * zevap(:,:) - ( tprecip(:,:) - sprecip(:,:) ) - sprecip(:,:) * (1._wp - zsnw )
      emp_ice(:,:) = SUM( a_i_b(:,:,:) * evap_ice(:,:,:), dim=3 ) - sprecip(:,:) * zsnw
      emp_tot(:,:) = emp_oce(:,:) + emp_ice(:,:)

      ! --- heat flux associated with emp --- !
      qemp_oce(:,:) = - ( 1._wp - at_i_b(:,:) ) * zevap(:,:) * zcptn(:,:)         & ! evap at sst
         &          + ( tprecip(:,:) - sprecip(:,:) )   *   zcptrain(:,:)         & ! liquid precip at Tair
         &          +   sprecip(:,:) * ( 1._wp - zsnw ) * ( zcptsnw (:,:) - rLfus ) ! solid precip at min(Tair,Tsnow)
      qemp_ice(:,:) =   sprecip(:,:) *           zsnw   * ( zcptsnw (:,:) - rLfus ) ! solid precip (only)

      ! --- total solar and non solar fluxes --- !
      qns_tot(:,:) = ( 1._wp - at_i_b(:,:) ) * qns_oce(:,:) + SUM( a_i_b(:,:,:) * qns_ice(:,:,:), dim=3 )  &
         &           + qemp_ice(:,:) + qemp_oce(:,:)
      qsr_tot(:,:) = ( 1._wp - at_i_b(:,:) ) * qsr_oce(:,:) + SUM( a_i_b(:,:,:) * qsr_ice(:,:,:), dim=3 )

      ! --- heat content of precip over ice in J/m3 (to be used in 1D-thermo) --- !
      qprec_ice(:,:) = rhos * ( zcptsnw(:,:) - rLfus )

      ! --- heat content of evap over ice in W/m2 (to be used in 1D-thermo) ---
      DO jl = 1, jpl
         qevap_ice(:,:,jl) = 0._wp ! should be -evap_ice(:,:,jl)*( ( Tice - rt0 ) * rcpi * tmask(:,:,1) )
         !                         ! But we do not have Tice => consider it at 0degC => evap=0
      END DO

      ! --- shortwave radiation transmitted thru the surface scattering layer (W/m2) --- !
      IF( nn_qtrice == 0 ) THEN
         ! formulation derived from Grenfell and Maykut (1977), where transmission rate
         !    1) depends on cloudiness
         !    2) is 0 when there is any snow
         !    3) tends to 1 for thin ice
         ztri(:,:) = 0.18 * ( 1.0 - cloud_fra(:,:) ) + 0.35 * cloud_fra(:,:)  ! surface transmission when hi>10cm
         DO jl = 1, jpl
            WHERE    ( phs(:,:,jl) <= 0._wp .AND. phi(:,:,jl) <  0.1_wp )     ! linear decrease from hi=0 to 10cm
               qtr_ice_top(:,:,jl) = qsr_ice(:,:,jl) * ( ztri(:,:) + ( 1._wp - ztri(:,:) ) * ( 1._wp - phi(:,:,jl) * 10._wp ) )
            ELSEWHERE( phs(:,:,jl) <= 0._wp .AND. phi(:,:,jl) >= 0.1_wp )     ! constant (ztri) when hi>10cm
               qtr_ice_top(:,:,jl) = qsr_ice(:,:,jl) * ztri(:,:)
            ELSEWHERE                                                         ! zero when hs>0
               qtr_ice_top(:,:,jl) = 0._wp
            END WHERE
         ENDDO
      ELSEIF( nn_qtrice == 1 ) THEN
         ! formulation is derived from the thesis of M. Lebrun (2019).
         !    It represents the best fit using several sets of observations
         !    It comes with snow conductivities adapted to freezing/melting conditions (see icethd_zdf_bl99.F90)
         qtr_ice_top(:,:,:) = 0.3_wp * qsr_ice(:,:,:)
      ENDIF
      !
      IF( iom_use('evap_ao_cea') .OR. iom_use('hflx_evap_cea') ) THEN
         CALL iom_put( 'evap_ao_cea'  , zevap(:,:) * ( 1._wp - at_i_b(:,:) ) * tmask(:,:,1)              )   ! ice-free oce evap (cell average)
         CALL iom_put( 'hflx_evap_cea', zevap(:,:) * ( 1._wp - at_i_b(:,:) ) * tmask(:,:,1) * zcptn(:,:) )   ! heat flux from evap (cell average)
      ENDIF
      IF( iom_use('rain') .OR. iom_use('rain_ao_cea') .OR. iom_use('hflx_rain_cea') ) THEN
         CALL iom_put( 'rain'         ,   tprecip(:,:) - sprecip(:,:)                             )          ! liquid precipitation 
         CALL iom_put( 'rain_ao_cea'  , ( tprecip(:,:) - sprecip(:,:) ) * ( 1._wp - at_i_b(:,:) ) )          ! liquid precipitation over ocean (cell average)
         CALL iom_put( 'hflx_rain_cea', ( tprecip(:,:) - sprecip(:,:) ) * zcptrain(:,:) )                    ! heat flux from rain (cell average)
      ENDIF
      IF(  iom_use('snow_ao_cea')   .OR. iom_use('snow_ai_cea')      .OR. &
         & iom_use('hflx_snow_cea') .OR. iom_use('hflx_snow_ao_cea') .OR. iom_use('hflx_snow_ai_cea')  )  THEN
         CALL iom_put( 'snow_ao_cea'     , sprecip(:,:)                            * ( 1._wp - zsnw(:,:) ) ) ! Snow over ice-free ocean  (cell average)
         CALL iom_put( 'snow_ai_cea'     , sprecip(:,:)                            *           zsnw(:,:)   ) ! Snow over sea-ice         (cell average)
         CALL iom_put( 'hflx_snow_cea'   , sprecip(:,:) * ( zcptsnw(:,:) - rLfus ) )                         ! heat flux from snow (cell average)
         CALL iom_put( 'hflx_snow_ao_cea', sprecip(:,:) * ( zcptsnw(:,:) - rLfus ) * ( 1._wp - zsnw(:,:) ) ) ! heat flux from snow (over ocean)
         CALL iom_put( 'hflx_snow_ai_cea', sprecip(:,:) * ( zcptsnw(:,:) - rLfus ) *           zsnw(:,:)   ) ! heat flux from snow (over ice)
      ENDIF
      IF( iom_use('hflx_prec_cea') ) THEN                                                                    ! heat flux from precip (cell average)
         CALL iom_put('hflx_prec_cea' ,    sprecip(:,:)                  * ( zcptsnw (:,:) - rLfus )  &
            &                          + ( tprecip(:,:) - sprecip(:,:) ) *   zcptrain(:,:) )
      ENDIF
      IF( iom_use('subl_ai_cea') .OR. iom_use('hflx_subl_cea') ) THEN
         CALL iom_put( 'subl_ai_cea'  , SUM( a_i_b(:,:,:) *  evap_ice(:,:,:), dim=3 ) * tmask(:,:,1) ) ! Sublimation over sea-ice (cell average)
         CALL iom_put( 'hflx_subl_cea', SUM( a_i_b(:,:,:) * qevap_ice(:,:,:), dim=3 ) * tmask(:,:,1) ) ! Heat flux from sublimation (cell average)
      ENDIF
      !
      IF(sn_cfctl%l_prtctl) THEN
         ALLOCATE(zmsk(jpi,jpj,jpl))
         DO jl = 1, jpl
            zmsk(:,:,jl) = tmask(:,:,1)
         END DO
         CALL prt_ctl(tab3d_1=qla_ice , clinfo1=' blk_ice: qla_ice  : ', mask1=zmsk,   &
            &         tab3d_2=z_qsb   , clinfo2=' z_qsb    : '         , mask2=zmsk, kdim=jpl)
         CALL prt_ctl(tab3d_1=z_qlw   , clinfo1=' blk_ice: z_qlw    : ', mask1=zmsk,   &
            &         tab3d_2=dqla_ice, clinfo2=' dqla_ice : '         , mask2=zmsk, kdim=jpl)
         CALL prt_ctl(tab3d_1=z_dqsb  , clinfo1=' blk_ice: z_dqsb   : ', mask1=zmsk,   &
            &         tab3d_2=z_dqlw  , clinfo2=' z_dqlw   : '         , mask2=zmsk, kdim=jpl)
         CALL prt_ctl(tab3d_1=dqns_ice, clinfo1=' blk_ice: dqns_ice : ', mask1=zmsk,   &
            &         tab3d_2=qsr_ice , clinfo2=' qsr_ice  : '         , mask2=zmsk, kdim=jpl)
         CALL prt_ctl(tab3d_1=ptsu    , clinfo1=' blk_ice: ptsu     : ', mask1=zmsk,   &
            &         tab3d_2=qns_ice , clinfo2=' qns_ice  : '         , mask2=zmsk, kdim=jpl)
         CALL prt_ctl(tab2d_1=tprecip , clinfo1=' blk_ice: tprecip  : ', mask1=tmask,   &
            &         tab2d_2=sprecip , clinfo2=' sprecip  : '         , mask2=tmask         )
         DEALLOCATE(zmsk)
      ENDIF

      !#LB:
      ! air-ice heat flux components that are not written from ice_stp()@icestp.F90:
      IF( iom_use('qla_ice') )  CALL iom_put( 'qla_ice', SUM( - qla_ice * a_i_b, dim=3 ) ) !#LB: sign consistent with what's done for ocean
      IF( iom_use('qsb_ice') )  CALL iom_put( 'qsb_ice', SUM( -   z_qsb * a_i_b, dim=3 ) ) !#LB:     ==> negative => loss of heat for sea-ice
      IF( iom_use('qlw_ice') )  CALL iom_put( 'qlw_ice', SUM(     z_qlw * a_i_b, dim=3 ) )
      !#LB.

   END SUBROUTINE blk_ice_2


   SUBROUTINE blk_ice_qcn( ld_virtual_itd, ptsu, ptb, phs, phi )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE blk_ice_qcn  ***
      !!
      !! ** Purpose :   Compute surface temperature and snow/ice conduction flux
      !!                to force sea ice / snow thermodynamics
      !!                in the case conduction flux is emulated
      !!
      !! ** Method  :   compute surface energy balance assuming neglecting heat storage
      !!                following the 0-layer Semtner (1976) approach
      !!
      !! ** Outputs : - ptsu    : sea-ice / snow surface temperature (K)
      !!              - qcn_ice : surface inner conduction flux (W/m2)
      !!
      !!---------------------------------------------------------------------
      LOGICAL                   , INTENT(in   ) ::   ld_virtual_itd  ! single-category option
      REAL(wp), DIMENSION(:,:,:), INTENT(inout) ::   ptsu            ! sea ice / snow surface temperature
      REAL(wp), DIMENSION(:,:)  , INTENT(in   ) ::   ptb             ! sea ice base temperature
      REAL(wp), DIMENSION(:,:,:), INTENT(in   ) ::   phs             ! snow thickness
      REAL(wp), DIMENSION(:,:,:), INTENT(in   ) ::   phi             ! sea ice thickness
      !
      INTEGER , PARAMETER ::   nit = 10                  ! number of iterations
      REAL(wp), PARAMETER ::   zepsilon = 0.1_wp         ! characteristic thickness for enhanced conduction
      !
      INTEGER  ::   ji, jj, jl           ! dummy loop indices
      INTEGER  ::   iter                 ! local integer
      REAL(wp) ::   zfac, zfac2, zfac3   ! local scalars
      REAL(wp) ::   zkeff_h, ztsu, ztsu0 !
      REAL(wp) ::   zqc, zqnet           !
      REAL(wp) ::   zhe, zqa0            !
      REAL(wp), DIMENSION(jpi,jpj,jpl) ::   zgfac   ! enhanced conduction factor
      !!---------------------------------------------------------------------

      ! -------------------------------------!
      !      I   Enhanced conduction factor  !
      ! -------------------------------------!
      ! Emulates the enhancement of conduction by unresolved thin ice (ld_virtual_itd = T)
      ! Fichefet and Morales Maqueda, JGR 1997
      !
      zgfac(:,:,:) = 1._wp

      IF( ld_virtual_itd ) THEN
         !
         zfac  = 1._wp /  ( rn_cnd_s + rcnd_i )
         zfac2 = EXP(1._wp) * 0.5_wp * zepsilon
         zfac3 = 2._wp / zepsilon
         !
         DO jl = 1, jpl
            DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
               zhe = ( rn_cnd_s * phi(ji,jj,jl) + rcnd_i * phs(ji,jj,jl) ) * zfac                            ! Effective thickness
               IF( zhe >=  zfac2 )   zgfac(ji,jj,jl) = MIN( 2._wp, 0.5_wp * ( 1._wp + LOG( zhe * zfac3 ) ) ) ! Enhanced conduction factor
            END_2D
         END DO
         !
      ENDIF

      ! -------------------------------------------------------------!
      !      II   Surface temperature and conduction flux            !
      ! -------------------------------------------------------------!
      !
      zfac = rcnd_i * rn_cnd_s
      !
      DO jl = 1, jpl
         DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
            !
            zkeff_h = zfac * zgfac(ji,jj,jl) / &                                    ! Effective conductivity of the snow-ice system divided by thickness
               &      ( rcnd_i * phs(ji,jj,jl) + rn_cnd_s * MAX( 0.01, phi(ji,jj,jl) ) )
            ztsu    = ptsu(ji,jj,jl)                                                ! Store current iteration temperature
            ztsu0   = ptsu(ji,jj,jl)                                                ! Store initial surface temperature
            zqa0    = qsr_ice(ji,jj,jl) - qtr_ice_top(ji,jj,jl) + qns_ice(ji,jj,jl) ! Net initial atmospheric heat flux
            !
            DO iter = 1, nit     ! --- Iterative loop
               zqc   = zkeff_h * ( ztsu - ptb(ji,jj) )                              ! Conduction heat flux through snow-ice system (>0 downwards)
               zqnet = zqa0 + dqns_ice(ji,jj,jl) * ( ztsu - ptsu(ji,jj,jl) ) - zqc  ! Surface energy budget
               ztsu  = ztsu - zqnet / ( dqns_ice(ji,jj,jl) - zkeff_h )              ! Temperature update
            END DO
            !
            ptsu   (ji,jj,jl) = MIN( rt0, ztsu )
            qcn_ice(ji,jj,jl) = zkeff_h * ( ptsu(ji,jj,jl) - ptb(ji,jj) )
            qns_ice(ji,jj,jl) = qns_ice(ji,jj,jl) + dqns_ice(ji,jj,jl) * ( ptsu(ji,jj,jl) - ztsu0 )
            qml_ice(ji,jj,jl) = ( qsr_ice(ji,jj,jl) - qtr_ice_top(ji,jj,jl) + qns_ice(ji,jj,jl) - qcn_ice(ji,jj,jl) )  &
               &   * MAX( 0._wp , SIGN( 1._wp, ptsu(ji,jj,jl) - rt0 ) )

            ! --- Diagnose the heat loss due to changing non-solar flux (as in icethd_zdf_bl99) --- !
            hfx_err_dif(ji,jj) = hfx_err_dif(ji,jj) - ( dqns_ice(ji,jj,jl) * ( ptsu(ji,jj,jl) - ztsu0 ) ) * a_i_b(ji,jj,jl)

         END_2D
         !
      END DO
      !
   END SUBROUTINE blk_ice_qcn

#endif

   !!======================================================================
END MODULE sbcblk
