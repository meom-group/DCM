MODULE icedia
   !!======================================================================
   !!                       ***  MODULE icedia  ***
   !!  Sea-Ice:   global budgets
   !!======================================================================
   !! History :  3.4  !  2012-10  (C. Rousset)       original code
   !!            4.0  !  2018     (many people)      SI3 [aka Sea Ice cube]
   !!----------------------------------------------------------------------
#if defined key_si3
   !!----------------------------------------------------------------------
   !!   'key_si3'                                       SI3 sea-ice model
   !!----------------------------------------------------------------------
   !!    ice_dia      : diagnostic of the sea-ice global heat content, salt content and volume conservation
   !!    ice_dia_init : initialization of budget calculation
   !!    ice_dia_rst  : read/write budgets restart
   !!----------------------------------------------------------------------
   USE dom_oce        ! ocean domain
   USE phycst         ! physical constant
   USE daymod         ! model calendar
   USE sbc_oce , ONLY : sfx, nn_fsbc   ! surface boundary condition: ocean fields
   USE ice            ! sea-ice: variables
   USE icerst         ! sea-ice: restart
   !
   USE in_out_manager ! I/O manager
   USE iom            ! I/O manager library
   USE lib_mpp        ! MPP library
   USE lib_fortran    ! fortran utilities (glob_sum + no signed zero)
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   ice_dia        ! called by icestp.F90
   PUBLIC   ice_dia_init   ! called in icestp.F90

   REAL(wp), SAVE ::   r1_area  ! inverse of the ocean area
   REAL(wp), DIMENSION(:,:), ALLOCATABLE ::   vol_loc_ini, sal_loc_ini, tem_loc_ini                    ! initial volume, salt and heat contents
   REAL(wp)                              ::   frc_sal, frc_voltop, frc_volbot, frc_temtop, frc_tembot  ! global forcing trends

   !!----------------------------------------------------------------------
   !! NEMO/ICE 4.0 , NEMO Consortium (2018)
   !! $Id: icedia.F90 15048 2021-06-23 16:02:14Z clem $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION ice_dia_alloc()
      !!---------------------------------------------------------------------!
      !!                ***  ROUTINE ice_dia_alloc ***
      !!---------------------------------------------------------------------!
      ALLOCATE( vol_loc_ini(jpi,jpj), sal_loc_ini(jpi,jpj), tem_loc_ini(jpi,jpj), STAT=ice_dia_alloc )

      CALL mpp_sum ( 'icedia', ice_dia_alloc )
      IF( ice_dia_alloc /= 0 )   CALL ctl_stop( 'STOP',  'ice_dia_alloc: failed to allocate arrays'  )
      !
   END FUNCTION ice_dia_alloc

   SUBROUTINE ice_dia( kt )
      !!---------------------------------------------------------------------------
      !!                  ***  ROUTINE ice_dia  ***
      !!
      !! ** Purpose:   Compute the sea-ice global heat content, salt content
      !!             and volume conservation
      !!---------------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time step
      !!
      REAL(wp), DIMENSION(jpi,jpj,16) ::   ztmp
      REAL(wp), DIMENSION(16)         ::   zbg          
      !!---------------------------------------------------------------------------
      IF( ln_timing )   CALL timing_start('ice_dia')

      IF( kt == nit000 .AND. lwp ) THEN
         WRITE(numout,*)
         WRITE(numout,*)'icedia: output ice diagnostics (integrated over the domain)'
         WRITE(numout,*)'~~~~~~'
      ENDIF

      IF( kt == nit000 ) THEN
         r1_area = 1._wp / glob_sum( 'icedia', e1e2t(:,:) )
      ENDIF

      ztmp(:,:,:) = 0._wp ! should be better coded
      
      ! ---------------------------!
      ! 1 - Trends due to forcing  !
      ! ---------------------------!
      ! they must be kept outside an IF(iom_use) because of the call to dia_rst below
      ztmp(:,:,1) = - ( wfx_ice(:,:) + wfx_snw(:,:) + wfx_err_sub(:,:) ) * e1e2t(:,:) ! freshwater flux ice/snow-ocean
      ztmp(:,:,2) = - ( wfx_sub(:,:) + wfx_spr(:,:) )                    * e1e2t(:,:) ! freshwater flux ice/snow-atm
      ztmp(:,:,3) = -   sfx    (:,:)                                     * e1e2t(:,:) ! salt fluxes ice/snow-ocean
      ztmp(:,:,4) =   qt_atm_oi(:,:)                                     * e1e2t(:,:) ! heat on top of ice-ocean
      ztmp(:,:,5) =   qt_oce_ai(:,:)                                     * e1e2t(:,:) ! heat on top of ocean (and below ice)
      
      ! ----------------------- !
      ! 2 -  Contents           !
      ! ----------------------- !
      IF( iom_use('ibgvol_tot' ) )   ztmp(:,:,6 ) = vt_i (:,:) * e1e2t(:,:) ! ice volume
      IF( iom_use('sbgvol_tot' ) )   ztmp(:,:,7 ) = vt_s (:,:) * e1e2t(:,:) ! snow volume
      IF( iom_use('ibgarea_tot') )   ztmp(:,:,8 ) = at_i (:,:) * e1e2t(:,:) ! area
      IF( iom_use('ibgsalt_tot') )   ztmp(:,:,9 ) = st_i (:,:) * e1e2t(:,:) ! salt content
      IF( iom_use('ibgheat_tot') )   ztmp(:,:,10) = et_i (:,:) * e1e2t(:,:) ! heat content
      IF( iom_use('sbgheat_tot') )   ztmp(:,:,11) = et_s (:,:) * e1e2t(:,:) ! heat content
      IF( iom_use('ipbgvol_tot') )   ztmp(:,:,12) = vt_ip(:,:) * e1e2t(:,:) ! ice pond volume
      IF( iom_use('ilbgvol_tot') )   ztmp(:,:,13) = vt_il(:,:) * e1e2t(:,:) ! ice pond lid volume

      ! ---------------------------------- !
      ! 3 -  Content variations and drifts !
      ! ---------------------------------- !
      IF( iom_use('ibgvolume') ) ztmp(:,:,14) = ( rhoi*vt_i(:,:) + rhos*vt_s(:,:) - vol_loc_ini(:,:) ) * e1e2t(:,:) ! freshwater trend
      IF( iom_use('ibgsaltco') ) ztmp(:,:,15) = ( rhoi*st_i(:,:)                  - sal_loc_ini(:,:) ) * e1e2t(:,:) ! salt content trend
      IF( iom_use('ibgheatco') .OR. iom_use('ibgheatfx') ) &
         &                       ztmp(:,:,16) = ( et_i(:,:) + et_s(:,:)           - tem_loc_ini(:,:) ) * e1e2t(:,:) ! heat content trend
      
      ! global sum
      zbg(1:16) = glob_sum_vec( 'icedia', ztmp(:,:,1:16) )

      ! change units for trends
      zbg(1) = zbg(1) * r1_rho0 * 1.e-9  * rDt_ice ! freshwater flux ice/snow-ocean (km3)
      zbg(2) = zbg(2) * r1_rho0 * 1.e-9  * rDt_ice ! freshwater flux ice/snow-atm (km3)
      zbg(3) = zbg(3) * r1_rho0 * 1.e-9  * rDt_ice ! salt fluxes ice/snow-ocean (km3*pss)
      zbg(4) = zbg(4)           * 1.e-20 * rDt_ice ! heat on top of ice-ocean (1.e20 J)
      zbg(5) = zbg(5)           * 1.e-20 * rDt_ice ! heat on top of ocean (and below ice) (1.e20 J)
      ! cumulative sum
      frc_voltop  = frc_voltop  + zbg(1)
      frc_volbot  = frc_volbot  + zbg(2)
      frc_sal     = frc_sal     + zbg(3)
      frc_temtop  = frc_temtop  + zbg(4)
      frc_tembot  = frc_tembot  + zbg(5)

      ! change units for contents
      zbg(6)  = zbg(6)  * 1.e-9  ! ice volume (km3)
      zbg(7)  = zbg(7)  * 1.e-9  ! snw volume (km3)
      zbg(8)  = zbg(8)  * 1.e-6  ! ice area (km2)
      zbg(9)  = zbg(9)  * 1.e-9  ! salt content (km3*pss)
      zbg(10) = zbg(10) * 1.e-20 ! ice heat content (1.e20 J)
      zbg(11) = zbg(11) * 1.e-20 ! snw heat content (1.e20 J)
      zbg(12) = zbg(12) * 1.e-9  ! pnd volume (km3)
      zbg(13) = zbg(13) * 1.e-9  ! pnd lid volume (km3)

      ! change units for trends
      zbg(14) = zbg(14) * r1_rho0 * 1.e-9  ! freshwater trend (km3)
      zbg(15) = zbg(15) * r1_rho0 * 1.e-9  ! salt content trend (km3*pss)
      zbg(16) = zbg(16)           * 1.e-20 ! heat content trend (1.e20 J)
      ! difference
      zbg(14) = zbg(14) - ( frc_voltop + frc_volbot )
      zbg(15) = zbg(15) -   frc_sal
      zbg(16) = zbg(16) - ( frc_tembot - frc_temtop )

      ! outputs
      CALL iom_put( 'ibgfrcvoltop' , frc_voltop )   ! vol  forcing ice/snw-atm          (km3 equivalent ocean water)
      CALL iom_put( 'ibgfrcvolbot' , frc_volbot )   ! vol  forcing ice/snw-ocean        (km3 equivalent ocean water)
      CALL iom_put( 'ibgfrcsal'    , frc_sal    )   ! sal  forcing                      (psu*km3 equivalent ocean water)
      CALL iom_put( 'ibgfrctemtop' , frc_temtop )   ! heat on top of ice/snw/ocean      (1.e20 J)
      CALL iom_put( 'ibgfrctembot' , frc_tembot )   ! heat on top of ocean(below ice)   (1.e20 J)
      CALL iom_put( 'ibgfrchfxtop' , frc_temtop * r1_area * 1.e-20 * kt*rn_Dt ) ! heat on top of ice/snw/ocean      (W/m2)
      CALL iom_put( 'ibgfrchfxbot' , frc_tembot * r1_area * 1.e-20 * kt*rn_Dt ) ! heat on top of ocean(below ice)   (W/m2)

      CALL iom_put( 'ibgvol_tot'  , zbg(6)  )
      CALL iom_put( 'sbgvol_tot'  , zbg(7)  )
      CALL iom_put( 'ibgarea_tot' , zbg(8)  )
      CALL iom_put( 'ibgsalt_tot' , zbg(9)  )
      CALL iom_put( 'ibgheat_tot' , zbg(10) )
      CALL iom_put( 'sbgheat_tot' , zbg(11) )
      CALL iom_put( 'ipbgvol_tot' , zbg(12) )
      CALL iom_put( 'ilbgvol_tot' , zbg(13) )
     
      CALL iom_put( 'ibgvolume' , zbg(14) )   ! ice/snow volume  drift            (km3 equivalent ocean water)
      CALL iom_put( 'ibgsaltco' , zbg(15) )   ! ice salt content drift            (psu*km3 equivalent ocean water)
      CALL iom_put( 'ibgheatco' , zbg(16) )   ! ice/snow heat content drift       (1.e20 J)
      !
      ! restarts
      IF( lrst_ice )   CALL ice_dia_rst( 'WRITE', kt_ice )
      !
      IF( ln_timing )   CALL timing_stop('ice_dia')
      !
   END SUBROUTINE ice_dia


   SUBROUTINE ice_dia_init
      !!---------------------------------------------------------------------------
      !!                  ***  ROUTINE ice_dia_init  ***
      !!
      !! ** Purpose: Initialization for the heat salt volume budgets
      !!
      !! ** Method : Compute initial heat content, salt content and volume
      !!
      !! ** Action : - Compute initial heat content, salt content and volume
      !!             - Initialize forcing trends
      !!             - Compute coefficients for conversion
      !!---------------------------------------------------------------------------
      INTEGER            ::   ios, ierror   ! local integer
      !!
      NAMELIST/namdia/ ln_icediachk, rn_icechk_cel, rn_icechk_glo, ln_icediahsb, ln_icectl, iiceprt, jiceprt
      !!----------------------------------------------------------------------
      !
      READ  ( numnam_ice_ref, namdia, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namdia in reference namelist' )
      READ  ( numnam_ice_cfg, namdia, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namdia in configuration namelist' )
      IF(lwm) WRITE ( numoni, namdia )
      !
      IF(lwp) THEN                  ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'ice_dia_init: ice diagnostics'
         WRITE(numout,*) ' ~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namdia:'
         WRITE(numout,*) '      Diagnose online heat/mass/salt conservation ln_icediachk  = ', ln_icediachk
         WRITE(numout,*) '         threshold for conservation (gridcell)    rn_icechk_cel = ', rn_icechk_cel
         WRITE(numout,*) '         threshold for conservation (global)      rn_icechk_glo = ', rn_icechk_glo
         WRITE(numout,*) '      Output          heat/mass/salt budget       ln_icediahsb  = ', ln_icediahsb
         WRITE(numout,*) '      control prints for a given grid point       ln_icectl     = ', ln_icectl
         WRITE(numout,*) '         chosen grid point position          (iiceprt,jiceprt)  = (', iiceprt,',', jiceprt,')'
      ENDIF
      !
      IF( ln_icediahsb ) THEN
         IF( ice_dia_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'ice_dia_init : unable to allocate arrays' )   ! allocate tke arrays
         CALL ice_dia_rst( 'READ' )   ! read or initialize all required files
      ENDIF
      !
   END SUBROUTINE ice_dia_init


   SUBROUTINE ice_dia_rst( cdrw, kt )
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE icedia_rst  ***
      !!
      !! ** Purpose :   Read or write DIA file in restart file
      !!
      !! ** Method  :   use of IOM library
      !!----------------------------------------------------------------------
      CHARACTER(len=*) , INTENT(in) ::   cdrw   ! "READ"/"WRITE" flag
      INTEGER, OPTIONAL, INTENT(in) ::   kt     ! ice time-step
      !
      INTEGER  ::   iter    ! local integer
      REAL(wp) ::   ziter   ! local scalar
      !!----------------------------------------------------------------------
      !
      IF( TRIM(cdrw) == 'READ' ) THEN        ! Read/initialise
         IF( ln_rstart ) THEN                   !* Read the restart file
            !
            CALL iom_get( numrir, 'kt_ice' , ziter )
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'ice_dia_rst read at time step = ', ziter
            IF(lwp) WRITE(numout,*) '~~~~~~~~~~'
#if defined key_drakkar
            IF( iom_varid( numrir, 'frc_voltop', ldstop = .FALSE. ) <= 0 )   THEN
              ! set trends to 0 if not found in restart file
              frc_voltop  = 0._wp                                          
              frc_volbot  = 0._wp                                          
              frc_temtop  = 0._wp                                                 
              frc_tembot  = 0._wp                                                 
              frc_sal     = 0._wp                                                 
              ! record initial ice volume, salt and temp
              vol_loc_ini(:,:) = rhoi * vt_i(:,:) + rhos * vt_s(:,:)  ! ice/snow volume (kg/m2)
              tem_loc_ini(:,:) = et_i(:,:) + et_s(:,:)                ! ice/snow heat content (J)
              sal_loc_ini(:,:) = rhoi * st_i(:,:)                     ! ice salt content (pss*kg/m2)
            ELSE
#else
            CALL iom_get( numrir, 'frc_voltop' , frc_voltop  )
            CALL iom_get( numrir, 'frc_volbot' , frc_volbot  )
            CALL iom_get( numrir, 'frc_temtop' , frc_temtop  )
            CALL iom_get( numrir, 'frc_tembot' , frc_tembot  )
            CALL iom_get( numrir, 'frc_sal'    , frc_sal     )
            CALL iom_get( numrir, jpdom_auto, 'vol_loc_ini', vol_loc_ini )
            CALL iom_get( numrir, jpdom_auto, 'tem_loc_ini', tem_loc_ini )
            CALL iom_get( numrir, jpdom_auto, 'sal_loc_ini', sal_loc_ini )
#endif
#if defined key_drakkar
            ENDIF
#endif
         ELSE
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) ' ice_dia at initial state '
            IF(lwp) WRITE(numout,*) '~~~~~~~'
            ! set trends to 0
            frc_voltop  = 0._wp
            frc_volbot  = 0._wp
            frc_temtop  = 0._wp
            frc_tembot  = 0._wp
            frc_sal     = 0._wp
            ! record initial ice volume, salt and temp
            vol_loc_ini(:,:) = rhoi * vt_i(:,:) + rhos * vt_s(:,:)  ! ice/snow volume (kg/m2)
            tem_loc_ini(:,:) = et_i(:,:) + et_s(:,:)                ! ice/snow heat content (J)
            sal_loc_ini(:,:) = rhoi * st_i(:,:)                     ! ice salt content (pss*kg/m2)
         ENDIF
         !
      ELSEIF( TRIM(cdrw) == 'WRITE' ) THEN   ! Create restart file
         !                                   ! -------------------
         iter = kt + nn_fsbc - 1   ! ice restarts are written at kt == nitrst - nn_fsbc + 1
         !
         IF( iter == nitrst ) THEN
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'ice_dia_rst write at time step = ', kt
            IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
         ENDIF
         !
         ! Write in numriw (if iter == nitrst)
         ! ------------------
         CALL iom_rstput( iter, nitrst, numriw, 'frc_voltop' , frc_voltop )
         CALL iom_rstput( iter, nitrst, numriw, 'frc_volbot' , frc_volbot )
         CALL iom_rstput( iter, nitrst, numriw, 'frc_temtop' , frc_temtop )
         CALL iom_rstput( iter, nitrst, numriw, 'frc_tembot' , frc_tembot )
         CALL iom_rstput( iter, nitrst, numriw, 'frc_sal'    , frc_sal    )
         CALL iom_rstput( iter, nitrst, numriw, 'vol_loc_ini', vol_loc_ini )
         CALL iom_rstput( iter, nitrst, numriw, 'tem_loc_ini', tem_loc_ini )
         CALL iom_rstput( iter, nitrst, numriw, 'sal_loc_ini', sal_loc_ini )
         !
      ENDIF
      !
   END SUBROUTINE ice_dia_rst

#else
   !!----------------------------------------------------------------------
   !!   Default option :         Empty module         NO SI3 sea-ice model
   !!----------------------------------------------------------------------
#endif

   !!======================================================================
END MODULE icedia
