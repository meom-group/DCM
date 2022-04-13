MODULE sbcfwb
   !!======================================================================
   !!                       ***  MODULE  sbcfwb  ***
   !! Ocean fluxes   : domain averaged freshwater budget
   !!======================================================================
   !! History :  OPA  ! 2001-02  (E. Durand)  Original code
   !!   NEMO     1.0  ! 2002-06  (G. Madec)  F90: Free form and module
   !!            3.0  ! 2006-08  (G. Madec)  Surface module
   !!            3.2  ! 2009-07  (C. Talandier) emp mean s spread over erp area 
   !!            3.6  ! 2014-11  (P. Mathiot  ) add ice shelf melting
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   sbc_fwb       : freshwater budget for global ocean configurations (free surface & forced mode)
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers
   USE dom_oce        ! ocean space and time domain
   USE sbc_oce        ! surface ocean boundary condition
   USE isf_oce , ONLY : fwfisf_cav, fwfisf_par                    ! ice shelf melting contribution
   USE sbc_ice , ONLY : snwice_mass, snwice_mass_b, snwice_fmass
   USE phycst         ! physical constants
   USE sbcrnf         ! ocean runoffs
   USE sbcssr         ! Sea-Surface damping terms
   !
   USE in_out_manager ! I/O manager
   USE iom            ! IOM
   USE lib_mpp        ! distribued memory computing library
   USE timing         ! Timing
   USE lbclnk         ! ocean lateral boundary conditions
   USE lib_fortran    ! 

   IMPLICIT NONE
   PRIVATE

   PUBLIC   sbc_fwb    ! routine called by step

   REAL(wp) ::   rn_fwb0   ! initial freshwater adjustment flux [kg/m2/s] (nn_fwb = 2 only)
   REAL(wp) ::   a_fwb     ! annual domain averaged freshwater budget from the previous year
   REAL(wp) ::   a_fwb_b   ! annual domain averaged freshwater budget from the year before or at initial state
   REAL(wp) ::   a_fwb_ini ! initial domain averaged freshwater budget
   REAL(wp) ::   area      ! global mean ocean surface (interior domain)

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: sbcfwb.F90 15439 2021-10-22 17:53:09Z clem $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE sbc_fwb( kt, kn_fwb, kn_fsbc, Kmm )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE sbc_fwb  ***
      !!
      !! ** Purpose :   Control the mean sea surface drift
      !!
      !! ** Method  :   several ways  depending on kn_fwb
      !!                =0 no control 
      !!                =1 global mean of emp set to zero at each nn_fsbc time step
      !!                =2 annual global mean corrected from previous year
      !!                =3 global mean of emp set to zero at each nn_fsbc time step
      !!                   & spread out over erp area depending its sign
      !! Note: if sea ice is embedded it is taken into account when computing the budget 
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt       ! ocean time-step index
      INTEGER, INTENT( in ) ::   kn_fsbc  ! 
      INTEGER, INTENT( in ) ::   kn_fwb   ! ocean time-step index
      INTEGER, INTENT( in ) ::   Kmm      ! ocean time level index
      !
      INTEGER  ::   ios, inum, ikty       ! local integers
      REAL(wp) ::   z_fwf, z_fwf_nsrf, zsum_fwf, zsum_erp                ! local scalars
      REAL(wp) ::   zsurf_neg, zsurf_pos, zsurf_tospread, zcoef          !   -      -
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) ::   ztmsk_neg, ztmsk_pos, z_wgt ! 2D workspaces
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) ::   ztmsk_tospread, zerp_cor    !   -      -
      REAL(wp)   ,DIMENSION(1) ::   z_fwfprv  
      COMPLEX(dp),DIMENSION(1) ::   y_fwfnow  
      !
      NAMELIST/namsbc_fwb/rn_fwb0
      !!----------------------------------------------------------------------
      !
      IF( kt == nit000 ) THEN
         READ( numnam_ref, namsbc_fwb, IOSTAT = ios, ERR = 901 )
901      IF( ios /= 0 ) CALL ctl_nam( ios, 'namsbc_fwb in reference namelist'     )
         READ( numnam_cfg, namsbc_fwb, IOSTAT = ios, ERR = 902 )
902      IF( ios >  0 ) CALL ctl_nam( ios, 'namsbc_fwb in configuration namelist' )
         IF(lwm) WRITE( numond, namsbc_fwb )
         IF(lwp) THEN
            WRITE(numout,*)
            WRITE(numout,*) 'sbc_fwb : FreshWater Budget correction'
            WRITE(numout,*) '~~~~~~~'
            IF( kn_fwb == 1 )   WRITE(numout,*) '          instantaneously set to zero'
            IF( kn_fwb == 3 )   WRITE(numout,*) '          fwf set to zero and spread out over erp area'
            IF( kn_fwb == 2 ) THEN
               WRITE(numout,*) '          adjusted from previous year budget'
               WRITE(numout,*)
               WRITE(numout,*) '   Namelist namsbc_fwb'
               WRITE(numout,*) '      Initial freshwater adjustment flux [kg/m2/s] = ', rn_fwb0
            END IF
         ENDIF
         !
         IF( kn_fwb == 3 .AND. nn_sssr /= 2 )   CALL ctl_stop( 'sbc_fwb: nn_fwb = 3 requires nn_sssr = 2, we stop ' )
         IF( kn_fwb == 3 .AND. ln_isfcav    )   CALL ctl_stop( 'sbc_fwb: nn_fwb = 3 with ln_isfcav = .TRUE. not working, we stop ' )
         !
         area = glob_sum( 'sbcfwb', e1e2t(:,:) * tmask(:,:,1))           ! interior global domain surface
         ! isf cavities are excluded because it can feedback to the melting with generation of inhibition of plumes
         ! and in case of no melt, it can generate HSSW.
         !
         IF( nn_ice == 0 ) THEN
            snwice_mass_b(:,:) = 0.e0               ! no sea-ice model is being used : no snow+ice mass
            snwice_mass  (:,:) = 0.e0
            snwice_fmass (:,:) = 0.e0
         ENDIF
         !
      ENDIF

      SELECT CASE ( kn_fwb )
      !
      CASE ( 1 )                             !==  global mean fwf set to zero  ==!
         !
         IF( MOD( kt-1, kn_fsbc ) == 0 ) THEN
            y_fwfnow(1) = local_sum( e1e2t(:,:) * ( emp(:,:) - rnf(:,:) - fwfisf_cav(:,:) - fwfisf_par(:,:) - snwice_fmass(:,:) ) )
            CALL mpp_delay_sum( 'sbcfwb', 'fwb', y_fwfnow(:), z_fwfprv(:), kt == nitend - nn_fsbc + 1 )
            z_fwfprv(1) = z_fwfprv(1) / area
            zcoef = z_fwfprv(1) * rcp
            emp(:,:) = emp(:,:) - z_fwfprv(1)        * tmask(:,:,1)
            qns(:,:) = qns(:,:) + zcoef * sst_m(:,:) * tmask(:,:,1) ! account for change to the heat budget due to fw correction
            ! outputs
            IF( iom_use('hflx_fwb_cea') )  CALL iom_put( 'hflx_fwb_cea', zcoef * sst_m(:,:) * tmask(:,:,1) )
            IF( iom_use('vflx_fwb_cea') )  CALL iom_put( 'vflx_fwb_cea', z_fwfprv(1)        * tmask(:,:,1) )
         ENDIF
#if defined key_drakkar
          CALL iom_put( "fwprv", z_fwfprv(1))
#endif
         !
      CASE ( 2 )                             !==  fw adjustment based on fw budget at the end of the previous year  ==!
         !                                                simulation is supposed to start 1st of January
         IF( kt == nit000 ) THEN                                                                 ! initialisation
            !                                                                                    ! set the fw adjustment (a_fwb)
            IF ( ln_rstart .AND. iom_varid( numror, 'a_fwb_b', ldstop = .FALSE. ) > 0     &      !    as read from restart file
               &           .AND. iom_varid( numror, 'a_fwb',   ldstop = .FALSE. ) > 0 ) THEN
               IF(lwp)   WRITE(numout,*) 'sbc_fwb : reading freshwater-budget from restart file'
               CALL iom_get( numror, 'a_fwb_b', a_fwb_b )
               CALL iom_get( numror, 'a_fwb'  , a_fwb )
               !
               a_fwb_ini = a_fwb_b
            ELSE                                                                                 !    as specified in namelist
               IF(lwp)   WRITE(numout,*) 'sbc_fwb : setting freshwater-budget from namelist rn_fwb0'
               a_fwb   = rn_fwb0
               a_fwb_b = 0._wp   ! used only the first year then it is replaced by a_fwb_ini
               !
               a_fwb_ini = glob_sum( 'sbcfwb', e1e2t(:,:) * ( ssh(:,:,Kmm) + snwice_mass(:,:) * r1_rho0 ) ) &
                  &      * rho0 / ( area * rday * REAL(nyear_len(1), wp) )
            END IF
            !
            IF(lwp)   WRITE(numout,*)
            IF(lwp)   WRITE(numout,*)'sbc_fwb : freshwater-budget at the end of previous year = ', a_fwb    , 'kg/m2/s'
            IF(lwp)   WRITE(numout,*)'          freshwater-budget at initial state            = ', a_fwb_ini, 'kg/m2/s'
            !
         ELSE
            ! at the end of year n:
            ikty = nyear_len(1) * 86400 / NINT(rn_Dt)
            IF( MOD( kt, ikty ) == 0 ) THEN   ! Update a_fwb at the last time step of a year
               !                                It should be the first time step of a year MOD(kt-1,ikty) but then the restart would be wrong
               !                                Hence, we make a small error here but the code is restartable
               a_fwb_b = a_fwb_ini
               ! mean sea level taking into account ice+snow
               a_fwb   = glob_sum( 'sbcfwb', e1e2t(:,:) * ( ssh(:,:,Kmm) + snwice_mass(:,:) * r1_rho0 ) )
               a_fwb   = a_fwb * rho0 / ( area * rday * REAL(nyear_len(1), wp) )   ! convert in kg/m2/s
            ENDIF
            !
         ENDIF
         !
         IF( MOD( kt-1, kn_fsbc ) == 0 ) THEN         ! correct the freshwater fluxes using previous year budget minus initial state
            zcoef = ( a_fwb - a_fwb_b )
            emp(:,:) = emp(:,:) + zcoef * tmask(:,:,1)
            qns(:,:) = qns(:,:) - zcoef * rcp * sst_m(:,:) * tmask(:,:,1) ! account for change to the heat budget due to fw correction
            ! outputs
            IF( iom_use('hflx_fwb_cea') )  CALL iom_put( 'hflx_fwb_cea', -zcoef * rcp * sst_m(:,:) * tmask(:,:,1) )
            IF( iom_use('vflx_fwb_cea') )  CALL iom_put( 'vflx_fwb_cea', -zcoef * tmask(:,:,1) )
         ENDIF
         ! Output restart information
         IF( lrst_oce ) THEN
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'sbc_fwb : writing FW-budget adjustment to ocean restart file at it = ', kt
            IF(lwp) WRITE(numout,*) '~~~~'
            CALL iom_rstput( kt, nitrst, numrow, 'a_fwb_b', a_fwb_b )
            CALL iom_rstput( kt, nitrst, numrow, 'a_fwb',   a_fwb   )
         END IF
         !
         IF( kt == nitend .AND. lwp ) THEN
            WRITE(numout,*) 'sbc_fwb : freshwater-budget at the end of simulation (year now) = ', a_fwb  , 'kg/m2/s'
            WRITE(numout,*) '          freshwater-budget at initial state                    = ', a_fwb_b, 'kg/m2/s'
         ENDIF
         !
      CASE ( 3 )                             !==  global fwf set to zero and spread out over erp area  ==!
         !
         ALLOCATE( ztmsk_neg(jpi,jpj) , ztmsk_pos(jpi,jpj) , ztmsk_tospread(jpi,jpj) , z_wgt(jpi,jpj) , zerp_cor(jpi,jpj) )
         !
         IF( MOD( kt-1, kn_fsbc ) == 0 ) THEN
            ztmsk_pos(:,:) = tmask_i(:,:)                      ! Select <0 and >0 area of erp
            WHERE( erp < 0._wp )   ztmsk_pos = 0._wp
            ztmsk_neg(:,:) = tmask_i(:,:) - ztmsk_pos(:,:)
            !                                                  ! fwf global mean (excluding ocean to ice/snow exchanges) 
            z_fwf     = glob_sum( 'sbcfwb', e1e2t(:,:) * ( emp(:,:) - rnf(:,:) - fwfisf_cav(:,:) - fwfisf_par(:,:) - snwice_fmass(:,:) ) ) / area
            !            
            IF( z_fwf < 0._wp ) THEN         ! spread out over >0 erp area to increase evaporation
               zsurf_pos = glob_sum( 'sbcfwb', e1e2t(:,:)*ztmsk_pos(:,:) )
               zsurf_tospread      = zsurf_pos
               ztmsk_tospread(:,:) = ztmsk_pos(:,:)
            ELSE                             ! spread out over <0 erp area to increase precipitation
               zsurf_neg = glob_sum( 'sbcfwb', e1e2t(:,:)*ztmsk_neg(:,:) )  ! Area filled by <0 and >0 erp 
               zsurf_tospread      = zsurf_neg
               ztmsk_tospread(:,:) = ztmsk_neg(:,:)
            ENDIF
            !
            zsum_fwf   = glob_sum( 'sbcfwb', e1e2t(:,:) * z_fwf )         ! fwf global mean over <0 or >0 erp area
!!gm :  zsum_fwf   = z_fwf * area   ???  it is right?  I think so....
            z_fwf_nsrf =  zsum_fwf / ( zsurf_tospread + rsmall )
            !                                                  ! weight to respect erp field 2D structure 
            zsum_erp   = glob_sum( 'sbcfwb', ztmsk_tospread(:,:) * erp(:,:) * e1e2t(:,:) )
            z_wgt(:,:) = ztmsk_tospread(:,:) * erp(:,:) / ( zsum_erp + rsmall )
            !                                                  ! final correction term to apply
            zerp_cor(:,:) = -1. * z_fwf_nsrf * zsurf_tospread * z_wgt(:,:)
            !
!!gm   ===>>>>  lbc_lnk should be useless as all the computation is done over the whole domain !
            CALL lbc_lnk( 'sbcfwb', zerp_cor, 'T', 1.0_wp )
            !
            emp(:,:) = emp(:,:) + zerp_cor(:,:)
            qns(:,:) = qns(:,:) - zerp_cor(:,:) * rcp * sst_m(:,:)  ! account for change to the heat budget due to fw correction
            erp(:,:) = erp(:,:) + zerp_cor(:,:)
            ! outputs
            IF( iom_use('hflx_fwb_cea') )  CALL iom_put( 'hflx_fwb_cea', -zerp_cor(:,:) * rcp * sst_m(:,:) )
            IF( iom_use('vflx_fwb_cea') )  CALL iom_put( 'vflx_fwb_cea', -zerp_cor(:,:) )
            !
            IF( lwp ) THEN                   ! control print
               IF( z_fwf < 0._wp ) THEN
                  WRITE(numout,*)'   z_fwf < 0'
                  WRITE(numout,*)'   SUM(erp+)     = ', SUM( ztmsk_tospread(:,:)*erp(:,:)*e1e2t(:,:) )*1.e-9,' Sv'
               ELSE
                  WRITE(numout,*)'   z_fwf >= 0'
                  WRITE(numout,*)'   SUM(erp-)     = ', SUM( ztmsk_tospread(:,:)*erp(:,:)*e1e2t(:,:) )*1.e-9,' Sv'
               ENDIF
               WRITE(numout,*)'   SUM(empG)     = ', SUM( z_fwf*e1e2t(:,:) )*1.e-9,' Sv'
               WRITE(numout,*)'   z_fwf         = ', z_fwf      ,' Kg/m2/s'
               WRITE(numout,*)'   z_fwf_nsrf    = ', z_fwf_nsrf ,' Kg/m2/s'
               WRITE(numout,*)'   MIN(zerp_cor) = ', MINVAL(zerp_cor) 
               WRITE(numout,*)'   MAX(zerp_cor) = ', MAXVAL(zerp_cor) 
            ENDIF
         ENDIF
         DEALLOCATE( ztmsk_neg , ztmsk_pos , ztmsk_tospread , z_wgt , zerp_cor )
         !
      CASE DEFAULT                           !==  you should never be there  ==!
         CALL ctl_stop( 'sbc_fwb : wrong nn_fwb value for the FreshWater Budget correction, choose either 1, 2 or 3' )
         !
      END SELECT
      !
   END SUBROUTINE sbc_fwb

   !!======================================================================
END MODULE sbcfwb
