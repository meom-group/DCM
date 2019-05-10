MODULE icestp
   !!======================================================================
   !!                       ***  MODULE  icestp  ***
   !! sea ice : Master routine for all the sea ice model
   !!=====================================================================
   !!
   !! The sea ice model SI3 (Sea Ice modelling Integrated Initiative),
   !!                        aka Sea Ice cube for its nickname
   !!
   !!    is originally based on LIM3, developed in Louvain-la-Neuve by: 
   !!       * Martin Vancoppenolle (UCL-ASTR, Belgium)
   !!       * Sylvain Bouillon (UCL-ASTR, Belgium)
   !!       * Miguel Angel Morales Maqueda (NOC-L, UK)
   !!      thanks to valuable earlier work by
   !!       * Thierry Fichefet
   !!       * Hugues Goosse
   !!      thanks also to the following persons who contributed
   !!       * Gurvan Madec, Claude Talandier, Christian Ethe (LOCEAN, France)
   !!       * Xavier Fettweis (UCL-ASTR), Ralph Timmermann (AWI, Germany)
   !!       * Bill Lipscomb (LANL), Cecilia Bitz (UWa) and Elisabeth Hunke (LANL), USA.
   !!
   !! SI3 has been made possible by a handful of persons who met as working group
   !!      (from France, Belgium, UK and Italy)
   !!    * Clement Rousset, Martin Vancoppenolle & Gurvan Madec (LOCEAN, France)
   !!    * Matthieu Chevalier & David Salas (Meteo France, France)
   !!    * Gilles Garric (Mercator Ocean, France)
   !!    * Thierry Fichefet & Francois Massonnet (UCL, Belgium)
   !!    * Ed Blockley & Jeff Ridley (Met Office, UK)
   !!    * Danny Feltham & David Schroeder (CPOM, UK)
   !!    * Yevgeny Aksenov (NOC, UK)
   !!    * Paul Holland (BAS, UK)
   !!    * Dorotea Iovino (CMCC, Italy)
   !!======================================================================
   !! History :  4.0  !  2018     (C. Rousset)      Original code SI3
   !!----------------------------------------------------------------------
#if defined key_si3
   !!----------------------------------------------------------------------
   !!   'key_si3'                                       SI3 sea-ice model
   !!----------------------------------------------------------------------
   !!   ice_stp       : sea-ice model time-stepping and update ocean SBC over ice-covered area
   !!   ice_init      : initialize sea-ice
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers
   USE dom_oce        ! ocean space and time domain
   USE c1d            ! 1D vertical configuration
   USE ice            ! sea-ice: variables
   USE ice1D          ! sea-ice: thermodynamical 1D variables
   !
   USE phycst         ! Define parameters for the routines
   USE eosbn2         ! equation of state
   USE sbc_oce        ! Surface boundary condition: ocean fields
   USE sbc_ice        ! Surface boundary condition: ice   fields
   !
   USE icesbc         ! sea-ice: Surface boundary conditions
   USE icedyn         ! sea-ice: dynamics
   USE icethd         ! sea-ice: thermodynamics
   USE icecor         ! sea-ice: corrections
   USE iceupdate      ! sea-ice: sea surface boundary condition update
   USE icedia         ! sea-ice: budget diagnostics
   USE icewri         ! sea-ice: outputs
   USE icerst         ! sea-ice: restarts
   USE icevar         ! sea-ice: operations
   USE icectl         ! sea-ice: control
   USE iceistate      ! sea-ice: initial state
   USE iceitd         ! sea-ice: remapping thickness distribution
   USE icealb         ! sea-ice: albedo
   !
   USE bdy_oce , ONLY : ln_bdy   ! flag for bdy
   USE bdyice         ! unstructured open boundary data for sea-ice
# if defined key_agrif
   USE agrif_ice
   USE agrif_ice_interp
# endif
   !
   USE in_out_manager ! I/O manager
   USE iom            ! I/O manager library
   USE lib_mpp        ! MPP library
   USE lib_fortran    ! fortran utilities (glob_sum + no signed zero)
   USE timing         ! Timing
   USE prtctl         ! Print control

   IMPLICIT NONE
   PRIVATE

   PUBLIC   ice_stp    ! called by sbcmod.F90
   PUBLIC   ice_init   ! called by sbcmod.F90

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/ICE 4.0 , NEMO Consortium (2018)
   !! $Id: icestp.F90 10535 2019-01-16 17:36:47Z clem $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE ice_stp( kt, ksbc )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE ice_stp  ***
      !!
      !! ** Purpose :   sea-ice model time-stepping and update ocean surface
      !!              boundary condition over ice-covered area
      !!
      !! ** Method  :   ice model time stepping
      !!              - call the ice dynamics routine
      !!              - call the ice advection/diffusion routine
      !!              - call the ice thermodynamics routine
      !!              - call the routine that computes mass and
      !!                heat fluxes at the ice/ocean interface
      !!              - save the outputs
      !!              - save the outputs for restart when necessary
      !!
      !! ** Action  : - time evolution of the LIM sea-ice model
      !!              - update all sbc variables below sea-ice:
      !!                utau, vtau, taum, wndm, qns , qsr, emp , sfx
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt      ! ocean time step
      INTEGER, INTENT(in) ::   ksbc    ! flux formulation (user defined, bulk, or Pure Coupled)
      !
      INTEGER ::   jl   ! dummy loop index
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('ice_stp')
      !
      !                                      !-----------------------!
      IF( MOD( kt-1, nn_fsbc ) == 0 ) THEN   ! --- Ice time step --- !
         !                                   !-----------------------!
         !
         kt_ice = kt                              ! -- Ice model time step
         !
         u_oce(:,:) = ssu_m(:,:)                  ! -- mean surface ocean current
         v_oce(:,:) = ssv_m(:,:)
         !
         CALL eos_fzp( sss_m(:,:) , t_bo(:,:) )   ! -- freezing temperature [Kelvin] (set to rt0 over land)
         t_bo(:,:) = ( t_bo(:,:) + rt0 ) * tmask(:,:,1) + rt0 * ( 1._wp - tmask(:,:,1) )
         !
         !                          !==  AGRIF Parent to Child  ==!
#if defined key_agrif
         !                              ! nbstep_ice ranges from 1 to the nb of child ocean steps inside one parent ice step
         IF( .NOT. Agrif_Root() )       nbstep_ice = MOD( nbstep_ice, Agrif_irhot() * Agrif_Parent(nn_fsbc) / nn_fsbc ) + 1
         !                              ! these calls must remain here for restartability purposes
                                        CALL agrif_interp_ice( 'T' ) 
                                        CALL agrif_interp_ice( 'U' )
                                        CALL agrif_interp_ice( 'V' )
#endif
                                        CALL store_fields             ! Store now ice values
         !
         !------------------------------------------------!
         ! --- Dynamical coupling with the atmosphere --- !
         !------------------------------------------------!
         ! It provides the following fields used in sea ice model:
         !    utau_ice, vtau_ice = surface ice stress [N/m2]
         !------------------------------------------------!
                                        CALL ice_sbc_tau( kt, ksbc, utau_ice, vtau_ice )          
         !-------------------------------------!
         ! --- ice dynamics and advection  --- !
         !-------------------------------------!
                                        CALL diag_set0                ! set diag of mass, heat and salt fluxes to 0
!                                       CALL ice_rst_opn( kt )        ! Open Ice restart file (if necessary) 
         !
         IF( ln_icedyn .AND. .NOT.lk_c1d )   &
            &                           CALL ice_dyn( kt )            ! -- Ice dynamics
         !
         !                          !==  lateral boundary conditions  ==!
         IF( ln_icethd .AND. ln_bdy )   CALL bdy_ice( kt )            ! -- bdy ice thermo
         !
         !                          !==  previous lead fraction and ice volume for flux calculations
                                        CALL ice_var_glo2eqv          ! h_i and h_s for ice albedo calculation
                                        CALL ice_var_agg(1)           ! at_i for coupling 
                                        CALL store_fields             ! Store now ice values
         !
         !------------------------------------------------------!
         ! --- Thermodynamical coupling with the atmosphere --- !
         !------------------------------------------------------!
         ! It provides the following fields used in sea ice model:
         !    emp_oce , emp_ice    = E-P over ocean and sea ice                    [Kg/m2/s]
         !    sprecip              = solid precipitation                           [Kg/m2/s]
         !    evap_ice             = sublimation                                   [Kg/m2/s]
         !    qsr_tot , qns_tot    = solar & non solar heat flux (total)           [W/m2]
         !    qsr_ice , qns_ice    = solar & non solar heat flux over ice          [W/m2]
         !    dqns_ice             = non solar  heat sensistivity                  [W/m2]
         !    qemp_oce, qemp_ice,  = sensible heat (associated with evap & precip) [W/m2]
         !    qprec_ice, qevap_ice
         !------------------------------------------------------!
                                        CALL ice_sbc_flx( kt, ksbc )
         !----------------------------!
         ! --- ice thermodynamics --- !
         !----------------------------!
         IF( ln_icethd )                CALL ice_thd( kt )            ! -- Ice thermodynamics      
         !
         IF( ln_icethd )                CALL ice_cor( kt , 2 )        ! -- Corrections
         !
                                        CALL ice_var_glo2eqv          ! necessary calls (at least for coupling)
                                        CALL ice_var_agg( 2 )         ! necessary calls (at least for coupling)
         !
                                        CALL ice_update_flx( kt )     ! -- Update ocean surface mass, heat and salt fluxes
         !
         IF( ln_icediahsb )             CALL ice_dia( kt )            ! -- Diagnostics outputs 
         !
                                        CALL ice_wri( kt )            ! -- Ice outputs 
         !
!        IF( lrst_ice )                 CALL ice_rst_write( kt )      ! -- Ice restart file 
         !
         IF( ln_icectl )                CALL ice_ctl( kt )            ! -- alerts in case of model crash
         !
      ENDIF   ! End sea-ice time step only

      !-------------------------!
      ! --- Ocean time step --- !
      !-------------------------!
      IF( ln_icedyn )                   CALL ice_update_tau( kt, ub(:,:,1), vb(:,:,1) )   ! -- update surface ocean stresses
!!gm   remark, the ocean-ice stress is not saved in ice diag call above .....  find a solution!!!
      !
      IF( ln_timing )   CALL timing_stop('ice_stp')
      !
   END SUBROUTINE ice_stp


   SUBROUTINE ice_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ice_init  ***
      !!
      !! ** purpose :   Initialize sea-ice parameters
      !!----------------------------------------------------------------------
      INTEGER :: ji, jj, ierr
      !!----------------------------------------------------------------------
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'Sea Ice Model: SI3 (Sea Ice modelling Integrated Initiative)' 
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~'
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'ice_init: Arrays allocation & Initialization of all routines & init state' 
      IF(lwp) WRITE(numout,*) '~~~~~~~~'
      !
      !                                ! Open the reference and configuration namelist files and namelist output file
      CALL ctl_opn( numnam_ice_ref, 'namelist_ice_ref',    'OLD',     'FORMATTED', 'SEQUENTIAL', -1, numout, lwp )
      CALL ctl_opn( numnam_ice_cfg, 'namelist_ice_cfg',    'OLD',     'FORMATTED', 'SEQUENTIAL', -1, numout, lwp )
      IF(lwm) CALL ctl_opn( numoni, 'output.namelist.ice', 'UNKNOWN', 'FORMATTED', 'SEQUENTIAL', -1, numout, lwp, 1 )
      !
      CALL par_init                ! set some ice run parameters
      !
      !                                ! Allocate the ice arrays (sbc_ice already allocated in sbc_init)
      ierr =        ice_alloc        ()      ! ice variables
      ierr = ierr + sbc_ice_alloc    ()      ! surface boundary conditions 
      ierr = ierr + ice1D_alloc      ()      ! thermodynamics
      !
      CALL mpp_sum( 'icestp', ierr )
      IF( ierr /= 0 )   CALL ctl_stop('STOP', 'ice_init : unable to allocate ice arrays')
      !
      CALL ice_itd_init                ! ice thickness distribution initialization
      !
      CALL ice_thd_init                ! set ice thermodynics parameters (clem: important to call it first for melt ponds)
      !
      !                                ! Initial sea-ice state
      IF( .NOT. ln_rstart ) THEN              ! start from rest: sea-ice deduced from sst
         CALL ice_istate_init
         CALL ice_istate
      ELSE                                    ! start from a restart file
         CALL ice_rst_read
      ENDIF
      CALL ice_var_glo2eqv
      CALL ice_var_agg(1)
      !
      CALL ice_sbc_init                ! set ice-ocean and ice-atm. coupling parameters
      !
      CALL ice_dyn_init                ! set ice dynamics parameters
      !
      CALL ice_update_init             ! ice surface boundary condition
      !
      CALL ice_alb_init                ! ice surface albedo
      !
      CALL ice_dia_init                ! initialization for diags
      !
      fr_i  (:,:)   = at_i(:,:)        ! initialisation of sea-ice fraction
      tn_ice(:,:,:) = t_su(:,:,:)      ! initialisation of surface temp for coupled simu
      !
      !                                ! set max concentration in both hemispheres
      WHERE( gphit(:,:) > 0._wp )   ;   rn_amax_2d(:,:) = rn_amax_n  ! NH
      ELSEWHERE                     ;   rn_amax_2d(:,:) = rn_amax_s  ! SH
      END WHERE

      IF( ln_rstart )   CALL iom_close( numrir )  ! close input ice restart file
      !
   END SUBROUTINE ice_init


   SUBROUTINE par_init
      !!-------------------------------------------------------------------
      !!                  ***  ROUTINE par_init ***
      !!
      !! ** Purpose :   Definition generic parameters for ice model
      !!
      !! ** Method  :   Read namelist and check the parameter
      !!                values called at the first timestep (nit000)
      !!
      !! ** input   :   Namelist nampar
      !!-------------------------------------------------------------------
      INTEGER  ::   ios                 ! Local integer
      !!
#if defined key_drakkar
      CHARACTER(LEN=20) :: cl_no     ! provide space for segment number
#endif
      NAMELIST/nampar/ jpl, nlay_i, nlay_s, ln_virtual_itd, ln_icedyn, ln_icethd, rn_amax_n, rn_amax_s,  &
         &             cn_icerst_in, cn_icerst_indir, cn_icerst_out, cn_icerst_outdir
      !!-------------------------------------------------------------------
      !
      REWIND( numnam_ice_ref )      ! Namelist nampar in reference namelist : Parameters for ice
      READ  ( numnam_ice_ref, nampar, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'nampar in reference namelist', lwp )
      REWIND( numnam_ice_cfg )      ! Namelist nampar in configuration namelist : Parameters for ice
      READ  ( numnam_ice_cfg, nampar, IOSTAT = ios, ERR = 902 )
902   IF( ios > 0 )   CALL ctl_nam ( ios , 'nampar in configuration namelist', lwp )
      IF(lwm) WRITE( numoni, nampar )

#if defined key_drakkar
      ! Add extension (job number to the restart dir. Differ for restart input and restart output
!{ DRAKKAR modification : NEMO reads restart files :<CN_ICERST_INDIR>.<<nn_no-1>>/<CN_ICERST_IN>-<<nn_no -1 >>_<RANK>.nc
      WRITE(cl_no,*) nn_no-1 ; cl_no = TRIM(ADJUSTL(cl_no) )
      cn_icerst_indir=TRIM(cn_icerst_indir)//'.'//TRIM(cl_no)
      cn_icerst_in= TRIM(cn_icerst_in)//'-'//TRIM(cl_no)

!  DRAKKAR modification : NEMO write restart files :<CN_ICERST_OUTDIR>.<<nn_no>>/<CN_ICERST_OUT>-<<nn_no >>_<RANK>.nc
      WRITE(cl_no,*) nn_no   ; cl_no = TRIM(ADJUSTL(cl_no) )
      cn_icerst_outdir=TRIM(cn_icerst_outdir)//'.'//TRIM(cl_no)
      cn_icerst_out= TRIM(cn_icerst_out)//'-'//TRIM(cl_no)
#endif
      !
      IF(lwp) THEN                  ! control print
         WRITE(numout,*)
         WRITE(numout,*) '   par_init: ice parameters shared among all the routines'
         WRITE(numout,*) '   ~~~~~~~~'
         WRITE(numout,*) '      Namelist nampar: '
         WRITE(numout,*) '         number of ice  categories                           jpl       = ', jpl
         WRITE(numout,*) '         number of ice  layers                               nlay_i    = ', nlay_i
         WRITE(numout,*) '         number of snow layers                               nlay_s    = ', nlay_s
         WRITE(numout,*) '         virtual ITD param for jpl=1 (T) or not (F)     ln_virtual_itd = ', ln_virtual_itd
         WRITE(numout,*) '         Ice dynamics       (T) or not (F)                   ln_icedyn = ', ln_icedyn
         WRITE(numout,*) '         Ice thermodynamics (T) or not (F)                   ln_icethd = ', ln_icethd
         WRITE(numout,*) '         maximum ice concentration for NH                              = ', rn_amax_n 
         WRITE(numout,*) '         maximum ice concentration for SH                              = ', rn_amax_s
      ENDIF
      !                                        !--- check consistency
      IF ( jpl > 1 .AND. ln_virtual_itd ) THEN
         ln_virtual_itd = .FALSE.
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) '   ln_virtual_itd forced to false as jpl>1, no need with multiple categories to emulate them'
      ENDIF
      !
      IF( ln_cpl .AND. nn_cats_cpl /= 1 .AND. nn_cats_cpl /= jpl ) THEN
         CALL ctl_stop( 'STOP', 'par_init: in coupled mode, nn_cats_cpl should be either 1 or jpl' )
      ENDIF
      !
      IF( ln_bdy .AND. ln_icediachk )   CALL ctl_warn('par_init: online conservation check does not work with BDY')
      !
      rdt_ice   = REAL(nn_fsbc) * rdt          !--- sea-ice timestep and its inverse
      r1_rdtice = 1._wp / rdt_ice
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) '      ice timestep rdt_ice = nn_fsbc*rdt = ', rdt_ice
      !
      r1_nlay_i = 1._wp / REAL( nlay_i, wp )   !--- inverse of nlay_i and nlay_s
      r1_nlay_s = 1._wp / REAL( nlay_s, wp )
      !
   END SUBROUTINE par_init


   SUBROUTINE store_fields
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE store_fields  ***
      !!
      !! ** purpose :  store ice variables at "before" time step
      !!----------------------------------------------------------------------
      INTEGER  ::   ji, jj, jl      ! dummy loop index
      !!----------------------------------------------------------------------
      !
      a_i_b (:,:,:)   = a_i (:,:,:)     ! ice area
      v_i_b (:,:,:)   = v_i (:,:,:)     ! ice volume
      v_s_b (:,:,:)   = v_s (:,:,:)     ! snow volume
      sv_i_b(:,:,:)   = sv_i(:,:,:)     ! salt content
      oa_i_b(:,:,:)   = oa_i(:,:,:)     ! areal age content
      e_s_b (:,:,:,:) = e_s (:,:,:,:)   ! snow thermal energy
      e_i_b (:,:,:,:) = e_i (:,:,:,:)   ! ice thermal energy
      WHERE( a_i_b(:,:,:) >= epsi20 )
         h_i_b(:,:,:) = v_i_b(:,:,:) / a_i_b(:,:,:)   ! ice thickness
         h_s_b(:,:,:) = v_s_b(:,:,:) / a_i_b(:,:,:)   ! snw thickness
      ELSEWHERE
         h_i_b(:,:,:) = 0._wp
         h_s_b(:,:,:) = 0._wp
      END WHERE
      
      WHERE( a_ip(:,:,:) >= epsi20 )
         h_ip_b(:,:,:) = v_ip(:,:,:) / a_ip(:,:,:)   ! ice pond thickness
      ELSEWHERE
         h_ip_b(:,:,:) = 0._wp
      END WHERE
      !
      ! ice velocities & total concentration
      at_i_b(:,:)  = SUM( a_i_b(:,:,:), dim=3 )
      u_ice_b(:,:) = u_ice(:,:)
      v_ice_b(:,:) = v_ice(:,:)
      !
   END SUBROUTINE store_fields


   SUBROUTINE diag_set0
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE diag_set0  ***
      !!
      !! ** purpose :  set ice-ocean and ice-atm. fluxes to zeros at the beggining
      !!               of the time step
      !!----------------------------------------------------------------------
      INTEGER  ::   ji, jj      ! dummy loop index
      !!----------------------------------------------------------------------
      sfx    (:,:) = 0._wp   ;
      sfx_bri(:,:) = 0._wp   ;   sfx_lam(:,:) = 0._wp
      sfx_sni(:,:) = 0._wp   ;   sfx_opw(:,:) = 0._wp
      sfx_bog(:,:) = 0._wp   ;   sfx_dyn(:,:) = 0._wp
      sfx_bom(:,:) = 0._wp   ;   sfx_sum(:,:) = 0._wp
      sfx_res(:,:) = 0._wp   ;   sfx_sub(:,:) = 0._wp
      !
      wfx_snw(:,:) = 0._wp   ;   wfx_ice(:,:) = 0._wp
      wfx_sni(:,:) = 0._wp   ;   wfx_opw(:,:) = 0._wp
      wfx_bog(:,:) = 0._wp   ;   wfx_dyn(:,:) = 0._wp
      wfx_bom(:,:) = 0._wp   ;   wfx_sum(:,:) = 0._wp
      wfx_res(:,:) = 0._wp   ;   wfx_sub(:,:) = 0._wp
      wfx_spr(:,:) = 0._wp   ;   wfx_lam(:,:) = 0._wp  
      wfx_snw_dyn(:,:) = 0._wp ; wfx_snw_sum(:,:) = 0._wp
      wfx_snw_sub(:,:) = 0._wp ; wfx_ice_sub(:,:) = 0._wp
      wfx_snw_sni(:,:) = 0._wp 
      wfx_pnd(:,:) = 0._wp

      hfx_thd(:,:) = 0._wp   ;
      hfx_snw(:,:) = 0._wp   ;   hfx_opw(:,:) = 0._wp
      hfx_bog(:,:) = 0._wp   ;   hfx_dyn(:,:) = 0._wp
      hfx_bom(:,:) = 0._wp   ;   hfx_sum(:,:) = 0._wp
      hfx_res(:,:) = 0._wp   ;   hfx_sub(:,:) = 0._wp
      hfx_spr(:,:) = 0._wp   ;   hfx_dif(:,:) = 0._wp
      hfx_err_rem(:,:) = 0._wp
      hfx_err_dif(:,:) = 0._wp
      wfx_err_sub(:,:) = 0._wp
      !
      afx_tot(:,:) = 0._wp   ;
      !
      diag_heat(:,:) = 0._wp ;   diag_sice(:,:) = 0._wp
      diag_vice(:,:) = 0._wp ;   diag_vsnw(:,:) = 0._wp

      ! SIMIP diagnostics
      qcn_ice_bot(:,:,:) = 0._wp ; qcn_ice_top(:,:,:) = 0._wp ! conductive fluxes
      t_si       (:,:,:) = rt0   ! temp at the ice-snow interface

      tau_icebfr(:,:)   = 0._wp   ! landfast ice param only (clem: important to keep the init here)
      cnd_ice   (:,:,:) = 0._wp   ! initialisation: effective conductivity at the top of ice/snow (ln_cndflx=T)
      qtr_ice_bot(:,:,:) = 0._wp  ! initialization: part of solar radiation transmitted through the ice needed at least for outputs
      !
      ! for control checks (ln_icediachk)
      diag_trp_vi(:,:) = 0._wp   ;   diag_trp_vs(:,:) = 0._wp
      diag_trp_ei(:,:) = 0._wp   ;   diag_trp_es(:,:) = 0._wp
      diag_trp_sv(:,:) = 0._wp

   END SUBROUTINE diag_set0

#else
   !!----------------------------------------------------------------------
   !!   Default option           Dummy module         NO SI3 sea-ice model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE ice_stp ( kt, ksbc )     ! Dummy routine
      INTEGER, INTENT(in) ::   kt, ksbc
      WRITE(*,*) 'ice_stp: You should not have seen this print! error?', kt
   END SUBROUTINE ice_stp
   SUBROUTINE ice_init                 ! Dummy routine
   END SUBROUTINE ice_init
#endif

   !!======================================================================
END MODULE icestp
