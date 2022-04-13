MODULE trabbl
   !!==============================================================================
   !!                       ***  MODULE  trabbl  ***
   !! Ocean physics :  advective and/or diffusive bottom boundary layer scheme
   !!==============================================================================
   !! History :  OPA  ! 1996-06  (L. Mortier)  Original code
   !!            8.0  ! 1997-11  (G. Madec)    Optimization
   !!   NEMO     1.0  ! 2002-08  (G. Madec)  free form + modules
   !!             -   ! 2004-01  (A. de Miranda, G. Madec, J.M. Molines ) add advective bbl
   !!            3.3  ! 2009-11  (G. Madec)  merge trabbl and trabbl_adv + style + optimization
   !!             -   ! 2010-04  (G. Madec)  Campin & Goosse advective bbl
   !!             -   ! 2010-06  (C. Ethe, G. Madec)  merge TRA-TRC
   !!             -   ! 2010-11  (G. Madec) add mbk. arrays associated to the deepest ocean level
   !!             -   ! 2013-04  (F. Roquet, G. Madec)  use of eosbn2 instead of local hard coded alpha and beta
   !!            4.0  ! 2017-04  (G. Madec)  ln_trabbl namelist variable instead of a CPP key
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   tra_bbl_alloc : allocate trabbl arrays
   !!   tra_bbl       : update the tracer trends due to the bottom boundary layer (advective and/or diffusive)
   !!   tra_bbl_dif   : generic routine to compute bbl diffusive trend
   !!   tra_bbl_adv   : generic routine to compute bbl advective trend
   !!   bbl           : computation of bbl diffu. flux coef. & transport in bottom boundary layer
   !!   tra_bbl_init  : initialization, namelist read, parameters control
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and active tracers
   USE dom_oce        ! ocean space and time domain
   USE phycst         ! physical constant
   USE eosbn2         ! equation of state
   USE trd_oce        ! trends: ocean variables
   USE trdtra         ! trends: active tracers
   !
   USE iom            ! IOM library
   USE in_out_manager ! I/O manager
   USE lbclnk         ! ocean lateral boundary conditions
   USE prtctl         ! Print control
   USE timing         ! Timing
   USE lib_fortran    ! Fortran utilities (allows no signed zero when 'key_nosignedzero' defined)

   IMPLICIT NONE
   PRIVATE

   PUBLIC   tra_bbl       !  routine called by step.F90
   PUBLIC   tra_bbl_init  !  routine called by nemogcm.F90
   PUBLIC   tra_bbl_dif   !  routine called by trcbbl.F90
   PUBLIC   tra_bbl_adv   !     -      -          -
   PUBLIC   bbl           !  routine called by trcbbl.F90 and dtadyn.F90

   !                                !!* Namelist nambbl *
   LOGICAL , PUBLIC ::   ln_trabbl   !: bottom boundary layer flag
   INTEGER , PUBLIC ::   nn_bbl_ldf  !: =1   : diffusive bbl or not (=0)
   INTEGER , PUBLIC ::   nn_bbl_adv  !: =1/2 : advective bbl or not (=0)
   !                                            !  =1 : advective bbl using the bottom ocean velocity
   !                                            !  =2 :     -      -  using utr_bbl proportional to grad(rho)
   REAL(wp), PUBLIC ::   rn_ahtbbl   !: along slope bbl diffusive coefficient [m2/s]
   REAL(wp), PUBLIC ::   rn_gambbl   !: lateral coeff. for bottom boundary layer scheme [s]
#if defined key_drakkar
   !                                !!* Namelist nambbl_drk *
   LOGICAL , PUBLIC ::   ln_kriteria !: =T use k- criteria instead of depth criteria (=F)
#endif

   LOGICAL , PUBLIC ::   l_bbl       !: flag to compute bbl diffu. flux coef and transport

   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:), PUBLIC ::   utr_bbl  , vtr_bbl   ! u- (v-) transport in the bottom boundary layer
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:), PUBLIC ::   ahu_bbl  , ahv_bbl   ! masked diffusive bbl coeff. at u & v-pts

   INTEGER , ALLOCATABLE, SAVE, DIMENSION(:,:), PUBLIC ::   mbku_d   , mbkv_d      ! vertical index of the "lower" bottom ocean U/V-level (PUBLIC for TAM)
   INTEGER , ALLOCATABLE, SAVE, DIMENSION(:,:), PUBLIC ::   mgrhu    , mgrhv       ! = +/-1, sign of grad(H) in u-(v-)direction (PUBLIC for TAM)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)         ::   ahu_bbl_0, ahv_bbl_0   ! diffusive bbl flux coefficients at u and v-points
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:), PUBLIC ::   e3u_bbl_0, e3v_bbl_0   ! thichness of the bbl (e3) at u and v-points (PUBLIC for TAM)

   !! * Substitutions
#  include "do_loop_substitute.h90"
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: trabbl.F90 15053 2021-06-24 15:39:38Z clem $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION tra_bbl_alloc()
      !!----------------------------------------------------------------------
      !!                ***  FUNCTION tra_bbl_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( utr_bbl  (jpi,jpj) , ahu_bbl  (jpi,jpj) , mbku_d(jpi,jpj) , mgrhu(jpi,jpj) ,     &
         &      vtr_bbl  (jpi,jpj) , ahv_bbl  (jpi,jpj) , mbkv_d(jpi,jpj) , mgrhv(jpi,jpj) ,     &
         &      ahu_bbl_0(jpi,jpj) , ahv_bbl_0(jpi,jpj) ,                                        &
         &      e3u_bbl_0(jpi,jpj) , e3v_bbl_0(jpi,jpj) ,                                    STAT=tra_bbl_alloc )
         !
      CALL mpp_sum ( 'trabbl', tra_bbl_alloc )
      IF( tra_bbl_alloc > 0 )   CALL ctl_warn('tra_bbl_alloc: allocation of arrays failed.')
   END FUNCTION tra_bbl_alloc


   SUBROUTINE tra_bbl( kt, Kbb, Kmm, pts, Krhs )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE bbl  ***
      !!
      !! ** Purpose :   Compute the before tracer (t & s) trend associated
      !!              with the bottom boundary layer and add it to the general
      !!              trend of tracer equations.
      !!
      !! ** Method  :   Depending on namtra_bbl namelist parameters the bbl
      !!              diffusive and/or advective contribution to the tracer trend
      !!              is added to the general tracer trend
      !!----------------------------------------------------------------------
      INTEGER,                                   INTENT(in   ) :: kt              ! ocean time-step
      INTEGER,                                   INTENT(in   ) :: Kbb, Kmm, Krhs  ! time level indices
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpts,jpt), INTENT(inout) :: pts             ! active tracers and RHS of tracer equation
      !
      INTEGER  ::   ji, jj, jk   ! Dummy loop indices
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::   ztrdt, ztrds
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start( 'tra_bbl')
      !
      IF( l_trdtra )   THEN                         !* Save the T-S input trends
         ALLOCATE( ztrdt(jpi,jpj,jpk), ztrds(jpi,jpj,jpk) )
         ztrdt(:,:,:) = pts(:,:,:,jp_tem,Krhs)
         ztrds(:,:,:) = pts(:,:,:,jp_sal,Krhs)
      ENDIF

      IF( l_bbl )   CALL bbl( kt, nit000, 'TRA', Kbb, Kmm )   !* bbl coef. and transport (only if not already done in trcbbl)

      IF( nn_bbl_ldf == 1 ) THEN                    !* Diffusive bbl
         !
         CALL tra_bbl_dif( pts(:,:,:,:,Kbb), pts(:,:,:,:,Krhs), jpts, Kmm )
         IF( sn_cfctl%l_prtctl )  &
         CALL prt_ctl( tab3d_1=pts(:,:,:,jp_tem,Krhs), clinfo1=' bbl_ldf  - Ta: ', mask1=tmask, &
            &          tab3d_2=pts(:,:,:,jp_sal,Krhs), clinfo2=           ' Sa: ', mask2=tmask, clinfo3='tra' )
         CALL iom_put( "ahu_bbl", ahu_bbl )   ! bbl diffusive flux i-coef
         CALL iom_put( "ahv_bbl", ahv_bbl )   ! bbl diffusive flux j-coef
         !
      ENDIF
      !
      IF( nn_bbl_adv /= 0 ) THEN                    !* Advective bbl
         !
         CALL tra_bbl_adv( pts(:,:,:,:,Kbb), pts(:,:,:,:,Krhs), jpts, Kmm )
         IF(sn_cfctl%l_prtctl)   &
         CALL prt_ctl( tab3d_1=pts(:,:,:,jp_tem,Krhs), clinfo1=' bbl_adv  - Ta: ', mask1=tmask, &
            &          tab3d_2=pts(:,:,:,jp_sal,Krhs), clinfo2=           ' Sa: ', mask2=tmask, clinfo3='tra' )
         CALL iom_put( "uoce_bbl", utr_bbl )  ! bbl i-transport
         CALL iom_put( "voce_bbl", vtr_bbl )  ! bbl j-transport
         !
      ENDIF

      IF( l_trdtra )   THEN                      ! send the trends for further diagnostics
         ztrdt(:,:,:) = pts(:,:,:,jp_tem,Krhs) - ztrdt(:,:,:)
         ztrds(:,:,:) = pts(:,:,:,jp_sal,Krhs) - ztrds(:,:,:)
         CALL trd_tra( kt, Kmm, Krhs, 'TRA', jp_tem, jptra_bbl, ztrdt )
         CALL trd_tra( kt, Kmm, Krhs, 'TRA', jp_sal, jptra_bbl, ztrds )
         DEALLOCATE( ztrdt, ztrds )
      ENDIF
      !
      IF( ln_timing )  CALL timing_stop( 'tra_bbl')
      !
   END SUBROUTINE tra_bbl


   SUBROUTINE tra_bbl_dif( pt, pt_rhs, kjpt, Kmm )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_bbl_dif  ***
      !!
      !! ** Purpose :   Computes the bottom boundary horizontal and vertical
      !!                advection terms.
      !!
      !! ** Method  : * diffusive bbl only (nn_bbl_ldf=1) :
      !!        When the product grad( rho) * grad(h) < 0 (where grad is an
      !!      along bottom slope gradient) an additional lateral 2nd order
      !!      diffusion along the bottom slope is added to the general
      !!      tracer trend, otherwise the additional trend is set to 0.
      !!      A typical value of ahbt is 2000 m2/s (equivalent to
      !!      a downslope velocity of 20 cm/s if the condition for slope
      !!      convection is satified)
      !!
      !! ** Action  :   pt_rhs   increased by the bbl diffusive trend
      !!
      !! References : Beckmann, A., and R. Doscher, 1997, J. Phys.Oceanogr., 581-591.
      !!              Campin, J.-M., and H. Goosse, 1999, Tellus, 412-430.
      !!----------------------------------------------------------------------
      INTEGER                              , INTENT(in   ) ::   kjpt   ! number of tracers
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(in   ) ::   pt     ! before and now tracer fields
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(inout) ::   pt_rhs ! tracer trend
      INTEGER                              , INTENT(in   ) ::   Kmm    ! time level indices
      !
      INTEGER  ::   ji, jj, jn   ! dummy loop indices
      INTEGER  ::   ik           ! local integers
      REAL(wp) ::   zbtr         ! local scalars
      REAL(wp), DIMENSION(A2D(nn_hls)) ::   zptb   ! workspace
      !!----------------------------------------------------------------------
      !
      DO jn = 1, kjpt                                     ! tracer loop
         !                                                ! ===========
         DO_2D( 1, 1, 1, 1 )
            ik = mbkt(ji,jj)                             ! bottom T-level index
            zptb(ji,jj) = pt(ji,jj,ik,jn)                ! bottom before T and S
         END_2D
         !
         DO_2D( 0, 0, 0, 0 )                               ! Compute the trend
            ik = mbkt(ji,jj)                            ! bottom T-level index
            pt_rhs(ji,jj,ik,jn) = pt_rhs(ji,jj,ik,jn)                                                  &
               &                + (  ahu_bbl(ji  ,jj  ) * ( zptb(ji+1,jj  ) - zptb(ji  ,jj  ) )     &
               &                   - ahu_bbl(ji-1,jj  ) * ( zptb(ji  ,jj  ) - zptb(ji-1,jj  ) )     &
               &                   + ahv_bbl(ji  ,jj  ) * ( zptb(ji  ,jj+1) - zptb(ji  ,jj  ) )     &
               &                   - ahv_bbl(ji  ,jj-1) * ( zptb(ji  ,jj  ) - zptb(ji  ,jj-1) )  )  &
               &                * r1_e1e2t(ji,jj) / e3t(ji,jj,ik,Kmm)
         END_2D
         !                                                  ! ===========
      END DO                                                ! end tracer
      !                                                     ! ===========
   END SUBROUTINE tra_bbl_dif


   ! NOTE: [tiling] tiling changes the results, but only the order of floating point operations is different
   SUBROUTINE tra_bbl_adv( pt, pt_rhs, kjpt, Kmm )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trc_bbl  ***
      !!
      !! ** Purpose :   Compute the before passive tracer trend associated
      !!     with the bottom boundary layer and add it to the general trend
      !!     of tracer equations.
      !! ** Method  :   advective bbl (nn_bbl_adv = 1 or 2) :
      !!      nn_bbl_adv = 1   use of the ocean near bottom velocity as bbl velocity
      !!      nn_bbl_adv = 2   follow Campin and Goosse (1999) implentation i.e.
      !!                       transport proportional to the along-slope density gradient
      !!
      !! References : Beckmann, A., and R. Doscher, 1997, J. Phys.Oceanogr., 581-591.
      !!              Campin, J.-M., and H. Goosse, 1999, Tellus, 412-430.
      !!----------------------------------------------------------------------
      INTEGER                              , INTENT(in   ) ::   kjpt   ! number of tracers
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(in   ) ::   pt     ! before and now tracer fields
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(inout) ::   pt_rhs ! tracer trend
      INTEGER                              , INTENT(in   ) ::   Kmm    ! time level indices
      !
      INTEGER  ::   ji, jj, jk, jn           ! dummy loop indices
      INTEGER  ::   iis , iid , ijs , ijd    ! local integers
      INTEGER  ::   ikus, ikud, ikvs, ikvd   !   -       -
      REAL(wp) ::   zbtr, ztra               ! local scalars
      REAL(wp) ::   zu_bbl, zv_bbl           !   -      -
      !!----------------------------------------------------------------------
      !                                                          ! ===========
      DO jn = 1, kjpt                                            ! tracer loop
         !                                                       ! ===========
         DO_2D_OVR( 1, 0, 1, 0 )            ! CAUTION start from i=1 to update i=2 when cyclic east-west
            IF( utr_bbl(ji,jj) /= 0.e0 ) THEN            ! non-zero i-direction bbl advection
               ! down-slope i/k-indices (deep)      &   up-slope i/k indices (shelf)
               iid  = ji + MAX( 0, mgrhu(ji,jj) )   ;   iis  = ji + 1 - MAX( 0, mgrhu(ji,jj) )
               ikud = mbku_d(ji,jj)                 ;   ikus = mbku(ji,jj)
               zu_bbl = ABS( utr_bbl(ji,jj) )
               !
               !                                               ! up  -slope T-point (shelf bottom point)
               zbtr = r1_e1e2t(iis,jj) / e3t(iis,jj,ikus,Kmm)
               ztra = zu_bbl * ( pt(iid,jj,ikus,jn) - pt(iis,jj,ikus,jn) ) * zbtr
               pt_rhs(iis,jj,ikus,jn) = pt_rhs(iis,jj,ikus,jn) + ztra
               !
               DO jk = ikus, ikud-1                            ! down-slope upper to down T-point (deep column)
                  zbtr = r1_e1e2t(iid,jj) / e3t(iid,jj,jk,Kmm)
                  ztra = zu_bbl * ( pt(iid,jj,jk+1,jn) - pt(iid,jj,jk,jn) ) * zbtr
                  pt_rhs(iid,jj,jk,jn) = pt_rhs(iid,jj,jk,jn) + ztra
               END DO
               !
               zbtr = r1_e1e2t(iid,jj) / e3t(iid,jj,ikud,Kmm)
               ztra = zu_bbl * ( pt(iis,jj,ikus,jn) - pt(iid,jj,ikud,jn) ) * zbtr
               pt_rhs(iid,jj,ikud,jn) = pt_rhs(iid,jj,ikud,jn) + ztra
            ENDIF
            !
            IF( vtr_bbl(ji,jj) /= 0.e0 ) THEN            ! non-zero j-direction bbl advection
               ! down-slope j/k-indices (deep)        &   up-slope j/k indices (shelf)
               ijd  = jj + MAX( 0, mgrhv(ji,jj) )     ;   ijs  = jj + 1 - MAX( 0, mgrhv(ji,jj) )
               ikvd = mbkv_d(ji,jj)                   ;   ikvs = mbkv(ji,jj)
               zv_bbl = ABS( vtr_bbl(ji,jj) )
               !
               ! up  -slope T-point (shelf bottom point)
               zbtr = r1_e1e2t(ji,ijs) / e3t(ji,ijs,ikvs,Kmm)
               ztra = zv_bbl * ( pt(ji,ijd,ikvs,jn) - pt(ji,ijs,ikvs,jn) ) * zbtr
               pt_rhs(ji,ijs,ikvs,jn) = pt_rhs(ji,ijs,ikvs,jn) + ztra
               !
               DO jk = ikvs, ikvd-1                            ! down-slope upper to down T-point (deep column)
                  zbtr = r1_e1e2t(ji,ijd) / e3t(ji,ijd,jk,Kmm)
                  ztra = zv_bbl * ( pt(ji,ijd,jk+1,jn) - pt(ji,ijd,jk,jn) ) * zbtr
                  pt_rhs(ji,ijd,jk,jn) = pt_rhs(ji,ijd,jk,jn)  + ztra
               END DO
               !                                               ! down-slope T-point (deep bottom point)
               zbtr = r1_e1e2t(ji,ijd) / e3t(ji,ijd,ikvd,Kmm)
               ztra = zv_bbl * ( pt(ji,ijs,ikvs,jn) - pt(ji,ijd,ikvd,jn) ) * zbtr
               pt_rhs(ji,ijd,ikvd,jn) = pt_rhs(ji,ijd,ikvd,jn) + ztra
            ENDIF
         END_2D
         !                                                       ! ===========
      END DO                                                     ! end tracer
      !                                                          ! ===========
   END SUBROUTINE tra_bbl_adv


   SUBROUTINE bbl( kt, kit000, cdtype, Kbb, Kmm )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE bbl  ***
      !!
      !! ** Purpose :   Computes the bottom boundary horizontal and vertical
      !!                advection terms.
      !!
      !! ** Method  : * diffusive bbl (nn_bbl_ldf=1) :
      !!        When the product grad( rho) * grad(h) < 0 (where grad is an
      !!      along bottom slope gradient) an additional lateral 2nd order
      !!      diffusion along the bottom slope is added to the general
      !!      tracer trend, otherwise the additional trend is set to 0.
      !!      A typical value of ahbt is 2000 m2/s (equivalent to
      !!      a downslope velocity of 20 cm/s if the condition for slope
      !!      convection is satified)
      !!              * advective bbl (nn_bbl_adv=1 or 2) :
      !!      nn_bbl_adv = 1   use of the ocean velocity as bbl velocity
      !!      nn_bbl_adv = 2   follow Campin and Goosse (1999) implentation
      !!        i.e. transport proportional to the along-slope density gradient
      !!
      !!      NB: the along slope density gradient is evaluated using the
      !!      local density (i.e. referenced at a common local depth).
      !!
      !! References : Beckmann, A., and R. Doscher, 1997, J. Phys.Oceanogr., 581-591.
      !!              Campin, J.-M., and H. Goosse, 1999, Tellus, 412-430.
      !!----------------------------------------------------------------------
      INTEGER         , INTENT(in   ) ::   kt       ! ocean time-step index
      INTEGER         , INTENT(in   ) ::   kit000   ! first time step index
      CHARACTER(len=3), INTENT(in   ) ::   cdtype   ! =TRA or TRC (tracer indicator)
      INTEGER         , INTENT(in   ) ::   Kbb, Kmm ! ocean time level index
      !
      INTEGER  ::   ji, jj                    ! dummy loop indices
      INTEGER  ::   ik                        ! local integers
      INTEGER  ::   iis, iid, ikus, ikud      !   -       -
      INTEGER  ::   ijs, ijd, ikvs, ikvd      !   -       -
      REAL(wp) ::   za, zb, zgdrho            ! local scalars
      REAL(wp) ::   zsign, zsigna, zgbbl      !   -      -
      REAL(wp), DIMENSION(A2D(nn_hls),jpts)   :: zts, zab         ! 3D workspace
      REAL(wp), DIMENSION(A2D(nn_hls))        :: zub, zvb, zdep   ! 2D workspace
      !!----------------------------------------------------------------------
      !
      IF( .NOT. l_istiled .OR. ntile == 1 )  THEN                       ! Do only on the first tile
         IF( kt == kit000 )  THEN
            IF(lwp)  WRITE(numout,*)
            IF(lwp)  WRITE(numout,*) 'trabbl:bbl : Compute bbl velocities and diffusive coefficients in ', cdtype
            IF(lwp)  WRITE(numout,*) '~~~~~~~~~~'
         ENDIF
      ENDIF
      !                                        !* bottom variables (T, S, alpha, beta, depth, velocity)
      DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
         ik = mbkt(ji,jj)                             ! bottom T-level index
         zts (ji,jj,jp_tem) = ts(ji,jj,ik,jp_tem,Kbb) ! bottom before T and S
         zts (ji,jj,jp_sal) = ts(ji,jj,ik,jp_sal,Kbb)
         !
         zdep(ji,jj) = gdept(ji,jj,ik,Kmm)            ! bottom T-level reference depth
         zub (ji,jj) = uu(ji,jj,mbku(ji,jj),Kmm)      ! bottom velocity
         zvb (ji,jj) = vv(ji,jj,mbkv(ji,jj),Kmm)
      END_2D
      !
      CALL eos_rab( zts, zdep, zab, Kmm )
      !
      !                                   !-------------------!
      IF( nn_bbl_ldf == 1 ) THEN          !   diffusive bbl   !
         !                                !-------------------!
         DO_2D_OVR( 1, 0, 1, 0 )                   ! (criteria for non zero flux: grad(rho).grad(h) < 0 )
            !                                                   ! i-direction
            za = zab(ji+1,jj,jp_tem) + zab(ji,jj,jp_tem)              ! 2*(alpha,beta) at u-point
            zb = zab(ji+1,jj,jp_sal) + zab(ji,jj,jp_sal)
            !                                                         ! 2*masked bottom density gradient
            zgdrho = (  za * ( zts(ji+1,jj,jp_tem) - zts(ji,jj,jp_tem) )    &
               &      - zb * ( zts(ji+1,jj,jp_sal) - zts(ji,jj,jp_sal) )  ) * umask(ji,jj,1)
            !
            zsign  = SIGN(  0.5_wp, -zgdrho * REAL( mgrhu(ji,jj) )  )    ! sign of ( i-gradient * i-slope )
            ahu_bbl(ji,jj) = ( 0.5 - zsign ) * ahu_bbl_0(ji,jj)       ! masked diffusive flux coeff.
            !
            !                                                   ! j-direction
            za = zab(ji,jj+1,jp_tem) + zab(ji,jj,jp_tem)              ! 2*(alpha,beta) at v-point
            zb = zab(ji,jj+1,jp_sal) + zab(ji,jj,jp_sal)
            !                                                         ! 2*masked bottom density gradient
            zgdrho = (  za * ( zts(ji,jj+1,jp_tem) - zts(ji,jj,jp_tem) )    &
               &      - zb * ( zts(ji,jj+1,jp_sal) - zts(ji,jj,jp_sal) )  ) * vmask(ji,jj,1)
            !
            zsign = SIGN(  0.5_wp, -zgdrho * REAL( mgrhv(ji,jj) )  )     ! sign of ( j-gradient * j-slope )
            ahv_bbl(ji,jj) = ( 0.5 - zsign ) * ahv_bbl_0(ji,jj)
         END_2D
         !
      ENDIF
      !
      !                                   !-------------------!
      IF( nn_bbl_adv /= 0 ) THEN          !   advective bbl   !
         !                                !-------------------!
         SELECT CASE ( nn_bbl_adv )             !* bbl transport type
         !
         CASE( 1 )                                   != use of upper velocity
            DO_2D_OVR( 1, 0, 1, 0 )                              ! criteria: grad(rho).grad(h)<0  and grad(rho).grad(h)<0
               !                                                  ! i-direction
               za = zab(ji+1,jj,jp_tem) + zab(ji,jj,jp_tem)               ! 2*(alpha,beta) at u-point
               zb = zab(ji+1,jj,jp_sal) + zab(ji,jj,jp_sal)
               !                                                          ! 2*masked bottom density gradient
               zgdrho = (  za * ( zts(ji+1,jj,jp_tem) - zts(ji,jj,jp_tem) )    &
                         - zb * ( zts(ji+1,jj,jp_sal) - zts(ji,jj,jp_sal) )  ) * umask(ji,jj,1)
               !
               zsign = SIGN(  0.5_wp, - zgdrho   * REAL( mgrhu(ji,jj) )  )   ! sign of i-gradient * i-slope
               zsigna= SIGN(  0.5_wp, zub(ji,jj) * REAL( mgrhu(ji,jj) )  )   ! sign of u * i-slope
               !
               !                                                          ! bbl velocity
               utr_bbl(ji,jj) = ( 0.5 + zsigna ) * ( 0.5 - zsign ) * e2u(ji,jj) * e3u_bbl_0(ji,jj) * zub(ji,jj)
               !
               !                                                  ! j-direction
               za = zab(ji,jj+1,jp_tem) + zab(ji,jj,jp_tem)               ! 2*(alpha,beta) at v-point
               zb = zab(ji,jj+1,jp_sal) + zab(ji,jj,jp_sal)
               !                                                          ! 2*masked bottom density gradient
               zgdrho = (  za * ( zts(ji,jj+1,jp_tem) - zts(ji,jj,jp_tem) )    &
                  &      - zb * ( zts(ji,jj+1,jp_sal) - zts(ji,jj,jp_sal) )  ) * vmask(ji,jj,1)
               zsign = SIGN(  0.5_wp, - zgdrho   * REAL( mgrhv(ji,jj) )  )   ! sign of j-gradient * j-slope
               zsigna= SIGN(  0.5_wp, zvb(ji,jj) * REAL( mgrhv(ji,jj) )  )   ! sign of u * i-slope
               !
               !                                                          ! bbl transport
               vtr_bbl(ji,jj) = ( 0.5 + zsigna ) * ( 0.5 - zsign ) * e1v(ji,jj) * e3v_bbl_0(ji,jj) * zvb(ji,jj)
            END_2D
            !
         CASE( 2 )                                 != bbl velocity = F( delta rho )
            zgbbl = grav * rn_gambbl
            DO_2D_OVR( 1, 0, 1, 0 )                         ! criteria: rho_up > rho_down
               !                                                  ! i-direction
               ! down-slope T-point i/k-index (deep)  &   up-slope T-point i/k-index (shelf)
               iid  = ji + MAX( 0, mgrhu(ji,jj) )
               iis  = ji + 1 - MAX( 0, mgrhu(ji,jj) )
               !
               ikud = mbku_d(ji,jj)
               ikus = mbku(ji,jj)
               !
               za = zab(ji+1,jj,jp_tem) + zab(ji,jj,jp_tem)               ! 2*(alpha,beta) at u-point
               zb = zab(ji+1,jj,jp_sal) + zab(ji,jj,jp_sal)
               !                                                          !   masked bottom density gradient
               zgdrho = 0.5 * (  za * ( zts(iid,jj,jp_tem) - zts(iis,jj,jp_tem) )    &
                  &            - zb * ( zts(iid,jj,jp_sal) - zts(iis,jj,jp_sal) )  ) * umask(ji,jj,1)
               zgdrho = MAX( 0.e0, zgdrho )                               ! only if shelf is denser than deep
               !
               !                                                          ! bbl transport (down-slope direction)
               utr_bbl(ji,jj) = e2u(ji,jj) * e3u_bbl_0(ji,jj) * zgbbl * zgdrho * REAL( mgrhu(ji,jj) )
               !
               !                                                  ! j-direction
               !  down-slope T-point j/k-index (deep)  &   of the up  -slope T-point j/k-index (shelf)
               ijd  = jj + MAX( 0, mgrhv(ji,jj) )
               ijs  = jj + 1 - MAX( 0, mgrhv(ji,jj) )
               !
               ikvd = mbkv_d(ji,jj)
               ikvs = mbkv(ji,jj)
               !
               za = zab(ji,jj+1,jp_tem) + zab(ji,jj,jp_tem)               ! 2*(alpha,beta) at v-point
               zb = zab(ji,jj+1,jp_sal) + zab(ji,jj,jp_sal)
               !                                                          !   masked bottom density gradient
               zgdrho = 0.5 * (  za * ( zts(ji,ijd,jp_tem) - zts(ji,ijs,jp_tem) )    &
                  &            - zb * ( zts(ji,ijd,jp_sal) - zts(ji,ijs,jp_sal) )  ) * vmask(ji,jj,1)
               zgdrho = MAX( 0.e0, zgdrho )                               ! only if shelf is denser than deep
               !
               !                                                          ! bbl transport (down-slope direction)
               vtr_bbl(ji,jj) = e1v(ji,jj) * e3v_bbl_0(ji,jj) * zgbbl * zgdrho * REAL( mgrhv(ji,jj) )
            END_2D
         END SELECT
         !
      ENDIF
      !
   END SUBROUTINE bbl


   SUBROUTINE tra_bbl_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_bbl_init  ***
      !!
      !! ** Purpose :   Initialization for the bottom boundary layer scheme.
      !!
      !! ** Method  :   Read the nambbl namelist and check the parameters
      !!              called by nemo_init at the first timestep (kit000)
      !!----------------------------------------------------------------------
      INTEGER ::   ji, jj                      ! dummy loop indices
      INTEGER ::   ii0, ii1, ij0, ij1, ios     ! local integer
      REAL(wp), DIMENSION(jpi,jpj) ::   zmbku, zmbkv   ! workspace
      !!
      NAMELIST/nambbl/ ln_trabbl, nn_bbl_ldf, nn_bbl_adv, rn_ahtbbl, rn_gambbl
#if defined key_drakkar
      NAMELIST/nambbl_drk/ ln_kriteria
#endif
      !!----------------------------------------------------------------------
      !
      READ  ( numnam_ref, nambbl, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'nambbl in reference namelist' )
      !
      READ  ( numnam_cfg, nambbl, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'nambbl in configuration namelist' )
      IF(lwm) WRITE ( numond, nambbl )
      !
#if defined key_drakkar
      READ  ( numnam_ref, nambbl_drk, IOSTAT = ios, ERR = 903)
903   IF( ios /= 0 )   CALL ctl_nam ( ios , 'nambbl_drk in reference namelist' )
      !
      READ  ( numnam_cfg, nambbl_drk, IOSTAT = ios, ERR = 904 )
904   IF( ios >  0 )   CALL ctl_nam ( ios , 'nambbl in configuration namelist' )
      IF(lwm) WRITE ( numond, nambbl_drk )
      !
#endif
      l_bbl = .TRUE.                 !* flag to compute bbl coef and transport
      !
      IF(lwp) THEN                   !* Parameter control and print
         WRITE(numout,*)
         WRITE(numout,*) 'tra_bbl_init : bottom boundary layer initialisation'
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '       Namelist nambbl : set bbl parameters'
         WRITE(numout,*) '          bottom boundary layer flag          ln_trabbl  = ', ln_trabbl
      ENDIF
      IF( .NOT.ln_trabbl )   RETURN
      !
      IF(lwp) THEN
         WRITE(numout,*) '          diffusive bbl (=1)   or not (=0)    nn_bbl_ldf = ', nn_bbl_ldf
         WRITE(numout,*) '          advective bbl (=1/2) or not (=0)    nn_bbl_adv = ', nn_bbl_adv
         WRITE(numout,*) '          diffusive bbl coefficient           rn_ahtbbl  = ', rn_ahtbbl, ' m2/s'
         WRITE(numout,*) '          advective bbl coefficient           rn_gambbl  = ', rn_gambbl, ' s'
      ENDIF
      !
      !                              ! allocate trabbl arrays
      IF( tra_bbl_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'tra_bbl_init : unable to allocate arrays' )
      !
      IF(lwp) THEN
         IF( nn_bbl_adv == 1 )    WRITE(numout,*) '       * Advective BBL using upper velocity'
         IF( nn_bbl_adv == 2 )    WRITE(numout,*) '       * Advective BBL using velocity = F( delta rho)'
      ENDIF
      !
      !                             !* vertical index of  "deep" bottom u- and v-points
      DO_2D( 1, 0, 1, 0 )                 ! (the "shelf" bottom k-indices are mbku and mbkv)
         mbku_d(ji,jj) = MAX(  mbkt(ji+1,jj  ) , mbkt(ji,jj)  )   ! >= 1 as mbkt=1 over land
         mbkv_d(ji,jj) = MAX(  mbkt(ji  ,jj+1) , mbkt(ji,jj)  )
      END_2D
      ! converte into REAL to use lbc_lnk ; impose a min value of 1 as a zero can be set in lbclnk
      zmbku(:,:) = REAL( mbku_d(:,:), wp )   ;     zmbkv(:,:) = REAL( mbkv_d(:,:), wp )
      CALL lbc_lnk( 'trabbl', zmbku,'U',1.0_wp, zmbkv,'V',1.0_wp)
      mbku_d(:,:) = MAX( INT( zmbku(:,:) ), 1 ) ;  mbkv_d(:,:) = MAX( NINT( zmbkv(:,:) ), 1 )
      !
      !                             !* sign of grad(H) at u- and v-points; zero if grad(H) = 0
      mgrhu(:,:) = 0   ;   mgrhv(:,:) = 0
#if defined key_drakkar
          IF ( ln_kriteria ) THEN
            DO_2D( 1, 0, 1, 0 )
              mgrhu(ji,jj) = mbkt(ji+1,jj) - mbkt(ji,jj)
              mgrhv(ji,jj) = mbkt(ji,jj+1) - mbkt(ji,jj)
              IF ( mgrhu(ji,jj) /= 0 ) THEN 
                  mgrhu(ji,jj) = INT(  SIGN( 1.e0, FLOAT(mbkt(ji+1,jj) - mbkt(ji,jj)) )) 
              ENDIF
              IF ( mgrhv(ji,jj) /= 0 ) THEN 
                  mgrhv(ji,jj) = INT(  SIGN( 1.e0, FLOAT(mbkt(ji,jj+1) - mbkt(ji,jj)) ))
              ENDIF  
            END_2D
          ELSE  ! ln_kriteria
#endif
      DO_2D( 1, 0, 1, 0 )
         IF( gdept_0(ji+1,jj,mbkt(ji+1,jj)) - gdept_0(ji,jj,mbkt(ji,jj)) /= 0._wp ) THEN
            mgrhu(ji,jj) = INT(  SIGN( 1.0_wp, gdept_0(ji+1,jj,mbkt(ji+1,jj)) - gdept_0(ji,jj,mbkt(ji,jj)) )  )
         ENDIF
         !
         IF( gdept_0(ji,jj+1,mbkt(ji,jj+1)) - gdept_0(ji,jj,mbkt(ji,jj)) /= 0._wp ) THEN
            mgrhv(ji,jj) = INT(  SIGN( 1.0_wp, gdept_0(ji,jj+1,mbkt(ji,jj+1)) - gdept_0(ji,jj,mbkt(ji,jj)) )  )
         ENDIF
      END_2D
#if defined key_drakkar
          ENDIF  ! ln_kriteria
#endif
      !
      DO_2D( 1, 0, 1, 0 )           !* bbl thickness at u- (v-) point; minimum of top & bottom e3u_0 (e3v_0)
         e3u_bbl_0(ji,jj) = MIN( e3u_0(ji,jj,mbkt(ji+1,jj  )), e3u_0(ji,jj,mbkt(ji,jj)) )
         e3v_bbl_0(ji,jj) = MIN( e3v_0(ji,jj,mbkt(ji  ,jj+1)), e3v_0(ji,jj,mbkt(ji,jj)) )
      END_2D
      CALL lbc_lnk( 'trabbl', e3u_bbl_0, 'U', 1.0_wp , e3v_bbl_0, 'V', 1.0_wp )      ! lateral boundary conditions
      !
      !                             !* masked diffusive flux coefficients
      ahu_bbl_0(:,:) = rn_ahtbbl * e2_e1u(:,:) * e3u_bbl_0(:,:) * umask(:,:,1)
      ahv_bbl_0(:,:) = rn_ahtbbl * e1_e2v(:,:) * e3v_bbl_0(:,:) * vmask(:,:,1)
      !
   END SUBROUTINE tra_bbl_init

   !!======================================================================
END MODULE trabbl
