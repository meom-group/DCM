MODULE diaprod
! Requires key_drakkar
# if defined key_drakkar
! Requires key_iom_put
# if defined key_xios
   !!======================================================================
   !!                     ***  MODULE  diaprod  ***
   !! Ocean diagnostics :  write ocean product diagnostics
   !!=====================================================================
   !! History :  3.4  ! 2012  (D. Storkey)  Original code
   !!            4.0  ! 2019  (D. Storkey)
   !!            4.2  ! 2022  (J.M. Molines)
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dia_prod      : calculate and write out product diagnostics
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers 
   USE dom_oce         ! ocean space and time domain
   USE iom             
   USE timing          ! performance summary

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dia_prod                 ! routines called by step.F90

   !! * Substitutions
#  include "do_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.2 , NEMO Consortium (2022)
   !! $Id: diaprod.F90 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dia_prod( kt ,Kmm)
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE dia_prod  ***
      !!                   
      !! ** Purpose :   Write out product diagnostics (uT, vS etc.)
      !!
      !! ** Method  :  use iom_put
      !!               Product diagnostics are not thickness-weighted in
      !this routine.
      !!               They should be thickness-weighted using XIOS if
      !key_vvl is set. 
      !!----------------------------------------------------------------------
      !!
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
      INTEGER, INTENT( in ) ::   Kmm     ! Ocean time-level index
      !!
      !!----------------------------------------------------------------------
      INTEGER                      ::   ji, jj, jk              ! dummy loop indices
      REAL(wp)                     ::   zztmp, zztmpx, zztmpy   ! 
      !!
      REAL(wp), POINTER, DIMENSION(:,:)   :: z2d      ! 2D workspace
      REAL(wp), POINTER, DIMENSION(:,:,:) :: z3d      ! 3D workspace
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zrhop    ! potential density
      !!----------------------------------------------------------------------
      ! 
      IF( ln_timing )   CALL timing_start('dia_prod')
      ! 
      ALLOCATE( z2d(jpi,jpj), z3d(jpi,jpj,jpk), zrhop(jpi,jpj,jpk) )
      !

      IF( iom_use("urhop") .OR. iom_use("vrhop") .OR. iom_use("wrhop") ) THEN 
         DO_3D(nn_hls, nn_hls, nn_hls, nn_hls,1,jpk)
           zrhop(ji,jj,jk) = rhop(ji,jj,jk)-1000.e0*tmask(ji,jj,jk)         ! reference potential density to 1000 to avoid precision issues in rhop2 calculation
         END_3D
      ENDIF

      IF( iom_use("ut") ) THEN
         z3d(:,:,:) = 0.e0 
            DO_3D(0,0,0,0,1,jpkm1)
                  z3d(ji,jj,jk) = uu(ji,jj,jk,Kmm) * 0.5 * ( ts(ji,jj,jk,jp_tem,Kmm) + ts(ji+1,jj,jk,jp_tem,Kmm) )
            END_3D
         CALL iom_put( "ut", z3d )                  ! product of temperature and zonal velocity at U points
      ENDIF

      IF( iom_use("vt") ) THEN
         z3d(:,:,:) = 0.e0
            DO_3D(0,0,0,0,1,jpkm1)
                  z3d(ji,jj,jk) = vv(ji,jj,jk,Kmm) * 0.5 * ( ts(ji,jj,jk,jp_tem,Kmm) + ts(ji,jj+1,jk,jp_tem,Kmm) )
            END_3D
         CALL iom_put( "vt", z3d )                  ! product of temperature and zonal velocity at U points
      ENDIF

      IF( iom_use("wt") ) THEN
         z3d(:,:,:) = 0.e0
         DO_2D(0,0,0,0)
             z3d(ji,jj,1) = ww(ji,jj,1) * ts(ji,jj,1,jp_tem,Kmm)
         END_2D
            DO_3D(0,0,0,0,2,jpkm1)
                  z3d(ji,jj,jk) = ww(ji,jj,jk) * 0.5 * ( ts(ji,jj,jk-1,jp_tem,Kmm) + ts(ji,jj,jk,jp_tem,Kmm) )
            END_3D
         CALL iom_put( "wt", z3d )                  ! product of temperature and zonal velocity at U points
      ENDIF

      IF( iom_use("us") ) THEN
         z3d(:,:,:) = 0.e0
            DO_3D(0,0,0,0,1,jpkm1)
                  z3d(ji,jj,jk) = uu(ji,jj,jk,Kmm) * 0.5 * ( ts(ji,jj,jk,jp_sal,Kmm) + ts(ji+1,jj,jk,jp_sal,Kmm) )
            END_3D
         CALL iom_put( "us", z3d )                  ! product of temperature and zonal velocity at U points
      ENDIF

      IF( iom_use("vs") ) THEN
         z3d(:,:,:) = 0.e0
            DO_3D(0,0,0,0,1,jpkm1)
                  z3d(ji,jj,jk) = vv(ji,jj,jk,Kmm) * 0.5 * ( ts(ji,jj,jk,jp_sal,Kmm) + ts(ji,jj+1,jk,jp_sal,Kmm) )
            END_3D
         CALL iom_put( "vs", z3d )                  ! product of temperature and zonal velocity at U points
      ENDIF

      IF( iom_use("ws") ) THEN
         z3d(:,:,:) = 0.e0
         DO_2D(0,0,0,0)
             z3d(ji,jj,1) = ww(ji,jj,1) * ts(ji,jj,1,jp_sal,Kmm)
         END_2D
            DO_3D(0,0,0,0,2,jpkm1)
                  z3d(ji,jj,jk) = ww(ji,jj,jk) * 0.5 * ( ts(ji,jj,jk-1,jp_sal,Kmm) + ts(ji,jj,jk,jp_sal,Kmm) )
            END_3D
         CALL iom_put( "ws", z3d )                  ! product of temperature and zonal velocity at U points
      ENDIF

      IF( iom_use("uu") ) THEN
         z3d(:,:,:) = 0.e0
            DO_3D(nn_hls,nn_hls,nn_hls,nn_hls,1,jpkm1)
                  z3d(ji,jj,jk) = uu(ji,jj,jk,Kmm) * uu(ji,jj,jk,Kmm)
            END_3D
         CALL iom_put( "uu", z3d )                  ! product of temperature and zonal velocity at U points
      ENDIF

      IF( iom_use("uv") ) THEN
         z3d(:,:,:) = 0.e0
            DO_3D(0,0,0,0,1,jpkm1)
                  z3d(ji,jj,jk) = 0.25 *( uu(ji-1,jj,jk,Kmm)+ uu(ji,jj,jk,Kmm) ) * ( vv(ji,jj-1,jk,Kmm)+ vv(ji,jj,jk,Kmm) )
            END_3D
         CALL iom_put( "uv", z3d )                  ! product of temperature and zonal velocity at U points
      ENDIF

      IF( iom_use("uw") ) THEN
         z3d(:,:,:) = 0.e0 
         DO_2D(0,0,0,0)
             z3d(ji,jj,1) = 0.5 * ( ww(ji,jj,1) + ww(ji+1,jj,1) ) * uu(ji,jj,1,Kmm)
         END_2D

         DO_3D(0,0,0,0,2,jpkm1)
               z3d(ji,jj,jk) = 0.25 * ( ww(ji,jj,jk) + ww(ji+1,jj,jk) ) * ( uu(ji,jj,jk-1,Kmm) + uu(ji,jj,jk,Kmm) )
         END_3D
         CALL iom_put( "uw", z3d )                  ! product of velocity and vertical velocity at UW points
      ENDIF

      IF( iom_use("vv") ) THEN
         z3d(:,:,:) = 0.e0
            DO_3D(nn_hls,nn_hls,nn_hls,nn_hls,1,jpkm1)
                  z3d(ji,jj,jk) = vv(ji,jj,jk,Kmm) * vv(ji,jj,jk,Kmm)
            END_3D
         CALL iom_put( "vv", z3d )                  ! product of temperature and zonal velocity at U points
      ENDIF

      IF( iom_use("vw") ) THEN
         z3d(:,:,:) = 0.e0
         DO_2D(0,0,0,0)
             z3d(ji,jj,1) = 0.5 * ( ww(ji,jj,1) + ww(ji,jj+1,1) ) * vv(ji,jj,1,Kmm)
         END_2D

         DO_3D(0,0,0,0,2,jpkm1)
               z3d(ji,jj,jk) = 0.25 * ( ww(ji,jj,jk) + ww(ji,jj+1,jk) ) * ( vv(ji,jj,jk-1,Kmm) + vv(ji,jj,jk,Kmm) )
         END_3D
         CALL iom_put( "vw", z3d )                  ! product of velocity and vertical velocity at UW points
      ENDIF

      IF( iom_use("urhop") ) THEN
         z3d(:,:,:) = 0.e0 
         DO_3D(0,0,0,0,1,jpkm1)
                  z3d(ji,jj,jk) = uu(ji,jj,jk,Kmm) * 0.5 * ( zrhop(ji,jj,jk) + zrhop(ji+1,jj,jk) )
         END_3D
         CALL iom_put( "urhop", z3d )                  ! product density and zonal velocity at U points
      ENDIF

      IF( iom_use("vrhop") ) THEN
         z3d(:,:,:) = 0.e0 
         DO_3D(0,0,0,0,1,jpkm1)
                  z3d(ji,jj,jk) = vv(ji,jj,jk,Kmm) * 0.5 * ( zrhop(ji,jj,jk) + zrhop(ji,jj+1,jk) )
         END_3D
         CALL iom_put( "vrhop", z3d )                  ! product density and zonal velocity at U points
      ENDIF

      IF( iom_use("wrhop") ) THEN
         z3d(:,:,:) = 0.e0 
         DO_2D(0,0,0,0)
             z3d(ji,jj,1) = ww(ji,jj,1) * zrhop(ji,jj,1)
         END_2D
        
         DO_3D(0,0,0,0,2,jpkm1)
                  z3d(ji,jj,jk) = ww(ji,jj,jk) * 0.5 * ( zrhop(ji,jj,jk-1) + zrhop(ji,jj,jk) )
         END_3D
         CALL iom_put( "wrhop", z3d )                  ! product density and vertical velocity at W points
      ENDIF

      !
      DEALLOCATE( z2d, z3d, zrhop )
      !
      IF( ln_timing )   CALL timing_stop('dia_prod')
      !
   END SUBROUTINE dia_prod
#else
   !!----------------------------------------------------------------------
   !!   Default option :                                         NO
   !diaprod
   !!----------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER :: lk_diaprod = .FALSE.   ! coupled flag
CONTAINS
   SUBROUTINE dia_prod( kt )   ! Empty routine
      INTEGER ::   kt
      WRITE(*,*) 'dia_prod: You should not have seen this print! error?', kt
   END SUBROUTINE dia_prod
#endif
#else
   !!----------------------------------------------------------------------
   !!   Default option :                                         NO
   !diaprod : key_drakkar is required for using diaprod
   !!----------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER :: lk_diaprod = .FALSE.   ! coupled flag
CONTAINS
   SUBROUTINE dia_prod( kt )   ! Empty routine
      INTEGER ::   kt
      WRITE(*,*) 'dia_prod: You should not have seen this print! error?', kt
   END SUBROUTINE dia_prod
#endif

   !!======================================================================
END MODULE diaprod

