MODULE diaprod
! Requires key_drakkar
# if defined key_drakkar
! Requires key_iom_put
# if defined key_iomput
   !!----------------------------------------------------------------------
   !!   Default option :                                         NO
   !diaprod
   !!----------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER :: lk_diaprod = .FALSE.   ! coupled flag
   PUBLIC   dia_prod                 ! routines called by step.F90

CONTAINS
   SUBROUTINE dia_prod( kt )   ! Empty routine
      INTEGER ::   kt
      !  Diaprod not yet implemented in 4.2RC
      ! this empty routine is there to solve  compilation issues.
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

