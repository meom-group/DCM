PROGRAM datfin
  !!----------------------------------------------------------------------
  !!               ***    PROGRAM DATFIN   ***
  !!
  !!  ** Purpose:  This PROGRAM takes a date and a number of days, it returns
  !!              the date after the specified number of days.
  !!              Date FORMAT is following the OPA ndastp variable, eg 851028
  !!              stands for October 28th 1985.
  !!
  !!  ** Method:
  !!      1- Decode the OPA FORMAT and build the date as iyyy mm dd
  !!      2- Compute the corresponding julian day
  !!      3- Add the  specified number of day
  !!      4- Convert the final date in julian day in iyyy mm dd FORMAT
  !!      5- PRINT the final date in the OPA FORMAT
  !!
  !!  ** Credit:
  !!      The original sources for julday and caldat were found in the
  !!      book Numerical Recipes.
  !!      
  !!  history
  !!        1998 : porting to OPA form by J.M. Molines
  !!        2007 : in fortran 90 by J.M. Molines
  !!
  !!-----------------------------------------------------------------------
  ! *  Declarations
  !
  IMPLICIT NONE
  !
  INTEGER :: ndastp, njour, iyyy, imm, id, juldeb, julfin
  INTEGER :: jfin,mfin,ianfin, mm
  INTEGER :: ichar, iargc,  narg, ILEN, ierr
  CHARACTER(LEN=80) :: line
  !
  ! *  red  arguments 
  narg=iargc()
  IF (narg /=  2 ) THEN
     PRINT *,' usage: datfin <ndastp> <ndays>'
     PRINT *,' Return the date (in ndastp format) of the'
     PRINT *,' day, ndays after the original ndastp date'
     STOP
  END IF
  CALL getarg(1,line)
  READ(line,*) ndastp
  CALL getarg(2,line)
  READ(line,*) njour
  !
  iyyy = ndastp/10000 
  mm = (ndastp - iyyy*10000)/100
  id = (ndastp - iyyy*10000 - mm*100)
  !        iyyy = iyyy + 1900
  !
  juldeb = julday(mm,id,iyyy) 
  julfin= juldeb + njour
  CALL caldat (julfin,mfin,jfin,ianfin)
  PRINT '(i4.2,2i2.2)',ianfin,mfin,jfin
  !	PRINT '(3i2.2)',ianfin-1900,mfin,jfin

CONTAINS
  FUNCTION julday(kmm,kid,kiyyy)
  !! ------------------------------------------------------------------
  !!          ***        FUNCTION JULDAY    ***
  !!
  !!   Purpose:   This routine returns the julian day number which begins at noon
  !!         of the calendar date specified by month kmm, day kid, and year kiyyy.
  !!         positive year signifies a.d.; negative, b.c.  (remember that the
  !!         year after 1 b.c. was 1 a.d.)
  !!         routine handles changeover to gregorian calendar on oct. 15, 1582.
  !!
  !!   Method:  This routine comes directly from the Numerical Recipe Book,
  !!           press et al., numerical recipes, cambridge univ. press, 1986.
  !!
  !!   Arguments:
  !!     kmm     : input, corresponding month
  !!     kid     : input, corresponding day
  !!     kiyyy   : input, corresponding year, positive IF a.d, negative b.c.
  !!      
  !!     
  !!   history
  !!     1998: J.M. Molines for the Doctor form. 
  !!     2007 : J.M. Molines in F90
  !! -----------------------------------------------------------------
   !  *  Declarations
   !
   INTEGER :: julday, kiyyy
   INTEGER, INTENT(in)  ::kmm, kid
   !  * Local 
   INTEGER, PARAMETER ::jpgreg=15+31*(10+12*1582)
   INTEGER  :: iy, im, ia
    ! ... Year 0 never existed ...
    IF (kiyyy == 0) STOP 101
    !
    IF (kiyyy < 0) kiyyy = kiyyy + 1
    IF (kmm > 2) THEN
       iy = kiyyy
       im = kmm + 1
    ELSE
       iy = kiyyy - 1
       im = kmm + 13
    END IF
    !
    julday = INT(365.25*iy) + INT(30.6001*im) + kid + 1720995 
    IF (kid+31*(kmm+12*kiyyy).GE.jpgreg) THEN
       ia = INT(0.01*iy)
       julday = julday + 2 - ia + INT(0.25*ia) 
    END IF
  END FUNCTION JULDAY

  SUBROUTINE caldat(kjulian,kmm,kid,kiyyy)
   !! -------------------------------------------------------------------
   !!             ***      SUBROUTINE CALDAT  ***
   !!
   !!   Purpose:  This routine convert a julian day in calendar date.
   !!             It is the inverse of the FUNCTION julday.  
   !!
   !!   Method:   This routine comes directly from the Numerical Recipe Book,
   !!           press et al., numerical recipes, cambridge univ. press, 1986.
   !!
   !!   Arguments
   !!     kjulian : input julian day number
   !!     kmm     : output, corresponding month
   !!     kid     : output, corresponding day
   !!     kiyyy   : output, corresponding year, positive IF a.d, negative b.c.
   !!      
   !!    
   !!   history
   !!     1998: J.M. Molines for the Doctor form.
   !!     2007: J.M. Molines in F90
   !!------------------------------------------------------------------------
    ! *  Declarations:
    !
    IMPLICIT NONE
    INTEGER,INTENT(in) ::  kjulian
    INTEGER,INTENT(out) ::  kmm, kid, kiyyy
    ! * Local
    INTEGER, PARAMETER ::jpgreg = 2299161
    INTEGER  :: ia,ialpha,ib,ic,id,ie
    !
    ! ... Cross over to Gregorian Calendar produces this correction:
    !
    IF ( kjulian >= jpgreg) THEN
       ialpha = INT ((( kjulian - 1867216) - 0.25)/36524.25 )
       ia     = kjulian +1 + ialpha -INT (0.25*ialpha)
    ELSE
       ia = kjulian
    END IF
    !
    ib = ia + 1524
    ic = INT (6680. + (( ib -2439870) - 122.1)/365.25 )
    id = 365* ic + INT (0.25*ic)
    ie = INT (( ib - id )/30.6001)
    !
    kid = ib - id - INT (30.6001*ie)
    kmm = ie -1
    IF ( kmm > 12 ) kmm = kmm - 12
    kiyyy = ic - 4715
    IF ( kmm   >  2 ) kiyyy = kiyyy - 1
    IF ( kiyyy <= 0 ) kiyyy = kiyyy - 1
  END SUBROUTINE caldat
END PROGRAM datfin
