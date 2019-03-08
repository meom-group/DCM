PROGRAM linreg
  ! This program return the slope of linear regression from X Y points given in a zrate file 
  ! $Date: 2011-12-07 21:57:43 +0100 (Wed, 07 Dec 2011) $
  ! $Rev: 353 $
  ! $Id: linreg.f90 353 2011-12-07 20:57:43Z molines $
  !-----------------------------------------------
  IMPLICIT NONE

  CHARACTER(LEN=80) :: cfile='zrate', cdum
  INTEGER :: narg, ji, np, ibeg, iend, numin=10
  INTEGER :: iargc
  REAL(KIND=4), DIMENSION(:), ALLOCATABLE :: x, y
  REAL(KIND=8) :: xybar, xbar, ybar, x2bar, y2bar, sxy, sx2, sy2
  REAL(KIND=8) :: rate
  
  narg=iargc()

  IF ( narg /= 2 ) THEN
    PRINT *,' USAGE: linreg ibeg iend [filename] '
    PRINT *,'    This program return the slope of linear regression '
    PRINT *,'    from X Y points given in a zrate file '
    PRINT *,'    The input filename (x y file) can be provided as 3rd argument'
    STOP
  ENDIF
  
  CALL getarg(1, cdum) ; READ(cdum,*) ibeg
  CALL getarg(2, cdum) ; READ(cdum,*) iend
  IF ( narg == 3 ) CALL getarg(3, cfile)

  OPEN(numin, FILE=cfile)
  np=-1
  DO 
    READ(numin,*, END=999)
    np=np+1
  ENDDO
999  print *, np,' lines in ', TRIM(cfile)
  ALLOCATE(x(np), y(np) )
  REWIND(numin)

  DO ji=1,np
    READ(numin,*) x(ji), y(ji)
  ENDDO
  xybar=0.d0
  xbar=0.d0
  ybar=0.d0
  x2bar=0.d0
  y2bar=0.d0

  np=iend-ibeg+1
  DO ji=ibeg, iend
    xybar=xybar+x(ji)*y(ji)/np
    xbar=xbar+x(ji)/np
    ybar=ybar+y(ji)/np
    x2bar=x2bar+x(ji)*x(ji)/np
    y2bar=y2bar+y(ji)*y(ji)/np
  ENDDO
    sxy=xybar-xbar*ybar
    sx2=x2bar-xbar*xbar
    sy2=y2bar-ybar*ybar
   
    rate=sxy/sx2
    PRINT *, 'RATE :', rate , 'step/min'
    
END PROGRAM linreg
