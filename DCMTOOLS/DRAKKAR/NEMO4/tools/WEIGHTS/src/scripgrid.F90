! ==============================================================================
# define ARGC

PROGRAM scripgrid

  USE scripgrid_mod

  CHARACTER(char_len) :: infile

#if defined ARGC
  INTEGER :: IARGC
  EXTERNAL IARGC

  if (IARGC() == 1) then
    CALL GETARG(1, infile)
    CALL convert( infile )
  ELSE
    write(6,*) 'need to supply a namelist file'
  ENDIF
#else
  write(6,*) 'enter name of namelist file'
  read(5,*) infile

  CALL convert( infile )
#endif

END PROGRAM scripgrid
