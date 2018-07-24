#define ARGC
! ==========================================================================

program scripinterp

  use scripinterp_mod

  character (char_len) ::  nm_in
#if defined ARGC
  integer :: iargc

  if (iargc() == 1) then
    call getarg(1, nm_in)
    call process_grid(nm_in)
  else
    write(6,*) 'need the name of an input namelist'
  endif
#else
  write(6,*) 'enter the name of an input namelist'
  read(5,*) nm_in
  call process_grid(nm_in)
#endif

end program scripinterp
