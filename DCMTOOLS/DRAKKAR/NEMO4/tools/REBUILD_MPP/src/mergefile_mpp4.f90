PROGRAM mergefile_mpp4
  !!======================================================================
  !!                     ***  PROGRAM  mergefile_mpp4 ***
  !!=====================================================================
  !!  ** Purpose : Merge a bunch of splitted files into a single
  !!               one. The merge function work with MPP nemo output.
  !!               This version must be run in // : each task merge one or
  !!               several files.
  !!
  !!  ** Method  :  Look for the description of the file set into
  !!                rank 0000 file. Then build a clone of the 0000
  !!                file except for the dimensions. Fill the merged
  !!                file by reading individual subdomains.
  !!
  !! History : 1.0 : 07-30-2014 : J.M. Molines Original code
  !!----------------------------------------------------------------------
  !!----------------------------------------------------------------------
  !!   routines                   : description
  !!   CreateMergedFile             Create the merged file
  !!   FillMergedFile               Fill the merged file
  !!
  !!   functions                  : description
  !!   GetNcFile                  : Get meta data from an ncfile into a structure
  !!   GetCoord                     Read coordinate.nc file
  !!   PrintNcFile                : Utility for dumping file structure
  !!   RenameFile                 : Utility for renaming file to DRAKKAR standards (-r option)
  !!----------------------------------------------------------------------
  USE netcdf
  USE modncfile
  USE modutil
  USE mpi
  !!----------------------------------------------------------------------
  !! DCM, MEOM 2014
  !! $Id$
  !! Copyright (c) 2014, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !!----------------------------------------------------------------------
  IMPLICIT NONE
  INTEGER(KIND=4) :: jf                     ! loop index
  INTEGER(KIND=4) :: narg, iargc, ijarg     ! browse command line variable
  INTEGER(KIND=4) :: iargfile, nfile        ! for deciphering file list
  INTEGER(KIND=4) :: mpierr                 ! error status for MPI calls
  INTEGER(KIND=4) :: mpisize, mpirank       ! number of core and current rank
  INTEGER(KIND=4) :: ncurrent, iask, nstop  ! 
  INTEGER(KIND=4) :: ndone                  ! number of rebuilt file for each mpitask
  INTEGER(KIND=4) :: numlst=10              ! logical unit for the list file (-F)
  INTEGER(KIND=4) :: integer_size           ! 
  INTEGER(KIND=4) :: i_bsend_overhead       !
  INTEGER(KIND=4) :: i_msg_size = 1         ! number of integer to exchange with bsend
  INTEGER(KIND=4) :: iask_number = 1        ! dummy variable send to trigger answer from #0
  INTEGER(KIND=4) :: istopflag = 9999       ! sent by slave mpitask when finished
  INTEGER(KIND=4), DIMENSION(MPI_STATUS_SIZE)  :: mpistatus               !
  INTEGER(KIND=4), DIMENSION(:),   ALLOCATABLE :: i_bsend_buffer

  CHARACTER(LEN=80 ) :: cdum
  CHARACTER(LEN=20 ) :: cl_format

  !!----------------------------------------------------------------------
  narg=iargc()

  IF ( narg == 0 ) THEN
     PRINT *,' usage :  mergefile4 -F | -f <file*>_0000.nc  [-c coordinate_file ] ...'
     PRINT *,'          [-d output_directory] [-r ] [-nc3] [-v] [-kmax K-max] [-agrif]...'
     PRINT *,'          [-w imin imax jmin jmax] [-b BLOCK-number] [-encoding ndigit]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :' 
     PRINT *,'         This program recombine splitted netcd files into a single global file.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       Use either one of the following ( not both) :'
     PRINT *,'        -F : takes all possible _0000.nc files in current directory'
     PRINT *,'           This option has been introduced to overrid shell limitations '
     PRINT *,'           regarding the length of the arguments.'
     PRINT *,'       or :'
     PRINT *,'        -f <file*>_0000.nc : list of full name of rank 0 <file*> to be rebuild'
     PRINT *,' '
     PRINT *,'       If the number of digit used for coding the rank number in the file'
     PRINT *,'       names is not 4, condider the use of -encoding option.'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'        -c coordinates file : Use coordinate file to patch nav_lon,nav_lat'
     PRINT *,'              where there are land processors'
     PRINT *,'        -d output_dir       : Specify the output directory'
     PRINT *,'        -r  : Rename output file following DRAKKAR rules. Imply that the output'
     PRINT *,'              directory for each file is <freq>_OUTPUT (<freq> is 1d, 5d ... )'
     PRINT *,'              (-r option supersedes the -d option.) '
     PRINT *,'        -nc3: Create netcdf3 files instead of default netcdf4 '
     PRINT *,'        -v  : Verbose, more informations on output'
     PRINT *,'        -w imin imax jmin jmax : Only rebuild the sub-domain that are within '
     PRINT *,'              the specified windows (in model coordinates).  This is usefull'
     PRINT *,'              when rebuilding HUGE restart file having in mind to use only'
     PRINT *,'              a sub zone of the total domain. Using this option reduces the'
     PRINT *,'              time of the rebuild. However, the resulting file is still on'
     PRINT *,'              the total domain and will need some croping afterward.'
     PRINT *,'              Default behaviour is to rebuild the total domain.'
     PRINT *,'        -kmax K-max : define the maximum number of level to be processed in'
     PRINT *,'              the merged output file. (Usefull with -w option).'
     PRINT *,'        -agrif  : indicates to the program that your files are prefixed by '
     PRINT *,'              a number (ex: 1_xxx) '
     PRINT *,'        -b BLOCK-number : This option was written in order to reduce the '
     PRINT *,'              number of files opened simultaneously during the rebuild.'
     PRINT *,'              The total number of (mpp) files will be processed by blocks.'
     PRINT *,'              BLOCK-number is the number of blocks you choose. Default is 1.'
     PRINT *,'              This option can considerably reduce the memory imprint, linked'
     PRINT *,'              with the massive opening of all the mpp files together.'
     PRINT *,'              As a rule of thumb, BLOCK-number can be choosen in order to'
     PRINT *,'              block size ranging from 1000 to 2000. It of course depends on'
     PRINT *,'              the available memory.'
     PRINT *,'        -encoding ndigit : indicates the number of digits used for coding'
     PRINT *,'              the rank number in the mpp file''s names. Default is 4'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'        none'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf files : <file>_merg.nc, unless -r option is used.'
     PRINT *,'       variables    : same than in <file>_0000.nc file'
     PRINT *,'      '
     PRINT *,'     SEE ALSO :'
     PRINT *,'      rebuild_nemo, splitfile '
     PRINT *,'      '
     STOP
  ENDIF

  ! MPI STUFF initialisation.
  CALL MPI_INIT(mpierr)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, mpisize, mpierr)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, mpirank, mpierr)
  mmpirank = mpirank  ! mmpirank is a variable declared in modncfile even in case of mono proc program
  !
  ! Prepare BUFFER for BSEND, allowing space for ... 1 integer !
  CALL MPI_TYPE_SIZE ( MPI_INTEGER, integer_size, mpierr)
  i_bsend_overhead = INT(( MPI_BSEND_OVERHEAD *1.)/integer_size + 1 )
  ALLOCATE(i_bsend_buffer (i_msg_size + i_bsend_overhead) )
  CALL MPI_BUFFER_ATTACH ( i_bsend_buffer, integer_size*(i_msg_size + i_bsend_overhead), mpierr)

  ijarg=1
  DO WHILE ( ijarg <= narg ) 
     CALL getarg(ijarg, cdum) ; ijarg = ijarg + 1
     SELECT CASE (cdum)
     CASE ( '-f' )
        cdum=' '
        nfile = 0
        iargfile = ijarg      ! count number of files in input
        DO WHILE ( iargfile <= narg )
           CALL getarg (iargfile, cdum ) ; iargfile=iargfile+1
           IF ( cdum(1:1) /= '-' ) THEN
              nfile = nfile + 1
           ELSE
              EXIT
           ENDIF
        ENDDO

        ! allocate file list and re-read args
        ALLOCATE ( cf_list(nfile) )
        DO jf = 1, nfile
           CALL getarg(ijarg, cf_list(jf) ) ; ijarg = ijarg+1
        ENDDO
     CASE ('-F') ! internal build of the file list with ALL _0000.nc file in current dir
        IF ( mpirank == 0 ) THEN   ! only rank 0 is creating the .zmergfile_list.tmp file
          CALL SYSTEM( "find . -maxdepth 1 -name ""*"// TRIM(c_pattern) // "*"" -printf ""%f\n"" | sort -r  > .zmergfile_list.tmp")
        ENDIF
        ! Synchro to be sure that zmergfile_list.tmp has been written
        CALL MPI_BARRIER(MPI_COMM_WORLD, mpierr )
        nfile=0
        OPEN(numlst, FILE='.zmergfile_list.tmp')
        DO 
           READ(numlst,*,END=999)
           nfile=nfile+1
        ENDDO
999     PRINT *, 'Files to process : ', nfile
        ALLOCATE( cf_list(nfile) )
        REWIND(numlst)
        DO jf=1,nfile
           READ(numlst,"(a)")  cf_list(jf)
        ENDDO
        CLOSE(numlst)
     CASE ('-c')
        lg_coord      = .TRUE.
        lg_coord_read = .TRUE.
        CALL getarg(ijarg, cf_coor) ; ijarg = ijarg + 1
     CASE ('-d')
        CALL getarg(ijarg, c_dirout) ; ijarg = ijarg + 1
        ! create directory if necessary
        CALL SYSTEM ("mkdir -p "//TRIM(c_dirout) )
     CASE ('-r')
        lg_rename     = .TRUE.
     CASE ('-nc3')
        lg_nc3        = .TRUE.
     CASE ('-v')
        lg_verbose    = .TRUE.
     CASE ('-w')
        lg_win        = .TRUE. 
        CALL getarg(ijarg, cdum) ; ijarg = ijarg + 1 ; READ(cdum,*) mini
        CALL getarg(ijarg, cdum) ; ijarg = ijarg + 1 ; READ(cdum,*) maxi
        CALL getarg(ijarg, cdum) ; ijarg = ijarg + 1 ; READ(cdum,*) minj
        CALL getarg(ijarg, cdum) ; ijarg = ijarg + 1 ; READ(cdum,*) maxj
     CASE ('-kmax') 
        lg_kmax       = .TRUE. 
        CALL getarg(ijarg, cdum) ; ijarg = ijarg + 1 ; READ(cdum,*) maxk
     CASE ('-agrif') 
        lg_agrif     = .TRUE.
     CASE ('-b')
        CALL getarg(ijarg, cdum) ; ijarg = ijarg + 1 ; READ(cdum,*) nbloc
     CASE ('-encoding')
        CALL getarg(ijarg, cdum) ; ijarg = ijarg + 1 ; READ(cdum,*) ndigit
        WRITE(cl_format,'(a,i1,a,i1,a)') '("_",I',ndigit,'.',ndigit,',".nc")'
        WRITE(c_pattern, cl_format) 0
     CASE DEFAULT
        PRINT *, TRIM(cdum),' : Unknown option '
        STOP
     END SELECT
  ENDDO

  ncurrent = 1  ; nstop = 0 ; ndone = 0
  IF ( mpirank == 0) THEN  ! master
     DO WHILE ( nstop /= mpisize -1 )   ! loop until all mpi task exited ( except task 0)
        CALL MPI_RECV (iask, 1, MPI_INTEGER , MPI_ANY_SOURCE , MPI_ANY_TAG , MPI_COMM_WORLD ,mpistatus, mpierr)
        IF (iask /= istopflag ) THEN ! return current number
           iask=ncurrent
           CALL MPI_BSEND (iask,1, MPI_INTEGER ,mpistatus( MPI_SOURCE ), mpistatus( MPI_TAG ), MPI_COMM_WORLD ,mpierr)
           ncurrent=ncurrent+1
        ELSE
           nstop=nstop+1
           IF (lg_verbose) PRINT '(i4.4,1x,a, i4.4)', mpirank,' receive stop from task ', mpistatus( MPI_SOURCE )
        ENDIF
     ENDDO
     PRINT *, '### 0 ###  done '
  ELSE
     DO WHILE ( ncurrent <   nfile  )
        ! Current mpirank ask for a number to rank #0
        CALL MPI_BSEND (iask_number, 1, MPI_INTEGER ,0 , 1           , MPI_COMM_WORLD ,                    mpierr)
        CALL MPI_RECV  (ncurrent,    1, MPI_INTEGER ,0 , MPI_ANY_TAG , MPI_COMM_WORLD , MPI_STATUS_IGNORE, mpierr)
        IF ( ncurrent <= nfile ) THEN
           ndone = ndone + 1 
           nndone = ndone  ! nndone is a variable declared in modncfile even in case of mono proc program
           IF (lg_verbose) PRINT '(I4.4, I4, x,a)', mpirank, ndone,  ' in charge of '// TRIM(cf_list( ncurrent))
           ! perform  action ncurrent 
           CALL RebuildNc ( ncurrent )
        ENDIF
     ENDDO
     PRINT *, '+++ ',mpirank,' +++ done'
     ! when mpi task finished, send istopflag to master task
     CALL MPI_BSEND (istopflag, 1, MPI_INTEGER ,0 ,  mpirank ,        MPI_COMM_WORLD ,mpierr)
  END IF

  CALL MPI_FINALIZE (mpierr)

END PROGRAM mergefile_mpp4
