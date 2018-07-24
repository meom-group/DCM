PROGRAM mergefile4
  !!======================================================================
  !!                     ***  PROGRAM  mergefile4 ***
  !!=====================================================================
  !!  ** Purpose : Merge a bunch of splitted files into a single
  !!               one. The merge function work with MPP nemo output.
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
  INTEGER(KIND=4) :: numlst=10              ! logical unit for the list file (-F)

  CHARACTER(LEN=80 ) :: cdum
  CHARACTER(LEN=80 ) :: cpattern='_0000.nc'

  !!----------------------------------------------------------------------
  narg=iargc()

  IF ( narg == 0 ) THEN
     PRINT *,' usage :  mergefile4 -F | -f <file*>_0000.nc  [-c coordinate_file ] ...'
     PRINT *,'          [ -d output_directory] [ -r ] [-nc3] [-v]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :' 
     PRINT *,'         This program recombine splitted netcd files into a single global file.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'        use either one of the follwing ( not both) :'
     PRINT *,'        -F : takes all possible _0000.nc files in current directory'
     PRINT *,'           This option has been introduced to overrid shell limitations '
     PRINT *,'           regarding the length of the arguments.'
     PRINT *,'       or :'
     PRINT *,'        -f <file*>_0000.nc : list of full name of rank 0 <file*> to be rebuild'
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
       CALL SYSTEM( "find . -maxdepth 1 -name ""*"// TRIM(cpattern) // "*"" -printf ""%f\n"" | sort -r  > .zmergfile_list.tmp")
       nfile=0
       OPEN(numlst, FILE='.zmergfile_list.tmp')
       DO 
       READ(numlst,*,END=999)
       nfile=nfile+1
       ENDDO
999    PRINT *, 'Files to process : ', nfile
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
     CASE DEFAULT
        PRINT *, TRIM(cdum),' : Unknown option '
        STOP
     END SELECT
  ENDDO

  ! rebuild merged file
  DO jf = 1, nfile
     CALL RebuildNC ( jf )
  ENDDO   ! file list

END PROGRAM mergefile4
