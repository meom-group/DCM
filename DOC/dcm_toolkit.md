# DCM toolkit
## Forewords
  After many years of run production, many 'small' scripts or programs were 
developped in order to help the run production. Even if 'small' they appear 
to be quite usefull ! We collect them into DCM toolkit for sharing and easier 
maintenance. In the following, a brief description of the toolkit is given, 
with an intent of classification by category. All scripts belonging to DCM toolkit
have a name starting with dcmtk_ . They are located in RUNTOOLS/toolkit directory.
Before use, there are some f90 programs to be compiled and installed in toolkit/src/.
(A Makefile is provided and no particular library is needed). `make install` copy the
binaries in RUNTOOLS/bin which is in the PATH, when using DCM.

  All dcmtk tools have the `-h `option to display usage instructions. The document
that you are actually readind was created with `dcmtk_mk_doc -o HOMEDCM/../DOC/dcm_toolkit.md`
## File management tools
### [dcmtk_chkmissing](../RUNTOOLS/toolkit/dcmtk_chkmissing)
  ```

USAGE: chkmissing [-h] dir1 dir2

  PURPOSE:
     Try to infer missing files in dir2 from files found in dir1 and vice-versa

  ARGUMENTS:
     dir1 dir2 refer to year directories of -S or -MEAN data baser

  OPTIONS:
     -h : display this help message

  ```
### [dcmtk_chkunlim](../RUNTOOLS/toolkit/dcmtk_chkunlim)
  ```

USAGE : dcmtk_chkunlim [-h ] [-v] [-p pattern] 

  PURPOSE:
     Scan the nc files in the current directory, and check if UNLIMITED 
     dimension is different from 0.
     With no option, gives the name of corrupted files (UNLIMITED dim = 0 ) 
 
  OPTIONS :
     -h : print this usage message
     -v : print the UNLIMITED line for all files
     -p pattern : use pattern for filtering the files to scan 
 
  ```
### [dcmtk_chkxios](../RUNTOOLS/toolkit/dcmtk_chkxios)
  ```

USAGE : dcmtk_chkxios -c [ -h ] [ -m ]

  PURPOSE : 
      Check files produced when using XIOS. Check that the number of
      subdomain files agrees with the number indicated within the file.
      This script must be run into the directory where all splitted
      files are, in particular the _0000.nc files. 
      This script is supporting ensemble runs.
       
  ARGUMENTS :
      -c : launch check. (With no arguments the help message is printed).
       
  OPTIONS : 
      -h : print this help message
      -m : migrate faulty set of files to local SAVE directory

  ```
### [dcmtk_cleano](../RUNTOOLS/toolkit/dcmtk_cleano)
  ```

USAGE: dcmtk_cleano [-h]

  PURPOSE:
     The goal is to reduce the size of this job.o erasing useless 
     information. Restore the date of the original file.  MACHINE is 
     supposed to be in the environment.

  OPTIONS:
     -h : Display this help message.
 
  ```
### [dcmtk_cpdd](../RUNTOOLS/toolkit/dcmtk_cpdd)
  ```

USAGE: dcmtk_cpdd [-h] 

  PURPOSE:
     This script is used to copy files from either a -S, -MEAN or -R 
     directory, located on /data/molines/WORKDIR to /data/molines/WORKDIR/SDIR (mirroring).
     For -S and -MEAN, it must be used in a 'year' subdirectory. 
     The same directory will be copied on /data/molines/WORKDIR/SDIR. 
     If files already exist, they are just skipped.
     Note that depending on the archiving system, mirroring is not
     necessarily the best choice as it does'nt reduce the number of
     used inodes on /data/molines/WORKDIR/SDIR file system.
   
  ARGUMENTS:
     This script is used without any arguments, CONFIG and CASE are deduced
     from the directory name.  mcp script (multiple cp) uses this script 
     (and its symetrical ddcp).  Look at mcp -h for more details 

  OPTIONS:
      -h : Display this help page 

  ```
### [dcmtk_ddcp](../RUNTOOLS/toolkit/dcmtk_ddcp)
  ```

USAGE: dcmtk_ddcp [-h] 

  PURPOSE:
     This script is used to copy files from either a -S, -MEAN or -R 
     directory, located on /data/molines/WORKDIR/SDIR to /data/molines/WORKDIR (mirroring).
     For -S and -MEAN, it must be used in a 'year' subdirectory. 
     The same directory will be copied on /data/molines/WORKDIR. 
     If files already exist, they are just skipped.
     Note that depending on the archiving system, mirroring is not
     necessarily the best choice as it does'nt reduce the number of
     used inodes on /data/molines/WORKDIR/SDIR file system.
   
  ARGUMENTS:
     This script is used without any arguments, CONFIG and CASE are deduced
     from the directory name.  mcp script (multiple cp) uses this script 
     (and its symetrical cpdd).  Look at mcp -h for more details 

  OPTIONS:
      -h : Display this help page 

  ```
### [dcmtk_inode](../RUNTOOLS/toolkit/dcmtk_inode)
  ```

USAGE: dcmtk_inode [-h] [list of directories]

  PURPOSE:
     Count the number of inodes in all the directories below the actual
     working directory. (All the sub directory tree is take into account.)
     This tool is very usefull to analyse places in your tree where there
     are too many files :) that may cause quota overflow.

  OPTIONS:
     -h : print this message 
     list_of_directories: restrict the inodes counting to the directories
                          in the list.

  ```
### [dcmtk_mcp](../RUNTOOLS/toolkit/dcmtk_mcp)
  ```

USAGE: dcmtk_mcp [-h] [-f firstyear] [-l lastyear ]

  PURPOSE:
     This script is used to perform file transfert between SDIR and WORKDIR.
     It is clever enough to determine the direction of the tranfert
     It uses dcmtk_cpdd or dcmtk_ddcp tools in an iterative way.
     It is now rather obsolete as data on SDIR are no more a mirror of WORKDIR.
     dcmtk_mcp must be used in the CONFCASE-S or CONFCASE-MEAN directory

  OPTIONS:
      -h : this help page 
      -f firstyear : specify the first year to transfer
                     If not given take the first existing year
      -l lastyear : specify the last year to transfer
                     If not given take the last existing year
   
  ```
### [dcmtk_mcpdd](../RUNTOOLS/toolkit/dcmtk_mcpdd)
  ```

USAGE: dcmtk_mcpdd [-h] 

  PURPOSE:
     Perform a mirroring of all the subdirectories of the current directory
     on SDIR. This script must be used in the CONFCASE-S or CONFCASE-MEAN 
     directory. All subdirectories are visited and files copied to SDIR/CONFIG/...

  ARGUMENTS:
     No arguments for this script.

  OPTIONS:
     -h : Display this help page 

  ```
### [dcmtk_mddcp](../RUNTOOLS/toolkit/dcmtk_mddcp)
  ```

USAGE: dcmtk_mddcp [-h] 

  PURPOSE:
     Perform a mirroring of all the subdirectories of the current SDIR directory
     on WORKDIR. This script must be used in the SDIR's CONFCASE-S or CONFCASE-MEAN 
     directory. All subdirectories are visited and files copied to WORKDIR/CONFIG/...

  ARGUMENTS:
     No arguments for this script.

  OPTIONS:
     -h : Display this help page 

  ```
### [dcmtk_mkordre](../RUNTOOLS/toolkit/dcmtk_mkordre)
  ```

USAGE: dcmtk_mkordre [-h ] [-a] [ extra file types ]

  PURPOSE:
     When having a bunch of files in a single directory, try to organize them
     by year directories.
     This script is used without arguments in the CONFCASE-S directory or in 
     CONFCASE-R directory.
     For restart files in -R, they are moved to v2.xx subdirs.
     For output files in -S, it moves the files in their respective year directory
     following file types are moved :
     gridT gridU gridV gridW icemod ptrcT diadT trends dynT flxT gridLOM gridSST gridSSU gridSSV 
     You can specify list of extra file type as arguments to this script
     0 fine agrif grid specified

  ARGUMENTS:
     If necessary, give a list of extra file type to deal with.

  OPTIONS:
      -h : this help message 
      -a number of subgrid : in case of AGRIF runs, specify the number of subgrids 

  ```
### [dcmtk_mvmo2s](../RUNTOOLS/toolkit/dcmtk_mvmo2s)
  ```

USAGE : dcmtk_mvmo2s  segment_number [-h] 

  PURPOSE:
     This script rename the mooring files to the corresponding -S directory
     Must be used in a CTL directory 

  ARGUMENTS:
     segment_number : the number of the segment to migrate.

  OPTIONS:
     -h : print this help message

  ```
### [dcmtk_mvnc2s](../RUNTOOLS/toolkit/dcmtk_mvnc2s)
  ```

USAGE : dcmtk_mvnc2s [-h ]

  PURPOSE:
     This script is used to move rebuilt nc file to the CONFCASE-S directory
     It is used in <CONFCASE>-XIOS.<seg> output directory after a run.
     This script does not store the files on the archive but only move then
     to 
     You can use cpdd of mcp to perform the archiving (as a mirror -deprecated-)
     This script manage ensemble run, each member is considered as a particular 
     case, with extension .NNN

  OPTIONS :
     [ -h ] : show this help message 

  ```
### [dcmtk_superinode](../RUNTOOLS/toolkit/dcmtk_superinode)
  ```

USAGE: dcmtk_superinode [-h] [list of directories]

  PURPOSE:
     Give the number of files in the directories of the list (or all directories),
     the size of each directories (Ko), mean file size (Ko) as well as the 
     directory name.

  OPTIONS:
     -h : print this message 
     list of directories: restrict action to the listed directories.

  ```
## Miscellaneous tools
### [dcmtk_mk_doc](../RUNTOOLS/toolkit/dcmtk_mk_doc)
  ```

USAGE: dcmtk_mk_doc [-h] -o OUT-md_file 

   PURPOSE:
      Scan the dcmtk_* scripts and sort then regarding the class it belongs to.
      Then output the usage message and build a md file.

   ARGUMENTS:
      -o OUT-md_file: Specify the output markdown file.

   OPTIONS:
      -h : Display this help message.

  ```
## Run management tools
### [dcmtk_journal_cpu](../RUNTOOLS/toolkit/dcmtk_journal_cpu)
  ```

USAGE : dcmtk_journal_cpu -w Journal_wiki.txt -c cpu_used.txt [-h] [-o output file] [-m machine ]

  PURPOSE:
     Use accounting information hold in cpu_used.txt in order to finish the run's journal
     created by dcmtk_journal_make.

  ARGUMENTS: (mandatory)
     -w Journal_wiki.txt  :  the name of the text file created from the wiki page
     -c cpu_used.txt      :  the name of the acounting file on HPC computer (see dcmtk_getcpu) 

  OPTIONS:
     -h                   :  this help message 
                          :  It can be obtained at https://reser.cines.fr/ with detailed option
     -o output file       :  specify the output file name ( optional, default is Journal_cpu.txt )
                          :  This file is then to be paste into the wiki, in place of the original table
     -m machine           :  Specify the machine you are working on [occigen]. ada is also suported
                             On ada, cpu_used_ada.txt file can be obtained with idrjar -d 01/01-nn/mm

  ```
### [dcmtk_journal_make](../RUNTOOLS/toolkit/dcmtk_journal_make)
  ```

USAGE: dcmtk_journal_make [-h ] [-f ] [-c confcase ] [-n name ] [-u user] [-o fileout]

  PURPOSE:
     This script create a wiki table (ReSTructured text trac format) with the segments of
     run present in the local directory (job output). It uses the namelist in  order
     to retrieve important parameters of the run, to put in the table. Hence,
     this script must be run in CTL directory for regular usage. It prepares the column
     for the CPU usage. This column can be completed afterward with accounting information
     (see dcmtk_journal_cpu) 

  OPTIONS:
     -h : help 
     -f : force : do not check that you are in a CTL dir 
     -n name : use name instead of nemo_<MACHINE>. <MACHINE> is inferred from hostname
            for irene, ada, occigen 
     -c confcase : use confcase instead of CONFCASE deduced from CTL
                 : this is usefull with option -f 
     -u user [default is molines ]
     -o output file [ default is journal.wiki 

  ```
### [dcmtk_split_logfile](../RUNTOOLS/toolkit/dcmtk_split_logfile)
  ```

USAGE:  dcmtk_split_logfile [-h] nemo_error_output_file

  PURPOSE:
     Split job.e and job.a when having a loop on segments within a single job.
     Ends up with individual job files (o and e) for each segment.
     Keep track of the date of the file.

  ARGUMENTS:
     nemo_error_output_file : the std error output file of a multi segment run.

  OPTIONS:
     -h : Display this usage message.

  ```
## Run performance tools
### [dcmtk_chkrate](../RUNTOOLS/toolkit/dcmtk_chkrate)
  ```

USAGE: dcmtk_chkrate [-h ] [-n job.o rootname]

  PURPOSE:
      Scans all the nemo_occigen.o* files in the current directory and 
      computes the running speed of the model during the first 180 timesteps 

  ARGUMENTS:
      No arguments, scan default nemo_occigen.o* file

  OPTIONS: 
    -h                : Display this help message
    -n job.o_rootname : Use rootname instead of nemo_occigen for the name of
                        the job output file.
  ```
### [dcmtk_elapsed](../RUNTOOLS/toolkit/dcmtk_elapsed)
  ```

USAGE:  dcmtk_elapsed [-h] nemo_machine.o 

  PURPOSE:
      This script Will scan the job outputfile and use the timing information to infer 
      a prognostic about the elapsed time required for 1 year of run.
      This command must be used in CTL where both logfiles and namelists are.

  ARGUMENT:
      nemo_machine.o : A std output job file produce by DRAKKAR nemo

  OPTIONS:
      -h : Display this help message

  ```
### [dcmtk_rate](../RUNTOOLS/toolkit/dcmtk_rate)
  ```

USAGE: dcmtk_rate [-h ] [-s istep ] [ -b nblock ] [ -f filename ] [ -m ] 
       [ -g ] [-t graphtype] [ -v ] [ -l ]

  PURPOSE: 
     Display the integration rate (step/minutes) computed from filename 
     (drakkar nemo job standard output). The analysis is done by chunks 
     of time steps, that are defined with the options.

  OPTIONS:
     -h : Display this help message. 
     -s istep : number of time step for the chunk of data to analyse [ 450 ]
     -b nblock : number of chunk of data to analyse [ 70 ]
     -f filename : nemo_machine.oxxx file name to consider [ none ]
     -m : indicate mean and std of the computed performance rates
     -g : pipe results to graph 
     -t graph output :  either X gif png ...[ X ]
     -v : verbose extra informations (do not use with -g option)
     -y ymin ymax :fix y axis scale to be between ymin and ymax 
     -l : label title above verbose informations (do use with -v options)

  ```
## Run progression tools
### [dcmtk_chkstp](../RUNTOOLS/toolkit/dcmtk_chkstp)
  ```

USAGE: dcmtk_chkstp [-h] [-d delay ] 

  PURPOSE:
      Display the evolution of time step to check if the run progresses
      smoothly. If run stop or too slow, a RED warning is dipslayed.

  ARGUMENTS:
      No arguments. Default delay is 5 seconds.

  OPTIONS:
     -h : this help page 
     -d delay : specify a delay (s) between each check 
                default is  5  seconds.
                decimal numbers can be used (eg 0.5 for 1/2 seconds)

  ```
### [dcmtk_eta](../RUNTOOLS/toolkit/dcmtk_eta)
  ```

USAGE: dcmtk_eta  [ -h | -help] [-j jobid] [-n name] [-u user] [-q queue] [-f ]

  PURPOSE:
     This tools aims at producing the Estimated Time of Arrival (eta) for a running
     job. It also provides usefull (or wacky!) real time information.

  OPTIONS:
     -h or -help : print this message 
     -j jobid    : take jobid to look in the running table
     -s sessid   : take sessid  to look in the running table
                 : This is usefull if many instances of name are running
     -n name     : take name instead of nemo_machine
     -q queue    : take queue instead of nhm for queue name
     -u user     : specify login name of the owner of the job
     -f          : fast : do not display progress bar
   
  ```
### [dcmtk_show_progression](../RUNTOOLS/toolkit/dcmtk_show_progression)
  ```

USAGE : dcmtk_show_progression [-h] [-t plotdevice] [-p prefix] [-s suffix] [-o png file].. 
                       [-m ] 

  PURPOSE : 
     Plot the number of files vs the time in the current directory.
     This is usefull when data transfert is slow, for estimating
     the real rate of transmission of the files.

  METHOD : 
     The basis of the command is ls -ltr --full-time to get the precise
     date and time of creation/access of each files in the directory.

  OPTIONS : 
     -h  : show this usage message
     -t plotdevice : plotdevice can be X [D] or png. X correspond to the
                     actual display while png create a png file show.png [D]
                     or the name specified with -o option
     -p prefix  : indicate a prefix for filtering the files in the directory
               The default value for the prefix is [ ORCA ]
     -s suffix  : indicate a suffix for filtering the files in the directory
               The default value for the suffix is [ nc ]
     -o png_file : Use png_file name instead of default show.png when using
               the -t png option
     -m  : Time axis in minutes, instead of default hours.

  ```
### [dcmtk_timeevol](../RUNTOOLS/toolkit/dcmtk_timeevol)
  ```

USAGE: dcmtk_timeevol  [-h] nemo.o1 nemo.o2 ..... minutes_max step_max 

  PURPOSE:
     This tools plot the number of time step as a function of time. It uses
     the timing information produced by DRAKKAR nemogcm at each time step, and
     available in the standard output job file. As far as job.o is accessible
     during the execution of a run, this tools is usefull to monitor real time
     run performances. It often shows the problems that may occur !
     This tools can plot many experiments on the same frame, allowing for easy
     performance comparison.

  ARGUMENTS: 
     The arguments consist of a list of std nemo job output files, followed
     by the maximum number of minutes to plot (x axis) and the maximum time
     steps to consider (y axis). 

  HISTORY: this tool was formely called timeevol4.ksh (for oldies !)

  ```
### [dcmtk_tourne](../RUNTOOLS/toolkit/dcmtk_tourne)
  ```

USAGE: dcmtk_tourne [-h] 

   PURPOSE:
      Display a spinning star according to the progress of time.step
      This script must be run in the running directory, where time.step is.

   OPTIONS:
      -h : Display this help message

  ```
### [dcmtk_visu](../RUNTOOLS/toolkit/dcmtk_visu)
  ```

USAGE: dcmtk_visu  [-h]  ocean.output 

  PURPOSE:
     This script decrypt an ocean.output file for max velocities, min SSH
     and produce a graph of these quantities vs step performed.
     This is helpfull for monitoring model instability.
     This tool is designed for NEMO release < 4.0. It is kept in dcmtk for memory
     and must be updated soon for NEMO4.

  ARGUMENT:
     ocean.output: the name of a particular ocean.ouput to check. It can be the
     the file corresponding to a running job.

  OPTIONS:
     -h : display this help message.

  ```
## Scalability experiment tools
### [dcmtk_chkrate_scal](../RUNTOOLS/toolkit/dcmtk_chkrate_scal)
  ```

USAGE : dcmtk_chkrate_scal [-h] [-n rootname of job.o] go 

  PURPOSE:
     Scan all the job output in current directory produced by scalability 
     experiment (e.g. nemo_occigen_10_20_180.o1234) and compute the rate (stp/mn)
     for every job. In the present script, 3 windows of 100 step are used.

  ARGUMENTS :
     go or anything ! : without arguments, display the usage message.

  OPTION : 
     [-n rootname of job.o] : default is nemo_occigen.

  ```
### [dcmtk_scal_mkjournal](../RUNTOOLS/toolkit/dcmtk_scal_mkjournal)
  ```

USAGE: dcmtk_scal_mkjournal [-h ] [-f ] [-c confcase ] [-n name ]

   PURPOSE:
      This tool is a variant of dcmtk_journal_make, where additional
      information is saved, in order to build a relevant wiki page for
      DRAKKAR scalability experiment. It is normally uses in a CTL directory
      from which a scalability experiment was performed.
      The output of this tool is using ReSTructured text wiki formating, working
      with trac. In the near future, this tool will be extended to produce markdown
      file.

   OPTIONS:
      -h : Display this help message.
      -f : force execution without checking  that you are in a CTL dir 
      -n name : pass alternative rootname instead of nemo_occigen
      -c confcase : use confcase instead of CONFCASE deduced from CTL
                  : this is usefull with option -f 
      This script must be run in CTL directory

  ```
### [dcmtk_scal_mkrunnemo.py](../RUNTOOLS/toolkit/dcmtk_scal_mkrunnemo.py)
  ```
 
USAGE: dcmtk_scal_mkrunnemo.py -h  -n <script_name> -c  <cores>  -x <nxios_min>
 
  PURPOSE:
     This python script ends up with a bash script (run.sh) holding
     the ad-hoc command lines for submitting a series of simulations
     with different domain decomposition, hence total number of cores.
  
     It assumes that you have already in the local directory an elementary
     script launching NEMO and taking 4 online arguments: jpni, jpnj, jpnij, nxios
     e.g.  ./run_nemo_occigen_scal.sh  20 30 50 3 
    
     This script is called from dcmtk_scal_prep which preprocess the processor.layout
     file produced by MPP_PREP, writing temporary log files (one for a given domain
     decomposition) named log_<totalcore>
 
  OPTIONS:
     -h : Display this help message
     -n <script_name> : define nemo scalability script name. Default:  ./run_nemo_occigen_scal.sh
     -c <cores> : gives number of core per compute note. Default : 28
     -x <cores> : gives the minimum number of cores dedicated to
           xios_server.exe. Default :  5
 
  ```
### [dcmtk_scal_prep](../RUNTOOLS/toolkit/dcmtk_scal_prep)
  ```

USAGE : dcmtk_scal_prep [-h] [-l layout_file] [-h] [-m minproc] [-M maxproc] [-s proc_step]

  PURPOSE:
     This tool is part of the process for setting up a scalability experiment.
     It assumes that you have already run the MPP_PREP tool for defining all the
     file, narrowing the results to a specific range of cores. In addition, it will
     call a python script that produces a job file to be submitted for all the 
     decomposition to test. At this point, the python tools assumes that the target
     machine is occigen (SLURM batch system). Hence, the job file may require some
     editing (batch header part) to fit your HPC system.
     Due to the interaction with MPP_PREP it is likely a good choice to run this tool
     where you have the MPP_PREP results, in particular with the screen.sh script.

  OPTIONS:
     -h : Display this help message.
     -l layout_file : Pass the name of the layout file produced by MPP_PREP  procedure.
           Default is processor.layout
     -m minproc : Gives the minimum number of core you want to test. Default=50 
     -M maxproc : Gives the maximum number of core you want to test. Default=600
     -s proc_step : Gives the core step for scaning the core range.  Default=50

  ```
## XIOS browse tools
### [dcmtk_lookforzoom](../RUNTOOLS/toolkit/dcmtk_lookforzoom)
  ```

USAGE: dcmtk_lookforzoom [-h] [iodef.xml]

  PURPOSE:
     This script scan the iodef.xml file (or the one given in argument, look for 
     enabled zoom and then scan the domain_ref.xml file in order to retrieve the 
     position imin imax jmin jmax for the zoomed box  output cat be made as a set
     of ncks commands which can be used to extract boxes from the global files.
     It can also provide an overlay file, that can be used by chart in order to 
     show the geographical position of the boxes with chart :) 
     This version of the script works for XIOS-1 only (an update will come soon 
     for XIOS-2.

  OPTIONS:
     -h : Display this help message
     iodef.xml : the name of a particular iodef.xml file, if not iodef.xml

  ```
### [dcmtk_mkdomain](../RUNTOOLS/toolkit/dcmtk_mkdomain)
  ```

USAGE: dcmtk_mkdomain -d DOM [-c coordinates_file] [ -h ] xmin xmax ymin ymax

  PURPOSE:
     Print the xml lines do be inserted in domain_ref.xml file in order to create 
     a new domain for XIOS output.
     This version is for XIOS1 : deprecated and obsolete
     New version to be written for XIOS2 

  ARGUMENTS:
     -d DOM : indicate a domain name for the xml file.
     xmin xmax ymin ymax : position of the domain in geographical coordinates.
                This 4 variable MUST be the last on the line.

  OPTIONS:
     -h : Display this help message.
     -c coordinates_file: Use this coordinate file instead of default coordinates.nc.

  ```
