#!/bin/bash

# $Id: function_4_ada.sh 1192 2012-11-08 12:43:19Z molines $
######################################################################################
#    bash script functions used below:
######################################################################################
# copy()
# rcopy()
# copyfor()
# rapatrie()
# core_rapatrie()
# rapatrie_res()
# expatrie()
# expatrie_res()
# clean_res()
# chkdir()
# chkdirg()
# chkfile()
# mkordre()
# submit() 
# runcode()
# runcode_mpmd()
# runcodescalasca()
# lsrestart()
#-------------------------------------------------------
# Source machine independant functions (already done in the script)
# . ./function_4_all.sh
# MACHINE DEPENDANT functions
# copy : a wrap up for scp or cp 
copy()    { \cp $1 $2 ; }
rcopy()   { \cp $1 $2 ; }
copyfor() { ln -sf $1 $2 ;}
# ---

# rapatrie is a shell function which tries to copy $1 on the local dir with the name $4
#         if the name is not defined or is none, it just skip
#         it then first tries to find the file on disk directory $2
#         if it fails it tries to find it from $3 directory on the storage machine
#         if it is not found anywhere, stop the script and issue a message
#         login_node is set in includefile.sh
rapatrie() {
      core_rapatrie $1 $2 $3 $4
      if [ $AGRIF = 1 ] ; then
         for idx in ${agrif_pref[@]} ; do
           core_rapatrie ${idx}_$1 $2 $3 ${idx}_$4
         done
      fi
           }
# ---

# The core of the rapatrie function.: called by rapatrie which also deals with agrif
core_rapatrie() {
     echo ${1:-none} | grep -iv none && \
    { if [ -f $4 ] ; then      # file already there
        echo $4 found
      elif [ -f $2/$1 ] ; then #copy from local file system
        cp $2/$1 $4
      else                     # copy from remote file system
        rcopy $3/$1 $4  ||  { echo $1 not found anywhere &&  exit ; }
      fi  ; }
           }
# ---

# rapatrie_res is normaly used for restart. On ada it is just the same as rapatrie
#  
rapatrie_res() { 
      rapatrie  $1 $2 $3 $4
               }
# ---

# expatrie  file $1 on the directory $2, with filename $3
expatrie() {
#      copy $1 $2/$3
       cp $1 $2/$3
           }
# ---

# expatrie_res  file $1 on the directory $2, with filename $3; copy $1 on local disk $4
expatrie_res() {
      mfput  $1 $2/$3
               }
# --- 

# remove old restart files from the working directory
# clean_res() { \rm  $P_R_DIR/*.$1.tar.* ; }
 clean_res() {
      echo no restart cleaning on ada
             }
# --- 

# check existence and eventually create directory
chkdir() { if [ ! -d $1 ] ; then mkdir $1 ; fi ; }
# ---

# chkdirg  : Usage: chkdirg DATA space_directory
#    check the existence of a directory on DATA space . Create it if not present
chkdirg() { if [ ! -d $SDIR/$1 ] ; then mkdir $1 ; fi  ; }
# --- 

# chkfile  : Usage: chkfile   dir_of_file/file  #CINES ajout 22-07-08
#    check the existence of a directory on dir_of_file on login_node .
chkfile() { if [ ! -f $1 ] ; then exit 1 ; fi  ; }
# --- 

# mkordre
mkordre() { cd $SDIR/${CONFIG}/${CONFIG_CASE}-S/ ; ~/bin/mkordre  ; }
# ---

# function for submitting jobs; modified for JADE
#submit() { rsh $login_node "cd ${P_CTL_DIR}; pwd; llsubmit $1 > $TMPDIR/logsubmit " ;}
submit() { cd ${P_CTL_DIR};  llsubmit $1 > $TMPDIR/logsubmit ; cd $TMPDIR  ; }
# ---

# function that source the .bashrc file if necessary (depends on the batch scheduler)
srcbash() {  echo > /dev/null  ; }  # this funcion is useless on ada
# ---

# function for running OPA : it takes the number of procs and name of program as argument
runcode() {
          poe ${*:2}
          }
# ---
# function for running mpp pogram with umpopulated  core computation AND node binding
runcode_u() {
          poe ${*:2}
          }
# ---
# function for running OPA ans XIOS : it takes the number of procs and name of programs as argument
#    runcode_mpmd  nproc1  prog1 nproc2 prog2
runcode_mpmd_ibm() {
         rm -f ./zrunfile
         for n in $(seq 1 $1 ) ; do
          echo $2 >> ./zrunfile
         done
         for n in $(seq 1 $3 ) ; do
          echo $4 >> ./zrunfile
         done

          poe -pgmmodel mpmd -cmdfile ./zrunfile
          }
# ---
# function for running OPA ans XIOS : it takes the number of procs and name of programs as argument
#    runcode_mpmd  nproc1  prog1 nproc2 prog2
runcode_mpmd() {
     mpirun --map-by node -np $3 $4 : -np $1 $2
          }
# ---


# function for running OPA with scalasca : it takes the number of procs and name of program as argument
runcodescalasca() {
         scan -t mpiexec_mpt -n $1 $2
          }
# ---
# ---
# function save_nc : put files on ergon
save_nc () {
     CONFIG=$1
     CONFIG_CASE=$2
     shift 2
     for tag in $@ ; do
       year=${tag:1:4}
       rsh ergon " mkdir -p ${CONFIG}/${CONFIG_CASE}-S/$year "
       for f in ${pref}${CONFIG_CASE}_${tag}_*.nc ; do
          mfput $f ${CONFIG}/${CONFIG_CASE}-S/$year
       done
     done
           }
# ---

# save_arch_file   file  archive_dir  # used when XIOS in use
save_arch_file () {
      rsh ergon "mkdir -p $2"
      mfput $1 $2
                  }

# ---

#--------------------------------------------------------------------------------------
# Restart file management
# function to perform a ls on the remote restart directory
lsrestart() { 
       if [ -f $P_R_DIR/$1.1 ] ; then # restarts we are looking for are on P_R_DIR
          ls $P_R_DIR/ | grep $1
       else                           # look on F_R_DIR
          ls $SDIR/$F_R_DIR/ | grep $1
       fi
            }
# ---
#-----------------------------------------------------------------------------------------
# Make batch header for submitted scripts
# mk_batch_hdr  --name name --wallclock wallclock --account account --nodes nodes --cores cores --par --seq --option "options line" --help --memory
mk_batch_hdr() {
   # initialization for ada
   name='none'
   account='none'
   wallclock=01:00:00
   nodes=1
   cores=1
   jobtype='serial'
   cluster='none'
   queue='none'
   adapp=0
   memory='3.5gb'

   mk_batch_hdr_core $@


# ADA  batch header creation according to options

if [ $adapp = 1 ] ; then
  prepro='# @ requirements = (Feature == "prepost")'
else
  prepro=''
fi

if [ $jobtype = 'serial' ] ; then
cat << eof 
#!/bin/bash
# @ job_name = $name
# @ output = \$(job_name).\$(jobid)
# @ error  = \$(output)
# @ job_type = $jobtype
$prepro
# @ wall_clock_limit=$wallclock
# @ as_limit = $memory
# @ queue
eof
else
cat << eof 
#!/bin/bash
# @ job_name = $name
# @ output = \$(job_name).\$(jobid)
# @ error  = \$(output)
# @ job_type = $jobtype
$prepro
# @ wall_clock_limit=$wallclock
# @ total_tasks = $cores
# @ environment = MP_PROCS=\$(total_tasks)
# @ environment = NPROC=\$(total_tasks)
# @ as_limit = $memory
# @ queue
eof
fi

               }
######################################################################################
