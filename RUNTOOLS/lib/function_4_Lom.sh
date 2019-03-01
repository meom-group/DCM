#!/bin/bash

# $Id: function_4_Lom.sh $
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
# runcodescalasca()
# lsrestart()  
# mk_batch_hdr
#-------------------------------------------------------
# Source machine independant functions (already done in the script)
# . ./function_4_all.sh
# MACHINE DEPENDANT functions
# copy : a wrap up for scp or cp 
copy()    { cp $1 $2 ; }
rcopy()   { cp $1 $2 ; }
copyfor() { ln -sf $1 $2 ;}
# ---

# rapatrie is a shell function which tries to copy $1 on the local dir with the name $4
#         if the name is not defined or is none, it just skip
#         it then first tries to find the file on disk directory $2
#         if it fails it tries to mfget from $3 directory 
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
#---

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

# rapatrie_res is normaly used for restart. On curie it is just the same as rapatrie
#  
rapatrie_res() { 
      rapatrie  $1 $2 $3 $4
               }
# ---

# expatrie  file $1 on the directory $2, with filename $3
expatrie() {
      copy $1 $2/$3
           }
# --- 

# expatrie_res  file $1 on the directory $2, with filename $3; copy $1 on local disk $4
expatrie_res() {
      echo expatrie_res does nothing on occigen2
               }
# --- 

# remove old restart files from the working directory
# clean_res() { \rm  $P_R_DIR/*.$1.tar.* ; }
 clean_res() {
      echo no restart cleaning on occigen2
             }
# --- 

# check existence and eventually create directory
chkdir() { if [ ! -d $1 ] ; then mkdir $1 ; fi ; }
# --- 

# chkdirg  : Usage: chkdirg DATA space_directory
#    check the existence of a directory on DATA space . Create it if not present
chkdirg() { if [ ! -d $SDIR/$1 ] ; then mkdir $SDIR/$1 ; fi  ; }
# --- 

# chkfile  : Usage: chkfile   dir_of_file/file  #CINES ajout 22-07-08
#    check the existence of a directory on dir_of_file on login_node .
chkfile() { if [ ! -f $1 ] ; then exit 1 ; fi  ; }
# --- 

# mkordre
mkordre() { cd $SDIR/${CONFIG}/${CONFIG_CASE}-S/ ; ~/bin/mkordre  ; }
# ---

# function for submitting jobs; modified for JADE

submit() { cd ${P_CTL_DIR} 
           if [ -f ~/.bad_node ] ; then 
           sbatch -x $(cat ~/.bad_node) $1 > $TMPDIR/logsubmit 
           else
           sbatch $1 > $TMPDIR/logsubmit 
           fi
           cd $TMPDIR 
         }
# ---

# function for running OPA : it takes the number of procs and name of program as argument
runcode() {
          srun  --mpi=openmpi  -m cyclic  -n $*
#         ccc_mprun -n $*
#         mpirun -np $*
          }

# function for running OPA : it takes the number of procs and name of program as argument
# unpopulated 
runcode_u() {
         mpirun --bynode -np $*
          }
# ---
# function for running OPA ans XIOS : it takes the number of procs and name of programs as argument
#    runcode_mpmd  nproc1 prog1    nproc2 prog2   ... nprocn progn
runcode_mpmd() {

# build a task file in the local directory, according to input parameters.
# in the main script, prog1 is nemo, prog2 is xios. Using -m cyclic, force a bind by node, so it is
# more clever to put xios first in the file in order to place XIOS on the first socket of various nodes.
         rm -f ./ztask_file.conf
         narg=$#
         npair=$(( narg / 2 ))
         if [ $npair = 0 ] ; then  # assume the argumet is already a zapp.conf file
          cp $1 ./ztask_file.conf
         else

n1=$1
n2=$3
prog1=$2
prog2=$4

rm -f ./ztask_file.conf
echo 0-$(( n2 - 1 ))         " " $prog2 > ./ztask_file.conf
echo ${n2}-$(( n1 + n2 -1 )) " " $prog1 >> ./ztask_file.conf
        fi

srun --mpi=openmpi -m cyclic --multi-prog ./ztask_file.conf
#srun --mpi=pmi2  -m cyclic -K1 \
#    --cpu_bind=map_cpu:0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27\
#    --multi-prog  ./ztask_file.conf
               }
# ---
         
# ---
# function for running OPA ans XIOS : it takes the number of procs and name of programs as argument
#    runcode_mpmdOLD  nproc1 prog1    nproc2 prog2   ... nprocn progn
runcode_mpmdOLD() {
         narg=$#
         npair=$(( narg / 2 ))
         nrest=$(( narg % 2 ))

         if [ $nrest != 0 ] ; then
           echo " You should give an even number of arguments : pair \( nproc  exec \)"
           return 1
         fi

         rm -f ./zapp.conf
         for ipair in $(seq 1 $npair) ; do
            echo $1 $2 >> zapp.conf
            shift 2
         done
#for zn in $(seq 1 56) ; do
#      cat << eof >> zapp.conf
#      7  $2
#      1   $4
#eof
#done


       ccc_mprun -f zapp.conf
               }
# ---


# function for running OPA with scalasca : it takes the number of procs and name of program as argument
runcodescalasca() {
         scan -t mpiexec_mpt -n $1 $2
          }

# function that source the .bashrc file if necessary (depends on the batch scheduler)
srcbash() {  
      source $HOME/.bashrc
          }  
# ---

# function save_nc : dummy function on jade
save_nc () {
     echo nothing to do ! > /dev/null
           }
# ---

# save_arch_file   file  archive_dir  # used when XIOS in use
save_arch_file () {
      echo nothing to do on curie ! > /dev/null
                  }
# --

#--------------------------------------------------------------------------------------
# Restart file management
# function to perform a ls on the remote restart directory
lsrestart() { 
       ls $F_R_DIR/ | grep $1
            }
# ---

#-----------------------------------------------------------------------------------------
# Make batch header for submitted scripts
# mk_batch_hdr  --name name --wallclock wallclock --account account --nodes nodes --cores cores --par --seq --option "options line" --help
mk_batch_hdr() {
   # initialization of variables on occigen2
   name=''
   account=''
   wallclock=01:00:00
   nodes=1
   cores=1
   jobtype='serial'
   cluster='nhm'
   queue='test'
   constraint=BDW28
   mk_batch_hdr_core $@     # pass all input argument to the core of the function (in function_all)

# on curie wall clock must be passed in seconds ( need to translate argument given as hh:mm:ss )
wallclock_second=$( echo $wallclock | awk -F: '{print $1*3600 +$2*60 +$3}')

#  Build header for Lomonosov

cat << eof 
#!/bin/bash
#SBATCH -J $name
#SBATCH -p gpu
#SBATCH --nodes=$nodes
#SBATCH --ntasks=$cores
#SBATCH --ntasks-per-node=8
#SBATCH --time=$wallclock
#SBATCH -e $name.e%j
#SBATCH -o $name.o%j
#SBATCH --constraint=$constraint
#SBATCH --exclusive
eof
               }
######################################################################################
