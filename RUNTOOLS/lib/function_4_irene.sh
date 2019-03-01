#!/bin/bash

# $Id: function_4_irene 301 2011-08-26 10:03:13Z molines $
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

# rapatrie_res is normaly used for restart. On irene it is just the same as rapatrie
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
      echo expatrie_res does nothing on irene
               }
# --- 

# remove old restart files from the working directory
# clean_res() { \rm  $P_R_DIR/*.$1.tar.* ; }
 clean_res() {
      echo no restart cleaning on irene 
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
submit() {  
      cd $P_CTL_DIR 
      ccc_msub $1 > $TMPDIR/logsubmit 
      cd $TMPDIR  
         }
# ---

# function for running OPA : it takes the number of procs and name of program as argument
runcode() {
         ccc_mprun -n $*
#         mpirun -np $*
          }

# function for running OPA : it takes the number of procs and name of program as argument
runcode_u() {
#        ccc_mprun -n $*
         mpirun --bynode -n  $*
          }
# ---
# function for running OPA and XIOS : it takes the number of procs and name of programs as argument
#    runcode_mpmd nproc1 prog1 nproc2 prog2  [...] nprocn progn
#
#
runcode_mpmd() { 
#         mpirun -bynode  -np $3 $4 : -np $1 $2
         rm -f ./zapp.conf
         narg=$#
         npair=$(( narg / 2 ))
         if [ $npair = 0 ] ; then  # assume the argumet is already a zapp.conf file
          cp $1 ./zapp.conf
         else
           for ipair in $(seq 1 $npair) ; do
              echo $1 $2 >> zapp.conf
              shift 2
           done
         fi
         ccc_mprun -E '-m cyclic' -f zapp.conf
                  }

# ---
# function for running OPA and XIOS : it takes the number of procs and name of programs as argument
#    This function runs XIOS on depopulated nodes ( 1 XIOS/ node) 
#    runcode_mpmd_dp -cpn 16 -dp 2 nproc1 prog1 nproc2 prog2  [...] nprocn progn
# in nemo_4.sh : runcode_mpmd_dp -cpn 16 -dp $RUN_DP $NB_NPROC ./opa $NB_NPROC_IOS ./xios_server.exe

#
#-------------------- !!!!! runcode_mpmd_dp HAS TO BE TESTED !!!!!!!!!!!! ------------------------
#
runcode_mpmd_dp() {
      #hosts list   : build file and set list_node array
      rm -f ${CTL_DIR}/hosts
      list_node=( $(nodeset -e ${SLURM_NODELIST}  ) )
      set -x
      for node in ${list_node[@]} ;do  
         echo $node >> ${CTL_DIR}/hosts 
      done
      set +x

      # browse arguments of the function
      narg=$#
      iarg=1 ; ipair=1
      core_per_node=${SLURM_CPUS_ON_NODE} # can me overrid using -cpn option

      while [ $iarg -le $narg ] ; do
        case $1 in
        ( -cpn ) shift 1 ; core_per_node=$1       ; shift 1 ; iarg=$(( iarg + 2 )) ;;
        ( -dp  ) shift 1 ; xios_core_per_node=$1  ; shift 1 ; iarg=$(( iarg + 2 )) ;;
        ( *    )  
            case  $ipair  in 
              ( 1 ) nb_nproc=$1      ; shift 1 ; iarg=$(( iarg + 1 )) 
                    nemo_prog=$1     ; shift 1 ; iarg=$(( iarg + 1 ))
                    ipair=$(( ipair + 1 )) ;;
              ( 2 ) nb_nproc_xios=$1 ; shift 1 ; iarg=$(( iarg + 1 )) 
                    xios_prog=$1     ; shift 1 ; iarg=$(( iarg + 1 ))
                    ipair=$(( ipair + 1 )) ;;
            esac
        esac
      done

      nnodes_nemo=$(( nb_nproc          / core_per_node      ))
      nnodes_xios=$(( nb_nproc_xios     / xios_core_per_node ))
#      nnodes_tot=$(( nnodes_nemo        + nnodes_xios       ))   # not used ?
#      ncores_per_slot=$(( core_per_node / xios_core_per_node ))  # not used

      #--- rankfile.txt
      rm -f ${CTL_DIR}/rankfile.txt
      touch  ${CTL_DIR}/rankfile.txt
      rank=0

      set +x
      #--- NEMO nodes ---#
      #------------------#
      for i in `seq 0 $(( core_per_node -1 ))` ; do
         for inode in $(seq 0 $((nnodes_nemo - 1)) ) ; do
            node=${list_node[$inode]}
            echo "rank $rank=$node slot=$i" >> ${CTL_DIR}/rankfile.txt
            rank=$(($rank+1))
         done
      done  

      #--- XIOS nodes ---#
      #------------------#
      for inode in $(seq $(( $nnodes_nemo  )) $(( $nnodes_nemo + $nnodes_xios -1 ))) ; do
          node=${list_node[$inode]}
          if   [ $xios_core_per_node = 1 ] ; then    # Depopulated xios nodes with 1 tasks
             echo "rank $rank=$node slot=0:0" >> ${CTL_DIR}/rankfile.txt  # task on core 0 of socket 0
             rank=$(($rank+1))                                                    
          elif [ $xios_core_per_node = 2 ] ; then    # Depopulated xios nodes with 2 tasks
             echo "rank $rank=$node slot=0:0" >> ${CTL_DIR}/rankfile.txt  # 1st task on core 0 of socket 0
             rank=$(($rank+1))                                                    
             echo "rank $rank=$node slot=1:0" >> ${CTL_DIR}/rankfile.txt  # 2nd task on core 0 of socket 1
             rank=$(($rank+1))                                                     
          elif [ $xios_core_per_node = 3 ] ; then    # Depopulated xios nodes with 3 tasks
             echo "rank $rank=$node slot=0:0" >> ${CTL_DIR}/rankfile.txt  # 1st task on core 0 of socket 0
             rank=$(($rank+1))                                                    
             echo "rank $rank=$node slot=1:0" >> ${CTL_DIR}/rankfile.txt  # 2nd task on core 0 of socket 1
             rank=$(($rank+1))                                                     
             echo "rank $rank=$node slot=1:4" >> ${CTL_DIR}/rankfile.txt  # 3nd task on core 4 of socket 1 (arbitrary)
             rank=$(($rank+1))                                                     
          else                                  
             echo "Warning : xios_core_per_node = $xios_core_per_node, this case of xios depopulated nodes has not been written yet."
             exit
          fi  ; 
      done
      set -x

         #--- Submit Nemo on full nodes & XIOS on depopulated nodes ---#
         #-------------------------------------------------------------#
         mpirun -hostfile ${CTL_DIR}/hosts -rankfile ${CTL_DIR}/rankfile.txt \
                 -np $nb_nproc $nemo_prog : -np $nb_nproc_xios $xios_prog
         #mpirun -hostfile ${CTL_DIR}/hosts -rankfile ${CTL_DIR}/rankfile.txt -np $5 $6 : -np $7 $8 -nw
               }

# ---
# function for running OPA with scalasca : it takes the number of procs and name of program as argument
runcodescalasca() {
         scan -t mpiexec_mpt -n $1 $2
          }
# ---
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
      echo nothing to do on irene ! > /dev/null
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
   # initialization of variables on irene
   name=''
   account=''
   wallclock=01:00:00
   nodes=1
   cores=1
   jobtype='serial'
   cluster='nhm'
   queue='test'
   option=''
   mk_batch_hdr_core $@     # pass all input argument to the core of the function (in function_all)

# on irene wall clock must be passed in seconds ( need to translate argument given as hh:mm:ss )
wallclock_second=$( echo $wallclock | awk -F: '{print $1*3600 +$2*60 +$3}')

#  Build header for irene 
cat << eof 
#!/bin/bash
#MSUB -r $name
#MSUB -n $cores
eof
if [ $nodes != 1 ] ; then

cat << eof
#MSUB -N $nodes
eof

fi
cat << eof
#MSUB -T $wallclock_second
#MSUB -q $queue
#MSUB -o $name.o%I
#MSUB -e $name.e%I
#MSUB -A $account
eof
# add option if any
if [ $option ] ; then
cat << eof
#MSUB -E "$option"
eof
fi
               }
######################################################################################
