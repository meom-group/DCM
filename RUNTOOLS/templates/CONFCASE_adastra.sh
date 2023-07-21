#!/bin/bash
#SBATCH --account=gda2307
#SBATCH --job-name=JOB_
#SBATCH --constraint=GENOA
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=192 # --cpus-per-task=1
#SBATCH --threads-per-core=1 # --hint=nomultithread
#SBATCH --exclusive
#SBATCH --output=z%x.o%j
#SBATCH --error=z%x.e%j
#SBATCH --time=02:00:00

set -x
ulimit -s 
ulimit -s unlimited

CONFIG=<CONFIG>
CASE=<CASE>

CONFCASE=${CONFIG}-${CASE}
CTL_DIR=$PDIR/RUN_${CONFIG}/${CONFCASE}/CTL

##################################################################
#  WARNING: On Jean-Zay, It seems that having NEMO and XIOS_SERVER
#        running on the same node, lead to freezing in iom_init
#        A workaround is to put all the xios on a separate node 
#        using NB_NCORE_DP /= 0 
#################################################################
# Following numbers must be consistant with the header of this job
export NB_NPROC=39     # number of cores used for NEMO
export NB_NPROC_IOS=1  # number of cores used for xios (number of xios_server.exe)
export NB_NCORE_DP=0   # activate depopulated core computation for XIOS. If not 0, RUN_DP is
                       # the number of cores used by XIOS on each exclusive node.
# Rebuild process 
export MERGE=0         # 1 = on the fly rebuild, 0 = dedicated job
export NB_NPROC_MER=15 # number of cores used for rebuild on the fly  (1/node is a good choice)
export NB_NNODE_MER=1  # number of nodes used for rebuild in dedicated job (MERGE=0). One instance of rebuild per node will be used.
export WALL_CLK_MER=3:00:00   # wall clock time for batch rebuild
export ACCOUNT=gda2307 # account to be used

date
#
echo " Read corresponding include file on the HOMEWORK "
.  ${CTL_DIR}/includefile.sh

. $RUNTOOLS/lib/function_4_all.sh
. $RUNTOOLS/lib/function_4.sh
#  you can eventually include function redefinitions here (for testing purpose, for instance).
. $RUNTOOLS/lib/nemo4.sh
