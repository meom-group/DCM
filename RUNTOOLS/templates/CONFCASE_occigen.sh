#!/bin/bash
#SBATCH --nodes=600
#SBATCH --ntasks=14400
#SBATCH --ntasks-per-node=24
#SBATCH --constraint=HSW24
##SBATCH -p all 
#SBATCH --threads-per-core=1
#SBATCH -J nemo_occigen
#SBATCH -e nemo_occigen.e%j
#SBATCH -o nemo_occigen.o%j
#SBATCH --time=24:00:00
#SBATCH --exclusive

set -x
ulimit -s 
ulimit -s unlimited

CONFIG=<CONFIG>
CASE=<CASE>

CONFCASE=${CONFIG}-${CASE}
CTL_DIR=$PDIR/RUN_${CONFIG}/${CONFCASE}/CTL

# Following numbers must be consistant with the header of this job
export NB_NPROC=3150    # number of cores used for NEMO
export NB_NPROC_IOS=34  # number of cores used for xios (number of xios_server.exe)
export NB_NCORE_DP=0    # activate depopulated core computation for XIOS. If not 0, RUN_DP is
                        # the number of cores used by XIOS on each exclusive node.
# Rebuild process 
export MERGE=0          # 0 = on the fly rebuild, 1 = dedicated job
export NB_NPROC_MER=231 # number of cores used for rebuild on the fly  (1/node is a good choice)
export NB_NNODE_MER=20  # number of nodes used for rebuild in dedicated job (MERGE=0). One instance of rebuild per node will be used.
export WALL_CLK_MER=3:00:00   # wall clock time for batch rebuild
export CONSTRAI_MER=HSW24  # For occigen, define either HSW24 or BDW28

date
#
echo " Read corresponding include file on the HOMEWORK "
.  ${CTL_DIR}/includefile.sh

. $RUNTOOLS/lib/function_4_all.sh
. $RUNTOOLS/lib/function_4.sh
#  you can eventually include function redefinitions here (for testing purpose, for instance).
. $RUNTOOLS/lib/nemo4.sh
