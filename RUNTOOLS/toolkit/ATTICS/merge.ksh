#!/bin/ksh
#MSUB -r merg_off_line             # Request name                     
#MSUB -n  50                      # Number of tasks to run            
#MSUB -N  10                       # Number of nodes to use          
#MSUB -T 36000                     # Elapsed time limit in seconds   
#MSUB -o zmergol.o%I            # Standard output. %I is the job id
#MSUB -e zmergol.e%I            # Error output. %I is the job id   
#MSUB -q standard                                                     
##MSUB -A gen0727                   # Project ID       
#MSUB -A ra2531                   # Project ID  
#MSUB -E "--exclusive"

CONFIG=ORCA025.L75
CASE=OCCITENS
ENSEMBLE_START=1
ENSEMBLE_END=50
ext=15
MERGE_EXEC=$WORKDIR/bin/mergefile_mpp4.exe
RUNTOOLS=$WORKDIR/DEV/RUN_TOOLS
NB_NPROC_MER=50

CONFCASE=${CONFIG}-${CASE}
DDIR=${DDIR:-$CDIR}
TMPDIR=$DDIR/TMPDIR_${CONFCASE}
XIOS=$DDIR/${CONFCASE}-XIOS.$ext
WKDIR=$XIOS/WRK.$$


. $RUNTOOLS/function_3.2.ksh
. $RUNTOOLS/function_3.2_all.ksh


mkdir -p $WKDIR
cd $WKDIR
mergeprog=$(basename $MERGE_EXEC )
# link all files in all member in a single dir
for member in $(seq  $ENSEMBLE_START $ENSEMBLE_END ) ; do
      nnn=$(getmember_extension $member nodot)
      mmm=.$nnn
      ln -sf $XIOS/$nnn/${CONFCASE}*.nc ./
done

ln -sf  $XIOS/coordinates.nc ./
ln -sf $MERGE_EXEC ./
getlst0000 24000
idxmax=$(( ${#lst0000[@]} - 1 )) #  index max in the list, starting from 0

for idx in $( seq 0 $idxmax ) ; do
     runcode_u $NB_NPROC_MER ./$mergeprog -f ${lst0000[$idx]} -c coordinates.nc -r
done
rm *.nc
