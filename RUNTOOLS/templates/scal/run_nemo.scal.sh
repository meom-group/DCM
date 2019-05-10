#!/bin/bash
# run_nemo_scal.sh
######################################################
#  wrapper for submitting nemo4 scalability experiment
######################################################

usage() {
   echo
   echo "USAGE: $(basename $0) [-h] [-c CORES-per-node] jpni jpnj jpnij nxios"
   echo "  "
   echo "  PURPOSE:"
   echo "    Launch scalability experiment corresponding to domain decomposition given"
   echo "    in the arguments."
   echo "  "
   echo "  ARGUMENTS:"
   echo "    jpni : number of subdomains in the I-direction"
   echo "    jpnj : number of subdomains in the J-direction"
   echo "    jpnij : Total number of ocean only subdomains."
   echo "    nxios : Number of xios_server.exe to be launched."
   echo "  "
   echo "  OPTIONS:"
   echo "    -h : print this usage message."
   echo "    -c CORES-per-node : set number of cores per computing node. Default :" $ncpn
   echo "  "
   exit 0
        }

ncpn=28

if [ ! $PDIR ] ; then
   echo "ERROR : You must set up your environment for DCM before using the RUN_TOOLS."
   echo "        PDIR environment variable not set."
   usage
fi

while getopts :hc: opt ; do
   case $opt in
     (h) usage ;;
     (c) ncpn=${OPTARG} ;;
     (*) ;;
    esac
done

shift $(($OPTIND-1))

if [ $# != 4 ] ; then
   echo "ERROR: incorrect number of arguments."
   usage
fi

set -x
jpni=$1
jpnj=$2
jpnij=$3
nxios=$4

ntasks=$(( jpnij + nxios ))
nodes=$(( ntasks / ncpn ))
if [ $(( ntasks % ncpn )) != 0 ] ; then
  nodes=$(( nodes + 1 ))
fi


if [ $ncpn = 24 ] ; then
   constraint='HSW24'
else
   constraint='BDW28'
fi

# create includefile.sh from template :
cat ./includefile.sh.tmp | sed -e "s/<JPNI>/$jpni/g"   -e "s/<JPNJ>/$jpnj/g"   \
          -e "s/<JPNIJ>/$jpnij/g" -e "s/<NXIOS>/$nxios/g" -e "s/<NODES>/$nodes/g" \
          -e "s/<NTASKS>/$ntasks/g" -e "s/<NCPN>/$ncpn/g" > includefile.sh_${jpni}_${jpnj}_${jpnij}

. ./includefile.sh_${jpni}_${jpnj}_${jpnij}

# Create namelist from template:
cat ./namelist.${CONFIG_CASE}.tmp  | sed -e "s/<JPNI>/$jpni/g" -e "s/<JPNJ>/$jpnj/g" \
          -e "s/<JPNIJ>/$jpnij/g" > namelist.${CONFIG_CASE}_${jpni}_${jpnj}_${jpnij}

# Create submit script from template:
cat ./${CONFIG_CASE}_${MACHINE}.sh.tmp  | sed -e "s/<JPNI>/$jpni/g" -e "s/<JPNJ>/$jpnj/g" \
          -e "s/<JPNIJ>/$jpnij/g" -e "s/<NXIOS>/$nxios/g" -e "s/<NCPN>/$ncpn/g" -e "s/<CONSTRAINT>/$constraint/g" \
          -e "s/<NODES>/$nodes/g" -e "s/<NTASKS>/$ntasks/g" > ${SUBMIT_SCRIPT}


date
set -x
echo " submitting ${SUBMIT_SCRIPT} ..."

$SUBMIT  ./${SUBMIT_SCRIPT}
