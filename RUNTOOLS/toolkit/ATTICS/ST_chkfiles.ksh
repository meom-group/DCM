#!/bin/bash
# This script report the number of files for a given config case on both /scratch and /data
# Contribution of Gildas Mainsant to the JADE_TOOLS
#
# $Id: ST_chkfiles.ksh 363 2011-12-19 13:15:01Z molines $

#************************************************
# input arguments :
#^^^^^^^^^^^^^^^^^^
CONFIG=$1
CASE=$2
CONFCASE=${CONFIG}-${CASE}

if [ $# = 2 ] ; then
 ybeg=1980
 yend=2004
elif [ $# = 3 ] ; then
 ybeg=$3
 yend=$3
elif [ $# = 4 ] ; then
 ybeg=$3
 yend=$4
else
   echo USAGE: ST_chkfiles.ksh  CONFIG CASE
   echo "      or"
   echo "      ST_chkfiles.ksh CONFIG CASE year-init year-end"
   echo "      or "
   echo "      ST_chkfile.ksh CONFIG CASE year "
   exit 0
fi
#************************************************

CURRENTDIR=$( pwd )
datas="/data/$USER/${CONFIG}/${CONFIG}-${CASE}-S"
scras="/scratch/$USER/${CONFIG}/${CONFIG}-${CASE}-S"

echo "****************************************"
echo "In ${CONFIG}-${CASE}-S there are : "
echo " "
echo "for year :  in /scratch    |  in /data"
echo "            ^^^^^^^^^^^       ^^^^^^^^"
echo "  "

for year in $( seq $ybeg $yend ) ; do
        cd $scras
        if [ -d $year ] ; then
          cd ${year}
          scrainfo="for ${year} : $( ls -1 | wc -l ) 'files'"
          cd $datas
          if [ -d $year ] ; then
            cd ${year}
            datainfo="$( ls -1 | wc -l ) 'files'"
          else
            datainfo="No ${year}/ directory"
            datainfo="${datainfo} '... What are you waiting for ! ;-)'"
          fi
        else
          scrainfo="No ${year}/ directory  "
          cd $datas
          if [ -d $year ] ; then
            cd ${year}
            datainfo="$( ls -1 | wc -l ) 'files'"
          else
            datainfo="No ${year}/ directory"
          fi
        fi 
        echo "$scrainfo     |  $datainfo"
done

cd $CURRENTDIR

scram="/scratch/$USER/${CONFIG}/${CONFIG}-${CASE}-MEAN/"
datam="/data/$USER/${CONFIG}/${CONFIG}-${CASE}-MEAN/"

echo "****************************************"
echo "In ${CONFIG}-${CASE}-MEAN there are : "
echo " "
echo "for year :  in /scratch    |  in /data"
echo "            ^^^^^^^^^^^       ^^^^^^^^"
echo "  "

for year in $( seq $ybeg $yend ) ; do
        cd $scram
        if [ -d $year ] ; then
          cd ${year}
          scrainfo="for ${year} : $( ls -1 | wc -l ) 'files'"
          cd $datam
          if [ -d $year ] ; then
            cd ${year}
            datainfo="$( ls -1 | wc -l ) 'files'"
          else
            datainfo="No ${year}/ directory"
            datainfo="${datainfo} '... What are you waiting for ! ;-)'"
          fi
        else
          scrainfo="No ${year}/ directory  "
          cd $datam
          if [ -d $year ] ; then
            cd ${year}
            datainfo="$( ls -1 | wc -l ) 'files'"
          else
            datainfo="No ${year}/ directory"
          fi
        fi
        echo "$scrainfo     |  $datainfo"
done

cd $CURRENTDIR
echo ''



