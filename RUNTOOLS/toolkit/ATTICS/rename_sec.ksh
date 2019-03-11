#!/bin/ksh
#MSUB -r renam_sec             # Request name                     
#MSUB -n  1                      # Number of tasks to run            
##MSUB -N  1                       # Number of nodes to use          
#MSUB -T 36000                     # Elapsed time limit in seconds   
#MSUB -o zrenamsec.o%I            # Standard output. %I is the job id
#MSUB -e zrenamsec.e%I            # Error output. %I is the job id   
#MSUB -q standard                                                     
##MSUB -A gen0727                   # Project ID       
#MSUB -A ra2531                   # Project ID  
##MSUB -E "--exclusive"


CONFIG=ORCA025.L75
CASE=OCCITENS

RUNTOOLS=$WORKDIR/DEV/RUN_TOOLS/

CONFIG_CASE=${CONFIG}-${CASE}
ENSEMBLE_START=1
ENSEMBLE_END=50
ext=15

TMPDIR=$DDIR/TMPDIR_${CONFIG_CASE}

  cd $TMPDIR

#. ./includefile.ksh
. $RUNTOOLS/function_3.2.ksh
. $RUNTOOLS/function_3.2_all.ksh

# look for zoom (in one_file mode). Assume domain_ref without 'grid' keyword
  if [ -f 06-secme.xml ] ; then
      cat 0*sec*xml > ziodef.xml
  else
      ln -sf iodef.xml ziodef.xml
  fi

  iodef=ziodef.xml
  lis_zoom_tmp=$( for f in $(  grep -i domain_ref $iodef | grep -vi grid_ ) ; do \
  cmd=$(echo $f | grep domain_ref ) ; \
     if [ $cmd ] ; then eval $cmd ; echo $domain_ref ;fi ;\
           done )
     lis_zoom=$( for zoom in $lis_zoom_tmp ; do \
           echo $zoom ; done | sort -u )
# check enabled zoom
  lis_zoom_enabled=''
  zoom_enable=''
  for zoom in $lis_zoom ; do
      grep -wB 1 $zoom $iodef | head -1 | grep -iq 'enabled=".true."'
      if [ $? = 0 ] ; then
            lis_zoom_enabled="$lis_zoom_enabled $zoom"
            zoom_enable=1
      fi
  done
  lis_zoomid=$( for f in $lis_zoom_enabled ; do n=$(( ${#f} - 1 )) ; echo ${f:0:$n} ; done | sort -u )

  if [ $zoom_enable ] ; then
       for member in $(seq  $ENSEMBLE_START $ENSEMBLE_END ) ; do
          echo MEMBER $member
          nnn=$(getmember_extension $member nodot)
          mmm=$(getmember_extension $member      )

          cd $DDIR/${CONFIG_CASE}-XIOS.$ext/$nnn
          if [ $nnn ] ; then
             ln -sf ../coordinates.nc ./
             ln -sf ../*xml ./
          fi
          for zoomid in $lis_zoomid ; do
             post_process_one_file $zoomid
          done
       done
       cd $TMPDIR
  fi
