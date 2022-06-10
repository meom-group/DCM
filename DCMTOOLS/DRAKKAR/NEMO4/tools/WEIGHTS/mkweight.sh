#!/bin/bash
  #  mkweight.sh is a wrapper to nocs-weight generator.
  ##-------------------------------------------------------------
  ##  $Rev: 41 $
  ##  $Date: 2010-12-10 12:12:38 +0100 (Fri, 10 Dec 2010) $
  ##  $Id: mkweight.sh 41 2010-12-10 11:12:38Z molines $
  ##--------------------------------------------------------------

# this function indicate how to use mkweight.sh
usage() {
      echo
      echo
      echo
      echo USAGE : $(basename $0 ) [-h] [-d] [-x] -c coordinates_file -M land_sea_mask  -m method 
      echo "                 [-h]       : this help"
      echo "                 [-x]       : only compile the code"
      echo "                 [-d]       : debug mode : do not erase intermediate files"
      echo "                              Also usefull if you want to use scripinterp"
      echo "        -c coordinates_file : indicate the coordinate file name of"
      echo "                              the target model"
      echo "        -M Land_sea_mask    : indicate the land sea mask file name of "
      echo "                              the atmospheric grid (input) "
      echo "        -m method           : The method of interpolation : "
      echo "                              bilinear  or  bicubic "
      exit 0
         }

sep()    {
      echo "########################################################################"
         }

# this function will compile the code from the source. You need to have envrironment variable NCDF set
compile() { 
          cd ../
          ./maketools -m $MACHINE -n WEIGHTS
          cd WEIGHTS
          }
###############
BINDIR=BLD/bin

if [ $# = 0 ] ; then usage ; fi          #  stop if no argument

while getopts :hdxc:M:m: opt ; do         # parse command line
  case $opt in 
    (h) usage ;;
    (d) debug=1 ;;
    (x) mkexe=1 ;;
    (c) COORD=$OPTARG ;;
    (M) LSM=$OPTARG   ;;
    (m) METHOD=$OPTARG ; if [[ $METHOD != 'bilinear' && $METHOD != 'bicubic' ]] ; then 
       sep
       echo " ERROR : the method is either 'bilinear' or 'bicubic' but not '"$METHOD"'"
       sep 
       usage ; fi  ;;
    (:)  echo $( basename $0 )": -"${OPTARG}" option : missing value" ; sep ; usage ;;
    (\?) echo $( basename $0 ): option -$OPTARG not valid. ; sep ; usage ;;
  esac 
done

if [ $mkexe ] ; then   # compile only
   compile 
   exit 0
fi

if [[ ! $COORD || ! $LSM  || ! $METHOD ]] ; then   # stop if missing informations to proceed
   sep
   echo "You must specify the 3 options -c ..., -M ... and -m ... !"
   sep
   usage
fi

if [ ! -f $COORD ] ; then                         # stop if files are missing
   echo " Coordinates file " $COORD " not found "
   sep
   exit 1
fi

if [ ! -f $LSM ] ; then
   echo " Land sea mask file  " $LSM " not found "
   sep
   exit 1
fi

if [[ ! -x $BINDIR/scripgrid.exe || ! -x $BINDIR/scrip.exe || ! -x $BINDIR/scripshape.exe ]] ; then   # recompile the code if necessary
   echo " Wait a moment, the code is to be recompiled ..."
   if [ $debug ] ; then 
     compile
   else
     compile > /dev/null 2>&1 
   fi
   echo " Now proceed to weight generation. "
fi

echo COORD = $COORD
echo LSM   = $LSM
echo METHOD = $METHOD

fid=$(basename $COORD)
fid=${fid/.nc}  ; fid=${fid/_coordinates}  ;  fid=${fid/coordinates_} ;  fid=${fid/coordinates}

echo  weight_${METHOD}_$fid.nc

ln -sf $COORD coordinates.nc
ln -sf $LSM lsm.nc

cat namelist.skel | sed -e "s/<MASK>/lsm.nc/" -e "s/<METHOD>/$METHOD/" -e "s/<FID>/$fid/" > namelist


$BINDIR/scripgrid.exe namelist
$BINDIR/scrip.exe namelist
$BINDIR/scripshape.exe  namelist

if [ ! $debug ] ; then   # clean the intermediate files
  \rm -r namelist lsm.nc coordinates.nc remap_nemo_grid.nc 
  \rm -r remap_data_grid.nc data_nemo.nc
fi
