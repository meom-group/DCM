#!/bin/bash
# This script produce a map of the bathymetry as defined in the namelist, with the mpp sub domains
# superimposed. Land only domains are marked with a large cross.
#  This script is a post-processing of the files produced with mpp_optimize.exe -keep i j 


usage () {
     echo USAGE : $(basename $0) -f proc_file [-n namelist ] [ -h ]
     
     echo "     -f proc_file : specify the name of the text file with the sub-domain layout"
     echo "     [-n namelist ] : indicate the name of the namelist file if not namelist "
     echo "     [-h ] : print this help message "
     echo "     Example : $(basename $0) -f ORCA2-zahir-007x011_069 -n namelist_orca2 "
     exit
         }

LookInNamelist()    { 
     eval grep -e $1 $namelist  | tr -d \' | tr -d \"  | sed -e 's/=/  = /' | \
           awk ' {if ( $1 == str ) print $3 }' str=$1 
                    }

if [ $# = 0 ] ; then usage ; fi
namelist=namelist
CHART=chart

while  getopts :hn:f: V  ; do
   case $V in
     (h)  usage ;;
     (n) namelist=${OPTARG} ;;
     (f) overdata=${OPTARG} ;;
     (c) confcase=${OPTARG} ;;
     (:)  echo ${b_n}" : -"${OPTARG}" option : missing value" 1>&2;
        exit 2;;
     (\?) echo ${b_n}" : -"${OPTARG}" option : not supported" 1>&2;
        exit 2;;
   esac
done

bathy=$(LookInNamelist cn_fbathy )
bathyvar=$(LookInNamelist cn_var )

$CHART -clrdata  $bathy -pixel  -jetpal 128  -clrmark bathy.clrmark \
   -format PALETTE i4 -xstep 80 -ystep 85  -spval 0 -clrvar $bathyvar\
   -overdata $overdata -o $overdata.cgm -string 0.5 0.95 1.3 0 "$overdata"  -prev

echo plot done in $overdata.cgm . Display it with idt $overdata.cgm
