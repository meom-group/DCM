#!/bin/bash

# functions
##  give usage of the script
usage() {
   echo "USAGE   :  dcmtk_inforun.sh"
   echo "PURPOSE : Give informations about the run corresponding to this CTL or CDF directory"
   exit 0
        }

## determine CONFCASE from the actual directory (CTL or CDF)
get_confcase() {
here=$(pwd)
b_name=$( basename $here )

case $b_name in
(CTL) tmp=$(dirname  $here )
      CONFCASE=$( basename $tmp )
      CONFIG=${CONFCASE%-*}
      CASE=${CONFCASE#*-} ;;

(CDF ) tmp=$(dirname $here  )
       tmp=$(dirname $tmp  )
      CONFCASE=$( basename $tmp )
      CONFIG=${CONFCASE%-*}
      CASE=${CONFCASE#*-} ;;
( * ) echo " cannot infer CONFASE !"
      echo " This is not a CTL dir nor a CTL/CDF dir "
      exit 1
esac
echo "     CONFIG   : " $CONFIG
echo "     CASE     : " $CASE
echo " ==> CONFCASE : " $CONFCASE
               }



# Main

get_confcase

# ensure to be in CTL dir
CTLDIR=$PDIR/RUN_${CONFIG}/${CONFCASE}/CTL
SSDIR=$DDIR/${CONFIG}/${CONFCASE}-S
STDIR=$SDIR/${CONFIG}/${CONFCASE}-S
MWDIR=$WORK/${CONFIG}/${CONFCASE}-MEAN
DIADIR=$WORK/${CONFIG}/${CONFCASE}-DIAGS/NC

DBFILE=${CONFCASE}.db
cd $CTLDIR

echo "      Last completed segment:  $( tail -2 $DBFILE | head -1 | awk '{print $1}' ) "
            tmp=$(tail -2 $DBFILE  | head -1 | awk '{print $NF}' )
echo "      Last date             : ${tmp:0:4}"-"${tmp:4:2}"-"${tmp:6:2}"
echo
echo " status on scratch :"
cd $SSDIR/1d
tmp=$( ls -1d *-concat | tail -1 )
echo "      Last concatenate year : " ${tmp%-concat}
echo
echo " status in WORK (-MEAN) "
cd $MWDIR/1y
tmp=$( ls -1d * | tail -1 )
echo "      Last year in -MEAN    : " ${tmp}
echo 
echo " Monitoring :"
cd $DIADIR
tmp=$(ls -1 *TSMEAN*nc | tail -1)
tmp=$( echo $tmp | sed -e "s/${CONFCASE}_y//" )
echo "      Last monitored year  : " ${tmp%%.*}

echo


#eORCA025.L75-IMHOTEP.GAI_y1996.1y_1y_TSMEAN.nc

