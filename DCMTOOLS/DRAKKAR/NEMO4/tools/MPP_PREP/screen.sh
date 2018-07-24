#!/bin/bash
# screen.sh : script used for screening the processor.layout file
###########################################################
# $Date: 2010-11-17 22:59:43 +0100 (Wed, 17 Nov 2010) $
# $Rev: 216 $
# $Id: screen.ksh 216 2010-11-17 21:59:43Z forge $
###########################################################

usage() {
     echo USAGE: $(basename $0 ) [ -h ] [-f file ] nproc
     echo "       -h : this help page "
     echo "       -f file : take file instead of processor.layout "
     echo "       nproc : number of sea proc to look for "
     exit 0
        }

 file='processor.layout'

while getopts :hf: opt ; do
   case $opt in
     (h) usage ;;
     (f) file=$OPTARG ;;
     (\?) echo $(basename $0 )": -"${OPTARG}" option : not supported " 1>&2 ; usage ;;
   esac
done

shift $(( $OPTIND - 1 ))

n=$1

if [ ! $n ] ; then echo you must specify nproc ; echo ; usage ; fi
if [ ! -f $file ] ; then echo $file is missing ; echo ;  exit 1 ; fi

s=$(printf "ocean processors%19d\n" $1)
lst=$(grep -A 10 -B 3 -e "$s" $file)

echo " Searching for $n sea processors in $file ...."

echo $lst | sed -e 's/--/,/g' | awk -F',' '{ for (i=1 ; i<=NF ; i++) print $i }' |\
    awk ' BEGIN { printf " jpni  jpnj   jpi   jpj  jpi x jpj  proc  elim   sup\n" } {printf "%5d %5d %5d %5d %10d %5d %5d %8.5f\n" , $2, $4,  $6 , $8 ,$8*$6,  $13 ,  $18 , $81}'
