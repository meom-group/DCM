#!/bin/bash
# class = @File management tools@

usage() {
   echo
   echo "USAGE : $(basename $0 ) directory* [-h] "
   echo
   echo "  PURPOSE :"
   echo "     Give the size, number of inodes and mean file size in the directoried given"
   echo "     as arguments."
   echo
   echo "  ARGUMENTS: "
   echo "     A list of directories to look for  "
   echo 
   echo "  OPTIONS:"
   echo "     -h : Display this help page "
   echo
   exit 0
        }

if [ $# =  0 ] ; then
   usage
fi

while getopts :h opt ; do
   case $opt in
     (h) usage ;;
     (*) usage ;;
   esac
done

for d in $* ; do
   if [ -d $d ] ; then
     size=$( du -sk --apparent-size $d | awk '{print $1/1024/1024.}')
     inod=$( inode.ksh $d | awk '{print $1}')
     if [ $inod != 0 ] ; then
        mfsz=$( echo $size $inod | awk '{ print $1/$2 }' )
  
        echo -n $d " : "
        echo    $inod $size Gb  $mfsz Gb
     fi
   fi
done
  
   
