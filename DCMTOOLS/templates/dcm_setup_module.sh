#!/bin/bash
### This is a template file to copy and edit in your.bashrc file in order
# to set up your environement for using DCM4.
# This version assumes that you are using modules and that you already
# set up the module file for DCM ( see an example in templates/module_example
####################################################
#
#       Drakkar environment
####################################################
export MODULEPATH=$HOME/modules:$MODULEPATH
module load DCM
#
# UDIR is the user directory (usually home)
#  where the CONFIG_xxx directories and source files reside.
#  on rhodes or zahir
export UDIR=$HOME

# CDIR is the directory for compilation.
#   On zahir, we choose CDIR  as $WORKDIR.
export CDIR=$WORKDIR

#[ If DDIR not defined, CDIR is used in place 
# DDIR is the directory where data ( -S -MEAN -R -I ) will be stored
#   On zahir, we choose DDIR  as $WORKDIR.
export DDIR=$WORKDIR
#]

# PDIR is the HOME directory of the Production Machine.
#   On zahir, it is set to something like /homegpfs/rech/cli/rcli002
export PDIR=$HOME

# SDIR is the HOME directory of the Storage machine.
#  On zahir, it is set to something like /u/rech/cli/rcli002
export SDIR=/u/rech/cli/rcli099

#
##############################################

