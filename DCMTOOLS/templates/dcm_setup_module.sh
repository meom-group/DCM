#!/bin/bash
### This is a template file to copy and edit in your.bashrc file in order
# to set up your environement for using DCM4.
# This version assumes that you are using modules and that you already
# set up the module file for DCM ( see an example in templates/module_example

# for memory:
# export WORKDIR=????   # if not set by the system
# export STOREDIR=????  # Access to archiving system
####################################################
#
#       Drakkar environment
####################################################
# Module set up :
export MODULEPATH=$HOME/modules:$MODULEPATH
module load DCM
#
# UDIR is the user directory (usually under $HOME)
#  where the CONFIG_xxx directories and source files reside.
export UDIR=$HOME/CONFIGS

# CDIR is the directory for compilation.
#   for a given CONFIG-CASE, compilation will use $CDIR/W<CONFIG>-<CASE>/
export CDIR=$WORKDIR

#[ If DDIR not defined, WORKDIR is used in place 
# DDIR is the directory where data ( -S -MEAN -R -I ) will be stored
# 
export DDIR=$WORKDIR
#]

# PDIR is the HOME directory of the Production Machine.
#   
export PDIR=$HOME/RUNS

# SDIR is the HOME directory of the Storage machine (independent machine) or
#       an NFS link to your archiving system
#  
export SDIR=$STOREDIR

#
##############################################

