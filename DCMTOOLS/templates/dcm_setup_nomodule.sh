#!/bin/bash
### This is a template file to copy and edit in your.bashrc file in order
# to set up your environement for using DCM4.
# This version assumes that you are not using modules and mimic the module 
# action in .bashrc/.
# Of course, without the modules, the file must be edited if you change the version
# of DCM ( and this is very tedious). I recommend the use of modules !
# for memory:
# export WORKDIR=????   # if not set by the system
# export STOREDIR=????  # Access to archiving system
####################################################
#
#       Drakkar environment
####################################################
# HOMEDCM is the directory where the DCM resides
# in general  it points to a unique directory for each machine
#
export HOMEDCM=/some/git/repository/DCM_4.0/DCMTOOLS
# RUNTOOLS point to the RUNTOOLS directory assiociated with this level of DCM
export RUNTOOLS=/some/git/repository/DCM_4.0/RUNTOOLS

# REFDIR is the reference NEMO  directory
# 
export REFDIR=$HOMEDCM/NEMOREF/NEMO4

# CUSDIR is the DRAKKAR customized directory
#  
export CUSDIR=$HOMEDCM/DRAKKAR/NEMO4

export PATH=${PATH}:$HOMEDCM/bin
#
########################################################################
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
