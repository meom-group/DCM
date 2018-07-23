#!/bin/bash
### This is a template file to copy and edit in your.bashrc file in order
# to set up your environement for using DCM4.
# This version assumes that you are not using modules and mimic the module 
# action in .bashrc/.
# Of course, without the modules, the file must be edited if you change the version
# of DCM ( and this is very tedious). I recommend the use of modules !
####################################################
#
#       Drakkar environment
####################################################
# HOMEDCM is the directory where the DCM resides
# in general  it points to a unique directory for each machine
#  on rhodes
export HOMEDCM=/some/git/repository/DCM_4.0/DCMTOOLS
# RUNTOOLS point to the RUNTOOLS directory assiociated with this level of DCM
export RUNTOOLS=/some/git/repository/DCM_4.0/RUNTOOLS

# REFDIR is the reference NEMO  directory
#  on zahir
export REFDIR=$HOMEDCM/NEMOREF/NEMO4

# CUSDIR is the DRAKKAR customized directory
#     holding the permanent customized files not yet under cvs, but valid
#     for all the users (in the modipsl tree) will overwrite the REF during install
#  on zahir
export CUSDIR=$HOMEDCM/DRAKKAR/NEMO4

export PATH=${PATH}:$HOMEDCM/TOOLS
#
##############################################

## The following is to be done  always ( module/no modules)
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

