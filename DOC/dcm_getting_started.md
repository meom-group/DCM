# Getting Started

## Pre-requisite :
  You need to have an operational XIOS library and ``xios_server.exe`` executable.  
You can refer to [this file](../DCMTOOLS/NEMOREF/xios_revision.md) for recomended XIOS revision for this DCM release.

## Road-map:
  Assume that your git repositories are located in some ``$DEVGIT`` (up to you!)
  1. Clone the DCM GitHub repository in ``$DEVGIT``  
``
   cd $DEVGIT
   git clone https://github.com/meom-group/DCM.git 
``
  1. Get the official NEMO code from NEMO-ST IPSL server.  
``
   cd $DEVGIT/DCM/DCMTOOLS/NEMOREF/
   ./getnemoref
``
  1. Set up DCM module  
   It is highly recommended to set DCM as a module, in particular if you are working with different release of NEMO at the same time. Likely, ``module`` is available on any HPC center, and is very easy to install on a linux system. However, you still have the possibility to set up the environment variables in your ``.bashrc`` file.


 
