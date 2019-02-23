# Getting Started

## Pre-requisite :
  You need to have an operational XIOS library and ``xios_server.exe`` executable.  
You can refer to [this file](../DCMTOOLS/NEMOREF/xios_revision.md) for recomended XIOS revision for this DCM release.

## Road-map:
  Assume that your git repositories are located in some ``$DEVGIT`` (up to you!)
  1. Clone the DCM GitHub repository in ``$DEVGIT``    
```
    cd $DEVGIT
    git clone https://github.com/meom-group/DCM.git DCM_4.0
```
  2. Get the official NEMO code from NEMO-ST IPSL server.  
```
   cd $DEVGIT/DCM/DCMTOOLS/NEMOREF/      
   ./getnemoref
```
  3. Set up DCM module  
   It is highly recommended to set DCM as a module, in particular if you are working with different release of NEMO at the same time. Likely, ``module`` is available on any HPC center, and is very easy to install on a linux system. However, you still have the possibility to set up the environment variables in your ``.bashrc`` file.
   The steps for setting up DCM modules are:  
    a) create a ``modules`` directory in your ``HOME`` (for instance).  
    b) append this directory name to ``MODULEPATH`` in your ``.bashrc`` file:  
```
    export MODULEPATH=$MODULEPATH:$HOME/modules/
```
    c) create there a DCM directory:  
```
    mkdir $HOME/modules/DCM
```
    d) copy the [template module file](../DCMTOOLS/templates/module_example) into ``$HOME/modules/DCM/``  
```
   cp $DEVGIT/DCM/DCMTOOLS/templates/module_example $HOME/modules/DCM/4.0
```
    e) edit ``$HOME/modules/DCM/4.0`` to fit your settings. Basically only few lines have to be modified :  
```
set             version         4.0
set             alter_version   4.0
set             base_path       $::env(DEVGIT)/DCM_$version/DCMTOOLS
set             alter_path      $::env(DEVGIT)/DCM_$alter_version/DCMTOOLS
```

> Note that _version_ refer to the actual version you are working with (DCM_4.0 in this example).  
> _alter_version_ points to another DCM version, and is used in some comparison tools (see advanced usage)  
> For first time use, no problems to set _version_ and _alter_version_ to the same thing.  
> IMPORTANT : version must be coherent with the name you give to your DCM repository in ```$DEVGIT```  

