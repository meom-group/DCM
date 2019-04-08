# DCM Compilation manual

  This document explains how to use DCM in order to create a model configuration,
from code installation to compilation. It assumes that you have already
installed an operational version of DCM, following instructions given in the
['Getting Started ' document](./dcm_getting_started.md).

## Important vocabulary and syntax:
#### CONFIG and CASE

  In DCM, we aim at minimizing the keywords for identifying a model configuration. In this sense, 2 keywords are essential :

 1. `CONFIG` : This is the name of a configuration, in general corresponding to a specific domain (like ORCA2, ORCA12, NATL025 etc...). In order to be even more specific, the number of levels can be specified as part of the `CONFIG` name, as in `ORCA025.L46` or `ORCA025.L75`. In this case, note the syntax with `.Lxx`, which is part of the DCM convention.
 1. `CASE` : This is the name of a numerical experiment, of course associated to a particular `CONFIG`. Hence, the sequence `CONFIG-CASE` clearly indentify a model realisation. (*e.g.*, ORCA12.L46-MAL90, ORCA025-KAB01 etc...). Note that you should use the `-` symbol between CONFIG and CASE.  In the DRAKKAR project we use an extra convention for coding the `CASE` name : The first letter identify the group or the machine where the configuration was produced, and possibly followed by the initial of the person in charge of the experiment, thus KAB stands for Kiel group, run performed by Arne, BAMT is a Brest run performed by Anne Marie etc.... Although it seems too simple, it turns out that this convention is quite usefull, especially when you have to revisit an old run !

  When using an HPC system, in some centers, there are up to 3 different machines for specific usage. For instance, on vectorial computer, we use to have cross compiler, meaning that the code was compiled on a different machine than the machine where it was run later. Also, it happens that the long term archiving system  is associated to a dedicated front-end computer. Although on more recent systems, those 3 different 'logical' machines corresponds to the same hardware, DCM allows for this diferenciation. In DCM jargon we speak about :
 1. `File machine` : This is the machine associated with long term archiving. Directories on this machine have an `F_` prefix in the running scripts.
 1. `Compile machine` : This is the machine used for compiling NEMO. The file system associated with this machine is `$CDIR` and is used for temporary files produced during the compilation process. It is typically located on a `WORKING directory`.
 1. `Production machine` : This is the machine where the code will be run, and where the running script will be copied. Also, run time data and model output will use the file system associated with this production machine. Directories on this file system have a `P_`prefix in the running scipts.

#### File name convention:

  `CONFIG` and `CASE` are used throughout all the production process, in particular for naming model output files. The rule is as follows :

```
    <CONFIG>-<CASE>_<TAG>_<GRID>.nc
```

  `TAG` is a string with the following format :`yYYYYmMMdDD.freq` (*e.g.* y2010m12d24.1d ). `freq` can be one of the id's used by XIOS for defining the output frequency: 1d, 1h, 1m etc.... 

 `GRID` is an indication of the content of the file, possibly linked with its position on the C-grid. (*e.g* gridT, gridU ... icemod etc...). 

Note that using DCM and DRAKKAR convention, special characters such as `. - _` are used as field separators. So it is very important not to use them for other purpose. This is why, in DRAKKAR, we have gridT, gridU files and not grid_T, grid_U as in standard NEMO output.

## Creating a brand new configuration (starting from nothing!)
### 1. Creating a new set of directories.

  This task will be performed by specific dcm_tools, installed with DCM. At first, you need to choose a `CONFIG` and `CASE` name. Then :

```
   dcm_mkconfdir_local <CONFIG> <CASE>
```

  This command wil create a set of directories/subdirectories that will be used for installing the code and for run production later on.

  Root names for the set of directories are those defined when installing DCM on your system. The created directories are :

 * on `$UDIR`:  
   `CONFIG_<CONFIG>/<CONFIG>-<CASE>/`
 * on `$PDIR`:  
   `RUN_<CONFIG>/<CONFIG>-<CASE>/EXE`  
   `RUN_<CONFIG>/<CONFIG>-<CASE>/CTL`  
 * on `$WORKDIR`:  
   `<CONFIG>/<CONFIG>-I`  
   `<CONFIG>/<CONFIG>-<CASE>-S`  
   `<CONFIG>/<CONFIG>-<CASE>-MEAN`  
 * on `$SDIR` :  
   `<CONFIG>/<CONFIG>-I`  
   `<CONFIG>/<CONFIG>-<CASE>-S`  
   `<CONFIG>/<CONFIG>-<CASE>-MEAN`  

    > At this point, just note that if `$SDIR` is located on a remote system and is not visible (*e.g.* NFS mount) on the login nodes, you may use the alternative command `dcm_mkconfdir_remote <CONFIG> <CASE>`

  For compilation purposes (this manual) only the directories created on `$UDIR` will be adressed. The other ones will be used and described in the [runtools manual](dcm_rt_manual.md).

### 2. Editing templates in `$UDIR/CONFIG_<CONFIG>/<CONFIG>-<CASE>/` (aka *confcase* directory)

  **From now onwards, all actions must be performed while staying in the *confcase* directory.**

```
   cd $UDIR/CONFIG_<CONFIG>/<CONFIG>-<CASE>/
```

  After step 1, some template files are copied in *confcase* directory and require some editing. Going to this directory, you see that a full NEMO tree hierachy of directories  (empty) has been created, (arch, cfgs, ext, src) and 2 template files, `CPP.keys` and `makefile`, are present. 
Another empty file whose name corresponds to the version name of the actual DCM (*e.g.* DCM_4.0) has been touched; it is there just for information.  You need to edit the following files:
 * `CPP.keys` : In this file you should define the cpp keys used for your configuration. The format is free. Any line starting with a # is considered as comment. The template CPP.keys, contains all possible keys for NEMO/DCM. You should clean the file in order to keep only the keys you need. You can refer to the NEMO book for relevant information about the cpp keys and their meaning. In the making process, this file will be used to create an *ad-hoc* file understandable by the compiling system (fcm).
 * `makefile` :  This file uses functionalities of GNU make. It controls all the processes of code installation and compilation. The header part of this file (clearly identified) need editing to match your requirements. All customisable variables have a corresponding self explanatory comment to help the user. For the basic usage, it is mandatory to have :  
    * `CONFIG = <CONFIG>` and  `CASE = <CASE>`  
 This is normally done automatically at step 1.  
    * `REFONLY = no/yes` : This fix the use of pure NEMO code ('yes') or DRAKKAR tainted code ('no'). Note that in order to use DRAKKAR code functionalities, you should have both `REFONLY = no` and `key_drakkar` defined in `CPP.keys` file.  
    * `MACHINE = your-machine-arch` 
 This is linked with the FCM compiling system used by NEMO (and XIOS). It requires a so called architecture file that contains usefull information about compiler name and options, libraries to be used etc... (See below for more details). All the architecture files are named such as `arch_your-machine-arch.fcm`. A full list of already available architecture can be displayed with the command `dcm_lsarch`. At the `makefile` level, you need to put `your-machine-arch` in the `MACHINE` variables.

    * set used NEMO components according to your configuration. This is done by setting either 'use' or 'notused' to the variables:  

      ```
      OCE = 'use'     # Ocean component
      ICE = 'use'     # Ice component
      TOP = 'notused' # Passive tracer component
      OFF = 'notused' # Passive tracer off-line
      NST = 'notused' # Nesting component (AGRIF, requires cpp key key_agrif
      SAS = 'notused' # Stand Alone Surface module
      MY_SRC='notused' # Personal code (likely useless when using DCM)
      ```
   Usefull but not mandatory:  
    * `NCOMPIL_PROC` defines the number of processors/cores used for parallel compilation. If greater than 1 it reduces the elapsed time needed for compilation.  
    * `GIT = no/check` : If set to check, an history file with information on the revision of NEMO and DCM is updated evry time the code in installed for this particular CONFIG-CASE. It is very usefull for tracability.   
    * Other variables set by default to 'none' will be discussed later on, when needed for more advanced usage.  

 * `arch_MACHINE.fcm`:  It is very likely than you will have to edit your architecture file to fit the architecture of the machine you are using. `dcm_lsarch` shows all the available architecture files available with NEMO and DCM. There are so many that it is not easy to choose. A good practice is to copy one of those (for instance `arch_X64-OCCIGEN.fcm`) into the arch/ sub directory of your *confcase* directory changing its name to avoid conflicts with standard arch files. Then you can edit it for your own customization. For example the content of a customized arch file can be :  

   ```
   %NCDF_HOME           $NETCDFHOME
   %NCDF_HOME_FORTRAN   $NETCDFFHOME
   %HDF5_HOME           $NETCDFFHOME
   %XIOS_HOME           $WORKDIR/DEV/xios-2.5-dev_olga
   %OASIS_HOME          $WORKDIR/now/models/oa3mct

   %NCDF_INC            -I$NETCDF_INCDIR -I$NETCDFF_INCDIR
   %NCDF_LIB            $NETCDF_LDFLAGS $NETCDFF_LDFLAGS
   %XIOS_INC            -I%XIOS_HOME/inc
   %XIOS_LIB            -L%XIOS_HOME/lib -lxios -lstdc++
   %OASIS_INC           -I%OASIS_HOME/build/lib/mct -I%OASIS_HOME/build/lib/psmile.MPI1
   %OASIS_LIB           -L%OASIS_HOME/lib -lpsmile.MPI1 -lmct -lmpeu -lscrip
   
   %CPP                 cpp
   %FC                  mpif90 -c -cpp
   %FCFLAGS             -i4 -r8 -O3 -fp-model precise -xAVX -fno-alias -init=zero -init=arrays
   %FFLAGS              %FCFLAGS
   %LD                  mpif90
   %FPPFLAGS            -P  -traditional
   %LDFLAGS             -lstdc++
   %AR                  ar
   %ARFLAGS             rs
   %MK                  gmake
   %USER_INC            %XIOS_INC %OASIS_INC %NCDF_INC
   %USER_LIB            %XIOS_LIB %OASIS_LIB %NCDF_LIB

   %CC                  cc
   %CFLAGS              -O0
   ```

  In this particular case, HOME directories for netcdf, hdf5, as well as include and lib directories are set using environment variables coming from the corresponding module ( `module show netcdf` for instance can help to know the (machine dependent) name of the environment variables. The attention is particularly drawn on `%XIOS_HOME` that must point to the directory where you installed XIOS as a pre-requisite.
 
 In general with few trial/error the good choice is found... It is not that difficult !

### 3. Installing the code.

  At the end of this step, you will have a `WORK` directory with the full fortran code and the tools for compilation. Being still in *confcase* this step is achieved by the single command:

```
   make install
```

 It is during this step that the NEMOREF, DRAKKAR and CONFIG-CASE layers of code are hierarchically copied to build the `WORK`. [Note that for now, there is no fortran code in *confcase*]. After this step, a link to `WORK` is available in *confcase* (so that you can have a global view on the code).

 This is also during this step that the intallation history is appended to `install_history` file in *confcase* (if you have `GIT = check` in `makefile`).

### 4. Importing pieces of code to *confcase*

  One strength of DCM is the faculty of importing NEMO fortran files into *confcase*, previous your own modifications. You always work in your *confcase* (private configuration). To do so, use the command :

```
   dcm_getfile  <code.F90>
```

Where \<code.F90\> is a NEMO fortran module (for instance `nemogcm.F90`). This command copy the target code file within the local NEMO tree (in src/) and creates a link for this file into *confcase*. Then you can edit/modify this code for your own purposes. It will be taken into account at the next compilation.

> Note that if for some reason you do not need anymore this code file into your *confcase*, you **must** use the `dcm_rmfile <code.F90>` command in order to erase the code from the local directory **AND** restore the corresponding file from NEMOREF+DRAKKAR into `WORK`

### 5. Compile the code :

 This is just the result of :  

```
   make
```

 At the end of this step, if no error are encountered at compile time, a link called `nemo4.exe` is created in `$PDIR/RUN_<CONFIG>/<CONFIG>-<CASE>/EXE/` directory, and a copy of `CPP.keys` is also put there (for the runtime processing), and you are done with code compilation.

## Compiling a NEMO REFERENCE CASE:

  NEMO system team maintains a reduced number of reference configurations covering almost all the NEMO capabilities.  In NEMO4(@rev 10650 ! ) there are 10 reference configurations:

```
AGRIF_DEMO
AMM12
C1D_PAPA
GYRE_BFM
GYRE_PISCES
ORCA2_ICE_PISCES
ORCA2_OFF_PISCES
ORCA2_OFF_TRC
ORCA2_SAS_ICE
SPITZ12
```

 In order to build a reference configuration just follow the steps (*e.g.* build GYRE_PISCES in GYRE-tst) :
   1. dcm_mkcondir_local GYRE-tst  
   1. ``` cd $UDIR/CONFIG_GYRE/GYRE-tst```
   1. edit `makefile` as previously for the basic things and :
       * set NEMO_REF_CASE to one of the official references (see above)
   1. prepare the reference configuration:
       * `make nemo_ref_case` : This will copy the official setup at the right  place (CPP.keys, cfg files etc ... )
   1. Compile the code
       * `make install `
       * `make`
  

## Compiling a NEMO TEST CASE:

  NEMO system team maintains a reduced number of test cases, coresponding to simplified model configuration. These test cases are intented to illustrate theoretical concept or a specific parameterization.  In NEMO4(@rev 10650 ! ) there are 10 test cases.

```
BENCH
CANAL
ICE_ADV1D
ICE_ADV2D
ICE_AGRIF
ISOMIP
LOCK_EXCHANGE
OVERFLOW
VORTEX
WAD
```

 In order to build a NEMO test case, just follow the steps (*e.g.* build WAD  in WAD-tst) :
   1. dcm_mkcondir_local WAD-tst  
   1. ``` cd $UDIR/CONFIG_WAS/WAD-tst```
   1. edit `makefile` as previously for the basic things and :
       * set TEST_CASE to one of the official NEMO test cases (see above)
   1. prepare the reference configuration:
       * `make test_case` : This will copy the official setup at the right  place (CPP.keys, cfg files etc ... )
   1. Compile the code
       * `make install `
       * `make`
 

## Cloning and modifying an already existing configuration:

  Actually, it very often happens that you already have a working `CONFIG-CASE` experiment and that you just want to clone as a starting point for further experiments.  DCM allows to do so very easily. 
  * As in the general case, you need to create the directory structure for the new `CONFIG-CASE`, and as previously this is done with :  

    ```
    dcm_mkconfdir_local <CONFIG> <CASE>
    ```

  * Then in the *confcase* directory, when editing `makefile`, you just need to set:  

    ```
    PREV_CONFIG = /full/path/to/CONFIG-CASEtoClone
    ```

    and use the 'magic' copyconfigall makefile target:  

    ```
    make copyconfigall
    ```

    All the fortran files from `CONFIG-CASEtoClone` are now in your new config, as well as the architecture file, the `CPP.keys` file and `makefile.prev`, exactly corresponding to the makefile used in `CONFIG-CASEtoClone`. You might rename it to `makefile` but remember to fix CONFIG and CASE to your new configuration names.

    > You can set PREV_CONFIG to any *confcase* directory (even from another user) provided you have read access to this directory, and provided that this *confcase* was created with the same version of DCM.

  * Then installation and compilation are achieved as in the general case:  

    ```
    make install && make 
    ```

## Extra tricks :

  * In some cases, you need to re-do an installation of the code in order to be sure of what is really in your WORK... This must be done by the sequence :  

    ```
    make cleaninst && make install
    ```

    This is usefull in particular if for some reasons you mixed up 2 different versions of DCM or NEMO ...
  * If you manually add files in the src/ tree, or if you add arch files into arch/ you can rebuild the links in *confcase* by using :  

    ```
    make links
    ```

  * In general, all makefile targets are listed with :

    ```
    make help
    ```






