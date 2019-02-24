# DCM Compilation manual

  This document explains how to use DCM in order to create a model configuration,
from code installation to compilation. It assumes that you have already
installed and operational version of DCM, following instruction given in the
['Getting Started ' document](./dcm_getting_started.md).

## Important vocabulary and syntax:
#### CONFIG and CASE

  In DCM, we aim at minimizing the keywords for identifying a model configuration. In this sense, 2 keywords are essential :

 1. `CONFIG` : This is the name of a configuration, in general corresponding to a specific domain (like ORCA2, ORCA12, NATL025 etc...). In order to be even more specific, the number of levels can be specified as part of the `CONFIG` name, as in `ORCA025.L46` or `ORCA025.L75`. In this case, note the syntax with `.Lxx`, which is part of the DCM convention.
 1. `CASE` : This is the name of a numerical experiment, of course associated to a particular `CONFIG`. Hence, the sequence `CONFIG-CASE` clearly indentify a model realisation. (*e.g.*, ORCA12.L46-MAL90, ORCA025-KAB01 etc...). Note the you should use the `-` symbol between CONFIG and CASE.  In the DRAKKAR project we use an extra convention for coding the `CASE` name : The first letter identify the group or the machine where the configuration was produced, and possibly followed by the initial of the person in charge of the experiment, thus KAB stands for Kiel group, run performed by Arne, BAMT is a Brest run performed by Anne Marie etc.... Although it seems too simple, it turns out that this convention is quite usefull, especially when you have to revisit an old run !

#### File name convention:

  Once the concept of `CONFIG` and `CASE` is well understood, all files produced by a model experiment will use the following rule for the name :

```
    <CONFIG>-<CASE>_<TAG>_<GRID>.nc
```

  `TAG` is a string with the following format :`yYYYYmMMdDD.freq` (*e.g.* y2010m12d24.1d ). `freq` can be one of the id's used by XIOS for defining the output frequency: 1d, 1h, 1m etc.... 

 `GRID` is an indication of the content of the file, possibly linked with its position on the C-grid. (*e.g* gridT, gridU ... icemod etc...). 

Note that using DCM and DRAKKAR convention, special characters such as `. - _` are used as field separators. So it is very important not to use them for other purpose. This is why, in DRAKKAR, we have gridT, gridU files and not grid_T, grid_U as in standard NEMO output.

## Creating a new configuration
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

### 2. Editing templates in `$UDIR/CONFIG_<CONFIG>/<CONFIG>-<CASE>/`
