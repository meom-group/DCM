# RUNTOOLS user manual
## Forewords:

  If you read this manual is that you probably have successfully installed and compiled a NEMO configuration, using DCM ! (to make it short you have a `nemo4.exe` in the `$PDIR/RUN_<CONFIG>/<CONFIG>-<CASE>/EXE/` directory, corresponding to your settings (code and cpp keys)). Ideally you are now only missing some runtime files to make your first run (and independently of using DCM or not ):  
 * the **namelist** files (for ocean, ice, passive tracers etc...)
    - Since NEMO_3.6, NEMO uses 2 levels of namelists: `namelist_ref` where all variables are defined with default values (in general corresponding to ORCA2 setting), and `namelist_cfg` which only re-define values changed with respect to the reference. Although it make sense, we found that having those 2 levels files (though never have the full set of parameters in a single file) is not convenient and probably a source of error. In DCM, we choose to maintain a single namelist (having all the variables in it, like `namelist_ref`) with the right parameters for the working configuration. This choice is transparent to NEMO as we just copy the *ad-hoc* `namelist_ref` to  `namelist_cfg`. 
    - Reference namelist (usefull to know about all the parameters) can be found in the `cfgs/SHARED` sub directory.
    - Namelists are a very clever and comfortable means to set parameters in NEMO, but it is really critical to double check the namelist parameters, as it is quite easy to make mistakes !
    - When using DRAKKAR modifications of the code, and in case a namelist block is envolved, DRAKKAR code redefine an additional namelist block for its own parameters. The DRAKKAR namelist block has the same name as the standard one, but with `_drk` appended (*e.g* `&namlbc` (std) and `&namlbc_drk` (DRAKKAR parameters). Doing so, a namelist setup can be run with standard NEMO code as well as with DRAKKAR modificated code.
    - Refer to the [NEMO book](https://www.nemo-ocean.eu/wp-content/uploads/NEMO_book.pdf), for detailed description of the namelist parameters (**caution: link point to NEMO_3.6 !!**)
 * the **xml** files used by XIOS for model output.
    - These files allow the link between NEMO and XIOS. At the end XIOS just need a file called `iodef.xml` where all the informations on both NEMO XIOS implementation and user requirements for output are described. Recent versions of XIOS (> 2) allow the use of separate 'sub' files which are used as sources in `iodef.xml`, which really improves the readibilty of the files.
    - In NEMO4 implementation, the xml sub files are :
      * `domain_def_nemo.xml`
      * `grid_def_nemo.xml`  
      * `field_def_nemo-oce.xml`  
      * `field_def_nemo-ice.xml`  
      * `field_def_nemo-pisces.xml`  
      * `file_def_nemo-oce.xml`  
      * `file_def_nemo-ice.xml`  
      * `file_def_nemo-pisces.xml`  
      * `context_nemo.xml`  
      * `iodef.xml`
    - Actually, domain_def, grid_def, field_def and iodef are distributed with NEMO and standard users might not change them.  The xml files requiring editing by the user are :
      * `file_def_nemo-xxx.xml` is where you define which fields/variables and at which frequency you want the model to output. It is really the manager of your model output. 
      * `context_nemo.xml` is a kind of menu where you set the file_def components (ocean, ice, tracers ...) you are using.
    - More informations will be given for advanced users willing non standard variables or defining sub domains.
 * **data** files such as initial conditions, domain configuration, atmospheric forcing etc...
    - In NEMO4, most of the data file are defined as a namelist parameter. 


 In order to run NEMO, the principle is to have all the runtime files in a common directory, together with `nemo4.exe` and `xios_server.exe` and issue (for instance -- there are variant according to the HPC system--) a command like :

```
   mpirun -np $XIOS_CORES xios_server.exe : -np $NEMO4_CORES nemo4.exe 
```

  Where `XIOS_CORES` is the number of cores dedicated to XIOS, and `NEMO4_CORES` is the number of cores dedicated to NEMO4. This statement launch 2 executables with the MPMD (Multiple Program Multiple Data) paradigm.

  DCM's runtools consist of a collection of scripts, handling all the machinery required to produce a long simulation (chaining elementary segments of run), dealing with the model output, restart files etc...  Although it works even for simple configurations (such as test cases or idealized cases), it is primarily designed for complex realistic cases (which explains the relative complexity of the tools).


  This document explains how to run a model configuration, using DCM runtools.
