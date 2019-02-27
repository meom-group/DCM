# RUNTOOLS user manual
## Forewords:

  If you read this manual is that you probably have successfully installed and compiled a NEMO configuration, using DCM ! (to make it short you have a `nemo4.exe` in the `$PDIR/RUN_<CONFIG>/<CONFIG>-<CASE>/EXE/` directory, corresponding to your settings (code and cpp keys)). Ideally you are now only missing some runtime files to make your first run:  
 * the **namelist** files (for ocean, ice, passive tracers etc...)
    - Since NEMO_3.6, NEMO uses 2 levels of namelists: `namelist_ref` where all variables are defined with default values (in general corresponding to ORCA2 setting), and `namelist_cfg` which only re-define values changed with respect to the reference. Although it make sense, we found that having those 2 levels files (though never have the full set of parameters in a single file) is not convenient and probably a source of error. In DCM, we choose to maintain a single namelist (having all the variables in it, like `namelist_ref`) with the right parameters for the working configuration. This choice is transparent to NEMO as we just copy the *ad-hoc* `namelist_ref` to  `namelist_cfg`. 
    - Reference namelist (usefull to know about all the parameters) can be found in the `cfgs/SHARED` sub directory.
    - Namelists are a very clever and comfortable means to set parameters in NEMO, but it is really critical to double check the namelist parameters, as it is quite easy to make mistakes !
    - When using DRAKKAR modifications of the code, and in case a namelist block is envolved, DRAKKAR code redefine an additional namelist block for its own parameters. The DRAKKAR namelist block has the same name as the standard one, but with `_drk` appended (*e.g* `&namlbc` (std) and `&namlbc_drk` (DRAKKAR parameters). Doing so, a namelist setup can be run with standard NEMO code as well as with DRAKKAR modificated code.
    - Refer to the [NEMO book](https://www.nemo-ocean.eu/wp-content/uploads/NEMO_book.pdf), for detailed description of the namelist parameters (**caution: link point to NEMO_3.6 !!**)
 * the **xml** files used by XIOS for model output.
 * **data** files such as initial conditions, domain configuration, atmospheric forcing etc...


  This document explains how to run a model configuration, installed and compiled with DCM.
