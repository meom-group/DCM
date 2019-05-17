# Porting ORCA025 from NEMO-3.6 to NEMO-4
### About the test configuration
We choose to port the very last version of ORCA025 ran in MEOM : ORCA025-GJM189.
 For the testing purpose the code is compiled with debug options on occigen:

     -i4 -r8 -O  -CB -fpe0 -ftrapuv -traceback -g  -fp-model precise -xAVX -fno-alias

### About creating ORCA025_domain_cfg.nc file : use of DOMAINcfg tool
### About running a simulation : use of RUNTOOLS
#### Problems of stability
  * explosion after 10 steps in the southern ocean under ice
#### Tests without ice
  * Compile w/o key_si3, and set nn_ice=0 in namelist.
    * compiles OK
    * blows up a 1rst step with crasy salinities (very high). sowaflup ~34 kg/m2/s ! 
  * Compile with key_si3, and set nn_ice=0 in namelist.

