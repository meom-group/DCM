# Porting ORCA12 from NEMO-3.6 to NEMO-4
### About the test configuration
We choose to port the very last version of ORCA12 ran in MEOM : ORCA12-MJM189
 For the testing purpose the code is compiled with debug options on occigen:

     -i4 -r8 -O  -CB -fpe0 -ftrapuv -traceback -g  -fp-model precise -xAVX -fno-alias

### About creating ORCA12_domain_cfg.nc file : use of DOMAINcfg tool
### About running a simulation : use of RUNTOOLS
#### Problems of stability
#### Tests without ice
