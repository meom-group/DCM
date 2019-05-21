# Porting ORCA025 from NEMO-3.6 to NEMO-4: test on ADA
### About the test configuration
We choose to port the very last version of ORCA025 ran in MEOM : ORCA025-GJM189.
 For the testing purpose the code is compiled with debug options on occigen:

     -i4 -r8 -O  -CB -fpe0 -ftrapuv -traceback -g  -fp-model precise -xAVX -fno-alias

### About creating ORCA025_domain_cfg.nc file : use of DOMAINcfg tool
   domain_cfg.nc file was created using the followong configuration files

   and the namelist 
### About running a simulation : use of RUNTOOLS
#### Problems of stability
  * explosion after 10 steps in the southern ocean under ice
  * reduce time step to 720 sec (!!) and explosion occurs at step 20.
#### Tests without ice
  * Compile w/o `key_si3`, and set `nn_ice=0` in namelist.
    * compiles OK
    * blows up a 1rst step with crasy salinities (very high). sowaflup ~34 kg/m2/s ! 
  * Compile with `key_si3`, and set `nn_ice=0` in namelist.
    * In this case, as far as `key_si3` is defined, `nn_ice` is reset to 2 despitre its value in the namelist.
    * I try to force `nn_ice=0` even with `key_si3` defined : Floating point execption in fwb
  * So : no easy way to run the code without ice (need extra debuging).
### Playing with namelist parameters :
  * using FCT second order (instead of 4th order) : explosion as well
  * using horizontal diffusion instead of isopycnal : better but still explosion
  * using all the above + UBS adv scheme : 180 steps OK !
> It seems that the instability comes from TS condition in the southern part of the southern ocean, where no data are in levitus and where the drowning procedure produce spurious fronts. However, in older version of NEMO it was not a problem !
### Working on the initial condition for TS
  * Note that when nn_istate=1 in the namelist, the initial condition which is output is eventually the model state after the 1rst time step. (In particular, velocities are not zero).
  * using the very last version of [sosie3](https://github.com/brodeau/sosie).
    

