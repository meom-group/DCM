# Porting BSAS12 from NEMO-3.6 to NEMO-4
### About the test configuration
 With the purpose of testing and evaluating NEMO4, we choose a small configuration (Black Sea and Azov Sea, aka BSAS12) at 1/12 deg. This configuration
was recently developped and ran in the frame of NEMO_3.6.  
 For the testing purpose the code is compiled with debug options on occigen:

     -i4 -r8 -O  -CB -fpe0 -ftrapuv -traceback -g  -fp-model precise -xAVX -fno-alias

### About creating BSAS12_domain_cfg.nc file : use of DOMAINcfg tool
 With DCM version of DOMAINcfg, I  took the files `bathy_meter.nc`, `coordinates.nc` and nemo namelist used with 3.6 version. I just adapt nammpp namelist
block in order to run the tool on a single node. I submit a job on 28 cores (occigen BDW28) and obtained domain_cfg file splitted on 28 domains. 
I recombine the pieces to a full domain file with `nemo_rebuild` tool, and finally get `domain_cfg.nc`, which was renamed `BSAS12_domain_cfg.nc`. 

> Although all was OK with DCM version of DOMAINcfg, there are some harmless error messages appearing in the ocean.output file after the job. They can be
disregarded.
>
> I think that in domain_cfg.nc file, there is a lack of information for tracability. I recommend to add some global attributes to describe the input files used when creating the file.  
>
> NEMO by itself is looking to `cn_cfg` and `nn_cfg` global attribute. If not found they are respectively replaced by `UNKNOWN` and `-999999`.  
>       
>         ncatted -h -a cn_cfg,global,c,c,bsas  -a nn_cfg,global,c,l,12 BSAS12_domain_cfg.nc
>
>  will add those two attributes to the domain_cfg.nc file.
> For the tracability, I also think that at least information regarding the full names of the bathymetry and coordinates file (not just bathy_meter and 
coordinates) should be indicated as a global attribute. 
>
> But this is not enough to rebuild correctly the domain_cfg.nc file : We also need the namelist, or at least some key values used in the setting of the partial steps vertical grid (and for sigma coordinates as well). However, giving just the name of the namelist in a global attribute is probably not sufficient
but will be better than nothing (namelist are too easy to erase or change ! ).
> I have coded a small fortran program that write the namelist in the netcdf variables. It might be an option.

### About running a simulation : use of RUNTOOLS
#### _Fixes for out-of-bound error_
 After some adjustment, the RUNTOOLS part of DCM is working enough to launch a job from the CTL directory. However, with the debug options of the compiler, 
I found an out-of-bound error in `traadv_fct.F90`.  It happens in these loops :

       DO jj = 2, jpjm1                 ! 2nd order centered at top & bottom
         DO ji = fs_2, fs_jpim1
            ikt = mikt(ji,jj) + 1            ! w-point below the 1st  wet point
            ikb = mbkt(ji,jj)                !     -   above the last wet point
            !
            zwd (ji,jj,ikt) = 1._wp          ! top
            zwi (ji,jj,ikt) = 0._wp
            zws (ji,jj,ikt) = 0._wp
            zwrm(ji,jj,ikt) = 0.5_wp * ( pt_in(ji,jj,ikt-1) + pt_in(ji,jj,ikt) )
            !
            zwd (ji,jj,ikb) = 1._wp          ! bottom
            zwi (ji,jj,ikb) = 0._wp
            zws (ji,jj,ikb) = 0._wp
            zwrm(ji,jj,ikb) = 0.5_wp * ( pt_in(ji,jj,ikb-1) + pt_in(ji,jj,ikb) )  ! <=========== HERE
         END DO
      END DO
 
The index ikb-1 turns to be 0 on land, because mkbt(:,:) is set to 1 on land. I think that this is probably harmless as the results on land will be masked.
(It might be a bug if there are places in the model domain where there are only 1 wet vertical cell).
But, as it is very practical to be able to check the out-of-bound errors, I decided to make a small modification setting `ikb = MAX(mbkt(ji,jj),2)`.
And with this fix/tweak, there no more out-of-bound error at run-time.

#### _Problems when using `ln_bdy=.true.`_ 

A run with no-bdy, 55 cores for nemo, 1 xios (2 occigen nodes) ends up normally.   
The same configuration, now with bdy, just freezes (seems to be in dynspg_ts::mpp_lnk_bdy, according to pstack on compute node). 

The very same one with only 27 nemo and 1 xios (1 occigen node) ends up normally !
> **Investigations under progress !**

   | ln_bdy  |   cores  | land elim  |  run   |  Compil | Comment |
   |:-------:|:--------:|:----------:|:------:|:-------:|---------|
   |  Y      |  1 x 27  |    N       |   OK   |  Debug  |         |
   |  Y      |  1 x 55  |    N       | **KO** |  Debug  |  Freeze in dynspg_ts : in mpp_max |
   |  Y      |  1 x 55  |    N       | **KO** |   O3    |  Frezze in dynspg_ts : in mpp_ma |
   |  Y      |  11 x 5  |    N       |   OK   |   O3    |         |
   |  Y      |  8 x 9   |    Y       | **KO** |  Debug  |  Freeze in dynspg_ts : in lbc_lnk_bdy |
   |  N      |  8 x 9   |    Y       |   OK   |  Debug  |       |

 Case domain decomp 1x55 seems to have too small jpj ? With land proc elim, (8x9) it freezes differently.  
 With the NEMO4 version of lib_mpp, the readilility of the mpp_link routines is degraded :(. 

 After a long search ( 2 days ! ) I found that the problems araised from the fact that the BDY location is just besides a limit 
betwee two domains.  Changing the domain decomposition (not optimal ) makes it work... 

Back to the problem for finding a fix, I found a problem with the values of nbondi_bdy between 2 adjacent processor : In my case
proc 0 has nbondi_bdy = -1 ( meaning exchange with the east proc, but the eastern proc ( proc 1 ) has a nbondi_bdy= 2 meaning no exchange.
This form a deadlock because proc 0 is waiting on a mpprcv from proc 1 that never comes.  Setting manually (!!!) nbondi_bdy=1 for proc 1 fix
the problem ( in bdyinit.F90) 

The idea is to debug this busy routine in order to end up with correct values for nbondi_bdy ...



#### _Introducing the ice model (SI3)_


#### _Comments_ 

 At the end of a segment there is a run.stat.nc file holding Smin, Smax, Tmin, Tmax and ABS(SSH)max, ABS(U)max

 Test convection with npc : it works !
