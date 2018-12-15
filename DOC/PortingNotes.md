# From NEMO_3.6 to NEMO_4.0,  from NEMODRAK to DCM
>  by J.-M. Molines  [MEOM/IGE](http://meom-group.github.io/) (Grenoble)
______________________

## **_Overview :_**
 Version 4.0 of NEMO is a major release.  It is the result of a simplification process which was decided by the NEMO consortium
in order to help the maintenance of the code, cleaning obsolete features, revisiting the user interface, adding/replacing
some new numerical schemes, etc... 

 Not surprisingly, the porting from 3.6 to 4.0 requires some work, at different levels.  This document was written 
while performing this porting, together with NEMO4 discovery. After a reminder on DCM concept, it holds basically 
two parts :     
   *  (A)  corresponds to all the changes in the code itself, and the compilation procedure. 
   *  (B) is dedicated to the impacts of the changes at run-time.

This work is done in the frame of the development of a fully new version of the Drakkar Config Manager (DCM). 

##  **_Reminder of the DCM concept:_**
 DCM is an integrated software tool allowing an easy deployment of a model configuration based on NEMO. For the model code,
it allows an easy management and tracking of code modifications, that may arise for instance when a new feature is being 
developped, or for a particular model configuration. It also provides (independently) usefull environnement for the production 
of simulations on HPC centers. Therefore, DCM is divided in two parts: **DCMTOOLS** and **RUNTOOLS**. This document only refer to the DCMTOOLS part. See [user manual for RUNTOOLS](rt_user_manual.md)

When building the code for a configuration, we use 3 layers, which are piled up to produce the final code:

The first layer
is the NEMO reference version, corresponding to a svn revision on the [NEMO system team repository](http://forge.ipsl.jussieu.fr/nemo/svn) . 
A DCM user must register as a declared user of NEMO before uploading the reference version. This layer is called
**NEMOREF** in the DCM jargon.

The second layer is the DRAKKAR layer, where only modified fortran modules are kept on the same layout as NEMOREF. This
layer may eventually be empty, but in practical it holds the DRAKKAR enhancement, not yet return to the NEMO system team, or
pieces of code that are just for DRAKKAR configurations. This layer is called **DRAKKAR** in the DCM jargon.

The third layer is user layer, where only modified fortran modules with respect to **NEMOREF** and/or **DRAKKAR** are kept,
a private directory of the user. Here again, the same NEMO layout is concerved in order to ease the building of the target
code. In the DCM jargon this layer is called the **CONFIG-CASE**, where CONFIG define the name of the configuration, 
and CASE is the name of a particular experiment performed with this configuration. (*e.g*  ORCA12-MJM189).

NEMOREF layer is fixed at a certain NEMO revision, and is the exact checkout of this revision (allowing tracability). DRAKKAR layer is maintained by DCM local administrator ( *e.g* J.M. Molines for the MEOM team at IGE). The CONFIG-CASE layer is at the responsability of the proper end user of the system.  It is initally an empty layer, holding only the tools for compilation.
If necessary, the user populate this layer from NEMOREF or DRAKKAR, then proceed to the modifications she/he wants to bring
to the system. 

DCM proposes a set of tools (mainly bash scripts) in order to deal with the file layout automatically. Since the first version,
we maintain almost the same syntax for these tools, maintaining a friendly user interface through the NEMO versions. See the [user manual](dcm_user_manual.md).


## **_A: Fortran code, compilation etc_**
 In DCM_4.0, the DRAKKAR modifications on fortran modules are all isolated from the standard code by a sequence of

    #if defined key_drakkar
      drakkar code ...
    #endif
 So that, it is very easy to identify the DRAKKAR modifications (using `meld` for instance). When a modification of a
namelist block is introduced ( *e.g*  add extra variables), we took the convention to put use a new namelist block whose
name is just the same with just `_drk` appended at the end of the block name. For example :

    &namtra_dmp
    ... NEMOREF std namelist
    /
    &namtra_dmp_drk
    ... DRAKKAR extension
    /
 Doing so, a user configuration can be run very easily switching off the DRAKKAR/CONFIG-CASE modifications, at least 
for testing purposes. This is a clear improvement with respect to previous versions of DCM (NEMODRAK).

### Code
 Compared with previous releases, NEMO_4.0 uses a new, completly different directory tree for the file layout. Although it 
changes a lot the old-users habits, it is much more coherent and profesional !

    arch/   --> architecture files for various computer/compilers
    cfgs/   --> configurations directory
        /SHARED  --> common namelist, xml files etc to all configs
        /<supported configurations>  --> for example : ORCA2_ICE_PISCES or AGRIF_DEMO
    doc/    --> latex sources and script for the NEMO_book
    ext/    --> external packages used by NEMO 
       /AGRIF    --> core of the AGRIF system ( not NEMO dependent)
       /FCM      --> perl scripts for FCM
       /IOIPSL   --> IOIPSL library
    mk/    -->  shell scripts and data files using FCM command for code compilation
    src/    -->  NEMO source code
       /ICE      --> source code for ICE engine : SI3 ( new ice model) pronounce : Sea Ice Cube
       /NST      --> source code for  NeSTing nemo  agrif interface 
       /OCE      --> source code for OCEan engine
       /OFF      --> source code for OFF-line passive tracers
       /SAO      --> source code for StandAloneObs module
       /SAS      --> source code for StandAloneSurface module
       /TOP      --> source code for TOP (passive tracer engine).
    tests/      --> configurations for tests cases ( new feature in NEMO 4)
         /<tests cases> : CANAL  ICEDYN  ISOMIP  LOCK_EXCHANGE  OVERFLOW   VORTEX  WAD 
    tools/        --> source code for various companion tools
         /<tools> : BDY_TOOLS/  DOMAINcfg/  MISCELLANEOUS/  NESTING/  REBUILD_NEMO/
                  SECTIONS_DIADCT/  WEIGHTS/ DMP_TOOLS/  GRIDGEN/    MPP_PREP/
                  OBSTOOLS/  REBUILD/ SCOORD_GEN/  SIREN/
         /maketools : compiling scripts (fcm based) for these tools

 As of December 2018, vith regard to `src/` **DRAKKAR** only contains modifications concerning 
 `OCE/` and `ICE/` as detailed below.

#### _src/OCE/ DRAKKAR content:_

  directory | module | reason|
  :----------|:--------|:-------|
  OCE/    | timing.F90 |  use more digits in format (for large simulations). |
  | "     | nemogcm.F90 | print time stamp at every step in std out <br> use flag BBBBBB and :( :( :( in ocean output <br> instead of AAAAAAAAA  |  
   OCE/ASM  |    |   |
   OCE/BDY  |    |   |
   OCE/C1D  |    |   |
   OCE/CRS  |    |   |
   OCE/DIA  |    |   |
   OCE/DIU  |    |   |
  OCE/DOM | domain.F90 | use segment number in the restart directory in order to ease production of big simulations <br> and avoid further renaming of restart files.    |
  |  "   | dtatsd.F90 | differentiate initial field  and damping fields |
  |  "    | dommsk.F90 | Allow for a 2D shlat value, read in external file. Add namlbc_drk.   | 
   OCE/DYN |     |   |
   OCE/FLO |     |   |
   OCE/ICB |     |   |
  OCE/IOM | restart.F90 | implement simple restart file names using segment numbers, avoiding renaming of restart files. |
  |  "    | iom.F90 | Implement new possible file_id in xml file ( extra iom_update_file_name ) |
   OCE/LBC |     |   |
   OCE/LDF |     |   |
   OCE/OBS |     |   |
  OCE/SBC | sbcblk.F90 | implement P. Mathiot phd katabatic winds enhancement. Use namsbc_blk_drk. |
  |  "    | _shapiro.F90_| New code used in sbcssr for filtering model Surface fields.  |
  |  "    | sbcssr.F90 | Add 3 new features in the sea surface restoring: <br>  1.  Permit local enhancement of the ssr strength. <br> 2. Add an option for limiting the near coast restoring.  <br> 3. Add an option for filtering SS fields before computing the mismatch with smooth observations or climatologies.  |
   OCE/STO |     |   |
  OCE/TRA | trabbl.F90 |  Implement H. Hervieux  K-criteria to replace (option ) the H-criteria for activation of BBL. |
  | "     | tradmp.F90 |  Allow for code modification of the restoring coeficient. Ideally this should be coded in a usr_tra_dmp ...  <br> Use namtra_dp_drk namelist |
   OCE/TRD |     |   |
  OCE/USR | usrdef_fmask.F90 |  Port change of lateral friction is some straits for ORCA025. <br> not that shlat2d is also implemented in drakkar modifications in dommsk|
  OCE/ZDF | zdftke.F90 |  Changes regarding influence of ice coverage on differents term not coded yet. Need more evalutation of pro and cons. |
  | "     | zdfdrg.F90 |  Implement namelist definition for the boost drag files. (instead of hard coded names). |
  

#### _src/ICE/ DRAKKAR content:_

  directory | module | reason|
  :----------|:--------|:-------|
  ICE/    | icerst.F90 |  Use drakkar restart file name allowing <br> restarting the code whithout renaming. |
  | "     | icestp.F90 |  Differentiate input and output restart <br> by their name, using segment number <br> and also different directories. |


#### _Comments:_
 In NEMO4 a driving idea was to eliminate as much as possible the CPP keys, and replace them by logical keys set in the namelist. 

Another idea was to eliminate from the code all particular cases regarding some configuration (for instance tweak for orca2 etc...). This latter point means that all the tweaks are reported in the input files (domain_cfg which replace `coordinates.nc`, `bathy_meter.nc`, and also provide the vertical metrics. In NEMO4, the depths are computed from the vertical metrics read in the domain_cfg.nc file. When really not easy to put tweaks in an external file, a USR/ directory is ready for saving user modules that may perform the user-wanted modifications. In the present version only usrdef_fmask.F90 is used in this sense for changing lateral friction is particular points of the grid ( straits etc).


## **_B: Run time changes_**

### Namelists

### domain_cfg.nc file

### xml files

#### _Comments:_
 At the end of a segment there is a run.stat.nc file holding Smin, Smax, Tmin, Tmax and ABS(SSH)max, ABS(U)max

 Test convection with npc : it works !
