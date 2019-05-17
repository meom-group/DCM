# From NEMO_3.6 to NEMO_4.0,  from NEMODRAK to DCM
>  by J.-M. Molines  [MEOM/IGE](http://meom-group.github.io/) (Grenoble)
______________________

## **_Overview :_**
 Version 4.0 of NEMO is a major release.  It is the result of a simplification process which was decided by the NEMO consortium
in order to help the maintenance of the code, cleaning obsolete features, revisiting the user interface, adding/replacing
some new numerical schemes, etc... 

 Not surprisingly, the porting from 3.6 to 4.0 requires some work, at different levels.  This document was written 
while performing this porting, together with NEMO4 discovery. After a reminder on DCM concept, it holds basically 
three parts :     
   *  (A)  corresponds to all the changes in the code itself, and the compilation procedure. 
   *  (B) is dedicated to the impacts of the changes at run-time.
   *  (C) Real tests : evaluation of NEMO4

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
 The major change at run time is the introduction of the domain_cfg.nc file (see below).  
 Namelist are also changed but it is not too difficult to figure out the values from a NEMO_3.6 namelist, except for some particular
namelist blocks that will be discussed in more details on a dedicated paragraph.  
 Finally, the xml files used with xios are also completly revisited and restructured. More explanation are given then.


### domain_cfg.nc file

 The major change for users is the  the fact that NEMO4 uses an input file called (`domain_cfg.nc`, name set in the namelist ) which 
hold the equivalent of mesh_hgr and mesh_zgr files. A NEMO tool is provided by the system in order to build this configuration file 
from the classical files used up to NEMO_3.6. (DOMAINcfg).  It is claimed in the documentation that it takes NEMO_3.6-like namelist 
as input, but this is not true (a hack of DOMAINcfg which works fine is provided in DRAKKAR).

It is important to note that the vertical metrics (defining the vertical grid) is now 
read from the domain_cfg.nc file. Thus, surprisingly, NEMO4 does not requires the bathymetry information as input, even if the bathymetry 
is still there in the domain_cfg file ( named `bathy_metry` ? ). 

When using partial steps or sigma coordinates any change in the 
bathymetry (for tunning purpose for instance) imply a re-computation of the domain_cfg file, which can be quite a subtantial piece of work 
for big configurations (although it can be donne in parallel on HPC). This latter point advocate for the use of HPC even when preparing a
new configuration.

In order to compile and use the tool DOMAINcfg, follow the instructions, once in your  NEMO4 CONFIG-CASE directory:

      dcm_mktools -n DOMAINcfg -m <machine> -c <CONFIG-CASE>

 This will compile the DOMAINcfg, including the DRAKKAR modifications. At the end of the compilation you will be told where the tool is available.  
Go there and copy in this directory all the required files from a 3.6 version of the configuration :

     coordinates.nc
     bathy_meter.nc
     namelist_<CONFIG-CASE>  # from a 3.6 configuration

As the tool requires both namelist_ref and namelist_cfg, just make the corresponding links:

     ln -sf namelist_<CONFIG-CASE> namelist_ref
     ln -sf namelist_<CONFIG-CASE> namelist_cfg

Edit namelist_**CONFIG-CASE** in order to add a line in the namcfg block:

     ln_e3_dep   = .true.   ! =T : e3=dk[depth] in discret sens.
                            !      ===>>> will become the only possibility in v4.0
                            ! =F : e3 analytical derivative of depth function
                            !      only there for backward compatibility test with v3.6

(Note that action is always the same as we want to follow the new NEMO4 standards. It is now defined as the default in the DRAKKAR version, so
that the modification above is not mandatory if you stick to ln_e3_dep=T ! ).

 Then you just have to run `make_domain_cfg.exe` (probably in parallel as the code was compiles in the same framework as CONFIG-CASE). You can adjust
the number of subdomains in the nammpp namelist_cfg block. Of course, after a parallel run you need to recombine the pieces of domain_cfg.nc into a single file.
(Using rebuild_nemo works fine).


### Namelists
 As stated in part A of this document, all DRAKKAR changes involving a namelist block, does use a new block, whose name is build from the main NEMO block name appending `_drk` at the end. The the present version of DCM the new blocks are : 

   |  code module|Namelist block| Comments |
   |:--------------|:-------------|:-------|
   |  dommsk.F90:   |  NAMELIST/**namlbc_drk**/ ln_shlat2d, cn_dir, sn_shlat2d | Implentation of 2D shlat |
   | dtatsd.F90:      | NAMELIST/**namtsd_drk**/  ln_tsd_init, ln_tsd_dmp, <br>  cn_dir, sn_tem_ini,sn_sal_ini, sn_tem_dmp, sn_sal_dmp | Differenciation of init and damping fields. |
   | sbcblk.F90:      | NAMELIST/**namsbc_blk_drk**/ ln_kata, sn_kati, sn_katj | Introduction of katabatic winds parametrisations |
   | sbcssr.F90:      | NAMELIST/**namsbc_ssr_drk**/ ln_sssr_flt, ln_sssr_msk, <br> sn_coast, rn_dist, nn_shap_iter | Introduction of coast distance, filtering |
   | trabbl.F90:      | NAMELIST/**nambbl_drk**/ ln_kriteria | Introduction of K-criteria for BBL |
   | tradmp.F90:      | NAMELIST/**namtra_dmp_drk**/ nn_hdmp , nn_file, ln_dmpmask,<br> rn_timsk, cn_dir, sn_dmp | DRAKKAR stype for 3D damping |
   | zdfdrg.F90:      | NAMELIST/**namdrg_top_drk**/ cn_dir, sn_boost | Specifiy information on boost coeficient file for top drag |
   | zdfdrg.F90:      | NAMELIST/**namdrg_bot_drk**/ cn_dir, sn_boost | Specifiy information on boost coeficient file for bottom drag |

 The only new feature with respect to NEMODRAK, is the namelist specification of the files in case of ln_dmpmask, and for the enhanced 
top/bottom friction. It used to be hard coded in NEMO.

#### _Comments:_
 For lateral mixing, diffusivity and viscosity coefficient are **no more** specified as values in the namelist. Instead, these coefficients are computed in the code as:

  aht =  fac * U<sub>d</sub> * L<sub>d</sub><sup>n</sup>, where U<sub>d</sub> is a typical diffusion velocity and 
 L<sub>d</sub> is a diffusion lenght scale, n is 1 for laplacian operator, 3 for bi-laplacian operator, fac is a coefficient
 (1/2 for laplacian operator, 1/12 for bilaplacian operator).

 In the namelist many choices are now possible through the `nn_aht_ijk_t` namelist parameter. According to its value we have the following possibilities :

  | nn_aht_ijk_t |  Formula for aht  | Use rn_Ud, rn_Ld |
  |:------------:|-----------|:------------------------:|
  |   -30        | read in eddy_diffusivity_3D.nc in file | no |
  |   -20        | read in eddy_diffusivity_2D.nc in file | no |
  |    0         | constant value = fac * U<sub>d</sub> * L<sub>d</sub><sup>n</sup> <br>  U<sub>d</sub>  and L<sub>d</sub> from the namelist | rn_Ud  rn_Ld |
  |   10         | aht=F(k) as specified in ldf_c1d |  rn_Ud  rn_Ld |
  |   20         | aht=F(i,j) as specified in ldf_c2d | rn_Ud  |
  |   21         | aht=F(i,j,t)  Treguier et al. JPO 1997 formulation | ? |
  |   30         | aht=F(i,j,k)  ldf_c2d * ldf_c1d | rn_Ud |
  |   31         | aht=F(i,j,k,t) F(local velocity and grid-spacing) | ? | 


 In most of the DRAKKAR run we used a 2D spatial variation of the eddy diffusivity, according to the local gridsize. This is what is 
achieved with `nn_aht_ijk_t = 20 `. In this case, only rn_Ud (U<sub>d</sub>) needs to be specified. In order to use the 
same diffusivity than in previous 3.6 run ( aht_0 in the 3.6 namelist), we may consider the following :

     aht[uv](:,:) = fac*Ud * max(e1,e2)(:,:) ^ n    [ fac=1/2, n=1 ; fac=1/12, n=3]  (4.0)
     aht[uv](:,:) = rn_aht_0 * [ max(e1,e2)(:,:)/ Max(e1,e2) ] ^ n    (3.6)

     ==>  rn_Ud = rn_aht_0/fac/MAX(e1,e2)^n

 In the above formulae, max(e1,e2)(:,:) stands for the local maximum between e1(i,j) and e2(i,j), while Max(e1,e2) stand for the
domain wide maximum for (e1,e2).

    Example for BSAS12 configuration :
       rn_aht_0= 67.52 m2/s (laplacian)
       MAX(e1,e2) = 7012 m 
       ==> rn_Ud = 67.52 * 2 / (7012)^1 = 0.0193  m/s

The advantage of this rigorous way of setting the diffusivity is that we see a link in the coefficents used with laplacian 
and bi-laplacian operators (for a given  U<sub>d</sub>).

       rn_aht0 = rn_Ud * fac * MAX(e1,e2)^n 

### XIOS and xml files
 NEMO4 requires the use of a recent XIOS library (2.5). Although it compiles fine with the 
[trunk](http://forge.ipsl.jussieu.fr/ioserver/svn/XIOS/branchs/xios-2.5) of xios_2.5, some problems have been encountered at run-time
with regard to domain NEMO domain decomposition and land processor elimination. This problem should be solve soon, but in the
mean time, using this particular [branch](http://forge.ipsl.jussieu.fr/ioserver/svn/XIOS/dev/dev_olga) (called 2.5_dev_olga) 
 of XIOS works very well. 

 xml files have been completly reorganized when using NEMO4. At the end, the file  which is used by nemo/xios is still `iodef.xml`,
but it is now a very short file where all the 'ingredients' are sourced:

    <?xml version="1.0"?>
    <simulation>

    <!-- ============================================================================================ -->
    <!-- XIOS context                                                                                 -->
    <!-- ============================================================================================ -->

    <context id="xios" >

      <variable_definition>

          <variable id="info_level"                type="int">10</variable>
          <variable id="using_server"              type="bool">true</variable>
          <variable id="using_oasis"               type="bool">false</variable>
          <variable id="oasis_codes_id"            type="string" >oceanx</variable>

      </variable_definition>
    </context>

    <!-- ============================================================================================ -->
    <!-- NEMO  CONTEXT add and suppress the components you need                                       -->
    <!-- ============================================================================================ -->

    <context id="nemo" src="./context_nemo.xml"/>       <!--  NEMO       -->

    </simulation>


So, all the nemo dependent information is within the `context_nemo.xml` file, which in turn is also very short :


     <!--
      ============================================================================================== 
           NEMO context
     ============================================================================================== 
     -->
     <context id="nemo">
     <!-- $id$ -->
     <!-- Fields definition -->
         <field_definition src="./field_def_nemo-oce.xml"/>    <!--  NEMO ocean dynamics     -->
         <field_definition src="./field_def_nemo-ice.xml"/>    <!--  NEMO sea-ice model      -->
         <field_definition src="./field_def_nemo-pisces.xml"/> --> <!--  NEMO ocean biology      -->
     
     <!-- Files definition -->
         <file_definition src="./file_def_nemo-oce.xml"/>     <!--  NEMO ocean dynamics      -->
         <file_definition src="./file_def_nemo-ice.xml"/>     <!--  NEMO sea-ice model       -->
         <file_definition src="./file_def_nemo-pisces.xml"/>  <!--  NEMO ocean biology       -->
         <!-- 
     ============================================================================================================
     = grid definition = = DO NOT CHANGE =
     ============================================================================================================
         -->
     
         <axis_definition>
           <axis id="deptht" long_name="Vertical T levels" unit="m" positive="down" />
           <axis id="depthu" long_name="Vertical U levels" unit="m" positive="down" />
           <axis id="depthv" long_name="Vertical V levels" unit="m" positive="down" />
           <axis id="depthw" long_name="Vertical W levels" unit="m" positive="down" />
           <axis id="profsed" long_name="Vertical S levels" unit="cm" positive="down" />
           <axis id="nfloat" long_name="Float number"      unit="-"                 />
           <axis id="icbcla"  long_name="Iceberg class"      unit="1"               />
           <axis id="ncatice" long_name="Ice category"       unit="1"               />
           <axis id="iax_20C" long_name="20 degC isotherm"   unit="degC"            />
           <axis id="iax_28C" long_name="28 degC isotherm"   unit="degC"            />
         </axis_definition>
     
         <domain_definition src="./domain_def_nemo.xml"/>
     
         <grid_definition src="./grid_def_nemo.xml"/>
     
     </context> 


In this context file, all the field definitions are sourced, as well as the file definitions for different NEMO
components : oce, ice, pisces.  In most of the cases, the end-user only have to change file_def_*_xml, to fit his
proper willings. 

When using DCM, we recommend to direct the output in a specific directory `$WORKDIR/<CONFIG>-<CASE>-XIOS.<seg>`. To do so,
we change  the `<file_definition>` line to 

         <file_definition type="multiple_file" name="<OUTDIR>/@expname@_@freq@" sync_freq="1d" min_digits="4">

The `<OUTDIR>` word is a keyword that the running script (`nemo4.sh`) will automatically update to the correct value.   
Also, in order to be fully compliant with  DCM specifications regarding the output files, each file may have some specific
global attributes. Example :

     // global attributes:
                ......
		:start_date = 19920101 ;
		:output_frequency = "1d" ;
		:CONFIG = "BSAS12" ;
		:CASE = "MJM151" ;
                .....

This is achieved by adding global attribute specification at the end of each file defined in the file_...xml, as chown
on the example ( it corresponds to the `<variable>` settings ):

     <!-- T FILES -->
        <file id="file6" name_suffix="_gridT_" description="ocean T grid variables" >
            <field field_ref="toce"         name="votemper"   />
            <field field_ref="soce"         name="vosaline"   />
            <field field_ref="ssh"          name="sossheig"   />

           <variable name="start_date"       type="int"><NDATE0>    </variable>
           <variable name="output_frequency" type="string">1d       </variable>
           <variable name="CONFIG"           type="string"><CONFIG> </variable>
           <variable name="CASE"             type="string"><CASE>   </variable>
        </file>

Here also `<NDATE0>`, `<CONFIG>` and `<CASE>` are keywords that the  running script (`nemo4.sh`) will automatically update.

## **_(C) Real tests, evaluation of NEMO4_**
###   Testing with [BSAS12](Porting_BSAS12-test.md)
###   Testing with [ORCA05](Porting_ORCA05-test.md)
###   Testing with [ORCA025](Porting_ORCA025-test.md)
###   Testing with [ORCA12](Porting_ORCA12-test.md)
