# RUNTOOLS user manual
## Forewords:

  If you read this manual is that you probably have successfully installed and compiled a NEMO configuration, using DCM ! (to make it short you have a `nemo4.exe` in the `$PDIR/RUN_<CONFIG>/<CONFIG>-<CASE>/EXE/` directory, corresponding to your settings (code and cpp keys)). Ideally you are now only missing some runtime files to make your first run (and independently of using DCM or not ):  
 * the **namelist** files (for ocean, ice, passive tracers etc...)
    - Since NEMO_3.6, NEMO uses 2 levels of namelists: `namelist_ref` where all variables are defined with default values (in general corresponding to ORCA2 setting), and `namelist_cfg` which only re-define values changed with respect to the reference. Although it make sense, we found that having those 2 levels files (though never have the full set of parameters in a single file) is not convenient and probably a source of error. In DCM, we choose to maintain a single namelist (having all the variables in it, like `namelist_ref`) with the right parameters for the working configuration. This choice is transparent to NEMO as we just copy the *ad hoc* `namelist_ref` to  `namelist_cfg`. 
    - Reference namelist (usefull to know about all the parameters) can be found in the `cfgs/SHARED` sub directory.
    - Namelists are a very clever and comfortable means to set parameters in NEMO, but it is really critical to double check the namelist parameters, as it is quite easy to make mistakes !
    - When using DRAKKAR modifications of the code, and in case a namelist block is envolved, DRAKKAR code redefine an additional namelist block for its own parameters. The DRAKKAR namelist block has the same name as the standard one, but with `_drk` appended (*e.g* `&namlbc` (std) and `&namlbc_drk` (DRAKKAR parameters). Doing so, a namelist setup can be run with standard NEMO code as well as with DRAKKAR modificated code.
    - Refer to the [NEMO book](https://www.nemo-ocean.eu/wp-content/uploads/NEMO_book.pdf), for detailed description of the namelist parameters (**caution: link point to NEMO_3.6 !!**)
 * the **xml** files used by XIOS for model output.
    - These files allow the link between NEMO and XIOS. At the end XIOS just need a file called `iodef.xml` where all the informations on both NEMO XIOS implementation and user requirements for output are described. Recent versions of XIOS (> 2) allow the use of separate 'sub' files which are used as sources in `iodef.xml`, which really improves the readibilty of the files.
    - In NEMO4 implementation, the xml sub files are :
      * `domain_def_nemo.xml`
      * `axis_def_nemo.xml`
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
    - The `domain_cfg.nc` configuration file is a new feature in NEMO4. It replaces the `coordinates.nc` and `bathy_meter.nc` files and includes the definition of the vertical grid
(no more initialized within NEMO, as it was the case in previous versions). A specific tool (`DOMAIN_cfg`) is dedicated to the building of the `domain_cfg.nc` file from files and namelist used 
in NEMO-3.6. See the [appendix](#domain_cfg-tool-) for instruction how to prepare the `domain_cfg.nc` file.


 In order to run NEMO, the principle is to have all the runtime files in a common directory, together with `nemo4.exe` and `xios_server.exe` and issue (for instance -- there are variant according to the HPC system--) a command like :

```
   mpirun -np $XIOS_CORES xios_server.exe : -np $NEMO4_CORES nemo4.exe 
```

  Where `XIOS_CORES` is the number of cores dedicated to XIOS, and `NEMO4_CORES` is the number of cores dedicated to NEMO4. This statement launch 2 executables with the MPMD (Multiple Program Multiple Data) paradigm. 

  This general overview of how to run a NEMO configuration is valid even if you do not use DCM. DCM's runtools offer an environment for automatization of most of the task required for run production.  It is  a collection of scripts, handling all the machinery required to produce a long simulation (chaining elementary segments of run), dealing with the model output, restart files etc...  Although it works even for simple configurations (such as test cases or idealized cases), it is primarily designed for complex realistic cases (which explains the relative complexity of the tools).

## Starting with DCM's runtools:
  When you did `dcm_mkconfdir_local` for preparing your configuration, both `EXE` and `CTL` directories were created in `$PDIR/RUN_<CONFIG>/<CONFIG>-<CASE>`. In `EXE` you already know that you have `nemo4.exe` and `CPP.keys` for your configuration. All the management of the run will be done from `CTL` (standing for 'control'). 

 1. **Prepare the RUNTOOLS** ( to be done once ):

    The main production script [`nemo4.sh`](../RUNTOOLS/lib/nemo4.sh), is valid for any system. The portability through different systems is achieved by using bash functions instead of machine dependent command. Let take an example to make it clear: for instance, the command used to submit a job on a batch system depends on the scheduler you are using; it may be `qsub`, `llsubmit`, `sbatch` .... In the main script a generic function `submit` is used in place of all the variant. Then `submit` is defined within a function file depending on the machine. Before using the main script, this specific function file is sourced so that the *ad hoc* command for submission will be used. 

    In order to be fully generic, the name of the functions file is hard coded as `function_4.sh`. Hence, preparing the RUNTOOLS is limited to make a link pointing to the *ad hoc* functions file for your machine:
    
    ```
    cd $RUNTOOLS/lib
    ln -sf function-<MACHINE>.sh function_4.sh
    ```

    In the actual DCM, functions for HPC machines **ada**, **irene** and **occigen** are provided. 

    For your purpose you may need to create the function file for your machine. Starting from one of the existing file is a good option. Most of the tricks are in `submit`, `runcode` and other job control statements. You can then contribute to the RUNTOOLS by sending your own `function_<MACHINE>.sh` file !

 1. **Create your run time environment** in `CTL`

    As already mentionned, all the management of the runs must be done in the corresponding `CTL` directory. Hence, the very first step is to populate it with relevant files. This can be done easily from the *confcase* directory (`$UDIR/CONFIF_<CONFIG>/<CONFIG>-<CASE>/`) with the single command:

    ```
    make ctl
    ```

    This will copy template files to your `CTL`, with some automatic editing to fit the CONFIG and CASE names :
    * `includefile.sh`
    * `<CONFIG>-<CASE>_<MACHINE>.sh`
    * `run_nemo.sh`
    * `<CONFIG>-<CASE>.db`
    * `NAMELIST/namelist*`
    * `XML/*.xml`

 1. Edit template files :
    * **`includefile.sh`**: This file is essentially used for setting up the path of some key directories, through environment variables. These variables have a name starting with F_ for directories on the archiving system (`SDIR`), and with P_ for the directories on the production system. At present, it is recommended to work mainly with the production machine (because the archiving system on HPC, is likely not a mirror of the production machine, as it used to be). Nevertheless, up to now we maintain this double path system: The running script uses special functions for retrieving files on the `TMPDIR` area, looking first on the production machine, and then on the archiving system if not found.  Most of the default setting may work if you follow a strict building of DCM. However, you still have some freedom to chose paths for forcing files, weight files, XIOS lib etc... Comments in the script are hopefully clear enough to describe the different environment variables.

      Additionally, it allows the correspondance between hard coded file names within NEMO and the corresponding file in the real world. Actually, only a few hard coded names remains in NEMO4. Most of the file names are now passed to NEMO through a namelist variable. In the `includefile.sh` on a same line you have the world name on the left and the corresponding NEMO name on the right.  Script comments may help for understanding who is who...

      Last but not least, the variable `MAXSUB` defined at the end of the file is there for automatic re-submission: Segments of job will be chained until `MAXSUB` segments are performed. ( See `<CONFIG>-<CASE>.db` below).
    * **`<CONFIG>-<CASE>_<MACHINE>.sh`** : This is the 'job' file that will be submitted to the HPC batch system. Therefore it starts with a header corresponding to the batch system directives. This is why it is 'machine' dependent. 

      The remaining part of this script is the same on any system. It needs editing for the following:

      ```
      # Following numbers must be consistant with the header of this job
        export NB_NPROC=3150    # number of cores used for NEMO
        export NB_NPROC_IOS=34  # number of cores used for xios (number of xios_server.exe)
        export NB_NCORE_DP=0    # activate depopulated core computation for XIOS. If not 0, RUN_DP is
                                # the number of cores used by XIOS on each exclusive node.
      # Rebuild process 
        export MERGE=0          # 0 = on the fly rebuild, 1 = dedicated job
        export NB_NPROC_MER=231 # number of cores used for rebuild on the fly  (1/node is a good choice)
        export NB_NNODE_MER=20  # number of nodes used for rebuild in dedicated job (MERGE=0). One instance of rebuild per node will be used.
        export WALL_CLK_MER=3:00:00   # wall clock time for batch rebuild
      ```

      Although the comments are self explanatory, some additional information on the rebuild process may be of interest ! In DRAKKAR, we advocate for using `XIOS` in detached mode (*i.e.* using the `xios_server.exe`) and with the `multiple_file` protocol (in oposition to `one_file` protocol), hence requiring a rebuild process. The reasons for this choice are two-fold: The rebuild process is used (i) to create netcdf4/hdf5 files with chunking and deflation (file size can be divided by 3 or 4 !) and (ii) to 'fill the holes' in `nav_lon nav_lat` corresponding to eliminated land processors in the NEMO computation; this latter point avoid having incohent mapping variables (on land).  Depending on the size of the domain, on the available computer resources etc... it may be of interest to launch the rebuild 'on the fly', meaning in the same job than NEMO, or to launch a dedicated job for the rebuild process. The choice is made by the value of the `MERGE` variables, and the requested number of core for the rebuild process must be defined, setting up the *ad hoc* variables. In both cases, DCM uses a DRAKKAR made rebuild tool (`REBUILD_MPP tool`) which work in parallel. More information about this tool, in the appendix of this manual.
    * **`run_nemo.sh`** : This script requires very few editing if any. It is the model launcher that submit the job. It is usefull to have it at hand in CTL, as in some cases or tests, you may want to add extra statements in the workflow. (Such as copying a non standard file to the `TMPDIR` for instance). Once all is correctly setup, `./run_nemo.sh` is the command that you will perform in order to start a production segment.
    * **`<CONFIG>-<CASE>.db`** : This file is automatically updated by the main running script. It holds information about the initial and last step of a job segment, as well as the date of the end of a segment. It just requires initialization (the first line of the file). In the template there are 3 values : Segment number, first model step, last model step of the segment. It is likely that the last step of the segment must be adjusted (depending on the model time step and length of the segment to perform). 
    * **`namelists`** : The template namelists are copied in CTL/NAMELIST. Remember that with DCM, we manage only a full namelist, *i.e.* a namelist similar to namelist_ref, where all the parameters correspond to the actual `<CONFIG>-<CASE>`. Once carefully checked, the operational namelists (ocean, ice, top etc... ) should be copied or moved to CTL/, with the name `namelist_<CONFIG>-<CASE>`. The main running script will duplicate it into namelist_ref and namelist_cfg for NEMO compliance.

      DRAKKAR code possibly uses extra variables in the namelists. In this case, a new namelist block is created, instead of modifying a standard namelist block. The DRAKKAR related block use the standard name with `_drk` appended. When using DCM, some namelist variables are updated during the run process. They appear in the namelists as `<KEY_WORD>`. Keep them carefully when editing the template namelist. (*i.e.* `<NN_NO>`,`<NIT000>`,`<NITEND>` etc...; for these particular 3 variables, the running script use the `<CONFIG>-<CASE>.db` information to set them correctly). 
    * **`xml files`** : Templates xml files are located in CTL/XML directory. The ones you need should be copied or moved to CTL/.  As already explained in the forewords of this manual, XIOS only needs `iodef.xml` file. However, all the other xml files are included at different levels, as shown in the following tree:

      ```
      iodef.xml                                            __|  check using_server : true
         |                                                 __
         |___ context_nemo.xml                             __|  only keep the needed file_def !
                          |                                __
                          |_______ field_def_nemo-oce.xml    |
                          |_______ field_def_nemo-ice.xml    | (no edit required for
                          |_______ field_def_nemo-pisces.xml |     standard use)
                          |_______ field_def_nemo-.....xml __|
                          |                                 
                          |                                __
                          |_______ file_def_nemo-oce.xml     |
                          |_______ file_def_nemo-ice.xml     |  NEED editing for your choices
                          |_______ file_def_nemo-pisces.xml  |  
                          |_______ file_def_nemo-.....xml  __|
                          |
                          |
                          |_______ axis_def_nemo.xml       __|  no edit required for std use (used from 4.0.2 onward)
                          |                                __
                          |_______ domain_def_nemo.xml     __|  no edit required for std use
                          |
                          |                                __
                          |_______ grid_def_nemo.xml       __|  no edit required for std use
      ```

    As shown on the graph above, major editing is only required in the `file_def_nemo-*.xml`. And in `context_nemo.xml` you should keep only the `file_def_nemo-*.xml` required by your setting (*e.g.* : if you do not use the ice model, you must not have `file_def_nemo-ice.xml` sourced in `context_nemo.xml`).

      Note the `file_def*xml` files used in DRAKKAR implements some changes with respect to the standard NEMO:
      *  `<file_definition\>` xml tag  is changed to (for example):

           ```
           <file_definition type="multiple_file" name="<OUTDIR>/@expname@_@freq@" sync_freq="1d" min_digits="4">
           ```

           `<OUTDIR>` is a `<KEY_WORD>` that will be changed to `$DDIR/<CONFIG>-<CASE>-XIOS.<seg>` in order to have the XIOS output in a separate directory for each production segment.

      * Global attributes in the output netcdf files are used, having for each file definition:

          ```
             <variable name="start_date"       type="int"><NDATE0>    </variable>
             <variable name="output_frequency" type="string">1h       </variable>
             <variable name="CONFIG"           type="string"><CONFIG> </variable>
             <variable name="CASE"             type="string"><CASE>   </variable>
          ```
          Here again, the `<KEY_WORD>` will be replaced by their value at run-time. These global attributes are used in the [DRAKKAR monitoring tools (DMONTOOLS)](https://github.com/meom-group/DMONTOOLS).

      `domain_def_nemo.xml` requires editing only if you add new sub-domains for output.  
      `field_def_nemo-*.xml` requires editing if you want to add extra variables in the already long list of possible variables to output. New variables should have a corresponding call to iom_put in the NEMO code.  

 1. Run the code

    ```
    ./run_nemo.sh  
    ```

    And that's it ! 

    In order to help the management of a run, small tools were developped. They are now collected in the DCM toolkit. See the [DCM toolkit manual](./dcm_toolkit.md) for description and use of the toolkit.
 1. Post processing the output : see [dedicated manual](./dcm_post_process.md) describing hints for post processing

## Cloning an existing configuration:

It happens very often that a new configuration is built to make sensitivity experiments with respect to some references.  For this particular case, DCM offers a very easy procedure to clone the CTL of an existing running configuratin (as it exists for the code itself). 
> Note this procedure is not limited to sensitivity experiment, you can clone a completly different configuration, but of course then you need to adjust the parameters in the namelists and in the xml files.

```
   cd <new empty CTL>
   dcm_clone_ctl -c <CONFIG>-<CASE>
```

With this command you will populate the new empty CTL with a valid set of files identical to the ones for \<CONFIG\>-\<CASE\>, but the the correct names.  Then you need to adjust the namelists, and possibly the xml files. 
> Note that if the CTL where you want to clone is not empty, no cloning will be done (in order to preserve possibly important settings!).

## Setting up a scalability experiment.
It is always very usefull to perform a scalability experiment before running a long simulation, in order to optimize the computer ressources. It 
is also required when you prepare a proposal to computing centers. 

With DCM, when using DRAKKAR modified code, we know the time (up to ms) when each step starts. This allow to easily draw a graph of the progression of the run (number of steps
performed vs time). In general, if all works smoothly, this graph is almost a perfect straight line, which slope (in stp/min, for instance) is an indicator of the
code performance.  On the other hand, this progression graph is also showing hardware problems very efficiently when they occur.  Due to the very linear shape of the progression graph,
only a short experiment (says 100 steps) is able to tell us the performance of the code for a given domain decomposition, hence number of cores.  Note that when evalutating the slope (stp/mn) we disregard the first and last 10 steps which are slow due to one-time initialisation or to closing files. The slope is evaluated with linear regression (specific dcmtk tool). 

  In order to make it easy, a special procedure is proposed to set up a scalability experiment.
 1. Create a new configuration  
     It is likely a good idea to use a CASE name refering somehow to scalability, but nothing mandatory.

    ```
    dcm_mkconfdir_local <CONFIG>-<CASE>
    ```

 1. Clone some differences for in the code (eliminating writing of restart files) and compile:  
    In `$HOMEDCM/DRAKKAR/NEMO4/cfgs`, `CONFIG-CASE.scal` configuration directory hold some 
    modified source code. You may use this config as a `PREV_CONFIG` in `makefile` in order to import
    those sources.

    ```
    cd $UDIR/CONFIG_<CONFIG>-<CASE>
    #  edit makefile there by setting:
    PREV_CONFIG=$(HOMEDCM)/DRAKKAR/NEMO4/cfgs/CONFIG-CASE.scal
    # then
    make copyconfig
    # check the makefile for `MACHINE` and other classical points ...
    # check the CPP.key file for the keys you need.
    make install && make
    ```

    > You may have some warning or error message when doing `make copyconfig`. You can ignore them safely.    
    > At this point the NEMO code is ready for use.

 1. Prepare the CTL directory  
    Preparing the CTL directory for scalability experiment is done by:

    ```
    cd  $PDIR/RUN_<CONFIG>/<CONFIG>-<CASE>/CTL/
    # use -s option in mkctl for scalability experiments
    dcm_mkctl -m occigen -s -a 
    ```

    Once this is done, there are templates for namelist, includefile and main running script. You still need to have a look at those file in order to fix some PATH (includefile), or some parameters (namelist). You also need to set up a set of xml files for the output (despite the fact that output will not be processed
in scalability experiment). You also need to setup the `<CONFIG>-<CASE>.db` in order to set the number of steps you want to perform foreach experiment. (A few hundreds is OK).
 1. Determine the domain decomposition you will test (using MPP_PREP tool)  
    This requires the instalation of the `MPP_PREP` tool

    ```
    cd $UDIR/CONFIG_<CONFIG>-<CASE>
    dcm_mktools -n MPP_PREP -m <MACHINE> -c <CONFIG>-<CASE>
    ```
    > Now the `mpp_optimize.exe` program is ready for use in `$WORKDIR/W<CONFIG>-<CASE>/tools/MPP_PREP`

    ```
    cd $WORKDIR/W<CONFIG>-<CASE>/tools/MPP_PREP
    # edit the namelist in order to 
    #   set cn_fbathy=<CONFIG>-<CASE>_domain_cfg.nc. Variable name is bottom_level by default.
    #   set the maximum domains you are looking for (nn_procmax)
    # Run the code (mpp_optimize.exe -h gives you USAGE information).
    ./mpp_optimize.exe -n namelist 
    ```

    > Now you have a `processor.layout` file in the `MPP_PREP` directory. You can find optimal domain decomposition by using `./screen.sh  <CORES>`. But for the scalability experiment you can even prepare a script that will launch all possible cases... 

    ```
    dcmtk_scal_prep -h
    USAGE : dcmtk_scal_prep [-h] [-l layout_file] [-m minproc] [-M maxproc] [-s proc_step]
  
      PURPOSE:
       This tool is part of the process for setting up a scalability experiment.
       It assumes that you have already run the MPP_PREP tool for defining all the
       file, narrowing the results to a specific range of cores. In addition, it will
       call a python script that produces a job file to be submitted for all the 
       decomposition to test. At this point, the python tools assumes that the target
       machine is occigen (SLURM batch system). Hence, the job file may require some
       editing (batch header part) to fit your HPC system.
       Due to the interaction with MPP_PREP it is likely a good choice to run this tool
       where you have the MPP_PREP results, in particular with the screen.sh script.
  
      OPTIONS:
       -h : Display this help message.
       -l layout_file : Pass the name of the layout file produced by MPP_PREP  procedure.
             Default is processor.layout
       -m minproc : Gives the minimum number of core you want to test. Default=50 
       -M maxproc : Gives the maximum number of core you want to test. Default=600
       -s proc_step : Gives the core step for scaning the core range.  Default=50
    ```

    So doing (for instance) :

    ```
    module load python  # the dcmtk script launch a python program...
    dcmtk_scal_prep -l processor.layout -m 100 -M 2000 -s 50 
    ```

    will produce `run.sh` metascript which is the driver for the scalability experiment.

    ```
    chmod 755 run.sh
    cp run.sh $PDIR/RUN_<CONFIG>/<CONFIG>-<CASE>/CTL/
    cd $PDIR/RUN_<CONFIG>/<CONFIG>-<CASE>/CTL/
    ```

    > You are now ready with the domain decomposition to test.
 1. Run the metascript
    The metascript `run.sh` prepared during the `MPP_PREP` phase is now ready to be run.

    ```
    ./run.sh
    ```
    This metascript launch quite a lot of jobs, each differing by the number of cores they are testing. When the jobs are all completed, the raw scaling results are with the `nemo_occigen*.o*` files (job output). 
 1. Analyse the results.  
    You may be able to plot a very primary scalability graph with the following commands (`dcmtk_scal_plot`):

    ```
    for f in nemo_occigen_*.o* ; do echo -n $f " : "  ; dcmtk_rate -s 180 -b 1 -f $f; done > perf.log
    cat perf.log |sort -t_ -k4n  | sed -e "s/\.o/\ \.o/g" -e "s/_/ /g" | awk '{print $5 " " $11}' | graph -TX -S 2 -m -1
    ```

    > Scalability experiments are an iterative process where you need to tune some parameters, such as the number of xios server you are requiring etc... This first and rather quick view of the scalability performance for a particular configuration helps you to choose an optimal number of cores. Then you may refine the experiment around this sweet point, in order to achieve the best performance before going to a production run.

## Appendix
### DOMAIN_cfg tool:
 This tool (proposed by NEMO ST) is aiming at creating the required `domain_cfg.nc` file from information used in NEMO3.6 ( namelist, coordinates.nc and bathy_meter.nc). In DCM, some fixes where added in order to be fully compliant with DRAKKAR namelist. In addition (DCM goodies!), `dcmtk_dom_doc.exe` can be used to encode the namelist_cfg and names of bathymetry and coordinates used when building the domain_cfg.nc file, in order to have a full tracability of the domain_cfg file (in particular, the namdom coefficient for the vertical grid are important to record if the same vertical grid need to be reproduced.)  
#### Compiling DOMAIN_cfg tool:
 In order to compile the tool, the best way in the following:

```
    cd $UDIR/CONFIG_<CONFIG>-<CASE>
    dcm_mktools -n DOMAIN_cfg -m <ARCH> -c <CONFIG>-<CASE>
```

where `<ARCH>` is the same used for compiling NEMO (see MACHINE in the makefile *e.g.* X64_OCCIGEN2 ). At the end of the compilation, a message indicates where you can find the executables. Go there and copy/move all the .exe files and template namelists to the directory where you build the domain_cfg.nc file.

#### Using make_domain_cfg.exe

#### Using dcmtk_dom_doc.exe 


### REBUILD_MPP tool:
  In order to compile this rebuild tool (allowing parallel rebuild), the best way is  the following:

```
    cd $UDIR/CONFIG_<CONFIG>-<CASE>
    dcm_mktools -n REBUILD_MPP -m <ARCH> -c <CONFIG>-<CASE>
```

where `<ARCH>` is the same used for compiling NEMO (see MACHINE in the makefile *e.g.* X64_OCCIGEN2 ). At the end of the compilation, a message indicates where you can find the executables. Go there and copy/move all the .exe files into `$WORKDIR/bin` (*i.e.* a directory in your PATH). There are many executables according to the scheme of parallelization, but they have all the same user interface, and in the DCM running script, `mergefile_mpp4.exe` is used as standard (and user can use it blindly !). Just for information, usage of mergefile series is :

   ```
  usage :  mergefile4 -F | -f <file*>_0000.nc  [-c coordinate_file ] ...
           [-d output_directory] [-r ] [-nc3] [-v] [-kmax K-max] ...
           [-w imin imax jmin jmax] [-b BLOCK-number] [-encoding ndigit]
       
      PURPOSE :
          This program recombine splitted netcd files into a single global file.
       
      ARGUMENTS :
        Use either one of the following ( not both) :
         -F : takes all possible _0000.nc files in current directory
            This option has been introduced to overrid shell limitations 
            regarding the length of the arguments.
        or :
         -f <file*>_0000.nc : list of full name of rank 0 <file*> to be rebuild
  
        If the number of digit used for coding the rank number in the file
        names is not 4, condider the use of -encoding option.
       
      OPTIONS :
         -c coordinates file : Use coordinate file to patch nav_lon,nav_lat
               where there are land processors
         -d output_dir       : Specify the output directory
         -r  : Rename output file following DRAKKAR rules. Imply that the output
               directory for each file is <freq>_OUTPUT (<freq> is 1d, 5d ... )
               (-r option supersedes the -d option.) 
         -nc3: Create netcdf3 files instead of default netcdf4 
         -v  : Verbose, more informations on output
         -w imin imax jmin jmax : Only rebuild the sub-domain that are within 
               the specified windows (in model coordinates).  This is usefull
               when rebuilding HUGE restart file having in mind to use only
               a sub zone of the total domain. Using this option reduces the
               time of the rebuild. However, the resulting file is still on
               the total domain and will need some croping afterward.
               Default behaviour is to rebuild the total domain.
         -kmax K-max : define the maximum number of level to be processed in
               the merged output file. (Usefull with -w option).
         -b BLOCK-number : This option was written in order to reduce the 
               number of files opened simultaneously during the rebuild.
               The total number of (mpp) files will be processed by blocks.
               BLOCK-number is the number of blocks you choose. Default is 1.
               This option can considerably reduce the memory imprint, linked
               with the massive opening of all the mpp files together.
               As a rule of thumb, BLOCK-number can be choosen in order to
               block size ranging from 1000 to 2000. It of course depends on
               the available memory.
         -encoding ndigit : indicates the number of digits used for coding
               the rank number in the mpp file's names. Default is 4
       
      REQUIRED FILES :
         none
       
      OUTPUT : 
        netcdf files : <file>_merg.nc, unless -r option is used.
        variables    : same than in <file>_0000.nc file
       
      SEE ALSO :
       rebuild_nemo, splitfile
   ```
