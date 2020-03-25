# Getting Started

## Pre-requisite :
  You need to have an operational XIOS library and `xios_server.exe` executable.  
You can refer to [this file](../DCMTOOLS/NEMOREF/xios_revision.md) for recomended XIOS revision for this DCM release.

## Road-map:
  Assume that your git repositories are located in some `$DEVGIT` (up to you!)
  1. Clone the DCM GitHub repository in `$DEVGIT`

      ```
      cd $DEVGIT
      git clone https://github.com/meom-group/DCM.git DCM_4.0
      ```

  1. If you want to use a specific DCM tag, please read carefully the [tag](Using_tags.md) document before proceeding.

  1. Get the official NEMO code from NEMO-ST IPSL server.

      ```
      cd $DEVGIT/DCM/DCMTOOLS/NEMOREF/      
      ./getnemoref
      ```

  1. Set up DCM module

      It is highly recommended to set DCM as a module, in particular if you are working with different release of NEMO at the same time. Likely, `module` is available on any HPC center, and is very easy to install on a linux system. However, you still have the possibility to set up the environment variables in your `.bashrc` file (see below _Configuring your environment_ paragraph).  
   The steps for setting up DCM modules are:  
      **a)** create a `modules` directory in your `HOME` (for instance).  
      **b)** append this directory name to `MODULEPATH` in your `.bashrc` file:  

        ```
        export MODULEPATH=$MODULEPATH:$HOME/modules/
        ```

      **c)** create  there a DCM directory:  

        ```
        mkdir $HOME/modules/DCM
        ```

      **d)** copy the [template module file](../DCMTOOLS/templates/module_example) into `$HOME/modules/DCM/`  

        ```
        cp $DEVGIT/DCM/DCMTOOLS/templates/module_example $HOME/modules/DCM/4.0
        ```

      **e)** edit `$HOME/modules/DCM/4.0` to fit your settings. Basically only few lines have to be modified :  
        ```
        set             version         4.0
        set             alter_version   4.0
        set             base_path       $::env(DEVGIT)/DCM_$version/DCMTOOLS
        set             alter_path      $::env(DEVGIT)/DCM_$alter_version/DCMTOOLS
        ```

        > Note that _version_ refer to the actual version you are working with (DCM_4.0 in this example).  
        > _alter_version_ points to another DCM version, and is used in some comparison tools (see advanced usage)  
        > For first time use, no problems to set _version_ and _alter_version_ to the same thing.  
        > IMPORTANT : version must be coherent with the name you give to your DCM repository in ```$DEVGIT```  


      Once these steps done, you can load the DCM environment by using :

       ```
       module load DCM/4.0
       ```

  1. Configuring your environment for DCM

      If you set up the DCM module, you already defined some environment variables directly linked with the DCM version you are using.

      You still have to define some environment variables which depends on your own working policy, defining important path for the DCM to work.   
      * `WORKDIR`: On some HPC system, this environment variable is set by the system. It points to a large and fast access file system, where data files, model output can be saved temporarly.  
      * `UDIR` : This is the place where the code of all configurations created with DCM will be saved. It does not require lot of space, but contains precious code. It is advisable to put it on a file system with frequent backup. Default is `$HOME/CONFIGS/`  
      * `PDIR`: This is the place where the running scripts will be copied, from where the jobs will be launched, and where the job's return files will be written. It is advisable to put it on a file system with frequent backup. Default is `$HOME/RUNS/`  
      * `CDIR`: This is where the compilation system will write the temporary files, the produced executable etc... Nothing really precious in there. If lost, it can be reconstructed by compiling the code again. Default is `$WORKDIR`  
      * `SDIR`: This is the root of the archiving system. (Permanent storage). In general it is located on very large file system, which common name is `$STOREDIR`. This is where the long term archiving of model results will take place. In most of the HPC centers, this `STOREDIR` is associated with automatic backup on tapes, and there are often restrictions on the number of files (inodes) to put on this space. It is much more efficient to put 'big' files on it (tar files or concatenated files, whose size is of the order of 100 Gb). Although this limitation looks like a strong constraint, it really helps when you need to get back a large amount of data from the tapes. Using big files, can reduce the required elapsed time of retrieval by a factor of 10 or 50 !  
      * `DDIR` : At some HPC center (*e.g.* TGCC/irene), there are 2 levels of working file system. One called WORKDIR with a relatively small quota (2Tb), but permanent (no automatic removal), and another one (*e.g.* `/scratch`) where the quota can be very high (~20 to 100 Tb), but with automatic removal of old files. In this case, `DDIR` point to this file sytem. The default is `$WORKDIR`.  

      A template file is provided ([dcm_setup_module.sh](../DCMTOOLS/templates/dcm_setup_module.sh)) to be appended to your `.bashrc` file or equivalent. Note that if you are not using the DCM modules, another template file is also provided ([dcm_setup_nomodule.sh](../DCMTOOLS/templates/dcm_setup_nomodule.sh)) for setting both DCM related environnement variables (normally set in the module), and the 4 environment variables already mentioned.

      The template files contains comments that describes the variables. Please refer to these comments for your own case.

## Next :

 **After completing these step you are ready to use DCM**  
  Please read the corresponding manual for [compilation and code maintenance](dcm_compil_manual.md) or for [run time scripts](dcm_rt_manual.md). 
 



