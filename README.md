# DCM
Drakkar Config Manager

## Purpose, history :
Drakkar Config Manager (DCM) is an environment tool with the main goal of helping the deployment a model configuration, based on the [NEMO](https://www.nemo-ocean.eu/) modeling system. This tool uses different layers of the code: the deepest one being the official NEMO release (NEMOREF release), the intermediate one (or DRAKKAR layer)  and the upper one being the user layer. All layers share the same code layout and building the code for a particular configuration consist of copying hierarchically (from the deepest to the uppermost layer), so that at the end of this code building process, the DRAKKAR modules and the user modules (differing from  NEMOREF) replace the original module.  A series of comparing tools are also provided for assessing the differences between layers. 

This tool has been initiated by our modeling group, in the mid 2000's in the frame of the [DRAKKAR project](https://www.drakkar-ocean.eu/). Since then, the tool followed the evolution of NEMO, and was maintained in a specific repository of the [LEGI forge](https://servforge.legi.grenoble-inp.fr/projects/DCM). One strength of this system is the ability to trace the source code, with reference to both the revision number of the NEMOREF code as well as the revision number of the DRAKKAR code, both saved into the user configuration directory. 

In July 2018, the anouncement of the major release NEMO_4.0 with deep modifications on the structure of the code implied also deep modifications in the DCM tool.  We decided to switch from the trac/svn forge to GitHub, as a new project.  The original forge will be maintained for tracability of previous DCM based configurations. 

9th, February, 2021 : v4.0.2 tag created, corrsponding to NEMO4 release r4.0.2  
9th, February, 2021 : HEAD (or master)  is in phase with NEMO4 release r4.0.5
27th, March, 2021 : HEAD (or master)  is in phase with NEMO4 release r4.0.6



## License :
   DCM is distributed under the CeCILL [license](License/DCMCeCILL.md).

## Using DCM :
  Users are invited to read the [user manual](DOC/dcm_user_manual.md) to learn how to use DCM !

