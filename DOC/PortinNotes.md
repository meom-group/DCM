# Porting notes : from NEMO_3.6 to NEMO_4.0, from NEMODRAK to DCM

## Overview :
 Version 4.0 of NEMO is a major release.  It is the result of a simplification process which was decided by the NEMO consortium
in order to help the maintenance of the code, cleaning obsolete features, revisiting the user interface, adding/replacing
some new numerical schemes, etc... 

 Not surprisingly, the porting from 3.6 to 4.0 requires some work, at different levels.  This document was written 
while performing this porting, together with NEMO4 discovery. After a reminder on DCM concept, it holds basically 
two parts :  (A)  corresponds to all the changes in the code itself, and the compilation procedure.  (B) is dedicated 
to the impacts of the changes at run-time.

This work is done in the frame of the development of a fully new version of the Drakkar Config Manager (DCM). 

##  Reminder of the DCM concept:
 DCM is an integrated software tool allowing an easy deployment of a model configuration based on NEMO. For the model code,
it allows an easy management and tracking of code modifications, that may arise for instance when a new feature is being 
developped, or for a particular model configuration. 

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

DCM propose a set of tools (mainly bash scripts) in order to deal with the file layout automatically. Since the first version,
we maintain almost the same syntax for these tools, maintaining a friendly user interface through the NEMO versions. See the [user manual](user_manual.md).


## A: Fortran code, compilation etc
 In DCM_4.0, the DRAKKAR modifications on fortran modules are all isolated from the standard code by a sequence of

    #if defined key_drakkar
      drakkar code ...
    #endif
 So that, it is very easy to identify the DRAKKAR modifications (using *meld* for instance). When a modification of a
namelist block is introduced ( *e.g*  add extra variables), we took the convention to put use a new namelist block whose
name is just the same with just *_drk* appended at the end of the block name. For example :

    &namtra_dmp
    ... NEMOREF std namelist
    /
    &namtra_dmp_drk
    ... DRAKKAR extension
    /
 Doing so, a user configuration can be run very easily switching off the DRAKKAR/CONFIG-CASE modifications, at least 
for testing purpose. This is a clear improvement with respect to previous versions of DCM (NEMODRAK).

### Code

### Namelists


## B:  Run time modifications
### Domain_cfg concept

### xml files

