# Using tags with DCM
## Policy of the NEMO revision and release.
A recent evolution in the way of naming NEMO revision has been published in March, 2020. Details are given on this [changeset message](https://forge.ipsl.jussieu.fr/nemo/changeset/12578 ).

The main idea is that all releases now correspond to a stable code.  For instance in the svn tree under releases/ 

 ```
   .../releases/r4.0
      r4.0.0
      r4.0.1
      r4.0.2
      r4.0-HEAD
 ```

In fact, for r4.0.0 and  r4.0.1 the old names (release-4.0, release-4.0.1) are still in use in order not to change existing things. Therefore the actual layout is:


 ```
   .../releases/
         release-3.4
         release-3.6
         release-4.0
         release-4.0.1
         r4.0/
           r4.0.2
           r4.0-HEAD
 ```

The important thing to note is that r4.0.2, is a stable tested version where some bugfixes where added to r4.0.1. r4.0-HEAD is the ongoing revision (trunk) where new bugfixes are being tested for a while.  After some time this r4.0-HEAD revision will be tagged r4.0.3 and HEAD will go on. Once tagged with a minor revision number, __code is not changed__.

## Implementation of the NEMO policy into DCM
DCM will follow the same rules, except that DCM versioning system is `git` instead of `subversion`. 

Git tags are created foreach NEMO revision (v4.0.0, v4.0.1 etc).  The master  or main branch will correspond to the r.0-HEAD, maybe with some delay.

Remember that in the NEMOREF directory, there is a bash script [getnemoref.sh](../DCMTOOLS/NEMOREF/getnemoref.sh), which will get the NEMO revision in phase with the actual DRAKKAR code in DCM. In case of mismatch, the comand `dcm_version` will indicate that something is not coherent between DRAKKAR and NEMOREF.

## Some usefull git command dealing with tags

   - In order to list the tag available after cloning a working copy :

    ```
     git tag
     or
     git tag -l
    ```
    For instance, 

    ```
    $ git tag
    v4.0.0
    v4.0.1
    ```

   -  In order to work with the code corresponding to *e.g.* v4.0.0, then do a git checkout of the tag:

   ```
    $ git checkout v4.0.0
    Note: checking out 'v4.0.0'.

    You are in 'detached HEAD' state. You can look around, make experimental
    changes and commit them, and you can discard any commits you make in this
    state without impacting any branches by performing another checkout.

    If you want to create a new branch to retain commits you create, you may
    do so (now or later) by using -b with the checkout command again. Example:

    git checkout -b new_branch_name

    HEAD is now at d26cd32... add capability for clim forcing and Renault parametrization
    $ git br
    * (no branch)
    master

   ```

The warning message indicated that then you are in detached HEAD state, meaning that no (easy) modifications can be done. In fact, this is the full philosophy of using a tag. It is a snapshot of the code at a given revision.  In this sense, it is not a branch. 
 
## Important warning using tags.
One particularity  of the system is the fact that NEMOREF is not part of the git repository.  NEMOREF is different from one DCM tag to the other. As stated above, [getnemoref.sh](../DCMTOOLS/NEMOREF/getnemoref.sh) script allows you to get the right version. But when using tags, you need to checkout the tag **before**  getting the NEMOREF code. And if for some reason you want to be back to the master, you need to  get NEMOREF again. 

In fact, I strongly recommend to use different working copy when working with tags. It is safer and you do not have temptation on going back to master ! This recommendation is coherent with the use of the modules.  
