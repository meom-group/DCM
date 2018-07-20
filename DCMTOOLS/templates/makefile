# DCM4.0
# ------------------------------------------
# --- Variables to be customized by users ---

# CONF = generic name of the configuration
# CASE   = case id for this variant
# CASEREF = reference case for this configuration (a former case ID or 'none' )


CONF = CONF
CASE = CASE
CASEREF = 'none'

# PREV_CONFIG  (type the full path or 'none' )
PREV_CONFIG = 'none'

# REFONLY ( yes or no : if yes, only reference is used, not the drakkar customs)
REFONLY = 'no'

# MACHINE : name of architecture for compilation. Assume that a arch_${MACHINE}.fcm exists in ARCH 
# to know available arch_xxx.fcm, use dcm_lsarch.ksh. You can edit/add your own arch_xxx.fcm in
# CONF-CASE/ARCH/ directory
MACHINE = ifort_meol

# NCOMPIL_PROC : number of procs to use for the compilation of the code.
NCOMPIL_PROC = 1

# set to 'use' the NEMO component required for this config, to anything else
# in the actual state of NEMO, OPA and LIM must be set to use for correct compilation
#           if notused
OCE = 'use'
ICE = 'use'
TOP = 'notused'
OFF = 'notused'
NST = 'notused'
SAS = 'notused'
# if set to 'check' maintain an 'install_history' file giving 
#  the date of install and the svn info of $HOMEDCM
SVN = 'no'
GIT = 'check'

# --- All customizations by users is done above this line ... Only experts touch things below ... ---
#----------------------------------------------------------------------------------------------------

# REF is the directory name holding the reference NEMOGCM tree from a NEMO-TEAM download
# CUS is the directory name holding the permanent customized files not yet under the reference, but valid
#     for all the users. Will overwrite the REF during install

REF = $(REFDIR)
CUS = $(CUSDIR)
LOC = $(PWD)

CONFIG_DIR = $(UDIR)/CONFIG_$(CONF)/$(CONF)-$(CASE)
COMPIL_DIR = $(CDIR)/W$(CONF)-$(CASE)
WK =  $(COMPIL_DIR)/CONFIG/$(CONF)-$(CASE)/WORK

CFILE := $(CONF)-$(CASE)

# CFILE is used below for setting the cfg.txt file in TOOLS/COMPILE
ifeq ($(OCE),'use')
  CFILE:=$(CFILE) OCE
endif

ifeq ($(ICE),'use')
  CFILE:=$(CFILE) ICE
endif

ifeq ($(TOP),'use')
  CFILE:=$(CFILE) TOP
endif

ifeq ($(OFF),'use')
  CFILE:=$(CFILE) OFF
endif

ifeq ($(NST),'use')
  CFILE:=$(CFILE) NST
endif

ifeq ($(SAS),'use')
  CFILE:=$(CFILE) SAS
endif


#####
# define the target that must always be done, in any condition, never up to date !
.PHONY: restore help 

compile: copyfile copynamelist
	@echo ' '
	@echo '... compile opa  with fcm ...'
	@cp CPP.keys $(PDIR)/RUN_$(CONF)/$(CONF)-$(CASE)/EXE/
	@dcm_cpp2fcm > $(COMPIL_DIR)/cfgs/$(CONF)-$(CASE)/cpp_$(CONF)-$(CASE).fcm
	@cd $(COMPIL_DIR) ; ./makenemo -n $(CONF)-$(CASE)  -j $(NCOMPIL_PROC)
	@ln -sf $(COMPIL_DIR)/CONFIG/$(CONF)-$(CASE)/BLD/bin/nemo.exe $(PDIR)/RUN_$(CONF)/$(CONF)-$(CASE)/EXE/opa


copyfile:
ifeq ($(CASEREF),'none')
	@echo "no CASEREF to take into account"
else
	@echo " CASEREF=$(CASEREF)"
  ifeq ($(OCE),'use') 
	@-cp -f -p -R ../$(CONF)-$(CASEREF)/OCE/  $(COMPIL_DIR)/NEMO
	@-cp -f -p -R ../$(CONF)-$(CASEREF)/IOIPSL/*  $(COMPIL_DIR)/EXTERNAL/IOIPSL/src/.
  endif
  ifeq ($(ICE),'use')
	@-cp -f -p -R ../$(CONF)-$(CASEREF)/ICE/  $(COMPIL_DIR)/NEMO
  endif
  ifeq ($(TOP),'use')
	@-cp -f -p -R ../$(CONF)-$(CASEREF)/TOP/  $(COMPIL_DIR)/NEMO
  endif
  ifeq ($(OFF),'use')
	@-cp -f -p -R ../$(CONF)-$(CASEREF)/OFF/  $(COMPIL_DIR)/NEMO
  endif
  ifeq ($(NST),'use')
	@-cp -f -p -R ../$(CONF)-$(CASEREF)/NST/  $(COMPIL_DIR)/NEMO
  endif
  ifeq ($(SAS),'use')
	@-cp -f -p -R ../$(CONF)-$(CASEREF)/SAS/  $(COMPIL_DIR)/NEMO
  endif
endif
ifeq ($(OPA),'use') 
	@-cp -f -p -R  $(CONFIG_DIR)/OPA/ $(COMPIL_DIR)/NEMO
	@-cp -f -p -R  $(CONFIG_DIR)/IOIPSL/*  $(COMPIL_DIR)/EXTERNAL/IOIPSL/src/.
endif
ifeq ($(ICE),'use')
	@-cp -f -p -R  $(CONFIG_DIR)/ICE/ $(COMPIL_DIR)/NEMO
endif
ifeq ($(TOP),'use')
	@-cp -f -p -R  $(CONFIG_DIR)/TOP/ $(COMPIL_DIR)/NEMO
endif
ifeq ($(OFF),'use')
	@-cp -f -p -R  $(CONFIG_DIR)/OFF/ $(COMPIL_DIR)/NEMO
endif
ifeq ($(NST),'use')
	@-cp -f -p -R  $(CONFIG_DIR)/NST/ $(COMPIL_DIR)/NEMO
endif
ifeq ($(SAS),'use')
	@-cp -f -p -R  $(CONFIG_DIR)/SAS/ $(COMPIL_DIR)/NEMO
endif
	@-cp -f -p -R $(CONFIG_DIR)/arch $(COMPIL_DIR)
	@gmake links
#	@-rm -f $(COMPIL_DIR)/CONFIG/$(CONF)-$(CASE)/MY_SRC/*.[Ffh]90
#	@-cp -f -p *.[Ffh]90 $(COMPIL_DIR)/CONFIG/$(CONF)-$(CASE)/MY_SRC/

copynamelist:
#	@rcp $(CONFIG_DIR)/CTL/namelist_ice.NATL4 uqbar:/u/rech/cli/rcli099/$(CONF)-$(CASE)/namelist_ice.NATL4
#	@rcp $(CONFIG_DIR)/CTL/namelist.NATL4 uqbar:/u/rech/cli/rcli099/$(CONF)-$(CASE)/namelist.NATL4


##############

install:   copyref update_history  WORK cleansvn cleangit

reinstall: copyref  

# Target copyconfig to duplicate an existant configuration to local directory
copyconfig:
ifeq ($(PREV_CONFIG),'none')
	@echo "no previous configuration defined"
else
	@echo "   copy a previous configuration $(PREV_CONFIG) to local directory"
	-cp -R $(PREV_CONFIG)/OCE .
	-cp -R $(PREV_CONFIG)/ICE . 
	-cp -R $(PREV_CONFIG)/TOP . 
	-cp -R $(PREV_CONFIG)/NST . 
	-cp -R $(PREV_CONFIG)/IOIPSL . 
	-cp -R $(PREV_CONFIG)/SAS . 
endif
	@gmake links


# Target copyconfigall do copyconfig + copy CPP.keys and template namelists
copyconfigall:  copyconfig 
ifeq ($(PREV_CONFIG),'none')
	@echo "no previous configuration defined"
else
	@echo "   copy a previous configuration $(PREV_CONFIG) to local directory"
	-cp  $(PREV_CONFIG)/arch/* arch/
	-cp  $(PREV_CONFIG)/CPP.keys .
	-cp  $(PREV_CONFIG)/makefile makefile.prev
	-cp  $(PREV_CONFIG)/nam* $(PDIR)/RUN_$(CONF)/$(CONF)-$(CASE)/CTL/
endif
	@gmake links


# Target copyref is to be done each time we need to rebuild the
# source code from the reference sources.
copyref:
	@echo "    copy reference version to $(COMPIL_DIR) "
	@cp -f -R $(REF)/NEMOGCM/NEMO $(COMPIL_DIR)/.
	@cp -f -R $(REF)/NEMOGCM/EXTERNAL $(COMPIL_DIR)/.
	@if [ ! -d $(COMPIL_DIR)/CONFIG ] ; then mkdir $(COMPIL_DIR)/CONFIG ; fi
	@if [ ! -d $(COMPIL_DIR)/CONFIG/$(CONF)-$(CASE) ] ; then mkdir $(COMPIL_DIR)/CONFIG/$(CONF)-$(CASE) ; fi
	@if [ ! -d $(COMPIL_DIR)/CONFIG/$(CONF)-$(CASE)/EXP00 ] ; then mkdir $(COMPIL_DIR)/CONFIG/$(CONF)-$(CASE)/EXP00 ; fi
	@if [ ! -d $(COMPIL_DIR)/CONFIG/$(CONF)-$(CASE)/MY_SRC ] ; then mkdir $(COMPIL_DIR)/CONFIG/$(CONF)-$(CASE)/MY_SRC ; fi
	@cp -f -R $(REF)/NEMOGCM/CONFIG/makenemo $(COMPIL_DIR)/CONFIG/.
	@cp -f -R $(REF)/NEMOGCM/ARCH $(COMPIL_DIR)/.
	@if [ ! -d $(COMPIL_DIR)/TOOLS ] ; then mkdir $(COMPIL_DIR)/TOOLS ; fi
	@cp -f -R $(REF)/NEMOGCM/TOOLS/COMPILE $(COMPIL_DIR)/TOOLS/.
ifeq ($(REFONLY),'yes')
	@echo "    WARNING : DRAKKAR version is not used "
else
	@echo "    copy DRAKKAR version to $(COMPIL_DIR) "
	@cp -f -R $(CUS)/NEMOGCM/NEMO $(COMPIL_DIR)/.
	@cp -f -R $(CUS)/NEMOGCM/EXTERNAL $(COMPIL_DIR)/.
	@cp -f -R $(CUS)/NEMOGCM/CONFIG/makenemo $(COMPIL_DIR)/CONFIG/.
	@cp -f -R $(CUS)/NEMOGCM/ARCH $(COMPIL_DIR)/.
	@cp -f -R $(CUS)/NEMOGCM/TOOLS/COMPILE $(COMPIL_DIR)/TOOLS/.
endif
	@gmake copyfile
	@if [ ! -d $(COMPIL_DIR)/CONFIG/$(CONF)-$(CASE) ] ; \
	then mkdir $(COMPIL_DIR)/CONFIG/$(CONF)-$(CASE) ; fi
	@echo $(CFILE) > $(COMPIL_DIR)/CONFIG/cfg.txt
	@dcm_cpp2fcm > $(COMPIL_DIR)/CONFIG/$(CONF)-$(CASE)/cpp_$(CONF)-$(CASE).fcm
	@cd $(COMPIL_DIR)/CONFIG/ ; ./makenemo -n $(CONF)-$(CASE) -m $(MACHINE) -j0

WORK:
	@ln -s $(WK) WORK

cleaninst:
	@-\rm -rf $(CDIR)/W$(CONF)-$(CASE)/*

clean:
	@cd $(COMPIL_DIR)/CONFIG/ ; ./makenemo -n $(CONF)-$(CASE) clean 

# restore target is used to restore a SINGLE file from the 3levels references
# must be invoked like : make restore RESTORE=OPA_SRC/SBC/ocesbc.F90  (for instance, with the full path relative to NEMO)
# it is used automatically with rmfile.

restore:
ifeq ($(RESTORE),)
	@echo "    Nothing to restore, guy !"
else
	@echo "    Restoring  $(RESTORE) from REFERENCES "
	@-cp  $(REF)/NEMOGCM/NEMO/$(RESTORE)  $(COMPIL_DIR)/NEMO/$(RESTORE)
  ifeq ($(REFONLY),'yes') 
	@echo "    WARNING : nothing restored from DRAKKAR "
  else
	@echo "    Restoring  $(RESTORE) from DRAKKAR "
	@-cp  $(CUS)/NEMOGCM/NEMO/$(RESTORE)  $(COMPIL_DIR)/NEMO/$(RESTORE)
  endif
  ifeq ($(CASEREF),'none')
	@echo "no CASEREF to take into account"
  else
	@echo " CASEREF=$(CASEREF) "
	@-cp  ../$(CONF)-$(CASEREF)/$(RESTORE)  $(COMPIL_DIR)/NEMO/$(RESTORE)
  endif
endif


restoreioipsl:
ifeq ($(RESTORE),)
	@echo "    Nothing to restore, guy !"
else
	@echo "    Restoring  $(RESTORE) from REFERENCES "
	@-cp  $(REF)/NEMOGCM/EXTERNAL/IOIPSL/src/$(RESTORE)  $(COMPIL_DIR)/EXTERNAL/IOIPSL/src/$(RESTORE)
  ifeq ($(REFONLY),'yes') 
	@echo "    WARNING : nothing restored from DRAKKAR "
  else
	@echo "    Restoring  $(RESTORE) from DRAKKAR "
	@-cp  $(CUS)/NEMOGCM/EXTERNAL/IOIPSL/src/$(RESTORE)  $(COMPIL_DIR)/EXTERNAL/IOIPSL/src/$(RESTORE)
  endif
  ifeq ($(CASEREF),'none')
	@echo "no CASEREF to take into account"
  else
	@echo " CASEREF=$(CASEREF) "
	@-cp  ../$(CONF)-$(CASEREF)/IOIPSL/$(RESTORE)  $(COMPIL_DIR)/EXTERNAL/IOIPSL/src/$(RESTORE)
  endif
endif
links:
	@-ln -sf OPA_SRC/*.*0 .
	@-ln -sf OPA_SRC/*/*.*0 .
	@-ln -sf OPATAM_SRC/*.*0 .
	@-ln -sf OPATAM_SRC/*/*.*0 .
	@-ln -sf LIM_SRC_2/*.*0 .
	@-ln -sf LIM_SRC_3/*.*0 .
	@-ln -sf TOP_SRC/*.*0 .
	@-ln -sf TOP_SRC/*/*.*0 .
	@-ln -sf NST_SRC/*.*0 .
	@-ln -sf SAS_SRC/*.*0 .
	@-ln -sf IOIPSL/*.*0 .
	@-ln -sf ARCH/*.fcm .

update_history:
ifeq ($(SVN),'check')
	@date >> install_history
	@svn info $(HOMEDCM) >> install_history
	@echo install_history updated
endif

ifeq ($(GIT),'check')
	@date >> $(LOC)/install_history
	@cd $(HOMEDCM) ; git remote -v >> $(LOC)/install_history ; git log  | head -4  >> $(LOC)/install_history ; cd -
	@echo install_history updated
endif


cleansvn:
	@-find . -name ".svn" -exec \rm -rf {} \;

cleangit:
	@-find . -name ".git" -exec \rm -rf {} \;

help:
	@echo Available targets  \(mostly used\)
	@echo  
	@echo  - make install       : install this configuration on the $(COMPIL_DIR)
	@echo  - make               : compile the code with actuall arch file : arch_$(MACHINE).fcm
	@echo  - make copyconfig    : clone PREV_CONFIG fortran files in the current directory
	@echo  - make copyconfigall : clone PREV_CONFIG fortran files, CPP_KEYS, and fcm files un the current directory
	@echo  - make cleaninst     : remove previous installation completly
	@echo  - make clean         : perform a makenemo clean in $(COMPIL_DIR)/CONFIG
	@echo  - make links   : redo links in the local directory
	@echo  - make help    : display this short help message
