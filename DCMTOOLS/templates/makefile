# DCM4.0
# ------------------------------------------
# --- Variables to be customized by users ---

# CONFIG = generic name of the configuration
# CASE   = case id for this variant
# CASEREF = reference case for this configuration (a former case ID or 'none' )


CONFIG = <CONFIG>
CASE = <CASE>
CASEREF = 'none'

# PREV_CONFIG  (type the full path or 'none' )
PREV_CONFIG = 'none'

# TEST_CASE ( one of CANAL  ICEDYN  ISOMIP  LOCK_EXCHANGE  OVERFLOW VORTEX  WAD)
#  NO QUOTE ' ' AROUND TEST_CASE
TEST_CASE = none

# NEMO REFERENC CASE ( one of AGRIF_DEMO AMM12 C1D_PAPA GYRE_BFM GYRE_PISCES ORCA2_ICE_PISCES
#                      ORCA2_OFF_PISCES ORCA2_OFF_TRC ORCA2_SAS_ICE  SPITZ12 )
#  NO QUOTE ' ' AROUND NEMO_REF_CASE
NEMO_REF_CASE = none


# REFONLY ( yes or no : if yes, only reference is used, not the drakkar customs)
REFONLY = 'no'

# MACHINE : name of architecture for compilation. Assume that a arch_${MACHINE}.fcm exists in ARCH 
# to know available arch_xxx.fcm, use dcm_lsarch. You can edit/add your own arch_xxx.fcm in
# CONFIG-CASE/ARCH/ directory
MACHINE = GCC_MEOM

# COMPUTER is the name of the machine on which the experiement will be run. Used for make ctl only
COMPUTER = jean-zay  # one of jean-zay irene occigen ...

# NCOMPIL_PROC : number of procs to use for the compilation of the code.
NCOMPIL_PROC = 1

# set to 'use' the NEMO component required for this config, to anything else  if notused
OCE = 'use'
ICE = 'use'
TOP = 'notused'
OFF = 'notused'
NST = 'notused'
SAS = 'notused'
MY_SRC='notused'
# if set to 'check' maintain an 'install_history' file giving 
#  the date of install and the git hash of $HOMEDCM
GIT = 'check'

# --- All customizations by users is done above this line ... Only experts touch things below ... ---
#----------------------------------------------------------------------------------------------------

# REF is the directory name holding the reference NEMOGCM tree from a NEMO-TEAM download
# CUS is the directory name holding the permanent customized files not yet under the reference, but valid
#     for all the users. Will overwrite the REF during install

REF = $(REFDIR)
CUS = $(CUSDIR)
LOC = $(PWD)

CONFIG_DIR = $(UDIR)/CONFIG_$(CONFIG)/$(CONFIG)-$(CASE)
COMPIL_DIR = $(CDIR)/W$(CONFIG)-$(CASE)
WK =  $(COMPIL_DIR)/cfgs/$(CONFIG)-$(CASE)/WORK

CFILE := $(CONFIG)-$(CASE)

# CFILE is used below for setting the ref_cfgs.txt file in TOOLS/COMPILE
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
	@echo '... compile nemo  with fcm ...'
	@cp CPP.keys $(PDIR)/RUN_$(CONFIG)/$(CONFIG)-$(CASE)/EXE/
	@dcm_cpp2fcm > $(COMPIL_DIR)/cfgs/$(CONFIG)-$(CASE)/cpp_$(CONFIG)-$(CASE).fcm
	@cd $(COMPIL_DIR) ; ./makenemo -r $(CONFIG)-$(CASE)  -j $(NCOMPIL_PROC)
	@ln -sf $(COMPIL_DIR)/cfgs/$(CONFIG)-$(CASE)/BLD/bin/nemo.exe $(PDIR)/RUN_$(CONFIG)/$(CONFIG)-$(CASE)/EXE/nemo4.exe


copyfile:
ifeq ($(CASEREF),'none')
	@echo "no CASEREF to take into account"
else
	@echo " CASEREF=$(CASEREF)"
  ifeq ($(OCE),'use') 
	@-cp -f -p -R ../$(CONFIG)-$(CASEREF)/src/OCE/  $(COMPIL_DIR)/src/
	@-cp -f -p -R ../$(CONFIG)-$(CASEREF)/ext/IOIPSL/*  $(COMPIL_DIR)/ext/IOIPSL/src/.
  endif
  ifeq ($(ICE),'use')
	@-cp -f -p -R ../$(CONFIG)-$(CASEREF)/src/ICE/  $(COMPIL_DIR)/src
  endif
  ifeq ($(TOP),'use')
	@-cp -f -p -R ../$(CONFIG)-$(CASEREF)/src/TOP/  $(COMPIL_DIR)/src
  endif
  ifeq ($(OFF),'use')
	@-cp -f -p -R ../$(CONFIG)-$(CASEREF)/src/OFF/  $(COMPIL_DIR)/src
  endif
  ifeq ($(NST),'use')
	@-cp -f -p -R ../$(CONFIG)-$(CASEREF)/src/NST/  $(COMPIL_DIR)/src
  endif
  ifeq ($(SAS),'use')
	@-cp -f -p -R ../$(CONFIG)-$(CASEREF)/src/SAS/  $(COMPIL_DIR)/src
  endif
endif
ifeq ($(OCE),'use') 
	@-cp -f -p -R  $(CONFIG_DIR)/src/OCE/ $(COMPIL_DIR)/src
	@-cp -f -p -R  $(CONFIG_DIR)/ext/IOIPSL/*  $(COMPIL_DIR)/ext/IOIPSL/src/.
endif
ifeq ($(ICE),'use')
	@-cp -f -p -R  $(CONFIG_DIR)/src/ICE/ $(COMPIL_DIR)/src
endif
ifeq ($(TOP),'use')
	@-cp -f -p -R  $(CONFIG_DIR)/src/TOP/ $(COMPIL_DIR)/src
endif
ifeq ($(OFF),'use')
	@-cp -f -p -R  $(CONFIG_DIR)/src/OFF/ $(COMPIL_DIR)/src
endif
ifeq ($(NST),'use')
	@-cp -f -p -R  $(CONFIG_DIR)/src/NST/ $(COMPIL_DIR)/src
endif
ifeq ($(SAS),'use')
	@-cp -f -p -R  $(CONFIG_DIR)/src/SAS/ $(COMPIL_DIR)/src
endif
# MY_SRC, if used must be the last one to copy
ifeq ($(MY_SRC),'use')
	@-cp -f -p -R  $(CONFIG_DIR)/src/MY_SRC/ $(COMPIL_DIR)/cfgs/$(CONFIG)-$(CASE)/
endif
	@-cp -f -p -R $(CONFIG_DIR)/arch $(COMPIL_DIR)
	@gmake links
#	@-rm -f $(COMPIL_DIR)/CONFIG/$(CONFIG)-$(CASE)/MY_SRC/*.[Ffh]90
#	@-cp -f -p *.[Ffh]90 $(COMPIL_DIR)/CONFIG/$(CONFIG)-$(CASE)/MY_SRC/

copynamelist:
#	@rcp $(CONFIG_DIR)/CTL/namelist_ice.NATL4 uqbar:/u/rech/cli/rcli099/$(CONFIG)-$(CASE)/namelist_ice.NATL4
#	@rcp $(CONFIG_DIR)/CTL/namelist.NATL4 uqbar:/u/rech/cli/rcli099/$(CONFIG)-$(CASE)/namelist.NATL4


##############

install:   copyref update_history  WORK cleansvn cleangit

reinstall: copyref  

# Target copyconfig to duplicate an existant configuration to local directory
copyconfig:
ifeq ($(PREV_CONFIG),'none')
	@echo "no previous configuration defined"
else
	@echo "   copy a previous configuration $(PREV_CONFIG) to local directory"
	-cp -R $(PREV_CONFIG)/src/OCE ./src/
	-cp -R $(PREV_CONFIG)/src/ICE ./src/
	-cp -R $(PREV_CONFIG)/src/TOP ./src/ 
	-cp -R $(PREV_CONFIG)/src/NST ./src/ 
	-cp -R $(PREV_CONFIG)/src/SAS ./src/ 
	-cp -R $(PREV_CONFIG)/ext/IOIPSL ./ext/
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
	-cp  $(PREV_CONFIG)/nam* $(PDIR)/RUN_$(CONFIG)/$(CONFIG)-$(CASE)/CTL/
endif
	@gmake links

# Target test_case : import respective test_case into the current config
test_case:
ifeq ($(TEST_CASE),'none')
	@echo "no test_case configuration defined"
else
	@echo "   copy a test case  configuration $(TEST_CASE) to local directory"
	-cp -f -R $(REFDIR)/tests/$(TEST_CASE)/MY_SRC ./src/
	-cat $(REFDIR)/tests/$(TEST_CASE)/cpp_$(TEST_CASE).fcm | awk '{for (i=2 ; i <= NF ; i++ ) print $$i }' > CPP.keys
	-cat $(REFDIR)/tests/demo_cfgs.txt | grep $(TEST_CASE) | sed -e "s/$(TEST_CASE)/$(CONFIG)-$(CASE)/" > $(CONFIG)-$(CASE)_cfgs.txt
endif
	@gmake links

# Target nemo_ref_case : import respective nemo reference case into the current config
nemo_ref_case:
ifeq ($(NEMO_REF_CASE),'none')
	@echo "no NEMO reference case configuration defined"
else
	@echo "   copy a nemo reference  case  configuration $(NEMO_REF_CASE) to local directory"
	-cp -f -R $(REFDIR)/cfgs/$(NEMO_REF_CASE)/MY_SRC ./src/
	-cat $(REFDIR)/cfgs/$(NEMO_REF_CASE)/cpp_$(NEMO_REF_CASE).fcm | awk '{for (i=2 ; i <= NF ; i++ ) print $$i }' > CPP.keys
	-cat $(REFDIR)/cfgs/ref_cfgs.txt | grep $(NEMO_REF_CASE) | sed -e "s/$(NEMO_REF_CASE)/$(CONFIG)-$(CASE)/" > $(CONFIG)-$(CASE)_cfgs.txt
endif
	@gmake links


# Target copyref is to be done each time we need to rebuild the
# source code from the reference sources.
copyref:
	@echo "    copy reference version to $(COMPIL_DIR) "
	@cp -f -R $(REF)/src $(COMPIL_DIR)/.
	@cp -f -R $(REF)/ext $(COMPIL_DIR)/.
	@cp -f -R $(REF)/mk  $(COMPIL_DIR)/.
	@cp -f -R $(REF)/arch  $(COMPIL_DIR)/.
	@mkdir -p $(COMPIL_DIR)/cfgs/$(CONFIG)-$(CASE)/EXPREF
	@mkdir -p $(COMPIL_DIR)/cfgs/$(CONFIG)-$(CASE)/MY_SRC
	@cp -f -R $(REF)/makenemo $(COMPIL_DIR)/
ifeq ($(REFONLY),'yes')
	@echo "    WARNING : DRAKKAR version is not used "
else
	@echo "    copy DRAKKAR version to $(COMPIL_DIR) "
	@cp -f -R $(CUS)/src $(COMPIL_DIR)/.
	@cp -f -R $(CUS)/ext $(COMPIL_DIR)/.
	@cp -f -R $(CUS)/mk  $(COMPIL_DIR)/.
	@cp -f -R $(CUS)/arch  $(COMPIL_DIR)/.
endif
	@gmake copyfile
ifeq ("$(wildcard $(CONFIG)-$(CASE)_cfgs.txt)","")
	@echo $(CFILE) > $(COMPIL_DIR)/cfgs/ref_cfgs.txt
else
	@cp $(CONFIG)-$(CASE)_cfgs.txt $(COMPIL_DIR)/cfgs/ref_cfgs.txt
endif
	@dcm_cpp2fcm > $(COMPIL_DIR)/cfgs/$(CONFIG)-$(CASE)/cpp_$(CONFIG)-$(CASE).fcm
	@echo  ./makenemo -r $(CONFIG)-$(CASE) -m $(MACHINE)  -j0
	@cd $(COMPIL_DIR) ;  ./makenemo -r $(CONFIG)-$(CASE) -m $(MACHINE)  -j0

WORK:
	@ln -s $(WK) WORK

cleaninst:
	@-\rm -rf $(CDIR)/W$(CONFIG)-$(CASE)/*

clean:
	@cd $(COMPIL_DIR) ; ./makenemo -r $(CONFIG)-$(CASE) clean 

# restore target is used to restore a SINGLE file from the 3levels references
# must be invoked like : make restore RESTORE=OPA_SRC/SBC/ocesbc.F90  (for instance, with the full path relative to NEMO)
# it is used automatically with rmfile.

restore:
ifeq ($(RESTORE),)
	@echo "    Nothing to restore, guy !"
else
	@echo "    Restoring  $(RESTORE) from REFERENCES "
	@-cp  $(REF)/$(RESTORE)  $(COMPIL_DIR)/$(RESTORE)
  ifeq ($(REFONLY),'yes') 
	@echo "    WARNING : nothing restored from DRAKKAR "
  else
    ifneq ("$(wildcard $(CUS)/$(RESTORE))","")
	@echo "    Restoring  $(RESTORE) from DRAKKAR "
	@-cp  $(CUS)/$(RESTORE)  $(COMPIL_DIR)/$(RESTORE)
    else
	@echo "         ===>  $(RESTORE) not in  DRAKKAR "
    endif
  endif
  ifeq ($(CASEREF),'none')
	@echo "         ===> no CASEREF to take into account"
  else
	@echo "    Restoring  $(RESTORE) from CASEREF=$(CASEREF) "
	@-cp  ../$(CONFIG)-$(CASEREF)/$(RESTORE)  $(COMPIL_DIR)/$(RESTORE)
  endif
endif


restoreioipsl:
# likely not working. Fix later if necessary.
ifeq ($(RESTORE),)
	@echo "    Nothing to restore, guy !"
else
	@echo "    Restoring  $(RESTORE) from REFERENCES "
	@-cp  $(REF)/ext/IOIPSL/src/$(RESTORE)  $(COMPIL_DIR)/ext/IOIPSL/src/$(RESTORE)
  ifeq ($(REFONLY),'yes') 
	@echo "    WARNING : nothing restored from DRAKKAR "
  else
	@echo "    Restoring  $(RESTORE) from DRAKKAR "
	@-cp  $(CUS)/ext/IOIPSL/src/$(RESTORE)  $(COMPIL_DIR)/ext/IOIPSL/src/$(RESTORE)
  endif
  ifeq ($(CASEREF),'none')
	@echo "no CASEREF to take into account"
  else
	@echo " CASEREF=$(CASEREF) "
	@-cp  ../$(CONFIG)-$(CASEREF)/ext/$(RESTORE)  $(COMPIL_DIR)/ext/IOIPSL/src/$(RESTORE)
  endif
endif
links:
	@-ln -sf src/OCE/*.*0 .
	@-ln -sf src/OCE/*/*.*0 .
	@-ln -sf src/ICE/*.*0 .
	@-ln -sf src/TOP/*.*0 .
	@-ln -sf src/TOP/*/*.*0 .
	@-ln -sf src/NST/*.*0 .
	@-ln -sf src/SAS/*.*0 .
	@-ln -sf ext/src/IOIPSL/*.*0 .
	@-ln -sf arch/*.fcm .
ifneq ("$(wildcard src/MY_SRC/*)","")
	@-ln -sf src/MY_SRC/*.*0 .
endif

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

ctl:
	@dcm_mkctl -c $(CONFIG)-$(CASE) -m $(COMPUTER)
	@mkdir -p $(PDIR)/RUN_$(CONFIG)/$(CONFIG)-$(CASE)/CTL/XML
	@mkdir -p $(PDIR)/RUN_$(CONFIG)/$(CONFIG)-$(CASE)/CTL/NAMELIST
ifeq ($(REFONLY),'yes')
	@cd $(PDIR)/RUN_$(CONFIG)/$(CONFIG)-$(CASE)/CTL/XML ; cp $(HOMEDCM)/NEMOREF/NEMO4/cfgs/SHARED/*xml ./
	@cd $(PDIR)/RUN_$(CONFIG)/$(CONFIG)-$(CASE)/CTL/NAMELIST ; cp $(HOMEDCM)/NEMOREF/NEMO4/cfgs/SHARED/namelist* ./
else
	@cd $(PDIR)/RUN_$(CONFIG)/$(CONFIG)-$(CASE)/CTL/XML ; cp $(HOMEDCM)/NEMOREF/NEMO4/cfgs/SHARED/*.xml ./
	@cd $(PDIR)/RUN_$(CONFIG)/$(CONFIG)-$(CASE)/CTL/XML ; cp $(HOMEDCM)/DRAKKAR/NEMO4/cfgs/CONFIG-CASE/EXPREF/*.xml ./
	@cd $(PDIR)/RUN_$(CONFIG)/$(CONFIG)-$(CASE)/CTL/NAMELIST ; cp $(HOMEDCM)/NEMOREF/NEMO4/cfgs/SHARED/namelist* ./
	@cd $(PDIR)/RUN_$(CONFIG)/$(CONFIG)-$(CASE)/CTL/NAMELIST ; cp $(HOMEDCM)/DRAKKAR/NEMO4/cfgs/CONFIG-CASE/EXPREF/namelist* ./
	@cd $(PDIR)/RUN_$(CONFIG)/$(CONFIG)-$(CASE)/CTL/NAMELIST ; \
           sed -e "s/<CONFCASE>/$(CONFIG)-$(CASE)/"  namelist_ref > ztmp ; \
           mv ztmp namelist.$(CONFIG)-$(CASE)
	@-cd $(PDIR)/RUN_$(CONFIG)/$(CONFIG)-$(CASE)/CTL/NAMELIST ; \
           sed -e "s/<CONFCASE>/$(CONFIG)-$(CASE)/"  namelist_ice_ref > ztmp ; \
           mv ztmp namelist_ice.$(CONFIG)-$(CASE)
	@cd $(PDIR)/RUN_$(CONFIG)/$(CONFIG)-$(CASE)/CTL/NAMELIST ; rm -f namelist*_ref 
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
	@echo  - make copyconfigall : clone PREV_CONFIG fortran files, CPP_KEYS, and fcm files in the current directory
	@echo  - make test_case     : clone TEST_CASE fortran files, CPP_KEYS, and fcm files in the current directory
	@echo  - make nemo_ref_case : clone NEMO_REF_CASE fortran files, CPP_KEYS, and fcm files in the current directory
	@echo  - make cleaninst     : remove previous installation completly
	@echo  - make clean         : perform a makenemo clean in $(COMPIL_DIR)/CONFIG
	@echo  - make links   : redo links in the local directory
	@echo  - make ctl     : create and populates the $(PDIR)/RUN_$(CONFIG)/$(CONFIG)-$(CASE)/CTL directory
	@echo "                   " - copy reference namelists in CTL/NAMELIST
	@echo "                   " - copy reference xml files in CTL/XML
	@echo  - make help    : display this short help message
