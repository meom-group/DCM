#!/bin/bash
date
set -x
########################################################################
#       2. PATHNAME   AND  VARIABLES INITIALISATION                    #
##^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
# Some FLAGS (formely deduced from cpp.options) 1= yes, 0= no
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# non standard features (even in DRAKKAR) ( no namelist nor cpp keys defined for that ! ) 
 UBAR_TIDE=0                          # 2D tidal bottom friction
 WAFDMP=0                             # Use WAter Flux DaMPing ( read previous SSS damping climatology in a file)

 RST_SKIP=1                           # if set, checking of the existence of the full set of restart files is disable (save time !)
 # next flags should be set to 1 if using DCM rev > 1674, to 0 otherwise.
 RST_DIRS=1                           # if set, assumes that restart files are written on multiple directories.
 RST_READY=1                          # if set assumes that restart file are ready to be read by NEMO (no links).

#########################################################################

 CONFIG=<CONFIG>
 CASE=<CASE>
 CONFIG_CASE=${CONFIG}-${CASE}

# Environmemt and miscelaneous
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
login_node=node    # usefull on jade or an any machines requiring scp or ssh to access remote data
MAILTO=<MAILTO>
ACCOUNT=none       # account number for project submission (e.g curie, vayu ...)
QUEUE=none         # queue name (e.g. curie )

# Directory names
#~~~~~~~~~~~~~~~~
# 
WORKDIR=/scratch/$USER
TMPDIR=$WORKDIR/TMPDIR_${CONFIG_CASE}
MACHINE=$<MACHINE>

case  $MACHINE  in
( occigen) SUBMIT=sbatch  ;;
( irene  ) SUBMIT=ccc_msub ;;
( ada    ) SUBMIT=SUBMIT=llsubmit ;;
( *      )  echo $MACHINE not yet supported for SUBMIT definition
esac

SUBMIT_SCRIPT=${CONFIG_CASE}_<MACHINE>.sh   # name of the script to be launched by run_nemo in CTL

if [ ! -d ${TMPDIR} ] ; then mkdir $TMPDIR ; fi

#
# Directory on the storage file system (F_xxx)
F_S_DIR=${SDIR}/${CONFIG}/${CONFIG_CASE}-S       # Stockage
F_R_DIR=${SDIR}/${CONFIG}/${CONFIG_CASE}-R       # Restarts
F_I_DIR=${SDIR}/${CONFIG}/${CONFIG}-I            # Initial + data
F_DTA_DIR=${SDIR}/${CONFIG}/${CONFIG}-I          # data dir
F_FOR_DIR=${SDIR}/DATA_FORCING/ERAinterim/ALL    # in function 3.2
F_OBC_DIR=${SDIR}/${CONFIG}/${CONFIG}-I/OBC      # OBC files
F_BDY_DIR=${SDIR}/${CONFIG}/${CONFIG}-I/BDY      # BDY files
F_MASK_DIR=${SDIR}/${CONFIG}/${CONFIG}-I/MASK    # AABW damping , Katabatic winds
F_INI_DIR=${SDIR}/${CONFIG}/${CONFIG}-I/          
F_WEI_DIR=$SDIR/DATA_FORCING/ERAinterim/ALL

F_OBS_DIR=/ccc/work/cont003/drakkar/drakkar      # for OBS operator
  F_ENA_DIR=${P_OBS_DIR}/ENACT-ENS
  F_SLA_DIR=${P_OBS_DIR}/j2

# Directories on the production machine (P_xxx)
P_S_DIR=$WORKDIR/${CONFIG}/${CONFIG_CASE}-S
P_R_DIR=$WORKDIR/${CONFIG}/${CONFIG_CASE}-R
P_I_DIR=$WORKDIR/${CONFIG}/${CONFIG}-I                  # mirror on the production machine of the F_I_DIR
P_DTA_DIR=$WORKDIR/${CONFIG}/${CONFIG}-I                # mirror on the production machine of the F_I_DIR
P_FOR_DIR=${WORKDIR}/DATA_FORCING/ERAinterim/ALL        # forcing files
P_OBC_DIR=${WORKDIR}/${CONFIG}/${CONFIG}-I/OBC          # OBC files
P_BDY_DIR=${WORKDIR}/${CONFIG}/${CONFIG}-I/BDY          # BDY files
P_WEI_DIR=${WORKDIR}/DATA_FORCING/ERAinterim/ALL        # weight files

P_CTL_DIR=${PDIR}/RUN_${CONFIG}/${CONFIG_CASE}/CTL      # directory from which the job is  launched
P_CDF_DIR=${PDIR}/RUN_${CONFIG}/${CONFIG_CASE}/CTL/CDF  # directory from which the diags are launched
P_EXE_DIR=${PDIR}/RUN_${CONFIG}/${CONFIG_CASE}/EXE      # directory where to find opa
P_UTL_DIR=${UDIR}/UTILS                                 # root directory of the build_nc programs (under bin )
P_XIOS_DIR=${WORKDIR}/XIOS                              # root directory of the XIOS library and xios_server.exe

P_OBS_DIR=/ccc/work/cont003/drakkar/drakkar     # for OBS operation
  P_ENA_DIR=${P_OBS_DIR}/ENACT-ENS
  P_SLA_DIR=${P_OBS_DIR}/j2

RUNTOOLS=${WORKDIR}/DEV/RUN_TOOLS                       # RUNTOOLS directory

# Executable code
#~~~~~~~~~~~~~~~~
EXEC=$P_EXE_DIR/nemo4.exe                              # nemo ...
XIOS_EXEC=$P_XIOS_DIR/bin/xios_server.exe              # xios server (used if code compiled with key_iomput
MERGE_EXEC=$P_UTL_DIR/bin/mergefile_mpp2.exe           # rebuild program (REBUILD_MPP TOOL)  either on the fly (MERGE=1) 
                                                       # or in specific job (MERGE=0). MERGE and corresponding cores number
                                                       # are set in CTL/${SUBMIT_SCRIPT}
                                                       # if you want netcdf4 output use mergefile_mpp4.exe

# In the following, set the name of some files that have a hard coded name in NEMO. Files with variable names
# are directly set up in the corresponding namelist, the script take care of them.
# For the following files, if not relevant set the 'world' name to ''
# set specific file names (data )(world name )                 ;   and their name in NEMO
#--------------------------------------------------------------------------------------------------------
# Bathymetry
BATFILE_LEVEL=                                                              ; OPA_BATFILE_LEVEL=bathy_level.nc
BATFILE_METER=ORCA025_bathy_etopo1_gebco1_smoothed_coast_corrected_mar10.nc ; OPA_BATFILE_METER=bathy_meter.nc

# Coordinates
COORDINATES=coordinates_ORCA_R025_lombok+ombai_v2.nc                        ; OPA_COORDINATES=coordinates.nc

# Ice initialiazation/damping
ICEINI=Init_Ice_GLORYS1V1_NSIDC_BOOTSTRAP_y1989m01_new.nc                   ; OPA_ICEINI=Ice_initialization.nc
ICEDMP=                                                                     ; OPA_ICEDMP=ice_damping.nc

# Bottom friction enhancement
BFR=orca025_bfr_coef_G45.nc                                                 ; OPA_BFR=bfr_coef.nc  # enhanced bottom coef for Torres

# 3D damping mask (aka AABW stuff)
WDMP=ORCA025.L75_dmp_mask.nc                                                ; OPA_WDMP=dmp_mask.nc


# AHM coef file LDF ( orca2 basically)
AHM2D=                                                                      ; OPA_AHM2D=ahmcoef

# Geothermal flux
GEO=                                                           ; OPA_GEO=geothermal_heating.nc

# TRACER new CFC file ends in 2005  ( probably obsolete or not up to date )
CFC=                                                           ; OPA_CFC=cfc1112.atm
CO2=                                                           ; OPA_CO2=splco2.dat
C14=                                                           ; OPA_C14=c14.dat

# --- not standard but already used in some config --
# Water flux damping 
WAFDMP_CLIM=ORCA025_wdmp_from_MJM95.nc                         ; OPA_WAFDMP_CLIM=wdmp_from_MJM95.nc

# ------------------------------------------------------


#OBC
OBCFILE=obc_ORCA025.L75-G70           # root name for OBC file (eg obc_SOSMOD12.L46-MAL95_north_TS_y1999m00.nc )
  # west open boundary file for T  S U    ( typically ${OBCFILE}_west . If set to 'xxx' this boundary is not open  )
  WESTOBC=xxx
  # east open boundary file for T  S U    ( typically ${OBCFILE}_east . If set to 'xxx' this boundary is not open  )
  EASTOBC=xxx
  # north open boundary file for T  S U   ( typically ${OBCFILE}_north . If set to 'xxx' this boundary is not open  )
  NORTHOBC=${OBCFILE}_north
  # south open boundary file for T  S U   ( typically ${OBCFILE}_south . If set to 'xxx' this boundary is not open  )
  SOUTHOBC=xxx

# Agrif 
# ======
AGRIF_FIXED_GRID=AGRIF_FixedGrids.in                  ; OPA_AGRIF_FIXED_GRID=AGRIF_FixedGrids.in

# Control parameters
# -----------------
MAXSUB=0                # resubmit job till job $MAXSUB
