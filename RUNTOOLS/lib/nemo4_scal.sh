# $Id: nemo4.sh$
# function.sh has already been sources in the main script
####### UNDER THIS LINE, YOU DON'T HAVE TO CHANGE ANYTHING #####################

echo Running on $( hostname )
echo NB_NPROC = $NB_NPROC
echo NB_NODES = $NB_NODES
## check existence of directories. Create them if they do'nt exist
 # TMPDIR can be set (if necessary in includefile.sh)
mkdir -p $TMPDIR

## Since curie machine, there might be different kind of WORKING directory:
#    traditional WORKDIR with relatively small quota, but no automatic cleaning : CDIR
#    scratchdir with large  quota, but eventually  cleaned automatically : DDIR
# CDIR is always set in DCM setup, DDIR is set to CDIR if not previously set

DDIR=${DDIR:-$CDIR}
RST_DIR=${RST_DIR:-0}
RST_READY=${RST_READY:-0}

mkdir -p  $P_I_DIR
mkdir -p  $P_R_DIR
mkdir -p  $P_S_DIR
mkdir -p  $P_S_DIR/ANNEX

## Generic name for some directories
CN_DIAOBS=${CONFIG_CASE}-DIAOBS     # receive files from diaobs functionality, if used
CN_DIRRST=${CONFIG_CASE}-RST        # receive restart files
CN_DIRICB=${CONFIG_CASE}-ICB        # receive Iceberg Output files

## -----------------------------------------------------
echo '(1) get all the working tools on the TMPDIR directory'
echo '-----------------------------------------------------'
cd $TMPDIR

#  This file is still in use in NEMO4. (in sbcfwd)
cat << eof > EMPave_old.dat
0 0 0
eof

## clean eventual (?) old files
\rm -f OK*    # remove all OK* file ( including those produced by ensemble runs)
\rm -f STOP*  # remove all STOP* file ( including those produced by ensemble runs)
\rm -f nitend.txt
\rm -f znitend.txt
\rm -f damping*

# usefull when TMPDIR is set dynamically by the system
pwd > waytmp
copy waytmp $P_CTL_DIR/

## copy of system and script tools: from P_CTL_DIR to TMPDIR
echo " [1.1]  copy script and get number ($$) usefull for run"
echo " ======================================================"

rcopy $P_UTL_DIR/bin/datfinyyyy ./
rcopy $P_CTL_DIR/includefile.sh_${JPNI}_${JPNJ}_${JPNIJ}  includefile.sh 

## copy the executable NEMO
set -x
chkfile $EXEC  

if [ $? = 0 ] ; then
    rcopy $EXEC ./nemo4.exe
else
    echo "   ===  ERROR: NEMO4 must me recompiled. Deleted from workdir"
    exit 1
fi

# copy list of CPPkeys ( available in EXE directory
chkfile $P_EXE_DIR/CPP.keys
if [ $? = 0 ] ; then
    rcopy $P_EXE_DIR/CPP.keys ./CPP.keys
else
    echo "   ===  ERROR: CPP.keys not found in EXE directory (weird ? ) "
    exit 1
fi

echo " [1.2]  set flags according to CPP_keys"
echo " ======================================"
# Reminder : NEMO4 CPP keys are :
# key_agrif
# key_asminc
# key_c1d
# key_cice
# key_cice4
# key_cyclone
# key_diadct
# key_diaharm
# key_diahth
# key_diainstant
# key_floats
# key_iomput
# key_mpp_mpi
# key_nemocice_decomp
# key_netcdf4
# key_nosignedzero
# key_oa3mct_v3
# key_oasis3
# key_si3
# key_top
#   key_sed_off
#   key_trdmxl_trc
#   key_trdtrc
# key_vectopt_loop

#  DRAKKAR CPP key
# key_drakkar


IOIPSL=0  # probably never used except with test cases ( if no XIOS ? )

AGRIF=0   ;  if [ $(keychk key_agrif  )  ] ; then AGRIF=1   ; fi
XIOS=0    ;  if [ $(keychk key_iomput )  ] ; then XIOS=1    ; DIROUText='XIOS'    ; fi
ICE=0     ;  if [ $(keychk key_si3    )  ] ; then ICE=1     ; fi
TOP=0     ;  if [ $(keychk key_top )     ] ; then TOP=1     ; fi
FLOAT=0   ;  if [ $(keychk key_floats)   ] ; then FLOAT=1   ; fi
CYCL=0    ;  if [ $(keychk key_cyclone)  ] ; then CYCL=1    ; fi
DIAHARM=0 ;  if [ $(keychk key_diaharm)  ] ; then DIAHARM=1 ; fi

# JMM fix for the time being
CFC=0 ; C14=0 ; MYTRC=0

#@@@@@@@@ XIOS2=0   ;  if [ $(keychk key_xios2  )  ] ; then XIOS2=1   ; fi
#@@@@@@@@ TRDMLD=0  ;  if [ $(keychk key_trdmld )  ] ; then TRDMLD=1  ; fi
#@@@@@@@@ CFC=0     ;  if [ $(keychk key_cfc )     ] ; then CFC=1     ; fi
#@@@@@@@@ FLXISH=0  ;  if [ $(keychk key_iceshelf) ] ; then FLXISH=1  ; fi   ===> ln_isf
#@@@@@@@@ DIAOBS=0  ;  if [ $(keychk key_diaobs )  ] ; then DIAOBS=1  ; fi

# 
## check if we are using new xml layout (ie with files like 04-files.xml
NEWXML=0
if [ $XIOS = 1 ] ; then
   if [ -f $P_CTL_DIR/04-file.xml ] ; then NEWXML=1 ; fi
fi
   
echo "   *** XIOS   = " $XIOS
echo "   *** NEWXML = " $NEWXML

## copy of the control files ( .db and and template namelist )
rcopy $P_CTL_DIR/namelist.${CONFIG_CASE}_${JPNI}_${JPNJ}_${JPNIJ} namelist
rcopy $P_CTL_DIR/$CONFIG_CASE.db ./

if [ $AGRIF = 1 ] ; then
    initagrif 
    for idx in ${agrif_pref[@]} ; do
        rcopy $P_CTL_DIR/${idx}_namelist.${CONFIG_CASE}_${JPNI}_${JPNJ}_${JPNIJ} ${idx}_namelist
        rcopy $P_CTL_DIR/${idx}_namelistio ${idx}_namelistio
    done
fi

## -------------------------------------
echo '(2) Set up the namelist for this run from template'
echo '--------------------------------------------------'
echo " [2.1]  ocean namelist"
echo " ====================="

## exchange <wildcards>  with the correct info from db
no=`tail -1 $CONFIG_CASE.db | awk '{print $1}' `
nit000=`tail -1 $CONFIG_CASE.db | awk '{print $2}' `
nitend=`tail -1 $CONFIG_CASE.db | awk '{print $3}' `

if [ $no != 1 ] ; then
    restart_flag=true
else
    restart_flag=false
fi

sed -e "s/<NN_NO>/$no/" \
    -e "s/<CONFCASE>/$CONFIG_CASE/" \
    -e "s/<NIT000>/$nit000/" \
    -e "s/<NITEND>/$nitend/" \
    -e "s/<RESTART>/$restart_flag/" \
    -e "s@<CN_DIAOBS>@$DDIR/${CN_DIAOBS}.$no@"   \
    -e "s@<CN_DIRICB>@$DDIR/${CN_DIRICB}.$no@"   \
    -e "s@<CN_DIRRST>@$DDIR/${CN_DIRRST}@"   namelist > znamelist1
\cp znamelist1 namelist

if [ $DIAOBS = 1 ] ; then  # modify namelist block namobs (and get data files)
    # Observation operator
    ENACT=0  
    tmp=$(LookInNamelist ln_ena) ; tmp=$(normalize $tmp)
    if [ $tmp = T ] ; then ENACT=1 ; fi

    SLA=0
    tmp=$(LookInNamelist ln_sla) ; tmp=$(normalize $tmp)
    if [ $tmp = T ] ; then SLA=1 ; fi
    echo "   ***  ENACT  = $ENACT"
    echo "   ***  SLA    = $SLA"

    getobs
fi

\cp namelist  namelist_ref
\cp namelist  namelist_cfg

## Agrif namelist update if any
if [ $AGRIF = 1 ] ; then
    # for each agrif subgrid, determine nit000 and nitend according to values for G0
    # in the following arrays nit0[@] and nite[@] are the values to be replaced in the namelists
    # note that index 0 correspond to G0. index idx corresponds to G{idx}
    nit0[0]=$nit000
    nite[0]=$nitend

    for idx in ${agrif_pref[@]} ; do
        nit0[idx]=$(( (${nit0[idx-1]} -1)*${timeref[idx]} + 1 ))
        nite[idx]=$((  ${nit0[idx]} - 1 + ( ${nite[idx-1]} - ${nit0[idx-1]} +1 ) *${timeref[idx]} ))

        sed -e "s/<NN_NO>/$no/"   \
            -e "s/<NIT000>/${nit0[idx]}/" \
            -e "s/<NITEND>/${nite[idx]}/" \
            -e "s@<CN_DIRRST>@$DDIR/${CN_DIRRST}@" \
            -e "s/<RESTART>/$restart_flag/" ${idx}_namelist > z${idx}_namelist1

        \cp z${idx}_namelist1 ${idx}_namelist
        \cp ${idx}_namelist  ${idx}_namelist_ref
        \cp ${idx}_namelist  ${idx}_namelist_cfg
    done
fi

echo "   ***  Check/Create directory : ${CONFIG_CASE}-${DIROUText}.$no"
mkdir -p  $DDIR/${CONFIG_CASE}-${DIROUText}_${JPNI}_${JPNJ}_${JPNIJ}.$no

echo "   ***  Check/Create directory : ${CONFIG_CASE}-${MOOROUText}.$no"
MOOROUText=MOORINGS
#JMM : eliminate creation of this dir ... dirty
#mkdir -p  $DDIR/${CONFIG_CASE}-${MOOROUText}.$no

if [ $DIAOBS = 1 ] ; then 
    echo "   ***  Check/Create directory : ${CN_DIAOBS}.$no"
    mkdir -p $DDIR/${CN_DIAOBS}.$no 
fi


if [ $RST_DIRS = 1 ] ; then 
    echo "   ***  Check/Create directory : ${CN_DIRRST}.$no"
    mkdir -p $DDIR/${CN_DIRRST}.$no
fi

rdt=$(LookInNamelist rn_rdt)

## place holder for time manager (eventually)
if [ $no != 1 ] ; then
    ndastpdeb=`tail -2 $CONFIG_CASE.db | head -1 | awk '{print $4}' `
else
    ndastpdeb=$(LookInNamelist nn_date0)
fi
echo "   ***  Intial date for this run : $ndastpdeb"

year=$(( ndastpdeb / 10000 ))
mmdd=$(( ndastpdeb - year * 10000 ))

if [ $mmdd = 1231 ] ; then
    year=$(( year + 1 ))
fi
year=$( printf "%04d" $year)

rdt=`echo 1 | awk "{ rdt=int($rdt); print rdt}" `
echo "   ***  Time step is: $rdt "

ndays=` echo 1 | awk "{ a=int( ($nitend - $nit000 +1)*$rdt /86400.) ; print a }" `
ndastpfin=`./datfinyyyy $ndastpdeb $ndays `
echo "   ***  $ndays days to run, starting $ndastpdeb ending $ndastpfin"

if [ $TOP = 1 ] ; then 
    echo ' [2.2]  Tracer namelist(s)'
    echo " ========================="
    rcopy $P_CTL_DIR/namelist_top ./
    sed -e "s@<CN_DIRRST>@$DDIR/${CN_DIRRST}@"   namelist_top > ztmpnmtop
    mv ztmpnmtop namelist_top
    cp namelist_top namelist_top_ref
    cp namelist_top namelist_top_cfg
    if [ $CFC = 1    ] ; then rapatrie $CFCATM $P_I_DIR $F_DTA_DIR $NEMO_CFCATM ; fi
  # if [ $C14 = 1    ] ; then rapatrie $CO2 $P_I_DIR $F_DTA_DIR $NEMO_CO2 ; fi
  # if [ $C14 = 1    ] ; then rapatrie $C14 $P_I_DIR $F_DTA_DIR $NEMO_C14 ; fi
    if [ $CFC = 1    ] ; then rcopy $P_CTL_DIR/namelist_cfc ./ ; fi
  # if [ $C14 = 1    ] ; then rcopy $P_CTL_DIR/namelist_c14 ./ ; fi
  # if [ $MYTRC = 1  ] ; then rcopy $P_CTL_DIR/namelist_mytrc ./ ; fi
fi

if [ $ICE != 0 ] ; then
    echo ' [2.3]  Ice namelist'
    echo " ========================="
    rcopy $P_CTL_DIR/namelist_ice.${CONFIG_CASE} namelist_ice
    sed -e "s@<CN_DIRRST>@$DDIR/${CN_DIRRST}@"   namelist_ice > ztmpnmice
    mv ztmpnmice namelist_ice
    cp namelist_ice namelist_ice_ref
    cp namelist_ice namelist_ice_cfg
    if [ $AGRIF = 1 ] ; then
        for idx in ${agrif_pref[@]} ; do
            rcopy $P_CTL_DIR/${idx}_namelist_ice ${idx}_namelist_ice
            cp ${idx}_namelist_ice ${idx}_namelist_ice_ref
            cp ${idx}_namelist_ice ${idx}_namelist_ice_cfg
        done
    fi
fi

# XIOS stuff migrated after ensemble run check !

echo ' [2.5]   Set flags according to namelists'
echo " ========================================"

# Domain cfg file
DOMAINcfg=0
  tmp=$( LookInNamelist ln_read_cfg namelist namcfg ) ; tmp=$( normalize $tmp )
  if [ $tmp = T ] ; then DOMAINcfg=1 ; fi

echo "   ***  DOMAINcfg = $DOMAINcfg"

BDY=0
  tmp=$(LookInNamelist ln_bdy namelist nambdy  ) ; tmp=$(normalize $tmp)
  if [ $tmp = T ] ; then BDY=1 ; fi

# Ice model
ICE_INI=0 ; ICE_DMP=0
if [ $ICE = 1 ] ; then   # SI3
    tmp=$(LookInNamelist ln_iceini namelist_ice namini ) ; tmp=$(normalize $tmp)
    if [ $tmp = T ] ; then ICE_INI=1 ; fi

#   No ice damping so far in SI3/NEMO4 ....
#    tmp=$(LookInNamelist ln_limdmp namelist_ice) ; tmp=$(normalize $tmp)
#    if [ $tmp = T ] ; then ICE_DMP=1 ; fi
fi
echo "   ***  ICE = $ICE"
echo "   ***  ICE_INI = $ICE_INI"
echo "   ***  ICE_DMP = $ICE_DMP"

# Tidal mixing ( Delavergne )
ZDFIWM=0
tmp=$( LookInNamelist ln_zdfiwm  ) ; tmp=$( normalize $tmp ) 
if [ $tmp = T ] ; then ZDFIWM=1  ; fi

# Enhanced bottom friction 
BOOST_DRG_BOT=0
tmp=$( LookInNamelist ln_boost namelist  namdrg_bot ) ; tmp=$( normalize $tmp ) 
if [ $tmp = T ] ; then BOOST_DRG_BOT=1  ; fi
echo "   ***  BOOST_DRG_BOT = $BOOST_DRG_BOT"

# Enhanced top friction 
BOOST_DRG_TOP=0
tmp=$( LookInNamelist ln_boost namelist namdrg_top ) ; tmp=$( normalize $tmp ) 
if [ $tmp = T ] ; then BOOST_DRG_TOP=1  ; fi
echo "   ***  BOOST_DRG_TOP = $BOOST_DRG_TOP"

# tracer damping 
TRADMP=0 
tmp=$(LookInNamelist ln_tradmp) ; tmp=$(normalize $tmp )
if [ $tmp = T ] ; then TRADMP=1 ; fi
echo "   ***  TRADMP = $TRADMP"

# damping mask (e.g. AABW )
AABW_DMP=0  
if [ $TRADMP = 1 ] ; then
    tmp=$(LookInNamelist ln_dmpmask namelist namtra_dmp_drk ) ; tmp=$(normalize $tmp)
    if [ $tmp = T ] ; then AABW_DMP=1 ; fi
fi
echo "   ***  AABW_DMP = $AABW_DMP"

# Lagrangian floats
RFLOAT=0 ; IFLOAT=$FLOAT
if [ $FLOAT = 1 ] ; then
    tmp=$(LookInNamelist ln_rstflo) ; tmp=$(normalize $tmp)
    if [ $tmp = T ] ; then RFLOAT=1 ; fi
fi
echo "   ***  FLOAT  = $FLOAT"
echo "   ***  IFLOAT = $IFLOAT"
echo "   ***  RFLOAT = $RFLOAT"

# Geothermal heating
GEOTH=0
tmp=$(LookInNamelist ln_trabbc) ; tmp=$(normalize $tmp)
if [ $tmp = T ] ; then 
    GEOTH=1
    nn_geoflx=$(LookInNamelist nn_geoflx)
    if [ $nn_geoflx = 2 ] ; then GEOTH=2 ; fi
fi
echo "   ***  GEOTH  = $GEOTH"

# Iceberg calving
ICB=0
tmp=$(LookInNamelist ln_icebergs) ; tmp=$(normalize $tmp)
if [ $tmp = T ] ; then
    ICB=1
    echo "   ***  Check/Create directory : ${CN_DIRICB}.$no"
    mkdir -p $DDIR/${CN_DIRICB}.$no 
fi
echo "   ***  ICB  = $ICB"

# Ice Shelves 
ISF=0
tmp=$(LookInNamelist ln_isf namelist namsbc) ; tmp=$(normalize $tmp)
if [ $tmp = T ] ; then ISF=1   ; fi
echo "   ***  ISF  = $ISF"

# Poleward Transport diagnostics
DIAPTR=0
SUBBAS=0
tmp=$(LookInNamelist ln_diaptr namelist namptr ) ; tmp=$(normalize $tmp)
if [ $tmp = T ] ; then
    DIAPTR=1
    tmp2=$(LookInNamelist ln_subbas namelist namptr) ; tmp2=$(normalize $tmp2)
    if [ $tmp2 = T ] ; then
        SUBBAS=1
    fi
fi
echo "   ***  DIAPTR  = $DIAPTR"
echo "   ***  SUBBAS  = $SUBBAS"

# Diffusion/viscosity read in file 
TRACOEF2D=0 ; TRACOEF3D=0
   tmp=$(LookInNamelist ln_traldf_OFF ) ; tmp=$(normalize $tmp)
   if [ $tmp = F ] ; then   #  Diffusion is ON !
     nn_aht_ijk_t=$(LookInNamelist nn_aht_ijk_t )
     if [ $nn_aht_ijk_t = -20 ] ; then TRACOEF2D=1 ; fi
     if [ $nn_aht_ijk_t = -30 ] ; then TRACOEF3D=1 ; fi
   fi

DYNCOEF2D=0  ; DYNCOEF3D=0
   tmp=$(LookInNamelist ln_dynldf_OFF ) ; tmp=$(normalize $tmp)
   if [ $tmp = F ] ; then   #  Diffusion is ON !
     nn_ahm_ijk_t=$(LookInNamelist nn_ahm_ijk_t )
     if [ $nn_ahm_ijk_t = -20 ] ; then DYNCOEF2D=1 ; fi
     if [ $nn_ahm_ijk_t = -30 ] ; then DYNCOEF3D=1 ; fi
   fi


# Ensemble run [ if not an ensemble run it is mandatory to have both ENSEMBLE_START and ENSEMBLE_END set to -1 ]
ENSEMBLE=0 ; ENSEMBLE_SIZE=1 ; ENSEMBLE_START=-1 ; ENSEMBLE_END=-1 ; ENSEMBLE_RST_IN=0
ENSDIAGS=0

tmp=$(LookInNamelist ln_ensemble) ; tmp=$(normalize $tmp)
if [ $tmp = T ] ; then  ENSEMBLE=1 ; fi

if [ $ENSEMBLE = 1 ] ; then
      # check if inter members diags are requested
    tmp=$(LookInNamelist ln_ens_diag) ; tmp=$(normalize $tmp)
    if [ $tmp = T ] ; then  ENSDIAGS=1 ; fi

      # set ensemble parameters from namelist
    ENSEMBLE_SIZE=$(LookInNamelist nn_ens_size)   
    ENSEMBLE_START=$(LookInNamelist nn_ens_start)
    ENSEMBLE_END=$(( ENSEMBLE_START + ENSEMBLE_SIZE -1 ))
    # next bloc is when members have different namelists (in test and not fully working )
    tmp=$(LookInNamelist ln_ens_namlist) ; tmp=$(normalize $tmp)
    if [ $tmp = T ] ; then  # each member has its own namelist.mmm  ( ocean only so far )
       for member in $(seq $ENSEMBLE_START $ENSEMBLE_END) ; do
          nnn=$(getmember_extension $member  nodot )  # number of the member without .
          rcopy $P_CTL_DIR/namelist.${CONFIG_CASE}.${nnn} namelist.${nnn}
          sed -e "s/NUMERO_DE_RUN/$no/" \
              -e "s/NIT000/$nit000/" \
              -e "s/NITEND/$nitend/" \
              -e "s/RESTART/$restart_flag/" \
              -e "s@CN_DIROUT@$DDIR/${CONFIG_CASE}-${DIROUText}.$no@" \
              -e "s@CN_DIROUT2@$DDIR/${CONFIG_CASE}-${DIROUText}.$no/SSF@" \
              -e "s@CN_DIROUT3@$DDIR/${CONFIG_CASE}-${DIROUText}.$no/5D@" \
              -e "s@CN_DIAOBS@$DDIR/${CN_DIAOBS}.$no@"   \
              -e "s@CN_DIRRST@$DDIR/${CN_DIRRST}@"   namelist.$nnn > znamelist1
              \cp znamelist1 namelist.$nnn
       done
    fi

      # check if members require individual restart file.
      # if not, all members restart from the same file 
      # ie, link will be made to the proper file
    tmp=$(LookInNamelist ln_ens_rst_in) ; tmp=$(normalize $tmp)
    if [ $tmp = T ] ; then ENSEMBLE_RST_IN=1 ; fi
    
      # check the number of available cores vs the required number of cores
    jpnij=$(LookInNamelist jpnij)  # core per member
    core_required=$(( jpnij * ENSEMBLE_SIZE ))
    if [ $core_required != $NB_NPROC ] ; then
        echo 'ERROR: inconsistent number of processors ...'
        exit
    fi
    for member in $(seq $ENSEMBLE_START $ENSEMBLE_END) ; do
        nnn=$(getmember_extension $member  nodot )  # number of the member without .
        mkdir -p  $DDIR/${CONFIG_CASE}-${DIROUText}.${no}/$nnn
        if [ $RST_DIRS = 1 ] ; then
            mkdir -p  $DDIR/${CN_DIRRST}.${no}/$nnn
        fi
    done
fi
echo "   ***  ENSEMBLE  = $ENSEMBLE"
echo "   ***  ENSEMBLE_SIZE  = $ENSEMBLE_SIZE"
echo "   ***  ENSEMBLE_START = $ENSEMBLE_START"
echo "   ***  ENSEMBLE_END   = $ENSEMBLE_END"

# Sochastic parameterization : set STO=1 if at least one of the ln_sto_xxx flag is true
STO=0 ; RSTO=0
tmp=$(LookInNamelist ln_sto_ldf) ; tmp=$(normalize $tmp)
if [ $tmp = T ] ; then STO=1 ; fi
tmp=$(LookInNamelist ln_sto_hpg) ; tmp=$(normalize $tmp)
if [ $tmp = T ] ; then STO=1 ; fi
tmp=$(LookInNamelist ln_sto_eos) ; tmp=$(normalize $tmp)
if [ $tmp = T ] ; then STO=1 ; fi

if [ $STO = 1 ] ; then
    tmp=$(LookInNamelist ln_rststo) ; tmp=$(normalize $tmp)
    if [ $tmp = T ] ; then RSTO=1 ; fi
fi
echo "   ***  STO  = $STO"
echo "   ***  RSTO = $RSTO"

# Use of the observation operator.
DIAOBS=0
 tmp=$(LookInNamelist ln_diaobs namelist_cfg namobs ) ; tmp=$(normalize $tmp)
 if [ $tmp = T ] ; then DIAOBS=1 ; fi
 echo "   *** DIAOBS  = $DIAOBS"

if [ $DIAOBS = 1 ] ; then
    missing_err=0
    # check if datfinyyyy and fbcomb are available
    chkfile $P_UTL_DIR/bin/datfinyyyy
    if [ $? = 0 ] ; then
        rcopy $P_UTL_DIR/bin/datfinyyyy ./datfinyyyy
    else
        echo "   === ERROR : missing datfinyyyy with diaobs in use "
        echo "       (Must be in $P_UTL_DIR/bin. Sources are in $RUN_TOOLS/UTILS )"
        missing_err=$(( missing_err + 1 ))
    fi

    chkfile $P_UTL_DIR/bin/fbcomb.exe
    if [ $? = 0 ] ; then
        rcopy $P_UTL_DIR/bin/fbcomb.exe ./fbcomb.exe
    else
        echo "   === ERROR : missing fbcomb with diaobs in use "
        echo "       (Must be in $P_UTL_DIR/bin. Sources are in TOOLS/OBSTOOLS ) "
        missing_err=$(( missing_err + 1 ))
    fi

    if [ $missing_err != 0 ] ; then exit 1 ; fi

fi


if [ $XIOS = 1 ] ; then
    echo ' [2.6]  XML files for XIOS'
    echo " ========================="
    cp $P_CTL_DIR/*xml ./
#   rcopy $P_CTL_DIR/field_def.xml field_def.xml
#   rcopy $P_CTL_DIR/domain_def.xml domain_def.xml
#   if [ $XIOS2 = 1 ] ; then
#      rcopy $P_CTL_DIR/file_def.xml file_def.xml
#      rcopy $P_CTL_DIR/grid_def.xml grid_def.xml
#   fi
    rcopy $XIOS_EXEC xios_server.exe

    if [ $NEWXML = 1 ] ; then
        if [ $ENSDIAGS = 1 ] ; then
            rcopy $P_CTL_DIR/01-ifile.xml  .
            rcopy $P_CTL_DIR/02-iseczo.xml .
            rcopy $P_CTL_DIR/03-isecme.xml .
        fi
        rcopy $P_CTL_DIR/04-file.xml   .
        rcopy $P_CTL_DIR/05-seczo.xml  .
        rcopy $P_CTL_DIR/06-secme.xml  .

       # build a iodef.xml from template according to the members
        rcopy $P_CTL_DIR/iodef_nemo_context.xml .
        rcopy $P_CTL_DIR/iodef_xios_context.xml .
        cat << eof > iodef.xml
<?xml version="1.0"?>
<simulation>

eof
       # first deal with 1rst member ( eventual inter member diags )
        member=$ENSEMBLE_START
        mmm=$(getmember_extension $member)  # if no ensemble, this function return empty string
        nnn=$(getmember_extension $member nodot )

        if [ $ENSDIAGS = 1 ] ; then
            cat iodef_nemo_context.xml | sed -e "/01-ifile.xml/r  01-ifile.xml"  \
                -e "/02-iseczo.xml/r 02-iseczo.xml" \
                -e "/03-isecme.xml/r 03-isecme.xml" \
                -e "/04-file.xml/r   04-file.xml"   \
                -e "/05-seczo.xml/r  05-seczo.xml"  \
                -e "/06-secme.xml/r  06-secme.xml" > ztemp.xml
        else
            cat iodef_nemo_context.xml | sed -e "/04-file.xml/r   04-file.xml"   \
                -e "/05-seczo.xml/r  05-seczo.xml"  \
                -e "/06-secme.xml/r  06-secme.xml" > ztemp.xml
        fi
        
        cat ztemp.xml | sed -e "s/<NEMO.MEMBER>/nemo$mmm/"  \
            -e "s@<OUTDIR>@$DDIR/${CONFIG_CASE}-XIOS_${JPNI}_${JPNJ}_${JPNIJ}.$no/$nnn@" \
            -e "s/<CASE>/$CASE${mmm}/" >>  iodef.xml
       # Then deal with other members (each identical so far)
        cat iodef_nemo_context.xml | sed -e "/04-file.xml/r   04-file.xml"   \
            -e "/05-seczo.xml/r  05-seczo.xml"  \
            -e "/06-secme.xml/r  06-secme.xml" > ztemp.xml

        memberp1=$(( ENSEMBLE_START + 1 ))
        for member in $(seq $memberp1 $ENSEMBLE_END) ; do
            mmm=$(getmember_extension $member)  # if no ensemble, this function return empty string
            nnn=$(getmember_extension $member nodot )
            cat ztemp.xml | sed -e "s/<NEMO.MEMBER>/nemo$mmm/"  \
                -e "s@<OUTDIR>@$DDIR/${CONFIG_CASE}-XIOS_${JPNI}_${JPNJ}_${JPNIJ}.$no/$nnn@" \
                -e "s@<MOORDIR>@$DDIR/${CONFIG_CASE}-MOORINGS.$no/$nnn@" \
                -e "s/<CASE>/$CASE${mmm}/" >>  iodef.xml
        done

        cat iodef_xios_context.xml >> iodef.xml
        cat << eof >> iodef.xml

</simulation>
eof
        rcopy $P_CTL_DIR/iodef.xml iodef.xml
    fi
    echo "  ***   Customize iodef.xml from template"
    # set <OUTDIR> in iodef.xml
    ndate0=$(LookInNamelist nn_date0)
   for  xml_fil in *.xml ; do
    cat $xml_fil | sed -e "s@<OUTDIR>@$DDIR/${CONFIG_CASE}-XIOS_${JPNI}_${JPNJ}_${JPNIJ}.$no@"  \
        -e "s@<MOORDIR>@$DDIR/${CONFIG_CASE}-MOORINGS.$no@" \
        -e "s/<CONFIG>/$CONFIG/" -e "s/<CASE>/$CASE/" \
        -e "s/<NDATE0>/$ndate0/" > ztmpxml
    mv ztmpxml $xml_fil
   done
#    if [ $XIOS2 = 1 ] ; then
#       cat file_def.xml | sed -e "s@<OUTDIR>@$DDIR/${CONFIG_CASE}-XIOS.$no@"  \
#        -e "s@<MOORDIR>@$DDIR/${CONFIG_CASE}-MOORINGS.$no@" \
#        -e "s/<CONFIG>/$CONFIG/" -e "s/<CASE>/$CASE/" \
#        -e "s/<NDATE0>/$ndate0/" > ztmp
#       mv ztmp file_def.xml
#    fi

fi
#--------------------------------------
echo '(3) Look for input files (According to flags)'
echo '---------------------------------------------'

echo ' [3.1] : configuration files'
echo ' =========================='
## DOMAINcfg
if [ $DOMAINcfg = 1 ] ; then  
    CN_DOMCFG=$(LookInNamelist cn_domcfg ).nc
    rapatrie $CN_DOMCFG $P_I_DIR $F_DTA_DIR $CN_DOMCFG ; fi

## bottom friction file
if [ $BOOST_DRG_BOT = 1 ] ; then             # enhance bottom friction used
    CN_BOOST=$(LookInNamelist sn_boost namelist namdrg_bot_drk ).nc
    rapatrie $CN_BOOST $P_I_DIR $F_DTA_DIR $CN_BOOST
fi

## top friction file
if [ $BOOST_DRG_TOP = 1 ] ; then             # enhance top friction used
    CN_BOOST=$(LookInNamelist sn_boost namelist namdrg_top_drk ).nc
    rapatrie $CN_BOOST $P_I_DIR $F_DTA_DIR $CN_BOOST
fi

## Tidal mixing (Delavergne)
if [ $ZDFIWM = 1 ] ; then
    getzdfiwm
fi

## TIDAL FRICTION ( use a 2D tide velocity in the quadratic friction law)
if [ $UBAR_TIDE  = 1 ] ; then
    rapatrie $BFR_TIDE  $P_I_DIR $F_DTA_DIR $NEMO_BFR_TIDE
fi

## geothermal heating
if [ $GEOTH = 2 ] ; then
    getgeo
fi

## iceberg calving
if [ $ICB = 1 ] ; then
    getcalving
fi

## iceshelve fluxes and/or circulation
if [ $ISF = 1 ] ; then
    getisf
fi

## diaobs : for memory : needed files are already copied (when updating namelist)

## diaptr subbasins mask
if [ $SUBBAS = 1 ] ; then
    rapatrie $SUBBAS  $P_I_DIR $F_INI_DIR $NEMO_SUBBAS
fi

##  Diffusivity/viscosity  prescribed on either 2D or 3D files
if [ $TRACOEF2D = 1 ] ; then
    rapatrie $AHT2D  $P_I_DIR $F_DTA_DIR $NEMO_AHT2D
fi

if [ $TRACOEF3D = 1 ] ; then
    rapatrie $AHT3D  $P_I_DIR $F_DTA_DIR $NEMO_AHT3D
fi

if [ $DYNCOEF2D = 1 ] ; then
    rapatrie $AHM2D  $P_I_DIR $F_DTA_DIR $NEMO_AHM2D
fi

if [ $DYNCOEF3D = 1 ] ; then
    rapatrie $AHM3D  $P_I_DIR $F_DTA_DIR $NEMO_AHM3D
fi



echo ' [3.2] : Initial conditions, damping files'
echo ' =========================================='

getinitdmp 

## Float initial position and restart float
if [ $IFLOAT = 1 ] ;  then
    rapatrie $FLOATFIL  $P_I_DIR $F_DTA_DIR  init_float
fi

## Ice initial condition only if no = 1
if [ $ICE_INI = 1 -a $no -eq 1 ] ; then
    tmp=$(LookInNamelist ln_iceini_file namelist_ice namini ) ; tmp=$(normalize $tmp)
    if [ $tmp = T ] ; then 
       geticeini
    fi
fi

## Ice damping file 
if [ $ICE_DMP = 1  ] ; then
    geticedmp
fi

echo ' [3.3] : Forcing fields'
echo ' ======================'
getforcing  # this function read the namelist and look for the name of the
            # forcing files to be fetched. 
            # It also get the weight files and/or katabatic mask (if any)
            # the runoff files and SSS/SST restoring files (as appropriate)

## SHLAT2D file
getshlat2d  # check namelist for ln_shlat2d. If true, get the file specified
            # in the namelist block namlbc, sn_shlat2s%clname 

## Use a climatology of SSS damping in order to correct E-P
if [ $WAFDMP = 1 ] ; then
    rapatrie $WAFDMP_CLIM $P_I_DIR $F_DTA_DIR $NEMO_WAFDMP_CLIM
fi

## Open boundaries files : there are no files for agrif nest -> core_rapatrie

echo ' [3.4] : BDY files '
echo '  ===================='
if [ $BDY = 1 ] ; then
    getbdy   # This function look for information in the namelist on how many bdy there are.
             # Then it looks for the name of the file to be imported and rapatrie the set of files
             # needed by the current segment, in the TMPDIR directory.
fi

echo ' [3.5] : restart files'
echo '  ===================='
prev_ext=$(( $no - 1 ))      # file extension of previous run

## model restarts
if [ $no -eq  1 ] ; then
    echo '   ***  Cold start, no restart to fetch.'
else
    proc0=0
    procn=$(( (NB_NPROC / ENSEMBLE_SIZE) - 1  ))
    filext='nc'
    echo "   ***  $filext restart files from $proc0 to $procn"
    format='%04d'
    if [ $procn -ge 10000  ] ; then format='%05d'  ; fi
    if [ $procn -ge 100000 ] ; then format='%06d'  ; fi  # in order to run nemo on more than 100 k cores
    
    tmp=$(LookInNamelist ln_ens_rst_in namelist ) ; ln_ens_rst_in=$( normalize $tmp )
        if [ $ln_ens_rst_in = F ] ; then
            ZENSEMBLE_START=-1 ; ZENSEMBLE_END=-1    # reduce RESTART loop to 1
        else
            ZENSEMBLE_START=$ENSEMBLE_START ; ZENSEMBLE_END=$ENSEMBLE_END
        fi

        for member in $(seq $ZENSEMBLE_START  $ZENSEMBLE_END ) ; do
            mmm=$(getmember_extension $member)          # if no ensemble, this function return empty string
            nnn=$(getmember_extension $member  nodot )  # number of the member without .

            zrstdir='./'   # set zrstdir to RST directory if used
            if [ $RST_DIRS = 1 ] ; then zrstdir=$DDIR/${CN_DIRRST}.$prev_ext/$nnn/ ; fi

            if [  $ln_ens_rst_in = F ] ; then mmm='' ; fi # force mmm to empty string if ln_ens_rst_in = false

            if [ $RST_READY = 1 ] ; then
                OCE_RST_IN=$(LookInNamelist cn_ocerst_in )-${prev_ext}$mmm
                ICE_RST_IN=$(LookInNamelist cn_icerst_in namelist_ice)-${prev_ext}$mmm
                TRC_RST_IN=$(LookInNamelist cn_trcrst_in namelist_top)-${prev_ext}$mmm
                TRD_RST_IN=$(LookInNamelist cn_trdrst_in )-${prev_ext}$mmm
                STO_RST_IN=$(LookInNamelist cn_storst_in )-${prev_ext}$mmm
                ICB_RST_IN=$(LookInNamelist cn_icbrst_in )-${prev_ext}$mmm
            else
##### O C E A N
###############
                OCE_RST_IN=$(LookInNamelist cn_ocerst_in )$mmm
                OCE_RST_OUT=$(LookInNamelist cn_ocerst_out )$mmm
                ok=1

                if [ ! $RST_SKIP  ] ; then 
         # ** check if the nc  restart file are in the current dir 
         #   Look for restart_xxxx.nc.$prev_ext  xxx=1,$procn
                    for (( proc=$proc0 ; proc <= $procn ; proc++ )) ; do
                        rest=$(printf "${OCE_RST_IN}_${format}.${filext}.$prev_ext\n" $proc)
                        if [ ! -f ${zrstdir}$rest ] ; then ok=0 ; fi
                    done
                fi
                
                if [ $ok -eq 1 ] ; then
                    echo "   ***  All ocean restart files are here for member  $mmm " 
                else
      # look for tar restart files
                    mkdir -p $zrstdir   # just in case zrstdir can be ./ : no harm !

                    for rest in $( lsrestart ${OCE_RST_OUT}_oce_v2.$prev_ext.tar. ) ; do
                        rapatrie $rest   $P_R_DIR  $F_R_DIR  ${zrstdir}$rest
                        cd ${zrstdir} ; tar xf $rest ; cd -   #  also work with zrstdir = ./
                    done
                fi

      # link  the xxx.${filext}.$prev_ext files to xxx.${filext}  (no extension).
                cd $zrstdir
                for rest in ${OCE_RST_IN}_[[:digit:]]*[[:digit:]].${filext}.$prev_ext ; do
                    ln -sf $rest  ${rest%.$prev_ext} 
                done

                if [ $AGRIF = 1 ] ; then
                    for idx in ${agrif_pref[@]} ; do
                        for rest in ${idx}_${OCE_RST_IN}_[[:digit:]]*[[:digit:]].${filext}.$prev_ext ; do
                            ln -sf $rest  ${rest%.$prev_ext} 
                        done
                    done
                fi
                cd -
            
##### I C E
###########
            if [ $ICE = 1 ] ; then
        # ** check if the ${filext} restart file are in the current dir 
                ICE_RST_IN=$(LookInNamelist cn_icerst_in namelist_ice)$mmm
                ICE_RST_OUT=$(LookInNamelist cn_icerst_out namelist_ice )$mmm
                ok=1
                if [ ! $RST_SKIP  ] ; then 
          #   Look for restart_ice_inxxxx.${filext}.$prev_ext  xxx=1,$procn
                     for (( proc=$proc0 ; proc <= $procn ; proc++ )) ; do
                         rest=$(printf "${ICE_RST_IN}_${format}.${filext}.$prev_ext\n" $proc)
                         if [ ! -f ${zrstdir}$rest ] ; then ok=0 ; fi
                     done
                 fi

                 if [ $ok -eq 1 ] ; then
                     echo "   ***  All ice restart files are here for member  $mmm " 
                 else
          # look for tar restart files
                     mkdir -p $zrstdir   # just in case zrstdir can be ./ : no harm !

                     for rest in $( lsrestart ${ICE_RST_OUT}_v2.$prev_ext.tar. ) ; do
                         rapatrie $rest   $P_R_DIR  $F_R_DIR  ${zrstdir}$rest
                         cd ${zrstdir} ; tar xf $rest ; cd -   #  also work with zrstdir = ./
                     done
                 fi

      # link  the xxx.${filext}.$prev_ext files to xxx.${filext}  (no extension).
                 cd $zrstdir
                 for rest in ${ICE_RST_IN}_[[:digit:]]*[[:digit:]].${filext}.$prev_ext ; do
                     ln -sf $rest  ${rest%.$prev_ext} 
                 done

                 if [ $AGRIF = 1 ] ; then
                     for idx in ${agrif_pref[@]} ; do
                         for rest in ${idx}_${ICE_RST_IN}_[[:digit:]]*[[:digit:]].${filext}.$prev_ext ; do
                             ln -sf $rest  ${rest%.$prev_ext} 
                         done
                     done
                 fi
                 cd -
            fi   # ICE = 1
            
##### T R A C E R S 
################### 
# Missing option RST_DIRS for TRACERS in case of global restart ...
            if [ $TOP = 1 ] ; then
                rsttrc=$(LookInNamelist lrsttr namelist_top) 

       # test if lrsttr is true or false. If true, then next line return a 1 status ($? )
                echo $rsttrc | grep -q false
                if [ $? = 1 ] ; then
                    TRC_RST_IN=$(LookInNamelist cn_trcrst_in namelist_top)$mmm
                    TRC_RST_OUT=$(LookInNamelist cn_trcrst_out namelist_top)$mmm
          # 1) look for restart.nc in P_R_DIR ( case of restart from nc file -- change number of proc, for example-- )
                    if [ -f $P_R_DIR/${TRC_RST_IN}.nc.$prev_ext ] ; then
                        ln -s $P_R_DIR/${TRC_RST_IN}.nc.$prev_ext $TRC_RST_IN.nc
                    else
                        ok=1
                        if [ ! $RST_SKIP  ] ; then 
                            for (( proc=$proc0 ; proc <= $procn ; proc++ )) ; do
                                rest=$(printf "${TRC_RST_IN}_${format}.${filext}.$prev_ext\n" $proc)
                                if [ ! -f ${zrstdir}$rest ] ; then ok=0 ; fi
                            done
                        fi

                        if [ $ok -eq 1 ] ; then
                            echo "   ***  All top restart files are here for member  $mmm " 
                        else

           # 2) standard procedure: look for tar files
                            for rest in $(lsrestart  ${TRC_RST_OUT}_v2.$prev_ext.tar. )  ; do
                                rapatrie_res $rest   $P_R_DIR  $F_R_DIR  ${zrstdir}restart.tar
                                cd $zrstdir ; tar xf restart.tar ; cd -
                            done
                        fi

           # restart files are archived with the correct TRC_RST_IN prefix !
                        cd $zrstdir
                        for rest in ${TRC_RST_IN}_[[:digit:]]*[[:digit:]].${filext}.$prev_ext ; do
                            ln -sf $rest  ${rest%.$prev_ext} 
                        done
                        cd -
                    fi
                fi
            fi
            
##### T R D  M L D
##################
            if [ $TRDMLD = 1  ] ; then
                trdmldrst=$(LookInNamelist ln_trdmld_restart )

         # test if trdmldrst  is true or false. If true, then next line return a 1 status ($? )
                echo $trdmldrst | grep -q false
                if [ $? = 1 ] ; then
                    TRD_RST_IN=$(LookInNamelist cn_trdrst_in )$mmm
                    TRD_RST_OUT=$(LookInNamelist cn_trdrst_out )$mmm
         # ** check if the ${filext} restart file are in the current dir
                    ok=1
                    if [ ! $RST_SKIP  ] ; then 
         #   Look for restart_mld_xxxx.${filext}.$prev_ext  xxx=1,$procn
                        for (( proc=$proc0 ; proc <= $procn ; proc++ )) ; do
                            rest=$(printf "${TRD_RST_IN}_${format}.${filext}.$prev_ext\n" $proc)
                            if [ ! -f ${zrstdir}$rest ] ; then ok=0 ; fi
                        done
                    fi

                    if [ $ok -eq 1 ] ; then
                        echo "   ***  All trdmld restart files are here for member  $mmm " 
                    else
                        for rest in $( lsrestart ${TRD_RST_OUT}_v2.$prev_ext.tar.) ; do
                            rapatrie $rest   $P_R_DIR  $F_R_DIR  ${zrstdir}$rest
                            cd $zrstdir ; tar xf $rest ; cd -
                        done
                    fi
         # remove the prev_ext from the name file
                    cd $zrstdir
                    for rest in ${TRD_RST_IN}_[[:digit:]]*[[:digit:]].${filext}.$prev_ext ; do
                        ln -sf $rest ${rest%.$prev_ext}
                    done
                    cd -
                fi
            fi

##### S T O 
###########
            if [ $RSTO = 1 ] ; then
                STO_RST_IN=$(LookInNamelist cn_storst_in)$mmm 
                STO_RST_OUT=$(LookInNamelist cn_storst_out)$mmm 

        # ** check if the ${filext} restart file are in the current dir
                ok=1
                if [ ! $RST_SKIP  ] ; then 
           #   Look for restart_mld_xxxx.${filext}.$prev_ext  xxx=1,$procn
                    for (( proc=$proc0 ; proc <= $procn ; proc++ )) ; do
                        rest=$(printf "${STO_RST_IN}_${format}.${filext}.$prev_ext\n" $proc)
                        if [ ! -f ${zrstdir}$rest ] ; then ok=0 ; fi
                    done
                fi

                if [ $ok -eq 1 ] ; then
                    echo "   ***  All sto restart files are here for member  $mmm " 
                else
                    for rest in $( lsrestart ${STO_RST_OUT}_v2.$prev_ext.tar.) ; do
                        rapatrie $rest   $P_R_DIR  $F_R_DIR  ${zrstdir}$rest
                        cd $zrstdir ; tar xf $rest ; cd -
                    done
                fi
        # remove the prev_ext from the name file
                cd $zrstdir
                for rest in ${STO_RST_IN}_[[:digit:]]*[[:digit:]].${filext}.$prev_ext ; do
                    ln -sf $rest ${rest%.$prev_ext}
                done
                cd -
            fi

##### I C B
###########
            if [ $ICB = 1 ] ; then  # need iceberg restart
                ICB_RST_IN=$(LookInNamelist cn_iscbrst_in  )$mmm
                ICB_RST_OUT=$(LookInNamelist cn_icbrst_out )$mmm


        # ** check if the ${filext} restart file are in the current dir
                ok=1
                if [ ! $RST_SKIP  ] ; then 
          #   Look for restart_icebergs_xxxx.${filext}.$prev_ext  xxx=1,$procn
                    for (( proc=$proc0 ; proc <= $procn ; proc++ )) ; do
                        rest=$(printf "${ICB_RST_IN}_${format}.${filext}.$prev_ext\n" $proc)
                        if [ ! -f ${zrstdir}$rest ] ; then ok=0 ; fi
                    done
                fi

                if [ $ok -eq 1 ] ; then
                    echo "   ***  All icb restart files are here for member  $mmm " 
                else
                    for rest in $( lsrestart ${ICB_RST_OUT}_v2.$prev_ext.tar.) ; do
                        rapatrie $rest   $P_R_DIR  $F_R_DIR  ${zrstdir}$rest
                        cd $zrstdir ; tar xf $rest ; cd -
                    done
                fi
        # remove the prev_ext from the name file
                cd $zrstdir
                for rest in ${ICB_RST_IN}_[[:digit:]]*[[:digit:]].${filext}.$prev_ext ; do
                    ln -sf $rest ${rest%.$prev_ext}
                done
                cd -
            fi

            fi  # RST_READY

        done   # loop on member
fi       # test on restart

## Float initial position and restart float
if [ $RFLOAT = 1 -a  $IFLOAT = 1 ] ;  then
    rapatrie restart_float.$prev_ext $P_R_DIR  $F_R_DIR restart_float 
fi

pwd
touch donecopy

echo '(4) Run the code'
echo '----------------'
date
if [ $XIOS = 1 -a $NB_NPROC_IOS != 0 ] ; then
   NB_NCORE_DP=${NB_NCORE_DP:=0}
   if [ $NB_NCORE_DP != 0 ] ; then
      runcode_mpmd_dp  -dp $NB_NCORE_DP $NB_NPROC ./nemo4.exe $NB_NPROC_IOS ./xios_server.exe
   else
      runcode_mpmd  $NB_NPROC ./nemo4.exe $NB_NPROC_IOS ./xios_server.exe
   fi
else
    runcode  $NB_NPROC ./nemo4.exe
fi
date
exit  # Scalability tests stop here
#--------------------------------------------------------
#########################################################
#--------------------------------------------------------
echo '(5) Post processing of the run'
echo '------------------------------'
 # gives the rights r to go
chmod -R go+r  $TMPDIR
cp layout.dat $P_S_DIR/ANNEX/
  ntiming=0 ; nmsh=0
  # check flags in namelist for further processing (before changing namelist)
  tmp=$(LookInNamelist ln_timing ) ; tmp=$( normalize $tmp )
  if [ $tmp = 'T' ] ; then ntiming=1 ; fi  
  tmp=$(LookInNamelist ln_meshmask ) ; tmp=$( normalize $tmp )
  if [ $tmp = 'T' ] ; then nmsh=1 ; fi  

echo ' [5.1] check the status of the run'
echo ' ================================'
  # touch OK.member to mark successfull members
nok=0
for member in $(seq $ENSEMBLE_START $ENSEMBLE_END) ; do
    mmm=$(getmember_extension $member)
    outfile=ocean.output$mmm
    # should be modified in NEMO to make it safer ( eg. write(numout,*) ' Run completed successfully'  )
    if [ "$(tail -20 $outfile | grep AAAAAAAA)" = 'AAAAAAAA' ] ; then touch OK$mmm ; nok=$(( nok + 1 )); fi

done

  # check nok vs ensemble size
STOP_FLAG=0    # 0 : run is clean,  1 : at least 1 member failed,   2: at least one member failed and no step done
if [[ $nok < $ENSEMBLE_SIZE ]] ; then 
    echo "   ---  WARNING: Some members failed, we set a stop flag"
    STOP_FLAG=1
fi

for member in $(seq $ENSEMBLE_START $ENSEMBLE_END) ; do
    mmm=$(getmember_extension $member)
    if [ ! -f time.step$mmm ] ; then
        echo "   ---  WARNING: No time-step are made by NEMO, for member  $member  Earlier crash."   
        touch STOP$mmm
        STOP_FLAG=2
    fi 
done

  # Post process the run according to the STOP_FLAG
case $STOP_FLAG in
    ( 0 )
    date
    ext=$no
    echo "   ***  Run OK"
    echo ' [5.2] rename  the restart files'
    echo ' ==============================='
    if [ $RST_READY = 1 ] ; then
       echo "   *** Restart files are ready from NEMO ..."
    else
    for member in  $(seq $ENSEMBLE_START $ENSEMBLE_END) ; do
        mmm=$(getmember_extension $member)
        nnn=$(getmember_extension $member nodot)
       # clean previous restart file from the DDIR/xxx-R directory
        clean_res $prev_ext

       # O C E A N 
       # *********
        OCE_RST_IN=$(LookInNamelist cn_ocerst_in )
        OCE_RST_OUT=$(LookInNamelist cn_ocerst_out )

        renamerst  $OCE_RST_IN $OCE_RST_OUT 

       # I C E 
       # *****
        if [ $ICE = 1 ] ; then
            ICE_RST_IN=$(LookInNamelist cn_icerst_in namelist_ice)
            ICE_RST_OUT=$(LookInNamelist cn_icerst_out namelist_ice )

          renamerst  $ICE_RST_IN $ICE_RST_OUT 
        fi

       # P A S S I V E   T R A C E R S
       # *****************************
        if [ $TOP = 1 ] ; then
            TRC_RST_IN=$(LookInNamelist cn_trcrst_in namelist_top)
            TRC_RST_OUT=$(LookInNamelist cn_trcrst_out namelist_top)

          renamerst  $TRC_RST_IN $TRC_RST_OUT 
        fi

       # T R D  M L D 
       # ************
        if [ $TRDMLD = 1 ] ; then
            TRD_RST_IN=$(LookInNamelist cn_trdrst_in )
            TRD_RST_OUT=$(LookInNamelist cn_trdrst_out )

            renamerst  $TRD_RST_IN $TRD_RST_OUT 
        fi
        
       # S T O 
       # *****
        if [ $STO = 1 ] ; then
            STO_RST_IN=$(LookInNamelist cn_storst_in)
            STO_RST_OUT=$(LookInNamelist cn_storst_out)

            renamerst  $STO_RST_IN $STO_RST_OUT 
        fi

       # I C B
       # *****
        if [ $ICB = 1 ] ; then
            ICB_RST_IN=restart_icebergs
            ICB_RST_OUT=icebergs_restart

#       ICB restart file are not named as other restart files ( with nitend in the name). Som renamerst does not work
#       till some update in the icb code.  
#        renamerst  $ICB_RST_IN $ICB_RST_OUT
            for f in ${ICB_RST_OUT}* ; do 
                g=$( echo $f | sed -e "s/$ICB_RST_OUT/$ICB_RST_IN/").$ext
                mv $f $g
            done
        fi

    done      # loop on members
    fi   # RST_READY

    date
    echo ' [5.3] Update the CONFIG_CASE.db file'
    echo ' ===================================='
    mmm=$(getmember_extension $ENSEMBLE_START) # take the first member(if any) for reference
    output_ref=ocean.output$mmm  # used in update_db_file function to infer the ending date
    update_db_file 
    cat $CONFIG_CASE.db
    copy $CONFIG_CASE.db $P_CTL_DIR/ 

    date
    echo ' [5.4] Rename namelists, ocean.output and other text files. Copy to P_S_DIR'
    echo ' =========================================================================='
    rename_txt_files  $ext
      # from now, take care of using the correct namelist name in the calls !

    date
    echo ' [5.5] Make restart tar files '
    echo ' ============================='
     # Build a script (to be submitted) for saving the individual ${filext}.$ext 
     # restart files into a set of tar files and expatrie_res them.
    mksavrst  zsrst.$ext.sh   

     # Submit the save-restart script 
     # When this script is finished ( asynchronously), there is a touch statement on file RST_DONE$mmm.$ext,
     # that need to be checked before cleaning. 
    submit ${P_CTL_DIR}/zsrst.$ext.sh
    cd $TMPDIR   # back in TMPDIR for sure

    date
    echo ' [5.6] Ready to re-submit the job NOW (to take place in the queue)'
    echo ' ================================================================='
    TESTSUB=$( wc $CONFIG_CASE.db | awk '{print $1}' )
    if [ $TESTSUB -le  $MAXSUB -o -f  FORCE_RESUB ] ; then
        submit  ${P_CTL_DIR}/${SUBMIT_SCRIPT}
        cd $TMPDIR
        cat $TMPDIR/logsubmit
    else
        echo "   --- WARNING: Maximum auto re-submit reached."
    fi ;;

    ( 1 | 2 )   # The run crashed :( . Send back some infos on the CTL directory
    date
    echo "   ---  WARNING: Run crashed."
    ext=$$.'ABORT'
    rename_txt_files $ext $P_CTL_DIR  # rename text files and copy ocean.output and run.stat to P_CTL_DIR
      # from now, take care of using the correct namelist name in the calls !

     # rename the output directory extension in order not to mix the output with correct run
    echo "   ***  ${CONFIG_CASE}-${DIROUText}.$no renamed to ${CONFIG_CASE}-${DIROUText}.$ext"
    mv $DDIR/${CONFIG_CASE}-${DIROUText}.$no $DDIR/${CONFIG_CASE}-${DIROUText}.$ext 
    if [ $DIAOBS = 1 ] ; then
        echo "   ***  ${CN_DIAOBS}.$no renamed to ${CN_DIAOBS}.$ext"
        mv $DDIR/${CN_DIAOBS}.$no $DDIR/${CN_DIAOBS}.$ext 
    fi ;;
esac

cd $TMPDIR
echo '(6) : Final Post Processing after next run is queued'
echo '----------------------------------------------------'
date

case $STOP_FLAG in
    ( 0 | 1 )  # run is OK or at least some step were performed, build the nc output files
    if [ $STOP_FLAG = 0 ] ; then ext=$no     ; fi
    if [ $STOP_FLAG = 1 ] ; then ext=$$.'ABORT' ; fi


    if [ $XIOS = 1 ] ; then
        echo ' [6.1] Process the rebuild of nc file from XIOS files '
        echo ' ========================================================='

        cp $CN_DOMCFG  $DDIR/${CONFIG_CASE}-${DIROUText}.$ext/
        cp iodef.xml domain_def.xml $DDIR/${CONFIG_CASE}-${DIROUText}.$ext/
   
        if [ $DIAPTR = 1 ] ; then  # process diaptr files (one_file mode). Only rename
        date
        echo ' [6.1.1] Process the rename of diaptr nc file from XIOS files '
        echo ' ========================================================='
            for member in $(seq  $ENSEMBLE_START $ENSEMBLE_END ) ; do
              nnn=$(getmember_extension $member nodot)
              mmm=$(getmember_extension $member      )
              cd  $DDIR/${CONFIG_CASE}-${DIROUText}.$ext/$nnn
              for f in ${CONFIG_CASE}${mmm}*diaptr_*.nc ; do
                 if [ -f $f ] ; then  # check if $f exist ( case of ln_diaptr=T but no i/o in iodef)
                   freq=$( echo $f | awk -F_ '{print $2}' )
                   tmp=${f##*-} ; indastp=${tmp%.nc}
                   tag=y${indastp:0:4}m${indastp:4:2}d${indastp:6:2}
                   g=${CONFIG_CASE}${mmm}_${tag}.${freq}_diaptr.nc
                   mkdir -p ${freq}_OUTPUT
                   mv $f ${freq}_OUTPUT/$g
                 fi
              done
            done
        fi
        cd $TMPDIR

        date
        echo ' [6.1.2] Process the rename of zoomed  nc file from XIOS files '
        echo ' ========================================================='
#      # look for zoom (in one_file mode). Assume domain_ref without 'grid' keyword
        if [ -f 06-secme.xml ] ; then    # case of ensemble runs : no need to scan the full iodef.xml (huge!)
           cat 0*sec*xml > ziodef.xml
        else
           ln -sf iodef.xml ziodef.xml
        fi
        iodef=ziodef.xml
        lis_zoom_tmp=$( for f in $(  grep -i domain_ref $iodef | grep -vi grid_ ) ; do \
            cmd=$(echo $f | grep domain_ref ) ; \
            if [ $cmd ] ; then eval $cmd ; echo $domain_ref ;fi ;\
                   done )
        lis_zoom=$( for zoom in $lis_zoom_tmp ; do \
            echo $zoom ; done | sort -u )
       # check enabled zoom
        lis_zoom_enabled=''
        zoom_enable=''
        for zoom in $lis_zoom ; do
            grep -wB 1 $zoom $iodef | head -1 | grep -iq 'enabled=".true."'
            if [ $? = 0 ] ; then
                lis_zoom_enabled="$lis_zoom_enabled $zoom"
                zoom_enable=1
            fi
        done
        lis_zoomid=$( for f in $lis_zoom_enabled ; do n=$(( ${#f} - 1 )) ; echo ${f:0:$n} ; done | sort -u )

        if [ $zoom_enable ] ; then
          if [ ! -f $CTL_DIR/post_process_zoom.sh ] ; then
             mk_post_process post_process_zoom.sh
             cp post_process_zoom.sh $CTL_DIR
           else
             cp $CTL_DIR/post_process_zoom.sh ./
          fi
         
          if [ $ENSEMBLE = 1 ] ; then
# JMM CAUTION not a generic function for all machines ...
            cmd='mpirun --bynode'
            for n in $(seq $ENSEMBLE_START $ENSEMBLE_END ) ; do
              cmd="$cmd -np 1 ./post_process_zoom.sh $n $ext :"
            done
            cmd=${cmd%:}
            $cmd
        #   runcode_mpmd ./zpp.conf.$ext
          else
             ./post_process_zoom.sh -1 $ext
          fi
          cd $TMPDIR
        fi
        # look for one-file or multiple_file mode (for global files)
        # file_definition is defined only once (must be) in iodef.xml file
        # if zoom are used, one_file is forced at the level of the file object
        MULTIPLE_FILE=$( grep '<file_definition' iodef.xml | \
            awk '{ i=1 ; while (i <= NF ) { print $i ; i++ } }' | \
            grep type | grep -q 'one_file' ; echo $? )
        # also check if there is only 1 XIOS server ==> ONE_FILE
        if [ $NB_NPROC_IOS = 1 ] ; then MULTIPLE_FILE=0 ; fi
        if [ $MULTIPLE_FILE = 1 ] ; then
            if [ $MERGE = 0 ] ; then
                echo "   ***  Recombine for XIOS using rebuild_nemo in a batch"
                mkbuild_merge zmergxios.$ext.sh  
                submit ${P_CTL_DIR}/zmergxios.$ext.sh
            else  # MERGE on the fly 
                echo "   ***  Recombine for XIOS using mergefile_mpp2 on the fly"
                if [ $ENSEMBLE = 1 ] ; then 
                    mergefiles_ens
                else
                    mergefiles 
                fi
            fi
        else   # one_file output 
         # almost nothing to do except renaming and correction of nav_lon,nav_lat if necessary
            cd $DDIR/${CONFIG_CASE}-XIOS.$ext/
            post_process_one_file 
            cd $TMPDIR
        fi
        
    fi 

    if [ $nmsh  != 0 ] ; then 
        date
        echo ' [6.2] Process the rebuild of mesh_mask files'
        echo ' ============================================'
            echo '  ***  netcdf meshmask files'
          # To be done
            echo "  ---  WARNING: Missing script for this case ..."
    fi 

    if [ $STOP_FLAG = 0 -a $DIAOBS = 1 ] ; then
        date
        echo '   [6.2.2] Recombine and save OBS fdbk files '
        echo '   ------------------------------------------'
#        cp ./fbcomb.exe $DDIR/${CN_DIAOBS}.$ext
        cd $DDIR/${CN_DIAOBS}.$ext

        if [ $ENACT = 1 ] ; then
            mkdir -p  $P_S_DIR/OBS/
            cat << eof > ztarobs.$ext.sh
#!/bin/bash
mmm=\$1
tar cf $P_S_DIR/OBS/${CONFIG_CASE}\${mmm}_enact_fdbk.$ext.tar enact\${mmm}_fdbk_*.nc
eof
            chmod 755 ./ztarobs.$ext.sh
            cmd="mpirun --bynode "
            for member in  $(seq $ENSEMBLE_START $ENSEMBLE_END) ; do
                mmm=$(getmember_extension $member)
# JMM : can be // with runcode_mpmd (mergefile_mpp2 )
#               ./fbcomb.exe  enact${mmm}_fdbk.nc.$ext enact${mmm}_fdbk_*.nc
#               expatrie enact${mmm}_fdbk.nc.$ext $F_S_DIR ${CONFIG_CASE}${mmm}_enact_fdbk.nc.$ext
#               mv enact${mmm}_fdbk.nc.$ext $P_S_DIR/OBS/${CONFIG_CASE}${mmm}_enact_fdbk.nc.$ext
#                tar cf $P_S_DIR/OBS/${CONFIG_CASE}${mmm}_enact_fdbk.$ext.tar enact${mmm}_fdbk_*.nc
             cmd="$cmd -np 1 ./ztarobs.$ext.sh $mmm :"
            done
            cmd=${cmd%:}  # eliminate last :
            $cmd
        fi

        if [ $SLA = 1 ] ; then
            mkdir -p  $P_S_DIR/OBS/
            cat << eof > ztarobssla.$ext.sh
#!/bin/bash
mmm=\$1
tar cf $P_S_DIR/OBS/${CONFIG_CASE}\${mmm}_sla_fdbk.$ext.tar enact\${mmm}_fdbk_*.nc
eof
            chmod 755 ./ztarobssla.$ext.sh
            cmd="mpirun --bynode "
            for member in  $(seq $ENSEMBLE_START $ENSEMBLE_END) ; do
                mmm=$(getmember_extension $member)
# JMM : can be // with runcode_mpmd (mergefile_mpp2 )
#               ./fbcomb.exe  sla${mmm}_fdbk.nc.$ext sla${mmm}_fdbk_*.nc
#               expatrie sla${mmm}_fdbk.nc.$ext $F_S_DIR ${CONFIG_CASE}${mmm}_sla_fdbk.nc.$ext
#               mv sla${mmm}_fdbk.nc.$ext $P_S_DIR/OBS/${CONFIG_CASE}${mmm}_sla_fdbk.nc.$ext
#             tar cf $P_S_DIR/OBS/${CONFIG_CASE}${mmm}_sla_fdbk.$ext.tar sla${mmm}_fdbk_*.nc
             cmd="$cmd -np 1 ./ztarobssla.$ext.sh $mmm :"
            done
            cmd=${cmd%:}  # eliminate last :
            $cmd
        fi
        cd $TMPDIR  # back to TMPDIR
    fi

    date
    echo ' [6.3] Pack some files in annex tar file for archiving on F machine'
    echo ' =================================================================='
     # note : next tar command takes into account all members of an ensemble run and all AGRIF zoom envolved
    tar cf tarfile.${CONFIG_CASE}_annex.$ext *run.stat*.$ext *ocean.outpu*.$ext *namelist_oce.$ext
    if [ $ntiming = 1 ] ; then tar rf tarfile.${CONFIG_CASE}_annex.$ext *timing.outpu*.$ext ; fi
    if [ $ICE = 1     ] ; then tar rf tarfile.${CONFIG_CASE}_annex.$ext *namelist_ice.$ext  ; fi
    
    expatrie tarfile.${CONFIG_CASE}_annex.$ext  $F_S_DIR tarfile.${CONFIG_CASE}_annex.$ext

    date
    echo ' [6.4] Extra diags output (for memory) '
    echo ' ========================='
    if [ $IFLOAT = 1 ] ; then
        echo '   ***  Floats'
        expatrie trajec_float  $F_S_DIR ${CONFIG_CASE}_trajfloat.$ext
        expatrie restart_float.out  $F_R_DIR  restart_float.$ext
        if [ $RFLOAT = 0 ] ; then
            echo "Ne pas effacer ce fichier.MERCI." > $P_CTL_DIR/float
            chmod a-wx $P_CTL_DIR/float
        fi
    fi ;;

    ( 2 )
    date
    echo ' [6.0] No step performed exit then'
    echo '================================='
    echo "   ===  ERROR : final exit "
    exit ;;
esac

#########################################################################
##                                END                                  ##
#########################################################################
