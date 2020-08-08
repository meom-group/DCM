#!/bin/bash
# $Id: function_4_all.sh $
######################################################################################
#    bash script functions used below: valid for any machines 
######################################################################################
# VALID FOR ANY MACHINES
#########################
# LookInNamelist()
# keychk() 
# normalize()
# getblock()
# getinitdmp()
# geticedmp()
# getshlat2d()
# SetYears()
# getforcing()
# gettmx()
# getbdy()
# getweight()
# getfiles()
# getmember_extension()
# number_of_file()
# renamerst()
# mktarrst()
# mksavrst()
# mkbuildnc()
# mkbuild_5d()
# mkbuildnc_daily()
# mkbuild_ssf()
# mkbuild_mesh_mask()
# initagrif()
# update_db_file()
#-------------------------------------------------------

# LookInNamelist returns the value of a variable in the namelist
#        examples: aht0=$(LookInNamelist aht0 )  <=> aht0=$(LookInNamelist aht0 namelist )
#                  ln_limdmp=$(LookInNamelist ln_limdmp namelist_ice )
#                  nit000=$(LookInNamelist nn_it000 namelist_oce.10 ) 
#        If there is a third argument it is used as a namelist block and the search is
#        limited to this block :
#                  ln_tsd_init=$(LookInNamelist ln_tsd_init namelist_cfg namtsd_drk )
LookInNamelist()    {
         if [ $# -ge 2 ] ; then znamelist=$2 ; else znamelist=namelist ; fi
         if [ $# = 3   ] ; then zblk=$3      ; else zblk=''            ; fi
         if [ ! $zblk ] ; then
           eval grep -e $1 $znamelist      | tr -d \' | tr -d \"  | sed -e 's/=/  = /' | awk ' {if ( $1 == str ) print $3 }' str=$1
         else
          getblock $zblk $znamelist | eval grep -e $1  | tr -d \' | tr -d \"  | sed -e 's/=/  = /' | awk ' {if ( $1 == str ) print $3 }' str=$1
         fi
                    }
# --- 

# keychk look for the presence of cpp key in CPP.keys file. Lines starting with #
#    are not take into account. If key found, return OK else keychk is unset
keychk() {
       cat CPP.keys | grep -v '^#' | grep -q $1
       if [ $? = 0 ] ; then
         echo OK
       fi
          }
# ---

# For logical value in namelist always return T or F despite the namelist format ( TRUE,true, true etc ...)
normalize()         {
               tmp=$(echo $1 | tr 'a-z' 'A-Z' )
               echo $tmp  | grep -q 'T'
               if [ $? = 0 ] ; then echo T ; else echo F ; fi
                    }
# ---

# Get a namelist block from its name in namelist
getblock()          { 
            # if a 2nd argument is passed, it is a namelist name. Default to 'namelist'
            if [ $2 ] ; then namelist=$2 ; else namelist=namelist ; fi
            cat $namelist | awk 'BEGIN {flip=0} { \
            if ( $1 == "&"blk && flip == 0 ) { flip=1   }  \
            if ( $1 != "/"  && flip == 1   ) { print $0 }  \
            if ( $1 == "/"  && flip == 1   ) { print $0 ; flip=0 }    \
                                    }' blk=$1
                    }
# ---
# Get initial forcing and/or damping files if required.
getinitdmp()        {
     # initial condition
     SetYears
     tmp=$(LookInNamelist ln_tsd_init namelist namtsd_drk) ; tmp=$(normalize $tmp )
     if [ $tmp = T ] ; then 
       filter='| grep -v sn_tem_dmp | grep -v sn_sal_dmp'  # at this level always filter dmp files
       blk=namtsd_drk ;  getfiles $blk $P_DTA_DIR $F_DTA_DIR   
                         getweight $blk $P_WEI_DIR $F_WEI_DIR
     fi
     # damping files (TRADMP is set in nemo4.sh)
     if [ $TRADMP = 1 ] ; then
        tmp=$(LookInNamelist ln_tsd_dmp namelist namtsd_drk) ; tmp=$(normalize $tmp )
        if [ $tmp = T ] ; then 
          filter='| grep -v sn_tem_ini | grep -v sn_sal_ini'  # at this level always filter ini files
          blk=namtsd_drk ;  getfiles  $blk $P_DTA_DIR $F_DTA_DIR 
                            getweight $blk $P_WEI_DIR $F_WEI_DIR
          # get also resto file ( in case of std NEMO stuff since 3.6 )
          nn_hdmp=$(LookInNamelist nn_hdmp namelist namtra_dmp_drk) 
          if [ $nn_hdmp != -2 ] ; then
             cn_resto=$(LookInNamelist cn_resto namelist namtra_dmp )
             rapatrie $cn_resto  $P_DTA_DIR $F_DTA_DIR $cn_resto
          else
             # look for ln_dmpmask
             tmp=$(LookInNamelist ln_dmpmask namelist namtra_dmp_drk) ; tmp=$(normalize $tmp )
             if [ $tmp = T ] ; then
                blk=namtra_dmp_drk ; getfiles  $blk $P_DTA_DIR $F_DTA_DIR   
             fi
          fi
        fi
     fi
                    }
# ---
# Get ice initialisation if required
geticeini()        {
    tmp=$(LookInNamelist ln_iceini_file namelist_ice namini ) ; tmp=$(normalize $tmp)
    if [ $tmp = T ] ; then
      filter=''
      blk=namini ;  getfiles $blk $P_DTA_DIR $F_DTA_DIR namelist_ice
    fi
                   }
# ---
# Get ice damping files if required.
geticedmp()        {
     tmp=$(LookInNamelist ln_limdmp namelist_ice) ; tmp=$(normalize $tmp )
     if [ $tmp = T ] ; then
       SetYears
       filter=''  # take all sn_xxx of the block
       blk=namice_dmp ;  getfiles  $blk $P_OBC_DIR $F_OBC_DIR namelist_ice
                         getweight $blk $P_WEI_DIR $F_WEI_DIR namelist_ice

     fi
                    }

# ---
# Get shlat2d file if required in namlbc namelist block
getshlat2d()        {
     filter=''
     tmp=$(LookInNamelist ln_shlat2d namelist namlbc_drk) ; tmp=$(normalize $tmp )
     if [ $tmp = T ] ; then
       blk=namlbc_drk ;  getfiles $blk $P_DTA_DIR $F_DTA_DIR
     fi
                    }
# ---
SetYears()          {
     dbfile=${CONFIG_CASE}.db
     # look for years being used :
     job_number=$( tail -1 $dbfile | awk '{print $1}' )
     # determine the number of years to get according to the lenght of the run
     zn1=$(LookInNamelist nn_it000 namelist )
     zn2=$(LookInNamelist nn_itend namelist )
     zrdt=$(LookInNamelist rn_rdt  namelist )
     zstpday=$( echo $zrdt | awk '{print 86400./$1 }' )
     znyear=$( echo $zn1 $zn2 $zstpday | awk '{ print int(( $2 - $1 +1)/$3/365+0.5 )}')
     if [ $znyear = 0 ] ; then znyear=1 ; fi   # force znyear to be at least 1

     if [ $job_number = 1 ] ; then    # case of cold start
       year=$( LookInNamelist nn_date0 namelist | awk '{print int($1/10000) }' )
       yearm1=$(( year - 1 )) 
       yearp1=$(( year + $znyear ))
     else
       yearm1=$( tail -2 $dbfile | head -1 | awk '{ print int($NF/10000) }' )
       yearm1=$(( yearm1 - 1 ))
       yearp1=$(( yearm1 + $znyear + 1 ))
     fi

       yearm1=$( printf "%04d" $yearm1 )
       yearp1=$( printf "%04d" $yearp1 )
                    }
# ---
# Getforcing : this function read information from the namelist and db file in order to know the 
#  forcing files to get, and get them.  No argument used
getforcing()        {
     SetYears
     
     # Look for type of forcing
     forcing=none
     tmp=$( LookInNamelist ln_usr      ) ; tmp=$( normalize $tmp ) ; if [ $tmp = T ] ; then forcing=ANA  ; blk=namsbc_ana  ; fi
### @@@@@@@@@  WIP for ANA : need to figura out usrdef_fbc
     tmp=$( LookInNamelist ln_flx      ) ; tmp=$( normalize $tmp ) ; if [ $tmp = T ] ; then forcing=FLX  ; blk=namsbc_flx  ; fi
     tmp=$( LookInNamelist ln_blk      ) ; tmp=$( normalize $tmp ) ; if [ $tmp = T ] ; then forcing=BLK  ; blk=namsbc_blk  ; fi
     tmp=$( LookInNamelist ln_cpl      ) ; tmp=$( normalize $tmp ) ; if [ $tmp = T ] ; then forcing=CPL  ; blk=namsbc_cpl  ; fi

     # get files corresponding the type of forcing
     case $forcing in 
     BLK ) 
#         getweight $blk $P_WEI_DIR $F_WEI_DIR

         echo  required forcing files :
         echo  =======================
         # check for optional files and set filter
         filter=''  
         tmp=$( LookInNamelist ln_taudif );  tmp=$(normalize $tmp )
         if [ $tmp = F ] ; then filter="$filter | grep -v sn_tdif " ; fi
         tmp=$( LookInNamelist ln_clim_forcing );  tmp=$(normalize $tmp )
         if [ $tmp = T ] ; then filter="$filter | grep -v sn_wndi  | grep -v sn_wndj " ; fi
         ln_clim_forcing=$tmp

         getweight $blk $P_WEI_DIR $F_WEI_DIR
         getfiles $blk  $P_FOR_DIR $F_FOR_DIR

         if [ $ln_clim_forcing = T ] ; then
           filter="$filter | grep -v sn_kati | grep -v sn_katj"
           blk_clim=namsbc_blk_drk
           getweight $blk_clim $P_WEI_DIR $F_WEI_DIR
           getfiles $blk_clim  $P_FOR_DIR $F_FOR_DIR
         fi ;;
    

      FLX )  
         getweight $blk $P_WEI_DIR $F_WEI_DIR
         getfiles $blk  $P_FOR_DIR $F_FOR_DIR ;;

     CPL ) ;;

     ANA ) ;;

     none )
         echo function getforcing only available for BLK FLX CPL or ANA  forcing
         exit 1 ;;
     esac

     # Extra forcing files 
     if [ $forcing = BLK ] ; then
        # KATABATIC mask
        filter=''
        tmp=$( LookInNamelist ln_kata namelist namsbc_blk_drk)  ;  tmp=$(normalize $tmp )
        if [ $tmp = T ] ; then 
          filter="$filter | grep -v sn_wmod | grep -v sn_uw | grep -v sn_vw"
          getfiles namsbc_blk_drk  $P_DTA_DIR $F_DTA_DIR 
        fi
     fi
   
     # Atmospheric pressure  # in BLK sea level pressure is  now required. So why 2 entries for atm pressure ?
     filter=''
     tmp=$(LookInNamelist ln_apr_dyn namelist) ; tmp=$(normalize $tmp )
     if [ $tmp = T ] ; then blk=namsbc_apr ;  getfiles $blk $P_DTA_DIR $F_DTA_DIR ;  fi

     # Drag coefficient from wave model 
     filter=''
     tmp=$(LookInNamelist ln_cdgw namelist) ; tmp=$(normalize $tmp )
     if [ $tmp = T ] ; then blk=namsbc_wave ;  getfiles $blk $P_DTA_DIR $F_DTA_DIR ;  fi

     # RUN_OFF 
     filter='| grep -v sn_s_rnf | grep -v sn_t_rnf | grep -v sn_dep_rnf '
     tmp=$(LookInNamelist ln_rnf namelist) ; tmp=$(normalize $tmp )
     if [ $tmp = T ] ; then blk=namsbc_rnf ;  getfiles $blk $P_DTA_DIR $F_DTA_DIR ;  fi

     # extra files 
     filter='| grep -v sn_rnf | grep -v sn_cnf '
     extra=0
     #   rnf depth file
     tmp=$(LookInNamelist ln_rnf_depth namelist) ; tmp=$(normalize $tmp )
     if [ $tmp = F ] ; then 
       filter="$filter | grep -v sn_dep_rnf "
     else 
       extra=1
     fi 
     #   rnf temperature file
     tmp=$(LookInNamelist ln_rnf_tem namelist) ; tmp=$(normalize $tmp )
     if [ $tmp = F ] ; then 
       filter="$filter | grep -v sn_t_rnf "
     else 
       extra=1
     fi 
     #   rnf salinity file
     tmp=$(LookInNamelist ln_rnf_sal namelist) ; tmp=$(normalize $tmp )
     if [ $tmp = F ] ; then 
       filter="$filter | grep -v sn_s_rnf "
     else 
       extra=1
     fi 

     if [ $extra = 1 ] ; then
       blk=namsbc_rnf ;  getfiles $blk $P_DTA_DIR $F_DTA_DIR 
     fi

     # Chlorophyl file
     tmp=$(LookInNamelist ln_traqsr namelist) ; tmp=$(normalize $tmp )
     if [ $tmp = T ] ; then   # use light penetration
         filter=''
         tmp=$(LookInNamelist ln_qsr_rgb namelist) ; tmp=$(normalize $tmp )   # use RGB parametrization
         if [ $tmp = T ] ; then 
           tmp=$(LookInNamelist nn_chldta namelist)    # use data on file ?
             if [ $tmp = 1 ] ; then 
              blk=namtra_qsr ;  getfiles $blk $P_DTA_DIR $F_DTA_DIR 
                                getweight $blk $P_WEI_DIR $F_WEI_DIR
             fi
         fi
     fi

     # Sea Surface Restoring files
     tmp=$(LookInNamelist ln_ssr namelist) ; tmp=$(normalize $tmp )   # use sea surface restoring
     filter=''
     if [ $tmp = T ] ; then   # use sea surface restoring
        tmp=$(LookInNamelist nn_sstr namelist )       # use SST damping ?
        if [ $tmp = 0 ] ; then filter="$filter | grep -v sn_sst " ; fi
        tmp=$(LookInNamelist nn_sssr namelist )       # use SSS damping ?
        if [ $tmp = 0 ] ; then filter="$filter | grep -v sn_sss " ; fi
        blk=namsbc_ssr ;  getfiles $blk  $P_DTA_DIR $F_DTA_DIR
                          getweight $blk $P_WEI_DIR $F_WEI_DIR
        filter=''
        tmp=$(LookInNamelist ln_sssr_msk namelist namsbc_ssr_drk ) ; tmp=$(normalize $tmp )   # use distance to coast file
        if [ $tmp = T ] ; then  
          blk=namsbc_ssr_drk ;  getfiles $blk  $P_DTA_DIR $F_DTA_DIR
                                getweight $blk $P_WEI_DIR $F_WEI_DIR
        fi

     fi
          }
# ---
# get tidal mixing files (osolete now)
gettmx()  {
        filter=''
        tmp=$(LookInNamelist ln_tmx_itf namelist ) ; tmp=$(normalize $tmp )
        if [ $tmp = F ] ; then filter="$filter | grep -v sn_mskitf " ; fi
        blk=namzdf_tmx ; getfiles $blk  $P_DTA_DIR $F_DTA_DIR
          }
# ---
#  mixing_power_bot mixing_power_pyc mixing_power_cri decay_scale_bot decay_scale_cri
getzdfiwm() { 
        filter=''
        blk=namzdf_iwm_drk ; getfiles $blk  $P_DTA_DIR $F_DTA_DIR
            }
# ---

# get BDY files 
getbdy()  {
        nb_bdy=$(LookInNamelist nb_bdy namelist)
        # (0) bdy_coordinates
        ln_coords_file=( $(LookInNamelist ln_coords_file namelist | sed -e 's/,/ /g' ) )
        cn_coords_file=( $(LookInNamelist cn_coords_file namelist | sed -e 's/,/ /g' ) )
        for i in $( seq 0 $(( nb_bdy -1 )) ) ; do
            tmp=${ln_coords_file[$i]}  ;  tmp=$(normalize $tmp )         
            if [ $tmp = T ] ; then 
               rapatrie ${cn_coords_file[$i]}  $P_BDY_DIR $F_BDY_DIR ${cn_coords_file[$i]}
            fi
        done
        # (0-1) bdy_mask
        ln_mask_file=$(LookInNamelist ln_mask_file namelist ) 
        cn_mask_file=$(LookInNamelist cn_mask_file namelist )
        tmp=$(normalize $ln_mask_file)
        if [ $tmp = T ] ; then
           rapatrie $cn_mask_file  $P_BDY_DIR $F_BDY_DIR $cn_mask_file
        fi
        # (1) data set for dyn2d:
        # get a list of comma separated integer (nb_bdy long)
        nn_dyn2d_dta_lst=($( LookInNamelist nn_dyn2d_dta  namelist | sed -e 's/,/ /g') ) 
        # get a list of comma separated integer (nb_bdy long)
        cn_dyn2d_lst=($( LookInNamelist cn_dyn2d namelist | sed -e 's/,/ /g') ) 
        for bdyset in $( seq 0 $(( nb_bdy-1)) ); do
          if [ ${cn_dyn2d_lst[$bdyset]} != 'none' ] ; then
             nn_dyn2d_dta=${nn_dyn2d_dta_lst[$bdyset]}
             if [ $nn_dyn2d_dta != 0 ] ; then
               echo "  ***  get bdy files for dyn2d, bdy set $bdyset"
               filter='| grep -v bn_u3d | grep -v bn_v3d '                               # skip dyn3d_dta
               filter="$filter | grep -v bn_tem | grep -v bn_sal"                        # skip tra_dta
               filter="$filter | grep -v bn_a_i | grep -v bn_h_i | grep -v bn_h_s"  # skip lim2_dta
#               filter="$filter | grep -v bn_a_i | grep -v bn_ht_i | grep -v bn_ht_s"    # skip lim3_dta
               ln_full_vel=$(LookInNamelist ln_full_vel namelist) ;  ln_full_vel=$(normalize $ln_full_vel)
               if [ $ln_full_vel = T ] ; then
                 filter="$filter | grep -v bn_u2d | grep -v bn_v2d "  # skip barotropic stuff
               fi
               blk=nambdy_dta ; getfiles $blk  $P_BDY_DIR $F_BDY_DIR BDY${bdyset}${nb_bdy}
#              if [ $nn_dyn2d_dta = 3 ] ; then
#
#                ln_bdytide_2ddta=$(LookInNamelist ln_bdytide_2ddta namelist) ; tmp=$(normalize $ln_bdytide_2ddta )
#                if [ $tmp = T ] 
#                    il faut ramener les fichiers <fileroot>_grid[TUV].nc
#                else
#                    file set for each boundary
#                fi
#              fi
             fi
          fi
        done

        # (2) data set for dyn3d:
        # get a list of comma separated integer (nb_bdy long)
        nn_dyn3d_dta_lst=($( LookInNamelist nn_dyn3d_dta  namelist | sed -e 's/,/ /g') )
        # get a list of comma separated integer (nb_bdy long)
        cn_dyn3d_lst=($( LookInNamelist cn_dyn3d namelist | sed -e 's/,/ /g') )
        for bdyset in $( seq 0 $(( nb_bdy-1)) ); do
          if [ ${cn_dyn3d_lst[$bdyset]} != 'none' ] ; then
             nn_dyn3d_dta=${nn_dyn3d_dta_lst[$bdyset]}
             if [ $nn_dyn3d_dta != 0 ] ; then
               echo "  ***  get bdy files for dyn3d, bdy set $bdyset"
               filter='| grep -v bn_ssh | grep -v bn_u2d | grep -v bn_v2d '              # skip dyn2d_dta
               filter="$filter | grep -v bn_tem | grep -v bn_sal"                        # skip tra_dta
               filter="$filter | grep -v bn_a_i | grep -v bn_h_i | grep -v bn_h_s"  # skip lim2_dta
#               filter="$filter | grep -v bn_a_i | grep -v bn_ht_i | grep -v bn_ht_s"    # skip lim3_dta
               blk=nambdy_dta ; getfiles $blk  $P_BDY_DIR $F_BDY_DIR BDY${bdyset}${nb_bdy}
             fi
          fi
        done

        # (3) data set for tracers:
        # get a list of comma separated integer (nb_bdy long)
        nn_tra_dta_lst=($( LookInNamelist nn_tra_dta  namelist | sed -e 's/,/ /g') )
        # get a list of comma separated integer (nb_bdy long)
        cn_tra_lst=($( LookInNamelist cn_tra namelist | sed -e 's/,/ /g') )
        for bdyset in $( seq 0 $(( nb_bdy-1)) ); do
          if [ ${cn_tra_lst[$bdyset]} != 'none' ] ; then
             nn_tra_dta=${nn_tra_dta_lst[$bdyset]}
             if [ $nn_tra_dta != 0 ] ; then
               echo "  ***  get bdy files for tracer, bdy set $bdyset"
               filter='| grep -v bn_ssh | grep -v bn_u2d | grep -v bn_v2d '              # skip dyn2d_dta
               filter="$filter | grep -v bn_u3d | grep -v bn_v3d "                       # skip dyn3d_dta
               filter="$filter | grep -v bn_a_i | grep -v bn_h_i | grep -v bn_h_s"       # skip si3_dta
               blk=nambdy_dta ; getfiles $blk  $P_BDY_DIR $F_BDY_DIR BDY${bdyset}${nb_bdy}
             fi
          fi
        done

        # (4) data set for si3
        # get a list of comma separated integer (nb_bdy long)
        nn_ice_dta_lst=($( LookInNamelist nn_ice_dta  namelist | sed -e 's/,/ /g') )
        # get a list of comma separated integer (nb_bdy long)
        cn_ice_lst=($( LookInNamelist cn_ice namelist | sed -e 's/,/ /g') )
        for bdyset in $( seq 0 $(( nb_bdy-1)) ); do
          if [ ${cn_ice_lst[$bdyset]} != 'none' ] ; then
             nn_ice_dta=${nn_ice_dta_lst[$bdyset]}
             if [ $nn_ice_dta != 0 ] ; then
               echo "  ***  get bdy files for ice model, bdy set $bdyset"
               filter='| grep -v bn_ssh | grep -v bn_u2d | grep -v bn_v2d '              # skip dyn2d_dta
               filter="$filter | grep -v bn_u3d | grep -v bn_v3d"                        # skip dyn3d_dta
               filter="$filter | grep -v bn_tem | grep -v bn_sal"                        # skip tra_dta
               blk=nambdy_dta ; getfiles $blk  $P_BDY_DIR $F_BDY_DIR BDY${bdyset}${nb_bdy}
             fi
          fi
        done

          }
# ---
# get geothermal heating files
getgeo()  {
        filter=''
        blk=nambbc    ; getfiles $blk  $P_DTA_DIR $F_DTA_DIR
                        getweight $blk $P_WEI_DIR $F_WEI_DIR
          }

# ---
# get calving file
getcalving()  {
        filter=''
        nn_test_icebergs=$(LookInNamelist nn_test_icebergs  namelist)
        if [ $nn_test_icebergs = -1 ] ; then
           blk=namberg ; getfiles $blk  $P_DTA_DIR $F_DTA_DIR
        fi
          }
# ---
# get isf files
getisf () {
       filter=''
       nn_isf=$(LookInNamelist nn_isf  namelist)
       blk=namsbc_isf
       # need to get files only for nn_isf = 2 3 or 4
       if [ $nn_isf = 2 ] ; then  # 
         filter='| grep -v sn_fwfisf | grep -v sn_rnfisf' 
         getfiles $blk $P_DTA_DIR $F_DTA_DIR
       fi
       if [ $nn_isf = 3 ] ; then  # 
         filter='| grep -v sn_fwfisf | grep -v sn_Leff_isf ' 
         getfiles $blk $P_DTA_DIR $F_DTA_DIR
       fi
       if [ $nn_isf = 4 ] ; then  # 
         filter='| grep -v sn_rnfisf | grep -v sn_Leff_isf | grep -v sn_depmax_isf | grep -v sn_depmin_isf' 
         getfiles $blk $P_DTA_DIR $F_DTA_DIR
       fi
          }
# ---
# get diaptr subbasins mask obsolete in Nemo4 so far ... (harcoded names)
getdiaptr()  {
        filter=''
        blk=namptr_drk  ; getfiles $blk $P_DTA_DIR $F_DTA_DIR
             }
# ---
# get obs data file
getobs () {
  if [ $no != 1 ] ; then
    ndastpdeb=$( tail -2 $CONFIG_CASE.db | head -1 | awk '{print $4}' )
  else
    ndastpdeb=$(LookInNamelist nn_date0)
  fi

  # hard coded file names
  root_enact=EN3_v2a_Profiles_
  root_sla=fdbk_j2_
  slaRefLevel='slaReferenceLevel.nc'

  rdt=$(LookInNamelist rn_rdt)
  rdt=$(echo 1 | awk "{ rdt=int($rdt); print rdt}" )

  ndays=$( echo 1 | awk "{ a=int( ($nitend - $nit000 +1)*$rdt /86400.) ; print a }" )
  ndastpfin=$( ./datfinyyyy $ndastpdeb $ndays )

  yyyy1=${ndastpdeb:0:4}
  mm1=${ndastpdeb:4:2}

  yyyy2=${ndastpfin:0:4}
  mm2=${ndastpfin:4:2}

 if [ $ENACT = 1 ] ; then
   # ENACT
   flist=''
   for y in $(seq $yyyy1 $yyyy2) ; do
     for m in $(seq -f '%02g' 1 12 ) ; do
       f=$root_enact$y$m.nc
       if [ $y = $yyyy1 ] ; then
         if [ $m -ge $mm1 ] ; then 
           if [ -f $P_ENA_DIR/$f ] ; then
             flist="$flist '$f' "
             rapatrie $f  $P_ENA_DIR $F_ENA_DIR $f
           fi
         fi
       elif [ $y = $yyyy2 ] ; then
         if [ $m -le $mm2 ] ; then 
           if [ -f $P_ENA_DIR/$f ] ; then
             flist="$flist '$f' "
             rapatrie $f  $P_ENA_DIR $F_ENA_DIR $f
           fi
         fi
       else
         if [ -f $P_ENA_DIR/$f ] ; then
           flist="$flist '$f' "
           rapatrie $f  $P_ENA_DIR $F_ENA_DIR $f
         fi
       fi
     done
   done

   cat namelist | sed -e "s/ENACTFILES_LIST/$flist/" > znamelist1
   mv znamelist1 namelist
 fi

 if [ $SLA = 1 ] ; then
   rapatrie $slaRefLevel $P_SLA_DIR $F_SLA_DIR $slaRefLevel
   flist=''
   for y in $(seq $yyyy1 $yyyy2) ; do
      f=$root_sla$y.nc
      if [ -f  $P_SLA_DIR/$f ] ; then
        flist="$flist '$f' "
        rapatrie $f  $P_SLA_DIR $F_SLA_DIR $f
      fi
   done

   cat namelist | sed -e "s/SLAFBFILES_LIST/$flist/" > znamelist1
   mv znamelist1 namelist
 fi
        }
#---

# getweight files associated with block blk in the namelist
getweight() {
         ZPDIR=$2
         ZFDIR=$3
         # if a 4th argument is passed, it is a namelist name. Default to 'namelist'
         if [ $4 ] ; then namelist=$4 ; else namelist=namelist ; fi
          cmd="getblock $1 $namelist |  grep sn_  $filter   |  awk -F, '{ print \$7}'  "
          lstw=$(eval $cmd | tr -d "'" | tr -d "," | sort -u )

         echo required weight files :
         echo =======================
         for f in $lstw ; do
            echo $f
#           copyfor $P_WEI_DIR/$f $f
#           copyfor $ZPDIR/$f $f
            rapatrie  $f $ZPDIR $ZFDIR  $f
         done
            }
# ---

# get data files associated with the block namelist ( eventually filtered )
getfiles()  {
         # This function is now accepting bdy request
         # the key "sn_" is thus now a variable that can take either sn_ or bn_ 9for BDY)
         zstr="sn_"
         zbdyset=""
         bdyflag=0
         # if a 4th argument is passed it can be either BDY{bdyset}{nbdy} or a namelist filename
         if [ $4 ] ; then 
             if [ $( echo $4 | awk '{ print index( $1,"BDY")}')  = 1 ] ; then
                # case of a BDY request. There will be repeated namelist block for nambdy_dta. Need to keep track
                # we do know localy the current bdy set making the call and the total number of bdy by decrypting
                # the BDYxx keyword
                zstr="bn_"
                ztmp=${4#BDY}
                zbdyset=${ztmp:0:1}   # assume key like BDY03 BDY13 BDY23 
                znbdy=${ztmp:1:1}      # may be improve with a key like BDY-00-03 ...( for bdy>9 ! )
                namelist=namelist 
                bdyflag=1
             else
                namelist=$4 
             fi
         else 
             namelist=namelist 
         fi
         ZPDIR=$2
         ZFDIR=$3

         cmd="getblock $1 $namelist |  awk '{if ( index(\$1,\"$zstr\") != 0 ) print \$0}'   $filter |  awk '{ print \$3 }'"
         lst=$(eval $cmd | tr -d "'" | tr -d ",")
         # check for eventual weight in the namelist block ( WARNING : problem if mixed ... )
         lstw=$(getblock $1 $namelist | grep $zstr  |  awk -F, '{ print $7 }' | tr -d "'" | tr -d " " | sort -u )

         if [ $zbdyset ] ; then  # dealing with bdy data
            # bdy data are not using weight files ( at leat for now)
            lstw=''
            lst_arr=($lst)
            lst=''  # reset lst to ''
            nfld=$(( ${#lst_arr[@]} / znbdy ))   # number of fields per bdy set
            for ifld in $(seq 1 $nfld) ; do
                ii=$(( (zbdyset  )*nfld +ifld -1 ))
                lst="$lst  ${lst_arr[$ii ]}"     # rebuild lst for bdy set zbdyset
                echo  $zbdyset $ifld $ii  ${lst_arr[$ii ]}
            done
         fi

         for f in $lst ; do
           if [ $f != none ] ; then   # non used bn_xxx are set to none
           ln_clim=$(getblock $1 $namelist | awk '{if ( index($1,zstr) != 0 ) print $0}' zstr=$zstr | grep $f  |  awk -F, '{ print $5 }' | sort -u )
           ln_clim=$(normalize $ln_clim )
           # look for file_type (ie yearly monthly weekly ... )
           file_type=$(getblock $1 $namelist | awk '{if ( index($1,zstr) != 0 ) print $0}' zstr=$zstr | grep $f  |  awk -F, '{ print $6 }' | tr -d "'" | tr -d " " | sort -u )

           case $file_type in 
           ( yearly )  # get files for 3 years 
           if [ $ln_clim = F ] ; then
             for y in $(seq -f "%04g" $yearm1 $yearp1) ; do
               echo ${f}_y${y}.nc
#              copyfor  $P_FOR_DIR/${f}_y${y}.nc ${f}_y${y}.nc
               if [ ${#lstw} != 0 -o $bdyflag = 1 ] ; then   # there are weight defined, so we do not try to copy nested files (??? AGRIF ??? )
                 copyfor  $ZPDIR/${f}_y${y}.nc ${f}_y${y}.nc
               else
                 rapatrie  ${f}_y${y}.nc $ZPDIR $ZFDIR  ${f}_y${y}.nc
               fi
             done
           else
             echo ${f}.nc
#            copyfor  $P_FOR_DIR/${f}.nc ${f}.nc
             if [ ${f} != 'NOT' ] ; then
             if [ ${#lstw} != 0   -o $bdyflag = 1 ] ; then   # there are weight defined, so we do not try to copy nested files ( ??? AGRIF ??? )
                copyfor  $ZPDIR/${f}.nc ${f}.nc
             else
                rapatrie  ${f}.nc $ZPDIR $ZFDIR  ${f}.nc
             fi
             fi
           fi ;;
           ( monthly ) # get file m12 for yearm1, then m01 -m12 for current year , m01 for yearp1
           if [ $ln_clim = F ] ; then
             for y in $(seq -f "%04g" $yearm1 $yearp1) ; do
               if [ $y = $yearm1 ] ; then
                 echo ${f}_y${y}m12.nc
                 if [ ${#lstw} != 0  -o $bdyflag = 1 ] ; then   # there are weight defined, so we do not try to copy nested files (??? AGRIF ??? )
                   copyfor  $ZPDIR/${f}_y${y}m12.nc ${f}_y${y}m12.nc
                 else
                   rapatrie  ${f}_y${y}m12.nc $ZPDIR $ZFDIR  ${f}_y${y}m12.nc
                 fi
               elif [ $y = $yearp1 ] ; then
                 echo ${f}_y${y}m01.nc
                 if [ ${#lstw} != 0 ] ; then   # there are weight defined, so we do not try to copy nested files (??? AGRIF ??? )
                   copyfor  $ZPDIR/${f}_y${y}m01.nc ${f}_y${y}m01.nc
                 else
                   rapatrie  ${f}_y${y}m01.nc $ZPDIR $ZFDIR  ${f}_y${y}m01.nc
                 fi
               else
                 for  mm in $(seq -w 01 12 ) ; do
                    echo ${f}_y${y}m${mm}.nc
                    if [ ${#lstw} != 0  -o $bdyflag = 1 ] ; then   # there are weight defined, so we do not try to copy nested files (??? AGRIF ??? )
                      copyfor  $ZPDIR/${f}_y${y}m${mm}.nc ${f}_y${y}m${mm}.nc
                    else
                      rapatrie  ${f}_y${y}m${mm}.nc $ZPDIR $ZFDIR  ${f}_y${y}m${mm}.nc
                    fi
                 done
                   
               fi
             done
           else
             for  mm in $(seq -w 01 12 ) ; do
                echo ${f}_m${mm}.nc
                if [ ${#lstw} != 0   -o $bdyflag = 1 ] ; then   # there are weight defined, so we do not try to copy nested files ( ??? AGRIF ??? )
                  copyfor  $ZPDIR/${f}_m${mm}.nc ${f}_m${mm}.nc
                else
                  rapatrie  ${f}_m${mm}.nc $ZPDIR $ZFDIR  ${f}_m${mm}.nc
                fi
             done
           fi ;;
           esac
           fi
         done   
            }
# ---

# Function returning the extension to be use for a given member of 
#     an ensemble run or nothing [ eg : .002 ]
getmember_extension()  {
      if [ $1 = -1 ] ; then
        echo ""  
      else
        if [ $# = 2 ] ; then
           printf "%03d" $1
        else
           printf ".%03d" $1
        fi
      fi
                       }
# ---

# Give the number of file containing $1 in its name
number_of_file() { ls -1 *$1* 2> /dev/null | wc -l  ; }
# ---

# function for renaming restart files according to current run
#       example : renamerst restart_ice_in restart_ice 
renamerst() {
            echo '   *** renaming restartfile to ' ${1}_xxxx.${filext}.$ext
 # look for nitend in namelist. Only the last restarts are to be renamed
 # in case of multiple restarts during the year
 znitend=$(LookInNamelist nn_itend  namelist)
 znitend=$( printf "%08d" $znitend )
 filext='nc'

 # Loop on cores for current member [ renamerst is called in the member loop, mmm is known ]
 set -x   # avoid very long list of statements in the log file
 zrstdir='./'
 if [ $RST_DIRS = 1 ] ; then zrstdir=$DDIR/${CN_DIRRST}.$ext/$nnn/ ; fi

#          clname = TRIM(cexper)//"_"//TRIM(ADJUSTL(clkt))//"_"//TRIM(cn_ocerst_out)
#          rstout = <CONFIG_CASE><.MBR>_<nitend>_<CN_RST_OUT>_<CORE>.<FILEXT>

 cd $zrstdir
 if [ -f ${CONFIG_CASE}${mmm}_${znitend}_${2}_${zcore0}.${filext} ] ; then
 for rest in ${CONFIG_CASE}${mmm}_${znitend}_${2}_[[:digit:]]*[[:digit:]].${filext} ; do
     if [ $RST_READY = 1 ] ; then
       CORE=${rest%.*} ; CORE=${CORE##*_}
       rest_in=${1}-${ext}${mmm}_${CORE}.${filext}
       mv $rest $rest_in
     else
       rest_in=${1}${mmm}_$( tmp=${rest##*_} ; echo ${tmp%.${filext}} ).${filext}
       mv $rest $rest_in.$ext
     fi
 done
 fi

 # Agrif case if necessary
 set -x
 if [ $AGRIF = 1 ] ; then
   for idx in  ${agrif_pref[@]} ; do
       znitend=$(LookInNamelist nn_itend ${idx}_namelist  )
       znitend=$( printf "%08d" $znitend )
       set +x   # avoid very long list of statements in the log file
       if [ -f ${idx}_${CONFIG_CASE}${mmm}_${znitend}_${2}_${zcore0}.${filext} ] ; then
       for rest in ${idx}_${CONFIG_CASE}${mmm}_${znitend}_${2}_[[:digit:]]*[[:digit:]].${filext} ; do
          if [ $RST_READY = 1 ] ; then
             CORE=${rest%.*} ; CORE=${CORE##*_}
             rest_in=${idx}_${1}-${ext}${mmm}_${CORE}.${filext}
             mv $rest $rest_in
          else
             rest_in=${idx}_${1}${mmm}_$( tmp=${rest##*_} ; echo ${tmp%.${filext}} ).${filext}
             mv $rest $rest_in.$ext
          fi
       done
       fi
       set -x
   done
 fi
 cd -
            }
# ---

# function that rename  and copy to apropriate P_S_DIR all text like files after a run 
#( ocean.ouput, namelist etc... ) taking care of members and AGRIF.
# It takes the extension as argument to be more generic (error case for instance)
rename_txt_files() {
     exten=$1
     EXTRA_COPY=
     if [ $# = 2 ] ; then  EXTRA_COPY="yes" ; fi
     # files common to all members
                              mv namelist      namelist_oce.$exten
      if [ $ICE = 1  ] ; then mv namelist_ice  namelist_ice.$exten     ; fi
      if [ $TOP = 1  ] ; then mv namelist_top  namelist_top.$exten     ; fi
      if [ $CFC = 1  ] ; then mv namelist_cfc  namelist_cfc.$exten     ; fi

                              cp namelist_oce.$exten $P_S_DIR/ANNEX
      if [ $ICE = 1  ] ; then cp namelist_ice.$exten $P_S_DIR/ANNEX    ; fi
      if [ $TOP = 1  ] ; then cp namelist_top.$exten $P_S_DIR/ANNEX    ; fi
      if [ $CFC = 1  ] ; then cp namelist_cfc.$exten $P_S_DIR/ANNEX    ; fi

      # member dependent files
      for member in  $(seq $ENSEMBLE_START $ENSEMBLE_END) ; do
         mmm=$(getmember_extension $member)
                                    mv ocean.output$mmm  ocean.output$mmm.$exten
         if [ $RUNSTAT = 1 ] ; then mv run.stat$mmm run.stat$mmm.$exten ; fi
         if [ $RUNSTAT = 1 ] ; then mv run.stat$mmm.nc run.stat$mmm.$exten.nc ; fi
         if [ $ntiming = 1 ] ; then mv timing.output$mmm timing.output$mmm.$exten ; fi

         ZP_S_DIR=$( echo $P_S_DIR | sed -e "s;$CONFIG_CASE;${CONFIG_CASE}$mmm;")
         mkdir -p $ZP_S_DIR/ANNEX
                                    cp ocean.output$mmm.$exten  $ZP_S_DIR/ANNEX
         if [ $RUNSTAT = 1 ] ; then cp run.stat$mmm.$exten  $ZP_S_DIR/ANNEX ; fi
         if [ $RUNSTAT = 1 ] ; then cp run.stat$mmm.$exten.nc  $ZP_S_DIR/ANNEX ; fi
         if [ $ntiming = 1 ] ; then cp timing.output$mmm.$exten $ZP_S_DIR/ANNEX ; fi
         if [ $EXTRA_COPY  ] ; then
                                         copy ocean.output$mmm.$exten  $2
            if [ $RUNSTAT = 1 ] ; then   copy run.stat$mmm.$exten  $2 ; fi
         fi

      done

      # treat case of AGRIF runs
      if [ $AGRIF = 1 ] ; then
         for idx in ${agrif_pref[@]} ; do
            # files common to all members
                                    mv ${idx}_namelist      ${idx}_namelist_oce.$exten
            if [ $ICE = 1  ] ; then mv ${idx}_namelist_ice  ${idx}_namelist_ice.$exten     ; fi
            if [ $TOP = 1  ] ; then mv ${idx}_namelist_top  ${idx}_namelist_top.$exten     ; fi
            if [ $CFC = 1  ] ; then mv ${idx}_namelist_cfc  ${idx}_namelist_cfc.$exten     ; fi

                                    cp ${idx}_namelist_oce.$exten $P_S_DIR/ANNEX
            if [ $ICE = 1  ] ; then cp ${idx}_namelist_ice.$exten $P_S_DIR/ANNEX    ; fi
            if [ $TOP = 1  ] ; then cp ${idx}_namelist_top.$exten $P_S_DIR/ANNEX    ; fi
            if [ $CFC = 1  ] ; then cp ${idx}_namelist_cfc.$exten $P_S_DIR/ANNEX    ; fi

            # member dependent files
            for member in  $(seq $ENSEMBLE_START $ENSEMBLE_END) ; do
               mmm=$(getmember_extension $member)
                                          mv ${idx}_ocean.output$mmm  ${idx}_ocean.output$mmm.$exten
                                          mv ${idx}_run.stat$mmm   ${idx}_run.stat$mmm.$exten
               if [ $ntiming = 1 ] ; then mv ${idx}_timing.output$mmm ${idx}_timing.output$mmm.$exten ; fi

               ZP_S_DIR=$( echo $P_S_DIR | sed -e "$;$CONFIG_CASE;${CONFIG_CASE}$mmm;")
               mkdir -p $ZP_S_DIR/ANNEX
                                          cp ${idx}_ocean.output$mmm.$exten  $ZP_S_DIR/ANNEX
                                          cp ${idx}_run.stat$mmm.$exten  $ZP_S_DIR/ANNEX
               if [ $ntiming = 1 ] ; then cp ${idx}_timing.output$mmm.$exten $ZP_S_DIR/ANNEX ; fi
               if [ $EXTRA_COPY  ] ; then
                                          copy ${idx}_ocean.output$mmm.$exten  $2
                                          copy ${idx}_run.stat$mmm.$exten  $2
               fi
            done
         done
      fi
                   }

# function for making restart tar file of about 1 Gb. Take 3 arguments: root_name_in root_name_out option
#        example : mktarrst restart restart _oce.v2 --> restart_oce.v2.1.tar.ext
mktarrst() {  echo "   *** making tar restart file for ${1}$3 "

 # ... build tar files of size < MAXTARSIZ (in bytes)
 MAXTARSIZ=40000000000
 set +x   # avoid very long list of file
 filext='nc'
 cd $zrstdir
 if [ $RST_READY = 1 ] ; then
    lscmd="/bin/ls -l  ${1%.*}-${ext}${mmm}_[[:digit:]]*[[:digit:]].${filext}"
 else
    lscmd="/bin/ls -l  ${1}_[[:digit:]]*[[:digit:]].${filext}.$ext"
 fi

 $lscmd | awk ' BEGIN {s=0 ; arch=0 ; list= ""}    \
       { if ( s < maxtar ) { s=s+$5 ; list=list " " $NF } \
         else { arch=arch+1 ;
         cmd="tar cf " dir name option "." ext ".tar." arch " " list ; \
          system (cmd ); s= $5; list=$NF } \
       } \
       END { arch=arch+1 ; cmd="tar cf  " dir name option "." ext ".tar." arch " " list ; \
             system (cmd ) } ' maxtar=$MAXTARSIZ ext=$ext name=$2 option=$3 dir=${P_R_DIR}/

 echo " end of restart tarfiles ." ; set -x 
 if [ $AGRIF = 1 ] ; then
    for idx in $nst_lst ; do
       set +x   # avoid very long list of file
       if [ $RST_READY = 1 ] ; then
         lscmd="/bin/ls -l  ${idx}_${1%.*}-${ext}${mmm}_[[:digit:]]*[[:digit:]].${filext}"
       else
         lscmd="/bin/ls -l  ${idx}_${1}_[[:digit:]]*[[:digit:]].${filext}.$ext"
       fi

       $lscmd  | awk ' BEGIN {s=0 ; arch=0 ; list= ""}    \
              { if ( s < maxtar ) { s=s+$5 ; list=list " " $NF } \
                else { arch=arch+1 ;
                      cmd="tar cf " dir name option "." ext ".tar." arch " " list ; \
                      system (cmd ); s= $5; list=$NF } \
              } \
               END { arch=arch+1 ; cmd="tar cf  " dir name option "." ext ".tar." arch " " list ; \
             system (cmd ) } ' maxtar=$MAXTARSIZ ext=$ext name=${idx}_$2 option=$3 dir=${P_R_DIR}/

         echo " end of restart tarfiles ." ; set -x 
    done
 fi
cd -
           }
# ---

# function for building script to save restart: 
#   This function buildt a submit script (with header) and as many secondary scripts as members in an ensemble run.
#   The submit script, then launch the execution in parallel of the secondary scripts (as many cores as members).
mksavrst() {

mk_batch_hdr --name ${1%%.*} --par --wallclock 0:20:00 --cores $ENSEMBLE_SIZE --cluster hpt --account $ACCOUNT --adapp --queue $QUEUE > $1

cat << eof >> $1    # Submit script name given as argument
 set -x
 cd $TMPDIR
 . ./includefile.sh
 . $RUNTOOLS/lib/function_4_all.sh
 . $RUNTOOLS/lib/function_4.sh
 srcbash # just in case
 # set local variables as in nemo4.sh (calling script)
 ext=$ext
 AGRIF=$AGRIF
 XIOS=$XIOS
 RST_DIRS=$RST_DIRS
 nst_lst="${agrif_pref[@]} "

 mpmd_arg=""
 for member in \$(seq $ENSEMBLE_START $ENSEMBLE_END) ; do
   mmm=\$(getmember_extension \$member)
   nnn=\$(getmember_extension \$member nodot )
   zrstdir='./'
   if [ \$RST_DIRS = 1 ] ; then zrstdir=$DDIR/${CN_DIRRST}.\$ext/\$nnn ; fi
   # create secondary scripts to be submitted in // trhough the members
   # $ to be maintained in the final script are replaces by @, then automatic edition
   # replace the @ by $ [ this is necessary because we are at the second level of script
   # creation !
   cat << eof1 > ztmp
#!/bin/bash
   set -x
   . ./includefile.sh
   . $RUNTOOLS/lib/function_4_all.sh
   . $RUNTOOLS/lib/function_4.sh

   # set local variables as in nemo3.4 (calling script)
   ext=$ext
   AGRIF=$AGRIF
   XIOS=$XIOS
   nst_lst="${agrif_pref[@]} "
   mmm=\$mmm
   zrstdir=\$zrstdir
   cd $DDIR
   tar cf $F_R_DIR/${CONFIG_CASE}\${mmm}-RST.$ext.tar ${CONFIG_CASE}-RST.$ext/\$mmm
eof1
   cat ztmp | sed -e 's/@/\$/g' > ./$1\${mmm}.sh    # change @ into \$ and create script for current member
   chmod 755 ./$1\${mmm}.sh                         # made it executable

   mpmd_arg="\$mpmd_arg 1 ./$1\${mmm}.sh"           # prepare the command line for runcode function
 done
  
  pwd
  if [ ! \$mmm ] ; then
     ./$1\${mmm}.sh                                 # not an ensemble run : serial process
  else
     runcode_mpmd  \$mpmd_arg                        # launch the scripts in parallele (mpmd mode)
  fi
eof

  # copy the script to P_CTL_DIR from where it will be launched by submit
  copy $1 $P_CTL_DIR                                 # for reference
              }
# --- 


# function for building script to save restart: 
#   This function buildt a submit script (with header) and as many secondary scripts as members in an ensemble run.
#   The submit script, then launch the execution in parallel of the secondary scripts (as many cores as members).
mksavrst2() {

mk_batch_hdr --name ${1%%.*} --par --wallclock 2:00:00 --cores $ENSEMBLE_SIZE --cluster hpt --account $ACCOUNT --adapp --queue $QUEUE > $1

cat << eof >> $1    # Submit script name given as argument
 set -x
 cd $TMPDIR
 . ./includefile.sh
 . $RUNTOOLS/lib/function_4_all.sh
 . $RUNTOOLS/lib/function_4.sh
 srcbash # just in case
 # set local variables as in nemo3.4 (calling script)
 ext=$ext
 AGRIF=$AGRIF
 XIOS=$XIOS
 RST_DIRS=$RST_DIRS
 nst_lst="${agrif_pref[@]} "

 mpmd_arg=""
 for member in \$(seq $ENSEMBLE_START $ENSEMBLE_END) ; do
   mmm=\$(getmember_extension \$member)
   nnn=\$(getmember_extension \$member nodot )
   zrstdir='./'
   if [ \$RST_DIRS = 1 ] ; then zrstdir=$DDIR/${CN_DIRRST}.\$ext/\$nnn ; fi
   # create secondary scripts to be submitted in // trhough the members
   # $ to be maintained in the final script are replaces by @, then automatic edition
   # replace the @ by $ [ this is necessary because we are at the second level od script
   # creation !
   cat << eof1 > ztmp
#!/bin/bash
   set -x
   . ./includefile.sh
   . $RUNTOOLS/lib/function_4_all.sh
   . $RUNTOOLS/lib/function_4.sh

   # set local variables as in nemo3.4 (calling script)
   ext=$ext
   AGRIF=$AGRIF
   XIOS=$XIOS
   nst_lst="${agrif_pref[@]} "
   mmm=\$mmm
   zrstdir=\$zrstdir

   # O C E A N
   # *********
   OCE_RST_IN=@(LookInNamelist cn_ocerst_in namelist_oce.$ext)\$mmm
   OCE_RST_OUT=@(LookInNamelist cn_ocerst_out namelist_oce.$ext)\$mmm

   mktarrst @OCE_RST_IN @OCE_RST_OUT _oce_v2

   # send them on DATA space
   cd \$P_R_DIR
   for f in @{OCE_RST_OUT}_oce_v2.\$ext.tar.*  ; do
      expatrie_res @f \$F_R_DIR @f \$P_R_DIR
   done
   if [ $AGRIF = 1 ] ; then
     for idx in \$nst_lst ; do 
       for f in @{idx}_@{OCE_RST_OUT}_oce_v2.\$ext.tar.*  ; do
          expatrie_res @f \$F_R_DIR @f \$P_R_DIR
       done
     done
   fi
   cd $TMPDIR

   # I C E
   # *****
   if [ $ICE = 1 ] ; then
     ICE_RST_IN=@(LookInNamelist cn_icerst_in namelist_ice.$ext)\$mmm
     ICE_RST_OUT=@(LookInNamelist cn_icerst_out namelist_ice.$ext)\$mmm

     mktarrst @ICE_RST_IN @ICE_RST_OUT _v2

   # send them on DATA space
    cd \$P_R_DIR
    for f in @{ICE_RST_OUT}_v2.\$ext.tar.*  ; do
       expatrie_res @f \$F_R_DIR @f \$P_R_DIR
    done
    if [ $AGRIF = 1 ] ; then
    for idx in \$nst_lst ; do 
      for f in @{idx}_@{ICE_RST_OUT}_v2.\$ext.tar.*  ; do
         expatrie_res @f \$F_R_DIR @f \$P_R_DIR
      done
    done
    fi
   fi
   cd $TMPDIR

   # P A S S I V E   T R A C E R S
   # *****************************
   if [ $TOP = 1 ] ; then
     TRC_RST_IN=@(LookInNamelist   cn_trcrst_in namelist_top.$ext)\$mmm
     TRC_RST_OUT=@(LookInNamelist cn_trcrst_out namelist_top.$ext)\$mmm

     mktarrst @TRC_RST_IN @TRC_RST_OUT _v2
  
    # send them on  DATA space
    cd \$P_R_DIR
    for f in @{TRC_RST_OUT}_v2.\$ext.tar.*  ; do
      expatrie_res @f \$F_R_DIR @f \$P_R_DIR
    done
    if [ $AGRIF = 1 ] ; then
      for idx in \$nst_lst ; do 
        for f in @{idx}_@{TRC_RST_OUT}_v2.\$ext.tar.*  ; do
           expatrie_res @f \$F_R_DIR @f \$P_R_DIR
        done
      done
    fi
   fi
   cd $TMPDIR
  
   # T R D  M L D
   # ************
   if [ $TRDMLD = 1 ] ; then
      TRD_RST_IN=@(LookInNamelist cn_trdrst_in namelist_oce.$ext)\$mmm
      TRD_RST_OUT=@(LookInNamelist cn_trdrst_out namelist_oce.$ext)\$mmm
  
      mktarrst @TRD_RST_IN @TRD_RST_OUT _v2
  
      # send them on  DATA space
      cd \$P_R_DIR
      for f in @{TRD_RST_OUT}_v2.\$ext.tar.*  ; do
         expatrie_res @f \$F_R_DIR @f \$P_R_DIR
      done
      if [ $AGRIF = 1 ] ; then
         for idx in \$nst_lst ; do 
            for f in @{idx}_@{TRD_RST_OUT}_v2.\$ext.tar.*  ; do
               expatrie_res @f \$F_R_DIR @f \$P_R_DIR
            done
         done
     fi
     cd $TMPDIR
   fi

   # S T O 
   # *****
   if [ $STO = 1 ] ; then
      STO_RST_IN=@(LookInNamelist cn_storst_in namelist_oce.$ext)\$mmm
      STO_RST_OUT=@(LookInNamelist cn_storst_out namelist_oce.$ext)\$mmm
      mktarrst @STO_RST_IN @STO_RST_OUT _v2
  
      # send them on  DATA space
      cd \$P_R_DIR
      for f in @{STO_RST_OUT}_v2.\$ext.tar.*  ; do
         expatrie_res @f \$F_R_DIR @f \$P_R_DIR
      done
     if [ $AGRIF = 1 ] ; then
        for idx in \$nst_lst ; do 
           for f in @{idx}_@{STO_RST_OUT}_v2.\$ext.tar.*  ; do
              expatrie_res @f \$F_R_DIR @f \$P_R_DIR
           done
        done
     fi
     cd $TMPDIR
   fi

   # I C B 
   # *****
   if [ $ICB = 1 ] ; then
      ICB_RST_IN=restart_icebergs\$mmm
      ICB_RST_OUT=icebergs_restart\$mmm
      mktarrst @ICB_RST_IN @ICB_RST_OUT _v2
  
      # send them on  DATA space
      cd \$P_R_DIR
      for f in @{ICB_RST_OUT}_v2.\$ext.tar.*  ; do
         expatrie_res @f \$F_R_DIR @f \$P_R_DIR
      done
     if [ $AGRIF = 1 ] ; then
        for idx in \$nst_lst ; do 
           for f in @{idx}_@{ICB_RST_OUT}_v2.\$ext.tar.*  ; do
              expatrie_res @f \$F_R_DIR @f \$P_R_DIR
           done
        done
     fi
     cd $TMPDIR
   fi

   touch RST_DONE\${mmm}.\$ext
eof1
   cat ztmp | sed -e 's/@/\$/g' > ./$1\${mmm}.sh    # change @ into \$ and create script for current member
   chmod 755 ./$1\${mmm}.sh                         # made it executable

   mpmd_arg="\$mpmd_arg 1 ./$1\${mmm}.sh"           # prepare the command line for runcode function
 done
  
  pwd
  if [ ! \$mmm ] ; then
     ./$1\${mmm}.sh                                 # not an ensemble run : serial process
  else
     runcode_mpmd  \$mpmd_arg                        # launch the scripts in parallele (mpmd mode)
  fi
eof

  # copy the script to P_CTL_DIR from where it will be launched by submit
  copy $1 $P_CTL_DIR                                 # for reference
              }
# --- 

mkbuild_mesh_mask() { 
mk_batch_hdr --name ${1%%.*} --seq --cores 1 --nodes 1 --wallclock 1:00:00  --account $ACCOUNT --cluster hpt  --queue $QUEUE > $1
cat << eof >> $1

set -x
 . ./includefile.sh
DDIR=\${DDIR:-\$CDIR}
cd \$DDIR/TMPDIR_${CONFIG_CASE}
. $RUNTOOLS/lib/function_4.sh
. $RUNTOOLS/lib/function_4_all.sh
srcbash

# if second argument present represent agrif nest number
AGRIF=$AGRIF # set if code compiled with key_agrif
                            agopt=''
if [ \$AGRIF = 1 ] ; then   agopt='-A'  ;  fi

idx=${2-0}
                          pref=''    
if [ \$idx != 0  ] ; then pref=\${idx}_ ; fi

# copy executable in the TMPDIR_xxx directory
if [ ! -f build_nc_iom  ] ; then rcopy \$P_UTL_DIR/bin/build_nc_iom ./  ; fi
if [ ! -f correct_holes ] ; then rcopy \$P_UTL_DIR/bin/correct_holes ./ ; fi

# BUILDNC_2.0
./build_nc_iom -r \${pref}mask     -p by \$agopt
./build_nc_iom -r \${pref}mesh_hgr -p r4 \$agopt
./build_nc_iom -r \${pref}mesh_zgr -p r4 \$agopt

# now correct holes
if [ \$idx = 0 ] ; then 
  ./correct_holes
else
  ./correct_holes -agrif \$idx
fi

# rename files:
  mv \${pref}mask.nc     \${pref}${CONFIG_CASE}_byte_mask.nc
  mv \${pref}mesh_hgr.nc \${pref}${CONFIG_CASE}_mesh_hgr.nc
  mv \${pref}mesh_zgr.nc \${pref}${CONFIG_CASE}_mesh_zgr.nc

# send then to $F_I_DIR and $P_I_DIR
  expatrie \${pref}${CONFIG_CASE}_mesh_hgr.nc   $F_I_DIR \${pref}${CONFIG}-${CASE}_mesh_hgr.nc
  expatrie \${pref}${CONFIG_CASE}_mesh_zgr.nc   $F_I_DIR \${pref}${CONFIG}-${CASE}_mesh_zgr.nc
  expatrie \${pref}${CONFIG_CASE}_byte_mask.nc  $F_I_DIR \${pref}${CONFIG}-${CASE}_byte_mask.nc

  mv  \${pref}${CONFIG_CASE}_mesh_hgr.nc  $P_I_DIR
  mv  \${pref}${CONFIG_CASE}_mesh_zgr.nc  $P_I_DIR
  mv  \${pref}${CONFIG_CASE}_byte_mask.nc $P_I_DIR

eof

  # copy the script to P_CTL_DIR from where it will be launched by submit
  copy $1 $P_CTL_DIR
                    }
# ---

mkbuild_xios() {
mk_batch_hdr --name ${1%%.*} --seq --wallclock 16:00:00  --account $ACCOUNT --cluster hpt  --adapp --queue $QUEUE > $1
cat << eof >> $1

set -x
CONFIG=$CONFIG
CONFIG_CASE=${CONFIG_CASE}
DDIR=\${DDIR:-\$CDIR}

cd $DDIR/\${CONFIG_CASE}-XIOS.$ext

ln -sf $P_UTL_DIR/bin/rebuild_nemo .
ln -sf $P_UTL_DIR/bin/rebuild_nemo.exe .

. $RUNTOOLS/lib/function_4_all.sh
. $RUNTOOLS/lib/function_4.sh


for freq in 1ts 1h 1d 5d ; do
mkdir \${freq}_OUTPUT
  typlst=\$( ls -1 \${CONFIG_CASE}_\${freq}_*nc | awk -F_ '{print \$3}' | sort -u )
#  for grid in gridT gridU gridV gridW icemod flxT icb; do
  for grid in \$typlst ; do
     # look for the number of server
     # check 0000.nc file always existing
     for f in \$(ls \${CONFIG_CASE}_\${freq}_\${grid}_*-*_0000.nc) ; do
        fbase=\${f%_0000.nc}
        tag=\${fbase##*-}
        ndastp=y\${tag:0:4}m\${tag:4:2}d\${tag:6:2}.\$freq
        echo \$fbase \$tag \$ndastp
        nio=\$( ls -1 \${fbase}_[[:digit:]]*[[:digit:]].nc  | wc -l )
        echo \$nio for  \${CONFIG_CASE}_\${freq} \${grid}
        if [ \$nio != 0 ] ; then
         ./rebuild_nemo \$fbase \$nio
         mv \$fbase.nc \${freq}_OUTPUT/\${CONFIG_CASE}_\${ndastp}_\$grid.nc
         # save file on the archiving system using save_arch_file function (machine dependent).
         # save_arch_file   file arch_directory 
         save_arch_file \${freq}_OUTPUT/\${CONFIG_CASE}_\${ndastp}_\$grid.nc \${CONFIG}/\${CONFIG_CASE}-S/\${freq}/\${tag:0:4}
        fi
    done
  done
done

eof
  # copy the script to P_CTL_DIR from where it will be launched by submit
  copy $1 $P_CTL_DIR
          }
# ---

mkbuild_xios_ens() {
mk_batch_hdr --name ${1%%.*} --cores $ENSEMBLE_SIZE --wallclock 4:00:00  \
             --account $ACCOUNT --cluster hpt  --adapp --queue $QUEUE > $1
cat << eof >> $1

set -x
CONFIG=$CONFIG
CONFIG_CASE=${CONFIG_CASE}
DDIR=\${DDIR:-\$CDIR}

. $RUNTOOLS/lib/function_4.sh
. $RUNTOOLS/lib/function_4_all.sh

cd $TMPDIR
mpmd_arg=""
for member in \$(seq  $ENSEMBLE_START $ENSEMBLE_END ) ; do
   nnn=\$(getmember_extension \$member nodot)
   mmm=.\$nnn
   cat << eof1 > zbxiosens
#!/bin/bash

   cd $DDIR/${CONFIG_CASE}-XIOS.$ext/\$nnn
   ln -sf $P_UTL_DIR/bin/rebuild_nemo .
   ln -sf $P_UTL_DIR/bin/rebuild_nemo.exe .
   ln -sf ../$CN_DOMCFG ./

   for freq in 1ts 1h 1d 5d ; do
      mkdir @{freq}_OUTPUT
      for grid in gridT gridU gridV gridW icemod flxT; do
         # look for the number of server
         # check 0000.nc file alwats existing
         for f in @(ls ${CONFIG_CASE}\${mmm}_@{freq}_@{grid}_*-*_0000.nc) ; do
            fbase=@{f%_0000.nc}
            tag=@{fbase##*-}
            ndastp=y@{tag:0:4}m@{tag:4:2}d@{tag:6:2}.@freq
            echo @fbase @tag @ndastp
            nio=@( ls -1 @{fbase}_[[:digit:]]*[[:digit:]].nc  | wc -l )
            echo @nio for  ${CONFIG_CASE}\${mmm}_@{freq} @{grid}
            if [ @nio != 0 ] ; then
              ./rebuild_nemo @fbase @nio
              mv @fbase.nc @{freq}_OUTPUT/${CONFIG_CASE}\${mmm}_@{ndastp}_@grid.nc
            fi
         done
      done
   done
eof1
   cat zbxiosens | sed -e 's/@/\$/g' > ./$1\${mmm}.sh
   chmod 755 ./$1\${mmm}.sh                         # made it executabl  
   mpmd_arg="\$mpmd_arg 1 ./$1\${mmm}.sh"           # prepare the command line for runcode function

done # members
   pwd
   runcode_mpmd \$mpmd_arg
eof
  # copy the script to P_CTL_DIR from where it will be launched by submit
  copy $1 $P_CTL_DIR
          }
# ---

# Prepare a script for merging files (using mergeproc- off-line).
# this script is valid for both ensemble and non ensemble run
  mkbuild_merge() {
  mk_batch_hdr --name ${1%%.*} --cores $NB_NPROC_MER --wallclock $WALL_CLK_MER \
             --account $ACCOUNT --cluster hpt  --adapp --queue $QUEUE \
             --nodes $NB_NNODE_MER --constraint $CONSTRAI_MER --option '--exclusive' > $1
       echo "  *** building merging script "
       cat  << eof >> $1
        set -x 
        ulimit -s unlimited
      . $RUNTOOLS/lib/function_4.sh
      . $RUNTOOLS/lib/function_4_all.sh
         DDIR=${DDIR:-$CDIR}
         zXIOS=$DDIR/${CONFIG_CASE}-XIOS.$ext
         mergeprog=$(basename $MERGE_EXEC )
         cd \$zXIOS
         # deal with scalar files
         ls *scalar*0000.nc > /dev/null  2>&1
         if [ \$? = 0 ] ; then
            mkdir -p SCALAR
            mv *scalar*.nc SCALAR
            cd SCALAR
              for f in *scalar*_0000.nc ; do
                 CONFCASE=\$( echo \$f | awk -F_ '{print \$1}' )
                 freq=\$( echo \$f | awk -F_ '{print \$2}' )
                 tag=\$( echo \$f | awk -F_ '{print \$5}' | awk -F- '{print \$1}' )
                 date=y\${tag:0:4}m\${tag:4:2}d\${tag:6:2}

                 g=\${CONFCASE}_\${date}.\${freq}_icescalar.nc
                 OUTDIR=../\${freq}_OUTPUT
                 mkdir -p \$OUTDIR
                 cp \$f \$OUTDIR/\$g

              done
            cd  \$zXIOS

         # end scalar file
         fi
         ln -sf $MERGE_EXEC ./
             runcode $NB_NPROC_MER ./\$mergeprog -F -c $CN_DOMCFG -r
eof
  copy $1 $P_CTL_DIR
            }
# ---


# Merge splitted netcdf files (produced by XIOS or standard NEMO output, on the fly
#  Take care that this program must be run on depopulated nodes ( intensive I/O) 
#   one or 2 by computing node. So far, we are on the fly in a NEMO production and we have
#   plenty of computing node available, but be sure to place them correctly
mergefiles()  {
         DDIR=${DDIR:-$CDIR}
         cd $DDIR/${CONFIG_CASE}-XIOS.$ext/
         ln -sf $MERGE_EXEC ./
         mergeprog=$(basename $MERGE_EXEC )
         date
         # list of file can be very long and some machine have limitation on the 
         # length of an argument (eg, IBM has a 24kb limit). A workaround for this
         # problem is to fraction the list of file to deal with in smaller chunks.
         # nlenmax is set to 24000 ( IBM example, but can be adjusted according to 
         # the limitation of the system you are using

         getlst0000 24000
         idxmax=$(( ${#lst0000[@]} - 1 )) #  index max in the list, starting from 0
         
         for idx in $( seq 0 $idxmax ) ; do 
            runcode_u $NB_NPROC_MER ./$mergeprog -f ${lst0000[$idx]} -c $CN_DOMCFG -r 
         done

         date
         cd $TMPDIR  # for the remnant of the script
              }
# ---

# Merge splitted netcdf files (produced by XIOS or standard NEMO output, on the fly
#  (idem as mergfile but for ensemble runs )
mergefiles_ens()  {
         DDIR=${DDIR:-$CDIR}
         zXIOS=$DDIR/${CONFIG_CASE}-XIOS.$ext
         WKDIR=$zXIOS/WRK.$$
         mkdir -p $WKDIR
         cd $WKDIR
         mergeprog=$(basename $MERGE_EXEC )
         # link all files in all member in a single dir
         for member in $(seq  $ENSEMBLE_START $ENSEMBLE_END ) ; do
            nnn=$(getmember_extension $member nodot)
            mmm=.$nnn
            ln -sf $zXIOS/$nnn/${CONFIG_CASE}*.nc ./
         done
         ln -sf  $zXIOS/$CN_DOMCFG ./
         ln -sf $MERGE_EXEC ./
         getlst0000 24000
         idxmax=$(( ${#lst0000[@]} - 1 )) #  index max in the list, starting from 0

         for idx in $( seq 0 $idxmax ) ; do
             runcode_u $NB_NPROC_MER ./$mergeprog -f ${lst0000[$idx]} -c $CN_DOMCFG -r
         done
         \rm *.nc  # in WKDIR
         cd $TMPDIR  # for the remnant of the script
              }
# ---
# mk_post_process : function that build a script to launch zoomed area post processing in parallel
mk_post_process()  {
    cat << eof > $1
#!/bin/bash
# automatically created post_processing script from
# function mk_post_process in function_4_all.sh
          . $RUNTOOLS/lib/function_4_all.sh
          CONFIG=$CONFIG
          DDIR=$DDIR
          CONFIG_CASE=$CONFIG_CASE
          CASE=$CASE
          member=\$1
          ext=\$2

          echo MEMBER \$member
          lis_zoomid='$lis_zoomid'
          
          nnn=\$(getmember_extension \$member nodot)
          mmm=\$(getmember_extension \$member      )

          cd \$DDIR/\${CONFIG_CASE}-XIOS.\$ext/\$nnn
          if [ \$nnn ] ; then
             ln -sf ../$CN_DOMCFG ./
             ln -sf ../*xml ./
          fi
          for zoomid in \$lis_zoomid ; do
             post_process_one_file \$zoomid
          done
eof
chmod 755 $1
                   }


# ---
# post_process_one_file : when using XIOS in one_file mode, need to rename the file DRAKKAR style, and to
#   correct nav_lon,nav_lat for masked land domain ( ini_mpp2)
# also used for one_file output for regional zoom. In this latter case, the call is made with the zoomid
post_process_one_file()  {
         if [ $# = 1 ] ; then zoomid=$1 ; zcoord=${zoomid}_coordinates.nc ; else zoomid=  ; zcoord="$CN_DOMCFG"  ; fi

         DDIR=${DDIR:-$CDIR}
         ls -ld WRK.* > /dev/null 2>&1
         if [ $? = 0 ] ; then
            ztmp=../WRK.*
         else
            ztmp=.
         fi
         # mkdir <freq>_OUTPUT directories according to existing files
         for freq in 1ts 1h 3h 1d 3d 5d 1m 1mo ; do
#            ls *${freq}*_????????-????????.nc  > /dev/null 2>&1 
            ls *${freq}*_*-*.nc  > /dev/null 2>&1 
            if [ $? = 0 ] ; then 
               mkdir -p ${ztmp}/${freq}_OUTPUT
            fi
         done
         # check if zoom coordinates are there in case of zoom
         mkdir -p $WORKDIR/${CONFIG}_Sections  
         if [ $zoomid ] ; then
               if [ ! -f $WORKDIR/${CONFIG}_Sections/${zoomid}_coordinates.nc ] ; then
                  mk_zoom_coord $zoomid  
                  mv ${zoomid}_coordinates.nc $WORKDIR/${CONFIG}_Sections/
               fi
               ln -sf $WORKDIR/${CONFIG}_Sections/${zoomid}_coordinates.nc ./
         fi
         # check that nav_lonlat.nc files are here...
         for ztyp in T U V ; do
           if [ ! -f  $WORKDIR/${CONFIG}_Sections/${CONFIG}_${zoomid}nav_lonlat_$ztyp.nc ] ; then
             # build then from $CN_DOMCFG file
             mk_nav_lonlat $zcoord  $zoomid
             mv ${CONFIG}_${zoomid}nav_lonlat_$ztyp.nc  $WORKDIR/${CONFIG}_Sections/
           fi 
         done

         ln -sf $WORKDIR/${CONFIG}_Sections/${CONFIG}_${zoomid}nav_lonlat_*.nc ./

         # correct nav_lon, nav_lat for gridU files
         for f in ${CONFIG_CASE}*_${zoomid}grid*U_*-*.nc ; do
            ncks -h -A -v nav_lon,nav_lat ${CONFIG}_${zoomid}nav_lonlat_U.nc $f
            rename_out $f
         done
         for f in ${CONFIG_CASE}*_${zoomid}grid*Usurf_*-*.nc ; do
            ncks -h -A -v nav_lon,nav_lat ${CONFIG}_${zoomid}nav_lonlat_U.nc $f
            rename_out $f
         done


         # correct nav_lon, nav_lat for gridV files
         for f in ${CONFIG_CASE}*_${zoomid}gridV_*-*.nc ; do
            ncks -h -A -v nav_lon,nav_lat ${CONFIG}_${zoomid}nav_lonlat_V.nc $f
            rename_out $f
         done
         for f in ${CONFIG_CASE}*_${zoomid}gridVsurf_*-*.nc ; do
            ncks -h -A -v nav_lon,nav_lat ${CONFIG}_${zoomid}nav_lonlat_U.nc $f
            rename_out $f
         done
         # correct nav_lon, nav_lat for gridiU files (instantaneous output)
         for f in ${CONFIG_CASE}*_${zoomid}gridiU_*-*.nc ; do
            if [ -f $f ] ; then
              ncks -h -A -v nav_lon,nav_lat ${CONFIG}_${zoomid}nav_lonlat_U.nc $f
              rename_out $f
            fi
         done
         # correct nav_lon, nav_lat for gridiV files (instantaneous output)
         for f in ${CONFIG_CASE}*_${zoomid}gridiV_*-*.nc ; do
            if [ -f $f ] ; then
              ncks -h -A -v nav_lon,nav_lat ${CONFIG}_${zoomid}nav_lonlat_V.nc $f
              rename_out $f
            fi
         done

         # correct nav_lon, nav_lat for igridV files (intermember files)
         for f in ${CONFIG_CASE}*_${zoomid}igridV_*-*.nc ; do
            if [ -f $f ] ; then
              ncks -h -A -v nav_lon,nav_lat ${CONFIG}_${zoomid}nav_lonlat_V.nc $f
              rename_out $f
            fi
         done

         # after moving gridU and gridV files all remnant files are on gridT
         for f in ${CONFIG_CASE}*_${zoomid}gridT_*-*.nc  \
                  ${CONFIG_CASE}*_${zoomid}gridTsurf_*-*.nc  \
                  ${CONFIG_CASE}*_${zoomid}gridiT_*-*.nc  \
                  ${CONFIG_CASE}*_${zoomid}ptrcT_*-*.nc  \
                  ${CONFIG_CASE}*_${zoomid}icemod_*-*.nc \
                  ${CONFIG_CASE}*_${zoomid}flxT_*-*.nc \
                  ${CONFIG_CASE}*_${zoomid}gridW_*-*.nc ; do
            ncks -h -A -v nav_lon,nav_lat ${CONFIG}_${zoomid}nav_lonlat_T.nc $f
            rename_out $f
         done
         # cd TMPDIR required in the calling program
                         }
# ---
# rename_out : rename NEMO output file to <ztmp>/<freq>_OUTPUT/<CONFIG><zoomid>-<CASE>_<tag>_<type>.nc
#    ztmp is either . in case of standard output or the name of a WRK directory (ensemble run)
#    it takes the NEMO name as input
rename_out() {
            nemo_file=$1
            zfreq=$( echo $nemo_file       | awk -F_ '{ print $2}')
            if [ $zoomid ] ; then
              ztype=$( echo $nemo_file       | sed -e "s/$zoomid//"  | awk -F_ '{print $3}' )
            else
              ztype=$( echo $nemo_file       | awk -F_ '{print $3}' )
            fi
            zndastp=$(  echo ${nemo_file%.nc} | awk -F_ '{ print $NF}'| awk -F- '{print $2}' )
            zyy=${zndastp:0:4} ; zyy=${zyy:=0000}
            zmm=${zndastp:4:2} ; zmm=${zmm:=00}
            zdd=${zndastp:6:2} ; zdd=${zdd:=00}
            ztag=y${zyy}m${zmm}d${zdd}.$zfreq

            drak_file=${ztmp}/${zfreq}_OUTPUT/${CONFIG}${zoomid}-${CASE}${mmm}_${ztag}_$ztype.nc
            cp $nemo_file $drak_file
#            mv $nemo_file $drak_file
             }
# ---
# mk_nav_lonlat is a function that build 3 files having nav_lon, and nav_lat (2D files) for U V T points
#    It takes the information from $CN_DOMCFG file, and need to rename glamx, gphix to nav_lon,nav_lat
#    it also reduce the nav_lon, nav_lat variables to 2D (x,y) and transform then to float as in XIOS files.
mk_nav_lonlat()   {
    if [ $# -gt 0 ] ; then coord=$1  ; else  coord="$CN_DOMCFG" ; fi
    if [ $# = 2   ] ; then zoomid=$2 ; else  zoomid=                ; fi
    # some coordinates file have an unlimited dimensions other no ...
    unlim=''
    ncdump -h $coord | grep -iq UNLIMITED
    if [ $? = 0 ] ; then unlim=1 ; fi
    zvar='time'
    ncdump -h $coord | grep -iq 'z ='
    if [ $? = 0 ] ; then zvar='z,time' ; fi

    # T points
    ncks -O -h -v glamt,gphit $coord zzz.nc
    ncrename -h -v glamt,nav_lon -v gphit,nav_lat zzz.nc
    if [ $unlim ] ; then 
       ncwa -O -a $zvar zzz.nc ${CONFIG}_${zoomid}nav_lonlat_T.nc
    else
       mv zzz.nc ${CONFIG}_${zoomid}nav_lonlat_T.nc
    fi
    ncap2 -h -s 'nav_lon=float(nav_lon)' -s 'nav_lat=float(nav_lat)' ${CONFIG}_${zoomid}nav_lonlat_T.nc zozo 
    mv zozo ${CONFIG}_${zoomid}nav_lonlat_T.nc

    # U points
    ncks -O -h -v glamu,gphiu $coord zzz.nc
    ncrename -h -v glamu,nav_lon -v gphiu,nav_lat zzz.nc
    if [ $unlim ] ; then 
       ncwa -O -a $zvar zzz.nc ${CONFIG}_${zoomid}nav_lonlat_U.nc
    else
       mv zzz.nc ${CONFIG}_${zoomid}nav_lonlat_U.nc
    fi
    ncap2 -h -s 'nav_lon=float(nav_lon)' -s 'nav_lat=float(nav_lat)' ${CONFIG}_${zoomid}nav_lonlat_U.nc zozo 
    mv zozo ${CONFIG}_${zoomid}nav_lonlat_U.nc

    # U points
    ncks -O -h -v glamv,gphiv $coord zzz.nc
    ncrename -h -v glamv,nav_lon -v gphiv,nav_lat zzz.nc
    if [  $unlim ] ; then 
       ncwa -O -a $zvar zzz.nc ${CONFIG}_${zoomid}nav_lonlat_V.nc
    else
       mv zzz.nc ${CONFIG}_${zoomid}nav_lonlat_V.nc
    fi
    ncap2 -h -s 'nav_lon=float(nav_lon)' -s 'nav_lat=float(nav_lat)' ${CONFIG}_${zoomid}nav_lonlat_V.nc zozo 
    mv zozo ${CONFIG}_${zoomid}nav_lonlat_V.nc
                  }
# ---
#  extract the zoom coordinated from the global coordinates file
#
mk_zoom_coord() {
    zoomid=$1
    domain="domain_def.xml"

    # now look for position of zoom in domain_ref.xml
    zoom=${zoomid}T   # JM Trick : suppose that the domain_ref ${zoomid}T exists (likely to be true)
    tmp=$(grep -i $zoom $domain)
    for f in $tmp ; do
       cmd=$(echo $f | grep -v '<' | grep -v '>')
       if [ $cmd ] ; then eval $cmd ; fi
    done
    zoom_iend=$(( zoom_ibegin + zoom_ni - 1 ))
    zoom_jend=$(( zoom_jbegin + zoom_nj - 1 ))

    zoomcoord=${zoomid}_coordinates.nc
    if [ ! -f $zoomcoord ] ; then
        ncks -F -O -d x,$zoom_ibegin,$zoom_iend -d y,$zoom_jbegin,$zoom_jend $CN_DOMCFG ${zoomcoord}
    fi
                }

# ---
# get the list of *0000.nc files in the current directory fractionned in sub-list whose
# length does not exceed the lencght (bytes) passed as argument
getlst0000()  {
         nlenmax=$1  # proxy for 24kb IBM limitation
         zlst=''
         idx=0   # index in the list
         for f in *${CONFCASE}*0000.nc ; do
            zlst="$zlst $f"
            len=$( echo $zlst | wc -c )
            if [ $len -gt $nlenmax ] ; then
              lst0000[$idx]=$zlst
              idx=$(( idx + 1 ))
              zlst=''
            fi
         done
         lst0000[$idx]=$zlst   # last chunk
              }
# ---

# Make batch header for submitted scripts ( common part to all machines). Called by specific machine function in function_4_machine.sh (mk_batch_hdr) with same arguments
# mk_batch_hdr_core  --name name --wallclock wallclock --account account --nodes nodes --cores cores --adapp --queue qname --cluster cname --par --seq --par  --option "options line" --help --memory memory
mk_batch_hdr_core () {

   while (( $# > 0 )) ; do
     case $1 in
     ("--name"      ) shift ; name=$1        ; shift ;;
     ("--wallclock" ) shift ; wallclock=$1   ; shift ;;
     ("--account"   ) shift ; account=$1     ; shift ;;
     ("--nodes"     ) shift ; nodes=$1       ; shift ;;
     ("--cores"     ) shift ; cores=$1       ; shift ;;
     ("--cluster"   ) shift ; cluster=$1     ; shift ;;
     ("--memory"    ) shift ; memory=$1      ; shift ;;
     ("--par"       ) shift ; jobtype='parallel' ;;
     ("--seq"       ) shift ; jobtype='serial'   ;;
     ("--option"    ) shift ; option=$1      ; shift ;;
     ("--queue"     ) shift ; queue=$1       ; shift ;;
     ("--adapp"     ) shift ; adapp=1            ;;
     ("--constraint") shift ; constraint=$1  ; shift ;;
     ("--help"      ) shift ;
         echo USAGE : mk_batch_hdr_core  --name name --wallclock wallclock --account account --nodes nodes  ... ;
         echo "       ... " --cores cores --par --seq --adapp --queue qname --option "options line" --constraint=constaint --help ; return ;;
     ( * )           shift ;;   # silently skip unknown options
     esac
   done
                    } 
#--------------------------------------------------------------------------------------
# AGRIF STUFF
# Init Agrif : set usefull variables in case of AGRIF run ( ARRAY agrif_pref )
initagrif() { 
   nbr_child=0 ;   
   rcopy $P_CTL_DIR/AGRIF_FixedGrids.in AGRIF_FixedGrids.in
   nbr_child=$(cat AGRIF_FixedGrids.in | awk 'BEGIN{s=0} NF==1 {s=s+$1} END{ print s}')
   n=1 ; lst=''
   while (( n <=  nbr_child )) ; do
     lst="$lst ${n}"
      n=$(( n + 1 ))
   done
   agrif_pref=( $lst )

   refinement=$( cat AGRIF_FixedGrids.in | awk ' NF==7 { print $7}' )
   timeref=( '  ' $refinement )
            } 
# ---
update_db_file()  {
     nop1=$(( $no + 1 ))
     # Procedure to change nn_itend during a run : If NEMO version supports it ( DRAKKAR enhancement)
     # creating a file, named nitend.txt with new desired value for nitend in the TMPDIR_CONFCASE directory,
     # will produce a modification in the stream of the run, changing nn_itend to the next possible choice
     # according to nn_write.  Doing so one can make a run shorter or longer ( In this later case, however,
     # there might be missing forcing files or OBC files ... )
     # check if 'nitend.txt' was used :
     if [ -f nitend.txt ] ; then  # change nitend by its real value
       echo modified nitend on the fly, correct db file
       newval=$(cat nitend.txt)
       nit000=`tail -1 $CONFIG_CASE.db | awk '{print $2}' `
       nitend=`tail -1 $CONFIG_CASE.db | awk '{print $3}' `
       nwrite=$(LookInNamelist nn_write namelist)
       # newnitend is computed as in NEMO
       newnitend=$( echo $nit000 $newval $nwrite | awk '{ print $1 - 1 + int( ( $2 - $1 ) / $3 +1 ) * $3 }')
       sed -e "s/ $nitend\$/ $newnitend/" $CONFIG_CASE.db > tmpdb
       mv tmpdb $CONFIG_CASE.db
       # also update the current namelist for restart renaming to work
       sed  -e "/nn_itend/s/$nitend/ $newnitend/" namelist > znamelist
       mv znamelist namelist
       \rm nitend.txt  # safe remove for next run 
     fi

     # add last date at the current line
     nline=$(wc $CONFIG_CASE.db | awk '{print $1}')

     # aammdd is the ndastp of the last day of the run ...
     # where can we get it ???? : in the ocean.output for sure !!
     aammdd=$( cat $output_ref | grep date | tail -1 | awk '{print $NF}' )

    # Look for line in db file  with only 3 columns, keep this line in last
    last=$( cat $CONFIG_CASE.db | awk ' NF == 3 ' )
    ncol=$( echo $last | wc -w )  # ncol = 0 means no 3 colums line found
    if [ $ncol = 0 ] ; then
      echo "db file is up to date with respect to date $aammdd"
    else
      sed -e "s/$last/$last\ $aammdd/" $CONFIG_CASE.db > tmpdb
      mv -f tmpdb $CONFIG_CASE.db
    fi

    # add a new last line for the next run
    nstep_per_day=$(( 86400 / $rdt ))

   if [ $ndays = 185 ] ; then
     dif=$(( 180 * $nstep_per_day ))
   elif [ $ndays = 180 ] ; then
     dif=$(( 185 * $nstep_per_day ))
   else
  nit000=`tail -1 $CONFIG_CASE.db | awk '{print $2}' `
  nitend=`tail -1 $CONFIG_CASE.db | awk '{print $3}' `
#     dif=$(( 365 * $nstep_per_day ))
     dif=$((  $nitend - $nit000  + 1  ))
   fi

    nit000=$(( $nitend + 1 ))
    nitend=$(( $nitend + $dif ))

    # offer the opportunity to modify last line of db file on the fly: use newdb last line if any
    if [ -f newdb ] ; then
      line=$( cat newdb )
      echo $line >> $CONFIG_CASE.db
      \rm newdb
    elif [ -f newdb.$nop1 ] ; then
      line=$( cat newdb.$nop1 )
      echo $line >> $CONFIG_CASE.db
    else
      echo $nop1 $nit000 $nitend >> $CONFIG_CASE.db
    fi
                  }

######################################################################################
