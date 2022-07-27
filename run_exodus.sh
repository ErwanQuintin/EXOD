#!/bin/bash

########################################################################
#                                                                      #
# EXODUS - EPIC XMM-Newton Outburst Detector Ultimate System           #
#                                                                      #
# Main Shell script to run the program                                 #
#                                                                      #
# Maitrayee Gupta (2022) - maitrayee.gupta@irap.omp.eu                 #
#                                                                      #
########################################################################


###
# Parsing arguments                                                            
###

# Default variables
LC_ENABLE=false	# by default do not generate light curves
ITERATIVE_DET_LEVEL=false

# Input variables
while [[ $# -gt 0 ]]; do
case "$1" in
  --obs_id_file|-obs)          OBS_ID_FILE=${2}
  shift; shift ;;
  --lc_enable|-lc)            LC_ENABLE=${2:-$LC_ENABLE}
  shift; shift ;;
  --iterative_det_level|-it)  ITERATIVE_DET_LEVEL=${2:-$ITERATIVE_DET_LEVEL}
esac
done

echo -e "\tOBS_ID_FILE = ${OBS_ID_FILE}"
echo -e "\tLC_ENABLE  = ${LC_ENABLE}"
echo -e "\tITERATIVE_DET_LEVEL = ${ITERATIVE_DET_LEVEL}"

# clean up any residual files from previous run
rm -rf  Master_Catalogue.csv match_countM1.txt match_countM2.txt  match_count_pn.txt  match_countPN.txt  separation_file.csv triple_match.csv  triple_match_obs.csv  

PN_START_LEVEL=6
M_START_LEVEL=4
if $ITERATIVE_DET_LEVEL
then
    PN_START_LEVEL=8
    M_START_LEVEL=8
fi



TW=100
names=$(<$OBS_ID_FILE)
path_prefix="/mnt/data/Maitrayee/EXOD/data/"
scripts_path="/mnt/data/Maitrayee/EXOD/scripts/"
for obs in $names
do
    echo $obs
    # clean up any residual files from previous run
    rm -rf  Master_Catalogue.csv match_countM1.txt match_countM2.txt  match_count_pn.txt  match_countPN.txt  separation_file.csv
    obs_path="$path_prefix$obs"
    echo $obs_path
    bash $SCRIPTS/download_observation.sh $FOLDER $obs PN MOS1 MOS2
    bash $SCRIPTS/filtering.sh -f $FOLDER -o $obs -s $SCRIPTS -i PN
    bash $SCRIPTS/filtering.sh -f $FOLDER -o $obs -s $SCRIPTS -i M1
    bash $SCRIPTS/filtering.sh -f $FOLDER -o $obs -s $SCRIPTS -i M2
    
    pn_det_level=$PN_START_LEVEL
    final_pn_det_level=$PN_START_LEVEL
    while [ $pn_det_level -gt 5 ]
    do
       echo $pn_det_leve
    
       python3 $SCRIPTS/detector.py -path $obs_path -bs 3 -dl $pn_det_level -tw $TW -gtr 1.0 -mta 16 -obs $obs --render -inst PN
       python3 $SCRIPTS/detect_pn_matches.py -path $obs_path  -bs 3 -pn_dl $pn_det_level -tw $TW -gtr 1.0 
       match_count=$(<match_count_pn.txt)
       echo $match_count
       if [ $match_count -gt 0 ]
       then
    	  final_pn_det_level=$pn_det_level
          pn_det_level=0
       else
    	       final_pn_det_level=$pn_det_level
	       let "pn_det_level=pn_det_level-1" 
       fi
       echo $pn_det_level

    done

  # remove failed directories
    for (( c= $final_pn_det_level+1 ; c<=$PN_START_LEVEL; c++ ))
	do  
		prefix=$obs_path
		slash="/"
		under_score="_"
		tw_name=$TW
		suffix="_3_1.0_PN"
                f_name="$prefix$slash$c$under_score$TW$suffix"
		rm -rf $f_name
	done

    
    m1_det_level=$M_START_LEVEL
    final_m1_det_level=$M_START_LEVEL
    while [ $m1_det_level -gt 3 ]
    do
       echo $m1_det_level
       python3 $SCRIPTS/detector.py -path  $obs_path  -bs 13 -dl $m1_det_level -tw $TW -gtr 1.0 -mta 16 -obs $obs --render -inst M1
       python3 $SCRIPTS/detect_corr.py -path $obs_path   -bs 13 -pn_dl $final_pn_det_level -m_dl $m1_det_level -tw $TW -gtr 1.0 -inst M1
       match_count=$(<match_countM1.txt)
       echo $match_count
       if [ $match_count -gt 0 ]
       then
    	  final_m1_det_level=$m1_det_level
          m1_det_level=0
       else
	       final_m1_det_level=$m1_det_level
	       let "m1_det_level=m1_det_level-1" 
       fi
       echo $m1_det_level

    done	

  # remove failed directories
    for (( c= $final_m1_det_level+1 ; c<=$M_START_LEVEL; c++ ))
	do
    		prefix=$obs_path
                slash="/"
		under_score="_"
                tw_name=$TW
                suffix="_13_1.0_M1"
                f_name="$prefix$slash$c$under_score$TW$suffix"
                rm -rf $f_name
	done


    m2_det_level=$M_START_LEVEL
    final_m2_det_level=$M_START_LEVEL
    while [ $m2_det_level -gt 3 ]
    do
       echo $m2_det_level
       python3 $SCRIPTS/detector.py -path  $obs_path  -bs 13 -dl $m2_det_level -tw $TW -gtr 1.0 -mta 16 -obs $obs --render -inst M2
       python3 $SCRIPTS/detect_corr.py -path $obs_path  -bs 13 -pn_dl $final_pn_det_level -m_dl $m2_det_level -tw $TW -gtr 1.0 -inst M2
       match_count=$(<match_countM2.txt)
       echo $match_count
       if [ $match_count -gt 0 ]
       then
    	  final_m2_det_level=$m2_det_level
          m2_det_level=0
       else
	       final_m2_det_level=$m2_det_level
	       let "m2_det_level=m2_det_level-1" 
       fi
       echo $m2_det_level

    done

  # remove failed directories
  for (( c= $final_m2_det_level+1 ; c<=$M_START_LEVEL; c++ ))
	do  
		prefix=$obs_path
		slash="/"
		under_score="_"
                tw_name=$TW
		suffix="_13_1.0_M2"
                f_name="$prefix$slash$c$under_score$TW$suffix"
		rm -rf $f_name
	done


    echo $final_pn_det_level 
    echo $final_m1_det_level
    echo $final_m2_det_level
    python3 $SCRIPTS/calc_all_separations.py -path $obs_path  -bs 13 -pn_dl $final_pn_det_level  -m1_dl $final_m1_det_level -m2_dl  $final_m2_det_level -obs $obs -tw $TW -gtr 1.0
    python3 $SCRIPTS/detect_final_corr.py -path $obs_path  -bs 13 -pn_dl $final_pn_det_level  -m1_dl $final_m1_det_level -m2_dl  $final_m2_det_level -obs $obs -tw $TW -gtr 1.0
   
    # clean up any residual files for future run
    rm -rf  Master_Catalogue.csv match_countM1.txt match_countM2.txt  match_count_pn.txt  match_countPN.txt  separation_file.csv

    if $LC_ENABLE
    then

    	python3 $SCRIPTS/pre_light_curve.py -path  $obs_path -bs 13 -pn_dl $final_pn_det_level -m1_dl $final_m1_det_level -m2_dl $final_m2_det_level -obs $obs -tw $TW -gtr 1.0
	pn_match_count=$(<match_countPN.txt)
	m1_match_count=$(<match_countM1.txt)
	m2_match_count=$(<match_countM2.txt)
	id=1
	while [ $id -lt $pn_match_count  ]
	do
		bash $SCRIPTS/lightcurve.sh   $obs_path  $SCRIPTS PN $id  $final_pn_det_level $TW 1.0 3  $obs_path/LC_LOG_$TW_PN_$id
		let "id=id+1"
	done

	id=1
	while [ $id -lt $m1_match_count  ]
	do
		bash $SCRIPTS/lightcurve.sh   $obs_path  $SCRIPTS M1 $id  $final_m1_det_level $TW 1.0 13  $obs_path/LC_LOG_$TW_M1_$id
		let "id=id+1"
	done

	id=1
	while [ $id -lt $m2_match_count  ]
	do
		bash $SCRIPTS/lightcurve.sh   $obs_path  $SCRIPTS M2 $id  $final_m2_det_level $TW 1.0 13  $obs_path/LC_LOG_100_M2_$id
		let "id=id+1"
	done
	# clean up any residual files for future run
        rm -rf  Master_Catalogue.csv match_countM1.txt match_countM2.txt  match_count_pn.txt  match_countPN.txt  separation_file.csv
 
     fi

done
echo Run done

python3 $SCRIPTS/match_with_xmm_catalog.py -path /mnt/data/Maitrayee/EXOD/
python3 $SCRIPTS/remove_bright_sources.py -path /mnt/data/Maitrayee/EXOD/
python3 $SCRIPTS/match_with_simbad.py
python3 $SCRIPTS/simbad_subclass.py


# clean up any residual files for future run
rm -rf  Master_Catalogue.csv match_countM1.txt match_countM2.txt  match_count_pn.txt  match_countPN.txt  separation_file.csv triple_match.csv  triple_match_obs.csv  
echo Post process done
   
