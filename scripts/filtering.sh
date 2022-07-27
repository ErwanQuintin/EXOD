#!/bin/bash

########################################################################
#                                                                      #
# EXODUS - EPIC XMM-Newton Outburst Detector Ultimate System           #
#                                                                      #
# Events file filtering                                                #
#                                                                      #
# Maitrayee Gupta (2022) - maitrayee.gupta@irap.omp.eu                 #
# Inés Pastor Marazuela (2019) - ines.pastor.marazuela@gmail.com      #
#                                                                      #
########################################################################

###
# Parsing arguments                                                            
###

# Default variables
RATE=0.5	# (PN Rate = 0.5 and MOS Rate = 0.4)
FOLDER=/home/mike/sas/xmmsas_20210317_1624/EXOD/data
SCRIPTS=/home/mike/sas/xmmsas_20210317_1624/EXOD/scripts
INST=PN

# Input variables
while [[ $# -gt 0 ]]; do
case "$1" in
  -o|-obs|--observation) OBS=${2}
  shift; shift ;;
  -r|--rate)             RATE=${2:-$RATE}
  shift; shift ;;
  -i|--instrument)       INST=${2:-$INST}
  shift; shift ;;
  -f|--folder)           FOLDER=${2:-$FOLDER}
  shift; shift ;;
  -s|--scripts)          SCRIPTS=${2:-$SCRIPTS}
  shift; shift ;;
esac
done

echo -e "\tFOLDER      = ${FOLDER}"
echo -e "\tOBSERVATION = ${OBS}"
echo -e "\tINSTRUMENT  = ${INST}"
echo -e "\tRATE        = ${RATE}"


path=$FOLDER/$OBS

###
# Defining functions
###

Title(){
  message=$1; i=0; x='===='
  while [[ i -lt ${#message} ]]; do x='='$x; ((++i)); done
  echo -e "\n\t  $message \n\t$x"
}

title(){
  message=$1; i=0; x=----
  while [[ i -lt ${#message} ]]; do x=-$x; ((++i)); done
  echo -e "\n\t  $message \n\t$x"
}

# Useful
########################################################################

var(){
  x=$1
  out=$(cat $SCRIPTS/file_names.py | grep ^$x | awk '{print $3}' | sed 's/"//g')
  echo $out
}

########################################################################
#                                                                      #
# Main programme                                                       #
#                                                                      #
########################################################################

Title "Filtering observation $OBS"

###
# Preliminaries
###

title "Preliminaries"

# Setting up SAS
cd $path
export SAS_ODF=$path
export SAS_CCF=$path/ccf.cif
export HEADAS=$(var HEADAS)
. $HEADAS/headas-init.sh
. $(var SAS)

if [ ! -f $path/ccf.cif ]; then cifbuild; fi

###
# Filtering
###

title "Cleaning events file"

# File names

org_file=$(ls $path/*${OBS}${INST}*IEVLI*)
clean_file=$path/${INST}_clean.fits
gti_file=$path/${INST}_gti.fits
img_file=$path/${INST}_image.fits
rate_file=$path/${INST}_rate.fits


echo -e "\tRAW FILE   = ${org_file}"
echo -e "\tCLEAN FILE = ${clean_file}"
echo -e "\tGTI FILE   = ${gti_file}"
echo -e "\tIMAGE FILE = ${img_file}"
echo -e "\tRATE FILE  = ${rate_file}"

# Creating GTI
title "Creating GTI"

if [ "$INST" == "PN" ]; then

  evselect table=$org_file withrateset=Y rateset=$rate_file maketimecolumn=Y timebinsize=100 makeratecolumn=Y expression='#XMMEA_EP && (PI in [10000:12000]) && (PATTERN==0)' -V 0
  RATE=0.5

elif [ "$INST" == "M1" ] || [ "$INST" == "M2" ]; then

  evselect table=$org_file withrateset=Y rateset=$rate_file maketimecolumn=Y timebinsize=100 makeratecolumn=Y expression='#XMMEA_EM && (PI>10000) && (PATTERN==0)' -V 0
  RATE=0.4

fi

if [[ $RATE != [0-9]* ]]; then
    echo "Opening ${INST}_rate.fits" 
    fv $rate_file &
    read -p "Choose the GTI cut rate : " RATE
fi
echo "Creating Good Time Intervals with threshold RATE=$RATE"

tabgtigen table=$rate_file expression="RATE<=$RATE" gtiset=$gti_file -V 0

# Cleaning events file

if [ "$INST" == "PN" ]; then

  evselect table=$org_file withfilteredset=Y filteredset=$clean_file destruct=Y keepfilteroutput=T expression="#XMMEA_EP && gti($gti_file,TIME) && (PATTERN<=4) && (PI in [500:12000])" -V 0

elif [ "$INST" == "M1" ] || [ "$INST" == "M2" ]; then

  evselect table=$org_file withfilteredset=Y filteredset=$clean_file destruct=Y keepfilteroutput=T expression="#XMMEA_EM && gti($gti_file,TIME) && (PATTERN<=12) && (PI in [500:10000])" -V 0

fi

# Creating image file

evselect table=$clean_file imagebinning=binSize imageset=$img_file withimageset=yes xcolumn=X ycolumn=Y ximagebinsize=80 yimagebinsize=80 -V 0


echo "Rate = $RATE" >> $path/${INST}_processing.log

echo "The end" 
date 
