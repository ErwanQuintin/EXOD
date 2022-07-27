#!/bin/bash

########################################################################
#                                                                      #
# EXODUS - EPIC XMM-Newton Outburst Detector Ultimate System           #
#                                                                      #
# Download Observauion files for OBSID                                 #
#                                                                      #
# Maitrayee Gupta (2022) - maitrayee.gupta@irap.omp.eu                 #
# Inés Pastor Marazuela (2019) - ines.pastor.marazuela@gmail.com      #
#                                                                      #
########################################################################

# Usage
# bash $SCRIPTS/download_observation.sh $FOLDER $obs PN MOS1 MOS2
#
#FOLDER       - where to download the observation data
#obsid        - XMM OBSID to download
#PN MOS1 MOS2 - One or more intrument names

folder=$1
obs=$2
mode=(${3:-"PN"} $4 $5) # Mode event list
path=$folder/$obs
if [ ! -d $path ]; then mkdir $path; fi

cd $path
echo ${mode[@]}
###
# Observation
###

#PN
if [[ "${mode[@]}" =~ "PN" ]]; then
  echo "Downloading PN data ..."
  curl -o $path/P${obs}PNS001PIEVLI.FTZ "http://nxsa.esac.esa.int/nxsa-sl/servlet/data-action-aio?obsno=${obs}&instname=PN&level=PPS&name=PIEVLI"
  curl -o $path/P${obs}PNS001FBKTSR0000.FTZ "http://nxsa.esac.esa.int/nxsa-sl/servlet/data-action-aio?obsno=${obs}&name=FBKTSR&instname=PN&level=PPS&extension=FTZ"
fi

# MOS1
if [[ "${mode[@]}" =~ "M1" ]] || [[ "${mode[@]}" =~ "MOS" ]] ; then
  echo "Downloading M1 data ..."
  curl -o $path/P${obs}M1S002MIEVLI.FTZ "http://nxsa.esac.esa.int/nxsa-sl/servlet/data-action-aio?obsno=${obs}&instname=M1&level=PPS&name=MIEVLI"
  curl -o $path/P${obs}M1S002FBKTSR0000.FTZ "http://nxsa.esac.esa.int/nxsa-sl/servlet/data-action-aio?obsno=${obs}&name=FBKTSR&instname=M1&level=PPS&extension=FTZ"
fi

# MOS2
if [[ "${mode[@]}" =~ "M2" ]] || [[ "${mode[@]}" =~ "MOS" ]] ; then
  echo "Downloading M2 data ..."
  curl -o $path/P${obs}M2S003MIEVLI.FTZ "http://nxsa.esac.esa.int/nxsa-sl/servlet/data-action-aio?obsno=${obs}&instname=M2&level=PPS&name=MIEVLI"
  curl -o $path/P${obs}M2S003FBKTSR0000.FTZ "http://nxsa.esac.esa.int/nxsa-sl/servlet/data-action-aio?obsno=${obs}&name=FBKTSR&instname=M2&level=PPS&extension=FTZ"
fi

# Summary file
curl -o $path/sas.TAR "http://nxsa.esac.esa.int/nxsa-sl/servlet/data-action-aio?obsno=${obs}&level=ODF&extension=SAS"
tar -xf $path/sas.TAR -C $path
rm *ATS.FIT *TCS.FIT *RAS.ASC *ROS.ASC MANIFEST* *.TAR

