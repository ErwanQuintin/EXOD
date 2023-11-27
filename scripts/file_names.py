"""
EXODUS - EPIC XMM-Newton Outburst Detector Ultimate System           
Declaration of file names                                            
Maitrayee Gupta (2022) - maitrayee.gupta@irap.omp.eu                 
Declaration of the file names handled both by the detector and the renderer
"""
import os
from logger import logger


# Created files
LOG = "log.txt"
VARIABILITY = "variability_file.fits"
REGION = "ds9_variable_sources.reg"
OUTPUT_IMAGE = "variability.pdf"
OUTPUT_IMAGE_SRCS = "sources.pdf"
OUTPUT_IMAGE_ALL = "variability_whole.pdf"
OUTPUT_EXODUS = "variability_all_inst.pdf"
BEST_MATCH = "best_match.csv"

# Observation files
# PN
CLEAN_FILE_PN = "PN_clean.fits"
GTI_FILE_PN = "PN_gti.fits"
IMG_FILE_PN = "PN_image.fits"
RATE_FILE_PN = "PN_rate.fits"

# MOS 1
CLEAN_FILE_M1 = "M1_clean.fits"
GTI_FILE_M1 = "M1_gti.fits"
IMG_FILE_M1 = "M1_image.fits"
RATE_FILE_M1 = "M1_rate.fits"

# MOS 2
CLEAN_FILE_M2 = "M2_clean.fits"
GTI_FILE_M2 = "M2_gti.fits"
IMG_FILE_M2 = "M2_image.fits"
RATE_FILE_M2 = "M2_rate.fits"

# software installation paths

HEADAS = os.environ["HEADAS"]
SAS = "/home/erwan/Documents/Softwares/sas_19.0.0-Ubuntu18.04/xmmsas_20201028_0905/setsas.sh"  # "/mnt/data/Maitrayee/EXOD/scripts/setsas.sh"

# scripts folder
FOLDER = "/mnt/data/Maitrayee/EXOD/data"
SCRIPTS = "/mnt/data/Maitrayee/EXOD/scripts"
