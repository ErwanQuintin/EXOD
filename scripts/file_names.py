"""
EXODUS - EPIC XMM-Newton Outburst Detector Ultimate System
Maitrayee Gupta (2022) - maitrayee.gupta@irap.omp.eu

Declaration of the file names handled both by the detector and the renderer.
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

# Software installation paths
# The SAS_PATH enviroment should be set to sas_19.0.0-Ubuntu18.04/xmmsas_20201028_0905
# or equivilent (NOTE NO TRAILING /)
HEADAS = os.environ["HEADAS"]
EXOD   = os.environ["EXOD"]
SAS_PATH = os.environ["SAS_PATH"]
SAS = f"{SAS_PATH}/setsas.sh"

# scripts folder
FOLDER = f"{EXOD}/data"  # this points to the data folder but it called FOLDER?!
SCRIPTS = f"{EXOD}/scripts"
