"""
EXODUS - EPIC XMM-Newton Outburst Detector Ultimate System
Maitrayee Gupta (2022) - maitrayee.gupta@irap.omp.eu      

Various resources for both detector and renderer.
"""
import os
import file_names as FileNames
from logger import logger

def open_files(folder_name):
    """
    Function opening files and writing their legend.
    @param  folder_name:  The directory to create the files.
    @return: info_file, variability_file, counter_per_tw,
    detected_variable_areas_file, time_windows_file
    """

    # Fixing the name of the folder
    if folder_name[-1] != "/":
        folder_name += "/"

    os.makedirs(folder_name, exist_ok=True)

    var_file = folder_name + FileNames.VARIABILITY
    reg_file = folder_name + FileNames.REGION
    best_match_file = folder_name + FileNames.BEST_MATCH
    return var_file, reg_file, best_match_file

