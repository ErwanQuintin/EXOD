"""
EXODUS - EPIC XMM-Newton Outburst Detector Ultimate System
Maitrayee Gupta (2022) - maitrayee.gupta@irap.omp.eu      

Various resources for both detector and renderer.
"""
import os
import subprocess
import sys
from math import *

import numpy as np
import scipy.ndimage as nd
import skimage.transform

import file_names as FileNames
from logger import logger

########################################################################
#                                                                      #
# Function and procedures                                              #
#                                                                      #
########################################################################


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





def read_from_file(file_path, counter=False, comment_token="#", separator=";"):
    """
    Function returning the content
    @param counter: True if it is the counters that are loaded, False if it is the variability
    @return: A list
    """
    data = []
    i = 0
    nb_lignes = -1

    with open(file_path) as f:
        for line in f:
            # if len(line) > 2 and comment_token not in line :
            if comment_token not in line:
                if i % 200 == 0:
                    data.append([])
                    nb_lignes += 1
                if not counter:
                    data[nb_lignes].append(float(line))
                else:
                    data[nb_lignes].append(
                        [float(tok) for tok in line.split(separator)]
                    )
                i += 1

    return data


########################################################################


def read_tws_from_file(file_path, comment_token="#", separator=";"):
    """
    Reads the list of time windows from its file.
    @return: A list of couples (ID of the TW, t0 of the TW)
    """
    tws = []

    with open(file_path) as f:
        for line in f:
            if len(line) > 2 and comment_token not in line:
                line_toks = line.split(separator)
                tws.append((int(line_toks[0]), float(line_toks[1])))

    return tws


########################################################################


########################################################################


def read_sources_from_file(file_path, comment_token="#", separator=";"):
    """
    Reads the source from their file
    @return: A list of Source objects
    """

    sources = []

    with open(file_path) as f:
        for line in f:
            if len(line) > 2 and comment_token not in line:
                toks = line.split(separator)
                sources.append(
                    Source(
                        toks[0],
                        int(toks[1]),
                        float(toks[2]),
                        float(toks[3]),
                        float(toks[4]),
                    )
                )

    return sources


########################################################################


def PN_config(data_matrix):
    """
    Arranges the variability data for EPIC_PN
    """
    data_v = []

    # XMM-Newton EPIC-pn CCD arrangement
    ccds = [[8, 7, 6, 9, 10, 11], [5, 4, 3, 0, 1, 2]]

    # Building data matrix
    for c in ccds[0]:
        data_v.extend(np.flipud(data_matrix[c]))
    i = 0
    for c in ccds[1]:
        m = np.flip(data_matrix[c])
        for j in range(64):
            data_v[i] = np.append(data_v[i], m[63 - j])
            i += 1
    return data_v


########################################################################


def M1_config(data_matrix):
    """
    Arranges the variability data for EPIC_MOS_1
    """
    corner = np.zeros((300, 600))

    # First part
    data_1 = np.concatenate(
        (corner, data_matrix[1].T, np.flip(data_matrix[6].T), corner), axis=0
    )

    # Second part
    data_2 = np.concatenate(
        (data_matrix[2].T, np.flipud(data_matrix[0]), np.flip(data_matrix[5].T)), axis=0
    )

    # Third part
    data_3 = np.concatenate(
        (corner, data_matrix[3].T, np.flip(data_matrix[4].T), corner), axis=0
    )

    # Building data matrix
    data_v = np.rot90(np.concatenate((data_1, data_2, data_3), axis=1))

    return data_v


########################################################################


def M2_config(data_matrix):
    """
    Arranges the variability data for EPIC_MOS_2
    """
    corner = np.zeros((300, 600))

    data_1 = np.concatenate(
        (corner, data_matrix[1].T, np.flip(data_matrix[6].T), corner), axis=0
    )

    data_2 = np.concatenate(
        (data_matrix[2].T, np.flipud(data_matrix[0]), np.flip(data_matrix[5].T)), axis=0
    )

    data_3 = np.concatenate(
        (corner, data_matrix[3].T, np.flip(data_matrix[4].T), corner), axis=0
    )

    # Building data matrix
    data_v = np.flip(np.concatenate((data_1, data_2, data_3), axis=1))

    return data_v


########################################################################
#                                                                      #
# Geometrical transformations                                          #
#                                                                      #
########################################################################


def data_transformation_PN(data, header):
    """
    Performing geometrical transformations from raw coordinates to sky coordinates
    @param data: variability matrix
    @param header: header of the clean events file
    @return: transformed variability data
    """

    # Header information
    angle = header["PA_PNT"]

    xproj = [float(header["TDMIN6"]), float(header["TDMAX6"])]  # projected x limits
    yproj = [float(header["TDMIN7"]), float(header["TDMAX7"])]  # projected y limits
    xlims = [float(header["TLMIN6"]), float(header["TLMAX6"])]  # legal x limits
    ylims = [float(header["TLMIN7"]), float(header["TLMAX7"])]  # legal y limits

    # scaling factor
    sx = 648 / (xlims[1] - xlims[0])
    sy = 648 / (ylims[1] - ylims[0])
    # pads (padding)
    padX = (int((xproj[0] - xlims[0]) * sx), int((xlims[1] - xproj[1]) * sx))
    padY = (int((yproj[0] - ylims[0]) * sy), int((ylims[1] - yproj[1]) * sy))
    # shape (resizing)
    pixX = 648 - (padX[0] + padX[1])
    pixY = 648 - (padY[0] + padY[1])

    # Transformations
    ## Rotation
    dataR = np.flipud(nd.rotate(data, angle, reshape=True))
    ## Resizing
    dataT = skimage.transform.resize(
        dataR, (pixY, pixX), mode="constant", cval=0
    )  # xy reversed
    ## Padding
    dataP = np.pad(dataT, (padY, padX), "constant", constant_values=0)  # xy reversed

    return dataP


########################################################################


def data_transformation_MOS(data, header):
    """
    Performing geometrical transformations from raw coordinates to sky coordinates
    @param data: variability matrix
    @param header: header of the clean events file
    @return: transformed variability data
    """

    # Header information
    angle = header["PA_PNT"]

    xproj = [float(header["TDMIN6"]), float(header["TDMAX6"])]  # projected x limits
    yproj = [float(header["TDMIN7"]), float(header["TDMAX7"])]  # projected y limits
    xlims = [float(header["TLMIN6"]), float(header["TLMAX6"])]  # legal x limits
    ylims = [float(header["TLMIN7"]), float(header["TLMAX7"])]  # legal y limits

    # scaling factor
    sx = 648 / (xlims[1] - xlims[0])
    sy = 648 / (ylims[1] - ylims[0])

    # pads (padding)
    interX = (int((xproj[0] - xlims[0]) * sx), int((xlims[1] - xproj[1]) * sx))
    interY = (int((yproj[0] - ylims[0]) * sy), int((ylims[1] - yproj[1]) * sy))

    # adding pad according to MOS image pix (i.e: 500x500 pix)
    numX = int((148 - (interX[0] + interX[1])) / 2)
    numY = int((148 - (interY[0] + interY[1])) / 2)

    padX = (interX[0] + numX, 148 - (interX[0] + numX))
    padY = (interY[0] + numY, 148 - (interY[0] + numY))

    # Transformations
    ## Rotation
    dataR = np.flipud(nd.rotate(data, angle, reshape=False))
    ## Resizing (MOS image are based on 500x500 pix)
    dataT = skimage.transform.resize(dataR, (500, 500), mode="constant", cval=0)
    ## Padding
    dataP = np.pad(dataT, (padY, padX), "constant", constant_values=0)  # xy reversed

    return dataP
