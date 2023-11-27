"""
EXODUS - EPIC XMM-Newton Outburst Detector Ultimate System
Maitrayee Gupta (2022) - maitrayee.gupta@irap.omp.eu      

Various resources for both detector and renderer.
"""
import sys
import os
from math import *
import subprocess

import numpy as np
import scipy.ndimage as nd
import skimage.transform

from logger import logger
import file_names as FileNames



########################################################################
#                                                                      #
# Function and procedures                                              #
#                                                                      #
########################################################################


def open_files(folder_name):
    """
    Function opening files and writing their legend.
    @param  folder_name:  The directory to create the files.
    @return: log_file, info_file, variability_file, counter_per_tw,
    detected_variable_areas_file, time_windows_file
    """

    # Fixing the name of the folder
    if folder_name[-1] != "/":
        folder_name += "/"

    # Creating the folder if needed
    if not os.path.exists(folder_name):
        try:
            os.makedirs(folder_name)
        except:
            print("Error in creating output directory.\nABORTING", file=sys.stderr)
            exit(-1)

    # Declaring the files
    log_file = None
    var_file = None
    reg_file = None
    best_match_file = None

    # Creating the log file
    try:
        log_file = open(folder_name + FileNames.LOG, "w+")

    except IOError as e:
        print("Error in creating log.txt.\nABORTING", file=sys.stderr)
        print(e, file=sys.stderr)
        close_files(log_file, var_file, reg_file, best_match_file)
        print_help()
        exit(-1)

    # Creating the file to store variability per pixel
    try:
        var_file = folder_name + FileNames.VARIABILITY

    except IOError as e:
        print(
            "Error in creating {0}.\nABORTING".format(FileNames.VARIABILITY),
            file=sys.stderr,
        )
        print(e, file=sys.stderr)
        close_files(log_file, var_file, reg_file, best_match_file)
        print_help()
        exit(-1)

    # Creating the region file to store the position of variable sources
    try:
        reg_file = folder_name + FileNames.REGION

    except IOError as e:
        print(
            "Error in creating {0}.\nABORTING".format(FileNames.REGION), file=sys.stderr
        )
        print(e, file=sys.stderr)
        close_files(log_file, var_file, reg_file, best_match_file)
        print_help()
        exit(-1)

        # Creating the best match file to store the best matches from SIMBAD
    try:
        best_match_file = folder_name + FileNames.BEST_MATCH

    except IOError as e:
        print(
            "Error in creating {0}.\nABORTING".format(FileNames.BEST_MATCH),
            file=sys.stderr,
        )
        print(e, file=sys.stderr)
        close_files(log_file, var_file, reg_file, best_match_file)
        print_help()
        exit(-1)

    return log_file, var_file, reg_file, best_match_file


########################################################################


def close_files(log_f, var_f, reg_f, bes_match_f):
    """
    Function closing all files.
    """

    if log_f:
        log_f.close()

    # if var_f :
    #    var_f.close()

    # if reg_f :
    #    reg_f.close()


########################################################################


class Tee(object):
    """
    Class object that will print the output to the log_f file
    and to terminal
    """

    def __init__(self, *files):
        self.files = files

    def write(self, obj):
        for f in self.files:
            f.write(obj)
            f.flush()  # If you want the output to be visible immediately

    def flush(self):
        pass


########################################################################


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


class Source(object):
    """
    Datastructure providing easy storage for detected sources.

    Attributes:
    id_src:  The identifier number of the source
    inst:    The type of CCD
    ccd:     The CCD where the source was detected at
    rawx:    The x coordinate on the CCD
    rawy:    The y coordinate on the CCD
    r:       The radius of the variable area
    x:       The x coordinate on the output image
    y:       The y coordinate on the output image
    """

    def __init__(self, src):
        """
        Constructor for Source class. Computes the x and y attributes.
        @param src : source, output of variable_sources_position
            id_src:  The identifier number of the source
            inst:    The type of CCD
            ccd:     The CCD where the source was detected at
            rawx:    The x coordinate on the CCD
            rawy:    The y coordinate on the CCD
            r:       The raw pixel radius of the detected variable area
        """
        super(Source, self).__init__()

        self.id_src = src[0]
        self.inst = src[1]
        self.ccd = src[2]
        self.rawx = src[3]
        self.rawy = src[4]
        self.rawr = src[5]
        self.vcount = src[6]
        self.x = None
        self.y = None
        self.skyr = self.rawr * 64
        self.ra = None
        self.dec = None
        self.r = self.skyr * 0.05  # arcseconds
        self.var_rawx = src[3] + 3  # 1.5 # 3
        self.var_rawy = src[4] + 3  # 1.5 # 3
        self.var_rawr = src[5]
        self.var_x = None
        self.var_y = None
        self.var_skyr = self.rawr * 64
        self.var_ra = None
        self.var_dec = None
        self.var_r = self.skyr * 0.05  # arcseconds

    def sky_coord(self, path, img, log_f):
        """
        Calculate sky coordinates with the sas task edet2sky.
        Return x, y, ra, dec
        """
        print(
            "file util call x= ", self.rawx, " Y = ", self.rawy, " CCD num = ", self.ccd
        )
        # Launching SAS commands
        command = f"""
        export SAS_ODF={path};
        export SAS_CCF={path}ccf.cif;
        export HEADAS={FileNames.HEADAS};
        . $HEADAS/headas-init.sh;
        . {FileNames.SAS};
        echo "# Variable source {self.id_src}";
        edet2sky datastyle=user inputunit=raw X={self.rawx} Y={self.rawy} ccd={self.ccd} calinfoset={img} -V 0
        """

        # Running command
        process = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)

        # Extracting the output
        try:
            outs, errs = process.communicate(timeout=15)
        except TimeoutExpired:
            process.kill()
            outs, errs = process.communicate()

        # Converting output in utf-8
        textout = outs.decode("utf8")
        # Splitting for each line
        txt = textout.split("\n")
        # Converting in numpy array
        det2sky = np.array(txt)
        print(det2sky)

        # Writing the results in log file
        # Finding the beginning of the text to write in log
        deb = (
            np.where(
                det2sky == "Do not forget to define SAS_CCFPATH, SAS_CCF and SAS_ODF"
            )[0][0]
            + 2
        )
        # Writing in log file
        for line in txt[deb:]:
            log_f.write(line + "\n")

        # Equatorial coordinates
        self.ra, self.dec = det2sky[
            np.where(det2sky == "# RA (deg)   DEC (deg)")[0][0] + 1
        ].split()

        print("RA = ", self.ra, " DEC = ", self.dec)

        # Sky pixel coordinates
        self.x, self.y = det2sky[
            np.where(det2sky == "# Sky X        Y pixel")[0][0] + 1
        ].split()
        print("x = ", self.x, " Y= ", self.y)
        print("raw x = ", self.rawx, " raw Y= ", self.rawy)


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
