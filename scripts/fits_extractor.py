"""
EXODUS - EPIC XMM-Newton Outburst Detector Ultimate System
Maitrayee Gupta (2022) - maitrayee.gupta@irap.omp.eu

Variability-related procedures specified into the documentation
"""
import time
from os.path import sys

import numpy as np
from astropy import wcs
from astropy.io import fits
from astropy.table import Table

from logger import logger

########################################################################


def extraction_photons(events_file):
    """
    Function extracting the E round list from its FITS events file.
    It alse returns the header information
    @param events_file: The events FITS file
    @return: The E round list
    @return: The events file header
    @return maximum, minimum
    @raise Exception: An exception from astropy if something went wrong
    """
    logger.debug('calling extraction_photons()')
    logger.debug(f'events_file={events_file}')

    hdulist = fits.open(events_file)
    events = hdulist[1].data
    header = hdulist[1].header

    # We will make use of the fact that the events file is already sorted by CCDs
    # First we split the events file by CCDs by retrieving the indices where the CCD# changes
    indices_ccd_change = np.where(np.diff(events["CCDNR"]))[0]
    events_filtered = np.split(np.array(events), indices_ccd_change)
    events_filtered = [np.array(ccd_events) for ccd_events in events_filtered]

    # Then we sort the events of each CCDs in chronological order
    events_filtered_sorted = [ccd_events[np.argsort(ccd_events["TIME"])] for ccd_events in events_filtered]
    return events_filtered_sorted, header


########################################################################


# Deprecated
def extraction_info(events_file):
    hdulist = fits.open(events_file)
    header = hdulist[1].header
    data = fits.getdata(events_file)
    dmin = [min(data["X"]), min(data["Y"])]
    dmax = [max(data["X"]), max(data["Y"])]
    return header, dmin, dmax


########################################################################


def get_gti_from_file(gti_file):
    """
    Function extracting the G round list from its Fits file.
    @param gti_file: The gti file
    @return: The G round list
    @raise Exception: An exception from astropy if something went wrong
    """
    logger.debug(f'Getting GTIs from file={gti_file}')
    hdulist = fits.open(gti_file)
    return hdulist[1].data


########################################################################


def fits_writer(data, sources, image, pars, file):
    """
    Function writing the variability and sources to a fits file
    @param data: image of the variability data
    @param sources: detected variable sources
    @param image: image obtained with evselect, needed for the header
    @param pars: variability parameters used in the variability computation
    @param file: output file name
    """

    hdulist = fits.open(image)
    header_img = hdulist[0].header
    hdulist.close()
    head_var_f = header_img
    head_var_f.append(card=("CREATOR", pars["CREATOR"], "[s] EXOD Time window"))
    head_var_f.append(card=("DATE", pars["DATE"], "[s] EXOD Time window"))
    head_var_f.append(card=("TW", pars["TW"], "[s] EXOD Time window"))
    head_var_f.append(card=("GTR", pars["GTR"], "EXOD Good time ratio"))
    head_var_f.append(card=("DL", pars["DL"], "EXOD Detection level"))
    head_var_f.append(card=("BS", pars["BS"], "[pix] EXOD Box size"))

    # Creating fits file
    hdul_f = fits.HDUList()
    hdul_var = fits.ImageHDU(data=data, header=head_var_f)
    hdul_src = fits.BinTableHDU(data=sources)
    hdul_f.append(hdul_var)
    hdul_f.append(hdul_src)

    # Writing to file
    hdul_f.writeto(file, overwrite=True)

    return True
