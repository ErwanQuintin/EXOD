import numpy as np
import subprocess

import file_names
from logger import logger

class Source:
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

    def __init__(self, id_src, inst, ccd, rawx, rawy, rawr):
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
        self.id_src = id_src
        self.inst = inst
        self.ccd = ccd
        self.rawx = rawx
        self.rawy = rawy 
        self.rawr = rawr
        self.vcount = None
        self.x = None
        self.y = None
        self.skyr = self.rawr * 64
        self.ra = None
        self.dec = None
        self.r = self.skyr * 0.05  # arcseconds
        self.var_rawx = rawx + 3  # 1.5 # 3
        self.var_rawy = rawy + 3  # 1.5 # 3
        self.var_rawr = rawr
        self.var_x = None
        self.var_y = None
        self.var_skyr = self.rawr * 64
        self.var_ra = None
        self.var_dec = None
        self.var_r = self.skyr * 0.05  # arcseconds


    def sky_coord(self, path, img):
        """
        Calculate sky coordinates with the sas task edet2sky.
        Return x, y, ra, dec
        """
        logger.debug(f"file util call rawx = {self.rawx} rawy = {self.rawy} ccd = {self.ccd}")

        # Launching SAS commands

        command = f"""
        export SAS_ODF={path};
        export SAS_CCF={path}ccf.cif;
        . {file_names.SAS};
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
        logger.debug(det2sky)

        # Writing the results in log file
        # Finding the beginning of the text to write in log
        deb = (np.where(det2sky == "Do not forget to define SAS_CCFPATH, SAS_CCF and SAS_ODF")[0][0]+ 2)

        # Equatorial coordinates
        self.ra, self.dec = det2sky[np.where(det2sky == "# RA (deg)   DEC (deg)")[0][0] + 1].split()

        logger.debug(f"RA = {self.ra} DEC = {self.dec}")

        # Sky pixel coordinates
        self.x, self.y = det2sky[np.where(det2sky == "# Sky X        Y pixel")[0][0] + 1].split()
        logger.debug(f'x = {self.x}, y = {self.y} rawx = {self.rawx}, rawy = {self.rawy}')

