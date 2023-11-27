# EXODUS
## EPIC XMM-Newton Outburst Detector Ultimate System

The purpose of EXODUS is to detect faint rapid transients in XMM-Newton data that may not have been detected automatically using the XMM-Newton pipeline. It works for all EPIC detectors and all read-out modes. The technique used in EXODUS is to break a long observation into short time windows and compare the number of counts in regions of the detector for each time window, to the average number of counts per region. EXODUS uses three sub-routines: the first one calculates the counts per region per window and determines the variability by comparing to the average rate for that region; the second one provides an image showing the level of variability on the whole detector, and a third one performs post-processing and collates the results.

We encourage potential users to read the user's guide (`EXODUS_users_guide.pdf`), and especially follow the tutorial presented in section 6.

This document provides both user and technical documentation for EXODUS. If you use EXODUS for your research, please acknowledge it by citing Gupta et al (to be submitted). EXODUS is built upon two projects: the original project 

"[Variabilitectron](https://framagit.org/DWojtowicz/Variabilitectron)" created by Damien Wojtowicz, and "[EXOD](https://www.aanda.org/articles/aa/full_html/2020/08/aa36869-19/aa36869-19.html)" created by [Pastor-Marazuela et al (2020)](https://www.aanda.org/articles/aa/full_html/2020/08/aa36869-19/aa36869-19.html).

Version History:
- [2017 - Variabilitectron - Damien Wojtowicz](https://framagit.org/DWojtowicz/Variabilitectron)
- [2018 - Variabilitectron2 - Inés Pastor Marazuela](https://framagit.org/InesPM/Variabilitectron)
- [2019 - EXOD - Inés Pastor Marazuela](https://github.com/InesPM/EXOD)
- [2020 - EXOD2 - Maitrayee Gupta & Florent Castellani](https://github.com/Monrillo/EXOD)
- [2022 - EXODUS - Maitrayee Gupta & Erwan Quintin](https://github.com/ErwanQuintin/EXOD)
- [2023 - ? - Erwan Quintin & Norman Khan](https://github.com/ErwanQuintin/EXOD)

## Tutorial

Let's set the path to where the scripts are located and where we want to store our data, for instance

```
SCRIPTS=/path/EXOD/scripts
FOLDER=/path/data
```

One can use the `run_exodus.sh` script to download, filter and compute the variability for a list of observations. Let's try, for instance, with observation 0831790701, and (see the obs_ids.txt in the examples folder):

```
bash run_exodus.sh -o obs_ids.txt -l true -i true
```
An example of the output of these commands can be found in the folder `examples`.

Belowe we see the results of run for OBSID 0831790701.


![variability](../master/example/merged_sources.png)
![variability](../master/example/merged_lc.png)
