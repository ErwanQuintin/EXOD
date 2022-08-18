# EXODUS
## EPIC XMM-Newton Outburst Detector Ultimate System

The purpose of EXODUS is to detect faint rapid transients in XMM-Newton data, that may not have been detected automatically using the XMM-Newton pipeline. It works for all EPIC detectors and all read-out modes. The technique used in EXODUS is to break a long observation into short time windows and compare the number of counts in regions of the detector for each time window, to the average number of counts per region. EXODUS uses three sub-routines: the first one calculates the counts per region per window and determines the variability by comparing to the average rate for that region; the second one provides an image showing the level of variability on the whole detector; and a third, performs post-processing and collates the results.

We encourage the potential users to read the users guide (EXODUS_users_guide.pdf), and especially follow the tutorial presented in section 6.

This document provides both user and technical documentation for EXODUS.
If you use EXODUS for your research, please acknowledge it by citing Gupta et al (to be submitted).
EXODUS is built upon two projects: the oOriginal project "<a href="https://framagit.org/DWojtowicz/Variabilitectron" target="_blank">Variabilitectron</a>", which was created by Damien Wojtowicz, and "EXOD" which was created by  <a href="https://www.aanda.org/articles/aa/full_html/2020/08/aa36869-19/aa36869-19.html" target="_blank">Pastor-Marazuela et al (2020)</a>.


See previous versions:
<a href="https://framagit.org/DWojtowicz/Variabilitectron" target="_blank">Variabilitectron v1</a>.
<a href="https://framagit.org/InesPM/Variabilitectron" target="_blank">Variabilitectron v2</a>. 
<a href="https://github.com/InesPM/EXOD" target="_blank">EXOD v1</a>.
<a href="https://github.com/Monrillo/EXOD" target="_blank">EXOD v2</a>.

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
