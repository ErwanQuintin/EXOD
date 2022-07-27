#!/usr/bin/env python3
# coding=utf-8


########################################################################
#                                                                      #
# EXODUS - EPIC XMM-Newton Outburst Detector Ultimate System           #
#                                                                      #
#  Classify the SIMBAD matches into one of 5 classes                   #
#                                                                      #
# Maitrayee Gupta (2022) - maitrayee.gupta@irap.omp.eu                 #
#                                                                      #
########################################################################
"""
Given the SIMBAD object type, define a dictionary and classify 
the objects into one of the following:
NA, Unspecified, Galaxy, AGN, Star, CompactObject or X-Ray source
"""

from math import *
import csv
import pandas as pd
import math
import numpy as np

########################################################################
#                                                                      #
# Define the dictionary.                                               #
# Object is classified into one of the following:                      #
# NA, Unspecified, Galaxy, AGN, Star, CompactObject or X-Ray source    #
#                                                                      #
########################################################################

dic_classifier = {'': '',

                  'NA': 'NA',
                  'Unknown': 'Unspecified', 'IR': 'Unspecified', 'NearIR': 'Unspecified','Infrared': 'Unspecified','FarIR': 'Unspecified', 'Radio': 'Unspecified', 'MIR': 'Unspecified',
                  'NIR': 'Unspecified', 'HH': 'Unspecified', 'HI': 'Unspecified', 'HII': 'Unspecified',
                  'LensedImage': 'Unspecified', 'LensingEv': 'Unspecified', 'Maser': 'Unspecified',
                  'MolCld': 'Unspecified', 'PartofCloud': 'Unspecified', 'Radio(sub-mm)': 'Unspecified',
                  'Blue': 'Unspecified', 'Possible_lensImage': 'Unspecified', 'Unspecified': 'Unspecified',
                  'Radio(mm)': 'Unspecified', 'denseCore': 'Unspecified', 'Radio(cm)': 'Unspecified',
                  'UV': 'Unspecified', 'PN': 'Unspecified', 'PN?': 'Unspecified', "EmObj": 'Unspecified',
                  'DkNeb': 'Unspecified', 'Transient': 'Unspecified', 'Candidate_LensSystem': 'Unspecified',
                  'FIR': 'Unspecified', 'multiple_object': 'Unspecified', 'GravLensSystem': 'Unspecified',
                  'Bubble': 'Unspecified', 'Cloud': 'Unspecified', 'SFregion': 'Unspecified', 'Galaxy_Candidate': 'Unspecified',
                  'Inexistent': 'Unspecified', 'gamma': 'Unspecified', 'GravLens': 'Unspecified',
                  'HVCld': 'Unspecified', 'Candidate_Lens': 'Unspecified', 'ISM': 'Unspecified', 'Void': 'Unspecified',
                  'RfNeb': 'Unspecified', 'HIshell': 'Unspecified', 'Outflow': 'Unspecified', 'Region': 'Unspecified',   'cmRad': 'Unspecified',
                  'Globule': 'Unspecified', 'outflow?': 'Unspecified', 'ComGlob': 'Unspecified', 'ULX_Candidate': 'Unspecified', 'PlanetaryNeb_Candidate': 'Unspecified', 'smmRad': 'Unspecified', 'radioBurst': 'Unspecified', 'DarkNeb' : 'Unspecified',  'Planet_Candidate' : 'Unspecified',  'MidIR' : 'Unspecified', 'mmRad' : 'Unspecified', 'Blend' : 'Unspecified', 'GravLens_Candidate' : 'Unspecified', 'HIIReg': 'Unspecified', 

                  'GinCl': 'Galaxy', 'Galaxy': 'Galaxy', 'GiC': 'Galaxy', 'BClG': 'Galaxy', 'LSB_G': 'Galaxy',
                  'LensedG': 'Galaxy', 'GroupG': 'Galaxy', 'PartOfG': 'Galaxy', 'GinPair': 'Galaxy',
                  'Possible_ClG': 'Galaxy', 'Possible_G': 'Galaxy', 'Possible_GrG': 'Galaxy', 'GinGroup': 'Galaxy',
                  'HII_G': 'Galaxy', 'ClG': 'Galaxy', 'StarburstG': 'Galaxy', 'IG': 'Galaxy', 'SuperClG': 'Galaxy', 'GtowardsGroup' : 'Galaxy',
                  'PartofG': 'Galaxy', 'Compact_Gr_G': 'Galaxy', 'PairG': 'Galaxy', 'BlueCompG': 'Galaxy', 'GtowardsCl' : 'Galaxy', 'BrightestCG' : 'Galaxy', 'Cluster*_Candidate': 'Galaxy', 'GlobCluster':'Galaxy',  'EmissionG': 'Galaxy', 'GlobCluster_Candidate' : 'Galaxy', 'BlueCompactG': 'Galaxy', 'LowSurfBrghtG': 'Galaxy', 'ClG_Candidate': 'Galaxy', 'HIIG' : 'Galaxy',

                  'AGN': 'AGN', 'Sy1': 'AGN','Seyfert1': 'AGN', 'Sy2': 'AGN', 'Seyfert2': 'AGN', 'AGN_Candidate': 'AGN', 'QSO': 'AGN', 'Seyfert_1': 'AGN',
                  'Seyfert_2': 'AGN', 'LINER': 'AGN', 'EmG': 'AGN', 'RadioG': 'AGN', 'LensedQ': 'AGN', 'BLLac': 'AGN',
                  'Blazar': 'AGN', 'QSO_Candidate': 'AGN', 'Seyfert': 'AGN', 'Blazar_Candidate': 'AGN',
                  'BLLac_Candidate': 'AGN',

                  'OrionV*': 'Star', 'bCepV*': 'Star', 'Orion_V*': 'Star', 'TTau*': 'Star', 'EB*': 'Star', 'YSO': 'Star', 'SB*': 'Star', '**': 'Star',
                  'Star': 'Star', 'RotV*': 'Star', 'Candidate_RGB*': 'Star', 'low-mass*': 'Star', 'V*': 'Star',
                  'PulsV*': 'Star', 'AGB*': 'Star', 'S*': 'Star', 'Candidate_YSO': 'Star', 'PM*': 'Star',
                  'Irregular_V*': 'Star', 'Em*': 'Star', 'LPV*': 'Star', 'Mira': 'Star', 'WR*': 'Star', 'Pec*': 'Star',
                  'Planet?': 'Star', 'Planet': 'Star', 'Eruptive*': 'Star', 'Cl*': 'Star', 'OpCl': 'Star',
                  'Assoc*': 'Star', 'PulsV*WVir': 'Star', 'PulsV*bCep': 'Star', 'RRLyr': 'Star', 'C*': 'Star',
                  'EllipVar': 'Star', 'Candidate_EB*': 'Star', 'Candidate_PulsV*WVir': 'Star', 'Candidate_LP*': 'Star',
                  'pulsV*SX': 'Star', 'Candidate_RSG*': 'Star', 'BYDra': 'Star', 'Be*': 'Star',
                  'Candidate_RRLyr': 'Star', 'BlueSG*': 'Star', 'Erupt*RCrB': 'Star', 'RGB*': 'Star', 'RSCVn': 'Star',
                  'gammaDor': 'Star', 'Cl*?': 'Star', 'Candidate_C*': 'Star', 'HB*': 'Star', 'Cepheid': 'Star',
                  'Ae*': 'Star', 'Candidate_TTau*': 'Star', 'deltaCep': 'Star', 'HotSubdwarf': 'Star',
                  'Candidate_AGB*': 'Star', 'YellowSG*': 'Star', 'Symbiotic*': 'Star', 'PulsV*delSct': 'Star',
                  'BlueStraggler': 'Star', 'Candidate_post-AGB*': 'Star', 'RotV*alf2CVn': 'Star', 'OH/IR': 'Star',
                  'V*?': 'Star', 'Candidate_BSG*': 'Star', 'RedSG*': 'Star', 'Candidate_brownD*': 'Star',
                  'Candidate_Mi*': 'Star', 'Candidate_HB*': 'Star', 'Candidate_Be*': 'Star', 'Candidate_SN*': 'Star',
                  'brownD*': 'Star', 'SG*': 'Star', 'PulsV*RVTau': 'Star', 'Candidate_WR*': 'Star', 'HV*': 'Star',
                  'Candidate_Hsd': 'Star', 'Candidate_Ae*': 'Star', 'Candidate_Cepheid': 'Star', 'post-AGB*': 'Star',
                  'Candidate_**': 'Star', 'Candidate_Symb*': 'Star', 'Candidate_S*': 'Star', 'Candidate_SG*': 'Star',
                  'Candidate_low-mass*': 'Star','HighPM*' : 'Star','WolfRayet*' : 'Star' ,'Cluster*' : 'Star' , 'Variable*' : 'Star' ,
                  'GlCl': 'CompactObject', 'GlCl?': 'CompactObject', 'Pulsar': 'CompactObject', 'ULX': 'CompactObject', 'Low-Mass*' : 'Star' ,
                  'TTauri*': 'star', 'AGB*_Candidate' : 'star', 'YSO_Candidate': 'star', 'Supernova':'star', 'BYDraV*':'star', 'RSCVnV*':'star', 'alf2CVnV*':'star', 'RCrBV*':'star', 'RGB*_Candidate':'star', 'RRLyrae_Candidate':'star', 'RedSG' :'star', 'HorBranch*':'star', 'IrregularV*':'star', 'LongPeriodV*':'star', 'LongPeriodV*_Candidate':'star',  'BrownD*_Candidate':'star', 'CataclyV*_Candidate':'star', 'ChemPec*':'star', 'delSctV*':'star', 'gammaDorV*':'star', 'SNRemnant':'star', 'Supernova_Candidate':'star', 'Be*_Candidate':'star',  'ClassicalCep':'star', 'EmLine*':'star',  'Mira_Candidate':'star',  'OH/IR*':'star','BlueSG' :'star','PlanetaryNeb':'star', 'RefNeb' :'star', 'Association' :'star',
                  
                  
                  
                  
                  'ULX?': 'CompactObject', 'HMXB': 'CompactObject', 'Candidate_HMXB': 'CompactObject',
                  'LMXB': 'CompactObject', 'Candidate_LMXB': 'CompactObject', 'Nova': 'CompactObject',
                  'CataclyV*': 'CompactObject', 'XB': 'CompactObject', 'SNR': 'CompactObject', 'SNR?': 'CompactObject',
                  'Candidate_WD*': 'CompactObject', 'WD*': 'CompactObject', 'SN': 'CompactObject',
                  'gammaBurst': 'CompactObject', 'Candidate_XB*': 'CompactObject', 'Candidate_BH': 'CompactObject',
                  'NS': 'CompactObject', 'Candidate_NS': 'CompactObject', 'Neutron*': 'CompactObject',
                  'Candidate_CV*': 'CompactObject', 'Candidate_Nova': 'CompactObject', 'WhiteDwarf': 'CompactObject', 'WhiteDwarf_Candidate': 'CompactObject', 'XrayBin':'CompactObject', 'HighMassXBin':'CompactObject', 'LowMassXBin':'CompactObject', 'LowMassXBin_Candidate':'CompactObject', 'Neutron*_Candidate':'CompactObject','EclBin':'CompactObject', 'EclBin_Candidate':'CompactObject', 'blue':'CompactObject',

                   'X': 'X_ray_source',

                  }

simbad_col_no=29


# Read in existing table and append additional column with this class

with open('triple_match_obs_sim.csv','r') as csvinput:
    with open('triple_match_obs_sim_classified.csv', 'w') as csvoutput:
        writer = csv.writer(csvoutput, lineterminator='\n')
        reader = csv.reader(csvinput)

        all = []
        row = next(reader)
        row.append('SIMBAD_CLASS')
        all.append(row)

        for row in reader:
            simbad_class= row[simbad_col_no]
            simbad_group = '-'
            simbad_group = dic_classifier[simbad_class]
            row.append( str(simbad_group))
            all.append(row)

        writer.writerows(all)










