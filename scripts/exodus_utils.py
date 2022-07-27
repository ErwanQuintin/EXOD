#!/usr/bin/env python3
# coding=utf-8

########################################################################
#                                                                      #
# EXODUS - EPIC XMM-Newton Outburst Detector Ultimate System           #
#                                                                      #
# Rendering for exodus program                                         #
#                                                                      #
# Maitrayee Gupta (2022) - maitrayee.gupta@irap.omp.eu                 #
#                                                                      #
########################################################################
"""
Various functions for EXODUS program
"""

# Third-party imports

import numpy as np
from astropy.table import Table, Column
from astropy.coordinates import SkyCoord
from astropy import units as u

# Internal imports


########################################################################
#                                                                      #
# Check if the radoius of one source overlaps with the othe            #
#                                                                      #
########################################################################

def check_correlation(src_1, src_2, corr_tab) :
    """
	Function checking the correlation sources between two source lists
    
	@param  src_1:      The source list of first detector
	@param  src_2:      The source list of second detector
	@param  corr_table: The correlation table

	@return: The correlation table appended
	"""
    
    for i in range(len(src_1)):
        for j in range(len(src_2)):
            c1 = SkyCoord(src_1['RA'][i], src_1['DEC'][i], frame='fk5', unit='deg')
            c2 = SkyCoord(src_2['RA'][j], src_2['DEC'][j], frame='fk5', unit='deg')
            sep = c1.separation(c2)
            if sep.arcsecond < (src_1['R'][i]+src_2['R'][j]):
                corr_tab.add_row([src_1['ID'][i], src_1['INST'][i], src_1['RA'][i], src_1['DEC'][i], src_1['R'][i],\
                                src_2['ID'][j], src_2['INST'][j], src_2['RA'][j], src_2['DEC'][j], src_2['R'][j]])
    
    return corr_tab


########################################################################
#                                                                      #
#  Take 3 sources across the 3 EPIC cameras and check if               #
#   they are they same source.                                         #
#                                                                      #
########################################################################
    
def check_triple(corr_1, corr_2, corr_3) :
    """
	Function checking correlations between 3 correlation tables
    
	@param  corr_1: The first correlation table
	@param  corr_2: The second correlation table
	@param  corr_3: The third correlation table

	@return: A list with ID and INST of triple correlation sources
	"""
    
    src_cand = []
    triple = []

    for i in range(len(corr_1)):
        src_cand.append([corr_1['ID_1'][i],corr_1['INST_1'][i],corr_1['ID_2'][i],corr_1['INST_2'][i]])
        
    for cand in src_cand:
        t1 = corr_3[np.where((corr_3['ID_1']==cand[0]) & (corr_3['INST_1']==cand[1]))]
        t2 = corr_2[np.where((corr_2['ID_1']==cand[2]) & (corr_2['INST_1']==cand[3]))]
        
    for line1 in t1:
        for line2 in t2:
            if line1['ID_2']==line2['ID_2'] and line1['INST_2']==line2['INST_2']:
                triple.append([line1['ID_1'],line1['INST_1'],line2['ID_1'],line2['INST_1'],line1['ID_2'],line1['INST_2']])
    
    return triple

########################################################################
#                                                                      #
# Check if one source falls within the radius of the other             #
#                                                                      #
########################################################################
    
def check_multiple_sources(src) :
    """
	Function checking if a source is already included in another one
	@param  src:   The astropy source table
	@return: The astropy source table without multiple sources
	"""
    print("here at check_multiple_sources")
    src.sort(['R'])     # table sorted in radius, from the largest to the shortest
    src.reverse()       # in order to eliminate sources with shorter radius

    multiple = []
    for i in range(len(src)):
        for j in range(i+1,len(src)):
            c1 = SkyCoord(src['RA'][i], src['DEC'][i], frame='fk5', unit='deg')
            c2 = SkyCoord(src['RA'][j], src['DEC'][j], frame='fk5', unit='deg')
            sep = c1.separation(c2)
            #if sep.arcsecond < (src_M2['R'][i]-(src_M2['R'][j]/2)):
            if sep.arcsecond < (src['R'][i]):
                multiple.append(j)
                
    src.remove_rows(multiple)   # removing multiple sources
    src.sort(['ID'])
    for i in range(len(src)):
        src['ID'][i]=i+1        # Putting sources in the same order
    
    return src

########################################################################
#                                                                      #
# Cross corelate sources between different detectors.                  #
#                                                                      #
########################################################################
   
def correl_flag(src, corr_table, triple_l) :
    """
	Function creating a column with correlation flag
        @param  src: The soure list from FITS record
	@param  corr_table: The correlation table
	@return: An astropy Table with correl flag column added
	"""
    # Changing in astropy table
    src1 = Table(src)
    # Adding an empty column
    A=np.empty((len(src1),), dtype='S25')
    col_c = Column(name='correl', data=A)
    src1.add_column(col_c)
    # Initializing empty
    src1['correl']=''
    
    # Sorting the table
    corr_1 = corr_table[np.where((corr_table['INST_1']=='PN') & (corr_table['INST_2']=='M1'))]
    corr_2 = corr_table[np.where((corr_table['INST_1']=='M1') & (corr_table['INST_2']=='M2'))]
    corr_3 = corr_table[np.where((corr_table['INST_1']=='PN') & (corr_table['INST_2']=='M2'))]

    # Flagging the correlation column for 3 EPIC
    if src1['INST'][0]=='PN':
        for s in src1:
            if s['ID'] in corr_1['ID_1']:
                s['correl']+='M1 '
            elif s['ID'] in corr_3['ID_1']:
                s['correl']+='M2 '
            for j in range(len(triple_l)):
                if s['ID']==triple_l[j][0]:
                    s['correl']='Triple'
    
    if src1['INST'][0]=='M1':
        for s in src1:
            if s['ID'] in corr_1['ID_2']:
                s['correl']+='PN '
            elif s['ID'] in corr_2['ID_1']:
                s['correl']+='M2 '
            for j in range(len(triple_l)):
                if s['ID']==triple_l[j][2]:
                    s['correl']='Triple'
                    
    if src1['INST'][0]=='M2':
        for s in src1:
            if s['ID'] in corr_2['ID_2']:
                s['correl']+='M1 '
            elif s['ID'] in corr_3['ID_2']:
                s['correl']+='PN '
            for j in range(len(triple_l)):
                if s['ID']==triple_l[j][4]:
                    s['correl']='Triple'
    
    return src1

########################################################################
