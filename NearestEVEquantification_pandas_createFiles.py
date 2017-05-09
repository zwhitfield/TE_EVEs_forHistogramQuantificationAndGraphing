# -*- coding: utf-8 -*-
"""
@author: zwhitfield

Made to do the stats portion only for TEs near and/or overlapping EVEs
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Oct  3 11:04:24 2016

@author: zwhitfield
This script takes in a table of EVEs with nearest upstream neighbor and nearest downstream neighbor (in two different files)
Outputs file ("TEsClosestToEVEs_overlapOrNearest.txt"), which has all TEs which overlap an EVE. If an EVE has no
overlapping TE, it will give the closest upstream AND downstream TE for that EVE.
"""

import sys
import pandas as pd
import numpy as np
import scipy.stats as sp
#scipy.__version__
#np.__version__
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
plt.ion()
#matplotlib.__version__
plt.style.use('ggplot')

inputdir = str(sys.argv[1])
outputdir = str(sys.argv[2])
analysisType = str(sys.argv[3])

# inputdir = "/home/zwhitfield/Desktop/ForMarkGenomePaper/FrozenData/Aag2_assembly/"
# outputdir="/home/zwhitfield/Desktop/ForMarkGenomePaper/FrozenData/Aag2_assembly/publicReleaseOutputCheck/"

#---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------
#------------------------------FUNCTIONS------------------------------------------------------
#---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------


def loadAndOrganizeData(TEfilepath,location,overlapStatus):
    TEfile =  pd.read_csv(TEfilepath,
                          names = ["ContigEVE","EVEstart","EVEend","EVEdescription","EVEscore","EVEstrand","EVEsomething","ContigTE","TEstart","TEend","TEdescription","TEscore","TEstrand","TEfamily","Distance"],
                          sep = "\t")
    
    #First, rename all LTR/Gypsy as LTR/Ty3_gyspsy. THere are both in the file, so need a consensus
    #First, rename all LTR/Copia as LTR/Ty1_copia. THere are both in the file, so need a consensus
    #First, rename all LTR/Pao as LTR/Pao_Bel. THere are both in the file, so need a consensus
    TEfile.loc[TEfile.TEfamily == 'LTR/Gypsy', 'TEfamily'] = "LTR/Ty3_gypsy"
    TEfile.loc[TEfile.TEfamily == 'LTR/Copia', 'TEfamily'] = "LTR/Ty1_copia"
    TEfile.loc[TEfile.TEfamily == 'LTR/Pao', 'TEfamily'] = "LTR/Pao_Bel"
                      
    #Split last column into two (there must be a better way to do this)
    TEfile["TEclass"]=TEfile["TEfamily"].str.split('/').str[0]
    TEfile["TEfamily"]=TEfile["TEfamily"].str.split('/').str[1]

    TEfile = TEfile[(TEfile['TEstart']>0)]

    TEfile["OverlapCategory"] = "NotApplicable"

    # Filter by EVEstart relative to TEstart.
    # With Overlaps allowed,distance of 0 is counted for upstream and downstream,
    # regardless of where other end of TE is relative to EVE, so need to manually categorize.
    # upstream is if TE starts before EVE, and end overlaps with EVE sequence. Assigned OverlapCategory of -1.9
    # downstream is if TE starts overlaps with EVE sequence, and ends after EVE sequence. Assigned OverlapCategory of 1.9
    # surround is if TE starts before, and ends after, EVE sequence. Assigned OverlapCategory of 0.9
    # contained is if TE start and stops within border of EVE. Assigned OverlapCategory of -0.9
    if overlapStatus == 'yes':
        TEfile = TEfile[(TEfile['Distance'] == 0)]

        if location == 'upstream':
            TEfile = TEfile[(TEfile['EVEstart']>TEfile['TEstart']) & (TEfile['EVEend']>TEfile['TEend'])]
            TEfile["OverlapCategory"] = -1.9

        elif location == 'downstream':
            TEfile = TEfile[(TEfile['EVEstart']<TEfile['TEstart']) & (TEfile['EVEend']<TEfile['TEend'])]
            TEfile["OverlapCategory"] = 1.9

        elif location == 'surround':
            TEfile = TEfile[(TEfile['EVEstart']>TEfile['TEstart']) & (TEfile['EVEend']<TEfile['TEend'])]
            TEfile["OverlapCategory"] = 0.9
        elif location == 'contained':
            TEfile = TEfile[(TEfile['EVEstart']<TEfile['TEstart']) & (TEfile['EVEend']>TEfile['TEend'])]
            TEfile["OverlapCategory"] = -0.9

    if overlapStatus == 'no':
        TEfile = TEfile[(TEfile['Distance'] != 0)]

        if location == 'upstream':
            TEfile = TEfile[(TEfile['Distance'] < 0)]

        elif location == 'downstream':
            TEfile = TEfile[(TEfile['Distance'] > 0)]

    return TEfile


def loadAndOrganizeData_withTaxonomy(TEfilepath, location, overlapStatus):
    TEfile = pd.read_csv(TEfilepath,
                         sep="\t")

    # First, rename all LTR/Gypsy as LTR/Ty3_gyspsy. THere are both in the file, so need a consensus
    # First, rename all LTR/Copia as LTR/Ty1_copia. THere are both in the file, so need a consensus
    # First, rename all LTR/Pao as LTR/Pao_Bel. THere are both in the file, so need a consensus
    TEfile.loc[TEfile.TEfamily == 'LTR/Gypsy', 'TEfamily'] = "LTR/Ty3_gypsy"
    TEfile.loc[TEfile.TEfamily == 'LTR/Copia', 'TEfamily'] = "LTR/Ty1_copia"
    TEfile.loc[TEfile.TEfamily == 'LTR/Pao', 'TEfamily'] = "LTR/Pao_Bel"

    # Split last column into two (there must be a better way to do this)
    TEfile["TEclass"] = TEfile["TEfamily"].str.split('/').str[0]
    TEfile["TEfamily"] = TEfile["TEfamily"].str.split('/').str[1]

    TEfile = TEfile[(TEfile['TEstart'] > 0)]

    TEfile["OverlapCategory"] = "NotApplicable"
    # Filter by EVEstart relative to TEstart. With Overlaps allowed,distance of
    # With Overlaps allowed,distance of 0 is counted for upstream and downstream,
    # regardless of where other end of TE is relative to EVE, so need to manually categorize.
    # upstream is if TE starts before EVE, and end overlaps with EVE sequence. Assigned OverlapCategory of -1.9
    # downstream is if TE starts overlaps with EVE sequence, and ends after EVE sequence. Assigned OverlapCategory of 1.9
    # surround is if TE starts before, and ends after, EVE sequence. Assigned OverlapCategory of 0.9
    # contained is if TE start and stops within border of EVE. Assigned OverlapCategory of -0.9
    if overlapStatus == 'yes':
        TEfile = TEfile[(TEfile['Distance'] == 0)]

        if location == 'upstream':
            TEfile = TEfile[(TEfile['EVEstart'] > TEfile['TEstart']) & (TEfile['EVEend'] > TEfile['TEend'])]
            TEfile["OverlapCategory"] = -1.9

        elif location == 'downstream':
            TEfile = TEfile[(TEfile['EVEstart'] < TEfile['TEstart']) & (TEfile['EVEend'] < TEfile['TEend'])]
            TEfile["OverlapCategory"] = 1.9

        elif location == 'surround':
            TEfile = TEfile[(TEfile['EVEstart'] > TEfile['TEstart']) & (TEfile['EVEend'] < TEfile['TEend'])]
            TEfile["OverlapCategory"] = 0.9
        elif location == 'contained':
            TEfile = TEfile[(TEfile['EVEstart'] < TEfile['TEstart']) & (TEfile['EVEend'] > TEfile['TEend'])]
            TEfile["OverlapCategory"] = -0.9

    if overlapStatus == 'no':
        TEfile = TEfile[(TEfile['Distance'] != 0)]

        if location == 'upstream':
            TEfile = TEfile[(TEfile['Distance'] < 0)]

        elif location == 'downstream':
            TEfile = TEfile[(TEfile['Distance'] > 0)]

    return TEfile

def concatAndsave (fileToSave, typeOfAnalysis, outputdirectory, taxStatus):
    # Turn into one big dataframe, resetting the index
    # This ends up treating the up and downstream TEs independently. It all becomes one big dataset.
    upANDdownstreamMERGED_Master = pd.concat(fileToSave, ignore_index=True)

    if taxStatus == "noTaxonomy":
        upANDdownstreamMERGED_Master.to_csv(outputdirectory + "TEsClosestToEVEs_" +  typeOfAnalysis + ".txt",
                                            sep="\t", index=False, quoting=False)
    if taxStatus == "withTaxonomy":
        upANDdownstreamMERGED_Master.to_csv(outputdirectory + "TEsClosestToEVEs_" +  typeOfAnalysis + "_withEVEtaxonomy.txt",
                                            sep="\t", index=False, quoting=False)



#---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------
#------------------------------END_OF_FUNCTIONS-----------------------------------------------
#---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------


#---------------------------------------------------------------------------------------------
#------------------------------Generate _overlapOrNearest files-------------------------------
#---------------------------------------------------------------------------------------------

if analysisType == "overlapOrNearest":
    # Get actaul categories of nearest TEs to EVEs. Overlapping TEs are same for upstream and downstream, so just pick one
    # to use as input file
    upstreamOVERLAP = loadAndOrganizeData(inputdir + "closestTEtoEVEs_UPSTREAM_allowOverlap.txt", "upstream", "yes")
    downstreamOVERLAP = loadAndOrganizeData(inputdir + "closestTEtoEVEs_DOWNSTREAM_allowOverlap.txt", "downstream", "yes")
    surroundOVERLAP = loadAndOrganizeData(inputdir + "closestTEtoEVEs_DOWNSTREAM_allowOverlap.txt", "surround", "yes")
    containedOVERLAP = loadAndOrganizeData(inputdir + "closestTEtoEVEs_DOWNSTREAM_allowOverlap.txt", "contained", "yes")

    upstreamNEARESTONLY = loadAndOrganizeData(inputdir + "closestTEtoEVEs_UPSTREAM_allowOverlap.txt", "upstream", "no")
    downstreamNEARESTONLY = loadAndOrganizeData(inputdir + "closestTEtoEVEs_DOWNSTREAM_allowOverlap.txt", "downstream", "no")

    # Combine and filter datasets
    # Stack/combine the dataframes
    # By combining OVERLAP and NEARESTONLY, end up with all EVEs with overlapping TEs (ie distance of 0), or if no overlap,
    # the nearest upstream AND downstream TE.
    dframes = [upstreamOVERLAP, downstreamOVERLAP, surroundOVERLAP, containedOVERLAP,
               upstreamNEARESTONLY, downstreamNEARESTONLY]

    # Do not include the "contained" category
    # dframes = [upstreamOVERLAP, downstreamOVERLAP, surroundOVERLAP,
    #            upstreamNEARESTONLY,downstreamNEARESTONLY]

    concatAndsave(dframes, "overlapOrNearest", outputdir, "noTaxonomy")

    #Now, same as above but for files containing EVE taxonomy assignments as well.

    upstreamOVERLAP = loadAndOrganizeData_withTaxonomy(inputdir + "closestTEtoEVEs_UPSTREAM_allowOverlap.txt_withTaxonomy.txt", "upstream", "yes")
    downstreamOVERLAP = loadAndOrganizeData_withTaxonomy(inputdir + "closestTEtoEVEs_DOWNSTREAM_allowOverlap.txt_withTaxonomy.txt", "downstream", "yes")
    surroundOVERLAP = loadAndOrganizeData_withTaxonomy(inputdir + "closestTEtoEVEs_DOWNSTREAM_allowOverlap.txt_withTaxonomy.txt", "surround", "yes")
    containedOVERLAP = loadAndOrganizeData_withTaxonomy(inputdir + "closestTEtoEVEs_DOWNSTREAM_allowOverlap.txt_withTaxonomy.txt", "contained", "yes")

    upstreamNEARESTONLY = loadAndOrganizeData_withTaxonomy(inputdir + "closestTEtoEVEs_UPSTREAM_allowOverlap.txt_withTaxonomy.txt", "upstream", "no")
    downstreamNEARESTONLY = loadAndOrganizeData_withTaxonomy(inputdir + "closestTEtoEVEs_DOWNSTREAM_allowOverlap.txt_withTaxonomy.txt", "downstream", "no")

    #Combine and filter datasets
    #Stack/combine the dataframes
    #By combining OVERLAP and NEARESTONLY, end up with all EVEs with overlapping TEs (ie distance of 0), or if no overlap,
    # the nearest upstream AND downstream TE.
    dframes = [upstreamOVERLAP, downstreamOVERLAP, surroundOVERLAP, containedOVERLAP,
               upstreamNEARESTONLY,downstreamNEARESTONLY]

    #Do not include the "contained" category
    # dframes = [upstreamOVERLAP, downstreamOVERLAP, surroundOVERLAP,
    #            upstreamNEARESTONLY,downstreamNEARESTONLY]

    concatAndsave(dframes, "overlapOrNearest",outputdir, "withTaxonomy")


#---------------------------------------------------------------------------------------------
#------------------------------Generate nearest TE only files---------------------------------
#---------------------------------------------------------------------------------------------

if analysisType == "nearestOnly":
    #Read in files that did NOT allow TEs which overlapped EVEs
    upstreamNEARESTONLY = loadAndOrganizeData(inputdir + "closestTEtoEVEs_UPSTREAM.txt", "upstream", "no")
    downstreamNEARESTONLY = loadAndOrganizeData(inputdir + "closestTEtoEVEs_DOWNSTREAM.txt", "downstream", "no")

    # Create vector of dataframes
    dframes = [upstreamNEARESTONLY, downstreamNEARESTONLY]

    concatAndsave(dframes, "nearestOnly", outputdir, "noTaxonomy")

    #Now, same as above but for files containing EVE taxonomy assignments as well.

    upstreamNEARESTONLY = loadAndOrganizeData_withTaxonomy(inputdir + "closestTEtoEVEs_UPSTREAM.txt_withTaxonomy.txt", "upstream", "no")
    downstreamNEARESTONLY = loadAndOrganizeData_withTaxonomy(inputdir + "closestTEtoEVEs_DOWNSTREAM.txt_withTaxonomy.txt", "downstream", "no")

    #Create vector of dataframes
    dframes = [upstreamNEARESTONLY,downstreamNEARESTONLY]

    concatAndsave(dframes, "nearestOnly",outputdir, "withTaxonomy")


#---------------------------------------------------------------------------------------------
#------------------------------Generate overlap TE only files---------------------------------
#---------------------------------------------------------------------------------------------
if analysisType == "overlapOnly":
    upstreamOVERLAP = loadAndOrganizeData(inputdir + "closestTEtoEVEs_UPSTREAM_allowOverlap.txt", "upstream", "yes")
    downstreamOVERLAP = loadAndOrganizeData(inputdir + "closestTEtoEVEs_DOWNSTREAM_allowOverlap.txt", "downstream", "yes")
    surroundOVERLAP = loadAndOrganizeData(inputdir + "closestTEtoEVEs_DOWNSTREAM_allowOverlap.txt", "surround", "yes")
    containedOVERLAP = loadAndOrganizeData(inputdir + "closestTEtoEVEs_DOWNSTREAM_allowOverlap.txt", "contained", "yes")

    dframes = [upstreamOVERLAP, downstreamOVERLAP, surroundOVERLAP, containedOVERLAP]

    # Do not include the "contained" category
    # dframes = [upstreamOVERLAP, downstreamOVERLAP, surroundOVERLAP,
    #            upstreamNEARESTONLY,downstreamNEARESTONLY]

    concatAndsave(dframes, "overlapOnly", outputdir, "noTaxonomy")

    #Now, same as above but for files containing EVE taxonomy assignments as well.

    upstreamOVERLAP = loadAndOrganizeData_withTaxonomy(inputdir + "closestTEtoEVEs_UPSTREAM_allowOverlap.txt_withTaxonomy.txt", "upstream", "yes")
    downstreamOVERLAP = loadAndOrganizeData_withTaxonomy(inputdir + "closestTEtoEVEs_DOWNSTREAM_allowOverlap.txt_withTaxonomy.txt", "downstream", "yes")
    surroundOVERLAP = loadAndOrganizeData_withTaxonomy(inputdir + "closestTEtoEVEs_DOWNSTREAM_allowOverlap.txt_withTaxonomy.txt", "surround", "yes")
    containedOVERLAP = loadAndOrganizeData_withTaxonomy(inputdir + "closestTEtoEVEs_DOWNSTREAM_allowOverlap.txt_withTaxonomy.txt", "contained", "yes")

    dframes = [upstreamOVERLAP, downstreamOVERLAP, surroundOVERLAP, containedOVERLAP]

    #Do not include the "contained" category
    # dframes = [upstreamOVERLAP, downstreamOVERLAP, surroundOVERLAP,
    #            upstreamNEARESTONLY,downstreamNEARESTONLY]

    concatAndsave(dframes, "overlapOnly",outputdir, "withTaxonomy")
