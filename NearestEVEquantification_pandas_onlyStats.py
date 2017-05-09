# -*- coding: utf-8 -*-
"""
@author: zwhitfield

Made to do the stats portion only for TEs near and/or overlapping EVEs
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Oct  3 11:04:24 2016

@author: zwhitfield
This script takes output of NearestEVEquantification_pandas_overlapOrNearest_createFiles.py and
NearestEVEquantification_GenomeWide_pandasBash_NowWithStats_FrozenDataNoTEfam.py

Tests significance of enrichment of certain TE types near EVEs compared to their prevelance in the entire genome

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

#This will be used to only select specific entries in dataset with TEs of the given type.
#Can be at any taxonomy level (class, subclass, family, element, etc...)
#The specific level is specified by filteredByCategory (ie name of column to subset using filteredBy)
filteredBy = str(sys.argv[3])
filteredByCategory = str(sys.argv[4])
groupedCategory = str(sys.argv[5])

analysisType = str(sys.argv[6])

# inputdir="/home/zwhitfield/Desktop/ForMarkGenomePaper/FrozenData/Aag2_assembly/publicReleaseOutputCheck/"
# outputdir="/home/zwhitfield/Desktop/ForMarkGenomePaper/FrozenData/Aag2_assembly/publicReleaseOutputCheck/forStats/"
# filteredBy="NONE" #Enter 'NONE' if don't want any filtering.
# filteredByCategory="TEclass"
# groupedCategory="TEclass"

possibleOrientations = ['same','opposite']
currentClassification = groupedCategory

#will not be graphing, but want to be consistent with other scripts
distanceCutoffs = [-20000,20000]
print "Quantifying upstream and downstream TEs closest to EVEs"

#---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------
#------------------------------FUNCTIONS------------------------------------------------------
#---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------

def getHistData (NearestTEdata, orientation, typeOfAnalysis):
    if orientation == 'same':
        upANDdownstreamMERGED = NearestTEdata[NearestTEdata['EVEstrand'] == NearestTEdata['TEstrand']]
    if orientation == 'opposite':
        upANDdownstreamMERGED = NearestTEdata[NearestTEdata['EVEstrand'] != NearestTEdata['TEstrand']]
    
    #upANDdownstreamMERGED will be by a particular taxonomy level of EVEs (i.e. LTRs)
    
    print ("Filtering TEs by " + filteredBy + " in the " + orientation + " orientation as EVEs" + " and grouping by " + currentClassification)
    
    #will become list of arrays of distances (one array per category of currentClassification)
    distances=list()
    
    #will become list of names of each category of currentClassification
    names = list()
    
    #Use pandas groupby command to group data frame by shared currentClassification (ie group by entries with same 'TEdescription')
    groupedDF = upANDdownstreamMERGED.groupby(currentClassification)

    #Obtain counts of each transposon in the entire genome, as produced by NearestEVEquantification_GenomeWide_pandas_NowWithStats.py
    wholeGenomeCounts = pd.read_csv(inputdir + '/ClassifiedBy_' + currentClassification+ '_FilteredBy_'+ filteredBy + '_WholeGenome.txt', sep='\t')

    #Get TOTAL number of TEs in genome
    wholeGenomeTotalCounts = sum(wholeGenomeCounts['Counts'])

    #Create outfile to hold the statistics of enrichment of TEs nearest EVEs versus whole genome. These comparisions have been filtered in the same way (i.e. filteredBy and currentClassification)
    if orientation == 'same':
        outfile= open(outputdir + 'ClassifiedBy_' + currentClassification+ '_FilteredBy_'+ filteredBy + '_StatsComparedToWholeGenome_SameStrand_' + typeOfAnalysis + '.txt','w')
    if orientation == 'opposite':
        outfile= open(outputdir + 'ClassifiedBy_' + currentClassification+ '_FilteredBy_'+ filteredBy + '_StatsComparedToWholeGenome_DifferentStrand_' + typeOfAnalysis + '.txt','w')

    outfile.write(currentClassification + '\t' + "specificTEcounts" + '\t' +"totalTEcountsOfType_" + filteredBy + '\t' +"pValueOneSidedBinom" + '\t' + "pValueOneSidedFEtest" + "\n")

    #Loop through groupedDF, getting the name (ie specific entry in currentClassification) and group (data in dataframe)
    #From there, extract distances from each group, append them to variable 'distance'. Do same for name
    #At the end, have 2 lists: each entry in 'distances' is all TE distances from EVEs of that currentCategory and each entry in 'names' is name of respective category in 'distances'
    #Example: If filteredBy was set to 'LTR', and classifications was set to 'TEfamily'
    #names=['Ty3_gypsy', 'Ty1_copia','Pao_Bel']
    #distances = [[all distances from EVEs of Ty3_gypsy elements], [all distances from EVEs of Ty1_copia elements], [all distances from EVEs of Pao Bel elements]]
    for name,group in groupedDF:
        distances.append(group['Distance'])
        names.append(name)
        #In order to perform enrichment tests, get number of TEs of given category ('names') in whole genome.
        specificTEcounts = wholeGenomeCounts[wholeGenomeCounts[currentClassification]==name].Counts.item()
        bTest = sp.binom_test(x=len(group['Distance']),n=len(upANDdownstreamMERGED['Distance']),p=float(specificTEcounts)/float(wholeGenomeTotalCounts),alternative='greater')
        oddsratio, FEpValue = sp.fisher_exact([[len(group['Distance']), specificTEcounts],
                                   [len(upANDdownstreamMERGED['Distance'])-len(group['Distance']), wholeGenomeTotalCounts-specificTEcounts]],
                                    alternative = 'greater')
        outfile.write(name + '\t' + str(len(group['Distance'])) + '\t' + str(len(upANDdownstreamMERGED['Distance'])) + '\t' + str(bTest) + '\t' + str(FEpValue) + "\n")
    outfile.close()

#---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------
#------------------------------END_OF_FUNCTIONS-----------------------------------------------
#---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
#------------------------------Read in EVE-TE pairs-----------------------------
#-------------------------------------------------------------------------------

filePath = inputdir + "TEsClosestToEVEs_" + analysisType + ".txt"
TEsNearestEVEs = pd.read_csv(filePath,
                             sep="\t")

#-------------------------------------------------------------------------------
#------------------------------Filter and calculate stats------------------
#-------------------------------------------------------------------------------

#Filter based on distance cutoffs (set at beginning of script)
TEsNearestEVEs = TEsNearestEVEs[(TEsNearestEVEs['Distance']>=distanceCutoffs[0]) & (TEsNearestEVEs['Distance']<=distanceCutoffs[1])]

#If desired, filter by desired class, established by filteredBy variable at top of script.
if filteredBy!= 'NONE':
    TEsNearestEVEs = TEsNearestEVEs[TEsNearestEVEs[filteredByCategory]==filteredBy]

#Perform enrichment tests
getHistData(TEsNearestEVEs, 'same', analysisType)
getHistData(TEsNearestEVEs, 'opposite', analysisType)