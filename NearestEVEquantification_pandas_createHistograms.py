# -*- coding: utf-8 -*-
"""
Created on Mon Nov 14 09:57:39 2016

@author: zwhitfield
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Oct  3 11:04:24 2016

@author: zwhitfield
This script takes output of NearestEVEquantification_pandas_overlapOrNearest_createFiles.py

Plots histogram depicting makeup of types of TEs nearest EVEs.

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
filteredBy = str(sys.argv[3]) #This will be used to only select/analyze specific entries in dataset with TEs of the given type.
filteredByCategory = str(sys.argv[4])#The specific level is specified by filteredByCategory (ie name of column to subset using filteredBy)
groupedCategory = str(sys.argv[5])#With TEs of only type 'filteredBy', plot at what taxonomy level?
analysisType = str(sys.argv[6])

# inputdir="/home/zwhitfield/Desktop/ForMarkGenomePaper/FrozenData/Aag2_assembly/"
# outputdir="/home/zwhitfield/Desktop/ForMarkGenomePaper/FrozenData/Aag2_assembly/Figures/Histograms/"
# filteredBy="LTR" #Enter 'NONE' if don't want any filtering.
# filteredByCategory="TEclass"
# groupedCategory="TEfamily"
# analysisType = "nearestOnly"

possibleOrientations = ['same','opposite']
currentClassification = groupedCategory

print "Quantifying upstream and downstream TEs closest to EVEs"

#---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------
#------------------------------FUNCTIONS------------------------------------------------------
#---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------



def getHistData (NearestTEdata, orientation, overlapStatus):
    if orientation == 'same':
        upANDdownstreamMERGED = NearestTEdata[NearestTEdata['EVEstrand'] == NearestTEdata['TEstrand']]
    if orientation == 'opposite':
        upANDdownstreamMERGED = NearestTEdata[NearestTEdata['EVEstrand'] != NearestTEdata['TEstrand']]

    if overlapStatus == 'yes':
        numBins = 4
        upANDdownstreamMERGED["OverlapCategory"] = pd.to_numeric(upANDdownstreamMERGED["OverlapCategory"])
    if overlapStatus == 'no':
        numBins = 20

    print ("Graphing TEs by " + filteredBy + " in the " + orientation + " orientation as EVEs" + " and grouping by " + currentClassification)
    
    #will become list of arrays of distances (one per category of currentClassification)
    distances=list()
    
    #will become list of names of each category of currentClassification
    names = list()
    
    #Use pandas groupby command to group data frame by shared currentClassification (ie group by entries with same 'TEdescription')
    groupedDF = upANDdownstreamMERGED.groupby(currentClassification)

    #Loop through groupedDF, getting the name (ie specific entry in currentClassification) and group (data in dataframe)
    #From there, extract distances from each group, append them to variable 'distance'. Do same for name
    #At the end, have 2 lists: each entry in 'distances' is all distances from EVEs of that currentCategory and each entry in 'names' is name of respective category in 'distances'
    #Example: If filteredBy was set to 'LTR', and classifications was set to 'TEfamily'
    #names=['Ty3_gypsy', 'Ty1_copia','Pao_Bel']
    #distances = [[all distances from EVEs of Ty3_gypsy elements], [all distances from EVEs of Ty1_copia elements], [all distances from EVEs of Pao Bel elements]]
    if overlapStatus == 'yes':
        for name,group in groupedDF:
            distances.append(group['OverlapCategory'])
            names.append(name)
    if overlapStatus == 'no':
        for name,group in groupedDF:
            distances.append(group['Distance'])
            names.append(name)

    #Order by number of hits per category. So can plot stacked histogram with bigest bar on bottom. Copied from stack overflow.
    distances_ordered,names_ordered=zip(*sorted(zip(distances,names), key = lambda count:len(count[0]),reverse=True))
    try: #I needed to add this try/except part because when EVERY enrtry in distances_ordered only has one entry, I get a key error and I don't know how to get rid of it.
        counts, bins, patches = plt.hist(distances_ordered, bins=numBins, label=names_ordered, range=(distanceCutoffs[0],distanceCutoffs[1]))
    except KeyError: #If the Keyerror mentioned above occurs, append a single extra value to the first category of distances. This distance is greater than the upper boundry extablished by 'distanceCutoffs', and so shouldn't affect the histogram
        distances[0] = distances[0].append(pd.Series([distanceCutoffs[1]+5000]))
        distances_ordered,names_ordered=zip(*sorted(zip(distances,names), key = lambda count:len(count[0]),reverse=True))
        counts, bins, patches = plt.hist(distances_ordered, label=names_ordered,bins=numBins, range=(distanceCutoffs[0],distanceCutoffs[1]))

    return [counts,bins,names_ordered]


def graphBarHist (counts, bins,orientation,names,overlapStatus):
    
    for i in range(0,len(names)):
        if orientation == 'opposite' and len(names)>1:
            counts[i] = counts[i]*-1 #Plot 'histogram' upside-down for opposite oriented TE/EVE pairs
        if orientation == 'opposite' and len(names)==1:#Need all of these extra if statements b/c if only 1 value, 'counts' is a list of each count, rather than a list of arrays of each count. So count[0] gave count in the first bin, rather than array of counts for first name in 'names'
                counts = counts * -1
        if overlapStatus == "yes":
            if len(names)==1:
               ax1.bar(bins[1:],counts, align = 'center', color=colorDict[names[i]], label = names[i])
            if len(names)>1:
                ax1.bar(bins[1:],counts[i], align = 'center', bottom = sum(counts[0:i]), color=colorDict[names[i]], label = names[i]) #To make a stacked plot, the bottom of current part, is all previous counts added up
            ax1.set_xlim((-5, 5))
        if overlapStatus == "no":
            if len(names) == 1:
                ax1.bar(bins[1:], counts, width=distanceCutoffs[0] / 10, align='edge', color=colorDict[names[i]], label=names[i])
            if len(names) > 1:
                ax1.bar(bins[1:], counts[i], width=distanceCutoffs[0] / 10, align='edge', bottom=sum(counts[0:i]), color=colorDict[names[i]], label=names[i])
            ax1.set_xlim((distanceCutoffs[0],distanceCutoffs[1]))

    if overlapStatus == "yes":
        ax1.set_xlabel('TE-EVE Overlap Type')
        plt.xticks(bins[1:], ['Upstream','EVEsurroundTE','TEsurroundEVE','Downstream'], rotation = 45)
    if overlapStatus == "no":
        ax1.set_xlabel('Distance from EVE')

    ax1.set_ylim((histYlims[0], histYlims[1]))
    ax1.set_ylabel('Counts')
    ax1.grid(False)


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

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#-------------------Graph overlap only-------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------

if analysisType == "overlapOnly" or analysisType == "overlapOrNearest":
    distanceCutoffs = [-2,2]
    histYlims = [-150,150]

    #Filter by removing non-overlapping hits
    upANDdownstreamMERGED_Master = TEsNearestEVEs[TEsNearestEVEs["Distance"] == 0]

    #Filter based on distance cutoffs
    upANDdownstreamMERGED_Master = upANDdownstreamMERGED_Master[(upANDdownstreamMERGED_Master['Distance']>=distanceCutoffs[0]) & (upANDdownstreamMERGED_Master['Distance']<=distanceCutoffs[1])]

    #--------------------------------------------------------------------------------
    #If desired, filter by desired class, establisehd by filteredBy variable at top of script.
    if filteredBy!= 'NONE':
        upANDdownstreamMERGED_Master = upANDdownstreamMERGED_Master[upANDdownstreamMERGED_Master[filteredByCategory]==filteredBy]



    sameData = getHistData(upANDdownstreamMERGED_Master, 'same','yes')
    sameData_counts = sameData[0]
    sameData_bins = sameData[1]
    sameData_names = sameData[2]

    oppData = getHistData(upANDdownstreamMERGED_Master, 'opposite', 'yes')
    oppData_counts = oppData[0]
    oppData_bins= oppData[1]
    oppData_names= oppData[2]

    #Get a defined color for each element in selected category
    N = len(set(upANDdownstreamMERGED_Master[currentClassification]))
    if filteredBy=="NONE" and currentClassification=="TEclass":
        # For TE class specific coloring (same colorscheme as in Patrick's figures)
        new_colors = [(0.68899655751153521, 0.8681737867056154, 0.54376011946622071),
                      (0.12572087695201239, 0.47323337360924367, 0.707327968232772),
                      (0.65098041296005249, 0.80784314870834351, 0.89019608497619629),
                      (0.98320646005518297, 0.5980161709820524, 0.59423301088459368),
                      (0.21171857311445125, 0.63326415104024547, 0.1812226118410335),
                      (0.89059593116535862, 0.10449827132271793, 0.11108035462744099),
                      (0.78329874347238004, 0.68724338552531095, 0.8336793640080622),
                      (0.99175701702342312, 0.74648213716698619, 0.43401768935077328),
                      (0.99990772780250103, 0.50099192647372981, 0.0051211073118098693),
                      (0.42485198495434734, 0.2511495584950722, 0.60386007743723258)]
        colorDict = dict(zip(['LTR', 'LINE', 'DNA', 'MITEs', 'Unknown', 'UD', 'Penelope', 'Helitrons', 'SINE', 'RC'], new_colors))
    else:
    # For using normal color brewer Paired palette
        sample_colors = sns.color_palette("Paired", N)
        colorDict = dict(zip(pd.unique(sameData_names + oppData_names), sample_colors))

    fig=plt.figure(figsize=(8.5,11), facecolor='white')
    ax1= fig.add_subplot(1,1,1)  #says use 1 row,1 column for plotting area, and insert current graph into position 1
    graphBarHist(sameData_counts,sameData_bins,'same',sameData_names, 'yes')
    graphBarHist(oppData_counts,oppData_bins,'opposite',oppData_names, 'yes')
    # Make legend only reflect up to top 10 Concordant TEs
    plt.legend(sameData_names[0:10], fontsize=11)

    ax1.axhline(linewidth=1, color="black")
    ax1.yaxis.set_ticks_position('left')
    ax1.xaxis.set_ticks_position('bottom')
    ax1.xaxis.set_tick_params(width=2)
    ax1.yaxis.set_tick_params(width=2)

    plt.show()
    fig.savefig(
        outputdir + 'ClassifiedBy_' + currentClassification + '_FilteredBy_' + filteredBy + '_BothStrands_' + analysisType + '_Overlap.pdf',
        dpi=600, facecolor=fig.get_facecolor(), transparent=True)



#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#-------------------Graph nearest only-------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------

if analysisType == "nearestOnly" or analysisType == "overlapOrNearest":
    distanceCutoffs = [-20000,20000]
    histYlims = [-400,400]

    #Filter based on distance cutoffs
    upANDdownstreamMERGED_Master = TEsNearestEVEs[(TEsNearestEVEs['Distance']>=distanceCutoffs[0]) & (TEsNearestEVEs['Distance']<=distanceCutoffs[1])]

    #Filter by removing all overlapping hits
    upANDdownstreamMERGED_Master = upANDdownstreamMERGED_Master[upANDdownstreamMERGED_Master["Distance"] != 0]

    #--------------------------------------------------------------------------------
    #If desired, filter by desired class, establisehd by filteredBy variable at top of script.
    if filteredBy!= 'NONE':
        upANDdownstreamMERGED_Master = upANDdownstreamMERGED_Master[upANDdownstreamMERGED_Master[filteredByCategory]==filteredBy]



    sameData = getHistData(upANDdownstreamMERGED_Master, 'same', 'no')
    sameData_counts = sameData[0]
    sameData_bins = sameData[1]
    sameData_names = sameData[2]

    oppData = getHistData(upANDdownstreamMERGED_Master, 'opposite', 'no')
    oppData_counts = oppData[0]
    oppData_bins= oppData[1]
    oppData_names= oppData[2]

    #Get a defined color for each element in selected category
    N = len(set(upANDdownstreamMERGED_Master[currentClassification]))

    if filteredBy=="NONE" and currentClassification=="TEclass":
        # For TE class specific coloring (same colorscheme as in Patrick's figures)
        new_colors = [(0.68899655751153521, 0.8681737867056154, 0.54376011946622071),
                      (0.12572087695201239, 0.47323337360924367, 0.707327968232772),
                      (0.65098041296005249, 0.80784314870834351, 0.89019608497619629),
                      (0.98320646005518297, 0.5980161709820524, 0.59423301088459368),
                      (0.21171857311445125, 0.63326415104024547, 0.1812226118410335),
                      (0.89059593116535862, 0.10449827132271793, 0.11108035462744099),
                      (0.78329874347238004, 0.68724338552531095, 0.8336793640080622),
                      (0.99175701702342312, 0.74648213716698619, 0.43401768935077328),
                      (0.99990772780250103, 0.50099192647372981, 0.0051211073118098693),
                      (0.42485198495434734, 0.2511495584950722, 0.60386007743723258)]
        colorDict = dict(zip(['LTR', 'LINE', 'DNA', 'MITEs', 'Unknown', 'UD', 'Penelope', 'Helitrons', 'SINE', 'RC'], new_colors))
    else:
    # For using normal color brewer Paired palette
        sample_colors = sns.color_palette("Paired", N)
        colorDict = dict(zip(pd.unique(sameData_names + oppData_names), sample_colors))

    fig=plt.figure(figsize=(8.5,11), facecolor='white')
    ax1= fig.add_subplot(1,1,1)  #says use 1 row,1 column for plotting area, and insert current graph into position 1
    graphBarHist(sameData_counts,sameData_bins,'same',sameData_names, 'no')
    # print(oppData_names)
    # print(oppData_bins)
    # print(oppData_counts)

    graphBarHist(oppData_counts,oppData_bins,'opposite',oppData_names, 'no')
    # Make legend only reflect up to top 10 Concordant TEs
    plt.legend(sameData_names[0:10], fontsize=11)

    ax1.axhline(linewidth=1, color="black")
    ax1.yaxis.set_ticks_position('left')
    ax1.xaxis.set_ticks_position('bottom')
    ax1.xaxis.set_tick_params(width=2)
    ax1.yaxis.set_tick_params(width=2)

    plt.show()
    fig.savefig(
        outputdir + 'ClassifiedBy_' + currentClassification + '_FilteredBy_' + filteredBy + '_BothStrands' + analysisType + '_NoOverlap.pdf',
        dpi=600, facecolor=fig.get_facecolor(), transparent=True)