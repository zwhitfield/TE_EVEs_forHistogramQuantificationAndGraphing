# -*- coding: utf-8 -*-
"""
Created on Wed Oct  5 16:35:23 2016

@author: zwhitfield
"""
import sys
import pandas as pd
import numpy as np
#np.__version__
import matplotlib.pyplot as plt
import matplotlib
#matplotlib.__version__
matplotlib.style.use('ggplot')

inputdir = str(sys.argv[1])
outputdir = str(sys.argv[2])

#This will be used to only select specific entries in dataset with TEs of the given type.
#Can be at any taxonomy level (class, subclass, family, element, etc...)
#The specific level is specified by filteredByCategory (ie name of column to subset using filteredBy)
filteredBy = str(sys.argv[3])
filteredByCategory = str(sys.argv[4])
groupedCategory =  str(sys.argv[5])

#outputdir="/home/zwhitfield/Desktop/ForMarkGenomePaper/FrozenData/Aag2_assembly/"
#filteredBy="NONE" #Enter 'NONE' if don't want any filtering.
#filteredByCategory="TEfamily"
#groupedCategory="TEclass"

#-----------------------------------------------------------------------------------------------------------
#Load ALL EVEs in genome. To have a baseline proportion for comparison
allTEsInGenome = pd.read_csv(outputdir + "Aag2_Contigs_TEs_sorted.bed", 
                              sep="\t",
                              names = ["TEcontig","TEstart","TEend","TEdescription","TEvalue","TEstrand", "TEfamily"])
#First, rename all LTR/Gypsy as LTR/Ty3_gyspsy. THere are both in the file, so need a consensus
#First, rename all LTR/Copia as LTR/Ty1_copia. THere are both in the file, so need a consensus
#First, rename all LTR/Pao as LTR/Pao_Bel. THere are both in the file, so need a consensus
allTEsInGenome.loc[allTEsInGenome.TEfamily == 'LTR/Gypsy', 'TEfamily'] = "LTR/Ty3_gypsy"
allTEsInGenome.loc[allTEsInGenome.TEfamily == 'LTR/Copia', 'TEfamily'] = "LTR/Ty1_copia"
allTEsInGenome.loc[allTEsInGenome.TEfamily == 'LTR/Pao', 'TEfamily'] = "LTR/Pao_Bel"


allTEsInGenome["TEclass"]=allTEsInGenome["TEfamily"].str.split('/').str[0]
allTEsInGenome["TEfamily"]=allTEsInGenome["TEfamily"].str.split('/').str[1]
allTEsInGenome['Distance'] = 0

#allTEsInGenome = allTEsInGenome.head(500)
allTEsInGenome['CombinedGroup'] = 'Unclassified'
allTEsInGenome.loc[(allTEsInGenome.TEclass == 'DNA') | (allTEsInGenome.TEclass == 'Helitrons')| (allTEsInGenome.TEclass == 'MITEs')| (allTEsInGenome.TEclass == 'RC'), 'CombinedGroup'] = 'DNAall'
allTEsInGenome.loc[(allTEsInGenome.TEclass == 'LTR') | (allTEsInGenome.TEclass == 'LINE')| (allTEsInGenome.TEclass == 'SINE')| (allTEsInGenome.TEclass == 'Penelope'), 'CombinedGroup'] = 'RNAall'


if filteredBy!= 'NONE':
    allTEsInGenome = allTEsInGenome[allTEsInGenome[filteredByCategory]==filteredBy]

classifications = [groupedCategory]

for currentClassification in classifications:
    print ("Filtering by " + filteredBy + " and grouping by " + currentClassification)
    outfile= open(outputdir + 'ClassifiedBy_' + currentClassification+ '_FilteredBy_'+ filteredBy + '_WholeGenome.txt','w')
    outfile.write(currentClassification + '\t' + "Counts" + "\n")
    distances=list()
    names = list()
    groupedDF = allTEsInGenome.groupby(currentClassification)
    
    for name,group in groupedDF:
        distances.append(group['Distance'])
        names.append(name)
        outfile.write(name + '\t' + str(len(group['Distance'])) + "\n")
    # distances_ordered,names_ordered=zip(*sorted(zip(distances,names), key = lambda count:len(count[0]),reverse=True))
    # fig=plt.figure()
    # plt.hist(distances_ordered, label=names_ordered,bins=5, stacked=True)
    # plt.suptitle(currentClassification)
    # plt.legend()
    #fig.savefig(outputdir + 'Figures/ClassifiedBy_' + currentClassification+ '_FilteredBy_'+ filteredBy + '_WholeGenome.png', dpi=600)
    outfile.close()
