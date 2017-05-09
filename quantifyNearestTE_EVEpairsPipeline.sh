#!/bin/bash

SCRIPTDIRECTORY="/home/zwhitfield/Desktop/ForMarkGenomePaper/ScriptsFinalVersions/WithArguments/publicRelease/forHistogramQuantificationAndGraphing"

#INPUTDIRECTORY="/home/zwhitfield/Desktop/ForMarkGenomePaper/FrozenData/Aag2_assembly/"
INPUTDIRECTORY="/home/zwhitfield/Desktop/ForMarkGenomePaper/FrozenData/Aag2_assembly/OLD_USING_PREV_EVE_COORD_FILE/"

FILTEREDBY="Ty3_gypsy" #Enter 'NONE' if don't want any filtering. Specifies category to filter TEs by. Will only keep TEs of this type for analysis.
FILTERDBYCATEGORY="TEfamily" #So script knows what 'column' of data frame to filter by
GROUPBY="TEdescription" # Specifies what level/category to quantify TEs by (ie what goes into groupBy function of pandas) for statistics and plotting on resulting histogram. Use 'CombinedGroup' if grouping by combined categories. Need to specify various categories to be combined in actual python script.

#Should the resulting EVE-TE pairs be...
# "overlapOrNearest": If an EVE has an overlapping TE(s), then use those pairs. If no overlapping TE, then look at 'nearest neighbor' TE (both upstream AND downstream).
# "nearestOnly": Ignore all TEs which overlap EVEs, and find 'nearest neighbor' TEs for all EVEs, both upstream and downstream.
# "overlapOnly": Only look at EVEs with a TE whose coordinates overlap. If no overlapping TE, that EVE is not used in the analysis.
ANALYSISTYPE="overlapOrNearest" # "overlapOrNearest" OR "nearestOnly" OR "overlapOnly".

#----------------------------------------------------------------------------------------------------------------------
#Get genome-wide counts of TE type of interest to use as comparision for statistical enrichment. 
#----------------------------------------------------------------------------------------------------------------------
printf "Quantifying TEs genome-wide..."
printf "\n"
python ${SCRIPTDIRECTORY}/"NearestEVEquantification_GenomeWide_pandasBash_NowWithStats_FrozenDataNoTEfam.py" ${INPUTDIRECTORY} ${INPUTDIRECTORY} ${FILTEREDBY} ${FILTERDBYCATEGORY} ${GROUPBY}
printf "Done"
printf "\n"

#----------------------------------------------------------------------------------------------------------------------
#Generate files of EVE-TE pairs type of analysis defined by $ANALYSISTYPE
#----------------------------------------------------------------------------------------------------------------------
printf "Quantifying EVE-TE pairs by specified analysis type..."
printf "\n"
python ${SCRIPTDIRECTORY}/"NearestEVEquantification_pandas_createFiles.py" ${INPUTDIRECTORY} ${INPUTDIRECTORY} ${ANALYSISTYPE}
printf "Done"
printf "\n"

#----------------------------------------------------------------------------------------------------------------------
#Stats only. This combines stats for overlap and nearest TE(up and downstream). Difference from before is nearest TE calculation is only for EVEs WITHOUT an overlapping TE.
#----------------------------------------------------------------------------------------------------------------------
OUTPUTDIRECTORY="/home/zwhitfield/Desktop/ForMarkGenomePaper/FrozenData/Aag2_assembly/OLD_USING_PREV_EVE_COORD_FILE/ForStats/"

printf "Quantifying enrichment of EVE-TE pairs compared to genome-wide background..."
printf "\n"
python ${SCRIPTDIRECTORY}/"NearestEVEquantification_pandas_onlyStats.py" ${INPUTDIRECTORY} ${OUTPUTDIRECTORY} ${FILTEREDBY} ${FILTERDBYCATEGORY} ${GROUPBY} ${ANALYSISTYPE}
printf "Done"
printf "\n"

#----------------------------------------------------------------------------------------------------------------------
#Histograms only. This graphs both overlap and nearest TE(up and downstream; if applicable). NEAREST TEs are graphed by distance and overlapping TEs are plotted categoricaly.
#----------------------------------------------------------------------------------------------------------------------
OUTPUTDIRECTORY="/home/zwhitfield/Desktop/ForMarkGenomePaper/FrozenData/Aag2_assembly/OLD_USING_PREV_EVE_COORD_FILE/Figures/"

printf "Graphing categories of TEs nearest EVEs..."
printf "\n"
python ${SCRIPTDIRECTORY}/"NearestEVEquantification_pandas_createHistograms.py" ${INPUTDIRECTORY} ${OUTPUTDIRECTORY} ${FILTEREDBY} ${FILTERDBYCATEGORY} ${GROUPBY} ${ANALYSISTYPE}
printf "Done"
printf "\n"

