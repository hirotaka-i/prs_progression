# Pathway specific PRS analysis for PD progression

 - **Project:** Pathway specific PRS analysis for PD progression. Collaboration with Manuela Tan, Lasse Pihlstrom (Oslo University Hospital)
 - **Author(s):** Hirotaka Iwaki 
 - **Date Last Updated:** August 20 2021
  - **Status:** Incomplete
    - **Update Description:** Starting analysis



## Overview
In this analysis, we create PD-PRS with the sliced SNPs close to the genes specific to PD-related pathways. We call them as *pathwas-specific-prs*. 
Then we test whether these pathway-specific-prs are associated with the PD progression using various PD progression markers. 

## Cohorts
cohorts=['DATATOP', 'PreCEPT_PostCEPT', 'DIGPD', 'PICNICS', 'PROPARK', 'NET_PD_LS1', 'PARKFIT', 'PARKWEST']


## Pathways
pathways=['endocytosis', 'adaptive_immune', 'innate_immune', 'lysosome', 'mitochondria', 'microglia', 'monocytes', 'alpha_synuclein']

##  Target gwas to create PRS
Using Chang 2017 GWAS discovery to create PRSs: http://research-pub.gene.com/chang_et_al_2017/

P-value thresholds are set to be the maximum_P < 0.05. The PRS will be normalized (z-transformed) per cohort

## Outcomes
outcomes=['age_at_onset', 'UPDRS2', 'UPDRS3', 'HY3', 'MOCA', 'DEMENTIA', 'PDQ8', 'RBD']

## Script
Please refer to the Main.ipynb