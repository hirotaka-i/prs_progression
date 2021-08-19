# Rscript --vanilla analyzePRS.R {dataset} {cohort} {pathway} {outcome} {model}
# This script test the association between the pathway-specific-prs and the outcome
# The result will be stored at report/{model}.csv

## libraries
suppressMessages(library(data.table))
suppressMessages(library(tidyverse))
suppressMessages(library(lmerTest))
suppressMessages(library(survival))


# arguments
args = commandArgs(trailingOnly=TRUE)
dataset=args[1]
cohort =args[2]
pathway=args[3]
outcome=args[4]
model=args[5] # cs, lt, surv

## dataset='DIGPD_chip'
## cohort = 'DIGPD'
## pathway='alpha_synuclein'
## outcome='HY3'
## model='surv' # cs, lt, surv



# parameters
clin_path='/data/LNG/iwakih2/MaleFemale/analysisSet.csv'
resdir="/data/CARD/projects/prs_prog/report"
prsfile=paste0('/data/CARD/projects/prs_prog/prs/', pathway, '.', dataset, '.csv')
pcafile=paste0('/data/CARD/PD/genotype_data/', dataset,'/', dataset,'_pc.eigenvec')
outfile=paste0(resdir, '/',  model, '.csv')

# create an analysis dataset
dclin = fread(clin_path) %>% filter(STUDY_NAME==cohort)
dprs = fread(prsfile)
dpca = fread(pcafile,  sep='\t')
dclin$IID = tolower(paste0(dclin$ID, '_', dclin$ID)) 
dpca$IID = tolower(paste0(dpca$IID, '_', dpca$IID))
dprs$IID= tolower(dprs$IID)
d = inner_join(dclin, dprs, by='IID') %>% 
    inner_join(., dpca, by='IID')
d$y = d[,outcome,with=FALSE]



if(model=='surv'){ 
    # create an output file
    if(!file.exists(outfile)){
        cols=c('dataset', 'outcome', 'pathway', 'Coeff', 'se', 'Pvalue', 'Cox.zphPVal', 'N', 'ovlik.ratio', 'logrank', 'r2')
        write(paste(cols, collapse=','),outfile, append=F)
    }

    
    # reduce to the complete set
    t = d[,c('y', 'IID', 'PRS_z', 'FEMALE','DiseaseDuration','AAO', 
             'PC1', 'PC2', 'PC3', 'PC4', 'PC5')]
    dt = t[complete.cases(t)]
    dt_last0  = dt %>% filter(y==0) %>% arrange(desc(DiseaseDuration)) %>% distinct(IID, .keep_all=T)
    dt_first1 = dt %>% filter(y==1) %>% arrange(DiseaseDuration) %>% distinct(IID, .keep_all=T)
    dfsurv=bind_rows(dt_first1, dt_last0) %>% distinct(IID, .keep_all=T)
    N = dim(dfsurv)[1]
    
    
    
    if(N == 0){
        
        print(paste0(outcome, ' is NA for ', dataset))
        
    }else{
        
        # analysis
        m1 = coxph(Surv(time=DiseaseDuration, event=y) ~ PRS_z + FEMALE + AAO 
                  + PC1 + PC2 + PC3 + PC4 + PC5, data=dfsurv)
        kmz <- cox.zph(m1, transform = "km")


        # resulut
        summ=summary(m1)
        res = coef(summ)
        resline=paste(paste(c(
            dataset, outcome, pathway, 
            res['PRS_z',c('coef', 'se(coef)', 'Pr(>|z|)')], 
            kmz$table['PRS_z', 'p'], N, 
            summ$logtest[[1]], summ$sctest[[1]], summ$rsq[[1]]
        ), collapse=','))
        write(resline, outfile, append=T)
        print(paste(model, resline))
    }
}


if(model=='lt'){ 
    # create an output file
    if(!file.exists(outfile)){
        cols=c('dataset', 'outcome', 'pathway', 'Coeff', 'se', 'df', 't_value', 'Pvalue', 'N_obs', 'N_inds')
        write(paste(cols, collapse=','),outfile, append=F)
    }

    
    # reduce to the complete set
    t = d[,c('y', 'IID', 'PRS_z', 'FEMALE','DiseaseDuration','AAO',
             'PC1', 'PC2', 'PC3', 'PC4', 'PC5')]
    dt = t[complete.cases(t)]
    N_obs = dim(dt)[1]
    N_inds= length(unique(dt$IID))
    
    if(N_inds==0){
        
        print(paste0(outcome, ' is NA for ', dataset))
        
    }else{
        
        # analysis
        m1 = lmer(y ~ PRS_z*DiseaseDuration + FEMALE*DiseaseDuration + AAO 
                  + PC1 + PC2 + PC3 + PC4 + PC5 + (DiseaseDuration|IID), data=dt)
        
        # resulut
        res = coef(summary(m1))
        resline=paste(paste(c(
            dataset, outcome, pathway, 
            res['PRS_z:DiseaseDuration',c('Estimate', 'Std. Error', 'df', 't value', 'Pr(>|t|)')], N_obs, N_inds
        ), collapse=','))
        write(resline, outfile, append=T)
        print(paste(model, resline))
    }
}


if(model=='cs'){ 
    # create an output file
    if(!file.exists(outfile)){
        cols=c('dataset', 'outcome', 'pathway', 'Coeff', 'se','Pvalue', 'N')
        write(paste(cols, collapse=','),outfile, append=F)
    }
    
    # reduce to the complete set
    t = d[,c('y', 'IID', 'PRS_z', 'FEMALE', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5')]
    t = t[complete.cases(t)]
    
    if(dim(t)[1] == 0){
        
        print(paste0(outcome, ' is NA for ', dataset))
        
    }else{
        
        # analysis
        dt = t %>% arrange(desc(y)) %>% distinct(IID, .keep_all=T)
        m1 = glm(y ~ PRS_z + FEMALE* + PC1 + PC2 + PC3 + PC4 + PC5, data=dt)
        
        # resulut
        N = dim(dt)[1]
        res = coef(summary(m1))
        resline=paste(paste(c(
            dataset, outcome, pathway, res['PRS_z',c('Estimate', 'Std. Error', 'Pr(>|t|)')], N
        ), collapse=','))
        write(resline, outfile, append=T)
        print(paste(model, resline))
    }
}