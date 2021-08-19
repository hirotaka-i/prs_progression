import sys
import os
import pandas as pd
import numpy as np
import argparse
import subprocess
import glob

"""
This fuction creates the prs scores from MINIMAC imputed output
"""

parser = argparse.ArgumentParser(description='Arguments')
parser.add_argument('--dataset', type=str, default='nope', help='dataset name')
parser.add_argument('--pathway', type=str, default='nope', help='pathway name')
args = parser.parse_args()



# parameter
dataset=args.dataset
pathway=args.pathway
gfolder=f'/data/CARD/PD/imputed_data/{dataset}'
exclusion_regions='/data/CARD/GENERAL/exclusion_regions_hg19.txt' # stored in the common resource

wkdir="/data/CARD/projects/prs_prog"
prs_ref=f'{wkdir}/PD_meta_analysis_PDGENE_PDWBS/PD_meta_analysis_PDGENE_PDWBS.txt' # downloaded from http://research-pub.gene.com/chang_et_al_2017/PD_meta_analysis_PDGENE_PDWBS.tar.gz
pathway_folder=f'{wkdir}/pathways_extract_hg19' # provided by manuela
temp=f'{wkdir}/temp/{pathway}/{dataset}'

plink2='/data/CARD/PD/AMP-PD/Plink/plink2_dev_June8/plink2'
prsice2='/data/CARD/projects/prs_prog/PRSice_linux'


# shell submitting functiono

def shell_do(command, log=False, return_log=False):
    print(f'Executing: {(" ").join(command.split())}', file=sys.stderr)

    res=subprocess.run(command.split(), stdout=subprocess.PIPE)

    if log:
        print(res.stdout.decode('utf-8'))
    if return_log:
        return(res.stdout.decode('utf-8'))

# Script Body

## creat a temp folder
shell_do(f'mkdir -p {temp}')


## slice relevant regions from choromosome-level vcfs
for num in range(1,23):
    chrnum=f'chr{num}'

    cmd = f"""\
{plink2} --vcf {gfolder}/{chrnum}.dose.vcf.gz dosage=DS \
--var-filter \
--extract-if-info R2 >= 0.3 \
--extract range {pathway_folder}/{pathway}_extract.txt \
--exclude range {exclusion_regions} --make-pgen --out {temp}/{chrnum} \
--rm-dup force-first \
--threads 1 \
--geno 0.05"""

    shell_do(cmd)

        
## combine chromosome level vcfs to bgen (which can use dosages with PRSice)
pvarlist=glob.glob(f'{temp}/chr*.pvar')
pfiles=[i.replace('.pvar', '') for i in pvarlist]

with open(f'{temp}/merge_list.txt', 'w') as f:
    f.write('\n'.join(pfiles))

cmd = f"{plink2} --pmerge-list {temp}/merge_list.txt --make-bed --geno 0.05 --export bgen-1.2 --threads 1 --out {temp}/all"

shell_do(cmd)

## coduct PRSice
cmd = f"""\
Rscript {prsice2}/PRSice.R \
--prsice {prsice2}/PRSice_linux \
--base {wkdir}/prs_base.txt \
--clump-kb 250 \
--clump-r2 0.1 \
--or \
--target {temp}/all,{temp}/all.sample \
--type bgen \
--no-regress \
--lower 5e-08 \
--upper 0.05 \
--interval 5e-08 \
--out {temp}/res"""

shell_do(cmd)

print(f'FINISH: The scores are stored at {temp}/res-')