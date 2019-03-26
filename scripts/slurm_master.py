import os
import subprocess
from subprocess import call
import glob
import pandas as pd
import argparse
import sys

if __name__ == "__main__":
        
    parser = argparse.ArgumentParser()
    parser.add_argument('--raw_files','-R',required=True,help='raw files directory')
    # /scratch/cgsb/gencore/out/Ghedin/2018-09-18_HMLHWBGX7/merged
    args = parser.parse_args()
    n1list = glob.glob(os.path.join(args.raw_files,'*_n01_*.fastq.gz' ))
    n2list = glob.glob(os.path.join(args.raw_files,'*_n02_*.fastq.gz' ))
    n1list.sort()
    n2list.sort()
    
    h3n2list = ['662','669','671','672','674','688','689','689','689','692','692','695','703','707','709','710','714','720','720','720','720','724','724','725','726','727','729','731','733','734','734','734','736','737','739','739','739','739','741','747','747','747','750','752','755','755','756','763','763','764','768','769','770','781','783']

    h1n1list = ['661','663','665','667','667','678','681','681','681','683','683','684','684','693','693','698','698','701','704','708','712','712','712','712','715','721','722','722','738','742','742','745','745','746','746','751','751','751','751','751','751','757','758','759','761','762','765','770','771','771','772','774','779','779','779','779','779','785','787','793','793']

    df = pd.read_csv('../files/HK_metadata_v2.csv')
    rootdict = {}
    for idx,row in df.iterrows():
        if row['Household'] in h3n2list:
            useref = '../reference/h3n2_cds'
            strain = 'H3N2'
        elif row['Household'] in h1n1list:
            useref = '../reference/h1n1_cds'
            strain = 'H1N1'
        else:
            continue
        rootid = row['realid'].split('_')[-1]
        rootdict[rootid] = [row['Fastq_Name'],strain,useref]
        

    for n1,n2 in zip(n1list,n2list):
        rootkey = n1.split('/')[-1].split('.')[0].split('_')[-1]
        try:
            a = rootdict[rootkey]
            os.system('sbatch --export=a="%s",b="%s",c="%s",d="%s",e="%s" bowtie.sh' % (n1,n2,a[0],a[1],a[2]))
            print('sbatch --export=a="%s",b="%s",c="%s",d="%s",e="%s" bowtie.sh' % (n1,n2,a[0],a[1],a[2]))
        except KeyError:
            pass