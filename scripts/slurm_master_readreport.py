import os
import subprocess
from subprocess import call
import glob
import pandas as pd
import argparse
import sys

if __name__ == "__main__":
        
    parser = argparse.ArgumentParser()
    parser.add_argument('--bam_dir','-b',required=True,help='bam files directory')
    args = parser.parse_args()
    bamlist = glob.glob(os.path.join(args.bam_dir,'*bam' ))


    for bam in bamlist:
        print bam
        # bamfiles/rmdup_bams/698-V2_0_rep1.H1N1.rmd.bam
        rootfile = bam.split('/')[-1]
        strain = rootfile.split('.')[1]
        if strain.upper() == 'H3N2':
            refpath = '../reference/h3n2_cds.fasta'
        elif strain.upper() == 'H1N1':
            refpath = '../reference/h1n1_cds.fasta'
        print('sbatch --export=strain="%s",ref="%s",infile="%s" readreport_runner.sh' % (strain,refpath,bam))
        os.system('sbatch --export=strain="%s",ref="%s",infile="%s" readreport_runner.sh' % (strain,refpath,bam))
