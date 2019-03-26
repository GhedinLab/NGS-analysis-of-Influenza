#!/bin/sh
#
#SBATCH --workdir=/scratch/ts2742/hk/scripts
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --time=32:00:00
#SBATCH --mem=8GB
#SBATCH --job-name=rr
#SBATCH --mail-type=END
#SBATCH --mail-user=ts2742@nyu.edu
#SBATCH --error=rr.error
#SBATCH --output=rr.out

module purge
module load pysam/intel/0.10.0

python readreport_v4_2.py --ref $ref --strain $strain --infile $infile 
