#!/bin/sh
#
#SBATCH --workdir=/scratch/ts2742/hk/scripts
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --time=32:00:00
#SBATCH --mem=32GB
#SBATCH --job-name=BOWTIE
#SBATCH --mail-type=END
#SBATCH --mail-user=ts2742@nyu.edu
#SBATCH --error=er.error
#SBATCH --output=out.out

module purge
module load pysam/intel/0.10.0
module load bowtie2/intel/2.3.2
module load samtools/intel/1.6
module load picard/2.8.2

python bowtie2_align.py -m1 $a -m2 $b --name $c --strain $d --refpath $e
