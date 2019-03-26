
import os
import subprocess
from subprocess import call
import glob
import pandas as pd
import argparse
import sys


if __name__ == "__main__":
        # python bowtie2_align.py --ref_dir $refpath --ref $strain -m1 $m1 -m2 $m2 --name $name 
    parser = argparse.ArgumentParser()

    parser.add_argument('--mate1','-m1',required=True,help='Give directory where read 1 is located')
    parser.add_argument('--mate2','-m2',required =True,help='Give directory where read 2 is located')
    parser.add_argument('--refpath',required =True,help='Give directory where read 2 is located')
    parser.add_argument('--strain',required =True,help='Give directory where read 2 is located')
    parser.add_argument('--name',required =True,help='Give directory where read 2 is located')
    args = parser.parse_args()

    mate1 = args.mate1
    mate2 = args.mate2
    refpath = args.refpath
    strain = args.strain
    name = args.name
    newname = name+'.'+strain


    directory_list = ['bamfiles','bamfiles/sorted_bams','bamfiles/rmdup_bams','MetricFiles',
                        'rawfiles/trimmed_mate1','rawfiles/trimmed_mate2','rawfiles/unpair','unmapped']

    for direc in directory_list:
        if not os.path.exists(direc):
            os.makedirs(direc)

    script = [
    "java -jar /share/apps/trimmomatic/0.36/trimmomatic-0.36.jar PE -phred33 -threads 20 {0} {1} ./rawfiles/trimmed_mate1/{2}_trimmed_1.fq ./rawfiles/unpair/{2}.unpair_trimmed_1.fq ./rawfiles/trimmed_mate2/{2}_trimmed_2.fq ./rawfiles/unpair/{2}.unpair_trimmed_2.fq ILLUMINACLIP:/share/apps/trimmomatic/0.36/adapters/NexteraPE-PE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:20".format(mate1, mate2, newname),
    "bowtie2 -x {0} -1 rawfiles/trimmed_mate1/{1}_trimmed_1.fq -2 rawfiles/trimmed_mate2/{1}_trimmed_2.fq -S bamfiles/{1}.sam --no-mixed --local --very-sensitive --un-conc unmapped/unmapped.{1}.fastq".format(refpath,newname),
    "samtools view -bSq 20 bamfiles/{0}.sam > bamfiles/{0}.bam".format(newname),
    "samtools sort -T bamfiles/{0}.sorted -o bamfiles/sorted_bams/{0}.sorted.bam bamfiles/{0}.bam".format(newname),
    "samtools index bamfiles/sorted_bams/{0}.sorted.bam".format(newname),
    "java -jar $PICARD_JAR MarkDuplicates I=bamfiles/sorted_bams/{0}.sorted.bam O=bamfiles/rmdup_bams/{0}.rmd.bam M=MetricFiles/{0}.met.star.txt REMOVE_DUPLICATES=true".format(newname),
    "samtools index bamfiles/rmdup_bams/{0}.rmd.bam".format(newname)]

    for command in script:
        print(command)
        os.system(command)