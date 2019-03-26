# NGS Scripts
These scripts take paired end fastqs, aligns them, and discovers minor variants.

## Method Overview
Paired end fastqs are trimmed using Trimmomatic[1] annd then aligned to their respective reference (H1N1pdm or seasonal H3N2) with bowtie2[2]. Samtools[3] was used to sort and index the bams which had their PCR duplicates removed with Picard MarkDuplicates. Python scripts using pysam parsed the bams to determine the frequencies of major and minor variants of each position along the reference sequence, as well check for strand bias. These are used to calculate the genetic distance (L1-norm) between every other sample to create a dissimilarity matrix. Samples were categorized by part of the same family (household) or not (other). R ggplot was used to graph boxplots of the genetic distances. 


1. Bolger, A. M., Lohse, M., & Usadel, B. (2014). Trimmomatic: A flexible trimmer for Illumina Sequence Data. Bioinformatics, btu170.
2. Langmead B, Salzberg S. Fast gapped-read alignment with Bowtie 2. Nature Methods. 2012, 9:357-359.
3. Li H.*, Handsaker B.*, Wysoker A., Fennell T., Ruan J., Homer N., Marth G., Abecasis G., Durbin R. and 1000 Genome Project Data Processing Subgroup (2009) The Sequence alignment/map (SAM) format and SAMtools. Bioinformatics, 25, 2078-9. [PMID: 19505943]

### Running alignment
1. `slurm_master.py` will kick off jobs to the cluster. This will submit a job for each mate1 and mate2 pair of a sample to expedite the process.
2. `bowtie.sh` is the shell script being kicked off by the above.
3. `bowtie2_align.py` will trim, perform alignemnt, and remove duplicates to create bamfiles.

### Creating minor variant lists
1. `slurm_master_readrepoert.py` will kick off jobs to the cluster so jobs can be run concurrently.
2. `readreport_runner.sh` is the shell script taking in the above.
3. `readreport_v4_2.py` parses bamfiles with pysam and the frequency of major and minor alleles, and performs a binomial distribution check to determine if there is strand bias. Defaults are coverage of 200 and minor variant frequency > 3%

### Plotting genetic distance
1. `l1norm.py` creates a dissimilarity matrix of the samples using the L1-norm.
2. `parse_dissim.py` prepares the output of above to be used by the R scripts
3. `plotter.R` plots boxplots of genetic distances.
