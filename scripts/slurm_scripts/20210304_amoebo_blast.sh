#!/bin/bash
## Job Name
#SBATCH --job-name=afcoyle_amoeboblast
## Allocation Definition
#SBATCH --account=srlab
#SBATCH --partition=srlab
## Resources
## Nodes
#SBATCH --nodes=1
## Walltime (days-hours:minutes:seconds format)
#SBATCH --time=1-12:00:00
## Memory per node
#SBATCH --mem=120G
##turn on e-mail notification
#SBATCH --mail-type=ALL
#SBATCH --mail-user=afcoyle@uw.edu
## Specify the working directory for this job
#SBATCH --chdir=/gscratch/scrubbed/afcoyle



/gscratch/srlab/programs/ncbi-blast-2.8.1+/bin/blastn \
-task="blastn" \
-query /gscratch/srlab/afcoyle/projects/21_ncbi_genome_blasts/data/cbai_hemat_transcriptomev2.0.fasta \
-db /gscratch/srlab/afcoyle/projects/21_ncbi_genome_blasts/output/blastdbs/amoebo/amoebo_blastdb \
-out /gscratch/srlab/afcoyle/projects/21_ncbi_genome_blasts/output/blastres/amoebo_eval10_2_blastres.tab \
-evalue 1E-02 \
-num_threads 40 \
-max_target_seqs 1 \
-outfmt 6

