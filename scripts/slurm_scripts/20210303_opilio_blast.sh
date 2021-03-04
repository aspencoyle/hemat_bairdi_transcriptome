#!/bin/bash
## Job Name
#SBATCH --job-name=afcoyle_opilioblast
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


/gscratch/srlab/programs/ncbi-blast-2.8.1+/bin/blastx \
-task="blastx" \
-query /gscratch/srlab/afcoyle/projects/21_ncbi_genome_blasts/data/cbai_hemat_transcriptomev2.0.fasta \
-db /gscratch/srlab/afcoyle/projects/21_ncbi_genome_blasts/output/blastdbs/opilio/opilio_blastdb \
-out /gscratch/srlab/afcoyle/projects/21_ncbi_genome_blasts/output/blastres/opilio_highereval_blastres.tab \
-evalue 1E-03 \
-num_threads 40 \
-max_target_seqs 1 \
-outfmt 6
