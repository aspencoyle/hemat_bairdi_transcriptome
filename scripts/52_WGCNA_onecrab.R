#############################
# Aidan Coyle, afcoyle@uw.edu
# Roberts lab, UW-SAFS
# 2021-04-12
#############################

# This file is an initial attempt to run WGCNA on a single crab over time.

# We will look solely at C. bairdi genes for Crab A

# This corresponds to Library IDs 178, 359, and 463 (Days 0, 2, and 17) for that crab

# We will first extract the TPM (transcripts per million) counts from the kallisto 
# libraries created earlier in the pipeline. We will then change those to logTPM 
# counts, and then begin the WGCNA analysis


# Load all libraries


