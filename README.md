# Germline mutation rate in vertebrates

Method used in the study on germline mutation rate estimation in 68 species of vertebrates.

These codes were used to estimate the germline mutation rates of 68 species of vertebrates from NGS data in pedigree samples (mother, father and offspring). 
The primary analysis was presented in https://github.com/lucieabergeron/germline_mutation_rate. 
We present here the difference in methodology and the codes to produce the figures and downstream analysis. 
The raw sequences can be found on NCBI under the project number PRJNA767781.

# Analysis

0. germline_mutation_rate discribe the additional step that has been conducted to produce the per generation and yearly rates.
1. phylogeny contains the analysis conducted to create the UCE phylogeny and to calibrate the tree.
2. mutations_analysis contains the script to analysis the DNM characteristics such as the mutational spectrum, the sex bias or the impact of the parental ages.

# Requirements

 - samtools 1.2
 - bcftools 1.2
 - bcftools 1.9
 - R 3.5.1
 - IQTREE
 - Phyluce


