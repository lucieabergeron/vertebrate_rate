# Heterozygosity

 - heterozygosity.py --> find the heterozygosity for each individual each species (uses bcftools 1.9)

# Nucleotide diversity

 - pi.py --> for each species, infer a global estimate of the site frequency spectrum (SFS) using a maximum likelihood method, estimate the observed nucleotide diversity per site, and finally, estimate the nucleotide diversity of a species on a slidding window of 50kb and a step size of 10kb. An R script average the pairwise estimation of the nucleotide diversity to have an estimate of π. 

# Effective population size

 - psmc.py --> estimate an effective population size for all species. 
 - psmc_plot.py --> uses plot_results_all.py and estime and plot the effective population size
