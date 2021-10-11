# Mutation spectrum

 - mutation_spectrum.py --> uses CpG.r to found for all species the mutations that are located in CpG sites.
 - type.r --> uses sp_all_type.txt (from the above script) to plot the type per group of species

# Mutation bias

 - Poo_0.py --> found the parents of origin for all samples, all SNPs, per chromosme. Need the sp.txt to find in which species uses POOHA.py
 - Poo_1.py --> catenate all the chromosomes in one file, grep the mutation in a file
 - Poo_2.py --> get the de novo mutations phased, uses poo.r
 - bias.r --> uses data.txt (from the above scripts) to get the mutation bias per group

# Parental ages

 - age.r --> correlation of mutation rate per generation and parental ages
