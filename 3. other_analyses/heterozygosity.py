# -*- coding: utf-8 -*-
"""
This script is to find the heterozygosity per individual.

Lucie Bergeron
03.12.18
"""
##################################################
# What you need ##################################
##################################################

# Packages:
import subprocess
import os
import pandas as pd

# Import species names:
species_sample = pd.read_csv('sp_sa2.txt', sep='\t', index_col=None, header=None)
# Species maybe import in ../ directory the one filtered

# The functions:
def nb(sp_target,sa_target):
    file = open('{}_stat_nb_snps.sh'.format(sp_target),'w')
    file.write('#!/bin/bash \n')
    file.write('#SBATCH --partition express,normal \n')
    file.write('#SBATCH --mem 10G \n')
    file.write('#SBATCH -c 1 \n')
    file.write('#SBATCH --time=00:58:00 \n')
    file.write('bcftools stats -s {} /home/lucie/MammalianMutation/faststorage/{}/vcf_files/genotype_genomicDBI_gather.g.vcf > {}_stats.txt \n'.format(sa_target, sp_target, sp_target))
    file.write('grep \'PSC\' {}_stats.txt >> {}_stats_sub.txt \n'.format(sp_target, sp_target))
    file.write('grep \'TSTV\' {}_stats.txt >> {}_stats_sub.txt \n'.format(sp_target, sp_target))
    file.close()
    ##"""Submit the .sh to the server"""
    sub_cmd = "sbatch -A MammalianMutation -o {}_stat_nb_snps.out {}_stat_nb_snps.sh".format(sp_target, sp_target)
    subprocess.call(sub_cmd, shell=True)


##################################################
# What you run  ##################################
##################################################
# For each species:
for sp in range(0,len(species_sample)):
    sp_target=species_sample.loc[sp,0]
    sa_target=species_sample.loc[sp,1]
    nb(sp_target=sp_target, sa_target=sa_target)

