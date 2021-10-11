# -*- coding: utf-8 -*-
"""
This script plot the psmc.

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
species = pd.read_csv('parameters.txt', sep='\t', index_col=None, header=None)

##################################################
# What you run  ##################################
##################################################
# For each chromosome one function:
for sp in range(0,len(species)):
    sp_target=species.loc[sp,0]
    gt=species.loc[sp,1]
    mu=species.loc[sp,2]
    sample=species.loc[sp,3]
    ## Script to run python
    file = open('{}/psmc_plot.sh'.format(sp_target),'a')
    file.write('#!/bin/bash \n')
    file.write('#SBATCH --mem 10G \n')
    file.write('#SBATCH -c 1 \n')
    file.write('#SBATCH --time=00:30:00 \n')
    file.write('## get seq \n')
    file.write('source activate matplot \n')
    file.write('python ~/MammalianMutation/pipeline_NGS/playing_around/PSMC/{}/plot_results.py \n'.format(sp_target))
    file.close()
    ## Write top of python scipt
    file = open('{}/plot_results.py'.format(sp_target),'a')
    file.write('PSMC_RESULTS = \"{}.psmc\" \n'.format(sample))
    file.write('MUTATION_RATE = {} \n'.format(mu))
    file.write('GENERAITON_TIME = {} \n'.format(gt))
    file.close()
    # Catenate the rest of the file
    cat = "cat plot_results_all.py >>{}/plot_results.py".format(sp_target, sp_target)
    subprocess.call(cat, shell=True)
    # Run
    run = "sbatch -A MammalianMutation -o {}/psmc_plot.out {}/psmc_plot.sh".format(sp_target, sp_target)
    subprocess.call(run, shell=True)



