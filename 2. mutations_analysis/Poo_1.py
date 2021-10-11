# -*- coding: utf-8 -*-
"""
This script catenate all chromosome in one file and grep the mutations.

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
species = pd.read_csv('sp.txt', sep=' ', index_col=None, header=None)
# Species maybe import in ../ directory the one filtered

##################################################
# What you run  ##################################
##################################################
# For each species:
for sp in range(0,len(species)):
    sp_target=species.loc[sp,0]
    print(sp_target)
    # Create directories
    mkdir = "mkdir {}/denovo_poo".format(sp_target)
    subprocess.call(mkdir, shell=True)
    # Import scaffold list
    scaff = pd.read_csv('../../../faststorage/{}/chromosomes.txt'.format(sp_target), sep=' ', index_col=None, header=None)
    scaff_list = scaff[1]
    # Import the name and number of offspring
    samples = pd.read_csv('../../../faststorage/{}/depth.txt'.format(sp_target), sep=' ', index_col=None, header=None)
    # For each sample
    for nb_off in range(0,len(samples)):
        sa_target=samples.loc[nb_off,0]
        denovo=pd.read_csv('../re_run_denovo/{}/data_denovo_{}.tab'.format(sp_target, sa_target), sep='\t', index_col=None)
        #Get catenate command
        cat_cmd='cat '
        for sc in scaff_list:
            cat_cmd+='/home/lucie/MammalianMutation/pipeline_NGS/playing_around/parental_origin/{}/find_poo/{}_{}_{}_POO.txt '.format(sp_target, sp_target, sa_target, sc)
        cat_cmd+='> /home/lucie/MammalianMutation/pipeline_NGS/playing_around/parental_origin/{}/find_poo/{}_{}_POO.txt'.format(sp_target, sp_target, sa_target)
        # Write a file with catenate and get de novo commands
        file = open('/home/lucie/MammalianMutation/pipeline_NGS/playing_around/parental_origin/{}/denovo_poo/{}_{}_POO.sh'.format(sp_target,sp_target,sa_target),'w')
        file.write('#!/bin/bash \n')
        file.write('#SBATCH --partition normal \n')
        file.write('#SBATCH --mem 6G \n')
        file.write('#SBATCH -c 1 \n')
        file.write('#SBATCH --time=02:00:00 \n')
        file.write('## \n')
        file.write('## Catenate \n')
        file.write(cat_cmd)
        file.write('\n')
        file.write('## Get denovo \n')
        for mu in range(0,len(denovo)):
            chrom0 = denovo.iloc[[mu],[1]]
            pos0 = denovo.iloc[[mu],[2]]
            chrom1 = chrom0.to_string(index=False, header=False)
            pos1 = pos0.to_string(index=False, header=False)
            chrom = chrom1.strip()
            pos = pos1.strip()
            file.write('grep \'{} {}\' /home/lucie/MammalianMutation/pipeline_NGS/playing_around/parental_origin/{}/find_poo/{}_{}_POO.txt >> /home/lucie/MammalianMutation/pipeline_NGS/playing_around/parental_origin/{}/denovo_poo/{}_{}_denovo_POO.txt \n'.format(chrom, pos, sp_target, sp_target, sa_target, sp_target, sp_target, sa_target))
        file.close()
        sub_cmd = "sbatch -A MammalianMutation -o /home/lucie/MammalianMutation/pipeline_NGS/playing_around/parental_origin/{}/denovo_poo/{}_{}_POO.out /home/lucie/MammalianMutation/pipeline_NGS/playing_around/parental_origin/{}/denovo_poo/{}_{}_POO.sh".format(sp_target,sp_target,sa_target,sp_target,sp_target,sa_target)
        subprocess.call(sub_cmd, shell=True)

