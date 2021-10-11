# -*- coding: utf-8 -*-
"""
This script found the parental origin in the de novo mutation.

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
    # Import the name and number of offspring
    samples = pd.read_csv('../../../faststorage/{}/depth.txt'.format(sp_target), sep=' ', index_col=None, header=None)
    # For each sample
    file = open('/home/lucie/MammalianMutation/pipeline_NGS/playing_around/parental_origin/{}/denovo_poo/poo.r'.format(sp_target),'w')
    file.write('nb_sample={} \n'.format(len(samples)))
    file.write('path="/home/lucie/MammalianMutation/pipeline_NGS/playing_around/parental_origin/{}/denovo_poo/" \n'.format(sp_target))
    for nb_off in range(0,len(samples)):
        sa_target=samples.loc[nb_off,0]
        file.write('data_{} = read.csv("/home/lucie/MammalianMutation/pipeline_NGS/playing_around/parental_origin/{}/denovo_poo/{}_{}_denovo_POO.txt", sep =" ", header=FALSE) \n'.format(nb_off+1, sp_target, sp_target, sa_target))
        file.write('data_full_{} = read.csv("/home/lucie/MammalianMutation/pipeline_NGS/playing_around/parental_origin/{}/find_poo/{}_{}_POO.txt", sep =" ", header=FALSE) \n'.format(nb_off+1, sp_target, sp_target, sa_target)) 
    file.close()
    # Write the sh file
    file = open('/home/lucie/MammalianMutation/pipeline_NGS/playing_around/parental_origin/{}/denovo_poo/{}_POO_R.sh'.format(sp_target,sp_target),'w')
    file.write('#!/bin/bash \n')
    file.write('#SBATCH --partition express,normal \n')
    file.write('#SBATCH --mem 6G \n')
    file.write('#SBATCH -c 1 \n')
    file.write('#SBATCH --time=00:30:00 \n')
    file.write('## \n')
    file.write('## Catenate \n')
    file.write('cat /home/lucie/MammalianMutation/pipeline_NGS/playing_around/parental_origin/poo.r >>/home/lucie/MammalianMutation/pipeline_NGS/playing_around/parental_origin/{}/denovo_poo/poo.r'.format(sp_target))
    file.write('\n')
    file.write('## Run R \n')
    file.write('Rscript /home/lucie/MammalianMutation/pipeline_NGS/playing_around/parental_origin/{}/denovo_poo/poo.r'.format(sp_target))
    file.close()
    sub_cmd = "sbatch -A MammalianMutation -o /home/lucie/MammalianMutation/pipeline_NGS/playing_around/parental_origin/{}/denovo_poo/{}_POO_R.out /home/lucie/MammalianMutation/pipeline_NGS/playing_around/parental_origin/{}/denovo_poo/{}_POO_R.sh".format(sp_target,sp_target,sp_target,sp_target)
    subprocess.call(sub_cmd, shell=True)

