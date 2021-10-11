# -*- coding: utf-8 -*-
"""
This script calculate the PSMC for all species.

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
species = pd.read_csv('sp.txt', sep='\t', index_col=None, header=None)

##################################################
# What you run  ##################################
##################################################
# For each chromosome one function:
for sp in range(0,len(species)):
    sp_target=species.loc[sp,0]
    ref_target=species.loc[sp,1]
    sample_target=species.loc[sp,2]
    dp_min = species.loc[sp,3]
    dp_max = species.loc[sp,4]
    # Create directory
    mkdir = "mkdir {}".format(sp_target)
    subprocess.call(mkdir, shell=True)
    file = open('{}/psmc.sh'.format(sp_target),'a')
    file.write('#!/bin/bash \n')
    file.write('#SBATCH --mem 30G \n')
    file.write('#SBATCH -c 3 \n')
    file.write('#SBATCH --time=70:00:00 \n')
    file.write('## get seq \n')
    file.write('/com/extra/samtools/0.1.19/bin/samtools mpileup -C50 -uf /home/lucie/MammalianMutation/faststorage/{}/ref_fasta/{} /home/lucie/MammalianMutation/faststorage/{}/bam_files/BACKUP/{}_sorted.merged.addg.uniq.rmdup.bam | /com/extra/samtools/0.1.19/bin/bcftools view -c - | /com/extra/samtools/0.1.19/bin/vcfutils.pl vcf2fq -d {} -D {} | gzip > {}.fq.gz \n'.format(sp_target, ref_target, sp_target, sample_target, dp_min, dp_max, sample_target))
    file.write('## Start psmc \n')
    file.write('~/bin/psmc/psmc-master/utils/fq2psmcfa -q20 {}.fq.gz > {}.psmcfa \n'.format(sample_target, sample_target))
    file.write('## Second step \n')
    file.write('~/bin/psmc/psmc-master/psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o {}.psmc {}.psmcfa \n'.format(sample_target, sample_target))
    file.write('## 3rd step \n')
    file.write('~/bin/psmc/psmc-master/utils/psmc2history.pl {}.psmc | ~/bin/psmc/psmc-master/utils/history2ms.pl > ms-cmd.sh \n'.format(sample_target))
    file.write('## Plot \n')
    file.write('source activate matplot \n')
    file.write('python plot_results.py \n')
    file.close()

