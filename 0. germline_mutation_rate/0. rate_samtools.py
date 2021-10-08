# -*- coding: utf-8 -*-
"""
This script calls variants with mpilup for each 
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
species = pd.read_csv('sp_here.txt', sep=' ', index_col=None, header=None)
pedigree = pd.read_csv('pedigree.txt', sep='\t', index_col=None, header=None)
genome = pd.read_csv('ref.txt', sep='\t', index_col=None, header=None)

##################################################
# What you run  ##################################
##################################################
# For each chromosome one function:
for sp in range(0,len(species)):
    sp_target=species.loc[sp,0]
    genome_target = genome.loc[genome[0] == sp_target][1].to_list()[0]
    # Create directory
    mkdir = "mkdir {}".format(sp_target)
    subprocess.call(mkdir, shell=True)
    # Directories:
    direct_denovo="/home/lucie/MammalianMutation/faststorage/{}/de_novo_mutation/".format(sp_target)
    direct_bam="/home/lucie/MammalianMutation/faststorage/{}/bam_files/".format(sp_target)
    direct_ref="/home/lucie/MammalianMutation/faststorage/{}/ref_fasta/".format(sp_target)
    direct = "/home/lucie/MammalianMutation/faststorage/{}/".format(sp_target)
    # Dictionary:
    pedigree_sp=pd.read_csv('{}/pedigree.ped'.format(direct), sep=' ', index_col=None, header=None)
    for sample in range(0,len(pedigree_sp)):
        fa = pedigree_sp.iloc[sample,4]
        mo = pedigree_sp.iloc[sample,6]
        off = pedigree_sp.iloc[sample,2]
        denovo_to_check=pd.read_csv('{}/data_denovo_{}.tab'.format(direct_denovo, off), sep='\t', index_col=None)
        file = open('{}/{}_samtools.sh'.format(sp_target, off),'w')
        file.write('#!/bin/bash \n')
        file.write('#SBATCH --partition express,normal,short \n')
        file.write('#SBATCH --mem 10G \n')
        file.write('#SBATCH -c 1 \n')
        file.write('#SBATCH --time=10:00:00 \n')
        for line in range(0,len(denovo_to_check)):
            chrom=denovo_to_check.iloc[line,0]
            pos=denovo_to_check.iloc[line,1]
            file.write('/com/extra/samtools/1.2/bin/samtools mpileup -ugf {}{} -r {}:{}-{} {}{}_sorted.merged.addg.uniq.rmdup.bam {}{}_sorted.merged.addg.uniq.rmdup.bam {}{}_sorted.merged.addg.uniq.rmdup.bam | /com/extra/bcftools/1.2/bin/bcftools call -m | tail -1 | cut -f1,2,4,5,10,11,12 >> /home/lucie/MammalianMutation/pipeline_NGS/playing_around/re_run_denovo/{}/{}_samtools.txt \n'.format(direct_ref, genome_target, chrom, pos, pos, direct_bam, fa, direct_bam, mo, direct_bam, off, sp_target, off))
        file.close()
        sub_cmd = "sbatch -A MammalianMutation -o {}/{}_samtools.out {}/{}_samtools.sh".format(sp_target, off, sp_target, off)
        subprocess.call(sub_cmd, shell=True)



