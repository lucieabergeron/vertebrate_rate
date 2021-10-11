"""
This script found the parents of origin of all SNPs for each species.

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

# The functions:
def Poo(name, chrom, species, father, mother, offspring):
    file = open('{}/find_poo/{}_{}_POO.sh'.format(species, name, chrom),'w')
    file.write('#!/bin/bash \n')
    file.write('#SBATCH --partition normal \n')
    file.write('#SBATCH --mem 10G \n')
    file.write('#SBATCH -c 2 \n')
    file.write('#SBATCH --time=78:00:00 \n')
    file.write('## \n')
    file.write('source activate pooha \n')
##    file.write('python /home/lucie/MammalianMutation/pipeline_NGS/playing_around/parental_origin/POOHA.py /home/lucie/MammalianMutation/faststorage/{}/vcf_files/genotype_genomicDBI_{}.g.vcf {} {} {} /home/lucie/MammalianMutation/faststorage/{}/bam_files/BACKUP/{}_sorted.merged.addg.uniq.rmdup.bam > /home/lucie/MammalianMutation/pipeline_NGS/playing_around/parental_origin/{}/find_poo/{}_{}_POO.txt'.format(species, chrom, father, mother, offspring, species, offspring, species, name, chrom))
##    file.write('python /home/lucie/MammalianMutation/pipeline_NGS/playing_around/parental_origin/POOHA.py /home/lucie/MammalianMutation/faststorage/{}/vcf_files/inter_vcf/genotype_genomicDBI_{}.g.vcf {} {} {} /home/lucie/MammalianMutation/faststorage/{}/bam_files/BACKUP/{}_sorted.merged.addg.uniq.rmdup.bam > /home/lucie/MammalianMutation/pipeline_NGS/playing_around/parental_origin/{}/find_poo/{}_{}_POO.txt'.format(species, chrom, father, mother, offspring, species, offspring, species, name, chrom))
    file.write('python /home/lucie/MammalianMutation/pipeline_NGS/playing_around/parental_origin/POOHA.py /home/lucie/MammalianMutation/faststorage/{}/vcf_files/inter_vcf/genotype_genomicDBI_{}.g.vcf {} {} {} /home/lucie/MammalianMutation/faststorage/{}/bam_files/{}_sorted.merged.addg.uniq.rmdup.bam > /home/lucie/MammalianMutation/pipeline_NGS/playing_around/parental_origin/{}/find_poo/{}_{}_POO.txt'.format(species, chrom, father, mother, offspring, species, offspring, species, name, chrom))
    file.write('\n')
    file.close()
    ##"""Submit the .sh to the server"""
    sub_cmd = "sbatch -A MammalianMutation -o {}/find_poo/{}_{}_POO.out {}/find_poo/{}_{}_POO.sh".format(species, name, chrom, species, name, chrom)
    subprocess.call(sub_cmd, shell=True)


##################################################
# What you run  ##################################
##################################################
# For each species:
for sp in range(0,len(species)):
    sp_target=species.loc[sp,0]
    print(sp_target)
    # Create directories
    mkdir = "mkdir {}".format(sp_target)
    subprocess.call(mkdir, shell=True)
    mkdir = "mkdir {}/find_poo".format(sp_target)
    subprocess.call(mkdir, shell=True)
    # Import scaffold list
    scaff = pd.read_csv('../../../faststorage/{}/chromosomes.txt'.format(sp_target), sep=' ', index_col=None, header=None)
    scaff_list = scaff[1]
    # Import the name and number of offspring
    samples = pd.read_csv('../../../faststorage/{}/depth.txt'.format(sp_target), sep=' ', index_col=None, header=None)
    # Import the pedigree
    pedigree = pd.read_csv('../../../faststorage/{}/pedigree.ped'.format(sp_target), sep=' ', index_col=None, header=None)
    pedigree_better = pedigree.dropna(axis=1, how='all')
    pedigree_final = pedigree_better.T.reset_index(drop=True).T
    # Create a file for each sample, each scaffold or chromosomes
    for nb_off in range(0,len(samples)):
        sa_target=samples.loc[nb_off,0]
        for sc in scaff_list:
            sc_target=sc
            fa0 = pedigree_final.loc[pedigree_final[1] == sa_target][2]
            fa1 = fa0.to_string(index=False)
            fa = fa1.strip()
            mo0 = pedigree_final.loc[pedigree_final[1] == sa_target][3]
            mo1 = mo0.to_string(index=False)
            mo = mo1.strip()
            Poo(name='{}_{}'.format(sp_target, sa_target), chrom=sc_target, species=sp_target ,father=fa, mother=mo, offspring=sa_target)


