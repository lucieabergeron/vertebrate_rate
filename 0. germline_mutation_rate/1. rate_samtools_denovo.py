# -*- coding: utf-8 -*-
"""
This script identifies in the new variant calling the true de novo candidates and the false positive ones.
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

##################################################
# What you run  ##################################
##################################################
# For each chromosome one function:
for sp in range(0,len(species)):
    sp_target=species.loc[sp,0]
    print(sp_target)
    # Directories:
    direct_denovo="/home/lucie/MammalianMutation/faststorage/{}/de_novo_mutation/".format(sp_target)
    direct_here="/home/lucie/MammalianMutation/pipeline_NGS/playing_around/re_run_denovo/{}/".format(sp_target)
    direct = "/home/lucie/MammalianMutation/faststorage/{}/".format(sp_target)
    # Dictionary:
    pedigree_sp=pd.read_csv('{}/pedigree.ped'.format(direct), sep=' ', index_col=None, header=None)
    for sample in range(0,len(pedigree_sp)):
        off = pedigree_sp.iloc[sample,2]
        print(off)
        denovo_to_check=pd.read_csv('{}data_denovo_{}.tab'.format(direct_denovo, off), sep='\t', index_col=None)
        denovo_samtools=pd.read_csv('{}{}_samtools.txt'.format(direct_here, off), sep='\t', index_col=None, header=None)
        idx_pass = []
        idx_nov = []
        idx_parentv = []
        idx_idk = []
        for line in range(0,len(denovo_samtools)):
            if denovo_samtools.iloc[line,4].split(':')[0]=='0/0' and denovo_samtools.iloc[line,5].split(':')[0]=='0/0' and denovo_samtools.iloc[line,6].split(':')[0]=='0/1':
                idx_pass.append(line)
            elif denovo_samtools.iloc[line,4].split(':')[0]=='0/0' and denovo_samtools.iloc[line,5].split(':')[0]=='0/0' and denovo_samtools.iloc[line,6].split(':')[0]=='0/0':
                idx_nov.append(line)
            elif (denovo_samtools.iloc[line,4].split(':')[0]=='0/0' and denovo_samtools.iloc[line,5].split(':')[0]=='0/1' and denovo_samtools.iloc[line,6].split(':')[0]=='0/1') or (denovo_samtools.iloc[line,4].split(':')[0]=='0/1' and denovo_samtools.iloc[line,5].split(':')[0]=='0/0' and denovo_samtools.iloc[line,6].split(':')[0]=='0/1'):
                idx_parentv.append(line)
            else:
                idx_idk.append(line)
        table_pass=denovo_samtools.iloc[idx_pass]
        table_nov=denovo_samtools.iloc[idx_nov]
        table_parentv=denovo_samtools.iloc[idx_parentv]
        table_idk=denovo_samtools.iloc[idx_idk]
        data_new= pd.DataFrame(columns=denovo_to_check.columns)
        data_nov= pd.DataFrame(columns=denovo_to_check.columns)
        data_parentv= pd.DataFrame(columns=denovo_to_check.columns)
        data_idk= pd.DataFrame(columns=denovo_to_check.columns)
        for denovo in range(0,len(denovo_to_check)):
            chrom = denovo_to_check.iloc[denovo,0]
            pos = denovo_to_check.iloc[denovo,1]
            if len(table_pass.loc[(table_pass[0] ==chrom) & (table_pass[1] == pos)])>0:
                data_new=data_new.append(denovo_to_check.iloc[denovo])
                print('pass')
            elif len(table_nov.loc[(table_nov[0] ==chrom) & (table_nov[1] == pos)])>0:
                data_nov=data_nov.append(denovo_to_check.iloc[denovo])
                print('nov')
            elif len(table_parentv.loc[(table_parentv[0] ==chrom) & (table_parentv[1] == pos)])>0:
                data_parentv=data_parentv.append(denovo_to_check.iloc[denovo])
                print('parentv')
            elif len(table_idk.loc[(table_idk[0] ==chrom) & (table_idk[1] == pos)])>0:
                data_idk=data_idk.append(denovo_to_check.iloc[denovo])
                print('idk')
        nb_denovo=len(data_new)
        file = open('{}denovo.txt'.format(direct_here),'a')
        file.write('{} {} \n'.format(off,nb_denovo))
        file.close()
        data_new.to_csv('{}data_denovo_{}.tab'.format(direct_here, off), sep='\t')
        data_nov.to_csv('{}data_denovo_{}_nov.tab'.format(direct_here, off), sep='\t')
        data_parentv.to_csv('{}data_denovo_{}_parentv.tab'.format(direct_here, off), sep='\t')
        data_idk.to_csv('{}data_denovo_{}_idk.tab'.format(direct_here, off), sep='\t')


