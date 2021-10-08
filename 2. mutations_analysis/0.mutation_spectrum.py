# -*- coding: utf-8 -*-
"""
This script explore the context of the mutation that are from a C or a G base.

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

##################################################
# What you run  ##################################
##################################################
# For each chromosome one function:
for sp in range(0,len(species)):
    sp_target=species.loc[sp,0]
    # Create directory
    mkdir = "mkdir {}".format(sp_target)
    subprocess.call(mkdir, shell=True)
    # Import the name and number of offspring
    samples = pd.read_csv('../re_run_denovo/{}/denovo.txt'.format(sp_target), sep=' ', index_col=None, header=None)
    # Import the mutation for ALL offspring and catenate
    appended_mut = []
    for nb_off in range(0,len(samples)):
        samples_target = samples.loc[nb_off,0]
        mut_sample = pd.read_csv('../re_run_denovo/{}/data_denovo_{}.tab'.format(sp_target, samples_target), sep='\t', index_col=None)
        mut_sample = mut_sample.drop(columns=['Unnamed: 0'])
        mut_sample.columns = ['CHROM', 'POS', 'TYPE', 'denovo', 'typedenovo', 'REF', 'ALT', 'FILTER', 'offspring.GT', 'offspring.AD', 'offspring.DP', 'offspring.GQ', 'offspring.PL', 'offspring.SAC', 'father.GT', 'father.AD', 'father.DP', 'father.GQ', 'father.PL', 'father.SAC', 'mother.GT', 'mother.AD','mother.DP', 'mother.GQ', 'mother.PL', 'mother.SAC']
        appended_mut.append(mut_sample)
    appended_mut = pd.concat(appended_mut, sort=False)
    # Save it
    pd.DataFrame.to_csv(appended_mut, '{}/all_denovo.tab'.format(sp_target), sep='\t')
#################################################################
    # Get the C or G to something
    awk_cmd1 = "awk '$7==\"C\" && $8==\"T\" || $7==\"G\" && $8==\"A\"' {}/all_denovo.tab | cut -f2,3 | sed 's/\"//g' >>{}/pos_CtoT.txt \n".format(sp_target, sp_target)
    awk_cmd2 = "awk '$7==\"C\" && $8==\"G\" || $7==\"G\" && $8==\"C\"' {}/all_denovo.tab | cut -f2,3 | sed 's/\"//g' >>{}/pos_CtoG.txt \n".format(sp_target, sp_target)
    awk_cmd3 = "awk '$7==\"C\" && $8==\"A\" || $7==\"G\" && $8==\"T\"' {}/all_denovo.tab | cut -f2,3 | sed 's/\"//g' >>{}/pos_CtoA.txt \n".format(sp_target, sp_target)
    pos_cmd1= "less {}/pos_CtoT.txt | while read a b; do echo -e \"$a \\t $(expr $b - 1) \\n$a \\t $b \\n$a \\t $(expr $b + 1)\"; done >>{}/pos_all_CtoT.txt \n".format(sp_target, sp_target)
    pos_cmd2= "less {}/pos_CtoG.txt | while read a b; do echo -e \"$a \\t $(expr $b - 1) \\n$a \\t $b \\n$a \\t $(expr $b + 1)\"; done >>{}/pos_all_CtoG.txt \n".format(sp_target, sp_target)
    pos_cmd3= "less {}/pos_CtoA.txt | while read a b; do echo -e \"$a \\t $(expr $b - 1) \\n$a \\t $b \\n$a \\t $(expr $b + 1)\"; done >>{}/pos_all_CtoA.txt \n".format(sp_target, sp_target)
    subprocess.call(awk_cmd1, shell=True)
    subprocess.call(awk_cmd2, shell=True)
    subprocess.call(awk_cmd3, shell=True)
    subprocess.call(pos_cmd1, shell=True)
    subprocess.call(pos_cmd2, shell=True)
    subprocess.call(pos_cmd3, shell=True)
##################################################################
    # CpG
    get_var1="less {}/pos_all_CtoT.txt | while read a b; do echo \"samtools faidx /home/lucie/MammalianMutation/faststorage/{}/ref_fasta/*.fa '$a:$b-$b' >> /home/lucie/MammalianMutation/pipeline_NGS/playing_around/type_of_mutation/{}/CpG_all_CtoT.txt\" >> {}/CpG_all_CtoT.sh ; done \n".format(sp_target,sp_target,sp_target,sp_target)
    get_var2="less {}/pos_all_CtoG.txt | while read a b; do echo \"samtools faidx /home/lucie/MammalianMutation/faststorage/{}/ref_fasta/*.fa '$a:$b-$b' >> /home/lucie/MammalianMutation/pipeline_NGS/playing_around/type_of_mutation/{}/CpG_all_CtoG.txt\" >> {}/CpG_all_CtoG.sh ; done \n".format(sp_target,sp_target,sp_target,sp_target)
    get_var3="less {}/pos_all_CtoA.txt | while read a b; do echo \"samtools faidx /home/lucie/MammalianMutation/faststorage/{}/ref_fasta/*.fa '$a:$b-$b' >> /home/lucie/MammalianMutation/pipeline_NGS/playing_around/type_of_mutation/{}/CpG_all_CtoA.txt\" >> {}/CpG_all_CtoA.sh ; done \n".format(sp_target,sp_target,sp_target,sp_target)
    subprocess.call(get_var1, shell=True)
    subprocess.call(get_var2, shell=True)
    subprocess.call(get_var3, shell=True)
    # format
    grep1 = "grep -v '>' /home/lucie/MammalianMutation/pipeline_NGS/playing_around/type_of_mutation/{}/CpG_all_CtoT.txt > /home/lucie/MammalianMutation/pipeline_NGS/playing_around/type_of_mutation/{}/CpG_all_CtoT1.txt".format(sp_target,sp_target)
    grep2 = "grep -v '>' /home/lucie/MammalianMutation/pipeline_NGS/playing_around/type_of_mutation/{}/CpG_all_CtoG.txt > /home/lucie/MammalianMutation/pipeline_NGS/playing_around/type_of_mutation/{}/CpG_all_CtoG1.txt".format(sp_target,sp_target)
    grep3 = "grep -v '>' /home/lucie/MammalianMutation/pipeline_NGS/playing_around/type_of_mutation/{}/CpG_all_CtoA.txt > /home/lucie/MammalianMutation/pipeline_NGS/playing_around/type_of_mutation/{}/CpG_all_CtoA1.txt".format(sp_target,sp_target)
    format1="sed 's/[a-z]/\\U&/g' /home/lucie/MammalianMutation/pipeline_NGS/playing_around/type_of_mutation/{}/CpG_all_CtoT1.txt > /home/lucie/MammalianMutation/pipeline_NGS/playing_around/type_of_mutation/{}/CpG_all_CtoT2.txt".format(sp_target,sp_target)
    format2="sed 's/[a-z]/\\U&/g' /home/lucie/MammalianMutation/pipeline_NGS/playing_around/type_of_mutation/{}/CpG_all_CtoG1.txt > /home/lucie/MammalianMutation/pipeline_NGS/playing_around/type_of_mutation/{}/CpG_all_CtoG2.txt".format(sp_target,sp_target)
    format3="sed 's/[a-z]/\\U&/g' /home/lucie/MammalianMutation/pipeline_NGS/playing_around/type_of_mutation/{}/CpG_all_CtoA1.txt >/home/lucie/MammalianMutation/pipeline_NGS/playing_around/type_of_mutation/{}/CpG_all_CtoA2.txt".format(sp_target, sp_target)
    paste1 = "paste /home/lucie/MammalianMutation/pipeline_NGS/playing_around/type_of_mutation/{}/pos_all_CtoT.txt /home/lucie/MammalianMutation/pipeline_NGS/playing_around/type_of_mutation/{}/CpG_all_CtoT2.txt >/home/lucie/MammalianMutation/pipeline_NGS/playing_around/type_of_mutation/{}/CpG_all_CtoT_used.txt".format(sp_target,sp_target,sp_target)
    paste2 = "paste /home/lucie/MammalianMutation/pipeline_NGS/playing_around/type_of_mutation/{}/pos_all_CtoG.txt /home/lucie/MammalianMutation/pipeline_NGS/playing_around/type_of_mutation/{}/CpG_all_CtoG2.txt >/home/lucie/MammalianMutation/pipeline_NGS/playing_around/type_of_mutation/{}/CpG_all_CtoG_used.txt".format(sp_target,sp_target,sp_target)
    paste3 = "paste /home/lucie/MammalianMutation/pipeline_NGS/playing_around/type_of_mutation/{}/pos_all_CtoA.txt /home/lucie/MammalianMutation/pipeline_NGS/playing_around/type_of_mutation/{}/CpG_all_CtoA2.txt >/home/lucie/MammalianMutation/pipeline_NGS/playing_around/type_of_mutation/{}/CpG_all_CtoA_used.txt".format(sp_target,sp_target,sp_target)
    file = open('{}/CpG_all_CtoT.sh'.format(sp_target),'a')
    file.write(grep1)
    file.write('\n')
    file.write(format1)
    file.write('\n')
    file.write(paste1)
    file.write('\n')
    file.close()
    file = open('{}/CpG_all_CtoG.sh'.format(sp_target),'a')
    file.write(grep2)
    file.write('\n')
    file.write(format2)
    file.write('\n')
    file.write(paste2)
    file.write('\n')
    file.close()
    file = open('{}/CpG_all_CtoA.sh'.format(sp_target),'a')
    file.write(grep3)
    file.write('\n')
    file.write(format3)
    file.write('\n')
    file.write(paste3)
    file.write('\n')
    file.close()
    # run the file
    bash_get_var1="bash /home/lucie/MammalianMutation/pipeline_NGS/playing_around/type_of_mutation/{}/CpG_all_CtoT.sh".format(sp_target)
    bash_get_var2="bash /home/lucie/MammalianMutation/pipeline_NGS/playing_around/type_of_mutation/{}/CpG_all_CtoG.sh".format(sp_target)
    bash_get_var3="bash /home/lucie/MammalianMutation/pipeline_NGS/playing_around/type_of_mutation/{}/CpG_all_CtoA.sh".format(sp_target)
    subprocess.call(bash_get_var1, shell=True)
    subprocess.call(bash_get_var2, shell=True)
    subprocess.call(bash_get_var3, shell=True)
##########################################################
    # Copy Rscript and run them
    file = open('{}/CpG.r'.format(sp_target),'w')
    file.write('path=\"/home/lucie/MammalianMutation/pipeline_NGS/playing_around/type_of_mutation/{}/\" \n'.format(sp_target))
    file.close()
    file = open('{}/type.r'.format(sp_target),'w')
    file.write('path=\"/home/lucie/MammalianMutation/pipeline_NGS/playing_around/type_of_mutation/{}/\" \n'.format(sp_target))
    file.close()
    cp1="cat CpG.r >>{}/CpG.r".format(sp_target)
    cp2="cat type.r >>{}/type.r".format(sp_target)
    subprocess.call(cp1, shell=True)
    subprocess.call(cp2, shell=True)
    file = open('{}/Rscript_run.sh'.format(sp_target),'w')
    file.write('#!/bin/bash \n')
    file.write('#SBATCH --mem 10G \n')
    file.write('#SBATCH -c 1 \n')
    file.write('#SBATCH --time=00:01:00 \n')
    file.write('Rscript /home/lucie/MammalianMutation/pipeline_NGS/playing_around/type_of_mutation/{}/CpG.r \n'.format(sp_target))
    file.write('sleep 20 \n')
    file.write('Rscript /home/lucie/MammalianMutation/pipeline_NGS/playing_around/type_of_mutation/{}/type.r \n'.format(sp_target))
    file.close()
    sub_cmd = "sbatch -A MammalianMutation -o {}/Rscript_run.out {}/Rscript_run.sh".format(sp_target,sp_target)
    subprocess.call(sub_cmd, shell=True)

