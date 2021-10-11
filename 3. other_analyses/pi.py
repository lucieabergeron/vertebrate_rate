# -*- coding: utf-8 -*-
"""
This script calculate the nucleotide diversity for the individuals of a species.

Lucie Bergeron
23.08.18
"""
##################################################
# What you need ##################################
##################################################

# Packages:
import subprocess
import os
##from variable import *
import pandas as pd

# The functions:
def pi(sp, ref, directory_pi):
    """Estimate SFS"""
    sfs_cmd1 = "~/bin/angsd/angsd -bam {}bamlist.txt -doSaf 1 -anc ~/MammalianMutation/faststorage/{}/ref_fasta/{} -GL 1 -P 4 -out {}thetasStep1-out".format(directory_pi, sp, ref, directory_pi)
    sfs_cmd2 = "~/bin/angsd/misc/realSFS {}thetasStep1-out.saf.idx -P 4 > {}thetasStep1-out.sfs".format(directory_pi, directory_pi)
    """Thetra per site"""
    thetra_cmd1 = "~/bin/angsd/angsd -bam {}bamlist.txt -out thetasStep1-out -doThetas 1 -doSaf 1 -pest thetasStep1-out.sfs -anc ~/MammalianMutation/faststorage/{}/ref_fasta/{} -GL 1".format(directory_pi, sp, ref, directory_pi)
    thetra_cmd2 = "~/bin/angsd/misc/thetaStat do_stat {}thetasStep1-out.thetas.idx".format(directory_pi)
    """Window estimate"""
    win_cmd = "~/bin/angsd/misc/thetaStat do_stat {}thetasStep1-out.thetas.idx -win 50000 -step 10000 -outnames theta.thetasWin-10kb.gz".format(directory_pi)
    """R command"""
    R_cmd = "Rscript {}r_sum.r".format(directory_pi)
    ##
    file = open('{}pi1.sh'.format(directory_pi),'w')
    file.write('#!/bin/bash \n')
    file.write('#SBATCH --partition normal \n')
    file.write('#SBATCH --mem 250G \n')
    file.write('#SBATCH -c 4 \n')
    file.write('#SBATCH --time=12:00:00 \n')
    file.write('#Step 1: global estimate of the SFS \n')
    file.write(sfs_cmd1)
    file.write('\n')
    file.write(sfs_cmd2)
    file.write('\n')
    file.close()
    file = open('{}pi2.sh'.format(directory_pi),'w')
    file.write('#!/bin/bash \n')
    file.write('#SBATCH --partition normal \n')
    file.write('#SBATCH --mem 250G \n')
    file.write('#SBATCH -c 4 \n')
    file.write('#SBATCH --time=12:00:00 \n')
    file.write('#Step2: Calculate thetas per site \n')
    file.write(thetra_cmd1)
    file.write('\n')
    file.write(thetra_cmd2)
    file.write('\n')
    file.close()
    file = open('{}pi3.sh'.format(directory_pi),'w')
    file.write('#!/bin/bash \n')
    file.write('#SBATCH --partition normal \n')
    file.write('#SBATCH --mem 250G \n')
    file.write('#SBATCH -c 4 \n')
    file.write('#SBATCH --time=12:00:00 \n')
    file.write('#Step3: Estimate over windows \n')
    file.write(win_cmd)
    file.write('\n')
    file.write('#Step 4: Summary R script\n')
    file.write(R_cmd)
    file.write('\n')
    file.close()


##################################################
# What you run  ##################################
##################################################
# For species in the file: create bamfile list, then run the thing
species = pd.read_csv('sp.txt', sep='\t', index_col=None, header=None)
for s in range(65):
  sp = species.iloc[s,0]
  ref = species.iloc[s,2]
  directory_pi= "/home/lucie/MammalianMutation/pipeline_NGS/playing_around/nucleotide_diversity/my_data/{}/".format(sp)
  # Write R script:
  file = open('{}r_sum.r'.format(directory_pi),'w')
  file.write('data=read.csv("{}theta.thetasWin-10kb.gz.pestPG", sep=\'\\t\') \n'.format(directory_pi))
  file.write('write.table(c(mean(na.omit(data$tP/data$nSites)), sd(na.omit(data$tP/data$nSites))), file="{}res.txt",row.names = FALSE, col.names = FALSE)'.format(directory_pi))
  file.close()
  # Write pi.sh
  pi(sp = sp, ref = ref, directory_pi= directory_pi)
