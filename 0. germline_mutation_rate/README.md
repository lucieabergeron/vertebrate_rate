The germline mutation rate analysis from raw sequenced to DNMs detection uses the same method than https://github.com/lucieabergeron/germline_mutation_rate.

Additional filters on the false positive calls was conducted on the candidate DNMs.
After runing the 3 first filtering step of the DNMs detection in lucieabergeron/germline_mutation_rate (0.callability.py, 1.negative_rate.py, and .nb_de_novo.py),
the following script were used:
- python 0.rate_samtools.py --> for all samples, produce variant calling with samtools mpileup and bcftools call (samtools and bcftools version 1.2)
- python 1.rate_samtools_denovo.py --> re calculate the number of DNMs (without the appearant FPs)
- python 2.rate_samtools_rate.py --> calculate a mutation rate as: DNM_corrected/(2xCGx(1-FNR))
- Rscript yearly_rate.r --> calculate a yearly rate per samples and then average per sample

