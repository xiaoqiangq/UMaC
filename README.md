# *uMaC*
uMaC is a new computational method that unifies two orthogonal pieces of information, i.e. methylation and CNAs, derived from whole-genome bisulfite sequencing (WGBS) data to estimating enriched ctDNA fragments' tumor fraction from cell-free DNA (cfDNA).

## Description
uMaC is a two-component framework, in which the first component models methylation changes, and the second one models CNAs, sequentially. Briefly, the methylation-component implements a Bayesian model to enrich ctDNA fragments from cfDNA data, based on distinct haplotypic methylation patterns in tumor vs. normal cfDNA; the CNA-component implements a hidden Markov model (HMM) to detect CNAs based on sequencing depth profiles in data with enriched ctDNA fragments.

**1. Methylation component**

    java -Xmx4g -jar uMaC-1.0.jar extract -I HCS-1070.dedupped.bam -O HCS-1070.bed 

**2. CNA component**

    Rscript R/MethyCNA.R --id HCS-1070  --BED HCS-1070.bed  --outDir normalization

## Contacts
If you have any questions or feedback, please contact us at:  
**Email:** <***@vanderbilt.edu>  
**Google Group:** <https://groups.google.com/****>

## Acknowledgements
uMaC is developed and maintained by Qiang Wei and Bingshan Li. 

## Software License
uMaC

