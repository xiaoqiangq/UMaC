# *uMaC*
uMaC is a new computational method that unifies two orthogonal pieces of information, i.e. methylation and CNAs, derived from whole-genome bisulfite sequencing (WGBS) data to estimate enriched ctDNA fragments' tumor fraction and CNA status from cell-free DNA (cfDNA).

## uMaC Wiki Page
**For more details on usage and outputs, please visit the [GitHub Wiki page for uMaC]
(https://github.com/xiaoqiangq/uMaC/wiki)**

## Description
uMaC is a two-component framework, in which the first component models methylation changes, and the second one models CNAs, sequentially. Briefly, the methylation-component implements a Bayesian model to enrich ctDNA fragments from cfDNA data, based on distinct haplotypic methylation patterns in tumor vs. normal cfDNA; the CNA-component implements a hidden Markov model (HMM) to detect CNAs based on sequencing depth profiles in data with enriched ctDNA fragments.

The Bayes model and HMM model are described in: A computational framework to unify orthogonal information in DNA methylation and copy number aberrations in cell-free DNA for cancer detection. 

**1. Methylation component**

    java -Xmx4g -jar uMaC-1.0.jar extract -I HCS-1070.dedupped.bam -O HCS-1070.bed 

**2. CNA component**

    Rscript R/uMaC.R --id HCS-1070  --BED HCS-1070.bed  --outDir normalization

## Contacts
If you have any questions or feedback, please contact us at:  
**Email:** <***@vanderbilt.edu>  
**Google Group:** <https://groups.google.com/****>

## Acknowledgements
uMaC is developed and maintained by Qiang Wei and Bingshan Li. 

## Software License
uMaC

