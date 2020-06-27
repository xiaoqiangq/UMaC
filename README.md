# *uMaC*
uMaC is a new computational method that unifies two orthogonal pieces of information, i.e. methylation and CNAs, derived from whole-genome bisulfite sequencing (WGBS) data to achieve improved detection accuracy. It implements a Bayes model to enrich circulating tumor DNA (ctDNA) from WGBS data based on hypo-methylation haplotypes, and subsequently, models copy number aberrations (CNAs) in the enriched ctDNA data for cancer detection. 



1. java -Xmx4g -jar uMaC-1.0.jar extract -I HCS-1070.dedupped.bam -O HCS-1070.bed 


2. Rscript R/MethyCNA.R --id HCS-1070  --BED HCS-1070.bed  --outDir normalization

