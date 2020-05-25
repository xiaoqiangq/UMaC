# MethyCNA
1. java -Xmx4g -jar MethyCNA-1.0.jar extract -I HCS-1070.dedupped.bam -O HCS-1070.bed 


2. Rscript R/MethyCNA.R --id HCS-1070  --BED HCS-1070.bed  --outDir normalization

