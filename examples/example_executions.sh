## update as required
rpsblastdb=/mnt/chaelab/shared/blastdb/Cdd.v3.18/Cdd
rpsblast=/usr/bin/rpsblast+

## example 1: generate gRNA for NB-ARC domain of genes in RPS6 cluster in reference genome
../minorg/main.py --dir ./example1 --prefix example1 \
                  --assembly ./subset_ref_TAIR10.fasta --annotation ./subset_ref_TAIR10.gff \
                  --domain 366375 --rps-db ${rpsblastdb} --rpsblast ${rpsblast} \
                  --cluster-set ./subset_cluster_mapping.txt \
                  --indv ref --cluster RPS6 --minlen 400 --screen-ref

## example 2: generate gRNA for homologue in non-reference genome of gene not in reference genome by extending reference genome using --extend-genome and --extend-cds. Keep files generated if programme crashes.
../minorg/main.py --dir ./example2 --prefix example2 \
                  --reference TAIR10 --reference-set ./athaliana_genomes.txt \
                  --extend-genome ./sample_gene.fasta --extend-cds ./sample_CDS.fasta \
                  --query ./subset_9944.fasta \
                  --gene AT1G10920 --check-recip --minid 95 --minlen 400 --keep-on-crash

## example 3: generate gRNA for orthologues in 3 species
../minorg/main.py --dir ./example3 --prefix example3 \
                  --reference TAIR10,araly2,araha1 --reference-set ./arabidopsis_genomes.txt \
                  --gene AT1G33560,AL1G47950.v2.1,Araha.3012s0003.v1.1 \
                  --indv ref --screen-ref


## export config file to environment
export MINORG_CONFIG=/path/to/sample_config.ini

## example 1 using config file
../minorg/main.py --dir ./example1_config --prefix example1 \
                  --assembly TAIR10 --annotation TAIR10 \
                  --domain 366375 --reference - \
                  --cluster-set arabidopsisnlr \
                  --indv ref --cluster RPS6 --minlen 400 --screen-ref

## example 2 using config file
../minorg/main.py --dir ./example2_config --prefix example2 \
                  --reference TAIR10 --reference-set athaliana \
                  --extend-genome ./sample_gene.fasta --extend-cds ./sample_CDS.fasta \
                  --query ./subset_9944.fasta \
                  --gene AT1G10920 --check-recip --minid 95 --minlen 400 --keep-on-crash

## example 3 using config file
../minorg/main.py --dir ./example3_config --prefix example3 \
                  --reference TAIR10,araly2,araha1 --reference-set arabidopsis \
                  --gene AT1G33560,AL1G47950.v2.1,Araha.3012s0003.v1.1 \
                  --indv ref --screen-ref

