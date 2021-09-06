rpsblastdb=/mnt/chaelab/shared/blastdb/Cdd.v3.18/Cdd

../minorg.py --dir ./example1 --prefix example1 \
             --reference ./subset_ref.fasta --bed ./subset_ref.bed --domain 366375 --rps-db ${rpsblastdb} \
             --cluster-set ./subset_cluster_mapping.txt \
             --indv ref --cluster RPS6 --check-recip --minid 95 --minlen 400

../minorg.py --dir ./example2 --prefix example2 \
             --reference ./subset_ref.fasta --bed ./subset_ref.bed \
             --extend-genome ./sample_gene.fasta --extend-cds ./sample_CDS.fasta \
             --query ./subset_9944.fasta \
             --gene AT1G10920 --check-recip --minid 95 --minlen 400 --keep-on-crash
