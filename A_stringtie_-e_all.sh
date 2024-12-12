#!/bin/bash
#SBATCH --cpus-per-task 8
#SBATCH --mem 32G
#SBATCH --time 3:00:00
#SBATCH --output finch_ST.out
#SBATCH --error finch_ST.err
#SBATCH --job-name finch_ST

module load stringtie/1.3.4d

cd /scratch/cchen635/GW_finch_color/Hisat2_opt

#for i in {553,557,561,555,559,563,649};
#for i in {554,558,562};
#for i in {651788,652791,653792,654787};
for i in {553,554,555,556,557,558,559,560,561,562,563,564,649,651788,652791,653792,654787,655789,656790};

do

stringtie Chuong${i}_ht2_opt.bam -e -p 8 \
-A /scratch/cchen635/GW_finch_color/String_opt/Chuong${i}_Gene_abundances.tab \
-G /home1/cchen635/genome_ref/bTaeGut1.4.pri/GCF_003957565.2_bTaeGut1.4.pri_genomic_rm_gene.gtf \
-o /scratch/cchen635/GW_finch_color/String_opt/Chuong${i}_ST_e_opt.gtf

done

