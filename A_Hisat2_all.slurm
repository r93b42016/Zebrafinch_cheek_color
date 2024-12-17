#!/bin/bash
#SBATCH --cpus-per-task 8
#SBATCH --mem 96G
#SBATCH --time 12:00:00
#SBATCH --output finch_color_HT2.out
#SBATCH --error finch_color_HT2.err
#SBATCH --job-name finch_color_HT2

module load gcc/11.3.0 hisat2/2.2.0 samtools/1.17

cd /scratch1/cchen635/GW_finchcolor/GW_trim_reads

#selected libraries
for i in {555,557,558,559,560,561,562,564,649,651788,652791,653792,654787};

do

hisat2 -p 8 -q -x /project/cmchuong_688/genome_ref/bTaeGut1.4.pri/Hisat2_index/bTaeGut1.4.pri_genomic \
 -U Chuong${i}_merged.trimmed.fastq.gz \
| samtools view -u - \
| samtools sort -o /scratch1/cchen635/GW_finchcolor/GW_Hisat2_opt/Chuong${i}_Hisat2_opt.bTaeGut1.4.pri.bam -

done

