#!/bin/bash
#PBS -P MST107119
#PBS -W group_list=MST107119
#PBS -q ngs16G
#PBS -N trimm_all
#PBS -l select=1:ncpus=4
#PBS -o ./trim_all_95-02.out
#PBS -e ./trim_all_95-02.err

cd /work1/u6585244/Gway_collaboration/new_finch_reads
mkdir QC
for i in {Chuong_787,Chuong_788,Chuong_789,Chuong_790,Chuong_791,Chuong_792};

do

java -jar /pkg/biology/Trimmomatic/Trimmomatic_v0.36/trimmomatic-0.36.jar SE ${i}_R1.fastq.gz QC/${i}_merged.trimmed.fastq.gz ILLUMINACLIP:TruSeq2and3_SE.fa:2:30:10 LEADING:18 TRAILING:18 SLIDINGWINDOW:4:15 MINLEN:36

done
