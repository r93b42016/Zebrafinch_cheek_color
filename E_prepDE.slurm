#!/bin/bash
#SBATCH --cpus-per-task 2
#SBATCH --mem 16G
#SBATCH --time 0:20:00
#SBATCH --output prepDE.out
#SBATCH --error prepDE.err
#SBATCH --job-name prepDE

module load python/3.9.2

####NOTE: no any space line in the sample list file!!!!!

#python prepDE.py3 -c -i sample_lst_2nd_refUCSC -g gene_count_stringtie_prepDE_2nd_-c_refUCSC.csv -t transcript_count_stringtie_prepDE_2nd_-c_refUCSC.csv
python prepDE.py3 -i sample_bTaeGut1.4.pri -g gene_count_stringtie_-e_bTaeGut1.4.pri_ZFcolor.csv -t transcript_count_stringtie_-e_bTaeGut1.4.pri_ZFcolor.csv

