#!/bin/bash
#SBATCH -J find_CpG   
#SBATCH -p rna             
#SBATCH -N 1                 
#SBATCH --mem=100G                 
#SBATCH --cpus-per-task=20    
#SBATCH -t 3-00:00:00        
#SBATCH -o find_CpG.out      
#SBATCH --mail-type=ALL
#SBATCH --mail-user=24211510018@m.fudan.edu.cn



# 生成 hg38_CpG_dyad.bed 文件，包含 hg38 基因组中所有 CpG 位点的位置信息
set -eo pipefail
source ~/.bashrc

conda activate ucsc_tools

seqkit locate --bed -P -i -j 20 -p CG ../TAPS_DEMO/reference_genome_annotation/hg38.fa \
| awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "CpG_"NR}' \
> hg38_CpG_dyad.bed


# because  hg38_CpG_dyad.bed contains some unstandard chromosome, so these codes are used to filter out the unstandard chromosome and rename the CpG sites in a standard way.
# hg38_CpG_dyad_clean.bed contains sex and mitochondial chromosomes
awk 'BEGIN{OFS="\t"}
$1 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y|M)$/ {
    n++
    print $1, $2, $3, "CpG_" n
}' hg38_CpG_dyad.bed > hg38_CpG_dyad_clean.bed



sort -k1,1 -k2,2n hg38_500bp_four_col_blocks.bed > hg38_500bp_four_col_blocks.sorted.bed
sort -k1,1 -k2,2n hg38_CpG_dyad_clean.bed > hg38_CpG_dyad_clean.sorted.bed

# hg38_500bp_four_col_blocks.bed have removed sex and mito chromosomes, so hg38_500bp_blocks_CpG_dyad_counts.bed only contains CpG dyad counts in autosomal chromosomes.
bedtools intersect \
  -a hg38_500bp_four_col_blocks.sorted.bed \
  -b hg38_CpG_dyad_clean.sorted.bed \
  -F 1.0 \
  -c \
> hg38_500bp_blocks_CpG_dyad_counts.bed

awk '$5 > 0' hg38_500bp_blocks_CpG_dyad_counts.bed > hg38_500bp_blocks_with_CpG_dyad.bed