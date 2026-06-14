#!/bin/bash
#SBATCH -J MeDen
#SBATCH -p rna             
#SBATCH -N 1                 
#SBATCH --mem=150G                 
#SBATCH --cpus-per-task=30    
#SBATCH -t 3-00:00:00        
#SBATCH -o MeDen.out      
#SBATCH --mail-type=ALL
#SBATCH --mail-user=24211510018@m.fudan.edu.cn


set -euo pipefail

input_dir="./"
output_dir="../../tissue_cell_type_methylation_density"

CpG_clean="/home/u24211510018/workspace/Atlas_WGBS/hg38_CpG_dyad_clean.sorted.bed"
block_file="/home/u24211510018/workspace/Atlas_WGBS/hg38_500bp_blocks_with_CpG_dyad.bed"


mkdir -p "$output_dir"

# ========= 全局 tmp =========
tmpdir_global="$(mktemp -d)"
trap 'rm -rf "$tmpdir_global"' EXIT

# ========= 只排序一次 =========
tmp_CpG_sorted="$tmpdir_global/CpG_clean.sorted.bed"
tmp_block_sorted="$tmpdir_global/block.sorted.bed"

echo "准备全局排序文件..."

sort -k1,1 -k2,2n "$CpG_clean" > "$tmp_CpG_sorted"
sort -k1,1 -k2,2n "$block_file" > "$tmp_block_sorted"

echo "全局排序完成"

# ============================

for bed_file in "$input_dir"/*hg38.sorted.bed; do
  [ -e "$bed_file" ] || continue

  base="$(basename "$bed_file" .sorted.bed)"

  methylation_density_file="${output_dir}/${base}.methylation_density.bed"

  tmpdir="$(mktemp -d)"

  tmp_bed_sorted="$tmpdir/${base}.sorted.bed"
  tmp_clean_sorted="$tmpdir/${base}.CpG_clean.sorted.bed"
  tmp_block_sum="$tmpdir/${base}.block_sum.bed"
  tmp_density_raw="$tmpdir/${base}.density.raw.bed"

  echo "==== 处理: $base ===="

  # ====== 1. sort bed_file（关键！）======
  echo "[1/5] 排序输入 bed..."
  sort -k1,1 -k2,2n "$bed_file" > "$tmp_bed_sorted"

  # ====== 2. CpG 精确匹配 ======
  echo "[2/5] CpG 精确匹配 + 过滤无效值..."
  bedtools map \
    -a "$tmp_CpG_sorted" \
    -b "$tmp_bed_sorted" \
    -c 4 \
    -o distinct \
    -null -1 \
    -f 1.0 \
    -F 1.0 \
  | awk '$5 >= 0 {print $1, $2, $3, $5}' OFS='\t' \
  > "$tmpdir/clean.raw.bed"

  # ====== 3. sort clean ======
  echo "[3/5] 排序 clean..."
  sort -k1,1 -k2,2n "$tmpdir/clean.raw.bed" > "$tmp_clean_sorted"

  # ====== 4. block map ======
  echo "[4/5] 计算 block methylation..."
  bedtools map \
    -a "$tmp_block_sorted" \
    -b "$tmp_clean_sorted" \
    -F 1.0 \
    -c 4 \
    -o sum \
    > "$tmp_block_sum"

  # ====== 5. 计算 density ======
  echo "[5/5] 计算 density..."
  awk 'BEGIN{OFS="\t"}
  {
      chrom=$1
      start=$2
      end=$3
      cpg_count=$5
      sum=$6

      if (sum=="." || sum=="") sum=-1

      if (cpg_count==0) {
          density=0
      } else {
          density=(sum / cpg_count) * 100
      }

      printf "%s\t%s\t%s\t%.3f\n", chrom, start, end, density
  }' "$tmp_block_sum" > "$tmp_density_raw"

  awk '$4 >= 0 {print $1, $2, $3, $4}' OFS='\t' \
    "$tmp_density_raw" > "$methylation_density_file"

  echo "完成: $methylation_density_file"

  rm -rf "$tmpdir"
done
