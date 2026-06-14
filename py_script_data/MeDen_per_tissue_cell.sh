#!/bin/bash

set -euo pipefail

# =========================
# help function
# =========================

show_help() {
cat << 'EOF'
Usage:
  MeDen_per_tissue_cell.sh CpG_clean block_file input_file output_dir file_prefix

Description:
  Convert methylation BED file into methylation density per genomic block.

Arguments:
  CpG_clean     CpG dyad BED file (sorted or unsorted, will be sorted internally)
  block_file    500bp genomic block file with CpG counts
  input_file    methylation BED file (chr start end beta)
  output_dir    output directory
  file_prefix   output prefix (final file: prefix.methylation_density.bed)

Output:
  output_dir/file_prefix.methylation_density.bed

Example:
  bash MeDen_per_tissue_cell.sh \
    hg38_CpG_dyad_clean.sorted.bed \
    hg38_500bp_blocks_with_CpG_dyad.bed \
    sample.hg38.sorted.bed \
    ./output \
    pancreas_sample1
EOF
}

# =========================
# help trigger
# =========================
if [[ $# -eq 0 || "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    show_help
    exit 0
fi

# =========================
# 参数输入
# =========================
CpG_clean="$1"
block_file="$2"
input_file="$3"
output_dir="$4"
file_prefix="$5"

mkdir -p "$output_dir"

# =========================
# 全局 tmp
# =========================
tmpdir_global="$(mktemp -d)"
trap 'rm -rf "$tmpdir_global"' EXIT

tmp_CpG_sorted="$tmpdir_global/CpG_clean.sorted.bed"
tmp_block_sorted="$tmpdir_global/block.sorted.bed"

echo "准备全局排序文件..."

sort -k1,1 -k2,2n "$CpG_clean" > "$tmp_CpG_sorted"
sort -k1,1 -k2,2n "$block_file" > "$tmp_block_sorted"

echo "全局排序完成"

# =========================
# 单样本处理
# =========================

bed_file="$input_file"
methylation_density_file="${output_dir}/${file_prefix}.methylation_density.bed"

tmpdir="$(mktemp -d)"

tmp_bed_sorted="$tmpdir/input.sorted.bed"
tmp_clean_sorted="$tmpdir/clean.sorted.bed"
tmp_block_sum="$tmpdir/block_sum.bed"
tmp_density_raw="$tmpdir/density.raw.bed"

echo "==== 处理: $file_prefix ===="

echo "[1/5] 排序输入 bed..."
sort -k1,1 -k2,2n "$bed_file" > "$tmp_bed_sorted"

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

echo "[3/5] 排序 clean..."
sort -k1,1 -k2,2n "$tmpdir/clean.raw.bed" > "$tmp_clean_sorted"

echo "[4/5] 计算 block methylation..."
bedtools map \
  -a "$tmp_block_sorted" \
  -b "$tmp_clean_sorted" \
  -F 1.0 \
  -c 4 \
  -o sum \
> "$tmp_block_sum"

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