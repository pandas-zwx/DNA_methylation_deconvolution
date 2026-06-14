#!/bin/bash
#SBATCH -J MeDen_MK
#SBATCH -p rna             
#SBATCH -N 1                 
#SBATCH --mem=150G                 
#SBATCH --cpus-per-task=30    
#SBATCH -t 3-00:00:00        
#SBATCH -o MeDen_MK.out      
#SBATCH --mail-type=ALL
#SBATCH --mail-user=24211510018@m.fudan.edu.cn


set -euo pipefail

input_dir="./"
output_dir="./"
mktemp_dir=$(mktemp -d)

CpG_clean="/home/u24211510018/workspace/Atlas_WGBS/hg38_CpG_dyad_clean.sorted.bed"
block_file="/home/u24211510018/workspace/Atlas_WGBS/hg38_500bp_blocks_with_CpG_dyad.bed"

# Define the cleanup function for temporary files
cleanup() {
    rm -rf "$mktemp_dir"
}

# Register the cleanup function to run on EXIT, INT, or TERM signals
trap cleanup EXIT INT TERM

for txt_file in "$input_dir"/*txt; do
  [ -e "$txt_file" ] || continue

  base="$(basename "$txt_file" .txt)"
  output_file="${output_dir}/${base}.bed"
  temp_file="${mktemp_dir}/${base}.temp.bed"

  # 清空临时文件
  > "$temp_file"

  # 读取txt文件，处理每一行
  while IFS=$'\t' read -r chr start_1_based me_read total_read; do
    # 将1-based转换为0-based，start减1，end是start+1
    start=$((start_1_based - 1))
    end=$((start_1_based + 1))

    # 输出 BED 格式的五列数据：chr, start, end, me_read, total_read
    echo -e "$chr\t$start\t$end\t$me_read\t$total_read" >> "$temp_file"
  done < "$txt_file"

  # 排序临时文件并保存到最终输出文件
  sort -k1,1 -k2,2n "$temp_file" > "$output_file"

done






set -euo pipefail

input_dir="./"
output_dir="./"

mktemp_dir="./temp"

mkdir -p "$mktemp_dir"

CpG_clean="/home/u24211510018/workspace/Atlas_WGBS/hg38_CpG_dyad_clean.sorted.bed"
block_file="/home/u24211510018/workspace/Atlas_WGBS/hg38_500bp_blocks_with_CpG_dyad.bed"

# Define the cleanup function for temporary files
cleanup() {
    # 这里不删除任何临时文件，以便调试
    echo "临时文件保持不变，方便后续调试"
}

# Register the cleanup function to run on EXIT, INT, or TERM signals
trap cleanup EXIT INT TERM

# Step 1: Sort the CpG file
sort -k1,1 -k2,2n "$CpG_clean" > "${mktemp_dir}/CpG_clean.sorted.bed"

base_name="MK_D3"
output_file="${output_dir}/${base_name}.CpG_dyad.bed"

# Temporary files to store results
tmp_mod_sum="${mktemp_dir}/mod_sum.bed"
tmp_depth_sum="${mktemp_dir}/depth_sum.bed"
tmp_out_raw="${mktemp_dir}/out_raw.bed"

# Step 2: Combine all bed files in the input directory
echo "[1/6] 汇聚所有 bed 文件中的 CpG 位点..."

# 清空临时文件
> "$tmp_mod_sum"
> "$tmp_depth_sum"

# Step 3: 合并所有 bed 文件的内容
cat "$input_dir"/*bed > "${mktemp_dir}/merged_temp.bed"

# Step 4: 对合并后的文件进行排序
sort -k1,1 -k2,2n "${mktemp_dir}/merged_temp.bed" > "${mktemp_dir}/merged_temp_sorted.bed"

# Step 5: 计算每个 CpG dyad 的甲基化读数总和
echo "[4/6] 计算每个 CpG dyad 的甲基化读数总和..."
bedtools map \
  -a "${mktemp_dir}/CpG_clean.sorted.bed" \
  -b "${mktemp_dir}/merged_temp_sorted.bed" \
  -c 4 \
  -f 1.0 \
  -F 1.0 \
  -o sum \
  > "$tmp_mod_sum"

head -10 "$tmp_mod_sum"

# Step 6: 计算每个 CpG dyad 的 depth 总和
echo "[5/6] 计算每个 CpG dyad 的 depth 总和..."
bedtools map \
  -a "${mktemp_dir}/CpG_clean.sorted.bed" \
  -b "${mktemp_dir}/merged_temp_sorted.bed" \
  -c 5 \
  -f 1.0 \
  -F 1.0 \
  -o sum \
  > "$tmp_depth_sum"

head -10 "$tmp_depth_sum"

# Step 7: 计算每个 CpG dyad 的 beta 值
echo "[6/6] 计算每个 CpG dyad 的 beta 值..."
awk 'BEGIN{OFS="\t"}
FNR==NR {
    modsum[$1"\t"$2"\t"$3] = $5
    next
}
{
    key = $1"\t"$2"\t"$3
    m = modsum[key]
    d = $5

    if (m=="." || m=="") m=-1
    if (d=="." || d=="") d=-1

    if (d <= 0 || m < 0) {
        print $1, $2, $3, -1
    } else {
        printf "%s\t%s\t%s\t%.3f\n", $1, $2, $3, m/d
    }
}' "$tmp_mod_sum" "$tmp_depth_sum" > "$tmp_out_raw"

# Step 8: 输出结果，筛选出 beta >= 0 的结果
awk '$4 >= 0 {print $1, $2, $3, $4}' OFS='\t' "$tmp_out_raw" > "$output_file"

echo "完成: $output_file"



set -euo pipefail

bed_file=$(find . -name "*CpG_dyad.bed")

block_file="/home/u24211510018/workspace/Atlas_WGBS/hg38_500bp_blocks_with_CpG_dyad.bed"

base="$(basename "$bed_file" .bed)"
methylation_density_file="${base}.methylation_density.bed"

tmpdir="$(mktemp -d)"
trap 'rm -rf "$tmpdir"' EXIT

tmp_clean_sorted="$tmpdir/${base}.sorted.bed"
tmp_block_sorted="$tmpdir/$(basename "$block_file").sorted.bed"
tmp_block_sum="$tmpdir/${base}.block_sum.bed"
tmp_density_raw="$tmpdir/${base}.methylation_density.raw.bed"

echo "[1/3] 排序 clean bed..."
sort -k1,1 -k2,2n "$bed_file" > "$tmp_clean_sorted"

echo "[2/3] 排序 block bed..."
sort -k1,1 -k2,2n "$block_file" > "$tmp_block_sorted"

echo "[3/3] 计算每个 500bp block 的 methylation density..."
bedtools map \
  -a "$tmp_block_sorted" \
  -b "$tmp_clean_sorted" \
  -F 1.0 \
  -c 4 \
  -o sum \
  > "$tmp_block_sum"

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

awk '$4 >= 0 {print $1, $2, $3, $4}' OFS='\t' "$tmp_density_raw" > "$methylation_density_file"

echo "完成"
echo "methylation density: $methylation_density_file"