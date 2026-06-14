#!/bin/bash
#SBATCH -J ToBed   
#SBATCH -p rna             
#SBATCH -N 1                 
#SBATCH --mem=100G                 
#SBATCH --cpus-per-task=20    
#SBATCH -t 3-00:00:00        
#SBATCH -o ToBed.out      
#SBATCH --mail-type=ALL
#SBATCH --mail-user=24211510018@m.fudan.edu.cn


set -eo pipefail

source ~/.bashrc
conda activate ucsc_tools

input_dir="./"
output_dir="../all_bed"
fail_file="../fail.txt"

mkdir -p "$output_dir"
: > "$fail_file"

tmpdir="$(mktemp -d)"
trap 'rm -rf "$tmpdir"' EXIT

for i in "$input_dir"/*hg38.bigwig; do
    [ -e "$i" ] || continue

    filename="$(basename "$i" .bigwig)"
    tmp_bedgraph="$(mktemp "$tmpdir/${filename}.XXXXXX.bedGraph")"
    out_file="${output_dir}/${filename}.sorted.bed"

    echo "正在处理: $i"

    if bigWigToBedGraph "$i" "$tmp_bedgraph"; then
        if sort -k1,1 -k2,2n "$tmp_bedgraph" > "$out_file"; then
            echo "已处理: $i -> $out_file"
        else
            echo "$i" >> "$fail_file"
            echo "排序失败，已跳过: $i"
            rm -f "$out_file"
        fi
    else
        echo "$i" >> "$fail_file"
        echo "转换失败，已跳过: $i"
    fi

    rm -f "$tmp_bedgraph"
done
