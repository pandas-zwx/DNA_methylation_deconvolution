#!/bin/bash

# 输入文件
chrm_size="hg38.chrom.sizes"

# 输出文件
autochrm="hg38.autosomes.sizes"
block_file="hg38_500bp_blocks.bed"
four_col_block_file="hg38_500bp_four_col_blocks.bed"

# 日志文件
log="generate_blocks.log"

# 清空或初始化日志文件
> "$log"

# 1. 提取常染色体
grep -E '^chr[0-9]+[[:space:]]' "$chrm_size" > "$autochrm" 2>> "$log"

# 2. 生成500bp窗口
bedtools makewindows -g "$autochrm" -w 500 -s 500 > "$block_file" 2>> "$log"

# 3. 添加block编号
awk '{print $1"\t"$2"\t"$3"\tblock_"NR}' "$block_file" > "$four_col_block_file" 2>> "$log"

echo "处理完成！"
echo "输出文件: $autochrm, $block_file, $four_col_block_file"
echo "日志文件: $log"