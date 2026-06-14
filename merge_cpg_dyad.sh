#!/usr/bin/env bash
set -euo pipefail

infile="${1:-hg38_CpG.bed}"
outfile="${2:-hg38_CpG_dyad.bed}"
unpaired="${3:-hg38_CpG_unpaired.bed}"

awk '
BEGIN {
    OFS = "\t"
    has_prev = 0
}

# 跳过注释和表头
/^#/ { next }
$1 == "CHROM" || $1 == "#CHROM" { next }

{
    chrom   = $1
    start   = $2
    end     = $3
    strand  = $4
    sctx    = $5
    ctx     = $6

    # 只处理 CpG
    if (ctx != "CpG") {
        print $0 >> unpaired_file
        next
    }

    # 当前没有缓存上一行，就先存起来
    if (!has_prev) {
        prev_chrom  = chrom
        prev_start  = start
        prev_end    = end
        prev_strand = strand
        prev_sctx   = sctx
        prev_ctx    = ctx
        prev_line   = $0
        has_prev = 1
        next
    }

    # 判断 prev 和当前行能不能组成一个 CpG dyad
    # 规则：
    # 1) 同一条染色体
    # 2) 前一条是 +，当前是 -
    # 3) 前一条区间是 [x, x+1)，当前是 [x+1, x+2)
    if (prev_chrom == chrom &&
        prev_strand == "+" &&
        strand == "-" &&
        prev_end == start) {

        # 合并输出为 dyad: [prev_start, start+1)
        print chrom, prev_start, end >> out_file
        has_prev = 0
    }
    else {
        # 前一条配不上，先记到 unpaired
        print prev_line >> unpaired_file

        # 当前行变成新的待配对记录
        prev_chrom  = chrom
        prev_start  = start
        prev_end    = end
        prev_strand = strand
        prev_sctx   = sctx
        prev_ctx    = ctx
        prev_line   = $0
        has_prev = 1
    }
}

END {
    # 最后一条如果还没配对，也记到 unpaired
    if (has_prev) {
        print prev_line >> unpaired_file
    }
}
' out_file="$outfile" unpaired_file="$unpaired" "$infile"

