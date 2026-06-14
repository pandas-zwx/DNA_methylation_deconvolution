import os
import re
import glob
import pandas as pd

os.chdir("/home/u24211510018/workspace/Atlas_WGBS/tissue_cell_type_methylation_density")

file_list = sorted(glob.glob("*.hg38.methylation_density.bed"))

series_list = []

for f in file_list:
    # 读入4列：chr, start, end, methylation_density
    df = pd.read_csv(
        f,
        sep="\t",
        header=None,
        names=["chr", "start", "end", "value"]
    )

    # 构造行索引 chr1_10000_10500
    region_id = (
        df["chr"].astype(str) + "_" +
        df["start"].astype(str) + "_" +
        df["end"].astype(str)
    )

    # 提取细胞类型
    # GSM5652185_Kidney-Glomerular-Endothel-Z00000443.hg38.methylation_density.bed
    # -> Kidney-Glomerular-Endothel
    m = re.match(r"^GSM\d+_(.*)-Z[0-9A-Z]+\.hg38\.methylation_density\.bed$", f)
    if m is None:
        raise ValueError(f"文件名格式不符合预期: {f}")
    cell_type = m.group(1)

    # 做成一个 Series，索引是 region_id，name 是列名
    s = pd.Series(df["value"].values, index=region_id, name=cell_type)
    series_list.append(s)

# 按索引外连接，缺失自动补 NaN
combined_df = pd.concat(series_list, axis=1, join="outer")

# 重复列名改成 .1 .2 .3 ...
counts = {}
new_cols = []
for col in combined_df.columns:
    if col not in counts:
        counts[col] = 0
        new_cols.append(col)
    else:
        counts[col] += 1
        new_cols.append(f"{col}.{counts[col]}")

combined_df.columns = new_cols

# 看看结果
print(combined_df.shape)
print(combined_df.head())

# 保存为Python 保存成parquet
combined_df.to_parquet("../R_script_data/combined_methylation_density_matrix.parquet")

# combined_df = pd.read_parquet("../R_script_data/combined_methylation_density_matrix.parquet")
