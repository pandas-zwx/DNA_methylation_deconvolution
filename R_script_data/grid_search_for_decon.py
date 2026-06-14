import os
import numpy as np
import pandas as pd
from scipy.optimize import nnls
from sklearn.metrics import mean_squared_error
from itertools import product
import warnings

warnings.filterwarnings('ignore')

# ================================
# 0. 设置工作目录
# ================================
os.chdir("/home/u24211510018/workspace/Atlas_WGBS/tissue_cell_type_methylation_density")
print(f"工作目录: {os.getcwd()}")

# ================================
# 1. 加载数据
# ================================
print("加载数据...")

combined_clean_mean_df = pd.read_parquet("../R_script_data/combined_clean_mean_df.parquet")
ref_df = combined_clean_mean_df.copy()
print(f"原始参考矩阵大小: {ref_df.shape}")

cpg_bed = pd.read_csv(
    "../hg38_500bp_blocks_with_CpG_dyad.bed",
    sep="\t",
    header=None,
    names=["chr", "start", "end", "block_id", "cpg_count"]
)

cpg_bed["region_id"] = (
    cpg_bed["chr"].astype(str) + "_" +
    cpg_bed["start"].astype(str) + "_" +
    cpg_bed["end"].astype(str)
)

cpg_info = cpg_bed[["region_id", "cpg_count"]].copy()
cpg_info = cpg_info[cpg_info["region_id"].isin(ref_df.index)].copy()
print(f"匹配的CpG block数: {cpg_info.shape[0]}")

# ================================
# 2. 构建参考矩阵（增加每个组织是否达标）
# ================================
def build_reference_matrix(ref_df, cpg_info, min_cpg_count, cv_threshold, top_n):
    """构建参考矩阵，同时返回每个组织是否达到 top_n*0.9"""

    # 1. CpG数量过滤
    cpg_pass_blocks = cpg_info.loc[cpg_info["cpg_count"] >= min_cpg_count, "region_id"]
    ref_df_cpg = ref_df.loc[ref_df.index.intersection(cpg_pass_blocks)].copy()

    if ref_df_cpg.empty:
        return None, pd.DataFrame(), {}

    # 2. CV过滤
    block_mean = ref_df_cpg.mean(axis=1)
    block_sd = ref_df_cpg.std(axis=1)
    block_cv = block_sd / block_mean.replace(0, np.nan)

    cv_pass_blocks = block_cv.index[block_cv >= cv_threshold]
    filtered_df = ref_df_cpg.loc[cv_pass_blocks].copy()

    if filtered_df.empty:
        return None, pd.DataFrame(), {}

    # 3. 筛选 top_n marker
    tissues = filtered_df.columns.tolist()
    all_markers = []
    tissue_sufficient_dict = {}  # 记录每个组织是否达到 top_n*0.9

    for tissue in tissues:
        this_values = filtered_df[tissue]
        other_tissues = [x for x in tissues if x != tissue]
        other_df = filtered_df[other_tissues]

        hyper_score = this_values - other_df.max(axis=1)
        hypo_score = other_df.min(axis=1) - this_values

        # hyper 候选
        hyper_mask = hyper_score > 0
        hyper_blocks = pd.DataFrame()
        if hyper_mask.any():
            hyper_blocks = pd.DataFrame({
                "block": filtered_df.index[hyper_mask],
                "tissue": tissue,
                "score": hyper_score[hyper_mask],
                "type": "hyper"
            })
            hyper_blocks["score_abs"] = hyper_blocks["score"].abs()

        # hypo 候选
        hypo_mask = hypo_score > 0
        hypo_blocks = pd.DataFrame()
        if hypo_mask.any():
            hypo_blocks = pd.DataFrame({
                "block": filtered_df.index[hypo_mask],
                "tissue": tissue,
                "score": hypo_score[hypo_mask],
                "type": "hypo"
            })
            hypo_blocks["score_abs"] = hypo_blocks["score"].abs()

        # 合并后排序取前 top_n
        merged_blocks = pd.concat([hyper_blocks, hypo_blocks], ignore_index=True)
        if not merged_blocks.empty:
            merged_blocks = merged_blocks.sort_values("score_abs", ascending=False).head(top_n)
            all_markers.append(merged_blocks)

        # 判断组织是否达到 top_n*0.9
        tissue_sufficient_dict[tissue] = len(merged_blocks) >= (top_n * 0.9)

    if not all_markers:
        return None, pd.DataFrame(), tissue_sufficient_dict

    marker_df = pd.concat(all_markers, ignore_index=True)
    marker_blocks = marker_df["block"].unique()

    # 构建最终参考矩阵
    final_ref_matrix = filtered_df.loc[filtered_df.index.intersection(marker_blocks)].copy()

    return final_ref_matrix, marker_df, tissue_sufficient_dict

# ================================
# 3. 评估函数
# ================================
def evaluate_reference_matrix(ref_matrix, n_simulations=50):
    if ref_matrix is None or ref_matrix.empty:
        return {
            "avg_corr": -1,
            "avg_mse": np.inf,
            "avg_residual": np.inf,
            "n_markers": 0
        }

    n_cell_types = ref_matrix.shape[1]
    n_regions = ref_matrix.shape[0]

    if n_regions < 10 or n_cell_types < 2:
        return {
            "avg_corr": -1,
            "avg_mse": np.inf,
            "avg_residual": np.inf,
            "n_markers": n_regions
        }

    X = ref_matrix.values.T
    all_correlations, all_mse, all_residuals = [], [], []

    np.random.seed(42)

    for _ in range(n_simulations):
        true_proportions = np.random.dirichlet(np.ones(n_cell_types) * 0.5, 1)[0]
        y_true = true_proportions @ X
        noise = np.random.normal(0, 0.05 * np.std(y_true), y_true.shape)
        y_observed = y_true + noise

        coef, residual = nnls(X.T, y_observed)
        coef = coef / (coef.sum() + 1e-10)

        all_residuals.append(residual)
        all_mse.append(mean_squared_error([true_proportions], [coef]))
        corr = np.corrcoef(true_proportions, coef)[0, 1]
        if not np.isnan(corr):
            all_correlations.append(corr)

    if not all_correlations:
        return {"avg_corr": -1, "avg_mse": np.inf, "avg_residual": np.inf, "n_markers": n_regions}

    return {
        "avg_corr": np.mean(all_correlations),
        "avg_mse": np.mean(all_mse),
        "avg_residual": np.mean(all_residuals),
        "std_corr": np.std(all_correlations),
        "std_mse": np.std(all_mse),
        "std_residual": np.std(all_residuals),
        "n_markers": n_regions
    }

# ================================
# 4. 综合评分函数
# ================================
def calculate_composite_score(avg_corr, avg_mse, avg_residual, n_markers,
                              mse_max=0.1, residual_max=2000, marker_good=8500, marker_max=20000):
    corr_norm = (avg_corr + 1) / 2
    mse_norm = 0 if avg_mse >= mse_max else 1 - (avg_mse / mse_max)
    residual_norm = 0 if avg_residual >= residual_max else 1 - (avg_residual / residual_max)
    if n_markers <= marker_good:
        marker_norm = 1.0
    elif n_markers >= marker_max:
        marker_norm = 0.0
    else:
        marker_norm = 1 - (n_markers - marker_good) / (marker_max - marker_good)
    composite = 0.40 * corr_norm + 0.30 * mse_norm + 0.15 * residual_norm + 0.15 * marker_norm
    return composite, corr_norm, mse_norm, residual_norm, marker_norm

# ================================
# 5. 网格搜索
# ================================
def grid_search_optimization(ref_df, cpg_info, param_grid, verbose=True,
                             output_file="../R_script_data/grid_search_results.csv"):
    results = []
    param_names = list(param_grid.keys())
    param_combinations = list(product(*param_grid.values()))

    print(f"开始网格搜索，共 {len(param_combinations)} 个参数组合...")

    for i, params in enumerate(param_combinations):
        param_dict = dict(zip(param_names, params))
        try:
            ref_matrix, marker_df, tissue_sufficient_dict = build_reference_matrix(
                ref_df, cpg_info,
                param_dict["min_cpg_count"],
                param_dict["cv_threshold"],
                param_dict["top_n"]
            )

            if ref_matrix is None or ref_matrix.empty:
                if verbose:
                    print(f"[{i+1:3d}/{len(param_combinations)}] {param_dict} -> 无效矩阵")
                continue

            metrics = evaluate_reference_matrix(ref_matrix, n_simulations=35)
            if metrics["avg_corr"] < 0 or metrics["avg_mse"] > 1:
                if verbose:
                    print(f"[{i+1:3d}/{len(param_combinations)}] {param_dict} -> 无效评估")
                continue

            composite_score, corr_norm, mse_norm, residual_norm, marker_norm = calculate_composite_score(
                metrics["avg_corr"], metrics["avg_mse"], metrics["avg_residual"], metrics["n_markers"]
            )

            all_sufficient = all(tissue_sufficient_dict.values())
            failed_tissues_str = ",".join([t for t, ok in tissue_sufficient_dict.items() if not ok])

            result = {
                **param_dict,
                "composite_score": composite_score,
                "corr_norm": corr_norm,
                "mse_norm": mse_norm,
                "residual_norm": residual_norm,
                "marker_norm": marker_norm,
                **metrics,
                "all_tissues_sufficient": all_sufficient,
                "failed_tissues": failed_tissues_str
            }

            results.append(result)

            # 实时写出
            if len(results) == 1:
                pd.DataFrame([result]).to_csv(output_file, index=False, mode="w")
            else:
                pd.DataFrame([result]).to_csv(output_file, index=False, mode="a", header=False)

            if verbose:
                print(f"[{i+1:3d}/{len(param_combinations)}] {param_dict} -> 综合评分: {composite_score:.4f}, 所有组织满足: {all_sufficient}")

        except Exception as e:
            if verbose:
                print(f"[{i+1:3d}/{len(param_combinations)}] {param_dict} -> 错误: {e}")
            continue

    if not results:
        print("警告: 没有有效的参数组合")
        return pd.DataFrame(), {}

    results_df = pd.DataFrame(results)
    best_idx = results_df["composite_score"].idxmax()
    best_params = results_df.loc[best_idx, param_names].to_dict()

    print(f"\n网格搜索完成! 最优参数: {best_params}, 综合评分: {results_df.loc[best_idx, 'composite_score']:.4f}")
    return results_df, best_params

# ================================
# 6. 执行网格搜索
# ================================
param_grid = {
    "min_cpg_count": [5, 8, 10, 12, 15],
    "cv_threshold": [0.15, 0.2, 0.25, 0.3, 0.35],
    "top_n": [25, 50, 75, 100, 125, 150]
}

results_df, best_params = grid_search_optimization(
    ref_df, cpg_info, param_grid,
    verbose=True,
    output_file="grid_search_results_round1_with_marker_penalty.csv"
)
