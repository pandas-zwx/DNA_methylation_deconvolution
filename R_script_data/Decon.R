library(data.table)

setwd("~/workspace/Atlas_WGBS/tissue_cell_type_methylation_density/")

file_list <- list.files(pattern = "\\.methylation_density\\.bed$")

# 先给每个文件生成唯一列名
cell_types <- sub(
  "^GSM[0-9]+_(.*)-Z[0-9A-Z]+\\.hg38\\.methylation_density\\.bed$",
  "\\1",
  file_list
)
cell_types_unique <- make.unique(cell_types)

dt_list <- vector("list", length(file_list))

for (i in seq_along(file_list)) {
  f <- file_list[i]
  
  dt <- fread(
    f,
    header = FALSE,
    sep = "\t",
    col.names = c("chr", "start", "end", "value"),
    showProgress = FALSE
  )
  
  dt[, region_id := paste(chr, start, end, sep = "_")]
  dt[, sample := cell_types_unique[i]]
  
  dt_list[[i]] <- dt[, .(region_id, sample, value)]
}

long_dt <- rbindlist(dt_list, use.names = TRUE)

wide_dt <- dcast(
  long_dt,
  region_id ~ sample,
  value.var = "value"
)

wide_df <- as.data.frame(wide_dt)
rownames(wide_df) <- wide_df$region_id
wide_df$region_id <- NULL

saveRDS(wide_df, file = "../R_script_data/combined_methylation_density_matrix.rds")
