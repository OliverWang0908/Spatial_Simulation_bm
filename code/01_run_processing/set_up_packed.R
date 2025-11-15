source("/Users/oukanyou/Desktop/SpatialSimBench/code/01_run_processing/S1setup.R")

# ============================================================
# Batch process all .rds files in a folder using the functions defined above
# ============================================================

library(SummarizedExperiment)
library(SingleCellExperiment)

# Path to the folder containing .rds files
data_dir <- "/Users/oukanyou/Desktop/SpatialSimBench/data"

# List all .rds files under the directory
rds_files <- list.files(data_dir, pattern = "\\.rds$", full.names = TRUE)

# Create a folder for processed results
output_dir <- file.path(data_dir, "processed_results")
if (!dir.exists(output_dir)) dir.create(output_dir)

# Loop over each .rds file and process it
for (file_path in rds_files) {
  cat("Processing:", file_path, "\n")
  
  # Extract a short name for saving results
  sample_name <- extract_name(file_path)
  
  # ------------------------------------------------------------
  # Step 1. Read the RDS file
  # ------------------------------------------------------------
  sce_obj <- readRDS(file_path)
  real sce <- reorder_matrix_real(sce_obj)
  print('reorder finish')
  
  
  # ------------------------------------------------------------
  # Step 2. Rebuild the SCE object (standardize spatial info)
  # ------------------------------------------------------------
  sce_rebuilt <- rebuilt_sce(sce_obj)
  
  # ------------------------------------------------------------
  # Step 3. Compute logcounts if not already present
  # ------------------------------------------------------------
  sce_rebuilt <- compute_logcounts(sce_rebuilt)
  
  # ------------------------------------------------------------
  # Step 4. Compute PCA for dimensionality reduction
  # ------------------------------------------------------------
  sce_pca <- computePCA(sce_rebuilt)
  
  # ------------------------------------------------------------
  # Step 5. (Optional) Run SimBench evaluation
  # ------------------------------------------------------------
  # Here, you would compare a simulated dataset with a real one.
  # For example, if you have both real and simulated versions:
  # fig <- simbench_result(sce_real = real_sce, sce_sim = sce_pca)
  # ggsave(paste0(output_dir, "/", sample_name, "_simbench_plot.pdf"), fig)
  
  # ------------------------------------------------------------
  # Step 6. Save the processed SCE object
  # ------------------------------------------------------------
  saveRDS(sce_pca, file = file.path(output_dir, paste0(sample_name, "_processed.rds")))
  
  cat("Finished:", sample_name, "\n\n")
}

cat("✅ All files processed and saved to:", output_dir, "\n")





##########
rebuilt_sce <- function(sce_obj) {
  counts_single <- as.matrix(assay(sce_obj, "counts"))
  
  # Extract colData and ensure rownames match counts
  col_data_df <- as.data.frame(colData(sce_obj))
  rownames(col_data_df) <- colnames(counts_single)  # Fix the rowname mismatch
  
  # Assign spatial coordinates
  col_data_df$spatial1 <- col_data_df$row
  col_data_df$spatial2 <- col_data_df$col
  
  # Optionally rename columns as "spatial1xspatial2"
  colnames(counts_single) <- paste0(col_data_df$spatial1, "x", col_data_df$spatial2)
  
  # Build new SCE object
  sce_new <- SingleCellExperiment(
    assays = list(counts = counts_single),
    colData = col_data_df[, c("spatial1", "spatial2")]
  )
  
  return(sce_new)
}



for (file_path in rds_files) {
  cat("Processing:", file_path, "\n")
  
  # Extract a short name for saving results
  sample_name <- extract_name(file_path)
  
  # ------------------------------------------------------------
  # Step 1. Read the RDS file
  # ------------------------------------------------------------
  sce_obj <- readRDS(file_path)
  counts_matrix <- assay(sce_obj, "counts")
  real_sce <- reorder_matrix(counts_matrix)
  real_sce <- real_sce[, !colnames(real_sce) %in% "Gene"]
  print('reorder finish')
  
  }
counts_df_real <- reorder_matrix_real_new(real_sce)
counts_df_real <- counts_df_real[, !colnames(real_sce) %in% "Gene"]

sce_obj <- compute_logcounts(sce_obj)
