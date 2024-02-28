#' pseudoBulk
#'
#' @title
#' Aggregate single-cell data into pseudo-bulk counts
#'
#' @description
#' Aaggregates single-cell counts into pseudo-bulk counts per sample, 
#' and stores the pseudo-bulk counts in a DGEList object along with any metadata specified as groups for further differential expression analysis. 
#' If no barcode is specified, it aggregates all cells.
#' It enforces a minimum sample size per group as specified by the threshold.
#'
#' @param sc.obj A single-cell object in either Seurat or SingleCellExperiment format, containing count data and sample IDs.
#' @param assay The name of the assay to pull counts from, with "counts" as the default.
#' @param sample A string specifying the column in the metadata that contains sample information.
#' @param group Optional. A string specifying the column in the metadata that contains group information for differential expression.
#' @param barcode Optional. A vector of two elements where the first element is the name of the barcode column from the metadata,
#' and the second element is the specific barcode ID to filter cells by. This argument is used to subset cells based on barcode information.
#' @param threshold The minimum number of samples for each group to be included in the analysis.
#' A threshold less than 2 will cause the function to stop with an error.
#'
#' @return Returns a DGE object with normalised counts for a sample.
#' @export
#' 

pseudoBulk <- function(sc.obj, 
                       assay = "counts", 
                       sample, 
                       group = NULL,
                       barcode = NULL,
                       threshold = 2) {
  
  if (threshold < 2) {
    stop("Threshold cannot be smaller than 2.")
  }
  
  # Get metadata and identify class
  if (class(sc.obj)[1] == "Seurat") {
    meta <- sc.obj@meta.data
    counts <- GetAssayData(sc.obj, assay = assay, slot = "counts")
    type <- "Seurat"
  } else {
    if (class(sc.obj)[1] == "SingleCellExperiment") {
      meta <- sc.obj@colData
      counts <- assay(sc.obj, assay)
      type <- "SingleCellExperiment"
    } else {
      stop("A single cell object must be supplied in Seurat or SingleCellExperiment format")
    }  
  }
  
  valid_cells <- !is.na(meta[[sample]])
  meta <- meta[valid_cells, ]
  counts <- counts[, valid_cells]
  
  # Check inputs
  if (!sample %in% colnames(meta)) {
    stop("sample column not found in sc.obj metadata")
  }
  if (!is.null(group) && !group %in% colnames(meta)) {
    stop("group column not found in sc.obj metadata")
  }
  
  # Check for consistency within each sample regarding the group
  if (!is.null(group)) {
    inconsistentSamples <- character()
    for (samp in unique(meta[[sample]])) {
      groups <- meta[meta[[sample]] == samp, group]
      if (length(unique(na.omit(groups))) != 1 || any(is.na(groups))) {
        inconsistentSamples <- c(inconsistentSamples, samp)
      }
    }
    if (length(inconsistentSamples) > 0) {
      stop(sprintf("Inconsistencies found in group values within samples: %s. Each sample should correspond to only one group, including NA as inconsistency.", paste(inconsistentSamples, collapse=", ")))
    }
  }
  
  # Filter cells by barcode if provided
  if (!is.null(barcode)) {
    if (!barcode[1] %in% colnames(meta)) {
      stop("Barcode column not found in sc.obj metadata.")
    }
    meta <- meta[!is.na(meta[[barcode[1]]]) & meta[[barcode[1]]] == barcode[2],]
    if (nrow(meta) == 0) {
      stop("No cells match the specified barcode.")
    }
    counts <- counts[, rownames(meta)]
  }
  
  if (!is.null(group)) {
    samples_per_group <- table(meta[[group]])
    if (any(samples_per_group < threshold)) {
      stop("After filtering, not all groups have the minimum number of samples specified by the threshold.")
    }
  }
  
  # Aggregate counts by sample
  unique_samples <- names(table(meta[[sample]]))[table(meta[[sample]]) > 1]
  if(length(unique_samples) < length(table(meta[[sample]]))){
    print("Removing samples with less than 2 cells.")
  }
  agg_counts <- sapply(unique_samples, function(samp) {
    rowSums(counts[, meta[[sample]] == samp])
  })
  
  # Filter based on threshold
  if (!is.null(group)) {
    group_counts <- table(meta[[sample]], meta[[group]])
    valid_samples <- unique_samples[colSums(group_counts >= threshold) > 0]
    agg_counts <- agg_counts[, valid_samples]
  }
  
  # Prepare sample data for DGEList
  sample_data <- data.frame(sample = colnames(agg_counts))
  if (!is.null(group)) {
    # Extract group for each sample
    sample_data$group <- sapply(sample_data$sample, function(s) {
      unique(meta[meta[[sample]] == s, group])[1]
    })
  }
  
  # Create DGEList object
  dge <- DGEList(counts = agg_counts)
  dge$samples <- sample_data
  return(dge)
}