# Build the SCP plot containers. The ones in the counts layer are layout
# placeholders required by the Seurat constructor, never measurements.
plot_shell <- function(features, ids, metadata) {
  stub <- matrix(1, length(features), length(ids), dimnames = list(features, ids))
  CreateSeuratObject(stub, meta.data = metadata)
}

# Use the observed group means and detection fractions directly. SCP 0.5.6
# calculates both from an expression matrix by default. Keep its original
# legend, geometry, label sizing and layout; locally substitute only the
# detection matrix (the means are already one observation per plotted group).
dot_from_summary <- function(values, genes, groups, ...) {
  stopifnot(as.character(packageVersion("SCP")) == "0.5.6",
            !anyDuplicated(values[c("gene","label")]),
            nrow(values) == length(genes) * length(groups))
  ix <- cbind(match(values$gene, genes), match(values$label, groups))
  mu <- detected <- matrix(NA_real_, length(genes), length(groups),
                           dimnames = list(genes, groups))
  mu[ix] <- values$mean_log_RNA
  detected[ix] <- values$percent_detected / 100
  stopifnot(all(is.finite(mu)), all(detected >= 0 & detected <= 1))
  ids <- sprintf("group%02d", seq_along(groups))
  obj <- plot_shell(genes, ids, data.frame(All.groups = factor(groups, levels = groups), row.names = ids))
  stored <- mu
  colnames(stored) <- ids
  obj <- SetAssayData(obj, assay = "RNA", layer = "data", new.data = stored)
  f <- SCP::GroupHeatmap
  original <- quote(mat_perc_list[[cell_group]] <- mat_perc)
  updated <- quote({
    stopifnot(identical(dimnames(mat_perc), dimnames(t(detection_fraction))))
    mat_perc_list[[cell_group]] <- t(detection_fraction)
  })
  hits <- 0L
  adapt <- function(e) {
    if (identical(e, original)) { hits <<- hits + 1L; return(updated) }
    if (is.call(e) && identical(e[[1]], as.name("for"))) {
      loop <- e[[4]]
      if (is.call(loop) && identical(loop[[1]], as.name("{"))) {
        for (k in seq_along(loop)[-1]) if (identical(loop[[k]], original)) loop[[k]] <- adapt(loop[[k]])
        e[[4]] <- loop
      }
    }
    e
  }
  code <- body(f)
  for (k in seq_along(code)[-1]) code[[k]] <- adapt(code[[k]])
  stopifnot(hits == 1L)
  body(f) <- code
  environment(f) <- list2env(list(detection_fraction = detected), parent = environment(SCP::GroupHeatmap))
  f(obj, features = genes, group.by = "All.groups", ...)
}
