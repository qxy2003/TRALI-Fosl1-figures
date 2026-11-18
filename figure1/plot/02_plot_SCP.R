# Run from the Figure 1 package root. The only data inputs are in 分析结果.
if (.Platform$OS.type == "windows") invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
options(Seurat.object.assay.version = "v3", digits = 15, device = function(...) grDevices::pdf(NULL))
suppressPackageStartupMessages({library(Seurat); library(SCP); library(ggplot2); library(cowplot)})
out <- "figure"
dir.create(out, showWarnings = FALSE)
source("../plot_helpers.R", encoding = "UTF-8")
point <- read.csv("分析结果/lung_points.csv.gz")
order_table <- read.csv("分析结果/marker_order.csv")
lineages <- unique(order_table$lineage)
r <- list(metadata = point[c("Sample","broad")],
          markers = split(order_table$gene, factor(order_table$lineage, levels=lineages)))
r$metadata$Sample <- factor(r$metadata$Sample, levels=c(paste0("C",1:3),paste0("L",1:4),paste0("T",1:4)))
r$metadata$broad <- factor(r$metadata$broad, levels=c(lineages,"Unresolved"))
ids <- sprintf("point%06d",seq_len(nrow(point)))
rownames(r$metadata) <- ids
xy <- as.matrix(point[c("UMAP_1","UMAP_2")]); rownames(xy) <- ids
lung <- plot_shell("layout",ids,r$metadata)
lung[["umap"]] <- CreateDimReducObject(embeddings=xy,key="UMAP_",assay="RNA")
short <- c("Neutrophils", "Myeloid", "T cells / ILCs", "NK cells", "B cells",
           "Endothelial", "Epithelial", "Fibroblasts", "Mural", "Mesothelial",
           "Cardiomyocytes", "Unresolved")
lung$Lineage <- factor(lung$broad, levels = c(names(r$markers), "Unresolved"), labels = short)
pal <- setNames(c("#C68325", "#A9553A", "#695596", "#AC6799", "#427DA5",
                  "#167A78", "#4FAFC8", "#78923D", "#A38665", "#B689AF",
                  "#CD5568", "#A9B0B3"), short)
stopifnot(!anyNA(lung$Lineage))
base <- theme(text = element_text(family = "sans", colour = "#20282A"),
              plot.margin = margin(2, 4, 2, 4), legend.title = element_blank())

# A and B share a lineage scale and legend.
a <- SCP::CellDimPlot(lung, group.by = "Lineage", reduction = "umap", palcolor = pal,
    raster = TRUE, raster.dpi = c(1600, 1600), pt.size = 1.5, show_stat = FALSE,
    label = FALSE, title = NULL, theme_use = "theme_void", legend.position = "none")
a <- a + labs(title = NULL, subtitle = NULL, x = "UMAP 1", y = "UMAP 2") +
  theme_classic(base_size = 9) + coord_fixed() + base +
  theme(legend.position = "none", axis.text = element_blank(), axis.ticks = element_blank())
b <- SCP::CellStatPlot(lung, stat.by = "Lineage", group.by = "Sample", plot_type = "bar",
    stat_type = "percent", position = "stack", palcolor = pal, bg_alpha = 0,
    label = FALSE, NA_stat = FALSE, title = NULL, xlab = NULL, ylab = "Recovered cells (%)")
# Verify every plotted stack against the independently tabulated percentages.
tab <- read.csv("分析结果/composition_by_mouse.csv")
tab$Lineage <- short[match(tab$lineage, levels(r$metadata$broad))]
matched <- merge(b$data, tab, by = c("Sample", "Lineage"))
stopifnot(nrow(matched) == nrow(tab), max(abs(100 * matched$value - matched$percent)) < 1e-8)
b$data$Condition <- factor(tab$Group[match(b$data$Sample, tab$Sample)], c("Control", "LPS", "TRALI"))
# SCP uses black outlines on every segment. Thin white dividers preserve small populations.
for (i in seq_along(b$layers)) {
  if (inherits(b$layers[[i]]$geom, "GeomCol")) {
    b$layers[[i]]$aes_params$colour <- "white"
    b$layers[[i]]$aes_params$linewidth <- 0.13
  }
}
b <- b + facet_grid(~Condition, scales = "free_x", space = "free_x") +
  labs(title = NULL, subtitle = NULL, x = NULL) + theme_classic(base_size = 9) + base +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold"),
        legend.position = "none", axis.ticks.x = element_blank(), panel.spacing.x = grid::unit(2, "mm"))
shared_legend <- cowplot::get_legend(b + guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  theme(legend.position = "bottom", legend.text = element_text(size = 8.5),
        legend.key.size = grid::unit(3, "mm"), legend.margin = margin(0,0,0,0)))

# C: dot diameter indicates the fraction detected. Color shows mean log-normalized expression.
selected <- unlist(r$markers, use.names = FALSE)
ComplexHeatmap::ht_opt(legend_title_gp = grid::gpar(fontsize = 8, fontface = "bold"),
                       legend_labels_gp = grid::gpar(fontsize = 8))
# SCP suppresses its extra class strip and legend for a grouping named All.groups.
# Preserve the lineage factor so the rows and their labels still identify every class.
lung$All.groups <- lung$Lineage
ref <- read.csv("分析结果/lineage_marker_summary.csv")
ref$label <- short[match(ref$identity, levels(r$metadata$broad))]
hm <- dot_from_summary(ref, selected, short, assay = "RNA",
    slot = "data", aggregate_fun = mean, exp_method = "raw", exp_cutoff = 0, lib_normalize = FALSE,
    add_dot = TRUE, dot_size = grid::unit(4, "mm"), flip = TRUE, border = FALSE,
    nlabel = 0,
    row_title = "", column_title = "", row_names_side = "left", column_names_side = "bottom",
    show_row_names = TRUE, show_column_names = TRUE, column_names_rot = 50,
    cluster_rows = FALSE, cluster_columns = FALSE, use_raster = FALSE,
    exp_legend_title = "Mean log RNA", heatmap_palcolor = c("#FFF5F0", "#FC9272", "#DE2D26", "#80000F"),
    limits = c(0, ceiling(max(read.csv("分析结果/lineage_marker_summary.csv")$mean_log_RNA))),
    ht_params = list(row_names_gp = grid::gpar(fontsize = 8.5),
                     column_names_gp = grid::gpar(fontsize = 8.5, fontface = "italic")))
cplot <- hm$plot
ref <- read.csv("分析结果/lineage_marker_summary.csv")
mat <- hm$matrix_list[[1]]
ix <- cbind(match(ref$gene, rownames(mat)), match(short[match(ref$identity, levels(r$metadata$broad))], colnames(mat)))
stopifnot(!anyNA(ix), max(abs(mat[ix] - ref$mean_log_RNA)) < 1e-8)

# D: one column per endothelial pseudobulk library, ordered by condition and mouse.
rows <- read.csv("分析结果/heatmap_z.csv",row.names="gene")
samples <- read.csv("分析结果/heatmap_samples.csv",row.names="Sample")
animals <- rownames(samples); de_genes <- rownames(rows)
bulk <- list(z=as.matrix(rows[,animals]),samples=samples,selected=rows["direction"])
mouse_obj <- plot_shell(de_genes,animals,samples)
mouse_obj <- SetAssayData(mouse_obj, assay="RNA", layer="data", new.data=bulk$z)
conditions <- c("Control","LPS","TRALI")
mouse_obj$All.groups <- factor(mouse_obj$Group, levels=conditions)
dfit <- SCP::FeatureHeatmap(mouse_obj, features=de_genes, group.by="All.groups",
    cells=animals, cell_order=animals, assay="RNA", slot="data", exp_method="raw",
    lib_normalize=FALSE, limits=c(-2,2), heatmap_palcolor=c("#2166AC","#F7F7F7","#B2182B"),
    exp_legend_title="Row z", nlabel=0, cluster_rows=FALSE, cluster_columns=FALSE,
    cluster_row_slices=FALSE, cluster_column_slices=FALSE, use_raster=FALSE, border=FALSE,
    show_row_names=TRUE, show_column_names=TRUE, row_names_side="left",
    column_names_side="bottom", column_names_rot=0, column_title=conditions, row_title="",
    ht_params=list(row_names_gp=grid::gpar(fontsize=8,fontface="italic"),
                   column_names_gp=grid::gpar(fontsize=7.5),
                   column_title_gp=grid::gpar(fontsize=9,fontface="bold")))
displayed <- dfit$matrix_list[[1]]
stopifnot(identical(dimnames(displayed),dimnames(bulk$z)),
          max(abs(displayed-bulk$z)) < 1e-10)
# E: six literature-defined transcriptional scores, with one value per mouse in each plot.
scores <- read.csv("分析结果/endothelial_scores/scores_by_mouse.csv", row.names = "Sample")
tests <- read.csv("分析结果/endothelial_scores/TRALI_vs_LPS.csv")
stopifnot(all(tests$method=="Two-sided unpaired Welch t-test"),
          all(tests$n_LPS==4L),all(tests$n_TRALI==4L),
          all(!tests$paired),all(!tests$var_equal),all(!tests$batch_adjusted),
          all(tests$multiplicity_adjustment=="None"),all(tests$significance_basis=="unadjusted P"))
axis_labels <- c(TightJunction = "Tight junction score", LeukocyteAdhesion = "Leukocyte adhesion score",
                 ActinMyosinSliding = "Actin-myosin\nfilament sliding score", Complement = "Complement score",
                 Coagulation = "Coagulation score", Fibrinolysis = "Fibrinolysis score")
stopifnot(identical(rownames(scores), animals), identical(tests$module, names(axis_labels)),
          identical(as.character(scores$Group), as.character(mouse_obj$Group)),
          identical(as.integer(scores$n_cells), as.integer(bulk$samples[animals, "n_cells"])))
for (module in names(axis_labels)) mouse_obj[[module]] <- scores[[module]]
mouse_obj$Group <- factor(mouse_obj$Group, levels = conditions)
group_colors <- c(Control = "#D7DADC", LPS = "#A7C6DE", TRALI = "#E4A79B")
score_panels <- lapply(names(axis_labels), function(module) {
  p <- SCP::FeatureStatPlot(mouse_obj, stat.by = module, group.by = "Group", plot_type = "box",
    palcolor = group_colors, bg_alpha = 0, add_point = TRUE, pt.color = "#22343A", pt.size = 1.45,
    jitter.width = .13, jitter.height = 0, seed = 41, combine = FALSE)[[1]]
  # Keep the box median and the eleven individual points.
  p$layers <- Filter(function(l) !inherits(l$stat, "StatSummary"), p$layers)
  point_layer <- which(vapply(p$layers, function(l) inherits(l$geom, "GeomPoint") &&
                             inherits(l$stat, "StatIdentity"), logical(1)))
  # Spread animals within each group in sample order so close scores remain visible.
  offsets <- ave(seq_len(nrow(p$data)), p$data$group.by, FUN = function(j) seq(-.16, .16, length.out = length(j)))
  p$layers[[point_layer]]$position <- position_nudge(x = offsets)
  for (i in seq_along(p$layers)) {
    if (inherits(p$layers[[i]]$geom, "GeomBoxplot")) {
      p$layers[[i]]$aes_params$width <- .48
      p$layers[[i]]$aes_params$linewidth <- .4
    }
  }
  stopifnot(setequal(p$data$cells, animals), nrow(p$data) == length(animals),
            max(abs(p$data$value - scores[p$data$cells, module])) < 1e-10)
  value_range <- range(p$data$value)
  span <- diff(value_range)
  top <- value_range[2] + .105 * span
  p_value <- tests$p_value[match(module, tests$module)]
  stopifnot(is.finite(p_value), p_value >= 0, p_value <= 1)
  n_stars <- sum(p_value < c(.05, .01, .001, .0001))
  stopifnot(tests$label[match(module,tests$module)]==strrep("*",n_stars))
  p <- p + labs(title = NULL, subtitle = NULL, x = NULL, y = axis_labels[[module]]) +
    theme_classic(base_size = 9) +
    theme(legend.position = "none", strip.text = element_blank(), strip.background = element_blank(),
          plot.margin = margin(5, 5, 2, 4), axis.title.y = element_text(size = 8.1, lineheight = .9), axis.text.x = element_text(size = 7.5),
          axis.text = element_text(colour = "#20282A"), axis.ticks.x = element_blank()) +
    scale_y_continuous(breaks = pretty(c(max(0, value_range[1] - .10 * span), top), n = 4),
                       expand = expansion(mult = c(.09, .13))) +
    coord_cartesian(xlim = c(.5, 3.5), ylim = c(max(0, value_range[1] - .10 * span), top + .27 * span),
                    expand = FALSE)
  if (n_stars > 0L) {
    p <- p +
      annotate("segment", x = 2, xend = 3, y = top, yend = top, linewidth = .32) +
      annotate("segment", x = c(2, 3), xend = c(2, 3), y = top, yend = top - .035 * span, linewidth = .32) +
      annotate("text", x = 2.5, y = top + .055 * span, label = strrep("*", n_stars), size = 3.5)
  }
  plotted <- ggplot_build(p)$data[[point_layer]]
  stopifnot(nrow(plotted) == length(animals), !anyDuplicated(plotted$x),
            isTRUE(all.equal(sort(plotted$y), sort(scores[[module]]), tolerance = 1e-10)))
  text_layers <- which(vapply(p$layers,function(l) inherits(l$geom,"GeomText"),logical(1)))
  segment_layers <- which(vapply(p$layers,function(l) inherits(l$geom,"GeomSegment"),logical(1)))
  stopifnot(length(text_layers)==as.integer(n_stars>0L),length(segment_layers)==2L*as.integer(n_stars>0L))
  if (n_stars>0L) stopifnot(identical(as.character(ggplot_build(p)$data[[text_layers]]$label),strrep("*",n_stars)))
  p
})
letter <- function(label) {
  ggdraw() + draw_label(label, x = .007, y = .5, hjust = 0, fontface = "bold", size = 14)
}
block <- function(label, body) plot_grid(letter(label), body, ncol = 1, rel_heights = c(.26, 2.6))
ab <- plot_grid(block("A", a), block("B", b), nrow = 1, rel_widths = c(.88, 1.12))
de_titles <- plot_grid(letter("D"), letter("E"), nrow = 1, rel_widths = c(.78, 1.22))
de_panels <- plot_grid(dfit$plot, plot_grid(plotlist = score_panels, ncol = 2, align = "hv", axis = "tblr"),
                       nrow = 1, rel_widths = c(.78, 1.22))
fig <- plot_grid(ab, shared_legend, letter("C"), cplot, de_titles, de_panels,
                 ncol = 1, rel_heights = c(2.32, .38, .23, 2.95, .23, 4.38))
fig <- ggdraw() + draw_plot(fig, x = .012, y = .012, width = .976, height = .976)
ggsave(file.path(out, "Figure1.pdf"), fig, width = 8.5, height = 10.8,
       device = cairo_pdf, bg = "white")
ggsave(file.path(out, "Figure1.tiff"), fig, width = 8.5, height = 10.8,
       device = ragg::agg_tiff, dpi = 600, compression = "lzw", bg = "white")
ggsave(file.path(out, "Figure1.png"), fig, width = 8.5, height = 10.8,
       device = ragg::agg_png, dpi = 200, bg = "white")
writeLines(capture.output(sessionInfo()), file.path(out, "plot_session.txt"))
