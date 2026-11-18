# Run from the repository root with Rscript --vanilla redraw.R 1 (or 2, 3).
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L || !args %in% c("1", "2", "3")) {
  stop("Usage: Rscript --vanilla redraw.R <1|2|3>")
}
if (.Platform$OS.type == "windows") {
  invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
}
setwd(paste0("figure", args))
script <- switch(args, "1" = "plot/02_plot_SCP.R",
                 "2" = "作图代码/02_plot.R", "3" = "作图代码/05_plot.R")
source(script, encoding = "UTF-8")

# Figures 1 and 2 originally saved the combined plate. Also export their plot objects.
save_letter <- function(letter, plot, width, height) {
  folder <- file.path("figure", "panels")
  dir.create(folder, showWarnings = FALSE)
  ggplot2::ggsave(file.path(folder, paste0("Figure", args, letter, ".pdf")),
                 plot = plot, width = width, height = height, device = grDevices::cairo_pdf,
                 bg = "white")
}
with_legend <- function(panel, key) {
  cowplot::plot_grid(panel, key, ncol = 1, rel_heights = c(1, .22))
}
if (args == "1") {
  save_letter("A", with_legend(a, shared_legend), 7, 5)
  save_letter("B", with_legend(b, shared_legend), 7, 4)
  save_letter("C", cplot, 8.5, 3.6)
  save_letter("D", dfit$plot, 4.2, 5)
  save_letter("E", cowplot::plot_grid(plotlist = score_panels, ncol = 2), 7, 6)
}
if (args == "2") {
  save_letter("A", with_legend(a, shared), 8, 3.6)
  save_letter("B", with_legend(b, shared), 8, 3.6)
  save_letter("C", cplot, 5.1, 3.7)
  save_letter("D", d, 6, 4.3)
  save_letter("E", e, 6.5, 3)
  save_letter("F", f, 5.5, 3)
  save_letter("G", g, 10.8, 3.5)
  save_letter("H", h, 6, 5)
  save_letter("I", i_panel, 6, 5)
}
