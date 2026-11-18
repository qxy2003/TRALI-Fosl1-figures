# Check the statistics plotted in Figures 1E, 2E–G and 4E.
if (.Platform$OS.type == "windows") invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
read_data <- function(dir, file) read.csv(file.path(dir, file), check.names = FALSE)
verify_welch <- function(values, estimate, p, low, high) {
  mice <- subset(values, Group %in% c("TRALI","LPS"))
  stopifnot(nrow(mice) == 8L, !anyDuplicated(mice$Sample), all(table(mice$Group) == 4L))
  treated <- mice$value[mice$Group == "TRALI"]
  control <- mice$value[mice$Group == "LPS"]
  fit <- t.test(treated, control, var.equal = FALSE)
  stopifnot(max(abs(c(mean(treated)-mean(control)-estimate, fit$p.value-p,
                      fit$conf.int-c(low,high)))) < 1e-10)
}
dir1 <- "figure1/分析结果/endothelial_scores"
values <- read_data(dir1, "test_mouse_values.csv")
tests <- read_data(dir1, "TRALI_vs_LPS.csv")
stopifnot(nrow(tests) == 6L)
for (i in seq_len(nrow(tests))) {
  r <- tests[i,]
  verify_welch(values[values$module == r$module,], r$mean_difference,
               r$p_value, r$CI95_low, r$CI95_high)
}
dir2 <- "figure2/分析结果"
values <- read_data(dir2, "boxplot_test_mouse_values.csv")
tests <- read_data(dir2, "boxplot_significance.csv")
stopifnot(nrow(tests) == 9L)
for (i in seq_len(nrow(tests))) {
  r <- tests[i,]
  keep <- values$panel == r$panel & values$population == r$population
  verify_welch(values[keep,], r$estimate, r$p, r$CI_low, r$CI_high)
}
cat("Figures 1 and 2: all 15 Welch contrasts and confidence intervals match.\n")
dir4 <- "figure4/分析结果"
s <- read.csv(file.path(dir4,"4E_pathways_by_mouse.csv"))
t <- read.csv(file.path(dir4,"4E_paired_tests.csv"))
stopifnot(!anyDuplicated(s[c("module","mouse")]),nrow(t)==6L,
          max(abs(s$model_KO - s$model_WT - s$delta)) < 1e-12)
paired_p <- vapply(t$module, function(m) {
  d4 <- s[s$module == m, ]
  stopifnot(setequal(d4$mouse, paste0("T",1:4)))
  wilcox.test(d4$delta, mu=0, alternative="two.sided", exact=TRUE)$p.value
}, numeric(1))
stopifnot(max(abs(paired_p-t$p_exact_two_sided),abs(p.adjust(paired_p,"BH")-t$q_BH_six_modules))<1e-12)
cat("Figure 4E: four paired mouse-level predictions per module; six exact p and BH q values match.\n")
