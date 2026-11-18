# Figure source data and plotting code for the TRALI lung endothelial study

This release redraws Figures 1–4 from processed figure data. Figure 4 is the version with **both gCap and aCap** in panels B–D. It includes mouse-level values and selected statistical results for the figures. It does not contain reads, 10x matrices, gene-by-cell UMI counts, original cell barcodes, Seurat/AnnData objects, whole-transcriptome expression, or trained networks. No pre-rendered panels or arrow layers are inputs.

| Figure | What is shared |
| --- | --- |
| 1 (A–E) | Lung UMAP coordinates with mouse and lineage labels, marker means and detection percentages by lineage, composition per mouse, z scores for the 16 displayed heatmap genes, and six endothelial scores with mouse-level contrasts. |
| 2 (A–I) | Endothelial UMAP coordinates and subtype labels, the normalized expression values for the six displayed feature maps, marker summaries, one Fosl1 AUC per scored cell, the 185 scored members (Fosl1 plus 184 candidate targets), and mouse/subtype summary results. Plot-specific row indices link the activity scores to the public coordinates. |
| 3 (A–D) | Displayed enrichment tests, set membership for the Venn diagram, six selected targets with their per-mouse heatmap values and correlations, and corresponding statistical summaries. Published adjusted values retain their original testing families. |
| 4 (A–E) | Perturbation ranks, paired model predictions by mouse, two subtype-specific sets of UMAP coordinates and pseudotime, grid vectors and gradients, distances, and pathways. Each observed pseudotime file contains only two coordinates and the displayed time value. Arrows are drawn from numbers. |

The per-cell values are derived results and can be analyzed within their released scope. The release supports redrawing the figures and checking selected mouse-level comparisons. QC, clustering, normalization from counts, differential expression from counts, network inference, perturbation and pseudotime fitting were performed upstream. Plotting functions that require a Seurat container use an explicitly artificial constant layer for layout; biological values come from the source tables. Marker dot plots use the supplied groups' observed means and detection percentages.

## Repository layout

Each figure has a folder named `figure1` through `figure4`. Its `data/` directory contains the processed tables and its `plot/` directory contains the plotting scripts. Files ending in `.csv.gz` are gzip-compressed CSV tables. The two `figure4/data/4C_*` directories contain the general capillary (`4C_gCap`) and aerocyte (`4C_aCap`) pseudotime data.

## Download and run

[Download the repository ZIP](https://github.com/qxy2003/TRALI-Fosl1-figures/archive/refs/heads/main.zip) or clone the repository. Extract the ZIP and open a terminal in the folder containing `redraw.R` and `verify.py`.

The commands below assume R and Python with the packages already installed. The verified environment was Windows, using R 4.5.2, SCP 0.5.6, Seurat 5.5.0, ggplot2 4.0.1 and Cairo. Full R package versions are in `environment/figure1.txt`, `figure2.txt`, and `figure3.txt`. The marker-dot helper requires SCP 0.5.6.

~~~sh
python verify.py
Rscript --vanilla redraw.R 1
Rscript --vanilla redraw.R 2
Rscript --vanilla redraw.R 3
python redraw4.py
Rscript --vanilla check_mouse_results.R
python verify.py --figures
~~~

Figure 4 layout was verified with Python 3.12.13, Matplotlib 3.11.2, pandas 3.0.6 and NumPy 2.5.3. PDF assembly and pixel checks used Python 3.14.2, PyMuPDF 2.1 (MuPDF 1.28.6) and Pillow 13.0.0. If using separate environments, run `python redraw4.py --assembly-python PATH_TO_ASSEMBLY_PYTHON` with the plotting Python, and run `verify.py --figures` with the assembly Python.

Combined PDFs are written to `figureN/figure/FigureN.pdf`. Figures 1–3 also export separate panels under `figureN/figure/panels/`; Figure 4 writes A–E beside the combined PDF. The R check reproduces the unadjusted two-sided Welch tests, mean differences and 95% confidence intervals for Figure 1E (six scores) and Figure 2E–G (seven displayed contrasts and two sensitivity contrasts), using four LPS and four TRALI mice. It also checks the six exact paired Wilcoxon tests in Figure 4E and their BH correction. The actin–myosin filament sliding score uses the frozen 13-member GO:0033275 set. Its exploratory screening family contains 115 sets (BH q = 0.1069). Figure 2I uses a separate batch-adjusted limma analysis.

Each of the four combined PDFs rendered identically to the corresponding manuscript reference at 150 dpi in the recorded environment. Verification compares RGB pixels rather than binary PDF files; creation timestamps, fonts, graphics libraries or the renderer can affect comparisons across installations.

## Study and software

Single-cell transcriptomics implicates Fosl1 in the early pulmonary capillary endothelial response to antibody-mediated transfusion-related acute lung injury. The plotting scripts use SCP (https://github.com/zhanghao-njmu/SCP), Seurat, SCENIC, AUCell, ggplot2 and Matplotlib. Figure 3D uses a boxplot and beeswarm adapted from the van der Wijst et al. [plotting code](https://github.com/molgenis/scRNA-seq/blob/35faef87d416d3912a70da95407b498ceba646bc/R-scripts/interaction_plots.R#L317-L333).
