# Regenerate the four panels from the analysis tables, at the Figure 3 root.
if(.Platform$OS.type=="windows") invisible(Sys.setlocale("LC_ALL","English_United States.utf8"))
suppressPackageStartupMessages({library(ggplot2);library(cowplot);library(enrichplot);library(Seurat);library(SCP)})
options(Seurat.object.assay.version="v3",device=function(...) grDevices::pdf(NULL))
out <- "figure"; singles <- file.path(out,"panels")
dir.create(singles,recursive=TRUE,showWarnings=FALSE)
read_result <- function(f) read.csv(file.path("分析结果",f),check.names=FALSE)
d <- read_result("key_targets_heatmap.csv"); reps <- read_result("key_targets.csv")
terms <- read_result("GO_displayed.csv")
ink <- "#20282A"; muted <- "#50676C"
# Figure 1D pseudobulk z scale; treatment and red palettes shared by Figures 1 and 2.
zcols <- c("#2166AC","#F7F7F7","#B2182B")
group_colors <- c(LPS="#A7C6DE",TRALI="#E4A79B")
red <- c("#FFF5F0","#FC9272","#DE2D26","#80000F")
tag_panel <- function(plot, letter) ggdraw(plot)+
 draw_label(letter,x=.003,y=.995,hjust=0,vjust=1,fontface="bold",size=15,colour=ink)
wid <- 8.75
make_page <- function(p) ggdraw()+draw_plot(p,x=.018,y=.02,width=.962,height=.96)
save_panel <- function(plot,name,width,height){
 grDevices::cairo_pdf(paste0(name,".pdf"),width=width,height=height,family="sans");print(plot);invisible(dev.off())
 ragg::agg_png(paste0(name,".png"),width=width,height=height,units="in",res=200,background="white");print(plot);invisible(dev.off())
 ragg::agg_tiff(paste0(name,".tiff"),width=width,height=height,units="in",res=600,compression="lzw",background="white");print(plot);invisible(dev.off())
}
# One pseudobulk library per animal and subtype, matching Figure 1D.
ids <- c(paste0("L",1:4),paste0("T",1:4))
pops <- c("gCap","aCap")
md <- expand.grid(Sample=ids,subtype=pops,stringsAsFactors=FALSE)
md$Group <- factor(ifelse(startsWith(md$Sample,"L"),"LPS","TRALI"),levels=names(group_colors))
md$All.groups <- factor(md$subtype,levels=pops)
rownames(md) <- paste(md$subtype,md$Sample,sep=".")
key <- paste(d$subtype,d$Sample,sep=".")
indices <- cbind(match(d$gene,reps$gene),match(key,rownames(md)))
stopifnot(nrow(d)==nrow(reps)*nrow(md),!anyNA(indices),!anyDuplicated(data.frame(indices)),
          all(as.character(md$Group[indices[,2]])==d$Group))
zmat <- matrix(NA_real_,nrow(reps),nrow(md),dimnames=list(reps$gene,rownames(md)))
source("../plot_helpers.R",encoding="UTF-8")
zmat[indices] <- d$within_subtype_z
stopifnot(all(is.finite(zmat)))
mouse <- plot_shell(rownames(zmat),colnames(zmat),md)
mouse <- SetAssayData(mouse,assay="RNA",layer="data",new.data=zmat)
ComplexHeatmap::ht_opt(legend_title_gp=grid::gpar(fontsize=9,fontface="bold"),
                       legend_labels_gp=grid::gpar(fontsize=8))
hm <- SCP::FeatureHeatmap(mouse,features=rownames(zmat),cells=colnames(zmat),cell_order=colnames(zmat),
  group.by="All.groups",assay="RNA",slot="data",exp_method="raw",lib_normalize=FALSE,
  limits=c(-2,2),heatmap_palcolor=zcols,exp_legend_title="RNA z score",
  cell_annotation="Group",cell_annotation_palcolor=list(group_colors),
  cell_annotation_params=list(height=grid::unit(3,"mm"),annotation_name_gp=grid::gpar(fontsize=8)),
  nlabel=0,cluster_rows=FALSE,cluster_columns=FALSE,cluster_row_slices=FALSE,cluster_column_slices=FALSE,
  use_raster=FALSE,border=FALSE,show_row_names=TRUE,row_names_side="left",
  show_column_names=TRUE,column_names_side="bottom",column_names_rot=0,row_title="",column_title=pops,
  ht_params=list(column_labels=md$Sample,column_gap=grid::unit(4,"mm"),
     row_names_gp=grid::gpar(fontsize=10,fontface="italic"),column_names_gp=grid::gpar(fontsize=8.5),
     column_title_gp=grid::gpar(fontsize=10,fontface="bold")))
actual <- hm$matrix_list[[1]]
stopifnot(identical(dimnames(actual),dimnames(zmat)),max(abs(actual-zmat))<1e-10)
C <- ggdraw()+draw_plot(hm$plot,x=.022,y=.025,width=.973,height=.90)+
 draw_label("C",x=.003,y=.995,hjust=0,vjust=1,fontface="bold",size=15,colour=ink)

# Use the installed enrichplot barplot method for each subtype, retaining the joint BH q values.
up <- terms[terms$direction=="Up",]
ordered <- unique(up[order(up$display_order),c("ID","term")])
description <- function(term)paste0(toupper(substr(term,1,1)),substring(term,2))
levels_y <- rev(description(ordered$term))
stopifnot(setequal(unique(up$subtype),pops),all(up$q<.05),!anyDuplicated(up[c("ID","subtype")]))
native_bars <- function(pop){
 t <- up[up$subtype==pop,]
 tab <- data.frame(ID=t$ID,Description=description(t$term),GeneRatio=paste(t$hits,t$annotated_query,sep="/"),
                   BgRatio=paste(t$set_size,t$background,sep="/"),pvalue=t$p,p.adjust=t$q,
                   qvalue=NA_real_,geneID=gsub(";","/",t$genes,fixed=TRUE),Count=t$hits,row.names=t$ID)
 obj <- methods::new("enrichResult",result=tab,pvalueCutoff=1,qvalueCutoff=1,pAdjustMethod="BH",
                     organism="Mus musculus",ontology="BP",keytype="SYMBOL")
 p <- graphics::barplot(obj,x="Count",color="p.adjust",showCategory=nrow(tab),font.size=10)
 ix <- match(p$data$ID,t$ID)
 stopifnot(inherits(p$layers[[1]]$geom,"GeomCol"),nrow(p$data)==nrow(t),!anyNA(ix),
           all(p$data$Count==t$hits[ix]),max(abs(p$data$p.adjust-t$q[ix]))<1e-12)
 p$data$subtype <- factor(pop,levels=pops)
 p$data$Description <- factor(as.character(p$data$Description),levels=levels_y)
 p
}
bar_objects <- lapply(pops,native_bars)
bars <- bar_objects[[1]]
bars$layers[[1]]$geom_params$width <- .70
bars$data <- do.call(rbind,lapply(bar_objects,function(p)p$data))
xmax <- ceiling(max(up$hits)/5)*5
label_term <- function(x)vapply(x,function(a)paste(strwrap(a,width=43),collapse="\n"),character(1))
bars <- suppressMessages(bars+
 facet_grid(cols=vars(subtype))+
 scale_y_discrete(labels=label_term)+
 scale_x_continuous(limits=c(0,xmax),breaks=seq(0,xmax,5),expand=expansion(mult=c(0,.03)))+
 scale_fill_gradientn(colours=rev(red),trans="log10",limits=c(min(up$q),.05),
                     breaks=c(.0001,.001,.01),labels=c("0.0001","0.001","0.01"),name="BH q")+
 labs(x=NULL,y=NULL,title=NULL)+
 theme_minimal(base_family="sans",base_size=10)+
 theme(panel.grid.major.y=element_blank(),panel.grid.minor=element_blank(),
       panel.grid.major.x=element_line(colour="#E4EAE8",linewidth=.28),
       panel.spacing.x=grid::unit(12,"pt"),axis.text.y=element_text(size=8.7,colour=ink,lineheight=.93),
       axis.text.x=element_text(size=8.2,colour=muted),axis.title.x=element_text(size=8.8,colour=ink),
       strip.text.x=element_text(size=10.4,face="bold",colour=ink,margin=margin(0,0,5,0)),
       legend.title=element_text(size=8.7),legend.text=element_text(size=8),
       legend.key.width=grid::unit(.36,"cm"),legend.key.height=grid::unit(.53,"cm"),
       plot.margin=margin(5,2,2,13)))
A <- tag_panel(bars,"A")

# The two subtype-specific DE sets and the correlation set share the Figure 2 target universe.
membership <- read_result("Figure3B_venn_membership.csv")
sets <- list(gCap=with(membership,gene[up_in_gCap]),
             aCap=with(membership,gene[up_in_aCap]),
             positive=with(membership,gene[positive_RNA_and_LOO_in_gCap_and_aCap]))
stopifnot(!anyDuplicated(membership$gene),setequal(Reduce(intersect,sets),reps$gene))
labels <- c(sprintf("Upregulated in gCap\n(n = %d)",length(sets$gCap)),
            sprintf("Upregulated in aCap\n(n = %d)",length(sets$aCap)),
            sprintf("Positive correlation with Fosl1 RNA and regulon activity\nin both subtypes (n = %d)",length(sets$positive)))
diagram <- VennDiagram::venn.diagram(
  x=sets,filename=NULL,disable.logging=TRUE,category.names=labels,
  euler.d=FALSE,scaled=FALSE,
  fill=c("#A8C5DA","#E6B8A7","#94C7B3"),alpha=rep(.52,3),
  col=c("#6C8EA4","#B58270","#649283"),lwd=1.05,
  cex=c(1,1,1,.86,1.34,.86,1),fontface=c(rep("plain",4),"bold","plain","plain"),
  fontfamily="sans",label.col=c(rep(ink,3),muted,ink,muted,ink),
  cat.pos=c(-80,80,180),cat.just=list(c(1,.5),c(0,.5),c(.5,1)),cat.dist=c(.055,.055,.035),
  cat.cex=.88,cat.col=ink,cat.fontfamily="sans",margin=.035)
circle_frame <- grid::viewport(y=.55,width=grid::unit(2.35,"in"),height=grid::unit(2.35,"in"))
B <- tag_panel(ggdraw()+draw_grob(grid::grobTree(children=diagram,vp=circle_frame)),"B")
# Each point is a within-mouse, within-subtype partial Spearman correlation.
# activity_LOO excludes the measured target from the Figure 2 AUCell signature.
panel_genes <- as.character(reps$gene)
mouse_rho <- read_result("correlations_each_mouse.csv")
points_c <- subset(mouse_rho, gene %in% panel_genes & measure == "activity_LOO")
points_c$gene <- factor(points_c$gene, levels=panel_genes)
points_c$subtype <- factor(points_c$subtype, levels=pops)
points_c$Sample <- factor(points_c$Sample, levels=ids)
points_c$Group <- factor(points_c$Group, levels=names(group_colors))
points_c <- points_c[order(points_c$gene,points_c$subtype,points_c$Sample),]
stopifnot(!anyNA(points_c), all(table(points_c$gene,points_c$subtype,points_c$Sample)==1),
          all(table(points_c$gene,points_c$subtype,points_c$Group)==4),
          all(as.character(points_c$Group)==ifelse(startsWith(as.character(points_c$Sample),"L"),"LPS","TRALI")),
          all(abs(points_c$rho)<=1), all(points_c$n_cells>0))

# Check mouse-level means against the full results before rendering.
summ <- subset(read_result("correlations_all_targets.csv"),gene %in% panel_genes & measure=="activity_LOO")
means <- aggregate(rho~gene+subtype,data=points_c,FUN=mean)
group_idx <- match(paste(means$gene,means$subtype),paste(summ$gene,summ$subtype))
stopifnot(nrow(summ)==length(panel_genes)*length(pops),!anyNA(group_idx),max(abs(means$rho-summ$mean_rho[group_idx]))<1e-12)
write.csv(points_c,file.path(singles,"Figure3D_source_data.csv"),row.names=FALSE)

# Adapted from van der Wijst et al., Nat Genet (2018), Fig. 2b, lines 317–333:
# https://github.com/molgenis/scRNA-seq/blob/35faef87d416d3912a70da95407b498ceba646bc/R-scripts/interaction_plots.R#L317-L333
# Give each treatment its own box and swarm within each subtype.
offset <- c(LPS=-.17, TRALI=.17)
points_c$x <- as.integer(points_c$subtype) + unname(offset[as.character(points_c$Group)])
rho_plot <- ggplot(points_c,aes(x=x,y=rho,group=x))+
 geom_hline(yintercept=0,linetype="dashed",linewidth=.42,colour="#73828A")+
 geom_boxplot(aes(fill=Group),width=.26,orientation="x",position="identity",
              outlier.shape=NA,notch=FALSE,colour=ink,linewidth=.42)+
 ggbeeswarm::geom_beeswarm(cex=4.2,corral="gutter",corral.width=.30,size=1.9,shape=21,stroke=.39,
                           colour=ink,fill=scales::alpha("white",.72),show.legend=FALSE)+
 facet_wrap(vars(gene),ncol=2,scales="fixed",axes="all_x",axis.labels="all_x")+
 scale_fill_manual(values=group_colors,drop=FALSE)+
 scale_x_continuous(breaks=seq_along(pops),labels=pops,limits=c(.57,2.43),expand=c(0,0))+
 scale_y_continuous(limits=c(min(-.08,floor(min(points_c$rho)*10)/10),ceiling(max(points_c$rho)*5)/5),
                    breaks=scales::breaks_width(.2),expand=expansion(mult=c(0,.012)))+
 labs(x=NULL,y="Within-mouse partial Spearman ρ",fill="Group")+
 theme_classic(base_family="sans",base_size=10)+
 theme(axis.text=element_text(size=9,colour=ink),axis.title.y=element_text(size=9.3,colour=ink),
       axis.line=element_line(linewidth=.35,colour=ink),axis.ticks=element_line(linewidth=.3,colour=ink),
       strip.background=element_blank(),strip.text.x=element_text(size=11,face="italic",colour=ink),
       panel.spacing.x=grid::unit(11,"mm"),panel.spacing.y=grid::unit(3.5,"mm"),legend.position="right",legend.title=element_text(size=9),
       legend.text=element_text(size=8.6),legend.key.height=grid::unit(.43,"cm"),
       legend.key.width=grid::unit(.42,"cm"),plot.margin=margin(2,7,1,8))
built_c <- ggplot_build(rho_plot)
boxes <- built_c$data[[2]]
dots <- built_c$data[[3]]
stopifnot(nrow(boxes)==24L,nrow(dots)==nrow(points_c))
for (i in seq_len(nrow(boxes))) {
  box <- boxes[i,]
  gene <- built_c$layout$layout$gene[match(box$PANEL,built_c$layout$layout$PANEL)]
  src <- points_c[points_c$gene==gene & abs(points_c$x-box$x)<1e-10,]
  pts <- dots[dots$PANEL==box$PANEL & dots$group==box$group,]
  quartiles <- quantile(src$rho,c(.25,.5,.75),names=FALSE)
  stopifnot(nrow(src)==4L,nrow(pts)==4L,
            max(abs(sort(pts$y)-sort(src$rho)))<1e-12,all(abs(pts$x-box$x)<.17),
            identical(box$fill,unname(group_colors[as.character(src$Group[1])])),
            max(abs(unlist(box[c("lower","middle","upper")],use.names=FALSE)-quartiles))<1e-12)
}
D <- tag_panel(rho_plot,"D")

combined <- plot_grid(A,B,C,D,ncol=1,rel_heights=c(3.70,2.68,2.62,4.95))
save_panel(make_page(combined),file.path(out,"Figure3"),wid,14.55)
save_panel(make_page(B),file.path(singles,"Figure3B_target_venn"),wid,2.8)
save_panel(make_page(D),file.path(singles,"Figure3D_mouse_correlations"),wid,5.15)
save_panel(make_page(A),file.path(singles,"Figure3A_GO_barplots"),wid,3.85)
save_panel(make_page(C),file.path(singles,"Figure3C_target_heatmap"),wid,2.86)
# Each native bar chart is also available independently with its own GO labels.
for (pop in pops){
 single <- bars
 single$data <- droplevels(bars$data[bars$data$subtype==pop,,drop=FALSE])
 save_panel(make_page(single),file.path(singles,paste0("Figure3A_",pop,"_up_GO")),6.5,3.55)
}
writeLines(capture.output(sessionInfo()),file.path(out,"plot_session.txt"))
