# Run at the Figure 2 package root. All scientific inputs are in 分析结果/.
if (.Platform$OS.type == "windows") invisible(Sys.setlocale("LC_ALL","English_United States.utf8"))
options(Seurat.object.assay.version="v3", digits=15, device=function(...) grDevices::pdf(NULL))
suppressPackageStartupMessages({library(Seurat); library(SCP); library(ggplot2); library(cowplot);
  library(SCENIC); library(AUCell); library(SummarizedExperiment)})
input <- "分析结果"
out <- "figure"
dir.create(out,showWarnings=FALSE)
source("../plot_helpers.R",encoding="UTF-8")
point <- read.csv(file.path(input,"endothelial_points.csv.gz"))
ids <- sprintf("point%05d",seq_len(nrow(point)))
md <- point[c("Sample","Group","Subtype")]; rownames(md) <- ids
md$Group <- factor(md$Group,c("Control","LPS","TRALI"))
md$Sample <- factor(md$Sample,c(paste0("C",1:3),paste0("L",1:4),paste0("T",1:4)))
md$Subtype <- factor(md$Subtype,read.csv(file.path(input,"subtype_order.csv"))$Subtype)
xy0 <- as.matrix(point[c("UMAP_1","UMAP_2")]); rownames(xy0) <- ids
result <- list(umap=xy0,meta=md)
features <- c("Cdh5","Gpihbp1","Car4","Gja5","Slc6a2","Ccl21a")
shown <- t(as.matrix(point[features])); colnames(shown) <- ids

table_in <- function(filename) read.csv(file.path(input,filename),check.names=FALSE)
composition <- table_in("composition_by_mouse.csv")
markers <- table_in("marker_summary.csv")
rna <- table_in("Fosl1_RNA_by_mouse.csv")
lineage <- table_in("lineage_AUC_by_mouse.csv")
stopifnot(nrow(lineage)==22L,identical(unique(lineage$signature),185L))
cap <- table_in("capillary_AUC_by_mouse.csv")
effects <- table_in("capillary_effects.csv")
annotations <- table_in("boxplot_significance.csv")
stopifnot(nrow(annotations)==9L,!anyDuplicated(annotations[c("panel","population")]))
plot_labels <- character(0)
box_label <- function(panel_name,population_name) {
  row <- annotations[annotations$panel==panel_name & annotations$population==population_name,]
  stopifnot(nrow(row)==1L,is.finite(row$p),row$contrast=="TRALI vs LPS",
            row$n_ref==4L,row$n_case==4L,
            row$method=="Two-sided unpaired Welch t-test",!row$paired,!row$var_equal,
            row$significance_basis=="unadjusted P",!row$batch_adjusted,
            row$multiplicity_adjustment=="None")
  label <- if (row$p>=.05) "ns" else strrep("*",sum(row$p<c(.05,.01,.001,.0001)))
  stopifnot(identical(row$label,label))
  if (label == "ns") "" else label
}
md <- result$meta
groups <- c("Control","LPS","TRALI")
red <- c("#FFF5F0","#FC9272","#DE2D26","#80000F")
feature_red <- c("#FEE5D9","#FB6A4A","#CB181D","#67000D")
sub_colors <- setNames(c("#B2182B","#167A78","#C68325","#695596","#427DA5","#A9B0B3"),levels(md$Subtype))
group_colors <- c(Control="#D7DADC",LPS="#A7C6DE",TRALI="#E4A79B")
ec <- plot_shell(features,ids,md)
ec <- SetAssayData(ec,assay="RNA",layer="data",new.data=shown)
ec[["umap"]] <- CreateDimReducObject(embeddings=result$umap,key="UMAP_",assay="RNA")
xr <- range(result$umap[,1]); yr <- range(result$umap[,2])
text_theme <- theme_classic(base_size=9) + theme(text=element_text(family="sans",colour="#20282A"),
  axis.text=element_text(colour="#20282A"),strip.background=element_blank(),strip.text=element_text(face="bold"),
  legend.background=element_blank(),plot.margin=margin(3,4,3,4))
map_theme <- theme_void(base_size=9) + theme(plot.title=element_text(size=9,hjust=.5,face="bold"),
  plot.margin=margin(2,3,2,3),strip.text=element_blank(),strip.background=element_blank(),
  legend.position="none")
letter <- function(tag) ggdraw()+draw_label(tag,x=.014,y=.44,hjust=0,size=14,fontface="bold")
block <- function(tag,body,h) plot_grid(letter(tag),body,ncol=1,rel_heights=c(.25,h))
note <- function(label,size=8) ggdraw()+draw_label(label,x=.5,y=.53,size=size,hjust=.5)
# A and B use the same identity colours and one legend.
maps <- lapply(groups,function(gr) {
  p <- SCP::CellDimPlot(ec,group.by="Subtype",cells=rownames(md)[md$Group==gr],reduction="umap",
    show_stat=FALSE,label=FALSE,palcolor=sub_colors,raster=TRUE,raster.dpi=c(850,850),pt.size=2.4,
    title=NULL,legend.position="none",theme_use="theme_void")
  p + coord_fixed(xlim=xr,ylim=yr,expand=TRUE) + map_theme + labs(title=gr,x=NULL,y=NULL)
})
a <- plot_grid(plotlist=maps,nrow=1)
b <- SCP::CellStatPlot(ec,stat.by="Subtype",group.by="Sample",stat_type="percent",plot_type="bar",
  position="stack",palcolor=sub_colors,bg_alpha=0,label=FALSE,NA_stat=FALSE,title=NULL,
  xlab=NULL,ylab="Recovered ECs (%)")
actual <- merge(b$data,composition,by=c("Sample","Subtype"))
stopifnot(nrow(actual)==nrow(composition),max(abs(actual$value*100-actual$percent))<1e-8)
b$data$Condition <- factor(md$Group[match(b$data$Sample,md$Sample)],groups)
for (i in seq_along(b$layers)) if (inherits(b$layers[[i]]$geom,"GeomCol")) {
  b$layers[[i]]$aes_params$colour <- "white"; b$layers[[i]]$aes_params$linewidth <- .13
}
b <- b+facet_grid(~Condition,scales="free_x",space="free_x")+text_theme+
  theme(legend.position="none",axis.ticks.x=element_blank(),panel.spacing.x=grid::unit(1.8,"mm"))+
  labs(title=NULL,subtitle=NULL,x=NULL,fill=NULL)
shared <- get_legend(b+guides(fill=guide_legend(nrow=1,byrow=TRUE))+
  theme(legend.position="bottom",legend.title=element_blank(),legend.text=element_text(size=8),
        legend.key.size=grid::unit(3,"mm"),legend.margin=margin(0,0,0,0)))
# C: dot diameter is detection; colour is mean log-normalized RNA including zeroes.
ec$All.groups <- md$Subtype
ComplexHeatmap::ht_opt(legend_title_gp=grid::gpar(fontsize=8,fontface="bold"),
                       legend_labels_gp=grid::gpar(fontsize=7.8))
markers$label <- markers$Subtype
hm <- dot_from_summary(markers,read.csv(file.path(input,"marker_order.csv"))$gene,levels(md$Subtype),assay="RNA",
  slot="data",aggregate_fun=mean,exp_method="raw",exp_cutoff=0,lib_normalize=FALSE,
  add_dot=TRUE,dot_size=grid::unit(4,"mm"),flip=TRUE,border=FALSE,nlabel=0,row_title="",column_title="",
  row_names_side="left",column_names_side="bottom",show_row_names=TRUE,show_column_names=TRUE,
  column_names_rot=52,cluster_rows=FALSE,cluster_columns=FALSE,use_raster=FALSE,exp_legend_title="Mean log RNA",
  heatmap_palcolor=red,limits=c(0,ceiling(max(markers$mean_log_RNA))),
  ht_params=list(row_names_gp=grid::gpar(fontsize=8),column_names_gp=grid::gpar(fontsize=8,fontface="italic")))
mat <- hm$matrix_list[[1]]
stopifnot(max(abs(mat[cbind(match(markers$gene,rownames(mat)),match(markers$Subtype,colnames(mat)))]-markers$mean_log_RNA))<1e-8)
# Give the SCP heatmap's two absolute-size legends their full device height
# before fitting it to the plate. Capture at print resolution.
capture_marker <- ragg::agg_capture(width=5.1,height=3.7,units="in",res=650,background="white")
print(hm$plot)
marker_img <- capture_marker()
invisible(grDevices::dev.off())
cplot <- ggdraw()+draw_grob(grid::rasterGrob(marker_img,interpolate=TRUE))
# D: the six reserved Figure 1 feature genes on the Figure 1 endothelial coordinates.
features <- c("Cdh5","Gpihbp1","Car4","Gja5","Slc6a2","Ccl21a")
fp <- SCP::FeatureDimPlot(ec,features=features,reduction="umap",assay="RNA",slot="data",raster=TRUE,
  raster.dpi=c(850,850),pt.size=2.5,bg_cutoff=0,bg_color="#E5E5E5",lower_cutoff=0,upper_cutoff=4,
  keep_scale="all",show_stat=FALSE,palcolor=feature_red,theme_use="theme_void",combine=FALSE)
fp <- lapply(seq_along(fp),function(i) {
  stopifnot(identical(unique(as.character(fp[[i]]$data$features)),features[i]))
  fp[[i]]+coord_fixed(xlim=xr,ylim=yr)+map_theme+labs(title=features[i],x=NULL,y=NULL)+
    theme(plot.title=element_text(face="italic",size=8,hjust=.5))
})
feat_key <- get_legend(fp[[1]]+guides(colour=guide_colourbar(title="log RNA",direction="horizontal",
  theme=theme(legend.key.width=grid::unit(19,"mm"),legend.key.height=grid::unit(2,"mm"))))+
  theme(legend.position="bottom",legend.title=element_text(size=8),legend.text=element_text(size=7),
        legend.margin=margin(0,0,0,0)))
d <- plot_grid(plot_grid(plotlist=fp,ncol=3),feat_key,ncol=1,rel_heights=c(1,.12))
# SCP mouse-level boxes for E, F and G. Each point keeps its observed y value.
mouse_box <- function(tab,panel_name,axis_label,upper,show_axis,stars=NULL,value_column="AUC") {
  tab$Group <- factor(tab$Group,groups)
  rownames(tab) <- as.character(tab$Sample)
  obs <- CreateSeuratObject(matrix(1,1,nrow(tab),dimnames=list("observation",rownames(tab))),meta.data=tab)
  p <- SCP::FeatureStatPlot(obs,stat.by=value_column,group.by="Group",plot_type="box",palcolor=group_colors,
    bg_alpha=0,add_point=TRUE,pt.color="#28363B",pt.size=1.45,jitter.height=0,jitter.width=.11,combine=FALSE)[[1]]
  p$layers <- Filter(function(l) inherits(l$geom,"GeomBoxplot") ||
    (inherits(l$geom,"GeomPoint") && inherits(l$stat,"StatIdentity")),p$layers)
  stopifnot(length(p$layers)==2L,nrow(p$data)==11L,setequal(p$data$cells,tab$Sample))
  for (i in seq_along(p$layers)) {
    if (inherits(p$layers[[i]]$geom,"GeomPoint")) {
      offsets <- ave(seq_len(nrow(p$data)),p$data$group.by,FUN=function(x) seq(-.13,.13,length.out=length(x)))
      p$layers[[i]]$position <- position_nudge(x=offsets)
    } else {
      p$layers[[i]]$aes_params$width <- .5
      p$layers[[i]]$aes_params$linewidth <- .36
      p$layers[[i]]$aes_params$outlier.shape <- NA
    }
  }
  p <- p+text_theme+coord_cartesian(ylim=c(0,upper))+labs(title=panel_name,x=NULL,y=if(show_axis)axis_label else NULL)+
    theme(legend.position="none",strip.text=element_blank(),strip.background=element_blank(),
      plot.title=element_text(hjust=.5,size=8.8,face="bold"),axis.ticks.x=element_blank(),axis.text.x=element_text(size=7.4),
      axis.title.y=element_text(family="sans",size=9,margin=margin(r=7)))
  if (!is.null(stars) && nzchar(stars)) {
    h <- max(tab[[value_column]])*1.095
    stopifnot(h+.08*upper<upper)
    plot_labels <<- c(plot_labels,stars)
    p <- p+annotate("segment",x=2,xend=3,y=h,yend=h,linewidth=.34)+
      annotate("segment",x=c(2,3),xend=c(2,3),y=h,yend=h-.025*upper,linewidth=.34)+
      annotate("text",x=2.5,y=h+.045*upper,label=stars,size=3.6)
  }
  b <- ggplot_build(p)$data
  j <- which(vapply(p$layers,function(l) inherits(l$geom,"GeomPoint"),logical(1)))
  expected <- tab[p$data$cells,,drop=FALSE]
  boxes <- which(vapply(p$layers,function(l) inherits(l$geom,"GeomBoxplot"),logical(1)))
  stopifnot(identical(as.character(p$data$group.by),as.character(expected$Group)),
            max(abs(p$data$value-expected[[value_column]]))<1e-12,
            length(boxes)==1L,nrow(b[[boxes]])==3L,nrow(b[[j]])==11L,
            max(abs(sort(b[[j]]$y)-sort(tab[[value_column]])))<1e-12)
  if (!is.null(stars) && nzchar(stars)) {
    label_layer <- which(vapply(p$layers,function(l) inherits(l$geom,"GeomText"),logical(1)))
    stopifnot(length(label_layer)==1L,identical(as.character(b[[label_layer]]$label),stars),
              all(b[[label_layer]]$y<upper))
  }
  p
}
# E. Each mouse contributes one original DESeq2-normalized pseudobulk value per population.
populations <- c("Total EC","gCap","aCap")
design_by_mouse <- unique(md[,c("Sample","Group")])
stopifnot(nrow(rna)==33L,!anyDuplicated(rna[,c("Sample","population")]),
          setequal(rna$population,populations),setequal(rna$Sample,design_by_mouse$Sample),
          all(rna$Group==as.character(design_by_mouse$Group[match(rna$Sample,design_by_mouse$Sample)])),
          all(rna$n_cells>0),all(rna$size_factor>0),
          max(abs(log2(1+rna$Fosl1_UMI/rna$size_factor)-rna$log2_normalized))<1e-10)
expression_boxes <- lapply(populations,function(pop) {
  tab <- rna[rna$population==pop,]
  mouse_box(tab,pop,"Fosl1 expression",
            max(rna$log2_normalized)*1.30,pop=="Total EC",stars=box_label("E",pop),value_column="log2_normalized")
})
e <- plot_grid(plotlist=expression_boxes,nrow=1,align="h",axis="tb")

lineage_boxes <- function(tab,axis_label,panel_name) {
  upper <- max(tab$AUC)*1.30
  plot_grid(mouse_box(subset(tab,compartment=="NonEC"),"Non-endothelium",axis_label,upper,TRUE,
                      stars=box_label(panel_name,"NonEC")),
            mouse_box(subset(tab,compartment=="EC"),"Endothelium",axis_label,upper,FALSE,
                      stars=box_label(panel_name,"EC")),
            nrow=1,align="h",axis="tb")
}
f <- lineage_boxes(lineage,"Fosl1 regulon activity","F")
# G. One 185-member signature, evaluated within each subtype.
# SCENIC returns a native blue channel. Display that channel on a pale-to-dark blue scale.
focal <- table_in("Fosl1_activity_points.csv.gz")
stopifnot(!anyDuplicated(focal$plot_row), all(focal$plot_row %in% seq_along(ids)))
au <- matrix(focal$AUC,1,dimnames=list("Fosl1",ids[focal$plot_row]))
rownames(au) <- "Fosl1 (185g)"
native <- new("aucellResults",SummarizedExperiment(assays=list(AUC=au)),nGenesDetected=numeric(0))
tmp_auc <- tempfile(fileext=".rds")
tmp_xy <- tempfile(fileext=".rds")
saveRDS(native,tmp_auc)
saveRDS(list(Y=result$umap[colnames(au),,drop=FALSE]),tmp_xy)
opt <- new("ScenicOptions",fileNames=list(int=data.frame(path=tmp_auc,row.names="aucell_regulonAUC")))
max_auc <- as.numeric(quantile(au,.99))
channel <- SCENIC::plotEmb_rgb(opt,regulonNames=list(blue="Fosl1"),aucType="AUC",
  aucMaxContrast=max_auc/max(au),offColor=NULL,showPlot=FALSE,showLegend=FALSE,tSNE_fileName=tmp_xy)
unlink(c(tmp_auc,tmp_xy))
native_intensity <- as.numeric(grDevices::col2rgb(channel)["blue",])/255
stopifnot(identical(names(channel),colnames(au)),max(abs(native_intensity-pmin(au[1,]/max_auc,1)))<=1/255)
xy <- data.frame(cell=colnames(au),X=result$umap[colnames(au),1],Y=result$umap[colnames(au),2],
              Group=md[colnames(au),"Group"],Subtype=md[colnames(au),"Subtype"],display_AUC=native_intensity*max_auc)
blue <- c("#E5F0F9","#A3CBE5","#4285BC","#13558F","#082D60")
xlim_g <- range(xy$X)+c(-.4,.4); ylim_g <- range(xy$Y)+c(-.4,.4)
cap_types <- c("gCap","aCap")
activity_map <- function(type,gr) {
  z <- xy[xy$Group==gr,]
  marked <- z[z$Subtype==type,]
  marked <- marked[order(marked$display_AUC),]
  ggplot(z,aes(X,Y))+geom_point(colour="#DEE0E2",size=.19,stroke=0)+
    geom_point(data=marked,aes(colour=display_AUC),size=.27,stroke=0)+
    scale_colour_gradientn(colours=blue,limits=c(0,max_auc),name="AUCell")+
    coord_fixed(xlim=xlim_g,ylim=ylim_g,expand=FALSE)+theme_void()+theme(legend.position="none",plot.margin=margin(0,0,0,0))
}
g_header <- plot_grid(NULL,plot_grid(plotlist=lapply(groups,note,size=8.3),nrow=1),nrow=1,rel_widths=c(.12,1))
g_rows <- lapply(cap_types,function(type) plot_grid(note(type,8.4),
  plot_grid(plotlist=lapply(groups,function(gr)activity_map(type,gr)),nrow=1),nrow=1,rel_widths=c(.12,1)))
bar <- seq(0,max_auc,length.out=256)
g_key <- ggplot(data.frame(x=bar,y=1),aes(x,y,fill=x))+geom_raster()+
  scale_fill_gradientn(colours=blue,limits=c(0,max_auc),guide="none")+
  scale_x_continuous(breaks=c(0,max_auc/2,max_auc),labels=c("0",sprintf("%.2f",max_auc/2),paste0("\u2265",sprintf("%.2f",max_auc))),expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+theme_void()+theme(axis.text.x=element_text(size=6.5))
g_maps <- plot_grid(g_header,g_rows[[1]],g_rows[[2]],plot_grid(NULL,note("AUCell",7.5),g_key,NULL,nrow=1,rel_widths=c(.9,.35,1,.9)),
  ncol=1,rel_heights=c(.14,1,1,.22))
gb <- lapply(cap_types,function(type) {
  q <- effects$BH_within[effects$Subtype==type]
  stopifnot(length(q)==1,is.finite(q))
  mouse_box(subset(cap,subtype==type),type,"Fosl1 regulon activity",max(cap$AUC)*1.30,type=="gCap",box_label("G",type))
})
g <- plot_grid(g_maps,plot_grid(plotlist=gb,nrow=1,align="h",axis="tb"),nrow=1,rel_widths=c(.56,.44))

# H. Native report selection: five largest RSS ranks per subtype, plus Fosl1.
ht <- table_in("de_novo_RSS_activity_RNA.csv")
types <- c("gCap","aCap","Arterial","Pulmonary venous","Lymphatic")
top5 <- lapply(types,function(ct) {v <- ht[ht$Subtype==ct,];head(v$TF[order(v$Rank,v$TF)],5)})
picked <- unique(c(unlist(top5),"Fosl1"))
ht <- ht[ht$TF %in% picked,]
ht$TF <- factor(ht$TF,rev(picked))
ht$Subtype <- factor(ht$Subtype,types)
tf_lab <- function(x) as.expression(lapply(x,function(g) if(g=="Fosl1") bquote(bolditalic(.(g))) else bquote(italic(.(g)))))
heatmap <- function(metric,title,with_dots) {
  p <- ggplot(ht,aes(Subtype,TF))+geom_tile(aes(fill=.data[[metric]]),colour="white",linewidth=.18)+
    scale_fill_gradient2(low="#3273A4",mid="#F8F7F4",high="#BA513E",limits=c(-2,2),breaks=c(-2,0,2),name="Z score")+
    scale_y_discrete(labels=tf_lab,expand=c(0,0))+
    scale_x_discrete(labels=c("gCap","aCap","Arterial","P. venous","Lymphatic"),expand=c(0,0))+
    text_theme+labs(title=title,x=NULL,y=NULL)+theme(axis.line=element_blank(),axis.ticks=element_blank(),
      axis.text.x=element_text(size=7,angle=40,hjust=1),axis.text.y=element_text(size=7),
      plot.title=element_text(hjust=.5,size=8.6,face="bold"),legend.position="none")
  if(with_dots) p <- p+geom_point(aes(size=RSS),colour="#25303A")+
    scale_size_area(name="RSS",max_size=2.25,limits=c(0,max(ht$RSS)),breaks=c(.2,.4,.6))
  p
}
ha <- heatmap("Z_AUC","Regulon activity",TRUE)
hr <- heatmap("Z_RNA","TF expression",FALSE)
hk <- get_legend(ha+guides(fill=guide_colourbar(barwidth=grid::unit(13,"mm"),barheight=grid::unit(2,"mm")),
  size=guide_legend(nrow=1))+theme(legend.position="bottom",legend.direction="horizontal",legend.text=element_text(size=6.8),
  legend.title=element_text(size=7),legend.key.width=grid::unit(3.2,"mm"),legend.spacing.x=grid::unit(1,"mm"),legend.margin=margin(0,0,0,0)))
h <- plot_grid(plot_grid(ha,hr,nrow=1,align="h",axis="tb"),hk,ncol=1,rel_heights=c(1,.085))

# I. Use identical x and size scales in the two top/bottom-six lollipops.
tt <- table_in("de_novo_TRALI_vs_LPS_all_TFs.csv")
bound <- max(abs(tt$DeltaAUC))*1.07
evidence_max <- max(-log10(tt$BH_within))
lollipop <- function(type) {
  v <- tt[tt$Subtype==type,]
  up <- head(v$TF[v$DeltaAUC>0][order(-v$DeltaAUC[v$DeltaAUC>0])],6)
  down <- head(v$TF[v$DeltaAUC<0][order(v$DeltaAUC[v$DeltaAUC<0])],6)
  d <- v[v$TF %in% unique(c(up,down,"Fosl1")),]
  d$TF <- factor(d$TF,d$TF[order(d$DeltaAUC)])
  d$direction <- ifelse(d$DeltaAUC>0,"Up","Down")
  d$evidence <- -log10(d$BH_within)
  ggplot(d,aes(DeltaAUC,TF))+geom_vline(xintercept=0,colour="#8C979A",linewidth=.34)+
    geom_segment(aes(x=0,xend=DeltaAUC,yend=TF,colour=direction),linewidth=.65)+
    geom_point(aes(size=evidence,fill=direction,shape=BH_within<.05),colour="#25303A",stroke=.42)+
    scale_shape_manual(values=c("TRUE"=21,"FALSE"=1),guide="none")+
    scale_colour_manual(values=c(Up="#B44936",Down="#3F739F"),guide="none")+
    scale_fill_manual(values=c(Up="#B44936",Down="#3F739F"),guide="none")+
    scale_size_continuous(name=expression(-log[10]("BH q")),limits=c(0,evidence_max),range=c(1.4,3.9))+
    scale_y_discrete(labels=tf_lab)+scale_x_continuous(limits=c(-bound,bound),breaks=scales::breaks_pretty(3))+
    labs(title=type,x="Change in regulon activity",y=NULL)+text_theme+
    theme(axis.ticks.y=element_blank(),axis.text.y=element_text(size=7.6),axis.text.x=element_text(size=7),
      axis.title.x=element_text(size=7.4),legend.position="none",plot.title=element_text(hjust=.5,size=8.6,face="bold"))
}
lp <- lapply(cap_types,lollipop)
lk <- get_legend(lp[[1]]+theme(legend.position="bottom",legend.title=element_text(size=7),
  legend.text=element_text(size=7),legend.margin=margin(0,0,0,0)))
i_panel <- plot_grid(plot_grid(plotlist=lp,nrow=1),lk,ncol=1,rel_heights=c(1,.085))
top <- plot_grid(block("A",a,1.55),block("B",b,1.55),nrow=1,rel_widths=c(.54,.46))
middle <- plot_grid(block("C",cplot,2.25),block("D",d,2.25),nrow=1,rel_widths=c(.48,.52))
ef <- plot_grid(block("E",e,1.66),block("F",f,1.66),nrow=1,rel_widths=c(.58,.42))
overview <- plot_grid(block("H",h,3.28),block("I",i_panel,3.28),nrow=1,rel_widths=c(.50,.50))
fig <- plot_grid(top,shared,middle,ef,block("G",g,2.30),overview,ncol=1,
                 rel_heights=c(1.75,.27,2.30,1.91,2.55,3.53))
fig <- ggdraw()+draw_plot(fig,.012,.008,.976,.984)
for (fmt in c("pdf","png","tiff")) {
  dest <- file.path(out,paste0("Figure2.",fmt))
  if (fmt=="pdf") ggsave(dest,fig,width=10.8,height=14.5,device=cairo_pdf,bg="white")
  if (fmt=="png") ggsave(dest,fig,width=10.8,height=14.5,device=ragg::agg_png,dpi=200,bg="white")
  if (fmt=="tiff") ggsave(dest,fig,width=10.8,height=14.5,device=ragg::agg_tiff,dpi=600,compression="lzw",bg="white")
}
sensitivity <- table_in("sensitivity_171/lineage_AUC_by_mouse.csv")
stopifnot(nrow(sensitivity)==22L,identical(unique(sensitivity$signature),171L))
sensitivity_plot <- plot_grid(note("171-gene sensitivity analysis",11),
                              lineage_boxes(sensitivity,"Mean 171-gene AUCell","sensitivity_171"),
                              ncol=1,rel_heights=c(.12,1))
dir.create(file.path(out,"sensitivity_171"),showWarnings=FALSE)
ggsave(file.path(out,"sensitivity_171","Fosl1_171_by_mouse.png"),sensitivity_plot,width=5,height=2.8,
       device=ragg::agg_png,dpi=400,bg="white")
expected_order <- c(box_label("E","Total EC"),box_label("E","gCap"),box_label("E","aCap"),
                    box_label("F","NonEC"),box_label("F","EC"),box_label("G","gCap"),box_label("G","aCap"),
                    box_label("sensitivity_171","NonEC"),box_label("sensitivity_171","EC"))
stopifnot(identical(plot_labels,expected_order[nzchar(expected_order)]))
annotation_order <- c(paste("E",populations),paste("F",c("NonEC","EC")),
                      paste("G",cap_types),paste("sensitivity_171",c("NonEC","EC")))
annotation_audit <- annotations[match(annotation_order,paste(annotations$panel,annotations$population)),]
annotation_audit$rendered_label <- expected_order
annotation_audit$bracket_drawn <- nzchar(expected_order)
write.csv(annotation_audit,file.path(out,"boxplot_annotation_audit.csv"),row.names=FALSE)
