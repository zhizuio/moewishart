rm(list = ls())

# user-defined directory to save simulation results as data files
file_directory <- "../moewishart/"
filename_save <- file_directory


# NOTE: manually copy the cluster indices from four methods into the CSV file `metaData_selected.csv` and save it as `metaData_selected2.csv`

# load CTRP drug sensitivity's AUC data
aucData <- read.csv(paste0(filename_save, "v20.data.curves_post_qc.txt"), header=T, sep="\t", fill=T)

load(paste0(filename_save, "CTRP_viability.RData"))
length(unique(viability_selected$master_cpd_id))
auc <- aucData[aucData$master_cpd_id %in% unique(viability_selected$master_cpd_id), ]
auc <- auc[auc$experiment_id %in% unique(viability_selected$experiment_id), ]

auc2 <- auc[, c(1,17,16)]
auc3 <- reshape(auc2, v.names="area_under_curve", timevar="master_cpd_id", idvar="experiment_id", direction="wide")
auc3[1:5,1:4]
length(unique(auc3$experiment_id))
sum(is.na(auc3))
colSums(is.na(auc3))


#===============
# select nonmissing AUC data
#===============
y0 <- auc3[, -1]
m0 <- dim(y0)[2]-1
eps <- 0.05
# r1.na is better to be not smaller than r2.na
r1.na <- 0.3
r2.na <- 0.2
k <- 1
while(sum(is.na(y0[,2:(1+m0)]))>0){
  r1.na <- r1.na - eps/k
  r2.na <- r1.na - eps/k
  k <- k + 1
  ## select drugs with <30% (decreasing with k) missing data overall cell lines
  na.y <- apply(y0[,2:(1+m0)], 2, function(xx) sum(is.na(xx))/length(xx))
  while(sum(na.y<r1.na)<m0){
    y0 <- y0[,-c(1+which(na.y>=r1.na))]
    m0 <- sum(na.y<r1.na)
    na.y <- apply(y0[,2:(1+m0)], 2, function(xx) sum(is.na(xx))/length(xx))
  }
  
  ## select cell lines with treatment of at least 80% (increasing with k) drugs
  na.y0 <- apply(y0[,2:(1+m0)], 1, function(xx) sum(is.na(xx))/length(xx))
  while(sum(na.y0<r2.na)<(dim(y0)[1])){
    y0 <- y0[na.y0<r2.na,]
    na.y0 <- apply(y0[,2:(1+m0)], 1, function(xx) sum(is.na(xx))/length(xx))
  }
  num.na <- sum(is.na(y0[,2:(1+m0)]))
  message("#{NA}=", num.na, "\n", "r1.na =", r1.na, ", r2.na =", r2.na, "\n")
}

dim(y0)
master_cpd_id172 <- substr(colnames(y0), 18, nchar(colnames(y0)))
y1 <- y0
colnames(y1) <- master_cpd_id172

# load meta data together with four methods' grouping results
metaData_selected_grouped <- read.csv(paste0(filename_save, "metaData_selected2.csv"), header=T, sep=";", fill=T)
metaData_selected_grouped172 <- metaData_selected_grouped[metaData_selected_grouped$master_cpd_id %in% colnames(y1), ]

y2 <- y1[, match(metaData_selected_grouped172$master_cpd_id, colnames(y1))]
sum(colnames(y2) == metaData_selected_grouped172$master_cpd_id)


metaData_grouped172 <- metaData_selected_grouped172[, c(2, 3, 13:18)]
length(unique(metaData_grouped172$EM_K8))

names(metaData_grouped172)[names(metaData_grouped172) == "combined_group"] <- "MoA"
metaData_grouped172$MoA[metaData_grouped172$MoA == "Miscellaneous or screening hit/natural product"] <- "Miscellaneous/screening hit"
metaData_grouped172$MoA <- factor(metaData_grouped172$MoA)

## plot heatmap with annotations
library(ComplexHeatmap)
anno <- metaData_grouped172
anno$status <- factor(anno$status, levels = c("approved_or_clinical", "experimental"),
  labels = c("FDA approved or clinical", "experimental"))
names(anno)[3:7] <- c("Status", "EM", "EM-MoE", "Bayes", "Bayes-MoE")
rownames(anno) <- anno$master_cpd_id
anno <- anno[, -c(1:2)]

## normalize AUC data
y3 <- scale(y2)
rownames(y3) <- paste0("cell_line", 1:nrow(y3))

# Define distinct color sets
status_colors = c("FDA approved or clinical" = "#E41A1C", "experimental" = "#377EB8")

combined_levels <- levels(as.factor(anno$MoA)) 
combined_colors <- setNames(ggsci::pal_igv()(11), combined_levels)
names(combined_colors) <- combined_levels

# Define Shared Colors 
shared_levels <- 1:10 
shared_colors <- structure(RColorBrewer::brewer.pal(10, "Paired"), names = shared_levels)

ordered_cols <- c("Status", "MoA", "Bayes", "Bayes-MoE", "EM", "EM-MoE")
anno <- anno[, ordered_cols]
ann_cols <- list(Status = status_colors, MoA = combined_colors)
for(name in setdiff(colnames(anno), c("Status", "MoA"))) {
    ann_cols[[name]] = shared_colors
}

# Define heatmap
col_ann <- HeatmapAnnotation(df = anno, col = ann_cols, show_legend = FALSE)
col_fun = colorRamp2(c(min(y3), 0, max(y3)), c("blue", "white", "red"))

ht <- Heatmap(y3, 
              name = "Z-Score", 
              col = col_fun,
              top_annotation = col_ann,
              show_row_names = FALSE, 
              show_column_names = FALSE,
              cluster_rows = FALSE,      
              show_heatmap_legend = FALSE, 
              column_title = "Compounds", 
              column_title_side = "bottom",
              row_title = "Cell Lines")

# Legends
lgd_zscore = Legend(title = "Standardized\nAUC", col_fun = col_fun, 
                    title_gp = gpar(fontsize = 8, fontface = "bold"),
                    legend_height = unit(2.5, "cm"))
lgd_status = Legend(title = "Status", at = names(status_colors), 
                    legend_gp = gpar(fill = status_colors),
                    grid_height = unit(3, "mm"), grid_width = unit(3, "mm"),
                    title_gp = gpar(fontsize = 11, fontface = "bold"),
                    labels_gp = gpar(fontsize = 10))
lgd_moa = Legend(title = "MoA", at = names(combined_colors), 
                      legend_gp = gpar(fill = combined_colors),
                      ncol = 1, # Split into two columns of 7
                      grid_height = unit(3, "mm"), grid_width = unit(3, "mm"),
                      title_gp = gpar(fontsize = 9, fontface = "bold"),
                      labels_gp = gpar(fontsize = 8))
lgd_shared = Legend(title = "Component index", at = names(shared_colors), 
                    legend_gp = gpar(fill = shared_colors),
                    grid_height = unit(3, "mm"), grid_width = unit(3, "mm"),
                    title_gp = gpar(fontsize = 11, fontface = "bold"),
                    labels_gp = gpar(fontsize = 10),
                    ncol = 1) # Set to 2 if too tall


pdf("moewishart_CTRP_heatmap.pdf",width=6, height=6)

# Draw heatmap with wide right padding for 3 columns of legends
draw(ht, padding = unit(c(2, 2, 2, 53), "mm")) 
draw(lgd_zscore, x = unit(0.52, "npc"), y = unit(0.55, "npc"), just = "left")
draw(lgd_status, x = unit(0.7, "npc"), y = unit(0.83, "npc"), just = "left")
draw(lgd_moa, x = unit(0.52, "npc"), y = unit(0.23, "npc"), just = "left")
draw(lgd_shared, x = unit(0.7, "npc"), y = unit(0.59, "npc"), just = "left")

dev.off()