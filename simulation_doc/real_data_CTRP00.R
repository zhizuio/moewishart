##===================
## Analyze data for paper by The Tien Mai & Zhi Zhao
##===================

#===============
# The user needs load the datasets from https://www.cancerrxgene.org/downloads/bulk_download archived data folder 'release-5.0'. 
# Downloading the three datasets used for our analysis.
#===============

auc0 <- read.table("/Users/zhiz/Downloads/IPFStructPenalty/CTRPv2.0_2015_ctd2_ExpandedDataset/v20.data.curves_post_qc.txt", header=T)
auc0 <- auc0[,c(1,17,16)]

compound <- read.table("/Users/zhiz/Downloads/IPFStructPenalty/CTRPv2.0_2015_ctd2_ExpandedDataset/v20.meta.per_compound.txt", fill=T)
compound <- compound[,1:2]
colnames(compound) <- compound[1,]
auc1 <- merge(auc0, compound, by="master_cpd_id")

experiment <- read.table("/Users/zhiz/Downloads/IPFStructPenalty/CTRPv2.0_2015_ctd2_ExpandedDataset/v20.meta.per_experiment.txt", header=T)
experiment <- experiment[,c(1,9)]
auc2 <- merge(auc1, experiment, by="experiment_id")


cell_line <- read.table("/Users/zhiz/Downloads/IPFStructPenalty/CTRPv2.0_2015_ctd2_ExpandedDataset/v20.meta.per_cell_line.txt", fill=T)
cell_line <- cell_line[,c(1,2)]
colnames(cell_line) <- cell_line[1,]
cell_line <- cell_line[-1,]
auc3 <- merge(auc2, cell_line, by="master_ccl_id")
