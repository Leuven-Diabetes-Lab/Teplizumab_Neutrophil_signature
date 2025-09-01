setwd("D:/aCD3_signature/Test_folder")

#Loading libraries:
library("DESeq2")
library("dplyr")
library("readxl")
library("Hmisc")
library("ggpubr")

#####Preparing data for correlation calculation########

#Download metadata: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA338797
#reading file with all AUC:
auc_file <- read_excel("SraRunTable_AbaTE_ITN027AI.xlsx", sheet= "SraRunTable")


#reading file with normalized gene expression all months

gene_expr_all <- read.csv("global.normalized_counts_with_controls_0_filtered.csv")



#making X column character
gene_expr_all$X <- as.character(gene_expr_all$X)
#renaming first column
names(gene_expr_all)[names(gene_expr_all) == 'X'] <- 'Symbol'
# making first column the index
rownames(gene_expr_all) <- gene_expr_all$Symbol
#remove original ID column from data frame
gene_expr_all$Symbol <- NULL

#selecting timepoint 0 participants
auc_0 <- subset(auc_file, auc_file$visitmonth == "0", select = c("Run","status","auc_percent_of_baseline","itn_name", "visitmonth"))
auc_0 = data.frame(auc_0)

#The Sra run table is missign some auc from visit 12month.
#You can download the full table at https://www.itntrialshare.org/study/Studies/ITN027AIPUBLIC/Study%20Data/dataset.view?datasetId=1201
#Here this file is named "Extra_info.xlsx"
extra_info <- read_excel("Extra_info.xlsx")

#selecting timepoint 12 from new extra info file to extract the auc values
extra_info_12 <- subset(extra_info, extra_info$`Visit Month` == c("12"), select = c("Participant ID","Visit Month","AUC Percent of Baseline"))
extra_info_12 <- data.frame(extra_info_12)

#left join old file with new file subsets
merge_auc_extra <-merge(x = auc_0, y = extra_info_12, by.x = "itn_name",  by.y = "Participant.ID", all.x = TRUE)

#make run column index (this is run from timepoint 0)
rownames(merge_auc_extra) <- merge_auc_extra$Run

#Seems like AbATE_693516 is missing AUC info for month 12, it is present in SraRunTable_AbaTE_ITN027AI.xlsx table.
#Manually adding the missing value:
merge_auc_extra[merge_auc_extra$itn_name == 'AbATE_693516', 'AUC.Percent.of.Baseline'] = 0

#####only run this line if want to delete controls##############
merge_auc_extra <- merge_auc_extra[!(merge_auc_extra$status %in% "C"),]
################################################################

#Make subsets for gene expression
gene_expr_0 <- subset(gene_expr_all, select = c(row.names(merge_auc_extra)))

#FILTER GENE EXPRESSION: they have to be expressed at least 10 times in minimum 60% of the patients
#transposing it too so that observations are rows and variables are columns

# smallest_count <- 10
# keep <- rowSums(gene_expr_0 >= smallest_count) >= ceiling(0.6 * ncol(gene_expr_0))
# gene_expr_subset_0 <- gene_expr_0[keep,]
# tgene_expr_0 <- data.frame(t(gene_expr_subset_0))

tgene_expr_0 <- data.frame(t(gene_expr_0))

#merge gex and auc
tgene_expr_all.auc_0_12 = merge(merge_auc_extra, tgene_expr_0, by = "row.names")
row.names(tgene_expr_all.auc_0_12) = tgene_expr_all.auc_0_12[,1]
tgene_expr_all.auc_0_12 = tgene_expr_all.auc_0_12[,-(1:7)]
tgene_expr_all.auc_0_12$itn_name <- NULL

#correlation with p-values
#making correlation matrix
cor.mat_0_12_with_p <- rcorr(as.matrix(tgene_expr_all.auc_0_12), type = c("pearson"))


#Accessing the two layers in the object and saving them separately then merging it into one readable table
d_p_corr_0_12 = data.frame(row.names(cor.mat_0_12_with_p$r), AUC = cor.mat_0_12_with_p$r[,"AUC.Percent.of.Baseline"])
d_p_sign_0_12 = data.frame(row.names(cor.mat_0_12_with_p$P), AUC = cor.mat_0_12_with_p$P[,"AUC.Percent.of.Baseline"])
merge_corr_sign_0_12 <-  merge(d_p_corr_0_12, d_p_sign_0_12, by = "row.names")
row.names(merge_corr_sign_0_12) = merge_corr_sign_0_12[,1]
merge_corr_sign_0_12 <- merge_corr_sign_0_12[,c(1,3,5)]
# assigning new names to the columns of the data frame 
colnames(merge_corr_sign_0_12) <- c('Symbol','Pearson_corr','p_value')
merge_corr_sign_0_12 = merge_corr_sign_0_12[order(-merge_corr_sign_0_12$Pearson_corr),]

#fixing gene names
merge_corr_sign_0_12 <- merge_corr_sign_0_12[-1,-1]
rownames(merge_corr_sign_0_12) <- sub("[.]", "-", rownames(merge_corr_sign_0_12))
merge_corr_sign_0_12$Symbol <- rownames(merge_corr_sign_0_12)

#Pearson correlation table for Supplemental Table 6
write.csv(merge_corr_sign_0_12, "AbATE_corr_sign_Baseline.csv")

#selecting only the significant genes to make the AbATE Baseline gene signature
Abate_signature<-merge_corr_sign_0_12[merge_corr_sign_0_12$p_value<(0.01),] 

#R signature - Baseline AbATE table for Supplemental Table 6
Abate_signature_R<-Abate_signature[Abate_signature$Pearson_corr >0,] 
write.csv(Abate_signature_R, "AbATE_Baseline_R_signature.csv")

#NR signature - Baseline AbATE table for Supplemental Table 6
Abate_signature_NR<-Abate_signature[Abate_signature$Pearson_corr <0,] 
write.csv(Abate_signature_NR, "AbATE_Baseline_NR_signature.csv")

#####Loading Human scRNAseq data########
library(Seurat)
library(dplyr)
library(ggplot2)
library(future)
library(SCpubr)

# check the current active plan
plan()
# change the current plan to access parallelization
plan("multisession", workers = 4)
options(future.globals.maxSize = 2e10)
future.seed=TRUE
plan()

WB_NO <-readRDS("./scRNAseq/WB_scRNAseq.rds")

Idents(WB_NO)<-"Annotation"


#Calculate the average expression levels of AbATE Baseline gene signature for each cell
AbATE_Baseline_Signature <- list(AbATE_Baseline_R_sign = na.omit(Abate_signature_R$Symbol),
                               AbATE_Baseline_NR_sign = na.omit(Abate_signature_NR$Symbol)
)

WB_NO <- AddModuleScore(WB_NO, AbATE_Baseline_Signature[1], name = names(AbATE_Baseline_Signature[1]))
WB_NO <- AddModuleScore(WB_NO, AbATE_Baseline_Signature[2], name = names(AbATE_Baseline_Signature[2]))


#Feature plot of the  average expression levels of AbATE Baseline gene signature with SCpubr
p_AbATE_Baseline_R_sign <- do_FeaturePlot(sample = WB_NO, features = paste0(names(AbATE_Baseline_Signature[1]),"1"), enforce_symmetry = F, diverging.palette =  "RdYlBu",
                                          min.cutoff = -0.06,
                                          max.cutoff = 0.13,
                                          plot_cell_borders = FALSE, use_viridis = TRUE, viridis.palette ="viridis", pt.size = 2,
                                          raster = TRUE, raster.dpi = 1024)

p_AbATE_Baseline_NR_sign <-do_FeaturePlot(sample = WB_NO, features = paste0(names(AbATE_Baseline_Signature[2]),"1"), enforce_symmetry = F, 
                                        min.cutoff = -0.06, 
                                        max.cutoff = 0.5,
                                        plot_cell_borders = FALSE, use_viridis = TRUE, viridis.palette ="viridis", pt.size = 2,
                                        raster = TRUE, raster.dpi = 1024)

ggsave(file="Fig5A.svg", plot=p_AbATE_Baseline_R_sign, width=8, height=7)
ggsave(file="Fig5B.svg", plot=p_AbATE_Baseline_NR_sign, width=8, height=7)

#We will now calculate the Average expression of the R and NR signature for each cell and plot in a Violin plot
#R signature
sub_WB_NO <- subset(WB_NO, features = AbATE_Baseline_Signature[[1]])
expression_matrix_subset_R <- GetAssayData(object = sub_WB_NO, layer = "data")
expression_matrix_subset_R_linear <- expm1(x = expression_matrix_subset_R) 
per_cell_average_R_score_linear <- colMeans(expression_matrix_subset_R_linear) # Average on linear scale
WB_NO$AbATE_Baseline_R_sign <- per_cell_average_R_score_linear[colnames(WB_NO)]

#NR signature
sub_WB_NO <- subset(WB_NO, features = AbATE_Baseline_Signature[[2]])
expression_matrix_subset_NR <- GetAssayData(object = sub_WB_NO, layer = "data")
expression_matrix_subset_NR_linear <- expm1(x = expression_matrix_subset_NR) 
per_cell_average_NR_score_linear <- colMeans(expression_matrix_subset_NR_linear) # Average on linear scale
WB_NO$AbATE_Baseline_NR_sign <- per_cell_average_NR_score_linear[colnames(WB_NO)]

#Violin plots
#Setup the color to match the cluster colors from Fig4B
library(forcats)
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
n = 13
cols = gg_color_hue(n)

#Violin plot of the R signature
feature_to_plot <- "AbATE_Baseline_R_sign"
grouping_variable <- "Annotation"

#Get the needed data
plot_data <- FetchData(WB_NO, vars = c(feature_to_plot, grouping_variable))
colnames(plot_data) <- c("Expression", "Group")
names(cols)<- levels(plot_data$Group)

#Manually reorganized to match fig5E order
desired_order <- c("NK/CD8 T cells","Plasma cells","DCs","T cells","Monocytes","B cells","Neutrophils 0","Basophils","Neutrophils 1", "Megakaryocytes","Neutrophils 2","Neutrophils 3","Neutrophils 4"  )
plot_data$Group <- fct_relevel(plot_data$Group, desired_order)

vln_plot_expr_R_sign<-ggplot(plot_data, aes(x = Group, y = Expression, fill = Group)) +
  geom_violin(trim = FALSE, # Set to TRUE to trim tails at data extremes
              alpha = 0.7,   # Transparency of the violin
              scale = "width" # "area", "count", or "width" - "width" makes all violins same max width
  ) +
  geom_point(position = position_jitter(width = 0.2, height = 0), # Jitter points horizontally
             size = 0.1,                                         # **Crucially, set a smaller dot size here**
             alpha = 0.1,                                        # Transparency of the points
             color = "black"                                     # Color of the points
  ) +
  # geom_boxplot(width = 0.1, # Add a small boxplot inside the violin
  #              outlier.shape = NA, # Hide outlier points from the boxplot as jittered points are shown
  #              color = "white",   # Color of the boxplot lines
  #              alpha = 0.5         # Transparency of the boxplot
  # ) +
  stat_summary(fun = mean,
               geom = "crossbar", # Horizontal line with optional width
               width = 0.2,
               color = "red",
               lwd = 0.5) +
  scale_fill_manual(
    values = cols, 
    breaks = desired_order
  )  + 
  labs(title = NULL,
       x = NULL,
       y = "Average expression") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 19, color = "black"),
        plot.title = element_text(size = 25),
        axis.title.y = element_text(size = 25, margin = margin(r = 33)),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 15),
        legend.position = "none")


ggsave(file="Fig5C.png", plot=vln_plot_expr_R_sign, width=6, height=6, dpi = 300)

#Violin plot of the NR signature
feature_to_plot <- "AbATE_Baseline_NR_sign"
grouping_variable <- "Annotation"

#Get the needed data
plot_data <- FetchData(WB_NO, vars = c(feature_to_plot, grouping_variable))
colnames(plot_data) <- c("Expression", "Group")
names(cols)<- levels(plot_data$Group)

#Manually reorganized to match fig6F order
desired_order <- c("Neutrophils 3","Neutrophils 1","Neutrophils 4","Neutrophils 2", "Megakaryocytes", "Basophils", "Neutrophils 0", "Monocytes", "NK/CD8 T cells","B cells","T cells","DCs","Plasma cells")
plot_data$Group <- fct_relevel(plot_data$Group, desired_order)

vln_plot_expr_NR_sign<-ggplot(plot_data, aes(x = Group, y = Expression, fill = Group)) +
  geom_violin(trim = FALSE, # Set to TRUE to trim tails at data extremes
              alpha = 0.7,   # Transparency of the violin
              scale = "width" # "area", "count", or "width" - "width" makes all violins same max width
  ) +
  geom_point(position = position_jitter(width = 0.2, height = 0), # Jitter points horizontally
             size = 0.1,                                         # **Crucially, set a smaller dot size here**
             alpha = 0.1,                                        # Transparency of the points
             color = "black"                                     # Color of the points
  ) +
  # geom_boxplot(width = 0.1, # Add a small boxplot inside the violin
  #              outlier.shape = NA, # Hide outlier points from the boxplot as jittered points are shown
  #              color = "white",   # Color of the boxplot lines
  #              alpha = 0.5         # Transparency of the boxplot
  # ) +
  stat_summary(fun = mean,
               geom = "crossbar", # Horizontal line with optional width
               width = 0.2,
               color = "red",
               lwd = 0.5) +
  scale_fill_manual(
    values = cols, 
    breaks = desired_order
  )  + 
  labs(title = NULL,
       x = NULL,
       y = "Average expression") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 19, color = "black"),
        plot.title = element_text(size = 25),
        axis.title.y = element_text(size = 25, margin = margin(r = 33)),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 15),
        legend.position = "none")

ggsave(file="Fig5D.png", plot=vln_plot_expr_NR_sign, width=6, height=6, dpi = 300)



#Let's plot pseudobulk average expression per cell type

#Plot average expression value per celltype for Positive correlation
R_avg_expression <- AverageExpression(WB_NO, features = AbATE_Baseline_Signature[1],group.by = "Annotation" )
R_avg_expression <- R_avg_expression$RNA

mean_expression_per_celltype <- colMeans(R_avg_expression)
mean_expression_per_celltype <- as.data.frame(mean_expression_per_celltype)
mean_expression_per_celltype$cell_type <- rownames(mean_expression_per_celltype)
colnames(mean_expression_per_celltype)[colnames(mean_expression_per_celltype) == "mean_expression_per_celltype"] <- "Mean_Pearson_R_sign"


ranked_cell_types <- mean_expression_per_celltype %>%
  arrange(desc(Mean_Pearson_R_sign))
ranked_cell_types$cell_type <- factor(ranked_cell_types$cell_type, 
                                      levels = ranked_cell_types$cell_type)


plot_expr_R_sign <- ggplot(ranked_cell_types, aes(x = cell_type, y = Mean_Pearson_R_sign)) +
  geom_point(size = 4) +
  labs(x = "Cell Type", y = "Average expression") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 20, color = "black"),
        plot.title = element_text(size = 25),
        axis.title.y = element_text(size = 25, margin = margin(r = 20)),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 15)) +
  labs(title = "R signature - AbATE Baseline")

plot_expr_R_sign

ggsave(file="Fig5E.svg", plot=plot_expr_R_sign, width=6, height=6)

#Plot average expression value per celltype for Positive correlation
NR_avg_expression <- AverageExpression(WB_NO, features = AbATE_Baseline_Signature[2],group.by = "Annotation" )
NR_avg_expression <- NR_avg_expression$RNA

mean_expression_per_celltype <- colMeans(NR_avg_expression)
mean_expression_per_celltype <- as.data.frame(mean_expression_per_celltype)
mean_expression_per_celltype$cell_type <- rownames(mean_expression_per_celltype)
colnames(mean_expression_per_celltype)[colnames(mean_expression_per_celltype) == "mean_expression_per_celltype"] <- "Mean_Pearson_NR_sign"


ranked_cell_types <- mean_expression_per_celltype %>%
  arrange(desc(Mean_Pearson_NR_sign))
ranked_cell_types$cell_type <- factor(ranked_cell_types$cell_type, 
                                      levels = ranked_cell_types$cell_type)


plot_expr_NR_sign <- ggplot(ranked_cell_types, aes(x = cell_type, y = Mean_Pearson_NR_sign)) +
  geom_point(size = 4) +
  labs(x = "Cell Type", y = "Average expression") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 20, color = "black"),
        plot.title = element_text(size = 25),
        axis.title.y = element_text(size = 25, margin = margin(r = 20)),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 15)) +
  labs(title = "NR signature - AbATE Baseline")

plot_expr_NR_sign

ggsave(file="Fig5F.svg", plot=plot_expr_NR_sign, width=6, height=6)


