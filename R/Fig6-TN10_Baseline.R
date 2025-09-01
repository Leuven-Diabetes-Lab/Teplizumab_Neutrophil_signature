setwd("D:/aCD3_signature/Test_folder")

#Loading libraries:
library("DESeq2")
library("dplyr")
library("readxl")
library("Hmisc")
library("ggpubr")

dds<-readRDS("TN10_WB_RNAseq.rds")

#Select for Baseline gene expression
to_keep<-which(dds$Time.Point %in% "Baseline")
dds_baseline <- dds[,to_keep]

#Separate Placebo from Teplizumab
to_keep<-which(dds_baseline$TREATMENTDESC %in% "Placebo")
dds_baseline_Placebo <- dds_baseline[,to_keep]

to_keep<-which(dds_baseline$TREATMENTDESC %in% "Teplizumab")
dds_baseline_Teplizumab <- dds_baseline[,to_keep]


#Make Matrix with 1st column is survival next to genes
genes_survival_matrix_baseline_Placebo<-as.data.frame(dds_baseline_Placebo$Survival)
colnames(genes_survival_matrix_baseline_Placebo)[1] <- "Survival"
cts_baseline_Placebo<-t(counts(dds_baseline_Placebo, normalized=TRUE))
genes_survival_matrix_baseline_Placebo<-cbind(genes_survival_matrix_baseline_Placebo, cts_baseline_Placebo)

genes_survival_matrix_baseline_Teplizumab<-as.data.frame(dds_baseline_Teplizumab$Survival)
colnames(genes_survival_matrix_baseline_Teplizumab)[1] <- "Survival"
cts_baseline_Teplizumab<-t(counts(dds_baseline_Teplizumab, normalized=TRUE))
genes_survival_matrix_baseline_Teplizumab<-cbind(genes_survival_matrix_baseline_Teplizumab, cts_baseline_Teplizumab)

######Correlation calculations########

#First Teplizumab treated patient
#making correlation matrix
cor.mat_baseline_Teplizumab <- rcorr(as.matrix(genes_survival_matrix_baseline_Teplizumab), type = c("pearson"))


#Accessing the two layers in the object and saving them separately then merging it into one readable table
d_p_corr_baseline_Teplizumab = data.frame(row.names(cor.mat_baseline_Teplizumab$r), Survival = cor.mat_baseline_Teplizumab$r[,"Survival"])
d_p_sign_baseline_Teplizumab = data.frame(row.names(cor.mat_baseline_Teplizumab$P), Survival = cor.mat_baseline_Teplizumab$P[,"Survival"])
merge_corr_sign_baseline_Teplizumab <-  merge(d_p_corr_baseline_Teplizumab, d_p_sign_baseline_Teplizumab, by = "row.names")
row.names(merge_corr_sign_baseline_Teplizumab) = merge_corr_sign_baseline_Teplizumab[,1]
merge_corr_sign_baseline_Teplizumab <- merge_corr_sign_baseline_Teplizumab[,c(1,3,5)]
# assigning new names to the columns of the data frame 
colnames(merge_corr_sign_baseline_Teplizumab) <- c('Symbol','Pearson_corr','p_value')
merge_corr_sign_baseline_Teplizumab = merge_corr_sign_baseline_Teplizumab[order(-merge_corr_sign_baseline_Teplizumab$Pearson_corr),]

#The same for the placebo
#making correlation matrix
cor.mat_baseline_Placebo <- rcorr(as.matrix(genes_survival_matrix_baseline_Placebo), type = c("pearson"))


#Accessing the two layers in the object and saving them separately then merging it into one readable table
d_p_corr_baseline_Placebo = data.frame(row.names(cor.mat_baseline_Placebo$r), Survival = cor.mat_baseline_Placebo$r[,"Survival"])
d_p_sign_baseline_Placebo = data.frame(row.names(cor.mat_baseline_Placebo$P), Survival = cor.mat_baseline_Placebo$P[,"Survival"])
merge_corr_sign_baseline_Placebo <-  merge(d_p_corr_baseline_Placebo, d_p_sign_baseline_Placebo, by = "row.names")
row.names(merge_corr_sign_baseline_Placebo) = merge_corr_sign_baseline_Placebo[,1]
merge_corr_sign_baseline_Placebo <- merge_corr_sign_baseline_Placebo[,c(1,3,5)]
# assigning new names to the columns of the data frame 
colnames(merge_corr_sign_baseline_Placebo) <- c('Symbol','Pearson_corr','p_value')
merge_corr_sign_baseline_Placebo = merge_corr_sign_baseline_Placebo[order(-merge_corr_sign_baseline_Placebo$Pearson_corr),]


#Select only for significant
corr_baseline_Placebo_signif <-  merge_corr_sign_baseline_Placebo[merge_corr_sign_baseline_Placebo$p_value<(0.05),] 
corr_baseline_Teplizumab_signif <-  merge_corr_sign_baseline_Teplizumab[merge_corr_sign_baseline_Teplizumab$p_value<(0.05),] 


Full_Tepli_signature <- corr_baseline_Teplizumab_signif[ -which(corr_baseline_Teplizumab_signif$Symbol %in% corr_baseline_Placebo_signif$Symbol),]
Teplizumab_positive_corr_only <-Full_Tepli_signature[Full_Tepli_signature$Pearson_corr>0,]
Teplizumab_negative_corr_only<-Full_Tepli_signature[Full_Tepli_signature$Pearson_corr<0,]


#Pearson correlation table for Supplemental Table 7
write.csv(merge_corr_sign_baseline_Teplizumab, "TN10_corr_sign_Baseline.csv")

#R signature - Baseline TN10 table for Supplemental Table 7
write.csv(Teplizumab_positive_corr_only, "TN10_Baseline_R_signature.csv")

#NR signature - Baseline TN10 table for Supplemental Table 7
write.csv(Teplizumab_negative_corr_only, "TN10_Baseline_NR_signature.csv")

#####Loading Human scRNAseq data########
library(Seurat)
library(dplyr)
library(ggplot2)
library(future)
library(SCpubr)
library(matrixStats)

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
TN10_Baseline_Signature <- list(TN10_Baseline_R_sign = na.omit(Teplizumab_positive_corr_only$Symbol),
                                 TN10_Baseline_NR_sign = na.omit(Teplizumab_negative_corr_only$Symbol)
)

WB_NO <- AddModuleScore(WB_NO, TN10_Baseline_Signature[1], name = names(TN10_Baseline_Signature[1]))
WB_NO <- AddModuleScore(WB_NO, TN10_Baseline_Signature[2], name = names(TN10_Baseline_Signature[2]), slot = "counts")


WB_NO_sub <- subset(WB_NO, features = TN10_Baseline_Signature[[1]])
expression_matrix_R <- GetAssayData(object = WB_NO_sub, layer = "data")
expression_matrix_R_linear <- expm1(x = expression_matrix_R)
per_cell_average_R_score_linear <- colMeans(expression_matrix_R_linear) 
WB_NO$TN10_Baseline_R_sign <- per_cell_average_R_score_linear[colnames(WB_NO)]

VlnPlot(WB_NO, features = "TN10_Baseline_R_sign")

WB_NO_sub <- subset(WB_NO, features = TN10_Baseline_Signature[[2]])
expression_matrix_NR <- GetAssayData(object = WB_NO_sub, layer = "data")
expression_matrix_NR_linear <- expm1(x = expression_matrix_NR)
per_cell_average_NR_score_linear <- colMeans(expression_matrix_NR_linear) 
WB_NO$TN10_Baseline_NR_sign <- per_cell_average_NR_score_linear[colnames(WB_NO)]

VlnPlot(WB_NO, features = "TN10_Baseline_NR_sign", pt.size = 0.01)



expression_matrix <- GetAssayData(object = subset(WB_NO, features = TN10_Baseline_Signature[[1]]), layer = "data")
mean_expression_per_gene <- rowMeans(expression_matrix)
sd_expression_per_gene <- rowSds(expression_matrix)
standardized_expression_matrix <- (expression_matrix - mean_expression_per_gene) / sd_expression_per_gene
standardized_expression_matrix[is.na(standardized_expression_matrix)] <- 0
standardized_expression_matrix[is.infinite(standardized_expression_matrix)] <- 0
z_score_per_cell <- colMeans(standardized_expression_matrix)
WB_NO$TN10_Baseline_R_sign <- z_score_per_cell


#expression_matrix <- GetAssayData(object = subset(WB_NO, features = TN10_Baseline_Signature[[2]]), layer = "data")
mean_expression_per_gene <- rowMeans(expression_matrix_NR_linear)
sd_expression_per_gene <- rowSds(expression_matrix_NR_linear)
standardized_expression_matrix <- (expression_matrix_NR_linear - mean_expression_per_gene) / sd_expression_per_gene
standardized_expression_matrix[is.na(standardized_expression_matrix)] <- 0
standardized_expression_matrix[is.infinite(standardized_expression_matrix)] <- 0
#z_score_per_cell <- colMeans(expression_matrix)
z_score_per_cell <- colMeans(standardized_expression_matrix)
WB_NO$TN10_Baseline_NR_sign_zscore <- z_score_per_cell
VlnPlot(WB_NO, features = "TN10_Baseline_NR_sign_zscore")

###Test Debug
test_matrix <- subset(WB_NO, features = TN10_Baseline_Signature[[2]])
test_matrix <- subset(WB_NO, features = c("TN10_Baseline_NR_sign"))
NR_avg_expression <- AverageExpression(test_matrix, group.by = "Annotation", layer = "data" )





#Feature plot of the  average expression levels of AbATE Baseline gene signature with SCpubr
p_TN10_Baseline_R_sign <- do_FeaturePlot(sample = WB_NO, features = paste0(names(TN10_Baseline_Signature[1])), enforce_symmetry = F, diverging.palette =  "RdYlBu",
                                         min.cutoff = 0.1,
                                         plot_cell_borders = FALSE, use_viridis = TRUE, viridis.palette ="viridis", pt.size = 2,
                                         raster = TRUE, raster.dpi = 1024)

p_TN10_Baseline_NR_sign <-do_FeaturePlot(sample = WB_NO, features = paste0(names(TN10_Baseline_Signature[2])), enforce_symmetry = F, 
                                          #min.cutoff = -0.05, 
                                          #max.cutoff = 0.08,
                                          plot_cell_borders = FALSE, use_viridis = TRUE, viridis.palette ="viridis", pt.size = 2,
                                          raster = TRUE, raster.dpi = 1024)

ggsave(file="Fig6A.svg", plot=p_TN10_Baseline_R_sign, width=8, height=7)
ggsave(file="Fig6B.svg", plot=p_TN10_Baseline_NR_sign, width=8, height=7)

#We will now calculate the Average expression of the R and NR signature for each cell and plot in a Violin plot
#R signature
sub_WB_NO <- subset(WB_NO, features = TN10_Baseline_Signature[[1]])
expression_matrix_subset_R <- GetAssayData(object = sub_WB_NO, layer = "data")
expression_matrix_subset_R_linear <- expm1(x = expression_matrix_subset_R) 
per_cell_average_R_score_linear <- colMeans(expression_matrix_subset_R_linear) # Average on linear scale
WB_NO$TN10_Baseline_R_sign <- per_cell_average_R_score_linear[colnames(WB_NO)]

#NR signature
sub_WB_NO <- subset(WB_NO, features = TN10_Baseline_Signature[[2]])
expression_matrix_subset_NR <- GetAssayData(object = sub_WB_NO, layer = "data")
expression_matrix_subset_NR_linear <- expm1(x = expression_matrix_subset_NR) 
per_cell_average_NR_score_linear <- colMeans(expression_matrix_subset_NR_linear) # Average on linear scale
WB_NO$TN10_Baseline_NR_sign <- per_cell_average_NR_score_linear[colnames(WB_NO)]

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
feature_to_plot <- "TN10_Baseline_R_sign"
grouping_variable <- "Annotation"

#Get the needed data
plot_data <- FetchData(WB_NO, vars = c(feature_to_plot, grouping_variable))
colnames(plot_data) <- c("Expression", "Group")
names(cols)<- levels(plot_data$Group)

#Manually reorganized to match fig6E order
desired_order <- c("T cells", "DCs","Plasma cells","Neutrophils 0","Monocytes","B cells","NK/CD8 T cells","Neutrophils 2","Neutrophils 4", "Neutrophils 3", "Megakaryocytes", "Basophils","Neutrophils 1" )
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
        axis.title.y = element_text(size = 25, margin = margin(r = 20)),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 15),
        legend.position = "none")


ggsave(file="Fig6C.png", plot=vln_plot_expr_R_sign, width=6, height=6, dpi = 300)

#Violin plot of the NR signature
feature_to_plot <- "TN10_Baseline_NR_sign"
grouping_variable <- "Annotation"

#Get the needed data
plot_data <- FetchData(WB_NO, vars = c(feature_to_plot, grouping_variable))
colnames(plot_data) <- c("Expression", "Group")
names(cols)<- levels(plot_data$Group)

#Manually reorganized to match fig6F order
desired_order <- c("Neutrophils 1", "B cells", "Basophils", "Megakaryocytes", "NK/CD8 T cells", "Plasma cells", "T cells", "DCs", "Neutrophils 3", "Neutrophils 0", "Monocytes", "Neutrophils 4", "Neutrophils 2" )
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
        axis.title.y = element_text(size = 25, margin = margin(r = 20)),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 15),
        legend.position = "none")

ggsave(file="Fig6D.png", plot=vln_plot_expr_NR_sign, width=6, height=6, dpi = 300)

#Let's plot pseudobulk average expression per cell type

#Plot average expression value per celltype for Positive correlation
R_avg_expression <- AverageExpression(WB_NO, features = TN10_Baseline_Signature[[1]],group.by = "Annotation" )
R_avg_expression <- R_avg_expression$RNA

mean_expression_per_celltype <- colMeans(R_avg_expression)
mean_expression_per_celltype <- as.data.frame(mean_expression_per_celltype)
mean_expression_per_celltype$cell_type <- rownames(mean_expression_per_celltype)
colnames(mean_expression_per_celltype)[colnames(mean_expression_per_celltype) == "mean_expression_per_celltype"] <- "Mean_Pearson_R_sign"


ranked_cell_types_R <- mean_expression_per_celltype %>%
  arrange(desc(Mean_Pearson_R_sign))
ranked_cell_types_R$cell_type <- factor(ranked_cell_types_R$cell_type, 
                                      levels = ranked_cell_types_R$cell_type)


plot_expr_R_sign <- ggplot(ranked_cell_types_R, aes(x = cell_type, y = Mean_Pearson_R_sign)) +
  geom_point(size = 4) +
  labs(x = "Cell Type", y = "Average expression") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 20, color = "black"),
        plot.title = element_text(size = 25),
        axis.title.y = element_text(size = 25, margin = margin(r = 20)),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 15)) +
  labs(title = "R signature - TN10 Baseline")

plot_expr_R_sign

ggsave(file="Fig5E.svg", plot=plot_expr_R_sign, width=6, height=6)

#Plot average expression value per celltype for Positive correlation
NR_avg_expression <- AverageExpression(WB_NO, features = TN10_Baseline_Signature[[2]],group.by = "Annotation", layer = "data" )
NR_avg_expression <- NR_avg_expression$RNA

mean_expression_per_celltype <- colMeans(NR_avg_expression)
mean_expression_per_celltype <- as.data.frame(mean_expression_per_celltype)
mean_expression_per_celltype$cell_type <- rownames(mean_expression_per_celltype)
colnames(mean_expression_per_celltype)[colnames(mean_expression_per_celltype) == "mean_expression_per_celltype"] <- "Mean_Pearson_NR_sign"


ranked_cell_types_NR <- mean_expression_per_celltype %>%
  arrange(desc(Mean_Pearson_NR_sign))
ranked_cell_types_NR$cell_type <- factor(ranked_cell_types_NR$cell_type, 
                                      levels = ranked_cell_types_NR$cell_type)


plot_expr_NR_sign <- ggplot(ranked_cell_types_NR, aes(x = cell_type, y = Mean_Pearson_NR_sign)) +
  geom_point(size = 4) +
  labs(x = "Cell Type", y = "Average expression") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 20, color = "black"),
        plot.title = element_text(size = 25),
        axis.title.y = element_text(size = 25, margin = margin(r = 20)),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 15)) +
  labs(title = "NR signature - TN10 Baseline")

plot_expr_NR_sign

ggsave(file="Fig5F.svg", plot=plot_expr_NR_sign, width=6, height=6)

