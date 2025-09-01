#Loading libraries:
library("DESeq2")
library("dplyr")
library("readxl")
library("Hmisc")
library("ggpubr")


#Load the dataset and previously generated Baseline signatures
dds_merge <- readRDS("TN10_AbATE_merged_Full_Annotation.rds") 


Abate_Baseline_full <- read.csv("AbATE_corr_sign_Baseline.csv")
TN10_Baseline_sign_R <- read.csv("TN10_Baseline_R_signature.csv")
TN10_Baseline_sign_NR <- read.csv("TN10_Baseline_NR_signature.csv")
Full_Tepli_signature<-rbind(TN10_Baseline_sign_R,TN10_Baseline_sign_NR)

#Select genes from the AbATE baseline Pearson correlation with p<0.05 and generate positive and negative lists
Abate_corr_signif<-Abate_Baseline_full[Abate_Baseline_full$p_value<(0.05),] 

#Compare both list
TN10_AbATE_both <-merge(x = Full_Tepli_signature, y = Abate_corr_signif, by.x='Symbol', by.y='Symbol')
Combined_Baseline_sign_R<-TN10_AbATE_both[which(TN10_AbATE_both$Pearson_corr.x >0 & TN10_AbATE_both$Pearson_corr.y >0),]
Combined_Baseline_sign_NR<-TN10_AbATE_both[which(TN10_AbATE_both$Pearson_corr.x <0 & TN10_AbATE_both$Pearson_corr.y <0),]



#Pearson correlation table for Supplemental Table 7
#write.csv(merge_corr_sign_baseline_Teplizumab, "TN10_corr_sign_Baseline.csv")

#R signature - Baseline TN10 table for Supplemental Table 7
#write.csv(Teplizumab_positive_corr_only, "TN10_Baseline_R_signature.csv")

#NR signature - Baseline TN10 table for Supplemental Table 7
#write.csv(Teplizumab_negative_corr_only, "TN10_Baseline_NR_signature.csv")

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
Combined_Baseline_Signature <- list(Combined_Baseline_R_sign = na.omit(Combined_Baseline_sign_R$Symbol),
                                    Combined_Baseline_NR_sign = na.omit(Combined_Baseline_sign_NR$Symbol)
)

WB_NO <- AddModuleScore(WB_NO, Combined_Baseline_Signature[1], name = names(Combined_Baseline_Signature[1]))
WB_NO <- AddModuleScore(WB_NO, Combined_Baseline_Signature[2], name = names(Combined_Baseline_Signature[2]))


#Feature plot of the  average expression levels of AbATE Baseline gene signature with SCpubr
p_Combined_Baseline_R_sign <- do_FeaturePlot(sample = WB_NO, features = paste0(names(Combined_Baseline_Signature[1]),"1"), enforce_symmetry = F, diverging.palette =  "RdYlBu",
                                         min.cutoff = -0.2,
                                         plot_cell_borders = FALSE, use_viridis = TRUE, viridis.palette ="viridis", pt.size = 2,
                                         raster = TRUE, raster.dpi = 1024)

p_Combined_Baseline_NR_sign <-do_FeaturePlot(sample = WB_NO, features = paste0(names(Combined_Baseline_Signature[2]),"1"), enforce_symmetry = F, 
                                             min.cutoff = -0.13, 
                                             max.cutoff = 0.3,
                                         plot_cell_borders = FALSE, use_viridis = TRUE, viridis.palette ="viridis", pt.size = 2,
                                         raster = TRUE, raster.dpi = 1024)

ggsave(file="Fig7A.svg", plot=p_TN10_Baseline_R_sign, width=8, height=7)
ggsave(file="Fig6B.svg", plot=p_TN10_Baseline_NR_sign, width=8, height=7)

#We will now calculate the Average expression of the R and NR signature for each cell and plot in a Violin plot
#R signature
sub_WB_NO <- subset(WB_NO, features = Combined_Baseline_Signature[[1]])
expression_matrix_subset_R <- GetAssayData(object = sub_WB_NO, layer = "data")
expression_matrix_subset_R_linear <- expm1(x = expression_matrix_subset_R) 
per_cell_average_R_score_linear <- colMeans(expression_matrix_subset_R_linear) # Average on linear scale
WB_NO$Combined_Baseline_R_sign <- per_cell_average_R_score_linear[colnames(WB_NO)]

#NR signature
sub_WB_NO <- subset(WB_NO, features = Combined_Baseline_Signature[[2]])
expression_matrix_subset_NR <- GetAssayData(object = sub_WB_NO, layer = "data")
expression_matrix_subset_NR_linear <- expm1(x = expression_matrix_subset_NR) 
per_cell_average_NR_score_linear <- colMeans(expression_matrix_subset_NR_linear) # Average on linear scale
WB_NO$Combined_Baseline_NR_sign <- per_cell_average_NR_score_linear[colnames(WB_NO)]

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
feature_to_plot <- "Combined_Baseline_R_sign"
grouping_variable <- "Annotation"

#Get the needed data
plot_data <- FetchData(WB_NO, vars = c(feature_to_plot, grouping_variable))
colnames(plot_data) <- c("Expression", "Group")
names(cols)<- levels(plot_data$Group)

#Manually reorganized to match fig6E order
desired_order <- c("T cells","B cells","Plasma cells","DCs","NK/CD8 T cells","Monocytes","Neutrophils 0","Basophils","Megakaryocytes","Neutrophils 1","Neutrophils 2","Neutrophils 3", "Neutrophils 4")
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


ggsave(file="Fig7C.png", plot=vln_plot_expr_R_sign, width=6, height=6, dpi = 300)

#Violin plot of the NR signature
feature_to_plot <- "Combined_Baseline_NR_sign"
grouping_variable <- "Annotation"

#Get the needed data
plot_data <- FetchData(WB_NO, vars = c(feature_to_plot, grouping_variable))
colnames(plot_data) <- c("Expression", "Group")
names(cols)<- levels(plot_data$Group)

#Manually reorganized to match fig6F order
desired_order <- c("Neutrophils 1","Basophils","Neutrophils 4","Megakaryocytes","Neutrophils 3", "Neutrophils 2", "NK/CD8 T cells", "T cells","B cells","Neutrophils 0","DCs","Monocytes","Plasma cells")
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

ggsave(file="Fig7D.png", plot=vln_plot_expr_NR_sign, width=6, height=6, dpi = 300)

#Let's plot pseudobulk average expression per cell type

#Plot average expression value per celltype for Positive correlation
R_avg_expression <- AverageExpression(WB_NO, features = Combined_Baseline_Signature[[1]],group.by = "Annotation" )
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
  labs(title = "R signature - Combined Baseline")

plot_expr_R_sign

ggsave(file="Fig7E.svg", plot=plot_expr_R_sign, width=6, height=6)

#Plot average expression value per celltype for Positive correlation
NR_avg_expression <- AverageExpression(WB_NO, features = Combined_Baseline_Signature[[2]],group.by = "Annotation", layer = "data" )
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
  labs(title = "NR signature - Combined Baseline")

plot_expr_NR_sign

ggsave(file="Fig7F.svg", plot=plot_expr_NR_sign, width=6, height=6)
