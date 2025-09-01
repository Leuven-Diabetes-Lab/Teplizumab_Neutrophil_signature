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

gene_expr_all <- read.csv("global.normalized_counts.csv")

#making X column character
gene_expr_all$X <- as.character(gene_expr_all$X)
#renaming first column
names(gene_expr_all)[names(gene_expr_all) == 'X'] <- 'Symbol'
# making first column the index
rownames(gene_expr_all) <- gene_expr_all$Symbol
#remove original ID column from data frame
gene_expr_all$Symbol <- NULL

#Separate auc file
auc_6 <- subset(auc_file, auc_file$visitmonth == "6", select = c("Run","status","auc_percent_of_baseline"))
auc_6 = data.frame(auc_6)
#make first column index
rownames(auc_6) <- auc_6$Run

#Make subsets for gene expression
gene_expr_6 <- subset(gene_expr_all, select = c(row.names(auc_6)))

#FILTER GENE EXPRESSION: they have to be expressed at least 10 times in minimum 60% of the patients
#transposing it too so that observations are rows and variables are columns

smallest_count <- 10
keep <- rowSums(gene_expr_6 >= smallest_count) >= ceiling(0.6 * ncol(gene_expr_6))
gene_expr_subset_6 <- gene_expr_6[keep,]
tgene_expr_6 <- data.frame(t(gene_expr_subset_6))

#merge with auc
tgene_expr_all.auc_6 = merge(auc_6, tgene_expr_6, by = "row.names")
row.names(tgene_expr_all.auc_6) = tgene_expr_all.auc_6[,1]
tgene_expr_all.auc_6 = tgene_expr_all.auc_6[,-(1:3)]

#correlation with p-values
#making correlation matrix
cor.mat_6_with_p <- rcorr(as.matrix(tgene_expr_all.auc_6), type = c("pearson"))

#Accessing the two layers in the object and saving them separately then merging it into one readable table
d_p_corr_6 = data.frame(row.names(cor.mat_6_with_p$r), AUC = cor.mat_6_with_p$r[,"auc_percent_of_baseline"])
d_p_sign_6 = data.frame(row.names(cor.mat_6_with_p$P), AUC = cor.mat_6_with_p$P[,"auc_percent_of_baseline"])
merge_corr_sign_6 <-  merge(d_p_corr_6, d_p_sign_6, by = "row.names")
row.names(merge_corr_sign_6) = merge_corr_sign_6[,1]
merge_corr_sign_6 <- merge_corr_sign_6[,c(1,3,5)]
# assigning new names to the columns of the data frame 
colnames(merge_corr_sign_6) <- c('Symbol','Pearson_corr','p_value')
merge_corr_sign_6 = merge_corr_sign_6[order(-merge_corr_sign_6$Pearson_corr),]
merge_corr_sign_6<-merge_corr_sign_6[!row.names(merge_corr_sign_6) %in% "auc_percent_of_baseline",]

#Pearson correlation table for Supplemental Table 4
write.csv(merge_corr_sign_6, "AbATE_corr_sign_6.csv")

#selecting only the significant genes to make the AbATE Month 6 gene signature
Abate_signature<-merge_corr_sign_6[merge_corr_sign_6$p_value<(0.05),] 

#R signature - Month 6 AbATE table for Supplemental Table 4
Abate_signature_R<-Abate_signature[Abate_signature$Pearson_corr >0,] 
write.csv(Abate_signature_R, "AbATE_Month6_R_signature.csv")

#NR signature - Month 6 AbATE table for Supplemental Table 4
Abate_signature_NR<-Abate_signature[Abate_signature$Pearson_corr <0,] 
write.csv(Abate_signature_NR, "AbATE_Month6_NR_signature.csv")

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

#Cluster plot for Figure 4B
DimPlot(WB_NO, reduction = "umap")

#Calculate the average expression levels of AbATE Month 6 gene signature for each cell
AbATE_Month6_Signature <- list(AbATE_Month6_R_sign = na.omit(Abate_signature_R$Symbol),
                               AbATE_Month6_NR_sign = na.omit(Abate_signature_NR$Symbol)
)

WB_NO <- AddModuleScore(WB_NO, AbATE_Month6_Signature[1], name = names(AbATE_Month6_Signature[1]))
WB_NO <- AddModuleScore(WB_NO, AbATE_Month6_Signature[2], name = names(AbATE_Month6_Signature[2]))


#Feature plot of the  average expression levels of AbATE Month 6 gene signature with SCpubr
p_AbATE_Month6_R_sign <- do_FeaturePlot(sample = WB_NO, features = paste0(names(AbATE_Month6_Signature[1]),"1"), enforce_symmetry = F, diverging.palette =  "RdYlBu",
                              min.cutoff = -0.1,
                              plot_cell_borders = FALSE, use_viridis = TRUE, viridis.palette ="viridis", pt.size = 2,
                              raster = TRUE, raster.dpi = 1024)

p_AbATE_Month6_NR_sign <-do_FeaturePlot(sample = WB_NO, features = paste0(names(AbATE_Month6_Signature[2]),"1"), enforce_symmetry = F, 
                              min.cutoff = -0.05, 
                              max.cutoff = 0.05,
                              plot_cell_borders = FALSE, use_viridis = TRUE, viridis.palette ="viridis", pt.size = 2,
                              raster = TRUE, raster.dpi = 1024)

ggsave(file="Fig4C.svg", plot=p_AbATE_Month6_R_sign, width=8, height=7)
ggsave(file="Fig4D.svg", plot=p_AbATE_Month6_NR_sign, width=8, height=7)

#We will now calculate the Average expression of the R and NR signature for each cell and plot in a Violin plot
#R signature
sub_WB_NO <- subset(WB_NO, features = AbATE_Month6_Signature[[1]])
expression_matrix_subset_R <- GetAssayData(object = sub_WB_NO, layer = "data")
expression_matrix_subset_R_linear <- expm1(x = expression_matrix_subset_R) 
per_cell_average_R_score_linear <- colMeans(expression_matrix_subset_R_linear) # Average on linear scale
WB_NO$AbATE_Month6_R_sign <- per_cell_average_R_score_linear[colnames(WB_NO)]

#NR signature
sub_WB_NO <- subset(WB_NO, features = AbATE_Month6_Signature[[2]])
expression_matrix_subset_NR <- GetAssayData(object = sub_WB_NO, layer = "data")
expression_matrix_subset_NR_linear <- expm1(x = expression_matrix_subset_NR) 
per_cell_average_NR_score_linear <- colMeans(expression_matrix_subset_NR_linear) # Average on linear scale
WB_NO$AbATE_Month6_NR_sign <- per_cell_average_NR_score_linear[colnames(WB_NO)]

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
feature_to_plot <- "AbATE_Month6_R_sign"
grouping_variable <- "Annotation"

#Get the needed data
plot_data <- FetchData(WB_NO, vars = c(feature_to_plot, grouping_variable))
colnames(plot_data) <- c("Expression", "Group")
names(cols)<- levels(plot_data$Group)

#Manually reorganized to match fig5E order
# desired_order <- c("NK/CD8 T cells","Plasma cells","DCs","T cells","Monocytes","B cells","Neutrophils 0","Basophils","Neutrophils 1", "Megakaryocytes","Neutrophils 2","Neutrophils 3","Neutrophils 4"  )
# plot_data$Group <- fct_relevel(plot_data$Group, desired_order)

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
    #breaks = desired_order
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


ggsave(file="Fig4E.png", plot=vln_plot_expr_R_sign, width=6, height=6, dpi = 300)

#Violin plot of the NR signature
feature_to_plot <- "AbATE_Month6_NR_sign"
grouping_variable <- "Annotation"

#Get the needed data
plot_data <- FetchData(WB_NO, vars = c(feature_to_plot, grouping_variable))
colnames(plot_data) <- c("Expression", "Group")
names(cols)<- levels(plot_data$Group)

#Manually reorganized to match fig6F order
# desired_order <- c("Neutrophils 3","Neutrophils 1","Neutrophils 4","Neutrophils 2", "Megakaryocytes", "Basophils", "Neutrophils 0", "Monocytes", "NK/CD8 T cells","B cells","T cells","DCs","Plasma cells")
# plot_data$Group <- fct_relevel(plot_data$Group, desired_order)

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
    #breaks = desired_order
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

ggsave(file="Fig4F.png", plot=vln_plot_expr_NR_sign, width=6, height=6, dpi = 300)

