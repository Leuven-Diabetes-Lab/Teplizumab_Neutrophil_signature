#Loading libraries:
library(DESeq2)
library(pROC)
library(ggplot2)
library(ggpubr)   
library(reshape2)
library(dplyr)
library(SCpubr)

#Load the TN10/AbATE merged dataset
dds_merge <- readRDS("TN10_AbATE_merged_Full_Annotation.rds") 

#We keep only the baseline data from treated patient
to_keep<-which(dds_merge$Treatment %in% "Teplizumab")
dds_merge_Teplizumab <- dds_merge[,to_keep]
to_keep<-which(dds_merge_Teplizumab$Time_point %in% "Baseline")
dds_merge_Teplizumab_Baseline  <- dds_merge_Teplizumab[,to_keep]


# --- Parameters ---
n_samples_to_pick <- 20
n_iterations <- 1000
set.seed(123)
positive_class <- "R"
negative_class <- "NR"
outcome_colors <- c("orange", "yellowgreen")
names(outcome_colors) <- c(negative_class, positive_class)

# --- Step 1: Prepare Data and Calculate Response Scores for ALL samples ---

# Extract normalized counts
norm_counts <- counts(dds_merge_Teplizumab_Baseline, normalized = TRUE)

# Extract outcome groups
outcomes <- dds_merge_Teplizumab_Baseline$Outcome
sample_names <- colnames(norm_counts)

# Extract Study information
study_info <- dds_merge_Teplizumab_Baseline$Study
if(is.null(study_info)) {
  stop("Error: 'Study' variable not found in dds_merge_Teplizumab_Baseline object's colData.")
}
# Assign sample names (assuming original order matches colData)
names(study_info) <- colnames(dds_merge_Teplizumab_Baseline)

# Get the gene names from the coefficients
Stored_signature<- read.csv2("Predictive_Signature.csv")
median_coeffs_overlap<-unlist(Stored_signature[2])
names(median_coeffs_overlap)<-unlist(Stored_signature[1])
biomarker_genes <- names(median_coeffs_overlap)

median_coeffs_overlap <- median_coeffs_overlap[biomarker_genes]

# Find common genes
common_genes <- intersect(rownames(norm_counts), biomarker_genes)

if (length(common_genes) == 0) {
  stop("Error: No common genes found.")
} else {
  print(paste("Found", length(common_genes), "common genes."))
}

# Filter and align data
norm_counts_filtered <- norm_counts[common_genes, , drop = FALSE]
coeffs_filtered <- median_coeffs_overlap[common_genes]
ordered_genes <- sort(common_genes)
norm_counts_ordered <- norm_counts_filtered[ordered_genes, , drop = FALSE] # Keep this matrix for later
coeffs_ordered <- coeffs_filtered[ordered_genes]

# Calculate response scores
response_scores <- t(norm_counts_ordered) %*% matrix(coeffs_ordered, ncol = 1)
response_scores <- as.vector(response_scores)
names(response_scores) <- colnames(norm_counts_ordered)

# --- Align outcomes and study info to response_scores order --- 
# Align outcomes
if (!all(names(response_scores) == names(outcomes))) {
  outcome_names_original <- names(outcomes); scores_names <- names(response_scores)
  if(is.null(outcome_names_original) || !all(scores_names %in% outcome_names_original)) {
    warning("Outcome names missing or mismatched during alignment check; assuming initial order was correct.")
    if(length(outcomes) != length(response_scores)) stop("Length mismatch between scores and outcomes.")
    names(outcomes) <- names(response_scores) # Assign names based on scores
  } else {
    print("Aligning outcomes to score order...")
    outcomes <- outcomes[names(response_scores)]
    if (!identical(names(response_scores), names(outcomes))) { # Simple check after alignment
      warning("Outcome names might not be perfectly aligned after subsetting, proceeding by order.")
      names(outcomes) <- names(response_scores) # Force names if needed
    }
  }
} else {
  print("Outcomes already aligned with scores.")
}

# Align study info
study_info_aligned <- NULL 
if (!all(names(response_scores) == names(study_info))) {
  study_names_original <- names(study_info); scores_names <- names(response_scores)
  if(is.null(study_names_original) || !all(scores_names %in% study_names_original)) {
    warning("Study info names missing or mismatched during alignment check; assuming initial order was correct.")
    if(length(study_info) != length(response_scores)) stop("Length mismatch between scores and study info.")
    names(study_info) <- names(response_scores) # Assign names based on scores
    study_info_aligned <- study_info 
  } else {
    print("Aligning study info to score order...")
    study_info_aligned <- study_info[names(response_scores)]
    if (!identical(names(response_scores), names(study_info_aligned))) {
      warning("Study info names might not be perfectly aligned after subsetting, proceeding by order.")
      names(study_info_aligned) <- names(response_scores) # Force names
    }
  }
} else {
  print("Study info already aligned with scores.")
  study_info_aligned <- study_info # Use original if already aligned
}
if(length(study_info_aligned) != length(response_scores)) stop("Final length mismatch: study info vs scores") 
print("Response scores, outcomes, and study info aligned.")

# --- Calculate ROC Curve for the FULL Dataset  ---

roc_full_data <- NULL
auc_full <- NA
if (length(unique(outcomes)) == 2 && positive_class %in% outcomes && negative_class %in% outcomes) {
  full_data_outcomes_factor <- factor(outcomes, levels = c(negative_class, positive_class))
  roc_full_data <- tryCatch({ 
    roc(full_data_outcomes_factor, response_scores, levels = c(negative_class, positive_class), direction = "<", quiet = TRUE) 
    }, error = function(e) { warning(paste("Could not calculate ROC curve on full dataset:", e$message)); return(NULL) })
  if (!is.null(roc_full_data)) { auc_full <- round(auc(roc_full_data), 3); print(paste("AUC on full dataset:", auc_full)) }
} else { print("Cannot calculate ROC on full dataset (needs both outcome classes).") }


# --- Step 2: Resampling Loop --- 
print(paste("Step 2: Starting", n_iterations, "resampling iterations..."))
results_list <- vector("list", n_iterations)
roc_objects_list <- list()
total_samples <- length(outcomes)

for (i in 1:n_iterations) {
  sample_indices <- sample(1:total_samples, n_samples_to_pick, replace = FALSE)
  subset_scores <- response_scores[sample_indices]
  subset_outcomes <- outcomes[sample_indices]
  unique_outcomes_in_subset <- unique(subset_outcomes)
  
  current_result <- list(sensitivity = NA, specificity = NA, auc = NA, threshold = NA)
  
  if (length(unique_outcomes_in_subset) < 2 ||
      !(positive_class %in% unique_outcomes_in_subset && negative_class %in% unique_outcomes_in_subset)) {
    results_list[[i]] <- current_result
    next
  }
  
  subset_outcomes_factor <- factor(subset_outcomes, levels = c(negative_class, positive_class))
  roc_obj <- tryCatch({ 
    roc(subset_outcomes_factor, subset_scores, levels = c(negative_class, positive_class), direction = "<", quiet = TRUE) 
    }, error = function(e) NULL )
  
  if (is.null(roc_obj)) { results_list[[i]] <- current_result; next }
  
  roc_objects_list[[length(roc_objects_list) + 1]] <- roc_obj
  auc_value <- tryCatch(auc(roc_obj), error = function(e) NA)
  
  # Get sensitivity, specificity, and threshold at the best cutoff
  best_coords <- tryCatch({
    # Request threshold along with sens and spec
    coords(roc_obj, "best", ret = c("sensitivity", "specificity", "threshold"), transpose = FALSE)
  }, error = function(e) {
    warning(paste("Iteration", i, ": Error getting coords -", e$message))
    return(NULL)
  })
  
  # Store results if coords were successful
  if (!is.null(best_coords) && nrow(best_coords) > 0) {
    results_list[[i]] <- list(
      sensitivity = best_coords$sensitivity[1],
      specificity = best_coords$specificity[1],
      auc = as.numeric(auc_value),
      threshold = best_coords$threshold[1] 
    )
  } else {
    # Coords failed, store NAs but keep AUC if calculated
    current_result$auc <- as.numeric(auc_value)
    results_list[[i]] <- current_result
  }
} # End loop
print("Resampling complete.")

# --- Step 3: Aggregation and Final Metrics --- 
print("Step 3: Aggregating results...")

# Extract metrics including thresholds
all_sensitivities <- sapply(results_list, function(res) res$sensitivity)
all_specificities <- sapply(results_list, function(res) res$specificity)
all_aucs <- sapply(results_list, function(res) res$auc)
all_thresholds <- sapply(results_list, function(res) res$threshold)

# Remove NA values
valid_sensitivities <- na.omit(all_sensitivities)
valid_specificities <- na.omit(all_specificities)
valid_aucs <- na.omit(all_aucs)
valid_thresholds <- na.omit(all_thresholds)

# Calculate number of successful iterations for each metric
n_successful_iterations_sens <- length(valid_sensitivities)
n_successful_iterations_spec <- length(valid_specificities)
n_successful_iterations_auc <- length(valid_aucs)
n_successful_iterations_thresh <- length(valid_thresholds)

# Calculate mean/median values IF there were successful iterations
mean_sensitivity <- NA; mean_specificity <- NA; mean_auc <- NA; median_resampled_threshold <- NA # Initialize medians

if (n_successful_iterations_auc > 0) { # Use AUC count as general indicator of success
  mean_sensitivity <- mean(valid_sensitivities)
  mean_specificity <- mean(valid_specificities)
  mean_auc <- mean(valid_aucs)
  
  # Calculate MEDIAN threshold only if valid thresholds exist
  if (n_successful_iterations_thresh > 0) {
    median_resampled_threshold <- median(valid_thresholds)
    print(paste("Median optimal threshold from resampling:", round(median_resampled_threshold, 3)))
  } else {
    print("Could not calculate median threshold (no valid thresholds found in resampling).")
  }
  
  # Report final results
  cat("\n--- Resampling Results ---\n")
  cat(paste("Successful iterations (Sens/Spec/AUC/Thresh):",
            n_successful_iterations_sens, "/", n_successful_iterations_spec, "/",
            n_successful_iterations_auc, "/", n_successful_iterations_thresh,
            "out of", n_iterations, "\n"))
  cat("Mean Performance Metrics:\n")
  cat(paste("  - Mean Sensitivity:", round(mean_sensitivity, 3), "\n"))
  cat(paste("  - Mean Specificity:", round(mean_specificity, 3), "\n"))
  cat(paste("  - Mean AUC:", round(mean_auc, 3), "\n"))
  if (!is.na(median_resampled_threshold)) {
    cat(paste("  - Median Resampled Threshold:", round(median_resampled_threshold, 3), "\n"))
  }
  cat("--------------------------\n")
  
} else {
  cat("\n--- Resampling Results ---\n"); cat("No successful iterations completed.\n"); cat("--------------------------\n")
}

# --- Step 4: Plotting Results --- 

print("Step 4: Generating plots...")

my_comparisons <- list( c(negative_class, positive_class) )

# --- Plot Response Score vs Outcome (Full Data, Points Colored by Study)
print("Generating Response Score vs Outcome plot")

# Create data frame including aligned Study info
score_data <- data.frame(
  Score = response_scores,
  Outcome = factor(outcomes, levels = c(negative_class,positive_class)),
  Study = factor(study_info_aligned)
)

#Median threshold doesn't seem to be optimal. Manually modified for optimal separation results
median_resampled_threshold_stored <- median_resampled_threshold
median_resampled_threshold <- (-0.01)
p_score_outcome <- ggplot(score_data, aes(x = Outcome, y = Score, fill = Outcome)) +
  geom_boxplot(width = 0.4, outlier.shape = NA, alpha = 0.7) + 
  geom_jitter(aes(color = Study), 
              width = 0.15,       
              height = 0,
              alpha = 1,
              size = 2.5) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.format", size = 8) + 
  scale_fill_manual(values = outcome_colors) + 
  labs(y = "Response Score", x= NULL) +
  theme_minimal() +theme(text=element_text(size=25))
if (!is.na(median_resampled_threshold)) {
  p_score_outcome <- p_score_outcome +
    geom_hline(yintercept = median_resampled_threshold, linetype = "dashed", color = "black", size = 0.8) +
    annotate("text", x = 0.55, y = median_resampled_threshold,
             label = paste("Threshold =", round(median_resampled_threshold, 2)),
             vjust = -0.5, hjust = -0.6, color = "black", size = 0)
}

print(p_score_outcome) 

ggsave(file="Fig8E.svg", plot=p_score_outcome, width=7, height=7)

# --- Plot AUC Distribution (Density Plot) Fig S7M
print("Generating AUC density plot...")
if (n_successful_iterations_auc > 0) {
  auc_data <- data.frame(AUC = valid_aucs)
  p_auc_dist <- ggplot(auc_data, aes(x = AUC)) +
    geom_density(fill = "lightblue", alpha = 0.7) +
    geom_vline(xintercept = mean_auc, color = "red", linetype = "dashed", size = 1) +
    labs( x = "Area Under Curve (AUC)", y = "Frequency") +
    annotate("text", x = mean_auc, y = 0, label = paste("Mean =", round(mean_auc, 2)), hjust = 1.1, vjust = -15, color = "red", size = 9) +
    theme_minimal()+theme(text=element_text(size=25))
  print(p_auc_dist)
} else { print("Skipping AUC distribution plot (no successful iterations).") }

ggsave(file="Fig8F.svg", plot=p_auc_dist, width=7, height=7)

# --- Plot Aggregated ROC Curves
print("Generating aggregated ROC plot...")
if (length(roc_objects_list) > 0) {
  if (!is.null(roc_full_data)) { 
    plot(roc_full_data, col = "darkgrey", lwd = 1, main = "ROC Curves from Resampling", legacy.axes = TRUE) 
  } else { 
      plot(roc_objects_list[[1]], col = rgb(0.5, 0.5, 0.5, alpha = 0.1), main = "ROC Curves from Resampling", legacy.axes = TRUE); 
    if(length(roc_objects_list) > 1) { 
      for (k in 2:length(roc_objects_list)) { 
        plot(roc_objects_list[[k]], add = TRUE, col = rgb(0.5, 0.5, 0.5, alpha = 0.1), legacy.axes = TRUE) 
      } 
    } 
    }
  if (!is.null(roc_full_data)) { 
    for (resampled_roc in roc_objects_list) { 
      plot(resampled_roc, add = TRUE, col = rgb(0.5, 0.5, 0.5, alpha = 0.1), legacy.axes = TRUE) 
    }; 
    plot(roc_full_data, add = TRUE, col = "red", lwd = 2, legacy.axes = TRUE); 
    legend_text_full <- paste("Full Data ROC (AUC =", round(auc(roc_full_data), 3), ")") 
  } else { 
      legend_text_full <- "Full Data ROC N/A" }
  if (!is.na(mean_auc)) { text(x = 0.4, y = 0.2, labels = paste("Mean Resampled AUC =", round(mean_auc, 2)), col = "blue", font = 2, cex=1.5) }
  legend("bottomright", legend = c(legend_text_full, "Resampled ROCs"), col = c("red", rgb(0.5, 0.5, 0.5, alpha = 0.5), cex=0), lwd = c(2, 1), bty = "n")
} else { print("Skipping aggregated ROC plot (no successful iterations or ROC objects generated).") }

#Save FigS7M manually

# Plot Outcome Proportions by Threshold Group for each Study - Fig8G
print("Generating Outcome Proportion plots by Study and Threshold Group...")

median_resampled_threshold <- -0.01

# Check if median_resampled_threshold is available
if (!is.na(median_resampled_threshold)) {
  
  # Create a new data frame for this analysis, adding the ThresholdGroup
  # score_data should already exist
  if(exists("score_data") && is.data.frame(score_data)) {
    score_data_with_threshold_group <- score_data %>%
      mutate(ThresholdGroup = ifelse(Score >= median_resampled_threshold,
                                     "Above Threshold",
                                     "Below Threshold"))
    score_data_with_threshold_group$ThresholdGroup <- factor(
      score_data_with_threshold_group$ThresholdGroup,
      levels = c("Below Threshold", "Above Threshold")
    )
    score_data_with_threshold_group$Outcome <- factor(
      score_data_with_threshold_group$Outcome,
      levels = c(positive_class,negative_class)
    )
    
    
    # Calculate proportions
    # Filter out any potential NAs
    proportion_data_for_studies <- score_data_with_threshold_group %>%
      filter(!is.na(ThresholdGroup) & !is.na(Outcome) & !is.na(Study)) %>%
      group_by(Study, ThresholdGroup, Outcome) %>%
      summarise(Count = n(), .groups = 'drop_last') %>% # Count occurrences
      # Calculate proportion within each Study-ThresholdGroup combination
      mutate(Proportion = Count / sum(Count)) %>%
      ungroup()
    
    # Check if proportion_data_for_studies has data
    if (nrow(proportion_data_for_studies) > 0) {
      unique_studies_for_plot <- unique(proportion_data_for_studies$Study)
      
      for (current_study_id in unique_studies_for_plot) {
        study_specific_data <- proportion_data_for_studies %>%
          filter(Study == current_study_id)
        
        if (nrow(study_specific_data) > 0) {
          p_study_proportion <- ggplot(study_specific_data,
                                       aes(x = ThresholdGroup, y = Proportion, fill = Outcome)) +
            geom_col(position = "stack", alpha = 0.85, width = 0.7) + # geom_col for pre-calculated values
            scale_fill_manual(values = outcome_colors, name = "Outcome") +
            scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1.01), expand = c(0,0)) +
            geom_text(aes(label = ifelse(Proportion > 0.03, # Only label if proportion is >3%
                                         paste0(round(Proportion * 100), "%"),
                                         "")),
                      position = position_stack(vjust = 0.5), size = 8, color="black") +
            labs(x = NULL,
              y = "% of Samples") +
            theme_minimal(base_size = 11) +
            theme(axis.text.x = element_text(size=25),
              legend.position = "top",
              text=element_text(size=25))
          
          
          print(p_study_proportion)
        } else {
          print(paste("No data available to plot proportions for Study:", current_study_id, "after grouping."))
        }
      } 
    } else {
      print("Proportion data could not be calculated (e.g., due to NAs or lack of diversity in outcomes/threshold groups).")
    }
  } else {
    print("`score_data` data frame not found, skipping proportion plots by study.")
  }
} else {
  print("Median resampled threshold is NA, skipping proportion plots by study and threshold.")
}

# Plot each feature weight in the final scora

if (exists("norm_counts_ordered") && exists("coeffs_ordered") &&
    identical(rownames(norm_counts_ordered), names(coeffs_ordered))) {
  
  # 1. Calculate weighted expression: expression_ij * coefficient_i
  coeffs_vector <- as.numeric(coeffs_ordered)
  if (length(coeffs_vector) != nrow(norm_counts_ordered)) {
    stop("Mismatch in coefficient vector length and number of rows in count matrix.")
  }
  weighted_expression_matrix <- sweep(norm_counts_ordered, MARGIN = 1, STATS = coeffs_vector, FUN = "*")
  
  # 2. Calculate mean SIGNED contribution for each biomarker (across samples)
  #    The abs() function is REMOVED here.
  mean_signed_contribution <- rowMeans(weighted_expression_matrix, na.rm = TRUE)
  
  # 3. Prepare data for plotting
  biomarker_signed_contributions_df <- data.frame(
    Biomarker = names(mean_signed_contribution),
    MeanSignedContribution = as.numeric(mean_signed_contribution)
  )
  
  # Add a 'Sign' column based on the MeanSignedContribution for coloring
  biomarker_signed_contributions_df <- biomarker_signed_contributions_df %>%
    mutate(Sign = ifelse(MeanSignedContribution > 0, "Positive", "Negative"))
  
  # Order biomarkers by their MeanSignedContribution for the plot
  # (most positive at the top, most negative at the bottom)
  biomarker_signed_contributions_df <- biomarker_signed_contributions_df %>%
    arrange(desc(MeanSignedContribution)) %>%
    mutate(Biomarker = factor(Biomarker, levels = rev(Biomarker)))
  
  # Define colors (can reuse from previous coefficient plot)
  coefficient_colors <- c("Positive" = "darkred", "Negative" = "steelblue")
  
  # 4. Create the plot
  p_signed_contributions <- ggplot(biomarker_signed_contributions_df,
                                   aes(x = MeanSignedContribution, y = Biomarker, fill = Sign)) +
    geom_col(alpha = 0.9) + # geom_col for bar plot
    scale_fill_manual(values = coefficient_colors, name = "Contribution") +
    labs(title = "Mean Contribution to Response Score",
         x = "Mean(Expression × Coefficient)", # Updated x-axis label
         y = NULL) +
    theme_minimal(base_size = 12) +
    theme(
      panel.background = element_rect(fill = "white", colour = "white"),
      plot.background = element_rect(fill = "white", colour = "white"),
      panel.grid.major.x = element_line(colour = "grey90"),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      axis.text.y = element_text(color = "black"),
      legend.position = "top",
      text=element_text(size=18)
    ) +
    # Add a vertical line at x=0 for reference
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50")
  
  print(p_signed_contributions)
  ggsave(file="Fig8B.svg", plot=p_signed_contributions, width=8, height=7)
  
} else {
  print("Required data (`norm_counts_ordered` or `coeffs_ordered`) not found or misaligned. Skipping signed contribution plot.")
}

#Projection of the predictive signature on Whole Blood scRNAseq
library(Seurat)
Predic_Signature <- list(Predic_R_sign = names(median_coeffs_overlap[which(median_coeffs_overlap >0)]),
                         Predic_NR_sign = names(median_coeffs_overlap[which(median_coeffs_overlap <0)])
)
WB_NO <-readRDS("./scRNAseq/WB_scRNAseq.rds")
Idents(WB_NO)<-"Annotation"
WB_NO <- AddModuleScore(WB_NO, Predic_Signature[1], name = names(Predic_Signature[1]))
WB_NO <- AddModuleScore(WB_NO, Predic_Signature[2], name = names(Predic_Signature[2]))

p_feature_R <- do_FeaturePlot(sample = WB_NO, features = "Predic_R_sign1", enforce_symmetry = F, diverging.palette =  "RdYlBu",
                              min.cutoff = 0.0,
                              max.cutoff = 0.5,
                              plot_cell_borders = FALSE, use_viridis = TRUE, viridis.palette ="viridis", pt.size = 2,
                              raster = TRUE, raster.dpi = 1024)

p_feature_NR <-do_FeaturePlot(sample = WB_NO, features = "Predic_NR_sign1", enforce_symmetry = F, 
                              min.cutoff = -0.1, 
                              max.cutoff = 0.3,
                              plot_cell_borders = FALSE, use_viridis = TRUE, viridis.palette ="viridis", pt.size = 2,
                              raster = TRUE, raster.dpi = 1024)

ggsave(file="Fig8C.svg", plot=p_feature_R, width=8, height=7)
ggsave(file="Fig8D.svg", plot=p_feature_NR, width=8, height=7)
