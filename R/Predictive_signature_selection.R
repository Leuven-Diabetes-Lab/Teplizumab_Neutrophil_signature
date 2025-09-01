setwd("D:/aCD3_signature/Test_folder")

#Loading libraries:
library(glmnet)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(pROC)
library(parallel)
library(doParallel)

#Load the TN10/AbATE merged dataset
dds_merge <- readRDS("TN10_AbATE_merged_Full_Annotation.rds") 

#We keep only the baseline data from treated patient
to_keep<-which(dds_merge$Treatment %in% "Teplizumab")
dds_merge_Teplizumab <- dds_merge[,to_keep]
to_keep<-which(dds_merge_Teplizumab$Time_point %in% "Baseline")
dds_merge_Teplizumab_Baseline  <- dds_merge_Teplizumab[,to_keep]

# --- Parameters ---
n_resamples <- 200
alpha_sequence <- seq(0, 1, by = 0.1)
train_prop <- 0.7
n_bootstrap_iterations <- 100
cv_nfolds <- 5
lambda_choice <- "lambda.min"
feature_threshold <- 0.0000
top_n_features_within_alpha <- 33 # N for defining top features per alpha/resample step
top_n_for_alpha_signature <- 33 # N for final signature

# Detect cores for parallel CV
num_cores <- parallel::detectCores() - 1 # Use parallel:: explicitly if needed
if (is.na(num_cores) || num_cores < 1) {
  num_cores <- 1 # Default to 1 if detection fails
}
cat(sprintf("Planning to use %d cores for parallel cv.glmnet if applicable.\n", num_cores))


# --- Data Preparation ---
cat("Preparing full dataset...\n")
if (!exists("dds_merge_Teplizumab_Baseline") || !inherits(dds_merge_Teplizumab_Baseline, "DESeqDataSet")) {
  stop("DESeqDataSet object 'dds_merge_Teplizumab_Baseline' not found or is not the correct type.")
}
norm_counts <- counts(dds_merge_Teplizumab_Baseline, normalized = TRUE)

X_all <- t(norm_counts)
if (!is.matrix(X_all)) { X_all <- as.matrix(X_all) }

#Binary Outcome for Evaluation
Outcome_factor <- factor(dds_merge_Teplizumab_Baseline$Outcome)
if (nlevels(Outcome_factor) != 2) {
  stop("Outcome variable 'dds_merge_Teplizumab_Baseline$Outcome' must have exactly two levels (e.g., 'R', 'NR') for evaluation.")
}
# IMPORTANT: Confirm this mapping is correct the data (e.g., R=positive class)
# Create numeric version (0/1) for pROC
Y_all_binary_numeric <- ifelse(Outcome_factor == "R", 1, 0) # Define R=1, NR=0

sample_ids_all <- colnames(dds_merge_Teplizumab_Baseline)
cat("Full data dimensions (Samples x Features):", dim(X_all), "\n")
cat("Using binary 'Outcome' for training.\n")
cat("Using binary 'Outcome' (", levels(Outcome_factor)[1],"=0,", levels(Outcome_factor)[2], "=1) for evaluation.\n")

# --- Storage for Overall Results ---
all_resample_results_list <- list()
overall_selected_features <- list()
overall_coefficients <- list()
collected_cv_data_list <- list()

# --- SETUP PARALLEL BACKEND ---
cl <- NULL; if (num_cores > 1) {cat("Registering parallel...\n"); cl <- makeCluster(num_cores); registerDoParallel(cl); cat("Registered.\n")} else {cat("Running sequentially.\n")}

# --- Main Loops ---
cat(sprintf("Starting %d resampling iterations...\n", n_resamples))
set.seed(42); resample_seeds <- sample.int(1e6, n_resamples)

tryCatch({ 
  for (resample_iter in 1:n_resamples) {
    current_resample_seed<-resample_seeds[resample_iter]; set.seed(current_resample_seed); cat(sprintf("\n--- Resample %d / %d ---\n", resample_iter, n_resamples))
    
    # Train split
    train_indices<-sample(1:nrow(X_all),size=floor(nrow(X_all)*train_prop));
    y_train_binary_numeric<-Y_all_binary_numeric[train_indices];
    y_test_binary_numeric<-Y_all_binary_numeric[-train_indices];
    x_train<-X_all[train_indices, , drop=F];
    x_test<-X_all[-train_indices, , drop=F];
    
    # NA handling
    valid_train_idx<-!is.na(y_train_binary_numeric);
    x_train<-x_train[valid_train_idx, , drop=F];
    y_train_binary_numeric<-y_train_binary_numeric[valid_train_idx];
    
    valid_test_idx<-!is.na(y_test_binary_numeric);
    x_test<-x_test[valid_test_idx, , drop=F];
    y_test_binary_numeric<-y_test_binary_numeric[valid_test_idx];
    
    if(nrow(x_train)<2||nrow(x_test)<2||length(unique(na.omit(y_test_binary_numeric)))<2||nrow(x_train)<=cv_nfolds){cat(" Skip resample due to insufficient data after NA handling or for CV.\n");next}
    
    bootstrap_sample_size<-nrow(x_train); 
    eval_results_this_resample<-numeric(length(alpha_sequence)); 
    names(eval_results_this_resample)<-as.character(alpha_sequence); 
    feature_count_this_resample<-numeric(length(alpha_sequence)); 
    names(feature_count_this_resample)<-as.character(alpha_sequence); 
    optimal_lambda_this_resample<-numeric(length(alpha_sequence)); 
    names(optimal_lambda_this_resample)<-as.character(alpha_sequence); 
    selected_features_this_resample<-list(); 
    coefficients_this_resample<-list()
    
    cat(sprintf(" Alpha: ")); 
    for(alpha_idx in 1:length(alpha_sequence)){ 
      current_alpha<-alpha_sequence[alpha_idx]; 
      cat(sprintf("%.1f ", current_alpha)); 
      optimal_lambda<-NA; 
      
      cv_lambda_fit<-tryCatch({
        cv.glmnet(x_train, y_train_binary_numeric, 
                  alpha=current_alpha,
                  nfolds=cv_nfolds,
                  parallel=(num_cores>1),
                  type.measure="deviance", 
                  family = "binomial")     # Ensure family is binomial
      },error=function(e){
        # cat(sprintf("\nError in cv.glmnet for alpha %.2f: %s\n", current_alpha, conditionMessage(e))) # Optional detailed error
        NULL
      }); 
      
      if(is.null(cv_lambda_fit)){
        eval_results_this_resample[as.character(current_alpha)]<-NA;
        feature_count_this_resample[as.character(current_alpha)]<-0;
        optimal_lambda_this_resample[as.character(current_alpha)]<-NA;
        next # continue to next alpha
      }
      
      optimal_lambda<-ifelse(lambda_choice=="lambda.1se",cv_lambda_fit$lambda.1se,cv_lambda_fit$lambda.min); 
      optimal_lambda_this_resample[as.character(current_alpha)]<-optimal_lambda; 
      
      # Store detailed CV curve data for alpha >= 0.1 for aggregated plot
      if (current_alpha >= 0.1) { 
        if(!is.null(cv_lambda_fit$lambda) && !is.null(cv_lambda_fit$cvm) && !is.null(cv_lambda_fit$nzero)){
          cv_df_this_run <- data.frame(
            Lambda = cv_lambda_fit$lambda,
            LogLambda = log(cv_lambda_fit$lambda),
            MeanDeviance = cv_lambda_fit$cvm,
            UpperDeviance = cv_lambda_fit$cvup, # For plotting error bands if desired
            LowerDeviance = cv_lambda_fit$cvlo, # For plotting error bands if desired
            NZero = cv_lambda_fit$nzero,
            Alpha = current_alpha,
            ResampleID = resample_iter
          )
          collected_cv_data_list[[length(collected_cv_data_list) + 1]] <- cv_df_this_run
        } else {
          cat(sprintf("\nWarning: cv.glmnet output for alpha %.2f, resample %d missing some components (lambda, cvm, nzero).\n", current_alpha, resample_iter))
        }
      }

      all_coefficients_list_alpha<-list(); all_selected_features_names_alpha<-list(); 
      bootstrap_seeds<-sample.int(1e6,n_bootstrap_iterations); 
      for(i in 1:n_bootstrap_iterations){ 
        set.seed(current_resample_seed+i+alpha_idx*1000); 
        bootstrap_indices<-sample(1:nrow(x_train),size=bootstrap_sample_size,replace=T); 
        X_sample<-x_train[bootstrap_indices,,drop=F]; 
        Y_sample<- y_train_binary_numeric[bootstrap_indices];
        
        if(length(unique(na.omit(Y_sample)))<=1){next}; 
        constant_features<-apply(X_sample,2,function(x)length(unique(x))<=1); 
        if(any(constant_features)){next}
        
        lasso_model_fit<-tryCatch({
          glmnet(X_sample,Y_sample,alpha=current_alpha,lambda=optimal_lambda, family = "binomial")
        },error=function(e){NULL}); 
        if(is.null(lasso_model_fit)){next}; 
        coef_lasso<-coef(lasso_model_fit,s=optimal_lambda)[,1]; 
        all_run_coeffs<-coef_lasso[-1]; 
        significant_coeffs<-all_run_coeffs[abs(all_run_coeffs) > feature_threshold];
        if(length(significant_coeffs)>0){ 
          all_coefficients_list_alpha[[length(all_coefficients_list_alpha)+1]]<-significant_coeffs; 
          all_selected_features_names_alpha[[length(all_selected_features_names_alpha)+1]]<-names(significant_coeffs)
        }
      }
      
      current_eval_metric<-NA; 
      current_n_features<-0; 
      top_features_alpha<-character(0); 
      if(length(all_selected_features_names_alpha)>0){ 
        all_features_flat_alpha_boot<-unlist(all_selected_features_names_alpha); 
        feature_counts_alpha_boot<-table(all_features_flat_alpha_boot); 
        feature_counts_sorted_alpha_boot<-sort(feature_counts_alpha_boot,decreasing=T); 
        current_top_n<-min(top_n_features_within_alpha,length(feature_counts_sorted_alpha_boot)); 
        top_features_alpha<-names(head(feature_counts_sorted_alpha_boot,current_top_n)); 
        current_n_features<-length(top_features_alpha); 
        median_coeffs_alpha<-numeric(length(top_features_alpha)); 
        names(median_coeffs_alpha)<-top_features_alpha
        
        if(length(all_coefficients_list_alpha)>0){
          for(feature in top_features_alpha){
            coeffs_for_feature<-sapply(all_coefficients_list_alpha,function(rc){
              if(!is.null(rc)&&feature%in%names(rc))rc[feature] else NA
            }); 
            median_coeffs_alpha[feature]<-median(coeffs_for_feature,na.rm=T)
          }; 
          median_coeffs_alpha[is.na(median_coeffs_alpha)]<-0
        }
        else{
          median_coeffs_alpha[top_features_alpha]<-0}; 
        final_coeffs_vector<-numeric(ncol(x_test)); 
        names(final_coeffs_vector)<-colnames(x_test); 
        common_features_test<-intersect(names(median_coeffs_alpha),names(final_coeffs_vector)); 
        final_coeffs_vector[common_features_test]<-median_coeffs_alpha[common_features_test]; 
        test_scores<-tryCatch({as.numeric(x_test %*% final_coeffs_vector)},error=function(e){NA})
        
        if(!anyNA(test_scores)&&length(unique(y_test_binary_numeric))>1&&length(unique(test_scores))>1){
          roc_obj<-tryCatch({pROC::roc(response=y_test_binary_numeric,predictor=test_scores,quiet=T,direction="<")},error=function(e){NULL}); 
          if(!is.null(roc_obj)){current_eval_metric<-as.numeric(pROC::auc(roc_obj))} else {}} else {}} else {}
      
      eval_results_this_resample[as.character(current_alpha)]<-current_eval_metric; 
      feature_count_this_resample[as.character(current_alpha)]<-current_n_features; 
      selected_features_this_resample[[as.character(current_alpha)]]<-top_features_alpha; 
      coefficients_this_resample[[as.character(current_alpha)]]<-all_coefficients_list_alpha
    } # End alpha loop
    
    cat("\n"); 
    resample_df<-data.frame(ResampleID=resample_iter,Alpha=alpha_sequence,Optimal_Lambda=optimal_lambda_this_resample,AUC=eval_results_this_resample,N_Features=feature_count_this_resample); 
    all_resample_results_list[[resample_iter]]<-resample_df; 
    overall_selected_features[[resample_iter]]<-selected_features_this_resample; 
    overall_coefficients[[resample_iter]]<-coefficients_this_resample # Store coeffs
  } # End outer resampling loop
}, finally = { # --- STOP PARALLEL BACKEND ---
  if (!is.null(cl)) {cat("\nStopping parallel...\n"); stopCluster(cl); registerDoSEQ(); cat("Stopped.\n")}
})

# --- Aggregate and Analyze Performance Results ---

cat("\n--- Aggregating and Analyzing Performance Results ---\n"); # ... (Print summary_stats, plots) ...
results_full_df <- bind_rows(all_resample_results_list); 
results_plot_df <- results_full_df %>% filter(!is.na(AUC))

if(nrow(results_plot_df)>0){ 
  summary_stats<-results_plot_df%>%group_by(Alpha)%>%
    dplyr::summarize(Mean_AUC=mean(AUC, na.rm=T), Median_AUC=median(AUC, na.rm=T), SD_AUC=sd(AUC, na.rm=T), Mean_Optimal_Lambda=mean(Optimal_Lambda, na.rm=T), Median_Optimal_Lambda=median(Optimal_Lambda, na.rm=T), N_Resamples=n(), Mean_N_Features=mean(N_Features, na.rm=T), .groups='drop')%>%
    arrange(desc(Median_AUC)); 
  cat("\nSummary Statistics per Alpha (Sorted by Median AUC):\n"); 
  print(summary_stats)
} else { cat("No valid AUC results.\n")}

# --- Analyze Feature Selection Frequency ---
cat("\n--- Analyzing Feature Selection Frequency ---\n")
all_features_lists <- unlist(overall_selected_features, recursive = TRUE)
if(length(all_features_lists)>0){
  overall_feature_counts <- sort(table(all_features_lists), decreasing=TRUE); 
  cat("\nOverall Top 33 Features:\n"); 
  print(head(overall_feature_counts,33))
} else {
    cat("No features selected.\n"); 
  overall_feature_counts <- NULL
  }
feature_df_long_list<-list(); 
for(r_id in seq_along(overall_selected_features)){
  for(alpha_val in names(overall_selected_features[[r_id]])){
    features<-overall_selected_features[[r_id]][[alpha_val]]; 
    if(length(features)>0){
      feature_df_long_list[[length(feature_df_long_list)+1]]<-data.frame(ResampleID=r_id, Alpha=as.numeric(alpha_val), Feature=features)
    }
  }
  }
if(length(feature_df_long_list)>0){
  feature_df_long<-bind_rows(feature_df_long_list); 
  per_alpha_freq <- feature_df_long %>% 
    group_by(Alpha, Feature) %>% 
    summarise(Frequency=n(), .groups='drop') %>% 
    arrange(Alpha, desc(Frequency)); 
  cat("\nTop 5 Features per Alpha:\n"); 
  print(per_alpha_freq%>%group_by(Alpha)%>%top_n(5, Frequency))
} else {
    cat("No features selected to calc per-alpha freq.\n"); 
  per_alpha_freq <- NULL
  }


overlapping_features <- character(0) # Default

if (!is.null(per_alpha_freq) && nrow(per_alpha_freq) > 0) {
  
  # 1. Identify Top N per alpha and filter
  cat(sprintf("Identifying Top %d features per alpha & filtering...\n", top_n_features_within_alpha))
  top_features_per_alpha_list <- list()
  eligible_top_feature_lists <- list()
  included_alphas <- c()
  ###We remove the Ridge regression which gives results, but very different from the rest
  for(alpha_val in alpha_sequence[-1]) {
    top_n_this_alpha <- per_alpha_freq %>% filter(Alpha == alpha_val) %>%
      top_n(top_n_features_within_alpha, Frequency) %>% pull(Feature)
    top_features_per_alpha_list[[as.character(alpha_val)]] <- top_n_this_alpha # Keep for reference if needed
    
    # Filter for intersection eligibility
    if (length(top_n_this_alpha) >= top_n_features_within_alpha) {
      eligible_top_feature_lists[[as.character(alpha_val)]] <- top_n_this_alpha
      included_alphas <- c(included_alphas, alpha_val)
      # cat(sprintf("  Including Alpha %s (%d features)\n", alpha_val, length(top_n_this_alpha))) # Optional print
    } # else { cat(sprintf("  Excluding Alpha %s (%d features)\n", alpha_val, length(top_n_this_alpha))) }
  }
  cat(sprintf(" Found %d alphas meeting feature count threshold (>=%d).\n", length(eligible_top_feature_lists), top_n_features_within_alpha))
  
  
  # 2. Calculate Intersection
  if (length(eligible_top_feature_lists) >= 2) {
    overlapping_features <- Reduce(intersect, eligible_top_feature_lists)
    cat(sprintf(" Intersection yields %d features.\n", length(overlapping_features)))
  } else if (length(eligible_top_feature_lists) == 1) {
    cat(" Only 1 alpha eligible. Using its features directly.\n")
    overlapping_features <- eligible_top_feature_lists[[1]]
  } else {
    cat(" Fewer than 1 alpha eligible. No overlapping features.\n")
  }
  
  # --- Proceed only if overlapping features are found ---
  if(length(overlapping_features) > 0) {
    cat(" Found overlapping features:\n"); print(overlapping_features)
    
    # 3. Calculate Median Coefficients for overlapping features
    cat(" Calculating median coefficients for overlapping features...\n")
    median_coeffs_overlap<-numeric(length(overlapping_features)); 
    names(median_coeffs_overlap)<-overlapping_features
    for(feature in overlapping_features){ 
      all_coeffs_for_feature<-list(); 
      for(r_coeffs_list in overall_coefficients){ 
        for(alpha_val_char in names(r_coeffs_list)[as.numeric(names(r_coeffs_list)) %in% included_alphas]) { 
          a_coeffs_boot_list<-r_coeffs_list[[alpha_val_char]]; 
          for(run_coeffs in a_coeffs_boot_list){ 
            if(!is.null(run_coeffs)&&feature%in%names(run_coeffs)){ 
              all_coeffs_for_feature[[length(all_coeffs_for_feature)+1]]<-run_coeffs[feature]
            }
          }
        }
        }
    if(length(all_coeffs_for_feature)>0){
      median_coeffs_overlap[feature]<-median(unlist(all_coeffs_for_feature),na.rm=T)
      }else{
        median_coeffs_overlap[feature]<-0
        }
      }; 
    median_coeffs_overlap[is.na(median_coeffs_overlap)]<-0
    print(median_coeffs_overlap) 
    
    # 4. Calculate Scores for ALL Samples
    final_coeffs_vector_B<-numeric(ncol(X_all)); 
    names(final_coeffs_vector_B)<-colnames(X_all); 
    common_B<-intersect(names(median_coeffs_overlap), names(final_coeffs_vector_B)); 
    final_coeffs_vector_B[common_B]<-median_coeffs_overlap[common_B]; 
    scores_B<-tryCatch({as.numeric(X_all %*% final_coeffs_vector_B)}, error=function(e){NA})
    
    
    # 5. Plot Scores vs Binary Outcome
    if (!anyNA(scores_B)) {
      df_plot_B <- data.frame(SampleID = sample_ids_all, Score = scores_B, Outcome = Outcome_factor)
      subtitle_b <- sprintf("%d features overlapping from %d eligible Alphas", length(overlapping_features), length(eligible_top_feature_lists))
      p_score_B <- ggplot(df_plot_B, aes(x = Outcome, y = Score, fill = Outcome)) + 
        geom_boxplot(outlier.shape=NA, alpha=0.7)+
        geom_jitter(width=0.2, height=0, alpha=0.5, size=1.5)+
        stat_compare_means(method="wilcox.test", label.x.npc="center", label.y.npc=0.9)+
        labs(title="Scenario B: Scores using Overlapping Top Features (Filtered Alphas)", subtitle=subtitle_b, x="Outcome Outcome (R vs NR)", y="Calculated Score")+
        theme_bw()+theme(legend.position="none")
      print(p_score_B)
    } else { cat(" Could not calculate scores.\n") }
    
  } else { cat("No overlapping features found.\n") }
  
} else { 
  cat("No per-alpha feature frequency data available.\n")
}

merge_Signature_All_genes<-median_coeffs_overlap
write.csv2(merge_Signature_All_genes, file = "Predictive_Signature.csv")

# --- Plot CV Curve and Coefficient Paths for Overlapping Features  ---
# Check if overlapping_features were found and are not empty
if (exists("overlapping_features") && !is.null(overlapping_features) && length(overlapping_features) > 0) {
  
  cat(sprintf("Found %d overlapping features for path plot: %s\n",
              length(overlapping_features), paste(overlapping_features, collapse=", ")))
  
  # 1. Subset the full data to include only the overlapping features
  if(!all(overlapping_features %in% colnames(X_all))) {
    cat("ERROR: Some overlapping features are not in X_all column names. Aborting path plot.\n")
    missing_in_X_all <- overlapping_features[!overlapping_features %in% colnames(X_all)]
    cat("Missing features:", paste(missing_in_X_all, collapse=", "), "\n")
  } else {
    X_overlap_for_path <- X_all[, overlapping_features, drop = FALSE]
    Y_binary_for_path_plot <- Y_all_binary_numeric # From your main script's data prep
    
    # *** CRUCIAL: Ensure X_overlap_for_path has correct column names ***
    # Standard subsetting X_all[, overlapping_features] SHOULD preserve colnames if X_all has them
    # and overlapping_features are valid column names of X_all. This is an extra check/assignment.
    if (ncol(X_overlap_for_path) == length(overlapping_features)) {
      # Assign names from overlapping_features to ensure order and presence
      colnames(X_overlap_for_path) <- overlapping_features
      cat("Confirmed/assigned column names to X_overlap_for_path using 'overlapping_features'.\n")
    } else {
      cat("ERROR: Mismatch between ncol(X_overlap_for_path) and length(overlapping_features). Cannot reliably assign column names.\n")
      # Proceeding, but plot labels might be numeric.
    }
    
    # Remove samples with NA in Y_binary_for_path_plot if any, and align X
    valid_samples_idx <- !is.na(Y_binary_for_path_plot)
    # Create final matrices to be used by glmnet
    X_overlap_for_path_final <- X_overlap_for_path[valid_samples_idx, , drop = FALSE]
    Y_binary_for_path_plot_final <- Y_binary_for_path_plot[valid_samples_idx]
    
    cat(sprintf("Dimensions of matrix for final path model: %d samples, %d features.\n",
                nrow(X_overlap_for_path_final), ncol(X_overlap_for_path_final)))
    
    
    if (nrow(X_overlap_for_path_final) < max(2, cv_nfolds) || # Need enough samples for CV if used, and for glmnet
        ncol(X_overlap_for_path_final) == 0 ||
        length(unique(Y_binary_for_path_plot_final)) < 2) {
      cat("Not enough data, features, or outcome levels to fit a binomial model for path plot.\n")
    } else {
      # 2. Choose an alpha for this final model visualization
      best_performing_alpha <- 0.5 # Default alpha
      if (exists("summary_stats") && is.data.frame(summary_stats) && nrow(summary_stats) > 0 &&
          "Alpha" %in% names(summary_stats) && "Median_AUC" %in% names(summary_stats)) {
        # Assuming summary_stats is sorted by Median_AUC descending
        best_performing_alpha <- summary_stats$Alpha[1]
        cat(sprintf("Using best performing alpha (%.2f based on Median_AUC) for path plot.\n", best_performing_alpha))
      } else {
        cat(sprintf("Using default alpha (%.2f) for path plot (summary_stats not found or not definitive).\n", best_performing_alpha))
      }
      
      # 3. Fit a new glmnet model on these features across a lambda sequence
      cat("Fitting final glmnet binomial model on overlapping features to generate paths...\n")
      final_path_model <- tryCatch({
        glmnet(X_overlap_for_path_final, Y_binary_for_path_plot_final,
               alpha = best_performing_alpha,
               family = "binomial") # Ensure family is binomial
      }, error = function(e) {
        cat("ERROR fitting final binomial model for path plot:", conditionMessage(e), "\n")
        NULL
      })
      
    }
  }
} else {
  cat("Skipping coefficient path plot: No overlapping features were identified in Scenario B or 'overlapping_features' not defined.\n")
}

if (!is.null(final_path_model)) {
  library(tidyr)
  lambda_seq <- final_path_model$lambda
  coef_matrix <- as.matrix(final_path_model$beta) # Feature coefficients
  
  # Ensure coef_matrix has rownames
  if(is.null(rownames(coef_matrix)) && ncol(X_overlap_for_path_final) == nrow(coef_matrix)){
    rownames(coef_matrix) <- colnames(X_overlap_for_path_final)
  }
  
  plot_df <- as.data.frame(t(coef_matrix)) %>%
    mutate(Lambda = lambda_seq, Log_Lambda = log(lambda_seq)) %>%
    pivot_longer(cols = -c(Lambda, Log_Lambda), names_to = "Feature", values_to = "Coefficient") %>%
    # Calculate L1 norm for each lambda if xvar = "norm" was intended
    group_by(Lambda) %>%
    mutate(L1_Norm = sum(abs(Coefficient))) %>%
    ungroup()
  
  # Get last non-zero point for labeling
  label_df <- plot_df %>%
    filter(Coefficient != 0) %>% # Or some small threshold
    group_by(Feature) %>%
    filter(L1_Norm == max(L1_Norm)) %>% # Label at largest L1 norm where coeff is non-zero
    slice(1) # Take one if multiple lambdas have max L1 norm
  
  p_manual <- ggplot(plot_df, aes(x = L1_Norm, y = Coefficient, group = Feature, color = Feature)) +
    geom_line(alpha=0.7) +
    ggrepel::geom_text_repel(data = label_df, aes(label = Feature),
                             size = 4, nudge_x = 0.1, segment.alpha=0.3,
                             direction="y", hjust=0) + # Added direction and hjust
    theme_bw() +
    labs(
      #title = plot_title,
      x = "L1 Norm", y = "Coefficients") +
    theme(legend.position = "none",text=element_text(size=20))
  print(p_manual)
}

ggsave(file="FigS7A.svg", plot=p_manual, width=8, height=8)

# --- Aggregated Cross-Validation Plot (Alphas 0.1 to 1.0) *** ---
cat("\n--- Generating Aggregated Cross-Validation Plot ---\n")

if (length(collected_cv_data_list) > 0) {
  full_cv_data <- bind_rows(collected_cv_data_list)
  
  # 1. Create a common LogLambda grid for interpolation
  min_log_lambda <- min(full_cv_data$LogLambda, na.rm = TRUE)
  max_log_lambda <- max(full_cv_data$LogLambda, na.rm = TRUE)
  common_log_lambda_grid <- seq(min_log_lambda, max_log_lambda, length.out = 100) # 100 points
  
  # 2. Interpolate and average for each alpha
  averaged_cv_curves_list <- list()
  alphas_to_plot <- unique(full_cv_data$Alpha) # Should be 0.1 to 1.0
  
  cat("Interpolating and averaging CV curves for each alpha (0.1-1.0)...\n")
  for (a_val in alphas_to_plot) {
    cat(sprintf("  Processing Alpha = %.2f\n", a_val))
    alpha_subset_data <- full_cv_data %>% dplyr::filter(Alpha == a_val)
    
    interpolated_deviances_this_alpha <- list()
    interpolated_nzeros_this_alpha <- list()
    
    for (r_id in unique(alpha_subset_data$ResampleID)) {
      resample_alpha_data <- alpha_subset_data %>%
        dplyr::filter(ResampleID == r_id) %>%
        arrange(LogLambda) # Ensure data is sorted by LogLambda for approx
      
      # Interpolate MeanDeviance
      interp_dev <- approx(resample_alpha_data$LogLambda,
                           resample_alpha_data$MeanDeviance,
                           xout = common_log_lambda_grid,
                           rule = 2)$y # rule=2: use endpoint values for extrapolation
      interpolated_deviances_this_alpha[[as.character(r_id)]] <- interp_dev
      
      # Interpolate NZero
      interp_nzero <- approx(resample_alpha_data$LogLambda,
                             resample_alpha_data$NZero,
                             xout = common_log_lambda_grid,
                             rule = 2)$y
      interpolated_nzeros_this_alpha[[as.character(r_id)]] <- interp_nzero
    }
    
    if (length(interpolated_deviances_this_alpha) > 0) {
      # Convert lists of vectors to matrices and calculate rowMeans
      matrix_deviances <- do.call(cbind, interpolated_deviances_this_alpha)
      mean_dev_for_alpha <- rowMeans(matrix_deviances, na.rm = TRUE)
      
      matrix_nzeros <- do.call(cbind, interpolated_nzeros_this_alpha)
      mean_nzero_for_alpha <- rowMeans(matrix_nzeros, na.rm = TRUE)
      
      averaged_cv_curves_list[[as.character(a_val)]] <- data.frame(
        LogLambda = common_log_lambda_grid,
        MeanDeviance = mean_dev_for_alpha,
        MeanNZero = round(mean_nzero_for_alpha), # Round for display
        Alpha = factor(a_val)
      )
    }
  }
  
  if (length(averaged_cv_curves_list) > 0) {
    averaged_cv_df <- bind_rows(averaged_cv_curves_list)
    
    # 3. Create the plot
    p_aggregated_cv <- ggplot(averaged_cv_df, aes(x = LogLambda, y = MeanDeviance, color = Alpha)) +
      geom_line(linewidth = 1) + # Use linewidth instead of size for lines
      labs(
        #title = "Aggregated Cross-Validated Binomial Deviance vs. Log(Lambda)",
        #subtitle = paste("Averaged across", n_resamples, "resamples. Curves for Alphas 0.1 to 1.0."),
        x = "Log(Lambda)",
        y = "Average Cross-Validated Binomial Deviance",
        color = "Alpha"
      ) +
      theme_bw() +
      theme(legend.position = "top",text=element_text(size=25))
    
    ggsave(file="Fig8A.svg", plot=p_aggregated_cv, width=7, height=7)

    print(p_aggregated_cv)
    
    cat("\nSummary of points at minimum average deviance for each alpha curve:\n")
    print(avg_lambda_min_points %>% select(Alpha, LogLambda, MeanDeviance, MeanNZero) %>% arrange(Alpha))
    
  } else {
    cat("Could not generate averaged CV curves. No data processed.\n")
  }
  
} else {
  cat("No detailed CV data was collected. Skipping aggregated CV plot.\n")
}
