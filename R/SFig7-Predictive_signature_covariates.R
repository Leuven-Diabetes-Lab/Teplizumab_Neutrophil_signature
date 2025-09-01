###First, run Fig8-Predictive_signature_validation.R###

# Create data frame including aligned Study info
score_data <- data.frame(
  Score = response_scores,
  Outcome = dds_merge_Teplizumab_Baseline$Outcome,
  Gender = dds_merge_Teplizumab_Baseline$Gender,
  Age = dds_merge_Teplizumab_Baseline$age_group,
  DR3= dds_merge_Teplizumab_Baseline$DR3,
  DR4= dds_merge_Teplizumab_Baseline$DR4 ,
  DQ8= dds_merge_Teplizumab_Baseline$DQ8 ,
  GAD65 = dds_merge_Teplizumab_Baseline$GAD65,
  IA2 = dds_merge_Teplizumab_Baseline$IA2,
  MIAA = as.factor( dds_merge_Teplizumab_Baseline$MIAA),
  ZNT8 = dds_merge_Teplizumab_Baseline$ZNT8,
  Study = factor(study_info_aligned) 
)

#Gender
my_comparisons <- list( c('Female', 'Male') )
outcome_colors <- c("pink", "blue")
names(outcome_colors) <- unlist(my_comparisons)

p_score_outcome <- ggplot(score_data, aes(x = Gender, y = Score, fill = Gender)) +
  geom_boxplot(width = 0.4, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(color = Study), 
              width = 0.15,       
              height = 0,
              alpha = 1,
              size = 2.5) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.format", size = 6) + 
  scale_fill_manual(values = outcome_colors) + 
  labs(y = "Response Score") +
  theme_minimal() +theme(text=element_text(size=25))

print(p_score_outcome) 
ggsave(file="FigS7D.svg", plot=p_score_outcome, width=7, height=7)

#Age
my_comparisons <- list( c("Children","Adolescent","Adult") )
outcome_colors <- c("#009E73", "#56B4E9","#F0E442")
names(outcome_colors) <- unlist(my_comparisons)
my_comparison1 <-list( c("Children","Adolescent") )
my_comparison2 <-list( c("Adolescent","Adult") )
my_comparison3 <-list( c("Children","Adult") )


p_score_outcome <- ggplot(score_data, aes(x = Age, y = Score, fill = Age)) +
  geom_boxplot(width = 0.4, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(color = Study), 
              width = 0.15,       
              height = 0,
              alpha = 1,
              size = 2.5) +
  stat_compare_means(comparisons = my_comparison1, method = "wilcox.test", label = "p.format", size = 6) + 
  stat_compare_means(comparisons = my_comparison2, method = "wilcox.test", label = "p.format", size = 6) + 
  stat_compare_means(comparisons = my_comparison3, method = "wilcox.test", label = "p.format", size = 6, label.y = 0.9) + 
  scale_fill_manual(values = outcome_colors) + 
  labs(y = "Response Score") +
  theme_minimal() +theme(text=element_text(size=25))

print(p_score_outcome) 
ggsave(file="FigS7E.svg", plot=p_score_outcome, width=7, height=7)

#HLA DR3 haplotype
my_comparisons <- list( c('ABSENT', 'PRESENT') )
outcome_colors <- c("grey50", "red")
names(outcome_colors) <- unlist(my_comparisons)

p_score_outcome <- ggplot(score_data[complete.cases(score_data),], aes(x = DR3, y = Score, fill = DR3)) +
  geom_boxplot(width = 0.4, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(color = Study), 
              width = 0.15,       
              height = 0,
              alpha = 1,
              size = 2.5) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.format", size = 6) + 
  scale_fill_manual(values = outcome_colors) + 
  labs(y = "Response Score") +
  theme_minimal() +theme(text=element_text(size=25))

print(p_score_outcome) 
ggsave(file="FigS7F.svg", plot=p_score_outcome, width=7, height=7)

#HLA DR4 haplotype
my_comparisons <- list( c('ABSENT', 'PRESENT') )
outcome_colors <- c("grey50", "red")
names(outcome_colors) <- unlist(my_comparisons)

p_score_outcome <- ggplot(score_data[complete.cases(score_data),], aes(x = DR4, y = Score, fill = DR4)) +
  geom_boxplot(width = 0.4, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(color = Study), 
              width = 0.15,       
              height = 0,
              alpha = 1,
              size = 2.5) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.format", size = 6) + 
  scale_fill_manual(values = outcome_colors) + 
  labs(y = "Response Score") +
  theme_minimal() +theme(text=element_text(size=25))


print(p_score_outcome) 
ggsave(file="FigS7G.svg", plot=p_score_outcome, width=7, height=7)

#HLA DQ8 haplotype
my_comparisons <- list( c('ABSENT', 'PRESENT') )
outcome_colors <- c("grey50", "red")
names(outcome_colors) <- unlist(my_comparisons)

p_score_outcome <- ggplot(score_data[complete.cases(score_data),], aes(x = DQ8, y = Score, fill = DQ8)) +
  geom_boxplot(width = 0.4, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(color = Study), 
              width = 0.15,       
              height = 0,
              alpha = 1,
              size = 2.5) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.format", size = 6) + 
  scale_fill_manual(values = outcome_colors) + 
  labs(y = "Response Score") +
  theme_minimal() +theme(text=element_text(size=25))

print(p_score_outcome) 
ggsave(file="FigS7H.svg", plot=p_score_outcome, width=7, height=7)

###And now, the Auto Antibodies###

#GAD65
my_comparisons <- list( c('0', '1') )
outcome_colors <- c("grey50", "red")
names(outcome_colors) <- unlist(my_comparisons)

p_score_outcome <- ggplot(score_data[complete.cases(score_data),], aes(x = GAD65, y = Score, fill = GAD65)) +
  geom_boxplot(width = 0.4, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(color = Study), 
              width = 0.15,     
              height = 0,
              alpha = 1,
              size = 2.5) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.format", size = 6) + 
  scale_fill_manual(values = outcome_colors) + 
  labs(y = "Response Score") +
  theme_minimal() +theme(text=element_text(size=25))

print(p_score_outcome) 
ggsave(file="FigS7I.svg", plot=p_score_outcome, width=7, height=7)


#IA2
my_comparisons <- list( c('0', '1') )
outcome_colors <- c("grey50", "red")
names(outcome_colors) <- unlist(my_comparisons)

p_score_outcome <- ggplot(score_data[complete.cases(score_data),], aes(x = IA2, y = Score, fill = IA2)) +
  geom_boxplot(width = 0.4, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(color = Study),
              width = 0.15,   
              height = 0,
              alpha = 1,
              size = 2.5) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.format", size = 6) + 
  scale_fill_manual(values = outcome_colors) + 
  labs(y = "Response Score") +
  theme_minimal() +theme(text=element_text(size=25))

print(p_score_outcome) 
ggsave(file="FigS7J.svg", plot=p_score_outcome, width=7, height=7)


#MIAA
my_comparisons <- list( c('0', '1') )
outcome_colors <- c("grey50", "red")
names(outcome_colors) <- unlist(my_comparisons)

p_score_outcome <- ggplot(score_data[complete.cases(score_data),], aes(x = MIAA, y = Score, fill = MIAA)) +
  geom_boxplot(width = 0.4, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(color = Study), 
              width = 0.15,     
              height = 0,
              alpha = 1,
              size = 2.5) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.format", size = 6) + 
  scale_fill_manual(values = outcome_colors) + 
  labs(y = "Response Score") +
  theme_minimal() +theme(text=element_text(size=25))

print(p_score_outcome) 
ggsave(file="FigS7K.svg", plot=p_score_outcome, width=7, height=7)


#ZNT8
my_comparisons <- list( c('0', '1') )
outcome_colors <- c("grey50", "red")
names(outcome_colors) <- unlist(my_comparisons)

p_score_outcome <- ggplot(score_data[complete.cases(score_data),], aes(x = ZNT8, y = Score, fill = ZNT8)) +
  geom_boxplot(width = 0.4, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(color = Study), 
              width = 0.15,   
              height = 0,
              alpha = 1,
              size = 2.5) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.format", size = 6) + 
  scale_fill_manual(values = outcome_colors) + 
  labs(y = "Response Score") +
  theme_minimal() +theme(text=element_text(size=25))

print(p_score_outcome) 
ggsave(file="FigS7L.svg", plot=p_score_outcome, width=7, height=7)


