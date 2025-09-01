library(forcats)
feature_to_plot <- "TN10_Baseline_NR_sign"
grouping_variable <- "Annotation"

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

n = 13
cols = gg_color_hue(n)

# dev.new(width = 4, height = 4)
# plot(1:n, pch = 16, cex = 2, col = cols)

plot_data <- FetchData(WB_NO, vars = c(feature_to_plot, grouping_variable))
colnames(plot_data) <- c("Expression", "Group")

names(cols)<- levels(plot_data$Group)

desired_order <- c("Neutrophils 1", "B cells", "Basophils", "Megakaryocytes", "NK/CD8 T cells", "Plasma cells", "T cells", "DCs", "Neutrophils 3", "Neutrophils 0", "Monocytes", "Neutrophils 4", "Neutrophils 2" )


plot_data$Group <- fct_relevel(plot_data$Group, desired_order)

#plot_data$Group <- factor(plot_data$Group, levels = desired_order)

# Verify the new order
print(paste("Reordered Cell Types:", paste(levels(plot_data$Group), collapse = ", ")))

# --- 3. Create the custom ggplot2 violin plot with manual color scale ---
ggplot(plot_data, aes(x = Group, y = Expression, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.7, scale = "width") +
  geom_point(position = position_jitter(width = 0.2, height = 0), size = 0.5, alpha = 0.5, color = "black") +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "grey30", alpha = 0.5) +
  stat_summary(fun = mean,
               geom = "point",
               shape = 23,
               size = 3,
               fill = "white",
               color = "red") +
  # --- Manual Color Scale - Using desired_order consistently ---
  scale_fill_manual(
    values = cols, # Use only as many colors as you have groups in your desired order
    breaks = desired_order                   # Ensure the order of breaks matches your desired order
  ) +
  labs(title = paste0("Expression of ", feature_to_plot, " by ", grouping_variable, " (Custom Colors)"),
       x = grouping_variable,
       y = paste0(feature_to_plot, " Expression")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "none" # Hide legend if x-axis labels are sufficient
  )



# Rename columns for clarity in ggplot


# --- 4. Create the custom ggplot2 violin plot ---
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


