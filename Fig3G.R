##### Fig3G
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(patchwork)

# Define region of interest
chr <- "8"
start_pos <- 127000000
end_pos <- 128200000

# List of cell lines and file paths
samples <- c("A549", "CAKI2", "HCT116", "HepG2", "PANC1", "T47D", "G401", "RPMI_7951")
files <- c("A549_4DNFID68JQY9_score.bedgraph",
           "CAKI2_4DNFIJIPWD63_score.bedgraph",
           "HCT116_4DNFIXTAS6EE_score.bedgraph",
           "HepG2_4DNFICSTCJQZ_score.bedgraph",
           "PANC1_4DNFIFDGVWLU_score.bedgraph",
           "T47D_4DNFITSIPCSK_score.bedgraph", 
           "G401_4DNFI5VO3E1W_score.bedgraph", 
           "RPMI_7951_4DNFIOAA3ZEQ_score.bedgraph")

# Read, filter, and label each dataset
all_insulation <- data.frame()

for (i in seq_along(files)) {
  df <- read.table(files[i], header = FALSE)
  colnames(df) <- c("chr", "start", "end", "score")
  df <- df[df$chr == chr & df$start >= start_pos & df$end <= end_pos, ]
  df$cell_line <- samples[i]
  all_insulation <- rbind(all_insulation, df)
}

# Define gene annotations as a data frame
genes <- data.frame(
  gene = c("CCAT1", "MYC", "PVT1"),
  start = c(127207382, 127735434, 127794524),
  end = c(127219268, 127736623, 128101256),
  y = max(all_insulation$score) + 0.3,    
  ymin = -1,                              
  ymax = 2                              
)

# Define your custom colors
custom_colors <- setNames(brewer.pal(length(samples), "Dark2"), samples)

# Plot
ggplot(all_insulation, aes(x = start, y = score, color = cell_line)) +
  geom_line(size = 1) +
  
  # Add gene rectangles
  geom_rect(data = genes, inherit.aes = FALSE,
            aes(xmin = start, xmax = end, ymin = ymin, ymax = ymax, fill = gene),
            alpha = 0.2, color = "black")+
  # Add gene labels
  geom_text(data = genes, inherit.aes = FALSE,
            aes(x = (start + end) / 2, y = y, label = gene),
            angle = 0, size = 4, vjust = -0.5) +
  theme_minimal() +
  scale_color_manual(values = custom_colors)+
  labs(title = "Insulation scores from Hi-C",
       x = "Genomic position on Chr8",
       y = "Insulation Index",
       color = "Cell Line",
       fill = "Gene") +
  theme(axis.text.x = element_text(), 
        axis.title = element_text(size = 13), 
        plot.title = element_text(face = 'bold', size = 15, hjust =0.5))

# Summarize by genomic position
summary_df <- all_insulation %>%
  group_by(start) %>%
  summarise(
    mean_score = mean(score, na.rm = TRUE),
    sd_score = sd(score, na.rm = TRUE)
  )
mean_plot <- ggplot(summary_df, aes(x = start)) +
  geom_line(aes(y = mean_score), color = "black", size = 1) +
  geom_ribbon(aes(ymin = mean_score - sd_score, ymax = mean_score + sd_score), 
              fill = "grey70", alpha = 0.5) +
  theme_minimal() +
  labs(x = "Genomic position on Chr8", 
       y = "Mean ± SD", 
       title = "Mean Insulation Index with Variability") +
  theme(axis.title = element_text(size = 13),
        plot.title = element_text(face = 'bold', size = 15, hjust = 0.5))

# Original plot assigned to a variable
original_plot <- ggplot(all_insulation, aes(x = start, y = score, color = cell_line)) +
  geom_line(size = 1) +
  geom_rect(data = genes, inherit.aes = FALSE,
            aes(xmin = start, xmax = end, ymin = ymin, ymax = ymax, fill = gene),
            alpha = 0.2, color = "black") +
  geom_text(data = genes, inherit.aes = FALSE,
            aes(x = (start + end) / 2, y = y, label = gene),
            angle = 0, size = 4, vjust = -0.5) +
  theme_minimal() +
  scale_color_manual(values = custom_colors) +
  labs(title = "Insulation scores from Hi-C",
       x = "Genomic position on Chr8",
       y = "Insulation Index",
       color = "Cell Line",
       fill = "Gene") +
  theme(axis.title = element_text(size = 13),
        plot.title = element_text(face = 'bold', size = 15, hjust = 0.5))

# Combine
combined_plot <- original_plot / mean_plot + plot_layout(heights = c(2, 1))
print(combined_plot)
