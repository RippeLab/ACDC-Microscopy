######### Calculate Intensity ratios from NFkB Immunostaining ###########

# This script quantifies the intensities of the target protein NFkB. It computes the intensity across the cytoplasmic and nuclear cell mask.

#1 Before using the script: cellpose masks are generated for the total cell, cytoplasm, nucleus and background area. 
#2 Intensity signal of NFkB ("_Fused_637nm_corr.tif) is quantified in each mask in FIJI and data is saved in a CSV file (Source Data for Extended Data Fig. 6g: CSV file, Images and masks)
#3 CSV files are used as input in the following script for quantification of NFkB-p65 signal (Note: CST indicates the antibody company, Cell Signaling Technology)

# load libraries and set working directory
library(ggplot2)
library(dplyr)
library(tidyr)

setwd("/filepath/inputfolder")

# Define input files grouped by time point and antibody type
infiles <- list()
infiles[['0min_CST']] <- list(
  cyto = "CST/000min_p65_CST_cytoplasm.csv",
  cell = "CST/000min_p65_CST_cell.csv",
  bg = "CST/000min_p65_CST_background.csv",
  time = 0, 
  antibody = "CST"
)
infiles[['30min_CST']] <- list(
  cyto = "CST/030min_p65_CST_cytoplasm.csv",
  cell = "CST/030min_p65_CST_cell.csv",
  bg = "CST/030min_p65_CST_background.csv",
  time = 30, 
  antibody = "CST"
)
infiles[['100min_CST']] <- list(
  cyto = "CST/100min_p65_CST_cytoplasm.csv",
  cell = "CST/100min_p65_CST_cell.csv",
  bg = "CST/100min_p65_CST_background.csv",
  time = 100, 
  antibody = "CST"
)
infiles[['170min_CST']] <- list(
  cyto = "CST/170min_p65_CST_cytoplasm.csv",
  cell = "CST/170min_p65_CST_cell.csv",
  bg = "CST/170min_p65_CST_background.csv",
  time = 170, 
  antibody = "CST"
)
infiles[['240min_CST']] <- list(
  cyto = "CST/240min_p65_CST_cytoplasm.csv",
  cell = "CST/240min_p65_CST_cell.csv",
  bg = "CST/240min_p65_CST_background.csv",
  time = 240, 
  antibody = "CST"
)

# create empty data frame
data <- data.frame()

# loop each time point for analysis
for (infile in names(infiles)){
  data_in_cyto <- read.table(infiles[[infile]][['cyto']], header = T, sep = ",") #read cyto mask
  data_in_cell <- read.table(infiles[[infile]][['cell']], header = T, sep = ",") #read cell mask
  data_in_bg <- read.table(infiles[[infile]][['bg']], header = T, sep = ",") #read background mask
  bg <- data_in_bg$Mean
  antibody <- infiles[[infile]][['antibody']]
  time <- infiles[[infile]][['time']]
  data_in <- left_join(data_in_cell, data_in_cyto, by = "X.1", suffix = c("_cell", "_cyto")) #merge by id
  data_in <- data_in %>%
    filter((Area_cell - Area_cyto) > 0,
           Min_cell > 0,
           Min_cyto > 0) %>%
    select(-X_cell, -XM_cell, -Y_cell, -YM_cell, -X_cyto, -Y_cyto, -XM_cyto, -YM_cyto) %>%
    select(-Min_cyto, -Max_cyto, -Min_cell, -Max_cell, -Median_cell, -Median_cyto, -StdDev_cell, -StdDev_cyto)
  
  data_in <- data_in %>% 
    mutate(
      sample = infile,
      Cumulative_cyto = Mean_cyto*Area_cyto,
      Cumulative_cell = Mean_cell*Area_cell,
      Cumulative_nuc = Cumulative_cell - Cumulative_cyto
    )
  # add info 
  data_in$bg <- bg
  data_in$antibody <- antibody
  data_in$time <- time
  # combine
  data <- rbind(data, data_in)  
}

# Convert data to long format for faceting
data_long <- data %>%
  pivot_longer(cols = starts_with("Cumulative_"), 
               names_to = "Region", 
               values_to = "Cumulative_Value") %>%
  mutate(Region = recode(Region, 
                         "Cumulative_cyto" = "Cytoplasm", 
                         "Cumulative_cell" = "Cell", 
                         "Cumulative_nuc" = "Nucleus"))

# Plot violin + boxplot
ggplot(data_long, aes(x = factor(time, levels = c(0, 30, 100, 170, 240)), y = Cumulative_Value, fill = Region)) +
  geom_violin(alpha = 0.6) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +
  facet_wrap(Region ~ antibody, scales = "fixed", nrow = 3) +
  labs(
    y = "Cumulative Intensity", 
    x = "Time (min)", 
    title = "Cumulative Intensity Distribution per Sample"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(),
    strip.text = element_text(face = "bold")
  )

# Scatter plot
ggplot(data, aes(x = Cumulative_cyto, y = Cumulative_nuc, color = sample)) +
  geom_point(alpha = 0.6) +
  labs(x = "Cumulative Cytoplasm", y = "Cumulative Nucleus", title = "Cytoplasm vs. Nucleus Intensity") +
  theme_minimal()

# Add metrics
data <- data %>%
  mutate(log2FC_nuc_cyto = log2(Cumulative_nuc / Cumulative_cyto),
         #sample = factor(sample, levels = c("0min", "30min", "240min")),
         Area_nuc = Area_cell - Area_cyto,
         Mean_nuc = Cumulative_nuc/Area_nuc,
         log2FC_nuc_cyto_mean = log2(Mean_nuc / Mean_cyto),
         log2FC_nuc_cyto_mean_bg_corrected = log2((Mean_nuc - bg) / (Mean_cyto - bg)),
         Cumulative_nuc_bg_corrected = Cumulative_nuc - (bg*Area_nuc),
         log_Cumulative_nuc_bg_corrected = log10(Cumulative_nuc_bg_corrected),
         log_mean_nuc_bg_corrected = log10(Mean_nuc - bg)
         )

# Violin + boxplot: log2 mean intensity ratio (nucleus vs cytoplasm)
ggplot(data, aes(x = factor(time, levels = c(0, 30, 100, 170, 240)), y = log2FC_nuc_cyto_mean, fill = factor(time))) +
  geom_violin(alpha = 0.6) +  # Violin plot
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +  # Boxplot inside violin
  facet_wrap(~ antibody, nrow = 2) +
  labs(y = "Log2 Ratio (Nucleus / Cytoplasm)", x = "Sample", 
       title = "Log2 Ratio of Mean Intensity") +
  theme_minimal() +
  theme(axis.text.x = element_text())  # Rotate x-axis labels for clarity

# Violin + boxplot: log2 cumulative intensity ratio
ggplot(data, aes(x = factor(time, levels = c(0, 30, 100, 170, 240)), y = log2FC_nuc_cyto, fill = factor(time))) +
  geom_violin(alpha = 0.6) +  # Violin plot
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +  # Boxplot inside violin
  facet_wrap(~ antibody, nrow = 2) +
  labs(y = "Log2 Ratio (Nucleus / Cytoplasm)", x = "Sample", 
       title = "Log2 Ratio of Cumulative Intensity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for clarity

# Violin + boxplot: background-corrected cumulative nuclear intensity (log)
ggplot(data, aes(x = factor(time, levels = c(0, 30, 100, 170, 240)), y = log_Cumulative_nuc_bg_corrected, fill = factor(time))) +
  geom_violin(alpha = 0.6) +  # Violin plot
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +  # Boxplot inside violin
  facet_wrap(~ antibody, nrow = 2) +
  scale_y_log10() +
  labs(y = "log10 Cumulative Nuclear intensity", x = "Sample", 
       title = "Cumulative Nuclear Intensity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for clarity

### Plot for Figure ###

# Define your custom color palette
custom_colors <- c(
  "0"   = "#B0B0B0",
  "30"  = "#BDE3F5",
  "100" = "#FF0000",
  "170" = "#B3D1B3",
  "240" = "#877EB9"
)

# Save final plot 
pdf(file = "FigE6g_ratios_bg_corrected.pdf")
# Plot
ggplot(data %>% filter(antibody == "CST"), 
       aes(x = factor(time, levels = c(0, 30, 100, 170, 240)), 
           y = log2FC_nuc_cyto_mean_bg_corrected, fill = factor(time))) +
  geom_violin(alpha = 0.6, trim = FALSE) +
  geom_boxplot(width = 0.2, outlier.shape = NA, 
               color = "black", fill = "white", alpha = 1) +
  scale_fill_manual(values = custom_colors) +
  scale_y_continuous(
    breaks = seq(-2, 4, by = 1),
    expand = expansion(mult = c(0.05, 0.1))
  ) +
  labs(y = "Log2 Ratio (Nucleus - background / Cytoplasm - background)", x = "Sample", 
       title = "Log2 Ratio of Mean Intensity - CST") +
  theme_classic()
dev.off()
