######### CALCULATE NUCLEAR INTENSITY and SIGNAL DISTRBUTION from IMMUNOFLUORESCENCE DATA ###########

# This script quantifies the intensities of target gene/protein. It computes nuclear intensity and distribution across conditions 
# It uses the "quantNuclei" custom function to compute the sum of intensity in the pre-segmented cellpose nuclear mask
# See Source Data for Extended Data Fig. 6f (CSV file with intensities, Microscopy images with NFkB channel and nuclear masks)

# Input data (stitched images):
#1 cellpose nuclear masks segmented from DAPI channel (_corr_cp_masks.png)
#2 (optional) ROIs masks to exclude from the analysis (e.g. laser overexposure), that. have been pre-selected, binarized and inverted in FIJI (_ROIs_masks_inverted.png)
#3 channel of interest in TIF format for target gene/protein (e.g., "_Fused_637nm_corr.tif")

########################### Load required packages and functions ##########################
library(EBImage)
library(tidyr)
library(dplyr)
library(ggplot2)
source("/filepath/quantNuclei_v01.R")
###########################################################################################

####### PART 1: retrieve data #######

# Assign working directory (without "/" in the end)
input_dir <- "/filepath/R_input/"
output_dir <- "/filepath/R_output/"

# Assign experiment condition 
time_point <- "240min"
experiment <- "p65_CST"   #target and antibody used
exp_date <- "ymd"         #assign date of experiment in ymd format

# Retrieve files in the directory 
all_files <- list.files(input_dir, recursive = TRUE, full.names = TRUE)

# Select files for nuclear masks (ROIs = areas to exclude from the original nuclear mask)
nucmasks_image <- all_files[grep("_corr_cp_masks.png$",all_files)]
#nucmasks_ROIs <- all_files[grep("mask_inverted.tif$", all_files)]
# Load images
nucmasks_ori <- readImage(nucmasks_image)
#nucmasks_exclude <- readImage(nucmasks_ROIs)

# Reformat nuclear masks (necessary for cellpose mask)
nucmasks_ori <- Image(as.numeric(as.factor(nucmasks_ori))-1,dim = dim(nucmasks_ori))
#nucmasks_exclude <- Image(as.numeric(as.factor(nucmasks_exclude))-1,dim = dim(nucmasks_exclude))

display(nucmasks_ori)
#display(nucmasks_exclude)

###### PART 2: FILTERING STEP 1 FOR NUCLEAR MASKS. Exclude regions of interest of the nuclear mask for the analysis (overexposure, stitching errors...) ######

# 2.1 Filter out masks in "nucmasks_exclude" from the original "nukmasks_ori" image (only if nucmasks_exclude is present in this condition)
#nucmasks_filt1 <- nucmasks_ori * nucmasks_exclude
#display(nucmasks_filt1)

# 2.2 Rename filtered image according to experiment and time point. Save in output folder
#nucmasks_filt1_rename <- paste0("00_", exp_date, "_HUVEC_", experiment, "_", time_point, "_nucmasks_filt1.tiff")
#writeImage(nucmasks_filt1, file.path(output_dir, nucmasks_filt1_rename),  type = "tiff", bits.per.sample = 16)

## NB: If Filtering Step 1 was applied, re-assign filtered image to nucmasks
#nucmasks <- nucmasks_filt1

## NB: If Filtering Step 1 was NOT applied, re-assign original image to nucmasks
nucmasks <- nucmasks_ori


###### PART 3: FILTERING STEP 2 FOR NUCLEAR MASKS.  ######
# 3.1 Remove nuclei that are cut at the edges of the image
dims<-dim(nucmasks)
border_size <- 10    # Set the border size (border of 1 pixel is NOT enough) - use 10 or more
top_row <- nucmasks[1:border_size,]
left_column <- nucmasks[,1:border_size]
bottom_row <- nucmasks[(dims[1]-border_size):dims[1],]
right_column <- nucmasks[,(dims[2]-border_size):dims[2]]
border <- c(top_row,left_column,bottom_row,right_column)    # Assign border
ids <- unique(border[which(border != 0)])
# Transform image as integer before removing objects, then change back to double
storage.mode(nucmasks) <- 'integer' 
nucmasks_filt2 <- rmObjects((nucmasks), ids)
storage.mode(nucmasks_filt2) <- 'double'

# 3.2 Inspect filtered image. 
display(nucmasks_filt2)   # If nuclei at borders are still included, repeat from: nucmasks <- nucmasks_filt1 OR nucmask <- nucmask_ori and change border_size (e.g., 15-20)

# 3.3 Rename and save final image (include "/2^16 max pixel value so that the nuclei IDs will be displayed in FIJI; include +1 to start cell count from ID=1)
nucmasks_filt2_rename <- paste0("01_", exp_date, "_HUVEC_", experiment, "_", time_point, "_nucmasks_final.tiff")    #rename image according to experiment conditions
writeImage((nucmasks_filt2 + 1) / 2^16, file.path(output_dir, nucmasks_filt2_rename),  type = "tiff", bits.per.sample = 16)

#OPTIONAL: check number of cells in the final mask
Ncell_initial <- length(unique(as.vector(nucmasks_filt2)))
print(Ncell_initial)

# 3.4 Assign filtered image to nuclei_masks before proceeding 
nuclei_masks <- nucmasks_filt2

###### PART 4: QUANTIFICATION OF SIGNAL IN NUCLEAR MASKS ######

# 4.1 Assign channel for the quantification 
IF_image <- all_files[grep("_Fused_637nm_corr.tif$",all_files)]
IF_channel <- readImage(IF_image)

# 4.2 use quantNuclei function to calculate sum intensity across nuclear masks 
p65_IntegrInt <- quantNuclei(nuclei_masks, IF_channel, sum) %>% data.frame() %>% setNames('IntegratedIntensity')

# 4.3 create data frame with x and y coordinates of nuclear masks (and area) 
nucfeat <- cbind(computeFeatures.moment(nuclei_masks),computeFeatures.shape(nuclei_masks))
nucfeat <- data.frame("x"=round(nucfeat[,"m.cx"]),"y"=round(nucfeat[,"m.cy"]),"nucleus_area"=nucfeat[,"s.area"])
nucfeat$maskID <- as.integer(rownames(nucfeat))
p65_IntegrInt$maskID <- as.integer(rownames(p65_IntegrInt))
fulldf <- full_join(nucfeat, p65_IntegrInt, by="maskID")  # bind nucfeat and gene_IntegrInt by the maskID

# 4.4 Rename according to experiment condition
IntegrInt_rename <- paste0("02_", exp_date, "_HUVEC_", experiment, "_", time_point, "_sum_intensity_table.csv")
fulldf_rename <- paste0("03_", exp_date, "_HUVEC_", experiment, "_", time_point, "_sum_intensity_table_XY_full.csv")
# 4.5 save CSV files
write.csv(p65_IntegrInt, file = file.path(output_dir, IntegrInt_rename), row.names = TRUE)
write.csv(fulldf, file = file.path(output_dir, fulldf_rename), row.names = TRUE)


###### PART 5: INSPECT AND PLOT SIGNAL INTENSITY ######

library(ggplot2)

# 5.1 Assign data (csv file) to be inspected
int_table <- p65_IntegrInt
#int_table <- read.csv("/filepath/output/02_date_HUVEC_p65_CST_240min_sum_int_table.csv")
hist(int_table$'IntegratedIntensity')

# Optional: rename and save plot
int_plot_rename <- paste0("04_", exp_date, "_HUVEC_", experiment, "_", time_point, "_intensity_plot.pdf")
plot_path <- file.path(output_dir, int_plot_rename)
pdf(plot_path)
hist(int_table$IntegratedIntensity)
dev.off()

#### repeat steps above for each condition/timepoint (generate all CSV files before proceeding with the next steps)

# 5.2 Retrieve intensity table for each condition
int_table_0h <- read.csv("/filepath/02_HUVEC_p65_CST_0h_intensity_table.csv")
int_table_30min <- read.csv("/filepath/02_HUVEC_p65_CST_30min_intensity_table.csv")
int_table_100min <- read.csv("/filepath/02_HUVEC_p65_CST_100min_intensity_table.csv")
int_table_170min <- read.csv("/filepath/02_HUVEC_p65_CST_170min_intensity_table.csv")
int_table_240min <- read.csv("/filepath/02_HUVEC_p65_CST_240min_intensity_table.csv")

# 5.3 Create data frame with all conditions
p65_0h_df <- int_table_0h %>% select(IntegratedIntensity)
colnames(p65_0h_df)<- "IntegratedIntensity"
p65_0h_df$condition <- "1.WT" 

p65_30min_df <- int_table_30min %>% select(IntegratedIntensity)
colnames(p65_30min_df)<- "IntegratedIntensity"
p65_30min_df$condition <- "2.TNFa+30min" 

p65_100min_df <- int_table_100min %>% select(IntegratedIntensity)
colnames(p65_100min_df)<- "IntegratedIntensity"
p65_100min_df$condition <- "3.TNFa+100min" 

p65_170min_df <- int_table_170min %>% select(IntegratedIntensity)
colnames(p65_170min_df)<- "IntegratedIntensity"
p65_170min_df$condition <- "4.TNFa+170min" 

p65_240min_df <- int_table_240min %>% select(IntegratedIntensity)
colnames(p65_240min_df)<- "IntegratedIntensity"
p65_240min_df$condition <- "5.TNFa+240min" 

df_conditions <- rbind(p65_0h_df, p65_30min_df, p65_100min_df, p65_170min_df, p65_240min_df)

# 5.4 Plot and inspect conditions as preferred
boxplot <- ggplot(df_conditions, aes(x=condition, y=IntegratedIntensity))+geom_boxplot()+labs(y= "Integrated Intensity")
boxplot

violin_boxplot <- ggplot(df_conditions, aes(x = condition, y = IntegratedIntensity)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.2, fill = "white") +
  labs(y = "Integrated Intensity") +
  theme_minimal()
violin_boxplot

log10_boxplot <- ggplot(df_conditions, aes(x=condition, y=log(IntegratedIntensity, 10)))+geom_boxplot()+labs(y= "Integrated Intensity")
log10_boxplot

log10_violin_boxplot <- ggplot(df_conditions, aes(x = condition, y = log10(IntegratedIntensity))) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.3, fill = "white") +
  labs(y = "Integrated Intensity (log10)") +
  theme_minimal()
log10_violin_boxplot

###### PART 6: PLOT INTENSITY DISTRIBUTION ######

library(ggplot2)
library(viridis) 
library(viridisLite)

# 6.1 Inspect intensity distributions
hist(int_table_0h$'IntegratedIntensity')
hist(int_table_30min$'IntegratedIntensity')
hist(int_table_100min$'IntegratedIntensity')
hist(int_table_170min$'IntegratedIntensity')
hist(int_table_240min$'IntegratedIntensity')

# 6.2 Calculate bins and bin width for each condition. Calculate sqrt(n) = square root on the number of observation (cells) to determine the bin to be used
calculate_binwidth <- function(data) {
  n_bins <- round(sqrt(nrow(data)))
  bin_width <- (max(data$IntegratedIntensity) - min(data$IntegratedIntensity)) / n_bins
  return(list(n_bins = n_bins, bin_width = bin_width))
}

bins_0h <- calculate_binwidth(int_table_0h)
bins_30min <- calculate_binwidth(int_table_30min)
bins_100min <- calculate_binwidth(int_table_100min)
bins_170min <- calculate_binwidth(int_table_170min)
bins_240min <- calculate_binwidth(int_table_240min)

# 6.3 Define function to plot histogram with density curve overlay
plot_histogram_density <- function(data, bins, title, color) {
  ggplot(data, aes(x = IntegratedIntensity)) +
    geom_histogram(aes(y = ..density..), binwidth = bins$bin_width, fill = color, color = "black", alpha = 0.5) +
    geom_density(color = "black", size = 0.8) +  # Add density curve
    labs(title = title, x = "Nuclear Intensity", y = "Density") +
    theme_minimal() }
colors <- viridis(5)

# 6.4 Generate and inspect individual plots
bin_plot_0h <- plot_histogram_density(int_table_0h, bins_0h, "Untreated", colors[1])
bin_plot_30min <- plot_histogram_density(int_table_30min, bins_30min, "30 min TNFα", colors[2])
bin_plot_100min <- plot_histogram_density(int_table_100min, bins_100min, "100 min TNFα", colors[3])
bin_plot_170min <- plot_histogram_density(int_table_170min, bins_170min, "170 min TNFα", colors[4])
bin_plot_240min <- plot_histogram_density(int_table_240min, bins_240min, "240 min TNFα", colors[5])

bin_plot_0h
bin_plot_30min
bin_plot_100min
bin_plot_170min
bin_plot_240min

# 6.5 Combine all conditions into a final density plot
combine_conditions <- rbind(
  transform(int_table_0h, Condition = "Untreated"),
  transform(int_table_30min, Condition = "30 min TNFα"),
  transform(int_table_100min, Condition = "100 min TNFα"),
  transform(int_table_170min, Condition = "170 min TNFα"),
  transform(int_table_240min, Condition = "240 min TNFα")
)

# 6.6 Generate density plot
density_plot <- ggplot(combine_conditions, aes(x = IntegratedIntensity, color = Condition, fill = Condition)) +
  geom_density(alpha = 0.3, size = 0.8) +  # Density curves
  scale_color_viridis_d() +  # Better color scheme
  scale_fill_viridis_d() +
  labs(title = "NFkB-p65 Intensity Distribution", x = "Nuclear Intensity", y = "Density") +
  theme_minimal() +
  coord_cartesian(xlim = c(0, 2000))  #Set x axis limit accordingly
print(binplot_density)

# 6.7 Rename and Save density plot
rename_density <- paste0("05_", exp_date, "_HUVEC_", experiment, "_final_density_plot.pdf")
ggsave(paste0("/filepath/", rename_density), density_plot, dpi=300)

