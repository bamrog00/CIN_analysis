library(ggplot2)
library(stats)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(RColorBrewer)
library(tidyr)


###### The following functions are to load the data, it is quite convoluted

loadCancerData = function(version, cancertype, mode){
  patient_level = read.csv(paste('../data/CIN_Output_version_',version,'/',cancertype,'/',cancertype,'_level_patient.csv', sep = ''), sep = ',', row.names = 1)
  chromosome_level = read.csv(paste('../data/CIN_Output_version_',version,'/',cancertype,'/',cancertype,'_level_chromosome.csv', sep = ''), sep = ',', row.names = 1)
  if (mode == 'Arm' | mode == 'All'){
    arm_level = read.csv(paste('../data/CIN_Output_version_',version,'/',cancertype,'/',cancertype,'_level_arms.csv', sep = ''), sep = ',', row.names = 1)
  }
  if (mode == 'Cytoband' | mode == 'All'){
    cytoband_level = read.csv(paste('../data/CIN_Output_version_',version,'/',cancertype,'/',cancertype,'_level_cytobands.csv', sep = ''), sep = ',', row.names = 1)
  }
  
  if (mode == 'Arm'){
    return (list(patient_level, chromosome_level, arm_level))
  }
  if (mode == 'Cytoband'){
    return (list(patient_level, chromosome_level, cytoband_level))
  }
  if (mode == 'All'){
    return (list(patient_level, chromosome_level, arm_level, cytoband_level))
  }
  
  return (list(patient_level, chromosome_level))
}


loadCytobandHRR_Data = function(version, cancertype){
  patient_HRR = read.csv(paste('../data/CIN_Output_version_',version,'/',cancertype,'/',cancertype,'_level_patient_HRR.csv', sep = ''), sep = ',', row.names = 1)
  cytoband_HRR = read.csv(paste('../data/CIN_Output_version_',version,'/',cancertype,'/',cancertype,'_level_cytobands_HRR.csv', sep = ''), sep = ',', row.names = 1)
  return (list(patient_HRR, cytoband_HRR))
}

rename_file_column <- function(dataframe, old, new) {
  if (old %in% names(dataframe)) {
    names(dataframe)[names(dataframe) == old] <- new
  }
  return(dataframe)
}

getCancerData = function(version, cancertype, mode, HRR){
  
  results = loadCancerData(version, cancertype, mode)
  
  patient_level = results[[1]]
  chromosome_level = results[[2]]
  
  ## Only needed for test_All
  if ('File.name' %in% names(chromosome_level)){
    chromosome_level = rename_file_column(chromosome_level, 'File.name', 'File.Name')
  }
  
  merged_data = merge(patient_level, chromosome_level, by = 'File.Name')
  
  if (mode == 'Arm'){
    arm_level = results[[3]]
    if ('File.name' %in% names(arm_level)){
      arm_level = rename_file_column(arm_level, 'File.name', 'File.Name')
    }
    merged_data = merge(merged_data, arm_level, by = 'File.Name')
  }
  
  if (mode == 'Cytoband'){
    cytoband_level = results[[4]]
    if ('File.name' %in% names(cytoband_level)){
      cytoband_level = rename_file_column(cytoband_level, 'File.name', 'File.Name')
    }
    merged_data = merge(merged_data, cytoband_level, by = 'File.Name')
  }
  
  if (mode == 'All'){
    arm_level = results[[3]]
    if ('File.name' %in% names(arm_level)){
      arm_level = rename_file_column(arm_level, 'File.name', 'File.Name')
    }
    merged_data = merge(merged_data, arm_level, by = 'File.Name')
    cytoband_level = results[[4]]
    if ('File.name' %in% names(cytoband_level)){
      cytoband_level = rename_file_column(cytoband_level, 'File.name', 'File.Name')
    }
    merged_data = merge(merged_data, cytoband_level, by = 'File.Name')
  }
  
  if ((mode == 'All' | mode == 'Cytoband') & HRR){
    results_HRR = loadCytobandHRR_Data(version, cancertype)
    pateint_hrr = results_HRR[[1]]
    cytoband_hrr = results_HRR[[2]]
    if ('File.name' %in% names(pateint_hrr)){
      pateint_hrr = rename_file_column(pateint_hrr, 'File.name', 'File.Name')
    }
    if ('File.name' %in% names(cytoband_hrr)){
      cytoband_hrr = rename_file_column(cytoband_hrr, 'File.name', 'File.Name')
    }
    merged_data = merge(merged_data, pateint_hrr, by = 'File.Name')
    merged_data = merge(merged_data, cytoband_hrr, by = 'File.Name')
  }
  
  return (merged_data)
}

loadAllCohorts = function(version, mode){
  cancertypes = getSubfolders(paste('../data/CIN_Output_version_', version, sep = ''))
  no = c('TARGET-CCSK','TARGET-ALL-P2', 'TARGET-AML','TCGA-CESC', 'TCGA-CHOL', 'TCGA-COAD', 'TCGA-DLBC', 'TCGA-ESCA', 'TCGA-GBM', 'TCGA-LAML', 'TCGA-MESO', 'TCGA-PCPG', 'TCGA-READ', 'TCGA-SKCM', 'TCGA-TGCT', 'TCGA-THCA', 'TCGA-THYM', 'TCGA-UCS', 'TCGA-UVM')
  cancertypes = cancertypes[!cancertypes %in% no]
  print(length(cancertypes))
  for (i in 1:length(cancertypes)){
    cancertype = cancertypes[i]
    
    if (i == 1){
      data = getCancerData(version, cancertype, mode, FALSE)
      data['ProjectID'] = cancertype
    }else{
      other_data = getCancerData(version, cancertype, mode, FALSE)
      other_data['ProjectID'] = cancertype
      data = rbind(data,other_data)
    }
    print(paste('Done cohort ', cancertype, sep = ''))
  }
  return (data)
  
}

getSubfolders = function(folder_path) {
  subfolders = list.dirs(folder_path, recursive = FALSE)
  subfolder_names <- basename(subfolders)
  return(subfolder_names)
}

annotateHRDType = function(data){
  cutoff_data = read.csv('../../GMM_Cutoff/data/kmeans_cutoffs_pancancer.csv', sep = ',')
  
  data$HRDtype = ifelse(data$HRD_sum >= cutoff_data$n2_cutoff[match(data$ProjectID, cutoff_data$Project.ID)], "High", "Low")
  return (data)
}

loadAllCohortsCytobands = function(version){
  cancertypes = getSubfolders(paste('../data/CIN_Output_version_', version, sep = ''))

  for (i in 1:length(cancertypes)){
    cancertype = cancertypes[i]
    
    if (i == 1){
      results_HRR = loadCytobandHRR_Data(version, cancertype)
      patient_level = read.csv(paste('../data/CIN_Output_version_',version,'/',cancertype,'/',cancertype,'_level_patient.csv', sep = ''), sep = ',', row.names = 1)
      pateint_hrr = results_HRR[[1]]
      cytoband_hrr = results_HRR[[2]]
      merged = merge(patient_level, pateint_hrr, by = 'File.Name')
      data = merge(merged, cytoband_hrr, by = 'File.Name')
      data['ProjectID'] = cancertype
    }else{
      results_HRR = loadCytobandHRR_Data(version, cancertype)
      patient_level = read.csv(paste('../data/CIN_Output_version_',version,'/',cancertype,'/',cancertype,'_level_patient.csv', sep = ''), sep = ',', row.names = 1)
      pateint_hrr = results_HRR[[1]]
      cytoband_hrr = results_HRR[[2]]
      merged = merge(patient_level, pateint_hrr, by = 'File.Name')
      other_data = merge(merged, cytoband_hrr, by = 'File.Name')
      other_data['ProjectID'] = cancertype
      data = rbind(data,other_data)
    }
    print(paste('Done cohort ', cancertype, sep = ''))
  }
  return (data)
  
}


##### Load the data used for thesis and presentation
data_all = loadAllCohorts('2_0', 'All')
data_test = data.frame(data_all)
data_test_type = annotateHRDType(data_test)




no = c('TARGET-CCSK','TARGET-ALL-P2', 'TARGET-AML','TCGA-CESC', 'TCGA-CHOL', 'TCGA-COAD', 'TCGA-DLBC', 'TCGA-ESCA', 'TCGA-GBM', 'TCGA-LAML', 'TCGA-MESO', 'TCGA-PCPG', 'TCGA-READ', 'TCGA-SKCM', 'TCGA-TGCT', 'TCGA-THCA', 'TCGA-THYM', 'TCGA-UCS', 'TCGA-UVM')
data_17_cohorts = data_test_type[!(data_test_type$ProjectID %in% no), ]

data_17_cohorts$totalCIN = data_17_cohorts$general_cin + data_17_cohorts$cn_cin




##### Load cytoband data
cytoband_all = loadAllCohortsCytobands('2_0')
cytoband_test = data.frame(cytoband_all)
cytoband_18 = cytoband_test[!(cytoband_test$ProjectID %in% no), ]
cytoband_18 = annotateHRDType(cytoband_18)

cytoband_18$totalCIN_HRR = cytoband_18$general_cin_HRR + cytoband_18$cn_cin_HRR




##### Some of the figures such as boxplot HRD high and HRD low boxplot for CIN results
general_cin_all_type = ggplot(data_test_type, aes(x = ProjectID, y = general_cin, fill = HRDtype)) +
  geom_boxplot() +
  labs(x = "Project ID", y = "General CIN", fill = "HRD-Type")  +
  ggtitle("General CIN, K-means") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  stat_compare_means(aes(group = HRDtype), method = 't.test', label = 'p.signif') +
  scale_fill_manual(values = c("High" = "#1F77B4", "Low" = "#E377C2"), 
                    labels = c("HRD-high", "HRD-low"))
general_cin_all_type
ggsave('../data/figures/boxplot_general_cin_all_type.png',general_cin_all_type,width = 8.63, height = 5.71)

cn_cin_all_type = ggplot(data_test_type, aes(x = ProjectID, y = cn_cin, fill = HRDtype)) +
  geom_boxplot() +
  labs(x = "Project ID", y = "CN-CIN", fill = "HRD-Type")  +
  ggtitle("CN-CIN, K-means") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  stat_compare_means(aes(group = HRDtype), method = 't.test', label ="p.signif") +
  scale_fill_manual(values = c("High" = "#1F77B4", "Low" = "#E377C2"), 
                    labels = c("HRD-high", "HRD-low"))

ggsave('../data/figures/boxplot_cn_cin_all_type.png',cn_cin_all_type,width = 8.63, height = 5.71)


general_cin_all = ggplot(data = data_test_type, aes(x = ProjectID, y = general_cin))+
  geom_boxplot()+
  ggtitle("General CIN score pan-cancer") +
  labs(x = "Project ID", y = "CIN score")  +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

ggsave('../data/figures/boxplot_general_cin_all.png',general_cin_all,width = 8.63, height = 5.71)


cn_cin_all = ggplot(data = data_test_type, aes(x = ProjectID, y = cn_cin))+
  geom_boxplot()+
  ggtitle("CN-CIN score pan-cancer") +
  labs(x = "Project ID", y = "CIN score")  +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

ggsave('../data/figures/boxplot_cn_cin_all.png',cn_cin_all,width = 8.63, height = 5.71)




general_cin_17_type = ggplot(data_17_cohorts, aes(x = ProjectID, y = general_cin, fill = HRDtype)) +
  geom_boxplot() +
  labs(x = "Project ID", y = "General CIN", fill = "HRD-Type")  +
  ggtitle("General CIN, K-means") +
  theme_bw()+
  ylim(0, 500)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  stat_compare_means(aes(group = HRDtype), method = 'wilcox.test', label ="p.signif") +
  scale_fill_manual(values = c("High" = "#1F77B4", "Low" = "#E377C2"), 
                    labels = c("HRD-high", "HRD-low"))

ggsave('../data/figures/boxplot_general_cin_17_type_scaled.png',general_cin_17_type,width = 8.63, height = 5.71)

cn_cin_17_type = ggplot(data_17_cohorts, aes(x = ProjectID, y = cn_cin, fill = HRDtype)) +
  geom_boxplot() +
  labs(x = "Project ID", y = "CN-CIN", fill = "HRD-Type")  +
  ggtitle("CN-CIN, K-means") +
  theme_bw()+
  ylim(0, 500)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  stat_compare_means(aes(group = HRDtype), method = 'wilcox.test', label ="p.signif") +
  scale_fill_manual(values = c("High" = "#1F77B4", "Low" = "#E377C2"), 
                    labels = c("HRD-high", "HRD-low"))

ggsave('../data/figures/boxplot_cn_cin_17_type_scaled.png',cn_cin_17_type,width = 8.63, height = 5.71)

total_cin_17_type = ggplot(data_17_cohorts, aes(x = ProjectID, y = totalCIN, fill = HRDtype)) +
  geom_boxplot() +
  labs(x = "Project ID", y = "Total CIN", fill = "HRD-Type")  +
  ggtitle("Total CIN, K-means") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  stat_compare_means(aes(group = HRDtype), method = 't.test', label ="p.signif") +
  scale_fill_manual(values = c("High" = "#1F77B4", "Low" = "#E377C2"), 
                    labels = c("HRD-high", "HRD-low"))

ggsave('../data/figures/boxplot_total_cin_17_type.png',total_cin_17_type,width = 8.63, height = 5.71)



general_cin_17 = ggplot(data = data_17_cohorts, aes(x = ProjectID, y = cn_cin))+
  geom_boxplot()+
  ggtitle("General CIN score pan-cancer") +
  labs(x = "Project ID", y = "CIN score")  +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

ggsave('../data/figures/boxplot_general_cin_17.png',general_cin_17,width = 8.63, height = 5.71)


cn_cin_17 = ggplot(data = data_17_cohorts, aes(x = ProjectID, y = cn_cin))+
  geom_boxplot()+
  ggtitle("General CIN score pan-cancer") +
  labs(x = "Project ID", y = "CIN score")  +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

ggsave('../data/figures/boxplot_cn_cin_17.png',cn_cin_17,width = 8.63, height = 5.71)


#######################################
#######################################
#######################################
#######################################



pattern <- "^chr[1-9][0-9]?_average_cn$"


subset_df <- data_17_cohorts[, c("ProjectID", "HRD_sum", "TAI", "LST", "LOH", grep(pattern, names(data_17_cohorts), value = TRUE))]

subsubset_df = subset_df[,c(grep(pattern, names(data_17_cohorts), value = TRUE))]

column_numbers <- sub("^chr([1-9][0-9]?).*", "\\1",  grep(pattern, names(data_17_cohorts), value = TRUE))

# Rename the matching columns


####### All copy number


new_column_names <- paste0("chr", column_numbers)
names(subsubset_df)[1:length(subsubset_df)] <- new_column_names

numbers <- 1:22
strings <- paste0("chr", numbers)

df_reorder <- subset(subsubset_df, select = strings)


z_scored = sapply(df_reorder, function(df_reorder) (df_reorder-mean(df_reorder))/sd(df_reorder))
matrix = as.matrix(z_scored)



matrix_t = t(matrix)

max_avg_cn = max(matrix)

col_fun = colorRamp2(c(0, 2, 7), c("blue", "white", "red"))

ha = HeatmapAnnotation(Cancertype = subset_df$ProjectID)

Heatmap(matrix_t, cluster_rows = FALSE, cluster_columns = FALSE, top_annotation  = ha, show_column_names  = FALSE, raster_quality = 10, column_split = subset_df$ProjectID, column_title=NULL, border = TRUE,  column_gap = unit(0, "mm"), name = 'Average copy number')


##### Mean average copy number


new_df <- data_17_cohorts %>%
  group_by(ProjectID) %>%
  summarise(
    mean_chr1 = mean(chr1_average_cn, na.rm = TRUE),
    mean_chr2 = mean(chr2_average_cn, na.rm = TRUE),
    mean_chr3 = mean(chr3_average_cn, na.rm = TRUE),
    mean_chr4 = mean(chr4_average_cn, na.rm = TRUE),
    mean_chr5 = mean(chr5_average_cn, na.rm = TRUE),
    mean_chr6 = mean(chr6_average_cn, na.rm = TRUE),
    mean_chr7 = mean(chr7_average_cn, na.rm = TRUE),
    mean_chr8 = mean(chr8_average_cn, na.rm = TRUE),
    mean_chr9 = mean(chr9_average_cn, na.rm = TRUE),
    mean_chr10 = mean(chr10_average_cn, na.rm = TRUE),
    mean_chr11 = mean(chr11_average_cn, na.rm = TRUE),
    mean_chr12 = mean(chr12_average_cn, na.rm = TRUE),
    mean_chr13 = mean(chr13_average_cn, na.rm = TRUE),
    mean_chr14 = mean(chr14_average_cn, na.rm = TRUE),
    mean_chr15 = mean(chr15_average_cn, na.rm = TRUE),
    mean_chr16 = mean(chr16_average_cn, na.rm = TRUE),
    mean_chr17 = mean(chr17_average_cn, na.rm = TRUE),
    mean_chr18 = mean(chr18_average_cn, na.rm = TRUE),
    mean_chr19 = mean(chr19_average_cn, na.rm = TRUE),
    mean_chr20 = mean(chr20_average_cn, na.rm = TRUE),
    mean_chr21 = mean(chr21_average_cn, na.rm = TRUE),
    mean_chr22 = mean(chr22_average_cn, na.rm = TRUE),
  )

numbers <- 1:22
strings <- paste0("chr", numbers)
names(new_df)[2:length(new_df)] = strings


subsetnew_df = new_df[,strings]

matrix_2 = as.matrix(subsetnew_df)

matrix_2_t = t(matrix_2)

max_avg_cn_2 = max(matrix_2)

col_fun = colorRamp2(c(0, 2, max_avg_cn_2), c("blue", "white", "red"))

ha = HeatmapAnnotation(Cancertype = new_df$ProjectID)

Heatmap(matrix_2_t, cluster_rows = FALSE, cluster_columns = FALSE, top_annotation  = ha, show_column_names  = FALSE, raster_quality = 10, column_split = new_df$ProjectID, column_title=NULL, border = TRUE,  column_gap = unit(0, "mm"), col = col_fun, name = 'Mean average CopyNum')
Heatmap(matrix_2, cluster_rows = FALSE, cluster_columns = FALSE, col = col_fun)








data_17_cohorts$Type_binary <- ifelse(data_17_cohorts$HRDtype == 'High', 1, 0)

scores = c('general_cin', 'cn_cin', 'n_cn_loh', 'n_cn_gain', 'n_cn_amp', 'n_homo_del', 'n_hemi_del', 'n_gain', 'n_amp')

ids = unique(data_17_cohorts$ProjectID)

patient_df <- data_17_cohorts[, c("ProjectID", "HRD_sum", 'Type_binary', 'general_cin', 'cn_cin', 'n_cn_loh',  'n_cn_gain', 'n_cn_amp', 'n_homo_del', 'n_hemi_del', 'n_gain', 'n_amp' )]

empty_df <- data.frame(matrix(ncol = length(ids), nrow = length(scores)))
colnames(empty_df) = ids
rownames(empty_df) = scores

for (score in scores){
  for (cancer in ids){
    subdata = patient_df[patient_df['ProjectID'] == cancer,]
    correlation = cor.test(subdata$Type_binary,subdata[,score], method = 'pearson')
    empty_df[score,cancer] = correlation$estimate
  }
}

corrmatrix = as.matrix(empty_df)


Heatmap(corrmatrix, show_column_dend = FALSE, show_row_dend = FALSE, name = 'Point-biserial correlation', column_title = 'Correlation CIN measrues level Patients and HRD-Type')





pattern_avg <- "^chr[1-9][0-9]?_average_cn$"
pattern_general <- "^chr[1-9][0-9]?_general_cin$"
pattern_cn <- "^chr[1-9][0-9]?_cn_cin$"



chromosome_subdf = data_17_cohorts[, c("ProjectID", "HRD_sum", grep(pattern_avg, names(data_17_cohorts), value = TRUE), grep(pattern_general, names(data_17_cohorts), value = TRUE), grep(pattern_cn, names(data_17_cohorts), value = TRUE))]

ids = unique(chromosome_subdf$ProjectID)

measures = c(grep(pattern_avg, names(chromosome_subdf), value = TRUE), grep(pattern_general, names(chromosome_subdf), value = TRUE), grep(pattern_cn, names(chromosome_subdf), value = TRUE))

empty_df_chr <- data.frame(matrix(ncol = length(ids), nrow = length(measures)))
colnames(empty_df_chr) = ids
rownames(empty_df_chr) = measures

for (measure in measures){
  for (cancer in ids){
    subdata = chromosome_subdf[chromosome_subdf['ProjectID'] == cancer,]
    correlation = cor(subdata$HRD_sum,subdata[,measure])
    empty_df_chr[measure,cancer] = correlation
  }
}

corrmatrix_chromosome = as.matrix(empty_df_chr)


Heatmap(corrmatrix_chromosome,cluster_rows = FALSE, row_names_gp = gpar(fontsize = 6), show_column_dend = FALSE, show_row_dend = FALSE, name = 'Pearson', column_title = 'Correlation CIN measrues level Chromosome and HRD-score')





scores = c('general_cin', 'cn_cin', 'n_cn_loh', 'n_cn_gain', 'n_cn_amp', 'n_homo_del', 'n_hemi_del', 'n_gain', 'n_amp')

ids = unique(data_test_type$ProjectID)

patient_df <- data_test_type[, c("ProjectID", "HRD_sum", 'general_cin', 'cn_cin', 'n_cn_loh',  'n_cn_gain', 'n_cn_amp', 'n_homo_del', 'n_hemi_del', 'n_gain', 'n_amp' )]

empty_df <- data.frame(matrix(ncol = length(ids), nrow = length(scores)))
colnames(empty_df) = ids
rownames(empty_df) = scores

for (score in scores){
  for (cancer in ids){
    subdata = patient_df[patient_df['ProjectID'] == cancer,]
    correlation = cor(subdata$HRD_sum,subdata[,score], method = 'spearman')
    empty_df[score,cancer] = correlation
  }
}

corrmatrix = as.matrix(empty_df)


Heatmap(corrmatrix, show_column_dend = FALSE, show_row_dend = FALSE, name = 'Spearman correlation', column_title = 'Correlation CIN measrues level Sample and HRDsum')






###### Colors for figures (Heatmaps)
num_colors <- 18

# Generate a palette of visually distinct colors
color_palette <- colorRampPalette(brewer.pal(12, "Paired"))(num_colors)

# Create a named color palette
named_palette <- setNames(color_palette, unique(data_test_type$ProjectID))

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
sub_col_vector = col_vector[1:num_colors]
named_palette <- setNames(sub_col_vector, unique(data_17_cohorts$ProjectID))

named_palette <- setNames(sub_col_vector, unique(cytoband_18$ProjectID))



################# Level: Chromosome ###### Plotting heatmaps



plotHeatmapChromosome = function(data_measurement, variable_list, color_vector, level, measure_full_name, measure_short_name, save = FALSE){
  
  subset_data = data_measurement[,variable_list]
  
  
  colnames(subset_data) <- sapply(colnames(subset_data), function(col) {
    parts <- unlist(strsplit(col, '_'))
    
    if (level == 'Arm') {
      # Keep first two parts and join with a '.'
      new_col <- paste(parts[1:2], collapse = '.')
    } else {
      # Keep only the first part
      new_col <- parts[1]
    }
    
    return(new_col)
  })
  

  
  chr_order <- paste0("chr", 1:22)
  if (level == 'Arm') {
    chr_order <- rep(chr_order, each = 2)
    chr_order <- paste0(chr_order, c(".p", ".q"))
  }
  
  subset_data <- subset_data[, colnames(subset_data)[order(match(colnames(subset_data), chr_order))]]
  
  z_scored = sapply(subset_data, function(subset_data) (subset_data-mean(subset_data))/sd(subset_data))
  matrix = as.matrix(z_scored)
  rownames(matrix) = data_measurement$ProjectID
  matrix_t = t(matrix)
  
  
  ha = HeatmapAnnotation(Cancertype = data_measurement$ProjectID, HRDType = data_measurement$HRDtype, col = list(Cancertype = color_vector,HRDType = c('High' = 'orange', 'Low' = 'blue')))
  if (level == 'Arm'){
    size_row = 7
    size_col = 8
  }else{
    size_row = 12
    size_col = 9
  }
  
  hm = Heatmap(matrix_t, cluster_rows = FALSE, cluster_columns = TRUE, show_column_dend = FALSE,
          top_annotation  = ha, show_column_names  = FALSE, raster_quality = 10,
          column_split = data_measurement$ProjectID, column_title= paste(measure_full_name,', Level ',level, sep = ''),
          column_gap = unit(0, "mm"),row_names_gp = gpar(fontsize = size_row),column_names_gp = gpar(fontsize = size_col), name = 'Z-Score')
  
  if (save){
    png(file=paste('../data/figures/heatmap_',measure_short_name,'_',level, '.png', sep = ''))
  }
  
  draw(hm)
  
  if (save){
    dev.off()
  }
}


plotmeanHeatmapChromosome = function(data_measurement, variable_list, color_vector, level, measure_full_name, measure_short_name, save = FALSE){
  
  chr_order <- paste0("chr", 1:22)
  chr_arm_order = chr_order
  if (level == 'Arm') {
    chr_order <- rep(chr_order, each = 2)
    chr_arm_order <- rep(chr_order, each = 2)
    chr_arm_order <- paste0(chr_order, c("_p", "_q"))
    chr_order <- paste0(chr_order, c(".p", ".q"))
    
  }
  
  mean_df <- data.frame(ProjectID = character(),
                          HRDtype = character(),
                          stringsAsFactors = FALSE)
  
  
  for (i in 1:length(chr_order)){
    
    specific_level = chr_order[[i]]
    if (level == 'Arm'){
      variable = variable_list[startsWith(variable_list, paste(chr_arm_order[[i]],'_', sep = ''))]
    }else{
      variable = variable_list[startsWith(variable_list, paste(specific_level,'_', sep = ''))]
    }
    
    mean_value = data_measurement %>%
      group_by(ProjectID, HRDtype) %>%
      summarise(!!specific_level := mean(!!sym(variable[[1]]), na.rm = TRUE), .groups = 'drop')
    if (i == 1){
      mean_df = mean_value
    }else{
      mean_df <- merge(mean_df, mean_value, by = c("ProjectID", "HRDtype"), all.x = TRUE)
    }
    
  }
  
  
  
  subset_mean_df = mean_df[,3:length(mean_df)]
  
  z_scored = sapply(subset_mean_df, function(subset_mean_df) (subset_mean_df-mean(subset_mean_df))/sd(subset_mean_df))
  
  matrix = as.matrix(z_scored)
  
  rownames(matrix) = mean_df$ProjectID
  
  matrix_t = t(matrix)
  
  ha = HeatmapAnnotation(Cancertype = mean_df$ProjectID, HRDType = mean_df$HRDtype, col = list(Cancertype = color_vector,HRDType = c('High' = 'orange', 'Low' = 'blue')))
  
  if (level == 'Arm'){
    size_row = 7
    size_col = 7
  }else{
    size_row = 12
    size_col = 7
  }
  
  hm = Heatmap(matrix_t, cluster_rows = FALSE, cluster_columns = TRUE, top_annotation  = ha,
          show_column_names  = TRUE, raster_quality = 10,row_names_gp = gpar(fontsize = size_row),column_names_gp = gpar(fontsize = size_col), column_title=paste('Mean ',measure_full_name,', Level ',level, sep = ''),
          column_gap = unit(0, "mm"), name = 'Z-score')
  
  if (save){
    png(file=paste('../data/figures/heatmap_',measure_short_name,'_mean_',level, '.png', sep = ''))
  }
  
  draw(hm)
  
  if (save){
    dev.off()
  }
}


patterns = c("^chr[1-9][0-9]?_general_cin$", "^chr[1-9][0-9]?_cn_cin$", "^chr[1-9][0-9]?_number_cn_loh$",
             "^chr[1-9][0-9]?_number_cn_gain$", "^chr[1-9][0-9]?_number_cn_amp$", "^chr[1-9][0-9]?_number_homo_del$",
             "^chr[1-9][0-9]?_number_hemi_del$", "^chr[1-9][0-9]?_number_gain$", "^chr[1-9][0-9]?_number_amp$")

measures_long = c('General CIN', 'CN-CIN', 'Number of CN-LOH', 'Number of CN-Gain', 'Number of CN-Amp',
                  'Number of Homo-Del', 'Number of Hemi-Del', 'Number of Gain', 'Number of Amp')

measures_short = c('general_cin', 'cn_cin', 'n_cn_loh', 'n_cn_gain', 'n_cn_amp', 'n_homo_del',
                   'n_hemi_del', 'n_gain', 'n_amp')


for (i in 1:length(patterns)){
  pattern = patterns[i]
  variables = grep(pattern, names(data_test_type), value = TRUE)
  subset = data_17_cohorts[, c("ProjectID", "HRDtype", variables)]
  plotmeanHeatmapChromosome(subset, variables, named_palette, 'Chromosome', measures_long[i], measures_short[i], save = TRUE)
  plotHeatmapChromosome(subset, variables, named_palette, 'Chromosome', measures_long[i], measures_short[i], save = TRUE)
}

######## Single measure use this here
i = 1
pattern = patterns[i]
variables = grep(pattern, names(data_test_type), value = TRUE)
subset = data_17_cohorts[, c("ProjectID", "HRDtype", variables)]
#plotmeanHeatmapChromosome(subset, variables, named_palette, 'Chromosome', measures_long[i], measures_short[i], save = FALSE)
plotHeatmapChromosome(subset, variables, named_palette, 'Chromosome', measures_long[i], measures_short[i], save = FALSE)





########################### Events ############################ PLOTTING heatmaps

plotHeatmapChromosomeEvents = function(binary_data, variable_list, level, event_name_long, event_name_short, save =FALSE){
  
  chr_order <- paste0("chr", 1:22)
  chr_arm_order = chr_order
  if (level == 'Arm') {
    chr_order <- rep(chr_order, each = 2)
    chr_arm_order <- rep(chr_order, each = 2)
    chr_arm_order <- paste0(chr_order, c("_p", "_q"))
    chr_order <- paste0(chr_order, c(".p", ".q"))
    
  }
  
  mean_df <- data.frame(ProjectID = character(),
                        HRDtype = character(),
                        stringsAsFactors = FALSE)
  
  
  for (i in 1:length(chr_order)){
    
    specific_level = chr_order[[i]]
    if (level == 'Arm'){
      variable = variable_list[startsWith(variable_list, paste(chr_arm_order[[i]],'_', sep = ''))]
    }else{
      variable = variable_list[startsWith(variable_list, paste(specific_level,'_', sep = ''))]
    }
    
    mean_value = binary_data %>%
      group_by(ProjectID, HRDtype) %>%
      summarise(!!specific_level := mean(!!sym(variable[[1]]), na.rm = TRUE), .groups = 'drop')
    if (i == 1){
      mean_df = mean_value
    }else{
      mean_df <- merge(mean_df, mean_value, by = c("ProjectID", "HRDtype"), all.x = TRUE)
    }
    
  }
  
  subset_mean_df <- mean_df[,3:ncol(mean_df)]

  matrix = as.matrix(subset_mean_df)
  rownames(matrix) = mean_df$ProjectID
  
  matrix_t = t(matrix)
  
  max_value = max(matrix_t)
  
  col_fun = colorRamp2(c(0, max_value), c("antiquewhite1", "red"))
  
  ha = HeatmapAnnotation(Cancertype = mean_df$ProjectID, HRDType = mean_df$HRDtype, col = list(Cancertype = named_palette,HRDType = c('High' = 'orange', 'Low' = 'blue')))
  
  if (level == 'Arm'){
    size_row = 7
    size_col = 7
  }else{
    size_row = 12
    size_col = 7
  }
  
  hm = Heatmap(matrix_t, cluster_rows = FALSE, cluster_columns = TRUE, show_column_dend = FALSE,
          top_annotation  = ha, show_column_names  = TRUE, raster_quality = 10, column_title= paste('Mean of whole ', level,' event, ', event_name_long, sep = ''),
          column_gap = unit(0, "mm"), name = 'Mean appearance', row_names_gp = gpar(fontsize = size_row),column_names_gp = gpar(fontsize = size_col), col = col_fun)
  
  if (save){
    png(file=paste('../data/figures/heatmap_event_',event_name_short,'_mean_',level, '.png', sep = ''))
  }
  
  draw(hm)
  
  if (save){
    dev.off()
  }
}


pattern_status = "^chr[1-9][0-9]?_status$"
variable_list = grep(pattern_status, names(data_test_type), value = TRUE)
subset_chr_status = data_17_cohorts[, c("ProjectID", "HRDtype", variable_list)]
events = c('CN-LOH', 'CN-GAIN', 'CN-AMP', 'Homo-Del', 'Hemi-Del', 'Neutral', 'Gain', 'Amp')
events_shortname = c('cn_loh', 'cn_gain', 'cn_amp', 'homo_del', 'hemi_del', 'neutral', 'gain', 'amp')
columns_to_convert <- names(subset_chr_status)[3:ncol(subset_chr_status)]

for (i in 1:length(events)){
  df_binary = subset_chr_status
  df_binary[columns_to_convert] <- lapply(subset_chr_status[columns_to_convert], function(x) as.integer(x == events[i]))
  plotHeatmapChromosomeEvents(df_binary,variable_list, 'Chromosome', events[i], events_shortname[i], save =TRUE)
}



################################## ARMS ###########################
patterns <- c("^chr[1-9][0-9]?_(q|p)_general_cin$", "^chr[1-9][0-9]?_(q|p)_cn_cin$", "^chr[1-9][0-9]?_(q|p)_number_cn_loh$",
              "^chr[1-9][0-9]?_(q|p)_number_cn_gain$", "^chr[1-9][0-9]?_(q|p)_number_cn_amp$", "^chr[1-9][0-9]?_(q|p)_number_homo_del$",
              "^chr[1-9][0-9]?_(q|p)_number_hemi_del$", "^chr[1-9][0-9]?_(q|p)_number_gain$", "^chr[1-9][0-9]?_(q|p)_number_amp$")



measures_long = c('General CIN', 'CN-CIN', 'Number of CN-LOH', 'Number of CN-Gain', 'Number of CN-Amp',
                  'Number of Homo-Del', 'Number of Hemi-Del', 'Number of Gain', 'Number of Amp')

measures_short = c('general_cin', 'cn_cin', 'n_cn_loh', 'n_cn_gain', 'n_cn_amp', 'n_homo_del',
                   'n_hemi_del', 'n_gain', 'n_amp')


for (i in 1:length(patterns)){
  pattern = patterns[i]
  variables = grep(pattern, names(data_test_type), value = TRUE)
  subset = data_17_cohorts[, c("ProjectID", "HRDtype", variables)]
  plotmeanHeatmapChromosome(subset, variables, named_palette, 'Arm', measures_long[i], measures_short[i], save = TRUE)
  plotHeatmapChromosome(subset, variables, named_palette, 'Arm', measures_long[i], measures_short[i], save = TRUE)
}

i = 2
pattern <- "^chr[1-9][0-9]?_(q|p)_cn_cin$"
variables = grep(pattern, names(data_test_type), value = TRUE)
subset = data_17_cohorts[, c("ProjectID", "HRDtype", variables)]
plotmeanHeatmapChromosome(subset, variables, named_palette, 'Arm', measures_long[i], measures_short[i], save = FALSE)
#plotHeatmapChromosome(subset, variables, named_palette, 'Arm', measures_long[i], measures_short[i], save = FALSE)


pattern_status = "^chr[1-9][0-9]?_(q|p)_status$"
variable_list = grep(pattern_status, names(data_test_type), value = TRUE)
subset_chr_status = data_17_cohorts[, c("ProjectID", "HRDtype", variable_list)]
events = c('CN-LOH', 'CN-GAIN', 'CN-AMP', 'Homo-Del', 'Hemi-Del', 'Neutral', 'Gain', 'Amp')
events_shortname = c('cn_loh', 'cn_gain', 'cn_amp', 'homo_del', 'hemi_del', 'neutral', 'gain', 'amp')
columns_to_convert <- names(subset_chr_status)[3:ncol(subset_chr_status)]

for (i in 1:length(events)){
  df_binary = subset_chr_status
  df_binary[columns_to_convert] <- lapply(subset_chr_status[columns_to_convert], function(x) as.integer(x == events[i]))
  plotHeatmapChromosomeEvents(df_binary,variable_list, 'Arm', events[i], events_shortname[i], save =TRUE)
}




#############################
# Statistic analysis

patterns <- c("^chr[1-9][0-9]?_general_cin$", "^chr[1-9][0-9]?_cn_cin$", "^chr[1-9][0-9]?_number_cn_loh$",
              "^chr[1-9][0-9]?_number_cn_gain$", "^chr[1-9][0-9]?_number_cn_amp$", "^chr[1-9][0-9]?_number_homo_del$",
              "^chr[1-9][0-9]?_number_hemi_del$", "^chr[1-9][0-9]?_number_gain$", "^chr[1-9][0-9]?_number_amp$", "^chr[1-9][0-9]?_(q|p)_general_cin$", "^chr[1-9][0-9]?_(q|p)_cn_cin$", "^chr[1-9][0-9]?_(q|p)_number_cn_loh$",
              "^chr[1-9][0-9]?_(q|p)_number_cn_gain$", "^chr[1-9][0-9]?_(q|p)_number_cn_amp$", "^chr[1-9][0-9]?_(q|p)_number_homo_del$",
              "^chr[1-9][0-9]?_(q|p)_number_hemi_del$", "^chr[1-9][0-9]?_(q|p)_number_gain$", "^chr[1-9][0-9]?_(q|p)_number_amp$")


measures = c('general_cin', 'cn_cin', 'totalCIN', 'n_cn_loh', 'n_cn_gain', 'n_cn_amp', 'n_homo_del', 'n_hemi_del', 'n_gain', 'n_amp')



##### Question to answer: Is there any chromosome that shows significant higher CIN in many cancer types or are most of them singifincat high across all cancer types

pattern = "^chr[1-9][0-9]?_general_cin$"
features = grep(pattern, names(data_17_cohorts), value = TRUE)
chromosome_data = data_17_cohorts[, c("ProjectID", "HRD_sum",'HRDtype', features )]
empty_df <- data.frame(matrix(nrow = length(features), ncol = 1))
colnames(empty_df) = 'p-value'
rownames(empty_df) = features

i = 1
for (feature in features){
  formula <- as.formula(paste(feature, "~ HRDtype"))
  m1 = wilcox.test(formula, data = chromosome_data,exact = FALSE, paired = FALSE)
  print(feature)
  print(m1)
  empty_df[i,'p-value'] = m1$p.value
  i = i + 1
}
write.csv(empty_df,'../data/general_cin_chromomsome_hrdh_l_p_value.csv', row.names = TRUE)
### Conclusion -> every chromosome shows signfificant higher general CIN in HRD-high over all cancertypes

#### Same question for each cancer type
pattern = "^chr[1-9][0-9]?_general_cin$"
features = grep(pattern, names(data_17_cohorts), value = TRUE)
ids = unique(data_17_cohorts$ProjectID)
chromosome_data = data_17_cohorts[, c("ProjectID", "HRD_sum",'HRDtype', features )]
empty_df <- data.frame(matrix(nrow = length(features)*length(ids), ncol = 3))
colnames(empty_df) = c('Feature','p-value', 'Cancertype')
non_significant <- data.frame(Cancertype = character(0), Feature = character(0), p_value = numeric(0), stringsAsFactors = FALSE)
colnames(non_significant) = c('Cancertype', 'Feature', 'p-value')
i = 1
num_non_significant = 0
for (cancertype in ids){
  subset = chromosome_data[chromosome_data$ProjectID == cancertype,]
  for (feature in features){
    formula <- as.formula(paste(feature, "~ HRDtype"))
    m1 = wilcox.test(formula, data = subset, exact = FALSE, paired = FALSE)
    empty_df[i,'p-value'] = m1$p.value
    empty_df[i,'Feature'] = feature
    empty_df[i,'Cancertype'] = cancertype
    i = i + 1
    if (m1$p.value >= 0.05){
      new_row = c(cancertype, feature, m1$p.value)
      non_significant[num_non_significant+1,'Cancertype'] = cancertype
      non_significant[num_non_significant+1,'Feature'] = feature
      non_significant[num_non_significant+1,'p-value'] = m1$p.value
      num_non_significant = num_non_significant + 1
    }
  }
}

## Conclusion: OS shows no significant difference in 12 chromosome, ACC in three, KICH in 9, KIRC in 1, KIRP in 5, LGG in 1, LUSC in 1 and OV in 6
write.csv(empty_df,'../data/general_cin_chromomsome_hrdh_l_p_value_by_cancertype.csv', row.names = FALSE)
write.csv(non_significant,'../data/general_cin_chromomsome_hrdh_l_p_value_by_cancertype_non_significant.csv', row.names = FALSE)




pattern = "^chr[1-9][0-9]?_(q|p)_general_cin$"
features = grep(pattern, names(data_17_cohorts), value = TRUE)
features <- features[features != 'chr13_p_general_cin']
features <- features[features != 'chr14_p_general_cin']
features <- features[features != 'chr15_p_general_cin']
features <- features[features != 'chr21_p_general_cin']
features <- features[features != 'chr22_p_general_cin']
chromosome_data = data_17_cohorts[, c("ProjectID", "HRD_sum",'HRDtype', features )]
empty_df <- data.frame(matrix(nrow = length(features), ncol = 1))
colnames(empty_df) = 'p-value'
rownames(empty_df) = features

i = 1
for (feature in features){
  formula <- as.formula(paste(feature, "~ HRDtype"))
  m1 = wilcox.test(formula, data = chromosome_data,exact = FALSE, paired = FALSE)
  print(feature)
  print(m1)
  empty_df[i,'p-value'] = m1$p.value
  i = i + 1
}
write.csv(empty_df,'../data/general_cin_arm_hrdh_l_p_value.csv', row.names = TRUE)


#### General CIN arms, do most of the arms show higher general CIN in the HRD-high than HRD-low in most cancertypes?
pattern = "^chr[1-9][0-9]?_(q|p)_general_cin$"
features = grep(pattern, names(data_17_cohorts), value = TRUE)
features <- features[features != 'chr13_p_general_cin']
features <- features[features != 'chr14_p_general_cin']
features <- features[features != 'chr15_p_general_cin']
features <- features[features != 'chr21_p_general_cin']
features <- features[features != 'chr22_p_general_cin']
ids = unique(data_17_cohorts$ProjectID)
chromosome_data = data_17_cohorts[, c("ProjectID", "HRD_sum",'HRDtype', features )]
my_list <- my_list[my_list != value_to_remove]
empty_df <- data.frame(matrix(nrow = length(features)*length(ids), ncol = 3))
colnames(empty_df) = c('Feature','p-value', 'Cancertype')
non_significant <- data.frame(Cancertype = character(0), Feature = character(0), p_value = numeric(0), stringsAsFactors = FALSE)
colnames(non_significant) = c('Cancertype', 'Feature', 'p-value')
i = 1
num_non_significant = 0
for (cancertype in ids){
  subset = chromosome_data[chromosome_data$ProjectID == cancertype,]
  for (feature in features){
    formula <- as.formula(paste(feature, "~ HRDtype"))
    m1 = wilcox.test(formula, data = subset, exact = FALSE, paired = FALSE)
    empty_df[i,'p-value'] = m1$p.value
    empty_df[i,'Feature'] = feature
    empty_df[i,'Cancertype'] = cancertype
    i = i + 1
    if (m1$p.value >= 0.05){
      new_row = c(cancertype, feature, m1$p.value)
      non_significant[num_non_significant+1,'Cancertype'] = cancertype
      non_significant[num_non_significant+1,'Feature'] = feature
      non_significant[num_non_significant+1,'p-value'] = m1$p.value
      num_non_significant = num_non_significant + 1
    }
  }
}

write.csv(empty_df,'../data/general_cin_arm_hrdh_l_p_value_by_cancertype.csv', row.names = FALSE)
write.csv(non_significant,'../data/general_cin_arm_hrdh_l_p_value_by_cancertype_non_significant.csv', row.names = FALSE)


##### General CIn cytoband
pattern = "^chr[1-9][0-9]?_.*_general_cin$"
features = grep(pattern, names(cytoband_18), value = TRUE)
chromosome_data = cytoband_18[, c("ProjectID", "HRD_sum",'HRDtype', features )]
empty_df <- data.frame(matrix(nrow = length(features), ncol = 1))
colnames(empty_df) = 'p-value'
rownames(empty_df) = features

i = 1
for (feature in features){
  formula <- as.formula(paste(feature, "~ HRDtype"))
  m1 = wilcox.test(formula, data = chromosome_data,exact = FALSE, paired = FALSE)
  print(feature)
  print(m1)
  empty_df[i,'p-value'] = m1$p.value
  i = i + 1
}
write.csv(empty_df,'../data/general_cin_cytoband_hrdh_l_p_value.csv', row.names = TRUE)


pattern = "^chr[1-9][0-9]?_.*_cn_cin$"
features = grep(pattern, names(cytoband_18), value = TRUE)
ids = unique(cytoband_18$ProjectID)
chromosome_data = cytoband_18[, c("ProjectID", "HRD_sum",'HRDtype', features )]
empty_df <- data.frame(matrix(nrow = length(features)*length(ids), ncol = 3))
colnames(empty_df) = c('Feature','p-value', 'Cancertype')
non_significant <- data.frame(Cancertype = character(0), Feature = character(0), p_value = numeric(0), stringsAsFactors = FALSE)
colnames(non_significant) = c('Cancertype', 'Feature', 'p-value')
i = 1
num_non_significant = 0
for (cancertype in ids){
  subset = chromosome_data[chromosome_data$ProjectID == cancertype,]
  for (feature in features){
    formula <- as.formula(paste(feature, "~ HRDtype"))
    m1 = wilcox.test(formula, data = subset, exact = FALSE, paired = FALSE)
    empty_df[i,'p-value'] = m1$p.value
    empty_df[i,'Feature'] = feature
    empty_df[i,'Cancertype'] = cancertype
    i = i + 1
    if (m1$p.value >= 0.05){
      new_row = c(cancertype, feature, m1$p.value)
      non_significant[num_non_significant+1,'Cancertype'] = cancertype
      non_significant[num_non_significant+1,'Feature'] = feature
      non_significant[num_non_significant+1,'p-value'] = m1$p.value
      num_non_significant = num_non_significant + 1
    }
  }
}

write.csv(empty_df,'../data/general_cin_cytoband_hrdh_l_p_value_by_cancertype.csv', row.names = FALSE)
write.csv(non_significant,'../data/general_cin_cytoband_hrdh_l_p_value_by_cancertype_non_significant.csv', row.names = FALSE)


pattern = "^chr[1-9][0-9]?_general_cin$"
features = grep(pattern, names(data_17_cohorts), value = TRUE)
feauter_data = data_17_cohorts[, c("ProjectID", "HRD_sum",'HRDtype', features )]
ids = unique(data_17_cohorts$ProjectID)
ids = c('All',ids)
#empty_df <- data.frame(Cancertype = character(0), Feature = character(0), p_value = numeric(0), stringsAsFactors = FALSE)
empty_df = data.frame(matrix(ncol=4,nrow=0, dimnames=list(NULL, c("Cancertype", "Feature", "p_value", 'significant'))))

i = 1
for (cancertype in ids){
  if (cancertype == 'All'){
    subdata = data.frame(feauter_data)
  }else{
    subdata = feauter_data[feauter_data$ProjectID == cancertype,]
  }
  for (feature in features){
    formula <- as.formula(paste(feature, "~ HRDtype"))
    m1 = wilcox.test(formula, data = subdata,exact = FALSE, paired = FALSE)
    empty_df[i,'Cancertype'] = cancertype
    empty_df[i,'Feature'] = feature
    empty_df[i,'p_value'] = m1$p.value
    if (m1$p.value >= 0.05){
      empty_df[i,'significant'] = 'No'
    }else{
      empty_df[i,'significant'] = 'Yes'
    }
    i = i + 1
  }
}




analysisFeatureAllLevels = function(feature, data, data_cytoband){
  
  pattern_chr = paste('^chr[1-9][0-9]?_',feature,'$',sep = '')
  pattern_arm = paste('^chr[1-9][0-9]?_(q|p)_',feature,'$',sep = '')
  pattern_cyto = paste('^chr[1-9][0-9]?_.*_',feature,'$',sep = '')
  
  analysisFeature(pattern_chr, data, 'chromosome', feature)
  analysisFeature(pattern_arm, data, 'arm', feature)
  analysisFeature(pattern_cyto, data_cytoband, 'cytoband', feature)
}

analysisFeature = function(pattern, data, level, feature_name){
  
  features = grep(pattern, names(data), value = TRUE)
  if (level == 'arm'){
    features <- features[features != paste('chr13_p_',feature_name,sep = '')]
    features <- features[features != paste('chr14_p_',feature_name,sep = '')]
    features <- features[features != paste('chr15_p_',feature_name,sep = '')]
    features <- features[features != paste('chr21_p_',feature_name,sep = '')]
    features <- features[features != paste('chr22_p_',feature_name,sep = '')]
  }
  feature_data = data[, c("ProjectID", "HRD_sum",'HRDtype', features )]
  empty_df = data.frame(matrix(ncol=4,nrow=0, dimnames=list(NULL, c("Cancertype", "Feature", "p_value", 'significant'))))
  ids = unique(data$ProjectID)
  ids = c('All',ids)
  i = 1
  for (cancertype in ids){
    if (cancertype == 'All'){
      subdata = data.frame(feature_data)
    }else{
      subdata = feature_data[feature_data$ProjectID == cancertype,]
    }
    for (feature in features){
      formula <- as.formula(paste(feature, "~ HRDtype"))
      m1 = wilcox.test(formula, data = subdata,exact = FALSE, paired = FALSE)
      empty_df[i,'Cancertype'] = cancertype
      empty_df[i,'Feature'] = feature
      empty_df[i,'p_value'] = m1$p.value
      if (m1$p.value >= 0.05 | is.nan(m1$p.value)){
        empty_df[i,'significant'] = 'No'
      }else{
        empty_df[i,'significant'] = 'Yes'
      }
      i = i + 1
    }
  }
  write.csv(empty_df,paste('../data/',feature_name,'_',level,'_wilcoxon_pvalues.csv', sep = ''), row.names = FALSE)
}


analysisFeatureAllLevels('general_cin', data_17_cohorts, cytoband_18)
analysisFeatureAllLevels('cn_cin', data_17_cohorts, cytoband_18)
analysisFeatureAllLevels('number_cn_loh', data_17_cohorts, cytoband_18)
analysisFeatureAllLevels('number_cn_gain', data_17_cohorts, cytoband_18)
analysisFeatureAllLevels('number_cn_amp', data_17_cohorts, cytoband_18)
analysisFeatureAllLevels('number_homo_del', data_17_cohorts, cytoband_18)
analysisFeatureAllLevels('number_hemi_del', data_17_cohorts, cytoband_18)
analysisFeatureAllLevels('number_gain', data_17_cohorts, cytoband_18)
analysisFeatureAllLevels('number_amp', data_17_cohorts, cytoband_18)

#### Was Nan because all is 0
pattern_cyto = paste('^chr[1-9][0-9]?_.*_','cn_cin','$',sep = '')
analysisFeature(pattern_cyto, cytoband_18, 'cytoband', 'cn_cin')
cytoband_paad = cytoband_18[cytoband_18$ProjectID == 'TCGA-PAAD',]
cytoband_chr20 = cytoband_paad[, c('HRDtype', 'chr20_q11.21_cn_cin')]

analysisResultsAllLevel = function(feature, cyto_annotation = NULL){
  analysisResults(feature, 'chromosome',cyto_annotation)
  analysisResults(feature, 'arm',cyto_annotation)
  analysisResults(feature, 'cytoband', cyto_annotation)
}

analysisResults = function(feature, level, cyto_annotation = NULL){
  file_path = paste('../data/analysis_',feature,'_',level,'.txt',sep = '')
  file_conn = file(file_path, "w")
  data = read.csv(paste('../data/',feature,'_',level,'_wilcoxon_pvalues.csv', sep = ''), sep = ',')
  cat('Wilcoxon rank sum test was used for calculating the p-value, wilcox.test(formula, data = subdata, exact = FALSE, paired = FALSE)', '\n', file = file_conn)
  cat('Significant are p-values < 0.05', '\n\n', file = file_conn)
  subset_all = data[data$Cancertype == 'All',]
  
  if(all(subset_all$significant == 'Yes')){
    cat(paste('All ',feature,' on ',level, ' level are significant higher in HRD-high (independent from cancertype)', sep = ''), "\n", file = file_conn)
  }else{
    cat(paste('Not all ',feature,' on ',level, ' level are significant higher in HRD-high (independent from cancertype)', sep = ''), "\n", file = file_conn)
    no_subset = subset_all[subset_all$significant == 'No',]
    cat('Not significant are: ', file = file_conn)
    for(no_feature in no_subset$Feature){
      cat(paste(no_feature,', ',sep = ''), file = file_conn)
    }
    cat('\n', file = file_conn)
  }
  
  cat('\n', file = file_conn)
  
  subset_byCancertype = data[data$Cancertype != 'All',]
  
  if (all(subset_byCancertype$significant == 'Yes')){
    cat(paste('For every cancer type the feature ',feature,' on ',level, ' level are significant higher in HRD-high', sep = ''), "\n", file = file_conn)
  }else{
    cat(paste('In these cancer types every ',level,' shows significant higher ' ,feature,' in HRD-high: ', sep = ''), file = file_conn)
    for (cancertype in unique(subset_byCancertype$Cancertype)){
      cancertype_subset = subset_byCancertype[subset_byCancertype$Cancertype == cancertype,]
      if(all(cancertype_subset$significant == 'Yes')){
        cat(paste(cancertype,', ',sep = ''), file = file_conn)
      }
    }
    cat('\n', file = file_conn)
    subset_not_significant = subset_byCancertype[subset_byCancertype$significant == 'No',]
    cat('\n', file = file_conn)
    cat(paste('The following cancer types showed at least one ',level,' with not significant difference:', sep=''), '\n', file = file_conn)
    cat('\n', file = file_conn)
    
    
    
    for(cancertype in unique(subset_not_significant$Cancertype)){
      cancertype_subset = subset_byCancertype[subset_byCancertype$Cancertype == cancertype,]
      
      not_signi = cancertype_subset[cancertype_subset$significant == 'No',]
      signi = cancertype_subset[cancertype_subset$significant == 'Yes',]
      
      cat(paste('-',cancertype, sep = ''),'\n', file = file_conn)
      cat(paste('     Significant (',nrow(signi),'/',nrow(cancertype_subset),'): ', sep = ''), file = file_conn)
      
      for (feature_level in signi$Feature){
        p_value = signi[signi$Feature == feature_level,'p_value']
        if (level == 'cytoband'){
          info = gsub(paste0("_", feature), "", feature_level)
          gene = cyto_annotation[cyto_annotation$cytoband_chr == info, 'symbol']
          if (length(gene) > 1){
            gene = paste(gene, collapse = ", ")
          }
          info = paste(info, '(', gene,')', sep = '')
        }else{
          info = gsub(paste0("_", feature), "", feature_level)
        }
        if (round(p_value,4) == 0){
          cat(paste(info,' (<0.0001), ', sep = ''), file = file_conn)
        }else{
          cat(paste(info,' (',round(p_value,4),'), ', sep = ''), file = file_conn)
        }
        
      }
      cat('\n', file = file_conn)
      cat(paste('     Not Significant (',nrow(not_signi),'/',nrow(cancertype_subset),'): ', sep = ''), file = file_conn)
      
      for (feature_level in not_signi$Feature){
        if (level == 'cytoband'){
          info = gsub(paste0("_", feature), "", feature_level)
          gene = cyto_annotation[cyto_annotation$cytoband_chr == info, 'symbol']
          if (length(gene) > 1){
            gene = paste(gene, collapse = ", ")
          }
          info = paste(info, '(', gene,')', sep = '')
        }else{
          info = gsub(paste0("_", feature), "", feature_level)
        }
        p_value = not_signi[not_signi$Feature == feature_level,'p_value']
        cat(paste(info,' (',round(p_value,4),'), ', sep = ''), file = file_conn)
      }
      cat('\n', file = file_conn)
      cat('\n', file = file_conn)
    }
    cat(paste('Here a list of the appearances of ',level,' that are not signifincat',sep = ''),'\n' , file = file_conn)
    for (feature_level in unique(subset_not_significant$Feature)){
      subset_feauture = subset_not_significant[subset_not_significant$Feature == feature_level,]
      if (level == 'cytoband'){
        info = gsub(paste0("_", feature), "", feature_level)
        gene = cyto_annotation[cyto_annotation$cytoband_chr == info, 'symbol']
        if (length(gene) > 1){
          gene = paste(gene, collapse = ", ")
        }
        info = paste(info,' ',nrow(subset_feauture),'/17 : ',sep = '')
        info = paste(info, ' (', gene,')', sep = '')
      }else{
        info = gsub(paste0("_", feature), "", feature_level)
        if (level == 'chromosome'){
          subset_genes = cyto_annotation[grepl(paste(info,'_',sep = ''), cyto_annotation$cytoband_chr), ]
        }else{
          subset_genes = cyto_annotation[grepl(info, cyto_annotation$cytoband_chr), ]
        }
        
        info = paste(info,' ',nrow(subset_feauture),'/18 : ',sep = '')
        genes = subset_genes$symbol
        genes = paste(genes, collapse = ", ")
        info = paste(info,' (', genes,')', sep = '')
      }
      info = paste(info,' ; Cancertypes: (', sep = '')
      for (cancertype in subset_feauture$Cancertype){
        info = paste(info,cancertype,', ', sep = '')
      }
      info = paste(info,')', sep = '')
      cat(info,'\n', file = file_conn)
    }
  }
  
  close(file_conn)
}

cytoband_HRR_annotation =  read.csv('../data/cytoband_HRRgenes_annotated.csv',row.names = 1)

test = cytoband_HRR_annotation[cytoband_HRR_annotation$cytoband_chr == 'chr15_q26.1', 'symbol']

test = cytoband_HRR_annotation[grepl('chr15', cytoband_HRR_annotation$cytoband_chr), ]
test_2 = test$symbol
print(paste(test_2, collapse = ", "))
combined_string <- paste(test, collapse = ", ")
analysisResults('general_cin', 'chromosome',cytoband_HRR_annotation)

analysisResultsAllLevel('general_cin', cytoband_HRR_annotation)
analysisResultsAllLevel('cn_cin', cytoband_HRR_annotation)
analysisResultsAllLevel('number_cn_loh', cytoband_HRR_annotation)
analysisResultsAllLevel('number_cn_gain', cytoband_HRR_annotation)
analysisResultsAllLevel('number_cn_amp', cytoband_HRR_annotation)
analysisResultsAllLevel('number_homo_del', cytoband_HRR_annotation)
analysisResultsAllLevel('number_hemi_del', cytoband_HRR_annotation)
analysisResultsAllLevel('number_gain', cytoband_HRR_annotation)
analysisResultsAllLevel('number_amp', cytoband_HRR_annotation)


feature = 'cn_cin'
level = 'arm'
data = read.csv(paste('../data/',feature,'_',level,'_wilcoxon_pvalues.csv', sep = ''), sep = ',')
subset_byCancertype = data[data$Cancertype != 'All',]

for (feature_level in unique(subset_byCancertype$Feature)){
  feature_subset = subset_byCancertype[subset_byCancertype$Feature == feature_level,]
  if(all(feature_subset$significant == 'Yes')){
    print(feature_level)
  }
}

################## General CIN percentages ############################
chromosomes = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22)
data = read.csv('../data/general_cin_chromosome_wilcoxon_pvalues.csv', sep = ',')
subset_byCancertype = data[data$Cancertype != 'All',]
for (chr in chromosomes){
  
  general_cin = paste('chr',chr,'_general_cin',sep = '')
  number_homo_del = paste('chr',chr,'_number_homo_del',sep = '')
  number_hemi_del = paste('chr',chr,'_number_hemi_del',sep = '')
  number_gain = paste('chr',chr,'_number_gain',sep = '')
  number_amp = paste('chr',chr,'_number_amp',sep = '')
  
  subdata = data_17_cohorts[, c('ProjectID', 'HRDtype', general_cin, number_homo_del, number_hemi_del,
                                number_gain, number_amp)]
  result <- subdata %>%
    group_by(ProjectID, HRDtype) %>%
    summarize(
      mean_general_cin = mean(!!sym(general_cin)),
      mean_number_homo_del = mean(!!sym(number_homo_del)),
      mean_number_hemi_del = mean(!!sym(number_hemi_del)),
      mean_number_gain = mean(!!sym(number_gain)),
      mean_number_amp = mean(!!sym(number_amp))
    )%>%
    # Calculate the percentage for each score relative to the sum of Score1
    mutate(
      Percentage_number_homo_del = round((mean_number_homo_del / mean_general_cin) * 100, 2),
      Percentage_number_hemi_del = round((mean_number_hemi_del / mean_general_cin) * 100, 2),
      Percentage_number_gain = round((mean_number_gain / mean_general_cin) * 100, 2),
      Percentage_number_amp = round((mean_number_amp / mean_general_cin) * 100, 2)
    )
  
  write.csv(result,paste('../data/general_cin_specific_chromosomes_precentages/general_cin_prec_chr',chr,'.csv', sep = ''), row.names = FALSE)

  df_long <- result %>%
    pivot_longer(cols = starts_with("Percentage"), names_to = "NumericColumn", values_to = "Value")
  
  df_long$significants = 'Not significant'
  
  significant_cancers <- unique(subset_byCancertype$Cancertype[subset_byCancertype$significant == 'Yes' & subset_byCancertype$Feature == general_cin])
  
  df_long[df_long$ProjectID %in% significant_cancers, 'significants'] = 'significant'
  
  df_long$Cancertype <- sapply(strsplit(df_long$ProjectID, "-"), function(x) x[2])
  df_long <- df_long[order(df_long$significants == 'Not significant', decreasing = TRUE), ]
  df_long$NumericColumn <- gsub("Percentage_number_homo_del", "Number of homo-del", df_long$NumericColumn)
  df_long$NumericColumn <- gsub("Percentage_number_hemi_del", "Number of hemi-del", df_long$NumericColumn)
  df_long$NumericColumn <- gsub("Percentage_number_gain", "Number of gain", df_long$NumericColumn)
  df_long$NumericColumn <- gsub("Percentage_number_amp", "Number of amp", df_long$NumericColumn)
  plot = ggplot(df_long) +
    geom_bar(aes(x = HRDtype, y = Value, fill = NumericColumn),
             position = "stack",
             stat = "identity") +
    facet_grid(~factor(Cancertype, levels=unique(df_long$Cancertype)),  switch = "x") +
    labs(title = paste("Precentage of events between HRD-high and HRD-low in Chromosome ",chr, sep =''),
         x = "Cancer Type", y = "Precentages",
         fill = "Events")+
    theme(strip.placement = "outside",
          strip.background = element_rect(fill = NA, color = "white"),
          strip.text.x = element_text(angle = 90, hjust = 1),  # Rotate facet labels
          axis.text.x = element_text(angle = 90, hjust = 1) )
    
  ggsave(paste('../data/figures/generalCIN_events_precentage_chromosome_',chr,'.png',sep = ''),plot,width = 8.63, height = 5.71)
  
}


chromosomes = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22)
#chromosomes = c(17)
data = read.csv('../data/cn_cin_chromosome_wilcoxon_pvalues.csv', sep = ',')
subset_byCancertype = data[data$Cancertype != 'All',]
for (chr in chromosomes){
  
  cn_cin = paste('chr',chr,'_cn_cin',sep = '')
  number_cn_loh = paste('chr',chr,'_number_cn_loh',sep = '')
  number_cn_gain = paste('chr',chr,'_number_cn_gain',sep = '')
  number_cn_amp = paste('chr',chr,'_number_cn_amp',sep = '')
  
  subdata = data_17_cohorts[, c('ProjectID', 'HRDtype', cn_cin, number_cn_loh, number_cn_gain,
                                number_cn_amp)]
  result <- subdata %>%
    group_by(ProjectID, HRDtype) %>%
    summarize(
      mean_cn_cin = mean(!!sym(cn_cin)),
      mean_number_cn_loh = mean(!!sym(number_cn_loh)),
      mean_number_cn_gain = mean(!!sym(number_cn_gain)),
      mean_number_cn_amp = mean(!!sym(number_cn_amp))
    )%>%
    # Calculate the percentage for each score relative to the sum of Score1
    mutate(
      Percentage_number_cn_loh = round((mean_number_cn_loh / mean_cn_cin) * 100, 2),
      Percentage_number_cn_gain = round((mean_number_cn_gain / mean_cn_cin) * 100, 2),
      Percentage_number_cn_amp = round((mean_number_cn_amp / mean_cn_cin) * 100, 2)
    )
  
  write.csv(result,paste('../data/cn_cin_specific_chromosomes_precentages/cn_cin_prec_chr',chr,'.csv', sep = ''), row.names = FALSE)
  
  df_long <- result %>%
    pivot_longer(cols = starts_with("Percentage"), names_to = "NumericColumn", values_to = "Value")
  
  df_long$significants = 'Not significant'
  
  significant_cancers <- unique(subset_byCancertype$Cancertype[subset_byCancertype$significant == 'Yes' & subset_byCancertype$Feature == cn_cin])
  
  df_long[df_long$ProjectID %in% significant_cancers, 'significants'] = 'significant'
  
  df_long$Cancertype <- sapply(strsplit(df_long$ProjectID, "-"), function(x) x[2])
  df_long <- df_long[order(df_long$significants == 'Not significant', decreasing = TRUE), ]
  df_long$NumericColumn <- gsub("Percentage_number_cn_loh", "Number of CN-LOH", df_long$NumericColumn)
  df_long$NumericColumn <- gsub("Percentage_number_cn_gain", "Number of CN-Gain", df_long$NumericColumn)
  df_long$NumericColumn <- gsub("Percentage_number_cn_amp", "Number of CN-Amp", df_long$NumericColumn)
  
  plot = ggplot(df_long) +
    geom_bar(aes(x = HRDtype, y = Value, fill = NumericColumn),
             position = "stack",
             stat = "identity") +
    facet_grid(~factor(Cancertype, levels=unique(df_long$Cancertype)),  switch = "x") +
    labs(title = paste("Precentage of events between HRD-high and HRD-low in Chromosome ",chr, sep =''),
         x = "Cancer Type", y = "Precentages",
         fill = "Events")+
    theme(strip.placement = "outside",
          strip.background = element_rect(fill = NA, color = "white"),
          strip.text.x = element_text(angle = 90, hjust = 1),  # Rotate facet labels
          axis.text.x = element_text(angle = 90, hjust = 1) )
  
  ggsave(paste('../data/figures/CN_CIN_events_precentage_chromosome_',chr,'.png',sep = ''),plot,width = 8.63, height = 5.71)
  
}



plot = ggplot(result, aes(x = ProjectID, y = Percentage_number_homo_del + Percentage_number_hemi_del + Percentage_number_gain + Percentage_number_amp, fill = HRDtype)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Stacked Barplot by ProjectID and HRDtype",
       x = "ProjectID", y = "Stacked Value",
       fill = "HRDtype") +
  theme_minimal()
print(plot)

break


data_17_cohorts$new_general_cin = data_17_cohorts$chr11_general_cin + data_17_cohorts$chr4_general_cin + data_17_cohorts$chr5_general_cin + data_17_cohorts$chr14_general_cin + data_17_cohorts$chr22_general_cin

ggplot(data_17_cohorts, aes(x = ProjectID, y = chr4_general_cin, fill = HRDtype)) +
  geom_boxplot() +
  labs(x = "Project ID", y = "General CIN HRR", fill = "HRD-Type")  +
  ggtitle("General CIN HRR Cytobands, K-means") +
  theme_bw()+
  #ylim(0, 500)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  stat_compare_means(aes(group = HRDtype), method = 'wilcox.test', label ="p.signif") +
  scale_fill_manual(values = c("High" = "#1F77B4", "Low" = "#E377C2"), 
                    labels = c("HRD-high", "HRD-low"))

###################################
###################################
# Cytoband analysis ###### Figures for general cin only using cytoband containing HRR genes + heatmaps

general_cin_17_type_HRR = ggplot(cytoband_18, aes(x = ProjectID, y = general_cin_HRR, fill = HRDtype)) +
  geom_boxplot() +
  labs(x = "Project ID", y = "General CIN HRR", fill = "HRD-Type")  +
  ggtitle("General CIN HRR Cytobands, K-means") +
  theme_bw()+
  #ylim(0, 500)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  stat_compare_means(aes(group = HRDtype), method = 't.test', label ="p.signif") +
  scale_fill_manual(values = c("High" = "#1F77B4", "Low" = "#E377C2"), 
                    labels = c("HRD-high", "HRD-low"))

print(general_cin_17_type_HRR)
ggsave('../data/figures/boxplot_general_cin_17_type_HRR.png',general_cin_17_type_HRR,width = 8.63, height = 5.71)

cn_cin_17_type_HRR = ggplot(cytoband_18, aes(x = ProjectID, y = cn_cin_HRR, fill = HRDtype)) +
  geom_boxplot() +
  labs(x = "Project ID", y = "CN-CIN HRR", fill = "HRD-Type")  +
  ggtitle("CN-CIN HRR Cytobands, K-means") +
  theme_bw()+
  #ylim(0, 500)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  stat_compare_means(aes(group = HRDtype), method = 't.test', label ="p.signif") +
  scale_fill_manual(values = c("High" = "#1F77B4", "Low" = "#E377C2"), 
                    labels = c("HRD-high", "HRD-low"))

print(cn_cin_17_type_HRR)
ggsave('../data/figures/boxplot_cn_cin_17_type_HRR.png',cn_cin_17_type_HRR,width = 8.63, height = 5.71)

total_cin_17_type_HRR = ggplot(cytoband_18, aes(x = ProjectID, y = totalCIN_HRR, fill = HRDtype)) +
  geom_boxplot() +
  labs(x = "Project ID", y = "Total CIN HRR", fill = "HRD-Type")  +
  ggtitle("Total CIN HRR Cytobands, K-means") +
  theme_bw()+
  #ylim(0, 500)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  stat_compare_means(aes(group = HRDtype), method = 't.test', label ="p.signif") +
  scale_fill_manual(values = c("High" = "#1F77B4", "Low" = "#E377C2"), 
                    labels = c("HRD-high", "HRD-low"))

print(total_cin_17_type_HRR)
ggsave('../data/figures/boxplot_total_cin_17_type_HRR.png',total_cin_17_type_HRR,width = 8.63, height = 5.71)



plotmeanHeatmapCytoband = function(data_measurement, variable_list, color_vector, measure_full_name, measure_short_name, save = FALSE){
  
  
  mean_df <- data.frame(ProjectID = character(),
                        HRDtype = character(),
                        stringsAsFactors = FALSE)
  
  
  for (i in 1:length(variable_list)){
    
    variable = variable_list[[i]]
    
    mean_value = data_measurement %>%
      group_by(ProjectID, HRDtype) %>%
      summarise(!!variable := mean(!!sym(variable), na.rm = TRUE), .groups = 'drop')
    if (i == 1){
      mean_df = mean_value
    }else{
      mean_df <- merge(mean_df, mean_value, by = c("ProjectID", "HRDtype"), all.x = TRUE)
    }
    
  }

  subset_mean_df = mean_df[,3:length(mean_df)]
  
  z_scored = sapply(subset_mean_df, function(subset_mean_df) (subset_mean_df-mean(subset_mean_df))/sd(subset_mean_df))
  
  matrix = as.matrix(z_scored)
  
  rownames(matrix) = mean_df$ProjectID
  
  matrix_t = t(matrix)
  
  ha = HeatmapAnnotation(Cancertype = mean_df$ProjectID, HRDType = mean_df$HRDtype, col = list(Cancertype = color_vector,HRDType = c('High' = 'orange', 'Low' = 'blue')))
  
  size_row = 7
  size_col = 7
  
  hm = Heatmap(matrix_t, cluster_rows = FALSE, cluster_columns = TRUE, top_annotation  = ha,
               show_column_names  = TRUE, raster_quality = 10,row_names_gp = gpar(fontsize = size_row),column_names_gp = gpar(fontsize = size_col), column_title=paste('Mean ',measure_full_name,', Level Cytoband', sep = ''),
               column_gap = unit(0, "mm"), name = 'Z-score')
  
  if (save){
    png(file=paste('../data/figures/heatmap_',measure_short_name,'_mean_cytoband.png', sep = ''))
  }
  
  draw(hm)
  
  if (save){
    dev.off()
  }
}

patterns <- c("^chr[1-9][0-9]?_.*_general_cin$", "^chr[1-9][0-9]?_.*_cn_cin$", "^chr[1-9][0-9]?_.*_number_cn_loh$",
              "^chr[1-9][0-9]?_.*_number_cn_gain$", "^chr[1-9][0-9]?_.*_number_cn_amp$", "^chr[1-9][0-9]?_.*_number_homo_del$",
              "^chr[1-9][0-9]?_.*_number_hemi_del$", "^chr[1-9][0-9]?_.*_number_gain$", "^chr[1-9][0-9]?_.*_number_amp$")



measures_long = c('General CIN HRR', 'CN-CIN HRR', 'Number of CN-LOH HRR', 'Number of CN-Gain HRR', 'Number of CN-Amp HRR',
                  'Number of Homo-Del HRR', 'Number of Hemi-Del HRR', 'Number of Gain HRR', 'Number of Amp HRR')

measures_short = c('general_cin_HRR', 'cn_cin_HRR', 'n_cn_loh_HRR', 'n_cn_gain_HRR', 'n_cn_amp_HRR', 'n_homo_del_HRR',
                   'n_hemi_del_HRR', 'n_gain_HRR', 'n_amp_HRR')

for (i in 1:length(patterns)){
  pattern = patterns[i]
  variables = grep(pattern, names(cytoband_18), value = TRUE)
  subset = cytoband_18[, c("ProjectID", "HRDtype", variables)]
  plotmeanHeatmapCytoband(subset, variables, named_palette, measures_long[i], measures_short[i], save = TRUE)
}

i = 1
pattern = patterns[i]
variables = grep(pattern, names(cytoband_18), value = TRUE)
subset = cytoband_18[, c("ProjectID", "HRDtype", variables)]
#plotmeanHeatmapChromosome(subset, variables, named_palette, 'Chromosome', measures_long[i], measures_short[i], save = FALSE)
plotmeanHeatmapCytoband(subset, variables, named_palette, measures_long[i], measures_short[i], save = FALSE)






########################################
########################################
########################################


########################################## OV merging and analysis (mutation was for checking out ovarian)

ov_data = data_17_cohorts[data_17_cohorts$ProjectID == 'TCGA-OV',]
hrd_results = data.frame(read.csv('../../HRD_score/data/HRD_scores_pan_cancer_annotated_v2.csv'))
hrd_results_ov = hrd_results[hrd_results['Project.ID'] == 'TCGA-OV',]
hrd_results_primary = hrd_results_ov[hrd_results_ov['Type'] == 'Primary',]
mutation_data = read.csv('../data/mutations.txt', sep = '\t', header = TRUE)
methylation_data = read.csv('../data/Methylation (HM27).txt', sep = '\t', header = TRUE)
colData_ov = read.csv('../../RNAseq_pancancer/RNAseq_OV/results_Kmeans/colData_scores.csv', sep = ',', row.names = 1)

hrd_results_primary$SAMPLE_ID <- substr(hrd_results_primary$Sample.ID, 1, nchar(hrd_results_primary$Sample.ID) - 1)
mut_unique <- subset(mutation_data, !duplicated(SAMPLE_ID))
methylation_data_unique <- subset(methylation_data, !duplicated(SAMPLE_ID))


####### Prepare HRD_results
hrd_results_primary_subset = hrd_results_primary[,c('File.Name','SAMPLE_ID')]

### Prepare CIN Results
pattern <- "^chr[1-9][0-9]?_mean_seg_length$"
variables = grep(pattern, names(ov_data), value = TRUE)
ov_data$total_mean_seg_length <- rowSums(ov_data[variables])
ov_data$HRDtype = ifelse(ov_data$HRDtype == 'High', 1, 0)

subset_cin_data = ov_data[, c('File.Name', "HRDtype",'HRD_sum','TAI','LST','LOH','general_cin','cn_cin','totalCIN',
                     'n_cn_loh', 'n_cn_gain', 'n_cn_amp', 'n_homo_del', 'n_hemi_del', 'n_gain', 'n_amp','total_mean_seg_length')]

##### Prepare colData
colData_ov_subset = colData_ov[,c('File.Name','scoreCINSARC','scoreHRR','scoreMMEJ','scoreSSA','scoreHRRv2')]

##### Prepare Mutation data
mut_unique$HRR_genes_mut <- as.integer(rowSums(mut_unique[, 3:ncol(mut_unique)] != "WT") > 0)
mut_unique$BRCA1_mut = ifelse(mut_unique$BRCA1 != 'WT', 1, 0)
mut_unique$BRCA2_mut = ifelse(mut_unique$BRCA2 != 'WT', 1, 0)
mut_unique$BRCA1_2_mut = ifelse((mut_unique$BRCA2 != 'WT' & mut_unique$BRCA1 != 'WT'), 1, 0)

mut_unique_subset = mut_unique[,c('SAMPLE_ID','HRR_genes_mut','BRCA1_mut','BRCA2_mut','BRCA1_2_mut')]

methylation_data_unique = methylation_data_unique[,2:ncol(methylation_data_unique)]

merged_hrd_mut <- merge(hrd_results_primary_subset, mut_unique_subset, by = "SAMPLE_ID")
merged_hrd_mut_meth <- merge(merged_hrd_mut, methylation_data_unique, by = "SAMPLE_ID")
merge_hrd_mut_meth_cin <- merge(merged_hrd_mut_meth, subset_cin_data, by = "File.Name")
merge_hrd_mut_meth_cin_colData <- merge(merge_hrd_mut_meth_cin, colData_ov_subset, by = "File.Name")


write.csv(merge_hrd_mut_meth_cin_colData,'../data/ov_merged_data.csv', row.names = TRUE)

brca_mut_mut = mut_unique_subset[mut_unique_subset$BRCA1_2_mut == 1,]

cor.test(merge_hrd_mut_meth_cin_colData$scoreCINSARC, merge_hrd_mut_meth_cin_colData$HRD_sum)

#####################################
#####################################
#####################################
#####################################
#####################################

























