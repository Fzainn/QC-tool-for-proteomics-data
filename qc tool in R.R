
#QC tool for limited proteolysis-coupled mass-spectrometric data (LIP-MS)
install.packages("PTXQC")

#to know where was PTX installed in my pc
cat(paste0("\nPTXQC was installed to '", .libPaths()[1], "'.\n\n"))

#####packages required for analysis

install.packages("protti")
install.packages("janitor")
install.packages("plotly")
install.packages("ggrepel")

library(conjurer)
library(protti)
library(tidyverse)
library(magrittr)
library(dplyr)

#we will create a random synthetic data set contains data similar to data obtained from a treatment exprermint with protein,
#metabolite or small-molecule
# random data set with 100 different proteins, out of which 5 % are significantly changing upon treatment.
#The data set includes 3 replicates for 2 different conditions (treated and untreated)
#QC should be conducted on raw unfiltered data (the direct output of your search engine of choice)

#to make sure random numbers can be reproduced or generated

set.seed(123) 

data <- create_synthetic_data(n_proteins = 100,
                              frac_change = 0.05,
                              n_replicates = 3,
                              n_conditions = 2,
                              method = "effect_random",
                              additional_metadata = TRUE)

#check QC of coefficient of variation between replicats r in reasonable range, it should be bellow 0.05
#if u see higher CVs that might be sample preparation error or the instrument performance wasn't ideal 

#qc_cvs() function that calculate coefficient variation of our data, this function works only with raw values
#so we will backtransform our log2 transformed data and create a new column called (raw_intensity)
#from the (peptide_inyinsity_missing) column which contains peptide_intensities.
# mutate() function to create a new column  

input <- data %>%
  # as the data is log2 transformed, we need to transform it back before calculating the CVs
  mutate(raw_intensity = 2^peptide_intensity_missin)

qc_cvs(data = input,
       grouping = peptide,
       condition = condition,
       intensity = raw_intensity, 
       plot = FALSE)

qc_cvs(data = input,
       grouping = peptide,
       condition = condition,
       intensity = raw_intensity, 
       plot = TRUE, 
       plot_style = "violin")
# in this plot condition-_1 and condition_2 has a higher cv than combined that means, 
#there is a problem with one or more sample or its noisy becuase of the treatment in this condition


#qc_ids function used to analysis of the number of identifications of precursor for peptides or protein
#the no.of protein should be the same for different samples. if if there r less observations in one sample, 
#this might a sample preparation or measurement error.

qc_ids(data = input,
       sample = sample,
       grouping = protein,
       intensity = peptide_intensity_missing,
       condition = condition, 
       plot = FALSE)

qc_ids(data = input,
       sample = sample,
       grouping = protein,
       intensity = peptide_intensity_missing,
       condition = condition, 
       title = "Protein identifications per sample",
       plot = TRUE)


#qc_peptide_type function to determine the type of peptides in our data
#which helps us in reproducibility of cleavage events in the experiment and if our treatment influenced protease activity or the digest in general
qc_peptide_type(data = input,
                sample = sample, 
                peptide = peptide, 
                pep_type = pep_type, 
                method = "intensity", 
                intensity = raw_intensity, 
                plot = TRUE, 
                interactive = FALSE)

qc_peptide_type(data = input,
                sample = sample, 
                peptide = peptide, 
                pep_type = pep_type, 
                method = "count", 
                plot = TRUE, 
                interactive = FALSE)


#run intensity with boxplot to help us to know if there r sample losses or measurment issues.

qc_intensity_distribution(data = input,
                          sample = sample,
                          grouping = peptide,
                          intensity_log2 = peptide_intensity_missing,
                          plot_style = "boxplot")
# lineplot to helps u to know if there r any trends in ur data
qc_median_intensities(data = input,
                      sample = sample,
                      grouping = peptide,
                      intensity = peptide_intensity_missing)
#charge state for ion should be the same for all samples
qc_charge_states(data = input, 
                 sample = sample, 
                 grouping = peptide, 
                 charge_states = charge, 
                 method = "intensity",
                 intensity = raw_intensity, 
                 plot = TRUE)



#check the numbers of missed cleavages in our dataset by using the function qc_missed_cleavages()
#this should be low

qc_missed_cleavages(data = input, 
                    sample = sample, 
                    grouping = peptide, 
                    missed_cleavages = n_missed_cleavage, 
                    method = "intensity",
                    intensity = raw_intensity, 
                    plot = TRUE)



#for sequence coverage, what is the precentage of the protein sequence is conered by identified peptides?
# it should be high
qc_sequence_coverage(data = input, 
                     protein_identifier = protein, 
                     coverage = coverage)


#the function qc_peak_width() In order to identify potential chromatographic issues that might have occurred during the measurement
#it should be the same in all samples

qc_peak_width(data = input,
              sample = sample, 
              intensity = peptide_intensity_missing,
              retention_time = retention_time, 
              peak_width = peak_width)

#The function qc_data_completeness() checks how many of all detected precursors, peptides or proteins were identified in each sample
qc_data_completeness(data = input, 
                     sample = sample, 
                     grouping = peptide, 
                     intensity = peptide_intensity_missing, 
                     plot = TRUE)
#For different kinds of analyses (e.g. t-tests) it is important that your data intensity follows a normal distribution. To ensure 
#that this is the case, we are going to use the function qc_intensity_distriubution()
qc_intensity_distribution(data = input,
                          grouping = peptide,
                          intensity_log2 = peptide_intensity_missing,
                          plot_style = "histogram")

install.packages("dendextend")
install.packages("pheatmap")
install.packages("seriation")
#Another approach to quality control is to check the correlation of your samples. Ideally,
#replicates should cluster together and different treatment conditions should be separated. 
#We are now going to check if this is the case for our data by using the function qc_sample_correlation()
qc_sample_correlation(data = input,
                      sample = sample, 
                      grouping = peptide, 
                      intensity_log2 = peptide_intensity_missing, 
                      condition = condition, 
                      interactive = FALSE)



#Another popular quality control method (which could also be considered part of data analysis) 
#is the principal component analysis (PCA). Generally, PCA is a method that reduces dimensionality of large datasets. 
#In our case this helps us to quickly assess how similar or different our replicates and conditions are
qc_pca(
  data = data,
  sample = sample,
  grouping = peptide,
  intensity = peptide_intensity_missing,
  condition = condition,
  digestion = NULL,
  plot_style = "scree"
)


qc_pca(
  data = data,
  sample = sample,
  grouping = peptide,
  intensity = peptide_intensity_missing,
  condition = condition,
  components = c("PC1", "PC2"), 
  plot_style = "pca"
)
