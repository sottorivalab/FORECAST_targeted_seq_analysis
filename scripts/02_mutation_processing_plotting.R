# Calculate SNVs we consider real plus some plotting and data outputting

library(ggplot2)
library(cowplot)
library(vcfR)
library(ComplexHeatmap)
library(circlize)
library(reshape2)
library(dplyr)
theme_set(theme_cowplot())

source("scripts/00-general_functions.R")

#################################################################################################################
##### 1 - Parameter definitions #################################################################################
#################################################################################################################

# Where is the pipeline?
pipe_root = "targeted_snakemake/"

# Which patients are blacklisted for mutation analysis?
blacklist_patients = c("", "")
# Which sample has posthoc been considered to be not tumour?
forced_out  = c("")
forced_norm = c("")
# Minimum number of samples
min_samples = 3
# Cohort frequency filtering, make it Inf to not subset
max_cohort_freq = 3
# Minimum mappability
min_map = 0.6
# Min VAF
min_vaf = 0.05
# Minimum qual filter
min_qual = 10

# Do we assess the SSARs filtered files?
ssars_filter = TRUE

# Recalculate patients with attempted normal?
redo_patients_with_normal = F

#################################################################################################################
##### 2 - Read in data and begin processing #####################################################################
#################################################################################################################

# Depends on filter
if(ssars_filter) {
  
  sub_dir = "filterSSARs"
  file_en = "_platypus_genotyped_annotated_ssars.vcf.gz"
  mut_map = "_mutation_mappability_ssars.out"
  
} else {
  
  sub_dir = "platypus_genotyping_calling"
  file_en = "_platypus_genotyped_annotated.vcf.gz"
  mut_map = "_mutation_mappability.out"
  
}

# Get the cohort patients
all_patients = list.files(pipe_root)
all_patients = c(all_patients[grep("^DNT", all_patients)], all_patients[grep("^IM", all_patients)])

# Find the mutation files
files = findFilesTES(file_pattern = paste0(file_en,"$"), sub_dir = sub_dir, pipeline_root = pipe_root, recursive_in_subdir = F)

# Find the deepSNV files
deepSNV_files = findFilesTES(file_pattern = "_deepSNV_mutations_chr", sub_dir = "deepSNV/mutations", 
                             pipeline_root = pipe_root, recursive_in_subdir = F)

# Read in the mutation files
mutations = lapply(files, function(f) {
  
  extract.indels(read.vcfR(f, verbose = F), return.indels = F)
  
})

# Get the deepSNV results per chromosome, takes a minute
deepSNV_per_chr_mutations = lapply(deepSNV_files, function(f) {
  
  muts = read.table(f, header = T)
  
  muts = paste0(muts$chr,"_",muts$pos)
  
  return(muts)
  
})

# Read in the coverage per mutation
coverage = lapply(mutations, function(s) {
  
  extract.gt(s, element = "NR", as.numeric = T)
  
})

# Read in the number of variants per mutation
number_variants = lapply(mutations, function(s) {
  
  extract.gt(s, element = "NV", as.numeric = T)
  
})

# Read in the filter flags
mutation_filters = lapply(mutations, function(s) {
  
  getFILTER(s)
  
})

# Read in the mutation quality scores
mutation_qual = lapply(mutations, function(s) {
  
  getQUAL(s)
  
})

# Create vafs
variant_allele_freq = lapply(1:length(mutations), function(i) {
  
  number_variants[[i]] / coverage[[i]]
  
})

# Get mappabilities per mutation
mutation_mappability = lapply(files, function(f) {
  
  patient = gsub(file_en, "", basename(f))
  
  map_data = read.table(paste0(pipe_root,
                               patient,"/",sub_dir,"/",patient,mut_map))
  
})

# Read in the metrics for passing QC or not
qc_metrics = read.table("refs/postUMI_qc_metrics.txt", 
                        header = T, stringsAsFactors = F)

# Summarise the occurence of different mutations across the cohort for filtering
cohort_pos_freq = table(unlist(lapply(variant_allele_freq, rownames)))

# Read the supp data
mut_data = read.table("refs/41588_2018_78_MOESM4_ESM.txt", 
                      header = T, sep = "\t", skip = 1, fill = NA, stringsAsFactors = F)

# Subset just the panel
panel_data = mut_data[mut_data$Hugo_Symbol %in% genes,]

# We don't want mets
panel_data = panel_data[panel_data$Primary_Met=="Primary",]

# Only point mutations
panel_data = panel_data[panel_data$type %in% c("Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation"),]

# Strip out the protein position
panel_data$protein_site = unlist(lapply(strsplit(panel_data$Protein_Change, split = ""), 
                                        function(i) as.numeric(paste0(i[-c(1:3,length(i))], collapse = ""))))

# Get all the mutated positions reported
prostate_mut_pos = unique(paste0(panel_data$Hugo_Symbol,"_",panel_data$protein_site))

#################################################################################################################
##### 3 - Make dataframes of the data summary ###################################################################
#################################################################################################################

# Collect mutation dataframes
mutation_data = lapply(1:length(mutations), function(s) {
  
  # Get data from lists we read in
  vafs = variant_allele_freq[[s]]
  maps = mutation_mappability[[s]]
  muts = mutations[[s]]
  filt = mutation_filters[[s]]
  
  # Annotation data
  gene           = unlist(lapply(strsplit(getINFO(muts), split = "[|]"), function(i) paste0(i[4])))
  protein_change = unlist(lapply(strsplit(getINFO(muts), split = "[|]"), function(i) paste0(i[4]," ",i[15]," ",i[16])))
  protein_pos    = unlist(lapply(strsplit(getINFO(muts), split = "[|]"), function(i) paste0(i[4],"_",i[15])))
  known_prostate = as.numeric(protein_pos %in% prostate_mut_pos)
  protein_coding = unlist(lapply(strsplit(getINFO(muts), split = "[|]"), function(i) i[16]))!=""
  
  # Get patient 
  patient     = gsub(file_en, "", basename(files[s]))
  
  # Which ones were called by deepSNV
  whitelist   = unlist(deepSNV_per_chr_mutations[grep(patient, deepSNV_files)])
  
  # Vaf filter  
  one_more_than_5_perc = apply(vafs, 1, function(r) any(r >= min_vaf))
  
  # Get data
  plot_df = data.frame(vafs,
                       Patient = patient,
                       Mappability = maps[,6][match(rownames(vafs), maps[,1])],
                       Cohort_Frequency = as.numeric(cohort_pos_freq[rownames(vafs)]),
                       QUAL = mutation_qual[[s]],
                       QD = as.numeric(unlist(lapply(strsplit(filt, split = ";"), function(i) "QD" %in% i))),
                       Known_Mut_Site = known_prostate,
                       deepSNV = rownames(vafs) %in% whitelist,
                       vaf_filt = one_more_than_5_perc,
                       Gene = gene,
                       Protein = protein_change,
                       protein_coding, check.names = F)
  
  return(plot_df)
  
})

#################################################################################################################
##### 4 - Filter for samples we can't use #######################################################################
#################################################################################################################

# Add patient names for processing
names(mutation_data) = unlist(lapply(mutation_data, function(i) unique(i$Patient)))

# To be excluded
exc_patient = c(names(which(!unlist(lapply(mutation_data, function(i) {
  
  # Remove normal samples
  i = i[,!(grepl("N[0-9]{1,2}_DW1", colnames(i)) | grepl("BC1_DNA1", colnames(i)))]
  if(i$Patient[1] %in% forced_out) {i = i[,!colnames(i) %in% forced_norm]}
  
  # Count number of samples
  ans = ncol(i)-11
  
})) >= min_samples)))

# Write that out
write.table(exc_patient, file = "results/mutation_calling/Patients_excluded_TES_too_few_samples.txt", 
            quote = F, col.names = F, row.names = F)
write.table(blacklist_patients, file = "results/mutation_calling/Patients_excluded_TES_mismatch.txt", 
            quote = F, col.names = F, row.names = F)

# Remove those on the blacklist
mutation_data = mutation_data[!names(mutation_data) %in% blacklist_patients]

#################################################################################################################
##### 5 - Collapse per mutation summaries and plot per mutation bar plot ########################################
#################################################################################################################

# Extract the summary columns for each mutation
cohort = do.call(rbind, lapply(mutation_data, function(m) m[,c("Patient", "Cohort_Frequency", "Gene", "Protein", 
                                                               "protein_coding", "Known_Mut_Site", "deepSNV",
                                                               "vaf_filt", "Mappability", "QUAL")]))

# Mutations that are overly recurrent in this cohort which on manual inspection are SSARs reads
ssars_artefacts = cohort[cohort$deepSNV & cohort$vaf_filt & cohort$Known_Mut_Site==0 & cohort$Mappability > min_map & 
                   cohort$Cohort_Frequency > max_cohort_freq & grepl("/", cohort$Protein) & cohort$QUAL >= min_qual,]

# Get pos again
ssars_artefacts$Pos = unlist(lapply(strsplit(rownames(ssars_artefacts), split = "[.]"), function(i) i[2]))

# Get unique ones
ssars_artefacts = ssars_artefacts[!duplicated(ssars_artefacts$Pos),]

# Subset
cohort = cohort[grepl("/", cohort$Protein) &
                cohort$protein_coding & 
                cohort$vaf_filt &
                cohort$Cohort_Frequency <= max_cohort_freq & 
                cohort$Mappability > min_map &
                cohort$deepSNV &
                cohort$QUAL >= min_qual,]

# Barplot df
bar_df = melt(table(cohort$Gene, cohort$Known_Mut_Site))
bar_df$Var2 = as.factor(bar_df$Var2)
bar_df$Var1 = factor(bar_df$Var1, levels = names(sort(table(cohort$Gene), decreasing = T)))

# Read in the expected frequencies
mut_freq = read.table("results/mutation_calling/Expected_cohort_frequencies_from_single_sample_cohorts.txt",
                      header = T)

# Create the barplot and show expected number of mutations
ggplot(data=bar_df, 
       aes(x=Var1, y = value, fill=Var2)) +
  geom_bar(stat="identity") +
  ggtitle(paste0("Frequency observed (SSAR filter=",ssars_filter,")")) +
  geom_point(data = mut_freq, aes(x = Var1, y = Freq*length(mutation_data)), fill = "black") +
  ylab("Number of mutations") +
  xlab("") +
  scale_fill_manual(values = c("0" = "#ff474c", "1" = "#0485d1")) +
  guides(fill=guide_legend(title="Reported\nMutations")) +
  theme_cowplot() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.title = element_text(hjust = 0.5))

#################################################################################################################
##### 6 - Assess subclonality of mutations in the cohort ########################################################
#################################################################################################################

# Apply just deepSNV and prevelance filters
mut_data_deepSNV     = lapply(mutation_data, function(i) i[i$deepSNV,])
mut_data_deepSNV     = lapply(mut_data_deepSNV, function(i) i[i$Cohort_Frequency <= max_cohort_freq,])

# Scroll through and summarise the status of clonality 
subclonality_summary = lapply(genes, function(g) {
  
  # Assess subclonality
  gene = lapply(mut_data_deepSNV, function(i) {
    
    # Which ones are normal?
    i = i[,!(grepl("N[0-9]{1,2}_DW1", colnames(i)) | grepl("BC1_DNA1", colnames(i)))]
    if(i$Patient[1] %in% forced_out) {i = i[,!colnames(i) %in% forced_norm]}
    
    # as.data.frame to capture single sample cases
    ans = as.data.frame(i[i$Gene==g & i$protein_coding & grepl("/", i$Protein) & i$vaf_filt & i$Mappability > min_map & 
                          i$QUAL >= min_qual,1:(ncol(i)-11)])
    
    # correct R's annoying conversion to numeric
    if(ncol(ans) == 1) {
      rownames(ans) = rownames(i[i$Gene==g & i$protein_coding & grepl("/", i$Protein) & i$vaf_filt & i$Mappability > min_map & 
                                 i$QUAL >= min_qual,])
      colnames(ans) = colnames(i)[1]
    }
    
    return(ans)
    
  })
  
  # Remove those with no mutation
  gene = gene[unlist(lapply(gene, function(i) nrow(i) > 0))]

  # Is it in fewer than every sample?
  gene_res = lapply(gene, function(i) apply(i, 1, function(j) {length(which(j >=min_vaf)) < (length(j))}))
  
  return(gene_res)
  
})

# Add names of genes
names(subclonality_summary) = genes

# Make this the snv result
subclonality_summary_snv = subclonality_summary

# Gather the indel result
subclonality_summary_indel = readRDS("results/mutation_calling/subclonality_assessment_indels.rds")

# Add in the indel mutations
for(g in names(subclonality_summary_indel)) {
  
  for(p in names(subclonality_summary_indel[[g]])) {
    
    subclonality_summary[[g]][[p]] = c(subclonality_summary[[g]][[p]], subclonality_summary_indel[[g]][[p]])
    
  }
  
}

# Final result with both indel and snv
saveRDS(subclonality_summary, file = "results/mutation_calling/subclonality_assessment.rds")

#################################################################################################################
##### 7 - Summarise in a heatmap the mutation status per patient and show missing data ##########################
#################################################################################################################

# Use the pipeline to find those patients were mutation calling was attempted
if(redo_patients_with_normal) {

  # Get files with attempted mutation calling
  files_first = findFilesTES(file_pattern = paste0("_platypus_genotyping_all_samples_ssars.vcf.gz","$"), 
                             sub_dir = sub_dir, pipeline_root = pipe_root)
  
  # Patients with normals taken from mutation calling output
  patient_w_normals = extractPatientName(basename(files_first))
  
  # Write that out
  write.table(patient_w_normals, file = "results/mutation_calling/Patients_with_normal_attempted_mut.txt", 
              quote = F, col.names = F, row.names = F)
  
} else {patient_w_normals = read.table("results/mutation_calling/Patients_with_normal_attempted_mut.txt", stringsAsFactors = F)$V1}

# Create empty vector for collection of results
res = NULL

# Run through each patient for making data for plotting
for(p in all_patients) {
  
  # Collection vector
  p_col = NULL
  
  # Run through the genes
  for(g in genes) {
    
    # Is the patient and gene there?
    pg = ifelse(length(subclonality_summary[[g]][[p]])>0, 1, 0)
    
    # If the patient isn't there then record NA
    pg = ifelse(p %in% patient_w_normals, pg, NA)
    
    # If the patient is blacklist, then record that
    pg = ifelse(!p %in% blacklist_patients, pg, NA)
    
    # Add something for subclonality
    if(!is.null(subclonality_summary[[g]][[p]])) {
      
      if(all(subclonality_summary[[g]][[p]])) {pg = 0.5}
      
    }
    
    # Add something for subclonality
    if(!is.null(subclonality_summary_indel[[g]][[p]])) {
      
      pg = -1*pg
      
      if(!is.null(subclonality_summary_snv[[g]][[p]])) {
        
        pg = pg-2
        
      }
      
    }
    
    # Add it on
    p_col = c(p_col, pg)
    
  }
  
  # Final collection
  res = cbind(res, p_col)
  
}

# Add names for heatmap
colnames(res) = all_patients
rownames(res) = genes

# Make res to plot
res_plot = res

# Convert res back to original style
res = ifelse(res < -1, res+2, res)

# Collapse frequencies
obs = data.frame(Var1 = names(apply(na.omit(t(res)), 2, function(i) sum(ceiling(abs(i)))) / nrow(t(res))),
                 Observed = apply(na.omit(t(res)), 2, function(i) sum(ceiling(abs(i)))) / nrow(t(res)))

# Compare to published data
freq_compare = merge(mut_freq, obs, by = "Var1")

# Make a plot
ggplot(freq_compare, aes(x = Freq, y = Observed)) + geom_point() + 
  ggtitle("Mutations observed per patient") +
  ylab("Observed Frequency (multiregion data)") +
  xlab("Expected Frequency (single sample)") + 
  geom_abline(slope = 1, intercept = 0) + theme(plot.title = element_text(hjust = 0.5))
ggsave("results/mutation_calling/Expected_per_patient_frequency_vs_FORECAST.pdf", width = 5, height = 5)

# Order by the frequency
res_plot = res_plot[order(apply(res, 1, function(i) sum(ceiling(abs(i)), na.rm = T)), decreasing = T),]
res      = res[order(apply(res, 1, function(i) sum(ceiling(abs(i)), na.rm = T)), decreasing = T),]

# Dplyr to sort
col_order = colnames(t(data.frame(t(abs(res))) %>% arrange(across(1:nrow(res), desc))))

# Add on frequency information
rownames(res_plot) = paste0(rownames(res)," (n=",apply(res, 1, function(i) sum(ceiling(abs(i)), na.rm = T)),")")

# Write out res because it is useful
write.csv(res_plot, file = "results/mutation_calling/Per_patient_mutation_status_heatmap_data.csv", quote = F)

# Make heatmap
pdf("results/mutation_calling/Patient_summary_heatmap.pdf", width = 16, height = 7)
hm = Heatmap(res_plot[,col_order], cluster_rows = F, cluster_columns = F, rect_gp = gpar(col = "white", lwd = 2),
        col = colorRamp2(c(-3, -2, -1, 0, 1), c("orange", "white", "#3f9b0b", "white", "#e50000")), border = T, show_heatmap_legend = F)
print(hm)
dev.off()

#################################################################################################################
##### 8 - Collect the mutation status of each sample individually ###############################################
#################################################################################################################

# Filter the mutations according to the filters
passing_mutations = lapply(mutation_data, function(i) i[i$deepSNV & i$Cohort_Frequency <= max_cohort_freq & 
                                                           i$protein_coding & grepl("/", i$Protein) & 
                                                           i$vaf_filt & i$Mappability > min_map & i$QUAL >= min_qual,])

# Take only sample
samples_of_interest = qc_metrics[qc_metrics$MEAN_TARGET_COVERAGE >= 10,c("Patient", "Sample")]

# Get sample id only
samples_of_interest$Sample = apply(samples_of_interest, 1, function(i) gsub(paste0(i[1],"_"), "", i[2]))

# Collect the data
per_sample_mutation = NULL

# Per row samples
for(sample in 1:nrow(samples_of_interest)) {
  
  # Get patient and sample ids
  pat = samples_of_interest[sample,"Patient"]
  sam = samples_of_interest[sample,"Sample"]
  
  # Collect per gene result
  gene_res = NULL
  
  for(g in genes) {
    
    # Get the sample and gene from patient
    mut = passing_mutations[[pat]][passing_mutations[[pat]]$Gene==g,sam]
    
    # Is the patient and gene there?
    pg = ifelse(length(mut)>0, ifelse(any(mut>=min_vaf), T, F), F)
    
    # If the patient isn't there then record NA
    pg = ifelse(pat %in% patient_w_normals, pg, NA)
    
    # If the patient is blacklist then record NA
    pg = ifelse(!pat %in% blacklist_patients, pg, NA)
    
    # Collect result for this gene
    gene_res = c(gene_res, pg)
    
  }
  
  # Name the results
  names(gene_res) = genes
  
  # Combine it
  per_sample_mutation = rbind(per_sample_mutation, gene_res)
  
}

# Associate it with the sample names
per_sample_mutation = cbind(samples_of_interest, per_sample_mutation)

# Read in the indel result
per_sample_mutation_indel = read.table("results/mutation_calling/Per_sample_gene_mutation_status_indels.txt", 
                                       header = T, sep = "\t", stringsAsFactors = F)

# Add on the indel status
per_sample_mutation = cbind(per_sample_mutation[,1:2], 
                            per_sample_mutation_indel[,3:ncol(per_sample_mutation_indel)] | 
                              per_sample_mutation[,3:ncol(per_sample_mutation)])

# Write it out for use in the main script
write.table(per_sample_mutation, file = "results/mutation_calling/Per_sample_gene_mutation_status.txt", 
            sep = "\t", quote = F, row.names = F)

# Gather those patients that have a subclonal driver (not KMTs)
subclonal_driver = NULL

# Run through patients
for(p in all_patients) {
  
  # Collection point
  sc_muts = NULL
  
  # We won't use KMT subclonal mutants
  for(g in genes[!grepl("KMT", genes)]) {
    
    # Are any listed as subclonal for this gene and patient?
    sc_mut = any(subclonality_summary[[g]][[p]])
    
    # Collect for patient 
    sc_muts = c(sc_mut, sc_muts)
  }
  
  # Did the patient have any subclonal mutations?
  sc_muts = any(sc_muts)
  
  # If the patient isn't there then record NA
  sc_muts = ifelse(p %in% patient_w_normals, sc_muts, NA)
  
  # If the patient is blacklist, then record that
  sc_muts = ifelse(!p %in% blacklist_patients, sc_muts, NA)
  
  # Collect results
  subclonal_driver = rbind(subclonal_driver, 
                            data.frame(Patient = p, Subclonal_Mut_Status = sc_muts, stringsAsFactors = F))
  
}

# Write out as a result
write.table(subclonal_driver, file = "results/mutation_calling/Subclonal_driver_status_summary.txt",
            quote = F, row.names = F, sep = "\t")

# Write out the passing mutations object for LOH analysis
saveRDS(passing_mutations, file = "results/mutation_calling/Passing_mutations_list_snvs.rds")

# Write out a melted version of the passing mutations that contains the mutation itself
per_sample_mutation_with_change = do.call(rbind, 
                                          lapply(passing_mutations, 
                                                 function(i) melt(i[,!colnames(i) %in% 
                                                                     c("Mappability", "Cohort_Frequency", "QUAL", "QD", "Known_Mut_Site")])))

# Format
rownames(per_sample_mutation_with_change) = NULL
colnames(per_sample_mutation_with_change)[7:8] = c("Sample", "VAF")

# Output this result too
write.table(per_sample_mutation_with_change, file = "results/mutation_calling/Per_sample_gene_mutation_status_with_change.txt",
            sep = "\t", quote = F, row.names = F)
