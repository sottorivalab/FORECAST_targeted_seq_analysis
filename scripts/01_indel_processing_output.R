# Calculate indels that we consider to be called

library(vcfR)
library(reshape2)
library(dplyr)

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
max_cohort_freq = 1
# Minimum mappability
min_map = 1
# Min QUAL score
min_qual = 100
# Min VAF
min_vaf = 0.05

# Do we assess the SSARs filtered files?
ssars_filter = TRUE

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
files = findFilesTES(file_pattern = paste0(file_en,"$"), sub_dir = sub_dir, pipeline_root = pipe_root,
                     recursive_in_subdir = F)

# Find the deepSNV files
deepSNV_files = findFilesTES(file_pattern = "_deepSNV_mutations_chr", sub_dir = "deepSNV/mutations", 
                             pipeline_root = pipe_root, recursive_in_subdir = F)

# Read in the mutation files
mutations = lapply(files, function(f) {
  
  extract.indels(read.vcfR(f, verbose = F), return.indels = T)
  
})

# Get the deepSNV results per chromosome, takes a minute
deepSNV_per_chr_mutations = lapply(deepSNV_files, function(f) {
  
  muts = read.table(f, header = T)
  
  muts = paste0(muts$chr,"_",muts$pos)
  
  return(muts)
  
})

# Read in the coverage data
coverage = lapply(mutations, function(s) {
  
  extract.gt(s, element = "NR", as.numeric = T)
  
})

# Read in the number of variant alleles
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

# Create vafs table
variant_allele_freq = lapply(1:length(mutations), function(i) {
  
  number_variants[[i]] / coverage[[i]]
  
})

# Get mappabilities
mutation_mappability = lapply(files, function(f) {
  
  patient = gsub(file_en, "", basename(f))
  
  map_data = read.table(paste0(pipe_root,
                               patient,"/",sub_dir,"/",patient,mut_map))
  
})

# Read in the metrics for passing QC or not
qc_metrics = read.table("refs/postUMI_qc_metrics.txt", 
                        header = T, stringsAsFactors = F)

# Summarise the occurence of different mutations across the cohort for filtering SSARs reads
cohort_pos_freq = table(unlist(lapply(variant_allele_freq, rownames)))

# Read the supp data
mut_data = read.table("refs/41588_2018_78_MOESM4_ESM.txt", 
                      header = T, sep = "\t", skip = 1, fill = NA, stringsAsFactors = F)

# Subset just the panel
panel_data = mut_data[mut_data$Hugo_Symbol %in% genes,]

# We don't want mets
panel_data = panel_data[panel_data$Primary_Met=="Primary",]

# Strip out the protein position
panel_data$protein_site = unlist(lapply(strsplit(panel_data$Protein_Change, split = ""), 
                                        function(i) as.numeric(paste0(i[-c(1:3,length(i))], collapse = ""))))

# Get all the mutated positions reported
prostate_mut_pos = unique(paste0(panel_data$Hugo_Symbol,"_",panel_data$protein_site))

# Gather mutation data
mutation_data = lapply(1:length(mutations), function(s) {
  
    # Get vafs from list
    vafs = variant_allele_freq[[s]]
    
    if(nrow(vafs) > 0) {
    
      # Get other features
      maps = mutation_mappability[[s]]
      muts = mutations[[s]]
      filt = mutation_filters[[s]]
      
      # Get alts
      alts = getALT(mutations[[s]])
      
      # Which just have one alt (like they are supposed to)
      single_alt = !grepl(",", alts)
      
      # Get other annotations
      gene           = unlist(lapply(strsplit(getINFO(muts), split = "[|]"), function(i) paste0(i[4])))
      protein_change = unlist(lapply(strsplit(getINFO(muts), split = "[|]"), function(i) paste0(i[4]," ",i[15]," ",i[16])))
      protein_pos    = unlist(lapply(strsplit(getINFO(muts), split = "[|]"), function(i) paste0(i[4],"_",i[15])))
      known_prostate = as.numeric(protein_pos %in% prostate_mut_pos)
      protein_coding = unlist(lapply(strsplit(getINFO(muts), split = "[|]"), function(i) i[16]))!=""
      
      # What's the patient name?
      patient     = gsub(file_en, "", basename(files[s]))
      
      # Which have deepSNV calls too
      whitelist   = unlist(deepSNV_per_chr_mutations[grep(patient, deepSNV_files)])
      
      # Vaf filter
      one_more_than_5_perc = apply(vafs, 1, function(r) any(r >= min_vaf))
      
      # Convert name
      white_search_name = unlist(lapply(strsplit(rownames(vafs), "_"), function(i) paste0(i[1],"_",as.numeric(i[2])+1)))
      
      # Output information for filtering indels
      plot_df = data.frame(vafs,
                           Patient = patient,
                           Mappability = maps[,6][match(rownames(vafs), maps[,1])],
                           Cohort_Frequency = as.numeric(cohort_pos_freq[rownames(vafs)]),
                           QUAL = mutation_qual[[s]],
                           QD = as.numeric(unlist(lapply(strsplit(filt, split = ";"), function(i) "QD" %in% i))),
                           Known_Mut_Site = known_prostate,
                           deepSNV = white_search_name %in% whitelist,
                           vaf_filt = one_more_than_5_perc,
                           Gene = gene,
                           Protein = protein_change,
                           protein_coding, 
                           single_alt = single_alt, check.names = F)
    
  } else {plot_df = NULL}
  
  return(plot_df)
  
})

# Add patient names to mutation data list
names(mutation_data) = unlist(lapply(strsplit(basename(files), split = "_"), function(i) i[1]))

# Remove those on the blacklist
mutation_data = mutation_data[!names(mutation_data) %in% blacklist_patients]

# Show the user which mutations have survived filtering
do.call(rbind, lapply(mutation_data, function(i) {
  
  if(!is.null(i)) {
    
    # Apply filter
    res = i[i$Mappability >= min_map & i$Cohort_Frequency <= max_cohort_freq & i$vaf_filt & i$protein_coding & 
            i$QUAL >= min_qual & i$single_alt,c("Patient","Mappability","Cohort_Frequency", "QUAL", "QD", 
                                                "Known_Mut_Site", "deepSNV", "vaf_filt", "Gene", "Protein", 
                                                "protein_coding", "single_alt")]
    
  } else {res = NULL}
  
}))

#################################################################################################################
##### 3 - Assess subclonality of mutations in the cohort ########################################################
#################################################################################################################

# Filter data by cohort frequency to avoid SSARs reads
mut_data_deepSNV     = lapply(mutation_data, function(i) i[i$Cohort_Frequency <= max_cohort_freq,])

# Remove those with no mutations called
mut_data_deepSNV     = mut_data_deepSNV[!unlist(lapply(mut_data_deepSNV, is.null))]

# Scroll through and summarise the status of clonality 
subclonality_summary = lapply(genes, function(g) {
  
  # Assess subclonality
  gene = lapply(mut_data_deepSNV, function(i) {
    
    # Remove those considered normal
    i = i[,!(grepl("N[0-9]{1,2}_DW1", colnames(i)) | grepl("BC1_DNA1", colnames(i)))]
    if(i$Patient[1] %in% forced_out) {i = i[,!colnames(i) %in% forced_norm]}
    
    # as.data.frame to capture single sample cases
    ans = as.data.frame(i[i$Gene==g & i$protein_coding & grepl("/", i$Protein) & i$vaf_filt & i$Mappability >= min_map & 
                          i$QUAL >= min_qual & i$single_alt,1:(ncol(i)-12)])
    
    # correct R's annoying conversion to numeric
    if(ncol(ans) == 1) {
      rownames(ans) = rownames(i[i$Gene==g & i$protein_coding & grepl("/", i$Protein) & i$vaf_filt & i$Mappability >= min_map & 
                                 i$QUAL >= min_qual & i$single_alt,])
      colnames(ans) = colnames(i)[1]
    }
    
    return(ans)
    
  })
  
  # Remove those without any mutations
  gene = gene[unlist(lapply(gene, function(i) nrow(i) > 0))]

  # Is it subclonal?
  gene_res = lapply(gene, function(i) apply(i, 1, function(j) {length(which(j >=min_vaf)) < (length(j))}))
  
  return(gene_res)
  
})

# Name the output by gene
names(subclonality_summary) = genes

# Save that result
saveRDS(subclonality_summary, file = "results/mutation_calling/subclonality_assessment_indels.rds")

#################################################################################################################
##### 4 - Per sample annotation of mutation status per gene #####################################################
#################################################################################################################

# Read in patients where there was a normal
patient_w_normals = read.table("results/mutation_calling/Patients_with_normal_attempted_mut.txt", stringsAsFactors = F)$V1

# Filter the mutations according to the filters
passing_mutations = lapply(mutation_data, function(i) i[i$Cohort_Frequency <= max_cohort_freq & 
                                                        i$protein_coding & grepl("/", i$Protein) & 
                                                        i$vaf_filt & i$Mappability >= min_map & i$QUAL >= min_qual & i$single_alt,])

# Take only samples that passed the coverage filter
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
    
    # If the patient has mismatched normal then record NA
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

# Write it out for use in the main script
write.table(per_sample_mutation, file = "results/mutation_calling/Per_sample_gene_mutation_status_indels.txt", 
            sep = "\t", quote = F, row.names = F)

# Write out the passing mutations object for LOH analysis
saveRDS(passing_mutations, file = "results/mutation_calling/Passing_mutations_list_indels.rds")

# Write out a melted version of the passing mutations that contains the mutation itself
per_sample_mutation_with_change = do.call(rbind, lapply(passing_mutations, function(i) {
  
  if(is.null(i)) {
    
    output = NULL
    
  } else {output = melt(i[,!colnames(i) %in% c("Mappability", "Cohort_Frequency", "QUAL", "QD", "Known_Mut_Site", "single_alt")])}
  
  return(output)
  
}))

# Format the rows and columns
rownames(per_sample_mutation_with_change) = NULL
colnames(per_sample_mutation_with_change)[7:8] = c("Sample", "VAF")

# Output this result too
write.table(per_sample_mutation_with_change, file = "results/mutation_calling/Per_sample_gene_mutation_status_with_change_indels.txt",
            sep = "\t", quote = F, row.names = F)
