source("scripts/00_general_functions.R")

library(reshape2)
library(dplyr)

inputs = c("results/mutation_calling/Passing_mutations_list_snvs.rds",
           "results/mutation_calling/Passing_mutations_list_indels.rds")

# Manipulate per patient the mutations
passing_mutations = lapply(inputs, function(f) {
  
  passing_mutations = readRDS(f)
  
  # Write out a melted version of the passing mutations that contains the mutation itself
  per_sample_mutation_with_change = do.call(rbind, lapply(passing_mutations, function(i) {
    
    if(is.null(i)) {
      
      output = NULL
      
    } else {output = melt(cbind(rownames(i), i[,!colnames(i) %in% c("Mappability", "Cohort_Frequency", "QUAL", "QD", "Known_Mut_Site", "single_alt")]))}
    
    return(output)
    
  }))
  
  return(per_sample_mutation_with_change)

})

# Combine outputs
per_sample_mutation_with_change = do.call(rbind, passing_mutations)

# Rename the location names as row names
rownames(per_sample_mutation_with_change) = NULL

# Appropriately change column names
colnames(per_sample_mutation_with_change)[1] = "Pos"
colnames(per_sample_mutation_with_change)[8:9] = c("Sample", "VAF")

# Class change
per_sample_mutation_with_change$Pos = as.character(per_sample_mutation_with_change$Pos)

# Extract chromosome and position
per_sample_mutation_with_change = cbind(do.call(rbind, strsplit(per_sample_mutation_with_change$Pos, split = "_")), 
                                        per_sample_mutation_with_change)

# Remove the combo column now
per_sample_mutation_with_change$Pos = NULL

# Give proper names
colnames(per_sample_mutation_with_change)[1:2] = c("chr", "position")

# Gather the purity estimates from the low coverage WGS data
per_sample_mutation_with_change$Purity = apply(per_sample_mutation_with_change, 1, function(i) {
  
  # What's the file?
  file = paste0("results/0_min_het_best_fit_ploidy_search/calls/",
                i["Patient"],"_",i["Sample"],"_cna_ploidy_search_metrics.txt")
  
  # Does the file exist?
  output = ifelse(file.exists(file), read.table(file, header = T, sep = "\t")$Purity, NA)
  
  return(output)
  
})

# Gather the position CN
per_sample_mutation_with_change$CN = apply(per_sample_mutation_with_change, 1, function(i) {
  
  # What file?
  file = paste0("results/0_min_het_best_fit_ploidy_search/calls/",
                i["Patient"],"_",i["Sample"],"_cna_ploidy_search_calls.txt")
  
  # If the file exists, get the CN
  if(file.exists(file)) {
    
    # Read in the cnas
    cnas = read.table(file, header = T)
    
    # Remove the text
    mut_chr = gsub("chr", "", i["chr"])
    
    # Get mutation position
    mut_pos = as.numeric(i["position"])
    
    # What the bin locations?
    cna_locs = extractChrRows(rownames(cnas))
    
    # Which bin does it lie in?
    hit = cnas[cna_locs$chr == mut_chr & cna_locs$start < mut_pos & cna_locs$end > mut_pos,]
    
  # Else just return NA  
  } else {hit = NA}
  
  # If there is no overlapping bin return NA
  if(length(hit) == 0) {hit = NA}
  
  return(hit)
  
})


# Calculate copies using formula
per_sample_mutation_with_change$Copies = calculateCopies(V = per_sample_mutation_with_change$VAF, 
                                                         C = per_sample_mutation_with_change$CN, 
                                                         r = per_sample_mutation_with_change$Purity,
                                                         exp_norm_cn = ifelse(per_sample_mutation_with_change$chr %in% c("chrX", "chrY"), 1, 2))

# Calculate LOH using extended function
per_sample_mutation_with_change$LOH = isitLOH(V = per_sample_mutation_with_change$VAF, 
                                              C = per_sample_mutation_with_change$CN, 
                                              r = per_sample_mutation_with_change$Purity,
                                              exp_norm_cn = ifelse(per_sample_mutation_with_change$chr %in% c("chrX", "chrY"), 1, 2))


# Delete columns
per_sample_mutation_with_change$deepSNV  = NULL
per_sample_mutation_with_change$vaf_filt = NULL
per_sample_mutation_with_change$protein_coding = NULL
per_sample_mutation_with_change$Purity   = NULL

per_sample_mutation_with_change$VAF    = signif(per_sample_mutation_with_change$VAF, digits = 3)
per_sample_mutation_with_change$Copies = signif(per_sample_mutation_with_change$Copies, digits = 3)

# Reorder
per_sample_mutation_with_change = per_sample_mutation_with_change[,c("Patient", "Sample", "chr", "position", 
                                                                     "Gene", "Protein", "VAF", "CN", "Copies", "LOH")]
colnames(per_sample_mutation_with_change)[3:4] = c("Chr", "Position")

# Order it by patient
per_sample_mutation_with_change = per_sample_mutation_with_change[order(per_sample_mutation_with_change$Patient),]

# Write it out
write.csv(per_sample_mutation_with_change, 
          file = "results/mutation_calling/Per_sample_gene_mutation_status_with_change_zygosity.csv", 
          quote = F, row.names = F)
