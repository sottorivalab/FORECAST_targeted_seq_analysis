# Script to perform dNdS analysis for the cohort

library(dndscv)

#################################################################################################################
############################ Run SNV section 'mutation_processing_plotting.R' ###################################
#################################################################################################################

# Collect mutation dataframes
dnds_data = lapply(1:length(mutations), function(s) {

  # Get the vafs
  vafs = variant_allele_freq[[s]]
  
  # Remove sample with an unusual enrichment of synonymous SNPs
  if(s==2) {vafs = vafs[,!colnames(vafs)%in%c("B11_DW1")]}
  
  # Get other features
  maps = mutation_mappability[[s]]
  muts = mutations[[s]]
  filt = mutation_filters[[s]]
  
  # Get ref and alts
  refs = getREF(mutations[[s]])
  alts = getALT(mutations[[s]])
  
  # Get annotation data
  gene           = unlist(lapply(strsplit(getINFO(muts), split = "[|]"), function(i) paste0(i[4])))
  protein_change = unlist(lapply(strsplit(getINFO(muts), split = "[|]"), function(i) paste0(i[4]," ",i[15]," ",i[16])))
  protein_pos    = unlist(lapply(strsplit(getINFO(muts), split = "[|]"), function(i) paste0(i[4],"_",i[15])))
  protein_coding = unlist(lapply(strsplit(getINFO(muts), split = "[|]"), function(i) i[16]))!=""
  
  # Extract patient names
  patient     = gsub(file_en, "", basename(files[s]))
  
  # Which mutations were called by deepSNV?
  whitelist   = unlist(deepSNV_per_chr_mutations[grep(patient, deepSNV_files)])

  # VAF filter
  one_more_than_5_perc = apply(vafs, 1, function(r) any(r >= min_vaf))
  
  # Remove those with normal namings
  tumour_cols = !(grepl("N[0-9]{1,2}_DW1", colnames(vafs)) | grepl("BC1_DNA1", colnames(vafs)))
  
  # Subset them
  tumour_vafs = vafs[,tumour_cols]
  
  # Also remove sample which has post-hoc been considered normal
  if(patient %in% forced_out) {tumour_vafs = tumour_vafs[,!colnames(tumour_vafs) %in% forced_norm]}
  
  # Capture for n=1 samples
  if(length(colnames(vafs)[tumour_cols])==1) {
    
    tumour_vafs = cbind(tumour_vafs)
    colnames(tumour_vafs) = colnames(vafs)[tumour_cols]
    
  }
  
  if(nrow(vafs)==1) {
    tumour_vafs = rbind(tumour_vafs)
    rownames(tumour_vafs) = rownames(vafs)
  }
  
  # Create large data frame for extracting data
  plot_df = data.frame(sampleID = patient,
                       chr = unlist(lapply(strsplit(rownames(vafs), split = "_"), function(i) i[1])),
                       pos = unlist(lapply(strsplit(rownames(vafs), split = "_"), function(i) i[2])),
                       ref = refs,
                       mut = alts,
                       gene,
                       protein_change,
                       Mappability = maps[,6][match(rownames(vafs), maps[,1])],
                       Cohort_Frequency = as.numeric(cohort_pos_freq[rownames(vafs)]),
                       QUAL = mutation_qual[[s]],
                       deepSNV = rownames(vafs) %in% whitelist,
                       vaf_filt = one_more_than_5_perc,
                       protein_coding, 
                       subclonal = apply(tumour_vafs, 1, function(j) length(which(j >=min_vaf)) < length(j)),
                       check.names = F)
  
  # Apply relevant filtering for SNVs
  plot_df = plot_df[plot_df$Mappability > min_map &
                    plot_df$Cohort_Frequency <= max_cohort_freq &
                    plot_df$QUAL >= min_qual &
                    plot_df$deepSNV &
                    plot_df$vaf_filt,]
  
  return(plot_df)
  
})

# Collapse the data into a single table
dnds_collected = do.call(rbind, dnds_data)
rownames(dnds_collected) = NULL

# Remove patients to exclude!
dnds_collected = dnds_collected[!dnds_collected$sampleID %in% blacklist_patients,]

# Get only the information the dndscv wants
dnds_input = dnds_collected[,1:5]

#################################################################################################################
############################ Run indel section 'indel_processing_output.R' ######################################
#################################################################################################################

# Collect mutation dataframes
dnds_indel_data = lapply(1:length(mutations), function(s) {
  
  # Get the vafs from the list elements
  vafs = variant_allele_freq[[s]]
  
  # Remove sample with an unusual enrichment of synonymous SNPs
  if(s==2) {vafs = vafs[,!colnames(vafs)%in%c("B11_DW1")]}
  
  # Sometimes there were no indels called
  if(nrow(vafs) > 0) {
    
    # Gather other features
    maps = mutation_mappability[[s]]
    muts = mutations[[s]]
    filt = mutation_filters[[s]]
    
    # Get the DNA change
    refs = getREF(mutations[[s]])
    alts = getALT(mutations[[s]])
    
    # Get annotation data
    gene           = unlist(lapply(strsplit(getINFO(muts), split = "[|]"), function(i) paste0(i[4])))
    protein_change = unlist(lapply(strsplit(getINFO(muts), split = "[|]"), function(i) paste0(i[4]," ",i[15]," ",i[16])))
    protein_pos    = unlist(lapply(strsplit(getINFO(muts), split = "[|]"), function(i) paste0(i[4],"_",i[15])))
    protein_coding = unlist(lapply(strsplit(getINFO(muts), split = "[|]"), function(i) i[16]))!=""
    
    # What's the patient name?
    patient     = gsub(file_en, "", basename(files[s]))
    
    # Which mutations were called by deepSNV?
    whitelist   = unlist(deepSNV_per_chr_mutations[grep(patient, deepSNV_files)])
    
    # Vaf filter
    one_more_than_5_perc = apply(vafs, 1, function(r) any(r >= min_vaf))
    
    # Alter the name according to deepSNV style annotation
    white_search_name = unlist(lapply(strsplit(rownames(vafs), "_"), function(i) paste0(i[1],"_",as.numeric(i[2])+1)))
    
    # Remove those with normal naming    
    tumour_cols = !(grepl("N[0-9]{1,2}_DW1", colnames(vafs)) | grepl("BC1_DNA1", colnames(vafs)))
    
    # Extract only those
    tumour_vafs = vafs[,tumour_cols]
    
    # Remove the sample considered to be normal posthoc
    if(patient %in% forced_out) {tumour_vafs = tumour_vafs[,!colnames(tumour_vafs) %in% forced_norm]}
    
    # Capture for n=1 samples
    if(length(colnames(vafs)[tumour_cols])==1) {
      
      tumour_vafs = cbind(tumour_vafs)
      colnames(tumour_vafs) = colnames(vafs)[tumour_cols]
      
    }
    
    if(nrow(vafs)==1) {
      tumour_vafs = rbind(tumour_vafs)
      rownames(tumour_vafs) = rownames(vafs)
    }
    
    # Create large dataframe for analysis
    plot_df = data.frame(sampleID = patient,
                         chr = unlist(lapply(strsplit(rownames(vafs), split = "_"), function(i) i[1])),
                         pos = unlist(lapply(strsplit(rownames(vafs), split = "_"), function(i) i[2])),
                         ref = refs,
                         mut = alts,
                         gene,
                         protein_change,
                         Mappability = maps[,6][match(rownames(vafs), maps[,1])],
                         Cohort_Frequency = as.numeric(cohort_pos_freq[rownames(vafs)]),
                         QUAL = mutation_qual[[s]],
                         deepSNV = white_search_name %in% whitelist,
                         vaf_filt = one_more_than_5_perc,
                         protein_coding, 
                         subclonal = apply(tumour_vafs, 1, function(j) length(which(j >=min_vaf)) < length(j)),
                         check.names = F)
    
    # Apply indel specific filtering
    plot_df = plot_df[plot_df$Mappability >= min_map &
                      plot_df$Cohort_Frequency <= max_cohort_freq &
                      plot_df$QUAL >= min_qual &
                      plot_df$vaf_filt,]
    
  } else {plot_df = NULL}
  
  return(plot_df)
  
})

# Collapse the data into a single table
dnds_indel_collected = do.call(rbind, dnds_indel_data)
rownames(dnds_indel_collected) = NULL

# Remove patients to exclude!
dnds_indel_collected = dnds_indel_collected[!dnds_indel_collected$sampleID %in% blacklist_patients,]

# Make into characters/numeric
dnds_indel_collected$ref = as.character(dnds_indel_collected$ref)
dnds_indel_collected$mut = as.character(dnds_indel_collected$mut)
dnds_indel_collected$pos = as.numeric(as.character(dnds_indel_collected$pos))

# Remove multiple alts because we are also filtering them out
dnds_indel_collected = dnds_indel_collected[unlist(lapply(strsplit(dnds_indel_collected$mut, split = ","), function(i) length(i)==1)),]

# Change annotation style
dnds_indel_collected$ref = unlist(lapply(strsplit(dnds_indel_collected$ref, split = ""), function(i) paste0(i[-1], collapse = "")))
dnds_indel_collected$mut = unlist(lapply(strsplit(dnds_indel_collected$mut, split = ""), function(i) paste0(i[-1], collapse = "")))

dnds_indel_collected$ref = ifelse(dnds_indel_collected$ref=="", "-", dnds_indel_collected$ref)
dnds_indel_collected$mut = ifelse(dnds_indel_collected$mut=="", "-", dnds_indel_collected$mut)

dnds_indel_collected$pos = dnds_indel_collected$pos+1

# Get only the information the dndscv wants
dnds_indel_input = dnds_indel_collected[,1:5]

#################################################################################################################
######################################### Now start dNdS proper #################################################
#################################################################################################################

# Combine SNVs and indels
dnds_collected = rbind(dnds_indel_collected, dnds_collected)
dnds_input = rbind(dnds_indel_input, dnds_input)

# Run dNdScv on all the data to start
res = dndscv(dnds_input, 
             refdb = "refs/RefCDS_human_GRCh38.p12.rda",
             gene_list=genes[!genes %in% c("NEAT1")], outmats=T,
             max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf)

# Here we produce the global dNdS values
ggplot(res$globaldnds) +
  geom_bar( aes(x=name, y=mle), stat="identity", alpha=0.7) +
  ylim(0,6) +
  ylab("dN/dS") +
  xlab("") +
  scale_x_discrete(labels=c("Missense", "Nonsense", "Splice Site", "All", "Truncating")) +
  geom_hline(yintercept = 1) +
  geom_errorbar( aes(x=name, ymin=cilow, ymax=cihigh), width=0.4, colour="grey", alpha=0.9, size=1.3)
ggsave("results/dnds/Global_dNdS_FORECAST.pdf", width = 7, height = 7)

# Save it
saveRDS(res, file = "results/dnds/Global_dNdS_FORECAST.rds")

# We would also like to calculate per gene confidence intervals
ci = geneci(res)

# Get that part of the list and print anything with a qvalue < 0.1
sel_cv = res$sel_cv
print(sel_cv[sel_cv$qallsubs_cv<0.1,c(1:10,ncol(sel_cv))], digits = 3)

# Sort them by missense mle
ci$gene = factor(ci$gene, levels = ci$gene[order(ci$mis_mle, decreasing = T)])

# Plot per gene missense mles
ggplot(ci) +
  geom_bar( aes(x=gene, y=mis_mle), stat="identity", alpha=0.7) +
  ylab("Missense dN/dS") +
  xlab("") +
  geom_hline(yintercept = 1, lty = "dotted") +
  geom_errorbar(aes(x=gene, ymin=mis_low, ymax=mis_high), width=0.4, colour="grey", alpha=0.9, size=1.3) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("results/dnds/Per_gene_missense_dNdS_FORECAST.pdf", width = 7, height = 7)

# Save it
saveRDS(ci, file = "results/dnds/Per_gene_missense_dNdS_FORECAST.rds")

# Sort them by truncating mle
ci$gene = factor(ci$gene, levels = ci$gene[order(ci$tru_mle, decreasing = T)])

# Plot per gene truncating mles
ggplot(ci) +
  geom_bar( aes(x=gene, y=tru_mle), stat="identity", alpha=0.7) +
  ylab("Truncating dN/dS") +
  xlab("") +
  geom_hline(yintercept = 1, lty = "dotted") +
  geom_errorbar(aes(x=gene, ymin=tru_low, ymax=tru_high), width=0.4, colour="grey", alpha=0.9, size=1.3) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("results/dnds/Per_gene_truncating_dNdS_FORECAST.pdf", width = 7, height = 7)

# Save it
saveRDS(ci, file = "results/dnds/Per_gene_truncating_dNdS_FORECAST.rds")

#################################################################################################################
############################################## Subclonal Mutations ##############################################
#################################################################################################################

# Remove those where there is not enough samples to comment on subclonality!
dnds_collected = dnds_collected[!dnds_collected$sampleID %in% exc_patient,]

# Get input of only the mutations that are subclonal
dnds_input = dnds_collected[dnds_collected$subclonal,1:5]

# Run dNdS on the subclonal mutations
res = dndscv(dnds_input, 
             refdb = "refs/RefCDS_human_GRCh38.p12.rda",
             gene_list=genes[!genes %in% c("NEAT1")], outmats=T,
             max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf)

# Plot results of subclonal analysis
ggplot(res$globaldnds) +
  geom_bar( aes(x=name, y=mle), stat="identity", alpha=0.7) +
  ylab("Subclonal dN/dS") +
  xlab("") +
  scale_y_continuous(breaks=seq(0,11,1)) +
  scale_x_discrete(labels=c("Missense", "Nonsense", "Splice Site", "All", "Truncating")) +
  geom_hline(yintercept = 1) +
  geom_errorbar( aes(x=name, ymin=cilow, ymax=cihigh), width=0.4, colour="grey", alpha=0.9, size=1.3)
ggsave("results/dnds/Subclonal_dNdS_FORECAST.pdf", width = 7, height = 7)

# Save it
saveRDS(res, file = "results/dnds/Subclonal_dNdS_FORECAST.rds")

#################################################################################################################
############################################### Clonal Mutations ################################################
#################################################################################################################

# Now take non-subclonal results
dnds_input = dnds_collected[!dnds_collected$subclonal,1:5]

# Run dNdS
res = dndscv(dnds_input, 
             refdb = "refs/RefCDS_human_GRCh38.p12.rda",
             gene_list=genes[!genes %in% c("NEAT1")], outmats=T,
             max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf)

# Plot those results
ggplot(res$globaldnds) +
  geom_bar( aes(x=name, y=mle), stat="identity", alpha=0.7) +
  ylab("Clonal dN/dS") +
  xlab("") +
  scale_y_continuous(breaks=seq(0,11,1)) +
  scale_x_discrete(labels=c("Missense", "Nonsense", "Splice Site", "All", "Truncating")) +
  geom_hline(yintercept = 1) +
  geom_errorbar( aes(x=name, ymin=cilow, ymax=cihigh), width=0.4, colour="grey", alpha=0.9, size=1.3)
ggsave("results/dnds/Clonal_dNdS_FORECAST.pdf", width = 7, height = 7)

# Save it
saveRDS(res, file = "results/dnds/Clonal_dNdS_FORECAST.rds")

#################################################################################################################
############################################### Summary for paper ###############################################
#################################################################################################################

# Publication plot
all = readRDS("results/dnds/Global_dNdS_FORECAST.rds")
clonal = readRDS("results/dnds/Clonal_dNdS_FORECAST.rds")
subclonal = readRDS("results/dnds/Subclonal_dNdS_FORECAST.rds")

# Combine results for ggplot2
plot_df = rbind(all$globaldnds, clonal$globaldnds, subclonal$globaldnds)

# Readable names
plot_df$Type = c(rep("All", times = 5), rep("Clonal", times = 5), rep("Subclonal", times = 5))
plot_df$name = as.character(plot_df$name)

# Get only missense and truncating
plot_df = plot_df[plot_df$name %in% c("wmis", "wtru"),]

# Replace with readable names
plot_df$name[plot_df$name=="wmis"] = "Missense"
plot_df$name[plot_df$name=="wtru"] = "Truncating"

# Which ones do we consider significant?
signif_rows = plot_df[plot_df$cilow > 1,]

# Make a star to annotate those
ann_text = data.frame(Type = signif_rows$Type, mle = signif_rows$cihigh + 1, 
                      label = "*",
                      name = factor(signif_rows$name, levels = unique(plot_df$name)))

# Make and save the plot for the manuscript
ggplot(plot_df, aes(x = Type, y = mle)) +
  geom_errorbar( aes(x=Type, ymin=cilow, ymax=cihigh), width=0.4, colour="grey", alpha=0.9, size=1.3) + 
  geom_point() + geom_hline(yintercept = 1) + xlab("") + ylab("dN/dS") + facet_grid(. ~ name) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  geom_text(data = ann_text, mapping = aes(x = Type, y = mle, label = label), size = 6)
ggsave("results/dnds/Paper_plot.png", width = 4, height = 3)
