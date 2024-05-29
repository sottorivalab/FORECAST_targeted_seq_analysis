library(vcfR)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(plyr)
library(reshape2)
library(ComplexHeatmap)
library(circlize)

# Function for expected VAF in tetraploid tumour
vaf_2_4 = function(r) (2*r) / (2*r+2)

# Parameters
indel_map    = 0.88
snv_map      = 0.6
indel_qual   = 100
snv_qual     = 10
max_rec      = 1
max_norm_vaf = 0

# Read in the driver gene list
driver_genes = read.csv("refs/Driver_genes.csv")
driver_genes = driver_genes[!duplicated(driver_genes$Gene),]

# Which patients will we analyse?
pats = c("")

# Colours for impact type
cols = c("#15b01a", "#0485d1", "#e50000")
names(cols) = c("LOW", "MODERATE", "HIGH")

# Colours for cohort 
cohort_cols = c("#bf77f6", "#95d0fc", "#d8dcd6")
names(cohort_cols) = c("FORECAST", "IMMUNE", "PANCAN")

# Read in the vcf files
muts = lapply(pats, function(pat) {
  
  # Read it in
  muts = read.vcfR(paste0("~/remotes/forecast_tes/nf-core-mutectplatypus/results/Processed/call_variants/platypus_min_reads_3/",
                          pat,"/",pat,"_platypus.ann.vcf.gz"))
  
})

# Create data frames for plotting and filtering
dfs = lapply(muts, function(mut) {
  
  # Collect all necessary metrics
  df = data.frame(chr = getCHROM(mut),
                  pos = getPOS(mut),
                  ref = getREF(mut),
                  alt = getALT(mut),
                  qual = getQUAL(mut),
                  do.call(rbind, lapply(getINFO(mut), function(i) unlist(strsplit(i, split = "[|]"))[c(4:2,16)])),
                  signif(extract.gt(mut, element = "NV", as.numeric = T) / extract.gt(mut, element = "NR", 
                                                                                      as.numeric = T), digits = 3),
                  stringsAsFactors = F)
  
  # Rename the columns
  colnames(df)[6:9] = c("gene", "Impact", "type", "protein_change")
  
  # Make the rows equal the mutation id
  df$id = rownames(df)
  
  return(df)

})

# Name for patient identification
names(dfs) = pats

# Count number of observations of each mutation
multimuts = table(unlist(lapply(dfs, rownames)))

# Read in the mappability score
maps = lapply(pats, function(pat) {
  
  # Read as a dataframe
  df = read.table(paste0("~/remotes/forecast_tes/nf-core-mutectplatypus/results/Reports/mappability/",
                         pat,"/",pat,"_mappability.out"), 
                  stringsAsFactors = F)
  
  return(df)
  
})

# Name for patient identification
names(maps) = pats

# Tag the filter passing result for mappability, indel, quality, if it is multimutated, protein change and multialt
dfs = lapply(pats, function(pat) {
  
  # Match the mappability
  maps = maps[[pat]]$V6[match(dfs[[pat]]$id, maps[[pat]]$V1)]
  
  # Get the dataframe from the list
  out = dfs[[pat]]
  
  # What's the associated mappability
  out$mappability = maps
  
  # If replacement length differs, call it a indel
  indel = unlist(lapply(strsplit(out$ref, split = ""), length)) != unlist(lapply(strsplit(out$alt, split = ""), 
                                                                                 length))
  
  # Record this
  out$indel = indel
  
  # Does it pass the different mappability filters depending on mutation type?
  out$map_filter = ifelse(out$indel, out$mappability >= indel_map, out$mappability >= snv_map)
  
  # Does it pass the different quality filters depending on mutation type?
  out$qual_filter = ifelse(out$indel, out$qual >= indel_qual, out$qual >= snv_qual)
  
  # Does it appear multiple times
  out$multimut    = multimuts[match(dfs[[pat]]$id, names(multimuts))]
  
  # Does it change the protein sequence?
  out$protein_change = grepl("/", dfs[[pat]]$protein_change)
  
  # Is there just one alternative?
  out$not_multi_alt  = !grepl(",", dfs[[pat]]$alt)
  
  return(out)
  
})

# Make a dataframe version with filters applied
dfs_filt = lapply(dfs, function(i) i[i$map_filter & i$qual_filter & i$multimut<=max_rec & 
                                       (i$protein_change | grepl("synonymous",i$type)) & i$not_multi_alt & 
                                       i$BC1_DNA1 <= max_norm_vaf,])

# Save this for other uses
saveRDS(dfs_filt, file = "results/cfDNA/FORECAST_cfDNA_mutations_filtered_nonsyn_and_syn.rds")

# Melt it to play ball with ggplot 
plt_df = lapply(dfs_filt, function(i) {
  
  i$indel          = NULL
  i$map_filter     = NULL
  i$qual_filter    = NULL
  i$multimut       = NULL
  i$protein_change = NULL
  i$not_multi_alt  = NULL
  
  res = melt(i[,9:(ncol(i)-1)], id = "id")
  res$value  = as.numeric(res$value)
  res$type   = i$type[match(res$id, i$id)]
  res$Impact = i$Impact[match(res$id, i$id)]
  res$map    = i$mappability[match(res$id, i$id)]
  res$qual   = i$qual[match(res$id, i$id)]
  res$gene   = i$gene[match(res$id, i$id)]
  res$type[grep("&", res$type)] = "mixed"
  return(res)
  
})

#################################################################################################################

# What overlaps with driver mutation?
key_genes = lapply(dfs_filt, function(i) i[i$gene %in% driver_genes$Gene,])

# Plot out the vaf distributions
plt_list = list()
loop = 0
for(i in c(1,2,3,"etc")) {
      
  loop = loop+1

  # Rearrange order, update manually
  if(i == "replace_with_number") {plt_df[[i]]$variable = factor(plt_df[[i]]$variable, levels = levels(plt_df[[i]]$variable)[c(1,2,3)])}

  # Grab the driver genes with some VAF
  driver_gene_pat = plt_df[[i]][plt_df[[i]]$gene %in% driver_genes$Gene & plt_df[[i]]$variable!="BC1_DNA1" & plt_df[[i]]$type!="synonymous_variant",]
  driver_gene_pat$list = driver_genes[match(driver_gene_pat$gene, driver_genes$Gene),2]
  driver_gene_pat = driver_gene_pat[driver_gene_pat$value!=0,]
  
  driver_gene_pat = driver_gene_pat[driver_gene_pat$list!="IMMUNE",]
    
  if(i=="code_to_add_legend_if_needed") {lg_pos = c(0.1,0.7)} else {lg_pos = "none"}
  
  # Plot it out
  p = ggplot(plt_df[[i]], aes(x = variable, y = value, group = id, col = Impact)) + geom_line(alpha = 0.1) + 
    geom_point() + 
    scale_colour_manual(values = cols) +
    annotate("label", 
             x = driver_gene_pat$variable,
             y = driver_gene_pat$value,
             label = driver_gene_pat$gene, 
             fill = cohort_cols[driver_gene_pat$list], cex = 2) +
    ylab("VAF") + xlab("") + 
    ggtitle(pats[i]) + 
    scale_y_continuous(trans='log10') + 
    theme_cowplot() + theme(plot.title = element_text(hjust = 0.5), strip.text.x = element_blank(),
                            strip.background = element_rect(colour="white", fill="white"),
                            legend.position=lg_pos)
  plt_list[[loop]] = p
  
}
ggarrange(plt_list[[1]], plt_list[[2]], plt_list[[3]], 
          plt_list[[4]], plt_list[[5]],
          ncol = 1, nrow = 5)
ggsave("results/cfDNA/FORECAST_cfDNA_vaf_distributions_logscale.png", width = 8.27, height = 11.69)
