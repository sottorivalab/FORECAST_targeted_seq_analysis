# A script to create the supplementary table
ind = read.table("results/mutation_calling/Per_sample_gene_mutation_status_with_change_indels.txt", 
                 stringsAsFactors = F, header = T, sep = "\t")

snv = read.table("results/mutation_calling/Per_sample_gene_mutation_status_with_change.txt", 
                 stringsAsFactors = F, header = T, sep = "\t")

# Merge and write out
df = rbind(snv, ind)

# Delete columns
df$deepSNV  = NULL
df$vaf_filt = NULL
df$protein_coding = NULL

# So the numbers aren't so long
df$VAF = signif(df$VAF, digits = 3)

# Reorder
df = df[,c("Patient", "Sample", "Gene", "Protein", "VAF")]

# Write it out
write.csv(df, file = "results/mutation_calling/Per_sample_gene_mutation_status_with_change_publication.csv", 
          quote = F, row.names = F)