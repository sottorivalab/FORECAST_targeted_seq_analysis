# Speed up finding files in the TES pipeline
findFilesTES = function(file_pattern, sub_dir, pipeline_root, patients = "auto", recursive_in_subdir = T) {
  
  if(patients == "auto") {
  
    dirs     = list.files(pipeline_root)
    
    patients = c(dirs[grep("^DNT", dirs)], dirs[grep("^IM", dirs)])
    
    print(patients)
    
  }
  
  per_patient = lapply(patients, function(p) {
    
    list.files(paste0(pipeline_root,p,"/",sub_dir), 
               pattern = file_pattern, recursive = recursive_in_subdir, full.names = T)
  
  })
  
  files = unlist(per_patient)
  
  return(files)
  
}

# Calculate mean from frequency table
meanFromFreqTable = function(df, value_col = 1, count_col = 2) {
  
  sum(apply(df, 1, function(i) i[value_col]*i[count_col])) / sum(df[,count_col])
  
}

# Extract patient name
extractPatientName = function(x) {unlist(lapply(strsplit(x, split = "_"), function(i) i[1]))}

# Genes in this study that are on the panel
genes = c("AKT1",
          "CDKN1B",
          "PTEN",
          "APC",
          "CDK12",
          "RB1",
          "AR",
          "CTNNB1",
          "SPOP",
          "ATM",
          "FOXA1",
          "TP53",
          "ARID1A",
          "ARID1B",
          "IDH1",
          "ZFHX3",
          "ASXL1",
          "KDM6A",
          "NEAT1",
          "BRAF",
          "KMT2C",
          "KMT2D",
          "CHD1",
          "BRCA1",
          "BRCA2",
          "PIK3CA",
          "PIK3R1",
          "PALB2")

# Calculate number of mutant copies
calculateCopies = function(V, C, r, exp_norm_cn = 2) {
  # Assumes CCF is 1
  (V * (((C - exp_norm_cn)*r)+exp_norm_cn)) / r
}

# Get simple clonal LOH assessment
isitLOH = function(V, C, r, exp_norm_cn = 2) {
  copies = calculateCopies(V = V, C = C, r = r, exp_norm_cn = exp_norm_cn)
  copies > (C - 0.5)
}

# Make location dataframe from chr pos rownames
extractChrRows = function(x) {
  m = do.call(rbind, strsplit(x, split = "[:]|-"))
  df = data.frame(chr = m[,1],
                  start = as.integer(m[,2]),
                  end = as.integer(m[,3]))
  df$chr = factor(df$chr, levels = c(1:22, "X", "Y"))
  return(df)
}
