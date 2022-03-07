# A script for making the heatmap for the patients after running heatmap script prior

# Where are our GitHub repos?
low_pass_repo = "low_pass_repo/"
targeted_repo = "targeted_repo/"

# First go to low pass analysis repo
setwd(low_pass_repo)

# Source general scripts
source("scripts/00_general_functions.R")

# Read in the outcome data
outcome = read.csv("refs/FORECAST_outcome_data_v2.csv", stringsAsFactors = F)

# Read in the data characteristics 
characteristics = read.csv("refs/FORECAST_patient_characteristics.csv", stringsAsFactors = F)

# Read in the rescoring pathologist scores
gleason_rescore = read.csv("refs/GleasonRescore_patients.csv", 
                           header = T, stringsAsFactors = F)

# Read in the metrics collected
metrics = read.table("results/7_collect_metrics/FORECAST_genomic_metrics.txt", sep = "\t", header = T)

# Read in GM data too
gleason_morisita = read.csv("refs/GleasonMorisita_Patients.csv", 
                            header = F, stringsAsFactors = F)

# Name GM columns
colnames(gleason_morisita) = c("Patient", "Gleason_Morisita")

# Make columns more readable
colnames(outcome) = c("Patient", "Recurrence", "Recurrence_Months", "Death", "Death_Months", "Mets", 
                      "Mets_Months", "Wide_Mets", "Wide_Met_Months")

# Rename the characteristics column
colnames(characteristics)[1] = "Patient"

# Add IM to beginning 
outcome$Patient         = ifelse(grepl("DNT", outcome$Patient), 
                                 yes = outcome$Patient, no = paste0("IM",outcome$Patient))
characteristics$Patient = ifelse(grepl("DNT", characteristics$Patient), 
                                 yes = characteristics$Patient, no = paste0("IM",characteristics$Patient))

# Merge by patient id
collected = merge(outcome, characteristics, by = "Patient")
collected = merge(collected, gleason_rescore[,c("Patient", "Gleason_Grade_Rescore")], 
                  by = "Patient", all.x = T)
collected = merge(collected, metrics, by = "Patient")

# Make into categories
collected$High_PSA               = collected$Presenting_PSA > 20
collected$Is_T3                  = collected$T_stage == "T3"
not_rescored                     = is.na(collected$Gleason_Grade_Rescore)
collected$Gleason_Grade_Rescore[not_rescored] = collected$Gleason_grade_group[not_rescored]

# Format grade groups
collected$Gleason_grade_group    = as.factor(collected$Gleason_Grade_Rescore)

# Return to targeted panel repo
setwd(targeted_repo)

# Get data in column order
col_anno = collected[match(col_order, collected$Patient),c("Patient", "Prefilter_samples", "mPGA", 
                                                           "Gleason_Grade_Rescore", "Presenting_PSA", 
                                                           "T_stage", "Recurrence", "Mets", "Death", 
                                                           "Recurrence_Months")]

# Make the heatmap annotations
column_ha = HeatmapAnnotation(Samples = col_anno$Prefilter_samples,
                              Gleason = col_anno$Gleason_Grade_Rescore,
                              Tstage  = col_anno$T_stage,
                              Outcome = col_anno$Recurrence | col_anno$Death,
                              col = list(Samples = colorRamp2(c(1, 2, 3, 10), c("red", "red", "white", "#0343df")),
                                         Gleason = c("2" = "#F0E442", "3" = "#D55E00", "4" = "#5D3A9B", "5" = "#000000"),
                                         Tstage = c("T1" = "#e6daa6", "T2" = "#8e82fe", "T3" = "#ffb16d"),
                                         Outcome = c("FALSE" = "#3f9b0b", "TRUE" = "#06470c")#,
                                         ),
                              show_legend = c(T,T,T,T,T,T,T),
                              annotation_legend_param = list(Outcome = list(title = "Recurrence/Death"),
                                                             Tstage = list(title="T-stage"),
                                                             Samples = list(at = c(1, 3, 10), labels = c("1", "3", "10+"))))

# Remove NEAT1 because we basically can't call coding mutation there
res_plot = res_plot[!grepl("NEAT1", rownames(res_plot)),]

# We use this multiple times now
col_fun = colorRamp2(c(-3, -2, -1, 0, 1), c("#601A4A", "white", "#63ACBE", "white", "#EE442F"))

# Make heatmap
pdf("results/mutation_calling/Patient_summary_heatmap_paper.pdf", width = 16, height = 7)
hm = Heatmap(res_plot[,col_order], cluster_rows = F, cluster_columns = F, rect_gp = gpar(col = "white", lwd = 2),
             top_annotation = column_ha,
             column_title = "Patient",
             column_title_side = "bottom",
             row_title = "Gene",
             column_names_gp = grid::gpar(fontsize = 8),
             col = col_fun, 
             border = T, show_heatmap_legend = F)
lg = Legend(at = c(1,-1,-3), title = "Mutation", 
            legend_gp = gpar(fill = col_fun(c(1,-1,-3))), 
            labels = c("SNV", "INDEL", "SNV/INDEL"))
lg2 = Legend(at = c(1,0.5), title = "Clonality", 
             legend_gp = gpar(fill = col_fun(c(1,0.5))), 
             labels = c("Clonal", "Subclonal"))
lgs = packLegend(lg, lg2)
draw(hm)
draw(lgs, x = unit(0.95, "npc"), y = unit(0.2, "npc"))
dev.off()
