# Pipeline Overview: Single-cell MR Discovery, Replication, and Colocalization

# ----------------------------
# Step 1: Discovery Mendelian Randomization (MR)
# ----------------------------

## 1.1 Load necessary libraries and data
library(readxl)
library(dplyr)
library(TwoSampleMR)
library(MRInstruments)
library(writexl)

# 1.1 Load raw single-cell eQTL data and format
dataset1K1K_raw <- read_excel("file_path/dataset1K1K.xlsx")
dataset1K1K <- dataset1K1K_raw %>%
  dplyr::filter(pvalue < 1e-5) %>%
  dplyr::mutate(
    se            = abs(`rho correlation coefficient` / qnorm(pvalue/2)),
    beta          = `rho correlation coefficient`,
    pval          = pvalue,
    cell          = `Cell type`,     
    phenotype     = `Gene ID`,
    effect_allele = `SNP assessed allele`
  ) %>%
  dplyr::select(cell, phenotype, SNP, effect_allele, beta, se, pval)

# Save formatted dataset1K1K file for reuse
write_xlsx(dataset1K1K, "file_path/dataset1K1K_formatted.xlsx")

## 1.2 Load and format outcome GWAS data
outcome_df <- data.table::fread(file = "file_path/outcome.tsv", header = TRUE)
outcome_df <- as.data.frame(outcome_df)
outcome <- format_data(
  outcome_df,
  type                 = "outcome",
  snp_col              = "SNP",
  beta_col             = "beta",
  se_col               = "standard_error",
  effect_allele_col    = "effect_allele",
  other_allele_col     = "other_allele",
  pval_col             = "p_value"
)

# Save formatted outcome for reuse
save(outcome, file = "file_path/outcome.RData")


## 1.3 Run MR for a cell type
cell_type <- "B IN"   # manual: set cell_type (`"B IN", "B Mem", "CD4 ET", 
#"CD4 NC", "CD4 SOX4", "CD8 ET", "CD8 NC", "CD8 S100B", "DC", "Mono C", "Mono NC", 
#"NK", "NK R", "Plasma")

exp_df <- dataset1K1K %>%
  dplyr::filter(cell == cell_type)

exp_obj <- format_data(
  exp_df,
  type              = "exposure",
  snp_col           = "SNP",
  beta_col          = "beta",
  se_col            = "se",
  effect_allele_col = "effect_allele",
  phenotype_col     = "phenotype",
  pval_col          = "pval"
)

har <- harmonise_data(exp_obj, outcome, action = 1)
mr_res_cell_type <- mr(har) %>%
  dplyr::filter(method %in% c("Inverse variance weighted", "Wald ratio")) %>%
  dplyr::mutate(qval = p.adjust(pval, method = "BH")) %>%
  mutate(cell = cell_type)

# (Repeat the above block for each cell type)

# After running for all desired cell types, bind and save combined results:
combined <- dplyr::bind_rows(mget(ls(pattern = "^mr_res_")))
write_xlsx(combined, "file_path/discMR_outcome.xlsx")

# ----------------------------
# Step 2: Replication Mendelian Randomization (MR)
# ----------------------------

## 2.1 Load previously formatted outcome for replication
load("file_path/outcome.RData")

## 2.1 Load libraries and replication exposure dataset
library(readr)
library(readxl)
library(dplyr)
library(meta)
library(genetics.binaRies)
library(TwoSampleMR)

scBloodNL_cell_type <- read_delim(
  "file_path/cell_type_expression_eQTLsFDR-ProbeLevel.txt", delim = "\t") %>% 
  dplyr::filter(PValue < 1e-5) # manual: set cell_type ("B", "CD4", "CD8", "DC", 
#"Mono", "NK")

## 2.2 Identify significant genes from discovery MR
discMR <- read_excel("file_path/discMR_outcome.xlsx")
sig_genes <- discMR %>% filter(cell %in% cell_types, qval < 0.05) %>%
  pull(exposure) # manual: set cell_types ("B IN", "B Mem", "Plasma" for "B", 
#"CD4 ET", "CD4 NC", "CD4 SOX4" for "CD4", "CD8 ET", "CD8 NC", "CD8 S100B" for "CD8", 
#"DC" for "DC", "Mono C", "Mono NC" for "Mono", "NK", "NK R" for "NK")

## 2.3 Subset replication exposure and export for manual formatting
subset_cell_type <- scBloodNL_cell_type %>% dplyr::filter(HGNCName %in% sig_genes)
write_xlsx(subset_cell_type, "file_path/subset_cell_type.xlsx")
# manual: keep required columns, format column names and add beta_i, se_i, beta_ii, se_ii and save

rep_cell_type <- read_excel("file_path/subset_cell_type.xlsx")

df1 <- rep_cell_type %>% dplyr::select(HGNCName, SNP, beta_i, se_i) %>% 
  dplyr::mutate(obs = row_number())
df2 <- rep_cell_type %>% dplyr::select(HGNCName, SNP, beta_ii, se_ii) %>%
  rename("beta_i" = "beta_ii", "se_i" = "se_ii") %>% dplyr::mutate(obs = row_number())
meta_df <- bind_rows(df1, df2) %>% arrange(obs)

ma_rep <- metagen(TE = beta_i, seTE = se_i, data = meta_df, common = TRUE, 
                  random = FALSE, sm = "MD", backtransf = FALSE, subgroup = obs)
rep_cell_type <- bind_cols(rep_cell_type,
  data.frame(beta = ma_rep$TE.fixed.w, se = ma_rep$seTE.fixed.w))

clumped <- ld_clump_local(
  tibble(rsid = rep_cell_type$SNP, pval = rep_cell_type$PValue),
  plink_bin = "file_path/plink_win64_20231211/plink.exe",
  bfile     = "file_path/Eur_ldref/EUR",
  clump_r2  = 0.1, clump_p = 1, clump_kb = 10000
)
rep_cell_type <- rep_cell_type %>% dplyr::filter(SNP %in% clumped$rsid)

## 2.4 Format and run replication MR
rep_cell_type <- format_data(
  rep_cell_type,
  type              = "exposure",
  snp_col           = "SNP",
  beta_col          = "beta",
  se_col            = "se",
  effect_allele_col = "effect_allele",
  phenotype_col     = "HGNCName",
  pval_col          = "PValue"
)
har_rep <- harmonise_data(rep_cell_type, outcome, action = 1)
mr_cell_type <- mr(har_rep) %>%
  dplyr::filter(method %in% c("Inverse variance weighted", "Wald ratio")) %>%
  dplyr::mutate(qval = p.adjust(pval, method = "BH")) %>%
  mutate(cell = cell_type)

# (Repeat the above block for each cell type)

# After running for all desired cell types, bind and save combined results:
combined <- dplyr::bind_rows(mget(ls(pattern = "^mr_")))
write_xlsx(combined, "file_path/repMR_outcome.xlsx")

# ----------------------------
# Step 3: Colocalization Analysis
# ----------------------------

## 3.1 Format outcome for colocalisation
outcome_df <- data.table::fread(file = "file_path/outcome.tsv", header = TRUE) %>%
  dplyr::rename("snp" = "SNP", "se" = "standard_error", "pval" = "p_value") %>%
  dplyr::mutate(varbeta = se^2) %>%
  dplyr::select(snp, beta, varbeta, pval)
outcome <- as.list(outcome_df)
outcome$type <- "cc"
outcome$s    <- p #manual: add proportion of cases in outcome population

## 3.2 Subset eQTLs for finer cell and matching genes
disc_hits    <- read_excel("file_path/discMR_outcome.xlsx") %>%
  filter(cell == cell_type, qval < 0.05) # manual: set cell_type (`"B IN", "B Mem", 
#"CD4 ET", "CD4 NC", "CD4 SOX4", "CD8 ET", "CD8 NC", "CD8 S100B", "DC", "Mono C", 
#"Mono NC", "NK", "NK R", "Plasma")
rep_hits     <- read_excel("file_path/repMR_outcome.xlsx") %>% 
  filter(cell == cell_type, qval < 0.05) # manual: set cell_types 
#("B" for "B IN","B Mem","Plasma","CD4" for "CD4 ET","CD4 NC","CD4 SOX4",
#  "CD8" for "CD8 ET","CD8 NC","CD8 S100B","DC" for "DC",
#  "Mono" for "Mono C","Mono NC","NK" for "NK","NK R")
matching_genes <- intersect(disc_hits$exposure, rep_hits$exposure)

dataset1K1K <- read_excel("file_path/dataset1K1K_formatted.xlsx", col_names = TRUE)
eqtl_df <- dataset1K1K %>% filter(cell == cell_type, phenotype %in% matching_genes)
# manual: set cell_type (`"B IN", "B Mem", #"CD4 ET", "CD4 NC", "CD4 SOX4", "CD8 ET", 
#"CD8 NC", "CD8 S100B", "DC", "Mono C", "Mono NC", "NK", "NK R", "Plasma")

eqtl_df$varbeta <- eqtl_df$se^2
eqtl_df$MAF <- MAF #MAF values of SNPs can be obtained from Ensembl
eqtl_df <- eqtl_df %>%
  rename("gene" = "phenotype", "snp" = "SNP", "position" = "Position")


## 3.3a Colocalization using coloc.abf
library(coloc)

# Get unique gene values
unique_genes <- unique(eqtl_df$gene)

# Initialize a list to store colocalization results
results_list <- list()

# Initialize an empty dataframe to store the posterior probabilities
pp_df_abf <- data.frame(
  Gene     = character(),
  PP_H0    = numeric(),
  PP_H1    = numeric(),
  PP_H2    = numeric(),
  PP_H3    = numeric(),
  PP_H4    = numeric(),
  stringsAsFactors = FALSE
)

for (gene in unique_genes) {
  # Subset dataframe for the current gene and convert to list
  gene_subset <- as.list(eqtl_df[eqtl_df$gene == gene, ])
  
  # Add scalars to the list
  gene_subset$type <- "quant"
  gene_subset$N    <- s  # manual: set sample size (CD4 NC=463528, 
  #CD4 ET=61786, CD4 SOX4=4065, CD8 ET=205077, CD8 NC=133482, CD8 S100B=34528, 
  #DC=8690, Plasma=3625, Mono C=38233, Mono NC=15166, B Mem=48023, B IN=82068, 
  #NK=159820, NK R=9677)
  
  # Run colocalization analysis using coloc.abf()
  coloc_result <- coloc.abf(gene_subset, outcome, MAF = NULL,
                            p1 = 1e-04, p2 = 1e-04, p12 = 1e-05)
  
  # Store the coloc_result for that gene in results_list
  results_list[[gene]] <- coloc_result
  
  # Check if the coloc_result object has row names
  if (length(rownames(coloc_result)) == 0) {
    print(paste("No colocalization results found for gene:", gene))
  } else {
    # Extract the posterior probabilities from the coloc_result object
    pp_h0 <- coloc_result$PP.H0.abf
    pp_h1 <- coloc_result$PP.H1.abf
    pp_h2 <- coloc_result$PP.H2.abf
    pp_h3 <- coloc_result$PP.H3.abf
    pp_h4 <- coloc_result$PP.H4.abf
    
    # Create a row for the current gene in the pp_df_abf dataframe
    gene_row <- data.frame(
      Gene  = gene,
      PP_H0 = pp_h0,
      PP_H1 = pp_h1,
      PP_H2 = pp_h2,
      PP_H3 = pp_h3,
      PP_H4 = pp_h4,
      stringsAsFactors = FALSE
    )
    
    # Append the gene_row to the pp_df_abf dataframe
    pp_df_abf <- rbind(pp_df_abf, gene_row)
  }
}

pp_df_abf$cell <- cell_type
pp_df_abf_cell_type <- pp_df_abf

# After running for all desired cell types, bind and save combined results:
combined <- dplyr::bind_rows(mget(ls(pattern = "^pp_df_abf_")))
write_xlsx(combined, "file_path/colabf_outcome.xlsx")


## 3.3b Colocalization using coloc.susie for Gene IDs with >1 eQTL as IV
remotes::install_github("chr1swallace/coloc")

library(coloc)
library(data.table)
library(TwoSampleMR)

# Identify overlapping SNPs with outcome
eqtl_df <- eqtl_df %>% rename(SNP = snp)
outcome <- data.table::fread(file = "file_path/outcome.tsv", header = TRUE)
common_snps <- intersect(eqtl_df$SNP, outcome$SNP)
outcome <- outcome[outcome$SNP %in% common_snps, ]

# Mismatch correction
# Identify allele mismatches and correct them
mismatched_snps <- eqtl_df %>%
  inner_join(outcome, by = "SNP") %>%
  filter(effect_allele.x != effect_allele.y)

if (nrow(mismatched_snps) > 0) {
  for (i in 1:nrow(mismatched_snps)) {
    snp <- mismatched_snps$SNP[i]
    
    # Update effect allele to match outcome and flip beta
    eqtl_df <- eqtl_df %>%
      mutate(
        effect_allele = ifelse(SNP == snp, outcome$effect_allele[outcome$SNP == snp], effect_allele),
        beta = ifelse(SNP == snp, -beta, beta)
      )
    
    cat("Corrected allele mismatch for SNP:", snp, "\n")
  }
}

# Write the list of common SNPs
write.table(common_snps, "locus_snps.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
system("file_path/plink_win64_20231211/plink.exe --bfile file_path/Eur_ldref/EUR --extract locus_snps.txt --keep-allele-order --out locus_ld --r square --write-snplist")

# Read LD matrix and SNP names
ld <- read.table("locus_ld.ld", header = FALSE)
ldnames <- fread("locus_ld.snplist", header = FALSE)

# Create LD matrix
ld.mat <- matrix(as.vector(data.matrix(ld)), nrow = nrow(ldnames), ncol = nrow(ldnames))

# Name rows and columns
colnames(ld.mat) <- ldnames$V1
rownames(ld.mat) <- ldnames$V1

# Prepare Dataset1 from eqtl_df
Dataset1 <- list(
  beta = eqtl_df$beta,
  varbeta = eqtl_df$se^2,
  snp = eqtl_df$SNP,
  type = "quant",
  MAF = eqtl_df$MAF,
  N = s, # manual: set sample size (CD4 NC=463528, 
  #CD4 ET=61786, CD4 SOX4=4065, CD8 ET=205077, CD8 NC=133482, CD8 S100B=34528, 
  #DC=8690, Plasma=3625, Mono C=38233, Mono NC=15166, B Mem=48023, B IN=82068, 
  #NK=159820, NK R=9677)
  LD = ld.mat
)

# Prepare Dataset2 from outcome_sub
Dataset2 <- list(
  beta = outcome$beta,
  varbeta = outcome$se^2,
  snp = outcome$SNP,
  type = "cc", 
  N = s, # manual: set sample size
  LD = ld.mat
)

# Run SuSiE fine-mapping on both traits
S1 <- runsusie(Dataset1)
S2 <- runsusie(Dataset2)

# Perform coloc.susie
susie_res <- coloc.susie(S1, S2)

# Extract rsid and PP.H4.abf
susie_df <- data.frame(
  rsid      = susie_res$summary$hit1,
  PP_H4_abf = susie_res$summary$`PP.H4.abf`,
  stringsAsFactors = FALSE
)

susie_df$cell <- cell_type
susie_df_cell_type <- susie_df

# After running for all desired cell types, bind and save combined results:
combined <- dplyr::bind_rows(mget(ls(pattern = "^susie_df_")))
write_xlsx(combined, "file_path/colsus_outcome.xlsx")

