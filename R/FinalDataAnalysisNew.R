################################################################################
# COMPREHENSIVE PARASITE PREVALENCE AND RISK FACTOR ANALYSIS
# Project: Zoonotic Endoparasites in Free-Roaming Cats (NYC)
#
# OVERVIEW:
# This script processes the raw data ('dataset1.csv') to generate all the 
# statistical results reported in the manuscript, including:
# 1. Prevalence calculations with 95% Confidence Intervals (Table 1).
# 2. Risk factor analysis (Fisher's Exact Tests) for infection presence.
# 3. Infection intensity analysis (Negative Binomial GLM / Mann-Whitney U) 
#    for egg shedding (EPG).
# 4. Generation of Table 1 and Table 2 as CSV files.
#
# INPUT: "dataset1.csv"
# OUTPUT: Console output and CSV files for Tables 1 & 2.
################################################################################

# ==============================================================================
# 1. SETUP AND LIBRARIES
# ==============================================================================
# We use standard libraries to ensure robust statistical methods.
# MASS: Required for the Negative Binomial GLM (for overdispersed egg counts).

required_packages <- c("MASS", "here", "dplyr", "broom", "gt")

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

output_tables  <- here("output", "tables")
output_figures <- here("output", "figures")

if (!dir.exists(output_tables)) {
  dir.create(output_tables, recursive = TRUE)
}

if (!dir.exists(output_figures)) {
  dir.create(output_figures, recursive = TRUE)
}

# ==============================================================================
# 2. DATA LOADING AND PREPARATION
# ==============================================================================
# Load the dataset. We treat empty strings "", "NA", and single spaces " " as missing data.
df <- read.csv(
  here("Data", "clean", "dataset1.csv"),
  na.strings = c("", "NA", " ")
)

# ------------------------------------------------------------------------------
# 2.1 DEFINE INFECTION STATUS (BINARY OUTCOMES)
# ------------------------------------------------------------------------------
# Goal: Create 0/1 variables for prevalence analysis.
# Logic: If an EPG count exists and is > 0, the cat is Positive (1).
#        If the EPG count is NA (missing), we assume it is Negative (0) 
#        because these columns usually only contain data for positive cases.
#        (Note: The denominators n=87 etc. handle the "Total Tested" logic later).

df$rounds_pos <- ifelse(!is.na(df$roundsepg) & df$roundsepg > 0, 1, 0)
df$hooks_pos  <- ifelse(!is.na(df$hooksepg) & df$hooksepg > 0, 1, 0)
df$coccid_pos <- ifelse(!is.na(df$coccidepg) & df$coccidepg > 0, 1, 0)

# "Any Parasite": Flag if a cat has AT LEAST ONE of the flotation targets (Helminths or Coccidia).
df$any_parasite <- ifelse(df$rounds_pos == 1 | df$hooks_pos == 1 | df$coccid_pos == 1, 1, 0)

# Co-infections: Calculate how many concurrent infections each cat has.
df$coinfection_count <- df$rounds_pos + df$hooks_pos + df$coccid_pos
# Define "Co-infected" as having 2 or more parasites.
df$is_coinfected <- ifelse(df$coinfection_count >= 2, 1, 0)
# Define the specific "Roundworm + Hookworm" coinfection (most common).
df$toxo_ancy_coinfection <- ifelse(df$rounds_pos == 1 & df$hooks_pos == 1, 1, 0)

# Antigen & Serology Outcomes:
# These columns contain "Y" for Positive. Convert to 1 (Pos) / 0 (Neg).
df$giardia_pos <- ifelse(df$giardia == "Y", 1, 0)
df$crypto_pos  <- ifelse(df$crypto == "Y", 1, 0)
df$toxo_Ab_pos <- ifelse(df$toxo == "Y", 1, 0) # Serology antibody positive

# ------------------------------------------------------------------------------
# 2.2 CLEAN DEMOGRAPHICS (RISK FACTORS)
# ------------------------------------------------------------------------------
# Age Grouping:
# Your study defines "Young" as <1 year (Kittens 'K' + Juveniles 'J').
# "Adult" is >=1 year (Adults 'A').
df$Age_Group <- NA # Initialize
df$Age_Group[df$Age %in% c("K", "J")] <- "Young"
df$Age_Group[df$Age == "A"] <- "Adult"

# Sex:
# Standardize values to "Male" and "Female" for clean table output.
df$Sex_Clean <- NA
df$Sex_Clean[df$Sex == "M"] <- "Male"
df$Sex_Clean[df$Sex == "F"] <- "Female"

# Set Reference Levels for Statistical Modeling:
# We set "Adult" and "Female" as the baseline (Reference) levels.
# This means Odds Ratios will represent risk relative to Adult Females.
df$Age_Group <- factor(df$Age_Group, levels = c("Adult", "Young"))
df$Sex_Clean <- factor(df$Sex_Clean, levels = c("Female", "Male"))

# ==============================================================================
# 3. ANALYSIS PART 1: PREVALENCE (TABLE 1 GENERATION)
# ==============================================================================
# Define the exact denominators (Total Tested) for each method as per your Methods section.
n_float   <- 87  # Fecal flotation total
n_antigen <- 43  # Antigen ELISA total
n_sero    <- 45  # Serology total

# ------------------------------------------------------------------------------
# FUNCTION: calc_prev_row
# Purpose: Calculates Prevalence % and 95% Confidence Intervals (Wilson Score)
# Input: Pathogen Name, Method Name, Count of Positives, Total Tested
# Output: A dataframe row formatted for Table 1
# ------------------------------------------------------------------------------
calc_prev_row <- function(name, method, positives, n_total) {
  # We use prop.test with correct=FALSE to get the Wilson Score Interval,
  # which is preferred for proportions (especially when prevalence is low/high).
  test <- suppressWarnings(prop.test(positives, n_total, conf.level = 0.95, correct = FALSE))
  
  data.frame(
    Parasites = name,
    Diagnostic_Method = method,
    Pos_Total = paste0(positives, " / ", n_total),
    Prevalence_Pct = round((positives/n_total)*100, 1),
    # Format CI as "Lower - Upper"
    CI_95 = paste0(round(test$conf.int[1]*100, 1), " - ", round(test$conf.int[2]*100, 1))
  )
}

# Generate rows for every parasite reported in Table 1
r1 <- calc_prev_row("Toxocara spp.", "Fecal Flotation", sum(df$rounds_pos, na.rm=T), n_float)
r2 <- calc_prev_row("Ancylostoma spp.", "Fecal Flotation", sum(df$hooks_pos, na.rm=T), n_float)
r3 <- calc_prev_row("Trichuris spp.", "Fecal Flotation", 0, n_float) # Explicitly 0 detected
r4 <- calc_prev_row("Coccidia", "Fecal Flotation", sum(df$coccid_pos, na.rm=T), n_float)
r5 <- calc_prev_row("Giardia spp.", "Coproantigen ELISA", sum(df$giardia_pos, na.rm=T), n_antigen)
r6 <- calc_prev_row("Cryptosporidium spp.", "Coproantigen ELISA", sum(df$crypto_pos, na.rm=T), n_antigen)
r7 <- calc_prev_row("Toxoplasma gondii", "Serology (MAT)", sum(df$toxo_Ab_pos, na.rm=T), n_sero)

# Combine rows into final Table 1
table1 <- rbind(r1, r2, r3, r4, r5, r6, r7)

# Print to console for checking
print("--- TABLE 1: PREVALENCE ---")
print(table1)

# Additional Coinfection Stats for the Results Text (Not in Table 1)
coinf <- calc_prev_row("Coinfection (2+)", "Flotation", sum(df$is_coinfected, na.rm=T), n_float)
spec_coinf <- calc_prev_row("Toxocara + Ancylostoma", "Flotation", sum(df$toxo_ancy_coinfection, na.rm=T), n_float)
print("--- COINFECTION TEXT STATS ---")
print(coinf)
print(spec_coinf)

# ==============================================================================
# 4. ANALYSIS PART 2: RISK FACTORS & INTENSITY (TABLE 2 GENERATION)
# ==============================================================================
# This section produces Table 2, which combines "Who has it?" (Prevalence/Risk)
# with "How bad is it?" (Intensity/EPG).

# ------------------------------------------------------------------------------
# FUNCTION: get_risk_stats
# Purpose: Calculates Prevalence %, Odds Ratio (OR), CI, and P-value.
# Method: Fisher's Exact Test (Robust for small sample sizes).
# ------------------------------------------------------------------------------
get_risk_stats <- function(data, group_col, outcome_col, target_group) {
  # 1. Subset data to rows where both Group and Outcome are known (no NAs)
  sub <- data[!is.na(data[[group_col]]) & !is.na(data[[outcome_col]]), ]
  
  # 2. Create a 2x2 contingency table (e.g., Young/Adult vs Pos/Neg)
  tbl <- table(sub[[group_col]], sub[[outcome_col]])
  
  # 3. Calculate raw prevalence % for the specific target group (e.g., % of Young cats positive)
  n_target <- sum(sub[[group_col]] == target_group)
  n_pos    <- sum(sub[[group_col]] == target_group & sub[[outcome_col]] == 1)
  prev_pct <- ifelse(n_target > 0, round((n_pos / n_target) * 100, 1), 0)
  
  # 4. Run Fisher's Exact Test
  # We check if the table is valid (at least 2x2) to avoid errors.
  if(nrow(tbl) < 2 || ncol(tbl) < 2) return(list(prev=prev_pct, or=NA, ci="NC", p=1))
  
  ft <- fisher.test(tbl)
  
  # Return list of results
  list(
    prev = prev_pct,
    or   = round(ft$estimate, 1), # Odds Ratio
    ci   = paste0(round(ft$conf.int[1], 1), "-", round(ft$conf.int[2], 1)), # 95% CI
    p    = ft$p.value # Raw P-value
  )
}

# ------------------------------------------------------------------------------
# FUNCTION: get_int_stats
# Purpose: Calculates Median EPG and statistical difference in intensity.
# Methods: 
#   - Negative Binomial GLM ("nb"): Used for Toxocara (highly overdispersed).
#   - Mann-Whitney U ("mw"): Used for Ancylostoma (standard non-parametric test).
# ------------------------------------------------------------------------------
get_int_stats <- function(data, group_col, epg_col, target_group, test_method="mw") {
  # 1. Filter: We only analyze intensity among POSITIVE animals (EPG > 0)
  sub_pos <- data[!is.na(data[[group_col]]) & !is.na(data[[epg_col]]) & data[[epg_col]] > 0, ]
  
  # 2. Calculate Median EPG for the target group
  vals <- sub_pos[[epg_col]][sub_pos[[group_col]] == target_group]
  median_val <- ifelse(length(vals) > 0, round(median(vals), 1), 0)
  
  # 3. Statistical Test (Compares the groups)
  p_val <- NA
  
  if (test_method == "nb") {
    # --- NEGATIVE BINOMIAL GLM ---
    # Appropriate for Toxocara because "Variance >> Mean" (Super-shedders)
    f <- as.formula(paste(epg_col, "~", group_col))
    try({
      model <- glm.nb(f, data = data) # Note: Runs on full data to model counts correctly
      p_val <- summary(model)$coefficients[2, 4] # Extract P-value for the group coefficient
    }, silent=TRUE)
    
  } else {
    # --- MANN-WHITNEY U TEST ---
    # Non-parametric rank sum test, standard for skewed data
    f <- as.formula(paste(epg_col, "~", group_col))
    try({
      test <- wilcox.test(f, data = sub_pos)
      p_val <- test$p.value
    }, silent=TRUE)
  }
  
  list(median=median_val, p=p_val)
}

# ------------------------------------------------------------------------------
# RUN CALCULATIONS FOR TABLE 2
# ------------------------------------------------------------------------------

# --- 1. TOXOCARA (Roundworm) ---
# Risk (Prevalence)
tox_age_risk <- get_risk_stats(df, "Age_Group", "rounds_pos", "Young")
tox_sex_risk <- get_risk_stats(df, "Sex_Clean", "rounds_pos", "Male")
# Reference group prevalence (Adult / Female)
tox_age_ref  <- get_risk_stats(df, "Age_Group", "rounds_pos", "Adult")
tox_sex_ref  <- get_risk_stats(df, "Sex_Clean", "rounds_pos", "Female")

# Intensity (Using GLM "nb" because Toxocara is overdispersed)
tox_age_int <- get_int_stats(df, "Age_Group", "roundsepg", "Young", "nb")
tox_sex_int <- get_int_stats(df, "Sex_Clean", "roundsepg", "Male", "nb")
# Reference group medians
tox_age_int_ref <- get_int_stats(df, "Age_Group", "roundsepg", "Adult", "nb")
tox_sex_int_ref <- get_int_stats(df, "Sex_Clean", "roundsepg", "Female", "nb")


# --- 2. ANCYLOSTOMA (Hookworm) ---
# Risk (Prevalence)
anc_age_risk <- get_risk_stats(df, "Age_Group", "hooks_pos", "Young")
anc_sex_risk <- get_risk_stats(df, "Sex_Clean", "hooks_pos", "Male")
# Reference group prevalence
anc_age_ref  <- get_risk_stats(df, "Age_Group", "hooks_pos", "Adult")
anc_sex_ref  <- get_risk_stats(df, "Sex_Clean", "hooks_pos", "Female")

# Intensity (Using Mann-Whitney "mw" - standard test)
anc_age_int <- get_int_stats(df, "Age_Group", "hooksepg", "Young", "mw")
anc_sex_int <- get_int_stats(df, "Sex_Clean", "hooksepg", "Male", "mw")
# Reference group medians
anc_age_int_ref <- get_int_stats(df, "Age_Group", "hooksepg", "Adult", "mw")
anc_sex_int_ref <- get_int_stats(df, "Sex_Clean", "hooksepg", "Female", "mw")

# ------------------------------------------------------------------------------
# CONSTRUCT FINAL TABLE 2 DATAFRAME
# ------------------------------------------------------------------------------
# Formatting helper for P-values
fmt_p <- function(p) {
  if (is.na(p)) return("-")
  if (p < 0.001) return("<0.001")
  return(format(round(p, 3), nsmall=3))
}

# Assemble the table column by column
table2 <- data.frame(
  Demographic_Group = c("Age", "Young (<1 yr)", "Adult (>=1 yr)", "Sex", "Male", "Female"),
  
  # Sample Sizes (N)
  N = c("", sum(df$Age_Group=="Young", na.rm=T), sum(df$Age_Group=="Adult", na.rm=T),
        "", sum(df$Sex_Clean=="Male", na.rm=T), sum(df$Sex_Clean=="Female", na.rm=T)),
  
  # --- Toxocara Columns ---
  # Column 1: Prevalence % (Odds Ratio; 95% CI)
  Tox_Prev_OR_CI = c(
    "",
    paste0(tox_age_risk$prev, "%\n(", tox_age_risk$or, "; ", tox_age_risk$ci, ")"),
    paste0(tox_age_ref$prev, "%\n(Reference)"),
    "",
    paste0(tox_sex_risk$prev, "%\n(", tox_sex_risk$or, "; ", tox_sex_risk$ci, ")"),
    paste0(tox_sex_ref$prev, "%\n(Reference)")
  ),
  # Column 2: Median EPG (P-value)
  Tox_Med_EPG_P = c(
    "",
    paste0(tox_age_int$median, "\n(", fmt_p(tox_age_int$p), ")"),
    paste0(tox_age_int_ref$median),
    "",
    paste0(tox_sex_int$median, "\n(", fmt_p(tox_sex_int$p), ")"),
    paste0(tox_sex_int_ref$median)
  ),
  
  # --- Ancylostoma Columns ---
  # Column 3: Prevalence % (Odds Ratio; 95% CI)
  Anc_Prev_OR_CI = c(
    "",
    paste0(anc_age_risk$prev, "%\n(", anc_age_risk$or, "; ", anc_age_risk$ci, ")"),
    paste0(anc_age_ref$prev, "%\n(Reference)"),
    "",
    paste0(anc_sex_risk$prev, "%\n(", anc_sex_risk$or, "; ", anc_sex_risk$ci, ")"),
    paste0(anc_sex_ref$prev, "%\n(Reference)")
  ),
  # Column 4: Median EPG (P-value)
  Anc_Med_EPG_P = c(
    "",
    paste0(anc_age_int$median, "\n(", fmt_p(anc_age_int$p), ")"),
    paste0(anc_age_int_ref$median),
    "",
    paste0(anc_sex_int$median, "\n(", fmt_p(anc_sex_int$p), ")"),
    paste0(anc_sex_int_ref$median)
  )
)

print("--- TABLE 2: DEMOGRAPHICS AND BURDEN ---")
print(table2)

# ==============================================================================
# 5. EXPORT TO CSV
# ==============================================================================
write.csv(table1,
          file = file.path(output_tables, "Table1.csv"),
          row.names = FALSE)

write.csv(table2,
          file = file.path(output_tables, "Table2.csv"),
          row.names = FALSE)

