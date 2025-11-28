# -------------------------------------------------------------------------------
# Description:
#   Build systemic-therapy lines from MASTER clinical data and compute
#   Time-to-Progression (TTP) and TTPr (line k+1 / line k) for LMS patients.
#   Progression is assigned within each line’s window; TTPr > 1.3 → Responder (R).
#
# Required Inputs:
#   - MASTER clinical exports (XLSX/CSV) under data/rawdata/validation/MASTER
#   - OncoTree mapping: sts_bone_oncotree_v2.csv
#
# User-defined parameters:
#   - dir_in, dir_out       : input/output roots
#   - cancer subtype filter : LMS (default; can change to other STS)
#
# Analyses performed:
#   - Parse rows → therapy/progression; build per-line chronology
#   - Tag systemic lines by prefix (CHT/CTH/ST/HCHT/EIA/DST) and drug text
#   - Assign earliest progression in [start, next_start); compute TTP
#   - Clean TTP (drop missing/negative); compute TTPr and R/NR flags
#   - Produce strict flags where both lines use confirmed progression
#
# Notes:
#   - Update non-systemic prefix blacklist if new local codes appear.
#   - Keep intermediates for clinical audit/review.
# -------------------------------------------------------------------------------
######################################################
## load libraries
######################################################
library(readxl)
library(dplyr)
library(tidyr)
library(readr)
library(tidyverse)
library(lubridate)
library(stringr)

######################################################
## set up directory
######################################################

dir_in <- 'data/rawdata'
dir_out <- 'data/procdata'

######################################################
## functions to extract drug names 
######################################################

extract_prefix <- function(txt){
  txt <- str_to_upper(str_squish(coalesce(txt,"")))
  m <- str_match(txt, "^([A-ZÄÖÜ]+)\\s*:")
  m[,2]
}

# Keep drug names exactly as they appear (no German→English, no title-casing)
extract_drugs <- function(x, treat_CT_as_imaging = TRUE){
  if (is.na(x) || !nzchar(x)) return(character(0))

  upper <- str_to_upper(str_squish(x))

  # Exclude clearly non-systemic rows by prefix
  if (str_detect(upper, "^(OP|RT|HSZT)\\s*:")) return(character(0))

  # Choose which prefixes to strip off the front of the string
  prefix_pattern <- if (treat_CT_as_imaging){
    "^(ST|CHT|CTH|HCHT|DST)\\s*:\\s*"
  } else {
    "^(ST|CHT|CTH|HCHT|DST|CT)\\s*:\\s*"
  }

  # Remove prefix and obvious parentheses/protocol notes
  txt <- x %>%
    str_replace(regex(prefix_pattern, ignore_case=TRUE), "") %>%
    str_replace_all("\\([^)]*\\)", " ") %>%
    str_replace_all("(?i)protokoll|protocol|schema|regime|im\\s+rahmen.*", " ")

  # Normalize delimiters to "|"
  txt <- txt %>%
    str_replace_all("\\s*[+;/]\\s*", "|") %>%
    str_replace_all("\\s*,\\s*", "|") %>%
    str_replace_all("\\|{2,}", "|") %>%
    str_squish()

  parts <- str_split(txt, "\\|") %>% unlist() %>% str_squish()
  parts <- parts[nzchar(parts)]
  if (!length(parts)) return(character(0))

  # Keep as-is (no normalization), but drop blatantly non-drug tokens if any
  parts <- parts[ str_detect(parts, "[A-Za-z]") ]

  unique(parts)
}

# extract RT dose
extract_rt_dose <- function(text) {
  # Return numeric dose from "GD XX Gy", e.g. 59.4 or 60
  str_extract(text, "(?i)GD\\s*[:\\s]?\\d+[.,]?\\d*") %>%  # match 'GD 59.4' or 'GD:60'
    str_replace_all(",", ".") %>%
    str_extract("\\d+\\.?\\d*") %>%
    as.numeric()
}

######################################################
## load NCT-MASTER data
######################################################
file_path <- file.path(dir_in, 'validation/MASTER', 'DKTK-P141_Clinical_Info.xlsx')
sheet_names <- excel_sheets(file_path)
sheet_names <- sheet_names[c(1,7,6)]

data_list <- lapply(1:length(sheet_names), function(k){

  df <- read_excel(file_path, sheet = sheet_names[k])
  colnames(df)[1] <- 'patientid'
  if(sheet_names[k] == 'Clinical Information'){
   
   df <- df[!duplicated(df$patientid), ]

  }

 df

})

# merge all sheets using Reduce and merge --> sheets 1 and 7
merged_data <- Reduce(function(x, y) full_join(x, y, by = "patientid"), data_list[1:2]) # 2020 patients

######################################################
## Add overall survival outcome
#####################################################

merged_data <- merged_data %>%
  mutate(
    date_diagnosis = as.Date(`Date of Diagnosis`, format = "%d.%m.%Y"),
    date_death     = as.Date(Deathdate, format = "%d.%m.%Y"),

    os_status = if_else(!is.na(date_death), 1, 0),
    os_event_date = date_death,

    os_days = as.numeric(os_event_date - date_diagnosis),
    os_months = os_days / 30.44
  ) %>%
  filter(!is.na(date_diagnosis))

merged_data$os_months <- round(merged_data$os_months) 
write.csv(merged_data, file = file.path(dir_out, 'validation', 'master_cohort_part1.csv'), row.names =FALSE)

######################################################
## Keep STS patients
######################################################
# load OncoTree data (Level 2)
annot <- read.csv(file.path(dir_in, 'validation/MASTER', 'sts_bone_oncotree_v2.csv'))
annot <- annot[annot$tissue == 'Soft Tissue', ] # 38 syubtypes

# Keep main subtypes 
annot_filtered <- annot %>%
  filter(oncotree_name %in% c(
    "Leiomyosarcoma",
    "Liposarcoma",
    "Myxofibrosarcoma",
    "Synovial Sarcoma",
    "Malignant Peripheral Nerve Sheath Tumor",
    "Undifferentiated Pleomorphic Sarcoma/Malignant Fibrous Histiocytoma/High-Grade Spindle Cell Sarcoma"
  ))

######################################################
## Let's focus on main STS subtypes
######################################################
merged_data <- read.csv(file.path(dir_out, 'validation', 'master_cohort_part1.csv'))

# keep sts patients
dat_sts <- merged_data[merged_data$'Oncotree.Code' %in% annot$oncotree_code, ] 
length(unique(dat_sts$patientid)) # 824 sts patients 
write.csv(dat_sts, file = file.path(dir_out, 'validation', 'master_cohort_sts.csv'), row.names = FALSE)

# number of patients per subtypes
group <- unique(dat_sts$Oncotree.Code) # 29 sts subtypes
subtype_dist <- lapply(1:length(group), function(k){
 
 df <- dat_sts[dat_sts$Oncotree.Code == group[k], ]
 data.frame(Oncotree_Code = group[k],
            oncotree_name = annot[annot$oncotree_code == group[k], 'oncotree_name'],
            n_patient = length(unique(df$patientid)))

})

subtype_dist <- do.call(rbind, subtype_dist)
subtype_dist <- subtype_dist[order(subtype_dist$n, decreasing=TRUE), ]
write.csv(subtype_dist, file = file.path(dir_out, 'validation', 'distribution_sts_subtypes.csv'), row.names = FALSE)

# keep main sts patients
dat_sts <- merged_data[merged_data$'Oncotree.Code' %in% annot_filtered$oncotree_code, ] 
length(unique(dat_sts$patientid)) # 390 main sts patients 
write.csv(dat_sts, file = file.path(dir_out, 'validation', 'master_cohort_main_sts_subtypes.csv'), row.names = FALSE)

######################################################
## Let's focus on LMS subtypes
######################################################
subtype <- 'LIPO'
dat_sts <- read.csv(file.path(dir_out, 'validation', 'master_cohort_main_sts_subtypes.csv'))
#lms_dat <- dat_sts[dat_sts$Oncotree.Code == 'LMS', ] # 171 patients
#ss_dat <- dat_sts[dat_sts$Oncotree.Code == 'SYNS', ] # 102 patients
#mfh_dat <- dat_sts[dat_sts$Oncotree.Code == 'MFH', ] # 62 patients
#mfs_dat <- dat_sts[dat_sts$Oncotree.Code == 'MFS', ] # 36 patients
lipo_dat <- dat_sts[dat_sts$Oncotree.Code == 'LIPO', ] # 19 patients
write.csv(lipo_dat, file = file.path(dir_out, 'validation', subtype, 'master_cohort.csv'), row.names = FALSE)

raw <- read_csv(
  file = file.path(dir_out, 'validation', subtype, "master_cohort.csv"),
  na = c("", "NA", "Na", "na", "N/A", "n/a", "NULL", "null", "#N/A"),
  guess_max = 100000
)

cn <- names(raw)

# Robust column picking
pick_col <- function(cn, patterns) {
  for (p in patterns) {
    hit <- grep(p, cn, ignore.case = TRUE, value = TRUE)
    if (length(hit) > 0) return(hit[1])
  }
  return(NA_character_)
}

pid_col   <- pick_col(cn, c("patientid", "^MASTER\\s*pid$", "^pid$"))
sex_col   <- pick_col(cn, c("Sex", "geschlecht"))
diagnosis_col <- pick_col(cn, c("Date.of.Diagnosis", "diagnosis"))
death_col     <- pick_col(cn, c("Deathdate", "verstorben", "tod"))
os_status <- pick_col(cn, c("os_status"))
os_days <- pick_col(cn, c("os_days"))
os_months <- pick_col(cn, c("os_months"))
start_col <- pick_col(cn, c("Therapy.start.date", "start", "beginn", "therapy_start"))
end_col   <- pick_col(cn, c("Therapy.end.date", "end", "therapy_end", "ende"))
stage_col <- pick_col(cn, c("Staging.method", "therapy.type", "treatment", "stage"))
remis_col <- pick_col(cn, c("Grouped.Remissionstatus", "remission", "status"))

picked <- c(pid_col, sex_col, diagnosis_col, death_col, os_status, os_days, os_months, start_col, end_col, stage_col, remis_col)
names(picked) <- c("pid", "sex", "diagnosis","death",  "os_status", "os_days", "os_months", "start", "end", "stage_type", "remiss")
print(picked)

if (any(is.na(picked))) {
  stop("Column detection failed. Missing: ", paste(names(picked)[is.na(picked)], collapse = ", "))
}

# ---- 3) parse dates & basic cleanup ----

parse_dt <- function(x) parse_date_time(
  x,
  orders = c("ymd","Y-m-d","d.m.Y","d.m.y","m/d/Y","m/d/y")
)

dat <- raw %>%
  transmute(
    pid            = .data[[pid_col]],
    sex            = .data[[picked["sex"]]],
    date_diagnosis = parse_dt(.data[[picked["diagnosis"]]]),
    date_death     = parse_dt(.data[[picked["death"]]]),
    
    # Keep these columns as-is (already numeric/logical, no need to parse dates)
    os_status      = .data[[picked["os_status"]]],
    os_days        = .data[[picked["os_days"]]],
    os_months      = .data[[picked["os_months"]]],
    
    start_raw      = .data[[start_col]],
    end_raw        = .data[[end_col]],
    stage_type     = .data[[stage_col]],   # therapy description / staging info
    remiss         = .data[[remis_col]]    # remission / progression status
  ) %>%
  mutate(
    start      = parse_dt(start_raw),
    end        = parse_dt(end_raw),
    stage_type = str_squish(as.character(stage_type)),
    remiss     = str_squish(as.character(remiss))
  ) %>%
  filter(!is.na(start)) %>%           # cannot use rows without start date
  arrange(pid, start) %>%
  mutate(
    start = as_date(start),
    end   = as_date(end)
  ) %>%
  # classify events
  mutate(
    is_progress = str_detect(tolower(remiss), "progress"),
    is_deceased = str_detect(tolower(remiss), "deceased"),
    is_therapy  = !is_progress & !is_deceased
  ) %>%
  # drop rows where we have no idea what they are
  filter(!(is.na(is_progress) & is.na(is_deceased) & is.na(is_therapy)))

dat %>% dplyr::count(is_progress, is_deceased, is_therapy)

# ---- 5) build therapy lines (one row per therapy line) ----
therapy <- dat %>%
  filter(is_therapy) %>%
  group_by(pid) %>%
  arrange(start, .by_group = TRUE) %>%
  mutate(
    line_no    = row_number(),
    next_start = lead(start),
    therapy_text = remiss  # optional renaming
  ) %>%
  ungroup()

# (Optional) If you want only systemic therapies for TTPr:
# therapy <- therapy %>%
#  filter(!str_detect(tolower(stage_type), "op|surgery|rt|radiation"))

# ---- 6) Prefix lists ----
systemic_prefix_keep    <- c("CHT","ST","CTH","HCHT","EIA","DST")
nonsystemic_prefix_drop <- c(
  "OP","OPR","S","ALG","KT",
  "RT","IMRT","IORT","RCT","RCTH","RCHT","RTCT",
  "PRRT","SIRT","TACE","HIPEC","ILP",
  "RHT","LITT","MWA","TECC"
)

therapy_sys <- therapy %>%
  mutate(
    therapy_text_clean = stringr::str_squish(coalesce(therapy_text, "")),
  # stage_type_clean   = stringr::str_squish(coalesce(stage_type, "")),

    prefix = extract_prefix(therapy_text_clean),

    # Categorize by prefix
    modality_by_prefix = dplyr::case_when(
      prefix %in% systemic_prefix_keep     ~ "systemic",
      prefix %in% nonsystemic_prefix_drop  ~ "nonsystemic",
      TRUE                                 ~ NA_character_
    ),

    # Parse drugs (optional: treat "CT" as imaging, not systemic)
    drugs_vec = purrr::map(therapy_text_clean, ~extract_drugs(.x, treat_CT_as_imaging = TRUE)),
    drug_list = purrr::map_chr(drugs_vec, ~ if (length(.x)) paste(.x, collapse = "; ") else NA_character_),
    has_drug_term = !is.na(drug_list),

    # Final systemic flag logic
    is_systemic = dplyr::case_when(
      modality_by_prefix == "systemic"          ~ TRUE,
      modality_by_prefix == "nonsystemic"       ~ FALSE,
      is.na(modality_by_prefix) & has_drug_term ~ TRUE,
      TRUE                                      ~ FALSE
    ),

    # Combine text for OP/RT flag detection
    full_text = tolower(therapy_text_clean),

    # Use word boundaries to avoid false positives
    op_flag = str_detect(full_text, "\\bop\\b|resect|surgery|operation|chirurg"),
    rt_flag = str_detect(full_text, "\\brt\\b|radiation|radiotherapie|bestrahlung"),

    # Label therapy phase
    therapy_phase = dplyr::case_when(
      op_flag & rt_flag ~ "Post-OP-RT",
      op_flag           ~ "Post-OP",
      rt_flag           ~ "Post-RT",
      TRUE              ~ "No-OP-RT"
    )
  ) 

write.csv(therapy_sys, file = file.path(dir_out, 'validation', subtype, 'step1.csv'), row.names = FALSE)

### 4) Extract RT dose

therapy_sys <- therapy_sys %>%
  mutate(
    rt_dose_gy = if_else(
      rt_flag & !is.na(therapy_text),
      extract_rt_dose(therapy_text),
      NA_real_
    )
  )

write.csv(therapy_sys,
          file.path(dir_out, 'validation', subtype, "step1_therapy_sys.csv"),
          row.names = FALSE)

### 5) Patient-level exclusion: ANY RT < 25 Gy

patients_lowdose_rt <- therapy_sys %>%
  filter(rt_flag, !is.na(rt_dose_gy), rt_dose_gy < 25) %>%
  distinct(pid)

write.csv(patients_lowdose_rt,
          file.path(dir_out, 'validation', subtype, "patients_excluded_rt_lt25.csv"),
          row.names = FALSE)

eligible_pids <- setdiff(unique(therapy_sys$pid), patients_lowdose_rt$pid)

therapy_sys_elig <- therapy_sys %>% filter(pid %in% eligible_pids)
dat_elig <- dat %>% filter(pid %in% eligible_pids)

# ---- 7) progress table ----
prog_tbl <- dat_elig %>%
  filter(is_progress) %>%
  transmute(pid, prog_date = start)  # progress date is stored in the "start" column for progress rows

therapy_with_prog <- therapy_sys %>%
  rowwise() %>%
  mutate(
    progress_date = {
      cand <- prog_tbl$prog_date[
        prog_tbl$pid == pid &
        prog_tbl$prog_date >= start &
        (is.na(next_start) | prog_tbl$prog_date < next_start)
      ]
      if (length(cand) == 0) NA_Date_ else min(cand, na.rm = TRUE)
    }
  ) %>%
  ungroup()

# ---- 8) compute TTP fields ----
therapy_ttp <- therapy_with_prog %>%
  mutate(
    # anchor: prefer confirmed progress; else use end as proxy; else missing
    ttp_anchor = coalesce(progress_date, end),
    ttp_source = case_when(
      !is.na(progress_date)               ~ "confirmed_progress",
      is.na(progress_date) & !is.na(end)  ~ "proxy_end",
      TRUE                                ~ "missing"
    ),
    # days from start to anchor
    ttp_days = as.numeric(ttp_anchor - start)
  ) 

# unlist drugs_vec column
therapy_ttp <- therapy_ttp %>%
  mutate(
    drugs_vec = sapply(drugs_vec, function(x) {
      if (length(x) == 0 || all(is.na(x))) {
        NA_character_
      } else {
        paste(x, collapse = "; ")
      }
    })
  )

write.csv(therapy_ttp, file = file.path(dir_out, 'validation', subtype, 'step2.csv'), row.names = FALSE)

therapy_ttp_all <- therapy_ttp %>%
  dplyr::mutate(
    exclusion_reason = dplyr::case_when(
      is.na(ttp_days) ~ "no_anchor_date",
      ttp_days < 0    ~ "negative_ttp",
      TRUE            ~ NA_character_
    ),
    is_valid_ttp = is.na(exclusion_reason)
  )

write.csv(therapy_ttp_all, file = file.path(dir_out, 'validation', subtype, 'step3.csv'), row.names = FALSE)

# Optional: Save filtered clean version
# therapy_ttp_clean <- therapy_ttp_all %>% filter(is_valid_ttp)
# write.csv(therapy_ttp_clean, file = file.path(dir_out, 'validation/LMS', 'lms_step3_clean.csv'), row.names = FALSE)

# Optional: Save just the excluded rows (for QA)
# ttp_exclusions <- therapy_ttp_all %>% filter(!is_valid_ttp)
# write.csv(ttp_exclusions, file = file.path(dir_out, "validation/LMS", "ttp_exclusions_log.csv"), row.names = FALSE)

therapy_ttpr <- therapy_ttp_all %>%
  dplyr::group_by(pid) %>%
  dplyr::arrange(line_no, .by_group = TRUE) %>%
  dplyr::mutate(
    prev_ttp_days       = dplyr::lag(ttp_days),
    prev_ttp_source     = dplyr::lag(ttp_source),

    # Compute ratio if both TTPs exist
    ttpr_vs_prev = ifelse(!is.na(ttp_days) & !is.na(prev_ttp_days),
                          ttp_days / prev_ttp_days, NA_real_),

    # Define valid pair only when both TTPs are known
    pair_ok_for_ttpr = !is.na(prev_ttp_days) & !is.na(ttp_days),

    # Both lines confirmed by progression (for strict evaluation)
    pair_both_confirmed =
      dplyr::lag(ttp_source == "confirmed_progress") &
      (ttp_source == "confirmed_progress"),

    # Define TTPr flags (standard)
    ttpr_flag = dplyr::case_when(
      pair_ok_for_ttpr & ttpr_vs_prev > 1.3 ~ "R",     # responder
      pair_ok_for_ttpr                      ~ "NR",    # non‑responder
      TRUE                                  ~ NA_character_  # not evaluable
    ),

    # Define TTPr flags (strict version — both confirmed progressions)
    ttpr_flag_strict = dplyr::case_when(
      pair_ok_for_ttpr & pair_both_confirmed & ttpr_vs_prev > 1.3 ~ "R",
      pair_ok_for_ttpr & pair_both_confirmed                      ~ "NR",
      TRUE                                                        ~ NA_character_
    )
  ) %>%
  dplyr::ungroup()

write.csv(therapy_ttpr, file.path(dir_out, 'validation', subtype, "step4.csv"), row.names=FALSE)

patient_summary_flags <- therapy_ttpr %>%
  group_by(pid) %>%
  summarise(
    n_total_lines     = n(),
    n_valid_ttps      = sum(!is.na(ttp_days)),
    n_ttpr_computed   = sum(!is.na(ttpr_flag)),
    had_responder     = any(ttpr_flag == "R", na.rm = TRUE),
    had_op            = any(op_flag, na.rm = TRUE),
    had_rt            = any(rt_flag, na.rm = TRUE),
    n_systemic_lines  = sum(is_systemic, na.rm = TRUE)
  ) %>%
  mutate(
    ttpr_eligible = n_ttpr_computed >= 1  # Can compute at least one TTPr
  )

patient_summary_flags <- patient_summary_flags %>%
  mutate(
    no_op_rt = !had_op & !had_rt  # TRUE if neither surgery nor RT was given
  )
  
write.csv(patient_summary_flags,
          file = file.path(dir_out, 'validation', subtype, "patient_ttpr_summary.csv"),
          row.names = FALSE)

#therapy_ttpr <- therapy_ttpr %>%
#  mutate(
#    # Optional: flag if systemic therapy only
#    keep_for_ttpr = is_systemic & is_valid_ttp
#  )

#therapy_ttpr %>% filter(keep_for_ttpr, !op_flag, !rt_flag)

## Include flags
therapy_ttpr <- therapy_ttpr %>%
  arrange(pid, line_no) %>%
  group_by(pid) %>%
  mutate(
    prior_op      = lag(op_flag),
    prior_rt      = lag(rt_flag),
    prior_rt_dose = lag(rt_dose_gy),

    ever_had_op = any(op_flag, na.rm = TRUE),
    ever_had_rt = any(rt_flag, na.rm = TRUE),

    include_strict = is_systemic == TRUE & (
      prior_op == TRUE |
      (prior_rt == TRUE & !is.na(prior_rt_dose) & prior_rt_dose >= 25) |
      (!ever_had_op & !ever_had_rt)
    ),

    include_relaxed = is_systemic == TRUE & (
      prior_op == TRUE |
      prior_rt == TRUE |
      (!ever_had_op & !ever_had_rt)
    )
  ) %>%
  ungroup()

write.csv(therapy_ttpr, file.path(dir_out, 'validation', subtype, "step5.csv"), row.names=FALSE)

######################################################################################################
############################################ Filtration Steps ########################################
######################################################################################################

non_drug_keywords <- c(
  "disease control","therapiepause","keine therapie",
  "resektion","operation","op","chirurgie","rfa","sirt","tae",
  "hyperthermie","mw[aä]","iort","hipec","ire","mistel","rt",
  "gd","bestrahlung","photonen","schilddrüsenmet","ed",
  "radiotherapie","neoadjuvant","unclear","unklar",
  "immunbiologische"
)

# Step 1: Filter valid systemic lines and exclude non-drug entries
strict_lines <- therapy_ttpr %>%
  filter(
    include_strict == TRUE,
    is_valid_ttp  == TRUE,
    is_systemic   == TRUE,
    !is.na(drug_list),
    !str_detect(tolower(drug_list), paste(non_drug_keywords, collapse="|"))
  ) %>%
  arrange(pid, line_no)

strict_ttpr_patients <- strict_lines %>%
  group_by(pid) %>%
  filter(!is.na(ttpr_vs_prev)) %>%
  pull(pid) %>% unique()

strict_final <- strict_lines %>%
  filter(pid %in% strict_ttpr_patients)

write.csv(strict_final, file.path(dir_out, 'validation', subtype, "strict_ttpr.csv"), row.names=FALSE)

# Filtered dataset: relaxed inclusion

relaxed_lines <- therapy_ttpr %>%
  filter(
    include_relaxed == TRUE,
    is_valid_ttp    == TRUE,
    is_systemic     == TRUE,
    !is.na(drug_list),
    !str_detect(tolower(drug_list), paste(non_drug_keywords, collapse="|"))
  ) %>%
  arrange(pid, line_no)

relaxed_ttpr_patients <- relaxed_lines %>%
  group_by(pid) %>%
  filter(!is.na(ttpr_vs_prev)) %>%
  pull(pid) %>% unique()

relaxed_final <- relaxed_lines %>%
  filter(pid %in% relaxed_ttpr_patients)

write.csv(relaxed_final, file.path(dir_out, 'validation', subtype, "relaxed_ttpr.csv"), row.names=FALSE)
