# =============================================================================
# Heatwave mortality among older adults in Region Stockholm, 2011-2024
# Primary analysis: heatwave (HW) and optimal temperature (OT) days defined
# from the cumulative relative risk (cumRR) of a distributed lag non-linear model
#
# The register data are not publicly available. Paths are relative to the
# project folder:
#   work/data_250307/  register extracts (SAS files)
#   documents/         lookup tables (ICD-10 blocks, ATC groups)
#   data/              DLNM output, intermediate files and fitted models
#   download/          result tables
#
# ssb0 = community dwellers, ssb1 = nursing home residents (SSB)
# Register variable names are kept as in the source data (lopnr = person ID,
# KOEN = sex, FMANAD = birth year-month, mosaicgrupp_sthlm = Mosaic group)
# =============================================================================

library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(lubridate)
library(haven)
library(readxl)
library(writexl)
library(caret)
library(gbm)
library(pROC)
library(ggplot2)

# 1. Daily cumulative relative risk (cumRR) over lags 0-7 from the stratum-specific DLNM
{
# Nursing home residents
load(file = "data/timeserie_ssb1_max.RData")
load(file = "data/pred_ssb1_sum.RData")

timeserie_ssb1_max <- timeserie_ssb1_max %>%
  mutate(year = substr(date, 1,4))

summer_data_ssb1 <- subset(
  timeserie_ssb1_max,
  format(date, "%m-%d") >= "03-25" &
    format(date, "%m-%d") <= "09-30"
)

summer_data_ssb1 <- summer_data_ssb1 %>%
  arrange(date) %>%
  mutate(
    temp_lag0 = temp,
    temp_lag1 = lag(temp,1),
    temp_lag2 = lag(temp,2),
    temp_lag3 = lag(temp,3),
    temp_lag4 = lag(temp,4),
    temp_lag5 = lag(temp,5),
    temp_lag6 = lag(temp,6),
    temp_lag7 = lag(temp,7)
  )

# RR for a given lag at the nearest temperature in the DLNM prediction grid
getRR <- function(temp, lagcol){
  if(is.na(temp)) return(NA_real_)

  temps_pred <- as.numeric(rownames(pred_ssb1_sum$matRRfit))

  idx <- which.min(abs(temps_pred - temp))

  pred_ssb1_sum$matRRfit[idx, lagcol]
}

summer_data_ssb1 <- summer_data_ssb1 %>%
  rowwise() %>%
  mutate(
    RR0 = getRR(temp_lag0,"lag0"),
    RR1 = getRR(temp_lag1,"lag1"),
    RR2 = getRR(temp_lag2,"lag2"),
    RR3 = getRR(temp_lag3,"lag3"),
    RR4 = getRR(temp_lag4,"lag4"),
    RR5 = getRR(temp_lag5,"lag5"),
    RR6 = getRR(temp_lag6,"lag6"),
    RR7 = getRR(temp_lag7,"lag7")
  ) %>%
  ungroup()

summer_data_ssb1 <- summer_data_ssb1 %>%
  mutate(
    cumRR_ssb1 = RR0 * RR1 * RR2 * RR3 *
      RR4 * RR5 * RR6 * RR7
  )

save(summer_data_ssb1, file = "data/summer_data_ssb1.RData")

rm(list = ls())
gc()

# Community dwellers
load(file = "data/timeserie_ssb0_max.RData")
load(file = "data/pred_ssb0_sum.RData")

timeserie_ssb0_max <- timeserie_ssb0_max %>%
  mutate(year = substr(date, 1,4))

summer_data_ssb0 <- subset(
  timeserie_ssb0_max,
  format(date, "%m-%d") >= "03-25" &
    format(date, "%m-%d") <= "09-30"
)

summer_data_ssb0 <- summer_data_ssb0 %>%
  arrange(date) %>%
  mutate(
    temp_lag0 = temp,
    temp_lag1 = lag(temp,1),
    temp_lag2 = lag(temp,2),
    temp_lag3 = lag(temp,3),
    temp_lag4 = lag(temp,4),
    temp_lag5 = lag(temp,5),
    temp_lag6 = lag(temp,6),
    temp_lag7 = lag(temp,7)
  )

getRR <- function(temp, lagcol){
  if(is.na(temp)) return(NA_real_)

  temps_pred <- as.numeric(rownames(pred_ssb0_sum$matRRfit))

  idx <- which.min(abs(temps_pred - temp))

  pred_ssb0_sum$matRRfit[idx, lagcol]
}

summer_data_ssb0 <- summer_data_ssb0 %>%
  rowwise() %>%
  mutate(
    RR0 = getRR(temp_lag0,"lag0"),
    RR1 = getRR(temp_lag1,"lag1"),
    RR2 = getRR(temp_lag2,"lag2"),
    RR3 = getRR(temp_lag3,"lag3"),
    RR4 = getRR(temp_lag4,"lag4"),
    RR5 = getRR(temp_lag5,"lag5"),
    RR6 = getRR(temp_lag6,"lag6"),
    RR7 = getRR(temp_lag7,"lag7")
  ) %>%
  ungroup()

summer_data_ssb0 <- summer_data_ssb0 %>%
  mutate(
    cumRR_ssb0 = RR0 * RR1 * RR2 * RR3 *
      RR4 * RR5 * RR6 * RR7
  )

save(summer_data_ssb0, file = "data/summer_data_ssb0.RData")

rm(list = ls())
gc()

# Merge both strata
load(file = "data/summer_data_ssb0.RData")
load(file = "data/summer_data_ssb1.RData")

cumRR_data <- inner_join(
  summer_data_ssb0 %>% select(date, temp, cumRR_ssb0),
  summer_data_ssb1 %>% select(date, cumRR_ssb1),
  by = "date"
) %>%
  filter(!is.na(cumRR_ssb1), !is.na(cumRR_ssb0))

save(cumRR_data, file = "data/cumRR_data.RData")

rm(list = ls())
gc()
}

# 2. HW and OT days per stratum (1 = HW, 0 = OT, NA = neither)
{
load(file = "data/cumRR_data.RData")

# Nursing home residents: MMT 17.8 °C
cumRR_data <- cumRR_data %>%
  mutate(
    HW_ssb1 = as.integer(cumRR_ssb1 > 1.1 & temp > 17.8),
    CTRL_ssb1 = as.integer(cumRR_ssb1 < 1.01 & temp >= 10.0 & temp <= 17.8)
  )

# Community dwellers: MMT 26.3 °C
cumRR_data <- cumRR_data %>%
  mutate(
    HW_ssb0 = as.integer(cumRR_ssb0 > 1.1 & temp > 26.3),
    CTRL_ssb0 = as.integer(cumRR_ssb0 < 1.01 & temp >= 10 & temp <= 26.3)
  )

cumRR_data <- cumRR_data %>%
  mutate(
    casecontrol_ssb1 = case_when(
      HW_ssb1 == 1   ~ 1L,
      CTRL_ssb1 == 1 ~ 0L,
      TRUE           ~ NA_integer_
    )
  )

cumRR_data <- cumRR_data %>%
  mutate(
    casecontrol_ssb0 = case_when(
      HW_ssb0 == 1   ~ 1L,
      CTRL_ssb0 == 1 ~ 0L,
      TRUE           ~ NA_integer_
    )
  )

cumRR_data <- cumRR_data %>%
  filter(!(is.na(casecontrol_ssb1) & is.na(casecontrol_ssb0)))

cumRR_data_all <- cumRR_data %>%
  select(date, temp, casecontrol_ssb1, casecontrol_ssb0)

save(cumRR_data_all, file = "data/cumRR_data_all.RData")
}

# 3. Study population: one valid death record, no emigration directly before death,
#    at least 1 year of residence in the region; linked to HW/OT days
{
  flytt <- read_sas("work/data_250307/multisjuka_flytt.sas7bdat")
  diseased_data <- flytt %>%
    filter(HKOD == "D") %>%
    group_by(lopnr) %>%
    summarize(count = n()) %>%
    filter(count == 1)

  diseased <- flytt %>%
    filter(lopnr %in% diseased_data$lopnr)

  diseased <- diseased %>%
    mutate(date = as.Date(HDAT, format = "%Y%m%d")) %>%
    arrange(lopnr, date) %>%
    select(-HDAT)

  na_dates <- diseased %>%
    filter(is.na(date)) %>%
    pull(lopnr) %>%
    unique()

  diseased <- diseased %>%
    filter(!lopnr %in% na_dates)

  # Remove individuals who emigrated (U) directly before death
  multiple_rows <- diseased %>%
    group_by(lopnr) %>%
    filter(n() > 1) %>%
    ungroup()

  lopnr_U_before_D <- multiple_rows %>%
    group_by(lopnr) %>%
    filter(any(HKOD == "U" & lead(HKOD) == "D")) %>%
    pull(lopnr) %>%
    unique()

  diseased <- diseased %>%
    filter(!lopnr %in% lopnr_U_before_D)

  # Residence: compare the latest immigration (I) date with the death date
  diseased <- diseased %>%
    arrange(lopnr, date) %>%
    mutate(
      I_date = ifelse(HKOD == "I", date, NA),
      D_date = ifelse(HKOD == "D", date, NA)
    ) %>%
    mutate(
      I_date = as.Date(I_date, origin = "1970-01-01"),
      D_date = as.Date(D_date, origin = "1970-01-01")
    )

  diseased <- diseased %>%
    filter(HKOD !="U")

  diseased <- diseased %>%
    group_by(lopnr) %>%
    mutate(
      D_criteria = ifelse(all(HKOD == "D"), TRUE, FALSE)) %>% ungroup()

  diseased1 <- diseased %>%
    filter(D_criteria == TRUE)

  diseased <- diseased %>%
    anti_join(diseased1, by = "lopnr")

  diseased <- diseased %>%
    arrange(lopnr, date) %>%
    group_by(lopnr) %>%
    slice_tail(n = 2) %>%
    ungroup()

  diseased <- diseased %>%
    group_by(lopnr) %>%
    mutate(
      I_date_lag = lag(I_date),
      D_date_lag = lag(D_date)
    ) %>%
    ungroup()

  diseased <- diseased %>%
    mutate(
      I_date_combined = ifelse(!is.na(I_date), I_date, I_date_lag),
      D_date_combined = ifelse(!is.na(D_date), D_date, D_date_lag)
    ) %>%
    mutate(
      I_date_combined = as.Date(I_date_combined, origin = "1970-01-01"),
      D_date_combined = as.Date(D_date_combined, origin = "1970-01-01")
    )

  diseased <- diseased %>%
    mutate(
      meets_1_years_criteria = ifelse((D_date_combined - I_date_combined) >= 365.25, TRUE, FALSE)
    )

  diseased <- diseased %>%
    filter(meets_1_years_criteria == TRUE)

  diseased <- diseased %>%
    select(lopnr, HKOD, date)

  diseased1 <- diseased1 %>%
    select(lopnr, HKOD, date)

  diseased_model_1 <- rbind(diseased, diseased1)

  # Link deaths to HW/OT days in each stratum
  load(file = "data/cumRR_data_all.RData")

  tempdata <- cumRR_data_all

  tempdata <- tempdata %>%
    mutate(date = as.Date(date))

  diseased_model_1 <- diseased_model_1 %>%
    mutate(date = as.Date(date))

  case.dates_ssb1 <- tempdata %>% filter(casecontrol_ssb1 == 1)
  control.dates_ssb1 <- tempdata %>% filter(casecontrol_ssb1 == 0)

  dead.case_ssb1 <- diseased_model_1 %>%
    filter(date %in% case.dates_ssb1$date) %>%
    left_join(tempdata, by = "date") %>%
    mutate(strata = "ssb1", case = 1)

  dead.control_ssb1 <- diseased_model_1 %>%
    filter(date %in% control.dates_ssb1$date) %>%
    left_join(tempdata, by = "date") %>%
    mutate(strata = "ssb1", case = 0)

  case.dates_ssb0 <- tempdata %>% filter(casecontrol_ssb0 == 1)
  control.dates_ssb0 <- tempdata %>% filter(casecontrol_ssb0 == 0)

  dead.case_ssb0 <- diseased_model_1 %>%
    filter(date %in% case.dates_ssb0$date) %>%
    left_join(tempdata, by = "date") %>%
    mutate(strata = "ssb0", case = 1)

  dead.control_ssb0 <- diseased_model_1 %>%
    filter(date %in% control.dates_ssb0$date) %>%
    left_join(tempdata, by = "date") %>%
    mutate(strata = "ssb0", case = 0)

  dead_model_1 <- bind_rows(
    dead.case_ssb1,
    dead.control_ssb1,
    dead.case_ssb0,
    dead.control_ssb0
  ) %>%
    group_by(lopnr) %>%
    summarise(
      death_date = first(date),
      case_ssb1 = if(all(is.na(casecontrol_ssb1))) NA_integer_ else max(casecontrol_ssb1, na.rm = TRUE),
      case_ssb0 = if(all(is.na(casecontrol_ssb0))) NA_integer_ else max(casecontrol_ssb0, na.rm = TRUE),
      .groups = "drop"
    )

  save(dead_model_1, file = "data/dead_model_1.RData")
  rm(list = ls())
  gc()
}

# 4. ICD-10 block grouping
{
  rank_2025 <- read_excel("documents/icd_group_se.xlsx")
  rank_2000 <- read_excel("documents/rank_2000.xlsx") %>%
    select(-PERCENT, -pvrank)

  unique_chapter <- rank_2025 %>%
    mutate(
      start_code = str_sub(chapter, 1, 3),
      end_code = str_sub(chapter, 5, 7)
    )

  match_chapter <- function(code, chapters) {
    for (i in 1:nrow(chapters)) {
      if (code >= chapters$start_code[i] && code <= chapters$end_code[i]) {
        return(chapters$chapter[i])
      }
    }
    return(NA)
  }

  rank_2000_2 <- rank_2000 %>%
    rowwise() %>%
    mutate(chapter = match_chapter(d3, unique_chapter)) %>%
    ungroup()

  icd2 <- rank_2000_2 %>%
    filter(!is.na(chapter))

  icd <- icd2 %>%
    group_by(chapter) %>%
    summarise(total_count = sum(COUNT)) %>%
    mutate(
      percentage = total_count / sum(total_count) * 100
    ) %>%
    arrange(desc(percentage)) %>%
    mutate(rank = row_number())

  # Expand a block range (e.g. A00-A09) into its three-character codes
  generate_codes <- function(chapter_range) {
    parts <- str_split(chapter_range, "-", simplify = TRUE)
    start <- parts[1]
    end <- parts[2]

    prefix_start <- substr(start, 1, 1)
    prefix_end <- substr(end, 1, 1)

    num_start <- as.integer(substr(start, 2, 3))
    num_end <- as.integer(substr(end, 2, 3))

    codes <- c()
    for (prefix in seq(from = utf8ToInt(prefix_start), to = utf8ToInt(prefix_end))) {
      current_prefix <- intToUtf8(prefix)
      start_num <- if (current_prefix == prefix_start) num_start else 0
      end_num <- if (current_prefix == prefix_end) num_end else 99

      for (i in start_num:end_num) {
        code <- paste0(current_prefix, str_pad(i, 2, pad = "0"))
        codes <- c(codes, code)
      }
    }

    paste(codes, collapse = "_")
  }

  icd <- icd %>%
    mutate(codes = sapply(chapter, generate_codes))

  icd <- icd %>%
    arrange(chapter)

  save(icd, file = "data/icd.RData")
}

# 5. Healthcare registers
# Inpatient care (SLV): diagnoses within 3 years before death
{
  load(file = "data/dead_model_1.RData")
  load(file = "data/icd.RData")
  rank_2000 <- icd
  slv <- read_sas("work/data_250307/multisjuka_slv.sas7bdat")

  slv <- slv %>%
    filter(lopnr %in% dead_model_1$lopnr) %>%
    mutate(UTDAT = ymd(UTDAT))

  slv <- slv %>%
    left_join(dead_model_1, slv, by = "lopnr") %>%
    relocate(death_date,  .after = "lopnr")

  slv <- slv %>%
    mutate(DIAG = paste(DIAG1,DIAG2,DIAG3,DIAG4,DIAG5,DIAG6,DIAG7,DIAG8,DIAG9,DIAG10, sep ="-")) %>%
    select(lopnr, KOEN, FMANAD, death_date, UTDAT, death_date, DIAG, AR, mosaicgrupp_sthlm)

  for (j in 1:nrow(rank_2000)) {
    chapter_name <- rank_2000$chapter[j]
    code_list <- unlist(strsplit(rank_2000$codes[j], "_"))

    slv[[chapter_name]] <- sapply(slv$DIAG, function(diag_string) {
      diag_codes <- unlist(strsplit(diag_string, "-"))
      diag_prefixes <- substr(diag_codes, 1, 3)
      any(diag_prefixes %in% code_list)
    }) * 1
  }

  slv_sum <- slv %>%
    group_by(lopnr) %>%
    filter((death_date - UTDAT) >= 0 & (death_date - UTDAT <= 365.25 * 3)) %>%
    summarise(across(all_of(rank_2000$chapter), max, .names = "{.col}"),
              across(c("KOEN","FMANAD", "death_date", "UTDAT", "mosaicgrupp_sthlm"), first)) %>%
    mutate(across("A00-A09":"Z80-Z99", ~ as.numeric(. > 0))) %>%
    relocate(KOEN, death_date, UTDAT, .after = lopnr)

  # Age at death (birth day set to the 15th of the birth month)
  slv_sum$FMANAD <- as.character(slv_sum$FMANAD)
  slv_sum$FMANAD <- paste0(substr(slv_sum$FMANAD, 1, 4), "-", substr(slv_sum$FMANAD, 5, 6), "-15")
  slv_sum$FMANAD <- ymd(slv_sum$FMANAD)
  slv_sum$death_age <- year(slv_sum$death_date) - year(slv_sum$FMANAD)
  slv_sum <- slv_sum %>%
    relocate(death_age, mosaicgrupp_sthlm, .after = death_date)

  slv_sum <- slv_sum %>%
    select(-"UTDAT", -"FMANAD")

  save(slv_sum, file = "data/slv_sum.RData")
  rm(list = ls())
  gc()
}

# Nursing home residence (HUSL): latest admission date per person
{
  load(file = "data/dead_model_1.RData")

  files <- list.files("work/data_250307", pattern = "^multisjuka_husl_\\d+\\.sas7bdat$", full.names = TRUE)

  husl_combined <- files %>%
    lapply(function(file_path) {
      message("Processing: ", basename(file_path))
      df <- read_sas(file_path)
      semi_join(df, dead_model_1, by = "lopnr")
    }) %>%
    bind_rows()

  husl_cleaned <- husl_combined %>%
    mutate(DATSSB = na_if(trimws(DATSSB), "")) %>%
    mutate(DATSSB = as.Date(DATSSB, format = "%Y%m%d")) %>%
    filter(!is.na(DATSSB)) %>%
    arrange(lopnr, DATSSB) %>%
    group_by(lopnr) %>%
    slice_tail(n = 1) %>%
    ungroup() %>%
    arrange(lopnr)

  husl_distinct_lopnr <- husl_combined %>%
    distinct(lopnr, .keep_all = TRUE)

  husl_remaining <- husl_distinct_lopnr %>%
    filter(!lopnr %in% husl_cleaned$lopnr)

  husl_remaining <- husl_remaining %>%
    mutate(DATSSB = na_if(trimws(DATSSB), "")) %>%
    mutate(DATSSB = as.Date(DATSSB, format = "%Y%m%d"))

  husl_all <- bind_rows(husl_remaining, husl_cleaned) %>%
    arrange(lopnr)

  husl_all <- left_join(husl_all, dead_model_1, by = "lopnr")

  husl_all$FMANAD <- as.character(husl_all$FMANAD)
  husl_all$FMANAD <- paste0(substr(husl_all$FMANAD, 1, 4), "-", substr(husl_all$FMANAD, 5, 6), "-15")
  husl_all$FMANAD <- ymd(husl_all$FMANAD)
  husl_all$death_age <- year(husl_all$death_date) - year(husl_all$FMANAD)
  husl_all <- husl_all %>%
    relocate(death_age, .after = death_date) %>%
    select(-AVDHL, -AVDSSB, -DATHL)

  save(husl_all, file = "data/husl_all.RData")
  rm(list = ls())
  gc()
}

# Primary care (KON), 2010-2013
{
  kon2010 <- read_sas("work/data_250307/multisjuka_kon_2010.sas7bdat")
  load(file = "data/dead_model_1.RData")
  kon2010 <- kon2010 %>% filter(lopnr %in% dead_model_1$lopnr)
  save(kon2010, file = "data/kon2010.RData")
  rm(list = ls())
  gc()

  kon2011 <- read_sas("work/data_250307/multisjuka_kon_2011.sas7bdat")
  load(file = "data/dead_model_1.RData")
  kon2011 <- kon2011 %>% filter(lopnr %in% dead_model_1$lopnr)
  save(kon2011, file = "data/kon2011.RData")
  rm(list = ls())
  gc()

  load(file = "data/dead_model_1.RData")
  kon2012 <- read_sas("work/data_250307/multisjuka_kon_2012.sas7bdat")
  kon2012 <- kon2012 %>%
    semi_join(dead_model_1, by = "lopnr")
  save(kon2012, file = "data/kon2012.RData")
  rm(list = ls())
  gc()

  load(file = "data/dead_model_1.RData")
  kon2013 <- read_sas("work/data_250307/multisjuka_kon_2013.sas7bdat")
  kon2013 <- kon2013 %>%
    semi_join(dead_model_1, by = "lopnr")
  save(kon2013, file = "data/kon2013.RData")
  rm(list = ls())
  gc()

  load(file = "data/kon2010.RData")
  load(file = "data/kon2011.RData")
  load(file = "data/kon2012.RData")
  load(file = "data/kon2013.RData")
  load(file = "data/dead_model_1.RData")
  load(file = "data/icd.RData")
  rank_2000 <- icd

  kon <- bind_rows(kon2010, kon2011, kon2012, kon2013) %>%
    arrange(lopnr)
  data.all <- kon %>% left_join(dead_model_1 %>% select(lopnr, death_date), by = "lopnr") %>%
    mutate(BDAT = ymd(BDAT)) %>%
    relocate(death_date, .after = BDAT)

  data.all.unique <- data.all %>%
    select(lopnr, KOEN, FMANAD, death_date, mosaicgrupp_sthlm) %>%
    distinct(lopnr, .keep_all = TRUE)

  data.all.unique$FMANAD <- as.character(data.all.unique$FMANAD)
  data.all.unique$FMANAD <- paste0(substr(data.all.unique$FMANAD, 1, 4), "-", substr(data.all.unique$FMANAD, 5, 6), "-15")
  data.all.unique$FMANAD <- ymd(data.all.unique$FMANAD)
  data.all.unique$death_age <- year(data.all.unique$death_date) - year(data.all.unique$FMANAD)
  data.all.unique <- data.all.unique %>%
    relocate(death_age, .after = death_date) %>%
    select(-FMANAD)

  # Care-giver category from VDG1: 1 = physician, 2 = nurse, 3 = assistant nurse, 4 = paramedical
  data.all <- data.all %>%
    mutate(DIAG = paste(DIAG1, DIAG2, DIAG3, DIAG4, DIAG5, DIAG6, DIAG7, DIAG8, sep = "-")) %>%
    select(lopnr, FMANAD, KOEN, BDAT, DIAG, VDG1, death_date, mosaicgrupp_sthlm) %>%
    mutate(VDG1_numeric = as.numeric(VDG1)) %>%
    mutate(VDG1_group = case_when(
      VDG1_numeric >= 1 & VDG1_numeric <= 69 ~ 1,
      VDG1_numeric %in% c(70, 89) ~ 2,
      VDG1_numeric == 91 ~ 3,
      VDG1_numeric %in% c(74, 75, 77, 78) ~ 4,
      TRUE ~ NA_real_
    ))

  data.all$FMANAD <- as.character(data.all$FMANAD)
  data.all$FMANAD <- paste0(substr(data.all$FMANAD, 1, 4), "-", substr(data.all$FMANAD, 5, 6), "-15")
  data.all$FMANAD <- ymd(data.all$FMANAD)
  data.all$death_age <- year(data.all$death_date) - year(data.all$FMANAD)
  data.all <- data.all %>%
    relocate(death_age, .after = death_date)

  # Diagnoses within 3 years before death, one binary variable per ICD-10 block
  data.all <- data.all %>%
    filter(death_date - BDAT <= years(3))

  for (j in 1:nrow(rank_2000)) {
    chapter_name <- rank_2000$chapter[j]

    code_list <- substr(unlist(strsplit(rank_2000$codes[j], "_")), 1, 3)

    data.all[[chapter_name]] <- sapply(data.all$DIAG, function(diag_string) {
      diag_codes <- strsplit(diag_string, "-")[[1]]
      diag_codes <- trimws(diag_codes[diag_codes != ""])
      diag_codes_truncated <- substr(diag_codes, 1, 3)
      any(diag_codes_truncated %in% code_list)
    }) * 1
  }

  summary_by_lopnr <- data.all %>%
    group_by(lopnr) %>%
    summarise(across(all_of(rank_2000$chapter), max, .names = "{.col}"),
              across(c(death_date, KOEN, death_age, mosaicgrupp_sthlm), first)
    ) %>%
    relocate(KOEN, death_date, death_age, mosaicgrupp_sthlm, .after = lopnr)

  # Visits within 1 year before death, by profession and time window (last month, months 2-12)
  data.all <- data.all %>%
    filter(VDG1_group %in% c(1, 2, 3, 4),
           (death_date - BDAT) <= (365.25 * 1)) %>%
    mutate(role = case_when(
      VDG1_group == 1 ~ "dr",
      VDG1_group == 2 ~ "nurse",
      VDG1_group == 3 ~ "as_nurse",
      VDG1_group == 4 ~ "param",
      TRUE ~ NA_character_
    ))

  data.all <- data.all %>%
    mutate(
      visits_year_1_dr = 0,
      visits_1_months_dr = 0,
      visits_year_1_nurse = 0,
      visits_1_months_nurse = 0,
      visits_year_1_asnurse = 0,
      visits_1_months_asnurse = 0,
      visits_year_1_param = 0,
      visits_1_months_param = 0,
    )

  data.all <- data.all %>%
    mutate(months_before_death = as.numeric((death_date - BDAT) / 30.44)) %>%
    relocate(months_before_death, .after = BDAT)

  data.all <- data.all %>%
    mutate(time_window = case_when(
      months_before_death <= 1 ~ "1_months",
      months_before_death > 1 & months_before_death <= 12 ~ "year_1",
      TRUE ~ NA_character_
    )) %>%
    relocate(time_window, .after = months_before_death)

  data.all <- data.all %>%
    relocate(role, .after = time_window)

  data.all <- data.all %>%
    mutate(
      visits_1_months_nurse = ifelse(role == "nurse" & time_window == "1_months", 1, 0),
      visits_3_months_nurse = ifelse(role == "nurse" & time_window == "3_months", 1, 0),
      visits_year_1_nurse = ifelse(role == "nurse" & time_window == "year_1", 1, 0),

      visits_1_months_dr = ifelse(role == "dr" & time_window == "1_months", 1, 0),
      visits_3_months_dr = ifelse(role == "dr" & time_window == "3_months", 1, 0),
      visits_year_1_dr = ifelse(role == "dr" & time_window == "year_1", 1, 0),

      visits_1_months_asnurse = ifelse(role == "as_nurse" & time_window == "1_months", 1, 0),
      visits_3_months_asnurse = ifelse(role == "as_nurse" & time_window == "3_months", 1, 0),
      visits_year_1_asnurse = ifelse(role == "as_nurse" & time_window == "year_1", 1, 0),

      visits_1_months_param = ifelse(role == "param" & time_window == "1_months", 1, 0),
      visits_3_months_param = ifelse(role == "param" & time_window == "3_months", 1, 0),
      visits_year_1_param = ifelse(role == "param" & time_window == "year_1", 1, 0)
    ) %>%
    group_by(lopnr) %>%
    summarise(
      across(starts_with("visits_"), sum, na.rm = TRUE),
      across(c(death_date, KOEN, death_age, mosaicgrupp_sthlm), first)
    ) %>%
    relocate(KOEN, death_date, death_age, mosaicgrupp_sthlm, .after = lopnr)

  missing_cols3 <- setdiff(names(data.all), names(summary_by_lopnr))
  for (col in missing_cols3) {
    summary_by_lopnr[[col]] <- 0
  }

  data.complete <- bind_rows(data.all, summary_by_lopnr) %>%
    arrange(lopnr) %>%
    group_by(lopnr) %>%
    summarise(
      # Visits from KON set to 0 to avoid double counting with OVR
      across(starts_with("visits"), ~ 0L),
      across(`A00-A09`:`Z80-Z99`, ~ sum(as.integer(.x > 0), na.rm = TRUE))
    )

  merged_data <- data.all.unique %>%
    left_join(data.complete, by = "lopnr")

  merged_data <- merged_data %>%
    mutate(across(setdiff(names(data.complete), "lopnr"), ~replace_na(., 0)))

  kon_complete <- merged_data

  save(kon_complete, file = "data/kon_complete.RData")

  rm(list = ls())
  gc()
}

# Outpatient care (OVR): restrict each year to the study population
{
  ovr2010 <- read_sas("work/data_250307/multisjuka_ovr_2010.sas7bdat")
  load(file = "data/dead_model_1.RData")
  ovr2010 <- ovr2010 %>%
    semi_join(dead_model_1, by = "lopnr")
  save(ovr2010, file = "data/ovr2010.RData")
  rm(list = ls())
  gc()

  ovr2011 <- read_sas("work/data_250307/multisjuka_ovr_2011.sas7bdat")
  load(file = "data/dead_model_1.RData")
  ovr2011 <- ovr2011 %>%
    semi_join(dead_model_1, by = "lopnr")
  save(ovr2011, file = "data/ovr2011.RData")
  rm(list = ls())
  gc()

  load(file = "data/dead_model_1.RData")
  ovr2012 <- read_sas("work/data_250307/multisjuka_ovr_2012.sas7bdat")
  ovr2012 <- ovr2012 %>%
    semi_join(dead_model_1, by = "lopnr")
  save(ovr2012, file = "data/ovr2012.RData")
  rm(list = ls())
  gc()

  load(file = "data/dead_model_1.RData")
  ovr2013 <- read_sas("work/data_250307/multisjuka_ovr_2013.sas7bdat")
  ovr2013 <- ovr2013 %>%
    semi_join(dead_model_1, by = "lopnr")
  save(ovr2013, file = "data/ovr2013.RData")
  rm(list = ls())
  gc()

  load(file = "data/dead_model_1.RData")
  ovr2014 <- read_sas("work/data_250307/multisjuka_ovr_2014.sas7bdat")
  ovr2014 <- ovr2014 %>%
    semi_join(dead_model_1, by = "lopnr")
  save(ovr2014, file = "data/ovr2014.RData")
  rm(list = ls())
  gc()
}

# OVR 2010-2012 (years processed in batches); same steps as for KON
{
  load(file = "data/ovr2010.RData")
  load(file = "data/ovr2011.RData")
  load(file = "data/ovr2012.RData")
  load(file = "data/dead_model_1.RData")
  load(file = "data/icd.RData")
  rank_2000 <- icd

  ovr <- bind_rows(ovr2010, ovr2011, ovr2012) %>%
    arrange(lopnr)

  data.all <- ovr %>% left_join(dead_model_1 %>% select(lopnr, death_date), by = "lopnr")

  data.all.unique <- data.all %>%
    filter(lopnr %in% dead_model_1$lopnr) %>%
    distinct(lopnr, .keep_all = TRUE)

  data.all.unique <- data.all.unique %>%
    select(lopnr, KOEN, FMANAD, death_date, mosaicgrupp_sthlm)

  data.all.unique$FMANAD <- as.character(data.all.unique$FMANAD)
  data.all.unique$FMANAD <- paste0(substr(data.all.unique$FMANAD, 1, 4), "-", substr(data.all.unique$FMANAD, 5, 6), "-15")
  data.all.unique$FMANAD <- ymd(data.all.unique$FMANAD)
  data.all.unique$death_age <- year(data.all.unique$death_date) - year(data.all.unique$FMANAD)
  data.all.unique <- data.all.unique %>%
    relocate(death_age, .after = death_date) %>%
    select(-FMANAD)

  data.all <- data.all %>%
    filter(lopnr %in% dead_model_1$lopnr) %>%
    mutate(BDAT = ymd(BDAT)) %>%
    relocate(death_date, .after = BDAT)

  data.all <- data.all %>%
    mutate(DIAG = paste(DIAG1, DIAG2, DIAG3, DIAG4, DIAG5, DIAG6, DIAG7, DIAG8, sep = "-")) %>%
    select(lopnr, FMANAD, KOEN, BDAT, DIAG, VDG1, death_date, mosaicgrupp_sthlm) %>%
    mutate(VDG1_numeric = as.numeric(VDG1)) %>%
    mutate(VDG1_group = case_when(
      VDG1_numeric >= 1 & VDG1_numeric <= 69 ~ 1,
      VDG1_numeric %in% c(70, 89) ~ 2,
      VDG1_numeric == 91 ~ 3,
      VDG1_numeric %in% c(74, 75, 77, 78) ~ 4,
      TRUE ~ NA_real_
    ))

  data.all$FMANAD <- as.character(data.all$FMANAD)
  data.all$FMANAD <- paste0(substr(data.all$FMANAD, 1, 4), "-", substr(data.all$FMANAD, 5, 6), "-15")
  data.all$FMANAD <- as.Date(data.all$FMANAD)
  data.all$death_date <- as.Date(data.all$death_date)
  data.all$death_age <- year(data.all$death_date) - year(data.all$FMANAD)
  data.all <- data.all %>%
    relocate(death_age, .after = death_date)

  data.all <- data.all %>%
    filter(death_date - BDAT <= years(3))

  for (j in 1:nrow(rank_2000)) {
    chapter_name <- rank_2000$chapter[j]

    code_list <- substr(unlist(strsplit(rank_2000$codes[j], "_")), 1, 3)

    data.all[[chapter_name]] <- sapply(data.all$DIAG, function(diag_string) {
      diag_codes <- strsplit(diag_string, "-")[[1]]
      diag_codes <- trimws(diag_codes[diag_codes != ""])
      diag_codes_truncated <- substr(diag_codes, 1, 3)
      any(diag_codes_truncated %in% code_list)
    }) * 1
  }

  summary_by_lopnr <- data.all %>%
    group_by(lopnr) %>%
    summarise(across(all_of(rank_2000$chapter), max, .names = "{.col}"),
              across(c(death_date, KOEN, death_age, mosaicgrupp_sthlm), first)
    ) %>%
    relocate(KOEN, death_date, death_age, mosaicgrupp_sthlm, .after = lopnr)

  data.all <- data.all %>%
    filter(VDG1_group %in% c(1, 2, 3, 4),
           (death_date - BDAT) <= (365.25 * 1)) %>%
    mutate(role = case_when(
      VDG1_group == 1 ~ "dr",
      VDG1_group == 2 ~ "nurse",
      VDG1_group == 3 ~ "as_nurse",
      VDG1_group == 4 ~ "param",
      TRUE ~ NA_character_
    ))

  data.all <- data.all %>%
    mutate(
      visits_year_1_dr = 0,
      visits_1_months_dr = 0,
      visits_year_1_nurse = 0,
      visits_1_months_nurse = 0,
      visits_year_1_asnurse = 0,
      visits_1_months_asnurse = 0,
      visits_year_1_param = 0,
      visits_1_months_param = 0,
    )

  data.all <- data.all %>%
    mutate(months_before_death = as.numeric((death_date - BDAT) / 30.44)) %>%
    relocate(months_before_death, .after = BDAT)

  data.all <- data.all %>%
    mutate(time_window = case_when(
      months_before_death <= 1 ~ "1_months",
      months_before_death > 1 & months_before_death <= 12 ~ "year_1",
      TRUE ~ NA_character_
    )) %>%
    relocate(time_window, .after = months_before_death)

  data.all <- data.all %>%
    relocate(role, .after = time_window)

  data.all <- data.all %>%
    mutate(
      visits_1_months_nurse = ifelse(role == "nurse" & time_window == "1_months", 1, 0),
      visits_3_months_nurse = ifelse(role == "nurse" & time_window == "3_months", 1, 0),
      visits_year_1_nurse = ifelse(role == "nurse" & time_window == "year_1", 1, 0),

      visits_1_months_dr = ifelse(role == "dr" & time_window == "1_months", 1, 0),
      visits_3_months_dr = ifelse(role == "dr" & time_window == "3_months", 1, 0),
      visits_year_1_dr = ifelse(role == "dr" & time_window == "year_1", 1, 0),

      visits_1_months_asnurse = ifelse(role == "as_nurse" & time_window == "1_months", 1, 0),
      visits_3_months_asnurse = ifelse(role == "as_nurse" & time_window == "3_months", 1, 0),
      visits_year_1_asnurse = ifelse(role == "as_nurse" & time_window == "year_1", 1, 0),

      visits_1_months_param = ifelse(role == "param" & time_window == "1_months", 1, 0),
      visits_3_months_param = ifelse(role == "param" & time_window == "3_months", 1, 0),
      visits_year_1_param = ifelse(role == "param" & time_window == "year_1", 1, 0)
    ) %>%
    group_by(lopnr) %>%
    summarise(
      across(starts_with("visits_"), sum, na.rm = TRUE),
      across(c(death_date, KOEN, death_age, mosaicgrupp_sthlm), first)
    ) %>%
    relocate(KOEN, death_date, death_age, mosaicgrupp_sthlm, .after = lopnr)

  missing_cols3 <- setdiff(names(data.all), names(summary_by_lopnr))
  for (col in missing_cols3) {
    summary_by_lopnr[[col]] <- 0
  }

  data.complete <- bind_rows(data.all, summary_by_lopnr) %>%
    arrange(lopnr) %>%
    group_by(lopnr) %>%
    summarise(
      across(starts_with("visits_"), ~ sum(.x, na.rm = TRUE)),
      across(`A00-A09`:`Z80-Z99`, ~ sum(as.numeric(.x > 0), na.rm = TRUE)))

  merged_data_ovr1 <- data.all.unique %>%
    left_join(data.complete, by = "lopnr")

  merged_data_ovr1 <- merged_data_ovr1 %>%
    mutate(across(setdiff(names(data.complete), "lopnr"), ~replace_na(., 0)))

  save(merged_data_ovr1, file = "data/merged_data_ovr1.RData")

  rm(list = ls())
  gc()
}

# OVR 2013-2014
{
  load(file = "data/ovr2013.RData")
  load(file = "data/ovr2014.RData")
  load(file = "data/dead_model_1.RData")
  load(file = "data/icd.RData")
  rank_2000 <- icd

  ovr <- bind_rows(ovr2013, ovr2014) %>%
    arrange(lopnr)

  data.all <- ovr %>% left_join(dead_model_1 %>% select(lopnr, death_date), by = "lopnr")

  data.all.unique <- data.all %>%
    filter(lopnr %in% dead_model_1$lopnr) %>%
    distinct(lopnr, .keep_all = TRUE)

  data.all.unique <- data.all.unique %>%
    select(lopnr, KOEN, FMANAD, death_date, mosaicgrupp_sthlm)

  data.all.unique$FMANAD <- as.character(data.all.unique$FMANAD)
  data.all.unique$FMANAD <- paste0(substr(data.all.unique$FMANAD, 1, 4), "-", substr(data.all.unique$FMANAD, 5, 6), "-15")
  data.all.unique$FMANAD <- ymd(data.all.unique$FMANAD)
  data.all.unique$death_age <- year(data.all.unique$death_date) - year(data.all.unique$FMANAD)
  data.all.unique <- data.all.unique %>%
    relocate(death_age, .after = death_date) %>%
    select(-FMANAD)

  data.all <- data.all %>%
    filter(lopnr %in% dead_model_1$lopnr) %>%
    mutate(BDAT = ymd(BDAT)) %>%
    relocate(death_date, .after = BDAT)

  data.all <- data.all %>%
    mutate(DIAG = paste(DIAG1, DIAG2, DIAG3, DIAG4, DIAG5, DIAG6, DIAG7, DIAG8, sep = "-")) %>%
    select(lopnr, FMANAD, KOEN, BDAT, DIAG, VDG1, death_date, mosaicgrupp_sthlm) %>%
    mutate(VDG1_numeric = as.numeric(VDG1)) %>%
    mutate(VDG1_group = case_when(
      VDG1_numeric >= 1 & VDG1_numeric <= 69 ~ 1,
      VDG1_numeric %in% c(70, 89) ~ 2,
      VDG1_numeric == 91 ~ 3,
      VDG1_numeric %in% c(74, 75, 77, 78) ~ 4,
      TRUE ~ NA_real_
    ))

  data.all$FMANAD <- as.character(data.all$FMANAD)
  data.all$FMANAD <- paste0(substr(data.all$FMANAD, 1, 4), "-", substr(data.all$FMANAD, 5, 6), "-15")
  data.all$FMANAD <- as.Date(data.all$FMANAD)
  data.all$death_date <- as.Date(data.all$death_date)
  data.all$death_age <- year(data.all$death_date) - year(data.all$FMANAD)
  data.all <- data.all %>%
    relocate(death_age, .after = death_date)

  data.all <- data.all %>%
    filter(death_date - BDAT <= years(3))

  for (j in 1:nrow(rank_2000)) {
    chapter_name <- rank_2000$chapter[j]

    code_list <- substr(unlist(strsplit(rank_2000$codes[j], "_")), 1, 3)

    data.all[[chapter_name]] <- sapply(data.all$DIAG, function(diag_string) {
      diag_codes <- strsplit(diag_string, "-")[[1]]
      diag_codes <- trimws(diag_codes[diag_codes != ""])
      diag_codes_truncated <- substr(diag_codes, 1, 3)
      any(diag_codes_truncated %in% code_list)
    }) * 1
  }

  summary_by_lopnr <- data.all %>%
    group_by(lopnr) %>%
    summarise(across(all_of(rank_2000$chapter), max, .names = "{.col}"),
              across(c(death_date, KOEN, death_age, mosaicgrupp_sthlm), first)
    ) %>%
    relocate(KOEN, death_date, death_age, mosaicgrupp_sthlm, .after = lopnr)

  data.all <- data.all %>%
    filter(VDG1_group %in% c(1, 2, 3, 4),
           (death_date - BDAT) <= (365.25 * 1)) %>%
    mutate(role = case_when(
      VDG1_group == 1 ~ "dr",
      VDG1_group == 2 ~ "nurse",
      VDG1_group == 3 ~ "as_nurse",
      VDG1_group == 4 ~ "param",
      TRUE ~ NA_character_
    ))

  data.all <- data.all %>%
    mutate(
      visits_year_1_dr = 0,
      visits_1_months_dr = 0,
      visits_year_1_nurse = 0,
      visits_1_months_nurse = 0,
      visits_year_1_asnurse = 0,
      visits_1_months_asnurse = 0,
      visits_year_1_param = 0,
      visits_1_months_param = 0,
    )

  data.all <- data.all %>%
    mutate(months_before_death = as.numeric((death_date - BDAT) / 30.44)) %>%
    relocate(months_before_death, .after = BDAT)

  data.all <- data.all %>%
    mutate(time_window = case_when(
      months_before_death <= 1 ~ "1_months",
      months_before_death > 1 & months_before_death <= 12 ~ "year_1",
      TRUE ~ NA_character_
    )) %>%
    relocate(time_window, .after = months_before_death)

  data.all <- data.all %>%
    relocate(role, .after = time_window)

  data.all <- data.all %>%
    mutate(
      visits_1_months_nurse = ifelse(role == "nurse" & time_window == "1_months", 1, 0),
      visits_3_months_nurse = ifelse(role == "nurse" & time_window == "3_months", 1, 0),
      visits_year_1_nurse = ifelse(role == "nurse" & time_window == "year_1", 1, 0),

      visits_1_months_dr = ifelse(role == "dr" & time_window == "1_months", 1, 0),
      visits_3_months_dr = ifelse(role == "dr" & time_window == "3_months", 1, 0),
      visits_year_1_dr = ifelse(role == "dr" & time_window == "year_1", 1, 0),

      visits_1_months_asnurse = ifelse(role == "as_nurse" & time_window == "1_months", 1, 0),
      visits_3_months_asnurse = ifelse(role == "as_nurse" & time_window == "3_months", 1, 0),
      visits_year_1_asnurse = ifelse(role == "as_nurse" & time_window == "year_1", 1, 0),

      visits_1_months_param = ifelse(role == "param" & time_window == "1_months", 1, 0),
      visits_3_months_param = ifelse(role == "param" & time_window == "3_months", 1, 0),
      visits_year_1_param = ifelse(role == "param" & time_window == "year_1", 1, 0)
    ) %>%
    group_by(lopnr) %>%
    summarise(
      across(starts_with("visits_"), sum, na.rm = TRUE),
      across(c(death_date, KOEN, death_age, mosaicgrupp_sthlm), first)
    ) %>%
    relocate(KOEN, death_date, death_age, mosaicgrupp_sthlm, .after = lopnr)

  missing_cols3 <- setdiff(names(data.all), names(summary_by_lopnr))
  for (col in missing_cols3) {
    summary_by_lopnr[[col]] <- 0
  }

  data.complete <- bind_rows(data.all, summary_by_lopnr) %>%
    arrange(lopnr) %>%
    group_by(lopnr) %>%
    summarise(
      across(starts_with("visits_"), ~ sum(.x, na.rm = TRUE)),
      across(`A00-A09`:`Z80-Z99`, ~ sum(as.numeric(.x > 0), na.rm = TRUE)))

  merged_data_ovr1b <- data.all.unique %>%
    left_join(data.complete, by = "lopnr")

  merged_data_ovr1b <- merged_data_ovr1b %>%
    mutate(across(setdiff(names(data.complete), "lopnr"), ~replace_na(., 0)))

  save(merged_data_ovr1b, file = "data/merged_data_ovr1b.RData")

  rm(list = ls())
  gc()
}

# OVR 2015-2019: restrict to the study population
{
  load(file = "data/dead_model_1.RData")
  ovr2015 <- read_sas("work/data_250307/multisjuka_ovr_2015.sas7bdat")
  ovr2015 <- ovr2015 %>%
    semi_join(dead_model_1, by = "lopnr")
  save(ovr2015, file = "data/ovr2015.RData")
  rm(list = ls())
  gc()

  load(file = "data/dead_model_1.RData")
  ovr2016 <- read_sas("work/data_250307/multisjuka_ovr_2016.sas7bdat")
  ovr2016 <- ovr2016 %>%
    semi_join(dead_model_1, by = "lopnr")
  save(ovr2016, file = "data/ovr2016.RData")
  rm(list = ls())
  gc()

  load(file = "data/dead_model_1.RData")
  ovr2017 <- read_sas("work/data_250307/multisjuka_ovr_2017.sas7bdat")
  ovr2017 <- ovr2017 %>%
    semi_join(dead_model_1, by = "lopnr")
  save(ovr2017, file = "data/ovr2017.RData")
  rm(list = ls())
  gc()

  load(file = "data/dead_model_1.RData")
  ovr2018 <- read_sas("work/data_250307/multisjuka_ovr_2018.sas7bdat")
  ovr2018 <- ovr2018 %>%
    semi_join(dead_model_1, by = "lopnr")
  save(ovr2018, file = "data/ovr2018.RData")
  rm(list = ls())
  gc()

  load(file = "data/dead_model_1.RData")
  ovr2019 <- read_sas("work/data_250307/multisjuka_ovr_2019.sas7bdat")
  ovr2019 <- ovr2019 %>%
    semi_join(dead_model_1, by = "lopnr")
  save(ovr2019, file = "data/ovr2019.RData")
  rm(list = ls())
  gc()
}

# OVR 2015-2017
{
  load(file = "data/ovr2015.RData")
  load(file = "data/ovr2016.RData")
  load(file = "data/ovr2017.RData")
  load(file = "data/dead_model_1.RData")
  load(file = "data/icd.RData")
  rank_2000 <- icd

  ovr <- bind_rows(ovr2015, ovr2016, ovr2017) %>%
    arrange(lopnr)

  data.all <- ovr %>% left_join(dead_model_1 %>% select(lopnr, death_date), by = "lopnr")

  data.all.unique <- data.all %>%
    filter(lopnr %in% dead_model_1$lopnr) %>%
    distinct(lopnr, .keep_all = TRUE)

  data.all.unique <- data.all.unique %>%
    select(lopnr, KOEN, FMANAD, death_date, mosaicgrupp_sthlm)

  data.all.unique$FMANAD <- as.character(data.all.unique$FMANAD)
  data.all.unique$FMANAD <- paste0(substr(data.all.unique$FMANAD, 1, 4), "-", substr(data.all.unique$FMANAD, 5, 6), "-15")
  data.all.unique$FMANAD <- ymd(data.all.unique$FMANAD)
  data.all.unique$death_age <- year(data.all.unique$death_date) - year(data.all.unique$FMANAD)
  data.all.unique <- data.all.unique %>%
    relocate(death_age, .after = death_date) %>%
    select(-FMANAD)

  data.all <- data.all %>%
    filter(lopnr %in% dead_model_1$lopnr) %>%
    mutate(BDAT = ymd(BDAT)) %>%
    relocate(death_date, .after = BDAT)

  data.all <- data.all %>%
    mutate(DIAG = paste(DIAG1, DIAG2, DIAG3, DIAG4, DIAG5, DIAG6, DIAG7, DIAG8, sep = "-")) %>%
    select(lopnr, FMANAD, KOEN, BDAT, DIAG, VDG1, death_date, mosaicgrupp_sthlm) %>%
    mutate(VDG1_numeric = as.numeric(VDG1)) %>%
    mutate(VDG1_group = case_when(
      VDG1_numeric >= 1 & VDG1_numeric <= 69 ~ 1,
      VDG1_numeric %in% c(70, 89) ~ 2,
      VDG1_numeric == 91 ~ 3,
      VDG1_numeric %in% c(74, 75, 77, 78) ~ 4,
      TRUE ~ NA_real_
    ))

  data.all$FMANAD <- as.character(data.all$FMANAD)
  data.all$FMANAD <- paste0(substr(data.all$FMANAD, 1, 4), "-", substr(data.all$FMANAD, 5, 6), "-15")
  data.all$FMANAD <- as.Date(data.all$FMANAD)
  data.all$death_date <- as.Date(data.all$death_date)
  data.all$death_age <- year(data.all$death_date) - year(data.all$FMANAD)
  data.all <- data.all %>%
    relocate(death_age, .after = death_date)

  data.all <- data.all %>%
    filter(death_date - BDAT <= years(3))

  for (j in 1:nrow(rank_2000)) {
    chapter_name <- rank_2000$chapter[j]

    code_list <- substr(unlist(strsplit(rank_2000$codes[j], "_")), 1, 3)

    data.all[[chapter_name]] <- sapply(data.all$DIAG, function(diag_string) {
      diag_codes <- strsplit(diag_string, "-")[[1]]
      diag_codes <- trimws(diag_codes[diag_codes != ""])
      diag_codes_truncated <- substr(diag_codes, 1, 3)
      any(diag_codes_truncated %in% code_list)
    }) * 1
  }

  summary_by_lopnr <- data.all %>%
    group_by(lopnr) %>%
    summarise(across(all_of(rank_2000$chapter), max, .names = "{.col}"),
              across(c(death_date, KOEN, death_age, mosaicgrupp_sthlm), first)
    ) %>%
    relocate(KOEN, death_date, death_age, mosaicgrupp_sthlm, .after = lopnr)

  data.all <- data.all %>%
    filter(VDG1_group %in% c(1, 2, 3, 4),
           (death_date - BDAT) <= (365.25 * 1)) %>%
    mutate(role = case_when(
      VDG1_group == 1 ~ "dr",
      VDG1_group == 2 ~ "nurse",
      VDG1_group == 3 ~ "as_nurse",
      VDG1_group == 4 ~ "param",
      TRUE ~ NA_character_
    ))

  data.all <- data.all %>%
    mutate(
      visits_year_1_dr = 0,
      visits_1_months_dr = 0,
      visits_year_1_nurse = 0,
      visits_1_months_nurse = 0,
      visits_year_1_asnurse = 0,
      visits_1_months_asnurse = 0,
      visits_year_1_param = 0,
      visits_1_months_param = 0,
    )

  data.all <- data.all %>%
    mutate(months_before_death = as.numeric((death_date - BDAT) / 30.44)) %>%
    relocate(months_before_death, .after = BDAT)

  data.all <- data.all %>%
    mutate(time_window = case_when(
      months_before_death <= 1 ~ "1_months",
      months_before_death > 1 & months_before_death <= 12 ~ "year_1",
      TRUE ~ NA_character_
    )) %>%
    relocate(time_window, .after = months_before_death)

  data.all <- data.all %>%
    relocate(role, .after = time_window)

  data.all <- data.all %>%
    mutate(
      visits_1_months_nurse = ifelse(role == "nurse" & time_window == "1_months", 1, 0),
      visits_3_months_nurse = ifelse(role == "nurse" & time_window == "3_months", 1, 0),
      visits_year_1_nurse = ifelse(role == "nurse" & time_window == "year_1", 1, 0),

      visits_1_months_dr = ifelse(role == "dr" & time_window == "1_months", 1, 0),
      visits_3_months_dr = ifelse(role == "dr" & time_window == "3_months", 1, 0),
      visits_year_1_dr = ifelse(role == "dr" & time_window == "year_1", 1, 0),

      visits_1_months_asnurse = ifelse(role == "as_nurse" & time_window == "1_months", 1, 0),
      visits_3_months_asnurse = ifelse(role == "as_nurse" & time_window == "3_months", 1, 0),
      visits_year_1_asnurse = ifelse(role == "as_nurse" & time_window == "year_1", 1, 0),

      visits_1_months_param = ifelse(role == "param" & time_window == "1_months", 1, 0),
      visits_3_months_param = ifelse(role == "param" & time_window == "3_months", 1, 0),
      visits_year_1_param = ifelse(role == "param" & time_window == "year_1", 1, 0)
    ) %>%
    group_by(lopnr) %>%
    summarise(
      across(starts_with("visits_"), sum, na.rm = TRUE),
      across(c(death_date, KOEN, death_age, mosaicgrupp_sthlm), first)
    ) %>%
    relocate(KOEN, death_date, death_age, mosaicgrupp_sthlm, .after = lopnr)

  missing_cols3 <- setdiff(names(data.all), names(summary_by_lopnr))
  for (col in missing_cols3) {
    summary_by_lopnr[[col]] <- 0
  }

  data.complete <- bind_rows(data.all, summary_by_lopnr) %>%
    arrange(lopnr) %>%
    group_by(lopnr) %>%
    summarise(
      across(starts_with("visits_"), ~ sum(.x, na.rm = TRUE)),
      across(`A00-A09`:`Z80-Z99`, ~ sum(as.numeric(.x > 0), na.rm = TRUE)))

  merged_data_ovr2 <- data.all.unique %>%
    left_join(data.complete, by = "lopnr")

  merged_data_ovr2 <- merged_data_ovr2 %>%
    mutate(across(setdiff(names(data.complete), "lopnr"), ~replace_na(., 0)))

  save(merged_data_ovr2, file = "data/merged_data_ovr2.RData")

  rm(list = ls())
  gc()
}

# OVR 2018-2019
{
  load(file = "data/ovr2018.RData")
  load(file = "data/ovr2019.RData")
  load(file = "data/dead_model_1.RData")
  load(file = "data/icd.RData")
  rank_2000 <- icd

  ovr <- bind_rows(ovr2018, ovr2019) %>%
    arrange(lopnr)

  data.all <- ovr %>% left_join(dead_model_1 %>% select(lopnr, death_date), by = "lopnr")

  data.all.unique <- data.all %>%
    filter(lopnr %in% dead_model_1$lopnr) %>%
    distinct(lopnr, .keep_all = TRUE)

  data.all.unique <- data.all.unique %>%
    select(lopnr, KOEN, FMANAD, death_date, mosaicgrupp_sthlm)

  data.all.unique$FMANAD <- as.character(data.all.unique$FMANAD)
  data.all.unique$FMANAD <- paste0(substr(data.all.unique$FMANAD, 1, 4), "-", substr(data.all.unique$FMANAD, 5, 6), "-15")
  data.all.unique$FMANAD <- ymd(data.all.unique$FMANAD)
  data.all.unique$death_age <- year(data.all.unique$death_date) - year(data.all.unique$FMANAD)
  data.all.unique <- data.all.unique %>%
    relocate(death_age, .after = death_date) %>%
    select(-FMANAD)

  data.all <- data.all %>%
    filter(lopnr %in% dead_model_1$lopnr) %>%
    mutate(BDAT = ymd(BDAT)) %>%
    relocate(death_date, .after = BDAT)

  data.all <- data.all %>%
    mutate(DIAG = paste(DIAG1, DIAG2, DIAG3, DIAG4, DIAG5, DIAG6, DIAG7, DIAG8, sep = "-")) %>%
    select(lopnr, FMANAD, KOEN, BDAT, DIAG, VDG1, death_date, mosaicgrupp_sthlm) %>%
    mutate(VDG1_numeric = as.numeric(VDG1)) %>%
    mutate(VDG1_group = case_when(
      VDG1_numeric >= 1 & VDG1_numeric <= 69 ~ 1,
      VDG1_numeric %in% c(70, 89) ~ 2,
      VDG1_numeric == 91 ~ 3,
      VDG1_numeric %in% c(74, 75, 77, 78) ~ 4,
      TRUE ~ NA_real_
    ))

  data.all$FMANAD <- as.character(data.all$FMANAD)
  data.all$FMANAD <- paste0(substr(data.all$FMANAD, 1, 4), "-", substr(data.all$FMANAD, 5, 6), "-15")
  data.all$FMANAD <- as.Date(data.all$FMANAD)
  data.all$death_date <- as.Date(data.all$death_date)
  data.all$death_age <- year(data.all$death_date) - year(data.all$FMANAD)
  data.all <- data.all %>%
    relocate(death_age, .after = death_date)

  data.all <- data.all %>%
    filter(death_date - BDAT <= years(3))

  for (j in 1:nrow(rank_2000)) {
    chapter_name <- rank_2000$chapter[j]

    code_list <- substr(unlist(strsplit(rank_2000$codes[j], "_")), 1, 3)

    data.all[[chapter_name]] <- sapply(data.all$DIAG, function(diag_string) {
      diag_codes <- strsplit(diag_string, "-")[[1]]
      diag_codes <- trimws(diag_codes[diag_codes != ""])
      diag_codes_truncated <- substr(diag_codes, 1, 3)
      any(diag_codes_truncated %in% code_list)
    }) * 1
  }

  summary_by_lopnr <- data.all %>%
    group_by(lopnr) %>%
    summarise(across(all_of(rank_2000$chapter), max, .names = "{.col}"),
              across(c(death_date, KOEN, death_age, mosaicgrupp_sthlm), first)
    ) %>%
    relocate(KOEN, death_date, death_age, mosaicgrupp_sthlm, .after = lopnr)

  data.all <- data.all %>%
    filter(VDG1_group %in% c(1, 2, 3, 4),
           (death_date - BDAT) <= (365.25 * 1)) %>%
    mutate(role = case_when(
      VDG1_group == 1 ~ "dr",
      VDG1_group == 2 ~ "nurse",
      VDG1_group == 3 ~ "as_nurse",
      VDG1_group == 4 ~ "param",
      TRUE ~ NA_character_
    ))

  data.all <- data.all %>%
    mutate(
      visits_year_1_dr = 0,
      visits_1_months_dr = 0,
      visits_year_1_nurse = 0,
      visits_1_months_nurse = 0,
      visits_year_1_asnurse = 0,
      visits_1_months_asnurse = 0,
      visits_year_1_param = 0,
      visits_1_months_param = 0,
    )

  data.all <- data.all %>%
    mutate(months_before_death = as.numeric((death_date - BDAT) / 30.44)) %>%
    relocate(months_before_death, .after = BDAT)

  data.all <- data.all %>%
    mutate(time_window = case_when(
      months_before_death <= 1 ~ "1_months",
      months_before_death > 1 & months_before_death <= 12 ~ "year_1",
      TRUE ~ NA_character_
    )) %>%
    relocate(time_window, .after = months_before_death)

  data.all <- data.all %>%
    relocate(role, .after = time_window)

  data.all <- data.all %>%
    mutate(
      visits_1_months_nurse = ifelse(role == "nurse" & time_window == "1_months", 1, 0),
      visits_3_months_nurse = ifelse(role == "nurse" & time_window == "3_months", 1, 0),
      visits_year_1_nurse = ifelse(role == "nurse" & time_window == "year_1", 1, 0),

      visits_1_months_dr = ifelse(role == "dr" & time_window == "1_months", 1, 0),
      visits_3_months_dr = ifelse(role == "dr" & time_window == "3_months", 1, 0),
      visits_year_1_dr = ifelse(role == "dr" & time_window == "year_1", 1, 0),

      visits_1_months_asnurse = ifelse(role == "as_nurse" & time_window == "1_months", 1, 0),
      visits_3_months_asnurse = ifelse(role == "as_nurse" & time_window == "3_months", 1, 0),
      visits_year_1_asnurse = ifelse(role == "as_nurse" & time_window == "year_1", 1, 0),

      visits_1_months_param = ifelse(role == "param" & time_window == "1_months", 1, 0),
      visits_3_months_param = ifelse(role == "param" & time_window == "3_months", 1, 0),
      visits_year_1_param = ifelse(role == "param" & time_window == "year_1", 1, 0)
    ) %>%
    group_by(lopnr) %>%
    summarise(
      across(starts_with("visits_"), sum, na.rm = TRUE),
      across(c(death_date, KOEN, death_age, mosaicgrupp_sthlm), first)
    ) %>%
    relocate(KOEN, death_date, death_age, mosaicgrupp_sthlm, .after = lopnr)

  missing_cols3 <- setdiff(names(data.all), names(summary_by_lopnr))
  for (col in missing_cols3) {
    summary_by_lopnr[[col]] <- 0
  }

  data.complete <- bind_rows(data.all, summary_by_lopnr) %>%
    arrange(lopnr) %>%
    group_by(lopnr) %>%
    summarise(
      across(starts_with("visits_"), ~ sum(.x, na.rm = TRUE)),
      across(`A00-A09`:`Z80-Z99`, ~ sum(as.numeric(.x > 0), na.rm = TRUE)))

  merged_data_ovr2b <- data.all.unique %>%
    left_join(data.complete, by = "lopnr")

  merged_data_ovr2b <- merged_data_ovr2b %>%
    mutate(across(setdiff(names(data.complete), "lopnr"), ~replace_na(., 0)))

  save(merged_data_ovr2b, file = "data/merged_data_ovr2b.RData")

  rm(list = ls())
  gc()
}

# OVR 2020-2024: restrict to the study population, then process 2020-2022
{
  {
    load(file = "data/dead_model_1.RData")
    ovr2020 <- read_sas("work/data_250307/multisjuka_ovr_2020.sas7bdat")
    ovr2020 <- ovr2020 %>%
      semi_join(dead_model_1, by = "lopnr")
    save(ovr2020, file = "data/ovr2020.RData")
    rm(list = ls())
    gc()

    load(file = "data/dead_model_1.RData")
    ovr2021 <- read_sas("work/data_250307/multisjuka_ovr_2021.sas7bdat")
    ovr2021 <- ovr2021 %>%
      semi_join(dead_model_1, by = "lopnr")
    save(ovr2021, file = "data/ovr2021.RData")
    rm(list = ls())
    gc()

    load(file = "data/dead_model_1.RData")
    ovr2022 <- read_sas("work/data_250307/multisjuka_ovr_2022.sas7bdat")
    ovr2022 <- ovr2022 %>%
      semi_join(dead_model_1, by = "lopnr")
    save(ovr2022, file = "data/ovr2022.RData")
    rm(list = ls())
    gc()

    load(file = "data/dead_model_1.RData")
    ovr2023 <- read_sas("work/data_250307/multisjuka_ovr_2023.sas7bdat")
    ovr2023 <- ovr2023 %>%
      semi_join(dead_model_1, by = "lopnr")
    save(ovr2023, file = "data/ovr2023.RData")
    rm(list = ls())
    gc()

    load(file = "data/dead_model_1.RData")
    ovr2024 <- read_sas("work/data_250307/multisjuka_ovr_2024.sas7bdat")
    ovr2024 <- ovr2024 %>%
      semi_join(dead_model_1, by = "lopnr")
    save(ovr2024, file = "data/ovr2024.RData")
    rm(list = ls())
    gc()
  }
  load(file = "data/ovr2020.RData")
  load(file = "data/ovr2021.RData")
  load(file = "data/ovr2022.RData")
  load(file = "data/dead_model_1.RData")
  load(file = "data/icd.RData")
  rank_2000 <- icd

  ovr <- bind_rows(ovr2020, ovr2021, ovr2022) %>%
    arrange(lopnr)
  data.all <- ovr %>% left_join(dead_model_1 %>% select(lopnr, death_date), by = "lopnr")

  data.all.unique <- data.all %>%
    filter(lopnr %in% dead_model_1$lopnr) %>%
    distinct(lopnr, .keep_all = TRUE)

  data.all.unique <- data.all.unique %>%
    select(lopnr, KOEN, FMANAD, death_date, mosaicgrupp_sthlm)

  data.all.unique$FMANAD <- as.character(data.all.unique$FMANAD)
  data.all.unique$FMANAD <- paste0(substr(data.all.unique$FMANAD, 1, 4), "-", substr(data.all.unique$FMANAD, 5, 6), "-15")
  data.all.unique$FMANAD <- ymd(data.all.unique$FMANAD)
  data.all.unique$death_age <- year(data.all.unique$death_date) - year(data.all.unique$FMANAD)
  data.all.unique <- data.all.unique %>%
    relocate(death_age, .after = death_date) %>%
    select(-FMANAD)

  data.all <- data.all %>%
    filter(lopnr %in% dead_model_1$lopnr) %>%
    mutate(BDAT = ymd(BDAT)) %>%
    relocate(death_date, .after = BDAT)

  data.all <- data.all %>%
    mutate(DIAG = paste(DIAG1, DIAG2, DIAG3, DIAG4, DIAG5, DIAG6, DIAG7, DIAG8, sep = "-")) %>%
    select(lopnr, FMANAD, KOEN, BDAT, DIAG, VDG1, death_date, mosaicgrupp_sthlm) %>%
    mutate(VDG1_numeric = as.numeric(VDG1)) %>%
    mutate(VDG1_group = case_when(
      VDG1_numeric >= 1 & VDG1_numeric <= 69 ~ 1,
      VDG1_numeric %in% c(70, 89) ~ 2,
      VDG1_numeric == 91 ~ 3,
      VDG1_numeric %in% c(74, 75, 77, 78) ~ 4,
      TRUE ~ NA_real_
    ))

  data.all$FMANAD <- as.character(data.all$FMANAD)
  data.all$FMANAD <- paste0(substr(data.all$FMANAD, 1, 4), "-", substr(data.all$FMANAD, 5, 6), "-15")
  data.all$FMANAD <- as.Date(data.all$FMANAD)
  data.all$death_date <- as.Date(data.all$death_date)
  data.all$death_age <- year(data.all$death_date) - year(data.all$FMANAD)
  data.all <- data.all %>%
    relocate(death_age, .after = death_date)

  data.all <- data.all %>%
    filter(death_date - BDAT <= years(3))

  for (j in 1:nrow(rank_2000)) {
    chapter_name <- rank_2000$chapter[j]

    code_list <- substr(unlist(strsplit(rank_2000$codes[j], "_")), 1, 3)

    data.all[[chapter_name]] <- sapply(data.all$DIAG, function(diag_string) {
      diag_codes <- strsplit(diag_string, "-")[[1]]
      diag_codes <- trimws(diag_codes[diag_codes != ""])
      diag_codes_truncated <- substr(diag_codes, 1, 3)
      any(diag_codes_truncated %in% code_list)
    }) * 1
  }

  summary_by_lopnr <- data.all %>%
    group_by(lopnr) %>%
    summarise(across(all_of(rank_2000$chapter), max, .names = "{.col}"),
              across(c(death_date, KOEN, death_age, mosaicgrupp_sthlm), first)
    ) %>%
    relocate(KOEN, death_date, death_age, mosaicgrupp_sthlm, .after = lopnr)

  data.all <- data.all %>%
    filter(VDG1_group %in% c(1, 2, 3, 4),
           (death_date - BDAT) <= (365.25 * 1)) %>%
    mutate(role = case_when(
      VDG1_group == 1 ~ "dr",
      VDG1_group == 2 ~ "nurse",
      VDG1_group == 3 ~ "as_nurse",
      VDG1_group == 4 ~ "param",
      TRUE ~ NA_character_
    ))

  data.all <- data.all %>%
    mutate(
      visits_year_1_dr = 0,
      visits_1_months_dr = 0,
      visits_year_1_nurse = 0,
      visits_1_months_nurse = 0,
      visits_year_1_asnurse = 0,
      visits_1_months_asnurse = 0,
      visits_year_1_param = 0,
      visits_1_months_param = 0,
    )

  data.all <- data.all %>%
    mutate(months_before_death = as.numeric((death_date - BDAT) / 30.44)) %>%
    relocate(months_before_death, .after = BDAT)

  data.all <- data.all %>%
    mutate(time_window = case_when(
      months_before_death <= 1 ~ "1_months",
      months_before_death > 1 & months_before_death <= 12 ~ "year_1",
      TRUE ~ NA_character_
    )) %>%
    relocate(time_window, .after = months_before_death)

  data.all <- data.all %>%
    relocate(role, .after = time_window)

  data.all <- data.all %>%
    mutate(
      visits_1_months_nurse = ifelse(role == "nurse" & time_window == "1_months", 1, 0),
      visits_3_months_nurse = ifelse(role == "nurse" & time_window == "3_months", 1, 0),
      visits_year_1_nurse = ifelse(role == "nurse" & time_window == "year_1", 1, 0),

      visits_1_months_dr = ifelse(role == "dr" & time_window == "1_months", 1, 0),
      visits_3_months_dr = ifelse(role == "dr" & time_window == "3_months", 1, 0),
      visits_year_1_dr = ifelse(role == "dr" & time_window == "year_1", 1, 0),

      visits_1_months_asnurse = ifelse(role == "as_nurse" & time_window == "1_months", 1, 0),
      visits_3_months_asnurse = ifelse(role == "as_nurse" & time_window == "3_months", 1, 0),
      visits_year_1_asnurse = ifelse(role == "as_nurse" & time_window == "year_1", 1, 0),

      visits_1_months_param = ifelse(role == "param" & time_window == "1_months", 1, 0),
      visits_3_months_param = ifelse(role == "param" & time_window == "3_months", 1, 0),
      visits_year_1_param = ifelse(role == "param" & time_window == "year_1", 1, 0)
    ) %>%
    group_by(lopnr) %>%
    summarise(
      across(starts_with("visits_"), sum, na.rm = TRUE),
      across(c(death_date, KOEN, death_age, mosaicgrupp_sthlm), first)
    ) %>%
    relocate(KOEN, death_date, death_age, mosaicgrupp_sthlm, .after = lopnr)

  missing_cols3 <- setdiff(names(data.all), names(summary_by_lopnr))
  for (col in missing_cols3) {
    summary_by_lopnr[[col]] <- 0
  }

  data.complete <- bind_rows(data.all, summary_by_lopnr) %>%
    arrange(lopnr) %>%
    group_by(lopnr) %>%
    summarise(
      across(starts_with("visits_"), ~ sum(.x, na.rm = TRUE)),
      across(`A00-A09`:`Z80-Z99`, ~ sum(as.numeric(.x > 0), na.rm = TRUE)))

  merged_data_ovr3 <- data.all.unique %>%
    left_join(data.complete, by = "lopnr")

  merged_data_ovr3 <- merged_data_ovr3 %>%
    mutate(across(setdiff(names(data.complete), "lopnr"), ~replace_na(., 0)))

  save(merged_data_ovr3, file = "data/merged_data_ovr3.RData")

  rm(list = ls())
  gc()
}

# OVR 2023-2024
{
  load(file = "data/ovr2023.RData")
  load(file = "data/ovr2024.RData")
  load(file = "data/dead_model_1.RData")
  load(file = "data/icd.RData")
  rank_2000 <- icd

  ovr <- bind_rows(ovr2023, ovr2024) %>%
    arrange(lopnr)
  data.all <- ovr %>% left_join(dead_model_1 %>% select(lopnr, death_date), by = "lopnr")

  data.all.unique <- data.all %>%
    filter(lopnr %in% dead_model_1$lopnr) %>%
    distinct(lopnr, .keep_all = TRUE)

  data.all.unique <- data.all.unique %>%
    select(lopnr, KOEN, FMANAD, death_date, mosaicgrupp_sthlm)

  data.all.unique$FMANAD <- as.character(data.all.unique$FMANAD)
  data.all.unique$FMANAD <- paste0(substr(data.all.unique$FMANAD, 1, 4), "-", substr(data.all.unique$FMANAD, 5, 6), "-15")
  data.all.unique$FMANAD <- ymd(data.all.unique$FMANAD)
  data.all.unique$death_age <- year(data.all.unique$death_date) - year(data.all.unique$FMANAD)
  data.all.unique <- data.all.unique %>%
    relocate(death_age, .after = death_date) %>%
    select(-FMANAD)

  data.all <- data.all %>%
    filter(lopnr %in% dead_model_1$lopnr) %>%
    mutate(BDAT = ymd(BDAT)) %>%
    relocate(death_date, .after = BDAT)

  data.all <- data.all %>%
    mutate(DIAG = paste(DIAG1, DIAG2, DIAG3, DIAG4, DIAG5, DIAG6, DIAG7, DIAG8, sep = "-")) %>%
    select(lopnr, FMANAD, KOEN, BDAT, DIAG, VDG1, death_date, mosaicgrupp_sthlm) %>%
    mutate(VDG1_numeric = as.numeric(VDG1)) %>%
    mutate(VDG1_group = case_when(
      VDG1_numeric >= 1 & VDG1_numeric <= 69 ~ 1,
      VDG1_numeric %in% c(70, 89) ~ 2,
      VDG1_numeric == 91 ~ 3,
      VDG1_numeric %in% c(74, 75, 77, 78) ~ 4,
      TRUE ~ NA_real_
    ))

  data.all$FMANAD <- as.character(data.all$FMANAD)
  data.all$FMANAD <- paste0(substr(data.all$FMANAD, 1, 4), "-", substr(data.all$FMANAD, 5, 6), "-15")
  data.all$FMANAD <- as.Date(data.all$FMANAD)
  data.all$death_date <- as.Date(data.all$death_date)
  data.all$death_age <- year(data.all$death_date) - year(data.all$FMANAD)
  data.all <- data.all %>%
    relocate(death_age, .after = death_date)

  data.all <- data.all %>%
    filter(death_date - BDAT <= years(3))

  for (j in 1:nrow(rank_2000)) {
    chapter_name <- rank_2000$chapter[j]

    code_list <- substr(unlist(strsplit(rank_2000$codes[j], "_")), 1, 3)

    data.all[[chapter_name]] <- sapply(data.all$DIAG, function(diag_string) {
      diag_codes <- strsplit(diag_string, "-")[[1]]
      diag_codes <- trimws(diag_codes[diag_codes != ""])
      diag_codes_truncated <- substr(diag_codes, 1, 3)
      any(diag_codes_truncated %in% code_list)
    }) * 1
  }

  summary_by_lopnr <- data.all %>%
    group_by(lopnr) %>%
    summarise(across(all_of(rank_2000$chapter), max, .names = "{.col}"),
              across(c(death_date, KOEN, death_age, mosaicgrupp_sthlm), first)
    ) %>%
    relocate(KOEN, death_date, death_age, mosaicgrupp_sthlm, .after = lopnr)

  data.all <- data.all %>%
    filter(VDG1_group %in% c(1, 2, 3, 4),
           (death_date - BDAT) <= (365.25 * 1)) %>%
    mutate(role = case_when(
      VDG1_group == 1 ~ "dr",
      VDG1_group == 2 ~ "nurse",
      VDG1_group == 3 ~ "as_nurse",
      VDG1_group == 4 ~ "param",
      TRUE ~ NA_character_
    ))

  data.all <- data.all %>%
    mutate(
      visits_year_1_dr = 0,
      visits_1_months_dr = 0,
      visits_year_1_nurse = 0,
      visits_1_months_nurse = 0,
      visits_year_1_asnurse = 0,
      visits_1_months_asnurse = 0,
      visits_year_1_param = 0,
      visits_1_months_param = 0,
    )

  data.all <- data.all %>%
    mutate(months_before_death = as.numeric((death_date - BDAT) / 30.44)) %>%
    relocate(months_before_death, .after = BDAT)

  data.all <- data.all %>%
    mutate(time_window = case_when(
      months_before_death <= 1 ~ "1_months",
      months_before_death > 1 & months_before_death <= 12 ~ "year_1",
      TRUE ~ NA_character_
    )) %>%
    relocate(time_window, .after = months_before_death)

  data.all <- data.all %>%
    relocate(role, .after = time_window)

  data.all <- data.all %>%
    mutate(
      visits_1_months_nurse = ifelse(role == "nurse" & time_window == "1_months", 1, 0),
      visits_3_months_nurse = ifelse(role == "nurse" & time_window == "3_months", 1, 0),
      visits_year_1_nurse = ifelse(role == "nurse" & time_window == "year_1", 1, 0),

      visits_1_months_dr = ifelse(role == "dr" & time_window == "1_months", 1, 0),
      visits_3_months_dr = ifelse(role == "dr" & time_window == "3_months", 1, 0),
      visits_year_1_dr = ifelse(role == "dr" & time_window == "year_1", 1, 0),

      visits_1_months_asnurse = ifelse(role == "as_nurse" & time_window == "1_months", 1, 0),
      visits_3_months_asnurse = ifelse(role == "as_nurse" & time_window == "3_months", 1, 0),
      visits_year_1_asnurse = ifelse(role == "as_nurse" & time_window == "year_1", 1, 0),

      visits_1_months_param = ifelse(role == "param" & time_window == "1_months", 1, 0),
      visits_3_months_param = ifelse(role == "param" & time_window == "3_months", 1, 0),
      visits_year_1_param = ifelse(role == "param" & time_window == "year_1", 1, 0)
    ) %>%
    group_by(lopnr) %>%
    summarise(
      across(starts_with("visits_"), sum, na.rm = TRUE),
      across(c(death_date, KOEN, death_age, mosaicgrupp_sthlm), first)
    ) %>%
    relocate(KOEN, death_date, death_age, mosaicgrupp_sthlm, .after = lopnr)

  missing_cols3 <- setdiff(names(data.all), names(summary_by_lopnr))
  for (col in missing_cols3) {
    summary_by_lopnr[[col]] <- 0
  }

  data.complete <- bind_rows(data.all, summary_by_lopnr) %>%
    arrange(lopnr) %>%
    group_by(lopnr) %>%
    summarise(
      across(starts_with("visits_"), ~ sum(.x, na.rm = TRUE)),
      across(`A00-A09`:`Z80-Z99`, ~ sum(as.numeric(.x > 0), na.rm = TRUE)))

  merged_data_ovr3b <- data.all.unique %>%
    left_join(data.complete, by = "lopnr")

  merged_data_ovr3b <- merged_data_ovr3b %>%
    mutate(across(setdiff(names(data.complete), "lopnr"), ~replace_na(., 0)))

  save(merged_data_ovr3b, file = "data/merged_data_ovr3b.RData")

  rm(list = ls())
  gc()
}

# Combine OVR and KON per person: visits summed, diagnoses as present/absent
{
  load(file = "data/merged_data_ovr1.RData")
  load(file = "data/merged_data_ovr1b.RData")
  load(file = "data/merged_data_ovr2.RData")
  load(file = "data/merged_data_ovr2b.RData")
  load(file = "data/merged_data_ovr3.RData")
  load(file = "data/merged_data_ovr3b.RData")
  load(file = "data/kon_complete.RData")

  kon_ovr <- bind_rows(merged_data_ovr1, merged_data_ovr1b, merged_data_ovr2, merged_data_ovr2b, merged_data_ovr3, merged_data_ovr3b, kon_complete) %>%
    arrange(lopnr) %>%
    group_by(lopnr) %>%
    summarise(
      across(c(death_date, KOEN, death_age, mosaicgrupp_sthlm), first),
      across(starts_with("visits_"), ~ sum(.x, na.rm = TRUE)),
      across(`A00-A09`:`Z80-Z99`, ~ as.integer(any(.x > 0, na.rm = TRUE)))
    )

  save(kon_ovr, file = "data/kon_ovr.RData")
  rm(list = ls())
  gc()
}

# Add inpatient diagnoses (SLV)
{
  load(file = "data/kon_ovr.RData")
  load(file = "data/slv_sum.RData")

  missing_cols <- setdiff(names(kon_ovr), names(slv_sum))

  for (col in missing_cols) {
    slv_sum[[col]] <- NA
  }
  kon_ovr_slv <- bind_rows(kon_ovr, slv_sum) %>%
    arrange(lopnr) %>%
    group_by(lopnr) %>%
    summarise(
      across(c(death_date, KOEN, death_age, mosaicgrupp_sthlm), first),
      across(starts_with("visits_"), ~ sum(.x, na.rm = TRUE)),
      across(`A00-A09`:`Z80-Z99`, ~ as.integer(any(.x > 0, na.rm = TRUE)))
    )

  save(kon_ovr_slv, file = "data/kon_ovr_slv.RData")
}

# 6. Dispensed medications
# ATC grouping
{
  atc_group <- read_xlsx("documents/atc-grouping.xlsx")

  atc_group$ATC <- gsub("\\+", "|", atc_group$ATC)
  save(atc_group, file = "data/atc_group.RData")
}

# Sex, birth month, death date and Mosaic group per person
{
  load(file = "data/ovr2010.RData")
  load(file = "data/ovr2011.RData")
  load(file = "data/ovr2012.RData")
  load(file = "data/ovr2013.RData")
  load(file = "data/ovr2014.RData")
  load(file = "data/ovr2015.RData")
  load(file = "data/ovr2016.RData")
  load(file = "data/ovr2017.RData")
  load(file = "data/ovr2018.RData")
  load(file = "data/ovr2019.RData")
  load(file = "data/ovr2020.RData")
  load(file = "data/ovr2021.RData")
  load(file = "data/ovr2022.RData")
  load(file = "data/ovr2023.RData")
  load(file = "data/ovr2024.RData")
  load(file = "data/kon2010.RData")
  load(file = "data/kon2011.RData")
  load(file = "data/kon2012.RData")
  load(file = "data/kon2013.RData")

  ovr <- bind_rows(ovr2010, ovr2011, ovr2012, ovr2013, ovr2014, ovr2015, ovr2016, ovr2017, ovr2018, ovr2019, ovr2020, ovr2021, ovr2022, ovr2023, ovr2024, kon2010, kon2011, kon2012, kon2013) %>%
    arrange(lopnr) %>%
    select(lopnr, FMANAD, KOEN, mosaicgrupp_sthlm) %>%
    distinct(lopnr, .keep_all = TRUE)

  load(file = "data/dead_model_1.RData")

  fmanad_df <- ovr %>% left_join(dead_model_1, by = "lopnr") %>%
    select(lopnr, FMANAD, KOEN, death_date, mosaicgrupp_sthlm)

  save(fmanad_df, file = "data/fmanad_df.RData")
  rm(list = ls())
  gc()
}

# Prescribed drugs: restrict each year to the study population
{
  load(file = "data/dead_model_1.RData")
  recept2010 <- read_sas("work/data_250307/multisjuka_recept_2010.sas7bdat")
  recept2010 <- recept2010 %>%
    semi_join(dead_model_1, by = "lopnr")
  save(recept2010, file = "data/recept2010.RData")
  rm(list = ls())
  gc()

  load(file = "data/dead_model_1.RData")
  recept2011 <- read_sas("work/data_250307/multisjuka_recept_2011.sas7bdat")
  recept2011 <- recept2011 %>%
    semi_join(dead_model_1, by = "lopnr")
  save(recept2011, file = "data/recept2011.RData")
  rm(list = ls())
  gc()

  load(file = "data/dead_model_1.RData")
  recept2012 <- read_sas("work/data_250307/multisjuka_recept_2012.sas7bdat")
  recept2012 <- recept2012 %>%
    semi_join(dead_model_1, by = "lopnr")
  save(recept2012, file = "data/recept2012.RData")
  rm(list = ls())
  gc()

  load(file = "data/dead_model_1.RData")
  recept2013 <- read_sas("work/data_250307/multisjuka_recept_2013.sas7bdat")
  recept2013 <- recept2013 %>%
    semi_join(dead_model_1, by = "lopnr")
  save(recept2013, file = "data/recept2013.RData")
  rm(list = ls())
  gc()

  load(file = "data/dead_model_1.RData")
  recept2014 <- read_sas("work/data_250307/multisjuka_recept_2014.sas7bdat")
  recept2014 <- recept2014 %>%
    semi_join(dead_model_1, by = "lopnr")
  save(recept2014, file = "data/recept2014.RData")
  rm(list = ls())
  gc()

  recept2015 <- read_sas("work/data_250307/multisjuka_recept_2015.sas7bdat")
  load(file = "data/dead_model_1.RData")
  recept2015 <- recept2015 %>%
    semi_join(dead_model_1, by = "lopnr")
  save(recept2015, file = "data/recept2015.RData")
  rm(list = ls())
  gc()

  recept2016 <- read_sas("work/data_250307/multisjuka_recept_2016.sas7bdat")
  load(file = "data/dead_model_1.RData")
  recept2016 <- recept2016 %>%
    semi_join(dead_model_1, by = "lopnr")
  save(recept2016, file = "data/recept2016.RData")
  rm(list = ls())
  gc()

  recept2017 <- read_sas("work/data_250307/multisjuka_recept_2017.sas7bdat")
  load(file = "data/dead_model_1.RData")
  recept2017 <- recept2017 %>%
    semi_join(dead_model_1, by = "lopnr")
  save(recept2017, file = "data/recept2017.RData")
  rm(list = ls())
  gc()

  recept2018 <- read_sas("work/data_250307/multisjuka_recept_2018.sas7bdat")
  load(file = "data/dead_model_1.RData")
  recept2018 <- recept2018 %>%
    semi_join(dead_model_1, by = "lopnr")
  save(recept2018, file = "data/recept2018.RData")
  rm(list = ls())
  gc()

  recept2019 <- read_sas("work/data_250307/multisjuka_recept_2019.sas7bdat")
  load(file = "data/dead_model_1.RData")
  recept2019 <- recept2019 %>%
    semi_join(dead_model_1, by = "lopnr")
  save(recept2019, file = "data/recept2019.RData")
  rm(list = ls())
  gc()

  recept2020 <- read_sas("work/data_250307/multisjuka_recept_2020.sas7bdat")
  load(file = "data/dead_model_1.RData")
  recept2020 <- recept2020 %>%
    semi_join(dead_model_1, by = "lopnr")
  save(recept2020, file = "data/recept2020.RData")
  rm(list = ls())
  gc()

  recept2021 <- read_sas("work/data_250307/multisjuka_recept_2021.sas7bdat")
  load(file = "data/dead_model_1.RData")
  recept2021 <- recept2021 %>%
    semi_join(dead_model_1, by = "lopnr")
  save(recept2021, file = "data/recept2021.RData")
  rm(list = ls())
  gc()

  recept2022 <- read_sas("work/data_250307/multisjuka_recept_2022.sas7bdat")
  load(file = "data/dead_model_1.RData")
  recept2022 <- recept2022 %>%
    semi_join(dead_model_1, by = "lopnr")
  save(recept2022, file = "data/recept2022.RData")
  rm(list = ls())
  gc()

  recept2023 <- read_sas("work/data_250307/multisjuka_recept_2023.sas7bdat")
  load(file = "data/dead_model_1.RData")
  recept2023 <- recept2023 %>%
    semi_join(dead_model_1, by = "lopnr")
  save(recept2023, file = "data/recept2023.RData")
  rm(list = ls())
  gc()

  recept2024 <- read_sas("work/data_250307/multisjuka_recept_2024.sas7bdat")
  load(file = "data/dead_model_1.RData")
  recept2024 <- recept2024 %>%
    semi_join(dead_model_1, by = "lopnr")
  save(recept2024, file = "data/recept2024.RData")
  rm(list = ls())
  gc()

  load(file = "data/recept2010.RData")
  load(file = "data/recept2011.RData")
  load(file = "data/recept2012.RData")
  load(file = "data/recept2013.RData")
  load(file = "data/recept2014.RData")
  load(file = "data/recept2015.RData")
  load(file = "data/recept2016.RData")
  load(file = "data/recept2017.RData")
  load(file = "data/recept2018.RData")
  load(file = "data/recept2019.RData")
  load(file = "data/recept2020.RData")
  load(file = "data/recept2021.RData")
  load(file = "data/recept2022.RData")
  load(file = "data/recept2023.RData")
  load(file = "data/recept2024.RData")

  data.recept.all <- bind_rows(recept2010, recept2011, recept2012, recept2013, recept2014, recept2015,
                               recept2016, recept2017, recept2018, recept2019, recept2020, recept2021,
                               recept2022, recept2023, recept2024) %>% arrange(lopnr)

  save(data.recept.all, file = "data/data.recept.all.RData")
  rm(list = ls())
  gc()
}

# Dispensations within 1 year before death. Per ATC group, DDD summed for drugs dispensed
# at least twice (number of dispensations when no DDD is available)
{
  load(file = "data/dead_model_1.RData")
  load(file = "data/data.recept.all.RData")
  varureg <- read_sas("work/data_250307/multisjuka_varureg.sas7bdat")
  load(file = "data/atc_group.RData")

  data.recept.all <- data.recept.all %>%
    select(-FODELSEARMANAD) %>%
    arrange(lopnr) %>%
    mutate(expdat = ymd(expdat))

  data.recept.all <- data.recept.all %>%
    left_join(varureg, by = "VARUNR") %>%
    select(lopnr, VARUNR, expdat, ATCKOD, ATCTEXT, DDD, STYRKA_KLARTEXT, ENHET_DDD, DDDFORP, ant)

  load(file = "data/fmanad_df.RData")

  data.recept.all <- data.recept.all %>%
    left_join(fmanad_df, by = "lopnr")

  data.recept.all <- data.recept.all %>%
    relocate(FMANAD, KOEN, death_date, mosaicgrupp_sthlm, .after = lopnr)

  data.recept.all <- data.recept.all %>%
    filter(!is.na(FMANAD))

  data.recept.all <- data.recept.all %>%
    mutate(
      FMANAD = ymd(paste0(FMANAD, "15")),
      death_date = as.Date(death_date, format = "%Y%m%d"),
      death_age = as.numeric(format(death_date, "%Y")) - as.numeric(format(FMANAD, "%Y"))
    )

  data.recept.all <- data.recept.all %>%
    relocate(death_age, .after = death_date)

  data.recept.all <- data.recept.all %>%
    filter(death_date - expdat <= years(1))

  data.recept.all <- data.recept.all %>%
    mutate(atc_short = substr(ATCKOD, 1, 5)) %>%
    relocate(atc_short, .after = ATCKOD)

  summary_list <- list()

  data.recept.all <- data.recept.all %>%
    mutate(DDDADJ = round(DDDFORP*ant))

  data.recept.all$DDDADJ <- as.numeric(data.recept.all$DDDADJ)

  filtered_data <- data.recept.all %>%
    group_by(lopnr, atc_short) %>%
    filter(n() > 1) %>%
    ungroup()

  for (j in 1:nrow(atc_group)) {
    atc_name <- atc_group$atc_namn[j]
    atc_pattern <- atc_group$ATC[j]

    filtered <- filtered_data %>%
      filter(str_detect(atc_short, atc_pattern))

    summary_j <- filtered %>%
      group_by(lopnr) %>%
      summarise(
        !!atc_name := if (all(is.na(DDDADJ) | DDDADJ == 0)) {
          n()
        } else {
          sum(DDDADJ, na.rm = TRUE)
        },
        .groups = "drop"
      )

    summary_list[[j]] <- summary_j
  }

  final_summary <- reduce(summary_list, full_join, by = "lopnr") %>%
    left_join(fmanad_df, by = "lopnr") %>%
    arrange(lopnr)

  final_summary <- final_summary %>%
    relocate(KOEN, FMANAD, death_date, mosaicgrupp_sthlm, .after = "lopnr")

  final_summary <- final_summary %>%
    mutate(across(where(is.numeric) & !matches("lopnr"), ~replace_na(., 0)))

  summary_data <- final_summary %>%
    mutate(
      FMANAD = ymd(paste0(FMANAD, "15")),
      death_date = as.Date(death_date, format = "%Y%m%d"),
      death_age = as.numeric(format(death_date, "%Y")) - as.numeric(format(FMANAD, "%Y"))
    )

  summary_data <- summary_data %>%
    relocate(death_age, .after = death_date)

  fmanad_df <- fmanad_df %>%
    mutate(
      FMANAD = ymd(paste0(FMANAD, "15")),
      death_date = as.Date(death_date, format = "%Y%m%d"),
      death_age = as.numeric(format(death_date, "%Y")) - as.numeric(format(FMANAD, "%Y"))
    )

  sum_data <- summary_data %>%
    select(-KOEN, -FMANAD, -death_date, -death_age, -mosaicgrupp_sthlm)

  final <- left_join(fmanad_df, sum_data, by = "lopnr")

  final <- final %>%
    mutate(across(where(is.numeric), ~replace_na(., 0)))

  save(final, file = "data/final.RData")

  rm(list = ls())
  gc()
}

# 7. Analysis dataset: nursing home status (SSB), medications and HW/OT status; age >= 65
{
  load(file = "data/kon_ovr_slv.RData")
  load(file = "data/final.RData")

  load(file = "data/husl_all.RData")

  husl_ssb <- husl_all %>%
    filter(!is.na(DATSSB)) %>%
    mutate(SSB = "1") %>%
    select(lopnr, SSB)

  kon_ovr_slv_mosaic_ssb <- left_join(kon_ovr_slv, husl_ssb, by = "lopnr") %>%
    relocate(SSB, .after = mosaicgrupp_sthlm)

  kon_ovr_slv_mosaic_ssb <- kon_ovr_slv_mosaic_ssb %>%
    mutate(SSB = as.numeric(SSB)) %>%
    mutate(SSB = replace_na(SSB, 0))

  save(kon_ovr_slv_mosaic_ssb, file = "data/kon_ovr_slv_mosaic_ssb.RData")

  recept <- final %>%
    select(-KOEN, -death_date, -death_age, -mosaicgrupp_sthlm,)

  final_summary_all <- left_join(kon_ovr_slv_mosaic_ssb, recept, by = "lopnr")

  final_summary_all <- final_summary_all %>%
    relocate(death_date, KOEN, death_age, mosaicgrupp_sthlm, .after = lopnr)

  load(file = "data/cumRR_data_all.RData")

  tempdata_final <- cumRR_data_all %>%
    mutate(date = as.Date(date)) %>%
    rename(case_ssb1 = casecontrol_ssb1) %>%
    rename(case_ssb0 = casecontrol_ssb0) %>%
    mutate(year = substr(date, 1,4))

  tempdata <- tempdata_final %>%
    rename(death_date = date) %>%
    mutate(death_date = ymd(death_date)) %>%
    relocate(year, .after = temp)

  final_summary <- final_summary_all %>%
    left_join(tempdata, by = "death_date") %>%
    as_tibble() %>%
    relocate(case_ssb0, case_ssb1, temp, year, .after = SSB)

  data.comb.65 <- final_summary %>%
    filter(death_age >= 65) %>%
    select(-FMANAD)

  save(data.comb.65, file = "data/data.comb.65.RData")

  rm(list = ls())
  gc()
}

# Mosaic group as an ordered factor
{
  load(file = "data/data.comb.65.RData")

  data.comb.65$mosaicgrupp_sthlm <- factor(
    recode(data.comb.65$mosaicgrupp_sthlm,
           "Saknas" = "Missing",
           "3.Sämre bemedlade" = "Low SES",
           "2.Medelgott bemedlade" = "Medium SES",
           "1.Bättre bemedlade" = "High SES"),
    levels = c("Missing", "Low SES", "Medium SES", "High SES"),
    ordered = TRUE
  )

  save(data.comb.65, file = "data/data.comb.65.RData")
  rm(list = ls())
  gc()
}

# 8. Stochastic gradient boosting
# Community dwellers (SSB == 0)
{
  load(file = "data/data.comb.65.RData")

  data.comb.65 <- data.comb.65 %>%
    filter(SSB == 0) %>%
    select(-SSB) %>%
    filter(!is.na(case_ssb0))

  data.comb.65$KOEN <- as.factor(data.comb.65$KOEN)

  set.seed(1234)

  trainIndex_comb_all <- createDataPartition(data.comb.65$case_ssb0, p = .70,
                                             list = FALSE,
                                             times = 1)

  # Predictors: sex, age, Mosaic group, visits, and diagnoses/medications present in > 25 individuals
  Train1 <- data.comb.65 %>%
    select(all_of(c(
      "lopnr",
      "KOEN", "death_age", "mosaicgrupp_sthlm", "case_ssb0",
      "visits_year_1_dr", "visits_1_months_dr",
      "visits_year_1_nurse", "visits_1_months_nurse",
      "visits_year_1_asnurse", "visits_1_months_asnurse",
      "visits_year_1_param", "visits_1_months_param"
    )))

  Train2 <- data.comb.65 %>%
    select(lopnr, `A00-A09`:`Z80-Z99`) %>%
    select(where(~ sum(. > 0, na.rm = TRUE) > 25))

  Train3 <- data.comb.65 %>%
    select(lopnr, A01:V01)

  Train3 <- Train3[, colSums(Train3 > 0, na.rm = TRUE) > 25]

  Train <- Train1 %>%
    left_join(Train2, by = "lopnr") %>%
    left_join(Train3, by = "lopnr")

  Train_comb_all_ssb0 <- Train[ trainIndex_comb_all,] %>%
    select(-lopnr)

  Test_comb_all_ssb0  <- Train[-trainIndex_comb_all,] %>%
    select(-lopnr)

  save(Test_comb_all_ssb0, file = "data/Test_comb_all_ssb0.RData")
  save(Train_comb_all_ssb0, file = "data/Train_comb_all_ssb0.RData")

  rm(list = ls())
  gc()

  load(file = "data/Train_comb_all_ssb0.RData")

  gbm.fit_all_comb_ssb0 <- gbm(
    formula = case_ssb0 ~ .,
    distribution = "bernoulli",
    data = Train_comb_all_ssb0,
    n.trees = 20000,
    interaction.depth = 5,
    bag.fraction = 0.5,
    shrinkage = 0.001,
    cv.folds = 10,
    n.minobsinnode = 10,
    class.stratify.cv = TRUE
  )

  save(gbm.fit_all_comb_ssb0, file = "data/gbm.fit_all_comb_ssb0.RData")
  rm(list = ls())
  gc()
}

# Nursing home residents (SSB == 1)
{
  load(file = "data/data.comb.65.RData")

  data.comb.65 <- data.comb.65 %>%
    filter(SSB == 1) %>%
    select(-case_ssb0) %>%
    filter(!is.na(case_ssb1))

  data.comb.65$KOEN <- as.factor(data.comb.65$KOEN)

  set.seed(1234)

  trainIndex_comb_all <- createDataPartition(data.comb.65$case_ssb1, p = .70,
                                             list = FALSE,
                                             times = 1)

  # Predictors: sex, age, Mosaic group, visits, and diagnoses/medications present in > 25 individuals
  Train1 <- data.comb.65 %>%
    select(all_of(c(
      "lopnr",
      "KOEN", "death_age", "mosaicgrupp_sthlm", "case_ssb1",
      "visits_year_1_dr", "visits_1_months_dr",
      "visits_year_1_nurse", "visits_1_months_nurse",
      "visits_year_1_asnurse", "visits_1_months_asnurse",
      "visits_year_1_param", "visits_1_months_param"
    )))

  Train2 <- data.comb.65 %>%
    select(lopnr, `A00-A09`:`Z80-Z99`) %>%
    select(where(~ sum(. > 0, na.rm = TRUE) > 25))

  Train3 <- data.comb.65 %>%
    select(lopnr, A01:V01)

  Train3 <- Train3[, colSums(Train3 > 0, na.rm = TRUE) > 25]

  Train <- Train1 %>%
    left_join(Train2, by = "lopnr") %>%
    left_join(Train3, by = "lopnr")

  Train_comb_all_ssb1 <- Train[ trainIndex_comb_all,] %>%
    select(-lopnr)

  Test_comb_all_ssb1  <- Train[-trainIndex_comb_all,] %>%
    select(-lopnr)

  save(Test_comb_all_ssb1, file = "data/Test_comb_all_ssb1.RData")
  save(Train_comb_all_ssb1, file = "data/Train_comb_all_ssb1.RData")

  rm(list = ls())
  gc()

  load(file = "data/Train_comb_all_ssb1.RData")

  gbm.fit_all_comb_ssb1 <- gbm(
    formula = case_ssb1 ~ .,
    distribution = "bernoulli",
    data = Train_comb_all_ssb1,
    n.trees = 20000,
    interaction.depth = 5,
    bag.fraction = 0.5,
    shrinkage = 0.001,
    cv.folds = 10,
    n.minobsinnode = 10,
    class.stratify.cv = TRUE
  )

  save(gbm.fit_all_comb_ssb1, file = "data/gbm.fit_all_comb_ssb1.RData")
  rm(list = ls())
  gc()
}

# 9. Test-set performance, relative influence and ORME
# Community dwellers
{
  load(file = "data/gbm.fit_all_comb_ssb0.RData")
  load(file = "data/Test_comb_all_ssb0.RData")

  best.iter = gbm.perf(gbm.fit_all_comb_ssb0, method = "cv")

  Effects_all <- as_tibble(gbm::summary.gbm(gbm.fit_all_comb_ssb0,
                                            plotit = FALSE))%>%
    filter(rel.inf>0)

  Effects_all %>%
    dplyr::arrange(desc(rel.inf)) %>%
    dplyr::top_n(10) %>%
    ggplot(aes(x = forcats::fct_reorder(.f = var,
                                        .x = rel.inf),
               y = rel.inf,
               fill = rel.inf)) +
    geom_col() +
    coord_flip() +
    scale_color_brewer(palette = "Dark2") +
    ggthemes::theme_fivethirtyeight() +
    theme(axis.title = element_text()) +
    xlab('Features') +
    ylab('Relative Influence') +
    ggtitle("Top 10 predictors by relative influence")

  pred_test_all <- predict(gbm.fit_all_comb_ssb0, Test_comb_all_ssb0, type = "response")

  roc_all <- pROC::roc(Test_comb_all_ssb0$case_ssb0, pred_test_all, ci = TRUE, auc = TRUE)
  print(roc_all)

  plot(roc_all, legacy.axes = TRUE, print.auc = TRUE)

  # Threshold maximising sensitivity + specificity
  max.sens.spec <- data.frame(roc_all$thresholds, roc_all$sensitivities, roc_all$specificities)

  max.sens.spec$summa <- max.sens.spec$roc_all.sensitivities + max.sens.spec$roc_all.specificities

  max_summa <- max(max.sens.spec$summa)

  cut_off_all <- max.sens.spec %>%
    filter(summa == max_summa) %>%
    pull(roc_all.thresholds)

  pred_test_all_res <- as.factor(data.table::fifelse(pred_test_all > cut_off_all, 1, 0))

  confusionMatrix(pred_test_all_res, as.factor(Test_comb_all_ssb0$case_ssb0), positive = "1")

  # ORME: odds at the highest vs lowest value of the partial dependence grid
  Odds_ratios <- NULL
  colnames <- as.vector(Effects_all$var)

  for (i in 1:nrow(Effects_all)){
    pdp.out <- plot(gbm.fit_all_comb_ssb0, i.var=colnames[[i]], return.grid = T, type = "response")

    odds.0 <- (pdp.out[1, 2]) / (1 - pdp.out[1, 2])

    odds.1 <- (pdp.out[100, 2]) / (1 - pdp.out[100, 2])

    OR.ME <- odds.1 / odds.0

    Odds_ratios[[i]] <- round(OR.ME, 2)
  }

  Odds <- as.data.frame(data.table::transpose(Odds_ratios))
  colnames(Odds) <- "OddsRatio"

  all_comb_ssb0 <- cbind(Effects_all, Odds)

  save(all_comb_ssb0, file = "data/all_comb_ssb0.RData")
  write_xlsx(all_comb_ssb0, "download/Relative_Influence_all_comb_ssb0.xlsx")
}

# ORME for sex (men vs women) and Mosaic group (high vs low)
{
  pdp.sex <- plot(
    gbm.fit_all_comb_ssb0,
    i.var = "KOEN",
    return.grid = TRUE,
    type = "response"
  )

  odds.K <- pdp.sex[1,2] / (1 - pdp.sex[1,2])
  odds.M <- pdp.sex[2,2] / (1 - pdp.sex[2,2])

  ORME.sex <- odds.M / odds.K

  ORME.sex

  pdp.mosaic <- plot(
    gbm.fit_all_comb_ssb0,
    i.var = "mosaicgrupp_sthlm",
    return.grid = TRUE,
    type = "response"
  )

  low.row  <- which(pdp.mosaic[,1] == "Low SES")
  high.row <- which(pdp.mosaic[,1] == "High SES")

  odds.low <- pdp.mosaic[low.row,2] /
    (1 - pdp.mosaic[low.row,2])

  odds.high <- pdp.mosaic[high.row,2] /
    (1 - pdp.mosaic[high.row,2])

  ORME.mosaic <- odds.high / odds.low

  ORME.mosaic

  pdp.sex
  pdp.mosaic

  rm(list = ls())
  gc()
}

# Nursing home residents
{
  load(file = "data/gbm.fit_all_comb_ssb1.RData")
  load(file = "data/Test_comb_all_ssb1.RData")

  best.iter = gbm.perf(gbm.fit_all_comb_ssb1, method = "cv")

  Effects_all <- as_tibble(gbm::summary.gbm(gbm.fit_all_comb_ssb1,
                                            plotit = FALSE))%>%
    filter(rel.inf>0)

  Effects_all %>%
    dplyr::arrange(desc(rel.inf)) %>%
    dplyr::top_n(10) %>%
    ggplot(aes(x = forcats::fct_reorder(.f = var,
                                        .x = rel.inf),
               y = rel.inf,
               fill = rel.inf)) +
    geom_col() +
    coord_flip() +
    scale_color_brewer(palette = "Dark2") +
    ggthemes::theme_fivethirtyeight() +
    theme(axis.title = element_text()) +
    xlab('Features') +
    ylab('Relative Influence') +
    ggtitle("Top 10 predictors by relative influence")

  pred_test_all <- predict(gbm.fit_all_comb_ssb1, Test_comb_all_ssb1, type = "response")

  roc_all <- pROC::roc(Test_comb_all_ssb1$case_ssb1, pred_test_all, ci = TRUE, auc = TRUE)
  print(roc_all)

  plot(roc_all, legacy.axes = TRUE, print.auc = TRUE)

  # Threshold maximising sensitivity + specificity
  max.sens.spec <- data.frame(roc_all$thresholds, roc_all$sensitivities, roc_all$specificities)

  max.sens.spec$summa <- max.sens.spec$roc_all.sensitivities + max.sens.spec$roc_all.specificities

  max_summa <- max(max.sens.spec$summa)

  cut_off_all <- max.sens.spec %>%
    filter(summa == max_summa) %>%
    pull(roc_all.thresholds)

  pred_test_all_res <- as.factor(data.table::fifelse(pred_test_all > cut_off_all, 1, 0))

  confusionMatrix(pred_test_all_res, as.factor(Test_comb_all_ssb1$case_ssb1), positive = "1")

  # ORME: odds at the highest vs lowest value of the partial dependence grid
  Odds_ratios <- NULL
  colnames <- as.vector(Effects_all$var)

  for (i in 1:nrow(Effects_all)){
    pdp.out <- plot(gbm.fit_all_comb_ssb1, i.var=colnames[[i]], return.grid = T, type = "response")

    odds.0 <- (pdp.out[1, 2]) / (1 - pdp.out[1, 2])

    odds.1 <- (pdp.out[100, 2]) / (1 - pdp.out[100, 2])

    OR.ME <- odds.1 / odds.0

    Odds_ratios[[i]] <- round(OR.ME, 2)
  }

  Odds <- as.data.frame(data.table::transpose(Odds_ratios))
  colnames(Odds) <- "OddsRatio"

  all_comb_ssb1 <- cbind(Effects_all, Odds)
  save(all_comb_ssb1, file = "data/all_comb_ssb1.RData")

  write_xlsx(all_comb_ssb1 ,"download/Relative_Influence_all_comb_ssb1.xlsx")
}

# ORME for sex (men vs women) and Mosaic group (high vs low)
{
  pdp.sex <- plot(
    gbm.fit_all_comb_ssb1,
    i.var = "KOEN",
    return.grid = TRUE,
    type = "response"
  )

  odds.K <- pdp.sex[1,2] / (1 - pdp.sex[1,2])
  odds.M <- pdp.sex[2,2] / (1 - pdp.sex[2,2])

  ORME.sex <- odds.M / odds.K

  ORME.sex

  pdp.mosaic <- plot(
    gbm.fit_all_comb_ssb1,
    i.var = "mosaicgrupp_sthlm",
    return.grid = TRUE,
    type = "response"
  )

  low.row  <- which(pdp.mosaic[,1] == "Low SES")
  high.row <- which(pdp.mosaic[,1] == "High SES")

  odds.low <- pdp.mosaic[low.row,2] /
    (1 - pdp.mosaic[low.row,2])

  odds.high <- pdp.mosaic[high.row,2] /
    (1 - pdp.mosaic[high.row,2])

  ORME.mosaic <- odds.high / odds.low

  ORME.mosaic

  pdp.sex
  pdp.mosaic

  rm(list = ls())
  gc()
}

# 10. ORME at the 25th, 50th, 75th and 95th percentiles of each variable's non-zero values
#     (daily scale, /365.25), relative to 0
{
  load(file = "data/data.comb.65.RData")
  load(file = "data/all_comb_ssb0.RData")
  load(file = "data/all_comb_ssb1.RData")
  load(file = "data/gbm.fit_all_comb_ssb0.RData")
  load(file = "data/gbm.fit_all_comb_ssb1.RData")
  quantile_res <- list()

  all_var <- rbind(all_comb_ssb0, all_comb_ssb1) %>%
    select(var) %>%
    distinct(var, .keep_all = TRUE)

  all_var <- all_var %>%
    filter(!startsWith(var, "`")) %>%
    filter(!startsWith(var, "mosaic")) %>%
    filter(!startsWith(var, "KOE"))

  inf0 <- gbm::summary.gbm(gbm.fit_all_comb_ssb0, plotit = FALSE) |>
    as_tibble() |>
    filter(rel.inf > 0)

  inf1 <- gbm::summary.gbm(gbm.fit_all_comb_ssb1, plotit = FALSE) |>
    as_tibble() |>
    filter(rel.inf > 0)

  vars <- unique(all_var$var)

  props <- c(0.25, 0.50, 0.75, 0.95)

  quantile_res_list <- list()

  for (v in vars) {
    message("Processing variable: ", v)

    pdp0 <- try(plot(gbm.fit_all_comb_ssb0, i.var = v, return.grid = TRUE, type = "response"), silent = TRUE)
    pdp1 <- try(plot(gbm.fit_all_comb_ssb1, i.var = v, return.grid = TRUE, type = "response"), silent = TRUE)

    if (inherits(pdp0, "try-error") || inherits(pdp1, "try-error")) next

    colnames(pdp0)[1:2] <- c("dose_annual", "prob")
    colnames(pdp1)[1:2] <- c("dose_annual", "prob")

    pdp0$dose <- pdp0$dose_annual / 365.25
    pdp1$dose <- pdp1$dose_annual / 365.25

    odds0 <- pdp0$prob / (1 - pdp0$prob)
    odds1 <- pdp1$prob / (1 - pdp1$prob)

    baseline_odds0 <- approx(x = pdp0$dose, y = odds0, xout = 0, rule = 2)$y
    baseline_odds1 <- approx(x = pdp1$dose, y = odds1, xout = 0, rule = 2)$y

    if (is.na(baseline_odds0) || is.na(baseline_odds1) || baseline_odds0 == 0 || baseline_odds1 == 0) next

    rel_odds0 <- odds0 / baseline_odds0
    rel_odds1 <- odds1 / baseline_odds1

    combined_odds <- rbind(
      data.frame(dose = pdp0$dose, odds = rel_odds0, model = "Community dwellers"),
      data.frame(dose = pdp1$dose, odds = rel_odds1, model = "Nursing home residents")
    )

    med <- data.comb.65[[v]]
    med <- med[!is.na(med) & med != 0] / 365.25

    if (length(med) == 0) next

    dose_cut <- stats::quantile(med, probs = 0.975, na.rm = TRUE)
    combined_odds_filt <- combined_odds[combined_odds$dose <= dose_cut, , drop = FALSE]

    if (nrow(combined_odds_filt) == 0) next

    qs <- stats::quantile(med, probs = props, na.rm = TRUE)

    models <- unique(combined_odds_filt$model)
    for (m in models) {
      df_m <- combined_odds_filt[combined_odds_filt$model == m, c("dose", "odds")]
      df_m <- df_m[order(df_m$dose), ]
      df_m <- df_m[!duplicated(df_m$dose), ]

      if (nrow(df_m) < 2) next

      interp <- approx(
        x = df_m$dose,
        y = df_m$odds,
        xout = qs,
        rule = 2
      )

      rel_inf <- if (m == "Community dwellers") {
        inf0$rel.inf[match(v, inf0$var)]
      } else {
        inf1$rel.inf[match(v, inf1$var)]
      }

      out_row <- data.frame(
        variable = v,
        model    = m,
        odds25   = unname(interp$y[1]),
        odds50   = unname(interp$y[2]),
        odds75   = unname(interp$y[3]),
        odds95   = unname(interp$y[4]),
        rel_inf  = ifelse(length(rel_inf) == 0, NA_real_, rel_inf),
        stringsAsFactors = FALSE
      )

      quantile_res_list[[length(quantile_res_list) + 1]] <- out_row
    }
  }

  quantile_res_RR <- dplyr::bind_rows(quantile_res_list)
  save(quantile_res_RR, file = "data/quantile_res_RR.RData")
  write_xlsx(quantile_res_RR, "download/quantile_res_RR.xlsx")

  rm(list = ls())
  gc()
}

# 11. Variables with the same ORME direction in both strata
{
  load(file = "data/quantile_res_RR.RData")

  quantile_res <- quantile_res_RR %>%
    mutate(rel_inf = round(rel_inf, 2)) %>%
    mutate(odds25 = round(odds25, 2)) %>%
    mutate(odds50 = round(odds50, 2)) %>%
    mutate(odds75 = round(odds75, 2)) %>%
    mutate(odds95 = round(odds95, 2))

  relinf_both <- quantile_res %>%
    filter(rel_inf >= 1) %>%
    group_by(variable) %>%
    filter(n_distinct(model) == 2) %>%
    ungroup()

  vars_below_1_both_models <- relinf_both %>%
    group_by(variable) %>%
    summarise(across(starts_with("odds"), ~all(.x <= 1, na.rm = TRUE))) %>%
    filter(if_all(starts_with("odds"), ~ .x == TRUE)) %>%
    pull(variable)

  vars_below1 <- quantile_res %>%
    filter(variable %in% vars_below_1_both_models) %>%
    filter(rel_inf >= 0.75)
  vars_above_1_both_models <- relinf_both %>%
    group_by(variable) %>%
    summarise(across(starts_with("odds"), ~all(.x >= 1, na.rm = TRUE))) %>%
    filter(if_all(starts_with("odds"), ~ .x == TRUE)) %>%
    pull(variable)

  vars_above1 <- quantile_res %>%
    filter(variable %in% vars_above_1_both_models) %>%
    filter(rel_inf >= 0.5)

  vars_above50_both <- relinf_both %>%
    group_by(variable) %>%
    summarise(include = all(odds50 >= 1 & odds75 >= 1 & odds95 >= 1, na.rm = TRUE)) %>%
    filter(include) %>%
    pull(variable)

  vars_above50 <- quantile_res %>%
    filter(variable %in% vars_above50_both) %>%
    filter(rel_inf >= 0.8)
  vars_above75_both <- relinf_both %>%
    group_by(variable) %>%
    summarise(include = all(odds75 >= 1 & odds95 >= 1, na.rm = TRUE)) %>%
    filter(include) %>%
    pull(variable)

  vars_above75 <- quantile_res %>%
    filter(variable %in% vars_above75_both) %>%
    filter(rel_inf >= 0.8)

  vars_above95_both <- relinf_both %>%
    group_by(variable) %>%
    summarise(include = all(odds95 >= 1, na.rm = TRUE)) %>%
    filter(include) %>%
    pull(variable)

  vars_above95 <- quantile_res %>%
    filter(variable %in% vars_above95_both) %>%
    filter(rel_inf >= 0.9)

  vars_below_1_both_models_all <- quantile_res %>%
    group_by(variable) %>%
    summarise(across(starts_with("odds"), ~all(.x < 1, na.rm = TRUE))) %>%
    filter(if_all(starts_with("odds"), ~ .x == TRUE)) %>%
    pull(variable)

  vars_below1_all <- quantile_res %>%
    filter(variable %in% vars_below_1_both_models_all) %>%
    filter(rel_inf > 0.1)

  vars_above_1_both_models_all <- quantile_res %>%
    group_by(variable) %>%
    summarise(across(starts_with("odds"), ~all(.x > 1, na.rm = TRUE))) %>%
    filter(if_all(starts_with("odds"), ~ .x == TRUE)) %>%
    pull(variable)

  vars_above1_all <- quantile_res %>%
    filter(variable %in% vars_above_1_both_models_all) %>%
    filter(rel_inf >= 0.1)

  vars_above50_both_all <- quantile_res %>%
    group_by(variable) %>%
    summarise(include = all(odds50 > 1 & odds75 > 1 & odds95 > 1, na.rm = TRUE)) %>%
    filter(include) %>%
    pull(variable)

  vars_above50_all <- quantile_res %>%
    filter(variable %in% vars_above50_both_all) %>%
    filter(rel_inf >= 0.1)

  vars_above75_both_all <- quantile_res %>%
    group_by(variable) %>%
    summarise(include = all(odds75 >= 1 & odds95 >= 1, na.rm = TRUE)) %>%
    filter(include) %>%
    pull(variable)

  vars_above75_all <- quantile_res %>%
    filter(variable %in% vars_above75_both_all) %>%
    filter(rel_inf >= 0.1)

  vars_above95_both_all <- quantile_res %>%
    group_by(variable) %>%
    summarise(include = all(odds95 >= 1, na.rm = TRUE)) %>%
    filter(include) %>%
    pull(variable)

  vars_above95_all <- quantile_res %>%
    filter(variable %in% vars_above95_both_all) %>%
    filter(rel_inf >= 0.1)

  vars_below50_both_all <- quantile_res %>%
    group_by(variable) %>%
    summarise(include = all(odds50 <1 & odds75 < 1 & odds95 < 1, na.rm = TRUE)) %>%
    filter(include) %>%
    pull(variable)

  vars_below50_all <- quantile_res %>%
    filter(variable %in% vars_below50_both_all) %>%
    filter(rel_inf >= 0.1)

  vars_below_lowdosage_both_all <- quantile_res %>%
    group_by(variable) %>%
    summarise(include = all(odds25 <1 & odds50 <1, na.rm = TRUE)) %>%
    filter(include) %>%
    pull(variable)

  vars_below_lowdosage_all <- quantile_res %>%
    filter(variable %in% vars_below_lowdosage_both_all) %>%
    filter(rel_inf >= 0.1)

  vars_below75_both <- quantile_res %>%
    group_by(variable) %>%
    summarise(include = all(odds75 <= 1 & odds95 <= 1, na.rm = TRUE)) %>%
    filter(include) %>%
    pull(variable)

  vars_below75 <- quantile_res %>%
    filter(variable %in% vars_below75_both) %>%
    filter(rel_inf >= 0.9)

  vars_protect_below95_both <- quantile_res %>%
    group_by(variable) %>%
    summarise(include = all(odds25 <= 1 & odds50 <= 1 & odds75 <= 1, na.rm = TRUE)) %>%
    filter(include) %>%
    pull(variable)

  vars_protect_below95_both <- quantile_res %>%
    filter(variable %in% vars_protect_below95_both) %>%
    filter(rel_inf >= 0.9)

  vars_protect_below75_both <- quantile_res %>%
    group_by(variable) %>%
    summarise(include = all(odds25 <= 1 & odds50 <= 1, na.rm = TRUE)) %>%
    filter(include) %>%
    pull(variable)

  vars_protect_below75_both <- quantile_res %>%
    filter(variable %in% vars_protect_below75_both) %>%
    filter(rel_inf >= 0.9)

  vars_mxied_both <- quantile_res %>%
    group_by(variable) %>%
    summarise(include = all(odds25 <= 1 & odds50 <= 1 & odds95 >= 1, na.rm = TRUE)) %>%
    filter(include) %>%
    pull(variable)

  vars_mixed_both <- quantile_res %>%
    filter(variable %in% vars_mxied_both) %>%
    filter(rel_inf >= 0.9)
}

# 12. Variables with ORME > 1 (or < 1) and relative influence > 0.1 in both strata
{
  load(file = "data/all_comb_ssb0.RData")
  load(file = "data/all_comb_ssb1.RData")

  consistent_positive_primary <- inner_join(
    all_comb_ssb0 %>%
      select(var, rel.inf, OddsRatio),
    all_comb_ssb1 %>%
      select(var, rel.inf, OddsRatio),
    by = "var",
    suffix = c("_ssb0", "_ssb1")
  ) %>%
    filter(
      rel.inf_ssb0 > 0.1,
      rel.inf_ssb1 > 0.1,
      OddsRatio_ssb0 > 1,
      OddsRatio_ssb1 > 1
    ) %>%
    arrange(desc((rel.inf_ssb0 + rel.inf_ssb1)/2))

  write_xlsx(consistent_positive_primary, "download/consistent_positive_primary.xlsx")

  rm(list = ls())
  gc()

  load(file = "data/all_comb_ssb0.RData")
  load(file = "data/all_comb_ssb1.RData")

  consistent_negative_primary <- inner_join(
    all_comb_ssb0 %>%
      select(var, rel.inf, OddsRatio),
    all_comb_ssb1 %>%
      select(var, rel.inf, OddsRatio),
    by = "var",
    suffix = c("_ssb0", "_ssb1")
  ) %>%
    filter(
      rel.inf_ssb0 > 0.1,
      rel.inf_ssb1 > 0.1,
      OddsRatio_ssb0 < 1,
      OddsRatio_ssb1 < 1
    ) %>%
    arrange(desc((rel.inf_ssb0 + rel.inf_ssb1)/2))

  write_xlsx(consistent_negative_primary, "download/consistent_negative_primary.xlsx")
}

# 13. Deaths per HW and OT day
{
  load(file = "data/data.comb.65.RData")

  ssb0 <- data.comb.65 %>% filter(SSB == 0)
  ssb1 <- data.comb.65 %>% filter(SSB == 1)

  ssb0 <- ssb0 %>%
    filter(!is.na(case_ssb0))

  ssb1 <- ssb1 %>%
    filter(!is.na(case_ssb1))
}

{
  death_summary_ssb0 <- ssb0 %>%
    group_by(death_date, case_ssb0) %>%
    summarise(n_deaths = n(), .groups = "drop")

  death_stats_ssb0 <- death_summary_ssb0 %>%
    group_by(case_ssb0) %>%
    summarise(
      total_deaths = sum(n_deaths),
      mean_deaths = mean(n_deaths),
      median_deaths = median(n_deaths),
      sd_deaths = sd(n_deaths),
      min_deaths = min(n_deaths),
      max_deaths = max(n_deaths),
      n_days = n()
    )
  write_xlsx(death_stats_ssb0, "download/death_stats_ssb0.xlsx")
}

{
  death_summary_ssb1 <- ssb1 %>%
    group_by(death_date, case_ssb1) %>%
    summarise(n_deaths = n(), .groups = "drop")

  death_stats_ssb1 <- death_summary_ssb1 %>%
    group_by(case_ssb1) %>%
    summarise(
      total_deaths = sum(n_deaths),
      mean_deaths = mean(n_deaths),
      median_deaths = median(n_deaths),
      sd_deaths = sd(n_deaths),
      min_deaths = min(n_deaths),
      max_deaths = max(n_deaths),
      n_days = n()
    )

  write_xlsx(death_stats_ssb1, "download/death_stats_ssb1.xlsx")
}

# 14. Mortality rate ratio, HW vs OT days: Poisson regression of daily death counts
#     (quasi-Poisson also fitted for nursing home residents)
{
  load(file = "data/data.comb.65.RData")
  deaths <- data.comb.65 %>%
    filter(SSB == 0) %>%
    filter(case_ssb0 %in% c(0,1))

  daily_counts <- deaths %>%
    count(death_date, case_ssb0, name = "death_count")

  model_pois <- glm(death_count ~ factor(case_ssb0),
                    family = poisson(link = "log"),
                    data = daily_counts)
  summary(model_pois)

  exp_coef <- exp(coef(model_pois)[2])
  percent_diff <- (exp_coef - 1) * 100
  percent_diff

  overdispersion_ratio <- sum(residuals(model_pois, type = "pearson")^2) /
    model_pois$df.residual

  overdispersion_ratio
  exp(confint(model_pois))
  ci <- exp(confint(model_pois))
  ci

  deaths <- data.comb.65 %>%
    filter(SSB == 1) %>%
    filter(case_ssb1 %in% c(0,1))

  daily_counts <- deaths %>%
    count(death_date, case_ssb1, name = "death_count")

  model_pois <- glm(death_count ~ factor(case_ssb1),
                    family = poisson(link = "log"),
                    data = daily_counts)
  summary(model_pois)

  exp_coef <- exp(coef(model_pois)[2])
  percent_diff <- (exp_coef - 1) * 100
  percent_diff

  overdispersion_ratio <- sum(residuals(model_pois, type = "pearson")^2) /
    model_pois$df.residual

  overdispersion_ratio
  exp(confint(model_pois))
  ci <- exp(confint(model_pois))
  ci

  model_quasi <- glm(death_count ~ factor(case_ssb1),
                     family = quasipoisson(link = "log"),
                     data = daily_counts)

  summary(model_quasi)
  exp(confint(model_quasi))
}

# 15. Tmax on HW and OT days
{
  load(file = "data/cumRR_data_all.RData")

  tempdata_final <- cumRR_data_all
  t_ssb0_case <- tempdata_final %>% filter(casecontrol_ssb0 == 1)

  t_ssb0_case %>%
    summarise(
      mean_temp = mean(temp, na.rm = TRUE),
      sd_temp   = sd(temp, na.rm = TRUE),
      min_temp  = min(temp, na.rm = TRUE),
      max_temp  = max(temp, na.rm = TRUE)
    )

  t_ssb0_control <- tempdata_final %>% filter(casecontrol_ssb0 == 0)

  t_ssb0_control %>%
    summarise(
      mean_temp = mean(temp, na.rm = TRUE),
      sd_temp   = sd(temp, na.rm = TRUE),
      min_temp  = min(temp, na.rm = TRUE),
      max_temp  = max(temp, na.rm = TRUE)
    )

  t_ssb1_case <- tempdata_final %>% filter(casecontrol_ssb1 == 1)

  t_ssb1_case %>%
    summarise(
      mean_temp = mean(temp, na.rm = TRUE),
      sd_temp   = sd(temp, na.rm = TRUE),
      min_temp  = min(temp, na.rm = TRUE),
      max_temp  = max(temp, na.rm = TRUE)
    )

  t_ssb1_control <- tempdata_final %>% filter(casecontrol_ssb1 == 0)

  t_ssb1_control %>%
    summarise(
      mean_temp = mean(temp, na.rm = TRUE),
      sd_temp   = sd(temp, na.rm = TRUE),
      min_temp  = min(temp, na.rm = TRUE),
      max_temp  = max(temp, na.rm = TRUE)
    )
}

# 16. Flowchart counts
# All deaths 2011-2024 with residence status
{
  flytt <- read_sas("work/data_250307/multisjuka_flytt.sas7bdat")
  diseased_data <- flytt %>%
    filter(HKOD == "D") %>%
    group_by(lopnr) %>%
    summarize(count = n()) %>%
    filter(count == 1)

  diseased <- flytt %>%
    filter(lopnr %in% diseased_data$lopnr)

  diseased <- diseased %>%
    mutate(date = as.Date(HDAT, format = "%Y%m%d")) %>%
    arrange(lopnr, date) %>%
    select(-HDAT)

  na_dates <- diseased %>%
    filter(is.na(date)) %>%
    pull(lopnr) %>%
    unique()

  diseased <- diseased %>%
    filter(!lopnr %in% na_dates)

  multiple_rows <- diseased %>%
    group_by(lopnr) %>%
    filter(n() > 1) %>%
    ungroup()

  lopnr_U_before_D <- multiple_rows %>%
    group_by(lopnr) %>%
    filter(any(HKOD == "U" & lead(HKOD) == "D")) %>%
    pull(lopnr) %>%
    unique()

  diseased <- diseased %>%
    filter(!lopnr %in% lopnr_U_before_D)

  diseased <- diseased %>%
    arrange(lopnr, date) %>%
    mutate(
      I_date = ifelse(HKOD == "I", date, NA),
      D_date = ifelse(HKOD == "D", date, NA)
    ) %>%
    mutate(
      I_date = as.Date(I_date, origin = "1970-01-01"),
      D_date = as.Date(D_date, origin = "1970-01-01")
    )

  diseased <- diseased %>%
    filter(HKOD !="U")

  diseased <- diseased %>%
    group_by(lopnr) %>%
    mutate(
      D_criteria = ifelse(all(HKOD == "D"), TRUE, FALSE)) %>% ungroup()

  diseased1 <- diseased %>%
    filter(D_criteria == TRUE)

  diseased <- diseased %>%
    anti_join(diseased1, by = "lopnr")

  diseased <- diseased %>%
    arrange(lopnr, date) %>%
    group_by(lopnr) %>%
    slice_tail(n = 2) %>%
    ungroup()

  diseased <- diseased %>%
    group_by(lopnr) %>%
    mutate(
      I_date_lag = lag(I_date),
      D_date_lag = lag(D_date)
    ) %>%
    ungroup()

  diseased <- diseased %>%
    mutate(
      I_date_combined = ifelse(!is.na(I_date), I_date, I_date_lag),
      D_date_combined = ifelse(!is.na(D_date), D_date, D_date_lag)
    ) %>%
    mutate(
      I_date_combined = as.Date(I_date_combined, origin = "1970-01-01"),
      D_date_combined = as.Date(D_date_combined, origin = "1970-01-01")
    )

  diseased <- diseased %>%
    mutate(
      meets_1_years_criteria = ifelse((D_date_combined - I_date_combined) >= 365.25, TRUE, FALSE)
    )

  diseased <- diseased %>%
    select(lopnr, HKOD, date, meets_1_years_criteria)

  diseased1 <- diseased1 %>%
    select(lopnr, HKOD, date)

  diseased1 <- diseased1 %>%
    mutate(meets_1_years_criteria = NA)

  diseased_model_1 <- rbind(diseased, diseased1)

  diseased_model_3 <- diseased_model_1 %>%
    filter(HKOD != "I") %>%
    distinct(lopnr, .keep_all = TRUE)

  diseased_model_3 <- diseased_model_3 %>%
    rename(death_date = date)

  dead_figure_1 <- diseased_model_3 %>%
    filter(between(death_date,
                   as.Date("2011-01-01"),
                   as.Date("2024-12-31")))

  save(dead_figure_1, file = "data/dead_figure_1.RData")
  rm(list = ls())
  gc()
}

# Birth month and nursing home status (HUSL)
{
  load("data/dead_figure_1.RData")

  files <- list.files(
    "work/data_250307",
    pattern = "^multisjuka_husl_\\d+\\.sas7bdat$",
    full.names = TRUE
  )

  ssb_info <- files %>%
    lapply(function(file_path) {
      message("Processing SSB: ", basename(file_path))

      read_sas(file_path) %>%
        semi_join(dead_figure_1, by = "lopnr") %>%
        mutate(DATSSB = na_if(trimws(DATSSB), "")) %>%
        filter(!is.na(DATSSB)) %>%
        mutate(SSB = 1L) %>%
        select(lopnr, DATSSB, SSB)
    }) %>%
    bind_rows() %>%
    distinct(lopnr, .keep_all = TRUE)

  husl_figure_1.1 <- files %>%
    lapply(function(file_path) {
      message("Processing: ", basename(file_path))

      read_sas(file_path) %>%
        semi_join(dead_figure_1, by = "lopnr") %>%
        left_join(
          dead_figure_1 %>%
            select(
              lopnr,
              death_date,
              meets_1_years_criteria
            ),
          by = "lopnr"
        ) %>%
        select(
          lopnr,
          FMANAD,
          death_date,
          meets_1_years_criteria
        )
    }) %>%
    bind_rows() %>%
    distinct(lopnr, .keep_all = TRUE) %>%
    left_join(ssb_info, by = "lopnr")

  save(husl_figure_1.1, file = "data/husl_figure_1.1.RData")
}

# Birth month (KON)
{
  load("data/dead_figure_1.RData")

  files <- list.files(
    "work/data_250307",
    pattern = "^multisjuka_kon_\\d+\\.sas7bdat$",
    full.names = TRUE
  )

  kon_figure_1 <- files %>%
    lapply(function(file_path) {
      message("Processing: ", basename(file_path))
      df <- read_sas(file_path)

      df %>%
        semi_join(dead_figure_1, by = "lopnr") %>%
        left_join(dead_figure_1 %>%
                    select(lopnr, death_date, meets_1_years_criteria),
                  by = "lopnr") %>%
        select(lopnr, FMANAD, death_date, meets_1_years_criteria)
    }) %>%
    bind_rows() %>%
    distinct(lopnr, .keep_all = TRUE)

  save(kon_figure_1, file = "data/kon_figure_1.RData")

  rm(files, df)
  gc()
}

# Birth month (OVR)
{
  load("data/dead_figure_1.RData")

  files <- list.files(
    "work/data_250307",
    pattern = "^multisjuka_ovr_\\d+\\.sas7bdat$",
    full.names = TRUE
  )

  ovr_figure_1 <- files %>%
    lapply(function(file_path) {
      message("Processing: ", basename(file_path))
      df <- read_sas(file_path)

      df %>%
        semi_join(dead_figure_1, by = "lopnr") %>%
        left_join(dead_figure_1 %>%
                    select(lopnr, death_date, meets_1_years_criteria),
                  by = "lopnr") %>%
        select(lopnr, FMANAD, death_date, meets_1_years_criteria)
    }) %>%
    bind_rows() %>%
    distinct(lopnr, .keep_all = TRUE)

  save(ovr_figure_1, file = "data/ovr_figure_1.RData")

  rm(list = ls())
  gc()
}

# Birth month (SLV)
{
  load("data/dead_figure_1.RData")

  slv <- read_sas("work/data_250307/multisjuka_slv.sas7bdat")

  slv_figure_1 <- slv %>%
    semi_join(dead_figure_1, by = "lopnr") %>%
    left_join(dead_figure_1 %>%
                select(lopnr, death_date, meets_1_years_criteria),
              by = "lopnr") %>%
    select(any_of(c("lopnr", "FMANAD", "death_date", "meets_1_years_criteria"))) %>%
    distinct(lopnr, .keep_all = TRUE)

  save(slv_figure_1, file = "data/slv_figure_1.RData")

  rm(slv)
  gc()
}

# Combine registers and apply exclusions: age >= 65, at least 1 year of residence
{
  load("data/husl_figure_1.1.RData")
  load("data/kon_figure_1.RData")
  load("data/ovr_figure_1.RData")
  load("data/slv_figure_1.RData")

  all_figure_1 <- bind_rows(husl_figure_1.1, kon_figure_1, ovr_figure_1, slv_figure_1) %>%
    distinct(lopnr, .keep_all = TRUE)

  all_figure_1 <- all_figure_1 %>%
    mutate(
      FMANAD = ymd(paste0(FMANAD, "15")),
      death_date = as.Date(death_date, format = "%Y%m%d"),
      death_age = as.numeric(format(death_date, "%Y")) - as.numeric(format(FMANAD, "%Y"))
    ) %>%
    arrange(lopnr)

  all_figure_2 <- all_figure_1 %>%
    filter(death_age >= 65)

  all_figure_3 <- all_figure_2 %>% filter(meets_1_years_criteria == FALSE)

  all_figure_4 <- all_figure_2 %>%
    filter(meets_1_years_criteria == TRUE | is.na(meets_1_years_criteria))

  all_figure_5 <- all_figure_4 %>%
    filter(between(death_date,
                   as.Date("2011-01-01"),
                   as.Date("2024-12-31")))

  save(all_figure_1, file = "data/all_figure_1.RData")

  rm(list = ls())
  gc()
}
