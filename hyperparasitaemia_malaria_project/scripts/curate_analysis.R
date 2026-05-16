pacman:: p_load(janitor, tidyverse, flextable, rio, here, gtsummary, epikit, gt, gtExtras, knitr, survfit) ### Upload packages

devtools::install_github("Infectious-Diseases-Data-Observatory/iddoverse") ### IDDOverse package
library(iddoverse)

###################################### Set up and prepare domains SDTM ##############################################################

#-----------------------------------------------
## 1. Import tables/datasets from SDTM_2 folder
#-----------------------------------------------

dm <- import("data/SDTM_2/DM 2025-11-20.csv")
ds <- import("data/SDTM_2/DS 2025-11-20.csv")
IN <- import("data/SDTM_2/IN 2025-11-20.csv")
mb <- import("data/SDTM_2/MB 2025-11-20.csv") 
pf <- import("data/SDTM_2/PF 2025-11-20.csv")
pt <- import("data/SDTM_2/PT 2025-11-20.csv")
rs <- import("data/SDTM_2/RS 2025-11-20.csv")
sa <- import("data/SDTM_2/SA 2025-11-20.csv")
vs <- import("data/SDTM_2/VS 2025-11-20.csv")


### 1. Microbiology domain 

prepare_domain(mb, "mb")
df_mb <- prepare_domain(mb, "mb", timing_variables = "VISITDY")

df_para <- create_malaria_parasitemia_table(
  mb,
  variables = c("PFALCIPA", "PFALCIPS", "PFALCIP", "PVIVAXA", "PVIVAXS", "PVIVAX"),
  include_method = TRUE,
  include_location = FALSE,
  timing_variables = c("MBHR", "MBDY", "MBSTDY", "VISITDY", "VISITNUM", "VISIT", "EPOCH",
                       "MBEVLINT", "MBEVINTX"),
  values_funct = first
)


### 2. Demographic domain 

dm <- dm %>% mutate(SEX = if_else(SEX == "U", NA_character_, SEX))


## Review information of certain domains. 
pf_2 <- prepare_domain(pf, "PF", timing_variables = "VISITDY")
  
df_part <- create_participant_table(dm_domain = dm,
                                    vs_domain = vs)

### 3. Treatment domain - Pharmacology features identified

treat <- prepare_domain(IN, "IN", timing_variables = "VISITDY")
treat_2 <- treat %>% 
  distinct(STUDYID, USUBJID, VISITDY, INCLAS, .keep_all = T) %>% 
  select(STUDYID, DOMAIN, USUBJID, INTRT, INDECOD, INCLAS, INDOSU, VISITDY, VISITNUM) %>% 
  mutate(VISITDY = as.numeric(VISITDY))

### 4. Data exploration of others domains according to variables of interes

pharmaco_res <- prepare_domain(pf, "PF", timing_variables = "VISITDY") %>% mutate(TIME = as.numeric(TIME))
classif <- prepare_domain(rs, "RS", timing_variables = "VISITDY") %>% mutate(TIME = as.numeric(TIME))
vital_s <- prepare_domain(vs, "VS", timing_variables = "VISITDY") %>% mutate(TIME = as.numeric(TIME))
final_status <- prepare_domain(ds %>%  select(-c(DSDECOD)), 
                                 "DS",  timing_variables = "VISITDY") %>% mutate(TIME = as.numeric(TIME))
final_status_2 <- prepare_domain(ds, "DS",  timing_variables = "VISITDY") %>% mutate(TIME = as.numeric(TIME))

final_status %>% 
  filter(STUDYID == "EFYGM") %>% 
  group_by(DISPOSITION, TIME) %>%
     summarise(n = n_distinct(USUBJID)) %>% 
     ungroup() %>% 
     pivot_wider(id_cols = DISPOSITION, names_from = TIME, values_from = n) %>% 
     adorn_totals("both") %>% 
     export(here("data", "DSDECOD.xlsx"))

final_status %>% 
  tabyl(DISPOSITION)
#-----------------------------------------------
## 2. Curated SDTM database and join domains
#-----------------------------------------------

sdtm_1 <- left_join(prepare_domain(mb, "mb", timing_variables = "MBDY"), 
          prepare_domain(dm, "dm", timing_variables = "VISITDY")) %>% 
  mutate(TIME = as.numeric(TIME))  ### Demographic and micriobiology domains

sdtm_2 <- left_join(sdtm_1, select(treat_2, STUDYID, USUBJID, INTRT, 
                                   INDECOD, INCLAS, INDOSU, VISITNUM, VISITDY),
                    by = c("TIME" = "VISITDY", "STUDYID", "USUBJID")) ### Treatment domain

sdtm_3 <- full_join(sdtm_2, pharmaco_res, by = c("TIME", "STUDYID", "USUBJID")) ### Pharmaco-domain
sdtm_4 <- full_join(sdtm_3, classif, by = c("TIME", "STUDYID", "USUBJID")) ### Classification (final cases' denotation)
sdtm_5 <- full_join(sdtm_4, vital_s, by = c("TIME", "STUDYID", "USUBJID")) ### Vital signals. 
sdtm_6 <- left_join(sdtm_5, select(final_status, TIME, STUDYID, USUBJID, DISPOSITION),
                    by = c("TIME", "STUDYID", "USUBJID")) ### Vital signals. 


### 2.1. Checking data and identify potential mistakes with the outocome's variables.
sdtm_5 %>% 
  group_by(STUDYID, WHOMAL01_NA) %>% 
  summarise(n = n_distinct(USUBJID)) %>% 
  ungroup() %>% 
  pivot_wider(id_cols = STUDYID, names_from = WHOMAL01_NA, values_from = n) %>% 
  flextable()

sdtm_5 %>% 
  group_by(STUDYID, WHOMAL01_NA, INTP_NA) %>% 
  summarise(n = n_distinct(USUBJID)) %>% 
  ungroup() %>% 
  pivot_wider(id_cols = WHOMAL01_NA, names_from = INTP_NA, values_from = n) %>% 
  flextable()

rs %>%
  group_by(STUDYID, RSORRES) %>% 
  summarise(n = n_distinct(USUBJID)) %>% 
  pivot_wider(id_cols = STUDYID, names_from = RSORRES, values_from = n) %>% 
  adorn_totals("col") %>% 
  flextable()


rs %>%
  filter(STUDYID == "GPXJK") %>% 
  group_by(RSSTRESC, VISITDY) %>% 
  summarise(n = n_distinct(USUBJID)) %>% 
  pivot_wider(id_cols = R, names_from = VISITDY, values_from = n) %>% 
  adorn_totals("row") %>% 
  flextable()

ds %>%
  #filter(STUDYID == "GPXJK") %>% 
  group_by(DSDECOD, STUDYID) %>% 
  summarise(n = n_distinct(USUBJID)) %>% 
  pivot_wider(id_cols = STUDYID, names_from = DSDECOD, values_from = n) %>% 
  adorn_totals("row") %>% 
  flextable()

#-----------------------------------------------------------
## 3. Adjust names of countries by standardised abbreviation
#-----------------------------------------------------------

sdtm_6 <- sdtm_6 %>% 
  #select(-c(TIME_SOURCE.x, TIME_SOURCE.y, TIME_SOURCE.x.x, TIME_SOURCE.y.y)) %>% 
  mutate(
    code_p = paste(STUDYID, USUBJID),
    COUNTRY_NAME = case_when(
      COUNTRY == "AGO" ~ "Angola",
      COUNTRY == "BEN" ~ "Benin",
      COUNTRY == "BFA" ~ "Burkina Faso",
      COUNTRY == "CIV" ~ "Côte d’Ivoire",
      COUNTRY == "CMR" ~ "Cameroon",
      COUNTRY == "GAB" ~ "Gabon",
      COUNTRY == "KEN" ~ "Kenya",
      COUNTRY == "MOZ" ~ "Mozambique",
      COUNTRY == "NER" ~ "Niger",
      COUNTRY == "NGA" ~ "Nigeria",
      COUNTRY == "RWA" ~ "Rwanda",
      COUNTRY == "SEN" ~ "Senegal",
      COUNTRY == "UGA" ~ "Uganda",
      COUNTRY == "ZMB" ~ "Zambia",
      TRUE ~ NA_character_
    )
  )

#-----------------------------------------------------------
## 4. Select patients treated with ACT in SDTM database
#-----------------------------------------------------------

filter_patients <- sdtm_6 %>% 
  filter(str_detect(ARM, regex(
      "ARTEMETHER LUMEFANTRINE|AMODIAQUINE ARTESUNATE|ARTHEMETER LUMEFANTRINE|ARTHEMETHER|ARTESUNATE AMODIAQUINE|DIHYDROARTEMISININ PIPERAQUINE|DIHYDROARTEMISININ-PIPERAQUINE|ARTENIMOL PIPERAQUINE|Artemether-Lumefantrine|Artesunate-Amodiaquine|Dihydroartemisinin-Piperaquine",
      ignore_case = TRUE)
  )) %>% 
  distinct(code_p, .keep_all = T)  
filter_patients_2 <- as.vector(filter_patients$code_p) ### Patient's vector filter based on treatment administrated

df <- sdtm_6 %>% 
  filter(code_p %in% c(filter_patients_2)) %>% 
  filter(!(TIME == 1 & is.na(`PFALCIPA_10^6/L`)))

### Create a new variable (column) with the name of treatment administrated regardless the time (day)
treatment_patients <- sdtm_6 %>% 
  select(STUDYID, USUBJID, ARM, COUNTRY_NAME) %>% 
  distinct(STUDYID, USUBJID, .keep_all = T) %>% 
  rename("ARM_2" = ARM,
         "COUNTRY_NAME_2" = COUNTRY_NAME)

df_2 <- left_join(df, treatment_patients, by = c("STUDYID", "USUBJID"))  ### Join treatment (new variable) 

final_status_2 <- final_status_2 %>% rename("DISPOSITION_2" = DISPOSITION) %>% select(-TIME_SOURCE)
df_2.1 <- left_join(df_2, final_status_2, by = c("STUDYID", "USUBJID", "TIME"))   ### Join treatment (new variable) 

df_2.1 <- df_2.1 %>% 
  select(-c(TIME_SOURCE.x, TIME_SOURCE.y, TIME_SOURCE.x.x, TIME_SOURCE.y.y))


df_3 <- df_2.1 %>% 
  #filter(!is.na(AGE_YEARS)) %>%  ### Filter cases without or non-age
  select(STUDYID, USUBJID, COUNTRY_NAME_2, TIME, AGE_YEARS, SEX, `PFALCIPA_10^6/L`, `PFALCIPS_10^6/L`, 
         ARM_2, INTP_NA, WHOMAL01_NA, DISPOSITION, DISPOSITION_2, WEIGHT_kg, HEIGHT_cm, TEMP_C) %>%  ### We need to define which ones would be the variables to include. 
  mutate(origin = "SDTM",
         outcome = case_when(!is.na(INTP_NA) & str_detect(DISPOSITION, "LATE") ~ INTP_NA,
                             str_detect(DISPOSITION, "LATE PARASITOLOGICAL") & str_detect(WHOMAL01_NA, "LATE PARASITOLOGICAL") ~ WHOMAL01_NA,
                             WHOMAL01_NA == "ACPR" & DISPOSITION == "ACPR" ~ "ACPR",
                             str_detect(INTP_NA, "CLINICAL FAILURE") & str_detect(DISPOSITION, "LATE TREATMENT") ~ INTP_NA,
                             str_detect(DISPOSITION,"ADEQUATE PARASITOLOGICAL AND CLINICAL RESPONSE") ~ "ACPR",
                             str_detect(DISPOSITION,"EARLY TREATMENT|ETF") ~ WHOMAL01_NA,
                             str_detect(WHOMAL01_NA, "LATE TREATMENT FAILURE") ~ INTP_NA,
                             str_detect(DISPOSITION, "ADVERSE") ~ DISPOSITION_2,
                             str_detect(DISPOSITION, "ANTI-MALARIAL") ~ DISPOSITION_2,
                             str_detect(DISPOSITION, "EXCLUDED") ~ DISPOSITION_2,
                             str_detect(DISPOSITION, "LATE CLINICAL FAILURE -|LATE CLINICAL FAILURE") & !is.na(INTP_NA) ~ INTP_NA,
                             is.na(INTP_NA) & !is.na(ARM_2) & str_detect(DISPOSITION, "LATE CLINICAL FAILURE") ~ "INDETERMINATE",
                             DISPOSITION == "COMPLETED" & WHOMAL01_NA == "ACPR" ~ WHOMAL01_NA,
                             DISPOSITION == "ETF" ~ "LTF",
                             DISPOSITION == "LCF" ~ INTP_NA,
                             DISPOSITION == "LPF" ~ INTP_NA,
                             str_detect(DISPOSITION, "LOST TO FOLLOW|LOST TO F/UP|LOSS") ~ DISPOSITION_2,
                             is.na(INTP_NA) & str_detect(WHOMAL01_NA, "LATE PARASITOLOGICAL FAILURE") ~ "LTF",
                             TRUE ~ DISPOSITION,
                             ),
         outcome = if_else(str_detect(outcome, "WITHDRAW|WITHDRE") , "WITHDRAWAL", outcome),
         outcome = if_else(str_detect(outcome, "LATE PARASITOLOGICAL FAILURE") & is.na(INTP_NA), "LTF", outcome),
         outcome = if_else(str_detect(outcome, "EARLY TREATMENT FAILURE") & is.na(INTP_NA), "ELTF", outcome),
         outcome = if_else(str_detect(outcome, "LATE TREATMENT FAILURE") & is.na(INTP_NA), "LTF", outcome),
         outcome = if_else(str_detect(outcome, "LATE TREATMENT FAILURE") & !is.na(INTP_NA), INTP_NA, outcome)
         ) 

df_3 %>% tabyl(outcome)


###################################### Set up and prepare domains MAPPER ##############################################################

#-----------------------------------------------
## 5. Import MAPPER data and wwarnset
#-----------------------------------------------

path <- here("data", "MAPPER_2", "wwarnset.csv _new.csv")  ### wwarnset
readLines(path, n = 5, encoding = "Latin1")
data_wwarnset <- read.csv(path, fileEncoding = "Latin1", stringsAsFactors = FALSE)

data_wwarnset <- data_wwarnset %>% 
  mutate(code_p = str_c(sid, pid)) 

### MAPPER data
data_selected <- import(here("data", "chosen_studies.xlsx"))
data_clinical_mapper_2 <- import(here("data", "MAPPER_2", "Clinical.csv _new.csv"), locale = locale(encoding = "Latin1"))
data_outcome_mapper_2 <- import(here("data", "MAPPER_2", "Malaria_parasitaemia.csv _new.csv"), locale = locale(encoding = "Latin1"))
data_pcr_mapper_2 <- import(here("data", "MAPPER_2", "Malaria_pcr.csv _new.csv"), locale = locale(encoding = "Latin1"))
data_subject_mapper_2 <- import(here("data", "MAPPER_2", "Subject.csv _new.csv"), locale = locale(encoding = "Latin1"))
data_treatment_mapper_2 <- import(here("data", "MAPPER_2", "treatment.csv _new.csv"), locale = locale(encoding = "Latin1"))
data_parasitaemia_mapper_2 <- import(here("data", "MAPPER_2", "Malaria_parasitaemia.csv _new.csv"), locale = locale(encoding = "Latin1"))
data_wwarnset_2 <- import(here("data", "MAPPER_2", "wwarnset.csv _new.csv"), locale = locale(encoding = "Latin1"))

data_treatment_mapper_3 <- data_treatment_mapper_2 %>% 
  rename("dayofobs" = dayoftreat,
         "hourofobs" = houroftreat)

data_subject_mapper_2 %>% 
  filter(sid == "AKABP") %>% 
  distinct(pid, .keep_all = T) %>% 
  count()

data_treatment_mapper_2 %>% 
  filter(sid == "UBTXH") %>% 
  distinct(pid, .keep_all = T) %>% 
  tabyl(treat)

data_clinical_mapper_2 %>% 
  tabyl(dayofobs)

#------------------------------------------------------
## 6. Joining MAPPER datasets and involve all variables
#------------------------------------------------------

data_cons_mapper <- full_join(data_subject_mapper_2, data_clinical_mapper_2, by = c("pid", "sid", "site", "dayofobs", "hourofobs"))
data_cons_mapper_2 <- full_join(data_cons_mapper, data_outcome_mapper_2, by = c("pid", "sid", "site", "dayofobs", "hourofobs"))
data_cons_mapper_3 <- full_join(data_cons_mapper_2, data_pcr_mapper_2, by = c("pid", "sid", "site", "dayofobs", "hourofobs"))
data_cons_mapper_4 <- full_join(data_cons_mapper_3, data_treatment_mapper_3, by = c("pid", "sid", "site", "dayofobs", "hourofobs"))

### Deleting duplicates
data_cons_mapper_5 <- data_cons_mapper_4 %>%
  distinct(pid, sid, site, dayofobs, hourofobs, .keep_all = T)

### Export final MAPPER database joined 
export(data_cons_mapper_5, here("data", "data_cons_mapper.rds")) ### Export MAPPER database consolidated according to Variables of interest

#------------------------------------------------------
## 7. Import consolidated MAPPER database
#------------------------------------------------------

df_mapper <- import(here("data", "data_cons_mapper.rds")) ### Mapper consolidated

df_mapper <- df_mapper %>% 
  select(-c(id.y, disease.y, id, disease)) %>% 
  rename("id" = id.x,
         "disease" = disease.x) %>% 
  mutate(code_p = str_c(sid, pid)) ### Cleaning and rename some variables

df_mapper_2 <- df_mapper %>% 
  filter(!dayofobs < 0) ### Remove registers which day of observation are less than 0 (zero)

#---------------------------------------------------------------------------------------------
## 8. Join MAPPER (consolidated) with wwarnset to carry variables translated according to 
## code book of wwarnset (accurate information about time, parasite load and treatment)
#------------------------------------------------------#--------------------------------------

data_mapper_2 <- left_join(df_mapper_2, select(data_wwarnset, site, outcome, 
                                                 pfmicl0, pfmicl1, pfmicl2, pfmicl3, pvmicl0,
                                                 pvmicl1, pvmicl2, pvmicl3, 
                                                 FlgDV1, FlgDV2, FlgDV3, FlgDV4, FlgDV13,
                                                 FlgDV11, FlgDV15, FlgDv17, FlgDV16, FlgDV18, 
                                                 FlgDV19, FlgDV20, FlgEF1, FlgEF2, FlgEF3, FlgEF4, 
                                                 FlgEF5, FlgPCR1, FlgPCR2, FlgPCR3, 
                                                 FlgPCR4, lastday, pardenslast, pardens, recurparday,
                                                 fspecies, StatAdjusted, StatUnadjusted, code_p),
                               by = c("code_p", "dayofobs" = "lastday")) %>% 
  rename("site" = site.x) %>% 
  mutate(dayofobs = dayofobs + 1) 

data_mapper_2 %>% 
  tabyl(dayofobs, outcome)


#----------------------------------------------------
## 9. Standardise names of treatment 
#----------------------------------------------------

data_mapper_3 <- data_mapper_2 %>% 
  mutate(
    # 1) Standardize WITHOUT dose/supervision/formulation/variant tags
    treat_2 = case_when(
      treat %in% c("AL", "AL3", "AL_6", "AL_supervised", "AL_adult",
                   "AL_4_child", "AL_4_adult", "AL_child", "AL_unsupervised") ~
        "Artemether–Lumefantrine",
      
      treat %in% c("AL_PQ", "AL3+PQh0") ~
        "Artemether–Lumefantrine + Primaquine",
      
      treat %in% c("AS") ~
        "Artesunate",
      
      treat %in% c("AQ") ~
        "Amodiaquine",
      
      # Collapse ALL AS+AQ / ASAQ / FDC / nFDC / formulation(1/2) into one
      treat %in% c("AS+AQ", "ASAQ", "AS+AQ-FDC", "AS+AQ(1)FDC", "AS+AQ(2)FDC",
                   "AS+AQ(1)-FDC", "AS+AQ(2)-FDC", "AS+AQ-nFDC",
                   "AS3+AQ-FDC", "AS3+AQ-nFDC") ~
        "Artesunate + Amodiaquine",
      
      # Collapse 3-day label into same regimen (same components)
      treat %in% c("AS+SP", "AS3+SP") ~
        "Artesunate + Sulfadoxine–Pyrimethamine",
      
      treat %in% c("AS+SMP") ~
        "Artesunate + Sulfamethoxypyrazine–Pyrimethamine",
      
      # Collapse adult/child into same regimen
      treat %in% c("ASMQ", "ASMQ_adult", "ASMQ_child") ~
        "Artesunate + Mefloquine",
      
      treat %in% c("AQ+SP") ~
        "Amodiaquine + Sulfadoxine–Pyrimethamine",
      
      # Collapse dose label
      treat %in% c("CQ", "CQ_double_dose", "CQ - double dose") ~
        "Chloroquine",
      
      treat %in% c("ART+NQ") ~
        "Artemisinin + Naphthoquine",
      
      # Collapse ALL DHA-PPQ / DHAPPQ / DP / formulation a/b / variant into one
      treat %in% c("DHA-PPQ", "DHA+PQP", "DHA+PQPa", "DHA+PQPb",
                   "DHAPPQ", "DHAPQP", "DHAPQPa", "DHAPPQb", "DP") ~
        "Dihydroartemisinin–Piperaquine",
      
      treat %in% c("SP") ~
        "Sulfadoxine–Pyrimethamine",
      
      treat %in% c("NA") ~ NA_character_,
      TRUE ~ treat
    ),
    components = case_when(
      treat_2 == "Artemether–Lumefantrine" ~ list(c("artemether","lumefantrine")),
      treat_2 == "Artemether–Lumefantrine + Primaquine" ~ list(c("artemether","lumefantrine","primaquine")),
      treat_2 == "Artesunate" ~ list(c("artesunate")),
      treat_2 == "Amodiaquine" ~ list(c("amodiaquine")),
      treat_2 == "Artesunate + Amodiaquine" ~ list(c("artesunate","amodiaquine")),
      treat_2 == "Artesunate + Sulfadoxine–Pyrimethamine" ~ list(c("artesunate","sulfadoxine","pyrimethamine")),
      treat_2 == "Artesunate + Sulfamethoxypyrazine–Pyrimethamine" ~ list(c("artesunate","sulfamethoxypyrazine","pyrimethamine")),
      treat_2 == "Artesunate + Mefloquine" ~ list(c("artesunate","mefloquine")),
      treat_2 == "Artemisinin + Naphthoquine" ~ list(c("artemisinin","naphthoquine")),
      treat_2 == "Dihydroartemisinin–Piperaquine" ~ list(c("dihydroartemisinin","piperaquine")),
      treat_2 == "Sulfadoxine–Pyrimethamine" ~ list(c("sulfadoxine","pyrimethamine")),
      treat_2 == "Chloroquine" ~ list(c("chloroquine")),
      treat_2 == "Amodiaquine + Sulfadoxine–Pyrimethamine" ~ list(c("amodiaquine","sulfadoxine","pyrimethamine")),
      TRUE ~ list(NA_character_)
    ),
    
    # 4) ACT / Non-ACT / Unknown computed from components
    act_category = {
      artemisinin_core <- c("artemether","artesunate","dihydroartemisinin","artemisinin")
      act_partners     <- c("lumefantrine","amodiaquine","piperaquine","mefloquine","naphthoquine")
      
      has_core    <- map_lgl(components, ~ any(.x %in% artemisinin_core, na.rm = TRUE))
      has_partner <- map_lgl(components, ~ any(.x %in% act_partners,     na.rm = TRUE))
      case_when(
        map_lgl(components, ~ all(is.na(.x))) ~ "Unknown",
        has_core & has_partner               ~ "ACT",
        TRUE                                 ~ "Non-ACT"
      )
      }
  )

### Check the names of the countries are standardize or they need to be adjusted
countries_names <- data_mapper_3 %>% 
  select(sid, site) %>% 
  distinct(sid, site) %>% 
  flextable()  ### Information sent to Hanna Jauncey to check it out!
  
filter_mappa_patients <- data_mapper_3 %>% 
  filter(act_category == "ACT") %>% 
  distinct() 

filter_mappa_patients <- as.vector(filter_mappa_patients$code)

data_mapper_4 <-data_mapper_3 %>% 
  filter(code_p %in% c(filter_mappa_patients))

data_mapper_4 %>% 
  count(n_distinct(code_p))

data_mapper_outcome <-  data_mapper_4 %>% 
  #group_by(code_p, pcr, outcome) %>% 
  #summarise(n = n()) %>% 
  filter(code_p == "JFPER16816")
  
  
data_mapper_5 <- data_mapper_4 %>% 
  select(sid, pid, site, dayofobs, ageyears, gender, pfmicl, gfmicl, 
         treat_2, outcome, temp, weight, height) %>% 
  rename("STUDYID" = sid,
         "USUBJID" = pid,
         "COUNTRY" = site,
         "TIME" = dayofobs,
         "AGE_YEARS" = ageyears,
         "SEX" = gender,
         `PFALCIPA_10^6/L` = pfmicl,
         `PFALCIPS_10^6/L` = gfmicl,
         "ARM" = treat_2,
         "WEIGHT_kg" = weight,
         "HEIGHT_cm" = height,
         "TEMP_C" = temp) %>% 
  mutate(TIME = if_else(TIME == 0, TIME + 1, TIME))


treatment_patients_2 <- data_mapper_5 %>%
  select(STUDYID, USUBJID, ARM, COUNTRY, AGE_YEARS) %>% 
  distinct(STUDYID, USUBJID, .keep_all = T) %>% 
  rename("ARM_2" = ARM,
         "COUNTRY_NAME_2" = COUNTRY,
         "AGE_2" = AGE_YEARS)

data_mapper_5 <- left_join(data_mapper_5, treatment_patients_2, by = c("STUDYID", "USUBJID")) %>% 
  select(-ARM, -COUNTRY, -AGE_YEARS) %>% 
  rename("AGE_YEARS" = AGE_2)

data_mapper_6 <- data_mapper_5 %>% 
  filter(AGE_YEARS < 5) %>% 
  mutate(origin = "MAPPER")

data_mapper_6 %>% 
  mutate(code_p = paste(STUDYID, USUBJID)) %>% 
  count(n_distinct(code_p)) ### Check it out then, and confirm the value is correct

df_mapper %>% 
  filter(sid == "AKABP") %>% 
  distinct(pid, .keep_all = T) %>% 
  tabyl(treat)

data_mapper_6 %>% 
  tabyl(outcome)

############################################### Final work's dataset ####################################################################  

#----------------------------------------------------
## 10. Join MAPPER data consolidated and SDTM
#----------------------------------------------------

df_4 <- df_3 %>%
  select(-c(WHOMAL01_NA, DISPOSITION, DISPOSITION_2, INTP_NA))

df_final <- rbind(df_4, data_mapper_6) %>% 
  rename("pfmicl" = `PFALCIPA_10^6/L`,
         "gfmicl" = `PFALCIPS_10^6/L`)%>% 
  mutate(code_p = paste(STUDYID, USUBJID))

df_final <- df_final %>% 
  distinct(USUBJID, STUDYID, COUNTRY_NAME_2, TIME, , .keep_all = T)

df_final %>% 
  filter(STUDYID == "UBTXH") %>% 
  distinct(USUBJID, .keep_all = T) %>% 
  tabyl(ARM_2)

df_final %>% 
  tabyl(COUNTRY_NAME_2) ### Check up names of the countries

#----------------------------------------------------
## 11. Identify treatments and number of patients
#----------------------------------------------------

df_final %>% 
  group_by(ARM_2) %>% 
  summarise(n = n_distinct(code_p))

data_characteristics <- df_final %>%
  select(STUDYID, USUBJID, AGE_YEARS, SEX) %>% 
  distinct(STUDYID, USUBJID,  .keep_all = T) %>% 
  rename("age" = AGE_YEARS,
         "gender" = SEX) %>% 
  mutate(code_p = paste(STUDYID, USUBJID))


### To add socio-demographic features for each participant
df_final_2 <- left_join(df_final, select(data_characteristics, code_p,
                                         age, gender), by = "code_p")



df_final_2 %>% 
  filter(!is.na(outcome)) %>% 
  group_by(STUDYID, outcome) %>% 
  summarise(n = n_distinct(USUBJID)) %>% 
  ungroup() %>% 
  pivot_wider(id_cols = STUDYID, names_from = outcome, values_from = n) %>% 
  adorn_totals("col") %>% 
  flextable()

#---------------------------------------------------------------------------------------
## 12. Cleaning and standardising main characteristics related to Outcome and age (group)
#---------------------------------------------------------------------------------------

df_final_3 <- df_final_2 %>%  
  clean_names() %>% 
  mutate(outcome_2 = case_when(outcome == "LTFNoPcr" ~ "LTF",
                               outcome == "LTFRC" ~ "RECRUDESCENCE",
                               outcome == "LTFRI" ~ "REINFECTION",
                               outcome %in% c("MissedD28", "MissedD42", "MissedD63") ~ "LOST TO FOLLOW-UP",
                               outcome == "ELTF" ~ "LTF",
                               outcome %in% c("NoParaAtInc", "NoPfAtInc") ~ "Non-PF-BASELINE",
                               #outcome == "NonPfMonoInf" ~ "NON_PF_MONO_INFX",
                               outcome %in% c("LowHb", "LowHt", "BSGap", "NonPfMonoInf", "UNKNOWN") ~ "OTHERS",
                             TRUE ~ outcome))

df_final_3  %>% 
  filter(!is.na(outcome_2)) %>% 
  group_by(studyid, outcome_2) %>% 
  summarise(n = n_distinct(usubjid)) %>% 
  ungroup() %>% 
  #mutate(time = as.numeric(time)) %>% 
  pivot_wider(id_cols = studyid, names_from = outcome_2, values_from = n) %>% 
  adorn_totals("col") %>% 
  flextable() %>% 
  align(align = "center", part = "all")
  
  group_by(sid, intp_na) %>% 
  distinct(usubjid, )


df_final_4 <- df_final_3 %>% 
  mutate(age_5 = if_else(age < 5, "< 5", "≥ 5" ),
         cod_id = paste(studyid, usubjid),
         age_years = as.numeric(age_years),
         age_group = case_when(age >= 0 & age < 1 ~ "0-1",
                               age >= 1 & age < 5 ~ "1-4",
                               age >= 5 & age < 10 ~ "5-10",
                               age >= 10 & age < 16 ~ "11-15",
                               age >= 16 & age < 46 ~ "16-45",
                               age >= 46 ~ ">46",
                               TRUE ~ NA),
         age_group = factor(age_group, levels = c("0-1", "1-4", "5-10",
                                                  "11-15", "16-45", ">46")),
         code_p = paste(studyid, usubjid),
         time = as.numeric(time), 
         week = case_when(time >= 1 & time <= 8 ~ 1,
                          time >= 9 & time <= 15 ~ 2,
                          time >= 16 & time <= 22 ~ 3,
                          time >= 23 & time <= 30 ~ 4,
                          time >= 31 & time <= 37 ~ 6,
                          time >= 38 & time <= 44 ~ 7,
                          time >= 45 & time <= 51 ~ 8,
                          time >= 52 & time <= 58 ~ 9,
                          time >= 59 & time <= 65 ~ 10,
                          TRUE ~ NA))


#---------------------------------------------------------------------------------------
## 13. Curating name of treatments and categories for everyone.
#---------------------------------------------------------------------------------------

df_final_5 <- df_final_4 %>%
  mutate(
    arm_raw = arm_2,
    arm_norm = arm_raw %>%
      str_squish() %>%
      str_to_upper() %>%
      str_replace_all("[\u2010\u2011\u2012\u2013\u2014\u2212]", "-"),
    act_name = case_when(
      # 1) AL + Primaquine (must come before plain AL)
      str_detect(arm_norm, "ARTEMETHER\\s*-?\\s*LUMEFANTRINE") &
        str_detect(arm_norm, "PRIMAQUINE") ~
        "Artemether–Lumefantrine + Primaquine (AL+PQ)",
      # 2) AL (includes: ARTEMETHER LUMEFANTRINE, ARTEMETHER-LUMEFANTRINE, TAB, SUSP)
      str_detect(arm_norm, "ARTEMETHER\\s*-?\\s*LUMEFANTRINE") |
        str_detect(arm_norm, "ARTHEMETHER\\s*-?\\s*LUMEFANTRINE") ~
        "Artemether–Lumefantrine (AL)",
      # 3) AS-AQ + Chlorpheniramine (must come before plain AS-AQ)
      (
        str_detect(arm_norm, "ARTESUNATE\\s*-?\\s*AMODIAQUINE") |
          str_detect(arm_norm, "AMODIAQUINE\\s+ARTESUNATE") |
          str_detect(arm_norm, "ARTESUNATE\\s+AMODIAQUINE")
      ) & str_detect(arm_norm, "CHLORPHENIRAMINE") ~
        "Artesunate–Amodiaquine + Chlorpheniramine (AS–AQ+CP)",
      # 4) AS-AQ (covers: with space, hyphen, reversed order)
      (
        str_detect(arm_norm, "ARTESUNATE\\s*-?\\s*AMODIAQUINE") |
          str_detect(arm_norm, "ARTESUNATE\\s+AMODIAQUINE") |
          str_detect(arm_norm, "ARTESUNATE + AMODIAQUINE") |
          str_detect(arm_norm, "AMODIAQUINE\\s+ARTESUNATE")
      ) ~
        "Artesunate–Amodiaquine (AS–AQ)",
      # 5) DHA-PPQ (space or hyphen)
    (
      str_detect(arm_norm, "DIHYDROARTEMISININ\\s*-?\\s*PIPERAQUINE")|
      str_detect(arm_norm, "ARTENIMOL PIPERAQUINE")
    ) ~
        "Dihydroartemisinin–Piperaquine (DHA–PPQ)",
      # 6) AS-MQ
      str_detect(arm_norm, "ARTESUNATE\\s*-?\\s*MEFLOQUINE") ~
        "Artesunate–Mefloquine (AS–MQ)",
      # 7) ART + NQ  (this is in your table and was missing)
      str_detect(arm_norm, "ARTEMISININ\\s*\\+\\s*NAPHTHOQUINE") ~
        "Artemisinin + Naphthoquine (ART+NQ)",
      TRUE ~ NA_character_
    ),
    ,
    act_name = if_else(arm_norm == "ARTESUNATE + AMODIAQUINE", 
                       "Artesunate–Amodiaquine (AS–AQ)", act_name),
    act_name = if_else(arm_norm == "ARTESUNATE + MEFLOQUINE", 
                       "Artesunate–Mefloquine (AS–MQ)", act_name),
    act_category = case_when(
      is.na(act_name) ~ "Non-ACT / Unknown",
      act_name %in% c(
        "Artemether–Lumefantrine (AL)",
        "Artemether–Lumefantrine + Primaquine (AL+PQ)",
        "Artesunate–Amodiaquine (AS–AQ)",
        "Artesunate–Amodiaquine + Chlorpheniramine (AS–AQ+CP)",
        "Dihydroartemisinin–Piperaquine (DHA–PPQ)",
        "Artesunate–Mefloquine (AS–MQ)",
        "Artemisinin + Naphthoquine (ART+NQ)"
      ) ~ "ACT",
      TRUE ~ "Other"
    ),
    act_code = case_when(
      act_name == "Artemether–Lumefantrine (AL)" ~ "AL",
      act_name == "Artemether–Lumefantrine + Primaquine (AL+PQ)" ~ "AL+PQ",
      act_name == "Artesunate–Amodiaquine (AS–AQ)" ~ "AS-AQ",
      act_name == "Artesunate–Amodiaquine + Chlorpheniramine (AS–AQ+CP)" ~ "AS-AQ+CP",
      act_name == "Dihydroartemisinin–Piperaquine (DHA–PPQ)" ~ "DHA-PPQ",
      act_name == "Artesunate–Mefloquine (AS–MQ)" ~ "AS-MQ",
      act_name == "Artemisinin + Naphthoquine (ART+NQ)" ~ "ART+NQ",
      TRUE ~ NA_character_
    )
  )

#---------------------------------------------------------------------------------------
## 15. Adjusting load of parasites and time such as numeric variables
#---------------------------------------------------------------------------------------

df_final_6 <- df_final_5 %>% 
  filter(!is.na(age) & !(time == 1 & is.na(pfmicl))) %>% 
  filter(age_years < 5) %>% 
  filter(!str_detect(act_name, "Chlorpheniramine|Naphthoquine")) 

### Preliminary table to identify treatment distribution among children under five years old
df_final_6 %>%
  group_by(act_name, age_5) %>% 
  summarise(n = n_distinct(code_p), .groups = "drop") %>% 
  rename("Name treatment scheme" = act_name,
         "Age condition" = age_5,
         "Participants" = n) %>% 
  arrange(desc(Participants)) %>% 
  flextable() %>% 
  autofit()

data <- df_final_6 %>% 
  filter(act_name %in% c("Artemether–Lumefantrine (AL)", 
                         "Artesunate–Amodiaquine (AS–AQ)",
                         "Dihydroartemisinin–Piperaquine (DHA–PPQ)",
                         "Artesunate–Mefloquine (AS–MQ)"))

### Import database with names of countries decoded
countries <- import(here("data", "countries_names.xlsx"))

countries_2 <- countries %>% 
  mutate(code_2 = paste0(studyid, usubjid)) %>% 
  clean_names() %>% 
  rename("country_name_1" = country_name) %>% 
  distinct(studyid, usubjid, country_name_1, .keep_all =T)
  
data_2 <- data %>% 
  mutate(code_2 = paste0(studyid, usubjid)) 

data_3 <- left_join(data_2, select(countries_2, code_2, country_name_1),
                               by = "code_2") 

data_3 <- data_3 %>% 
  mutate(country_name = if_else(is.na(country_name_1), country_name_2, country_name_1)) %>% 
  filter(!(origin == "MAPPER" & studyid == "YYDSM"))


#data_3 %>% 
#  select(studyid, usubjid, country_name) %>% 
#  distinct(studyid, usubjid, country_name, .keep_all = T) %>% 
#  export(here("hyperparasitaemia_malaria_project", "data", "names_countries.xlsx"))

#---------------------------------------------------------------------------------------
## 16. Final database exported to be processed
#---------------------------------------------------------------------------------------

export(data_3, here("data", "data_final.rds"))

######################################## Descriptive analysis evaluation ##################################################################

#---------------------------------------------------------------------------------------
## 16. Import final database (work's dataset)
#---------------------------------------------------------------------------------------

data <- import(here("data", "data_final.rds"))

#---------------------------------------------------------------------------------------
## 17. Distribution the load parasite table by treatment scheme 
#---------------------------------------------------------------------------------------

data_2 <- data %>% 
  mutate(pfmicl = as.numeric(pfmicl),
         parasitaemia_categories = case_when(pfmicl == 0 ~ "0",
                                             pfmicl >= 1 & pfmicl <= 99 ~ "1-99",
                                             pfmicl >= 100 & pfmicl <= 999 ~ "100-999",
                                             pfmicl >= 1000 & pfmicl <= 9999 ~ "1,000-9,999",
                                             pfmicl >= 10000 & pfmicl <= 49999 ~ "10,000-49,999",
                                             pfmicl >= 50000 & pfmicl <= 99999 ~ "50,000-99,999",
                                             pfmicl >= 100000 & pfmicl <= 499999 ~ "100,000-499,999",
                                             pfmicl >= 500000 & pfmicl <= 999999 ~ "500,000-499,999",
                                             pfmicl >= 1000000  ~ "> 1,000,000",
                                             TRUE ~ NA),
         parasitaemia_categories = factor(parasitaemia_categories, levels = c("0",
                                                                              "1-99",
                                                                              "100-999",
                                                                            "1,000-9,999",
                                                                            "10,000-49,999",
                                                                            "50,000-99,999",
                                                                            "100,000-499,999",
                                                                            "500,000-499,999",
                                                                            "> 1,000,000")),
         pfmicl = as.numeric(pfmicl),
         outcome_2 = na_if(str_squish(outcome_2), "")) 

condition_pcr <- data_2 %>%
  filter(!is.na(outcome_2)) %>% 
  select(studyid, usubjid, outcome_2) %>% 
  distinct(studyid, usubjid, outcome_2) %>% 
  rename("condition" = outcome_2)

#data_3 <- left_join(data_2, select(condition_pcr, studyid, usubjid, condition),
 #                   by = c("studyid", "usubjid"))

data %>% 
  tabyl(outcome_2)

#---------------------------------------------------------------------------------------
## 18. Table 1. Distribution load's ranges of parasite among treatments scheme
#---------------------------------------------------------------------------------------

### Day 0 (baseline)

data_2 %>% 
  filter(time == 1 & !is.na(pfmicl)) %>% 
  group_by(parasitaemia_categories, act_name) %>%
  summarise(n = n_distinct(code_p), .groups = "drop") %>% 
  pivot_wider(id_cols = parasitaemia_categories, names_from = act_name, values_from = n) %>% 
  mutate(across(-1, replace_na, 0)) %>% 
  adorn_totals() %>% 
  adorn_percentages(denominator = "col") %>% 
  adorn_pct_formatting(digits = 2) %>% 
  adorn_ns(position = "front") %>% 
  flextable() %>% 
  autofit() %>% 
  align(part = "all", j = 2:5, align = "center") %>% 
  hline(i = 9, part = "body", border = fp_border_default(width = 1.5)) %>% 
  set_header_labels(parasitaemia_categories = paste0("Parasitaemia ", "\U00B5", "L"))

### 2. Day 28

data_2 %>%  
  filter(time == 29 & !is.na(pfmicl)) %>% 
  group_by(parasitaemia_categories, act_name) %>%
  summarise(n = n_distinct(code_p), .groups = "drop") %>% 
  pivot_wider(id_cols = parasitaemia_categories, names_from = act_name, values_from = n) %>% 
  mutate(across(-1, replace_na, 0)) %>% 
  adorn_totals() %>% 
  adorn_percentages(denominator = "col") %>% 
  adorn_pct_formatting(digits = 2) %>% 
  adorn_ns(position = "front") %>% 
  flextable() %>% 
  autofit() %>% 
  align(part = "all", j = 2:5, align = "center") %>% 
  hline(i = 8, part = "body", border = fp_border_default(width = 1.5)) %>% 
  set_header_labels(parasitaemia_categories = paste0("Parasitaemia ", "\U00B5", "L"))

### 3. Day 42

data_2 %>% 
  filter(time == 43 & !is.na(pfmicl)) %>% 
  group_by(parasitaemia_categories, act_name) %>%
  summarise(n = n_distinct(code_p), .groups = "drop") %>% 
  pivot_wider(id_cols = parasitaemia_categories, names_from = act_name, values_from = n) %>% 
  mutate(across(-1, replace_na, 0)) %>% 
  adorn_totals() %>% 
  adorn_percentages(denominator = "col") %>% 
  adorn_pct_formatting(digits = 2) %>% 
  adorn_ns(position = "front") %>% 
  flextable() %>% 
  autofit() %>% 
  align(part = "all", j = 2:4, align = "center") %>% 
  hline(i = 7, part = "body", border = fp_border_default(width = 1.5)) %>% 
  set_header_labels(parasitaemia_categories = paste0("Parasitaemia ", "\U00B5", "L"))

#---------------------------------------------------------------------------------------
## 19. Evaluation of distribution studies by outcomes
#---------------------------------------------------------------------------------------

data_2 %>% 
  filter(time == 1) %>% 
  group_by(act_name) %>% 
  summarise(median = median(pfmicl, na.rm = T))

#---------------------------------------------------------------------------------------
## 21. Outcome distribution by study ID and outcome
#---------------------------------------------------------------------------------------

### a. Export data of ID studies and outcome (recrudescence, relapse, reinfection)

outcome <- data_2 %>% 
  filter(!is.na(outcome_2)) %>% 
  distinct(studyid, usubjid, outcome_2) %>% 
  group_by(studyid, outcome_2) %>% 
  summarise(n = n_distinct(usubjid, na.rm = T), groups = "drop") %>% 
  pivot_wider(id_cols = studyid, names_from = outcome_2, values_from = n) %>% 
  adorn_totals("col") %>% 
  arrange(desc(Total))

outcome

data_2 %>% 
  group_by(studyid, intp_na) %>% 
  summarise(n = n_distinct(cod_id), .groups = "drop") %>% 
  pivot_wider(id_cols = studyid, names_from = intp_na, values_from = n) %>% 
  mutate(across(, replace_na, 0)) %>% 
  arrange(desc(RECRUDESCENCE)) %>% 
  adorn_totals("both") %>% 
  flextable() %>% 
  autofit() %>% 
  save_as_docx(path = here("outcomes", "outcomes_distribution.docx"))


### b. Figure 1. Outcome distribution by study (recrudescence, relapse, reinfection)

outcome <- c("REINFECTION" = "salmon1", "RECRUDESCENCE" = "orangered", "INDETERMINATE" = "gold",
             "OTHERS" = "lightgray", "LTF" = "burlywood2", "ACPR" = "cadetblue2")

data_2 %>%
  filter(!is.na(outcome_2)) %>% 
  group_by(studyid, outcome_2) %>%  
  summarise(n = n_distinct(usubjid), .groups = "drop") %>% 
  pivot_wider(id_cols = studyid, names_from = outcome_2, values_from = n) %>% 
  adorn_totals("col") %>% 
  arrange(desc("Total")) 

data_2 %>%
  filter(!is.na(outcome_2)) %>% 
  filter(outcome_2 %in% c("RECRUDESCENCE", "REINFECTION", "INDETERMINATE", "ACPR", "LFT", "OTHERS")) %>% 
  group_by(studyid, outcome_2) %>% 
  summarise(n = n_distinct(usubjid)) %>% 
  mutate(
    total_n = sum(n),
    prop = n / total_n,
  ) %>% 
  arrange(desc(total_n)) %>% 
  mutate(studyid = fct_reorder(studyid, total_n)) %>% 
  ggplot(aes(x = reorder(studyid, total_n), y = n, fill = outcome_2)) +
  geom_col(position = "fill", colour = "black", na.rm = TRUE) +
  
  # NEW: Adds the count 'n' inside each colored segment
  geom_text(
    aes(label = n), 
    position = position_fill(vjust = 0.5), # Centers text vertically inside the fill
    size = 3,                              # Adjust size so it fits inside the bars
    color = "black",                       # Change to "white" if your fill colors are very dark
    fontface = "bold"
  ) +
  geom_text(
    aes(label = paste0("n = ", " ", total_n), y = 1), 
    hjust = -0.2,           # Pushes text to the right of the bar
    size = 3    
  ) +
  
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, .05))) +
  scale_fill_manual(values = c(outcome)) +
  coord_flip() + 
  theme_minimal() +
  labs(y = "Proportion", x = "Study ID", fill = "OUTCOME") +
  theme(
    axis.title = element_text(size = 10, color = "black", face = "bold"),
    legend.title = element_text(size = 10, color = "black", face = "bold"),
    axis.text.x = element_text(size = 9, color = "black"),
    legend.position = "bottom",
  ) 


  ggsave(path = here("figures and tables"),
         width = 12, height = 10, dpi = 300,
         file = "distribution_sid_outcome.png")
  

#---------------------------------------------------------------------------------------
## 22. Final database to share with Prabin and James. 
#---------------------------------------------------------------------------------------
  
data_final_2 <- data_2 %>% 
  filter(!is.na(outcome_2)) %>% 
  select(studyid, usubjid, age_years, sex, time, pfmicl, origin, outcome_2, act_name)


export(data_final_2, here("data", "data_outcome_hyperp.xlsx"))

#### Table 1. 

data_2 %>% 
  distinct(usubjid, .keep_all = T) %>% 
  tabyl(act_name)

sdtm_5 %>% 
  group_by(STUDYID) %>% 
  summarise(n = n_distinct(USUBJID))


match_study_table <- data_2 %>% 
  group_by(studyid) %>% 
  distinct(usubjid, .keep_all = T) %>% 
  summarise(n = n_distinct(usubjid),
            age = round(median(age_years), digits = 1),
            sex_F = sum(sex == "F"),
            sex_M = sum(sex == "M"), 
            act_name_AL = sum(act_name == "Artemether–Lumefantrine (AL)"),
            act_name_AS_AQ = sum(act_name == "Artesunate–Amodiaquine (AS–AQ)"),
            act_name_DHA_PPQ = sum(act_name == "Dihydroartemisinin–Piperaquine (DHA–PPQ)"),
            act_name_AS_MQ = sum(act_name == "Artesunate–Mefloquine (AS–MQ)")) %>% 
  mutate(across(c(sex_F, sex_M), ~ ifelse(. == 0, "", .))) %>% 
  arrange(desc(n)) %>% 
  mutate(num = seq_len(n())) %>% 
  flextable()

match_study_table

match_study_table %>% 
  save_as_docx(path = here("outcomes", "match_studies_table.docx"))




########## Baseline characteristics based on parasitaemia and distribution among main variables ##################################
data_3 %>%
  filter(time == 1) %>%
  #filter(time%in% c(1, 28, 42)) %>% 
  select(age, gender, act_name, pfmicl, condition) %>% 
  mutate(condition = str_to_title(condition),
         condition = if_else(is.na(condition), "NPI", condition),
         parasitaemia = case_when(pfmicl == 0 ~ "0 ~ Negative",
                                  pfmicl > 0 ~ "≥1 ~ Positive",
                                  is.na(pfmicl) ~ "No data/Not measured"),
         gender = if_else(is.na(gender)|gender=="", "No data", gender)) %>% 
  select(age, gender, parasitaemia, pfmicl, condition, act_name) %>% 
  tbl_summary(by = act_name, 
              missing = "no",
              label = list(age ~ "Age",
                           gender ~ "Gender",
                           parasitaemia ~ "Parasites detected",
                           pfmicl ~ paste0("Parasitaemia \n(baseline)"," \U00B5","L"),
                           condition ~ "Condition")) %>% 
  bold_labels()
  
data_3 %>% 
  filter(time == 29) %>% 
  mutate(parasitaemia = case_when(pfmicl == 0 ~ "0 ~ Negative",
                           pfmicl > 0 ~ "≥1 ~ Positive",
                           is.na(pfmicl) ~ "No data/Not measured")) %>% 
  tabyl(intp_na, parasitaemia, act_name)

data_3 <- data_3 %>% 
  mutate(parasitaemia = case_when(pfmicl == 0 ~ "0 ~ Negative",
                                  pfmicl > 0 ~ "≥1 ~ Positive",
                                  is.na(pfmicl) ~ "No data/Not measured")) 


ggplot(data = data_2 %>% group_by(code_p) %>% 
         filter(time == 1 & !is.na(pfmicl) & !is.na(intp_na)),
       aes(x = pfmicl)) +
  geom_histogram(aes(y = after_stat(count))) +
  scale_x_continuous(transform = "log10") +
  facet_wrap(~ act_name + intp_na)


ggplot(data = data_2 %>% filter(time %in% c(1,2,3,7,14,29, 43) & pfmicl > 100000) %>% 
         mutate(pfmicl = as.numeric(pfmicl)), aes(x = factor(time))) +
  geom_boxplot(aes(y = pfmicl, fill = act_name), na.rm = T) +
  scale_y_continuous(n.breaks = 8) +
  theme(
    legend.position = "bottom"
  )

ggplot(data = data_2 %>% filter(time %in% c(1,2,3,7,14,29, 43) & pfmicl > 100000) %>% 
         mutate(pfmicl = as.numeric(pfmicl)), aes(x = factor(time))) +
  geom_jitter(aes(y = pfmicl, fill = act_name), na.rm = T) + 
  scale_y_continuous(n.breaks = 10, transform = "log10") +
  facet_wrap(~intp_na)



ggplot(data = data %>% group_by(code_p) %>% 
         filter(time == 28) %>% 
         mutate(pfmicl = if_else((time == 28 & !is.na(pfmicl)) & is.na)), 
       aes(x = as.numeric(pfmicl))) +
  geom_histogram(aes(y = after_stat(count))) +
  scale_x_continuous(transform = "log10") +
  facet_wrap(~ act_name + intp_na)

########################################################## Survival analysis ##############################################################

pacman::p_load(ggsurvfit, survival)

data_4 <- data_2 %>% 
  mutate(recrudescence = if_else(outcome_2 == "RECRUDESCENCE", 1, 0)) %>% 
  group_by(code_p, time, recrudescence, act_name) %>% 
  summarise(n = n(), .groups = "drop") %>% 
  ungroup() %>%  
  filter(!is.na(recrudescence)) %>% 
  filter(!is.na(time))



data_4 %>% 
  tabyl(time, recrudescence)


prob <- survfit2(Surv(as.numeric(time), recrudescence) ~ 1, data = data_4 %>% filter(time < 50))
summary(prob)

ggsurvfit(prob, linewidth = 1, color = "pink2") +
  add_confidence_interval() +
  add_risktable(risktable_stats = "{n.risk} ({cum.event})") +
  add_censor_mark(color = "pink3") +
  #add_pvalue(caption="Log-rank {p.value}",
  #         location = "annotation", x = 50, y = 0.68) +
  scale_ggsurvfit() + 
  coord_cartesian(ylim = c(0.95, 1)) +
  scale_x_continuous(n.breaks = 10) +
  scale_y_continuous(n.breaks = 12, 
                     #label = scales::label_percent()
                     ) +
  add_quantile(x_value = 30, linetype = "dotted", 
               color = "grey30", linewidth = 0.8) +
  #geom_vline(xintercept = 30, linetype = "dashed", colour = "grey30") +
  labs(x = "Time (days after to start treatment)") +
  theme_bw() +
  theme(legend.position = "bottom", 
        axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(color = "black", size = 14, face = "bold"),
        legend.title = element_text(face = "bold", size = 14, color = "black"),
        legend.text = element_text(color = "black", size = 12),
  ) 



prob_2 <- survfit2(Surv(as.numeric(time), recrudescence) ~ act_name, data = data_4 %>% filter(time < 50))
summary(prob_2)

ggsurvfit(prob_2, linewidth = 1) +
  add_confidence_interval() +
  add_risktable(risktable_stats = "{n.risk} ({cum.event})") +
  add_censor_mark(color = "pink3") +
  #add_pvalue(caption="Log-rank {p.value}",
  #         location = "annotation", x = 50, y = 0.68) +
  scale_ggsurvfit() + 
  coord_cartesian(ylim = c(0.95, 1)) +
  scale_x_continuous(n.breaks = 10) +
  scale_y_continuous(n.breaks = 12, 
                     #label = scales::label_percent()
  ) 
  
  
 

model_1 <- glm(recurrence ~ sex + age_years, )
  
  

mortality_risk <- survfit2(Surv(dif_treat_outcome, death_vl) ~ 1, data = data_6) 
