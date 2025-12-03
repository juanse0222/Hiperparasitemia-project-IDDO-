pacman:: p_load(janitor, tidyverse, flextable, rio, here, gtsummary, epikit, gt, gtExtras, knitr) ### Upload packages

############################################ MAPPER DATABASES #####################################################################
data_selected <- import("chosen_studies.xlsx")
data_clinical_mapper <- import(here("MAPPER", "Clinical.csv _new.csv"))
data_outcome_mapper <- import(here("MAPPER", "Malaria_parasitaemia.csv _new.csv"))
data_pcr_mapper <- import(here("MAPPER", "Malaria_pcr.csv _new.csv"))
data_subject_mapper <- import(here("MAPPER", "Subject.csv _new.csv"))
data_treatment_mapper <- import(here("MAPPER", "treatment.csv _new.csv"))
data_parasitaemia_mapper <- import(here("MAPPER", "Malaria_parasitaemia.csv _new.csv"))

data_treatment_mapper_2 <- data_treatment_mapper %>% 
  rename("dayofobs" = dayoftreat,
         "hourofobs" = houroftreat)

data_subject_mapper %>% 
  summarise(n_distinct(sid))


### Joining databases and involve all variables
data_cons_mapper <- full_join(data_subject_mapper, data_clinical_mapper, by = c("pid", "sid", "site", "dayofobs", "hourofobs"))
data_cons_mapper_2 <- full_join(data_cons_mapper, data_outcome_mapper, by = c("pid", "sid", "site", "dayofobs", "hourofobs"))
data_cons_mapper_3 <- full_join(data_cons_mapper_2, data_pcr_mapper, by = c("pid", "sid", "site", "dayofobs", "hourofobs"))
data_cons_mapper_4 <- full_join(data_cons_mapper_3, data_treatment_mapper_2, by = c("pid", "sid", "site", "dayofobs", "hourofobs"))

### Review information available so far. 
data_cons_mapper_5 <- data_cons_mapper_4 %>% 
  mutate(dayofobs_1 = if_else(is.na(dayofobs.x)|dayofobs.x == 0, dayofobs.y, dayofobs.x),
         dayofobs_2 = if_else(is.na(dayofobs_1)|dayofobs_1 == 0, dayofobs.x.x, dayofobs_1),
         dayofobs_3 = if_else(is.na(dayofobs_2)|dayofobs_2 == 0, dayofobs.y.y, dayofobs_2),
         hourofobs_1 = if_else(is.na(hourofobs.x)|hourofobs.x == 0, hourofobs.y, hourofobs.x),
         hourofobs_2 = if_else(is.na(hourofobs_1)|hourofobs_1 == 0, hourofobs.x.x, hourofobs_1),
         hourofobs_3 = if_else(is.na(hourofobs_2)|hourofobs_2 == 0, hourofobs.y.y, hourofobs_2)) %>% 
  select(-c(dayofobs_1, dayofobs_2, dayofobs.x, dayofobs.x.x, dayofobs.y, dayofobs.y.y, hourofobs_1, 
            hourofobs_2, hourofobs.x, hourofobs.y, hourofobs.x.x, hourofobs.y.y)) %>% 
  rename("dayofobs" = dayofobs_3,
         "hourofobs" = hourofobs_3)

### Deleting duplicates
data_cons_mapper_5 <- data_cons_mapper_4 %>% 
  distinct(pid, sid, site, dayofobs, .keep_all = T)

export(data_cons_mapper_5, here("data", "data_cons_mapper.csv"))

### Evaluating and reviewing dataset.

data_8 <- data_cons_mapper_6 %>% 
  mutate(recrudescence = if_else(pcr == "RC", 1, 0))


data_8 %>% 
  group_by(sid) %>% 
  summarize(n = n_distinct(pid), 
            failure = sum(recrudescence, na.rm = T),
            .groups = "drop") %>% 
  arrange(desc(failure)) %>% 
  print(n = 100)
  

data_cons_mapper_5 %>% 
  summarise(n = n_distinct(pid))

### Joining and filtering based-on studies selected. 
data_cons_mapper_6 <- left_join(data_cons_mapper_5, select(data_selected, sid, sid_inv),
                                by = "sid") %>% 
  filter(!is.na(sid_inv))

data_cons_mapper_6 %>% 
  summarise(n = n_distinct(sid)) ### Number of studies 

data_cons_mapper_6 %>% 
  summarise(n = n_distinct(pid)) ### number of patients

data_cons_mapper_6 %>% 
  mutate(lastdayfup = as.numeric(lastdayfup)) %>% 
  group_by(lastdayfup) %>% 
  summarise(n = n(), .groups = "drop") %>%  ### last follow-up table
  print(n= 100)

data_cons_mapper_6 %>% 
  group_by(pid) %>% 
  summarise(sum(is.na(lastdayfup)))
  
data_cons_mapper_6 %>%
  mutate(pfmicl = as.numeric(pfmicl)) %>% 
  filter(dayofobs == 0) %>%
  group_by(pid) %>% 
  summarise(n_missing = sum(is.na(pfmicl)), .groups = "drop") ### Missing data parasites/microlitre

data_cons_mapper_6 %>%
  mutate(lastdayfup = as.numeric(lastdayfup),
         ageyears = as.numeric(ageyears)) %>% 
  filter(dayofobs > 27) %>%
  group_by(pid) %>% 
  summarise(n = n(), .groups = "drop") %>% 
  summarise(n = n_distinct(pid)) ### Patients with follows-up equal or more than 28 days. 

data_cons_mapper_6 %>%
  mutate(ageyears = as.numeric(ageyears)) %>% 
  filter(ageyears < 5) %>% 
  group_by(pid) %>% 
  summarise(n = n(), .groups = "drop") %>% 
  summarise(n = n_distinct(pid)) ### patients under five years old

########################################################### SDTM #####################################################################################################

data_sdtm_DM <- import(here("SDTM", "DM 2025-01-10.csv")) %>% select(-c(DOMAIN)) ### DEMOGRAPHIC
data_sdtm_DS <- import(here("SDTM", "DS 2025-01-10.csv")) %>% rename("SEQ" = "DSSEQ") %>% select(-c(DOMAIN)) ### DISPOSITION
data_sdtm_IN <- import(here("SDTM", "IN 2025-01-10.csv")) %>% rename("SEQ" = "INSEQ") %>% select(-c(DOMAIN)) ### INTERVENTION
data_sdtm_MB <- import(here("SDTM", "MB 2025-01-10.csv")) %>% rename("SEQ" = "MBSEQ") %>% select(-c(DOMAIN)) ### MICROBIOLOGY
data_sdtm_PF <- import(here("SDTM", "PF 2025-01-10.csv")) %>% rename("SEQ" = "PFSEQ") %>% select(-c(DOMAIN)) ### PHARMOCHOGENETIC
data_sdtm_PT <- import(here("SDTM", "PT 2025-01-10.csv")) %>% rename("SEQ" = "PTSEQ") %>% select(-c(DOMAIN))
data_sdtm_RS <- import(here("SDTM", "RS 2025-01-10.csv")) %>% rename("SEQ" = "RSSEQ") %>% select(-c(DOMAIN))
data_sdtm_SA <- import(here("SDTM", "SA 2025-01-10.csv")) %>% rename("SEQ" = "SASEQ") %>% select(-c(DOMAIN))
data_sdtm_VS <- import(here("SDTM", "VS 2025-01-10.csv")) %>% rename("SEQ" = "VSSEQ") %>% select(-c(DOMAIN))
data_sdtm_TI <- import(here("SDTM", "TI.csv")) %>% select(-c(DOMAIN)) ### IT DOESN'T INCLUDED IN THE FINAL DATABASE
data_sdtm_TS <- import(here("SDTM", "TS.csv")) %>% rename("SEQ" = "TSSEQ") %>% select(-c(DOMAIN)) ### IT DOESN'T INCLUDED IN THE FINAL DATABASE
data_sdtm_TV <- import(here("SDTM", "TV.csv")) %>% select(-c(DOMAIN)) ### IT DOESN'T INCLUDED IN THE FINAL DATABASE

data_sdtm_1 <- full_join(data_sdtm_DS, data_sdtm_DM, by = c("STUDYID", "USUBJID"))
data_sdtm_2 <- full_join(data_sdtm_1, data_sdtm_MB, by = c("STUDYID", "USUBJID", "SEQ", "VISIT", "VISITNUM", "VISITDY", "EPOCH"))
data_sdtm_3 <- full_join(data_sdtm_2, data_sdtm_VS, by = c("STUDYID", "USUBJID", "SEQ", "VISIT", "VISITNUM", "VISITDY", "EPOCH"))  
data_sdtm_4 <- full_join(data_sdtm_3, data_sdtm_PF, by = c("STUDYID", "USUBJID", "SEQ", "VISIT", "VISITNUM", "VISITDY", "EPOCH")) 
data_sdtm_5 <- full_join(data_sdtm_4, data_sdtm_PT, by = c("STUDYID", "USUBJID", "SEQ", "VISIT", "VISITNUM", "VISITDY", "EPOCH"))
data_sdtm_6 <- full_join(data_sdtm_5, data_sdtm_IN, by = c("STUDYID", "USUBJID", "SEQ", "VISIT", "VISITNUM", "VISITDY", "EPOCH", "DOTIND"))
data_sdtm_7 <- full_join(data_sdtm_6, data_sdtm_RS, by = c("STUDYID", "USUBJID", "SEQ", "VISIT", "VISITNUM", "VISITDY", "EPOCH")) 
data_sdtm_8 <- full_join(data_sdtm_7, data_sdtm_VS, by = c("STUDYID", "USUBJID", "SEQ", "VISIT", "VISITNUM", "VISITDY", "EPOCH"))

export(data_sdtm_8, here("data", "data_sdtm_joined.csv"))

### Studies that are part of SDTM (SID)
data_sdtm <- data_sdtm_8 %>% 
  select(STUDYID, COUNTRY) %>% 
  distinct(STUDYID) %>% 
  mutate(SDTM = TRUE)

### studies that aren't part of list selected
data_sdtm_selected <- left_join(data_sdtm_8, select(data_selected, sid, sid_inv, data_format),
                                                    by = c("STUDYID" = "sid"))

data_sdtm_selected %>% 
  filter(is.na(data_format)) %>%
  select(STUDYID, data_format) %>% 
  distinct(STUDYID, data_format) ### Total five studies found that aren't part the initial list (56)

########################################################### WWARNSET ####################################################################

### 1. Curated and cleaning database WWARNset
data_wwarn <- read_csv("hyperparasitimia_malaria_project/data/WWARNset.csv", locale = locale(encoding = "Latin1"))
data_wwarn_2 <- data_wwarn %>% 
  select(sid, site) %>% 
  distinct(sid) %>% 
  mutate(wwarn = T)

### Reviewing and evaluation of database for both  sources
data_selected %>% summarise(n = n_distinct(sid))

data_selected_2 <- left_join(data_selected, select(data_wwarn_2, sid, wwarn), by = "sid")

data_not_wwarn <- data_selected_2 %>% 
  filter(is.na(wwarn)) %>% 
  select(c(sid, sid_inv, data_format)) 

### Cross and filter with not wwarnset studies and identify missing 

data_not_warnn_yes_sdtm <- left_join(data_not_wwarn, select(data_sdtm, STUDYID, SDTM), by = c("sid" = "STUDYID")) %>% 
  filter(is.na(SDTM))

### Total studies without match among WWARNset and SDTM 
data_not_warnn_yes_sdtm %>% 
  summarise(n_distinct(sid))

###################################### Studies WWARNset that should be transformed in SDTM ######################################

data_wwarn %>% 
  filter(ageyears < 5) %>% 
  mutate(recrudescence = if_else(outcome == "LTFRC", 1, 0)) %>% 
  group_by(sid) %>% 
  summarise(n = n(),
            failure = sum(recrudescence, na.rm = T)) %>% 
  arrange(desc(failure)) %>% 
  gt()

