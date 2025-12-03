pacman:: p_load(janitor, tidyverse, flextable, rio, here, gtsummary, epikit, gt, gtExtras, knitr) ### Upload packages
###########################################################################################################################

### Import SDTM and MAPPER with the all fields merged

mapper <- import("data/data_cons_mapper.csv")
sdtm <- import("data/data_sdtm_joined.csv")

#### Selection of variables and variables transformation for merging with sdtm 

mapper_2 <- mapper %>% 
  select(id.x, sid, pid, site, dayofobs, hourofobs, ageyears, gender, height, lastdayfup, 
         weight, gfmicl, gfmicl2, pfmicl, molecular_method_used, treat, pcr, trt1, trt2, 
         trt3, trt4, trt5) %>% 
  rename("id" = id.x)


### Review and evaluation datasets sdtm 
dm <- import("data/SDTM_2/DM 2025-11-20.csv")
ds <- import("data/SDTM_2/DS 2025-11-20.csv") %>% rename("SEQ" = DSSEQ, "EPOCHDS" = EPOCH)
mb <- import("data/SDTM_2/MB 2025-11-20.csv") 
pf <- import("data/SDTM_2/PF 2025-11-20.csv")
pt <- import("data/SDTM_2/PT 2025-11-20.csv")


dm <- dm %>% #### Demographic characteristics
  select(STUDYID, USUBJID, AGE, AGEU, SEX, ARMCD, ARM, COUNTRY)
ds <- ds %>% ### Domain state the most relevant characteristics are: Visit number, day of outcome, arm of the treatmen, and general condition associated with the informed day           
  select(STUDYID, USUBJID, SEQ, VISITNUM, VISIT, DSCDSTDY, VISITDY, DSMODIFY, EPOCHDS)
mb <- mb %>% 
  select(STUDYID, USUBJID, MBSEQ, MBTEST, MBSTRESN, MBSTRESU, MBMETHOD, VISITNUM, 
         VISITDY, EPOCH, MBDY, MBCDSTDY)
pf <- pf %>% 
  select(STUDYID, USUBJID, PFCAT, PFORRES, PFSTRESC, PFSPEC, PFMETHOD, VISITDY)
pt <- pt %>% 
  select(STUDYID, USUBJID, PTTRT, PTCLAS, PTROUTE, PTPSTRGU, VISITDY) %>% 
  distinct(STUDYID, USUBJID, PTCLAS, .keep_all = T)

pt_2 <- pt %>% ### Patients chosen with ACT treatment according with the protocol
  filter(str_detect(PTTRT, "ARTEMETHER|LUMEFANTRIN|ARTESUNATE|Dihydroartem|AMODIAQUINE|
                    AMODIAQUINE|Artesunate|Coartem|Coarsucam|ARTENIMOL|DUO-|Lumefantrine")) 

tabyl(pt_2$PTTRT)


df <- left_join(mb, select(dm, STUDYID, USUBJID, COUNTRY, AGE, AGEU, SEX, ARM),
                  by = c("STUDYID", "USUBJID")) ### Sociodemographics variables

df$MBDY %>% tabyl() %>% adorn_totals()  

df %>% 
  filter(EPOCH == "BASELINE") %>% 
  group_by(STUDYID, USUBJID) %>% 
  summarise(n = sum(MBDY, na.rm = T)) %>% 
  filter(n == 0)

ds$VISITDY %>% tabyl() %>% adorn_totals()

table(df$VISITDY, df$MBDY)


df_2 <- left_join(df, select(ds, STUDYID, USUBJID, DSDECOD, DSMODIFY, VISITDY),
                by = c("STUDYID", "USUBJID", "VISITDY"))

hist(log10(df_2$MBSTRESN[df_2$EPOCH == "BASELINE" & df_2$AGE < 5])) ### Histogram distribution parasiteamia baseline

df_3 <- left_join(df_2, select(pf, STUDYID, USUBJID, PFCAT, PFORRES, 
                               PFSTRESC, PFSPEC, PFMETHOD, VISITDY),
                  by = c("STUDYID", "USUBJID", "VISITDY"))


df_3 %>%
  filter(MBTEST == "Plasmodium falciparum, Asexual" & AGE > 5) %>% 
  mutate(pf_miss = if_else(!is.na(PFORRES), 1, 0)) %>% 
  tabyl(PFSTRESC, pf_miss) %>% 
  adorn_totals()

df_3 %>%
  filter(MBTEST == "Plasmodium falciparum, Asexual" & AGE > 5) %>% 
  mutate(pf_miss = if_else(!is.na(PFORRES), 1, 0)) %>% 
  ggplot() + 
  geom_histogram(aes(x = log10(MBSTRESN)))


df_3 %>%
  filter(MBTEST == "Plasmodium falciparum, Asexual" & AGE > 5) 


df_4 <- left_join(df_3, select(pt, STUDYID, USUBJID, PTCLAS, PTROUTE, PTPSTRGU, VISITDY),
                  by = c("STUDYID", "USUBJID", "VISITDY"))




df_3 <- df_2 %>% 
  filter(str_detect(MBTEST, "Plasmodium falciparum, Asexual"))

ds_dm_mb_pf <- left_join(ds_dm_mb, select(pf, STUDYID, USUBJID,PFORRES, PFSTRESC,
                                          PFMETHOD, VISITDY, PFCDSTDY, EPOCH),
                         by = c("STUDYID", "USUBJID", "VISITDY", "EPOCH"))




df <- ds_dm_mb_pf
df %>% 
  filter(str_detect(MBTEST, "Plasmodium falciparum, Asexual")) %>% 
  tabyl(PFSTRESC)










sdtm_2 <- sdtm %>% 
  select(STUDYID, USUBJID, SEQ, AGE, SEX, VISITNUM, VISIT, VISITDY, EPOCH, DSCDSTDY, INDECOD, INTRT, MBTEST, 
         MBORRES, PFCAT, PFSTRESC, PFCDSTDY, RSSTRESC, RSDY, RSCDSTDY
         )

sdtm %>% 
  tabyl(DSTERM)
