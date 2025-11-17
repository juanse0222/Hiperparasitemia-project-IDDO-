require(pacman)
pacman::p_load(tidyverse, rio, lubridate, here, epikit, flextable, janitor, scales) ### Load packages
###################################################################################################################################

### 1. Curated and cleaning database WWARNset
data <- read_csv("data/WWARNset.csv", guess_max = 100000, locale = locale(encoding = "Latin1"))
data_chosen <- import(here("outcomes", "chosen_studies.xlsx"))

data_chosen %>% summarise(n = n_distinct(sid))

### 2. Identify how many patients corresponding to studies chosen with registers in WWARNset.

data_2 <- left_join(data, select(data_chosen, sid, pubmed), by = "sid")
data_2 <- data_2 %>% filter(!is.na(pubmed)) ### All patients that contain the wwarnset are part of the studies chosen before

data_2 %>%
  summarise(n = n_distinct(sid)) ### The total studies in wwarnset was 44 compared with 56 of studies chosen initially.

colnames(data_2) ### 56 variables,

### 2.1. Age distribution - histogram
ggplot(data_2 %>% filter(!is.na(ageyears)), aes(x = ageyears)) +
  geom_histogram(bins = 20) +
  geom_density(aes(y = after_stat(count))) +
  scale_x_continuous(n.breaks = 20)

### 3. Apply selection criteria to obtain the final number of patients

#### 3.1. Evaluation main variables and identifying behaviour and distribution

### 3.1.1. Table 2xN age categorized (< 5yr and 5yrs and older x last follow-up < 28d, 28d, >28d)
data_2 <- data_2 %>%
  mutate(
    age_5 = if_else(ageyears < 5, "<5y", "5y >"),
    day_28 = case_when(
      lastday < 28 ~ "<28",
      lastday == 28 ~ "28",
      lastday > 28 ~ ">28",
      TRUE ~ NA
    ),
    day_28 = factor(day_28, levels = c("<28", "28", ">28")),
    pfmicl0 = as.numeric(pfmicl0)
  ) ### Transformation sort of variables to evaluate forward

data_2 %>%
  tabyl(age_5, day_28) ##### Finding 818 registers without age and  4,688 register that manage age <5 and last follow-up 28 days
##### WWARNset doesn't include gender (sex)

### 3.1.2. Distribution of parasitema in day 0

median(data_2$pfmicl0, na.rm = T) ### median of load parasite visit 0  = 18,420

ggplot(data_2, aes(x = pfmicl0)) +
  geom_histogram() +
  geom_vline(xintercept = median(data_2$pfmicl0, na.rm = T), linetype = "dashed", color = "black") +
  labs(title = "Load of parasite visit 0 without inclusion criteria") +
  scale_x_continuous(transform = "log10", n.breaks = 10) +
  scale_y_continuous(n.breaks = 20) +
  theme(
    axis.text.x = element_text(angle = 90)
  )

### 3.1.3. Applying criteria under 5y and follow-up until 28 days

m <- median(data_2$pfmicl0[data_2$age_5 == "<5y" & data_2$day_28 == "28"], na.rm = TRUE) ### Median with criteria employed was = 1,5326

ggplot(data_2 %>% filter(age_5 == "<5y" & day_28 == "28"), aes(x = pfmicl0)) +
  geom_histogram() +
  geom_vline(xintercept = m, linetype = "dashed", color = "black") +
  labs(title = "Load of parasite visit 0 children under 5y and follow-up day 28") +
  scale_x_continuous(transform = "log10", n.breaks = 10) +
  scale_y_continuous(n.breaks = 20) +
  theme(
    axis.text.x = element_text(angle = 90)
  )


### 3.1.4 Characteristics and distribution of treatment.
#### Transformation of variables including names and establish those who are ACT
data_2 <- data_2 %>%
  mutate(
    treat_2 = case_when(
      treat == "AL" ~ "Artemether–Lumefantrine",
      treat == "AL_adult" ~ "Artemether–Lumefantrine (adult dose)",
      treat == "AL_child" ~ "Artemether–Lumefantrine (child dose)",
      treat == "AL_PQ" ~ "Artemether–Lumefantrine + Primaquine",
      treat == "AL_unsupervised" ~ "Artemether–Lumefantrine (unsupervised)",
      treat == "AS" ~ "Artesunate",
      treat == "AS+AQ" ~ "Artesunate + Amodiaquine",
      treat == "AS+AQ-FDC" ~ "Artesunate–Amodiaquine (fixed-dose combination)",
      treat == "AS+AQ(1)FDC" ~ "Artesunate–Amodiaquine FDC (formulation 1)",
      treat == "AS+AQ(2)FDC" ~ "Artesunate–Amodiaquine FDC (formulation 2)",
      treat == "AS+SP" ~ "Artesunate + Sulfadoxine–Pyrimethamine",
      treat == "AS3+AQ-FDC" ~ "Artesunate–Amodiaquine (3-day FDC)",
      treat == "AS3+SP" ~ "Artesunate + Sulfadoxine–Pyrimethamine (3-day)",
      treat == "ASAQ" ~ "Artesunate + Amodiaquine (alternate code)",
      treat == "ASMQ" ~ "Artesunate + Mefloquine",
      treat == "ASMQ_adult" ~ "Artesunate + Mefloquine (adult dose)",
      treat == "ASMQ_child" ~ "Artesunate + Mefloquine (child dose)",
      treat == "AQ+SP" ~ "Amodiaquine + Sulfadoxine–Pyrimethamine",
      treat == "CQ" ~ "Chloroquine",
      treat == "CQ_double_dose" ~ "Chloroquine (double dose)",
      treat == "ART+NQ" ~ "Artemisinin + Naphthoquine",
      treat == "DHA-PPQ" ~ "Dihydroartemisinin–Piperaquine",
      treat == "DHAPQP" ~ "Dihydroartemisinin–Piperaquine (variant)",
      treat == "DHAPQPa" ~ "Dihydroartemisinin–Piperaquine (formulation a)",
      treat == "DHAPPQb" ~ "Dihydroartemisinin–Piperaquine (formulation b)",
      treat == "DP" ~ "Dihydroartemisinin–Piperaquine (DP)",
      treat == "SP" ~ "Sulfadoxine–Pyrimethamine",
      treat == "NA" ~ NA_character_,
      TRUE ~ "Other / Unclassified"
    ),
    act_category = case_when(
      treat %in% c(
        "AL", "AL_adult", "AL_child", "AL_PQ", "AL_unsupervised",
        "AS+AQ", "AS+AQ-FDC", "AS+AQ(1)FDC", "AS+AQ(2)FDC",
        "AS+SP", "AS3+AQ-FDC", "AS3+SP", "ASAQ",
        "ASMQ", "ASMQ_adult", "ASMQ_child",
        "ART+NQ", "DHA-PPQ", "DHAPQP", "DHAPQPa", "DHAPPQb", "DP"
      ) ~ "ACT",
      treat %in% c("CQ", "CQ_double_dose", "SP", "AQ+SP") ~ "Non-ACT",
      TRUE ~ "Unknown"
    )
  )


### Bar-chart to evaluate treatment distribution applying inclusion criteria
ggplot(data_2 %>% filter(age_5 == "<5y" & day_28 == "28" & act_category == "ACT"), aes(x = treat)) +
  geom_bar() +
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5, size = 3) +
  theme(
    axis.text.x = element_text(angle = 90)
  )


### 4 Arrange studies by descendant order
sid_arr_desc <- data_3 %>%
  tabyl(sid) %>%
  arrange(desc(n))

sid_arr_desc_2 <- data_3 %>%
  group_by(sid, pubmed) %>%
  summarize(n = n(), .groups = "drop") %>%
  arrange(desc(n))

sid_arr_desc_2 %>% print(n = 1000)

### 5. Outcome distribution bar graph

ggplot(data_3, aes(x = outcome)) +
  geom_bar() +
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5)
