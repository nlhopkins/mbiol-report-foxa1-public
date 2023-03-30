#### Packages ####
library(tidyverse)
library(hrbrthemes)
library(ggthemes)
library(eulerr)
library(ggpubr)
library(purrr)
library(ggpmisc)
library(ggstatsplot)

options(scipen = 999)

#### Tidy ####
raw <- read.delim("data/raw/hg19.txt") %>%
    janitor::clean_names()

threshold <- raw %>%
    mutate(across(contains("ed"), ~ .x / ed_input)) %>%
    mutate(across(contains("dox"), ~ .x / dox_input)) %>%
    pivot_longer(cols = contains(c("_foxa1", "_h3k27ac")),
                 names_to = "treatment") %>%
    group_by(treatment) %>%
    summarise_at("value", median) %>%
    filter(treatment == "ed_h3k27ac") %>%
    pull(2)

data <- raw %>%
    mutate(across(contains("ed"), ~ .x / ed_input)) %>%
    mutate(across(contains("dox"), ~ .x / dox_input)) %>%
    mutate(fox_fold = (dox_foxa1_high - ed_foxa1) / ed_foxa1) %>%
    mutate(h3_fold = (dox_h3k27ac - ed_h3k27ac) / ed_h3k27ac)  %>%
    mutate(fox_abs = dox_foxa1_high - ed_foxa1) %>%
    mutate(h3_abs = dox_h3k27ac - ed_h3k27ac) %>%
    pivot_longer(cols = contains(c("_foxa1", "_h3k27ac")),
                 names_to = "treatment") %>%
    mutate(condition = str_extract(treatment, "[^_]+"))


# mutate(activity = ifelse(
#     ed_h3k27ac < threshold & dox_h3k27ac > threshold,
#     'gain',
#     ifelse(ed_h3k27ac > threshold &
#                dox_h3k27ac < threshold, 'loss', 'shared')
# )) %>%

means <- data %>% group_by(treatment) %>%
    summarise_at("value", mean)


#### fox volcano ####
volcano_fox <- data %>%
    pivot_wider(names_from = "treatment",
                values_from = "value",
                id_cols = "name") %>%
    select(c(name, ed_foxa1, dox_foxa1_high)) %>%
    mutate(p_value = -log10(pmap_dbl(select(., -1),  ~ chisq.test(c(
        ...
    ))$p.value))) %>%
    mutate(fox_fold = ((dox_foxa1_high - ed_foxa1) / ed_foxa1)) %>%
    filter(fox_fold != "Inf")

# add a column of NAs
volcano_fox$diffexpressed <- "NS"

# if log2Foldchange > 4 and pvalue < .000000001, set as "UP"
volcano_fox$diffexpressed[volcano_fox$fox_fold > 0.5 &
                              volcano_fox$p_value > .000000001] <-
    "UP"

# if log2Foldchange < -4 and pvalue < .000000001, set as "DOWN"
volcano_fox$diffexpressed[volcano_fox$fox_fold < -0.5 &
                              volcano_fox$p_value > .000000001] <-
    "DOWN"

#### h3 volcano ####
volcano_h3 <- data %>%
    pivot_wider(names_from = "treatment",
                values_from = "value",
                id_cols = "name") %>%
    select(c(name, ed_h3k27ac, dox_h3k27ac)) %>%
    mutate(p_value = -log10(pmap_dbl(select(., -1),  ~ chisq.test(c(
        ...
    ))$p.value))) %>%
    mutate(h3_fold = ((dox_h3k27ac - ed_h3k27ac) / ed_h3k27ac)) %>%
    filter(h3_fold != "Inf")


# add a column of NAs
volcano_h3$diffexpressed <- "NS"

# if log2Foldchange > 4 and pvalue < .000000001, set as "UP"
volcano_h3$diffexpressed[volcano_h3$h3_fold > 0.5 &
                             volcano_h3$p_value > .000000001] <-
    "UP"

# if log2Foldchange < -4 and pvalue < .000000001, set as "DOWN"
volcano_h3$diffexpressed[volcano_h3$h3_fold < -0.5 &
                             volcano_h3$p_value > .000000001] <-
    "DOWN"


diffexpressed <-
    merge(volcano_fox, volcano_h3, by = "name") %>%
    pivot_longer(cols = contains(c("diffexpressed.x", "diffexpressed.y")),
                 names_to = "comparison",
                 values_to = "diffexpressed") %>%
    select(c("name", "comparison", "diffexpressed")) %>%
    mutate(comparison = str_replace(comparison, "diffexpressed.x", "fox")) %>%
    mutate(comparison = str_replace(comparison, "diffexpressed.y", "h3")) %>%
    distinct()

data2 <- merge(diffexpressed, data, by = "name")


#### Venn ####

fox_loss <-
    length(grep("DOWN", volcano_fox$diffexpressed)) %>% as.numeric()
h3_loss <-
    length(grep("DOWN", volcano_h3$diffexpressed)) %>% as.numeric()
loss_loss <-
    sum(volcano_fox$diffexpressed == "DOWN" &
            grepl("DOWN", volcano_h3$diffexpressed))

fox_gain <-
    length(grep("UP", volcano_fox$diffexpressed)) %>% as.numeric()
h3_gain <-
    length(grep("UP", volcano_h3$diffexpressed)) %>% as.numeric()
gain_gain <-
    sum(volcano_fox$diffexpressed == "UP" &
            grepl("UP", volcano_h3$diffexpressed))


gain_loss <-
    sum(volcano_fox$diffexpressed == "UP" &
            grepl("DOWN", volcano_h3$diffexpressed))
loss_gain <-
    sum(volcano_fox$diffexpressed == "DOWN" &
            grepl("UP", volcano_h3$diffexpressed))


up <- merge(volcano_fox, volcano_h3, by = "name") %>%
    filter(diffexpressed.x == "UP" &
               diffexpressed.y == "UP")

dn <- merge(volcano_fox, volcano_h3, by = "name") %>%
    filter(diffexpressed.x == "DOWN" &
               diffexpressed.y == "DOWN")


#### relative position ####
upstream <- read.delim("data/raw/upstream.txt") %>%
    janitor::clean_names() %>%
    mutate(position = "upstream")

downstream <- read.delim("data/raw/downstream.txt") %>%
    janitor::clean_names() %>%
    mutate(position = "downstream")


total <- full_join(upstream,
                   downstream) %>%
    mutate(across(contains("ed"), ~ .x / ed_input)) %>%
    mutate(across(contains("dox"), ~ .x / dox_input)) %>%
    pivot_longer(cols = contains(c("_foxa1", "_h3k27ac")),
                 names_to = "treatment")

total <- merge(diffexpressed, total, by = "name")

total_threshold <- data %>%
    group_by(treatment) %>%
    summarise_at(vars(value), list(name = median)) %>%
    filter(treatment == "ed_h3k27ac") %>%
    pull(2)

save.image(file = 'environments/data.RData')
