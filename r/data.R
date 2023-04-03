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

activity_threshold <- raw %>%
    pivot_longer(cols = contains(c("_foxa1", "_h3k27ac")),
                 names_to = "treatment") %>%
    group_by(treatment) %>%
    summarise_at("value", median) %>%
    filter(treatment == "ed_h3k27ac") %>%
    pull(2)

ed_biding_threshold <- raw %>%
    pivot_longer(cols = contains(c("input")),
                 names_to = "input",
                 values_to = "input_value") %>%
    group_by(input) %>%
    summarise_at("input_value", median) %>%
    filter(input == "ed_input") %>%
    pull(2)

dox_biding_threshold <- raw %>%
    pivot_longer(cols = contains(c("input")),
                 names_to = "input",
                 values_to = "input_value") %>%
    group_by(input) %>%
    summarise_at("input_value", median) %>%
    filter(input == "dox_input") %>%
    pull(2)


data <- raw %>%
    mutate(dox_fox_bound = ifelse(dox_foxa1_high > dox_biding_threshold, '1', '0')) %>%
    mutate(ed_fox_bound = ifelse(ed_foxa1 > ed_biding_threshold, '1', '0')) %>%
    mutate(dox_h3k27ac_bound = ifelse(dox_h3k27ac > dox_biding_threshold, '1', '0')) %>%
    mutate(ed_h3k27ac_bound = ifelse(ed_h3k27ac > ed_biding_threshold, '1', '0')) %>%
    mutate(fox_fold = (dox_foxa1_high - ed_foxa1) / ed_foxa1) %>%
    mutate(h3_fold = (dox_h3k27ac - ed_h3k27ac) / ed_h3k27ac)  %>%
    mutate(fox_abs = dox_foxa1_high - ed_foxa1) %>%
    mutate(h3_abs = dox_h3k27ac - ed_h3k27ac) %>%
    mutate(
        activity = ifelse(
            ed_h3k27ac < activity_threshold & dox_h3k27ac > activity_threshold,
            'gain',
            ifelse(
                ed_h3k27ac > activity_threshold &
                    dox_h3k27ac < activity_threshold,
                'loss',
                'shared'
            )
        )
    ) %>%
    pivot_longer(cols = contains(c("_foxa1", "_h3k27ac")) &
                     !contains("bound"),
                 names_to = "treatment") %>%
    mutate(condition = str_extract(treatment, "[^_]+")) %>%
    mutate(amino = str_replace(name, "^(([^-]+-){1}[^-]+)-.*", "\\1")) %>%
    filter(treatment !="dox_foxa1_low")



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
volcano_fox$diffexpressed[volcano_fox$fox_fold > 0.25 &
                              volcano_fox$p_value > .01] <-
    "UP"

# if log2Foldchange < -4 and pvalue < .000000001, set as "DOWN"
volcano_fox$diffexpressed[volcano_fox$fox_fold < -0.25 &
                              volcano_fox$p_value > .01] <-
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
volcano_h3$diffexpressed[volcano_h3$h3_fold > 0.25 &
                             volcano_h3$p_value > .01] <-
    "UP"

# if log2Foldchange < -4 and pvalue < .000000001, set as "DOWN"
volcano_h3$diffexpressed[volcano_h3$h3_fold < -0.25 &
                             volcano_h3$p_value > .01] <-
    "DOWN"


diffexpressed <-
    merge(volcano_fox, volcano_h3, by = "name") %>%
    pivot_longer(cols = contains(c("diffexpressed.x", "diffexpressed.y")),
                 names_to = "comparison",
                 values_to = "diffexpressed") %>%
    select(c("name", "comparison", "diffexpressed")) %>%
    mutate(comparison = str_replace(comparison, "diffexpressed.x", "fox_diff")) %>%
    mutate(comparison = str_replace(comparison, "diffexpressed.y", "h3_diff")) %>%
    distinct()

data <- merge(diffexpressed, data, by = "name")


#### Venn ####

fox_loss <-
    length(grep("DOWN", volcano_fox$diffexpressed)) %>% as.numeric()
h3_loss <-
    length(grep("DOWN", volcano_h3$diffexpressed)) %>% as.numeric()

loss_loss <- fox_loss + h3_loss

fox_gain <-
    length(grep("UP", volcano_fox$diffexpressed)) %>% as.numeric()
h3_gain <-
    length(grep("UP", volcano_h3$diffexpressed)) %>% as.numeric()

gain_gain <- fox_gain + h3_gain


gain_loss <- fox_gain + h3_loss
loss_gain <- fox_loss + h3_gain


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
                   downstream)

pos_activity_threshold <- total %>%
    pivot_longer(cols = contains(c("_foxa1", "_h3k27ac")),
                 names_to = "treatment") %>%
    group_by(treatment) %>%
    summarise_at("value", median) %>%
    filter(treatment == "ed_h3k27ac")

pos_ed_biding_threshold <- total %>%
    pivot_longer(cols = contains(c("input")),
                 names_to = "input",
                 values_to = "input_value") %>%
    group_by(input) %>%
    summarise_at("input_value", median) %>%
    filter(input == "ed_input") %>%
    pull(2)

pos_dox_biding_threshold <- total %>%
    pivot_longer(cols = contains(c("input")),
                 names_to = "input",
                 values_to = "input_value") %>%
    group_by(input) %>%
    summarise_at("input_value", median) %>%
    filter(input == "ed_input") %>%
    pull(2)


total <- total %>%
    mutate(dox_fox_bound = ifelse(dox_foxa1_high > dox_biding_threshold, '1', '0')) %>%
    mutate(ed_fox_bound = ifelse(ed_foxa1 > ed_biding_threshold, '1', '0')) %>%
    mutate(dox_h3k27ac_bound = ifelse(dox_h3k27ac > pos_dox_biding_threshold, '1', '0')) %>%
    mutate(ed_h3k27ac_bound = ifelse(ed_h3k27ac > pos_ed_biding_threshold, '1', '0')) %>%
    mutate(fox_fold = (dox_foxa1_high - ed_foxa1) / ed_foxa1) %>%
    mutate(h3_fold = (dox_h3k27ac - ed_h3k27ac) / ed_h3k27ac)  %>%
    mutate(fox_abs = dox_foxa1_high - ed_foxa1) %>%
    mutate(h3_abs = dox_h3k27ac - ed_h3k27ac) %>%
    mutate(
        activity = ifelse(
            ed_h3k27ac < pos_activity_threshold &
                dox_h3k27ac > pos_activity_threshold,
            'gain',
            ifelse(
                ed_h3k27ac > pos_activity_threshold &
                    dox_h3k27ac < pos_activity_threshold,
                'loss',
                'shared'
            )
        )
    ) %>%
    mutate(h3_status = ifelse(
        dox_h3k27ac > pos_activity_threshold,
        'active',
        ifelse(ed_h3k27ac > pos_activity_threshold, 'active', 'inactive')
    )) %>%
    pivot_longer(cols = contains(c("_foxa1", "_h3k27ac")) &
                     !contains("bound"),
                 names_to = "treatment") %>%
    mutate(condition = str_extract(treatment, "[^_]+")) %>%
    mutate(amino = str_replace(name, "^(([^-]+-){1}[^-]+)-.*", "\\1"))

total <- merge(diffexpressed, total, by = "name")

save.image(file = 'environments/data.RData')
