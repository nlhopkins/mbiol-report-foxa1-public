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
    mutate(across(contains("ed"), ~ .x / ed_input)) %>%
    mutate(across(contains("dox"), ~ .x / dox_input)) %>%
    pivot_longer(cols = contains(c("_foxa1", "_h3k27ac")),
                 names_to = "treatment") %>%
    group_by(treatment) %>%
    summarise_at("value", median) %>%
    filter(treatment == "ed_h3k27ac") %>%
    pull(2)

# ed_biding_threshold <- raw %>%
#     pivot_longer(cols = contains(c("input")),
#                  names_to = "input",
#                  values_to = "input_value") %>%
#     group_by(input) %>%
#     summarise_at("input_value", median) %>%
#     filter(input == "ed_input") %>%
#     pull(2)
#
# dox_biding_threshold <- raw %>%
#     pivot_longer(cols = contains(c("input")),
#                  names_to = "input",
#                  values_to = "input_value") %>%
#     group_by(input) %>%
#     summarise_at("input_value", median) %>%
#     filter(input == "dox_input") %>%
#     pull(2)


data <- raw %>%
    mutate(across(contains("ed"), ~ .x / ed_input)) %>%
    mutate(across(contains("dox"), ~ .x / dox_input)) %>%
    mutate(ed_fox_bound = ifelse(ed_foxa1 > 1, '1', '0')) %>%
    mutate(dox_h3k27ac_bound = ifelse(dox_h3k27ac > 1, '1', '0')) %>%
    mutate(ed_h3k27ac_bound = ifelse(ed_h3k27ac > 1, '1', '0')) %>%
    mutate(dox_fox_bound = ifelse(dox_foxa1_high > 1, '1', '0')) %>%
    mutate(ed_overlaid_bound = ifelse(ed_foxa1 > 1 &
                                          ed_h3k27ac > 1, '1', '0')) %>%
    mutate(dox_overlaid_bound = ifelse(dox_foxa1_high > 1 &
                                           dox_h3k27ac > 1, '1', '0')) %>%
    mutate(fox_fold = dox_foxa1_high / ed_foxa1) %>%
    mutate(h3_fold = dox_h3k27ac / ed_h3k27ac)  %>%
    mutate(
        h3_status = ifelse(
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
    filter(treatment != "dox_foxa1_low") %>%
    mutate(condition = str_extract(treatment, "[^_]+")) %>%
    mutate(amino = str_replace(name, "^(([^-]+-){1}[^-]+)-.*", "\\1"))



means <- data %>% group_by(treatment) %>%
    summarise_at("value", mean)


#### fox volcano ####
volcano_fox <- data %>%
    pivot_wider(
        names_from = "treatment",
        values_from = "value",
        id_cols = c("name")
    ) %>%
    select(c(name, ed_foxa1, dox_foxa1_high)) %>%
    filter(ed_foxa1 > 1 |
               dox_foxa1_high > 1) %>%
    mutate(p_value = -log10(pmap_dbl(select(., -1),  ~ chisq.test(c(
        ...
    ))$p.value))) %>%
    mutate(fox_fold = log2(dox_foxa1_high / ed_foxa1)) %>%
    filter(fox_fold != "Inf") %>%
    unique()

# add a column of NAs
volcano_fox$diffexpressed <- "NS"

# if log2Foldchange > 4 and pvalue < .000000001, set as "GAIN"
volcano_fox$diffexpressed[volcano_fox$fox_fold > log2(1.5) &
                              volcano_fox$p_value > .001] <-
    "GAIN"

# if log2Foldchange < -4 and pvalue < .000000001, set as "LOSS"
volcano_fox$diffexpressed[volcano_fox$fox_fold < -log2(1.5) &
                              volcano_fox$p_value > .001] <-
    "LOSS"

#### h3 volcano ####
volcano_h3 <- data %>%
    pivot_wider(
        names_from = "treatment",
        values_from = "value",
        id_cols = c("name")
    ) %>%
    select(c(name, ed_h3k27ac, dox_h3k27ac)) %>%
    filter(ed_h3k27ac > 1 |
               dox_h3k27ac > 1) %>%
    mutate(p_value = -log10(pmap_dbl(select(., -1),  ~ chisq.test(c(
        ...
    ))$p.value))) %>%
    mutate(h3_fold = log2(dox_h3k27ac / ed_h3k27ac)) %>%
    filter(h3_fold != "Inf") %>%
    unique()


# add a column of NAs
volcano_h3$diffexpressed <- "NS"

# if log2Foldchange > 4 and pvalue < .000000001, set as "GAIN"
volcano_h3$diffexpressed[volcano_h3$h3_fold > log2(1.5) &
                             volcano_h3$p_value > .001] <-
    "GAIN"

# if log2Foldchange < -4 and pvalue < .000000001, set as "LOSS"
volcano_h3$diffexpressed[volcano_h3$h3_fold < -log2(1.5) &
                             volcano_h3$p_value > .001] <-
    "LOSS"


diffexpressed <-
    merge(volcano_fox, volcano_h3, by = "name") %>%
    pivot_longer(cols = contains(c("diffexpressed.x", "diffexpressed.y")),
                 names_to = "comparison",
                 values_to = "diffexpressed") %>%
    select(c("name", "comparison", "diffexpressed")) %>%
    mutate(comparison = str_replace(comparison, "diffexpressed.x", "fox_diff")) %>%
    mutate(comparison = str_replace(comparison, "diffexpressed.y", "h3_diff")) %>%
    distinct()


#### Venn ####

fox_loss <-
    length(grep("LOSS", volcano_fox$diffexpressed)) %>% as.numeric()

h3_loss <-
    length(grep("LOSS", volcano_h3$diffexpressed)) %>% as.numeric()

x <- volcano_fox %>%
    filter(diffexpressed == "LOSS") %>%
    select(name)

y <- volcano_h3 %>%
    filter(diffexpressed == "LOSS") %>%
    select(name)

loss_loss <- inner_join(x, y) %>% nrow() %>% as.numeric()




fox_gain <-
    length(grep("GAIN", volcano_fox$diffexpressed)) %>% as.numeric()
h3_gain <-
    length(grep("GAIN", volcano_h3$diffexpressed)) %>% as.numeric()

x <- volcano_fox %>%
    filter(diffexpressed == "GAIN") %>%
    select(name)

y <- volcano_h3 %>%
    filter(diffexpressed == "GAIN") %>%
    select(name)

gain_gain <- inner_join(x, y) %>% nrow() %>% as.numeric()



x <- volcano_fox %>%
    filter(diffexpressed == "GAIN") %>%
    select(name)

y <- volcano_h3 %>%
    filter(diffexpressed == "LOSS") %>%
    select(name)

gain_loss <- inner_join(x, y) %>% nrow() %>% as.numeric()


x <- volcano_fox %>%
    filter(diffexpressed == "LOSS") %>%
    select(name)

y <- volcano_h3 %>%
    filter(diffexpressed == "GAIN") %>%
    select(name)

loss_gain <- inner_join(x, y) %>% nrow() %>% as.numeric()




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


total <- total

total <- merge(diffexpressed, total, by = "name")

save.image(file = 'environments/data.RData')

