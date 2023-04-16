#### Packages ####
library(tidyverse)
library(hrbrthemes)
library(ggthemes)
library(eulerr)
library(ggpubr)
library(purrr)
library(ggpmisc)
library(ggstatsplot)
library(colorspace)

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
    mutate(h3_status = ifelse(
        ed_h3k27ac_bound == "0" & dox_h3k27ac_bound == "1",
        "UP",
        ifelse(ed_h3k27ac_bound == "1" &
                   dox_h3k27ac_bound == "0",
               "DN",
               "NC")
    )) %>%
    mutate(fox_status = ifelse(
        ed_fox_bound == "0" & dox_fox_bound == "1",
        "UP",
        ifelse(ed_fox_bound == "1" &
                   dox_fox_bound == "0",
               "DN",
               'shared')
    )) %>%
    mutate(co_status = ifelse(
        ed_overlaid_bound == "0" & dox_overlaid_bound == "1",
        "UP",
        ifelse(
            ed_overlaid_bound == "1" &
                dox_overlaid_bound == "0",
            "DN",
            'shared'
        )
    )) %>%
    mutate(
        activity_status = ifelse(
            ed_h3k27ac < activity_threshold & dox_h3k27ac > activity_threshold,
            "GAIN",
            ifelse(
                ed_h3k27ac > activity_threshold &
                    dox_h3k27ac < activity_threshold,
                "LOSS",
                "NC"
            )
        )
    ) %>%
    pivot_longer(cols = contains(c("_foxa1", "_h3k27ac")) &
                     !contains("bound"),
                 names_to = "treatment") %>%
    filter(treatment != "dox_foxa1_low") %>%
    mutate(condition = str_extract(treatment, "[^_]+")) %>%
    mutate(amino = str_replace(name, "^(([^-]+-){1}[^-]+)-.*", "\\1")) %>%
    mutate(activity = ifelse(value > activity_threshold, '1', '0'))



#### fox volcano ####
volcano_fox <- data %>%
    pivot_wider(
        names_from = "treatment",
        values_from = "value",
        id_cols = c("name")
    ) %>%
    select(c(name, ed_foxa1, dox_foxa1_high)) %>%
    mutate(p_value = -log10(pmap_dbl(select(., -1),  ~ chisq.test(c(
        ...
    ))$p.value))) %>%
    mutate(fox_fold = log2(dox_foxa1_high / ed_foxa1)) %>%
    filter(fox_fold != "Inf") %>%
    distinct()

# add a column of NAs
volcano_fox$diffexpressed <- "NS"

# if log2Foldchange > 4 and pvalue < .000000001, set as "UP"
volcano_fox$diffexpressed[volcano_fox$fox_fold > log2(1.5) &
                              volcano_fox$p_value > .001] <-
    "UP"

# if log2Foldchange < -4 and pvalue < .000000001, set as "DN"
volcano_fox$diffexpressed[volcano_fox$fox_fold < -log2(1.5) &
                              volcano_fox$p_value > .001] <-
    "DN"

#### h3 volcano ####
volcano_h3 <- data %>%
    pivot_wider(
        names_from = "treatment",
        values_from = "value",
        id_cols = c("name")
    ) %>%
    select(c(name, ed_h3k27ac, dox_h3k27ac)) %>%
    mutate(p_value = -log10(pmap_dbl(select(., -1),  ~ chisq.test(c(
        ...
    ))$p.value))) %>%
    mutate(h3_fold = log2(dox_h3k27ac / ed_h3k27ac)) %>%
    filter(h3_fold != "Inf") %>%
    distinct()


# add a column of NAs
volcano_h3$diffexpressed <- "NS"

# if log2Foldchange > 4 and pvalue < .000000001, set as "UP"
volcano_h3$diffexpressed[volcano_h3$h3_fold > log2(1.5) &
                             volcano_h3$p_value > .001] <-
    "UP"

# if log2Foldchange < -4 and pvalue < .000000001, set as "DN"
volcano_h3$diffexpressed[volcano_h3$h3_fold < -log2(1.5) &
                             volcano_h3$p_value > .001] <-
    "DN"


diffexpressed <-
    full_join(volcano_fox, volcano_h3, by = "name", copy = T) %>%
    pivot_longer(cols = contains(c("diffexpressed.x", "diffexpressed.y")),
                 names_to = "comparison",
                 values_to = "diffexpressed") %>%
    mutate(comparison = str_replace(comparison, "diffexpressed.x", "fox_diff")) %>%
    mutate(comparison = str_replace(comparison, "diffexpressed.y", "h3_diff")) %>%
    distinct() 


#### Venn ####

fox_loss <-
    length(grep("DN", volcano_fox$diffexpressed)) %>% as.numeric()

h3_loss <-
    length(grep("DN", volcano_h3$diffexpressed)) %>% as.numeric()

x <- volcano_fox %>%
    filter(diffexpressed == "DN") %>%
    select(name)

y <- volcano_h3 %>%
    filter(diffexpressed == "DN") %>%
    select(name)

loss_loss <- inner_join(x, y) %>% nrow() %>% as.numeric()




fox_gain <-
    length(grep("UP", volcano_fox$diffexpressed)) %>% as.numeric()
h3_gain <-
    length(grep("UP", volcano_h3$diffexpressed)) %>% as.numeric()

x <- volcano_fox %>%
    filter(diffexpressed == "UP") %>%
    select(name)

y <- volcano_h3 %>%
    filter(diffexpressed == "UP") %>%
    select(name)

gain_gain <- inner_join(x, y) %>% nrow() %>% as.numeric()



x <- volcano_fox %>%
    filter(diffexpressed == "UP") %>%
    select(name)

y <- volcano_h3 %>%
    filter(diffexpressed == "DN") %>%
    select(name)

gain_loss <- inner_join(x, y) %>% nrow() %>% as.numeric()


x <- volcano_fox %>%
    filter(diffexpressed == "DN") %>%
    select(name)

y <- volcano_h3 %>%
    filter(diffexpressed == "UP") %>%
    select(name)

loss_gain <- inner_join(x, y) %>% nrow() %>% as.numeric()




#### relative position ####
upstream <- read.delim("data/raw/upstream.txt") %>%
    janitor::clean_names() %>%
    mutate(position = "upstream")

downstream <- read.delim("data/raw/downstream.txt") %>%
    janitor::clean_names() %>%
    mutate(position = "downstream")

total_activity_threshold <- full_join(upstream,
                                      downstream) %>%
    mutate(across(contains("ed"), ~ .x / ed_input)) %>%
    mutate(across(contains("dox"), ~ .x / dox_input)) %>%
    pivot_longer(cols = contains(c("_foxa1", "_h3k27ac")),
                 names_to = "treatment") %>%
    group_by(treatment) %>%
    summarise_at("value", median, na.rm = TRUE) %>%
    filter(treatment == "ed_h3k27ac") %>%
    pull(2)


total <- full_join(upstream,
                   downstream) %>%
    mutate(across(contains("ed"), ~ .x / ed_input)) %>%
    mutate(across(contains("dox"), ~ .x / dox_input)) %>%
    filter_all(all_vars(!is.infinite(.))) %>%
    mutate(ed_fox_bound = ifelse(ed_foxa1 > 1, '1', '0')) %>%
    mutate(dox_h3k27ac_bound = ifelse(dox_h3k27ac > 1, '1', '0')) %>%
    mutate(ed_h3k27ac_bound = ifelse(ed_h3k27ac > 1, '1', '0')) %>%
    mutate(dox_fox_bound = ifelse(dox_foxa1_high > 1, '1', '0')) %>%
    mutate(ed_overlaid_bound = ifelse(ed_foxa1 > 1 &
                                          ed_h3k27ac > 1, '1', '0')) %>%
    mutate(dox_overlaid_bound = ifelse(dox_foxa1_high > 1 &
                                           dox_h3k27ac > 1, '1', '0')) %>%
    mutate(h3_status = ifelse(
        ed_h3k27ac_bound == "0" & dox_h3k27ac_bound == "1",
        "GAIN",
        ifelse(ed_h3k27ac_bound == "1" &
                   dox_h3k27ac_bound == "0",
               "LOSS",
               "NC")
    )) %>%
    mutate(fox_status = ifelse(
        ed_fox_bound == "0" & dox_fox_bound == "1",
        "UP",
        ifelse(ed_fox_bound == "1" &
                   dox_fox_bound == "0",
               "DN",
               'shared')
    )) %>%
    mutate(co_status = ifelse(
        ed_overlaid_bound == "0" & dox_overlaid_bound == "1",
        "UP",
        ifelse(
            ed_overlaid_bound == "1" &
                dox_overlaid_bound == "0",
            "DN",
            'shared'
        )
    )) %>%
    mutate(
        activity_status = ifelse(
            ed_h3k27ac < total_activity_threshold & dox_h3k27ac > total_activity_threshold,
            "GAIN",
            ifelse(
                ed_h3k27ac > total_activity_threshold &
                    dox_h3k27ac < total_activity_threshold,
                "LOSS",
                "NC"
            )
        )
    ) %>%
    pivot_longer(cols = contains(c("_foxa1", "_h3k27ac")) &
                     !contains("bound"),
                 names_to = "treatment") %>%
    mutate(condition = str_extract(treatment, "[^_]+")) %>%
    filter(treatment != "dox_foxa1_low") %>%
    distinct() %>%
    pivot_longer(cols = contains("bound"),
                 names_to = "target",
                 values_to = "bound") %>%
    drop_na(bound) %>%
    mutate(activity = ifelse(value > total_activity_threshold, '1', '0'))



save.image(file = 'environments/data.RData')

