library(tidyverse)
library(hrbrthemes)
library(eulerr)


raw <- read.delim("data/raw/hg19.txt") %>%
    janitor::clean_names()

threshold <- median(raw$ed_h3k27ac)

raw <- raw %>%
    mutate(across(contains("ed"), ~ .x / ed_input)) %>%
    mutate(across(contains("dox"), ~ .x / dox_input)) %>%
    mutate(fox_fold = (dox_foxa1_high - ed_foxa1) / ed_foxa1) %>%
    mutate(h3_fold = (dox_h3k27ac - ed_h3k27ac) / ed_h3k27ac) %>%
    mutate(fox_status = if_else(dox_foxa1_high > ed_foxa1, "gain", "loss")) %>%
    mutate(h3_status = ifelse(
        ed_h3k27ac < threshold & dox_h3k27ac > threshold,
        'gain',
        ifelse(ed_h3k27ac > threshold &
                   dox_h3k27ac < threshold, 'loss', 'shared')
    ))


upstream <- read.delim("data/raw/upstream.txt") %>%
    janitor::clean_names() %>%
    mutate(position = "upstream") %>%
    mutate(across(contains("ed"), ~ .x / ed_input)) %>%
    mutate(across(contains("dox"), ~ .x / dox_input)) %>%
    mutate(fox_fold = (dox_foxa1_high - ed_foxa1) / ed_foxa1) %>%
    mutate(h3_fold = (dox_h3k27ac - ed_h3k27ac) / ed_h3k27ac) %>%
    mutate(fox_status = if_else(dox_foxa1_high > ed_foxa1, "gain", "loss")) %>%
    mutate(h3_status = ifelse(
        ed_h3k27ac < threshold & dox_h3k27ac > threshold,
        'gain',
        ifelse(ed_h3k27ac > threshold &
                   dox_h3k27ac < threshold, 'loss', 'shared')
    ))

downstream <- read.delim("data/raw/downstream.txt") %>%
    janitor::clean_names() %>%
    mutate(position = "downstream")


total <- full_join(upstream,
                   downstream,
                   by = "name",
                   suffix = c(".up", ".dn"))


ggplot(raw, aes(x = log10(fox_fold), y = log10(h3_fold))) +
    geom_point(color = "black") +
    geom_smooth(
        method = lm ,
        color = "red",
        fill = "#69b3a2",
        se = TRUE
    ) +
    theme_ipsum()

ggplot(raw, aes(x = log10(dox_foxa1_high), y = log10(dox_h3k27ac))) +
    geom_point(color = "black") +
    geom_smooth(
        method = lm ,
        color = "red",
        fill = "#69b3a2",
        se = TRUE
    ) +
    theme_ipsum()


fox_gain <- length(grep("gain", raw$fox_status)) %>% as.numeric()
fox_loss <- length(grep("loss", raw$fox_status)) %>% as.numeric()

h3_gain <- length(grep("gain", raw$h3_status)) %>% as.numeric()
h3_loss <- length(grep("loss", raw$h3_status)) %>% as.numeric()

loss_loss <-
    raw %>% filter(fox_status == "loss" &
                       h3_status == "loss") %>% length() %>% as.numeric()
gain_gain <-
    raw %>% filter(fox_status == "gain" &
                       h3_status == "gain") %>% length() %>% as.numeric()




euler(c(A = fox_gain, B = h3_gain, "A&B" = gain_gain)) %>% plot


euler(c(A = fox_loss, B = h3_loss, "A&B" = loss_loss)) %>% plot





# foxa1
raw %>%
    pivot_longer(cols = contains(c("_foxa1", "_h3k27ac")),
                 names_to = "treatment") %>%
    filter(!grepl("h3k27ac|dox_foxa1_low", treatment)) %>%
    ggplot(aes(x = value, fill = treatment)) +
    geom_histogram(color = "#e9ecef",
                   alpha = 0.6,
                   position = 'identity') +
    scale_fill_manual(values = c("#69b3a2", "#404080")) +
    theme_ipsum() +
    labs(fill = "")


# h3k27ac
raw %>%
    pivot_longer(cols = contains(c("_foxa1", "_h3k27ac")),
                 names_to = "treatment") %>%
    filter(!grepl("foxa1", treatment)) %>%
    ggplot(aes(x = value, fill = treatment)) +
    geom_histogram(color = "#e9ecef",
                   alpha = 0.6,
                   position = 'identity') +
    scale_fill_manual(values = c("#69b3a2", "#404080")) +
    theme_ipsum() +
    labs(fill = "")
