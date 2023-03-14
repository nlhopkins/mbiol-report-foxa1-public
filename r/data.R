library(tidyverse)
library(hrbrthemes)
library(eulerr)

raw <- read.delim("data/raw/hg19.txt") %>%
    janitor::clean_names() %>%
    mutate(across(contains("ed"), ~ .x / ed_input)) %>%
    mutate(across(contains("dox"), ~ .x / dox_input)) %>%
    mutate(fox_fold = (dox_foxa1_high - ed_foxa1) / ed_foxa1) %>%
    mutate(h3_fold = (dox_h3k27ac - ed_h3k27ac) / ed_h3k27ac)

threshold <- median(raw$ed_h3k27ac)

raw <- raw %>%
    mutate(fox_status = if_else(dox_foxa1_high > ed_foxa1, "gain", "loss")) %>%
    mutate(h3_status = ifelse(
        ed_h3k27ac < threshold & dox_h3k27ac > threshold,
        'gain',
        ifelse(ed_h3k27ac > threshold &
                   dox_h3k27ac < threshold, 'loss', 'shared')
    ))

library(eulerr)
fit <- euler(c(A = 450, B = 1800, "A&B" = 230))
plot(fit)

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
