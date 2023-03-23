library(tidyverse)
library(hrbrthemes)
library(ggthemes)
library(eulerr)
library(ggpubr)

raw <- read.delim("data/raw/hg19.txt") %>%
    janitor::clean_names()

threshold <- median(raw$ed_h3k27ac, na.rm = TRUE)

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
    )) %>%
    mutate(fox_abs = dox_foxa1_high - ed_foxa1) %>%
    mutate(h3_abs = dox_h3k27ac - ed_h3k27ac)


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
    mutate(fox_fold = (dox_foxa1_high - ed_foxa1) / ed_foxa1) %>%
    mutate(h3_fold = (dox_h3k27ac - ed_h3k27ac) / ed_h3k27ac) %>%
    mutate(fox_status = if_else(dox_foxa1_high > ed_foxa1, "gain", "loss")) %>%
    mutate(h3_status = ifelse(
        ed_h3k27ac < threshold & dox_h3k27ac > threshold,
        'gain',
        ifelse(ed_h3k27ac > threshold &
                   dox_h3k27ac < threshold, 'loss', 'shared')
    ))


pos_threshold <- median(total$ed_h3k27ac, na.rm = TRUE)

# fold change
ggplot(raw, aes(x = log10(fox_fold), y = log10(h3_fold))) +
    geom_point(color = "black") +
    geom_smooth(
        method = lm ,
        color = "red",
        fill = "#69b3a2",
        se = TRUE
    ) +
    theme_classic()

# abs values
ggplot(raw, aes(x = log10(fox_abs), y = log10(h3_abs))) +
    geom_point(color = "black") +
    geom_smooth(
        method = lm ,
        color = "red",
        fill = "#69b3a2",
        se = TRUE
    ) +
    theme_classic()

# q values
ggplot(raw, aes(x = log10(dox_foxa1_high), y = log10(dox_h3k27ac))) +
    geom_point(color = "black") +
    geom_smooth(
        method = lm ,
        color = "red",
        fill = "#69b3a2",
        se = TRUE
    ) +
    theme_classic()


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
raw   %>%
    pivot_longer(cols = contains(c("_foxa1", "_h3k27ac")),
                 names_to = "treatment") %>%
    filter(!grepl("h3k27ac|dox_foxa1_low", treatment)) %>%
    ggplot(aes(x = value, fill = treatment)) +
    geom_histogram(
        color = "#e9ecef",
        alpha = 0.6,
        position = 'identity',
        bins = 50
    ) +
    scale_fill_manual(values = c("#69b3a2", "#404080")) +
    theme_ipsum() +
    labs(fill = "")


# h3k27ac
raw  %>%
    pivot_longer(cols = contains(c("_foxa1", "_h3k27ac")),
                 names_to = "treatment") %>%
    filter(!grepl("foxa1", treatment)) %>%
    ggplot(aes(x = value, fill = treatment)) +
    geom_histogram(
        color = "#e9ecef",
        alpha = 0.6,
        position = 'identity',
        bins = 50
    ) +
    scale_fill_manual(values = c("#69b3a2", "#404080")) +
    theme_ipsum() +
    labs(fill = "")


# position
# foxa1
total %>%
    pivot_longer(cols = contains(c("_foxa1", "_h3k27ac")),
                 names_to = "treatment") %>%
    filter(!grepl("h3k27ac|dox_foxa1_low", treatment)) %>%
    ggplot(aes(x = value, fill = position)) +
    geom_histogram(
        color = "#e9ecef",
        alpha = 0.6,
        position = 'identity',
        bins = 50
    ) +
    scale_fill_manual(values = c("#69b3a2", "#404080")) +
    theme_ipsum() +
    labs(fill = "")

# h3
total  %>%
    pivot_longer(cols = contains(c("_foxa1", "_h3k27ac")),
                 names_to = "treatment") %>%
    filter(!grepl("foxa1", treatment)) %>%
    ggplot(aes(x = value, fill = treatment)) +
    geom_histogram(
        color = "#e9ecef",
        alpha = 0.6,
        position = 'identity',
        bins = 50
    ) +
    scale_fill_manual(values = c("#69b3a2", "#404080")) +
    theme_ipsum() +
    labs(fill = "")

# fold
total  %>%
    pivot_longer(cols = contains(c("_foxa1", "_h3k27ac")),
                 names_to = "treatment") %>%
    filter(!grepl("h3k27ac|dox_foxa1_low", treatment)) %>%
    ggplot(aes(x = fox_fold, fill = position)) +
    geom_histogram(
        color = "#e9ecef",
        alpha = 0.6,
        position = 'identity',
        bins = 50
    ) +
    scale_fill_manual(values = c("#69b3a2", "#404080")) +
    theme_ipsum() +
    labs(fill = "")

total  %>%
    pivot_longer(cols = contains(c("_foxa1", "_h3k27ac")),
                 names_to = "treatment") %>%
    filter(!grepl("foxa1", treatment)) %>%
    ggplot(aes(x = h3_fold, fill = position)) +
    geom_histogram(
        color = "#e9ecef",
        alpha = 0.6,
        position = 'identity',
        bins = 50
    ) +
    scale_fill_manual(values = c("#69b3a2", "#404080")) +
    theme_ipsum() +
    labs(fill = "")


total  %>%
    pivot_longer(cols = contains(c("_foxa1", "_h3k27ac")),
                 names_to = "treatment") %>%
    pivot_longer(cols = contains(c("_foxa1", "_h3k27ac")),
                 names_to = "treatment") %>%
    filter(!grepl("h3k27ac|dox_foxa1_low", treatment)) %>%
    ggplot(aes(x = treatment, y = log10(value), fill = position)) +
    geom_boxplot() +
    theme_ipsum()


total %>%
    pivot_longer(cols = contains(c("_foxa1", "_h3k27ac")),
                 names_to = "treatment") %>%
    mutate(activity = if_else(value > pos_threshold, "active", "inactive")) %>%
    filter(!grepl("h3k27ac|dox_foxa1_low", treatment)) %>%
    filter(!grepl("inactive", activity)) %>%
    ggplot(aes(
        fill = factor(position, level = c("upstream", "downstream")),
        x = factor(treatment,
                   level = c("ed_foxa1", "dox_foxa1_high")),
        y = after_stat(count)
    )) +
    geom_bar(position = "dodge") +
    xlab("") +
    scale_fill_discrete(name = "") +
    theme_ipsum()



# whisker plots
x <- raw %>%
    pivot_longer(cols = contains(c("_foxa1", "_h3k27ac")),
                 names_to = "treatment") %>%
    mutate(target = case_when(str_detect(treatment, "fox") ~ "fox", TRUE ~ "h3")) %>%
    mutate(variable = case_when(str_detect(treatment, "dox") ~ "dox", TRUE ~ "ed"))

mod <- aov(data = x, value ~ target * variable)
summary(mod)
TukeyHSD(mod)
plot(mod, which = 2)
plot(mod, which = 1)


raw  %>%
    pivot_longer(cols = contains(c("_foxa1", "_h3k27ac")),
                 names_to = "treatment") %>%
    mutate(target = case_when(str_detect(treatment, "fox") ~ "fox", TRUE ~ "h3")) %>%
    mutate(variable = case_when(str_detect(treatment, "dox") ~ "dox", TRUE ~ "ed"))  %>%
    filter(!grepl("dox_foxa1_low", treatment)) %>%
    ggplot(aes(x = target,
               y = log10(value),
               fill = variable)) +
    geom_boxplot() +
    theme_ipsum() +
    stat_compare_means(method = "wilcox.test",
                       label = "p.signif",
                       paired = T)


raw  %>%
    pivot_longer(cols = contains(c("_foxa1", "_h3k27ac")),
                 names_to = "treatment") %>%
    mutate(target = case_when(str_detect(treatment, "fox") ~ "fox", TRUE ~ "h3")) %>%
    mutate(variable = case_when(str_detect(treatment, "dox") ~ "dox", TRUE ~ "ed"))  %>%
    filter(!grepl("dox_foxa1_low", treatment)) %>% 
    ggplot(aes(treatment, name, fill = log10(value))) +
    geom_tile()
