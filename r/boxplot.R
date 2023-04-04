load('environments/data.RData')

#### whisker plots ####
model <- data %>%
    mutate(target = case_when(str_detect(treatment, "fox") ~ "fox", TRUE ~ "h3")) %>%
    mutate(variable = case_when(str_detect(treatment, "dox") ~ "dox", TRUE ~ "ed"))

mod <- aov(data = model, value ~ target * variable)
summary(mod)
TukeyHSD(mod)
plot(mod, which = 2)
plot(mod, which = 1)

#### q values ####
data  %>%
    mutate(target = case_when(str_detect(treatment, "fox") ~ "fox", TRUE ~ "h3")) %>%
    mutate(variable = case_when(str_detect(treatment, "dox") ~ "dox", TRUE ~ "ed"))  %>%
    filter(!grepl("dox_foxa1_low", treatment)) %>%
    ggplot(aes(
        x = factor(variable, levels = c("ed", "dox")),
        y = log2(value),
        color = variable
    )) +
    geom_violin(trim = FALSE) +
    facet_wrap( ~ target) +
    stat_summary(
        fun.data = "mean_sdl",
        fun.args = list(mult = 1),
        geom = "pointrange",
        color = "black"
    ) +
    theme_ipsum() +
    stat_compare_means(
        comparisons = list(c("ed", "dox")),
        method = "wilcox.test",
        label = "p.signif",
        paired = T,
        #change to F to compare means
        label.x = 1.4,
        label.y = 5,
        hide.ns = T,
        bracket.size = 0
    ) +
    labs(subtitle = "Wilcoxon") +
    theme(legend.position = "none")



data  %>%
    mutate(target = case_when(str_detect(treatment, "fox") ~ "fox", TRUE ~ "h3")) %>%
    mutate(variable = case_when(str_detect(treatment, "dox") ~ "dox", TRUE ~ "ed"))  %>%
    filter(!grepl("dox_foxa1_low|h3k27ac", treatment)) %>%
    group_by(treatment) %>%
    summarise_at("value", median)



#### q values for down and up genes ####

data  %>%
    mutate(target = case_when(str_detect(treatment, "fox") ~ "fox", TRUE ~ "h3")) %>%
    mutate(variable = case_when(str_detect(treatment, "dox") ~ "dox", TRUE ~ "ed"))  %>%
    filter(grepl("h3k27ac", treatment)) %>%
    filter(!grepl("NS", diffexpressed)) %>%
    ggplot(aes(
        color = diffexpressed,
        x = variable,
        y = log2(value)
    )) +
    geom_violin(trim = FALSE) +
    facet_wrap(vars(diffexpressed)) +
    stat_summary(
        fun.data = "mean_sdl",
        fun.args = list(mult = 1),
        geom = "pointrange",
        color = "black"
    ) +
    theme_ipsum() +
    stat_compare_means(
        comparisons = list(c("ed", "dox")),
        method = "wilcox.test",
        label = "p.signif",
        paired = F,
        label.x = 1.5,
        hide.ns = T,
        bracket.size = 0
    ) +
    labs(subtitle = "Wilcoxon") +
    theme(legend.position = "none")




### save ####
save.image(file = 'environments/boxplot.RData')
