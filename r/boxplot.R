load('environments/data.RData')

# whisker plots
model <- data %>%
    mutate(target = case_when(str_detect(treatment, "fox") ~ "fox", TRUE ~ "h3")) %>%
    mutate(variable = case_when(str_detect(treatment, "dox") ~ "dox", TRUE ~ "ed"))

mod <- aov(data = model, value ~ target * variable)
summary(mod)
TukeyHSD(mod)
plot(mod, which = 2)
plot(mod, which = 1)

## q values
data  %>%
    mutate(target = case_when(str_detect(treatment, "fox") ~ "fox", TRUE ~ "h3")) %>%
    mutate(variable = case_when(str_detect(treatment, "dox") ~ "dox", TRUE ~ "ed"))  %>%
    filter(!grepl("dox_foxa1_low", treatment)) %>%
    ggplot(aes(
        color = variable,
        x = factor(variable, levels = c("ed", "dox")),
        y = log2(value)
    )) +
    geom_violin(trim = FALSE) +
    stat_summary(
        fun.data = "mean_sdl",
        fun.args = list(mult = 1),
        geom = "pointrange",
        color = "black"
    ) +
    theme_ipsum() +
    stat_compare_means(
        method = "wilcox.test",
        label = "p.signif",
        paired = T,
        label.x = 1.4,
        label.y = 5
    ) +
    labs(subtitle = "Wilcoxon") +
    facet_wrap(vars(target)) +
    theme(legend.position = "none")


## q values for down and up genes

data  %>%
    mutate(target = case_when(str_detect(treatment, "fox") ~ "fox", TRUE ~ "h3")) %>%
    mutate(variable = case_when(str_detect(treatment, "dox") ~ "dox", TRUE ~ "ed"))  %>%
    filter(grepl("h3k27ac", treatment)) %>%
    filter(!grepl("NS", diffexpressed)) %>%
    ggplot(aes(
        color = diffexpressed,
        x = diffexpressed,
        y = log2(value)
    )) +
    geom_violin(trim = FALSE) +
    stat_summary(
        fun.data = "mean_sdl",
        fun.args = list(mult = 1),
        geom = "pointrange",
        color = "black"
    ) +
    theme_ipsum() +
    stat_compare_means(
        method = "wilcox.test",
        label = "p.signif",
        paired = F,
        label.x = 1.5
    ) +
    labs(subtitle = "Wilcoxon") +
    facet_wrap(vars(variable)) + theme(legend.position = "none")




# relative position
total  %>%
    filter(!grepl("dox_foxa1_low", treatment)) %>%
    ggplot(aes(x = condition, 
               y = log2(value), fill = position)) +
    geom_boxplot() +
    theme_ipsum() 


save.image(file = 'environments/boxplot.RData')
