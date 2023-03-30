load('environments/data.RData')

# whisker plots
x <- data %>%
    mutate(target = case_when(str_detect(treatment, "fox") ~ "fox", TRUE ~ "h3")) %>%
    mutate(variable = case_when(str_detect(treatment, "dox") ~ "dox", TRUE ~ "ed"))

mod <- aov(data = x, value ~ target * variable)
summary(mod)
TukeyHSD(mod)
plot(mod, which = 2)
plot(mod, which = 1)

data  %>%
    mutate(target = case_when(str_detect(treatment, "fox") ~ "fox", TRUE ~ "h3")) %>%
    mutate(variable = case_when(str_detect(treatment, "dox") ~ "dox", TRUE ~ "ed"))  %>%
    filter(!grepl("dox_foxa1_low", treatment)) %>%
    ggplot(aes(
        x = target,
        y = log2(value),
        fill = variable
    )) +
    geom_boxplot() +
    theme_ipsum() +
    stat_compare_means(method = "wilcox.test",
                       label = "p",
                       paired = T)




data  %>%
    mutate(target = case_when(str_detect(treatment, "fox") ~ "fox", TRUE ~ "h3")) %>%
    mutate(variable = case_when(str_detect(treatment, "dox") ~ "dox", TRUE ~ "ed"))  %>%
    filter(!grepl("dox_foxa1_low", treatment)) %>%
    ggplot(aes(
        x = target,
        y = log2(value),
        fill = factor(variable, levels = c("ed", "dox"))
    )) +
    geom_boxplot() +
    theme_ipsum() +
    stat_compare_means(method = "wilcox.test",
                       label = "p.signif",
                       paired = T) +
    labs(subtitle = "Wilcoxon")


# relative position
total  %>%
    filter(!grepl("dox_foxa1_low", treatment)) %>%
    ggplot(aes(x = treatment, y = log10(value), fill = position)) +
    geom_boxplot() +
    theme_ipsum() +
    facet_wrap("position")


save.image(file = 'environments/boxplot.RData')

