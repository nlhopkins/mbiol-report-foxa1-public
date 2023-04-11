load('environments/data.RData')


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
    scale_colour_manual(values = qualitative_hcl(2, palette = "Cold")) +
    geom_violin(trim = FALSE) +
    facet_wrap(vars(factor(
        target,
        levels = c("fox", "h3") ,
        labels = c("FOXA1", "H3K27ac")
    )),
    strip.position = "bottom") +
    stat_summary(
        fun.data = "mean_sdl",
        fun.args = list(mult = 1),
        geom = "pointrange",
        color = "black"
    ) +
    theme_classic(base_size = 20) +
    xlab("") +
    ylab("log2(Q-Value)") +
    stat_compare_means(
        comparisons = list(c("ed", "dox")),
        method = "wilcox.test",
        label = "p.signif",
        paired = T,
        #change to F to compare means
        label.x = 1.4,
        label.y = 5,
        hide.ns = T,
        bracket.size = 0,
        size = 7
    ) +
    labs(subtitle = "Wilcoxon") +
    theme(
        legend.position = "none",
        strip.placement = "outside",
        strip.background = element_rect(color = NA),
        panel.spacing = unit(0, "lines")
    ) +
    scale_y_continuous(
        limits = c(-6, 9),
        expand = c(0, 0),
        breaks = seq(-6, 9, by = 3)
    ) +
    scale_x_discrete(labels = c("ed" = "-Dox", "dox" = "+Dox"))


#### medians ####
data  %>%
    mutate(target = case_when(str_detect(treatment, "fox") ~ "fox", TRUE ~ "h3")) %>%
    mutate(variable = case_when(str_detect(treatment, "dox") ~ "dox", TRUE ~ "ed"))  %>%
    filter(!grepl("dox_foxa1_low|h3k27ac", treatment)) %>%
    group_by(treatment) %>%
    summarise_at("value", median)



#### q values for down and up genes ####

data %>%
    merge(diffexpressed, by = "name") %>%
    mutate(target = case_when(str_detect(treatment, "fox") ~ "fox", TRUE ~ "h3")) %>%
    mutate(variable = case_when(str_detect(treatment, "dox") ~ "dox", TRUE ~ "ed"))  %>%
    filter(grepl("h3k27ac", treatment)) %>%
    filter(!grepl("NS", diffexpressed)) %>%
    ggplot(aes(
        color = variable,
        x = factor(variable, levels = c("ed", "dox")),
        y = log2(value)
    )) +
    scale_colour_manual(values = qualitative_hcl(2, palette = "Cold")) +
    geom_violin(trim = FALSE) +
    facet_wrap(vars(diffexpressed), strip.position = "bottom") +
    stat_summary(
        fun.data = "mean_sdl",
        fun.args = list(mult = 1),
        geom = "pointrange",
        color = "black"
    ) +
    xlab("") +
    ylab("Q-Value") +
    theme_classic(base_size = 20) +
    theme(
        legend.position = "none",
        strip.placement = "outside",
        strip.background = element_rect(color = NA),
        panel.spacing = unit(0, "lines")
    ) +
    stat_compare_means(
        comparisons = list(c("ed", "dox")),
        method = "wilcox.test",
        label = "p.signif",
        paired = F,
        label.x = 1.5,
        hide.ns = T,
        bracket.size = 0,
        size = 7
    ) +
    labs(subtitle = "Mann-Whitney U") +
    theme(legend.position = "none") +
    scale_y_continuous(
        limits = c(-6, 9),
        expand = c(0, 0),
        breaks = seq(-6, 9, by = 3)
    ) +
    scale_x_discrete(labels = c("ed" = "-Dox", "dox" = "+Dox"))




### save ####
save.image(file = 'environments/boxplot.RData')
