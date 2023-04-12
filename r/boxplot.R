load('environments/data.RData')


#### q values ####
data  %>%
    mutate(target = case_when(str_detect(treatment, "fox") ~ "fox", TRUE ~ "h3")) %>%
    mutate(variable = case_when(str_detect(treatment, "dox") ~ "dox", TRUE ~ "ed"))  %>%
    filter(!grepl("dox_foxa1_low", treatment)) %>%
    ggplot(aes(
        x = factor(variable, levels = c("ed", "dox")),
        y = log2(value),
        fill = treatment
    )) +
    scale_fill_manual(
        values = c(
            dox_foxa1_high = "#aca4e0",
            ed_foxa1 = "#3ec1b6",
            dox_h3k27ac = "#dc9c87",
            ed_h3k27ac = "#adb364"
        )
    ) +
    geom_violin(trim = FALSE,
                color = NA) +
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
    mutate(target = case_when(str_detect(comparison, "fox") ~ "FOXA1", TRUE ~ "H3K27ac")) %>%
    mutate(variable = case_when(str_detect(treatment, "dox") ~ "dox", TRUE ~ "ed")) %>%
    mutate(target = target, factor(target, levels = c("FOXA1", "H3K27ac"))) %>%
    filter(!grepl("NS", diffexpressed)) %>%
    unite('tag', c("comparison", "diffexpressed"), remove = F) %>%
    unite('colour', c("tag", "variable"), remove = F) %>%
    ggplot(aes(
        group = interaction(variable, condition),
        fill = colour,
        x = factor(variable, levels = c("ed", "dox")),
        y = log2(value)
    )) +
    geom_violin(trim = FALSE,
                color = NA)  +
    scale_fill_manual(
        values = c(
            "#aca4e0",
            "#3ec1b6",
            "#dc9c87",
            "#adb364",
            "#aca4e0",
            "#3ec1b6",
            "#dc9c87",
            "#adb364"
        )
    ) +
    facet_wrap(target ~ diffexpressed,
               nrow = 1,
               strip.position = "bottom") +
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
        aes(group = variable),
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
