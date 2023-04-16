load('environments/data.RData')


#### q values ####
data  %>%
    pivot_longer(cols = contains("bound"),
                 names_to = "target",
                 values_to = "bound") %>%
    mutate(target = case_when(str_detect(treatment, "fox") ~ "fox", TRUE ~ "h3")) %>%
    select(c("name", "condition", "value", "treatment", "target")) %>%
    distinct() %>%
    ggplot(aes(
        x = factor(condition, levels = c("ed", "dox")),
        y = log2(value),
        fill = treatment
    )) +
    scale_fill_manual(
        values = c(
            dox_foxa1_high = "#ACA4E1",
            ed_foxa1 = "#39BDB1",
            dox_h3k27ac = "#DB9D85",
            ed_h3k27ac = "#ABB064"
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
    stat_summary(fun.data = "median_iqr",
                 geom = "pointrange",
                 color = "black") +
    theme_classic(base_size = 40) +
    xlab("") +
    ylab("log2(Q-value)") +
    stat_compare_means(
        comparisons = list(c("ed", "dox")),
        method = "wilcox.test",
        label = "p",
        paired = F,
        #change to F to compare means
        label.x = 1.4,
        label.y = 5,
        hide.ns = T,
        bracket.size = 0,
        size = 14
    ) +
    labs(subtitle = "Mann-Whitney U") +
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
    pivot_longer(cols = contains("bound"),
                 names_to = "target",
                 values_to = "bound") %>%
    mutate(target = case_when(str_detect(treatment, "fox") ~ "fox", TRUE ~ "h3")) %>%
    select(c("name", "condition", "value", "treatment", "target")) %>%
    distinct() %>%
    group_by(condition, target) %>%
    summarise_at("value", median)



#### q values for down and up genes ####

data %>%
    merge(diffexpressed, by = "name") %>%
    filter(diffexpressed != "NS") %>%
    select(c(
        "name",
        "value",
        "diffexpressed",
        "treatment",
        "comparison",
        "condition"
    )) %>%
    distinct() %>%
    drop_na() %>%
    ggplot(aes(
        fill = factor(
            condition,
            levels = c("ed", "dox"),
            labels = c("-Dox", "+Dox")
        ),
        x = diffexpressed,
        y = log2(value)
    )) +
    geom_violin(trim = FALSE,
                color = NA)  +
    scale_fill_manual(name = "",
                      values = c("#ACA4E1",
                                 "#39BDB1")) +
    facet_wrap(vars(factor(
        comparison,
        levels = c("fox_diff", "h3_diff"),
        labels = c("FOXA1", "H3K27ac")
    )),
    nrow = 1,
    strip.position = "bottom") +
    stat_summary(
        fun.data = "median_iqr",
        geom = "pointrange",
        color = "black",
        position = position_dodge(0.9)
    ) +
    xlab("") +
    ylab("log2(Q-value") +
    theme_classic(base_size = 40) +
    theme(
        legend.position = "top",
        strip.placement = "outside",
        strip.background = element_rect(color = NA),
        panel.spacing = unit(0, "lines")
    ) +
    stat_compare_means(
        aes(group = condition),
        method = "wilcox.test",
        label = "p.signif",
        paired = F,
        label.x = 1.5,
        label.y = 6,
        hide.ns = T,
        size = 14,
        na.rm = T
    ) +
    labs(subtitle = "Mann-Whitney U") +
    scale_y_continuous(
        limits = c(-6, 9),
        expand = c(0, 0),
        breaks = seq(-6, 9, by = 3)
    )




data %>%
    merge(diffexpressed, by = "name") %>%
    filter(diffexpressed != "NS") %>%
    select(c(
        "name",
        "value",
        "diffexpressed",
        "treatment",
        "comparison",
        "condition"
    )) %>%
    distinct() %>%
    drop_na() %>%
    group_by(diffexpressed, comparison, condition) %>%
    summarise_at("value", median)





#### active genes ####
data %>%
    filter(grepl("h3k27ac", treatment)) %>%
    mutate(activity = ifelse(value > activity_threshold, '1', '0')) %>%
    ggplot(aes(
        fill = condition,
        x = factor(condition, levels = c("ed", "dox")),
        y = log2(value)
    )) +
    geom_violin(trim = FALSE,
                color = NA)  +
    scale_fill_manual(values = c("#adb364", "#dc9c87")) +
    facet_wrap(vars(factor(
        activity, labels = c("0" = "Inactive", "1" = "Active")
    )),
    nrow = 1,
    strip.position = "bottom") +
    
    stat_summary(fun.data = "median_iqr",
                 
                 geom = "pointrange",
                 color = "black") +
    xlab("") +
    ylab("log2(Q-value)") +
    theme_classic(base_size = 40) +
    theme(
        legend.position = "none",
        strip.placement = "outside",
        strip.background = element_rect(color = NA),
        panel.spacing = unit(0, "lines")
    ) +
    stat_compare_means(
        aes(group = condition),
        comparisons = list(c("ed", "dox")),
        method = "wilcox.test",
        label = "p.signif",
        paired = F,
        label.x = 1.5,
        hide.ns = T,
        bracket.size = 0,
        size = 14
    ) +
    labs(subtitle = "Mann-Whitney U") +
    theme(legend.position = "none") +
    scale_y_continuous(
        limits = c(-6, 9),
        expand = c(0, 0),
        breaks = seq(-6, 9, by = 3)
    ) +
    scale_x_discrete(labels = c("ed" = "-Dox", "dox" = "+Dox"))

data %>%
    filter(grepl("h3k27ac", treatment)) %>%
    mutate(activity = ifelse(value > activity_threshold, '1', '0')) %>%
    group_by(condition, activity) %>%
    summarise_at("value", median)


### inactive vs active active foxa1
data %>%
    drop_na() %>%
    filter(treatment == c("ed_foxa1", "dox_foxa1_high")) %>%
    ggplot(aes(
        fill = condition,
        x = factor(condition, levels = c("ed", "dox")),
        y = log2(value)
    )) +
    geom_violin(trim = FALSE,
                color = NA)  +
    scale_fill_manual(values = c("#adb364", "#dc9c87")) +
    facet_wrap(vars(factor(
        activity_status,
        levels = c("LOSS", "NC", "GAIN")
    )),
    nrow = 1,
    strip.position = "bottom") +
    
    stat_summary(fun.data = "median_iqr",
                 
                 geom = "pointrange",
                 color = "black") +
    xlab("") +
    ylab("log2(Q-value)") +
    theme_classic(base_size = 40) +
    theme(
        legend.position = "none",
        strip.placement = "outside",
        strip.background = element_rect(color = NA),
        panel.spacing = unit(0, "lines")
    ) +
    stat_compare_means(
        aes(group = condition),
        comparisons = list(c("ed", "dox")),
        method = "wilcox.test",
        label = "p.signif",
        paired = F,
        label.x = 1.5,
        hide.ns = T,
        bracket.size = 0,
        size = 14
    ) +
    labs(subtitle = "Mann-Whitney U") +
    theme(legend.position = "none") +
    scale_y_continuous(
        limits = c(-6, 9),
        expand = c(0, 0),
        breaks = seq(-6, 9, by = 3)
    ) +
    scale_x_discrete(labels = c("ed" = "-Dox", "dox" = "+Dox"))

data %>%
    drop_na() %>%
    filter(treatment == c("ed_foxa1", "dox_foxa1_high")) %>%
    group_by(activity_status, condition) %>%
    summarise_at("value", median)




#### upstream and downstream
total %>%
    mutate(
        target = case_when(
            grepl("ed_fox_bound", target) ~ "ed_foxa1",
            grepl("ed_h3k27ac_bound", target) ~ "ed_h3k27ac",
            grepl("dox_fox_bound", target) ~ "dox_foxa1_high",
            grepl("dox_h3k27ac_bound", target) ~ "dox_h3k27ac",
            grepl("ed_overlaid_bound", target) ~ "ed_co-bound",
            grepl("dox_overlaid_bound", target) ~ "dox_co-bound"
        )
    ) %>%
    drop_na() %>%
    filter(bound == "1") %>%
    mutate(target = case_when(
        grepl("foxa1", target) ~ "foxa1",
        grepl("h3k27ac", target) ~ "h3k27ac",
        grepl("co-bound", target) ~ "co-bound"
    )) %>%
    mutate(treatment = case_when(
        grepl("fox", treatment) ~ "foxa1",
        grepl("h3k27ac", treatment) ~ "h3k27ac"
    )) %>%
    filter(target != "co-bound") %>%
    filter(treatment == target) %>%
    distinct() %>%
    ggplot(aes(
        fill = factor(
            condition,
            levels = c("ed", "dox"),
            labels = c("-Dox", "+Dox")
        ),
        x = factor(
            position,
            levels = c("upstream", "downstream"),
            labels = c("Upstream", "Downstream")
        ),
        y = log2(value)
    )) +
    geom_violin(trim = FALSE,
                color = NA)  +
    scale_fill_manual(name = "",
                      values = c("#ACA4E1", "#39BDB1")) +
    facet_wrap(vars(treatment),
               nrow = 1,
               strip.position = "bottom") +
    stat_summary(
        fun.data = "median_iqr",
        
        geom = "pointrange",
        color = "black",
        position = position_dodge(0.9)
    ) +
    xlab("") +
    ylab("log2(Q-value)") +
    theme_classic(base_size = 40) +
    theme(
        legend.position = "top",
        strip.placement = "outside",
        strip.background = element_rect(color = NA),
        panel.spacing = unit(0, "lines")
    ) +
    stat_compare_means(
        aes(group = condition),
        method = "wilcox.test",
        label = "p.signif",
        paired = F,
        label.x = 1.5,
        label.y = 6,
        hide.ns = T,
        size = 14,
        na.rm = T
    ) +
    labs(subtitle = "Mann-Whitney U") +
    scale_y_continuous(
        limits = c(-6, 9),
        expand = c(0, 0),
        breaks = seq(-6, 9, by = 3)
    )


total %>%
    mutate(
        target = case_when(
            grepl("ed_fox_bound", target) ~ "ed_foxa1",
            grepl("ed_h3k27ac_bound", target) ~ "ed_h3k27ac",
            grepl("dox_fox_bound", target) ~ "dox_foxa1_high",
            grepl("dox_h3k27ac_bound", target) ~ "dox_h3k27ac"
        )
    ) %>%
    filter(target == treatment) %>%
    drop_na() %>%
    filter(bound == "1") %>%
    mutate(target = case_when(
        grepl("ed_foxa1", target) ~ "foxa1",
        grepl("ed_h3k27ac", target) ~ "h3k27ac",
        grepl("dox_foxa1_high", target) ~ "foxa1",
        grepl("dox_h3k27ac", target) ~ "h3k27ac"
    )) %>%
    distinct() %>%
    group_by(target, position, condition) %>%
    summarise_at("value", median)




total %>%
    mutate(
        target = case_when(
            grepl("ed_fox_bound", target) ~ "ed_foxa1",
            grepl("ed_h3k27ac_bound", target) ~ "ed_h3k27ac",
            grepl("dox_fox_bound", target) ~ "dox_foxa1_high",
            grepl("dox_h3k27ac_bound", target) ~ "dox_h3k27ac",
            grepl("ed_overlaid_bound", target) ~ "ed_co-bound",
            grepl("dox_overlaid_bound", target) ~ "dox_co-bound"
        )
    ) %>%
    drop_na() %>%
    filter(bound == "1") %>%
    mutate(target = case_when(
        grepl("foxa1", target) ~ "foxa1",
        grepl("h3k27ac", target) ~ "h3k27ac",
        grepl("co-bound", target) ~ "co-bound"
    )) %>%
    mutate(treatment = case_when(
        grepl("fox", treatment) ~ "foxa1",
        grepl("h3k27ac", treatment) ~ "h3k27ac"
    )) %>%
    filter(target != "co-bound") %>%
    filter(treatment == target) %>%
    distinct() %>%
    ggplot(aes(
        fill = factor(
            position,
            levels = c("upstream", "downstream"),
            labels = c("Upstream", "Downstream")
        ),
        x = factor(
            condition,
            levels = c("ed", "dox"),
            labels = c("-Dox", "+Dox")
        ),
        y = log2(value)
    )) +
    geom_violin(trim = FALSE,
                color = NA)  +
    scale_fill_manual(name = "",
                      values = c("#ACA4E1", "#39BDB1")) +
    facet_wrap(vars(treatment),
               nrow = 1,
               strip.position = "bottom") +
    stat_summary(
        fun.data = "median_iqr",
        
        geom = "pointrange",
        color = "black",
        position = position_dodge(0.9)
    ) +
    xlab("") +
    ylab("log2(Q-value)") +
    theme_classic(base_size = 40) +
    theme(
        legend.position = "top",
        strip.placement = "outside",
        strip.background = element_rect(color = NA),
        panel.spacing = unit(0, "lines")
    ) +
    stat_compare_means(
        aes(group = position),
        method = "wilcox.test",
        label = "p.signif",
        paired = F,
        label.x = 1.5,
        label.y = 6,
        hide.ns = T,
        size = 14,
        na.rm = T
    ) +
    labs(subtitle = "Mann-Whitney U") +
    scale_y_continuous(
        limits = c(-6, 9),
        expand = c(0, 0),
        breaks = seq(-6, 9, by = 3)
    )





### save ####
save.image(file = 'environments/boxplot.RData')