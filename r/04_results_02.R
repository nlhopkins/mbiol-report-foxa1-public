#### load environments ####
load(file = 'environments/03_results_01.RData')

#### FOXA1 volcano ####
# do chisq test on FOXA1 -dox vs +dox
volcano_fox <- data %>%
    pivot_wider(
        # pivot wider for readability
        names_from = "treatment",
        values_from = "value",
        id_cols = c("name")
    ) %>%
    select(c(name, ed_foxa1, dox_foxa1_high)) %>% # select name and colums to compare
    mutate(p_value = -log10(pmap_dbl(select(., -1),  ~ chisq.test(c(
        # get p value from chisq
        ...
    ))$p.value))) %>%
    mutate(fox_fold = log2(dox_foxa1_high / ed_foxa1)) %>% # calculate fold change
    filter(fox_fold != "Inf") %>% # remove infinite values (from dividing by 0)
    distinct() %>%  # remove duplicates
    mutate(diffexpressed = ifelse(
        fox_fold > log2(1.5) & p_value > .001,
        "UP",
        # if log2 > 1.5 and p > 0.001 = gain fox enrichment
        ifelse(fox_fold < -log2(1.5) & p_value > .001,
               "DN", # if log2 < 1.5 and p > 0.001 = gain fox enrichment
               "NS") # otherwise not significant
    ))



# plot
volcano_fox %>%
    add_count(diffexpressed) %>% # add counts
    mutate(diffexpressed = factor(diffexpressed,
                                  levels = c("DN", "NS", "UP"))) %>% # levels = order displayed in plots
    mutate(diffexpressed = paste0(diffexpressed, ' (', n, ')')) %>% # add counts to expression status for plot
    ggplot(aes(x = fox_fold, # scatter plot
               y = p_value,
               colour = diffexpressed)) +
    geom_point(size = 3) +
    geom_vline(xintercept = c(-log2(1.5), log2(1.5)), col = "grey") + # threshold lines
    geom_hline(yintercept = (.01), col = "grey") +
    scale_x_continuous(limits = c(-4, 4),
                       breaks = seq(-4, 4, by = 2)) + # adjust scales
    scale_y_continuous(
        limits = c(-0.2, 2),
        expand = c(0, 0),
        breaks = seq(0, 2, by = 0.5)
    ) +
    theme_classic(base_size = 40) + # theme
    theme(
        legend.position = "top",
        strip.placement = "outside",
        strip.background = element_rect(color = NA),
        panel.spacing = unit(0, "lines")
    ) +
    scale_color_manual(name = "", # pretty colours :)
                       # colours
                       values =
                           c("#3ec1b6", "grey", "#aca4e0")) +
    xlab("log2(FOXA1 Fold Change)") + # x title
    ylab("-log10(P-value)") + # y title
    guides(colour = guide_legend(override.aes = list(size = # legend dot size
                                                         10)))








#### H3k27ac volcano ####
volcano_h3 <- data %>%
    pivot_wider(
        # pivot wider for readability
        names_from = "treatment",
        values_from = "value",
        id_cols = c("name")
    ) %>%
    select(c(name, ed_h3k27ac, dox_h3k27ac)) %>% # select name and colums to compare
    mutate(p_value = -log10(pmap_dbl(select(., -1),  ~ chisq.test(c(
        # get p value from chisq
        ...
    ))$p.value))) %>%
    mutate(h3_fold = log2(dox_h3k27ac / ed_h3k27ac)) %>% # calculate fold change
    filter(h3_fold != "Inf") %>% # remove infinite values (from dividing by 0)
    distinct() %>%  # remove duplicates
    mutate(diffexpressed = ifelse(
        h3_fold > log2(1.5) & p_value > .001,
        "UP",
        # if log2 > 1.5 and p > 0.001 = gain fox enrichment
        ifelse(h3_fold < -log2(1.5) & p_value > .001,
               "DN", # if log2 < 1.5 and p > 0.001 = gain fox enrichment
               "NS") # otherwise not significant
    ))


# plot
volcano_h3 %>%
    add_count(diffexpressed) %>% # add counts
    mutate(diffexpressed = factor(diffexpressed,
                                  levels = c("DN", "NS", "UP"))) %>% # levels = order displayed in plots
    mutate(diffexpressed = paste0(diffexpressed, ' (', n, ')')) %>%  # add counts to expression status for plot
    ggplot(aes(x = h3_fold, # scatter plot
               y = p_value,
               colour = diffexpressed)) +
    geom_point(size = 3) +
    geom_vline(xintercept = c(-log2(1.5), log2(1.5)), col = "grey") + # threshold lines
    geom_hline(yintercept = (.01), col = "grey") +
    scale_x_continuous(limits = c(-4, 4),
                       breaks = seq(-4, 4, by = 2)) + # adjust scales
    scale_y_continuous(
        limits = c(-0.2, 2),
        expand = c(0, 0),
        breaks = seq(0, 2, by = 0.5)
    ) +
    theme_classic(base_size = 40) + # theme
    theme(
        legend.position = "top",
        strip.placement = "outside",
        strip.background = element_rect(color = NA),
        panel.spacing = unit(0, "lines")
    ) +
    scale_color_manual(name = "",
                       # colours
                       values =
                           c("#adb364", "grey", "#dc9c87")) +
    xlab("log2(FOXA1 Fold Change)") + # x title
    ylab("-log10(P-value)") + # y title
    guides(colour = guide_legend(override.aes = list(size = # legend dot size
                                                         10)))


#### Table of both fox and h3k27ac change ####

diffexpressed <-
    full_join(volcano_fox, volcano_h3, by = "name", copy = T) %>% # merge fox_diff and h3_diff
    pivot_longer(cols = contains(c("diffexpressed.x", "diffexpressed.y")),
                 names_to = "comparison",
                 values_to = "diffexpressed") %>% # pivot catagories to combine them
    mutate(comparison = str_replace(comparison, "diffexpressed.x", "fox_diff")) %>% # replace x with foxa1
    mutate(comparison = str_replace(comparison, "diffexpressed.y", "h3_diff")) %>% # replace y with h3k27ac
    distinct() %>%
    merge(data, by = "name") %>% # merge with data to get other catagories
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
    filter(
        grepl("fox", comparison) &
            grepl("fox", treatment) |
            grepl("h3", comparison) &
            grepl("h3", treatment)
    ) # make sure diff matched treament

diffexpressed %>% ggplot(aes(
    # plot
    fill = factor(
        condition,
        levels = c("ed", "dox"),
        labels = c("-Dox", "+Dox")
    ),
    x = diffexpressed,
    y = log2(value)
)) +
    geom_violin(trim = FALSE, # violin
                color = NA) +
    facet_wrap(vars(factor(
        # facet by diff catagories
        comparison,
        levels = c("fox_diff", "h3_diff"),
        labels = c("FOXA1", "H3K27ac")
    )),
    nrow = 1,
    strip.position = "bottom") +
    stat_summary(
        # box plot
        fun.data = "median_iqr",
        geom = "pointrange",
        color = "black",
        position = position_dodge(0.9)
    ) +
    stat_compare_means(
        # compare means
        aes(group = condition),
        method = "wilcox.test",
        label = "p.signif",
        # change to "p" for no
        paired = F,
        label.x = 1.5,
        label.y = 6,
        hide.ns = T,
        size = 14,
        na.rm = T
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
    labs(subtitle = "Mann-Whitney U") +
    scale_y_continuous(
        limits = c(-6, 9),
        expand = c(0, 0),
        breaks = seq(-6, 9, by = 3)
    ) +
    scale_fill_manual(name = "",
                      values = c("#ACA4E1",
                                 "#39BDB1"))


#### Venn ####

fox_diff <-
    diffexpressed %>% filter(diffexpressed == "UP" &
                                 treatment == "dox_foxa1_high") %>% select(name) # tDNAs with UP FOXA1

h3_diff <-
    diffexpressed %>% filter(diffexpressed == "UP" &
                                 treatment == "dox_h3k27ac") %>% select(name) # tDNAs with UP H3K27ac

both_up <- inner_join(fox_diff, h3_diff) # tDNAs with UP both



eulerr_options(labels = list(cex = 3)) # venn text size

euler(c(
    # venn
    Fox_Gain = nrow(fox_diff),
    H3K27ac_Gain = nrow(h3_diff),
    "Fox_Gain&H3K27ac_Gain" = nrow(both_up)
)) %>% plot(
    quantities = list(type = c("counts", "percent"), cex = 3),
    lty = 0,
    alpha = 0.75,
    fills = c(
        Fox_Gain = "#ACA4E1",
        H3K27ac_Gain = "#DB9D85"
    ),
    labels = c(Fox_Gain = "FOXA1 UP",
               H3K27ac_Gain = "H3K27ac UP")
)

save.image(file = "environments/04_results_02.RData")
