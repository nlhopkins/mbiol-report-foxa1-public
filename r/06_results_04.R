#### load environment ####
load(file = 'environments/02_data.RData')


#### total positional data ####
total_binding <- total %>%
    pivot_longer(cols = contains("bound"),
                 # pivot for readability
                 names_to = "target",
                 values_to = "bound") %>%
    mutate(condition = case_when(grepl("ed", target) ~ "ed", # extract condition
                                 grepl("dox", target) ~ "dox")) %>% # remove treatment column
    mutate(binding = case_when(
        # extract categories
        grepl("fox", target) ~ "fox",
        grepl("h3k27ac", target) ~ "h3k27ac"
    )) %>%
    select(c("name", "condition", "bound", "binding", "target", "position")) %>% # select columns we need
    distinct() %>%
    group_by(target, position) %>% # grouping for counts
    filter(bound == "1") %>%
    mutate(bound_count = n()) %>% # get count of bound tdnas
    ungroup() %>%
    mutate(gene_count = n_distinct(name)) # get number of tdnas (for percentages)


#### barplot -dox vs +dox ####
total_binding  %>% ggplot(aes(
    # plot
    x = factor(
        condition,
        levels = c("ed", "dox"),
        labels = c("-Dox", "+Dox")
    ),
    y = bound_count / gene_count * 100,
    # percentage
    fill = factor( # fill = position
        position,
        levels = c("upstream", "downstream"),
        labels = c("Upstream", "Downstream")
    )
)) +
    geom_bar(position = "dodge", stat = "identity") + # barplot
    facet_wrap(vars(factor( # facet by binding target
        binding,
        levels = c("fox", "h3k27ac"),
        labels = c("FOXA1", "H3K27ac")
    )),
    nrow = 1,
    strip.position = "bottom") +
    xlab("") +
    ylab("% of Bound tDNAs") +
    scale_fill_manual(name = "",
                      values = c("#ACA4E1", "#39BDB1")) +
    theme_classic(base_size = 40) +
    theme(
        legend.position = "top",
        strip.placement = "outside",
        strip.background = element_rect(color = NA),
        panel.spacing = unit(0, "lines")
    ) +
    scale_y_continuous(
        limits = c(0, 100),
        expand = c(0, 0),
        breaks = seq(0, 100, by = 20)
    ) +
    geom_text( # text showing the counts and percentage
        position = position_dodge(width = 1),
        aes(
            label = paste(
                bound_count,
                "\n",
                "(",
                round(bound_count / gene_count * 100, digits = 1),
                "%)",
                sep = ""
            ),
            y = (bound_count / gene_count * 100 / 2)
        ),
        size = 7,
        colour = "white",
        check_overlap = T
    )


#### barplot up vs down ####
total_binding  %>% ggplot(aes( #plot
    x = factor(
        position,
        levels = c("upstream", "downstream"),
        labels = c("Upstream", "Downstream")
    ),
    y = bound_count / gene_count * 100, # count and percentage
    fill = factor(
        condition,
        levels = c("ed", "dox"),
        labels = c("-Dox", "+Dox")
    )
)) +
    geom_bar(position = "dodge", stat = "identity") + # barplot
    facet_wrap(vars(factor( # facet by binding target
        binding,
        levels = c("fox", "h3k27ac"),
        labels = c("FOXA1", "H3K27ac")
    )),
    nrow = 1,
    strip.position = "bottom") +
    xlab("") + # x title
    ylab("% of Bound tDNAs") + # y title
    scale_fill_manual(name = "", # colours
                      values = c("#ACA4E1", "#39BDB1")) +
    theme_classic(base_size = 40) + # theme + font size
    theme(
        legend.position = "top",
        strip.placement = "outside",
        strip.background = element_rect(color = NA),
        panel.spacing = unit(0, "lines")
    ) +
    scale_y_continuous(
        limits = c(0, 100),
        expand = c(0, 0),
        breaks = seq(0, 100, by = 20)
    ) +
    geom_text( # add counts and pecerntage to the bars
        position = position_dodge(width = 1),
        aes(
            label = paste(
                bound_count,
                "\n",
                "(",
                round(bound_count / gene_count * 100, digits = 1),
                "%)",
                sep = ""
            ),
            y = (bound_count / gene_count * 100 / 2)
        ),
        size = 7,
        colour = "white",
        check_overlap = T
    )



#### boxplot of bound tdnas by position ####
total %>%
    pivot_longer(cols = contains("bound"), # pivot bound catargories into 2 columns
                 # pivot for readability
                 names_to = "target",
                 values_to = "bound") %>%
    mutate(target = case_when(str_detect(treatment, "fox") ~ "fox", TRUE ~ "h3")) %>% # extract target
    select(c("name", "condition", "value", "treatment", "target", "position")) %>% # select these columns
    distinct() %>% # remove repeated rows
    ggplot(aes( # plot
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
    facet_wrap(vars(factor( # facet by binding target
        target,
        levels = c("fox", "h3"),
        labels = c("FOXA1", "H3K27ac")
    )),
    nrow = 1,
    strip.position = "bottom") +
    scale_fill_manual(name = "", # colours
                      values = c("#ACA4E1", "#39BDB1")) +
    stat_summary( # boxplot
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
    stat_compare_means( # compare means
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
