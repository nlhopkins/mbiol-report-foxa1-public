#### load data file ####
load('environments/02_data.RData')


#### Violin plot of Q-values of bound tDNAs (Results 1C) ####
data  %>%
    pivot_longer(cols = contains("bound"),
                 # pivot for readability
                 names_to = "target",
                 values_to = "bound") %>%
    mutate(target = case_when(str_detect(treatment, "fox") ~ "fox", TRUE ~ "h3")) %>% # extract target
    select(c("name", "condition", "value", "treatment", "target")) %>% # select these columns
    distinct() %>% # remove repeated rows
    ggplot(aes(
        # create plot
        x = factor(condition, levels = c("ed", "dox")),
        y = log2(value),
        fill = treatment
    )) +
    geom_violin(trim = FALSE, # violin plot
                color = NA) +
    facet_wrap(vars(factor(
        target,
        levels = c("fox", "h3") ,
        labels = c("FOXA1", "H3K27ac")
    )),
    strip.position = "bottom") +
    stat_summary(fun.data = "median_iqr",
                 # add boxplot
                 geom = "pointrange",
                 color = "black") +
    theme_classic(base_size = 40) +
    stat_compare_means(
        # compare means
        comparisons = list(c("ed", "dox")),
        method = "wilcox.test",
        label = "p",
        # change to "p.signif" to show asterisks
        paired = F,
        label.x = 1.4,
        label.y = 5,
        hide.ns = T,
        bracket.size = 0,
        size = 14
    ) +
    scale_y_continuous( # adjust scale
        limits = c(-6, 9),
        expand = c(0, 0),
        breaks = seq(-6, 9, by = 3)
    ) + # make pretty (can ignore everything below if needed)
    xlab("") + # remove x label
    ylab("log2(Q-value)") + # rename y
    labs(subtitle = "Mann-Whitney U") + # subtitle of test used
    theme(
        legend.position = "none",
        strip.placement = "outside",
        strip.background = element_rect(color = NA),
        panel.spacing = unit(0, "lines")
    ) +
    scale_x_discrete(labels = c("ed" = "-Dox", "dox" = "+Dox")) + # rename x scale
    scale_fill_manual( # nice colours :)
        values = c(
            dox_foxa1_high = "#ACA4E1",
            ed_foxa1 = "#39BDB1",
            dox_h3k27ac = "#DB9D85",
            ed_h3k27ac = "#ABB064"
        )
    )



#### Counts of bound tDNAs ####
# add a column with the sum of bound tdnas
binding <- data %>%
    pivot_longer(cols = contains("bound"), # pivot for readability
                 names_to = "target",
                 values_to = "bound") %>%
    mutate(condition = case_when(grepl("ed", target) ~ "ed", # extract condition
                                 grepl("dox", target) ~ "dox")) %>% # remove treatment column
    mutate(binding = case_when( # extract categories
        grepl("fox", target) ~ "fox",
        grepl("h3k27ac", target) ~ "h3k27ac",
        grepl("cobound", target) ~ "co-bound"
    )) %>%
    select(c("name", "condition", "bound", "binding", "target")) %>% # select columns we need 
    distinct() %>% 
    group_by(bound, condition, binding) %>% # grouping for counts
    filter(bound == "1") %>% 
    mutate(bound_count = n()) %>% # get count of bound tdnas
    ungroup() %>% 
    mutate(gene_count = n_distinct(name)) # get number of tdnas (for percentages)




binding %>% ggplot(aes(
    x = factor( # barplot of bound tdnas counts
        condition,
        levels = c("ed", "dox"),
        labels = c("-Dox", "+Dox")
    ),
    y = bound_count / gene_count * 100,
    fill = target
)) +
    geom_bar(position = "dodge", stat = "identity") + # barplot
    facet_wrap(vars(factor( # facet by catagories
        binding,
        levels = c("fox", "h3k27ac", "co-bound") ,
        labels = c("FOXA1", "H3K27ac", "Co-bound")
    )),
    strip.position = "bottom") +
    xlab("") + # remove x title
    ylab("% of Bound tDNAs") + # rename y title
    theme_classic(base_size = 40) + # plot theme, base size is text size
    theme(
        legend.position = "none",
        strip.placement = "outside",
        strip.background = element_rect(color = NA),
        panel.spacing = unit(0, "lines")
    ) +
    scale_y_continuous( # adjust y axis scale
        limits = c(0, 100),
        expand = c(0, 0),
        breaks = seq(0, 100, by = 20)
    ) +
    geom_text( # add counts and percentages to the bars
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
            y = bound_count / gene_count * 100 / 2
        ),
        size = 10,
        colour = "white",
        check_overlap = T
    ) +
    scale_fill_manual( # pretty colours :)
        name = "",
        values = c(
            dox_fox_bound = "#ACA4E1",
            ed_fox_bound = "#39BDB1",
            dox_h3k27ac_bound = "#DB9D85",
            ed_h3k27ac_bound = "#ABB064",
            ed_cobound_bound = "#E092C3",
            dox_cobound_bound = "#86B875"
        )
    )


save.image(file = 'environments/03_results_01.RData')
