load('environments/data.RData')

#### bound genes ####
binding <- data %>%
    pivot_longer(cols = contains("bound"),
                 names_to = "target",
                 values_to = "bound") %>%
    mutate(condition = case_when(grepl("ed", target) ~ "ed",
                                 grepl("dox", target) ~ "dox")) %>%
    select(!"treatment") %>%
    unique() %>%
    mutate(binding = case_when(
        grepl("fox", target) ~ "fox",
        grepl("h3k27ac", target) ~ "h3k27ac",
        grepl("overlaid", target) ~ "co-bound"
    )) %>%
    group_by(target, bound) %>%
    mutate(bound_count = n()) %>%
    ungroup() %>%
    mutate(bound_count = bound_count / (n_distinct(condition) + n_distinct(bound))) %>%
    mutate(gene_count = n_distinct(name)) %>%
    filter(bound == "1")




binding %>% ggplot(aes(
    x = factor(
        condition,
        levels = c("ed", "dox"),
        labels = c("-Dox", "+Dox")
    ),
    y = bound_count / gene_count * 100,
    fill = condition
)) +
    geom_bar(position = "dodge", stat = "identity") +
    xlab("") +
    ylab("% of Bound High Confidence hg19 tDNAs") +
    scale_fill_discrete(name = "") +
    theme_classic(base_size = 15) +
    theme(
        legend.position = "none",
        strip.placement = "outside",
        strip.background = element_rect(color = NA),
        panel.spacing = unit(0, "lines")
    ) +
    scale_y_continuous(
        limits = c(0, 100),
        expand = c(0, 0),
        breaks = seq(0, 100, by = 20)
    ) +
    geom_text(
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
        size = 4,
        colour = "white",
        check_overlap = T
    ) +
    facet_wrap(vars(factor(
        binding,
        levels = c("fox", "h3k27ac", "co-bound") ,
        labels = c("FOXA1", "H3K27ac", "Co-bound")
    )),
    strip.position = "bottom")



x <- binding %>%
    filter(binding == "h3k27ac") %>%
    filter(bound == "1") %>%
    filter(condition == "ed") %>%
    select(name, h3_status) %>% unique()

y <- binding %>%
    filter(binding == "fox") %>%
    filter(bound == "1") %>%
    filter(condition == "ed") %>%
    select(name, fox_status) %>% unique()

z <- inner_join(x, y)


x <- binding %>%
    filter(binding == "h3k27ac") %>%
    filter(bound == "1") %>%
    filter(condition == "dox") %>%
    select(name, h3_status) %>% unique()

y <- binding %>%
    filter(binding == "fox") %>%
    filter(bound == "1") %>%
    filter(condition == "dox") %>%
    select(name, fox_status) %>% unique()

a <- inner_join(x, y)

b <- anti_join(a, z)

#### active genes ####

data %>%
    filter(grepl("h3k27ac", treatment)) %>%
    mutate(activity = ifelse(value > activity_threshold, '1', '0')) %>%
    group_by(condition, activity) %>%
    mutate(active_count = n()) %>%
    ungroup() %>%
    mutate(gene_count = n_distinct(name)) %>%
    ggplot(aes(
        x = factor(
            condition,
            levels = c("ed", "dox"),
            labels = c("-Dox", "+Dox")
        ),
        y = active_count / gene_count * 100,
        fill = activity
    )) +
    geom_bar(position = "dodge", stat = "identity") +
    xlab("") +
    ylab("% of High Confidence hg19 tDNAs") +
    scale_fill_discrete(name = "",
                        labels = c('Inactive', 'Active')) +
    theme_classic(base_size = 15) +
    theme(legend.title = element_blank(),
          legend.position = "top") +
    scale_y_continuous(
        limits = c(0, 60),
        expand = c(0, 0),
        breaks = seq(0, 60, by = 20)
    ) +
    geom_text(
        position = position_dodge(width = 1),
        aes(
            label = paste(
                active_count,
                "\n",
                "(",
                round(active_count / gene_count * 100, digits = 1),
                "%)",
                sep = ""
            ),
            y = active_count / gene_count * 100 / 2
        ),
        size = 4,
        colour = "white",
        check_overlap = T
    )


data %>%
    filter(grepl("h3k27ac", treatment)) %>%
    mutate(activity = ifelse(value > activity_threshold, '1', '0')) %>%
    group_by(condition, h3_status) %>%
    mutate(active_count = n()) %>%
    ungroup() %>%
    mutate(gene_count = n_distinct(name)) %>%
    ggplot(aes(
        x = factor(
            h3_status,
            levels = c("loss", "shared", "gain"),
            labels = c("Loss", "Shared", "Gain")
        ),
        y = active_count,
        fill = h3_status
    )) +
    geom_bar(position = "dodge", stat = "identity") +
    xlab("") +
    scale_fill_discrete(name = "") +
    theme_classic(base_size = 15) +
    theme(legend.position = "none") +
    scale_y_continuous(
        limits = c(0, 400),
        expand = c(0, 0),
        breaks = seq(0, 400, by = 100)
    ) +
    geom_text(
        position = position_dodge(width = 1),
        aes(
            label = paste(
                active_count,
                " ",
                "(",
                round(active_count / gene_count * 100, digits = 1),
                "%)",
                sep = ""
            ),
            y = active_count / 2
        ),
        size = 4,
        colour = "white",
        check_overlap = T
    )



#### save ####
save.image(file = 'environments/barchart.RData')


#### bound genes by position ####
binding <- total %>%
    select(!"treatment") %>%
    mutate(binding = case_when(
        grepl("fox", target) ~ "fox",
        grepl("h3k27ac", target) ~ "h3k27ac",
        grepl("overlaid", target) ~ "co-bound"
    )) %>%
    group_by(condition,target, bound, position) %>%
    mutate(bound_count = n()) %>%
    ungroup() %>%
    mutate(bound_count = bound_count / n_distinct(condition)) %>%
    mutate(gene_count = n_distinct(name)) %>%
    filter(bound == "1")




binding %>% ggplot(aes(
    x = factor(
        condition,
        levels = c("ed", "dox"),
        labels = c("-Dox", "+Dox")
    ),
    y = bound_count / gene_count * 100,
    fill = position
)) +
    geom_bar(position = "dodge", stat = "identity") +
    xlab("") +
    ylab("% of Bound High Confidence hg19 tDNAs") +
    scale_fill_discrete(name = "") +
    theme_classic(base_size = 15) +
    theme(
        legend.position = "none",
        strip.placement = "outside",
        strip.background = element_rect(color = NA),
        panel.spacing = unit(0, "lines")
    ) +
    scale_y_continuous(
        limits = c(0, 100),
        expand = c(0, 0),
        breaks = seq(0, 100, by = 20)
    ) +
    geom_text(
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
        size = 4,
        colour = "white",
        check_overlap = T
    ) +
    facet_wrap(vars(factor(
        binding,
        levels = c("fox", "h3k27ac", "co-bound") ,
        labels = c("FOXA1", "H3K27ac", "Co-bound")
    )),
    strip.position = "bottom")
