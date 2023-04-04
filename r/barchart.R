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
        grepl("h3k27ac", target) ~ "h3k27ac"
    )) %>%
    group_by(target, bound) %>%
    mutate(bound_count = n()) %>%
    ungroup() %>%
    mutate(bound_count = bound_count / n_distinct(target)) %>%
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
    ylab("% of High Confidence hg19 tDNAs") +
    scale_fill_discrete(name = "") +
    theme_ipsum(base_size = 15,
                axis_title_size = 15,
                strip_text_size = 15) +
    theme(legend.position = "none") +
    scale_y_continuous() +
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
    facet_wrap(~ factor(binding, labels = c("FOXA1", "H3K27ac")))



x <- binding %>%
    filter(binding == "h3k27ac") %>%
    filter(bound == "1") %>%
    filter(condition == "ed") %>%
    select(name) %>% unique()

y <- binding %>%
    filter(binding == "fox") %>%
    filter(bound == "1") %>%
    filter(condition == "ed") %>%
    select(name) %>% unique()

z <- inner_join(x, y)


x <- binding %>%
    filter(binding == "h3k27ac") %>%
    filter(bound == "1") %>%
    filter(condition == "dox") %>%
    select(name) %>% unique()

y <- binding %>%
    filter(binding == "fox") %>%
    filter(bound == "1") %>%
    filter(condition == "dox") %>%
    select(name) %>% unique()

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
    theme_ipsum(base_size = 15,
                axis_title_size = 15) +
    theme(legend.title = element_blank()) +
    scale_y_continuous() +
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
    theme_ipsum(base_size = 15,
                axis_title_size = 15) +
    theme(legend.position = "none") +
    scale_y_continuous() +
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
