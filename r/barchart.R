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
    fill = target
)) +
    geom_bar(position = "dodge", stat = "identity") +
    xlab("") +
    ylab("% of Bound tDNAs") +
    theme_classic(base_size = 40) +
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
        size = 10,
        colour = "white",
        check_overlap = T
    ) +
    facet_wrap(vars(factor(
        binding,
        levels = c("fox", "h3k27ac", "co-bound") ,
        labels = c("FOXA1", "H3K27ac", "Co-bound")
    )),
    strip.position = "bottom") +
    scale_fill_manual(
        name = "",
        values = c(
            dox_fox_bound = "#ACA4E1",
            ed_fox_bound = "#39BDB1",
            dox_h3k27ac_bound = "#DB9D85",
            ed_h3k27ac_bound = "#ABB064",
            ed_overlaid_bound = "#E092C3",
            dox_overlaid_bound = "#86B875"
        )
    )



binding %>%
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
    ylab("% of tDNAs") +
    scale_fill_manual(
        name = "",
        labels = c('Inactive', 'Active'),
        values = c("#DB9D85", "#86B875")
    ) +
    theme_classic(base_size = 40) +
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
        size = 10,
        colour = "white",
        check_overlap = T
    )


data %>%
    filter(grepl("h3k27ac", treatment)) %>%
    group_by(condition, activity_status) %>%
    mutate(active_count = n()) %>%
    ungroup() %>%
    mutate(gene_count = n_distinct(name)) %>%
    ggplot(aes(
        x = factor(activity_status,
                   levels = c("LOSS", "NC", "GAIN")),
        y = active_count,
        fill = activity_status
    )) +
    geom_bar(position = "dodge", stat = "identity") +
    xlab("") +
    ylab("No. Active tDNAs") +
    scale_fill_manual(name = "",
                      values = c("#86B875", "#E092C3", "grey")) +
    theme_classic(base_size = 40) +
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
        size = 10,
        colour = "white",
        check_overlap = T
    )



#### bound genes by position ####
binding <- total %>%
    mutate(condition = case_when(grepl("ed", target) ~ "ed",
                                 grepl("dox", target) ~ "dox")) %>%
    unique() %>%
    filter(bound == "1") %>%
    select(!"bound") %>%
    distinct() %>%
    mutate(binding = case_when(
        grepl("fox", target) ~ "fox",
        grepl("h3k27ac", target) ~ "h3k27ac",
        grepl("overlaid", target) ~ "co-bound"
    )) %>%
    group_by(position, treatment, target, condition, binding) %>%
    mutate(bound_count = n()) %>%
    ungroup() %>%
    group_by(binding) %>%
    mutate(gene_count = n_distinct(name)) %>%
    mutate(binding = factor(binding,
                            levels = c("fox",
                                       "h3k27ac",
                                       "co-bound"))) %>%
    filter(
        grepl("ed_fox", treatment) & grepl("ed_fox", target) |
            grepl("ed_h3k27ac", treatment) &
            grepl("ed_h3k27ac", target) |
            grepl("dox_fox", treatment) & grepl("dox_fox", target) |
            grepl("dox_h3k27ac", treatment) &
            grepl("dox_h3k27ac", target)
    ) %>% mutate(star = case_when(grepl("h3k", target) ~ "***"))


binding %>% ggplot(aes(
    x = factor(
        position,
        levels = c("upstream", "downstream"),
        labels = c("Upstream", "Downstream")
    ),
    y = bound_count / gene_count * 100,
    fill = factor(
        condition,
        levels = c("ed", "dox"),
        labels = c("-Dox", "+Dox")
    )
)) +
    geom_bar(position = "dodge", stat = "identity") +
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
        size = 7,
        colour = "white",
        check_overlap = T
    ) +
    facet_wrap(vars(factor(
        binding,
        levels = c("fox", "h3k27ac", "co-bound"),
        labels = c("FOXA1", "H3K27ac", "Co-bound")
    )),
    nrow = 1,
    strip.position = "bottom") +
    geom_text(aes(label = star, y = 90),
              colour = "black",
              size = 10)

binding %>% ggplot(aes(
    x = factor(
        condition,
        levels = c("ed", "dox"),
        labels = c("-Dox", "+Dox")
    ),
    y = bound_count / gene_count * 100,
    fill = factor(
        position,
        levels = c("upstream", "downstream"),
        labels = c("Upstream", "Downstream")
    )
)) +
    geom_bar(position = "dodge", stat = "identity") +
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
        size = 6,
        colour = "white",
        check_overlap = T
    ) +
    facet_wrap(vars(factor(
        binding,
        levels = c("fox", "h3k27ac", "co-bound"),
        labels = c("FOXA1", "H3K27ac", "Co-bound")
    )),
    nrow = 1,
    strip.position = "bottom") +
    geom_text(aes(label = "****", y = 90),
              colour = "black",
              size = 10)






x <- upstream  %>%
    mutate(across(contains("ed"), ~ .x / ed_input)) %>%
    mutate(across(contains("dox"), ~ .x / dox_input)) %>%
    filter_all(all_vars(!is.infinite(.))) %>%
    mutate(ed_overlaid = ifelse(ed_foxa1 > 1 &
                                    ed_h3k27ac > 1, '1', '0')) %>%
    mutate(dox_overlaid = ifelse(dox_foxa1_high > 1 &
                                     dox_h3k27ac > 1, '1', '0')) %>%
    mutate(ed_foxa1 = ifelse(ed_foxa1 > 1, '1', '0')) %>%
    mutate(dox_h3k27ac = ifelse(dox_h3k27ac > 1, '1', '0')) %>%
    mutate(ed_h3k27ac = ifelse(ed_h3k27ac > 1, '1', '0')) %>%
    mutate(dox_foxa1_high = ifelse(dox_foxa1_high > 1, '1', '0')) %>%
    select(
        c(
            "ed_foxa1",
            "dox_h3k27ac",
            "ed_h3k27ac",
            "dox_foxa1_high",
            "ed_overlaid",
            "dox_overlaid",
        )
    ) %>% tibble::rownames_to_column()





y <- downstream %>%
    mutate(across(contains("ed"), ~ .x / ed_input)) %>%
    mutate(across(contains("dox"), ~ .x / dox_input)) %>%
    filter_all(all_vars(!is.infinite(.))) %>%
    mutate(ed_overlaid = ifelse(ed_foxa1 > 1 &
                                    ed_h3k27ac > 1, '1', '0')) %>%
    mutate(dox_overlaid = ifelse(dox_foxa1_high > 1 &
                                     dox_h3k27ac > 1, '1', '0')) %>%
    mutate(ed_foxa1 = ifelse(ed_foxa1 > 1, '1', '0')) %>%
    mutate(dox_h3k27ac = ifelse(dox_h3k27ac > 1, '1', '0')) %>%
    mutate(ed_h3k27ac = ifelse(ed_h3k27ac > 1, '1', '0')) %>%
    mutate(dox_foxa1_high = ifelse(dox_foxa1_high > 1, '1', '0')) %>%
    select(
        c(
            "ed_foxa1",
            "dox_h3k27ac",
            "ed_h3k27ac",
            "dox_foxa1_high",
            "ed_overlaid",
            "dox_overlaid"
        )
    ) %>% tibble::rownames_to_column()

z <- full_join(x,
               y,
               keep = T,
               suffix = c(".x", ".y"),
               by = "rowname") %>% select(!c("rowname.x", "rowname.y"))


library(rstatix)



chisq_test(z$ed_h3k27ac.x, z$ed_h3k27ac.y)
chisq_test(z$dox_h3k27ac.x, z$dox_h3k27ac.y)
chisq_test(z$ed_foxa1.x, z$ed_foxa1.y)
chisq_test(z$dox_foxa1_high.x, z$dox_foxa1_high.y)

chisq_test(z$ed_h3k27ac.x, z$dox_h3k27ac.x)
chisq_test(z$ed_foxa1.x, z$dox_foxa1_high.x)
chisq_test(z$ed_h3k27ac.x, z$dox_h3k27ac.x)
chisq_test(z$ed_foxa1.x, z$dox_foxa1_high.x)

chisq_test(z$ed_h3k27ac.y, z$dox_h3k27ac.y)
chisq_test(z$ed_foxa1.y, z$dox_foxa1_high.y)
chisq_test(z$ed_h3k27ac.y, z$dox_h3k27ac.y)
chisq_test(z$ed_foxa1.y, z$dox_foxa1_high.y)
