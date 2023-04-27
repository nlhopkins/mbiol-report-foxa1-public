#### load environment ####
load(file = "environments/04_results_02.RData")

#### Violin of active tdnas ####
data %>%
    filter(grepl("h3k27ac", treatment)) %>% # filter only h3k27ac values
    ggplot(aes( # plot
        fill = condition,
        x = factor(condition, levels = c("ed", "dox")),
        y = log2(value)
    )) +
    geom_violin(trim = FALSE, # violin
                color = NA)  +
    scale_fill_manual(values = c("#adb364", "#dc9c87")) +
    facet_wrap(vars(factor( # facet
        activity, labels = c("0" = "Inactive", "1" = "Active")
    )),
    nrow = 1,
    strip.position = "bottom") +
    
    stat_summary(fun.data = "median_iqr", # boxplot
                 
                 geom = "pointrange",
                 color = "black") +
    xlab("") + # x title
    ylab("log2(Q-value)") + # y title
    theme_classic(base_size = 40) + # theme font size
    theme(
        legend.position = "none",
        strip.placement = "outside",
        strip.background = element_rect(color = NA),
        panel.spacing = unit(0, "lines")
    ) +
    stat_compare_means( # compare means
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
    scale_x_discrete(labels = c("ed" = "-Dox", "dox" = "+Dox")) # x axis labels

#### barplot of counts of active tdnas ####
data %>%
    filter(grepl("h3k27ac", treatment)) %>% # filter to only h3k27ac
    group_by(condition, activity) %>% # group by to get counts based on these categories
    mutate(active_count = n()) %>% # count number tdnas
    ungroup() %>% # remove grouping
    mutate(gene_count = n_distinct(name)) %>% # count total number of tdnas in each catagory
    ggplot(aes( # plot
        x = factor(
            condition,
            levels = c("ed", "dox"),
            labels = c("-Dox", "+Dox")
        ),
        y = active_count / gene_count * 100, # y as percentage
        fill = activity
    )) +
    geom_bar(position = "dodge", stat = "identity") + # barplot
    xlab("") +
    ylab("% of tDNAs") +
    scale_fill_manual( # colour of bars
        name = "",
        labels = c('Inactive', 'Active'),
        values = c("#DB9D85", "#86B875")
    ) +
    theme_classic(base_size = 40) + # theme + font size
    theme(legend.title = element_blank(), 
          legend.position = "top") +
    scale_y_continuous( 
        limits = c(0, 60),
        expand = c(0, 0),
        breaks = seq(0, 60, by = 20)
    ) +
    geom_text( # add count + percentage to bars
        position = position_dodge(width = 1),
        aes(
            label = paste(
                active_count, #counts
                "\n",
                "(",
                round(active_count / gene_count * 100, digits = 1), # percentage
                "%)",
                sep = ""
            ),
            y = active_count / gene_count * 100 / 2 # label height = half of %
        ),
        size = 10,
        colour = "white",
        check_overlap = T
    )
#### barplot of counts of gain/loss tdnas ####
data %>%
filter(grepl("h3k27ac", treatment)) %>% # filter to only h3k27ac
    group_by(condition, activity_status) %>% # group by these categories
    mutate(active_count = n()) %>% # get counts
    ungroup() %>%
    mutate(gene_count = n_distinct(name)) %>% # total number of tdnas in each group
    ggplot(aes( # plot
        x = factor(activity_status,
                   levels = c("LOSS", "NC", "GAIN")),
        y = active_count,
        fill = activity_status
    )) +
    geom_bar(position = "dodge", stat = "identity") + # barplot
    xlab("") +
    ylab("No. Active tDNAs") +
    scale_fill_manual(name = "", # pretty colours :)
                      values = c("#86B875", "#E092C3", "grey")) +
    theme_classic(base_size = 40) +
    theme(legend.position = "none") +
    scale_y_continuous(
        limits = c(0, 400),
        expand = c(0, 0),
        breaks = seq(0, 400, by = 100)
    ) +
    geom_text( # add counts and percentages to bars
        position = position_dodge(width = 1), # position of text
        aes(
            label = paste(
                active_count, # count
                " ",
                "(",
                round(active_count / gene_count * 100, digits = 1), # pecrcent
                "%)",
                sep = ""
            ),
            y = active_count / 2 # text height = half of %
        ),
        size = 10,
        colour = "white",
        check_overlap = T
    )

#### Violin of FOXA1 by gain/loss ####
data %>%
    filter(grepl("foxa1", treatment)) %>% # 
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
        paired = F, # change to T for paired wilcoxon text (compares same tdnas rather than means)
        label.x = 1.5,
        hide.ns = T,
        bracket.size = 0,
        size = 14
    ) +
    labs(subtitle = "Mann-Whitney U") + # stat test name
    theme(legend.position = "none") +
    scale_y_continuous(
        limits = c(-6, 9),
        expand = c(0, 0),
        breaks = seq(-6, 9, by = 3)
    ) +
    scale_x_discrete(labels = c("ed" = "-Dox", "dox" = "+Dox")) # x axis names

#### Venn ####
# other groups taken from 04_results_02.R

activity_gain <- data %>% filter(activity_status == "GAIN") %>% select(name) %>% distinct() # tDNAs with UP FOXA1

gain_fox_up <- inner_join(activity_gain, fox_diff) # activated + gained fox

gain_h3_up <- inner_join(activity_gain, h3_diff) # activated + gained h3k27ac

gain_up_up <- inner_join(activity_gain, both_up) # activited + gained fox + gained h3k27ac


euler(
    c(
        GAIN = nrow(activity_gain),
        FOXUP = nrow(fox_diff),
        H3UP = nrow(h3_diff),
        "GAIN&FOXUP" = nrow(gain_fox_up),
        "GAIN&H3UP" = nrow(gain_h3_up),
        "GAIN&FOXUP&H3UP" = nrow(gain_up_up),
        "FOXUP&H3UP" = nrow(both_up)
    ),
    shape = "circle"
)  %>%
    plot(
        quantities = list(type = c("counts"), cex = 3),
        lty = 0,
        alpha = 0.75,
        fills = c(
            GAIN = "#86B875",
            FOXUP = "#ACA4E1",
            H3UP = "#DB9D85"
        ),
        labels = c(GAIN = "Activity GAIN",
                   FOXUP = "FOXA1 UP",
                   H3UP = "H3K27ac UP")
    )


save.image(file = "environments/05_results_03.RData")
