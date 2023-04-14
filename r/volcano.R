load('environments/data.RData')

options(scipen = 999)
options(ggrepel.max.overlaps = Inf)

top_fox = volcano_fox %>%
    filter(diffexpressed != "NS") %>%
    arrange(., desc(p_value)) %>%
    head(., 10)

# Add column label, containing the gene name for the top hits or nothing for all others
volcano_fox$label = if_else(volcano_fox$name %in% top_fox$name,
                            volcano_fox$name, NA)

volcano_fox %>%
    add_count(diffexpressed) %>%
    mutate(diffexpressed = factor(diffexpressed,
                                  levels = c("DN", "NS", "UP"))) %>%
    mutate(diffexpressed = paste0(diffexpressed, ' (', n, ')')) %>%
    ggplot(aes(x = fox_fold,
               y = p_value,
               colour = diffexpressed)) +
    geom_point(size = 3) +
    geom_vline(xintercept = c(-log2(1.5), log2(1.5)), col = "grey") +
    geom_hline(yintercept = (.01), col = "grey") +
    scale_color_manual(name = "",
                       values =
                           c("#3ec1b6", "grey", "#aca4e0")) +
    scale_x_continuous(limits = c(-4, 4),
                       breaks = seq(-4, 4, by = 2)) +
    scale_y_continuous(
        limits = c(-0.2, 2),
        expand = c(0, 0),
        breaks = seq(0, 2, by = 0.5)
    ) +
    theme_classic(base_size = 40) +
    theme(
        legend.position = "top",
        strip.placement = "outside",
        strip.background = element_rect(color = NA),
        panel.spacing = unit(0, "lines")
    ) +
    xlab("log2(FOXA1 Fold Change)") +
    ylab("-log10(P-value)") + guides(colour = guide_legend(override.aes = list(size =
                                                                                   10)))


top_h3 = volcano_h3 %>%
    filter(diffexpressed != "NS") %>%
    arrange(., desc(p_value)) %>%
    head(., 10)


# Add column label, containing the gene name for the top hits or nothing for all others
volcano_h3$label = if_else(volcano_h3$name %in% top_h3$name,
                           volcano_h3$name, NA)




volcano_h3 %>%
    add_count(diffexpressed) %>%
    mutate(diffexpressed = factor(diffexpressed,
                                  levels = c("DN", "NS", "UP"))) %>% 
    mutate(diffexpressed = paste0(diffexpressed, ' (', n, ')')) %>%
    ggplot(aes(x = h3_fold,
               y = p_value,
               colour = diffexpressed)) +
    geom_point(size = 3) +
    geom_vline(xintercept = c(-log2(1.5), log2(1.5)), col = "grey") +
    geom_hline(yintercept = (.01), col = "grey") +
    scale_color_manual(name = "",
                       values =
                           c("#adb364", "grey", "#dc9c87")) +
    scale_x_continuous(limits = c(-4, 4)) +
    scale_y_continuous(
        limits = c(-0.2, 2),
        expand = c(0, 0),
        breaks = seq(0, 2, by = 0.5)
    ) +
    theme_classic(base_size = 40) +
    theme(
        legend.position = "top",
        strip.placement = "outside",
        strip.background = element_rect(color = NA),
        panel.spacing = unit(0, "lines")
    ) +
    xlab("log2(H3K27ac Fold Change)") +
    ylab("-log10(P-value)") + guides(colour = guide_legend(override.aes = list(size =
                                                                                   10)))


top_fox = volcano_fox %>%
    filter(diffexpressed == "UP") %>%
    arrange(., desc(p_value)) %>%
    head(., 10)


save.image(file = 'environments/volcano.RData')

x <- full_join(volcano_fox, volcano_h3, by = "name")

x %>% filter(diffexpressed.x == "UP" & diffexpressed.y == "UP") %>%
    ggplot(aes(x = dox_foxa1_high, y = dox_h3k27ac)) +
    geom_point()

x %>% filter(diffexpressed.x == "UP" & diffexpressed.y == "UP") %>%
    ggplot(aes(x = fox_fold, y = h3_fold)) +
    geom_point()

cor.test(x$fox_fold, x$h3_fold, method = "pearson")

cor.test(x$dox_foxa1_high, x$dox_h3k27ac, method = "pearson")
