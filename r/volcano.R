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
    mutate(diffexpressed = paste0(diffexpressed, ' (', n, ')')) %>%
    ggplot(aes(x = fox_fold,
               y = p_value,
               colour = diffexpressed)) +
    geom_point(size = 1) +
    theme_minimal() +
    geom_vline(xintercept = c(-0.25, 0.25), col = "grey") +
    geom_hline(yintercept = (.01), col = "grey") +
    scale_color_manual(values =
                           c("#7B52AE", "grey", "#94C773")) +
    scale_x_continuous(limits = c(-5, 11)) +
    theme(legend.position = "none") +
    ggrepel::geom_text_repel(aes(label = label),
                             size = 3,
                             show.legend = FALSE) +
    theme_ipsum() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())



top_h3 = volcano_h3 %>%
    filter(diffexpressed == "UP") %>%
    arrange(., desc(p_value)) %>%
    head(., 10)

# Add column label, containing the gene name for the top hits or nothing for all others
volcano_h3$label = if_else(volcano_h3$name %in% top_h3$name,
                           volcano_h3$name, NA)




volcano_h3 %>%
    add_count(diffexpressed) %>%
    mutate(diffexpressed = paste0(diffexpressed, ' (', n, ')')) %>%
    ggplot(aes(x = h3_fold,
               y = p_value,
               col = diffexpressed)) +
    geom_point(size = 0.5) +
    theme_minimal() +
    geom_vline(xintercept = c(-0.25, 0.25), col = "grey") +
    geom_hline(yintercept = (.01), col = "grey") + scale_color_manual(values =
                                                                                 c("#7B52AE", "grey", "#94C773")) +
    scale_x_continuous(limits = c(-5, 22)) +
    ggrepel::geom_text_repel(aes(label = label),
                             size = 3, show.legend = FALSE) +
    theme_ipsum() + theme(panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank())


top_fox = volcano_fox %>%
    filter(diffexpressed == "UP") %>%
    arrange(., desc(p_value)) %>%
    head(., 10)


save.image(file = 'environments/volcano.RData')

