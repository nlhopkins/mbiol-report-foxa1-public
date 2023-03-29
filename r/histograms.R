load('environments/data.RData')


# foxa1
data %>%
    filter(!grepl("h3k27ac|dox_foxa1_low", treatment)) %>%
    ggplot(aes(x = value, fill = treatment)) +
    geom_histogram(
        color = "#e9ecef",
        alpha = 0.6,
        position = 'identity',
        bins = 50
    ) +
    scale_fill_manual(values = c("#69b3a2", "#404080")) +
    theme_ipsum() +
    labs(fill = "")


# h3k27ac
data  %>%
    filter(!grepl("foxa1", treatment)) %>%
    ggplot(aes(x = value, fill = treatment)) +
    geom_histogram(
        color = "#e9ecef",
        alpha = 0.6,
        position = 'identity',
        bins = 50
    ) +
    scale_fill_manual(values = c("#69b3a2", "#404080")) +
    theme_ipsum() +
    labs(fill = "")


# position
# foxa1
total  %>%
    ggplot(aes(x = value, fill = position)) +
    geom_histogram(
        color = "#e9ecef",
        alpha = 0.6,
        position = 'identity',
        bins = 50
    ) +
    scale_fill_manual(values = c("#69b3a2", "#404080")) +
    theme_ipsum() +
    labs(fill = "")

# h3
total %>%
    filter(!grepl("foxa1", treatment)) %>%
    ggplot(aes(x = value, fill = treatment)) +
    geom_histogram(
        color = "#e9ecef",
        alpha = 0.6,
        position = 'identity',
        bins = 50
    ) +
    scale_fill_manual(values = c("#69b3a2", "#404080")) +
    theme_ipsum() +
    labs(fill = "")

# fold
total  %>%
    filter(!grepl("h3k27ac|dox_foxa1_low", treatment)) %>%
    ggplot(aes(x = fox_fold, fill = position)) +
    geom_histogram(
        color = "#e9ecef",
        alpha = 0.6,
        position = 'identity',
        bins = 50
    ) +
    scale_fill_manual(values = c("#69b3a2", "#404080")) +
    theme_ipsum() +
    labs(fill = "")

total  %>%
    filter(!grepl("foxa1", treatment)) %>%
    ggplot(aes(x = h3_fold, fill = position)) +
    geom_histogram(
        color = "#e9ecef",
        alpha = 0.6,
        position = 'identity',
        bins = 50
    ) +
    scale_fill_manual(values = c("#69b3a2", "#404080")) +
    theme_ipsum() +
    labs(fill = "")


save.image(file = 'environments/histograms.RData')
