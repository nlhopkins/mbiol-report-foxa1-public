#### Load Data ####
raw <-
    read.delim("data/raw/hg19.txt") %>% # Read in quantified data from EaSeq
    janitor::clean_names()


#### Activity Threshold ####
# get median -Dox H3K27ac for 'active' threshold

activity_threshold <- raw %>%
    mutate(across(contains("ed"), ~ .x / ed_input)) %>% # normalise -Dox to -Dox input
    mutate(across(contains("dox"), ~ .x / dox_input)) %>% # normalise Dox to Dox input
    pivot_longer(cols = contains(c("_foxa1", "_h3k27ac")), # pivot for R readability
                 names_to = "treatment") %>%
    group_by(treatment) %>% # group by treatment
    summarise_at("value", median) %>% # get medians
    filter(treatment == "ed_h3k27ac") %>% # filter by -Dox
    pull(2) # extract value


data <- raw %>%
    mutate(across(contains("ed"), ~ .x / ed_input)) %>% # normalise -Dox to -Dox input
    mutate(across(contains("dox"), ~ .x / dox_input)) %>% # normalise Dox to Dox input
    mutate(ed_fox_bound = ifelse(ed_foxa1 > 1, '1', '0')) %>% # 'bound' if fold-change > 1, 0 = not bound, 1 = bound
    mutate(dox_h3k27ac_bound = ifelse(dox_h3k27ac > 1, '1', '0')) %>%
    mutate(ed_h3k27ac_bound = ifelse(ed_h3k27ac > 1, '1', '0')) %>%
    mutate(dox_fox_bound = ifelse(dox_foxa1_high > 1, '1', '0')) %>%
    mutate(ed_cobound_bound = ifelse(ed_foxa1 > 1 &
                                         ed_h3k27ac > 1, '1', '0')) %>% # cobound if both foxa1 + h3k27ac bound
    mutate(dox_cobound_bound = ifelse(dox_foxa1_high > 1 &
                                          dox_h3k27ac > 1, '1', '0')) %>%
    mutate(fox_fold = dox_foxa1_high / ed_foxa1) %>% # FOXA1 -dox/+dox fold change
    mutate(h3_fold = dox_h3k27ac / ed_h3k27ac)  %>% # FOXA1 -dox/+dox fold change
    mutate(
        activity_status = ifelse(
            ed_h3k27ac < activity_threshold & dox_h3k27ac > activity_threshold,
            "GAIN",
            ifelse(
                ed_h3k27ac > activity_threshold &
                    dox_h3k27ac < activity_threshold,
                "LOSS",
                "NC"
            )
        )
    ) %>%
    pivot_longer(cols = contains(c("_foxa1", "_h3k27ac")) &
                     !contains("bound"),
                 names_to = "treatment") %>%
    filter(treatment != "dox_foxa1_low") %>%
    mutate(condition = str_extract(treatment, "[^_]+")) %>%
    mutate(amino = str_replace(name, "^(([^-]+-){1}[^-]+)-.*", "\\1")) %>%
    mutate(activity = ifelse(value > activity_threshold, '1', '0'))



#### relative position ####
upstream <- read.delim("data/raw/upstream.txt") %>%
    janitor::clean_names() %>%
    mutate(position = "upstream")

downstream <- read.delim("data/raw/downstream.txt") %>%
    janitor::clean_names() %>%
    mutate(position = "downstream")

total_activity_threshold <- full_join(upstream,
                                      downstream) %>%
    mutate(across(contains("ed"), ~ .x / ed_input)) %>%
    mutate(across(contains("dox"), ~ .x / dox_input)) %>%
    pivot_longer(cols = contains(c("_foxa1", "_h3k27ac")),
                 names_to = "treatment") %>%
    group_by(treatment) %>%
    summarise_at("value", median, na.rm = TRUE) %>%
    filter(treatment == "ed_h3k27ac") %>%
    pull(2)

total <- full_join(upstream,
                   downstream) %>%
    mutate(across(contains("ed"), ~ .x / ed_input)) %>% # normalise -Dox to -Dox input
    mutate(across(contains("dox"), ~ .x / dox_input)) %>% # normalise Dox to Dox input
    mutate(ed_fox_bound = ifelse(ed_foxa1 > 1, '1', '0')) %>% # 'bound' if fold-change > 1, 0 = not bound, 1 = bound
    mutate(dox_h3k27ac_bound = ifelse(dox_h3k27ac > 1, '1', '0')) %>%
    mutate(ed_h3k27ac_bound = ifelse(ed_h3k27ac > 1, '1', '0')) %>%
    mutate(dox_fox_bound = ifelse(dox_foxa1_high > 1, '1', '0')) %>%
    mutate(fox_fold = dox_foxa1_high / ed_foxa1) %>% # FOXA1 -dox/+dox fold change
    mutate(h3_fold = dox_h3k27ac / ed_h3k27ac)  %>% # FOXA1 -dox/+dox fold change
    mutate(
        activity_status = ifelse(
            ed_h3k27ac < activity_threshold & dox_h3k27ac > activity_threshold,
            "GAIN",
            ifelse(
                ed_h3k27ac > activity_threshold &
                    dox_h3k27ac < activity_threshold,
                "LOSS",
                "NC"
            )
        )
    ) %>%
    pivot_longer(cols = contains(c("_foxa1", "_h3k27ac")) &
                     !contains("bound"),
                 names_to = "treatment") %>%
    filter(treatment != "dox_foxa1_low") %>%
    mutate(condition = str_extract(treatment, "[^_]+")) %>%
    mutate(amino = str_replace(name, "^(([^-]+-){1}[^-]+)-.*", "\\1")) %>%
    mutate(activity = ifelse(value > activity_threshold, '1', '0'))



save.image(file = 'environments/02_data.RData')
