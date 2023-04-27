#### table of clusters/families ####
tdnas <- data %>% select(c("name", "fox_fold", "h3_fold","treatment", "value")) %>% 
    pivot_wider(
    names_from = "treatment",
    values_from = "value") %>%
    mutate(ed_activity = ifelse(ed_h3k27ac > activity_threshold,
                                "1",
                                "0")) %>%
    
    mutate(dox_activity = ifelse(dox_h3k27ac > activity_threshold,
                                 "1",
                                 "0"))

imet <- tdnas %>% filter(
    name %in% c(
        "tRNA-iMet-CAT-1-1",
        "tRNA-iMet-CAT-1-2",
        "tRNA-iMet-CAT-1-3",
        "tRNA-iMet-CAT-1-4",
        "tRNA-iMet-CAT-1-5",
        "tRNA-iMet-CAT-1-6",
        "tRNA-iMet-CAT-1-7",
        "tRNA-iMet-CAT-1-8",
        "tRNA-iMet-CAT-2-1"
    )
) %>% mutate(group = "iMET",
             func = "Transcription Initiation")


met <- tdnas %>% filter(
    name %in% c(
        "tRNA-Met-CAT-1-1",
        "tRNA-Met-CAT-2-1",
        "tRNA-Met-CAT-3-1",
        "tRNA-Met-CAT-3-2",
        "tRNA-Met-CAT-4-1",
        "tRNA-Met-CAT-4-2",
        "tRNA-Met-CAT-4-3",
        "tRNA-Met-CAT-5-1",
        "tRNA-Met-CAT-6-1",
        "tRNA-Met-CAT-7-1"
    )
) %>% mutate(group = "MET",
             func = "iMET Control")



argccg <- tdnas %>% filter(
    name %in% c(
        "tRNA-Arg-CCG-1-1",
        "tRNA-Arg-CCG-1-2",
        "tRNA-Arg-CCG-1-3",
        "tRNA-Arg-CCG-2-1"
    )
) %>% mutate(group = "Arg-CCG",
             func = "Promotes Breast Cancer Metastasis")

gluttc <- tdnas %>% filter(
    name %in% c(
        "tRNA-Glu-TTC-1-1",
        "tRNA-Glu-TTC-1-2",
        "tRNA-Glu-TTC-2-1",
        "tRNA-Glu-TTC-2-2",
        "tRNA-Glu-TTC-3-1",
        "tRNA-Glu-TTC-4-1",
        "tRNA-Glu-TTC-4-2"
    )
) %>% mutate(group = "Glu-TTC",
             func = "Promotes Breast Cancer Metastasis")

sec <- tdnas %>% filter(name %in% c("tRNA-SeC-TCA-1-1",
                                "tRNA-SeC-TCA-2-1")) %>%
    mutate(group = "SeC",
           func = "Oxidation (REDOX) & Poor Survival")

ebersole <- tdnas %>% filter(
    name %in% c(
        "tRNA-Glu-CTC-1-5",
        "tRNA-Gly-TCC-2-5",
        "tRNA-Asp-GTC-2-5",
        "tRNA-Leu-CAG-1-5"
    )
) %>% mutate(group = "Ebersole",
             func = "Chromatin Insulator")



aloxe <- tdnas %>% filter(
    name %in% c(
        "tRNA-Lys-TTT-3-5",
        "tRNA-Gln-CTG-1-5",
        "tRNA-Leu-TAG-1-1",
        "tRNA-Arg-TCT-2-1"
    )
) %>% mutate(group = "ALOXE3",
             func = "Chromatin Insulator")

hes <- tdnas %>% filter(name %in% c("tRNA-Gly-GCC-2-6")) %>%
    mutate(group = "HES7",
           func = "Proximate to ALOXE")

per1 <- tdnas %>% filter(name %in% c("tRNA-Ser-CGA-1-1",
                                 "tRNA-Thr-AGT-5-1")) %>%
    mutate(group = "Per1",
           func = "Chromatin Insulator")

tmem <- tdnas %>% filter(
    name %in% c(
        "tRNA-Trp-CCA-3-3",
        "tRNA-Ser-GCT-4-3",
        "tRNA-Thr-AGT-1-1",
        "tRNA-Ile-AAT-5-5"
    )
) %>%
    mutate(group = "TMEM107",
           func = "Chromatin Insulator")

clusters <- list(aloxe,
           ebersole,
           hes,
           per1,
           tmem,
           argccg,
           gluttc,
           imet,
           met,
           sec) %>%
    reduce(full_join) %>%
    mutate(ed_activity = ifelse(ed_activity == "1", "Active", "Inactive")) %>%
    mutate(dox_activity = ifelse(dox_activity == "1", "Active", "Inactive")) %>%
    mutate_if(is.numeric, ~ round(., 2)) %>%
    rename(
        "-Dox Activity" = ed_activity,
        "+Dox Activity" = dox_activity,
        "FOXA1 Fold Change" = fox_fold,
        "H3K27ac Fold Change" = h3_fold
    ) %>% select(c("group", "name","func", "+Dox Activity", "-Dox Activity", "FOXA1 Fold Change","H3K27ac Fold Change"))

