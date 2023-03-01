library(tidyverse)
library(ggplot2)
library(ggvenn)

#### high_confidence ####
high_confidence <-
    read_delim(
        "data/raw/high_confidence_tRNA_genes.txt",
        delim = "\t",
        col_names = TRUE,
        col_types = "ccc",
    ) %>%
    janitor::clean_names()

#### easeq ####
easeq <- read_delim(
    "data/raw/easeq_quantified.txt",
    delim = "\t",
    col_names = TRUE,
    col_types = cols(
        `ED FOXA1` = col_number(),
        `Dox FOXA1 High` = col_number(),
        `Dox FOXA1 Low` = col_number(),
        `Dox H3K27ac` = col_number(),
        `ED H3K27ac` = col_number()
    ),
    
    show_col_types = TRUE
) %>%
    semi_join(high_confidence,
              by = c("name" = "high_confidence_t_rna_genes_hg19"))


write.table(easeq,"data/processed/hg19tRNA_genes_high_confidence.txt.txt", sep = "\t",row.names = FALSE, quote = FALSE)







#### normalised ####
easeq_normalised <- easeq %>%
    mutate_at(.vars = vars(starts_with("foxa1_ed")),
              .funs = funs(. - foxa1_ed_input)) %>%
    mutate_at(.vars = vars(starts_with("foxa1_dox")),
              .funs = funs(. - foxa1_dox_input))


normalised_fox <- easeq_normalised %>%
    select(!c("foxa1_ed_input", "foxa1_dox_input")) %>%
    pivot_longer(cols = starts_with("foxa1"),
                 names_to = "chip",
                 values_to = "quant") %>%
    filter(!str_detect(chip, 'h3k27ac'))


normalised_h3k27ac <- easeq_normalised %>%
    select(!c("foxa1_ed_input", "foxa1_dox_input")) %>%
    pivot_longer(cols = starts_with("foxa1"),
                 names_to = "chip",
                 values_to = "quant") %>%
    filter(str_detect(chip, 'h3k27ac'))




#### heatmaps ####
heatmap_fox <- normalised_fox %>%
    arrange(desc(chip), quant) %>% 
    ggplot(aes(chip, name, fill = quant)) +
    geom_tile() +
    scale_x_discrete(labels = c("FOXA1 Rep1 DOX+", "FOXA1 Rep1 DOX+", "H3K27ac DOX-")) +
theme(axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank())



heatmap_h3k27ac <- normalised_h3k27ac[order(normalised_h3k27ac$quant),] %>%
    ggplot(aes(chip, name, fill = quant)) +
    geom_tile() +
    scale_x_discrete(labels = c("H3K27ac DOX+", "H3K27ac DOX-")) +
    theme(axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank())


#### diff ####


diff_fox <- easeq_normalised %>%
    mutate_at(.vars = vars(starts_with("foxa1")),
              .funs = funs(. - foxa1_ed_foxa1_chipseq)) %>%
    select(
        c(
            "chromosome",
            "start",
            "end",
            "strand",
            "name",
            "foxa1_dox_foxa1_chipseq_rep1",
            "foxa1_dox_foxa1_chipseq_rep2"
        )
    ) %>%
    mutate(
        foxa1_dox_foxa1_chipseq_rep1 = if_else(foxa1_dox_foxa1_chipseq_rep1 > 0, "up", "dn", missing = "N/A")
    ) %>%
    mutate(
        foxa1_dox_foxa1_chipseq_rep2 = if_else(foxa1_dox_foxa1_chipseq_rep2 > 0, "up", "dn", missing = "N/A")
    )


diff_h3k27ac <- easeq_normalised %>%
    mutate_at(.vars = vars(starts_with("foxa1")),
              .funs = funs(. - foxa1_ed_h3k27ac_chipseq)) %>%
    select(c(
        "chromosome",
        "start",
        "end",
        "strand",
        "name",
        "foxa1_dox_h3k27ac_chipseq"
    )) %>%
    mutate(
        foxa1_dox_h3k27ac_chipseq = if_else(foxa1_dox_h3k27ac_chipseq > 0, "up", "dn", missing = "N/A")
    )

diff <-
    full_join(diff_fox,
              diff_h3k27ac,
              by = c("name",
                     "chromosome",
                     "start",
                     "end",
                     "strand"))

#### gain/loss ####

gain_loss_foxa1 <- normalised_fox %>%
    mutate(status = if_else(quant > median(quant), "gain", "loss", missing = "shared"))


gain_loss_h3k27ac <- normalised_h3k27ac %>%
    mutate(status = if_else(quant > median(quant), "gain", "loss", missing = "shared"))



gain_loss <-
    full_join(
        gain_loss_foxa1,
        gain_loss_h3k27ac,
        by = c(
            "name",
            "chip",
            "quant",
            "chromosome",
            "start",
            "end",
            "strand",
            "status"
        )
    )



#### venn ####

rep1 <- as.list(
    gain_loss %>%
        filter(!str_detect(chip, "ed|rep2")) %>%
        filter(str_detect(status, "gain")) %>%
        pivot_wider(names_from = chip, values_from = name) %>%
        select(
            c("foxa1_dox_foxa1_chipseq_rep1", "foxa1_dox_h3k27ac_chipseq")
        )
) %>%
    ggvenn(c("foxa1_dox_foxa1_chipseq_rep1", "foxa1_dox_h3k27ac_chipseq"))

