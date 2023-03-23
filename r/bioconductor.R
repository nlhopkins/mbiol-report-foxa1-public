library("tidyverse")
BiocManager::install("Rsamtools")
library(Rsamtools)
library(GenomicAlignments)
library(googledrive)



ed_fox <- scanBam(BamFile("~/Library/CloudStorage/GoogleDrive-nlh524@york.ac.uk/My Drive/Stage 4/Research Project/EaSeq/FOXA1/Datasets/FOXA1-Dox_FOXA1_ChIPSeq.bam"))

ed_h3 <- "~/Library/CloudStorage/GoogleDrive-nlh524@york.ac.uk/My Drive/Stage 4/Research Project/EaSeq/FOXA1/Datasets/FOXA1-Dox_H3K27ac_ChIPSeq.bam"

dox_fox <- "~/Library/CloudStorage/GoogleDrive-nlh524@york.ac.uk/My Drive/Stage 4/Research Project/EaSeq/FOXA1/Datasets/FOXA1+Dox_FOXA1_ChIPSeq rep2.bam"

dox_h3 <- "~/Library/CloudStorage/GoogleDrive-nlh524@york.ac.uk/My Drive/Stage 4/Research Project/EaSeq/FOXA1/Datasets/FOXA1+Dox_H3K27ac_ChIPSeq.bam"

ed_input <- "~/Library/CloudStorage/GoogleDrive-nlh524@york.ac.uk/My Drive/Stage 4/Research Project/EaSeq/FOXA1/Datasets/FOXA1-Dox_input.bam"

dox_input <- "~/Library/CloudStorage/GoogleDrive-nlh524@york.ac.uk/My Drive/Stage 4/Research Project/EaSeq/FOXA1/Datasets/FOXA1+Dox_input.bam"