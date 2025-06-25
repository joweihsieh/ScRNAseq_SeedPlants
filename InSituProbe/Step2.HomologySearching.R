# Homology and Unique Region Detection Pipeline for Probe Design
# Author: Jo-Wei Allison Hsieh
# Description: This script performs BLASTn analysis of candidate CDS sequences, filters homologous hits, checks for self-alignment, and extracts unique probe regions for in situ design

# ================================
# Run BLASTn
# ================================

#ln -s /home/f06b22037/SSD2/JW/genome/Cunninghamia_lanceolata/Chr_genome_final_gene.gtf ./
#ln -s /home/woodydrylab/HDD/GenomicsData/Cunninghamia_lanceolata/Cunninghamia_lanceolata_Shuai/Lachesis_assembly_changed.fa ./
#makeblastdb -in Lachesis_assembly_changed.fa -dbtype nucl -out Lachesis_db

#blastn \
#  -task blastn-short \
#  -query maker_genes_37_CDS.fa \
#  -db Lachesis_db \
#  -word_size 10 \
#  -dust no \
#  -perc_identity 70 \
#  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
#  -out blast_results_w10_s70.txt

# ================================
# Parse BLAST Results and Filter
# ================================
library(tidyverse)
library(Biostrings)

probe_length <- 45
blast <- read_tsv("blast_results_w10_s70.txt", 
                  col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"))

GTF <- read.table("Chr_genome_final_gene.gtf", sep = "\t")
filtered_GTF <- GTF %>%
  filter(V3 == "transcript" & grepl(paste(unique(blast$qseqid), collapse = "|"), V9))
filtered_GTF$transcript_id <- sub(".*transcript_id ([^;]+);.*", "\\1", filtered_GTF$V9)

merged_data <- merge(blast, filtered_GTF, by.x = "qseqid", by.y = "transcript_id")
merged_data$remark <- ifelse(
  merged_data$sstart >= merged_data$V4 & merged_data$send <= merged_data$V5,
  "itself", "not")

merged_data_v2 <- merged_data %>% 
  filter(pident >= 70 & length >= probe_length)

CDS_file <- "maker_genes_37_CDS.fa"
CDS_seq <- readDNAStringSet(CDS_file)
CDS_seq_df <- data.frame(names = names(CDS_seq), width = width(CDS_seq))
merged_data_v3 <- merge(merged_data_v2, CDS_seq_df, by.x = "qseqid", by.y = "names")
write.csv(merged_data_v3[, c(1:13, 16:17, 22, 23)], "maker_genes_37_CDS_homo.csv", quote = FALSE, row.names = FALSE)

# ================================
# Identify Unique (Non-Homologous) Regions
# ================================
seq_length_df <- data.frame(gene = names(CDS_seq), length = width(CDS_seq))

homology_regions <- merged_data_v3 %>% 
  filter(remark != "itself") %>%
  select(qseqid, qstart, qend) %>% 
  arrange(qseqid, qstart)

find_unique_regions <- function(gene, homology_df, seq_length_df){
  len <- seq_length_df$length[seq_length_df$gene == gene]
  full_region <- rep(TRUE, len)
  
  homo_sub <- homology_df %>% 
    filter(qseqid == gene, !is.na(qstart), !is.na(qend))
  
  if(nrow(homo_sub) == 0){
    return(data.frame(gene = gene, start = 1, end = len, unique = "full"))
  }
  
  for(i in 1:nrow(homo_sub)){
    region_start <- max(1, homo_sub$qstart[i])
    region_end <- min(len, homo_sub$qend[i])
    full_region[region_start:region_end] <- FALSE
  }
  
  rle_regions <- rle(full_region)
  ends <- cumsum(rle_regions$lengths)
  starts <- ends - rle_regions$lengths + 1
  
  unique_regions <- data.frame(start = starts[rle_regions$values],
                               end = ends[rle_regions$values])

  if(nrow(unique_regions) == 0){
    return(data.frame(gene = gene, start = NA, end = NA, unique = "none"))
  } else {
    unique_regions$gene <- gene
    unique_regions$unique <- "partial"

    if(nrow(unique_regions) == 1 && unique_regions$start[1] == 1 && unique_regions$end[1] == len){
      unique_regions$unique <- "full"
    }
  }
  
  return(unique_regions)
}

gene_list <- unique(seq_length_df$gene)

unique_regions_list <- lapply(gene_list, function(g){
  regions <- find_unique_regions(g, homology_regions, seq_length_df)
  return(regions)
})

unique_regions_df <- bind_rows(unique_regions_list) %>% select(gene, start, end, unique)

#table(unique_regions_df$unique)

unique_regions_df <- unique_regions_df %>%
  mutate(
    start = ifelse(is.na(start), 0, start),
    end = ifelse(is.na(end), 0, end)
  )

unique_regions_df$len <- unique_regions_df$end - unique_regions_df$start

unique_regions_df <- unique_regions_df %>%
  mutate(ID = sub("\\..*", "", gene))



FC <- read.csv("UMAP_log2FC_&_2others_clean.csv")
Table <- merge(FC, unique_regions_df, by = "ID")
Table2 <- merge(Table, seq_length_df, by = "gene")
colnames(Table2)[ncol(Table2)] <- "total_len"

Table_ordered <- Table2 %>% 
  mutate(ID = factor(ID, levels = unique(FC$ID))) %>% 
  arrange(ID)

# ================================
# Extract Probe Sequences
# ================================
fasta_seqs <- readDNAStringSet("maker_genes_37_CDS.fa")
names(fasta_seqs) <- sub("\\..*", "", names(fasta_seqs))
Table_ordered$sequence_fragment <- NA_character_

for (i in 1:nrow(Table_ordered)) {
  gene_id <- Table_ordered$ID[i]
  seq_start <- Table_ordered$start[i]
  seq_end <- Table_ordered$end[i]
  full_seq <- fasta_seqs[[gene_id]]
  if (seq_start >= 1 && seq_end <= length(full_seq)) {
    Table_ordered$sequence_fragment[i] <- as.character(subseq(full_seq, start = seq_start, end = seq_end))
  }
}

Table_ordered2 <- Table_ordered[Table_ordered$len >= probe_length, ]
write.csv(Table_ordered2, "maker_genes_CDS_unique_regions_w10_min45bp.csv", quote = FALSE, row.names = FALSE)

valid_seqs <- Table_ordered2 %>% filter(!is.na(sequence_fragment))
output_seqs <- DNAStringSet(valid_seqs$sequence_fragment)
names(output_seqs) <- paste0(valid_seqs$gene, "_", valid_seqs$start, "_", valid_seqs$end)
writeXStringSet(output_seqs, filepath = "maker_genes_unique_CDS_w10_min45bp.fa")

# ================================
# Summary Table of Gene Counts by Cluster
# ================================
gene_count_by_type <- Table_ordered2 %>% group_by(Type) %>% summarize(n_genes = n_distinct(gene))
