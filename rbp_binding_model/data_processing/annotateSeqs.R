library(GenomicRanges)
library(data.table)
library(annotatr)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library("org.Hs.eg.db")
library(dplyr)

file_to_annotate <- 'neg_rbp_seqs_150.csv'
rbp_seqs_df <- read.csv(file_to_annotate)

############## annotate gene expressin (lnTPM) #################
get_gene_expression <- function(cell_line) {
  cell_line_gene_df <- read.csv(paste0(cell_line, '_targetGene_all.csv'))
  # remove genes with missing lnTPM values
  cell_line_gene_df <- cell_line_gene_df[is.finite(cell_line_gene_df$lnTPM), ]
  cell_line_gene_GR <- makeGRangesFromDataFrame(cell_line_gene_df)
  
  sequence_df <- rbp_seqs_df[grep(cell_line, rbp_seqs_df$RBP), ]
  sequence_GR <- makeGRangesFromDataFrame(sequence_df, keep.extra.columns = FALSE)
  # find midpoint of sequence to annotate with gene expression
  start(sequence_GR) <- ((start(sequence_GR) + end(sequence_GR)) %/% 2)
  end(sequence_GR) <- start(sequence_GR)
  
  ix_list <- nearest(sequence_GR, cell_line_gene_GR, ignore.strand=TRUE)
  TPM_list <- cell_line_gene_df$lnTPM[ix_list]
  sequence_df$lnTPM <- TPM_list
  
  return(sequence_df)
}

cell_lines <- c('adrenal_gland','HepG2','K562')
annot_gene_seq_list <- list()
for (i in 1:length(cell_lines)) {
  annot_gene_seq_list[[i]] <- get_gene_expression(cell_lines[i])
}
annot_gene_seqs <- rbindlist(annot_gene_seq_list)

############# annotate sequence region ################
# Build the annotations (a single GRanges object) using annotatr
annotations_GR <- build_annotations(genome = 'hg19', annotations = 'hg19_basicgenes')
annot_gene_seqs_GR <- makeGRangesFromDataFrame(annot_gene_seqs, keep.extra.columns = TRUE)
start(annot_gene_seqs_GR) <- ((start(annot_gene_seqs_GR) + end(annot_gene_seqs_GR)) %/% 2)
end(annot_gene_seqs_GR) <- start(annot_gene_seqs_GR)

# Intersect the regions we read in with the annotations
dm_annotated <- annotate_regions(
  regions = annot_gene_seqs_GR,
  annotations = annotations_GR,
  ignore.strand = TRUE,
  quiet = FALSE)
# A GRanges object is returned
dm_annotated_df <- data.frame(dm_annotated)
dm_annotated_df <- dm_annotated_df[, c('RBP','seqnames','start','end','strand','lnTPM','annot.width','annot.type')]

# select annotation type from longest transcript for a given sequence
df_temp <- dm_annotated_df %>% group_by(RBP, seqnames, start, end, strand, lnTPM) %>% mutate(the_rank  = rank(-annot.width, ties.method = "random")) %>% filter(the_rank == 1) %>% dplyr::select(-the_rank)
df_temp <- df_temp[, c('RBP','seqnames','start','end','strand','lnTPM','annot.type')]

colnames(annot_gene_seqs) <- c('RBP','seqnames','true_start','true_end','strand','sequence','lnTPM')
annot_gene_seqs$start <- (annot_gene_seqs$true_start + annot_gene_seqs$true_end) %/% 2
annot_gene_seqs$end <- annot_gene_seqs$start

full_annot_seqs <- merge(x=annot_gene_seqs, y=df_temp, by=c('RBP','seqnames','start','end','strand','lnTPM'), all.x=TRUE)
names(full_annot_seqs)[names(full_annot_seqs) == "annot.type"] <- "annot"
full_annot_seqs <- full_annot_seqs[, c('RBP','seqnames','true_start','true_end','strand','sequence','lnTPM','annot')]
colnames(full_annot_seqs) <- c('RBP','seqnames','start','end','strand','sequence','lnTPM', 'annot')
write.csv(full_annot_seqs, file_to_annotate, row.names = FALSE)

# get only seq ranges and annotations (without actual sequence)
full_annot_seqs <- full_annot_seqs[, c('RBP','seqnames','start','end','strand','lnTPM','annot')]
write.csv(full_annot_seqs, '1kb_seqs/neg_rbp_seqs_1kb.csv', row.names = FALSE)
