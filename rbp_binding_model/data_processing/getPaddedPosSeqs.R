library(GenomicRanges)
library("BSgenome.Hsapiens.UCSC.hg19")

# get seqs from 150bp ranges file
pos_ranges_df <- read.csv('pos_ranges_150bp.csv')
pos_ranges_GR <- makeGRangesFromDataFrame(pos_ranges_df, keep.extra.columns=TRUE)

# get midpoint, then subtract 500bp for start, and scale to 1kb sequence length
start(pos_ranges_GR) <- ((start(pos_ranges_GR) + end(pos_ranges_GR)) %/% 2) - 500
end(pos_ranges_GR) <- start(pos_ranges_GR) + 999
pos_ranges_1kb_df <- as.data.frame(pos_ranges_GR)
write.csv(pos_ranges_1kb_df, 'pos_rbp_seqs_1kb.csv', row.names = FALSE)

# combine pos ranges with sequences, and reformat dataframe
seqs_1kb <- as.data.frame(getSeq(Hsapiens, pos_ranges_GR))
pos_ranges_1kb_df$sequence <- seqs_1kb$x #x is column name of seq column in df
pos_ranges_1kb_df <- pos_ranges_1kb_df[, c('RBP','seqnames','start','end','strand','sequence')]
write.csv(pos_ranges_1kb_df, 'pos_rbp_seqs_1kb.csv', row.names = FALSE)
