library(GenomicRanges)
library(IRanges)
library("BSgenome.Hsapiens.UCSC.hg19")
library(dplyr)
library(data.table)

######### get eligible regions to sample for negatives ##########
# start with list of transcribed regions downloaded from GENCODE
gencode_file = 'gencode.v34lift37.annotation.gff3'
read.delim(gencode_file, header=F, comment.char="#") -> gff
trans <- gff[gff$V3 == 'transcript', ]
trans <- trans[, c(1,4,5,7)]
trans <- trans[trans$seqnames != 'chrM', ]
names(trans) <- c('seqnames','start','end','strand')
trans_GR <- makeGRangesFromDataFrame(trans)
trans_GR <- intersect(trans_GR, trans_GR)

# allow 250bp overlap with nontranscribed regions (max < 500bp if padding overlaps)
start(trans_GR) <- start(trans_GR) - 250
end(trans_GR) <- end(trans_GR) + 250
trans_GR <- intersect(trans_GR, trans_GR)
trans_GR <- trans_GR[width(trans_GR) > 999]

# set fixed truncated windows of genome
neg_windows <- slidingWindows(trans_GR, width = 1000, step = 100)
windows_good <- unlist(neg_windows[width(neg_windows) == 1000])
windows_bad <- unlist(neg_windows[width(neg_windows) < 1000])
end(windows_bad) <- start(windows_bad) + 999
fixed_windows <- c(windows_bad, windows_good)
rm(windows_bad)
rm(windows_good)

# calculate gc values for fixed windows
window_seqs <- getSeq(Hsapiens, fixed_windows)
ff = alphabetFrequency(window_seqs, baseOnly=TRUE)
gc_array = (ff[,"C"]+ff[,"G"])/rowSums(ff)

# get partitioned dfs by gc to sample negative sequences from
neg_seqs_0 = fixed_windows[gc_array < 0.2]
neg_seqs_2 = fixed_windows[(gc_array >= 0.2) & (gc_array < 0.4)]
neg_seqs_4 = fixed_windows[(gc_array >= 0.4) & (gc_array < 0.6)]
neg_seqs_6 = fixed_windows[(gc_array >= 0.6) & (gc_array < 0.8)]
neg_seqs_8 = fixed_windows[gc_array >= 0.8]


######### gen neg seqs that are gc matched with pos seqs #########
# function to generate negative sequences for a given RBP
gen_neg_seqs <- function(RBP) {
  # get pos seqs for specified RBP and calculate gc distribution
  pos_RBP <- pos_rbp_seqs[pos_rbp_seqs$RBP == RBP ,]
  pos_RBP_GR <- makeGRangesFromDataFrame(pos_RBP)
  pos_seqs <- getSeq(Hsapiens, pos_RBP_GR)
  pos_ff = alphabetFrequency(pos_seqs, baseOnly=TRUE)
  pos_gc_array = (pos_ff[,"C"]+pos_ff[,"G"])/rowSums(pos_ff)
  
  # remove pos ranges from neg ranges with gc. allow for 200bp overlap with pos sequences
  good_ranges_0 = subsetByOverlaps(neg_seqs_0, pos_RBP_GR, minoverlap=200, invert = TRUE, ignore.strand=TRUE)
  good_ranges_2 = subsetByOverlaps(neg_seqs_2, pos_RBP_GR, minoverlap=200, invert = TRUE, ignore.strand=TRUE)
  good_ranges_4 = subsetByOverlaps(neg_seqs_4, pos_RBP_GR, minoverlap=200, invert = TRUE, ignore.strand=TRUE)
  good_ranges_6 = subsetByOverlaps(neg_seqs_6, pos_RBP_GR, minoverlap=200, invert = TRUE, ignore.strand=TRUE)
  good_ranges_8 = subsetByOverlaps(neg_seqs_8, pos_RBP_GR, minoverlap=200, invert = TRUE, ignore.strand=TRUE)
  
  # sample for random neg ranges in proper gc bucket
  out_0 = sample(good_ranges_0, sum(pos_gc_array < 0.2), replace = TRUE)
  out_2 = sample(good_ranges_2, sum((pos_gc_array >= 0.2) & (pos_gc_array < 0.4)), replace = TRUE)
  out_4 = sample(good_ranges_4, sum((pos_gc_array >= 0.4) & (pos_gc_array < 0.6)), replace = TRUE)
  out_6 = sample(good_ranges_6, sum((pos_gc_array >= 0.6) & (pos_gc_array < 0.8)), replace = TRUE)
  if (length(good_ranges_8) == 0) {
    out_8 = sample(good_ranges_6, sum(pos_gc_array >= 0.8), replace = TRUE)
  } else {
    out_8 = sample(good_ranges_8, sum(pos_gc_array >= 0.8), replace = TRUE)
  }
  out_all = c(out_0, out_2, out_4, out_6, out_8)
  
  # set random padding for neg seqs within window
  pad <- sample(0:99, length(out_all), replace = TRUE)
  start(out_all) <- start(out_all) + pad
  end(out_all) <- end(out_all) + pad
  
  # get actual 1kb neg seqs from granges object
  strand(out_all) <- sample(c('-','+'), size = length(out_all), replace = TRUE)
  out_seqs <- getSeq(Hsapiens, out_all)
  df_out <- as.data.frame(out_all)
  df_out$sequence <- as.data.frame(out_seqs)$x
  df_out$RBP <- RBP
  df_out <- df_out[, c(7,1,2,3,5,6)]
  return(df_out)
}

pos_rbp_seqs <- read.csv('pos_rbp_seqs_1kb.csv')
RBP_list <- as.character(unique(pos_rbp_seqs[['RBP']]))

neg_seqs_per_rbp_list <- list()
for (i in 1:length(RBP_list)) {
  neg_seqs_per_rbp_list[[i]] <- gen_neg_seqs(RBP_list[i])
}
neg_rbp_seqs <- rbindlist(neg_seqs_per_rbp_list)
write.csv(neg_rbp_seqs, 'neg_rbp_seqs_1kb_SEQ.csv', row.names = FALSE)

# get only neg seq ranges
neg_rbp_seqs <- neg_rbp_seqs[ , c('RBP','seqnames','start','end','strand')]
write.csv(neg_rbp_seqs, 'neg_rbp_seqs_1kb.csv', row.names = FALSE)
