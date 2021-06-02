library(gtools)
library(seqinr)
library("stats")
library("xlsx")
#library("calibrate")
#library("plotrix")
library("RColorBrewer")
library(igraph)
library("readxl")
source("alex_suite.R")

#################################################################################################################
# PARAMETERS
#################################################################################################################
# Delimiter between folders/files in filepaths. Should be "/" for Linux and "\\" for Windows.
FILEPATH_DELIM <- "/"
# Path to data folder, allowing simpler loading of data with the alex_suite.R data_path function.
DATA_FOLDER <- "/home/local/ARCS/ak3792/Documents/Research/data"
# Path to output folder, allowing simpler writing of output with the alex_suite.R output_path function.
OUTPUT_FOLDER <- "/home/local/ARCS/ak3792/Documents/Research/WGS/output"

# Load Ensembl to gene name mappings.
ensembl_genes <- get_ensembl_name_mappings()

#################################################################################################################
# Define top (HIGH_EXPRESSION_PERCENTILE)% HDE and HHE genes, and store the respective genebodies.
# Forces number of HDE and HHE genes to be the same (NUM_HE_GENES) for better cross-disease result comparison.
#################################################################################################################
if(TRUE) {
    genebody <- read.csv(data_path("refGene\\refGene_hg19_genebody_fixed.bed"), sep="\t", stringsAsFactors=FALSE)
    genebody_granges <- to_genomic_regions(genebody, chr_colname="chromosome", start_colname="tss", end_colname="tes", strand_colname="strand", label_colname="gene", order_coordinates=TRUE, remove_duplicate_labels=TRUE)
    
    genes_by_expression <- get_genes_by_expression()
    HDE_genes <- genes_by_expression[["HDE"]]; HHE_genes <- genes_by_expression[["HHE"]]
    HMHDE_genes <- genes_by_expression[["HMHDE"]]; HMHHE_genes <- genes_by_expression[["HMHHE"]]
    MHDE_genes <- genes_by_expression[["MHDE"]]; MHHE_genes <- genes_by_expression[["MHHE"]]
    MLDE_genes <- genes_by_expression[["MLDE"]]; MLHE_genes <- genes_by_expression[["MLHE"]]
    LDE_genes <- genes_by_expression[["LDE"]]; LHE_genes <- genes_by_expression[["LHE"]]
    constrained_genes <- ls(get_constrained_genes("pLI>0.5"))
    unconstrained_genes <- ls(get_constrained_genes("pLI<=0.5"))
    bivalent_genes <- get_bivalent_genes()
    candidate_CDH_genes <- unique(read.csv(data_path("Table S5 CDH training genes with mouse model.csv"), stringsAsFactors=FALSE, header=TRUE, skip=1)[,1]) # 61 candidate CDH genes, from Lan's file
    candidate_CHD_genes <- unique(read.csv(data_path("mouse_CHD_gene0516.txt"), sep="\t", stringsAsFactors=FALSE)[,1]) # Latest list of 727 candidate CHD genes.
    known_mouse_ko_genes <- unique(read.csv(data_path("Curated known Human-Mouse CHD genes.txt"), sep="\t", stringsAsFactors=FALSE)[,1]) # 253 known CHD mouse KO genes, from Lan's file
    
    genebody_HDE <- genebody[genebody$gene %in% HDE_genes,]; genebody_HHE <- genebody[genebody$gene %in% HHE_genes,]
    genebody_HMHDE <- genebody[genebody$gene %in% HMHDE_genes,]; genebody_HMHHE <- genebody[genebody$gene %in% HMHHE_genes,]
    genebody_MHDE <- genebody[genebody$gene %in% MHDE_genes,]; genebody_MHHE <- genebody[genebody$gene %in% MHHE_genes,]
    genebody_MLDE <- genebody[genebody$gene %in% MLDE_genes,]; genebody_MLHE <- genebody[genebody$gene %in% MLHE_genes,]
    genebody_LDE <- genebody[genebody$gene %in% LDE_genes,]; genebody_LHE <- genebody[genebody$gene %in% LHE_genes,]
    genebody_constrained <- genebody[genebody$gene %in% constrained_genes,]
    genebody_unconstrained <- genebody[genebody$gene %in% unconstrained_genes,]
    genebody_bivalent <- genebody[genebody$gene %in% bivalent_genes,]
    genebody_candidate_CDH <- genebody[genebody$gene %in% candidate_CDH_genes,]
    genebody_candidate_CHD <- genebody[genebody$gene %in% candidate_CHD_genes,]
    genebody_known_mouse_ko <- genebody[genebody$gene %in% known_mouse_ko_genes,]
}

#################################################################################################################
# Find HDE, MHDE, MLDE, LDE variants in all of the data sets.
# Define all (ex: HDE, HHE, constrained, etc) TSS/TES(3'UTR) regions.
#################################################################################################################
# Get particularly defined TS regions, given the genebody and left and right pillow sizes around the specified "TSS" or "TES" sites.
get_TS_regions <- function(genebody, sites, left, right) {
    strands <- rep(1, nrow(genebody)); strands[genebody[,6]=="-"] <- -1
    if (sites=="TSS") { TS_index = 2 } else if (sites=="TES") { TS_index = 3 } else { print("ERROR: Invalid sites parameter."); return(0) }
    TS_regions <- unique(data.frame(gsub("chr", "", genebody[,1]), genebody[,TS_index]-(left*strands), genebody[,TS_index]+(right*strands), genebody[,4]))
    colnames(TS_regions) <- c("chromosome", "start", "end", "gene")
    ordered_coordinates <- t(apply(TS_regions[,c("start", "end")], 1, sort))
    TS_regions$start <- ordered_coordinates[,1]
    TS_regions$end <- ordered_coordinates[,2]
    TS_regions <- TS_regions[!duplicated(TS_regions$gene),]
    return(TS_regions)
}
# Define all TES/TSS region variables.
if(TRUE) {
    TSS40_regions <- get_TS_regions(genebody, sites="TSS", left=40000, right=40000); TSS_regions <- get_TS_regions(genebody, sites="TSS", left=20000, right=20000); TES_regions <- get_TS_regions(genebody, sites="TES", left=5000, right=20000)
    Felix_TSS_regions <- get_TS_regions(genebody, sites="TSS", left=10000, right=0) 
    HDE_TSS40_regions <- TSS40_regions[TSS40_regions$gene %in% HDE_genes,]; HDE_TSS_regions <- TSS_regions[TSS_regions$gene %in% HDE_genes,]; HDE_TES_regions <- TES_regions[TES_regions$gene %in% HDE_genes,]
    HHE_TSS40_regions <- TSS40_regions[TSS40_regions$gene %in% HHE_genes,]; HHE_TSS_regions <- TSS_regions[TSS_regions$gene %in% HHE_genes,]; HHE_TES_regions <- TES_regions[TES_regions$gene %in% HHE_genes,]
    HMHDE_TSS40_regions <- TSS40_regions[TSS40_regions$gene %in% HMHDE_genes,]; HMHDE_TSS_regions <- TSS_regions[TSS_regions$gene %in% HMHDE_genes,]; HMHDE_TES_regions <- TES_regions[TES_regions$gene %in% HMHDE_genes,]
    HMHHE_TSS40_regions <- TSS40_regions[TSS40_regions$gene %in% HMHHE_genes,]; HMHHE_TSS_regions <- TSS_regions[TSS_regions$gene %in% HMHHE_genes,]; HMHHE_TES_regions <- TES_regions[TES_regions$gene %in% HMHHE_genes,]
    MHDE_TSS40_regions <- TSS40_regions[TSS40_regions$gene %in% MHDE_genes,]; MHDE_TSS_regions <- TSS_regions[TSS_regions$gene %in% MHDE_genes,]; MHDE_TES_regions <- TES_regions[TES_regions$gene %in% MHDE_genes,]
    MHHE_TSS40_regions <- TSS40_regions[TSS40_regions$gene %in% MHHE_genes,]; MHHE_TSS_regions <- TSS_regions[TSS_regions$gene %in% MHHE_genes,]; MHHE_TES_regions <- TES_regions[TES_regions$gene %in% MHHE_genes,]
    MLDE_TSS40_regions <- TSS40_regions[TSS40_regions$gene %in% MLDE_genes,]; MLDE_TSS_regions <- TSS_regions[TSS_regions$gene %in% MLDE_genes,]; MLDE_TES_regions <- TES_regions[TES_regions$gene %in% MLDE_genes,]
    MLHE_TSS40_regions <- TSS40_regions[TSS40_regions$gene %in% MLHE_genes,]; MLHE_TSS_regions <- TSS_regions[TSS_regions$gene %in% MLHE_genes,]; MLHE_TES_regions <- TES_regions[TES_regions$gene %in% MLHE_genes,]
    LDE_TSS40_regions <- TSS40_regions[TSS40_regions$gene %in% LDE_genes,]; LDE_TSS_regions <- TSS_regions[TSS_regions$gene %in% LDE_genes,]; LDE_TES_regions <- TES_regions[TES_regions$gene %in% LDE_genes,]
    LHE_TSS40_regions <- TSS40_regions[TSS40_regions$gene %in% LHE_genes,]; LHE_TSS_regions <- TSS_regions[TSS_regions$gene %in% LHE_genes,]; LHE_TES_regions <- TES_regions[TES_regions$gene %in% LHE_genes,]
    constrained_TSS40_regions <- TSS40_regions[TSS40_regions$gene %in% constrained_genes,]; constrained_TSS_regions <- TSS_regions[TSS_regions$gene %in% constrained_genes,]; constrained_TES_regions <- TES_regions[TES_regions$gene %in% constrained_genes,];
    unconstrained_TSS40_regions <- TSS40_regions[TSS40_regions$gene %in% unconstrained_genes,]; unconstrained_TSS_regions <- TSS_regions[TSS_regions$gene %in% unconstrained_genes,]; unconstrained_TES_regions <- TES_regions[TES_regions$gene %in% unconstrained_genes,];
    candidate_CHD_TSS40_regions <- TSS40_regions[TSS40_regions$gene %in% candidate_CHD_genes,]; candidate_CHD_TSS_regions <- TSS_regions[TSS_regions$gene %in% candidate_CHD_genes,]; candidate_CHD_TES_regions <- TES_regions[TES_regions$gene %in% candidate_CHD_genes,];
    known_mouse_ko_TSS40_regions <- TSS40_regions[TSS40_regions$gene %in% known_mouse_ko_genes,]; known_mouse_ko_TSS_regions <- TSS_regions[TSS_regions$gene %in% known_mouse_ko_genes,]; known_mouse_ko_TES_regions <- TES_regions[TES_regions$gene %in% known_mouse_ko_genes,];
}
# Define all TES/TSS GRanges objects, for more efficient overlap calculation.
if(TRUE) {
    TSS40_granges <- to_genomic_regions(TSS40_regions, labels=TSS40_regions$gene); TSS_granges <- to_genomic_regions(TSS_regions, labels=TSS_regions$gene); TES_granges <- to_genomic_regions(TES_regions, labels=TES_regions$gene);
    start(TSS40_granges)[start(TSS40_granges) < 1] <- 1; start(TSS_granges)[start(TSS_granges) < 1] <- 1; start(TES_granges)[start(TES_granges) < 1] <- 1;
    Felix_TSS_granges <- to_genomic_regions(Felix_TSS_regions, labels=Felix_TSS_regions$gene)
    HDE_TSS40_granges <- TSS40_granges[names(TSS40_granges) %in% HDE_genes]; HDE_TSS_granges <- TSS_granges[names(TSS_granges) %in% HDE_genes]; HDE_TES_granges <- TES_granges[names(TES_granges) %in% HDE_genes]
    HHE_TSS40_granges <- TSS40_granges[names(TSS40_granges) %in% HHE_genes]; HHE_TSS_granges <- TSS_granges[names(TSS_granges) %in% HHE_genes,]; HHE_TES_granges <- TES_granges[names(TES_granges) %in% HHE_genes,]
    HMHDE_TSS40_granges <- TSS40_granges[names(TSS40_granges) %in% HMHDE_genes]; HMHDE_TSS_granges <- TSS_granges[names(TSS_granges) %in% HMHDE_genes,]; HMHDE_TES_granges <- TES_granges[names(TES_granges) %in% HMHDE_genes,]
    HMHHE_TSS40_granges <- TSS40_granges[names(TSS40_granges) %in% HMHHE_genes]; HMHHE_TSS_granges <- TSS_granges[names(TSS_granges) %in% HMHHE_genes,]; HMHHE_TES_granges <- TES_granges[names(TES_granges) %in% HMHHE_genes,]
    MHDE_TSS40_granges <- TSS40_granges[names(TSS40_granges) %in% MHDE_genes]; MHDE_TSS_granges <- TSS_granges[names(TSS_granges) %in% MHDE_genes,]; MHDE_TES_granges <- TES_granges[names(TES_granges) %in% MHDE_genes,]
    MHHE_TSS40_granges <- TSS40_granges[names(TSS40_granges) %in% MHHE_genes]; MHHE_TSS_granges <- TSS_granges[names(TSS_granges) %in% MHHE_genes,]; MHHE_TES_granges <- TES_granges[names(TES_granges) %in% MHHE_genes,]
    MLDE_TSS40_granges <- TSS40_granges[names(TSS40_granges) %in% MLDE_genes]; MLDE_TSS_granges <- TSS_granges[names(TSS_granges) %in% MLDE_genes,]; MLDE_TES_granges <- TES_granges[names(TES_granges) %in% MLDE_genes,]
    MLHE_TSS40_granges <- TSS40_granges[names(TSS40_granges) %in% MLHE_genes]; MLHE_TSS_granges <- TSS_granges[names(TSS_granges) %in% MLHE_genes,]; MLHE_TES_granges <- TES_granges[names(TES_granges) %in% MLHE_genes,]
    LDE_TSS40_granges <- TSS40_granges[names(TSS40_granges) %in% LDE_genes]; LDE_TSS_granges <- TSS_granges[names(TSS_granges) %in% LDE_genes,]; LDE_TES_granges <- TES_granges[names(TES_granges) %in% LDE_genes,]
    LHE_TSS40_granges <- TSS40_granges[names(TSS40_granges) %in% LHE_genes]; LHE_TSS_granges <- TSS_granges[names(TSS_granges) %in% LHE_genes,]; LHE_TES_granges <- TES_granges[names(TES_granges) %in% LHE_genes,]
    constrained_TSS40_granges <- TSS40_granges[names(TSS40_granges) %in% constrained_genes]; constrained_TSS_granges <- TSS_granges[names(TSS_granges) %in% constrained_genes,]; constrained_TES_granges <- TES_granges[names(TES_granges) %in% constrained_genes,];
    unconstrained_TSS40_granges <- TSS40_granges[names(TSS40_granges) %in% unconstrained_genes]; unconstrained_TSS_granges <- TSS_granges[names(TSS_granges) %in% unconstrained_genes,]; unconstrained_TES_granges <- TES_granges[names(TES_granges) %in% unconstrained_genes,];
    bivalent_TSS40_granges <- TSS40_granges[names(TSS40_granges) %in% bivalent_genes]; bivalent_TSS_granges <- TSS_granges[names(TSS_granges) %in% bivalent_genes,]; bivalent_TES_granges <- TES_granges[names(TES_granges) %in% bivalent_genes,];
    candidate_CDH_TSS40_granges <- TSS40_granges[names(TSS40_granges) %in% candidate_CDH_genes]; candidate_CDH_TSS_granges <- TSS_granges[names(TSS_granges) %in% candidate_CDH_genes,]; candidate_CDH_TES_granges <- TES_granges[names(TES_granges) %in% candidate_CDH_genes,];
    candidate_CHD_TSS40_granges <- TSS40_granges[names(TSS40_granges) %in% candidate_CHD_genes]; candidate_CHD_TSS_granges <- TSS_granges[names(TSS_granges) %in% candidate_CHD_genes,]; candidate_CHD_TES_granges <- TES_granges[names(TES_granges) %in% candidate_CHD_genes,];
    known_mouse_ko_TSS40_granges <- TSS40_granges[names(TSS40_granges) %in% known_mouse_ko_genes]; known_mouse_ko_TSS_granges <- TSS_granges[names(TSS_granges) %in% known_mouse_ko_genes,]; known_mouse_ko_TES_granges <- TES_granges[names(TES_granges) %in% known_mouse_ko_genes,];
}
# Define Roadmap eid names
store_roadmap_eid_names()
# Define chromosome lengths
chr_lengths <- get_global("chr_lengths")
if(is.null(chr_lengths)) { chr_lengths <- hash_chromosome_lengths() }
# Define standard chromosomes
standard_chromosomes <- get_global("standard_chromosomes")
if(is.null(standard_chromosomes)) { standard_chromosomes <- c(1:22,"X","Y") }
# Define hg38 versions of mappability/region filters
all_filters_granges_hg38 <- get_global("all_filters_granges_hg38")
if(is.null(all_filters_granges_hg38)) { 
    if(file.exists(output_path("all_filters_hg38_granges.rds"))) {
        all_filters_granges_hg38 <- readRDS(output_path("all_filters_hg38_granges.rds")) #to_genomic_regions(genomic_coordinates_from_strings(read.csv(output_path("all_filters_hg38_granges.txt"),header=FALSE)[,1]))
    } else {
        # Define mappability
        if(TRUE) {
            mappability_hg38_dat <- read.table(data_path("PCGC/mappability/hg38_300bp.bed"))
            colnames(mappability_hg38_dat) <- c("chromosome", "start", "end", "mappability")
            mappability_hg38_dat <- mappability_hg38_dat[mappability_hg38_dat$chromosome %in% paste0("chr",standard_chromosomes),] 
            mappability_hg38_dat <- mappability_hg38_dat[mappability_hg38_dat$mappability != 1,]
            mappability_hg38_dat$chromosome <- gsub("chr", "", mappability_hg38_dat$chromosome)
            mappability_hg38_granges <- makeGRangesFromDataFrame(mappability_hg38_dat, starts.in.df.are.0based=TRUE)
            mappability_hg38_granges <- intersect(mappability_hg38_granges, mappability_hg38_granges)
            rm(mappability_hg38_dat)
            write(paste0("chr",mappability_hg38_granges), file=output_path("mappability_hg38_granges.txt"))
        }
        print(sum(as.numeric(width(mappability_hg38_granges))))
        
        # Define LCR
        if(TRUE) {
            lcr_hg38_dat <- read.table(data_path("PCGC/LCR/LCR-hs38.bed"))
            colnames(lcr_hg38_dat) <- c("chromosome", "start", "end")
            lcr_hg38_dat <- lcr_hg38_dat[lcr_hg38_dat$chromosome %in% paste0("chr",standard_chromosomes),] 
            lcr_hg38_dat$chromosome <- gsub("chr", "", lcr_hg38_dat$chromosome)
            lcr_hg38_granges <- makeGRangesFromDataFrame(lcr_hg38_dat, starts.in.df.are.0based=TRUE)
            lcr_hg38_granges <- intersect(lcr_hg38_granges, lcr_hg38_granges)
            rm(lcr_hg38_dat)
            write(paste0("chr",lcr_hg38_granges), file=output_path("lcr_hg38_granges.txt"))
        }
        print(sum(as.numeric(width(lcr_hg38_granges))))
        
        # Define segdups
        if(TRUE) {
            segdups_hg38_dat <- read.table(data_path("PCGC/segdups/hg38_genomicSuperDups.txt"))[,c(2:4,27)]
            colnames(segdups_hg38_dat) <- c("chromosome", "start", "end", "similarity")
            segdups_hg38_dat <- segdups_hg38_dat[segdups_hg38_dat$chromosome %in% paste0("chr",standard_chromosomes),] 
            segdups_hg38_dat <- segdups_hg38_dat[segdups_hg38_dat$similarity > 0.99,]
            segdups_hg38_dat$chromosome <- gsub("chr", "", segdups_hg38_dat$chromosome)
            segdups_hg38_granges <- makeGRangesFromDataFrame(segdups_hg38_dat, starts.in.df.are.0based=TRUE)
            segdups_hg38_granges <- intersect(segdups_hg38_granges, segdups_hg38_granges)
            rm(segdups_hg38_dat)
            write(paste0("chr",segdups_hg38_granges), file=output_path("segdups_hg38_granges.txt"))
        }
        print(sum(as.numeric(width(segdups_hg38_granges))))
        
        # Define ENCODE Duke Blacklist
        if(TRUE) {
            duke_blacklist_hg38_dat <- read.table(data_path("PCGC/ENCODE_blacklists/encode_duke_blacklist_hg38.bed"))[,c(1:3)]
            colnames(duke_blacklist_hg38_dat) <- c("chromosome", "start", "end")
            duke_blacklist_hg38_dat <- duke_blacklist_hg38_dat[duke_blacklist_hg38_dat$chromosome %in% paste0("chr",standard_chromosomes),] 
            duke_blacklist_hg38_dat$chromosome <- gsub("chr", "", duke_blacklist_hg38_dat$chromosome)
            duke_blacklist_hg38_granges <- makeGRangesFromDataFrame(duke_blacklist_hg38_dat, starts.in.df.are.0based=TRUE)
            duke_blacklist_hg38_granges <- intersect(duke_blacklist_hg38_granges, duke_blacklist_hg38_granges)
            rm(duke_blacklist_hg38_dat)
            write(paste0("chr",duke_blacklist_hg38_granges), file=output_path("duke_blacklist_hg38_granges.txt"))
        }
        print(sum(as.numeric(width(duke_blacklist_hg38_granges))))
        
        # Define ENCODE DAC Blacklist
        if(TRUE) {
            dac_blacklist_hg38_dat <- read.table(data_path("PCGC/ENCODE_blacklists/encode_dac_blacklist_hg38.bed"))[,c(1:3)]
            colnames(dac_blacklist_hg38_dat) <- c("chromosome", "start", "end")
            dac_blacklist_hg38_dat <- dac_blacklist_hg38_dat[dac_blacklist_hg38_dat$chromosome %in% paste0("chr",standard_chromosomes),] 
            dac_blacklist_hg38_dat$chromosome <- gsub("chr", "", dac_blacklist_hg38_dat$chromosome)
            dac_blacklist_hg38_granges <- makeGRangesFromDataFrame(dac_blacklist_hg38_dat, starts.in.df.are.0based=TRUE)
            dac_blacklist_hg38_granges <- intersect(dac_blacklist_hg38_granges, dac_blacklist_hg38_granges)
            rm(dac_blacklist_hg38_dat)
            write(paste0("chr",dac_blacklist_hg38_granges), file=output_path("dac_blacklist_hg38_granges.txt"))
        }
        print(sum(as.numeric(width(dac_blacklist_hg38_granges))))
        
        # Define combination of hg38 filters
        if(TRUE) {
            all_filters_granges_hg38 <- c(mappability_hg38_granges, lcr_hg38_granges, segdups_hg38_granges, duke_blacklist_hg38_granges, dac_blacklist_hg38_granges)
            all_filters_granges_hg38 <- intersect(all_filters_granges_hg38, all_filters_granges_hg38)
            write(paste0("chr",all_filters_granges_hg38), file=output_path("all_filters_hg38_granges.txt"))
            saveRDS(all_filters_granges_hg38, file=output_path("all_filters_hg38_granges.rds"))
        }
    }
}
print(paste0(sum(as.numeric(chr_lengths))-sum(as.numeric(width(all_filters_granges_hg38)))," bp (",round((sum(as.numeric(chr_lengths))-sum(as.numeric(width(all_filters_granges_hg38))))/sum(as.numeric(chr_lengths)),3)," of genome) callable regions remaining after ",sum(as.numeric(width(all_filters_granges_hg38)))," total bp filtered out."))
# Define hg19 versions of mappability/region filters
all_filters_granges_hg19 <- get_global("all_filters_granges_hg19")
if(is.null(all_filters_granges_hg19)) { 
    if(file.exists(output_path("all_filters_hg19_granges.rds"))) {
        all_filters_granges_hg19 <- readRDS(output_path("all_filters_hg19_granges.rds")) #to_genomic_regions(genomic_coordinates_from_strings(read.csv(output_path("all_filters_hg19_granges.txt"),header=FALSE)[,1]))
    } else {
        # Define mappability
        if(TRUE) {
            mappability_hg19_dat <- read.table(data_path("PCGC/mappability/hg19_300bp.bed"))
            colnames(mappability_hg19_dat) <- c("chromosome", "start", "end", "mappability")
            mappability_hg19_dat <- mappability_hg19_dat[mappability_hg19_dat$chromosome %in% paste0("chr",standard_chromosomes),] 
            mappability_hg19_dat <- mappability_hg19_dat[mappability_hg19_dat$mappability != 1,]
            mappability_hg19_dat$chromosome <- gsub("chr", "", mappability_hg19_dat$chromosome)
            mappability_hg19_granges <- makeGRangesFromDataFrame(mappability_hg19_dat, starts.in.df.are.0based=TRUE)
            mappability_hg19_granges <- intersect(mappability_hg19_granges, mappability_hg19_granges)
            rm(mappability_hg19_dat)
            write(paste0("chr",mappability_hg19_granges), file=output_path("mappability_hg19_granges.txt"))
        }
        print(sum(as.numeric(width(mappability_hg19_granges))))
        
        # Define LCR
        if(TRUE) {
            lcr_hg19_dat <- read.table(data_path("PCGC/LCR/LCR-hs37d5.bed"))
            colnames(lcr_hg19_dat) <- c("chromosome", "start", "end")
            lcr_hg19_dat <- lcr_hg19_dat[lcr_hg19_dat$chromosome %in% standard_chromosomes,] 
            lcr_hg19_dat$chromosome <- gsub("chr", "", lcr_hg19_dat$chromosome)
            lcr_hg19_granges <- makeGRangesFromDataFrame(lcr_hg19_dat, starts.in.df.are.0based=TRUE)
            lcr_hg19_granges <- intersect(lcr_hg19_granges, lcr_hg19_granges)
            rm(lcr_hg19_dat)
            write(paste0("chr",lcr_hg19_granges), file=output_path("lcr_hg19_granges.txt"))
        }
        print(sum(as.numeric(width(lcr_hg19_granges))))
        
        # Define segdups
        if(TRUE) {
            segdups_hg19_dat <- read.table(data_path("PCGC/segdups/hg19_genomicSuperDups.txt"))[,c(2:4,27)]
            colnames(segdups_hg19_dat) <- c("chromosome", "start", "end", "similarity")
            segdups_hg19_dat <- segdups_hg19_dat[segdups_hg19_dat$chromosome %in% paste0("chr",standard_chromosomes),] 
            segdups_hg19_dat <- segdups_hg19_dat[segdups_hg19_dat$similarity > 0.99,]
            segdups_hg19_dat$chromosome <- gsub("chr", "", segdups_hg19_dat$chromosome)
            segdups_hg19_granges <- makeGRangesFromDataFrame(segdups_hg19_dat, starts.in.df.are.0based=TRUE)
            segdups_hg19_granges <- intersect(segdups_hg19_granges, segdups_hg19_granges)
            rm(segdups_hg19_dat)
            write(paste0("chr",segdups_hg19_granges), file=output_path("segdups_hg19_granges.txt"))
        }
        print(sum(as.numeric(width(segdups_hg19_granges))))
        
        # Define ENCODE Duke Blacklist
        if(TRUE) {
            duke_blacklist_hg19_dat <- read.table(data_path("PCGC/ENCODE_blacklists/encode_duke_blacklist_hg19.bed"))[,c(1:3)]
            colnames(duke_blacklist_hg19_dat) <- c("chromosome", "start", "end")
            duke_blacklist_hg19_dat <- duke_blacklist_hg19_dat[duke_blacklist_hg19_dat$chromosome %in% paste0("chr",standard_chromosomes),] 
            duke_blacklist_hg19_dat$chromosome <- gsub("chr", "", duke_blacklist_hg19_dat$chromosome)
            duke_blacklist_hg19_granges <- makeGRangesFromDataFrame(duke_blacklist_hg19_dat, starts.in.df.are.0based=TRUE)
            duke_blacklist_hg19_granges <- intersect(duke_blacklist_hg19_granges, duke_blacklist_hg19_granges)
            rm(duke_blacklist_hg19_dat)
            write(paste0("chr",duke_blacklist_hg19_granges), file=output_path("duke_blacklist_hg19_granges.txt"))
        }
        print(sum(as.numeric(width(duke_blacklist_hg19_granges))))
        
        # Define ENCODE DAC Blacklist
        if(TRUE) {
            dac_blacklist_hg19_dat <- read.table(data_path("PCGC/ENCODE_blacklists/encode_dac_blacklist_hg19.bed"))[,c(1:3)]
            colnames(dac_blacklist_hg19_dat) <- c("chromosome", "start", "end")
            dac_blacklist_hg19_dat <- dac_blacklist_hg19_dat[dac_blacklist_hg19_dat$chromosome %in% paste0("chr",standard_chromosomes),] 
            dac_blacklist_hg19_dat$chromosome <- gsub("chr", "", dac_blacklist_hg19_dat$chromosome)
            dac_blacklist_hg19_granges <- makeGRangesFromDataFrame(dac_blacklist_hg19_dat, starts.in.df.are.0based=TRUE)
            dac_blacklist_hg19_granges <- intersect(dac_blacklist_hg19_granges, dac_blacklist_hg19_granges)
            rm(dac_blacklist_hg19_dat)
            write(paste0("chr",dac_blacklist_hg19_granges), file=output_path("dac_blacklist_hg19_granges.txt"))
        }
        print(sum(as.numeric(width(dac_blacklist_hg19_granges))))
        
        # Define combination of hg19 filters
        if(TRUE) {
            all_filters_granges_hg19 <- c(mappability_hg19_granges, lcr_hg19_granges, segdups_hg19_granges, duke_blacklist_hg19_granges, dac_blacklist_hg19_granges)
            all_filters_granges_hg19 <- intersect(all_filters_granges_hg19, all_filters_granges_hg19)
            write(paste0("chr",all_filters_granges_hg19), file=output_path("all_filters_hg19_granges.txt"))
            saveRDS(all_filters_granges_hg19, file=output_path("all_filters_hg19_granges.rds"))
        }
    }
}
print(paste0(sum(as.numeric(chr_lengths))-sum(as.numeric(width(all_filters_granges_hg19)))," bp (",round((sum(as.numeric(chr_lengths))-sum(as.numeric(width(all_filters_granges_hg19))))/sum(as.numeric(chr_lengths)),3)," of genome) callable regions remaining after ",sum(as.numeric(width(all_filters_granges_hg19)))," total bp filtered out."))

#################################################################################################################
# Determine which variants fall in the defined HDE/HHE TSS/TES regions for cases and controls.
#################################################################################################################
# Returns true for variants that are near either transcription start sites (TSS) or transcription end sites (TES), and false otherwise for each index.
# Chromosomes and positions must be parallel vectors.
# If samples is specified, considers all variants/sample in one go. This helps overcome cases where an individual may have multiple variants in a single gene.
is_TS <- function(variant_granges, TS_granges=NULL, sites=NULL, expression=NULL, return_genes=FALSE, samples=NULL) {
    if(class(variant_granges)[1] != "GRanges") { print("ERROR: Parameter variant_granges is not a GRanges object!"); return() }
    if(is.null(TS_granges)) {
        TS_granges_name = paste(c(expression, sites, "granges"), collapse="_")
        TS_granges <- get(TS_granges_name)
    }
    
    hits_indices <- data.frame(findOverlaps(variant_granges, TS_granges))
    hits <- cbind(hits_indices[,1], names(TS_granges)[hits_indices[,2]])
    hits_genes <- sapply(seq(1:length(variant_granges)), function(i) { return(paste0(unique(hits[(hits[,1] == i),2]),collapse=",")) })
    
    if(return_genes) { return(hits_genes)
    } else { return(hits_genes != "") }
}

#################################################################################################################
# Define case and control data for both CDH and CHD, using standardized column names.
# Note that "ssc" refers to 438 SSC control variants called by Hongjian with GATK, as in the CDH calls.
# On the other hand, "sscfb" refers to the full set of SSC control variants, including Harvard FreeBayes calls, as in the CHD calls.
#################################################################################################################
# Returns separate WGSA and gene-annotated data tables for snps, indels, and all mutations for specifically the TES and TSS regions.
clean_wgsa_data <- function(snps, indels, snp_granges=NULL, indel_granges=NULL, handle_indel_brackets=TRUE) {
    # If snp_granges or indel_granges parallel GRanges object not passed in as paremeters, build them for use in this function.
    if(class(snp_granges)[1] != "GRanges") { snp_granges <- to_genomic_regions(snps, chr_colname="X.chr", start_colname="pos", end_colname="pos", labels=rep("snv",nrow(snps))) }
    if(class(indel_granges)[1] != "GRanges") { indel_granges <- to_genomic_regions(indels, chr_colname="X.chr", start_colname="pos", end_colname="pos", labels=rep("indel",nrow(indels))) }
    
    wgsa_data <- new.env()
    if(handle_indel_brackets) {
        indels <- unfactorize(indels)
        for(i in 1:nrow(indels)) {
            print(i)
            for(j in 6:ncol(indels)) { #355:ncol(indels)
                indels[i,j] <- gsub("{[0-9]+}", "", indels[i,j], perl=TRUE)
            }
        }
    }
    snps_TSS40_genes <- is_TS(snp_granges, sites="TSS40", return_genes=TRUE)
    snps_TSS_genes <- is_TS(snp_granges, sites="TSS", return_genes=TRUE)
    snps_TES_genes <- is_TS(snp_granges, sites="TES", return_genes=TRUE)
    indels_TSS40_genes <- is_TS(indel_granges, sites="TSS40", return_genes=TRUE)
    indels_TSS_genes <- is_TS(indel_granges, sites="TSS", return_genes=TRUE)
    indels_TES_genes <- is_TS(indel_granges, sites="TES", return_genes=TRUE)
    
    snps_TSS40 <- cbind(snps_TSS40_genes, snps)[snps_TSS40_genes != "",]; colnames(snps_TSS40)[1] <- "TS_gene"
    snps_TSS <- cbind(snps_TSS_genes, snps)[snps_TSS_genes != "",]; colnames(snps_TSS)[1] <- "TS_gene"
    snps_TES <- cbind(snps_TES_genes, snps)[snps_TES_genes != "",]; colnames(snps_TES)[1] <- "TS_gene"
    indels_TSS40 <- cbind(indels_TSS40_genes, indels)[indels_TSS40_genes != "",]; colnames(indels_TSS40)[1] <- "TS_gene"
    indels_TSS <- cbind(indels_TSS_genes, indels)[indels_TSS_genes != "",]; colnames(indels_TSS)[1] <- "TS_gene"
    indels_TES <- cbind(indels_TES_genes, indels)[indels_TES_genes != "",]; colnames(indels_TES)[1] <- "TS_gene"
    
    shared_features <- intersect(colnames(snps), colnames(indels))
    muts <- rbind(snps[,shared_features], indels[,shared_features])
    shared_features <- intersect(colnames(snps_TSS40), colnames(indels_TSS40))
    muts_TSS40 <- rbind(snps_TSS40[,shared_features], indels_TSS40[,shared_features]) # Should be same as TSS40, but re-running for TES annotations just in case.
    shared_features <- intersect(colnames(snps_TSS), colnames(indels_TSS))
    muts_TSS <- rbind(snps_TSS[,shared_features], indels_TSS[,shared_features])
    shared_features <- intersect(colnames(snps_TES), colnames(indels_TES)) # Should be same as TSS, but re-running for TES annotations just in case.
    muts_TES <- rbind(snps_TES[,shared_features], indels_TES[,shared_features])
    
    wgsa_data[["snps"]] <- snps
    wgsa_data[["snps_TSS40"]] <- snps_TSS40
    wgsa_data[["snps_TSS"]] <- snps_TSS
    wgsa_data[["snps_TES"]] <- snps_TES
    wgsa_data[["indels"]] <- indels
    wgsa_data[["indels_TSS40"]] <- indels_TSS40
    wgsa_data[["indels_TSS"]] <- indels_TSS
    wgsa_data[["indels_TES"]] <- indels_TES
    wgsa_data[["muts"]] <- muts
    wgsa_data[["muts_TSS40"]] <- muts_TSS40
    wgsa_data[["muts_TSS"]] <- muts_TSS
    wgsa_data[["muts_TES"]] <- muts_TES
    return(wgsa_data)
}

# Load unprocessed WGS dataset, liftover to hg19 genomic coordinates if on hg38, apply mappability + LCR filter, and write processed dataset to file.
process_wgs_dataset <- function(dat_name="dat", exclude_samples=c(), input_file, output_file=NULL, data_already_annotated=FALSE, from="hg38", to="hg38", DV_threshold=0, unpack_variant_info=TRUE, remove_low_qual_variants=TRUE, annotations=c("refGene", "gnomad_genome", "exac03", "caddgt10"), qc_plots=TRUE, apply_region_filters=FALSE, skip_DV=FALSE, confirm_refseq=FALSE, dat_sep=NULL, skip=0, header=TRUE) { #annotations=c("ucsc_refgene", "exac03")
    cat(paste0("\nLoading WGS dataset from input file ",input_file,"..."))
    return_env <- new.env()
    if (is.null(dat_sep) && grepl("csv", input_file)) { dat_sep = "," } else { dat_sep = "\t" }
    dat <- read.csv(input_file, sep=dat_sep, skip=skip, header=header)
    if(data_already_annotated) {
        cat("Done.\nData already formatted and annotated.\n")
    } else {
        cat("Done.\nFormatting and standardizing data...")
        # Format sample column and standardize column names.
        sample_col <- which(colnames(dat) == "IID")[1]
        if(is.na(sample_col)) { sample <- unlist(lapply(strsplit(paste0(dat$proband), "\\("), function(x) x[[1]][1])); dat <- cbind(dat, sample)
        } else { colnames(dat)[sample_col] <- "sample" }
        # Remove samples designated for exclusion.
        dat <- dat[which(!(dat$sample %in% exclude_samples)),]
        # Standardize column names for easier use.
        dat <- standardize_colnames(dat, re_order=TRUE)
        # Add snv_indel annotation for quick lookup of whether a variant is an snv or an indel
        snv_indel <- rep("snv", nrow(dat)); snv_indel[is_indel(dat$Ref, dat$Alt)] <- "indel"
        dat <- cbind(dat, snv_indel)
        dat_colnames <- colnames(dat)
        # Remove "chromosomes" other than chromosomes 1 through 22 and X and Y; ie: remove black-listed GL* chromosomes.
        dat <- dat[which(dat$Chrom %in% c(paste(1:22),"X","Y", paste0("chr",c(paste(1:22),"X","Y")))),]
        cat("Done.\n")
    
        # Unpack GATK variant info, based on the FORMAT column
        if(unpack_variant_info) {
            if(sum(c("FORMAT", "proband", "parents") %in% dat_colnames) == 3) {  
                variant_info_colnames <- apply(expand.grid(sort(unique(unlist(strsplit(paste0(dat$FORMAT), ":")))), c("proband", "father", "mother")), 1, function(x) paste(x[c(2,1)], collapse="_"))
                variant_info_colnames_length <- length(variant_info_colnames)
                variant_info_blank <- rep(".", variant_info_colnames_length); names(variant_info_blank) <- variant_info_colnames
                num_variants = nrow(dat)
                dat_variant_info <- data.frame(t(data.frame(sapply(1:nrow(dat), function(i) {
                    cat(paste0("Unpacking variant info...[",i," / ",num_variants,"]\n"))
                    variant_info_keys <- strsplit(paste0(dat$FORMAT[i]), ":")[[1]]
                    variant_info <- variant_info_blank
                    variant_info[paste0("proband_",variant_info_keys)] <- strsplit(gsub("^.*\\(|\\).*$", "", paste0(dat$proband[i])), ":")[[1]]
                    parents_variant_info <- strsplit(paste0(dat$parents[i]), "\\),")[[1]]
                    variant_info[paste0("father_",variant_info_keys)] <- strsplit(gsub("^.*\\(|\\).*$", "", parents_variant_info[1]), ":")[[1]]
                    variant_info[paste0("mother_",variant_info_keys)] <- strsplit(gsub("^.*\\(|\\).*$", "", parents_variant_info[2]), ":")[[1]]
                    return(variant_info)
                }))))
                dat <- cbind(dat, dat_variant_info)
                cat("Done.\n")
                dat_colnames <- colnames(dat)
            } 
        }
        
        # Use AnnovarR to annotate the requested annotations, if available.
        if(!is.null(annotations)) { 
            cat(paste0("Annotating variants with AnnovarR...[",from,": { ",paste(annotations,collapse=", ")," }]\n"))
            dat <- run_annovar(dat, annotations, buildver=from)
            cat("Done.\n")
            dat_colnames <- colnames(dat)
        }
    }
    
    # Remove low quality variants, double-checking results of the first-pass denovos script.
    if(remove_low_qual_variants) {
        snv_indel <- dat$snv_indel
        num_variants = nrow(dat)
        ### Reduce false positives in cases. ###
        cat("Removing low quality variants...[Reducing false positives in cases]")
        cat("\n",nrow(dat))
        # Remove genotype quality proband_PL(0/0) < 70
        dat <- dat[which(as.numeric(unlist(lapply(strsplit(paste0(dat$proband_PL), ","), function(x) x[1]))) >= 70),]
        cat(" ->",nrow(dat))
        # Remove highly recurrent within cohort AC >= 3 (in cohorts with fewer than 1000 cases)
        if(length(unique(dat$sample)) < 1000) { dat_strings <- genomic_dat_to_strings(dat); dat_strings_counts <- table(dat_strings); dat <- dat[which(dat_strings_counts[dat_strings] < 3),] }
        cat(" ->",nrow(dat))
        # Remove variants that do not pass GATK best practice: strand bias FS > 25, quality over depth QD < 2 (SNV) or QD < 1 (indel), Read position bias ReadPosRankSum < -3 (indels)
        dat <- dat[which(as.numeric(paste0(dat$FS)) <= 25 & ((snv_indel == "snv" & as.numeric(paste0(dat$QD)) >= 2) | (snv_indel == "indel" & as.numeric(paste0(dat$QD)) >= 1 & as.numeric(paste0(dat$ReadPosRankSum)) >= -3))),]
        cat(" ->",nrow(dat))
        # Remove number of reads supporting alternative allele proband_AD < 6 OR proband's alternative allele fraction < 20%
        dat <- dat[which(unlist(lapply(strsplit(paste0(dat$proband_AD), ","), function(x) { x <- as.numeric(x); AAD <- max(x[-c(1)]); return(AAD >= 6 && AAD/sum(x) >= 0.2) }))),]
        cat(" ->",nrow(dat))
        # Remove black-listed MUC or HLA genes
        if("Gene.refGene" %in% dat_colnames) { dat <- dat[which(!grepl("(^|/|;)(MUC|HLA)", dat$Gene.refGene)),]
        } else if("Gene" %in% dat_colnames) { dat <- dat[which(!grepl("(^|/|;)(MUC|HLA)", dat$Gene)),] }
        cat(" ->",nrow(dat))
        cat(paste0("\n",num_variants-nrow(dat)," variants removed.\n"))
        num_variants = nrow(dat)
        ### Reduce false negatives in parents. ###
        cat("Removing low quality variants...[Reducing false negatives in parents]")
        cat("\n",nrow(dat))
        # Remove genomic population allele frequency gnomAD_genome_ALL > 0.1%
        if("gnomAD_genome_ALL" %in% dat_colnames) { dat <- dat[which(is.na(dat$gnomAD_genome_ALL) | dat$gnomAD_genome_ALL == "." | as.numeric(paste0(dat$gnomAD_genome_ALL)) <= 0.0001),] }
        cat(" ->",nrow(dat))
        # Remove coding population allele frequency ExAC_ALL > 0.1%
        if("ExAC_ALL" %in% dat_colnames) { dat <- dat[which(is.na(dat$ExAC_ALL) | dat$ExAC_ALL == "." | as.numeric(paste0(dat$ExAC_ALL)) <= 0.001),] }
        cat(" ->",nrow(dat))
        # Remove depth of coverage in parents father_DP < 10 OR mother_DP < 10
        dat <- dat[which(as.numeric(paste0(dat$father_DP)) >= 10 & as.numeric(paste0(dat$mother_DP)) >= 10),]
        cat(" ->",nrow(dat))
        # Remove genotype quality in parents father_GQ < 30 OR mother_GQ < 30
        dat <- dat[which(as.numeric(paste0(dat$father_GQ)) >= 30 & as.numeric(paste0(dat$mother_GQ)) >= 30),]
        cat(" ->",nrow(dat))
        # Remove alternative allele fraction in parents > 3.5%
        dat <- dat[which(unlist(lapply(strsplit(paste0(dat$father_AD), ","), function(x) { x <- as.numeric(x); AAD <- max(x[-c(1)]); return(AAD/sum(x) <= 0.035) })) 
                         & unlist(lapply(strsplit(paste0(dat$mother_AD), ","), function(x) { x <- as.numeric(x); AAD <- max(x[-c(1)]); return(AAD/sum(x) <= 0.035) }))),]
        cat(" ->",nrow(dat))
        cat(paste0("\n",num_variants-nrow(dat)," variants removed.\nDone.\n"))
        num_variants = nrow(dat)
    }
    
    # Perform liftover (either hg19->hg38, or hg38->hg19), if "from" and "to" assemblies are specified.
    if(!is.null(from) && !is.null(to) && from!=to) { 
        cat(paste0("Performing ",from,"->",to," liftover with UCSC tool..."))
        dat <- liftover(dat, from=from, to=to, chr_colname="Chrom", start_colname="Position", ref_colname="Ref", alt_colname="Alt", mismatches_pause=FALSE, confirm_refseq=confirm_refseq) 
        cat("Done.\n")
    }
    
    # Create GRanges object from dat
    dat <- data.frame(dat)
    dat_granges <- to_genomic_regions(dat, chr_colname="Chrom", start_colname="Position", end_colname="Position", labels=dat$snv_indel)
    dat <- unfactorize(dat)
    dat_colnames <- colnames(dat)
    DAT_SAMPLE_COUNT = length(unique(dat$sample))
    
    # Apply all region filters (300bp mappability, segdups, low complexity regions, DAC and Duke Blacklists)
    if(apply_region_filters) {
        cat("Applying mappability and region filters (300bp mappability, segdups, low complexity regions, DAC and Duke Blacklists)...")
        all_filters_granges <- get_global(paste0("all_filters_granges_",to))
        if(is.null(all_filters_granges)) {
            cat(paste0("ERROR: Did not find all_filters_granges_",to," in workspace!\n"))
        } else {
            hits <- unique(queryHits(findOverlaps(dat_granges, all_filters_granges)))
            dat <- dat[-c(hits),]; dat_granges <- dat_granges[-c(hits)]
            cat("Done.\n")
        }
    }
    
    # Draw QC plots
    if(qc_plots) {
        # Effect of DV Quality threshold on variants/sample
        if(!skip_DV) {
            dat_include <- dat$DeepvarFilter == "PASS" & dat$FamMembers.DeepvarFilterPASS == "."
            dat <- dat[dat_include,]; dat_granges <- dat_granges[dat_include,]
            dat$DeepvarQual <- as.numeric(paste0(dat$DeepvarQual))
            filename = output_path(paste0(dat_name,"_effect_of_dv_quality_threshold.pdf"))
            pdf(file=filename)
            dv_cutoffs <- seq(0, 50, by=1)
            plot(dv_cutoffs, sapply(dv_cutoffs, function(x) sum(dat$DeepvarQual >= x))/DAT_SAMPLE_COUNT, main=paste0("Effect of DV quality threshold on ",dat_name," variants/sample"), xlab="DV_Qual", ylab="variants/sample", type="l", lwd=2, cex.axis=1.4, cex.lab=1.4, cex.main=1.3)
            dev.off()
            pdf_to_png(filename)
            # Filter by DeepVariant, with DV_Quality >= DV_threshold (default 30)
            dat_include <- dat$DeepvarQual >= DV_threshold
            dat <- dat[dat_include,]; dat_granges <- dat_granges[dat_include,]
        }
    
        # Histograms of the number of de novo SNVs and indels per sample
        incomplete_processing = c()
        if(!("gnomAD_genome_ALL" %in% dat_colnames)) { incomplete_processing = c("gnomAD") }
        if(!("ExAC_ALL" %in% dat_colnames)) { incomplete_processing = c(incomplete_processing, "ExAC") }
        if(sum(c("Gene.refGene","Gene") %in% dat_colnames)==0) { incomplete_processing = c(incomplete_processing, "MUC/HLA genes") }
        filename = output_path(paste0(dat_name,"_denovos_snv_count_density.pdf"))
        pdf(file=filename)
        p1 <- hist(table(dat[dat$snv_indel == "snv",]$sample), plot=FALSE, breaks=30)
        plot(p1, main=paste0(dat_name," denovo SNV count density"), xlab="# variants", col="grey", xaxt="n") 
        axis(1, at = seq(0, (round(max(p1$mids)/10)+1)*10, by=10))
        info_text = paste0("mean = ",round(mean(table(dat[dat$snv_indel == "snv",]$sample)), 2),", median = ",round(median(table(dat[dat$snv_indel == "snv",]$sample)), 2))
        if(length(incomplete_processing) > 0) { info_text = paste0(info_text,", variants not yet filtered with ",paste(incomplete_processing, collapse=" or ")) }
        mtext(info_text)
        dev.off()
        pdf_to_png(filename)
        filename = output_path(paste0(dat_name,"_denovos_indel_count_density.pdf"))
        pdf(file=filename)
        p2 <- hist(table(dat[dat$snv_indel == "indel",]$sample), plot=FALSE, breaks=30)
        plot(p2, main=paste0(dat_name," denovo indel count density"), xlab="# variants", col="grey", xaxt="n")
        axis(1, at = seq(0, (round(max(p2$mids)/10)+1)*10, by=10))
        info_text = paste0("mean = ",round(mean(table(dat[dat$snv_indel == "indel",]$sample)), 2),", median = ",round(median(table(dat[dat$snv_indel == "indel",]$sample)), 2))
        if(length(incomplete_processing) > 0) { info_text = paste0(info_text,", variants not yet filtered with ",paste(incomplete_processing, collapse=" or ")) }
        mtext(info_text)
        dev.off()
        pdf_to_png(filename)
        
        # Average (per sample) de novos per chromosome
        chromosomes <- c(1:22,"X","Y")
        dat_chrom_counts <- table(dat$Chrom)[paste0("chr",chromosomes)] / DAT_SAMPLE_COUNT
        filename = output_path(paste0(dat_name,"_denovos_per_chrom.pdf"))
        pdf(file=filename)
        barplot(dat_chrom_counts, names.arg=chromosomes, main=paste0("Average ",dat_name," de novos per chromosome"), xlab="chromosome", ylab="variants/sample", cex.main=1.3, cex.lab=1.4)
        dev.off()
        pdf_to_png(filename)
    }
    
    # Write processed dataset to file.
    if(!is.null(output_file)) {
        cat(paste0("Writing processed dataset to output file ",output_file,"..."))
        if(grepl("\\.csv$", output_file)) { write.csv(dat, file=output_file, quote=FALSE, row.names=FALSE)
        } else { write.table(dat, file=output_file, quote=FALSE, row.names=FALSE, sep="\t") }
        cat("Done.\n")
    }
    
    #cat("\nCommands to run WGSA (on Poweredge) with 32 GB memory and 20 CPU threads:\nEdit dat_WGSA065_config.txt in WGSA directory (/home/local/ARCS/yshen/software/WGSA) to point to correct SNP and indel files for this dataset.\n:> java WGSA065 dat_WGSA065_config.txt 32 20\n:> bash dat_WGSA065_config.txt.sh\n")
    
    return_env[["dat"]] <- dat
    return_env[["granges"]] <- dat_granges
    return_env[["sample_count"]] <- DAT_SAMPLE_COUNT
    return_env[["output_file"]] <- output_file
    return(return_env)
}

##################################################################################################################################################################
# Re-called SSC and CDH 
##################################################################################################################################################################
# SSC
# Create SSC controls ped file, using one unique unaffected sibling per family to avoid multi-counting variants
ssc_id_mappings <- read.csv(data_path("WGS/SSC/SSC_postQC_IDMap.txt"), sep="\t", header=FALSE)
ssc_controls_ped <- data.frame(t(data.frame(sapply(sort(unique(gsub("\\..*$", "", ssc_id_mappings$V2))), function(x) return(c(x, paste0(x,".s1"), paste0(x,".fa"), paste0(x,".mo")))))))
colnames(ssc_controls_ped) <- c("family_id", "sample_id", "father_id", "mother_id")
write.table(ssc_controls_ped, file=output_path("SSC_controls.ped"), sep="\t", row.names=FALSE)
# Raw denovos -> DeepVariant input
raw_denovos_file = data_path("WGS/re_called/vcf_ssc_all_filtered_2020-01-18.csv")
denovos_raw <- read.csv(raw_denovos_file, sep=",")
IID <- unlist(lapply(strsplit(paste0(denovos_raw$proband), "\\("), function(x) x[[1]][1]))
denovos_raw <- cbind(denovos_raw[,c(1:8)], IID, denovos_raw[,-c(1:8)])
denovos_raw <- denovos_raw[!(denovos_raw$ALT == "*"),] #| paste0(denovos_raw$IID) %in% c("1-15237", "1-09567")
sort(table(denovos_raw$IID))
denovos_raw <- denovos_raw[!(paste0(denovos_raw$IID) %in% c("11022.s1", "13976.s1", "14655.s1", "12970.s1", "13949.s1")),]
sort(table(denovos_raw$IID))
write.table(denovos_raw, data_path("WGS/re_called/vcf_ssc_all_filtered_2020-01-18.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

# CDH1
# Raw denovos -> DeepVariant input
raw_denovos_file = data_path("WGS/re_called/vcf_cdh1_filtered_2019-09-10.csv")
denovos_raw <- read.csv(raw_denovos_file, sep=",")
IID <- unlist(lapply(strsplit(paste0(denovos_raw$proband), "\\("), function(x) x[[1]][1]))
denovos_raw <- cbind(denovos_raw[,c(1:8)], IID, denovos_raw[,-c(1:8)])
denovos_raw <- denovos_raw[!(denovos_raw$ALT == "*"),] #| paste0(denovos_raw$IID) %in% c("1-15237", "1-09567")
denovos_raw[1,]
denovos_raw <- liftover(denovos_raw, from="hg19", to="hg38", chr_colname="CHROM", start_colname="POS", end_colname=NULL, ref_colname="REF", alt_colname="ALT", work_folder=OUTPUT_FOLDER, program_path="/home/local/ARCS/ak3792/programs/ucsc/liftOver", confirm_refseq=FALSE, mismatches_pause=FALSE)
denovos_raw$CHROM <- paste0("chr",denovos_raw$CHROM)
sort(table(denovos_raw$IID))
write.table(denovos_raw, data_path("WGS/re_called/vcf_cdh1_filtered_2019-09-10_hg38.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

# CDH2
# Raw denovos -> DeepVariant input
raw_denovos_file = data_path("WGS/re_called/vcf_cdh2_filtered_2019-08-13.csv")
denovos_raw <- read.csv(raw_denovos_file, sep=",")
IID <- unlist(lapply(strsplit(paste0(denovos_raw$proband), "\\("), function(x) x[[1]][1]))
denovos_raw <- cbind(denovos_raw[,c(1:8)], IID, denovos_raw[,-c(1:8)])
denovos_raw <- denovos_raw[!(denovos_raw$ALT == "*"),] #| paste0(denovos_raw$IID) %in% c("1-15237", "1-09567")
denovos_raw[1,]
sample_info <- read.csv(data_path("CDH2_Apr2018_sampleinfo.csv"))
a <- read.table(data_path("WGS/re_called/WGS_cdh2_Broad.ped"))
b <- unique(denovos_raw$IID[!(denovos_raw$IID %in% c(paste0(sample_info$IID), paste0(sample_info$DAD), paste0(sample_info$MOM)))])
a[a$V2 %in% b,]
haha <- read.table(data_path("WGS/re_called/CDH2_sample.txt"), header=TRUE)

sort(table(denovos_raw$IID))
denovos_raw <- denovos_raw[!(paste0(denovos_raw$IID) %in% c("151760", "CDH05-0025", "CDH06-0012")),]
sort(table(denovos_raw$IID))
write.table(denovos_raw, data_path("WGS/re_called/vcf_cdh2_filtered_2019-08-13.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

# CDH3
# Create CDH3 ped file
cdh3_pedigree <- read.csv(data_path("WGS/CDH/batch3/WGS_PCRFree.ped"), sep="\t", header=FALSE)
colnames(cdh3_pedigree) <- c("family_id", "sample_id", "father_id", "mother_id", "gender", "phenotype")
cdh3_pedigree$gender[cdh3_pedigree$gender == 1] <- "M"; cdh3_pedigree$gender[cdh3_pedigree$gender == 2] <- "F"
cdh3_pedigree$phenotype[cdh3_pedigree$phenotype == 1] <- "healthy"; cdh3_pedigree$phenotype[cdh3_pedigree$phenotype == 2] <- "CDH"
for(i in 1:ncol(cdh3_pedigree)) { cdh3_pedigree[,i] <- gsub("\\s+", "", cdh3_pedigree[,i]) }
father_info <- cdh3_pedigree[cdh3_pedigree$father_id == "0" & cdh3_pedigree$mother_id == "0" & grepl("M",cdh3_pedigree$gender), c("family_id", "sample_id", "phenotype")]; colnames(father_info)[-c(1:2)] <- paste0("father_", colnames(father_info)[-c(1:2)])
mother_info <- cdh3_pedigree[cdh3_pedigree$father_id == "0" & cdh3_pedigree$mother_id == "0" & grepl("F",cdh3_pedigree$gender), c("family_id", "sample_id", "phenotype")]; colnames(mother_info)[-c(1:2)] <- paste0("mother_", colnames(mother_info)[-c(1:2)])
cdh3_ped <- merge(cdh3_pedigree[cdh3_pedigree$father_id != "0" & cdh3_pedigree$mother_id != "0",], father_info, by.x=c("family_id", "father_id"), by.y=c("family_id", "sample_id"), all.x=TRUE)
cdh3_ped <- merge(cdh3_ped, mother_info, by.x=c("family_id", "mother_id"), by.y=c("family_id", "sample_id"), all.x=TRUE)
cdh3_ped <- cdh3_ped[,c("family_id", "sample_id", "father_id", "mother_id", "gender", "phenotype", "father_phenotype", "mother_phenotype")]
write.table(cdh3_ped, file=data_path("WGS/CDH/batch3/CDH3.ped"), sep="\t", row.names=FALSE)

cdh3_raw_denovos_file = data_path("WGS/CDH/CDH3_vcf_filtered_2019-05-17_hg38_norm_DeepVariant.txt")
cdh3_raw <- read.csv(cdh3_raw_denovos_file, sep="\t")
cdh3_raw[grepl("CDH2-75", cdh3_raw$proband),][1:10,]
cdh3_raw[grepl("CDH640", cdh3_raw$proband),][1:10,]
dat <- read.csv(data_path("WGS/re_called/vcf_cdh3_filtered_2019-05-17_hg38_norm_DeepVariant.txt"), sep="\t")
sort(table(dat$IID))

# CHD3 (TOPMed)
cram_list <- read.csv(data_path("WGS/CHD/topmed810_cram_list.txt"), sep=" ", header=FALSE) 
nrow(cram_list)
crams <- cram_list[,2]
crams_counts <- table(gsub("-0[12]$", "", crams))
length(unique(topmed810_raw$IID))
sum(crams_counts == 3)
length(unique(topmed810$sample))
topmed810_ped <- read.csv(data_path("WGS/CHD/topmed810.ped"), sep="\t", header=FALSE) 
chd3 <- chd3[!(chd3$sample %in% c("1-15237", "1-09567")),]
# Save GATK-only TOPMed result
topmed810_raw <- read.csv(data_path("WGS/CHD/vcf_topmed810_filtered_2019-12-08_hg38_norm.txt"), sep="\t")
#sort(table(topmed810_raw$IID))
denovos_processed <- process_wgs_dataset("TOPMed", input_file=data_path("WGS/CHD/vcf_topmed810_filtered_2019-12-08_hg38_norm.txt"), output_file=data_path("WGS/re_called/chd3_gatk_only_final.tsv"), skip_DV=TRUE, qc_plots=FALSE)
denovos_processed <- process_wgs_dataset("TOPMed_filtered", input_file=data_path("WGS/re_called/chd3_gatk_only_final.tsv"), data_already_annotated=TRUE, remove_low_qual_variants=FALSE, apply_region_filters=TRUE)
chd3 <- denovos_processed$dat
chd3_granges <- denovos_processed$granges
CHD3_SAMPLE_COUNT <- denovos_processed$sample_count
TOPMED_SAMPLE_COUNT <- CHD3_SAMPLE_COUNT

# CHD4 (TOPMed second batch!)
# Raw denovos -> DeepVariant input
denovos_raw <- read.csv(data_path("WGS/re_called/vcf_topmed259_filtered_2020-05-07.csv"))
IID <- unlist(lapply(strsplit(paste0(denovos_raw$proband), "\\("), function(x) x[[1]][1]))
sort(table(IID))
denovos_raw <- cbind(denovos_raw[,c(1:8)], IID, denovos_raw[,-c(1:8)])
denovos_raw <- denovos_raw[!(denovos_raw$ALT == "*" | paste0(denovos_raw$IID) %in% c("1-13142", "1-07385")),] #| paste0(denovos_raw$IID) %in% c("1-15237", "1-09567")
write.table(denovos_raw, data_path("WGS/re_called/vcf_topmed259_filtered_2020-05-07.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

denovos_processed <- process_wgs_dataset("TOPMed_259", input_file=data_path("WGS/re_called/vcf_topmed259_filtered_2020-05-07.tsv"), output_file=data_path("WGS/re_called/topmed259_hg38_gatk_only_before_region_filters.tsv"), skip_DV=TRUE, qc_plots=FALSE)
denovos_processed <- process_wgs_dataset("TOPMed_259_filtered", input_file=data_path("WGS/re_called/topmed259_hg38_gatk_only_before_region_filters.tsv"), output_file=data_path("WGS/re_called/topmed259_hg38_gatk_only.tsv"), skip_DV=TRUE, qc_plots=FALSE, data_already_annotated=TRUE, remove_low_qual_variants=FALSE, apply_region_filters=TRUE)
chd3 <- denovos_processed$dat

##################################################################################################################################################################
# Analysis of DeepVariant output, and processing to final de novo variant files.
##################################################################################################################################################################
re_called <- combine_datasets(combine_datasets(combine_datasets(combine_datasets(
        read.csv(data_path("WGS/re_called/vcf_ssc_all_filtered_2020-01-18_hg38_norm_DeepVariant.txt"), sep="\t"),
        read.csv(data_path("WGS/re_called/vcf_cdh1_filtered_2019-09-10_hg38_norm_DeepVariant.txt"), sep="\t")),
        read.csv(data_path("WGS/re_called/vcf_cdh2_filtered_2019-08-13_hg38_norm_DeepVariant.txt"), sep="\t")),
        read.csv(data_path("WGS/re_called/vcf_cdh3_filtered_2019-05-17_hg38_norm_DeepVariant.txt"), sep="\t")),
        read.csv(data_path("WGS/re_called/vcf_topmed810_filtered_2020-01-07_hg38_norm_DeepVariant.txt"), sep="\t"))
re_called <- re_called[which(re_called$DeepvarFilter == "PASS" & re_called$FamMembers.DeepvarFilterPASS == "."),]
colnames(re_called)[which(colnames(re_called) == "IID")[1]] <- "sample"
snv_indel <- rep("snv", nrow(re_called)); snv_indel[is_indel(re_called$REF, re_called$ALT)] <- "indel"
re_called <- cbind(re_called, snv_indel)
re_called <- standardize_colnames(re_called, re_order=TRUE, remove_chr_prefix=TRUE)
write.table(re_called[,c("Chrom", "Position", "Ref", "Alt", "sample", "snv_indel")], file=data_path("WGS/re_called/re_called.tsv"), row.names=FALSE, quote=FALSE)
re_called_hg19 <- liftover(re_called, from="hg38", to="hg19", chr_colname="Chrom", start_colname="Position", end_colname=NULL, ref_colname="Ref", alt_colname="Alt", confirm_refseq=FALSE, mismatches_pause=FALSE)#[,c("Chrom","Position","Ref","Alt","sample","snv_indel","Chrom_hg38","Position_hg38")]
write.table(re_called_hg19[,c("Chrom", "Position", "Position", "Ref", "Alt", "sample", "snv_indel", "Chrom_hg38","Position_hg38")], file=data_path("WGS/re_called/re_called_hg19.tsv"), row.names=FALSE, quote=FALSE)
#perl annotate_variation.pl '/home/local/ARCS/ak3792/Documents/Research/data/WGS/re_called/re_called_hg19.tsv' humandb/ -filter -dbtype gnomad_genome -buildver hg19 -out ex1 -otherinfo
#perl table_annovar.pl '/home/local/ARCS/ak3792/Documents/Research/data/WGS/re_called/re_called_hg19.tsv' humandb/ -buildver hg19 -out myanno -remove -protocol refGene,gnomad_genome,exac03,caddgt10 -operation g,f,f,f -nastring . -csvout -polish


#cat("\nCommands to run WGSA (on Poweredge) with 32 GB memory and 40 CPU threads:\nEdit dat_WGSA065_config.txt in WGSA directory (/home/local/ARCS/yshen/software/WGSA) to point to correct SNP and indel files for this dataset.\n:> java WGSA065 re_called_WGSA065_config.txt 32 40\n:> bash re_called_WGSA065_config.txt.sh\n")

# SSC
#sort(table(read.csv(data_path("WGS/re_called/vcf_ssc_all_filtered_2020-01-18_hg38_norm_DeepVariant.txt"), sep="\t")$IID))
#denovos_processed <- process_wgs_dataset("SSC", input_file=data_path("WGS/re_called/vcf_ssc_all_filtered_2020-01-18_hg38_norm_DeepVariant.txt"), output_file=data_path("WGS/re_called/ssc_final.tsv"))
denovos_processed <- process_wgs_dataset("SSC_filtered", input_file=data_path("WGS/re_called/ssc_final.tsv"), data_already_annotated=TRUE, remove_low_qual_variants=FALSE, apply_region_filters=TRUE, DV_threshold=0)
ssc <- denovos_processed$dat
ssc_granges <- denovos_processed$granges
SSC_SAMPLE_COUNT <- denovos_processed$sample_count

# CDH1
#sort(table(read.csv(data_path("WGS/re_called/vcf_cdh1_filtered_2019-09-10_hg38_norm_DeepVariant.txt"), sep="\t")$IID))
#denovos_processed <- process_wgs_dataset("CDH1", exclude_samples=c(), input_file=data_path("WGS/re_called/vcf_cdh1_filtered_2019-09-10_hg38_norm_DeepVariant.txt"), output_file=data_path("WGS/re_called/cdh1_final.tsv"))
denovos_processed <- process_wgs_dataset("CDH1_filtered", input_file=data_path("WGS/re_called/cdh1_final.tsv"), data_already_annotated=TRUE, remove_low_qual_variants=FALSE, apply_region_filters=TRUE, DV_threshold=35)
cdh1 <- denovos_processed$dat
cdh1_granges <- denovos_processed$granges
CDH1_SAMPLE_COUNT <- denovos_processed$sample_count

# CDH2
#sort(table(read.csv(data_path("WGS/re_called/vcf_cdh2_filtered_2019-08-13_hg38_norm_DeepVariant.txt"), sep="\t")$IID))
#denovos_processed <- process_wgs_dataset("CDH2", exclude_samples=c(), input_file=data_path("WGS/re_called/vcf_cdh2_filtered_2019-08-13_hg38_norm_DeepVariant.txt"), output_file=data_path("WGS/re_called/cdh2_final.tsv"))
denovos_processed <- process_wgs_dataset("CDH2_filtered", input_file=data_path("WGS/re_called/cdh2_final.tsv"), data_already_annotated=TRUE, remove_low_qual_variants=FALSE, apply_region_filters=TRUE, DV_threshold=35)
cdh2 <- denovos_processed$dat
cdh2_granges <- denovos_processed$granges
CDH2_SAMPLE_COUNT <- denovos_processed$sample_count

# CDH3
#sort(table(read.csv(data_path("WGS/re_called/vcf_cdh3_filtered_2019-05-17_hg38_norm_DeepVariant.txt"), sep="\t")$IID))
#denovos_processed <- process_wgs_dataset("CDH3", exclude_samples=c("CDH472", "CDH1162","CDH1560","CDH1206"), input_file=data_path("WGS/re_called/vcf_cdh3_filtered_2019-05-17_hg38_norm_DeepVariant.txt"), output_file=data_path("WGS/re_called/cdh3_final.tsv"))
denovos_processed <- process_wgs_dataset("CDH3_filtered", input_file=data_path("WGS/re_called/cdh3_final.tsv"), data_already_annotated=TRUE, remove_low_qual_variants=FALSE, apply_region_filters=TRUE, DV_threshold=35)
cdh3 <- denovos_processed$dat
cdh3_granges <- denovos_processed$granges
CDH3_SAMPLE_COUNT <- denovos_processed$sample_count
#annotations=c("caddgt10")

# All CDH
cdh <- combine_datasets(combine_datasets(cdh1, cdh2),cdh3); cdh_granges <- c(cdh1_granges, cdh2_granges, cdh3_granges); CDH_SAMPLE_COUNT = CDH1_SAMPLE_COUNT + CDH2_SAMPLE_COUNT + CDH3_SAMPLE_COUNT


# CHD3 (TOPMed)
#sort(table(read.csv(data_path("WGS/re_called/vcf_topmed810_filtered_2020-01-07_hg38_norm_DeepVariant.txt"), sep="\t")$IID))
#denovos_processed <- process_wgs_dataset("TOPMed", input_file=data_path("WGS/re_called/vcf_topmed810_filtered_2020-01-07_hg38_norm_DeepVariant.txt"), output_file=data_path("WGS/re_called/chd3_final.tsv"))
denovos_processed <- process_wgs_dataset("TOPMed_filtered", input_file=data_path("WGS/re_called/chd3_final.tsv"), data_already_annotated=TRUE, remove_low_qual_variants=FALSE, apply_region_filters=TRUE, DV_threshold=35)
chd3 <- denovos_processed$dat
chd3_granges <- denovos_processed$granges
CHD3_SAMPLE_COUNT <- denovos_processed$sample_count
TOPMED_SAMPLE_COUNT <- CHD3_SAMPLE_COUNT

# CHD4 (TOPMed second batch!)
sort(table(read.csv(data_path("WGS/re_called/vcf_topmed259_filtered_2020-05-20_hg38_norm_DeepVariant.txt"), sep="\t")$IID))
denovos_processed <- process_wgs_dataset("TOPMed_259", input_file=data_path("WGS/re_called/vcf_topmed259_filtered_2020-05-20_hg38_norm_DeepVariant.txt"), output_file=data_path("WGS/re_called/topmed259_hg38_dvqual0_before_region_filters.tsv"), remove_low_qual_variants=TRUE, apply_region_filters=FALSE, DV_threshold=0)
denovos_processed <- process_wgs_dataset("TOPMed_259_filtered", input_file=data_path("WGS/re_called/chd3_final.tsv"), output_file=data_path("WGS/re_called/topmed259_hg38_dvqual0_final.tsv"), data_already_annotated=TRUE, remove_low_qual_variants=FALSE, apply_region_filters=TRUE, DV_threshold=0)
chd4_dv <- denovos_processed$dat
sort(table(chd4_dv$sample))
sort(table(chd4_dv$sample[chd4_dv$snv_indel == "snv"]))
sort(table(chd4_dv$sample[chd4_dv$snv_indel == "indel"]))
chd4_granges <- denovos_processed$granges
CHD4_SAMPLE_COUNT <- denovos_processed$sample_count
# GATK only below:
denovos_processed <- process_wgs_dataset("TOPMed_259", input_file=data_path("WGS/re_called/vcf_topmed259_filtered_2020-05-07.tsv"), output_file=data_path("WGS/re_called/topmed259_hg38_gatk_only_before_region_filters.tsv"), skip_DV=TRUE, qc_plots=FALSE)
denovos_processed <- process_wgs_dataset("TOPMed_259_filtered", input_file=data_path("WGS/re_called/topmed259_hg38_gatk_only_before_region_filters.tsv"), output_file=data_path("WGS/re_called/topmed259_hg38_gatk_only.tsv"), skip_DV=TRUE, qc_plots=FALSE, data_already_annotated=TRUE, remove_low_qual_variants=FALSE, apply_region_filters=TRUE)
chd4_gatk <- denovos_processed$dat
sort(table(chd4_gatk$sample))
sort(table(chd4_gatk$sample[chd4_gatk$snv_indel == "snv"]))
sort(table(chd4_gatk$sample[chd4_gatk$snv_indel == "indel"]))

# Process GATK+FreeBayes files sent by Sarah
a <- load(data_path("WGS/re_called/topmed1069_gatk_freebayes_dnvs.rda"))
all_gatk[1,]
all_gatk_fb[1,]
denovos_raw <- all_gatk
# Raw denovos -> DeepVariant input
denovos_raw <- read.csv(data_path("WGS/re_called/vcf_topmed259_filtered_2020-05-07.csv"))
IID <- unlist(lapply(strsplit(paste0(denovos_raw$proband), "\\("), function(x) x[[1]][1]))
sort(table(IID))
denovos_raw <- cbind(denovos_raw[,c(1:8)], IID, denovos_raw[,-c(1:8)])
denovos_raw <- denovos_raw[!(denovos_raw$ALT == "*" | paste0(denovos_raw$IID) %in% c("1-13142", "1-07385")),] #| paste0(denovos_raw$IID) %in% c("1-15237", "1-09567")
write.table(denovos_raw, data_path("WGS/re_called/vcf_topmed259_filtered_2020-05-07.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
# Annotate with DV output
sort(table(read.csv(data_path("WGS/re_called/vcf_topmed1069_filtered_2020-05-27_hg38_norm_DeepVariant.txt"), sep="\t")$IID))
denovos_processed <- read.csv(data_path("WGS/re_called/vcf_topmed1069_filtered_2020-05-27_hg38_norm_DeepVariant.txt"), sep="\t")
denovos_raw <- read.csv(data_path("WGS/re_called/vcf_topmed1069_filtered_2020-05-22.tsv"), sep="\t")
all_gatk[1,]
fb_variants <- unique(all_gatk_fb$uniqueVarID)
a <- merge(cbind(rbind(all_gatk_fb[,colnames(all_gatk)],all_gatk),c(rep(TRUE,nrow(all_gatk_fb)),rep(FALSE,nrow(all_gatk)))), denovos_processed[,c("uniqueVarID","IID","DeepvarFilter","DeepvarQual","NearVars","FamMembers.DeepvarFilterPASS","FamMembers.DeepvarFilterRefCall","FamMembers.DeepvarQual")], all.x=TRUE)
colnames(a) <- c("uniqueVarID","Chrom","Position","End","Ref","Alt","AC","FS","MQ","QD","FreeBayes_PASS","sample","DeepvarFilter","DeepVariant_QUAL","NearVars","FamMembers.DeepvarFilterPASS","FamMembers.DeepvarFilterRefCall","FamMembers.DeepvarQual")   
a <- a[!duplicated(a$uniqueVarID),]
DeepVariant_PASS <- a$DeepvarFilter == "PASS" & a$FamMembers.DeepvarFilterPASS == "."
a <- cbind(a, DeepVariant_PASS)
a <- a[a$uniqueVarID %in% unique(all_gatk$uniqueVarID),c("Chrom","Position","Ref","Alt","sample","uniqueVarID","End","AC","FS","MQ","QD","FreeBayes_PASS","DeepVariant_PASS","DeepVariant_QUAL")]
a$DeepVariant_QUAL <- as.numeric(paste0(a$DeepVariant_QUAL)); a$DeepVariant_QUAL[is.na(a$DeepVariant_QUAL)] <- 0
a$DeepVariant_PASS[is.na(a$DeepVariant_PASS)] <- FALSE; a$FreeBayes_PASS[is.na(a$FreeBayes_PASS)] <- FALSE
nrow(a)
sum(a$DeepVariant_PASS)
sum(a$DeepVariant_QUAL > 0)
sum(a$FreeBayes_PASS)
write.table(a, data_path("WGS/re_called/topmed1069_HMS_MSSM_FreeBayes_DeepVariant.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

##################################################################################################################################################################
# WGS Batch comparisons
##################################################################################################################################################################
# gnomAD densities for Average de novos per chrom
filename = output_path("gnomad_densities.pdf")
pdf(file=filename)
ssc_gnomad_densities <- paste0(ssc$gnomAD_genome_ALL); ssc_gnomad_densities[is.na(ssc_gnomad_densities) | ssc_gnomad_densities == "."] <- "0"; ssc_gnomad_densities <- as.numeric(ssc_gnomad_densities)
plot(density(log10(ssc_gnomad_densities), to=log10(1e-3)), col="blue", lwd=2, main="gnomAD_ALL Density", xlab="log10(gnomAD_ALL)", xaxs="i", yaxs="i", cex.main=1.3, cex.lab=1.4)
cdh_gnomad_densities <- paste0(cdh$gnomAD_genome_ALL); cdh_gnomad_densities[is.na(cdh_gnomad_densities) | cdh_gnomad_densities == "."] <- "0"; cdh_gnomad_densities <- as.numeric(cdh_gnomad_densities)
lines(density(log10(cdh_gnomad_densities), to=log10(1e-3)), col="red", lwd=2)
#plot(density(log10(as.numeric(ssc$gnomAD_genome_ALL[!is.na(ssc$gnomAD_genome_ALL) & ssc$gnomAD_genome_ALL != "."])), to=log10(1e-3)), col="blue", lwd=2, main="gnomAD_ALL Density", xlab="log10(gnomAD_ALL)", xaxs="i", yaxs="i", cex.main=1.3, cex.lab=1.4)
#lines(density(log10(as.numeric(cdh$gnomAD_genome_ALL[!is.na(cdh$gnomAD_genome_ALL) & cdh$gnomAD_genome_ALL != "."])), to=log10(1e-3)), col="red", lwd=2)
legend("topright", legend=c(paste0("CDH (",sum(cdh_gnomad_densities>=1e-4)," / ",length(cdh_gnomad_densities)," >= 1e-4)"), paste0("SSC (",sum(ssc_gnomad_densities>=1e-4)," / ",length(ssc_gnomad_densities)," >= 1e-4)")), col=c("red", "blue"), pch=15, cex=1.2)
dev.off()
pdf_to_png(filename)

# AF densities
ssc_an <- sapply(paste0(ssc$AN), function(x) max(as.numeric(strsplit(x,",")[[1]])))
ssc_af <- sapply(paste0(ssc$AF), function(x) max(as.numeric(strsplit(x,",")[[1]])))
ssc_ac <- ssc_an * ssc_af
cdh1_an <- sapply(paste0(cdh1$AN), function(x) max(as.numeric(strsplit(x,",")[[1]])))
cdh1_af <- sapply(paste0(cdh1$AF), function(x) max(as.numeric(strsplit(x,",")[[1]])))
cdh1_ac <- cdh1_an * cdh1_af
cdh2_an <- sapply(paste0(cdh2$AN), function(x) max(as.numeric(strsplit(x,",")[[1]])))
cdh2_af <- sapply(paste0(cdh2$AF), function(x) max(as.numeric(strsplit(x,",")[[1]])))
cdh2_ac <- cdh2_an * cdh2_af
cdh3_an <- sapply(paste0(cdh3$AN), function(x) max(as.numeric(strsplit(x,",")[[1]])))
cdh3_af <- sapply(paste0(cdh3$AF), function(x) max(as.numeric(strsplit(x,",")[[1]])))
cdh3_ac <- cdh3_an * cdh3_af

plot(density(log10(ssc_gnomad_densities), to=log10(1e-3)), col="blue", lwd=2, main="gnomAD_ALL Density", xlab="log10(gnomAD_ALL)", xaxs="i", yaxs="i", cex.main=1.3, cex.lab=1.4)

# Find GATA4/GATA6 variants
cdh[which(grepl("(^|/|;)(GATA4|GATA6)", cdh$Gene.refGene)),]

# Average de novos per chromosome, for different batches
chromosomes <- c(1:22,"X")
ssc_chrom_counts <- table(ssc$Chrom)[paste0("chr",chromosomes)] / SSC_SAMPLE_COUNT
cdh1_chrom_counts <- table(cdh1$Chrom)[paste0("chr",chromosomes)] / CDH1_SAMPLE_COUNT
cdh2_chrom_counts <- table(cdh2$Chrom)[paste0("chr",chromosomes)] / CDH2_SAMPLE_COUNT
cdh3_chrom_counts <- table(cdh3$Chrom)[paste0("chr",chromosomes)] / CDH3_SAMPLE_COUNT
chd3_chrom_counts <- table(chd3$Chrom)[paste0("chr",chromosomes)] / CHD3_SAMPLE_COUNT
filename = output_path("denovos_per_chrom.pdf")
pdf(file=filename, width=14)
barplot(rbind(ssc_chrom_counts,cdh1_chrom_counts,cdh2_chrom_counts,cdh3_chrom_counts,chd3_chrom_counts), names.arg=chromosomes, beside=TRUE, col=rainbow(5), legend.text=c("SSC", "CDH1", "CDH2", "CDH3","TOPMed"), main="Average de novos per chromosome", xlab="chromosome", ylab="variants/sample", cex.main=1.3, cex.lab=1.4)
dev.off()
pdf_to_png(filename)

# Effect of paternal age on variants/sample, for different batches
plot_variants_vs_paternal_age <- function(dats, paternal_ages, batches=NULL, id_match_prefixes=c("^"), cols=NULL, col_alpha=0.4, filename=NULL) {
    if(class(dats) != "list") { dats <- list(dats) }
    if(class(paternal_ages) != "list") { paternal_ages <- list(paternal_ages) }
    if(class(batches) != "list") { batches <- list(batches) }
    if(class(id_match_prefixes) != "list") { id_match_prefixes <- list(id_match_prefixes) }
    num_dats <- length(dats)
    all_dats_paternal_age <- lapply(1:num_dats, function(i) {
        dat <- dats[[i]]
        if(i > length(paternal_ages)) { paternal_age <- paternal_ages[[1]] } else { paternal_age <- paternal_ages[[i]] } 
        if(i > length(batches)) { batch <- batches[[1]] } else { batch = batches[[i]] } 
        if(i > length(id_match_prefixes)) { id_match_prefix <- id_match_prefixes[[1]] } else { id_match_prefix = id_match_prefixes[[i]] } 
        dat_samples <- dat$sample
        dat_snv_counts <- data.frame(table(dat$sample[dat$snv_indel=="snv"])); dat_indel_counts <- data.frame(table(dat$sample[dat$snv_indel=="indel"]))
        dat_variant_counts <- merge(dat_snv_counts, dat_indel_counts, by="Var1", all.x=TRUE, all.y=TRUE); dat_variant_counts[is.na(dat_variant_counts)] <- 0; colnames(dat_variant_counts) <- c("sample", "snv_count", "indel_count")
        if(is.null(batch)) { batch <- rep(paste0("dat",i),length(dat_samples)) 
        } else if(length(batch) < length(dat_samples)) { batch <- rep(batch,length(dat_samples)) }
        all_samples <- data.frame(dat_samples, batch)
        all_samples <- all_samples[!duplicated(dat_samples),]
        colnames(all_samples)[1] <- "sample"
        all_samples <- unfactorize(merge(all_samples, dat_variant_counts, by="sample"))
        dat_paternal_age_samples <- unfactorize(data.frame(t(sapply(paste0(paternal_age$Blinded.ID), function(id) {
            return(unlist(all_samples[grepl(paste0(paste0(id_match_prefix,id),collapse="|"), all_samples$sample),][1,]))
        }))))
        table(dat_paternal_age_samples$batch)
        dat_paternal_age <- cbind(paternal_age, dat_paternal_age_samples)
        dat_paternal_age <- dat_paternal_age[!is.na(dat_paternal_age$sample),]
        dat_paternal_age <- dat_paternal_age[!(is.na(dat_paternal_age$Paternal.Age.at.Proband.Birth) | dat_paternal_age$Paternal.Age.at.Proband.Birth == ""),]
        table(dat_paternal_age$batch)
        if(!is.numeric(dat_paternal_age$Paternal.Age.at.Proband.Birth)) {
            dat_paternal_age$Paternal.Age.at.Proband.Birth <- sapply(dat_paternal_age$Paternal.Age.at.Proband.Birth, function(x) { years <- as.numeric(gsub("(y[, ]+[0-9]+m)|(.*m[, ]+)|y","", x)); months <-as.numeric(gsub("(m[, ]+[0-9]+y)|(.*y[, ]+)|m","", x)); if(!is.na(months)) { years = years + (months/12) }; return(years) })
        }
        dat_paternal_age$snv_count <- as.numeric(paste0(dat_paternal_age$snv_count)); dat_paternal_age$indel_count <- as.numeric(paste0(dat_paternal_age$indel_count))
        return(dat_paternal_age)
    })
    # Draw combined plot
    combined_paternal_age <- all_dats_paternal_age[[1]]
    if(num_dats > 1) { for(i in 2:num_dats) { combined_paternal_age <- combine_datasets(combined_paternal_age, all_dats_paternal_age[[i]]) } }
    save_plot_to_file = !is.null(filename)
    if(save_plot_to_file) { pdf(file=filename) }
    batch <- as.factor(combined_paternal_age$batch)
    batches <- unique(batch)
    num_batches <- length(batches)
    if(is.null(cols)) { cols <- palette()[4:(3+num_batches)] }
    cols <- adjustcolor(cols, alpha.f=col_alpha)
    scatterplot(combined_paternal_age$Paternal.Age.at.Proband.Birth, combined_paternal_age$snv_count + combined_paternal_age$indel_count, legend.plot=FALSE, groups=batch, col=cols, pch=rep(1,num_batches), boxplots="xy", main="", xlab="paternal age (years)", ylab="variants/sample", cex.main=1.5, cex.axis=1.4, cex.lab=1.3, spread=FALSE, smooth=FALSE, legend.coords="topleft", legend.columns=2)
    #mtext(paste0("Pearson = ", round(cor(combined_paternal_age$Paternal.Age.at.Proband.Birth, combined_paternal_age$snv_count + combined_paternal_age$indel_count, method="pearson"),3)))
    if(save_plot_to_file) {
        dev.off(); pdf_to_png(filename) 
        cdh_batch_means <- rbind(sapply(batches, function(batch) get_global(paste0(toupper(batch),"_SAMPLE_COUNT")))
                                 ,sapply(batches, function(batch) mean((combined_paternal_age$Paternal.Age.at.Proband.Birth)[combined_paternal_age$batch == batch]))
                                 ,sapply(batches, function(batch) mean((combined_paternal_age$snv_count + combined_paternal_age$indel_count)[combined_paternal_age$batch == batch]))
                                 ,sapply(batches, function(batch) mean((combined_paternal_age$snv_count)[combined_paternal_age$batch == batch]))
                                 ,sapply(batches, function(batch) mean((combined_paternal_age$indel_count)[combined_paternal_age$batch == batch]))
        ); colnames(cdh_batch_means) <- batches; rownames(cdh_batch_means) <- c("num_samples", "paternal_age", "variants", "SNVs", "indels")
        write.csv(cdh_batch_means, file=gsub(".pdf","_batch_means.csv",filename))
    }
}
# Load SSC paternal age information.
ssc_parental_age <- read.csv(data_path("parental_age_at_proband_birth_2018_06_27.txt"), sep="\t")
ssc_id_mappings <- read.csv(data_path("WGS/SSC/SSC_postQC_IDMap.txt"), sep="\t", header=FALSE); ssc_id_mappings$V1 <- paste0(ssc_id_mappings$V1); ssc_id_mappings$V2 <- paste0(ssc_id_mappings$V2)
ssc_parental_age$Blinded.ID <- sapply(paste0(ssc_parental_age$Blinded.ID), function(x) ssc_id_mappings$V2[which(ssc_id_mappings$V1 == x)[1]])
ssc_parental_age <- ssc_parental_age[which(!is.na(ssc_parental_age$Blinded.ID)),]
# Load CHD parental age information.
chd_parental_age <- read.csv(data_path("WGS/CHD/all_chd_wgs_ages.csv"))
# Load CDH parental age information.
cdh_parental_age <- read.csv(data_path("CDH_info_20190624DHREAMSDB.csv")) #read.csv(data_path("WGS/CDH/cdh_paternal_ages.csv"))
colnames(cdh_parental_age) <- c("Blinded.ID", "Paternal.Age.at.Proband.Birth", "Maternal.Age.at.Proband.Birth") #c("Blinded.ID", "Paternal.Age.at.Proband.Birth")
cdh_parental_age <- cdh_parental_age[which(cdh_parental_age$Paternal.Age.at.Proband.Birth < 100),]
# Draw combined plot
plot_variants_vs_paternal_age(list(ssc, combine_datasets(combine_datasets(cdh1,cdh2),cdh3), chd3), list(ssc_parental_age, cdh_parental_age, chd_parental_age), batches=list("SSC", c(rep("CDH1",nrow(cdh1)), rep("CDH2",nrow(cdh2)), rep("CDH3",nrow(cdh3))), "TOPMed"), id_match_prefixes=list("^", c("^","CDH","-"), "^"), filename=output_path(paste0("variants_vs_paternal_age.pdf"))) 

# Load CDH phenotype (isolated vs. complex, plus gender) information.
#cdh1_phenotype_info <- merge(cdh[!duplicated(cdh$sample),], read.csv(data_path("CDH1_sampleinfo.csv")), by.x="sample", by.y="StudyID")[,c("sample","Gender","CaseClass")]; colnames(cdh1_phenotype_info) <- c("sample", "gender", "class")
#cdh2_phenotype_info <- merge(cdh2[!duplicated(cdh2$sample),], read.csv(data_path("CDH2_Apr2018_sampleinfo.csv")), by.x="sample", by.y="IID")[,c("sample","SEX","TYPE")]; colnames(cdh2_phenotype_info) <- c("sample", "gender", "class")
#cdh_phenotype_info <- rbind(cdh1_phenotype_info, cdh2_phenotype_info)
#rm(cdh1_phenotype_info); rm(cdh2_phenotype_info)
#cdh_phenotype_info_dhream <- read.csv(data_path("WGS/CDH/CDH_WGS_clinical information_DHREAM.csv"))
#cdh_phenotype_info_mgh <- read.csv(data_path("WGS/CDH/CDH_WGS_clinical information_MGH.csv"))

#blue <- rgb(0, 0, 1, alpha=0.7)
#red <- rgb(1, 0, 0, alpha=0.7)
#barplot(table(denovos[denovos$snv_indel == "snv",]$sample), main="denovos variant counts", xlab="sample", ylab="# variants", col=blue, breaks=seq(1, max(table(denovos[denovos$snv_indel == "snv",]$sample)), by=1), axisnames = FALSE)
#barplot(table(denovos[denovos$snv_indel == "indel",]$sample), main="denovos variant counts", xlab="sample", ylab="# variants", col=red, breaks=seq(1, max(table(denovos[cdh3$snv_indel == "snv",]$sample)), by=1), add=TRUE, axisnames = FALSE)
#legend("topleft", legend=c("SNV", "indel"), col=c(blue, red), pch=15)

# Annotate CADD scores
cdh1_hg19 <- liftover(cdh1, from="hg38", to="hg19", chr_colname="Chrom", start_colname="Position", end_colname=NULL, ref_colname="Ref", alt_colname="Alt", confirm_refseq=FALSE, mismatches_pause=FALSE)[,c("Chrom","Position","Ref","Alt","sample","Chrom_hg38","Position_hg38")]
cdh1_hg19 <- run_annovar(cdh1_hg19, annotations="caddgt10", buildver = "hg19")
x <- merge(cdh1, cdh1_hg19[,c("Chrom_hg38", "Position_hg38", "CADD")], by.x=, all.x=TRUE)
denovos_raw$CHROM <- paste0("chr",denovos_raw$CHROM)


BiocManager::install("GenomicScores")
library(GenomicScores)
availableGScores()
cadd_hg19 <- getGScores("cadd.v1.3.hg19")
gr <- GRanges(seqnames="chr22", IRanges(50967020:50967025, width=1))
gscores(cadd_hg19, gr)
gscores(cadd_hg19, cdh1_granges)

##################################################################################################################################################################
# Write all WGS dataframes for the specified batch to file (called name.tsv, for each respective dataframe name)
##################################################################################################################################################################
write_wgs_dataframes_to_file <- function(batch, directory=output_path("wgs_dataframes")) {
    dir.create(file.path(directory), showWarnings = FALSE)
    print(paste0("Writing ",batch," WGS dataframes to ",directory))
    dat_names <- paste0(tolower(batch), "_", sapply(c("snps", "indels", "muts"), function(x) { paste0(x, c("_wgsa", "_TSS40", "_TSS", "_TES")) }))
    for(dat_name in dat_names) {
        if(dat_name %in% ls(envir=.GlobalEnv)) {
            print(paste0("Writing ",dat_name))
            dat <- get(dat_name)
            write.table(dat, file=full_path(directory, paste0(dat_name,".tsv")), sep="\t", row.names=FALSE)
        }
    }
}
for(batch in c("CDH","CDH2","CDH_combined","SSC_518","SSC_1088","SSC_all","CHD","CHD2","CHD_combined","CHDFB","CHDFB2","CHDFB_combined","SSCFB","SSCFB2","SSCFB_combined","regulatory")) {
    write_wgs_dataframes_to_file(batch)
}
##################################################################################################################################################################
# Load all WGS dataframes for the specified batch from file, if they were previously written.
##################################################################################################################################################################
load_wgs_dataframes_from_file <- function(batch, filter=TRUE, filter_preview=FALSE, directory=output_path("wgs_dataframes")) {
    if(!file.exists(directory)) { print(paste0("WGS dataframe directory ",directory," does not exist!")); return() }
    print(paste0("Loading ",batch," WGS dataframes from ",directory))
    #dat_names <- paste0(tolower(batch), "_", sapply(c("snps", "indels", "muts"), function(x) { paste0(x, c("_wgsa", "_TSS40", "_TSS", "_TES")) }))
    dat_names <- paste0(tolower(batch), "_", sapply(c("muts"), function(x) { paste0(x, c("_wgsa")) }))
    for(dat_name in dat_names) {
        filename = full_path(directory, paste0(dat_name,".tsv"))
        if(file.exists(filename)) {
            print(paste0("Loading ",dat_name))
            dat <- read.csv(file=filename, sep="\t")
            set_global(dat, name=dat_name)
            #assign(dat_name, dat, pos=globalenv())
            
            if(grepl("muts_wgsa", dat_name)) { 
                sample_size_name = paste0(toupper(batch),"_SAMPLE_COUNT")
                print(paste0("Loading ",sample_size_name))
                SAMPLE_SIZE = nrow(dat); 
                set_global(SAMPLE_SIZE, name=sample_size_name) 
            }
        }
    }
    
    if(filter) { filter_dataset(batch, preview=filter_preview) }
}


#################################################################################################################
# RBP relevance ranking by 
#################################################################################################################
rank_RBP_tissue_relevance <- function() {
    eids <- c("E080", "E123", "E118", "E003", "E008", "E083")
    eid_names <- get_roadmap_epigenome_names(eids)
    marks <- c("H3K36me3", "H3K4me1", "H3K4me2", "H3K4me3", "H3K27ac", "H3K27me3", "H3K9me3", "H3K79me2")
    #histone_marks <- apply(expand.grid(eids, paste0(marks,".broadPeak")), 1, paste0, collapse=".")
    tissue_information_table <- matrix(0, nrow=length(eids), ncol=length(marks), dimnames=list(eids, marks))
    genome_length = sum(chr_lengths)
    sapply_out1 <- sapply(eids, function(eid) {
        sapply_out2 <- sapply(marks, function(mark) {
            histone_mark_name <- paste0(eid,".",mark,".broadPeak")
            print(paste0("Processing ",histone_mark_name,"..."))
            histone_mark <- load_annotation(histone_mark_name)
            tissue_information_table[eid, mark] <<- footprint(histone_mark)/genome_length
            print(tissue_information_table)
        })
    })
    print(paste0("Filling in data gaps with mean values..."))
    sapply_out1 <- sapply(1:length(marks), function(i) { tissue_information_table[tissue_information_table[,i] == 0, i] <<- sum(tissue_information_table[,i])/sum(tissue_information_table[,i] > 0) })
    print(tissue_information_table)
    
    rbps <- get_features_by_group("RBP")
    unique_rbps <- sort(unique(paste0(unlist(data.frame(strsplit(rbps, "\\."))[2,]))))
    rbp_relevance_table <- matrix(0, nrow=length(eids), ncol=length(unique_rbps), dimnames=list(eids, unique_rbps))
    sapply_out1 <- sapply(rbps, function(rbp_feature) {
        print(paste0("Processing ",rbp_feature,"..."))
        rbp_split <- strsplit(rbp_feature, "\\.")[[1]]
        rbp = rbp_split[2]
        tissue = rbp_split[1]
        if(tissue == "HepG2") { eid = "E118" 
        } else if (tissue == "K562") { eid = "E123" 
        } else if (tissue == "adrenal_gland") { eid = "E080"
        } else { return(0) }
        histone_mark_name <- paste0(eid,".H3K36me3.broadPeak")
        histone_mark <- load_annotation(histone_mark_name)
        rbp_granges <- load_annotation(rbp_feature)
        rbp_footprint <- footprint(rbp_granges)
        rbp_relevance_table[eid, rbp] <<- footprint(intersect(rbp_granges, histone_mark))/rbp_footprint
        print(rbp_relevance_table)
    })
    rbp_relevance_table <- rbp_relevance_table[,c(colSums(rbp_relevance_table[c("E080", "E123","E118"),] > 0)) > 1]
    
    #tissue_information_table <- cbind(tissue_information_table, rowSums(tissue_information_table)/ncol(tissue_information_table)); colnames(tissue_information_table)[ncol(tissue_information_table)] <- "mean_footprint"
    plot(prcomp(t(tissue_information_table))$rotation[,"PC1"], prcomp(t(tissue_information_table))$rotation[,"PC2"], xlab="PC1", ylab="PC2")
    library(factoextra)
    res.pca <- prcomp(t(tissue_information_table))
    fviz_eig(res.pca, main="Histone mark footprints tissue PCA contributions")
    dev.copy2pdf(file=output_path("tissue_pca_contributions.pdf"))
    fviz_pca_ind(res.pca, col.ind = "cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE, title="Histone mark footprints tissue PCA")
    dev.copy2pdf(file=output_path("tissue_pca.pdf"))
    fviz_pca_biplot(res.pca, repel = TRUE, col.var = "#2E9FDF", col.ind = "#696969", title="Histone mark footprints tissue PCA biplot") 
    dev.copy2pdf(file=output_path("tissue_pca_biplot.pdf"))
    #legend("topleft", legend=paste0(eids,": ",eid_names), cex=0.7, pch=15, col=3)
    write.csv(tissue_information_table, file=output_path("tissue_information_table.csv"))
    write.csv(rbp_relevance_table, file=output_path("rbp_relevance_table.csv"))
    
    
    A <- prcomp(t(tissue_information_table))$rotation[,"PC1"]
    #A <- A / sqrt(sum(A^2)) 
    Y <- rbp_relevance_table[c("E080","E123","E118"),]; #rownames(Y) <- c("E123", "E118")                   # DV matrix
    my_model <- lm(Y ~ A[c("E080","E123","E118")])            # the multivariate model
    my_model
    num_rbps_to_plot = min(c(20, ncol(Y)))
    rbp_cols <- rainbow(num_rbps_to_plot)
    #pdf(output_path("rbp_relevance_tissue_regression.pdf"))
    plot(-10, -10, main="RBP relevance tissue regression", xlab="Tissue PC1 of histone mark footprints", ylab="Fraction of RBP footprint in H3K36me3 regions", xlim=c(min(A), max(A)+0.2), ylim=c(0,1), yaxs="i", type="o", cex.lab=1.4, cex.main=1.2)
    mtext(paste0(ncol(Y)," RBPs with multiple tissue eCLIP data"))
    sapply_out1 <- sapply(1:length(eid_names), function(i) {
        abline(v=A[i], col="black", lty=1, lwd=1)
        text(A[i]+0.005, 0.01+(0.02*i), adj=c(0,0), paste0(eid_names[i]), col="black", cex=0.6) 
    })
    sapply_out1 <- sapply(1:num_rbps_to_plot, function(i) { 
        available_tissues <- c("E080","E123","E118")[Y[c("E080","E123","E118"),i] > 0]
        abline(lm(Y[available_tissues, i] ~ A[available_tissues]), col=rbp_cols[i], lty=3)
        lines(A[available_tissues], Y[available_tissues, i], col=rbp_cols[i], type="o")
        
        # Plot points (with x) of RBP relevance in remaining tissues inferred from mean result using eCLIP data from available tissues.
        rbp = colnames(Y)[i]
        other_tissues <- rownames(rbp_relevance_table)[!(rownames(rbp_relevance_table) %in% available_tissues)]
        sapply_out2 <- sapply(other_tissues, function(eid) {
            print(paste0("Plotting inferred ",eid,".",rbp," relevance derived from available tissue eCLIP data..."))
            histone_mark_name <- paste0(eid,".H3K36me3.broadPeak")
            histone_mark <- load_annotation(histone_mark_name)
            histone_mark_footprint <- footprint(histone_mark)
            
            estimated_rbp_relevance_from_available_eclip <- sapply(available_tissues, function(available_tissue_eid) {
                if(available_tissue_eid == "E118") { available_tissue = "HepG2"
                } else if (available_tissue_eid == "E123") { available_tissue = "K562"
                } else if (available_tissue_eid == "E080") { available_tissue = "adrenal_gland"
                } else { return(0) }
                
                rbp_feature <- paste0(available_tissue,".",rbp)
                rbp_granges <- load_annotation(rbp_feature)
                rbp_footprint <- footprint(rbp_granges)
                
                return(footprint(intersect(rbp_granges, histone_mark))/histone_mark_footprint)  # used to be rbp_footprint in denominator  
            })
            
            rbp_relevance_table[eid,rbp] <<- mean(estimated_rbp_relevance_from_available_eclip)
            
            points(rep(A[eid], length(estimated_rbp_relevance_from_available_eclip)), estimated_rbp_relevance_from_available_eclip, col=rbp_cols[i], pch=4)
        })
    })
    legend("topright", legend=colnames(Y)[1:num_rbps_to_plot], col=rbp_cols, pch=15, cex=0.75, bg="white")
    #dev.off()
    
    # For each tisue, rank RBPs by binding site overlap with H3K36me3!
    ################################### WORKING ON THIS!!!
    rbp_relevance_rankings <- apply(rbp_relevance_table, 1, function(x) { 
        rbp_ranking <- order(x, decreasing=TRUE)
        #x <- formatC(x,format="e",digits=2)
        return(paste0(paste0(colnames(rbp_relevance_table)[rbp_ranking]," (",x[rbp_ranking],")"), collapse=", "))
    })
    write.csv(rbp_relevance_rankings, file=output_path("tissue_rbp_relevance_rankings.csv"))
}


#################################################################################################################
# Look at recurrent H3K36me3_RBP variants
#################################################################################################################
find_recurrent_genes <- function(cases, controls, constraints=c(), filename=NULL, gene_colname="nearest_gene") {
    a <- annotate(cases, c(gene_colname, constraints))
    b <- annotate(controls, c(gene_colname, constraints))
    a_genes <- a$nearest_gene#[a$RBP == "Y" & a$E003.H3K36me3.fullPeak == "Y"]
    b_genes <- b$nearest_gene#[b$RBP == "Y" & b$E003.H3K36me3.fullPeak == "Y"]
    a_genes <- table(a_genes)
    b_genes <- table(b_genes)
    a_genes <- a_genes[a_genes > 1]
    a_genes <- cbind(names(a_genes), a_genes); colnames(a_genes) <- c("gene", "hits_cases")
    b_genes <- cbind(names(b_genes), b_genes); colnames(b_genes) <- c("gene", "count")
    
    a_genes <- cbind(a_genes, unlist(sapply(a_genes[,"gene"], function(x) { pli <- pLI_scores[[paste0(x)]]; if(is.null(pli)) { pli = "." } else { pli = round(pli, 3) }; return(pli) }))); colnames(a_genes)[ncol(a_genes)] <- "pLI"
    a_genes <- cbind(a_genes, unlist(sapply(a_genes[,"gene"], function(x) { mis_z <- mis_z_scores[[paste0(x)]]; if(is.null(mis_z)) { mis_z = "." } else { mis_z = round(mis_z, 3) }; return(mis_z) }))); colnames(a_genes)[ncol(a_genes)] <- "mis_z"
    a_genes <- merge(a_genes, heart_expression_ranks_dat, by=("gene"), all.x=TRUE)
    a_genes[,"heart_rank"] <- paste0(a_genes[,"heart_rank"])
    a_genes[which(is.na(a_genes[,"heart_rank"]) | a_genes[,"heart_rank"] == "NA"), "heart_rank"] <- "."
    
    a_genes <- merge(a_genes, b_genes, by=("gene"), all.x=TRUE); colnames(a_genes)[ncol(a_genes)] <- "hits_controls"
    a_genes[,"hits_controls"] <- paste0(a_genes[,"hits_controls"])
    a_genes[which(is.na(a_genes[,"hits_controls"]) | a_genes[,"hits_controls"] == "NA"), "hits_controls"] <- 0
    
    a_genes <- data.frame(a_genes)
    a_genes$heart_rank <- paste0(a_genes$heart_rank); a_genes$heart_rank[is.na(a_genes$heart_rank) | a_genes$heart_rank == "NA"] <- "."
    
    a_genes <- cbind(a_genes, unlist(sapply(a_genes[,"gene"], function(x) { 
        a_hits <- a[a$nearest_gene == paste0(x),]# & a$RBP == "Y" & a$E003.H3K36me3.fullPeak == "Y",]
        return(paste(paste0(a_hits$sample,"_",a_hits$Chrom_hg38,":",a_hits$Position_hg38,"_",a_hits$Ref,">",a_hits$Alt), collapse=", "))
    })))
    colnames(a_genes)[ncol(a_genes)] <- "case_variants"
    
    a_genes <- cbind(a_genes, unlist(sapply(a_genes[,"gene"], function(x) { 
        b_hits <- b[b$nearest_gene == paste0(x),]# & a$RBP == "Y" & a$E003.H3K36me3.fullPeak == "Y",]
        if(nrow(b_hits) > 0) {
            return(paste(paste0(b_hits$sample,"_",b_hits$Chrom_hg38,":",b_hits$Position_hg38,"_",b_hits$Ref,">",b_hits$Alt), collapse=", "))
        } else { return(".") }
    })))
    colnames(a_genes)[ncol(a_genes)] <- "control_variants"
    
    a_genes <- a_genes[order(a_genes$hits_cases, decreasing=TRUE),]
    write.csv(a_genes, file=filename, row.names=FALSE)
    return(a_genes)
}
find_recurrent_genes(cases=chdfb_combined, controls=sscfb_combined, constraints=c("RBP", "E003.H3K36me3.fullPeak"), filename=output_path("chdfb_rbp_h3k36me3_recurrent_genes.csv"))

# Recurrent gene analysis for union of PCGC manuscript approaches
filename_suffixes = c("", "_rbp_top_hit", "_rbp_constrained_hit")
for(j in 1:length(filename_suffixes)) {
    rbp_variants_dat <- data.frame(read_excel(data_path("chd_manuscript_supplementary_data/chd_manuscript_supplementary_data_2019_02_22.xlsx"), sheet="S15 RBP enrichment"))
    if(j > 1) { rbp_variants_dat <- rbp_variants_dat[rbp_variants_dat[,(5+j-1)] == "T",] }
    filename_suffix = filename_suffixes[j]
    print(paste0(j,": ",filename_suffix))
    
    colnames(rbp_variants_dat)[1:2] <- c("Chrom_hg38", "Position_hg38")
    heartenn_variants_dat <- data.frame(read_excel(data_path("chd_manuscript_supplementary_data/chd_manuscript_supplementary_data_2019_02_22.xlsx"), sheet="S8 Top HeartENN score")) # S14 Multi-gene-hit DNVs # S15 RBP enrichment
    heartenn_variants_dat <- heartenn_variants_dat[heartenn_variants_dat$HeartENN_score >= 0.1,]; colnames(heartenn_variants_dat)[1:2] <- c("Chrom_hg38", "Position_hg38")
    multigene_variants_dat <- data.frame(read_excel(data_path("chd_manuscript_supplementary_data/chd_manuscript_supplementary_data_2019_02_22.xlsx"), sheet="S14 Multi-gene-hit DNVs"))
    
    heartenn_variants_dat <- standardize_colnames(heartenn_variants_dat[,1:5])
    multigene_variants_dat <- standardize_colnames(multigene_variants_dat[,1:5])
    rbp_variants_dat <- standardize_colnames(rbp_variants_dat[,1:5])
    heartenn_variants_dat[1:5,]
    multigene_variants_dat[1:5,]
    rbp_variants_dat[1:5,]
    nrow(heartenn_variants_dat)
    nrow(multigene_variants_dat)
    nrow(rbp_variants_dat)
    
    union_variants_dat <- rbind(heartenn_variants_dat, multigene_variants_dat, rbp_variants_dat)
    nrow(union_variants_dat)
    union_variants_dat <- unique(union_variants_dat)
    nrow(union_variants_dat)
    
    union_variants <- unlist(strsplit(genomic_dat_to_strings(rbp_variants_dat, chr_colname="Chrom_hg38", pos_colname="Position_hg38", ref_colname="Ref", alt_colname="Alt"), ";"))
    chdfb_combined_variants <- unlist(strsplit(genomic_dat_to_strings(chdfb_combined_muts_wgsa, chr_colname="Chrom_hg38", pos_colname="Position_hg38", ref_colname="ref", alt_colname="alt"), ";"))
    chdfb_union_variants <- chdfb_combined_muts_wgsa[chdfb_combined_variants %in% union_variants,]
    nrow(chdfb_union_variants)
    sscfb_combined_variants <- unlist(strsplit(genomic_dat_to_strings(sscfb_combined_muts_wgsa, chr_colname="Chrom_hg38", pos_colname="Position_hg38", ref_colname="ref", alt_colname="alt"), ";"))
    sscfb_union_variants <- sscfb_combined_muts_wgsa[sscfb_combined_variants %in% union_variants,]
    nrow(sscfb_union_variants)
    
    recurrent_union_genes <- find_recurrent_genes(cases=chdfb_union_variants, controls=sscfb_union_variants, filename=output_path(paste0("chdfb_union_analyses_recurrent_genes",filename_suffix,".csv")))
    
    chd_coding_mutations <- data.frame(read_excel(data_path("ng.3970-S3.xlsx"), sheet="S9", skip=1))
    chd_coding_mutations <- standardize_colnames(chd_coding_mutations[chd_coding_mutations$Variant_Class != "syn",])
    chd_coding_mutations <- cbind(chd_coding_mutations, unlist(strsplit(genomic_dat_to_strings(chd_coding_mutations, chr_colname="Chrom", pos_colname="Position", ref_colname="Ref", alt_colname="Alt", notes_colname="Variant_Class"), ";"))); colnames(chd_coding_mutations)[ncol(chd_coding_mutations)] <- "variant"
    coding_mutations <- unlist(sapply(recurrent_union_genes$gene, function(x) { coding_muts <- chd_coding_mutations$variant[chd_coding_mutations$Gene == x]; if(length(coding_muts) == 0) { coding_muts = "." }; return(paste0(sort(coding_muts), collapse=", ")) }))
    
    chd_genes <- data.frame(read_excel(data_path("chd_manuscript_supplementary_data/chd_manuscript_supplementary_data_2019_02_22.xlsx"), sheet="S5 CHD gene lists"))
    human_chd_genes <- unique(paste0(chd_genes[!is.na(chd_genes[,2]),2]))
    mouse_chd_genes <- unique(paste0(chd_genes[!is.na(chd_genes[,3]),3]))
    
    human_chd_gene <- rep("N", nrow(recurrent_union_genes)); human_chd_gene[recurrent_union_genes$gene %in% human_chd_genes] <- "Y"
    mouse_chd_gene <- rep("N", nrow(recurrent_union_genes)); mouse_chd_gene[recurrent_union_genes$gene %in% mouse_chd_genes] <- "Y"
    chd_pathways <- unlist(sapply(recurrent_union_genes$gene, function(x) { pathways <- unique(chd_pathway_genes[grepl(x, chd_pathway_genes$gene),]$pathway); if(is.null(pathways)) { pathways = "." }; return(paste0(sort(pathways), collapse=", ")) }))
    binom_results <- data.frame(t(sapply(1:nrow(recurrent_union_genes), function(i) {
        binom_test_result <- binomial_test(as.numeric(paste0(recurrent_union_genes$hits_cases[i])), as.numeric(paste0(recurrent_union_genes$hits_controls[i])), length(unique(chdfb_combined_muts_wgsa$sample)), length(unique(sscfb_combined_muts_wgsa$sample)))
        return(c(binom_test_result[["estimate"]], binom_test_result[["p.value"]]))
    }))); colnames(binom_results) <- c("binom_enrichment", "binom_p.value")
    # Look at phenotypes
    chd_phenotypes <- standardize_colnames(data.frame(read_excel(data_path("chd_manuscript_supplementary_data/chd_manuscript_supplementary_data_2019_02_22.xlsx"), sheet="S1 Phenotypes")))
    phenotypes <- unlist(sapply(recurrent_union_genes$case_variants, function(x) {
        samples <- sort(unique(sapply(unlist(strsplit(paste0(x), ", ")), function(case_variant) strsplit(case_variant, "_")[[1]][1])))
        sample_phenotypes <- sapply(samples, function(sample) {
            chd_phenotypes_entry <- chd_phenotypes[chd_phenotypes$sample == sample,]
            paste0(chd_phenotypes_entry[1],": [ ",paste0(names(chd_phenotypes_entry[-1]),": ",chd_phenotypes_entry[-1], collapse="; ")," ]")
        })
        return(paste0(sample_phenotypes, collapse=", "))
    }))
    recurrent_union_genes <- cbind(recurrent_union_genes, coding_mutations, human_chd_gene, mouse_chd_gene, chd_pathways, binom_results, phenotypes)
    recurrent_union_genes <- recurrent_union_genes[,c("gene", "hits_cases", "pLI", "mis_z", "heart_rank", "human_chd_gene", "mouse_chd_gene", "chd_pathways", "coding_mutations", "hits_controls", "binom_enrichment", "binom_p.value", "case_variants", "control_variants", "phenotypes")]
    # Define score metric for internal gene ranking
    score <- rep(0, nrow(recurrent_union_genes))
    recurrent_union_genes$pLI <- paste0(recurrent_union_genes$pLI); recurrent_union_genes$mis_z <- paste0(recurrent_union_genes$mis_z)
    recurrent_union_genes$pLI[recurrent_union_genes$pLI == "."] <- "0"; recurrent_union_genes$mis_z[recurrent_union_genes$mis_z == "."] <- "0"
    score[as.numeric(recurrent_union_genes$pLI) > 0.5 | as.numeric(recurrent_union_genes$mis_z) > 3] <- score[as.numeric(recurrent_union_genes$pLI) > 0.5 | as.numeric(recurrent_union_genes$mis_z) > 3] + 1
    score[recurrent_union_genes$heart_rank != "." & as.numeric(recurrent_union_genes$heart_rank) < 0.25] <- score[recurrent_union_genes$heart_rank != "." & as.numeric(recurrent_union_genes$heart_rank) < 0.25] + 1
    score[recurrent_union_genes$human_chd_gene == "Y" | recurrent_union_genes$mouse_chd_gene == "Y" | recurrent_union_genes$chd_pathways != ""] <- score[recurrent_union_genes$human_chd_gene == "Y" | recurrent_union_genes$mouse_chd_gene == "Y" | recurrent_union_genes$chd_pathways != ""] + 3
    score[recurrent_union_genes$coding_mutations != "."] <- score[recurrent_union_genes$coding_mutations != "."] + 2
    recurrent_union_genes <- cbind(recurrent_union_genes, score)
    recurrent_union_genes <- recurrent_union_genes[order(recurrent_union_genes$score, decreasing=TRUE),]
    # Select only marginally significantly enriched (case vs. control) genes
    recurrent_union_genes <- recurrent_union_genes[recurrent_union_genes$binom_p.value < 0.05,]
    write.csv(recurrent_union_genes, file=output_path(paste0("chdfb_union_analyses_recurrent_genes",filename_suffix,".csv")), row.names=FALSE)
}


# Look at recurrent H3K36me3_RBP or H3K4me1_indels variants
H3K36me3_RBP_or_H3K4me1_indels_variants <- paste0(read.table(output_path("recurrent_genes_H3K36me3_RBP_or_H3K4me1_indels_variants.txt"))[,1])

H3K36me3_RBP_or_H3K4me1_indels_variants_dat <- t(sapply(H3K36me3_RBP_or_H3K4me1_indels_variants, function(variant_gene) { 
    variant_gene_strsplit <- strsplit(variant_gene, "_")[[1]]; variant = variant_gene_strsplit[1]; gene = variant_gene_strsplit[2]
    mutation = gsub("[^ACTG>]", "", variant)
    location = unlist(strsplit(gsub("[chrACTG>]", "", variant), ":")) # "[ACTG]+>[ACTG]" also works for regex
    chromosome_hg19 = location[1]; position_hg19 = location[2]
    match_index = chdfb_combined_muts_wgsa$X.chr == chromosome_hg19 & chdfb_combined_muts_wgsa$pos == position_hg19
    chromosome_hg38 = paste0(chdfb_combined_muts_wgsa$Chrom_hg38[match_index])[1]; position_hg38 = paste0(chdfb_combined_muts_wgsa$Position_hg38[match_index])[1]
    return(c(paste0(chromosome_hg38,":",position_hg38,mutation), variant, gene))
}))
chdfb_all_variants <- paste0("chr",chdfb_combined_muts_wgsa$X.chr,":",chdfb_combined_muts_wgsa$pos,chdfb_combined_muts_wgsa$ref,">",chdfb_combined_muts_wgsa$alt)
names(chdfb_all_variants) <- paste0(chdfb_combined_muts_wgsa$sample)
H3K36me3_RBP_or_H3K4me1_indels_variants_dat <- cbind(sapply(H3K36me3_RBP_or_H3K4me1_indels_variants_dat[,2], function(x) { names(chdfb_all_variants)[grepl(x, chdfb_all_variants)] }), H3K36me3_RBP_or_H3K4me1_indels_variants_dat)
colnames(H3K36me3_RBP_or_H3K4me1_indels_variants_dat) <- c("sample", "hg38", "hg19", "gene")
rownames(H3K36me3_RBP_or_H3K4me1_indels_variants_dat) <- c()

sarah_annotation_variable_names <- load(data_path("Sarah_Morton_annotated_wgs_31Jul18.rda"))
denovo_all_annotated <- cbind(paste0(denovo_all_annotated$seqnames,":",denovo_all_annotated$start,denovo_all_annotated$Ref,">",denovo_all_annotated$Alt), denovo_all_annotated)
colnames(denovo_all_annotated)[1] <- "variant"

H3K36me3_RBP_or_H3K4me1_indels_variants_dat <- merge(H3K36me3_RBP_or_H3K4me1_indels_variants_dat, denovo_all_annotated, by.x=c("hg38","sample"), by.y=c("variant","Blinded.ID"))[,c("sample", "hg38", "hg19", "gene", "nearest_gene", "upstream_gene", "downstream_gene")]
colnames(H3K36me3_RBP_or_H3K4me1_indels_variants_dat)[(ncol(H3K36me3_RBP_or_H3K4me1_indels_variants_dat)-2):ncol(H3K36me3_RBP_or_H3K4me1_indels_variants_dat)] <- c("Sarah_nearest_gene", "Sarah_upstream_gene", "Sarah_downstream_gene")
write.csv(H3K36me3_RBP_or_H3K4me1_indels_variants_dat, file=output_path("H3K36me3_RBP_or_H3K4me1_indels_variants.csv"), row.names=FALSE)


# CHD1 AC
chd1_ac <- chd1_ac[chd1_ac$Chrom %in% standard_chromosomes & (as.numeric(unlist(lapply(strsplit(paste0(chd1_ac$AC),","), function(x) max(x)))) > 1), c("Chrom", "Position", "Ref", "Alt", "sample", "AC")]; chd1_ac$sample <- unlist(lapply(strsplit(paste0(chd1_ac$sample), "\\("), function(x) x[[1]][1])); chd1_ac <- chd1_ac[chd1_ac$sample %in% unique(chdfb$sample),]
chd1_ac_granges <- to_genomic_regions(chd1_ac, chr_colname="Chrom", start_colname="Position", end_colname="Position", labels=is_indel(chd1_ac$Ref, chd1_ac$Alt)); names(chd1_ac_granges)[names(chd1_ac_granges) == "TRUE"] <- "indel"; names(chd1_ac_granges)[names(chd1_ac_granges) == "FALSE"] <- "snv"
chd1_ac <- chd1_ac[-c(unique(c(queryHits(findOverlaps(chd1_ac_granges, mappability_hg19_granges)), queryHits(findOverlaps(chd1_ac_granges, lcr_hg19_granges)), queryHits(findOverlaps(chd1_ac_granges, segdups_hg19_granges)),queryHits(findOverlaps(chd1_ac_granges, duke_blacklist_hg19_granges)),queryHits(findOverlaps(chd1_ac_granges, dac_blacklist_hg19_granges))))),]
chd1_ac_granges <- to_genomic_regions(chd1_ac, chr_colname="Chrom", start_colname="Position", end_colname="Position", labels=is_indel(chd1_ac$Ref, chd1_ac$Alt)); names(chd1_ac_granges)[names(chd1_ac_granges) == "TRUE"] <- "indel"; names(chd1_ac_granges)[names(chd1_ac_granges) == "FALSE"] <- "snv"
# CHD2 AC
chd2_ac <- chd2_ac[chd2_ac$Chrom %in% paste0("chr",standard_chromosomes) & (as.numeric(unlist(lapply(strsplit(paste0(chd2_ac$AC),","), function(x) max(x)))) > 1), c("Chrom", "Position", "Ref", "Alt", "sample", "AC")]; chd2_ac$Chrom <- gsub("chr","",chd2_ac$Chrom); chd2_ac$sample <- unlist(lapply(strsplit(paste0(chd2_ac$sample), "\\("), function(x) x[[1]][1])); chd2_ac <- chd2_ac[chd2_ac$sample %in% unique(chdfb2$sample),]
chd2_ac_granges <- to_genomic_regions(chd2_ac, chr_colname="Chrom", start_colname="Position", end_colname="Position", labels=is_indel(chd2_ac$Ref, chd2_ac$Alt)); names(chd2_ac_granges)[names(chd2_ac_granges) == "TRUE"] <- "indel"; names(chd2_ac_granges)[names(chd2_ac_granges) == "FALSE"] <- "snv"
chd2_ac <- chd2_ac[-c(unique(c(queryHits(findOverlaps(chd2_ac_granges, mappability_hg38_granges)), queryHits(findOverlaps(chd2_ac_granges, lcr_hg38_granges)), queryHits(findOverlaps(chd2_ac_granges, segdups_hg38_granges)),queryHits(findOverlaps(chd2_ac_granges, duke_blacklist_hg38_granges)),queryHits(findOverlaps(chd2_ac_granges, dac_blacklist_hg38_granges))))),]
chd2_ac_granges <- to_genomic_regions(chd2_ac, chr_colname="Chrom", start_colname="Position", end_colname="Position", labels=is_indel(chd2_ac$Ref, chd2_ac$Alt)); names(chd2_ac_granges)[names(chd2_ac_granges) == "TRUE"] <- "indel"; names(chd2_ac_granges)[names(chd2_ac_granges) == "FALSE"] <- "snv"
# SSC1 AC
ssc1_ac <- ssc1_ac[ssc1_ac$Chrom %in% standard_chromosomes & (as.numeric(unlist(lapply(strsplit(paste0(ssc1_ac$AC),","), function(x) max(x)))) > 1), c("Chrom", "Position", "Ref", "Alt", "sample", "AC")]; ssc1_ac$sample <- unlist(lapply(strsplit(paste0(ssc1_ac$sample), "\\("), function(x) x[[1]][1])); ssc1_ac <- ssc1_ac[ssc1_ac$sample %in% unique(sscfb$sample),]
ssc1_ac_granges <- to_genomic_regions(ssc1_ac, chr_colname="Chrom", start_colname="Position", end_colname="Position", labels=is_indel(ssc1_ac$Ref, ssc1_ac$Alt)); names(ssc1_ac_granges)[names(ssc1_ac_granges) == "TRUE"] <- "indel"; names(ssc1_ac_granges)[names(ssc1_ac_granges) == "FALSE"] <- "snv"
ssc1_ac <- ssc1_ac[-c(unique(c(queryHits(findOverlaps(ssc1_ac_granges, mappability_hg19_granges)), queryHits(findOverlaps(ssc1_ac_granges, lcr_hg19_granges)), queryHits(findOverlaps(ssc1_ac_granges, segdups_hg19_granges)),queryHits(findOverlaps(ssc1_ac_granges, duke_blacklist_hg19_granges)),queryHits(findOverlaps(ssc1_ac_granges, dac_blacklist_hg19_granges))))),]
ssc1_ac_granges <- to_genomic_regions(ssc1_ac, chr_colname="Chrom", start_colname="Position", end_colname="Position", labels=is_indel(ssc1_ac$Ref, ssc1_ac$Alt)); names(ssc1_ac_granges)[names(ssc1_ac_granges) == "TRUE"] <- "indel"; names(ssc1_ac_granges)[names(ssc1_ac_granges) == "FALSE"] <- "snv"
# SSC2 AC
ssc2_ac <- ssc2_ac[ssc2_ac$Chrom %in% paste0("chr",standard_chromosomes) & (as.numeric(unlist(lapply(strsplit(paste0(ssc2_ac$AC),","), function(x) max(x)))) > 1), c("Chrom", "Position", "Ref", "Alt", "sample", "AC")]; ssc2_ac$Chrom <- gsub("chr","",ssc2_ac$Chrom)
ssc2_ac$sample <- unlist(lapply(strsplit(paste0(ssc2_ac$sample), "\\("), function(x) x[[1]][1])); ssc2_ac$sample <- sapply(ssc2_ac$sample, function(sample) { mapped_sample <- ssc_id_conversion$V1[ssc_id_conversion$V2 == sample]; if(length(mapped_sample)>0) { return(mapped_sample[1]) } else { return(sample) } }); ssc2_ac <- ssc2_ac[ssc2_ac$sample %in% unique(sscfb2$sample),]
ssc2_ac_granges <- to_genomic_regions(ssc2_ac, chr_colname="Chrom", start_colname="Position", end_colname="Position", labels=is_indel(ssc2_ac$Ref, ssc2_ac$Alt)); names(ssc2_ac_granges)[names(ssc2_ac_granges) == "TRUE"] <- "indel"; names(ssc2_ac_granges)[names(ssc2_ac_granges) == "FALSE"] <- "snv"
ssc2_ac <- ssc2_ac[-c(unique(c(queryHits(findOverlaps(ssc2_ac_granges, mappability_hg38_granges)), queryHits(findOverlaps(ssc2_ac_granges, lcr_hg38_granges)), queryHits(findOverlaps(ssc2_ac_granges, segdups_hg38_granges)),queryHits(findOverlaps(ssc2_ac_granges, duke_blacklist_hg38_granges)),queryHits(findOverlaps(ssc2_ac_granges, dac_blacklist_hg38_granges))))),]
ssc2_ac_granges <- to_genomic_regions(ssc2_ac, chr_colname="Chrom", start_colname="Position", end_colname="Position", labels=is_indel(ssc2_ac$Ref, ssc2_ac$Alt)); names(ssc2_ac_granges)[names(ssc2_ac_granges) == "TRUE"] <- "indel"; names(ssc2_ac_granges)[names(ssc2_ac_granges) == "FALSE"] <- "snv"

chd1_ac <- chd1_ac[-c(unique(queryHits(findOverlaps(chd1_ac_granges, chdfb1_hg19_granges)))),]; chd1_ac$Chrom <- paste0("chr", chd1_ac$Chrom)
chd2_ac <- chd2_ac[-c(unique(queryHits(findOverlaps(chd2_ac_granges, chdfb2_hg38_granges)))),]; chd2_ac$Chrom <- paste0("chr", chd2_ac$Chrom)
ssc1_ac <- ssc1_ac[-c(unique(queryHits(findOverlaps(ssc1_ac_granges, sscfb1_hg19_granges)))),]; ssc1_ac$Chrom <- paste0("chr", ssc1_ac$Chrom)
ssc2_ac <- ssc2_ac[-c(unique(queryHits(findOverlaps(ssc2_ac_granges, sscfb2_hg38_granges)))),]; ssc2_ac$Chrom <- paste0("chr", ssc2_ac$Chrom)
write.csv(chd1_ac, output_path("GMKF1_candidate_AC_variants_hg19.csv"), row.names=FALSE)
write.csv(chd2_ac, output_path("GMKF2_candidate_AC_variants_hg38.csv"), row.names=FALSE)
write.csv(ssc1_ac, output_path("SSC1_candidate_AC_variants_hg19.csv"), row.names=FALSE)
write.csv(ssc2_ac, output_path("SSC2_candidate_AC_variants_hg38.csv"), row.names=FALSE)

length(sscfb2_granges) - length(unique(queryHits(findOverlaps(sscfb2_granges, ssc_1088_granges))))
length(ssc_1088_granges) - length(unique(subjectHits(findOverlaps(sscfb2_granges, ssc_1088_granges))))
sscfb2_unique <- sscfb2_muts_wgsa[-c(unique(queryHits(findOverlaps(sscfb2_granges, ssc_1088_granges)))),]
chr_temp <- sscfb2_unique$X.chr; pos_temp <- sscfb2_unique$pos; sscfb2_unique$X.chr <- sscfb2_unique$Chrom_hg38; sscfb2_unique$pos <- sscfb2_unique$Position_hg38; sscfb2_unique$Chrom_hg38 <- chr_temp; sscfb2_unique$Position_hg38 <- pos_temp; colnames(sscfb2_unique)[c(1,2)] <- c("X.chr", "pos"); colnames(sscfb2_unique)[which(colnames(sscfb2_unique) %in% c("Chrom_hg38","Position_hg38"))] <- c("X.chr_hg19", "pos_hg19")
ssc2_unique <- ssc_1088_muts_wgsa[-c(unique(subjectHits(findOverlaps(sscfb2_granges, ssc_1088_granges)))),]
ssc2_unique <- liftover(ssc2_unique, from="hg19", to="hg38", chr_colname="X.chr", start_colname="pos")
sscfb2_unique <- sscfb2_unique[,c("X.chr", "pos", "ref", "alt", "sample", "E008.H3K36me3.gappedPeak", "RBP", "RBP0", "X.chr_hg19", "pos_hg19", "ANNOVAR_ensembl_summary", "Simons.FamID")]
ssc2_unique <- ssc2_unique[,c("X.chr", "pos", "ref", "alt", "sample", "E008.H3K36me3.gappedPeak", "RBP", "RBP0", "X.chr_hg19", "pos_hg19")]; ssc2_unique$X.chr <- paste0("chr",ssc2_unique$X.chr)
ssc2_unique <- cbind(ssc2_unique, sapply(ssc2_unique$sample, function(x) ssc_id_conversion$V2[ssc_id_conversion$V1 == x])); colnames(ssc2_unique)[ncol(ssc2_unique)] <- "Simons.FamID"

ssc_columbia_unique_call_filtering_results <- read.csv(data_path("WGS/columbia_control_call_info/ssc_columbia_unique_call_filtering_results.txt"), sep="\t")
ssc_hs_unique_call_reasons <- read.csv(data_path("WGS/columbia_control_call_info/ssc_hs_unique_call_reason.txt"), sep="\t")

sum(paste0(sscfb2_unique$Simons.FamID,"_",sscfb2_unique$X.chr,"_",sscfb2_unique$pos,"_",sscfb2_unique$ref,"_",sscfb2_unique$alt) %in% ssc_hs_unique_call_reasons$var_id)
sum(paste0(ssc2_unique$Simons.FamID,"_",ssc2_unique$X.chr,"_",ssc2_unique$pos,"_",ssc2_unique$ref,"_",ssc2_unique$alt) %in% ssc_columbia_unique_call_filtering_results$cmb_var_id)
ssc_hs_unique_call_reasons[ssc_hs_unique_call_reasons$var_id %in% paste0(sscfb2_unique$Simons.FamID,"_",sscfb2_unique$X.chr,"_",sscfb2_unique$pos,"_",sscfb2_unique$ref,"_",sscfb2_unique$alt),]
a <- ssc_hs_unique_call_reasons[ssc_hs_unique_call_reasons$var_id %in% paste0(sscfb2_unique$Simons.FamID,"_",sscfb2_unique$X.chr,"_",sscfb2_unique$pos,"_",sscfb2_unique$ref,"_",sscfb2_unique$alt),]
a <- a[!(grepl("coding_variants", a$columbia_exclude_reason)),]
write.csv(a, file=output_path("ssc_hs_unique_call_variants_failed_for_other_reasons.csv"), row.names=FALSE)
a_annotated <- lapply(strsplit(paste0(a$var_id), "_"), function(x) {  
    dat <- sscfb_combined_muts_wgsa[sscfb_combined_muts_wgsa$Chrom_hg38 == x[2] & sscfb_combined_muts_wgsa$Position_hg38 == x[3],]
    return(dat[,c("Chrom_hg38", "Position_hg38", "ref", "alt", "snv_indel", "E008.H3K36me3.gappedPeak","RBP")])
})
a_annotated <- data.frame(Reduce(rbind, a_annotated))
a_annotated[a_annotated$E008.H3K36me3.gappedPeak=="Y" & !is.na(a_annotated$RBP),]

sscfb_combined_muts_wgsa[sscfb_combined_muts_wgsa$E008.H3K36me3.gappedPeak=="Y" & !is.na(sscfb_combined_muts_wgsa$RBP),]
ssc_all_muts_wgsa[ssc_all_muts_wgsa$E008.H3K36me3.gappedPeak=="Y" & !is.na(ssc_all_muts_wgsa$RBP),]
b <- sscfb_combined_muts_wgsa[sscfb_combined_muts_wgsa$E008.H3K36me3.gappedPeak=="Y" & !is.na(sscfb_combined_muts_wgsa$RBP),][!(paste(sscfb_combined_muts_wgsa$sample[sscfb_combined_muts_wgsa$E008.H3K36me3.gappedPeak=="Y" & !is.na(sscfb_combined_muts_wgsa$RBP)], sscfb_combined_muts_wgsa$X.chr[sscfb_combined_muts_wgsa$E008.H3K36me3.gappedPeak=="Y" & !is.na(sscfb_combined_muts_wgsa$RBP)], sscfb_combined_muts_wgsa$pos[sscfb_combined_muts_wgsa$E008.H3K36me3.gappedPeak=="Y" & !is.na(sscfb_combined_muts_wgsa$RBP)]) %in% paste(ssc_all_muts_wgsa$sample[ssc_all_muts_wgsa$E008.H3K36me3.gappedPeak=="Y" & !is.na(ssc_all_muts_wgsa$RBP)], ssc_all_muts_wgsa$X.chr[ssc_all_muts_wgsa$E008.H3K36me3.gappedPeak=="Y" & !is.na(ssc_all_muts_wgsa$RBP)], ssc_all_muts_wgsa$pos[ssc_all_muts_wgsa$E008.H3K36me3.gappedPeak=="Y" & !is.na(ssc_all_muts_wgsa$RBP)])),]
b_reasons <- ssc_hs_unique_call_reasons[lapply(strsplit(paste0(ssc_hs_unique_call_reasons$var_id),"_"), function(x) paste0(x[1],"_",x[2],"_",x[3])) %in% paste0(sapply(b$sample, function(x) ssc_id_conversion$V2[ssc_id_conversion$V1 == x]),"_",b$Chrom_hg38,"_",b$Position_hg38),] 
b_reasons$columbia_exclude_reason <- paste0(b_reasons$columbia_exclude_reason)
b_reasons$columbia_exclude_reason[which(grepl("High FS", b_reasons$columbia_exclude_reason))] <- "High FS"
b_reasons$columbia_exclude_reason[which(grepl("Low QD", b_reasons$columbia_exclude_reason))] <- "Low QD"
table(b_reasons$columbia_exclude_reason)
write.csv(b_reasons[!(b_reasons$var_id %in% a$var_id) & !(grepl("coding_variants", b_reasons$columbia_exclude_reason)),], file=output_path("ssc_hs_unique_call_variants_failed_for_other_reasons.csv"), row.names=FALSE)


ssc_columbia_unique_call_filtering_results[ssc_columbia_unique_call_filtering_results$cmb_var_id %in% paste0(ssc2_unique$Simons.FamID,"_",ssc2_unique$X.chr,"_",ssc2_unique$pos,"_",ssc2_unique$ref,"_",ssc2_unique$alt),]

ssc2_unique_granges <- to_genomic_regions(ssc2_unique, chr_colname="X.chr", start_colname="pos", end_colname="pos")
sscfb2_unique_granges <- to_genomic_regions(sscfb2_unique, chr_colname="X.chr", start_colname="pos", end_colname="pos")

# CDH samples to select for RNAseq 
cdh_lan <- cdh_muts_wgsa[,!grepl("narrowPeak", colnames(cdh_muts_wgsa)) & !(grepl("Peak", colnames(cdh_muts_wgsa)) & !grepl(paste0(get_relevant_roadmap_eids("CDH"), collapse="|"), colnames(cdh_muts_wgsa)))]
cdh_lan2 <- cdh2_muts_wgsa[,!grepl("narrowPeak", colnames(cdh2_muts_wgsa)) & !(grepl("Peak", colnames(cdh2_muts_wgsa)) & !grepl(paste0(get_relevant_roadmap_eids("CDH"), collapse="|"), colnames(cdh2_muts_wgsa)))]
write.csv(cdh_lan, output_path("cdh1_noncoding_variants.csv"), row.names=FALSE)
write.csv(cdh_lan2, output_path("cdh2_noncoding_variants.csv"), row.names=FALSE)

TFRBP <- read.csv(data_path("TFRBP.txt"), sep="\t"); TFRBP$gene_symbol <- paste0(TFRBP$gene_symbol)
transcription_factor_genes <- TFRBP$gene_symbol[includes_at_least_one(TFRBP$entrez_human, "TF")]
rbp_genes <- TFRBP$gene_symbol[includes_at_least_one(TFRBP$entrez_human, "RBP")]
splicing_regulator_genes <- TFRBP$gene_symbol[includes_at_least_one(TFRBP$entrez_human, "splicing regulator")]

get_interesting_variants <- function(dat) {
    dat <- dat[!is.na(dat$RBP) & apply(dat[,grepl("H3K36me3",colnames(dat))], 1, function(x) sum(x=="Y")>0),c("X.chr","pos","ref","alt","sample","ANNOVAR_refseq_summary","ANNOVAR_refseq_Closest_gene.intergenic_only.","phastCons46way_placental","RBP",colnames(dat)[grepl("H3K36me3",colnames(dat))]),]
    dat$phastCons46way_placental[dat$phastCons46way_placental=="."] <- 0; dat$phastCons46way_placental <- as.numeric(paste0(dat$phastCons46way_placental)); evolutionarily_conserved <- rep("N",nrow(dat)); evolutionarily_conserved[dat$phastCons46way_placental > 0.5] <- "Y"; dat$phastCons46way_placental <- evolutionarily_conserved 
    dat$ANNOVAR_refseq_summary <- gsub("\\)[^\\|]+", ")", gsub("\\(\\d+\\)", ")", gsub("\\(\\d+\\):", "(", dat$ANNOVAR_refseq_summary))); dat$ANNOVAR_refseq_summary[grepl("intergenic", dat$ANNOVAR_refseq_summary)] <- "."
    dat$ANNOVAR_refseq_Closest_gene.intergenic_only. <- gsub(",*NONE.*\\),*", "", gsub(":N[M|R]_\\d+|\\.|\\|", "", dat$ANNOVAR_refseq_Closest_gene.intergenic_only.))
    dat$ANNOVAR_refseq_Closest_gene.intergenic_only.[dat$ANNOVAR_refseq_Closest_gene.intergenic_only. == ""] <- dat$ANNOVAR_refseq_summary[dat$ANNOVAR_refseq_Closest_gene.intergenic_only. == ""] # dat[,-which(colnames(dat) == "ANNOVAR_refseq_summary")]
    variant <- paste0("chr",dat[,1],":",dat[,2],"_",dat[,3],">",dat[,4]); dat <- cbind(variant, dat)
    colnames(dat)[2:9] <- c("Chrom", "Position", "Ref", "Alt", "sample", "genes", "genes_with_trans_effect", "evolutionary_conservation")
    colnames(dat)[grepl("H3K36me3",colnames(dat))] <- sapply(colnames(dat)[grepl("H3K36me3",colnames(dat))], function(x) {
        histone_feature_split <- strsplit(x, "\\.")[[1]]
        eid <- histone_feature_split[1]
        tissue <- get_roadmap_epigenome_names(eid)
        if(tissue == "-") { tissue <- eid; eid <- paste0(relevant_eids, collapse=",") } # special case for combined tissue features, where the eid is initially "relevant"
        histone_mark <- histone_feature_split[2]
        peak_type <- histone_feature_split[3]
        if (!(grepl("me", histone_mark) | grepl("ac", histone_mark))) { next }
        test_name <- paste0(histone_mark, " in ", tissue, " (", eid, ", ", peak_type, ")")
    })
    #dat$genes_with_trans_effect <- rep(".", nrow(dat))
    dat$genes_with_trans_effect <- paste0(unlist(lapply(strsplit(gsub("\\([^\\|]*\\)", "", dat$genes), "\\|"), function(x) { return(paste0(sapply(x[x %in% TFRBP$gene_symbol], function(x_TFRBP) paste0(x_TFRBP,"(",TFRBP$entrez_human[TFRBP$gene_symbol == x_TFRBP],")") ), collapse=",") ) })))
    dat$genes_with_trans_effect[dat$genes_with_trans_effect == ""] <- "."
    return(dat)
}
cdh_lan <- get_interesting_variants(cdh_lan)
cdh_lan2 <- get_interesting_variants(cdh_lan2)
write.csv(cdh_lan, output_path("cdh1_h3k36me3_rbp_disturbing_variants.csv"), row.names=FALSE)
write.csv(cdh_lan2, output_path("cdh2_h3k36me3_rbp_disturbing_variants.csv"), row.names=FALSE)

get_interesting_samples <- function(dat) {
    dat_samples <- aggregate(dat[-c(2:6)], by=list(dat$sample), FUN=paste, collapse=";"); colnames(dat_samples)[1:6] <- c("sample", "variants", "genes", "genes_with_trans_effect", "evolutionary_conservation", "RBP_binding_sites_disrupted")
    disrupted_rbps <- t(sapply(dat_samples$RBP_binding_sites_disrupted, function(x) {
        RBP_counts <- table(strsplit(x, ",")[[1]])
        return(c(paste(paste0(names(RBP_counts),"(",RBP_counts,")"), collapse=","), sum(RBP_counts)))
    }))
    dat_samples$RBP_binding_sites_disrupted <- disrupted_rbps[,1]
    dat_samples <- cbind(dat_samples[,1], as.numeric(sapply(dat_samples[,2], function(x) length(strsplit(x,";")[[1]]))), dat_samples[,2:5], as.numeric(disrupted_rbps[,2]), dat_samples[,6:ncol(dat_samples)]); colnames(dat_samples)[1:7] <- c("sample", "num_variants", "variants", "genes", "genes_with_trans_effect", "evolutionary_conservation", "num_RBP_binding_sites_disrupted")
    dat_samples$genes_with_trans_effect <- unlist(lapply(strsplit(dat_samples$genes_with_trans_effect, ";|\\."), function(x) paste0(x[x!=""], collapse=","))); dat_samples$genes_with_trans_effect[dat_samples$genes_with_trans_effect == ""] <- "."
    dat_samples <- dat_samples[order(dat_samples$num_RBP_binding_sites_disrupted, dat_samples$num_variants, dat_samples$sample, decreasing=TRUE),]
    return(dat_samples)
}
cdh_lan_samples <- get_interesting_samples(cdh_lan)
cdh_lan2_samples <- get_interesting_samples(cdh_lan2)
write.csv(cdh_lan_samples, output_path("cdh1_h3k36me3_rbp_disturbing_samples.csv"), row.names=FALSE)
write.csv(cdh_lan2_samples, output_path("cdh2_h3k36me3_rbp_disturbing_samples.csv"), row.names=FALSE)


# Read coding variants
# Return coding variants data frame from the specified vcf directory.
get_coding_variants <- function(vcf_directory, filter=TRUE) {
    vcf_files <- list.files(vcf_directory)
    vcf_file <- vcf_files[1]
    b <- sapply(1:length(vcf_files), function(vcf_file_num) {
        vcf_file = vcf_files[vcf_file_num]
        print(paste0(vcf_file_num, " ", vcf_file))
        a <- read.csv(gzfile(paste0(vcf_directory, "/", vcf_file)), sep="\t", comment.char = "#", header=FALSE)
        a <- cbind(a[,c(1:2,4:5)], rep(vcf_file, nrow(a)), a[,c(6:ncol(a))])
        for(j in 1:ncol(a)) { a[,j] <- paste0(a[,j]) } #a$Chrom <- paste0(a$Chrom); a$Position <- paste0(a$Position); a$Ref <- paste0(a$Ref); a$Alt <- paste0(a$Alt); a$file <- paste0(a$file)
        colnames(a)[1:5] <- c("Chrom", "Position", "Ref", "Alt", "file")
        a$Chrom <- gsub("chr", "", a$Chrom)
        a <- a[a$Chrom %in% standard_chromosomes,]
        if(filter) {
            # WGS filters
            print("Applying WGS quality filters...")
            parse_info <- function(query, info) { result <- unlist(strsplit(info[which(grepl(paste0("^",query,"="),info))],"=|,")); if(length(result) < 2) { return(NA) }; result <- result[2]; if(is.na(as.numeric(result))) { return(result) } else { return(as.numeric(result)) } }
            passes_wgs_filters <- sapply(1:nrow(a), function(j) {
                info <- strsplit(a$V8[j], ";")[[1]]
                #AC = parse_info("AC",info); FS = parse_info("FS",info); QD = parse_info("QD",info); 
                if(is_indel(a$Ref[j],a$Alt[j])) {
                    ReadPosRankSum = parse_info("ReadPosRankSum",info)
                    return(parse_info("AC",info) < 3 & parse_info("FS",info) <= 25 & parse_info("QD",info) < 1 & parse_info("ReadPosRankSum",info) < -3)
                } else {
                    return(parse_info("AC",info) < 3 & parse_info("FS",info) <= 25 & parse_info("QD",info) < 2)
                }
            }); passes_wgs_filters[is.na(passes_wgs_filters)] <- FALSE
            a <- a[passes_wgs_filters,]
            # GATK filters
            print("Applying GATK filters...")
            parse_gatk_info <- function(dat) {
                gatk_info_names <- unlist(strsplit(paste0(dat$V9),":"))
                proband_gatk_info <- unlist(strsplit(paste0(dat$V10),":")); mother_gatk_info <- unlist(strsplit(paste0(dat$V11),":")); father_gatk_info <- unlist(strsplit(paste0(dat$V12),":"))
                gatk_info <- new.env()
                for(j in 1:length(gatk_info_names)) { gatk_info[[paste0("proband_",gatk_info_names[j])]] <- unlist(strsplit(paste0(proband_gatk_info[j]),",")); gatk_info[[paste0("mother_",gatk_info_names[j])]] <- unlist(strsplit(paste0(mother_gatk_info[j]),",")); gatk_info[[paste0("father_",gatk_info_names[j])]] <- unlist(strsplit(paste0(father_gatk_info[j]),",")) }
                return(gatk_info)
            }
            passes_gatk_filters <- sapply(1:nrow(a), function(j) {
                gatk_info <- parse_gatk_info(a[j,])
                used_info <- c(gatk_info[["proband_PL"]], gatk_info[["proband_AD"]], gatk_info[["mother_DP"]], gatk_info[["father_DP"]], gatk_info[["mother_GQ"]], gatk_info[["father_GQ"]], gatk_info[["mother_AD"]], gatk_info[["father_AD"]])
                #print(used_info)
                if(is.null(used_info) || sum(lengths(used_info) == 0) > 0) { return(NA)
                } else {
                    return(as.numeric(gatk_info[["proband_PL"]][1]) >= 70 & as.numeric(gatk_info[["proband_AD"]][2]) >= 6 & # as.numeric(gatk_info[["proband_AD"]][2])/sum(as.numeric(gatk_info[["proband_AD"]])) >= 0.2 &
                               as.numeric(gatk_info[["mother_DP"]]) >= 10 & as.numeric(gatk_info[["father_DP"]]) >= 10 & as.numeric(gatk_info[["mother_GQ"]]) >= 30 & as.numeric(gatk_info[["father_GQ"]] >= 30)) #&
                    # as.numeric(gatk_info[["mother_AD"]][2])/sum(as.numeric(gatk_info[["mother_AD"]])) <= 0.035 & as.numeric(gatk_info[["father_AD"]][2])/sum(as.numeric(gatk_info[["father_AD"]])) <= 0.035))
                }
            }); passes_gatk_filters[is.na(passes_gatk_filters)] <- FALSE
            a <- a[passes_gatk_filters,]
            print(nrow(a))
        }
        return(a)
    })
    dat <- data.frame()
    for(i in 1:ncol(b)) {
        print(paste0(i, " ", colnames(b)[i]))
        a <- b[,i]
        a <- data.frame(Reduce(cbind, a))
        for(j in 1:ncol(a)) { a[,j] <- paste0(a[,j]) } # a$Chrom <- paste0(a$Chrom); a$Position <- paste0(a$Position); a$Ref <- paste0(a$Ref); a$Alt <- paste0(a$Alt); a$file <- paste0(a$file)
        colnames(a)[1:5] <- c("Chrom", "Position", "Ref", "Alt", "file")
        dat <- rbind(dat, a)
    }
    return(dat)
}
chd1_coding_hg19 <- get_coding_variants(data_path("WGS/coding/CHD1_vcf_denovo_1%_hg19"))
#chd1_pilot_coding_hg19 <- get_coding_variants(data_path("WGS/coding/CHD1_pilot_vcf_denovo_1%_hg19"))
chd2_coding_hg38 <- get_coding_variants(data_path("WGS/coding/CHD2_vcf_denovo_1%coding_hg38"))
ssc_all_coding_hg38 <- get_coding_variants(data_path("WGS/coding/SSC_all_denovo_control_1%coding"))

chd2_coding_hg19 <- liftover(chd2_coding_hg38, from="hg38", to="hg19", chr_colname="Chrom", start_colname="Position", ref_colname="Ref", alt_colname="Alt", mismatches_pause=FALSE, confirm_refseq=FALSE); chd2_coding_hg19 <- chd2_coding_hg19[chd2_coding_hg19$Chrom %in% standard_chromosomes,]
ssc_all_coding_hg19 <- liftover(ssc_all_coding_hg38, from="hg38", to="hg19", chr_colname="Chrom", start_colname="Position", ref_colname="Ref", alt_colname="Alt", mismatches_pause=FALSE, confirm_refseq=FALSE); ssc_all_coding_hg38 <- ssc_all_coding_hg38[ssc_all_coding_hg38$Chrom %in% standard_chromosomes,]


chd1_coding_hg19 <- cbind(chd1_coding_hg19, rep("",nrow(chd1_coding_hg19)), rep("",nrow(chd1_coding_hg19))); colnames(chd1_coding_hg19)[(ncol(chd1_coding_hg19)-1):ncol(chd1_coding_hg19)] <- c("Chrom_hg38", "Position_hg38")
hg19_coding_variants <- combine_datasets(combine_datasets(chd1_coding_hg19, chd2_coding_hg19), ssc_all_coding_hg19)
hg19_coding_variants <- hg19_coding_variants[hg19_coding_variants$Chrom %in% standard_chromosomes,]
write.table(hg19_coding_variants, file=output_path("hg19_coding_variants.tsv"), quote=FALSE, row.names=FALSE, sep="\t")

hg19_coding_snps_annovar <- read.csv(output_path("hg19_coding_variants_annotated2.tsv.snp"), sep="\t")
hg19_coding_snps_features <- read.csv(output_path("hg19_coding_variants_annotated.tsv.snp"), sep="\t")
hg19_coding_snps <- merge(hg19_coding_snps_annovar, hg19_coding_snps_features, by=c("Chrom","Position","Ref","Alt","file",paste0("V",6:12),"Chrom_hg38","Position_hg38"), all.x=TRUE, all.y=TRUE)
hg19_coding_indels_annovar <- read.csv(output_path("hg19_coding_variants_annotated2.tsv.indel"), sep="\t")
hg19_coding_indels_features <- read.csv(output_path("hg19_coding_variants_annotated.tsv.indel"), sep="\t")


all_filters_granges_hg19


chd_burden_hg19[grepl("^RBP$", chdfb_burden_hg19$burden), 1:8]
chdfb_burden_hg19[grepl("^RBP$", chdfb_burden_hg19$burden), 1:8]
best_result_chd_all[grepl("^H3K36me3 in H9 Cells", best_result_chd_all$label), 1:8]
best_result_chdfb_hg19[grepl("^H3K36me3 in H9 Cells", best_result_chdfb_hg19$label), 1:8]
best_result_chd_all[grepl("^RBP H3K36me3 in H9 Cells", best_result_chd_all$label), 1:8]
best_result_chdfb_hg19[grepl("^RBP H3K36me3 in H9 Cells", best_result_chdfb_hg19$label), 1:8]


old_variants_controls_unmapped <- old_variants_controls
old_variants_controls <- sapply(old_variants_controls, function(x) { aba <- strsplit(x,"_")[[1]]; return(paste0(sapply(aba[1], function(sample) { mapped_sample <- ssc_id_conversion$V1[ssc_id_conversion$V2 == sample]; if(length(mapped_sample)>0) { return(mapped_sample[1]) } else { return(sample) } }),"_",aba[2])) })
names(old_variants_controls) <- c()
table(unlist(lapply(strsplit(new_variants_controls[!(new_variants_controls %in% old_variants_controls)],"_"), function(x) x[1])))

plot(as.numeric(paste0(table(chd$sample)))[names(table(chd$sample)) %in% unique(chdfb$sample)], as.numeric(paste0(table(chdfb$sample)))[names(table(chdfb$sample)) %in% unique(chd$sample)], xlab="GATK-only", ylab="GATK+FB", main="CHD Sample Variant Counts")
plot(as.numeric(paste0(table(ssc_518$sample)))[names(table(ssc_518$sample)) %in% unique(sscfb$sample)], as.numeric(paste0(table(sscfb$sample)))[names(table(sscfb$sample)) %in% unique(ssc_518$sample)], xlab="GATK-only", ylab="GATK+FB", main="SSC Sample Variant Counts")
plot(as.numeric(paste0(table(chd2$sample)))[names(table(chd2$sample)) %in% unique(chdfb2$sample)], as.numeric(paste0(table(chdfb2$sample)))[names(table(chdfb2$sample)) %in% unique(chd2$sample)], xlab="GATK-only", ylab="GATK+FB", main="CHD2 Sample Variant Counts")
plot(as.numeric(paste0(table(ssc_1088$sample)))[names(table(ssc_1088$sample)) %in% unique(sscfb2$sample)], as.numeric(paste0(table(sscfb2$sample)))[names(table(sscfb2$sample)) %in% unique(ssc_1088$sample)], xlab="GATK-only", ylab="GATK+FB", main="SSC2 Sample Variant Counts")
plot(as.numeric(paste0(table(ssc_1088$sample[!is.na(ssc_1088_muts_wgsa$RBP)])))[names(table(ssc_1088$sample[!is.na(ssc_1088_muts_wgsa$RBP)])) %in% unique(sscfb2$sample[!is.na(sscfb2_muts_wgsa$RBP)])], as.numeric(paste0(table(sscfb2$sample[!is.na(sscfb2_muts_wgsa$RBP)])))[names(table(sscfb2$sample[!is.na(sscfb2_muts_wgsa$RBP)])) %in% unique(ssc_1088$sample[!is.na(ssc_1088_muts_wgsa$RBP)])], xlab="GATK-only", ylab="GATK+FB", main="SSC2 Sample RBP Variant Counts")


all_variants_wgsa <- combine_datasets(cdh_combined_muts_wgsa, ssc_all_muts_wgsa)
all_variants_wgsa <- cbind(rep(0, nrow(all_variants_wgsa)), all_variants_wgsa); colnames(all_variants_wgsa)[1] <- "CDH"; all_variants_wgsa$CDH[1:nrow(cdh_combined_muts_wgsa)] <- 1 
all_variants_wgsa <- all_variants_wgsa[,colnames(all_variants_wgsa) %in% c("CDH","phastCons46way_placental","phastCons100way_vertebrate","CADD_phred","Eigen.pred","FANTOM5_enhancer_robust","funseq2_noncoding_rankscore",colnames(all_variants_wgsa)[grepl("gappedPeak",colnames(all_variants_wgsa))])]
#all_variants_wgsa[is.na(all_variants_wgsa)] <- 0
for(i in 1:ncol(all_variants_wgsa)) {
    print(i)
    all_variants_wgsa[,i] <- paste0(all_variants_wgsa[,i])
    all_variants_wgsa[,i][is.na(all_variants_wgsa[,i])] <- 0 
    all_variants_wgsa[,i][all_variants_wgsa[,i] == "N"] <- 0; all_variants_wgsa[,i][all_variants_wgsa[,i] == "Y"] <- 1
    all_variants_wgsa[,i] <- as.numeric(all_variants_wgsa[,i])
    all_variants_wgsa[,i][is.na(all_variants_wgsa[,i])] <- 0
}
all_variants_wgsa[1:4,]

all_variantsfb_wgsa <- combine_datasets(chdfb_combined_muts_wgsa, sscfb_combined_muts_wgsa)
all_variantsfb_wgsa <- cbind(rep(0, nrow(all_variantsfb_wgsa)), all_variantsfb_wgsa); colnames(all_variantsfb_wgsa)[1] <- "CHD"; all_variantsfb_wgsa$CDH[1:nrow(cdh_combined_muts_wgsa)] <- 1 
all_variantsfb_wgsa <- all_variantsfb_wgsa[,colnames(all_variantsfb_wgsa) %in% c("CHD","phastCons46way_placental","phastCons100way_vertebrate","CADD_phred","Eigen.pred","FANTOM5_enhancer_robust","funseq2_noncoding_rankscore",colnames(all_variantsfb_wgsa)[grepl("gappedPeak",colnames(all_variantsfb_wgsa))])]
#all_variantsfb_wgsa[is.na(all_variantsfb_wgsa)] <- 0
for(i in 1:ncol(all_variantsfb_wgsa)) {
    print(i)
    all_variantsfb_wgsa[,i] <- paste0(all_variantsfb_wgsa[,i])
    all_variantsfb_wgsa[,i][is.na(all_variantsfb_wgsa[,i])] <- 0 
    all_variantsfb_wgsa[,i][all_variantsfb_wgsa[,i] == "N"] <- 0; all_variantsfb_wgsa[,i][all_variantsfb_wgsa[,i] == "Y"] <- 1
    all_variantsfb_wgsa[,i] <- as.numeric(all_variantsfb_wgsa[,i])
    all_variantsfb_wgsa[,i][is.na(all_variantsfb_wgsa[,i])] <- 0
}
all_variantsfb_wgsa[1:4,]

library("autoencoder")
## Set up the autoencoder architecture:
nl=3 # number of layers (default is 3: input, hidden, output)
unit.type = "logistic" # specify the network unit type, i.e., the unit's activation function ("logistic" or "tanh")
Nx.patch=10 # width of training image patches, in pixels
Ny.patch=10 # height of training image patches, in pixels
N.input = Nx.patch*Ny.patch # number of units (neurons) in the input layer (one unit per pixel)
N.hidden = 5*5 # number of units in the hidden layer
lambda = 0.0002 # weight decay parameter
beta = 6 # weight of sparsity penalty term
rho = 0.01 # desired sparsity parameter
epsilon <- 0.001 # a small parameter for initialization of weights as small gaussian random numbers sampled from N(0,epsilon^2)
max.iterations = 2000 ## number of iterations in optimizer

a <- autoencode(as.matrix(all_variants_wgsa[,-c(1:5)]), X.test = NULL, nl=nl, N.hidden=2, unit.type=unit.type,
                lambda=lambda, beta=beta, rho=rho, epsilon=epsilon, optim.method="BFGS", #optim.method can be one of "BFGS", "L-BFGS-B", or "CG".
                rel.tol=sqrt(.Machine$double.eps), max.iterations=max.iterations,
                rescale.flag=T, rescaling.offset = 0.001)
a_feature_cols <- rep("red", ncol(all_variants_wgsa[,-c(1:5)]));
plot(a$W[[1]][1,], a$W[[1]][2,])
visualize.hidden.units(a,1, 238)

colnames(all_variants_wgsa[,-c(1:5)])[order(a$W[[1]][1,])]


library("biomartr") 
library("PWMEnrich")
library("PWMEnrich.Hsapiens.background")
library("TFBSTools")
library("seqLogo")
library("msa")
#library("ggplot2")

# load the pre-compiled lognormal background
data(PWMLogn.hg19.MotifDb.Hsap)
# load the stripe2 sequences from a FASTA file for motif enrichment
#sequence = DNAString("TGCATCAAGTGTGTAGTGCAAGTGAGTGATGAGTAGAAGTTGAGTGAGGTAGATGC") # readDNAStringSet(system.file(package="PWMEnrich", dir="extdata", file="stripe2.fa"))
#refseq <- get_refseq("1", "hg19") #cat(paste0("Reading chr", chromosome, "...")); #get_refseq(chromosome, version); #cat("Done.\nProcessing positions...\n")sequence

#a <- BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg19::Hsapiens, "chr1", 10000, 20000, as.character = T)
#hg <- read_genome(file=data_path("HS.genome.38/Homo_sapiens_genomic_refseq.fna.gz"), obj.type="Biostrings")


get_aligned_sequences <- function(gr, min_sequence_length=NULL, max_sequence_length=NULL, bucket_size=1000, bucket=1, shuffle=FALSE, filename=NULL) {
    if(shuffle) { gr <- sample(gr) }
    if(!is.null(min_sequence_length)) { gr <- gr[width(gr) >= min_sequence_length] }
    sequence_lengths <- width(gr)
    if(!is.null(max_sequence_length)) { sequence_lengths <- sapply(sequence_lengths, function(sequence_length) max(c(sequence_length, max_sequence_length))) }
    padding_length = max(sequence_lengths) * 2
    gr_sequences <- granges_to_DNAStringSet(gr, length=padding_length*2+sequence_lengths)
    
    num_buckets = ceiling(length(gr_sequences)/bucket_size)
    num_buckets = 1
    gr_sequences <- sapply(bucket, function(bucket) { #sapply(1:num_buckets, function(bucket) {
        bucket_start_index = ((bucket-1)*bucket_size+1); bucket_end_index = min(c(bucket*bucket_size,length(gr_sequences)))
        print(paste0(bucket_start_index,":",bucket_end_index))
        padded_sequences <- gr_sequences[bucket_start_index:bucket_end_index]
        padded_sequence_lengths <- width(padded_sequences)
        sequences_left_padding <- subseq(padded_sequences,start=1, end=padding_length)
        sequences <- subseq(padded_sequences,start=padding_length+1, end=padded_sequence_lengths-padding_length)
        sequences_right_padding <- subseq(padded_sequences,start=padded_sequence_lengths-padding_length+1, end=padded_sequence_lengths)
        
        sequences_msa <- msa(sequences, "ClustalOmega", cluster=8, order="input")
        m_consensus <- msaConsensusSequence(sequences_msa)
        print(m_consensus)
        m <- sequences_msa@unmasked
        num_start_gaps = width(m) - nchar(gsub("^-*","",m))
        num_end_gaps = width(m) - nchar(gsub("-*$","",m))
        
        #at <- matrix(c(TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE), nrow=length(probes), ncol=width(probes)[1], byrow=TRUE)
        
        # Fill in gaps at the starts of the sequences
        print("Filling in gaps at the starts of sequences...")
        m_filled <- DNAStringSet(sapply(1:length(m), function(i) { m_sequence <- unlist(m[i]); if(num_start_gaps[i] > 0) { replaceLetterAt(m_sequence, c(1:num_start_gaps[i]), subseq(unlist(sequences_left_padding[i]), start=padding_length-num_start_gaps[i]+1, end=padding_length)) } else { m_sequence } }))
        # Fill in gaps at the ends of the sequences
        print("Filling in gaps at the ends of sequences...")
        m_filled <- DNAStringSet(sapply(1:length(m_filled), function(i) { m_sequence <- unlist(m_filled[i]); if(num_end_gaps[i] > 0) { replaceLetterAt(m_sequence, c((length(m_sequence)-num_end_gaps[i]+1):length(m_sequence)), subseq(unlist(sequences_right_padding[i]), start=1, end=num_end_gaps[i])) } else { m_sequence } }))
        
        if(!is.null(filename)) {
            gr_subset <- gr[bucket_start_index:bucket_end_index] 
            log_table <- data.frame(paste0("chr",seqnames(gr_subset)), start(gr_subset), end(gr_subset), width(gr_subset), strand(gr_subset), names(gr_subset), paste0(m), paste0(m_filled), width(m_filled))
            colnames(log_table) <-  c("chrom", "start", "end", "width", "strand", "p.value", "msa", "msa_filled", "msa_width")
            write.csv(log_table, file=filename, row.names=FALSE)
        }
        
        return(m_filled)
    })[[1]]
    
    if(!is.null(filename)) {
        
        
    }
    
    return(gr_sequences)
}

pwm_analysis <- function(rbp, motif_length=8, min_sequence_length=NULL, max_sequence_length=NULL, padding=0, bucket_size=1000, bucket=1, composite=FALSE) {
    all_rbp_features <- get_features_by_group("RBP_50bp")
    rbp_features <- all_rbp_features[grepl(rbp, all_rbp_features)]
    rbp_granges <- sapply(rbp_features, function(rbp_feature) { load_annotation(rbp_feature, padding=padding) })
    rbp_granges <- Reduce(c, rbp_granges)
    rbp_granges <- rbp_granges[order(as.numeric(names(rbp_granges)))] # order by p.value/significance of peak
    
    rbp_sequences <- get_aligned_sequences(rbp_granges, min_sequence_length=min_sequence_length, max_sequence_length=max_sequence_length, bucket_size=bucket_size, bucket=bucket, shuffle=FALSE, filename=output_path(paste0(rbp,"_eclip_msa",bucket,".csv")))
    rbp_sequences <- rbp_sequences[!grepl("[^ATCG]", rbp_sequences)]
    if(composite) { rbp_sequences <- sequence_composite(rbp_sequences) }
    rbp_pwm <- TFBSTools::toPWM(rbp_sequences, type="prob", pseudocounts=0.8,  bg=genomic.acgt)
    
    rbp_seqlogo_pwm <- makePWM(rbp_pwm)
    rbp_seqlogo_pwm@consensus; rbp_seqlogo_pwm@ic
    
    motif_start_index = which.max(rollapply(rbp_seqlogo_pwm@ic, width=motif_length, mean))
    motif_end_index = motif_start_index + motif_length - 1
    trimmed_rbp_pwm <- rbp_pwm[,motif_start_index:motif_end_index]
    
    seqlogo_title <- paste0(rbp," PWM logo"); if(composite) { seqlogo_title <- paste0(seqlogo_title," composite") }
    
    seqLogo(trimmed_rbp_pwm, ic.scale=TRUE)
    grid.text(seqlogo_title, x=0.5, y=0.95, gp=gpar(fontsize=18, fontface="bold", col="black"))
    dev.copy2pdf(file=output_path(paste0(gsub(" ", "_", tolower(seqlogo_title)),bucket)))
    
    seqlogo_title = gsub("logo", "extended logo", seqlogo_title)
    seqLogo(rbp_pwm, ic.scale=TRUE)
    grid.text(seqlogo_title, x=0.5, y=0.95, gp=gpar(fontsize=18, fontface="bold", col="black"))
    dev.copy2pdf(file=output_path(paste0(gsub(" ", "_", tolower(seqlogo_title)),bucket)))
    
    return(trimmed_rbp_pwm)
}
sapply(1:3, function(j) pwm_analysis("RBFOX2", padding=50, bucket_size=5000, bucket=j)) # motif should be TGCATG
sum(grepl("TGCATG", paste0(sequences)))
sum(grepl("TGCATC", paste0(sequences)))
grepl("TGCATG", paste0(sequences))[1:10]
bases <- c("A","C","G","T")
k = 5
all_kmers <- apply(expand.grid(lapply(1:k, function(i) return(bases))), 1, paste0, collapse="")

all_kmer_exact_searches <- lapply(1:length(all_kmers), function(i) { print(i); candidate_motif <- all_kmers[i]; sequences[grepl(candidate_motif, paste0(sequences))] })
all_kmer_fuzzy_searches <- lapply(1:length(all_kmers), function(i) { print(i); candidate_motif <- all_kmers[i]; sequences[agrepl(candidate_motif, paste0(sequences), max.distance=list(insertions=0, deletions=0, substitutions=1))] })
all_kmer_fuzzy_searches_lengths <- sapply(all_kmer_fuzzy_searches, FUN=length)
plot(density(all_kmer_fuzzy_searches_lengths))
all_kmers <- all_kmers[order(all_kmer_fuzzy_searches_lengths, decreasing=TRUE)]
all_kmer_fuzzy_searches <- all_kmer_fuzzy_searches[order(all_kmer_fuzzy_searches_lengths, decreasing=TRUE)]
msa(unique(Reduce(c, all_kmer_fuzzy_searches[1:10])), "ClustalOmega")



sapply(1:3, function(j) pwm_analysis("QKI", bucket_size=1000, bucket=j)) # motif should be ACTAAYN{1-20}TAAY

library("DescTools")
Entropy(as.matrix(rep(1/8, 8)))


#1 - sum(rbfox2_sequences %in% qki_sequences)/length(rbfox2_sequences)
#1 - sum(qki_sequences %in% rbfox2_sequences)/length(qki_sequences)
#genomic.acgt <- getBackgroundFrequencies("hg19")

rbfox2_granges <- c(load_annotation("HepG2.RBFOX2"), load_annotation("K562.RBFOX2"))
rbfox2_sequences <- get_aligned_sequences(rbfox2_granges, sequence_length=7, bucket_size=5000)
rbfox2_sequences <- sequence_composite(rbfox2_sequences[!grepl("[^ATCG]", rbfox2_sequences)])       
rbfox2_pwm <- TFBSTools::toPWM(rbfox2_sequences[!(rbfox2_sequences %in% qki_sequences)], type="prob", pseudocounts=0.8,  bg=genomic.acgt)
seqLogo(rbfox2_pwm, ic.scale=TRUE)
makePWM(rbfox2_pwm)@consensus; makePWM(rbfox2_pwm)@ic

motif_length = 6
motif_start_index = which.max(rollapply(makePWM(rbfox2_pwm)@ic, width=motif_length, mean))
motif_end_index = motif_start_index + motif_length - 1
trimmed_rbfox2_pwm <- rbfox2_pwm[,motif_start_index:motif_end_index]
seqLogo(trimmed_rbfox2_pwm, ic.scale=TRUE)


qki_granges <- c(load_annotation("HepG2.QKI"), load_annotation("K562.QKI"))
qki_sequences <- get_aligned_sequences(qki_granges, sequence_length=7, bucket_size=5000)
qki_sequences <- qki_sequences[!grepl("[^ATCG]", qki_sequences)]
qki_pwm <- TFBSTools::toPWM(qki_sequences[!(qki_sequences %in% rbfox2_sequences)], type="prob", pseudocounts=0.8,  bg=genomic.acgt)
seqLogo(qki_pwm, ic.scale=TRUE)
makePWM(qki_pwm)@consensus; makePWM(qki_pwm)@ic

motif_length = 6
motif_start_index = which.max(rollapply(makePWM(qki_pwm)@ic, width=motif_length, mean))
motif_end_index = motif_start_index + motif_length - 1
trimmed_qki_pwm <- qki_pwm[,motif_start_index:motif_end_index]
seqLogo(trimmed_qki_pwm, ic.scale=TRUE)


rbfox2_m2 <- apply(floor(rbfox2_pwm*100), c(1,2), function (x) { return(as.integer(x)) })
qki_m2 <- apply(floor(qki_pwm*100), c(1,2), function (x) { return(as.integer(x)) })
motifs <- list(RBFOX2=rbfox2_m2, QKI=qki_m2)
motifs_bg = makeBackground(motifs, organism="hg19", type="logn", quick=TRUE)

# perform TSS regions motif enrichment!
TSS_sequences <- granges_to_DNAStringSet(TSS_granges)
TSS_motif_enrichment = motifEnrichment(TSS_sequences, motifs_bg)
report = sequenceReport(TSS_motif_enrichment, 1)
report
groupReport(TSS_motif_enrichment)

all_tss_sequences_report <- data.frame(t(sapply(1:length(TSS_sequences), function(i) { return(as.data.frame(sequenceReport(TSS_motif_enrichment, i))) })))
all_tss_sequences_report$p.value <- as.numeric(paste0(all_tss_sequences_report$p.value)); all_tss_sequences_report$raw.score <- as.numeric(paste0(all_tss_sequences_report$raw.score))
all_tss_sequences_report[all_tss_sequences_report$p.value < 0.05,]
nrow(all_tss_sequences_report[all_tss_sequences_report$p.value < 0.05,])

# perform TES/3'UTR regions motif enrichment!
TES_sequences <- granges_to_DNAStringSet(TES_granges)
TES_motif_enrichment = motifEnrichment(TES_sequences, motifs_bg)
report = sequenceReport(TES_motif_enrichment, 1)
report
groupReport(TES_motif_enrichment)

all_tes_sequences_report <- data.frame(t(sapply(1:length(TES_sequences), function(i) { return(as.data.frame(sequenceReport(TES_motif_enrichment, i))) })))
all_tes_sequences_report$p.value <- as.numeric(paste0(all_tes_sequences_report$p.value)); all_tes_sequences_report$raw.score <- as.numeric(paste0(all_tes_sequences_report$raw.score))
all_tes_sequences_report[all_tes_sequences_report$p.value < 0.05,]
nrow(all_tes_sequences_report[all_tes_sequences_report$p.value < 0.05,])

# perform case vs. control regions motif enrichment!
# Add sequence context features:
add_sequence_context_feature <- function(dat, chr_colname="Chrom", pos_colname="Position", alt_colname="Alt", width=51, version="hg19") {
    #sequences <- get_trimers(dat[,chr_colname], dat[,pos_colname], width=width, allow_BSgenome=TRUE, version=version)
    arm_width = floor(width/2)
    all_starts <- as.numeric(paste0(dat[,pos_colname])) - arm_width; all_ends <- as.numeric(paste0(dat[,pos_colname])) + arm_width
    all_chromosomes <- paste0(dat[,chr_colname])
    all_alts <- paste0(dat[,alt_colname])
    
    sequences <- lapply(unique(all_chromosomes), function(chromosome) {
        if(grepl("Y",chromosome)) { return(matrix(nrow=0, ncol=ncol(dat)+2)) }
        cat(paste0("Reading chr", gsub("chr","",chromosome), "..."))
        refseq <- get_refseq(chromosome, version="hg19", allow_BSgenome=TRUE)[[1]]
        cat("Done.\nGrabbing sequences...")
        curr_chrom_indices <- which(all_chromosomes == chromosome)
        starts <- all_starts[curr_chrom_indices]; ends <- all_ends[curr_chrom_indices]; widths <- ends - starts + 1
        split_indices <- cumsum(widths)
        curr_chrom_rbp_sequences <- eval(parse(text=paste0("substring(paste0(refseq[c(",paste0(starts,":",ends, collapse=","),")]), c(0,split_indices[-length(split_indices)])+1, split_indices)")))
        cat("Done.\n")
        
        mutated_sequences <- curr_chrom_rbp_sequences
        substr(mutated_sequences, arm_width+1, arm_width+1) <- all_alts[curr_chrom_indices]
        return(cbind(dat[curr_chrom_indices,], c(curr_chrom_rbp_sequences), c(mutated_sequences)))
    })
    sequences <- data.frame(rbindlist(lapply(sequences, function(x) return(data.frame(x)))))
    colnames(sequences)[(ncol(sequences)-1):ncol(sequences)] <- c("ref_sequence", "alt_sequence")
    sequences$ref_sequence <- paste0(sequences$ref_sequence); sequences$alt_sequence <- paste0(sequences$alt_sequence)
    
    return(sequences)
}
chd_variant_dat <- add_sequence_context_feature(chdfb_combined_muts_wgsa, chr_colname="X.chr", pos_colname="pos", alt_colname="alt", width=51)
ssc_variant_dat <- add_sequence_context_feature(sscfb_combined_muts_wgsa, chr_colname="X.chr", pos_colname="pos", alt_colname="alt", width=51)

chd_ref_sequences <- DNAStringSet(chd_variant_dat$ref_sequence); chd_alt_sequences <- DNAStringSet(chd_variant_dat$alt_sequence) 
ssc_ref_sequences <- DNAStringSet(ssc_variant_dat$ref_sequence); ssc_alt_sequences <- DNAStringSet(ssc_variant_dat$alt_sequence) 

chd_ref_motif_enrichment = motifEnrichment(chd_ref_sequences, motifs_bg)
groupReport(chd_ref_motif_enrichment)

plot(sequenceReport(chd_ref_motif_enrichment, 1), fontsize=7, id.fontsize=6)

chd_alt_motif_enrichment = motifEnrichment(chd_alt_sequences, motifs_bg)
groupReport(chd_alt_motif_enrichment)
ssc_ref_motif_enrichment = motifEnrichment(ssc_ref_sequences, motifs_bg)
groupReport(ssc_ref_motif_enrichment)
ssc_alt_motif_enrichment = motifEnrichment(ssc_alt_sequences, motifs_bg)
groupReport(ssc_alt_motif_enrichment)

#length(chd_ref_sequences)
chd_sequences_report_rbfox2 <- data.frame(t(sapply(1:length(chd_ref_sequences), function(i) { 
    print(i)
    ref_report <- as.data.frame(sequenceReport(chd_ref_motif_enrichment, i))[,c("target", "raw.score", "p.value")]
    alt_report <- as.data.frame(sequenceReport(chd_alt_motif_enrichment, i))[,c("target", "raw.score", "p.value")]
    report <- merge(ref_report, alt_report, by=c("target"))
    colnames(report) <- gsub("\\.x", "_ref", colnames(report)); colnames(report) <- gsub("\\.y", "_alt", colnames(report))
    return(report[report$target == "RBFOX2",])
})))
chd_sequences_report_rbfox2$raw.score_ref <- as.numeric(paste0(chd_sequences_report_rbfox2$raw.score_ref)); chd_sequences_report_rbfox2$raw.score_alt <- as.numeric(paste0(chd_sequences_report_rbfox2$raw.score_alt))

chd_sequences_report_qki <- data.frame(t(sapply(1:length(chd_ref_sequences), function(i) { 
    print(i)
    ref_report <- as.data.frame(sequenceReport(chd_ref_motif_enrichment, i))[,c("target", "raw.score", "p.value")]
    alt_report <- as.data.frame(sequenceReport(chd_alt_motif_enrichment, i))[,c("target", "raw.score", "p.value")]
    report <- merge(ref_report, alt_report, by=c("target"))
    colnames(report) <- gsub("\\.x", "_ref", colnames(report)); colnames(report) <- gsub("\\.y", "_alt", colnames(report))
    return(report[report$target == "QKI",])
})))
chd_sequences_report_qki$raw.score_ref <- as.numeric(paste0(chd_sequences_report_qki$raw.score_ref)); chd_sequences_report_qki$raw.score_alt <- as.numeric(paste0(chd_sequences_report_qki$raw.score_alt))

ssc_sequences_report_rbfox2 <- data.frame(t(sapply(1:length(ssc_ref_sequences), function(i) { 
    print(i)
    ref_report <- as.data.frame(sequenceReport(ssc_ref_motif_enrichment, i))[,c("target", "raw.score", "p.value")]
    alt_report <- as.data.frame(sequenceReport(ssc_alt_motif_enrichment, i))[,c("target", "raw.score", "p.value")]
    report <- merge(ref_report, alt_report, by=c("target"))
    colnames(report) <- gsub("\\.x", "_ref", colnames(report)); colnames(report) <- gsub("\\.y", "_alt", colnames(report))
    return(report[report$target == "RBFOX2",])
})))
ssc_sequences_report_rbfox2$raw.score_ref <- as.numeric(paste0(ssc_sequences_report_rbfox2$raw.score_ref)); ssc_sequences_report_rbfox2$raw.score_alt <- as.numeric(paste0(ssc_sequences_report_rbfox2$raw.score_alt))

ssc_sequences_report_qki <- data.frame(t(sapply(1:length(chd_ref_sequences), function(i) { 
    print(i)
    ref_report <- as.data.frame(sequenceReport(chd_ref_motif_enrichment, i))[,c("target", "raw.score", "p.value")]
    alt_report <- as.data.frame(sequenceReport(chd_alt_motif_enrichment, i))[,c("target", "raw.score", "p.value")]
    report <- merge(ref_report, alt_report, by=c("target"))
    colnames(report) <- gsub("\\.x", "_ref", colnames(report)); colnames(report) <- gsub("\\.y", "_alt", colnames(report))
    return(report[report$target == "QKI",])
})))
chd_sequences_report_qki$raw.score_ref <- as.numeric(paste0(chd_sequences_report_qki$raw.score_ref)); chd_sequences_report_qki$raw.score_alt <- as.numeric(paste0(chd_sequences_report_qki$raw.score_alt))



plot(log(chd_sequences_report_rbfox2$raw.score_ref), log(chd_sequences_report_rbfox2$raw.score_alt), col="red", main="RBFOX2 Alt vs Ref score")
lines(c(-10000,10000), c(-10000,10000))
dev.copy2pdf(file=output_path("RBFOX2_alt_vs_ref_score.pdf"))
#points(log(chd_sequences_report_qki$raw.score_ref), log(chd_sequences_report_qki$raw.score_alt), col="blue")
#legend("topleft", legend=c("RBFOX2", "QKI"), col=c("red", "blue"), pch=15)

sum(chd_sequences_report_rbfox2$raw.score_alt > chd_sequences_report_rbfox2$raw.score_ref)/nrow(chd_sequences_report_rbfox2)
sum(chd_sequences_report_qki$raw.score_alt > chd_sequences_report_qki$raw.score_ref)/nrow(chd_sequences_report_qki)
sum(ssc_sequences_report_rbfox2$raw.score_alt > ssc_sequences_report_rbfox2$raw.score_ref)/nrow(ssc_sequences_report_rbfox2)
sum(ssc_sequences_report_qki$raw.score_alt > ssc_sequences_report_qki$raw.score_ref)/nrow(ssc_sequences_report_qki)

plot(density(log(chd_sequences_report_rbfox2$raw.score_ref)), col="blue", ylim=c(0,0.16), main="PWM Score Densities", xlab="log(raw_pwm_enrichment_score)")
lines(density(log(chd_sequences_report_rbfox2$raw.score_alt)), col="green")
lines(density(log(chd_sequences_report_qki$raw.score_ref)), col="red")
lines(density(log(chd_sequences_report_qki$raw.score_alt)), col="orange")
legend("topleft", legend=c("RBFOX2 Ref", "RBFOX2 Alt", "QKI Ref", "QKI Alt"), col=c("blue", "green", "red", "orange"), pch=15)
mtext("51 bp motif lengths")
dev.copy2pdf(file=output_path("pwm_score_densities.pdf"))

plot(density(chd_sequences_report_rbfox2$raw.score_alt / chd_sequences_report_rbfox2$raw.score_ref), col="red", xlab="RBFOX2 alt_score / ref_score", main="")
dev.copy2pdf(file=output_path("disruption_density.pdf"))
#lines(density(chd_sequences_report_qki$raw.score_alt / chd_sequences_report_qki$raw.score_ref), col="blue")
legend("topright", legend=c("RBFOX2", "QKI"), col=c("red", "blue"), pch=15)

plot(log(chd_sequences_report_rbfox2$raw.score_ref), log(chd_sequences_report_rbfox2$raw.score_alt / chd_sequences_report_rbfox2$raw.score_ref), col="red", xlab="log(ref_score)", ylab="log(alt_score / ref_score)", main="RBFOX2")
dev.copy2pdf(file=output_path("RBFOX2_delta_vs_ref_score.pdf"))

all_tes_sequences_report$p.value <- as.numeric(paste0(all_tes_sequences_report$p.value)); all_tes_sequences_report$raw.score <- as.numeric(paste0(all_tes_sequences_report$raw.score))
all_tes_sequences_report[all_tes_sequences_report$p.value < 0.05,]
nrow(all_tes_sequences_report[all_tes_sequences_report$p.value < 0.05,])


sequences <- granges_to_DNAStringSet(TSS_granges[1:50])
motif_enrichment = motifEnrichment(sequences, rbfox2_pwm)

# perform motif enrichment!
sequences <- granges_to_DNAStringSet(TSS_granges[1:50]) #DNAStringSet(list("chr1"=BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg19::Hsapiens, "chr1", 10000, 10050), "chr1"=BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg19::Hsapiens, "chr1", 20000, 20050)))
motif_enrichment = motifEnrichment(sequences, PWMLogn.hg19.MotifDb.Hsap)
## Calculating motif enrichment scores ...
report = sequenceReport(motif_enrichment, 1)[c(which(sequenceReport(motif_enrichment, 1)$target == "GATA4"), 1:20)]
report
sequenceReport(motif_enrichment, 2)[c(which(sequenceReport(motif_enrichment, 2)$target == "GATA4"), 1:20)]
# plot the motif with P-value < 0.05
plot(report[report$target == "GATA4" | report$p.value < 0.05], fontsize=7, id.fontsize=6)


# Analyze motif disturbance in CDH TSS and 3'UTR.
analyze_motif_disturbance <- function(dat, proximity=12, filename=NULL) {
    ref_sequences <- get_trimers(dat$X.chr, dat$pos, width=(2*proximity+1), allow_BSgenome=TRUE)$trimers
    mut_sequences <- sapply(1:length(ref_sequences), function(i) { seq <- ref_sequences[i]; substr(seq, proximity+1, proximity+1) <- paste0(cdh_muts_TES$alt[i]); return(seq) } )
    ref_sequences <- sequences_to_DNAStringSet(ref_sequences, dat$X.chr)
    mut_sequences <- sequences_to_DNAStringSet(mut_sequences, dat$X.chr)
    ref_motif_enrichment = motifEnrichment(ref_sequences, PWMLogn.hg19.MotifDb.Hsap)
    mut_motif_enrichment = motifEnrichment(mut_sequences, PWMLogn.hg19.MotifDb.Hsap)
    motif_disturbance_dats <- new.env(); motif_disturbance_dats_group_size = 200; num_motif_disturbance_groups <- ceiling(length(ref_sequences)/motif_disturbance_dats_group_size)
    for(i in 1:num_motif_disturbance_groups) { motif_disturbance_dats[[paste0(i)]] <- data.frame() }
    motif_disturbance_per_variant_dat <- dat[,c("X.chr","pos","ref","alt","sample","TS_gene")]; colnames(motif_disturbance_per_variant_dat)[which(colnames(motif_disturbance_per_variant_dat) == "TS_gene")] <- "gene"
    motif_disturbances <- sapply(1:length(ref_sequences), function(i) {
        print(paste0(i, " / ", length(ref_sequences)))
        ref_report <- sequenceReport(ref_motif_enrichment, i)#[c(which(sequenceReport(ref_motif_enrichment, i)$target == "GATA4"), 1:20)]
        mut_report <- sequenceReport(mut_motif_enrichment, i)#[c(which(sequenceReport(mut_motif_enrichment, i)$target == "GATA4"), 1:20)]
        ref_report <- cbind(ref_report$rank, ref_report$target, ref_report$id, ref_report$raw.score, ref_report$p.value); colnames(ref_report) <- c("rank", "target", "id", "raw.score", "p.value")
        mut_report <- cbind(mut_report$rank, mut_report$target, mut_report$id, mut_report$raw.score, mut_report$p.value); colnames(mut_report) <- c("rank", "target", "id", "raw.score", "p.value")
        report <- merge(ref_report, mut_report, by=c("target","id"))
        report$p.value.x <- as.numeric(paste0(report$p.value.x)); report$p.value.y <- as.numeric(paste0(report$p.value.y))
        harmed_motif_bindings <- report[report$p.value.x < 0.05 & report$p.value.y > 0.05 & !grepl("UW\\.Motif",report$target),]
        if(nrow(harmed_motif_bindings) == 0) { return(".") }
        harmed_motif_bindings <- harmed_motif_bindings[!duplicated(harmed_motif_bindings$target),]
        motif_disturbance_dat_group_index = paste0(ceiling(i/motif_disturbance_dats_group_size))
        motif_disturbance_dats[[motif_disturbance_dat_group_index]] <- rbind(motif_disturbance_dats[[motif_disturbance_dat_group_index]], cbind(motif_disturbance_per_variant_dat[i,], harmed_motif_bindings))
        return(paste0(paste0(harmed_motif_bindings$target," (",formatC(harmed_motif_bindings$p.value.x,format="e",digits=2)," -> ",formatC(harmed_motif_bindings$p.value.y,format="e",digits=2),")"), collapse=", "))
    })
    if(!is.null(filename)) {
        motif_disturbance_per_variant_dat <- cbind(motif_disturbance_per_variant_dat, motif_disturbances)
        write.csv(motif_disturbance_per_variant_dat, file=filename, row.names=FALSE)
    }
    motif_disturbance_dat <- data.frame()
    for(i in 1:num_motif_disturbance_groups) { 
        print(paste0(i, " / ", num_motif_disturbance_groups))
        motif_disturbance_dat <- rbind(motif_disturbance_dat, motif_disturbance_dats[[paste0(i)]]) 
    }
    
    
    return(motif_disturbance_dat)
}
cdh_tes_motif_disturbance <- analyze_motif_disturbance(cdh_combined_muts_TES, filename=output_path("cdh_tes_motif_disturbance_per_variant.csv"))
write.csv(cdh_tes_motif_disturbance, file=output_path("cdh_tes_motif_disturbance.csv"), row.names=FALSE)
cdh_tss_motif_disturbance <- analyze_motif_disturbance(cdh_combined_muts_TSS, filename=output_path("cdh_tss_motif_disturbance_per_variant.csv"))
write.csv(cdh_tss_motif_disturbance, file=output_path("cdh_tss_motif_disturbance.csv"), row.names=FALSE)

source("alex_suite.R")
library("biomartr"); library("PWMEnrich"); library("PWMEnrich.Hsapiens.background")
data(PWMLogn.hg19.MotifDb.Hsap)

ssc1_tes_motif_disturbance <- analyze_motif_disturbance(ssc_518_muts_TES, filename=output_path("ssc1_tes_motif_disturbance_per_variant.csv"))
write.csv(ssc1_tes_motif_disturbance, file=output_path("ssc1_tes_motif_disturbance.csv"), row.names=FALSE)
ssc1_tss_motif_disturbance <- analyze_motif_disturbance(ssc_518_muts_TSS, filename=output_path("ssc1_tss_motif_disturbance_per_variant.csv"))
write.csv(ssc1_tss_motif_disturbance, file=output_path("ssc1_tss_motif_disturbance.csv"), row.names=FALSE)
ssc2_tes_motif_disturbance <- analyze_motif_disturbance(ssc_1088_muts_TES, filename=output_path("ssc2_tes_motif_disturbance_per_variant.csv"))
write.csv(ssc2_tes_motif_disturbance, file=output_path("ssc2_tes_motif_disturbance.csv"), row.names=FALSE)
ssc2_tss_motif_disturbance <- analyze_motif_disturbance(ssc_1088_muts_TSS, filename=output_path("ssc2_tss_motif_disturbance_per_variant.csv"))
write.csv(ssc2_tss_motif_disturbance, file=output_path("ssc2_tss_motif_disturbance.csv"), row.names=FALSE)

chd_tes_motif_disturbance <- analyze_motif_disturbance(chd_combined_muts_TES, filename=output_path("chd_tes_motif_disturbance_per_variant.csv"))
write.csv(cdh_tes_motif_disturbance, file=output_path("chd_tes_motif_disturbance.csv"), row.names=FALSE)
chd_tss_motif_disturbance <- analyze_motif_disturbance(chd_combined_muts_TSS, filename=output_path("chd_tss_motif_disturbance_per_variant.csv"))
write.csv(cdh_tss_motif_disturbance, file=output_path("chd_tss_motif_disturbance.csv"), row.names=FALSE)


# Analyze motif disturbance in cases vs. controls
cdh_tes_motif_disturbance[which(abs(log2(cdh_tes_motif_disturbance$raw.score.y/cdh_tes_motif_disturbance$raw.score.x)) > 3),] # about 1/3 of variants that strongly disturb a TF binding site motif

pdf(output_path("delta_log2_motif_score_distribution.pdf"))
plot(density(log2(cdh_tes_motif_disturbance$raw.score.y/cdh_tes_motif_disturbance$raw.score.x)), main="delta_log2(score) distribution", col="red", lty=2, xlab="delta_log2(score)", ylim=c(0,0.7))
lines(density(log2(cdh_tss_motif_disturbance$raw.score.y/cdh_tss_motif_disturbance$raw.score.x)), col="orange", lty=2)
lines(density(log2(ssc1_tes_motif_disturbance$raw.score.y/ssc1_tes_motif_disturbance$raw.score.x)), col="blue", lty=2)
lines(density(log2(ssc1_tss_motif_disturbance$raw.score.y/ssc1_tss_motif_disturbance$raw.score.x)), col="green", lty=2)
legend("topleft", legend=c("CDH 3'UTR", "CDH TSS", "SSC 3'UTR", "SSC TSS"), col=c(c("red", "orange", "blue", "green")), lty=2)
dev.off()

#cdh_tss_motif_disturbance[cdh_tss_motif_disturbance$target == "GATA4",]
#nrow(unique(cdh_tss_motif_disturbance[cdh_tss_motif_disturbance$target == "GATA4",c("X.chr","pos")]))
write_motif_counts_tables <- function(tss_motif_disturbance, tes_motif_disturbance, name) {
    tss_motif_target_table <- table(tss_motif_disturbance$target); tss_motif_target_table <- sort(tss_motif_target_table[tss_motif_target_table != 0], decreasing=TRUE)
    write.csv(data.frame(tss_motif_target_table), file=output_path(paste0(name,"_tss_motif_disturbance_counts.csv")), row.names=FALSE)
    tes_motif_target_table <- table(tes_motif_disturbance$target); tes_motif_target_table <- sort(tes_motif_target_table[tes_motif_target_table != 0], decreasing=TRUE)
    write.csv(data.frame(tes_motif_target_table), file=output_path(paste0(name,"_tes_motif_disturbance_counts.csv")), row.names=FALSE)
    combined_motif_target_table <- merge(tss_motif_target_table, tes_motif_target_table, by="Var1", all.x=TRUE, all.y=TRUE)
    combined_motif_target_table <- cbind(combined_motif_target_table, combined_motif_target_table$Freq.x + combined_motif_target_table$Freq.y)
    colnames(combined_motif_target_table) <- c("disrupted_motif", "TSS_count", "TES_count", "total_count")
    combined_motif_target_table <- combined_motif_target_table[order(combined_motif_target_table$total_count, decreasing=TRUE),]
    write.csv(data.frame(combined_motif_target_table), file=output_path(paste0(name,"_motif_disturbance_counts.csv")), row.names=FALSE)
}
write_motif_counts_tables(cdh_tss_motif_disturbance, cdh_tes_motif_disturbance, "cdh")
write_motif_counts_tables(ssc1_tss_motif_disturbance, ssc1_tes_motif_disturbance, "ssc1")
write_motif_counts_tables(cdh_tss_motif_disturbance[which(abs(log2(cdh_tss_motif_disturbance$raw.score.y/cdh_tss_motif_disturbance$raw.score.x)) > 3),], cdh_tes_motif_disturbance[which(abs(log2(cdh_tes_motif_disturbance$raw.score.y/cdh_tes_motif_disturbance$raw.score.x)) > 3),], "cdh_strong")
write_motif_counts_tables(ssc1_tss_motif_disturbance[which(abs(log2(ssc1_tss_motif_disturbance$raw.score.y/ssc1_tss_motif_disturbance$raw.score.x)) > 3),], ssc1_tes_motif_disturbance[which(abs(log2(ssc1_tes_motif_disturbance$raw.score.y/ssc1_tes_motif_disturbance$raw.score.x)) > 3),], "ssc1_strong")

cdh_strong_motif_disturbance <- read.csv(output_path("cdh_strong_motif_disturbance_counts.csv"))
ssc1_strong_motif_disturbance <- read.csv(output_path("ssc1_strong_motif_disturbance_counts.csv"))
strong_motif_disturbance <- merge(cdh_strong_motif_disturbance, ssc1_strong_motif_disturbance, by="disrupted_motif", all.x=TRUE, all.y=TRUE)
strong_motif_disturbance[is.na(strong_motif_disturbance)] <- 0

pdf(output_path("cdh_vs_ssc_motif_disruption_comparison.pdf"))
strong_motif_disturbance[c(nrow(strong_motif_disturbance), which(strong_motif_disturbance$disrupted_motif == "GATA4")),] <- strong_motif_disturbance[c(which(strong_motif_disturbance$disrupted_motif == "GATA4"), nrow(strong_motif_disturbance)),]
colors <- rep("blue", nrow(strong_motif_disturbance)); colors[strong_motif_disturbance$disrupted_motif == "GATA4"] <- "red"
plot(strong_motif_disturbance$total_count.y, strong_motif_disturbance$total_count.x/1.116, main="Motif Disruption Comparison", xlab="count/sample in SSC", ylab="count/sample in CDH", col=colors, xlim=c(0,1000), ylim=c(0,1000))
lines(seq(-100,1100), seq(-100,1100))
#text(strong_motif_disturbance$total_count.y[strong_motif_disturbance$disrupted_motif == "E2F1"], strong_motif_disturbance$total_count.x[strong_motif_disturbance$disrupted_motif == "E2F1"]+35, "E2F1")
#text(strong_motif_disturbance$total_count.y[strong_motif_disturbance$disrupted_motif == "USF2"], strong_motif_disturbance$total_count.x[strong_motif_disturbance$disrupted_motif == "USF2"]+35, "USF2")
text(strong_motif_disturbance$total_count.y[strong_motif_disturbance$disrupted_motif == "GATA4"], strong_motif_disturbance$total_count.x[strong_motif_disturbance$disrupted_motif == "GATA4"]+80, "GATA4", col="red")
mtext("delta_log2(score) > 3")
dev.off()

pdf(output_path("cdh_vs_ssc_motif_disruption_comparison_candidate_genes.pdf"))
strong_motif_disturbance <- strong_motif_disturbance[strong_motif_disturbance$disrupted_motif %in% candidate_CDH_genes,]
strong_motif_disturbance[c(nrow(strong_motif_disturbance), which(strong_motif_disturbance$disrupted_motif == "GATA4")),] <- strong_motif_disturbance[c(which(strong_motif_disturbance$disrupted_motif == "GATA4"), nrow(strong_motif_disturbance)),]
colors <- rep("blue", nrow(strong_motif_disturbance)); colors[strong_motif_disturbance$disrupted_motif == "GATA4"] <- "red"
plot(strong_motif_disturbance$total_count.y[strong_motif_disturbance$disrupted_motif %in% candidate_CDH_genes], strong_motif_disturbance$total_count.x[strong_motif_disturbance$disrupted_motif %in% candidate_CDH_genes]/1.116, main="Motif Disruption Comparison", xlab="count/sample in SSC", ylab="count/sample in CDH", col=colors, xlim=c(0,400), ylim=c(0,400))
lines(seq(-100,1100), seq(-100,1100))
text(strong_motif_disturbance$total_count.y[strong_motif_disturbance$disrupted_motif == "GATA4"], strong_motif_disturbance$total_count.x[strong_motif_disturbance$disrupted_motif == "GATA4"]+40, "GATA4", col="red")
mtext("candidate mouse genes (n=61), delta_log2(score) > 3")
dev.off()


strong_motif_disturbance[which(strong_motif_disturbance$total_count.x > strong_motif_disturbance$total_count.y & strong_motif_disturbance$total_count.y > 300),]
sum(strong_motif_disturbance$total_count.x > strong_motif_disturbance$total_count.y)
nrow(strong_motif_disturbance)
strong_motif_disturbance[strong_motif_disturbance$disrupted_motif %in% c("GATA4", "GATA6"),]


# extract the 3 PWMs for the TFs we are interested in
ids = report$id[1:3]
sel.pwms = PWMLogn.dm3.MotifDb.Dmel$pwms[ids]
names(sel.pwms) = report$target[1:3]
# scan and get the raw scores
scores = motifScores(sequence, sel.pwms, raw.scores=TRUE)
# raw scores for the first (and only) input sequence
dim(scores[[1]])
## [1] 968 3
head(scores[[1]])
## bcd gt Kr
## [1,] 7.484914e-05 4.213929e-05 1.141957e-07
## [2,] 9.180413e-02 4.275114e-04 1.162378e-03
## [3,] 1.020698e-02 2.326263e+00 1.480311e-02
## [4,] 2.206202e-07 4.600757e-07 2.085725e-07
## [5,] 7.044890e-06 7.690586e-07 1.638103e-06
## [6,] 4.913950e-04 7.229475e-08 4.625971e-07
# score starting at position 1 of forward strand
scores[[1]][1, "bcd"]
## bcd
## 7.484914e-05
# score for the reverse complement of the motif, starting at the same position
scores[[1]][485, "bcd"]
## bcd
## 2.055192e-06
# plot
plotMotifScores(scores, cols=c("green", "red", "blue"))

hgmd <- read.csv(data_path("hgmd_annovar_input_2.txt"), sep="\t")
hgmd_granges <- to_genomic_regions(hgmd, chr_colname="Chrom", start_colname="Start", end_colname="End", label_colname="var_id")
hgmd_TSS40_granges <- hgmd_granges[unique(queryHits(findOverlaps(hgmd_granges, TSS40_granges)))]
hgmd_TSS_granges <- hgmd_granges[unique(queryHits(findOverlaps(hgmd_granges, TSS_granges)))]
hgmd_TES_granges <- hgmd_granges[unique(queryHits(findOverlaps(hgmd_granges, TES_granges)))]
print(paste0(length(hgmd_granges)," HGMD total regulatory and splicing variants:"))
print(paste0(length(hgmd_TSS40_granges)," in TSS40 regions (",sum(grepl("Regulatory", names(hgmd_TSS40_granges)))," regulatory, ",sum(grepl("Splicing", names(hgmd_TSS40_granges)))," splicing)"))
print(paste0(length(hgmd_TSS_granges)," in TSS regions (",sum(grepl("Regulatory", names(hgmd_TSS_granges)))," regulatory, ",sum(grepl("Splicing", names(hgmd_TSS_granges)))," splicing)"))
print(paste0(length(hgmd_TES_granges)," in 3'UTR regions (",sum(grepl("Regulatory", names(hgmd_TES_granges)))," regulatory, ",sum(grepl("Splicing", names(hgmd_TES_granges)))," splicing)"))

random_TES_variants <- generate_random_variants(sum(grepl("Regulatory", names(hgmd_TES_granges))), TES_granges)
random_TSS40_variants <- generate_random_variants(sum(grepl("Regulatory", names(hgmd_TSS40_granges))), TSS40_granges)

#random_gnomad_variants <- read.csv(data_path("random_gnomad_variants.txt"), sep="\t", header=FALSE)[,c(1,2,4,5)]
random_gnomad_variants <- read.csv(data_path("gnomAD/random_gnomad_variants_3m.txt"), sep="\t", header=FALSE)[,c(1,2,4,5)]
random_gnomad_variants <- cbind(random_gnomad_variants, rep("random",nrow(random_gnomad_variants))); colnames(random_gnomad_variants) <- c("Chrom", "Position", "Ref", "Alt", "sample")
random_gnomad_variants$Chrom <- paste0(random_gnomad_variants$Chrom); random_gnomad_variants$Position <- as.numeric(paste0(random_gnomad_variants$Position))
random_gnomad_variants <- random_gnomad_variants[!is.na(random_gnomad_variants$Chrom) & !is.na(random_gnomad_variants$Position),]
random_gnomad_granges <- to_genomic_regions(random_gnomad_variants, chr_colname="Chrom", start_colname="Position", label_colname="sample")
random_gnomad_regulatory_indices <- unique(c(queryHits(findOverlaps(random_gnomad_granges, TSS40_granges)), queryHits(findOverlaps(random_gnomad_granges, TES_granges))))
random_gnomad_granges <- random_gnomad_granges[random_gnomad_regulatory_indices]
random_gnomad_variants <- random_gnomad_variants[random_gnomad_regulatory_indices,]

regulatory_variant_dat <-  hgmd[grepl("Regulatory", hgmd$var_id),c(1,2,4,5,6)]; colnames(regulatory_variant_dat) <- c("Chrom", "Position", "Ref", "Alt", "sample"); regulatory_variant_dat$Chrom <- gsub("chr","",regulatory_variant_dat$Chrom)
regulatory_variant_dat <- rbind(regulatory_variant_dat, random_gnomad_variants)
#regulatory_variant_dat <- rbind(regulatory_variant_dat, random_TES_variants, random_TSS40_variants)
write.table(regulatory_variant_dat, file=output_path("regulatory_variants.tsv"), sep="\t", row.names=FALSE, quote=FALSE) # For annotation with WGSA

regulatory_variant_dat_snps <- read.csv(output_path("regulatory_variants_annotated.tsv.snp"), sep="\t"); regulatory_variant_dat_snps <- regulatory_variant_dat_snps[!grepl(",", regulatory_variant_dat_snps$Ref) & !grepl(",", regulatory_variant_dat_snps$Alt),]
regulatory_variant_dat_indels <- read.csv(output_path("regulatory_variants_annotated.tsv.indel"), sep="\t"); regulatory_variant_dat_indels <- regulatory_variant_dat_indels[!grepl(",", regulatory_variant_dat_indels$Ref) & !grepl(",", regulatory_variant_dat_indels$Alt),]
regulatory_variant_dat <- standardize_colnames(rbind(regulatory_variant_dat_snps[,1:5], regulatory_variant_dat_indels[,1:5]), remove_chr_prefix=TRUE)
snv_indel <- rep("snv", nrow(regulatory_variant_dat)); snv_indel[is_indel(regulatory_variant_dat$Ref, regulatory_variant_dat$Alt)] <- "indel"
regulatory_variant_dat <- cbind(regulatory_variant_dat, snv_indel)
regulatory_variant_granges <- to_genomic_regions(regulatory_variant_dat, chr_colname="Chrom", start_colname="Position", end_colname="Position", labels=regulatory_variant_dat$snv_indel)
regulatory_variant_dat_wgsa <- clean_wgsa_data(regulatory_variant_dat_snps, regulatory_variant_dat_indels, snp_granges=regulatory_variant_granges[names(regulatory_variant_granges)=="snv"], indel_granges=regulatory_variant_granges[names(regulatory_variant_granges)=="indel"])
regulatory_variant_dat_muts_wgsa <- regulatory_variant_dat_wgsa[["muts"]]; colnames(regulatory_variant_dat_muts_wgsa)[1:5] <- c("X.chr", "pos", "ref", "alt", "sample")
regulatory_muts_wgsa <- regulatory_variant_dat_muts_wgsa
annotate_with_RBP("regulatory")
regulatory_muts_wgsa$RBP0[!is.na(regulatory_muts_wgsa$RBP0)] <- 1; regulatory_muts_wgsa$RBP0[is.na(regulatory_muts_wgsa$RBP0)] <- 0; regulatory_muts_wgsa$RBP0 <- as.numeric(regulatory_muts_wgsa$RBP0)
regulatory_muts_wgsa$RBP[!is.na(regulatory_muts_wgsa$RBP)] <- 1; regulatory_muts_wgsa$RBP[is.na(regulatory_muts_wgsa$RBP)] <- 0; regulatory_muts_wgsa$RBP <- as.numeric(regulatory_muts_wgsa$RBP)
regulatory_variant_dat <- regulatory_muts_wgsa
regulatory_variant_dat <- cbind(rep(1, nrow(regulatory_variant_dat)), regulatory_variant_dat); colnames(regulatory_variant_dat)[1] <- "hgmd"; regulatory_variant_dat$hgmd[grepl("random",regulatory_variant_dat$sample)] <- 0

# Add sequence context features:
add_sequence_context_features <- function(dat, chr_colname="X.chr", pos_colname="pos", width=12, version="hg19") {
    sequences <- get_trimers(dat[,chr_colname], dat[,pos_colname], width=(2*width+1), allow_BSgenome=TRUE, version=version)$trimers
    sequence_context_feature_names <- paste0(width,"_bp")
    sapply(1:width, function(start_pos) substr(sequences,start_pos,start_pos+width-1))
    
    return(dat)
}
#haha <- add_sequence_context_features(regulatory_muts_wgsa)

regulatory_variant_dat <- regulatory_variant_dat[,-c(2,3,4,5,6)] #regulatory_variant_dat <- regulatory_variant_dat[,-c(2,3,4,5,6)]
regulatory_variant_dat <- regulatory_variant_dat[c(which(regulatory_variant_dat$hgmd == 1), sample(which(regulatory_variant_dat$hgmd == 0), sum(regulatory_variant_dat$hgmd), replace=FALSE)),]

regulatory_variant_dat$ANNOVAR_ensembl_Effect <- sapply(regulatory_variant_dat$ANNOVAR_ensembl_Effect, function(x) strsplit(paste0(x), "\\|")[[1]][1])
regulatory_variant_dat$RBP[!is.na(regulatory_variant_dat$RBP)] <- 1; regulatory_variant_dat$RBP[is.na(regulatory_variant_dat$RBP)] <- 0; regulatory_variant_dat$RBP <- as.numeric(regulatory_variant_dat$RBP)
regulatory_variant_dat$RBP0[!is.na(regulatory_variant_dat$RBP0)] <- 1; regulatory_variant_dat$RBP0[is.na(regulatory_variant_dat$RBP0)] <- 0; regulatory_variant_dat$RBP0 <- as.numeric(regulatory_variant_dat$RBP0)
# Remove unwanted WGSA columns
regulatory_variant_dat <- regulatory_variant_dat[,!(colnames(regulatory_variant_dat) %in% c("ANNOVAR_ensembl_Effect", "ANNOVAR_ensembl_Transcript_ID", "ANNOVAR_ensembl_Gene_ID", "ANNOVAR_ensembl_Closest_gene.intergenic_only.", "ANNOVAR_ensembl_HGVSc", "ANNOVAR_ensembl_HGVSp", "ANNOVAR_ensembl_Exon_Rank", 
                                                                                            "ANNOVAR_ensembl_summary", "ANNOVAR_refseq_Effect", "ANNOVAR_refseq_Transcript_ID", "ANNOVAR_refseq_Gene_ID", "ANNOVAR_refseq_Closest_gene.intergenic_only.", "ANNOVAR_refseq_HGVSc", 
                                                                                            "ANNOVAR_refseq_HGVSp", "ANNOVAR_refseq_Exon_Rank", "ANNOVAR_refseq_summary", "splicing_consensus_ada_score", "splicing_consensus_rf_score", "clinvar_rs", "clinvar_clnsig", "clinvar_trait", 
                                                                                            "clinvar_golden_stars", "COSMIC_ID", "ExAC_AC", "ExAC_AF", "ExAC_Adj_AC", "ExAC_Adj_AF", "ExAC_AFR_AC", "ExAC_AFR_AF", "ExAC_AMR_AC", "ExAC_AMR_AF", "ExAC_EAS_AC", "ExAC_EAS_AF", "ExAC_FIN_AC", 
                                                                                            "ExAC_FIN_AF", "ExAC_NFE_AC", "ExAC_NFE_AF", "ExAC_SAS_AC", "fathmm.MKL_non.coding_rankscore", "fathmm.MKL_non.coding_pred", "fathmm.MKL_non.coding_group", "fathmm.MKL_coding_rankscore", 
                                                                                            "fathmm.MKL_coding_pred", "fathmm.MKL_coding_group", "Eigen.phred", "Eigen_coding_or_noncoding", "SuperEnhancer_tissue_cell", "SuperEnhancer_RefSeq_id", "SuperEnhancer_Gene_symbol", "phastCons100way_vertebrate", "phastCons100way_vertebrate_rankscore"))]
# Make numeric columns numeric.
for(i in 1:ncol(regulatory_variant_dat)) {
    feature <- regulatory_variant_dat[,i]
    feature_numeric <- as.numeric(paste0(feature[feature != "."]))
    if(sum(!is.na(feature_numeric)) == 0) { next # not a numeric feature
    } else {
        regulatory_variant_dat[,i] <- rep(mean(feature_numeric), nrow(regulatory_variant_dat))
        regulatory_variant_dat[,i][feature != "."] <- feature_numeric
    }
}

#tests <- append_test("FANTOM5_enhancer_robust", "SNPs", cases_snps[which(cases_snps$FANTOM5_enhancer_robust == "Y"),], controls_snps[which(controls_snps$FANTOM5_enhancer_robust == "Y"),], "Enhancer elements found with CAGE (Cap Analysis of Gene Expression)")
#tests <- append_test("CADD > 10", "SNPs", cases_snps[which(as.numeric(paste0(cases_snps$CADD_phred)) > 10),], controls_snps[which(as.numeric(controls_snps$CADD_phred) > 10),], "")
#tests <- append_test("Eigen.phred > 10", "SNPs", cases_snps[which(as.numeric(paste0(cases_snps$Eigen.phred)) > 10),], controls_snps[which(as.numeric(paste0(controls_snps$Eigen.phred)) > 10),], "")
#tests <- append_test("GERP_RS > 3", "SNPs", cases_snps[which(as.numeric(paste0(cases_snps$GERP_RS)) > 3),], controls_snps[which(as.numeric(paste0(controls_snps$GERP_RS)) > 3),], "")
#tests <- append_test("funseq2_noncoding_rankscore > 0.9", "SNPs", cases_snps[which(as.numeric(paste0(cases_snps$funseq2_noncoding_rankscore)) > 0.9),], controls_snps[which(as.numeric(paste0(controls_snps$funseq2_noncoding_rankscore)) > 0.9),], "")



run_adaboost <- function(dat, mtext="", cv=10) {
    library("fastAdaboost")
    set.seed(9999)
    dat <- dat[sample(1:num_rows),]
    num_rows = nrow(dat); test_set_size = floor(num_rows/cv)
    labels <- c(); pred_scores <- c(); pred_scores_no_rbp <- c(); pred_scores_cadd <- c(); pred_scores_eigen <- c()
    allowed_columns <- which(!grepl("GERP|COSMIC|fathmm|Eigen|CADD|funseq2|FANTOM5", colnames(dat)))
    for(cv_i in 1:cv) {
        print(paste0("CV: ",cv_i," / ",cv))
        test_rows <- ((cv_i-1)*test_set_size+1):(cv_i*test_set_size) #(1:nrow(dat)) %in% sample(1:nrow(dat), floor(0.2*nrow(dat)), replace=FALSE)
        # RBP
        print("Training Adaboost with RBP...")
        train_dat <- dat[-c(test_rows),allowed_columns]; test_dat <- dat[c(test_rows),allowed_columns]
        a <- adaboost(hgmd ~ ., train_dat, 10)
        pred <- predict(a, newdata=test_dat)
        # No RBP
        print("Training Adaboost without RBP...")
        train_dat <- train_dat[,which(!grepl("RBP", colnames(train_dat)))]; test_dat <- test_dat[,which(!grepl("RBP", colnames(test_dat)))]
        a <- adaboost(hgmd ~ ., train_dat, 10)
        pred_no_rbp <- predict(a, newdata=test_dat)
        # CADD
        print("Training Adaboost with only CADD...")
        train_dat <- dat[-c(test_rows),c("hgmd","CADD_phred")]; test_dat <- dat[c(test_rows),c("hgmd","CADD_phred")]
        a <- adaboost(hgmd ~ ., train_dat, 10)
        pred_cadd <- predict(a, newdata=test_dat)
        # Eigen
        print("Training Adaboost with only Eigen...")
        train_dat <- dat[-c(test_rows),c("hgmd","Eigen.raw")]; test_dat <- dat[c(test_rows),c("hgmd","Eigen.raw")]
        a <- adaboost(hgmd ~ ., train_dat, 10)
        pred_eigen <- predict(a, newdata=test_dat)
        
        labels <- c(labels, test_dat$hgmd == 1)
        pred_scores <- c(pred_scores, pred$prob[,1])
        pred_scores_no_rbp <- c(pred_scores_no_rbp, pred_no_rbp$prob[,1])
        pred_scores_cadd <- c(pred_scores_cadd, pred_cadd$prob[,1])
        pred_scores_eigen <- c(pred_scores_eigen, pred_eigen$prob[,1])
        #pred <- predict(a, newdata=regulatory_variant_dat[!test_rows,])
        #print(pred$error)
        #print(table(pred$class, regulatory_variant_dat[!test_rows,]$hgmd))
        #print(pred$error)
        #print(table(pred$class, regulatory_variant_dat[test_rows,]$hgmd))
        #pdf(file=output_path("adaboost_trees.pdf"))
        #par(mfrow=c(2,2))
        #plot(get_tree(a, 1)[[1]], main="AdaBoost Tree 1"); mtext(paste0("weight: ",a$weights[1]))
        #plot(get_tree(a, 2)[[1]], main="AdaBoost Tree 2"); mtext(paste0("weight: ",a$weights[2]))
        #plot(get_tree(a, 3)[[1]], main="AdaBoost Tree 3"); mtext(paste0("weight: ",a$weights[3]))
        #plot(get_tree(a, 4)[[1]], main="AdaBoost Tree 4"); mtext(paste0("weight: ",a$weights[4]))
        #par(mfrow=c(1,1))
        #dev.off()
    }
    
    # Draw ROC curve
    # RBP
    roc_result <- sapply(seq(0,1.01,by=0.01), function(cutoff) { calls <- pred_scores < cutoff; sensitivity <- sum(calls[labels==1]==labels[labels==1])/sum(labels==1); one_minus_specificity <- 1 - (sum(calls[labels==0]==labels[labels==0])/sum(labels==0)); return(c(sensitivity, one_minus_specificity)) })
    sensitivity <- roc_result[1,]
    one_minus_specificity <- roc_result[2,]
    pdf(file=output_path("roc_curve.pdf"))
    plot(c(0,1), c(0,1), type="l", col="grey", xaxs="i", yaxs="i", xlab="1 - Specificity", ylab="Sensitivity", main="ROC Curve", cex.axis=1.3, cex.lab=1.3, cex.main=1.2)
    lines(one_minus_specificity, sensitivity, col="red")
    auc_all = sum(diff(one_minus_specificity) * rollmean(sensitivity, 2))
    # No RBP
    roc_result <- sapply(seq(0,1.01,by=0.01), function(cutoff) { calls <- pred_scores_no_rbp < cutoff; sensitivity <- sum(calls[labels==1]==labels[labels==1])/sum(labels==1); one_minus_specificity <- 1 - (sum(calls[labels==0]==labels[labels==0])/sum(labels==0)); return(c(sensitivity, one_minus_specificity)) })
    sensitivity <- roc_result[1,]
    one_minus_specificity <- roc_result[2,]
    lines(one_minus_specificity, sensitivity, col="green")
    auc_no_rbp = sum(diff(one_minus_specificity) * rollmean(sensitivity, 2))
    # CADD
    roc_result <- sapply(seq(0,1.01,by=0.01), function(cutoff) { calls <- pred_scores_cadd < cutoff; sensitivity <- sum(calls[labels==1]==labels[labels==1])/sum(labels==1); one_minus_specificity <- 1 - (sum(calls[labels==0]==labels[labels==0])/sum(labels==0)); return(c(sensitivity, one_minus_specificity)) })
    sensitivity <- roc_result[1,]
    one_minus_specificity <- roc_result[2,]
    lines(one_minus_specificity, sensitivity, col="blue")
    auc_cadd = sum(diff(one_minus_specificity) * rollmean(sensitivity, 2))
    # No RBP
    roc_result <- sapply(seq(0,1.01,by=0.01), function(cutoff) { calls <- pred_scores_eigen < cutoff; sensitivity <- sum(calls[labels==1]==labels[labels==1])/sum(labels==1); one_minus_specificity <- 1 - (sum(calls[labels==0]==labels[labels==0])/sum(labels==0)); return(c(sensitivity, one_minus_specificity)) })
    sensitivity <- roc_result[1,]
    one_minus_specificity <- roc_result[2,]
    lines(one_minus_specificity, sensitivity, col="cyan")
    auc_eigen = sum(diff(one_minus_specificity) * rollmean(sensitivity, 2))
    
    mtext(paste0(cv,"-fold CV, AUC_RBP = ", round(auc_all,3), ", AUC_no_RBP = ", round(auc_no_rbp,3)), cex=1.1)
    legend("bottomright", legend=c("RBP", "No RBP", "CADD", "Eigen"), col=c("red", "green", "blue", "cyan"), pch=15)
    dev.off()
}
run_adaboost(regulatory_variant_dat, cv=10)

# Find adaboost score distribution on real WGS data, and compare case distribution to control distribution.
get_adaboost_scores <- function(train_dat, new_dat) {
    library("fastAdaboost")
    # Train Adaboost classifier on train_dat (regulatory_variant_dat, in our case).
    allowed_columns <- which(!grepl("GERP|COSMIC|fathmm|Eigen|CADD|funseq2|FANTOM5", colnames(train_dat)))
    train_dat <- train_dat[,allowed_columns]
    print("Training AdaBoost classifier on training data...")
    a <- adaboost(hgmd ~ ., train_dat, 10)
    
    # Prepare new dataset, for which to find the prediction scores for.
    new_dat <- new_dat[,colnames(new_dat) %in% colnames(train_dat)]
    #new_dat$ANNOVAR_ensembl_Effect <- sapply(new_dat$ANNOVAR_ensembl_Effect, function(x) strsplit(paste0(x), "\\|")[[1]][1])
    new_dat$RBP[!is.na(new_dat$RBP)] <- 1; new_dat$RBP[is.na(new_dat$RBP)] <- 0; new_dat$RBP <- as.numeric(new_dat$RBP)
    new_dat$RBP0[!is.na(new_dat$RBP0)] <- 1; new_dat$RBP0[is.na(new_dat$RBP0)] <- 0; new_dat$RBP0 <- as.numeric(new_dat$RBP0)
    
    # Make new_dat's numeric columns numeric.
    print("Making numeric columns numeric...")
    for(i in 1:ncol(new_dat)) {
        feature <- new_dat[,i]
        feature_numeric <- as.numeric(paste0(feature[feature != "."]))
        if(sum(!is.na(feature_numeric)) == 0) { next # not a numeric feature
        } else {
            new_dat[,i] <- rep(mean(feature_numeric), nrow(new_dat))
            new_dat[,i][feature != "."] <- feature_numeric
        }
    }
    # Remove "NY" histone mark levels.
    new_dat <- new_dat[apply(new_dat,1, function(x) !("NY" %in% paste0(x))),]
    
    # Predict scores for new_dat
    print("Predicting scores...")
    new_dat_pred <- predict(a, newdata=new_dat)$prob[,1]
    return(new_dat_pred)
}
chdfb_combined_pred <- get_adaboost_scores(regulatory_variant_dat, chdfb_combined_muts_wgsa[which(!is.na(chdfb_combined_muts_wgsa$RBP)),])
sscfb_combined_pred <- get_adaboost_scores(regulatory_variant_dat, sscfb_combined_muts_wgsa[which(!is.na(sscfb_combined_muts_wgsa$RBP)),])
chdfb_combined_pred <- get_adaboost_scores(regulatory_variant_dat, chdfb_combined_muts_wgsa[which(apply(chdfb_combined_muts_wgsa[,grepl("H3K36me3",colnames(chdfb_combined_muts_wgsa))], 1, function(x) sum(grepl("Y",x))>0)),])
sscfb_combined_pred <- get_adaboost_scores(regulatory_variant_dat, sscfb_combined_muts_wgsa[which(apply(sscfb_combined_muts_wgsa[,grepl("H3K36me3",colnames(sscfb_combined_muts_wgsa))], 1, function(x) sum(grepl("Y",x))>0)),])# Draw plot comparing score densities
pdf(file=output_path("wgs_adaboost_score_distribution_fb.pdf"))
plot(density(chdfb_combined_pred, bw=0.1), col="red", cex.axis=1.3, cex.lab=1.3, cex.main=1.2)
lines(density(sscfb_combined_pred, bw=0.1), col="green")
legend("topleft", legend=c("CHDFB cases", "SSCFB controls"), col=c("red", "green"), pch=15)
dev.off()
# cdh_combined_pred <- get_adaboost_scores(regulatory_variant_dat, cdh_combined_muts[which(!is.na(cdh_combined_muts$RBP)),])
# chd_combined_pred <- get_adaboost_scores(regulatory_variant_dat, chd_combined_muts[which(!is.na(cdh_combined_muts$RBP)),])
# ssc_all_pred <- get_adaboost_scores(regulatory_variant_dat, ssc_all_muts_TSS[which(!is.na(cdh_combined_muts$RBP)),])
# pdf(file=output_path("wgs_adaboost_score_distribution.pdf"))
# plot(density(cdh_combined_pred, bw=0.1), col="blue", cex.axis=1.3, cex.lab=1.3, cex.main=1.2)
# lines(density(chd_combined_pred, bw=0.1), col="red")
# lines(density(ssc_all_pred, bw=0.1), col="green")
# legend("topleft", legend=c("CDH cases", "CHD cases", "SSC controls"), col=c("blue", "red", "green"), pch=15)
# dev.off()


# Send Felix and Sarah all gene annotations
# CHDFB
chdfb_gene_annotations <- merge(chdfb_combined_muts_wgsa[,c("X.chr", "pos", "ref", "alt", "sample", "Chrom_hg38", "Position_hg38")], chdfb_combined_muts_TES[,c("X.chr", "pos", "ref", "alt", "sample", "TS_gene")], all.x=TRUE)
colnames(chdfb_gene_annotations)[ncol(chdfb_gene_annotations)] <- "3'UTR_genes"
chdfb_gene_annotations <- merge(chdfb_gene_annotations, chdfb_combined_muts_TSS[,c("X.chr", "pos", "ref", "alt", "sample", "TS_gene")], all.x=TRUE)
colnames(chdfb_gene_annotations)[ncol(chdfb_gene_annotations)] <- "TSS_genes"
chdfb_gene_annotations <- merge(chdfb_gene_annotations, chdfb_combined_muts_TSS40[,c("X.chr", "pos", "ref", "alt", "sample", "TS_gene")], all.x=TRUE)
colnames(chdfb_gene_annotations)[ncol(chdfb_gene_annotations)] <- "TSS40_genes"
colnames(chdfb_gene_annotations)[1:4] <- c("Chrom_hg19", "Position_hg19", "Ref", "Alt")
chdfb_gene_annotations <- chdfb_gene_annotations[!(is.na(chdfb_gene_annotations[,"3'UTR_genes"]) & is.na(chdfb_gene_annotations[,"TSS_genes"]) & is.na(chdfb_gene_annotations[,"TSS40_genes"])),]
nrow(chdfb_combined_muts_wgsa)
nrow(chdfb_gene_annotations)
chdfb_gene_annotations$Chrom_hg19 <- paste0("chr",chdfb_gene_annotations$Chrom_hg19)
chdfb_gene_annotations <- chdfb_gene_annotations[,c("Chrom_hg38", "Position_hg38", "Ref", "Alt", "sample", "3'UTR_genes", "TSS_genes", "TSS40_genes", "Chrom_hg19", "Position_hg19")]
chdfb_gene_annotations[1:5,]
write.table(chdfb_gene_annotations, file=output_path("chdfb_gene_annotations.tsv"), sep="\t", row.names=FALSE)
# SSCFB
sscfb_gene_annotations <- merge(sscfb_combined_muts_wgsa[,c("X.chr", "pos", "ref", "alt", "sample", "Chrom_hg38", "Position_hg38")], sscfb_combined_muts_TES[,c("X.chr", "pos", "ref", "alt", "sample", "TS_gene")], all.x=TRUE)
colnames(sscfb_gene_annotations)[ncol(sscfb_gene_annotations)] <- "3'UTR_genes"
sscfb_gene_annotations <- merge(sscfb_gene_annotations, sscfb_combined_muts_TSS[,c("X.chr", "pos", "ref", "alt", "sample", "TS_gene")], all.x=TRUE)
colnames(sscfb_gene_annotations)[ncol(sscfb_gene_annotations)] <- "TSS_genes"
sscfb_gene_annotations <- merge(sscfb_gene_annotations, sscfb_combined_muts_TSS40[,c("X.chr", "pos", "ref", "alt", "sample", "TS_gene")], all.x=TRUE)
colnames(sscfb_gene_annotations)[ncol(sscfb_gene_annotations)] <- "TSS40_genes"
colnames(sscfb_gene_annotations)[1:4] <- c("Chrom_hg19", "Position_hg19", "Ref", "Alt")
sscfb_gene_annotations <- sscfb_gene_annotations[!(is.na(sscfb_gene_annotations[,"3'UTR_genes"]) & is.na(sscfb_gene_annotations[,"TSS_genes"]) & is.na(sscfb_gene_annotations[,"TSS40_genes"])),]
nrow(sscfb_combined_muts_wgsa)
nrow(sscfb_gene_annotations)
sscfb_gene_annotations$Chrom_hg19 <- paste0("chr",sscfb_gene_annotations$Chrom_hg19)
sscfb_gene_annotations <- sscfb_gene_annotations[,c("Chrom_hg38", "Position_hg38", "Ref", "Alt", "sample", "3'UTR_genes", "TSS_genes", "TSS40_genes", "Chrom_hg19", "Position_hg19")]
sscfb_gene_annotations[1:5,]
write.table(sscfb_gene_annotations, file=output_path("sscfb_gene_annotations.tsv"), sep="\t", row.names=FALSE)

# Phastcons cutoff analysis
a_info <- strsplit(readLines(data_path("WGSA/phastCons/chr16.phastCons46way.placental.wigFix.gz.rank0"), n=1), " ")[[1]]; a_start = strsplit(a_info[3], "=")[[1]][2]
a <- read.csv(data_path("WGSA/phastCons/chr16.phastCons46way.placental.wigFix.gz.rank0"), sep="\t", header=FALSE, skip=1)
a$V1 <- as.numeric(paste0(a$V1)); a$V2 <- as.numeric(paste0(a$V2)); sum(is.na(a$V1)); sum(is.na(a$V2))
a <- a[!(is.na(a$V1) | is.na(a$V2)),]
density(a$V1)
density(a$V2)
sum(a$V1 > 0.5) / nrow(a)
sum(a$V2 > 0.5) / nrow(a)

# CDH phenotype analysis (including isolated vs. complex)
cdh_info <- read.csv(data_path("CDH_baylor_extrainfo.csv"))

a <- read.csv(data_path("felix_chdfb_hg19_vars_to_rm.csv"))
nrow(a)
sum(a$segmentalDuplications)
sum(a$LCR.hs37d5_chr)
sum(!a$mappability1_300)
sum(a$genes.MUC.HLA)
sum(a$encode_duke_blacklist)
sum(a$exonic_by_refseq_or_ensembl)
sum(a$dn_within_500bp)

a <- read.csv(data_path("felix_dnvs_filtered_hg38.csv"))


chdfb_combined_granges_expanded <- chdfb_combined_granges; start(chdfb_combined_granges_expanded) <- start(chdfb_combined_granges_expanded) - 250; end(chdfb_combined_granges_expanded) <- end(chdfb_combined_granges_expanded) + 250
chdfb_combined_granges_clusters <- findOverlaps(chdfb_combined_granges_expanded, chdfb_combined_granges_expanded); chdfb_combined_granges_clusters <- chdfb_combined_granges_clusters[queryHits(chdfb_combined_granges_clusters) != subjectHits(chdfb_combined_granges_clusters)]
chdfb_denovo_clusters <- cluster_overlapping_genes(paste0(paste0(chdfb_combined_muts_wgsa$sample,"_",chdfb_combined_muts_wgsa$Chrom_hg38,":",chdfb_combined_muts_wgsa$Position_hg38,chdfb_combined_muts_wgsa$ref,">",chdfb_combined_muts_wgsa$alt)[queryHits(chdfb_combined_granges_clusters)] , ",", paste0(chdfb_combined_muts_wgsa$sample,"_",chdfb_combined_muts_wgsa$Chrom_hg38,":",chdfb_combined_muts_wgsa$Position_hg38,chdfb_combined_muts_wgsa$ref,">",chdfb_combined_muts_wgsa$alt)[subjectHits(chdfb_combined_granges_clusters)] ))
chdfb_denovo_clusters <- unique(chdfb_denovo_clusters)
chdfb_denovo_clusters <- unlist(lapply(strsplit(chdfb_denovo_clusters, ","), function(x) { 
    samples <- unlist(strsplit(x, "_")); samples <- samples[(1:length(samples)) %% 2 == 1]
    if(sum(duplicated(samples))>0) { return(paste0(x[samples %in% samples[duplicated(samples)]], collapse=",")) 
    } else { return(".") }
})); chdfb_denovo_clusters <- chdfb_denovo_clusters[chdfb_denovo_clusters != "."]
write(chdfb_denovo_clusters, file=output_path("chdfb_denovo_clusters_hg38.txt"))

chd_combined_granges_expanded <- chd_combined_granges; start(chd_combined_granges_expanded) <- start(chd_combined_granges_expanded) - 250; end(chd_combined_granges_expanded) <- end(chd_combined_granges_expanded) + 250
chd_combined_granges_clusters <- findOverlaps(chd_combined_granges_expanded, chd_combined_granges_expanded); chd_combined_granges_clusters <- chd_combined_granges_clusters[queryHits(chd_combined_granges_clusters) != subjectHits(chd_combined_granges_clusters)]
chd_denovo_clusters <- cluster_overlapping_genes(paste0(paste0(chd_combined$sample,"_",chd_combined$Chrom,":",chd_combined$Position,chd_combined$Ref,">",chd_combined$Alt)[queryHits(chd_combined_granges_clusters)] , ",", paste0(chd_combined$sample,"_",chd_combined$Chrom,":",chd_combined$Position,chd_combined$Ref,">",chd_combined$Alt)[subjectHits(chd_combined_granges_clusters)] ))
chd_denovo_clusters <- unique(chd_denovo_clusters)
chd_denovo_clusters <-  unlist(lapply(strsplit(chd_denovo_clusters, ","), function(x) { 
    samples <- unlist(strsplit(x, "_")); samples <- samples[(1:length(samples)) %% 2 == 1]
    if(sum(duplicated(samples))>0) { return(paste0(x[samples %in% samples[duplicated(samples)]], collapse=",")) 
    } else { return(".") }
})); chd_denovo_clusters <- chd_denovo_clusters[chd_denovo_clusters != "."]
write(chd_denovo_clusters, file=output_path("chd_denovo_clusters_hg19.txt"))

sscfb_combined_granges_expanded <- sscfb_combined_granges; start(sscfb_combined_granges_expanded) <- start(sscfb_combined_granges_expanded) - 250; end(sscfb_combined_granges_expanded) <- end(sscfb_combined_granges_expanded) + 250
sscfb_combined_granges_clusters <- findOverlaps(sscfb_combined_granges_expanded, sscfb_combined_granges_expanded); sscfb_combined_granges_clusters <- sscfb_combined_granges_clusters[queryHits(sscfb_combined_granges_clusters) != subjectHits(sscfb_combined_granges_clusters)]
sscfb_denovo_clusters <- cluster_overlapping_genes(paste0(paste0(sscfb_combined_muts_wgsa$sample,"_",sscfb_combined_muts_wgsa$Chrom_hg38,":",sscfb_combined_muts_wgsa$Position_hg38,sscfb_combined_muts_wgsa$ref,">",sscfb_combined_muts_wgsa$alt)[queryHits(sscfb_combined_granges_clusters)] , ",", paste0(sscfb_combined_muts_wgsa$sample,"_",sscfb_combined_muts_wgsa$Chrom_hg38,":",sscfb_combined_muts_wgsa$Position_hg38,sscfb_combined_muts_wgsa$ref,">",sscfb_combined_muts_wgsa$alt)[subjectHits(sscfb_combined_granges_clusters)] ))
sscfb_denovo_clusters <- unique(sscfb_denovo_clusters)
sscfb_denovo_clusters <-  unlist(lapply(strsplit(sscfb_denovo_clusters, ","), function(x) { 
    samples <- unlist(strsplit(x, "_")); samples <- samples[(1:length(samples)) %% 2 == 1]
    if(sum(duplicated(samples))>0) { return(paste0(x[samples %in% samples[duplicated(samples)]], collapse=",")) 
    } else { return(".") }
})); sscfb_denovo_clusters <- sscfb_denovo_clusters[sscfb_denovo_clusters != "."]
write(sscfb_denovo_clusters, file=output_path("sscfb_denovo_clusters_hg38.txt"))

ssc_all_granges_expanded <- ssc_all_granges; start(ssc_all_granges_expanded) <- start(ssc_all_granges_expanded) - 250; end(ssc_all_granges_expanded) <- end(ssc_all_granges_expanded) + 250
ssc_all_granges_clusters <- findOverlaps(ssc_all_granges_expanded, ssc_all_granges_expanded); ssc_all_granges_clusters <- ssc_all_granges_clusters[queryHits(ssc_all_granges_clusters) != subjectHits(ssc_all_granges_clusters)]
ssc_denovo_clusters <- cluster_overlapping_genes(paste0(paste0(ssc_all$sample,"_",ssc_all$Chrom,":",ssc_all$Position,ssc_all$Ref,">",ssc_all$Alt)[queryHits(ssc_all_granges_clusters)] , ",", paste0(ssc_all$sample,"_",ssc_all$Chrom,":",ssc_all$Position,ssc_all$Ref,">",ssc_all$Alt)[subjectHits(ssc_all_granges_clusters)] ))
ssc_denovo_clusters <- unique(ssc_denovo_clusters)
ssc_denovo_clusters <-  unlist(lapply(strsplit(ssc_denovo_clusters, ","), function(x) { 
    samples <- unlist(strsplit(x, "_")); samples <- samples[(1:length(samples)) %% 2 == 1]
    if(sum(duplicated(samples))>0) { return(paste0(x[samples %in% samples[duplicated(samples)]], collapse=",")) 
    } else { return(".") }
})); ssc_denovo_clusters <- ssc_denovo_clusters[ssc_denovo_clusters != "."]
write(ssc_denovo_clusters, file=output_path("ssc_denovo_clusters_hg19.txt"))

a <- function() {
    cat(paste0("CHD: # clusters (500bp): ", length(chd_denovo_clusters), ", # clustered variants/sample in CHD: ", length(unlist(strsplit(chd_denovo_clusters,",")))/CHD_COMBINED_SAMPLE_COUNT),"\n")
    cat(paste0("SSC: # clusters (500bp): ", length(ssc_denovo_clusters), ", # clustered variants/sample in SSC: ", length(unlist(strsplit(ssc_denovo_clusters,",")))/SSC_ALL_SAMPLE_COUNT),"\n")
    cat(paste0("CHDFB: # clusters (500bp): ", length(chdfb_denovo_clusters), ", # clustered variants/sample in CHDFB: ", length(unlist(strsplit(chdfb_denovo_clusters,",")))/CHDFB_COMBINED_SAMPLE_COUNT),"\n")
    cat(paste0("SSCFB: # clusters (500bp): ", length(sscfb_denovo_clusters), ", # clustered variants/sample in SSCFB: ", length(unlist(strsplit(sscfb_denovo_clusters,",")))/SSCFB_COMBINED_SAMPLE_COUNT),"\n")
}
a()

#################################################################################################################
# Load and process rare inherited variants.
#################################################################################################################
if(FALSE) {
    # Load rare inherited CDH noncoding data.
    if(file.exists(data_path("WGS\\CDH\\cdh_rare_noncoding.txt"))) {
        cdh_rare_nc <- read.csv(data_path("WGS\\CDH\\cdh_rare_noncoding.txt"), sep="\t")
    } else {
        cdh_rare_nc_unprocessed <- read.csv(data_path("WGS\\CDH\\cdh_rare_noncoding_unprocessed.txt"), sep="\t")
        cdh_rare_nc_unprocessed_genotype_split <- strsplit(paste0(cdh_rare_nc_unprocessed$genotype), ":|/")
        cdh_rare_nc_unprocessed_alt_split <- strsplit(paste0(cdh_rare_nc_unprocessed$Alt), ",")
        cdh_rare_nc_alt <- sapply(1:length(cdh_rare_nc_unprocessed_genotype_split), function(i) { gtype = as.numeric(cdh_rare_nc_unprocessed_genotype_split[[i]][2]); return(cdh_rare_nc_unprocessed_alt_split[[i]][gtype]) } )
        #cbind(cdh_rare_nc_unprocessed, cdh_rare_nc_alt) 
        cdh_rare_nc <- cdh_rare_nc_unprocessed[,-which(colnames(cdh_rare_nc_unprocessed) == "genotype")]
        cdh_rare_nc$Alt <- cdh_rare_nc_alt
        cdh_rare_nc <- cdh_rare_nc[paste0(cdh_rare_nc$Alt) != "*",]
        write.table(cdh_rare_nc, file=data_path("WGS\\CDH\\cdh_rare_noncoding.txt"), sep="\t", row.names=FALSE)
    }
    cdh_rare_nc_is_indel <- is_indel(cdh_rare_nc$Ref, cdh_rare_nc$Alt)
    cdh_rare_nc <- cbind(cdh_rare_nc, cdh_rare_nc_is_indel); colnames(cdh_rare_nc)[ncol(cdh_rare_nc)] <- "snv_indel"
    cdh_rare_nc$snv_indel[cdh_rare_nc$snv_indel] <- "indel"; cdh_rare_nc$snv_indel[cdh_rare_nc$snv_indel != "indel"] <- "snv"
    sum(cdh_rare_nc$snv_indel == "snv")
    sum(cdh_rare_nc$snv_indel == "indel")   
    
    # Load rare inherited SSC noncoding data.
    if(file.exists(data_path("WGS\\SSC\\ssc_rare_noncoding.txt"))) {
        ssc_rare_nc <- read.csv(data_path("WGS\\SSC\\ssc_rare_noncoding.txt"), sep="\t")
    } else {
        ssc_rare_nc_unprocessed <- read.csv(data_path("WGS\\SSC\\ssc_rare_noncoding_unprocessed.txt"), sep="\t")
        ssc_rare_nc_unprocessed_genotype_split <- strsplit(paste0(ssc_rare_nc_unprocessed$genotype), ":|/")
        ssc_rare_nc_unprocessed_alt_split <- strsplit(paste0(ssc_rare_nc_unprocessed$Alt), ",")
        ssc_rare_nc_alt <- sapply(1:length(ssc_rare_nc_unprocessed_genotype_split), function(i) { gtype = as.numeric(ssc_rare_nc_unprocessed_genotype_split[[i]][2]); return(ssc_rare_nc_unprocessed_alt_split[[i]][gtype]) } )
        #cbind(ssc_rare_nc_unprocessed, ssc_rare_nc_alt) 
        ssc_rare_nc <- ssc_rare_nc_unprocessed[,-which(colnames(ssc_rare_nc_unprocessed) == "genotype")]
        ssc_rare_nc$Alt <- ssc_rare_nc_alt
        ssc_rare_nc <- ssc_rare_nc[paste0(ssc_rare_nc$Alt) != "*",]
        write.table(ssc_rare_nc, file=data_path("WGS\\SSC\\ssc_rare_noncoding.txt"), sep="\t", row.names=FALSE)
    }
    ssc_rare_nc_is_indel <- is_indel(ssc_rare_nc$Ref, ssc_rare_nc$Alt)
    ssc_rare_nc <- cbind(ssc_rare_nc, ssc_rare_nc_is_indel); colnames(ssc_rare_nc)[ncol(ssc_rare_nc)] <- "snv_indel"
    ssc_rare_nc$snv_indel[ssc_rare_nc$snv_indel] <- "indel"; ssc_rare_nc$snv_indel[ssc_rare_nc$snv_indel != "indel"] <- "snv"
    sum(ssc_rare_nc$snv_indel == "snv")
    sum(ssc_rare_nc$snv_indel == "indel")
}

#################################################################################################################
# Download desired histone modification peaks from Roadmap Epigenomics
#################################################################################################################
download_roadmap_histone_peaks <- function(marks, eids, peak_types) {
    if(peak_types[1] == "all") { peak_types <- c("narrowPeak", "gappedPeak", "broadPeak") } # gappedPeak is the same as broadPeak, but selecting only those that contain at least one significant narrowPeak inside them, increasing confidence.
    if(eids[1] %in% c("CDH","CHD")) { eids <- get_relevant_roadmap_eids(eids) }
    if(marks[1] == "all") { marks = "*" }
    
    HISTONE_MODIFICATIONS_ANNOTATIONS_PATH = data_path("Roadmap")
    already_owned_histone_modification_files <- list.files(path=HISTONE_MODIFICATIONS_ANNOTATIONS_PATH)
    
    desired_histone_marks <- apply(expand.grid(eids, "-", marks, ".", peak_types), 1, paste0, collapse="")
    desired_histone_marks <- desired_histone_marks[!(desired_histone_marks %in% already_owned_histone_modification_files)]
    
    sapply_out <- sapply(peak_types, function(peak_type) {
        roadmap_url_base <- paste0("https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/",peak_type,"/")
        for(histone_mark in desired_histone_marks[grepl(peak_type, desired_histone_marks)]) {
            #roadmap_url <- paste0(roadmap_url_base,"/",histone_mark,".gz")
            #print(paste0("Fetching ",roadmap_url,"..."))
            if(marks[1] == "*") { wget_command = paste0("wget -nc -P ",HISTONE_MODIFICATIONS_ANNOTATIONS_PATH,"/ -r -l1 -np -nd -A '",histone_mark,".gz' ",roadmap_url_base)
            } else { wget_command = paste0("wget -nc -P ",HISTONE_MODIFICATIONS_ANNOTATIONS_PATH,"/ ",roadmap_url_base,histone_mark,".gz") }
            print(wget_command)
            system(wget_command)
        }
        return()
    })
    print(paste0("Unzipping files..."))
    system(paste0("gunzip -f ",HISTONE_MODIFICATIONS_ANNOTATIONS_PATH,"/*.gz"))
    print(paste0("Done."))
}
#get_relevant_roadmap_eids("ASD")[!(get_relevant_roadmap_eids("ASD") %in% get_relevant_roadmap_eids("CHD"))]
#download_roadmap_histone_peaks("H3K36me3", get_relevant_roadmap_eids("ASD")[!(get_relevant_roadmap_eids("ASD") %in% get_relevant_roadmap_eids("CHD"))], "all")
download_roadmap_histone_peaks("all", c("E080"), "broadPeak")

#################################################################################################################
# Unpack RBP eCLIP files downloaded from ENCODE
#################################################################################################################
unpack_RBP_eCLIP_peaks <- function(eclip_path, version="hg19", filter_quality=TRUE, p_filter=0.01, q_filter=0.2) {
    eclip_file_info <- read.csv(full_path(eclip_path, "metadata.tsv"), sep="\t")
    eclip_file_info <- eclip_file_info[sapply(strsplit(paste0(eclip_file_info$Biological.replicate.s.),","), length) > 1 & paste0(eclip_file_info$Assembly) == version,]
    
    heart_expression_ranks <- read.csv(file=output_path("heart_rank.csv"))
    # Unpack RBP eCLIP .gz files.
    print(paste0("Unpacking RBP eCLIP files..."))
    eclip_files = list.files(path=eclip_path, pattern="\\.bed.gz$")
    num_eclip_files = length(eclip_files)
    for(i in 1:num_eclip_files) {
        rbp_eclip_file = eclip_files[i]
        print(paste0(rbp_eclip_file,"...[",i," / ",num_eclip_files,"]"))
        eclip_accession = strsplit(rbp_eclip_file, "\\.")[[1]][1]
        rbp_info <- eclip_file_info[eclip_file_info$File.accession == eclip_accession, c("Biosample.term.name", "Experiment.target")]
        if(nrow(rbp_info) == 0) { next }
        rbp_info <- gsub(" ", "_", paste0(unlist(rbp_info)))
        rbp_target = strsplit(rbp_info[2], "-")[[1]][1]
        heart_rank = heart_expression_ranks[heart_expression_ranks$gene == rbp_target,"heart_rank"]
        if(length(heart_rank) < 1 || heart_rank >= 0.5) { next }
        unpacked_filepath = full_path(eclip_path, paste0(rbp_info[1],".",rbp_target,".bed"))
        packed_filepath = full_path(eclip_path, rbp_eclip_file)
        system(paste0("zcat ",packed_filepath," > ",unpacked_filepath)) # can use gunzip -c as well instead of zcat
        
        if(filter_quality) {
            rbp_dat <- read.csv(unpacked_filepath, sep="\t", header=FALSE)
            rbp_dat <- rbp_dat[exp(-rbp_dat$V8) < p_filter,] # p-value threshold
            if(nrow(rbp_dat)>0) { write.table(rbp_dat, file=unpacked_filepath, sep="\t", col.names=FALSE, row.names=FALSE)
            } else { system(paste0("rm ",unpacked_filepath)) }
        }
    }
    
    print(paste0("Done."))
}
unpack_RBP_eCLIP_peaks(eclip_path="/home/local/ARCS/ak3792/Documents/Research/data/Yeo_eCLIP/RBP_hg19", version="hg19", p_filter=0.05)
rm(features_env); gc()
setup_features_env()
# fullPeak CHDFB
chdfb_burdens_all_2s_full <- burden_analysis2("CHDFB_combined", ssc="SSCFB_combined", desired_tests=desired_tests, include_individual_RBPs=FALSE, peak_types=c("fullPeak"), statistical_test="fisher_exact_2s")
write.csv(chdfb_burdens_all_2s_full[,c(1:10)], file=output_path("chdfb_burdens_all_2s_full.csv"), row.names=FALSE)
best_result_chdfb_all_2s_full <- volcano_plot2("CHDFB", chdfb_burdens_all_2s_full[-c(1:3, which(grepl("TSS40", "overall", chdfb_burdens_all_2s_full$burden))),], bonferroni_p.value_cutoff=0.05/133, highlight_best=FALSE, num_best_hits=20, pval_relative_importance=1000, number_labels=TRUE, x_min=-0.5, x_max=0.5, label_gap=0.5, label_points_max_iterations=20)


#################################################################################################################
# Set up features_env environment that will store named features (filepaths, data frames, or GRanges objects)
#################################################################################################################
setup_features_env <- function(make_fullPeak=FALSE) {
    features_env <- new.env()
    # Add Roadmap histone modification features
    HISTONE_MODIFICATIONS_ANNOTATIONS_PATH = data_path("Roadmap")
    for(histone_modification_file in list.files(path=HISTONE_MODIFICATIONS_ANNOTATIONS_PATH)) { 
        features_env_key = paste0(gsub("-", ".", histone_modification_file))
        features_env[[features_env_key]] <- new.env()
        features_env[[features_env_key]][["peaks"]] <- full_path(HISTONE_MODIFICATIONS_ANNOTATIONS_PATH, histone_modification_file)
        histone_mark_info <- strsplit(features_env_key, "\\.")[[1]]
        if(histone_mark_info[1] %in% get_relevant_roadmap_eids("CHD")) { features_env[[features_env_key]][["notes"]] <- paste0("histone_mark; ",histone_mark_info[2]) } else { features_env[[features_env_key]][["notes"]] <- "" }
    }
    # Add RBP eCLIP features
    RBP_ANNOTATIONS_PATH = data_path("ENCODE_eCLIP_peak") #data_path("Yeo_eCLIP/RBP_hg19") #data_path("ENCODE_eCLIP_peak")
    for(rbp_eclip_file in list.files(path=RBP_ANNOTATIONS_PATH, pattern="\\.bed$")) { 
        features_env_key = paste0(strsplit(rbp_eclip_file, "\\.")[[1]][1:2], collapse=".")
        features_env[[features_env_key]] <- new.env()
        features_env[[features_env_key]][["peaks"]] <- full_path(RBP_ANNOTATIONS_PATH, rbp_eclip_file)
        features_env[[features_env_key]][["notes"]] <- "RBP; default_padding=50"
    }
    features_env[["RBP"]] <- new.env(); features_env[["RBP"]][["notes"]] <- "default_padding=50"
    # Add variant_type features. Function will be executed with eval(parse(text=c)), like in the following example: eval(parse(text="!is_indel(variants$Ref, variants$Alt)"))
    features_env[["is_snv"]] <- new.env(); features_env[["is_snv"]][["function"]] <- "!is_indel(variants$Ref, variants$Alt, verbose=TRUE)"; features_env[["is_snv"]][["notes"]] <- "variant_type"
    features_env[["is_indel"]] <- new.env(); features_env[["is_indel"]][["function"]] <- "is_indel(variants$Ref, variants$Alt, verbose=TRUE)"; features_env[["is_indel"]][["notes"]] <- "variant_type"
    features_env[["is_variant"]] <- new.env(); features_env[["is_variant"]][["function"]] <- "rep('Y', nrow(variants))"; features_env[["is_variant"]][["notes"]] <- "variant_type"
    features_env[["snv_indel"]] <- new.env(); features_env[["snv_indel"]][["function"]] <- "snv_indel <- rep('snv', nrow(variants)); snv_indel[is_indel(variants$Ref, variants$Alt, verbose=TRUE)] <- 'indel'; return(snv_indel)"; features_env[["snv_indel"]][["notes"]] <- ""
    # Add region_type features.
    features_env[["3'UTR"]] <- new.env(); features_env[["3'UTR"]][["peaks"]] <- TES_granges; features_env[["3'UTR"]][["notes"]] <- "region_type"
    features_env[["TSS"]] <- new.env(); features_env[["TSS"]][["peaks"]] <- TSS_granges; features_env[["TSS"]][["notes"]] <- "region_type"
    features_env[["TSS40"]] <- new.env(); features_env[["TSS40"]][["peaks"]] <- TSS40_granges; features_env[["TSS40"]][["notes"]] <- "region_type"
    # Add nearby_genes features.
    features_env[["nearest_gene"]] <- new.env(); features_env[["nearest_gene"]][["peaks"]] <- genebody_granges; features_env[["nearest_gene"]][["notes"]] <- "nearby_gene; annotation_type=nearest"
    # Add gene_set features, using function.
    features_env[["HDE"]] <- new.env(); features_env[["HDE"]][["function"]] <- "variants$nearest_gene %in% HDE_genes"; features_env[["HDE"]][["notes"]] <- "gene_set"
    features_env[["HHE"]] <- new.env(); features_env[["HHE"]][["function"]] <- "variants$nearest_gene %in% HHE_genes"; features_env[["HHE"]][["notes"]] <- "gene_set"
    features_env[["HE"]] <- new.env(); features_env[["HE"]][["function"]] <- "variants$nearest_gene %in% HMHHE_genes"; features_env[["HE"]][["notes"]] <- "gene_set"
    features_env[["constrained"]] <- new.env(); features_env[["constrained"]][["function"]] <- "variants$nearest_gene %in% constrained_genes"; features_env[["constrained"]][["notes"]] <- "gene_set"
    features_env[["constrained HE"]] <- new.env(); features_env[["constrained HE"]][["function"]] <- "variants$nearest_gene %in% HMHHE_genes & variants$nearest_gene %in% constrained_genes"; features_env[["constrained HE"]][["notes"]] <- "gene_set"
    features_env[["bivalent genes"]] <- new.env(); features_env[["bivalent genes"]][["function"]] <- "variants$nearest_gene %in% bivalent_genes"; features_env[["bivalent genes"]][["notes"]] <- "gene_set"
    # Add CADD features, using function.
    features_env[["CADD10"]] <- new.env(); features_env[["CADD10"]][["function"]] <- "!is.na(variants$CADDgt10_Phred) & variants$CADDgt10_Phred != '.' & as.numeric(paste0(variants$CADDgt10_Phred)) > 10"; features_env[["CADD10"]][["notes"]] <- "CADD"
    features_env[["CADD15"]] <- new.env(); features_env[["CADD15"]][["function"]] <- "!is.na(variants$CADDgt10_Phred) & variants$CADDgt10_Phred != '.' & as.numeric(paste0(variants$CADDgt10_Phred)) > 15"; features_env[["CADD15"]][["notes"]] <- "CADD"
    features_env[["CADD20"]] <- new.env(); features_env[["CADD20"]][["function"]] <- "!is.na(variants$CADDgt10_Phred) & variants$CADDgt10_Phred != '.' & as.numeric(paste0(variants$CADDgt10_Phred)) > 20"; features_env[["CADD20"]][["notes"]] <- "CADD"
    features_env[["CADD25"]] <- new.env(); features_env[["CADD25"]][["function"]] <- "!is.na(variants$CADDgt10_Phred) & variants$CADDgt10_Phred != '.' & as.numeric(paste0(variants$CADDgt10_Phred)) > 25"; features_env[["CADD25"]][["notes"]] <- "CADD"
    features_env[["CADD30"]] <- new.env(); features_env[["CADD30"]][["function"]] <- "!is.na(variants$CADDgt10_Phred) & variants$CADDgt10_Phred != '.' & as.numeric(paste0(variants$CADDgt10_Phred)) > 30"; features_env[["CADD30"]][["notes"]] <- "CADD"
    
    set_global(features_env)
    
    # Append .fullPeak features
    if(make_fullPeak) {
        features_to_make_fullPeak <- c(get_features_by_group("H3K36me3"), get_features_by_group("H3K79me2"))
        features_to_make_fullPeak <- features_to_make_fullPeak[grepl("broadPeak",features_to_make_fullPeak)]
        gap_size_to_fill_in = 1000
        padding = gap_size_to_fill_in/2
        sapply_out <- sapply(features_to_make_fullPeak, function(x) {
            print(paste0("Creating fullPeak annotation from ",x,"..."))
            annotations <- features_env[[x]][["peaks"]]
            if(class(annotations)!="GRanges") { # If annotations is not already a GRanges class, it must be a data frame or a filepath from which we will construct a new GRanges object.
                if(class(annotations)!="data.frame") { # If annotations is not already a data frame, it must be a filepath from which we will construct a new data frame and then GRanges object.
                    if (grepl("csv", annotations)) { dat_sep = "," } else { dat_sep = "\t" }
                    annotations <- read.csv(annotations, sep=dat_sep, header=FALSE); colnames(annotations)[1:3] <- c("chromosome", "start", "end")
                }
                annotations <- to_genomic_regions(annotations, labels=genomic_coordinates_to_strings(annotations$chromosome, annotations$start, annotations$end))
            }
            
            annotations_full <- annotations
            start(annotations_full) <- start(annotations_full) - padding
            end(annotations_full) <- end(annotations_full) + padding
            annotations_full <- intersect(annotations_full, genebody_granges)
            annotations_full <- union(annotations_full, annotations) # Same as original broadPeak, but with reasonably sized gaps in genebody filled in
            
            features_env_key <- gsub("broad", "full", x)
            histone_mark_info <- strsplit(features_env_key, "\\.")[[1]]
            features_env[[features_env_key]] <- new.env()
            features_env[[features_env_key]][["peaks"]] <- annotations_full
            features_env[[features_env_key]][["notes"]] <- paste0("histone_mark; ",histone_mark_info[2])
            set_global(features_env)
        })
    }
    
    return(features_env)
}
setup_features_env()

#################################################################################################################
# Returns all available features that belong to the specified group, such as "RBP" (across types of RBPs), "H3K36me3" (across tissue types), etc...Supports padding specification, just like the annotate function.
#################################################################################################################
get_features_by_group <- function(group, features_env=NULL) {
    if (is.null(features_env)) { features_env <- get_global("features_env") }; if (is.null(features_env)) { features_env <- setup_features_env() }
    available_features <- ls(features_env)
    num_available_features = length(available_features)
    available_features_notes <- sapply(available_features, function(x) { return(unlist(strsplit(features_env[[x]][["notes"]], "; *"))) })
    
    padding_specified = FALSE
    padding = suppressWarnings(as.numeric(gsub(".*_([0-9]+)(bp)$","\\1", group))) # Capture padding size in bp, if it was appended to the group name. For example, gsub(".*_([0-9]+)(bp)$","\\1", "RBP_50bp") returns 50
    if(!is.na(padding)) { group = gsub("_[0-9]+bp$", "", group); padding_specified = TRUE }
    
    features_group <- available_features[sapply(1:num_available_features, function(i) { group %in% available_features_notes[[i]] })]
    
    if(padding_specified) { features_group <- paste0(features_group,"_",padding,"bp") }
    return(features_group)
}

#################################################################################################################
# Annotate the given variants (passed in as either data frame, GRanges, or variant_list format) with the specified features (saved in features_env as either GRanges, data frames, or filepaths)
#################################################################################################################
annotate <- function(variants, features, features_env=NULL, variants_granges=NULL, save_globally=TRUE) {
    if(is.null(features_env)) { features_env <- get_global("features_env") }; if (is.null(features_env)) { features_env <- setup_features_env() }
    if(class(variants) == "GRanges") { variants_granges <- variants; variants <- genomic_regions_to_dat(variants) } # handles case where variants are passed in directly as GRanges
    variants <- standardize_colnames(variants)
    if (is.null(variants_granges)) { variants_granges <- to_genomic_regions(variants, chr_colname="Chrom", start_colname="Position", end_colname="Position") }
    num_features = length(features)
    available_features <- ls(features_env)
    available_features_notes <- sapply(available_features, function(x) { return(unlist(strsplit(features_env[[x]][["notes"]], "; *"))) })
    # Function to load annotation GRanges for a given feature and padding size.
    load_annotation <- function(feature, padding=0, save_globally=save_globally) {
        annotations <- features_env[[gsub("_[0-9]+bp$", "", feature)]][["peaks"]]
        if(class(annotations)!="GRanges") { # If annotations is not already a GRanges class, it must be a data frame or a filepath from which we will construct a new GRanges object.
            if(class(annotations)!="data.frame") { # If annotations is not already a data frame, it must be a filepath from which we will construct a new data frame and then GRanges object.
                if (grepl("csv", annotations)) { dat_sep = "," } else { dat_sep = "\t" }
                annotations <- read.csv(annotations, sep=dat_sep, header=FALSE); colnames(annotations)[1:3] <- c("chromosome", "start", "end")
                strand_col <- which(apply(annotations, 2, function(annotations_col) sum(paste0(annotations_col) %in% c("+", "-")) == nrow(annotations)))[1]
                if(!is.na(strand_col)) { colnames(annotations)[strand_col] <- "strand" }
            }
            annotations_labels <- genomic_coordinates_to_strings(annotations$chromosome, annotations$start, annotations$end)
            p.value_col <- which(tolower(colnames(annotations)) %in% c("p.value", "p"))[1]
            if(is.na(p.value_col)) { 
                p.value_col <- which(apply(annotations, 2, function(annotations_col) sum(grepl("[pP]=[0-9]", annotations_col)) == nrow(annotations)))[1] 
                if(!is.na(p.value_col)) { annotations_labels <- gsub(".*\\[[pP]=([0-9\\.eE\\-]+)\\].*","\\1", annotations[,p.value_col]) }
            } else { annotations_labels <- paste0(annotations[,p.value_col]) }
            annotations <- to_genomic_regions(annotations, labels=annotations_labels)
        }
        #print(paste0("Feature: ",feature,", Padding in load_annotation: ",padding))
        if(padding != 0) { start(annotations) <- start(annotations) - padding; end(annotations) <- end(annotations) + padding }
        #print(annotations)
        if(save_globally) { 
            if(is.null(features_env[[paste0(feature)]])) { features_env[[paste0(feature)]] <- new.env(); features_env[[paste0(feature)]][["notes"]] <- ""; }
            features_env[[paste0(feature)]][["peaks"]] <- annotations # Save annotations in features_env as a GRanges object for future efficiency
            set_global(features_env)
        }
        return(annotations)
    }
    # Annotate all desired features
    annotation_successes <- sapply(1:num_features, function(i) {
        feature = features[i]
        print(paste0("Annotating feature ",feature,"...[",i," / ",num_features,"]"))
        if(features[i] %in% colnames(variants)) { variants <<- variants[,-c(feature)] }
        padding_specified = FALSE
        padding = suppressWarnings(as.numeric(gsub(".*_([0-9]+)(bp)$","\\1", feature))) # Capture padding size in bp, if it was appended to the feature name. For example, gsub(".*_([0-9]+)(bp)$","\\1", "RBP_50bp") returns 50
        if(!is.na(padding)) { feature = gsub("_[0-9]+bp$", "", feature); padding_specified = TRUE }
        feature_notes <- unlist(available_features_notes[[feature]])
        # Sets annotation type other than "hit", if it is specified in the feature notes.
        annotation_type = gsub("annotation_type=([a-zA-Z]+)","\\1", feature_notes[grepl("annotation_type=", feature_notes)])[1]
        if(is.null(annotation_type) || is.na(annotation_type)) { annotation_type = "hit" }
        # Sets default padding other than zero, if it is specified in the feature notes.
        default_padding = as.numeric(gsub("default_padding=([0-9]+)","\\1", feature_notes[grepl("default_padding=", feature_notes)]))[1]
        if(is.null(default_padding) || is.na(default_padding)) { default_padding = 0 }
        if(!padding_specified) { padding = default_padding }
        #print(paste0("Padding specified: ",padding_specified))
        #print(paste0("Using padding of ",padding,"bp"))
        # Get annotations for the single feature or feature group.
        if(feature %in% available_features && "function" %in% ls(features_env[[feature]])) {
            annotation_vector <- eval(parse(text=features_env[[feature]][["function"]]))
            if(class(annotation_vector) == "logical") { 
                annotation_vector_char <- rep('N', length(annotation_vector))  
                annotation_vector_char[annotation_vector] <- 'Y'
                annotation_vector <- annotation_vector_char
            }
        } else {
            if(feature %in% available_features && "peaks" %in% ls(features_env[[feature]])) { 
                annotations <- load_annotation(features[i], padding-default_padding, save_globally)
            } else {
                # Scan for possible feature group
                features_group <- available_features[unlist(lapply(available_features_notes, function(x) { return(feature %in% unlist(x)) }))]
                num_features_in_group = length(features_group)
                if(num_features_in_group > 0) {
                    annotations <- GRanges()
                    sapply_out <- sapply(1:num_features_in_group, function(features_group_index) { 
                        features_group_feature = features_group[features_group_index]
                        if(padding_specified) { features_group_feature = paste0(features_group_feature,"_",padding,"bp") } 
                        print(paste0("    Grabbing ",feature," group feature ",features_group_feature," [",features_group_index,"/",num_features_in_group,"]..."))
                        #annotations <<- c(annotations, load_annotation(features_group_feature, padding, save_globally))
                        annotations <<- c(annotations, load_annotation(features_group_feature, padding, FALSE))
                        annotations <<- intersect(annotations, annotations)
                    })
                    if(save_globally) { 
                        if(is.null(features_env[[features[i]]])) { features_env[[features[i]]] <- new.env(); features_env[[features[i]]][["notes"]] <- ""; }
                        features_env[[features[i]]][["peaks"]] <- annotations # Save annotations in features_env as a GRanges object for future efficiency
                        set_global(features_env)
                    }
                } else { print(paste0("Skipping feature ",feature,": not in available features.")); return(FALSE) }
            }
            
            if(annotation_type == "nearest") {
                annotation_vector <- names(annotations)[nearest(variants_granges, annotations)]
            }
            else { # simple Y/N "hit" annotation
                overlaps_with_annotation <- unique(queryHits(findOverlaps(variants_granges, annotations)))
                annotation_vector <- rep("N", nrow(variants))
                annotation_vector[overlaps_with_annotation] <- "Y"
            } 
        }
        variants <<- cbind(variants, annotation_vector)
        colnames(variants)[ncol(variants)] <<- features[i]
        return(TRUE)
    })
    print(paste0("Annotation complete! ",sum(annotation_successes)," out of ",num_features," features successfully annotated."))
    
    return(variants)
}

#################################################################################################################
# Annotates all available RBPs with their gene targets, based on nearest gene to each peak of every RBP
#################################################################################################################
get_RBP_targets <- function(nearby_peaks_required=1, return_target_counts_matrix=FALSE) {
    load_annotation <- function(feature) {
        annotations <- features_env[[gsub("_[0-9]+bp$", "", feature)]][["peaks"]]
        if(class(annotations)!="GRanges") { # If annotations is not already a GRanges class, it must be a data frame or a filepath from which we will construct a new GRanges object.
            if(class(annotations)!="data.frame") { # If annotations is not already a data frame, it must be a filepath from which we will construct a new data frame and then GRanges object.
                if (grepl("csv", annotations)) { dat_sep = "," } else { dat_sep = "\t" }
                annotations <- read.csv(annotations, sep=dat_sep, header=FALSE); colnames(annotations)[1:3] <- c("chromosome", "start", "end")
                strand_col <- which(apply(annotations, 2, function(annotations_col) sum(paste0(annotations_col) %in% c("+", "-")) == nrow(annotations)))[1]
                if(!is.na(strand_col)) { colnames(annotations)[strand_col] <- "strand" }
            }
            annotations_labels <- genomic_coordinates_to_strings(annotations$chromosome, annotations$start, annotations$end)
            p.value_col <- which(tolower(colnames(annotations)) %in% c("p.value", "p"))[1]
            if(is.na(p.value_col)) { 
                p.value_col <- which(apply(annotations, 2, function(annotations_col) sum(grepl("[pP]=[0-9]", annotations_col)) == nrow(annotations)))[1] 
                if(!is.na(p.value_col)) { annotations_labels <- gsub(".*\\[[pP]=([0-9\\.eE\\-]+)\\].*","\\1", annotations[,p.value_col]) }
            } else { annotations_labels <- paste0(annotations[,p.value_col]) }
            annotations <- to_genomic_regions(annotations, labels=annotations_labels)
        }
        return(annotations)
    }
    
    rbps <- get_features_by_group("RBP")
    num_rbps <- length(rbps)
    rbp_targets <- data.frame(t(sapply(1:num_rbps, function(i) {
        rbp <- rbps[i]
        print(paste0("Finding RBP targets for ",rbp, " (",nearby_peaks_required, " nearby peaks required)... [",i," / ",num_rbps,"]"))
        rbp_peaks <- load_annotation(rbp)
        rbp_annotated <- annotate(rbp_peaks, "nearest_gene")
        rbp_target_gene_counts <- table(rbp_annotated[,ncol(rbp_annotated)])
        if(return_target_counts_matrix) {
            all_gene_counts <- rep(0, length(genebody$gene)); names(all_gene_counts) <- sort(genebody$gene)
            all_gene_counts[names(rbp_target_gene_counts)] <- rbp_target_gene_counts
            return(c(rbp, all_gene_counts))
        } else {
            rbp_target_genes <- sort(names(rbp_target_gene_counts)[rbp_target_gene_counts >= nearby_peaks_required])
            return(c(rbp, footprint(rbp_peaks), length(rbp_target_genes), paste0(rbp_target_genes, collapse=",")))
        }
    })))
    if(return_target_counts_matrix) {
        colnames(rbp_targets) <- c("RBP", sort(genebody$gene))
        rbp_targets <- unfactorize(rbp_targets)
    } else {
        colnames(rbp_targets) <- c("RBP", "footprint", "num_targets", "targets")
        rbp_targets$footprint <- as.numeric(paste0(rbp_targets$footprint)); rbp_targets$num_targets <- as.numeric(paste0(rbp_targets$num_targets))
    }
    return(rbp_targets)
}
write.csv(get_RBP_targets(return_target_counts_matrix=TRUE), file=output_path("rbp_target_peak_counts.csv"), row.names=FALSE)
write.csv(get_RBP_targets(), file=output_path("rbp_targets.csv"), row.names=FALSE)
write.csv(get_RBP_targets(2), file=output_path("rbp_targets_2_nearby_peaks_required.csv"), row.names=FALSE)
write.csv(get_RBP_targets(25), file=output_path("rbp_targets_25_nearby_peaks_required.csv"), row.names=FALSE)

nearby_peaks_required_to_test <- c(1:5,10,20,25,50)
sapply_out <- sapply(nearby_peaks_required_to_test, function(nearby_peaks_required) {
    rbp_targets_info <- get_RBP_targets(nearby_peaks_required)[,c(1,2,3)]
    if(nearby_peaks_required == 1) { 
        pdf(file=output_path("rbp_footprint_size_density.pdf"))
        plot(density(rbp_targets_info$footprint), main="", xlab="footprint size (bp)")
        dev.off()
        mtext_msg = "1 nearby peak required to implicate target"; filename_suffix = "_1_peak_required.pdf"
    } else { mtext_msg = paste0(nearby_peaks_required," nearby peaks required to implicate target"); filename_suffix = paste0("_",nearby_peaks_required,"_peaks_required.pdf") }
    pdf(file=output_path(paste0("rbp_target_count_density",filename_suffix)))
    plot(density(rbp_targets_info$num_targets), main="", xlab="# targets")
    mtext(mtext_msg)
    dev.off()
    pdf(file=output_path(paste0("rbp_target_count_vs_footprint_size",filename_suffix)))
    plot(rbp_targets_info$footprint, rbp_targets_info$num_targets, xlab="footprint size (bp)", ylab="# targets")
    mtext(mtext_msg)
    dev.off()
})

load_annotation <- function(feature, padding=0, save_globally=FALSE) {
    annotations <- features_env[[gsub("_[0-9]+bp$", "", feature)]][["peaks"]]
    if(is.null(annotations)) { return(NULL) }
    if(class(annotations)!="GRanges") { # If annotations is not already a GRanges class, it must be a data frame or a filepath from which we will construct a new GRanges object.
        if(class(annotations)!="data.frame") { # If annotations is not already a data frame, it must be a filepath from which we will construct a new data frame and then GRanges object.
            if (grepl("csv", annotations)) { dat_sep = "," } else { dat_sep = "\t" }
            annotations <- read.csv(annotations, sep=dat_sep, header=FALSE); colnames(annotations)[1:3] <- c("chromosome", "start", "end")
            strand_col <- which(apply(annotations, 2, function(annotations_col) sum(paste0(annotations_col) %in% c("+", "-")) == nrow(annotations)))[1]
            if(!is.na(strand_col)) { colnames(annotations)[strand_col] <- "strand" }
        }
        annotations_labels <- genomic_coordinates_to_strings(annotations$chromosome, annotations$start, annotations$end)
        p.value_col <- which(tolower(colnames(annotations)) %in% c("p.value", "p"))[1]
        if(is.na(p.value_col)) { 
            p.value_col <- which(apply(annotations, 2, function(annotations_col) sum(grepl("[pP]=[0-9]", annotations_col)) == nrow(annotations)))[1] 
            if(!is.na(p.value_col)) { annotations_labels <- gsub(".*\\[[pP]=([0-9\\.eE\\-]+)\\].*","\\1", annotations[,p.value_col]) }
        } else { annotations_labels <- paste0(annotations[,p.value_col]) }
        annotations <- to_genomic_regions(annotations, labels=annotations_labels)
    }
    #print(paste0("Feature: ",feature,", Padding in load_annotation: ",padding))
    if(padding != 0) { start(annotations) <- start(annotations) - padding; end(annotations) <- end(annotations) + padding }
    #print(annotations)
    if(save_globally) { 
        if(is.null(features_env[[paste0(feature)]])) { features_env[[paste0(feature)]] <- new.env(); features_env[[paste0(feature)]][["notes"]] <- ""; }
        features_env[[paste0(feature)]][["peaks"]] <- annotations # Save annotations in features_env as a GRanges object for future efficiency
        set_global(features_env)
    }
    return(annotations)
}



#random_gnomad_variants <- get_random_gnomad_variants(1000000)

find_test_hits <- function(annotated_variants, tests) {
    num_tests = length(tests)
    test_elements <- lapply(1:num_tests, function(i) unlist(strsplit(tests[i], "\\+")))
    features = unique(unlist(test_elements))
    num_features = length(features)
    
    annotations <- apply(annotated_variants, 2, function(x) which(x == "Y"))
    
    stored_results <- new.env()
    fill_in_test <- function(test_elements) {
        stored_results_key <- paste0(test_elements, collapse="+")
        #print(stored_results_key)
        stored_result <<- stored_results[[stored_results_key]]
        if(is.null(stored_result)) {
            indices <- unlist(annotations[test_elements[1]])
            if(length(test_elements) > 1) { indices <- intersect(indices, fill_in_test(test_elements[-1])) }
            stored_results[[stored_results_key]] <<- indices
            #print(indices)
            #readline()
            #print(length(indices))
            return(indices)
        } else {
            #print(length(stored_result))
            return(stored_result)
        }
    }
    annotations_hits <- sapply(1:num_tests, function(i) { print(i); fill_in_test(unlist(test_elements[i])) })
    names(annotations_hits) <- tests
    #annotations_hits <- sapply(1:num_tests, function(i) { print(i); Reduce(intersect, annotations[unlist(test_elements[i])]) })
    #annotations_hits <- sapply(1:num_tests, function(i) { print(i); sum(rowSums(annotations[,unlist(test_elements[i]),drop=F]) == 0) })
    return(annotations_hits)
}

find_test_hits <- function(annotated_variants, tests) {
    num_tests = length(tests)
    test_elements <- lapply(1:num_tests, function(i) unlist(strsplit(tests[i], "\\+")))
    features = unique(unlist(test_elements))
    num_features = length(features)
    
    annotations <- apply(annotated_variants, 2, function(x) which(x == "Y"))
    annotations_hits <- sapply(1:num_tests, function(i) { print(i); Reduce(intersect, annotations[unlist(test_elements[i])]) })
    #annotations_hits <- sapply(1:num_tests, function(i) { print(i); sum(rowSums(annotations[,unlist(test_elements[i]),drop=F]) == 0) })
    return(annotations_hits)
}

# Uses Werling eigenvector decomposition/PCA to calculate the number of independent tests.
calculate_num_independent_tests <- function(tests, num_variants_per_iteration=10000, num_cases_per_iteration=NULL, num_controls_per_iteration=NULL, num_iterations=100, random_gnomad_variants=random_gnomad_variants) {
    num_tests <- length(tests)
    #num_variants_half = floor(num_variants_per_iteration/2)
    if(is.null(num_cases_per_iteration)) { num_cases_per_iteration = floor(num_variants_per_iteration/2) }
    if(is.null(num_controls_per_iteration)) { num_controls_per_iteration = floor(num_variants_per_iteration/2) }
    num_variants_per_iteration = num_cases_per_iteration + num_controls_per_iteration
    
    all_features <- unique(unlist(strsplit(tests, "\\+")))
    num_random_gnomad_variants <- nrow(random_gnomad_variants)
    random_gnomad_variants_granges <- to_genomic_regions(random_gnomad_variants, chr_colname="Chrom", start_colname="Position", end_colname="Position")
    #random_gnomad_variants <- annotate(random_gnomad_variants, "nearest_gene")
    random_gnomad_variants_annotated <- annotate(random_gnomad_variants, all_features, variants_granges=random_gnomad_variants_granges)
    random_gnomad_variants_annotated_hits <- find_test_hits(random_gnomad_variants_annotated, tests)
    
    #random_gnomad_samples_all <- sample(1:num_random_gnomad_variants, num_variants_per_iteration*num_iterations, replace=TRUE)
    #random_gnomad_samples_cases_all <- random_gnomad_samples_all[1:(num_cases_per_iteration*num_iterations)]
    #random_gnomad_samples_controls_all <- random_gnomad_samples_all[(length(random_gnomad_samples_all) - (num_controls_per_iteration*num_iterations) + 1):length(random_gnomad_samples_all)]
    random_gnomad_samples_cases_all <- sample(1:num_random_gnomad_variants, num_cases_per_iteration*num_iterations, replace=TRUE)
    random_gnomad_samples_controls_all <- sample(1:num_random_gnomad_variants, num_controls_per_iteration*num_iterations, replace=TRUE)
    
    # Calculate and populate the (#_iterations x #_tests) Z-scores matrix
    z_scores_matrix <- data.frame(t(sapply(1:num_iterations, function(iteration) {
        print(paste0("Iteration ",iteration," / ",num_iterations))
        #annotated_variants <- random_gnomad_variants_annotated[((iteration-1)*num_variants_per_iteration):(iteration*num_variants_per_iteration),]
        #annotated_variants_hits <- find_test_hits(annotated_variants, tests)
        #p.values <- sapply(1:num_tests, function(i) {
        #    print(paste0("Test ",i," / ",num_tests))
        #    test_elements <- unlist(strsplit(tests[i], "\\+"))
        #    num_test_elements = length(test_elements)
        #    annotations <- annotated_variants[,which(colnames(annotated_variants) %in% test_elements)]
        #    if(num_test_elements > 1) { hits <- rowSums(annotations == "Y") == ncol(annotations) } else { hits <- annotations == "Y" }
        #    num_variants_half = num_variants_per_iteration/2
        #    num_hits_cases <- sum(hits[1:num_variants_half])
        #    num_hits_controls <- sum(hits[(num_variants_half + 1):num_variants_per_iteration])
        #    return(fisher_exact_test(num_hits_cases, num_hits_controls, num_variants_half, num_variants_half)$p.value)
        #})
        
        random_gnomad_samples_cases <- table(random_gnomad_samples_cases_all[(((iteration-1)*num_cases_per_iteration)+1):(iteration*num_cases_per_iteration)])
        random_gnomad_samples_controls <- table(random_gnomad_samples_controls_all[(((iteration-1)*num_controls_per_iteration)+1):(iteration*num_controls_per_iteration)])
        
        p.values <- sapply(random_gnomad_variants_annotated_hits, function(x) {
            hits <- unlist(x)
            num_hits_cases <- sum(random_gnomad_samples_cases[hits])
            num_hits_controls <- sum(random_gnomad_samples_controls[hits])
            #print(paste(c(num_hits_cases, num_hits_controls, num_variants_half, num_variants_half), collapse=", "))
            return(fisher_exact_test(num_hits_cases, num_hits_controls, num_cases_per_iteration, num_controls_per_iteration)$p.value)
        })
        z_scores <- (p.values-mean(p.values)) / sd(p.values)
        
        return(z_scores)
    })))
    #dat[] <- lapply(dat, factor)
    #a <- sapply(1:ncol(dat), function(i) { if(length(levels(dat[,i])) > 1) { return(as.numeric(dat[,i])) } else { return(NULL) } } )
    #if(is.list(a)) { 
    #    testable_feature <- sapply(1:length(a), function(i) { return(length(unlist(a[i])) > 0) } ) 
    #    a <- a[testable_feature]
    #} else {
    #    testable_feature <- sapply(1:ncol(a), function(i) { return(length(a[,i]) > 0) } ) 
    #    a <- a[,testable_feature]
    #}
    #a1 <- data.frame(a)
    #rownames(z_scores_matrix) <- NULL
    #colnames(z_scores_matrix) <- NULL
    min_num_iterations = 10
    num_components_to_explain_99_iteration_analysis <- sapply(min_num_iterations:num_iterations, function(iteration_count) {
        z_scores_matrix_subset <- z_scores_matrix[1:iteration_count,,drop=FALSE]
        zero_variance_cols <- which(apply(z_scores_matrix_subset, 2, var)==0)
        if (length(zero_variance_cols) > 0) { z_scores_matrix_subset <- z_scores_matrix_subset[,-c(zero_variance_cols)] } # remove zero-variance columns
        if(ncol(z_scores_matrix_subset) > 1) { 
            pca <- prcomp(z_scores_matrix_subset, center=TRUE, scale.=TRUE)
            percent_variance_explained <- pca$sdev**2/sum(pca$sdev**2)
            num_components_to_explain_99 <- min(which(cumsum(percent_variance_explained) > 0.99))
        } else { num_components_to_explain_99 = 1 }
        return(num_components_to_explain_99)
    })
    pdf(output_path(paste0("num_independent_components_vs_num_iterations_with_",num_variants_per_iteration,"_variants_per_iteration.pdf")))
    plot(min_num_iterations:num_iterations, num_components_to_explain_99_iteration_analysis, type="l", main="Effect of num_iterations on num_independent_components", xlab="num_iterations", ylab="num_independent_components")
    mtext(paste0("Analysis run for ",num_tests," total tests (made from ",length(all_features)," unique features); ",num_variants_per_iteration," variants per iteration"))
    dev.off()
    
    zero_variance_cols <- which(apply(z_scores_matrix, 2, var)==0)
    if (length(zero_variance_cols) > 0) { z_scores_matrix <- z_scores_matrix[,-c(zero_variance_cols)] } # remove zero-variance columns
    if(ncol(z_scores_matrix) > 1) { 
        pca <- prcomp(z_scores_matrix, center=TRUE, scale.=TRUE)
        percent_variance_explained <- pca$sdev**2/sum(pca$sdev**2)
        num_components_to_explain_99 <- min(which(cumsum(percent_variance_explained) > 0.99))
    } else { num_components_to_explain_99 = 1 }
    
    print(paste0("Number of unique features: ",length(all_features)))
    print(paste0("Number of tests (total): ",num_tests))
    print(paste0("Number of independent tests: ",num_components_to_explain_99))
    print(paste0("Bonferonni correction threshold: ",0.05/num_components_to_explain_99))
    
    return_env <- new.env()
    return_env[["num_independent_tests"]] <- num_components_to_explain_99
    return_env[["num_components_to_explain_99_iteration_analysis"]] <- num_components_to_explain_99_iteration_analysis
    return_env[["z_scores"]] <- z_scores_matrix
    return(return_env)
}
#desired_tests <- c("RBP", "E008.H3K36me3", "E003.H3K36me3", "RBP+E008.H3K36me3", "RBP+E003.H3K36me3", "is_indel+RBP+E008.H3K36me3") # test
#num_independent_tests = calculate_num_independent_tests(desired_tests, num_variants_per_iteration=10000, num_iterations=100, random_gnomad_variants=random_gnomad_variants)

####################################################################################################################
# DEFINE DESIRED TESTS!
####################################################################################################################
# Get footprint info of desired histone marks
histone_mark_footprints <- data.frame(t(sapply(c(get_features_by_group("H3K79me2"), get_features_by_group("H3K36me3")), function(desired_histone_mark_feature) {
    print(desired_histone_mark_feature)
    histone_mark_feature <- load_annotation(desired_histone_mark_feature, save_globally=FALSE)
    histone_feature_split <- strsplit(paste0(desired_histone_mark_feature), "\\.")[[1]]
    eid = histone_feature_split[1]
    tissue = get_roadmap_epigenome_names(eid)
    return(c(desired_histone_mark_feature, tissue, sum(width(histone_mark_feature))))
})))
colnames(histone_mark_footprints) <- c("histone_mark", "tissue", "footprint"); rownames(histone_mark_footprints) <- NULL
histone_mark_footprints$footprint <- as.numeric(paste0(histone_mark_footprints$footprint))
write.csv(histone_mark_footprints[grepl("fullPeak", histone_mark_footprints$histone_mark),], file=output_path("histone_mark_footprints.csv"), row.names=FALSE)

pdf(output_path("histone_mark_footprint_distribution.pdf"))
hist(histone_mark_footprints[grepl("fullPeak", histone_mark_footprints$histone_mark),]$footprint, breaks=8, main="Histone mark footprint distribution", xlab="footprint (bp)")
dev.off()

desired_variant_type_features <- get_features_by_group("variant_type") # is_snv, is_indel, and is_variant
desired_region_type_features <- c("", "3'UTR", "TSS") # 3'UTR, TSS, and TSS40
desired_gene_set_features <- c("", "constrained") #c("", "HHE", "HE", "constrained", "constrained HE")  #c("", get_features_by_group("gene_set")) # HDE and HHE
desired_RBP_features <- c("RBP") # combined RBP peaks      # plus each individual RBP feature, get_features_by_group("RBP")
#desired_RBP_features <- c(desired_RBP_features, get_features_by_group("RBP")) # Add individual RBPs
desired_histone_mark_features <- c("", get_features_by_group("H3K79me2"), get_features_by_group("H3K36me3")) #, get_features_by_group("H3K4me1"), get_features_by_group("H3K27ac"), get_features_by_group("H3K9me3")) # H3K36me3 in all relevant tissue types
desired_histone_mark_features <- c(desired_histone_mark_features[grepl("H3K36me3.broadPeak", desired_histone_mark_features)])
desired_CADD_features <- c("", get_features_by_group("CADD"))  

desired_tests <- apply(expand.grid(desired_variant_type_features, desired_region_type_features, desired_gene_set_features, desired_RBP_features, desired_histone_mark_features), 1, function(x) paste0(x[x != ""], collapse="+"))
#desired_tests <- unique(c(apply(expand.grid(desired_variant_type_features, desired_region_type_features, desired_gene_set_features, desired_RBP_features, desired_CADD_features), 1, function(x) paste0(x[x != ""], collapse="+")), desired_tests))
num_tests = length(desired_tests)
all_features <- unique(unlist(strsplit(desired_tests, "\\+")))
print(paste0(num_tests," tests, composed of ",length(all_features)," total unique features."))

library("corrplot")
source("http://www.sthda.com/upload/rquery_cormat.r") # load rquery.cormat(mydata) function
#random_gnomad_variants <- annotate(random_gnomad_variants, "nearest_gene")

num_variants_per_iteration_to_try <- c(100, 500, 1000, 2000, 5000, 10000, 20000, 50000)
# Get analysis results
full_num_independent_tests_analysis <- sapply(1:length(num_variants_per_iteration_to_try), function(i) {
    independent_tests <- calculate_num_independent_tests(desired_tests, num_variants_per_iteration=num_variants_per_iteration_to_try[i], num_iterations=1000, random_gnomad_variants=random_gnomad_variants)
    return(independent_tests[["num_components_to_explain_99_iteration_analysis"]])
})

# Plot analysis
pdf(output_path(paste0("num_independent_components_full_analysis.pdf")))
colors <- rainbow(length(num_variants_per_iteration_to_try))
# Create an empty plot.
plot(1, type = 'n', xlim=c(0, 1000), ylim=c(0, max(full_num_independent_tests_analysis)), main="Effect of num_iterations on num_independent_components", xlab="num_iterations", ylab="num_independent_components", cex.lab = 1.7)
mtext(paste0("Analysis run for ",num_tests," total tests (made from ",length(all_features)," unique features)"))
# Plot lines
sapply_out <- sapply(1:length(num_variants_per_iteration_to_try), function(i) {
    lines(10:1000, full_num_independent_tests_analysis[,i], col=colors[i], lwd=2)
})
legend("topleft", legend = paste0("variants/iteration = ",num_variants_per_iteration_to_try), col = colors, lwd = 1, cex = 0.8)
dev.off()


independent_tests <- 104 #calculate_num_independent_tests(desired_tests, num_cases_per_iteration=nrow(chdfb_combined_muts_wgsa), num_controls_per_iteration=nrow(sscfb_combined_muts_wgsa), num_iterations=1000, random_gnomad_variants=random_gnomad_variants[1:1000000,])

#num_independent_tests
#z_scores <- independent_tests[["z_scores"]]




cormat <- function(mydata) {
    mydata <- unfactorize(mydata[, -which(colnames(mydata) %in% c("Chrom","Position","Ref","Alt","nearest_gene"))])
    mydata[mydata == "N"] <- 0
    mydata[mydata == "Y"] <- 1
    for(i in 1:ncol(mydata)) { mydata[,i] <- as.numeric(mydata[,i]) }
    zero_variance_cols <- which(apply(mydata, 2, var)==0)
    mydata <- mydata[,-zero_variance_cols]
    
    rquery.cormat(mydata[,1:20])
}
col<- colorRampPalette(c("blue", "white", "red"))(20)
heatmap(x = as.matrix(mydata[,c("HHE", "HDE", "TSS", "3'UTR", "RBP", "E008.H3K36me3.broadPeak", "E003.H3K36me3.broadPeak")]), col = col)




volcano_plot2 <- function(disease, burdens, burden_colors=NULL, num_best_hits=0, strong_effect_cutoff=1, bonferroni_p.value_cutoff=-1, highlight_best=TRUE, requested_labels=c(), requested_labels_description="Requested results", filename=NULL, copy_to_images=TRUE, number_labels=TRUE, pval_relative_importance=300, label_cex = 0.7, label_gap=1.1, label_points_max_iterations=200, x_min=NULL, x_max=NULL, y_max=NULL) {
    cat(paste0("Making volcano plot for ", disease, "..."))
    x_max_provided = !is.null(x_max); x_min_provided = !is.null(x_min)
    marginal_p.value_cutoff = 0.05
    stronger_p.value_cutoff = 0.001
    burdens <- burdens[!grepl("RBP0|candidate|HMHDE|HMHHE", burdens$burden),]
    burdens$enrichment <- as.numeric(paste0(burdens$enrichment))
    burdens$p.value <- as.numeric(paste0(burdens$p.value))
    num_total_tests = nrow(burdens)
    num_marginally_significant = sum(burdens$p.value < marginal_p.value_cutoff)
    burdens <- burdens[burdens$enrichment < Inf & burdens$enrichment > 0,] 
    #burdens <- burdens[burdens$p.value < marginal_p.value_cutoff,]
    label <- paste0(gsub(" \\(.+\\)","",burdens$burden, perl=TRUE), " (", burdens$variants, ")")
    burdens <- cbind(burdens, label)
    if (grepl("CDH", disease)) { requested_labels <- gsub("HE ", "DE ", requested_labels); requested_labels <- gsub("CHD", "CDH", requested_labels)
    if(!is.null(burden_colors)) { burden_colors$label <- gsub("HE ", "DE ", burden_colors$label); burden_colors$label <- gsub("CHD", "CDH", burden_colors$label) }
    } else if (grepl("CHD", disease)) { requested_labels <- gsub("DE ", "HE ", requested_labels); requested_labels <- gsub("CDH", "CHD", requested_labels)
    if(!is.null(burden_colors)) { burden_colors$label <- gsub("DE ", "HE ", burden_colors$label); burden_colors$label <- gsub("CDH", "CHD", burden_colors$label) } 
    }
    if (length(requested_labels) > 0) { burdens <- burdens[c(which(!(burdens$label %in% requested_labels)), which((burdens$label %in% requested_labels))),] }
    if(length(burden_colors) < 1) { burden_colors <- NULL }
    
    cols <- rep("gray50", nrow(burdens))
    if (is.null(burden_colors)) {
        significant_bonferroni <- (burdens$p.value < bonferroni_p.value_cutoff)
        significant_stronger <- (burdens$p.value < stronger_p.value_cutoff & !significant_bonferroni)
        significant_marginal <- (burdens$p.value < marginal_p.value_cutoff & !significant_stronger & !significant_bonferroni)
        requested <- burdens$label %in% requested_labels
        strong_effect <- (burdens$enrichment > strong_effect_cutoff)
        
        ##cols[requested] <- "black"
        #cols[strong_effect] <- "gray50"
        #cols[requested & strong_effect] <- "gray50"
        #cols[significant_marginal] <- "gray50"
        #cols[requested & significant_marginal] <- "gray50"
        #cols[significant_marginal & strong_effect] <- "gray50" #"orange"
        ##cols[requested & significant_marginal & strong_effect] <- "orange4"
        cols[significant_stronger & strong_effect] <- "gray50" #"darkorange1"
        cols[significant_bonferroni & strong_effect] <- "red"
    } else {
        burden_colors <- data.frame(cbind(paste0(burden_colors$label), 1:nrow(burden_colors))); colnames(burden_colors) <- c("label","rank"); rownames(burden_colors) <- paste0(burden_colors$label)
        color_rank <- burden_colors[paste0(burdens$label),]; color_rank$rank <- as.numeric(paste0(color_rank$rank))
        col_ramp <- colorRampPalette(brewer.pal(11,"RdBu"), bias=3)(sum(!is.na(color_rank$rank)))
        cols[!is.na(color_rank$rank)] <- col_ramp[color_rank$rank[!is.na(color_rank$rank)]]
    }
    
    if(is.null(filename)) { 
        filename = paste0(tolower(disease),"_burden_volcano_plot.pdf"); filename = gsub(" *\\( *", "_", filename); filename = gsub(" *\\) *", "", filename);
        if(!is.null(burden_colors)) { filename <- gsub("plot", "heat_plot", filename) }
    }
    
    if(number_labels) { pdf(file=output_path(filename), width=14)
    } else { pdf(file=output_path(filename)) }
    
    if(!x_min_provided) { x_min = min(log2(burdens$enrichment))-0.5 } #0 # Min value for x axis
    if(!x_max_provided) { x_max = max(log2(burdens$enrichment))+0.5 } #18 # Max value for x axis
    if(is.null(y_max)) { y_max = ceiling(-log10(min(as.numeric(paste0(burdens$p.value))))) } #min(c(5, ceiling(-log10(min(as.numeric(paste0(burdens$p.value)))))))
    cex.axis=1.6; cex.lab=1.8; cex.main=1.4; cex.mtext=1.1
    par_mar <- par()$mar
    par_pin <- par()$pin
    if(number_labels) {
        par(mfrow=c(1,2), mar=c(par_mar[1]*2, par_mar[2]*2, par_mar[2]/1.5, 0))
        if(!x_max_provided) { x_max = max(log2(burdens$enrichment))+0.5 }
    }
    y_min = 0 #-log10(0.05)
    
    col_alpha = 0.4
    cols <- adjustcolor(cols, alpha.f=col_alpha)
    plot(log2(burdens$enrichment), -log10(burdens$p.value), pch=19, col=cols, main=paste0(disease, " Burden Volcano Plot"), xlab="log2(Odds Ratio)", ylab="-log10(p.value)", xlim=c(x_min,x_max), ylim=c(y_min,y_max), xaxt="n", xaxs="i", yaxs="i", cex.lab=cex.lab, cex.axis=cex.axis, cex.main=cex.main)
    axis(1, at=seq((floor(2*min(log2(burdens$enrichment)))/2), (ceiling(x_max*2)/2), by=0.5), cex.axis=cex.axis, cex.lab=cex.lab)
    library(Hmisc)
    minor.tick(nx=5, ny=5, tick.ratio=0.5)
    left_edge = par("usr")[1]; bottom_edge = par("usr")[3]
    abline(v=0, col="black", lty=2, lwd=1.5); 
    #abline(v=log2(strong_effect_cutoff), col="red", lty=2); #text(log2(strong_effect_cutoff)+0.05, bottom_edge+0.02, adj=c(0,0), paste0(strong_effect_cutoff,"x"), col="red", cex=0.6)
    abline(h=-log10(bonferroni_p.value_cutoff), col="red", lty=2, lwd=1.5); text(left_edge+0.05, -log10(bonferroni_p.value_cutoff)+0.03, adj=c(0,0), paste0("p=", formatC(bonferroni_p.value_cutoff,format="e",digits=2),""), col="red", cex=1) 
    #abline(h=-log10(stronger_p.value_cutoff), col="red", lty=2); text(left_edge+0.05, -log10(stronger_p.value_cutoff)+0.03, adj=c(0,0), paste0("p=",stronger_p.value_cutoff), col="red", cex=0.6)
    #abline(h=-log10(marginal_p.value_cutoff), col="red", lty=2); text(left_edge+0.05, -log10(marginal_p.value_cutoff)+0.03, adj=c(0,0), paste0("p=",marginal_p.value_cutoff), col="red", cex=0.6)
    #mtext("Marginally significant (p < 0.05) tests") #mtext(paste0(num_total_tests, " total tests, ", num_marginally_significant, " (", round(100*num_marginally_significant/num_total_tests,1),"%) pass marginal significance"), cex=cex.mtext)
    if (is.null(burden_colors)) { best_results_col = "gray20" # "red4"
    } else { best_results_col = "black" }
    if(length(requested_labels) > 0 && is.null(burden_colors)) {
        #legend("topleft", legend=c(paste0("Top ", num_best_hits, " results"), requested_labels_description), col=c(best_results_col,"black"), pch=15, cex=0.8)
    } else { #legend("topleft", legend=c(paste0("Top ", num_best_hits, " results")), col=c(best_results_col), pch=15, cex=0.8) }
    }
    burdens <- burdens[order(pval_relative_importance*(-log10(burdens$p.value)/mean(-log10(burdens$p.value))) + (log2(burdens$enrichment)/mean(log2(burdens$enrichment))), decreasing=TRUE),]
    
    cat("Done.\n")
    if (length(requested_labels) > 0) {
        requested_results <- burdens[burdens$label %in% requested_labels,]
    } else { requested_results <- c() }
    if(num_best_hits > 0) {
        best_results <- burdens[0:num_best_hits,] #[burdens$enrichment > strong_effect_cutoff & burdens$p.value < marginal_p.value_cutoff,][0:num_best_hits,]
        cat("Top results: \n"); cat(paste0(1:num_best_hits, ". ", best_results$label, ":   ", round(best_results$enrichment,1), "x enrichment, p=", round(best_results$p.value,4), "\n"))
        best_results_cols <- rep(best_results_col,num_best_hits); best_results_cols[best_results$label %in% requested_labels] <- "red3"
        requested_results <- requested_results[which(!(requested_results$label %in% best_results$label)),]; requested_labels <- requested_labels[!(requested_labels %in% best_results$label)]
        #number_labels_legend <- label_points(log2(c(best_results$enrichment,requested_results$enrichment)), -log10(c(best_results$p.value,requested_results$p.value)), cex=label_cex, c(paste0(best_results$label),paste0(requested_results$label)), col=c(best_results_cols,rep("black",length(requested_results))), number_labels=number_labels, gap=label_gap, max_iterations=label_points_max_iterations)
    } else {
        #number_labels_legend <- label_points(log2(requested_results$enrichment), -log10(requested_results$p.value), cex=label_cex, paste0(requested_results$label), col="black", number_labels=number_labels, gap=label_gap, max_iterations=label_points_max_iterations)
    }
    
    if(number_labels) {
        par(mar=c(par_mar[1], par_mar[1]/7, par_mar[3], 0))
        #par(pin=c(par_pin[1]*1.5, par_pin[2]))
        plot.new()
        #legend("topleft", legend=number_labels_legend, col=c("white"), pch=26, cex=0.8, bty="n")
        library("plotrix")
        highlighted_tests <- burdens[burdens$enrichment > strong_effect_cutoff & burdens$p.value < stronger_p.value_cutoff, c("burden", "variants", "m1", "m0", "enrichment", "p.value")]
        if(nrow(highlighted_tests)>0) {  
            colnames(highlighted_tests) <- c("feature_combination", "variants", "hits_cases", "hits_controls", "OR", "p.value")
            highlighted_tests$OR <- round(highlighted_tests$OR, 2); highlighted_tests$p.value <- formatC(highlighted_tests$p.value, format="e", digits=2)
            highlighted_tests$variants <- gsub("SNP", "SNV", highlighted_tests$variants)
            highlighted_tests$feature_combination <- gsub(" \\(E.+Peak\\)", "", highlighted_tests$feature_combination)
            
            best_rbp_h3k36me3_hit <- which(grepl("RBP H3K36me3", highlighted_tests$test))[1]; if(is.na(best_rbp_h3k36me3_hit)) { best_rbp_h3k36me3_hit = -1 }
            light_dark <- c("gray87", "gray95")
            table_bg_cols <- t(sapply(1:nrow(highlighted_tests), function(i) { if(i == best_rbp_h3k36me3_hit && highlight_best) { bg_col = "lightgoldenrodyellow" } else { bg_col = light_dark[(i %% 2) + 1] }; return(rep(bg_col, ncol(highlighted_tests))) }))
            #table_text_cols <- t(sapply(1:nrow(highlighted_tests), function(i) { if(i %in% best_rbp_h3k36me3_hit) { text_col = "red" } else { text_col = "black" }; return(rep(text_col, ncol(highlighted_tests))) }))
            #matrix(rep(rep("gray90",ncol=ncol(highlighted_tests)),nrow(highlighted_tests)), nrow=nrow(highlighted_tests), ncol=ncol(highlighted_tests))
            max_rows = 20
            addtable2plot("top", table=highlighted_tests[1:min(c(nrow(highlighted_tests),max_rows)),], hlines=FALSE, vlines=TRUE, bg=table_bg_cols, text.col="black", cex=0.85, yjust=-1, xpad=0.15, ypad=1)
            
            #if(!is.na(best_rbp_h3k36me3_hit)) {
            #    highlighted_tests[-c(best_rbp_h3k36me3_hit),] <- ""
            #    addtable2plot("topleft", table=highlighted_tests, hlines=FALSE, vlines=TRUE, bg=table_bg_cols, text.col="red", cex=0.7, yjust=-1, xpad=0.15, ypad=1)
            #}
            
            mtext("Highlighted Tests", font=2)
            
            highlighted_tests_filename = gsub("volcano_plot.pdf", "highlighted_tests.csv", filename)
            write.csv(highlighted_tests, file=output_path(highlighted_tests_filename), row.names=FALSE)
            
            #library(gridExtra)
            #library(grid)
            #g <- tableGrob(iris[1:4, 1:3])
            #find_cell <- function(table, row, col, name="core-fg"){
            #    l <- table$layout
            #    which(l$t==row & l$l==col & l$name==name)
            #}
            
            #ind <- find_cell(g, 3, 2, "core-fg")
            #ind2 <- find_cell(g, 2, 3, "core-bg")
            #g$grobs[ind][[1]][["gp"]] <- gpar(fontsize=15, fontface="bold")
            #g$grobs[ind2][[1]][["gp"]] <- gpar(fill="darkolivegreen1", col = "darkolivegreen4", lwd=5)
            #grid.draw(g)
            ##grid.table(highlighted_tests, base_size = 8)
        }
    }
    
    dev.off()
    par(mfrow=c(1,1), mar=c(par_mar), pin=c(par_pin))
    
    if(copy_to_images) { pdf_to_png(output_path(filename)) }
    
    return(burdens)
}
best_result_cdh3_2s_full <- volcano_plot2("CDH3", cdh3_burdens[-c(1:3, which(grepl("TSS40", "overall", cdh3_burdens$burden))),], bonferroni_p.value_cutoff=0.05/105, highlight_best=FALSE, num_best_hits=20, pval_relative_importance=1000, number_labels=TRUE, x_min=-0.5, x_max=0.8, label_gap=0.5, label_points_max_iterations=20)
best_result_cdh1_cdh2_2s_full <- volcano_plot2("CDH1+CDH2", cdh1_cdh2_burdens[-c(1:3, which(grepl("TSS40", "overall", cdh1_cdh2_burdens$burden))),], bonferroni_p.value_cutoff=0.05/105, highlight_best=FALSE, num_best_hits=20, pval_relative_importance=1000, number_labels=TRUE, x_min=-0.5, x_max=0.5, label_gap=0.5, label_points_max_iterations=20)

write.csv(cdh3_burdens, file=output_path("cdh3_burdens.csv"), row.names=FALSE)
write.csv(cdh1_cdh2_burdens, file=output_path("cdh1_cdh2_burdens.csv"), row.names=FALSE)
write.csv(cdh_burdens_all, file=output_path("cdh_burdens_all.csv"), row.names=FALSE)

best_result_chdfb_all_2s_full <- volcano_plot2("CHDFB", chdfb_burdens_all_2s_full[-c(1:3, which(grepl("TSS40", "overall", chdfb_burdens_all_2s_full$burden))),], bonferroni_p.value_cutoff=0.05/105, highlight_best=FALSE, num_best_hits=20, pval_relative_importance=1000, number_labels=TRUE, x_min=-0.5, x_max=0.5, y_max=4.2, label_gap=0.5, label_points_max_iterations=20)

#chdfb_burdens_all_2s <- burden_analysis2("CHDFB_combined", ssc="SSCFB_combined", include_individual_RBPs=TRUE, statistical_test="fisher_exact_2s")
best_result_chdfb_all_2s <- volcano_plot2("CHDFB", chdfb_burdens_all_2s[-c(1:3),], bonferroni_p.value_cutoff=0.05/97, num_best_hits=20, pval_relative_importance=1000, number_labels=TRUE, x_min=-4.5, x_max=4.5, label_gap=0.5, label_points_max_iterations=20)

#cdh_burdens_all_2s <- burden_analysis2("CDH_combined", ssc="SSC_all", include_individual_RBPs=TRUE, statistical_test="fisher_exact_2s")
best_result_cdh_all_2s <- volcano_plot2("CDH", cdh_burdens_all_2s[-c(1:3),], bonferroni_p.value_cutoff=0.05/97, num_best_hits=20, pval_relative_importance=1000, number_labels=TRUE, x_min=-3.5, x_max=3.5, label_gap=0.5, label_points_max_iterations=20)

#cdh_chd_burdens_all_2s <- burden_analysis2("CDH_CHD", ssc="SSC_all", include_individual_RBPs=TRUE, statistical_test="fisher_exact_2s")
best_result_cdh_chd_all_2s <- volcano_plot2("CDH+CHD", cdh_chd_burdens_all_2s[-c(1:3),], bonferroni_p.value_cutoff=0.05/97, num_best_hits=20, pval_relative_importance=1000, number_labels=TRUE, x_min=-5, x_max=5, label_gap=0.5, label_points_max_iterations=20)


#chdfb_burdens_all_2s_broad <- burden_analysis2("CHDFB_combined", ssc="SSCFB_combined", include_individual_RBPs=TRUE, peak_type="broadPeak", statistical_test="fisher_exact_2s")
#write.csv(chdfb_burdens_all_2s_broad, file=output_path("chdfb_burdens_all_2s_broad.csv"), row.names=FALSE)
best_result_chdfb_all_2s_broad <- volcano_plot2("CHD", chdfb_burdens_all_2s_broad[-c(1:3, which(grepl("TSS40", "overall", chdfb_burdens_all_2s_broad$burden))),], bonferroni_p.value_cutoff=0.05/97, num_best_hits=20, pval_relative_importance=1000, number_labels=TRUE, x_min=-4.5, x_max=4.5, label_gap=0.5, label_points_max_iterations=20)

cdh_burdens_all_2s_broad <- burden_analysis2("CDH_combined", ssc="SSC_all", include_individual_RBPs=TRUE, peak_type="broadPeak", statistical_test="fisher_exact_2s")
write.csv(cdh_burdens_all_2s_broad, file=output_path("cdh_burdens_all_2s_broad.csv"), row.names=FALSE)
best_result_cdh_all_2s_broad <- volcano_plot2("CDH", cdh_burdens_all_2s_broad[-c(1:3, which(grepl("TSS40", chdfb_burdens_all_2s_broad$burden))),], bonferroni_p.value_cutoff=0.05/97, num_best_hits=20, pval_relative_importance=1000, number_labels=TRUE, x_min=-3.5, x_max=3.5, label_gap=0.5, label_points_max_iterations=20)

cdh_chd_burdens_all_2s_broad <- burden_analysis2("CDH_CHD", ssc="SSC_all", include_individual_RBPs=TRUE, peak_type="broadPeak", statistical_test="fisher_exact_2s")
write.csv(cdh_chd_burdens_all_2s_broad, file=output_path("cdh_chd_burdens_all_2s_broad.csv"), row.names=FALSE)
best_result_cdh_chd_all_2s_broad <- volcano_plot2("CDH+CHD", cdh_chd_burdens_all_2s_broad[-c(1:3, which(grepl("TSS40", chdfb_burdens_all_2s_broad$burden))),], bonferroni_p.value_cutoff=0.05/97, num_best_hits=20, pval_relative_importance=1000, number_labels=TRUE, x_min=-5, x_max=5, label_gap=0.5, label_points_max_iterations=20)

# CDH3
best_result_chd3_2s_full <- volcano_plot2("CDH3", cdh3_burdens[-c(1:3, which(grepl("TSS40", "overall", chdfb_burdens_all_2s_full$burden))),], bonferroni_p.value_cutoff=0.05/105, highlight_best=FALSE, num_best_hits=20, pval_relative_importance=1000, number_labels=TRUE, x_min=-0.5, x_max=0.5, label_gap=0.5, label_points_max_iterations=20)

# fullPeak CHDFB
chdfb_burdens_all_2s_full <- burden_analysis2("CHDFB_combined", ssc="SSCFB_combined", variant_type_testing=TRUE, desired_tests=desired_tests, include_individual_RBPs=FALSE, peak_types=c("fullPeak"), statistical_test="fisher_exact_2s")
write.csv(chdfb_burdens_all_2s_full[,c(1:10,16,17)], file=output_path("chdfb_burdens_all_2s_full.csv"), row.names=FALSE)
best_result_chdfb_all_2s_full <- volcano_plot2("CHDFB", chdfb_burdens_all_2s_full[-c(1:3, which(grepl("TSS40", "overall", chdfb_burdens_all_2s_full$burden))),], bonferroni_p.value_cutoff=0.05/105, highlight_best=FALSE, num_best_hits=20, pval_relative_importance=1000, number_labels=TRUE, x_min=-0.5, x_max=0.5, y_max=4.2, label_gap=0.5, label_points_max_iterations=20)

# fullPeak CHD w/o FreeBayes
chd_burdens_all_2s_full <- burden_analysis2("CHD_combined", ssc="SSC_all", desired_tests=desired_tests, include_individual_RBPs=FALSE, peak_types=c("fullPeak"), statistical_test="fisher_exact_2s")
write.csv(chd_burdens_all_2s_full[,c(1:10)], file=output_path("chd_burdens_all_2s_full.csv"), row.names=FALSE)
best_result_chd_all_2s_full <- volcano_plot2("CHD", chd_burdens_all_2s_full[-c(1:3, which(grepl("TSS40", "overall", chd_burdens_all_2s_full$burden))),], bonferroni_p.value_cutoff=0.05/133, highlight_best=FALSE, num_best_hits=20, pval_relative_importance=1000, number_labels=TRUE, x_min=-0.5, x_max=0.5, label_gap=0.5, label_points_max_iterations=20)

# fullPeak CDH w/o FreeBayes
cdh_burdens_all_2s_full <- burden_analysis2("CDH_combined", ssc="SSC_all", desired_tests=desired_tests, include_individual_RBPs=FALSE, peak_types=c("fullPeak"), statistical_test="fisher_exact_2s")
write.csv(cdh_burdens_all_2s_full[,c(1:10)], file=output_path("cdh_burdens_all_2s_full.csv"), row.names=FALSE)
best_result_cdh_all_2s_full <- volcano_plot2("CDH", cdh_burdens_all_2s_full[-c(1:3, which(grepl("TSS40", "overall", cdh_burdens_all_2s_full$burden))),], bonferroni_p.value_cutoff=0.05/133, highlight_best=FALSE, num_best_hits=20, pval_relative_importance=1000, number_labels=TRUE, x_min=-0.5, x_max=0.5, label_gap=0.5, label_points_max_iterations=20)

# fullPeak CDH+CHD w/o FreeBayes
cdh_burdens_all_2s_full <- burden_analysis2("CDH_combined", ssc="SSC_all", desired_tests=desired_tests, include_individual_RBPs=FALSE, peak_types=c("fullPeak"), statistical_test="fisher_exact_2s")
write.csv(cdh_burdens_all_2s_full[,c(1:10)], file=output_path("cdh_burdens_all_2s_full.csv"), row.names=FALSE)
best_result_cdh_all_2s_full <- volcano_plot2("CDH", cdh_burdens_all_2s_full[-c(1:3, which(grepl("TSS40", "overall", cdh_burdens_all_2s_full$burden))),], bonferroni_p.value_cutoff=0.05/133, highlight_best=FALSE, num_best_hits=20, pval_relative_importance=1000, number_labels=TRUE, x_min=-0.5, x_max=0.5, label_gap=0.5, label_points_max_iterations=20)

# Get target expression of the given peaks
target_expression <- function(peaks) {
    targets <- unique(annotate(genomic_regions_to_dat(peaks, chr_colname="Chrom", pos_colname="Position"), "nearest_gene")$nearest_gene)
    expressions <- heart_expression_ranks_dat[heart_expression_ranks_dat[,"gene"] %in% targets,"heart_rank"]
    expressions <- as.numeric(expressions[expressions != "."])
    return(expressions)
}

# Print top_burden_result_variants to file, for both CHDFB and SSCFB
rbp_info <- unfactorize(data.frame(t(sapply(get_features_by_group("RBP"), function(rbp) {
    print(rbp)
    rbp_name_split <- unlist(strsplit(rbp,"\\."))
    gene = rbp_name_split[2]
    rbp_pli <- pLI_scores[[gene]]; if(is.null(rbp_pli)) { rbp_pli = "." } else { rbp_pli = round(rbp_pli, 3) }
    rbp_heart_rank <- heart_expression_ranks_dat[heart_expression_ranks_dat[,"gene"] == gene,"heart_rank"][1]
    if(is.na(rbp_heart_rank)) { rbp_heart_rank= "." }
    rbp_granges <- load_annotation(rbp, save_globally=FALSE)
    target_expressions <- target_expression(rbp_granges)
    rbp_targets_heart_expression <- median(target_expressions)
    fraction_hhe_targets <- sum(target_expressions < 0.25)/length(target_expressions)
    return(c(gene, rbp_name_split[1], footprint(rbp_granges), rbp_pli, rbp_heart_rank, rbp_targets_heart_expression, fraction_hhe_targets))
})))); rownames(rbp_info) <- NULL; colnames(rbp_info) <- c("RBP", "cell_line", "footprint", "pLI", "heart_rank", "median_target_heart_rank", "fraction_hhe_targets")

rbp_info_no_blanks <- rbp_info[rbp_info$heart_rank != ".",]
rbp_heart_ranks <- as.numeric(rbp_info_no_blanks[,"heart_rank"])
rbp_median_target_heart_ranks <- as.numeric(rbp_info_no_blanks[,"median_target_heart_rank"])
rbp_fraction_hhe_targets <- as.numeric(rbp_info_no_blanks[,"fraction_hhe_targets"])

pdf(output_path("rbp_heart_expression_densities.pdf"))
plot(density(as.numeric(heart_expression_ranks_dat[,"heart_rank"]), from=0, to=1), xlim=c(0,1), ylim=c(0,11), xaxs="i", yaxs="i", xlab="E14.5 mouse heart rank", main="RBP heart expression densities", cex.main=1.3, cex.lab=1.4, cex.axis=1.4, lwd=2)
lines(density(rbp_heart_ranks, from=0, to=1), col="blue", lwd=2)
abline(v=as.numeric(rbp_info_no_blanks[rbp_info_no_blanks[,"RBP"] == "RBFOX2","heart_rank"]), lty=2, col="blue")
abline(v=as.numeric(rbp_info_no_blanks[rbp_info_no_blanks[,"RBP"] == "QKI","heart_rank"]), lty=3, col="blue")
lines(density(rbp_median_target_heart_ranks, from=0, to=1), col="red", lwd=2)
abline(v=mean(as.numeric(rbp_info_no_blanks[rbp_info_no_blanks[,"RBP"] == "RBFOX2","median_target_heart_rank"])), lty=2, col="red")
abline(v=mean(as.numeric(rbp_info_no_blanks[rbp_info_no_blanks[,"RBP"] == "QKI","median_target_heart_rank"])), lty=3, col="red")
lines(density(rbp_fraction_hhe_targets, from=0, to=1), col="green", lwd=2)
abline(v=mean(as.numeric(rbp_info_no_blanks[rbp_info_no_blanks[,"RBP"] == "RBFOX2","fraction_hhe_targets"])), lty=2, col="green")
abline(v=mean(as.numeric(rbp_info_no_blanks[rbp_info_no_blanks[,"RBP"] == "QKI","fraction_hhe_targets"])), lty=3, col="green")
legend("topright", legend=c(paste0("all genes (N=",nrow(heart_expression_ranks_dat),")"), paste0("RBPs (N=",nrow(rbp_info_no_blanks),")"), "RBP median targets", "RBP fraction HHE targets", "density", "RBFOX2", "QKI"), col=c("black", "blue", "red", "green", "black", "black", "black"), pch=c(15, 15, 15, 15, NA, NA, NA), lty=c(NA, NA, NA, NA, 1, 2, 3))
dev.off()

rbp_subset <- read.csv(output_path("rbp_contribution_1. RBP H3K36me3 in ES-UCSF4  Cells (SNVs).csv")); rbp_subset <- rbp_subset[rbp_subset$hits_cases > 0,]
rbp_subset <- rbp_subset[1:floor(nrow(rbp_subset)/4),]
rbp_info_no_blanks <- rbp_info_no_blanks[rbp_info_no_blanks$RBP %in% paste0(rbp_subset$RBP),]
rbp_heart_ranks_result1 <- as.numeric(rbp_info_no_blanks[,"heart_rank"])
rbp_median_target_heart_ranks_result1 <- as.numeric(rbp_info_no_blanks[,"median_target_heart_rank"])

rbp_info_no_blanks <- rbp_info[rbp_info$heart_rank != ".",]
rbp_subset <- read.csv(output_path("rbp_contribution_2. constrained TSS RBP H3K36me3 in H1 Cells (SNVs).csv")); rbp_subset <- rbp_subset[rbp_subset$hits_cases > 0,]
rbp_subset <- rbp_subset[1:floor(nrow(rbp_subset)/4),]
rbp_info_no_blanks <- rbp_info_no_blanks[rbp_info_no_blanks$RBP %in% paste0(rbp_subset$RBP),]
rbp_heart_ranks_result2 <- as.numeric(rbp_info_no_blanks[,"heart_rank"])
rbp_median_target_heart_ranks_result2 <- as.numeric(rbp_info_no_blanks[,"median_target_heart_rank"])

cor(as.numeric(rbp_info[rbp_info$heart_rank != ".", "heart_rank"]), as.numeric(rbp_info[rbp_info$heart_rank != ".", "fraction_hhe_targets"]))
pdf(file=output_path(paste0("rbp_target_expression_vs_rbp_expression.pdf")))
#batch <- as.factor(sapply(ssc_parental_age$Blinded.ID, function(sample) { if (paste0(sample) %in% paste0(sscfb2$sample)) { return("SSC2") } else { return("SSC1") } }))
#scatterplot(rbp_heart_ranks, rbp_median_target_heart_ranks, legend.plot=FALSE, groups=batch, pch=c(1,1,1), boxplots="xy", main="SSCFB Samples Summary", xlab="paternal age (years)", ylab="variants/sample", cex.main=1.5, cex.axis=1.4, cex.lab=1.3, spread=FALSE, smooth=FALSE)
scatterplot(rbp_heart_ranks, rbp_median_target_heart_ranks, legend.plot=FALSE, pch=c(1,1), boxplots="xy", main="RBP target expression vs RBP expression", xlab="RBP heart rank", ylab="RBP median target heart rank", cex.main=1.5, cex.axis=1.4, cex.lab=1.3, spread=FALSE, smooth=FALSE)
#legend("topleft", legend=c("1st batch","2nd batch"), col=c("blue","magenta"), pch=15)
#mtext(paste0("Pearson = ", round(cor(rbp_heart_ranks, rbp_median_target_heart_ranks, method="pearson"),3)))
dev.off()

cor(as.numeric(rbp_info[rbp_info$heart_rank != ".", "heart_rank"]), as.numeric(rbp_info[rbp_info$heart_rank != ".", "median_target_heart_rank"]))
pdf(file=output_path(paste0("rbp_target_expression_vs_rbp_expression.pdf")))
#batch <- as.factor(sapply(ssc_parental_age$Blinded.ID, function(sample) { if (paste0(sample) %in% paste0(sscfb2$sample)) { return("SSC2") } else { return("SSC1") } }))
#scatterplot(rbp_heart_ranks, rbp_median_target_heart_ranks, legend.plot=FALSE, groups=batch, pch=c(1,1,1), boxplots="xy", main="SSCFB Samples Summary", xlab="paternal age (years)", ylab="variants/sample", cex.main=1.5, cex.axis=1.4, cex.lab=1.3, spread=FALSE, smooth=FALSE)
scatterplot(rbp_heart_ranks, rbp_fraction_hhe_targets, legend.plot=FALSE, pch=c(1,1), boxplots="xy", main="RBP target expression vs RBP expression", xlab="RBP heart rank", ylab="Fraction of RBP targets in HHE", cex.main=1.5, cex.axis=1.4, cex.lab=1.3, spread=FALSE, smooth=FALSE)
#legend("topleft", legend=c("1st batch","2nd batch"), col=c("blue","magenta"), pch=15)
#mtext(paste0("Pearson = ", round(cor(rbp_heart_ranks, rbp_median_target_heart_ranks, method="pearson"),3)))
dev.off()


exonic_variants <- read.csv(output_path("WGS_exonic_variants_2019_02_20.txt"), sep="\t")
chdfb_combined_muts_wgsa <- chdfb_combined_muts_wgsa[!(paste(chdfb_combined_muts_wgsa$Chrom_hg38, chdfb_combined_muts_wgsa$Position_hg38) %in% paste(exonic_variants$Chrom_hg38, exonic_variants$Pos_hg38)),]
sscfb_combined_muts_wgsa <- sscfb_combined_muts_wgsa[!(paste(sscfb_combined_muts_wgsa$Chrom_hg38, sscfb_combined_muts_wgsa$Position_hg38) %in% paste(exonic_variants$Chrom_hg38, exonic_variants$Pos_hg38)),]


num_best_results = 7
best_result_indices = order(as.numeric(paste0(chdfb_burdens_all_2s_full$p.value)))[1:num_best_results]
sapply_out <- sapply(1:length(best_result_indices), function(i) {
    best_result_index = best_result_indices[i]
    burden <- chdfb_burdens_all_2s_full[best_result_index,]
    name = paste0(i,". ",gsub(" \\(E.+Peak\\)", gsub("SNP", "SNV", paste0(" (",burden$variants,")")), paste0(burden$burden)))
    a <- unlist(sapply(unlist(strsplit(paste0(burden$cases_gene_hits), ",")), function(target) {
        target_split <- unlist(strsplit(target, "\\(|\\)"))
        if(length(target_split) > 1) { return(rep(target_split[1], as.numeric(target_split[2])))
        } else { return(target) }   
    })); names(a) <- NULL
    
    b <- unlist(sapply(unlist(strsplit(paste0(burden$controls_gene_hits), ",")), function(target) {
        target_split <- unlist(strsplit(target, "\\(|\\)"))
        if(length(target_split) > 1) { return(rep(target_split[1], as.numeric(target_split[2])))
        } else { return(target) }   
    })); names(b) <- NULL
    
    pdf(output_path(paste0("rbp_target_distribution_",name,".pdf")))
    barplot(table(sort(table(a))), space=0, main="Distribution of RBP target disturbance", xlab="variants per target", ylab="# targets", cex.axis=1.4, cex.names=1.4, cex.lab=1.4, cex.main=1.4) #cex.mtext=1.1
    mtext(name, cex=1.1)
    dev.off()
    
    a_genes <- sort(table(a), decreasing=TRUE)
    a_genes <- data.frame(cbind(names(a_genes), a_genes))
    rownames(a_genes) <- NULL
    colnames(a_genes) <- c("gene", "hits_cases")
    a_genes$hits_cases <- as.numeric(paste0(a_genes$hits_cases))
    
    b_genes <- sort(table(b), decreasing=TRUE)
    b_genes <- data.frame(cbind(names(b_genes), b_genes))
    rownames(b_genes) <- NULL
    colnames(b_genes) <- c("gene", "hits_controls")
    b_genes$hits_controls <- as.numeric(paste0(b_genes$hits_controls))
    
    a_genes <- cbind(a_genes, unlist(sapply(a_genes[,"gene"], function(x) { pli <- pLI_scores[[paste0(x)]]; if(is.null(pli)) { pli = "." } else { pli = round(pli, 3) }; return(pli) }))); colnames(a_genes)[ncol(a_genes)] <- "pLI"
    a_genes <- merge(a_genes, heart_expression_ranks_dat, by=("gene"), all.x=TRUE)
    a_genes[,"heart_rank"] <- paste0(a_genes[,"heart_rank"])
    a_genes[which(is.na(a_genes[,"heart_rank"]) | a_genes[,"heart_rank"] == "NA"), "heart_rank"] <- "."
    
    a_genes <- merge(a_genes, b_genes, by=("gene"), all.x=TRUE); colnames(a_genes)[ncol(a_genes)] <- "hits_controls"
    a_genes[is.na(a_genes[,"hits_controls"]), "hits_controls"] <- 0
    
    a_genes <- data.frame(a_genes)
    a_genes$heart_rank <- paste0(a_genes$heart_rank); a_genes$heart_rank[is.na(a_genes$heart_rank) | a_genes$heart_rank == "NA"] <- "."
    
    # Annotate case variants
    a_hits <- t(sapply(unlist(strsplit(paste0(burden$cases_hg38_variant_hits), ";")), function(x) {
        x_split <- unlist(strsplit(x, "_\\(|\\)"))
        return(x_split[c(2,1)])
        #if(x_split[2] == gene) { return(x_split[1]) } else { return(NULL) }
    })); rownames(a_hits) <- NULL; colnames(a_hits) <- c("gene", "case_variants")
    a_hits <- aggregate(a_hits[,"case_variants"], by=list(a_hits[,"gene"]), FUN=paste, collapse=", ")
    colnames(a_hits) <- c("gene", "case_variants")
    
    a_genes <- merge(a_genes, a_hits, by=("gene"), all.x=TRUE); colnames(a_genes)[ncol(a_genes)] <- "case_variants_hg38"
    
    a_genes <- a_genes[order(a_genes$hits_cases, decreasing=TRUE),]
    
    write.csv(a_genes, file=output_path(paste0("rbp_target_counts_",name,".csv")), row.names=FALSE)
    
    ##########################################################
    # Find contribution of individual RBPs to signal
    ##########################################################
    case_variants_hg38 <- genomic_strings_to_dat(burden[,"cases_hg38_variant_hits"], column_names=c("sample", "Chrom_hg38", "Position_hg38", "Ref", "Alt"), column_order=c(2,3,4,5,1))
    write.csv(case_variants_hg38 , file=output_path(paste0("chd_variants_",name,".csv")), row.names=FALSE)
    control_variants_hg38  <- genomic_strings_to_dat(burden[,"controls_hg38_variant_hits"], column_names=c("sample", "Chrom_hg38", "Position_hg38", "Ref", "Alt"), column_order=c(2,3,4,5,1))
    write.csv(control_variants_hg38 , file=output_path(paste0("ssc_variants_",name,".csv")), row.names=FALSE)
    
    case_variants_hg19 <- genomic_strings_to_dat(burden[,"cases_hg19_variant_hits"], column_names=c("sample", "Chrom", "Position", "Ref", "Alt"), column_order=c(2,3,4,5,1))
    control_variants_hg19  <- genomic_strings_to_dat(burden[,"controls_hg19_variant_hits"], column_names=c("sample", "Chrom", "Position", "Ref", "Alt"), column_order=c(2,3,4,5,1))
    
    case_variants_hg19 <- annotate(case_variants_hg19, get_features_by_group("RBP"))
    control_variants_hg19 <- annotate(control_variants_hg19, get_features_by_group("RBP"))
    
    a_genes <- sort(apply(case_variants_hg19[,-c(1:5)], 2, function(x) { sum(x == "Y") }), decreasing=TRUE)
    
    pdf(output_path(paste0("rbp_distribution_",name,".pdf")))
    barplot(table(a_genes), space=0, main="Distribution of RBP disturbance", xlab="variants per RBP", ylab="# RBPs", cex.axis=1.4, cex.names=1.4, cex.lab=1.4, cex.main=1.4) #cex.mtext=1.1
    mtext(name, cex=1.1)
    dev.off()
    
    a_genes_names_split <- t(data.frame(strsplit(names(a_genes),"\\.")))
    a_genes <- data.frame(cbind(a_genes_names_split[,2], a_genes_names_split[,1], a_genes))
    rownames(a_genes) <- NULL
    colnames(a_genes) <- c("gene", "cell_line", "hits_cases")
    a_genes$hits_cases <- as.numeric(paste0(a_genes$hits_cases))
    
    b_genes <- sort(apply(control_variants_hg19[,-c(1:5)], 2, function(x) { sum(x == "Y") }), decreasing=TRUE)
    b_genes_names_split <- t(data.frame(strsplit(names(b_genes),"\\.")))
    b_genes <- data.frame(cbind(b_genes_names_split[,2], b_genes_names_split[,1], b_genes))
    rownames(b_genes) <- NULL
    colnames(b_genes) <- c("gene", "cell_line", "hits_controls")
    b_genes$hits_controls <- as.numeric(paste0(b_genes$hits_controls))
    
    a_genes <- cbind(a_genes, unlist(sapply(a_genes[,"gene"], function(x) { pli <- pLI_scores[[paste0(x)]]; if(is.null(pli)) { pli = "." } else { pli = round(pli, 3) }; return(pli) }))); colnames(a_genes)[ncol(a_genes)] <- "pLI"
    a_genes <- merge(a_genes, heart_expression_ranks_dat, by=("gene"), all.x=TRUE)
    a_genes[,"heart_rank"] <- paste0(a_genes[,"heart_rank"])
    a_genes[which(is.na(a_genes[,"heart_rank"]) | a_genes[,"heart_rank"] == "NA"), "heart_rank"] <- "."
    
    a_genes <- merge(a_genes, b_genes, by=c("gene","cell_line"), all.x=TRUE); colnames(a_genes)[ncol(a_genes)] <- "hits_controls"
    a_genes[is.na(a_genes[,"hits_controls"]), "hits_controls"] <- 0
    
    a_genes <- data.frame(a_genes)
    a_genes$heart_rank <- paste0(a_genes$heart_rank); a_genes$heart_rank[is.na(a_genes$heart_rank) | a_genes$heart_rank == "NA"] <- "."
    
    a_genes <- a_genes[order(a_genes$hits_cases, decreasing=TRUE),]
    colnames(a_genes)[1] <- "RBP"
    
    a_genes <- merge(a_genes, rbp_info[,c("RBP","cell_line", "footprint", "median_target_heart_rank")], by=c("RBP","cell_line"), all.x=TRUE); colnames(a_genes)[(ncol(a_genes)-1):ncol(a_genes)] <- c("footprint", "median_target_heart_rank")
    normalized_contribution <- a_genes$hits_cases / a_genes$footprint; normalized_contribution <- normalized_contribution / sum(normalized_contribution)
    a_genes <- cbind(a_genes, normalized_contribution)
    
    fet_results <- t(sapply(1:nrow(a_genes), function(i) {
        fet_result <- fisher_exact_test(a_genes$hits_cases[i], a_genes$hits_controls[i], sum(a_genes$hits_cases), sum(a_genes$hits_controls))
        return(c(fet_result[["estimate"]], fet_result[["p.value"]]))
    })); rownames(fet_results) <- NULL; colnames(fet_results) <- c("enrichment", "p.value")
    a_genes <- cbind(a_genes, fet_results)
    
    a_genes <- a_genes[,c("RBP","cell_line", "hits_cases", "pLI", "heart_rank", "median_target_heart_rank", "hits_controls", "footprint", "normalized_contribution", "enrichment", "p.value")]
    a_genes <- a_genes[order(a_genes$normalized_contribution, a_genes$hits_cases, decreasing=TRUE),]
    write.csv(a_genes, file=output_path(paste0("rbp_contribution_",name,".csv")), row.names=FALSE)
})




cases <- standardize_colnames(chdfb_combined_muts_wgsa[,1:5])
controls <- standardize_colnames(sscfb_combined_muts_wgsa[,1:5])
controls$sample <- paste0("random_",controls$sample)
regulatory_variant_dat <- rbind(cases, controls)


library(keras)
library(kerasR)
library(stringr)
library(pbapply)
library(EBImage)
model <- load_model_hdf5("../ML/output/k562.hnrnpu_model.h5")

# Convert genomic sequences into a tensor
genomic_sequences_to_tensor <- function(sequences, sequence_annotations=NULL, sequence_length=NULL) {
    bases <- c("A","C","G","T")
    if(is.null(sequence_length)) { sequence_length = max(nchar(sequences)) }
    num_sequences = length(sequences)
    tensor <- array(0, dim=c(num_sequences, sequence_length, 4, 1), dimnames = list(paste0("sequence_",1:num_sequences), 1:sequence_length, bases, "channel")); #colnames(mat) <- apply(expand.grid(bases, 1:sequence_length)[,c(2,1)], 1, paste0, collapse="")
    sequences_split <- strsplit(sequences,"")
    for(i in 1:num_sequences) { 
        sequence <- sequences_split[[i]]
        for(base_position in 1:min(c(sequence_length, length(sequence)))) { tensor[i, base_position, which(bases == sequence[base_position]), 1] <- 1 }
    }
    
    return(tensor)
}

# Add sequence context features:
add_sequence_context_feature <- function(dat, chr_colname="Chrom", pos_colname="Position", ref_colname="Ref", alt_colname="Alt", width=51, version="hg19") {
    #sequences <- get_trimers(dat[,chr_colname], dat[,pos_colname], width=width, allow_BSgenome=TRUE, version=version)
    arm_width = floor(width/2)
    all_starts <- as.numeric(paste0(dat[,pos_colname])) - arm_width; all_ends <- as.numeric(paste0(dat[,pos_colname])) + arm_width
    all_chromosomes <- paste0(dat[,chr_colname])
    all_refs <- paste0(dat[,ref_colname])
    all_refs_lengths <- nchar(all_refs)
    all_alts <- paste0(dat[,alt_colname])
    
    sequences <- lapply(unique(all_chromosomes), function(chromosome) {
        if(grepl("Y",chromosome)) { return(matrix(nrow=0, ncol=ncol(dat)+2)) }
        cat(paste0("Reading chr", gsub("chr","",chromosome), "..."))
        refseq <- get_refseq(chromosome, version="hg19", allow_BSgenome=TRUE)[[1]]
        cat("Done.\nGrabbing sequences...")
        curr_chrom_indices <- which(all_chromosomes == chromosome)
        starts <- all_starts[curr_chrom_indices]; ends <- all_ends[curr_chrom_indices]; widths <- ends - starts + 1
        split_indices <- cumsum(widths)
        curr_chrom_rbp_sequences <- eval(parse(text=paste0("substring(paste0(refseq[c(",paste0(starts,":",ends, collapse=","),")]), c(0,split_indices[-length(split_indices)])+1, split_indices)")))
        cat("Done.\n")
        
        mutated_sequences <- paste0(substr(curr_chrom_rbp_sequences, 1, arm_width), all_alts[curr_chrom_indices],substr(curr_chrom_rbp_sequences, arm_width+1+all_refs_lengths[curr_chrom_indices], nchar(curr_chrom_rbp_sequences)))  
        curr_chrom_rbp_sequences <- substring(curr_chrom_rbp_sequences, 1, width)
        mutated_sequences <- substring(mutated_sequences, 1, width)
        #substr(mutated_sequences, arm_width+1, arm_width+1) <- all_alts[curr_chrom_indices]
        #mutated_sequences <- curr_chrom_rbp_sequences
        #substr(mutated_sequences, arm_width+1, arm_width+1) <- all_alts[curr_chrom_indices]
        return(cbind(dat[curr_chrom_indices,], c(curr_chrom_rbp_sequences), c(mutated_sequences)))
    })
    sequences <- data.frame(rbindlist(lapply(sequences, function(x) return(data.frame(x)))))
    colnames(sequences)[(ncol(sequences)-1):ncol(sequences)] <- c("ref_sequence", "alt_sequence")
    sequences$ref_sequence <- paste0(sequences$ref_sequence); sequences$alt_sequence <- paste0(sequences$alt_sequence)
    
    return(sequences)
}
regulatory_variant_dat <- add_sequence_context_feature(regulatory_variant_dat, width=151)

regulatory_variant_refs_tensor <- genomic_sequences_to_tensor(regulatory_variant_dat$ref_sequence, sequence_length=151)
regulatory_variant_alts_tensor <- genomic_sequences_to_tensor(regulatory_variant_dat$alt_sequence, sequence_length=151)
control_indices <- grepl("random", regulatory_variant_dat$sample)

prediction_score_to_likelihood_mapping_table <- read.csv(file="../ML/output/prediction_score_to_likelihood_mapping_table.csv")
rbps_to_focus_on <- get_features_by_group("RBP") #c("K562.RBFOX2", "K562.EFTUD2", "K562.HNRNPU")
#rbps_to_focus_on <- rbps_to_focus_on[gsub("^.+\\.","",rbps_to_focus_on) %in% ls(get_constrained_genes("pLI>0.5"))]
for(rbp in rbps_to_focus_on[c(1:160)]) {
    print(rbp)
    model <- load_model_hdf5(paste0("../ML/output/",tolower(rbp),"_model2.h5"))
    ref_pred_scores <- model %>% predict(regulatory_variant_refs_tensor[,,,])
    alt_pred_scores <- model %>% predict(regulatory_variant_alts_tensor[,,,])
    ref_LR <- prediction_score_to_likelihood_mapping_table[round_to_nearest(ref_pred_scores/0.01)+1,rbp]
    alt_LR <- prediction_score_to_likelihood_mapping_table[round_to_nearest(alt_pred_scores/0.01)+1,rbp]
    delta_pred_scores <- alt_LR - ref_LR
    #delta_pred_scores <- -log(alt_pred_scores/ref_pred_scores)
    #delta_pred_scores[(alt_pred_scores < ref_pred_scores & (ref_pred_scores < 0.5 | alt_pred_scores > 0.5)) | (alt_pred_scores > ref_pred_scores & (ref_pred_scores > -0.5 | alt_pred_scores < -0.5))] <- 0
    return_dat <- cbind(ref_pred_scores, alt_pred_scores, ref_LR, alt_LR, delta_pred_scores)
    colnames(return_dat) <- c("ref_pred_score", "alt_pred_score", "ref_LR", "alt_LR", "delta_LR")
    write.csv(return_dat, file=output_path(paste0(rbp,"_CHD_scores.csv")), row.names=FALSE)
    #return(return_dat)
    rm(model); gc()
}
#saveRDS(variant_pred_scores, file = output_path("variant_pred_scores.rds"))
delta_pred_scores <- unfactorize(data.frame(sapply(rbps_to_focus_on, function(rbp) {
    print(rbp)
    return(read.csv(output_path(paste0(rbp,"_CHD_scores.csv")))[,5])
})))

delta_pred_scores <- cbind(delta_pred_scores[,1:160], t(apply(delta_pred_scores[,gsub("^.+\\.","",colnames(delta_pred_scores)) %in% ls(get_constrained_genes("pLI>0.5"))], 1, FUN=range)))
#delta_pred_scores <- cbind(delta_pred_scores[,1:160], t(apply(delta_pred_scores, 1, FUN=range)))
colnames(delta_pred_scores)[(ncol(delta_pred_scores)-1):ncol(delta_pred_scores)] <- c("GOF", "LOF")

saveRDS(delta_pred_scores, file=output_path("delta_pred_scores_CHD.rds"))

delta_pred_scores <- readRDS(file="../WGS/output/delta_pred_scores_CHD.rds")


#delta_pred_scores <- apply(delta_pred_scores, 2, function(x) (x - mean(x)) / sd(x))
sapply_out <- sapply(1:ncol(delta_pred_scores), function(i) print(range(delta_pred_scores[,i])))
#filename = output_path(paste0("delta_prediction_score_densities.pdf"))
#pdf(file=filename)
cols <- rainbow(ncol(delta_pred_scores))
#plot_start = -0.9; plot_end = -0.5; y_max = 0.05
#plot(density(delta_pred_scores[,1], from=plot_start), col="white", main="Delta prediction score densities", lwd=2, lty=1, xlab="-log(alt_pred_scores/ref_pred_scores)", xlim=c(plot_start, plot_end), ylim=c(0, y_max), xaxs="i", cex.axis=1.4, cex.lab=1.4, cex.main=1.4) # xlim
control_indices <- grepl("random", regulatory_variant_dat$sample)
sapply_out <- sapply(1:ncol(delta_pred_scores), function(i) print(range(delta_pred_scores[,i])))
sapply_out <- sapply(1:ncol(delta_pred_scores), function(i) print(range(delta_pred_scores[!control_indices,i])))
sapply_out <- sapply(1:ncol(delta_pred_scores), function(i) print(range(delta_pred_scores[control_indices,i])))
starts <- rep(0.25, ncol(delta_pred_scores)) #c(0.25, 0.25, 0.25)
bandwidths = rep(0.01, ncol(delta_pred_scores)) #c(0.01, 0.01, 0.002)
regulatory_variant_dat_indels <- nchar(paste0(regulatory_variant_dat$Ref)) > 1 | nchar(paste0(regulatory_variant_dat$Alt)) > 1 
all_limit_variants = c("indels_only", "SNVs_only", "")
for(limit_variants in all_limit_variants) {
    sapply_out <- sapply(1:ncol(delta_pred_scores), function(i) {
        rbp = colnames(delta_pred_scores)[i]
        print(rbp)
        if(limit_variants == "SNVs_only") {
            delta_pred_scores_cases <- delta_pred_scores[(!control_indices) & (!regulatory_variant_dat_indels),i]
            delta_pred_scores_controls <- delta_pred_scores[control_indices & (!regulatory_variant_dat_indels),i]
            rbp = paste0(rbp," SNVs")
        } else if(limit_variants == "indels_only") {
            delta_pred_scores_cases <- delta_pred_scores[(!control_indices) & (regulatory_variant_dat_indels),i]
            delta_pred_scores_controls <- delta_pred_scores[control_indices & (regulatory_variant_dat_indels),i]
            rbp = paste0(rbp," indels")
        } else {
            delta_pred_scores_cases <- delta_pred_scores[!control_indices,i]
            delta_pred_scores_controls <- delta_pred_scores[control_indices,i]
        }
        tail_width = 0.5
        
        #filename = output_path(paste0(gsub(" ","_",rbp),"_delta_prediction_score_densities.pdf"))
        #pdf(file=filename, width=14)
        #par(mfrow=c(1,2))
        
        delta_pred_scores_cases_density <- density(delta_pred_scores_cases, bw=bandwidths[i])
        delta_pred_scores_controls_density <- density(delta_pred_scores_controls, bw=bandwidths[i])
        
        #plot_start = min(delta_pred_scores[,i]); plot_end = -starts[i] #plot_start * (1-tail_width)
        ##delta_pred_scores_cases_density <- density(delta_pred_scores_cases, to=plot_end, bw=bandwidth)
        ##delta_pred_scores_controls_density <- density(delta_pred_scores_controls, to=plot_end, bw=bandwidth)
        ##print("Cases density: "); print(delta_pred_scores_cases_density)
        ##print("Controls density: "); print(delta_pred_scores_controls_density)
        #plot(delta_pred_scores_controls_density$x[delta_pred_scores_controls_density$x <= plot_end], delta_pred_scores_controls_density$y[delta_pred_scores_controls_density$x <= plot_end], col="blue", main=paste0(rbp," delta prediction score densities"), lwd=2, lty=1, xlab="-log(alt_pred_scores/ref_pred_scores)", ylab="Density", xlim=c(1.1*plot_start, plot_end), xaxs="i", cex.axis=1.4, cex.lab=1.4, cex.main=1.4, type="l")    
        #lines(delta_pred_scores_cases_density$x[delta_pred_scores_cases_density$x <= plot_end], delta_pred_scores_cases_density$y[delta_pred_scores_cases_density$x <= plot_end], col="red", lwd=2, lty=1)
        #legend("topleft", legend=c("CHD", "SSC"), col=c("red", "blue"), pch=15, cex=1.3)
        #mtext(paste0("score deltas below threshold ",plot_end), cex=1.2) #mtext(paste0("score deltas (lowest ",tail_width*100,"% tail)"), cex=1.2)
        #
        #plot_end = max(delta_pred_scores[,i]); plot_start = starts[i] #plot_start = plot_end * (1-tail_width)
        ##delta_pred_scores_cases_density <- density(delta_pred_scores_cases, from=plot_start, bw=bandwidth)
        ##delta_pred_scores_controls_density <- density(delta_pred_scores_controls, from=plot_start, bw=bandwidth)
        ##print("Cases density: "); print(delta_pred_scores_cases_density)
        ##print("Controls density: "); print(delta_pred_scores_controls_density)
        #plot(delta_pred_scores_controls_density$x[delta_pred_scores_controls_density$x >= plot_start], delta_pred_scores_controls_density$y[delta_pred_scores_controls_density$x >= plot_start], col="blue", main=paste0(rbp," delta prediction score densities"), lwd=2, lty=1, xlab="-log(alt_pred_scores/ref_pred_scores)", ylab="Density", xlim=c(plot_start, 1.1*plot_end), xaxs="i", cex.axis=1.4, cex.lab=1.4, cex.main=1.4, type="l")
        #lines(delta_pred_scores_cases_density$x[delta_pred_scores_cases_density$x >= plot_start], delta_pred_scores_cases_density$y[delta_pred_scores_cases_density$x >= plot_start], col="red", lwd=2, lty=1)
        #legend("topright", legend=c("CHD", "SSC"), col=c("red", "blue"), pch=15, cex=1.3)
        #mtext(paste0("score deltas above threshold ",plot_start), cex=1.2) #mtext(paste0("score deltas (highest ",tail_width*100,"% tail)"), cex=1.2)
        #dev.off()
        #pdf_to_png(filename)
        #par(mfrow=c(1,1)) 
        
        thresholds <- seq(-1.5, 1.5, by=0.1)
        greater_than_threshold_fets <- unfactorize(data.frame(t(sapply(thresholds, function(threshold) { 
            fet_result <- fisher.test(matrix(c(sum(delta_pred_scores_cases > threshold), sum(delta_pred_scores_controls > threshold), sum(delta_pred_scores_cases <= threshold), sum(delta_pred_scores_controls <= threshold)), nrow = 2, dimnames = list(hits = c("Y", "N"), status = c("case", "control"))), alternative = "two.sided")
            return(c(threshold, fet_result$estimate, fet_result$conf.int, fet_result$p.value)) 
        })))); colnames(greater_than_threshold_fets) <- c("threshold", "estimate", "conf.int_lower", "conf.int_higher", "p.value")
        greater_than_threshold_fets[greater_than_threshold_fets == Inf] <- 0
        greater_than_threshold_fets[,c("estimate", "conf.int_lower", "conf.int_higher")] <- log2(greater_than_threshold_fets[,c("estimate", "conf.int_lower", "conf.int_higher")])
        greater_than_threshold_fets[greater_than_threshold_fets == -Inf] <- 0
        less_than_threshold_fets <- unfactorize(data.frame(t(sapply(thresholds, function(threshold) { 
            fet_result <- fisher.test(matrix(c(sum(delta_pred_scores_cases < threshold), sum(delta_pred_scores_controls < threshold), sum(delta_pred_scores_cases >= threshold), sum(delta_pred_scores_controls >= threshold)), nrow = 2, dimnames = list(hits = c("Y", "N"), status = c("case", "control"))), alternative = "two.sided")
            return(c(threshold, fet_result$estimate, fet_result$conf.int, fet_result$p.value)) 
        })))); colnames(less_than_threshold_fets) <- c("threshold", "estimate", "conf.int_lower", "conf.int_higher", "p.value")
        less_than_threshold_fets[less_than_threshold_fets == Inf] <- 0
        less_than_threshold_fets[,c("estimate", "conf.int_lower", "conf.int_higher")] <- log2(less_than_threshold_fets[,c("estimate", "conf.int_lower", "conf.int_higher")])
        less_than_threshold_fets[less_than_threshold_fets == -Inf] <- 0
        
        filename = output_path(paste0(gsub(" ","_",rbp),"_delta_prediction_score_enrichments.pdf"))
        pdf(file=filename, width=14)
        par(mfrow=c(1,2)) 
        
        cols <- c("black","red")[as.numeric(greater_than_threshold_fets$p.value < 0.05)+1]
        plot(greater_than_threshold_fets$threshold, greater_than_threshold_fets$estimate, main=paste0("CHD ",rbp," LOF Enrichment"), xlab="RBP binding delta threshold", ylab="log2(Odds Ratio)", col=cols, pch=19, ylim=c(min(greater_than_threshold_fets$conf.int_lower), max(greater_than_threshold_fets$conf.int_higher)), cex.axis=1.4, cex.lab=1.4, cex.main=1.4)
        abline(h=0, col="blue")
        mtext("Binding delta > threshold", cex=1.2)
        segments(x0=greater_than_threshold_fets$threshold, y0=greater_than_threshold_fets$conf.int_lower, x1=greater_than_threshold_fets$threshold, y1=greater_than_threshold_fets$conf.int_higher, col=cols)
        legend("topleft", legend=c("NS", "p<0.5"), col=c("black", "red"), pch=19, cex=1.2)
        
        cols <- c("black","red")[as.numeric(less_than_threshold_fets$p.value < 0.05)+1]
        plot(less_than_threshold_fets$threshold, less_than_threshold_fets$estimate, main=paste0("CHD ",rbp," GOF Enrichment"), xlab="RBP binding delta threshold", ylab="log2(Odds Ratio)", col=cols, pch=19, ylim=c(min(less_than_threshold_fets$conf.int_lower), max(less_than_threshold_fets$conf.int_higher)), cex.axis=1.4, cex.lab=1.4, cex.main=1.4)
        abline(h=0, col="blue")
        mtext("Binding delta < threshold", cex=1.2)
        segments(x0=less_than_threshold_fets$threshold, y0=less_than_threshold_fets$conf.int_lower, x1=less_than_threshold_fets$threshold, y1=less_than_threshold_fets$conf.int_higher, col=cols)
        legend("topright", legend=c("NS", "p<0.5"), col=c("black", "red"), pch=19, cex=1.2)
        
        dev.off()
        pdf_to_png(filename)
        par(mfrow=c(1,1))
    })
}

variant_annotations <- annotate(regulatory_variant_dat, c("autism_genes_20000bp", "constrained_genes_20000bp", "H3K36me3"))[,-c(1:ncol(regulatory_variant_dat))] == "Y"
autism_gene_indices <- variant_annotations[,1]; constrained_gene_indices <- variant_annotations[,2]; H3K36me3_indices <- constrained_gene_indices <- variant_annotations[,3]

# Run variant threshold test to get real p.values
variant_types <- c("SNV", "indel")
variant_constraints <- c("autism_gene", "constrained_gene", "H3K36me3", "")
thresholds <- seq(4, 150, by=2)
#thresholds <- seq(0.5, 1, by=0.01)
do_gof = TRUE
sapply_out <- sapply(161:161, function(i) { #1:ncol(delta_pred_scores)
    variable_threshold_results_lof <- data.frame()
    variable_threshold_results_gof <- data.frame()
    for(variant_constraint in variant_constraints) {
        rbp = paste0("",colnames(delta_pred_scores)[i])
        if(variant_constraint == "") { constrained_variant_indices <- rep(TRUE, nrow(delta_pred_scores))
        } else { constrained_variant_indices <- get(paste0(variant_constraint,"_indices")); rbp = paste0(rbp," ",gsub("_gene","",variant_constraint)) }
        
        for(variant_type in variant_types) {
            # print(paste0(c(rbp,variant_type,variant_constraint),collapse=", "))
            if(variant_type == "SNV") {
                delta_pred_scores_cases <- delta_pred_scores[(!control_indices) & (!regulatory_variant_dat_indels) & constrained_variant_indices,i]
                delta_pred_scores_controls <- delta_pred_scores[control_indices & (!regulatory_variant_dat_indels) & constrained_variant_indices,i]
            } else if(variant_type == "indel") {
                delta_pred_scores_cases <- delta_pred_scores[(!control_indices) & (regulatory_variant_dat_indels) & constrained_variant_indices,i]
                delta_pred_scores_controls <- delta_pred_scores[control_indices & (regulatory_variant_dat_indels) & constrained_variant_indices,i]
            } else {
                delta_pred_scores_cases <- delta_pred_scores[!control_indices,i]
                delta_pred_scores_controls <- delta_pred_scores[control_indices,i]
            }
            n1 = length(delta_pred_scores_cases); n0 = length(delta_pred_scores_controls)
            
            num_permutations = 100
            variable_threshold_results <- t(sapply(1:(num_permutations+1), function(j) {
                print(j)
                
                greater_than_threshold_fets <- unfactorize(data.frame(t(sapply(thresholds, function(threshold) {
                    m1 = sum(delta_pred_scores_cases > threshold); m0 = sum(delta_pred_scores_controls > threshold)
                    fet_result <- fisher.test(matrix(c(m1, n1-m1, m0, n0-m0), nrow = 2, dimnames = list(hits = c("Y", "N"), status = c("case", "control"))), alternative = "two.sided")
                    return(c(threshold, fet_result$estimate, fet_result$conf.int, fet_result$p.value, m1, m0, n1, n0)) 
                })))); colnames(greater_than_threshold_fets) <- c("threshold", "estimate", "conf.int_lower", "conf.int_higher", "p.value", "m1", "m0", "n1", "n0")
                greater_than_threshold_fets[greater_than_threshold_fets == Inf] <- 0
                #greater_than_threshold_fets[,c("estimate", "conf.int_lower", "conf.int_higher")] <- log2(greater_than_threshold_fets[,c("estimate", "conf.int_lower", "conf.int_higher")])
                #greater_than_threshold_fets[greater_than_threshold_fets == -Inf] <- 0
                
                if(do_gof) {
                    thresholds <- -thresholds
                    less_than_threshold_fets <- unfactorize(data.frame(t(sapply(thresholds, function(threshold) { 
                        m1 = sum(delta_pred_scores_cases < threshold); m0 = sum(delta_pred_scores_controls < threshold)
                        fet_result <- fisher.test(matrix(c(m1, n1-m1, m0, n0-m0), nrow = 2, dimnames = list(hits = c("Y", "N"), status = c("case", "control"))), alternative = "two.sided")
                        return(c(threshold, fet_result$estimate, fet_result$conf.int, fet_result$p.value, m1, m0, n1, n0)) 
                    })))); colnames(less_than_threshold_fets) <- c("threshold", "estimate", "conf.int_lower", "conf.int_higher", "p.value", "m1", "m0", "n1", "n0")
                    less_than_threshold_fets[less_than_threshold_fets == Inf] <- 0
                    #less_than_threshold_fets[,c("estimate", "conf.int_lower", "conf.int_higher")] <- log2(less_than_threshold_fets[,c("estimate", "conf.int_lower", "conf.int_higher")])
                    #less_than_threshold_fets[less_than_threshold_fets == -Inf] <- 0
                }
                
                #more_extreme_than_threshold_fets <- rbind(less_than_threshold_fets, greater_than_threshold_fets)
                #write.csv(more_extreme_than_threshold_fets, file=output_path(paste0(gsub(" ","_",rbp),"_",variant_type,"_delta_prediction_score_variable_threshold_test.csv")), row.names=FALSE)
                
                # shuffle for next permutation
                delta_pred_scores_all <- c(delta_pred_scores_cases, delta_pred_scores_controls)
                delta_pred_scores_cases_indices <- sample(1:length(delta_pred_scores_all), length(delta_pred_scores_cases), replace=FALSE)
                delta_pred_scores_cases <<- delta_pred_scores_all[delta_pred_scores_cases_indices]
                delta_pred_scores_controls <<- delta_pred_scores_all[-c(delta_pred_scores_cases_indices)]
                
                if(do_gof) { return(unlist(c(greater_than_threshold_fets[which.min(greater_than_threshold_fets$p.value),], less_than_threshold_fets[which.min(less_than_threshold_fets$p.value),]))) 
                } else { return(unlist(greater_than_threshold_fets[which.min(greater_than_threshold_fets$p.value),])) }
            }))
            
            variable_threshold_results_lof_entry <- c(colnames(delta_pred_scores)[i], variant_type, variant_constraint, variable_threshold_results[1,1:9], (sum(variable_threshold_results[-1,5] <= variable_threshold_results[1,5])+1)/(num_permutations+1))
            names(variable_threshold_results_lof_entry)[c(1:3,13)] <- c("RBP", "variant_type", "variant_constraint", "real_p.value")
            variable_threshold_results_lof <- rbind(unfactorize(variable_threshold_results_lof), unfactorize(variable_threshold_results_lof_entry))
            
            if(do_gof) { 
                variable_threshold_results_gof_entry <- c(colnames(delta_pred_scores)[i], variant_type, variant_constraint, variable_threshold_results[1,10:18], (sum(variable_threshold_results[-1,14] <= variable_threshold_results[1,14])+1)/(num_permutations+1))
                names(variable_threshold_results_gof_entry)[c(1:3,13)] <- c("RBP", "variant_type", "variant_constraint", "real_p.value")
                variable_threshold_results_gof <- rbind(unfactorize(variable_threshold_results_gof), unfactorize(variable_threshold_results_gof_entry))
            }
        }
    }
    colnames(variable_threshold_results_lof) <- c("RBP", "variant_type", "variant_constraint", "threshold", "estimate", "conf.int_lower", "conf.int_higher", "p.value", "m1", "m0", "n1", "n0", "real_p.value")
    write.csv(variable_threshold_results_lof, file=output_path(paste0(gsub(" ","_",rbp),"_lof_prediction_score_variable_threshold_enrichment.csv")), row.names=FALSE)
    if(do_gof) {
        colnames(variable_threshold_results_gof) <- c("RBP", "variant_type", "variant_constraint", "threshold", "estimate", "conf.int_lower", "conf.int_higher", "p.value", "m1", "m0", "n1", "n0", "real_p.value")
        write.csv(variable_threshold_results_gof, file=output_path(paste0(gsub(" ","_",rbp),"_gof_prediction_score_variable_threshold_enrichment.csv")), row.names=FALSE)
    }
})

legend("topright", legend=c(gsub("^.*\\.", "", colnames(delta_pred_scores)), "Positives", "Negatives"), col=c(cols,"black","black"), pch=c(rep(15,ncol(delta_pred_scores)),NA,NA), lty=c(rep(NA,ncol(delta_pred_scores)),1,3))
#mtext(paste0("Validation set of ",length(test_indices)," length-",sequence_length," sequences"))
dev.off()
pdf_to_png(filename)

#mean(ref_pred_scores)
#mean(alt_pred_scores)
#pdf(output_path("prediction_score_densities.pdf"))
#plot(density(ref_pred_scores), col="blue", main="Prediction score densities")
#lines(density(alt_pred_scores), col="red")
#legend("topleft", legend=c("Refs", "Alts"), col=c("blue", "red"), pch=15)
#dev.off()

#delta_pred_scores <- -log(alt_pred_scores/ref_pred_scores)
#plot(density(delta_pred_scores))

control_indices <- grepl("random", regulatory_variant_dat$sample)
delta_pred_scores_cases <- delta_pred_scores[!control_indices]
delta_pred_scores_controls <- delta_pred_scores[control_indices]
pdf(output_path(paste0(rbp,"_delta_prediction_score_densities.pdf")))
plot(density(delta_pred_scores_controls, from=0.5), col="blue", main="Delta prediction score densities", lwd=2, lty=1, xlab="-log(alt_pred_scores/ref_pred_scores)", xlim=c(1.15, 1.275), xaxs="i", cex.axis=1.4, cex.lab=1.4, cex.main=1.4) # xlim
lines(density(delta_pred_scores_cases, from=0.5), col="red", lwd=2, lty=1)
#abline(v=1.0015, lty=3)
#abline(v=1.00325, lty=3)
legend("topright", legend=c("CDH+CHD", "SSC"), col=c("red", "blue"), pch=15)
mtext("score deltas above threshold 1.15")
dev.off()


density(delta_pred_scores_cases)
density(delta_pred_scores_controls)
pdf(output_path("delta_prediction_score_densities.pdf"))
plot(density(delta_pred_scores_controls), col="blue", main="Delta prediction score densities", lwd=2, lty=1, xlab="-log(alt_pred_scores/ref_pred_scores)", xaxs="i") #, xlim=c(0.995, 1.005)) #, xlim=c(0.995, 1.005)) # xlim
lines(density(delta_pred_scores_cases), col="red", lwd=2, lty=1)
abline(v=1.0015, lty=3)
abline(v=1.00325, lty=3)
legend("topright", legend=c("CDH+CHD", "SSC"), col=c("red", "blue"), pch=15)
dev.off()

pdf(output_path("delta_prediction_score_densities.pdf"))
plot(density(delta_pred_scores_controls, from=1.15), col="blue", main="Delta prediction score densities", lwd=2, lty=1, xlab="-log(alt_pred_scores/ref_pred_scores)", xlim=c(1.15, 1.275), xaxs="i", cex.axis=1.2, cex.lab=1.2) # xlim
lines(density(delta_pred_scores_cases, from=1.15), col="red", lwd=2, lty=1)
#abline(v=1.0015, lty=3)
#abline(v=1.00325, lty=3)
legend("topright", legend=c("CDH+CHD", "SSC"), col=c("red", "blue"), pch=15)
mtext("score deltas above threshold 1.15")
dev.off()

plot(density(delta_pred_scores_controls, to=-1.15), col="blue", main="Delta prediction score densities", lwd=2, lty=1, xlab="-log(alt_pred_scores/ref_pred_scores)", xlim=c(-1.15, -1.275), xaxs="i", cex.axis=1.2, cex.lab=1.2) # xlim
lines(density(delta_pred_scores_cases, to=-1.15), col="red", lwd=2, lty=1)
#abline(v=1.0015, lty=3)
#abline(v=1.00325, lty=3)
legend("topright", legend=c("CDH+CHD", "SSC"), col=c("red", "blue"), pch=15)
mtext("score deltas below threshold -1.15")

# Load ASD
#asd_processed <- process_wgs_dataset(input_file=data_path("WGS/An_2018_Table_S2_denovos.csv"), output_file=data_path("WGS/ASD/ASD_final.tsv"), dat_sep="\t", skip=1, header=TRUE, from="hg38", to="hg19", confirm_refseq=FALSE)
asd_processed <- process_wgs_dataset(input_file=data_path("ASD_OlgaT.tsv"), output_file=data_path("WGS/ASD/ASD_final.tsv"), confirm_refseq=FALSE)
asd_variant_dat <- asd_processed$dat
asd_variant_granges <- asd_processed$granges
asd_control_indices <- !asd_variant_dat$Proband #asd_variant_dat$Pheno == "control"
asd_ssc_variant_dat <- asd_variant_dat[asd_control_indices,]
asd_ssc_variant_granges <- asd_variant_granges[asd_control_indices]
asd_variant_dat <- asd_variant_dat[!asd_control_indices,]
asd_variant_granges <- asd_variant_granges[!asd_control_indices]
asd_muts_wgsa <- asd_variant_dat
asd_ssc_muts_wgsa <- asd_ssc_variant_dat

setup_features_env()
# Burden analysis on fullPeak ASD
asd_burdens <- burden_analysis2("ASD", ssc="ASD_SSC", desired_tests=desired_tests, include_individual_RBPs=FALSE, statistical_test="fisher_exact_2s")
write.csv(asd_burdens[,c(1:10)], file=output_path("asd_burdens.csv"), row.names=FALSE)
best_result_asd <- volcano_plot2("ASD", asd_burdens[-c(1:3, which(grepl("TSS40", "overall", asd_burdens$burden))),], bonferroni_p.value_cutoff=0.05/133, highlight_best=FALSE, num_best_hits=20, pval_relative_importance=1000, number_labels=TRUE, x_min=-0.5, x_max=0.5, label_gap=0.5, label_points_max_iterations=20)



cdh1_burdens <- burden_analysis2("CDH1", desired_tests=desired_tests, include_individual_RBPs=FALSE, statistical_test="fisher_exact_2s")
write.csv(burden_analysis2("CDH1", desired_tests=NULL, include_individual_RBPs=FALSE, statistical_test="fisher_exact_2s")[,c(1:10)], file=output_path("cdh1_burdens_all.csv"), row.names=FALSE)
write.csv(burden_analysis2("CDH1", desired_tests=NULL, include_individual_RBPs=FALSE, statistical_test="binomial")[,c(1:10)], file=output_path("cdh1_burdens_binom_all.csv"), row.names=FALSE)
write.csv(cdh1_burdens[,c(1:10)], file=output_path("cdh1_burdens.csv"), row.names=FALSE)
best_result_cdh1 <- volcano_plot2("CDH1", cdh1_burdens[-c(1:3, which(grepl("TSS40", "overall", cdh1_burdens$burden))),], bonferroni_p.value_cutoff=0.05/133, highlight_best=FALSE, num_best_hits=20, pval_relative_importance=1000, number_labels=TRUE, x_min=-0.5, x_max=0.5, label_gap=0.5, label_points_max_iterations=20)

cdh2_burdens <- burden_analysis2("CDH2", desired_tests=desired_tests, include_individual_RBPs=FALSE, statistical_test="fisher_exact_2s")
write.csv(burden_analysis2("CDH2", desired_tests=NULL, include_individual_RBPs=FALSE, statistical_test="fisher_exact_2s")[,c(1:10)], file=output_path("cdh2_burdens_all.csv"), row.names=FALSE)
write.csv(burden_analysis2("CDH2", desired_tests=NULL, include_individual_RBPs=FALSE, statistical_test="binomial")[,c(1:10)], file=output_path("cdh2_burdens_binom_all.csv"), row.names=FALSE)
write.csv(cdh2_burdens[,c(1:10)], file=output_path("cdh2_burdens.csv"), row.names=FALSE)
best_result_cdh2 <- volcano_plot2("CDH2", cdh2_burdens[-c(1:3, which(grepl("TSS40", "overall", cdh2_burdens$burden))),], bonferroni_p.value_cutoff=0.05/133, highlight_best=FALSE, num_best_hits=20, pval_relative_importance=1000, number_labels=TRUE, x_min=-0.5, x_max=0.5, label_gap=0.5, label_points_max_iterations=20)

cdh3_burdens <- burden_analysis2("CDH3", desired_tests=desired_tests, include_individual_RBPs=FALSE, statistical_test="fisher_exact_2s")
write.csv(burden_analysis2("CDH3", desired_tests=NULL, include_individual_RBPs=FALSE, statistical_test="fisher_exact_2s")[,c(1:10)], file=output_path("cdh3_burdens_all.csv"), row.names=FALSE)
write.csv(burden_analysis2("CDH3", desired_tests=NULL, include_individual_RBPs=FALSE, statistical_test="binomial")[,c(1:10)], file=output_path("cdh3_burdens_binom_all.csv"), row.names=FALSE)
write.csv(cdh3_burdens[,c(1:10)], file=output_path("cdh3_burdens.csv"), row.names=FALSE)
best_result_cdh3 <- volcano_plot2("CDH3", cdh3_burdens[-c(1:3, which(grepl("TSS40", "overall", cdh3_burdens$burden))),], bonferroni_p.value_cutoff=0.05/133, highlight_best=FALSE, num_best_hits=20, pval_relative_importance=1000, number_labels=TRUE, x_min=-0.5, x_max=0.5, label_gap=0.5, label_points_max_iterations=20)

cdh_burdens <- burden_analysis2("CDH", desired_tests=desired_tests, include_individual_RBPs=FALSE, statistical_test="fisher_exact_2s")
write.csv(burden_analysis2("CDH", desired_tests=NULL, include_individual_RBPs=FALSE, statistical_test="fisher_exact_2s")[,c(1:10)], file=output_path("cdh_burdens_all.csv"), row.names=FALSE)
write.csv(burden_analysis2("CDH", desired_tests=NULL, include_individual_RBPs=FALSE, statistical_test="binomial")[,c(1:10)], file=output_path("cdh_burdens_binom_all.csv"), row.names=FALSE)
write.csv(cdh_burdens[,c(1:10)], file=output_path("cdh_burdens.csv"), row.names=FALSE)
best_result_cdh <- volcano_plot2("CDH", cdh_burdens[-c(1:3, which(grepl("TSS40", "overall", cdh_burdens$burden))),], bonferroni_p.value_cutoff=0.05/133, highlight_best=FALSE, num_best_hits=20, pval_relative_importance=1000, number_labels=TRUE, x_min=-0.5, x_max=0.5, label_gap=0.5, label_points_max_iterations=20)

chd3_burdens <- burden_analysis2("CHD3", desired_tests=desired_tests, include_individual_RBPs=FALSE, statistical_test="fisher_exact_2s")
write.csv(burden_analysis2("CHD3", desired_tests=NULL, include_individual_RBPs=FALSE, statistical_test="fisher_exact_2s")[,c(1:10)], file=output_path("chd3_burdens_all.csv"), row.names=FALSE)
write.csv(burden_analysis2("CHD3", desired_tests=NULL, include_individual_RBPs=FALSE, statistical_test="binomial")[,c(1:10)], file=output_path("chd3_burdens_binom_all.csv"), row.names=FALSE)
write.csv(chd3_burdens[,c(1:10)], file=output_path("chd3_burdens.csv"), row.names=FALSE)
best_result_chd3 <- volcano_plot2("TOPMed", chd3_burdens[-c(1:3, which(grepl("TSS40", "overall", chd3_burdens$burden))),], bonferroni_p.value_cutoff=0.05/133, highlight_best=FALSE, num_best_hits=20, pval_relative_importance=1000, number_labels=TRUE, x_min=-0.5, x_max=0.5, label_gap=0.5, label_points_max_iterations=20)

burden_analysis2 <- function(disease, ssc_controls=NULL, desired_tests=NULL, annotate_phenotypes=FALSE, statistical_test="fisher_exact", variant_type_testing=TRUE, peak_types="broadPeak", include_individual_RBPs=FALSE, include_hg38_variant_annotation=TRUE, subsample_snvs=FALSE) {
    if(grepl("SSC", disease)) { control_name <- "SSC"; expression_string_token <- "DE"; expression_tissue_name = "E11.5 mouse diaphragm"
    } else if(grepl("CDH", disease)) { control_name <- "SSC"; expression_string_token <- "DE"; expression_tissue_name = "E11.5 mouse diaphragm"
    } else if(grepl("ASD", disease)) { control_name <- "ASD_SSC"; expression_string_token <- "HE"; expression_tissue_name = "E14.5 mouse heart"  
    } else if(grepl("CHD|TOPMed", disease)) { control_name <- "SSC"; expression_string_token <- "HE"; expression_tissue_name = "E14.5 mouse heart" } else { return() }
    high_expressed_name <- paste0(paste0("H", expression_string_token)); high_expressed_genes <- get(paste0(high_expressed_name, "_genes"))
    high_medium_high_expressed_name <- paste0(paste0("HMH", expression_string_token)); high_medium_high_expressed_genes <- get(paste0(high_medium_high_expressed_name, "_genes"))
    medium_high_expressed_name <- paste0(paste0("MH", expression_string_token)); medium_high_expressed_genes <- get(paste0(medium_high_expressed_name, "_genes"))
    low_expressed_name <- paste0(paste0("L", expression_string_token)); low_expressed_genes <- get(paste0(low_expressed_name, "_genes"))
    medium_low_expressed_name <- paste0(paste0("ML", expression_string_token)); medium_low_expressed_genes <- get(paste0(medium_low_expressed_name, "_genes"))
    if (!is.null(ssc_controls)) { control_name <- ssc_controls }
    
    # Get case variables for the given disease
    cases <- standardize_colnames(get(paste0(tolower(disease)))); cases <- cases[,c(1:5,which(colnames(cases) %in% c("CADDgt10_Phred", "Chrom_hg19", "Position_hg19", "Chrom_hg38", "Position_hg38")))]
    #case_granges <- standardize_colnames(get(paste0(tolower(disease),"_granges")))
    if (statistical_test != "binomial") { n1 = nrow(cases)
    } else { n1 = get(paste0(toupper(disease), "_SAMPLE_COUNT")) } # statistical_test="binomial"
    # Get control variables for the given disease
    controls <- standardize_colnames(get(paste0(tolower(control_name)))); controls <- controls[,c(1:5,which(colnames(controls) %in% c("CADDgt10_Phred", "Chrom_hg19", "Position_hg19", "Chrom_hg38", "Position_hg38")))]
    #control_granges <- standardize_colnames(get(paste0(tolower(control_name),"_granges")))
    if (statistical_test != "binomial") { n0 = nrow(controls)
    } else { n0 = get(paste0(toupper(control_name), "_SAMPLE_COUNT")) } # statistical_test="binomial"
    
    shared_colnames <- intersect(colnames(cases), colnames(controls))
    cases <- cases[,shared_colnames]
    controls <- controls[,shared_colnames]
    
    ####################################################################################################################
    # FIGURE OUT DESIRED TESTS!
    ####################################################################################################################
    if(is.null(desired_tests)) {
        desired_variant_type_features <- get_features_by_group("variant_type") # is_snv, is_indel, and is_variant
        desired_region_type_features <- c("", "3'UTR", "TSS") # 3'UTR, TSS, and TSS40
        desired_gene_set_features <- c("", "constrained") #c("", "HE", "constrained", "constrained HE") #get_features_by_group("gene_set")) # HDE (25% high diaphragm expressed), HHE (25% high heart expressed), HE (50% heart expressed)
        desired_RBP_features <- c("", "RBP") # combined RBP peaks      # plus each individual RBP feature, get_features_by_group("RBP")
        if(include_individual_RBPs) { desired_RBP_features <- c(desired_RBP_features, get_features_by_group("RBP")) }
        #desired_histone_mark_features <- c("", get_features_by_group("H3K79me2"), get_features_by_group("H3K36me3"), get_features_by_group("H3K4me1"), get_features_by_group("H3K27ac"), get_features_by_group("H3K9me3")) # H3K36me3 in all relevant tissue types
        desired_histone_mark_features <- c("", get_features_by_group("H3K79me2"), get_features_by_group("H3K36me3"))
        desired_histone_mark_features <- c("", desired_histone_mark_features[rowSums(sapply(peak_types, function(peak_type) grepl(peak_type, desired_histone_mark_features))) > 0])
        desired_CADD_features <- c("", get_features_by_group("CADD"))
        
        desired_tests <- apply(expand.grid(desired_variant_type_features, desired_region_type_features, desired_gene_set_features, desired_RBP_features, desired_histone_mark_features), 1, function(x) paste0(x[x != ""], collapse="+")) #desired_gene_set_features
        desired_tests <- unique(c(desired_tests, apply(expand.grid("is_snv", desired_CADD_features, desired_region_type_features, desired_gene_set_features, desired_RBP_features), 1, function(x) paste0(x[x != ""], collapse="+"))))
    }
    num_tests = length(desired_tests)
    all_features <- unique(unlist(strsplit(desired_tests, "\\+")))
    #num_independent_tests = calculate_num_independent_tests(desired_tests, num_variants_per_iteration=10000, num_iterations=100, random_gnomad_variants=random_gnomad_variants)
    
    ####################################################################################################################
    
    # Annotate cases and controls
    print("Annotating cases and controls...")
    cases <- annotate(cases, "nearest_gene"); controls <- annotate(controls, "nearest_gene"); 
    print("Annotating variants with features...")
    all_annotations <- annotate(rbind(cases, controls), all_features) #, variants_granges=c(case_granges, control_granges))
    cases <- all_annotations[1:nrow(cases),]
    controls <- all_annotations[(nrow(cases)+1):nrow(all_annotations),]
    #cases <- unfactorize(cases)
    #controls <- unfactorize(controls)
    
    # Returns the formatted test to append to the tests data frame.
    append_test <- function(test_name, variant_type, d1, d0, notes="", phenotype_annotation=annotate_phenotypes, gene_hits_annotation=FALSE, variant_hits_annotation=FALSE, hg38_variant_hits_annotation=FALSE) {
        gene_colname = "nearest_gene"; chr_colname = "Chrom"; pos_colname = "Position"; ref_colname = "Ref"; alt_colname = "Alt"
        cases_have_gene_annotation = (gene_colname %in% colnames(d1)); controls_have_gene_annotation = (gene_colname %in% colnames(d0))
        notable_genes_annotation = (cases_have_gene_annotation && controls_have_gene_annotation) # If this is true, does annotation for notable genes.
        phenotype_burdens = ""; notable_genes = ""; cases_gene_hits = ""; controls_gene_hits = ""; cases_variant_hits = ""; controls_variant_hits = ""; cases_hg38_variant_hits = ""; controls_hg38_variant_hits = ""
        # Annotate variants
        if(variant_hits_annotation) {
            cases_variant_hits <- genomic_dat_to_strings(d1, chr_colname=chr_colname, pos_colname=pos_colname, ref_colname=ref_colname, alt_colname=alt_colname, notes_colname=gene_colname, add_chr_prefix=TRUE, sep=";", empty_string="-")
            controls_variant_hits <- genomic_dat_to_strings(d0, chr_colname=chr_colname, pos_colname=pos_colname, ref_colname=ref_colname, alt_colname=alt_colname, notes_colname=gene_colname, add_chr_prefix=TRUE, sep=";", empty_string="-")
        }
        if(hg38_variant_hits_annotation && "Chrom_hg38" %in% colnames(d1) && "Position_hg38" %in% colnames(d1) && "Chrom_hg38" %in% colnames(d0) && "Position_hg38" %in% colnames(d0)) {
            cases_hg38_variant_hits <- genomic_dat_to_strings(d1, chr_colname="Chrom_hg38", pos_colname="Position_hg38", ref_colname=ref_colname, alt_colname=alt_colname, notes_colname=gene_colname, add_chr_prefix=TRUE, sep=";", empty_string="-")
            controls_hg38_variant_hits <- genomic_dat_to_strings(d0, chr_colname="Chrom_hg38", pos_colname="Position_hg38", ref_colname=ref_colname, alt_colname=alt_colname, notes_colname=gene_colname, add_chr_prefix=TRUE, sep=";", empty_string="-")
        }
        # Annotate genes
        if(notable_genes_annotation || (cases_have_gene_annotation && gene_hits_annotation)) { 
            cases_gene_counts <- table(unlist(strsplit(paste0(d1[,gene_colname]), ",")))
            cases_gene_hits_names <- paste0(names(cases_gene_counts))
            cases_recurrent_genes_dummy <- (cases_gene_counts > 1)
            cases_gene_hits <- cases_gene_hits_names
            cases_gene_hits[cases_recurrent_genes_dummy] <- paste0(cases_gene_hits[cases_recurrent_genes_dummy], "(", cases_gene_counts[cases_recurrent_genes_dummy], ")")
        }
        if(notable_genes_annotation || (controls_have_gene_annotation && gene_hits_annotation)) { 
            controls_gene_counts <- table(unlist(strsplit(paste0(d0[,gene_colname]), ",")))
            controls_gene_hits_names <- paste0(names(controls_gene_counts))
            controls_recurrent_genes_dummy <- (controls_gene_counts > 1)
            controls_gene_hits <- controls_gene_hits_names
            controls_gene_hits[controls_recurrent_genes_dummy] <- paste0(controls_gene_hits[controls_recurrent_genes_dummy], "(", controls_gene_counts[controls_recurrent_genes_dummy], ")")
        }
        # Notable genes are those that show up in cases but not in controls, AND are either recurrent, constrained, known candidates, or highly expressed.
        if(notable_genes_annotation) {
            if(grepl("CDH", disease)) { known_candidates <- candidate_CDH_genes } else if(grepl("CHD", disease)) { known_candidates <- candidate_CHD_genes } else { known_candidates <- candidate_CDH_genes }
            notable_genes <- paste0(cases_gene_hits[!(cases_gene_hits_names %in% controls_gene_hits_names) & (cases_recurrent_genes_dummy | cases_gene_hits_names %in% unique(c(constrained_genes, known_candidates, high_expressed_genes)))], collapse=",")
        }
        
        if(phenotype_annotation) {
            if (grepl("CDH", disease)) { 
                all_phenotype_info <- cdh_phenotype_info; phenotype_info <- all_phenotype_info[all_phenotype_info$sample %in% d1$sample,]
                female_complex = sum(phenotype_info$gender == "F" & phenotype_info$class == "Complex")
                female_isolated = sum(phenotype_info$gender == "F" & phenotype_info$class == "Isolated")
                male_complex = sum(phenotype_info$gender == "M" & phenotype_info$class == "Complex")
                male_isolated = sum(phenotype_info$gender == "M" & phenotype_info$class == "Isolated")
                all_female_complex = sum(all_phenotype_info$gender == "F" & all_phenotype_info$class == "Complex")
                all_female_isolated = sum(all_phenotype_info$gender == "F" & all_phenotype_info$class == "Isolated")
                all_male_complex = sum(all_phenotype_info$gender == "M" & all_phenotype_info$class == "Complex")
                all_male_isolated = sum(all_phenotype_info$gender == "M" & all_phenotype_info$class == "Isolated")
                
                female_binom_result <- binomial_test(female_complex+female_isolated, male_complex+male_isolated, all_female_complex+all_female_isolated, all_male_complex+all_male_isolated, alternative="two.sided")
                isolated_binom_result <- binomial_test(female_isolated+male_isolated, female_complex+male_complex, all_female_isolated+all_male_isolated, all_female_complex+all_male_complex, alternative="two.sided")
                female_isolated_binom_result <- binomial_test(female_isolated, female_complex+male_complex+male_isolated, all_female_isolated, all_female_complex+all_male_complex+all_male_isolated, alternative="two.sided")
                female_complex_binom_result <- binomial_test(female_complex, female_isolated+male_complex+male_isolated, all_female_complex, all_female_isolated+all_male_complex+all_male_isolated, alternative="two.sided")
                male_isolated_binom_result <- binomial_test(male_isolated, male_complex+female_complex+female_isolated, all_male_isolated, all_male_complex+all_female_complex+all_female_isolated, alternative="two.sided")
                male_complex_binom_result <- binomial_test(male_complex, male_isolated+female_complex+female_isolated, all_male_complex, all_male_isolated+all_female_complex+all_female_isolated, alternative="two.sided")
                phenotype_annotation_elements <- c()
                for(phenotype_test_name in c("female", "isolated", "female_isolated", "female_complex", "male_isolated", "male_complex")) {
                    phenotype_test <- get(paste0(phenotype_test_name,"_binom_result"))
                    phenotype_count = phenotype_test[["m1"]]
                    phenotype_enrichment = phenotype_test[["estimate"]]
                    phenotype_p.value = phenotype_test[["p.value"]]
                    phenotype_annotation_elements <- c(phenotype_annotation_elements, paste0(phenotype_test_name, ": ", phenotype_count, " samples, ", round(phenotype_enrichment,2), "x (p=", formatC(phenotype_p.value,format="e",digits=2), ")"))
                }
            }
            phenotype_burdens <- paste0(phenotype_annotation_elements, collapse=", ")
        }
        
        if(gene_hits_annotation) { 
            if(length(cases_gene_hits)>0) { cases_gene_hits <- paste0(cases_gene_hits, collapse=",")
            } else { cases_gene_hits <- "-" }
            if(length(controls_gene_hits)>0) { controls_gene_hits <- paste0(controls_gene_hits, collapse=",") 
            } else { controls_gene_hits <- "-" }
        } else { cases_gene_hits = ""; controls_gene_hits = "" }
        
        return(c(test_name, variant_type, nrow(d1), nrow(d0), phenotype_burdens, notes, notable_genes, cases_gene_hits, controls_gene_hits, cases_variant_hits, controls_variant_hits, cases_hg38_variant_hits, controls_hg38_variant_hits))
    }
    
    # Add remaining tests to the data frame
    tests <- unfactorize(data.frame(t(sapply(1:num_tests, function(i) {
        print(paste0("Appending test ",desired_tests[i],"...[",i," / ",num_tests,"]"))
        test_elements <- unlist(strsplit(desired_tests[i], "\\+"))
        num_test_elements = length(test_elements)
        annotations_cases <- cases[,which(colnames(cases) %in% test_elements)]; annotations_controls <- controls[,which(colnames(controls) %in% test_elements)]
        if(num_test_elements > 1) { 
            hits_cases <- rowSums(annotations_cases == "Y") == ncol(annotations_cases)
            hits_controls <- rowSums(annotations_controls == "Y") == ncol(annotations_controls)
        } else { 
            hits_cases <- annotations_cases == "Y"
            hits_controls <- annotations_controls == "Y"
        }
        
        # Determine variant type
        if("is_indel" %in% test_elements) { variant_type = "indels"
        } else if("is_snv" %in% test_elements) { 
            variant_type = "SNPs";
            if(subsample_snvs) { 
                hits_cases <- hits_cases & (1:length(hits_cases) %in% sample(1:length(hits_cases), sum(cases$is_indel == "Y"), replace=FALSE))
                hits_controls <- hits_controls & (1:length(hits_controls) %in% sample(1:length(hits_controls), sum(controls$is_indel == "Y"), replace=FALSE)) 
            }
        } else { variant_type = "SNP+indel" }
        
        # Determine test name
        gene_set_feature = intersect(test_elements, desired_gene_set_features); region_feature = intersect(test_elements, desired_region_type_features); RBP_feature = intersect(test_elements, desired_RBP_features); histone_mark_feature = intersect(test_elements, desired_histone_mark_features); CADD_feature = intersect(test_elements, desired_CADD_features)
        if(length(gene_set_feature)>0) { test_name = gene_set_feature } else { test_name = NULL }
        if(length(CADD_feature)>0) { test_name = paste(c(test_name, CADD_feature), collapse=" ") }
        if(length(region_feature)>0) { test_name = paste(c(test_name, region_feature), collapse=" ") }
        if(length(RBP_feature)>0) { test_name = paste(c(test_name, RBP_feature), collapse=" ") }
        
        if(length(histone_mark_feature)>0) {
            histone_feature_split <- strsplit(histone_mark_feature, "\\.")[[1]]
            eid = histone_feature_split[1]
            if (eid == "H3K27ac_Dickel_heart_enhancer") { tissue = "heart"; histone_mark = "H3K27ac"; peak_type = "Dickel" 
            } else if(histone_feature_split[2] == "tissue_majority") { eid = "tissue_majority"; tissue = "tissue_majority"; histone_mark = histone_feature_split[1]; peak_type = histone_feature_split[3]
            } else { tissue = get_roadmap_epigenome_names(eid); histone_mark = histone_feature_split[2]; peak_type = histone_feature_split[3] }
            if(tissue == "-") { tissue <- eid; eid <- paste0(relevant_eids, collapse=",") } # special case for combined tissue features, where the eid is initially "relevant"
            test_name = paste(c(test_name, paste0(histone_mark, " in ", tissue, " (", eid, ", ", peak_type, ")")), collapse=" ")
            histone_info = get_histone_mark_function(histone_mark)
        } else { histone_info = "" }
        if(is.null(test_name) || test_name == "") {
            test_name = "overall" 
            histone_info = paste0(round(sum(hits_cases)/get(paste0(toupper(disease),"_SAMPLE_COUNT")),2), " variants/case, ", round(sum(hits_controls)/get(paste0(toupper(control_name),"_SAMPLE_COUNT")),2), " variants/control, two-sided binomial test")
        }
        
        return(append_test(test_name, variant_type, cases[hits_cases,], controls[hits_controls,], histone_info, gene_hits_annotation=TRUE, variant_hits_annotation=TRUE, hg38_variant_hits_annotation=include_hg38_variant_annotation))
    }))))
    colnames(tests) <- c("test", "variants", "m1", "m0", "phenotype_burdens", "notes", "notable_genes", "cases_gene_hits", "controls_gene_hits", "cases_hg19_variant_hits", "controls_hg19_variant_hits", "cases_hg38_variant_hits", "controls_hg38_variant_hits")
    overall_burden_tests <- tests$test == "overall"
    if(sum(overall_burden_tests)>0) { tests <- tests[c(which(overall_burden_tests),which(!overall_burden_tests)),] }
    
    tests$m1 <- as.numeric(tests$m1); tests$m0 <- as.numeric(tests$m0)
    binom_enrichments <- c()
    binom_p.values <- c()
    overall_norm_factors <- c(1, 1, 1); names(overall_norm_factors) <- c("SNP+indel", "SNPs", "indels")
    #if (normalize_by_overall_burden) { 
    #    overall_norm_factors[1] <- (tests$m1[1]/n1)/(tests$m0[1]/n0) 
    #    overall_norm_factors[2] <- (tests$m1[2]/n1)/(tests$m0[2]/n0) 
    #    overall_norm_factors[3] <- (tests$m1[3]/n1)/(tests$m0[3]/n0)
    #}
    for(i in 1:nrow(tests)) {
        if (statistical_test != "poisson regression" && tests[i,1] == "overall") { 
            binom_result <- binomial_test(tests$m1[i], tests$m0[i], get(paste0(toupper(disease),"_SAMPLE_COUNT")), get(paste0(toupper(control_name),"_SAMPLE_COUNT")), alternative="two.sided")
        } else if (statistical_test == "binomial") { norm_factor = overall_norm_factors[tests[i,2]]; binom_result <- binomial_test(floor(tests$m1[i]/norm_factor), tests$m0[i], n1, n0, alternative="greater") 
        } else if (statistical_test == "poisson regression") { 
            if (tests[i,1] == "overall") {
                tests[i,6] <- gsub(", two-sided binomial test", "", tests[i,6])
                variants = tolower(tests[i,2]); if(variants == "snp+indel") { variants = "muts" }
                case_samples <- table(get(paste0("cases_",variants))$sample)
                control_samples <- table(get(paste0("controls_",variants))$sample)
            } else {
                case_samples <- table(paste0(sapply(unlist(strsplit(tests$cases_variant_hits[i],";")), function(x) strsplit(x,"_")[[1]][1])))
                control_samples <- table(paste0(sapply(unlist(strsplit(tests$controls_variant_hits[i],";")), function(x) strsplit(x,"_")[[1]][1])))
            }
            case_samples <- cbind(names(case_samples), case_samples); colnames(case_samples) <- c("sample", "variant_count")
            control_samples <- cbind(names(control_samples), control_samples);  colnames(control_samples) <- c("sample", "variant_count")
            if(grepl("CHD",disease)) { case_paternal_ages <- get("chd_parental_age") } else { case_paternal_ages <- get("cdh_parental_age") }
            case_paternal_ages <- case_paternal_ages[case_paternal_ages$Blinded.ID %in% cases_muts$sample,]
            control_paternal_ages <- get("ssc_parental_age")[,1:2]; control_paternal_ages <- control_paternal_ages[control_paternal_ages$Blinded.ID %in% controls_muts$sample,]
            case_dat <- merge(case_paternal_ages, case_samples, by.x="Blinded.ID", by.y="sample", all.x=TRUE, all.y=TRUE)
            case_dat$variant_count <- as.numeric(paste0(case_dat$variant_count)); case_dat$Paternal.Age.at.Proband.Birth <- as.numeric(paste0(case_dat$Paternal.Age.at.Proband.Birth))
            case_dat$variant_count[which(is.na(case_dat$variant_count))] <- 0; case_dat$Paternal.Age.at.Proband.Birth[which(is.na(case_dat$Paternal.Age.at.Proband.Birth))] <- mean(case_dat$Paternal.Age.at.Proband.Birth[!is.na(case_dat$Paternal.Age.at.Proband.Birth)])
            control_dat <- merge(control_paternal_ages, control_samples, by.x="Blinded.ID", by.y="sample", all.x=TRUE, all.y=TRUE)
            control_dat$variant_count <- as.numeric(paste0(control_dat$variant_count)); control_dat$Paternal.Age.at.Proband.Birth <- as.numeric(paste0(control_dat$Paternal.Age.at.Proband.Birth))
            control_dat$variant_count[which(is.na(control_dat$variant_count))] <- 0; control_dat$Paternal.Age.at.Proband.Birth[which(is.na(control_dat$Paternal.Age.at.Proband.Birth))] <- mean(control_dat$Paternal.Age.at.Proband.Birth[!is.na(control_dat$Paternal.Age.at.Proband.Birth)])
            variant_count_dat <- rbind(case_dat, control_dat)
            case_control <- rep(0, nrow(variant_count_dat)); case_control[1:nrow(case_dat)] <- 1; variant_count_dat <- cbind(variant_count_dat, case_control)
            
            m1 <- glm(variant_count ~ case_control + Paternal.Age.at.Proband.Birth, family="poisson", data=variant_count_dat) # control = list(maxit = 50)
            #binom_enrichments <- c(binom_enrichments, coef(summary(m1))[2,1])
            #binom_p.values <- c(binom_p.values, coef(summary(m1))[2,4])
            binom_result <- new.env()
            binom_result[["estimate"]] <- coef(summary(m1))[2,1]
            binom_result[["p.value"]] <- coef(summary(m1))[2,4]
            # paste0(c("Controls model: log(variant_count) = ", " * paternal_age + "), rev(m1$coefficients), collapse="")
            # #plot(control_dat$Paternal.Age.at.Proband.Birth, control_dat$variant_count)
            # #sum(predict(m1, type="response"))
            # predicted_case_variant_count = round(sum(predict(m1, newdata=case_dat, type="response")))
            # m2 <- glm(variant_count ~ Paternal.Age.at.Proband.Birth, family="poisson", data=case_dat)
            # #summary(m2)
            # paste0(c("Cases model: log(variant_count) = ", " * paternal_age + "), rev(m2$coefficients), collapse="")
            # #plot(case_dat$Paternal.Age.at.Proband.Birth, case_dat$variant_count)
            # #sum(predict(m2, type="response"))
            # predicted_control_variant_count = round(sum(predict(m2, newdata=control_dat, type="response")))
            
            #binom_result <- fisher_exact_test(tests$m1[i], predicted_case_variant_count, n1, n1, alternative = c("greater"))
        } else if (statistical_test == "fisher_exact_2s") {
            n1_i = n1; n0_i = n0
            if(variant_type_testing && tests$variants[i] != "SNP+indel") {
                if(tests$variants[i] == "SNPs") { n1_i = sum(cases$is_snv == "Y"); n0_i = sum(controls$is_snv == "Y")
                } else { n1_i = sum(cases$is_indel == "Y"); n0_i = sum(controls$is_indel == "Y") }
            }
            binom_result <- fisher_exact_test(tests$m1[i], tests$m0[i], n1_i, n0_i, alternative="two.sided")
        } else { binom_result <- fisher_exact_test(tests$m1[i], tests$m0[i], n1, n0, alternative="greater") } # fisher_exact
        binom_enrichments <- c(binom_enrichments, binom_result[["estimate"]])
        binom_p.values <- c(binom_p.values, binom_result[["p.value"]])
    }
    burdens <- cbind(tests$test, tests$variants, tests$m1, tests$m0, binom_enrichments, binom_p.values, n1, n0, tests$phenotype_burdens, tests$notes, tests$notable_genes, tests$cases_gene_hits, tests$controls_gene_hits, tests$cases_hg19_variant_hits, tests$controls_hg19_variant_hits, tests$cases_hg38_variant_hits, tests$controls_hg38_variant_hits)
    colnames(burdens) <- c("burden", "variants", "m1", "m0", "enrichment", "p.value", "n1", "n0", "phenotype_burdens", "notes", "notable_genes", "cases_gene_hits", "controls_gene_hits", "cases_hg19_variant_hits", "controls_hg19_variant_hits", "cases_hg38_variant_hits", "controls_hg38_variant_hits")
    if(statistical_test == "poisson regression") { colnames(burdens)[5] <- "z_score" }
    
    # Return burden results
    return(data.frame(burdens))
}
#chdfb_burdens[chdfb_burdens$burden == "RBP",1:8]
#chdfb_burdens[chdfb_burdens$burden == "H3K36me3 in H9 Cells (E008, narrowPeak)",1:8]
#chdfb_burdens[chdfb_burdens$burden == "RBP H3K36me3 in H9 Cells (E008, narrowPeak)",1:8]
#chdfb_burdens[chdfb_burdens$burden == "TSS RBP H3K36me3 in H9 Cells (E008, narrowPeak)",1:8]
#chdfb_burdens[grepl("H3K79me2", chdfb_burdens$burden),1:8]
chdfb_burdens_all <- burden_analysis2("CHDFB_combined", ssc="SSCFB_combined", include_individual_RBPs=TRUE)
best_result_chdfb_all <- volcano_plot("CHDFB", chdfb_burdens_all[-c(1:3),], strong_effect_cutoff=0.9, bonferroni_p.value_cutoff=0.05/97, num_best_hits=20, pval_relative_importance=1000, number_labels=TRUE, x_max=4.5, label_gap=0.5, label_points_max_iterations=20)
chdfb_burdens <- burden_analysis2("CHDFB_combined", ssc="SSCFB_combined")
best_result_chdfb <- volcano_plot("CHDFB", chdfb_burdens[-c(1:3),], strong_effect_cutoff=0.9, bonferroni_p.value_cutoff=0.05/93, num_best_hits=10, pval_relative_importance=1000, number_labels=FALSE)

cdh_burdens_all <- burden_analysis2("CDH_combined", ssc="SSC_all", include_individual_RBPs=FALSE)
best_result_cdh_all <- volcano_plot("CDH", cdh_burdens_all[-c(1:3),], strong_effect_cutoff=0.9, bonferroni_p.value_cutoff=0.05/97, num_best_hits=20, pval_relative_importance=1000, number_labels=TRUE, x_max=3.5, label_gap=0.5, label_points_max_iterations=20)

cdh3_burdens <- burden_analysis2("CDH3", ssc="SSC_all", include_individual_RBPs=FALSE)
best_result_cdh3 <- volcano_plot("CDH3", cdh3_burdens[-c(1:3),], strong_effect_cutoff=0.9, bonferroni_p.value_cutoff=0.05/97, num_best_hits=20, pval_relative_importance=1000, number_labels=TRUE, x_max=3.5, label_gap=0.5, label_points_max_iterations=20)

cdh_chd_burdens_all <- burden_analysis2("CDH_CHD", ssc="SSC_all", include_individual_RBPs=TRUE)
best_result_cdh_chd_all <- volcano_plot("CDH+CHD", cdh_chd_burdens_all[-c(1:3),], strong_effect_cutoff=0.9, bonferroni_p.value_cutoff=0.05/97, num_best_hits=20, pval_relative_importance=1000, number_labels=TRUE, x_max=5, label_gap=0.5, label_points_max_iterations=20)


#################################################################################################################
# Analyze RBP binding site disruption by denovo nc mutations. 
# RBP eCLIP data from Chaolin has peaks for 160 RBPs, out of ~1500 (600 structurally distinct) total RBPs in human genome.

# Confirmed: ~500 RBPs in human genome that bind to single-stranded RNA binding domains, plus an additional 200 RBPs that bind to dsRNA.
# ~1900 RBPs total taking predicted ones in GO, though this includes many false positives.
# 1542 tested RBPs -> 600 structurally distinct ones (Gerstberger S, Hafner M, & Tuschl T, 2014. http://www.nature.com/nrg/journal/v15/n12/full/nrg3813.html)
# 160/600 ~ 26.7% is the eCLIP coverage we have.
#################################################################################################################
RBP_ANNOTATIONS_FOLDER <- data_path("RBP_annotations")
dir.create(file.path(RBP_ANNOTATIONS_FOLDER), showWarnings = FALSE)
RBP_OUTPUT_FOLDER <- output_path("RBP_analysis")
dir.create(file.path(RBP_OUTPUT_FOLDER), showWarnings = FALSE)

# Find the hits/overlaps of the given case and control variants with all 160 RBP eCLIP peaks files.
# The start parameter, if specified, indicates which file (of the 160) to begin processing at. It is for use if run needs to be restarted partway through.
# The pillows parameter is a vector of all pillow/padding values to try adding onto original eCLIP peaks and running overlaps for, in efficient nested fashion.
find_rbp_hits <- function(cases, controls, disease="", rbp_output_dir=NULL, pillows=0, cases_granges=NULL, controls_granges=NULL, start=1) {
    rbp_dir = data_path("ENCODE_eCLIP_peak")
    rbp_output_dir <- full_path(RBP_ANNOTATIONS_FOLDER, paste0(c("RBP",disease),collapse="_"))
    if (is.null(cases_granges)) { cases_granges <- to_genomic_regions(cases, chr_colname="Chrom", start_colname="Position", end_colname="Position", labels=genomic_coordinates_to_strings(cases$Chrom, cases$Position)) }
    if (is.null(controls_granges)) { controls_granges <- to_genomic_regions(controls, chr_colname="Chrom", start_colname="Position", end_colname="Position", labels=genomic_coordinates_to_strings(controls$Chrom, controls$Position)) }
    cases_granges_backup <- cases_granges
    controls_granges_backup <- controls_granges
    
    # combined_rbp_granges <- c()
    # rbp_footprints <- c()
    
    dir.create(file.path(rbp_output_dir), showWarnings = FALSE)
    
    for(rbp_file in list.files(path=rbp_dir)[start:160]) { # Example file with some negative strands: "HepG2.CSTF2T.R2.tag.uniq.peak.sig.bed"
        print(paste0(rbp_dir, FILEPATH_DELIM, rbp_file))
        rbp_dat <- read.csv(full_path(rbp_dir, rbp_file), sep="\t", header=FALSE)[,c(1:3,6)]
        #colnames(rbp_dat)[1:3] <- c("chromosome", "start", "end")
        colnames(rbp_dat)[1:4] <- c("chromosome", "start", "end", "strand")
        rbp_dat_backup <- rbp_dat
        cases_granges <- cases_granges_backup
        controls_granges <- controls_granges_backup
        
        #rbp_dat <- rbp_dat[rbp_dat$strand=="-",1:3] # Consider only minus strand!
        sapply_out <- sapply(1:nrow(rbp_dat), function(i) { rbp_dat[i,2:3] <- sort(c(rbp_dat$end[i], rbp_dat$start[i])) } )
        
        pillows <- sort(pillows, decreasing=TRUE)
        for(pillow_index in 1:length(pillows)) {
            pillow <- pillows[pillow_index]
            print(paste0("Pillow size: ", pillow))
            rbp_dat$start <- rbp_dat$start - pillow; rbp_dat$end <- rbp_dat$end + pillow
            rbp <- genomic_coordinates_to_strings(rbp_dat[,1], rbp_dat[,2], rbp_dat[,3])
            duplicate_rbp <- duplicated(rbp) | is.na(rbp_dat[,1])
            rbp <- rbp[!duplicate_rbp]
            rbp_dat <- rbp_dat[!duplicate_rbp,]
            rbp_granges <- to_genomic_regions(rbp_dat, labels=rbp)
            
            # 			if(pillow_index == 1) {
            #     			combined_rbp_granges <- c(combined_rbp_granges, rbp_granges)
            #     			combined_rbp_granges <- intersect(combined_rbp_granges, combined_rbp_granges)
            #     			rbp_footprints <- c(rbp_footprints, sum(end(rbp_granges) - start(rbp_granges)))
            #     			names(rbp_footprints)[length(rbp_footprints)] = strsplit(rbp_file,"\\.")[[1]][2]
            # 			}
            
            hits_cases <- unique(data.frame(olRanges(cases_granges, rbp_granges)))
            if(nrow(hits_cases) > 0) { 
                hits_cases <- cbind(hits_cases, genomic_coordinates_to_strings(hits_cases$seqnames, hits_cases$Sstart, hits_cases$Send))
                colnames(hits_cases)[ncol(hits_cases)] <- "RBP_peak"
                cases_rbp <- merge(cases, hits_cases, by.x=c("Chrom", "Position"), by.y=c("seqnames", "start"))[c("Chrom", "Position", "Ref", "Alt", "sample", "RBP_peak")]
                cases_rbp <- aggregate(cases_rbp, by=list(cases_rbp$Chrom, cases_rbp$Position), FUN=function(x) { paste0(x[!duplicated(x)], collapse=",") })[,c("Chrom", "Position", "Ref", "Alt", "sample", "RBP_peak")]  
            } else { cases_rbp <- data.frame("-", "-", "-", "-", "-", "-") }
            colnames(cases_rbp) <- c("Chrom", "Position", "Ref", "Alt", "sample", "RBP_peak")
            hits_controls <- unique(data.frame(olRanges(controls_granges, rbp_granges)))
            if(nrow(hits_controls) > 0) {
                hits_controls <- cbind(hits_controls, genomic_coordinates_to_strings(hits_controls$seqnames, hits_controls$Sstart, hits_controls$Send))
                colnames(hits_controls)[ncol(hits_controls)] <- "RBP_peak"
                controls_rbp <- merge(controls, hits_controls, by.x=c("Chrom", "Position"), by.y=c("seqnames", "start"))[c("Chrom", "Position", "Ref", "Alt", "sample", "RBP_peak")]
                controls_rbp <- aggregate(controls_rbp, by=list(controls_rbp$Chrom, controls_rbp$Position), FUN=function(x) { paste0(x[!duplicated(x)], collapse=",") })[,c("Chrom", "Position", "Ref", "Alt", "sample", "RBP_peak")]  
            } else { controls_rbp <- data.frame("-", "-", "-", "-", "-", "-") }
            colnames(controls_rbp) <- c("Chrom", "Position", "Ref", "Alt", "sample", "RBP_peak")
            
            rbp_output_sub_dir = paste0(pillow,"bp")
            dir.create(file.path(rbp_output_dir, rbp_output_sub_dir), showWarnings = FALSE)
            rbp_output_file = paste0(rbp_output_dir, FILEPATH_DELIM, rbp_output_sub_dir, FILEPATH_DELIM, strsplit(rbp_file, ".R2")[[1]][1], "_overlaps_cases.csv")
            write.csv(cases_rbp, file=rbp_output_file, row.names=FALSE)
            rbp_output_file <- gsub("cases", "controls", rbp_output_file)
            write.csv(controls_rbp, file=rbp_output_file, row.names=FALSE)
            
            hits_cases_strings <- cbind(genomic_coordinates_to_strings(hits_cases$seqnames, hits_cases$start), genomic_coordinates_to_strings(hits_cases$seqnames, hits_cases$Sstart, hits_cases$Send))
            hits_controls_strings <- cbind(genomic_coordinates_to_strings(hits_controls$seqnames, hits_controls$start), genomic_coordinates_to_strings(hits_controls$seqnames, hits_controls$Sstart, hits_controls$Send))
            cases_granges <- cases_granges[paste(cases_granges) %in% hits_cases_strings[,1],]
            controls_granges <- controls_granges[paste(controls_granges) %in% hits_controls_strings[,1],]
            rbp_dat <- rbp_dat[rbp %in% union(hits_cases_strings[,2], hits_controls_strings[,2]),]
            rbp_dat$start <- rbp_dat$start + pillow; rbp_dat$end <- rbp_dat$end - pillow
        }
    }
    # rbp_footprints <- c(rbp_footprints, sum(end(combined_rbp_granges) - start(combined_rbp_granges)))
    # names(rbp_footprints)[length(rbp_footprints)] = "combined"
    # return(rbp_footprints)
}
#find_rbp_hits(cdh, ssc, disease="CDH", cases_granges=cdh_granges, controls_granges=ssc_granges, pillows=seq(0,100,by=10))
#find_rbp_hits(chd, ssc, disease="CHD", cases_granges=chd_granges, controls_granges=ssc_granges, pillows=seq(0,100,by=10))
#find_rbp_hits(cdh2, ssc_1088, disease="CDH2", cases_granges=cdh2_granges, controls_granges=ssc_1088_granges, pillows=seq(0,100,by=10))
#find_rbp_hits(chd2, ssc_1088, disease="CHD2", cases_granges=chd2_granges, controls_granges=ssc_1088_granges, pillows=seq(0,100,by=10))
#find_rbp_hits(cdh, ssc_518, disease="CDH_SSC518", cases_granges=cdh_granges, controls_granges=ssc_518_granges, pillows=seq(0,100,by=10))
#find_rbp_hits(chdfb, sscfb, disease="CHDFB", cases_granges=chdfb_granges, controls_granges=sscfb_granges, pillows=seq(0,100,by=10))
#find_rbp_hits(chdfb2, sscfb2, disease="CHDFB2", cases_granges=chdfb2_granges, controls_granges=sscfb2_granges, pillows=seq(0,100,by=10))
#find_rbp_hits(chd2, ssc_1088, disease="CHD2", cases_granges=chd2_granges, controls_granges=ssc_1088_granges, pillows=seq(0,100,by=10))
#find_rbp_hits(regulatory_variant_dat, ssc_518, disease="regulatory", cases_granges=regulatory_variant_granges, controls_granges=ssc_518_granges, pillows=seq(0,100,by=10))
chdfb_combined <- chdfb_combined_muts_wgsa[1,1:8]; colnames(chdfb_combined)[1:4] <- c("Chrom", "Position", "Ref", "Alt")
sscfb_combined <- sscfb_combined_muts_wgsa[1,1:8]; colnames(sscfb_combined)[1:4] <- c("Chrom", "Position", "Ref", "Alt")
rbp_footprints <- find_rbp_hits(chdfb_combined, sscfb_combined, disease="test")



# Other RBP analysis functions 
if (TRUE) {
    # Plot RBP eCLIP peak count and length and footprint/coverage distribution statistics.
    make_rbp_stats <- function() {
        lengths <- c(); counts <- c(); rbp_footprints <- c(); combined_rbp_granges <- c()
        rbp_dir = data_path("ENCODE_eCLIP_peak")
        rbp_files <- list.files(path=rbp_dir)[1:160]
        for(i in 1:length(rbp_files)) {
            rbp_file <- rbp_files[i]
            rbp_dat <- read.csv(full_path(rbp_dir, rbp_file), sep="\t", header=FALSE)[,c(1:3,6)]
            colnames(rbp_dat)[1:4] <- c("chromosome", "start", "end", "strand")
            
            print(paste0(i, "    ", rbp_file, " : ", nrow(rbp_dat), " eCLIP peaks"))
            #sapply_out <- sapply(1:nrow(rbp_dat), function(i) { rbp_dat[i,2:3] <- sort(c(rbp_dat$end[i], rbp_dat$start[i])) } )
            length <- abs(rbp_dat$end-rbp_dat$start)
            length_density <- density(length)
            lengths <- c(lengths, length)
            counts <- c(counts, nrow(rbp_dat))
            
            rbp <- genomic_coordinates_to_strings(rbp_dat[,1], rbp_dat[,2], rbp_dat[,3])
            duplicate_rbp <- duplicated(rbp) | is.na(rbp_dat[,1])
            rbp <- rbp[!duplicate_rbp]
            rbp_dat <- rbp_dat[!duplicate_rbp,]
            rbp_granges <- to_genomic_regions(rbp_dat, labels=rbp)
            
            if(length(combined_rbp_granges) == 0) { combined_rbp_granges <- rbp_granges } else { combined_rbp_granges <- c(combined_rbp_granges, rbp_granges) }
            combined_rbp_granges <- intersect(combined_rbp_granges, combined_rbp_granges)
            rbp_footprints <- c(rbp_footprints, sum(end(rbp_granges) - start(rbp_granges)))
            names(rbp_footprints)[length(rbp_footprints)] = strsplit(rbp_file,"\\.")[[1]][2]
        }
        # Plot peak counts distribution
        pdf(output_path("rbp_peak_count_distribution.pdf"))
        plot(density(counts, from=0), main="RBP eCLIP Peak Count Distribution", xlab="peak count", xaxs="i")
        mtext(paste0("160 RBPs, ", length(lengths), " peaks, mean(peak_count) = ", round(mean(counts),0)))
        dev.off()
        # Plot peak lengths distribution
        pdf(output_path("rbp_peak_length_distribution.pdf"))
        plot(density(lengths, from=0, adjust=4), main="RBP eCLIP Peak Length Distribution", xlab="bp", xaxs="i")
        mtext(paste0("160 RBPs, ", length(lengths), " peaks, mean(peak_length) = ", round(mean(lengths),3)))
        dev.off()
        
        rbp_footprints <- c(rbp_footprints, sum(end(combined_rbp_granges) - start(combined_rbp_granges)))
        names(rbp_footprints)[length(rbp_footprints)] = "combined"
        # Plot peak footprint distribution
        pdf(output_path("rbp_peak_footprint_distribution.pdf"))
        plot(density(rbp_footprints[-length(rbp_footprints)], from=0), main="RBP eCLIP Peak Footprint Distribution", xlab="bp", xaxs="i")
        mtext(paste0("160 RBPs, mean(peak_footprint) = ", round(mean(rbp_footprints[-length(rbp_footprints)])), ", combined_peak_footprint = ", rbp_footprints[length(rbp_footprints)]))
        dev.off()
        return(rbp_footprints)
    }
    rbp_footprints <- make_rbp_stats()
    
    # Load RBP hits/annotations for the given disease and pillow size, and return an environment containing them for both cases and controls.
    load_rbp_hits <- function(disease, pillow=NULL, verbose=FALSE) {
        if(is.null(pillow)) { rbp_output_dir <- paste0(RBP_ANNOTATIONS_FOLDER, FILEPATH_DELIM, paste0(c("RBP",disease),collapse="_"))
        } else { rbp_output_dir <- paste0(RBP_ANNOTATIONS_FOLDER, FILEPATH_DELIM, paste0(c("RBP",disease),collapse="_"), FILEPATH_DELIM, pillow, "bp") }
        
        rbp_results_cases <- data.frame()
        rbp_results_controls <- data.frame()
        for(rbp_file in list.files(path=rbp_output_dir)) {
            if(verbose) { print(rbp_file) }
            rbp_dat <- cbind(read.csv(paste0(rbp_output_dir, FILEPATH_DELIM, rbp_file), stringsAsFactors=FALSE), rbp_file)
            if(grepl("cases", rbp_file)) {
                rbp_results_cases <- rbind(rbp_results_cases, rbp_dat)
            } else if(grepl("controls", rbp_file)) {
                rbp_results_controls <- rbind(rbp_results_controls, rbp_dat)
            }
        }
        # Remove empty entries.
        rbp_results_cases <- rbp_results_cases[rbp_results_cases$Chrom != "-",]
        rbp_results_controls <- rbp_results_controls[rbp_results_controls$Chrom != "-",]
        # Turn "T" nucleotides, incorrectly converted to "TRUE" by read.csv, back into "T".
        rbp_results_cases$Ref[rbp_results_cases$Ref == "TRUE"] <- "T"
        rbp_results_cases$Alt[rbp_results_cases$Alt == "TRUE"] <- "T"
        rbp_results_controls$Ref[rbp_results_controls$Ref == "TRUE"] <- "T"
        rbp_results_controls$Alt[rbp_results_controls$Alt == "TRUE"] <- "T"
        # Append RBP annotation to results data frames.
        rbp_results_cases <- cbind(rbp_results_cases, sapply(strsplit(paste0(rbp_results_cases$rbp_file), "\\.|_overlaps", perl=TRUE), function(x) { print(x[2]) }) )
        rbp_results_controls <- cbind(rbp_results_controls , sapply(strsplit(paste0(rbp_results_controls$rbp_file), "\\.|_overlaps", perl=TRUE), function(x) { print(x[2]) }) )
        colnames(rbp_results_cases)[ncol(rbp_results_cases)] <- "RBP"
        colnames(rbp_results_controls)[ncol(rbp_results_controls)] <- "RBP"
        # Save results to new environment and return it.
        return_env <- new.env()
        return_env[["cases"]] <- rbp_results_cases
        return_env[["controls"]] <- rbp_results_controls
        return(return_env)
    }
    
    # Binomial test: Null is binom(m, p), where m = m1 + m0 (m0 is the number of such variants in controls),  p = n1/(n1+n0),   
    # n1 is the number of case samples in total, n0 is the number of control samples in total
    # Run this test for SNV+Indel and then separately for the two groups, and return all results in a table.
    run_rbp_binomial_tests <- function(rbp_results_cases, rbp_results_controls, cases, controls, n0=0, n1=0) {
        rbp_results_cases <- unique_variants(rbp_results_cases, aggregate=FALSE)
        rbp_results_controls <- unique_variants(rbp_results_controls, aggregate=FALSE)
        cases <- unique_variants(cases, aggregate=FALSE)
        controls <- unique_variants(controls, aggregate=FALSE)
        
        results_dat <- data.frame()
        if(n0 == 0) { n0 <- length(unique(controls$sample)) }
        if(n1 == 0) { n1 <- length(unique(cases$sample)) }
        binom_prob = n1/(n1+n0) # cases_sample_count/(cases_sample_count+controls_sample_count)
        
        # Binomial test for all variants
        m0 = nrow(rbp_results_controls)
        m1 = nrow(rbp_results_cases)
        controls_variant_count = nrow(controls)	
        cases_variant_count = nrow(cases)
        binom_test_result <- binomial_test(m1, m0, n1, n0, alternative="greater")
        all_binom_test_result <- binomial_test(cases_variant_count, controls_variant_count, n1, n0, alternative="greater")
        results_dat <- rbind(results_dat, cbind("SNV + indel", binom_test_result[["estimate"]], binom_test_result[["p.value"]], m1, cases_variant_count, n1, m0, controls_variant_count, n0, all_binom_test_result[["estimate"]], all_binom_test_result[["p.value"]]))
        
        # Binomial test for just indels
        rbp_indels_cases <- (nchar(rbp_results_cases$Ref)>1 | nchar(rbp_results_cases$Alt)>1) 
        rbp_indels_controls <- (nchar(rbp_results_controls$Ref)>1 | nchar(rbp_results_controls$Alt)>1) 
        indels_cases <- (nchar(paste(cases$Ref))>1 | nchar(paste(cases$Alt))>1) 
        indels_controls <- (nchar(paste(controls$Ref))>1 | nchar(paste(controls$Alt))>1) 
        m0 = nrow(rbp_results_controls[rbp_indels_controls,])
        m1 = nrow(rbp_results_cases[rbp_indels_cases,])
        controls_variant_count = nrow(controls[indels_controls,])	
        cases_variant_count = nrow(cases[indels_cases,])
        binom_test_result <- binomial_test(m1, m0, n1, n0, alternative="greater")
        all_binom_test_result <- binomial_test(cases_variant_count, controls_variant_count, n1, n0, alternative="greater")
        results_dat <- rbind(results_dat, cbind("indels only", binom_test_result[["estimate"]], binom_test_result[["p.value"]], m1, cases_variant_count, n1, m0, controls_variant_count, n0, all_binom_test_result[["estimate"]], all_binom_test_result[["p.value"]]))
        
        # Binomial test for just SNVs
        m0 = nrow(rbp_results_controls[!rbp_indels_controls,])
        m1 = nrow(rbp_results_cases[!rbp_indels_cases,])
        controls_variant_count = nrow(controls[!indels_controls,])
        cases_variant_count = nrow(cases[!indels_cases,])
        binom_test_result <- binomial_test(m1, m0, n1, n0, alternative="greater")
        all_binom_test_result <- binomial_test(cases_variant_count, controls_variant_count, n1, n0, alternative="greater")
        results_dat <- rbind(results_dat, cbind("SNVs only", binom_test_result[["estimate"]], binom_test_result[["p.value"]], m1, cases_variant_count, n1, m0, controls_variant_count, n0, all_binom_test_result[["estimate"]], all_binom_test_result[["p.value"]]))
        
        colnames(results_dat) <- c("variant_type", "enrichment", "p.value", "cases_RBP", "num_case_variants", "num_case_samples", "controls_RBP", "num_control_variants", "num_control_samples", "all_enrichment", "all_p.value")
        return(results_dat)
    }
    
    # Analyze RBP binding peak disruption for the specified disease and peak pillow/padding size.
    # The default analysis is to calculate "quick"/overall RBP nc burden between cases and controls.
    # The analysis parameter can also be specified as "expression", "constrained", and "all", to further do breakdown by expression category, pLI, and both, respectively. 
    analyze_RBPs <- function(disease, analysis="quick", pillow=NULL) {
        expressed_TS <- new.env()
        if (grepl("CDH_new", disease)) { 
            cases <- cdh_new; controls <- ssc; cases_granges <- cdh_new_granges; controls_granges <- ssc_granges; expression_string_token <- "DE"; n0 <- SSC_SAMPLE_COUNT; n1 <- CDH_NEW_SAMPLE_COUNT
            expressed_TS[["L_TSS"]] <- LDE_TSS_granges; expressed_TS[["ML_TSS"]] <- MLDE_TSS_granges; expressed_TS[["MH_TSS"]] <- MHDE_TSS_granges; expressed_TS[["H_TSS"]] <- HDE_TSS_granges
            expressed_TS[["L_TES"]] <- LDE_TES_granges; expressed_TS[["ML_TES"]] <- MLDE_TES_granges; expressed_TS[["MH_TES"]] <- MHDE_TES_granges; expressed_TS[["H_TES"]] <- HDE_TES_granges
        } else if (grepl("CDH", disease)) { 
            cases <- cdh; controls <- ssc; cases_granges <- cdh_granges; controls_granges <- ssc_granges; expression_string_token <- "DE"; n0 <- SSC_SAMPLE_COUNT; n1 <- CDH_SAMPLE_COUNT
            expressed_TS[["L_TSS"]] <- LDE_TSS_granges; expressed_TS[["ML_TSS"]] <- MLDE_TSS_granges; expressed_TS[["MH_TSS"]] <- MHDE_TSS_granges; expressed_TS[["H_TSS"]] <- HDE_TSS_granges
            expressed_TS[["L_TES"]] <- LDE_TES_granges; expressed_TS[["ML_TES"]] <- MLDE_TES_granges; expressed_TS[["MH_TES"]] <- MHDE_TES_granges; expressed_TS[["H_TES"]] <- HDE_TES_granges
        } else if (grepl("CHD", disease)) { 
            cases <- chd; controls <- sscfb; cases_granges <- chd_granges; controls_granges <- sscfb_granges; expression_string_token <- "HE"; n0 <- SSCFB_SAMPLE_COUNT; n1 <- CHD_SAMPLE_COUNT
            expressed_TS[["L_TSS"]] <- LHE_TSS_granges; expressed_TS[["ML_TSS"]] <- MLHE_TSS_granges; expressed_TS[["MH_TSS"]] <- MHDE_TSS_granges; expressed_TS[["H_TSS"]] <- HHE_TSS_granges
            expressed_TS[["L_TES"]] <- LHE_TES_granges; expressed_TS[["ML_TES"]] <- MLHE_TES_granges; expressed_TS[["MH_TES"]] <- MHDE_TES_granges; expressed_TS[["H_TES"]] <- HHE_TES_granges
        } else { return(NULL) }
        
        if(is.null(pillow)) { pillow_string <- ""
        } else { pillow_string <- paste0("_", pillow, "bp") }
        
        rbp_annotations <- load_rbp_hits(disease, pillow)
        rbp_results_cases <- rbp_annotations[["cases"]]
        rbp_results_controls <- rbp_annotations[["controls"]]
        burden_results <- run_rbp_binomial_tests(rbp_results_cases, rbp_results_controls, cases, controls, n0=n0, n1=n1)
        write.csv(burden_results, file=full_path(RBP_OUTPUT_FOLDER, paste0("rbp_burden_", disease, pillow_string, ".csv")), row.names=FALSE)
        
        if(analysis!="quick") {
            dir.create(file.path(full_path(RBP_OUTPUT_FOLDER, "rbp_hits")), showWarnings = FALSE)
            rbp_results_cases_granges <- to_genomic_regions(rbp_results_cases, chr_colname="Chrom", start_colname="Position", label_colname="RBP")
            rbp_results_controls_granges <- to_genomic_regions(rbp_results_controls, chr_colname="Chrom", start_colname="Position", label_colname="RBP")
            if(analysis=="all" || "expression" %in% analysis) {
                for(expressed_TS_type in ls(expressed_TS)) {
                    expressed_TS_name <- gsub("_", paste0(expression_string_token, "_"), expressed_TS_type)
                    print(expressed_TS_name)
                    rbp_results_cases_subset <- rbp_results_cases[is_TS(rbp_results_cases_granges, TS_granges=expressed_TS[[expressed_TS_type]]),]
                    rbp_results_controls_subset <- rbp_results_controls[is_TS(rbp_results_controls_granges, TS_granges=expressed_TS[[expressed_TS_type]]),]
                    cases_subset <- cases[is_TS(cases_granges, TS_granges=expressed_TS[[expressed_TS_type]]),]
                    controls_subset <- controls[is_TS(controls_granges, TS_granges=expressed_TS[[expressed_TS_type]]),]
                    burden_results <- run_rbp_binomial_tests(rbp_results_cases_subset, rbp_results_controls_subset, cases_subset, controls_subset, n0=n0, n1=n1)
                    write.csv(burden_results, file=full_path(RBP_OUTPUT_FOLDER, paste0("rbp_burden_", disease, pillow_string, "_", expressed_TS_name, ".csv")), row.names=FALSE)
                    write.csv(rbp_results_cases_subset, file=full_path(RBP_OUTPUT_FOLDER, paste0("rbp_hits/rbp_case_hits_", disease, pillow_string, "_", expressed_TS_name, ".csv")), row.names=FALSE)
                    write.csv(rbp_results_controls_subset, file=full_path(RBP_OUTPUT_FOLDER, paste0("rbp_hits/rbp_control_hits_", disease, pillow_string, "_", expressed_TS_name, ".csv")), row.names=FALSE)
                }
            }
            if(analysis=="all" || "constrained" %in% analysis) {
                constraint_TS <- new.env()
                constraint_TS[["constrained_TSS"]] <- constrained_TSS_granges; constraint_TS[["constrained_TES"]] <- constrained_TES_granges;
                constraint_TS[["unconstrained_TSS"]] <- unconstrained_TSS_granges; constraint_TS[["unconstrained_TES"]] <- unconstrained_TES_granges;
                for(constraint_TS_type in ls(constraint_TS)) {
                    print(constraint_TS_type)
                    rbp_results_cases_subset <- rbp_results_cases[is_TS(rbp_results_cases_granges, TS_granges=constraint_TS[[constraint_TS_type]]),]
                    rbp_results_controls_subset <- rbp_results_controls[is_TS(rbp_results_controls_granges, TS_granges=constraint_TS[[constraint_TS_type]]),]
                    cases_subset <- cases[is_TS(cases_granges, TS_granges=constraint_TS[[constraint_TS_type]]),]
                    controls_subset <- controls[is_TS(controls_granges, TS_granges=constraint_TS[[constraint_TS_type]]),]
                    burden_results <- run_rbp_binomial_tests(rbp_results_cases_subset, rbp_results_controls_subset, cases_subset, controls_subset, n0=n0, n1=n1)
                    write.csv(burden_results, file=full_path(RBP_OUTPUT_FOLDER, paste0("rbp_burden_", disease, pillow_string, "_", constraint_TS_type, ".csv")), row.names=FALSE)
                    write.csv(rbp_results_cases_subset, file=full_path(RBP_OUTPUT_FOLDER, paste0("rbp_hits/rbp_case_hits_", disease, pillow_string, "_", constraint_TS_type, ".csv")), row.names=FALSE)
                    write.csv(rbp_results_controls_subset, file=full_path(RBP_OUTPUT_FOLDER, paste0("rbp_hits/rbp_control_hits_", disease, pillow_string, "_", constraint_TS_type, ".csv")), row.names=FALSE)
                }
            }
        }
    }
    
    # Make stats table showing RBP result comparison between 10bp pillow and regular peak runs, in a single table.
    make_RBP_comparison_table <- function(disease) {
        if (grepl("CDH", disease)) { 
            expression_string_token <- "DE"
        } else if (grepl("CHD", disease)) { 
            expression_string_token <- "HE"
        } else { return(NULL) }
        
        expression_bucket_names <- paste0(c("L", "ML", "MH", "H"), expression_string_token)
        for(region in c("TSS", "TES")) {
            summary_all_expression <- data.frame()
            for(expression_bucket_name in expression_bucket_names) {
                filename <- paste0("rbp_burden_", disease, "_", expression_bucket_name, "_", region, ".csv")
                filename_10bp <- gsub(disease, paste0(disease, "_10bp"), filename)
                output_file <- gsub(".csv", "_comparison.csv", filename)
                #print(output_file)
                dat <- read.csv(filename)
                dat_10bp <- read.csv(filename_10bp)
                if (dat$num_case_variants != dat_10bp$num_case_variants) {
                    print("Number case variants not equal!!!"); return()	
                }
                if (dat$num_control_variants != dat_10bp$num_control_variants) {
                    print("Number control variants not equal!!!"); return()	
                }
                if (dat$num_case_samples != dat_10bp$num_case_samples) {
                    print("Number case samples not equal!!!"); return()	
                }
                if (dat$num_control_samples != dat_10bp$num_control_samples) {
                    print("Number control samples not equal!!!"); return()	
                }
                if (sum(dat$cases_RBP > dat_10bp$cases_RBP) > 0) {
                    print(paste0("More case hits than case_10bp hits in ", output_file, "!!!"))	
                }
                if (sum(dat$controls_RBP > dat_10bp$controls_RBP) > 0) {
                    print(paste0("More control hits than control_10bp hits in ", output_file, "!!!"))	
                }
                summary <- cbind(paste0(dat$variant_type), dat$enrichment, dat$p.value, dat_10bp$enrichment, dat_10bp$p.value, dat$cases_RBP, dat_10bp$cases_RBP, dat$controls_RBP, dat_10bp$controls_RBP, dat$num_case_variants, dat$num_control_variants, dat$num_case_samples, dat$num_control_samples)
                colnames(summary) <- c("variant_type", "enrichment", "p.value", "enrichment_10bp", "p.value_10bp", "cases_RBP", "cases_10bp_RBP", "controls_RBP", "controls_10bp_RBP", "num_case_variants", "num_control_variants", "num_case_samples", "num_control_samples")
                write.csv(summary, file=output_file, row.names=FALSE)
                summary_all_expression <- rbind(summary_all_expression, cbind(rep(expression_bucket_name,3), summary))
            }
            colnames(summary_all_expression)[1] <- "expression"
            write.csv(summary_all_expression, file=paste0("rbp_burden_", disease, "_", region, "_comparison.csv"), row.names=FALSE)
        }
        
    }
    
    # Plot RBP burden analysis results for the specified disease and peak pillow/padding size. 
    # The plot argument can be "expression", "constrained", or "all" to make all possible plots.
    plot_RBP_burden_results <- function(disease, plot="all", pillow=NULL) {
        #if (is.null(pillows)) { pillow = NULL }
        
        if(!is.null(pillow)) { disease = paste0(disease, "_", pillow, "bp") }
        
        # Expression breakdown burden plots
        if(plot=="all" || plot=="expression") {
            if (grepl("CDH", disease)) { 
                expression_string_token <- "DE"
            } else if (grepl("CHD", disease)) { 
                expression_string_token <- "HE"
            } else { return(NULL) }
            
            expression_bucket_names <- paste0(c("L", "ML", "MH", "H"), expression_string_token)
            for(region in c("TSS", "TES")) {
                results <- data.frame()
                single_star <- c()
                double_star <- c()
                triple_star <- c()
                variant_counts <- c()
                for(expression_bucket_name in expression_bucket_names) {
                    filename <- full_path(RBP_OUTPUT_FOLDER, paste0("rbp_burden_", disease, "_", expression_bucket_name, "_", region, ".csv"))
                    dat <- read.csv(filename)
                    single_star <- c(single_star, dat$enrichment*(2*(dat$p.value < 0.05 & dat$p.value >= 0.01)-1)+0.05-2*(dat$enrichment==0), -2) 
                    double_star <- c(double_star, dat$enrichment*(2*(dat$p.value < 0.01 & dat$p.value >= 0.001)-1)+0.05-2*(dat$enrichment==0), -2) 
                    triple_star <- c(triple_star, dat$enrichment*(2*(dat$p.value < 0.001)-1)+0.05-2*(dat$enrichment==0), -2)	
                    variant_counts <- c(variant_counts, dat$cases_RBP, -10)
                    if(ncol(results) == 0) {
                        results <- cbind(dat$enrichment)
                    } else { results <- cbind(results, dat$enrichment) }
                }
                colnames(results) <- expression_bucket_names
                rownames(results) <- c("SNV + indel", "indels only", "SNVs only")
                
                file <- output_path(paste0("rbp_burden_", disease, "_", region, ".pdf"))
                if(!is.null(file)) { pdf(file) }
                mar=1.1*par("mar")
                mar[4] = 5
                if(!is.null(mar)) { default_mar = par("mar"); par(mar=mar) }
                
                barplot(results, beside=TRUE, main=paste0(disease, " ", region, " RBP Burden Analysis"), xlab="Expression", ylab="Enrichment", ylim=c(0, 1.4*max(results)), col=c(2,3,4), cex.lab=1.8, cex.axis=1.5, cex.names=1.5, cex.main=1.5)
                legend("topleft", legend=c(rownames(results), "count"), col=c(2,3,4,1), pch=c(15, 15, 15, 16), cex=1.2)
                mtext("p.values: * < 0.05, ** < 0.01, *** < 0.001", cex=1.25)
                
                # Add significance stars
                #text(c(2,2,2,0, 1.5, 1.5, 1.5, 0, 1,1,1), "    *", cex=1.8) # test
                text(single_star, "     *", col="orange", cex=1.4) # single star
                text(double_star, "     **", col="orange", cex=1.4) # double star
                text(triple_star, "     ***", col="orange", cex=1.4) # triple star
                
                # Plot RBP_cases counts using right axis.
                par(new = T)
                plot(seq(0.5, length(variant_counts)), variant_counts, pch=16, axes=F, xlab=NA, ylab=NA, yaxs="i", xlim=c(0, length(variant_counts)-1), ylim=c(0, 1.4*max(variant_counts)))
                axis(side=4, cex.axis=1.5)
                mtext(side=4, line=3, "Case Variant Count", cex=1.7)
                
                if(!is.null(file)) { dev.off() }
                if(!is.null(mar)) { par(mar=default_mar); dev.off() }
            }
        }
        
        # Constrained gene breakdown burden plots
        if(plot=="all" || plot=="constrained") {
            for(region in c("TSS", "TES")) {
                results <- data.frame()
                single_star <- c()
                double_star <- c()
                triple_star <- c()
                variant_counts <- c()
                for(bucket_name in c("unconstrained", "constrained")) {
                    filename <- full_path(RBP_OUTPUT_FOLDER, paste0("rbp_burden_", disease, "_", bucket_name, "_", region, ".csv"))
                    dat <- read.csv(filename)
                    single_star <- c(single_star, dat$enrichment*(2*(dat$p.value < 0.05 & dat$p.value >= 0.01)-1)+0.05-2*(dat$enrichment==0), -2) 
                    double_star <- c(double_star, dat$enrichment*(2*(dat$p.value < 0.01 & dat$p.value >= 0.001)-1)+0.05-2*(dat$enrichment==0), -2) 
                    triple_star <- c(triple_star, dat$enrichment*(2*(dat$p.value < 0.001)-1)+0.05-2*(dat$enrichment==0), -2)	
                    variant_counts <- c(variant_counts, dat$cases_RBP, -10)
                    if(ncol(results) == 0) {
                        results <- cbind(dat$enrichment)
                    } else { results <- cbind(results, dat$enrichment) }
                }
                colnames(results) <- c("pLI<=0.5", "pLI>0.5")
                rownames(results) <- c("SNV + indel", "indels only", "SNVs only")
                
                file <- output_path(paste0("rbp_burden_", disease, "_", region, "_pli_breakdown.pdf"))
                if(!is.null(file)) { pdf(file) }
                mar=1.1*par("mar")
                mar[4] = 5
                if(!is.null(mar)) { default_mar = par("mar"); par(mar=mar) }
                
                barplot(results, beside=TRUE, main=paste0(disease, " ", region, " RBP Burden Analysis"), xlab="Constraint", ylab="Enrichment", ylim=c(0, 1.4*max(results)), col=c(2,3,4), cex.lab=1.8, cex.axis=1.5, cex.names=1.5, cex.main=1.5)
                legend("topleft", legend=c(rownames(results), "count"), col=c(2,3,4,1), pch=c(15, 15, 15, 16), cex=1.2)
                mtext("p.values: * < 0.05, ** < 0.01, *** < 0.001", cex=1.25)
                
                # Add significance stars
                text(single_star, "     *", col="orange", cex=1.4) # single star
                text(double_star, "     **", col="orange", cex=1.4) # double star
                text(triple_star, "     ***", col="orange", cex=1.4) # triple star
                
                # Plot RBP_cases counts using right axis.
                par(new = T)
                plot(seq(0.5, length(variant_counts)), variant_counts, pch=16, axes=F, xlab=NA, ylab=NA, yaxs="i", xlim=c(0, length(variant_counts)-1), ylim=c(0, 1.4*max(variant_counts)))
                axis(side=4, cex.axis=1.5)
                mtext(side=4, line=3, "Case Variant Count", cex=1.7)
                
                if(!is.null(file)) { dev.off() }
                if(!is.null(mar)) { par(mar=default_mar); dev.off() }
            }
            
        }
    }
    
    # Plot RBP burden analysis results for multiple pillow values.
    plot_RBP_burden_results_multi_pillow <- function(disease) {
        if (grepl("CDH", disease)) { 
            expression_string_token <- "DE"
        } else if (grepl("CHD", disease)) { 
            expression_string_token <- "HE"
        } else { return(NULL) }
        expression_bucket_names <- paste0(c("L", "ML", "MH", "H"), expression_string_token)
        pillows <- seq(0, 100, by=10)
        variant_types <- c("SNV + indel", "indels only", "SNVs only")
        
        for(region in c("TSS", "TES")) {
            file <- output_path(paste0("rbp_burden_", disease, "_", region, "_multi_pillow.pdf"))
            if(!is.null(file)) { pdf(file, height=7, width=21); par(mfrow=c(1,3)) }
            for(variant_type_index in 1:length(variant_types)) {
                variant_type = variant_types[variant_type_index]
                results <- data.frame()
                variant_counts <- data.frame()
                single_star <- c()
                double_star <- c()
                triple_star <- c()
                for(expression_bucket_name in expression_bucket_names) {
                    enrichment_vector <- c()
                    counts_vector <- c()
                    for(pillow in pillows) {
                        filename <- full_path(RBP_OUTPUT_FOLDER, paste0("rbp_burden_", disease, "_", pillow, "bp_", expression_bucket_name, "_", region, ".csv"))
                        dat <- read.csv(filename)
                        #single_star <- c(single_star, dat$enrichment*(2*(dat$p.value < 0.05 & dat$p.value >= 0.01)-1)+0.05-2*(dat$enrichment==0), -2) 
                        #double_star <- c(double_star, dat$enrichment*(2*(dat$p.value < 0.01 & dat$p.value >= 0.001)-1)+0.05-2*(dat$enrichment==0), -2) 
                        #triple_star <- c(triple_star, dat$enrichment*(2*(dat$p.value < 0.001)-1)+0.05-2*(dat$enrichment==0), -2)	
                        enrichment_vector <- c(enrichment_vector, dat$enrichment[variant_type_index])
                        counts_vector <- c(counts_vector, dat$cases_RBP[variant_type_index])
                    }
                    if(ncol(results) == 0) {
                        results <- cbind(enrichment_vector)
                        variant_counts <- cbind(counts_vector)
                    } else { results <- cbind(results, enrichment_vector); variant_counts <- cbind(variant_counts, counts_vector) }
                }
                colnames(results) <- expression_bucket_names
                rownames(results) <- pillows
                colnames(variant_counts) <- expression_bucket_names
                rownames(variant_counts) <- pillows
                
                #mar=1.1*par("mar")
                #mar[4] = 5
                #if(!is.null(mar)) { default_mar = par("mar"); par(mar=mar) }
                
                expression_line_colors <- c("blue", "green", "orange", "red")
                for(expression_bucket_index in 1:ncol(results)) {
                    expression_bucket_name <- expression_bucket_names[expression_bucket_index]
                    expression_line_color <- expression_line_colors[expression_bucket_index]
                    if(expression_bucket_index==1) { plot(pillows, results[,expression_bucket_index], col=expression_line_color, type="l", main=paste0(disease, " ", region, " RBP Burden Analysis"), xlab="eCLIP Peak Padding (bp)", ylab="Enrichment", ylim=c(0,2.5))
                    } else { lines(pillows, results[,expression_bucket_index], col=expression_line_color, type="l") }
                }
                #lines(pillows, )
                
                legend("topleft", legend=c(rev(expression_bucket_names),"count"), col=c(rev(expression_line_colors),"grey"), pch=c(16, 16, 16, 16, 16), cex=1.2)
                mtext(paste0(variant_type, "    p.values: * < 0.05, ** < 0.01, *** < 0.001"), cex=1.25)
                
                # Add significance stars
                ##text(c(2,2,2,0, 1.5, 1.5, 1.5, 0, 1,1,1), "    *", cex=1.8) # test
                #text(single_star, "     *", col="orange", cex=1.4) # single star
                #text(double_star, "     **", col="orange", cex=1.4) # double star
                #text(triple_star, "     ***", col="orange", cex=1.4) # triple star
                
                # Plot RBP_cases counts using right axis.
                #par(new = T)
                #plot(seq(0.5, length(variant_counts)), variant_counts, pch=16, axes=F, xlab=NA, ylab=NA, yaxs="i", xlim=c(0, length(variant_counts)-1), ylim=c(0, 1.4*max(variant_counts)))
                #axis(side=4, cex.axis=1.5)
                #mtext(side=4, line=3, "Case Variant Count", cex=1.7)
                
                #if(!is.null(mar)) { par(mar=default_mar); dev.off() }
            }
            if(!is.null(file)) { dev.off() }
            par(mfrow=c(1,1)) 
        }
    }
    
    # Run RBP analyis for diseases with different pillows for RBP peaks.
    if(FALSE) {
        # CDH RBP analysis
        pillows <- seq(0,100,by=10)
        for(pillow in pillows) {
            analyze_RBPs("CDH", pillow=pillow, analysis="all")
            plot_RBP_burden_results("CDH", pillow=pillow, plot="all")
            #make_RBP_comparison_table("CDH")
        }
        plot_RBP_burden_results_multi_pillow("CDH")
        
        # CHD RBP analysis
        for(pillow in pillows) {
            analyze_RBPs("CHD", pillow=pillow, analysis="all")
            plot_RBP_burden_results(paste0("CHD_", pillow, "bp"))
        }
        plot_RBP_burden_results_multi_pillow("CHD")
        
        # New CDH batch RBP analysis
        pillows <- seq(0,100,by=10)
        for(pillow in pillows) {
            analyze_RBPs("CDH_new", pillow=pillow, analysis="all")
            plot_RBP_burden_results("CDH_new", pillow=pillow, plot="all")
        }
        plot_RBP_burden_results_multi_pillow("CDH_new")
    }
    
}

# Add RBP annotations with the specified pillow sizes (and respective RBP feature names) to WGSA annotated data frames.
annotate_with_RBP <- function(disease, controls_name=NULL, controls_only=FALSE, pillows=c(0,50), feature_names=c("RBP0","RBP")) {
    for(pillow_index in 1:length(pillows)) {
        pillow = pillows[pillow_index]
        feature_name = feature_names[pillow_index]
        rbp_annotations <- load_rbp_hits(disease, pillow)
        
        if(!controls_only) {
            rbp_results <- rbp_annotations[["cases"]]
            rbp_results <- aggregate(rbp_results, by=list(rbp_results$Chrom, rbp_results$Position), FUN=function(x) { paste0(x[!duplicated(x)], collapse=",") })[,c("Chrom", "Position", "Ref", "Alt", "sample", "RBP")]  
            colnames(rbp_results)[which(colnames(rbp_results) == "RBP")] <- feature_name
            
            dat_names <- paste0(tolower(disease), "_", sapply(c("snps", "indels", "muts"), function(x) { paste0(x, c("_wgsa", "_TSS40", "_TSS", "_TES")) }))
            for(dat_name in dat_names) {
                if(dat_name %in% ls(envir=.GlobalEnv)) {
                    print(dat_name)
                    dat <- get(dat_name)
                    dat <- merge(dat[,which(colnames(dat)!=feature_name)], rbp_results[c("Chrom", "Position", "Ref", "Alt", "sample", feature_name)], by.x=c("X.chr", "pos", "ref", "alt", "sample"), by.y=c("Chrom", "Position", "Ref", "Alt", "sample"), all.x = TRUE)
                    assign(dat_name, dat, pos=globalenv())
                }
            }
        }
        if(!is.null(controls_name)) {
            rbp_results <- rbp_annotations[["controls"]]
            rbp_results <- aggregate(rbp_results, by=list(rbp_results$Chrom, rbp_results$Position), FUN=function(x) { paste0(x[!duplicated(x)], collapse=",") })[,c("Chrom", "Position", "Ref", "Alt", "sample", "RBP")]  
            colnames(rbp_results)[which(colnames(rbp_results) == "RBP")] <- feature_name
            
            dat_names <- paste0(tolower(controls_name), "_", sapply(c("snps", "indels", "muts"), function(x) { paste0(x, c("_wgsa", "_TSS40", "_TSS", "_TES")) }))
            for(dat_name in dat_names) {
                if(dat_name %in% ls(envir=.GlobalEnv)) {
                    print(dat_name)
                    dat <- get(dat_name)
                    dat <- merge(dat[,which(colnames(dat)!=feature_name)], rbp_results[c("Chrom", "Position", "Ref", "Alt", "sample", feature_name)], by.x=c("X.chr", "pos", "ref", "alt", "sample"), by.y=c("Chrom", "Position", "Ref", "Alt", "sample"), all.x = TRUE)
                    assign(dat_name, dat, pos=globalenv())
                }
            }
        }
    }
}
annotate_with_RBP("CDH", "SSC")
annotate_with_RBP("CHD")
annotate_with_RBP("CDH2", "SSC_1088")
annotate_with_RBP("CHD2")
annotate_with_RBP("CDH_SSC518", "SSC_518", controls_only=TRUE)
annotate_with_RBP("CHDFB", "SSCFB")
annotate_with_RBP("CHDFB2", "SSCFB2")
annotate_with_RBP("regulatory")

# Find histone mark footprint sizes
# Get combined H3K36me3 GRanges
roadmap_dir = data_path("Roadmap")
histone_mark_files <- list.files(path=roadmap_dir)
combined_H3K36me3_granges <- GRanges()
for(i in 1:length(histone_mark_files)[1]) {
    print(i)
    histone_mark_file = histone_mark_files[i]
    if(!grepl("H3K36me3", histone_mark_file)) { next }
    histone_peaks <- read.csv(full_path(roadmap_dir, histone_mark_file), sep="\t", header=FALSE)
    colnames(histone_peaks)[1:3] <- c("chromosome", "start", "end")
    histone_peak_granges <- to_genomic_regions(histone_peaks)
    combined_H3K36me3_granges <- c(combined_H3K36me3_granges, histone_peak_granges)
    combined_H3K36me3_granges <- intersect(combined_H3K36me3_granges, combined_H3K36me3_granges)
}
sum(end(combined_H3K36me3_granges) - start(combined_H3K36me3_granges))
# Get combined RBP GRanges
rbp_dir = data_path("ENCODE_eCLIP_peak")
rbp_files <- list.files(path=rbp_dir)
combined_rbp_granges <- GRanges()
for(i in 1:length(rbp_files)[1]) {
    print(i)
    rbp_file = rbp_files[i]
    eclip_peaks <- read.csv(full_path(rbp_dir, rbp_file), sep="\t", header=FALSE)
    colnames(eclip_peaks)[1:3] <- c("chromosome", "start", "end")
    eclip_peak_granges <- to_genomic_regions(eclip_peaks)
    combined_rbp_granges <- c(combined_rbp_granges, eclip_peak_granges)
    combined_rbp_granges <- intersect(combined_rbp_granges, combined_rbp_granges)
}
sum(end(combined_rbp_granges) - start(combined_rbp_granges))
# Combine H3K36me3 with RBP and calculate combined footprint
rbp_H3K36me3_granges <- intersect(combined_rbp_granges, combined_H3K36me3_granges)
sum(end(rbp_H3K36me3_granges) - start(rbp_H3K36me3_granges))

CTTNB1_dat <- data.frame(t(c("3", 41236328-20000, 41301587+20000))); colnames(CTTNB1_dat) <- c("chromosome", "start", "end")
CTTNB1_Grange <- to_genomic_regions(CTTNB1_dat)
rbp_H3K36me3_cttnb1_granges <- intersect(CTTNB1_Grange, rbp_H3K36me3_granges)
sum(end(rbp_H3K36me3_cttnb1_granges) - start(rbp_H3K36me3_cttnb1_granges))

annotate_with_manual_features <- function(disease) {
    if(grepl("SSC",disease)) { relevant_tissues <- get_relevant_roadmap_eids("CHD") } else { relevant_tissues <- get_relevant_roadmap_eids(disease) }
    num_relevant_tissues = length(relevant_tissues)
    dat_names <- paste0(tolower(disease), "_", sapply(c("snps", "indels", "muts"), function(x) { paste0(x, c("_wgsa", "_TSS40", "_TSS", "_TES")) }))
    for(dat_name in dat_names) {
        if(dat_name %in% ls(envir=.GlobalEnv)) {
            print(dat_name)
            dat <- get(dat_name)
            # Annotate majority_vote histone marks
            histone_mark_types <- unique(sapply(colnames(dat)[which(unlist(lapply(strsplit(colnames(dat),"\\."), function(x) length(x) == 3)))], function(full_mark_name) strsplit(full_mark_name,"\\.")[[1]][2]))
            histone_mark_types <- histone_mark_types[!grepl("tissue_majority",histone_mark_types)]
            histone_marks_env <- new.env()
            for(histone_mark_type in histone_mark_types) {
                for(tissue in relevant_tissues) {
                    histone_mark = paste0(c(tissue, histone_mark_type, "gappedPeak"), collapse=".")
                    if(!(histone_mark %in% colnames(dat))) { next }
                    histone_marks_env[[histone_mark_type]] <- c(histone_marks_env[[histone_mark_type]], histone_mark)
                }
            }
            histone_mark_types = ls(histone_marks_env)
            num_histone_mark_types = length(histone_mark_types)
            for(i in 1:num_histone_mark_types) {
                histone_mark_type = histone_mark_types[i]
                histone_marks <- histone_marks_env[[histone_mark_type]] 
                num_histone_marks = length(histone_marks); majority_vote_threshold = min(c(3, num_histone_marks/2))
                print(paste0("Finding tissue majority vote for ",histone_mark_type," histone mark...[", i, " / ", num_histone_mark_types, "]"))
                histone_mark_tissue_majority <- rep("N", nrow(dat))
                if(num_histone_marks == 1) { histone_mark_tissue_majority[sapply(dat[,colnames(dat) %in% histone_marks], function(x) sum(x == "Y") >= majority_vote_threshold)] <- "Y"
                } else { histone_mark_tissue_majority[apply(dat[,colnames(dat) %in% histone_marks], 1, function(x) sum(x == "Y") >= majority_vote_threshold)] <- "Y" }
                #plot(density(apply(dat[,colnames(dat) %in% histone_marks], 1, function(x) sum(x == "Y"))))
                #sum(apply(dat[,colnames(dat) %in% histone_marks], 1, function(x) sum(x == "Y")) >= 5 & apply(dat[,colnames(dat) %in% histone_marks[grepl("E003|E008",histone_marks)]], 1, function(x) sum(x == "Y")) >= 1)/sum(apply(dat[,colnames(dat) %in% histone_marks], 1, function(x) sum(x == "Y")) >= 5)
                histone_mark_tissue_majority_colname = paste0(histone_mark_type,".tissue_majority.gappedPeak")
                if(histone_mark_tissue_majority_colname %in% colnames(dat)) {
                    dat[,histone_mark_tissue_majority_colname] <- histone_mark_tissue_majority
                } else {
                    dat <- cbind(dat, histone_mark_tissue_majority)
                    colnames(dat)[ncol(dat)] <- histone_mark_tissue_majority_colname
                }
            }
            
            # Annotate with Felix's heart enhancers.
            print("Annotating H3K27ac Dickel heart enhancer (data from Felix)")
            H3K27ac_heart_enhancer_regions <- read.table(data_path("enh_heart_Dickel_2015_sort.bed"))
            colnames(H3K27ac_heart_enhancer_regions) <- c("chromosome", "start", "end", "enhancer_hit_score")
            H3K27ac_heart_enhancer_granges <- to_genomic_regions(H3K27ac_heart_enhancer_regions, labels=H3K27ac_heart_enhancer_regions$enhancer_hit_score)
            H3K27ac_Dickel_heart_enhancer <- rep("N", nrow(dat))
            dat_granges <- to_genomic_regions(dat, chr_colname="X.chr", start_colname="pos")
            H3K27ac_Dickel_heart_enhancer[unique(queryHits(findOverlaps(dat_granges, H3K27ac_heart_enhancer_granges)))] <- "Y"
            if("H3K27ac_Dickel_heart_enhancer" %in% colnames(dat)) {
                dat[,"H3K27ac_Dickel_heart_enhancer"] <- H3K27ac_Dickel_heart_enhancer
            } else { dat <- cbind(dat, H3K27ac_Dickel_heart_enhancer) }
            
            assign(dat_name, dat, pos=globalenv())
        }
    }
}
annotate_with_manual_features("CHDFB_combined")
annotate_with_manual_features("SSCFB_combined")

filter_dataset <- function(dataset, version="hg19", remove_synonymous=TRUE, remove_missense=TRUE, preview=FALSE) {
    all_filters_granges <- get(paste0("all_filters_granges_",version))
    
    dat_names <- paste0(tolower(dataset), "_", sapply(c("snps", "indels", "muts"), function(x) { paste0(x, c("_wgsa", "_TSS40", "_TSS", "_TES")) }))
    for(dat_name in dat_names) {
        if(dat_name %in% ls(envir=.GlobalEnv)) {
            print(dat_name)
            dat <- get(dat_name)
            num_variants <- nrow(dat)
            if(remove_missense || remove_synonymous) { 
                dat_annovar_ensembl <- gsub("\\(.*", "", gsub("^.*:", "", dat$ANNOVAR_ensembl_summary)); dat$ANNOVAR_ensembl_summary <- dat_annovar_ensembl
                if(remove_synonymous) { coding_to_remove <- c("nonsynonymous", "synonymous") } else { coding_to_remove <- c("nonsynonymous") }
                dat <- dat[!(dat_annovar_ensembl %in% coding_to_remove),]
                dat_annovar_ensembl <- table(dat_annovar_ensembl[dat_annovar_ensembl %in% coding_to_remove])
                print(paste0(num_variants-nrow(dat)," variants removed (", paste(dat_annovar_ensembl, names(dat_annovar_ensembl), collapse=", "), ") due to ANNOVAR ensembl coding region."))
                num_variants <- nrow(dat)
            }
            if ("X.chr" %in% colnames(dat)) { chr_colname = "X.chr"; pos_colname = "pos" } else { chr_colname = "Chrom"; pos_colname = "Position" }
            dat_granges <- to_genomic_regions(dat, chr_colname=chr_colname, start_colname=pos_colname, end_colname=pos_colname)
            dat_indices_to_filter_out <- c(unique(queryHits(findOverlaps(dat_granges, all_filters_granges))))
            if(length(dat_indices_to_filter_out) > 0) { dat <- dat[-c(dat_indices_to_filter_out),]; dat_granges <- dat_granges[-c(dat_indices_to_filter_out)] }
            print(paste0(num_variants-length(dat_granges)," variants removed due to region filters."))
            if(preview) {
                print(dat_granges)
            } else {
                assign(dat_name, dat, pos=globalenv())
            }
        }
    }
}
filter_dataset("CDH")
filter_dataset("CDH2")
filter_dataset("CHD")
filter_dataset("CHD2")
filter_dataset("SSC_518")
filter_dataset("SSC_1088")
filter_dataset("CHDFB")
filter_dataset("CHDFB2")
filter_dataset("SSCFB")
filter_dataset("SSCFB2")

filter_dataset("CHDFB_combined", remove_synonymous=FALSE)
filter_dataset("SSCFB_combined", remove_synonymous=FALSE)

# Combine all cases (CDH + CHD) together for burden analysis purposes.
combine_cdh_with_chd <- function() {
    cdh_dat_names <- paste0("cdh_combined_", sapply(c("snps", "indels", "muts"), function(x) { paste0(x, c("_wgsa", "_TSS40", "_TSS", "_TES")) }))
    chd_dat_names <- gsub("cdh", "chd", cdh_dat_names)
    cdh_chd_dat_names <- gsub("cdh_combined", "cdh_chd", cdh_dat_names)
    for(i in 1:length(cdh_dat_names)) {
        cdh_chd_dat_name <- cdh_chd_dat_names[i]
        print(cdh_chd_dat_name)
        cdh_dat <- get(cdh_dat_names[i])
        chd_dat <- get(chd_dat_names[i])
        assign(cdh_chd_dat_name, combine_datasets(cdh_dat, chd_dat), pos=globalenv())
    }
    assign("CDH_CHD_SAMPLE_COUNT", length(unique(cdh_chd_muts_wgsa$sample)), pos=globalenv())
}
combine_cdh_with_chd()

# More fine-grain (ie: integrating important histone marks) RBP burden checks.
# SORT OF LEGACY CODE!
if (TRUE) {
    cdh_muts_TSS_rbp <- merge(cdh_muts_TSS, rbp_results_cases, by.x=c("X.chr", "pos", "sample"), by.y=c("Chrom", "Position", "sample"))
    cdh_muts_TES_rbp <- merge(cdh_muts_TES, rbp_results_cases, by.x=c("X.chr", "pos", "sample"), by.y=c("Chrom", "Position", "sample"))
    cdh_muts_TSS_rbp <- aggregate(cdh_muts_TSS_rbp$RBP, by=cdh_muts_TSS_rbp[,c("TS_gene", "X.chr", "pos", "ref", "alt", "sample", colnames(cdh_muts_TSS_rbp)[grepl("gappedPeak", colnames(cdh_muts_TSS_rbp))])], paste0, collapse=","); colnames(cdh_muts_TSS_rbp)[ncol(cdh_muts_TSS_rbp)] <- "RBP"
    cdh_muts_TES_rbp <- aggregate(cdh_muts_TES_rbp$RBP, by=cdh_muts_TES_rbp[,c("TS_gene", "X.chr", "pos", "ref", "alt", "sample", colnames(cdh_muts_TES_rbp)[grepl("gappedPeak", colnames(cdh_muts_TES_rbp))])], paste0, collapse=","); colnames(cdh_muts_TES_rbp)[ncol(cdh_muts_TES_rbp)] <- "RBP"
    # H3K36me3 in H1 stem cells of TES of HDE genes
    rbp_H3K36me3_H1_tes_cases  <- cdh_muts_TES_rbp[cdh_muts_TES_rbp$E003.H3K36me3.gappedPeak == "Y",c("TS_gene", "X.chr", "pos", "ref", "alt", "sample", "RBP")]
    # H3K36me3 in H9 stem cells of TES of HDE genes
    rbp_H3K36me3_H9_tes_cases <- cdh_muts_TES_rbp[cdh_muts_TES_rbp$E008.H3K36me3.gappedPeak == "Y",c("TS_gene", "X.chr", "pos", "ref", "alt", "sample", "RBP")]
    # H3K18ac in H1 stem cells of TSS of HDE genes
    rbp_H3K18ac_H1_tss_cases <- cdh_muts_TSS_rbp[cdh_muts_TSS_rbp$E003.H3K18ac.gappedPeak == "Y",c("TS_gene", "X.chr", "pos", "ref", "alt", "sample", "RBP")]
    # all TES of HDE genes
    rbp_all_tes_cases <- cdh_muts_TES_rbp[,c("TS_gene", "X.chr", "pos", "ref", "alt", "sample", "RBP")]
    # all TSS of HDE genes
    rbp_all_tss_cases <- cdh_muts_TSS_rbp[,c("TS_gene", "X.chr", "pos", "ref", "alt", "sample", "RBP")]
    
    ssc_muts_TSS_rbp <- merge(ssc_muts_TSS, rbp_results_controls, by.x=c("X.chr", "pos"), by.y=c("Chrom", "Position")); colnames(ssc_muts_TSS_rbp)[which(colnames(ssc_muts_TSS_rbp)=="sample.x")] <- "sample"
    ssc_muts_TES_rbp <- merge(ssc_muts_TES, rbp_results_controls, by.x=c("X.chr", "pos"), by.y=c("Chrom", "Position")); colnames(ssc_muts_TES_rbp)[which(colnames(ssc_muts_TES_rbp)=="sample.x")] <- "sample"
    ssc_muts_TSS_rbp <- aggregate(ssc_muts_TSS_rbp$RBP, by=ssc_muts_TSS_rbp[,c("TS_gene", "X.chr", "pos", "ref", "alt", "sample", colnames(ssc_muts_TSS_rbp)[grepl("gappedPeak", colnames(ssc_muts_TSS_rbp))])], paste0, collapse=","); colnames(ssc_muts_TSS_rbp)[ncol(ssc_muts_TSS_rbp)] <- "RBP"
    ssc_muts_TES_rbp <- aggregate(ssc_muts_TES_rbp$RBP, by=ssc_muts_TES_rbp[,c("TS_gene", "X.chr", "pos", "ref", "alt", "sample", colnames(ssc_muts_TES_rbp)[grepl("gappedPeak", colnames(ssc_muts_TES_rbp))])], paste0, collapse=","); colnames(ssc_muts_TES_rbp)[ncol(ssc_muts_TES_rbp)] <- "RBP"
    # SSC H3K36me3 in H1 stem cells of TES of HDE genes
    rbp_H3K36me3_H1_tes_controls <- ssc_muts_TES_rbp[ssc_muts_TES_rbp$E003.H3K36me3.gappedPeak == "Y",c("TS_gene", "X.chr", "pos", "ref", "alt", "sample", "RBP")]
    # SSC H3K36me3 in H9 stem cells of TES of HDE genes
    rbp_H3K36me3_H9_tes_controls <- ssc_muts_TES_rbp[ssc_muts_TES_rbp$E008.H3K36me3.gappedPeak == "Y",c("TS_gene", "X.chr", "pos", "ref", "alt", "sample", "RBP")]
    # SSC H3K18ac in H1 stem cells of TSS of HDE genes
    rbp_H3K18ac_H1_tss_controls <- ssc_muts_TSS_rbp[ssc_muts_TSS_rbp$E003.H3K18ac.gappedPeak == "Y",c("TS_gene", "X.chr", "pos", "ref", "alt", "sample", "RBP")]
    #write.csv(file="out.csv", ssc_muts_TSS_rbp[ssc_muts_TSS_rbp$E003.H3K18ac.gappedPeak == "Y",c("TS_gene", "X.chr", "pos", "ref", "alt", "sample", "RBP")], row.names=FALSE)
    # all TES of HDE genes
    rbp_all_tes_controls <- ssc_muts_TES_rbp[,c("TS_gene", "X.chr", "pos", "ref", "alt", "sample", "RBP")]
    # all TSS of HDE genes
    rbp_all_tss_controls <- ssc_muts_TSS_rbp[,c("TS_gene", "X.chr", "pos", "ref", "alt", "sample", "RBP")]
    
    rbp_genes_H3K36me3_H1_tes <- sort(unique(unlist(strsplit(paste(rbp_H3K36me3_H1_tes_cases$TS_gene), ",")))[! unique(unlist(strsplit(paste(rbp_H3K36me3_H1_tes_cases$TS_gene), ","))) %in% unlist(strsplit(paste(rbp_H3K36me3_H1_tes_controls$TS_gene), ","))])
    rbp_genes_H3K36me3_H9_tes <- sort(unique(unlist(strsplit(paste(rbp_H3K36me3_H9_tes_cases$TS_gene), ",")))[! unique(unlist(strsplit(paste(rbp_H3K36me3_H9_tes_cases$TS_gene), ","))) %in% unlist(strsplit(paste(rbp_H3K36me3_H9_tes_controls$TS_gene), ","))])
    rbp_genes_H3K18ac_H1_tss <- sort(unique(unlist(strsplit(paste(rbp_H3K18ac_H1_tss_cases$TS_gene), ",")))[! unique(unlist(strsplit(paste(rbp_H3K18ac_H1_tss_cases$TS_gene), ","))) %in% unlist(strsplit(paste(rbp_H3K18ac_H1_tss_controls$TS_gene), ","))])
    rbp_genes_all_tes <- sort(unique(unlist(strsplit(paste(rbp_all_tes_cases$TS_gene), ",")))[! unique(unlist(strsplit(paste(rbp_all_tes_cases$TS_gene), ","))) %in% unlist(strsplit(paste(rbp_all_tes_controls$TS_gene), ","))])
    rbp_genes_all_tss <- sort(unique(unlist(strsplit(paste(rbp_all_tss_cases$TS_gene), ",")))[! unique(unlist(strsplit(paste(rbp_all_tss_cases$TS_gene), ","))) %in% unlist(strsplit(paste(rbp_all_tss_controls$TS_gene), ","))])
    
    RBP_gene_results <- rbind(cbind("3'UTR of HDE genes", paste0(rbp_genes_all_tes, collapse=",")),
                              cbind("3'UTR of HDE genes with H3K36me3 in H1 stem cells", paste0(rbp_genes_H3K36me3_H1_tes, collapse=",")),
                              cbind("3'UTR of HDE genes with H3K36me3 in H9 stem cells", paste0(rbp_genes_H3K36me3_H9_tes, collapse=",")),
                              cbind("TSS of HDE genes", paste0(rbp_genes_all_tss, collapse=",")),
                              cbind("TSS of HDE genes with H3K18ac in H1 stem cells", "-"))
    colnames(RBP_gene_results) <- c("region", "genes")
    write.csv(file="rbp_gene_results.csv", RBP_gene_results, row.names=FALSE)
    
    
    cstf2_dat <- read.csv(paste0(rbp_dir, "\\HepG2.CSTF2.R2.tag.uniq.peak.sig.bed"), sep="\t", header=FALSE)
    colnames(cstf2_dat)[1:3] <- c("chromosome", "start", "end")
    cstf2 <- genomic_coordinates_to_strings(cstf2_dat[,1], cstf2_dat[,2], cstf2_dat[,3])
    duplicate_cstf2 <- duplicated(cstf2)
    cstf2_granges <- to_genomic_regions(cstf2_dat[!duplicate_cstf2,], labels=cstf2[!duplicate_cstf2])
    hits <- unique(data.frame(olRanges(cases_granges, cstf2_granges)))
    
    rac1_dat <- read.csv(paste0(rbp_dir, "\\HepG2.CSTF2.R2.tag.uniq.peak.sig.bed"), sep="\t", header=FALSE)
    colnames(cstf2_dat)[1:3] <- c("chromosome", "start", "end")
    cstf2 <- genomic_coordinates_to_strings(cstf2_dat[,1], cstf2_dat[,2], cstf2_dat[,3])
    duplicate_cstf2 <- duplicated(cstf2)
    cstf2_granges <- to_genomic_regions(cstf2_dat[!duplicate_cstf2,], labels=cstf2[!duplicate_cstf2])
    hits <- unique(data.frame(olRanges(cases_granges, cstf2_granges)))
    
    #genes <- genomic_coordinates_to_strings(genebody[,1], min(genebody[,2], genebody[,3]), max(genebody[,2], genebody[,3]))
    TES_regions <- get_TS_regions(genebody, sites="TES", left=5000, right=20000)
    TES_granges <- to_genomic_regions(TES_regions, labels=genomic_coordinates_to_strings(TES_regions[,1], TES_regions[,2], TES_regions[,3]))
    hits <- unique(data.frame(olRanges(TES_granges, cstf2_granges)))
    cstf2_genes <- TES_regions[genomic_coordinates_to_strings(TES_regions$chromosome, TES_regions$start, TES_regions$end) %in% genomic_coordinates_to_strings(hits$seqnames, hits$start, hits$end),]
    write.csv(file="cstf2_genes.csv", cstf2_genes, row.names=FALSE)
    
    # Annotate CDH2 with histone marks:
    if (TRUE) {
        # E001.H3K36me3.gappedPeak
        E001_H3K36me3_regions <- read.table(data_path("Roadmap/E001-H3K36me3.gappedPeak"))[,1:3]
        colnames(E001_H3K36me3_regions) <- c("chromosome", "start", "end")
        E001_H3K36me3_granges <- to_genomic_regions(E001_H3K36me3_regions)
        #hits_variants <- paste0(apply(cdh_new[unique(queryHits(findOverlaps(cdh_new_granges, E001_H3K36me3_granges))),1:2], 1, paste, collapse=":"))
        E001.H3K36me3.gappedPeak <- rep("N", nrow(cdh_new))
        E001.H3K36me3.gappedPeak[unique(queryHits(findOverlaps(cdh_new_granges, E001_H3K36me3_granges)))] <- "Y"
        cdh_new <- cbind(cdh_new, E001.H3K36me3.gappedPeak)
        # E003.H3K36me3.gappedPeak
        E003_H3K36me3_regions <- read.table(data_path("Roadmap/E003-H3K36me3.gappedPeak"))[,1:3]
        colnames(E003_H3K36me3_regions) <- c("chromosome", "start", "end")
        E003_H3K36me3_granges <- to_genomic_regions(E003_H3K36me3_regions)
        E003.H3K36me3.gappedPeak <- rep("N", nrow(cdh_new))
        E003.H3K36me3.gappedPeak[unique(queryHits(findOverlaps(cdh_new_granges, E003_H3K36me3_granges)))] <- "Y"
        cdh_new <- cbind(cdh_new, E003.H3K36me3.gappedPeak)
        # E008.H3K36me3.gappedPeak
        E008_H3K36me3_regions <- read.table(data_path("Roadmap/E008-H3K36me3.gappedPeak"))[,1:3]
        colnames(E008_H3K36me3_regions) <- c("chromosome", "start", "end")
        E008_H3K36me3_granges <- to_genomic_regions(E008_H3K36me3_regions)
        E008.H3K36me3.gappedPeak <- rep("N", nrow(cdh_new))
        E008.H3K36me3.gappedPeak[unique(queryHits(findOverlaps(cdh_new_granges, E008_H3K36me3_granges)))] <- "Y"
        cdh_new <- cbind(cdh_new, E008.H3K36me3.gappedPeak)
        # E015.H3K36me3.gappedPeak
        E015_H3K36me3_regions <- read.table(data_path("Roadmap/E015-H3K36me3.gappedPeak"))[,1:3]
        colnames(E015_H3K36me3_regions) <- c("chromosome", "start", "end")
        E015_H3K36me3_granges <- to_genomic_regions(E015_H3K36me3_regions)
        E015.H3K36me3.gappedPeak <- rep("N", nrow(cdh_new))
        E015.H3K36me3.gappedPeak[unique(queryHits(findOverlaps(cdh_new_granges, E015_H3K36me3_granges)))] <- "Y"
        cdh_new <- cbind(cdh_new, E015.H3K36me3.gappedPeak)
        # E083.H3K4me1.gappedPeak
        E083_H3K4me1_regions <- read.table(data_path("Roadmap/E083-H3K4me1.gappedPeak"))[,1:3]
        colnames(E083_H3K4me1_regions) <- c("chromosome", "start", "end")
        E083_H3K4me1_granges <- to_genomic_regions(E083_H3K4me1_regions)
        E083.H3K4me1.gappedPeak <- rep("N", nrow(cdh_new))
        E083.H3K4me1.gappedPeak[unique(queryHits(findOverlaps(cdh_new_granges, E083_H3K4me1_granges)))] <- "Y"
        cdh_new <- cbind(cdh_new, E083.H3K4me1.gappedPeak)
        
    }
}

# Annotate with Felix's H3K27ac/p300 enhancer data:
H3K27ac_heart_enhancer_regions <- read.table(data_path("enh_heart_Dickel_2015_sort.bed"))
colnames(H3K27ac_heart_enhancer_regions) <- c("chromosome", "start", "end", "enhancer_hit_score")
H3K27ac_heart_enhancer_granges <- to_genomic_regions(H3K27ac_heart_enhancer_regions, labels=H3K27ac_heart_enhancer_regions$enhancer_hit_score)

chd_in_mouse_candidates <- sapply(strsplit(paste0(chd_muts_wgsa$SnpEff_ensembl_Gene_name), "\\|"), function(x) { sum(paste0(unique(x)) %in% candidate_CHD_genes) > 0 })
sscfb_in_mouse_candidates <- sapply(strsplit(paste0(sscfb_muts_wgsa$SnpEff_ensembl_Gene_name), "\\|"), function(x) { sum(paste0(unique(x)) %in% candidate_CHD_genes) > 0 })

chd_muts_felix_TSS_H3K27ac_heart_enhancer_Y <- sum(includes_at_least_one(is_TS(chd_granges[is_TS(chd_granges, TS_granges=H3K27ac_heart_enhancer_granges)], TS_granges=Felix_TSS_granges, return_genes=TRUE), candidate_CHD_genes))
chd_muts_felix_TSS_H3K27ac_heart_enhancer_Y
sscfb_felix_muts_TSS_H3K27ac_heart_enhancer_Y <- sum(includes_at_least_one(is_TS(sscfb_granges[is_TS(sscfb_granges, TS_granges=H3K27ac_heart_enhancer_granges)], TS_granges=Felix_TSS_granges, return_genes=TRUE), candidate_CHD_genes))
sscfb_felix_muts_TSS_H3K27ac_heart_enhancer_Y
chd_muts_TSS_H3K27ac_heart_enhancer_Y <- sum(includes_at_least_one(is_TS(chd_granges[is_TS(chd_granges, TS_granges=H3K27ac_heart_enhancer_granges)], sites="TSS", return_genes=TRUE), candidate_CHD_genes))
chd_muts_TSS_H3K27ac_heart_enhancer_Y
sscfb_muts_TSS_H3K27ac_heart_enhancer_Y <- sum(includes_at_least_one(is_TS(sscfb_granges[is_TS(sscfb_granges, TS_granges=H3K27ac_heart_enhancer_granges)], sites="TSS", return_genes=TRUE), candidate_CHD_genes))
sscfb_muts_TSS_H3K27ac_heart_enhancer_Y
chd_muts_TES_H3K27ac_heart_enhancer_Y <- sum(includes_at_least_one(is_TS(chd_granges[is_TS(chd_granges, TS_granges=H3K27ac_heart_enhancer_granges)], sites="TES", return_genes=TRUE), candidate_CHD_genes))
chd_muts_TES_H3K27ac_heart_enhancer_Y
sscfb_muts_TES_H3K27ac_heart_enhancer_Y <- sum(includes_at_least_one(is_TS(sscfb_granges[is_TS(sscfb_granges, TS_granges=H3K27ac_heart_enhancer_granges)], sites="TES", return_genes=TRUE), candidate_CHD_genes))
sscfb_muts_TES_H3K27ac_heart_enhancer_Y


test_result <- binomial_test(chd_muts_felix_TSS_H3K27ac_heart_enhancer_Y, sscfb_felix_muts_TSS_H3K27ac_heart_enhancer_Y, sample_count=CHD_SAMPLE_COUNT, control_sample_count=SSCFB_SAMPLE_COUNT, alternative=c("greater"))
print(paste0("Felix TSS result: ", test_result[["estimate"]], ", ", test_result[["p.value"]]))
test_result <- binomial_test(chd_muts_TSS_H3K27ac_heart_enhancer_Y, sscfb_muts_TSS_H3K27ac_heart_enhancer_Y, sample_count=CHD_SAMPLE_COUNT, control_sample_count=SSCFB_SAMPLE_COUNT, alternative=c("greater"))
print(paste0("TSS result: ", test_result[["estimate"]], ", ", test_result[["p.value"]]))
test_result <- binomial_test(chd_muts_TES_H3K27ac_heart_enhancer_Y, sscfb_muts_TES_H3K27ac_heart_enhancer_Y, sample_count=CHD_SAMPLE_COUNT, control_sample_count=SSCFB_SAMPLE_COUNT, alternative=c("greater"))
print(paste0("TES result: ", test_result[["estimate"]], ", ", test_result[["p.value"]]))

#################################################################################################################
# Annotate HDE and HHE TSS/TES regions with expected mutation counts, from trimer background mutation rates. 
# Note that this is a very time-consuming step with current version of get_expected_mut_counts method.
#################################################################################################################
ANNOTATE_EXPECTED_MUT_COUNTS = FALSE
if (ANNOTATE_EXPECTED_MUT_COUNTS) {
    TES_regions_exp <- get_expected_mut_counts(TES_regions, precision=1000)
    TSS_regions_exp <- get_expected_mut_counts(TSS_regions, precision=1000)
}


rbp_results_controls_10bp_subset <- rbp_results_controls_10bp[is_TS(rbp_results_controls_10bp$Chrom, rbp_results_controls_10bp$Position, disease="SSC", TS_regions=HHE_TSS_regions, clear_hash=TRUE),]

rbp_results_controls_subset <- rbp_results_controls_subset[order(rbp_results_controls_subset$Chrom, rbp_results_controls_subset$Position),]
rbp_results_controls_10bp_subset <- rbp_results_controls_10bp_subset[order(rbp_results_controls_10bp_subset$Chrom, rbp_results_controls_10bp_subset$Position),]
rbp_results_controls_subset[is_indel(rbp_results_controls_subset$Ref, rbp_results_controls_subset$Alt),]	
rbp_results_controls_subset_pasted <- sapply(1:nrow(rbp_results_controls_subset), function(i) { paste0(rbp_results_controls_subset[i,c(1:5)], collapse="_") } ) 
rbp_results_controls_10bp_subset_pasted <- sapply(1:nrow(rbp_results_controls_10bp_subset), function(i) { paste0(rbp_results_controls_10bp_subset[i,c(1:5)], collapse="_") } ) 
rbp_results_controls_subset_pasted[!(rbp_results_controls_subset_pasted %in% rbp_results_controls_10bp_subset_pasted)] 
rbp_results_controls_subset[1:10,]
rbp_results_controls_10bp_subset[1:10,]
sum(rbp_results_controls_subset_pasted %in% rbp_results_controls_10bp_subset_pasted)
nrow(rbp_results_controls_subset)
nrow(rbp_results_controls_10bp_subset)


TS_regions_none <- get_TS_regions(genebody, sites="TES", left=0, right=0)[1:1000,]
TS_regions_left <- get_TS_regions(genebody, sites="TES", left=5000, right=0)[1:1000,]
TS_regions_right <- get_TS_regions(genebody, sites="TES", left=0, right=20000)[1:1000,]
TS_regions_all <- get_TS_regions(genebody, sites="TES", left=5000, right=20000)[1:1000,]
cases_all <- cases[is_TS(cases$Chrom, cases$Position, disease="CDH", sites="TES", TS_regions=TS_regions_all, print_counts=FALSE),]
cases_left <- cases[is_TS(cases$Chrom, cases$Position, disease="CDH", sites="TES", TS_regions=TS_regions_left, print_counts=FALSE),]
cases_right <- cases[is_TS(cases$Chrom, cases$Position, disease="CDH", sites="TES", TS_regions=TS_regions_right, print_counts=FALSE),]
cases_none <- cases[is_TS(cases$Chrom, cases$Position, disease="CDH", sites="TES", TS_regions=TS_regions_none, print_counts=FALSE),]

apply(cases_left[,1:2], 1, paste0, collapse=",") %in% apply(cases_all[,1:2], 1, paste0, collapse=",")
apply(cases_right[,1:2], 1, paste0, collapse=",") %in% apply(cases_all[,1:2], 1, paste0, collapse=",")
cases_left[apply(cases_left[,1:2], 1, paste0, collapse=",") %in% apply(cases_right[,1:2], 1, paste0, collapse=","),1:2]
apply(cases_right[,1:2], 1, paste0, collapse=",") %in% apply(cases_left[,1:2], 1, paste0, collapse=",")
c(apply(cases_left[,1:2], 1, paste0, collapse=","), apply(cases_right[,1:2], 1, paste0, collapse=",")) %in% apply(cases_all[,1:2], 1, paste0, collapse=",")


genebody[1:2,]
TS_regions_none[1:2,]
TS_regions_left[1:2,]
TS_regions_right[1:2,]	
TS_regions_all[1:2,]


test <- is_TS(cdh$Chrom, cdh$Position, disease="CDH", TS_regions=MLDE_TSS_regions, return_genes=TRUE)

#################################################################################################################
# Compare raw nc burden for the specified sites ("TSS" or TES"), at different offset levels.
#################################################################################################################
test_TS_burden <- function(cases, controls, genebody, disease, sites, left, right, by=1000, p_val_bin_size=0.01, filename_string) {
    num_control_samples = length(unique(controls$sample)); num_case_samples = length(unique(cases$sample)) 
    snv_indel_defs <- new.env(); snv_indel_defs[["snv"]] <- "snv"; snv_indel_defs[["indel"]] <- "indel"; snv_indel_defs[["snv+indel"]] <- c("snv", "indel")
    offsets <- seq(-left, right, by=by); offsets <- offsets[offsets!=0]
    odds_ratios <- c()
    p_values <- c()
    dat <- data.frame(stringsAsFactors=FALSE)
    for(offset in offsets) {
        if(offset == 0) { odds_ratios <- c(odds_ratios, 0); p_values <- c(p_values, 1); next }
        print(paste0("Processing ", sites, " Definition: [", min(offset, 0), ", ", max(offset, 0), "]"))
        TS_regions_curr <- get_TS_regions(genebody, sites=sites, left=-min(offset, 0), right=max(offset, 0))
        TS_granges_curr <- to_genomic_regions(TS_regions_curr)
        cases_TS <- cases[is_TS(cases_granges, TS_granges=TS_granges_curr),]
        controls_TS <- controls[is_TS(controls_granges, TS_granges=TS_granges_curr),]
        
        for(snv_indel_def_index in ls(snv_indel_defs)) {
            cases_selected <- cases_HE_TS[cases_HE_TS$snv_indel %in% snv_indel_defs[[snv_indel_def_index]],]
            controls_selected <- controls_HE_TS[controls_HE_TS$snv_indel %in% snv_indel_defs[[snv_indel_def_index]],]
            cases_Y <- nrow(cases_selected)
            cases_N <- sum(cases$snv_indel %in% snv_indel_defs[[snv_indel_def_index]]) - cases_Y
            controls_Y <- nrow(controls_selected)
            controls_N <- sum(controls$snv_indel %in% snv_indel_defs[[snv_indel_def_index]]) - controls_Y
            mat <- matrix(c(cases_Y, controls_Y, cases_N, controls_N), nrow=2, ncol=2)
            ft <- fisher.test(mat, alternative="greater")
            print(mat)
            print(ft)
            
            # Binomial test: Null is binom(m, p), where m = m1 + m0 (m0 is the number of such variants in controls),  p = n1/(n1+n0),   
            # n1 is the number of cases in total, n0 is the number of controls in total
            m0 = controls_Y
            m1 = cases_Y
            n0 = num_control_samples	
            n1 = num_case_samples
            binom_size = m1+m0
            binom_prob = n1/(n1+n0)
            binom_enrichment = m1/(binom_size*binom_prob)
            #binom.test(m1, binom_size, p = binom_prob, alternative = "greater", conf.level = 0.95)
            binom_p.value = pbinom(q=m1-1, size=binom_size, prob=binom_prob, lower.tail=FALSE)	
            
            odds_ratios <- c(odds_ratios, as.numeric(paste(ft$estimate)))
            p_values <- c(p_values, as.numeric(paste(ft$p.value)))
            new_dat_row <- paste(c(paste0("[", min(offset, 0), ", ", max(offset, 0), "]"), snv_indel_def_index, cases_Y, cases_N, controls_Y, controls_N, binom_enrichment, binom_p.value))
            #new_dat_row <- paste(c(paste0("[", min(offset, 0), ", ", max(offset, 0), "]"), snv_indel_def_index, cases_Y, cases_N, controls_Y, controls_N, ft$estimate, ft$p.value))
            dat <- rbind(data.frame(dat, stringsAsFactors=FALSE), t(data.frame(c(new_dat_row))))
        }
    }
    # Do one more test for our TES/TSS definition.
    if (sites == "TSS") {
        left = 20000
        right = 20000
    } else { # sites == "TES"
        left = 5000
        right = 20000
    }
    print(paste0("Processing ", sites, " Definition: [", -left, ", ", right, "]"))
    TS_regions <- get_TS_regions(genebody, sites=sites, left=left, right=right)
    cases_HE_TS <- cases[is_TS(cases$Chrom, cases$Position, disease=disease, sites=sites, TS_regions=TS_regions, print_counts=FALSE),]
    controls_HE_TS <- controls[is_TS(controls$Chrom, controls$Position, disease=disease, sites=sites, TS_regions=TS_regions, print_counts=FALSE),]
    
    for(snv_indel_def_index in ls(snv_indel_defs)) {
        cases_selected <- cases_HE_TS[cases_HE_TS$snv_indel %in% snv_indel_defs[[snv_indel_def_index]],]
        controls_selected <- controls_HE_TS[controls_HE_TS$snv_indel %in% snv_indel_defs[[snv_indel_def_index]],]
        cases_Y <- nrow(cases_selected)
        cases_N <- sum(cases$snv_indel %in% snv_indel_defs[[snv_indel_def_index]]) - cases_Y
        controls_Y <- nrow(controls_selected)
        controls_N <- sum(controls$snv_indel %in% snv_indel_defs[[snv_indel_def_index]]) - controls_Y
        mat <- matrix(c(cases_Y, controls_Y, cases_N, controls_N), nrow=2, ncol=2)
        ft <- fisher.test(mat, alternative="greater")
        print(mat)
        print(ft)
        
        # Binomial test: Null is binom(m, p), where m = m1 + m0 (m0 is the number of such variants in controls),  p = n1/(n1+n0),   
        # n1 is the number of cases in total, n0 is the number of controls in total
        m0 = controls_Y
        m1 = cases_Y
        n0 = num_control_samples	
        n1 = num_case_samples	
        binom_size = m1+m0
        binom_prob = 195/(195+438) ######### n1/(n1+n0)
        binom_enrichment = m1/(binom_size*binom_prob)
        binom_p.value = pbinom(q=m1-1, size=binom_size, prob=binom_prob, lower.tail=FALSE)	
        
        new_dat_row <- paste0(c(paste0("[", -left, ", ", right, "]"), snv_indel_def_index, cases_Y, cases_N, controls_Y, controls_N, binom_enrichment, binom_p.value))
        #new_dat_row <- paste0(c(paste0("[", -left, ", ", right, "]"), snv_indel_def_index, cases_Y, cases_N, controls_Y, controls_N, ft$estimate, ft$p.value))
        dat <- rbind(data.frame(dat, stringsAsFactors=FALSE), t(data.frame(c(new_dat_row))))
    }
    colnames(dat) <- c("region", "mut_type", "cases_Y", "cases_N", "controls_Y", "controls_N", "enrichment", "p.value")
    dat <- dat[order(dat$mut_type),] # order by mutation type
    
    # Print to file
    write.csv(dat, file=paste0(filename_string, ".csv"), row.names=FALSE)
    
    # Plot results
    #num_p_val_bins <- round(1/p_val_bin_size)
    #colfunc <- colorRampPalette(c("royalblue", "red"))
    #colors <- colfunc(num_p_val_bins)
    #p_value_colors <- data.frame(seq(1,p_val_bin_size,by=-p_val_bin_size), colfunc(num_p_val_bins))
    #colnames(p_value_colors) <- c("p_val", "col")
    #rounded_p_values <- ceiling(p_values*num_p_val_bins)/num_p_val_bins
    #rounded_p_values_indices <- sapply(rounded_p_values, function(x) which(abs(p_value_colors$p_val - x) < p_val_bin_size/100))
    #plot(offsets, odds_ratios, main=paste0(sites, " nc Enrichment"), xlab="Offset", ylab="Enrichment", col=paste(p_value_colors$col[rounded_p_values_indices]), pch=19, cex=2)
    #dev.copy2pdf(file=paste0(filename_string, "_plot.pdf"))
    
    return(dat)
}
# CDH TS burden analysis
blah <- test_TS_burden(cdh, ssc, genebody_HDE, disease="CDH", sites="TSS", left=20000, right=20000, by=5000, p_val_bin_size=0.01, filename_string="CDH_HDE_TSS_region_burdens")
blah <- test_TS_burden(cdh, ssc, genebody_HDE, disease="CDH", sites="TES", left=20000, right=20000, by=5000, p_val_bin_size=0.01, filename_string="CDH_HDE_TES_region_burdens")
blah <- test_TS_burden(cdh, ssc, genebody, disease="CDH", sites="TSS", left=20000, right=20000, by=5000, p_val_bin_size=0.01, filename_string="CDH_TSS_region_burdens")
blah <- test_TS_burden(cdh, ssc, genebody, disease="CDH", sites="TES", left=20000, right=20000, by=5000, p_val_bin_size=0.01, filename_string="CDH_TES_region_burdens")
# CHD TS burden analysis
blah <- test_TS_burden(chd, sscfb, genebody_HHE, disease="CHD", sites="TSS", left=20000, right=20000, by=5000, p_val_bin_size=0.01, filename_string="CHD_HHE_TSS_region_burdens")
blah <- test_TS_burden(chd, sscfb, genebody_HHE, disease="CHD", sites="TES", left=20000, right=20000, by=5000, p_val_bin_size=0.01, filename_string="CHD_HHE_TES_region_burdens")
blah <- test_TS_burden(chd, sscfb, genebody, disease="CHD", sites="TSS", left=20000, right=20000, by=5000, p_val_bin_size=0.01, filename_string="CHD_TSS_region_burdens")
blah <- test_TS_burden(chd, sscfb, genebody, disease="CHD", sites="TES", left=20000, right=20000, by=5000, p_val_bin_size=0.01, filename_string="CHD_TES_region_burdens")

#################################################################################################################
# Builds and saves gene-focused mutation count data frames.
#################################################################################################################
# Reads ExAC data.
exac_dat <- read.csv("fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt", sep="\t")

# Builds mutation count data frame 
build_mutations_dat <- function(cases_snps, cases_indels, controls_snps, controls_indels) {
    CDH_HE_TS_20k_snps_genes <- paste0(cases_snps$TS_gene); CDH_HE_TS_20k_snps_genes <- CDH_HE_TS_20k_snps_genes[CDH_HE_TS_20k_snps_genes != "."]
    CDH_HE_TS_20k_indels_genes <- paste0(cases_indels$TS_gene); CDH_HE_TS_20k_indels_genes <- CDH_HE_TS_20k_indels_genes[CDH_HE_TS_20k_indels_genes != "."]
    SSC_HE_TS_20k_snps_genes <- paste0(controls_snps$TS_gene); SSC_HE_TS_20k_snps_genes <- SSC_HE_TS_20k_snps_genes[SSC_HE_TS_20k_snps_genes != "."]
    SSC_HE_TS_20k_indels_genes <- paste0(controls_indels$TS_gene); SSC_HE_TS_20k_indels_genes <- SSC_HE_TS_20k_indels_genes[SSC_HE_TS_20k_indels_genes != "."]
    
    # Store (CDH_snp_count, CDH_indel_count, SSC_snp_count, SSC_indel_count) pairs for all implicated genes for the specified sites ("TES" or "TSS").
    TS_genes <- new.env()
    for(genes in CDH_HE_TS_20k_snps_genes) {
        for(gene in strsplit(genes, ",")[[1]]) {
            if(!is.null(TS_genes[[paste(gene)]])) {
                TS_genes[[paste(gene)]][1] <- as.numeric(TS_genes[[paste(gene)]][1]) + 1
            } else {
                TS_genes[[paste(gene)]] <- c(1, 0, 0, 0)
            }
        }
    }
    for(genes in CDH_HE_TS_20k_indels_genes) {
        for(gene in strsplit(genes, ",")[[1]]) {
            if(!is.null(TS_genes[[paste(gene)]])) {
                TS_genes[[paste(gene)]][2] <- as.numeric(TS_genes[[paste(gene)]][2]) + 1
            } else {
                TS_genes[[paste(gene)]] <- c(0, 1, 0, 0)
            }
        }
    }
    for(genes in SSC_HE_TS_20k_snps_genes) {
        for(gene in strsplit(genes, ",")[[1]]) {
            if(!is.null(TS_genes[[paste(gene)]])) {
                TS_genes[[paste(gene)]][3] <- as.numeric(TS_genes[[paste(gene)]][3]) + 1
            } else {
                TS_genes[[paste(gene)]] <- c(0, 0, 1, 0)
            }
        }
    }
    for(genes in SSC_HE_TS_20k_indels_genes) {
        for(gene in strsplit(genes, ",")[[1]]) {
            if(!is.null(TS_genes[[paste(gene)]])) {
                TS_genes[[paste(gene)]][4] <- as.numeric(TS_genes[[paste(gene)]][4]) + 1
            } else {
                TS_genes[[paste(gene)]] <- c(0, 0, 0, 1)
            }
        }
    }
    TS_genes_dat <- data.frame()
    for(gene in ls(TS_genes)) {
        counts <- TS_genes[[paste(gene)]]
        gene_row <- c(gene, as.numeric(counts))
        #print(paste(gene_row))
        #readline()
        TS_genes_dat <- rbind(TS_genes_dat, gene_row, stringsAsFactors=FALSE)
    }
    colnames(TS_genes_dat) <- c("gene", "cases_snv_count", "cases_indel_count", "controls_snv_count", "controls_indel_count")
    TS_genes_dat$cases_snv_count <- as.numeric(TS_genes_dat$cases_snv_count)
    TS_genes_dat$cases_indel_count <- as.numeric(TS_genes_dat$cases_indel_count)
    TS_genes_dat$controls_snv_count <- as.numeric(TS_genes_dat$controls_snv_count)
    TS_genes_dat$controls_indel_count <- as.numeric(TS_genes_dat$controls_indel_count)
    #TS_genes_dat <- merge(TS_genes_dat, HE_TS_20k_regions[,c("gene", "exp_mutations")], by="gene")
    #TS_genes_dat$exp_mutations <- TS_genes_dat$exp_mutations * SAMPLE_COUNT
    #TS_genes_dat <- cbind(TS_genes_dat, dpois(TS_genes_dat$cases_snv_count, lambda=TS_genes_dat$exp_mutations))
    #colnames(TS_genes_dat)[(ncol(TS_genes_dat)-1):ncol(TS_genes_dat)] <- c("exp_snvs", "dpois(cases_snvs, exp_snvs)")
    #write.csv(TS_genes_dat[as.numeric(TS_genes_dat$CDH_snv_count) + as.numeric(TS_genes_dat$CDH_indel_count) > 1,], file=paste0(DISEASE, "_HE_TSS_genes.csv"), row.names=FALSE)
    
    # Annotate genes with pLI from ExAC
    TS_genes_dat <- cbind(TS_genes_dat, sapply(1:nrow(TS_genes_dat), function(i) return(round(exac_dat$pLI[exac_dat$gene==TS_genes_dat$gene[i]][1], 3))))
    colnames(TS_genes_dat)[ncol(TS_genes_dat)] <- "pLI"
    
    return(TS_genes_dat)
}
# Build CDH recurrence mutation data tables
cdh_TSS_genes_dat <- build_mutations_dat(cdh_snps_TSS, cdh_indels_TSS, ssc_snps_TSS, ssc_indels_TSS)
cdh_TES_genes_dat <- build_mutations_dat(cdh_snps_TES, cdh_indels_TES, ssc_snps_TES, ssc_indels_TES)
write.csv(cdh_TSS_genes_dat, file="CDH_HDE_TSS_genes_all.csv", row.names=FALSE)
write.csv(cdh_TES_genes_dat, file="CDH_HDE_TES_genes_all.csv", row.names=FALSE)
write.csv(cdh_TSS_genes_dat[cdh_TSS_genes_dat$cases_snv_count + cdh_TSS_genes_dat$cases_indel_count >1 & cdh_TSS_genes_dat$controls_snv_count + cdh_TSS_genes_dat$controls_indel_count <1,], file="CDH_HDE_TSS_genes_recurrent.csv", row.names=FALSE)
write.csv(cdh_TES_genes_dat[cdh_TES_genes_dat$cases_snv_count + cdh_TES_genes_dat$cases_indel_count >1 & cdh_TES_genes_dat$controls_snv_count + cdh_TES_genes_dat$controls_indel_count <1,], file="CDH_HDE_TES_genes_recurrent.csv", row.names=FALSE)
# Build CHD recurrence mutation data tables
chd_TSS_genes_dat <- build_mutations_dat(chd_snps_TSS, chd_indels_TSS, sscfb_snps_TSS, sscfb_indels_TSS)
chd_TES_genes_dat <- build_mutations_dat(chd_snps_TES, chd_indels_TES, sscfb_snps_TES, sscfb_indels_TES)
write.csv(chd_TSS_genes_dat, file="CHD_HHE_TSS_genes_all.csv", row.names=FALSE)
write.csv(chd_TES_genes_dat, file="CHD_HHE_TES_genes_all.csv", row.names=FALSE)
write.csv(chd_TSS_genes_dat[chd_TSS_genes_dat$cases_snv_count + chd_TSS_genes_dat$cases_indel_count >1 & chd_TSS_genes_dat$controls_snv_count + chd_TSS_genes_dat$controls_indel_count <1,], file="CHD_HHE_TSS_genes_recurrent.csv", row.names=FALSE)
write.csv(chd_TES_genes_dat[chd_TES_genes_dat$cases_snv_count + chd_TES_genes_dat$cases_indel_count >1 & chd_TES_genes_dat$controls_snv_count + chd_TES_genes_dat$controls_indel_count <1,], file="CHD_HHE_TES_genes_recurrent.csv", row.names=FALSE)

#################################################################################################################
# Run permutation recurrence analysis for null model.
#################################################################################################################
# Permutation approach
permutation <- function(cases_dat, controls_dat, sites, disease, n_greater_than=1, num_iterations=10) {
    #cases_dat_size = nrow(cases_dat)
    #controls_dat_size = nrow(controls_dat)
    case_samples <- unique(cases_dat$sample)
    control_samples <- unique(controls_dat$sample)
    case_num_samples <- length(case_samples)
    control_num_samples <- length(control_samples)
    
    results <- c()
    for(i in 1:num_iterations) {
        print(paste0("Iteration ", i))
        #indices <- sample(seq(1,cases_dat_size+controls_dat_size), cases_dat_size)
        #cases_dat_indices <- indices[indices <= cases_dat_size]
        #controls_dat_indices <- indices[indices > cases_dat_size] - cases_dat_size
        #dat <- rbind(cases_dat[cases_dat_indices,], controls_dat[controls_dat_indices,])
        indices <- sample(seq(1,case_num_samples+control_num_samples), case_num_samples)
        cases_sample_indices <- indices[indices <= case_num_samples]
        controls_sample_indices <- indices[indices > case_num_samples] - case_num_samples
        dat <- rbind(cases_dat[cases_dat$sample %in% case_samples[cases_sample_indices],], controls_dat[controls_dat$sample %in% control_samples[controls_sample_indices],])
        #snps_genes <- is_TS(dat$Chrom, dat$Position, sites=sites, disease=disease, return_genes=TRUE, print_counts=FALSE, samples=dat$sample)
        snps_genes <- is_TS(dat$Chrom, dat$Position, sites=sites, disease=disease, return_genes=TRUE, print_counts=FALSE)
        snps_genes <- snps_genes[snps_genes != ""]
        sites_genes <- new.env()
        for(genes in snps_genes) {
            for(gene in strsplit(genes, ",")[[1]]) {
                if(!is.null(sites_genes[[paste(gene)]])) {
                    sites_genes[[paste(gene)]] <- as.numeric(sites_genes[[paste(gene)]]) + 1
                } else {
                    sites_genes[[paste(gene)]] <- 1
                }
            }
        }
        sites_genes_dat <- data.frame()
        for(gene in ls(sites_genes)) {
            count <- as.numeric(sites_genes[[paste(gene)]])
            gene_row <- c(gene, count)
            sites_genes_dat <- rbind(sites_genes_dat, gene_row, stringsAsFactors=FALSE)
        }
        if(!(nrow(sites_genes_dat) > 0)) { sites_genes_dat <- data.frame(x=numeric(0), y=integer(0)) }
        colnames(sites_genes_dat) <- c("gene", "snv_count")
        sites_genes_dat$snv_count <- as.numeric(sites_genes_dat$snv_count)
        results <- c(results, sum(sites_genes_dat$snv_count > n_greater_than))
    }
    return(results)
}
# CDH permutations
CDH_permutations_TSS_snv <- permutation(cdh[cdh$snv_indel == "snv",], ssc[ssc$snv_indel == "snv",], sites="TSS", disease="CDH", num_iterations=1000)
CDH_permutations_TES_snv <- permutation(cdh[cdh$snv_indel == "snv",], ssc[ssc$snv_indel == "snv",], sites="TES", disease="CDH", num_iterations=1000)
CDH_permutations_TSS_mut <- permutation(cdh, ssc, sites="TSS", disease="CDH", num_iterations=1000)
CDH_permutations_TES_mut <- permutation(cdh, ssc, sites="TES", disease="CDH", num_iterations=1000)
# CHD permutations
CHD_permutations_TSS_snv <- permutation(chd[chd$snv_indel == "snv",], sscfb[sscfb$snv_indel == "snv",], sites="TSS", disease="CHD", num_iterations=1000)
CHD_permutations_TES_snv <- permutation(chd[chd$snv_indel == "snv",], sscfb[sscfb$snv_indel == "snv",], sites="TES", disease="CHD", num_iterations=1000)
CHD_permutations_TSS_mut <- permutation(chd, sscfb, sites="TSS", disease="CHD", num_iterations=1000)
CHD_permutations_TES_mut <- permutation(chd, sscfb, sites="TES", disease="CHD", num_iterations=1000)
CHD_permutations_TES_snv <- permutation(chd[chd$snv_indel == "snv",], sscfb[sscfb$snv_indel == "snv",], sites="TES", disease="CHD", num_iterations=1000)

#################################################################################################################
# Draw plots
#################################################################################################################
# CDH plots
DISEASE="CDH"
SAMPLE_COUNT=195
plot_recurrance(cases_dat=cdh_TSS_genes_dat, mut_type="snv", main=paste0("CDH TSS Region SNV Recurrence Analysis"), n_greater_than=1, cases_permutations=CDH_permutations_TSS_snv, filename=paste0("CDH_TSS_over_1_nc_snv.pdf"))
plot_recurrance(cases_dat=cdh_TES_genes_dat, mut_type="snv", main=paste0("CDH 3'UTR Region SNV Recurrence Analysis"), n_greater_than=1, cases_permutations=CDH_permutations_TES_snv, filename=paste0("CDH_TES_over_1_nc_snv.pdf"))
plot_recurrance(cases_dat=cdh_TSS_genes_dat, mut_type="mut", main=paste0("CDH TSS Region SNV+Indel Recurrence Analysis"), n_greater_than=1, cases_permutations=CDH_permutations_TSS_mut, filename=paste0("CDH_TSS_over_1_nc_mut.pdf"))
plot_recurrance(cases_dat=cdh_TES_genes_dat, mut_type="mut", main=paste0("CDH 3'UTR Region SNV+Indel Recurrence Analysis"), n_greater_than=1, cases_permutations=CDH_permutations_TES_mut, filename=paste0("CDH_TES_over_1_nc_mut.pdf"))

# CHD plots
DISEASE="CHD"
SAMPLE_COUNT=350
plot_recurrance(cases_dat=chd_TSS_genes_dat, mut_type="snv", main=paste0("CHD TSS Region SNV Recurrence Analysis"), n_greater_than=1, cases_permutations=CHD_permutations_TSS_snv, filename=paste0("CHD_TSS_over_1_nc_snv.pdf"))
plot_recurrance(cases_dat=chd_TES_genes_dat, mut_type="snv", main=paste0("CHD 3'UTR Region SNV Recurrence Analysis"), n_greater_than=1, cases_permutations=CHD_permutations_TES_snv, filename=paste0("CHD_TES_over_1_nc_snv.pdf"))
plot_recurrance(cases_dat=chd_TSS_genes_dat, mut_type="mut", main=paste0("CHD TSS Region SNV+Indel Recurrence Analysis"), n_greater_than=1, cases_permutations=CHD_permutations_TSS_mut, filename=paste0("CHD_TSS_over_1_nc_mut.pdf"))
plot_recurrance(cases_dat=chd_TES_genes_dat, mut_type="mut", main=paste0("CHD 3'UTR Region SNV+Indel Recurrence Analysis"), n_greater_than=1, cases_permutations=CHD_permutations_TES_mut, filename=paste0("CHD_TES_over_1_nc_mut.pdf"))






ssc_snps <- read.csv("SSC_WGSA_518_ColumbiaCall.annotated.snp", sep="\t") #"SSC_noncoding.WGSA.annotated.snp"
ssc_indels <- read.csv("SSC_WGSA_518_ColumbiaCall.annotated.indel", sep="\t") #"SSC_noncoding.WGSA.annotated.indel"
ssc_wgsa <- get_wgsa_data(ssc_snps, ssc_indels, disease="CDH")
ssc_snps_TSS <- ssc_wgsa[["snps_TSS"]]
ssc_snps_TES <- ssc_wgsa[["snps_TES"]]
ssc_indels_TSS <- ssc_wgsa[["indels_TSS"]]
ssc_indels_TES <- ssc_wgsa[["indels_TES"]]
ssc_muts_TSS <- ssc_wgsa[["muts_TSS"]]
ssc_muts_TES <- ssc_wgsa[["muts_TES"]]


ssc_indels <- unfactorize(ssc_ind)
for(i in 1:nrow(ssc_indels)) {
    print(i)
    for(j in 355:ncol(ssc_indels)) {
        ssc_indels[i,j] <- gsub("{[0-9]+}", "", ssc_indels[i,j], perl=TRUE)
    }
}
ssc_snps_TSS <- ssc_snps[is_TS(ssc_snps$X.chr, ssc_snps$pos, sites="TSS"),]
ssc_snps_TES <- ssc_snps[is_TS(ssc_snps$X.chr, ssc_snps$pos, sites="TES"),]
ssc_indels_TSS <- ssc_indels[is_TS(ssc_indels$X.chr, ssc_indels$pos, sites="TSS"),]
ssc_indels_TES <- ssc_indels[is_TS(ssc_indels$X.chr, ssc_indels$pos, sites="TES"),]
ssc_snps_TSS <- cbind(is_TS(ssc_snps_TSS$X.chr, ssc_snps_TSS$pos, sites="TSS", return_genes=TRUE), ssc_snps_TSS)
ssc_snps_TES <- cbind(is_TS(ssc_snps_TES$X.chr, ssc_snps_TES$pos, sites="TES", return_genes=TRUE), ssc_snps_TES)
ssc_indels_TSS <- cbind(is_TS(ssc_indels_TSS$X.chr, ssc_indels_TSS$pos, sites="TSS", return_genes=TRUE), ssc_indels_TSS)
ssc_indels_TES <- cbind(is_TS(ssc_indels_TES$X.chr, ssc_indels_TES$pos, sites="TES", return_genes=TRUE), ssc_indels_TES)
colnames(ssc_snps_TSS)[1] <- "TS_gene"
colnames(ssc_snps_TES)[1] <- "TS_gene"
colnames(ssc_indels_TSS)[1] <- "TS_gene"
colnames(ssc_indels_TES)[1] <- "TS_gene"

cdh_snps <- read.csv("CDH_noncoding.WGSA_annotated.snp", sep="\t")
cdh_ind <- read.csv("CDH_noncoding.WGSA_annotated.indel", sep="\t")
cdh_indels <- unfactorize(cdh_ind)
for(i in 1:nrow(cdh_indels)) {
    print(i)
    for(j in 355:ncol(cdh_indels)) {
        cdh_indels[i,j] <- gsub("{[0-9]+}", "", cdh_indels[i,j], perl=TRUE)
    }
}
cdh_snps_TSS <- cdh_snps[is_TS(cdh_snps$X.chr, cdh_snps$pos, sites="TSS"),]
cdh_snps_TES <- cdh_snps[is_TS(cdh_snps$X.chr, cdh_snps$pos, sites="TES"),]
cdh_indels_TSS <- cdh_indels[is_TS(cdh_indels$X.chr, cdh_indels$pos, sites="TSS"),]
cdh_indels_TES <- cdh_indels[is_TS(cdh_indels$X.chr, cdh_indels$pos, sites="TES"),]
cdh_snps_TSS <- cbind(is_TS(cdh_snps_TSS$X.chr, cdh_snps_TSS$pos, sites="TSS", return_genes=TRUE), cdh_snps_TSS)
cdh_snps_TES <- cbind(is_TS(cdh_snps_TES$X.chr, cdh_snps_TES$pos, sites="TES", return_genes=TRUE), cdh_snps_TES)
cdh_indels_TSS <- cbind(is_TS(cdh_indels_TSS$X.chr, cdh_indels_TSS$pos, sites="TSS", return_genes=TRUE), cdh_indels_TSS)
cdh_indels_TES <- cbind(is_TS(cdh_indels_TES$X.chr, cdh_indels_TES$pos, sites="TES", return_genes=TRUE), cdh_indels_TES)
colnames(cdh_snps_TSS)[1] <- "TS_gene"
colnames(cdh_snps_TES)[1] <- "TS_gene"
colnames(cdh_indels_TSS)[1] <- "TS_gene"
colnames(cdh_indels_TES)[1] <- "TS_gene"

shared_features <- intersect(colnames(cdh_snps), colnames(cdh_indels))
cdh_muts <- rbind(cdh_snps[,shared_features], cdh_indels[,shared_features])
cdh_muts_TSS <- cdh_muts[is_TS(cdh_muts$X.chr, cdh_muts$pos, sites="TSS"),]
cdh_muts_TES <- cdh_muts[is_TS(cdh_muts$X.chr, cdh_muts$pos, sites="TES"),]
ssc_muts <- rbind(ssc_snps[,shared_features], ssc_indels[,shared_features])
ssc_muts_TSS <- ssc_muts[is_TS(ssc_muts$X.chr, ssc_muts$pos, sites="TSS"),]
ssc_muts_TES <- ssc_muts[is_TS(ssc_muts$X.chr, ssc_muts$pos, sites="TES"),]
cdh_muts_TSS <- cbind(is_TS(cdh_muts_TSS$X.chr, cdh_muts_TSS$pos, sites="TSS", return_genes=TRUE), cdh_muts_TSS)
cdh_muts_TES <- cbind(is_TS(cdh_muts_TES$X.chr, cdh_muts_TES$pos, sites="TES", return_genes=TRUE), cdh_muts_TES)
ssc_muts_TSS <- cbind(is_TS(ssc_muts_TSS$X.chr, ssc_muts_TSS$pos, sites="TSS", return_genes=TRUE), ssc_muts_TSS)
ssc_muts_TES <- cbind(is_TS(ssc_muts_TES$X.chr, ssc_muts_TES$pos, sites="TES", return_genes=TRUE), ssc_muts_TES)
colnames(cdh_muts_TSS)[1] <- "TS_gene"
colnames(cdh_muts_TES)[1] <- "TS_gene"
colnames(ssc_muts_TSS)[1] <- "TS_gene"
colnames(ssc_muts_TES)[1] <- "TS_gene"

# SNV Burdens
TSS_burdens_env <- analyze_burden(cdh_snps_TSS, ssc_snps_TSS, CDH_SAMPLE_COUNT, SSC_SAMPLE_COUNT)
TSS_burdens <- TSS_burdens_env[["qualitative"]]
epigenome_names <- get_roadmap_epigenome_names(unlist(lapply(seq(1:length(TSS_burdens$feature)), function(i) { return(strsplit(paste(TSS_burdens$feature[i]), "\\.")[[1]][1]) } )))
TSS_burdens <- cbind(TSS_burdens, epigenome_names)
write.csv(TSS_burdens, file="CDH_TSS_snv_burden_analysis.csv", row.names=FALSE)
write.csv(TSS_burdens[which(TSS_burdens$p.value < 0.1),], file="CDH_TSS_snv_burden_analysis_filtered.csv", row.names=FALSE)
TES_burdens_env <- analyze_burden(cdh_snps_TES, ssc_snps_TES)
TES_burdens <- TES_burdens_env[["qualitative"]]
epigenome_names <- get_roadmap_epigenome_names(unlist(lapply(seq(1:length(TES_burdens$feature)), function(i) { return(strsplit(paste(TES_burdens$feature[i]), "\\.")[[1]][1]) } )))
TES_burdens <- cbind(TES_burdens, epigenome_names)
write.csv(TES_burdens, file="CDH_TES_snv_burden_analysis.csv", row.names=FALSE)
write.csv(TES_burdens[which(TES_burdens$p.value < 0.1),], file="CDH_TES_snv_burden_analysis_filtered.csv", row.names=FALSE)
#num_ind_tests_results <- rbind(num_ind_tests_results, cbind("TSS SNV", num_independent_tests(cdh_snps_TSS[,colnames(cdh_snps_TSS) %in% paste0(TSS_burdens[TSS_burdens$cases_Y > 0,]$feature)]), num_independent_tests(ssc_muts_TSS[,colnames(ssc_muts_TSS) %in% paste0(TSS_burdens[TSS_burdens$cases_Y > 0,]$feature)])))
#num_ind_tests_results <- rbind(num_ind_tests_results, cbind("TES SNV", num_independent_tests(cdh_snps_TES[,colnames(cdh_snps_TES) %in% paste0(TES_burdens[TES_burdens$cases_Y > 0,]$feature)]), num_independent_tests(ssc_muts_TES[,colnames(ssc_muts_TES) %in% paste0(TES_burdens[TES_burdens$cases_Y > 0,]$feature)])))

# Indel Burdens
TSS_burdens_env <- analyze_burden(cdh_indels_TSS, ssc_indels_TSS)
TSS_burdens <- TSS_burdens_env[["qualitative"]]
epigenome_names <- get_roadmap_epigenome_names(unlist(lapply(seq(1:length(TSS_burdens$feature)), function(i) { return(strsplit(paste(TSS_burdens$feature[i]), "\\.")[[1]][1]) } )))
TSS_burdens <- cbind(TSS_burdens, epigenome_names)
write.csv(TSS_burdens, file="CDH_TSS_indel_burden_analysis.csv", row.names=FALSE)
write.csv(TSS_burdens[which(TSS_burdens$p.value < 0.1),], file="CDH_TSS_indel_burden_analysis_filtered.csv", row.names=FALSE)
TES_burdens_env <- analyze_burden(cdh_indels_TES, ssc_indels_TES)
TES_burdens <- TES_burdens_env[["qualitative"]]
epigenome_names <- get_roadmap_epigenome_names(unlist(lapply(seq(1:length(TES_burdens$feature)), function(i) { return(strsplit(paste(TES_burdens$feature[i]), "\\.")[[1]][1]) } )))
TES_burdens <- cbind(TES_burdens, epigenome_names)
write.csv(TES_burdens, file="CDH_TES_indel_burden_analysis.csv", row.names=FALSE)
write.csv(TES_burdens[which(TES_burdens$p.value < 0.1),], file="CDH_TES_indel_burden_analysis_filtered.csv", row.names=FALSE)
#num_ind_tests_results <- rbind(num_ind_tests_results, cbind("TSS Indel", num_independent_tests(cdh_indels_TSS[,colnames(cdh_indels_TSS) %in% paste0(TSS_burdens[TSS_burdens$cases_Y > 0,]$feature)]), num_independent_tests(ssc_indels_TSS[,colnames(ssc_indels_TSS) %in% paste0(TSS_burdens[TSS_burdens$cases_Y > 0,]$feature)])))
#num_ind_tests_results <- rbind(num_ind_tests_results, cbind("TES Indel", num_independent_tests(cdh_indels_TES[,colnames(cdh_indels_TES) %in% paste0(TES_burdens[TES_burdens$cases_Y > 0,]$feature)]), num_independent_tests(ssc_indels_TES[,colnames(ssc_indels_TES) %in% paste0(TES_burdens[TES_burdens$cases_Y > 0,]$feature)])))

# SNV+Indel Burdens
TSS_burdens_env <- analyze_burden(cdh_muts_TSS, ssc_muts_TSS, CDH_SAMPLE_COUNT, SSC_SAMPLE_COUNT)
TSS_burdens <- TSS_burdens_env[["qualitative"]]
epigenome_names <- get_roadmap_epigenome_names(unlist(lapply(seq(1:length(TSS_burdens$feature)), function(i) { return(strsplit(paste(TSS_burdens$feature[i]), "\\.")[[1]][1]) } )))
TSS_burdens <- cbind(TSS_burdens, epigenome_names)
write.csv(TSS_burdens, file="CDH_TSS_mut_burden_analysis.csv", row.names=FALSE)
write.csv(TSS_burdens[which(TSS_burdens$p.value < 0.1),], file="CDH_TSS_mut_burden_analysis_filtered.csv", row.names=FALSE)
TES_burdens_env <- analyze_burden(cdh_muts_TES, ssc_muts_TES, CDH_SAMPLE_COUNT, SSC_SAMPLE_COUNT)
TES_burdens <- TES_burdens_env[["qualitative"]]
epigenome_names <- get_roadmap_epigenome_names(unlist(lapply(seq(1:length(TES_burdens$feature)), function(i) { return(strsplit(paste(TES_burdens$feature[i]), "\\.")[[1]][1]) } )))
TES_burdens <- cbind(TES_burdens, epigenome_names)
write.csv(TES_burdens, file="CDH_TES_mut_burden_analysis.csv", row.names=FALSE)
write.csv(TES_burdens[which(TES_burdens$p.value < 0.1),], file="CDH_TES_mut_burden_analysis_filtered.csv", row.names=FALSE)
#num_ind_tests_results <- rbind(num_ind_tests_results, cbind("TSS SNV+Indel", num_independent_tests(cdh_muts_TSS[,colnames(cdh_muts_TSS) %in% paste0(TSS_burdens[TSS_burdens$cases_Y > 0,]$feature)]), num_independent_tests(ssc_muts_TSS[,colnames(ssc_muts_TSS) %in% paste0(TSS_burdens[TSS_burdens$cases_Y > 0,]$feature)])))
#num_ind_tests_results <- rbind(num_ind_tests_results, cbind("TES SNV+Indel", num_independent_tests(cdh_muts_TES[,colnames(cdh_muts_TES) %in% paste0(TES_burdens[TES_burdens$cases_Y > 0,]$feature)]), num_independent_tests(ssc_muts_TES[,colnames(ssc_muts_TES) %in% paste0(TES_burdens[TES_burdens$cases_Y > 0,]$feature)])))

colnames(num_ind_tests_results) <- c("feature set", "# independent tests CDH", "# independent tests SSC")
write.csv(num_ind_tests_results, file="independent_tests_results.csv", row.names=FALSE)

############################### FIXED!!!!!!!!!!!!!!!############## # E008 is H9 stem cell, E003 is H1 stem cell

volcano_plot <- function(disease, burdens, burden_colors=NULL, num_best_hits=10, strong_effect_cutoff=1.25, bonferroni_p.value_cutoff=-1, requested_labels=c(), requested_labels_description="Requested results", filename=NULL, number_labels=TRUE, pval_relative_importance=300, label_cex = 0.7, label_gap=1.1, label_points_max_iterations=200, x_max=NULL, y_max=NULL) {
    cat(paste0("Making volcano plot for ", disease, "..."))
    x_max_provided = !is.null(x_max)
    marginal_p.value_cutoff = 0.05
    stronger_p.value_cutoff = 0.001
    burdens <- burdens[!grepl("RBP0|candidate|HMHDE|HMHHE", burdens$burden),]
    burdens$enrichment <- as.numeric(paste0(burdens$enrichment))
    burdens$p.value <- as.numeric(paste0(burdens$p.value))
    num_total_tests = nrow(burdens)
    num_marginally_significant = sum(burdens$p.value < marginal_p.value_cutoff)
    burdens <- burdens[burdens$enrichment < Inf & burdens$enrichment > 0 & burdens$p.value < marginal_p.value_cutoff,]
    label <- paste0(gsub(" \\(.+\\)","",burdens$burden, perl=TRUE), " (", burdens$variants, ")")
    burdens <- cbind(burdens, label)
    if (grepl("CDH", disease)) { requested_labels <- gsub("HE ", "DE ", requested_labels); requested_labels <- gsub("CHD", "CDH", requested_labels)
    if(!is.null(burden_colors)) { burden_colors$label <- gsub("HE ", "DE ", burden_colors$label); burden_colors$label <- gsub("CHD", "CDH", burden_colors$label) }
    } else if (grepl("CHD", disease)) { requested_labels <- gsub("DE ", "HE ", requested_labels); requested_labels <- gsub("CDH", "CHD", requested_labels)
    if(!is.null(burden_colors)) { burden_colors$label <- gsub("DE ", "HE ", burden_colors$label); burden_colors$label <- gsub("CDH", "CHD", burden_colors$label) } 
    }
    if (length(requested_labels) > 0) { burdens <- burdens[c(which(!(burdens$label %in% requested_labels)), which((burdens$label %in% requested_labels))),] }
    if(length(burden_colors) < 1) { burden_colors <- NULL }
    
    cols <- rep("grey", nrow(burdens))
    if (is.null(burden_colors)) {
        significant_bonferroni <- (burdens$p.value < bonferroni_p.value_cutoff)
        significant_stronger <- (burdens$p.value < stronger_p.value_cutoff & !significant_bonferroni)
        significant_marginal <- (burdens$p.value < marginal_p.value_cutoff & !significant_stronger & !significant_bonferroni)
        requested <- burdens$label %in% requested_labels
        strong_effect <- (burdens$enrichment > strong_effect_cutoff)
        
        cols[requested] <- "black"
        cols[strong_effect] <- "lightblue1"
        cols[requested & strong_effect] <- "blue"
        cols[significant_marginal] <- "green"
        cols[requested & significant_marginal] <- "green4"
        cols[significant_marginal & strong_effect] <- "gray50" #"orange"
        cols[requested & significant_marginal & strong_effect] <- "orange4"
        cols[significant_stronger & strong_effect] <- "darkorange1"
        cols[significant_bonferroni] <- "red"
    } else {
        burden_colors <- data.frame(cbind(paste0(burden_colors$label), 1:nrow(burden_colors))); colnames(burden_colors) <- c("label","rank"); rownames(burden_colors) <- paste0(burden_colors$label)
        color_rank <- burden_colors[paste0(burdens$label),]; color_rank$rank <- as.numeric(paste0(color_rank$rank))
        col_ramp <- colorRampPalette(brewer.pal(11,"RdBu"), bias=3)(sum(!is.na(color_rank$rank)))
        cols[!is.na(color_rank$rank)] <- col_ramp[color_rank$rank[!is.na(color_rank$rank)]]
    }
    
    if(is.null(filename)) { 
        filename = paste0(tolower(disease),"_burden_volcano_plot.pdf"); filename = gsub(" *\\( *", "_", filename); filename = gsub(" *\\) *", "", filename);
        if(!is.null(burden_colors)) { filename <- gsub("plot", "heat_plot", filename) }
    }
    
    if(number_labels) { pdf(file=output_path(filename), width=14)
    } else { pdf(file=output_path(filename)) }
    
    if(!x_max_provided) { x_max = max(log2(burdens$enrichment))+0.5 } #18 # Max value for x axis
    if(is.null(y_max)) { y_max = ceiling(-log10(min(as.numeric(paste0(burdens$p.value))))) } #min(c(5, ceiling(-log10(min(as.numeric(paste0(burdens$p.value)))))))
    cex.axis=1.3; cex.lab=1.3; cex.main=1.2; cex.mtext=1.1
    par_mar <- par()$mar
    if(number_labels) {
        par(mfrow=c(1,2), mar=c(par_mar[1], par_mar[2], par_mar[3], 0))
        if(!x_max_provided) { x_max = max(log2(burdens$enrichment))+0.5 }
    }
    x_min = 0 #min(log2(burdens$enrichment))
    y_min = -log10(0.05) #0
    plot(log2(burdens$enrichment), -log10(burdens$p.value), col=cols, main=paste0(disease, " Burden Volcano Plot"), xlab="log2(Odds Ratio)", ylab="-log10(p.value)", xlim=c(x_min,x_max), ylim=c(y_min,y_max), xaxt="n", xaxs="i", yaxs="i", cex.lab=cex.lab, cex.axis=cex.axis, cex.main=cex.main)
    axis(1, at=floor(min(log2(burdens$enrichment))):x_max, cex.axis=cex.axis, cex.lab=cex.lab)
    left_edge = par("usr")[1]; bottom_edge = par("usr")[3]
    abline(v=0, col="black", lty=2); 
    abline(v=log2(strong_effect_cutoff), col="red", lty=2); #text(log2(strong_effect_cutoff)+0.05, bottom_edge+0.02, adj=c(0,0), paste0(strong_effect_cutoff,"x"), col="red", cex=0.6)
    abline(h=-log10(bonferroni_p.value_cutoff), col="red", lty=2); text(left_edge+0.05, -log10(bonferroni_p.value_cutoff)+0.03, adj=c(0,0), paste0("p=", formatC(bonferroni_p.value_cutoff,format="e",digits=2)," (Bonf.)"), col="red", cex=0.6) 
    abline(h=-log10(stronger_p.value_cutoff), col="red", lty=2); text(left_edge+0.05, -log10(stronger_p.value_cutoff)+0.03, adj=c(0,0), paste0("p=",stronger_p.value_cutoff), col="red", cex=0.6)
    #abline(h=-log10(marginal_p.value_cutoff), col="red", lty=2); text(left_edge+0.05, -log10(marginal_p.value_cutoff)+0.03, adj=c(0,0), paste0("p=",marginal_p.value_cutoff), col="red", cex=0.6)
    mtext("Marginally significant (p < 0.05) tests") #mtext(paste0(num_total_tests, " total tests, ", num_marginally_significant, " (", round(100*num_marginally_significant/num_total_tests,1),"%) pass marginal significance"), cex=cex.mtext)
    if (is.null(burden_colors)) { best_results_col = "gray20" # "red4"
    } else { best_results_col = "black" }
    if(length(requested_labels) > 0 && is.null(burden_colors)) {
        #legend("topleft", legend=c(paste0("Top ", num_best_hits, " results"), requested_labels_description), col=c(best_results_col,"black"), pch=15, cex=0.8)
    } else { #legend("topleft", legend=c(paste0("Top ", num_best_hits, " results")), col=c(best_results_col), pch=15, cex=0.8) }
    }
    burdens <- burdens[order(pval_relative_importance*(-log10(burdens$p.value)/mean(-log10(burdens$p.value))) + (log2(burdens$enrichment)/mean(log2(burdens$enrichment))), decreasing=TRUE),]
    
    cat("Done.\n")
    if (length(requested_labels) > 0) {
        requested_results <- burdens[burdens$label %in% requested_labels,]
    } else { requested_results <- c() }
    if(num_best_hits > 0) {
        best_results <- burdens[burdens$enrichment > strong_effect_cutoff & burdens$p.value < marginal_p.value_cutoff,][0:num_best_hits,]
        cat("Top results: \n"); cat(paste0(1:num_best_hits, ". ", best_results$label, ":   ", round(best_results$enrichment,1), "x enrichment, p=", round(best_results$p.value,4), "\n"))
        best_results_cols <- rep(best_results_col,num_best_hits); best_results_cols[best_results$label %in% requested_labels] <- "red3"
        requested_results <- requested_results[which(!(requested_results$label %in% best_results$label)),]; requested_labels <- requested_labels[!(requested_labels %in% best_results$label)]
        number_labels_legend <- label_points(log2(c(best_results$enrichment,requested_results$enrichment)), -log10(c(best_results$p.value,requested_results$p.value)), cex=label_cex, c(paste0(best_results$label),paste0(requested_results$label)), col=c(best_results_cols,rep("black",length(requested_results))), number_labels=number_labels, gap=label_gap, max_iterations=label_points_max_iterations)
    } else { number_labels_legend <- label_points(log2(requested_results$enrichment), -log10(requested_results$p.value), cex=label_cex, paste0(requested_results$label), col="black", number_labels=number_labels, gap=label_gap, max_iterations=label_points_max_iterations) } 
    
    if(number_labels) {
        par(mar=c(par_mar[1], par_mar[1]/20, par_mar[3], 0))
        plot.new()
        legend("topleft", legend=number_labels_legend, col=c("white"), pch=26, cex=0.8, bty="n")
    }
    
    dev.off()
    par(mfrow=c(1,1), mar=c(par_mar))
    
    return(burdens)
}
best_result_cdh_chd <- volcano_plot("All Cases (CDH+CHD)", cdh_chd_burden_hg19[-c(1:3),], strong_effect_cutoff=0.9, bonferroni_p.value_cutoff=2.19e-4, num_best_hits=10, pval_relative_importance=1000, number_labels=FALSE)
best_result_chdfb_hg19 <- volcano_plot("All Cases (CDH+CHD)", cdh_chd_burden_hg19[-c(1:3),], strong_effect_cutoff=0.9, bonferroni_p.value_cutoff=2.19e-4, num_best_hits=10, pval_relative_importance=1000, number_labels=FALSE)
# CHDFB vs SSCFB
best_result_chdfb_hg19 <- volcano_plot("All CHDFB", chdfb_burden_hg19[-c(1:3),], strong_effect_cutoff=0.9, bonferroni_p.value_cutoff=2.19e-4, num_best_hits=10, number_labels=FALSE)
# CHDFB vs SSCFB Conserved
best_result_chdfb_hg19_cons <- volcano_plot("Conserved CHDFB", chdfb_burden_hg19_cons[-c(1:3),], strong_effect_cutoff=0.9, bonferroni_p.value_cutoff=2.19e-4, num_best_hits=10, number_labels=FALSE)


best_results_chd2 <- volcano_plot("CHD2", chd2_burden, num_best_hits=10, requested_labels=best_results_chd$label[1:10], requested_labels_description="Top 10 in CHD1/SSC1 analysis")
best_results_chd <- volcano_plot("CHD", chd_burden, num_best_hits=10)
volcano_plot_output <- volcano_plot("CHD2", chd2_burden, burden_colors=best_results_chd, requested_labels_description="Rank in CHD1/SSC1 analysis", num_best_hits=20)
best_result_chd_all <- volcano_plot("All CHD (745) All SSC (1606)", chd_combined_burden, num_best_hits=10, requested_labels=best_results_chd$label[1:10], requested_labels_description="Top 10 in CHD1/SSC1 analysis")
volcano_plot_output <- volcano_plot("All CHD (745) All SSC (1606)", chd_combined_burden, burden_colors=best_results_chd, requested_labels_description="Rank in CHD1/SSC1 analysis", num_best_hits=10)
volcano_plot_output <- volcano_plot("All CDH (489) All SSC (1606)", cdh_combined_burden, burden_colors=best_results_cdh, requested_labels_description="Rank in CDH1/SSC1 analysis", num_best_hits=20)
volcano_plot_output <- volcano_plot("Conserved CDH (all 489 samples)", cdh_combined_burden_cons, burden_colors=best_results_cdh, requested_labels_description="Rank in CDH1/SSC1 analysis", num_best_hits=20)
volcano_plot_output <- volcano_plot("CHD", chd_burden, burden_colors=best_results_cdh, requested_labels_description="Rank in CDH1/SSC1 analysis", num_best_hits=20)
#volcano_plot_output <- volcano_plot("CHD2", chd2_burden_norm, burden_colors=best_results_chd, requested_labels_description="Rank in CHD1/SSC1 analysis", num_best_hits=20)

# The gap parameter is the vertical gap, in fraction of label height, to include between labels for better readability.
label_points <- function(x, y, labels, col="black", gap=1.1, x_offset = 0.1, cex=0.5, number_labels=FALSE, max_iterations=200) { # gap=0.8
    number_labels_legend = ""
    if(number_labels) { number_labels_legend <- paste0(seq(1:length(labels)), ": ", labels); labels <- seq(1:length(labels)) }
    
    precision = 10000
    label_widths <- strwidth(labels, font = 12, units = "user", cex = cex) * precision
    label_height <- strheight(labels, font = 12, units = "user", cex = cex)[1] * precision
    #print(paste0("(", label_widths, ", ", label_height, ")"))
    left_edge = par("usr")[1]*precision; right_edge = par("usr")[2]*precision; bottom_edge = par("usr")[3]*precision; top_edge = par("usr")[4]*precision
    dat <- data.frame(labels, round((x+x_offset)*precision), round(((x+x_offset)*precision)+label_widths), rep("",length(labels)), round(y*precision), round((y*precision)+label_height), rep("",length(labels)), rep(1,length(labels)))
    colnames(dat) <- c("labels", "x_start", "x_end", "x_grange", "y_start", "y_end", "y_grange", "seqname")
    solution_found = FALSE; iteration = 1
    while(!solution_found & iteration <= max_iterations) {
        dat$y_start <- dat$y_start - (label_height * gap/2); dat$y_end <- dat$y_end + (label_height * gap/2)
        dat$x_grange <- paste0("1:",dat$x_start,"-",dat$x_end); dat$y_grange <- paste0("1:",dat$y_start,"-",dat$y_end)
        find_conflicts <- function(dat) {
            # Convert full label coordinates into GRanges objects.
            x_granges <- makeGRangesFromDataFrame(dat, start.field="x_start", end.field="x_end")
            y_granges <- makeGRangesFromDataFrame(dat, start.field="y_start", end.field="y_end")
            # Find overlaps, and remove duplicates in hits by selecting the member of overlap pair with smaller label index.
            x_hits <- data.frame(findOverlaps(x_granges, x_granges)); x_hits <- x_hits[x_hits$queryHits < x_hits$subjectHits,]
            y_hits <- data.frame(findOverlaps(y_granges, y_granges)); y_hits <- y_hits[y_hits$queryHits < y_hits$subjectHits,]
            # Find conflicts by combining x and y hits, and return this conflicts data frame.
            conflicts <- x_hits[paste0(x_hits$queryHits,",",x_hits$subjectHits) %in% paste0(y_hits$queryHits,",",y_hits$subjectHits),]
            return(conflicts)
        }
        conflicts <- find_conflicts(dat)
        if (nrow(conflicts) == 0) { solution_found = TRUE
        } else {
            conflict_clusters <- unique(cluster_overlapping_genes(paste0(apply(conflicts, 1, paste0, collapse=",")), verbose=FALSE))
            for(conflict_cluster in conflict_clusters) {
                label_indices <- as.numeric(unlist(strsplit(conflict_cluster, ",")))
                label_indices <- label_indices[order(dat$y_start[label_indices])] # Sort label indices from bottom to top 
                num_labels <- length(label_indices)
                #middle = mean(dat$y_start[label_indices]) + ((1+gap)*label_height/2)
                middle = mean(y[label_indices]*precision) + ((1+gap)*label_height/2)
                bottom = max(c(middle - floor(num_labels/2)*(1+gap)*label_height, bottom_edge+(gap*label_height)))
                top = bottom + (num_labels-1)*(1+gap)*label_height
                if (top > top_edge-label_height) { top = top_edge-label_height; bottom = top - (num_labels-1)*(1+gap)*label_height } # Makes sure no conflict with top edge.
                # Spread out (along y-axis) labels in this conflict cluster from bottom to top, fixing the conflict. This can introduce new conflicts, so we run multiple iterations.
                for(i in 1:num_labels) {
                    label_index <- label_indices[i]
                    dat$y_start[label_index] = bottom + (i-1)*(1+gap)*label_height
                    dat$y_end[label_index] <- bottom + (i)*(1+gap)*label_height
                    dat$y_start[label_index] <- round(dat$y_start[label_index]); dat$y_end[label_index] <- round(dat$y_end[label_index])
                }
            }
            iteration = iteration + 1
        }
    }
    
    text(dat$x_start/precision, dat$y_start/precision, dat$labels, col=col, cex=cex, adj=c(0,0))
    segments(x0=x, y0=y, x1=dat$x_start/precision, y1=(dat$y_start+(label_height/2))/precision, col=col)
    
    return(number_labels_legend)
}
best_result_cdh_hg19 <- volcano_plot("All CDH", data.frame(cdh_burden_hg19_fet[-c(1:3),]), strong_effect_cutoff=0.9, bonferroni_p.value_cutoff=2.19e-4, num_best_hits=10, number_labels=FALSE, x_max=2.0, label_cex=0.75, pval_relative_importance=1000)



best_results_cdh <- volcano_plot("CDH1", cdh_burden, num_best_hits=10)
best_results_cdh_new <- volcano_plot("CDH2 (SSC2)", cdh_new_burden, num_best_hits=10, requested_labels=best_results_cdh$label[1:10], requested_labels_description="Top 10 in CDH1")
best_results_cdh_new <- volcano_plot("CDH2 (SSC1)", cdh_new_ssc_old_burden, num_best_hits=10, requested_labels=best_results_cdh$label[1:10], requested_labels_description="Top 10 in CDH1")
best_results_cdh <- volcano_plot("CDH1", cdh_burden, num_best_hits=10, requested_labels=best_results_cdh_new$label[1:10], requested_labels_description="Top 10 in CDH2 (SSC2)")
best_results_chd <- volcano_plot("CHD", chd_burden, num_best_hits=10, requested_labels=union(best_results_cdh$label[1:10],best_results_cdh_new$label[1:10]), requested_labels_description="Top 10 in CDH1/CDH2")

volcano_plot_output <- volcano_plot("CDH1", cdh_burden, burden_colors=best_results_chd, num_best_hits=10, requested_labels_description="Rank in CHD")
volcano_plot_output <- volcano_plot("CDH2 (SSC2)", cdh_new_burden, burden_colors=best_results_cdh, requested_labels_description="Rank in CDH1", num_best_hits=10)
volcano_plot_output <- volcano_plot("CDH2 (SSC1)", cdh_new_ssc_old_burden, burden_colors=best_results_cdh, requested_labels_description="Rank in CDH1", num_best_hits=10)
volcano_plot_output <- volcano_plot("CHD", chd_burden, burden_colors=best_results_cdh, requested_labels_description="Rank in CDH1", num_best_hits=10)


best_results_cdh <- volcano_plot("CDH (combined)", cdh_combined_burden, num_best_hits=10, requested_labels=union(best_results_cdh$label[1:10],best_results_cdh_new$label[1:10]), requested_labels_description="Top 10 in CDH1/CDH2")
best_results_ssc <- volcano_plot("SSC1 (SSC2)", ssc_burden, num_best_hits=10, requested_labels=union(best_results_cdh$label[1:10],best_results_cdh_new$label[1:10]), requested_labels_description="Top 10 in CDH1/CDH2")

best_results_ssc_1088 <- volcano_plot("SSC_1088 (SSC1)", ssc_1088_burden, num_best_hits=10)
best_results_cdh_1088 <- volcano_plot("CDH (SSC_1088)", cdh_1088_burden, num_best_hits=10)

best_result_cdh_all <- volcano_plot("All CDH (489) All SSC (1606)", cdh_all_burden, num_best_hits=10, requested_labels=best_results_cdh$label[1:10], requested_labels_description="Top 10 in original analysis")
volcano_plot_output <- volcano_plot("All CDH (489) All SSC (1606)", cdh_all_burden, burden_colors=best_results_cdh, requested_labels_description="Rank in original analysis", num_best_hits=10)

# Statistical test can be "binomial" (default), which is blind to total number of variants, or "fisher_exact", which controls for total number of variants and better handles cases where overall enrichments are not 1.
basic_burden_analysis <- function(disease, recurrent=FALSE, filter_result=FALSE, ssc=NULL, conserved_only=FALSE, process_RBP=TRUE, normalize_by_overall_burden=FALSE, annotate_phenotypes=FALSE, statistical_test="binomial") {
    if(grepl("SSC", disease)) { control_name <- "SSC"; expression_string_token <- "DE"; expression_tissue_name = "E11.5 mouse diaphragm"
    } else if(grepl("CDH", disease)) { control_name <- "SSC"; expression_string_token <- "DE"; expression_tissue_name = "E11.5 mouse diaphragm"
    } else if(grepl("CHDFB", disease)) { control_name <- "SSCFB"; expression_string_token <- "HE"; expression_tissue_name = "E14.5 mouse heart" 
    } else if(grepl("CHD", disease)) { control_name <- "SSC"; expression_string_token <- "HE"; expression_tissue_name = "E14.5 mouse heart" } else { return() }
    if(disease == "CDH2" || disease == "CHD2") { control_name <- "SSC_1088" } else if (disease == "CDH_combined" || disease == "CHD_combined") { control_name <- "SSC_all" } else if (disease == "CHDFB_combined") { control_name <- "SSCFB_combined" } 
    high_expressed_name <- paste0(paste0("H", expression_string_token)); high_expressed_genes <- get(paste0(high_expressed_name, "_genes"))
    high_medium_high_expressed_name <- paste0(paste0("HMH", expression_string_token)); high_medium_high_expressed_genes <- get(paste0(high_medium_high_expressed_name, "_genes"))
    medium_high_expressed_name <- paste0(paste0("MH", expression_string_token)); medium_high_expressed_genes <- get(paste0(medium_high_expressed_name, "_genes"))
    low_expressed_name <- paste0(paste0("L", expression_string_token)); low_expressed_genes <- get(paste0(low_expressed_name, "_genes"))
    medium_low_expressed_name <- paste0(paste0("ML", expression_string_token)); medium_low_expressed_genes <- get(paste0(medium_low_expressed_name, "_genes"))
    if (!is.null(ssc)) { control_name <- ssc }
    if (conserved_only) { TSS_data_suffix = "40" } else { TSS_data_suffix = "" }
    
    # Get case variables for the given disease
    cases_snps <- get(paste0(tolower(disease), "_", "snps_wgsa"))
    cases_indels <- get(paste0(tolower(disease), "_", "indels_wgsa"))
    cases_muts <- get(paste0(tolower(disease), "_", "muts_wgsa"))
    if (statistical_test != "binomial") { n1 = nrow(cases_muts)
    } else { n1 = get(paste0(toupper(disease), "_SAMPLE_COUNT")) } # statistical_test="binomial"
    cases_snps_TSS <- get(paste0(tolower(disease), "_", "snps_TSS", TSS_data_suffix)); cases_snps_TES <- get(paste0(tolower(disease), "_", "snps_TES"))
    cases_indels_TSS <- get(paste0(tolower(disease), "_", "indels_TSS", TSS_data_suffix)); cases_indels_TES <- get(paste0(tolower(disease), "_", "indels_TES"))
    cases_muts_TSS <- get(paste0(tolower(disease), "_", "muts_TSS", TSS_data_suffix)); cases_muts_TES <- get(paste0(tolower(disease), "_", "muts_TES"))
    # Get control variables for the given disease
    controls_snps <- get(paste0(tolower(control_name), "_", "snps_wgsa"))
    controls_indels <- get(paste0(tolower(control_name), "_", "indels_wgsa"))
    controls_muts <- get(paste0(tolower(control_name), "_", "muts_wgsa"))
    if (statistical_test != "binomial") { n0 = nrow(controls_muts)
    } else { n0 = get(paste0(toupper(control_name), "_SAMPLE_COUNT")) } # statistical_test="binomial"
    controls_snps_TSS <- get(paste0(tolower(control_name), "_", "snps_TSS", TSS_data_suffix)); controls_snps_TES <- get(paste0(tolower(control_name), "_", "snps_TES"))
    controls_indels_TSS <- get(paste0(tolower(control_name), "_", "indels_TSS", TSS_data_suffix)); controls_indels_TES <- get(paste0(tolower(control_name), "_", "indels_TES"))
    controls_muts_TSS <- get(paste0(tolower(control_name), "_", "muts_TSS", TSS_data_suffix)); controls_muts_TES <- get(paste0(tolower(control_name), "_", "muts_TES"))
    
    if (conserved_only) {
        cases_snps$phastCons46way_placental[cases_snps$phastCons46way_placental=="."] <- 0; cases_snps$phastCons46way_placental <- as.numeric(paste0(cases_snps$phastCons46way_placental)); cases_snps <- cases_snps[cases_snps$phastCons46way_placental > 0.5,]
        cases_indels$phastCons46way_placental[cases_indels$phastCons46way_placental=="."] <- 0; cases_indels$phastCons46way_placental <- as.numeric(paste0(cases_indels$phastCons46way_placental)); cases_indels <- cases_indels[cases_indels$phastCons46way_placental > 0.5,]
        cases_muts$phastCons46way_placental[cases_muts$phastCons46way_placental=="."] <- 0; cases_muts$phastCons46way_placental <- as.numeric(paste0(cases_muts$phastCons46way_placental)); cases_muts <- cases_muts[cases_muts$phastCons46way_placental > 0.5,]
        cases_snps_TSS$phastCons46way_placental[cases_snps_TSS$phastCons46way_placental=="."] <- 0; cases_snps_TSS$phastCons46way_placental <- as.numeric(paste0(cases_snps_TSS$phastCons46way_placental)); cases_snps_TSS <- cases_snps_TSS[cases_snps_TSS$phastCons46way_placental > 0.5,]
        cases_indels_TSS$phastCons46way_placental[cases_indels_TSS$phastCons46way_placental=="."] <- 0; cases_indels_TSS$phastCons46way_placental <- as.numeric(paste0(cases_indels_TSS$phastCons46way_placental)); cases_indels_TSS <- cases_indels_TSS[cases_indels_TSS$phastCons46way_placental > 0.5,]
        cases_muts_TSS$phastCons46way_placental[cases_muts_TSS$phastCons46way_placental=="."] <- 0; cases_muts_TSS$phastCons46way_placental <- as.numeric(paste0(cases_muts_TSS$phastCons46way_placental)); cases_muts_TSS <- cases_muts_TSS[cases_muts_TSS$phastCons46way_placental > 0.5,]
        cases_snps_TES$phastCons46way_placental[cases_snps_TES$phastCons46way_placental=="."] <- 0; cases_snps_TES$phastCons46way_placental <- as.numeric(paste0(cases_snps_TES$phastCons46way_placental)); cases_snps_TES <- cases_snps_TES[cases_snps_TES$phastCons46way_placental > 0.5,]
        cases_indels_TES$phastCons46way_placental[cases_indels_TES$phastCons46way_placental=="."] <- 0; cases_indels_TES$phastCons46way_placental <- as.numeric(paste0(cases_indels_TES$phastCons46way_placental)); cases_indels_TES <- cases_indels_TES[cases_indels_TES$phastCons46way_placental > 0.5,]
        cases_muts_TES$phastCons46way_placental[cases_muts_TES$phastCons46way_placental=="."] <- 0; cases_muts_TES$phastCons46way_placental <- as.numeric(paste0(cases_muts_TES$phastCons46way_placental)); cases_muts_TES <- cases_muts_TES[cases_muts_TES$phastCons46way_placental > 0.5,]
        controls_snps$phastCons46way_placental[controls_snps$phastCons46way_placental=="."] <- 0; controls_snps$phastCons46way_placental <- as.numeric(paste0(controls_snps$phastCons46way_placental)); controls_snps <- controls_snps[controls_snps$phastCons46way_placental > 0.5,]
        controls_indels$phastCons46way_placental[controls_indels$phastCons46way_placental=="."] <- 0; controls_indels$phastCons46way_placental <- as.numeric(paste0(controls_indels$phastCons46way_placental)); controls_indels <- controls_indels[controls_indels$phastCons46way_placental > 0.5,]
        controls_muts$phastCons46way_placental[controls_muts$phastCons46way_placental=="."] <- 0; controls_muts$phastCons46way_placental <- as.numeric(paste0(controls_muts$phastCons46way_placental)); controls_muts <- controls_muts[controls_muts$phastCons46way_placental > 0.5,]
        controls_snps_TSS$phastCons46way_placental[controls_snps_TSS$phastCons46way_placental=="."] <- 0; controls_snps_TSS$phastCons46way_placental <- as.numeric(paste0(controls_snps_TSS$phastCons46way_placental)); controls_snps_TSS <- controls_snps_TSS[controls_snps_TSS$phastCons46way_placental > 0.5,]
        controls_indels_TSS$phastCons46way_placental[controls_indels_TSS$phastCons46way_placental=="."] <- 0; controls_indels_TSS$phastCons46way_placental <- as.numeric(paste0(controls_indels_TSS$phastCons46way_placental)); controls_indels_TSS <- controls_indels_TSS[controls_indels_TSS$phastCons46way_placental > 0.5,]
        controls_muts_TSS$phastCons46way_placental[controls_muts_TSS$phastCons46way_placental=="."] <- 0; controls_muts_TSS$phastCons46way_placental <- as.numeric(paste0(controls_muts_TSS$phastCons46way_placental)); controls_muts_TSS <- controls_muts_TSS[controls_muts_TSS$phastCons46way_placental > 0.5,]
        controls_snps_TES$phastCons46way_placental[controls_snps_TES$phastCons46way_placental=="."] <- 0; controls_snps_TES$phastCons46way_placental <- as.numeric(paste0(controls_snps_TES$phastCons46way_placental)); controls_snps_TES <- controls_snps_TES[controls_snps_TES$phastCons46way_placental > 0.5,]
        controls_indels_TES$phastCons46way_placental[controls_indels_TES$phastCons46way_placental=="."] <- 0; controls_indels_TES$phastCons46way_placental <- as.numeric(paste0(controls_indels_TES$phastCons46way_placental)); controls_indels_TES <- controls_indels_TES[controls_indels_TES$phastCons46way_placental > 0.5,]
        controls_muts_TES$phastCons46way_placental[controls_muts_TES$phastCons46way_placental=="."] <- 0; controls_muts_TES$phastCons46way_placental <- as.numeric(paste0(controls_muts_TES$phastCons46way_placental)); controls_muts_TES <- controls_muts_TES[controls_muts_TES$phastCons46way_placental > 0.5,]
        #plot(density(as.numeric(paste0(ssc_all_snps_wgsa$phastCons46way_placental_rankscore[ssc_all_snps_wgsa$phastCons46way_placental_rankscore!="."]))), main="phastCons placental rankscore")
    }
    
    # If recurrent argument is TRUE, limit TSS/3'UTR variants to only those annotated with recurrently mutated genes.
    if (recurrent) {
        print("Limiting variants to recurrently mutated genes...")
        genes_TSS_cases <- table(unlist(strsplit(paste0(cases_muts_TSS$TS_gene), ",")))
        recurrent_genes_TSS_cases <- names(genes_TSS_cases[genes_TSS_cases > 1])
        genes_TES_cases <- table(unlist(strsplit(paste0(cases_muts_TES$TS_gene), ",")))
        recurrent_genes_TES_cases <- names(genes_TES_cases[genes_TES_cases > 1])
        genes_TSS_controls <- table(unlist(strsplit(paste0(controls_muts_TSS$TS_gene), ",")))
        recurrent_genes_TSS_controls <- names(genes_TSS_controls[genes_TSS_controls > 1])
        genes_TES_controls <- table(unlist(strsplit(paste0(controls_muts_TES$TS_gene), ",")))
        recurrent_genes_TES_controls <- names(genes_TES_controls[genes_TES_controls > 1])
        cases_snps_TSS <- cases_snps_TSS[includes_at_least_one(cases_snps_TSS$TS_gene, recurrent_genes_TSS_cases),]; cases_snps_TES <- cases_snps_TES[includes_at_least_one(cases_snps_TES$TS_gene, recurrent_genes_TES_cases),]
        cases_indels_TSS <- cases_indels_TSS[includes_at_least_one(cases_indels_TSS$TS_gene, recurrent_genes_TSS_cases),]; cases_indels_TES <- cases_indels_TES[includes_at_least_one(cases_indels_TES$TS_gene, recurrent_genes_TES_cases),]
        cases_muts_TSS <- cases_muts_TSS[includes_at_least_one(cases_muts_TSS$TS_gene, recurrent_genes_TSS_cases),]; cases_muts_TES <- cases_muts_TES[includes_at_least_one(cases_muts_TES$TS_gene, recurrent_genes_TES_cases),]
        controls_snps_TSS <- controls_snps_TSS[includes_at_least_one(controls_snps_TSS$TS_gene, recurrent_genes_TSS_controls),]; controls_snps_TES <- controls_snps_TES[includes_at_least_one(controls_snps_TES$TS_gene, recurrent_genes_TES_controls),]
        controls_indels_TSS <- controls_indels_TSS[includes_at_least_one(controls_indels_TSS$TS_gene, recurrent_genes_TSS_controls),]; controls_indels_TES <- controls_indels_TES[includes_at_least_one(controls_indels_TES$TS_gene, recurrent_genes_TES_controls),]
        controls_muts_TSS <- controls_muts_TSS[includes_at_least_one(controls_muts_TSS$TS_gene, recurrent_genes_TSS_controls),]; controls_muts_TES <- controls_muts_TES[includes_at_least_one(controls_muts_TES$TS_gene, recurrent_genes_TES_controls),]
    }
    # Append all tests to the tests data frame.
    print("Appending basic burden tests...")
    tests <- data.frame("overall", "SNP+indel", nrow(cases_muts), nrow(controls_muts), "", paste0(round(nrow(cases_muts)/get(paste0(toupper(disease),"_SAMPLE_COUNT")),2), " variants/case, ", round(nrow(controls_muts)/get(paste0(toupper(control_name),"_SAMPLE_COUNT")),2), " variants/control, two-sided binomial test"), "", "", "", "", "", stringsAsFactors=FALSE); colnames(tests) <- c("test", "variants", "m1", "m0", "phenotype_burdens", "notes", "notable_genes", "cases_gene_hits", "controls_gene_hits", "cases_variant_hits", "controls_variant_hits")
    tests <- rbind(tests, c("overall", "SNPs", nrow(cases_snps), nrow(controls_snps), "", paste0(round(nrow(cases_snps)/get(paste0(toupper(disease),"_SAMPLE_COUNT")),2), " variants/case, ", round(nrow(controls_snps)/get(paste0(toupper(control_name),"_SAMPLE_COUNT")),2), " variants/control, two-sided binomial test"), "", "", "", "", ""))
    tests <- rbind(tests, c("overall", "indels", nrow(cases_indels), nrow(controls_indels), "", paste0(round(nrow(cases_indels)/get(paste0(toupper(disease),"_SAMPLE_COUNT")),2), " variants/case, ", round(nrow(controls_indels)/get(paste0(toupper(control_name),"_SAMPLE_COUNT")),2), " variants/control, two-sided binomial test"), "", "", "", "", ""))
    
    if(TRUE) {
        # Append a test to the tests data frame.
        append_test <- function(test_name, variant_type, d1, d0, notes="", phenotype_annotation=annotate_phenotypes, gene_hits_annotation=FALSE, variant_hits_annotation=FALSE) {
            cases_have_gene_annotation = ("TS_gene" %in% colnames(d1)); controls_have_gene_annotation = ("TS_gene" %in% colnames(d0))
            notable_genes_annotation = (cases_have_gene_annotation && controls_have_gene_annotation) # If this is true, does annotation for notable genes.
            phenotype_burdens = ""; notable_genes = ""; cases_gene_hits = ""; controls_gene_hits = ""; cases_variant_hits = ""; controls_variant_hits = ""
            # Annotate variants
            if(variant_hits_annotation) {
                if(nrow(d1)>0) { 
                    cases_variant_hits <- paste0(d1$sample, "_chr", d1$X.chr, ":", d1$pos, d1$ref, ">", d1$alt)
                    if(cases_have_gene_annotation) { cases_variant_hits <- paste0(cases_variant_hits, "_(", d1$TS_gene, ")", collapse=";") } else { cases_variant_hits <- paste0(cases_variant_hits, collapse=";") }
                } else { cases_variant_hits <- "-" }
                if(nrow(d0)>0) {
                    controls_variant_hits <- paste0(d0$sample, "_chr", d0$X.chr, ":", d0$pos, d0$ref, ">", d0$alt)
                    if(controls_have_gene_annotation) { controls_variant_hits <- paste0(controls_variant_hits, "_(", d0$TS_gene, ")", collapse=";") } else { controls_variant_hits <- paste0(controls_variant_hits, collapse=";") }
                } else { controls_variant_hits <- "-" }
            }
            # Annotate genes
            if(notable_genes_annotation || (cases_have_gene_annotation && gene_hits_annotation)) { 
                cases_gene_counts <- table(unlist(strsplit(paste0(d1$TS_gene), ",")))
                cases_gene_hits_names <- paste0(names(cases_gene_counts))
                cases_recurrent_genes_dummy <- (cases_gene_counts > 1)
                cases_gene_hits <- cases_gene_hits_names
                cases_gene_hits[cases_recurrent_genes_dummy] <- paste0(cases_gene_hits[cases_recurrent_genes_dummy], "(", cases_gene_counts[cases_recurrent_genes_dummy], ")")
            }
            if(notable_genes_annotation || (controls_have_gene_annotation && gene_hits_annotation)) { 
                controls_gene_counts <- table(unlist(strsplit(paste0(d0$TS_gene), ",")))
                controls_gene_hits_names <- paste0(names(controls_gene_counts))
                controls_recurrent_genes_dummy <- (controls_gene_counts > 1)
                controls_gene_hits <- controls_gene_hits_names
                controls_gene_hits[controls_recurrent_genes_dummy] <- paste0(controls_gene_hits[controls_recurrent_genes_dummy], "(", controls_gene_counts[controls_recurrent_genes_dummy], ")")
            }
            # Notable genes are those that show up in cases but not in controls, AND are either recurrent, constrained, known candidates, or highly expressed.
            if(notable_genes_annotation) {
                if(grepl("CDH", disease)) { known_candidates <- candidate_CDH_genes } else if(grepl("CHD", disease)) { known_candidates <- candidate_CHD_genes } else { known_candidates <- candidate_CDH_genes }
                notable_genes <- paste0(cases_gene_hits[!(cases_gene_hits_names %in% controls_gene_hits_names) & (cases_recurrent_genes_dummy | cases_gene_hits_names %in% unique(c(constrained_genes, known_candidates, high_expressed_genes)))], collapse=",")
            }
            
            if (phenotype_annotation) {
                if (grepl("CDH", disease)) { 
                    all_phenotype_info <- cdh_phenotype_info; phenotype_info <- all_phenotype_info[all_phenotype_info$sample %in% d1$sample,]
                    female_complex = sum(phenotype_info$gender == "F" & phenotype_info$class == "Complex")
                    female_isolated = sum(phenotype_info$gender == "F" & phenotype_info$class == "Isolated")
                    male_complex = sum(phenotype_info$gender == "M" & phenotype_info$class == "Complex")
                    male_isolated = sum(phenotype_info$gender == "M" & phenotype_info$class == "Isolated")
                    all_female_complex = sum(all_phenotype_info$gender == "F" & all_phenotype_info$class == "Complex")
                    all_female_isolated = sum(all_phenotype_info$gender == "F" & all_phenotype_info$class == "Isolated")
                    all_male_complex = sum(all_phenotype_info$gender == "M" & all_phenotype_info$class == "Complex")
                    all_male_isolated = sum(all_phenotype_info$gender == "M" & all_phenotype_info$class == "Isolated")
                    
                    female_binom_result <- binomial_test(female_complex+female_isolated, male_complex+male_isolated, all_female_complex+all_female_isolated, all_male_complex+all_male_isolated, alternative="two.sided")
                    isolated_binom_result <- binomial_test(female_isolated+male_isolated, female_complex+male_complex, all_female_isolated+all_male_isolated, all_female_complex+all_male_complex, alternative="two.sided")
                    female_isolated_binom_result <- binomial_test(female_isolated, female_complex+male_complex+male_isolated, all_female_isolated, all_female_complex+all_male_complex+all_male_isolated, alternative="two.sided")
                    female_complex_binom_result <- binomial_test(female_complex, female_isolated+male_complex+male_isolated, all_female_complex, all_female_isolated+all_male_complex+all_male_isolated, alternative="two.sided")
                    male_isolated_binom_result <- binomial_test(male_isolated, male_complex+female_complex+female_isolated, all_male_isolated, all_male_complex+all_female_complex+all_female_isolated, alternative="two.sided")
                    male_complex_binom_result <- binomial_test(male_complex, male_isolated+female_complex+female_isolated, all_male_complex, all_male_isolated+all_female_complex+all_female_isolated, alternative="two.sided")
                    phenotype_annotation_elements <- c()
                    for(phenotype_test_name in c("female", "isolated", "female_isolated", "female_complex", "male_isolated", "male_complex")) {
                        phenotype_test <- get(paste0(phenotype_test_name,"_binom_result"))
                        phenotype_count = phenotype_test[["m1"]]
                        phenotype_enrichment = phenotype_test[["estimate"]]
                        phenotype_p.value = phenotype_test[["p.value"]]
                        phenotype_annotation_elements <- c(phenotype_annotation_elements, paste0(phenotype_test_name, ": ", phenotype_count, " samples, ", round(phenotype_enrichment,2), "x (p=", formatC(phenotype_p.value,format="e",digits=2), ")"))
                    }
                }
                phenotype_burdens <- paste0(phenotype_annotation_elements, collapse=", ")
            }
            
            if(gene_hits_annotation) { 
                if(length(cases_gene_hits)>0) { cases_gene_hits <- paste0(cases_gene_hits, collapse=",")
                } else { cases_gene_hits <- "-" }
                if(length(controls_gene_hits)>0) { controls_gene_hits <- paste0(controls_gene_hits, collapse=",") 
                } else { controls_gene_hits <- "-" }
            } else { cases_gene_hits = ""; controls_gene_hits = "" }
            
            return(rbind(tests, c(test_name, variant_type, nrow(d1), nrow(d0), phenotype_burdens, notes, notable_genes, cases_gene_hits, controls_gene_hits, cases_variant_hits, controls_variant_hits)))
        }
        
        print("Appending known candidate genes burden tests...")
        if (grepl("CDH", disease)) {
            candidate_genes <- candidate_CDH_genes
            candidate_genes_name = "candidate mouse CDH genes"
            tests <- append_test("candidate mouse CDH genes TSS", "SNP+indel", cases_muts_TSS[includes_at_least_one(cases_muts_TSS$TS_gene, candidate_CDH_genes),], controls_muts_TSS[includes_at_least_one(controls_muts_TSS$TS_gene, candidate_CDH_genes),], paste0(length(candidate_CDH_genes)," candidate CDH genes in mouse, from Lan's file"), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
            tests <- append_test("candidate mouse CDH genes TSS", "SNPs", cases_snps_TSS[includes_at_least_one(cases_snps_TSS$TS_gene, candidate_CDH_genes),], controls_snps_TSS[includes_at_least_one(controls_snps_TSS$TS_gene, candidate_CDH_genes),], paste0(length(candidate_CDH_genes)," candidate CDH genes in mouse, from Lan's file"), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
            tests <- append_test("candidate mouse CDH genes TSS", "indels", cases_indels_TSS[includes_at_least_one(cases_indels_TSS$TS_gene, candidate_CDH_genes),], controls_indels_TSS[includes_at_least_one(controls_indels_TSS$TS_gene, candidate_CDH_genes),], paste0(length(candidate_CDH_genes)," candidate CDH genes in mouse, from Lan's file"), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
            tests <- append_test("candidate mouse CDH genes 3'UTR", "SNP+indel", cases_muts_TES[includes_at_least_one(cases_muts_TES$TS_gene, candidate_CDH_genes),], controls_muts_TES[includes_at_least_one(controls_muts_TES$TS_gene, candidate_CDH_genes),], paste0(length(candidate_CDH_genes)," candidate CDH genes in mouse, from Lan's file"), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
            tests <- append_test("candidate mouse CDH genes 3'UTR", "SNPs", cases_snps_TES[includes_at_least_one(cases_snps_TES$TS_gene, candidate_CDH_genes),], controls_snps_TES[includes_at_least_one(controls_snps_TES$TS_gene, candidate_CDH_genes),], paste0(length(candidate_CDH_genes)," candidate CDH genes in mouse, from Lan's file"), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
            tests <- append_test("candidate mouse CDH genes 3'UTR", "indels", cases_indels_TES[includes_at_least_one(cases_indels_TES$TS_gene, candidate_CDH_genes),], controls_indels_TES[includes_at_least_one(controls_indels_TES$TS_gene, candidate_CDH_genes),], paste0(length(candidate_CDH_genes)," candidate CDH genes in mouse, from Lan's file"), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
        } else if (grepl("CHD", disease)) {
            candidate_genes <- candidate_CHD_genes
            candidate_genes_name = "candidate mouse CHD genes"
            tests <- append_test("candidate mouse CHD genes TSS", "SNPs", cases_snps_TSS[includes_at_least_one(cases_snps_TSS$TS_gene, candidate_CHD_genes),], controls_snps_TSS[includes_at_least_one(controls_snps_TSS$TS_gene, candidate_CHD_genes),], paste0(length(candidate_CHD_genes)," latest candidate CHD genes in mouse"), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
            tests <- append_test("candidate mouse CHD genes TSS", "indels", cases_indels_TSS[includes_at_least_one(cases_indels_TSS$TS_gene, candidate_CHD_genes),], controls_indels_TSS[includes_at_least_one(controls_indels_TSS$TS_gene, candidate_CHD_genes),], paste0(length(candidate_CHD_genes)," latest candidate CHD genes in mouse"), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
            tests <- append_test("candidate mouse CHD genes 3'UTR", "SNPs", cases_snps_TES[includes_at_least_one(cases_snps_TES$TS_gene, candidate_CHD_genes),], controls_snps_TES[includes_at_least_one(controls_snps_TES$TS_gene, candidate_CHD_genes),], paste0(length(candidate_CHD_genes)," latest candidate CHD genes in mouse"), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
            tests <- append_test("candidate mouse CHD genes 3'UTR", "indels", cases_indels_TES[includes_at_least_one(cases_indels_TES$TS_gene, candidate_CHD_genes),], controls_indels_TES[includes_at_least_one(controls_indels_TES$TS_gene, candidate_CHD_genes),], paste0(length(candidate_CHD_genes)," latest candidate CHD genes in mouse"), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
            tests <- append_test("known mouse KO genes TSS", "SNPs", cases_snps_TSS[includes_at_least_one(cases_snps_TSS$TS_gene, known_mouse_ko_genes),], controls_snps_TSS[includes_at_least_one(controls_snps_TSS$TS_gene, known_mouse_ko_genes),], paste0(length(known_mouse_ko_genes)," known KO genes in mouse, from Lan's file"), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
            tests <- append_test("known mouse KO genes TSS", "indels", cases_indels_TSS[includes_at_least_one(cases_indels_TSS$TS_gene, known_mouse_ko_genes),], controls_indels_TSS[includes_at_least_one(controls_indels_TSS$TS_gene, known_mouse_ko_genes),], paste0(length(known_mouse_ko_genes)," known KO genes in mouse, from Lan's file"), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
            tests <- append_test("known mouse KO genes 3'UTR", "SNPs", cases_snps_TES[includes_at_least_one(cases_snps_TES$TS_gene, known_mouse_ko_genes),], controls_snps_TES[includes_at_least_one(controls_snps_TES$TS_gene, known_mouse_ko_genes),], paste0(length(known_mouse_ko_genes)," known KO genes in mouse, from Lan's file"), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
            tests <- append_test("known mouse KO genes 3'UTR", "indels", cases_indels_TES[includes_at_least_one(cases_indels_TES$TS_gene, known_mouse_ko_genes),], controls_indels_TES[includes_at_least_one(controls_indels_TES$TS_gene, known_mouse_ko_genes),], paste0(length(known_mouse_ko_genes)," known KO genes in mouse, from Lan's file"), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
        }
        tests <- append_test("FANTOM5_enhancer_robust", "SNPs", cases_snps[which(cases_snps$FANTOM5_enhancer_robust == "Y"),], controls_snps[which(controls_snps$FANTOM5_enhancer_robust == "Y"),], "Enhancer elements found with CAGE (Cap Analysis of Gene Expression)")
        #tests <- append_test("FANTOM5_enhancer_robust", "indels", cases_indels[which(cases_indels$FANTOM5_enhancer_robust == "Y"),], controls_indels[which(controls_indels$FANTOM5_enhancer_robust == "Y"),], "Enhancer elements found with CAGE (Cap Analysis of Gene Expression)")
        tests <- append_test("CADD > 10", "SNPs", cases_snps[which(as.numeric(paste0(cases_snps$CADD_phred)) > 10),], controls_snps[which(as.numeric(controls_snps$CADD_phred) > 10),], "")
        #tests <- append_test("CADD > 10", "indels", cases_indels[which(as.numeric(paste0(cases_indels$CADD_phred)) > 10),], controls_indels[which(as.numeric(controls_indels$CADD_phred) > 10),], "")
        tests <- append_test("CADD > 15", "SNPs", cases_snps[which(as.numeric(paste0(cases_snps$CADD_phred)) > 15),], controls_snps[which(as.numeric(controls_snps$CADD_phred) > 15),], "")
        #tests <- append_test("CADD > 15", "indels", cases_indels[which(as.numeric(paste0(cases_indels$CADD_phred)) > 15),], controls_indels[which(as.numeric(controls_indels$CADD_phred) > 15),], "")
        tests <- append_test("CADD > 20", "SNPs", cases_snps[which(as.numeric(paste0(cases_snps$CADD_phred)) > 20),], controls_snps[which(as.numeric(controls_snps$CADD_phred) > 20),], "")
        #tests <- append_test("CADD > 20", "indels", cases_indels[which(as.numeric(paste0(cases_indels$CADD_phred)) > 20),], controls_indels[which(as.numeric(controls_indels$CADD_phred) > 20),], "")
        tests <- append_test("Eigen.phred > 10", "SNPs", cases_snps[which(as.numeric(paste0(cases_snps$Eigen.phred)) > 10),], controls_snps[which(as.numeric(paste0(controls_snps$Eigen.phred)) > 10),], "")
        #tests <- append_test("Eigen.phred > 10", "indels", cases_indels[which(as.numeric(paste0(cases_indels$Eigen.phred)) > 10),], controls_indels[which(as.numeric(paste0(controls_indels$Eigen.phred)) > 10),], "")
        tests <- append_test("GERP_RS > 3", "SNPs", cases_snps[which(as.numeric(paste0(cases_snps$GERP_RS)) > 3),], controls_snps[which(as.numeric(paste0(controls_snps$GERP_RS)) > 3),], "")
        #tests <- append_test("GERP_RS > 3", "indels", cases_indels[which(as.numeric(paste0(cases_indels$GERP_RS)) > 3),], controls_indels[which(as.numeric(paste0(controls_indels$GERP_RS)) > 3),], "")
        tests <- append_test("funseq2_noncoding_rankscore > 0.9", "SNPs", cases_snps[which(as.numeric(paste0(cases_snps$funseq2_noncoding_rankscore)) > 0.9),], controls_snps[which(as.numeric(paste0(controls_snps$funseq2_noncoding_rankscore)) > 0.9),], "")
        #tests <- append_test("funseq2_noncoding_rankscore > 0.9", "indels", cases_indels[which(as.numeric(paste0(cases_indels$funseq2_noncoding_rankscore)) > 0.9),], controls_indels[which(as.numeric(paste0(controls_indels$funseq2_noncoding_rankscore)) > 0.9),], "")
        tests <- append_test("TSS", "SNPs", cases_snps_TSS, controls_snps_TSS, gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
        tests <- append_test("TSS", "indels", cases_indels_TSS, controls_indels_TSS, gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
        tests <- append_test("3'UTR", "SNPs", cases_snps_TES, controls_snps_TES, gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
        tests <- append_test("3'UTR", "indels", cases_indels_TES, controls_indels_TES, gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
        print("Appending expression breakdown burden tests...")
        tests <- append_test(paste0(high_expressed_name," TSS"), "SNPs", cases_snps_TSS[includes_at_least_one(cases_snps_TSS$TS_gene, high_expressed_genes),], controls_snps_TSS[includes_at_least_one(controls_snps_TSS$TS_gene, high_expressed_genes),], paste0("High (top 25%) expression in ",expression_tissue_name), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
        tests <- append_test(paste0(high_expressed_name," TSS"), "indels", cases_indels_TSS[includes_at_least_one(cases_indels_TSS$TS_gene, high_expressed_genes),], controls_indels_TSS[includes_at_least_one(controls_indels_TSS$TS_gene, high_expressed_genes),], paste0("High (top 25%) expression in ",expression_tissue_name), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
        tests <- append_test(paste0(high_expressed_name," TSS"), "SNP+indel", cases_muts_TSS[includes_at_least_one(cases_muts_TSS$TS_gene, high_expressed_genes),], controls_muts_TSS[includes_at_least_one(controls_muts_TSS$TS_gene, high_expressed_genes),], paste0("High (top 25%) expression in ",expression_tissue_name), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
        tests <- append_test(paste0(high_expressed_name," 3'UTR"), "SNPs", cases_snps_TES[includes_at_least_one(cases_snps_TES$TS_gene, high_expressed_genes),], controls_snps_TES[includes_at_least_one(controls_snps_TES$TS_gene, high_expressed_genes),], paste0("High (top 25%) expression in ",expression_tissue_name), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
        tests <- append_test(paste0(high_expressed_name," 3'UTR"), "indels", cases_indels_TES[includes_at_least_one(cases_indels_TES$TS_gene, high_expressed_genes),], controls_indels_TES[includes_at_least_one(controls_indels_TES$TS_gene, high_expressed_genes),], paste0("High (top 25%) expression in ",expression_tissue_name), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
        tests <- append_test(paste0(high_expressed_name," 3'UTR"), "SNP+indel", cases_muts_TES[includes_at_least_one(cases_muts_TES$TS_gene, high_expressed_genes),], controls_muts_TES[includes_at_least_one(controls_muts_TES$TS_gene, high_expressed_genes),], paste0("High (top 25%) expression in ",expression_tissue_name), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
        tests <- append_test(paste0(medium_high_expressed_name," TSS"), "SNPs", cases_snps_TSS[includes_at_least_one(cases_snps_TSS$TS_gene, medium_high_expressed_genes),], controls_snps_TSS[includes_at_least_one(controls_snps_TSS$TS_gene, medium_high_expressed_genes),], paste0("Medium-high (50th-75th percentile) expression in ",expression_tissue_name), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
        tests <- append_test(paste0(medium_high_expressed_name," TSS"), "indels", cases_indels_TSS[includes_at_least_one(cases_indels_TSS$TS_gene, medium_high_expressed_genes),], controls_indels_TSS[includes_at_least_one(controls_indels_TSS$TS_gene, medium_high_expressed_genes),], paste0("Medium-high (50th-75th percentile) expression in ",expression_tissue_name), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
        tests <- append_test(paste0(medium_high_expressed_name," TSS"), "SNP+indel", cases_muts_TSS[includes_at_least_one(cases_muts_TSS$TS_gene, medium_high_expressed_genes),], controls_muts_TSS[includes_at_least_one(controls_muts_TSS$TS_gene, medium_high_expressed_genes),], paste0("Medium-high (50th-75th percentile) expression in ",expression_tissue_name), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
        tests <- append_test(paste0(medium_high_expressed_name," 3'UTR"), "SNPs", cases_snps_TES[includes_at_least_one(cases_snps_TES$TS_gene, medium_high_expressed_genes),], controls_snps_TES[includes_at_least_one(controls_snps_TES$TS_gene, medium_high_expressed_genes),], paste0("Medium-high (50th-75th percentile) expression in ",expression_tissue_name), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
        tests <- append_test(paste0(medium_high_expressed_name," 3'UTR"), "indels", cases_indels_TES[includes_at_least_one(cases_indels_TES$TS_gene, medium_high_expressed_genes),], controls_indels_TES[includes_at_least_one(controls_indels_TES$TS_gene, medium_high_expressed_genes),], paste0("Medium-high (50th-75th percentile) expression in ",expression_tissue_name), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
        tests <- append_test(paste0(medium_high_expressed_name," 3'UTR"), "SNP+indel", cases_muts_TES[includes_at_least_one(cases_muts_TES$TS_gene, medium_high_expressed_genes),], controls_muts_TES[includes_at_least_one(controls_muts_TES$TS_gene, medium_high_expressed_genes),], paste0("Medium-high (50th-75th percentile) expression in ",expression_tissue_name), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
        tests <- append_test(paste0(high_medium_high_expressed_name," TSS"), "SNPs", cases_snps_TSS[includes_at_least_one(cases_snps_TSS$TS_gene, high_medium_high_expressed_genes),], controls_snps_TSS[includes_at_least_one(controls_snps_TSS$TS_gene, high_medium_high_expressed_genes),], paste0("High and medium-high (top 50% percentile) expression in ",expression_tissue_name), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
        tests <- append_test(paste0(high_medium_high_expressed_name," TSS"), "indels", cases_indels_TSS[includes_at_least_one(cases_indels_TSS$TS_gene, high_medium_high_expressed_genes),], controls_indels_TSS[includes_at_least_one(controls_indels_TSS$TS_gene, high_medium_high_expressed_genes),], paste0("High and medium-high (top 50% percentile) expression in ",expression_tissue_name), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
        tests <- append_test(paste0(high_medium_high_expressed_name," TSS"), "SNP+indel", cases_muts_TSS[includes_at_least_one(cases_muts_TSS$TS_gene, high_medium_high_expressed_genes),], controls_muts_TSS[includes_at_least_one(controls_muts_TSS$TS_gene, high_medium_high_expressed_genes),], paste0("High and medium-high (top 50% percentile) expression in ",expression_tissue_name), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
        tests <- append_test(paste0(high_medium_high_expressed_name," 3'UTR"), "SNPs", cases_snps_TES[includes_at_least_one(cases_snps_TES$TS_gene, high_medium_high_expressed_genes),], controls_snps_TES[includes_at_least_one(controls_snps_TES$TS_gene, high_medium_high_expressed_genes),], paste0("High and medium-high (top 50% percentile) expression in ",expression_tissue_name), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
        tests <- append_test(paste0(high_medium_high_expressed_name," 3'UTR"), "indels", cases_indels_TES[includes_at_least_one(cases_indels_TES$TS_gene, high_medium_high_expressed_genes),], controls_indels_TES[includes_at_least_one(controls_indels_TES$TS_gene, high_medium_high_expressed_genes),], paste0("High and medium-high (top 50% percentile) expression in ",expression_tissue_name), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
        tests <- append_test(paste0(high_medium_high_expressed_name," 3'UTR"), "SNP+indel", cases_muts_TES[includes_at_least_one(cases_muts_TES$TS_gene, high_medium_high_expressed_genes),], controls_muts_TES[includes_at_least_one(controls_muts_TES$TS_gene, high_medium_high_expressed_genes),], paste0("High and medium-high (top 50% percentile) expression in ",expression_tissue_name), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
        tests <- append_test(paste0(medium_low_expressed_name," TSS"), "SNPs", cases_snps_TSS[includes_at_least_one(cases_snps_TSS$TS_gene, medium_low_expressed_genes),], controls_snps_TSS[includes_at_least_one(controls_snps_TSS$TS_gene, medium_low_expressed_genes),], paste0("Medium-low (25th-50th percentile) expression in ",expression_tissue_name), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
        tests <- append_test(paste0(medium_low_expressed_name," TSS"), "indels", cases_indels_TSS[includes_at_least_one(cases_indels_TSS$TS_gene, medium_low_expressed_genes),], controls_indels_TSS[includes_at_least_one(controls_indels_TSS$TS_gene, medium_low_expressed_genes),], paste0("Medium-low (25th-50th percentile) expression in ",expression_tissue_name), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
        tests <- append_test(paste0(medium_low_expressed_name," TSS"), "SNP+indel", cases_muts_TSS[includes_at_least_one(cases_muts_TSS$TS_gene, medium_low_expressed_genes),], controls_muts_TSS[includes_at_least_one(controls_muts_TSS$TS_gene, medium_low_expressed_genes),], paste0("Medium-low (25th-50th percentile) expression in ",expression_tissue_name), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
        tests <- append_test(paste0(medium_low_expressed_name," 3'UTR"), "SNPs", cases_snps_TES[includes_at_least_one(cases_snps_TES$TS_gene, medium_low_expressed_genes),], controls_snps_TES[includes_at_least_one(controls_snps_TES$TS_gene, medium_low_expressed_genes),], paste0("Medium-low (25th-50th percentile) expression in ",expression_tissue_name), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
        tests <- append_test(paste0(medium_low_expressed_name," 3'UTR"), "indels", cases_indels_TES[includes_at_least_one(cases_indels_TES$TS_gene, medium_low_expressed_genes),], controls_indels_TES[includes_at_least_one(controls_indels_TES$TS_gene, medium_low_expressed_genes),], paste0("Medium-low (25th-50th percentile) expression in ",expression_tissue_name), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
        tests <- append_test(paste0(medium_low_expressed_name," 3'UTR"), "SNP+indel", cases_muts_TES[includes_at_least_one(cases_muts_TES$TS_gene, medium_low_expressed_genes),], controls_muts_TES[includes_at_least_one(controls_muts_TES$TS_gene, medium_low_expressed_genes),], paste0("Medium-low (25th-50th percentile) expression in ",expression_tissue_name), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
        tests <- append_test(paste0(low_expressed_name," TSS"), "SNPs", cases_snps_TSS[includes_at_least_one(cases_snps_TSS$TS_gene, low_expressed_genes),], controls_snps_TSS[includes_at_least_one(controls_snps_TSS$TS_gene, low_expressed_genes),], paste0("Low (bottom 25%) expression in ",expression_tissue_name), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
        tests <- append_test(paste0(low_expressed_name," TSS"), "indels", cases_indels_TSS[includes_at_least_one(cases_indels_TSS$TS_gene, low_expressed_genes),], controls_indels_TSS[includes_at_least_one(controls_indels_TSS$TS_gene, low_expressed_genes),], paste0("Low (bottom 25%) expression in ",expression_tissue_name), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
        tests <- append_test(paste0(low_expressed_name," TSS"), "SNP+indel", cases_muts_TSS[includes_at_least_one(cases_muts_TSS$TS_gene, low_expressed_genes),], controls_muts_TSS[includes_at_least_one(controls_muts_TSS$TS_gene, low_expressed_genes),], paste0("Low (bottom 25%) expression in ",expression_tissue_name), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
        tests <- append_test(paste0(low_expressed_name," 3'UTR"), "SNPs", cases_snps_TES[includes_at_least_one(cases_snps_TES$TS_gene, low_expressed_genes),], controls_snps_TES[includes_at_least_one(controls_snps_TES$TS_gene, low_expressed_genes),], paste0("Low (bottom 25%) expression in ",expression_tissue_name), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
        tests <- append_test(paste0(low_expressed_name," 3'UTR"), "indels", cases_indels_TES[includes_at_least_one(cases_indels_TES$TS_gene, low_expressed_genes),], controls_indels_TES[includes_at_least_one(controls_indels_TES$TS_gene, low_expressed_genes),], paste0("Low (bottom 25%) expression in ",expression_tissue_name), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
        tests <- append_test(paste0(low_expressed_name," 3'UTR"), "SNP+indel", cases_muts_TES[includes_at_least_one(cases_muts_TES$TS_gene, low_expressed_genes),], controls_muts_TES[includes_at_least_one(controls_muts_TES$TS_gene, low_expressed_genes),], paste0("Low (bottom 25%) expression in ",expression_tissue_name), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
        tests <- append_test("constrained (pLI>0.5) genes TSS", "SNPs", cases_snps_TSS[includes_at_least_one(cases_snps_TSS$TS_gene, constrained_genes),], controls_snps_TSS[includes_at_least_one(controls_snps_TSS$TS_gene, constrained_genes),], gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
        tests <- append_test("constrained (pLI>0.5) genes TSS", "indels", cases_indels_TSS[includes_at_least_one(cases_indels_TSS$TS_gene, constrained_genes),], controls_indels_TSS[includes_at_least_one(controls_indels_TSS$TS_gene, constrained_genes),], gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
        tests <- append_test("constrained (pLI>0.5) genes TSS", "SNP+indel", cases_muts_TSS[includes_at_least_one(cases_muts_TSS$TS_gene, constrained_genes),], controls_muts_TSS[includes_at_least_one(controls_muts_TSS$TS_gene, constrained_genes),], gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
        tests <- append_test("constrained (pLI>0.5) genes 3'UTR", "SNPs", cases_snps_TES[includes_at_least_one(cases_snps_TES$TS_gene, constrained_genes),], controls_snps_TES[includes_at_least_one(controls_snps_TES$TS_gene, constrained_genes),], gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
        tests <- append_test("constrained (pLI>0.5) genes 3'UTR", "indels", cases_indels_TES[includes_at_least_one(cases_indels_TES$TS_gene, constrained_genes),], controls_indels_TES[includes_at_least_one(controls_indels_TES$TS_gene, constrained_genes),], gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
        tests <- append_test("constrained (pLI>0.5) genes 3'UTR", "SNP+indel", cases_muts_TES[includes_at_least_one(cases_muts_TES$TS_gene, constrained_genes),], controls_muts_TES[includes_at_least_one(controls_muts_TES$TS_gene, constrained_genes),], gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
        tests <- append_test("bivalent genes TSS", "SNPs", cases_snps_TSS[includes_at_least_one(cases_snps_TSS$TS_gene, bivalent_genes),], controls_snps_TSS[includes_at_least_one(controls_snps_TSS$TS_gene, bivalent_genes),], paste0(length(bivalent_genes)," genes containing both repressor and activator marks"), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
        tests <- append_test("bivalent genes TSS", "indels", cases_indels_TSS[includes_at_least_one(cases_indels_TSS$TS_gene, bivalent_genes),], controls_indels_TSS[includes_at_least_one(controls_indels_TSS$TS_gene, bivalent_genes),], paste0(length(bivalent_genes)," genes containing both repressor and activator marks"), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
        tests <- append_test("bivalent genes 3'UTR", "SNPs", cases_snps_TES[includes_at_least_one(cases_snps_TES$TS_gene, bivalent_genes),], controls_snps_TES[includes_at_least_one(controls_snps_TES$TS_gene, bivalent_genes),], paste0(length(bivalent_genes)," genes containing both repressor and activator marks"), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
        tests <- append_test("bivalent genes 3'UTR", "indels", cases_indels_TES[includes_at_least_one(cases_indels_TES$TS_gene, bivalent_genes),], controls_indels_TES[includes_at_least_one(controls_indels_TES$TS_gene, bivalent_genes),], paste0(length(bivalent_genes)," genes containing both repressor and activator marks"), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
        tests <- append_test(paste0(high_expressed_name," bivalent genes TSS"), "SNPs", cases_snps_TSS[includes_at_least_one(cases_snps_TSS$TS_gene, intersect(high_expressed_genes, bivalent_genes)),], controls_snps_TSS[includes_at_least_one(controls_snps_TSS$TS_gene, intersect(high_expressed_genes, bivalent_genes)),], gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
        tests <- append_test(paste0(high_expressed_name," bivalent genes TSS"), "indels", cases_indels_TSS[includes_at_least_one(cases_indels_TSS$TS_gene, intersect(high_expressed_genes, bivalent_genes)),], controls_indels_TSS[includes_at_least_one(controls_indels_TSS$TS_gene, intersect(high_expressed_genes, bivalent_genes)),], gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
        tests <- append_test(paste0(high_expressed_name," bivalent genes TSS"), "SNP+indel", cases_muts_TSS[includes_at_least_one(cases_muts_TSS$TS_gene, intersect(high_expressed_genes, bivalent_genes)),], controls_muts_TSS[includes_at_least_one(controls_muts_TSS$TS_gene, intersect(high_expressed_genes, bivalent_genes)),], gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
        tests <- append_test(paste0(high_expressed_name," bivalent genes 3'UTR"), "SNPs", cases_snps_TES[includes_at_least_one(cases_snps_TES$TS_gene, intersect(high_expressed_genes, bivalent_genes)),], controls_snps_TES[includes_at_least_one(controls_snps_TES$TS_gene, intersect(high_expressed_genes, bivalent_genes)),], gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
        tests <- append_test(paste0(high_expressed_name," bivalent genes 3'UTR"), "indels", cases_indels_TES[includes_at_least_one(cases_indels_TES$TS_gene, intersect(high_expressed_genes, bivalent_genes)),], controls_indels_TES[includes_at_least_one(controls_indels_TES$TS_gene, intersect(high_expressed_genes, bivalent_genes)),], gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
        tests <- append_test(paste0(high_expressed_name," bivalent genes 3'UTR"), "SNP+indel", cases_muts_TES[includes_at_least_one(cases_muts_TES$TS_gene, intersect(high_expressed_genes, bivalent_genes)),], controls_muts_TES[includes_at_least_one(controls_muts_TES$TS_gene, intersect(high_expressed_genes, bivalent_genes)),], gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
        tests <- append_test(paste0(high_medium_high_expressed_name," bivalent genes TSS"), "SNPs", cases_snps_TSS[includes_at_least_one(cases_snps_TSS$TS_gene, intersect(high_medium_high_expressed_genes, bivalent_genes)),], controls_snps_TSS[includes_at_least_one(controls_snps_TSS$TS_gene, intersect(high_medium_high_expressed_genes, bivalent_genes)),], gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
        tests <- append_test(paste0(high_medium_high_expressed_name," bivalent genes TSS"), "indels", cases_indels_TSS[includes_at_least_one(cases_indels_TSS$TS_gene, intersect(high_medium_high_expressed_genes, bivalent_genes)),], controls_indels_TSS[includes_at_least_one(controls_indels_TSS$TS_gene, intersect(high_medium_high_expressed_genes, bivalent_genes)),], gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
        tests <- append_test(paste0(high_medium_high_expressed_name," bivalent genes TSS"), "SNP+indel", cases_muts_TSS[includes_at_least_one(cases_muts_TSS$TS_gene, intersect(high_medium_high_expressed_genes, bivalent_genes)),], controls_muts_TSS[includes_at_least_one(controls_muts_TSS$TS_gene, intersect(high_medium_high_expressed_genes, bivalent_genes)),], gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
        tests <- append_test(paste0(high_medium_high_expressed_name," bivalent genes 3'UTR"), "SNPs", cases_snps_TES[includes_at_least_one(cases_snps_TES$TS_gene, intersect(high_medium_high_expressed_genes, bivalent_genes)),], controls_snps_TES[includes_at_least_one(controls_snps_TES$TS_gene, intersect(high_medium_high_expressed_genes, bivalent_genes)),], gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
        tests <- append_test(paste0(high_medium_high_expressed_name," bivalent genes 3'UTR"), "indels", cases_indels_TES[includes_at_least_one(cases_indels_TES$TS_gene, intersect(high_medium_high_expressed_genes, bivalent_genes)),], controls_indels_TES[includes_at_least_one(controls_indels_TES$TS_gene, intersect(high_medium_high_expressed_genes, bivalent_genes)),], gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
        tests <- append_test(paste0(high_medium_high_expressed_name," bivalent genes 3'UTR"), "SNP+indel", cases_muts_TES[includes_at_least_one(cases_muts_TES$TS_gene, intersect(high_medium_high_expressed_genes, bivalent_genes)),], controls_muts_TES[includes_at_least_one(controls_muts_TES$TS_gene, intersect(high_medium_high_expressed_genes, bivalent_genes)),], gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
        
        # Append tests for RBP binding site disruption.
        if (process_RBP) {
            print("Appending RBP burden tests...")
            for(RBP in c("RBP", "RBP0")) {
                if(RBP == "RBP") { pillow_info_string = "+/- 50bp " } else if(RBP == "RBP0") { pillow_info_string = "" } else { pillow_info_string = paste0("+/- ",gsub("RBP","",RBP),"bp ") }
                tests <- append_test(RBP, "SNPs", cases_snps[which(!is.na(cases_snps[,RBP])),], controls_snps[which(!is.na(controls_snps[,RBP])),], paste0("Variants disrupting a binding site ",pillow_info_string,"of at least one of the 160 RBPs."), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
                tests <- append_test(RBP, "indels", cases_indels[which(!is.na(cases_indels[,RBP])),], controls_indels[which(!is.na(controls_indels[,RBP])),], paste0("Variants disrupting a binding site ",pillow_info_string,"of at least one of the 160 RBPs."), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
                tests <- append_test(RBP, "SNP+indel", cases_muts[which(!is.na(cases_muts[,RBP])),], controls_muts[which(!is.na(controls_muts[,RBP])),], paste0("Variants disrupting a binding site ",pillow_info_string,"of at least one of the 160 RBPs."), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
                tests <- append_test(paste0("TSS ",RBP), "SNPs", cases_snps_TSS[which(!is.na(cases_snps_TSS[,RBP])),], controls_snps_TSS[which(!is.na(controls_snps_TSS[,RBP])),], paste0("Variants disrupting a binding site ",pillow_info_string,"in a TSS of at least one of the 160 RBPs."), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
                tests <- append_test(paste0("TSS ",RBP), "indels", cases_indels_TSS[which(!is.na(cases_indels_TSS[,RBP])),], controls_indels_TSS[which(!is.na(controls_indels_TSS[,RBP])),], paste0("Variants disrupting a binding site ",pillow_info_string,"in a TSS of at least one of the 160 RBPs."), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
                tests <- append_test(paste0("TSS ",RBP), "SNP+indel", cases_muts_TSS[which(!is.na(cases_muts_TSS[,RBP])),], controls_muts_TSS[which(!is.na(controls_muts_TSS[,RBP])),], paste0("Variants disrupting a binding site ",pillow_info_string,"in a TSS of at least one of the 160 RBPs."), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
                tests <- append_test(paste0("3'UTR ",RBP), "SNPs", cases_snps_TES[which(!is.na(cases_snps_TES[,RBP])),], controls_snps_TES[which(!is.na(controls_snps_TES[,RBP])),], paste0("Variants disrupting a binding site ",pillow_info_string,"in a 3'UTR of at least one of the 160 RBPs."), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
                tests <- append_test(paste0("3'UTR ",RBP), "indels", cases_indels_TES[which(!is.na(cases_indels_TES[,RBP])),], controls_indels_TES[which(!is.na(controls_indels_TES[,RBP])),], paste0("Variants disrupting a binding site ",pillow_info_string,"in a 3'UTR of at least one of the 160 RBPs."), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
                tests <- append_test(paste0("3'UTR ",RBP), "SNP+indel", cases_muts_TES[which(!is.na(cases_muts_TES[,RBP])),], controls_muts_TES[which(!is.na(controls_muts_TES[,RBP])),], paste0("Variants disrupting a binding site ",pillow_info_string,"in a 3'UTR of at least one of the 160 RBPs."), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
                tests <- append_test(paste0(high_expressed_name," TSS ",RBP), "SNPs", cases_snps_TSS[!is.na(cases_snps_TSS[,RBP]) & includes_at_least_one(cases_snps_TSS$TS_gene, high_expressed_genes),], controls_snps_TSS[!is.na(controls_snps_TSS[,RBP]) & includes_at_least_one(controls_snps_TSS$TS_gene, high_expressed_genes),], paste0("Variants disrupting a binding site ",pillow_info_string,"in a top 25% expressed TSS of at least one of the 160 RBPs."), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
                tests <- append_test(paste0(high_expressed_name," TSS ",RBP), "indels", cases_indels_TSS[!is.na(cases_indels_TSS[,RBP]) & includes_at_least_one(cases_indels_TSS$TS_gene, high_expressed_genes),], controls_indels_TSS[!is.na(controls_indels_TSS[,RBP]) & includes_at_least_one(controls_indels_TSS$TS_gene, high_expressed_genes),], paste0("Variants disrupting a binding site ",pillow_info_string,"in a top 25% expressed TSS of at least one of the 160 RBPs."), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
                tests <- append_test(paste0(high_expressed_name," TSS ",RBP), "SNP+indel", cases_muts_TSS[!is.na(cases_muts_TSS[,RBP]) & includes_at_least_one(cases_muts_TSS$TS_gene, high_expressed_genes),], controls_muts_TSS[!is.na(controls_muts_TSS[,RBP]) & includes_at_least_one(controls_muts_TSS$TS_gene, high_expressed_genes),], paste0("Variants disrupting a binding site ",pillow_info_string,"in a top 25% expressed TSS of at least one of the 160 RBPs."), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
                tests <- append_test(paste0(high_expressed_name," 3'UTR ",RBP), "SNPs", cases_snps_TES[!is.na(cases_snps_TES[,RBP]) & includes_at_least_one(cases_snps_TES$TS_gene, high_expressed_genes),], controls_snps_TES[!is.na(controls_snps_TES[,RBP]) & includes_at_least_one(controls_snps_TES$TS_gene, high_expressed_genes),], paste0("Variants disrupting a binding site ",pillow_info_string,"in a top 25% expressed 3'UTR of at least one of the 160 RBPs."), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
                tests <- append_test(paste0(high_expressed_name," 3'UTR ",RBP), "indels", cases_indels_TES[!is.na(cases_indels_TES[,RBP]) & includes_at_least_one(cases_indels_TES$TS_gene, high_expressed_genes),], controls_indels_TES[!is.na(controls_indels_TES[,RBP]) & includes_at_least_one(controls_indels_TES$TS_gene, high_expressed_genes),], paste0("Variants disrupting a binding site ",pillow_info_string,"in a top 25% expressed 3'UTR of at least one of the 160 RBPs."), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
                tests <- append_test(paste0(high_expressed_name," 3'UTR ",RBP), "SNP+indel", cases_muts_TES[!is.na(cases_muts_TES[,RBP]) & includes_at_least_one(cases_muts_TES$TS_gene, high_expressed_genes),], controls_muts_TES[!is.na(controls_muts_TES[,RBP]) & includes_at_least_one(controls_muts_TES$TS_gene, high_expressed_genes),], paste0("Variants disrupting a binding site ",pillow_info_string,"in a top 25% expressed 3'UTR of at least one of the 160 RBPs."), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
                tests <- append_test(paste0(high_medium_high_expressed_name," TSS ",RBP), "SNPs", cases_snps_TSS[!is.na(cases_snps_TSS[,RBP]) & includes_at_least_one(cases_snps_TSS$TS_gene, high_medium_high_expressed_genes),], controls_snps_TSS[!is.na(controls_snps_TSS[,RBP]) & includes_at_least_one(controls_snps_TSS$TS_gene, high_medium_high_expressed_genes),], paste0("Variants disrupting a binding site ",pillow_info_string,"in a top 50% expressed TSS of at least one of the 160 RBPs."), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
                tests <- append_test(paste0(high_medium_high_expressed_name," TSS ",RBP), "indels", cases_indels_TSS[!is.na(cases_indels_TSS[,RBP]) & includes_at_least_one(cases_indels_TSS$TS_gene, high_medium_high_expressed_genes),], controls_indels_TSS[!is.na(controls_indels_TSS[,RBP]) & includes_at_least_one(controls_indels_TSS$TS_gene, high_medium_high_expressed_genes),], paste0("Variants disrupting a binding site ",pillow_info_string,"in a top 50% expressed TSS of at least one of the 160 RBPs."), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
                tests <- append_test(paste0(high_medium_high_expressed_name," TSS ",RBP), "SNP+indel", cases_muts_TSS[!is.na(cases_muts_TSS[,RBP]) & includes_at_least_one(cases_muts_TSS$TS_gene, high_medium_high_expressed_genes),], controls_muts_TSS[!is.na(controls_muts_TSS[,RBP]) & includes_at_least_one(controls_muts_TSS$TS_gene, high_medium_high_expressed_genes),], paste0("Variants disrupting a binding site ",pillow_info_string,"in a top 50% expressed TSS of at least one of the 160 RBPs."), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
                tests <- append_test(paste0(high_medium_high_expressed_name," 3'UTR ",RBP), "SNPs", cases_snps_TES[!is.na(cases_snps_TES[,RBP]) & includes_at_least_one(cases_snps_TES$TS_gene, high_medium_high_expressed_genes),], controls_snps_TES[!is.na(controls_snps_TES[,RBP]) & includes_at_least_one(controls_snps_TES$TS_gene, high_medium_high_expressed_genes),], paste0("Variants disrupting a binding site ",pillow_info_string,"in a top 50% expressed 3'UTR of at least one of the 160 RBPs."), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
                tests <- append_test(paste0(high_medium_high_expressed_name," 3'UTR ",RBP), "indels", cases_indels_TES[!is.na(cases_indels_TES[,RBP]) & includes_at_least_one(cases_indels_TES$TS_gene, high_medium_high_expressed_genes),], controls_indels_TES[!is.na(controls_indels_TES[,RBP]) & includes_at_least_one(controls_indels_TES$TS_gene, high_medium_high_expressed_genes),], paste0("Variants disrupting a binding site ",pillow_info_string,"in a top 50% expressed 3'UTR of at least one of the 160 RBPs."), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
                tests <- append_test(paste0(high_medium_high_expressed_name," 3'UTR ",RBP), "SNP+indel", cases_muts_TES[!is.na(cases_muts_TES[,RBP]) & includes_at_least_one(cases_muts_TES$TS_gene, high_medium_high_expressed_genes),], controls_muts_TES[!is.na(controls_muts_TES[,RBP]) & includes_at_least_one(controls_muts_TES$TS_gene, high_medium_high_expressed_genes),], paste0("Variants disrupting a binding site ",pillow_info_string,"in a top 50% expressed 3'UTR of at least one of the 160 RBPs."), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
                tests <- append_test(paste0(low_expressed_name," TSS ",RBP), "SNPs", cases_snps_TSS[!is.na(cases_snps_TSS[,RBP]) & includes_at_least_one(cases_snps_TSS$TS_gene, low_expressed_genes),], controls_snps_TSS[!is.na(controls_snps_TSS[,RBP]) & includes_at_least_one(controls_snps_TSS$TS_gene, low_expressed_genes),], paste0("Variants disrupting a binding site ",pillow_info_string,"in a low expressed TSS of at least one of the 160 RBPs."), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
                tests <- append_test(paste0(low_expressed_name," TSS ",RBP), "indels", cases_indels_TSS[!is.na(cases_indels_TSS[,RBP]) & includes_at_least_one(cases_indels_TSS$TS_gene, low_expressed_genes),], controls_indels_TSS[!is.na(controls_indels_TSS[,RBP]) & includes_at_least_one(controls_indels_TSS$TS_gene, low_expressed_genes),], paste0("Variants disrupting a binding site ",pillow_info_string,"in a low expressed TSS of at least one of the 160 RBPs."), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
                tests <- append_test(paste0(low_expressed_name," TSS ",RBP), "SNP+indel", cases_muts_TSS[!is.na(cases_muts_TSS[,RBP]) & includes_at_least_one(cases_muts_TSS$TS_gene, low_expressed_genes),], controls_muts_TSS[!is.na(controls_muts_TSS[,RBP]) & includes_at_least_one(controls_muts_TSS$TS_gene, low_expressed_genes),], paste0("Variants disrupting a binding site ",pillow_info_string,"in a low expressed TSS of at least one of the 160 RBPs."), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
                tests <- append_test(paste0(low_expressed_name," 3'UTR ",RBP), "SNPs", cases_snps_TES[!is.na(cases_snps_TES[,RBP]) & includes_at_least_one(cases_snps_TES$TS_gene, low_expressed_genes),], controls_snps_TES[!is.na(controls_snps_TES[,RBP]) & includes_at_least_one(controls_snps_TES$TS_gene, low_expressed_genes),], paste0("Variants disrupting a binding site ",pillow_info_string,"in a low expressed 3'UTR of at least one of the 160 RBPs."), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
                tests <- append_test(paste0(low_expressed_name," 3'UTR ",RBP), "indels", cases_indels_TES[!is.na(cases_indels_TES[,RBP]) & includes_at_least_one(cases_indels_TES$TS_gene, low_expressed_genes),], controls_indels_TES[!is.na(controls_indels_TES[,RBP]) & includes_at_least_one(controls_indels_TES$TS_gene, low_expressed_genes),], paste0("Variants disrupting a binding site ",pillow_info_string,"in a low expressed 3'UTR of at least one of the 160 RBPs."), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
                tests <- append_test(paste0(low_expressed_name," 3'UTR ",RBP), "SNP+indel", cases_muts_TES[!is.na(cases_muts_TES[,RBP]) & includes_at_least_one(cases_muts_TES$TS_gene, low_expressed_genes),], controls_muts_TES[!is.na(controls_muts_TES[,RBP]) & includes_at_least_one(controls_muts_TES$TS_gene, low_expressed_genes),], paste0("Variants disrupting a binding site ",pillow_info_string,"in a low expressed 3'UTR of at least one of the 160 RBPs."), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
                tests <- append_test(paste0("constrained (pLI>0.5) genes TSS ",RBP), "SNPs", cases_snps_TSS[!is.na(cases_snps_TSS[,RBP]) & includes_at_least_one(cases_snps_TSS$TS_gene, constrained_genes),], controls_snps_TSS[!is.na(controls_snps_TSS[,RBP]) & includes_at_least_one(controls_snps_TSS$TS_gene, constrained_genes),], paste0("Variants disrupting a binding site ",pillow_info_string,"in a constrained TSS of at least one of the 160 RBPs."), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
                tests <- append_test(paste0("constrained (pLI>0.5) genes TSS ",RBP), "indels", cases_indels_TSS[!is.na(cases_indels_TSS[,RBP]) & includes_at_least_one(cases_indels_TSS$TS_gene, constrained_genes),], controls_indels_TSS[!is.na(controls_indels_TSS[,RBP]) & includes_at_least_one(controls_indels_TSS$TS_gene, constrained_genes),], paste0("Variants disrupting a binding site ",pillow_info_string,"in a constrained TSS of at least one of the 160 RBPs."), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
                tests <- append_test(paste0("constrained (pLI>0.5) genes TSS ",RBP), "SNP+indel", cases_muts_TSS[!is.na(cases_muts_TSS[,RBP]) & includes_at_least_one(cases_muts_TSS$TS_gene, constrained_genes),], controls_muts_TSS[!is.na(controls_muts_TSS[,RBP]) & includes_at_least_one(controls_muts_TSS$TS_gene, constrained_genes),], paste0("Variants disrupting a binding site ",pillow_info_string,"in a constrained TSS of at least one of the 160 RBPs."), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
                tests <- append_test(paste0("constrained (pLI>0.5) genes 3'UTR ",RBP), "SNPs", cases_snps_TES[!is.na(cases_snps_TES[,RBP]) & includes_at_least_one(cases_snps_TES$TS_gene, constrained_genes),], controls_snps_TES[!is.na(controls_snps_TES[,RBP]) & includes_at_least_one(controls_snps_TES$TS_gene, constrained_genes),], paste0("Variants disrupting a binding site ",pillow_info_string,"in a constrained 3'UTR of at least one of the 160 RBPs."), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
                tests <- append_test(paste0("constrained (pLI>0.5) genes 3'UTR ",RBP), "indels", cases_indels_TES[!is.na(cases_indels_TES[,RBP]) & includes_at_least_one(cases_indels_TES$TS_gene, constrained_genes),], controls_indels_TES[!is.na(controls_indels_TES[,RBP]) & includes_at_least_one(controls_indels_TES$TS_gene, constrained_genes),], paste0("Variants disrupting a binding site ",pillow_info_string,"in a constrained 3'UTR of at least one of the 160 RBPs."), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
                tests <- append_test(paste0("constrained (pLI>0.5) genes 3'UTR ",RBP), "SNP+indel", cases_muts_TES[!is.na(cases_muts_TES[,RBP]) & includes_at_least_one(cases_muts_TES$TS_gene, constrained_genes),], controls_muts_TES[!is.na(controls_muts_TES[,RBP]) & includes_at_least_one(controls_muts_TES$TS_gene, constrained_genes),], paste0("Variants disrupting a binding site ",pillow_info_string,"in a constrained 3'UTR of at least one of the 160 RBPs."), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
                tests <- append_test(paste0(candidate_genes_name," TSS ",RBP), "SNPs", cases_snps_TSS[!is.na(cases_snps_TSS[,RBP]) & includes_at_least_one(cases_snps_TSS$TS_gene, candidate_genes),], controls_snps_TSS[!is.na(controls_snps_TSS[,RBP]) & includes_at_least_one(controls_snps_TSS$TS_gene, candidate_genes),], paste0("Variants disrupting a binding site ",pillow_info_string,"in a candidate gene TSS of at least one of the 160 RBPs."), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
                tests <- append_test(paste0(candidate_genes_name," TSS ",RBP), "indels", cases_indels_TSS[!is.na(cases_indels_TSS[,RBP]) & includes_at_least_one(cases_indels_TSS$TS_gene, candidate_genes),], controls_indels_TSS[!is.na(controls_indels_TSS[,RBP]) & includes_at_least_one(controls_indels_TSS$TS_gene, candidate_genes),], paste0("Variants disrupting a binding site ",pillow_info_string,"in a candidate gene TSS of at least one of the 160 RBPs."), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
                tests <- append_test(paste0(candidate_genes_name," TSS ",RBP), "SNP+indel", cases_muts_TSS[!is.na(cases_muts_TSS[,RBP]) & includes_at_least_one(cases_muts_TSS$TS_gene, candidate_genes),], controls_muts_TSS[!is.na(controls_muts_TSS[,RBP]) & includes_at_least_one(controls_muts_TSS$TS_gene, candidate_genes),], paste0("Variants disrupting a binding site ",pillow_info_string,"in a candidate gene TSS of at least one of the 160 RBPs."), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
                tests <- append_test(paste0(candidate_genes_name," 3'UTR ",RBP), "SNPs", cases_snps_TES[!is.na(cases_snps_TES[,RBP]) & includes_at_least_one(cases_snps_TES$TS_gene, candidate_genes),], controls_snps_TES[!is.na(controls_snps_TES[,RBP]) & includes_at_least_one(controls_snps_TES$TS_gene, candidate_genes),], paste0("Variants disrupting a binding site ",pillow_info_string,"in a candidate gene 3'UTR of at least one of the 160 RBPs."), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
                tests <- append_test(paste0(candidate_genes_name," 3'UTR ",RBP), "indels", cases_indels_TES[!is.na(cases_indels_TES[,RBP]) & includes_at_least_one(cases_indels_TES$TS_gene, candidate_genes),], controls_indels_TES[!is.na(controls_indels_TES[,RBP]) & includes_at_least_one(controls_indels_TES$TS_gene, candidate_genes),], paste0("Variants disrupting a binding site ",pillow_info_string,"in a candidate gene 3'UTR of at least one of the 160 RBPs."), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
                tests <- append_test(paste0(candidate_genes_name," 3'UTR ",RBP), "SNP+indel", cases_muts_TES[!is.na(cases_muts_TES[,RBP]) & includes_at_least_one(cases_muts_TES$TS_gene, candidate_genes),], controls_muts_TES[!is.na(controls_muts_TES[,RBP]) & includes_at_least_one(controls_muts_TES$TS_gene, candidate_genes),], paste0("Variants disrupting a binding site ",pillow_info_string,"in a candidate gene 3'UTR of at least one of the 160 RBPs."), gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
            }
            # Append tests for individual RBP binding site disruption
            print("Appending individual RBP burden tests...")
            rbp_list <- paste0(sapply(list.files(path=data_path("ENCODE_eCLIP_peak")), function(x) { strsplit(x, "\\.")[[1]][2] })) 
        }
        
        # Append tests for histone marks features
        print("Appending histone mark burden tests...")
        histone_tests_start_index <- nrow(tests) + 1
        relevant_eids <- get_relevant_roadmap_eids(disease)
        relevant_histone_column_names <- search_colnames(relevant_eids, cases_muts)
        relevant_histone_column_names <- c(relevant_histone_column_names, colnames(cases_muts)[grepl("tissue_majority",colnames(cases_muts))], "H3K27ac_Dickel_heart_enhancer")
        
        # Combined tissues
        if (FALSE) {
            # Get combined tissue feature names for all histone marks.
            combined_tissue_features_env <- new.env()
            for(relevant_histone_column_name in relevant_histone_column_names) {
                histone_feature_split <- strsplit(relevant_histone_column_name, "\\.")[[1]]
                eid <- histone_feature_split[1]
                tissue <- get_roadmap_epigenome_names(eid)
                histone_mark <- histone_feature_split[2]
                peak_type <- histone_feature_split[3]
                if(is.null(combined_tissue_features_env[[histone_mark]])) { combined_tissue_features_env[[histone_mark]] <- relevant_histone_column_name
                } else { combined_tissue_features_env[[histone_mark]] <- c(combined_tissue_features_env[[histone_mark]], relevant_histone_column_name) }
            }
            
            append_combined_features_to_dat <- function(dat) {
                # Append combined tissue features to all data files and the name to relevant_histone_column_names.
                for(combined_tissue_feature_name in ls(combined_tissue_features_env)) {
                    combined_tissue_feature_dummy <- rowSums(dat[,combined_tissue_features_env[[combined_tissue_feature_name]]] == "Y") > 0
                    combined_tissue_feature <- rep("N", length(combined_tissue_feature_dummy)); combined_tissue_feature[combined_tissue_feature_dummy] <- "Y"
                    dat <- cbind(dat, combined_tissue_feature); colnames(dat)[ncol(dat)] <- paste0("relevant.", combined_tissue_feature_name, ".anyPeak")
                }
                return(dat)
            }
            # Append combined features to cases data frames.
            cases_snps <- append_combined_features_to_dat(cases_snps)
            cases_indels <- append_combined_features_to_dat(cases_indels)
            cases_muts <- append_combined_features_to_dat(cases_muts)
            cases_snps_TSS <- append_combined_features_to_dat(cases_snps_TSS)
            cases_indels_TSS <- append_combined_features_to_dat(cases_indels_TSS)
            cases_muts_TSS <- append_combined_features_to_dat(cases_muts_TSS)
            cases_snps_TES <- append_combined_features_to_dat(cases_snps_TES)
            cases_indels_TES <- append_combined_features_to_dat(cases_indels_TES)
            cases_muts_TES <- append_combined_features_to_dat(cases_muts_TES)
            # Append combined features to controls data frames.
            controls_snps <- append_combined_features_to_dat(controls_snps)
            controls_indels <- append_combined_features_to_dat(controls_indels)
            controls_muts <- append_combined_features_to_dat(controls_muts)
            controls_snps_TSS <- append_combined_features_to_dat(controls_snps_TSS)
            controls_indels_TSS <- append_combined_features_to_dat(controls_indels_TSS)
            controls_muts_TSS <- append_combined_features_to_dat(controls_muts_TSS)
            controls_snps_TES <- append_combined_features_to_dat(controls_snps_TES)
            controls_indels_TES <- append_combined_features_to_dat(controls_indels_TES)
            controls_muts_TES <- append_combined_features_to_dat(controls_muts_TES)
            # Append combined feature names to relevant_histone_column_names.
            relevant_histone_column_names <- c(relevant_histone_column_names, paste0("relevant.", ls(combined_tissue_features_env), ".anyPeak"))
        }
        for(relevant_histone_column_name in relevant_histone_column_names) {
            histone_feature_split <- strsplit(relevant_histone_column_name, "\\.")[[1]]
            eid = histone_feature_split[1]
            if (eid == "H3K27ac_Dickel_heart_enhancer") { tissue = "heart"; histone_mark = "H3K27ac"; peak_type = "Dickel" 
            } else if(histone_feature_split[2] == "tissue_majority") { eid = "tissue_majority"; tissue = "tissue_majority"; histone_mark = histone_feature_split[1]; peak_type = histone_feature_split[3]
            } else { tissue = get_roadmap_epigenome_names(eid); histone_mark = histone_feature_split[2]; peak_type = histone_feature_split[3] }
            if(tissue == "-") { tissue <- eid; eid <- paste0(relevant_eids, collapse=",") } # special case for combined tissue features, where the eid is initially "relevant"
            if(peak_type == "narrowPeak") { next } # Only process gappedPeak; results seem almost entirely identical anyways.
            if (!(grepl("me", histone_mark) | grepl("ac", histone_mark))) { next }
            test_name = paste0(histone_mark, " in ", tissue, " (", eid, ", ", peak_type, ")")
            histone_info = get_histone_mark_function(histone_mark)
            tests <- append_test(test_name, "SNPs", cases_snps[which(cases_snps[,relevant_histone_column_name] == "Y"),], controls_snps[which(controls_snps[,relevant_histone_column_name] == "Y"),], histone_info, gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
            tests <- append_test(test_name, "indels", cases_indels[which(cases_indels[,relevant_histone_column_name] == "Y"),], controls_indels[which(controls_indels[,relevant_histone_column_name] == "Y"),], histone_info, gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
            tests <- append_test(test_name, "SNP+indel", cases_muts[which(cases_muts[,relevant_histone_column_name] == "Y"),], controls_muts[which(controls_muts[,relevant_histone_column_name] == "Y"),], histone_info, gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
            tests <- append_test(paste0("TSS ", test_name), "SNPs", cases_snps_TSS[which(cases_snps_TSS[,relevant_histone_column_name] == "Y"),], controls_snps_TSS[which(controls_snps_TSS[,relevant_histone_column_name] == "Y"),], histone_info, gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
            tests <- append_test(paste0("TSS ", test_name), "indels", cases_indels_TSS[which(cases_indels_TSS[,relevant_histone_column_name] == "Y"),], controls_indels_TSS[which(controls_indels_TSS[,relevant_histone_column_name] == "Y"),], histone_info, gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
            tests <- append_test(paste0("TSS ", test_name), "SNP+indel", cases_muts_TSS[which(cases_muts_TSS[,relevant_histone_column_name] == "Y"),], controls_muts_TSS[which(controls_muts_TSS[,relevant_histone_column_name] == "Y"),], histone_info, gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
            tests <- append_test(paste0("3'UTR ", test_name), "SNPs", cases_snps_TES[which(cases_snps_TES[,relevant_histone_column_name] == "Y"),], controls_snps_TES[which(controls_snps_TES[,relevant_histone_column_name] == "Y"),], histone_info, gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
            tests <- append_test(paste0("3'UTR ", test_name), "indels", cases_indels_TES[which(cases_indels_TES[,relevant_histone_column_name] == "Y"),], controls_indels_TES[which(controls_indels_TES[,relevant_histone_column_name] == "Y"),], histone_info, gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
            tests <- append_test(paste0("3'UTR ", test_name), "SNP+indel", cases_muts_TES[which(cases_muts_TES[,relevant_histone_column_name] == "Y"),], controls_muts_TES[which(controls_muts_TES[,relevant_histone_column_name] == "Y"),], histone_info, gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
            tests <- append_test(paste0(high_expressed_name, " TSS ", test_name), "SNPs", cases_snps_TSS[which(cases_snps_TSS[,relevant_histone_column_name] == "Y"),][includes_at_least_one(cases_snps_TSS$TS_gene[which(cases_snps_TSS[,relevant_histone_column_name] == "Y")], high_expressed_genes),], controls_snps_TSS[which(controls_snps_TSS[,relevant_histone_column_name] == "Y"),][includes_at_least_one(controls_snps_TSS$TS_gene[which(controls_snps_TSS[,relevant_histone_column_name] == "Y")], high_expressed_genes),], histone_info, gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
            tests <- append_test(paste0(high_expressed_name, " TSS ", test_name), "indels", cases_indels_TSS[which(cases_indels_TSS[,relevant_histone_column_name] == "Y"),][includes_at_least_one(cases_indels_TSS$TS_gene[which(cases_indels_TSS[,relevant_histone_column_name] == "Y")], high_expressed_genes),], controls_indels_TSS[which(controls_indels_TSS[,relevant_histone_column_name] == "Y"),][includes_at_least_one(controls_indels_TSS$TS_gene[which(controls_indels_TSS[,relevant_histone_column_name] == "Y")], high_expressed_genes),], histone_info, gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
            tests <- append_test(paste0(high_expressed_name, " TSS ", test_name), "SNP+indel", cases_muts_TSS[which(cases_muts_TSS[,relevant_histone_column_name] == "Y"),][includes_at_least_one(cases_muts_TSS$TS_gene[which(cases_muts_TSS[,relevant_histone_column_name] == "Y")], high_expressed_genes),], controls_muts_TSS[which(controls_muts_TSS[,relevant_histone_column_name] == "Y"),][includes_at_least_one(controls_muts_TSS$TS_gene[which(controls_muts_TSS[,relevant_histone_column_name] == "Y")], high_expressed_genes),], histone_info, gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
            tests <- append_test(paste0(high_expressed_name, " 3'UTR ", test_name), "SNPs", cases_snps_TES[which(cases_snps_TES[,relevant_histone_column_name] == "Y"),][includes_at_least_one(cases_snps_TES$TS_gene[which(cases_snps_TES[,relevant_histone_column_name] == "Y")], high_expressed_genes),], controls_snps_TES[which(controls_snps_TES[,relevant_histone_column_name] == "Y"),][includes_at_least_one(controls_snps_TES$TS_gene[which(controls_snps_TES[,relevant_histone_column_name] == "Y")], high_expressed_genes),], histone_info, gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
            tests <- append_test(paste0(high_expressed_name, " 3'UTR ", test_name), "indels", cases_indels_TES[which(cases_indels_TES[,relevant_histone_column_name] == "Y"),][includes_at_least_one(cases_indels_TES$TS_gene[which(cases_indels_TES[,relevant_histone_column_name] == "Y")], high_expressed_genes),], controls_indels_TES[which(controls_indels_TES[,relevant_histone_column_name] == "Y"),][includes_at_least_one(controls_indels_TES$TS_gene[which(controls_indels_TES[,relevant_histone_column_name] == "Y")], high_expressed_genes),], histone_info, gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
            tests <- append_test(paste0(high_expressed_name, " 3'UTR ", test_name), "SNP+indel", cases_muts_TES[which(cases_muts_TES[,relevant_histone_column_name] == "Y"),][includes_at_least_one(cases_muts_TES$TS_gene[which(cases_muts_TES[,relevant_histone_column_name] == "Y")], high_expressed_genes),], controls_muts_TES[which(controls_muts_TES[,relevant_histone_column_name] == "Y"),][includes_at_least_one(controls_muts_TES$TS_gene[which(controls_muts_TES[,relevant_histone_column_name] == "Y")], high_expressed_genes),], histone_info, gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
            tests <- append_test(paste0(high_medium_high_expressed_name, " TSS ", test_name), "SNPs", cases_snps_TSS[which(cases_snps_TSS[,relevant_histone_column_name] == "Y"),][includes_at_least_one(cases_snps_TSS$TS_gene[which(cases_snps_TSS[,relevant_histone_column_name] == "Y")], high_medium_high_expressed_genes),], controls_snps_TSS[which(controls_snps_TSS[,relevant_histone_column_name] == "Y"),][includes_at_least_one(controls_snps_TSS$TS_gene[which(controls_snps_TSS[,relevant_histone_column_name] == "Y")], high_medium_high_expressed_genes),], histone_info, gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
            tests <- append_test(paste0(high_medium_high_expressed_name, " TSS ", test_name), "indels", cases_indels_TSS[which(cases_indels_TSS[,relevant_histone_column_name] == "Y"),][includes_at_least_one(cases_indels_TSS$TS_gene[which(cases_indels_TSS[,relevant_histone_column_name] == "Y")], high_medium_high_expressed_genes),], controls_indels_TSS[which(controls_indels_TSS[,relevant_histone_column_name] == "Y"),][includes_at_least_one(controls_indels_TSS$TS_gene[which(controls_indels_TSS[,relevant_histone_column_name] == "Y")], high_medium_high_expressed_genes),], histone_info, gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
            tests <- append_test(paste0(high_medium_high_expressed_name, " TSS ", test_name), "SNP+indel", cases_muts_TSS[which(cases_muts_TSS[,relevant_histone_column_name] == "Y"),][includes_at_least_one(cases_muts_TSS$TS_gene[which(cases_muts_TSS[,relevant_histone_column_name] == "Y")], high_medium_high_expressed_genes),], controls_muts_TSS[which(controls_muts_TSS[,relevant_histone_column_name] == "Y"),][includes_at_least_one(controls_muts_TSS$TS_gene[which(controls_muts_TSS[,relevant_histone_column_name] == "Y")], high_medium_high_expressed_genes),], histone_info, gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
            tests <- append_test(paste0(high_medium_high_expressed_name, " 3'UTR ", test_name), "SNPs", cases_snps_TES[which(cases_snps_TES[,relevant_histone_column_name] == "Y"),][includes_at_least_one(cases_snps_TES$TS_gene[which(cases_snps_TES[,relevant_histone_column_name] == "Y")], high_medium_high_expressed_genes),], controls_snps_TES[which(controls_snps_TES[,relevant_histone_column_name] == "Y"),][includes_at_least_one(controls_snps_TES$TS_gene[which(controls_snps_TES[,relevant_histone_column_name] == "Y")], high_medium_high_expressed_genes),], histone_info, gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
            tests <- append_test(paste0(high_medium_high_expressed_name, " 3'UTR ", test_name), "indels", cases_indels_TES[which(cases_indels_TES[,relevant_histone_column_name] == "Y"),][includes_at_least_one(cases_indels_TES$TS_gene[which(cases_indels_TES[,relevant_histone_column_name] == "Y")], high_medium_high_expressed_genes),], controls_indels_TES[which(controls_indels_TES[,relevant_histone_column_name] == "Y"),][includes_at_least_one(controls_indels_TES$TS_gene[which(controls_indels_TES[,relevant_histone_column_name] == "Y")], high_medium_high_expressed_genes),], histone_info, gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
            tests <- append_test(paste0(high_medium_high_expressed_name, " 3'UTR ", test_name), "SNP+indel", cases_muts_TES[which(cases_muts_TES[,relevant_histone_column_name] == "Y"),][includes_at_least_one(cases_muts_TES$TS_gene[which(cases_muts_TES[,relevant_histone_column_name] == "Y")], high_medium_high_expressed_genes),], controls_muts_TES[which(controls_muts_TES[,relevant_histone_column_name] == "Y"),][includes_at_least_one(controls_muts_TES$TS_gene[which(controls_muts_TES[,relevant_histone_column_name] == "Y")], high_medium_high_expressed_genes),], histone_info, gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
            tests <- append_test(paste0(candidate_genes_name," TSS ", test_name), "SNPs", cases_snps_TSS[which(cases_snps_TSS[,relevant_histone_column_name] == "Y"),][includes_at_least_one(cases_snps_TSS$TS_gene[which(cases_snps_TSS[,relevant_histone_column_name] == "Y")], candidate_genes),], controls_snps_TSS[which(controls_snps_TSS[,relevant_histone_column_name] == "Y"),][includes_at_least_one(controls_snps_TSS$TS_gene[which(controls_snps_TSS[,relevant_histone_column_name] == "Y")], candidate_genes),], histone_info, gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
            tests <- append_test(paste0(candidate_genes_name," TSS ", test_name), "indels", cases_indels_TSS[which(cases_indels_TSS[,relevant_histone_column_name] == "Y"),][includes_at_least_one(cases_indels_TSS$TS_gene[which(cases_indels_TSS[,relevant_histone_column_name] == "Y")], candidate_genes),], controls_indels_TSS[which(controls_indels_TSS[,relevant_histone_column_name] == "Y"),][includes_at_least_one(controls_indels_TSS$TS_gene[which(controls_indels_TSS[,relevant_histone_column_name] == "Y")], candidate_genes),], histone_info, gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
            tests <- append_test(paste0(candidate_genes_name," TSS ", test_name), "SNP+indel", cases_muts_TSS[which(cases_muts_TSS[,relevant_histone_column_name] == "Y"),][includes_at_least_one(cases_muts_TSS$TS_gene[which(cases_muts_TSS[,relevant_histone_column_name] == "Y")], candidate_genes),], controls_muts_TSS[which(controls_muts_TSS[,relevant_histone_column_name] == "Y"),][includes_at_least_one(controls_muts_TSS$TS_gene[which(controls_muts_TSS[,relevant_histone_column_name] == "Y")], candidate_genes),], histone_info, gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
            tests <- append_test(paste0(candidate_genes_name," 3'UTR ", test_name), "SNPs", cases_snps_TES[which(cases_snps_TES[,relevant_histone_column_name] == "Y"),][includes_at_least_one(cases_snps_TES$TS_gene[which(cases_snps_TES[,relevant_histone_column_name] == "Y")], candidate_genes),], controls_snps_TES[which(controls_snps_TES[,relevant_histone_column_name] == "Y"),][includes_at_least_one(controls_snps_TES$TS_gene[which(controls_snps_TES[,relevant_histone_column_name] == "Y")], candidate_genes),], histone_info, gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
            tests <- append_test(paste0(candidate_genes_name," 3'UTR ", test_name), "indels", cases_indels_TES[which(cases_indels_TES[,relevant_histone_column_name] == "Y"),][includes_at_least_one(cases_indels_TES$TS_gene[which(cases_indels_TES[,relevant_histone_column_name] == "Y")], candidate_genes),], controls_indels_TES[which(controls_indels_TES[,relevant_histone_column_name] == "Y"),][includes_at_least_one(controls_indels_TES$TS_gene[which(controls_indels_TES[,relevant_histone_column_name] == "Y")], candidate_genes),], histone_info, gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
            tests <- append_test(paste0(candidate_genes_name," 3'UTR ", test_name), "SNP+indel", cases_muts_TES[which(cases_muts_TES[,relevant_histone_column_name] == "Y"),][includes_at_least_one(cases_muts_TES$TS_gene[which(cases_muts_TES[,relevant_histone_column_name] == "Y")], candidate_genes),], controls_muts_TES[which(controls_muts_TES[,relevant_histone_column_name] == "Y"),][includes_at_least_one(controls_muts_TES$TS_gene[which(controls_muts_TES[,relevant_histone_column_name] == "Y")], candidate_genes),], histone_info, gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
            if (process_RBP) {
                for(RBP in c("RBP", "RBP0")) {
                    #if(RBP == "RBP") { pillow_info_string = "+/- 50bp " } else if(RBP == "RBP0") { pillow_info_string = "" } else { pillow_info_string = paste0("+/- ",gsub("RBP","",RBP),"bp ") }
                    tests <- append_test(paste0(RBP," ",test_name), "SNPs", cases_snps[which(cases_snps[,relevant_histone_column_name] == "Y" & !is.na(cases_snps[,RBP])),], controls_snps[which(controls_snps[,relevant_histone_column_name] == "Y" & !is.na(controls_snps[,RBP])),], histone_info, gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
                    tests <- append_test(paste0(RBP," ",test_name), "indels", cases_indels[which(cases_indels[,relevant_histone_column_name] == "Y" & !is.na(cases_indels[,RBP])),], controls_indels[which(controls_indels[,relevant_histone_column_name] == "Y" & !is.na(controls_indels[,RBP])),], histone_info, gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
                    tests <- append_test(paste0(RBP," ",test_name), "SNP+indel", cases_muts[which(cases_muts[,relevant_histone_column_name] == "Y" & !is.na(cases_muts[,RBP])),], controls_muts[which(controls_muts[,relevant_histone_column_name] == "Y" & !is.na(controls_muts[,RBP])),], histone_info, gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
                    tests <- append_test(paste0("TSS ",RBP," ",test_name), "SNPs", cases_snps_TSS[which(cases_snps_TSS[,relevant_histone_column_name] == "Y" & !is.na(cases_snps_TSS[,RBP])),], controls_snps_TSS[which(controls_snps_TSS[,relevant_histone_column_name] == "Y" & !is.na(controls_snps_TSS[,RBP])),], histone_info, gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
                    tests <- append_test(paste0("TSS ",RBP," ",test_name), "indels", cases_indels_TSS[which(cases_indels_TSS[,relevant_histone_column_name] == "Y" & !is.na(cases_indels_TSS[,RBP])),], controls_indels_TSS[which(controls_indels_TSS[,relevant_histone_column_name] == "Y" & !is.na(controls_indels_TSS[,RBP])),], histone_info, gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
                    tests <- append_test(paste0("TSS ",RBP," ",test_name), "SNP+indel", cases_muts_TSS[which(cases_muts_TSS[,relevant_histone_column_name] == "Y" & !is.na(cases_muts_TSS[,RBP])),], controls_muts_TSS[which(controls_muts_TSS[,relevant_histone_column_name] == "Y" & !is.na(controls_muts_TSS[,RBP])),], histone_info, gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
                    tests <- append_test(paste0("3'UTR ",RBP," ",test_name), "SNPs", cases_snps_TES[which(cases_snps_TES[,relevant_histone_column_name] == "Y" & !is.na(cases_snps_TES[,RBP])),], controls_snps_TES[which(controls_snps_TES[,relevant_histone_column_name] == "Y" & !is.na(controls_snps_TES[,RBP])),], histone_info, gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
                    tests <- append_test(paste0("3'UTR ",RBP," ",test_name), "indels", cases_indels_TES[which(cases_indels_TES[,relevant_histone_column_name] == "Y" & !is.na(cases_indels_TES[,RBP])),], controls_indels_TES[which(controls_indels_TES[,relevant_histone_column_name] == "Y" & !is.na(controls_indels_TES[,RBP])),], histone_info, gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)
                    tests <- append_test(paste0("3'UTR ",RBP," ",test_name), "SNP+indel", cases_muts_TES[which(cases_muts_TES[,relevant_histone_column_name] == "Y" & !is.na(cases_muts_TES[,RBP])),], controls_muts_TES[which(controls_muts_TES[,relevant_histone_column_name] == "Y" & !is.na(controls_muts_TES[,RBP])),], histone_info, gene_hits_annotation=TRUE, variant_hits_annotation=TRUE)}
            }
        }
    } # REMOVE THIS IF TO STOP SKIP
    histone_tests_end_index <- nrow(tests)
    
    tests$m1 <- as.numeric(tests$m1); tests$m0 <- as.numeric(tests$m0)
    binom_enrichments <- c()
    binom_p.values <- c()
    overall_norm_factors <- c(1, 1, 1); names(overall_norm_factors) <- c("SNP+indel", "SNPs", "indels")
    if (normalize_by_overall_burden) { 
        overall_norm_factors[1] <- (tests$m1[1]/n1)/(tests$m0[1]/n0) 
        overall_norm_factors[2] <- (tests$m1[2]/n1)/(tests$m0[2]/n0) 
        overall_norm_factors[3] <- (tests$m1[3]/n1)/(tests$m0[3]/n0)
    }
    for(i in 1:nrow(tests)) {
        if (statistical_test != "poisson regression" && tests[i,1] == "overall") { 
            binom_result <- binomial_test(tests$m1[i], tests$m0[i], get(paste0(toupper(disease),"_SAMPLE_COUNT")), get(paste0(toupper(control_name),"_SAMPLE_COUNT")), alternative="two.sided")
        } else if (statistical_test == "binomial") { norm_factor = overall_norm_factors[tests[i,2]]; binom_result <- binomial_test(floor(tests$m1[i]/norm_factor), tests$m0[i], n1, n0, alternative="greater") 
        } else if (statistical_test == "poisson regression") { 
            if (tests[i,1] == "overall") {
                tests[i,6] <- gsub(", two-sided binomial test", "", tests[i,6])
                variants = tolower(tests[i,2]); if(variants == "snp+indel") { variants = "muts" }
                case_samples <- table(get(paste0("cases_",variants))$sample)
                control_samples <- table(get(paste0("controls_",variants))$sample)
            } else {
                case_samples <- table(paste0(sapply(unlist(strsplit(tests$cases_variant_hits[i],";")), function(x) strsplit(x,"_")[[1]][1])))
                control_samples <- table(paste0(sapply(unlist(strsplit(tests$controls_variant_hits[i],";")), function(x) strsplit(x,"_")[[1]][1])))
            }
            case_samples <- cbind(names(case_samples), case_samples); colnames(case_samples) <- c("sample", "variant_count")
            control_samples <- cbind(names(control_samples), control_samples);  colnames(control_samples) <- c("sample", "variant_count")
            if(grepl("CHD",disease)) { case_paternal_ages <- get("chd_parental_age") } else { case_paternal_ages <- get("cdh_parental_age") }
            case_paternal_ages <- case_paternal_ages[case_paternal_ages$Blinded.ID %in% cases_muts$sample,]
            control_paternal_ages <- get("ssc_parental_age")[,1:2]; control_paternal_ages <- control_paternal_ages[control_paternal_ages$Blinded.ID %in% controls_muts$sample,]
            case_dat <- merge(case_paternal_ages, case_samples, by.x="Blinded.ID", by.y="sample", all.x=TRUE, all.y=TRUE)
            case_dat$variant_count <- as.numeric(paste0(case_dat$variant_count)); case_dat$Paternal.Age.at.Proband.Birth <- as.numeric(paste0(case_dat$Paternal.Age.at.Proband.Birth))
            case_dat$variant_count[which(is.na(case_dat$variant_count))] <- 0; case_dat$Paternal.Age.at.Proband.Birth[which(is.na(case_dat$Paternal.Age.at.Proband.Birth))] <- mean(case_dat$Paternal.Age.at.Proband.Birth[!is.na(case_dat$Paternal.Age.at.Proband.Birth)])
            control_dat <- merge(control_paternal_ages, control_samples, by.x="Blinded.ID", by.y="sample", all.x=TRUE, all.y=TRUE)
            control_dat$variant_count <- as.numeric(paste0(control_dat$variant_count)); control_dat$Paternal.Age.at.Proband.Birth <- as.numeric(paste0(control_dat$Paternal.Age.at.Proband.Birth))
            control_dat$variant_count[which(is.na(control_dat$variant_count))] <- 0; control_dat$Paternal.Age.at.Proband.Birth[which(is.na(control_dat$Paternal.Age.at.Proband.Birth))] <- mean(control_dat$Paternal.Age.at.Proband.Birth[!is.na(control_dat$Paternal.Age.at.Proband.Birth)])
            variant_count_dat <- rbind(case_dat, control_dat)
            case_control <- rep(0, nrow(variant_count_dat)); case_control[1:nrow(case_dat)] <- 1; variant_count_dat <- cbind(variant_count_dat, case_control)
            
            m1 <- glm(variant_count ~ case_control + Paternal.Age.at.Proband.Birth, family="poisson", data=variant_count_dat) # control = list(maxit = 50)
            #binom_enrichments <- c(binom_enrichments, coef(summary(m1))[2,1])
            #binom_p.values <- c(binom_p.values, coef(summary(m1))[2,4])
            binom_result <- new.env()
            binom_result[["estimate"]] <- coef(summary(m1))[2,1]
            binom_result[["p.value"]] <- coef(summary(m1))[2,4]
            # paste0(c("Controls model: log(variant_count) = ", " * paternal_age + "), rev(m1$coefficients), collapse="")
            # #plot(control_dat$Paternal.Age.at.Proband.Birth, control_dat$variant_count)
            # #sum(predict(m1, type="response"))
            # predicted_case_variant_count = round(sum(predict(m1, newdata=case_dat, type="response")))
            # m2 <- glm(variant_count ~ Paternal.Age.at.Proband.Birth, family="poisson", data=case_dat)
            # #summary(m2)
            # paste0(c("Cases model: log(variant_count) = ", " * paternal_age + "), rev(m2$coefficients), collapse="")
            # #plot(case_dat$Paternal.Age.at.Proband.Birth, case_dat$variant_count)
            # #sum(predict(m2, type="response"))
            # predicted_control_variant_count = round(sum(predict(m2, newdata=control_dat, type="response")))
            
            #binom_result <- fisher_exact_test(tests$m1[i], predicted_case_variant_count, n1, n1, alternative = c("greater"))
        } else { binom_result <- fisher_exact_test(tests$m1[i], tests$m0[i], n1, n0, alternative="greater") }
        binom_enrichments <- c(binom_enrichments, binom_result[["estimate"]])
        binom_p.values <- c(binom_p.values, binom_result[["p.value"]])
    }
    burdens <- cbind(tests$test, tests$variants, tests$m1, tests$m0, binom_enrichments, binom_p.values, n1, n0, tests$phenotype_burdens, tests$notes, tests$notable_genes, tests$cases_gene_hits, tests$controls_gene_hits, tests$cases_variant_hits, tests$controls_variant_hits)
    colnames(burdens) <- c("burden", "variants", "m1", "m0", "enrichment", "p.value", "n1", "n0", "phenotype_burdens", "notes", "notable_genes", "cases_gene_hits", "controls_gene_hits", "cases_variant_hits", "controls_variant_hits")
    if(statistical_test == "poisson regression") { colnames(burdens)[5] <- "z_score" }
    
    # filter out histone mark tests that are not even marginally sigificance (p < 0.05).
    if (filter_result) {
        rows_to_keep <- unique(c(1:(histone_tests_start_index-1), which(binom_p.values < 0.1)))
        if (histone_tests_end_index > nrow(burdens)) { rows_to_keep <- unique(c(rows_to_keep, (histone_tests_end_index+1):nrow(burdens))) }
        burdens <- burdens[rows_to_keep,]
    }
    # If recurrent, keep only TSS/3'UTR burden tests, since they are the data tables annotated with TS_gene.
    if (recurrent) { burdens <- burdens[grepl("TSS|3'UTR", burdens[,1]),] }
    
    # Return burden results
    return(data.frame(burdens))
}

compare_batch_burdens <- function(burden1, burden2, name1, name2) {
    burden1$burden <- gsub("mouse CHD", "mouse", burden1$burden); burden1$burden <- gsub("mouse CDH", "mouse", burden1$burden)
    burden1$burden <- gsub("HE", "E",burden1$burden); burden1$burden <- gsub("DE", "E", burden1$burden)
    burden2$burden <- gsub("mouse CHD", "mouse", burden2$burden); burden2$burden <- gsub("mouse CDH", "mouse", burden2$burden)
    burden2$burden <- gsub("HE", "E", burden2$burden); burden2$burden <- gsub("DE", "E", burden2$burden)
    
    comparison <- merge(burden1[1:6], burden2[1:6], by=c("burden","variants"));
    comparison$enrichment.x <- as.numeric(paste0(comparison$enrichment.x)); comparison$p.value.x <- as.numeric(paste0(comparison$p.value.x))
    comparison$enrichment.y <- as.numeric(paste0(comparison$enrichment.y)); comparison$p.value.y <- as.numeric(paste0(comparison$p.value.y))
    comparison <- comparison[which((!is.na(comparison$enrichment.x)) & (!is.na(comparison$enrichment.y)) & (comparison$enrichment.x != 0) & (comparison$enrichment.y != 0)),]
    #comparison$enrichment.x[is.na(comparison$enrichment.x)] <- 16; comparison$enrichment.y[is.na(comparison$enrichment.y)] <- 16
    #comparison$enrichment.x[comparison$enrichment.x == Inf] <- 16; comparison$enrichment.y[comparison$enrichment.y == Inf] <- 16
    #comparison$enrichment.x[comparison$enrichment.x == 0] <- 1/8; comparison$enrichment.y[comparison$enrichment.y == 0] <- 1/8
    comparison <- comparison[order(sapply(1:nrow(comparison), function(i) min(c(comparison$p.value.x[i], comparison$p.value.y[i]))), decreasing=TRUE),]
    colors <- rep("gray49", nrow(comparison))
    colors[comparison$p.value.x < 0.05] <- "green"
    colors[comparison$p.value.y < 0.05] <- "blue"
    colors[comparison$p.value.x < 0.05 & comparison$p.value.y < 0.05] <- "red"
    
    filename = output_path(paste0(name1,"_vs_",name2,"_burdens_scatterplot.pdf"))
    pdf(filename)
    plot(log2(comparison$enrichment.x), log2(comparison$enrichment.y), col=colors, main=paste0(name1," vs ",name2," Burden Comparison"), xlab=paste0("log2(",name1," enrichment)"), ylab=paste0("log2(",name2," enrichment)"), cex.axis=1.3, cex.lab=1.3, cex.main=1.2)
    mtext(paste0("Spearman: ", round(cor(log2(comparison$enrichment.x), log2(comparison$enrichment.y), method="spearman"),3), " for all tests, ", round(cor(log2(comparison$enrichment.x[colors!="gray49"]), log2(comparison$enrichment.y[colors!="gray49"]), method="spearman"),3), " for marginally sig., zero counts removed"), cex=1)
    legend("topleft", legend=c(paste0("marginal sig. in ",name1), paste0("marginal sig. in ",name2), paste0("marginal sig. in both")), col=c("green","blue","red"), pch=15, bty="n")
    dev.off()
    print(comparison[colors=="red",])
    return(comparison)
}
compare_batch_burdens(cdh_burden, chd_burden, "CDH1", "CHD1")
compare_batch_burdens(cdh_burden, cdh_combined_burden, "CDH1", "CDH_combined")
write.csv(compare_batch_burdens(cdh_combined_burden, chd_combined_burden, "CDH_all", "CHD_all"), file=output_path("cross_disease_both_marginal.csv"), row.names=FALSE)
compare_batch_burdens(cdh_burden, cdh2_burden, "CDH1", "CDH2")
write.csv(compare_batch_burdens(cdh_combined_burden, cdh_combined_burden_cons, "CDH", "CDH_conserved"), file=output_path("CDH_vs_CDH_conserved_comparison.csv"), row.names=FALSE)
write.csv(compare_batch_burdens(chd_combined_burden, chd_combined_burden_cons, "CHD", "CHD_conserved"), file=output_path("CHD_vs_CHD_conserved_comparison.csv"), row.names=FALSE)
write.csv(compare_batch_burdens(cdh_combined_burden, chd_combined_burden, "CDH", "CHD"), file=output_path("CDH_vs_CHD_comparison.csv"), row.names=FALSE)
write.csv(compare_batch_burdens(cdh_combined_burden_cons, chd_combined_burden_cons, "CDH_conserved", "CHD_conserved"), file=output_path("CDH_conserved_vs_CHD_conserved_comparison.csv"), row.names=FALSE)

#################################################################
# Run burden analysis with volcano plots!
#################################################################
# CDH vs SSC
cdh_burden_hg19 <- basic_burden_analysis("CDH_combined", ssc="SSC_all", annotate_phenotypes=FALSE, statistical_test="fisher_exact")[,-c(12,13)]
cdh_burden_hg19 <- basic_burden_analysis("CDH_combined", ssc="SSC_all", annotate_phenotypes=FALSE, statistical_test="binomial")[,-c(12,13)]
write.table(cdh_burden_hg19, file=output_path("CDH_all_burden_analysis_binomial.tsv"), sep="\t", row.names=FALSE)
best_result_cdh_hg19 <- volcano_plot("All CDH", data.frame(cdh_burden_hg19[-c(1:3),]), strong_effect_cutoff=1, bonferroni_p.value_cutoff=2.19e-4, num_best_hits=10, number_labels=FALSE, x_max=2.2, label_cex=0.75, pval_relative_importance=1000)
# CDH vs SSC Conserved
cdh_burden_hg19_cons <- basic_burden_analysis("CDH_combined", ssc="SSC_all", annotate_phenotypes=FALSE, statistical_test="fisher_exact", conserved_only=TRUE)[,-c(12,13)]
write.table(cdh_burden_hg19_cons, file=output_path("CDH_all_burden_analysis_cons.tsv"), sep="\t", row.names=FALSE)
best_result_cdh_hg19_cons <- volcano_plot("Conserved CDH", cdh_burden_hg19_cons[-c(1:3),], strong_effect_cutoff=1, bonferroni_p.value_cutoff=2.19e-4, num_best_hits=20)

# CHD vs SSC
chd_burden_hg19 <- basic_burden_analysis("CHD_combined", ssc="SSC_all", annotate_phenotypes=FALSE, statistical_test="fisher_exact")[,-c(12,13)]
write.table(chd_burden_hg19, file=output_path("CHD_all_burden_analysis.tsv"), sep="\t", row.names=FALSE)
best_result_chd_hg19 <- volcano_plot("All CHD", chd_burden_hg19[-c(1:3),], strong_effect_cutoff=1, bonferroni_p.value_cutoff=2.19e-4, num_best_hits=20, pval_relative_importance=1000)
# CHD vs SSC Conserved
chd_burden_hg19_cons <- basic_burden_analysis("CHD_combined", ssc="SSC_all", annotate_phenotypes=FALSE, statistical_test="fisher_exact", conserved_only=TRUE)[,-c(12,13)]
write.table(chd_burden_hg19_cons, file=output_path("CHD_all_burden_analysis_cons.tsv"), sep="\t", row.names=FALSE)
best_result_chd_hg19_cons <- volcano_plot("Conserved CHD", chd_burden_hg19_cons[-c(1:3),], strong_effect_cutoff=1, bonferroni_p.value_cutoff=2.19e-4, num_best_hits=20)

filter_dataset("CHDFB_combined", remove_synonymous=FALSE, remove_missense=FALSE)
filter_dataset("SSCFB_combined", remove_synonymous=FALSE, remove_missense=FALSE)
filter_dataset("CHDFB_combined", remove_synonymous=TRUE)
filter_dataset("SSCFB_combined", remove_synonymous=TRUE)
# CHDFB vs SSCFB
chdfb_burden_hg19 <- basic_burden_analysis("CHDFB_combined", ssc="SSCFB_combined", annotate_phenotypes=FALSE, statistical_test="fisher_exact")[,-c(12,13)]
write.table(chdfb_burden_hg19, file=output_path("CHDFB_all_burden_analysis.tsv"), sep="\t", row.names=FALSE)
chdfb_burden_hg19$burden <- sapply(chdfb_burden_hg19$burden, function(x) gsub("tissue_majority", "3tissue", x))
best_result_chdfb_hg19 <- volcano_plot("All CHDFB", chdfb_burden_hg19[-c(1:3),], strong_effect_cutoff=0.9, bonferroni_p.value_cutoff=2.19e-4, num_best_hits=10, number_labels=FALSE, x_max=3.15, label_cex=0.75)
# CHDFB vs SSCFB Conserved
chdfb_burden_hg19_cons <- basic_burden_analysis("CHDFB_combined", ssc="SSCFB_combined", annotate_phenotypes=FALSE, statistical_test="fisher_exact", conserved_only=TRUE)[,-c(12,13)]
write.table(chdfb_burden_hg19_cons, file=output_path("CHDFB_all_burden_analysis_cons.tsv"), sep="\t", row.names=FALSE)
best_result_chdfb_hg19_cons <- volcano_plot("Conserved CHDFB", chdfb_burden_hg19_cons[-c(1:3),], strong_effect_cutoff=1, bonferroni_p.value_cutoff=2.19e-4, num_best_hits=20, requested_labels=best_result_chd_hg19_cons$label[1:10], requested_labels_description="Top 10 in GATK-only analysis")
# CHDFB vs SSCFB Poisson Regression
chdfb_burden_hg19_poisson <- basic_burden_analysis("CHDFB_combined", ssc="SSCFB_combined", annotate_phenotypes=FALSE, statistical_test="poisson regression")[,-c(12,13)]
write.table(chdfb_burden_hg19_poisson, file=output_path("CHDFB_all_burden_analysis_poisson.tsv"), sep="\t", row.names=FALSE)
best_result_chdfb_hg19_poisson <- volcano_plot("All CHDFB", chdfb_burden_hg19_poisson[-c(1:3),], strong_effect_cutoff=0.9, bonferroni_p.value_cutoff=2.19e-4, num_best_hits=10, number_labels=FALSE)

# CDH+CHD vs SSC
cdh_chd_burden_hg19 <- basic_burden_analysis("CDH_CHD", ssc="SSC_all", annotate_phenotypes=FALSE, statistical_test="fisher_exact")[,-c(12,13)]
write.table(cdh_chd_burden_hg19, file=output_path("CDH_CHD_burden_analysis.tsv"), sep="\t", row.names=FALSE)
best_result_cdh_chd <- volcano_plot("All Cases (CDH+CHD)", cdh_chd_burden_hg19[-c(1:3),], strong_effect_cutoff=1, bonferroni_p.value_cutoff=2.19e-4, num_best_hits=20, pval_relative_importance=1000)
# CDH+CHD vs SSC Conserved
cdh_chd_burden_hg19_cons <- basic_burden_analysis("CDH_CHD", ssc="SSC_all", annotate_phenotypes=FALSE, statistical_test="fisher_exact", conserved_only=TRUE)[,-c(12,13)]
write.table(cdh_chd_burden_hg19_cons, file=output_path("CDH_CHD_burden_analysis_cons.tsv"), sep="\t", row.names=FALSE)
best_result_cdh_chd_cons <- volcano_plot("Conserved Cases (CDH+CHD)", cdh_chd_burden_hg19_cons[-c(1:3),], strong_effect_cutoff=1, bonferroni_p.value_cutoff=2.19e-4, num_best_hits=20)
#################################################################


# Comparison of standard vs. conserved burden analyses
a <- merge(cdh_combined_burden[,c(1:2,5)], cdh_combined_burden_cons[,c(1:2,5)], by=c("burden","variants"))
a$enrichment.x <- as.numeric(paste0(a$enrichment.x)); a$enrichment.y <- as.numeric(paste0(a$enrichment.y))
a <- a[order(a$enrichment.y/(a$enrichment.x+0.00001), decreasing=TRUE),]
nrow(a)
boosted_tokens <- unlist(strsplit(paste0(a$burden[a$enrichment.y > a$enrichment.x]), ",|\\)|\\(| ")); boosted_tokens <- table(boosted_tokens[boosted_tokens != ""])
unboosted_tokens <- unlist(strsplit(paste0(a$burden[a$enrichment.y <= a$enrichment.x]), ",|\\)|\\(| ")); unboosted_tokens <- table(unboosted_tokens[boosted_tokens != ""])
boosted_tokens 
unboosted_tokens
boosted_tokens <- cbind(names(boosted_tokens), boosted_tokens); colnames(boosted_tokens) <- c("token", "count"); unboosted_tokens <- cbind(names(unboosted_tokens), unboosted_tokens); colnames(unboosted_tokens) <- c("token", "count")
token_comparison <-merge(boosted_tokens, unboosted_tokens, by="token", all.x=TRUE, all.y=TRUE); token_comparison$count.x <- as.numeric(paste0(token_comparison$count.x)); token_comparison$count.y <- as.numeric(paste0(token_comparison$count.y))
token_comparison$count.x[is.na(token_comparison$count.x)] <- 0; token_comparison$count.y[is.na(token_comparison$count.y)] <- 0
token_comparison <- aggregate(token_comparison[,2:3], by=list(tolower(token_comparison[,1])), FUN=sum)
colnames(token_comparison) <- c("token", "conserved", "standard")
token_comparison <- token_comparison[token_comparison$conserved != token_comparison$standard,]; token_comparison <- token_comparison[order(token_comparison$conserved/token_comparison$standard, decreasing=TRUE),]
token_comparison <- token_comparison[!(token_comparison$token %in% c(">","0.9","10","15","20","3"," ","","line","cell","in","cdh","mouse","genes")),]
colnames(token_comparison) <- c("token", "enriched_in_conserved", "enriched_in_standard")
write.csv(token_comparison, file=output_path("conserved_vs_standard_burden_token_comparison.csv"), row.names=FALSE)

token_comparison <- rep(0, length(unique(c(names(boosted_tokens), names(unboosted_tokens))))); unique(c(names(boosted_tokens), names(unboosted_tokens)))

a1 <- best_results_cdh_518[order(best_results_cdh_518$burden),c("burden","variants","enrichment","p.value")]
a2 <- best_results_cdh2[order(best_results_cdh2$burden),c("burden","variants","enrichment","p.value")]
a_combined <- best_result_cdh_all[order(best_result_cdh_all$burden),c("burden","variants","enrichment","p.value")]
enrichment_comparison <- merge(a1, a2, by=c("burden","variants"))
enrichment_comparison <- merge(enrichment_comparison, a_combined, by=c("burden","variants"))
colnames(enrichment_comparison) <- c("test", "variants", "CDH1_SSC518_enrichment", "CDH1_SSC518_p.value", "CDH2_SSC1088_enrichment", "CDH2_SSC1088_p.value", "CDH_all_SSC_all_enrichment", "CDH_all_SSC_all_p.value")
enrichment_comparison <- enrichment_comparison[enrichment_comparison$CDH_all_SSC_all_p.value < 0.05,]
pdf(file=output_path("cdh_enrichment_comparison_with_cdh2.pdf"))
scatterplot(log10(enrichment_comparison$CDH1_SSC518_enrichment), log10(enrichment_comparison$CDH2_SSC1088_enrichment), main="CDH Enrichment Comparison", xlab="CDH195 vs. SSC518 log10(enrichment)", ylab="CDH294 vs. SSC1088 log10(enrichment)", cex.main=1.5, cex.axis=1.4, cex.lab=1.4, spread=FALSE, smoother=FALSE)
legend("topright", legend=paste0("Spearman = ", round(cor(enrichment_comparison$CDH1_SSC518_enrichment, enrichment_comparison$CDH2_SSC1088_enrichment, method="spearman"),3)), col=c("white"), lty=c(1), bty="n")
mtext(paste0(sum(enrichment_comparison$CDH_all_SSC_all_p.value < 0.05), " tests with at least marginal significance in combined analysis"), padj=-1.05, cex=1)
dev.off()
pdf(file=output_path("cdh_enrichment_comparison_with_cdh_all.pdf"))
scatterplot(log10(enrichment_comparison$CDH1_SSC518_enrichment), log10(enrichment_comparison$CDH_all_SSC_all_enrichment), main="CDH Enrichment Comparison", xlab="CDH195 vs. SSC518 log10(enrichment)", ylab="CDH_all vs. SSC_all log10(enrichment)", cex.main=1.5, cex.axis=1.4, cex.lab=1.4, spread=FALSE, smoother=FALSE)
legend("topright", legend=paste0("Spearman = ", round(cor(enrichment_comparison$CDH1_SSC518_enrichment, enrichment_comparison$CDH_all_SSC_all_enrichment, method="spearman"),3)), col=c("white"), lty=c(1), bty="n")
mtext(paste0(sum(enrichment_comparison$CDH_all_SSC_all_p.value < 0.05), " tests with at least marginal significance in combined analysis"), padj=-1.05, cex=1)
dev.off()
pdf(file=output_path("cdh2_enrichment_comparison_with_cdh_all.pdf"))
scatterplot(log10(enrichment_comparison$CDH2_SSC1088_enrichment), log10(enrichment_comparison$CDH_all_SSC_all_enrichment), main="CDH Enrichment Comparison", xlab="CDH294 vs. SSC1088 log10(enrichment)", ylab="CDH_all vs. SSC_all log10(enrichment)", cex.main=1.5, cex.axis=1.4, cex.lab=1.4, spread=FALSE, smoother=FALSE)
legend("topright", legend=paste0("Spearman = ", round(cor(enrichment_comparison$CDH2_SSC1088_enrichment, enrichment_comparison$CDH_all_SSC_all_enrichment, method="spearman"),3)), col=c("white"), lty=c(1), bty="n")
mtext(paste0(sum(enrichment_comparison$CDH_all_SSC_all_p.value < 0.05), " tests with at least marginal significance in combined analysis"), padj=-1.05, cex=1)
dev.off()

#mtext("200 random genes (675.2 Kbp, exonic); model samples every 30th trimer", padj=-1.05, cex=1)
legend("topright", legend=c("Smoothed Line", "Regression Line", paste0("Spearman = ", round(cor(a$model_rate, a$hongjian_rate, method="spearman"),3))), col=c("red","green","white"), lty=c(1,1), bty="n")


# Burden of recurrently mutated genes, to hopefully improve signal-to-noise.
cdh_burden_recurrent <- basic_burden_analysis("CDH", recurrent=TRUE)[,-c(11,12)]
write.table(cdh_burden_recurrent, file=output_path("CDH_recurrent_basic_burden_analysis.tsv"), sep="\t", row.names=FALSE)
chd_burden_recurrent <- basic_burden_analysis("CHD", recurrent=TRUE)[,-c(11,12)]
write.table(chd_burden_recurrent, file=output_path("CHD_recurrent_basic_burden_analysis.tsv"), sep="\t", row.names=FALSE)
cdh_new_burden_recurrent <- basic_burden_analysis("CDH_new", recurrent=TRUE)[,-c(11,12)]
write.table(cdh_new_burden_recurrent, file=output_path("CDH2_SSC2_recurrent_burden_analysis.tsv"), sep="\t", row.names=FALSE)



# Write variants to file.
if (TRUE) {
    # CDH
    variants <- cdh_muts_TES[,c("X.chr","pos","ref","alt","sample","TS_gene","RBP")]; colnames(variants) <- c("chromosome","pos","ref","alt","sample","gene","RBP")
    snv_indel <- rep("snv", nrow(variants)); snv_indel[is_indel(variants$ref, variants$alt)] <- "indel"; variants <- cbind(variants, snv_indel)
    write.table(variants, file=output_path("CDH_3UTR_variants.tsv"), sep="\t", row.names=FALSE)
    variants <- cdh_muts_TSS[,c("X.chr","pos","ref","alt","sample","TS_gene","RBP")]; colnames(variants) <- c("chromosome","pos","ref","alt","sample","gene","RBP")
    snv_indel <- rep("snv", nrow(variants)); snv_indel[is_indel(variants$ref, variants$alt)] <- "indel"; variants <- cbind(variants, snv_indel)
    write.table(variants, file=output_path("CDH_TSS_variants.tsv"), sep="\t", row.names=FALSE)
    #CHD
    variants <- chd_muts_TES[,c("X.chr","pos","ref","alt","sample","TS_gene","RBP")]; colnames(variants) <- c("chromosome","pos","ref","alt","sample","gene","RBP")
    snv_indel <- rep("snv", nrow(variants)); snv_indel[is_indel(variants$ref, variants$alt)] <- "indel"; variants <- cbind(variants, snv_indel)
    write.table(variants, file=output_path("CHD_3UTR_variants.tsv"), sep="\t", row.names=FALSE)
    variants <- chd_muts_TSS[,c("X.chr","pos","ref","alt","sample","TS_gene","RBP")]; colnames(variants) <- c("chromosome","pos","ref","alt","sample","gene","RBP")
    snv_indel <- rep("snv", nrow(variants)); snv_indel[is_indel(variants$ref, variants$alt)] <- "indel"; variants <- cbind(variants, snv_indel)
    write.table(variants, file=output_path("CHD_TSS_variants.tsv"), sep="\t", row.names=FALSE)
    #CDH2
    variants <- cdh_new_muts_TES[,c("X.chr","pos","ref","alt","sample","TS_gene","RBP")]; colnames(variants) <- c("chromosome","pos","ref","alt","sample","gene","RBP")
    snv_indel <- rep("snv", nrow(variants)); snv_indel[is_indel(variants$ref, variants$alt)] <- "indel"; variants <- cbind(variants, snv_indel)
    write.table(variants, file=output_path("CDH2_3UTR_variants.tsv"), sep="\t", row.names=FALSE)
    variants <- cdh_new_muts_TSS[,c("X.chr","pos","ref","alt","sample","TS_gene","RBP")]; colnames(variants) <- c("chromosome","pos","ref","alt","sample","gene","RBP")
    snv_indel <- rep("snv", nrow(variants)); snv_indel[is_indel(variants$ref, variants$alt)] <- "indel"; variants <- cbind(variants, snv_indel)
    write.table(variants, file=output_path("CDH2_TSS_variants.tsv"), sep="\t", row.names=FALSE)
}

# Load TAD regions data and make TAD Granges object.
TADs <- read.csv(data_path("TAD H1\\combined\\hg19_total.combined.domain"), sep="\t", header=FALSE)
TADs[,1]<- gsub("chr", "", TADs[,1])
TADs <- cbind(paste0(TADs[,1], ":", TADs[,2], "-", TADs[,3]), TADs)
colnames(TADs) <- c("name", "chromosome", "start", "end")
TAD_granges <- to_genomic_regions(TADs, label_colname="name")

# Returns which TADs (topologically-associated domains) the given regions (GRanges) fall in; in case a region intersects multiple TADs, choose the one with most overlap.
get_TADs <- function(region_granges) {
    TAD_results <- c()
    all_hits <- data.frame(unique(olRanges(region_granges, TAD_granges)))
    for(i in 1:length(region_granges)) {
        hits <- all_hits[all_hits$Qindex == i,]
        if(nrow(hits) > 0) {
            best_hit <- names(TAD_granges)[hits$Sindex[which.max(hits$OLpercQ)]]
        } else { best_hit <- "" }
        TAD_results <- c(TAD_results, best_hit)
    }
    return(TAD_results)
}

# Return clusterings for the given regions (GRanges), based on overlaps. Can optionally break up by TADs, to avoid nonsensical long chains forming.
# Requires igraph package
cluster_overlapping_regions <- function(region_granges, region_name="", break_by_TAD=FALSE, make_plots=TRUE) {
    regions <- cbind(names(region_granges), data.frame(region_granges)); colnames(regions)[1] <- "gene"
    hits <- data.frame(findOverlaps(region_granges, region_granges))
    gene_pairs <- unique(rbind(t(apply(hits[hits$queryHits != hits$subjectHits,], 1, function(x) { names(region_granges)[sort(x)] } )), cbind(names(region_granges),names(region_granges))))
    if(break_by_TAD) {
        print("Annotating regions with TADs...")
        gene_TADs <- new.env()
        region_TADs_vector <- get_TADs(region_granges)
        for(i in 1:nrow(regions)) {
            print(i)
            gene_TADs[[paste(regions$gene[i])]] <- region_TADs_vector[i]
        }
        same_TAD <- rep(TRUE, nrow(gene_pairs))
        print("Breaking gene pairs by TAD...")
        for(i in 1:nrow(gene_pairs)) {
            print(i)
            gene_pair <- gene_pairs[i,]
            same <- (gene_TADs[[gene_pair[1]]] == gene_TADs[[gene_pair[2]]])
            if(!same) { same_TAD[i] <- FALSE }
        }
        gene_pairs <- gene_pairs[same_TAD,]
    }
    g <- graph.data.frame(gene_pairs, directed=FALSE)
    comps <- components(g)
    num_comps <- comps$no
    if (make_plots) {
        if (break_by_TAD) { tad_string = "_TAD" } else { tad_string = "" }
        pdf(output_path(paste0("clustered_", gsub('[^0-9a-zA-Z]', '', region_name), tad_string, "_region_size_distribution.pdf")))
        hist(comps$csize, breaks=seq(0.5, (max(comps$csize)+0.5)), main=paste0("Clustered ", region_name, " Region Sizes"), xlab="# genes", cex.main=1.5, cex.axis=1.3, cex.lab=1.4, lwd=2, col="blue")
        if (break_by_TAD) { tad_string = "TAD, " } else { tad_string = "" }
        mtext(paste0(tad_string, length(V(g)), " genes, ", num_comps, " comps, largest with ", max(comps$csize), " genes"), cex=1.2)
        dev.off()
    }
    
    clustered_regions <- data.frame(stringsAsFactors=FALSE)
    print("Creating clustered regions...")
    for(comp in 1:num_comps) {
        print(comp)
        connected_genes <- paste0(names(V(g)[which(comps$membership == comp)]))
        connected_genes_string <- paste0(sort(connected_genes), collapse=",")
        comp_regions <- regions[regions$gene %in% connected_genes,]
        clustered_regions <- rbind(clustered_regions, cbind(connected_genes_string, paste(comp_regions$seqnames[1]), min(comp_regions$start), max(comp_regions$end)))
    }
    colnames(clustered_regions) <- c("gene", "chromosome", "start", "end")
    #clustered_regions$start <- as.numeric(paste0(clustered_regions$start)); clustered_regions$end <- as.numeric(paste0(clustered_regions$end))
    #colnames(clustered_regions) <- c("gene", "chromosome", "start", "end")
    clustered_region_granges <- c(to_genomic_regions(clustered_regions[!(clustered_regions$gene %in% names(region_granges)),]), region_granges[names(region_granges) %in% clustered_regions$gene])
    
    if (make_plots) {
        if (break_by_TAD) { tad_string = "_TAD" } else { tad_string = "" }
        pdf(output_path(paste0("clustered_", gsub('[^0-9a-zA-Z]', '', region_name), tad_string, "_region_length_distribution.pdf")))
        plot(density(log10(end(clustered_region_granges)-start(clustered_region_granges)), from=log10(min(end(clustered_region_granges)-start(clustered_region_granges)))), main=paste0("Clustered ", region_name, " Region Lengths"), xlab="log10(bp)", cex.main=1.5, cex.axis=1.3, cex.lab=1.4, lwd=3, col="blue", xaxs="i")
        if (break_by_TAD) { tad_string = "TAD, " } else { tad_string = "" }
        mtext(paste0(tad_string, length(V(g)), " genes, ", num_comps, " comps, lengths in [", min(end(clustered_region_granges)-start(clustered_region_granges)), " bp, ", max(end(clustered_region_granges)-start(clustered_region_granges)), " bp]"), cex=1.2)
        dev.off()
    }
    
    return(clustered_region_granges)
}
clustered_TES_TAD_granges <- cluster_overlapping_regions(TES_granges, region_name="3'UTR", break_by_TAD=TRUE)
clustered_TSS_TAD_granges <- cluster_overlapping_regions(TSS_granges, region_name="TSS", break_by_TAD=TRUE)
clustered_TES_granges <- cluster_overlapping_regions(TES_granges, region_name="3'UTR")
clustered_TSS_granges <- cluster_overlapping_regions(TSS_granges, region_name="TSS")

tad_string = "_TAD"
clustered_region_granges <- clustered_TES_TAD_granges; region_name = "3'UTR"
pdf(output_path(paste0("clustered_", gsub('[^0-9a-zA-Z]', '', region_name), tad_string, "_region_length_distribution.pdf")))
#plot(density(log10(end(clustered_region_granges)-start(clustered_region_granges)), from=log10(min(end(clustered_region_granges)-start(clustered_region_granges)))), main=paste0("Clustered ", region_name, " Region Lengths"), xlab="log10(bp)", xlim=c(0, log10(max(end(clustered_region_granges)-start(clustered_region_granges))*10)), cex.main=1.5, cex.axis=1.3, cex.lab=1.4, lwd=3, col="blue", xaxs="i")
plot(density(end(clustered_region_granges)-start(clustered_region_granges), from=min(end(clustered_region_granges)-start(clustered_region_granges))), main=paste0("Clustered ", region_name, " Region Lengths"), xlab="bp", xlim=c(0, max(end(clustered_region_granges)-start(clustered_region_granges))*1.1), cex.main=1.5, cex.axis=1.3, cex.lab=1.4, lwd=3, col="blue", xaxs="i")
tad_string <- gsub("_TAD", "TAD, ", tad_string)
mtext(paste0(tad_string, length(V(g)), " genes, ", length(clustered_region_granges), " comps, lengths in [", min(end(clustered_region_granges)-start(clustered_region_granges)), " bp, ", max(end(clustered_region_granges)-start(clustered_region_granges)), " bp]"), cex=1.2)
dev.off()
tad_string = "_TAD"
clustered_region_granges <- clustered_TSS_TAD_granges; region_name = "TSS"
pdf(output_path(paste0("clustered_", gsub('[^0-9a-zA-Z]', '', region_name), tad_string, "_region_length_distribution.pdf")))
plot(density(end(clustered_region_granges)-start(clustered_region_granges), from=min(end(clustered_region_granges)-start(clustered_region_granges))), main=paste0("Clustered ", region_name, " Region Lengths"), xlab="bp", xlim=c(0, max(end(clustered_region_granges)-start(clustered_region_granges))*1.1), cex.main=1.5, cex.axis=1.3, cex.lab=1.4, lwd=3, col="blue", xaxs="i")
tad_string <- gsub("_TAD", "TAD, ", tad_string)
mtext(paste0(tad_string, length(V(g)), " genes, ", length(clustered_region_granges), " comps, lengths in [", min(end(clustered_region_granges)-start(clustered_region_granges)), " bp, ", max(end(clustered_region_granges)-start(clustered_region_granges)), " bp]"), cex=1.2)
dev.off()
tad_string = ""
clustered_region_granges <- clustered_TES_granges; region_name = "3'UTR"
pdf(output_path(paste0("clustered_", gsub('[^0-9a-zA-Z]', '', region_name), tad_string, "_region_length_distribution.pdf")))
#plot(density(log10(end(clustered_region_granges)-start(clustered_region_granges)), from=log10(min(end(clustered_region_granges)-start(clustered_region_granges)))), main=paste0("Clustered ", region_name, " Region Lengths"), xlab="log10(bp)", xlim=c(0, log10(max(end(clustered_region_granges)-start(clustered_region_granges))*10)), cex.main=1.5, cex.axis=1.3, cex.lab=1.4, lwd=3, col="blue", xaxs="i")
plot(density(end(clustered_region_granges)-start(clustered_region_granges), from=min(end(clustered_region_granges)-start(clustered_region_granges))), main=paste0("Clustered ", region_name, " Region Lengths"), xlab="bp", xlim=c(0, max(end(clustered_region_granges)-start(clustered_region_granges))*1.1), cex.main=1.5, cex.axis=1.3, cex.lab=1.4, lwd=3, col="blue", xaxs="i")
tad_string <- gsub("_TAD", "TAD, ", tad_string)
mtext(paste0(tad_string, length(V(g)), " genes, ", length(clustered_region_granges), " comps, lengths in [", min(end(clustered_region_granges)-start(clustered_region_granges)), " bp, ", max(end(clustered_region_granges)-start(clustered_region_granges)), " bp]"), cex=1.2)
dev.off()
tad_string = ""
clustered_region_granges <- clustered_TSS_granges; region_name = "TSS"
pdf(output_path(paste0("clustered_", gsub('[^0-9a-zA-Z]', '', region_name), tad_string, "_region_length_distribution.pdf")))
plot(density(end(clustered_region_granges)-start(clustered_region_granges), from=min(end(clustered_region_granges)-start(clustered_region_granges))), main=paste0("Clustered ", region_name, " Region Lengths"), xlab="bp", xlim=c(0, max(end(clustered_region_granges)-start(clustered_region_granges))*1.1), cex.main=1.5, cex.axis=1.3, cex.lab=1.4, lwd=3, col="blue", xaxs="i")
tad_string <- gsub("_TAD", "TAD, ", tad_string)
mtext(paste0(tad_string, length(V(g)), " genes, ", length(clustered_region_granges), " comps, lengths in [", min(end(clustered_region_granges)-start(clustered_region_granges)), " bp, ", max(end(clustered_region_granges)-start(clustered_region_granges)), " bp]"), cex=1.2)
dev.off()

# Convert normal gene annotation of variants to the clusters that each represents. 
convert_genes_to_clusters <- function(dat, sites, TAD=FALSE) {
    genes <- paste0(dat$TS_gene)
    if(TAD) { TAD_string <- "TAD" } else { TAD_string <- NULL }
    cluster_granges_name = paste(c("clustered", sites, TAD_string, "granges"), collapse="_")
    cluster_granges <- try(get(cluster_granges_name))
    if(inherits(cluster_granges, "try-error")) { 
        region_granges <- get(paste0(c(sites, "_granges")))
        cluster_granges <- cluster_overlapping_regions(region_granges, region_name=sites, break_by_TAD=TAD, make_plots=FALSE)
    }
    
    print("Hashing gene-to-cluster mappings...")
    gene_cluster_mapping <- new.env()
    for(i in 1:length(cluster_granges)) {
        print(i)
        cluster <- paste0(names(cluster_granges)[i])
        all_genes_in_cluster <- unlist(strsplit(cluster, ","))
        for(cluster_gene in all_genes_in_cluster) { gene_cluster_mapping[[cluster_gene]] <- cluster }
    }
    
    print("Assigning genes to clusters...")
    clusters <- c()
    for(i in 1:length(genes)) {
        print(i)
        optimum_cluster <- NULL
        all_genes <- unlist(strsplit(paste0(genes[i]), ","))
        mapped_clusters <- unique(unlist(lapply(all_genes, function(x) { gene_cluster_mapping[[x]] } )))
        mapped_clusters <- mapped_clusters[order(start(cluster_granges[mapped_clusters]))] # order by starting position, to process from left to right.
        # Break any possible ties by using position and centralness of it in each cluster. Thus, each position will only map to 1 cluster, reproduceably.
        if (length(mapped_clusters) > 1) {  
            centralness_scores <- c() # Vector storing how central the mutation is in each mapped cluster.
            position = dat$pos[i]
            for(mapped_cluster_index in 1:length(mapped_clusters)) {
                mapped_cluster <- mapped_clusters[mapped_cluster_index]
                cluster_start = start(cluster_granges[mapped_cluster])
                cluster_end = end(cluster_granges[mapped_cluster])
                cluster_length = cluster_end - cluster_start
                if(mapped_cluster_index == 1) { cluster_edge_position = cluster_end # first (upstream) cluster, so we concerned with end
                } else { cluster_edge_position = cluster_start } # second (downstream) cluster, so we concerned with start
                centralness_score = abs(cluster_edge_position - position) / cluster_length
                centralness_scores <- c(centralness_scores, centralness_score)
            }
            optimum_cluster = mapped_clusters[which.max(centralness_scores)] # in case of tie, will pick leftmost (most upstream) cluster
        } else { optimum_cluster = mapped_clusters[1] }
        
        if(length(optimum_cluster) < 1) { optimum_cluster = "" } # No cluster could be mapped.
        clusters <- c(clusters, optimum_cluster)
    }
    
    return(clusters)
}
#a <- convert_genes_to_clusters(chd_indels_TES, "TES", TAD=TRUE)

# Run new version of recurrence analysis! If fix_variant_count is false (default), fixes sample count instead.
recurrence_analysis2 <- function(disease, sites, cases, controls, filename, main, num_iterations=25000, n_greater_than=1, case_samples=NULL, control_samples=NULL, TAD=TRUE, plot_sim=FALSE, plot_control_line=FALSE, fix_variant_count=FALSE) {
    all_regions <- unique(rbind(cases[,which(colnames(cases) %in% c("sample","TS_gene"))], controls[,which(colnames(controls) %in% c("sample","TS_gene"))]))
    if (is.null(case_samples)) { case_samples <- unique(paste0(cases$sample)) }
    if (is.null(control_samples)) { control_samples <- unique(paste0(controls$sample)) } 
    case_samples <- paste0(case_samples); control_samples <- paste0(control_samples)
    all_samples <- c(case_samples, control_samples)
    num_case_variants = nrow(cases); num_cases = length(case_samples); num_controls = length(control_samples); total_samples <- num_cases + num_controls
    
    #all_regions$TS_gene <- cluster_overlapping_genes(all_regions$TS_gene)
    all_regions$TS_gene <- convert_genes_to_clusters(all_regions, sites, TAD)
    all_regions <- all_regions[all_regions$TS_gene != "",]
    all_regions <- all_regions[order(all_regions$sample),]
    all_regions <- unique(all_regions)
    cases <- all_regions[all_regions$sample %in% case_samples,]
    controls <- all_regions[all_regions$sample %in% control_samples,]
    
    recurrence_counts <- c()
    for(i in 1:num_iterations) {
        print(i)
        if(fix_variant_count) {
            permuted_case_genes <- paste0(all_regions$TS_gene[sample(1:nrow(all_regions), num_case_variants)])
        } else {
            permuted_case_samples <- sample(all_samples, num_cases)
            permuted_case_genes <- paste0(all_regions$TS_gene[all_regions$sample %in% permuted_case_samples])
        }
        recurrence_count = length(unique(permuted_case_genes[duplicated(permuted_case_genes)]))
        recurrence_counts <- c(recurrence_counts, recurrence_count)
    }
    
    case_genes <- paste0(cases$TS_gene)
    original_recurrence_count <- length(unique(case_genes[duplicated(case_genes)]))
    control_genes1 <- paste0(controls$TS_gene[controls$sample %in% control_samples[1:num_cases]])
    control_recurrence_count1 <- length(unique(control_genes1[duplicated(control_genes1)]))
    control_genes2 <- paste0(controls$TS_gene[controls$sample %in% control_samples[(length(control_samples)-num_cases+1):length(control_samples)]])
    control_recurrence_count2 <- length(unique(control_genes2[duplicated(control_genes2)]))
    print(paste0(c(original_recurrence_count, mean(recurrence_counts), control_recurrence_count1, control_recurrence_count2), collapse=", "))
    
    # Plot
    cases_line <- original_recurrence_count
    cases_permutations <- recurrence_counts
    
    lines <- new.env()
    cols <- list()
    ltys <- list()
    mtext_label <- ""
    mtext_label <- paste0(mtext_label, num_cases, " case samples")
    
    # Cases line (must be plotted)
    label <- paste0(disease, " data (value=", original_recurrence_count, ")")
    lines[[label]] <- paste0("v=", cases_line)
    cols[[label]] <- "red"
    
    # Controls line
    if (plot_control_line) {
        label <- paste0("SSC data (value=", mean(c(control_recurrence_count1, control_recurrence_count2)), ")")
        lines[[label]] <- paste0("v=", mean(c(control_recurrence_count1, control_recurrence_count2)))
        cols[[label]] <- "black"
        ltys[[label]] <- 3
    }
    
    # Permutation
    label <- paste0("Permutations (N=", prettyNum(num_iterations, big.mark=",", scientific=FALSE), ")")
    lines[[label]] <- density(cases_permutations, from=0, bw=2)
    cols[[label]] <- "blue"
    perm_p_val = sum(cases_permutations > cases_line)/length(cases_permutations)
    if (perm_p_val == 0) { perm_p_val <- paste0("< ", prettyNum(1/num_iterations, scientific=TRUE)) } else { perm_p_val <- round(perm_p_val, 5) }
    mtext_label <- paste0(mtext_label, ", Perm. p-value: ", perm_p_val)
    
    # Background expectation for clustered regions.
    if(plot_sim) {
        num_sims = 100
        label <- paste0("Background Sim (N=",num_sims,")")
        #lines[[label]] <- density(background_for_clusters(sites, num_cases=num_cases, num_simulations=num_sims), from=0, bw=2)
        lines[[label]] <- density(get(paste0(tolower(sites),"_background_sim_for_clusters")), from=0, bw=2)
        cols[[label]] <- "green"
    }
    
    multi_plot(lines=lines, cols=cols, ltys=ltys, lwd=2, cex.lab=1.8, cex.axis=2, cex.main=1.5, mtext_cex=1.5, mar=1.1*par("mar"), main=main, mtext=mtext_label, xlab=paste0("n regions mutated in >", n_greater_than, " sample"), ylab="Density", legend_location="topleft", legend_cex=1.3, xaxs="i", yaxs="i", file=filename)
    
    return(recurrence_counts)
}

background_for_clusters <- function(sites, num_cases, num_simulations=100) {
    region_granges <- get(paste0("clustered_", sites, "_TAD_granges"))
    trimer_counts_folder = output_path(paste0("trimer_counts_for_clustered_",sites,"_h9"))
    background_expectation_per_individual <- sapply(1:length(region_granges), function(i) {
        print(paste0(i, " / ", length(region_granges)))
        region_grange <- region_granges[i]
        chromosome = paste0(seqnames(region_grange)); start = start(region_grange); end = end(region_grange)
        tri_counts <- read.csv(file=full_path(trimer_counts_folder, paste0("chr",chromosome,"_",start,"-",end,"_trimer_counts.csv")))
        tri_counts_names <- names(tri_counts); tri_counts <- as.numeric(tri_counts); names(tri_counts) <- tri_counts_names
        tri_mut_rates <- get(paste0(tolower(sites),"_trimer_mutation_rates_from_ssc_cpg_meth_aggregate"))
        expected_muts = 0
        for(tri in names(tri_counts[tri_counts > 0])) {
            expected_muts = expected_muts + as.numeric(tri_counts[tri]) * sum(tri_mut_rates[grepl(paste0(tri,"->"),names(tri_mut_rates))])
        }
        return(expected_muts)
    })
    total_mutations <- 0
    recurrence_sim <- sapply(1:num_simulations, function(i) { print(i); recurrence_count = sum(sapply(background_expectation_per_individual, function(lamb) { sim <- rpois(num_cases, lambda=lamb); total_mutations <<- total_mutations + sum(sim); return(sum(sim>=1)>1) })); return(recurrence_count) })
    print(paste0(total_mutations, " total mutations simulated over ", num_simulations, " simulations."))
    return(recurrence_sim)
}
haha <- cbind(names(region_granges), (end(region_granges)-start(region_granges)), background_expectation_per_individual)
colnames(haha) <- c("cluster", "length(bp)", "background_expectation_per_sample")
haha[sample(1:nrow(haha)),][1:10,]
haha <- haha[order(as.numeric(haha[,3])/as.numeric(haha[,2]), decreasing=TRUE),]

tes_background_sim_for_clusters <- background_for_clusters("TES", num_cases=750, num_simulations=100)
tss_background_sim_for_clusters <- background_for_clusters("TSS", num_cases=750, num_simulations=100)

recurrence_counts <- recurrence_analysis2("CHD", "TES", chd_indels_TES, sscfb_indels_TES, output_path("recurrence_CHD_selected_indels_TES.pdf"), main=paste0("CHD 3'UTR Region Indel Recurrence Analysis"))
recurrence_counts <- recurrence_analysis2("CHD", "TSS", chd_indels_TSS, sscfb_indels_TSS, output_path("recurrence_CHD_selected_indels_TSS.pdf"), main=paste0("CHD TSS Region Indel Recurrence Analysis"))
recurrence_counts <- recurrence_analysis2("CHD", "TES", chd_indels_TES, sscfb_indels_TES, case_samples=unique(chd$sample), control_samples=unique(sscfb$sample), output_path("recurrence_CHD_indels_TES.pdf"), main=paste0("CHD 3'UTR Region Indel Recurrence Analysis"))
recurrence_counts <- recurrence_analysis2("CHD", "TSS", chd_indels_TSS, sscfb_indels_TSS, case_samples=unique(chd$sample), control_samples=unique(sscfb$sample), output_path("recurrence_CHD_indels_TSS.pdf"), main=paste0("CHD TSS Region Indel Recurrence Analysis"))
recurrence_counts <- recurrence_analysis2("CHD", "TES", chd_snps_TES, sscfb_snps_TES, case_samples=unique(chd$sample), control_samples=unique(sscfb$sample), output_path("recurrence_CHD_snps_TES.pdf"), main=paste0("CHD 3'UTR Region SNV Recurrence Analysis"))
recurrence_counts <- recurrence_analysis2("CHD", "TSS", chd_snps_TSS, sscfb_snps_TSS, case_samples=unique(chd$sample), control_samples=unique(sscfb$sample), output_path("recurrence_CHD_snps_TSS.pdf"), main=paste0("CHD TSS Region SNV Recurrence Analysis"))
recurrence_counts <- recurrence_analysis2("CHD", "TES", chd_snps_TES[includes_at_least_one(chd_snps_TES$TS_gene, HHE_genes),], sscfb_snps_TES[includes_at_least_one(sscfb_snps_TES$TS_gene, HHE_genes),], case_samples=unique(chd$sample), control_samples=unique(sscfb$sample), output_path("recurrence_CHD_snps_HHE_TES.pdf"), main=paste0("CHD HHE 3'UTR Region SNV Recurrence Analysis"))
recurrence_counts <- recurrence_analysis2("CHD", "TSS", chd_snps_TSS[includes_at_least_one(chd_snps_TSS$TS_gene, HHE_genes),], sscfb_snps_TSS[includes_at_least_one(sscfb_snps_TSS$TS_gene, HHE_genes),], case_samples=unique(chd$sample), control_samples=unique(sscfb$sample), output_path("recurrence_CHD_snps_HHE_TSS.pdf"), main=paste0("CHD HHE TSS Region SNV Recurrence Analysis"))
recurrence_counts <- recurrence_analysis2("CHD", "TES", chd_muts_TES, sscfb_muts_TES, case_samples=unique(chd$sample), control_samples=unique(sscfb$sample), output_path("recurrence_CHD_muts_TES.pdf"), main=paste0("CHD 3'UTR Region SNV+Indel Recurrence Analysis"))
recurrence_counts <- recurrence_analysis2("CHD", "TSS", chd_muts_TSS, sscfb_muts_TSS, case_samples=unique(chd$sample), control_samples=unique(sscfb$sample), output_path("recurrence_CHD_muts_TSS.pdf"), main=paste0("CHD TSS Region SNV+Indel Recurrence Analysis"))
recurrence_counts <- recurrence_analysis2("CHD", "TES", chd_muts_TES[includes_at_least_one(chd_muts_TES$TS_gene, HHE_genes),], sscfb_muts_TES[includes_at_least_one(sscfb_muts_TES$TS_gene, HHE_genes),], case_samples=unique(chd$sample), control_samples=unique(sscfb$sample), output_path("recurrence_CHD_muts_HHE_TES.pdf"), main=paste0("CHD HHE 3'UTR Region SNV+Indel Recurrence Analysis"))
recurrence_counts <- recurrence_analysis2("CHD", "TSS", chd_muts_TSS[includes_at_least_one(chd_muts_TSS$TS_gene, HHE_genes),], sscfb_muts_TSS[includes_at_least_one(sscfb_muts_TSS$TS_gene, HHE_genes),], case_samples=unique(chd$sample), control_samples=unique(sscfb$sample), output_path("recurrence_CHD_muts_HHE_TSS.pdf"), main=paste0("CHD HHE TSS Region SNV+Indel Recurrence Analysis"))

recurrence_counts <- recurrence_analysis2("CDH", "TES", cdh_indels_TES, ssc_indels_TES, case_samples=unique(cdh$sample), control_samples=unique(ssc$sample), output_path("recurrence_CDH_indels_TES.pdf"), main=paste0("CDH 3'UTR Region Indel Recurrence Analysis"))
recurrence_counts <- recurrence_analysis2("CDH", "TSS", cdh_indels_TSS, ssc_indels_TSS, case_samples=unique(cdh$sample), control_samples=unique(ssc$sample), output_path("recurrence_CDH_indels_TSS.pdf"), main=paste0("CDH TSS Region Indel Recurrence Analysis"))
recurrence_counts <- recurrence_analysis2("CDH", "TES", cdh_snps_TES, ssc_snps_TES, case_samples=unique(cdh$sample), control_samples=unique(ssc$sample), output_path("recurrence_CDH_snps_TES.pdf"), main=paste0("CDH 3'UTR Region SNV Recurrence Analysis"))
recurrence_counts <- recurrence_analysis2("CDH", "TSS", cdh_snps_TSS, ssc_snps_TSS, case_samples=unique(cdh$sample), control_samples=unique(ssc$sample), output_path("recurrence_CDH_snps_TSS.pdf"), main=paste0("CDH TSS Region SNV Recurrence Analysis"))
recurrence_counts <- recurrence_analysis2("CDH", "TES", cdh_snps_TES[includes_at_least_one(cdh_snps_TES$TS_gene, HDE_genes),], ssc_snps_TES[includes_at_least_one(ssc_snps_TES$TS_gene, HDE_genes),], case_samples=unique(cdh$sample), control_samples=unique(ssc$sample), output_path("recurrence_CDH_snps_HDE_TES.pdf"), main=paste0("CDH HDE 3'UTR Region SNV Recurrence Analysis"), num_iterations=25000)
recurrence_counts <- recurrence_analysis2("CDH", "TSS", cdh_snps_TSS[includes_at_least_one(cdh_snps_TSS$TS_gene, HDE_genes),], ssc_snps_TSS[includes_at_least_one(ssc_snps_TSS$TS_gene, HDE_genes),], case_samples=unique(cdh$sample), control_samples=unique(ssc$sample), output_path("recurrence_CDH_snps_HDE_TSS.pdf"), main=paste0("CDH HDE TSS Region SNV Recurrence Analysis"), num_iterations=25000)
recurrence_counts <- recurrence_analysis2("CDH", "TES", cdh_muts_TES, ssc_muts_TES, case_samples=unique(cdh$sample), control_samples=unique(ssc$sample), output_path("recurrence_CDH_muts_TES.pdf"), main=paste0("CDH 3'UTR Region SNV+Indel Recurrence Analysis"))
recurrence_counts <- recurrence_analysis2("CDH", "TSS", cdh_muts_TSS, ssc_muts_TSS, case_samples=unique(cdh$sample), control_samples=unique(ssc$sample), output_path("recurrence_CDH_muts_TSS.pdf"), main=paste0("CDH TSS Region SNV+Indel Recurrence Analysis"))
recurrence_counts <- recurrence_analysis2("CDH", "TES", cdh_muts_TES[includes_at_least_one(cdh_muts_TES$TS_gene, HDE_genes),], ssc_muts_TES[includes_at_least_one(ssc_muts_TES$TS_gene, HDE_genes),], case_samples=unique(cdh$sample), control_samples=unique(ssc$sample), output_path("recurrence_CDH_muts_HDE_TES.pdf"), main=paste0("CDH HDE 3'UTR Region SNV+Indel Recurrence Analysis"), num_iterations=25000)
recurrence_counts <- recurrence_analysis2("CDH", "TSS", cdh_muts_TSS[includes_at_least_one(cdh_muts_TSS$TS_gene, HDE_genes),], ssc_muts_TSS[includes_at_least_one(ssc_muts_TSS$TS_gene, HDE_genes),], case_samples=unique(cdh$sample), control_samples=unique(ssc$sample), output_path("recurrence_CDH_muts_HDE_TSS.pdf"), main=paste0("CDH HDE TSS Region SNV+Indel Recurrence Analysis"), num_iterations=25000)

# Run these!
recurrence_counts <- recurrence_analysis2("CDH", "TES", cdh_muts_TES, ssc_518_muts_TES, case_samples=unique(cdh$sample), control_samples=unique(ssc_518$sample), output_path("recurrence_CDH1_muts_TES.pdf"), main=paste0("CDH1 3'UTR Region SNV+Indel Recurrence Analysis"), plot_sim=TRUE)
recurrence_counts <- recurrence_analysis2("CDH", "TSS", cdh_muts_TSS, ssc_518_muts_TSS, case_samples=unique(cdh$sample), control_samples=unique(ssc_518$sample), output_path("recurrence_CDH1_muts_TSS.pdf"), main=paste0("CDH1 TSS Region SNV+Indel Recurrence Analysis"), plot_sim=TRUE)
recurrence_counts <- recurrence_analysis2("CDH", "TES", cdh2_muts_TES, ssc_1088_muts_TES, case_samples=unique(cdh2$sample), control_samples=unique(ssc_1088$sample), output_path("recurrence_CDH2_muts_TES.pdf"), main=paste0("CDH2 3'UTR Region SNV+Indel Recurrence Analysis"), plot_sim=TRUE)
recurrence_counts <- recurrence_analysis2("CDH", "TSS", cdh2_muts_TSS, ssc_1088_muts_TSS, case_samples=unique(cdh2$sample), control_samples=unique(ssc_1088$sample), output_path("recurrence_CDH2_muts_TSS.pdf"), main=paste0("CDH2 TSS Region SNV+Indel Recurrence Analysis"), plot_sim=TRUE)

recurrence_counts <- recurrence_analysis2("CDH", "TES", cdh_combined_muts_TES, ssc_all_muts_TES, case_samples=unique(cdh_combined$sample), control_samples=unique(ssc_all$sample), output_path("recurrence_CDH_combined_muts_TES.pdf"), main=paste0("CDH 3'UTR Region SNV+Indel Recurrence Analysis"), plot_sim=FALSE, num_iterations=25000)
recurrence_counts <- recurrence_analysis2("CDH", "TSS", cdh_combined_muts_TSS, ssc_all_muts_TSS, case_samples=unique(cdh_combined$sample), control_samples=unique(ssc_all$sample), output_path("recurrence_CDH_combined_muts_TSS.pdf"), main=paste0("CDH TSS Region SNV+Indel Recurrence Analysis"), plot_sim=FALSE, num_iterations=25000)
recurrence_counts <- recurrence_analysis2("CDH", "TES", cdh_combined_muts_TES[includes_at_least_one(cdh_combined_muts_TES$TS_gene, HDE_genes),], ssc_all_muts_TES[includes_at_least_one(ssc_all_muts_TES$TS_gene, HDE_genes),], case_samples=unique(cdh_combined$sample), control_samples=unique(ssc_all$sample), output_path("recurrence_CDH_combined_muts_HDE_TES.pdf"), main=paste0("CDH HDE 3'UTR Region SNV+Indel Recurrence Analysis"), num_iterations=25000)
recurrence_counts <- recurrence_analysis2("CDH", "TSS", cdh_combined_muts_TSS[includes_at_least_one(cdh_combined_muts_TSS$TS_gene, HDE_genes),], ssc_all_muts_TSS[includes_at_least_one(ssc_all_muts_TSS$TS_gene, HDE_genes),], case_samples=unique(cdh_combined$sample), control_samples=unique(ssc_all$sample), output_path("recurrence_CDH_combined_muts_HDE_TSS.pdf"), main=paste0("CDH HDE TSS Region SNV+Indel Recurrence Analysis"), num_iterations=25000)
# combined 3'UTR and TSS
recurrence_counts <- recurrence_analysis2("CDH", "TSS", rbind(cdh_combined_muts_TSS,cdh_combined_muts_TES), rbind(ssc_all_muts_TSS,ssc_all_muts_TES), case_samples=unique(cdh_combined$sample), control_samples=unique(ssc_all$sample), output_path("recurrence_CDH_combined_muts_TES_plus_TSS.pdf"), main=paste0("CDH 3'UTR+TSS Region SNV+Indel Recurrence Analysis"), num_iterations=25000)
# CHDFB
recurrence_counts <- recurrence_analysis2("CHD", "TES", chdfb_muts_TES, sscfb_muts_TES, case_samples=unique(chdfb$sample), control_samples=unique(sscfb$sample), output_path("recurrence_CHDFB_muts_TES.pdf"), main=paste0("CHD 3'UTR Region SNV+Indel Recurrence Analysis"), plot_sim=FALSE, num_iterations=25000)
recurrence_counts <- recurrence_analysis2("CHD", "TSS", chdfb_muts_TSS, sscfb_muts_TSS, case_samples=unique(chdfb$sample), control_samples=unique(sscfb$sample), output_path("recurrence_CHDFB_muts_TSS.pdf"), main=paste0("CHD TSS Region SNV+Indel Recurrence Analysis"), plot_sim=FALSE, num_iterations=25000)
recurrence_counts <- recurrence_analysis2("CHD2", "TES", chdfb2_muts_TES, sscfb2_muts_TES, case_samples=unique(chdfb2$sample), control_samples=unique(sscfb2$sample), output_path("recurrence_CHDFB2_muts_TES.pdf"), main=paste0("CHD2 3'UTR Region SNV+Indel Recurrence Analysis"), plot_sim=FALSE, num_iterations=25000)
recurrence_counts <- recurrence_analysis2("CHD2", "TSS", chdfb2_muts_TSS, sscfb2_muts_TSS, case_samples=unique(chdfb2$sample), control_samples=unique(sscfb2$sample), output_path("recurrence_CHDFB2_muts_TSS.pdf"), main=paste0("CHD2 TSS Region SNV+Indel Recurrence Analysis"), plot_sim=FALSE, num_iterations=25000)
recurrence_counts <- recurrence_analysis2("CHD", "TES", chdfb_combined_muts_TES, sscfb_combined_muts_TES, case_samples=unique(chdfb_combined$sample), control_samples=unique(sscfb_combined$sample), output_path("recurrence_CHDFB_combined_muts_TES.pdf"), fix_variant_count=TRUE, main=paste0("CHD 3'UTR Region SNV+Indel Recurrence Analysis"), plot_sim=FALSE, num_iterations=25000)
recurrence_counts <- recurrence_analysis2("CHD", "TSS", chdfb_combined_muts_TSS, sscfb_combined_muts_TSS, case_samples=unique(chdfb_combined$sample), control_samples=unique(sscfb_combined$sample), output_path("recurrence_CHDFB_combined_muts_TSS.pdf"), fix_variant_count=TRUE, main=paste0("CHD TSS Region SNV+Indel Recurrence Analysis"), plot_sim=FALSE, num_iterations=25000)


dat <- read.table(output_path("CHD_3UTR_variants.tsv"), header=TRUE)
dat <- read.table(output_path("CDH_3UTR_variants.tsv"), header=TRUE)
query <- c("SHOC2","NF1","KALRN")
query <- c("GATA4")
dat[(! is.na(dat$RBP)),][includes_at_least_one(dat$gene[! is.na(dat$RBP)], query),]
dat[(is.na(dat$RBP)),][includes_at_least_one(dat$gene[is.na(dat$RBP)], query),]


expression_ranks <- get_genes_by_expression(return_type="heart_rank") # get expression rank scores for all possible genes
heart_expression_ranks_dat <- t(sapply(ls(expression_ranks), function(x) return(c(x, expression_ranks[[x]]))))
colnames(heart_expression_ranks_dat) <- c("gene", "heart_rank")
write.csv(heart_expression_ranks_dat, file=output_path("heart_rank.csv"), row.names=FALSE)
expression_ranks <- get_genes_by_expression(return_type="diaphragm_rank") # get expression rank scores for all possible genes
diaphragm_expression_ranks_dat <- t(sapply(ls(expression_ranks), function(x) return(c(x, expression_ranks[[x]]))))
colnames(diaphragm_expression_ranks_dat) <- c("gene", "diaphragm_rank")
write.csv(diaphragm_expression_ranks_dat, file=output_path("diaphragm_rank.csv"), row.names=FALSE)
pLI_scores <- get_constrained_genes("pLI>=0")
pLI_scores_dat <- t(sapply(ls(pLI_scores), function(x) return(c(x, pLI_scores[[x]]))))
colnames(pLI_scores_dat) <- c("gene", "pLI")
write.csv(pLI_scores_dat, file=output_path("pli_scores.csv"), row.names=FALSE)
mis_z_scores <- get_constrained_genes("mis_z>=-Inf")
mis_z_scores_dat <- t(sapply(ls(mis_z_scores), function(x) return(c(x, mis_z_scores[[x]]))))
colnames(mis_z_scores_dat) <- c("gene", "mis_z")
write.csv(mis_z_scores_dat, file=output_path("mis_z_scores.csv"), row.names=FALSE)
# Write tables of recurrently mutated genes (with no mutations in controls, if remove_genes_hit_in_controls=TRUE).
write_recurrent_genes_table <- function(disease, sites, cases_snps, cases_indels, controls, filename, cases2_snps=NULL, cases2_indels=NULL, controls2=NULL, include_expected=TRUE, remove_genes_hit_in_controls=FALSE, TAD=TRUE, return_table=FALSE) {
    if(grepl("CDH", disease)) { expression_query <- "diaphragm_rank" } else if(grepl("CHD", disease)) { expression_query <- "heart_rank" } else { print("ERROR: Invalid disease!"); return() }
    expression_ranks <- get_genes_by_expression(return_type=expression_query) # get expression rank scores for all possible genes
    pLI_scores <- get_constrained_genes("pLI>=0") # get pLI scores for all possible genes
    sample_count = length(unique(c(paste0(cases_snps$sample), paste0(cases_indels$sample), paste0(cases2_snps$sample), paste0(cases2_indels$sample))))
    if(remove_genes_hit_in_controls) { filename = gsub("\\.", "_filtered.", filename) }
    
    case_snps_genes <- unlist(strsplit(paste0(cases_snps$TS_gene), ","))
    case_indels_genes <- unlist(strsplit(paste0(cases_indels$TS_gene), ","))
    control_genes <- unlist(strsplit(paste0(controls$TS_gene), ","))
    if(!is.null(controls2)) { control_genes <- c(control_genes, unlist(strsplit(paste0(controls2$TS_gene), ","))) }
    if(remove_genes_hit_in_controls) { case_snps_genes <- case_snps_genes[!(case_snps_genes %in% control_genes)]; case_indels_genes <- case_indels_genes[!(case_indels_genes %in% control_genes)] }
    snps_genes <- data.frame(table(case_snps_genes))
    indels_genes <- data.frame(table(case_indels_genes))
    #genes <- merge(snps_genes, indels_genes, by="Var1", all=TRUE)
    genes <- merge(snps_genes, indels_genes, by.x=colnames(snps_genes)[1], by.y=colnames(indels_genes)[1], all=TRUE)
    colnames(genes) <- c("gene", paste0(disease,"_snps"), paste0(disease,"_indels"))
    desired_columns <- c("X.chr", "pos", "ref", "alt", "sample", "TS_gene", "RBP")
    all_variants <- rbind(cases_snps[,which(colnames(cases_snps) %in% desired_columns)], cases_indels[,which(colnames(cases_indels) %in% desired_columns)])
    if(!is.null(cases2_snps) && !is.null(cases2_indels)) {
        case2_snps_genes <- unlist(strsplit(paste0(cases2_snps$TS_gene), ","))
        case2_indels_genes <- unlist(strsplit(paste0(cases2_indels$TS_gene), ","))
        if(remove_genes_hit_in_controls) { case2_snps_genes <- case2_snps_genes[!(case2_snps_genes %in% control_genes)]; case2_indels_genes <- case2_indels_genes[!(case2_indels_genes %in% control_genes)] }
        snps2_genes <- data.frame(table(case2_snps_genes))
        indels2_genes <- data.frame(table(case2_indels_genes))
        #genes2 <- merge(snps2_genes, indels2_genes, by="Var1", all=TRUE)
        genes2 <- merge(snps2_genes, indels2_genes, by.x=colnames(snps2_genes)[1], by.y=colnames(indels2_genes)[1], all=TRUE)
        colnames(genes2) <- c("gene", paste0(disease,"2_snps"), paste0(disease,"2_indels"))
        genes2[is.na(genes2)] <- 0
        
        genes <- merge(genes, genes2, by="gene", all=TRUE)
        all_variants <- rbind(all_variants, rbind(cases2_snps[,which(colnames(cases2_snps) %in% desired_columns)], cases2_indels[,which(colnames(cases2_indels) %in% desired_columns)]))
    }
    genes[is.na(genes)] <- 0
    genes <- cbind(genes, rowSums(genes[,c(-1)]))
    colnames(genes)[ncol(genes)] <- "total_hits"
    all_variants <- cbind(all_variants, paste0(all_variants$sample, "_chr", all_variants$X.chr, ":", all_variants$pos, all_variants$ref, ">", all_variants$alt))
    colnames(all_variants)[ncol(all_variants)] <- "variant"
    
    recurrent_genes <- genes[genes$total_hits > 1,]
    pLI <- unlist(sapply(recurrent_genes$gene, function(gene) { annot <- pLI_scores[[paste(gene)]]; if(is.null(annot)) {return(NA)} else {return(paste0(round(annot,3)))} } ))
    expression_rank <- unlist(sapply(recurrent_genes$gene, function(gene) { annot <- expression_ranks[[paste(gene)]]; if(is.null(annot)) {return(NA)} else {return(paste0(round(annot,3)))} } ))
    recurrent_genes <- cbind(recurrent_genes, rep(0,nrow(recurrent_genes)), pLI, expression_rank, rep(NA,nrow(recurrent_genes)), rep(NA,nrow(recurrent_genes)))
    colnames(recurrent_genes)[c((ncol(recurrent_genes)-4),(ncol(recurrent_genes)-2):ncol(recurrent_genes))] <- c("total_samples", expression_query, "variants", "RBPs_disrupted")
    recurrent_genes$pLI <- paste0(recurrent_genes$pLI); recurrent_genes[,expression_query] <- paste0(recurrent_genes[,expression_query]) # Change factor columns to String columns, to be able to edit them in loop below.
    recurrent_genes$variants <- paste0(recurrent_genes$variants); recurrent_genes$RBPs_disrupted <- paste0(recurrent_genes$RBPs_disrupted) # Change factor columns to String columns, to be able to edit them in loop below.
    recurrent_genes[recurrent_genes == "NA"] <- NA
    for(i in 1:nrow(recurrent_genes)) {
        print(paste0("Doing variant-level annotation... [", i, " / ", nrow(recurrent_genes), "]"))
        variants <- all_variants[includes_at_least_one(all_variants$TS_gene, recurrent_genes$gene[i]),]
        recurrent_genes$total_samples[i] <- length(unique(variants$sample))
        recurrent_genes$variants[i] <- paste(sort(unique(variants$variant)),collapse=",")
        rbps <- unique(variants$RBP[!is.na(variants$RBP)])
        if(length(rbps)>0) { recurrent_genes$RBPs_disrupted[i] <- paste(sort(rbps),collapse=",") }
    }
    print("Done.")
    
    recurrent_genes <- recurrent_genes[recurrent_genes$total_samples > 1,]
    
    # Annotate clusters
    if(TAD) { TAD_string <- "TAD" } else { TAD_string <- NULL }
    cluster_granges_name = paste(c("clustered", sites, TAD_string, "granges"), collapse="_")
    cluster_granges <- try(get(cluster_granges_name))
    if(inherits(cluster_granges, "try-error")) { 
        region_granges <- get(paste0(c(sites, "_granges")))
        cluster_granges <- cluster_overlapping_regions(region_granges, region_name=sites, break_by_TAD=TAD, make_plots=FALSE)
    }
    print("Hashing gene-to-cluster mappings...")
    gene_cluster_mapping <- new.env()
    for(i in 1:length(cluster_granges)) {
        cluster <- paste0(names(cluster_granges)[i])
        all_genes_in_cluster <- unlist(strsplit(cluster, ","))
        for(cluster_gene in all_genes_in_cluster) { gene_cluster_mapping[[cluster_gene]] <- cluster }
    }
    print("Done.")
    cluster <- unlist(sapply(recurrent_genes$gene, function(gene) { annot <- gene_cluster_mapping[[paste(gene)]]; if(is.null(annot)) {return(NA)} else {return(paste0(annot))} } ))
    cluster_region <- paste0(cluster_granges[cluster])
    recurrent_genes <- cbind(recurrent_genes, cluster, cluster_region)
    
    if(include_expected) {  
        #if(bin_size > 1) { bin_file = full_path(work_folder, paste0("chr",chromosome,"_",start(g_chr[bin_start]),"-",end(g_chr[bin_end]),"_bin",bin,"_",bin_size,"_granges_trimer_counts.csv"))
        #} else { bin_file = full_path(work_folder, paste0("chr",chromosome,"_",start(g_chr[bin_start]),"-",end(g_chr[bin_end]),"_trimer_counts.csv")) }
        #bin_files <- c(bin_files, bin_file)
        #if (!restart && file.exists(bin_file)) {}
        #tri_counts_files <- list.files(output_path(paste0("trimer_counts_for_",sites,"_temp2")))
        trimer_counts_folder = output_path(paste0("trimer_counts_for_",sites,"_temp2"))
        region_granges <- get(paste0(sites,"_granges"))
        background_expectation_per_individual <- sapply(1:nrow(recurrent_genes), function(i) {
            region_grange <- region_granges[recurrent_genes$gene[i]]
            chromosome = paste0(seqnames(region_grange)); start = start(region_grange); end = end(region_grange)
            tri_counts <- read.csv(file=full_path(trimer_counts_folder, paste0("chr",chromosome,"_",start,"-",end,"_trimer_counts.csv")))
            tri_counts_names <- names(tri_counts); tri_counts <- as.numeric(tri_counts); names(tri_counts) <- tri_counts_names
            tri_mut_rates <- get(paste0(tolower(sites),"_trimer_mutation_rates_from_ssc_cpg_meth_aggregate"))
            expected_muts = 0
            for(tri in names(tri_counts[tri_counts > 0])) {
                expected_muts = expected_muts + as.numeric(tri_counts[tri]) * sum(tri_mut_rates[grepl(paste0(tri,"->"),names(tri_mut_rates))])
            }
            return(expected_muts)
        })
        background_expectation <- background_expectation_per_individual * sample_count
        expectation_dat <- cbind(ppois(recurrent_genes$total_hits - 1, lambda=background_expectation, lower.tail=FALSE), background_expectation, background_expectation_per_individual); colnames(expectation_dat)[1] <- "ppois"
        recurrent_genes <- cbind(recurrent_genes[,1:which(colnames(recurrent_genes) == expression_query)], expectation_dat, recurrent_genes[,(which(colnames(recurrent_genes) == expression_query)+1):ncol(recurrent_genes)])
        recurrent_genes <- recurrent_genes[order(recurrent_genes$ppois, decreasing=FALSE),]
    } else {
        recurrent_genes <- recurrent_genes[order(recurrent_genes$total_samples, recurrent_genes$total_hits, (recurrent_genes[,2] + recurrent_genes[,3]), recurrent_genes$pLI, recurrent_genes[,expression_query], recurrent_genes$RBPs_disrupted, decreasing=TRUE),]
    }
    
    recurrent_genes[is.na(recurrent_genes)] <- "."
    write.table(recurrent_genes, file=filename, sep="\t", row.names=FALSE)
    if (return_table) { return(recurrent_genes) } else { return() }
}
write_recurrent_genes_table("CHD", "TES", chd_snps_TES, chd_indels_TES, controls=sscfb_muts_TES, output_path("CHD_3UTR_recurrent_genes2.tsv"))
write_recurrent_genes_table("CHD", "TSS", chd_snps_TSS, chd_indels_TSS, controls=sscfb_muts_TSS, output_path("CHD_TSS_recurrent_genes.tsv"))
write_recurrent_genes_table("CDH", "TES", cdh_snps_TES, cdh_indels_TES, controls=ssc_muts_TES, output_path("CDH1_3UTR_recurrent_genes2.tsv"), return_table=TRUE)
write_recurrent_genes_table("CDH", "TSS", cdh_snps_TSS, cdh_indels_TSS, controls=ssc_muts_TSS, output_path("CDH1_TSS_recurrent_genes.tsv"))
write_recurrent_genes_table("CDH", "TES", cdh_new_snps_TES, cdh_new_indels_TES, controls=ssc_new_muts_TES, output_path("CDH2_3UTR_recurrent_genes.tsv"))
write_recurrent_genes_table("CDH", "TSS", cdh_new_snps_TSS, cdh_new_indels_TSS, controls=ssc_new_muts_TSS, output_path("CDH2_TSS_recurrent_genes.tsv"))
write_recurrent_genes_table("CDH", "TES", cdh_snps_TES, cdh_indels_TES, controls=ssc_muts_TES, cases2_snps=cdh_new_snps_TES, cases2_indels=cdh_new_indels_TES, controls2=ssc_new_muts_TES, output_path("CDH_combined_3UTR_recurrent_genes.tsv"))
write_recurrent_genes_table("CDH", "TSS", cdh_snps_TSS, cdh_indels_TSS, controls=ssc_muts_TSS, cases2_snps=cdh_new_snps_TSS, cases2_indels=cdh_new_indels_TSS, controls2=ssc_new_muts_TSS, output_path("CDH_combined_TSS_recurrent_genes.tsv"))

cdh_recurrent_tes <- write_recurrent_genes_table("CDH", "TES", cdh_snps_TES, cdh_indels_TES, controls=ssc_518_muts_TES, cases2_snps=cdh2_snps_TES, cases2_indels=cdh2_indels_TES, controls2=ssc_1088_muts_TES, output_path("CDH_489_3UTR_recurrent_genes.tsv"), return_table=TRUE)
cdh_recurrent_tss <- write_recurrent_genes_table("CDH", "TSS", cdh_snps_TSS, cdh_indels_TSS, controls=ssc_518_muts_TSS, cases2_snps=cdh2_snps_TSS, cases2_indels=cdh2_indels_TSS, controls2=ssc_1088_muts_TSS, output_path("CDH_489_TSS_recurrent_genes.tsv"), return_table=TRUE)
chd_recurrent_tes <- write_recurrent_genes_table("CHD", "TES", chd_snps_TES, chd_indels_TES, controls=ssc_518_muts_TES, cases2_snps=chd2_snps_TES, cases2_indels=chd2_indels_TES, controls2=ssc_1088_muts_TES, output_path("CHD_745_3UTR_recurrent_genes.tsv"), return_table=TRUE)
chd_recurrent_tss <- write_recurrent_genes_table("CHD", "TSS", chd_snps_TSS, chd_indels_TSS, controls=ssc_518_muts_TSS, cases2_snps=chd2_snps_TSS, cases2_indels=chd2_indels_TSS, controls2=ssc_1088_muts_TSS, output_path("CHD_745_TSS_recurrent_genes.tsv"), return_table=TRUE)
chdfb_recurrent_tes <- write_recurrent_genes_table("CHD", "TES", chdfb_snps_TES, chdfb_indels_TES, controls=sscfb_muts_TES, cases2_snps=chdfb2_snps_TES, cases2_indels=chdfb2_indels_TES, controls2=sscfb2_muts_TES, output_path("CHDFB_combined_3UTR_recurrent_genes.tsv"), return_table=TRUE, include_expected=FALSE)
chdfb_recurrent_tss <- write_recurrent_genes_table("CHD", "TSS", chdfb_snps_TSS, chdfb_indels_TSS, controls=sscfb_muts_TSS, cases2_snps=chdfb2_snps_TSS, cases2_indels=chdfb2_indels_TSS, controls2=sscfb2_muts_TSS, output_path("CHDFB_combined_TSS_recurrent_genes.tsv"), return_table=TRUE, include_expected=FALSE)


##############################################################################

TSS_burdens_env <- analyze_burden(chd, ssc)
TSS_burdens <- TSS_burdens_env[["qualitative"]]
epigenome_names <- get_roadmap_epigenome_names(unlist(lapply(seq(1:length(TSS_burdens$feature)), function(i) { return(strsplit(paste(TSS_burdens$feature[i]), "\\.")[[1]][1]) } )))
TSS_burdens <- cbind(TSS_burdens, epigenome_names)
write.csv(TSS_burdens, file="Felix_check_burden_analysis.csv", row.names=FALSE)
write.csv(TSS_burdens[TSS_burdens$enrichment > 1.5,], file="Felix_check_burden_analysis_strong.csv", row.names=FALSE)
felix_tissues <- c("E083", "E104", "E095", "E105", "E013") # Fetal Heart, Right Atrium, Left Ventricle, Right Ventricle, and Mesoderm, respectively.
felix_histone_marks <- c("H3K27ac.gappedPeak", "H3K4me1.gappedPeak")
felix_features_grid <- expand.grid(felix_tissues, felix_histone_marks)
felix_features <- paste0(felix_features_grid$Var1, ".", felix_features_grid$Var2)
felix_enhancer_burdens <- TSS_burdens[TSS_burdens$feature %in% felix_features,]
write.csv(felix_enhancer_burdens, file="Felix_check_burden_analysis_specific.csv", row.names=FALSE)
paste0(nrow(chd), "/", nrow(chd_wgsa))
paste0(nrow(ssc), "/", nrow(ssc_wgsa))
cases_Y <- nrow(chd)
cases_N <- nrow(chd_wgsa) - cases_Y
controls_Y <- nrow(ssc)
controls_N <- nrow(ssc_wgsa) - controls_Y
ft <- fisher.test(matrix(c(cases_Y, controls_Y, cases_N, controls_N), nrow=2, ncol=2), alternative="greater")

skip <- function() {
    
    CDH_HE_TSS_20k_snps <- cases[is_TS(cases$Chrom, cases$Position, sites="TSS"),]
    CDH_HE_TES_20k_snps <- cases[is_TS(cases$Chrom, cases$Position, sites="TES"),]
    SSC_HE_TSS_20k_snps <- controls[is_TS(controls$Chrom, controls$Position, sites="TSS"),]
    SSC_HE_TES_20k_snps <- controls[is_TS(controls$Chrom, controls$Position, sites="TES"),]
    
    #CDH_HE_TSS_20k_snps <- cases[unlist(lapply(seq(1:nrow(cases)), function(i) { print(i); return(sum(cases$X.chr[i] == paste0(HE_TSS_20k_regions$chromosome) & cases$pos[i] >= HE_TSS_20k_regions$start & cases$pos[i] <= HE_TSS_20k_regions$end) > 0) })),]
    #CDH_HE_TES_20k_snps <- cases[unlist(lapply(seq(1:nrow(cases)), function(i) { print(i); return(sum(cases$X.chr[i] == paste0(HE_TES_20k_regions$chromosome) & cases$pos[i] >= HE_TES_20k_regions$start & cases$pos[i] <= HE_TES_20k_regions$end) > 0) })),]
    #SSC_HE_TSS_20k_snps <- cases[unlist(lapply(seq(1:nrow(controls)), function(i) { print(i); return(sum(controls$X.chr[i] == paste0(HE_TSS_20k_regions$chromosome) & controls$pos[i] >= HE_TSS_20k_regions$start & controls$pos[i] <= HE_TSS_20k_regions$end) > 0) })),]
    #SSC_HE_TES_20k_snps <- cases[unlist(lapply(seq(1:nrow(controls)), function(i) { print(i); return(sum(controls$X.chr[i] == paste0(HE_TES_20k_regions$chromosome) & controls$pos[i] >= HE_TES_20k_regions$start & controls$pos[i] <= HE_TES_20k_regions$end) > 0) })),]
    
    bivalent_genes <- read.csv(data_path("E001_bivalent_wgenes.bed"), sep="\t", header=FALSE)
    colnames(bivalent_genes) <- c("chromosome", "start", "end", "something", "ensembl_id")
    #bivalent_genes_granges <- to_genomic_regions(bivalent_genes, labels=genomic_coordinates_to_strings(bivalent_genes$chromosome, bivalent_genes$start, bivalent_genes$end))
    bivalent_genes_names <- c()
    for(ensembl_id in bivalent_genes$ensembl_id) {
        name <- ensembl_genes[[paste(ensembl_id)]]
        if(!is.null(name)) { 
            bivalent_genes_names <- c(bivalent_genes_names, name) 
        } else { bivalent_genes_names <- c(bivalent_genes_names, NA) } 
    } 
    bivalent_genes <- cbind(bivalent_genes, bivalent_genes_names)
    colnames(bivalent_genes)[ncol(bivalent_genes)] <- "gene"
    #heart_brain_rank <- read.csv("heart_brain_rank.csv")
    #HE <- heart_brain_rank[!duplicated(HE$human.External.Gene.Name),]
    #HE <- HE[HE$e14.5_rank <= 25,]
    #HE_genes <- unique(HE$human.External.Gene.Name) 
    biv_genes <- unique(bivalent_genes$gene)
    
    fraction_bivalent_HE <- sum(HE_genes %in% biv_genes)/length(HE_genes)
    fraction_HE_bivalent <- sum(biv_genes %in% HE_genes)/length(biv_genes)
    
    #fraction_bivalent_HE <- sum(HE_genes %in% bivalent_genes$gene)/length(HE_genes)
    #fraction_HE_bivalent <- sum(bivalent_genes$gene %in% HE_genes)/length(unique(bivalent_genes$gene))
    #fraction_bivalent_human <- sum(human_genes %in% bivalent_genes$gene)/length(human_genes)
    
    cases_H3K27me3 <- cases_CDH_snps[,grepl("H3K27me3", colnames(cases_CDH_snps))]
    cases_H3K27me3_indices <- sapply(1:nrow(cases_H3K27me3), function(i) { print(i); row <- cases_H3K27me3[i,]; return("Y" %in% paste(row)) } ) #return(sum(row=="Y")>0)
    cases_H3K27me3 <- cases_CDH_snps[cases_H3K27me3_indices,]
    
    H3K27me3_indices <- c()
    H3K27me3_cols <- which(grepl("H3K27me3", colnames(cases_CDH_snps)))
    for(i in 1:nrow(cases_CDH_snps)) {
        print(cases_CDH_snps[i,])
        print(paste(cases_CDH_snps[i,c(H3K27me3_cols)]))
        readline()
    }
    
    # Calculate raw nc burden in TSS and TES!
    cases_HE_TSS <- cases[is_TS(cases$Chrom, cases$Position, sites="TSS"),]
    cases_HE_TES <- cases[is_TS(cases$Chrom, cases$Position, sites="TES"),]
    controls_HE_TSS <- controls[is_TS(controls$Chrom, controls$Position, sites="TSS"),]
    controls_HE_TES <- controls[is_TS(controls$Chrom, controls$Position, sites="TES"),]
    
    cases_Y <- nrow(cases_HE_TSS)
    cases_N <- nrow(cases) - cases_Y
    controls_Y <- nrow(controls_HE_TSS)
    controls_N <- nrow(controls) - controls_Y
    ft_nc_TSS <- fisher.test(matrix(c(cases_Y, controls_Y, cases_N, controls_N), nrow=2, ncol=2), alternative="greater")
    matrix(c(cases_Y, controls_Y, cases_N, controls_N), nrow=2, ncol=2)
    ft_nc_TSS
    
    cases_Y <- nrow(cases_HE_TES)
    cases_N <- nrow(cases) - cases_Y
    controls_Y <- nrow(controls_HE_TES)
    controls_N <- nrow(controls) - controls_Y
    ft_nc_TES <- fisher.test(matrix(c(cases_Y, controls_Y, cases_N, controls_N), nrow=2, ncol=2), alternative="greater")
    matrix(c(cases_Y, controls_Y, cases_N, controls_N), nrow=2, ncol=2)
    ft_nc_TES
    
    
    run_tests <- function(dat, SAMPLE_COUNT) {
        controls_350 <- dat[dat$sample %in% sample(unique(paste(dat$sample)), SAMPLE_COUNT+1),]
        over_1_snv_SSC_jackknife_TSS <- jackknife_snv_count(controls_350, sites="TSS", block_size=1)
        over_1_snv_SSC_jackknife_TES <- jackknife_snv_count(controls_350, sites="TES", block_size=1)
        #over_2_snv_SSC_jackknife_TSS <- jackknife_snv_count(controls_350, sites="TSS", block_size=1, n_greater_than=2)
        #over_2_snv_SSC_jackknife_TES <- jackknife_snv_count(controls_350, sites="TES", block_size=1, n_greater_than=2)
        #over_1_mut_SSC_jackknife_TSS <- jackknife_snv_count(ssc_controls_350, sites="TSS", block_size=1)
        #over_1_mut_SSC_jackknife_TES <- jackknife_snv_count(ssc_controls_350, sites="TES", block_size=1)
        #over_2_mut_SSC_jackknife_TSS <- jackknife_snv_count(ssc_controls_350, sites="TSS", block_size=1, n_greater_than=2)
        #over_2_mut_SSC_jackknife_TES <- jackknife_snv_count(ssc_controls_350, sites="TES", block_size=1, n_greater_than=2)
        
        results <- new.env()
        results[["over_1_snv_SSC_jackknife_TSS"]] <- over_1_snv_SSC_jackknife_TSS
        results[["over_1_snv_SSC_jackknife_TES"]] <- over_1_snv_SSC_jackknife_TES
        return(results)
    }
    hits <- unique(data.frame(olRanges(controls_granges, bivalent_genes_granges)))
    controls_biv <- controls[genomic_coordinates_to_strings(controls$Chrom, controls$Position) %in% genomic_coordinates_to_strings(hits$seqnames, hits$start),]
    hits <- unique(data.frame(olRanges(cases_granges, bivalent_genes_granges)))
    cases_biv <- cases[gsub("chr", "", genomic_coordinates_to_strings(cases$Chrom, cases$Position)) %in% genomic_coordinates_to_strings(hits$seqnames, hits$start),]
    #biv_results_cases <- run_tests(cases_biv, SAMPLE_COUNT=length(unique(controls_biv$sample)))
    biv_results <- run_tests(controls_biv, SAMPLE_COUNT=length(unique(cases_biv$sample)))
    
    controls_granges_TSS <- controls_granges[is_TS(controls$Chrom, controls$Position, sites="TSS")]
    controls_granges_TES <- controls_granges[is_TS(controls$Chrom, controls$Position, sites="TES")]
    cases_granges_TSS <- controls_granges[is_TS(cases$Chrom, cases$Position, sites="TSS")]
    cases_granges_TES <- controls_granges[is_TS(cases$Chrom, cases$Position, sites="TES")]
    hits <- unique(data.frame(olRanges(controls_granges_TSS, bivalent_genes_granges)))
    controls_biv_TSS <- controls[genomic_coordinates_to_strings(controls$Chrom, controls$Position) %in% genomic_coordinates_to_strings(hits$seqnames, hits$start),]
    hits <- unique(data.frame(olRanges(cases_granges_TSS, bivalent_genes_granges)))
    cases_biv_TSS <- cases[gsub("chr", "", genomic_coordinates_to_strings(cases$Chrom, cases$Position)) %in% genomic_coordinates_to_strings(hits$seqnames, hits$start),]
    hits <- unique(data.frame(olRanges(controls_granges_TES, bivalent_genes_granges)))
    controls_biv_TES <- controls[genomic_coordinates_to_strings(controls$Chrom, controls$Position) %in% genomic_coordinates_to_strings(hits$seqnames, hits$start),]
    hits <- unique(data.frame(olRanges(cases_granges_TES, bivalent_genes_granges)))
    cases_biv_TES <- cases[gsub("chr", "", genomic_coordinates_to_strings(cases$Chrom, cases$Position)) %in% genomic_coordinates_to_strings(hits$seqnames, hits$start),]
    
    GO_dat <- read.csv(file="GO_stemcell_annotation.txt", sep="\t", header=FALSE)
    GO_genes <- unique(GO_dat[,3])
    GO_genes_ranges <- hg19_genebody[hg19_genebody[,5] %in% GO_genes, c(1,2,3)]
    colnames(GO_genes_ranges) <- c("chromosome", "start", "end")
    GO_genes_granges <- to_genomic_regions(GO_genes_ranges, labels=genomic_coordinates_to_strings(GO_genes_ranges$chromosome, GO_genes_ranges$start, GO_genes_ranges$end))
    
    hits <- unique(data.frame(olRanges(controls_granges, GO_genes_granges)))
    controls_GO <- controls[genomic_coordinates_to_strings(controls$Chrom, controls$Position) %in% genomic_coordinates_to_strings(hits$seqnames, hits$start),]
    hits <- unique(data.frame(olRanges(cases_granges, GO_genes_granges)))
    cases_GO <- cases[gsub("chr", "", genomic_coordinates_to_strings(cases$Chrom, cases$Position)) %in% genomic_coordinates_to_strings(hits$seqnames, hits$start),]
    GO_results <- run_tests(controls_GO, SAMPLE_COUNT=length(unique(cases_GO$sample)))
    GO_TSS_genes_dat <- build_mutations_dat(cases_GO$Chrom, cases_GO$Position, is_indel(cases_GO$Ref,cases_GO$Alt), controls_GO$Chrom, controls_GO$Position, is_indel(controls_GO$Ref,controls_GO$Alt), sites="TSS")
    GO_TES_genes_dat <- build_mutations_dat(cases_GO$Chrom, cases_GO$Position, is_indel(cases_GO$Ref,cases_GO$Alt), controls_GO$Chrom, controls_GO$Position, is_indel(controls_GO$Ref,controls_GO$Alt), sites="TES")
    # GO HE TSS, with >1 SNVs
    plot_recurrance(cases_dat=GO_TSS_genes_dat, mut_type="mut", n_greater_than=1, controls_jackknife=GO_results[["over_1_snv_SSC_jackknife_TSS"]])
    # GO HE TES, with >1 SNVs
    plot_recurrance(cases_dat=GO_TES_genes_dat, mut_type="mut", n_greater_than=1, controls_jackknife=GO_results[["over_1_snv_SSC_jackknife_TES"]])
    
    
    test_biv <- function(cases_biv, cases, controls_biv, controls) {
        cases_Y <- nrow(cases_biv)
        cases_N <- nrow(cases) - cases_Y
        controls_Y <- nrow(controls_biv)
        controls_N <- nrow(controls) - controls_Y	
        print(matrix(c(cases_Y, controls_Y, cases_N, controls_N), nrow=2, ncol=2))
        ft_biv <- fisher.test(matrix(c(cases_Y, controls_Y, cases_N, controls_N), nrow=2, ncol=2), alternative="greater")
        return(ft_biv)
    }
    test_biv(cases_biv, cases, controls_biv, controls)
    test_biv(cases_biv_TSS, cases, controls_biv_TSS, controls)
    test_biv(cases_biv_TES, cases, controls_biv_TES, controls)
    
    cases_Y <- nrow(cases_biv)
    cases_N <- nrow(cases) - cases_Y
    controls_Y <- nrow(controls_biv)
    controls_N <- nrow(controls) - controls_Y
    ft_biv <- fisher.test(matrix(c(cases_Y, controls_Y, cases_N, controls_N), nrow=2, ncol=2), alternative="greater")
    
    cases_Y <- nrow(cases_GO)
    cases_N <- nrow(cases) - cases_Y
    controls_Y <- nrow(controls_GO)
    controls_N <- nrow(controls) - controls_Y
    ft_GO <- fisher.test(matrix(c(cases_Y, controls_Y, cases_N, controls_N), nrow=2, ncol=2), alternative="greater")
    
    
    cases_GO_biv <- 0
    for(i in 1:nrow(cases_GO)) {
        if((cases_GO$Chrom[i] %in% cases_biv$Chrom) & (cases_GO$Position[i] %in% cases_biv$Position)) { cases_GO_biv = cases_GO_biv + 1 }
    } 
    cases_Y <- cases_GO_biv
    controls_GO_biv <- 0
    for(i in 1:nrow(controls_GO)) {
        if((controls_GO$Chrom[i] %in% controls_biv$Chrom) & (controls_GO$Position[i] %in% controls_biv$Position)) { controls_GO_biv = controls_GO_biv + 1 }
    }
    controls_Y <- cases_GO_biv
    cases_N <- nrow(cases) - cases_Y
    controls_N <- nrow(controls) - controls_Y
    ft_GO_biv <- fisher.test(matrix(c(cases_Y, controls_Y, cases_N, controls_N), nrow=2, ncol=2), alternative="greater")
    
    cases_Y <- nrow(cases_GO) + nrow(cases_biv) - cases_GO_biv
    controls_Y <- nrow(controls_GO) + nrow(controls_biv) - controls_GO_biv
    cases_N <- nrow(cases) - cases_Y
    controls_N <- nrow(controls) - controls_Y
    ft_GO_or_biv <- fisher.test(matrix(c(cases_Y, controls_Y, cases_N, controls_N), nrow=2, ncol=2), alternative="greater")
    print(matrix(c(cases_Y, controls_Y, cases_N, controls_N), nrow=2, ncol=2))
    print(ft_GO_or_biv)
    
    
    ## All HE TSS, with >1 SNVs
    #plot_recurrance(cases_dat=all_TSS_genes_dat, mut_type="mut", n_greater_than=1, controls_jackknife=biv_results[["over_1_snv_SSC_jackknife_TSS"]])
    ## All HE TES, with >1 SNVs
    #plot_recurrance(cases_dat=all_TES_genes_dat, mut_type="mut", n_greater_than=1, controls_jackknife=biv_results[["over_1_snv_SSC_jackknife_TES"]])
    
    
    
    
    # Annotate mutations in dat with annotation in regions vector.
    annotate_with_feature <- function(dat, annotations, name="annotation", dat_granges=NULL) {
        if(is.null(dat_granges)) { dat_granges <- to_genomic_regions(dat, chr_colname="X.chr", start_colname="pos", end_colname="pos", labels=dat$sample) }
        if(class(annotations)!="GRanges") { annotations <- to_genomic_regions(annotations, labels=genomic_coordinates_to_strings(annotations$chromosome, annotations$start, annotations$end)) }
        hits <- data.frame(olRanges(dat_granges, annotations))
        hits <- unique(hits)
        overlaps_with_annotation <- genomic_coordinates_to_strings(dat$X.chr, dat$pos, dat$pos) %in% genomic_coordinates_to_strings(hits$seqnames, hits$start, hits$end)
        
        annotation_vector <- character(nrow(dat))
        annotation_vector[overlaps_with_annotation] <- "Y"
        annotation_vector[!overlaps_with_annotation] <- "N"
        dat <- cbind(dat, annotation_vector)
        colnames(dat)[ncol(dat)] <- name
        return(dat)
    }
    CDH_HE_TES_20k_snps <- annotate_with_feature(CDH_HE_TES_20k_snps, bivalent_genes_granges, name="bivalent_marks")
    CDH_HE_TSS_20k_snps <- annotate_with_feature(CDH_HE_TSS_20k_snps, bivalent_genes_granges, name="bivalent_marks")
    SSC_HE_TES_20k_snps <- annotate_with_feature(SSC_HE_TES_20k_snps, bivalent_genes_granges, name="bivalent_marks")
    SSC_HE_TSS_20k_snps <- annotate_with_feature(SSC_HE_TSS_20k_snps, bivalent_genes_granges, name="bivalent_marks")
    
    annotate_with_genes <- function(dat, genes, name="annotation") {
        genes_hash <- new.env()
        sapply_out <- sapply(genes, function(gene) { genes_hash[[paste(gene)]] <- 1 } )
        has_bivalent_gene <- sapply(seq(1:nrow(dat)), function(i) { genes <- strsplit(paste0(dat$ANNOVAR_ensembl_Gene_ID[i]), "\\|")[[1]]; for(gene in genes) { if (!is.null(genes_hash[[paste(gene)]])) { return(TRUE) }} ; return(FALSE) } )
        annotation_vector <- character(nrow(dat))
        annotation_vector[has_bivalent_gene] <- "Y"
        annotation_vector[!has_bivalent_gene] <- "N"
        dat <- cbind(dat, annotation_vector)
        colnames(dat)[ncol(dat)] <- name
        return(dat)
    }
    CDH_HE_TES_20k_snps <- annotate_with_genes(CDH_HE_TES_20k_snps, unique(bivalent_genes$ensembl_id), name="bivalent_gene")
    CDH_HE_TSS_20k_snps <- annotate_with_genes(CDH_HE_TSS_20k_snps, unique(bivalent_genes$ensembl_id), name="bivalent_gene")
    SSC_HE_TES_20k_snps <- annotate_with_genes(SSC_HE_TES_20k_snps, unique(bivalent_genes$ensembl_id), name="bivalent_gene")
    SSC_HE_TSS_20k_snps <- annotate_with_genes(SSC_HE_TSS_20k_snps, unique(bivalent_genes$ensembl_id), name="bivalent_gene")
    
    # For lasso.R
    #cases <- CDH_HE_TES_20k_snps
    #controls <- SSC_HE_TES_20k_snps
    
    #CDH_FANTOM <- CDH_HE_TES_20k_snps$FANTOM5_enhancer_robust
    #CDH_FANTOM <- CDH_FANTOM[!is.na(CDH_FANTOM)]
    #SSC_FANTOM <- SSC_HE_TES_20k_snps$FANTOM5_enhancer_robust
    #SSC_FANTOM <- SSC_FANTOM[!is.na(SSC_FANTOM)]
    #ft <- fisher.test(matrix(c(sum(CDH_FANTOM == "Y"), sum(SSC_FANTOM == "Y"), sum(CDH_FANTOM == "N"), sum(SSC_FANTOM == "N")), nrow=2, ncol=2), alternative="greater")
    #print(paste0(round(ft$estimate, 3), ", p-value: ", round(ft$p.value, 4)))
    
    # Use cdh_snps_TES$ANNOVAR_ensembl_Closest_gene to find only closest gene, not genes on both sides of the variant (predicted ANNOVAR consequence)!!!!!!
    # Allow for vector of multiple constraints!!!
    annotate_wgsa_gene_names <- function(dat, constraint="pLI>0.5", include_counts=TRUE) {
        return_env <- new.env()
        return_env[["genes"]] <- ""
        return_env[["constrained_genes"]] <- ""
        return_env[["mouse_candidate_genes"]] <- ""
        return_env[["recurrent_genes"]] <- ""
        return(return_env)
        
        #gene_annots <- strsplit(paste0(dat$ANNOVAR_ensembl_Gene_ID), "\\|")
        gene_annots <- strsplit(paste0(dat$TS_gene), ",")
        gene_names_env <- new.env()
        for(j in 1:length(gene_annots)) {
            for(gene_name in unlist(gene_annots[j])) { 
                if(is.null(gene_names_env[[paste(gene_name)]])) {
                    gene_names_env[[paste(gene_name)]] <- 1 
                } else { gene_names_env[[paste(gene_name)]] <- gene_names_env[[paste(gene_name)]] + 1 }
            }
            
            #already_processed <- new.env() 
            #for(ensembl_id in unlist(gene_annots[j])) { 
            #	if(!is.null(already_processed[[paste(ensembl_id)]]) || is.null(ensembl_genes[[paste(ensembl_id)]])) { 
            #		next 
            #	} else {
            #		if(is.null(gene_names_env[[ensembl_genes[[paste(ensembl_id)]]]])) { 
            #			gene_names_env[[ensembl_genes[[paste(ensembl_id)]]]] <- 1 
            #		} else { gene_names_env[[ensembl_genes[[paste(ensembl_id)]]]] <- gene_names_env[[ensembl_genes[[paste(ensembl_id)]]]] + 1 }
            #		already_processed[[paste(ensembl_id)]] <- 1
            #	}
            #}
        }
        
        constrained_genes <- get_constrained_genes(constraint)
        mouse_candidate_genes <- get_cdh_mouse_candidate_genes()
        
        genes_vector <- c()
        constrained_genes_vector <- c()
        mouse_candidate_genes_vector <- c()
        recurrent_genes_vector <- c()
        if (include_counts) {
            for(gene in ls(gene_names_env)) {
                genes_vector <- c(genes_vector, paste0(gene, "(", gene_names_env[[gene]], ")"))
                if(!is.null(constrained_genes[[gene]])) { constrained_genes_vector <- c(constrained_genes_vector, paste0(gene, "(", gene_names_env[[gene]], ")")) }
                if(!is.null(mouse_candidate_genes[[gene]])) { mouse_candidate_genes_vector <- c(mouse_candidate_genes_vector, paste0(gene, "(", gene_names_env[[gene]], ")")) }
                if(as.numeric(gene_names_env[[gene]]) > 1) { recurrent_genes_vector <- c(recurrent_genes_vector, paste0(gene, "(", gene_names_env[[gene]], ")")) }
            }
        } else {
            genes_vector <- ls(gene_names_env) 
            for(gene in ls(gene_names_env)) {
                if(!is.null(constrained_genes[[gene]])) { constrained_genes_vector <- c(constrained_genes_vector, gene) }
                if(!is.null(mouse_candidate_genes[[gene]])) { mouse_candidate_genes_vector <- c(mouse_candidate_genes_vector, gene) }
            }
        }
        gene_names <- paste0(genes_vector, collapse=",")
        constrained_gene_names <- paste0(constrained_genes_vector, collapse=",")
        mouse_candidate_gene_names <- paste0(mouse_candidate_genes_vector, collapse=",")
        recurrent_gene_names <- paste0(recurrent_genes_vector, collapse=",")
        if(gene_names == "") { gene_names <- "-" }
        if(constrained_gene_names == "") { constrained_gene_names <- "-" }
        if(mouse_candidate_gene_names == "") { mouse_candidate_gene_names <- "-" }
        if(recurrent_gene_names == "") { recurrent_gene_names <- "-" }
        
        return_env <- new.env()
        return_env[["genes"]] <- gene_names
        return_env[["constrained_genes"]] <- constrained_gene_names
        return_env[["mouse_candidate_genes"]] <- mouse_candidate_gene_names
        return_env[["recurrent_genes"]] <- recurrent_gene_names
        
        return(return_env)
    }
    
    # Analyze burden for all annotation variables.
    analyze_burden <- function(cases_annotation, controls_annotation, n1=0, n0=0) {
        if(n0 == 0) { n0 <- length(unique(controls_annotation$sample)) }
        if(n1 == 0) { n1 <- length(unique(cases_annotation$sample)) }
        
        burden_D_qual <- data.frame()
        burden_D_quant <- data.frame()
        shared_features <- intersect(colnames(cases_annotation), colnames(controls_annotation))
        constraint="pLI>0.5"
        
        for(i in 1:length(shared_features)) {
            feature_name <- shared_features[i]
            print(paste0(i, "    ", feature_name))
            CDH_feature_all <- cases_annotation[,c(feature_name)]
            CDH_feature <- CDH_feature_all[!is.na(CDH_feature_all)]
            SSC_feature_all <- controls_annotation[,c(feature_name)]
            SSC_feature <- SSC_feature_all[!is.na(SSC_feature_all)]
            
            CDH_feature_num <- get_numeric(CDH_feature)
            SSC_feature_num <- get_numeric(SSC_feature)
            
            if (length(CDH_feature_num) < 2 || length(SSC_feature_num) < 2) { # Qualitative feature
                #ft <- fisher.test(matrix(c(sum(CDH_feature == "Y"), sum(SSC_feature == "Y"), sum(CDH_feature == "N"), sum(SSC_feature == "N")), nrow=2, ncol=2), alternative="greater")
                ft <- binomial_test(sum(CDH_feature == "Y"), sum(SSC_feature == "Y"), n1, n0, alternative=c("greater")) # two.sided
                if (ft$p.value < 0.1) {
                    cases_genes_env <- annotate_wgsa_gene_names(cases_annotation[CDH_feature == "Y",], constraint)
                    cases_genes <- cases_genes_env[["genes"]]
                    cases_constrained_genes <- cases_genes_env[["constrained_genes"]]
                    cases_mouse_candidates <- cases_genes_env[["mouse_candidate_genes"]]
                    cases_recurrent_genes <- cases_genes_env[["recurrent_genes"]]
                    controls_genes_env <- annotate_wgsa_gene_names(controls_annotation[SSC_feature == "Y",], constraint)
                    controls_genes <- controls_genes_env[["genes"]]
                    controls_constrained_genes <- controls_genes_env[["constrained_genes"]]	
                    controls_mouse_candidates <- controls_genes_env[["mouse_candidate_genes"]]
                    controls_recurrent_genes <- controls_genes_env[["recurrent_genes"]]
                } else { cases_genes = "-"; cases_constrained_genes = "-"; cases_mouse_candidates = "-"; cases_recurrent_genes = "-"; controls_genes = "-"; controls_constrained_genes = "-"; controls_mouse_candidates = "-"; controls_recurrent_genes = "-" }
                burden_D_qual <- rbind(burden_D_qual, cbind(feature_name, sum(CDH_feature == "Y"), sum(CDH_feature == "N"), sum(SSC_feature == "Y"), sum(SSC_feature == "N"), ft$estimate, ft$p.value, cases_genes, cases_constrained_genes, cases_mouse_candidates, cases_recurrent_genes, controls_genes, controls_constrained_genes, controls_mouse_candidates, controls_recurrent_genes)) #sum(is.na(CDH_feature_all))
            } else { # Quantitative feature
                #enrichment = mean(CDH_feature_num) / mean(SSC_feature_num)
                #tt <- try(t.test(CDH_feature_num, SSC_feature_num))
                #if(inherits(tt, "try-error")) { next }
                ##p.value = pnorm(mean(CDH_feature_num), mean=mean(SSC_feature_num), sd=sd(SSC_feature_num), lower.tail=FALSE)
                #burden_D_quant <- rbind(burden_D_quant, cbind(feature_name, mean(CDH_feature_num), mean(SSC_feature_num), length(CDH_feature_num), length(SSC_feature_num), enrichment, tt$p.value))
            } 
        }
        colnames(burden_D_qual) <- c("feature", "cases_Y", "cases_N", "controls_Y", "controls_N", "enrichment", "p.value", "cases_genes", paste0("cases_constrained_genes(",constraint,")"), "cases_mouse_candidate_genes", "cases_recurrent_genes", "controls_genes", paste0("controls_constrained_genes(",constraint,")"), "controls_mouse_candidate_genes", "controls_recurrent_genes")
        rownames(burden_D_qual) <- seq(1:nrow(burden_D_qual))
        burden_D_qual$cases_Y <- as.numeric(paste(burden_D_qual$cases_Y))
        burden_D_qual$cases_N <- as.numeric(paste(burden_D_qual$cases_N))
        burden_D_qual$controls_Y <- as.numeric(paste(burden_D_qual$controls_Y))
        burden_D_qual$controls_N <- as.numeric(paste(burden_D_qual$controls_N))
        burden_D_qual$enrichment <- as.numeric(paste(burden_D_qual$enrichment))
        burden_D_qual$p.value <- as.numeric(paste(burden_D_qual$p.value))
        
        burden_D_quant
        colnames(burden_D_quant) <- c("feature", "cases_mean", "controls_mean", "cases_count", "controls_count", "enrichment", "p.value")
        rownames(burden_D_quant) <- seq(1:nrow(burden_D_quant))
        burden_D_quant$cases_mean <- as.numeric(paste(burden_D_quant$cases_mean))
        burden_D_quant$controls_mean <- as.numeric(paste(burden_D_quant$controls_mean))
        burden_D_quant$cases_count <- as.numeric(paste(burden_D_quant$cases_count))
        burden_D_quant$controls_count <- as.numeric(paste(burden_D_quant$controls_count))
        burden_D_quant$enrichment <- as.numeric(paste(burden_D_quant$enrichment))
        burden_D_quant$p.value <- as.numeric(paste(burden_D_quant$p.value))
        
        burdens <- new.env()
        burdens[["qualitative"]] <- burden_D_qual	
        burdens[["quantitative"]] <- burden_D_quant
        return(burdens)
    }
    a <- analyze_burden(cdh_muts_TES, ssc_muts_TES, n1=CDH_SAMPLE_COUNT, n0=SSC_SAMPLE_COUNT)
    
    
    CDH_HE_TSS_20k_snps
    CDH_HE_TES_20k_snps
    SSC_HE_TSS_20k_snps
    SSC_HE_TES_20k_snps
    
    TSS_burdens_env <- analyze_burden(CDH_HE_TSS_20k_snps, SSC_HE_TSS_20k_snps)
    TSS_burdens <- TSS_burdens_env[["qualitative"]]
    epigenome_names <- get_roadmap_epigenome_names(unlist(lapply(seq(1:length(TSS_burdens$feature)), function(i) { return(strsplit(paste(TSS_burdens$feature[i]), "\\.")[[1]][1]) } )))
    TSS_burdens <- cbind(TSS_burdens, epigenome_names)
    write.csv(TSS_burdens[TSS_burdens$enrichment > 2,], file=paste0(DISEASE, "_HHE_TSS_burden_analysis_strong.csv"), row.names=FALSE)
    write.csv(TSS_burdens[TSS_burdens$enrichment > 1,], file=paste0(DISEASE, "_HHE_TSS_burden_analysis.csv"), row.names=FALSE)
    TES_burdens_env <- analyze_burden(CDH_HE_TES_20k_snps, SSC_HE_TES_20k_snps)
    TES_burdens <- TES_burdens_env[["qualitative"]]
    epigenome_names <- get_roadmap_epigenome_names(unlist(lapply(seq(1:length(TES_burdens$feature)), function(i) { return(strsplit(paste(TES_burdens$feature[i]), "\\.")[[1]][1]) } )))
    TES_burdens <- cbind(TES_burdens, epigenome_names)
    write.csv(TES_burdens[TES_burdens$enrichment > 2,], file=paste0(DISEASE, "_HHE_TES_burden_analysis_strong.csv"), row.names=FALSE)
    write.csv(TES_burdens[TES_burdens$enrichment > 1,], file=paste0(DISEASE, "_HHE_TES_burden_analysis.csv"), row.names=FALSE)
    
    TES_burdens <- TES_burdens_env[["quantitative"]]
    TES_burdens[1:10,]
    TES_burdens[TES_burdens$p.value < 0.1,]
    
    columbia <- read.csv("CDH_denovo_IGV_noncoding_uniq.csv")
    columbia_fams <- paste(as.matrix(data.frame(strsplit(paste(columbia$proband), "\\("))[1,]))
    columbia <- columbia[!(columbia_fams %in% excluded_samples),]
    columbia_indels <- columbia[nchar(paste(columbia$REF)) > 1 | nchar(paste(columbia$ALT)) > 1,]
    ssc_controls <- read.csv("SSC_denovo_IGV_noncoding.csv")
    ssc_indels <- ssc_controls[nchar(paste(ssc_controls$REF)) > 1 | nchar(paste(ssc_controls$ALT)) > 1,]
    CDH_HE_TSS_20k_indels_genes <- is_TS(columbia_indels$CHROM, columbia_indels$POS, sites="TSS", return_genes=TRUE)
    CDH_HE_TSS_20k_indels <- columbia_indels[CDH_HE_TSS_20k_indels_genes != "",]
    CDH_HE_TSS_20k_indels_genes <- CDH_HE_TSS_20k_indels_genes[CDH_HE_TSS_20k_indels_genes != ""]
    CDH_HE_TES_20k_indels_genes <- is_TS(columbia_indels$CHROM, columbia_indels$POS, sites="TES", return_genes=TRUE)
    CDH_HE_TES_20k_indels <- columbia_indels[CDH_HE_TES_20k_indels_genes != "",]
    CDH_HE_TES_20k_indels_genes <- CDH_HE_TES_20k_indels_genes[CDH_HE_TES_20k_indels_genes != ""]
    SSC_HE_TSS_20k_indels_genes <- is_TS(ssc_indels$CHROM, ssc_indels$POS, sites="TSS", return_genes=TRUE)
    SSC_HE_TSS_20k_indels <- ssc_indels[SSC_HE_TSS_20k_indels_genes != "",]
    SSC_HE_TSS_20k_indels_genes <- SSC_HE_TSS_20k_indels_genes[SSC_HE_TSS_20k_indels_genes != ""]
    SSC_HE_TES_20k_indels_genes <- is_TS(ssc_indels$CHROM, ssc_indels$POS, sites="TES", return_genes=TRUE)
    SSC_HE_TES_20k_indels <- ssc_indels[SSC_HE_TES_20k_indels_genes != "",]
    SSC_HE_TES_20k_indels_genes <- SSC_HE_TES_20k_indels_genes[SSC_HE_TES_20k_indels_genes != ""]
    
}

# Testing enrichment of mutations in TSS and TES HE regions!!!
cases_Y <- apply(CDH_HE_TSS_20k_snps, 1, function(row) { paste0(gsub(" ", "", row[1]), ":", gsub(" ", "", row[2])) } )
#cases_Y <- c(cases_Y, apply(CDH_HE_TSS_20k_indels, 1, function(row) { paste0(gsub(" ", "", row[1]), ":", gsub(" ", "", row[2])) } ))
cases_N <- paste0(gsub(" ", "", columbia[,1]), ":", gsub(" ", "", columbia[,2]))
cases_N <- cases_N[!(cases_N %in% unique(cases_Y))]
controls_Y <- apply(SSC_HE_TSS_20k_snps, 1, function(row) { paste0(gsub(" ", "", row[1]), ":", gsub(" ", "", row[2])) } )
#controls_Y <- c(controls_Y, apply(SSC_HE_TSS_20k_indels, 1, function(row) { paste0(gsub(" ", "", row[1]), ":", gsub(" ", "", row[2])) } ))
controls_N <- paste0(gsub(" ", "", ssc_controls[,1]), ":", gsub(" ", "", ssc_controls[,2]))
controls_N <- controls_N[!(controls_N %in% unique(controls_Y))]
ft <- fisher.test(matrix(c(length(cases_Y), length(controls_Y), length(cases_N), length(controls_N)), nrow=2, ncol=2), alternative="greater")
print("TSS Results: ")
ft	

cases_Y <- apply(CDH_HE_TES_20k_snps, 1, function(row) { paste0(gsub(" ", "", row[1]), ":", gsub(" ", "", row[2])) } )
#cases_Y <- c(cases_Y, apply(CDH_HE_TES_20k_indels, 1, function(row) { paste0(gsub(" ", "", row[1]), ":", gsub(" ", "", row[2])) } ))
cases_N <- paste0(gsub(" ", "", columbia[,1]), ":", gsub(" ", "", columbia[,2]))
cases_N <- cases_N[!(cases_N %in% unique(cases_Y))]
controls_Y <- apply(SSC_HE_TES_20k_snps, 1, function(row) { paste0(gsub(" ", "", row[1]), ":", gsub(" ", "", row[2])) } )
#controls_Y <- c(controls_Y, apply(SSC_HE_TES_20k_indels, 1, function(row) { paste0(gsub(" ", "", row[1]), ":", gsub(" ", "", row[2])) } ))
controls_N <- paste0(gsub(" ", "", ssc_controls[,1]), ":", gsub(" ", "", ssc_controls[,2]))
controls_N <- controls_N[!(controls_N %in% unique(controls_Y))]
ft <- fisher.test(matrix(c(length(cases_Y), length(controls_Y), length(cases_N), length(controls_N)), nrow=2, ncol=2), alternative="greater")
print("TES Results: ")
ft

############## STARTING CHD RUN FROM HERE!!!!! ##############

columbia <- read.csv("CHD_denovo_IGV_filtered_science.csv")
columbia_fams <- paste(as.matrix(data.frame(strsplit(paste(columbia$proband), "\\("))[1,]))
columbia_indels <- columbia[nchar(paste(columbia$REF)) > 1 | nchar(paste(columbia$ALT)) > 1,]
ssc_controls <- read.csv("SSC_denovo_IGV_noncoding.csv")
ssc_indels <- ssc_controls[nchar(paste(ssc_controls$REF)) > 1 | nchar(paste(ssc_controls$ALT)) > 1,]
CDH_HE_TSS_20k_indels_genes <- is_TS(columbia_indels$CHROM, columbia_indels$POS, sites="TSS", return_genes=TRUE)
CDH_HE_TSS_20k_indels <- columbia_indels[CDH_HE_TSS_20k_indels_genes != "",]
CDH_HE_TSS_20k_indels_genes <- CDH_HE_TSS_20k_indels_genes[CDH_HE_TSS_20k_indels_genes != ""]
CDH_HE_TES_20k_indels_genes <- is_TS(columbia_indels$CHROM, columbia_indels$POS, sites="TES", return_genes=TRUE)
CDH_HE_TES_20k_indels <- columbia_indels[CDH_HE_TES_20k_indels_genes != "",]
CDH_HE_TES_20k_indels_genes <- CDH_HE_TES_20k_indels_genes[CDH_HE_TES_20k_indels_genes != ""]
SSC_HE_TSS_20k_indels_genes <- is_TS(ssc_indels$CHROM, ssc_indels$POS, sites="TSS", return_genes=TRUE)
SSC_HE_TSS_20k_indels <- ssc_indels[SSC_HE_TSS_20k_indels_genes != "",]
SSC_HE_TSS_20k_indels_genes <- SSC_HE_TSS_20k_indels_genes[SSC_HE_TSS_20k_indels_genes != ""]
SSC_HE_TES_20k_indels_genes <- is_TS(ssc_indels$CHROM, ssc_indels$POS, sites="TES", return_genes=TRUE)
SSC_HE_TES_20k_indels <- ssc_indels[SSC_HE_TES_20k_indels_genes != "",]
SSC_HE_TES_20k_indels_genes <- SSC_HE_TES_20k_indels_genes[SSC_HE_TES_20k_indels_genes != ""]

CDH_HE_TSS_20k_snps_genes <- is_TS(cases$X.chr, cases$pos, sites="TSS", return_genes=TRUE)
CDH_HE_TSS_20k_snps_genes <- CDH_HE_TSS_20k_snps_genes[CDH_HE_TSS_20k_snps_genes != ""]
CDH_HE_TES_20k_snps_genes <- is_TS(cases$X.chr, cases$pos, sites="TES", return_genes=TRUE)
CDH_HE_TES_20k_snps_genes <- CDH_HE_TES_20k_snps_genes[CDH_HE_TES_20k_snps_genes != ""]
SSC_HE_TSS_20k_snps_genes <- is_TS(controls$X.chr, controls$pos, sites="TSS", return_genes=TRUE)
SSC_HE_TSS_20k_snps_genes <- SSC_HE_TSS_20k_snps_genes[SSC_HE_TSS_20k_snps_genes != ""]
SSC_HE_TES_20k_snps_genes <- is_TS(controls$X.chr, controls$pos, sites="TES", return_genes=TRUE)
SSC_HE_TES_20k_snps_genes <- SSC_HE_TES_20k_snps_genes[SSC_HE_TES_20k_snps_genes != ""]

# Store (CDH_snp_count, CDH_indel_count, SSC_snp_count, SSC_indel_count) pairs for all implicated genes for TSS.
TSS_genes <- new.env()
for(genes in CDH_HE_TSS_20k_snps_genes) {
    for(gene in strsplit(genes, ",")[[1]]) {
        if(!is.null(TSS_genes[[paste(gene)]])) {
            TSS_genes[[paste(gene)]][1] <- as.numeric(TSS_genes[[paste(gene)]][1]) + 1
        } else {
            TSS_genes[[paste(gene)]] <- c(1, 0, 0, 0)
        }
    }
}
for(genes in CDH_HE_TSS_20k_indels_genes) {
    for(gene in strsplit(genes, ",")[[1]]) {
        if(!is.null(TSS_genes[[paste(gene)]])) {
            TSS_genes[[paste(gene)]][2] <- as.numeric(TSS_genes[[paste(gene)]][2]) + 1
        } else {
            TSS_genes[[paste(gene)]] <- c(0, 1, 0, 0)
        }
    }
}
for(genes in SSC_HE_TSS_20k_snps_genes) {
    for(gene in strsplit(genes, ",")[[1]]) {
        if(!is.null(TSS_genes[[paste(gene)]])) {
            TSS_genes[[paste(gene)]][3] <- as.numeric(TSS_genes[[paste(gene)]][3]) + 1
        } else {
            TSS_genes[[paste(gene)]] <- c(0, 0, 1, 0)
        }
    }
}
for(genes in SSC_HE_TSS_20k_indels_genes) {
    for(gene in strsplit(genes, ",")[[1]]) {
        if(!is.null(TSS_genes[[paste(gene)]])) {
            TSS_genes[[paste(gene)]][4] <- as.numeric(TSS_genes[[paste(gene)]][4]) + 1
        } else {
            TSS_genes[[paste(gene)]] <- c(0, 0, 0, 1)
        }
    }
}
TSS_genes_dat <- data.frame()
for(gene in ls(TSS_genes)) {
    counts <- TSS_genes[[paste(gene)]]
    gene_row <- c(gene, as.numeric(counts))
    print(paste(gene_row))
    #readline()
    TSS_genes_dat <- rbind(TSS_genes_dat, gene_row, stringsAsFactors=FALSE)
}
colnames(TSS_genes_dat) <- c("gene", "CDH_snv_count", "CDH_indel_count", "SSC_snv_count", "SSC_indel_count")
TSS_genes_dat$CDH_snv_count <- as.numeric(TSS_genes_dat$CDH_snv_count)
TSS_genes_dat$CDH_indel_count <- as.numeric(TSS_genes_dat$CDH_indel_count)
TSS_genes_dat$SSC_snv_count <- as.numeric(TSS_genes_dat$SSC_snv_count)
TSS_genes_dat$SSC_indel_count <- as.numeric(TSS_genes_dat$SSC_indel_count)
TSS_genes_dat <- merge(TSS_genes_dat, HE_TSS_20k_regions[,c("gene", "exp_mutations")], by="gene")
TSS_genes_dat$exp_mutations <- TSS_genes_dat$exp_mutations * SAMPLE_COUNT
TSS_genes_dat <- cbind(TSS_genes_dat, dpois(TSS_genes_dat$CDH_snv_count, lambda=TSS_genes_dat$exp_mutations))
colnames(TSS_genes_dat)[(ncol(TSS_genes_dat)-1):ncol(TSS_genes_dat)] <- c("exp_snvs", "dpois(CDH_snvs, exp_snvs)")
write.csv(TSS_genes_dat[as.numeric(TSS_genes_dat$CDH_snv_count) + as.numeric(TSS_genes_dat$CDH_indel_count) > 1,], file=paste0(DISEASE, "_HE_TSS_genes.csv"), row.names=FALSE)

# Store (CDH_snp_count, CDH_indel_count, SSC_snp_count, SSC_indel_count) pairs for all implicated genes for TES.
TES_genes <- new.env()
for(genes in CDH_HE_TES_20k_snps_genes) {
    for(gene in strsplit(genes, ",")[[1]]) {
        if(!is.null(TES_genes[[paste(gene)]])) {
            TES_genes[[paste(gene)]][1] <- as.numeric(TES_genes[[paste(gene)]][1]) + 1
        } else {
            TES_genes[[paste(gene)]] <- c(1, 0, 0, 0)
        }
    }
}
for(genes in CDH_HE_TES_20k_indels_genes) {
    for(gene in strsplit(genes, ",")[[1]]) {
        if(!is.null(TES_genes[[paste(gene)]])) {
            TES_genes[[paste(gene)]][2] <- as.numeric(TES_genes[[paste(gene)]][2]) + 1
        } else {
            TES_genes[[paste(gene)]] <- c(0, 1, 0, 0)
        }
    }
}
for(genes in SSC_HE_TES_20k_snps_genes) {
    for(gene in strsplit(genes, ",")[[1]]) {
        if(!is.null(TES_genes[[paste(gene)]])) {
            TES_genes[[paste(gene)]][3] <- as.numeric(TES_genes[[paste(gene)]][3]) + 1
        } else {
            TES_genes[[paste(gene)]] <- c(0, 0, 1, 0)
        }
    }
}
for(genes in SSC_HE_TES_20k_indels_genes) {
    for(gene in strsplit(genes, ",")[[1]]) {
        if(!is.null(TES_genes[[paste(gene)]])) {
            TES_genes[[paste(gene)]][4] <- as.numeric(TES_genes[[paste(gene)]][4]) + 1
        } else {
            TES_genes[[paste(gene)]] <- c(0, 0, 0, 1)
        }
    }
}
TES_genes_dat <- data.frame()
for(gene in ls(TES_genes)) {
    counts <- TES_genes[[paste(gene)]]
    gene_row <- c(gene, as.numeric(counts))
    TES_genes_dat <- rbind(TES_genes_dat, gene_row, stringsAsFactors=FALSE)
}
colnames(TES_genes_dat) <- c("gene", "CDH_snv_count", "CDH_indel_count", "SSC_snv_count", "SSC_indel_count")
TES_genes_dat$CDH_snv_count <- as.numeric(TES_genes_dat$CDH_snv_count)
TES_genes_dat$CDH_indel_count <- as.numeric(TES_genes_dat$CDH_indel_count)
TES_genes_dat$SSC_snv_count <- as.numeric(TES_genes_dat$SSC_snv_count)
TES_genes_dat$SSC_indel_count <- as.numeric(TES_genes_dat$SSC_indel_count)
TES_genes_dat <- merge(TES_genes_dat, HE_TES_20k_regions[,c("gene", "exp_mutations")], by="gene")
TES_genes_dat$exp_mutations <- TES_genes_dat$exp_mutations * SAMPLE_COUNT
TES_genes_dat <- cbind(TES_genes_dat, dpois(TES_genes_dat$CDH_snv_count, lambda=TES_genes_dat$exp_mutations))
colnames(TES_genes_dat)[(ncol(TES_genes_dat)-1):ncol(TES_genes_dat)] <- c("exp_snvs", "dpois(CDH_snvs, exp_snvs)")

write.csv(TES_genes_dat[as.numeric(TES_genes_dat$CDH_snv_count) + as.numeric(TES_genes_dat$CDH_indel_count) > 1,], file=paste0(DISEASE, "_HE_TES_genes.csv"), row.names=FALSE)

#ft <- fisher.test(matrix(c(sum(CDH_feature == "Y"), sum(SSC_feature == "Y"), sum(CDH_feature == "N"), sum(SSC_feature == "N")), nrow=2, ncol=2), alternative="greater")

write.csv(TSS_genes_dat, file=paste0(DISEASE, "_HE_TSS_all_genes.csv"), row.names=FALSE)
write.csv(TES_genes_dat, file=paste0(DISEASE, "_HE_TES_all_genes.csv"), row.names=FALSE)


# Return gene snv counts for random subsample
# If the parameter unique_genes_per_sample is TRUE (default), count multiple variants in a gene from same individual as a single occurence, not multiple occurences.
individual_variants_hash <- new.env()
subsample_snv_count <- function(dat, sites, sample_size=SAMPLE_COUNT, with_replacement=FALSE, unique_genes_per_sample=TRUE) {
    subsampled_individuals <- paste(sample(unique(dat$sample), size=sample_size, replace=with_replacement)) # subsample with replacement (bootstrap approach)
    
    #dat2 <- dat[unlist(sapply(subsampled_individuals[1:5], function(individual) { if (is.null(individual_variants_hash[[individual]])) { individual_variants_hash[[individual]] <- which(dat$sample == individual) } ;  return(individual_variants_hash[[individual]]) } )), c("X.chr", "pos", "sample")]
    ##dat2 <- dat[rowSums(is.na(dat)) == 0,]
    dat <- dat[dat$sample %in% subsampled_individuals, c("Chrom", "Position", "sample")]
    
    #dat <- dat[unlist(sapply(subsampled_individuals, function(individual) { if (is.null(individual_variants_hash[[individual]])) { individual_variants_hash[[individual]] <- which(dat$sample == individual) } ;  return(individual_variants_hash[[individual]]) } )), c("X.chr", "pos", "sample")]
    snps_genes <- is_TS(dat$Chrom, dat$Position, sites=sites, return_genes=TRUE, print_counts=FALSE, samples=dat$sample)
    snps_genes <- snps_genes[snps_genes != ""]
    sites_genes <- new.env()
    for(genes in snps_genes) {
        for(gene in strsplit(genes, ",")[[1]]) {
            if(!is.null(sites_genes[[paste(gene)]])) {
                sites_genes[[paste(gene)]] <- as.numeric(sites_genes[[paste(gene)]]) + 1
            } else {
                sites_genes[[paste(gene)]] <- 1
            }
        }
    }
    sites_genes_dat <- data.frame()
    for(gene in ls(sites_genes)) {
        count <- as.numeric(sites_genes[[paste(gene)]])
        gene_row <- c(gene, count)
        sites_genes_dat <- rbind(sites_genes_dat, gene_row, stringsAsFactors=FALSE)
    }
    colnames(sites_genes_dat) <- c("gene", "snv_count")
    sites_genes_dat$snv_count <- as.numeric(sites_genes_dat$snv_count)
    return(sites_genes_dat)
}
num_subsample_iterations <- 10
over_1_snv_SSC_subsamples_TSS <- sapply(seq(1:num_subsample_iterations), function(x) { print(paste0("Subsample iteration: ", x)); sub_snv_counts <- subsample_snv_count(controls, "TSS", with_replacement=FALSE); return(sum(sub_snv_counts$snv_count > 1)) } )
over_1_snv_SSC_subsamples_TES <- sapply(seq(1:num_subsample_iterations), function(x) { print(paste0("Subsample iteration: ", x)); sub_snv_counts <- subsample_snv_count(controls, "TES", with_replacement=FALSE); return(sum(sub_snv_counts$snv_count > 1)) } )
over_2_snv_SSC_subsamples_TSS <- sapply(seq(1:num_subsample_iterations), function(x) { print(paste0("Subsample iteration: ", x)); sub_snv_counts <- subsample_snv_count(controls, "TSS", with_replacement=FALSE); return(sum(sub_snv_counts$snv_count > 2)) } )
over_2_snv_SSC_subsamples_TES <- sapply(seq(1:num_subsample_iterations), function(x) { print(paste0("Subsample iteration: ", x)); sub_snv_counts <- subsample_snv_count(controls, "TES", with_replacement=FALSE); return(sum(sub_snv_counts$snv_count > 2)) } )


controls_individual_mut_counts <- sapply(unique(controls$sample), function(individual) { return(sum(controls$sample==individual)) } )
# Return gene snv counts for bootstrap, ie: random subsample with replacement
# If the parameter unique_genes_per_sample is TRUE (default), count multiple variants in a gene from same individual as a single occurence, not multiple occurences.
bootstrap_snv_count <- function(dat, sites, sample_size=SAMPLE_COUNT, unique_genes_per_sample=TRUE) {
    dat <- dat[unlist(sapply(seq(1:sample_size), function(i) { n <- sample(controls_individual_mut_counts, size=1); sample(seq(1:nrow(dat)), size=n, replace=TRUE) } )), c("X.chr", "pos")]
    snps_genes <- is_TS(dat$X.chr, dat$pos, sites=sites, return_genes=TRUE, print_counts=FALSE, samples=dat$sample)
    snps_genes <- snps_genes[snps_genes != ""]
    sites_genes <- new.env()
    for(genes in snps_genes) {
        for(gene in strsplit(genes, ",")[[1]]) {
            if(!is.null(sites_genes[[paste(gene)]])) {
                sites_genes[[paste(gene)]] <- as.numeric(sites_genes[[paste(gene)]]) + 1
            } else {
                sites_genes[[paste(gene)]] <- 1
            }
        }
    }
    sites_genes_dat <- data.frame()
    for(gene in ls(sites_genes)) {
        count <- as.numeric(sites_genes[[paste(gene)]])
        gene_row <- c(gene, count)
        sites_genes_dat <- rbind(sites_genes_dat, gene_row, stringsAsFactors=FALSE)
    }
    colnames(sites_genes_dat) <- c("gene", "snv_count")
    sites_genes_dat$snv_count <- as.numeric(sites_genes_dat$snv_count)
    return(sites_genes_dat)
}

sum(table(unlist(strsplit(paste0(cdh_muts_TSS$TS_gene), ",")))>1)
sum(table(unlist(strsplit(paste0(cdh_muts_TES$TS_gene), ",")))>1)
sum(table(unlist(strsplit(paste0(ssc_muts_TSS$TS_gene), ",")))>1)
sum(table(unlist(strsplit(paste0(ssc_muts_TES$TS_gene), ",")))>1)
sum(table(unlist(strsplit((c(paste0(cdh_muts_TSS$TS_gene), paste0(ssc_muts_TSS$TS_gene))), ",")))>1)
sum(table(unlist(strsplit((c(paste0(cdh_muts_TES$TS_gene), paste0(ssc_muts_TES$TS_gene))), ",")))>1)

# Simple subsampling (without replacement) approach
subsample <- function(full_dat, size, sites, n_greater_than=1, num_iterations=10) {
    samples <- unique(full_dat$sample)
    num_samples <- length(samples)
    
    results <- c()
    for(i in 1:num_iterations) {
        print(paste0("Iteration ", i))
        selected_samples <- sample(samples, size)
        
        dat <- full_dat[full_dat$sample %in% selected_samples,]
        #snps_genes <- is_TS(dat$Chrom, dat$Position, sites=sites, return_genes=TRUE, print_counts=FALSE, samples=dat$sample)
        snps_genes <- is_TS(dat$Chrom, dat$Position, sites=sites, return_genes=TRUE, print_counts=FALSE)
        snps_genes <- snps_genes[snps_genes != ""]
        sites_genes <- new.env()
        for(genes in snps_genes) {
            for(gene in strsplit(genes, ",")[[1]]) {
                if(!is.null(sites_genes[[paste(gene)]])) {
                    sites_genes[[paste(gene)]] <- as.numeric(sites_genes[[paste(gene)]]) + 1
                } else {
                    sites_genes[[paste(gene)]] <- 1
                }
            }
        }
        sites_genes_dat <- data.frame()
        for(gene in ls(sites_genes)) {
            count <- as.numeric(sites_genes[[paste(gene)]])
            gene_row <- c(gene, count)
            sites_genes_dat <- rbind(sites_genes_dat, gene_row, stringsAsFactors=FALSE)
        }
        if(!(nrow(sites_genes_dat) > 0)) { sites_genes_dat <- data.frame(x=numeric(0), y=integer(0)) }
        colnames(sites_genes_dat) <- c("gene", "snv_count")
        sites_genes_dat$snv_count <- as.numeric(sites_genes_dat$snv_count)
        results <- c(results, sum(sites_genes_dat$snv_count > n_greater_than))
    }
    return(results)
}
sizes_to_test <- seq(0, 400, by=25)
CDH_subsample_TSS_tests <- c(0)
CDH_subsample_TES_tests <- c(0)
SSC_subsample_TSS_tests <- c(0)
SSC_subsample_TES_tests <- c(0)
for(curr_size in sizes_to_test) {
    if(curr_size == 0) { next }
    print(paste0("Testing subsample size ", curr_size))
    if(curr_size < SAMPLE_COUNT) {
        CDH_subsample_TES_tests <- c(CDH_subsample_TES_tests, mean(subsample(cases[,c(1:9, which(colnames(cases)=="snv_indel"))], size=curr_size, sites="TES", num_iterations=100)))
        CDH_subsample_TSS_tests <- c(CDH_subsample_TSS_tests, mean(subsample(cases[,c(1:9, which(colnames(cases)=="snv_indel"))], size=curr_size, sites="TSS", num_iterations=100)))
    } else { CDH_subsample_TES_tests <- c(CDH_subsample_TES_tests, NA); CDH_subsample_TSS_tests <- c(CDH_subsample_TSS_tests, NA) }
    SSC_subsample_TES_tests <- c(SSC_subsample_TES_tests, mean(subsample(controls[,c(1:9, which(colnames(controls)=="snv_indel"))], size=curr_size, sites="TES", num_iterations=100)))
    SSC_subsample_TSS_tests <- c(SSC_subsample_TSS_tests, mean(subsample(controls[,c(1:9, which(colnames(controls)=="snv_indel"))], size=curr_size, sites="TSS", num_iterations=100)))
}
lines_x <- new.env(); lines_y <- new.env(); cols <- list(); ltys <- list(); exclude_from_legend <- new.env(); sizes_to_test2 <- sizes_to_test**2
quadratic.model <- lm(CDH_subsample_TES_tests ~ sizes_to_test2 + 0); label <- paste0("CDH TES (y = ", round(quadratic.model$coefficients[1],5), "x^2)"); lines_x[[label]] <- sizes_to_test; lines_y[[label]] <- CDH_subsample_TES_tests; cols[[label]] <- "red"; ltys[[label]] <- 1
lines_x[["quad1"]] <- sizes_to_test; lines_y[["quad1"]] <- predict(quadratic.model); cols[["quad1"]] <- 1; ltys[["quad1"]] <- 2; exclude_from_legend[["quad1"]] <- 1
quadratic.model <- lm(CDH_subsample_TSS_tests ~ sizes_to_test2 + 0); label <- paste0("CDH TSS (y = ", round(quadratic.model$coefficients[1],5), "x^2)"); lines_x[[label]] <- sizes_to_test; lines_y[[label]] <- CDH_subsample_TSS_tests; cols[[label]] <- "orange"; ltys[[label]] <- 1
lines_x[["quad2"]] <- sizes_to_test; lines_y[["quad2"]] <- predict(quadratic.model); cols[["quad2"]] <- 1; ltys[["quad2"]] <- 2; exclude_from_legend[["quad2"]] <- 1
quadratic.model <- lm(SSC_subsample_TES_tests ~ sizes_to_test2 + 0); label <- paste0("SSC TES (y = ", round(quadratic.model$coefficients[1],5), "x^2)"); lines_x[[label]] <- sizes_to_test; lines_y[[label]] <- SSC_subsample_TES_tests; cols[[label]] <- "blue"; ltys[[label]] <- 1
lines_x[["quad3"]] <- sizes_to_test; lines_y[["quad3"]] <- predict(quadratic.model); cols[["quad3"]] <- 1; ltys[["quad3"]] <- 2; exclude_from_legend[["quad3"]] <- 1
quadratic.model <- lm(SSC_subsample_TSS_tests ~ sizes_to_test2 + 0); label <- paste0("SSC TSS (y = ", round(quadratic.model$coefficients[1],5), "x^2)"); lines_x[[label]] <- sizes_to_test; lines_y[[label]] <- SSC_subsample_TSS_tests; cols[[label]] <- "green"; ltys[[label]] <- 1
lines_x[["quad4"]] <- sizes_to_test; lines_y[["quad4"]] <- predict(quadratic.model); cols[["quad4"]] <- 1; ltys[["quad4"]] <- 2; exclude_from_legend[["quad4"]] <- 1
#label <- paste0("y = ", round(quadratic.model$coefficients[1],3), "x^2"); lines_x[[label]] <- sizes_to_test; lines_y[[label]] <- predict(quadratic.model, type="response"); cols[[label]] <- "black"; ltys[[label]] <- 2
multi_plot(lines_x=lines_x, lines_y=lines_y, exclude_from_legend=exclude_from_legend, cols=cols, ltys=ltys, lwd=2, cex.lab=1.8, cex.axis=2, cex.main=1.5, mtext_cex=1.5, mar=1.1*par("mar"), main="Subsampling Recurrence Analysis", mtext="Mean of 100 iterations/size", xlab="subsample size", ylab="n regions mutated in >1 sample", legend_location="topleft", legend_cex=1.3, file="recurrence_subsample_analysis.pdf")


# Permutation approach
permutation <- function(cases_dat, controls_dat, sites, n_greater_than=1, num_iterations=10) {
    #cases_dat_size = nrow(cases_dat)
    #controls_dat_size = nrow(controls_dat)
    case_samples <- unique(cases_dat$sample)
    control_samples <- unique(controls_dat$sample)
    case_num_samples <- length(case_samples)
    control_num_samples <- length(control_samples)
    
    results <- c()
    for(i in 1:num_iterations) {
        print(paste0("Iteration ", i))
        #indices <- sample(seq(1,cases_dat_size+controls_dat_size), cases_dat_size)
        #cases_dat_indices <- indices[indices <= cases_dat_size]
        #controls_dat_indices <- indices[indices > cases_dat_size] - cases_dat_size
        #dat <- rbind(cases_dat[cases_dat_indices,], controls_dat[controls_dat_indices,])
        indices <- sample(seq(1,case_num_samples+control_num_samples), case_num_samples)
        cases_sample_indices <- indices[indices <= case_num_samples]
        controls_sample_indices <- indices[indices > case_num_samples] - case_num_samples
        dat <- rbind(cases_dat[cases_dat$sample %in% case_samples[cases_sample_indices],], controls_dat[controls_dat$sample %in% control_samples[controls_sample_indices],])
        #snps_genes <- is_TS(dat$Chrom, dat$Position, sites=sites, return_genes=TRUE, print_counts=FALSE, samples=dat$sample)
        snps_genes <- is_TS(dat$Chrom, dat$Position, sites=sites, return_genes=TRUE, print_counts=FALSE)
        snps_genes <- snps_genes[snps_genes != ""]
        sites_genes <- new.env()
        for(genes in snps_genes) {
            for(gene in strsplit(genes, ",")[[1]]) {
                if(!is.null(sites_genes[[paste(gene)]])) {
                    sites_genes[[paste(gene)]] <- as.numeric(sites_genes[[paste(gene)]]) + 1
                } else {
                    sites_genes[[paste(gene)]] <- 1
                }
            }
        }
        sites_genes_dat <- data.frame()
        for(gene in ls(sites_genes)) {
            count <- as.numeric(sites_genes[[paste(gene)]])
            gene_row <- c(gene, count)
            sites_genes_dat <- rbind(sites_genes_dat, gene_row, stringsAsFactors=FALSE)
        }
        if(!(nrow(sites_genes_dat) > 0)) { sites_genes_dat <- data.frame(x=numeric(0), y=integer(0)) }
        colnames(sites_genes_dat) <- c("gene", "snv_count")
        sites_genes_dat$snv_count <- as.numeric(sites_genes_dat$snv_count)
        results <- c(results, sum(sites_genes_dat$snv_count > n_greater_than))
    }
    return(results)
}
snv_indel_indices <- is_indel(cases$Ref, cases$Alt)
snv_indel <- rep("snv", nrow(cases))
snv_indel[snv_indel_indices] <- "indel"
cases <- cbind(cases, snv_indel)
snv_indel_indices <- is_indel(controls$Ref, controls$Alt)
snv_indel <- rep("snv", nrow(controls))
snv_indel[snv_indel_indices] <- "indel"
controls <- cbind(controls, snv_indel)
CDH_permutations_TSS_snv <- permutation(cases[cases$snv_indel == "snv",c(1:9, which(colnames(cases)=="snv_indel"))], controls[controls$snv_indel == "snv",c(1:9, which(colnames(controls)=="snv_indel"))], sites="TSS", num_iterations=1000)
CDH_permutations_TES_snv <- permutation(cases[cases$snv_indel == "snv",c(1:9, which(colnames(cases)=="snv_indel"))], controls[controls$snv_indel == "snv",c(1:9, which(colnames(controls)=="snv_indel"))], sites="TES", num_iterations=1000)
CDH_permutations_TSS_mut <- permutation(cases[,c(1:9, which(colnames(cases)=="snv_indel"))], controls[,c(1:9, which(colnames(controls)=="snv_indel"))], sites="TSS", num_iterations=10000)
CDH_permutations_TES_mut <- permutation(cases[,c(1:9, which(colnames(cases)=="snv_indel"))], controls[,c(1:9, which(colnames(controls)=="snv_indel"))], sites="TES", num_iterations=10000)

SSC_permutations_TSS_mut <- permutation(controls[,c(1:9, which(colnames(controls)=="snv_indel"))], cases[,c(1:9, which(colnames(cases)=="snv_indel"))], sites="TSS", num_iterations=1000)
SSC_permutations_TES_mut <- permutation(controls[,c(1:9, which(colnames(controls)=="snv_indel"))], cases[,c(1:9, which(colnames(cases)=="snv_indel"))], sites="TES", num_iterations=1000)
SSC_TSS_genes_dat <- all_TSS_genes_dat[,c(1,4,5,2,3,6)]
colnames(SSC_TSS_genes_dat) <- colnames(all_TSS_genes_dat)
SSC_TES_genes_dat <- all_TES_genes_dat[,c(1,4,5,2,3,6)]
colnames(SSC_TES_genes_dat) <- colnames(all_TES_genes_dat)
plot_recurrance(cases_dat=SSC_TSS_genes_dat, mut_type="mut", main=paste0("SSC TSS Region SNV+Indel Recurrence Analysis"), n_greater_than=1, cases_permutations=SSC_permutations_TSS_mut, filename=paste0("SSC_TSS_over_1_nc_mut.pdf"))
plot_recurrance(cases_dat=SSC_TES_genes_dat, mut_type="mut", main=paste0("SSC 3'UTR Region SNV+Indel Recurrence Analysis"), n_greater_than=1, cases_permutations=SSC_permutations_TES_mut, filename=paste0("SSC_TES_over_1_nc_mut.pdf"))


# Jackknife approach
# If the parameter unique_genes_per_sample is TRUE (default), count multiple variants in a gene from same individual as a single occurence, not multiple occurences.
jackknife_snv_count <- function(dat, sites, block_size=1, sample_size=NULL, n_greater_than=1, num_subsample_iterations_per_jack=1, unique_genes_per_sample=TRUE) {
    individuals <- unique(dat$sample)
    num_individuals <- length(individuals)
    if (!is.null(sample_size)) {
        block_size <- num_individuals - sample_size
    } else {
        sample_size <- num_individuals - block_size
    }
    num_blocks = ceiling(num_individuals/block_size)
    
    cat(paste0("Jackknife with ", num_blocks, " Blocks"), "\n")
    cat("Subsampling all individuals:\n")
    result_all <- mean(sapply(seq(1:num_subsample_iterations_per_jack), function(x) { print(paste0("Subsample iteration: ", x)); sub_snv_counts <- subsample_snv_count(dat, sites, sample_size=num_individuals, with_replacement=FALSE, unique_genes_per_sample= unique_genes_per_sample); return(sum(sub_snv_counts$snv_count > n_greater_than)) } ))
    
    pseudo <- numeric(num_blocks)
    for(block in 1:num_blocks) {
        block_start_index = (block-1)*block_size + 1
        block_end_index = block_start_index + block_size - 1
        extra_samples_in_block_indices <- c()
        if (block == num_blocks) { # Last block, so process individuals until the end and then randomly sample enough extra from rest of data to top off count for this block!
            extra_samples_in_block_indices <- sample(seq(1:(block_start_index-1)), (block_end_index - num_individuals))
            block_end_index = num_individuals
        }
        #cat(paste0("Subsampling individuals excluding Block ", block, " (", block_start_index, "-", block_end_index, "): "), "\n")
        subsampled_individuals <- c(paste0(individuals[-c(block_start_index:block_end_index, extra_samples_in_block_indices)]))
        result_jack <- mean(sapply(seq(1:num_subsample_iterations_per_jack), function(x) { print(paste0("Subsample iteration: ", x)); sub_snv_counts <- subsample_snv_count(dat[dat$sample %in% subsampled_individuals, c("Chrom", "Position", "sample")], sites, sample_size=sample_size, with_replacement=FALSE, unique_genes_per_sample= unique_genes_per_sample); return(sum(sub_snv_counts$snv_count > n_greater_than)) } ))
        #pseudo[block] <- num_blocks*result_all - (num_blocks-1)*result_jack
        pseudo[block] <- result_jack
        print(paste0("Block: ", block, ", Result Jack: ", result_jack, ", Pseudo: ", pseudo[block]))
    }
    jackknife_estimate <- mean(pseudo)
    jackknife_var <- (num_blocks-1)/num_blocks * sum((jackknife_estimate-pseudo)**2)
    #jackknife_ci <- c(mean(pseudo) - qt(0.975,num_blocks-1)*sqrt(var(pseudo)/num_blocks), mean(pseudo) + qt(0.975,num_blocks-1)*sqrt(var(pseudo)/num_blocks))
    
    jackknife_return <- new.env()
    jackknife_return[["pseudo"]] <- pseudo
    jackknife_return[["estimate"]] <- jackknife_estimate
    jackknife_return[["var"]] <- jackknife_var
    #jackknife_return[["ci"]] <- jackknife_ci
    return(jackknife_return)
}

many_jackknife <- function(dat, sites, n_greater_than=1, num_subsample_iterations=2, sample_count=SAMPLE_COUNT) {
    jack_estimates <- c()
    jack_variances <- c()
    for(y in 1:num_subsample_iterations) { 
        print(paste0("Iteration ", y))
        dat_subsampled <- dat[dat$sample %in% sample(unique(paste(dat$sample)), sample_count+1),]
        jack <- jackknife_snv_count(dat_subsampled, sites=sites, block_size=1, n_greater_than=n_greater_than)
        jack_estimates <- c(jack_estimates, jack[["estimate"]])
        jack_variances <- c(jack_variances, jack[["var"]])
    }
    many_jack <- new.env()
    many_jack[["estimates"]] <- jack_estimates
    many_jack[["vars"]] <- jack_variances
    return(many_jack)
    
}

# Many Jackknife
num_subsample_iterations = 10
controls_snvs <- controls[controls$snv_indel == "snv",]
over_1_snv_SSC_many_jackknife_TSS <- many_jackknife(controls_snvs, sites="TSS", n_greater_than=1, num_subsample_iterations=num_subsample_iterations)
over_1_snv_SSC_many_jackknife_TES <- many_jackknife(controls_snvs, sites="TES", n_greater_than=1, num_subsample_iterations=num_subsample_iterations)
over_1_mut_SSC_many_jackknife_TSS <- many_jackknife(controls, sites="TSS", n_greater_than=1, num_subsample_iterations=num_subsample_iterations)
over_1_mut_SSC_many_jackknife_TES <- many_jackknife(controls, sites="TES", n_greater_than=1, num_subsample_iterations=num_subsample_iterations)
over_2_snv_SSC_many_jackknife_TSS <- many_jackknife(controls_snvs, sites="TSS", n_greater_than=2, num_subsample_iterations=num_subsample_iterations)
over_2_snv_SSC_many_jackknife_TES <- many_jackknife(controls_snvs, sites="TES", n_greater_than=2, num_subsample_iterations=num_subsample_iterations)
over_2_mut_SSC_many_jackknife_TSS <- many_jackknife(controls, sites="TSS", n_greater_than=2, num_subsample_iterations=num_subsample_iterations)
over_2_mut_SSC_many_jackknife_TES <- many_jackknife(controls, sites="TES", n_greater_than=2, num_subsample_iterations=num_subsample_iterations)

# CDH 195 samples Many Jackknife
over_1_snv_SSC_many_jackknife_TSS_195 <- many_jackknife(controls_snvs, sites="TSS", n_greater_than=1, num_subsample_iterations=num_subsample_iterations, sample_count=195)
over_1_snv_SSC_many_jackknife_TES_195 <- many_jackknife(controls_snvs, sites="TES", n_greater_than=1, num_subsample_iterations=num_subsample_iterations, sample_count=195)
over_1_mut_SSC_many_jackknife_TSS_195 <- many_jackknife(controls, sites="TSS", n_greater_than=1, num_subsample_iterations=num_subsample_iterations, sample_count=195)
over_1_mut_SSC_many_jackknife_TES_195 <- many_jackknife(controls, sites="TES", n_greater_than=1, num_subsample_iterations=num_subsample_iterations, sample_count=195)
#over_2_snv_SSC_many_jackknife_TSS_195 <- many_jackknife(controls_snvs, sites="TSS", n_greater_than=2, num_subsample_iterations=num_subsample_iterations, sample_count=195)
#over_2_snv_SSC_many_jackknife_TES_195 <- many_jackknife(controls_snvs, sites="TES", n_greater_than=2, num_subsample_iterations=num_subsample_iterations, sample_count=195)
#over_2_mut_SSC_many_jackknife_TSS_195 <- many_jackknife(controls, sites="TSS", n_greater_than=2, num_subsample_iterations=num_subsample_iterations, sample_count=195)
#over_2_mut_SSC_many_jackknife_TES_195 <- many_jackknife(controls, sites="TES", n_greater_than=2, num_subsample_iterations=num_subsample_iterations, sample_count=195)

# Recurrance plots
#plot_recurrance(cases_dat=all_TSS_genes_dat, mut_type="snv", main=paste0(DISEASE, " TSS Region SNV Recurrence Analysis"), n_greater_than=1, cases_sim_dat=sim_dat_CDH_TSS, cases_permutations=CHD_permutations_TSS, controls_many_jackknife=TSS_many_jack_snv, filename=paste0(DISEASE, "_TSS_over_1_nc_snv.pdf"))
#plot_recurrance(cases_dat=all_TES_genes_dat, mut_type="snv", main=paste0(DISEASE, " 3'UTR Region SNV Recurrence Analysis"), n_greater_than=1, cases_sim_dat=sim_dat_CDH_TES, cases_permutations=CHD_permutations_TES, controls_many_jackknife=TES_many_jack_snv, filename=paste0(DISEASE, "_TES_over_1_nc_snv.pdf"))
plot_recurrance(cases_dat=all_TSS_genes_dat, mut_type="snv", main=paste0(DISEASE, " TSS Region SNV Recurrence Analysis"), n_greater_than=1, cases_permutations=CHD_permutations_TSS_snv, controls_many_jackknife=TSS_many_jack_snv, filename=paste0(DISEASE, "_TSS_over_1_nc_snv.pdf"))
plot_recurrance(cases_dat=all_TES_genes_dat, mut_type="snv", main=paste0(DISEASE, " 3'UTR Region SNV Recurrence Analysis"), n_greater_than=1, cases_permutations=CHD_permutations_TES_snv, controls_many_jackknife=TES_many_jack_snv, filename=paste0(DISEASE, "_TES_over_1_nc_snv.pdf"))
plot_recurrance(cases_dat=all_TSS_genes_dat, mut_type="mut", main=paste0(DISEASE, " TSS Region SNV+Indel Recurrence Analysis"), n_greater_than=1, cases_permutations=CHD_permutations_TSS_mut, controls_many_jackknife=TSS_many_jack_mut, filename=paste0(DISEASE, "_TSS_over_1_nc_mut.pdf"))
plot_recurrance(cases_dat=all_TES_genes_dat, mut_type="mut", main=paste0(DISEASE, " 3'UTR Region SNV+Indel Recurrence Analysis"), n_greater_than=1, cases_permutations=CHD_permutations_TES_mut, controls_many_jackknife=TES_many_jack_mut, filename=paste0(DISEASE, "_TES_over_1_nc_mut.pdf"))
#plot_recurrance(cases_dat=all_TSS_genes_dat, mut_type="mut", main=paste0(DISEASE, " TSS Region SNV+Indel Recurrence Analysis"), n_greater_than=1, controls_many_jackknife=TSS_many_jack_mut, filename=paste0(DISEASE, "_TSS_over_1_nc_mut.pdf"))
#plot_recurrance(cases_dat=all_TES_genes_dat, mut_type="mut", main=paste0(DISEASE, " 3'UTR Region SNV+Indel Recurrence Analysis"), n_greater_than=1, controls_many_jackknife=TES_many_jack_mut, filename=paste0(DISEASE, "_TES_over_1_nc_mut.pdf"))
#plot_recurrance(cases_dat=all_TSS_genes_dat, mut_type="snv", main=paste0(DISEASE, " TSS Region SNV Recurrence Analysis"), n_greater_than=2, cases_sim_dat=sim_dat_CDH_TSS, controls_many_jackknife=over_2_snv_SSC_many_jackknife_TSS, filename=paste0(DISEASE, "_TSS_over_2_nc_snv.pdf"))
#plot_recurrance(cases_dat=all_TES_genes_dat, mut_type="snv", main=paste0(DISEASE, " 3'UTR Region SNV Recurrence Analysis"), n_greater_than=2, cases_sim_dat=sim_dat_CDH_TES, controls_many_jackknife=over_2_snv_SSC_many_jackknife_TES, filename=paste0(DISEASE, "_TES_over_2_nc_snv.pdf"))
#plot_recurrance(cases_dat=all_TSS_genes_dat, mut_type="mut", main=paste0(DISEASE, " TSS Region SNV+Indel Recurrence Analysis"), n_greater_than=2, controls_many_jackknife=over_2_mut_SSC_many_jackknife_TSS, filename=paste0(DISEASE, "_TSS_over_2_nc_mut.pdf"))
#plot_recurrance(cases_dat=all_TES_genes_dat, mut_type="mut", main=paste0(DISEASE, " 3'UTR Region SNV+Indel Recurrence Analysis"), n_greater_than=2, controls_many_jackknife=over_2_mut_SSC_many_jackknife_TES, filename=paste0(DISEASE, "_TES_over_2_nc_mut.pdf"))

# Controls recurrence plots
controls_sim_dat_TSS <- simulate_mutations(HE_TSS_20k_regions_backup, num_individuals=438)
controls_sim_dat_TES <- simulate_mutations(HE_TES_20k_regions_backup, num_individuals=438)
DISEASE="SSC"
SAMPLE_COUNT=438
controls_TSS_genes_dat <- all_TSS_genes_dat[,c(1,4,5,2,3)]
controls_TES_genes_dat <- all_TES_genes_dat[,c(1,4,5,2,3)]
colnames(controls_TSS_genes_dat)[1:5] <- colnames(all_TSS_genes_dat)[1:5]
colnames(controls_TES_genes_dat)[1:5] <- colnames(all_TES_genes_dat)[1:5]
plot_recurrance(cases_dat=controls_TSS_genes_dat, mut_type="snv", main=paste0("SSC TSS Region SNV Recurrence Analysis"), n_greater_than=1, cases_sim_dat=controls_sim_dat_TSS, filename=paste0(DISEASE, "_TSS_over_1_nc_snv.pdf"))
plot_recurrance(cases_dat=controls_TES_genes_dat, mut_type="snv", main=paste0("SSC 3'UTR Region SNV Recurrence Analysis"), n_greater_than=1, cases_sim_dat=controls_sim_dat_TES, filename=paste0(DISEASE, "_TES_over_1_nc_snv.pdf"))
DISEASE="CHD"
SAMPLE_COUNT=350

# CDH plots
#cases_CDH <- read.csv("CDH_noncoding.txt", sep="\t")
CDH_TSS_genes_dat <- read.csv("CDH_HDE_TSS_genes_all.csv")
CDH_TES_genes_dat <- read.csv("CDH_HDE_TES_genes_all.csv")
#cases_CDH_snps <- read.csv("CDH_noncoding.WGSA_annotated.snp", sep="\t")   # WGSA ANNOTATED!!!
#cases_CDH_indels <- read.csv("CDH_noncoding.WGSA_annotated.indel", sep="\t") # WGSA ANNOTATED!!! CDH_noncoding.WGSA_annotated.snp
#snv_indel <- rep("snv", nrow(cases_CDH_snps)); cases_CDH_snps <- cbind(cases_CDH_snps, snv_indel)
#snv_indel <- rep("indel", nrow(cases_CDH_indels)); cases_CDH_indels <- cbind(cases_CDH_indels, snv_indel)
cases_CDH <- rbind(cases_CDH_snps[,c(1:5, which(colnames(cases_CDH_snps) == "snv_indel"))], cases_CDH_indels[,c(1:5, which(colnames(cases_CDH_indels) == "snv_indel"))]) 
cases_CDH <- standardize_colnames(cases_CDH, remove_chr_prefix=TRUE)
excluded_samples <- readLines("CDH_excluded_sample.txt")
cases_CDH <- cases_CDH[!(cases_CDH$sample %in% excluded_samples),]
#colnames(cases_CDH) <- c("Chrom", "Position", "Ref", "Alt", "sample")
CDH_TSS_genes_dat <- build_mutations_dat(cases_CDH$Chrom, cases_CDH$Position, is_indel(cases_CDH$Ref,cases_CDH$Alt), controls$Chrom, controls$Position, is_indel(controls$Ref,controls$Alt), sites="TSS")
CDH_TES_genes_dat <- build_mutations_dat(cases_CDH$Chrom, cases_CDH$Position, is_indel(cases_CDH$Ref,cases_CDH$Alt), controls$Chrom, controls$Position, is_indel(controls$Ref,controls$Alt), sites="TES")

CDH_sim_dat_TSS <- simulate_mutations(HE_TSS_20k_regions, num_individuals=195)
CDH_sim_dat_TES <- simulate_mutations(HE_TES_20k_regions, num_individuals=195)
DISEASE="CDH"
SAMPLE_COUNT=195
#plot_recurrance(cases_dat=CDH_TSS_genes_dat, mut_type="snv", main=paste0("CDH TSS Region SNV Recurrence Analysis"), n_greater_than=1, cases_sim_dat=CDH_sim_dat_TSS, cases_permutations=CDH_permutations_TSS_snv, controls_many_jackknife=over_1_snv_SSC_many_jackknife_TSS_195, filename=paste0("CDH2_TSS_over_1_nc_snv.pdf"))
#plot_recurrance(cases_dat=CDH_TES_genes_dat, mut_type="snv", main=paste0("CDH 3'UTR Region SNV Recurrence Analysis"), n_greater_than=1, cases_sim_dat=CDH_sim_dat_TES, cases_permutations=CDH_permutations_TES_snv, controls_many_jackknife=over_1_snv_SSC_many_jackknife_TES_195, filename=paste0("CDH2_TES_over_1_nc_snv.pdf"))
#plot_recurrance(cases_dat=CDH_TSS_genes_dat, mut_type="mut", main=paste0("CDH TSS Region SNV+Indel Recurrence Analysis"), n_greater_than=1, cases_permutations=CDH_permutations_TSS_mut, controls_many_jackknife=over_1_mut_SSC_many_jackknife_TSS_195, filename=paste0("CDH2_TSS_over_1_nc_mut.pdf"))
#plot_recurrance(cases_dat=CDH_TES_genes_dat, mut_type="mut", main=paste0("CDH 3'UTR Region SNV+Indel Recurrence Analysis"), n_greater_than=1, cases_permutations=CDH_permutations_TES_mut, controls_many_jackknife=over_1_mut_SSC_many_jackknife_TES_195, filename=paste0("CDH2_TES_over_1_nc_mut.pdf"))
plot_recurrance(cases_dat=CDH_TSS_genes_dat, mut_type="snv", main=paste0("CDH TSS Region SNV Recurrence Analysis"), n_greater_than=1, cases_permutations=CDH_permutations_TSS_snv, filename=paste0("CDH2_TSS_over_1_nc_snv.pdf"))
plot_recurrance(cases_dat=CDH_TES_genes_dat, mut_type="snv", main=paste0("CDH 3'UTR Region SNV Recurrence Analysis"), n_greater_than=1, cases_permutations=CDH_permutations_TES_snv, filename=paste0("CDH2_TES_over_1_nc_snv.pdf"))
plot_recurrance(cases_dat=CDH_TSS_genes_dat, mut_type="mut", main=paste0("CDH TSS Region SNV+Indel Recurrence Analysis"), n_greater_than=1, cases_permutations=CDH_permutations_TSS_mut, filename=paste0("CDH2_TSS_over_1_nc_mut.pdf"))
plot_recurrance(cases_dat=CDH_TES_genes_dat, mut_type="mut", main=paste0("CDH 3'UTR Region SNV+Indel Recurrence Analysis"), n_greater_than=1, cases_permutations=CDH_permutations_TES_mut, filename=paste0("CDH2_TES_over_1_nc_mut.pdf"))
SAMPLE_COUNT=350
DISEASE="CHD"


candidate_CHD_TSS_genes_dat <- build_mutations_dat(cases$Chrom, cases$Position, is_indel(cases$Ref,cases$Alt), controls$Chrom, controls$Position, is_indel(controls$Ref,controls$Alt), sites="TSS")
candidate_CHD_TES_genes_dat <- build_mutations_dat(cases$Chrom, cases$Position, is_indel(cases$Ref,cases$Alt), controls$Chrom, controls$Position, is_indel(controls$Ref,controls$Alt), sites="TES")
candidate_CHD_sim_dat_Felix_TSS <- simulate_mutations(HE_TSS_20k_regions, num_individuals=350)
candidate_CHD_sim_dat_TES <- simulate_mutations(HE_TES_20k_regions, num_individuals=350)
plot_recurrance(cases_dat=candidate_CHD_TSS_genes_dat, mut_type="snv", main=paste0("CHD (TSS-10k,TSS) Region SNV Recurrence Analysis"), n_greater_than=1, cases_sim_dat=candidate_CHD_sim_dat_Felix_TSS, filename=paste0("candidate_CHD_TSS_over_1_nc_snv.pdf"))
plot_recurrance(cases_dat=candidate_CHD_TES_genes_dat, mut_type="snv", main=paste0("CHD 3'UTR Region SNV Recurrence Analysis"), n_greater_than=1, cases_sim_dat=candidate_CHD_sim_dat_TES, filename=paste0("candidate_CHD_TES_over_1_nc_snv.pdf"))


# Compare counts!

CHD_TSS_snv_count_cases <- sum(CHD_TSS_genes_dat[,2])
CHD_TSS_snv_count_sim <- mean(apply(CHD_sim_dat_TSS, 1, sum))
cat(paste0("# of CHD TSS SNVs in Cases: ", CHD_TSS_snv_count_cases), "\n", paste0("Mean # of CHD TSS SNVs in Simulation: ", CHD_TSS_snv_count_sim), "\n", sep="")
CHD_TES_snv_count_cases <- sum(CHD_TES_genes_dat[,2])
CHD_TES_snv_count_sim <- mean(apply(CHD_sim_dat_TES, 1, sum))
cat(paste0("# of CHD TES SNVs in Cases: ", CHD_TES_snv_count_cases), "\n", paste0("Mean # of CHD TES SNVs in Simulation: ", CHD_TES_snv_count_sim), "\n", sep="")

CDH_TSS_snv_count_cases <- sum(CDH_TSS_genes_dat[,2])
CDH_TSS_snv_count_sim <- mean(apply(CDH_sim_dat_TSS, 1, sum))
cat(paste0("# of CDH TSS SNVs in Cases: ", CDH_TSS_snv_count_cases), "\n", paste0("Mean # of CDH TSS SNVs in Simulation: ", CDH_TSS_snv_count_sim), "\n", sep="")
CDH_TES_snv_count_cases <- sum(CDH_TES_genes_dat[,2])
CDH_TES_snv_count_sim <- mean(apply(CDH_sim_dat_TES, 1, sum))
cat(paste0("# of CDH TES SNVs in Cases: ", CDH_TES_snv_count_cases), "\n", paste0("Mean # of CDH TES SNVs in Simulation: ", CDH_TES_snv_count_sim), "\n", sep="")



# Find probabilities of number of genes with x number of mutations from simulation.
simulate_mutations <- function(regions, num_simulations=10000, num_individuals=SAMPLE_COUNT, fixed_mutation_count=NULL, print_sim_dat=FALSE) {
    print(paste0("Num individuals: ", num_individuals))
    num_regions <- nrow(regions)
    #print(paste0("Sampling: ", sample_size))
    #sample_indices <- sample(seq(1:num_individuals), size=sample_size)
    #sim_dat <- sapply(seq(1:nrow(regions)), function(i) { print(paste0(regions$gene[i], ": ", regions$exp_mutations[i])); readline(); rpois(num_simulations, lambda=regions$exp_mutations[i]) } )
    sim_dat <- sapply(seq(1:nrow(regions)), function(i) { print(paste0(i, "/", num_regions)); muts <- sapply(seq(1:num_individuals), function(j) { rpois(num_simulations, lambda=regions$exp_mutations[i]) }); return(rowSums(muts>0)) } )
    sim_dat_results <- table(sim_dat)
    ncells <- nrow(sim_dat)*ncol(sim_dat)
    cat(paste0("Num simulations: ", num_simulations, "\n"))
    cat(paste0("Mean SNV rate (SNVs/gene): ", sum(as.numeric(names(sim_dat_results))*sim_dat_results)/ncells, "\n"))
    cat(paste0("Num regions: ", nrow(regions), "\n"))
    sapply_out <- sapply(seq(1:length(sim_dat_results)), function(i) { cat(paste0("E(regions with ", names(sim_dat_results)[i], " SNVs): ", sim_dat_results[i]/num_simulations, "\n")) })
    return(sim_dat)
}
#sim_dat_CDH_TSS <- simulate_mutations(HE_TSS_20k_regions, num_individuals=SAMPLE_COUNT)
#mean(sapply(seq(1:nrow(sim_dat_CDH_TSS)), function(i) { return(sum(sim_dat_CDH_TSS[i,]>1)) } ))

num_simulations=10000
#HE_TSS_20k_regions$exp_mutations <- HE_TSS_20k_regions$exp_mutations / SAMPLE_COUNT
#HE_TES_20k_regions$exp_mutations <- HE_TES_20k_regions$exp_mutations / SAMPLE_COUNT
#HE_TSS_20k_regions$exp_mutations <- HE_TSS_20k_regions$exp_mutations / 350
#HE_TES_20k_regions$exp_mutations <- HE_TES_20k_regions$exp_mutations / 350
cat("Simulating mutations for CDH HE TSS 20k regions:\n")
sim_dat_CDH_TSS <- simulate_mutations(HE_TSS_20k_regions, num_individuals=SAMPLE_COUNT)
cat("Simulating mutations for CDH HE TES 20k regions:\n")
sim_dat_CDH_TES <- simulate_mutations(HE_TES_20k_regions, num_individuals=SAMPLE_COUNT)
#cat("Simulating mutations for SSC HE TSS 20k regions:\n")
#sim_dat_SSC_TSS <- simulate_mutations(HE_TSS_20k_regions, num_individuals=438)
#cat("Simulating mutations for SSC HE TES 20k regions:\n")
#sim_dat_SSC_TES <- simulate_mutations(HE_TES_20k_regions, num_individuals=438)



############## PLOTS!!! ##############

# The mut_type parameter can be "snv" or "mut"
plot_recurrance <- function(cases_dat, mut_type="snv", main="test", cases_sim_dat=NULL, cases_permutations, n_greater_than=1, controls_distrib=NULL, controls_bootstrap=NULL, controls_jackknife=NULL, controls_many_jackknife=NULL, controls_block_jackknife=NULL, filename=NULL) {
    lines <- new.env()
    cols <- list()
    ltys <- list()
    mtext_label <- ""
    mtext_label <- paste0(mtext_label, SAMPLE_COUNT, " case samples")
    
    
    # Cases line (must be plotted)
    if (mut_type == "snv") {
        cases_line <- sum(cases_dat$cases_snv_count > 1)
    } else if (mut_type == "mut") { 
        cases_line <- sum(cases_dat$cases_snv_count + cases_dat$cases_indel_count > 1) 
    } else { return(0) }
    label <- paste0(DISEASE, " data (value=", cases_line, ")")
    lines[[label]] <- paste0("v=", cases_line)
    cols[[label]] <- "red"
    
    # Subsampling (simple)
    if(!is.null(controls_distrib)) {
        empirical_p_val = sum(controls_distrib > cases_line)/length(controls_distrib)
        if (empirical_p_val == 0) { empirical_p_val <- paste0("< 1/", length(controls_distrib)) } else { empirical_p_val <- paste0("= " , round(empirical_p_val, 3)) }
        label <- paste0("SSC subsamples (N=", num_subsample_iterations, ")")
        lines[[label]] <- density(controls_distrib)
    }
    
    # Poisson
    if(!is.null(controls_distrib)) {
        poisson_mean = mean(controls_distrib)
        poisson_p_val = round(ppois(cases_line, lambda=poisson_mean, lower=FALSE), 3)
        label <- paste0("SSC Poisson (lambda=", round(poisson_mean,2), ")")
        lines[[label]] <- density(rpois(10000, poisson_mean))
        cols[[label]] <- "purple"; ltys[[label]] <- 2
    }
    
    # Simulation
    if(!is.null(cases_sim_dat)) {
        label <- paste0("Simulation (N=", num_simulations, ")")
        lines[[label]] <- density(sapply(seq(1:nrow(cases_sim_dat)), function(i) { return(sum(cases_sim_dat[i,]>n_greater_than)) } ))
        cols[[label]] <- "orange"; #ltys[[label]] <- 2
    }
    
    # Permutation
    if(!is.null(cases_permutations)) {
        label <- paste0("Permutations (N=", 10000, ")")
        lines[[label]] <- density(cases_permutations)
        cols[[label]] <- "blue"; "green"
        perm_p_val = sum(cases_permutations > cases_line)/length(cases_permutations)
        if (perm_p_val == 0) { perm_p_val <- paste0("< 1/", length(cases_permutations)) } else { perm_p_val <- round(perm_p_val, 5) }
        mtext_label <- paste0(mtext_label, ", Perm. p-value: ", perm_p_val)
    }
    
    # Bootstrap
    if(!is.null(controls_bootstrap)) {
        label <- paste0("SSC bootstrap (N=", num_subsample_iterations, ")")
        lines[[label]] <- density(controls_bootstrap)
        cols[[label]] <- "green"
    }
    
    # Many Jackknife
    if(!is.null(controls_many_jackknife)) {
        jack_estimates <- controls_many_jackknife[["estimates"]]
        jack_vars <- controls_many_jackknife[["vars"]]
        jack_sims <- c()
        for(j in 1:length(jack_estimates)) {
            jack_sims <- c(jack_sims, rnorm(10000, mean=jack_estimates[j], sd=sqrt(jack_vars[j]))) 
        }
        #many_jack_p_val = round(pnorm(cases_line, mean=controls_jackknife[["estimate"]], sd=sqrt(controls_jackknife[["var"]]), lower=FALSE), 3)
        label <- paste0("SSC Many Jackknife (N=10)")
        lines[[label]] <- density(jack_sims)
        cols[[label]] <- "purple"
        
        controls_jackknife <- new.env()
        controls_jackknife[["estimate"]] <- jack_estimates[1]
        controls_jackknife[["var"]] <- jack_vars[1]
    }
    
    # Jackknife
    if(!is.null(controls_jackknife)) {
        jack_p_val = round(pnorm(cases_line, mean=controls_jackknife[["estimate"]], sd=sqrt(controls_jackknife[["var"]]), lower=FALSE), 3)
        label <- paste0("SSC Jackknife (Var=", round(controls_jackknife[["var"]],2), ")")
        lines[[label]] <- density(rnorm(10000, mean=controls_jackknife[["estimate"]], sd=sqrt(controls_jackknife[["var"]])))
        cols[[label]] <- "cyan"
        
        mtext_label <- paste0(mtext_label, ", Jack p-value: ", jack_p_val)
    }
    
    # Block Jackknife
    if(!is.null(controls_block_jackknife)) {
        block_jack_p_val = round(pnorm(cases_line, mean=controls_block_jackknife[["estimate"]], sd=sqrt(controls_block_jackknife[["var"]]), lower=FALSE), 3)
        label <- paste0("SSC Jackknife (Var=", round(controls_block_jackknife[["var"]],2), ")")
        lines[[label]] <- density(rnorm(10000, mean=controls_block_jackknife[["estimate"]], sd=sqrt(controls_block_jackknife[["var"]])))
        cols[[label]] <- "pink"
    }
    
    multi_plot(lines=lines, cols=cols, ltys=ltys, lwd=2, cex.lab=1.8, cex.axis=2, cex.main=1.5, mtext_cex=1.5, mar=1.1*par("mar"), main=main, mtext=mtext_label, xlab=paste0("n regions mutated in >", n_greater_than, " sample"), ylab="Density", legend_location="topleft", legend_cex=1.3, file=filename)
}
plot_recurrance(cases_dat=CDH_TSS_genes_dat, mut_type="snv", main=paste0("CDH TSS Region SNV Recurrence Analysis"), n_greater_than=1, cases_permutations=CDH_permutations_TSS_snv, filename=paste0("CDH3_TSS_over_1_nc_snv.pdf"))
plot_recurrance(cases_dat=CDH_TES_genes_dat, mut_type="snv", main=paste0("CDH 3'UTR Region SNV Recurrence Analysis"), n_greater_than=1, cases_permutations=CDH_permutations_TES_snv, filename=paste0("CDH3_TES_over_1_nc_snv.pdf"))
plot_recurrance(cases_dat=CDH_TSS_genes_dat, mut_type="mut", main=paste0("CDH TSS Region SNV+Indel Recurrence Analysis"), n_greater_than=1, cases_permutations=CDH_permutations_TSS_mut, filename=paste0("CDH3_TSS_over_1_nc_mut.pdf"))
plot_recurrance(cases_dat=CDH_TES_genes_dat, mut_type="mut", main=paste0("CDH 3'UTR Region SNV+Indel Recurrence Analysis"), n_greater_than=1, cases_permutations=CDH_permutations_TES_mut, filename=paste0("CDH3_TES_over_1_nc_mut.pdf"))


main_title_prefix = "CDH HE"
if (DISEASE == "CHD") { main_title_prefix = "CHD HE" }

# CDH HE TES, with >1 SNVs
plot_recurrance(cases_dat=TES_genes_dat, mut_type="snv", cases_sim_dat=sim_dat_CDH_TES, n_greater_than=1, controls_distrib=over_1_snv_SSC_subsamples_TES, controls_bootstrap=over_1_snv_SSC_bootstrap_TES, controls_jackknife=over_1_snv_SSC_jackknife_TES, controls_block_jackknife=over_2_snv_SSC_block_jackknife_TES)
# CDH HE TSS, with >1 SNVs
plot_recurrance(cases_dat=TSS_genes_dat, mut_type="snv", cases_sim_dat=sim_dat_CDH_TSS, n_greater_than=1, controls_distrib=over_1_snv_SSC_subsamples_TSS, controls_bootstrap=over_1_snv_SSC_bootstrap_TSS, controls_jackknife=over_1_snv_SSC_jackknife_TSS, controls_block_jackknife=over_2_snv_SSC_block_jackknife_TSS)
# CDH HE TES, with >2 SNVs
plot_recurrance(cases_dat=TES_genes_dat, mut_type="snv", cases_sim_dat=sim_dat_CDH_TES, n_greater_than=2, controls_distrib=over_2_snv_SSC_subsamples_TES, controls_bootstrap=over_2_snv_SSC_bootstrap_TES, controls_jackknife=over_2_snv_SSC_jackknife_TES, controls_block_jackknife=over_2_snv_SSC_block_jackknife_TES)
# CDH HE TSS, with >2 SNVs
plot_recurrance(cases_dat=TSS_genes_dat, mut_type="snv", cases_sim_dat=sim_dat_CDH_TSS, n_greater_than=2, controls_distrib=over_2_snv_SSC_subsamples_TSS, controls_bootstrap=over_2_snv_SSC_bootstrap_TSS, controls_jackknife=over_2_snv_SSC_jackknife_TSS, controls_block_jackknife=over_2_snv_SSC_block_jackknife_TSS)
# CDH HE TES, with >1 SNVs+Indels
plot_recurrance(cases_dat=TES_genes_dat, mut_type="mut", cases_sim_dat=sim_dat_CDH_TES, n_greater_than=1, controls_distrib=over_1_mut_SSC_subsamples_TES, controls_bootstrap=over_1_mut_SSC_bootstrap_TES, controls_jackknife=over_1_mut_SSC_jackknife_TES, controls_block_jackknife=over_2_mut_SSC_block_jackknife_TES)
# CDH HE TSS, with >1 SNVs+Indels
plot_recurrance(cases_dat=TSS_genes_dat, mut_type="mut", cases_sim_dat=sim_dat_CDH_TSS, n_greater_than=1, controls_distrib=over_1_mut_SSC_subsamples_TSS, controls_bootstrap=over_1_mut_SSC_bootstrap_TSS, controls_jackknife=over_1_mut_SSC_jackknife_TSS, controls_block_jackknife=over_2_mut_SSC_block_jackknife_TSS)
# CDH HE TES, with >2 SNVs+Indels
plot_recurrance(cases_dat=TES_genes_dat, mut_type="mut", cases_sim_dat=sim_dat_CDH_TES, n_greater_than=2, controls_distrib=over_2_mut_SSC_subsamples_TES, controls_bootstrap=over_2_mut_SSC_bootstrap_TES, controls_jackknife=over_2_mut_SSC_jackknife_TES, controls_block_jackknife=over_2_mut_SSC_block_jackknife_TES)
# CDH HE TSS, with >2 SNVs+Indels
plot_recurrance(cases_dat=TSS_genes_dat, mut_type="mut", cases_sim_dat=sim_dat_CDH_TSS, n_greater_than=2, controls_distrib=over_2_mut_SSC_subsamples_TSS, controls_bootstrap=over_2_mut_SSC_bootstrap_TSS, controls_jackknife=over_2_mut_SSC_jackknife_TSS, controls_block_jackknife=over_2_mut_SSC_block_jackknife_TSS)

# CHD HE TES, with >1 SNVs
plot_recurrance(cases_dat=TES_genes_dat, mut_type="snv", cases_sim_dat=sim_dat_CDH_TES, n_greater_than=1, controls_distrib=over_1_snv_SSC_subsamples_TES, controls_bootstrap=over_1_snv_SSC_bootstrap_TES, controls_jackknife=over_1_snv_SSC_jackknife_TES, controls_block_jackknife=over_2_snv_SSC_block_jackknife_TES)
# CHD HE TSS, with >1 SNVs
plot_recurrance(cases_dat=TSS_genes_dat, mut_type="snv", cases_sim_dat=sim_dat_CDH_TSS, n_greater_than=1, controls_distrib=over_1_snv_SSC_subsamples_TSS, controls_bootstrap=over_1_snv_SSC_bootstrap_TSS, controls_jackknife=over_1_snv_SSC_jackknife_TSS, controls_block_jackknife=over_2_snv_SSC_block_jackknife_TSS)
# CHD HE TES, with >2 SNVs
plot_recurrance(cases_dat=TES_genes_dat, mut_type="snv", cases_sim_dat=sim_dat_CDH_TES, n_greater_than=2, controls_distrib=over_2_snv_SSC_subsamples_TES, controls_bootstrap=over_2_snv_SSC_bootstrap_TES, controls_jackknife=over_2_snv_SSC_jackknife_TES, controls_block_jackknife=over_2_snv_SSC_block_jackknife_TES)
# CHD HE TSS, with >2 SNVs
plot_recurrance(cases_dat=TSS_genes_dat, mut_type="snv", cases_sim_dat=sim_dat_CDH_TSS, n_greater_than=2, controls_distrib=over_2_snv_SSC_subsamples_TSS, controls_bootstrap=over_2_snv_SSC_bootstrap_TSS, controls_jackknife=over_2_snv_SSC_jackknife_TSS, controls_block_jackknife=over_2_snv_SSC_block_jackknife_TSS)
# CHD HE TES, with >1 SNVs+Indels
plot_recurrance(cases_dat=TES_genes_dat, mut_type="mut", cases_sim_dat=sim_dat_CDH_TES, n_greater_than=1, controls_distrib=over_1_mut_SSC_subsamples_TES, controls_bootstrap=over_1_mut_SSC_bootstrap_TES, controls_jackknife=over_1_mut_SSC_jackknife_TES, controls_block_jackknife=over_2_mut_SSC_block_jackknife_TES)
# CHD HE TSS, with >1 SNVs+Indels
plot_recurrance(cases_dat=TSS_genes_dat, mut_type="mut", cases_sim_dat=sim_dat_CDH_TSS, n_greater_than=1, controls_distrib=over_1_mut_SSC_subsamples_TSS, controls_bootstrap=over_1_mut_SSC_bootstrap_TSS, controls_jackknife=over_1_mut_SSC_jackknife_TSS, controls_block_jackknife=over_2_mut_SSC_block_jackknife_TSS)
# CHD HE TES, with >2 SNVs+Indels
plot_recurrance(cases_dat=TES_genes_dat, mut_type="mut", cases_sim_dat=sim_dat_CDH_TES, n_greater_than=2, controls_distrib=over_2_mut_SSC_subsamples_TES, controls_bootstrap=over_2_mut_SSC_bootstrap_TES, controls_jackknife=over_2_mut_SSC_jackknife_TES, controls_block_jackknife=over_2_mut_SSC_block_jackknife_TES)
# CHD HE TSS, with >2 SNVs+Indels
plot_recurrance(cases_dat=TSS_genes_dat, mut_type="mut", cases_sim_dat=sim_dat_CDH_TSS, n_greater_than=2, controls_distrib=over_2_mut_SSC_subsamples_TSS, controls_bootstrap=over_2_mut_SSC_bootstrap_TSS, controls_jackknife=over_2_mut_SSC_jackknife_TSS, controls_block_jackknife=over_2_mut_SSC_block_jackknife_TSS)






































SSC_distrib <- over_1_snv_SSC_subsamples_TES
SSC_bootstrap <- over_1_snv_SSC_bootstrap_TES
SSC_jackknife <- over_1_snv_SSC_jackknife_TES
SSC_block_jackknife <- over_1_snv_SSC_block_jackknife_TES
cases_line <- sum(TES_genes_dat$CDH_snv_count > 1)
empirical_p_val = sum(SSC_distrib >= cases_line)/length(SSC_distrib)
if (empirical_p_val == 0) { empirical_p_val <- paste0("< 1/", length(SSC_distrib)) } else { empirical_p_val <- paste0("= " , round(empirical_p_val, 3)) }
plot(density(SSC_distrib), main=paste0(main_title_prefix, " TES Non-Coding SNV Results"), xlab="n genes with > 1 de novo", col="blue", xlim=c(0, max(c(SSC_distrib, SSC_bootstrap))*1.2)) # col=4
lines(density(SSC_bootstrap), col=3)
abline(v=cases_line, col=2)
CDH_sim_distrib <- sapply(seq(1:nrow(sim_dat_CDH_TES)), function(i) { return(sum(sim_dat_CDH_TES[i,]>1)) } )
lines(density(CDH_sim_distrib), col="orange", lty=2)
poisson_mean = mean(SSC_distrib)
poisson_p_val = round(ppois(cases_line, lambda=poisson_mean, lower=FALSE), 3)
jack_p_val = round(pnorm(cases_line, mean=SSC_jackknife[["estimate"]], sd=sqrt(SSC_jackknife[["var"]]), lower=FALSE), 3)
#lines(density(rpois(10000, poisson_mean)), col="purple", lty=2)
lines(density(rnorm(10000, mean=SSC_jackknife[["estimate"]], sd=sqrt(SSC_jackknife[["var"]]))), col="cyan")
lines(density(rnorm(10000, mean=SSC_block_jackknife[["estimate"]], sd=sqrt(SSC_block_jackknife[["var"]]))), col="pink")
poisson_legend_label <- paste0("SSC Poisson (lambda=", round(poisson_mean,2), ")")
legend("topleft", legend=c(paste0("SSC subsamples (N=", num_subsample_iterations, ")"), paste0("SSC Jackknife (Var=", round(SSC_jackknife[["var"]],2), ")"), paste0("SSC bootstrap (N=", num_subsample_iterations, ")"), poisson_legend_label, paste0(DISEASE, " data (value=", cases_line, ")"), paste0("Simulation (N=", num_simulations, ")")), col=c(4,"cyan",3,"purple",2,"orange"), lty=c(1,1,1,2,1,2), lwd = 1, cex = 0.9)
mtext(paste0(SAMPLE_COUNT, " samples, Emp. p-value ", empirical_p_val, ", Poiss. p-value: ", poisson_p_val, ", Jack p-value: ", jack_p_val))
dev.copy2pdf(file=paste0(DISEASE, "_TES_over_1_nc_snv.pdf"))

# CHD HE TSS, with >1 SNVs
SSC_distrib <- over_1_snv_SSC_subsamples_TSS
SSC_bootstrap <- over_1_snv_SSC_bootstrap_TSS
SSC_jackknife <- over_1_snv_SSC_jackknife_TSS
SSC_block_jackknife <- over_1_snv_SSC_block_jackknife_TSS
cases_line <- sum(TSS_genes_dat$CDH_snv_count > 1)
empirical_p_val = sum(SSC_distrib >= cases_line)/length(SSC_distrib)
if (empirical_p_val == 0) { empirical_p_val <- paste0("< 1/", length(SSC_distrib)) } else { empirical_p_val <- paste0("= " , round(empirical_p_val, 3)) }
plot(density(SSC_distrib), main=paste0(main_title_prefix, " TSS Non-Coding SNV Results"), xlab="n genes with > 1 de novo", col="blue", xlim=c(0, max(c(SSC_distrib, SSC_bootstrap))*1.2))
lines(density(SSC_bootstrap), col=3)
abline(v=cases_line, col=2)
CDH_sim_distrib <- sapply(seq(1:nrow(sim_dat_CDH_TSS)), function(i) { return(sum(sim_dat_CDH_TSS[i,]>1)) } )
lines(density(CDH_sim_distrib), col="orange", lty=2)
poisson_mean = mean(SSC_distrib)
poisson_p_val = round(ppois(cases_line, lambda=poisson_mean, lower=FALSE), 3)
jack_p_val = round(pnorm(cases_line, mean=SSC_jackknife[["estimate"]], sd=sqrt(SSC_jackknife[["var"]]), lower=FALSE), 3)
#lines(density(rpois(10000, poisson_mean)), col="purple", lty=2)
lines(density(rnorm(10000, mean=SSC_jackknife[["estimate"]], sd=sqrt(SSC_jackknife[["var"]]))), col="cyan")
lines(density(rnorm(10000, mean=SSC_block_jackknife[["estimate"]], sd=sqrt(SSC_block_jackknife[["var"]]))), col="pink")
poisson_legend_label <- paste0("SSC Poisson (lambda=", round(poisson_mean,2), ")")
legend("topleft", legend=c(paste0("SSC subsamples (N=", num_subsample_iterations, ")"), paste0("SSC Jackknife (Var=", round(SSC_jackknife[["var"]],2), ")"), paste0("SSC bootstrap (N=", num_subsample_iterations, ")"), poisson_legend_label, paste0(DISEASE, " data (value=", cases_line, ")"), paste0("Simulation (N=", num_simulations, ")")), col=c(4,"cyan",3,"purple",2,"orange"), lty=c(1,1,1,2,1,2), lwd = 1, cex = 0.9)
mtext(paste0(SAMPLE_COUNT, " samples, Emp. p-value ", empirical_p_val, ", Poiss. p-value: ", poisson_p_val, ", Jack p-value: ", jack_p_val))
dev.copy2pdf(file=paste0(DISEASE, "_TSS_over_1_nc_snv.pdf"))

# CHD HE TES, with >2 SNVs
SSC_distrib <- over_2_snv_SSC_subsamples_TES
SSC_bootstrap <- over_2_snv_SSC_bootstrap_TES
SSC_jackknife <- over_2_snv_SSC_jackknife_TES
SSC_block_jackknife <- over_2_snv_SSC_block_jackknife_TES
cases_line <- sum(TES_genes_dat$CDH_snv_count > 2)
empirical_p_val = sum(SSC_distrib >= cases_line)/length(SSC_distrib)
if (empirical_p_val == 0) { empirical_p_val <- paste0("< 1/", length(SSC_distrib)) } else { empirical_p_val <- paste0("= " , round(empirical_p_val, 3)) }
plot(density(SSC_distrib), main=paste0(main_title_prefix, " TES Non-Coding SNV Results"), xlab="n genes with > 2 de novos", col="blue", xlim=c(0, max(c(SSC_distrib, SSC_bootstrap))*1.2))
lines(density(SSC_bootstrap), col=3)
abline(v=cases_line, col=2)
CDH_sim_distrib <- sapply(seq(1:nrow(sim_dat_CDH_TES)), function(i) { return(sum(sim_dat_CDH_TES[i,]>2)) } )
lines(density(CDH_sim_distrib), col="orange", lty=2)
poisson_mean = mean(SSC_distrib)
poisson_p_val = round(ppois(cases_line, lambda=poisson_mean, lower=FALSE), 3)
jack_p_val = round(pnorm(cases_line, mean=SSC_jackknife[["estimate"]], sd=sqrt(SSC_jackknife[["var"]]), lower=FALSE), 3)
#lines(density(rpois(10000, poisson_mean)), col="purple", lty=2)
lines(density(rnorm(10000, mean=SSC_jackknife[["estimate"]], sd=sqrt(SSC_jackknife[["var"]]))), col="cyan")
lines(density(rnorm(10000, mean=SSC_block_jackknife[["estimate"]], sd=sqrt(SSC_block_jackknife[["var"]]))), col="pink")
poisson_legend_label <- paste0("SSC Poisson (lambda=", round(poisson_mean,2), ")")
legend("topleft", legend=c(paste0("SSC subsamples (N=", num_subsample_iterations, ")"), paste0("SSC Jackknife (Var=", round(SSC_jackknife[["var"]],2), ")"), paste0("SSC bootstrap (N=", num_subsample_iterations, ")"), poisson_legend_label, paste0(DISEASE, " data (value=", cases_line, ")"), paste0("Simulation (N=", num_simulations, ")")), col=c(4,"cyan",3,"purple",2,"orange"), lty=c(1,1,1,2,1,2), lwd = 1, cex = 0.9)
mtext(paste0(SAMPLE_COUNT, " samples, Emp. p-value ", empirical_p_val, ", Poiss. p-value: ", poisson_p_val, ", Jack p-value: ", jack_p_val))
dev.copy2pdf(file=paste0(DISEASE, "_TES_over_2_nc_snv.pdf"))

# CHD HE TSS, with >2 SNVs
SSC_distrib <- over_2_snv_SSC_subsamples_TSS
SSC_bootstrap <- over_2_snv_SSC_bootstrap_TSS
SSC_jackknife <- over_2_snv_SSC_jackknife_TSS
SSC_block_jackknife <- over_2_snv_SSC_block_jackknife_TSS
cases_line <- sum(TSS_genes_dat$CDH_snv_count > 2)
empirical_p_val = sum(SSC_distrib >= cases_line)/length(SSC_distrib)
if (empirical_p_val == 0) { empirical_p_val <- paste0("< 1/", length(SSC_distrib)) } else { empirical_p_val <- paste0("= " , round(empirical_p_val, 3)) }
plot(density(SSC_distrib), main=paste0(main_title_prefix, " TSS Non-Coding SNV Results"), xlab="n genes with > 2 de novos", col="blue", xlim=c(0, max(c(SSC_distrib, SSC_bootstrap))*1.2))
lines(density(SSC_bootstrap), col=3)
abline(v=cases_line, col=2)
CDH_sim_distrib <- sapply(seq(1:nrow(sim_dat_CDH_TSS)), function(i) { return(sum(sim_dat_CDH_TSS[i,]>2)) } )
lines(density(CDH_sim_distrib), col="orange", lty=2)
poisson_mean = mean(SSC_distrib)
poisson_p_val = round(ppois(cases_line, lambda=poisson_mean, lower=FALSE), 3)
jack_p_val = round(pnorm(cases_line, mean=SSC_jackknife[["estimate"]], sd=sqrt(SSC_jackknife[["var"]]), lower=FALSE), 3)
#lines(density(rpois(10000, poisson_mean)), col="purple", lty=2)
lines(density(rnorm(10000, mean=SSC_jackknife[["estimate"]], sd=sqrt(SSC_jackknife[["var"]]))), col="cyan")
lines(density(rnorm(10000, mean=SSC_block_jackknife[["estimate"]], sd=sqrt(SSC_block_jackknife[["var"]]))), col="pink")
poisson_legend_label <- paste0("SSC Poisson (lambda=", round(poisson_mean,2), ")")
legend("topleft", legend=c(paste0("SSC subsamples (N=", num_subsample_iterations, ")"), paste0("SSC Jackknife (Var=", round(SSC_jackknife[["var"]],2), ")"), paste0("SSC bootstrap (N=", num_subsample_iterations, ")"), poisson_legend_label, paste0(DISEASE, " data (value=", cases_line, ")"), paste0("Simulation (N=", num_simulations, ")")), col=c(4,"cyan",3,"purple",2,"orange"), lty=c(1,1,1,2,1,2), lwd = 1, cex = 0.9)
mtext(paste0(SAMPLE_COUNT, " samples, Emp. p-value ", empirical_p_val, ", Poiss. p-value: ", poisson_p_val, ", Jack p-value: ", jack_p_val))
dev.copy2pdf(file=paste0(DISEASE, "_TSS_over_2_nc_snv.pdf"))

############## PLOTS FOR SNVs+INDELs ##############

# CHD HE TES, with >1 SNVs+Indels
SSC_distrib <- over_1_mut_SSC_subsamples_TES
SSC_bootstrap <- over_1_mut_SSC_bootstrap_TES
SSC_jackknife <- over_1_mut_SSC_jackknife_TES
SSC_block_jackknife <- over_1_mut_SSC_block_jackknife_TES
cases_line <- sum(TES_genes_dat$CDH_snv_count + TES_genes_dat$CDH_indel_count > 1)
empirical_p_val = sum(SSC_distrib >= cases_line)/length(SSC_distrib)
if (empirical_p_val == 0) { empirical_p_val <- paste0("< 1/", length(SSC_distrib)) } else { empirical_p_val <- paste0("= " , round(empirical_p_val, 3)) }
plot(density(SSC_distrib), main=paste0(main_title_prefix, " TES Non-Coding SNV+Indel Results"), xlab="n genes with > 1 de novo", col="blue", xlim=c(0, max(c(SSC_distrib, SSC_bootstrap))*1.2))
lines(density(SSC_bootstrap), col=3)
abline(v=cases_line, col=2)
poisson_mean = mean(SSC_distrib)
poisson_p_val = round(ppois(cases_line, lambda=poisson_mean, lower=FALSE), 3)
jack_p_val = round(pnorm(cases_line, mean=SSC_jackknife[["estimate"]], sd=sqrt(SSC_jackknife[["var"]]), lower=FALSE), 3)
#lines(density(rpois(10000, poisson_mean)), col="purple", lty=2)
lines(density(rnorm(10000, mean=SSC_jackknife[["estimate"]], sd=sqrt(SSC_jackknife[["var"]]))), col="cyan")
lines(density(rnorm(10000, mean=SSC_block_jackknife[["estimate"]], sd=sqrt(SSC_block_jackknife[["var"]]))), col="pink")
#abline(v=over_1_mut_SSC_jackknife_TES[["estimate"]], col="cyan")
poisson_legend_label <- paste0("SSC Poisson (lambda=", round(poisson_mean,2), ")")
legend("topleft", legend=c(paste0("SSC subsamples (N=", num_subsample_iterations, ")"), paste0("SSC Jackknife (Var=", round(SSC_jackknife[["var"]],2), ")"), paste0("SSC bootstrap (N=", num_subsample_iterations, ")"), poisson_legend_label, paste0(DISEASE, " data (value=", cases_line, ")")), col=c(4,"cyan",3,"purple",2), lty=c(1,1,1,2,1), lwd = 1, cex = 0.9)
mtext(paste0(SAMPLE_COUNT, " samples, Emp. p-value ", empirical_p_val, ", Poiss. p-value: ", poisson_p_val, ", Jack p-value: ", jack_p_val))
dev.copy2pdf(file=paste0(DISEASE, "_TES_over_1_nc_mut.pdf"))

# CHD HE TSS, with >1 SNVs+Indels
SSC_distrib <- over_1_mut_SSC_subsamples_TSS
SSC_bootstrap <- over_1_mut_SSC_bootstrap_TSS
SSC_jackknife <- over_1_mut_SSC_jackknife_TSS
SSC_block_jackknife <- over_1_mut_SSC_block_jackknife_TSS
cases_line <- sum(TSS_genes_dat$CDH_snv_count + TSS_genes_dat$CDH_indel_count > 1)
empirical_p_val = sum(SSC_distrib >= cases_line)/length(SSC_distrib)
if (empirical_p_val == 0) { empirical_p_val <- paste0("< 1/", length(SSC_distrib)) } else { empirical_p_val <- paste0("= " , round(empirical_p_val, 3)) }
plot(density(SSC_distrib), main=paste0(main_title_prefix, " TSS Non-Coding SNV+Indel Results"), xlab="n genes with > 1 de novo", col="blue", xlim=c(0, max(c(SSC_distrib, SSC_bootstrap))*1.2))
lines(density(SSC_bootstrap), col=3)
abline(v=cases_line, col=2)
poisson_mean = mean(SSC_distrib)
poisson_p_val = round(ppois(cases_line, lambda=poisson_mean, lower=FALSE), 3)
jack_p_val = round(pnorm(cases_line, mean=SSC_jackknife[["estimate"]], sd=sqrt(SSC_jackknife[["var"]]), lower=FALSE), 3)
#lines(density(rpois(10000, poisson_mean)), col="purple", lty=2)
lines(density(rnorm(10000, mean=SSC_jackknife[["estimate"]], sd=sqrt(SSC_jackknife[["var"]]))), col="cyan")
lines(density(rnorm(10000, mean=SSC_block_jackknife[["estimate"]], sd=sqrt(SSC_block_jackknife[["var"]]))), col="pink")
#abline(v=over_1_mut_SSC_jackknife_TSS[["estimate"]], col="cyan")
poisson_legend_label <- paste0("SSC Poisson (lambda=", round(poisson_mean,2), ")")
legend("topleft", legend=c(paste0("SSC subsamples (N=", num_subsample_iterations, ")"), paste0("SSC Jackknife (Var=", round(SSC_jackknife[["var"]],2), ")"), paste0("SSC bootstrap (N=", num_subsample_iterations, ")"), poisson_legend_label, paste0(DISEASE, " data (value=", cases_line, ")")), col=c(4,"cyan",3,"purple",2), lty=c(1,1,1,2,1), lwd = 1, cex = 0.9)
mtext(paste0(SAMPLE_COUNT, " samples, Emp. p-value ", empirical_p_val, ", Poiss. p-value: ", poisson_p_val, ", Jack p-value: ", jack_p_val))
dev.copy2pdf(file=paste0(DISEASE, "_TSS_over_1_nc_mut.pdf"))

# CHD HE TES, with >2 SNVs+Indels
SSC_distrib <- over_2_mut_SSC_subsamples_TES
SSC_bootstrap <- over_2_mut_SSC_bootstrap_TES
SSC_jackknife <- over_2_mut_SSC_jackknife_TES
SSC_block_jackknife <- over_2_mut_SSC_block_jackknife_TES
cases_line <- sum(TES_genes_dat$CDH_snv_count + TES_genes_dat$CDH_indel_count > 2)
empirical_p_val = sum(SSC_distrib >= cases_line)/length(SSC_distrib)
if (empirical_p_val == 0) { empirical_p_val <- paste0("< 1/", length(SSC_distrib)) } else { empirical_p_val <- paste0("= " , round(empirical_p_val, 3)) }
plot(density(SSC_distrib), main=paste0(main_title_prefix, " TES Non-Coding SNV+Indel Results"), xlab="n genes with > 2 de novos", col="blue", xlim=c(0, max(c(SSC_distrib, SSC_bootstrap))*1.2))
lines(density(SSC_bootstrap), col=3)
abline(v=cases_line, col=2)
poisson_mean = mean(SSC_distrib)
poisson_p_val = round(ppois(cases_line, lambda=poisson_mean, lower=FALSE), 3)
jack_p_val = round(pnorm(cases_line, mean=SSC_jackknife[["estimate"]], sd=sqrt(SSC_jackknife[["var"]]), lower=FALSE), 3)
#lines(density(rpois(10000, poisson_mean)), col="purple", lty=2)
lines(density(rnorm(10000, mean=SSC_jackknife[["estimate"]], sd=sqrt(SSC_jackknife[["var"]]))), col="cyan")
lines(density(rnorm(10000, mean=SSC_block_jackknife[["estimate"]], sd=sqrt(SSC_block_jackknife[["var"]]))), col="pink")
poisson_legend_label <- paste0("SSC Poisson (lambda=", round(poisson_mean,2), ")")
legend("topleft", legend=c(paste0("SSC subsamples (N=", num_subsample_iterations, ")"), paste0("SSC Jackknife (Var=", round(SSC_jackknife[["var"]],2), ")"), paste0("SSC bootstrap (N=", num_subsample_iterations, ")"), poisson_legend_label, paste0(DISEASE, " data (value=", cases_line, ")")), col=c(4,"cyan",3,"purple",2), lty=c(1,1,1,2,1), lwd = 1, cex = 0.9)
mtext(paste0(SAMPLE_COUNT, " samples, Emp. p-value ", empirical_p_val, ", Poiss. p-value: ", poisson_p_val, ", Jack p-value: ", jack_p_val))
dev.copy2pdf(file=paste0(DISEASE, "_TES_over_2_nc_mut.pdf"))

# CHD HE TSS, with >2 SNVs+Indels
SSC_distrib <- over_2_mut_SSC_subsamples_TSS
SSC_bootstrap <- over_2_mut_SSC_bootstrap_TSS
SSC_jackknife <- over_2_mut_SSC_jackknife_TSS
SSC_block_jackknife <- over_2_mut_SSC_block_jackknife_TSS
cases_line <- sum(TSS_genes_dat$CDH_snv_count + TSS_genes_dat$CDH_indel_count > 2)
empirical_p_val = sum(SSC_distrib >= cases_line)/length(SSC_distrib)
if (empirical_p_val == 0) { empirical_p_val <- paste0("< 1/", length(SSC_distrib)) } else { empirical_p_val <- paste0("= " , round(empirical_p_val, 3)) }
plot(density(SSC_distrib), main=paste0(main_title_prefix, " TSS Non-Coding SNV+Indel Results"), xlab="n genes with > 2 de novos", col="blue", xlim=c(0, max(c(SSC_distrib, SSC_bootstrap))*1.2))
lines(density(SSC_bootstrap), col=3)
abline(v=cases_line, col=2)
poisson_mean = mean(SSC_distrib)
poisson_p_val = round(ppois(cases_line, lambda=poisson_mean, lower=FALSE), 3)
jack_p_val = round(pnorm(cases_line, mean=SSC_jackknife[["estimate"]], sd=sqrt(SSC_jackknife[["var"]]), lower=FALSE), 3)
#lines(density(rpois(10000, poisson_mean)), col="purple", lty=2)
lines(density(rnorm(10000, mean=SSC_jackknife[["estimate"]], sd=sqrt(SSC_jackknife[["var"]]))), col="cyan")
lines(density(rnorm(10000, mean=SSC_block_jackknife[["estimate"]], sd=sqrt(SSC_block_jackknife[["var"]]))), col="pink")
poisson_legend_label <- paste0("SSC Poisson (lambda=", round(poisson_mean,2), ")")
legend("topleft", legend=c(paste0("SSC subsamples (N=", num_subsample_iterations, ")"), paste0("SSC Jackknife (Var=", round(SSC_jackknife[["var"]],2), ")"), paste0("SSC bootstrap (N=", num_subsample_iterations, ")"), poisson_legend_label, paste0(DISEASE, " data (value=", cases_line, ")")), col=c(4,"cyan",3,"purple",2), lty=c(1,1,1,2,1), lwd = 1, cex = 0.9)
mtext(paste0(SAMPLE_COUNT, " samples, Emp. p-value ", empirical_p_val, ", Poiss. p-value: ", poisson_p_val, ", Jack p-value: ", jack_p_val))
dev.copy2pdf(file=paste0(DISEASE, "_TSS_over_2_nc_mut.pdf"))

#over_1_snv_CDH <- sapply(seq(1:nrow(sim_dat_CDH_TES)), function(i) { return(sum(sim_dat_CDH_TES[i,]>1)) } )
##plot(density(over_1_snv_CDH), main="HE TES SNV Simulation", xlab="n genes with > 1 de novo", col=1, xlim=c(0, 80), ylim=c(0, 0.1))
#over_1_snv_CDH <- sapply(seq(1:nrow(sim_dat_CDH_TES)), function(i) { return(sum(sim_dat_CDH_TES[i,]>1)) } )
##plot(density(over_1_snv_CDH), main="HE TES SNV Simulation", xlab="n genes with > 1 de novo", col=1, xlim=c(0, 80), ylim=c(0, 0.1))

#ppois(cases_line, mean(SSC_distrib), lower=FALSE)


# Code from before March 14th, 2017
#over_1_snv_CDH <- sapply(seq(1:nrow(sim_dat_CDH_TES)), function(i) { return(sum(sim_dat_CDH_TES[i,]>1)) } )
##plot(density(over_1_snv_CDH), main="HE TES SNV Simulation", xlab="n genes with > 1 de novo", col=1, xlim=c(0, 80), ylim=c(0, 0.1))
#plot(density(over_1_snv_CDH), main="CHD HE TES SNV Simulation", xlab="n genes with > 1 de novo", col=1, xlim=c(0, 100), ylim=c(0, 0.1))
#mtext(paste0(SAMPLE_COUNT, " samples, ", nrow(HE_TES_20k_regions), " total genes, ", num_simulations, " sims, ", num_subsample_iterations, " SSC subsamples"))
#lines(density(over_1_snv_SSC_subsamples_TES), col=4)
#lines(density(rpois(100000, mean(over_1_snv_SSC_subsamples_TES))), col=3)
#abline(v=sum(TES_genes_dat$CDH_snv_count > 1), col=2)
#legend("topleft", legend=c("Simulations", paste0(DISEASE, " data"), "SSC subsamples", paste0("Poisson(lambda=", round(mean(over_1_snv_SSC_subsamples_TES),2), ")")), col=c(1,2,4,3), lty=c(1, 1, 1), lwd = 1, cex = 0.85)
#dev.copy2pdf(file=paste0(DISEASE, "_HE_TES_simulation.pdf"))



print(paste0(nrow(columbia)/length(unique(columbia_fams))))
print(paste0(nrow(ssc_controls)/length(unique(ssc_controls$sample))))

cases_freeze <- read.csv("CHD_WGS_DNMs_PCGC.txt", sep="\t")
controls_freeze <- read.csv("CHD_WGS_DNMs_Simons.txt", sep="\t")
cases_freeze_vars <- paste0(cases_freeze$Chrom, ":", cases_freeze$Position, cases_freeze$Ref, ">", cases_freeze$Alt)
cases_vars <- paste0("chr", columbia$CHROM, ":", columbia$POS, columbia$REF, ">", columbia$ALT)
sum(cases_vars %in% cases_freeze_vars)/length(cases_vars)
controls_freeze_vars <- paste0(controls_freeze$Chrom, ":", controls_freeze$Position, controls_freeze$Ref, ">", controls_freeze$Alt)
controls_vars <- paste0("chr", ssc_controls$X.chr, ":", ssc_controls$pos, ssc_controls$REF, ">", ssc_controls$ALT)
sum(controls_vars %in% controls_freeze_vars)/length(controls_vars)


samocha <- read.csv("Daly_rate.csv")
gene_ranges <- hg19_genebody_HE[,c(1,2,3,5)]
colnames(gene_ranges) <- c("chromosome", "start", "end", "gene")
gene_ranges <- gene_ranges[!duplicated(gene_ranges$gene),]
gene_ranges$chromosome <- gsub("chr", "", paste0(gene_ranges$chromosome))
gene_ranges_expected_muts <- get_expected_mut_counts(gene_ranges)

samocha_sim_comparison <- merge(gene_ranges_expected_muts, samocha, by="gene")
samocha_sim_comparison$all <- exp(samocha_sim_comparison$all)
samocha_sim_comparison$syn <- exp(samocha_sim_comparison$syn)
samocha_sim_comparison$mis <- exp(samocha_sim_comparison$mis)
samocha_sim_comparison$non <- exp(samocha_sim_comparison$non)
samocha_sim_comparison$splice_site <- exp(samocha_sim_comparison$splice_site)
samocha_sim_comparison$frameshift <- exp(samocha_sim_comparison$frameshift)
plot(samocha_sim_comparison$all, samocha_sim_comparison$exp_mutations, main="Simulated Gene Mutation Rate Comparison", xlab="Samocha", ylab="Simulation")
line_of_best_fit <- lm(samocha_sim_comparison$exp_mutations ~ samocha_sim_comparison$all)
abline(line_of_best_fit, col=2)
mtext(paste0("Corr: ", cor(samocha_sim_comparison$all, samocha_sim_comparison$exp_mutations), ", Best Fit Slope: ", line_of_best_fit$coefficients[2]))


# Finds region overlaps.
get_region_overlaps <- function(HE_TS_regions) {
    regions_vector <- c()
    for(i in 1:nrow(HE_TS_regions)) {
        if (length(regions_vector) == 0) {
            regions_vector <- GRanges(seqnames=HE_TS_regions$chromosome[i], ranges=IRanges(start=HE_TS_regions$start[i], end=HE_TS_regions$end[i], names=HE_TS_regions$gene[i]))
        } else {
            regions_vector <- c(regions_vector, GRanges(seqnames=HE_TS_regions$chromosome[i], ranges=IRanges(start=HE_TS_regions$start[i], end=HE_TS_regions$end[i], names=HE_TS_regions$gene[i])))
        }
    }
    hits <- data.frame(olRanges(regions_vector, regions_vector))
    hits_same <- hits[(hits$OLpercS == 100),]
    
    hits <- hits[(hits$OLpercS != 100),]
    cat(paste0("Fraction bases overlapped: ", sum(hits$OLlength), "/", sum(hits_same$OLlength), " = ", round(sum(hits$OLlength)/sum(hits_same$OLlength), 4)), "\n")
    return(hits)
}
#cat(paste0(DISEASE, " TSS overlaps analysis:\n"))
#ols <- get_region_overlaps(HE_TSS_20k_regions)
#cat(paste0(DISEASE, " TES overlaps analysis:\n"))
#ols <- get_region_overlaps(HE_TES_20k_regions)




over_1_snv_SSC_subsamples_TSS <- sapply(seq(1:num_subsample_iterations), function(x) { print(paste0("Subsample iteration: ", x)); sub_snv_counts <- subsample_snv_count(controls, "TSS", with_replacement=FALSE); return(sum(sub_snv_counts$snv_count > 1)) } )
over_1_snv_SSC_subsamples_TES <- sapply(seq(1:num_subsample_iterations), function(x) { print(paste0("Subsample iteration: ", x)); sub_snv_counts <- subsample_snv_count(controls, "TES", with_replacement=FALSE); return(sum(sub_snv_counts$snv_count > 1)) } )
over_2_snv_SSC_subsamples_TSS <- sapply(seq(1:num_subsample_iterations), function(x) { print(paste0("Subsample iteration: ", x)); sub_snv_counts <- subsample_snv_count(controls, "TSS", with_replacement=FALSE); return(sum(sub_snv_counts$snv_count > 2)) } )
over_2_snv_SSC_subsamples_TES <- sapply(seq(1:num_subsample_iterations), function(x) { print(paste0("Subsample iteration: ", x)); sub_snv_counts <- subsample_snv_count(controls, "TES", with_replacement=FALSE); return(sum(sub_snv_counts$snv_count > 2)) } )

ssc_controls <- read.csv("SSC_denovo_IGV_noncoding.csv")
colnames(ssc_controls)[c(1,2,9)] <- c("X.chr", "pos", "sample")
ssc_controls$sample <- sapply(1:nrow(ssc_controls), function(i) { strsplit(paste(ssc_controls$sample[i]), "\\(")[[1]][1] })
over_1_mut_SSC_subsamples_TSS <- sapply(seq(1:num_subsample_iterations), function(x) { print(paste0("Subsample iteration: ", x)); sub_snv_counts <- subsample_snv_count(ssc_controls, "TSS", with_replacement=FALSE); return(sum(sub_snv_counts$snv_count > 1)) } )
over_1_mut_SSC_subsamples_TES <- sapply(seq(1:num_subsample_iterations), function(x) { print(paste0("Subsample iteration: ", x)); sub_snv_counts <- subsample_snv_count(ssc_controls, "TES", with_replacement=FALSE); return(sum(sub_snv_counts$snv_count > 1)) } )
over_2_mut_SSC_subsamples_TSS <- sapply(seq(1:num_subsample_iterations), function(x) { print(paste0("Subsample iteration: ", x)); sub_snv_counts <- subsample_snv_count(ssc_controls, "TSS", with_replacement=FALSE); return(sum(sub_snv_counts$snv_count > 2)) } )
over_2_mut_SSC_subsamples_TES <- sapply(seq(1:num_subsample_iterations), function(x) { print(paste0("Subsample iteration: ", x)); sub_snv_counts <- subsample_snv_count(ssc_controls, "TES", with_replacement=FALSE); return(sum(sub_snv_counts$snv_count > 2)) } )

# Block Jackknife
over_1_snv_SSC_block_jackknife_TSS <- jackknife_snv_count(controls, sites="TSS", sample_size=SAMPLE_COUNT)
over_1_snv_SSC_block_jackknife_TES <- jackknife_snv_count(controls, sites="TES", sample_size=SAMPLE_COUNT)
over_2_snv_SSC_block_jackknife_TSS <- jackknife_snv_count(controls, sites="TSS", sample_size=SAMPLE_COUNT, n_greater_than=2)
over_2_snv_SSC_block_jackknife_TES <- jackknife_snv_count(controls, sites="TES", sample_size=SAMPLE_COUNT, n_greater_than=2)
over_1_mut_SSC_block_jackknife_TSS <- jackknife_snv_count(ssc_controls, sites="TSS", sample_size=SAMPLE_COUNT)
over_1_mut_SSC_block_jackknife_TES <- jackknife_snv_count(ssc_controls, sites="TES", sample_size=SAMPLE_COUNT)
over_2_mut_SSC_block_jackknife_TSS <- jackknife_snv_count(ssc_controls, sites="TSS", sample_size=SAMPLE_COUNT, n_greater_than=2)
over_2_mut_SSC_block_jackknife_TES <- jackknife_snv_count(ssc_controls, sites="TES", sample_size=SAMPLE_COUNT, n_greater_than=2)

# Bootstrap (subsamples with replacement)
num_subsample_iterations = 100
over_1_snv_SSC_bootstrap_TSS <- sapply(seq(1:num_subsample_iterations), function(x) { print(paste0("Subsample iteration: ", x)); sub_snv_counts <- bootstrap_snv_count(controls, "TSS"); return(sum(sub_snv_counts$snv_count > 1)) } )
over_1_snv_SSC_bootstrap_TES <- sapply(seq(1:num_subsample_iterations), function(x) { print(paste0("Subsample iteration: ", x)); sub_snv_counts <- bootstrap_snv_count(controls, "TES"); return(sum(sub_snv_counts$snv_count > 1)) } )
over_2_snv_SSC_bootstrap_TSS <- sapply(seq(1:num_subsample_iterations), function(x) { print(paste0("Subsample iteration: ", x)); sub_snv_counts <- bootstrap_snv_count(controls, "TSS"); return(sum(sub_snv_counts$snv_count > 2)) } )
over_2_snv_SSC_bootstrap_TES <- sapply(seq(1:num_subsample_iterations), function(x) { print(paste0("Subsample iteration: ", x)); sub_snv_counts <- bootstrap_snv_count(controls, "TES"); return(sum(sub_snv_counts$snv_count > 2)) } )
over_1_mut_SSC_bootstrap_TSS <- sapply(seq(1:num_subsample_iterations), function(x) { print(paste0("Subsample iteration: ", x)); sub_snv_counts <- bootstrap_snv_count(ssc_controls, "TSS"); return(sum(sub_snv_counts$snv_count > 1)) } )
over_1_mut_SSC_bootstrap_TES <- sapply(seq(1:num_subsample_iterations), function(x) { print(paste0("Subsample iteration: ", x)); sub_snv_counts <- bootstrap_snv_count(ssc_controls, "TES"); return(sum(sub_snv_counts$snv_count > 1)) } )
over_2_mut_SSC_bootstrap_TSS <- sapply(seq(1:num_subsample_iterations), function(x) { print(paste0("Subsample iteration: ", x)); sub_snv_counts <- bootstrap_snv_count(ssc_controls, "TSS"); return(sum(sub_snv_counts$snv_count > 2)) } )
over_2_mut_SSC_bootstrap_TES <- sapply(seq(1:num_subsample_iterations), function(x) { print(paste0("Subsample iteration: ", x)); sub_snv_counts <- bootstrap_snv_count(ssc_controls, "TES"); return(sum(sub_snv_counts$snv_count > 2)) } )

# Jackknife (leave-one-out), using random 351 individuals.
controls_350 <- controls[controls$sample %in% sample(unique(paste(controls$sample)), SAMPLE_COUNT+1),]
ssc_controls_350 <- ssc_controls[ssc_controls$sample %in% sample(unique(paste(ssc_controls$sample)), SAMPLE_COUNT+1),]
over_1_snv_SSC_jackknife_TSS <- jackknife_snv_count(controls_350, sites="TSS", block_size=1)
over_1_snv_SSC_jackknife_TES <- jackknife_snv_count(controls_350, sites="TES", block_size=1)
over_2_snv_SSC_jackknife_TSS <- jackknife_snv_count(controls_350, sites="TSS", block_size=1, n_greater_than=2)
over_2_snv_SSC_jackknife_TES <- jackknife_snv_count(controls_350, sites="TES", block_size=1, n_greater_than=2)
over_1_mut_SSC_jackknife_TSS <- jackknife_snv_count(ssc_controls_350, sites="TSS", block_size=1)
over_1_mut_SSC_jackknife_TES <- jackknife_snv_count(ssc_controls_350, sites="TES", block_size=1)
over_2_mut_SSC_jackknife_TSS <- jackknife_snv_count(ssc_controls_350, sites="TSS", block_size=1, n_greater_than=2)
over_2_mut_SSC_jackknife_TES <- jackknife_snv_count(ssc_controls_350, sites="TES", block_size=1, n_greater_than=2)
#print(paste0("TSS SNV Jackknife: ", over_1_snv_SSC_jackknife_TSS[["estimate"]], ", Var = ", over_1_snv_SSC_jackknife_TSS[["var"]]))
#print(paste0("TES SNV Jackknife: ", over_1_snv_SSC_jackknife_TES[["estimate"]], ", Var = ", over_1_snv_SSC_jackknife_TES[["var"]]))
#print(paste0("TSS SNV+Indel Jackknife: ", over_1_mut_SSC_jackknife_TSS[["estimate"]], ", Var = ", over_1_mut_SSC_jackknife_TSS[["var"]]))
#print(paste0("TES SNV+Indel Jackknife: ", over_1_mut_SSC_jackknife_TES[["estimate"]], ", Var = ", over_1_mut_SSC_jackknife_TES[["var"]]))

# Block Jackknife (subsamples without replacement but with block exclusion)
#num_subsample_iterations = 100
#over_1_snv_SSC_jackknife_TSS <- sapply(seq(1:num_subsample_iterations), function(x) { print(paste0("Subsample iteration: ", x)); sub_snv_counts <- jackknife_snv_count(controls, "TSS"); return(sum(sub_snv_counts$snv_count > 1)) } )
#over_1_snv_SSC_jackknife_TES <- sapply(seq(1:num_subsample_iterations), function(x) { print(paste0("Subsample iteration: ", x)); sub_snv_counts <- jackknife_snv_count(controls, "TES"); return(sum(sub_snv_counts$snv_count > 1)) } )
#over_2_snv_SSC_jackknife_TSS <- sapply(seq(1:num_subsample_iterations), function(x) { print(paste0("Subsample iteration: ", x)); sub_snv_counts <- jackknife_snv_count(controls, "TSS"); return(sum(sub_snv_counts$snv_count > 2)) } )
#over_2_snv_SSC_jackknife_TES <- sapply(seq(1:num_subsample_iterations), function(x) { print(paste0("Subsample iteration: ", x)); sub_snv_counts <- jackknife_snv_count(controls, "TES"); return(sum(sub_snv_counts$snv_count > 2)) } )
#over_1_mut_SSC_jackknife_TSS <- sapply(seq(1:num_subsample_iterations), function(x) { print(paste0("Subsample iteration: ", x)); sub_snv_counts <- jackknife_snv_count(ssc_controls, "TSS"); return(sum(sub_snv_counts$snv_count > 1)) } )
#over_1_mut_SSC_jackknife_TES <- sapply(seq(1:num_subsample_iterations), function(x) { print(paste0("Subsample iteration: ", x)); sub_snv_counts <- jackknife_snv_count(ssc_controls, "TES"); return(sum(sub_snv_counts$snv_count > 1)) } )
#over_2_mut_SSC_jackknife_TSS <- sapply(seq(1:num_subsample_iterations), function(x) { print(paste0("Subsample iteration: ", x)); sub_snv_counts <- jackknife_snv_count(ssc_controls, "TSS"); return(sum(sub_snv_counts$snv_count > 2)) } )
#over_2_mut_SSC_jackknife_TES <- sapply(seq(1:num_subsample_iterations), function(x) { print(paste0("Subsample iteration: ", x)); sub_snv_counts <- jackknife_snv_count(ssc_controls, "TES"); return(sum(sub_snv_counts$snv_count > 2)) } )



#@ Many Jackknife
#over_1_snv_SSC_subsamples_TSS_estimates <- c()
#over_1_snv_SSC_subsamples_TSS_variances <- c()
#for(y in 1:num_subsample_iterations) { 
#	print(paste0("Iteration ", y))
#	controls_350 <- controls[controls$sample %in% sample(unique(paste(controls$sample)), SAMPLE_COUNT+1),]
#	jack <- jackknife_snv_count(controls_350, sites="TSS", block_size=1)
#	over_1_snv_SSC_subsamples_TSS_estimates <- c(over_1_snv_SSC_subsamples_TSS_estimates, jack[["estimate"]])
#	over_1_snv_SSC_subsamples_TSS_variances <- c(over_1_snv_SSC_subsamples_TSS_variances, jack[["var"]])
#}
#TSS_many_jack_mut <- new.env()
#TSS_many_jack_mut[["estimates"]] <- over_1_snv_SSC_subsamples_TSS_estimates
#TSS_many_jack_mut[["vars"]] <- over_1_snv_SSC_subsamples_TSS_variances

