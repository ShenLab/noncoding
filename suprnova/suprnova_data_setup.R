library(keras)
library(kerasR)
library(stringr)
library(pbapply)
#library(EBImage)
library(rhdf5)
library(gtools)
library("PWMEnrich")
#library("TFBSTools")
library("seqLogo")
library("msa")
library("ComplexHeatmap")
library("hexbin")
library("RColorBrewer")
library("ggplot2")
library("tfprobability")
library("tidyverse")
library("stringi")
library("parallel")
library("data.table")
library("abind")
library("reticulate")
source("alex_suite.R")
# //TODO: Not all of the above are dependencies any more, and will be cleaned up.

#################################################################################################################
# PARAMETERS
#################################################################################################################
# Delimiter between folders/files in filepaths. Should be "/" for Linux and "\\" for Windows.
FILEPATH_DELIM <- "/"
# Path to data folder, allowing simpler loading of data with the alex_suite.R data_path function.
DATA_FOLDER <- "/mnt/data/ak3792"
# Path to big data folder, allowing simpler loading of big data with the alex_suite.R bigdata_path function.
BIGDATA_FOLDER <- "/mnt/data/ak3792"
# Path to output folder, allowing simpler writing of output with the alex_suite.R output_path function.
OUTPUT_FOLDER <- "/home/ak3792/Research/ML/output"

#################################################################################################################
# FUNCTIONS TO TRANSFORM BETWEEN GENOMIC SEQUENCES AND THEIR HOT-ENCODED TENSOR REPRESENTATION
#################################################################################################################

# Convert genomic sequences into matrix form.
genomic_sequences_to_matrix <- function(sequences, sequence_annotations=NULL) {
    bases <- c("A","C","G","T")
    sequence_length = max(nchar(sequences))
    num_sequences = length(sequences)
    mat <- matrix(0, nrow=num_sequences, ncol=4*sequence_length); rownames(mat) <- paste0("sequence_",1:num_sequences); colnames(mat) <- apply(expand.grid(bases, 1:sequence_length)[,c(2,1)], 1, paste0, collapse="")
    sequences_split <- strsplit(sequences,"")
    for(i in 1:num_sequences) { 
        sequence <- sequences_split[[i]]
        encoded_sequence <- sapply(1:length(sequence), function(base_position) { (base_position-1)*4 + which(bases == sequence[base_position]) })
        mat[i,c(encoded_sequence)] <- rep(1,length(encoded_sequence))
    }
    # Append sequence annotations to matrix
    if(!is.null(sequence_annotations)) {
        for(i in 1:length(sequence_annotations)) {
            mat <- cbind(mat, sequence_annotations[[i]]); colnames(mat)[ncol(mat)] <- names(sequence_annotations)[i]
        }
    }
    return(mat)
}

# Convert genomic sequences into tensor form.
genomic_sequences_to_tensor <- function(sequences, sequence_annotations=NULL, sequence_length=NULL, frame_width=NULL, num_frames=1, verbose=TRUE) {
    bases <- c("A","C","G","T"); base_index_lookup <- 1:4; names(base_index_lookup) <- bases
    if(is.null(sequence_length)) { sequence_length = max(nchar(sequences)) }
    frame_width_null = is.null(frame_width); num_frames_null = is.null(num_frames)
    if(frame_width_null) { frame_width = sequence_length }
    if(num_frames_null) { num_frames = sequence_length }
    if(frame_width_null && num_frames_null) { num_frames = 1; frame_width = sequence_length }
    num_sequences = length(sequences)
    if(num_frames > 1) { tensor <- array(0, dim=c(num_sequences, num_frames, frame_width, 4, 1), dimnames = list(paste0("sequence_",1:num_sequences), paste0("frame_",1:num_frames), 1:frame_width, bases, "channel")) 
    } else { tensor <- array(0, dim=c(num_sequences, frame_width, 4, 1), dimnames = list(paste0("sequence_",1:num_sequences), 1:frame_width, bases, "channel")) }
    #colnames(mat) <- apply(expand.grid(bases, 1:sequence_length)[,c(2,1)], 1, paste0, collapse="")
    is_biostring <- class(sequences) == "DNAStringSet"
    if(!is_biostring) { sequences <- DNAStringSet(sequences) }
    if(verbose) print("Converting sequences...")
    sequence_tensors <- lapply(1:num_sequences, function(i) { 
        if(verbose) print(paste0(i," / ",num_sequences))
        #if(verbose && i %% 1000 == 1) { print(paste0(i," / ",num_sequences)) }
        sequence <- sequences[[i]]
        if(num_frames > 1) {
            if(length(padded_sequence) < frame_width) {
                tensor[i,1,1:length(sequence),,1] <<- letterFrequencyInSlidingView(sequence, view.width=1, letters=bases)
            } else {
                rollapply(1:length(padded_sequence), width=frame_width, FUN=function(sequence_coords) {
                    sequence <- sequence[sequence_coords]
                    tensor[i,sequence_coords[1],1:length(sequence),,1] <<- letterFrequencyInSlidingView(sequence, view.width=1, letters=bases)
                    return(0)
                })
            }
        } else { tensor[i,1:length(sequence),,1] <<- letterFrequencyInSlidingView(sequence, view.width=1, letters=bases) }
        return(0)
    }) #, mc.cores=detectCores())
    
    return(tensor)
}

# Convert tensor into genomic sequences.
tensor_to_genomic_sequences <- function(tensor) {
    #bases <- c("A","C","G","T",)
    tensor_dimensions = dim(tensor)
    if(length(tensor_dimensions) > 3) { tensor <- tensor[,,,1] 
    } else if( length(tensor_dimensions) < 3) { 
        three_dimensional_tensor <- array(0, dim=c(1, tensor_dimensions))
        three_dimensional_tensor[1,,] <- tensor
        tensor <- three_dimensional_tensor
        tensor_dimensions = dim(tensor)
    }
    sequence_length = tensor_dimensions[2]
    num_sequences = tensor_dimensions[1]
    
    bases <- c("A","C","G","T")
    sequences <- unlist(sapply(1:num_sequences, function(i) {
        sequence_hot_encoding <- tensor[i,,]
        return(paste0(bases[unlist(apply(sequence_hot_encoding, 1, function(base) which(base == 1)))], collapse=""))
    }))
    
    return(sequences)
}

#################################################################################################################
# FUNCTIONS TO UNPACK (AND RE-PACK, FOR rSUPRNOVA) VARIANTS FOR PER-ALLELE CHANGE gnomAD ANNOTATION
#################################################################################################################

# Unpacks each mutation, given reference sequence context annotations, into a mutation per position in each region/sequence context. Useful for regional gnomAD annotation.
unpack_muts <- function(dat, all_alts=TRUE, remove_N_bases=TRUE) {
    dat_colnames <- colnames(dat)
    num_rows = nrow(dat)
    colnames_to_remove <- which(dat_colnames %in% c("Position", "Ref", "Alt"))
    bases <- c("A","C","G","T")
    dat_unpacked <- suppressWarnings(unfactorize(data.frame(rbindlist(lapply(1:num_rows, function(i) {
        print(paste0(i," / ",num_rows))
        region_ref <- unlist(strsplit(dat[i,"ref_sequence"],""))
        region_ref_length = length(region_ref)
        padding = floor(region_ref_length / 2)
        if(remove_N_bases && "N" %in% region_ref) { 
            return(cbind(dat[i,-colnames_to_remove], ".", "N", "N", dat[i,colnames_to_remove]))
        } else if(all_alts) {
            paddings <- floor(seq(-padding, padding+0.9, by=0.25))
            paddings_length = length(paddings)
            dat_expanded <- as.data.frame(lapply(dat[i,], rep, paddings_length))
            return(cbind(dat_expanded[,-colnames_to_remove], dat[i,"Position"] + paddings, unlist(lapply(region_ref, function(x) rep(x,4))), rep(bases,region_ref_length), dat_expanded[,colnames_to_remove]))
        } else {
            return(cbind(dat[i,-colnames_to_remove], dat[i,"Position"] + (-padding):padding, region_ref), dat[i,colnames_to_remove])
        }
    })))))
    colnames(dat_unpacked)[(ncol(dat_unpacked)-(2*length(colnames_to_remove))+1):ncol(dat_unpacked)] <- c("Position", "Ref", "Alt", "Position_original", "Ref_original", "Alt_original")
    dat_unpacked <- dat_unpacked[dat_unpacked$Ref != dat_unpacked$Alt,]
    dat_unpacked <- standardize_colnames(dat_unpacked, re_order=TRUE)
    dat_unpacked$Position <- as.numeric(dat_unpacked$Position)
    return(dat_unpacked)
}

# Groups sequences from unpacked ANNOVAR annotation into tensor objects
dat <- readRDS(output_path("transcribed100k_1_fully_unpacked_annotated.rds")); collapse_variants=FALSE; re_order=FALSE
group_seqs_into_tensors <- function(dat, collapse_variants=TRUE, re_order=TRUE) {
    if(re_order) {
        if("seq" %in% colnames(dat)) { dat <- dat[order(dat$seq, dat$Chrom, dat$Position),]
        } else { dat <- dat[order(dat$Chrom, dat$Position),] }
    }
    dat_num_rows = nrow(dat)
    position_difs <- diff(dat$Position)
    rows_per_position <- which(position_difs > 0)[1]
    frame_width <- which(position_difs > 1)[1] / rows_per_position
    rows_per_frame = frame_width * rows_per_position
    num_sequences = floor(dat_num_rows / rows_per_frame)
    if(num_sequences * rows_per_frame != dat_num_rows) {
        print(paste0("ERROR: ",rows_per_position," * ",frame_width," is not an even multiple of ",dat_num_rows))
        return(1)
    }
    seqs <- unlist(lapply(1:num_sequences, function(seq_i) rep(seq_i, rows_per_frame)))
    sequences <- dat$ref_sequence
    AFs <- dat$AF
    dat_seq_breakpoints <- c(0,which(diff(seqs)>0),nrow(dat))
    if(collapse_variants) { 
        af_tensor <- array(0, dim=c(num_sequences, frame_width), dimnames = list(paste0("sequence_",1:num_sequences), 1:frame_width)) #colnames(mat) <- apply(expand.grid(bases, 1:sequence_length)[,c(2,1)], 1, paste0, collapse="")
    } else { 
        alts <- 1:4; names(alts) <- c("A","C","G","T"); dat_alts <- dat$Alt
        empty_row <- rep(-1,length(alts))
        af_tensor <- array(0, dim=c(num_sequences, frame_width, length(alts)), dimnames = list(paste0("sequence_",1:num_sequences), 1:frame_width, names(alts)))
    }
    tensor <- abind(lapply(1:num_sequences, function(seq_i) {
        print(seq_i)
        seq_start = dat_seq_breakpoints[seq_i]+1
        seq_end = dat_seq_breakpoints[seq_i+1]
        
        if(collapse_variants) { 
            seq_AFs <- rollapply(AFs[seq_start:seq_end], width=rows_per_position, by=rows_per_position, FUN=sum)
            af_tensor[seq_i,1:length(seq_AFs)] <<- seq_AFs
        } else { 
            seq_AFs <- rollapply(seq_start:seq_end, width=rows_per_position, by=rows_per_position, FUN=function(x) {
                filled_row <- empty_row; filled_row[alts[dat_alts[x]]] <- AFs[x]
                return(filled_row)
            })
            af_tensor[seq_i,1:nrow(seq_AFs),] <<- seq_AFs
        }
        return(genomic_sequences_to_tensor(sequences[seq_start], sequence_length=frame_width, verbose=FALSE))
    }), along=4); tensor <- aperm(tensor, c(4,2,3,1))
    
    return_env <- new.env()
    return_env[["input"]] <- tensor
    return_env[["output"]] <- af_tensor
    return(return_env)
}

#################################################################################################################
# FUNCTIONS TO DEFINE REGIONAL ANNOTATION DATA USING GENCODE
#################################################################################################################

# Process GENCODE v34 data and use it to build and save GRanges objects for each type of region (5'UTR, 3'UTR, CDS, splice site, intronic, etc.) 
process_gencode <- function() {
    gencode <- read.csv(data_path("gencode.v34.annotation.gff3"), sep="\t", skip=7)[,-c(6,8)]
    colnames(gencode) <- c("chromosome", "source", "region_type", "start", "end", "strand", "info")
    gencode <- gencode[gencode$chromosome %in% paste0("chr",c(1:22,"X","Y")),] #gencode$source == "ENSEMBL"
    region_types <- c("five_prime_UTR", "three_prime_UTR", "transcript", "exon", "CDS")
    info_to_extract <- c("gene_name", "gene_id", "transcript_id", "exon_number", "remap_original_location")
    empty_info <- t(data.frame(rep(".",length(info_to_extract)))); colnames(empty_info) <- info_to_extract; rownames(empty_info) <- NULL
    info_to_extract_sorted <- sort(info_to_extract)
    b <- rbindlist(lapply(region_types, function(region_type) {
        print(region_type)
        region_type_indices <- which(gencode$region_type == region_type)
        a <- rbindlist(lapply(strsplit(paste0(gencode$info[region_type_indices]), ";"), function(region_info) {
            region_info <- unfactorize(data.frame(strsplit(region_info, "=")))
            colnames(region_info) <- region_info[1,]
            region_info <- cbind(region_info[2,intersect(info_to_extract,colnames(region_info))],empty_info); region_info <- region_info[,!duplicated(colnames(region_info))]
            return(region_info[,info_to_extract])
        }))
        a <- cbind(gencode[region_type_indices,],a)
        saveRDS(a, output_path(paste0("gencode_",region_type,"_dat_hg38.rds")))
        a_granges <- to_genomic_regions(a, label_colname="gene_name")
        saveRDS(a_granges, output_path(paste0("gencode_",region_type,"_granges_hg38.rds")))
        return(a)
    }))
    saveRDS(b, output_path(paste0("gencode_full_dat_hg38.rds")))
    
    cds_granges <- readRDS(output_path("gencode_CDS_granges_hg38.rds"))
    transcript_granges <- readRDS(output_path("gencode_transcript_granges_hg38.rds"))
    exon_granges <- readRDS(output_path("gencode_exon_granges_hg38.rds"))
    five_prime_UTR_granges <- readRDS(output_path("gencode_five_prime_UTR_granges_hg38.rds"))
    three_prime_UTR_granges <- readRDS(output_path("gencode_three_prime_UTR_granges_hg38.rds"))
    
    intron_granges <- GenomicRanges::setdiff(transcript_granges, GenomicRanges::union(cds_granges, GenomicRanges::union(five_prime_UTR_granges, three_prime_UTR_granges)))
    saveRDS(intron_granges, output_path(paste0("gencode_intron_granges_hg38.rds")))
    intron_widths <- width(intron_granges)
    splice_site_distance_to_edge = 100
    
    five_prime_ss_granges <- intron_granges
    end(five_prime_ss_granges[strand(five_prime_ss_granges) == "+"]) <- pmin(start(five_prime_ss_granges[strand(five_prime_ss_granges) == "+"]) + (splice_site_distance_to_edge - 1), end(five_prime_ss_granges[strand(five_prime_ss_granges) == "+"]))
    start(five_prime_ss_granges[strand(five_prime_ss_granges) == "-"]) <- pmax(end(five_prime_ss_granges[strand(five_prime_ss_granges) == "-"]) - (splice_site_distance_to_edge - 1), start(five_prime_ss_granges[strand(five_prime_ss_granges) == "-"]))
    saveRDS(five_prime_ss_granges, output_path(paste0("gencode_five_prime_ss_granges_hg38.rds")))
    three_prime_ss_granges <- intron_granges
    start(three_prime_ss_granges[strand(three_prime_ss_granges) == "+"]) <- pmax(end(three_prime_ss_granges[strand(three_prime_ss_granges) == "+"]) - (splice_site_distance_to_edge - 1), start(three_prime_ss_granges[strand(three_prime_ss_granges) == "+"]))
    end(three_prime_ss_granges[strand(three_prime_ss_granges) == "-"]) <- pmin(start(three_prime_ss_granges[strand(three_prime_ss_granges) == "-"]) + (splice_site_distance_to_edge - 1), end(three_prime_ss_granges[strand(three_prime_ss_granges) == "-"]))
    saveRDS(three_prime_ss_granges, output_path(paste0("gencode_three_prime_ss_granges_hg38.rds")))
    
    noncoding_exon_granges <- GenomicRanges::setdiff(exon_granges, cds_granges)
    saveRDS(noncoding_exon_granges, output_path(paste0("gencode_noncoding_exon_granges_hg38.rds")))
}

# Return a data frame of a specified number of randomly sampled variants from gnomAD, with optional limitation to a certain desired region type.
sample_gnomad_variants <- function(N, region_type=NULL, af_cutoff=1e-3) {
    print("Reading gnomAD data...")
    gnomad_trimers <- readRDS(data_path("hg38_gnomad3.0_genome_snvs_trimers.rds"))
    if(!is.null(af_cutoff)) {  gnomad_trimers <- gnomad_trimers[gnomad_trimers$AF < af_cutoff,] }
    if(is.null(region_type)) { num_gnomAD_variants_to_sample = N } else { num_gnomAD_variants_to_sample <- N * 1000 }
    
    print("Sampling gnomAD variants...")
    if(num_gnomAD_variants_to_sample > nrow(gnomad_trimers)) {
        num_gnomAD_variants_to_sample = nrow(gnomad_trimers)
        gnomad_dat <- unfactorize(data.frame(gnomad_trimers))
    } else { gnomad_dat <- unfactorize(data.frame(gnomad_trimers[sample(1:nrow(gnomad_trimers), num_gnomAD_variants_to_sample),])) }
    
    if(!is.null(region_type)) {
        gnomad_dat_annot <- annotate_region_types(gnomad_dat, "hg38", region_type)
        gnomad_dat <- gnomad_dat[sample(which(gnomad_dat_annot[,region_type] == TRUE), min(c(N,num_gnomAD_variants_to_sample)))]
    }
    return(gnomad_dat)
}

# Annotate the data with the specified region types. 
# CDS / coding sequence can be optionally broken up into syn/mis/nonsense if breakdown_CDS is set to TRUE, but note this part of the code is still in testing.
annotate_region_types <- function(dat, version="hg19", region_types=c("CDS","five_prime_UTR","three_prime_UTR","five_prime_ss","three_prime_ss","intron","intergenic"), breakdown_CDS=FALSE, num_variants_to_sample=NULL, label_colname="ref_sequence") {
    if(is.null(num_variants_to_sample) || num_variants_to_sample > nrow(dat_trimers)) {
        num_variants_to_sample = nrow(dat)
        a <- unfactorize(data.frame(dat))
    } else { a <- unfactorize(data.frame(dat[sample(1:nrow(dat), num_variants_to_sample),])) }
    a_indices <- 1:num_variants_to_sample
    if(version == "hg38") {
        a_granges <- to_genomic_regions(a, chr_colname="Chrom", start_colname="Position") #, label_colname=label_colname)
    } else if(sum(c("Chrom_hg38","Position_hg38") %in% colnames(a)) == 2) {
        a_granges <- to_genomic_regions(a, chr_colname="Chrom_hg38", start_colname="Position_hg38") #, label_colname=label_colname)
    } else { print("Version should be hg38 or have Chrom_hg38 and Position_hg38 colnames for this function (for now)!"); return(1) }
    
    region_annotation <- lapply(region_types, function(region_type) {
        print(region_type)
        if(region_type == "CDS") { reg_granges <- readRDS(output_path("gencode_CDS_granges_hg38.rds"))
        } else if(region_type == "transcript") { reg_granges <- readRDS(output_path("gencode_transcript_granges_hg38.rds"))
        } else if(region_type == "exon") { reg_granges <- readRDS(output_path("gencode_exon_granges_hg38.rds"))
        } else if(region_type == "five_prime_UTR") { reg_granges <- readRDS(output_path("gencode_five_prime_UTR_granges_hg38.rds"))
        } else if(region_type == "three_prime_UTR") { reg_granges <- readRDS(output_path("gencode_three_prime_UTR_granges_hg38.rds"))
        } else if(region_type == "five_prime_ss") { reg_granges <- readRDS(output_path(paste0("gencode_five_prime_ss_granges_hg38.rds")))
        } else if(region_type == "three_prime_ss") { reg_granges <- readRDS(output_path(paste0("gencode_three_prime_ss_granges_hg38.rds")))
        } else if(region_type == "intron") { reg_granges <- readRDS(output_path(paste0("gencode_intron_granges_hg38.rds")))
        } else if(region_type == "noncoding_exon") { reg_granges <- readRDS(output_path(paste0("gencode_noncoding_exon_granges_hg38.rds"))) # 5'UTR and 3'UTR
        } else { return(NULL) }
        
        is_region <- a_indices %in% queryHits(findOverlaps(a_granges, reg_granges))
        if(region_type == "CDS" && breakdown_CDS) {
            if(sum(is_region) > 0) {
                cds_refGene <- a[is_region,]; cds_refGene <- cbind(cds_refGene, paste0("CDS_",1:nrow(cds_refGene))); colnames(cds_refGene)[ncol(cds_refGene)] <- "sample"
                cds_refGene <- run_annovar(cds_refGene, "refGene", "hg38")
                is_cds_missense <- a_indices %in% which(is_cds)[cds_refGene$ExonicFunc.refGene == "nonsynonymous SNV"]
                is_cds_syn <- a_indices %in% which(is_cds)[cds_refGene$ExonicFunc.refGene == "synonymous SNV"]
                is_cds_nonsense <- a_indices %in% which(is_cds)[cds_refGene$ExonicFunc.refGene == "stopgain"]
                is_cds_startloss <- a_indices %in% which(is_cds)[cds_refGene$ExonicFunc.refGene == "startloss"]
                is_cds_stoploss <- a_indices %in% which(is_cds)[cds_refGene$ExonicFunc.refGene == "stoploss"]
            } else { is_cds_missense <- is_region; is_cds_syn <- is_region; is_cds_nonsense <- is_region; is_cds_startloss <- is_region; is_cds_stoploss <- is_region }
            is_region <- data.frame(is_region, is_cds_missense, is_cds_syn, is_cds_nonsense, is_cds_startloss, is_cds_stoploss)
            colnames(is_region) <- c(region_type, "missense", "synonymous", "nonsense", "startloss", "stoploss")
        }
        else { 
            is_region <- data.frame(is_region)
            colnames(is_region) <- c(region_type)
        }
        return(is_region)
    }); names(region_annotation) <- region_types
    if("intergenic" %in% region_types) { 
        is_intergenic <- !(region_annotation$CDS[,1] | region_annotation$five_prime_UTR | region_annotation$three_prime_UTR | region_annotation$intron) 
        colnames(is_intergenic) <- "intergenic"
        region_annotation[["intergenic"]] <- is_intergenic
    }
    
    region_annotation <- Reduce(cbind, region_annotation)
    return(region_annotation)
}

# Return a data frame of randomly sampled variants for each desired region type, with a specified number to sample per type.
sample_region_type_variants <- function(region_types = c("five_prime_UTR", "three_prime_UTR", "five_prime_ss", "three_prime_ss", "intron", "CDS"), regions_to_sample_per_type = 20000, width = 151) {
    sampled_region_type_variants <- unfactorize(rbindlist(lapply(region_types, function(region_type) {
        region_type_granges <- readRDS(output_path(paste0("gencode_",region_type,"_granges_hg38.rds")))
        sampled_indices <- 1:length(region_type_granges)
        sampled_indices <- sample(1:length(region_type_granges), regions_to_sample_per_type, replace=TRUE)
        sampled_regions_for_type <- rbindlist(lapply(1:length(sampled_indices), function(i) {
            if(i %% 1000 == 1) { print(paste0("Sampling ",region_type,": ",i," / ",length(sampled_indices))) }
            sampled_i = sampled_indices[i]
            region <- region_type_granges[sampled_i]
            region_chrom = seqnames(region); region_start = start(region); region_end = end(region)
            possible_starts <- region_start:(region_end-width+1)
            if(length(possible_starts) > 1 && possible_starts[2] > possible_starts[1]) { # Region is longer than width
                region_start <- sample(possible_starts, 1); region_end <- region_start + width - 1
                start <- region_start; end <- region_end
            } else { 
                midpoint <- region_start + floor((region_end - region_start + 1)/2)
                arm_width = floor(width/2)
                start <- midpoint - arm_width; end <- start + width - 1
            }
            random_region <- data.frame(unlist(c(region_type, i, region_chrom, start, end, region_start, region_end)))
            colnames(random_region) <- c("region_type", "index", "Chrom", "start", "end", "region_start", "region_end")
            return(random_region)
        }))
    })))
    return(sampled_region_type_variants)
}

# Get particularly defined TS regions, given the genebody and left and right pillow sizes around the specified "TSS" or "TES" sites.
# This function is deprecated in favor of the results from process_gencode() for overall use, but is still in the code for simple customizable TSS/TES definitions.
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
#TSS_regions <- get_TS_regions(genebody, sites="TSS", left=20000, right=20000); TES_regions <- get_TS_regions(genebody, sites="TES", left=5000, right=20000)
#TSS_granges <- to_genomic_regions(TSS_regions, labels=TSS_regions$gene); TES_granges <- to_genomic_regions(TES_regions, labels=TES_regions$gene)

#################################################################################################################
# FUNCTIONS TO CALCULATE GENOMIC ANNOTATION FEATURES
#################################################################################################################

# Finds gnomAD o/e (observed/expected) gene mutation constraints. This metric is more parallel to noncoding allele frequency constraint than pLI
get_obs_over_exp <- function(genebody) {
    genebody <- read.csv(data_path("refGene\\refGene_hg19_genebody_fixed.bed"), sep="\t", stringsAsFactors=FALSE)
    genebody <- unfactorize(genebody[!duplicated(genebody$gene),])
    #genebody_granges_hg19 <- to_genomic_regions(genebody_hg19, chr_colname="chromosome", start_colname="tss", end_colname="tes", strand_colname="strand", label_colname="gene", order_coordinates=TRUE, remove_duplicate_labels=TRUE)
    ordered_coordinates <- t(apply(genebody[,c("tss", "tes")], 1, sort))
    #genebody$tss <- ordered_coordinates[,1]
    #genebody$tes <- ordered_coordinates[,2]
    genes <- cbind(genebody[,1], apply(genebody, 1, function(x) {
        return(round(sum(as.numeric(x[2:3]))/2))
    }), rep("?",nrow(genebody)), ordered_coordinates, genebody[,-c(1:3)])
    colnames(genes) <- c("Chrom", "Position", "Ref", "start", "end", colnames(genebody)[-c(1:3)])
    genes <- standardize_colnames(genes, re_order=TRUE, remove_chr_prefix=TRUE)
    genes <- genes[genes$Chrom %in% c(1:22),]
    widths <- abs(genes$end - genes$start) + 1
    genes <- add_sequence_context_feature(genes, version="hg19", width=widths)
    saveRDS(genes, output_path("genebody.rds"))
    
    colnames(genes)[which(colnames(genes) %in% c("Position", "ref_sequence"))] <- c("Position_hg19", "ref_sequence_hg19")
    genes <- liftover(genes, from="hg19", to="hg38", chr_colname="Chrom", start_colname="start", end_colname="end", ref_colname="Ref", alt_colname="Alt", confirm_refseq=FALSE, mismatches_pause=FALSE)
    genes <- genes[genes$Chrom %in% c(1:22),] # & (genes$end - genes$start + 1 == width)
    genes <- cbind(genes, apply(genes, 1, function(x) {
        return(round(sum(as.numeric(x[c("start","end")]))/2))
    })); colnames(genes)[ncol(genes)] <- "Position"
    widths <- abs(genes$end - genes$start) + 1
    genes <- add_sequence_context_feature(genes, version="hg38", width=widths)
    genes <- genes[,c("Chrom", "Position", "start", "end", "gene", "transcript", "strand", "Ref", "ref_sequence", "Chrom_hg19", "Position_hg19", "start_hg19", "end_hg19")] #genes$ref_sequence == genes$ref_sequence_hg19
    nrow(genes)
    saveRDS(genes, output_path("genebody.rds"))
    
    genes <- readRDS(output_path("genebody.rds"))
    genes_granges <- to_genomic_regions(genes, chr_colname="Chrom", label_colname="gene")
    gnomad_exome_granges <- to_genomic_regions(read.csv("/data/annovar/humandb/hg38_gnomad_exome.txt", sep="\t"), chr_colname="X.Chr", start_colname="Start", end_colname="End", label_colname="gnomAD_exome_ALL")
    gnomad_exome_granges <- gnomad_exome_granges[width(gnomad_exome_granges) == 1]
    gnomad_exome_granges_AF <- as.numeric(names(gnomad_exome_granges)); gnomad_exome_granges_AF[is.na(gnomad_exome_granges_AF)] <- 0
    gnomad_gene_overlaps <- data.frame(findOverlaps(genes_granges, gnomad_exome_granges))
    observed <- aggregate(subjectHits ~ queryHits, data=gnomad_gene_overlaps, FUN=function(x) { AFs <- gnomad_exome_granges_AF[x]; return(mean(AFs[AFs < 1e-3])) })
    colnames(observed) <- c("gene", "AF")
    observed$gene <- names(genes_granges)[observed$gene]
    #sample_size = 71702; observed$AC <- round(observed$AC * sample_size * 2 + 0.01)
    genes <- merge(genes, observed, all.x=TRUE)
    saveRDS(genes, output_path("genes.rds"))
    
    expected <- get_expected_muts(genes, chr_colname="Chrom_hg19", start_colname="start_hg19", end_colname="end_hg19", aggregate_region=TRUE, version="hg19")
    genes <- cbind(genes, expected)
    saveRDS(genes, output_path("genes.rds"))
    
    genes$AF[is.na(genes$AF)] <- 0
    genes$AF <- genes$AF + (min(genes$AF[genes$AF > 0])/2)
    genes$expected[genes$expected < 1e-8] <- 0
    o_e_Z <- genes$AF / (genes$expected + (min(genes$expected[genes$expected > 0])/2))
    o_e_Z[o_e_Z > 2500] <- 2500
    plot(density(genes$expected))
    plot(density(genes$AF))
    plot(density(o_e_Z))
    o_e_Z <- norm_tensor(o_e_Z)
    plot(density(o_e_Z))
    genes <- cbind(genes, o_e_Z)
    #genes <- cbind(genes[,-which(colnames(genes)=="o_e_Z")], o_e_Z)
    saveRDS(genes, output_path("genes.rds")) 
    #genes <- readRDS(output_path("genes.rds")) 
    exac_dat <- read.csv(data_path("fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt"), sep="\t")
    genes <- merge(genes, exac_dat[!duplicated(exac_dat$gene),c("gene","pLI")], all.x=TRUE)
    #genes <- merge(genes[,-which(colnames(genes)=="pLI")], exac_dat[!duplicated(exac_dat$gene),c("gene","pLI")], all.x=TRUE)
    #draw_plot(data.frame(x=genes$pLI[genes$o_e_Z < 1.5 & !is.na(genes$pLI)], y=genes$o_e_Z[genes$o_e_Z < 1.5 & !is.na(genes$pLI)]), hex_density = 25, title="Obs/Exp Z-score vs. pLI", xlab="ExAC pLI", ylab="gnomAD 3.0 obs/exp Z", legend_text="# genes", linear_best_fit=TRUE, quadratic_best_fit=FALSE, filename="obs_exp_Z_vs_pLI.pdf", cor_method="Spearman")
    genes$pLI[is.na(genes$pLI)] <- mean(genes$pLI[!is.na(genes$pLI)])
    cor(genes$pLI, genes$o_e_Z, method="spearman")
    saveRDS(genes, output_path("genes.rds")) 
    
    # Convert genes data frame into one with all positions and mutation types fully unpacked and ready to be inputted for ANNOVAR gnomAD annotation.
    genes <- readRDS(output_path("genebody.rds"))
    genes_fully_unpacked <- unpack_muts(genes)
    num_seqs = length(unique(genes_fully_unpacked$gene)); num_seqs
    seq_row_counts <- table(table(genes_fully_unpacked$genes)); seq_row_counts
    rows_per_seq = as.numeric(names(seq_row_counts))
    rows_per_seq == width*3
    genes_fully_unpacked_midpoints <- seq((arm_width+1)*3,nrow(genes_fully_unpacked),by=width*3)-1
    genes_hg38_to_hg19_pos_shifts <- genes_fully_unpacked$Position_hg19[genes_fully_unpacked_midpoints]-genes_fully_unpacked$Position[genes_fully_unpacked_midpoints]
    genes_fully_unpacked$Position_hg19 <- genes_fully_unpacked$Position + unlist(lapply(genes_hg38_to_hg19_pos_shifts, function(x) rep(x,rows_per_seq)))
    saveRDS(genes_fully_unpacked, output_path("genes_fully_unpacked.rds"))
    genes_fully_unpacked <- readRDS(output_path("genes_fully_unpacked.rds"))
    gnomad_fully_unpacked_annotated <- run_annovar(gnomad_fully_unpacked, "gnomad30_genome", buildver="hg38")
    
    # Annotate fully unpacked gnomad variants data frame with gnomAD 3.0.0, using local ANNOVAR install.
    gnomad_fully_unpacked <- readRDS(output_path("gnomad_fully_unpacked.rds"))
    gnomad_fully_unpacked_annotated <- run_annovar(gnomad_fully_unpacked, "gnomad30_genome", buildver="hg38")
    gnomad_fully_unpacked_annotated$AF[is.na(gnomad_fully_unpacked_annotated$AF) | gnomad_fully_unpacked_annotated$AF == "."] <- 0
    gnomad_fully_unpacked_annotated$AF <- as.numeric(gnomad_fully_unpacked_annotated$AF)
    table(table(gnomad_fully_unpacked_annotated$seq))
    sum(gnomad_fully_unpacked_annotated$AF == 0)/nrow(gnomad_fully_unpacked_annotated)
    saveRDS(gnomad_fully_unpacked_annotated, output_path("gnomad_fully_unpacked_annotated.rds"))
    
    exac_dat <- read.csv(data_path("fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt"), sep="\t")
    nearest_genes_dat <- names(genebody_granges_hg19)[nearest(dat_variant_granges_hg19, genebody_granges_hg19)]
    pLIs_dat <- genes[unlist(sapply(nearest_genes_dat,function(x) { nearest_gene <- which(paste0(exac_dat$gene) == x); if(length(nearest_gene) < 1) { return(NA) } else { return(nearest_gene[1]) }  })),"o_e_Z"]
    pLIs_dat[is.na(pLIs_dat)] <- mean(pLIs_dat[!is.na(pLIs_dat)])
    pLIs_dat <- array(pLIs_dat, c(length(pLIs_dat),1))
    saveRDS(pLIs_dat, output_path("gnomad_pLIs.rds"))
    # Position-specific expected (background mutation rate)
    expected <- get_expected_muts(dat[dat_midpoints,], chr_colname="Chrom_hg19", pos_colname="Position_hg19", width=151, version="hg19")
    #expected <- t(array(expected, dim(dat_gradcams)[2:1]))
    sum(expected== 0)/prod(dim(expected))
    saveRDS(expected, output_path("gnomad_expected.rds"))
}

#################################################################################################################
# FUNCTIONS TO ANNOTATE DATA WITH USEFUL GENOMIC FEATURES
#################################################################################################################

# Set up features_env environment that will store named features (filepaths, data frames, or GRanges objects)
setup_features_env <- function(make_fullPeak=FALSE, features_to_make_fullPeak=NULL) {
    features_env <- new.env()
    # Add Roadmap histone modification features
    HISTONE_MODIFICATIONS_ANNOTATIONS_PATH = data_path("Roadmap")
    for(histone_modification_file in list.files(path=HISTONE_MODIFICATIONS_ANNOTATIONS_PATH)) { 
        features_env_key = paste0(gsub("-", ".", histone_modification_file))
        features_env[[features_env_key]] <- new.env()
        features_env[[features_env_key]][["peaks"]] <- full_path(HISTONE_MODIFICATIONS_ANNOTATIONS_PATH, histone_modification_file)
        histone_mark_info <- strsplit(features_env_key, "\\.")[[1]]
        if(histone_mark_info[1] %in% get_relevant_roadmap_eids("ASD")) { features_env[[features_env_key]][["notes"]] <- paste0("histone_mark; ",histone_mark_info[2]) } else { features_env[[features_env_key]][["notes"]] <- "" }
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
    
    # Add CDH, CHD, autism, HDE, HHE, HBE, constrained gene features
    genebody <- read.csv(data_path("refGene\\refGene_hg19_genebody_fixed.bed"), sep="\t", stringsAsFactors=FALSE)
    genebody[,c("tss","tes")] <- t(sapply(1:nrow(genebody), function(i) { sort(c(genebody$tss[i], genebody$tes[i])) }))
    # Liftover to hg38
    #genebody <- liftover(genebody, from="hg19", to="hg38", chr_colname="chromosome", start_colname="tss", end_colname="tes", confirm_refseq=FALSE)
    genebody_granges <- to_genomic_regions(genebody, chr_colname="chromosome", start_colname="tss", end_colname="tes", strand_colname="strand", label_colname="gene", order_coordinates=TRUE, remove_duplicate_labels=TRUE)
    # Grab gene sets
    candidate_cdh_genes <- unique(read.csv(data_path("Table S5 CDH training genes with mouse model.csv"), stringsAsFactors=FALSE, header=TRUE, skip=1)[,1]) # 61 candidate CDH genes, from Lan's file
    candidate_chd_genes <- unique(read.csv(data_path("mouse_CHD_gene0516.txt"), sep="\t", stringsAsFactors=FALSE)[,1]) # Latest list of 727 candidate CHD genes.
    known_chd_mouse_ko_genes <- unique(read.csv(data_path("Curated known Human-Mouse CHD genes.txt"), sep="\t", stringsAsFactors=FALSE)[,1]) # 253 known CHD mouse KO genes, from Lan's file
    autism_genes <- unique(c(paste0(read.csv(data_path("sfari_high_confidence_genes_112118.csv"))$gene.symbol), paste0(read.csv(data_path("satterstrom_fdr0.1_asd_genes.csv"))$gene))) # 86 autism genes with score 1 or 2 in SFARI database, compiled by Siying
    genes_by_expression <- get_genes_by_expression()
    HDE_genes <- genes_by_expression[["HDE"]]; HHE_genes <- genes_by_expression[["HHE"]]; HBE_genes <- genes_by_expression[["HBE"]]
    constrained_genes <- ls(get_constrained_genes("pLI>0.5"))
    #bivalent_genes <- get_bivalent_genes()
    # Store to features_env, as both hg19 and hg38
    for(features_env_key in c("candidate_cdh_genes", "candidate_chd_genes", "known_chd_mouse_ko_genes", "autism_genes", "HDE_genes", "HHE_genes", "HBE_genes", "constrained_genes")) {
        gene_set <- get(features_env_key)
        features_env[[features_env_key]] <- new.env()
        features_env[[features_env_key]][["peaks"]] <- genebody_granges[names(genebody_granges) %in% gene_set]
        features_env[[features_env_key]][["notes"]] <- "gene_set; default_padding=0"
    }
    
    set_global(features_env)
    
    # Append .fullPeak features
    if(make_fullPeak | !is.null(features_to_make_fullPeak)) {
        if(is.null(features_to_make_fullPeak)) {
            features_to_make_fullPeak <- c(get_features_by_group("H3K36me3"), get_features_by_group("H3K79me2"))
            features_to_make_fullPeak <- features_to_make_fullPeak[grepl("broadPeak",features_to_make_fullPeak)]
        }
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
            annotations_full <- intersect(annotations_full, genebody_granges, ignore.strand=TRUE)
            annotations_full <- union(annotations_full, annotations) # Same as original broadPeak, but with reasonably sized gaps in genebody filled in
            annotations_full <- intersect(annotations_full, annotations_full)
            
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

# Returns all available features that belong to the specified group, such as "RBP" (across types of RBPs), "H3K36me3" (across tissue types), etc...Supports padding specification, just like the annotate function.
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

# Loads the specified feature from features_env.
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

# Annotate the given variants (passed in as either data frame, GRanges, or variant_list format) with the specified features (saved in features_env as either GRanges, data frames, or filepaths)
annotate <- function(variants, features, features_env=NULL, variants_granges=NULL, save_globally=FALSE) {
    if(is.null(features_env)) { features_env <- get_global("features_env") }; if (is.null(features_env)) { features_env <- setup_features_env() }
    if(class(variants) == "GRanges") { variants_granges <- variants; variants <- genomic_regions_to_dat(variants) } # handles case where variants are passed in directly as GRanges
    variants <- standardize_colnames(variants)
    if (is.null(variants_granges)) { variants_granges <- to_genomic_regions(variants, chr_colname="Chrom", start_colname="Position", end_colname="Position") }
    num_features = length(features)
    available_features <- ls(features_env)
    available_features_notes <- sapply(available_features, function(x) { return(unlist(strsplit(features_env[[x]][["notes"]], "; *"))) })
    
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

# Setup allele frequency-annotated data for SUPRNOVA under the given name. 
# The function currently assumes that data is in hg38 but has a Position_hg19 column; typically, run setup_supermodel_data() function before this one.
# Run this function multiple times (with specified number of tasks to complete) to ensure that all setup is complete; this is designed to manage RAM for large datasets.
annotate_AFs <- function(dat_name, dat=NULL, version="hg38", type="r", num_tasks_to_do=1, rewrite=FALSE) {
    tasks_done = 0
    filename = output_path(paste0(dat_name,"_fully_unpacked.rds"))
    if(tasks_done < num_tasks_to_do && (!file.exists(filename) || rewrite)) {
        if(is.null(dat)) { dat <- readRDS(output_path(paste0(dat_name,"_dat.rds"))) }
        regional_fully_unpacked <- unpack_muts(dat)
        num_seqs = length(unique(paste0(regional_fully_unpacked$seq))); num_seqs
        seq_row_counts <- table(table(regional_fully_unpacked$sample)); seq_row_counts
        rows_per_seq = as.numeric(names(seq_row_counts))
        rows_per_seq = width*3
        rows_per_seq == width*3
        regional_fully_unpacked_midpoints <- seq((arm_width+1)*3,nrow(regional_fully_unpacked),by=width*3)-1
        regional_hg38_to_hg19_pos_shifts <- regional_fully_unpacked$Position_hg19[regional_fully_unpacked_midpoints]-regional_fully_unpacked$Position[regional_fully_unpacked_midpoints]
        regional_fully_unpacked$Position_hg19 <- regional_fully_unpacked$Position + unlist(lapply(regional_hg38_to_hg19_pos_shifts, function(x) rep(x,rows_per_seq)))
        saveRDS(regional_fully_unpacked, filename)
        tasks_done = tasks_done + 1
    } else if(tasks_done < num_tasks_to_do) { print(paste0(filename," already defined!")) }
    
    # Annotate fully unpacked variants data frame with gnomAD 3.0.0, using local ANNOVAR install.
    filename = output_path(paste0(dat_name,"_fully_unpacked_annotated.rds"))
    if(tasks_done < num_tasks_to_do && (!file.exists(filename) || rewrite)) {
        regional_fully_unpacked <- readRDS(output_path(paste0(dat_name,"_fully_unpacked.rds")))
        regional_fully_unpacked_annotated <- run_annovar(regional_fully_unpacked, "gnomad30_genome", buildver="hg38")
        regional_fully_unpacked_annotated$AF[is.na(regional_fully_unpacked_annotated$AF) | regional_fully_unpacked_annotated$AF == "."] <- 0
        regional_fully_unpacked_annotated$AF <- as.numeric(regional_fully_unpacked_annotated$AF)
        table(table(regional_fully_unpacked_annotated$seq))
        sum(regional_fully_unpacked_annotated$AF == 0)/nrow(regional_fully_unpacked_annotated)
        saveRDS(regional_fully_unpacked_annotated, filename)
        tasks_done = tasks_done + 1
    } else if(tasks_done < num_tasks_to_do) { print(paste0(filename," already defined!")) }
    
    # Repackage annotated gnomad data into tensor objects.
    collapse_variants = (type != "v")
    if(!collapse_variants) { filename = output_path(paste0(dat_name,"_AFs_variant.rds")) } else { filename = output_path(paste0(dat_name,"_AFs.rds")) }
    if(tasks_done < num_tasks_to_do && (!file.exists(filename) || rewrite)) {
        regional_fully_unpacked_annotated <- readRDS(output_path(paste0(dat_name,"_fully_unpacked_annotated.rds")))
        regional_tensors_env <- group_seqs_into_tensors(regional_fully_unpacked_annotated, collapse_variants=collapse_variants, re_order=FALSE)
        dat_tensor <- regional_tensors_env[["input"]]
        af_tensor <- regional_tensors_env[["output"]]
        rm(regional_tensors_env)
        # dat_tensor[1,,,]; sum(af_tensor== 0)/prod(dim(af_tensor)); sum(rowSums(af_tensor)== 0)/nrow(af_tensor)
        saveRDS(dat_tensor, output_path(paste0(dat_name,"_tensor2.rds")))
        saveRDS(af_tensor, filename)
        tasks_done = tasks_done + 1
    } else if(tasks_done < num_tasks_to_do) { print(paste0(filename," already defined!")) }
    
    print(paste0(tasks_done," tasks completed!"))
    if(tasks_done > 0) { print("Exiting for now. Reboot R and re-run setup function again.") } else { print("Setup already complete!") }
}

# Returns the expected number of mutations for each width-sized region around the specified positions.
get_expected_muts <- function(dat, chr_colname="Chrom", pos_colname="Position", width=151, start_colname=NULL, end_colname=NULL, version="hg19", aggregate_region=FALSE, background_rate_path=data_path("background_mut_rate")) {
    if(version!="hg19") { dat <- liftover(dat, from=version, to="hg19", chr_colname=chr_colname, start_colname=pos_colname, confirm_refseq=FALSE, mismatches_pause=FALSE) }
    if(aggregate_region) { expected_muts <- rep(0, nrow(dat))
    } else { expected_muts <- array(0, c(nrow(dat), width)) }
    if(is.null(start_colname) || is.null(end_colname)) {
        arm_width = floor(width/2)
        all_starts <- as.numeric(paste0(dat[,pos_colname])) - arm_width; all_ends <- as.numeric(paste0(dat[,pos_colname])) + arm_width 
    } else {
        all_starts <- as.numeric(paste0(dat[,start_colname])); all_ends <- as.numeric(paste0(dat[,end_colname]))
    }
    chromosomes = intersect(1:22, unique(dat[,chr_colname]))
    for(chromosome in chromosomes) {
        chr_background_path <- full_path(background_rate_path, paste0("chr",chromosome,"_n.wig"))
        cat(paste0("Reading chr", gsub("chr","",chromosome), "..."))
        chr_background <- readLines(chr_background_path)
        cat("Done.\n")
        chr_indices <- which(dat[,chr_colname] == chromosome)
        if(length(chr_indices)>0) { 
            chr_expected_muts <- sapply(chr_indices, function(i) { as.numeric(paste0("0",chr_background[all_starts[i]:all_ends[i]])) })
            if(aggregate_region) {
                expected_muts[chr_indices] <- unlist(lapply(chr_expected_muts, function(x) { mean(x[which(!is.na(x))]) }))
            } else {
                expected_muts[chr_indices,] <- t(chr_expected_muts)
            }
        }
    }
    return(expected_muts)
}

# Finds overlaps with RBP eCLIP data. Make sure the ranges given are in hg19 to match the eCLIP data!!!!
find_eclip_overlaps <- function(dat_granges=NULL, dat=NULL, chr_colname="chromosome", start_colname="start", end_colname="end", label_colname="gene", eclip_padding=0) {
    if(is.null(dat_granges)) {
        if(is.null(dat)) { print("ERROR: No genomic ranges or data frame provided!"); return(1) }
        dat_granges <- to_genomic_regions(dat, chr_colname=chr_colname, start_colname=start_colname, end_colname=end_colname, label_colname=label_colname)
        #dat_granges <- dat_granges_hg19[order(as.numeric(gsub("seq","",names(dat_granges))))]
    }
    eclip_granges_filename = output_path(paste0("all_rbp_eclip_peaks_",eclip_padding,"bp.rds"))
    if(file.exists(eclip_granges_filename)) { all_rbp_eclip_peaks <- readRDS(eclip_granges_filename) 
    } else {
        print("Loading RBP eCLIP peaks...")
        rbp_eclip_peaks <- sapply(get_features_by_group("RBP"), function(rbp) { print(rbp); x <- load_annotation(paste0(rbp,"_",eclip_padding,"bp")); return(intersect(x,x)) })
        rbp_names <- names(rbp_eclip_peaks)
        all_rbp_eclip_peaks <- Reduce(c, rbp_eclip_peaks)
        names(all_rbp_eclip_peaks) <- unlist(lapply(1:length(rbp_eclip_peaks), function(rbp_i) {
            rep(rbp_names[rbp_i], length(rbp_eclip_peaks[[rbp_i]]))
        }))
        saveRDS(all_rbp_eclip_peaks, eclip_granges_filename)
    }
    dat_eclip_overlaps <- data.frame(findOverlaps(dat_granges, all_rbp_eclip_peaks))
    dat_eclip_starts <- start(dat_granges[dat_eclip_overlaps$queryHits])
    dat_eclip_overlaps <- cbind(dat_eclip_overlaps, names(all_rbp_eclip_peaks)[dat_eclip_overlaps$subjectHits], start(all_rbp_eclip_peaks[dat_eclip_overlaps$subjectHits])-dat_eclip_starts, end(all_rbp_eclip_peaks[dat_eclip_overlaps$subjectHits])-dat_eclip_starts)
    colnames(dat_eclip_overlaps) <- c("query", "subject", "RBP", "start", "end")
    return(dat_eclip_overlaps)
}

#################################################################################################################################################
# FUNCTIONS FOR CALCULATING GradCAM LOCALIZATION RESULTS FROM INDIVIDUAL RBP BINDING MODELS, AND LOADING THEM INTO ONE TENSOR.
#################################################################################################################################################

# Calculate and store GradCAM heatmap results for all RBPs on the given input sequences tensor
calculate_gradcams <- function(input_tensor, work_folder, k_reset_freq=3, motif_length=8, motif_min_percent_distinct=0.5, motif_score_min_cutoff=0.1) {
    dir.create(work_folder, showWarnings = FALSE)
    tf$compat$v1$disable_eager_execution()
    rbps <- get_features_by_group("RBP")
    rbp_heatmaps <- lapply(1:length(rbps), function(rbp_i) {
        if(rbp_i %% k_reset_freq == 0) { k_clear_session(); gc() }
        rbp = rbps[rbp_i]
        print(paste0(rbp_i,". ",rbp))
        rbp_filename = full_path(work_folder,paste0(rbp,"_gnomad_gradcam.rds"))
        if(file.exists(rbp_filename)) { return(0) }#read.csv(rbp_filename, row.names=1)) }
        model <- load_model_hdf5(paste0("../ML/output/",tolower(rbp),"_model2.h5"))
        #input_tensor <- dat_tensor[,,,][1,,,drop=FALSE]
        #model %>% predict(input_tensor)
        rbp_output <- model$output[,1] 
        last_conv_layer <- model$get_layer("conv1") 
        grads <- k_gradients(rbp_output, last_conv_layer$output)[[1]]
        if(is.null(grads)) { return(NULL) }
        pooled_grads <- k_mean(grads, axis=c(1,2))
        iterate <- k_function(list(model$input), list(pooled_grads, last_conv_layer$output[1,,]))
        
        #sequence_to_display = test_indices[sequence_index_to_display]
        seq_heatmaps <- abind(lapply(1:nrow(input_tensor), function(sequence_index) {
            #print(sequence_index)
            iterate_result <- iterate(input_tensor[sequence_index,,,drop=FALSE])
            pooled_grads_value <- iterate_result[[1]]; conv_layer_output_value <- iterate_result[[2]]
            for (i in 1:100) {
                conv_layer_output_value[,i] <- conv_layer_output_value[,i] * pooled_grads_value[i]
            }
            # The channel-wise mean of the resulting feature map is our heatmap of class activation
            heatmap <- smooth.spline(spline(apply(conv_layer_output_value, 1, mean), n=ncol(input_tensor), method="fmm")$y)$y
            # Normalize the heatmap between 0 and 1, for better visualization.
            heatmap <- pmax(heatmap, 0)
            heatmap <- heatmap / max(heatmap)
            
            #################################################################################################################################################
            # CODE FOR DRAWING SEQUENCE HEATMAPS LIKE I SHOW IN PPT!
            #################################################################################################################################################
            # sequence <- strsplit(tensor_to_genomic_sequences(data[sequence_to_display,,,]), "")[[1]]
            # sequence_length <- length(sequence)
            # sequence_matrix <- suppressWarnings(matrix(sequence, ncol=ceiling(sqrt(ncol(data))), byrow=TRUE))
            # if((ncol(sequence_matrix) * nrow(sequence_matrix)) > sequence_length) { sequence_matrix[nrow(sequence_matrix), seq((sequence_length %% ncol(sequence_matrix))+1, ncol(sequence_matrix))] <- "" }
            # rownames(sequence_matrix) <- paste0(seq(1, sequence_length, by=ncol(sequence_matrix)),"...",sapply(seq(ncol(sequence_matrix), sequence_length+ncol(sequence_matrix)-1, by=ncol(sequence_matrix)), function(x) min(c(x, sequence_length))))
            # library(plotrix)
            # filename=output_path(paste0(model_name,"_seq",sequence_to_display,"_rbp_binding_heatmap.pdf"))
            # pdf(file=filename)
            # plot(-1, -1, xlim=c(0,1), ylim=c(0,1), main=paste0(colnames(labels)[labels[sequence_to_display,]==1],paste0(" seq",sequence_to_display),", pred_score=",round(pred_scores[sequence_index_to_display],3)), xlab="", ylab="", xaxt="n", yaxt="n", axes=FALSE, cex.main=1.6, cex.lab=1.5) #"RBP binding strength heatmap"
            # #mtext(paste0(colnames(labels)[labels[sequence_to_display,]==1]," sequence, pred_score=",round(pred_scores[sequence_index_to_display],3)))
            # sequence_matrix_cols <- matrix(c(adjustcolor(rgb(colorRamp(c("blue","white","red"))(heatmap[1:sequence_length]), maxColorValue=255), alpha.f = 0.8), rep("white", (ncol(sequence_matrix) * nrow(sequence_matrix)) - sequence_length)), ncol=ceiling(sqrt(ncol(data))), byrow=TRUE)
            # addtable2plot("topleft",table=sequence_matrix, bg=sequence_matrix_cols, display.rownames=TRUE, display.colnames=FALSE, cex=2.2)
            # 
            # # Determine the optimal bound region of fixed motif_length from the heatmap 
            # motif_start_index = which.max(rollapply(heatmap, width=motif_length, mean))
            # motif_end_index = motif_start_index + motif_length - 1
            # trimmed_binding_site <- paste0(sequence[motif_start_index:motif_end_index], collapse="")
            # text(0.5, 0.01, paste0("Optimal binding site: ",trimmed_binding_site), cex=1.7, font=2)
            # 
            # dev.off()
            # pdf_to_png(filename)
            # 
            # heatmap_smoothed <- rollapply(heatmap, width=motif_length, mean, align="left", partial=TRUE)
            # motif_start_indices <- which(heatmap_smoothed >= motif_score_min_cutoff)
            # motif_start_indices <- motif_start_indices[order(heatmap_smoothed[motif_start_indices], decreasing=TRUE)]
            # get_best_nonoverlapping <- function(motif_start_indices) { if(length(motif_start_indices) > 0) { return(c(motif_start_indices[1], get_best_nonoverlapping(motif_start_indices[abs(motif_start_indices - motif_start_indices[1]) >= ceiling(motif_length*motif_min_percent_distinct)]))) } else { return(c()) } }
            # motif_start_indices <- get_best_nonoverlapping(motif_start_indices)
            # motif_end_indices <- motif_start_indices + motif_length - 1
            # trimmed_binding_sites <- sapply(1:length(motif_start_indices), function(motif_i) paste0(sequence[motif_start_indices[motif_i]:motif_end_indices[motif_i]], collapse=""))
            # 
            # motifs_dat <- unfactorize(data.frame(rep(paste0("seq",sequence_to_display),length(motif_start_indices)), rep(model_name,length(motif_start_indices)), rep(pred_scores[sequence_index_to_display],length(motif_start_indices)), motif_start_indices, motif_end_indices, trimmed_binding_sites, heatmap_smoothed[motif_start_indices], heatmap_smoothed[motif_start_indices]/heatmap_smoothed[motif_start_indices][1]))
            # colnames(motifs_dat) <- c("seq", "model", "pred_score", "start", "end", "motif", "score", "score_norm")
            # return(motifs_dat)
            
            return(heatmap)
        }), along=2)
        saveRDS(seq_heatmaps, rbp_filename)
        #write.csv(seq_heatmaps, file=rbp_filename)
        rm(seq_heatmaps)
        return(0)
    })
    return(0)
}

# Get calculated GradCAMS from intermediate results folder, and stack them into one Tensor object.
get_gradcams <- function(work_folder, region_indices=NULL, rbp_pat=NULL, dim=NULL) {
    rbps <- get_features_by_group("RBP")
    subset_regions = !is.null(region_indices)
    subset_rbps = !is.null(rbp_pat) && !is.null(dim)
    rbp_heatmaps <- abind(lapply(1:length(rbps), function(rbp_i) {
        rbp = rbps[rbp_i]
        print(paste0(rbp_i,". ",rbp))
        if(subset_rbps && !grepl(rbp_pat,rbp)) { return(array(0,dim=dim)) }
        rbp_filename = full_path(work_folder,paste0(rbp,"_gnomad_gradcam.rds"))
        if(file.exists(rbp_filename)) {
            if(subset_regions) { return(readRDS(rbp_filename)[,region_indices]) 
            } else { return(readRDS(rbp_filename)) }
            #if(subset_regions) { return(read.csv(rbp_filename, row.names=1)[,region_indices]) 
            #} else { return(read.csv(rbp_filename, row.names=1)) }
        } else { return(NULL) }
    }), along=3)
    rbp_heatmaps <- aperm(rbp_heatmaps, c(2,1,3))
    return(rbp_heatmaps)
}

# Get the subset of the full GradCAM tensor, defined by the specified indices groups; not used in current workflow, but kept in code for now.
# For example:> split_gradcams("asd", indices_groups=seq(40001, 126700, by=28900))
split_gradcams <- function(full_dat_name, indices_groups) {
    full_dat <- readRDS(output_path(paste0(full_dat_name,"_dat_all.rds")))
    indices_groups <- rollapply(unique(c(indices_groups, nrow(full_dat))), width=2, FUN=function(x) { return(c(x[1],x[2])) })
    print(indices_groups)
    subdat_names <- sapply(1:nrow(indices_groups), function(i) paste0(full_dat_name,i))
    for(i in 1:length(subdat_names)) {
        print(paste0("Writing ",subdat_names[i],"_dat.rds regions subset..."))
        subdat_indices <- indices_groups[i,1]:indices_groups[i,2]
        saveRDS(full_dat[subdat_indices,], output_path(paste0(subdat_names[i],"_dat.rds")))
    }
    work_folders <- sapply(subdat_names, function(subdat_name) { work_folder = output_path(paste0(subdat_name,"_gradcams")); dir.create(work_folder, showWarnings = FALSE); return(work_folder) })
    full_gradcams_folder = output_path(paste0(full_dat_name,"_gradcams"))
    full_gradcams_files <- list.files(full_gradcams_folder)
    for(gradcam_filename in full_gradcams_files) {
        print(paste0("Loading ",gradcam_filename,"..."))
        rbp_full_gradcam <- read.csv(full_path(full_gradcams_folder, gradcam_filename), row.names=1)
        for(i in 1:length(subdat_names)) {
            subdat_indices <- indices_groups[i,1]:indices_groups[i,2]
            print(paste0("Writing GradCAMs subset to ",work_folders[i],"..."))
            write.csv(rbp_full_gradcam[,subdat_indices], file=full_path(work_folders[i],gradcam_filename))
        }
    }
    return(0)
}

#################################################################################################################
# MAIN: FUNCTIONS TO SETUP SUPERMODEL DATA TENSORS WITH A GIVEN NAME
#################################################################################################################

# Turn a regional type of tensor into four different tensors representing each variant-specific possibility for vSUPRNOVA.
# convert_tensor_to_variant_forms(paste0(c("transcribed100k_"),1:4))
convert_tensor_to_variant_forms <- function(dat_names, make_gradcams=TRUE) {
    alts <- c("A","C","G","T")
    for(dat_name in dat_names) {
        print(dat_name)
        ref_tensor_filename = output_path(paste0(dat_name,"_tensor.rds"))
        ref_tensor <- readRDS(ref_tensor_filename)
        for(alt_i in 1:4) {
            alt <- alts[alt_i]; print(alt)
            alt_tensor <- ref_tensor
            alt_tensor[,ceiling(ncol(alt_tensor)/2),alt_i,1] <- 1; alt_tensor[,ceiling(ncol(alt_tensor)/2),-c(alt_i),1] <- 0
            alt_tensor_filename = gsub("_tensor.rds$", paste0("_",alt,"_tensor.rds"), ref_tensor_filename)
            saveRDS(alt_tensor, alt_tensor_filename)
        }
    }
    
    if(make_gradcams) {
        dat_names <- c(sapply(dat_names, FUN=function(x) paste0(x,"_",alts)))
        for(dat_name in dat_names) {
            dat_gradcams_filename = output_path(paste0(dat_name,"_gradcams.rds"))
            if(file.exists(dat_gradcams_filename)) { next }
            
            dat_gradcams_folder= output_path(paste0(dat_name,"_gradcams"))
            print(paste0("Calculating GradCAMS for ",dat_gradcams_folder," folder..."))
            dat_tensor_filename = output_path(paste0(dat_name,"_tensor.rds"))
            dat_tensor <- readRDS(dat_tensor_filename)
            calculate_gradcams(dat_tensor[,,,], dat_gradcams_folder, k_reset_freq=3)
            gc()
            
            dat_gradcams <- get_gradcams(dat_gradcams_folder)
            saveRDS(dat_gradcams, dat_gradcams_filename)
            dat_gradcams <- array_reshape(dat_gradcams, c(dim(dat_gradcams),1))
            dat_gradcams[is.na(dat_gradcams)] <- 0
            saveRDS(dat_gradcams, dat_gradcams_filename)
        }
    }
}

# Setup and save all data tensors for SUPRNOVA under the given name. 
# Run this function multiple times (with specified number of tasks to complete) to ensure that all setup is complete; this is designed to manage RAM for large datasets.
setup_supermodel_data <- function(dat_name, dat=NULL, version="hg19", width=151, num_tasks_to_do=1, rewrite=FALSE) {
    arm_width = floor(width/2)
    dat_regions_filename = output_path(paste0(dat_name,"_dat.rds"))
    dat_granges_hg19_filename = output_path(paste0(dat_name,"_granges_hg19.rds"))
    dat_granges_hg38_filename = output_path(paste0(dat_name,"_granges_hg38.rds"))
    dat_trimers_filename = output_path(paste0(dat_name,"_trimers.rds"))
    dat_cpg_filename = output_path(paste0(dat_name,"_cpg.rds"))
    dat_tensor_filename = output_path(paste0(dat_name,"_tensor.rds"))
    dat_gradcams_folder= output_path(paste0(dat_name,"_gradcams"))
    dat_gradcams_filename = output_path(paste0(dat_name,"_gradcams.rds"))
    #dat_pLIs_filename = output_path(paste0(dat_name,"_pLIs.rds"))
    dat_pLIs_filename = output_path(paste0(dat_name,"_obs_exp.rds"))
    dat_expected_filename = output_path(paste0(dat_name,"_expected.rds"))
    
    tasks_done = 0
    # Processed data / genomic regions
    if(tasks_done < num_tasks_to_do && (!file.exists(dat_regions_filename) || rewrite)) {
        print(paste0("Setting up ",dat_regions_filename,"..."))
        dat <- standardize_colnames(data.frame(dat), remove_chr_prefix=TRUE, re_order=TRUE)
        dat <- cbind(paste0("seq",1:nrow(dat)), dat[,c(1:2,5,3:4)], dat[,2]-arm_width, dat[,2]+arm_width)
        colnames(dat) <- c("seq", "Chrom", "Position", "sample", "Ref", "Alt", "start", "end")
        dat <- dat[dat$Chrom %in% c(1:22),]
        dat <- add_sequence_context_feature(dat, version="hg19", width=width)
        colnames(dat)[which(colnames(dat) %in% c("Position", "ref_sequence"))] <- c("Position_hg19", "ref_sequence_hg19")
        dat <- liftover(dat, from="hg19", to="hg38", chr_colname="Chrom", start_colname="start", end_colname="end", ref_colname="Ref", alt_colname="Alt", confirm_refseq=FALSE, mismatches_pause=FALSE)#[,c("Chrom","Position","Ref","Alt","sample","snv_indel","Chrom_hg38","Position_hg38")]
        dat <- dat[dat$Chrom %in% c(1:22) & (dat$end - dat$start + 1 == width),]
        dat <- cbind(dat, dat$start + arm_width); colnames(dat)[ncol(dat)] <- "Position"
        dat <- add_sequence_context_feature(dat, version="hg38", width=width)
        dat <- dat[dat$ref_sequence == dat$ref_sequence_hg19, c("seq", "Chrom", "Position", "start", "end", "sample", "Ref", "Alt", "ref_sequence", "Chrom_hg19", "Position_hg19", "start_hg19", "end_hg19")]
        saveRDS(dat, dat_regions_filename)
        
        print(paste0("Setting up ",dat_granges_hg19_filename,"..."))
        dat_granges_hg19 <- to_genomic_regions(dat, chr_colname="Chrom_hg19", start_colname="start_hg19", end_colname="end_hg19", label_colname="seq")
        dat_granges_hg19 <- dat_granges_hg19[order(as.numeric(gsub("seq","",names(dat_granges_hg19))))]
        saveRDS(dat_granges_hg19, dat_granges_hg19_filename)
        print(paste0("Setting up ",dat_granges_hg38_filename,"..."))
        dat_granges_hg38 <- to_genomic_regions(dat, chr_colname="Chrom", start_colname="start", end_colname="end", label_colname="seq")
        dat_granges_hg38 <- dat_granges_hg38[order(as.numeric(gsub("seq","",names(dat_granges_hg38))))]
        saveRDS(dat_granges_hg38, dat_granges_hg38_filename)
        tasks_done = tasks_done + 1
    } else if(tasks_done < num_tasks_to_do) { print(paste0(dat_regions_filename," already defined!")) }
    # Trimers and CpG annotation
    if(tasks_done < num_tasks_to_do && (!file.exists(dat_cpg_filename) || rewrite)) {
        print(paste0("Setting up ",dat_trimers_filename,"..."))
        dat <- data.frame(readRDS(dat_regions_filename))
        dat_trimers <- rbindlist(lapply(1:nrow(dat), function(i) {
            if(i %% 1000 == 1) { print(i) }
            ref_seq <- paste0(dat$ref_sequence[i])
            return(data.frame(t(data.frame(rollapply(strsplit(ref_seq,"")[[1]], width=3, by=1, partial=2, FUN=function(x) paste0(x,collapse=""))))))
        }))
        saveRDS(dat_trimers, dat_trimers_filename)
        print(paste0("Setting up ",dat_cpg_filename,"..."))
        dat_cpg <- t(apply(dat_trimers, 1, function(x) grepl("CG", x)))
        saveRDS(dat_cpg, dat_cpg_filename)
        tasks_done = tasks_done + 1
    } else if(tasks_done < num_tasks_to_do) { print(paste0(dat_cpg_filename," already defined!")) }
    # Genomic tensor object
    if(tasks_done < num_tasks_to_do && (!file.exists(dat_tensor_filename) || rewrite)) {
        print(paste0("Setting up ",dat_tensor_filename,"..."))
        dat <- readRDS(dat_regions_filename)
        dat_tensor <- genomic_sequences_to_tensor(dat$ref_sequence, sequence_length=151)
        saveRDS(dat_tensor, dat_tensor_filename)
        tasks_done = tasks_done + 1
    } else if(tasks_done < num_tasks_to_do) { print(paste0(dat_tensor_filename," already defined!")) }
    # Calculate and stack GradCAMs
    if(tasks_done < num_tasks_to_do && (!file.exists(dat_gradcams_filename) || rewrite)) {
        if(length(list.files(dat_gradcams_folder)) < 160) {
            print(paste0("Calculating GradCAMS for ",dat_gradcams_folder," folder..."))
            dat_tensor <- readRDS(dat_tensor_filename)
            calculate_gradcams(dat_tensor[,,,], dat_gradcams_folder)
            tasks_done = tasks_done + 1
        } else { print(paste0("GradCAMS in ",dat_gradcams_folder," already all calculated!")) }
        if(tasks_done < num_tasks_to_do) {
            print(paste0("Setting up ",dat_gradcams_filename,"..."))
            dat_gradcams <- get_gradcams(dat_gradcams_folder) #, region_indices=(dat_start:dat_end))
            saveRDS(dat_gradcams, dat_gradcams_filename)
            dat_gradcams <- array_reshape(dat_gradcams, c(dim(dat_gradcams),1))
            dat_gradcams[is.na(dat_gradcams)] <- 0
            saveRDS(dat_gradcams, dat_gradcams_filename)
            tasks_done = tasks_done + 1
        }
    } else if(tasks_done < num_tasks_to_do) { print(paste0(dat_gradcams_filename," already defined!")) }
    # ExAC pLI
    if(tasks_done < num_tasks_to_do && (!file.exists(dat_pLIs_filename) || rewrite)) {
        print(paste0("Setting up ",dat_pLIs_filename,"..."))
        dat_granges_hg19 <- readRDS(dat_granges_hg19_filename)
        genebody_hg19 <- read.csv(data_path("refGene\\refGene_hg19_genebody_fixed.bed"), sep="\t", stringsAsFactors=FALSE)
        genebody_granges_hg19 <- to_genomic_regions(genebody_hg19, chr_colname="chromosome", start_colname="tss", end_colname="tes", strand_colname="strand", label_colname="gene", order_coordinates=TRUE, remove_duplicate_labels=TRUE)
        #exac_dat <- read.csv(data_path("fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt"), sep="\t")
        genes <- readRDS(output_path("genes.rds"))
        nearest_genes_dat <- names(genebody_granges_hg19)[nearest(dat_granges_hg19, genebody_granges_hg19)]
        #pLIs_dat <- exac_dat[unlist(sapply(nearest_genes_dat,function(x) { nearest_gene <- which(paste0(exac_dat$gene) == x); if(length(nearest_gene) < 1) { return(NA) } else { return(nearest_gene[1]) }  })),"pLI"]
        pLIs_dat <- genes[unlist(sapply(nearest_genes_dat,function(x) { nearest_gene <- which(paste0(exac_dat$gene) == x); if(length(nearest_gene) < 1) { return(NA) } else { return(nearest_gene[1]) }  })),"o_e_Z"]
        pLIs_dat[is.na(pLIs_dat)] <- mean(pLIs_dat[!is.na(pLIs_dat)])
        pLIs_dat <- array(pLIs_dat, c(length(pLIs_dat),1))
        saveRDS(pLIs_dat, dat_pLIs_filename)
        tasks_done = tasks_done + 1
    } else if(tasks_done < num_tasks_to_do) { print(paste0(dat_pLIs_filename," already defined!")) }
    # Position-specific expected (background mutation rate)
    if(tasks_done < num_tasks_to_do && (!file.exists(dat_expected_filename) || rewrite)) {
        print(paste0("Setting up ",dat_expected_filename,"..."))
        dat <- readRDS(dat_regions_filename)
        dat_expected <- get_expected_muts(dat, chr_colname="Chrom_hg19", pos_colname="Position_hg19", width=151, version="hg19")
        #dat_expected <- t(array(dat_expected, c(151,dim(dat)[1])))
        saveRDS(dat_expected, dat_expected_filename)
        tasks_done = tasks_done + 1
    } else if(tasks_done < num_tasks_to_do) { print(paste0(dat_expected_filename," already defined!")) }
    
    print(paste0(tasks_done," tasks completed!"))
    if(tasks_done > 0) { print("Exiting for now. Reboot R and re-run setup function again.") } else { print("Setup already complete!") }
}

# Combine supermodel data that was processed and setup in split batches.
# An example with ASD:> setup_supermodel_data("asd1"); setup_supermodel_data("asd2"); setup_supermodel_data("asd3"); combine_supermodel_data(dat_names=c("asd","asd1","asd2","asd3"), combined_dat_name="asd_all")
combine_supermodel_data <- function(dat_names, combined_dat_name) {
    dat_regions_filenames = output_path(paste0(c(combined_dat_name,dat_names),"_dat.rds"))
    dat_granges_hg19_filenames = output_path(paste0(c(combined_dat_name,dat_names),"_granges_hg19.rds"))
    dat_granges_hg38_filenames = output_path(paste0(c(combined_dat_name,dat_names),"_granges_hg38.rds"))
    dat_trimers_filenames = output_path(paste0(c(combined_dat_name,dat_names),"_trimers.rds"))
    dat_cpg_filenames = output_path(paste0(c(combined_dat_name,dat_names),"_cpg.rds"))
    dat_tensor_filenames = output_path(paste0(c(combined_dat_name,dat_names),"_tensor.rds"))
    #dat_gradcams_folder = output_path(paste0(c(combined_dat_name,dat_names),"_gradcams"))
    dat_gradcams_filenames = output_path(paste0(c(combined_dat_name,dat_names),"_gradcams.rds"))
    dat_pLIs_filenames = output_path(paste0(c(combined_dat_name,dat_names),"_pLIs.rds"))
    dat_expected_filenames = output_path(paste0(c(combined_dat_name,dat_names),"_expected.rds"))
    
    #saveRDS(rbindlist(lapply(dat_regions_filenames[-1], function(f) readRDS(f))), dat_regions_filenames[1])
    saveRDS(Reduce(c, lapply(dat_granges_hg19_filenames[-1], function(f) readRDS(f))), dat_granges_hg19_filenames[1])
    saveRDS(Reduce(c, lapply(dat_granges_hg38_filenames[-1], function(f) readRDS(f))), dat_granges_hg38_filenames[1])
    saveRDS(abind(lapply(dat_trimers_filenames[-1], function(f) readRDS(f)), along=1), dat_trimers_filenames[1])
    saveRDS(abind(lapply(dat_cpg_filenames[-1], function(f) readRDS(f)), along=1), dat_cpg_filenames[1])
    #saveRDS(abind(lapply(dat_tensor_filenames[-1], function(f) readRDS(f)), along=1), dat_tensor_filenames[1])
    saveRDS(abind(lapply(dat_gradcams_filenames[-1], function(f) readRDS(f)), along=1), dat_gradcams_filenames[1])
}

# Combine independently trained models, such as for different RBPs, into a single model with shared input that can be loaded and run much faster.
combine_models <- function(models, combined_model_name="combined_model") {
    all_models <- lapply(models, function(model_path) {
        print(paste0("Loading model from ",model_path))
        model_path <- paste0("../ML/output/",tolower("K562.RBFOX2"),"_model2.h5")
        model_name <- gsub("^.*\\/([a-zA-Z0-9\\.]*)_model.h5*","\\1", model_path)
        model <- load_model_hdf5(model_path)
        #for(layer in model$layers) { print(layer$name) } #layer$name <- paste0(model_name,"_",layer$name) }
        #sapply_out <- sapply(model$layers, function(layer) paste0(layer$name))
        return(model)
    })
    shared_input <- all_models[[1]]$input
    
    all_model_sections <- lapply(all_models, function(model) {
        conv1 <- get_layer(model, "conv1")
        return(shared_input %>% model)
    })
    combined_model <- keras_model(inputs=c(shared_input), outputs=c(all_model_sections))
    
    opt <- optimizer_sgd(lr = 0.1, decay = 1e-2, momentum=0.5, nesterov=TRUE)
    masked_loss_function <- function(y_true, y_pred, mask=-1) { mask_vector <- k_cast(k_not_equal(y_true, mask), k_floatx());
    return(k_binary_crossentropy(y_true * mask_vector, y_pred * mask_vector)) }
    #model %>% compile(loss = masked_loss_function, optimizer = opt, metrics = "accuracy")
    combined_model %>% compile(loss = "binary_crossentropy", optimizer = opt, metrics = "accuracy") #loss = "categorical_crossentropy"
    
    # Draw model network
    #plot_model(combined_model, to_file = output_path(paste0(combined_model_name,".png")), show_shapes = TRUE, show_layer_names = TRUE)
    #plot_model(combined_model, to_file = output_path(paste0("images/",combined_model_name,".png")), show_shapes = TRUE, show_layer_names = TRUE)
    
    combined_model %>% save_model_hdf5(output_path(paste0(combined_model_name,".h5")))
    return(combined_model)
}
#combined_model <- combine_models(models=paste0("../ML/output/",tolower(get_features_by_group("RBP")),"_model.h5"))

#################################################################################################################
# FUNCTIONS TO DO SEQUENCE ALIGNMENT ON DATA AND ANALYZE PWMs 
#################################################################################################################

# Return multiple sequence alignment of the given sequences, with gaps filled in with the correct neighboring sequences.
get_aligned_sequences <- function(gr, min_sequence_length=NULL, max_sequence_length=NULL, bucket_size=1000, bucket=1, shuffle=FALSE, filename=NULL, method="ClustalOmega") {
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
        
        sequences_msa <- msa(sequences, order="input", method=method)
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
        
        m_consensus <- msaConsensusSequence(sequences_msa)
        print(m_consensus)
        
        if(!is.null(filename)) {
            gr_subset <- gr[bucket_start_index:bucket_end_index] 
            log_table <- data.frame(paste0("chr",seqnames(gr_subset)), start(gr_subset), end(gr_subset), width(gr_subset), strand(gr_subset), names(gr_subset), paste0(m), paste0(m_filled), width(m_filled))
            colnames(log_table) <-  c("chrom", "start", "end", "width", "strand", "p.value", "msa", "msa_filled", "msa_width")
            has_p.value <- sum(is.na(as.numeric(paste0(log_table$p.value)))) == 0
            if(! has_p.value) { log_table <- log_table[,-which(colnames(log_table) == "p.value")] }
            write.csv(log_table, file=filename, row.names=FALSE)
        }
        
        return(m_filled)
    })[[1]]
    
    if(!is.null(filename)) {
        
        
    }
    
    return(gr_sequences)
}

# Analyze PWM, either given directly or derived from sequences, and draw sequence logo
pwm_analysis <- function(rbp=NULL, rbp_granges=NULL, rbp_pwm=NULL, rbp_sequences=NULL, motif_length=7, min_sequence_length=NULL, max_sequence_length=NULL, bucket_size=1000, bucket=1, composite=FALSE) {
    if(is.null(rbp_pwm)) {
        if(is.null(rbp_sequences)) {
            if(is.null(rbp_granges)) {
                all_rbp_features <- get_features_by_group("RBP_0bp")
                rbp_features <- all_rbp_features[grepl(rbp, all_rbp_features)]
                rbp_granges <- sapply(rbp_features, function(rbp_feature) { load_annotation(rbp_feature) })
                rbp_granges <- Reduce(c, rbp_granges)
            }
            names_numeric <- as.numeric(names(rbp_granges))
            has_p.value <- sum(is.na(names_numeric)) == 0
            if(has_p.value) { rbp_granges <- rbp_granges[order(names_numeric)] } # order by p.value/significance of peak
            
            rbp_sequences <- get_aligned_sequences(rbp_granges, min_sequence_length=min_sequence_length, max_sequence_length=max_sequence_length, bucket_size=bucket_size, bucket=bucket, shuffle=FALSE, filename=output_path(paste0(rbp,"_eclip_msa",bucket,".csv")))
            rbp_sequences <- rbp_sequences[!grepl("[^ATCG]", rbp_sequences)]
        }
        if(composite) { rbp_sequences <- sequence_composite(rbp_sequences) }
        rbp_pwm <- TFBSTools::toPWM(rbp_sequences, type="prob", pseudocounts=0.8,  bg=genomic.acgt)
    }
    rbp_seqlogo_pwm <- makePWM(rbp_pwm)
    rbp_seqlogo_pwm@consensus; rbp_seqlogo_pwm@ic
    
    motif_length <- min(c(motif_length, ncol(rbp_pwm)))
    motif_start_index = which.max(rollapply(rbp_seqlogo_pwm@ic, width=motif_length, mean))
    motif_end_index = motif_start_index + motif_length - 1
    trimmed_rbp_pwm <- rbp_pwm[,motif_start_index:motif_end_index]
    
    seqlogo_title <- paste0(rbp," PWM logo"); if(composite) { seqlogo_title <- paste0(seqlogo_title," composite") }
    
    seqLogo(trimmed_rbp_pwm, ic.scale=TRUE)
    grid.text(seqlogo_title, x=0.5, y=0.95, gp=gpar(fontsize=18, fontface="bold", col="black"))
    filename = output_path(paste0(gsub(" ", "_", tolower(seqlogo_title)),bucket,".pdf"))
    dev.copy2pdf(file=filename)
    pdf_to_png(filename)
    
    seqlogo_title = gsub("logo", "extended logo", seqlogo_title)
    seqLogo(rbp_pwm, ic.scale=TRUE)
    grid.text(seqlogo_title, x=0.5, y=0.95, gp=gpar(fontsize=18, fontface="bold", col="black"))
    filename = output_path(paste0(gsub(" ", "_", tolower(seqlogo_title)),bucket,".pdf"))
    dev.copy2pdf(file=filename)
    pdf_to_png(filename)
    
    return(trimmed_rbp_pwm)
}

#################################################################################################################
# OTHER FUNCTIONS
#################################################################################################################

# Generate construct random sequences from bases; note that this is very different from randomly selecting sequences from genome.
construct_random_sequences <- function(sequence_length, num_sequences) {
    random_sequences <- sapply(1:num_sequences, function(i) paste0(sample(c("A","C","G","T"), sequence_length, replace=TRUE), collapse=""))
    return(random_sequences)
}

# Make sure that the start and end coordinates in data are properly ordered. 
# //TODO: Modify this to take strand into account, which "flips" which of start/end should be smaller.
order_start_end_coordinates <- function(dat, start_colname = "start", end_colname = "end") {
    badly_ordered <- dat[,start_colname] > dat[,end_colname]
    if(sum(badly_ordered) > 0) {
        badly_ordered_starts <- dat[badly_ordered, start_colname]
        dat[badly_ordered, start_colname] <- dat[badly_ordered, end_colname]
        dat[badly_ordered, end_colname] <- badly_ordered_starts
    }
    return(dat)
}

# Run variant threshold test to get real p.values
# This is an older function, but left in the the code for now.
variant_threshold_test <- function(variant_types = c("SNV", "indel"), variant_constraints = c("autism_gene", "constrained_gene", "H3K36me3", ""), thresholds = seq(0.5, 1, by=0.01)) {
    sapply_out <- sapply(1:160, function(i) { #1:ncol(delta_pred_scores)
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
                    
                    if(FALSE) {
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
                    
                    if(FALSE) { return(unlist(c(greater_than_threshold_fets[which.min(greater_than_threshold_fets$p.value),], less_than_threshold_fets[which.min(less_than_threshold_fets$p.value),]))) 
                    } else { return(unlist(greater_than_threshold_fets[which.min(greater_than_threshold_fets$p.value),])) }
                }))
                
                variable_threshold_results_lof_entry <- c(colnames(delta_pred_scores)[i], variant_type, variant_constraint, variable_threshold_results[1,1:9], (sum(variable_threshold_results[-1,5] <= variable_threshold_results[1,5])+1)/(num_permutations+1))
                names(variable_threshold_results_lof_entry)[c(1:3,13)] <- c("RBP", "variant_type", "variant_constraint", "real_p.value")
                variable_threshold_results_lof <- rbind(unfactorize(variable_threshold_results_lof), unfactorize(variable_threshold_results_lof_entry))
                
                if(FALSE) { 
                    variable_threshold_results_gof_entry <- c(colnames(delta_pred_scores)[i], variant_type, variant_constraint, variable_threshold_results[1,10:18], (sum(variable_threshold_results[-1,14] <= variable_threshold_results[1,14])+1)/(num_permutations+1))
                    names(variable_threshold_results_gof_entry)[c(1:3,13)] <- c("RBP", "variant_type", "variant_constraint", "real_p.value")
                    variable_threshold_results_gof <- rbind(unfactorize(variable_threshold_results_gof), unfactorize(variable_threshold_results_gof_entry))
                }
            }
        }
        colnames(variable_threshold_results_lof) <- c("RBP", "variant_type", "variant_constraint", "threshold", "estimate", "conf.int_lower", "conf.int_higher", "p.value", "m1", "m0", "n1", "n0", "real_p.value")
        write.csv(variable_threshold_results_lof, file=output_path(paste0(gsub(" ","_",rbp),"_ref_prediction_score_variable_threshold_enrichment.csv")), row.names=FALSE)
        if(FALSE) {
            colnames(variable_threshold_results_gof) <- c("RBP", "variant_type", "variant_constraint", "threshold", "estimate", "conf.int_lower", "conf.int_higher", "p.value", "m1", "m0", "n1", "n0", "real_p.value")
            write.csv(variable_threshold_results_gof, file=output_path(paste0(gsub(" ","_",rbp),"_ref_prediction_score_variable_threshold_gof.csv")), row.names=FALSE)
        }
    })
    
    legend("topright", legend=c(gsub("^.*\\.", "", colnames(delta_pred_scores)), "Positives", "Negatives"), col=c(cols,"black","black"), pch=c(rep(15,ncol(delta_pred_scores)),NA,NA), lty=c(rep(NA,ncol(delta_pred_scores)),1,3))
    #mtext(paste0("Validation set of ",length(test_indices)," length-",sequence_length," sequences"))
    dev.off()
    pdf_to_png(filename)
}










