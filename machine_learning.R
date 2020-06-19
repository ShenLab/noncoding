library(keras)
library(kerasR)
library(stringr)
library(pbapply)
library(EBImage)
library(rhdf5)
library(gtools)
library("PWMEnrich")
library("TFBSTools")
library("seqLogo")
library("msa")
library("ComplexHeatmap")
library("abind")
library("hexbin")
library("RColorBrewer")
library("ggplot2")
library("tfprobability")
library("tidyverse")
source("../WGS/alex_suite.R")

#a <- read.csv(data_path("WGS/An_2018_Table_S2_denovos.csv"), skip=1, header=TRUE, sep="\t") #https://www.ncbi.nlm.nih.gov/pubmed/30545852

#################################################################################################################
# PARAMETERS
#################################################################################################################
# Delimiter between folders/files in filepaths. Should be "/" for Linux and "\\" for Windows.
FILEPATH_DELIM <- "/"
# Path to data folder, allowing simpler loading of data with the alex_suite.R data_path function.
DATA_FOLDER <- "/home/local/ARCS/ak3792/Documents/Research/data"
# Path to big data folder, allowing simpler loading of big data with the alex_suite.R bigdata_path function.
BIGDATA_FOLDER <- "/data"
# Path to output folder, allowing simpler writing of output with the alex_suite.R output_path function.
OUTPUT_FOLDER <- "/home/local/ARCS/ak3792/Documents/Research/ML/output"

#################################################################################################################
# Load regulatory variants (HGMD + gnomAD) data (processed in ../WGS/nc_analysis.R)
#################################################################################################################
regulatory_variant_dat <- read.csv(file="../WGS/output/regulatory_variants.tsv", sep="\t") 
regulatory_variant_dat <- standardize_colnames(regulatory_variant_dat)
#asd_variant_dat <- read.csv(data_path("ASD_OlgaT.tsv"), sep="\t") # a <- read.table(file="/home/local/ARCS/ak3792/Documents/Research/WGS/output/regulatory_variants.tsv", sep="\t", header=1) # For annotation with WGSA
asd_variant_dat <- read.csv(data_path("WGS/An_2018_Table_S2_denovos.csv"), skip=1, header=TRUE, sep="\t")
asd_variant_dat <- standardize_colnames(asd_variant_dat)
chd_variant_dat <- read.csv(data_path("WGS/CHDFB/chdfb_sscfb_pcgc_final_variants.csv"))[,-1]
chd_variant_dat <- standardize_colnames(chd_variant_dat, re_order=TRUE)

# Add sequence context features:
add_sequence_context_feature <- function(dat, width=51, version="hg19", chr_colname="Chrom", pos_colname="Position", ref_colname="Ref", alt_colname=NULL, split_multiple_alts=TRUE) {
    #dat <- dat[,-which(colnames(dat) %in% c("ref_sequence","alt_sequence"))]
    arm_width = floor(width/2)
    if(is.null(alt_colname)) {
        skip_alts = TRUE; all_alt_combinations = FALSE
    } else if (alt_colname == "all") {
        skip_alts = FALSE; all_alt_combinations = TRUE
    } else {
        skip_alts = !(alt_colname %in% colnames(dat)); all_alt_combinations = FALSE
    }
    if(!skip_alts && split_multiple_alts && !all_alt_combinations && sum(grepl(",", dat[,alt_colname]) > 0) & length(width)==1) { 
        cat("Splitting multiple alternate sequences...")
        library("tidyr")
        alt_colname_index = which(colnames(dat) == alt_colname)
        colnames(dat)[alt_colname_index] <- "Alt"
        dat <- separate_rows(dat, Alt)
        colnames(dat)[alt_colname_index] <- alt_colname
        cat("Done.\n")
    }
    all_starts <- as.numeric(paste0(dat[,pos_colname])) - arm_width; all_ends <- as.numeric(paste0(dat[,pos_colname])) + arm_width 
    all_chromosomes <- paste0(dat[,chr_colname])
    all_refs <- paste0(dat[,ref_colname])
    all_refs_lengths <- nchar(all_refs)
    if(all_alt_combinations) { 
        bases <- c("A","C","G","T")
        all_alts_list <- lapply(bases, function(alt_base) rep(alt_base,nrow(dat)))
        names(all_alts_list) <- bases
    } else { 
        if(skip_alts) { all_alts_list <- list("")
        } else { all_alts_list <- list(gsub(",.*$", "", paste0(dat[,alt_colname]))) }
        names(all_alts_list) <- ""
    }
    sequences_all <- lapply(1:length(all_alts_list), function(alt_base_i) {
        alt_base = names(all_alts_list)[alt_base_i]
        all_alts <- all_alts_list[[alt_base_i]]
        sequences <- mclapply(unique(all_chromosomes), function(chromosome) {
            print(chromosome)
            if(grepl("Y",chromosome)) { return(matrix(nrow=0, ncol=ncol(dat)+2)) }
            cat(paste0("Reading chr", gsub("chr","",chromosome), "..."))
            refseq <- get_refseq(chromosome, version=version, allow_BSgenome=TRUE)[[1]]
            cat("Done.\nGrabbing sequences...")
            curr_chrom_indices <- which(all_chromosomes == chromosome)
            starts <- all_starts[curr_chrom_indices]; ends <- all_ends[curr_chrom_indices]
            widths <- ends - starts + 1; arm_width = floor(widths/2)
            split_indices <- cumsum(widths)
            curr_chrom_rbp_sequences <- eval(parse(text=paste0("substring(paste0(refseq[c(",paste0(starts,":",ends, collapse=","),")]), c(0,split_indices[-length(split_indices)])+1, split_indices)")))
            cat("Done.\n")
            
            if(!skip_alts) {
                mutated_sequences <- paste0(substr(curr_chrom_rbp_sequences, 1, arm_width), all_alts[curr_chrom_indices],substr(curr_chrom_rbp_sequences, arm_width+1+all_refs_lengths[curr_chrom_indices], nchar(curr_chrom_rbp_sequences)))  
                mutated_sequences <- substring(mutated_sequences, 1, widths)
            }
            curr_chrom_rbp_sequences <- substring(curr_chrom_rbp_sequences, 1, widths)
            if(skip_alts) {
                return(cbind(dat[curr_chrom_indices,], c(curr_chrom_rbp_sequences)))
            } else {
                return(cbind(dat[curr_chrom_indices,], c(curr_chrom_rbp_sequences), c(mutated_sequences)))
            }
        }, mc.cores=detectCores())
        sequences <- data.frame(rbindlist(lapply(sequences, function(x) return(data.frame(x)))))
        if(skip_alts) {
            colnames(sequences)[ncol(sequences)] <- "ref_sequence"
        } else {
            alt_sequence_name = paste0("alt",alt_base,"_sequence")
            colnames(sequences)[(ncol(sequences)-1):ncol(sequences)] <- c("ref_sequence", alt_sequence_name)
            sequences[,alt_sequence_name] <- paste0(sequences[,alt_sequence_name])
        }
        sequences$ref_sequence <- paste0(sequences$ref_sequence)
        
        return(sequences)
    })
    if(length(sequences_all) > 1) { 
        sequences_all <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, all.x = TRUE, all.y=TRUE), sequences_all)
    } else { sequences_all <- sequences_all[[1]] } 
    return(sequences_all)
}

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
#full_length = sequence_length*2 - 1
#dat <- add_sequence_context_feature(dat[,-which(colnames(dat) %in% c("ref_sequence", "alt_sequence"))], width=full_length, version="hg19")
dat_refs_tensor <- genomic_sequences_to_tensor(dat$ref_sequence, sequence_length=full_length, frame_width=sequence_length)
dat_alts_tensor <- genomic_sequences_to_tensor(dat$alt_sequence, sequence_length=full_length, frame_width=sequence_length)


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

regulatory_variant_dat <- add_sequence_context_feature(regulatory_variant_dat, width=151, version="hg19")
regulatory_variant_dat[13527,]
control_indices <- grepl("random", regulatory_variant_dat$sample)
regulatory_variant_refs_tensor <- genomic_sequences_to_tensor(regulatory_variant_dat$ref_sequence, sequence_length=151)
regulatory_variant_alts_tensor <- genomic_sequences_to_tensor(regulatory_variant_dat$alt_sequence, sequence_length=151)

asd_variant_dat <- add_sequence_context_feature(asd_variant_dat, width=151, version="hg19")
asd_control_indices <- !asd_variant_dat$Proband #asd_variant_dat$Pheno == "control" #, for An et al dataset
asd_variant_refs_tensor <- genomic_sequences_to_tensor(asd_variant_dat$ref_sequence, sequence_length=151)
asd_variant_alts_tensor <- genomic_sequences_to_tensor(asd_variant_dat$alt_sequence, sequence_length=151)

chd_variant_dat <- add_sequence_context_feature(chd_variant_dat, width=151, version="hg19")
chd_control_indices <- chd_variant_dat$case_control == "control"
chd_variant_refs_tensor <- genomic_sequences_to_tensor(chd_variant_dat$ref_sequence, sequence_length=151)
chd_variant_alts_tensor <- genomic_sequences_to_tensor(chd_variant_dat$alt_sequence, sequence_length=151)


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
group_seqs_into_tensors <- function(dat, re_order=TRUE) {
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
    #if("seq" %in% colnames(dat)) { seqs <- dat$seq
    #} else {
        seqs <- unlist(lapply(1:num_sequences, function(seq_i) rep(seq_i, rows_per_frame)))
    #}
    sequences <- dat$ref_sequence
    AFs <- dat$AF
    af_tensor <- array(0, dim=c(num_sequences, frame_width), dimnames = list(paste0("sequence_",1:num_sequences), 1:frame_width)); #colnames(mat) <- apply(expand.grid(bases, 1:sequence_length)[,c(2,1)], 1, paste0, collapse="")
    
    dat_seq_breakpoints <- c(0,which(diff(seqs)>0),nrow(dat))
    tensor <- abind(lapply(1:num_sequences, function(seq_i) {
        print(seq_i)
        seq_start = dat_seq_breakpoints[seq_i]+1
        seq_end = dat_seq_breakpoints[seq_i+1]
        seq_AFs <- rollapply(AFs[seq_start:seq_end], width=rows_per_position, by=rows_per_position, FUN=sum)
        af_tensor[seq_i,1:length(seq_AFs)] <<- seq_AFs
        return(genomic_sequences_to_tensor(sequences[seq_start], sequence_length=frame_width, verbose=FALSE))
    }), along=4)
    tensor <- aperm(tensor, c(4,2,3,1))
    
    return_env <- new.env()
    return_env[["input"]] <- tensor
    return_env[["output"]] <- af_tensor
    return(return_env)
}

# Returns the expected number of mutations for each width-sized region around the specified positions.
get_expected_muts <- function(dat, chr_colname="Chrom", pos_colname="Position", width=151, version="hg19", aggregate_region=FALSE, background_rate_path="/home/local/ARCS/ak3792/Documents/Research/data/background_mut_rate") {
    if(version!="hg19") { dat <- liftover(dat, from=version, to="hg19", chr_colname=chr_colname, start_colname=pos_colname, confirm_refseq=FALSE, mismatches_pause=FALSE) }
    if(aggregate_region) { expected_muts <- rep(0, nrow(dat))
    } else { expected_muts <- array(0, c(nrow(dat), width)) }
    arm_width = floor(width/2)
    all_starts <- as.numeric(paste0(dat[,pos_colname])) - arm_width; all_ends <- as.numeric(paste0(dat[,pos_colname])) + arm_width 
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
                expected_muts[chr_indices] <- mean(chr_expected_muts)
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
    pLIs_dat <- exac_dat[unlist(sapply(nearest_genes_dat,function(x) { nearest_gene <- which(paste0(exac_dat$gene) == x); if(length(nearest_gene) < 1) { return(NA) } else { return(nearest_gene[1]) }  })),"pLI"]
    pLIs_dat[is.na(pLIs_dat)] <- mean(pLIs_dat[!is.na(pLIs_dat)])
    pLIs_dat <- array(pLIs_dat, c(length(pLIs_dat),1))
    saveRDS(pLIs_dat, output_path("gnomad_pLIs.rds"))
    # Position-specific expected (background mutation rate)
    expected <- get_expected_muts(dat[dat_midpoints,], chr_colname="Chrom_hg19", pos_colname="Position_hg19", width=151, version="hg19")
    #expected <- t(array(expected, dim(dat_gradcams)[2:1]))
    sum(expected== 0)/prod(dim(expected))
    saveRDS(expected, output_path("gnomad_expected.rds"))
}

######################################################
# Code to set up all tensors for machine learning!
######################################################
# Set up initial gnomad data frame with eligible random variants.
width=151
arm_width=floor(width/2)
gnomad <- read.csv(data_path("random_gnomad_variants.txt"), sep="\t", header=FALSE)[,1:5]
gnomad <- cbind(paste0("seq",rep(1:nrow(gnomad))), gnomad, gnomad[,2]-arm_width, gnomad[,2]+arm_width)
colnames(gnomad) <- c("seq", "Chrom", "Position", "sample", "Ref", "Alt", "start", "end")
gnomad <- gnomad[gnomad$Chrom %in% c(1:22),]
gnomad <- add_sequence_context_feature(gnomad, version="hg19", width=width)
colnames(gnomad)[which(colnames(gnomad) %in% c("Position", "ref_sequence"))] <- c("Position_hg19", "ref_sequence_hg19")
gnomad <- liftover(gnomad, from="hg19", to="hg38", chr_colname="Chrom", start_colname="start", end_colname="end", ref_colname="Ref", alt_colname="Alt", confirm_refseq=FALSE, mismatches_pause=FALSE)#[,c("Chrom","Position","Ref","Alt","sample","snv_indel","Chrom_hg38","Position_hg38")]
gnomad <- gnomad[gnomad$Chrom %in% c(1:22) & (gnomad$end - gnomad$start + 1 == width),]
gnomad <- cbind(gnomad, gnomad$start + arm_width); colnames(gnomad)[ncol(gnomad)] <- "Position"
gnomad <- add_sequence_context_feature(gnomad, version="hg38", width=width)
gnomad <- gnomad[gnomad$ref_sequence == gnomad$ref_sequence_hg19, c("seq", "Chrom", "Position", "start", "end", "sample", "Ref", "Alt", "ref_sequence", "Chrom_hg19", "Position_hg19", "start_hg19", "end_hg19")]
nrow(gnomad)
saveRDS(gnomad, output_path("gnomad.rds"))

# Convert gnomad data frame into one with all positions and mutation types fully unpacked and ready to be inputted for ANNOVAR gnomAD annotation.
gnomad_fully_unpacked <- unpack_muts(gnomad)
num_seqs = length(unique(gnomad_fully_unpacked$seq)); num_seqs
seq_row_counts <- table(table(gnomad_fully_unpacked$seq)); seq_row_counts
rows_per_seq = as.numeric(names(seq_row_counts))
rows_per_seq == width*3
gnomad_fully_unpacked_midpoints <- seq((arm_width+1)*3,nrow(gnomad_fully_unpacked),by=width*3)-1
gnomad_hg38_to_hg19_pos_shifts <- gnomad_fully_unpacked$Position_hg19[gnomad_fully_unpacked_midpoints]-gnomad_fully_unpacked$Position[gnomad_fully_unpacked_midpoints]
gnomad_fully_unpacked$Position_hg19 <- gnomad_fully_unpacked$Position + unlist(lapply(gnomad_hg38_to_hg19_pos_shifts, function(x) rep(x,rows_per_seq)))
saveRDS(gnomad_fully_unpacked, output_path("gnomad_fully_unpacked.rds"))

# Annotate fully unpacked gnomad variants data frame with gnomAD 3.0.0, using local ANNOVAR install.
gnomad_fully_unpacked <- readRDS(output_path("gnomad_fully_unpacked.rds"))
gnomad_fully_unpacked_annotated <- run_annovar(gnomad_fully_unpacked, "gnomad30_genome", buildver="hg38")
gnomad_fully_unpacked_annotated$AF[is.na(gnomad_fully_unpacked_annotated$AF) | gnomad_fully_unpacked_annotated$AF == "."] <- 0
gnomad_fully_unpacked_annotated$AF <- as.numeric(gnomad_fully_unpacked_annotated$AF)
table(table(gnomad_fully_unpacked_annotated$seq))
sum(gnomad_fully_unpacked_annotated$AF == 0)/nrow(gnomad_fully_unpacked_annotated)
saveRDS(gnomad_fully_unpacked_annotated, output_path("gnomad_fully_unpacked_annotated.rds"))

# Repackage annotated gnomad data into tensor objects.
gnomad_fully_unpacked_annotated <- readRDS(output_path("gnomad_fully_unpacked_annotated.rds"))
gnomad_tensors_env <- group_seqs_into_tensors(gnomad_fully_unpacked_annotated, re_order=FALSE)
dat_tensor <- gnomad_tensors_env[["input"]]
af_tensor <- gnomad_tensors_env[["output"]]
rm(gnomad_tensors_env)
dat_tensor[1,,,]
sum(af_tensor== 0)/prod(dim(af_tensor))
sum(rowSums(af_tensor)== 0)/nrow(af_tensor)
saveRDS(dat_tensor, output_path("gnomad_dat.rds"))
saveRDS(af_tensor, output_path("gnomad_AFs.rds"))

# Calculate GradCAMs from dat tensor.
dat_tensor <- readRDS(output_path("gnomad_dat.rds"))
calculate_gradcams(dat_tensor[,,,], output_path("gnomad_gradcams"))
dat_gradcams <- get_gradcams(output_path("gnomad_gradcams"))
dat_gradcams <- array_reshape(dat_gradcams, c(dim(dat_gradcams),1))
dat_gradcams[is.na(dat_gradcams)] <- 0
sum(dat_gradcams== 0)/prod(dim(dat_gradcams))
#sum(rowSums(af_tensor)== 0)/nrow(dat_gradcams)
saveRDS(dat_gradcams, output_path("gnomad_gradcams.rds"))

# Calculate expected (background mutation rate) and pLI additional input tensors.
dat <- readRDS(output_path("gnomad_fully_unpacked_annotated.rds"))
dat_midpoints <- seq((arm_width+1)*3,nrow(dat),by=width*3)-1
dat_variant_granges_hg19 <- to_genomic_regions(dat[dat_midpoints,], chr_colname="Chrom_hg19", start_colname="start_hg19", end_colname="end_hg19", label_colname="seq")
dat_variant_granges_hg19 <- dat_variant_granges_hg19[order(as.numeric(gsub("seq","",names(dat_variant_granges_hg19))))]
genebody_hg19 <- read.csv(data_path("refGene\\refGene_hg19_genebody_fixed.bed"), sep="\t", stringsAsFactors=FALSE)
genebody_granges_hg19 <- to_genomic_regions(genebody_hg19, chr_colname="chromosome", start_colname="tss", end_colname="tes", strand_colname="strand", label_colname="gene", order_coordinates=TRUE, remove_duplicate_labels=TRUE)
# Nearest gene pLI (haploinsufficiency); gnomAD obs/exp for nearest gene would be more parallel to AF constraint
# Load ExAC pLIs data
exac_dat <- read.csv(data_path("fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt"), sep="\t")
nearest_genes_dat <- names(genebody_granges_hg19)[nearest(dat_variant_granges_hg19, genebody_granges_hg19)]
pLIs_dat <- exac_dat[unlist(sapply(nearest_genes_dat,function(x) { nearest_gene <- which(paste0(exac_dat$gene) == x); if(length(nearest_gene) < 1) { return(NA) } else { return(nearest_gene[1]) }  })),"pLI"]
pLIs_dat[is.na(pLIs_dat)] <- mean(pLIs_dat[!is.na(pLIs_dat)])
pLIs_dat <- array(pLIs_dat, c(length(pLIs_dat),1))
saveRDS(pLIs_dat, output_path("gnomad_pLIs.rds"))
# Position-specific expected (background mutation rate)
expected <- get_expected_muts(dat[dat_midpoints,], chr_colname="Chrom_hg19", pos_colname="Position_hg19", width=151, version="hg19")
expected <- t(array(expected, dim(dat_gradcams)[2:1]))
sum(expected== 0)/prod(dim(expected))
saveRDS(expected, output_path("gnomad_expected.rds"))

# Load prediction to likelihood mapping table.
prediction_score_to_likelihood_mapping_table <- read.csv(output_path("prediction_score_to_likelihood_mapping_table.csv"))

#################################################################################################################
# Super-models code and analysis
#################################################################################################################
save_supermodel <- function(model, model_name) {
    filename = output_path(paste0(model_name))
    export_savedmodel(model, filename)
}
save_supermodel(model, "suprnova")
load_supermodel <- function(model_name) {
    filename = output_path(paste0(model_name))
    model <- load_model_tf(filename)
    return(model)
}

supermodel <- function(dat, name=NULL, version="hg19") {
    # gnomAD sample size, used for converting allele count to allele frequency
    sample_size = 71702
    # Load all input and output tensors, which were calculated and saved in the previous section.
    af_tensor <- readRDS(output_path("gnomad_AFs.rds"))
    dat_gradcams <- readRDS(output_path("gnomad_gradcams.rds"))
    pLIs_dat <- readRDS(output_path("gnomad_pLIs.rds"))
    expected <- readRDS(output_path("gnomad_expected.rds"))
    dat <- readRDS(output_path("gnomad_fully_unpacked_annotated.rds"))
    width=151; arm_width=floor(width/2); dat_midpoints <- seq((arm_width+1)*3,nrow(dat),by=width*3)-1
    #af_tensor <- (af_tensor - mean(af_tensor)) / sd(af_tensor)
    #dat_gradcams <- (dat_gradcams - mean(dat_gradcams)) / sd(dat_gradcams)
    #pLIs_dat <- (pLIs_dat - mean(pLIs_dat)) / sd(pLIs_dat)
    #expected <- (expected - mean(expected)) / sd(expected)

    # Draw some important input and output metrics/distributions to file.
    filename = output_path("gnomAD_distributions.pdf")
    pdf(file=filename, width=11)
    par(mfrow=c(2,3))
    plot(density(log10(expected)), type="l", lwd=1, col="blue", main="Background mutation rates", xlab="log10(mu)", cex.lab=1.2)
    mtext(paste0(formatC(100*sum(expected==0)/prod(dim(expected)),format="f",digits=2),"% of sites have 0 mu"))
    plot(density(log10(af_tensor)), type="l", lwd=1, col="blue", main="gnomAD allele frequencies", xlab=paste0("log10(AF), ",sample_size," WGS samples"), cex.lab=1.2)
    mtext(paste0(formatC(100*sum(af_tensor==0)/prod(dim(af_tensor)),format="f",digits=2),"% of sites have 0 AF"))
    plot(density(log10(af_tensor[rowSums(af_tensor)>0,])), type="l", lwd=1, col="blue", main="gnomAD allele frequencies", xlab=paste0("log10(AF), ",sample_size," WGS samples"), cex.lab=1.2)
    mtext(paste0(formatC(100*sum(af_tensor[rowSums(af_tensor)>0,]==0)/prod(dim(af_tensor[rowSums(af_tensor)>0,])),format="f",digits=2),"% of >0 region sites have 0 AF"))
    plot(density(pLIs_dat), type="l", lwd=1, col="blue", main="Nearest gene pLIs", xlab="pLI (haploinsufficiency)", cex.lab=1.2)
    plot(density(sample(dat_gradcams,100000,replace=TRUE)), type="l", lwd=1, col="blue", main="RBP GradCAM values", xlab="normalized activation", cex.lab=1.2)
    mtext(paste0(formatC(100*sum(dat_gradcams==0)/prod(dim(dat_gradcams)),format="f",digits=2),"% of sites have 0 activation"))
    plot(density(sample(dat_gradcams,100000,replace=TRUE)), type="l", lwd=1, col="blue", main="RBP GradCAM values", xlab="normalized activation", cex.lab=1.2)
    mtext(paste0(formatC(100*sum(dat_gradcams==0)/prod(dim(dat_gradcams)),format="f",digits=2),"% of sites have 0 activation"))
    par(mfrow=c(1,1))
    dev.off()
    pdf_to_png(filename)
    
    # Get dat trimers
    dat_trimers <- rbindlist(lapply(1:length(dat_midpoints), function(i) {
        print(i)
        ref_seq <- dat$ref_sequence[dat_midpoints[i]]
        return(data.frame(t(data.frame(rollapply(strsplit(ref_seq,"")[[1]], width=3, by=1, partial=2, FUN=function(x) paste0(x,collapse=""))))))
    }))
    # Find CpG sites
    dat_cpg <- t(apply(dat_trimers, 1, function(x) grepl("CG", x)))
    # Sanity check: Make sure that background mutation rate is correlated with CpG sites
    cor(c(expected), c(dat_cpg), method="spearman")
    
    maps <- round(af_tensor * sample_size * 2 + 0.01)
    maps <- t(apply(maps, 1, function(x) return(c(sum(x==0),sum(x==1),sum(x==2),sum(x==3),sum(x==4),sum(x==5),sum(x>5)))))[,-1]
    colnames(maps) <- c("Singletons", "Doubletons", "Tripletons", "AC = 4", "AC = 5", "AC > 5")
    barplot(colSums(maps)/sum(maps)*100, las=2, main="", ylab="Proportion (%)") #, ylim = c(0,5 + max(mtcars$qsec)), xlab = "", space = 1)
    end_point = 0.5 + ncol(maps) + ncol(maps) - 1 #this is the line which does the trick (together with barplot "space = 1" parameter)
    #rotate 60 degrees (srt = 60)
    text(seq(1.5, end_point, by = 2), par("usr")[3]-5, srt=60, adj=1, xpd=TRUE, labels=paste(colnames(maps)), cex=0.65)
    
    filename = output_path("proportion_of_singletons.pdf")
    pdf(file=filename)
    plot(density(maps), main="PS Distribution in Data", xlab="Proportion of Singletons", col="blue", lwd=2, cex.lab=1.4, cex.axis=1.4, cex.main=1.3)
    mtext("Proportion of singletons measured for each 151bp sequence", cex=1.2)
    dev.off()
    pdf_to_png(filename)
    allowed_indices = which(rowSums(expected)>0)
    cor(log10(rowMeans(expected)[allowed_indices]), maps[allowed_indices], method="spearman")
    #plot(log10(rowMeans(expected)[allowed_indices]), maps[allowed_indices], main="", xlab="", ylab="Proportion of Singletons")
    draw_plot(data.frame(x=log10(rowMeans(expected)[allowed_indices]), y=maps[allowed_indices]), hex_density = 25, title="PS vs. Mutability", xlab="log10(mean regional mutability)", ylab="Proportion of Singletons", legend_text="# regions", linear_best_fit=FALSE, quadratic_best_fit=FALSE, filename="PS_vs_mutability.pdf")
    
    ac_tensor <- round(af_tensor * sample_size * 2 + 0.01)
    contexts <- paste0(unique(unlist(dat_trimers))); contexts <- sort(contexts[sapply(contexts, nchar)==3])
    #contexts <- unique(sapply(contexts, function(x) paste0(sort(c(x,sequence_composite(x))),collapse="/")))
    context_sfs <- t(sapply(contexts, function(context) {
        context_indices <- matrix(unlist(dat_trimers == context | dat_trimers == sequence_composite(context)), nrow=nrow(dat_trimers), ncol=ncol(dat_trimers))
        #context_indices <- context_indices & ac_tensor[context_indices] < mean()
        background <- expected[context_indices]
        x <- ac_tensor[context_indices]
        return(c(mean(background), sum(x==0),sum(x==1),sum(x==2),sum(x==3),sum(x==4),sum(x==5),sum(x>5)))
    }))
    colnames(context_sfs) <- c("mutability", "Zerotons", "Singletons", "Doubletons", "Tripletons", "AC = 4", "AC = 5", "AC > 5")
    context_ps <- sort(context_sfs[,3]/rowSums(context_sfs[,-c(1:2)]))
    rownames(context_sfs) <- sapply(contexts, function(x) paste0(sort(c(x,sequence_composite(x))),collapse="/"))
    context_sfs <- context_sfs[which(!duplicated(rownames(context_sfs))),]
    context_weights <- rowSums(context_sfs[,-c(1:2)])
    sort(context_sfs[,3]/context_weights)
    maps_y <- context_sfs[,3]/context_weights
    maps_y <- unlist(sapply(1:length(maps_y), function(i) rep(maps_y[i],context_weights[i])))
    x <- log10(context_sfs[,1])
    x <- unlist(sapply(1:length(x), function(i) rep(x[i],context_weights[i])))
    maps_model <- lm(maps_y ~ x )
    filename = output_path("proportion_of_singletons_vs_mutability.pdf")
    pdf(file=filename)
    cols <- sapply(rownames(context_sfs), function(x) return(c("blue","green")[as.numeric(grepl("CG",x))+1]))
    plot(log10(context_sfs[,1]), context_sfs[,3]/rowSums(context_sfs[,-c(1:2)]), col=cols, main="Proportion of Singletons vs. Mutability", xlab="log10(mutability)", ylab="Proportion of Singletons", cex.lab=1.4, cex.lab=1.4, cex.main=1.3)
    mtext("Trimers pooled with their respective sequence composite", cex=1.1)
    maps_line_x <- seq(floor(min(log10(context_sfs[,1]))), ceiling(max(log10(context_sfs[,1]))), by=diff(range(log10(context_sfs[,1])))/100)
    lines(maps_line_x, predict(maps_model, new=data.frame(x=maps_line_x)), col="red", lwd=2)
    legend("bottomleft", legend=c("CpG trimers","Non-CpG trimers","Expected PS regression"), col=c("green","blue","red"), pch=c(15,15,NA), lty=c(NA,NA,1), lwd=c(NA,NA,2), cex=1.2)
    dev.off()
    pdf_to_png(filename)
        
    #ps <- t(apply(ps, 1, function(x) return(c(sum(x==0),sum(x==1),sum(x==2),sum(x==3),sum(x==4),sum(x==5),sum(x>5)))))[,-1]
    #colnames(ps) <- c("Singletons", "Doubletons", "Tripletons", "AC = 4", "AC = 5", "AC > 5")
    #ps <- ps[,1]/rowSums(ps); ps[is.nan(ps)] <- 0
    #mutability <- rowMeans(expected)
    #train_indices <- randomized_order_training_indices
    #maps_model <- lm(ps[intersect(train_indices,which(mutability>0))] ~ log10(mutability[intersect(train_indices,which(mutability>0))]))
    #draw_plot(data.frame(x=log10(mutability[intersect(train_indices,which(mutability>0))]), y=ps[intersect(train_indices,which(mutability>0))]), hex_density = 25, title="PS vs. Mutability", xlab="log10(mean regional mutability)", ylab="Proportion of Singletons", legend_text="# regions", linear_best_fit=T, quadratic_best_fit=FALSE)#, filename="PS_vs_mutability.pdf")
    
    # Calculate correlations between GradCAM input and allele frequency output.
    rbps <- get_features_by_group("RBP")
    gradcam_af_correlations <- unlist(lapply(1:length(rbps), function(rbp_i) {
        print(rbp_i)
        return(cor(c(dat_gradcams[,,rbp_i,1]), c(af_tensor), method="spearman"))
    }))
    range(gradcam_af_correlations)
    rbps[order(gradcam_af_correlations, decreasing=TRUE)][1:10]
    gradcam_af_correlations[order(gradcam_af_correlations, decreasing=TRUE)][1:10]
    
    fetal_heart_H3K36me3 <- load_annotation("E083.H3K36me3.broadPeak")
    fetal_heart_H3K36me3 <- intersect(fetal_heart_H3K36me3, fetal_heart_H3K36me3)
    sum(width(fetal_heart_H3K36me3))
    
    # Number of filters to use
    num_filters = 100
    ## Number of RBPs in matrix
    #num_rbps = dim(scores_tensor)[[3]] #dim(DF)[[2]]

    # Start with hidden 2D convolutional layer being fed (num_rbps + 2) x 2  pixel images
    #rbp_binding_input <- layer_input(name="rbp_binding_input", shape = c(dim(DF)[[2]], dim(DF)[[3]], 1))
    #DF <- aperm(DF, c(1,3,2))
    rbp_binding_input <- layer_input(name="rbp_binding_input", shape = c(dim(dat_gradcams)[-1])) #c(dim(DF)[2:3]))
    # First convolution layer
    conv1 <- layer_conv_2d(name="conv1", filters=num_filters, kernel_size=c(10,160), padding="valid", use_bias=FALSE)# %>% layer_activation_leaky_relu(alpha=0.3)
    a <- rbp_binding_input %>% conv1
    # Second convolution layer
    #conv2 <- layer_conv_1d(name="conv2", filters=num_filters, kernel_size=8, use_bias=FALSE) %>% layer_activation_leaky_relu(alpha=0.3)
    # Third convolution layer
    #conv3 <- layer_conv_1d(name="conv3", filters=num_filters, kernel_size=8, use_bias=FALSE) %>% layer_activation_leaky_relu(alpha=0.3)
    # First max pooling layer
    maxpool1 <- layer_max_pooling_2d(name="maxpool1", pool_size=c(4,num_filters))
    # Second max pooling layer
    #maxpool2 <- layer_max_pooling_1d(name="maxpool2", pool_size = 4)
    
    # Module responsible for processing tensor with CNN and learning relevant RBP binding patterns.
    rbp_disruption_module <- rbp_binding_input %>% 
        conv1 %>% layer_reshape(unlist(conv1$output_shape[c(2,4,3)])) %>%  layer_activation_leaky_relu(alpha=0.3) %>%
        maxpool1 %>% layer_dropout(0.1, name="dropout_weak1") %>% layer_flatten()
    # Upsamples compressed rbp disruption encoding back up to the input sequence length, and outputs it as a flattened vector.
    rbp_disruption_module_output <- rbp_disruption_module %>% 
        layer_reshape(c(rbp_disruption_module$shape[[2]],1)) %>%
        layer_upsampling_1d(ncol(dat_gradcams)) %>% layer_average_pooling_1d(rbp_disruption_module$shape[[2]]) %>%
        layer_flatten(name="rbp_disruption_output")

    # Auxiliary input for gene expression
    gene_expression_module <- layer_input(name="gene_expression", shape = c(1)) 
    gene_expression_module_activation <- gene_expression_module %>% layer_dense(1, name="gene_expression_activation", activation="sigmoid")
    # Auxiliary input for pLI/haploinsufficiency
    pli_module <- layer_input(name="pLI", shape = c(1))
    pli_module_activation <- pli_module %>% layer_dense(1, name="pLI_activation", activation="sigmoid", use_bias=FALSE)
    # Learning of selection coefficient using RBP gene regulation disruption output and gene-level features.
    gene_damagingness_layer <- layer_multiply(c(rbp_disruption_module_output, pli_module_activation)) %>%
        layer_activation_softmax(name="gene_damagingness")
    
    ## Secondary structure prediction output
    #ss_predictions <- cnn_module %>% 
    #    layer_dense(name="ss_predictions", sequence_length, activation="sigmoid")

    # Auxiliary input for background mutation rate
    background_mut_rate_module <- layer_input(name="background_mut_rate", shape = c(151))
    
    #pseudocount = min(expected[expected!=0])
    #expected <- expected + pseudocount
    #af_tensor <- af_tensor + pseudocount
    # Compile model, and draw model network.
    #model <- keras_model(inputs=c(rbp_binding_input, pli_module, background_mut_rate_module), outputs=c(AF_emission_layer))
    
    mu_scaling_factor = 1e4 # In effect, the starting point for s
    actual_s_upper_bound = 0.1 # 0.04
    s_lower_bound = 1/(actual_s_upper_bound*mu_scaling_factor) # Because s_real = 1/(mu_scaling_factor * s_calc) -> s_calc = 1/(s_real * mu_scaling_factor), and we bound 0.001 < s_real < 0.04
    s_upper_bound = 1/(0.00005*mu_scaling_factor) # 1/0.00005*mu_scaling_factor)
    print(paste0("Constraints: ",1/(mu_scaling_factor*s_upper_bound)," <= s <= ",1/(mu_scaling_factor*s_lower_bound)))
    s_constraint <- function(w) { 
        constraint_maxnorm(s_upper_bound)(w) * k_cast(k_greater_equal(w, s_lower_bound), k_floatx()) 
    }
    adjusted_selection_coef_layer <- layer_lambda(name="s", f = function(inputs) { 
        s <- inputs[[1]]
        return(s_constraint(s))
    }) (c(gene_damagingness_layer)) #(c(gene_damagingness_layer))
    
    positions <- lapply(rep(0,ncol(dat_gradcams)), function(x) x)
    for(i in 1:ncol(dat_gradcams)) {
        positions[[i]] <- layer_lambda(name=paste0("af",i), f = function(inputs) { 
            mu <- inputs[[1]]; s <- inputs[[2]]
            return(mu[,i,drop=FALSE] * s[,i,drop=FALSE])
        }) (c(background_mut_rate_module, adjusted_selection_coef_layer))
        #positions[[i]] <- layer_lambda(name=paste0("s",i), f = function(inputs) { 
        #    w <- inputs[[1]]
        #    return(constraint_maxnorm(s_upper_bound)(w) * k_cast(k_greater_equal(w, s_lower_bound), k_floatx()))
        #}) (positions[[i]])
        #positions[[i]] <- positions[[i]] %>% layer_dense(1, name=paste0("pos",i), activation="relu", kernel_constraint=s_constraint, use_bias=FALSE) #kernel_constraint=constraint_nonneg
    }
    # Allele frequency emission per position, given positional background rates and selection coefficients.
    AF_emission_layer <- layer_concatenate(positions, name="AF_emission")
    model <- keras_model(inputs=c(rbp_binding_input, pli_module, background_mut_rate_module), outputs=c(AF_emission_layer))
    
    # Custom loss function, that can optionally take in a mask arg and additional extra args.
    batch_size = 128
    custom_loss_function <- function(extra_args=NULL, mask=-1) {
        sample_size = array(extra_args[1], c(batch_size,ncol(dat_gradcams)))
        loss_function <- function(y_true, y_pred) { # mask=mask
            #pseudocount = k_variable(extra_args[2])
            #mask_vector <- k_cast(k_not_equal(y_true, mask), k_floatx());
            #y_pred <- y_pred * mask_vector; y_true <- y_true * mask_vector
            
            #y_pred <- expected[10130:10257,,drop=FALSE]; y_true <- af_tensor[10130:10257,,drop=FALSE]
            #y_pred <- expected[1,1:5,drop=FALSE]; y_true <- af_tensor[1,1:5,drop=FALSE]
            
            #return(k_random_normal(1) * k_sum(y_true %>% (tfp$distributions$Poisson(rate=c(0.1))$prob)))
            # obs <- 0.2 %>% (tfp$distributions$Poisson(rate=c(0.1,1e-8,1e-8/0.02))$prob)
            # 0 %>% (tfp$distributions$Poisson(rate=c(1e-9,1e-8))$prob)
            #obs %>% (tfp$distributions$Poisson(rate=c(0.1,1e-8,1e-8/0.02))$prob)
            #y_true <- y_true * sample_size # observed
            
            #k_print_tensor(y_true)
            #k_print_tensor(k_sum(y_pred))
            #k_print_tensor(y_pred$shape)
            #k_print_tensor(y_pred[1:2,1:5,drop=FALSE])
            #k_print_tensor((y_pred + pseudocount)[1:2,1:5,drop=FALSE])
            #k_print_tensor(k_sum(k_cast(k_equal(y_pred, nan)), k_floatx()))
            
            # YOU LEFT OFF HERE!!!!
            #sum(apply(expected, 1, function(x) sum(x==0)>0))
            #sum(apply(expected, 1, function(x) sum(x==0)>1))
            #nrow(expected)
            
            eligible_indices <- tf$where(k_flatten(k_less(y_true, 1e-2) & k_greater(y_pred, 0))) #which(as.matrix(k_flatten(k_less(y_true, 1e-3))))
            y_true <- k_gather(k_flatten(y_true), eligible_indices)
            y_pred <- k_gather(k_flatten(y_pred), eligible_indices)
            log_probs <- k_mean(y_true) %>% (tfp$distributions$Poisson(rate=(k_mean(y_pred)))$log_prob)
            #tf$gather(log_probs, mask_vector)
            #mask_vector <- k_cast(k_not_equal(y_pred, NaN) & k_not_equal(y_pred, -Inf), k_floatx())
            #k_print_tensor(log_probs)
            sum_log_probs <- -k_mean(log_probs) #* mask_vector)
            #k_print_tensor(sum_log_probs)
            
            #mask_vector <- k_cast(k_not_equal(y_true, mask), k_floatx());
            #y_pred <- y_pred * mask_vector; y_true <- y_true * mask_vector
            #k_not_equal(obs, -Inf)
            return(sum_log_probs)
            
            # KL divergence of two Poissons is lambda1 * log(lambda1/lambda2) + lambda2 - lambda1
            #return(k_log(y_pred * k_log(y_pred / y_true) + y_true - y_pred))
            
            #return(-(y_pred %>% tfd_log_prob(y_true))) 
            #y_true <- y_true * (1/sample_size)
            ## return log(KL(Poisson_pred, Poisson_true))
            
        }
        return(loss_function)
    }
    opt <- optimizer_sgd(lr = 50, decay = 1e-2, momentum=0.5, nesterov=TRUE, clipnorm=1) #optimizer_adam(clipnorm=1)
    model %>% compile(loss=custom_loss_function(sample_size, pseudocount), optimizer="rmsprop", metrics="mae") # custom_loss_function(sample_size)
    model_name = "supermodel"
    #plot_model(model, to_file = output_path(paste0("images/",model_name,".png")), show_shapes = TRUE, show_layer_names = TRUE)
    
    norm_tensor <- function(x) { return((x - mean(x)) / sd(x)) }
    #expected <- array(expected, c(dim(expected),1))
    # Fit model to training data
    #randomized_order_training_indices = sample(1:nrow(dat_gradcams), floor(0.8*nrow(dat_gradcams)))
    allowed_indices = which(rowSums(expected)>0 & rowSums(af_tensor)>0)
    randomized_order_training_indices = sample(allowed_indices, floor(0.8*length(allowed_indices)))
    test_indices = allowed_indices[!(allowed_indices %in% randomized_order_training_indices)]
    ##history <- model %>% fit(list(array(DF, c(dim(DF), 1))[randomized_order_training_indices,,,,drop=FALSE], pLIs_dat[randomized_order_training_indices]), list(labels[randomized_order_training_indices]), epochs=10, batch_size=32, validation_split = 0.2)
    #history <- model %>% fit(list(dat_gradcams[randomized_order_training_indices,,,,drop=FALSE], pLIs_dat[randomized_order_training_indices,,drop=FALSE], expected[randomized_order_training_indices,,drop=FALSE]), list(af_tensor[randomized_order_training_indices,,drop=FALSE]*(1/sample_size)), epochs=5, batch_size=128, validation_split=0.2)
    ##history <- model %>% fit(list(norm_tensor(dat_gradcams[randomized_order_training_indices,,,,drop=FALSE]), norm_tensor(pLIs_dat[randomized_order_training_indices,,drop=FALSE]), norm_tensor(expected[randomized_order_training_indices,,drop=FALSE])), list(norm_tensor(af_tensor[randomized_order_training_indices,,drop=FALSE]*(1/sample_size))), epochs=5, batch_size=128, validation_split=0.2)
    num_epochs=150
    history <- model %>% fit(list(dat_gradcams[randomized_order_training_indices,,,,drop=FALSE], pLIs_dat[randomized_order_training_indices,,drop=FALSE], expected[randomized_order_training_indices,,drop=FALSE]*mu_scaling_factor), list(af_tensor[randomized_order_training_indices,,drop=FALSE]), epochs=num_epochs, batch_size=128, validation_split=0.2)
    print(history$metrics$loss)
    print(history$metrics$val_loss)
    
    print(paste0("Constraints: ",1/(mu_scaling_factor*s_upper_bound)," <= s <= ",1/(mu_scaling_factor*s_lower_bound)))
    plot(1:length(history$metrics$val_loss), history$metrics$val_loss, type="l", col="blue", lwd=2, xlab="epochs", ylab="loss", main="Training loss vs. epochs")
    lines(1:length(history$metrics$val_loss), history$metrics$loss, type="l", col="red", lwd=2)
    
    filename=output_path("training_loss_vs_epochs.pdf")
    pdf(file=filename)
    plot(1:length(history$metrics$val_loss), history$metrics$val_loss, type="l", col="blue", lwd=2, xlab="epochs", ylab="loss", main="Training loss vs. epochs", cex.lab=1.4, cex.axis=1.4, cex.main=1.3)
    #lines(1:length(history$metrics$val_loss), history$metrics$loss, type="l", col="red", lwd=2)
    dev.off()
    pdf_to_png(filename)
    
    # Write model architecture to file.
    write(paste0(model), file=output_path("model_architecture.txt"))

    pLI_act <- as.numeric(model$get_layer("pLI_activation")$weights[[1]][1,1])
    pLI_act
    
    #s_vals <- 1/(mu_scaling_factor*unlist(lapply(1:ncol(dat_gradcams), function(i) as.numeric(model$get_layer(paste0("pos",i))$weights[[1]][1,1]))))
    #s_vals[s_vals == -Inf | s_vals == Inf] <- 0
    #s_vals
    
    # Directly look at weights in the trained model to determine learned motifs.
    conv1_weights <- get_weights(model$get_layer("conv1"))[[1]]
    conv1_weights
    conv1_kernel_size = as.numeric(gsub("[^0-9]", "", paste0(model$get_layer("conv1")$kernel_size)))
    conv1_kernel_size
    sapply_out <- sapply(1:num_filters, function(i) { 
        filter <- t(conv1_weights[,,,i])
        #filter <- -filter
        rownames(filter) <- c("A","C","G","T")
        print("TGCATG")
        print(sequence_composite("TGCATG"))
        filter
        filter <- apply(filter, 2, function(x) { if(sum(x < 0)>0) { x <- x - min(x) }; return(x / sum(x)) })
        print(paste0("Calculating PWM for filter c1f",i,"..."))
        #paste0(model_name,"_c1f",i))
    })
    pli_weight <- get_weights(model$get_layer("pLI_activation"))[[1]][1,1]
    
    # Calculate expected (background mutation rate) and pLI additional input tensors.
    dat <- readRDS(output_path("gnomad_fully_unpacked_annotated.rds"))
    dat_midpoints <- seq((arm_width+1)*3,nrow(dat),by=width*3)-1
    dat_variant_granges_hg19 <- to_genomic_regions(dat[dat_midpoints,], chr_colname="Chrom_hg19", start_colname="start_hg19", end_colname="end_hg19", label_colname="seq")
    dat_variant_granges_hg19 <- dat_variant_granges_hg19[order(as.numeric(gsub("seq","",names(dat_variant_granges_hg19))))]

    eclip_padding = 0
    rbp_eclip_peaks <- sapply(get_features_by_group("RBP"), function(rbp) { print(rbp); x <- load_annotation(paste0(rbp,"_",eclip_padding,"bp")); return(intersect(x,x)) })
    rbp_names <- names(rbp_eclip_peaks)
    all_rbp_eclip_peaks <- Reduce(c, rbp_eclip_peaks)
    names(all_rbp_eclip_peaks) <- unlist(lapply(1:length(rbp_eclip_peaks), function(rbp_i) {
        rep(rbp_names[rbp_i], length(rbp_eclip_peaks[[rbp_i]]))
    }))
    dat_eclip_overlaps <- data.frame(findOverlaps(dat_variant_granges_hg19, all_rbp_eclip_peaks))
    dat_eclip_starts <- start(dat_variant_granges_hg19[dat_eclip_overlaps$queryHits])
    dat_eclip_overlaps <- cbind(dat_eclip_overlaps, names(all_rbp_eclip_peaks)[dat_eclip_overlaps$subjectHits], start(all_rbp_eclip_peaks[dat_eclip_overlaps$subjectHits])-dat_eclip_starts, end(all_rbp_eclip_peaks[dat_eclip_overlaps$subjectHits])-dat_eclip_starts)
    colnames(dat_eclip_overlaps) <- c("query", "subject", "RBP", "start", "end")
    #dat_eclip_intersections <- merge(data.frame(queryHits=c(1:length(dat_variant_granges_hg19))), all.x=TRUE)

    # Investigate trained model layers by feeding in new sequences.
    layers_to_investigate = c("conv1", "maxpool1", "rbp_disruption_output", "gene_damagingness", "s", "AF_emission")
    model_layer <- model$get_layer(layers_to_investigate[1])
    model_layer_names <- sapply(model$layers, function(layer) paste0(layer$name))
    layers_to_investigate <- layers_to_investigate[layers_to_investigate %in% model_layer_names]
    activation_model <- keras_model(model$input, outputs=lapply(layers_to_investigate, function(layer_to_investigate) { return(get_output_at(model$get_layer(layer_to_investigate), 1)) })) #sapply(model$layers[2:3], function(layer) return(get_output_at(layer[[1]], 1)))) #layer[[1]]$output)
    output_folder = output_path("model_outputs")
    dir.create(output_folder, showWarnings = FALSE)

    test_sequence_index = test_indices
    labels_all <- af_tensor[test_sequence_index,,drop=FALSE]
    background_all <- expected[test_sequence_index,,drop=FALSE]
    activations <- activation_model %>% predict(list(dat_gradcams[test_sequence_index,,,,drop=FALSE], pLIs_dat[test_sequence_index,,drop=FALSE], expected[test_sequence_index,,drop=FALSE]*mu_scaling_factor))
    rbp_disruption_all <- activations[[which(layers_to_investigate == "rbp_disruption_output")]]
    rbp_disruption_range <- range(rbp_disruption_all)
    rbp_disruption_all <- (rbp_disruption_all - rbp_disruption_range[1])/diff(rbp_disruption_range)
    preds_all <- activations[[which(layers_to_investigate == "AF_emission")]]
    s_vals_all <- 1/(mu_scaling_factor*activations[[which(layers_to_investigate == "s")]])
    s_vals_all[s_vals_all == -Inf | s_vals_all == Inf] <- 0
    af_range <- c(0,7.5e-4)
    which(test_indices %in% dat_eclip_overlaps$query)
    test_indices[which(test_indices %in% dat_eclip_overlaps$query)]
    maps_result <- rbindlist(lapply(1:length(test_indices), function(i) {
        print(i)
        test_sequence_index = test_indices[i]
        #test_sequence_index = 28209
        single_input = length(test_sequence_index) == 1
        labels <- af_tensor[test_sequence_index,,drop=FALSE]
        background <- expected[test_sequence_index,,drop=FALSE]
        pli <- pLIs_dat[i]
        rbp_disruption <- rbp_disruption_all[i,]
        #rbp_disruption
        Ne = 10000
        num_simulated = 100000
        preds <- preds_all[i,]
        zero_preds <- preds == 0
        preds[zero_preds] <- rpois(sum(preds==0), lambda=num_simulated*Ne*background[zero_preds])/num_simulated
        preds_simulated <- preds; preds[zero_preds] <- NA #rep(NA, sum(zero_preds))
        #preds
        s_vals <- s_vals_all[i,]
        s_vals_zero <- s_vals; s_vals_zero[s_vals_zero == 0] <- runif(sum(s_vals_zero == 0))*0.01
        s_vals_nonzero <- s_vals; s_vals_nonzero[s_vals==0] <- NA
        #s_vals
        #cor(c(s_vals), c(preds), method="spearman")
        #cor(c(preds > 0), c(labels > 0), method="spearman")
        #cor(c(s_vals > 0), c(labels > 0), method="spearman")
        
        # Draw some important input and output metrics/distributions to file.
        filename = full_path(output_folder, paste0("model_outputs_seq",test_sequence_index,".pdf"))
        if(test_sequence_index %in% test_indices[which(test_indices %in% dat_eclip_overlaps$query)]) {
            print("Has eCLIP overlap!")
            filename <- gsub(".pdf", "_eclip.pdf", filename)
        }
        pdf(file=filename, height=12)
        par_mar <- par()$mar
        num_plots = 5
        cex_axis = 1.4
        cex_lab = 1.5
        total_mar_to_remove = 4.1
        bonus_mar_to_remove = c(0.75,-0.5,0,0)
        custom_par_mars <- lapply(seq(1,0,by=-(1/(num_plots-1))), function(mar_top_bottom_proportion) { 
            return(par_mar - total_mar_to_remove*c(mar_top_bottom_proportion, 0, 1-mar_top_bottom_proportion, 0) - bonus_mar_to_remove)
        })
        par(mfrow=c(num_plots,1), mar=custom_par_mars[[1]])
        plot(1:length(labels), labels, type="l", col="blue", xaxt="n", xlab="", main="Supermodel Activation Analysis", ylim=af_range, cex.axis=cex_lab, cex.lab=cex_lab, cex.main=1.4)
        #mtext(paste0("Sequence ",test_sequence_index,", pLI = ",formatC(pli,format="f",digits=2)), cex=1.1)
        mtext(paste0("Sequence ",test_sequence_index), cex=1.1)
        par(mar=custom_par_mars[[2]])
        plot(1:length(labels), background, type="l", col="green", xaxt="n", xlab="", ylim=c(range(background_all)), cex.axis=cex_lab, cex.lab=cex_lab)
        par(mar=custom_par_mars[[3]])
        plot(1:length(labels), rbp_disruption, type="l", col="black", xaxt="n", xlab="", yaxs="i", ylim=c(-0.08,1), cex.axis=cex_lab, cex.lab=cex_lab)
        abline(h=0, col="black", lty=1)
        eclip_annotation_loc = -0.04
        abline(h=eclip_annotation_loc, col="gray20", lty=3)
        eclip_overlaps <- dat_eclip_overlaps[dat_eclip_overlaps$query == test_sequence_index,]
        if(nrow(eclip_overlaps)>0) {
            overlapped_rbps <- unique(eclip_overlaps$RBP)
            cols <- rainbow(length(overlapped_rbps)); names(cols) <- overlapped_rbps
            #col_alpha = 0.4; cols <- adjustcolor(cols, alpha.f=col_alpha)
            segments(x0=eclip_overlaps$start, y0=eclip_annotation_loc, x1=eclip_overlaps$end, y1=eclip_annotation_loc, col=cols[paste0(eclip_overlaps$RBP)], lty=1, lwd=3)
            legend("topleft", legend=overlapped_rbps, col=cols, pch=15)
        }
        par(mar=custom_par_mars[[4]])
        plot(1:length(labels), s_vals_zero, type="l", col="grey", xaxt="n", xlab="", ylim=c(range(s_vals_all)), cex.axis=cex_lab, cex.lab=cex_lab) #ylim=c(0,1/(mu_scaling_factor*s_lower_bound))
        lines(1:length(labels), s_vals_nonzero, type="l", col="red")
        par(mar=custom_par_mars[[5]])
        plot(1:length(labels), preds_simulated, type="l", col="grey", xlab="Position", ylim=c(ylim=af_range), cex.axis=cex_lab, cex.lab=cex_lab)
        lines(1:length(labels), preds, type="l", col="red")
        par(mfrow=c(1,1), mar=par_mar)
        dev.off()
        pdf_to_png(filename, output_folder="model_outputs")
        
        maps <- round(labels * sample_size * 2 + 0.01)
        maps_res <- data.frame(t(rbind(maps, background, s_vals, rep(nrow(eclip_overlaps)>0,151))))
        colnames(maps_res) <- c("AC", "mutability", "s", "nearby_eclip")
        return(maps_res)
    }))
    
    #################################################################################################################
    # Apply SUPRNOVA model to WGS datasets
    #################################################################################################################
    # Set up ASD data frame.
    asd_variant_dat <- read.csv(data_path("ASD_OlgaT.tsv"), sep="\t") # a <- read.table(file="/home/local/ARCS/ak3792/Documents/Research/WGS/output/regulatory_variants.tsv", sep="\t", header=1) # For annotation with WGSA
    #asd_variant_dat <- read.csv(data_path("WGS/An_2018_Table_S2_denovos.csv"), skip=1, header=TRUE, sep="\t")
    #asd_variant_dat[,c(1:2,6,3:4)]
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
        dat_pLIs_filename = output_path(paste0(dat_name,"_pLIs.rds"))
        dat_expected_filename = output_path(paste0(dat_name,"_expected.rds"))
        
        tasks_done = 0
        # Processed data / genomic regions
        if(tasks_done < num_tasks_to_do && (!file.exists(dat_regions_filename) || rewrite)) {
            print(paste0("Setting up ",dat_regions_filename,"..."))
            dat <- standardize_colnames(dat, remove_chr_prefix=TRUE, re_order=TRUE)
            dat <- cbind(paste0("seq",rep(1:nrow(dat))), dat[,c(1:2,5,3:4)], dat[,2]-arm_width, dat[,2]+arm_width)
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
            dat <- readRDS(dat_regions_filename)
            dat_trimers <- rbindlist(lapply(1:nrow(dat), function(i) {
                if(i %% 1000 == 1) { print(i) }
                ref_seq <- dat$ref_sequence[i]
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
            exac_dat <- read.csv(data_path("fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt"), sep="\t")
            nearest_genes_dat <- names(genebody_granges_hg19)[nearest(dat_granges_hg19, genebody_granges_hg19)]
            pLIs_dat <- exac_dat[unlist(sapply(nearest_genes_dat,function(x) { nearest_gene <- which(paste0(exac_dat$gene) == x); if(length(nearest_gene) < 1) { return(NA) } else { return(nearest_gene[1]) }  })),"pLI"]
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
    
    # Setup CHD/SSC data frames (data freeze in accepted Nature paper, with FreeBayes)
    chd_variant_dat <- read.csv(data_path("WGS/CHDFB/chdfb_sscfb_pcgc_final_variants.csv"))
    chd_variant_dat$sample <- paste0(chd_variant_dat$case_control,"_",chd_variant_dat$sample)
    chd_indices <- grepl("case_",chd_variant_dat$sample); ssc_indices <- !chd_indices
    subsample_indices <- c(sample(which(chd_indices),5000), sample(which(ssc_indices),5000))
    setup_supermodel_data(chd_variant_dat[subsample_indices,], "chd", version="hg19")
    chd_dat <- readRDS(output_path("chd_dat.rds"))
    nrow(chd_dat)
    table(grepl("case_",chd_dat$sample))
    setup_supermodel_data(NULL, "chd", version="hg19", num_tasks_to_do=1)

    # Setup ASD data frames
    setup_supermodel_data("asd", num_tasks_to_do=1)
    chd_dat <- readRDS(output_path("asd_all_dat.rds"))

    # Interpet s prediction results for WGS datasets
    analyze_s_predictions <- function(dat_name, case_pat=NULL, case_indices=NULL) {
        dat <- readRDS(output_path(paste0(dat_name,"_dat.rds")))
        dat_granges_hg19 <- readRDS(output_path(paste0(dat_name,"_granges_hg19.rds")))
        if(is.null(case_indices)) { case_indices <- grepl(case_pat, dat$sample) }; control_indices <- !case_indices
        num_case_variants = sum(case_indices); num_control_variants = sum(control_indices)
        dat_s_preds_filename = output_path(paste0(dat_name,"_s_preds.rds"))
        if(file.exists(dat_s_preds_filename)) {
            s_vals_all <- readRDS(dat_s_preds_filename)
        } else {
            dat_gradcams <- readRDS(output_path(paste0(dat_name,"_gradcams.rds")))
            dat_expected <- readRDS(output_path(paste0(dat_name,"_expected.rds")))
            dat_pLIs <- readRDS(output_path(paste0(dat_name,"_pLIs.rds")))
            dat_cpg <- readRDS(output_path(paste0(dat_name,"_cpg.rds")))
            if(!exists("activation_model") || is.null(activation_model)) { 
                layers_to_investigate = c("s")
                model_layer <- model$get_layer(layers_to_investigate[1])
                model_layer_names <- sapply(model$layers, function(layer) paste0(layer$name))
                layers_to_investigate <- layers_to_investigate[layers_to_investigate %in% model_layer_names]
                activation_model <- keras_model(model$input, outputs=lapply(layers_to_investigate, function(layer_to_investigate) { return(get_output_at(model$get_layer(layer_to_investigate), 1)) })) #sapply(model$layers[2:3], function(layer) return(get_output_at(layer[[1]], 1)))) #layer[[1]]$output)
                set_global(activation_model)
            }
            activations <- activation_model %>% predict(list(dat_gradcams[,,,,drop=FALSE], dat_pLIs[,,drop=FALSE], dat_expected[,,drop=FALSE]*mu_scaling_factor))
            if(class(activations) == "list") {
                s_vals_all <- 1/(mu_scaling_factor*activations[[which(layers_to_investigate == "s")]])
            } else { # activations is a matrix instead of a list, caused by having just a single entity in layers_to_investigate
                s_vals_all <- 1/(mu_scaling_factor*activations)
            }
            s_vals_all[s_vals_all == -Inf | s_vals_all == Inf] <- 0
            saveRDS(s_vals_all, dat_s_preds_filename)
        }
        cutoff=0.02
        sum(s_vals_all[case_indices,76] > cutoff)/(num_case_variants*151); sum(s_vals_all[control_indices,76] > cutoff)/(sum(num_control_variants)*151)
        
        filename = output_path(paste0(dat_name,"_s_distributions.pdf"))
        pdf(filename)
        plot(density(s_vals_all[case_indices,76]), col="red", xlab="Selection Coefficient (s)", main=paste0(toupper(dat_name)," Selection Coef. Distributions"), cex.lab=1.4, cex.axis=1.4, cex.main=1.3)
        lines(density(s_vals_all[control_indices,76]), col="blue")
        mtext(paste0("All variants"), cex=1.2)
        variant_counts <- c(prod(dim(s_vals_all[case_indices,76,drop=FALSE])),prod(dim(s_vals_all[control_indices,76,drop=FALSE])))
        legend("topright", legend=c(paste0(c("case (N = ","control (N = "),variant_counts,")"), paste0("ratio = ",variant_counts[1]/variant_counts[2])), col=c("red","blue","white"), pch=15, cex=1.2)
        dev.off()
        pdf_to_png(filename)

        filename = output_path(paste0(dat_name,"_s_distributions_righttail.pdf"))
        pdf(filename)
        plot(density(s_vals_all[case_indices,76]), xlim=c(cutoff,0.11), xaxs="i", ylim=c(0,2), col="red", xlab="Selection Coefficient (s)", main=paste0(toupper(dat_name)," Selection Coef. Distributions"), cex.lab=1.4, cex.axis=1.4, cex.main=1.3)
        lines(density(s_vals_all[control_indices,76]), col="blue")
        mtext(paste0("Variants with s >= 0.02 only"), cex=1.2)
        variant_counts <- c(sum(s_vals_all[case_indices,76]>=0.02),sum(s_vals_all[control_indices,76]>=0.02))
        legend("topright", legend=c(paste0(c("case (N = ","control (N = "),variant_counts,")"), paste0("ratio = ",variant_counts[1]/variant_counts[2])), col=c("red","blue","white"), pch=15, cex=1.2)
        dev.off()
        pdf_to_png(filename)

        dat_eclip_overlaps <- find_eclip_overlaps(dat_granges_hg19)
        case_variants_eclip <- intersect(unique(dat_eclip_overlaps$query), which(case_indices))
        control_variants_eclip <- intersect(unique(dat_eclip_overlaps$query), which(control_indices))
        print(paste0("Freq. regions with nearby eCLIP in cases: ",length(case_variants_eclip) / num_case_variants))
        print(paste0("Freq. regions with nearby eCLIP in controls: ",length(control_variants_eclip) / num_control_variants))
        
        filename = output_path(paste0(dat_name,"_s_distributions_eclip.pdf"))
        pdf(filename)
        plot(density(s_vals_all[case_variants_eclip,76]), col="red", xlab="Selection Coefficient (s)", main=paste0(toupper(dat_name)," Selection Coef. Distributions"), cex.lab=1.4, cex.axis=1.4, cex.main=1.3)
        lines(density(s_vals_all[control_variants_eclip,76]), col="blue")
        mtext(paste0("Variants with nearby eCLIP only"), cex=1.2)
        variant_counts <- c(length(case_variants_eclip),length(control_variants_eclip))
        legend("topright", legend=c(paste0(c("case (N = ","control (N = "),variant_counts,")"), paste0("ratio = ",variant_counts[1]/variant_counts[2])), col=c("red","blue","white"), pch=15, cex=1.2)
        dev.off()
        pdf_to_png(filename)

        filename = output_path(paste0(dat_name,"_s_distributions_eclip_righttail.pdf"))
        pdf(filename)
        plot(density(s_vals_all[case_variants_eclip,76]), xlim=c(0.02,0.11), xaxs="i", ylim=c(0,2), col="red", xlab="Selection Coefficient (s)", main=paste0(toupper(dat_name)," Selection Coef. Distributions"), cex.lab=1.4, cex.axis=1.4, cex.main=1.3)
        lines(density(s_vals_all[control_variants_eclip,76]), col="blue")
        mtext(paste0("Variants with nearby eCLIP and s >= 0.02 only"), cex=1.2)
        variant_counts <- c(sum(s_vals_all[case_variants_eclip,76]>=0.02),sum(s_vals_all[control_variants_eclip,76]>=0.02))
        legend("topright", legend=c(paste0(c("case (N = ","control (N = "),variant_counts,")"), paste0("ratio = ",variant_counts[1]/variant_counts[2])), col=c("red","blue","white"), pch=15, cex=1.2)
        dev.off()
        pdf_to_png(filename)
    }
    analyze_s_predictions("chd", case_pat="^case_")
    analyze_s_predictions("asd", case_pat="-p")

    # Read ASD data frames
    asd_dat <- readRDS(output_path("ASD_OlgaT.rds"))[1:40000,]
    asd_dat_is_case <- grepl("-p", asd_dat$sample) # FALSE means controls / unaffected sibling rather than proband
    asd_num_case_variants = sum(asd_dat_is_case); asd_num_control_variants = length(asd_dat_is_case) - asd_num_case_variants
    asd_gradcams <- readRDS(output_path("asd_gradcams.rds"))
    asd_expected <- readRDS(output_path("asd_expected.rds"))
    asd_pLIs <- readRDS(output_path("asd_pLIs.rds"))
    # Find CpG sites
    asd_cpg <- readRDS(output_path("asd_cpg.rds"))
    # Sanity check: Make sure that background mutation rate is correlated with CpG sites
    cor(c(asd_expected), c(asd_cpg), method="spearman")
    mean(asd_expected[asd_cpg]); mean(asd_expected[!asd_cpg])
    # Interpet ASD s prediction results
    layers_to_investigate = c("s")
    model_layer <- model$get_layer(layers_to_investigate[1])
    model_layer_names <- sapply(model$layers, function(layer) paste0(layer$name))
    layers_to_investigate <- layers_to_investigate[layers_to_investigate %in% model_layer_names]
    activation_model <- keras_model(model$input, outputs=lapply(layers_to_investigate, function(layer_to_investigate) { return(get_output_at(model$get_layer(layer_to_investigate), 1)) })) #sapply(model$layers[2:3], function(layer) return(get_output_at(layer[[1]], 1)))) #layer[[1]]$output)
    
    asd_activations <- activation_model %>% predict(list(asd_gradcams[,,,,drop=FALSE], asd_pLIs[,,drop=FALSE], asd_expected[,,drop=FALSE]*mu_scaling_factor))
    s_vals_all <- 1/(mu_scaling_factor*asd_activations[[which(layers_to_investigate == "s")]])
    s_vals_all[s_vals_all == -Inf | s_vals_all == Inf] <- 0
    sort(s_vals_all[asd_dat_is_case,], decreasing=TRUE)[1:10]; sort(s_vals_all[!asd_dat_is_case,], decreasing=TRUE)[1:10]
    cutoff=0.09
    sum(s_vals_all[asd_dat_is_case,] > cutoff)/(asd_num_case_variants*151); sum(s_vals_all[!asd_dat_is_case,] > cutoff)/(sum(asd_num_control_variants)*151)
    plot(density(s_vals_all[asd_dat_is_case,76]), xlim=c(0.02,0.11), ylim=c(0,1), col="red", main="")
    lines(density(s_vals_all[!asd_dat_is_case,76]), col="blue")
    mtext(table(asd_dat_is_case))
    
    # ASD eCLIP overlaps
    asd_granges_hg19 <- to_genomic_regions(asd_dat, chr_colname="Chrom_hg19", start_colname="start_hg19", end_colname="end_hg19", label_colname="seq")
    asd_granges_hg19 <- asd_granges_hg19[order(as.numeric(gsub("seq","",names(asd_granges_hg19))))]
    dat_eclip_overlaps <- find_eclip_overlaps(asd_granges_hg19)
    paste0("Freq. regions with nearby eCLIP in cases: ",length(intersect(unique(dat_eclip_overlaps$query), which(asd_dat_is_case))) / asd_num_case_variants)
    paste0("Freq. regions with nearby eCLIP in controls: ",length(intersect(unique(dat_eclip_overlaps$query), which(!asd_dat_is_case))) / asd_num_control_variants)
    
    split_points <- seq(-0.5, 0, by=0.1)
    pred_score_enrichments <- lapply(split_points, function(split_point) {
        enrichment_test <- multi_threshold_test(pred_scores[labels[test_indices] < split_point,,drop=FALSE], pred_scores[labels[test_indices] >= split_point,,drop=FALSE], seq(-0.05, 0.05, by=0.005), threshold_dir="<", label=paste0("split=",split_point))
        filename = output_path(paste0("pred_score_distributions_",split_point,"split.pdf"))
        pdf(filename)
        plot(density(pred_scores[labels[test_indices] < split_point]), type="l", lwd=2, col="red", cex.axis=1.4, cex.lab=1.4, cex.main=1.3, main="AF_lognorm Prediction Distributions", xlab="normalized ln(AF+pseudocount)")
        lines(density(pred_scores[labels[test_indices] >= split_point], from=min_unconstrained_pred_score), lwd=2, col="blue")
        optimal_threshold = enrichment_test$threshold[which.min(enrichment_test[,which(grepl(":p.value", colnames(enrichment_test)))])]
        abline(v=optimal_threshold, lty=2)
        mtext(paste0("Optimal constraint threshold: ",optimal_threshold), cex=1.2)
        legend("topleft", legend=c(paste0("labels < ",split_point), paste0("labels >= ",split_point)), col=c("red", "blue"), pch=15, cex=1.3)
        dev.off()
        pdf_to_png(filename)
        return(enrichment_test)
    }) %>% reduce(full_join, by="threshold")
    pred_score_enrichments[1:5,1:9]
    apply(pred_score_enrichments[,c(1,which(grepl(":p.value", colnames(pred_score_enrichments))))], 2, min)
    pred_score_enrichments[,c(1,which(grepl(":estimate", colnames(pred_score_enrichments))))]
    # Print optimal cutoff and p.value for each split_point value.
    t(apply(pred_score_enrichments[,which(grepl(":p.value", colnames(pred_score_enrichments)))], 2, function(x) { optimal = which.min(x); return(c(pred_score_enrichments[optimal,c(1)], min(x))) }))
    write.csv(pred_score_enrichments, file=output_path(paste0("pred_scores_enrichments.csv")), row.names=FALSE)
    #write.csv(pred_score_enrichments[,c(1,which(grepl(":estimate", colnames(greater_than_threshold_fets_ss))),which(grepl(":m1", colnames(greater_than_threshold_fets_ss))))], file=output_path(paste0(rbp,"_gene_expression_enrichment_estimates_vs_secondary_structure.csv")), row.names=FALSE)
    
    filename = output_path(paste0("pred_score_enrichments.pdf"))
    ci_width = 0.001
    mtext_label = paste0(nrow(pred_scores)," test variants")
    cols <- rainbow(length(split_points))
    sapply_out <- sapply(1:length(split_points), function(ss_i) { #1:length(split_points)
        ss <- paste0("split=",split_points[ss_i])
        col <- cols[ss_i]
        if(ss_i == 1) {
            pdf(filename)
            plot(pred_score_enrichments$threshold, pred_score_enrichments[,paste0(ss,":estimate")], main=paste0("RBP Binding Enrichment (high vs. low expressed)"), xlab="Gap between high/low expression (ln(TPM+1), around median)", ylab="log2(Odds Ratio of RBP binding)", col=col, type="l", lwd=2, pch=19, xaxs="i", yaxs="i", ylim=c(min(c(0,min(unlist(pred_score_enrichments[,grepl(":conf.int_lower",colnames(pred_score_enrichments))])))), max(unlist(pred_score_enrichments[,grepl(":conf.int_higher",colnames(pred_score_enrichments))]))), cex.axis=1.4, cex.lab=1.4, cex.main=1.4)
            #for(quant in quantile(DF$gene_expression)[2:4]) { abline(v=quant, col="black", lty=3) }
            abline(h=0, col="black", lty=1)
            mtext(mtext_label, cex=1.1)
        } else {
            lines(pred_score_enrichments$threshold, pred_score_enrichments[,paste0(ss,":estimate")], col=cols[ss_i], type="l", lwd=2, pch=19)
        }
        ci_col <- adjustcolor(col, alpha.f=0.3)
        segments(x0=pred_score_enrichments$threshold, y0=pred_score_enrichments[,paste0(ss,":conf.int_lower")], x1=pred_score_enrichments$threshold, y1=pred_score_enrichments[,paste0(ss,":conf.int_higher")], col=ci_col)
        segments(x0=pred_score_enrichments$threshold-ci_width, y0=pred_score_enrichments[,paste0(ss,":conf.int_lower")], x1=pred_score_enrichments$threshold+ci_width, y1=pred_score_enrichments[,paste0(ss,":conf.int_lower")], col=ci_col)
        segments(x0=pred_score_enrichments$threshold-ci_width, y0=pred_score_enrichments[,paste0(ss,":conf.int_higher")], x1=pred_score_enrichments$threshold+ci_width, y1=pred_score_enrichments[,paste0(ss,":conf.int_higher")], col=ci_col)
    })
    legend("bottomleft", legend=c(split_points,"quantiles"), col=c(cols,"black"), pch=c(rep(19,length(split_points)),NA), lty=c(rep(NA,length(split_points)),3), cex=1.05)
    dev.off()
    pdf_to_png(filename) 
    
    
    
    cor(maps_result$AC, maps_result$mutability, method="spearman")
    cor(maps_result$AC, maps_result$s, method="spearman")
    draw_plot(data.frame(x=maps_result$AC[maps_result$AC > 0] == 1, y=maps_result$s[maps_result$AC > 0], hex_density = 25, title="PS vs. Mutability", xlab="AC", ylab="s", legend_text="# regions", linear_best_fit=FALSE, quadratic_best_fit=FALSE)) #, filename="PS_vs_mutability.pdf")
    range(maps_result$s[maps_result$AC == 1]); range(maps_result$s[maps_result$AC > 1]) 
    plot(density(maps_result$s[maps_result$AC == 1 & maps_result$s > 0]), col="red")
    lines(density(maps_result$s[maps_result$AC > 1 & maps_result$s > 0]), col="blue")
    
    s_by = 0.01
    #s_bucket_cuts <- c(0,0.01,0.04,1)
    s_bucket_cuts <- rbind(s_bucket_cuts[-length(s_bucket_cuts)], s_bucket_cuts[-1])
    s_bucket_cuts <- rbind(c(0,0,0,rep(s_by*3,1/s_by-3)), c(s_by,seq(s_by*2,1,by=s_by)))
    s_bucket_cuts <- rbind(c(0,0,0,0.03),c(0.01,0.02,0.03,1))
    s_bucket_names <- c("s = 0", paste0(s_bucket_cuts[1,]," < s <= ",s_bucket_cuts[2,]))
    sites = "all" # "all", "CpG", or "non-CpG"
    maps_cpg_vector <- unlist(dat_cpg[test_indices,])
    maps_sfs <- rbindlist(lapply(0:ncol(s_bucket_cuts), function(s_i) {
        if(s_i == 0) { curr_indices <- maps_result$s == 0 } else { curr_indices <- maps_result$s > s_bucket_cuts[1,s_i] & maps_result$s <= s_bucket_cuts[2,s_i]  }
        if(sites == "CpG") { curr_indices <- curr_indices & maps_cpg_vector 
        } else if (sites == "non-CpG") { curr_indices <- curr_indices & !maps_cpg_vector
        } else { sites = "all" }
        x <- maps_result$AC[curr_indices]
        ps <- sum(x == 1)/sum(x >= 1)
        #expected_ps <- 
        expected_ps <- predict(maps_model, new=data.frame(x=log10(maps_result$mutability[curr_indices])))
        maps_sampled <- ps - sort(sapply(1:1000, function(sample_i) mean(sample(expected_ps, length(expected_ps), replace=TRUE))), decreasing=TRUE)
        maps <- ps - mean(expected_ps)
        return(data.frame(t(data.frame(c(sum(x==0),sum(x==1),sum(x==2),sum(x==3),sum(x==4),sum(x==5),sum(x>5),ps,mean(expected_ps),maps,maps_sampled[25],maps_sampled[975])[-1]))))
    }))
    #maps_sfs <- t(apply(maps_sfs, 1, function(x) x/sum(x)))
    colnames(maps_sfs) <- c("Singletons", "Doubletons", "Tripletons", "AC = 4", "AC = 5", "AC > 5", "PS_observed", "PS_expected", "MAPS", "MAPS_lower", "MAPS_upper")
    rownames(maps_sfs) <- s_bucket_names
    s_last_filled = min(c(max(which(apply(maps_sfs[,1:6], 1, function(x) sum(!is.nan(x) & x > 0)>0))), max(which(diff(maps_sfs$MAPS)>0))+1))
    maps_sfs_rownames <- c(rownames(maps_sfs)[1:(s_last_filled-1)], gsub("<=.*", "<= 1", rownames(maps_sfs)[s_last_filled]))
    maps_sfs <- maps_sfs[1:s_last_filled,]; rownames(maps_sfs) <- maps_sfs_rownames
    maps_sfs
    maps_sfs_to_write <- cbind(rownames(maps_sfs),maps_sfs); colnames(maps_sfs_to_write)[1] <- "AC"
    write.table(maps_sfs_to_write, sep=",", output_path(paste0("model_sfs_",sites,"_table.csv")), row.names=FALSE)
    if(sites == "all") {
        filename = output_path(paste0("model_MAPS.pdf"))
        pdf(file=filename)
        maps_cols <- c("black", rev(rainbow(nrow(maps_sfs)-1)))
        plot(1:nrow(maps_sfs), maps_sfs$MAPS, col=maps_cols, main="MAPS vs. Predicted Selection Coef.", xlab="", ylab="MAPS", ylim=range(c(maps_sfs[,c("MAPS_lower","MAPS_upper")])), xaxt="n", cex.axis=1.3, cex.lab=1.3, cex.main=1.3)
        abline(h=0, col="gray50", lty=3)
        segments(x0=1:nrow(maps_sfs), y0=maps_sfs$MAPS_lower, x1=1:nrow(maps_sfs), y1=maps_sfs$MAPS_upper, col=maps_cols)
        text(1:nrow(maps_sfs)+0.2, par("usr")[3]-0.001, srt=45, adj=1, xpd=TRUE, labels=rownames(maps_sfs), cex=1)
        dev.off()
        pdf_to_png(filename)
    }
    
    maps_sfs_full <- maps_sfs
    maps_sfs <- maps_sfs_full[,1:6]
    rownames(maps_sfs) <- maps_sfs_rownames
    
    filename = output_path(paste0("model_sfs_",sites,".pdf"))
    pdf(file=filename)
    s_cols <- rainbow(nrow(maps_sfs))
    plot(1:ncol(maps_sfs), maps_sfs[1,]/sum(maps_sfs[1,])*100, col=s_cols[1], type="o", main="Site Frequency Spectrum", xlab="", ylab="Proportion (%)", xaxt="n", xlim=c(1,ncol(maps_sfs)), ylim=c(0,ceiling(max(maps_sfs_full$PS_observed*10))*10), cex.axis=1.3, cex.lab=1.3, cex.main=1.3)
    for(maps_i in 1:ncol(maps_sfs)) { abline(v=maps_i, lty=3, col="gray50") }
    for(s_i in 2:nrow(maps_sfs)) { lines(1:ncol(maps_sfs), maps_sfs[s_i,]/sum(maps_sfs[s_i,])*100, col=s_cols[s_i], type="o") }
    if(sites == "CpG") { mtext_text = paste0(sum(maps_cpg_vector)," supermodel test sequence CpG sites")
    } else if (sites == "non-CpG") { mtext_text = paste0(sum(!maps_cpg_vector)," supermodel test sequence non-CpG sites")
    } else { mtext_text = paste0(length(test_indices)*151," supermodel test sequence sites") }
    mtext(mtext_text, cex=1.1)
    legend("topright", legend=c("Predicted selection coef. s", paste0(rownames(maps_sfs)," (N = ",rowSums(maps_sfs[,1:6]),")")), col=c("white",s_cols), pch=15, cex=1)
    text(1:ncol(maps_sfs)+0.3, par("usr")[3]-1, srt=45, adj=1, xpd=TRUE, labels=colnames(maps_sfs), cex=1.3)
    dev.off()
    pdf_to_png(filename)
    
    filename = output_path(paste0("model_sfs_",sites,"_logscale.pdf"))
    pdf(file=filename)
    s_cols <- rainbow(nrow(maps_sfs))
    plot(log10(1:ncol(maps_sfs)), maps_sfs[1,]/sum(maps_sfs[1,])*100, col=s_cols[1], type="o", main="Site Frequency Spectrum", xlab="", ylab="Proportion (%)", xaxt="n", xlim=log10(c(1,ncol(maps_sfs))), ylim=c(0,ceiling(max(maps_sfs_full$PS_observed*10))*10), cex.axis=1.3, cex.lab=1.3, cex.main=1.3)
    for(maps_i in 1:ncol(maps_sfs)) { abline(v=log10(maps_i), lty=3, col="gray50") }
    for(s_i in 2:nrow(maps_sfs)) { lines(log10(1:ncol(maps_sfs)), maps_sfs[s_i,]/sum(maps_sfs[s_i,])*100, col=s_cols[s_i], type="o") }
    if(sites == "CpG") { mtext_text = paste0(sum(maps_cpg_vector)," supermodel test sequence CpG sites")
    } else if (sites == "non-CpG") { mtext_text = paste0(sum(!maps_cpg_vector)," supermodel test sequence non-CpG sites")
    } else { mtext_text = paste0(length(test_indices)*151," supermodel test sequence sites") }
    mtext(mtext_text, cex=1.1)
    legend("topright", legend=c("Predicted selection coef. s", paste0(rownames(maps_sfs)," (N = ",rowSums(maps_sfs[,1:6]),")")), col=c("white",s_cols), pch=15, cex=1)
    text(log10(1:ncol(maps_sfs))+0.04, par("usr")[3]-1, srt=45, adj=1, xpd=TRUE, labels=colnames(maps_sfs), cex=1.3)
    dev.off()
    pdf_to_png(filename)
    
    length(maps_result$s)
    length(c(expected[test_indices,]))
    ps <- sum(maps_result$AC == 1)/sum(maps_result$AC > 1)
    expected_ps <- predict(maps_model, new=data.frame(x=log10(maps_result$mutability)))
    cor(maps[c(expected[test_indices,])>0], maps_result$s[c(expected[test_indices,])>0], method="spearman")
    
    #filename = output_path("proportion_of_singletons.pdf")
    #pdf(file=filename)
    plot(density(maps), main="PS Distribution in Data", xlab="Proportion of Singletons", col="blue", lwd=2, cex.lab=1.4, cex.axis=1.4, cex.main=1.3)
    mtext("Proportion of singletons measured for each 151bp sequence", cex=1.2)
    #dev.off()
    #pdf_to_png(filename)
    allowed_indices = which(rowSums(expected)>0)
    cor(log10(rowMeans(expected)[allowed_indices]), maps[allowed_indices], method="spearman")
    #plot(log10(rowMeans(expected)[allowed_indices]), maps[allowed_indices], main="", xlab="", ylab="Proportion of Singletons")
    draw_plot(data.frame(x=log10(rowMeans(expected)[allowed_indices]), y=maps[allowed_indices]), hex_density = 25, title="PS vs. Mutability", xlab="log10(mean regional mutability)", ylab="Proportion of Singletons", legend_text="# regions", linear_best_fit=FALSE, quadratic_best_fit=FALSE, filename="PS_vs_mutability.pdf")
    
    
    
    maps <- t(apply(maps, 1, function(x) return(c(sum(x==0),sum(x==1),sum(x==2),sum(x==3),sum(x==4),sum(x==5),sum(x>5)))))[,-1]
    colnames(maps) <- c("Singletons", "Doubletons", "Tripletons", "AC = 4", "AC = 5", "AC > 5")
    
    try_to_align=TRUE
    num_best_sequence_candidates = 2
    motif_length = 8
    motif_length = min(c(motif_length, conv1_kernel_size))
    filter_contributions <- matrix(data=0, nrow=length(test_sequence_index), ncol=num_filters); colnames(filter_contributions) <- paste0("c1f",1:num_filters)
    activated_sequences <- lapply(1:length(test_sequence_index), function(i) {  #input_tensor_dims[length(input_tensor_dims)-1]
        print(paste0(i))
        #sequence <- input_sequence[i]
        conv1_windows <- unfactorize(data.frame(rollapply(1:input_tensor_dims[length(input_tensor_dims)-1], width=conv1_kernel_size, function(positions) return(c(tensor_to_genomic_sequences(input_tensor[i,positions,]), positions[1], positions[conv1_kernel_size])))))
        colnames(conv1_windows) <- c("motif", "start", "end")
        conv1_activation_sequences <- conv1_windows$motif #rollapply(1:input_tensor_dims[length(input_tensor_dims)-1], width=conv1_kernel_size, function(positions) tensor_to_genomic_sequences(input_tensor[i,positions,]))
        sequences_acceptable <- nchar(conv1_activation_sequences) >= motif_length
        conv1_activation_sequences <- conv1_activation_sequences[sequences_acceptable]
        conv1_windows <- conv1_windows[sequences_acceptable,]
        
        repped_grange <- rbp_granges[rows_to_pick][test_sequence_index[i]] # rbp_granges[which(sequence == gr_sequences)]
        
        conv1_activation_scores <- lapply(1:num_filters, function(filter_index) { activation_scores <- activations[[1]][i,sequences_acceptable,filter_index]; names(activation_scores) <- 1:length(conv1_activation_sequences); return(sort(activation_scores, decreasing=TRUE)) })
        names(conv1_activation_scores) <- paste0("c1f",1:length(conv1_activation_scores))
        best_sequences <- sapply(1:num_filters, function(filter_index) { 
            x <- conv1_activation_scores[[filter_index]]
            filter_contributions[i,filter_index] <<- max(x[1:num_best_sequence_candidates])
            best_sequence_indices <- as.numeric(names(x)[1:num_best_sequence_candidates])
            repped_granges <- rep(repped_grange, num_best_sequence_candidates)
            start(repped_granges) <- start(repped_granges) + conv1_windows$start[best_sequence_indices] - 1; end(repped_granges) <- start(repped_granges) + conv1_windows$end[best_sequence_indices] - conv1_windows$start[best_sequence_indices]
            return(repped_granges) 
        })
        names(best_sequences) <- paste0("c1f",1:num_filters)
        return(best_sequences)
    })
    activated_sequences <- lapply(1:num_filters, function(filter_index) { do.call(c, lapply(activated_sequences, function(x) x[[filter_index]])) })
    #activated_sequences <- data.frame(rbindlist(activated_sequences))
    filter_contributions <- colMeans(filter_contributions)
    top_k = 3
    top_k_filters <- order(filter_contributions, decreasing=TRUE)[1:top_k]
    sapply(top_k_filters, function(i) {
        print(paste0("Calculating PWM for filter c1f",i,"..."))
        pwm_analysis(paste0(model_name,"_",rbp_to_analyze,"_c1f",i), activated_sequences[[i]], bucket_size=5000, bucket=1)
    })
    print(paste0("Calculating combined PWM for ",rbp_to_analyze,"..."))
    pwm_analysis(paste0(model_name,"_",rbp_to_analyze,"_combined"), do.call(c, sapply(1:top_k, function(k) { activated_seqs <- activated_sequences[[top_k_filters[k]]]; return(sample(activated_seqs, floor(length(activated_seqs)*filter_contributions[top_k_filters[k]]))) })), bucket_size=5000, bucket=1) #do.call(c, activated_sequences[top_k_filters])
    print(sort(filter_contributions, decreasing=TRUE))
    print("TGCATG")
    print("CATGCA")
    
    pdf(output_path(paste0(model_name,"_",rbp_to_analyze,"_conv1_filter_contribution.pdf")))
    barplot_cols <- c(rainbow(top_k), rep("grey", num_filters-top_k))
    barplot(sort(filter_contributions, decreasing=TRUE), col=barplot_cols, main="CNN Conv1 Filter Contribution", ylab="Average sequence max activation", xlab=paste0(num_filters," Conv1 filters"), cex.main=1.4, cex.lab=1.4, cex.axis=1.4, names.arg=rep("",num_filters))
    legend("topright", title="Filter", legend=c(paste0("c1f",top_k_filters), "other"), col=barplot_cols[1:(top_k+1)], pch=15, cex=1.4)
    dev.off()
    pdf_to_png(output_path(paste0(model_name,"_",rbp_to_analyze,"_conv1_filter_contribution.pdf")))
} # end of supermodel code
    
    
    model <- keras_model(inputs=c(background_mut_rate_module), outputs=c(AF_emission_layer)) 
    model %>% compile(loss="mse", optimizer=optimizer_rmsprop(clipnorm=1), metrics="mse") # custom_loss_function(sample_size)
    history <- model %>% fit(list(norm_tensor(expected[randomized_order_training_indices,,drop=FALSE])), list(norm_tensor(af_tensor[randomized_order_training_indices,,drop=FALSE])), epochs=100, batch_size=128, validation_split=0.2)
    #history <- model %>% fit(list(expected[randomized_order_training_indices,,drop=FALSE]), list(af_tensor[randomized_order_training_indices,,drop=FALSE]), epochs=100, batch_size=128, validation_split=0.2)
    print(history$metrics)
    
    #    layer_independent_poisson(event_shape=c(151), convert_to_tensor_fn=tfp$distributions$Poisson(rate=0.1)$log_prob)
    # dpois(2, lambda=0.1); dpois(2, lambda=0.2)
    # 2 %>% (tfp$distributions$Poisson(rate=c(0.1,0.2))$prob)
    # background mu ~ 1e-8, s ~ 0.02, af ~ 1e-4, sample_size N ~ 15000
    # quantile(af_tensor[af_tensor!=0]); mean(af_tensor)
    # quantile(expected[expected!=0]); mean(expected)
    # obs = af * N ~ 0-2, which is the allele count in gnomAD
    # obs %>% (tfp$distributions$Poisson(rate=c(0.1,1e-8,1e-8/0.02))$prob)
    # af / N = mu / s  ->  af = (mu * N)/s  AND  mu = (af * s)/
    # mean(af_tensor)/sample_size is the expectation 
    layer_lambda(name="AF_emission", f = function(inputs) {
        mu <- inputs[[1]]; s <- inputs[[2]]
        obs %>% (tfp$distributions$Poisson(rate=mu*sample_size)$prob)
        return(k_log(s * (mu)))
        #return(k_log(s * mu + pseudocount))
        #return(k_log(k_sum(rpois(sample_size, lambda=mu))/sample_size + pseudocount))
        #return(k_log(k_sum(c(mu,s))+pseudocount))
        #if(s < s_cutoff) { return(log(sum(rpois(sample_size, lambda=mu))/sample_size + pseudocount))
        #} else { return(log(sum(rpois(sample_size, lambda=mu))/sample_size + pseudocount)) }
    }) (c(background_mut_rate_module, selection_coef_layer))
    
    s_tensor <- dpois(lambda=expected)

    model <- keras_model(inputs=c(rbp_binding_input), outputs=c(rbp_disruption_module_output))
    opt <- optimizer_sgd(lr = 0.1, decay = 1e-2, momentum=0.5, nesterov=TRUE)
    masked_loss_function <- function(y_true, y_pred, mask=-1) { mask_vector <- k_cast(k_not_equal(y_true, mask), k_floatx());
                                                                return(k_binary_crossentropy(y_true * mask_vector, y_pred * mask_vector)) }
    # Compile model, and draw model network.
    model %>% compile(loss = c("mse"), optimizer=optimizer_rmsprop(), metrics="mean_absolute_error")
    model_name = "supermodel"
    kerasR::plot_model(model, to_file = output_path(paste0("images/",model_name,".png")), show_shapes = TRUE, show_layer_names = TRUE)
    plot_model(model)
    kerasR::modules$keras.utils$plot_model(model = model, to_file = output_path(paste0("images/",model_name,".png")), show_shapes = TRUE, show_layer_names = TRUE)
    

    #all_models <- lapply(models, function(model_path) {
    #    print(paste0("Loading model from ",model_path))
    #    return(load_model_hdf5(model_path))
    #})
    #shared_input <- all_models[[1]]$input
    #
    #all_model_sections <- lapply(all_models, function(model) {
    #    conv1 <- get_layer(model, "conv1")
    #    return(shared_input %>% model)
    #})
    #combined_model <- keras_model(inputs=c(shared_input), outputs=c(all_model_sections))

    # Fit model to training data
    randomized_order_training_indices = sample(1:nrow(dat_gradcams), floor(0.8*nrow(dat_gradcams)))
    test_indices = (1:nrow(dat_gradcams))[!((1:nrow(dat_gradcams)) %in% randomized_order_training_indices)]
    #history <- model %>% fit(list(array(DF, c(dim(DF), 1))[randomized_order_training_indices,,,,drop=FALSE], pLIs_dat[randomized_order_training_indices]), list(labels[randomized_order_training_indices]), epochs=10, batch_size=32, validation_split = 0.2)
    history <- model %>% fit(list(dat_gradcams[randomized_order_training_indices,,,,drop=FALSE], pLIs_dat[randomized_order_training_indices], expected[randomized_order_training_indices]), list(af_tensor[randomized_order_training_indices,]), epochs=10, batch_size=16, validation_split=0.2)
    print(history$metrics)
    filename = output_path("supermodel_training.pdf")
    pdf(filename)
    plot(1:length(history$metrics$mean_absolute_error), history$metrics$mean_absolute_error, type="l", lty=1, lwd=2, col="blue", ylim=range(c(history$metrics$mean_absolute_error, history$metrics$val_mean_absolute_error)), main="Supermodel Training", xlab="epoch", ylab="Mean Absolute Error", cex.lab=1.4, cex.axis=1.4, cex.main=1.3)
    lines(1:length(history$metrics$val_mean_absolute_error), history$metrics$val_mean_absolute_error, lty=1, lwd=2, col="red")
    legend("topleft", legend=c("train", "validation"), col=c("blue","red"), pch=15, cex=1.2)
    dev.off()
    pdf_to_png(filename)
    #model %>% save_model_hdf5(output_path(paste0(model_name,".h5")))
    #model <- load_model_hdf5(output_path(paste0(model_name,".h5")))
    
    # Use model to predict labels for new data!
    #pred_scores <- model %>% predict(c(list(array(DF, c(dim(DF), 1))[test_indices,,,,drop=FALSE], pLIs_dat[test_indices])))
    pred_scores <- model %>% predict(c(list(DF[test_indices,,,drop=FALSE], pLIs_dat[test_indices], expected[test_indices])))
    non_zero_observed_indices <- labels[test_indices] != min(labels[test_indices])
    cor(labels[test_indices], pred_scores, method="pearson")
    cor(labels[test_indices], pred_scores, method="spearman")
    cor(labels[test_indices][non_zero_observed_indices], pred_scores[non_zero_observed_indices], method="pearson")
    cor(labels[test_indices][non_zero_observed_indices], pred_scores[non_zero_observed_indices], method="spearman")
    mean(pred_scores[labels[test_indices] <= 0]); mean(pred_scores[labels[test_indices] > 0])
    unconstrained_pred_scores <- pred_scores[labels[test_indices] >= 0]
    constrained_pred_scores <- pred_scores[labels[test_indices] < 0]
    min_unconstrained_pred_score = min(unconstrained_pred_scores)
    local_min_unconstrained_pred_score = min(unconstrained_pred_scores[unconstrained_pred_scores > -0.01 & unconstrained_pred_scores < 0])
    dat[test_indices[which(pred_scores < min_unconstrained_pred_score)],]
    unconstrained_dat <- dat[test_indices[labels[test_indices] >= 0],]
    constrained_dat <- dat[test_indices[labels[test_indices] < 0],]
    nrow(unconstrained_dat)
    sum(unconstrained_dat$sample != "random")
    sum(unconstrained_dat$sample != "random") / nrow(unconstrained_dat)
    nrow(constrained_dat)
    sum(constrained_dat$sample != "random")
    sum(constrained_dat$sample != "random") / nrow(constrained_dat)

    unconstrained_dat <- dat[test_indices[pred_scores >= 0],]
    constrained_dat <- dat[test_indices[pred_scores < 0],]
    nrow(unconstrained_dat)
    sum(unconstrained_dat$sample != "random")
    sum(unconstrained_dat$sample != "random") / nrow(unconstrained_dat)
    nrow(constrained_dat)
    sum(constrained_dat$sample != "random")
    sum(constrained_dat$sample != "random") / nrow(constrained_dat)

    filename = output_path("pred_score_distributions.pdf")
    pdf(filename)
    plot(density(pred_scores[labels[test_indices] < 0]), type="l", lwd=2, col="red", cex.axis=1.4, cex.lab=1.4, cex.main=1.3, main="AF_lognorm Prediction Distributions", xlab="normalized ln(AF+pseudocount)")
    lines(density(pred_scores[labels[test_indices] >= 0], from=min_unconstrained_pred_score), lwd=2, col="blue")
    abline(v=min_unconstrained_pred_score, lty=2)
    #abline(v=local_min_unconstrained_pred_score, lty=2)
    #mtext(sum(pred_scores < min_unconstrained_pred_score), cex=1.2)
    legend("topleft", legend=c("labels < 0", "labels >= 0"), col=c("red", "blue"), pch=15, cex=1.3)
    dev.off()
    pdf_to_png(filename)

    # Multiple-threshold Enrichment Test
    multi_threshold_test <- function(cases, controls, thresholds=seq(0, 1, by=0.1), threshold_dir=">", feature_name=NULL, label="") {
        feature = which(colnames(cases) == feature_name)[1]; if(is.na(feature)) { feature = 1 }
        threshold_dir = 2*((threshold_dir == ">") - 0.5)
        greater_than_threshold_fets <- unfactorize(data.frame(t(sapply(thresholds, function(threshold) {
            m1 = sum((threshold_dir*cases[,feature]) >= (threshold_dir*threshold))
            n1 = nrow(cases)
            m0 = sum(threshold_dir*(controls[,feature]) >= (threshold_dir*threshold))
            n0 = nrow(controls)
            fet_result <- fisher.test(matrix(c(m1, n1-m1, m0, n0-m0), nrow = 2, dimnames = list(hits = c("Y", "N"), status = c("case", "control"))), alternative = "two.sided")
            return(c(threshold, fet_result$estimate, fet_result$conf.int, fet_result$p.value, m1, m0, n1, n0)) 
        })))); colnames(greater_than_threshold_fets) <- c("threshold", "estimate", "conf.int_lower", "conf.int_higher", "p.value", "m1", "m0", "n1", "n0")
        greater_than_threshold_fets[greater_than_threshold_fets == Inf] <- 0
        greater_than_threshold_fets[,c("estimate", "conf.int_lower", "conf.int_higher")] <- log2(greater_than_threshold_fets[,c("estimate", "conf.int_lower", "conf.int_higher")])
        greater_than_threshold_fets[greater_than_threshold_fets == -Inf] <- 0
        colnames(greater_than_threshold_fets) <- c("threshold", paste0(label,":", c("estimate", "conf.int_lower", "conf.int_higher", "p.value", "m1", "m0", "n1", "n0")))
        return(greater_than_threshold_fets)
    }
    split_points <- seq(-0.5, 0, by=0.1)
    pred_score_enrichments <- lapply(split_points, function(split_point) {
        enrichment_test <- multi_threshold_test(pred_scores[labels[test_indices] < split_point,,drop=FALSE], pred_scores[labels[test_indices] >= split_point,,drop=FALSE], seq(-0.05, 0.05, by=0.005), threshold_dir="<", label=paste0("split=",split_point))
        filename = output_path(paste0("pred_score_distributions_",split_point,"split.pdf"))
        pdf(filename)
        plot(density(pred_scores[labels[test_indices] < split_point]), type="l", lwd=2, col="red", cex.axis=1.4, cex.lab=1.4, cex.main=1.3, main="AF_lognorm Prediction Distributions", xlab="normalized ln(AF+pseudocount)")
        lines(density(pred_scores[labels[test_indices] >= split_point], from=min_unconstrained_pred_score), lwd=2, col="blue")
        optimal_threshold = enrichment_test$threshold[which.min(enrichment_test[,which(grepl(":p.value", colnames(enrichment_test)))])]
        abline(v=optimal_threshold, lty=2)
        mtext(paste0("Optimal constraint threshold: ",optimal_threshold), cex=1.2)
        legend("topleft", legend=c(paste0("labels < ",split_point), paste0("labels >= ",split_point)), col=c("red", "blue"), pch=15, cex=1.3)
        dev.off()
        pdf_to_png(filename)
        return(enrichment_test)
    }) %>% reduce(full_join, by="threshold")
    pred_score_enrichments[1:5,1:9]
    apply(pred_score_enrichments[,c(1,which(grepl(":p.value", colnames(pred_score_enrichments))))], 2, min)
    pred_score_enrichments[,c(1,which(grepl(":estimate", colnames(pred_score_enrichments))))]
    # Print optimal cutoff and p.value for each split_point value.
    t(apply(pred_score_enrichments[,which(grepl(":p.value", colnames(pred_score_enrichments)))], 2, function(x) { optimal = which.min(x); return(c(pred_score_enrichments[optimal,c(1)], min(x))) }))
    write.csv(pred_score_enrichments, file=output_path(paste0("pred_scores_enrichments.csv")), row.names=FALSE)
    #write.csv(pred_score_enrichments[,c(1,which(grepl(":estimate", colnames(greater_than_threshold_fets_ss))),which(grepl(":m1", colnames(greater_than_threshold_fets_ss))))], file=output_path(paste0(rbp,"_gene_expression_enrichment_estimates_vs_secondary_structure.csv")), row.names=FALSE)
    
    filename = output_path(paste0("pred_score_enrichments.pdf"))
    ci_width = 0.001
    mtext_label = paste0(nrow(pred_scores)," test variants")
    cols <- rainbow(length(split_points))
    sapply_out <- sapply(1:length(split_points), function(ss_i) { #1:length(split_points)
        ss <- paste0("split=",split_points[ss_i])
        col <- cols[ss_i]
        if(ss_i == 1) {
            pdf(filename)
            plot(pred_score_enrichments$threshold, pred_score_enrichments[,paste0(ss,":estimate")], main=paste0("RBP Binding Enrichment (high vs. low expressed)"), xlab="Gap between high/low expression (ln(TPM+1), around median)", ylab="log2(Odds Ratio of RBP binding)", col=col, type="l", lwd=2, pch=19, xaxs="i", yaxs="i", ylim=c(min(c(0,min(unlist(pred_score_enrichments[,grepl(":conf.int_lower",colnames(pred_score_enrichments))])))), max(unlist(pred_score_enrichments[,grepl(":conf.int_higher",colnames(pred_score_enrichments))]))), cex.axis=1.4, cex.lab=1.4, cex.main=1.4)
            #for(quant in quantile(DF$gene_expression)[2:4]) { abline(v=quant, col="black", lty=3) }
            abline(h=0, col="black", lty=1)
            mtext(mtext_label, cex=1.1)
        } else {
            lines(pred_score_enrichments$threshold, pred_score_enrichments[,paste0(ss,":estimate")], col=cols[ss_i], type="l", lwd=2, pch=19)
        }
        ci_col <- adjustcolor(col, alpha.f=0.3)
        segments(x0=pred_score_enrichments$threshold, y0=pred_score_enrichments[,paste0(ss,":conf.int_lower")], x1=pred_score_enrichments$threshold, y1=pred_score_enrichments[,paste0(ss,":conf.int_higher")], col=ci_col)
        segments(x0=pred_score_enrichments$threshold-ci_width, y0=pred_score_enrichments[,paste0(ss,":conf.int_lower")], x1=pred_score_enrichments$threshold+ci_width, y1=pred_score_enrichments[,paste0(ss,":conf.int_lower")], col=ci_col)
        segments(x0=pred_score_enrichments$threshold-ci_width, y0=pred_score_enrichments[,paste0(ss,":conf.int_higher")], x1=pred_score_enrichments$threshold+ci_width, y1=pred_score_enrichments[,paste0(ss,":conf.int_higher")], col=ci_col)
    })
    legend("bottomleft", legend=c(split_points,"quantiles"), col=c(cols,"black"), pch=c(rep(19,length(split_points)),NA), lty=c(rep(NA,length(split_points)),3), cex=1.05)
    dev.off()
    pdf_to_png(filename) 
    #plot(density(DF_control$gene_expression, from=0), main=paste0("ln(TPM+1) Gene Expression Distributions"), xlab="Gene expression: ln(TPM+1)", lty=1, lwd=2, col="blue", xaxs="i", yaxs="i", cex.axis=1.4, cex.lab=1.4, cex.main=1.4)
    #mtext(mtext_label, cex=1.1)
    #lines(density(DF_eclip$gene_expression, from=0), lty=1, lwd=2, col="red")
    #for(quant in quantile(DF$gene_expression)[2:4]) { abline(v=quant, col="black", lty=3) }
    #legend("topright", legend=c("eCLIP", "control", "quantiles"), col=c("red","blue","black"), lty=c(1,1,3), cex=1.1)
    #filename = output_path(paste0(rbp,"_gene_expression_distributions.pdf"))
    #dev.copy2pdf(file=filename)
    #pdf_to_png(filename)
            

    pdf(output_path("preds_vs_labels.pdf"))
    plot(labels[test_indices], pred_scores)
    dev.off()
    
    # Directly look at weights in the trained model to determine learned motifs.
    filters <- t(get_weights(model$get_layer("conv1"))[[1]][,,])
    rownames(filters) <- paste0(c(1:nrow(filters))) #colnames(filters) <- c("REF", "ALT")
    #filters
    filter_cols <- c("green","orange")[apply(filters, 1, function(x) x[2] < x[1]) + 1]
    filename = output_path("rbp_disruption_filters.pdf")
    pdf(filename)
    print(Heatmap(filters, 
            show_heatmap_legend = TRUE, name = "weight", #title of legend
            row_title = "ref_binding_score : alt_binding_score interaction filters", column_title = "REF                                         ALT",
            cluster_rows=TRUE, cluster_columns=FALSE,
            row_dend_side="left", column_dend_side="bottom",
            row_names_side="left", column_names_side="top", #column_names_rot=0,
            row_names_gp = gpar(col=filter_cols, fontsize = 10), column_names_gp = gpar(fontsize = 20) # Text size for row and column names
    ))
    dev.off()
    pdf_to_png(filename)

    rbp_matrix_output_all_weights <- get_weights(model$get_layer("dense1"))[[1]]
    dim(rbp_matrix_output_all_weights)
    rbp_matrix_output_weights <- matrix(rowSums(rbp_matrix_output_all_weights), nrow=num_filters, ncol=num_rbps)
    rownames(rbp_matrix_output_weights) <- paste0(c(1:nrow(rbp_matrix_output_weights)))
    filename = output_path("rbp_disruption_weights.pdf")
    pdf(filename)
    print(Heatmap(rbp_matrix_output_weights, 
            show_heatmap_legend = TRUE, name = "weight", #title of legend
            row_title = "ref_binding_score : alt_binding_score interaction filters", column_title = "RBP",
            cluster_rows=TRUE, cluster_columns=TRUE,
            row_dend_side="left", column_dend_side="top",
            row_names_side="left", column_names_side="top",
            row_names_gp = gpar(col=filter_cols, fontsize = 10), column_names_gp = gpar(fontsize = 5) # Text size for row and column names
    ))
    dev.off()
    pdf_to_png(filename)


    
    # CNN + AdaBoost model using HGMD and gnomAD
    library("caret")
    library("xgboost")
    train_indices_vec <- sample(1:nrow(DF_gnomad), floor(nrow(DF_gnomad)*0.8))
    train_indices <- rep(FALSE, nrow(DF_gnomad)); train_indices[train_indices_vec] <- TRUE
    model <- train(harmful ~., data = data.frame(DF_gnomad[train_indices,]), method = "xgbTree", trControl = trainControl("cv", number = 10))
    model$bestTune
    predictions <- model %>% predict(data.frame(DF_gnomad[!train_indices,]))
    cor(predictions, DF_gnomad[!train_indices,"harmful"])
    RMSE(predictions, DF_gnomad[!train_indices,"harmful"]) # Compute the average prediction error RMSE
    varImp(model) # Rank variable importance
    pdf(output_path("xgboost_performance.pdf"))
    plot(DF_gnomad[!train_indices,"harmful"], predictions, main="", xlab="", ylab="", cex.axis=1.4, cex.lab=1.4, cex.main=1.3)
    mtext()
    dev.off()
    
    draw_plot(data.frame(x=DF_gnomad[!train_indices,"harmful"], y=predictions), hex_density = 25, title="gnomAD XGBoost Result", xlab="log(obs/exp)", ylab="Predicted log(obs/exp)", legend_text="# variants", filename="gnomAD_XGBoost_result.pdf")
}

supermodel(standardize_colnames(read.csv(file="../WGS/output/regulatory_variants.tsv", sep="\t")))


#################################################################################################################
# gnomAD-only
#################################################################################################################

gnomad <- read.csv(data_path("random_gnomad_variants.txt"), sep="\t", header=FALSE)
colnames(gnomad)[1:5] <- c("Chrom", "Position", "sample", "Ref", "Alt")
observed <- lapply(strsplit(paste0(gnomad$V8), ";"), function(x) { as.numeric(strsplit(x[which(grepl("AF=",x))[1]],"[=,]")[[1]][-1]) } )


expected <- get_expected_muts(gnomad, width=1, version="hg19")
gnomad <- cbind(gnomad[,c("Chrom", "Position", "Ref", "Alt", "sample")], observed, expected, observed/expected)
gnomad <- gnomad[gnomad$expected > 0 & !is.na(gnomad$observed),]
gnomad <- add_sequence_context_feature(gnomad, width=151, version="hg19")
gnomad_refs_tensor <- genomic_sequences_to_tensor(gnomad$ref_sequence, sequence_length=151)
gnomad_alts_tensor <- genomic_sequences_to_tensor(gnomad$alt_sequence, sequence_length=151)

prediction_score_to_likelihood_mapping_table <- read.csv(output_path("prediction_score_to_likelihood_mapping_table.csv"))
calculate_scores("gnomad_scores", gnomad_refs_tensor, gnomad_alts_tensor)
gnomad_scores <- get_scores("gnomad_scores", c("ref_pred_score", "alt_pred_score"))

# load genebody
genebody <- read.csv(data_path("refGene\\refGene_hg19_genebody_fixed.bed"), sep="\t", stringsAsFactors=FALSE)
genebody_granges <- to_genomic_regions(genebody, chr_colname="chromosome", start_colname="tss", end_colname="tes", strand_colname="strand", label_colname="gene", order_coordinates=TRUE, remove_duplicate_labels=TRUE)
# Load pLIs
exac_dat <- read.csv(data_path("fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt"), sep="\t")

# Set up DF and gradient boosting super-model, this time for regression task.
DF_gnomad <- cbind(log(gnomad$observed.expected/15708+1), gnomad_scores[["ref_pred_score"]], gnomad_scores[["alt_pred_score"]]) #data.frame(as.numeric(!control_indices), ref_pred_scores, alt_pred_scores);
colnames(DF_gnomad) <- c("harmful", gsub("$","_ref",colnames(gnomad_scores[["ref_pred_score"]])), gsub("$","_alt",colnames(gnomad_scores[["alt_pred_score"]]))) #colnames(DF) <- c("harmful","ref_pred_score","alt_pred_score")
gnomad_variant_granges <- to_genomic_regions(gnomad, chr_colname="Chrom", start_colname="Position", end_colname="Position", label_colname="SampleID")
nearest_genes_gnomad <- names(genebody_granges)[nearest(gnomad_variant_granges, genebody_granges)]
#nearest_genes_gnomad <- cbind(nearest_genes_gnomad, b[unlist(sapply(nearest_genes_gnomad,function(x) { nearest_gene <- which(b$Gene.name == x); if(length(nearest_gene) < 1) { return(NA) } else { return(nearest_gene[1]) }  })),c("K562_TPM","HepG2_TPM",colnames(b)[(ncol(b) - 53):ncol(b)])]); colnames(nearest_genes_gnomad)[1:3] <- c("gene", "K562_logTPM", "HepG2_logTPM")
#for(i in 2:ncol(nearest_genes_gnomad)) {
#    nearest_genes_gnomad[,i][is.na(nearest_genes_gnomad[,i])] <- mean(nearest_genes_gnomad[,i][!is.na(nearest_genes_gnomad[,i])])
#    nearest_genes_gnomad[,i] <- log(nearest_genes_gnomad[,i] + 1)
#}
#DF_gnomad <- cbind(DF_gnomad,nearest_genes_gnomad[,-1])
pLIs_gnomad <- exac_dat[unlist(sapply(nearest_genes_gnomad,function(x) { nearest_gene <- which(paste0(exac_dat$gene) == x); if(length(nearest_gene) < 1) { return(NA) } else { return(nearest_gene[1]) }  })),"pLI"]
pLIs_gnomad[is.na(pLIs_gnomad)] <- mean(pLIs_gnomad[!is.na(pLIs_gnomad)])
DF_gnomad <- cbind(DF_gnomad, pLIs_gnomad)
colnames(DF_gnomad)[ncol(DF_gnomad)] <- "pLI"
DF_gnomad <- as.matrix(DF_gnomad)

pdf(output_path("gnomad_obs_exp_density.pdf"))
plot(density(exp(DF_gnomad[,"harmful"])))
dev.off()
pdf_to_png(output_path("gnomad_obs_exp_density.pdf"))
pdf(output_path("gnomad_obs_exp_density_ln.pdf"))
plot(density(log(DF_gnomad[,"harmful"])))
dev.off()
pdf_to_png(output_path("gnomad_obs_exp_density_ln.pdf"))

# CNN + AdaBoost model using HGMD and gnomAD
library("caret")
library("xgboost")
train_indices_vec <- sample(1:nrow(DF_gnomad), floor(nrow(DF_gnomad)*0.8))
train_indices <- rep(FALSE, nrow(DF_gnomad)); train_indices[train_indices_vec] <- TRUE
model <- train(harmful ~., data = data.frame(DF_gnomad[train_indices,]), method = "xgbTree", trControl = trainControl("cv", number = 10))
model$bestTune
predictions <- model %>% predict(data.frame(DF_gnomad[!train_indices,]))
cor(predictions, DF_gnomad[!train_indices,"harmful"])
RMSE(predictions, DF_gnomad[!train_indices,"harmful"]) # Compute the average prediction error RMSE
varImp(model) # Rank variable importance
pdf(output_path("xgboost_performance.pdf"))
plot(DF_gnomad[!train_indices,"harmful"], predictions, main="", xlab="", ylab="", cex.axis=1.4, cex.lab=1.4, cex.main=,1.3)
mtext()
dev.off()

draw_plot(data.frame(x=DF_gnomad[!train_indices,"harmful"], y=predictions), hex_density = 25, title="gnomAD XGBoost Result", xlab="log(obs/exp)", ylab="Predicted log(obs/exp)", legend_text="# variants", filename="gnomAD_XGBoost_result.pdf")


#a <- c(0,0,0,0)
#haha <- mclapply(1:4, function(i) a[i] <<- i*2, mc.cores=detectCores())

selected_chromosome = 17
chr7_refseq <- get_refseq(selected_chromosome, version="hg19", allow_BSgenome=TRUE)[[1]]
motif_length = 8
chr_length = length(chr7_refseq)
max_length_per_batch = 100000000
num_batches = ceiling(chr_length/max_length_per_batch)
length_per_batch = ceiling(chr_length/num_batches)

start = 1; end = chr_length #floor(length(chr7_refseq)*0.1) 
discoverable_ranges <- data.frame(selected_chromosome, start, end); colnames(discoverable_ranges) <- c("chromosome", "start", "end")
discoverable_ranges <- to_genomic_regions(discoverable_ranges)
h3k36me3_peaks <- load_annotation("E123.H3K36me3.fullPeak") #"E003.H3K36me3.broadPeak")
h3k36me3_peaks <- h3k36me3_peaks[unique(queryHits(findOverlaps(h3k36me3_peaks, discoverable_ranges)))]

rbps_to_focus_on <- c("K562.RBFOX2", "K562.EFTUD2", "K562.HNRNPU", "K562.ILF3", "K562.QKI")
lapply_out <- lapply(rbps_to_focus_on[-c(2)], function(rbp) {
    print(paste0("Finding genomic PPV for ",rbp))
    model_name = tolower(paste0(tolower(rbp),"_model2"))
    model <- load_model_hdf5(output_path(paste0(model_name,".h5")))

    eclip_peaks <- load_annotation(rbp)
    eclip_peaks <- eclip_peaks[unique(queryHits(findOverlaps(eclip_peaks, discoverable_ranges)))]
    eclip_peaks <- eclip_peaks[unique(queryHits(findOverlaps(eclip_peaks, h3k36me3_peaks)))]

    a <- lapply(1:num_batches, function(batch) {
        print(paste0("Batch ",batch," / ",num_batches))
        start = (batch-1)*length_per_batch + 1; end = min(c(batch*length_per_batch, chr_length))
        chr7_slices <- rollapply(start:end, width=151, by=151-(motif_length-1), FUN=function(x) { return(paste0(chr7_refseq[x])) }) #print(max(x)); 
        chr7_slices_ranges <- data.frame(cbind(rep(selected_chromosome,length(chr7_slices)), rollapply(start:end, width=151, by=151-(motif_length-1), FUN=range))); colnames(chr7_slices_ranges) <- c("chromosome", "start", "end")
        chr7_slices_ranges <- to_genomic_regions(chr7_slices_ranges)
        
        # LIMIT SLICES ONLY TO THE EXPRESSED ONES!
        expressed_indices <- unique(queryHits(findOverlaps(chr7_slices_ranges, h3k36me3_peaks)))
        chr7_slices_ranges <- chr7_slices_ranges[expressed_indices]
        chr7_slices <- chr7_slices[expressed_indices]
        
        chr7_tensor <- genomic_sequences_to_tensor(chr7_slices, sequence_length=151)
        binding_sites <- find_binding_sites(chr7_tensor, model, binding_sites=FALSE, model_name=gsub("model.*","chr",selected_chromosome,model_name))
        rm(chr7_tensor); gc()

        chr7_slices_with_eclip <- 1:length(chr7_slices_ranges) %in% queryHits(findOverlaps(chr7_slices_ranges, eclip_peaks))

        return(data.frame(binding_sites, chr7_slices_with_eclip))
    })
    a <- rbindlist(a)
    binding_sites <- unlist(a[,1])
    chr7_slices_with_eclip <- a[,2] == TRUE

    plot(density(binding_sites))

    length(eclip_peaks)
    length(chr7_slices) # All ranges
    sum(chr7_slices_with_eclip) # Number of Positives in real genomic data
    sum(!chr7_slices_with_eclip) # Number of Negatives in real genomic data
    ppv_result <- get_roc_result(pred_scores=binding_sites, labels=chr7_slices_with_eclip, plot_roc=FALSE, plot_precision_recall=TRUE, filename_prefix=paste0(rbp,"_chr7"), mtext_text=paste0(rbp," on chr7"))
    ppv_result <- data.frame(ppv_result[["cutoffs"]], ppv_result[["precision"]], ppv_result[["sensitivity"]]); colnames(ppv_result) <- c("cutoff", "precision", "recall"); ppv_result <- ppv_result[!is.nan(ppv_result$precision),]
    ppv_result <- unfactorize(cbind(ppv_result, data.frame(t(sapply(as.numeric(paste0(ppv_result$cutoff)), function(cutoff) {
        predicted_positives = sum(binding_sites > cutoff) # Number of predicted Positives
        TP = sum(binding_sites[chr7_slices_with_eclip] > cutoff) # Number of True Positives
        FP = sum(binding_sites[!chr7_slices_with_eclip] > cutoff) # Number of False Positives
        TN = sum(binding_sites[!chr7_slices_with_eclip] <= cutoff) # Number of True Negatives
        FN = sum(binding_sites[chr7_slices_with_eclip] <= cutoff) # Number of False Negatives
        #print(paste0("PPV with cutoff ",cutoff,": ", sum(binding_sites[chr7_slices_with_eclip] > cutoff) / sum(binding_sites > cutoff) ))
        return(c(TP, FP, TN, FN))
    }))))); colnames(ppv_result) <- c("cutoff", "precision", "recall", "TP", "FP", "TN", "FN")
    ppv_result[1:10,]
    write.csv(ppv_result, file=output_path(paste0(rbp,"_chr7_PPV.csv")), row.names=FALSE)

    return(ppv_result)
})
saveRDS(chr7_tensor, file = output_path("chr17_tensor.rds"))
saveRDS(chr7_slices, file = output_path("chr17_slices.rds"))
saveRDS(chr7_slices_ranges, file = output_path("chr17_slices_ranges.rds"))
saveRDS(chr7_slices_with_eclip, file = output_path("chr17_slices_with_eclip.rds"))
saveRDS(binding_sites, file = output_path("chr17_binding_sites.rds"))

read.csv(data_path("neg_rbp_seqs_H1.csv"))
print(load(file=data_path("expr_log2_medianTPM_plus1.rda")))
exp_median_gtex[1:10,]
exp_median_gtex[rownames(exp_median_gtex) == "GATA4",]

genes_k562_tpms <- read.csv(data_path("ENCFF047WAI_K562_gene_quantifications.tsv"), sep="\t"); genes_k562_tpms <- cbind(gsub("\\..*$","",genes_k562_tpms$gene_id), genes_k562_tpms); colnames(genes_k562_tpms)[1] <- "Gene.stable.ID"; colnames(genes_k562_tpms)[which(colnames(genes_k562_tpms) == "TPM")] <- "K562_TPM"
genes_hepg2_tpms <- read.csv(data_path("ENCFF945LNB_HepG2_gene_quantifications.tsv"), sep="\t"); genes_hepg2_tpms <- cbind(gsub("\\..*$","",genes_hepg2_tpms$gene_id), genes_hepg2_tpms); colnames(genes_hepg2_tpms)[1] <- "Gene.stable.ID"; colnames(genes_hepg2_tpms)[which(colnames(genes_hepg2_tpms) == "TPM")] <- "HepG2_TPM"
ensembl_dat <- read.csv(data_path("EnsemblID_genename.txt"), sep="\t")
a <- merge(ensembl_dat, genes_k562_tpms); a <- merge(a[,c("Gene.stable.ID","Gene.Start..bp.","Gene.End..bp.","Chromosome.scaffold.name","Gene.name","gene_id","transcript_id.s.","length","effective_length","expected_count","K562_TPM")], genes_hepg2_tpms[,c("Gene.stable.ID","HepG2_TPM")])

plot(log(a$K562_TPM + 1), log(a$HepG2_TPM + 1), main="HepG2 vs. K562 log(TPM+1)")
mtext(paste0("Spearman Corr: ",cor(log(a$K562_TPM + 1), log(a$HepG2_TPM + 1), method="spearman")))

b <- merge(a,gene_eclip_overlaps, by.x="Gene.name", by.y="gene")
plot(log(b$K562_TPM + 1), b$K562.RBFOX2_peaks_within_20kb, main="# RBFOX2 peaks within 20kb of gene vs. K562 log(TPM+1)")
mtext(paste0("Spearman Corr: ",cor(log(b$K562_TPM + 1), b$K562.RBFOX2_peaks_within_20kb, method="spearman")))


get_gene_eclip_overlaps <- function(rbp) {
    rbp_cell_line = c("E118","E123")[as.numeric(grepl("K562\\.",rbp))+1]
    eclip_peaks <- load_annotation(rbp)
    h3k36me3_peaks <- load_annotation(paste0(rbp_cell_line,".H3K36me3.fullPeak"))
    genebody_padded <- genebody_granges; start(genebody_padded) <- start(genebody_padded) - 20000; end(genebody_padded) <- end(genebody_padded) + 20000
    gene_eclip_overlaps <- unfactorize(data.frame(table(factor(queryHits(findOverlaps(genebody_padded, eclip_peaks)), levels=1:length(genebody_padded)))))
    colnames(gene_eclip_overlaps) <- c("gene", paste0(rbp,"_peaks_within_20kb"))
    gene_eclip_overlaps$gene <- names(genebody_padded)[gene_eclip_overlaps$gene]
    return(gene_eclip_overlaps)
}
gene_eclip_overlaps <- get_gene_eclip_overlaps("K562.RBFOX2")
gene_eclip_overlaps <- get_gene_eclip_overlaps("HepG2.RBFOX2")


DF <- data.frame(x=log(b$K562_TPM + 1), y=log(b$K562.RBFOX2_peaks_within_20kb + 1))

draw_plot <- function(DF, hex_density = 10, title="", xlab="log(TPM+1)", ylab="log(# peaks + 1)", legend_text="# genes", cor_method="Pearson", linear_best_fit=TRUE, quadratic_best_fit=TRUE, ignore=NA, filename=NULL) {
    p = ggplot(DF, aes(x=x, y=y)) + stat_binhex(bins=hex_density) + scale_fill_gradient(name=paste0(cor_method,"\nCorr:\n", round(cor(DF[,1], DF[,2], method=tolower(cor_method)), 3),paste0("\n\n",legend_text)), trans="log", breaks=10^(0:4)) + labs(title=title, x=xlab, y=ylab) + 
        theme(
            text = element_text(size=22),
            plot.title = element_text(size=16, face="bold", hjust=0.5, vjust=1),
            axis.title.x = element_text(size=16, face="bold", hjust=0.5, vjust=0.5),
            axis.title.y = element_text(size=16, face="bold", hjust=0.5, vjust=1),
            legend.title=element_text(size=16, face="bold"),
            legend.text=element_text(size=16)
        )
    if(linear_best_fit) { p = p + geom_smooth(method="lm", colour="green", se=TRUE) }
    if(quadratic_best_fit) { p = p + geom_smooth(method="lm", formula = y ~ x + I(x^2), colour="red", se=TRUE) }
    
    #+ annotation_custom(grob = textGrob(paste0("Pearson Corr:\n", round(cor(DF[,1], DF[,2], method="pearson"), 2))), xmin=max(DF$x), xmax=max(DF$x), ymin=max(DF$y), ymax=max(DF$y))
    
    print(p)
    
    ## Code to override clipping
    #gt <- ggplotGrob(p)
    #gt$layout$clip[gt$layout$name=="panel"] <- "off"
    #grid.draw(gt)
    
    if(!is.null(filename)) {
        ggsave(filename)
        pdf_to_png(filename)
    }
}



#table(unlist(ss_at_binding_site))/length(unlist(ss_at_binding_site))
#table(dominant_ss_at_binding_site)/length(dominant_ss_at_binding_site)
#table(dominant_ss_at_binding_site)
#table(DF$label[dominant_ss_at_binding_site == "internal_loop"])
#ss_at_binding_site[[which(dominant_ss_at_binding_site == "external_region")[1]]]

#DF_eclip <- DF[DF$label == 1,]
#DF_control <- DF[DF$label == 0,]
#p <- ggplot(DF_eclip, aes(x=dominant_ss_at_binding_site, y=gene_expression)) + geom_violin() + geom_boxplot(width=0.1)
#print(p)
#
#p <- ggplot(DF, aes(x=label, y=gene_expression)) + geom_violin() + geom_boxplot(width=0.1)
#print(p)
#
#p <- ggplot(DF, aes(x=dominant_ss_at_binding_site, y=gene_expression, color=label)) + geom_violin(width=0.5, position=position_dodge()) + geom_boxplot(width=0.05, binaxis='y', stackdir='center', position=position_dodge(0.5)) + scale_fill_manual(breaks=c('control','eCLIP'),values=c('red','darkgreen')) # + geom_jitter(alpha=0.3)
#print(p)
#
#
#ss_freq_at_binding_site <- t(data.frame(lapply(ss_at_binding_site, function(x) { empty <- rep(0, 5); names(empty) <- dimnames(secondary_structure)[[3]]; x_freqs <- table(x)/length(x); empty[names(x_freqs)] <- x_freqs; return(empty) })))
#rownames(ss_freq_at_binding_site) <- NULL
#DF <- data.frame(labels[,1], gene_expressions[,1], ss_freq_at_binding_site); colnames(DF)[1:2] <- c("label", "gene_expression"); DF$label <- as.factor(DF$label)
#DF_eclip <- DF[DF$label == 1,]
#DF_control <- DF[DF$label == 0,]
#p <- ggplot(DF, aes(x=label, y=paired)) + geom_violin() + geom_boxplot(width=0.1)
#print(p)
#cor(DF$paired, DF$label)
#p <- ggplot(DF, aes(x=label, y=hairpin_loop)) + geom_violin() + geom_boxplot(width=0.1)
#print(p)
#p <- ggplot(DF, aes(x=label, y=internal_loop)) + geom_violin() + geom_boxplot(width=0.1)
#print(p)
#p <- ggplot(DF, aes(x=label, y=multi_loop)) + geom_violin() + geom_boxplot(width=0.1)
#print(p)
#p <- ggplot(DF, aes(x=label, y=external_region)) + geom_violin() + geom_boxplot(width=0.1)
#print(p)
#
#p <- ggplot(DF, aes(x=dominant_ss_at_binding_site, y=gene_expression, color=label)) + geom_boxplot(width=0.5, position=position_dodge()) + scale_fill_manual(breaks=c('control','eCLIP'),values=c('red','darkgreen')) # + geom_jitter(alpha=0.3)
#print(p)

#filename = output_path("gene_expression_enrichment_vs_secondary_structure.pdf")
#pdf(file=filename)

# Return data frame with secondary structure, gene expression, and label information.
get_ss_DF <- function(secondary_structure, gene_expressions, labels=NULL, padding_to_remove=50, allow_paired_to_dominate=FALSE) {
    if(is.null(labels)) { labels <- data.frame(c(rep(1,nrow(secondary_structure)/2),rep(0,nrow(secondary_structure)/2))) }
    ss_at_binding_site <- lapply(1:nrow(secondary_structure), function(s) { 
        print(paste0(s," / ",nrow(secondary_structure)))
        s_length = sum(rowSums(secondary_structure[s,1:151,1:5,1])); s_start = padding_to_remove + 1; s_end = s_length-padding_to_remove
        ss <- dimnames(secondary_structure)[[3]][apply(secondary_structure[s,s_start:s_end,1:5,1], 1, which.max)]
        names(ss) <- dimnames(data)[[3]][apply(data[s,s_start:s_end,1:4,1], 1, which.max)]
        return(ss)
    })
    dominant_ss_at_binding_site <- unlist(lapply(ss_at_binding_site, function(x) { if(allow_paired_to_dominate && "paired" %in% x[min(c(length(x),2)):max(c((length(x)-1),1))]) { return("paired") } else { x_table <- table(x); return(names(x_table)[which.max(x_table)]) }  }))
    DF <- data.frame(dominant_ss_at_binding_site, gene_expressions[,1], labels[,1])
    colnames(DF)[2:3] <- c("gene_expression","label"); DF$label <- as.factor(DF$label)
    return(DF)
}
# Multiple-threshold Enrichment Test
library("tidyverse")
midpoint_expression = median(DF$gene_expression)
thresholds <- seq(0, midpoint_expression, by=0.1) #seq(0, ceiling(max(DF$gene_expression)), by=0.1)
secondary_structures <- c("any", "paired", "hairpin_loop", "internal_loop", "multi_loop", "external_region") #dimnames(secondary_structure)[[3]])
#rbp = "K562.RBFOX2"
rbps = c("K562.RBFOX2", "K562.QKI", "K562.EFTUD2", "K562.HNRNPU")
remove_0_expression = FALSE
for(rbp in rbps) {
    print(rbp)
    DF <- get_ss_DF(readRDS(output_path(paste0("secondary_structures/",rbp,"_ss.rds"))), read.csv(output_path(paste0("gene_expressions/",tolower(rbp),"_gene_expressions.csv"))))
    if(remove_0_expression) { DF <- DF[DF$gene_expression > 0,]; } #thresholds <- thresholds[-c(1)]
    DF_eclip <- DF[DF$label == 1,]
    DF_control <- DF[DF$label == 0,]
    greater_than_threshold_fets_ss <- lapply(secondary_structures, function(ss) {
        if(ss == "any") { DF_ss <- DF } else { DF_ss <- DF[DF$dominant_ss_at_binding_site == ss,] }
        #DF_ss_eclip <- DF_ss[DF_ss$label == 1,]; DF_ss_control <- DF_ss[DF_ss$label == 0,]
        greater_than_threshold_fets <- unfactorize(data.frame(t(sapply(thresholds, function(threshold) {
            high_expressed_indices <- DF_ss$gene_expression > midpoint_expression + (threshold/2)
            low_expressed_indices <- DF_ss$gene_expression < midpoint_expression - (threshold/2)
            DF_ss_high_expressed <- DF_ss[high_expressed_indices,]; DF_ss_low_expressed <- DF_ss[low_expressed_indices,]
            m1 = sum(DF_ss_high_expressed$label == 1) #sum(DF_ss_eclip$gene_expression > threshold)
            n1 = nrow(DF_ss_high_expressed) #nrow(DF_ss_eclip)
            m0 = sum(DF_ss_low_expressed$label == 1) #sum(DF_ss_control$gene_expression > threshold)
            n0 = nrow(DF_ss_low_expressed) #nrow(DF_ss_control)
            fet_result <- fisher.test(matrix(c(m1, n1-m1, m0, n0-m0), nrow = 2, dimnames = list(hits = c("Y", "N"), status = c("case", "control"))), alternative = "two.sided")
            return(c(threshold, fet_result$estimate, fet_result$conf.int, fet_result$p.value, m1, m0, n1, n0)) 
        })))); colnames(greater_than_threshold_fets) <- c("high_expression_threshold", "estimate", "conf.int_lower", "conf.int_higher", "p.value", "m1", "m0", "n1", "n0")
        greater_than_threshold_fets[greater_than_threshold_fets == Inf] <- 0
        greater_than_threshold_fets[,c("estimate", "conf.int_lower", "conf.int_higher")] <- log2(greater_than_threshold_fets[,c("estimate", "conf.int_lower", "conf.int_higher")])
        greater_than_threshold_fets[greater_than_threshold_fets == -Inf] <- 0
        colnames(greater_than_threshold_fets) <- c("high_expression_threshold", paste0(ss,":", c("estimate", "conf.int_lower", "conf.int_higher", "p.value", "m1", "m0", "n1", "n0")))
        return(greater_than_threshold_fets)
    }) %>% reduce(full_join, by="high_expression_threshold")
    
    greater_than_threshold_fets_ss[,c(1,which(grepl(":estimate", colnames(greater_than_threshold_fets_ss))))]
    write.csv(greater_than_threshold_fets_ss, file=output_path(paste0(rbp,"_gene_expression_enrichment_vs_secondary_structure.csv")), row.names=FALSE)
    write.csv(greater_than_threshold_fets_ss[,c(1,which(grepl(":estimate", colnames(greater_than_threshold_fets_ss))),which(grepl(":m1", colnames(greater_than_threshold_fets_ss))))], file=output_path(paste0(rbp,"_gene_expression_enrichment_estimates_vs_secondary_structure.csv")), row.names=FALSE)
    
    mtext_label = paste0(rbp,", ",nrow(DF_eclip)," eCLIP peaks")
    cols <- rainbow(length(secondary_structures))
    sapply_out <- sapply(1:length(secondary_structures), function(ss_i) {
        ss <- secondary_structures[ss_i]
        col <- cols[ss_i]
        if(ss_i == 1) {
            plot(greater_than_threshold_fets_ss$high_expression_threshold, greater_than_threshold_fets_ss[,paste0(ss,":estimate")], main=paste0("RBP Binding Enrichment (high vs. low expressed)"), xlab="Gap between high/low expression (ln(TPM+1), around median)", ylab="log2(Odds Ratio of RBP binding)", col=col, type="l", lwd=2, pch=19, xaxs="i", yaxs="i", ylim=c(min(c(0,min(unlist(greater_than_threshold_fets_ss[,grepl(":conf.int_lower",colnames(greater_than_threshold_fets_ss))])))), max(unlist(greater_than_threshold_fets_ss[,grepl(":conf.int_higher",colnames(greater_than_threshold_fets_ss))]))), cex.axis=1.4, cex.lab=1.4, cex.main=1.4)
            for(quant in quantile(DF$gene_expression)[2:4]) { abline(v=quant, col="black", lty=3) }
            abline(h=0, col="black", lty=1)
            mtext(mtext_label, cex=1.1)
        } else {
            lines(greater_than_threshold_fets_ss$high_expression_threshold, greater_than_threshold_fets_ss[,paste0(ss,":estimate")], col=cols[ss_i], type="l", lwd=2, pch=19)
        }
        ci_col <- adjustcolor(col, alpha.f=0.3)
        segments(x0=greater_than_threshold_fets_ss$high_expression_threshold, y0=greater_than_threshold_fets_ss[,paste0(ss,":conf.int_lower")], x1=greater_than_threshold_fets_ss$high_expression_threshold, y1=greater_than_threshold_fets_ss[,paste0(ss,":conf.int_higher")], col=ci_col)
        segments(x0=greater_than_threshold_fets_ss$high_expression_threshold-0.025, y0=greater_than_threshold_fets_ss[,paste0(ss,":conf.int_lower")], x1=greater_than_threshold_fets_ss$high_expression_threshold+0.025, y1=greater_than_threshold_fets_ss[,paste0(ss,":conf.int_lower")], col=ci_col)
        segments(x0=greater_than_threshold_fets_ss$high_expression_threshold-0.025, y0=greater_than_threshold_fets_ss[,paste0(ss,":conf.int_higher")], x1=greater_than_threshold_fets_ss$high_expression_threshold+0.025, y1=greater_than_threshold_fets_ss[,paste0(ss,":conf.int_higher")], col=ci_col)
    })
    legend("topleft", legend=c(secondary_structures,"quantiles"), col=c(cols,"black"), pch=c(rep(19,length(secondary_structures)),NA), lty=c(rep(NA,length(secondary_structures)),3), cex=1.05)
    filename = output_path(paste0(rbp,"_gene_expression_enrichments.pdf"))
    dev.copy2pdf(file=filename)
    pdf_to_png(filename)
    
    plot(density(DF_control$gene_expression, from=0), main=paste0("ln(TPM+1) Gene Expression Distributions"), xlab="Gene expression: ln(TPM+1)", lty=1, lwd=2, col="blue", xaxs="i", yaxs="i", cex.axis=1.4, cex.lab=1.4, cex.main=1.4)
    mtext(mtext_label, cex=1.1)
    lines(density(DF_eclip$gene_expression, from=0), lty=1, lwd=2, col="red")
    for(quant in quantile(DF$gene_expression)[2:4]) { abline(v=quant, col="black", lty=3) }
    legend("topright", legend=c("eCLIP", "control", "quantiles"), col=c("red","blue","black"), lty=c(1,1,3), cex=1.1)
    filename = output_path(paste0(rbp,"_gene_expression_distributions.pdf"))
    dev.copy2pdf(file=filename)
    pdf_to_png(filename)
}

if(grepl("GOF",rbp)) { more_extreme_than_threshold_fets <- more_extreme_than_threshold_fets[more_extreme_than_threshold_fets$threshold <= 0,]
} else if(grepl("LOF",rbp)) { more_extreme_than_threshold_fets <- more_extreme_than_threshold_fets[more_extreme_than_threshold_fets$threshold >= 0,] }
cols <- c("black","red")[as.numeric(more_extreme_than_threshold_fets$p.value < 0.05)+1]

plot(more_extreme_than_threshold_fets$threshold, more_extreme_than_threshold_fets$estimate, main=paste0("",rbp," ",variant_type," case enrichment"), xlab="RBP binding delta threshold", ylab="log2(Odds Ratio)", col=cols, pch=19, ylim=c(min(more_extreme_than_threshold_fets$conf.int_lower), max(more_extreme_than_threshold_fets$conf.int_higher)), cex.axis=1.4, cex.lab=1.4, cex.main=1.4)
abline(h=0, col="blue")
abline(v=0, col="blue")
#if(!(grepl("LOF|GOF",rbp))) { text(x=c(-0.25, 0.25), y=rep(min(more_extreme_than_threshold_fets$conf.int_lower), 2), labels=c("GOF", "LOF"), font=2, adj=0.5, cex=1.2) }
mtext("Gene expression more extreme than threshold", cex=1.2)
segments(x0=more_extreme_than_threshold_fets$threshold, y0=more_extreme_than_threshold_fets$conf.int_lower, x1=more_extreme_than_threshold_fets$threshold, y1=more_extreme_than_threshold_fets$conf.int_higher, col=cols)
segments(x0=more_extreme_than_threshold_fets$threshold-0.025, y0=more_extreme_than_threshold_fets$conf.int_lower, x1=more_extreme_than_threshold_fets$threshold+0.025, y1=more_extreme_than_threshold_fets$conf.int_lower, col=cols)
segments(x0=more_extreme_than_threshold_fets$threshold-0.025, y0=more_extreme_than_threshold_fets$conf.int_higher, x1=more_extreme_than_threshold_fets$threshold+0.025, y1=more_extreme_than_threshold_fets$conf.int_higher, col=cols)
# Draw case variants involved line on same plot
par(new = T)
plot(more_extreme_than_threshold_fets$threshold, more_extreme_than_threshold_fets$m1, type="l", lty=3, col="black", pch=16, axes=F, xlab=NA, ylab=NA, cex=1.4)
axis(side = 4)
mtext(side = 4, line = 3, "eCLIP sites involved", cex=1.4)
# Draw legend
#if(grepl("LOF",rbp)) { legend_location = "topright"
#} else { legend_location = "topleft" }
#legend(legend_location, legend=c("NS", "p<0.5", "count"), col=c("black","red","black"), pch=c(19,19,NA), lty=c(NA,NA,3), cex=1.1)

dev.off()
pdf_to_png(filename)




DF <- data.frame(x=gene_expressions[,1], y=rowMax(secondary_structure[,1:151,5,1]))
draw_plot(data.frame(x=log(b$K562_TPM + 1), y=log(b$K562.RBFOX2_peaks_within_20kb + 1)), title="# K562.RBFOX2 peaks within 20kb of gene vs. log(TPM+1)", filename=output_path("K562.RBFOX2_peaks_vs_gene_expression.pdf"))


# K562
b <- merge(a,get_gene_eclip_overlaps("K562.RBFOX2"), by.x="Gene.name", by.y="gene")
draw_plot(data.frame(x=log(b$K562_TPM + 1), y=log(b$K562.RBFOX2_peaks_within_20kb + 1)), title="# K562.RBFOX2 peaks within 20kb of gene vs. log(TPM+1)", filename=output_path("K562.RBFOX2_peaks_vs_gene_expression.pdf"))

draw_plot(data.frame(x=nearest_genes$K562_logTPM, y=DF$K562.RBFOX2_alt), title="K562.RBFOX2_ref binding predictions vs. log(TPM+1)", xlab="log(TPM+1)", ylab="binding score prediction", filename=output_path("K562.RBFOX2_ref_binding_vs_gene_expression.pdf"))
draw_plot(data.frame(x=nearest_genes$HepG2_logTPM, y=DF$HepG2.RBFOX2_alt), title="HepG2.RBFOX2_ref binding predictions vs. log(TPM+1)", xlab="log(TPM+1)", ylab="binding score prediction", filename=output_path("HepG2.RBFOX2_ref_binding_vs_gene_expression.pdf"))
draw_plot(data.frame(x=nearest_genes$K562_logTPM, y=DF$K562.RBFOX2_alt), title="K562.RBFOX2_alt binding predictions vs. log(TPM+1)", xlab="log(TPM+1)", ylab="binding score prediction", filename=output_path("K562.RBFOX2_alt_binding_vs_gene_expression.pdf"))
draw_plot(data.frame(x=nearest_genes$HepG2_logTPM, y=DF$HepG2.RBFOX2_alt), title="HepG2.RBFOX2_alt binding predictions vs. log(TPM+1)", xlab="log(TPM+1)", ylab="binding score prediction", filename=output_path("HepG2.RBFOX2_alt_binding_vs_gene_expression.pdf"))
draw_plot(data.frame(x=nearest_genes$K562_logTPM, y=DF$LOF_ref), title="K562.LOF_ref binding predictions vs. log(TPM+1)", xlab="log(TPM+1)", ylab="binding score prediction", filename=output_path("K562.LOF_ref_binding_vs_gene_expression.pdf"))
draw_plot(data.frame(x=nearest_genes$HepG2_logTPM, y=DF$LOF_ref), title="HepG2.LOF_ref binding predictions vs. log(TPM+1)", xlab="log(TPM+1)", ylab="binding score prediction", filename=output_path("HepG2.LOF_ref_binding_vs_gene_expression.pdf"))
draw_plot(data.frame(x=nearest_genes$K562_logTPM, y=DF$LOF_alt), title="K562.LOF_alt binding predictions vs. log(TPM+1)", xlab="log(TPM+1)", ylab="binding score prediction", filename=output_path("K562.LOF_alt_binding_vs_gene_expression.pdf"))
draw_plot(data.frame(x=nearest_genes$HepG2_logTPM, y=DF$LOF_alt), title="HepG2.LOF_alt binding predictions vs. log(TPM+1)", xlab="log(TPM+1)", ylab="binding score prediction", filename=output_path("HepG2.LOF_alt_binding_vs_gene_expression.pdf"))


#HepG2
b <- merge(a,get_gene_eclip_overlaps("HepG2.RBFOX2"), by.x="Gene.name", by.y="gene")
draw_plot(data.frame(x=log(b$HepG2_TPM + 1), y=log(b$HepG2.RBFOX2_peaks_within_20kb + 1)), title="# HepG2.RBFOX2 peaks within 20kb of gene vs. log(TPM+1)", filename=output_path("HepG2.RBFOX2_peaks_vs_gene_expression.pdf"))

# GTEX tissue TPM data
a <- read.table(data_path("GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct"), skip=2, header=3, sep="\t")
b <- merge(b, a, by.x="Gene.name", by.y="Description")


# Annotate nearest gene log(TPM+1)
genebody <- read.csv(data_path("refGene\\refGene_hg19_genebody_fixed.bed"), sep="\t", stringsAsFactors=FALSE)
genebody_granges <- to_genomic_regions(genebody, chr_colname="chromosome", start_colname="tss", end_colname="tes", strand_colname="strand", label_colname="gene", order_coordinates=TRUE, remove_duplicate_labels=TRUE)
regulatory_variant_granges <- to_genomic_regions(regulatory_variant_dat, chr_colname="Chrom", start_colname="Position", end_colname="Position", label_colname="sample")
nearest_genes <- names(genebody_granges)[nearest(regulatory_variant_granges, genebody_granges)]
nearest_genes <- cbind(nearest_genes, b[unlist(sapply(nearest_genes,function(x) { nearest_gene <- which(b$Gene.name == x); if(length(nearest_gene) < 1) { return(NA) } else { return(nearest_gene[1]) }  })),c("K562_TPM","HepG2_TPM",colnames(b)[(ncol(b) - 53):ncol(b)])]); colnames(nearest_genes)[1:3] <- c("gene", "K562_logTPM", "HepG2_logTPM")
for(i in 2:ncol(nearest_genes)) {
    nearest_genes[,i][is.na(nearest_genes[,i])] <- mean(nearest_genes[,i][!is.na(nearest_genes[,i])])
    nearest_genes[,i] <- log(nearest_genes[,i] + 1)
}
DF <- cbind(DF,nearest_genes[,-1])

nearest_genes$gene
exac_dat <- read.csv(data_path("fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt"), sep="\t")

pLIs <- exac_dat[unlist(sapply(nearest_genes$gene,function(x) { nearest_gene <- which(paste0(exac_dat$gene) == x); if(length(nearest_gene) < 1) { return(NA) } else { return(nearest_gene[1]) }  })),"pLI"]
pLIs[is.na(pLIs)] <- mean(pLIs[!is.na(pLIs)])

DF <- cbind(DF, pLIs)

# CNN + logistic regression model using HGMD and gnomAD
rbp = "K562.RBFOX2"
model <- load_model_hdf5(paste0("../ML/output/",tolower(rbp),"_model3.h5"))
ref_pred_scores <- model %>% predict(regulatory_variant_refs_tensor[,,,])
alt_pred_scores <- model %>% predict(regulatory_variant_alts_tensor[,,,])
train_indices <- sample(c(TRUE,FALSE),floor(0.8*length(control_indices)),replace=TRUE)
library("caret")
library("glmnet")
DF <- cbind(as.numeric(!control_indices), all_scores[["ref_pred_score"]], all_scores[["alt_pred_score"]]) #data.frame(as.numeric(!control_indices), ref_pred_scores, alt_pred_scores);
colnames(DF) <- c("harmful", gsub("$","_ref",colnames(all_scores[["ref_pred_score"]])), gsub("$","_alt",colnames(all_scores[["alt_pred_score"]]))) #colnames(DF) <- c("harmful","ref_pred_score","alt_pred_score")
DF <- as.matrix(DF)
logitMod_alpha = 1
logitMod_lambda_fit <- cv.glmnet(x=DF[train_indices,-1], y=DF[train_indices,1], family="binomial", alpha=logitMod_alpha, type.measure="mse", nfolds=10, parallel=TRUE)
logitMod_lambda = logitMod_lambda_fit$lambda.min
logitMod = glmnet(x=DF[train_indices,-1], y=DF[train_indices,1], family="binomial", alpha=logitMod_alpha, lambda=logitMod_lambda)
predictions_alpha1_expr <- c(predict(logitMod, newx=DF[!train_indices,-1], type="response"))
#logitMod <- glm(harmful ~ ., data=DF, subset=train_indices, family=binomial(link="logit"))
#predictions <- predict(logitMod, newdata=DF[!train_indices,], type="response")
summary(logitMod)
get_roc_result(predictions_alpha0_expr, DF[!train_indices,"harmful"], curve_names="alpha-0 elastic net log regression + expr", filename_prefix=output_path(paste0("HGMD_logistic_regression")))

# Lasso regression feature selection
a <- coef(logitMod)[-1,1] * apply(DF[,-1], 2, sd)
a <- a[a != 0]
a[order(abs(a), decreasing=TRUE)]
a[order(names(a))]
feature_counts <- sort(table(gsub("_(ref|alt)$","",names(a))), decreasing=TRUE)
important_disrupted_rbps <- a[names(a) %in% c(paste0(names(feature_counts)[feature_counts == 2],"_ref"),paste0(names(feature_counts)[feature_counts == 2],"_alt"))]
important_disrupted_rbps[]
important_gene_features <- a[!(grepl("_(ref|alt)$",names(a)))]
important_gene_features[]


# CNN + AdaBoost model using HGMD and gnomAD
library("fastAdaboost")
train_indices_vec <- sample(1:nrow(DF), floor(nrow(DF)*0.8))
train_indices <- rep(FALSE, nrow(DF)); train_indices[train_indices_vec] <- TRUE
adaboostMod <- adaboost(harmful ~ ., data.frame(DF[train_indices,]), 50)
predictions <- predict(adaboostMod, newdata=data.frame(DF[!train_indices,]))$prob[,2]
summary(adaboostMod)
get_roc_result(predictions, DF[!train_indices,"harmful"], filename_prefix=output_path(paste0("HGMD_AdaBoost")))

predictions_100trees_expr_pli <- predictions

get_roc_result(data.frame(predictions_150trees_expr_pli, predictions_100trees_expr_pli, predictions_50trees_expr_pli, predictions_alpha1_expr, predictions_alpha0_expr, predictions_25trees_expr, predictions_50trees_expr, predictions_alpha1, predictions_alpha05, predictions_alpha0, predictions_25trees, predictions_50trees, predictions_100trees), DF[!train_indices,"harmful"], curve_names=c("150-tree AdaBoost + expr + pLI", "100-tree AdaBoost + expr + pLI", "50-tree AdaBoost + expr + pLI", "alpha-1 elastic net log regression + expr", "alpha-0 elastic net log regression + expr", "25-tree AdaBoost + expr", "50-tree AdaBoost + expr", "alpha-1 elastic net log regression", "alpha-0.5 elastic net log regression", "alpha-0 elastic net log regression", "25-tree AdaBoost", "50-tree AdaBoost", "100-tree AdaBoost"), filename_prefix=output_path(paste0("HGMD_AdaBoost")), mtext_text="HGMD vs. gnomAD, using only ref/alt RBP binding predictions", legend.cex=0.8)
get_roc_result(data.frame(predictions_100trees_expr_pli, predictions_50trees_expr_pli, predictions_50trees_expr, predictions_alpha0_expr, predictions_50trees, predictions_alpha05), DF[!train_indices,"harmful"], curve_names=c("100-tree AdaBoost + expr + pLI", "50-tree AdaBoost + expr + pLI", "50-tree AdaBoost + expr", "elastic net log regression + expr", "50-tree AdaBoost", "elastic net log regression"), filename_prefix=output_path(paste0("HGMD_AdaBoost_ashg")), mtext_text="HGMD vs. gnomAD", legend.cex=1.1)

get_roc_result(data.frame(predictions_100trees_expr_pli), DF[!train_indices,"harmful"], curve_names="100-tree AdaBoost + expr + pLI", filename_prefix=output_path(paste0("HGMD_AdaBoost_100tree")), mtext_text="HGMD vs. gnomAD", legend.cex=1.1)
get_roc_result(data.frame(predictions_100trees_expr_pli), DF[!train_indices,"harmful"], filename_prefix=output_path(paste0("HGMD_AdaBoost_100tree")), mtext_text="HGMD vs. gnomAD", legend.cex=1.1)


pdf(output_path("HGMD_score_distributions.pdf"))
plot(density(ref_pred_scores[control_indices]), lwd=2, lty=1, col="blue", cex.axis=1.4, cex.lab=1.4, cex.main=1.3, main="RBP binding score distributions", xlab="Predicted binding score (from CNN)")
lines(density(alt_pred_scores[control_indices]), lwd=2, lty=2, col="blue")
lines(density(ref_pred_scores[!control_indices]), lwd=2, lty=1, col="red")
lines(density(alt_pred_scores[!control_indices]), lwd=2, lty=2, col="red")
legend("bottomright", legend=c("gnomAD ref","gnomAD alt","HGMD ref","HGMD alt"), col=c("blue","blue","red","red"), lty=c(1,2,1,2), cex=1.2)
dev.off()
pdf_to_png(output_path("HGMD_score_distributions.pdf"))

# ASD
DF_asd <- cbind(as.numeric(!asd_control_indices), asd_scores[["ref_pred_score"]], asd_scores[["alt_pred_score"]]) #data.frame(as.numeric(!control_indices), ref_pred_scores, alt_pred_scores);
colnames(DF_asd) <- c("harmful", gsub("$","_ref",colnames(asd_scores[["ref_pred_score"]])), gsub("$","_alt",colnames(asd_scores[["alt_pred_score"]]))) #colnames(DF) <- c("harmful","ref_pred_score","alt_pred_score")
asd_variant_granges <- to_genomic_regions(asd_variant_dat, chr_colname="Chrom", start_colname="Position", end_colname="Position", label_colname="SampleID")
nearest_genes_asd <- names(genebody_granges)[nearest(asd_variant_granges, genebody_granges)]
nearest_genes_asd <- cbind(nearest_genes_asd, b[unlist(sapply(nearest_genes_asd,function(x) { nearest_gene <- which(b$Gene.name == x); if(length(nearest_gene) < 1) { return(NA) } else { return(nearest_gene[1]) }  })),c("K562_TPM","HepG2_TPM",colnames(b)[(ncol(b) - 53):ncol(b)])]); colnames(nearest_genes_asd)[1:3] <- c("gene", "K562_logTPM", "HepG2_logTPM")
for(i in 2:ncol(nearest_genes_asd)) {
    nearest_genes_asd[,i][is.na(nearest_genes_asd[,i])] <- mean(nearest_genes_asd[,i][!is.na(nearest_genes_asd[,i])])
    nearest_genes_asd[,i] <- log(nearest_genes_asd[,i] + 1)
}
DF_asd <- cbind(DF_asd,nearest_genes_asd[,-1])
pLIs_asd <- exac_dat[unlist(sapply(nearest_genes_asd$gene,function(x) { nearest_gene <- which(paste0(exac_dat$gene) == x); if(length(nearest_gene) < 1) { return(NA) } else { return(nearest_gene[1]) }  })),"pLI"]
pLIs_asd[is.na(pLIs_asd)] <- mean(pLIs_asd[!is.na(pLIs_asd)])
DF_asd <- cbind(DF_asd, pLIs_asd)
colnames(DF_asd) <- colnames(DF)
DF_asd <- as.matrix(DF_asd)

predictions_asd <- predict(adaboostMod, newdata=data.frame(DF_asd))$prob[,2]
asd_roc_result <- get_roc_result(predictions_asd, DF_asd[,"harmful"], filename_prefix=output_path(paste0("ASD_AdaBoost")))

asd_control_predictions <- predictions_asd[asd_control_indices]
asd_case_predictions <- predictions_asd[!asd_control_indices]

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
TSS_regions <- get_TS_regions(genebody, sites="TSS", left=20000, right=20000); TES_regions <- get_TS_regions(genebody, sites="TES", left=5000, right=20000)
TSS_granges <- to_genomic_regions(TSS_regions, labels=TSS_regions$gene); TES_granges <- to_genomic_regions(TES_regions, labels=TES_regions$gene)
asd_TSS_indices <- 1:length(asd_variant_granges) %in% unique(queryHits(findOverlaps(asd_variant_granges, TSS_granges)))
asd_TES_indices <- 1:length(asd_variant_granges) %in% unique(queryHits(findOverlaps(asd_variant_granges, TES_granges)))

pdf(output_path("ASD_variant_score_distributions.pdf"))
plot(density(predictions_asd[control_indices]), lwd=2, lty=1, col="blue", cex.axis=1.4, cex.lab=1.4, cex.main=1.3, main="Variant RBP disruption score distributions", xlab="Predicted RBP disruption score")
lines(density(predictions_asd[!control_indices]), lwd=2, lty=1, col="red")
legend("topright", legend=c("ASD","control"), col=c("red","blue"), lty=c(1,1), cex=1.2)
dev.off()
pdf_to_png(output_path("ASD_variant_score_distributions.pdf"))
pdf(output_path("ASD_variant_score_distributions_righttail.pdf"))
plot(density(predictions_asd[control_indices]), lwd=2, lty=1, col="blue", xlim=c(0.4,1), ylim=c(0,2), cex.axis=1.4, cex.lab=1.4, cex.main=1.3, main="Variant RBP disruption score distributions", xlab="Predicted RBP disruption score")
lines(density(predictions_asd[!control_indices]), lwd=2, lty=1, col="red")
legend("topright", legend=c("ASD","control"), col=c("red","blue"), lty=c(1,1), cex=1.2)
dev.off()
pdf_to_png(output_path("ASD_variant_score_distributions_righttail.pdf"))

fet_result <- fisher_exact_test(sum(!asd_control_indices & asd_TES_indices), sum(asd_control_indices & asd_TES_indices), length(asd_case_predictions), length(asd_control_predictions), alternative=c("two.sided"))
fet_result$estimate
fet_result$p.value
fet_result <- fisher_exact_test(sum(!asd_control_indices & asd_TSS_indices), sum(asd_control_indices & asd_TSS_indices), length(asd_case_predictions), length(asd_control_predictions), alternative=c("two.sided"))
fet_result$estimate
fet_result$p.value

for(regulatory_region in c("3'UTR", "TSS", "3'UTR+TSS")) {
    if(regulatory_region == "3'UTR") { regulatory_region_indices <- asd_TES_indices
    } else if(regulatory_region == "TSS") { regulatory_region_indices <- asd_TSS_indices
    } else { regulatory_region_indices <- asd_TES_indices | asd_TSS_indices }
    fet_result <- fisher_exact_test(sum(!asd_control_indices & regulatory_region_indices), sum(asd_control_indices & regulatory_region_indices), length(asd_case_predictions), length(asd_control_predictions), alternative=c("two.sided"))
    cat("\n")
    cat(paste0(regulatory_region," overall enrichment: ",fet_result$estimate," (p=",fet_result$p.value,")"))
    cat("\n")
    for(cutoff in seq(0.5, 0.9, by=0.05)) {
        m1 = sum(predictions_asd[!asd_control_indices & regulatory_region_indices]>cutoff)
        m0 = sum(predictions_asd[asd_control_indices & regulatory_region_indices]>cutoff)
        n1 = length(asd_case_predictions)
        n0 = length(asd_control_predictions)
        fet_result <- fisher.test(matrix(c(m1, n1-m1, m0, n0-m0), nrow = 2, dimnames = list(hits = c("Y", "N"), status = c("case", "control"))), alternative = "two.sided")
        cat(paste0(regulatory_region," pathogenic variant (score>",cutoff,") enrichment: ",fet_result$estimate," (p=",fet_result$p.value,", ",m1,"/",n1," case variants, ",m0,"/",n0," control variants)"))
        cat("\n")
    }
}


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


cor(DF[,1], DF[,2], method="pearson")
lm(DF, formula = y ~ x + I(x^2))

rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')), bias=500)
r <- rf(32)
hb <- hexbin(cbind(log(b$K562_TPM + 1), b$K562.RBFOX2_peaks_within_20kb), xbins=hex_density)
plot(hb, colramp=rf, main=paste0("# RBFOX2 peaks within 20kb of gene vs. K562 log(TPM+1)"), xlab="log(TPM+1)", ylab="log(# peaks)")

dev.copy2pdf(output_path(paste0("")))
dev.off()

load_annotation("E123.H3K36me3.fullPeak")

rbp = "HepG2.RBFOX2" #"K562.RBFOX2"
rbp_cell_line = c("E118","E123")[as.numeric(grepl("K562\\.",rbp))+1]
eclip_peaks <- load_annotation(rbp)
h3k36me3_peaks <- load_annotation(paste0(rbp_cell_line,".H3K36me3.fullPeak"))
genebody_padded <- genebody_granges; start(genebody_padded) <- start(genebody_padded) - 20000; end(genebody_padded) <- end(genebody_padded) + 20000
gene_eclip_overlaps <- unfactorize(data.frame(table(factor(queryHits(findOverlaps(genebody_padded, eclip_peaks)), levels=1:length(genebody_padded)))))
colnames(gene_eclip_overlaps) <- c("gene", paste0(rbp,"_peaks_within_20kb"))
gene_eclip_overlaps$gene <- names(genebody_padded)[gene_eclip_overlaps$gene]

h3k36me3_peaks <- intersect(h3k36me3_peaks, genebody_padded, ignore.strand=TRUE)
gene_h3k36me3_overlaps <- unfactorize(data.frame(findOverlaps(genebody_padded, h3k36me3_peaks)))
padded_gene_lengths <- width(genebody_padded)
gene_h3k36me3_overlaps <- aggregate(gene_h3k36me3_overlaps$subjectHits, by=list(gene_h3k36me3_overlaps$queryHits), FUN=function(x){
    combined_peaks <- h3k36me3_peaks[x]
    combined_peaks <- intersect(combined_peaks, combined_peaks)
    return(sum(width(combined_peaks)))
})
colnames(gene_h3k36me3_overlaps) <- c("gene", "H3K36me3_coverage")
gene_h3k36me3_overlaps$H3K36me3_coverage <- gene_h3k36me3_overlaps$H3K36me3_coverage / padded_gene_lengths[gene_h3k36me3_overlaps$gene]
gene_h3k36me3_overlaps$gene <- names(genebody_padded)[gene_h3k36me3_overlaps$gene]

gene_fraction_covered_by_h3k36me3 <- sapply(1:nrow(gene_eclip_overlaps), function(i) {
    print(i)
    gene_h3k36me3_overlaps$queryHits == 
    return(sum(width(intersect(gene_h3k36me3_overlaps, genebody_padded[i], ignore.strand=TRUE))) / padded_gene_lengths[i])
})



rbps_to_focus_on <- get_features_by_group("RBP") #c("K562.RBFOX2", "K562.EFTUD2", "K562.HNRNPU")
#rbps_to_focus_on <- rbps_to_focus_on[gsub("^.+\\.","",rbps_to_focus_on) %in% ls(get_constrained_genes("pLI>0.5"))]

calculate_scores_tensor <- function(scores_name, variant_refs_tensor, variant_alts_tensor, rbps_to_focus_on=NULL, combined_model=TRUE) {
    #if(!grepl("^/",scores_folder)) { scores_folder = output_path(scores_folder) }
    #dir.create(scores_folder, showWarnings = FALSE)
    if(combined_model) { rbps_to_focus_on = "combined"
    } else if(is.null(rbps_to_focus_on)) { rbps_to_focus_on <- get_features_by_group("RBP") }

    counter = 1
    for(rbp in rbps_to_focus_on) {
        print(paste0(counter,". ",rbp))
        model <- load_model_hdf5(paste0("../ML/output/",tolower(rbp),"_model3.h5"))
        variant_refs_tensor_dims <- dim(variant_refs_tensor); variant_alts_tensor_dims <- dim(variant_alts_tensor)
        if(sum(variant_refs_tensor_dims != variant_alts_tensor_dims)>0) {
            print("ERROR: Ref and Alt Tensor dimensions do not match!")
            return(1)
        } else if(length(variant_refs_tensor_dims)<5) { # Handle positional, rather than regional, case.
            variant_refs_tensor <- array(variant_refs_tensor, dim=c(variant_refs_tensor_dims[1], 1, variant_refs_tensor_dims[-c(1)]))
            variant_alts_tensor <- array(variant_alts_tensor, dim=c(variant_alts_tensor_dims[1], 1, variant_alts_tensor_dims[-c(1)]))
        }
        rbp_features <- get_features_by_group("RBP")
        region_width = dim(variant_refs_tensor)[2]
        scores_tensor <- abind(lapply(1:region_width, function(region_pos) {
            print(region_pos)
            ref_pred_scores <- model %>% predict(variant_refs_tensor[,region_pos,,,])
            alt_pred_scores <- model %>% predict(variant_alts_tensor[,region_pos,,,])
        
            if(!is.list(ref_pred_scores)) { ref_pred_scores <- list(rbp=ref_pred_scores) }
            if(!is.list(alt_pred_scores)) { alt_pred_scores <- list(rbp=alt_pred_scores) }
            if(length(ref_pred_scores) == 160) { names(ref_pred_scores) <- rbp_features }
            if(length(alt_pred_scores) == 160) { names(alt_pred_scores) <- rbp_features }

            return(abind(data.frame(ref_pred_scores), data.frame(alt_pred_scores), along=3))
        }), along=4) #, mc.cores=12)
        scores_tensor <- aperm(scores_tensor, c(1,4,2,3))

        print("Writing scores tensor to file...")
        saveRDS(scores_tensor, file=output_path(paste0(scores_name,"_",rbp,"_scores_tensor.rds")))
 
        rm(model); gc()
        counter = counter + 1
        if(length(rbps_to_focus_on) == 1) { return(scores_tensor) }
    }
}

calculate_scores <- function(scores_folder, variant_refs_tensor, variant_alts_tensor, rbps_to_focus_on=NULL, combined_model=TRUE) {
    if(!grepl("^/",scores_folder)) { scores_folder = output_path(scores_folder) }
    dir.create(scores_folder, showWarnings = FALSE)
    if(combined_model) { rbps_to_focus_on = "combined"
    } else if(is.null(rbps_to_focus_on)) { rbps_to_focus_on <- get_features_by_group("RBP") }

    multiple_alts = class(variant_alts_tensor) == "list"

    counter = 1
    for(rbp in rbps_to_focus_on) {
        print(paste0(counter,". ",rbp))
        model <- load_model_hdf5(paste0("../ML/output/",tolower(rbp),"_model3.h5"))
        if(length(dim(ref_pred_scores))>4) {
            
        }
        ref_pred_scores <- model %>% predict(variant_refs_tensor[,,,])
        alt_pred_scores <- model %>% predict(variant_alts_tensor[,,,])
        if(!is.list(ref_pred_scores)) { ref_pred_scores <- list(rbp=ref_pred_scores) }
        if(!is.list(alt_pred_scores)) { alt_pred_scores <- list(rbp=alt_pred_scores) }
        if(length(ref_pred_scores) == 160) { names(ref_pred_scores) <- get_features_by_group("RBP") }
        if(length(alt_pred_scores) == 160) { names(alt_pred_scores) <- get_features_by_group("RBP") }
        rbps <- intersect(names(ref_pred_scores), names(alt_pred_scores))
        lapply_out <- lapply(rbps, function(rbp) {
            print(rbp)
            ref_pred_scores <- ref_pred_scores[[rbp]]; alt_pred_scores <- alt_pred_scores[[rbp]]
            ref_LR <- prediction_score_to_likelihood_mapping_table[round_to_nearest(ref_pred_scores/0.01)+1,rbp]
            alt_LR <- prediction_score_to_likelihood_mapping_table[round_to_nearest(alt_pred_scores/0.01)+1,rbp]
            delta_pred_scores <- ref_LR - alt_LR
            #delta_pred_scores <- -log(alt_pred_scores/ref_pred_scores)
            #delta_pred_scores[(alt_pred_scores < ref_pred_scores & (ref_pred_scores < 0.5 | alt_pred_scores > 0.5)) | (alt_pred_scores > ref_pred_scores & (ref_pred_scores > -0.5 | alt_pred_scores < -0.5))] <- 0
            return_dat <- cbind(ref_pred_scores, alt_pred_scores, ref_LR, alt_LR, delta_pred_scores)
            colnames(return_dat) <- c("ref_pred_score", "alt_pred_score", "ref_LR", "alt_LR", "delta_LR")
            write.csv(return_dat, file=paste0(scores_folder,"/",rbp,"_binding_scores.csv"), row.names=FALSE)
            #return(return_dat)
        })
        rm(model); gc()
        counter = counter + 1
    }
}
#saveRDS(variant_pred_scores, file = output_path("variant_pred_scores.rds"))
get_scores <- function(scores_folder, score_names=c("ref_pred_score", "alt_pred_score", "ref_LR", "alt_LR", "delta_LR"), rbps_to_focus_on=NULL) {
    if(!grepl("^/",scores_folder)) { scores_folder = output_path(scores_folder) }
    if(is.null(rbps_to_focus_on)) { rbps_to_focus_on <- get_features_by_group("RBP") }

    per_rbp_scores <- lapply(rbps_to_focus_on, function(rbp) {
        print(rbp)
        return(read.csv(paste0(scores_folder,"/",rbp,"_binding_scores.csv"))[,score_names])
    }); names(per_rbp_scores) <- rbps_to_focus_on
    all_scores <- new.env()
    for(i in 1:length(score_names)) {
        all_scores_i <- unfactorize(data.frame(sapply(rbps_to_focus_on, function(rbp) {
            return(per_rbp_scores[[rbp]][,score_names[i]])
        })))
        all_scores_i <- cbind(all_scores_i, t(apply(all_scores_i, 1, FUN=range)))
        colnames(all_scores_i)[(ncol(all_scores_i)-1):ncol(all_scores_i)] <- c("GOF", "LOF")
        all_scores[[score_names[i]]] <- all_scores_i
    }
    return(all_scores)
}
#disruptive_score_matrix_axis_names <- c("ref_LR", "alt_LR")
calculate_scores("HGMD_scores", regulatory_variant_refs_tensor, regulatory_variant_alts_tensor)
all_scores <- get_scores("HGMD_scores", c("ref_pred_score", "alt_pred_score"))

calculate_scores("ASD_scores", asd_variant_refs_tensor, asd_variant_alts_tensor, rbps_to_focus_on=get_features_by_group("RBP"))
asd_scores <- get_scores("ASD_scores", c("ref_pred_score", "alt_pred_score"))

calculate_scores("CHD_scores", chd_variant_refs_tensor, chd_variant_alts_tensor)
chd_scores <- get_scores("CHD_scores", c("ref_pred_score", "alt_pred_score"))

# CHD
DF_chd <- cbind(as.numeric(!chd_control_indices), chd_scores[["ref_pred_score"]], chd_scores[["alt_pred_score"]]) #data.frame(as.numeric(!control_indices), ref_pred_scores, alt_pred_scores);
colnames(DF_chd) <- c("harmful", gsub("$","_ref",colnames(chd_scores[["ref_pred_score"]])), gsub("$","_alt",colnames(chd_scores[["alt_pred_score"]]))) #colnames(DF) <- c("harmful","ref_pred_score","alt_pred_score")
chd_variant_granges <- to_genomic_regions(chd_variant_dat, chr_colname="Chrom", start_colname="Position", end_colname="Position", label_colname="SampleID")
nearest_genes_chd <- names(genebody_granges)[nearest(chd_variant_granges, genebody_granges)]
nearest_genes_chd <- cbind(nearest_genes_chd, b[unlist(sapply(nearest_genes_chd,function(x) { nearest_gene <- which(b$Gene.name == x); if(length(nearest_gene) < 1) { return(NA) } else { return(nearest_gene[1]) }  })),c("K562_TPM","HepG2_TPM",colnames(b)[(ncol(b) - 53):ncol(b)])]); colnames(nearest_genes_chd)[1:3] <- c("gene", "K562_logTPM", "HepG2_logTPM")
for(i in 2:ncol(nearest_genes_chd)) {
    nearest_genes_chd[,i][is.na(nearest_genes_chd[,i])] <- mean(nearest_genes_chd[,i][!is.na(nearest_genes_chd[,i])])
    nearest_genes_chd[,i] <- log(nearest_genes_chd[,i] + 1)
}
DF_chd <- cbind(DF_chd,nearest_genes_chd[,-1])
pLIs_chd <- exac_dat[unlist(sapply(nearest_genes_chd$gene,function(x) { nearest_gene <- which(paste0(exac_dat$gene) == x); if(length(nearest_gene) < 1) { return(NA) } else { return(nearest_gene[1]) }  })),"pLI"]
pLIs_chd[is.na(pLIs_chd)] <- mean(pLIs_chd[!is.na(pLIs_chd)])
DF_chd <- cbind(DF_chd, pLIs_chd)
colnames(DF_chd) <- colnames(DF)
DF_chd <- as.matrix(DF_chd)

predictions_chd <- predict(adaboostMod, newdata=data.frame(DF_chd))$prob[,2]
chd_roc_result <- get_roc_result(predictions_chd, DF_chd[,"harmful"], filename_prefix=output_path(paste0("CHD_AdaBoost")))

chd_control_predictions <- predictions_chd[chd_control_indices]
chd_case_predictions <- predictions_chd[!chd_control_indices]

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
TSS_regions <- get_TS_regions(genebody, sites="TSS", left=20000, right=20000); TES_regions <- get_TS_regions(genebody, sites="TES", left=5000, right=20000)
TSS_granges <- to_genomic_regions(TSS_regions, labels=TSS_regions$gene); TES_granges <- to_genomic_regions(TES_regions, labels=TES_regions$gene)
chd_TSS_indices <- 1:length(chd_variant_granges) %in% unique(queryHits(findOverlaps(chd_variant_granges, TSS_granges)))
chd_TES_indices <- 1:length(chd_variant_granges) %in% unique(queryHits(findOverlaps(chd_variant_granges, TES_granges)))

pdf(output_path("CHD_variant_score_distributions.pdf"))
plot(density(predictions_chd[control_indices]), lwd=2, lty=1, col="blue", cex.axis=1.4, cex.lab=1.4, cex.main=1.3, main="Variant RBP disruption score distributions", xlab="Predicted RBP disruption score")
lines(density(predictions_chd[!control_indices]), lwd=2, lty=1, col="red")
legend("topright", legend=c("CHD","SSC"), col=c("red","blue"), lty=c(1,1), cex=1.2)
dev.off()
pdf_to_png(output_path("CHD_variant_score_distributions.pdf"))
pdf(output_path("CHD_variant_score_distributions_righttail.pdf"))
plot(density(predictions_chd[control_indices]), lwd=2, lty=1, col="blue", xlim=c(0.4,1), ylim=c(0,2), cex.axis=1.4, cex.lab=1.4, cex.main=1.3, main="Variant RBP disruption score distributions", xlab="Predicted RBP disruption score")
lines(density(predictions_chd[!control_indices]), lwd=2, lty=1, col="red")
legend("topright", legend=c("CHD","SSC"), col=c("red","blue"), lty=c(1,1), cex=1.2)
dev.off()
pdf_to_png(output_path("CHD_variant_score_distributions_righttail.pdf"))

fet_result <- fisher_exact_test(sum(!chd_control_indices & chd_TES_indices), sum(chd_control_indices & chd_TES_indices), length(chd_case_predictions), length(chd_control_predictions), alternative=c("two.sided"))
fet_result$estimate
fet_result$p.value
fet_result <- fisher_exact_test(sum(!chd_control_indices & chd_TSS_indices), sum(chd_control_indices & chd_TSS_indices), length(chd_case_predictions), length(chd_control_predictions), alternative=c("two.sided"))
fet_result$estimate
fet_result$p.value

for(regulatory_region in c("3'UTR", "TSS", "3'UTR+TSS")) {
    if(regulatory_region == "3'UTR") { regulatory_region_indices <- chd_TES_indices
    } else if(regulatory_region == "TSS") { regulatory_region_indices <- chd_TSS_indices
    } else { regulatory_region_indices <- chd_TES_indices | chd_TSS_indices }
    fet_result <- fisher_exact_test(sum(!chd_control_indices & regulatory_region_indices), sum(chd_control_indices & regulatory_region_indices), length(chd_case_predictions), length(chd_control_predictions), alternative=c("two.sided"))
    cat("\n")
    cat(paste0(regulatory_region," overall enrichment: ",fet_result$estimate," (p=",fet_result$p.value,")"))
    cat("\n")
    for(cutoff in seq(0.5, 0.9, by=0.05)) {
        m1 = sum(predictions_chd[!chd_control_indices & regulatory_region_indices]>cutoff)
        m0 = sum(predictions_chd[chd_control_indices & regulatory_region_indices]>cutoff)
        n1 = length(chd_case_predictions)
        n0 = length(chd_control_predictions)
        fet_result <- fisher.test(matrix(c(m1, n1-m1, m0, n0-m0), nrow = 2, dimnames = list(hits = c("Y", "N"), status = c("case", "control"))), alternative = "two.sided")
        cat(paste0(regulatory_region," pathogenic variant (score>",cutoff,") enrichment: ",fet_result$estimate," (p=",fet_result$p.value,", ",m1,"/",n1," case variants, ",m0,"/",n0," control variants)"))
        cat("\n")
    }
}

#################################################################################################################################################
# Combine independently trained models, such as for different RBPs, into a single model with shared input that can be loaded and run much faster.
#################################################################################################################################################
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
combined_model <- combine_models(models=paste0("../ML/output/",tolower(get_features_by_group("RBP")),"_model4.h5"))
combined_model <- combine_models(models=paste0("../ML/output/",tolower(get_features_by_group("RBP")),"_model3.h5"))

a1 <- model1 %>% predict(regulatory_variant_alts_tensor[,,,])
a2 <- model2 %>% predict(regulatory_variant_alts_tensor[,,,])
a <- combined_model %>% predict(regulatory_variant_alts_tensor[,,,])

calculate_gradcams <- function(input_tensor, work_folder, motif_length = 8, motif_min_percent_distinct=0.5, motif_score_min_cutoff=0.1) {
    # This is the "african elephant" entry in the prediction vector
    #model_name = "K562.RBFOX2"; model <- load_model_hdf5(paste0("../ML/output/",tolower(model_name),"_model2.h5"))
    #model_name = "combined"
    #model <- load_model_hdf5(paste0("../ML/output/",tolower(model_name),"_model3.h5"))
    #model_layer_names <- sapply(model$layers, function(layer) paste0(layer)) #$name))
    
    #sequence_to_display = test_indices[sequence_index_to_display]
    #input_tensor <- list(adrop(data[sequence_to_display,,,,drop=FALSE], 4))
    dir.create(work_folder, showWarnings = FALSE)
    tf$compat$v1$disable_eager_execution()
    rbps <- get_features_by_group("RBP")
    rbp_heatmaps <- lapply(1:length(rbps), function(rbp_i) {
        rbp = rbps[rbp_i]
        print(paste0(rbp_i,". ",rbp))
        rbp_filename = full_path(work_folder,paste0(rbp,"_gnomad_gradcam.csv"))
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
            return(heatmap)

            # # Determine the optimal bound region of fixed motif_length from the heatmap. 
            # motif_start_index = which.max(rollapply(heatmap, width=motif_length, mean))
            # motif_end_index = motif_start_index + motif_length - 1
            # trimmed_binding_site <- paste0(sequence[motif_start_index:motif_end_index], collapse="")
            
            # heatmap_smoothed <- rollapply(heatmap, width=motif_length, mean, align="left", partial=TRUE)
            # motif_start_indices <- which(heatmap_smoothed >= motif_score_min_cutoff)
            # motif_start_indices <- motif_start_indices[order(heatmap_smoothed[motif_start_indices], decreasing=TRUE)]
            # get_best_nonoverlapping <- function(motif_start_indices) { if(length(motif_start_indices) > 0) { return(c(motif_start_indices[1], get_best_nonoverlapping(motif_start_indices[abs(motif_start_indices - motif_start_indices[1]) >= ceiling(motif_length*motif_min_percent_distinct)]))) } else { return(c()) } }
            # motif_start_indices <- get_best_nonoverlapping(motif_start_indices)
            # motif_end_indices <- motif_start_indices + motif_length - 1
            # trimmed_binding_sites <- sapply(1:length(motif_start_indices), function(motif_i) paste0(sequence[motif_start_indices[motif_i]:motif_end_indices[motif_i]], collapse=""))
            # motifs_dat <- unfactorize(data.frame(rep(paste0("seq",sequence_to_display),length(motif_start_indices)), rep(model_name,length(motif_start_indices)), rep(pred_scores[sequence_index_to_display],length(motif_start_indices)), motif_start_indices, motif_end_indices, trimmed_binding_sites, heatmap_smoothed[motif_start_indices], heatmap_smoothed[motif_start_indices]/heatmap_smoothed[motif_start_indices][1]))
            # colnames(motifs_dat) <- c("seq", "model", "pred_score", "start", "end", "motif", "score", "score_norm")
        }), along=2)
        write.csv(seq_heatmaps, file=rbp_filename)
        rm(seq_heatmaps)
        gc()
        return(0)
    })
    return(0)
    # #model$layers[[1]]$submodules
    # #submodules <- model$layers[[2]]$submodules
    # #submodule_types <- unlist(lapply(submodules, function(x) gsub(">","",gsub("^.*\\.([^\\.]*)$","\\1", paste0(x)))))

    # rbp_index = 1
    # rbp_output <- model$output[,1] 
    # #rbp_output <- model$output[[rbp_index]][,1]
    # # The is the output feature map of the `block5_conv3` layer,
    # # the last convolutional layer in VGG16
    # last_conv_layer <- model$get_layer("conv1") 
    # #last_conv_layer <- submodules[[min(which(grepl("Conv",submodule_types)))]]
    # # This is the gradient of the "african elephant" class with regard to
    # # the output feature map of `block5_conv3`

    # #tape <- tf$GradientTape()
    # #grads <- k_gradients(submodules[[13]], last_conv_layer$output)
    # grads <- k_gradients(rbp_output, last_conv_layer$output)[[1]]
    # #grads <- lapply(1:160, function(rbp_index) k_gradients(model$output[[rbp_index]][,1], last_conv_layer$output)[[1]])
    # if(is.null(grads)) { return(NULL) }
    # # This is a vector of shape (512,), where each entry
    # # is the mean intensity of the gradient over a specific feature map channel
    # pooled_grads <- k_mean(grads, axis=c(1,2))
    # # This function allows us to access the values of the quantities we just defined:
    # # `pooled_grads` and the output feature map of `block5_conv3`,
    # # given a sample image
    # iterate <- k_function(list(model$input), list(pooled_grads, last_conv_layer$output[1,,]))
    
    # #pred_scores[sequence_index_to_display]
    # sequence_to_display = test_indices[sequence_index_to_display]
    # iterate_result <- iterate(list(adrop(data[sequence_to_display,,,,drop=FALSE], 4)))
    # pooled_grads_value <- iterate_result[[1]]; conv_layer_output_value <- iterate_result[[2]]
    # # We multiply each channel in the feature map array
    # # by "how important this channel is" with regard to the elephant class
    # for (i in 1:100) {
    #     conv_layer_output_value[,i] <- conv_layer_output_value[,i] * pooled_grads_value[i]
    # }
    # # The channel-wise mean of the resulting feature map is our heatmap of class activation
    # #heatmap <- smooth.spline(apply(conv_layer_output_value, 1, mean))$y
    # heatmap <- smooth.spline(spline(apply(conv_layer_output_value, 1, mean), n=ncol(data), method="fmm")$y)$y
    # # Normalize the heatmap between 0 and 1, for better visualization.
    # heatmap <- pmax(heatmap, 0)
    # heatmap <- heatmap / max(heatmap)
    
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
    
    # # Determine the optimal bound region of fixed motif_length from the heatmap 
    # motif_start_index = which.max(rollapply(heatmap, width=motif_length, mean))
    # motif_end_index = motif_start_index + motif_length - 1
    # trimmed_binding_site <- paste0(sequence[motif_start_index:motif_end_index], collapse="")
    # text(0.5, 0.01, paste0("Optimal binding site: ",trimmed_binding_site), cex=1.7, font=2)
    
    # dev.off()
    # pdf_to_png(filename)
    
    # heatmap_smoothed <- rollapply(heatmap, width=motif_length, mean, align="left", partial=TRUE)
    # motif_start_indices <- which(heatmap_smoothed >= motif_score_min_cutoff)
    # motif_start_indices <- motif_start_indices[order(heatmap_smoothed[motif_start_indices], decreasing=TRUE)]
    # get_best_nonoverlapping <- function(motif_start_indices) { if(length(motif_start_indices) > 0) { return(c(motif_start_indices[1], get_best_nonoverlapping(motif_start_indices[abs(motif_start_indices - motif_start_indices[1]) >= ceiling(motif_length*motif_min_percent_distinct)]))) } else { return(c()) } }
    # motif_start_indices <- get_best_nonoverlapping(motif_start_indices)
    # motif_end_indices <- motif_start_indices + motif_length - 1
    # trimmed_binding_sites <- sapply(1:length(motif_start_indices), function(motif_i) paste0(sequence[motif_start_indices[motif_i]:motif_end_indices[motif_i]], collapse=""))
    
    # motifs_dat <- unfactorize(data.frame(rep(paste0("seq",sequence_to_display),length(motif_start_indices)), rep(model_name,length(motif_start_indices)), rep(pred_scores[sequence_index_to_display],length(motif_start_indices)), motif_start_indices, motif_end_indices, trimmed_binding_sites, heatmap_smoothed[motif_start_indices], heatmap_smoothed[motif_start_indices]/heatmap_smoothed[motif_start_indices][1]))
    # colnames(motifs_dat) <- c("seq", "model", "pred_score", "start", "end", "motif", "score", "score_norm")
    # return(motifs_dat)
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
        rbp_filename = full_path(work_folder,paste0(rbp,"_gnomad_gradcam.csv"))
        if(file.exists(rbp_filename)) { 
            if(subset_regions) { return(read.csv(rbp_filename, row.names=1)[,region_indices]) 
            } else { return(read.csv(rbp_filename, row.names=1)) }
        } else { return(NULL) }
    }), along=3)
    rbp_heatmaps <- aperm(rbp_heatmaps, c(2,1,3))
    return(rbp_heatmaps)
}

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
split_gradcams("asd", indices_groups=seq(40001, 126700, by=28900))

rbp 
 "LOF"
score_types <- c("pred_score", "LR")
delta_cutoffs <- new.env(); delta_cutoffs[["pred_score"]] <- seq(0, 0.25, by=0.05); delta_cutoffs[["LR"]] <- c(0, 0.5, 1, 2, 5, 15)
disruptive_score_bandwiths <- new.env(); disruptive_score_bandwiths[["pred_score"]] = 0.02; disruptive_score_bandwiths[["LR"]] = 1
max_disruptive_scores <- new.env(); max_disruptive_scores[["pred_score"]] = 1; max_disruptive_scores[["LR"]] = 50
for(score_type in score_types) {
    max_disruptive_score = max_disruptive_scores[[score_type]]
    disruptive_score_bandwith = disruptive_score_bandwiths[[score_type]]
    disruptive_scores <- seq(0, max_disruptive_score, by=disruptive_score_bandwith)
    if(score_type == "LR") { disruptive_scores[length(disruptive_scores)] <- paste0(">=",disruptive_scores[length(disruptive_scores)]) }
    
    disruptive_score_matrix_axis_names <- paste0(c("ref_","alt_"),score_type)
    ref_scores <- all_scores[[disruptive_score_matrix_axis_names[1]]][,rbp]
    alt_scores <- all_scores[[disruptive_score_matrix_axis_names[2]]][,rbp]
    delta_pred_scores <- ref_scores - alt_scores #-all_scores[["delta_LR"]][,rbp]
    
    for(delta_cutoff in delta_cutoffs[[score_type]]) {
        print(paste0("Drawing heatmap for ",rbp,"_",score_type," with cutoff ",delta_cutoff))
        disruptive_score_counts <- table(paste0(round_to_nearest(ref_scores[abs(delta_pred_scores) >= abs(delta_cutoff)],disruptive_score_bandwith),"->",round_to_nearest(alt_scores[abs(delta_pred_scores) >= abs(delta_cutoff)],disruptive_score_bandwith)))
        disruptive_score_counts_split <- unfactorize(data.frame(strsplit(names(disruptive_score_counts), "->")))
        if(nrow(disruptive_score_counts_split) < 2) { next }
        disruptive_score_counts_split[disruptive_score_counts_split > max_disruptive_score/disruptive_score_bandwith] <- max_disruptive_score/disruptive_score_bandwith
        disruptive_score_refs <- unlist(disruptive_score_counts_split[1,])
        disruptive_score_alts <- unlist(disruptive_score_counts_split[2,])
        disruptive_score_matrix <- matrix(rep(0, length(disruptive_scores)**2), nrow=length(disruptive_scores))
        #diag(length(disruptive_scores))-(2*diag(length(disruptive_scores)))
        rownames(disruptive_score_matrix) <- disruptive_scores; colnames(disruptive_score_matrix) <- disruptive_scores
        for(j in 1:length(disruptive_score_counts)) { 
            if(score_type == "LR") { disruptive_score_matrix[disruptive_score_refs[j], disruptive_score_alts[j]] <- disruptive_score_counts[j]  
            } else { disruptive_score_matrix[paste0(disruptive_score_refs[j]), paste0(disruptive_score_alts[j])] <- disruptive_score_counts[j] }
        }
        num_significant_variants = sum(disruptive_score_matrix)
        print(paste0("# significant variants: ",num_significant_variants))
        pdf(file=output_path(paste0(rbp,"_disruptive_score_matrix_ASD_",score_type,"_delta",delta_cutoff,"_",num_significant_variants,"_variants.pdf")))
        print(Heatmap(disruptive_score_matrix, 
                show_heatmap_legend = TRUE, name = "variants", #title of legend
                row_title = disruptive_score_matrix_axis_names[1], column_title = disruptive_score_matrix_axis_names[2],
                cluster_rows=FALSE, cluster_columns=FALSE
                ,row_dend_side="left", column_dend_side="top"
                ,row_names_gp = gpar(fontsize = 5), column_names_gp = gpar(fontsize = 5) # Text size for row and column names
        ))
        print("Done.")
        dev.off()
    }
}

a <- list.files(output_path()); a <- a[grepl("disruptive_score_matrix", a)]
for(a_i in a) { pdf_to_png(output_path(a_i)) }

#delta_pred_scores <- cbind(delta_pred_scores[,1:160], t(apply(delta_pred_scores[,gsub("^.+\\.","",colnames(delta_pred_scores)) %in% ls(get_constrained_genes("pLI>0.5"))], 1, FUN=range)))
delta_pred_scores <- cbind(delta_pred_scores[,1:160], t(apply(delta_pred_scores, 1, FUN=range)))
colnames(delta_pred_scores)[(ncol(delta_pred_scores)-1):ncol(delta_pred_scores)] <- c("GOF", "LOF")

saveRDS(delta_pred_scores, file=output_path("delta_pred_scores_ASD.rds"))

delta_pred_scores <- readRDS(file="../ML/output/delta_pred_scores_ASD.rds")

#delta_pred_scores <- apply(delta_pred_scores, 2, function(x) (x - mean(x)) / sd(x))
sapply_out <- sapply(1:ncol(delta_pred_scores), function(i) print(range(delta_pred_scores[,i])))
#filename = output_path(paste0("delta_prediction_score_densities.pdf"))
#pdf(file=filename)
cols <- rainbow(ncol(delta_pred_scores))
#plot_start = -0.9; plot_end = -0.5; y_max = 0.05
#plot(density(delta_pred_scores[,1], from=plot_start), col="white", main="Delta prediction score densities", lwd=2, lty=1, xlab="-log(alt_pred_scores/ref_pred_scores)", xlim=c(plot_start, plot_end), ylim=c(0, y_max), xaxs="i", cex.axis=1.4, cex.lab=1.4, cex.main=1.4) # xlim
sapply_out <- sapply(1:ncol(delta_pred_scores), function(i) print(range(delta_pred_scores[,i])))
sapply_out <- sapply(1:ncol(delta_pred_scores), function(i) print(range(delta_pred_scores[!control_indices,i])))
sapply_out <- sapply(1:ncol(delta_pred_scores), function(i) print(range(delta_pred_scores[control_indices,i])))
delta_pred_scores_cases_density <- density(unlist(delta_pred_scores[!control_indices,1:160]))
delta_pred_scores_controls_density <- density(unlist(delta_pred_scores[control_indices,1:160]))
pdf(file=output_path("delta_prediction_score_global_density.pdf"))
plot(delta_pred_scores_cases_density, col="red", lty=2, main="Delta prediction score global density")
lines(delta_pred_scores_controls_density, col="blue", lty=2)
legend("topright", legend=c("case", "control"), col=c("red", "blue"), pch=15)
dev.off()
pdf(file=output_path("delta_prediction_score_global_density_right_tail.pdf"))
plot(delta_pred_scores_cases_density, col="red", lty=2, main="Delta prediction score global density", xlim=c(1, 3.5), ylim=c(0,0.001))
lines(delta_pred_scores_controls_density, col="blue", lty=2)
legend("topright", legend=c("case", "control"), col=c("red", "blue"), pch=15)
dev.off()
pdf(file=output_path("delta_prediction_score_global_density_left_tail.pdf"))
plot(delta_pred_scores_cases_density, col="red", lty=2, main="Delta prediction score global density", xlim=c(-3, -1), ylim=c(0,0.001))
lines(delta_pred_scores_controls_density, col="blue", lty=2)
legend("topleft", legend=c("case", "control"), col=c("red", "blue"), pch=15)
dev.off()

variant_annotations <- annotate(regulatory_variant_dat, c("autism_genes_20000bp", "constrained_genes_20000bp", "H3K36me3"))[,-c(1:ncol(regulatory_variant_dat))] == "Y"
autism_gene_indices <- variant_annotations[,1]; constrained_gene_indices <- variant_annotations[,2]; H3K36me3_indices <- constrained_gene_indices <- variant_annotations[,3]

starts <- rep(0.25, ncol(delta_pred_scores)) #c(0.25, 0.25, 0.25)
bandwidths = rep(0.01, ncol(delta_pred_scores)) #c(0.01, 0.01, 0.002)
regulatory_variant_dat_indels <- regulatory_variant_dat$Type == "Indel" #nchar(paste0(regulatory_variant_dat$Ref)) > 1 | nchar(paste0(regulatory_variant_dat$Alt)) > 1 
variant_types <- c("SNV", "indel")
variant_constraints <- c("autism_gene", "constrained_gene", "H3K36me3", "")
sapply_out <- sapply(161:162, function(i) { #1:ncol(delta_pred_scores)
    for(variant_constraint in variant_constraints) {
        rbp = paste0("",colnames(delta_pred_scores)[i])
        if(variant_constraint == "") { constrained_variant_indices <- rep(TRUE, nrow(delta_pred_scores))
        } else { constrained_variant_indices <- get(paste0(variant_constraint,"_indices")); rbp = paste0(rbp," ",gsub("_gene","",variant_constraint)) }
        filename = output_path(paste0(gsub(" ","_",rbp),"_delta_prediction_score_enrichments.pdf"))
        pdf(file=filename, width=14)
        par(mfrow=c(1,2)) 
        par(mar=c(5.1,4.1,4.1,5.1))
        for(variant_type in variant_types) {
            print(paste0(c(rbp,variant_type,variant_constraint),collapse=", "))
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
            tail_width = 0.5
            
            #filename = output_path(paste0(gsub(" ","_",rbp),"_delta_prediction_score_densities.pdf"))
            #pdf(file=filename, width=14)
            #par(mfrow=c(1,2))
            
            #delta_pred_scores_cases_density <- density(delta_pred_scores_cases, bw=bandwidths[i])
            #delta_pred_scores_controls_density <- density(delta_pred_scores_controls, bw=bandwidths[i])
            
            #plot_start = min(delta_pred_scores[,i]); plot_end = -starts[i] #plot_start * (1-tail_width)
            ##delta_pred_scores_cases_density <- density(delta_pred_scores_cases, to=plot_end, bw=bandwidth)
            ##delta_pred_scores_controls_density <- density(delta_pred_scores_controls, to=plot_end, bw=bandwidth)
            ##print("Cases density: "); print(delta_pred_scores_cases_density)
            ##print("Controls density: "); print(delta_pred_scores_controls_density)
            #plot(delta_pred_scores_controls_density$x[delta_pred_scores_controls_density$x <= plot_end], delta_pred_scores_controls_density$y[delta_pred_scores_controls_density$x <= plot_end], col="blue", main=paste0(rbp," delta prediction score densities"), lwd=2, lty=1, xlab="-log(alt_pred_scores/ref_pred_scores)", ylab="Density", xlim=c(1.1*plot_start, plot_end), xaxs="i", cex.axis=1.4, cex.lab=1.4, cex.main=1.4, type="l")    
            #lines(delta_pred_scores_cases_density$x[delta_pred_scores_cases_density$x <= plot_end], delta_pred_scores_cases_density$y[delta_pred_scores_cases_density$x <= plot_end], col="red", lwd=2, lty=1)
            #legend("topleft", legend=c("case", "control"), col=c("red", "blue"), pch=15, cex=1.3)
            #mtext(paste0("score deltas below threshold ",plot_end), cex=1.2) #mtext(paste0("score deltas (lowest ",tail_width*100,"% tail)"), cex=1.2)
            
            #plot_end = max(delta_pred_scores[,i]); plot_start = starts[i] #plot_start = plot_end * (1-tail_width)
            ##delta_pred_scores_cases_density <- density(delta_pred_scores_cases, from=plot_start, bw=bandwidth)
            ##delta_pred_scores_controls_density <- density(delta_pred_scores_controls, from=plot_start, bw=bandwidth)
            ##print("Cases density: "); print(delta_pred_scores_cases_density)
            ##print("Controls density: "); print(delta_pred_scores_controls_density)
            #plot(delta_pred_scores_controls_density$x[delta_pred_scores_controls_density$x >= plot_start], delta_pred_scores_controls_density$y[delta_pred_scores_controls_density$x >= plot_start], col="blue", main=paste0(rbp," delta prediction score densities"), lwd=2, lty=1, xlab="-log(alt_pred_scores/ref_pred_scores)", ylab="Density", xlim=c(plot_start, 1.1*plot_end), xaxs="i", cex.axis=1.4, cex.lab=1.4, cex.main=1.4, type="l")
            #lines(delta_pred_scores_cases_density$x[delta_pred_scores_cases_density$x >= plot_start], delta_pred_scores_cases_density$y[delta_pred_scores_cases_density$x >= plot_start], col="red", lwd=2, lty=1)
            #legend("topright", legend=c("case", "control"), col=c("red", "blue"), pch=15, cex=1.3)
            #mtext(paste0("score deltas above threshold ",plot_start), cex=1.2) #mtext(paste0("score deltas (highest ",tail_width*100,"% tail)"), cex=1.2)
            #dev.off()
            #pdf_to_png(filename)
            #par(mfrow=c(1,1))
            
            thresholds <- seq(-150, 150, by=2)
            greater_than_threshold_fets <- unfactorize(data.frame(t(sapply(thresholds, function(threshold) {
                m1 = sum(delta_pred_scores_cases > threshold); n1 = length(delta_pred_scores_cases); m0 = sum(delta_pred_scores_controls > threshold); n0 = length(delta_pred_scores_controls)
                fet_result <- fisher.test(matrix(c(m1, n1-m1, m0, n0-m0), nrow = 2, dimnames = list(hits = c("Y", "N"), status = c("case", "control"))), alternative = "two.sided")
                return(c(threshold, fet_result$estimate, fet_result$conf.int, fet_result$p.value, m1, m0, n1, n0)) 
            })))); colnames(greater_than_threshold_fets) <- c("threshold", "estimate", "conf.int_lower", "conf.int_higher", "p.value", "m1", "m0", "n1", "n0")
            greater_than_threshold_fets[greater_than_threshold_fets == Inf] <- 0
            greater_than_threshold_fets[,c("estimate", "conf.int_lower", "conf.int_higher")] <- log2(greater_than_threshold_fets[,c("estimate", "conf.int_lower", "conf.int_higher")])
            greater_than_threshold_fets[greater_than_threshold_fets == -Inf] <- 0
            less_than_threshold_fets <- unfactorize(data.frame(t(sapply(thresholds, function(threshold) { 
                m1 = sum(delta_pred_scores_cases < threshold); n1 = length(delta_pred_scores_cases); m0 = sum(delta_pred_scores_controls < threshold); n0 = length(delta_pred_scores_controls)
                fet_result <- fisher.test(matrix(c(m1, n1-m1, m0, n0-m0), nrow = 2, dimnames = list(hits = c("Y", "N"), status = c("case", "control"))), alternative = "two.sided")
                return(c(threshold, fet_result$estimate, fet_result$conf.int, fet_result$p.value, m1, m0, n1, n0)) 
            })))); colnames(less_than_threshold_fets) <- c("threshold", "estimate", "conf.int_lower", "conf.int_higher", "p.value", "m1", "m0", "n1", "n0")
            less_than_threshold_fets[less_than_threshold_fets == Inf] <- 0
            less_than_threshold_fets[,c("estimate", "conf.int_lower", "conf.int_higher")] <- log2(less_than_threshold_fets[,c("estimate", "conf.int_lower", "conf.int_higher")])
            less_than_threshold_fets[less_than_threshold_fets == -Inf] <- 0
            
            if(grepl("GOF",rbp)) { more_extreme_than_threshold_fets <- rbind(less_than_threshold_fets[less_than_threshold_fets$threshold <= 0,], greater_than_threshold_fets[greater_than_threshold_fets$threshold > 0,])
            } else { more_extreme_than_threshold_fets <- rbind(less_than_threshold_fets[less_than_threshold_fets$threshold < 0,], greater_than_threshold_fets[greater_than_threshold_fets$threshold >= 0,]) }
            write.csv(more_extreme_than_threshold_fets, file=output_path(paste0(gsub(" ","_",rbp),"_",variant_type,"_delta_prediction_score_enrichments.csv")), row.names=FALSE)
            
            if(grepl("GOF",rbp)) { more_extreme_than_threshold_fets <- more_extreme_than_threshold_fets[more_extreme_than_threshold_fets$threshold <= 0,]
            } else if(grepl("LOF",rbp)) { more_extreme_than_threshold_fets <- more_extreme_than_threshold_fets[more_extreme_than_threshold_fets$threshold >= 0,] }
            cols <- c("black","red")[as.numeric(more_extreme_than_threshold_fets$p.value < 0.05)+1]
            
            plot(more_extreme_than_threshold_fets$threshold, more_extreme_than_threshold_fets$estimate, main=paste0("",rbp," ",variant_type," case enrichment"), xlab="RBP binding delta threshold", ylab="log2(Odds Ratio)", col=cols, pch=19, ylim=c(min(more_extreme_than_threshold_fets$conf.int_lower), max(more_extreme_than_threshold_fets$conf.int_higher)), cex.axis=1.4, cex.lab=1.4, cex.main=1.4)
            abline(h=0, col="blue")
            abline(v=0, col="blue")
            if(!(grepl("LOF|GOF",rbp))) { text(x=c(-0.25, 0.25), y=rep(min(more_extreme_than_threshold_fets$conf.int_lower), 2), labels=c("GOF", "LOF"), font=2, adj=0.5, cex=1.2) }
            mtext("Binding delta more extreme than threshold", cex=1.2)
            segments(x0=more_extreme_than_threshold_fets$threshold, y0=more_extreme_than_threshold_fets$conf.int_lower, x1=more_extreme_than_threshold_fets$threshold, y1=more_extreme_than_threshold_fets$conf.int_higher, col=cols)
            segments(x0=more_extreme_than_threshold_fets$threshold-0.025, y0=more_extreme_than_threshold_fets$conf.int_lower, x1=more_extreme_than_threshold_fets$threshold+0.025, y1=more_extreme_than_threshold_fets$conf.int_lower, col=cols)
            segments(x0=more_extreme_than_threshold_fets$threshold-0.025, y0=more_extreme_than_threshold_fets$conf.int_higher, x1=more_extreme_than_threshold_fets$threshold+0.025, y1=more_extreme_than_threshold_fets$conf.int_higher, col=cols)
            # Draw case variants involved line on same plot
            par(new = T)
            plot(more_extreme_than_threshold_fets$threshold, more_extreme_than_threshold_fets$m1, type="l", lty=3, col="black", pch=16, axes=F, xlab=NA, ylab=NA, cex=1.4)
            axis(side = 4)
            mtext(side = 4, line = 3, "Case variants involved", cex=1.4)
            # Draw legend
            if(grepl("LOF",rbp)) { legend_location = "topright"
            } else { legend_location = "topleft" }
            legend(legend_location, legend=c("NS", "p<0.5", "count"), col=c("black","red","black"), pch=c(19,19,NA), lty=c(NA,NA,3), cex=1.1)
            
            #cols <- c("black","red")[as.numeric(less_than_threshold_fets$p.value < 0.05)+1]
            #plot(less_than_threshold_fets$threshold, less_than_threshold_fets$estimate, main=paste0("",rbp," GOF case enrichment"), xlab="RBP binding delta threshold", ylab="log2(Odds Ratio)", col=cols, pch=19, ylim=c(min(less_than_threshold_fets$conf.int_lower), max(less_than_threshold_fets$conf.int_higher)), cex.axis=1.4, cex.lab=1.4, cex.main=1.4)
            #abline(h=0, col="blue")
            #mtext("Binding delta < threshold", cex=1.2)
            #segments(x0=less_than_threshold_fets$threshold, y0=less_than_threshold_fets$conf.int_lower, x1=less_than_threshold_fets$threshold, y1=less_than_threshold_fets$conf.int_higher, col=cols)
            #legend("topright", legend=c("NS", "p<0.5"), col=c("black", "red"), pch=19, cex=1.2)
        }
        dev.off()
        pdf_to_png(filename)
        par(mfrow=c(1,1))
        par(mar=c(5.1,4.1,4.1,2.1))
    }
})

# Run variant threshold test to get real p.values
variant_types <- c("SNV", "indel")
variant_constraints <- c("autism_gene", "constrained_gene", "H3K36me3", "")
#thresholds <- seq(4, 150, by=2)
thresholds <- seq(0.5, 1, by=0.01)
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


for(rbp in rbps[45]) {
    print(rbp)
    model_name = paste0(tolower(rbp),"_model")
    model <- load_model_hdf5(output_path(paste0(model_name,".h5")))
    ref_pred_scores <- model %>% predict(regulatory_variant_refs_tensor[,,,])
    alt_pred_scores <- model %>% predict(regulatory_variant_alts_tensor[,,,])
    delta_pred_scores <- -log(alt_pred_scores/ref_pred_scores)
    
    for(variant_type in c("SNV", "Indel")) {
        print(variant_type)
        print(paste0("Mean: ",mean(delta_pred_scores[regulatory_variant_dat$Type == variant_type])))
        print("Quantiles:")
        print(quantile(delta_pred_scores[regulatory_variant_dat$Type == variant_type]))
    }
}


pdf(output_path("prediction_score_densities.pdf"))
plot(density(ref_pred_scores), col="blue", main="Prediction score densities")
lines(density(alt_pred_scores), col="red")
legend("topleft", legend=c("Refs", "Alts"), col=c("blue", "red"), pch=15)
dev.off()

delta_pred_scores <- -log(alt_pred_scores/ref_pred_scores)
plot(density(delta_pred_scores))

control_indices <- regulatory_variant_dat$Pheno == "control" #grepl("random", regulatory_variant_dat$sample)
delta_pred_scores_cases <- delta_pred_scores[!control_indices]
delta_pred_scores_controls <- delta_pred_scores[control_indices]
mean(delta_pred_scores_cases)
mean(delta_pred_scores_controls)
pdf(output_path("delta_prediction_score_densities.pdf"))
plot(density(delta_pred_scores_cases), col="red", main="Delta prediction score densities")
lines(density(delta_pred_scores_controls), col="blue")
legend("topleft", legend=c("case", "control"), col=c("red", "blue"), pch=15)
#legend("topleft", legend=c("HGMD", "gnomAD"), col=c("red", "blue"), pch=15)
dev.off()

#################################################################################################################################################
# Application of CNN to genomic data
#################################################################################################################################################

# Generate dummy data
library(gtools)
# Generate 1000 random sequences of length 7 each.
generate_random_sequences <- function(sequence_length, num_sequences) {
    random_sequences <- sapply(1:num_sequences, function(i) paste0(sample(c("A","C","G","T"), sequence_length, replace=TRUE), collapse=""))
    #H3K36me3_annotation <- rep(0,num_sequences); H3K36me3_annotation[sample(c(TRUE,FALSE), length(H3K36me3_annotation), prob=c(0.05,0.95), replace=TRUE)] <- 1
    #H3K4me1_annotation <- rep(0,num_sequences); H3K4me1_annotation[sample(c(TRUE,FALSE), length(H3K4me1_annotation), prob=c(0.4,0.6), replace=TRUE)] <- 1
    #phastCons_annotation <- rep(0,num_sequences); phastCons_annotation[sample(c(TRUE,FALSE), length(phastCons_annotation), prob=c(0.1,0.9), replace=TRUE)] <- 1
    #annotations <- list("H3K36me3"=H3K36me3_annotation, "H3K4me1"=H3K4me1_annotation, "evolutionary_conservation"=phastCons_annotation)
    #return_env <- new.env()
    #return_env[["sequences"]] <- random_sequences
    #return_env[["annotations"]] <- annotations
    #return(return_env)
    return(random_sequences)
}
#random_sequences <- generate_random_sequences(sequence_length=7, num_sequences=1000)
#data <- genomic_sequences_to_matrix(random_sequences, annotations)
#labels <- matrix(round(runif(num_sequences, min = 0, max = 1)), nrow=num_sequences, ncol=1) # { 1 = "RBP binding site", 0 = "not an RBP binding site" }
#num_features <- 4*sequence_length+length(annotations)

# Load real RBP eCLIP data (from Chaolin), just for one RBP: RBFOX2. We are going to predict additional binding sites for it across the genome! Or at least around specified variants of interest (like those in regulatory_variant_dat) for simpler computation time (represents predicted RBP feature).
rbp_dat <- read.csv("/home/local/ARCS/ak3792/Documents/Research/data/ENCODE_eCLIP_peak/K562.DROSHA.R2.tag.uniq.peak.sig.bed", sep="\t", header=FALSE)[,c(1:3,6)] #read.csv("/home/local/ARCS/ak3792/Documents/Research/data/ENCODE_eCLIP_peak/HepG2.RBFOX2.R2.tag.uniq.peak.sig.bed", sep="\t", header=FALSE)[,c(1:3,6)]

# Add RBP eCLIP features
RBP_ANNOTATIONS_PATH = data_path("ENCODE_eCLIP_peak")
rbps <- new.env()
rbp_files <- list.files(path=RBP_ANNOTATIONS_PATH)
num_rbps = length(rbp_files)
individual_rbp_dats <- sapply(1:num_rbps, function(i) { 
    rbp_eclip_file = rbp_files[i]
    rbp = paste0(strsplit(rbp_eclip_file, "\\.")[[1]][1:2], collapse=".")
    print(paste0("Storing RBP peaks for ",rbp, "... [",i," / ",num_rbps,"]"))
    rbp_dat <- read.csv(full_path(RBP_ANNOTATIONS_PATH, rbp_eclip_file), sep="\t", header=FALSE)[,c(1:3,6)]
    return(as.matrix(cbind(rbp_dat, rbp)))
})
rbp_dat <- data.frame(rbindlist(lapply(individual_rbp_dats, function(x) return(data.frame(x)))))
colnames(rbp_dat) <- c("chromosome", "start", "end", "strand", "RBP")
rbp_dat$start <- as.numeric(paste0(rbp_dat$start)); rbp_dat$end <- as.numeric(paste0(rbp_dat$end)); rbp_dat$RBP <- paste0(rbp_dat$RBP)

# Make sure start and end coordinates are properly ordered.
start_colname = "start"; end_colname = "end"
badly_ordered <- rbp_dat[,start_colname] > rbp_dat[,end_colname]
if(sum(badly_ordered) > 0) {
    badly_ordered_starts <- rbp_dat[badly_ordered, start_colname]
    rbp_dat[badly_ordered, start_colname] <- rbp_dat[badly_ordered, end_colname]
    rbp_dat[badly_ordered, end_colname] <- badly_ordered_starts
}
#sapply_out <- sapply(1:nrow(rbp_dat), function(i) { rbp_dat[i,2:3] <- sort(c(rbp_dat$end[i], rbp_dat$start[i])) } )

# Keep only unique rows
rbp_dat <- unique(rbp_dat)

# Add pillow/padding to the RBP peaks, if it is desired to loosen the RBP disruption label.
pillow = 0
rbp_dat$start <- rbp_dat$start - pillow; rbp_dat$end <- rbp_dat$end + pillow

rbp_sequences <- lapply(unique(rbp_dat$chromosome), function(chromosome) {
    if(grepl("Y",chromosome)) { return(matrix(nrow=0, ncol=2)) }
    cat(paste0("Reading chr", gsub("chr","",chromosome), "..."))
    refseq <- get_refseq(chromosome, version="hg19", allow_BSgenome=TRUE)[[1]]
    cat("Done.\nGrabbing sequences...")
    curr_chrom_indices <- which(rbp_dat$chromosome == chromosome)
    starts <- rbp_dat$start[curr_chrom_indices]; ends <- rbp_dat$end[curr_chrom_indices]; widths <- ends - starts + 1
    split_indices <- cumsum(widths)
    curr_chrom_rbp_sequences <- eval(parse(text=paste0("substring(paste0(refseq[c(",paste0(starts,":",ends, collapse=","),")]), c(0,split_indices[-length(split_indices)])+1, split_indices)")))
    cat("Done.\n")
    #num_curr_chrom_indices = length(curr_chrom_indices)
    #curr_chrom_rbp_sequences <- sapply(1:num_curr_chrom_indices, function(i) {
    #    cat(paste0("[",i," / ",num_curr_chrom_indices,"]"), "\n")
    #    curr_chrom_index = curr_chrom_indices[i]
    #    start = rbp_dat$start[curr_chrom_index]; end = rbp_dat$end[curr_chrom_index]
    #    return(paste0(refseq[start:end]))
    #})
    return(cbind(c(curr_chrom_rbp_sequences), c(rbp_dat$RBP[curr_chrom_indices]), chromosome, starts, ends))
})
rbp_sequences <- data.frame(rbindlist(lapply(rbp_sequences, function(x) return(data.frame(x)))))

rbp_sequences <- read.csv(data_path("rbp_sequences.csv"))

colnames(rbp_sequences) <- c("sequence", "RBP")
rbp_sequences$sequence <- paste0(rbp_sequences$sequence); rbp_sequences$RBP <- paste0(rbp_sequences$RBP)

rbp_sequences[1:10,]



#random_sequence_lengths <- table(nchar(rbp_sequences$sequence))
##random_sequences <- unlist(sapply(1:length(random_sequence_lengths), function(i) { print(random_sequence_lengths[i]); return(generate_random_sequences(sequence_length=as.numeric(paste0(names(random_sequence_lengths)))[i], num_sequences=random_sequence_lengths[i])) }))
##data <- genomic_sequences_to_tensor(c(rbp_sequences,random_sequences))
##labels <- c(rep(1,length(rbp_sequences)), rep(0,length(random_sequences)))

#rbps <- sort(unique(rbp_sequences$RBP))
#rbp_sequences_aggregated <- aggregate(rbp_sequences, by=list(rbp_sequences$sequence), FUN=function(x) { paste0(x[!duplicated(x)], collapse=",") } )[,c("sequence","RBP")]

#data_full <- genomic_sequences_to_tensor(rbp_sequences_aggregated$sequence)
#labels_full <- sapply(rbps, function(rbp) {
#    return(as.numeric(grepl(rbp, rbp_sequences_aggregated$RBP)))
#})

###### trinucleotideFrequency(refseq)  THIS IS WHAT TO USE FOR BACKGROUND MUTATION RATE!!!! VERY FAST AND EFFICIENT!
###### oligonucleotideFrequency(refseq, width=7)   THIS FUNCTION IS EVEN MORE GENERAL AND AMAZING FOR ANY k-MER SIZE!

#rbp_sequences_backup <- rbp_sequences

setup_features_env()
rbp_padding = 50000
all_rbp_features <- get_features_by_group("RBP")
rbps_to_focus_on <- c("K562.RBFOX2", "K562.QKI")
rbp_features <- all_rbp_features #[rowSums(sapply(rbps_to_focus_on, function(rbp) grepl(rbp, all_rbp_features))) == 1]
rbp_granges <- sapply(rbp_features, function(rbp_feature) { load_annotation(rbp_feature, padding=rbp_padding) })
rbp_granges <- lapply(rbp_granges, function(x) x[order(as.numeric(names(x)))]) # order by p.value/significance of peak 
sapply_out <- sapply(1:length(rbp_granges), function(i) { rbp_name <- gsub("_[^\\.]*$", "", names(rbp_granges)[i]); names(rbp_granges[[i]]) <<- rep(rbp_name, length(rbp_granges[[i]])) })
rbp_granges <- Reduce(c, rbp_granges)
rbp_granges_names <- names(rbp_granges)
rbps <- sort(unique(rbp_granges_names))

labels_full <- sapply(rbps, function(rbp) {
    print(rbp)
    rbp_indices <- rbp_granges_names == rbp
    rbp_granges_overlaps <- data.frame(findOverlaps(rbp_granges, rbp_granges[rbp_indices]))
    rbp_granges_starts <- start(rbp_granges); rbp_granges_ends <- end(rbp_granges)
    rbp_granges_overlaps_subjects <- rbp_granges_overlaps$subjectHits; rbp_granges_overlaps_queries <- rbp_granges_overlaps$queryHits
    rbp_granges_overlap_sizes <- pmin(abs(rbp_granges_starts[rbp_indices][rbp_granges_overlaps_subjects] - rbp_granges_ends[rbp_granges_overlaps_queries]), abs(rbp_granges_starts[rbp_granges_overlaps_queries] - rbp_granges_ends[rbp_indices][rbp_granges_overlaps_subjects])) + 1
    rbp_granges_overlaps <- cbind(rbp_granges_overlaps, rbp_granges_overlap_sizes); colnames(rbp_granges_overlaps)[ncol(rbp_granges_overlaps)] <- "overlap_size"
    
    overlap_sizes_vector <- rep(0, length(rbp_granges))
    overlap_sizes_per_query <- rbp_granges_overlaps[,c("queryHits", "overlap_size")]
    overlap_sizes_per_query <- aggregate(overlap_sizes_per_query$overlap_size, by=list(overlap_sizes_per_query$queryHits), FUN=max) # Get most overlapped/best score for each query RBP, for the current subject RBP
    indices <- overlap_sizes_per_query[,1]
    overlap_sizes_vector[indices]  <- overlap_sizes_per_query[,2]
    return(overlap_sizes_vector)
})
saveRDS(labels_full, file=data_path("all_rbp_overlap_size_labels_50kbp_padding.rds"))
write.csv(labels_full, file=data_path("all_rbp_overlap_size_labels_50kbp_padding.csv"), row.names=FALSE)

distances_full <- apply(labels_full, 2, function(x) pmax(0, (2 * rbp_padding) - x))
rownames(distances_full) <- rbp_granges_names
saveRDS(distances_full, file=data_path("all_rbp_distances.rds"))

rbp_eclip_counts <- table(rownames(distances_full))
rbp_closeness_threshold = 100 # bp
distances_means <- sapply(rbps, function(rbp) {
    print(rbp)
    rbp_indices <- rownames(distances_full) == rbp
    rbp_closeness_score <- colSums(distances_full[rbp_indices,] < rbp_closeness_threshold)/((rbp_eclip_counts+rbp_eclip_counts[rbp])/2) # percent of peaks for each RBP within the bp closeness threshold of this RBP's peaks #colMeans(distances_full[rbp_indices,])
    return(rbp_closeness_score)
})
write.csv(distances_means, file=data_path(paste0("all_rbp_closeness_",rbp_closeness_threshold,"bp.csv")))
tissues = c("K562", "HepG2")
for(tissue in tissues) {
    pdf(file=data_path(paste0("all_rbp_closeness_",rbp_closeness_threshold,"bp_heatmap_",tissue,".pdf")))
    #heatmap(distances_means, col=brewer.pal(9,"RdBu"), symm=TRUE, Rowv=NA, Colv=NA, cexRow=0.25, cexCol=0.25)
    heatmap(-(distances_means[grepl(tissue,rbps),grepl(tissue,rbps)])**0.25, col=brewer.pal(9,"RdBu"), symm=TRUE, cexRow=0.4, cexCol=0.4)
    dev.off()
}

#install.packages("rngtools", repos="http://R-Forge.R-project.org") # RNGtools compatible with R v3.5
library("NMF")
distances_nmf <- NULL
distances_nmf_measures <- data.frame()
ranks <- seq(5, 50, by=5)
nruns <- c(10, 30, 50, 80, 100)
for(nrun_curr in nruns) {
    print(nrun_curr)
    distances_nmf <- nmf(distances_means, rank=ranks, nrun=nrun_curr)
    distances_nmf_measures <- rbind(unfactorize(distances_nmf_measures), unfactorize(distances_nmf$measures))
}
distances_nmf$fit

#fit(distances_nmf)
#fitted(distances_nmf)
#cophcor(distances_nmf)
#dispersion(distances_nmf)


pdf(file=output_path(".pdf"))
par(mar = c(5,5,2,5))
plot(distances_nmf$measures$rank, distances_nmf$measures$dispersion, type="o", col="red", ylab=expression(-log[10](italic(p))), ylim=c(0,3))
par(new = T)
plot(distances_nmf$measures$rank, distances_nmf$measures$cophenetic, type="o", col="blue", pch=16, axes=F, xlab=NA, ylab=NA, cex=1.2)
axis(side = 4)
mtext(side = 4, line = 3, 'Number genes selected')
legend("topleft", legend=c(expression(-log[10](italic(p))), "N genes"), lty=c(1,0), pch=c(NA, 16), col=c("red3", "black"))
dev.off()

tissues = c("K562", "HepG2")
for(tissue in tissues) {
    pdf(file=data_path(paste0("all_rbp_closeness_",rbp_closeness_threshold,"bp_heatmap_",tissue,"_nmf.pdf")))
    #heatmap(distances_means, col=brewer.pal(9,"RdBu"), symm=TRUE, Rowv=NA, Colv=NA, cexRow=0.25, cexCol=0.25)
    #heatmap(apply(t(coef(distances_nmf)[,grepl(tissue,rbps)]), 2, function(x) -((x-mean(x))/sd(x))), col=brewer.pal(9,"RdBu"), Colv=NA, cexRow=0.4, cexCol=0.4)
    heatmap(apply(t(coef(distances_nmf)), 2, function(x) -((x-mean(x))/sd(x))), col=brewer.pal(9,"RdBu"), Colv=NA, cexRow=0.25, cexCol=0.25)
    dev.off()
}

#library("grid")
#library("ComplexHeatmap")
#library("circlize")
#heatmap_col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
#grid.newpage()
#pushViewport(viewport(layout = grid.layout(nr = 2, nc = 4, heights=c(0.15, 0.85))))
#
#pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
#grid.text(paste0("CNN layer activations for ",paste0(unique(input_label),collapse=",")," seq\n\"",input_sequence,"\"\n(",paste0(gsub("^.*\\.", "", colnames(labels)),"_score=",round(activations[[3]],3), collapse=", "),")"), x=0.1, y=0.55, gp=gpar(fontface = "bold", cex=1.5), just="left") #formatC(bonferroni_p.value_cutoff,format="e",digits=2)
#upViewport()
#
#lgd = Legend(at = c(-2, -1, 0, 1, 2), col_fun=heatmap_col_fun, title="Signal Strength", title_position="topleft")
#pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 4))
#grid.draw(lgd)
#upViewport()

#pushViewport(viewport(layout.pos.row = 2, layout.pos.col = c(1,2)))
#draw(Heatmap(activations[[1]], column_title = paste0("Conv1 (kernel_size=",conv1_kernel_size,")"), cluster_rows=FALSE, cluster_columns=FALSE, show_heatmap_legend = FALSE), newpage = FALSE) #row_order=
#upViewport()

#pushViewport(viewport(layout.pos.row = 2, layout.pos.col = c(3,4)))
#draw(Heatmap(activations[[2]], column_title = paste0("Conv2 (kernel_size=",conv2_kernel_size,")"), cluster_rows=FALSE, cluster_columns=FALSE, show_heatmap_legend = FALSE), newpage = FALSE) #row_order=
#upViewport()

#upViewport()
#dev.copy2pdf(file=output_path(paste0("cnn_layer_activations_",input_label,"_",input_sequence,".pdf")))

close_rbp_pairs <- paste0(unlist(sapply(1:nrow(distances_means), function(i) {
    distances_means_i <- distances_means[i,]
    close_rbps <- distances_means_i < 2000
    close_rbps[i] <- FALSE
    if(sum(close_rbps) > 0) { return(paste0(rbps[i]," - ",rbps[close_rbps],": ",round(distances_means_i[close_rbps])," bp")) } else { return(c()) }
})), collapse=", ")

a <- aggregate(distances_full, by=list(rownames(distances_full)), FUN=mean)


gr_sequences <- granges_to_DNAStringSet(rbp_granges)

#rbp_sequences <- data.frame(names(rbp_granges), seqnames(rbp_granges), start(rbp_granges), end(rbp_granges), strand(rbp_granges), paste0(gr_sequences))
#colnames(rbp_sequences) <- c("RBP", "chromosome", "start", "end", "strand", "sequence")
#write.csv(rbp_sequences, data_path("all_rbp_sequences.csv"), row.names=FALSE)
rbp_sequences <- read.csv(data_path("all_rbp_sequences.csv"))
rbp_sequences <- unfactorize(rbp_sequences)
rbps <- sort(unique(rbp_sequences$RBP))
##rbp_sequences_aggregated <- aggregate(rbp_sequences, by=list(rbp_sequences$sequence), FUN=function(x) { paste0(x[!duplicated(x)], collapse=",") } )[,c("sequence","RBP")]
#data_full <- genomic_sequences_to_tensor(gr_sequences, verbose=TRUE)
#saveRDS(data_full, file=data_path("all_rbp_tensor.rds"))


# rbp_granges_overlaps <- data.frame(findOverlaps(rbp_granges, rbp_granges))
# #rbp_granges_overlaps <- rbp_granges_overlaps[rbp_granges_overlaps$queryHits != rbp_granges_overlaps$subjectHits,]
# rbp_granges_starts <- start(rbp_granges); rbp_granges_ends <- end(rbp_granges)
# rbp_granges_overlaps_subjects <- rbp_granges_overlaps$subjectHits; rbp_granges_overlaps_queries <- rbp_granges_overlaps$queryHits
# rbp_granges_overlap_sizes <- pmin(abs(rbp_granges_starts[rbp_granges_overlaps_subjects] - rbp_granges_ends[rbp_granges_overlaps_queries]), abs(rbp_granges_starts[rbp_granges_overlaps_queries] - rbp_granges_ends[rbp_granges_overlaps_subjects])) + 1
# rbp_granges_overlaps <- cbind(rbp_granges_overlaps, rbp_granges_overlap_sizes, rbp_granges_names[rbp_granges_overlaps$subjectHits]); colnames(rbp_granges_overlaps)[(ncol(rbp_granges_overlaps)-1):ncol(rbp_granges_overlaps)] <- c("overlap_size", "subject_name")
# #a <- aggregate(rbp_granges_overlaps$queryHits, by=list(rbp_granges_overlaps$subject_name), FUN=table)
# rbps <- sort(unique(rbp_granges_names))
# labels_full <- sapply(rbps, function(rbp) {
#     print(rbp)
#     overlap_sizes_vector <- rep(0, length(rbp_granges))
#     overlap_sizes_per_query <- rbp_granges_overlaps[rbp_granges_overlaps$subject_name == rbp,c("queryHits", "overlap_size")]
#     overlap_sizes_per_query <- aggregate(overlap_sizes_per_query$overlap_size, by=list(overlap_sizes_per_query$queryHits), FUN=max) # Get most overlapped/best score for each query RBP, for the current subject RBP
#     indices <- overlap_sizes_per_query[,1]
#     overlap_sizes_vector[indices]  <- overlap_sizes_per_query[,2]
#     return(overlap_sizes_vector)
# }); colnames(labels_full) <- rbps
# saveRDS(labels_full, file=data_path("all_rbp_overlap_size_labels.rds"))
# write.csv(labels_full, file=data_path("all_rbp_overlap_size_labels.csv"), row.names=FALSE)
#saveRDS(labels_full, file=data_path("all_rbp_labels.rds"))


#h5createFile(data_path("all_rbp_tensor.h5"))
#h5write(data_full[,,,1], data_path("all_rbp_tensor.h5"),"df")
#data_full <- h5read(data_path("all_rbp_tensor.h5"), "S")

#labels_full <- sapply(rbps, function(rbp) {
#    return(as.numeric(grepl(rbp, rbp_sequences_aggregated$RBP)))
#}); colnames(labels_full) <- rbps

#data_full <- data
#labels_full <- labels
#sequence_length_full <- sequence_length

#################################################################################################################################################
# Simple CNN version! Includes motif/filter analysis in different layers.
#################################################################################################################################################
constrained_rbp_features <- all_rbp_features #[!(gsub("^.+\\.","",all_rbp_features) %in% ls(get_constrained_genes("pLI>0.5")))]
#rbp_neg <- read.csv(data_path("new_cnn_training_data/neg_rbp_seqs.csv"))
#rbp_pos <- read.csv(data_path("new_cnn_training_data/pos_rbp_seqs.csv"))
 
gene_expressions_matrix <- get_gene_expressions_matrix(c("K562","HepG2"))
#17:length(constrained_rbp_features)
training_metrics <- lapply(110:160, function(rbp_feature_i) { #1:length(constrained_rbp_features)
    rbps_to_focus_on <- constrained_rbp_features[rbp_feature_i] #QKI is 126, EFTUD2 is 87 c("K562.HNRNPU" is 104) #RBFOX2 is 127
    print(paste0(rbp_feature_i,". ",rbps_to_focus_on))
    multi_labels = FALSE
    multi_label_cutoff = 1
    process_labels_full = FALSE
    
    #if(multi_labels) {
    #    if(process_labels_full) {
    #        labels_full <- readRDS(file=data_path("all_rbp_overlap_size_labels.rds"))
    #        for(i in 1:nrow(labels_full)) { labels_full[i,] <- labels_full[i,]/max(labels_full[i,]) }
    #        filename = output_path(paste0("padded_eclip_rbp_overlap_percentage_density_",rbp_padding,"bp_padding.pdf"))
    #        pdf(file=filename)
    #        plot(density(labels_full[1:1000,], from=0, to=1), main="Padded eCLIP RBP overlap percentage density", xlab="Padded eCLIP RBP overlap percentage", xaxs="i", yaxs="i", lwd=2, col="red", cex.axis=1.4, cex.lab=1.4)
    #        mtext(paste0("+/- ",rbp_padding," bp-padded eCLIP sequences"), cex=1.2)
    #        dev.off()
    #        pdf_to_png(filename)
    #        for(i in 1:ncol(labels_full)) { print(i); very_high_overlapped <- labels_full[,i] >= multi_label_cutoff; very_low_overlapped <- labels_full[,i] <= (1 - multi_label_cutoff); labels_full[,i] <- rep(-1, nrow(labels_full)); labels_full[very_high_overlapped,i] <- 1; labels_full[very_low_overlapped,i] <- 0 }
    #        saveRDS(labels_full, file=data_path(paste0("all_rbp_overlap_size_labels_",multi_label_cutoff,".rds")))
    #    } else {
    #        labels_full <- readRDS(file=data_path(paste0("all_rbp_overlap_size_labels_",multi_label_cutoff,".rds")))
    #    }
    #    rows_to_pick <- sample(unique(c(which(rowSums(labels_full[,rbps_to_focus_on,drop=FALSE]) > 0))), 100000)
    #} else {
    #    labels_full <- readRDS(file=data_path("all_rbp_labels.rds"))
    #    rows_to_pick <- unique(c(which(rowSums(labels_full[,rbps_to_focus_on,drop=FALSE]) > 0)))
    #}
    #labels <- labels_full[rows_to_pick, rbps_to_focus_on, drop=FALSE]
    #if(!multi_labels) { labels[rowSums(labels) > 0,][labels[rowSums(labels) > 0,] == 0] <- -1 } # Sets 0 labels to unknown (-1) instead, if we have not yet run overlap analysis. 
    #rm(labels_full); gc()
    
    #data_full <- readRDS(file=data_path("all_rbp_tensor.rds"))
    #data <- data_full[rows_to_pick,,,,drop=FALSE]
    #rm(data_full); gc()
    
    rbp = rbps_to_focus_on
    rbp_cell_line_tpm = c("HepG2_logTPM","K562_logTPM")[as.numeric(grepl("K562\\.",rbp))+1]
    SEQUENCE_LENGTH = 151
    use_nearest_gene_expression = TRUE
    use_secondary_structure = TRUE
    
    gene_expression_path = output_path(paste0("gene_expressions/",tolower(rbps_to_focus_on),"_gene_expressions.csv"))
    if(file.exists(gene_expression_path)) {
        ge_already_stored = TRUE
        gene_expressions <- read.csv(gene_expression_path)
    } else { ge_already_stored = FALSE }
    secondary_structure_path = output_path(paste0("secondary_structures/",rbps_to_focus_on,"_ss.rds"))
    if(file.exists(secondary_structure_path)) {
        ss_already_stored = TRUE
        secondary_structure <- readRDS(secondary_structure_path)
    } else { ss_already_stored = FALSE }
    
    # Combine positive data set and labels with negative data set and labels
    data_positives_csv <- read.csv(data_path("RBP/pos_rbp_seqs_transcribed_regions.csv")) #new_cnn_training_data/pos_rbp_seqs.csv"))
    colnames(data_positives_csv)[2] <- "chromosome"
    indices_curr_rbp <- paste0(data_positives_csv$RBP) %in% rbps_to_focus_on
    data <- genomic_sequences_to_tensor(paste0(data_positives_csv$sequence[indices_curr_rbp]), sequence_length=SEQUENCE_LENGTH, verbose=TRUE)
    labels <- data.frame(rep(1, nrow(data))); colnames(labels) <- rbps_to_focus_on
    if(use_nearest_gene_expression && !ge_already_stored) { gene_expressions <- data.frame(gene_expressions_matrix[get_nearest_genes(data_positives_csv[indices_curr_rbp,]),rbp_cell_line_tpm]) }
    if(use_secondary_structure && !ss_already_stored) { secondary_structure <- annotate_rna_ss(data_positives_csv$sequence[indices_curr_rbp], as_tensor=TRUE, sequence_length=SEQUENCE_LENGTH) }
    rm(data_positives_csv)
    
    # Combine positive data set and labels with negative data set and labels
    data_negatives_csv <- read.csv(data_path("RBP/neg_rbp_seqs_transcribed_regions_11072019.csv"))[,-c(1)] #new_cnn_training_data/neg_rbp_seqs.csv")) #read.csv(data_path("negative_rbp_sequences/all_neg_seqs.csv"))
    colnames(data_negatives_csv)[2] <- "chromosome"
    indices_curr_rbp <- data_negatives_csv$RBP %in% rbps_to_focus_on
    data <- abind(data, genomic_sequences_to_tensor(paste0(data_negatives_csv$sequence[indices_curr_rbp]), sequence_length=ncol(data), verbose=TRUE), along=1)
    if(use_nearest_gene_expression && !ge_already_stored) {
        gene_expressions <- abind(gene_expressions, array(gene_expressions_matrix[get_nearest_genes(data_negatives_csv[indices_curr_rbp,]),rbp_cell_line_tpm], dim=c(nrow(data)- nrow(labels), 1)), along=1)
        colnames(gene_expressions) <- "logTPM"
        write.csv(gene_expressions, file=output_path(paste0(tolower(rbps_to_focus_on),"_gene_expressions.csv")), row.names=FALSE)
    }
    if(use_secondary_structure && !ss_already_stored) {
        secondary_structure <- abind(secondary_structure, annotate_rna_ss(data_negatives_csv$sequence[indices_curr_rbp], as_tensor=TRUE, sequence_length=ncol(data)), along=1)
        cat("Saving secondary structure to file...")
        saveRDS(secondary_structure, file=secondary_structure_path)
        cat("Done.\n")
    }
    rm(data_negatives_csv); gc()
    
    labels <- abind(labels, array(0, dim=c(nrow(data)- nrow(labels), ncol(labels))), along=1)
    
    data_without_ss <- data
    data <- abind(data, secondary_structure, along=3)
    
    dim(data)
    dim(labels)
    if(use_nearest_gene_expression) { print(dim(gene_expressions)) }
    if(use_secondary_structure) { print(dim(secondary_structure)) }
    num_bases = 4
    sequence_length = ncol(data) #max(nchar(rbp_sequences_aggregated$sequence)) #num_features/num_bases
    
    # Check and plot GC and sequence length distributions for positive and negative datasets.
    #i=1
    #pos_data_props <- data.frame(t(apply(data[labels[,i] == 1,,,], 1, FUN=function(x) { seq_length = sum(x); seq_gc = sum(x[,c("C","G")])/seq_length; return(c(paste0(tolower(rbps_to_focus_on[i]),"_pos"), seq_length, seq_gc)) })))
    #colnames(pos_data_props) <- c("label", "length", "GC"); pos_data_props$length <- as.numeric(paste0(pos_data_props$length)); pos_data_props$GC <- as.numeric(paste0(pos_data_props$GC))
    #neg_data_props <- data.frame(t(apply(data[labels[,i] == 0,,,], 1, FUN=function(x) { seq_length = sum(x); seq_gc = sum(x[,c("C","G")])/seq_length; return(c(paste0(tolower(rbps_to_focus_on[i]),"_neg"), seq_length, seq_gc)) })))
    #colnames(neg_data_props) <- c("label", "length", "GC"); neg_data_props$length <- as.numeric(paste0(neg_data_props$length)); neg_data_props$GC <- as.numeric(paste0(neg_data_props$GC))
    #
    #plot(density(pos_data_props$length), main="Sequence length distributions", xlab="length (bp)", lty=3, col="red", ylim=c(0, max(c(density(pos_data_props$length)$y, density(neg_data_props$length)$y))), cex.axis=1.4, cex.lab=1.4, cex.main=1.3)
    #lines(density(neg_data_props$length), lty=3, col="blue")
    #abline(v=mean(pos_data_props$length), col="red", lty=3); abline(v=mean(neg_data_props$length), col="blue", lty=3)
    #legend("topright", legend=c("pos", "neg"), col=c("red", "blue"), pch=15)
    #plot(density(pos_data_props$GC), main="Sequence GC distributions", xlab="GC", lty=3, col="red", ylim=c(0, max(c(density(pos_data_props$GC)$y, density(neg_data_props$GC)$y))), cex.axis=1.4, cex.lab=1.4, cex.main=1.3)
    #lines(density(neg_data_props$GC), lty=3, col="blue")
    #abline(v=mean(pos_data_props$GC), col="red", lty=3); abline(v=mean(neg_data_props$GC), col="blue", lty=3)
    #legend("topright", legend=c("pos", "neg"), col=c("red", "blue"), pch=15)
    
    sequence_length = ncol(data)
    # Number of filters to use
    num_filters = 100
    
    # Start with hidden 1D convolutional layer being fed 4 x sequence_length pixel images
    input <- layer_input(name="input", shape = c(sequence_length, dim(data)[[3]]))
    # First convolution layer
    conv1 <- layer_conv_1d(name="conv1", filters=num_filters, kernel_size=8, padding="valid", activation="relu", 
                           use_bias=FALSE, kernel_regularizer=regularizer_l2(l = 0.02))
    # Second convolution layer
    conv2 <- layer_conv_1d(name="conv2", filters=num_filters, kernel_size=4, activation="relu", use_bias=FALSE)
    # First max pooling layer
    maxpool1 <- layer_max_pooling_1d(name="maxpool1", pool_size = 8)
    # Second max pooling layer
    maxpool2 <- layer_max_pooling_1d(name="maxpool2", pool_size = 4)
    
    # Module responsible for processing tensor with CNN and learning matching sequence motifs.
    cnn_module <- input %>% 
        conv1 %>% 
        maxpool1 %>%
        layer_dropout(0.1, name="dropout_weak1") %>%
        conv2 %>% 
        maxpool2 %>% 
        layer_dropout(0.1, name="dropout_weak2") %>%
        # Flatten max filtered output into feature vector and feed into dense layer
        layer_flatten() %>%
        layer_dense(512, name="dense1", activation="relu") %>%
        layer_dense(512, name="dense2", activation="relu") %>%
        layer_dense(512, name="dense3", activation="relu")
        
    cnn_module_output <- cnn_module %>% layer_dropout(0.25, name="dropout_strong") %>%
        layer_dense(name="cnn_module_final_score", ncol(labels), activation="sigmoid")

    # Auxiliary input for nearest gene expression
    nearest_gene_expression_module <- layer_input(name="nearest_gene_expression", shape = c(1)) 
    
    # Final binding prediction output
    binding_predictions <- layer_concatenate(c(cnn_module_output, nearest_gene_expression_module)) %>% 
        layer_dense(name="final_score", ncol(labels), activation="sigmoid")
    
    # Secondary structure prediction output
    ss_predictions <- cnn_module %>% 
        layer_dense(name="ss_predictions", sequence_length, activation="sigmoid")

    #use_nearest_gene_expression = FALSE

    if(use_nearest_gene_expression) {
        if(use_secondary_structure) { 
            model <- keras_model(inputs=c(input, nearest_gene_expression_module), outputs=c(binding_predictions, ss_predictions))
        } else {
            model <- keras_model(inputs=c(input, nearest_gene_expression_module), outputs=c(binding_predictions))
        }
    } else {
        if(use_secondary_structure) { 
            model <- keras_model(inputs=c(input), outputs=c(cnn_module_output, ss_predictions))
        } else {
            model <- keras_model(inputs=c(input), outputs=c(cnn_module_output))
        }
    }
    opt <- optimizer_sgd(lr = 0.1, decay = 1e-2, momentum=0.5, nesterov=TRUE)
    masked_loss_function <- function(y_true, y_pred, mask=-1) { mask_vector <- k_cast(k_not_equal(y_true, mask), k_floatx());
                                                                return(k_binary_crossentropy(y_true * mask_vector, y_pred * mask_vector)) }
    #model %>% compile(loss = masked_loss_function, optimizer = opt, metrics = "accuracy")
    model %>% compile(loss = c("binary_crossentropy", "binary_crossentropy"), loss_weights=c(0.5, 0.5), optimizer = opt, metrics = "accuracy") #loss = "categorical_crossentropy"
    
    if(use_nearest_gene_expression) {
        if(use_secondary_structure) { 
            model_return <- analyze_model(model, data, labels, other_inputs=list(gene_expressions), other_outputs=list(secondary_structure[,,"paired",1]), model_name=paste0(tolower(rbps_to_focus_on),"_model3"), only_check_performance=TRUE, epochs=10)
        } else {
            model_return <- analyze_model(model, data, labels, other_inputs=list(gene_expressions), model_name=paste0(tolower(rbps_to_focus_on),"_model3"), only_check_performance=TRUE, epochs=10)
        }
    }
    else {
        if(use_secondary_structure) { 
            model_return <- analyze_model(model, data, labels, model_name=paste0(tolower(rbps_to_focus_on),"_model3"), only_check_performance=TRUE, epochs=20)
        } else {
            model_return <- analyze_model(model, data, labels, model_name=paste0(tolower(rbps_to_focus_on),"_model3"), only_check_performance=TRUE, epochs=20)
        }
    }
    model_return[[1]]
    
    return(model_return)
})

mean(unlist(training_metrics))
median(unlist(training_metrics))

score_breakpoint_size = 0.01
prediction_score_to_likelihood_mapping_table <- Reduce(function(df1, df2) merge(df1, df2, by="score", all.x=TRUE, all.y=TRUE), lapply(all_rbp_features, function(rbp) {
    rbp_score_densities <- lapply(c("positives", "negatives"), function(dataset) {
        rbp_score_density <- read.csv(output_path(paste0(tolower(rbp),"_model2_prediction_score_distribution_",dataset,".csv")))
        rbp_score_density$x <- round_to_nearest(rbp_score_density$x, score_breakpoint_size)
        rbp_score_density <-  rbp_score_density[rbp_score_density$x >= 0 & rbp_score_density$x <= 1,]
        rbp_score_density <- aggregate(rbp_score_density$y, by=list(rbp_score_density$x), FUN=mean)
        colnames(rbp_score_density) <- c("score", "density")
        rbp_score_density$density <- rbp_score_density$density / sum(rbp_score_density$density)
        return(rbp_score_density)
    })
    likelihood_ratio_table <- merge(rbp_score_densities[[1]], rbp_score_densities[[2]], by="score")
    likelihood_ratio_table <- cbind(likelihood_ratio_table, likelihood_ratio_table$density.x/likelihood_ratio_table$density.y)[,c(1,4)] #(likelihood_ratio_table$density.x + likelihood_ratio_table$density.y)
    colnames(likelihood_ratio_table) <- c("score", rbp)
    return(likelihood_ratio_table)
})); prediction_score_to_likelihood_mapping_table[is.na(prediction_score_to_likelihood_mapping_table)] <- 0
write.csv(prediction_score_to_likelihood_mapping_table, file=output_path("prediction_score_to_likelihood_mapping_table.csv"), row.names=FALSE)

pdf(file=output_path("likelihood_ratio_vs_prediction_score.pdf"))
rbps_to_plot <- c("RBFOX2", "EFTUD2", "QKI", "ILF3", "HNRNPU", "HNRNPA1")
cols <- rainbow(length(rbps_to_plot))
plot(prediction_score_to_likelihood_mapping_table$score, prediction_score_to_likelihood_mapping_table$HepG2.HNRNPU, type="l", ylim=c(0, 100), main="Likelihood Ratio vs Binding Score", xlab="RBP binding prediction score", ylab="Likelihood Ratio (positives/negatives)", cex.axis=1.4, cex.lab=1.4, cex.main=1.3)
for(i in 1:length(rbps_to_plot)) { lines(prediction_score_to_likelihood_mapping_table$score, prediction_score_to_likelihood_mapping_table[,paste0("HepG2.",rbps_to_plot[i])], col=cols[i]) }
legend("topleft", legend=rbps_to_plot, col=cols, lty=1, cex=1.2)
dev.off()

#################################################################################################################
# Train and save the compiled model, and write performance graphs and various other analysis plots to files. Includes motif/filter analysis in different layers.
#################################################################################################################
analyze_model <- function(model, data, labels, other_inputs=NULL, other_outputs=NULL, model_name="mymodel", skip_model_training=FALSE, only_check_performance=TRUE, epochs=5) {
    randomized_order_training_indices = sample(1:nrow(data), floor(0.8*nrow(data)))
    test_indices = (1:nrow(data))[!((1:nrow(data)) %in% randomized_order_training_indices)]
    if(!skip_model_training) {
        # Train the model, iterating on the data in batches of 32 samples
        if(is.null(other_inputs)) { 
            all_inputs <- data[randomized_order_training_indices,,,]
        } else {
            all_inputs <- c(list(data[randomized_order_training_indices,,,]), lapply(other_inputs, function(x) x[randomized_order_training_indices,]))
        }
        if(is.null(other_outputs)) { 
            all_outputs <- labels[randomized_order_training_indices,]
        } else {
            all_outputs <- c(list(labels[randomized_order_training_indices,]), lapply(other_outputs, function(x) x[randomized_order_training_indices,]))
        }
        history <- model %>% fit(all_inputs, all_outputs, epochs=epochs, batch_size=32, validation_split = 0.2)
        #plot.new()
        #plot(history)
        #dev.copy2pdf(file=output_path(paste0(model_name,"_training_history.pdf")))
        print(history$metrics)
        model %>% save_model_hdf5(output_path(paste0(model_name,".h5")))
        #model <- load_model_hdf5(output_path(paste0(model_name,".h5")))
    }
    # Draw model network
    plot_model(model, to_file = output_path(paste0(model_name,".png")), show_shapes = TRUE, show_layer_names = TRUE)
    plot_model(model, to_file = output_path(paste0("images/",model_name,".png")), show_shapes = TRUE, show_layer_names = TRUE)
    
    # Use model to predict labels for new data!
    #model %>% evaluate(data[test_indices,,,], labels[test_indices,])
    if(is.null(other_inputs)) { pred_scores <- model %>% predict(data[test_indices,,,]) 
    } else { pred_scores <- model %>% predict(c(list(data[test_indices,,,]), lapply(other_inputs, function(x) x[test_indices,]))) }
    pdf(output_path(paste0(model_name,"_prediction_score_distribution.pdf")))
    plot(density(pred_scores[,1]), main="RBP binding site prediction score distribution", col="white", xlim=c(0,1), ylim=c(0, max(sapply(1:ncol(pred_scores), function(i) max(c(density(pred_scores[,i][labels[test_indices,i] == 1])$y, density(pred_scores[,i][labels[test_indices,i] == 0])$y))))))
    cols <- rainbow(ncol(pred_scores))
    densities <- lapply(1:ncol(pred_scores), function(i) { 
        densities <- lapply(c(1,0), function(is_case) {
            ps <- pred_scores[,i][labels[test_indices,i] == is_case]
            print(ps)
            print(length(ps))
            return(density(ps))
        })
        lines(densities[[1]], col=cols[i], lwd=2, lty=1)
        lines(densities[[2]], col=cols[i], lwd=2, lty=3)
        return(densities)
    })
    write.csv(data.frame(densities[[1]][[1]][c("x", "y")]), output_path(paste0(model_name,"_prediction_score_distribution_positives.csv")), row.names=FALSE)
    write.csv(data.frame(densities[[1]][[2]][c("x", "y")]), output_path(paste0(model_name,"_prediction_score_distribution_negatives.csv")), row.names=FALSE)
    legend("topright", legend=c(gsub("^.*\\.", "", colnames(labels)), "Positives", "Negatives"), col=c(cols,"black","black"), pch=c(rep(15,ncol(labels)),NA,NA), lty=c(rep(NA,ncol(labels)),1,3))
    mtext(paste0("Validation set of ",length(test_indices)," length-",sequence_length," sequences"))
    dev.off()
    pdf_to_png(output_path(paste0(model_name,"_prediction_score_distribution.pdf")))
    
    pred_scores_old <- pred_scores
    labels_old <- labels
    roc_results <- c()
    for(i in 1:ncol(labels_old)) {
        rbp = colnames(labels_old)[i]
        print(paste0("Drawing ROC and PR curves for ",rbp))
        roc_results <- c(roc_results, get_roc_result(pred_scores_old[,i], labels_old[test_indices,i], plot_combined_only=TRUE, filename_prefix=output_path(paste0(model_name,"_",rbp)), mtext=paste0(rbp," binding sites"))[["auc_roc"]])
    }
    #plot(density(roc_results))
    
    if(only_check_performance) { return(roc_results) }
    
    if(FALSE) {
    # Directly look at weights in the trained model to determine learned motifs.
    conv1_weights <- get_weights(model$get_layer("conv1"))[[1]]
    get_weights(model$get_layer("conv1"))
    sapply_out <- sapply(1:num_filters, function(i) { 
        filter <- t(conv1_weights[,,i])
        #filter <- -filter
        rownames(filter) <- c("A","C","G","T")
        print("TGCATG")
        print(sequence_composite("TGCATG"))
        filter
        filter <- apply(filter, 2, function(x) { if(sum(x < 0)>0) { x <- x - min(x) }; return(x / sum(x)) })
        print(paste0("Calculating PWM for filter c1f",i,"..."))
        pwm_analysis(rbp_pwm=filter, rbp=paste0(model_name,"_c1f",i))
    })
    
    
    
    # Investigate trained model layers by feeding in new sequences.
    layers_to_investigate = c("conv1", "maxpool1", "final_score")
    model_layer <- model$get_layer(layers_to_investigate[1])
    model_layer_names <- sapply(model$layers, function(layer) paste0(layer$name))
    layers_to_investigate <- layers_to_investigate[layers_to_investigate %in% model_layer_names]
    activation_model <- keras_model(model$input, outputs=lapply(layers_to_investigate, function(layer_to_investigate) { return(get_output_at(model$get_layer(layer_to_investigate), 1)) })) #sapply(model$layers[2:3], function(layer) return(get_output_at(layer[[1]], 1)))) #layer[[1]]$output)
    rbp_to_analyze = "K562.RBFOX2"
    for(rbp_to_analyze in rbps_to_focus_on) {
        test_sequence_index = test_indices[labels[test_indices,rbp_to_analyze] == 1][1:200] #[127:128] #12 for RBFOX2, 128 for QKI
        single_input = length(test_sequence_index) == 1
        processed_labels <- gsub("^.*\\.", "", colnames(labels))
        input_label = apply(labels[test_sequence_index,,drop=FALSE], 1, function(x) paste0(sort(processed_labels[which(x == 1)]), collapse=","))
        #if(length(input_label) == 0) { input_label = "QKI" }
        input_label
        input_tensor <- data[test_sequence_index,,,]
        input_sequence = tensor_to_genomic_sequences(input_tensor)
        input_tensor_dims <- dim(input_tensor)
        #input_sequence
        conv1_kernel_size = as.numeric(gsub("[^0-9]", "", paste0(model$get_layer("conv1")$kernel_size)))
        conv2_kernel_size = as.numeric(gsub("[^0-9]", "", paste0(model$get_layer("conv2")$kernel_size)))
        if(single_input) {
            conv_names <- new.env()
            conv_names[["conv1"]] <- rollapply(1:nrow(input_tensor), width=conv1_kernel_size, function(positions) tensor_to_genomic_sequences(input_tensor[positions,]))
            #conv_names[["conv2"]] <- rollapply(conv_names[["conv1"]], width=conv2_kernel_size, function(seqs) paste0("")) #rollapply(conv_names[["conv1"]], width=conv2_kernel_size, function(seqs) paste0(seqs[seqs!=""], collapse="+"))
        }
        activation_model_outputs <- activation_model %>% predict(data[c(1,test_sequence_index),,,][-c(1),,,drop=FALSE])
        activations <- activation_model_outputs #lapply(1:length(activation_model_outputs), function(i) { if(grepl("conv1",layers_to_investigate[i])) { x <- activation_model_outputs[[i]][,,,drop=FALSE]; x_dims <- dim(x); if(single_input) { dimnames(x)[[length(x_dims)-1]] <- conv_names[[paste0("conv",i)]] }; dimnames(x)[[length(x_dims)]] <- paste0("c",i,"f",1:x_dims[length(x_dims)]) } else { x <- activation_model_outputs[[i]] }; return(x) })
        
        #heatmap(activations[[1]], Colv = NA, Rowv = NA, scale="column")
        
        try_to_align=TRUE
        num_best_sequence_candidates = 2
        motif_length = 8
        motif_length = min(c(motif_length, conv1_kernel_size))
        filter_contributions <- matrix(data=0, nrow=length(test_sequence_index), ncol=num_filters); colnames(filter_contributions) <- paste0("c1f",1:num_filters)
        activated_sequences <- lapply(1:length(test_sequence_index), function(i) {  #input_tensor_dims[length(input_tensor_dims)-1]
            print(paste0(i))
            #sequence <- input_sequence[i]
            conv1_windows <- unfactorize(data.frame(rollapply(1:input_tensor_dims[length(input_tensor_dims)-1], width=conv1_kernel_size, function(positions) return(c(tensor_to_genomic_sequences(input_tensor[i,positions,]), positions[1], positions[conv1_kernel_size])))))
            colnames(conv1_windows) <- c("motif", "start", "end")
            conv1_activation_sequences <- conv1_windows$motif #rollapply(1:input_tensor_dims[length(input_tensor_dims)-1], width=conv1_kernel_size, function(positions) tensor_to_genomic_sequences(input_tensor[i,positions,]))
            sequences_acceptable <- nchar(conv1_activation_sequences) >= motif_length
            conv1_activation_sequences <- conv1_activation_sequences[sequences_acceptable]
            conv1_windows <- conv1_windows[sequences_acceptable,]
            
            repped_grange <- rbp_granges[rows_to_pick][test_sequence_index[i]] # rbp_granges[which(sequence == gr_sequences)]
            
            conv1_activation_scores <- lapply(1:num_filters, function(filter_index) { activation_scores <- activations[[1]][i,sequences_acceptable,filter_index]; names(activation_scores) <- 1:length(conv1_activation_sequences); return(sort(activation_scores, decreasing=TRUE)) })
            names(conv1_activation_scores) <- paste0("c1f",1:length(conv1_activation_scores))
            best_sequences <- sapply(1:num_filters, function(filter_index) { 
                x <- conv1_activation_scores[[filter_index]]
                filter_contributions[i,filter_index] <<- max(x[1:num_best_sequence_candidates])
                best_sequence_indices <- as.numeric(names(x)[1:num_best_sequence_candidates])
                repped_granges <- rep(repped_grange, num_best_sequence_candidates)
                start(repped_granges) <- start(repped_granges) + conv1_windows$start[best_sequence_indices] - 1; end(repped_granges) <- start(repped_granges) + conv1_windows$end[best_sequence_indices] - conv1_windows$start[best_sequence_indices]
                return(repped_granges) 
            })
            names(best_sequences) <- paste0("c1f",1:num_filters)
            return(best_sequences)
        })
        activated_sequences <- lapply(1:num_filters, function(filter_index) { do.call(c, lapply(activated_sequences, function(x) x[[filter_index]])) })
        #activated_sequences <- data.frame(rbindlist(activated_sequences))
        filter_contributions <- colMeans(filter_contributions)
        top_k = 3
        top_k_filters <- order(filter_contributions, decreasing=TRUE)[1:top_k]
        sapply(top_k_filters, function(i) {
            print(paste0("Calculating PWM for filter c1f",i,"..."))
            pwm_analysis(paste0(model_name,"_",rbp_to_analyze,"_c1f",i), activated_sequences[[i]], bucket_size=5000, bucket=1)
        })
        print(paste0("Calculating combined PWM for ",rbp_to_analyze,"..."))
        pwm_analysis(paste0(model_name,"_",rbp_to_analyze,"_combined"), do.call(c, sapply(1:top_k, function(k) { activated_seqs <- activated_sequences[[top_k_filters[k]]]; return(sample(activated_seqs, floor(length(activated_seqs)*filter_contributions[top_k_filters[k]]))) })), bucket_size=5000, bucket=1) #do.call(c, activated_sequences[top_k_filters])
        print(sort(filter_contributions, decreasing=TRUE))
        print("TGCATG")
        print("CATGCA")
        
        pdf(output_path(paste0(model_name,"_",rbp_to_analyze,"_conv1_filter_contribution.pdf")))
        barplot_cols <- c(rainbow(top_k), rep("grey", num_filters-top_k))
        barplot(sort(filter_contributions, decreasing=TRUE), col=barplot_cols, main="CNN Conv1 Filter Contribution", ylab="Average sequence max activation", xlab=paste0(num_filters," Conv1 filters"), cex.main=1.4, cex.lab=1.4, cex.axis=1.4, names.arg=rep("",num_filters))
        legend("topright", title="Filter", legend=c(paste0("c1f",top_k_filters), "other"), col=barplot_cols[1:(top_k+1)], pch=15, cex=1.4)
        dev.off()
        pdf_to_png(output_path(paste0(model_name,"_",rbp_to_analyze,"_conv1_filter_contribution.pdf")))
    }
    #msa(granges_to_DNAStringSet(activated_sequences[[3]]), order="input", method="ClustalOmega", cluster=100)
    
    
    # This is the "african elephant" entry in the prediction vector
    african_elephant_output <- model$output[,1]
    # The is the output feature map of the `block5_conv3` layer,
    # the last convolutional layer in VGG16
    last_conv_layer <- model$get_layer("conv2")
    # This is the gradient of the "african elephant" class with regard to
    # the output feature map of `block5_conv3`
    grads <- k_gradients(african_elephant_output, last_conv_layer$output)[[1]]
    # This is a vector of shape (512,), where each entry
    # is the mean intensity of the gradient over a specific feature map channel
    pooled_grads <- k_mean(grads, axis=c(1,2))
    # This function allows us to access the values of the quantities we just defined:
    # `pooled_grads` and the output feature map of `block5_conv3`,
    # given a sample image
    iterate <- k_function(list(model$input), list(pooled_grads, last_conv_layer$output[1,,])) 
    # These are the values of these two quantities, as arrays,
    # given our sample image of two elephants
    library(abind)
    which(sapply(1:nrow(data), function(i) sum(data[i,,,])) > 144 & 1:nrow(data) %in% test_indices)
    which(test_indices == 29074)
    which(pred_scores > 0.9)
    sequence_index_to_display = 8
    determine_rbp_sequence_heatmap(8)
    }
    
    determine_rbp_sequence_heatmap <- function(sequence_index_to_display, motif_length = 8, motif_min_percent_distinct=0.5, motif_score_min_cutoff = 0.1) {
        # This is the "african elephant" entry in the prediction vector
        african_elephant_output <- model$output[,1]
        # The is the output feature map of the `block5_conv3` layer,
        # the last convolutional layer in VGG16
        last_conv_layer <- model$get_layer("conv1")
        # This is the gradient of the "african elephant" class with regard to
        # the output feature map of `block5_conv3`
        grads <- k_gradients(african_elephant_output, last_conv_layer$output)[[1]]
        if(is.null(grads)) { return(NULL) }
        # This is a vector of shape (512,), where each entry
        # is the mean intensity of the gradient over a specific feature map channel
        pooled_grads <- k_mean(grads, axis=c(1,2))
        # This function allows us to access the values of the quantities we just defined:
        # `pooled_grads` and the output feature map of `block5_conv3`,
        # given a sample image
        iterate <- k_function(list(model$input), list(pooled_grads, last_conv_layer$output[1,,]))
        
        #pred_scores[sequence_index_to_display]
        sequence_to_display = test_indices[sequence_index_to_display]
        iterate_result <- iterate(list(adrop(data[sequence_to_display,,,,drop=FALSE], 4)))
        pooled_grads_value <- iterate_result[[1]]; conv_layer_output_value <- iterate_result[[2]]
        # We multiply each channel in the feature map array
        # by "how important this channel is" with regard to the elephant class
        for (i in 1:100) {
            conv_layer_output_value[,i] <- conv_layer_output_value[,i] * pooled_grads_value[i]
        }
        # The channel-wise mean of the resulting feature map is our heatmap of class activation
        #heatmap <- smooth.spline(apply(conv_layer_output_value, 1, mean))$y
        heatmap <- smooth.spline(spline(apply(conv_layer_output_value, 1, mean), n=ncol(data), method="fmm")$y)$y
        # Normalize the heatmap between 0 and 1, for better visualization.
        heatmap <- pmax(heatmap, 0)
        heatmap <- heatmap / max(heatmap)
        
        sequence <- strsplit(tensor_to_genomic_sequences(data[sequence_to_display,,,]), "")[[1]]
        sequence_length <- length(sequence)
        sequence_matrix <- suppressWarnings(matrix(sequence, ncol=ceiling(sqrt(ncol(data))), byrow=TRUE))
        if((ncol(sequence_matrix) * nrow(sequence_matrix)) > sequence_length) { sequence_matrix[nrow(sequence_matrix), seq((sequence_length %% ncol(sequence_matrix))+1, ncol(sequence_matrix))] <- "" }
        rownames(sequence_matrix) <- paste0(seq(1, sequence_length, by=ncol(sequence_matrix)),"...",sapply(seq(ncol(sequence_matrix), sequence_length+ncol(sequence_matrix)-1, by=ncol(sequence_matrix)), function(x) min(c(x, sequence_length))))
        library(plotrix)
        filename=output_path(paste0(model_name,"_seq",sequence_to_display,"_rbp_binding_heatmap.pdf"))
        pdf(file=filename)
        plot(-1, -1, xlim=c(0,1), ylim=c(0,1), main=paste0(colnames(labels)[labels[sequence_to_display,]==1],paste0(" seq",sequence_to_display),", pred_score=",round(pred_scores[sequence_index_to_display],3)), xlab="", ylab="", xaxt="n", yaxt="n", axes=FALSE, cex.main=1.6, cex.lab=1.5) #"RBP binding strength heatmap"
        #mtext(paste0(colnames(labels)[labels[sequence_to_display,]==1]," sequence, pred_score=",round(pred_scores[sequence_index_to_display],3)))
        sequence_matrix_cols <- matrix(c(adjustcolor(rgb(colorRamp(c("blue","white","red"))(heatmap[1:sequence_length]), maxColorValue=255), alpha.f = 0.8), rep("white", (ncol(sequence_matrix) * nrow(sequence_matrix)) - sequence_length)), ncol=ceiling(sqrt(ncol(data))), byrow=TRUE)
        addtable2plot("topleft",table=sequence_matrix, bg=sequence_matrix_cols, display.rownames=TRUE, display.colnames=FALSE, cex=2.2)
        
        # Determine the optimal bound region of fixed motif_length from the heatmap 
        motif_start_index = which.max(rollapply(heatmap, width=motif_length, mean))
        motif_end_index = motif_start_index + motif_length - 1
        trimmed_binding_site <- paste0(sequence[motif_start_index:motif_end_index], collapse="")
        text(0.5, 0.01, paste0("Optimal binding site: ",trimmed_binding_site), cex=1.7, font=2)
        
        dev.off()
        pdf_to_png(filename)
        
        heatmap_smoothed <- rollapply(heatmap, width=motif_length, mean, align="left", partial=TRUE)
        motif_start_indices <- which(heatmap_smoothed >= motif_score_min_cutoff)
        motif_start_indices <- motif_start_indices[order(heatmap_smoothed[motif_start_indices], decreasing=TRUE)]
        get_best_nonoverlapping <- function(motif_start_indices) { if(length(motif_start_indices) > 0) { return(c(motif_start_indices[1], get_best_nonoverlapping(motif_start_indices[abs(motif_start_indices - motif_start_indices[1]) >= ceiling(motif_length*motif_min_percent_distinct)]))) } else { return(c()) } }
        motif_start_indices <- get_best_nonoverlapping(motif_start_indices)
        motif_end_indices <- motif_start_indices + motif_length - 1
        trimmed_binding_sites <- sapply(1:length(motif_start_indices), function(motif_i) paste0(sequence[motif_start_indices[motif_i]:motif_end_indices[motif_i]], collapse=""))
        
        motifs_dat <- unfactorize(data.frame(rep(paste0("seq",sequence_to_display),length(motif_start_indices)), rep(model_name,length(motif_start_indices)), rep(pred_scores[sequence_index_to_display],length(motif_start_indices)), motif_start_indices, motif_end_indices, trimmed_binding_sites, heatmap_smoothed[motif_start_indices], heatmap_smoothed[motif_start_indices]/heatmap_smoothed[motif_start_indices][1]))
        colnames(motifs_dat) <- c("seq", "model", "pred_score", "start", "end", "motif", "score", "score_norm")
        return(motifs_dat)
    }
    best_motifs <- rbindlist(lapply(which(pred_scores > 0.9)[1:10], function(j) { print(j); return(determine_rbp_sequence_heatmap(j)) }))
    best_motifs <- best_motifs[best_motifs$score > 0.33,]
    
    return_env <- new.env()
    return_env[["best_motifs"]] <- best_motifs
    return(best_motifs)
    
    #sort(unlist(activated_sequences), decreasing=TRUE)
    
    # if(!single_input) {
    #     #specific_gr_sequences <- paste0(gr_sequences[names(rbp_granges) == rbp_to_analyze])
    #     #specific_rbp_granges <- rbp_granges[names(rbp_granges) == rbp_to_analyze]
    #     filter_proportions <- colMeans(colMeans(activations[[1]])); filter_proportions <- filter_proportions / sum(filter_proportions)
    #     write.csv(data.frame(filter_proportions)[order(filter_proportions, decreasing=TRUE),,drop=FALSE], file=output_path(paste0(rbp_to_analyze,"_filter_proportions.csv")))
    #     
    #     calculate_filter_pwms = TRUE
    #     if(calculate_filter_pwms) {
    #         filter_sequences <- lapply(1:num_filters, function(i) {
    #             print(paste0("Grabbing activated sequences for filter c1f",i,"..."))
    #             filter_seqs <- unlist(lapply(1:length(input_sequence), function(j) {
    #                 sequence <- input_sequence[j]
    #                 conv1_reps <- floor(activations[[1]][j,,i]*5 + 0.5)
    #                 conv1_windows <- unfactorize(data.frame(rollapply(1:nrow(input_tensor[j,,]), width=conv1_kernel_size, function(positions) return(c(tensor_to_genomic_sequences(input_tensor[j,positions,]), positions[1], min(c(positions[conv1_kernel_size], nchar(sequence))))))))
    #                 colnames(conv1_windows) <- c("motif", "start", "end")
    #                 conv1_reps <- conv1_reps[conv1_windows$motif != ""]; conv1_windows <- conv1_windows[conv1_windows$motif != "",]
    #                 repped_sequences <- conv1_windows[unlist(sapply(1:nrow(conv1_windows), function(k) rep(k, conv1_reps[k]))),]
    #                 
    #                 repped_granges <- rbp_granges[as.numeric(sample(paste0(which(sequence == gr_sequences)), nrow(repped_sequences), replace=TRUE))]
    #                 start(repped_granges) <- start(repped_granges) + repped_sequences$start - 1; end(repped_granges) <- start(repped_granges) + repped_sequences$end - repped_sequences$start
    #                 return(repped_granges)
    #             }))
    #             return(Reduce(c, filter_seqs))
    #         })
    #         filter_sequences <- lapply(filter_sequences, function(x) x[width(x) >= 6])
    #         
    #         #genomic.acgt <- getBackgroundFrequencies("hg19")
    #         sapply(1:num_filters, function(i) {
    #             print(paste0("Calculating PWM for filter c1f",i,"..."))
    #             pwm_analysis(paste0(rbp_to_analyze,"_c1f",i), filter_sequences[[i]], bucket_size=5000, bucket=1)
    #         })
    #     }
    #     
    #     activations <- lapply(activations, function(x) colMeans(x))
    #     input_sequence <- "mean" 
    #     
    #     colMeans(activations[[1]])
    # }
    
    #library("grid")
    #library("ComplexHeatmap")
    #library("circlize")
    #heatmap_col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
    #grid.newpage()
    #pushViewport(viewport(layout = grid.layout(nr = 2, nc = 4, heights=c(0.15, 0.85))))
    #
    #pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
    #grid.text(paste0("CNN layer activations for ",paste0(unique(input_label),collapse=",")," seq\n\"",input_sequence,"\"\n(",paste0(gsub("^.*\\.", "", colnames(labels)),"_score=",round(activations[[3]],3), collapse=", "),")"), x=0.1, y=0.55, gp=gpar(fontface = "bold", cex=1.5), just="left") #formatC(bonferroni_p.value_cutoff,format="e",digits=2)
    #upViewport()
    #
    #lgd = Legend(at = c(-2, -1, 0, 1, 2), col_fun=heatmap_col_fun, title="Signal Strength", title_position="topleft")
    #pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 4))
    #grid.draw(lgd)
    #upViewport()
                      
    #pushViewport(viewport(layout.pos.row = 2, layout.pos.col = c(1,2)))
    #draw(Heatmap(activations[[1]], column_title = paste0("Conv1 (kernel_size=",conv1_kernel_size,")"), cluster_rows=FALSE, cluster_columns=FALSE, show_heatmap_legend = FALSE), newpage = FALSE) #row_order=
    #upViewport()
    
    #pushViewport(viewport(layout.pos.row = 2, layout.pos.col = c(3,4)))
    #draw(Heatmap(activations[[2]], column_title = paste0("Conv2 (kernel_size=",conv2_kernel_size,")"), cluster_rows=FALSE, cluster_columns=FALSE, show_heatmap_legend = FALSE), newpage = FALSE) #row_order=
    #upViewport()
    
    #upViewport()
    #dev.copy2pdf(file=output_path(paste0("cnn_layer_activations_",input_label,"_",input_sequence,".pdf")))
}

#################################################################################################################
# Draw ROC curve and/or Precision-Recall curve
#################################################################################################################
get_roc_result <- function(pred_scores_dat, labels, cols=NULL, curve_names=NULL, plot_roc=TRUE, plot_precision_recall=TRUE, plot_combined_only=TRUE, filename_prefix=NULL, mtext_text=NULL, mask=-1, legend.cex=1.3) {
    if(!is.null(filename_prefix) & !grepl("/",filename_prefix)) { filename_prefix <- output_path(filename_prefix) }
    if(is.null(ncol(pred_scores_dat))) { pred_scores_dat <- data.frame(pred_scores_dat) }
    num_curves = ncol(pred_scores_dat)
    if(is.null(curve_names)) { curve_names <- rep(" ", num_curves) }
    if(is.null(cols)) { cols <- rainbow(num_curves) }
    
    sensitivity <- new.env()
    one_minus_specificity <- new.env()
    precision <- new.env()
    auc_roc = new.env()
    auc_precision_recall = new.env()
    f1_score = new.env()
    
    known_labels <- labels != mask; labels <- labels[known_labels]
    
    print(paste0("Calculating ",num_curves," curves: "))
    for(i in 1:num_curves) {
        curve_name = curve_names[i]
        print(curve_name)
        pred_scores <- pred_scores_dat[known_labels,i]
        cutoffs <- seq(1.01,0,by=-0.01)
        roc_result <- sapply(cutoffs, function(cutoff) {
            calls <- pred_scores > cutoff
            sensitivity <- sum(calls[labels==1]==labels[labels==1])/sum(labels==1)
            one_minus_specificity <- 1 - (sum(calls[labels==0]==labels[labels==0])/sum(labels==0))
            precision <- sum(calls[labels==1]==labels[labels==1])/sum(calls==1)
            return(c(sensitivity, one_minus_specificity, precision))
        })
        sensitivity[[curve_name]] <- roc_result[1,]
        one_minus_specificity[[curve_name]] <- roc_result[2,]
        precision[[curve_name]] <- roc_result[3,]
        auc_roc[[curve_name]] = sum(diff(one_minus_specificity[[curve_name]]) * rollmean(sensitivity[[curve_name]], 2))
        auc_precision_recall[[curve_name]] = sum(diff(sensitivity[[curve_name]][!is.nan(precision[[curve_name]])]) * rollmean(precision[[curve_name]][!is.nan(precision[[curve_name]])], 2))
        f1_scores <- 2 * sensitivity[[curve_name]] * precision[[curve_name]] / (sensitivity[[curve_name]] + precision[[curve_name]])
        f1_score[[curve_name]] = max(f1_scores[!is.nan(f1_scores)])
    }
    
    draw_roc <- function() {
        plot(c(0,1), c(0,1), type="l", col="grey", xaxs="i", yaxs="i", xlab="1 - Specificity", ylab="Sensitivity", main="ROC Curve", cex.axis=1.3, cex.lab=1.3, cex.main=1.2)
        for(i in 1:num_curves) {
            lines(one_minus_specificity[[curve_names[i]]], sensitivity[[curve_names[i]]], col=cols[i], lwd=2)
        }
        if(num_curves == 1) {
            mtext(paste(c(mtext_text, paste0("F1 = ", round(f1_score[[curve_names[1]]],3)), paste0("AUC = ", round(auc_roc[[curve_names[1]]],3))), collapse=", "))
        } else {
            legend("bottomright", legend=sapply(1:num_curves, function(i) paste0(curve_names[i]," (","F1 = ", round(f1_score[[curve_names[i]]],3),", AUC = ", round(auc_roc[[curve_names[i]]],3),")")), col=cols, pch=rep(15,num_curves), cex=legend.cex)
            mtext(mtext_text)
        }
        #mtext(paste0(cv,"-fold CV, AUC_RBP = ", round(auc_all,3), ", AUC_no_RBP = ", round(auc_no_rbp,3)), cex=1.1)
        #legend("bottomright", legend=c("RBP", "No RBP", "CADD", "Eigen"), col=c("red", "green", "blue", "cyan"), pch=15)
    }
    draw_precision_recall <- function() {
        plot(c(0,1), c(0.5,0.5), type="l", col="grey", xaxs="i", yaxs="i", ylim=c(0,1), xlab="Recall", ylab="Precision", main="Precision-Recall Curve", cex.axis=1.3, cex.lab=1.3, cex.main=1.2)
        for(i in 1:num_curves) {
            lines(sensitivity[[curve_names[i]]][!is.nan(precision[[curve_names[i]]])], precision[[curve_names[i]]][!is.nan(precision[[curve_names[i]]])], col=cols[i], lwd=2)
        }
        if(num_curves == 1) {
            mtext(paste(c(mtext_text, paste0("F1 = ", round(f1_score[[curve_names[1]]],3)), paste0("PR AUC = ", round(auc_precision_recall[[curve_names[1]]],3))), collapse=", "))
        } else {
            legend("bottomright", legend=sapply(1:num_curves, function(i) paste0(curve_names[i]," (","F1 = ", round(f1_score[[curve_names[i]]],3),", PR AUC = ", round(auc_precision_recall[[curve_names[i]]],3),")")), col=cols, pch=rep(15,num_curves), cex=legend.cex)
            mtext(mtext_text)
        }
    }
    
    if(plot_roc && !plot_combined_only) { 
        # Draw ROC curve
        if(!is.null(filename_prefix)) { filename = paste0(filename_prefix,"_roc_curve.pdf"); pdf(file=filename) }
        draw_roc() 
        if(!is.null(filename_prefix)) { dev.off(); pdf_to_png(filename) }
    }
    if(plot_precision_recall) { 
        # Draw Precision-Recall curve
        if(!plot_combined_only) {
            if(!is.null(filename_prefix)) { filename = paste0(filename_prefix,"_precision_recall_curve.pdf"); pdf(file=filename) }
            draw_precision_recall() 
            if(!is.null(filename_prefix)) { dev.off(); pdf_to_png(filename) }
        }
        # Draw both performance curves side by side
        if(plot_roc && !is.null(filename_prefix)) {
            filename = paste0(filename_prefix,"_performance_metric_curves.pdf")
            pdf(file=filename, width=14)
            par(mfrow=c(1,2)) 
            draw_roc() 
            draw_precision_recall() 
            dev.off()
            pdf_to_png(filename) 
            par(mfrow=c(1,1)) 
        }
    }
    
    return_env <- new.env()
    return_env[["sensitivity"]] <- sensitivity; return_env[["one_minus_specificity"]] <- one_minus_specificity; return_env[["precision"]] <- precision; return_env[["cutoffs"]] <- cutoffs; return_env[["auc_roc"]] <- auc_roc; return_env[["auc_precision_recall"]] <- auc_precision_recall
    return(return_env)
}

#################################################################################################################
# Determine RBP binding sites along the input sequence
#################################################################################################################
find_binding_sites <- function(sequences, model, binding_sites=TRUE, conv_layer=NULL, model_name=NULL, draw=TRUE, motif_length = 8, motif_min_percent_distinct=0.5, motif_score_min_cutoff = 0.1) {
    #sequences_num_dimensions <- length(dim(sequences))
    sequences_class = class(sequences)
    if(sequences_class == "array") { #(sequences_num_dimensions == 4) {
        sequences_tensor <- sequences
        # Get prediction scores
        pred_scores <- model %>% predict(sequences_tensor[,,,])
        if(!binding_sites) { return(pred_scores) }
        sequences <- tensor_to_genomic_sequences(sequences)
    } else if (sequences_class == "character") { #(sequences_num_dimensions == 1) {
        sequences_tensor <- genomic_sequences_to_tensor(sequences, sequence_length=model_sequence_length)
        # Get prediction scores
        pred_scores <- model %>% predict(sequences_tensor[,,,])
        if(!binding_sites) { return(pred_scores) }
    }
    num_sequences <- length(sequences)
    
    # Get sequence length for this model.
    model_sequence_length = model$input$shape[1]
    # Get model name if one is not set
    if(is.null(model_name)) { model_name = model$name }
    # Find last convolutional layer
    if(is.null(conv_layer)) {
        model_layer_names <- sapply(model$layers, function(layer) paste0(layer$name))
        model_layer_names <- model_layer_names[grepl("conv",model_layer_names)]
        conv_layer = model_layer_names[length(model_layer_names)]
    }
    # This is the "african elephant" entry in the prediction vector
    african_elephant_output <- model$output[,1]
    # The is the output feature map of the `block5_conv3` layer,
    # the last convolutional layer in VGG16
    last_conv_layer <- model$get_layer(conv_layer)
    # This is the gradient of the "african elephant" class with regard to
    # the output feature map of `block5_conv3`
    grads <- k_gradients(african_elephant_output, last_conv_layer$output)[[1]]
    if(is.null(grads)) { return(NULL) }
    # This is a vector of shape (512,), where each entry
    # is the mean intensity of the gradient over a specific feature map channel
    pooled_grads <- k_mean(grads, axis=c(1,2))
    # This function allows us to access the values of the quantities we just defined:
    # `pooled_grads` and the output feature map of `block5_conv3`,
    # given a sample image
    iterate <- k_function(list(model$input), list(pooled_grads, last_conv_layer$output[1,,]))
        
    # Format sequences, and process each one at a time
    #if(length(dim(sequences)) > 3) { sequences_tensor <- adrop(sequences_tensor, 4) }
    sequences_tensor <- adrop(sequences_tensor, 4)
    binding_sites <- rbindlist(lapply(1:num_sequences, function(sequences_i) { 
        sequence = sequences[sequences_i]
        print(sequences_i)
        sequence_tensor = sequences_tensor[sequences_i,,,drop=FALSE]
        iterate_result <- iterate(list(sequence_tensor))
        pooled_grads_value <- iterate_result[[1]]; conv_layer_output_value <- iterate_result[[2]]
        # We multiply each channel in the feature map array by "how important this channel is" with regard to the elephant class
        for (i in 1:100) {
            conv_layer_output_value[,i] <- conv_layer_output_value[,i] * pooled_grads_value[i]
        }
        # The channel-wise mean of the resulting feature map is our heatmap of class activation
        heatmap <- smooth.spline(spline(apply(conv_layer_output_value, 1, mean), n=ncol(data), method="fmm")$y)$y
        # Normalize the heatmap between 0 and 1, for better visualization.
        heatmap <- pmax(heatmap, 0)
        heatmap <- heatmap / max(heatmap)
        
        sequence <- strsplit(sequence, "")[[1]]
        sequence_length <- length(sequence)
        
        heatmap_smoothed <- rollapply(heatmap, width=motif_length, mean, align="left", partial=TRUE)
        motif_start_indices <- which(heatmap_smoothed >= motif_score_min_cutoff)
        motif_start_indices <- motif_start_indices[order(heatmap_smoothed[motif_start_indices], decreasing=TRUE)]
        get_best_nonoverlapping <- function(motif_start_indices) { if(length(motif_start_indices) > 0) { return(c(motif_start_indices[1], get_best_nonoverlapping(motif_start_indices[abs(motif_start_indices - motif_start_indices[1]) >= ceiling(motif_length*motif_min_percent_distinct)]))) } else { return(c()) } }
        motif_start_indices <- get_best_nonoverlapping(motif_start_indices)
        motif_end_indices <- motif_start_indices + motif_length - 1
        trimmed_binding_sites <- sapply(1:length(motif_start_indices), function(motif_i) {
            sequence_to_paste <- sequence[motif_start_indices[motif_i]:motif_end_indices[motif_i]]
            sequence_to_paste <- sequence_to_paste[!is.na(sequence_to_paste)]
            return(paste0(sequence_to_paste, collapse=""))
        })
        motifs_dat <- unfactorize(data.frame(rep(paste0("seq",sequences_i),length(motif_start_indices)), rep(model_name,length(motif_start_indices)), rep(pred_scores[sequences_i],length(motif_start_indices)), motif_start_indices, motif_end_indices, trimmed_binding_sites, heatmap_smoothed[motif_start_indices], heatmap_smoothed[motif_start_indices]/heatmap_smoothed[motif_start_indices][1]))
        colnames(motifs_dat) <- c("seq", "model", "pred_score", "start", "end", "motif", "score", "score_norm")
        #motifs_dat <- motifs_dat[motifs_dat$motif != "",]
        motifs_dat <- motifs_dat[nchar(motifs_dat$motif) == motif_length,]
        
        if(draw) {
            sequence_matrix <- suppressWarnings(matrix(sequence, ncol=ceiling(sqrt(ncol(data))), byrow=TRUE))
            if((ncol(sequence_matrix) * nrow(sequence_matrix)) > sequence_length) { sequence_matrix[nrow(sequence_matrix), seq((sequence_length %% ncol(sequence_matrix))+1, ncol(sequence_matrix))] <- "" }
            rownames(sequence_matrix) <- paste0(seq(1, sequence_length, by=ncol(sequence_matrix)),"...",sapply(seq(ncol(sequence_matrix), sequence_length+ncol(sequence_matrix)-1, by=ncol(sequence_matrix)), function(x) min(c(x, sequence_length))))
            library(plotrix)
            filename=output_path(paste0(model_name,"_seq",sequences_i,"_rbp_binding_heatmap.pdf"))
            pdf(file=filename)
            plot(-1, -1, xlim=c(0,1), ylim=c(0,1), main=paste0(colnames(labels)[labels[sequences_i,]==1],paste0(" seq",sequences_i),", pred_score=",round(pred_scores[sequences_i],3)), xlab="", ylab="", xaxt="n", yaxt="n", axes=FALSE, cex.main=1.6, cex.lab=1.5) #"RBP binding strength heatmap"
            #mtext(paste0(colnames(labels)[labels[sequence_to_display,]==1]," sequence, pred_score=",round(pred_scores[sequence_index_to_display],3)))
            sequence_matrix_cols <- matrix(c(adjustcolor(rgb(colorRamp(c("blue","white","red"))(heatmap[1:sequence_length]), maxColorValue=255), alpha.f = 0.8), rep("white", (ncol(sequence_matrix) * nrow(sequence_matrix)) - sequence_length)), ncol=ceiling(sqrt(ncol(data))), byrow=TRUE)
            addtable2plot("topleft",table=sequence_matrix, bg=sequence_matrix_cols, display.rownames=TRUE, display.colnames=FALSE, cex=2.2)
            # Write the optimal bound region of fixed motif_length from the heatmap 
            text(0.5, 0.01, paste0("Optimal binding site: ",motifs_dat$motif[which.max(motifs_dat$score)]), cex=1.7, font=2)
            dev.off()
            pdf_to_png(filename)
        }
        
        return(motifs_dat)
    }))
    return(binding_sites)
}
best_motifs <- find_binding_sites(data[1:5,,,,drop=FALSE], model, model_name=model_name, binding_sites=TRUE)
best_motifs <- best_motifs[best_motifs$score > 0.33,]

#################################################################################################################################################
# More complex model with multiple inputs and outputs
#################################################################################################################################################

rbps_to_focus_on <- c("K562.RBFOX2", "K562.EFTUD2", "K562.HNRNPU", "K562.ILF3", "K562.QKI")
multi_labels = FALSE
multi_label_cutoff = 1
process_labels_full = FALSE

if(multi_labels) {
    if(process_labels_full) {
        labels_full <- readRDS(file=data_path("all_rbp_overlap_size_labels.rds"))
        for(i in 1:nrow(labels_full)) { labels_full[i,] <- labels_full[i,]/max(labels_full[i,]) }
        filename = output_path(paste0("padded_eclip_rbp_overlap_percentage_density_",rbp_padding,"bp_padding.pdf"))
        pdf(file=filename)
        plot(density(labels_full[1:1000,], from=0, to=1), main="Padded eCLIP RBP overlap percentage density", xlab="Padded eCLIP RBP overlap percentage", xaxs="i", yaxs="i", lwd=2, col="red", cex.axis=1.4, cex.lab=1.4)
        mtext(paste0("+/- ",rbp_padding," bp-padded eCLIP sequences"), cex=1.2)
        dev.off()
        pdf_to_png(filename)
        for(i in 1:ncol(labels_full)) { print(i); very_high_overlapped <- labels_full[,i] >= multi_label_cutoff; very_low_overlapped <- labels_full[,i] <= (1 - multi_label_cutoff); labels_full[,i] <- rep(-1, nrow(labels_full)); labels_full[very_high_overlapped,i] <- 1; labels_full[very_low_overlapped,i] <- 0 }
        saveRDS(labels_full, file=data_path(paste0("all_rbp_overlap_size_labels_",multi_label_cutoff,".rds")))
    } else {
        labels_full <- readRDS(file=data_path(paste0("all_rbp_overlap_size_labels_",multi_label_cutoff,".rds")))
    }
    rows_to_pick <- sample(unique(c(which(rowSums(labels_full[,rbps_to_focus_on,drop=FALSE]) > 0))), 100000)
} else {
    labels_full <- readRDS(file=data_path("all_rbp_labels.rds"))
    rows_to_pick <- unique(c(which(rowSums(labels_full[,rbps_to_focus_on,drop=FALSE]) > 0)))
}
labels <- labels_full[rows_to_pick, rbps_to_focus_on, drop=FALSE]
if(!multi_labels) { labels[rowSums(labels) > 0,][labels[rowSums(labels) > 0,] == 0] <- -1 } # Sets 0 labels to unknown (-1) instead, if we have not yet run overlap analysis. 
rm(labels_full); gc()

data_full <- readRDS(file=data_path("all_rbp_tensor.rds"))
data <- data_full[rows_to_pick,,,,drop=FALSE]
rm(data_full); gc()

# Combine positive data set and labels with negative data set and labels
#data_negatives <- genomic_sequences_to_tensor(paste0(read.csv(data_path("negative_rbp_sequences/RBFOX2_neg.csv"))$seq), sequence_length=ncol(data), verbose=TRUE)
data_negatives_csv <- read.csv(data_path("negative_rbp_sequences/all_neg_seqs.csv"))
data_negatives <- genomic_sequences_to_tensor(paste0(data_negatives_csv$sequence[data_negatives_csv$RBP %in% rbps_to_focus_on]), sequence_length=ncol(data), verbose=TRUE)
data <- abind(data, data_negatives, along=1)
rm(data_negatives); gc()
labels <- abind(labels, array(0, dim=c(nrow(data)- nrow(labels), ncol(labels))), along=1)

dim(data)
dim(labels)
num_bases = 4
sequence_length = ncol(data) #max(nchar(rbp_sequences_aggregated$sequence)) #num_features/num_bases
num_rbps = ncol(labels)


# Number of filters to use
num_filters = 100

# Start with hidden 1D convolutional layer being fed 4 x sequence_length pixel images
input <- layer_input(name="sequence_tensor", shape=c(sequence_length, 4))
# First convolution layer
conv1 <- layer_conv_1d(name="conv1", filters=num_filters, kernel_size=8, padding="valid", activation="relu", use_bias=FALSE, kernel_regularizer=regularizer_l2(l = 0.05))
# Second convolution layer
conv2 <- layer_conv_1d(name="conv2", filters=num_filters, kernel_size=2, activation="relu", use_bias=FALSE)
# First max pooling layer
maxpool1 <- layer_max_pooling_1d(name="maxpool1", pool_size = 8)
# Second max pooling layer
maxpool2 <- layer_max_pooling_1d(name="maxpool2", pool_size = 1)

# Auxiliary input, corresponding simply to the RBP indices, for the embedding table.
aux_input <- layer_input(name="input_to_embedding_layer", shape = c(1))
rbp_indices_tensor <- t(sapply(1:nrow(labels), function(i) c(0:(num_rbps-1)))) #diag(nrow=num_rbps) #rep(1, num_rbps)
colnames(rbp_indices_tensor) <- colnames(labels)
# RBP embedding table
rbp_embedding_dimensionality = 20
rbp_embeddings <- aux_input %>% layer_embedding(name="rbp_embeddings", input_dim=num_rbps, output_dim=rbp_embedding_dimensionality, input_length=num_rbps) %>% layer_flatten()

# Get motif features from sequence by convolution, and flatten max filtered output into feature vector
conv_output <- input %>% 
    conv1 %>% 
    maxpool1 %>%
    layer_dropout(0.1, name="dropout_weak1") %>%
    conv2 %>% 
    maxpool2 %>% 
    layer_dropout(0.1, name="dropout_weak2") %>%
    layer_flatten()
    
# Combine information from  
predictions <- conv_output %>% #layer_concatenate(c(conv_output, rbp_embeddings)) %>%  
    layer_dense(512, activation="relu", name="dense1") %>%
    layer_dense(64, activation="relu", name="dense2") %>%
    layer_dense(16, activation="relu", name="dense3") %>%
    layer_dropout(0.5, name="dropout_strong") %>%
    layer_dense(name="final_score", num_rbps, activation="sigmoid")

#model <- keras_model(inputs=c(input, aux_input), outputs=predictions)
model <- keras_model(inputs=c(input), outputs=predictions)

opt <- optimizer_sgd(lr = 0.1, decay = 1e-2, momentum=0.5, nesterov=TRUE)
masked_loss_function <- function(y_true, y_pred, mask=-1) { mask_vector <- k_cast(k_not_equal(y_true, mask), k_floatx()); return(k_binary_crossentropy(y_true * mask_vector, y_pred * mask_vector)) }
model %>% compile(loss = masked_loss_function, optimizer = "rmsprop", metrics = "accuracy")

analyze_model(model, data, labels, model_name="mymodel", skip_model_training=FALSE, only_check_performance=FALSE, epochs=5) #other_inputs=list(rbp_indices_tensor)


model_name = "mymodel"
plot_model(model, to_file = output_path(paste0(model_name,".png")), show_shapes = TRUE, show_layer_names = TRUE)
randomized_order_training_indices = sample(1:nrow(data), floor(0.8*nrow(data)))
test_indices = (1:nrow(data))[!((1:nrow(data)) %in% randomized_order_training_indices)]
history <- model %>% fit(x=list(data[randomized_order_training_indices,,,], rbp_indices_tensor[randomized_order_training_indices,]), labels[randomized_order_training_indices,], epochs=5, batch_size=32, validation_split = 0.2)
plot(history)
# Plot RBP embeddings
rbp_embedding_table <- get_weights(model$get_layer("rbp_embeddings"))[[1]]
rbp_embedding_table_pca <- prcomp(rbp_embedding_table)$x
cols <- rainbow(num_rbps)
filename=output_path("rbp_embedding_table_pca.pdf")
pdf(file=filename)
plot(rbp_embedding_table_pca[,1], rbp_embedding_table_pca[,2], main="RBP Embedding Table PCA", xlab="PC1", ylab="PC2", xlim=range(rbp_embedding_table_pca[,1])*1.2, ylim=range(rbp_embedding_table_pca[,2])*1.2, col=cols, cex.lab=1.4, cex.axis=1.4, cex.main=1.4)
text(rbp_embedding_table_pca[,1], rbp_embedding_table_pca[,2], labels=colnames(labels), pos=3, offset=0.75, cex=1, col=cols)
mtext(paste0(rbp_embedding_dimensionality,"-dimensional embedding for ",num_rbps," RBPs"), cex=1.2)
dev.off()
pdf_to_png(filename)

#################################################################################################################
# Set up features_env environment that will store named features (filepaths, data frames, or GRanges objects)
#################################################################################################################
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
setup_features_env(features_to_make_fullPeak=c("E118.H3K36me3.broadPeak", "E123.H3K36me3.broadPeak", "E080.H3K36me3.broadPeak")) # E118 is HepG2, E123 is K562, E080 is adrenal_gland
#write.table(file=output_path("E123.H3K36me3.fullPeak"), genomic_regions_to_dat(load_annotation("E123.H3K36me3.fullPeak")))
saveRDS(load_annotation("E123.H3K36me3.fullPeak"), file = output_path("K562_E123.H3K36me3.fullPeak.rds"))
saveRDS(load_annotation("E118.H3K36me3.fullPeak"), file = output_path("HepG2_E118.H3K36me3.fullPeak.rds"))
saveRDS(load_annotation("E080.H3K36me3.fullPeak"), file = output_path("adrenal_gland_E080.H3K36me3.fullPeak.rds"))

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
# Loads the specified feature from features_env.
#################################################################################################################
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

#################################################################################################################
# Annotate the given variants (passed in as either data frame, GRanges, or variant_list format) with the specified features (saved in features_env as either GRanges, data frames, or filepaths)
#################################################################################################################
annotate <- function(variants, features, features_env=NULL, variants_granges=NULL, save_globally=FALSE) {
    if(is.null(features_env)) { features_env <- get_global("features_env") }; if (is.null(features_env)) { features_env <- setup_features_env() }
    if(class(variants) == "GRanges") { variants_granges <- variants; variants <- genomic_regions_to_dat(variants) } # handles case where variants are passed in directly as GRanges
    variants <- standardize_colnames(variants)
    if (is.null(variants_granges)) { variants_granges <- to_genomic_regions(variants, chr_colname="Chrom", start_colname="Position", end_colname="Position") }
    num_features = length(features)
    available_features <- ls(features_env)
    available_features_notes <- sapply(available_features, function(x) { return(unlist(strsplit(features_env[[x]][["notes"]], "; *"))) })
    # Function to load annotation GRanges for a given feature and padding size.
    # load_annotation <- function(feature, padding=0, save_globally=save_globally) {
    #     annotations <- features_env[[gsub("_[0-9]+bp$", "", feature)]][["peaks"]]
    #     if(class(annotations)!="GRanges") { # If annotations is not already a GRanges class, it must be a data frame or a filepath from which we will construct a new GRanges object.
    #         if(class(annotations)!="data.frame") { # If annotations is not already a data frame, it must be a filepath from which we will construct a new data frame and then GRanges object.
    #             if (grepl("csv", annotations)) { dat_sep = "," } else { dat_sep = "\t" }
    #             annotations <- read.csv(annotations, sep=dat_sep, header=FALSE); colnames(annotations)[1:3] <- c("chromosome", "start", "end")
    #             strand_col <- which(apply(annotations, 2, function(annotations_col) sum(paste0(annotations_col) %in% c("+", "-")) == nrow(annotations)))[1]
    #             if(!is.na(strand_col)) { colnames(annotations)[strand_col] <- "strand" }
    #         }
    #         annotations_labels <- genomic_coordinates_to_strings(annotations$chromosome, annotations$start, annotations$end)
    #         p.value_col <- which(tolower(colnames(annotations)) %in% c("p.value", "p"))[1]
    #         if(is.na(p.value_col)) { 
    #             p.value_col <- which(apply(annotations, 2, function(annotations_col) sum(grepl("[pP]=[0-9]", annotations_col)) == nrow(annotations)))[1] 
    #             if(!is.na(p.value_col)) { annotations_labels <- gsub(".*\\[[pP]=([0-9\\.eE\\-]+)\\].*","\\1", annotations[,p.value_col]) }
    #         } else { annotations_labels <- paste0(annotations[,p.value_col]) }
    #         annotations <- to_genomic_regions(annotations, labels=annotations_labels)
    #     }
    #     #print(paste0("Feature: ",feature,", Padding in load_annotation: ",padding))
    #     if(padding != 0) { start(annotations) <- start(annotations) - padding; end(annotations) <- end(annotations) + padding }
    #     #print(annotations)
    #     if(save_globally) { 
    #         if(is.null(features_env[[paste0(feature)]])) { features_env[[paste0(feature)]] <- new.env(); features_env[[paste0(feature)]][["notes"]] <- ""; }
    #         features_env[[paste0(feature)]][["peaks"]] <- annotations # Save annotations in features_env as a GRanges object for future efficiency
    #         set_global(features_env)
    #     }
    #     return(annotations)
    # }
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
# Return multiple sequence alignment of the given sequences, with gaps filled in with the correct neighboring sequences.
#################################################################################################################
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

#################################################################################################################
# Analyze PWM, either given directly or derived from sequences, and draw sequence logo
#################################################################################################################
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












