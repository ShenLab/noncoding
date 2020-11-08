# This script is a suite of general R tools and functions written by Alexander Kitaygorodsky.

library(gtools)
library(seqinr)
library("GenomicRanges")
library("BSgenome.Hsapiens.UCSC.hg19")
library("rtracklayer")
library(igraph)
library("zoo")
library("car")
#library("rapportools")
library("deconstructSigs")
library("RColorBrewer")
#library("GenomicFeatures")
library("purrr")
library("data.table")
library("KEGGREST")
library("readxl")
library('annovarR')
library("AnnotationDbi")

#devtools::install_github("JhuangLab/annovarR")
#library(BioInstaller)
#install.bioinfo('annovar', '/home/local/ARCS/ak3792/Documents/Research/data/annovar')
# Use BioInstaller to install vcfanno easily in R
#install.bioinfo('vcfanno', '/home/local/ARCS/ak3792/Documents/Research/data/vcfanno')
                
#download.database(show.all.names = TRUE)
#download.database('db_annovar_exac03', database.dir = "/data/annovar/humandb", buildver = "hg19")
#download.database('db_annovar_refgene', database.dir = "/data/annovar/humandb", buildver = "hg19")
#download.database('db_ucsc_cytoband', database.dir = "/data/annovar/humandb", buildver = "hg19")
#download.database('db_annovar_cadd', database.dir = "/data/annovar/humandb", buildver = "hg19")
##All annovarR supported big annotation database required SQLite format
#download.database('db_annovar_avsnp147_sqlite', database.dir = "/data/annovar", buildver = "hg19")
#download.database('db_annovar_avsnp', database.dir = sprintf('%s/databases/', tempdir()), show.all.versions = TRUE)
#download.database('db_annovar_avsnp', database.dir = sprintf('%s/databases/', tempdir()), version="avsnp150")
#
#library(data.table)
#library("RSQLite")
#
#database.dir <- "/data/annovar/humandb"
#database.cfg <- system.file('extdata', 'config/databases.toml', package = "annovarR")
##is.toml.file(database.cfg)
##get.config.type(database.cfg)
##library("RcppTOML")
##parseTOML(database.cfg)
#get.annotation.names()
#download.name <- "db_annovar_refgene" #get.download.name("exac03nontcga")
#download.name
#download.database(download.name = download.name, buildver = "hg19", database.dir = database.dir, db.type="txt")
chr <- c("chr22", "chr2", "chr1")
start <- c("46615880", "100030", "200020")
end <- c("46615880", "100035", "200030")
ref <- c("T", "A", "A")
alt <- c("C", "T", "C")
dat <- data.table(Chrom = chr, start = start, end = end, ref = ref, alt = alt)
##x <- annotation.merge(dat = dat, anno.name = c("ucsc_refgene", "exac03nontcga"), database.dir = database.dir, db.type = 'txt')
#x <- annotation.merge(dat = dat, anno.name = c("ucsc_refgene", "exac03"), database.dir = database.dir, db.type = 'txt')
#
#colnames(dat) <- c("Chrom", "Position", "End", "Ref", "Alt")
#x <- run_annovar(dat, c("ucsc_refgene", "exac03", "1000g2015aug_all"), buildver="hg19")

# Use ANNOVAR to annotate the given dataframe with the desired. Download using perl script as follows: 
# /home/local/ARCS/ak3792/Documents/Research/data/annovar/table_annovar.pl example/ex1.avinput humandb/ -buildver hg19 -out myanno -remove -protocol refGene,cytoBand,exac03,avsnp147,dbnsfp30a -operation gx,r,f,f,f -nastring . -csvout -polish -xref example/gene_xref.txt
#ak3792@c2b2ysld5:/data/annovar$ perl ./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar exac03 humandb
#NOTICE: Web-based checking to see whether ANNOVAR new version is available ... Done
#NOTICE: Downloading annotation database http://www.openbioinformatics.org/annovar/download/hg19_exac03.txt.gz ... OK
#NOTICE: Downloading annotation database http://www.openbioinformatics.org/annovar/download/hg19_exac03.txt.idx.gz ... OK
#NOTICE: Uncompressing downloaded files
run_annovar <- function(dat, annotations=NULL, buildver="hg19", run_in_shell=TRUE, allow_recursion=TRUE, work_folder=output_path("annovar_temp"), annovar_path="/data/annovar", database.dir="/data/annovar/humandb", database.cfg=system.file('extdata', 'config/databases.toml', package="annovarR"), allow_sqlite=FALSE) {
    # ANNOVAR search command
    if(length(dat) == 1 && class(dat) == "character") {
        library('annovarR')
        if(dat == "databases") { return(download.database(show.all.names = TRUE))
        } else if (dat == "annotations") { return(get.annotation.names())
        } else {
            dat_split <- strsplit(dat, " ")[[1]]
            dat_name <- get.download.name(dat_split[length(dat_split)])
            if(!allow_sqlite) { dat_name <- gsub("_sqlite", "", dat_name) }
            if(dat_split[1] == "download") {
                return(download.database(download.name=dat_name, buildver=buildver, database.dir=database.dir, db.type="txt"))
            } else { return(dat_name) }
        }
    }
    # Run ANNOVAR, either in shell through the system command (default), or with the annovarR R package (susceptible to memory issues/crashing).
    if(run_in_shell) {
        # Prepare input file to ANNOVAR
        dir.create(work_folder, showWarnings = FALSE)
        annovar_input_filename = full_path(work_folder,"dat.tsv")
        annovar_output_name = gsub(".tsv$", "", annovar_input_filename)
        annovar_operations <- sapply(annotations, function(x) return(c("f","g")[(x=="refGene")+1]))
        
        # If not yet done, add snv_indel annotation for quick lookup of whether a variant is an snv or an indel
        if(!("snv_indel" %in% colnames(dat))) {
            snv_indel <- rep("snv", nrow(dat)); snv_indel[is_indel(dat$Ref, dat$Alt)] <- "indel"
            dat <- cbind(dat, snv_indel)
        }
        
        #dat_annovar <- dat[!duplicated(dat[c("Chrom","Position","Ref","Alt","sample"),]),]
        hg19_only_annotations <- grepl("cadd", annotations)
        hg38_only_annotations <- grepl("gnomad30_genome", annotations)
        if(buildver=="hg38" && sum(hg19_only_annotations)>0) {
            if(allow_recursion && sum(!hg19_only_annotations)>0) {
                dat1 <- run_annovar(dat, annotations=annotations[!hg19_only_annotations], buildver=buildver, run_in_shell=run_in_shell, allow_recursion=allow_recursion, work_folder=work_folder, annovar_path=annovar_path, database.dir=database.dir, database.cfg=database.cfg, allow_sqlite=allow_sqlite)
                dat2 <- run_annovar(dat, annotations=annotations[hg19_only_annotations], buildver=buildver, run_in_shell=run_in_shell, allow_recursion=allow_recursion, work_folder=work_folder, annovar_path=annovar_path, database.dir=database.dir, database.cfg=database.cfg, allow_sqlite=allow_sqlite)
                return(merge(dat1, dat2[!duplicated(dat2[,c("Chrom","Position","Ref","Alt","sample")]),], all.x=TRUE))
            }
            dat_annovar <- liftover(dat, from="hg38", to="hg19", chr_colname="Chrom", start_colname="Position", end_colname=NULL, ref_colname="Ref", alt_colname="Alt", confirm_refseq=FALSE, mismatches_pause=FALSE)#[,c("Chrom","Position","Ref","Alt","sample","snv_indel","Chrom_hg38","Position_hg38")]
            buildver = "hg19"
            hg38_to_hg19_liftover = TRUE
            hg19_to_hg38_liftover = FALSE
            liftover_suffix = "_hg38"
        } else if(buildver=="hg19" && sum(hg38_only_annotations)>0) {
            if(allow_recursion && sum(!hg38_only_annotations)>0) {
                dat1 <- run_annovar(dat, annotations=annotations[!hg38_only_annotations], buildver=buildver, run_in_shell=run_in_shell, allow_recursion=allow_recursion, work_folder=work_folder, annovar_path=annovar_path, database.dir=database.dir, database.cfg=database.cfg, allow_sqlite=allow_sqlite)
                dat2 <- run_annovar(dat, annotations=annotations[hg38_only_annotations], buildver=buildver, run_in_shell=run_in_shell, allow_recursion=allow_recursion, work_folder=work_folder, annovar_path=annovar_path, database.dir=database.dir, database.cfg=database.cfg, allow_sqlite=allow_sqlite)
                return(merge(dat1, dat2[!duplicated(dat2[,c("Chrom","Position","Ref","Alt","sample")]),], all.x=TRUE))
            }
            dat_annovar <- liftover(dat, from="hg19", to="hg38", chr_colname="Chrom", start_colname="Position", end_colname=NULL, ref_colname="Ref", alt_colname="Alt", confirm_refseq=FALSE, mismatches_pause=FALSE)
            buildver = "hg38"
            hg19_to_hg38_liftover = TRUE
            hg38_to_hg19_liftover = FALSE
            liftover_suffix = "_hg19"
        } else {
            dat_annovar <- dat
            hg38_to_hg19_liftover = FALSE
            hg19_to_hg38_liftover = FALSE
            liftover_suffix = ""
        }
        dat_annovar <- unfactorize(dat_annovar)[,c("Chrom","Position","Position","Ref","Alt","sample","snv_indel",paste0("Chrom",liftover_suffix),paste0("Position",liftover_suffix),"Ref","Alt")]
        colnames(dat_annovar)[c(3,8,9,10,11)] <- c("Position_end", "Chrom_original", "Position_original", "Ref_original", "Alt_original")
        indel_indices <- dat_annovar$snv_indel == "indel"
        dat_annovar$Ref[indel_indices] <- gsub("^.", "", dat_annovar$Ref[indel_indices]); dat_annovar$Alt[indel_indices] <- gsub("^.", "", dat_annovar$Alt[indel_indices])
        dat_annovar$Ref[indel_indices][dat_annovar$Ref[indel_indices]==""] <- "-"
        deletion_indices <- dat_annovar$Alt[indel_indices]==""
        dat_annovar$Alt[indel_indices][deletion_indices] <- "-"
        dat_annovar$Position_end[indel_indices][deletion_indices] <- format(dat_annovar$Position[indel_indices][deletion_indices] + nchar(dat_annovar$Ref[indel_indices][deletion_indices]) - 1, , scientific=FALSE)
        dat_annovar$Position <- format(dat_annovar$Position, scientific=FALSE)
        dat_annovar$Position_end <- format(dat_annovar$Position_end, scientific=FALSE)
        dat_annovar$Position_original <- format(dat_annovar$Position_original, scientific=FALSE)
        
        write.table(dat_annovar, file=annovar_input_filename, row.names=FALSE, col.names=FALSE, quote=FALSE)
        messages <- system(paste0("perl ",annovar_path,"/table_annovar.pl '",annovar_input_filename,"' ",database.dir," -buildver ",buildver," -out ",annovar_output_name," -remove -protocol ",paste(annotations,collapse=",")," -operation ",paste(annovar_operations,collapse=",")," -nastring ."), intern=TRUE)
        annovar_invalid_input_filename = paste0(annovar_output_name,".invalid_input")
        annovar_output_filename = paste0(annovar_output_name,".",buildver,"_multianno.txt")
        
        dat_annotations <- read.csv(annovar_output_filename, sep="\t")
        dat_annotations <- cbind(dat_annovar, dat_annotations[,-which(colnames(dat_annotations) %in% colnames(dat_annovar))])
        colnames_to_revert_to_original <- colnames(dat_annotations)[grepl("_original",colnames(dat_annotations))]
        dat_annotations[,gsub("_original","",colnames_to_revert_to_original)] <- dat_annotations[,colnames_to_revert_to_original]
        dat_annotations <- dat_annotations[,-which(colnames(dat_annotations) %in% c("Chr","Start","End","Position_end",colnames_to_revert_to_original))]
        dat_annotations <- merge(dat, dat_annotations[!duplicated(dat_annotations[,c("Chrom","Position","Ref","Alt","sample")]),], all.x=TRUE)
    } else {
        library('annovarR')
        dat_annotations <- cbind(dat, annotation.merge(dat=dat[,c("Chrom", "Position", "Position", "Ref", "Alt")], anno.name=annotations, buildver=buildver, database.dir=database.dir, db.type='txt'))
    }
    cat("Done annotation!\n")
    return(unfactorize(dat_annotations))
}
#dat <- read.csv("/home/local/ARCS/ak3792/Documents/Research/data/WGS/re_called/re_called.tsv", sep=" ")
#a <- run_annovar(dat, annotations=c("refGene", "gnomad_genome", "exac03", "caddgt10"), buildver="hg38")
#dat_annotations <- annotation.merge(dat=cdh3[,c("Chrom", "Position", "Position", "Ref", "Alt")], anno.name=c("ucsc_refgene"), buildver="hg19", database.dir="/data/annovar/humandb", db.type='txt')
#esp6500siv2_all: alternative allele frequency in All subjects in the NHLBI-ESP project with 6500 exomes, including the indel calls and the chrY calls.
#gnomad_genome
# 1000g2015aug_all; ALL.sites.2015_08 is 1000g2015aug
#"ucsc_refgene"
#cdh3 <- cdh3[,c(1:6)]
#for(annot in c("exac03", "1000g2015aug_all")) {
#    print(annot)
#    cdh3 <- run_annovar(cdh3, annot)
#}
#sum(!is.na(cdh3$ExAC_ALL))
#as.numeric(cdh3$ExAC_ALL[!is.na(cdh3$ExAC_ALL)]) > 0.001

# List all functions available in alex_suite.R, except for this one. The function alex_suite_functions is an alias of this alex_suite function, so either can be used.
alex_suite <- function(search=NULL) {
    alex_suite_lines <- readLines("alex_suite.R")
    function_line_indices <- grep("^[^\\s]+ <- function\\(.*\\)", alex_suite_lines, perl=TRUE)
    comment_line_indices <- grep("^#", alex_suite_lines, perl=TRUE)
    
    if (is.null(search)) { 
        cat("Welcome to Alex Suite!\n") 
        potential_starter_comment_line_index = 1
        while(potential_starter_comment_line_index %in% comment_line_indices) {
            cat(gsub("^# *", "", alex_suite_lines[potential_starter_comment_line_index]), "\n")
            potential_starter_comment_line_index = potential_starter_comment_line_index + 1
        }
        if (potential_starter_comment_line_index > 1) { cat("\n") } # Starter comment was found and printed, so make a newline space before printing functions.
        cat("Functions: \n\n")
    } else {
        cat(paste0("Alex Suite functions containing keyword \"", search,"\":\n\n")) 
    }
    
    for(function_line_index in function_line_indices) {
        f <- gsub(" *\\{.*", "", alex_suite_lines[function_line_index])
        if (f == "alex_suite <- function()") { next } # skip this alex_suite() information function
        if (!is.null(search) && !grepl(search, strsplit(f, " <- ")[[1]][1])) { next }
        
        function_comment_lines <- c()
        potential_comment_line_index = function_line_index - 1
        while(potential_comment_line_index %in% comment_line_indices) {
            function_comment_lines <- c(potential_comment_line_index, function_comment_lines)
            potential_comment_line_index = potential_comment_line_index - 1
        }
        if (length(function_comment_lines) > 0) { sapply_out <- sapply(alex_suite_lines[function_comment_lines], function(x) { cat(x, "\n") }) }
        cat(f, "\n\n")
    }
    #cat(paste0(gsub("\\{.*", "", alex_suite_lines[function_line_indices]), collapse="\n"))
}
alex_suite_functions <- alex_suite

# returns string w/o leading or trailing whitespace
trim <- function (x) {
    return(gsub("^\\s+|\\s+$", "", x))
}

# Return the genes closest to each of the respective specified genomic positions.
# Requires GenomicFeatures package, from Bioconductor
get_closest_genes <- function(chromosomes, positions) { 
    
}

# Return clusterings for the given gene annotation, based on annotations containing multiple genes.
# Requires igraph package
cluster_overlapping_genes <- function(genes, verbose=TRUE) { 
    genes <- paste0(genes)
    if(verbose) { print("Clustering overlapping genes...") }
    gene_pairs <- data.frame(stringsAsFactors=FALSE)
    for(i in 1:length(genes)) {
        if(verbose) { print(i) }
        all_genes <- unlist(strsplit(paste0(genes[i]), ","))
        gene_pairs <- rbind(gene_pairs, expand.grid(all_genes, all_genes)) 
    }
    g <- graph.data.frame(unique(gene_pairs), directed=FALSE)
    comps <- components(g) 
    num_comps <- comps$no
    if(verbose) { print("Relabeling gene components...") }
    for(comp in 1:num_comps) {
        if(verbose) { print(comp) }
        connected_genes <- paste0(names(V(g)[which(comps$membership == comp)]))
        connected_genes_string <- paste0(sort(connected_genes), collapse=",")
        genes_indices_to_relabel <- sapply(strsplit(genes, ","), function(x) { sum(x %in% connected_genes) > 0 })
        genes[genes_indices_to_relabel] <- connected_genes_string
    }
    return(genes)
}

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

# Return histone modification function for the given histone mark
get_histone_mark_function <- function(histone_mark) {
    info <- read.csv(data_path("histone_mark_information.csv"), skip=1)
    if(grepl("me", histone_mark)) { 
        type = "methylation"
        histone_mark_split <- strsplit(histone_mark, "me")[[1]]
        aa <- histone_mark_split[1] # amino acid
        modification <- paste0("me", histone_mark_split[2]) # modification
    } else if (grepl("ac", histone_mark)) { 
        type = "acetylation"
        histone_mark_split <- strsplit(histone_mark, "ac")[[1]]
        aa <- histone_mark_split[1] # amino acid
        modification <- paste0("ac", histone_mark_split[2]) # modification
    } else { return("unknown") }
    
    info <- info[info$modification_type == type,] # filter info data frame by appropriate modification type.
    functions <- unique(c(unlist(strsplit(paste0(info$proposed_function[paste0(info$histone, info$site) == aa]), ", ")), "unknown"))
    
    # Check if there is specific information on the methylation in any proposed functions; if there is, remove functions specified for NOT this modification.
    if (type == "methylation" ) {
        me1_hits <- unlist(sapply(functions, function(fun) grepl("mono-Me", fun))) 
        me2_hits <- unlist(sapply(functions, function(fun) grepl("di-Me", fun))) 
        me3_hits <- unlist(sapply(functions, function(fun) grepl("tri-Me", fun)))
        if (modification == "me1") { # if(sum(me1_hits) > 0) { functions <- gsub(" \\(mono-Me\\)", "", functions[me1_hits]) } 
            functions <- gsub(" \\(mono-Me\\)", "", functions[!(me2_hits | me3_hits)])
        } else if (modification == "me2") { 
            functions <- gsub(" \\(di-Me\\)", "", functions[!(me1_hits | me3_hits)])
        } else if (modification == "me3") { 
            functions <- gsub(" \\(tri-Me\\)", "", functions[!(me1_hits | me2_hits)])
        }
    }
    functions <- unique(functions)
    
    # Remove "unknown" from function, unless it is the only proposed function.
    if (length(functions) > 1) { functions <- functions[functions != "unknown"] }
        
    return(paste0(functions, collapse=", "))
}

# Get relevant Roadmap EIDs for the specified disease.
get_relevant_roadmap_eids <- function(disease) {
    roadmap_eids <- read.csv(data_path("roadmap_eids.csv"))
    if (grepl("CDH", disease) || grepl("SSC", disease)) { 
        return(paste0(roadmap_eids$EID)[roadmap_eids$GROUP %in% c("IMR90", "ESC") | roadmap_eids$standardized_name %in% c("Fetal Lung", "Fetal Heart")]) # "iPSC"
    } else if (grepl("CHD", disease)) { 
        return(paste0(roadmap_eids$EID)[roadmap_eids$GROUP %in% c("ESC") | roadmap_eids$standardized_name %in% c("Fetal Heart")]) # "iPSC"
    } else if (grepl("ASD|autism", disease)) { 
        return(paste0(roadmap_eids$EID)[roadmap_eids$GROUP %in% c("Brain", "ESC")]) # "iPSC"
    } else { return(c()) }
}

# Store Roadmap EID name mappings.
store_roadmap_eid_names <- function() {
    roadmap_eid_names <<- new.env()
    roadmap_eids <- read.csv(data_path("roadmap_eids.csv"))
    for(i in 1:nrow(roadmap_eids)) {
        roadmap_eid_names[[paste(roadmap_eids$EID[i])]] <<- paste(roadmap_eids$standardized_name[i])
    }
    set_global(roadmap_eid_names)
}
#store_roadmap_eid_names()

# Return Roadmap epigenome names for the given EIDs vector.
get_roadmap_epigenome_names <- function(eids) {
    epigenome_names <- c()
    for(i in 1:length(eids)) {
        epigenome_name <- roadmap_eid_names[[paste(eids[i])]]
        if(is.null(epigenome_name)) {
            epigenome_names <- c(epigenome_names, "-")
        } else {
            epigenome_names <- c(epigenome_names, epigenome_name)
        }
    }
    return(epigenome_names)
}
#epigenome_names <- get_roadmap_epigenome_names(unlist(lapply(seq(1:length(TSS_burdens$feature)), function(i) { return(strsplit(paste(TSS_burdens$feature[i]), "\\.")[[1]][1]) } )))

# Combine two datasets by shared feature names, and return the combined dataset.
combine_datasets <- function(dat1, dat2, use_temp_file=FALSE) {
    shared_features <- intersect(colnames(dat1), colnames(dat2))
    dat1 <- dat1[,shared_features]
    dat2 <- dat2[,shared_features]
    if (use_temp_file) {
        write.table(dat1, file="temp.tsv", sep="\t", row.names=FALSE)
        write.table(dat2, file="temp.tsv", sep="\t", row.names=FALSE, append=TRUE)
        dat_combined <- read.table("temp.tsv", sep="\t", header=TRUE)
    } else { dat_combined <- rbind(dat1, dat2) }
    return(dat_combined)
}

get_numeric <- function(vec) {
    nums <- as.numeric(paste0(vec[vec != "."]))
    #nums <- nums[!is.na(nums)]
    if (sum(is.na(nums)) > 0) { 
        return(NA)  
    } else {
        return(nums)
    }
}

# Get and return CDH mouse candidate genes.
get_cdh_mouse_candidate_genes <- function() {
    dat <- read.csv(data_path("Table S5 CDH training genes with mouse model.csv"), header=TRUE, skip=1)
    candidate_genes <- new.env()
    for(gene in dat$Gene) {
        candidate_genes[[paste(gene)]] <- 1
    }
    return(candidate_genes)
}

# Return results of word search over the column names of the given data frame. 
search_colnames <- function(query, dat) {
    if (length(query) > 0) {
        return(unique(paste0(unlist(sapply(1:length(query), function(i) { colnames(dat)[grepl(query[i], colnames(dat))] } )))))
    } else { (return(c())) }
}


# Returns whether or not each entry specified by first parameter, split by delim, intersects the second parameter.
includes_at_least_one <- function(query, v, delim=",") { 
    if (length(query) > 0) {
        return(sapply(1:length(query), function(i) { sum(strsplit(paste0(query[i]), delim)[[1]] %in% v) > 0 }))
    } else { return(0) }
}

# Updates OS being currently worked on, plus relevant environment variables. 
update_os <- function() {
    sysname = Sys.info()[["sysname"]]
    nodename = Sys.info()[["nodename"]]
    # Delimiter between folders/files in filepaths. Should be "/" for Linux and "\\" for Windows.
    if (sysname == "Windows") { FILEPATH_DELIM <<- "\\" } else { FILEPATH_DELIM <<- "/" } 
    if (nodename == "c2b2ysld5") { 
        # Path to data folder, allowing simpler loading of data with the alex_suite.R data_path function.
        DATA_FOLDER <<- "/home/local/ARCS/ak3792/Documents/Research/data"
        # Path to big data folder, allowing simpler loading of big data with the alex_suite.R bigdata_path function.
        BIGDATA_FOLDER <- "/mnt/data/ak3792"
        # Path to output folder, allowing simpler writing of output with the alex_suite.R output_path function.
        OUTPUT_FOLDER <<- "/home/local/ARCS/ak3792/Documents/Research/ML/output"
    } else if (nodename == "c2b2ysld1") { 
        # Path to data folder, allowing simpler loading of data with the alex_suite.R data_path function.
        DATA_FOLDER <<- "/home/local/ARCS/ak3792/machine/Research/data"
        # Path to output folder, allowing simpler writing of output with the alex_suite.R output_path function.
        OUTPUT_FOLDER <<- "/home/local/ARCS/ak3792/machine/Research/WGS/output"
    } else if (nodename == "GANYMEDE") { 
        # Path to data folder, allowing simpler loading of data with the alex_suite.R data_path function.
        DATA_FOLDER <<- "E:\\Documents\\Columbia\\Research\\data"
        # Path to output folder, allowing simpler writing of output with the alex_suite.R output_path function.
        OUTPUT_FOLDER <<- "E:\\Documents\\Columbia\\Research\\WGS\\output"
    }
    cat(paste0("Updated environment for ", nodename, " (", sysname, "):"), "\n")
    cat(paste0("    FILEPATH_DELIM: \"", gsub("\\\\", "\\\\\\\\", FILEPATH_DELIM), "\""), "\n")
    cat(paste0("    DATA_FOLDER: ", gsub("\\\\", "\\\\\\\\", DATA_FOLDER)), "\n")
    cat(paste0("    OUTPUT_FOLDER: ", gsub("\\\\", "\\\\\\\\", OUTPUT_FOLDER)), "\n")
}

# Return full path, given a directory and path from that directory. Simple utility method.
full_path <- function(directory, path_from_directory) {
    return_string <- paste0(directory, FILEPATH_DELIM, path_from_directory)
    return_string <- gsub("/|\\\\", gsub("\\\\", "\\\\\\\\", FILEPATH_DELIM), return_string)
    return(return_string) 
}

# Return full data path, using DATA_FOLDER global variable if it is specified in the current R workspace.
data_path <- function(path_from_data_folder="") { 
    path_from_data_folder <- gsub("/|\\\\", gsub("\\\\", "\\\\\\\\", FILEPATH_DELIM), path_from_data_folder)
    return_string <- path_from_data_folder
    if(exists("DATA_FOLDER")) { return_string <- paste0(DATA_FOLDER, FILEPATH_DELIM, path_from_data_folder) }
    return(return_string) 
}

# Return full output path, using OUTPUT_FOLDER global variable if it is specified in the current R workspace.
output_path <- function(path_from_output_folder="") { 
    path_from_output_folder <- gsub("/|\\\\", gsub("\\\\", "\\\\\\\\", FILEPATH_DELIM), path_from_output_folder)
    return_string <- path_from_output_folder
    if(exists("OUTPUT_FOLDER")) { return_string <- paste0(OUTPUT_FOLDER, FILEPATH_DELIM, path_from_output_folder) }
    return(return_string) 
}

# Draw density ggplot of the given data frame (with x and y columns)
draw_plot <- function(DF, hex_density = 10, title="", xlab="log(TPM+1)", ylab="log(# peaks + 1)", legend_text="# genes", cor_method="Pearson", linear_best_fit=TRUE, quadratic_best_fit=TRUE, identity_diagonal=FALSE, ignore=NA, filename=NULL) {
    p = ggplot(DF, aes(x=x, y=y)) + stat_binhex(bins=hex_density) + scale_fill_gradient(name=paste0(cor_method,"\nCorr:\n", round(cor(DF[,1], DF[,2], method=tolower(cor_method)), 3),paste0("\n\n",legend_text)), trans="log", breaks=10^(0:4)) + labs(title=title, x=xlab, y=ylab) + 
        theme(
            text = element_text(size=22),
            plot.title = element_text(size=16, face="bold", hjust=0.5, vjust=1),
            axis.title.x = element_text(size=16, face="bold", hjust=0.5, vjust=0.5),
            axis.title.y = element_text(size=16, face="bold", hjust=0.5, vjust=1),
            legend.title=element_text(size=16, face="bold"),
            legend.text=element_text(size=16)
        )
    if(identity_diagonal) { p = p + geom_segment(aes(x=max(DF), y=max(DF), xend=min(DF), yend=min(DF))) }
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

# Convert pdf to png, with the specified pixels per inch (ppi) resolution
pdf_to_png <- function(input_filepath, ppi=500, output_folder=NULL, output_filepath=NULL) {
    input_filename <- unlist(strsplit(input_filepath, FILEPATH_DELIM)); input_filename = input_filename[length(input_filename)]
    if(is.null(output_filepath)) { 
        filename = gsub("\\.pdf$",".pdf", input_filename)
        if(is.null(output_folder)) { output_filepath <- output_path(paste0("images/",filename))
        } else {
            output_folder <- output_path(paste0("images/",output_folder))
            dir.create(output_folder, showWarnings = FALSE)
            output_filepath = paste0(output_folder,"/",filename) 
        }
    }
    system(paste0("pdftoppm -r ",ppi," -png ",input_filepath," ",output_filepath))
    print(paste0("PDF converted and written to ",output_filepath))
}

# Return as string the name of the given variable.
varname <- function(variable) {
    return(deparse(substitute(variable)))
}

# If aggregate is true, aggregates other columns together for duplicates. Otherwise, just ignores duplicates and uses first instance of each variant.
# When aggregate is true, unique argument determines whether to delimit identical values, while delim argument is used to choose character for paste0 collapse. 
unique_variants <- function(variants, aggregate=TRUE, unique=TRUE, delim=",") {
	if(aggregate) {
		aggregated_variants <- aggregate(variants, by=list(variants$Chrom, variants$Position, variants$Ref, variants$Alt, variants$sample), FUN=function(x){ if(unique){x <- unique(x)}; paste0(x,collapse=delim)} )
		return(aggregated_variants[,colnames(variants)])
	} else {
		variant_strings <- apply(variants[,c("Chrom", "Position", "Ref", "Alt", "sample")], 1, paste0, collapse=";")
		variant_strings_duplicates_indices <- duplicated(variant_strings)
		return(variants[!variant_strings_duplicates_indices,])
	}
}

get_ensembl_name_mappings <- function() {
	ensembl_dat <- read.csv(data_path("EnsemblID_genename.txt"), sep="\t")
	ensembl_genes <- new.env()
	sapply_out <- sapply(seq(1:nrow(ensembl_dat)), function(i) { ensembl_genes[[paste(ensembl_dat$Gene.stable.ID[i])]] <- paste(ensembl_dat$Gene.name[i]) } )
	return(ensembl_genes)
}

# The constraint parameter should be in form like "y>x", where "y" is the colname and x is value. Can use "<=" and any other inequalities too.
# The return_type parameter can be set to "transcript", or the default "gene"
get_constrained_genes <- function(constraint, return_type="gene", exclude_transcript_suffix=TRUE) {
	exac_dat <- read.csv(data_path("fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt"), sep="\t")
	x <- paste0("exac_dat$", constraint)
	constrained_indices <- which(eval(parse(text=x)))
	constraint_variable <- strsplit(constraint, " |=|<|>", perl=TRUE)[[1]][1]

	constrained <- new.env()
	for(constrained_index in constrained_indices) {
		key <- paste0(exac_dat[constrained_index, return_type])
		if(return_type == "transcript" && exclude_transcript_suffix) { key <- strsplit(key, "\\.")[[1]][1] }
		constrained[[key]] <- exac_dat[constrained_index, constraint_variable]
	}

	return(constrained)
}

# Get bivalent genes
get_bivalent_genes <- function() {
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
    biv_genes <- unique(bivalent_genes$gene)
    return(biv_genes)
}

# Annotate nearest gene
get_nearest_genes <- function(gr) {
    if(class(gr) != "GRanges") { gr <- to_genomic_regions(gr) }
    genebody <- read.csv(data_path("refGene\\refGene_hg19_genebody_fixed.bed"), sep="\t", stringsAsFactors=FALSE)
    genebody_granges <- to_genomic_regions(genebody, chr_colname="chromosome", start_colname="tss", end_colname="tes", strand_colname="strand", label_colname="gene", order_coordinates=TRUE, remove_duplicate_labels=TRUE)
    nearest_genes <- names(genebody_granges)[nearest(gr, genebody_granges)]
    return(nearest_genes)
}

# Gets GTEX tissues and HepG2 and K562 log(TPM+1) gene expressions matrix
get_gene_expressions_matrix <- function(sources=c("K562","HepG2","GTEx"), log_scale=TRUE, fill_empty=TRUE, unique_genes=TRUE, include_refseq=TRUE) {
    a <- read.csv(data_path("EnsemblID_genename.txt"), sep="\t") # Ensembl genes
    if(include_refseq) {
        genebody <- read.csv(data_path("refGene\\refGene_hg19_genebody_fixed.bed"), sep="\t", stringsAsFactors=FALSE)
        genebody <- genebody[!(genebody$gene %in% a$Gene.name),]
        genebody$chromosome <- gsub("chr","",genebody$chromosome)
        genebody <- genebody[,c("tss","tes","transcript","chromosome","gene")]
        colnames(genebody) <- colnames(a)
        a <- rbind(a, genebody)
    }

    sources_queried <- c("K562","HEPG2","GTEX") %in% toupper(sources)

    if(sources_queried[1]) { # K562
        genes_k562_tpms <- read.csv(data_path("ENCFF047WAI_K562_gene_quantifications.tsv"), sep="\t")
        genes_k562_tpms <- cbind(gsub("\\..*$","",genes_k562_tpms$gene_id), genes_k562_tpms)
        colnames(genes_k562_tpms)[1] <- "Gene.stable.ID"
        colnames(genes_k562_tpms)[which(colnames(genes_k562_tpms) == "TPM")] <- "K562_TPM"
        a <- merge(a, genes_k562_tpms, all.x=TRUE)
    }
    if(sources_queried[2]) { # HepG2
        genes_hepg2_tpms <- read.csv(data_path("ENCFF945LNB_HepG2_gene_quantifications.tsv"), sep="\t")
        genes_hepg2_tpms <- cbind(gsub("\\..*$","",genes_hepg2_tpms$gene_id), genes_hepg2_tpms)
        colnames(genes_hepg2_tpms)[1] <- "Gene.stable.ID"
        colnames(genes_hepg2_tpms)[which(colnames(genes_hepg2_tpms) == "TPM")] <- "HepG2_TPM"
        a <- merge(a[,c("Gene.stable.ID","Gene.Start..bp.","Gene.End..bp.","Chromosome.scaffold.name","Gene.name","gene_id","transcript_id.s.","length","effective_length","expected_count","K562_TPM")], genes_hepg2_tpms[,c("Gene.stable.ID","HepG2_TPM")], all.x=TRUE)
    }
    if(sources_queried[3]) { # GTEx
        genes_gtex_tpms <- read.table(data_path("GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct"), skip=2, header=3, sep="\t")
        genes_gtex_tpms <- cbind(gsub("\\..*$","",genes_gtex_tpms$Name), genes_gtex_tpms)
        colnames(genes_gtex_tpms) <- c("Gene.stable.ID", "gene_id", "Gene.name", paste0(gsub("\\.$","",gsub("\\.+",".",colnames(genes_gtex_tpms)[-c(1:3)])),"_TPM"))
        a <- merge(a, genes_gtex_tpms, all.x=TRUE)
    }

    # Polish columns as specified, converting to logscale and/or filling empty cells with mean gene expression.
    a_tpm_columns <- grepl("TPM",colnames(a))
    if(sum(a_tpm_columns)>0 && (log_scale || fill_empty)) {
        for(i in which(a_tpm_columns)) {
            na_exp <- is.na(a[,i])
            a[,i][na_exp] <- suppressWarnings(mean(a[,i][!na_exp]))
            if(log_scale) { a[,i] <- suppressWarnings(log(a[,i] + 1)) }
            if(!fill_empty) { a[,i][na_exp] <- NA }
        }
    }

    a <- a[,c("gene_id","Gene.stable.ID","Gene.name","Chromosome.scaffold.name","Gene.Start..bp.","Gene.End..bp.","transcript_id.s.","length","effective_length","expected_count",colnames(a)[a_tpm_columns])]
    colnames(a)[1:7] <- c("ensembl_id","ensembl_stable_id","gene","chromosome","start","end","transcript_ids")
    if(log_scale) { colnames(a) <- gsub("TPM$", "logTPM", colnames(a)) }

    if(unique_genes) { 
        a <- a[!duplicated(a$gene),] 
        rownames(a) <- a$gene
    }

    return(a)
}


#################################################################################################################
# Define top (HIGH_EXPRESSION_PERCENTILE)% HDE and HHE genes, and store the respective genebodies.
# Forces number of HDE and HHE genes to be the same (NUM_HE_GENES) for better cross-disease result comparison.
# H = top 25%, L = bottom 25%, MH and ML are respective middle quartiles (50-75th, 25-50th), and HMH is top 50%.
# The default return_type of "genes" returns all final gene vectors. Also possible is return_type of "heart_rank" or "diaphragm_rank".
#################################################################################################################
get_genes_by_expression <- function(return_type="genes") {
    HIGH_EXPRESSION_PERCENTILE = 25
	human_genes <- read.csv(data_path("human_genesize.csv"))$gene # Can use to filter for only human genes.
	genebody <- read.csv(data_path("refGene\\refGene_hg19_genebody_fixed.bed"), sep="\t", stringsAsFactors=FALSE)
	HDE_rank_dat <- read.csv(data_path("diaphragm_rank.csv"), stringsAsFactors=FALSE)
	HDE_rank_dat <- HDE_rank_dat[!duplicated(HDE_rank_dat$gene) & HDE_rank_dat$gene %in% human_genes & HDE_rank_dat$gene %in% genebody$gene,]
	HDE_rank_dat <- HDE_rank_dat[order(HDE_rank_dat$rank),]
	HHE_rank_dat <- read.csv(data_path("heart_brain_rank.csv"), stringsAsFactors=FALSE)
	HHE_rank_dat <- HHE_rank_dat[!duplicated(HHE_rank_dat$human.External.Gene.Name) & HHE_rank_dat$human.External.Gene.Name %in% human_genes & HHE_rank_dat$human.External.Gene.Name %in% genebody$gene,]
	HHE_rank_dat$e14.5_rank <- 100 - HHE_rank_dat$e14.5_rank # HHE and HDE rankings are backwards!!!
	HHE_rank_dat <- HHE_rank_dat[order(HHE_rank_dat$e14.5_rank),]
	NUM_HE_GENES = floor(min(c(nrow(HDE_rank_dat)*HIGH_EXPRESSION_PERCENTILE/100, nrow(HHE_rank_dat)*HIGH_EXPRESSION_PERCENTILE/100)))
	#HDE_genes <- unique(HDE_rank_dat$gene[HDE_rank_dat$rank <= HIGH_EXPRESSION_PERCENTILE])
	#HHE_genes <- unique(HHE_rank_dat$human.External.Gene.Name[HHE_rank_dat$e14.5_rank <= HIGH_EXPRESSION_PERCENTILE])
	HDE_genes <- HDE_rank_dat$gene[1:NUM_HE_GENES] 
	HHE_genes <- HHE_rank_dat$human.External.Gene.Name[1:NUM_HE_GENES]
	HMHDE_genes <- HDE_rank_dat$gene[1:(2*NUM_HE_GENES)]
	HMHHE_genes <- HHE_rank_dat$human.External.Gene.Name[1:(2*NUM_HE_GENES)]
	MHDE_genes <- HDE_rank_dat$gene[(NUM_HE_GENES+1):(2*NUM_HE_GENES)] 
	MHHE_genes <- HHE_rank_dat$human.External.Gene.Name[(NUM_HE_GENES+1):(2*NUM_HE_GENES)]
	LDE_genes <- HDE_rank_dat$gene[(nrow(HDE_rank_dat)-NUM_HE_GENES+1):nrow(HDE_rank_dat)] 
	LHE_genes <- HHE_rank_dat$human.External.Gene.Name[(nrow(HHE_rank_dat)-NUM_HE_GENES+1):nrow(HHE_rank_dat)]
	MLDE_genes <- HDE_rank_dat$gene[(nrow(HDE_rank_dat)-(2*NUM_HE_GENES)+1):(nrow(HDE_rank_dat)-NUM_HE_GENES)] 
	MLHE_genes <- HHE_rank_dat$human.External.Gene.Name[(nrow(HHE_rank_dat)-(2*NUM_HE_GENES)+1):(nrow(HHE_rank_dat)-NUM_HE_GENES)]

	return_env <- new.env()
	if(return_type == "heart_rank") {
	    all_genes <- c(HHE_genes, MHHE_genes, MLHE_genes, LHE_genes)
	    for(i in 1:length(all_genes)) { return_env[[all_genes[i]]] <- i/length(all_genes) }
	} else if(return_type == "diaphragm_rank") {
	    all_genes <- c(HDE_genes, MHDE_genes, MLDE_genes, LDE_genes)
	    for(i in 1:length(all_genes)) { return_env[[all_genes[i]]] <- i/length(all_genes) }
	} else { # return_type == "genes"
    	return_env <- new.env()
    	return_env[["HDE"]] <- HDE_genes; return_env[["HHE"]] <- HHE_genes
    	return_env[["HMHDE"]] <- HMHDE_genes; return_env[["HMHHE"]] <- HMHHE_genes
    	return_env[["MHDE"]] <- MHDE_genes; return_env[["MHHE"]] <- MHHE_genes
    	return_env[["MLDE"]] <- MLDE_genes; return_env[["MLHE"]] <- MLHE_genes
    	return_env[["LDE"]] <- LDE_genes; return_env[["LHE"]] <- LHE_genes
	} 
	return(return_env)
}

# Binomial test: Null is binom(m, p), where m = m1 + m0 (m0 is the number of such variants in controls),  p = n1/(n1+n0),   
# n1 is the number of cases in total, n0 is the number of controls in total
binomial_test <- function(cases_Y, controls_Y, sample_count, control_sample_count, alternative=c("two.sided")) {
    m0 = controls_Y
	m1 = cases_Y
	n0 = control_sample_count	
	n1 = sample_count
	if (n0 < 1 || n1 < 1 || m1 < 1) { 
		binom_enrichment = 0
		binom_p.value = 1
	} else {
		binom_size = m1+m0
		binom_prob = n1/(n1+n0)
		binom_enrichment = (m1/n1)/(m0/n0)#m1/(binom_size*binom_prob)
		#binom_p.value = pbinom(q=m1-1, size=binom_size, prob=binom_prob, lower.tail=FALSE)
		binom_p.value = binom.test(m1, binom_size, binom_prob, alternative=alternative)$p.value
	} 

	result <- new.env()
	result[["p.value"]] <- binom_p.value
	result[["estimate"]] <- binom_enrichment
	result[["m0"]] <- m0
	result[["m1"]] <- m1
	result[["n0"]] <- n0
	result[["n1"]] <- n1
	return(result)
}

# Fisher's exact test: Null is m1/n1 = m0/n0, where m1 and m0 are the number of hits in cases and controls, respectively, and n1 and n0 are the total number of variants in cases and controls, respectively.
# n1 is the number of case variants in total, n0 is the number of control variants in total
fisher_exact_test <- function(cases_Y, controls_Y, variant_count, control_variant_count, alternative=c("two.sided")) {
    m0 = as.numeric(paste0(controls_Y))
    m1 = as.numeric(paste0(cases_Y))
    n0 = as.numeric(paste0(control_variant_count))	
    n1 = as.numeric(paste0(variant_count))
    if (n0 < 1 || n1 < 1 || m1 < 1) { 
        fet_enrichment = 0
        fet_ci = c(-Inf,Inf)
        fet_p.value = 1
    } else {
        fet_result <- fisher.test(matrix(c(m1, n1-m1, m0, n0-m0), nrow = 2, dimnames = list(hits = c("Y", "N"), status = c("case", "control"))), alternative = alternative)
        fet_enrichment = fet_result$estimate
        fet_ci = fet_result$conf.int
        fet_p.value = fet_result$p.value
    } 
    
    result <- new.env()
    result[["p.value"]] <- fet_p.value
    result[["estimate"]] <- fet_enrichment
    result[["ci"]] <- fet_ci
    result[["m0"]] <- m0
    result[["m1"]] <- m1
    result[["n0"]] <- n0
    result[["n1"]] <- n1
    return(result)
}

# Poisson exact test: Null is m1/n1 = m0/n0, where m1 and m0 are the number of hits in cases and controls, respectively, and n1 and n0 are the total number of variants in cases and controls, respectively.
# n1 is the number of case variants in total, n0 is the number of control variants in total
poisson_exact_test <- function(cases_Y, controls_Y, variant_count, control_variant_count, alternative=c("two.sided")) {
    m0 = controls_Y
    m1 = cases_Y
    n0 = control_variant_count	
    n1 = variant_count
    if (n0 < 1 || n1 < 1 || m1 < 1) { 
        fet_enrichment = 0
        fet_p.value = 1
    } else {
        fet_result <- poisson.test(matrix(c(m1, n1-m1, m0, n0-m0), nrow = 2, dimnames = list(cases = c("Y", "N"), control = c("Y", "N"))), alternative = alternative)
        fet_enrichment = fet_result$estimate
        fet_p.value = fet_result$p.value
    } 
    
    result <- new.env()
    result[["p.value"]] <- fet_p.value
    result[["estimate"]] <- fet_enrichment
    result[["m0"]] <- m0
    result[["m1"]] <- m1
    result[["n0"]] <- n0
    result[["n1"]] <- n1
    return(result)
}

genomic_regions_to_dat <- function(regions_granges, chr_colname="Chrom", pos_colname="Position", start_colname=NULL, end_colname=NULL, strand_colname="strand", label_colname="label", add_chr_prefix=FALSE) {
    dat_colnames <- colnames(dat)
    if(is.null(start_colname) || !(start_colname %in% dat_colnames) || is.null(end_colname) || !(end_colname %in% dat_colnames)) { 
        genomic_dat <- data.frame(seqnames(regions_granges), start(regions_granges), strand(regions_granges), names(regions_granges))
        colnames(genomic_dat) <- c(chr_colname, pos_colname, strand_colname, label_colname)
    } else {
        genomic_dat <- data.frame(seqnames(regions_granges), start(regions_granges), end(regions_granges), strand(regions_granges), names(regions_granges))
        colnames(genomic_dat) <- c(chr_colname, start_colname, end_colname, strand_colname, label_colname)
    }
    if(add_chr_prefix) { genomic_dat[,chr_colname] <- paste0("chr", gsub("chr", "", genomic_dat[,chr_colname])) }
    return(genomic_dat)
}

genomic_dat_to_strings <- function(dat, chr_colname="Chrom", pos_colname="Position", start_colname=NULL, end_colname=NULL, sample_colname=NULL, ref_colname="Ref", alt_colname="Alt", notes_colname=NULL, add_chr_prefix=TRUE, sep=NULL, empty_string="-") {
  if(!is.null(dat) && nrow(dat)>0) { 
    dat_colnames <- colnames(dat)
    if(is.null(start_colname) || !(start_colname %in% dat_colnames)) { start_colname <- pos_colname }
    dat_chromosomes <- gsub("chr", "", dat[,chr_colname]); if(add_chr_prefix) { dat_chromosomes <- paste0("chr", dat_chromosomes) }
    if(is.null(sample_colname)) { dat_variants <- paste0(dat_chromosomes,":",dat[,start_colname]) } else { dat_variants <- paste0(dat[,sample_colname],"_",dat_chromosomes,":",dat[,start_colname]) }
    if(!is.null(end_colname) && end_colname %in% dat_colnames) { dat_variants <- paste0(dat_variants,"-",dat[,end_colname]) }
    if(!is.null(ref_colname) && ref_colname %in% dat_colnames && !is.null(alt_colname) && alt_colname %in% dat_colnames) { dat_variants <- paste0(dat_variants,"_",dat[,ref_colname],">",dat[,alt_colname]) }
    if(!is.null(notes_colname)) { dat_variants <- paste0(dat_variants, "_(", dat[,notes_colname], ")", collapse=sep) } else { dat_variants <- paste0(dat_variants, collapse=sep) }
    return(dat_variants)
  } else {return(empty_string) }
}
#genomic_dat_to_strings(chdfb_combined_muts_wgsa[1:10,], chr_colname="X.chr", pos_colname="pos", ref_colname="ref", alt_colname="alt", notes_colname="snv_indel")

genomic_strings_to_dat <- function(genomic_strings, column_names=c("sample", "Chrom", "Position", "Ref", "Alt", "notes"), column_order=NULL, add_chr_prefix=TRUE, sep=";") {
    num_cols = length(column_names)
    dat <- t(data.frame(sapply(unlist(strsplit(gsub("\\(|\\)", "", paste0(genomic_strings)), sep)), function(x) unlist(strsplit(x, "_|:|>"))[1:num_cols]))); rownames(dat) <- NULL
    dat[,2] <- gsub("chr", "", dat[,2]); if(add_chr_prefix) { dat[,2] <- paste0("chr", dat[,2]) }
    colnames(dat) <- column_names
    if(!is.null(column_order)) { dat <- dat[,column_order[1:num_cols]] }
    dat <- unfactorize(data.frame(dat))
    return(dat)
}

genomic_coordinates_to_strings <- function(chromosomes, starts, ends=NULL, add_chr_prefix=FALSE) {
	if(!is.null(ends)) { ends <- paste0("-", ends) }
	chromosomes <- gsub("chr", "", chromosomes)
	if(add_chr_prefix) { chromosomes <- paste0("chr", chromosomes) }
	return(paste0(chromosomes, ":", starts, ends))
}

genomic_coordinates_from_strings <- function(genomic_coordinate_strings, keep_chr_prefix=FALSE) {
	if(!keep_chr_prefix) { genomic_coordinate_strings <- gsub("chr", "", genomic_coordinate_strings) }
	genomic_coordinates_dat <- data.frame(t(unname(data.frame(strsplit(genomic_coordinate_strings, ":|-")))))
	if(ncol(genomic_coordinates_dat) < 3) {
		colnames(genomic_coordinates_dat) <- c("chromosome", "pos")
		genomic_coordinates_dat$pos <- as.numeric(paste(genomic_coordinates_dat$pos))
	} else {
		colnames(genomic_coordinates_dat) <- c("chromosome", "start", "end")
		genomic_coordinates_dat$start <- as.numeric(paste(genomic_coordinates_dat$start))
		genomic_coordinates_dat$end <- as.numeric(paste(genomic_coordinates_dat$end))
	}
	rownames(genomic_coordinates_dat) <- seq(1:nrow(genomic_coordinates_dat))
	return(genomic_coordinates_dat)
}

to_genomic_regions <- function(regions_dat, chr_colname="chromosome", start_colname="start", end_colname="end", strand_colname="strand", labels=NULL, label_colname="gene", order_coordinates=FALSE, remove_duplicate_labels=FALSE, keep_chr_prefix=FALSE) {
    regions_dat[,chr_colname] <- paste0(regions_dat[,chr_colname])
    if (!keep_chr_prefix) { regions_dat[,chr_colname] <- gsub("chr", "", regions_dat[,chr_colname]) }
    regions_dat[,start_colname] <- as.numeric(paste0(regions_dat[,start_colname]))
    if (end_colname %in% colnames(regions_dat)) { regions_dat[,end_colname] <- as.numeric(paste0(regions_dat[,end_colname])) } else { end_colname <- start_colname }
    if (is.null(labels)) {
	    if (label_colname %in% colnames(regions_dat)) { labels <- paste0(regions_dat[,label_colname]) 
	    } else { labels <- genomic_coordinates_to_strings(regions_dat[,chr_colname], regions_dat[,start_colname], regions_dat[,end_colname]) }
    }
    
    if(order_coordinates) {
        badly_ordered <- regions_dat[,start_colname] > regions_dat[,end_colname]
        badly_ordered_starts <- regions_dat[badly_ordered, start_colname]
        regions_dat[badly_ordered, start_colname] <- regions_dat[badly_ordered, end_colname]
        regions_dat[badly_ordered, end_colname] <- badly_ordered_starts
        #ordered_coordinates <- t(apply(TS_regions[,c("start", "end")], 1, sort))
        #regions_dat[,start_colname] <- ordered_coordinates[,1]
        #regions_dat[,end_colname] <- ordered_coordinates[,2]
    }
    
    if (strand_colname %in% colnames(regions_dat)) { regions_vector <- makeGRangesFromDataFrame(regions_dat, seqnames.field=chr_colname,start.field=start_colname, end.field=end_colname, strand.field=strand_colname)
    } else { regions_vector <- makeGRangesFromDataFrame(regions_dat, seqnames.field=chr_colname,start.field=start_colname, end.field=end_colname) }
    names(regions_vector) <- labels
	
    if(remove_duplicate_labels) {
        duplicated_labels <- unique(labels[duplicated(labels)])
        duplicate_indices_to_remove <- c(unlist(sapply(duplicated_labels, function(duplicated_label) { relevant_indices <- which(labels== duplicated_label); return(relevant_indices[-c(which.max(width(regions_vector[relevant_indices])))]) })))
        regions_vector <- regions_vector[which(!(1:length(regions_vector) %in% duplicate_indices_to_remove))]
    }
    
	return(regions_vector)
}

# Returns the total footprint of the given GRanges object
footprint <- function(regions_granges) {
    if(is.null(regions_granges)) { return(0) }
    regions_granges <- intersect(regions_granges, regions_granges)
    return(as.numeric(sum(end(regions_granges) - start(regions_granges))))
}

# Returns boolean vector indicating which positions are indels and which are just snps, based on the ref and alt vectors.
is_indel <- function(ref, alt, verbose=FALSE) {
	return(nchar(paste0(ref)) > 1 | nchar(paste0(alt)) > 1)
}

# Gets and returns the requested entity from global context.
get_global <- function(entity, verbose=FALSE) {
    if(!is.character(entity)) { entity = deparse(substitute(entity)) }
    if(exists(entity, envir = globalenv())) { return(get(entity, envir=globalenv())) 
    } else { 
        if(verbose) { cat(paste0(entity," not found in global context!\n")) }
        return(NULL) 
    }
}
#haha <- function() { for(i in 1:10) { print(paste0(i,", ",get_global("i"),", ",get_global(i))) } }
#haha()

# Sets the given entity in global context, optionally with a different name.
set_global <- function(entity, name=NULL, verbose=FALSE) {
    if(is.null(name)) { name = deparse(substitute(entity)) }
    if(verbose) { cat(paste0(name, " saved in global context!\n")) }
    assign(name, entity, pos=globalenv())
    #return(0)
}

##########################################################################
## MULTI_PLOT - Plot multiple lines with automatic axis determination. ##
##########################################################################
# Plots multiple lines with automatic x-axis and y-axis determination.
# lines must be an environment where key is name of line and value is the line that should be plotted.
# Instead of lines, can specify individual lines_x and lines_y environments, where again key is name of line and x and y values split in parallel environments.
# cols and ltys must be lists where the entry names are the names of the lines and values are the values (of the color/line type). 
multi_plot <- function(lines=NULL, lines_x=NULL, lines_y=NULL, main, xlab="", ylab="", mtext=NULL, cols=list(), ltys=list(), exclude_from_legend=NULL, xmin=NULL, xmax=NULL, ymin=NULL, ymax=NULL, xaxs="r", yaxs="r", type="l", file=NULL, draw_legend=TRUE, legend_location="topright", legend_cex=1, buffer=0.1, force_xmin_0=TRUE, force_ymin_0=TRUE, verbose=FALSE, lwd=1, cex.lab=1, cex.axis=1, cex.main=1, mtext_cex=1, mar=NULL) {
	# Builds lines environment from lines_x and lines_y if lines is not specified as a parameter.
	if(is.null(lines)) {
		if(is.null(lines_x) || is.null(lines_y)) { return(0) }
		lines <- new.env()
		for(line_label in ls(lines_x)) {
			x <- lines_x[[line_label]]
			y <- lines_y[[line_label]]
			if(is.null(y) || length(y) > length(x)) { return(0) }
			if(length(y) < length(x)) { y <- c(y, rep(NA, length(x)-length(y))) }
			lines[[line_label]] <- data.frame(cbind(x, y), stringsAsFactors = FALSE)
		}
	}
	line_labels <- ls(lines)
	num_lines <- length(line_labels)

	# Parse vertical or horizontal abline from string
	parse_abline <- function(string) {
		parsed_string <- try(regmatches(string, regexec("^([v|h])=([0-9]+)$", string, perl=TRUE))[[1]], silent=TRUE)
		if(inherits(parsed_string, "try-error") || length(parsed_string) != 3 || sum(sapply(parsed_string, is.na)) != 0) { 
			return(NULL) 
		} else { return(parsed_string[2:3]) }
	}
	
	# Find optimal min and max x and y values for the plot, and store which indices have ablines
	horizontal_abline_indices <- new.env()
	vertical_abline_indices <- new.env()
	x_min <- NULL
	x_max <- NULL
	y_min <- NULL
	y_max <- NULL
	for(i in 1:num_lines) {
		line <- lines[[(line_labels[i])]]
		abline_def <- parse_abline(line)
		if(is.null(abline_def)) { # line is not an abline
			line_x <- line$x
			line_y <- line$y
			if(is.null(line_x) || is.null(line_y)) { return(0) }
			x_min <- min(c(x_min, line_x[!is.na(line_x)]))
			x_max <- max(c(x_max, line_x[!is.na(line_x)]))
			y_min <- min(c(y_min, line_y[!is.na(line_y)]))
			y_max <- max(c(y_max, line_y[!is.na(line_y)]))
		} else { # line is an abline
			abline_val = as.numeric(abline_def[2])
			if(abline_def[1] == "v") { # vertical line
				vertical_abline_indices[[paste(i)]] <- abline_val
				x_min <- min(c(x_min, abline_val))
				x_max <- max(c(x_max, abline_val))
			} else { # abline_def[1] == "h"; horizontal line
				horizontal_abline_indices[[paste(i)]] <- abline_val
				y_min <- min(c(y_min, abline_val))
				y_max <- max(c(y_max, abline_val))
			}
		}
	}
	if(!is.null(xmin)) { x_min <- xmin }
	if(!is.null(xmax)) { x_max <- xmax }
	if(!is.null(ymin)) { y_min <- ymin }
	if(!is.null(ymax)) { y_max <- ymax }
	x_min <- as.numeric(x_min)
	x_max <- as.numeric(x_max)
	y_min <- as.numeric(y_min)
	y_max <- as.numeric(y_max)
	x_axis_length <- x_max - x_min
	y_axis_length <- y_max - y_min
	if (x_min < 0 || force_xmin_0==FALSE) { x_min <- x_min-(buffer*x_axis_length) } else { x_min <- 0 }
	x_max <- x_max+(buffer*x_axis_length)
	if (y_min < 0 || force_ymin_0==FALSE) { y_min <- y_min-(buffer*y_axis_length) } else { y_min <- 0 }
	y_max <- y_max+(buffer*y_axis_length)

	# Sets colors and line types using the defined entries in cols and ltys (default col and lty for undefined entries are both 1).
	cols <- cols[line_labels] # order cols, generating NULL elements for missing col entries
	ltys <- ltys[line_labels] # order ltys, generating NULL elements for missing lty entries
	undefined_cols_indices <- sapply(seq(1:length(cols)), function(j) { return(is.null(cols[j][[1]])) }) # Fill in missing col entries with default col
	undefined_ltys_indices <- sapply(seq(1:length(ltys)), function(j) { return(is.null(ltys[j][[1]])) }) # Fill in missing lty entries with default lty 
	names(cols)[undefined_cols_indices] <- line_labels[undefined_cols_indices]
	names(ltys)[undefined_ltys_indices] <- line_labels[undefined_ltys_indices]
	cols[undefined_cols_indices] <- 1
	ltys[undefined_ltys_indices] <- 1
	col <- unlist(cols)
	lty <- unlist(ltys)
	
	# Opens pdf file connection if file parameter specified.
	if(!is.null(file)) { pdf(file) }
	
	# Plot each line
	if(!is.null(mar)) { default_mar = par("mar"); par(mar=mar) }
	plot(c(), c(), main=main, xlab=xlab, ylab=ylab, xlim=c(x_min,x_max), ylim=c(y_min,y_max), xaxs=xaxs, yaxs=yaxs, cex.lab=cex.lab, cex.axis=cex.axis, cex.main=cex.main)
	for(i in 1:num_lines) {
		if (verbose) { cat(paste0("Plotting ", line_labels[i], "...")) }
		if(!is.null(vertical_abline_indices[[paste(i)]])) { # vertical abline
			abline(v=vertical_abline_indices[[paste(i)]], col=col[i], lty=lty[i], lwd=lwd)
		} else if(!is.null(horizontal_abline_indices[[paste(i)]])) { # horizontal abline
			abline(h=horizontal_abline_indices[[paste(i)]], col=col[i], lty=lty[i], lwd=lwd)
		} else { # not abline
			lines(lines[[(line_labels[i])]], col=col[i], lty=lty[i], lwd=lwd, type=type)
		}
		if (verbose) { cat("Done.\n") }
	}
	if (verbose) { cat("Finished all plots.\n") }

	if(draw_legend) { 
		if(is.null(exclude_from_legend)) { legend(legend_location, legend=line_labels, col=col, lty=lty, lwd=lwd, cex=legend_cex) 
		} else { legend(legend_location, legend=line_labels[!(line_labels %in% ls(exclude_from_legend))], col=col[!(line_labels %in% ls(exclude_from_legend))], lty=lty[!(line_labels %in% ls(exclude_from_legend))], lwd=lwd, cex=legend_cex)  }
	}
	if(!is.null(mtext)) { mtext(mtext, cex=mtext_cex) }

	if(!is.null(file)) { dev.off() }
	if(!is.null(mar)) { par(mar=default_mar); dev.off() }
	return(1)
}
#lines <- new.env()
#lines[["CDH_sim_distrib"]] <- density(sapply(seq(1:nrow(sim_dat_CDH_TES)), function(i) { return(sum(sim_dat_CDH_TES[i,]>1)) } ))
#lines[["cases_line"]] <- paste0("v=", sum(TES_genes_dat$CDH_snv_count > 1))
#lines[["SSC_bootstrap"]] <- density(SSC_bootstrap)
#cols <- list("CDH_sim_distrib"=1, "cases_line"="red", "SSC_bootstrap"=3)
#ltys <- list("CDH_sim_distrib"=2, "SSC_bootstrap"=3)
#multi_plot(lines=lines, cols=cols, ltys=ltys, main="test", xlab="xlab", ylab="ylab")

##########################################################################
## olRanges - Identify Range Overlaps ##
##########################################################################
# Author: Thomas Girke
# Last update: 8-Feb-11
# Details on usage and use cases are available here:
# http://manuals.bioinformatics.ucr.edu/home/ht-seq#TOC-Analysis-Routines-with-IRanges-Geno
# Utility: identify overlaps in range data sets, such as annotation or alignment positions defined
# by two IRanges/GRanges objects.  
## Overlap types
	## olup: startup & endin
	## Q --------------
	## S       -------------

	## oldown: startin & enddown
	## Q       -------------
	## S --------------

	## inside: startin & endin 
	## Q     -----
	## S --------------

	## contained: startup & enddown
	## Q --------------
	## S     -----

###########################################################
## (A) olRanges Function for IRanges and GRanges Objects ##
###########################################################
olRanges <- function(query, subject, output="gr", ...) {
    require(GenomicRanges); require(IRanges)
        
    ## Input check
    if(!((class(query)=="GRanges" & class(subject)=="GRanges") | (class(query)=="IRanges" & class(subject)=="IRanges"))) {
        stop("Query and subject need to be of same class, either GRanges or IRanges!")
    }
        
    ## Find overlapping ranges
    if(class(query)=="GRanges") {
    	seqlengths(query) <- rep(NA, length(seqlengths(query)))
    	seqlengths(subject) <- rep(NA, length(seqlengths(subject)))
    }
	olindex <- as.matrix(findOverlaps(query, subject, ...))
    query <- query[olindex[,1]]
    subject <- subject[olindex[,2]]
    olma <- cbind(Qstart=start(query), Qend=end(query), Sstart=start(subject), Send=end(subject))
    
    ## Pre-queries for overlaps
	startup <- olma[,"Sstart"] < olma[,"Qstart"]
	enddown <- olma[,"Send"] > olma[,"Qend"]
	startin <- olma[,"Sstart"] >= olma[,"Qstart"] & olma[,"Sstart"] <= olma[,"Qend"]
	endin <- olma[,"Send"] >= olma[,"Qstart"] & olma[,"Send"] <=  olma[,"Qend"]

	## Overlap types
	olup <- startup & endin
	oldown <- startin & enddown
	inside <- startin & endin 
	contained <- startup & enddown
	
    ## Overlap types in one vector
	OLtype <- rep("", length(olma[,"Qstart"]))
	OLtype[olup] <- "olup"
	OLtype[oldown] <- "oldown"
	OLtype[inside] <- "inside" 
	OLtype[contained] <- "contained"
	
    ## Overlap positions
	OLstart <- rep(0, length(olma[,"Qstart"]))
	OLend <- rep(0, length(olma[,"Qstart"]))
	OLstart[olup] <- olma[,"Qstart"][olup]
	OLend[olup] <- olma[,"Send"][olup]
	OLstart[oldown] <- olma[,"Sstart"][oldown]
	OLend[oldown] <- olma[,"Qend"][oldown]
	OLstart[inside] <- olma[,"Sstart"][inside]
	OLend[inside] <- olma[,"Send"][inside]
	OLstart[contained] <- olma[,"Qstart"][contained]
	OLend[contained] <- olma[,"Qend"][contained]

	## Absolute and relative length of overlaps
	OLlength <- (OLend - OLstart) + 1
    OLpercQ <- OLlength/width(query)*100
    OLpercS <- OLlength/width(subject)*100
	
	## Output type
        oldf <- data.frame(Qindex=olindex[,1], Sindex=olindex[,2], olma, OLstart, OLend, OLlength, OLpercQ, OLpercS, OLtype)
        if(class(query) == "GRanges") {
		oldf <- cbind(space=as.character(seqnames(query)), oldf)
	}
        if(output=="df") {
                return(oldf)
        }
        if(output=="gr") {
                if(class(query)=="GRanges") {
                        elementMetadata(query) <- cbind(as.data.frame(elementMetadata(query)), oldf)
                }
                if(class(query)=="IRanges") {
                        query <- GRanges(seqnames = Rle(rep("dummy", length(query))), ranges = IRanges(start=oldf[,"Qstart"], end=oldf[,"Qend"]), strand = Rle(strand(rep("+", length(query)))), oldf)  
                }
                return(query)
        }
}

## Run olRanges function
## Sample Data Sets
# grq <- GRanges(seqnames = Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)), 
#                ranges = IRanges(seq(1, 100, by=10), end = seq(30, 120, by=10)), 
#                strand = Rle(strand(c("-", "+", "-")), c(1, 7, 2)))
# grs <- shift(grq[c(2,5,6)], 5)
# olRanges(query=grq, subject=grs, output="df") 
# olRanges(query=grq, subject=grs, output="gr") 


################################################
## (B) Older rangeOL Function for Data Frames ##
################################################
rangeOL <- function(featureDF=featureDF, label="ID1", start=25, end=35) {
	## Pre-queries for overlaps
	qstartup <- start < featureDF[,1]
	qenddown <- end > featureDF[,2]
	qstartin <- start >= featureDF[,1] & start <= featureDF[,2]
	qendin <- end >= featureDF[,1] & end <= featureDF[,2]

	## Overlap types
	qolup <- qstartup & qendin
	qoldown <- qstartin & qenddown
	qinside <- qstartin & qendin 
	qcontains <- qstartup & qenddown

	## Overlap types in one vector
	OLtype <- rep("", length(featureDF[,1]))
	OLtype[qolup] <- "olup"
	OLtype[qoldown] <- "oldown"
	OLtype[qinside] <- "inside" 
	OLtype[qcontains] <- "contained"

	## Overlap Positions
	OLstart <- rep("", length(featureDF[,1]))
	OLend <- rep("", length(featureDF[,1]))
	OLstart[qolup] <- featureDF[,1][qolup]
	OLend[qolup] <- end
	OLstart[qoldown] <- start
	OLend[qoldown] <- featureDF[,2][qoldown]
	OLstart[qinside] <- start
	OLend[qinside] <- end
	OLstart[qcontains] <- featureDF[,1][qcontains]
	OLend[qcontains] <- featureDF[,2][qcontains]

	## Length of overlaps
	OLlength <- (as.numeric(OLend) - as.numeric(OLstart)) + 1 
	
	## Feature label
	OLquery <- rep("", length(featureDF[,1]))
	OLquery[nchar(OLtype)>0] <- paste(label, " (", start, ":", end, ")", sep="")
	
	## Output
	featureDF <- cbind(featureDF, OLquery, OLtype, OLstart, OLend, OLlength)
	return(featureDF[nchar(OLtype)>0,])
}

## Run rangeOL function
# featureDF <- data.frame(Start=c(1, 20, 90, 150, 180), End=c(40, 30, 105, 160, 200), Feature="mRNA")
# rangeOL(featureDF=featureDF, label="ID1", start=25, end=35)
# lapply(seq(along=featureDF[,1]), function(x) rangeOL(featureDF=featureDF, label=featureDF[x,3], start=as.numeric(featureDF[x,1]), end=as.numeric(featureDF[x,2])))




#ensembl_names_dat <- read.csv("EnsemblID_genename.txt", sep="\t")
#ensembl_names <- new.env()
#for(i in 1:nrow(ensembl_names_dat)) {
#	ensembl_names[[paste(ensembl_names_dat$Gene.stable.ID[i])]] <- paste(ensembl_names_dat$Gene.name[i])
#}

# Get genes for the specified pathways
get_pathway_genes <- function(pathways=c("Wnt signaling pathway", "TGF-beta signaling pathway", "Notch signaling pathway", "mRNA surveillance pathway", "Ras signaling pathway")) {
    kegg_pathways <- keggList("pathway", "hsa") # get human pathway names and KEGG IDs
    kegg_pathways <-  kegg_pathways[kegg_pathways %in% paste0(pathways, " - Homo sapiens (human)")]
    pathway_genes <- data.frame(t(sapply(1:length(kegg_pathways), function(i) { 
        kegg_pathway_name <- gsub(" - Homo sapiens \\(human\\)", "", kegg_pathways[i])
        kegg_pathway_id <- names(kegg_pathways[i])
        kegg_pathway <- keggGet(kegg_pathway_id)[[1]]
        kegg_pathway_genes <- unlist(sapply(kegg_pathway$GENE, function(x) { gene <- unlist(strsplit(x, ";")); if(length(gene) > 1) { return(gene[1]) } else { return(NULL) } }))
        return(c(kegg_pathway_name, paste0(kegg_pathway_genes, collapse=",")))
    }))); colnames(pathway_genes) <- c("pathway", "genes")
    pathway_genes$pathway <- paste0(pathway_genes$pathway); pathway_genes$genes <- paste0(pathway_genes$genes)
    
    return(pathway_genes)
}
#chd_pathway_genes <- get_pathway_genes(c("Wnt signaling pathway", "TGF-beta signaling pathway", "Notch signaling pathway", "mRNA surveillance pathway", "Ras signaling pathway"))

# Hash chromosome lengths for all chromosomes.
hash_chromosome_lengths <- function(version="hg19") {
	chr_lengths <- sapply(c(1:22,"X","Y"), function(chr) length(get_refseq(chr, version=version)[[1]]))
	names(chr_lengths) <- c(1:22,"X","Y")
	set_global(chr_lengths)
	return(chr_lengths)
}
#chr_lengths <- hash_chromosome_lengths()

# Grab the specified n number of random variants from gnomAD, and returns it in the specified format: "data.frame", "GRanges", or "string"
get_random_gnomad_variants <- function(n=10000, format="data.frame") {
    chr_lengths <- get_global("chr_lengths"); if (is.null(chr_lengths)) { chr_lengths <- hash_chromosome_lengths() }
    chr_lengths <- chr_lengths[which(names(chr_lengths) != "Y")] # Remove chrY, since we do not have gnomAD data for it.
    print(chr_lengths)
    chromosomes <- names(chr_lengths)
    chr_proportions <- chr_lengths / sum(as.numeric(chr_lengths))
    variants_per_chrom <- ceiling(n * chr_proportions)
    GNOMAD_DIRECTORY = data_path("gnomAD")
    
    random_gnomad_variants <- unlist(sapply(1:length(chromosomes), function(i) { 
        print(paste0("Grabbing ",variants_per_chrom[i]," random variants for chromosome ",chromosomes[i],"..."))
        gnomad_file = paste0(GNOMAD_DIRECTORY,"/gnomad.genomes.r2.0.2.sites.chr",chromosomes[i],".vcf.bgz.2.bed.gz")
        if(file.exists(gnomad_file)) { return(system(paste0("zcat ",gnomad_file," | shuf -n ",variants_per_chrom[i]), intern=TRUE))#, function(variant) unlist(strsplit(variant, "\t"))#[c(1,2,4,5)])
        } else { print(paste0("Skipping chromosome ",chromosomes[i],": could not find gnomAD file ",gnomad_file)); return(c()) }
    }))
    total_variants_grabbed = length(random_gnomad_variants)
    random_gnomad_variants <- data.frame(t(sapply(1:total_variants_grabbed, function(i) {
        print(paste0("Parsing variants...[",i," / ",total_variants_grabbed,"]"))
        return(unlist(strsplit(random_gnomad_variants[i], "\t"))[c(1,2,4,5)])
    })))
    colnames(random_gnomad_variants) <- c("Chrom", "Position", "Ref", "Alt")
    random_gnomad_variants <- unfactorize(random_gnomad_variants[sample(1:nrow(random_gnomad_variants),n),])
    
    if(format == "GRanges") { random_gnomad_variants  <- to_genomic_regions(random_gnomad_variants, chr_colname="Chrom", start_colname="Position") 
    } else if (format == "string") { random_gnomad_variants  <- paste0(random_gnomad_variants$Chrom,":",random_gnomad_variants$Position,random_gnomad_variants$Ref,">",random_gnomad_variants$Alt) }
    
    print(paste0("Returning ",n," total random gnomAD variants!"))
    
    return(random_gnomad_variants)
}

generate_random_variants <- function(n, within_granges=NULL, version="hg19") {
    if(!("chr_lengths" %in% ls())) { chr_lengths <<- hash_chromosome_lengths() }
    chromosomes <- sort(sample(c(1:22,"X","Y"), n, replace=TRUE, prob=chr_lengths))
    chromosome_counts <- table(chromosomes)
    atcg <- c("A","T","C","G")
    variants_per_chromosome <- sapply(names(chromosome_counts), function(chr) { 
        cat(paste0("Generating random variants for chromosome chr", chr, "...\n"))
        chr_refseq <- get_refseq(chr, version=version)[[1]]
        if(is.null(within_granges)) { eligible_positions <- 1:length(chr_refseq)
        } else { cat("    Limiting to eligible positions..."); within_granges_chr <- within_granges[seqnames(within_granges) == chr]; eligible_positions <- c(unlist(sapply(1:length(within_granges_chr), function(i) { return(start(within_granges_chr[i]):end(within_granges_chr[i])) }))); cat("Done.\n") }
        chr_positions <- sample(eligible_positions, chromosome_counts[[chr]]) 
        chr_refs <- sapply(chr_positions, function(pos) return(paste0(chr_refseq[pos])))
        is_N <- !(chr_refs %in% atcg)
        while(sum(is_N) > 0) {
            chr_positions[is_N] <- sample(eligible_positions, sum(is_N))
            chr_refs[is_N] <- sapply(chr_positions[is_N], function(pos) return(paste0(chr_refseq[pos])))
            is_N[is_N] <- !(chr_refs[is_N] %in% atcg)
        }
        chr_alts <- sapply(chr_refs, function(chr_ref) return(sample(atcg[atcg != chr_ref], 1)))
        cat("    Done.\n")
        return(cbind(chr, chr_positions, chr_refs, chr_alts))
    })
    variants <- data.frame(rbindlist(lapply(variants_per_chromosome, function(x) return(data.frame(x)))))
    variants <- cbind(variants, rep("random", nrow(variants)))
    colnames(variants) <- c("Chrom", "Position", "Ref", "Alt", "sample")
    variants[,1] <- paste0(variants[,1]); variants[,2] <- as.numeric(paste0(variants[,2])); variants[,3] <- paste0(variants[,3]); variants[,4] <- paste0(variants[,4]); variants[,5] <- paste0(variants[,5]);
    return(variants)
}

# Generate random region in the specific chromosome (can be given either with or without prefix 'chr'), with the specified length.
gen_random_region <- function(chromosome, length) {
	chr = gsub("chr", "", paste(chromosome))
	start = floor(runif(1, min=1, max=(chr_lengths[[chr]]-length+2)))
	#start=(chr_lengths[[chr]]-length+1)
	end = start+length-1
	return(c(start, end))
}

## Generate random controls with same sample distribution (and count) as cases.
#random_TSS_length_regions <- c()
#for(i in 1:nrow(CDH_HDE_TSS_20k_snps)) {
#	
#}


# Turn data frame with character strings or numerics into data frame with factors in place of all character strings or numerics.
factorize <- function(df){
    return(data.frame(as.matrix(df)))
}

# Turn data frame with factors into data frame with character strings or numerics in place of all factors.
unfactorize <- function(df, numeric=TRUE){
    sapply_out <- sapply(which(sapply(df, class) == "factor"), function(i) { 
        df_i_char <- as.character(df[[i]])
        if(numeric) { 
            df_i_num <- as.numeric(df_i_char)
            if(sum(is.na(df_i_num)) == 0) { df[[i]] <<- df_i_num; return() } 
        }
        df[[i]] <<- df_i_char  
    })
    return(df)
}
#a <- data.frame(cbind(rbind(c(1,2,3),rbind(paste0(c(1,2,3)))))); a <- cbind(a, a[,1]); a$X2 <- paste0(a$X2); a$X3 <- as.numeric(paste0(a$X3)); colnames(a)[ncol(a)] <- "X4"; a <- cbind(a, paste0("chr",6:7)); colnames(a)[ncol(a)] <- "X5"
#a; a$X1+10; a$X2; a$X3+10; a$X4+10; a$X5
#a <- unfactorize(a)
#a; a$X1+10; a$X2; a$X3+10; a$X4+10; a$X5
#a <- factorize(a)
#a; a$X1+10; a$X2; a$X3+10; a$X4+10; a$X5
#a <- unfactorize(factorize(a))
#a; a$X1+10; a$X2; a$X3+10; a$X4+10; a$X5

# Hash codon to amino acid mappings
hash_codon_amino_acid_mappings <- function() {
    aa_codon_dat <- read.csv(data_path("amino_acid_codon_mappings.txt"), sep="\t")
    codon_aa_mappings <- new.env()
    sapply_out <- sapply(1:nrow(aa_codon_dat), function(i) { for(codon in unlist(strsplit(paste0(aa_codon_dat$dna_codons[i]), ","))) { codon_aa_mappings[[codon]] <- paste0(aa_codon_dat$amino_acid_code[i]) } })
    return(codon_aa_mappings)
}
#codon_aa_mappings <- hash_codon_amino_acid_mappings()

# Returns whether or not the mutation is synonymous, for each trimer mutation in trimer_mutation_strings parameter vector.
is_synonymous <- function(trimer_mutation_strings) {
    return(sapply(trimer_mutation_strings, function(trimer_mutation_string) {
        trimer_mutation_string_split <- unlist(strsplit(paste0(trimer_mutation_string), "->"))
        return(codon_aa_mappings[[trimer_mutation_string_split[1]]] == codon_aa_mappings[[trimer_mutation_string_split[2]]])
        }))
}

# Hash trimer mutation rates
hash_trimer_mut_rates <- function() {
	expected_mut_counts <- c()
	tri_mutation_rates_dat <- read.csv(data_path("3mer_table_from_felix.txt"), sep=" ")
	tri_mutation_rates_dat$mu_snp <- as.numeric(tri_mutation_rates_dat$mu_snp)
	tri_mutation_rates <- new.env()
	#sapply_out <- sapply(unique(tri_mutation_rates_dat$from), function(x) { tri_mutation_rates[[paste(x)]] <- 0 })
	#sapply_out <- sapply(seq(1:nrow(tri_mutation_rates_dat)), function(i) { tri_mutation_rates[[paste(tri_mutation_rates_dat$from[i])]] <- tri_mutation_rates[[paste(tri_mutation_rates_dat$from[i])]] + tri_mutation_rates_dat$mu_snp[i] })
	
	sapply_out <- sapply(1:nrow(tri_mutation_rates_dat), function(i) { tri_mutation_rates[[paste0(tri_mutation_rates_dat$from[i], "->", tri_mutation_rates_dat$to[i])]] <- tri_mutation_rates_dat$mu_snp[i] })
	
	#tri_mutation_rates 
	return(tri_mutation_rates)
}
#tri_mutation_rates <- hash_trimer_mut_rates()

# Returns DNA methylation Granges object loaded from BigWig, for the specified eID (data from Roadmap).
load_dna_methylation <- function(eid="E008", include_all_positions=FALSE) {
    filename = data_path(paste0(c("fractional_methylation", "bigwig", paste0(eid,"_WGBS_FractionalMethylation.bigwig")), collapse=FILEPATH_DELIM))
    if (!file.exists(filename)) { return(NULL) }
    meth_granges <- import(filename)
    if(include_all_positions) {
        non_gc_meth_granges <- invert_granges(meth_granges)
        non_gc_meth_granges$score <- 0
        meth_granges <- c(meth_granges, non_gc_meth_granges)
        meth_granges <- sortSeqlevels(meth_granges)
        meth_granges <- sort(meth_granges)
    }
    return(meth_granges)
}
#meth_granges <- load_dna_methylation("E008") # h9 stem cells
#meth_granges_all <- load_dna_methylation("E008", include_all_positions=TRUE) #h9 stem cells
#meth_granges <- load_dna_methylation("E017") #imr90
#meth_granges_all <- load_dna_methylation("E017", include_all_positions=TRUE) #imr90
if (FALSE) {
    cpg_coverage <- c()
    avg_methylation <- c()
    methylation_qc <- data.frame()
    for(eid in get_relevant_roadmap_eids("CDH")) {
        meth_granges <- load_dna_methylation(eid)
        if(!is.null(meth_granges)) { 
            print(paste0(eid, " added."))
            cpg_coverage <- c(cpg_coverage, length(meth_granges)) 
            avg_methylation <- c(avg_methylation, mean(score(meth_granges)))
            methylation_qc <- rbind(methylation_qc, cbind(eid, length(meth_granges), mean(score(meth_granges))))
        }
    }
    colnames(methylation_qc) <- c("eid", "cpg_coverage", "mean_cpg_meth")
    hist(cpg_coverage)
    hist(avg_methylation)
    plot(cpg_coverage, avg_methylation)
    text(methylation_qc$cpg_coverage, methylation_qc$mean_cpg_meth, methylation_qc$eid)
}

# Gets DNA methylation fractional score for the specified genomic positions, padded using the specified bw bandwidth around start if end is not specfied.
get_dna_methylation <- function(chromosome, start=NULL, end=NULL, bw=0) {
    chromosome <- paste0(chromosome)
    
    if (chromosome=="all") { return(score(meth_granges))
    } else if (is.null(start)) { 
        start = 1
        end = seqlengths(Hsapiens)[paste0("chr", gsub("chr", "", chromosome))]
        dna_meth <- score(meth_granges[seqnames(meth_granges) == paste0("chr", gsub("chr", "", chromosome))]) #detect_index(start(meth_granges_all), function(x) x >= 10470)
    } else if (is.null(end)) { start = (start-bw/2); end = (start+bw/2) } # dna_meth <- score(meth_granges[seqnames(meth_granges) == paste0("chr", gsub("chr", "", chromosome)) & start(meth_granges) <= (start-bw/2) & start(meth_granges) >= (start+bw/2)])
    
    if(length(start)>1) {
        query_dat <- data.frame(rep(paste0("chr", gsub("chr", "", chromosome)), length(start)), start)
    } else { query_dat <- data.frame(rep(paste0("chr", gsub("chr", "", chromosome)), end-start+1), seq(start, end)) }
    colnames(query_dat) <- c("chromosome", "start")
    query_granges <- to_genomic_regions(query_dat, keep_chr_prefix=TRUE)
    hits <- findOverlaps(query_granges, meth_granges_all)
    dna_meth <- rep(0, nrow(query_dat))
    dna_meth[queryHits(hits)] <- score(meth_granges_all[subjectHits(hits)])
    return(dna_meth)
}
#get_dna_methylation(1, 6174581, 6174590)
#get_dna_methylation("1", 6174581, bw=1000)
#get_dna_methylation("1", 11360767, end=11362767)

# Does the inverse of the R aggregate function, stretching the aggregated data frame dat into a larger data frame again by splitting the columns identified in the "by" parameter (with the "delim" parameter as delimiter).
unaggregate <- function(dat, by, delim=",", add_index=TRUE, verbose=FALSE) {
    if (verbose) { print(paste0("by: ", paste(by,collapse=","))) }
    by_env <- new.env()
    num_elements <- c()
    for(b in by) { by_env[[paste0(b)]] <- strsplit(paste0(dat[,paste0(b)]), delim); num_elements <- c(num_elements, length(by_env[[paste0(b)]])) }; names(num_elements) <- by
    if(length(unique(num_elements)) > 1) { print("ERROR: Desired split element lengths are not all identical:"); print(num_elements); return() }
    
    not_by <- colnames(dat)[!(colnames(dat) %in% by)]
    if (verbose) { print(paste0("not_by: ", paste(not_by,collapse=", "))) }
    new_dat <- data.frame()
    for(i in 1:nrow(dat)) { 
        if (verbose) { print(paste0(i, " / ", nrow(dat))) }
        new_addition <- suppressWarnings(cbind(dat[i,not_by], rbind(sapply(by, function(b) { unlist(by_env[[paste0(b)]][i]) }))))
        if (add_index) { index <- 1:nrow(new_addition); new_addition <- cbind(new_addition, index) }
        new_dat <- rbind(new_dat, new_addition) 
    }  
    
    column_order <- colnames(dat); if (add_index) { column_order <- c(column_order, "index") }
    new_dat <- new_dat[,column_order]
    return(new_dat)
}

# Returns exon regions in a combined data frame for the specified genes. The transcripts parameter for which transcipt to pick can be "first", "longest" (default), "all", or the specific NCBI protein/transcript IDs.
get_exon_regions <- function(genes=NULL, transcripts="longest", separate_exon_entries=TRUE) {
    exons_dat <- read.table(data_path("refseq_exons.txt"), sep="\t", header=TRUE)
    
    if (!("length" %in% colnames(exons_dat))) {
        exon_starts <- strsplit(paste0(exons_dat$exonStarts), ",")
        exon_ends <- strsplit(paste0(exons_dat$exonEnds), ",")
        length <- sapply(1:nrow(exons_dat), function(i) { return(sum(abs(as.numeric(unlist(exon_ends[i])) - as.numeric(unlist(exon_starts[i])))+1)) } )
        exons_dat <- cbind(exons_dat, length)
        write.table(exons_dat, file=data_path("refseq_exons.txt"), sep="\t", row.names=FALSE, quote=FALSE)
    }
    if (!is.null(genes)) { exons_dat <- exons_dat[exons_dat$gene %in% genes,] }
    genes <- sort(unique(exons_dat$gene))
    ids <- sort(unique(exons_dat$id))
    
    if (transcripts[1] == "first") { 
        indices_to_keep <- sapply(genes, function(gene) { return(which(exons_dat$gene == gene)[1]) })
    } else if (transcripts[1] == "longest") {
        indices_to_keep <- sapply(genes, function(gene) { return(which(exons_dat$gene == gene)[which.max(exons_dat[exons_dat$gene == gene, "length"])]) })
    } else if (transcripts[1] == "all") { 
        indices_to_keep <- 1:nrow(exons_dat)
    } else if (sum(transcripts %in% ids) > 0) {
        indices_to_keep <- which(exons_dat$id %in% transcripts)
    } else { cat(paste0("ERROR: Invalid transcripts parameter \"", paste0(transcript,collapse=",") ,"\".\nUse \"first\", \"longest\", \"all\", or a vector of valid NCBI protein/transcript IDs.")); return() }
    
    exons_dat <- exons_dat[indices_to_keep,]
    
    if (separate_exon_entries) {
        exons_dat <- unaggregate(exons_dat, by=c("exonStarts", "exonEnds", "exonFrames"))[,c("chrom", "exonStarts", "exonEnds", "gene", "id", "index", "exonCount", "txStart", "txEnd", "length", "strand", "exonFrames")]
        colnames(exons_dat) <- c("chromosome", "start", "end", "gene", "transcript_id", "exon_num", "exon_count", "transcript_start", "transcript_end", "transcript_length", "strand", "frame")
        exons_dat$start <- as.numeric(paste0(exons_dat$start)); exons_dat$end <- as.numeric(paste0(exons_dat$end)); exons_dat$transcript_start <- as.numeric(paste0(exons_dat$transcript_start)); exons_dat$transcript_end<- as.numeric(paste0(exons_dat$transcript_end)); exons_dat$transcript_length<- as.numeric(paste0(exons_dat$transcript_length)); 
    } else {
        exons_dat <- exons_dat[, c("chrom", "exonStarts", "exonEnds", "gene", "id", "exonCount", "txStart", "txEnd", "length", "strand", "exonFrames")]
        colnames(exons_dat) <- c("chromosome", "exon_starts", "exon_ends", "gene", "transcript_id", "exon_count", "transcript_start", "transcript_end", "transcript_length", "strand", "exon_frames")
    }
    rownames(exons_dat) <- c()
    
    return(exons_dat)
}

# Returns expected mutation counts per individual for each region. The regions data frame must have columns "chromosome", "start", and "end".
get_expected_mut_counts <- function(regions, precision=3) {
	REFERENCE_SEQ_FOLDER <<- data_path("hg19_chr_fa")
	expected_mut_dat <- data.frame()
	regions$chromosome <- gsub("chr", "", regions$chromosome)
	for(chromosome in unique(regions$chromosome)) {
		#if (chromosome != "18") { next }
		cat(paste0("Reading chr", chromosome, "..."))
		refseq <- read.fasta(full_path(REFERENCE_SEQ_FOLDER, paste0("chr", chromosome, ".fa")), seqtype="DNA")
		cat("Done.\nProcessing regions...\n")
		regs <- regions[regions$chromosome==chromosome,]
		regs_length = nrow(regs)
		average_dna_methylation <- get_dna_methylation("all") # genome-wide average DNA methylation score
	    # average_dna_methylation <- get_dna_methylation(chromosome, 1, length(refseq[[1]])) # chromosome-wide average DNA methylation score
		for(reg_index in 1:nrow(regs)) {
		    print(paste0(reg_index, " / ", regs_length))
			start = regs$start[reg_index]
			end = regs$end[reg_index]
			average_dna_methylation <- get_dna_methylation(chromosome, start, end)
			seq_positions <- (start-1):(end+1)
			seq <- toupper(refseq[[1]][seq_positions])
			tris <- new.env()
			apply_out <- apply(permutations(n=4, r=3, v=c("A", "T", "C", "G"), repeats.allowed=TRUE), 1, function(x) { tris[[paste(x, collapse="")]] <- 0 })
			# previous_base = seq[1]
			# curr_base = seq[2]
			# next_base = seq[3]
			# for(i in 2:(length(seq)-1)) {
			# 	if(!is.na(previous_base) && !is.na(curr_base) && !is.na(next_base)) {
			# 		tri = paste(c(previous_base, curr_base, next_base), collapse="")
			# 		tris[[tri]] <- tris[[tri]] + 1
			# 	}
			# 	previous_base = curr_base
			# 	curr_base = next_base
			# 	next_base = seq[i+2]
			# }

			rollapply_out <- rollapply(1:(length(seq_positions)), width=3, by=precision, FUN = function(x) { 
			    print(x); tri = paste(seq[x], collapse=""); 
			    norm_dna_methylation <- 1 #get_dna_methylation(chromosome, seq_positions[x[2]], bw=100)/average_dna_methylation;
			    if(is.nan(norm_dna_methylation)) { norm_dna_methylation = 1 }; tris[[tri]] <- tris[[tri]] + norm_dna_methylation }
			    , align = "left")
			expected_mut_count = 0
			for(tri in ls(tris)) {
				if (length(tris[[tri]]) == 0 || length(tri_mutation_rates[[tri]]) == 0) { print("SKIPPING!"); next }
				#print(tris[[tri]])
				#print(tris[[tri]]*tri_mutation_rates[[tri]]*2)
				expected_mut_count = expected_mut_count + tris[[tri]]*tri_mutation_rates[[tri]]*2 # since 2 chromosomes.
			}
			expected_mut_count = expected_mut_count * precision
			
			expected_mut_dat <- rbind(expected_mut_dat, cbind(regs[reg_index,], expected_mut_count))
		}
		cat("Done.\n")
	}
	colnames(expected_mut_dat) <- c(colnames(regions), "exp_mutations")
	return(expected_mut_dat)
}
if (FALSE) {
    num_genes_to_test = 25
    random_start = 600
    genes <- TSS_regions_exp$gene # unique(unlist(strsplit(paste0(unique(ssc_combined_muts_TES$TS_gene)[random_start:(random_start+num_genes_to_test-1)]), ",")))[1:num_genes_to_test]
    empirical_counts_table <- table(unlist(strsplit(paste0(ssc_all_muts_TSS$TS_gene), ",")))
    empirical_counts <- sapply(1:length(genes), function(i) { print(i); gene <- genes[i]; empirical_counts_table[gene] }); empirical_counts[is.na(empirical_counts)] <- 0
    cdh_empirical_counts_table <- table(unlist(strsplit(paste0(cdh_combined_muts_TSS$TS_gene), ",")))
    cdh_empirical_counts <- sapply(1:length(genes), function(i) { print(i); gene <- genes[i]; cdh_empirical_counts_table[gene] }); cdh_empirical_counts[is.na(cdh_empirical_counts)] <- 0
    expected_counts <- TSS_regions_exp$exp_mutations * SSC_ALL_SAMPLE_COUNT; names(expected_counts) <- TSS_regions_exp$gene # get_expected_mut_counts(TES_regions[TES_regions$gene %in% genes,], precision=480)$exp_mutations * SSC_COMBINED_SAMPLE_COUNT; names(expected_counts) <- TES_regions$gene[TES_regions$gene %in% genes]
    #TSS_sim <- simulate_mutations(TSS_regions_exp, num_individuals=SSC_COMBINED_SAMPLE_COUNT, num_simulations=1000)
    TSS_sim_dat_results <- table(TSS_sim)
    ncells <- nrow(TSS_sim)*ncol(TSS_sim)
    cat(paste0("Num simulations: ", 1000, "\n"))
    cat(paste0("Mean SNV rate (SNVs/gene): ", sum(as.numeric(names(TSS_sim_dat_results))*TSS_sim_dat_results)/ncells, "\n"))
    cat(paste0("Num regions: ", nrow(regions), "\n"))
    sapply_out <- sapply(seq(1:length(TSS_sim_dat_results)), function(i) { cat(paste0("E(regions with ", names(TSS_sim_dat_results)[i], " SNVs): ", TSS_sim_dat_results[i]/1000, "\n")) })
    
    empirical_counts <- empirical_counts[order(names(empirical_counts))]
    cdh_empirical_counts <- cdh_empirical_counts[order(names(cdh_empirical_counts))]
    expected_counts <- expected_counts[order(names(expected_counts))]
    
    plot(density(TSS_sim), col="black", main="", xlim=c(0,6))
    lines(density(cdh_empirical_counts), col="red")
    lines(density(empirical_counts), col="blue")
    abline(v=mean(TSS_sim), col="black", lty=2)
    abline(v=mean(cdh_empirical_counts), col="red", lty=2) #*SSC_COMBINED_SAMPLE_COUNT/CDH_COMBINED_SAMPLE_COUNT
    abline(v=mean(empirical_counts), col="blue", lty=2)
    legend("topright", legend=c(paste0("SSC (N=",SSC_COMBINED_SAMPLE_COUNT,")"), paste0("CDH (N=",CDH_COMBINED_SAMPLE_COUNT,")"), paste0("Expected (N=",SSC_COMBINED_SAMPLE_COUNT,")")), col=c("red", "blue", "black"), lty=c(1,1,1))
    
    pdf(output_path("recurrent_TSS_hits_vs_expectation.pdf"))
    barplot(table(cdh_empirical_counts)[1:7], main=paste0("Recurrent TSS Hits vs. Expectation"), xlab="# mutations", ylab="# genes", cex.main=1.5, cex.axis=1.3, cex.lab=1.4, lwd=2, xlim=c(0,8), col=adjustcolor("red",alpha.f=0.2))
    barplot(table(empirical_counts)[1:7], col=adjustcolor("blue",alpha.f=0.2), xaxt="n", yaxt="n", add=TRUE)
    barplot(round(TSS_sim_dat_results/1000)[1:7], col=adjustcolor("black",alpha.f=0.2), xaxt="n", yaxt="n", add=TRUE)
    legend("topright", legend=c(paste0("CDH (N=",CDH_COMBINED_SAMPLE_COUNT,")"), paste0("SSC (N=",SSC_COMBINED_SAMPLE_COUNT,")"), paste0("Expected (N=",SSC_COMBINED_SAMPLE_COUNT,")")), col=c(adjustcolor("red",alpha.f=0.5), adjustcolor("blue",alpha.f=0.5), adjustcolor("black",alpha.f=0.5)), pch=c(15, 15, 15))
    mtext(paste0(length(genes), " genes in total"))
    dev.off()
    
    mean(empirical_counts)
    mean(cdh_empirical_counts)*SSC_COMBINED_SAMPLE_COUNT/CDH_COMBINED_SAMPLE_COUNT
    mean(TSS_sim)
    
    print_results <- function() {
        cat(paste0("Control sample count: ", SSC_COMBINED_SAMPLE_COUNT, "    (CDH sample count: ", CDH_COMBINED_SAMPLE_COUNT, ")"), "\n")
        cat(paste0("Expected hits/TSS using 1.2e-8 mutations/site/generation genome-wide mutation rate: ", (1.2e-8)*2*40000*SSC_COMBINED_SAMPLE_COUNT), "\n") 
        cat(paste0("Empirical hits/TSS (SSC): ", mean(empirical_counts)), "\n")
        cat(paste0("Empirical hits/TSS (CDH, normalized sample count): ", mean(cdh_empirical_counts)*SSC_COMBINED_SAMPLE_COUNT/CDH_COMBINED_SAMPLE_COUNT), "\n")
        cat(paste0("Expected hits/TSS using trimer mutation rate model: ", mean(expected_counts)), "\n")
    }
    print_results()
    
    boxplot(expected_counts ~ empirical_counts, main="Trimer Mutation Rate Model Results", xlab="Empirical Count", ylab="Expected Count")
    mtext(paste0("Expected hits using 1.1e-8 mutations/site/generation genome-wide mutation rate: ", (1.1e-8)*2*40000*SSC_COMBINED_SAMPLE_COUNT))
    
    # Print trimer mutation rates
    sapply(ls(tri_mutation_rates), function(x) { print(paste0(x, "  ", round(tri_mutation_rates[[x]], 12))) })
    
    dnmr <- read.table(data_path("DNMR/DNMR-DM.txt"), header=TRUE)
    dnmr_avg <- read.table(data_path("DNMR/DNMR-average.txt"), header=TRUE)
    
    hongjian_background <- read.table(data_path("MutationRate_20170710_rate.txt"), comment.char = "", header=TRUE)
    dnmr <- cbind(hongjian_background[,c(3:8)], sapply(1:nrow(hongjian_background), function(i) { sum(hongjian_background[i,c("p_synonymous","p_misense","p_nonsense")]) } ))
    colnames(dnmr) <- c("transcript", "exon_count", "chromosome", "transcript_start", "transcript_end", "gene", "snp_rate")
    dnmr$snp_rate <- as.numeric(paste0(dnmr$snp_rate)); dnmr$exon_count <- as.numeric(paste0(dnmr$exon_count)); dnmr$transcript_start <- as.numeric(paste0(dnmr$transcript_start)); dnmr$transcript_end <- as.numeric(paste0(dnmr$transcript_end))
    
    num_genes_to_test = 200
    random_start = 300
    genes <- unique(unlist(strsplit(sample(paste0(unique(ssc_combined_muts_TSS$TS_gene))[random_start:(random_start+num_genes_to_test-1)]), ",")))[1:num_genes_to_test]
    dnmr[dnmr$gene %in% genes,]
    regions <- get_exon_regions(genes, transcripts="longest") # genebody[genebody$gene %in% genes,1:4]; colnames(regions) <- c("chromosome", "start", "end", "gene"); 
    regions[,c("start","end")] <- t(sapply(1:nrow(regions), function(i) { sort(c(regions$start[i], regions$end[i])) }))
    expected <- get_expected_mut_counts(regions, precision=30)
    expected_colnames_to_keep <- c("chromosome","transcript_start","transcript_end","gene","transcript_id","exon_count","transcript_length","strand","exp_mutations")
    expected <- aggregate(exp_mutations ~ ., data=expected[,expected_colnames_to_keep], FUN=sum) #aggregate(expected[,c(expected_constant_colnames,"exp_mutations")], by=sapply(expected_constant_colnames, function(colname) { list(expected[,colname]) }), FUN=sum) 
    expected
    dnmr[dnmr$gene %in% genes,]
    dnmr_avg[dnmr_avg$gene %in% genes,]
    a <- merge(expected, dnmr, by="gene"); colnames(a)[which(colnames(a) %in% c("exp_mutations","snp_rate"))] <- c("model_rate", "hongjian_rate"); a <- a[a$model_rate != Inf,] # a <- cbind(a, a$end - a$start); a$exp_mutations <- a$exp_mutations/a[,ncol(a)]; colnames(a) <- c("gene", "chromosome", "start", "end", "model_rate", "mirDNMR_rate", "length")
    mean(a$model_rate / a$hongjian_rate)
    cor(a$model_rate, a$hongjian_rate, method="spearman")
    pdf(file=output_path("background_model_vs_mirdnmr.pdf"))
    scatterplot(a$model_rate, a$hongjian_rate, main="Background Mutation Model vs. Hongjian Rates", xlab="model rate", ylab="mirDNMR rate", xlim=c(0,6e-4), ylim=c(0,6e-4), cex.main=1.5, cex.axis=1.4, cex.lab=1.4, spread=FALSE)
    mtext("200 random genes (675.2 Kbp, exonic); model samples every 30th trimer", padj=-1.05, cex=1)
    legend("topright", legend=c("Smoothed Line", "Regression Line", paste0("Spearman = ", round(cor(a$model_rate, a$hongjian_rate, method="spearman"),3))), col=c("red","green","white"), lty=c(1,1), bty="n")
    dev.off()
    write.csv(a, file=output_path("gene_mutation_rates.csv"), row.names=FALSE, quote=FALSE)
    
    dnmr_annot <- merge(genebody, dnmr_avg); dnmr_annot <- cbind(dnmr_annot, abs(dnmr_annot$tes - dnmr_annot$tss)); colnames(dnmr_annot)[ncol(dnmr_annot)] <- "length"

    
    a[1:10,]
    b <- get_expected_mut_counts(regions[regions$gene %in% unique(regions$gene)[1:10],], precision=240)
    b_colnames_to_keep <- c("chromosome","transcript_start","transcript_end","gene","transcript_id","exon_count","transcript_length","strand","exp_mutations")
    b <- aggregate(exp_mutations ~ ., data=b[,b_colnames_to_keep], FUN=sum) #aggregate(expected[,c(expected_constant_colnames,"exp_mutations")], by=sapply(expected_constant_colnames, function(colname) { list(expected[,colname]) }), FUN=sum) 
    b
}

# Normalizes tensor or vector
norm_tensor <- function(x) { 
    return((x - mean(x)) / sd(x)) 
}

# Returns Refseq fasta object for the specified chromosome and assemby version.
get_refseq  <- function(chromosome, version="hg19", allow_BSgenome=TRUE) {
    if(allow_BSgenome) {
        if(version=="hg19") { return(list(BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg19::Hsapiens, paste0("chr",gsub("chr","",chromosome)))))
        } else if(version=="hg38") { return(list(BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38::Hsapiens, paste0("chr",gsub("chr","",chromosome))))) }
        #return(list(chr=strsplit(paste0(BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg19::Hsapiens, paste0("chr",gsub("chr","",chromosome)))),"")))
    } else { return(read.fasta(full_path(gsub("hg[0-9]+", version, REFERENCE_SEQ_FOLDER), paste0("chr", gsub("chr","",chromosome), ".fa")), seqtype="DNA")) }
}
#refseq1 <- get_refseq(12,  allow_BSgenome=FALSE)
#refseq2 <- get_refseq(12)

# Gets Refseq trimer strings for the union of regions specified in the GRanges object g.
get_trimer_counts_for_granges <- function(g_all, split_by_gc_meth=TRUE, work_folder=output_path("trimer_counts_for_granges_temp"), restart=FALSE, bin_size=1, version="hg19") {
    dir.create(file.path(work_folder), showWarnings = FALSE)
    bin_files <- c()
    tri_counts_names <- apply(permutations(n=5, r=3, v=c("A", "T", "C", "G", "N"), repeats.allowed=TRUE), 1, function(x) { paste(x, collapse="") })
    if(split_by_gc_meth) { tri_counts_names <- c(tri_counts_names, paste0(c("A", "T", "C", "G", "N"), "CmG"), paste0("CmG", c("A", "T", "C", "G", "N"))) } # apply(permutations(n=5, r=2, v=c("Cm", "T", "C", "G", "N")), 1, function(x) { paste0(x[1], c("A", "T", "C", "G", "Cm", "Gm", "N"), x[2]) })) }
    tri_counts_empty <- rep(0, length(tri_counts_names))
    names(tri_counts_empty) <- tri_counts_names
    cpg_trimer_names <- tri_counts_names[grepl("CG", tri_counts_names)]
    for(chromosome in unique(seqnames(g_all))) {
        if (!(chromosome %in% c(1:22,"X","Y"))) { next }
        cat(paste0("Loading refseq for chromosome ", chromosome, "..."))
        refseq <- get_refseq(chromosome, version)
        #if(split_by_gc_meth) { gc_meths <- get_dna_methylation(chromosome) }
        cat("Done.\n")
        g_chr <- g_all[seqnames(g_all) == chromosome]
        num_bins = ceiling(length(g_chr)/bin_size)
        for(bin in 1:num_bins) {
            bin_start = ((bin-1)*bin_size+1)
            bin_end = min(c(((bin)*bin_size), length(g_chr)))
            if(bin_start > bin_end) { next }
            start_pos = start(g_chr[bin_start]); 
            if(bin_size > 1) { bin_file = full_path(work_folder, paste0("chr",chromosome,"_",start(g_chr[bin_start]),"-",end(g_chr[bin_end]),"_bin",bin,"_",bin_size,"_granges_trimer_counts.csv"))
            } else { bin_file = full_path(work_folder, paste0("chr",chromosome,"_",start(g_chr[bin_start]),"-",end(g_chr[bin_end]),"_trimer_counts.csv")) }
            bin_files <- c(bin_files, bin_file)
            if (!restart && file.exists(bin_file)) { next }
            
            g <- g_chr[bin_start:bin_end]
            tri_counts <- tri_counts_empty
            for(i in 1:length(g)) {
                print(paste0(bin_start," / ",length(g_chr),"    ",i," / ",length(g)))
                grange <- g[i]
                grange_refseq <- toupper(refseq[[1]][(start(grange)-1):(end(grange)+1)])
                trimers <- rollapply(grange_refseq, width=3, by=1, paste0, collapse="")
                b <- table(trimers)
                if(split_by_gc_meth) { 
                    gc_meths <- get_dna_methylation(chromosome, start(grange), end(grange))
                    cpg_trimer_meths <- sapply(cpg_trimer_names, function(cpg_trimer_name) { cpg_trimer_indices <- (trimers == cpg_trimer_name); if(sum(cpg_trimer_indices)>0) { return(mean(gc_meths[cpg_trimer_indices])) } else { return(0) } })
                    for(cpg_trimer_name in names(cpg_trimer_meths[cpg_trimer_meths > 0])) {
                        if(is.na(b[cpg_trimer_name])) { next }
                        methylated_cpg_trimer_name <- gsub("CG", "CmG", cpg_trimer_name) #gsub("^([A-Z]{2})([A-Z])$", "\\1m\\2", cpg_trimer_name)
                        b[methylated_cpg_trimer_name] <- b[cpg_trimer_name]*cpg_trimer_meths[cpg_trimer_name]
                        b[cpg_trimer_name] <- b[cpg_trimer_name]*(1-cpg_trimer_meths[cpg_trimer_name])
                    }
                }
                tri_counts[names(b)] <- tri_counts[names(b)] + b
            }
            #print(sum(as.numeric(width(g))))
            #print(sum(tri_counts))
            #print(tri_counts)
            write.csv(rbind(tri_counts), file=bin_file, row.names=FALSE)
        }
    }
    tri_counts_total <- c()
    for(i in 1:length(bin_files)) {
        print(paste0())
        bin_file = bin_files[i]
        #print(bin_file)
        tri_counts <- as.numeric(read.csv(file=bin_file))
        names(tri_counts) <- tri_counts_names
        #print(tri_counts)
        if(length(tri_counts_total) > 0) { tri_counts_total <- tri_counts_total + tri_counts
        } else { tri_counts_total <- tri_counts  }
    }
    return(tri_counts_total)
}
if(FALSE) {
    tri_counts_TES_imr90 <- get_trimer_counts_for_granges(TES_granges, work_folder=output_path("trimer_counts_for_TES_imr90"))
    tri_counts_TSS_imr90 <- get_trimer_counts_for_granges(TSS_granges, work_folder=output_path("trimer_counts_for_TSS_imr90"))
    #tri_counts_mappable <- get_trimer_counts_for_granges(mappability_granges, work_folder=output_path("trimer_counts_mappable_temp"), bin_size=2000)
    tri_counts_clustered_TES_h9 <- get_trimer_counts_for_granges(clustered_TES_TAD_granges, work_folder=output_path("trimer_counts_for_clustered_TES_h9"))
    tri_counts_clustered_TSS_h9 <- get_trimer_counts_for_granges(clustered_TSS_TAD_granges, work_folder=output_path("trimer_counts_for_clustered_TSS_h9"))
    
    print(sum(as.numeric(width(mappability_granges))))
    print(sum(tri_counts_mappable))
    tri_counts_mappable
    write.csv(rbind(tri_counts_mappable), file=output_path("mappable_trimer_counts.csv"), row.names=FALSE)
}

# Gets Refseq trimer strings for the specified chromosomes and positions. Refs can also be specified for sanity check. If Alts are specified, will return trimers mutations (ex: "AAA->ACA") for each element instead of just reference trimers.
# The width parameter defaults to 3, but can be optionally set to any odd number >=3 to return width-mers (ex: 5-mers or 7-mers) centered on the reference base, rather than trimers.
get_trimers <- function(chromosomes, positions, refs=NULL, alts=NULL, width=3, version="hg19", gc_meths=FALSE, allow_composite=FALSE, allow_1bp_misalignment=FALSE, mismatches_pause=TRUE, allow_BSgenome=FALSE) {
    chromosomes <- paste0(chromosomes); positions <- as.numeric(paste0(positions))
    return_env <- new.env()
    if (!(width %% 2 == 1 && width > 3)) { width=3 }
    cent = ceiling(width/2)

    # chromosomes <- rep(paste0("."), sum(as.numeric(width(g))))
    # sapply(1:length(g), function(i) { rep(seqnames(g[i]), width(g[i])) }) 
    # positions <- rep(paste0("."), sum(as.numeric(width(g))))
    # 
    #positions_dat <- unlist(sapply(1:length(g), function(i) { print(paste0(i," / ",length(g))); grange <- g[i]; return(rbind(paste0(seqnames(grange)), seq(start(grange), end(grange)))) }))
    #chromosomes <- positions_dat[seq(1,length(positions_dat),by=2)]
    #positions <- positions_dat[seq(2,length(positions_dat),by=2)]
    # trimers <- BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg19::Hsapiens, chromosomes, positions-floor(width/2), positions+floor(width/2), as.character = T)
    # 
    # a <- BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg19::Hsapiens, paste0("chr",seqnames(g[1])), start(g[1]), end(g[1]), as.character = T)
    # 
    # g <- mappability_granges[1:1000]
    # tri_counts <- rep(0, 65)
    # names(tri_counts) <- c(apply(permutations(n=4, r=3, v=c("A", "T", "C", "G"), repeats.allowed=TRUE), 1, function(x) { paste(x, collapse="") }), "NNN")
    # for(i in 1:length(g)) {
    #     print(paste0(i," / ",length(g))); grange <- g[i]
    #     b <- table(rollapply(toupper(refseq[[1]][(start(grange)-1):(end(grange)+1)]), width=3, by=1, paste0, collapse=""))
    #     tri_counts[names(b)] <- tri_counts[names(b)] + b
    # }
    # sum(as.numeric(width(g)))
    # sum(tri_counts)
    # tri_counts
    #BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg19::Hsapiens, chromosomes, positions-floor(width/2), positions+floor(width/2), as.character = T)
    
    trimers <- rep(paste0("."), length(chromosomes))
    if(gc_meths) { gc_meths_vector <- rep(0, length(chromosomes)) }
    for(chromosome in unique(chromosomes)) {
        cat(paste0("Reading chr", gsub("chr","",chromosome), "..."))
        refseq <- get_refseq(chromosome, version, allow_BSgenome=allow_BSgenome)
        cat("Done.\nProcessing positions...\n")
        curr_chrom_indices <- which(chromosomes == chromosome)
        curr_positions <- positions[curr_chrom_indices]
        #if(allow_BSgenome) { trimers_for_chrom <- data.frame(sapply(1:length(curr_chrom_indices), function(i) { toupper(unlist(strsplit(refseq[[1]][(curr_positions[i]-floor(width/2)):(curr_positions[i]+floor(width/2))], ""))) }),stringsAsFactors=FALSE)
        trimers_for_chrom <- data.frame(sapply(1:length(curr_chrom_indices), function(i) { toupper(refseq[[1]][(curr_positions[i]-floor(width/2)):(curr_positions[i]+floor(width/2))]) }),stringsAsFactors=FALSE)
        cat("Done.\n")
    
        if(!is.null(refs)) { 
            cat("Sanity checking Refseq reference sequence...")
            curr_refs <- paste0(refs[curr_chrom_indices])
            refseq_refs <- paste0(trimers_for_chrom[cent,]) # unlist(strsplit(trimers_for_chrom, ""))[seq(2, (3*length(trimers_for_chrom)), by=3)]
            insertion_indices <- (nchar(curr_refs) > 1)
            curr_refs_to_compare <- curr_refs
            refseq_refs_to_compare <- refseq_refs
            if(sum(insertion_indices) > 0) { curr_refs_to_compare[insertion_indices] <- substr(curr_refs_to_compare[insertion_indices], 1, 2); refseq_refs_to_compare[insertion_indices] <- apply(data.frame(trimers_for_chrom[(cent:(cent+1)),insertion_indices]), 2, paste0, collapse="") }
            
            all_equal = tryCatch({ all.equal(curr_refs_to_compare, refseq_refs_to_compare) })
            if (inherits(all_equal, "try-error") || !(is.boolean(all_equal) && all_equal)) { 
                ref_mismatches <- sapply(1:length(curr_refs), function(i) { curr_refs_to_compare[i] != refseq_refs_to_compare[i] })
                mismatches <- data.frame(cbind(rep(chromosome, sum(ref_mismatches)), curr_positions[ref_mismatches], curr_refs[ref_mismatches], apply(data.frame(trimers_for_chrom[,ref_mismatches]), 2, paste0, collapse="")))
                colnames(mismatches) <- c("Chrom", "Position", "Ref", "refseq_trimer")
                if(width > 3) { colnames(mismatches)[ncol(mismatches)] <- gsub("tri", width, colnames(mismatches)[ncol(mismatches)]) }
                if(allow_composite) { 
                    refseq_refs_to_compare[insertion_indices] <- apply(data.frame(trimers_for_chrom[((cent-1):cent),insertion_indices]), 2, paste0, collapse="")
                    composite_matches_indices <- sapply(which(ref_mismatches), function(i) { curr_refs_to_compare[i] == sequence_composite(refseq_refs_to_compare[i]) })
                    composite_matches <- mismatches[composite_matches_indices,]; mismatches <- mismatches[!composite_matches_indices,]
                    
                    if(nrow(composite_matches)>0) {
                        if(is.null(return_env[["composite_matches"]])) { return_env[["composite_matches"]] <- composite_matches
                        } else { return_env[["composite_matches"]] <- rbind(return_env[["composite_matches"]], composite_matches) }
                    }
                }
                if(allow_1bp_misalignment) { 
                    # refseq_refs_to_compare[insertion_indices] <- apply(data.frame(trimers_for_chrom[((cent-1):cent),insertion_indices]), 2, paste0, collapse="")
                    # misaligned_matches_indices <- sapply(which(ref_mismatches), function(i) { curr_refs_to_compare[i] == sequence_composite(refseq_refs_to_compare[i]) })
                    # misaligned_matches <- mismatches[misaligned_matches_indices,]; mismatches <- mismatches[!misaligned_matches_indices,]
                    # 
                    # if(nrow(composite_matches)>0) {
                    #     if(is.null(return_env[["misaligned_matches"]])) { return_env[["misaligned_matches"]] <- misaligned_matches
                    #     } else { return_env[["misaligned_matches"]] <- rbind(return_env[["misaligned_matches"]], misaligned_matches) }
                    # }
                }
                if(nrow(mismatches)>0) {
                    if(is.null(return_env[["mismatches"]])) { return_env[["mismatches"]] <- mismatches 
                    } else { return_env[["mismatches"]] <- rbind(return_env[["mismatches"]], mismatches) }
                }
            }
            cat("Done.\n")
        }
        if(is.null(alts)) {
            apply_direction = 2;
            trimers[curr_chrom_indices] <- apply(data.frame(trimers_for_chrom), apply_direction, paste0, collapse="")
        } else {
            cat("Converting trimers to trimer mutations using alts...")
            if(width == 3) {
                trimers[curr_chrom_indices] <- paste0(trimers_for_chrom[1,], trimers_for_chrom[2,], trimers_for_chrom[3,], "->", trimers_for_chrom[1,], alts[curr_chrom_indices], trimers_for_chrom[3,])
            } else { 
                trimers[curr_chrom_indices] <- paste0(apply(data.frame(trimers_for_chrom), 2, paste0, collapse=""), "->", apply(data.frame(trimers_for_chrom[(1:(cent-1)),]), 2, paste0, collapse=""), alts, apply(data.frame(trimers_for_chrom[((cent+1):width),]), 2, paste0, collapse=""))
            }
            cat("Done.\n")
        }
        
        if(gc_meths) { gc_meths_vector[curr_chrom_indices] <- get_dna_methylation(chromosome, curr_positions) }
    }
    
    return_env[["trimers"]] <- trimers
    if(gc_meths) { return_env[["gc_meths"]] <- gc_meths_vector }
    if(!is.null(return_env[["composite_matches"]])) { rownames(return_env[["composite_matches"]]) <- c() }
    if(!is.null(return_env[["mismatches"]])) { rownames(return_env[["mismatches"]]) <- c(); print("ERROR: Not all provided refs matched Refseq refs!!! Check 'mismatches' key in return environment for details."); if(mismatches_pause) { readline() } }
    return(return_env)
}
#get_trimers_test <- get_trimers(ssc$Chrom[ssc$Chrom %in% c(1,10)], ssc$Position[ssc$Chrom %in% c(1,10)], ssc$Ref[ssc$Chrom %in% c(1,10)], width=9, version="hg38")
#get_trimers_test[["mismatches"]]

# Returns the methylated CpG versions of the provided trimer strings. Works for both trimers and trimer mutations.
methylate_cpg_trimer <- function(tris) {
    if (grepl("/", tris)) { parts <- strsplit(tris, "/")[[1]]; return(paste0(methylate_cpg_trimer(parts[1]), "/", methylate_cpg_trimer(parts[2]))) }
    return(gsub("^(.?)CG", "\\1CmG", tris))
}

# Returns the base-paired composite of the given DNA sequences. Can handle "m" character denoting CpG methylation, as well as ignoring "->" for mutation entries.
sequence_composite <- function(sequences, flip_direction=TRUE, remove_na = TRUE) {
    base_composites <- c("T", "A", "G", "C", "N", "m" , "m", "-", "?", "*"); names(base_composites) <- c("A", "T", "C", "G", "N", "M", "m", "-", "?", "*")
    return(unname(sapply(sequences, function(sequence) { 
        if(grepl("->", sequence)) { sequence_parts <- strsplit(sequence, "->")[[1]]; return(paste0(sequence_composite(sequence_parts), collapse="->")) }
        is_upper = (toupper(gsub("m", "", sequence)) == gsub("m", "", sequence))
        sequence <- toupper(sequence)
        sequence_comp <- unname(base_composites[unlist(strsplit(sequence, ""))])
        if(remove_na) { sequence_comp <- sequence_comp[!is.na(sequence_comp)] } else { sequence_comp[is.na(sequence_comp)] <- "?" }
        if(flip_direction) { sequence_comp <- rev(sequence_comp) }
        sequence_comp <- paste0(sequence_comp, collapse="")
        if(!is_upper) { sequence_comp <- tolower(sequence_comp) }
        return(sequence_comp)
    })))
}
#genome_trimer_counts <- (tri.counts.genome-tri.counts.exome) # tri.counts.genome data frame is from deconstructSigs package (source("https://bioconductor.org/biocLite.R"); biocLite("BSgenome"); biocLite("BSgenome.Hsapiens.UCSC.hg19") to get Bioconductor dependencies)
#trimer_pairs_mappings <- c(row.names(genome_trimer_counts)); names(trimer_pairs_mappings) <- sequence_composite(trimer_pairs_mappings)
#expanded_trimer_pairs_mappings <- c(trimer_pairs_mappings, paste0(c("A", "T", "C", "G"), "CmG")); names(expanded_trimer_pairs_mappings) <- sequence_composite(expanded_trimer_pairs_mappings)

# Gets all trimer mutation pair strings mappings, mapping to strings using trimers present in trimer_pairs_mappings. For example, "AAA->ACA" maps to "TTT->ACA".
get_trimer_mutations_pairs_mappings <- function() {
    mutations <- permutations(n=4, r=2, v=c("A", "T", "C", "G"))
    contexts <- permutations(n=4, r=2, v=c("A", "T", "C", "G"), repeats.allowed=TRUE)
    ref_trimers <- c(sapply(1:nrow(mutations), function(i) { paste0(contexts[,1], mutations[i,1], contexts[,2]) }))
    ref_cpg_indices <- grepl("CG", ref_trimers)
    ref_trimers <- c(ref_trimers, gsub("^(.?)CG", "\\1CmG", ref_trimers[ref_cpg_indices]))
    ref_trimers_mapped <- ref_trimers; ref_trimers_mapped[!(ref_trimers %in% expanded_trimer_pairs_mappings)] <- expanded_trimer_pairs_mappings[ref_trimers[!(ref_trimers %in% expanded_trimer_pairs_mappings)]]
    alt_trimers <- c(sapply(1:nrow(mutations), function(i) { paste0(contexts[,1], mutations[i,2], contexts[,2]) })); 
    alt_trimers <- c(alt_trimers, alt_trimers[ref_cpg_indices])
    alt_trimers_mapped <- alt_trimers; alt_trimers_mapped[!(alt_trimers %in% expanded_trimer_pairs_mappings)] <- expanded_trimer_pairs_mappings[alt_trimers[!(alt_trimers %in% expanded_trimer_pairs_mappings)]]
    trimer_mutations_pairs_mappings <- new.env()
    trimer_mutations_pairs_mappings <- data.frame(cbind(paste0(ref_trimers, "->", alt_trimers), paste0(ref_trimers_mapped, "->", alt_trimers_mapped)), stringsAsFactors=FALSE)
    colnames(trimer_mutations_pairs_mappings) <- c("trimer_mutation", "mapped_trimer_mutation")
    row.names(trimer_mutations_pairs_mappings) <- trimer_mutations_pairs_mappings$trimer_mutation
    return(trimer_mutations_pairs_mappings)
}
#trimer_mutations_pairs_mappings <- get_trimer_mutations_pairs_mappings()
#write.csv(trimer_mutations_pairs_mappings, file=output_path("trimer_mutations_pairs_mappings.csv"), row.names=FALSE)

estimate_trimer_mutation_rates_from_data <- function(dat, tri_counts, sample_subset=NULL) {
    original_sample_count = length(unique(dat$sample))
    if (!is.null(sample_subset)) { dat <- dat[dat$sample %in% unique(dat$sample)[1:sample_subset],] }
    sample_count = length(unique(dat$sample))
    
    if("Chrom" %in% colnames(dat)) { dat <- dat[dat$snv_indel == "snv",]; chr_colname="Chrom"; pos_colname="Position"; ref_colname="Ref"; alt_colname="Alt"
    } else { chr_colname="X.chr"; pos_colname="pos"; ref_colname="ref"; alt_colname="alt" }
    
    trimer_muts <- get_trimers(dat[,chr_colname], dat[,pos_colname], dat[,ref_colname], dat[,alt_colname], gc_meths=TRUE)
    if(!is.null(trimer_muts[["mismatches"]])) { return(trimer_muts) }
    gc_meths <- trimer_muts[["gc_meths"]]
    trimer_muts <- to_both_strands(trimer_muts[["trimers"]])
    #trimer_muts <- trimer_mutations_pairs_mappings[trimer_muts,"mapped_trimer_mutation"]
    
    #if (nrow(genome_trimer_counts) == 32) { genome_trimer_counts <- rbind(genome_trimer_counts, genome_trimer_counts); row.names(genome_trimer_counts)[(nrow(genome_trimer_counts)/2+1):nrow(genome_trimer_counts)] <- sequence_composite(row.names(genome_trimer_counts)[(nrow(genome_trimer_counts)/2+1):nrow(genome_trimer_counts)]) } # If data does not include counts for composites.
    #trimer_pairs_mappings <- row.names(genome_trimer_counts); names(trimer_pairs_mappings) <- sequence_composite(row.names(genome_trimer_counts))
    #trimers_mapped <- trimer_pairs_mappings[trimers]; trimers[!is.na(trimers_mapped)] <- trimers_mapped[!is.na(trimers_mapped)]
    
    tri_muts <- new.env()
    apply_out <- sapply(to_both_strands(unique(trimer_mutations_pairs_mappings[,1])), function(x) { tri_muts[[x]] <- rep(0, length(unique(dat$sample))); names(tri_muts[[x]]) <- paste0(unique(dat$sample)) })
    
    cpg_trimer_names <- to_both_strands(gsub("m", "", trimer_mutations_pairs_mappings[grepl("m", trimer_mutations_pairs_mappings[,1]),1]))
    cpg_dat <- data.frame()
    sapply_out <- sapply(unique(dat$sample), function(samp) { 
        samp_indices <- (dat$sample == samp)
        mutated_trimer_counts <- table(trimer_muts[samp_indices])
        
        cpg_trimer_meths <- sapply(cpg_trimer_names, function(cpg_trimer_name) { cpg_trimer_indices <- samp_indices & (trimer_muts == cpg_trimer_name); if(sum(cpg_trimer_indices)>0) { return(mean(gc_meths[cpg_trimer_indices])) } else { return(0) } })
        for(cpg_trimer_name in names(cpg_trimer_meths[cpg_trimer_meths > 0])) {
            if(is.na(mutated_trimer_counts[cpg_trimer_name])) { next }
            methylated_cpg_trimer_name <- methylate_cpg_trimer(cpg_trimer_name) #trimer_mutations_pairs_mappings[methylate_cpg_trimer(cpg_trimer_name),"mapped_trimer_mutation"]
            mutated_trimer_counts[methylated_cpg_trimer_name] <- mutated_trimer_counts[cpg_trimer_name]*cpg_trimer_meths[cpg_trimer_name]
            mutated_trimer_counts[cpg_trimer_name] <- mutated_trimer_counts[cpg_trimer_name]*(1-cpg_trimer_meths[cpg_trimer_name])
            cpg_dat <<- rbind(cpg_dat, cbind(cpg_trimer_name, paste0(samp), cpg_trimer_meths[cpg_trimer_name]))
        }
        for(tri_mut in names(mutated_trimer_counts)) {
            #for(tri_mut in trimer_mutations_pairs_mappings$trimer_mutation[trimer_mutations_pairs_mappings$mapped_trimer_mutation == mapped_tri_mut]) {
                if(is.null(tri_muts[[tri_mut]])) { tri_muts[[tri_mut]] <- rep(0, length(unique(dat$sample))); names(tri_muts[[tri_mut]]) <- paste0(unique(dat$sample)) }
                tri_muts[[tri_mut]][paste0(samp)] <- tri_muts[[tri_mut]][paste0(samp)] + mutated_trimer_counts[tri_mut]
            #}
        }
    })
    colnames(cpg_dat) <- c("mutation", "sample", "gc_meth")
    cpg_dat$gc_meth <- as.numeric(paste0(cpg_dat$gc_meth))
    rownames(cpg_dat) <- c()
    cpg_dat <- cpg_dat[order(cpg_dat$mutation, cpg_dat$sample, cpg_dat$gc_meth),]
    
    #empirical_trimer_mutation_rates <- sapply(ls(tri_muts), function(tri_mut) { tri_muts[[tri_mut]]/genome_trimer_counts[strsplit(tri_mut, "->")[[1]][1],] })
    #empirical_trimer_mutation_rates <- sapply(ls(tri_muts), function(tri_mut) { tri_muts[[tri_mut]]/sum(c(tri_counts_TES[gsub("m", "", strsplit(tri_mut, "->")[[1]][1])], tri_counts_TES[sequence_composite(gsub("m", "", strsplit(tri_mut, "->")[[1]][1]))])) })
    empirical_trimer_mutation_rates <- sapply(ls(tri_muts), function(tri_mut) { tri_muts[[tri_mut]]/sum(c(tri_counts[strsplit(tri_mut, "->|/")[[1]][c(1,3)]])) })
    print(colMeans(empirical_trimer_mutation_rates[,gsub("m","",colnames(empirical_trimer_mutation_rates)) %in% cpg_trimer_names]))
    print(paste0("ACG/CGT trimer count: ", sum(tri_counts[c("ACG","CGT")])))
    print(paste0("ACmG/CmGT trimer count: ", sum(tri_counts[c("ACmG","CmGT")])))
    print(paste0("mean(tri_muts[['ACG->ATG/CGT->CAT']]): ", mean(tri_muts[["ACG->ATG/CGT->CAT"]])))
    print(paste0("mean(tri_muts[['ACmG->ATG/CmGT->CAT']]): ", mean(tri_muts[["ACmG->ATG/CmGT->CAT"]])))
    print(paste0("ACG->ATG/CGT->CAT count in SSC: ", sum(trimer_muts == "ACG->ATG/CGT->CAT")))
    print(paste0("ACG->ATG/CGT->CAT avg gc methylation: ", mean(gc_meths[trimer_muts == "ACG->ATG/CGT->CAT"])))
    print(paste0("ACG->ATG/CGT->CAT avg gc methylation from cpg_dat: ", mean(cpg_dat$gc_meth[cpg_dat$mutation == "ACG->ATG/CGT->CAT"])))
    print(colMeans(empirical_trimer_mutation_rates[,gsub("m","",colnames(empirical_trimer_mutation_rates)) == "ACG->ATG/CGT->CAT"]))
    print(paste0("ACG->ATG/CGT->CAT avg gc methylation from averaged cpg_dat samples (with zero-ing): ", mean(sapply(paste0(unique(cpg_dat$sample)), function(samp) { gc = mean(cpg_dat$gc_meth[cpg_dat$sample == samp & cpg_dat$mutation == "ACG->ATG/CGT->CAT"]); if(is.nan(gc)) { gc = 0 }; return(gc) } ))))
    print(paste0("ACG->ATG/CGT->CAT avg gc methylation from averaged cpg_dat samples (w/o zero-ing): ", mean(sapply(paste0(unique(cpg_dat$sample)), function(samp) { gcs = cpg_dat$gc_meth[cpg_dat$sample == samp & cpg_dat$mutation == "ACG->ATG/CGT->CAT"]; return(mean(gcs[!is.nan(gcs)])) } ))))
    
    # Re-calculating from raw data to avoid Simpson's Paradox!
    empirical_trimer_mutation_rates_aggregate <- sapply(colnames(empirical_trimer_mutation_rates), function(mutation) {
        methylated = grepl("m", mutation)
        original_trimer_mutation = gsub("m", "", mutation)
        if(methylated) { return(sum(gc_meths[trimer_muts == original_trimer_mutation]) / (sample_count * sum(c(tri_counts[strsplit(mutation, "->|/")[[1]][c(1,3)]]))))
        } else { return(sum(1 - gc_meths[trimer_muts == original_trimer_mutation]) / (sample_count * sum(c(tri_counts[strsplit(mutation, "->|/")[[1]][c(1,3)]])))) }
    })
    names(empirical_trimer_mutation_rates_aggregate) <- colnames(empirical_trimer_mutation_rates)
    
    #genome_trimer_counts
        
    #summary(m1 <- glm(num_awards ~ prog + math, data=p, family="poisson"))
    
    return_env <- new.env()
    return_env[["per_sample"]] <- empirical_trimer_mutation_rates
    return_env[["aggregate"]] <- empirical_trimer_mutation_rates_aggregate
    return(return_env)
}
if (FALSE) {
    mappable_trimer_counts <- read.csv(data_path("mappable_trimer_counts.csv"))
    b <- sapply(rownames(genome_trimer_counts), function(tri) { mappable_trimer_counts[,tri] + mappable_trimer_counts[,trimer_pairs_mappings[tri]] })
    mappable_trimer_counts <- cbind(mappable_trimer_counts, trimer_pairs_mappings[names(mappable_trimer_counts)])
    
    pdf(file=output_path("trimer_prevalence.pdf"))
    barplot(t(tri.counts.genome/sum(tri.counts.genome)), main="Genomic Trimer Prevalence", ylab="prevalence", xaxt="n", cex.lab=1.2, col="grey")
    #points(c(1:length(unique(a$trimer_mutation))), sqrt(sort(colMeans(trimer_mutation_rates_from_ssc))), pch = 1, col = "blue", cex=0.4)
    labs <- rownames(tri.counts.genome) #rbind(cbind(rownames(tri.counts.genome), tri.counts.genome), t(sapply(1:nrow(tri.counts.genome), function(i) { tri = rownames(tri.counts.genome)[i]; return(c(names(trimer_pairs_mappings)[trimer_pairs_mappings == tri], tri.counts.genome[i,1])) })))
    axis(1, at=1.2*c(1:length(labs))-0.5, las=2, family="mono", labels=labs, cex.axis=0.5)
    #points(c(1:length(labs)), sapply(labs, function(tri_mut) { sqrt(tri_mutation_rates[[tri_mut]]) }), pch = 1, col = "red", cex=0.4)
    #legend("topleft", legend=c("no_decomposition", "decomposition"), col=c("blue","red"), pch=15)
    dev.off()
    
    pdf(file=output_path("trimer_prevalence_cpg_decomposed.pdf"))
    trimer_prevalence <- tri_counts_TES[!grepl("N|m", names(tri_counts_TES))]/sum(tri_counts_TES)
    trimer_prevalence <- rbind(trimer_prevalence, trimer_prevalence); rownames(trimer_prevalence) <- c("unmethylated CpG", "methylated CpG")
    trimer_prevalence[2,][gsub("m", "", names(tri_counts_TES[grepl("m", names(tri_counts_TES)) & !grepl("N", names(tri_counts_TES))]))] <- tri_counts_TES[grepl("m", names(tri_counts_TES)) & !grepl("N", names(tri_counts_TES))]/sum(tri_counts_TES)
    trimer_prevalence[2,][!(colnames(trimer_prevalence) %in% gsub("m", "", names(tri_counts_TES[grepl("m", names(tri_counts_TES)) & !grepl("N", names(tri_counts_TES))])))] <- 0
    for(i in 1:length(trimer_pairs_mappings)) { 
        tri = names(trimer_pairs_mappings)[i]
        trimer_prevalence[,tri] <- rowSums(trimer_prevalence[,c(tri,trimer_pairs_mappings[tri])])
        colnames(trimer_prevalence)[colnames(trimer_prevalence) == tri] <- paste0(sort(c(tri,trimer_pairs_mappings[tri])),collapse="/")
    }
    trimer_prevalence <- trimer_prevalence[,grepl("/",colnames(trimer_prevalence))]
    trimer_prevalence <- trimer_prevalence[,order(colnames(trimer_prevalence))] #trimer_prevalence[,order(colSums(trimer_prevalence))]
    barplot(trimer_prevalence, main="Genomic Trimer Prevalence", ylab="prevalence", xaxt="n", cex.lab=1.2, col=c("lightblue", "coral1"))
    axis(1, at=1.2*c(1:length(labs))-0.5, las=2, family="mono", labels=colnames(trimer_prevalence), cex.axis=0.8)
    mtext("CpG methylation decomposed")
    legend("topright", legend=c("methylated CpG"), col=c("coral1"), pch=15)
    dev.off()
    
    pdf(file=output_path("trimer_mutation_prevalence_cpg_decomposed.pdf"))
    tri_mut_counts_TES <- c(trimer_mutation_rates_from_ssc_means[,2]); names(tri_mut_counts_TES) <- trimer_mutation_rates_from_ssc_means[,1]
    trimer_prevalence <- tri_mut_counts_TES[!grepl("N|m", names(tri_mut_counts_TES))]/sum(tri_mut_counts_TES)
    trimer_prevalence <- rbind(trimer_prevalence, trimer_prevalence); rownames(trimer_prevalence) <- c("unmethylated CpG", "methylated CpG")
    trimer_prevalence[2,][gsub("m", "", names(tri_mut_counts_TES[grepl("m", names(tri_mut_counts_TES)) & !grepl("N", names(tri_mut_counts_TES))]))] <- tri_mut_counts_TES[grepl("m", names(tri_mut_counts_TES)) & !grepl("N", names(tri_mut_counts_TES))]/sum(tri_mut_counts_TES)
    trimer_prevalence[2,][!(colnames(trimer_prevalence) %in% gsub("m", "", names(tri_mut_counts_TES[grepl("m", names(tri_mut_counts_TES)) & !grepl("N", names(tri_mut_counts_TES))])))] <- 0
    for(i in 1:nrow(trimer_mutations_pairs_mappings)) { 
        tri = trimer_mutations_pairs_mappings$trimer_mutation[i]
        mapped_tri = trimer_mutations_pairs_mappings$mapped_trimer_mutation[i]
        if(tri %in% colnames(trimer_prevalence) && mapped_tri %in% colnames(trimer_prevalence)) {
            trimer_prevalence[,tri] <- rowSums(trimer_prevalence[,c(tri,mapped_tri)])
            colnames(trimer_prevalence)[colnames(trimer_prevalence) == tri] <- mapped_tri #paste0(sort(c(tri,mapped_tri)),collapse="/")
        }
    }
    #trimer_prevalence <- trimer_prevalence[,grepl("/",colnames(trimer_prevalence))]
    trimer_prevalence <- trimer_prevalence[,order(colnames(trimer_prevalence))] #trimer_prevalence[,order(colSums(trimer_prevalence))]
    barplot(trimer_prevalence, main="Empirical Trimer Mutation Prevalence", ylab="prevalence", xaxt="n", cex.lab=1.2, col=c("lightblue", "coral1"))
    axis(1, at=1.2*c(1:ncol(trimer_prevalence))-0.5, las=2, family="mono", labels=colnames(trimer_prevalence), cex.axis=0.4)
    mtext("CpG methylation decomposed")
    legend("topright", legend=c("methylated CpG"), col=c("coral1"), pch=15)
    dev.off()
    
    
    # trimer_mutation_rates_from_ssc_means <- apply(trimer_mutation_rates_from_ssc, 2, mean) #[order(row.names(trimer_mutation_rates_from_ssc))]
    # trimer_mutation_rates_from_ssc_means <- trimer_mutation_rates_from_ssc_means / 2
    # trimer_mutation_rates_from_ssc_means <- merge(data.frame(names(trimer_mutation_rates_from_ssc_means), trimer_mutation_rates_from_ssc_means), trimer_mutations_pairs_mappings, by.x="names.trimer_mutation_rates_from_ssc_means.", by.y="mapped_trimer_mutation")[,c("trimer_mutation", "trimer_mutation_rates_from_ssc_means")]
    # 
    # #trimer_mutation_rates_from_felix <- sapply(names(trimer_mutation_rates_from_ssc_means), function(mapped_tri_mut) { original_tri_muts <- trimer_mutations_pairs_mappings$trimer_mutation[trimer_mutations_pairs_mappings$mapped_trimer_mutation == mapped_tri_mut]; return(mean(c(tri_mutation_rates[[original_tri_muts[1]]], tri_mutation_rates[[original_tri_muts[2]]]))) })
    # trimer_mutation_rates_from_felix <- sapply(paste0(trimer_mutation_rates_from_ssc_means[,1]), function(tri_mut) { tri_mutation_rates[[tri_mut]] })
    # 
    # trimer_mutation_rates_estimates <- cbind(trimer_mutation_rates_from_ssc_means, trimer_mutation_rates_from_felix)
    # #trimer_mutation_rates_estimates <- data.frame(cbind(names(trimer_mutation_rates_from_ssc_means)[order(names(trimer_mutation_rates_from_ssc_means))], trimer_mutation_rates_from_ssc_means[order(names(trimer_mutation_rates_from_ssc_means))], trimer_mutation_rates_from_felix[order(names(trimer_mutation_rates_from_felix))]))
    # colnames(trimer_mutation_rates_estimates) <- c("trimer", "model_rate", "samocha_rate"); trimer_mutation_rates_estimates$model_rate <- as.numeric(paste0(trimer_mutation_rates_estimates$model_rate)); trimer_mutation_rates_estimates$samocha_rate <- as.numeric(paste0(trimer_mutation_rates_estimates$samocha_rate))
    # #trimer_mutation_rates_estimates <- merge(trimer_mutation_rates_estimates, trimer_mutations_pairs_mappings, by.x="trimer", by.y="mapped_trimer_mutation")[,c("trimer_mutation", "model_rate", "samocha_rate")]
    # 
    # mean(trimer_mutation_rates_estimates$model_rate / trimer_mutation_rates_estimates$samocha_rate)
    # cor(trimer_mutation_rates_estimates$model_rate, trimer_mutation_rates_estimates$samocha_rate, method="spearman")
    # pdf(file=output_path("ssc_vs_samocha_trimer_rates.pdf"))
    # scatterplot(trimer_mutation_rates_estimates$model_rate, trimer_mutation_rates_estimates$samocha_rate, main="SSC vs. Samocha Trimer Rates", xlab="ssc rate", ylab="samocha rate", cex.main=1.5, cex.axis=1.4, cex.lab=1.4, spread=FALSE)
    # most_different <- order(abs(log(trimer_mutation_rates_estimates$model_rate / trimer_mutation_rates_estimates$samocha_rate)), decreasing=TRUE)[1:20]
    # text(trimer_mutation_rates_estimates$model_rate[most_different], trimer_mutation_rates_estimates$samocha_rate[most_different], paste0("                           ", trimer_mutation_rates_estimates$trimer[most_different]), cex=0.75) 
    # #mtext("200 random genes (675.2 Kbp, exonic); model samples every 30th trimer", padj=-1.05, cex=1)
    # legend("topleft", legend=c("Smoothed Line", "Regression Line", paste0("Spearman = ", round(cor(trimer_mutation_rates_estimates$model_rate, trimer_mutation_rates_estimates$samocha_rate, method="spearman"),3))), col=c("red","green","white"), lty=c(1,1), bty="n")
    # dev.off()
    # write.csv(trimer_mutation_rates_estimates, file=output_path("ssc_trimer_mutation_rates.csv"), row.names=FALSE, quote=FALSE)
    # write.csv(trimer_mutation_rates_from_ssc, file=output_path("ssc_trimer_mutation_rates_per_sample.csv"), row.names=TRUE, quote=FALSE)
    # 
    # apply(trimer_mutation_rates_estimates[c("AAA->AGA","AAA->ATA","AAA->ACA","TTT->TCT","TTT->TGT","TTT->TAT"),2:3], 2, mean) #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # trimer_mutation_rates_estimates_refs <- unlist(strsplit(trimer_mutation_rates_estimates$trimer, "->"))[seq(1,nrow(trimer_mutation_rates_estimates)*2,by=2)]
    # trimer_mutation_rates_estimates_refs_mapped <- trimer_pairs_mappings[trimer_mutation_rates_estimates_refs]
    # trimer_mutation_rates_estimates_refs[!is.na(trimer_mutation_rates_estimates_refs_mapped)] <- trimer_mutation_rates_estimates_refs_mapped[!is.na(trimer_mutation_rates_estimates_refs_mapped)]
    # ref_background_rates_estimates <- aggregate(trimer_mutation_rates_estimates[,2:3], by=list(trimer_mutation_rates_estimates_refs), FUN=sum); colnames(ref_background_rates_estimates)[1] <- c("trimer")
    # mean(ref_background_rates_estimates$model_rate / ref_background_rates_estimates$samocha_rate)
    # cor(ref_background_rates_estimates$model_rate, ref_background_rates_estimates$samocha_rate, method="spearman")
    # scatterplot(ref_background_rates_estimates$model_rate, ref_background_rates_estimates$samocha_rate, main="SSC vs. Samocha Trimer Rates", xlab="ssc rate", ylab="samocha rate", cex.main=1.5, cex.axis=1.4, cex.lab=1.4, spread=FALSE)
    # most_different <- order(abs(log(ref_background_rates_estimates$model_rate / ref_background_rates_estimates$samocha_rate)), decreasing=TRUE)[1:20]
    # text(ref_background_rates_estimates$model_rate[most_different], ref_background_rates_estimates$samocha_rate[most_different], paste0("                           ", ref_background_rates_estimates$trimer[most_different]), cex=0.75) 
    # #mtext("200 random genes (675.2 Kbp, exonic); model samples every 30th trimer", padj=-1.05, cex=1)
    # legend("topleft", legend=c("Smoothed Line", "Regression Line", paste0("Spearman = ", round(cor(ref_background_rates_estimates$model_rate, ref_background_rates_estimates$samocha_rate, method="spearman"),3))), col=c("red","green","white"), lty=c(1,1), bty="n")
    # write.csv(ref_background_rates_estimates, file=output_path("ssc_trimer_total_mutation_rates.csv"), row.names=FALSE, quote=FALSE)
    # 
    #trimer_mutation_rates_from_ssc <- trimer_mutation_rates_from_ssc[,colMeans(trimer_mutation_rates_from_ssc)>0]
    #colnames(trimer_mutation_rates_from_ssc) <- to_both_strands(colnames(trimer_mutation_rates_from_ssc))
    #trimer_mutation_rates_from_ssc <- trimer_mutation_rates_from_ssc[,!duplicated(colnames(trimer_mutation_rates_from_ssc))]
    
    #trimer_mutation_rates_from_ssc_test <- estimate_trimer_mutation_rates_from_data(ssc_all_snps_TES[1:100,], tri_counts_TES) #estimate_trimer_mutation_rates_from_data(ssc_all)
    #ssc_all_snps_TES[1:100,1:5]
    #trimer_mutation_rates_from_ssc_test[,colMeans(trimer_mutation_rates_from_ssc_test)>0]
    #table(paste0(ssc_all_snps_TES[1:100,3], "->", ssc_all_snps_TES[1:100,4]))
    #colMeans(trimer_mutation_rates_from_ssc_test)[colMeans(trimer_mutation_rates_from_ssc_test)>0]
}
if (FALSE) {
    # Estimate cpg_meth decomposed mutation rate (both aggregate and per sample) from TSS
    tss_trimer_mutation_rates_from_ssc_cpg_meth <- estimate_trimer_mutation_rates_from_data(ssc_all_snps_TSS, tri_counts_TSS_h9) #estimate_trimer_mutation_rates_from_data(ssc_all)
    tss_trimer_mutation_rates_from_ssc_cpg_meth_aggregate <- tss_trimer_mutation_rates_from_ssc_cpg_meth[["aggregate"]]; tss_trimer_mutation_rates_from_ssc_cpg_meth <- tss_trimer_mutation_rates_from_ssc_cpg_meth[["per_sample"]]
    table(paste0(ssc_all_snps_TSS[,3], "->", ssc_all_snps_TSS[,4]))
    colMeans(tss_trimer_mutation_rates_from_ssc_cpg_meth)[colMeans(tss_trimer_mutation_rates_from_ssc_cpg_meth)>0]
    
    # TSS basic, per sample
    tss_trimer_mutation_rates_from_ssc_basic <- tss_trimer_mutation_rates_from_ssc_cpg_meth
    tss_trimer_mutation_rates_from_ssc_basic[,gsub("m","",colnames(tss_trimer_mutation_rates_from_ssc_basic)[grepl("m", colnames(tss_trimer_mutation_rates_from_ssc_basic))])] <- tss_trimer_mutation_rates_from_ssc_basic[,grepl("m", colnames(tss_trimer_mutation_rates_from_ssc_basic))] + tss_trimer_mutation_rates_from_ssc_basic[,gsub("m","",colnames(tss_trimer_mutation_rates_from_ssc_basic)[grepl("m", colnames(tss_trimer_mutation_rates_from_ssc_basic))])]
    tss_trimer_mutation_rates_from_ssc_basic <- tss_trimer_mutation_rates_from_ssc_basic[,!grepl("m", colnames(tss_trimer_mutation_rates_from_ssc_basic))]
    
    # TSS basic, aggregate
    tss_trimer_mutation_rates_from_ssc_basic_aggregate <- tss_trimer_mutation_rates_from_ssc_cpg_meth_aggregate
    tss_trimer_mutation_rates_from_ssc_basic_aggregate[gsub("m","",names(tss_trimer_mutation_rates_from_ssc_basic_aggregate)[grepl("m", names(tss_trimer_mutation_rates_from_ssc_basic_aggregate))])] <- tss_trimer_mutation_rates_from_ssc_basic_aggregate[grepl("m", names(tss_trimer_mutation_rates_from_ssc_basic_aggregate))] + tss_trimer_mutation_rates_from_ssc_basic_aggregate[gsub("m","",names(tss_trimer_mutation_rates_from_ssc_basic_aggregate)[grepl("m", names(tss_trimer_mutation_rates_from_ssc_basic_aggregate))])]
    tss_trimer_mutation_rates_from_ssc_basic_aggregate <- tss_trimer_mutation_rates_from_ssc_basic_aggregate[!grepl("m", names(tss_trimer_mutation_rates_from_ssc_basic_aggregate))]
    
    # Write TSS-derived files
    write.csv(tss_trimer_mutation_rates_from_ssc_basic, file=output_path("trimer_mutation_rates_basic_from_ssc_tss_h9_per_sample.csv"))
    write.csv(tss_trimer_mutation_rates_from_ssc_cpg_meth, file=output_path("trimer_mutation_rates_cpg_meth_from_ssc_tss_h9_per_sample.csv"))
    write.csv(tss_trimer_mutation_rates_from_ssc_basic_aggregate, file=output_path("trimer_mutation_rates_basic_from_ssc_tss_h9.csv"))
    write.csv(tss_trimer_mutation_rates_from_ssc_cpg_meth_aggregate, file=output_path("trimer_mutation_rates_cpg_meth_from_ssc_tss_h9.csv"))
    
    # Estimate cpg_meth decomposed mutation rate (both aggregate and per sample) from TES
    tes_trimer_mutation_rates_from_ssc_cpg_meth <- estimate_trimer_mutation_rates_from_data(ssc_all_snps_TES, tri_counts_TES_h9) #estimate_trimer_mutation_rates_from_data(ssc_all)
    tes_trimer_mutation_rates_from_ssc_cpg_meth_aggregate <- tes_trimer_mutation_rates_from_ssc_cpg_meth[["aggregate"]]; tes_trimer_mutation_rates_from_ssc_cpg_meth <- tes_trimer_mutation_rates_from_ssc_cpg_meth[["per_sample"]]
    table(paste0(ssc_all_snps_TES[,3], "->", ssc_all_snps_TES[,4]))
    colMeans(tes_trimer_mutation_rates_from_ssc_cpg_meth)[colMeans(tes_trimer_mutation_rates_from_ssc_cpg_meth)>0]
    
    # TES basic, per sample
    tes_trimer_mutation_rates_from_ssc_basic <- tes_trimer_mutation_rates_from_ssc_cpg_meth
    tes_trimer_mutation_rates_from_ssc_basic[,gsub("m","",colnames(tes_trimer_mutation_rates_from_ssc_basic)[grepl("m", colnames(tes_trimer_mutation_rates_from_ssc_basic))])] <- tes_trimer_mutation_rates_from_ssc_basic[,grepl("m", colnames(tes_trimer_mutation_rates_from_ssc_basic))] + tes_trimer_mutation_rates_from_ssc_basic[,gsub("m","",colnames(tes_trimer_mutation_rates_from_ssc_basic)[grepl("m", colnames(tes_trimer_mutation_rates_from_ssc_basic))])]
    tes_trimer_mutation_rates_from_ssc_basic <- tes_trimer_mutation_rates_from_ssc_basic[,!grepl("m", colnames(tes_trimer_mutation_rates_from_ssc_basic))]
    
    # TES basic, aggregate
    tes_trimer_mutation_rates_from_ssc_basic_aggregate <- tes_trimer_mutation_rates_from_ssc_cpg_meth_aggregate
    tes_trimer_mutation_rates_from_ssc_basic_aggregate[gsub("m","",names(tes_trimer_mutation_rates_from_ssc_basic_aggregate)[grepl("m", names(tes_trimer_mutation_rates_from_ssc_basic_aggregate))])] <- tes_trimer_mutation_rates_from_ssc_basic_aggregate[grepl("m", names(tes_trimer_mutation_rates_from_ssc_basic_aggregate))] + tes_trimer_mutation_rates_from_ssc_basic_aggregate[gsub("m","",names(tes_trimer_mutation_rates_from_ssc_basic_aggregate)[grepl("m", names(tes_trimer_mutation_rates_from_ssc_basic_aggregate))])]
    tes_trimer_mutation_rates_from_ssc_basic_aggregate <- tes_trimer_mutation_rates_from_ssc_basic_aggregate[!grepl("m", names(tes_trimer_mutation_rates_from_ssc_basic_aggregate))]
    
    # Write TES-derived files
    write.csv(tes_trimer_mutation_rates_from_ssc_basic, file=output_path("trimer_mutation_rates_basic_from_ssc_tes_h9_per_sample.csv"))
    write.csv(tes_trimer_mutation_rates_from_ssc_cpg_meth, file=output_path("trimer_mutation_rates_cpg_meth_from_ssc_tes_h9_per_sample.csv"))
    write.csv(tes_trimer_mutation_rates_from_ssc_basic_aggregate, file=output_path("trimer_mutation_rates_basic_from_ssc_tes_h9.csv"))
    write.csv(tes_trimer_mutation_rates_from_ssc_cpg_meth_aggregate, file=output_path("trimer_mutation_rates_cpg_meth_from_ssc_tes_h9.csv"))
    
    # Read trimer mutation rates from file
    # TSS
    tissue = "imr90" #"h9"
    tss_trimer_mutation_rates_from_ssc_basic <- read.table(output_path(paste0("trimer_mutation_rates_basic_from_ssc_tss_", tissue, "_per_sample.csv")), sep=",", row.names=1, header=TRUE, check.names=FALSE)
    tss_trimer_mutation_rates_from_ssc_cpg_meth <- read.table(output_path(paste0("trimer_mutation_rates_cpg_meth_from_ssc_tss_", tissue, "_per_sample.csv")), sep=",", row.names=1, header=TRUE, check.names=FALSE)
    tss_trimer_mutation_rates_from_ssc_basic_aggregate <- read.table(output_path(paste0("trimer_mutation_rates_basic_from_ssc_tss_", tissue, ".csv")), sep=",", row.names=1, header=TRUE)
    basic_rownames <- rownames(tss_trimer_mutation_rates_from_ssc_basic_aggregate)
    tss_trimer_mutation_rates_from_ssc_basic_aggregate <- as.numeric(paste0(tss_trimer_mutation_rates_from_ssc_basic_aggregate[,1])); names(tss_trimer_mutation_rates_from_ssc_basic_aggregate) <- basic_rownames
    tss_trimer_mutation_rates_from_ssc_cpg_meth_aggregate <- read.table(output_path(paste0("trimer_mutation_rates_cpg_meth_from_ssc_tss_", tissue, ".csv")), sep=",", row.names=1, header=TRUE)
    cpg_meth_rownames <- rownames(tss_trimer_mutation_rates_from_ssc_cpg_meth_aggregate)
    tss_trimer_mutation_rates_from_ssc_cpg_meth_aggregate <- as.numeric(paste0(tss_trimer_mutation_rates_from_ssc_cpg_meth_aggregate[,1])); names(tss_trimer_mutation_rates_from_ssc_cpg_meth_aggregate) <- cpg_meth_rownames
    # TES
    tes_trimer_mutation_rates_from_ssc_basic <- read.table(output_path(paste0("trimer_mutation_rates_basic_from_ssc_tes_", tissue, "_per_sample.csv")), sep=",", row.names=1, header=TRUE, check.names=FALSE)
    tes_trimer_mutation_rates_from_ssc_cpg_meth <- read.table(output_path(paste0("trimer_mutation_rates_cpg_meth_from_ssc_tes_", tissue, "_per_sample.csv")), sep=",", row.names=1, header=TRUE, check.names=FALSE)
    tes_trimer_mutation_rates_from_ssc_basic_aggregate <- read.table(output_path(paste0("trimer_mutation_rates_basic_from_ssc_tes_", tissue, ".csv")), sep=",", row.names=1, header=TRUE)
    tes_trimer_mutation_rates_from_ssc_basic_aggregate <- as.numeric(paste0(tes_trimer_mutation_rates_from_ssc_basic_aggregate[,1])); names(tes_trimer_mutation_rates_from_ssc_basic_aggregate) <- basic_rownames
    tes_trimer_mutation_rates_from_ssc_cpg_meth_aggregate <- read.table(output_path(paste0("trimer_mutation_rates_cpg_meth_from_ssc_tes_", tissue, ".csv")), sep=",", row.names=1, header=TRUE)
    tes_trimer_mutation_rates_from_ssc_cpg_meth_aggregate <- as.numeric(paste0(tes_trimer_mutation_rates_from_ssc_cpg_meth_aggregate[,1])); names(tes_trimer_mutation_rates_from_ssc_cpg_meth_aggregate) <- cpg_meth_rownames

    for (region in c("tes", "tss")) {
        if (region == "tss") { object_prefix = "tss_" } else { object_prefix = "tes_" }
        for (decompose_cpg_meth in c(FALSE, TRUE)) {
            if (decompose_cpg_meth) { trimer_mutation_rates_from_ssc <- get(paste0(object_prefix,"trimer_mutation_rates_from_ssc_cpg_meth")); trimer_mutation_rates_from_ssc_agg <- get(paste0(object_prefix,"trimer_mutation_rates_from_ssc_cpg_meth_aggregate")); file_suffix = paste0("_cpg_meth_",region,"_imr90.pdf") 
            } else { trimer_mutation_rates_from_ssc <- get(paste0(object_prefix,"trimer_mutation_rates_from_ssc_basic")); trimer_mutation_rates_from_ssc_agg <- get(paste0(object_prefix,"trimer_mutation_rates_from_ssc_basic_aggregate")); file_suffix = paste0("_",region,".pdf") }
            a <- data.frame()
            for(lab in names(sort(trimer_mutation_rates_from_ssc_agg))) {
                i = which(colnames(trimer_mutation_rates_from_ssc) == lab)
                a <- rbind(a, cbind(colnames(trimer_mutation_rates_from_ssc)[i], trimer_mutation_rates_from_ssc[,i]))
            }
            colnames(a) <- c("trimer_mutation", "rate"); a$rate <- as.numeric(paste0(a$rate))
            
            pdf(file=output_path(paste0("ssc_trimer_mutation_signatures",file_suffix)))
            boxplot(a$rate~a$trimer_mutation, main="SSC Trimer Mutation Signature", ylab="mutation rate", xaxt="n", lty=1, cex.lab=1.2, cex=0.4, col="grey")
            labs <- unique(a$trimer_mutation) #sapply(unique(a$trimer_mutation), function(tri_mut) { trimer_mutations_pairs_mappings$trimer_mutation[trimer_mutations_pairs_mappings$mapped_trimer_mutation == tri_mut][1] })
            #points(c(1:length(unique(a$trimer_mutation))), sort(colMeans(trimer_mutation_rates_from_ssc)), pch = 1, col = "blue", cex=0.4)
            points(c(1:length(labs)), trimer_mutation_rates_from_ssc_agg[paste0(labs)], pch = 1, col = "blue", cex=0.4)
            labs_snp_type <- to_both_strands(sapply(labs, function(lab) { return(paste0(strsplit(gsub("m","",lab),"")[[1]][c(2,7)], collapse="->")) }))
            labs_cols <- rainbow(length(unique(labs_snp_type))); labs_cols[labs_cols == "#FFFF00FF"] <- "#FFA500FF" # Use rainbow for colors, but switch yellow to orange for better visibility
            for(i in 1:length(unique(labs_snp_type))) { snp_type = unique(labs_snp_type)[i]; axis(1, at=c(1:length(labs))[labs_snp_type == snp_type], las=2, family="mono", labels=labs[labs_snp_type == snp_type], col.axis=labs_cols[i], col.ticks=labs_cols[i], cex.axis=0.4) }
            points(c(1:length(labs)), sapply(labs, function(tri_mut) { samoch = tri_mutation_rates[[strsplit(gsub("m","",tri_mut),"/")[[1]][1]]]; if(is.null(samoch)) { return(-1) } else { return(samoch) }  }), pch = 1, col = "red", cex=0.4)
            legend("topleft", legend=c(paste0("ssc (N = ",SSC_ALL_SAMPLE_COUNT,")"), "ssc_mean", "samocha", unique(labs_snp_type)), col=c("darkgrey","blue","red",labs_cols), pch=c(15,1,1,rep(3,length(unique(labs_snp_type)))))
            if (decompose_cpg_meth) { mtext("CpG methylation decomposed", cex=1) }
            dev.off()
            
            pdf(file=output_path(paste0("ssc_trimer_mutation_signatures_sqrt",file_suffix)))
            boxplot(sqrt(a$rate)~a$trimer_mutation, main="SSC Trimer Mutation Signature", ylab="sqrt(mutation rate)", xaxt="n", lty=1, cex.lab=1.2, cex=0.4, col="grey")
            labs <- unique(a$trimer_mutation) #sapply(unique(a$trimer_mutation), function(tri_mut) { trimer_mutations_pairs_mappings$trimer_mutation[trimer_mutations_pairs_mappings$mapped_trimer_mutation == tri_mut][1] })
            #points(c(1:length(unique(a$trimer_mutation))), sqrt(sort(colMeans(trimer_mutation_rates_from_ssc))), pch = 1, col = "blue", cex=0.4)
            points(c(1:length(labs)), sqrt(trimer_mutation_rates_from_ssc_agg[paste0(labs)]), pch = 1, col = "blue", cex=0.4)
            labs_snp_type <- to_both_strands(sapply(labs, function(lab) { return(paste0(strsplit(gsub("m","",lab),"")[[1]][c(2,7)], collapse="->")) }))
            labs_cols <- rainbow(length(unique(labs_snp_type))); labs_cols[labs_cols == "#FFFF00FF"] <- "#FFA500FF" # Use rainbow for colors, but switch yellow to orange for better visibility
            for(i in 1:length(unique(labs_snp_type))) { snp_type = unique(labs_snp_type)[i]; axis(1, at=c(1:length(labs))[labs_snp_type == snp_type], las=2, family="mono", labels=labs[labs_snp_type == snp_type], col.axis=labs_cols[i], col.ticks=labs_cols[i], cex.axis=0.4) }
            points(c(1:length(labs)), sapply(labs, function(tri_mut) { samoch = tri_mutation_rates[[strsplit(gsub("m","",tri_mut),"/")[[1]][1]]]; if(is.null(samoch)) { return(-1) } else { return(sqrt(samoch)) }  }), pch = 1, col = "red", cex=0.4)
            legend("topleft", legend=c(paste0("ssc (N = ",SSC_ALL_SAMPLE_COUNT,")"), "ssc_mean", "samocha", unique(labs_snp_type)), col=c("darkgrey","blue","red",labs_cols), pch=c(15,1,1,rep(3,length(unique(labs_snp_type)))))
            if (decompose_cpg_meth) { mtext("CpG methylation decomposed", cex=1) }
            dev.off()
            
            pdf(file=output_path(paste0("ssc_vs_samocha_trimer_mutation",file_suffix)))
            labs <- unique(a$trimer_mutation) #sapply(unique(a$trimer_mutation), function(tri_mut) { trimer_mutations_pairs_mappings$trimer_mutation[trimer_mutations_pairs_mappings$mapped_trimer_mutation == tri_mut][1] })
            labs_snp_type <- to_both_strands(sapply(labs, function(lab) { return(paste0(strsplit(gsub("m","",lab),"")[[1]][c(2,7)], collapse="->")) }))
            empiric <- trimer_mutation_rates_from_ssc_agg[paste0(labs)]
            labs <- labs[empiric > 0]
            labs_snp_type <- labs_snp_type[empiric > 0]
            empiric <- empiric[empiric > 0]
            samoch <- sapply(labs, function(tri_mut) { samoch = tri_mutation_rates[[strsplit(gsub("m","",tri_mut),"/")[[1]][1]]]; if(is.null(samoch)) { return(-1) } else { return(samoch) }  })
            plot_width = 1.1*max(c(samoch, empiric))
            scatterplot(empiric, samoch, main="SSC vs. Samocha Trimer Rates", xlab="ssc rate", ylab="samocha rate", xlim=c(0,plot_width), ylim=c(0,plot_width), cex.main=1.5, cex.axis=1.4, cex.lab=1.4, spread=FALSE, smooth=FALSE)
            if (decompose_cpg_meth) { mtext("CpG methylation decomposed", padj=-1.05, cex=1) }
            most_different <- which(empiric > plot_width/4 | samoch > plot_width/4)
            most_different_cols <- labs_cols[labs_snp_type[most_different]]
            #most_different_cols[samoch[most_different] > 0.03 & empiric[most_different] < 0.03] <- "grey"
            text(empiric[most_different], samoch[most_different], paste0("                           ", labs[most_different]), col=most_different_cols, cex=0.5) 
            legend("topleft", legend=c("Regression Line", paste0("Spearman = ", round(cor(empiric, samoch, method="spearman"),3))), col=c("green","white"), lty=c(1,1), bty="n")
            lines(c(-1,1),c(-1,1),lty=3, col="grey")
            dev.off()
            
            cpg_labs <- (gsub("m","",labs) %in% cpg_trimer_names)
            empiric_sum_no_cpg = sum(sort(trimer_mutation_rates_from_ssc_agg)[!cpg_labs])
            samoch_sum_no_cpg = sum(sapply(labs[!cpg_labs], function(tri_mut) { samoch = tri_mutation_rates[[strsplit(gsub("m","",tri_mut),"/")[[1]][1]]]; if(is.null(samoch)) { return(-1) } else { return(samoch) }  }))
            norm_factor = empiric_sum_no_cpg / samoch_sum_no_cpg
            
            trimer_mutation_rates_from_ssc_agg <- trimer_mutation_rates_from_ssc_agg / norm_factor
            a$rate <- a$rate / norm_factor
            
            pdf(file=output_path(paste0("ssc_trimer_mutation_signatures_norm",file_suffix)))
            boxplot(a$rate~a$trimer_mutation, main="SSC Trimer Mutation Signature", ylab="normalized mutation rate", xaxt="n", lty=1, cex.lab=1.2, cex=0.4, col="grey")
            labs <- unique(a$trimer_mutation)
            points(c(1:length(labs)), trimer_mutation_rates_from_ssc_agg[paste0(labs)], pch = 1, col = "blue", cex=0.4)
            labs_snp_type <- to_both_strands(sapply(labs, function(lab) { return(paste0(strsplit(gsub("m","",lab),"")[[1]][c(2,7)], collapse="->")) }))
            labs_cols <- rainbow(length(unique(labs_snp_type))); labs_cols[labs_cols == "#FFFF00FF"] <- "#FFA500FF" # Use rainbow for colors, but switch yellow to orange for better visibility
            names(labs_cols) <- unique(labs_snp_type)
            for(i in 1:length(unique(labs_snp_type))) { snp_type = unique(labs_snp_type)[i]; axis(1, at=c(1:length(labs))[labs_snp_type == snp_type], las=2, family="mono", labels=labs[labs_snp_type == snp_type], col.axis=labs_cols[i], col.ticks=labs_cols[i], cex.axis=0.4) }
            points(c(1:length(labs)), sapply(labs, function(tri_mut) { samoch = tri_mutation_rates[[strsplit(gsub("m","",tri_mut),"/")[[1]][1]]]; if(is.null(samoch)) { return(-1) } else if (tri_mut %in% labs[cpg_labs]) { return(samoch/norm_factor) } else { return(samoch) }  }), pch = 1, col = "red", cex=0.4)
            legend("topleft", legend=c(paste0("ssc (N = ",SSC_ALL_SAMPLE_COUNT,")"), "ssc_mean", "samocha", unique(labs_snp_type)), col=c("darkgrey","blue","red",labs_cols), pch=c(15,1,1,rep(3,length(unique(labs_snp_type)))))
            if (decompose_cpg_meth) { mtext("CpG methylation decomposed", cex=1) }
            dev.off()
            
            pdf(file=output_path(paste0("ssc_vs_samocha_trimer_mutation_norm",file_suffix)))
            labs <- unique(a$trimer_mutation) #sapply(unique(a$trimer_mutation), function(tri_mut) { trimer_mutations_pairs_mappings$trimer_mutation[trimer_mutations_pairs_mappings$mapped_trimer_mutation == tri_mut][1] })
            labs_snp_type <- to_both_strands(sapply(labs, function(lab) { return(paste0(strsplit(gsub("m","",lab),"")[[1]][c(2,7)], collapse="->")) }))
            empiric <- trimer_mutation_rates_from_ssc_agg[paste0(labs)]
            labs <- labs[empiric > 0]
            labs_snp_type <- labs_snp_type[empiric > 0]
            empiric <- empiric[empiric > 0]
            samoch <- sapply(labs, function(tri_mut) { samoch = tri_mutation_rates[[strsplit(gsub("m","",tri_mut),"/")[[1]][1]]]; if(is.null(samoch)) { return(-1) } else if (tri_mut %in% labs[cpg_labs]) { return(samoch/norm_factor) } else { return(samoch) }  })
            plot_width = 1.1*max(c(samoch, empiric))
            scatterplot(empiric, samoch, main="SSC vs. Samocha Trimer Rates", xlab="normalized ssc rate", ylab="normalized samocha rate", xlim=c(0,plot_width), ylim=c(0,plot_width), cex.main=1.5, cex.axis=1.4, cex.lab=1.4, spread=FALSE, smooth=FALSE)
            if (decompose_cpg_meth) { mtext("CpG methylation decomposed", padj=-1.05, cex=1) }
            most_different <- which(empiric > plot_width/4 | samoch > plot_width/4)
            most_different_cols <- labs_cols[labs_snp_type[most_different]]
            #most_different_cols[samoch[most_different] > 0.03 & empiric[most_different] < 0.03] <- "grey"
            text(empiric[most_different], samoch[most_different], paste0("                           ", labs[most_different]), col=most_different_cols, cex=0.5) 
            legend("topleft", legend=c("Regression Line", paste0("Spearman = ", round(cor(empiric, samoch, method="spearman"),3))), col=c("green","white"), lty=c(1,1), bty="n")
            lines(c(-1,1),c(-1,1),lty=3, col="grey")
            dev.off()
            
            # pdf(file=output_path(paste0("ssc_vs_samocha_trimer_mutation_zoom",file_suffix)))
            # labs <- unique(a$trimer_mutation) #sapply(unique(a$trimer_mutation), function(tri_mut) { trimer_mutations_pairs_mappings$trimer_mutation[trimer_mutations_pairs_mappings$mapped_trimer_mutation == tri_mut][1] })
            # labs_snp_type <- to_both_strands(sapply(labs, function(lab) { return(paste0(strsplit(gsub("m","",lab),"")[[1]][c(2,7)], collapse="->")) }))
            # empiric <- trimer_mutation_rates_from_ssc_agg[paste0(labs)]
            # labs <- labs[empiric > 0]
            # labs_snp_type <- labs_snp_type[empiric > 0]
            # empiric <- empiric[empiric > 0]
            # samoch <- sapply(labs, function(tri_mut) { samoch = tri_mutation_rates[[strsplit(gsub("m","",tri_mut),"/")[[1]][1]]]; if(is.null(samoch)) { return(-1) } else { return(samoch) }  })
            # empiric <- empiric[-c(order(samoch, decreasing=TRUE)[1:4])]
            # labs <- labs[-c(order(samoch, decreasing=TRUE)[1:4])]
            # labs_snp_type <- labs_snp_type[-c(order(samoch, decreasing=TRUE)[1:4])]
            # samoch <- samoch[-c(order(samoch, decreasing=TRUE)[1:4])]
            # plot_width = 1.1*max(c(samoch, empiric))
            # scatterplot(empiric, samoch, main="SSC vs. Samocha Trimer Rates", xlab="ssc rate", ylab="samocha rate", xlim=c(0,plot_width), ylim=c(0,plot_width), cex.main=1.5, cex.axis=1.4, cex.lab=1.4, spread=FALSE, smooth=FALSE)
            # if (decompose_cpg_meth) { mtext("CpG methylation decomposed", padj=-1.05, cex=1) }
            # #most_different <- which(empiric > plot_width/4 | samoch > plot_width/4)
            # #most_different_cols <- labs_cols[labs_snp_type[most_different]]
            # #most_different_cols[samoch[most_different] > 2e-8 & empiric[most_different] < 2.5e-8] <- "grey"
            # #text(empiric[most_different], samoch[most_different], paste0("                           ", labs[most_different]), col=most_different_cols, cex=0.5) 
            # legend("topleft", legend=c("Regression Line", paste0("Spearman = ", round(cor(empiric, samoch, method="spearman"),3))), col=c("green","white"), lty=c(1,1), bty="n")
            # lines(c(-1,1),c(-1,1),lty=3, col="grey")
            # dev.off()
        }
    }
    # Do sanity check for SNP type to make sure mutation model worked as intended!
    trimer_mutation_rates_from_ssc_basic <- tes_trimer_mutation_rates_from_ssc_basic
    pdf(file=output_path("mutation_model_snp_frequency_sanity_check.pdf"))
    b <- aggregate(colMeans(trimer_mutation_rates_from_ssc_basic), by=list(unlist(lapply(gsub("m","",colnames(trimer_mutation_rates_from_ssc_basic)), function(x) trimer_mutations_to_snps(x)))), FUN=sum)
    b_names <- b[,1]
    b <- b[,2]; b <- b/sum(b); names(b) <- b_names
    b_raw <- table(paste0(ssc_all_snps_TES$ref, "->", ssc_all_snps_TES$alt))
    names(b_raw) <-to_both_strands(names(b_raw)); b_raw <- aggregate(as.numeric(b_raw), by=list(names(b_raw)), FUN=sum)
    b_raw_names <- b_raw[,1]
    b_raw <- b_raw[,2]; b_raw <- b_raw/sum(b_raw); names(b_raw) <- b_raw_names
    barplot(rbind(b, b_raw), cex.names=0.85, main="SNP Type Sanity Check", ylab="frequency", beside=TRUE, col=c("red","blue"))
    legend("topleft", legend=c("Derived Rate from Mutation Model","Raw SSC SNP counts"), col=c("red","blue"), pch=15)
    dev.off()
    
    # Find outlier sample in SSC data.
    sample_snp_counts <- sapply(unique(ssc_all$sample), function(samp) sum(ssc_all$snv_indel == "snv" & ssc_all$sample == samp)) 
    hist(sample_snp_counts, breaks=20)
    which(sample_snp_counts > 120)
    unique(ssc_all$sample)[which(sample_snp_counts > 120)]
    
    #library("ggplot2")
    library("vioplot")
    for (tri in c("ACG")) { #apply(permutations(n=4, r=3, v=c("A", "T", "C", "G"), repeats.allowed=TRUE), 1, function(x) { paste0(x, collapse="") })) {
        print(tri)
        tri_indices <- (paste0(unlist(data.frame(strsplit(labs,"-"))[1,])) == tri)
        if(sum(tri_indices) == 0) { next }
        tri_a <- a[a$trimer_mutation %in% sapply(tri_labs, function(tri_lab) trimer_mutations_pairs_mappings$mapped_trimer_mutation[trimer_mutations_pairs_mappings$trimer_mutation == tri_lab]),]
        #tri_a <- tri_a[!is.nan(tri_a$rate),]
        tri_mean_rates <- sort(colMeans(trimer_mutation_rates_from_ssc))[tri_indices]
        tri_labs <- labs[tri_indices]
        tri_mutations <- unique(tri_a$trimer_mutation)
        #ggplot(tri_a,aes(x=tri_a$trimer_mutation, tri_a$rate) + geom_violin())
        vioplot(tri_a$rate[tri_a$trimer_mutation == tri_mutations[1]], tri_a$rate[tri_a$trimer_mutation == tri_mutations[2]], tri_a$rate[tri_a$trimer_mutation == tri_mutations[3]], names=tri_labs, col="lightblue")
    }
    
    
    dat <- ssc
    dat$Chrom <- paste0("chr", dat$Chrom)
    deconstructSigs::mut.to.sigs.input(mut.ref = dat, sample.id = "sample", chr = "Chrom", pos = "Position", ref = "Ref", alt = "Alt")
    deconstructSigs::mut.to.sigs.input(mut.ref = mappability_dat, chr = "chromosome", pos = "Position", ref = "Ref", alt = "Alt")
    #TES_regions_exp <- get_expected_mut_counts(TES_regions, precision=120)
}

to_both_strands <- function(tris) {
    tris_comp <- sequence_composite(tris)
    return(sapply(1:length(tris), function(i) { return(paste0(sort(c(tris[i],tris_comp[i])), collapse="/")) }))
}

trimer_mutations_to_snps <- function(tris, sort_strands=TRUE) {
    snps <- c()
    for(tri in tris) {
        if (grepl("/", tri)) { parts <- strsplit(tri, "/")[[1]]; if (sort_strands) { snp = paste0(sort(trimer_mutations_to_snps(parts)), collapse="/") } else { snp = paste0(trimer_mutations_to_snps(parts), collapse="/") }
        } else { snp = gsub("^.(.).->.(.).", "\\1->\\2", tri) }
        snps <- c(snps, snp)
    }
    return(snps)
}

# Normalize a vector (so that direction the same but sum of its elements is 1).
normalize_vector <- function(x) {
    #return(x / sqrt(sum(x**2)))
    return(x / sum(x))
}

# Annotate RNA secondary structure using the RNAplfold tool, modified by Hilal Kazan (Quaid Morris lab, UofT). URL: http://www.cs.toronto.edu/~hilal/rnacontext/
# For each position in the input sequence, get probabilities of 1. pairedness, 2. being in a hairpin loop (H), 3. being in an internal loop (I), 4. being in a multi loop (M), and 5. being in an external region (E).
# Parameter W is the window length, L is the maximum span, and u is the width.
annotate_rna_ss <- function(sequences, sequence_names=NULL, programs=c("H","I","M","E"), W=240, L=160, u=1, as_tensor=FALSE, sequence_length=NULL, work_folder=output_path("RNAplfold_temp"), program_path="/home/local/ARCS/ak3792/programs/RNAplfold") {
    dir.create(work_folder, showWarnings = FALSE)
    sink_filename = full_path(work_folder,"messages.txt")
    fasta_filename = full_path(work_folder,"input.fasta")
    combined_letter_profiles_filename = full_path(work_folder,"combined_profile.txt")
    sequences <- paste0(sequences)
    num_sequences = length(sequences)
    if(is.null(sequence_names)) { sequence_names <- paste0("seq",1:num_sequences) }
    if(num_sequences > 1) { cat("Writing sequence to FASTA...") } else { cat("Writing sequences to FASTA...") }
    write(paste0(">",sequence_names,"\n",sequences), file=fasta_filename)
    cat("Done.\n")
    programs <- sort(intersect(gsub("_RNAplfold|.*/","",programs), c("H","I","M","E")))
    program_output_filenames <- full_path(work_folder, paste0(programs,"_profile.txt"))
    combined_profile_rownames <- c("paired", c("hairpin_loop", "internal_loop", "multi_loop", "external_region")[c("H","I","M","E") %in% programs])
    programs <- full_path(program_path, paste0(programs,"_RNAplfold"))
    cat("Running RNAplfold modules:\n")
    individual_outputs <- mclapply(1:length(programs), function(i) {
        cat(paste0("\tGenerating ",program_output_filenames[i],"...\n"))
        messages <- system(paste0(programs[i]," -W ",W," -L ",L," -u ",u," <",fasta_filename," >",program_output_filenames[i]," 2>",sink_filename), intern=TRUE)
        cat("Done.\n")
    }, mc.cores=detectCores())
    cat("Combining secondary structure outputs...")
    combined_output <- system(paste0("python ",full_path(program_path,"combine_letter_profiles.py "),paste(program_output_filenames, collapse=" ")," 1 ",combined_letter_profiles_filename), intern=TRUE)
    cat("Done.\n")
    combined_profile_height = length(programs) + 2
    combined_profiles_lines <- readLines(combined_letter_profiles_filename)
    combined_profiles <- lapply(1:num_sequences, function(i) {
        sequence_bases <- strsplit(sequences[i], "")[[1]]
        sequence_bases_length = length(sequence_bases)
        combined_profile <- unfactorize(data.frame(t(data.frame(strsplit(combined_profiles_lines[((i-1)*combined_profile_height+2):(i*combined_profile_height)], "\t"))[-c(1),])))
        rownames(combined_profile) <- combined_profile_rownames
        if(as_tensor) { 
            if(is.null(sequence_length)) { sequence_length = sequence_bases_length; padding = 0 } else { padding = sequence_length - sequence_bases_length } 
            if(padding > 0) {
                combined_profile <- cbind(combined_profile, matrix(rep(0,combined_profile_height-1), nrow=combined_profile_height-1, ncol=padding))
                sequence_bases <- c(sequence_bases, rep("",padding))
            } else { 
                combined_profile <- combined_profile[,1:sequence_length] 
                sequence_bases <- sequence_bases[1:sequence_length]
            }
        }
        colnames(combined_profile) <- sequence_bases
        return(combined_profile)
    })
    names(combined_profiles) <- sequence_names
    if(as_tensor) { 
        tensor <- array(0, dim=c(num_sequences, sequence_length, combined_profile_height-1, 1), dimnames = list(sequence_names, 1:sequence_length, combined_profile_rownames, "channel"))
        tensor[,,,1] <- aperm(abind(combined_profiles, along=0), c(1,3,2))
        combined_profiles <- tensor
    }

    return(combined_profiles)
}
#a <- annotate_rna_ss(c("CTGTACGTAGCAGTACA", "TTTCGTATGA", "CGTAGGTAGGAG"))

# Lift over genomic coordinates between from one assembly to another using the UCSC liftOver tool.
liftover <- function(dat, from, to, chr_colname="chromosome", start_colname="start", end_colname=NULL, ref_colname="Ref", alt_colname="Alt", work_folder=OUTPUT_FOLDER, program_path="/home/local/ARCS/ak3792/programs/ucsc/liftOver", confirm_refseq=TRUE, mismatches_pause=TRUE) {
    if(is.null(end_colname)) { end_colname = start_colname }
    original_chr_colname <- paste0(chr_colname,"_",from); original_start_colname <- paste0(start_colname,"_",from); original_end_colname <- paste0(end_colname,"_",from)
    current_assembly_variants_file = full_path(work_folder, paste0("liftOver_", from, "_variants.txt"))
    new_assembly_variants_file = full_path(work_folder, paste0("liftOver_", to, "_variants.txt")) 
    chain_file = gsub("liftOver$", paste0(from,"To",capitalize_string(to),".over.chain"), program_path)
    unmapped_file = gsub(".txt$", "_unmapped.txt", new_assembly_variants_file)
        
    current_assembly_variants <- genomic_coordinates_to_strings(dat[,chr_colname], dat[,start_colname], ends=dat[,end_colname], add_chr_prefix=TRUE)
    write(current_assembly_variants, file=current_assembly_variants_file)
    
    liftOver_output <- system(paste0(program_path, " ", current_assembly_variants_file, " ", chain_file, " ", new_assembly_variants_file, " ", unmapped_file, " -positions"), intern=TRUE)
    
    liftover_successful_variants <- try({ paste0(read.delim(file=new_assembly_variants_file, header=FALSE)[,1]) }, silent=TRUE)
    if(inherits(liftover_successful_variants, "try-error")) { liftover_successful_variants <- c() }
    liftover_failed_variants <- try({ paste0(read.delim(file=unmapped_file, comment.char='#', header=FALSE)[,1]) }, silent=TRUE)
    if(inherits(liftover_failed_variants, "try-error")) { liftover_failed_variants <- c() }
    
    if(length(liftover_failed_variants) + length(liftover_successful_variants) != length(current_assembly_variants) || length(liftover_failed_variants) != sum(current_assembly_variants %in% liftover_failed_variants)) {
        print("ERROR: liftOver successful and failed variant counts do not add up to total queried variants count!"); readline()
    }
    new_dat <- dat[!(current_assembly_variants %in% liftover_failed_variants),]
    new_dat <- cbind(new_dat, new_dat[,chr_colname], new_dat[,start_colname]); colnames(new_dat)[(ncol(new_dat)-1):ncol(new_dat)] <- c(original_chr_colname, original_start_colname)
    new_dat_liftover <- genomic_coordinates_from_strings(liftover_successful_variants)
    new_dat[,chr_colname] <- new_dat_liftover$chromosome 
    new_dat[,start_colname] <- new_dat_liftover$start
    if(end_colname != start_colname && end_colname %in% colnames(new_dat)) { 
        new_dat <- cbind(new_dat, new_dat[,end_colname]); colnames(new_dat)[ncol(new_dat)] <- original_end_colname
        new_dat[,end_colname] <- new_dat_liftover$end 
    }
    #if("proband" %in% colnames(new_dat)) { new_dat$proband <- unlist(lapply(strsplit(paste0(ssc_new_batch$proband), "\\("), function(x) x[[1]][1])) }
    
    if(confirm_refseq && (ref_colname %in% colnames(new_dat)) && (alt_colname %in% colnames(new_dat))) {
        new_dat_trimers <- get_trimers(new_dat[,chr_colname], new_dat[,start_colname], new_dat[,ref_colname], width=9, version=to, allow_composite=TRUE, allow_1bp_misalignment=TRUE, mismatches_pause=mismatches_pause)
        #b <<- new_dat_trimers
        if(!is.null(new_dat_trimers[["composite_matches"]])) { # Mapped to Refseq composite, so Ref should be changed!!!
            new_dat_trimers_composite <- new_dat_trimers[["composite_matches"]]
            new_dat_trimers_composite_variants <- genomic_coordinates_to_strings(new_dat_trimers_composite$Chrom, new_dat_trimers_composite$Position, ends=new_dat_trimers_composite$Position, add_chr_prefix=TRUE)
            #b1 <<- liftover_successful_variants
            #b2 <<- new_dat
            #b2_old <<- dat
            composite_variant_indices <- (liftover_successful_variants %in% new_dat_trimers_composite_variants)
            new_dat[,ref_colname] <- as.character(new_dat[,ref_colname])
            new_dat[,ref_colname] <- as.character(new_dat[,ref_colname])
            new_dat[composite_variant_indices, ref_colname] <- sequence_composite(new_dat[composite_variant_indices, ref_colname])
            new_dat[composite_variant_indices, alt_colname] <- sequence_composite(new_dat[composite_variant_indices, alt_colname])
        }
        if(!is.null(new_dat_trimers[["mismatches"]])) { # Refseq mismatch!!!
            new_dat_trimers_mismatches <- new_dat_trimers[["mismatches"]]
            new_dat_trimers_mismatches_variants <- genomic_coordinates_to_strings(new_dat_trimers_mismatches$Chrom, new_dat_trimers_mismatches$Position, ends=new_dat_trimers_mismatches$Position, add_chr_prefix=TRUE)
            new_dat_trimers_mismatches <- cbind(new_dat_trimers_mismatches, current_assembly_variants[!(current_assembly_variants %in% liftover_failed_variants)][liftover_successful_variants %in% new_dat_trimers_mismatches_variants], liftover_successful_variants[liftover_successful_variants %in% new_dat_trimers_mismatches_variants])
            colnames(new_dat_trimers_mismatches)[(ncol(new_dat_trimers_mismatches)-1):ncol(new_dat_trimers_mismatches)] <- c(paste0("original_", from, "_variant"), paste0("mapped_", to, "_variant"))
            print("Mismatches:")
            print(new_dat_trimers_mismatches)
            
            new_dat <- new_dat[!(liftover_successful_variants %in% new_dat_trimers_mismatches_variants),]
        }
    } else { print("Skipped Refseq confirmation.") }
    
    return(new_dat)
}

# Return, as a GRanges object, genomic ranges for the whole hg19 genome *except* for the regions specified in g. 
invert_granges <- function(g) {
    chr_names <- paste0("chr", c(1:22,"X","Y"))
    chr_lengths <- seqlengths(Hsapiens)[chr_names]
    if (!grepl("chr", paste0(seqnames(g)[1]))) {
        chr_names <- gsub("chr", "", chr_names)
        names(chr_lengths) <- gsub("chr", "", names(chr_lengths))
    }
    
    granges <- c()
    
    print("Processing chromosomes...")
    for(chr_name in chr_names) {
        print(chr_name)
        g_chr <- g[seqnames(g) == chr_name]
        if (length(g_chr) > 0) {
            starts <- c(1, end(g_chr) + 1)
            ends <- c(start(g_chr) - 1, chr_lengths[chr_name])
            dat <- data.frame(rep(chr_name, length(starts)), starts, ends)
            colnames(dat) <- c("chromosome", "start", "end")
            if (length(granges) > 0) { granges <- c(granges, to_genomic_regions(dat, keep_chr_prefix=TRUE))
            } else { granges <- to_genomic_regions(dat, keep_chr_prefix=TRUE) }
        }
    }
    print("Done.")
    print("Removing 0-length elements...")
    #dat <- dat[as.numeric(paste0(dat$start)) <= as.numeric(paste0(dat$end)),]
    granges <- granges[start(granges) <= end(granges)]
    print("Done.")
    return(granges)
}
#bad_mappability_granges <- invert_granges(mappability_granges)

# Return DNAStringSet object constrcted from the given GRanges.
granges_to_DNAStringSet <- function(gr, length=NULL) {
    chromosomes <- paste0("chr",gsub("chr","",seqnames(gr)))
    if(is.null(length)) {
        sequences <- BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg19::Hsapiens, chromosomes, start(gr), end(gr), strand=strand(gr))
    } else {
        padding <- floor((length - width(gr)) / 2)
        starts <- start(gr) - padding
        sequences <- BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg19::Hsapiens, chromosomes, starts, starts+length-1, strand=strand(gr))
    }
    return(sequences)
}
#gr <- load_annotation("K562.RBFOX2")[3925:3934]

# Return DNAStringSet object constrcted from the given sequences and chromosomes.
sequences_to_DNAStringSet <- function(sequences, chromosomes) {
    dna_string_set <- DNAStringSet(sequences)
    names(dna_string_set) <- paste0("chr",gsub("chr","",paste0(chromosomes)))
    return(dna_string_set)
}

rename_col <- function(dat, from, to) {
	#gsub("^Chrom$", "chr", colnames(cases), perl=TRUE)
	col_indices <- which(colnames(dat) %in% from)
	if (length(col_indices) > 0) { colnames(dat)[col_indices] <- to }
	return(dat)
}

split_multiple_alts <- function(dat, alt_colname="Alt", linked_colname=NULL) {
    if(sum(grepl(",", dat[,alt_colname]) > 0)) {
        library("tidyr")
        alt_colname_index = which(colnames(dat) == alt_colname)
        colnames(dat)[alt_colname_index] <- "Alt"
        if(!is.null(linked_colname)) {
            linked_colname_index = which(colnames(dat) == linked_colname)
            colnames(dat)[linked_colname_index] <- "linked_colname"
            dat <- separate_rows(dat, Alt, linked_colname)
            colnames(dat)[linked_colname_index] <- linked_colname
        }
        else { dat <- separate_rows(dat, Alt) }
        colnames(dat)[alt_colname_index] <- alt_colname
    }
    return(dat)
}
#dat <- cbind(gnomad, unlist(lapply(observed, function(x) paste0(x,collapse=","))))
#colnames(dat)[ncol(dat)] <- "observed"
#table(unlist(lapply(observed, length)))
#which(unlist(lapply(observed, length)) > 3)
#split_multiple_alts(gnomad)


standardize_colnames <- function(dat, re_order=FALSE, remove_chr_prefix=FALSE) {
	dat <- rename_col(dat, from=c("chromosome", "X.chr", "chr", "Chr", "CHROM"), to="Chrom")
	dat <- rename_col(dat, from=c("CHROM_hg19"), to="Chrom_hg19")
	dat <- rename_col(dat, from=c("CHROM_hg38"), to="Chrom_hg38")
	dat <- rename_col(dat, from=c("Pos_hg38"), to="Position_hg38")
	dat <- rename_col(dat, from=c("pos", "POSITION", "POS", "Pos", "Start"), to="Position")
	dat <- rename_col(dat, from=c("POS_hg19"), to="Position_hg19")
	dat <- rename_col(dat, from=c("POS_hg38"), to="Position_hg38")
	dat <- rename_col(dat, from=c("ref", "REF"), to="Ref")
	dat <- rename_col(dat, from=c("alt", "ALT"), to="Alt")
	dat <- rename_col(dat, from=c("Blinded.ID", "Blinded_ID", "SampleID", "Individual.ID"), to="sample")
	dat <- rename_col(dat, from=c("Type"), to="snv_indel")
	if(remove_chr_prefix && "Chrom" %in% colnames(dat)) { dat$Chrom <- gsub("chr", "", paste0(dat$Chrom)) } 
	if(re_order) { 
	    cols_to_reorder <- c("Chrom", "Position", "Ref", "Alt", "sample", "snv_indel")
	    dat_colnames <- colnames(dat)
	    cols_to_reorder <- cols_to_reorder[cols_to_reorder %in% dat_colnames]
	    dat_colnames_reorder_indices <- which(dat_colnames %in% cols_to_reorder)
	    if(length(dat_colnames_reorder_indices) == length(cols_to_reorder)) {
	        dat <- dat[,c(cols_to_reorder,dat_colnames[-c(dat_colnames_reorder_indices)])]
	    }
	}
	return(dat)
}

capitalize_string <- function(x) {
    s <- strsplit(x, " ")[[1]]
    return(paste(toupper(substring(s, 1,1)), substring(s, 2), sep="", collapse=" "))
}

round_to_nearest <- function(x, nearest=1, round_0.5_up=TRUE) {
    x <- x/nearest
    if(round_0.5_up) { x <- x + (nearest/10) }
    return(round(x)*nearest)
}

cat("Successfully loaded alex_suite.R!\n\nRun help function alex_suite() to get a list of all available functions in this suite.\nAlternatively, run alex_suite([QUERY]) to find only functions containing keyword [QUERY].\n\nRun update_os() to refresh environment variables for your current operating system.\n")


#try-catch structure
#num_pages <- try(as.numeric(paste0(xpathSApply(html_data, "((//*[@id='content']/div//td[@class='paging-list-index'])[1]/a)[last()]", xmlValue))))
#if(inherits(num_pages, "try-error") || is.na(num_pages) ) { num_pages = 1 }


