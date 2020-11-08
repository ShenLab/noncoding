library(keras)
library(kerasR)
library("tfprobability")
gpu <- tf$config$experimental$list_physical_devices("GPU")[[1]]
tf$config$experimental$set_memory_growth(gpu, TRUE) # Important step on some servers to uncap the GPU memory usage right away. 
library(stringr)
library(pbapply)
#library(EBImage)
library(rhdf5)
library(gtools)
#library("PWMEnrich")
#library("TFBSTools")
#library("seqLogo")
#library("msa")
library("ComplexHeatmap")
library("hexbin")
#library("colorspace")
library("RColorBrewer")
library("ggplot2")
library("tfprobability")
library("tidyverse")
library("stringi")
library("parallel")
library("data.table")
library("abind")
library("reticulate")
# //TODO: Not all of the above are dependencies any more, and will be cleaned up.

source("alex_suite.R")
source("suprnova_data_setup.R")
source("suprnova_design.R")

#################################################################################################################
# SIMPLE FUNCTIONS
#################################################################################################################

#################################################################################################################
# MAIN: Super-models code and analysis
#################################################################################################################

if(file.exists(output_path("rbp_order.rds"))) { rbp_order <- readRDS(output_path("rbp_order.rds"))
} else { rbp_order <- order_rbps(); saveRDS(rbp_order, output_path("rbp_order.rds")) }
load_supermodel_training_data("transcribed100k_1"); dat_obs_exp[,1] <- rep(1,nrow(dat_obs_exp))
source("suprnova_design.R")
s_anchor = log(5e-5) #log(mean((expected/af_tensor)[af_tensor > 0 & af_tensor < 1e-3]))
s_upper_bound = 0.005
s_scaling_factor = (log(s_upper_bound) - s_anchor)
cat(paste0("s = exp(s_model * s_scaling_factor + s_anchor)\ns_anchor = ",round(s_anchor,3),", s_scaling_factor = ",round(s_scaling_factor,3),"\n"),paste0("s_model = ",c(-1,0,1),"  ->  s = ",exp((c(-1,0,1)) * s_scaling_factor + s_anchor),"\n"))
gc(); model <- define_supermodel_architecture(L=15, num_filters=5, rbp_kernel_width=3, sequence_kernel_width=1, custom_conv1=TRUE, s_anchor=s_anchor, s_scaling_factor=s_scaling_factor)
model_layers <- lapply(model$layers, function(x) paste0("Output Dim: c(", gsub("[()]*","",gsub("list","",paste(gsub("NULL","N",x$get_output_shape_at(as.integer(0))),collapse=", "))),")"))
names(model_layers) <- lapply(model$layers, function(x) x$name)
model_layers[1:15]
training_results <- train_supermodel(num_epochs=50, batch_size=256, rbp_order=rbp_order, fraction_for_training=0.8)
#maps_result <- plot_model_outputs(model, test_indices[rowSums(af_tensor)[test_indices]>0][1:100])
#weighted_conv_altref <- unfactorize(data.frame(cbind(t(as.matrix(model$get_layer("conv_alt_ref")$weights[[1]][1,1,1,,])), as.matrix(model$get_layer("rbp_binding_dense")$weights[[1]])[,1])))
#rownames(weighted_conv_altref) <- c(paste0("filter_",1:nrow(weighted_conv_altref))); colnames(weighted_conv_altref) <- c("Alt-Ref Delta", "Ref", "weight")
#weighted_conv_altref <- weighted_conv_altref[order(abs(weighted_conv_altref$weight), decreasing=TRUE),]
#weighted_conv_altref
#model$get_layer("conv_ref_binding")$weights[[1]]
#model$get_layer("conv_delta")$weights[[1]]
#model$get_layer("conv_alt_ref")$weights[[1]]
dim(model$get_layer("rbp_binding_activation")$weights[[1]]); range(as.matrix(model$get_layer("rbp_binding_activation")$weights[[1]]))
dim(model$get_layer("rbp_complex_activation")$weights[[1]]); range(as.matrix(model$get_layer("rbp_complex_activation")$weights[[1]]))
dim(model$get_layer("gene_regulatory_disruption")$weights[[1]]); range(as.matrix(model$get_layer("gene_regulatory_disruption")$weights[[1]]))
#model$get_layer("rbp_complex_impacts_adjust")$weights[[1]]
#model$get_layer("gene_constraint_activation")$weights[[1]]

#save_supermodel(model, "suprnova")
#k_clear_session(); model <- load_supermodel("suprnova", L=15, num_filters=5, s_anchor=s_anchor, s_scaling_factor=s_scaling_factor); gc()
#model_s <- get_activation_model(model, c("d", "s"))
plot_model_distributions <- function(test_idx=NULL) {
    if(is.null(test_idx)) { test_idx = 1:nrow(altref_gradcams) }; set_global(test_idx)
    AFs_all <- af_tensor[test_idx,]; set_global(AFs_all)
    expected_all <- expected[test_idx,]; set_global(expected_all)
    model_s <- get_activation_model(model, c("rbp_binding_activation", "altref_binding_changes", "rbp_complex_activation", "central_site_select", "gene_regulatory_disruption", "d", "selection_coef_output", "s"))
    print("Finding model activations...")
    gc(); preds <- model_predict(model_s, test_idx); set_global(preds)
    s_vals_all <- preds[["s"]]; s_vals_all[is.nan(s_vals_all)] <- 0; set_global(s_vals_all)
    
    print("Plotting model distributions...")
    filename=output_path(paste0("model_distributions_preds.pdf")); pdf(file=filename, width=10)
    par(mfrow=c(2,3))
    mtext_cex = 0.75
    par_mar <- par()$mar; par(mar=par_mar+c(-1,0,0,0))
    print(paste0("range(rbp binding activations) = c(",paste(range(preds[["rbp_binding_activation"]]),collapse=", "),")"))
    #filename=output_path(paste0("model_distributions_rbp_binding.pdf")); pdf(file=filename)
    plot(density(preds[["rbp_binding_activation"]]), lwd=1, col="blue", main="", xlab="RBP binding activation", cex.main=1.3, cex.lab=1.4, cex.axis=1.4)
    mtext(paste0(round(4/3 * 100 * sum(preds[["rbp_binding_activation"]] > 0.5) / prod(dim(preds[["rbp_binding_activation"]])),2)," % of variants have binding activation > 0.5"), cex=mtext_cex, line=1)
    #dev.off(); pdf_to_png(filename)
    print(paste0("range(ref-alt binding changes) = c(",paste(range(preds[["altref_binding_changes"]]),collapse=", "),")"))
    #filename=output_path(paste0("model_distributions_binding_change.pdf")); pdf(file=filename)
    plot(density(preds[["altref_binding_changes"]]), lwd=1, col="blue", main="", xlab="ref-alt activated binding change", cex.main=1.8, cex.lab=1.4, cex.axis=1.4)
    title("SUPRNOVA Prediction Distributions", line=3, cex=1.8)
    mtext(paste0(round(4/3 * 100 * sum(preds[["altref_binding_changes"]] > 0.25) / prod(dim(preds[["altref_binding_changes"]])),2)," % of variants have delta > 0.25 (LOF), "), cex=mtext_cex, line=1.5)
    mtext(paste0(round(4/3 * 100 * sum(preds[["altref_binding_changes"]] < -0.25) / prod(dim(preds[["altref_binding_changes"]])),2)," % have < -0.25 (GOF)"), cex=mtext_cex, line=0.5)
    #dev.off(); pdf_to_png(filename)
    print(paste0("range(rbp complex disruptions) = c(",paste(range(preds[["central_site_select"]]),collapse=", "),")"))
    #filename=output_path(paste0("model_distributions_rbp_complexes.pdf")); pdf(file=filename)
    plot(density(sample(preds[["central_site_select"]], 1000000, replace=TRUE)), lwd=1, col="blue", main="", xlab="RBP complex disruption", cex.main=1.3, cex.lab=1.4, cex.axis=1.4)
    mtext(paste0(round(4/3 * 100 * sum(preds[["central_site_select"]] > 0.25) / prod(dim(preds[["central_site_select"]])),2)," % of complexes have delta > 0.25 (LOF), "), cex=mtext_cex, line=1.5)
    mtext(paste0(round(4/3 * 100 * sum(preds[["central_site_select"]] < -0.25) / prod(dim(preds[["central_site_select"]])),2)," % have < -0.25 (GOF)"), cex=mtext_cex, line=0.5)
    #dev.off(); pdf_to_png(filename)
    par(mar=par_mar+c(-0.25,0,-0.75,0))
    print(paste0("range(d) = c(",paste(range(preds[["d"]]),collapse=", "),")"))
    #filename=output_path(paste0("model_distributions_d.pdf")); pdf(file=filename)
    plot(density(preds[["d"]]), lwd=1, col="blue", main="", xlab="gene reg. damagingness d", cex.main=1.3, cex.lab=1.4, cex.axis=1.4)
    mtext(paste0(round(4/3 * 100 * sum(preds[["d"]] > 0.75) / prod(dim(preds[["d"]])),2)," % of variants have d > 0.75"), cex=mtext_cex, line=1)
    #dev.off(); pdf_to_png(filename)
    print(paste0("range(model_s) = c(",paste(range(preds[["selection_coef_output"]]),collapse=", ")))
    #filename=output_path(paste0("model_distributions_s_coef.pdf")); pdf(file=filename)
    plot(density(preds[["selection_coef_output"]]), lwd=1, col="blue", main="", xlab="in-model selection coef.", cex.main=1.3, cex.lab=1.4, cex.axis=1.4)
    mtext(paste0("s = exp(s_model * s_scaling_factor + s_anchor), "), cex=mtext_cex, line=1.5)
    mtext(paste0("s_anchor = ",round(s_anchor,3),", s_scaling_factor = ",round(s_scaling_factor,3)), cex=mtext_cex, line=0.5)
    #dev.off(); pdf_to_png(filename)
    print(paste0("range(s) = c(",paste(range(preds[["s"]]),collapse=", "),"), mean(s) = ",mean(preds[["s"]])))
    #filename=output_path(paste0("model_distributions_s.pdf")); pdf(file=filename)
    plot(density(log10(s_vals_all)), lwd=1, col="blue", main="", xlab="log10(selection coef. s)", cex.main=1.3, cex.lab=1.4, cex.axis=1.4)
    mtext(paste0("mean(s) = ",mean(s_vals_all),", "), cex=mtext_cex, line=1.5)
    mtext(paste0("mean(mu/AF) = ",mean((expected_all/AFs_all)[AFs_all > 0])), cex=mtext_cex, line=0.5)
    #dev.off(); pdf_to_png(filename)
    dev.off(); pdf_to_png(filename)
    par(mar=par_mar)
    par(mfrow=c(1,1))
    return(0)
}
plot_model_distributions(test_idx=1:10000)

rbp_binding_thresholds <- as.matrix(model$get_layer("rbp_binding_activation")$weights[[1]])[,1]; names(rbp_binding_thresholds) <- rbps[rbp_order]
sort(rbp_binding_thresholds); mean(rbp_binding_thresholds); range(rbp_binding_thresholds)
dat_trimers[8,75:77]
rbp_binding_thresholds[91:100]
dat_gradcams[8,7:9,rbp_order[91:100],]
preds[["rbp_binding_activation"]][8,7:9,91:100,c(1,4)]
altref_gradcams[8,7:9,rbp_order[91:100],]
round(preds[["altref_binding_changes"]][8,7:9,91:100,],3)
preds[["altref_binding_changes"]][8,7:9,91:100,c(1,4)]
sum(preds[["altref_binding_changes"]] > 0); prod(dim(preds[["delta_binding_activation"]]))
range(preds[["rbp_complex_activation"]])
preds[["rbp_complex_activation"]][8,7:9,91:100,,c(1,4)]

conv1_weights <- as.array(model$get_layer("rbp_complex_activation")$weights[[1]][1,1,,1,])
rbp_complex_names <- rollapply(c(NA,names(rbp_order),NA), width=3, by=1, FUN=function(x) paste0(x, collapse=", ")) #x[!is.na(x)]
rbp_complex_filter_weights <- as.matrix(model$get_layer("gene_regulatory_disruption")$weights[[1]] %>% k_reshape(c(160,5))) # Consider making a 2D Dense custom layer
rbp_complex_filter_names <-  paste0("filter",1:ncol(rbp_complex_filter_weights),"_",apply(conv1_weights, 2, paste0, collapse=""))
rownames(rbp_complex_filter_weights) <- rbp_complex_names; colnames(rbp_complex_filter_weights) <- rbp_complex_filter_names
colnames(conv1_weights) <- rbp_complex_filter_names
conv1_weights
rbp_complex_weights <- rbp_complex_filter_weights %*% t(conv1_weights)
range(rbp_complex_weights)
rbp_complex_filter_weights[1:3,]
t(conv1_weights)
rbp_complex_weights[1:3,]

rbp_complex_activations <- preds[["central_site_select"]]; dim(rbp_complex_activations); range(rbp_complex_activations)
rbp_complex_activation_names <- apply(expand.grid(gsub("^.*_","",rbp_complex_filter_names), rbp_complex_names), 1, function(x) {
    x <- strsplit(x[2],", ")[[1]][strsplit(x[1],"")[[1]] == "1"]
    return(paste0(x[x!="NA"], collapse=", "))
})
dimnames(rbp_complex_activations)[[3]] <- rbp_complex_activation_names
rbp_complex_activations[8,1,order(abs(rbp_complex_activations[8,1,]))]
rbp_complex_activations[8,2,order(abs(rbp_complex_activations[8,2,]))]
rbp_complex_activations[8,3,order(abs(rbp_complex_activations[8,3,]))]
rbp_complex_activations[8,4,order(abs(rbp_complex_activations[8,4,]))]
# NEXT, CONSIDER PUTTING A SIGMOID OVER THE ACTIVATIONS TO MODIFY THE DISTRIBUTION OF d. AND FIX THE CONNECTION TO s!
##################STOPPED HERE!!!###################
dim(preds[["gene_regulatory_disruption"]]); range(preds[["gene_regulatory_disruption"]])

# Investigate relationship of region types vs variant score
investigate_region_type_relationship <- function(region_types, scores, bucket_cuts=NULL, num_buckets=5, quant_start=0, quant_end=0.95, score_type="d", normalize_counts=TRUE, normalize_buckets=TRUE, include_class_any=NULL, tilted_labels=TRUE, scientific_notation=FALSE, draw_legend=TRUE, legend_loc="topleft", legend_ncol=1, legend_cex=0.8, class_label_cex=0.65, class_label_x_offset = 0.5, srt=60, verbose=FALSE) {
    scores_dims <- dim(scores); if(!is.null(scores_dims) && length(scores_dims) == 2) { scores <- apply(scores, 1, max) }
    if(is.null(include_class_any)) { include_class_any = normalize_counts }
    if(is.null(bucket_cuts)) { quant_by=(quant_end-quant_start)/(num_buckets-1) ; bucket_cuts <- c(-Inf,round(quantile(scores, seq(quant_start+quant_by,quant_end,by=quant_by)),2)) }
    if(is.null(dim(bucket_cuts))) { bucket_cuts <- rbind(bucket_cuts, rep(Inf,length(bucket_cuts))) }
    #bucket_names <- c(paste0(bucket_cuts[1,]," < ",score_type," <= ",bucket_cuts[2,]))
    bucket_names <- rep("",ncol(bucket_cuts))
    b <- Reduce(cbind, lapply(1:ncol(bucket_cuts), function(bucket_i) {
        if(verbose) { print(bucket_i) }
        #if(s_i == 0) { curr_indices_backup <- s_vals_all == 0 } else if(s_i > ncol(s_bucket_cuts)) { curr_indices_backup <- rep(TRUE, length(s_vals_all)) } else {
        lower_bound = bucket_cuts[1,bucket_i]; upper_bound = bucket_cuts[2,bucket_i]
        if(upper_bound < Inf) { bucket_name <- paste0(lower_bound," < ",score_type," <= ",upper_bound)
        } else { if(lower_bound > -Inf) { bucket_name <- paste0(score_type," > ",lower_bound) } else { bucket_name <- paste0(score_type," = any") } }
        bucket_names[bucket_i] <<- bucket_name
        curr_indices <- scores > lower_bound & scores <= upper_bound
        return(colSums(region_types[curr_indices,]))
    })); colnames(b) <- bucket_names
    if(include_class_any) { b <- rbind(colSums(b), b); rownames(b)[1] <- "any class" }
    b_counts <- b
    #[,!grepl("= any",colnames(b_counts))]
    region_variant_proportions <- rowSums(b_counts); region_variant_proportions <- region_variant_proportions / sum(region_variant_proportions)
    if(normalize_counts || normalize_buckets) {
        things_to_norm_by <- c()
        if(normalize_buckets) { b <- apply(b, 2, function(x) x / sum(x)); things_to_norm_by <- c(things_to_norm_by, paste0(score_type," bucket")) }
        if(normalize_counts) { b <- b / region_variant_proportions; things_to_norm_by <- c(things_to_norm_by, "class") }
        things_to_norm_by <- rev(things_to_norm_by)
        mtext_text = paste0("Normalization done by total variant count per ",paste(things_to_norm_by, collapse=" and ")); ylab_text = "Normalized Proportion"
    } else { mtext_text = "No normalization performed"; ylab_text = "Proportion" }
    #return(b)
    if(tilted_labels) { b_backup <- b; colnames(b) <- NULL }
    cols <- c(rep("gray",include_class_any), brewer.pal(nrow(b)-include_class_any, "Spectral"))
    filename=paste0(output_path(paste0(score_type,"_region_type_comparison.pdf")))
    if(normalize_counts) { filename <- gsub("\\.pdf$", "_normcounts.pdf", filename) }
    if(normalize_buckets) { filename <- gsub("\\.pdf$", "_normbuckets.pdf", filename) }
    pdf(filename)
    par_mar <- par()$mar; par(mar=par_mar+c(0,1.5,0,0))
    main_text = "Regional Class Distribution vs. "
    if(score_type == "d") { main_text = paste0(main_text,"Gene Reg. Damagingness d")
    } else if(score_type == "s") { main_text = paste0(main_text,"Selection Coef. s")
    } else { main_text = paste0(main_text,score_type) }
    barplot(b, beside=TRUE, space=c(0,1), las=2, ylim=c(0, max(b)*1.04), main=main_text, ylab="", col=cols) #, ylim = c(0,5 + max(mtcars$qsec)), xlab = "", space = 1)
    title(ylab=ylab_text, mgp=c(4,1,0), cex.lab=1.2)
    title(xlab=paste0(score_type," bucket"), mgp=c(3.15,1,0), cex.lab=1.3)
    if(tilted_labels) {
        class_label_starts = seq(nrow(b)/2 + class_label_x_offset, (nrow(b)+1)*ncol(b) + class_label_x_offset, by=nrow(b)+1)
        text(class_label_starts, par("usr")[3]-(max(b)*0.05), srt=srt, adj=1, xpd=TRUE, labels=bucket_names, cex=class_label_cex)  #rotate 60 degrees (srt = 60)
    }
    mtext(mtext_text)
    if(draw_legend) { legend(legend_loc, legend=c(rownames(b)), col=cols, ncol=legend_ncol, pch=15, cex=legend_cex) }
    dev.off()
    pdf_to_png(filename)
    par(mar=par_mar)
    if(tilted_labels) { b <- b_backup }
    write.csv(b_counts, output_path(paste0(score_type,"_region_type_comparison_counts.csv")))
    return(b)
}

# Investigate relationship of variant score densities vs. region types
investigate_score_vs_region_type_relationship <- function(region_types, scores, score_type="d", righttail_threshold=NULL, righttail_quantile=0.5, logscale=NULL, include_class_any=FALSE, verbose=FALSE) {
    scores_dims <- dim(scores); if(!is.null(scores_dims) && length(scores_dims) == 2) { scores <- apply(scores, 1, max) }
    cols <- c(rep("gray",include_class_any), brewer.pal(ncol(region_types), "Spectral"))
    include_class_any=FALSE
    verbose=FALSE
    scores_dims <- dim(scores); if(!is.null(scores_dims) && length(scores_dims) == 2) { scores <- apply(scores, 1, max) }
    cols <- c(rep("gray",include_class_any), brewer.pal(ncol(region_types), "Spectral"))
    righttail_quantile=0.5
    logscale=NULL
    if(is.null(logscale)) { logscale = score_type == "s" }
    xlab_text = score_type; if(logscale) { xlab_text = paste0("log10(",xlab_text,")") }
    main_text = "Regional Class Distribution vs. "
    if(score_type == "d") { main_text = paste0(main_text,"Gene Reg. Damagingness d")
    } else if(score_type == "s") { main_text = paste0(main_text,"Selection Coef. s")
    } else { main_text = paste0(main_text,score_type) }
    for(righttail in c(FALSE,TRUE)) {
        filename=paste0(output_path(paste0(score_type,"_distributions_by_region_type.pdf")))
        if(righttail) { if(is.null(righttail_threshold)) { righttail_threshold = quantile(scores_density$x,righttail_quantile) }; filename <- gsub("\\.pdf$", "_righttail.pdf", filename); mtext_text = paste0("Right tail of distributions only, ",score_type," > ",round(righttail_threshold,3))
        } else { mtext_text = "Full distributions" }
        pdf(filename)
        plot_i_start = 1 - include_class_any; plot_indices <- plot_i_start:ncol(region_types); legend_text <- colnames(region_types)[plot_indices]
        for(plot_i in plot_indices) {
            if(verbose) { print(plot_i) }
            if(plot_i == 0) { curr_indices <- rep(TRUE, nrow(region_types)); legend_text <- c("any class", legend_text)
            } else { curr_indices <- which(region_types[,plot_i] == TRUE) }
            scores_density <- density(scores[curr_indices], from=0)
            if(plot_i == plot_i_start) {
                domain = pmax(c(1.1*range(scores_density$x)),c(max(righttail*c(righttail_threshold,0)),0))
                plot(scores_density, col=cols[plot_i+include_class_any], lwd=2, ylim=c(0,max(scores_density$y[scores_density$x > domain[1] & scores_density$x < domain[2]])*1.1), xlim=domain, xaxs="i", main=paste0("Distributions of ",gsub("^.*vs\\. ", "", main_text)," by Region Class"), xlab=xlab_text, cex.main=1.3, cex.lab=1.4, cex.axis=1.4)
            } else { lines(scores_density, col=cols[plot_i+include_class_any], lwd=2) }
        }
        mtext(mtext_text, cex=1.2)
        legend("topright", legend=legend_text, col=cols, ncol=1, pch=15, cex=1.2)
        dev.off()
        pdf_to_png(filename)
    }
}

# Test canonical splicing sites and region types
if(file.exists(output_path("rbp_order.rds"))) { rbp_order <- readRDS(output_path("rbp_order.rds"))
} else { rbp_order <- order_rbps(); saveRDS(rbp_order, output_path("rbp_order.rds")) }
s_anchor = -8.365
s_upper_bound = 0.01
s_scaling_factor = (log(s_upper_bound) - s_anchor)
cat(paste0("s = exp(s_model * s_scaling_factor + s_anchor)\ns_anchor = ",round(s_anchor,3),", s_scaling_factor = ",round(s_scaling_factor,3),"\n"),paste0("s_model = ",c(-1,0,1),"  ->  s = ",exp(c(-1,0,1) * s_scaling_factor + s_anchor),"\n"))
model <- load_supermodel("suprnova", L=15, num_filters=5, s_anchor=s_anchor, s_scaling_factor=s_scaling_factor, custom_conv1=TRUE); gc()
model_s <- get_activation_model(model, c("d","s")); gc()
AFs_all <- array(0,dim=c(0,4)); expected_all <- array(0,dim=c(0,4)); s_vals_all <- array(0,dim=c(0,4)); d_all <- array(0,dim=c(0,4)); region_types_all <- array(0,dim=c(0,8))
dat_name = "transcribed100k_4"
for(dat_name in paste0("transcribed100k_",c(1,2))) {
    print(dat_name)
    load_supermodel_training_data(dat_name); dat_obs_exp[,1] <- rep(1,nrow(dat_obs_exp))
    refgene <- run_annovar(dat, "refGene", "hg38")
    is_canonical_splicing <- which(refgene$Func.refGene == "splicing")
    test_idx = c(is_canonical_splicing,1:7500)
    preds <- model_predict(model_s, test_idx)
    s_vals <- preds[["s"]]; s_vals[which(is.nan(s_vals_all))] <- 0
    s_vals_all <- abind(s_vals_all, s_vals, along=1)
    canon_splice <- rep(FALSE, length(test_idx)); if(length(is_canonical_splicing) > 0) { canon_splice[1:length(is_canonical_splicing)] <- TRUE }
    region_types_all <- abind(region_types_all, cbind(dat_region_types[test_idx,],canon_splice), along=1)
    d_all <- abind(d_all, preds[["d"]], along=1)
    AFs_all <- abind(AFs_all, af_tensor[test_idx,], along=1)
    expected_all <- abind(expected_all, expected[test_idx,], along=1)
}; region_types_all <- data.frame(region_types_all)
colnames(region_types_all) <- c("CDS","5'UTR","3'UTR","5'ss","3'ss","intron","intergenic","canon. splice")
region_types_all <- region_types_all[,which(!(colnames(region_types_all) %in% "intergenic"))]
# Plot region distributions for each score bucket.
a <- investigate_region_type_relationship(region_types_all, d_all, num_buckets=6, score_type="d", normalize_counts=FALSE, normalize_buckets=FALSE, tilted_labels=TRUE, class_label_cex=1.1, class_label_x_offset=3.5, srt=0, draw_legend=TRUE, legend_loc="topright", legend_ncol=1, scientific_notation=FALSE)
a <- investigate_region_type_relationship(region_types_all, d_all, num_buckets=6, score_type="d", normalize_counts=TRUE, normalize_buckets=FALSE, tilted_labels=TRUE, class_label_cex=1.1, class_label_x_offset=3.5, srt=0, draw_legend=TRUE, legend_loc="topright", legend_ncol=2, scientific_notation=FALSE)
a <- investigate_region_type_relationship(region_types_all, d_all, num_buckets=6, score_type="d", normalize_counts=FALSE, normalize_buckets=TRUE, tilted_labels=TRUE, class_label_cex=1.1, class_label_x_offset=3.5, srt=0, draw_legend=TRUE, legend_ncol=1, scientific_notation=FALSE)
a <- investigate_region_type_relationship(region_types_all, d_all, num_buckets=6, score_type="d", normalize_counts=TRUE, normalize_buckets=TRUE, tilted_labels=TRUE, class_label_cex=1.1, class_label_x_offset=3.5, srt=0, draw_legend=TRUE, legend_ncol=ceiling((ncol(dat_region_types)+1)/2), scientific_notation=FALSE)
a
# Plot score distributions for each region type.
investigate_score_vs_region_type_relationship(region_types_all, d_all, score_type="d", righttail_quantile=0.5)


range(preds[["conv_alt_ref"]])
dim(preds[["conv_alt_ref"]])
dim(preds[["rbp_binding_dense"]])
lapply(model$weights, function(x) return(c(x$name, paste0("Dim: c(",paste(dim(x),collapse=","),")"))))
model$layers[[23]]

conv1_weights <- as.array(model$get_layer("conv1")$weights[[1]][,,1,1,])
dimnames(conv1_weights) <- list(paste0("site",1:nrow(conv1_weights)), paste0("RBP",1:4), paste0("filter",1:dim(conv1_weights)[3]))
conv1_activations <- preds[["conv1"]]
rbp_complex_names <- rollapply(names(rbp_order), width=4, by=2, FUN=function(x) paste0(x, collapse=", "))
dimnames(conv1_activations) <- list(paste0("seq",1:nrow(conv1_activations)), paste0("pos",1:ncol(conv1_activations)), rbp_complex_names, c("A","C","G","T"), paste0("filter",1:dim(conv1_activations)[5]))
for(i in 1:(dim(conv1_weights)[3])) {
    filt <- t(conv1_weights[,,i]) #%>% layer_activation_leaky_relu(alpha=0.3))
    colnames(filt) <- NULL
    filename = output_path(paste0("suprnova_conv1_filter",i,".pdf"))
    pdf(filename)
    print(Heatmap(filt, 
                  show_heatmap_legend = TRUE, name = "weight", #title of legend
                  row_title = "", column_title = "Position",
                  cluster_rows=FALSE, cluster_columns=FALSE,
                  row_dend_side="left", column_dend_side="bottom",
                  row_names_side="left", column_names_side="top" #, column_names_rot=0#, #,
                  #row_names_gp = gpar(col=filter_cols, fontsize = 10), column_names_gp = gpar(fontsize = 20) # Text size for row and column names
    ))
    dev.off()
    pdf_to_png(filename)
}

conv1_weights_rbp_involvements <- t(apply(conv1_weights, 3, function(x) { apply(x, 2, max) }))    
rbp_complex_involvements <- abind(lapply(1:nrow(conv1_activations), function(seq_i) {
    abind(lapply(1:(dim(conv1_activations)[3]), function(rbp_complex_i) {
        conv1_activations[seq_i,1,rbp_complex_i,,] %*% conv1_weights_rbp_involvements
    }), along=3)
}), along=4)
rbp_complex_involvements <- aperm(rbp_complex_involvements, c(4,3,2,1)); rownames(rbp_complex_involvements) <- paste0("seq",1:nrow(conv1_activations)); colnames(rbp_complex_involvements) <- rbp_complex_names
dim(rbp_complex_involvements)
rbp_complex_involvements[1,1:3,,]

bases <- c("A","C","G","T")
ref_indices <- apply(af_tensor[test_idx,], 1, function(x) which(x < 0)); refs <- bases[ref_indices]
rbp_involvements[1,]
most_constrained_alt_indices <- unfactorize(data.frame(t(data.frame(apply(s_vals_all, 1, function(x) { max_i = which.max(x); return(c(max_i, bases[max_i], x[max_i])) }))))); colnames(most_constrained_alt_indices) <- c("alt_index", "alt", "s")
#paste0(refs,">",most_constrained_alt_indices$alt)

most_constrained_sites <- order(most_constrained_alt_indices$s, decreasing=TRUE)[1:10]
most_constrained_alt_indices[most_constrained_sites,]

site_i <- most_constrained_sites[2]
most_constrained_complex <- rbp_complex_involvements[site_i,,,most_constrained_alt_indices$alt_index[1]]
most_constrained_complex <- most_constrained_complex[order(rowSums(most_constrained_complex), decreasing=TRUE),]
most_constrained_complex[1:10,]

rbps <- sort(names(rbp_order)); midpoint = ceiling(ncol(dat_gradcams)/2)
most_constrained_gradcams <- dat_gradcams[site_i,midpoint,,most_constrained_alt_indices$alt_index[1],]
rownames(most_constrained_gradcams) <- rbps; colnames(most_constrained_gradcams) <- c("alt-ref delta", "ref_binding_score")
most_constrained_gradcams <- most_constrained_gradcams[order(abs(most_constrained_gradcams[,1]), decreasing=TRUE),]
most_constrained_gradcams[1:20,]

#range(preds[["gradcam_reformat"]]); preds[["gradcam_reformat"]][1,1:3,1:5,,]
range(preds[["conv_alt_ref"]]); preds[["conv_alt_ref"]][1,1:3,1:5,,]
range(preds[["rbp_binding_dense"]]); preds[["rbp_binding_dense"]][1,1:3,1:5,,] # plot(density(preds[["rbp_binding_dense"]]))
range(preds[["conv1"]]) # plot(density(preds[["conv1"]]))
range(preds[["d"]]) # plot(density(preds[["d"]]))
range(preds[["dense1"]]) # plot(density(preds[["dense1"]]))
range(preds[["dense2"]]) # plot(density(preds[["dense2"]]))

plot(density(log10((expected_all/s_vals_all)[AFs_all >= 0 & AFs_all < 1e-3])), col="red", main="AF vs. mu/s Comparison")
lines(density(log10(AFs_all[AFs_all >= 0 & AFs_all < 1e-3])), col="blue")
mtext(paste0("Spearman Corr = ",round(cor((expected_all/s_vals_all)[AFs_all > 0 & AFs_all < 1e-3], AFs_all[AFs_all > 0 & AFs_all < 1e-3], method="spearman"),5)))
legend("topleft", legend=c("AF","mu / s_pred"), col=c("blue","red"), pch=15, cex=0.7)
round(cor((expected_all/s_vals_all)[AFs_all > 0 & AFs_all < 1e-3], AFs_all[AFs_all > 0 & AFs_all < 1e-3], method="spearman"),5)
draw_plot(data.frame(x=log10(expected_all/s_vals_all)[AFs_all > 0 & AFs_all < 1e-3], y=log10(AFs_all)[AFs_all > 0 & AFs_all < 1e-3]), hex_density = 25, title="AF vs. mu/s", xlab="log10(mu/s)", ylab="log10(AF)", legend_text="# variants", linear_best_fit=TRUE, quadratic_best_fit=FALSE, filename="AF_vs_mu_over_s.pdf", cor_method="Spearman")


s_cutoff = 0.002
filename=output_path(paste0("model_distributions_test_set_s",s_cutoff,".pdf"))
pdf(file=filename)
AF_zero_value = 1e-6
AFs_all_backup <- AFs_all
AFs_all[AFs_all == 0] <- AF_zero_value
num_constrained = sum(AFs_all >=0 & AFs_all < 1e-3 & s_vals_all > s_cutoff)
plot(density(log10(expected_all[AFs_all >=0 & AFs_all < 1e-3 & s_vals_all > s_cutoff])), col="green", lty=3, main="SUPRNOVA Test Set Distributions", xlab="log10-scale", ylim=c(0,6), xlim=c(min(log10(expected[AFs_all >=0 & s_vals_all > 0.001])),max(AFs_all)), cex.main=1.3, cex.lab=1.4, cex.axis=1.4)
lines(density(sample(log10(expected_all[AFs_all >=0 & AFs_all < 1e-3 & s_vals_all <= s_cutoff]),num_constrained,replace=TRUE)), col="green", lty=1)
lines(density(log10(s_vals_all[AFs_all >=0 & AFs_all < 1e-3 & s_vals_all > s_cutoff])), col="red", lty=3)
lines(density(sample(log10(s_vals_all[AFs_all >=0 & AFs_all < 1e-3 & s_vals_all <= s_cutoff]),num_constrained,replace=TRUE)), col="red", lty=1)
lines(density(log10((expected_all/AFs_all)[AFs_all >=0 & AFs_all < 1e-3 & s_vals_all > s_cutoff])), col="orange", lty=3)
lines(density(sample(log10((expected_all/AFs_all)[AFs_all >=0 & AFs_all < 1e-3 & s_vals_all < s_cutoff]),num_constrained,replace=TRUE)), col="orange", lty=1)
lines(density(log10(AFs_all[AFs_all >=0 & AFs_all < 1e-3 & s_vals_all > s_cutoff]), from=log10(AF_zero_value), to=-3), col="blue", lty=3)
lines(density(sample(log10(AFs_all[AFs_all >=0 & AFs_all < 1e-3 & s_vals_all <= s_cutoff]),num_constrained,replace=TRUE), from=log10(AF_zero_value), to=-3), col="blue", lty=1)
mtext(paste0("Common variants (AF > ",1e-3,") excluded"), cex=1)
legend("topleft", legend=c(paste0("AF, Zerotons set to ",AF_zero_value),"mu","mu/AF","s",paste0("s > ",s_cutoff," (N=",num_constrained,")"),paste0("s <= ",s_cutoff," (N=",sum(AFs_all >=0 & AFs_all < 1e-3 & s_vals_all <= s_cutoff),")")), col=c("blue","green","orange","red","black","black"), pch=c(15,15,15,15,NA,NA), lty=c(NA,NA,NA,NA,3,1), cex=0.8)
dev.off()
pdf_to_png(filename)
quantile(log10((expected_all/AFs_all)[AFs_all >=0 & AFs_all < 1e-3 & s_vals_all > s_cutoff]), probs=c(0,0.05,0.1,0.2))
quantile(log10((expected_all/AFs_all)[AFs_all >=0 & AFs_all < 1e-3 & s_vals_all <= s_cutoff]), probs=c(0,0.05,0.1,0.2))
round(cor((expected_all/AFs_all)[AFs_all >=0 & AFs_all < 1e-3], s_vals_all[AFs_all >=0 & AFs_all < 1e-3], method="spearman"),5)
round(cor((expected_all/AFs_all)[AFs_all >=0 & AFs_all < 1e-3 & s_vals_all > s_cutoff], s_vals_all[AFs_all >=0 & AFs_all < 1e-3 & s_vals_all > s_cutoff], method="spearman"),5)
draw_plot(data.frame(x=log10((expected_all/AFs_all)[AFs_all >=0 & AFs_all < 1e-3]), y=log10(s_vals_all[AFs_all >=0 & AFs_all < 1e-3])), hex_density = 25, title="Predicted s vs. mu/AF", xlab=paste0("log10(mu/AF), Zerotons set to ",AF_zero_value," AF"), ylab="log10(s)", legend_text="# variants", linear_best_fit=TRUE, quadratic_best_fit=FALSE, filename="s_vs_mu_over_AF.pdf", cor_method="Spearman")
AFs_all <- AFs_all_backup

filename=output_path("s_distribution_test_set.pdf")
pdf(file=filename)
plot(density(log10(s_vals_all)), main="Selection Coef. Distribution in Test Data", xlab="log10(s)", col="blue", lwd=2, cex.axis=1.4, cex.lab=1.4, cex.main=1.3)
s_constrained_threshold = 0.005
mtext(paste0(sum(s_vals_all > s_constrained_threshold)," test variants have s > ",s_constrained_threshold), cex=1.2)
dev.off()
pdf_to_png(filename)

d <- preds[["d"]]; #d <- max(d) - d
range(preds[["rbp_binding_dense"]])
filename=output_path("d_distribution_test_set.pdf")
pdf(file=filename)
plot(density(d), main="Gene Reg. Damagingness Distribution in Test Data", xlab="d", col="blue", lwd=2, cex.axis=1.4, cex.lab=1.4, cex.main=1.3)
mtext(paste0("Representation of total RBP disruption, aggregated over 160 RBPs"), cex=1.2)
dev.off()
pdf_to_png(filename)

draw_plot(data.frame(x=c(d), y=c(s_vals_all)), hex_density = 25, title="Selection Coef. vs. Gene Reg. Damagingness", xlab="Gene Reg. Damagingness d", ylab="Selecion Coef. s", legend_text="# variants", linear_best_fit=FALSE, quadratic_best_fit=FALSE, filename="s_vs_d.pdf", cor_method="Spearman")

conv1 <- preds[["conv1"]]
cor(c(max(d)-d), c(s_vals_all), method="spearman")
plot(log10(s_vals_all)[AFs_all < 1e-3], d[AFs_all < 1e-3], main="d vs. s", xlab="s", ylab="d")
mtext(paste0("Spearman Corr = ",round(cor(log10(s_vals_all)[AFs_all < 1e-3], d[AFs_all < 1e-3], method="spearman"),5)))
dense1 <- preds[["dense1"]][,76,]
dense2 <- preds[["dense2"]][,76,]
dense3 <- preds[["dense3"]][,76,]
dense_final <- preds[["dense_final"]][,76,]

filename=output_path("s_distribution_test_set.pdf")
pdf(file=filename)
#xlim=c(log10(min(s_vals_all[AFs_all < 1e-3])),0), xaxs="i"
plot(density(log10((s_vals_all[AFs_all < 1e-3]))), main="Selection Coef. Distribution in Test Data", xlab="log10(s)", col="blue", lwd=2, cex.axis=1.4, cex.lab=1.4, cex.main=1.3)
s_constrained_threshold = 0.005
mtext(paste0("N = ",sum(s_vals_all[AFs_all < 1e-3] > s_constrained_threshold),", or ",round(100*sum(s_vals_all[AFs_all < 1e-3] > s_constrained_threshold)/length(c(s_vals_all[AFs_all < 1e-3])),3),"% of test set sites have s > ",s_constrained_threshold), cex=1.2)
dev.off()
pdf_to_png(filename)

plot(density(d))
range(d)
plot(density(s_vals_all))
range(s_vals_all)
plot(density(log10(s_vals_all)))
model$get_layer("conv1")$weights[[1]]
model$get_layer("conv1_relu")$weights[[1]]
as.matrix(model$get_layer("conv2")$weights[[1]][,,1,2])
as.matrix(model$get_layer("conv1_relu")$weights[[1]][,,1,2])
model$get_layer("dense1")$weights[[1]]
model$get_layer("dense2")$weights[[1]]
model$get_layer("per_site_sigmoids")$weights[[1]]
model$get_layer("pLI_activation")$weights[[1]]
range((expected_all/s_vals_all)[AFs_all < 1e-3]); range(AFs_all[AFs_all < 1e-3])
plot(density(log10((expected_all/s_vals_all)[AFs_all > 0 & AFs_all < 1e-3])), col="red", main="AF vs. mu/s Comparison")
lines(density(log10(AFs_all[AFs_all > 0 & AFs_all < 1e-3])), col="blue")
mtext(paste0("Spearman Corr = ",round(cor((expected_all/s_vals_all)[AFs_all < 1e-3], AFs_all[AFs_all < 1e-3], method="spearman"),5)))
legend("topleft", legend=c("AF","mu / s_pred"), col=c("blue","red"), pch=15, cex=0.7)
round(cor((expected_all/s_vals_all)[AFs_all < 1e-3 & AFs_all > 1e-4], AFs_all[AFs_all < 1e-3 & AFs_all > 1e-4], method="spearman"),5)
draw_plot(data.frame(x=log10(expected_all/s_vals_all)[AFs_all < 1e-3], y=log10(AFs_all)[AFs_all < 1e-3]), hex_density = 25, title="AF vs. mu/s", xlab="log10(mu/s)", ylab="log10(AF)", legend_text="# variants", linear_best_fit=TRUE, quadratic_best_fit=FALSE, filename="AF_vs_mu_over_s.pdf", cor_method="Spearman")

load_supermodel_training_data("transcribed100k_2")
test_indices2 <- train_supermodel(num_epochs=30, batch_size=256, rbp_order=rbp_order, fraction_for_training=0.5)
save_supermodel(model, "suprnova")
maps_result <- plot_model_outputs(model, test_indices[rowSums(af_tensor)[test_indices]>0][1:100])

suprnova <- load_supermodel("suprnova", L=15, num_filters=5, rbp_kernel_width=3, sequence_kernel_width=1, custom_conv1=TRUE, s_anchor=s_anchor, s_scaling_factor=s_scaling_factor)
model <- get_activation_model(suprnova, "s")
load_supermodel_training_data("transcribed100k_3")
trimers_all <- paste0(unlist(dat_trimers[,76])); rm(dat_trimers)
region_types <- dat_region_types; rm(dat_region_types)
plis <- dat_pLIs[,1]; rm(dat_pLIs)
AFs_all <- af_tensor[,76]; expected_all <- expected[,76]
s_vals_all <- (model %>% predict(list(dat_gradcams[,,rbp_order,,drop=FALSE], dat_obs_exp[,,drop=FALSE], expected[,,drop=FALSE]*mu_scaling_factor, af_tensor[,,drop=FALSE])))[,76]
s_vals_all[is.nan(s_vals_all)] <- 0
plot(density(log10(s_vals_all)))
gc(); load_supermodel_training_data("transcribed100k_4")
trimers_all <- c(trimers_all, paste0(unlist(dat_trimers[,76]))); rm(dat_trimers)
region_types <- rbind(region_types, dat_region_types); rm(dat_region_types)
plis <- c(plis, dat_pLIs[,1]); rm(dat_pLIs)
AFs_all <- c(AFs_all, af_tensor[,76]); expected_all <- c(expected_all, expected[,76])
s_vals_all <- c(s_vals_all, (model %>% predict(list(dat_gradcams[,,rbp_order,,drop=FALSE], dat_obs_exp[,,drop=FALSE], expected[,,drop=FALSE]*mu_scaling_factor, af_tensor[,,drop=FALSE])))[,76])
s_vals_all[is.nan(s_vals_all)] <- 0

plot(density(log10(s_vals_all)))
plot(density(expected_all / AFs_all))
plot(log10(expected_all/s_vals_all)[AFs_all < 1e-3], log10(AFs_all)[AFs_all < 1e-3])
mtext(paste0("Spearman Corr = ",round(cor((expected_all/s_vals_all)[AFs_all < 1e-3], AFs_all[AFs_all < 1e-3], method="spearman"),2)))
draw_plot(data.frame(x=log10(expected_all/s_vals_all)[AFs_all < 1e-3], y=log10(AFs_all)[AFs_all < 1e-3]), hex_density = 25, title="AF vs. mu/s", xlab="log10(mu/s)", ylab="log10(AF)", legend_text="# variants", linear_best_fit=TRUE, quadratic_best_fit=FALSE, filename="AF_vs_mu_over_s.pdf", cor_method="Spearman")
#draw_plot(data.frame(x=log10(expected_all/s_vals_all)[AFs_all < 1e-3], y=log10(AFs_all)[AFs_all < 1e-3]), hex_density = 10, title="AF vs. mu/s", xlab="log10(mu/s)", ylab="log10(AF)", legend_text="# variants", linear_best_fit=TRUE, quadratic_best_fit=FALSE, filename="AF_vs_mu_over_s.pdf", cor_method="Pearson")

plot(log10(expected_all), log10(AFs_all))
mtext(paste0("Spearman Corr = ",round(cor(expected_all, AFs_all, method="spearman"),2)))
plot(log10(expected_all/s_vals_all)[s_vals_all>0.01], log10(AFs_all)[s_vals_all>0.01])
mtext(paste0("Spearman Corr = ",round(cor((expected_all/s_vals_all)[s_vals_all>0.02], AFs_all[s_vals_all>0.02], method="spearman"),2)))
draw_plot(data.frame(x=log10(expected_all/s_vals_all), y=log10(AFs_all)), hex_density = 10, title="AF vs. mu/s", xlab="log10(mu/s)", ylab="log10(AF)", legend_text="# variants", linear_best_fit=TRUE, quadratic_best_fit=FALSE, filename="AF_vs_mu_over_s.pdf", cor_method="Spearman")
draw_plot(vdata.frame(x=log10(expected_all/s_vals_all)[s_vals_all>0.02], y=log10(AFs_all)[s_vals_all>0.02]), hex_density = 10, title="AF vs. mu/s", xlab="log10(mu/s)", ylab="log10(AF)", legend_text="# variants", linear_best_fit=TRUE, quadratic_best_fit=FALSE, filename="AF_vs_mu_over_s.pdf", cor_method="Spearman")
AFs_all[s_vals_all > 0.1]
expected_all[s_vals_all > 0.1]
expected_all[s_vals_all > 0.1] / s_vals_all[s_vals_all > 0.1]
plot(log10(expected_all/s_vals_all), log10(AFs_all))
mtext(paste0("Spearman Corr = ",round(cor(expected_all/s_vals_all, AFs_all, method="spearman"),2)))

cpg_all <- grepl("CG",trimers_all)

s_bucket_cuts <- rbind(c(0,0,0,0,0,0,0,0,0,0,0.003,0.004,0.005,0.006),c(0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.01,0.02,0.02,0.02,0.02,0.02))
#s_bucket_cuts <- rbind(c(0,0,0,0,0,0,0,0,0,0,0.005),c(0.001,0.0015,0.002,0.0025,0.003,0.0035,0.004,0.0045,0.005,0.01,0.01))
s_bucket_names <- c(paste0(s_bucket_cuts[1,]," < s <= ",s_bucket_cuts[2,]), "s = any")
library("fitdistrplus")
for(cpg in c(FALSE,TRUE)) {
    if(cpg) { cpg_suffix = "_CpG" } else { cpg_suffix = "_non-CpG" }
    haha <- Reduce(cbind, lapply(1:(ncol(s_bucket_cuts)+1), function(s_i) {
        print(s_i)
        if(s_i == 0) { curr_indices_backup <- s_vals_all == 0 } else if(s_i > ncol(s_bucket_cuts)) { curr_indices_backup <- rep(TRUE, length(s_vals_all)) } else { curr_indices_backup <- s_vals_all > s_bucket_cuts[1,s_i] & s_vals_all <= s_bucket_cuts[2,s_i]  }
        curr_indices <- curr_indices_backup & AFs_all < 1e-4
        if(cpg) { curr_indices <- curr_indices & cpg_all } else { curr_indices <- curr_indices & !cpg_all }
        x <- round(AFs_all[curr_indices] * sample_size + 0.01)
        x[is.nan(x)] <- 0
        if(length(x) > 1) { f <- fitdist(x, "nbinom"); f <- c(rev(f$estimate), rev(f$sd), f$loglik, f$aic, f$bic)
        } else { f <- c(0,0,0,0,0,0,0) }
        x <- c(sum(x==0),sum(x==1),sum(x==2),sum(x==3),sum(x==4),sum(x==5),sum(x>5),f)
        x <- x
        return(x)
    }))
    colnames(haha) <- s_bucket_names; rownames(haha) <- c("Zerotons", "Singletons", "Doubletons", "Tripletons", "AC = 4", "AC = 5", "AC > 5", "NB_mu", "NB_size", "NB_mu_sd", "NB_size_sd", "NB_LL", "NB_AIC", "NB_BIC")
    haha
    haha["NB_mu",] #colSums(haha[1:7,])
    haha["NB_size",]
    write.csv(haha, output_path(paste0("s_NB_fits",cpg_suffix,".csv")))
    
    filename=paste0(output_path("s_NB_comparison"),cpg_suffix)
    pdf(file=filename)
    fitD <- dnbinom(0:10, size= haha["NB_size",13], mu=haha["NB_mu",13]) #rnbinom(100000, size= haha["NB_size",6], mu=haha["NB_mu",6])
    plot(0:10, fitD, type="o", lwd=2, col="red", main="Comparison of Negative Binomial Fits", xlab="Allele Count (AC)", ylab="Density", cex.lab=1.4, cex.axis=1.4, cex.main=1.3)
    lines(0:10, dnbinom(0:10, size= haha["NB_size",4], mu=haha["NB_mu",4]), lwd=2, col="blue", type="o")
    legend("topright", legend=c(colnames(haha)[c(4,13)],paste0("NB(mu=",round(haha["NB_mu",c(4,13)],3),", size=",round(haha["NB_size",c(4,13)],3),")"))[c(1,3,2,4)], col=c("blue","white","red","white"), pch=15, cex=1.1)
    mtext(paste0(gsub("^_","",cpg_suffix)," only"), cex=1.2)
    dev.off()
    pdf_to_png(filename)
    
    haha[-c(1:7),]
    table(s_vals_all[s_vals_all > 0.006])
    apply(haha, 2, function(x) sum(x[1:2])/sum(x))
    apply(haha, 2, function(x) sum(x[2])/sum(x[-1]))
    apply(haha, 2, function(x) sum(x[1])/sum(x))
    
    haha <- haha[1:7,]
    for(leg in c(TRUE, FALSE)) {
        filename=paste0(output_path("s_sfs_comparison"),cpg_suffix)
        if(leg) { filename = paste0(filename, "_legend.pdf") } else { filename = paste0(filename, ".pdf") }
        pdf(file=filename)
        cols <- rainbow(nrow(haha))
        colnames(haha) <- NULL
        barplot(apply(haha[,-ncol(haha)], 2, function(x) x/sum(x)), beside=TRUE, las=2, main="Site Frequency Spectrum vs. Selection Coef.", ylab="Proportion (%)", col=cols) #, ylim = c(0,5 + max(mtcars$qsec)), xlab = "", space = 1)
        end_point = 0.5 + ncol(haha)-1 + ncol(haha)-1 - 1 #this is the line which does the trick (together with barplot "space = 1" parameter)
        #rotate 60 degrees (srt = 60)
        text(seq(1.5+5, end_point*4+5, by = 2*4), par("usr")[3]-0.03, srt=60, adj=1, xpd=TRUE, labels=s_bucket_names[-ncol(haha)], cex=0.65)
        if(leg) { legend("topleft", legend=c(rownames(haha)), col=cols, pch=15, cex=0.75) }
        mtext(paste0(gsub("^_","",cpg_suffix)," only"), cex=1.2)
        dev.off()
        pdf_to_png(filename)
    }
    colnames(haha) <- s_bucket_names
    write.csv(haha, output_path(paste0("s_sfs_comparison_counts",cpg_suffix,".csv")))
}

# Investigate relationship of gene constraint and s.
s_vals <- s_vals_all[(length(s_vals_all)-length(pLIs_dat)+1):length(s_vals_all)]
draw_plot(data.frame(x=plis, y=log10(s_vals_all)), hex_density = 10, title="s vs. Gene Constraint", xlab="pLI", ylab="log10(s)", legend_text="# variants", linear_best_fit=TRUE, quadratic_best_fit=FALSE, filename="s_vs_gene_constraint.pdf", cor_method="Spearman")
draw_plot(data.frame(x=pLIs_dat, y=log10(s_vals_all)[(length(s_vals_all)-length(pLIs_dat)+1):length(s_vals_all)]), hex_density = 10, title="s vs. Gene Constraint", xlab="pLI", ylab="log10(s)", legend_text="# variants", linear_best_fit=TRUE, quadratic_best_fit=FALSE, filename="s_vs_gene_constraint.pdf", cor_method="Spearman")
draw_plot(data.frame(x=pLIs_dat[log10(s_vals) > -2.4], y=log10(s_vals)[log10(s_vals) > -2.4]), hex_density = 10, title="AF vs. mu/s", xlab="pLIs_dat", ylab="log10(s)", legend_text="# variants", linear_best_fit=TRUE, quadratic_best_fit=FALSE, filename="s_vs_gene_constraint.pdf", cor_method="Spearman")

pli_threshold = 0.9
b <- cbind(mean(plis[AFs_all < 1e-3] > pli_threshold), mean(plis[AFs_all < 1e-3 & s_vals_all > 0.004] > pli_threshold), mean(plis[AFs_all < 1e-3 & s_vals_all > 0.005] > pli_threshold))
rownames(b) <- c(paste0("pLI > ",pli_threshold)); colnames(b) <- c("All variants", "s > 0.004", "s > 0.005")
b
filename=paste0(output_path("s_pLI_comparison.pdf"))
pdf(filename)
cols <- "gray"
barplot(b, main="Site Frequency Spectrum vs. Selection Coef.", ylab="Proportion (%)", col=cols) #, ylim = c(0,5 + max(mtcars$qsec)), xlab = "", space = 1)
legend("topright", legend=c(rownames(b)), col=cols, pch=15, cex=0.75)
dev.off()
pdf_to_png(filename)

dat_region_types <- dat_region_types[,which(!(colnames(dat_region_types) %in% "intergenic"))]
colnames(dat_region_types) <- c("CDS","5'UTR","3'UTR","5'ss","3'ss","intron")
a <- investigate_region_type_relationship(dat_region_types[test_idx,], preds[["d"]], num_buckets=6, score_type="d", normalize_counts=FALSE, normalize_buckets=FALSE, tilted_labels=TRUE, class_label_cex=1.1, class_label_x_offset=3.5, srt=0, draw_legend=TRUE, legend_ncol=1, scientific_notation=FALSE)
a <- investigate_region_type_relationship(dat_region_types[test_idx,], preds[["d"]], num_buckets=6, score_type="d", normalize_counts=TRUE, normalize_buckets=FALSE, tilted_labels=TRUE, class_label_cex=1.1, class_label_x_offset=3.5, srt=0, draw_legend=TRUE, legend_loc="topright", legend_ncol=ceiling((ncol(dat_region_types)+1)/2), scientific_notation=FALSE)
a <- investigate_region_type_relationship(dat_region_types[test_idx,], preds[["d"]], num_buckets=6, score_type="d", normalize_counts=FALSE, normalize_buckets=TRUE, tilted_labels=TRUE, class_label_cex=1.1, class_label_x_offset=3.5, srt=0, draw_legend=TRUE, legend_ncol=1, scientific_notation=FALSE)
a <- investigate_region_type_relationship(dat_region_types[test_idx,], preds[["d"]], num_buckets=6, score_type="d", normalize_counts=TRUE, normalize_buckets=TRUE, tilted_labels=TRUE, class_label_cex=1.1, class_label_x_offset=3.5, srt=0, draw_legend=TRUE, legend_ncol=ceiling((ncol(dat_region_types)+1)/2), scientific_notation=FALSE)
a
investigate_score_vs_region_type_relationship(dat_region_types[test_idx,], preds[["d"]], score_type="d", righttail_quantile=0.5)

refgene <- run_annovar(dat, "refGene", "hg38")
is_canonical_splicing <- which(refgene$Func.refGene == "splicing")


region_classes <- list(is_cds,is_cds_syn,is_cds_startloss,is_cds_missense,is_cds_stoploss,is_cds_nonsense,is_5utr,is_3utr,is_intron,is_5ss,is_3ss,is_intergenic)
region_class_names <- c("All CDS","Synonymous","Start lost","Missense","Stop lost","Nonsense","5'UTR","3'UTR","intron","5'ss","3'ss","intergenic")
num_bootstrap_samples = 10
b <- unfactorize(data.frame(rbindlist(lapply(1:length(region_classes), function(class_indices_i) {
    print(paste0(class_indices_i,". ",region_class_names[class_indices_i]))
    class_indices <- region_classes[[class_indices_i]]
    AC <- a_AC[class_indices]; trimers <- a_trimers[class_indices]; N = sum(class_indices)
    x_sampled <- t(sapply(1:num_bootstrap_samples, function(sample_i) {
        if(N > 100000000) { gc() }
        print(sample_i)
        sample_indices <- sample(1:N, replace=TRUE)
        ps_obs <- sum(AC[sample_indices] == 1) / sum(AC[sample_indices] > 0)
        ps_exp <- mean(na.omit(trimers_expected_ps[trimers[sample_indices]]))
        rm(sample_indices)
        return(c(ps_obs, ps_exp))
    })); x_sampled <- cbind(x_sampled, x_sampled[,1] - x_sampled[,2]); colnames(x_sampled) <- c("PS_obs", "PS_exp", "MAPS")
    return(data.frame(t(data.frame(c(region_class_names[class_indices_i], N, sum(grepl("CG",trimers))/N, colMeans(x_sampled), x_sampled[max(c(1,floor(0.025*num_bootstrap_samples))),], x_sampled[min(c(ceiling(0.975*num_bootstrap_samples),num_bootstrap_samples)),])))))
}))))
colnames(b) <- c("region_type","N","CpG_freq","PS_obs","PS_exp","MAPS","PS_obs_lower","PS_exp_lower","MAPS_lower","PS_obs_upper","PS_exp_upper","MAPS_upper")
b
saveRDS(b, output_path("b_100k.rds"))

filename = output_path(paste0("regional_s_comparison.pdf"))
pdf(filename)
cols = c("black","black","white",rainbow(nrow(b)-1),"gray40")
plot(1:length(b$region_type), rev(b$MAPS-MAPS_offset), pch=16, xlab="", ylab="MAPS", main="Regional MAPS Distributions", ylim=range(c(b[,c("MAPS_lower","MAPS_upper")]))-MAPS_offset, col=rev(cols[-c(1:3)]), xaxt="n", cex.axis=1.3, cex.lab=1.3, cex.main=1.3)
segments(x0=1:length(b$region_type), y0=rev(b$MAPS_lower)-MAPS_offset, x1=1:length(b$region_type), y1=rev(b$MAPS_upper)-MAPS_offset, col=rev(cols[-c(1:3)]))
text(1:length(b$region_type)+0.2, par("usr")[3]-0.001, srt=45, adj=1, xpd=TRUE, labels=rev(b$region_type), cex=0.8)
if(num_gnomAD_variants_to_sample < nrow(dat_trimers)) { mtext(paste0(num_gnomAD_variants_to_sample," random genomic variants"), cex=1.2)
} else { mtext(paste0(num_gnomAD_variants_to_sample," gnomAD 3.0 variants"), cex=1.2) }
dev.off()
pdf_to_png(filename)


filename = output_path(paste0("b.pdf"))
pdf(filename)
cols = c("black","black","white",rainbow(nrow(b)-1),"gray40")
plot(b$CpG_freq, b$PS_obs, pch=16, xlab="CpG trimer frequency", ylim=range(b$PS_exp+MAPS_offset)*c(0.96,1), ylab="Proportion Singletons", main="Regional PS Breakdown", col=cols[-c(1:3)], cex.lab=1.4, cex.axis=1.4, cex.main=1.3)
points(b$CpG_freq, b$PS_exp+MAPS_offset, pch=8, col=cols[-c(1:3)])
legend("bottomleft", legend=c("Observed PS in gnomAD","Expected PS from trimers","",paste0(b$region_type," (MAPS = ",round(b$MAPS-MAPS_offset,4),")")), col=cols, pch=c(16,8,8,rep(15,nrow(b))), cex=1.1)
if(num_gnomAD_variants_to_sample < nrow(dat_trimers)) { mtext(paste0(num_gnomAD_variants_to_sample," random genomic variants"), cex=1.2)
} else { mtext(paste0(num_gnomAD_variants_to_sample," gnomAD 3.0 variants"), cex=1.2) }
dev.off()
pdf_to_png(filename)

filename = output_path(paste0("regional_MAPS_comparison_gnomAD.pdf"))
pdf(filename)
cols = c("black","black","white",rainbow(nrow(b)-1),"gray40")
plot(1:length(b$region_type), rev(b$MAPS-MAPS_offset), pch=16, xlab="", ylab="MAPS", main="Regional MAPS Distributions", ylim=range(c(b[,c("MAPS_lower","MAPS_upper")]))-MAPS_offset, col=rev(cols[-c(1:3)]), xaxt="n", cex.axis=1.3, cex.lab=1.3, cex.main=1.3)
segments(x0=1:length(b$region_type), y0=rev(b$MAPS_lower)-MAPS_offset, x1=1:length(b$region_type), y1=rev(b$MAPS_upper)-MAPS_offset, col=rev(cols[-c(1:3)]))
text(1:length(b$region_type)+0.2, par("usr")[3]-0.001, srt=45, adj=1, xpd=TRUE, labels=rev(b$region_type), cex=0.8)
if(num_gnomAD_variants_to_sample < nrow(dat_trimers)) { mtext(paste0(num_gnomAD_variants_to_sample," random genomic variants"), cex=1.2)
} else { mtext(paste0(num_gnomAD_variants_to_sample," gnomAD 3.0 variants"), cex=1.2) }
dev.off()
pdf_to_png(filename)

maps_sfs <- t(haha[,-ncol(haha)])
filename = output_path(paste0("s_sfs_logscale.pdf"))
pdf(file=filename)
s_cols <- rainbow(nrow(maps_sfs))
#c(0,ceiling(max(maps_sfs_full$PS_observed*10))*10)
plot(log10(1:(ncol(maps_sfs)-1)), log10(maps_sfs[1,-1]/sum(maps_sfs[1,-1])), col=s_cols[1], type="o", main="Site Frequency Spectrum", xlab="", ylab="log10(Proportion freq)", xaxt="n", xlim=log10(c(1,ncol(maps_sfs))), cex.axis=1.3, cex.lab=1.3, cex.main=1.3)
for(maps_i in 1:(ncol(maps_sfs)-1)) { abline(v=log10(maps_i), lty=3, col="gray50") }
for(s_i in 2:nrow(maps_sfs)) { lines(log10(1:(ncol(maps_sfs)-1)), log10(maps_sfs[s_i,-1]/sum(maps_sfs[s_i,-1])), col=s_cols[s_i], type="o") }
#if(sites == "CpG") { mtext_text = paste0(sum(maps_cpg_vector)," supermodel test sequence CpG sites")
#} else if (sites == "non-CpG") { mtext_text = paste0(sum(!maps_cpg_vector)," supermodel test sequence non-CpG sites")
#} else { mtext_text = paste0("all"," supermodel test sequence sites") }
mtext_text = paste0(sites," sites")
mtext(mtext_text, cex=1.1)
legend("topright", legend=c("Predicted selection coef. s", paste0(rownames(maps_sfs)," (N = ",rowSums(maps_sfs[,1:6]),")")), col=c("white",s_cols), pch=15, cex=1)
text(log10(1:ncol(maps_sfs))+0.04, par("usr")[3]-0.01, srt=45, adj=1, xpd=TRUE, labels=colnames(maps_sfs), cex=1.3)
dev.off()
pdf_to_png(filename)
apply(haha[,-ncol(haha)], 2, function(x) x[2]/sum(x[-1]))


#################################################################################################################
# Analyis of background mutation rate and allele frequency on the gradient of the Negative Binomial (NB) Loss with respect to s.  
#################################################################################################################
analyze_NB_parameter_effects <- function() {
    AFs_all <- c(af_tensor[rowSums(af_tensor)>0]); expected_all <- c(expected[rowSums(af_tensor)>0])
    filename = output_path(paste0("suprnova_loss_function.pdf"))
    pdf(file=filename, width=10, height=10)
    mus <- c(quantile(expected_all)[3], quantile(expected_all)[5], quantile(expected_all)[3], quantile(expected_all)[5], quantile(expected_all)[3], quantile(expected_all)[5], quantile(expected_all)[3], quantile(expected_all)[5])
    y_trues <- c(min(AFs_all), min(AFs_all), quantile(AFs_all)[2], quantile(AFs_all)[2], 2*quantile(AFs_all)[2], 2*quantile(AFs_all)[2], quantile(AFs_all)[5], quantile(AFs_all)[5])
    ss <- c(10**seq(-8,-4),0.001,0.005,0.01,0.015,0.02,0.025,0.03,0.04,0.05,0.1,0.5,1)
    par(mfrow=c(3,3))
    for(i in 1:(prod(par()$mfrow))) {
        if(i < min(c(length(mus)+1, length(y_trues)+1, prod(par()$mfrow)))) {
            #print(i)
            mu = mus[i]; y_true <- y_trues[i]
            probs <- sapply(ss, function(s) {
                z = 4 * Ne * s / sample_size # Since s = 1/s_real #/ sample_size  #0.04 / 71702
                r = 4 * Ne * mu # mu must be >= 0
                probs <- (round(y_true*sample_size+0.01) %>% (tfp$distributions$NegativeBinomial(total_count=r, probs=(1/(z+1)))$prob))
                neg_log_probs <- - (round(y_true*sample_size+0.01) %>% (tfp$distributions$NegativeBinomial(total_count=r, probs=(1/(z+1)))$log_prob))
                return(c(as.matrix(probs)[,1], as.matrix(neg_log_probs)[,1]))
            })
            #print(probs)
            if(i == floor((par()$mfrow[2]+1)/2)) { main_text = "Objective: Minimize -log(prob)" } else { main_text = "" }
            if(i %% (par()$mfrow[2]) == 1 || par()$mfrow[2] == 1) { ylab_text = "NB prob (%)" } else { ylab_text = "" }
            xlab_text = "log10(s)"
            plot(log10(ss), formatC(probs[1,]*100,format="e",digits=4), col="blue", type="o", main=main_text, xlab=xlab_text, ylab=ylab_text, cex.main=1.4, cex.lab=1.2)
            par(new = TRUE)
            plot(log10(ss), probs[2,], col="red", type = "o", xaxt = "n", yaxt = "n", ylab = "", xlab = "")
            par(new = FALSE)
            axis(side = 4, cex=0.25)
            #mtext("-log(prob)", side = 4, line = 3)
            mtext(paste0("AF = ",formatC(y_true,format="e",digits=2),", mu = ",formatC(mu,format="e",digits=2),""), cex=0.8)
        } else {
            plot.new()
            if (i == prod(par()$mfrow)) { 
                legend("center", legend=c("NB prob","NB -log(prob)"), col = c("blue","red"), pch=15, cex=1.5) 
                mtext(paste0(sample_size," gnomAD 3.0 samples"), cex=0.9)
            }
        }
    }
    dev.off()
    pdf_to_png(filename)
    par(mfrow=c(1,1))
    z = 4 * Ne * s_vals_all / sample_size # Since s = 1/s_real #/ sample_size  #0.04 / 71702
    r = 4 * Ne * expected_all
    probs <- as.matrix(round(AFs_all*sample_size+0.01) %>% (tfp$distributions$NegativeBinomial(total_count=r, probs=(1/(z+1)))$prob))[,1]
    mean(probs)
}

#################################################################################################################
# Relationship of max s nearest each gene with the gene's s_het: latter should be upper bound of the former, with some semi-linear relationship.
#################################################################################################################
a <- read.csv(data_path("genes_shet_cassa_sunaev.csv")); a$s_het[is.na(a$s_het)] <- mean(a$s_het[!is.na(a$s_het)])
plot(density(a$s_het), main="Gene s_het Distribution", xlab="s_het", cex.main=1.3, cex.lab=1.4, cex.axis=1.4)
mtext("Gene s_het values calculated from obs/exp ExAC PTVs", cex=1.1)
plot(density(log10(a$s_het)), col="blue", lwd=2, main="All Genes s_het Distribution", xlab="log10(s_het)", cex.main=1.3, cex.lab=1.4, cex.axis=1.4)
mtext("Gene s_het values calculated from obs/exp ExAC PTVs", cex=1.1)

genes_shet <- a$s_het; names(genes_shet) <- a$gene_symbol
# MAKE A GENERAL ANNOTATE-GENES(a$genes) FUNCTION TO EASILY ANNOTATE PLI, OBS/EXP, S_HET, etc. USING JUST ONE LINE, AND RUN IT ON GENE NAMES.
# STOPPPED HERE!!!
genebody_hg19 <- read.csv(data_path("refGene\\refGene_hg19_genebody_fixed.bed"), sep="\t", stringsAsFactors=FALSE)
genebody_granges_hg19 <- to_genomic_regions(genebody_hg19, chr_colname="chromosome", start_colname="tss", end_colname="tes", strand_colname="strand", label_colname="gene", order_coordinates=TRUE, remove_duplicate_labels=TRUE)

nearest_genes_all <- Reduce(c,lapply(paste0("transcribed100k_",1:4), function(name) {
    print(name)
    dat_granges_hg19 <- readRDS(output_path(paste0(name,"_granges_hg19.rds")))
    dat_nearest_genes <- names(genebody_granges_hg19)[nearest(dat_granges_hg19, genebody_granges_hg19)]
    saveRDS(dat_nearest_genes, output_path(paste0(name,"_nearest_genes.rds")))
    shet <- genes_shet[dat_nearest_genes]; shet[is.na(shet)] <- mean(genes_shet)
    saveRDS(data.frame(shet), output_path(paste0(name,"_shet.rds")))
    dat_eclip_overlaps <- find_eclip_overlaps(dat_granges_hg19)
    saveRDS(data.frame(dat_eclip_overlaps), output_path(paste0(name,"_eclip_overlaps.rds")))
    
    
    expected <- readRDS(output_path(paste0(name,"_expected.rds")))
    to_keep <- apply(expected, 1, function(x) sum(x==0) == 0)
    return(dat_nearest_genes[to_keep])
}))
shet_all <- genes_shet[nearest_genes_all]; shet_all[is.na(shet_all)] <- mean(genes_shet)
saveRDS(nearest_genes_all, output_path("transcribed100k_nearest_genes.rds"))
saveRDS(shet_all, output_path("transcribed100k_shet.rds"))

k_reset_freq = 20
rbps <- get_features_by_group("RBP")
for(name in paste0("transcribed100k_",1:4)) {
    print(name)
    dat_tensor <- readRDS(output_path(paste0(name,"_tensor.rds")))
    
    rbp_binding_scores <- Reduce(cbind, lapply(1:length(rbps), function(rbp_i) {
        if(rbp_i %% k_reset_freq == 0) { k_clear_session(); gc() }
        rbp = rbps[rbp_i]
        print(paste0(rbp_i,". ",rbp))
        model <- load_model_hdf5(paste0("../ML/output/old/",tolower(rbp),"_model2.h5"))
        return(model %>% predict(dat_tensor[,,,1]))
    }))
    
    saveRDS(rbp_binding_scores, output_path(paste0(name,"_binding_scores.rds")))
}


filename = output_path("shet_distribution.pdf")
pdf(filename)
plot(density(log10(shet_all)), col="red", lwd=2, main="Gene s_het Distribution", xlab="log10(s_het)", cex.main=1.3, cex.lab=1.4, cex.axis=1.4)
lines(density(log10(a$s_het)), col="blue", lwd=2, main="Gene s_het Distribution in Data", xlab="log10(s_het)", cex.main=1.3, cex.lab=1.4, cex.axis=1.4)
mtext("Gene s_het values calculated from obs/exp ExAC PTVs", cex=1.1)
legend("topleft", legend=c("All genes", "Nearest genes to variants"), col=c("blue","red"), pch=15, cex=1.2)
dev.off()
pdf_to_png(filename)

max_s_gene <- sapply(unique(nearest_genes_all), function(nearest_gene) {
    return(max(s_vals_all[which(nearest_genes_all == nearest_gene)]))
})
shet_gene <- genes_shet[unique(nearest_genes_all)]; shet_gene[is.na(shet_gene)] <- mean(genes_shet)

num_bootstrap_samples = 1000
shet_bucket_cuts <- rbind(c(0,0,0,0,0,0,0,0,0,0,0.003,0.004,0.005,0.006),c(0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.01,0.02,0.02,0.02,0.02,0.02))
#s_bucket_cuts <- rbind(c(0,0,0,0,0,0,0,0,0,0,0.005),c(0.001,0.0015,0.002,0.0025,0.003,0.0035,0.004,0.0045,0.005,0.01,0.01))
shet_bucket_names <- c(paste0(shet_bucket_cuts[1,]," < s_het <= ",shet_bucket_cuts[2,]), "s_het = any")
shet_analysis <- Reduce(cbind, lapply(1:(ncol(shet_bucket_cuts)+1), function(s_i) {
    print(s_i)
    if(s_i == 0) { curr_indices_backup <- shet_all == 0 } else if(s_i > ncol(shet_bucket_cuts)) { curr_indices_backup <- rep(TRUE, length(shet_all)) } else { curr_indices_backup <- shet_all > shet_bucket_cuts[1,s_i] & shet_all <= shet_bucket_cuts[2,s_i]  }
    curr_indices <- curr_indices_backup
    x <- s_vals_all[curr_indices]
    s_sampled <- sort(sapply(1:num_bootstrap_samples, function(sample_i) { 
        return(sample(x, length(x), replace=TRUE))
    }), decreasing=FALSE)
    s <- mean(x)
    return(data.frame(t(data.frame(c(s,s_sampled[round(num_bootstrap_samples*0.025)],s_sampled[round(num_bootstrap_samples*0.975)],length(x))[-1]))))
}))
colnames(shet_analysis) <- c("s", "s_lower", "s_upper")#,"N_genes")
rownames(shet_analysis) <- shet_bucket_names


filename = output_path("max_s_vs_shet_distribution.pdf")
pdf(filename)
plot(density(log10(max_s_gene)), col="red", lwd=2, main="Gene s_het Distribution", xlab="log10(selection coef.)", cex.main=1.3, cex.lab=1.4, cex.axis=1.4)
lines(density(log10(shet_all)), col="blue", lwd=2, main="Gene s_het Distribution in Data", xlab="log10(s_het)", cex.main=1.3, cex.lab=1.4, cex.axis=1.4)
mtext("Gene s_het values calculated from obs/exp ExAC PTVs", cex=1.1)
legend("topleft", legend=c("gene s_het", "max(variant s)"), col=c("blue","red"), pch=15, cex=1.2)
dev.off()
pdf_to_png(filename)

draw_plot(data.frame(x=max_s_gene, y=shet_gene), hex_density = 25, title="Gene s_het vs. s", xlab="Variant s", ylab="Gene s_het", legend_text="# variants", linear_best_fit=TRUE, quadratic_best_fit=FALSE, filename="AF_vs_mu_over_s.pdf", cor_method="Spearman")
draw_plot(data.frame(x=log10(max_s_gene), y=log10(shet_gene)), hex_density = 25, title="Gene s_het vs. s", xlab="Variant s", ylab="Gene s_het", legend_text="# variants", linear_best_fit=TRUE, quadratic_best_fit=FALSE, filename="AF_vs_mu_over_s.pdf", cor_method="Spearman")
draw_plot(data.frame(y=log10(max_s_gene), x=log10(shet_gene)), hex_density = 25, title="Max s vs. Gene s_het", xlab="Gene log10(s_het)", ylab="Max variant log10(s)", legend_text="# variants", linear_best_fit=FALSE, quadratic_best_fit=FALSE, filename="AF_vs_mu_over_s.pdf", cor_method="Spearman")
draw_plot(data.frame(y=log10(max_s_gene)[log10(max_s_gene) > -2.35], x=log10(shet_gene)[log10(max_s_gene) > -2.35]), hex_density = 25, title="Max s vs. Gene s_het", xlab="Gene log10(s_het)", ylab="Max variant log10(s)", legend_text="# variants", linear_best_fit=TRUE, quadratic_best_fit=FALSE, identity_diagonal=TRUE, filename="AF_vs_mu_over_s.pdf", cor_method="Spearman")

shet_gene[log10(max_s_gene) > -2.35]
unique(nearest_genes_all)[log10(max_s_gene) > -2.35]

name = paste0("transcribed100k_",1)
dat <- readRDS(output_path(paste0(name,"_dat.rds")))
dat_eclip_overlaps <- readRDS(output_path(paste0(name,"_eclip_overlaps.rds")))

for(name in paste0("transcribed100k_",1:4)) {
    print(name)
    dat_sequences <- readRDS(output_path(paste0(name,"_dat.rds")))$ref_sequence  #[to_keep,site_indices,,,drop=FALSE]
    dat_gradcams <- readRDS(output_path(paste0(name,"_gradcams.rds")))
    rbp_binding_scores <- readRDS(output_path(paste0(name,"_binding_scores.rds")))
    dat_eclip_overlaps <- readRDS(output_path(paste0(name,"_eclip_overlaps.rds")))
    strong_binding_cutoff = 0.85
    idx = which(apply(dat_eclip_overlaps, 1, function(x) { return(rbp_binding_scores[as.numeric(x[1]),which(rbps == x[3])]) }) > strong_binding_cutoff)
    draw_gradcams(dat_sequences, dat_gradcams, rbp_binding_scores, labels=dat_eclip_overlaps$RBP[idx], eclip_peaks=lapply(idx, function(i) dat_eclip_overlaps[dat_eclip_overlaps$query == dat_eclip_overlaps$query[i] & dat_eclip_overlaps$RBP == dat_eclip_overlaps$RBP[i],4:5]), dat_name=name, sequence_indices_to_display=dat_eclip_overlaps$query[idx])
}

##################################################################################################
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
setup_supermodel_data("asd1", num_tasks_to_do=1)
setup_supermodel_data("asd2", num_tasks_to_do=1)
setup_supermodel_data("asd3", num_tasks_to_do=1)

all_asd_dat <- readRDS(output_path("asd_dat_all.rds"))
    
setup_supermodel_data("chd", num_tasks_to_do=1)

# Interpet s prediction results for WGS datasets
analyze_s_predictions <- function(dat_names, case_pat=NULL, case_indices=NULL, cutoff=0.002, transcribed_only=FALSE, regional=FALSE, rewrite=FALSE) {
    to_keep_all <- c()
    s_vals_all <- abind(lapply(dat_names, function(dat_name) {
        print(paste0("Getting s predictions for ",dat_name))
        dat_s_preds_filename = output_path(paste0(dat_name,"_s_preds.rds"))
        if(!rewrite && file.exists(dat_s_preds_filename) && FALSE) {
            s_vals_all <- readRDS(dat_s_preds_filename)
        } else {
            to_keep <- load_supermodel_training_data(dat_name); dat_obs_exp[,1] <- rep(1,nrow(dat_obs_exp))
            #dat_gradcams <- readRDS(output_path(paste0(dat_name,"_gradcams.rds")))
            #dat_expected <- readRDS(output_path(paste0(dat_name,"_expected.rds")))
            ##dat_pLIs <- readRDS(output_path(paste0(dat_name,"_pLIs.rds")))
            #dat_obs_exp <- readRDS(output_path(paste0(dat_name,"_obs_exp.rds")))
            #dat_cpg <- readRDS(output_path(paste0(dat_name,"_cpg.rds")))
            activation_model <- get_activation_model(model, "s")
            arm_width = floor(activation_model$input[[1]]$shape[[2]]/2)
            midpoint = ceiling(ncol(dat_gradcams)/2)
            activations <- activation_model %>% predict(list(altref_gradcams[,(midpoint-arm_width):(midpoint+arm_width),rbp_order,,drop=FALSE], dat_gradcams[,(midpoint-arm_width):(midpoint+arm_width),rbp_order,,drop=FALSE], dat_obs_exp[,,drop=FALSE], expected[,,drop=FALSE], expected[,,drop=FALSE]))
            if(class(activations) == "list") {
                s_vals_all <- activations[[which(layers_to_investigate == "s")]]
            } else { # activations is a matrix instead of a list, caused by having just a single entity in layers_to_investigate
                s_vals_all <- activations
            }
            s_vals_all[s_vals_all == -Inf | s_vals_all == Inf] <- 0
            bases <- 1:4; names(bases) <- c("A","C","G","T")
            s_vals_all <- cbind(sapply(1:nrow(dat), function(i) s_vals_all[i,bases[paste0(dat$Alt[i])]]))
            saveRDS(s_vals_all, dat_s_preds_filename)
            to_keep_all <<- c(to_keep_all, to_keep)
        }
        return(s_vals_all)
    }), along=1)
    dat_granges_hg19 <- suppressWarnings(Reduce(c, lapply(dat_names, function(dat_name) {
        print(paste0("Getting GRanges (hg19) for ",dat_name))
        return(readRDS(output_path(paste0(dat_name,"_granges_hg19.rds"))))
    })))[to_keep_all]
    case_indices <- Reduce(c, lapply(dat_names, function(dat_name) {
        print(paste0("Getting case indices for ",dat_name))
        dat <- readRDS(output_path(paste0(dat_name,"_dat.rds")))
        if(is.null(case_indices)) { case_indices <- grepl(case_pat, dat$sample) }
        return(case_indices)
    }))[to_keep_all]
    if(transcribed_only) {
        if(!exists("fetal_brain_transcribed") || is.null(fetal_brain_transcribed)) { 
            #store_roadmap_eid_names()
            #get_relevant_roadmap_eids("ASD")
            #get_roadmap_epigenome_names(get_relevant_roadmap_eids("ASD"))
            fetal_brain_transcribed <- c(load_annotation("E081.H3K36me3.broadPeak"), load_annotation("E082.H3K36me3.broadPeak"))
            fetal_brain_transcribed <- intersect(fetal_brain_transcribed, fetal_brain_transcribed)
        }
        transcribed_indices <- unique(queryHits(findOverlaps(dat_granges_hg19, fetal_brain_transcribed)))
        s_vals_all <- s_vals_all[transcribed_indices,,drop=FALSE]
        dat_granges_hg19 <- dat_granges_hg19[transcribed_indices]
        case_indices <- case_indices[transcribed_indices]
        transcribed_mtext = " (transcribed in fetal brain)"
        transcribed_filename_text = "_transcribed"
    } else { transcribed_mtext = ""; transcribed_filename_text = "" }
    if(regional) { reg_idx = 1:ncol(s_vals_all); regional_filename_text = "_regional" } else { reg_idx = round((ncol(s_vals_all)+1)/2); regional_filename_text = "" }
    control_indices <- !case_indices; num_case_variants = sum(case_indices); num_control_variants = sum(control_indices)
    
    dat_name = gsub("[0-9_]*$","",dat_names[1])
    sum(s_vals_all[case_indices,reg_idx] > cutoff)/(num_case_variants*length(reg_idx)); sum(s_vals_all[control_indices,reg_idx] > cutoff)/(sum(num_control_variants)*length(reg_idx))
    
    filename = output_path(paste0(dat_name,transcribed_filename_text,regional_filename_text,"_s_distributions.pdf"))
    pdf(filename)
    plot(density(log10(s_vals_all[case_indices,reg_idx])), col="red", xlab="Selection Coefficient log10(s)", main=paste0(toupper(dat_name)," Selection Coef. Distributions"), cex.lab=1.4, cex.axis=1.4, cex.main=1.3)
    lines(density(log10(s_vals_all[control_indices,reg_idx])), col="blue")
    mtext(paste0("All variants",transcribed_mtext), cex=1.2)
    variant_counts <- c(prod(dim(s_vals_all[case_indices,reg_idx,drop=FALSE])),prod(dim(s_vals_all[control_indices,reg_idx,drop=FALSE])))
    print(paste0("Ratio: ", variant_counts[1]/variant_counts[2]))
    total_variant_counts <- variant_counts
    legend("topright", legend=c(paste0(c("case (N = ","control (N = "),variant_counts,")")), col=c("red","blue"), pch=15, cex=1.2)
    dev.off()
    pdf_to_png(filename)
    
    cutoffs <- seq(sqrt(max(min(s_vals_all[control_indices,reg_idx]),min(s_vals_all[control_indices,reg_idx]))), sqrt(min(max(s_vals_all[control_indices,reg_idx]),max(s_vals_all[control_indices,reg_idx]))), by=0.0025)**2
    fet_results <- unfactorize(data.frame(t(data.frame(sapply(cutoffs, function(cutoff) {
        variant_counts <- c(sum(s_vals_all[case_indices,reg_idx]>=cutoff),sum(s_vals_all[control_indices,reg_idx]>=cutoff))
        fet_result <- fisher_exact_test(variant_counts[1], variant_counts[2], total_variant_counts[1], total_variant_counts[2], alternative=c("two.sided"))
        return(c(cutoff, fet_result[["estimate"]], fet_result[["ci"]], fet_result[["p.value"]]))
    }))))); colnames(fet_results) <- c("s_cutoff", "OR", "conf.int_lower", "conf.int_upper", "p.value")
    optimal_cutoff_index = which(fet_results$OR < Inf & fet_results$OR > 1)[which.min(fet_results$p.value[fet_results$OR < Inf & fet_results$OR > 1])]
    cutoff <- cutoffs[optimal_cutoff_index]
    fet_result <- fet_results[optimal_cutoff_index,,drop=FALSE]
    filename = output_path(paste0(dat_name,transcribed_filename_text,regional_filename_text,"_s_multiple_threshold_test.pdf"))
    pdf(filename)
    cols <- rep("black", length(cutoffs)); cols[fet_results$p.value < 0.05] <- "red"
    y_lim <- 1.05*c(min(log2(fet_results$conf.int_lower[fet_results$conf.int_lower!=-Inf & fet_results$conf.int_lower!=0])),max(log2(fet_results$conf.int_upper[fet_results$conf.int_upper!=Inf & fet_results$conf.int_upper!=0])))
    print(y_lim)
    plot(log10(cutoffs), log2(fet_results$OR), col=cols, cex.axis=1.4, cex.lab=1.4, cex.main=1.3, pch=16, main="Selection Coef. Multiple Threshold Test", xlab="log10(s) cutoff", ylab="log2(Odds Ratio)", ylim=y_lim)
    abline(h=0, lty=3, col="gray50")
    for(cutoff_index in 1:length(cutoffs)) {
        segments(x0=log10(cutoffs[cutoff_index]), y0=log2(fet_results$conf.int_lower[cutoff_index]), x1=log10(cutoffs[cutoff_index]), y1=log2(fet_results$conf.int_upper[cutoff_index]), col=cols[cutoff_index], lwd=1)
    }
    #abline(v=cutoff, lty=3, col="gray50")
    mtext(paste0("Optimal threshold: s >= ",round(cutoff,4)," (OR = ",round(fet_result$OR,2),", p = ",formatC(fet_result$p.value, format="e", digits=2),", two-sided FET)"), cex=1.2)
    #legend("topleft", legend=c(paste0("labels < ",split_point), paste0("labels >= ",split_point)), col=c("red", "blue"), pch=15, cex=1.3)
    dev.off()
    pdf_to_png(filename)
    
    cutoff = 0.0038
    filename = output_path(paste0(dat_name,transcribed_filename_text,regional_filename_text,"_s_distributions_righttail.pdf"))
    pdf(filename)
    plot(density(log10(s_vals_all[case_indices,reg_idx])), xlim=log10(c(cutoff,max(s_vals_all[,reg_idx])*1.05)), xaxs="i", ylim=c(0,20), col="red", xlab="Selection Coefficient (s)", main=paste0(toupper(dat_name)," Selection Coef. Distributions"), cex.lab=1.4, cex.axis=1.4, cex.main=1.3)
    lines(density(log10(s_vals_all[control_indices,reg_idx])), col="blue")
    mtext(paste0("Variants",transcribed_mtext," with s >= ",round(cutoff,4)," only"), cex=1.2)
    abline(v=cutoff, lty=3, col="gray50")
    variant_counts <- c(sum(s_vals_all[case_indices,reg_idx]>=cutoff),sum(s_vals_all[control_indices,reg_idx]>=cutoff))
    print(paste0("Ratio: ", variant_counts[1]/variant_counts[2]))
    #fet_result <- fisher_exact_test(variant_counts[1], variant_counts[2], total_variant_counts[1], total_variant_counts[2], alternative=c("two.sided"))
    legend("topright", legend=c(paste0(c("case (N = ","control (N = "),variant_counts,")"), paste0("Odds Ratio = ",round(fet_result$OR,2)), paste0("(p = ",formatC(fet_result$p.value, format="e", digits=2),", two-sided FET)")), col=c("red","blue","white","white"), pch=15, cex=1.2)
    dev.off()
    pdf_to_png(filename)
    
    if(FALSE) {
        dat_eclip_overlaps <- find_eclip_overlaps(dat_granges_hg19)
        case_variants_eclip <- intersect(unique(dat_eclip_overlaps$query), which(case_indices))
        control_variants_eclip <- intersect(unique(dat_eclip_overlaps$query), which(control_indices))
        print(paste0("Freq. regions with nearby eCLIP in cases: ",length(case_variants_eclip) / num_case_variants))
        print(paste0("Freq. regions with nearby eCLIP in controls: ",length(control_variants_eclip) / num_control_variants))
        
        filename = output_path(paste0(dat_name,transcribed_filename_text,regional_filename_text,"_s_distributions_eclip.pdf"))
        pdf(filename)
        plot(density(s_vals_all[case_variants_eclip,reg_idx]), col="red", xlab="Selection Coefficient (s)", main=paste0(toupper(dat_name)," Selection Coef. Distributions"), cex.lab=1.4, cex.axis=1.4, cex.main=1.3)
        lines(density(s_vals_all[control_variants_eclip,reg_idx]), col="blue")
        mtext(paste0("Variants",transcribed_mtext," with nearby eCLIP only"), cex=1.2)
        variant_counts <- c(length(case_variants_eclip),length(control_variants_eclip))
        print(paste0("Ratio: ", variant_counts[1]/variant_counts[2]))
        fet_result <- fisher_exact_test(variant_counts[1], variant_counts[2], total_variant_counts[1], total_variant_counts[2], alternative=c("two.sided"))
        legend("topright", legend=c(paste0(c("case (N = ","control (N = "),variant_counts,")"), paste0("Odds Ratio = ",round(fet_result[["estimate"]],2)), paste0("(p = ",formatC(fet_result[["p.value"]], format="e", digits=2),", two-sided FET)")), col=c("red","blue","white","white"), pch=15, cex=1.2)
        dev.off()
        pdf_to_png(filename)
        
        filename = output_path(paste0(dat_name,transcribed_filename_text,regional_filename_text,"_s_distributions_eclip_righttail.pdf"))
        pdf(filename)
        plot(density(s_vals_all[case_variants_eclip,reg_idx]), xlim=c(0.02,0.11), xaxs="i", ylim=c(0,2), col="red", xlab="Selection Coefficient (s)", main=paste0(toupper(dat_name)," Selection Coef. Distributions"), cex.lab=1.4, cex.axis=1.4, cex.main=1.3)
        lines(density(s_vals_all[control_variants_eclip,reg_idx]), col="blue")
        mtext(paste0("Variants",transcribed_mtext," with nearby eCLIP and s >= ",cutoff," only"), cex=1.2)
        abline(v=cutoff, lty=3, col="gray50")
        variant_counts <- c(sum(s_vals_all[case_variants_eclip,reg_idx]>=cutoff),sum(s_vals_all[control_variants_eclip,reg_idx]>=cutoff))
        print(paste0("Ratio: ", variant_counts[1]/variant_counts[2]))
        fet_result <- fisher_exact_test(variant_counts[1], variant_counts[2], total_variant_counts[1], total_variant_counts[2], alternative=c("two.sided"))
        legend("topright", legend=c(paste0(c("case (N = ","control (N = "),variant_counts,")"), paste0("Odds Ratio = ",round(fet_result[["estimate"]],2)), paste0("(p = ",formatC(fet_result[["p.value"]], format="e", digits=2),", two-sided FET)")), col=c("red","blue","white","white"), pch=15, cex=1.2)
        dev.off()
        pdf_to_png(filename)
    }
    
    return(fet_results)
}
#analyze_s_predictions("chd", case_pat="^case_")
#model <- load_supermodel("suprnova")
a <- analyze_s_predictions(c("asd1","asd2"), case_pat="-p", cutoff=0.002, transcribed_only=TRUE, regional=FALSE, rewrite=TRUE)
a <- analyze_s_predictions(c("asd1","asd2"), case_pat="-p", cutoff=0.002, transcribed_only=FALSE, regional=FALSE, rewrite=TRUE)
#a <- analyze_s_predictions("chd", case_pat="^case_", cutoff=0.002, transcribed_only=TRUE, regional=FALSE, rewrite=TRUE)
#a <- analyze_s_predictions("chd", case_pat="^case_", cutoff=0.002, transcribed_only=FALSE, regional=FALSE, rewrite=TRUE)
a <- analyze_s_predictions("hgmd", case_pat="_Regulatory", cutoff=0.002, transcribed_only=TRUE, regional=FALSE, rewrite=TRUE)
a <- analyze_s_predictions("hgmd", case_pat="_Regulatory", cutoff=0.002, transcribed_only=FALSE, regional=FALSE, rewrite=TRUE)
#a <- analyze_s_predictions(c("asd","asd1","asd2","asd3"), case_pat="-p", cutoff=0.04, transcribed_only=TRUE, regional=TRUE)

# STOPPED HERE!!! TRY HGMD vs. RANDOM GNOMAD DATASET AND SEE IF THIS IS AT LEAST AS GOOD AS OUR GRADIENT BOOSTING METHOD!



#################################################################################################################
# REGION TYPE ANALYSIS
#################################################################################################################
# Randomly sample variants across different region types for comparison.
# sampled_region_variants_filename = output_path("sampled_region_type_variants_dat.rds")
# saveRDS(sample_region_type_variants(region_types, regions_to_sample_per_type=20000, width=151), sampled_region_variants_filename)

sampled_gnomad_variants_filename = output_path("random_transcribed_gnomAD_variants_100k.rds")
saveRDS(sample_gnomad_variants(100000, "transcribed"), sampled_gnomad_variants_filename)
dat_100k <- readRDS(sampled_gnomad_variants_filename)
colnames(dat_100k)[1:5] <- c("Chrom", "Position", "Ref", "Alt", "sample")
dat_100k <- liftover(dat_100k, from="hg38", to="hg19", chr_colname="Chrom", start_colname="Position", ref_colname="Ref", alt_colname="Alt", confirm_refseq=FALSE, mismatches_pause=FALSE)
saveRDS(dat_100k, output_path("random_transcribed_gnomAD_variants_100k_hg19.rds"))
dat_100k <- readRDS(output_path("random_transcribed_gnomAD_variants_100k_hg19.rds"))
setup_supermodel_data("transcribed100k_1", dat=dat_100k[1:25000,c("Chrom", "Position", "Ref", "Alt", "sample")], version="hg19", num_tasks_to_do=1)
setup_supermodel_data("transcribed100k_1", num_tasks_to_do=20)
setup_supermodel_data("transcribed100k_2", dat=dat_100k[25001:50000,c("Chrom", "Position", "Ref", "Alt", "sample")], version="hg19", num_tasks_to_do=1)
setup_supermodel_data("transcribed100k_2", num_tasks_to_do=20)
setup_supermodel_data("transcribed100k_3", dat=dat_100k[50001:75000,c("Chrom", "Position", "Ref", "Alt", "sample")], version="hg19", num_tasks_to_do=1)
setup_supermodel_data("transcribed100k_3", num_tasks_to_do=20)
setup_supermodel_data("transcribed100k_4", dat=dat_100k[75001:100000,c("Chrom", "Position", "Ref", "Alt", "sample")], version="hg19", num_tasks_to_do=1)
setup_supermodel_data("transcribed100k_4", num_tasks_to_do=20)
# Annotate AFs and region types for the random variants
for(transcribed100k_i in c(1:4)) {
    transcribed100k_name = paste0("transcribed100k_",transcribed100k_i)
    #transcribed100k_name = "asd2"
    print(transcribed100k_name)
    dat <- readRDS(output_path(paste0(transcribed100k_name,"_dat.rds")))
    region_types <- annotate_region_types(dat, version="hg38", region_types=c("CDS","five_prime_UTR","three_prime_UTR","five_prime_ss","three_prime_ss","intron","intergenic"), breakdown_CDS=FALSE, num_variants_to_sample=NULL)
    saveRDS(region_types, output_path(paste0(transcribed100k_name,"_region_types.rds")))
    dat_AFs <- matrix(0, nrow=nrow(dat), ncol=151)
    dat_AFs[,76] <- dat$sample #dat_100k$AF[(25000*(transcribed100k_i-1)+1):(25000*transcribed100k_i)]
    saveRDS(dat_AFs, output_path(paste0(transcribed100k_name,"_AFs.rds")))
}


# Do an unbiased MAPS
sample_size = 71702
dat <- readRDS(output_path("random_transcribed_gnomAD_variants_100k.rds"))

trimers_expected_ps <- readRDS(output_path("trimers_expected_ps.rds"))
dat_trimers <- readRDS("/data/hg38_gnomad3.0_genome_snvs_trimers.rds")
dat_trimers <- dat_trimers[dat_trimers$AF < 1e-3,]

a_AC <- round(a$AF * sample_size * 2 + 0.01)
a_trimers <- a$ref_sequence
sort(table(a_AC), decreasing=TRUE)

            
b <- rbind(cbind(sum(a_AC[is_cds] == 1) / sum(a_AC[is_cds] > 0), sum(grepl("CG",a_trimers[is_cds]))/sum(is_cds), mean(na.omit(trimers_expected_ps[a_trimers[is_cds]]))),
        cbind(sum(a_AC[is_5utr] == 1) / sum(a_AC[is_5utr] > 0), sum(grepl("CG",a_trimers[is_5utr]))/sum(is_5utr), mean(na.omit(trimers_expected_ps[a_trimers[is_5utr]]))),
        cbind(sum(a_AC[is_3utr] == 1) / sum(a_AC[is_3utr] > 0), sum(grepl("CG",a_trimers[is_3utr]))/sum(is_3utr), mean(na.omit(trimers_expected_ps[a_trimers[is_3utr]]))),
        cbind(sum(a_AC[is_intron] == 1) / sum(a_AC[is_intron] > 0), sum(grepl("CG",a_trimers[is_intron]))/sum(is_intron), mean(na.omit(trimers_expected_ps[a_trimers[is_intron]]))),
        cbind(sum(a_AC[is_5ss] == 1) / sum(a_AC[is_5ss] > 0), sum(grepl("CG",a_trimers[is_5ss]))/sum(is_5ss), mean(na.omit(trimers_expected_ps[a_trimers[is_5ss]]))),
        cbind(sum(a_AC[is_3ss] == 1) / sum(a_AC[is_3ss] > 0), sum(grepl("CG",a_trimers[is_3ss]))/sum(is_3ss), mean(na.omit(trimers_expected_ps[a_trimers[is_3ss]]))),
        cbind(sum(a_AC[is_intergenic] == 1) / sum(a_AC[is_intergenic] > 0), sum(grepl("CG",a_trimers[is_intergenic]))/sum(is_intergenic), mean(na.omit(trimers_expected_ps[a_trimers[is_intergenic]]))))

region_classes <- list(is_cds,is_cds_syn,is_cds_startloss,is_cds_missense,is_cds_stoploss,is_cds_nonsense,is_5utr,is_3utr,is_intron,is_5ss,is_3ss,is_intergenic)
region_class_names <- c("All CDS","Synonymous","Start lost","Missense","Stop lost","Nonsense","5'UTR","3'UTR","intron","5'ss","3'ss","intergenic")
num_bootstrap_samples = 10
b <- unfactorize(data.frame(rbindlist(lapply(1:length(region_classes), function(class_indices_i) {
        print(paste0(class_indices_i,". ",region_class_names[class_indices_i]))
        class_indices <- region_classes[[class_indices_i]]
        AC <- a_AC[class_indices]; trimers <- a_trimers[class_indices]; N = sum(class_indices)
        x_sampled <- t(sapply(1:num_bootstrap_samples, function(sample_i) {
            if(N > 100000000) { gc() }
            print(sample_i)
            sample_indices <- sample(1:N, replace=TRUE)
            ps_obs <- sum(AC[sample_indices] == 1) / sum(AC[sample_indices] > 0)
            ps_exp <- mean(na.omit(trimers_expected_ps[trimers[sample_indices]]))
            rm(sample_indices)
            return(c(ps_obs, ps_exp))
        })); x_sampled <- cbind(x_sampled, x_sampled[,1] - x_sampled[,2]); colnames(x_sampled) <- c("PS_obs", "PS_exp", "MAPS")
        return(data.frame(t(data.frame(c(region_class_names[class_indices_i], N, sum(grepl("CG",trimers))/N, colMeans(x_sampled), x_sampled[max(c(1,floor(0.025*num_bootstrap_samples))),], x_sampled[min(c(ceiling(0.975*num_bootstrap_samples),num_bootstrap_samples)),])))))
}))))
colnames(b) <- c("region_type","N","CpG_freq","PS_obs","PS_exp","MAPS","PS_obs_lower","PS_exp_lower","MAPS_lower","PS_obs_upper","PS_exp_upper","MAPS_upper")
b
saveRDS(b, output_path("b.rds"))

MAPS_offset <- sum(b$MAPS*b$N)/sum(b$N)

filename = output_path(paste0("b.pdf"))
pdf(filename)
cols = c("black","black","white",rainbow(nrow(b)-1),"gray40")
plot(b$CpG_freq, b$PS_obs, pch=16, xlab="CpG trimer frequency", ylim=range(b$PS_exp+MAPS_offset)*c(0.96,1), ylab="Proportion Singletons", main="Regional PS Breakdown", col=cols[-c(1:3)], cex.lab=1.4, cex.axis=1.4, cex.main=1.3)
points(b$CpG_freq, b$PS_exp+MAPS_offset, pch=8, col=cols[-c(1:3)])
legend("bottomleft", legend=c("Observed PS in gnomAD","Expected PS from trimers","",paste0(b$region_type," (MAPS = ",round(b$MAPS-MAPS_offset,4),")")), col=cols, pch=c(16,8,8,rep(15,nrow(b))), cex=1.1)
if(num_gnomAD_variants_to_sample < nrow(dat_trimers)) { mtext(paste0(num_gnomAD_variants_to_sample," random genomic variants"), cex=1.2)
} else { mtext(paste0(num_gnomAD_variants_to_sample," gnomAD 3.0 variants"), cex=1.2) }
dev.off()
pdf_to_png(filename)
    
filename = output_path(paste0("regional_MAPS_comparison_gnomAD.pdf"))
pdf(filename)
cols = c("black","black","white",rainbow(nrow(b)-1),"gray40")
plot(1:length(b$region_type), rev(b$MAPS-MAPS_offset), pch=16, xlab="", ylab="MAPS", main="Regional MAPS Distributions", ylim=range(c(b[,c("MAPS_lower","MAPS_upper")]))-MAPS_offset, col=rev(cols[-c(1:3)]), xaxt="n", cex.axis=1.3, cex.lab=1.3, cex.main=1.3)
segments(x0=1:length(b$region_type), y0=rev(b$MAPS_lower)-MAPS_offset, x1=1:length(b$region_type), y1=rev(b$MAPS_upper)-MAPS_offset, col=rev(cols[-c(1:3)]))
text(1:length(b$region_type)+0.2, par("usr")[3]-0.001, srt=45, adj=1, xpd=TRUE, labels=rev(b$region_type), cex=0.8)
if(num_gnomAD_variants_to_sample < nrow(dat_trimers)) { mtext(paste0(num_gnomAD_variants_to_sample," random genomic variants"), cex=1.2)
} else { mtext(paste0(num_gnomAD_variants_to_sample," gnomAD 3.0 variants"), cex=1.2) }
dev.off()
pdf_to_png(filename)

#################################################################################################################
# Read WGS data frames
#################################################################################################################

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

#################################################################################################################
# Find case-control enrichment boost when predicted s or d are included.
#################################################################################################################
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

#################################################################################################################
# Analyze s distributions in different regional bins (3'UTR,  CDS, splice sites, etc.)
#################################################################################################################
dat_name = "regional2"
dat_s_preds_filename = output_path(paste0(dat_name,"_s_preds.rds"))
if(file.exists(dat_s_preds_filename)) {
    s_vals_all <- readRDS(dat_s_preds_filename)
} else {
    dat_gradcams <- readRDS(output_path(paste0(dat_name,"_gradcams.rds")))
    dat_expected <- readRDS(output_path(paste0(dat_name,"_expected.rds")))
    #dat_pLIs <- readRDS(output_path(paste0(dat_name,"_pLIs.rds")))
    dat_pLIs <- readRDS(output_path(paste0(dat_name,"_obs_exp.rds")))
    dat_cpg <- readRDS(output_path(paste0(dat_name,"_cpg.rds")))
    activation_model <- get_activation_model(model, "s")
    activations <- activation_model %>% predict(list(dat_gradcams[,,,,drop=FALSE], dat_pLIs[,,drop=FALSE], dat_expected[,,drop=FALSE]*mu_scaling_factor))
    if(class(activations) == "list") {
        s_vals_all <- 1/(mu_scaling_factor*activations[[which(layers_to_investigate == "s")]])
    } else { # activations is a matrix instead of a list, caused by having just a single entity in layers_to_investigate
        s_vals_all <- 1/(mu_scaling_factor*activations)
    }
    s_vals_all[s_vals_all == -Inf | s_vals_all == Inf] <- 0
    saveRDS(s_vals_all, dat_s_preds_filename)
}


regional_dat <- readRDS(output_path("regional2_dat.rds"))
regional_dat_labels <- data.frame(t(data.frame(strsplit(regional_dat$sample, "\\.")))); colnames(regional_dat_labels) <- c("region_type", "index"); rownames(regional_dat_labels) <- NULL
sampled_region_type_variants <- data.frame(readRDS(output_path("sampled_region_type_variants_dat.rds")))
sampled_region_type_variants <- sampled_region_type_variants[paste0(sampled_region_type_variants$region_type,".",sampled_region_type_variants$index) %in% regional_dat$sample,]
region_offsets_filename = output_path("regional2_offsets.rds")
if(file.exists(region_offsets_filename)) {
    region_offsets <- readRDS(region_offsets_filename)
} else {
    region_offsets <- t(sapply(1:nrow(regional_dat), function(i) { 
        if(i %% 1000 == 1) { print(paste0(i," / ",nrow(regional_dat))) }
        label = unlist(regional_dat_labels[i,])
        region_entry <- sampled_region_type_variants[which(sampled_region_type_variants$region_type == paste0(label[1]) & sampled_region_type_variants$index == paste0(label[2])),] 
        start_offset = region_entry$region_start - region_entry$start; end_offset = region_entry$end - region_entry$region_end
        return(c(start_offset, end_offset))
    })); colnames(region_offsets) <- c("start_offset", "end_offset")
    saveRDS(region_offsets, region_offsets_filename)
}

library("GenomicAlignments")
genebody_hg19 <- read.csv(data_path("refGene\\refGene_hg19_genebody_fixed.bed"), sep="\t", stringsAsFactors=FALSE)
genebody_granges_hg19 <- to_genomic_regions(genebody_hg19, chr_colname="chromosome", start_colname="tss", end_colname="tes", strand_colname="strand", label_colname="gene", order_coordinates=TRUE, remove_duplicate_labels=TRUE)
genes <- readRDS(output_path("genes.rds"))
#nearest_genes_dat <- names(genebody_granges_hg19)[nearest(dat_variant_granges_hg19, genebody_granges_hg19)]
#pLIs_dat <- genes[unlist(sapply(nearest_genes_dat,function(x) { nearest_gene <- which(paste0(exac_dat$gene) == x); if(length(nearest_gene) < 1) { return(NA) } else { return(nearest_gene[1]) }  })),"o_e_Z"]
regional_granges <- readRDS(output_path("regional2_granges_hg19.rds"))
regional_granges[grepl("five_prime_UTR",regional_dat$sample)]

a <- findOverlaps(regional_granges[grepl("five_prime_UTR",regional_dat$sample)], regional_granges[grepl("three_prime_UTR",regional_dat$sample)])
a
a <- findOverlaps(regional_granges[grepl("five_prime_UTR",regional_dat$sample)], regional_granges[grepl("intron",regional_dat$sample)])
a
a <- findOverlaps(regional_granges[grepl("five_prime_UTR",regional_dat$sample)], regional_granges[grepl("CDS",regional_dat$sample)])
a
a <- findOverlaps(regional_granges[grepl("CDS",regional_dat$sample)], regional_granges[grepl("intron",regional_dat$sample)])
a
a <- findOverlaps(regional_granges[grepl("CDS",regional_dat$sample)], genebody_granges_hg19)
a

gencode <- readRDS(output_path(paste0("gencode_full_dat.rds")))


region_types <- paste0(unlist(unique(regional_dat_labels$region_type)))
cols <- rainbow(length(region_types))
regional_AFs <- readRDS(output_path("regional2_AFs.rds"))
regional_expected <- readRDS(output_path("regional2_expected.rds"))
regional_trimers <- readRDS(output_path("regional2_trimers.rds"))
regional_cpg <- readRDS(output_path("regional2_cpg.rds"))
filename = output_path(paste0("regional2_s_distributions.pdf"))
sample_size = 71702
pdf(file=filename)
maps_result <- unfactorize(data.frame(rbindlist(lapply(1:length(region_types), function(region_i) { #data.frame(t(sapply
    region_type = region_types[region_i]
    print(paste0(region_type))
    regional_indices <- which(regional_dat_labels[,1] == region_type)
    s_vals_region <- unlist(sapply(regional_indices, function(i) { s_vals_all[i,region_offsets[i,1]:(151-region_offsets[i,2])] }))
    AC_region <- unlist(sapply(regional_indices, function(i) { regional_AFs[i,region_offsets[i,1]:(151-region_offsets[i,2])] }))
    background_region <- unlist(sapply(regional_indices, function(i) { regional_expected[i,region_offsets[i,1]:(151-region_offsets[i,2])] }))
    trimers_region <- unlist(sapply(regional_indices, function(i) { paste0(unlist(regional_trimers[i,region_offsets[i,1]:(151-region_offsets[i,2])])) }))
    cpg_region <- unlist(sapply(regional_indices, function(i) { regional_cpg[i,region_offsets[i,1]:(151-region_offsets[i,2])] }))
    
    if(region_type == region_types[1]) {
        plot(density(s_vals_region), xlim=c(0.01,0.1), ylim=c(0,4), col=cols[region_i], xlab="Selection Coef. s", main="Regional Selection. Coef Distributions")
    } else {
        lines(density(s_vals_region), col=cols[region_i])
    }
    #s_vals_region_bootstrap <- sort(sample(s_vals_region,10000,replace=TRUE))
    
    AC_region <- round(AC_region * sample_size * 2 + 0.01)
    maps_res <- data.frame(t(rbind(AC_region, background_region, s_vals_region, trimers_region, cpg_region, region_type)))
    colnames(maps_res) <- c("AC", "mutability", "s", "trimer", "cpg", "region_type")
    return(maps_res)
    #return(c(mean(s_vals_region_bootstrap), s_vals_region_bootstrap[250], s_vals_region_bootstrap[9750]))
}))));# colnames(s_vals_all_means) <- c("estimate", "conf.int_lower", "conf.int_higher")
legend("topright", legend=c(paste0(region_types)), col=cols, pch=15)
dev.off()
pdf_to_png(filename)

saveRDS(maps_result, output_path("regional2_maps_result.rds"))

#maps_model <- readRDS(output_path("maps_model.rds"))
trimers_expected_ps <- readRDS(output_path("trimers_expected_ps.rds"))
s_by = 0.01
#s_bucket_cuts <- c(0,0.01,0.04,1)
s_bucket_cuts <- rbind(s_bucket_cuts[-length(s_bucket_cuts)], s_bucket_cuts[-1])
s_bucket_cuts <- rbind(c(0,0,0,rep(s_by*3,1/s_by-3)), c(s_by,seq(s_by*2,1,by=s_by)))
s_bucket_cuts <- rbind(c(0,0,0,0,0,0.05),c(0.01,0.02,0.03,0.4,0.5,1))
s_bucket_names <- c(paste0(s_bucket_cuts[1,]," < s <= ",s_bucket_cuts[2,]), "s = any") #"s = 0", 
region_types <- region_types[region_types != "noncoding_exon"]
for(sites in c("non-CpG", "all")) { #CpG 
    #sites = "CpG" # "all", "CpG", or "non-CpG"
    print(sites)
    ##region_types <- list("noncoding_exon", c("five_prime_ss", "three_prime_ss"), "intron", "CDS")
    ##region_type_names <- c("UTRs", "splice sites", "introns", "CDS")
    maps_sfs <- unfactorize(data.frame(rbindlist(lapply(1:(ncol(s_bucket_cuts)+1), function(s_i) {
        print(s_i)
        if(s_i == 0) { curr_indices_backup <- maps_result$s == 0 } else if(s_i > ncol(s_bucket_cuts)) { curr_indices_backup <- rep(TRUE, nrow(maps_result)) } else { curr_indices_backup <- maps_result$s > s_bucket_cuts[1,s_i] & maps_result$s <= s_bucket_cuts[2,s_i]  }
        return(unfactorize(data.frame(rbindlist(lapply(0:length(region_types), function(region_i) {
            curr_indices <- curr_indices_backup
            if(region_i == 0) { region_type = "all"; regional_indices <- rep(TRUE, length(curr_indices)) 
            }else { region_type = paste0(region_types[region_i]); regional_indices <- maps_result$region_type == region_type }
            print(region_type)
            curr_indices <- curr_indices & regional_indices
            print(sum(curr_indices))
            maps_cpg_vector <- maps_result$cpg == "TRUE"
            if(sites == "CpG") { curr_indices <- curr_indices & maps_cpg_vector 
            } else if (sites == "non-CpG") { curr_indices <- curr_indices & !maps_cpg_vector
            } else { sites = "all" }
            x <- as.numeric(paste0(maps_result$AC))
            ps <- sum(x[curr_indices] == 1)/sum(x[curr_indices] >= 1)
            s <- as.numeric(paste0(maps_result$s))
            #expected_ps <- predict(maps_model, new=data.frame(x=log10(maps_result$mutability[curr_indices])))
            reg_trimers <- paste0(maps_result$trimer); #reg_trimers <- reg_trimers[nchar(reg_trimers)==3]
            expected_ps <- trimers_expected_ps[reg_trimers]; expected_ps[is.na(expected_ps)] <- 0
            curr_indices_which <- which(c(curr_indices)); curr_indices_sum <- sum(curr_indices)
            x_flat <- c(x); expected_ps_flat <- c(expected_ps); s_flat <- c(s)
            num_bootstrap_samples = 100
            maps_sampled <- sort(sapply(1:num_bootstrap_samples, function(sample_i) { 
                sample_indices <- sample(curr_indices_which, curr_indices_sum, replace=TRUE)
                return(mean((sum(x_flat[sample_indices] == 1)/sum(x_flat[sample_indices] >= 1))) - mean(expected_ps_flat[sample_indices]))
            }), decreasing=FALSE)
            maps <- ps - mean(expected_ps[curr_indices])
            x <- x[curr_indices]; expected_ps <- expected_ps[curr_indices]
            s_sampled <- sort(sapply(1:num_bootstrap_samples, function(sample_i) { 
                sample_indices <- sample(curr_indices_which, curr_indices_sum, replace=TRUE)
                return(mean(s_flat[sample_indices]))
            }), decreasing=FALSE)
            s <- mean(s[curr_indices])
            return(data.frame(t(data.frame(c(sum(x==0),sum(x==1),sum(x==2),sum(x==3),sum(x==4),sum(x==5),sum(x>5),ps,mean(expected_ps),maps,maps_sampled[round(num_bootstrap_samples*0.025)],maps_sampled[round(num_bootstrap_samples*0.975)],s,s_sampled[round(num_bootstrap_samples*0.025)],s_sampled[round(num_bootstrap_samples*0.975)],region_type,curr_indices_sum)[-1]))))
        })))))
    }))))
    #maps_sfs <- t(apply(maps_sfs, 1, function(x) x/sum(x)))
    colnames(maps_sfs) <- c("Singletons", "Doubletons", "Tripletons", "AC = 4", "AC = 5", "AC > 5", "PS_observed", "PS_expected", "MAPS", "MAPS_lower", "MAPS_upper", "s", "s_lower", "s_upper", "region_type","N_sites")
    rownames(maps_sfs) <- c(sapply(s_bucket_names, function(x) paste0(x,", ",unique(maps_sfs$region_type))))
    maps_sfs$MAPS_upper[is.na(maps_sfs$MAPS_upper)] <- 0; maps_sfs$MAPS_upper <- as.numeric(maps_sfs$MAPS_upper)
    maps_sfs$MAPS_lower[is.na(maps_sfs$MAPS_lower)] <- 0; maps_sfs$MAPS_lower <- as.numeric(maps_sfs$MAPS_lower)
    s_last_filled = min(c(max(which(apply(maps_sfs[-nrow(maps_sfs),1:6], 1, function(x) sum(!is.nan(x) & x > 0)>0))), max(which(diff(maps_sfs$MAPS)!=0))+1))
    maps_sfs_rownames <- c(rownames(maps_sfs)[1:(s_last_filled-1)], gsub("<=.*,", "<= 1,", rownames(maps_sfs)[s_last_filled]), rownames(maps_sfs)[nrow(maps_sfs)])
    maps_sfs <- maps_sfs[c(1:s_last_filled,nrow(maps_sfs)),]; rownames(maps_sfs) <- maps_sfs_rownames
    maps_sfs
    maps_sfs_to_write <- cbind(rownames(maps_sfs),maps_sfs); colnames(maps_sfs_to_write)[1] <- "AC"
    write.table(maps_sfs_to_write, sep=",", output_path(paste0("model_sfs_",sites,"_table.csv")), row.names=FALSE)
    
    filename = output_path(paste0("regional_MAPS_comparison_",sites,".pdf"))
    pdf(file=filename)
    #maps_sfs <- maps_sfs[maps_sfs$region_type != "noncoding_exon",]
    maps_sfs_names <- unfactorize(data.frame(strsplit(rownames(maps_sfs),", ")))
    region_buckets <- unique(paste0(maps_sfs_names[2,]))
    cols <- rainbow(length(region_buckets))
    any_s <- grepl(" = any",rownames(maps_sfs))
    plot(1:length(region_buckets), maps_sfs$MAPS[any_s], col=cols, main="Regional MAPS Distributions", xlab="", ylab="MAPS", ylim=range(c(maps_sfs[any_s,c("MAPS_lower","MAPS_upper")])), xaxt="n", cex.axis=1.3, cex.lab=1.3, cex.main=1.3)
    #abline(h=0, col="gray50", lty=3)
    segments(x0=1:length(region_buckets), y0=maps_sfs$MAPS_lower[any_s], x1=1:length(region_buckets), y1=maps_sfs$MAPS_upper[any_s], col=cols)
    text(1:length(region_buckets)+0.2, par("usr")[3]-0.001, srt=45, adj=1, xpd=TRUE, labels=region_buckets, cex=0.8)
    dev.off()
    pdf_to_png(filename)
    
    filename = output_path(paste0("regional_s_comparison_",sites,".pdf"))
    pdf(file=filename)
    plot(1:length(region_buckets), maps_sfs$s[any_s], col=cols, main="Regional Selection Coef. Distributions", xlab="", ylab="mean Selection Coef. s", ylim=range(c(maps_sfs[any_s,c("s_lower","s_upper")])), xaxt="n", cex.axis=1.3, cex.lab=1.3, cex.main=1.3)
    segments(x0=1:length(region_buckets), y0=maps_sfs$s_lower[any_s], x1=1:length(region_buckets), y1=maps_sfs$s_upper[any_s], col=cols)
    text(1:length(region_buckets)+0.2, par("usr")[3]-0.00001, srt=45, adj=1, xpd=TRUE, labels=region_buckets, cex=0.8)
    dev.off()
    pdf_to_png(filename)
    
    filename = output_path(paste0("model_MAPS_regional_",sites,".pdf"))
    pdf(file=filename)
    maps_sfs <- maps_sfs[maps_sfs$region_type != "noncoding_exon" & !grepl(" = ",rownames(maps_sfs)),]
    maps_sfs_names <- unfactorize(data.frame(strsplit(rownames(maps_sfs),", ")))
    s_buckets <- unique(paste0(maps_sfs_names[1,])); region_buckets <- unique(paste0(maps_sfs_names[2,]))
    #maps_cols <- c("black", rev(rainbow(nrow(maps_sfs)-1)))
    maps_x <- c(sapply(1:length(s_buckets), function(x) rep(x, length(region_buckets))))
    maps_cols <- rep(c("black",rainbow(length(region_buckets))), length(s_buckets))
    plot(maps_x, maps_sfs$MAPS, col=maps_cols, main="MAPS vs. Predicted Selection Coef.", xlab="", ylab="MAPS", ylim=range(c(maps_sfs[,c("MAPS_lower","MAPS_upper")])), xaxt="n", cex.axis=1.3, cex.lab=1.3, cex.main=1.3)
    abline(h=0, col="gray50", lty=3)
    #segments(x0=maps_x, y0=maps_sfs$MAPS_lower, x1=maps_x, y1=maps_sfs$MAPS_upper, col=maps_cols)
    segments(x0=maps_x-0.05, y0=maps_sfs$MAPS_lower, x1=maps_x+0.05, y1=maps_sfs$MAPS_lower, col=maps_cols, lty=2)
    segments(x0=maps_x-0.05, y0=maps_sfs$MAPS_upper, x1=maps_x+0.05, y1=maps_sfs$MAPS_upper, col=maps_cols, lty=1)
    text(1:length(s_buckets)+0.15, par("usr")[3]-0.0012, srt=45, adj=1, xpd=TRUE, labels=s_buckets, cex=1)
    legend("bottomleft", legend=c(region_buckets), col=maps_cols[1:length(region_buckets)], pch=15, cex=1)
    mtext(paste0(sites," sites"), cex=1.2)
    dev.off()
    pdf_to_png(filename)
    
    maps_sfs_full <- maps_sfs; rownames(maps_sfs_full) <- c(sapply(s_buckets, function(x) paste0(x,", ",unique(maps_sfs_full$region_type))))
    maps_sfs <- maps_sfs_full[maps_sfs_full$region_type == "all",1:6]
    rownames(maps_sfs) <- gsub(", all", "", rownames(maps_sfs_full)[maps_sfs_full$region_type == "all"])
    filename = output_path(paste0("model_sfs_",sites,"_logscale.pdf"))
    pdf(file=filename)
    s_cols <- rainbow(nrow(maps_sfs))
    #c(0,ceiling(max(maps_sfs_full$PS_observed*10))*10)
    plot(log10(1:ncol(maps_sfs)), log10(maps_sfs[1,]/sum(maps_sfs[1,])), col=s_cols[1], type="o", main="Site Frequency Spectrum", xlab="", ylab="log10(Proportion freq)", xaxt="n", xlim=log10(c(1,ncol(maps_sfs))), cex.axis=1.3, cex.lab=1.3, cex.main=1.3)
    for(maps_i in 1:ncol(maps_sfs)) { abline(v=log10(maps_i), lty=3, col="gray50") }
    for(s_i in 2:nrow(maps_sfs)) { lines(log10(1:ncol(maps_sfs)), log10(maps_sfs[s_i,]/sum(maps_sfs[s_i,])), col=s_cols[s_i], type="o") }
    #if(sites == "CpG") { mtext_text = paste0(sum(maps_cpg_vector)," supermodel test sequence CpG sites")
    #} else if (sites == "non-CpG") { mtext_text = paste0(sum(!maps_cpg_vector)," supermodel test sequence non-CpG sites")
    #} else { mtext_text = paste0("all"," supermodel test sequence sites") }
    mtext_text = paste0(sites," sites")
    mtext(mtext_text, cex=1.1)
    legend("topright", legend=c("Predicted selection coef. s", paste0(rownames(maps_sfs)," (N = ",rowSums(maps_sfs[,1:6]),")")), col=c("white",s_cols), pch=15, cex=1)
    text(log10(1:ncol(maps_sfs))+0.04, par("usr")[3]-0.01, srt=45, adj=1, xpd=TRUE, labels=colnames(maps_sfs), cex=1.3)
    dev.off()
    pdf_to_png(filename)
    
    maps_sfs <- aggregate(maps_sfs_full[,1:6], by=list(gsub("^.*, ","",rownames(maps_sfs_full))), FUN=sum) #grepl("s = any",rownames(maps_sfs_full)),1:6]
    rownames(maps_sfs) <- maps_sfs$Group.1; maps_sfs <- maps_sfs[,-1]
    colnames(maps_sfs) <- c("Singletons", "Doubletons", "Tripletons", "AC = 4", "AC = 5", "AC > 5")
    filename = output_path(paste0("model_sfs_regional_",sites,"_logscale.pdf"))
    pdf(file=filename)
    s_cols <- rainbow(nrow(maps_sfs))
    #c(0,ceiling(max(maps_sfs_full$PS_observed*10))*10)
    plot(log10(1:ncol(maps_sfs)), log10(maps_sfs[1,]/sum(maps_sfs[1,])), col=s_cols[1], type="o", main="Site Frequency Spectrum", xlab="", ylab="log10(Proportion freq)", xaxt="n", xlim=log10(c(1,ncol(maps_sfs))), cex.axis=1.3, cex.lab=1.3, cex.main=1.3)
    for(maps_i in 1:ncol(maps_sfs)) { abline(v=log10(maps_i), lty=3, col="gray50") }
    for(s_i in 2:nrow(maps_sfs)) { lines(log10(1:ncol(maps_sfs)), log10(maps_sfs[s_i,]/sum(maps_sfs[s_i,])), col=s_cols[s_i], type="o") }
    #if(sites == "CpG") { mtext_text = paste0(sum(maps_cpg_vector)," supermodel test sequence CpG sites")
    #} else if (sites == "non-CpG") { mtext_text = paste0(sum(!maps_cpg_vector)," supermodel test sequence non-CpG sites")
    #} else { mtext_text = paste0("all"," supermodel test sequence sites") }
    mtext_text = paste0(sites," sites")
    mtext(mtext_text, cex=1.1)
    legend("topright", legend=c("Region type", paste0(rownames(maps_sfs)," (N = ",rowSums(maps_sfs[,1:6]),")")), col=c("white",s_cols), pch=15, cex=1)
    text(log10(1:ncol(maps_sfs))+0.04, par("usr")[3]-0.01, srt=45, adj=1, xpd=TRUE, labels=colnames(maps_sfs), cex=1.3)
    dev.off()
    pdf_to_png(filename)
    
    #plot(maps_sfs_full$MAPS[maps_sfs_full$region_type == "three_prime_UTR"], maps_sfs_full$s[maps_sfs_full$region_type == "three_prime_UTR"])
}
maps_sfs <- read.csv(output_path("model_sfs_all_table.csv"))

# Find correlation of s with obs/exp
o_e_z <- rep(readRDS(output_path(paste0(dat_name,"_obs_exp.rds")))[,1], 151)
s <- c(1/s_vals_all)
cor(s, o_e_z)
plot(o_e_z, s, main="", xlab="Obs/Exp", ylab="s", cex.axis=1.4, cex.lab=1.4, cex.main=1.3)
draw_plot(data.frame(x=o_e_z, y=s), title="Predicted s vs. obs/exp input", xlab="obs/exp input", ylab="s", legend_text="# sites", cor_method="Spearman", linear_best_fit=FALSE, quadratic_best_fit=FALSE, filename=output_path(paste0("model_s_vs_oez.pdf")))

dev.off()
pdf_to_png(filename)
filename = output_path(paste0("model_s_vs_oez_logscale.pdf"))
pdf(file=filename)
plot(o_e_z, log10(s), main="", xlab="Obs/Exp", ylab="log10(s)", cex.axis=1.4, cex.lab=1.4, cex.main=1.3)
dev.off()
pdf_to_png(filename)
    

# Plot regional s comparison
# ADD x-labels that are slanted and equal to region_types
# YOU LEFT OFF HERE!!!!
filename = output_path(paste0("regional_s_comparison.pdf"))
pdf(file=filename)
plot(1:length(region_types), s_vals_all_means$estimate, col=cols, main="Regional Selection Coef. Distributions", xlab="", ylab="Selection Coef. s", ylim=range(c(s_vals_all_means[,c("conf.int_lower","conf.int_higher")])), xaxt="n", cex.axis=1.3, cex.lab=1.3, cex.main=1.3)
#abline(h=0, col="gray50", lty=3)
segments(x0=1:length(region_types), y0=s_vals_all_means$conf.int_lower, x1=1:length(region_types), y1=s_vals_all_means$conf.int_higher, col=cols)
text(1:length(region_types)+0.2, par("usr")[3]-0.001, srt=45, adj=1, xpd=TRUE, labels=region_types, cex=0.8)
dev.off()
pdf_to_png(filename)

cor(maps_result$AC, maps_result$mutability, method="spearman")
cor(maps_result$AC, maps_result$s, method="spearman")
draw_plot(data.frame(x=maps_result$AC[maps_result$AC > 0] == 1, y=maps_result$s[maps_result$AC > 0], hex_density = 25, title="PS vs. Mutability", xlab="AC", ylab="s", legend_text="# regions", linear_best_fit=FALSE, quadratic_best_fit=FALSE)) #, filename="PS_vs_mutability.pdf")
range(maps_result$s[maps_result$AC == 1]); range(maps_result$s[maps_result$AC > 1]) 
plot(density(maps_result$s[maps_result$AC == 1 & maps_result$s > 0]), col="red")
lines(density(maps_result$s[maps_result$AC > 1 & maps_result$s > 0]), col="blue")

maps_model <- readRDS(output_path("maps_model.rds"))
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
    text(1:nrow(maps_sfs)+0.2, par("usr")[3]-0.0002, srt=45, adj=1, xpd=TRUE, labels=rownames(maps_sfs), cex=1)
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
a <- get_gradcams(output_path("regional2_gradcams"))

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
        
        
        
# Investigate filter motifs and respective activations in a given dataset.
# //TODO: Either clean up or delete this function.
investigate_filter_motifs <- function() {
    layers_to_investigate = c("conv1", "maxpool1")
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
# OTHER FUNCTIONS
#################################################################################################################













