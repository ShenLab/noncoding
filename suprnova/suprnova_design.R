library(keras)
library(kerasR)
library("tfprobability")
gpu <- tf$config$experimental$list_physical_devices("GPU")[[1]]
tf$config$experimental$set_memory_growth(gpu, TRUE) # Important step on some servers to uncap the GPU memory usage right away. 
library(stringr)
library(pbapply)
library(EBImage)
library(rhdf5)
library(gtools)
#library("PWMEnrich")
library("TFBSTools")
#library("seqLogo")
#library("msa")
library("ComplexHeatmap")
library("abind")
library("hexbin")
library("RColorBrewer")
library("ggplot2")
library("tidyverse")
source("alex_suite.R")
# //TODO: Not all of the above are dependencies any more, and will be cleaned up.

#################################################################################################################
# PARAMETERS
#################################################################################################################
Ne = 10000 # Effective population size for humans; can try between 10,000 to 20,000.
s_scaling_factor = -1 # The neural network represents s in more natural scale and then converts it to real s with exp(s * s_scaling_factor); this significantly helps smooth/improve the training.

#################################################################################################################
# SIMPLE FUNCTIONS
#################################################################################################################
# Get all layer names of the given model
get_model_layer_names <- function(model) {
    model_layer_names <- sapply(model$layers, function(layer) paste0(layer$name))
    return(model_layer_names)
}

# Return an activation model of the given model, with outputs specified at every layer_to_investigate listed.
get_activation_model <- function(model, layers_to_investigate = c("conv1","maxpool1","d","s")) {
    model_layer <- model$get_layer(layers_to_investigate[1])
    model_layer_names <- get_model_layer_names(model)
    layers_to_investigate <- layers_to_investigate[layers_to_investigate %in% model_layer_names]
    activation_model <- keras_model(model$input, outputs=lapply(layers_to_investigate, function(layer_to_investigate) { return(get_output_at(model$get_layer(layer_to_investigate), 1)) })) #sapply(model$layers[2:3], function(layer) return(get_output_at(layer[[1]], 1)))) #layer[[1]]$output)
    return(activation_model)
}

# Saves the model (weights) under a given model name.
save_supermodel <- function(model, model_name) {
    filename = output_path(paste0(model_name,"_weights.hdf5"))
    save_model_weights_hdf5(model, filename)
}

# Loads the model with the given model name; because of custom architecture, this requires first re-defining the architecture in keras and then loading in the saved weights.
load_supermodel <- function(model_name, type="v", L=15, num_rbps=160, num_filters=5, rbp_kernel_width=4, sequence_kernel_width=6, s_anchor=2e-6, s_scaling_factor=-1) {
    filename = output_path(paste0(model_name,"_weights.hdf5"))
    model <- define_supermodel_architecture(type, L, num_rbps, num_filters, rbp_kernel_width, sequence_kernel_width, s_anchor, s_scaling_factor) # HAVE TO FIGURE OUT A WAY TO SAVE THE TRANSFER THE PARAMETERS AS WELL AND LOAD THEM FIRST!
    model <- load_model_weights_hdf5(model, filename)
    return(model)
}

#################################################################################################################
# Define the super-model architecture, given the parameters. 
# Type can be "r"/"regional" or "v"/"variant" for rSUPRNOVA and vSUPRNOVA, respectively.
#################################################################################################################
# type="v"; L=15; num_rbps=160; num_filters=5; rbp_kernel_width=4; sequence_kernel_width=6; s_scaling_factor=-1
define_supermodel_architecture <- function(type="v", L=15, num_rbps=160, num_filters=5, rbp_kernel_width=3, sequence_kernel_width=3, s_anchor=2e-6, s_scaling_factor=-1, custom_conv1=FALSE) {
    preserve_variants = (type == "v" || type == "variant")
    if(preserve_variants) {
        # Define the shape of the RBP alt-ref GradCAM tensor input to the neural network.
        delta_binding_input <- layer_input(name="delta_binding_input", shape = c(L,num_rbps,4))
        # Define the shape of the RBP GradCAM tensor input to the neural network.
        rbp_binding_input <- layer_input(name="rbp_binding_input", shape = c(L,num_rbps,1))
        # alt-ref convolution layer
        #w <- model$get_layer("conv_alt_ref")$weights[[1]][1,1,1,,]
        conv_altref_constraint <- function(w) { k_minimum(k_maximum(w, -1), 1) }
        conv_altref <- layer_conv_3d(name="conv_alt_ref", filters=1, kernel_size=c(1,1,1), strides=c(1,1,1), padding="valid", use_bias=FALSE) #kernel_constraint=conv_altref_constraint
        # Custom conv1 filter weights, to account for every RBP complex possibility
        if(custom_conv1) {
            row_weights <- expand.grid(lapply(1:rbp_kernel_width, function(i) c(0,1)))
            row_weights_combs <- expand.grid(lapply(1:sequence_kernel_width, function(i) 1:nrow(row_weights)))
            conv1_filters <- abind(lapply(1:nrow(row_weights_combs), function(i) { f <- as.matrix(row_weights[unlist(row_weights_combs[i,]),]); colnames(f) <- gsub("Var","RBP",colnames(f)); return(f) }), along=3)
            num_filters = dim(conv1_filters)[3]
            midpoint_RBP <- ceiling(rbp_kernel_width / 2)
            midpoint_seq <- ceiling(sequence_kernel_width / 2)
            conv1_filters <- conv1_filters[,,unlist(lapply(1:num_filters, function(i) { # 25 filters excluded when filters are 3x3 dimension
                f <- conv1_filters[,,i,drop=FALSE]; rbp_inv <- colSums(f); seq_inv <- rowSums(f)
                return(sum(!unlist(lapply(list(rbp_inv, seq_inv), function(x) {
                    midpoint = ceiling(length(x)/2)
                    return(x[midpoint] > 0 | ((sum(x[1:(midpoint-1)]) > 0) & (sum(x[(midpoint+1):rbp_kernel_width]) > 0)))
                }))) == 0)
            })),drop=FALSE]; conv1_filters <- conv1_filters[,,order(apply(conv1_filters, 3, sum)),drop=FALSE]
            num_filters = dim(conv1_filters)[3]
            conv1_filters <- k_cast(conv1_filters[,,1:num_filters,drop=FALSE], dtype="float32")
            my_filter <- function(shape=c(sequence_kernel_width,rbp_kernel_width,1,1,num_filters), dtype=NULL) { k_reshape(conv1_filters, shape) }
            # First convolution layer
            conv1 <- layer_conv_3d(name="rbp_complex_activation", filters=num_filters, kernel_size=c(1,sequence_kernel_width,rbp_kernel_width), strides=c(1,1,rbp_kernel_width-2), padding="same", use_bias=FALSE,
                                   kernel_initializer=my_filter, trainable=FALSE, activation="tanh")
        } else {
            # First convolution layer
            conv1 <- layer_conv_3d(name="rbp_complex_activation", filters=num_filters, kernel_size=c(1,sequence_kernel_width,rbp_kernel_width), strides=c(1,1,rbp_kernel_width-2), padding="same", use_bias=FALSE, activation="tanh")
        }
        #gradcam_reformat_layer <- layer_lambda(name="gradcam_reformat", f = function(inputs) {
        #    delta_gradcams <- inputs[[1]]; ref_gradcams <- inputs[[2]]
        #    return(k_stack(list(delta_gradcams, ref_gradcams), axis=5))
        #}) (c(delta_binding_input, rbp_binding_activation))
        
        #scaling_layer_function <- function(inputs) {}
        #scaling_layer <- layer_lambda(name="s", f = function(inputs) {
        #    s <- inputs[[1]]; #s <- s_constraint(s)
        #    s <- k_exp(s * s_scaling_factor + s_anchor)
        #    return(s)
        #}) (c(selection_coef_layer))
        
        # Custom Keras layer to weight the final dimension of the input by trainable weightings.
        # If the act parameter is set to "threshold", a compressed sigmoid is done around the kernel to push all weights towards 0 or 1.  
        TrainableWeightingLayer <- R6::R6Class("CustomLayer",      
            inherit = KerasLayer,
            public = list(
                output_dim = NULL,
                kernel = NULL,
                threshold = NULL,
                ndims = 2,
                along = 0,
                anchor = 0.5,
                initialize = function(threshold, along, anchor, ndims) { self$ndims <- ndims; self$threshold <- threshold; self$along <- along; self$anchor <- anchor }, #self$output_dim <- output_dim
                build = function(input_shape) {
                    self$ndims <- length(input_shape)
                    if(self$along < 1 || self$along > self$ndims) { self$along <- self$ndims }
                    self$kernel <- self$add_weight(
                        name = "kernel",
                        shape = lapply(2:self$ndims, function(i) { if(i == self$along) { return(input_shape[[i]]) } else { return(k_cast(1,dtype="int32")) }}),
                        #shape = list(input_shape[[self$along]]), #list(input_shape[[2]], self$output_dim),
                        #initializer = initializer_truncated_normal(mean=self$anchor, stddev=(1-self$anchor)/2)
                        initializer = initializer_truncated_normal(mean=self$anchor, stddev=min(c(1-self$anchor, self$anchor))/2), # so that both 0 and 1 are 2 standard devs away and will be redrawn.
                    )
                },
                call = function(x, mask = NULL) {
                    if(is.null(self$threshold)) { return(x * self$kernel)
                    } else {
                        if(!is.numeric(self$threshold)) { self$threshold <- 5/(1-self$anchor) }
                        return(k_sigmoid(self$threshold*(x - self$kernel)))
                    }
                },
                compute_output_shape = function(input_shape) { input_shape }
            )
          )
        ## define layer wrapper function
        layer_weighting <- function(object, name = NULL, along=0, threshold=NULL, anchor=0.5, trainable=TRUE) {
            create_layer(TrainableWeightingLayer, object, list(threshold=threshold, name=name, along=along, anchor=anchor, trainable=trainable, ndims=2))
        }
        #layer_weighting(name="delta_binding_activation")
        #layer_custom <- function(object, output_dim, name = NULL, trainable = TRUE) { create_layer(CustomLayer, object, list(output_dim=as.integer(unlist(output_dim)), name = name, trainable = trainable)) }
        #delta_binding_activation_layer <- layer_weighting(name="delta_binding_activation", anchor=0.4, along=2)
        #delta_binding_activation <- abs(w) %>% delta_binding_activation_layer
        #delta_binding_activation
        #mean(as.matrix(delta_binding_activation)); sd(as.matrix(delta_binding_activation))
        #delta_binding_activation_layer$weights
        #delta_binding_activation_layer <- layer_weighting(name="delta_binding_activation", threshold=TRUE)
        #delta_binding_activation <- abs(w) %>% delta_binding_activation_layer
        #delta_binding_activation
        #mean(as.matrix(delta_binding_activation)); sd(as.matrix(delta_binding_activation))
        #delta_binding_activation_layer$weights
        #delta_binding_activation_layer <- layer_weighting(name="delta_binding_activation", threshold=100)
        #delta_binding_activation <- abs(w) %>% delta_binding_activation_layer
        #delta_binding_activation
        #mean(as.matrix(delta_binding_activation)); sd(as.matrix(delta_binding_activation))
        #delta_binding_activation_layer$weights
        #delta_binding_activation_layer$weights
        #w
        #delta_binding_activation
        #delta_binding_activation_layer <- layer_weighting(name="delta_binding_activation", along=3, threshold=TRUE, anchor=0.75)
        #delta_binding_activation <- delta_binding_input %>% delta_binding_activation_layer
        #delta_binding_activation_layer$weights
        
        # Module responsible for processing tensor with CNN and learning relevant RBP binding patterns.
        #conv_delta <- layer_conv_2d(name="conv_delta", filters=1, kernel_size=c(1,1), strides=c(1,1), padding="valid", activation="sigmoid", use_bias=FALSE) #kernel_constraint=conv_altref_constraint
        #delta_binding_activation <- delta_binding_input %>% layer_permute(c(1,3,2)) %>% conv_delta %>% layer_permute(name="delta_binding_activation",c(1,3,2))
        #conv_rbp_binding <- layer_conv_2d(name="conv_ref_binding", filters=1, kernel_size=c(1,1), strides=c(1,1), padding="valid", activation="sigmoid", use_bias=FALSE) 
        rbp_activation_layer <- layer_weighting(name="rbp_binding_activation", along=3, threshold=TRUE, anchor=0.75)
        delta_binding_activation <- delta_binding_input %>% rbp_activation_layer
        rbp_binding_activation <- rbp_binding_input %>% rbp_activation_layer
        rbp_disruption_module <- layer_subtract(c(delta_binding_activation, rbp_binding_activation), name="altref_binding_changes") %>%
            layer_reshape(c(unlist(delta_binding_input$shape[c(4,2,3)]),1)) %>% #c(prod(unlist(rbp_binding_input$shape[2:3])),unlist(rbp_binding_input$shape[4:5]))) %>%
            #layer_dense(1, name="rbp_binding_dense", activation="linear", use_bias=FALSE) %>% #conv_alt_ref %>% layer_dense
            #layer_reshape(c(unlist(rbp_binding_input$shape[2:4]),1)) %>%      kernel_constraint=constraint_nonneg()
            conv1 %>% layer_spatial_dropout_3d(0.2, name="dropout_weak1")
            
        # Select the middle site; prior to this a conv kernel across positions can be done to pool nearby site information, but from now on we focus only on central one.
        site_select_layer <- layer_lambda(name="central_site_select", f = function(inputs) {
            x <- inputs[[1]][,,ceiling(L/2),,]
            return(x %>% layer_reshape(c(x$shape[[2]],prod(unlist(x$shape[3:4])))))
        }) (c(rbp_disruption_module))
        
        #layer_reshape(c(unlist(rbp_disruption_module$shape[2:3]),prod(unlist(rbp_disruption_module$shape[4:5])))) %>%
        #layer_max_pooling_3d(name="maxpool1", pool_size=c(1,2,1)) %>% conv2
        
        # Second convolution layer
        #conv2 <- layer_conv_3d(name="conv2", filters=num_filters, kernel_size=c(1,1,rbp_kernel_width), strides=c(1,1,rbp_kernel_width-2), padding="same", use_bias=FALSE, activation="relu")                  
        # Pool together all RBP disruptions into one final aggregate RBP disruption representation.
        rbp_disruption_module_output <- site_select_layer %>% layer_dense(1, name="gene_regulatory_disruption", activation="sigmoid", use_bias=TRUE)
            #layer_reshape(unlist(rbp_disruption_module$shape[2:4])) %>% layer_permute(c(1,3,2,4)) %>% 
            #layer_max_pooling_3d(name="maxpool_final", pool_size=c(rbp_disruption_module$shape[[3]],rbp_disruption_module$shape[[2]],num_filters)) %>% 
            #layer_flatten(name="rbp_disruption_output")
        
    } else {
        # Define the shape of the RBP GradCAM tensor input to the neural network.
        rbp_binding_input <- layer_input(name="rbp_binding_input", shape = c(L,num_rbps,1,1))
        # First convolution layer
        conv1 <- layer_conv_2d(name="rbp_complex_activation", filters=num_filters, kernel_size=c(sequence_kernel_width,rbp_kernel_width), strides=c(1,rbp_kernel_width-2), padding="valid", use_bias=FALSE, activation="relu")
        # Second convolution layer
        conv2 <- layer_conv_2d(name="conv2", filters=num_filters, kernel_size=c(sequence_kernel_width,rbp_kernel_width), strides=c(1,rbp_kernel_width-2), padding="valid", use_bias=FALSE, activation="relu")
    
        # Module responsible for processing tensor with CNN and learning relevant RBP binding patterns.
        rbp_disruption_module <- rbp_binding_input %>% layer_reshape(c(unlist(rbp_binding_input$shape[2:4]))) %>%
            conv1 %>% layer_spatial_dropout_2d(0.2, name="dropout_weak1") %>% 
            layer_max_pooling_2d(name="maxpool1", pool_size=c(1,2)) %>% conv2
        
        # Pool together all RBP disruptions into one final aggregate RBP disruption representation.
        rbp_disruption_module_filters_pooled <- rbp_disruption_module %>% layer_reshape(c(unlist(rbp_disruption_module$shape[2:4]),1)) %>% 
            layer_max_pooling_3d(name="maxpool_final", pool_size=c(1,rbp_disruption_module$shape[[3]],num_filters)) %>%
            layer_flatten() 
        # Upsample compressed RBP disruption encoding back up to the input sequence length, and outputs it as a flattened vector.
        rbp_disruption_module_output <- rbp_disruption_module_filters_pooled %>% 
            layer_reshape(c(rbp_disruption_module_filters_pooled$shape[[2]],1)) %>% layer_upsampling_1d(L) %>% 
            layer_average_pooling_1d(rbp_disruption_module_filters_pooled$shape[[2]]) %>%
            layer_flatten(name="gene_regulatory_disruption")
    }
    
    # Adjust gene damagingness layer to follow constraints on d.
    # Depending on how we want the representation of d to be constrained, we can adjust the following bounds:
    #d_lower_bound = -Inf; d_upper_bound = Inf
    #d_constraint <- function(w) { k_minimum(k_maximum(w, d_lower_bound), d_upper_bound) }
    #adjusted_gene_damagingness_layer <- layer_lambda(name="d", f = function(inputs) {
    #    d <- inputs[[1]]; #d <- d_constraint(d)
    #    return(d)
    #}) (c(rbp_disruption_module_output))

    adjusted_gene_damagingness_layer <- rbp_disruption_module_output %>% layer_reshape(c(rbp_disruption_module_output$shape[[2]],1)) %>%
        #layer_dense(1, name="d_dense1", activation="linear", use_bias=TRUE) %>%
        #layer_dense(1, name="d_dense2", activation="tanh", use_bias=TRUE) %>%
        #layer_dense(1, name="d_dense3", activation="tanh", use_bias=TRUE) %>%
        #layer_dense(1, name="d_dense4", activation="tanh", use_bias=TRUE) %>%
        #layer_dense(1, name="d_dense5", activation="sigmoid", use_bias=TRUE) %>%
        layer_flatten(name="d")
    
    # Add nearest gene and epigenomic features to the neural network, which are combined with gene damagingness d to make s!
    num_tissues = 7
    # Auxiliary input for gene expression
    gene_expression_module <- layer_input(name="gene_expression", shape = c(num_tissues)) 
    gene_expression_module_activation <- gene_expression_module %>% layer_dense(1, name="gene_expression_activation", activation="sigmoid")
    # Auxiliary input for observed/expected gene constraint
    gene_constraint_module <- layer_input(name="gene_constraint", shape = c(1))
    gene_constraint_module_activation <- gene_constraint_module %>% layer_dense(1, name="gene_constraint_activation", activation="sigmoid", use_bias=FALSE) #%>%
    # Auxiliary input for histone/epigenomic marks annotation
    epigenomic_module <- layer_input(name="epigenomic_marks", shape = c(num_tissues,1)) #adjusted_gene_damagingness_layer$shape[[2]]))
    epigenomic_module_activation <- epigenomic_module %>% layer_dense(1, name="epigenomic_site_aggregate", activation="sigmoid", use_bias=FALSE) %>%
        layer_flatten() %>% layer_dense(1, name="epigenomic_activation", activation="sigmoid", use_bias=FALSE) 
    
    # Learning of selection coefficient using RBP gene regulation disruption output and gene-level features.
    # layer_multiply(c(adjusted_gene_damagingness_layer, gene_expression_module_activation, gene_constraint_module_activation, epigenomic_module_activation)) %>% # %>% #gene_damagingness_stack_layer %>% 
    selection_coef_layer <- #layer_multiply(c(adjusted_gene_damagingness_layer, gene_expression_module_activation, gene_constraint_module_activation, epigenomic_module_activation)) %>% 
        layer_multiply(c(adjusted_gene_damagingness_layer, gene_constraint_module_activation)) %>%
        layer_reshape(c(adjusted_gene_damagingness_layer$shape[[2]],1)) %>%
        layer_dense(1, name="dense1", activation="tanh", use_bias=TRUE) %>% #sigmoid
        layer_dense(1, name="dense2", activation="tanh", use_bias=TRUE) %>% #layer_activation_leaky_relu(alpha=0.3) %>%
        layer_dense(1, name="dense3", activation="tanh", use_bias=TRUE) %>% #layer_activation_leaky_relu(alpha=0.3) %>%
        layer_dense(1, name="dense4", activation="tanh", use_bias=TRUE) %>%
        #layer_dense(1, name="dense5", activation="tanh", use_bias=TRUE) %>%
        #layer_dense(1, name="dense6", activation="tanh", use_bias=TRUE) %>%
        #layer_dense(1, name="dense7", activation="tanh", use_bias=TRUE) %>%
        #layer_dense(1, name="dense8", activation="tanh", use_bias=TRUE) %>%
        #layer_dense(1, name="dense9", activation="tanh", use_bias=TRUE) %>%
        #layer_dense(1, name="dense10", activation="tanh", use_bias=TRUE) %>%
        layer_dense(1, name="dense_final", activation="tanh", use_bias=TRUE) %>%
        layer_flatten(name="selection_coef_output")
    
    # Adjust selection coefficient layer to follow constraints on s.
    # exp(s * s_scaling_factor) has to be [0,1] in the end, so there need to be according bounds for the neural net's representation of s. 
    if(s_scaling_factor < 0) { s_lower_bound = 0; s_upper_bound = Inf } else { s_upper_bound = 0; s_lower_bound = -Inf }
    s_constraint <- function(w) { k_minimum(k_maximum(w, s_lower_bound), s_upper_bound) }
    adjusted_selection_coef_layer <- layer_lambda(name="s", f = function(inputs) {
        s <- inputs[[1]]; #s <- s_constraint(s)
        s <- k_exp(s * s_scaling_factor + s_anchor)
        return(s)
    }) (c(selection_coef_layer))
    
    # Auxiliary input for background mutation rate
    background_mut_rate_module <- layer_input(name="background_mut_rate", shape = c(adjusted_gene_damagingness_layer$shape[[2]]))
    
    # Auxiliary input for y_true (observed allele frequency)
    y_true_module <- layer_input(name="y_true", shape = c(adjusted_gene_damagingness_layer$shape[[2]]))
    
    # Define the Negative Binomial Loss layer, which calculates the negative logloss of the fit of our predicted s given the observed allele count (y_true * sample_size) and background mutation rate mu.  
    # Old Poisson method: log_probs <- k_mean(y_true) %>% (tfp$distributions$Poisson(rate=(k_mean(y_pred)))$log_prob)
    # New approach with Negative Binomial (allows overdispersion)
    # Probs p = 1 / (z + 1), where z = 4 * Ne * s / N 
    # Number of failures r = theta = 4 * Ne * mu
    # "Success" rate AF is the variable i term fed into the NB
    # Mean of NB will be pr/(1-p) = ( 1 * 4 * Ne * mu / (z+1) ) / ( z/(z+1) ) = sample_size * mu / s
    loss_layer <- layer_lambda(name="neg_log_probs", f = function(inputs) { 
        y_true <- inputs[[1]]; mu <- inputs[[2]]; s <- inputs[[3]]
        z = 4 * Ne * s / sample_size
        r = 4 * Ne * mu
        neg_log_probs <- - ((y_true * sample_size) %>% (tfp$distributions$NegativeBinomial(total_count=r, probs=(1/(z+1)))$log_prob)) #y_true * sample_size
        return(neg_log_probs)
    }) (c(y_true_module, background_mut_rate_module, adjusted_selection_coef_layer))
    
    # Define model with keras functional API
    model <- keras_model(inputs=c(delta_binding_input, rbp_binding_input, gene_constraint_module, background_mut_rate_module, y_true_module), outputs=c(loss_layer))
    #model <- keras_model(inputs=c(rbp_binding_input, gene_expression_module, gene_constraint_module, epigenomic_module, background_mut_rate_module, y_true_module), outputs=c(loss_layer))
    
    # Custom loss function, that can optionally take in a mask arg and additional extra args.
    custom_loss_function <- function(extra_args=NULL, mask=-1) {
        loss_function <- function(y_true, y_pred) { 
            # Focus on plausible de novos (y_true < 1e-3) and plausible constrained / non-random sites (0 < neg_log_probs < Inf)
            eligible_indices <- tf$where(k_flatten(k_greater_equal(y_true, 0) & k_less(y_true, 1e-3) & k_greater(y_pred, 0) & k_less(y_pred, Inf)))
            y_true <- k_gather(k_flatten(y_true), eligible_indices)
            y_pred <- k_gather(k_flatten(y_pred), eligible_indices)
            # Calculate total logloss of all eligible sites taken together in this training iteration.
            sum_neg_log_probs <- k_mean(y_pred) # Use mean instead of sum to control for different number of eligible (mu available) sites
            return((sum_neg_log_probs-1)*100) # Can print the result out before this (during runtime) with k_print_tensor(sum_neg_log_probs)
        }
        return(loss_function)
    }
    opt <- optimizer_adam(lr=0.01) # optimizer_rmsprop(lr=0.01)
    model %>% compile(loss=custom_loss_function(), optimizer=opt) #, metrics="mae")
    
    model_name = "supermodel"
    # plot_model(model, to_file = output_path(paste0("images/",model_name,".png")), show_shapes = TRUE, show_layer_names = TRUE)
    # Write model architecture to file, and return the defined model.
    write(paste0(model), file=output_path("model_architecture.txt"))
    return(model)
}

#################################################################################################################
# Load all of the diffrent tensor objects (for the given named dataset), required for SUPRNOVA training.
#################################################################################################################
load_supermodel_training_data <- function(name="gnomad", type="v", remove_zero_expected=TRUE, remove_common_variants=TRUE) {
    version="hg19"; set_global(version)
    # gnomAD sample size, used for converting allele count to allele frequency
    sample_size = 71702; set_global(sample_size)
    expected <- readRDS(output_path(paste0(name,"_expected.rds")))
    if(remove_zero_expected) { to_keep <- apply(expected, 1, function(x) sum(x==0) == 0) } else { to_keep <- rep(TRUE, nrow(expected)) }
    # Load all input and output tensors, which were calculated and saved in the previous section.
    if(type=="v") {
        midpoint = ceiling(ncol(expected)/2); arm_width=7; site_indices <- (midpoint-arm_width):(midpoint+arm_width)
        af_tensor <- readRDS(output_path(paste0(name,"_AFs_variant.rds")))[,midpoint,,drop=TRUE] * 2
        if(remove_common_variants) { to_keep <- to_keep & apply(af_tensor, 1, function(x) (sum(x >= 1e-3) == 0 & sum(x < 0) == 1)) }
        af_tensor <- af_tensor[to_keep,,drop=FALSE]
        dat_cpg <- readRDS(output_path(paste0(name,"_cpg.rds")))[to_keep,,drop=FALSE]; set_global(dat_cpg)
        
        expected <- expected[to_keep,rep(midpoint,4),drop=FALSE]
        cpg_sites <- dat_cpg[,midpoint]
        expected[cpg_sites,] <- expected[cpg_sites,] * c(1,1,1,10)/12
        expected[!cpg_sites,] <- expected[!cpg_sites,]/3
        dat_gradcams <- readRDS(output_path(paste0(name,"_gradcams.rds")))[to_keep,site_indices,,,drop=FALSE]
        dat_A_gradcams <- readRDS(output_path(paste0(name,"_A_gradcams.rds")))[to_keep,site_indices,,,drop=FALSE]
        dat_C_gradcams <- readRDS(output_path(paste0(name,"_C_gradcams.rds")))[to_keep,site_indices,,,drop=FALSE]
        dat_G_gradcams <- readRDS(output_path(paste0(name,"_G_gradcams.rds")))[to_keep,site_indices,,,drop=FALSE]
        dat_T_gradcams <- readRDS(output_path(paste0(name,"_T_gradcams.rds")))[to_keep,site_indices,,,drop=FALSE]
        #dat_gradcams <- abind(list(abind(list(dat_A_gradcams, dat_C_gradcams, dat_G_gradcams, dat_T_gradcams), along=4),
        #                            abind(list(dat_gradcams, dat_gradcams, dat_gradcams, dat_gradcams), along=4)), along=5); set_global(dat_gradcams)
        #dat_gradcams <- aperm(dat_gradcams, c(1,2,3,5,4))
        altref_gradcams <- abind(lapply(list(dat_A_gradcams, dat_C_gradcams, dat_G_gradcams, dat_T_gradcams), function(alt_gradcams) {
            return(alt_gradcams) # - dat_gradcams)
        }), along=4); set_global(altref_gradcams)
    } else {
        dat_gradcams <- abind(readRDS(output_path(paste0(name,"_gradcams.rds")))[to_keep,,,,drop=FALSE],5)
        af_tensor <- readRDS(output_path(paste0(name,"_AFs.rds")))[to_keep,,drop=FALSE] * 2
        dat_cpg <- readRDS(output_path(paste0(name,"_cpg.rds")))[to_keep,,drop=FALSE]; set_global(dat_cpg)
    }
    set_global(af_tensor); set_global(dat_gradcams); set_global(expected)
    
    dat_trimers <- readRDS(output_path(paste0(name,"_trimers.rds")))[to_keep,,drop=FALSE]; set_global(dat_trimers)
    dat_region_types <- readRDS(output_path(paste0(name,"_region_types.rds")))[to_keep,,drop=FALSE]; set_global(dat_region_types)
    dat_pLIs <- readRDS(output_path(paste0(name,"_pLIs.rds")))[to_keep,,drop=FALSE]; set_global(dat_pLIs)
    dat_obs_exp <- readRDS(output_path(paste0(name,"_obs_exp.rds")))[to_keep,,drop=FALSE]; set_global(dat_obs_exp)
    width=ncol(dat_gradcams); set_global(width)
    arm_width=floor(width/2); set_global(arm_width)
    if(name == "gnomad") { dat <- readRDS(output_path("gnomad_fully_unpacked_annotated.rds")); dat_midpoints <- seq((arm_width+1)*3,nrow(dat),by=width*3)-1
    } else { dat <- readRDS(output_path(paste0(name,"_dat.rds"))); dat_midpoints <- seq(1:nrow(dat)) }
    dat <- dat[to_keep,]; set_global(dat)
    dat_midpoints <- dat_midpoints[to_keep]; set_global(dat_midpoints)
}

#################################################################################################################
# Order RBPs by eCLIP binding site clustering
#################################################################################################################
order_rbps <- function(padding=50, cell_lines=c("K562","HepG2","all")) {
    setup_features_env()
    rbps <- get_features_by_group("RBP")
    all_rbp_features <- get_features_by_group(paste0("RBP_",padding,"bp"))
    rbp_granges <- lapply(all_rbp_features, function(rbp_feature) { print(rbp_feature); rbp_g <- load_annotation(rbp_feature); names(rbp_g) <- rep(gsub("_[0-9]*bp$","",rbp_feature),length(rbp_g)); return(rbp_g) })
    rbp_granges <- Reduce(c, rbp_granges)
    
    rbp_counts <- table(names(rbp_granges))
    rbp_similarities <- lapply(rbps, function(rbp) {
        print(rbp)
        #rbp_pad_name = paste0(rbp,"_",padding,"bp")
        rbp_g <- rbp_granges[names(rbp_granges) == rbp]
        hits <- table(names(rbp_granges)[unique(subjectHits(findOverlaps(rbp_g, rbp_granges)))])
        return(hits / sapply(rbp_counts, function(x) min(x,rbp_counts[rbp])))
    })
    rbp_similarities <- Reduce(rbind, rbp_similarities)
    rownames(rbp_similarities) <- colnames(rbp_similarities)
    colnames(rbp_similarities) <- NULL
    rbp_similarities[1:6,1:6]
    rbp_similarities_all <- rbp_similarities
    
    for(cell_line in cell_lines) {
        cell_line_indices <- grepl(cell_line, rownames(rbp_similarities_all))
        if(sum(cell_line_indices) == 0) { cell_line_indices <- rep(TRUE, nrow(rbp_similarities_all)) }
        rbp_similarities <- rbp_similarities_all[cell_line_indices,cell_line_indices]
        
        filename = output_path(paste0("rbp_binding_hclust_",cell_line,".pdf"))
        pdf(filename)
        a <- hclust(as.dist(1 - rbp_similarities))
        hlcust_cex = 0.3
        if(cell_line == "all") { hlcust_cex = 0.25 }
        plot(a, cex=hlcust_cex)
        dev.off()
        pdf_to_png(filename)
        
        rownames(rbp_similarities)[a$order]
        
        filename = output_path(paste0("rbp_binding_similarities_",cell_line,".pdf"))
        pdf(filename)
        print(Heatmap(sqrt(sqrt(rbp_similarities[a$order,a$order])), 
                      show_heatmap_legend = TRUE, name = "weight", #title of legend
                      row_title = "", column_title = "RBP Binding Similarity Matrix",
                      cluster_rows=FALSE, cluster_columns=FALSE,
                      row_dend_side="left", column_dend_side="bottom",
                      row_names_side="left", column_names_side="top",#, #column_names_rot=0,
                      row_names_gp = gpar(col="black", fontsize = 3.5), column_names_gp = gpar(fontsize = 16, fontface="bold") # Text size for row and column names
        ))
        dev.off()
        pdf_to_png(filename)
        
        if(cell_line == cell_lines[length(cell_lines)]) {
            return(a$order)
        }
    }
}
#if(file.exists(output_path("rbp_order.rds"))) { rbp_order <- readRDS(output_path("rbp_order.rds"))
#} else { rbp_order <- order_rbps(); saveRDS(rbp_order, output_path("rbp_order.rds")) }

#################################################################################################################
# Fit the supermodel to training data, returning testing indices and a vector of the validation losses each epoch.
# The number of epochs, batch size, RBP order, and fraction data for training can all be customized.
#################################################################################################################
train_supermodel <- function(num_epochs=150, batch_size=128, rbp_order=1:160, fraction_for_training=0.8) {
    allowed_indices = which(rowSums(expected)>0) # & rowSums(af_tensor)>0)
    randomized_order_training_indices = sample(allowed_indices, floor(fraction_for_training*length(allowed_indices)))
    test_indices = allowed_indices[!(allowed_indices %in% randomized_order_training_indices)]
    #history <- model %>% fit(list(dat_gradcams[randomized_order_training_indices,,rbp_order,,drop=FALSE], dat_obs_exp[randomized_order_training_indices,,drop=FALSE], expected[randomized_order_training_indices,,drop=FALSE], af_tensor[randomized_order_training_indices,,drop=FALSE]), list(af_tensor[randomized_order_training_indices,,drop=FALSE]), epochs=num_epochs, batch_size=batch_size, validation_split=0.2)
    #history <- model %>% fit(list(dat_gradcams[randomized_order_training_indices,,rbp_order,,,drop=FALSE], dat_obs_exp[randomized_order_training_indices,,drop=FALSE], expected[randomized_order_training_indices,,drop=FALSE], af_tensor[randomized_order_training_indices,,drop=FALSE]), list(af_tensor[randomized_order_training_indices,,drop=FALSE]), epochs=num_epochs, batch_size=batch_size, validation_split=0.2)
    history <- model %>% fit(list(altref_gradcams[randomized_order_training_indices,,rbp_order,,drop=FALSE], dat_gradcams[randomized_order_training_indices,,rbp_order,,drop=FALSE], dat_obs_exp[randomized_order_training_indices,,drop=FALSE], expected[randomized_order_training_indices,,drop=FALSE], af_tensor[randomized_order_training_indices,,drop=FALSE]), list(af_tensor[randomized_order_training_indices,,drop=FALSE]), epochs=num_epochs, batch_size=batch_size, validation_split=0.2)
    print(history$metrics$loss)
    print(history$metrics$val_loss)
    
    filename=output_path("training_loss_vs_epochs.pdf")
    pdf(file=filename)
    plot(1:length(history$metrics$val_loss), history$metrics$val_loss, type="l", col="blue", lwd=2, xlab="epochs", ylab="loss", main="Training loss vs. epochs", cex.lab=1.4, cex.axis=1.4, cex.main=1.3)
    #lines(1:length(history$metrics$val_loss), history$metrics$loss, type="l", col="red", lwd=2)
    dev.off()
    pdf_to_png(filename)
    
    return_env <- new.env()
    return_env[["test_indices"]] <- test_indices
    return_env[["val_loss"]] <- history$metrics$val_loss
    return(return_env)
}

#################################################################################################################
# Plot model outputs at several informative levels in the neural network.
#################################################################################################################
plot_model_outputs <- function(model, test_indices, dat_eclip_overlaps=NULL, output_folder=output_path("model_outputs")) {
    dir.create(output_folder, showWarnings = FALSE)
    activation_model <- get_activation_model(model, c("s","rbp_disruption_output"))
    activations <- activation_model %>% predict(list(dat_gradcams[test_indices,,rbp_order,,drop=FALSE], dat_obs_exp[test_indices,,drop=FALSE], expected[test_indices,,drop=FALSE]*mu_scaling_factor, af_tensor[test_indices,,drop=FALSE]))
    names(activations) <- activation_model$output_names
    labels_all <- af_tensor[test_indices,,drop=FALSE]
    background_all <- expected[test_indices,,drop=FALSE]
    s_vals_all <- activations[["s"]]; s_vals_all[is.nan(s_vals_all)] <- 0; s_vals_all <- exp(s_vals_all * s_scaling_factor)
    rbp_disruption_all <- activations[["rbp_disruption_output"]]
    
    rbp_disruption_range <- range(rbp_disruption_all)
    rbp_disruption_all <- (rbp_disruption_all - rbp_disruption_range[1])/diff(rbp_disruption_range)
    af_range <- c(0,7.5e-4)
    #which(test_indices %in% dat_eclip_overlaps$query); test_indices[which(test_indices %in% dat_eclip_overlaps$query)]
    maps_result <- rbindlist(lapply(1:length(test_indices), function(i) {
        print(i)
        test_sequence_index = test_indices[i]
        single_input = length(test_sequence_index) == 1
        labels <- labels_all[i,,drop=FALSE]
        background <- background_all[i,,drop=FALSE] #pli <- dat_obs_exp[i]
        rbp_disruption <- rbp_disruption_all[i,]
        num_simulated = sample_size
        preds <- rpois(length(background), lambda=num_simulated*Ne*background)/num_simulated
        zero_preds <- preds == 0
        preds[zero_preds] <- rpois(sum(preds==0), lambda=num_simulated*Ne*background[zero_preds])/num_simulated
        preds_simulated <- preds; preds[zero_preds] <- NA #rep(NA, sum(zero_preds))
        s_vals <- s_vals_all[i,]
        s_vals_zero <- s_vals; s_vals_zero[s_vals_zero == 0] <- runif(sum(s_vals_zero == 0))*0.01
        s_vals_nonzero <- s_vals; s_vals_nonzero[s_vals==0] <- NA
        
        # Draw some important input and output metrics/distributions to file.
        filename = full_path(output_folder, paste0("model_outputs_seq",test_sequence_index,".pdf"))
        if(!is.null(dat_eclip_overlaps) && test_sequence_index %in% test_indices[which(test_indices %in% dat_eclip_overlaps$query)]) {
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
        plot(1:length(labels), labels, type="l", col="blue", xaxt="n", xlab="", ylab="AF (obs)", main="Supermodel Activation Analysis", ylim=af_range, cex.axis=cex_lab, cex.lab=cex_lab, cex.main=1.4)
        #mtext(paste0("Sequence ",test_sequence_index,", pLI = ",formatC(pli,format="f",digits=2)), cex=1.1)
        mtext(paste0("Sequence ",test_sequence_index), cex=1.1)
        par(mar=custom_par_mars[[2]])
        plot(1:length(labels), background, type="l", col="green", xaxt="n", xlab="", ylab="mu", ylim=c(range(background_all)), cex.axis=cex_lab, cex.lab=cex_lab)
        par(mar=custom_par_mars[[3]])
        plot(1:length(labels), rbp_disruption, type="l", col="black", xaxt="n", xlab="", ylab="d", yaxs="i", ylim=c((!is.null(dat_eclip_overlaps))*(-0.08),1.05), cex.axis=cex_lab, cex.lab=cex_lab)
        abline(h=0, col="black", lty=1)
        if(!is.null(dat_eclip_overlaps)) {
            eclip_annotation_loc = -0.04
            abline(h=eclip_annotation_loc, col="gray20", lty=3)
            eclip_overlaps <- dat_eclip_overlaps[dat_eclip_overlaps$query == test_sequence_index,]
            has_nearby_eclip = nrow(eclip_overlaps)>0
            if(has_nearby_eclip) {
                overlapped_rbps <- unique(eclip_overlaps$RBP)
                cols <- rainbow(length(overlapped_rbps)); names(cols) <- overlapped_rbps
                #col_alpha = 0.4; cols <- adjustcolor(cols, alpha.f=col_alpha)
                segments(x0=eclip_overlaps$start, y0=eclip_annotation_loc, x1=eclip_overlaps$end, y1=eclip_annotation_loc, col=cols[paste0(eclip_overlaps$RBP)], lty=1, lwd=3)
                legend("topleft", legend=overlapped_rbps, col=cols, pch=15)
            }
        } else { has_nearby_eclip = FALSE }
        par(mar=custom_par_mars[[4]])
        plot(1:length(labels), s_vals_zero, type="l", col="grey", xaxt="n", xlab="", ylab="s", ylim=c(range(s_vals_all)), cex.axis=cex_lab, cex.lab=cex_lab) #ylim=c(0,1/(mu_scaling_factor*s_lower_bound))
        lines(1:length(labels), s_vals_nonzero, type="l", col="red")
        par(mar=custom_par_mars[[5]])
        plot(1:length(labels), preds_simulated, type="l", col="grey", xlab="Position", ylab="AF (sim)", ylim=c(ylim=af_range), cex.axis=cex_lab, cex.lab=cex_lab)
        lines(1:length(labels), preds, type="l", col="orange")
        par(mfrow=c(1,1), mar=par_mar)
        dev.off()
        pdf_to_png(filename, output_folder="model_outputs")
        
        maps <- round(labels * sample_size)
        maps_res <- data.frame(t(rbind(maps, background, s_vals, rep(has_nearby_eclip,151))))
        colnames(maps_res) <- c("AC", "mutability", "s", "nearby_eclip")
        return(maps_res)
    }))
    return(maps_result)
}

#################################################################################################################
# Draw visual representations of SUPRNOVA filters. 
#################################################################################################################
visualize_filters <- function(model) {
    # Directly look at filters in the trained model to determine learned motifs.
    filters <- aperm(get_weights(model$get_layer("conv1"))[[1]][,,,], c(2,1,3))
    num_filters = dim(filters)[3]
    #rownames(filters) <- paste0(c(1:nrow(filters))) #colnames(filters) <- c("REF", "ALT")
    #filters
    #filter_cols <- c("green","orange")[apply(filters, 1, function(x) x[2] < x[1]) + 1]
    for(i in 1:30) { #1:nrow(filters)
        filt <- as.matrix(filters[,,i] %>% layer_activation_leaky_relu(alpha=0.3))
        filename = output_path(paste0("suprnova_filter",i,".pdf"))
        pdf(filename)
        print(Heatmap(filt, 
                      show_heatmap_legend = TRUE, name = "weight", #title of legend
                      row_title = "RBP", column_title = "Position",
                      cluster_rows=FALSE, cluster_columns=FALSE,
                      row_dend_side="left", column_dend_side="bottom",
                      row_names_side="left", column_names_side="top"#, #column_names_rot=0,
                      #row_names_gp = gpar(col=filter_cols, fontsize = 10), column_names_gp = gpar(fontsize = 20) # Text size for row and column names
        ))
        dev.off()
        pdf_to_png(filename)
    }
    
    # Directly look at weights in the trained model
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
}

#################################################################################################################
# Draw some important input and output metrics/distributions to file.
#################################################################################################################
analyze_supermodel_data <- function() {
    filename = output_path("gnomAD_distributions.pdf")
    pdf(file=filename, width=11)
    par(mfrow=c(2,3))
    plot(density(log10(expected)), type="l", lwd=1, col="blue", main="Background mutation rates", xlab="log10(mu)", cex.lab=1.2)
    mtext(paste0(formatC(100*sum(expected==0)/prod(dim(expected)),format="f",digits=2),"% of sites have 0 mu"))
    plot(density(log10(af_tensor)), type="l", lwd=1, col="blue", main="gnomAD allele frequencies", xlab=paste0("log10(AF), ",sample_size," WGS samples"), cex.lab=1.2)
    mtext(paste0(formatC(100*sum(af_tensor==0)/prod(dim(af_tensor)),format="f",digits=2),"% of sites have 0 AF"))
    plot(density(log10(af_tensor[rowSums(af_tensor)>0,])), type="l", lwd=1, col="blue", main="gnomAD allele frequencies", xlab=paste0("log10(AF), ",sample_size," WGS samples"), cex.lab=1.2)
    mtext(paste0(formatC(100*sum(af_tensor[rowSums(af_tensor)>0,]==0)/prod(dim(af_tensor[rowSums(af_tensor)>0,])),format="f",digits=2),"% of >0 region sites have 0 AF"))
    plot(density(dat_obs_exp), type="l", lwd=1, col="blue", main="Nearest gene pLIs", xlab="pLI (haploinsufficiency)", cex.lab=1.2)
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
}

