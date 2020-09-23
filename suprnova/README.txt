FILE DESCRIPTIONS:

suprnova_data_setup.R: Responsible for calculating and saving all of the different input tensor objects that SUPRNOVA will use for training.

suprnova_design.R: Responsible for supermodel definition and training, as well as saving and loading models. Allows customizable generation of a model that can then be used for prediction and analysis. The specific keras/Tensorflow design is defined in the define_supermodel_architecture() function, with two versions available: 1. regional/rSUPRNOVA, which makes non-basepair specific predictions for s for every position in the input sequence region, and variant/vSUPRNOVA, which makes a prediction for each possible allele change/variant for the centre-most position in the input sequence region. 

suprnova_analysis.R: Contains all analysis code, including comparison of SUPRNOVA predicted s and d with other metrics such as MAPS, association of s with gene constraint, region type, and other sanity checks, and case vs. control enrichment boost in WGS variant data when s and/or d are included as features. This script calls both of the above R files to set up the data and trained model.

alex_suite.R: Contains many useful helper functions, used in all of the above R files.


WORKFLOW:

To make changes to the supermodel (experimentation with different network structure, layer types, regularization approaches, etc.), edit the define_supermodel_architecture() function in suprnova_design.R, load required tensors into the environment with load_supermodel_training_data(), and finally call train_supermodel() with the specified number of epochs, batch size, and fraction of data to use for training vs. testing. The R object "model" in the environment can then be used to make predictions given new data, as demonstrated in suprnova_analysis.R. The suprnova_design.R functions visualize_filters() and plot_model_ouputs() can then be used to visualize and understand the learned RBP convolution filters and intermediate network outputs, respectively. Finally, save_supermodel() and load_supermodel() help make the trained model persistent in file and share-able, while get_activation_model() returns a version of the supermodel that outputs at every queried layer stage for each input, rather than simply the final s prediction. To see which layers are available, run get_model_layer_names().


DEPENDENCIES: 

//TODO
