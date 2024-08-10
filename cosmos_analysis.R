library(decoupleR)
library(OmnipathR)
library(cosmosR)
library(org.Hs.eg.db)
library(CARNIVAL)
library(reshape2)
library(readr)

# Get the CollecTRI network
net <- get_collectri(organism = 'human', split_complexes = FALSE)

##### TRANSCRIPTOME
transcriptome <- read.csv("/workspaces/env/AFTER_PY/exp_input_editedit.csv")
t_matrix <- transcriptome$l2fc  
names(t_matrix) <- transcriptome$gene
t_matrix <- as.matrix(t_matrix)
colnames(t_matrix) <- "score"

##### METABOLOME
metabolic_data <- read.csv("/workspaces/env/AFTER_PY/metabolic_input_editedit.csv")
# Run decoupleR analysis
metabolic_acts <- run_ulm(mat = as.matrix(metabolic_data), 
                          net = net, 
                          .source = 'source', 
                          .target = 'target', 
                          .mor = 'mor', 
                          minsize = 5)

# Sort and view results
metabolic_acts_sorted <- metabolic_acts[order(metabolic_acts$p_value), ]
head(metabolic_acts_sorted, 10)



### PHOSPHOPROTEMICS
ph <- read.csv("/workspaces/env/ph_FILT.csv")
# Create the matrix
phospho_matrix <- -log10(ph$adj.p.value)
names(phospho_matrix) <- ph$gene
phospho_matrix <- as.matrix(phospho_matrix)
colnames(phospho_matrix) <- "score"


#### combined phosphopr. & transcr.
# num [1:29253, 1] -0.822 -0.822 -0.828 -0.841 -0.841 ...
# - attr(*, "dimnames")=List of 2
#  ..$ : chr [1:29253] "HDAC7" "LRPAP1" "CCNYL1" "IQSEC1" ...
#  ..$ : chr "score"
combined_matrix <- rbind(phospho_matrix, t_matrix)




####### RUN decoupler
tf_acts <- run_ulm(mat = combined_matrix, 
                         net = net, 
                         .source = 'source', 
                         .target = 'target', 
                         .mor = 'mor', 
                         minsize = 5
                    )


### visualize
tf_acts_sorted <- tf_acts[order(tf_acts$p_value), ]
head(tf_acts_sorted)
# A tibble: 6 × 5
# statistic source condition score  p_value
#  <chr>     <chr>  <chr>     <dbl>    <dbl>
#1 ulm       SP1    score     14.6  6.44e-48
#2 ulm       MYC    score     12.5  7.87e-36


#### SIGN Data
sig_data <- read.csv("/workspaces/env/AFTER_PY/sig_data_editedit.csv")

sig_data <- as.matrix(sig_data)
str(sig_data)
metabolomic_data <- read.csv("/workspaces/env/AFTER_PY/sig_data_editedit.csv")
metabolomic_data <- subset(metabolomic_data, select = -X)
head(metabolomic_data, n=5)
met_matrix <- as.matrix(metabolomic_data)
names(met_matrix) <- metabolic_data$HMDB_ID
str(met_matrix)


load_tf_regulon_dorothea(confidence = c("A", "B", "C"))

# Extract the 'score' column and convert to numeric
sig_data <- as.data.frame(sig_data)
scores <- as.numeric(trimws(sig_data$score))
names(scores) <- sig_data$signaling_nodes
is.vector(scores)

met_matrix <- as.data.frame(met_matrix)
str(met_matrix)
metab <- as.numeric(trimws(met_matrix$score))
names(metab) <- met_matrix$signaling_nodes
is.vector(metab)

tra <- as.data.frame(transcriptome)
str(tra)
tr_vector <- as.numeric(trimws(tra$l2fc))
names(tr_vector) <- tra$gene
is.vector(tr_vector)



net_meta <- as.data.frame(net)
str(net_meta)
# Rename a column by index
colnames(net_meta)[3] <- "interaction"  
head(net_meta)
net_meta <- net_meta[, c("source", "interaction", "target")]




test_for <- preprocess_COSMOS_signaling_to_metabolism(meta_network = net_meta,
                                        signaling_data = scores,
                                        metabolic_data = metab,
                                                      diff_expression_data = tr_vector,
                                                      maximum_network_depth = 15,
                                                      remove_unexpressed_nodes = TRUE,
                                                      CARNIVAL_options = CARNIVAL_options)




CARNIVAL_options <- default_CARNIVAL_options(solver = "lpSolve")
CARNIVAL_options$threads <- 2
CARNIVAL_options$mipGAP <- 0.05
CARNIVAL_options$timelimit <- 3600




"""
[1] "COSMOS: all 622 signaling nodes from data were found in the meta PKN"
[1] "COSMOS: all 622 metabolic nodes from data were found in the meta PKN"
[1] "COSMOS: 3710 of the 12086 genes in expression data were found as transcription factor target"
[1] "COSMOS: 3710 of the 5321 transcription factor targets were found in expression data"
[1] "COSMOS: removing unexpressed nodes from PKN..."
[1] "COSMOS: 14554 interactions removed"
[1] "COSMOS: 230 input/measured nodes are not in PKN any more: AHR, AIRE, AP1, ARX, ASCL1, ASCL2 and 224 more."
[1] "COSMOS: 230 input/measured nodes are not in PKN any more: AHR, AIRE, AP1, ARX, ASCL1, ASCL2 and 224 more."
[1] "COSMOS: removing nodes that are not reachable from inputs within 15 steps"
[1] "COSMOS: 239 from  16819 interactions are removed from the PKN"
[1] "COSMOS: removing nodes that are not observable by measurements within 15 steps"
[1] "COSMOS: 12984 from  16580 interactions are removed from the PKN"
[1] "COSMOS: 9 input/measured nodes are not in PKN any more: ARID4B, ELF2, FOXP2, IKZF4, KCNIP3, NFE2L3 and 3 more."
[1] "COSMOS: 9 input/measured nodes are not in PKN any more: ARID4B, ELF2, FOXP2, IKZF4, KCNIP3, NFE2L3 and 3 more."
[1] "COSMOS:  2374 interactions are removed from the PKN based on consistency check between TF activity and gene expression"
[1] "COSMOS: 43 input/measured nodes are not in PKN any more: BTG2, CEBPG, CREB3, CREB3L1, CUX1, ELK4 and 37 more."
[1] "COSMOS: 43 input/measured nodes are not in PKN any more: BTG2, CEBPG, CREB3, CREB3L1, CUX1, ELK4 and 37 more."
--- Start of the CARNIVAL pipeline ---
22:53:48 16.07.2024 Carnival flavour: vanilla
22:53:48 16.07.2024 Generating variables for lp problem
22:53:48 16.07.2024 Done: generating variables for lp problem
Saving preprocessed data.
Done: saving parsed data: /workspaces/env//parsedData_t22_53_48d16_07_2024n30.RData
22:53:48 16.07.2024 Generating formulation for LP problem
22:53:49 16.07.2024 Done: generating formulation for LP problem.
Saving LP file
Done: Saving LP file: /workspaces/env//lpFile_t22_53_48d16_07_2024n30.lp
22:53:50 16.07.2024 Solving LP problem
Writing cplex command file
Done: writing cplex command file
/workspaces/env/ampl.linux-intel64/cplex: bad option -f
usage: cplex [options] stub [-AMPL] [<assignment> ...]

Options:
        --  {end of options}
        -=  {show name= possibilities}
        -?  {show usage}
        -bf {read boundsfile f}
        -e  {suppress echoing of assignments}
        -of {write .sol file to file f}
        -s  {write .sol file (without -AMPL)}
        -v  {just show version}
Saving results...
Error in solversFunctions$solve(carnivalOptions) : 
  CPLEX solution file is not found. CPLEX was likely interrupted (exceeding memory limit is the usual cause).
           Try to increase the available resources (memory) or reducing the PKN. 
In addition: Warning messages:
1: In preprocessPriorKnowledgeNetwork(priorKnowledgeNetwork) :
  self loop(s) detected and removed from prior knowledge network.
2: In checkPerturbations(perturbations, nodesPriorKnowledgeNetwork) :
  These perturbation nodes are not in prior knowledge network and will be ignored: DEAF1 | HMGB2 | PLAGL2 | VEZF1
3: In checkMeasurements(measurements, nodesPriorKnowledgeNetwork) :
  These measurement nodes are not in prior knowledge network and will be ignored: DEAF1 | HMGB2 | PLAGL2 | VEZF1
"""


load('/workspaces/env/parsedData_t23_28_08d16_07_2024n34.RData')
test_result_for <- run_COSMOS_signaling_to_metabolism(data = test_for,
                                                      CARNIVAL_options = CARNIVAL_options)



formated_result_for <- format_COSMOS_res(test_result_for)
head(formated_result_for)
write.csv(formated_result_for, file = "formated_result_for.csv")
summary(formated_result_for)

#
#
#

####BACK
test_back <- preprocess_COSMOS_metabolism_to_signaling(meta_network = net_meta,
                                        signaling_data = scores,
                                        metabolic_data = metab,
                                                       diff_expression_data = tr_vector,
                                                       maximum_network_depth = 15,
                                                       remove_unexpressed_nodes = FALSE,
                                                       CARNIVAL_options = CARNIVAL_options)



#CARNIVAL_options$timelimit <- 28800
test_result_back <- run_COSMOS_metabolism_to_signaling(data = test_back,
                                                       CARNIVAL_options = CARNIVAL_options)



formated_result_back <- format_COSMOS_res(test_result_back)




full_sif <- as.data.frame(rbind(formated_result_for[[1]], formated_result_back[[1]]))
full_sif <- full_sif[full_sif$Weight>0,]
full_attributes <- as.data.frame(rbind(formated_result_for[[2]], formated_result_back[[2]]))

full_sif <- unique(full_sif)
full_attributes <- unique(full_attributes)




#
#
#
#
#
sig_input <- scores
metab_input <- metab
RNA_inpus <- tr_vector


#Choose which compartment to assign to the metabolic measurments
metab_input <- prepare_metab_inputs(metab_input, c("c","m"))

##Filter significant inputs
sig_input <- sig_input[abs(sig_input) > 2]
metab_input <- metab_input[abs(metab_input) > 2]



meta_network <- meta_network_cleanup(netM)

#Remove genes that are not expressed from the meta_network
meta_network <- cosmosR:::filter_pkn_expressed_genes(names(RNA_inpus), meta_pkn = meta_network)



#Filter inputs and prune the meta_network to only keep nodes that can be found downstream of the inputs
#The number of step is quite flexible, 7 steps already covers most of the network
n_steps <- 8


# in this step we prune the network to keep only the relevant part between upstream and downstream nodes
sig_input <- cosmosR:::filter_input_nodes_not_in_pkn(sig_input, meta_network)

#meta_network <- cosmosR:::keep_controllable_neighbours(meta_network, n_steps, names(sig_input))
metab_input <- cosmosR:::filter_input_nodes_not_in_pkn(metab_input, meta_network)
meta_network <- cosmosR:::keep_observable_neighbours(meta_network, n_steps, names(metab_input))
sig_input <- cosmosR:::filter_input_nodes_not_in_pkn(sig_input, meta_network)


#compress the network
meta_network_compressed_list <- compress_same_children(meta_network, sig_input = sig_input, metab_input = metab_input)
meta_network_compressed <- meta_network_compressed_list$compressed_network
node_signatures <- meta_network_compressed_list$node_signatures
duplicated_parents <- meta_network_compressed_list$duplicated_signatures
meta_network_compressed <- meta_network_cleanup(meta_network_compressed)

