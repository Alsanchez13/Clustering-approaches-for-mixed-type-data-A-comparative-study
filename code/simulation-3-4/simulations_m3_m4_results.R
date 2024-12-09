# TITLE: clustering mixture TAN BN mixed-type data

#-----------
# LIBRARIES
#----------
list.of.packages <- c(
  "foreach",
  "doParallel",
  "ranger",
  "palmerpenguins",
  # "tidyverse",
  "kableExtra",
  "iterators",
  "itertools",
  # clustering algorithms
  "kamila",
  "clustMixType",
  "FPDclustering",
  "bnlearn",
  'VarSelLCM',
  # clustering metrics
  "MixGHD",
  "mclust",
  "MASS",
  "dplyr"
)

#loading packages
for(package.i in list.of.packages){
  suppressPackageStartupMessages(
    library(
      package.i, 
      character.only = TRUE
    )
  )
}
# Functions
source('...\\clustering_bayesian_network.R')
source("...\\sample_sbn_bn.R")
source("...\\sample_mixture_bn.R")
#---------------------------------
# SETTING UP DISTRIBUTED COMPUTING
#---------------------------------
# Number of cores to use for parallel computing
n.cores <- parallel::detectCores() - 1
#create the cluster
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)
#check cluster definition (optional)
print(my.cluster)
#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)
#check if it is registered (optional)
foreach::getDoParRegistered()
#how many workers are available? (optional)
foreach::getDoParWorkers()
#------------------------------
# CLUSTERING FOR K = 2 CLUSTERS
#------------------------------
# Parameters
cl <- 6 # 2 or 4
size <- c(300, 1200)
# Prepare the output table
lab.method <- rep(c("KAMILA", "K Prototyes", "PDQ", "Convex K Means", "MBN", "LCM"),2)
lab.size <-rep(c("300", "1200"),each=6)
G <- rep(rep(c('6'), each=6),1)
tab<- data.frame(Method= lab.method, Size = lab.size, Cluster= G)

View(tab)
# initialize metrics
ari_mean <- c()
ari_sd <- c()
#-----------
# START LOOP
#-----------
# Start timer
start <- Sys.time()
# Import parameters of the simulation model
# sim_model <- readRDS("...\\m3_parameters2.rds")
# sim_model <- readRDS("...\\m3_parameters4.rds")
sim_model <- readRDS("...\\m3_parameters6.rds")
# sim_model <- readRDS("...\\m4_parameters2.rds")
# sim_model <- readRDS("...\\m4_parameters4.rds")
# sim_model <- readRDS("...\\m4_parameters6.rds")

bn_simulations_results <-
  foreach(s = size, .combine = "c") %do% {
    print(s)
        
    ari_results <-
      foreach(l = icount(10), .combine = "rbind") %dopar% {
      
        #----------------
        # DATA GENERATION
        #----------------
        # M3 model
        df <- bnlearn::rbn(sim_model$model, s)
        # M4 model
        # df <- sample_data_sim_model_mixture_bn(sim_model, s)
        
        #-------------------
        # Data preprocessing
        #-------------------
        # separate and rename continuous variables
        con_df <- df[, c("X3", "X4", "X6")]
        con_name <- c(sprintf("con%d",1:3))
        colnames(con_df) <- con_name

        # separate and rename categorical variables
        cat_df <- df[, c("X1", "X2", "X5")]
        cat_name <- c(sprintf("cat%d",1:3))
        colnames(cat_df) <- cat_name
        
        ### True label
        # M3 model
        true_label <- df[,1]
        # M4 model
        # true_label <- df[,7]
        
        #---------
        # ANALYSIS
        #---------
        # KAMILA
        kam_df <- kamila::kamila(
          conVar = as.data.frame(con_df),
          catFactor = as.data.frame(cat_df, stringsAsFactors = TRUE),
          numClust = cl,
          numInit = 20,
          maxIter = 25
          )
        
        ari_kam <- 0.1*(MixGHD::ARI(true_label, kam_df$finalMemb))
        
        # K PROTOTYPES
        kprot_df <- clustMixType::kproto(
          data.frame(con_df, cat_df, stringsAsFactors = TRUE),
          k = cl,
          lambda = clustMixType::lambdaest(
            data.frame(con_df, cat_df, stringsAsFactors = TRUE)
            ),
          nstart = 10
          )

        ari_kprot <- 0.1*(MixGHD::ARI(true_label, kprot_df$cluster))
        
        # PDQ
        pdq_df <- FPDclustering::PDQ(
          x=data.frame(con_df, cat_df, stringsAsFactors = TRUE),
          k=cl,
          ini='random',
          dist='gower',
          cont=1:3,
          cat=4:5,
          bin = 6
          )

        ari_pdq <- 0.1*(MixGHD::ARI(true_label, pdq_df$label))
        
        # CONVEX K-MEANS
        convk_df <- kamila::gmsClust(
          conData = as.data.frame(con_df),
          catData = kamila::dummyCodeFactorDf(
            data.frame(apply(cat_df, 2, factor), stringsAsFactors = TRUE)
            ),
          nclust = cl,
          searchDensity = 20
          )
        
        ari_convk <- 0.1*(MixGHD::ARI(true_label, convk_df$results$cluster))

        # BN
        bn_df <- clust_bn(
          df=data.frame(con_df, cat_df, stringsAsFactors = TRUE),
          K=cl,
          score_bn='bic-cg',
          multi_start = 10
          )
        
        ari_bn<- 0.1*(MixGHD::ARI(true_label, bn_df$class))
        
        # LCM
        lcm <- VarSelLCM::VarSelCluster(
          data.frame(con_df, cat_df),
          cl,
          vbleSelec = FALSE,
          crit.varsel = "BIC"
          )
        lcm_label <- VarSelLCM::fitted(lcm, type = "partition")
        
        ari_lcm <- 0.1*(MixGHD::ARI(true_label, lcm_label))

        return(c(ari_kam, ari_kprot, ari_pdq, ari_convk, ari_bn, ari_lcm))
      }
    
    ari_mean <- c(sum(ari_results[,1]), sum(ari_results[,2]), sum(ari_results[,3]),
                  sum(ari_results[,4]), sum(ari_results[,5]), sum(ari_results[,6]))
    
    ari_sd <- c(sd(ari_results[,1])*10, sd(ari_results[,2])*10, sd(ari_results[,3])*10, 
               sd(ari_results[,4])*10, sd(ari_results[,5])*10, sd(ari_results[,6])*10)
    
    return(list(ari_mean, ari_sd))
  }
# ARI for all the methods
tab$ari_mean <- unlist(bn_simulations_results[c(1,3)])
tab$ari_sd <- unlist(bn_simulations_results[c(2,4)])

# End timer
print( Sys.time() - start)

# STOP cluster
parallel::stopCluster(cl = my.cluster)
