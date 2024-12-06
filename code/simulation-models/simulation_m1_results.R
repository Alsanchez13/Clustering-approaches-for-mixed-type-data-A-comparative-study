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
# Import custom functions
source("C:\\Users\\a809716\\Desktop\\these_airbus\\article_1\\retour_journal_probability_statistics\\code\\sample_mixed_gaussian.R")
source("C:\\Users\\a809716\\Desktop\\these_airbus\\article_1\\codes\\R\\clustering_bayesian_network.R")

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
#-----------
# CLUSTERING
#-----------
# Prepare the output table
lab.method <- rep(c("KAMILA", "K Prototyes", "PDQ", "Convex KM", "MBN", "LCM"),18)
lab.overlap <- c( rep(c("60%"), times = 108))
lab.cluster <- rep(rep(c("2","5","10"),each=36),1)
lab.contprop <- rep(rep(c("0.33","0.5","0.66"),each=6),6)
lab.size <- rep(rep(c("700", "1400"),each=18),3)
tab<- data.frame(Method= lab.method, Over.=lab.overlap, Clust. = lab.cluster, Size = lab.size, Cont.Prop = lab.contprop)
View(tab)
####################################
# CLUSTERING VARYING CLUSTER NUMBERS
####################################
# Parameters
var.num <- 12
cont.ratio <- c(0.33, 0.5, 0.66)
sizes <- c(700, 1400)
clusters <- c(2, 5, 10)

ari_mean <- c()
ari_sd <- c()
# Start timer
start <- Sys.time()

m1_gaussian_results <-
  foreach(cl = clusters, .combine = "c") %:%
  foreach(size = sizes, .combine = "c") %do% {
    print(cl)
    print(size)
    
    ari_results <-
      foreach(i = icount(length(cont.ratio)), .combine = "cbind") %:%
      foreach(l = icount(10), .combine = "rbind") %dopar% {
        
        con_var <- round(var.num*cont.ratio[i], digits = 0)
        nom_var <- round(var.num*(1- cont.ratio[i]), digits = 0)
        
        #----------------
        # DATA GENERATION
        #----------------
        # replace 0.6 by 0.3 to get results from 30% overlap
        df <- mixt_data(cont_var = con_var, over_cont = rep(0.6, con_var),
                        nom_var = nom_var, over_nom = rep(0.6, nom_var),
                        clusters = cl, obser = size)
        
        # add cluster membership column
        mem.count <- floor(rep((size/cl), cl))
        mem.count[cl] <- size - sum(mem.count[-cl])
        membership <- rep(1:cl, times = mem.count)
        df <-cbind(rep(1:cl, times = mem.count), df$data) #membership + data
        
        #-------------------
        # Data preprocessing
        #-------------------
        # separate and rename continuous variables
        con_df <- df[, 2:(con_var + 1)]
        con_name <- c(sprintf("con%d",1:con_var))
        colnames(con_df) <- con_name
        
        # separate and rename categorical variables
        cat_df <- df[, (con_var + 2):(con_var + nom_var + 1)]
        cat_name <- c(sprintf("cat%d",1:nom_var))
        colnames(cat_df) <- cat_name
        mode(cat_df) <- "character" #make entries into strings while preserving the matrix
        
        # true label
        true_label <- df[,1]
        
        #---------
        # ANALYSIS
        #---------
        # KAMILA
        kam_df <- kamila::kamila(conVar = as.data.frame(con_df),
                                 catFactor = as.data.frame(cat_df, stringsAsFactors = TRUE),
                                 numClust = cl,
                                 numInit = 20,
                                 maxIter = 25)
        
        ari_kam <- 0.1*(MixGHD::ARI(true_label, kam_df$finalMemb))
        
        # K PROTOTYPES
        kprot_df <- clustMixType::kproto(data.frame(con_df, as.data.frame(cat_df, stringsAsFactors = TRUE)),
                                         k = cl,
                                         lambda = clustMixType::lambdaest(
                                           data.frame(con_df,
                                                      as.data.frame(cat_df, stringsAsFactors = TRUE))),
                                         nstart = 10)
        
        ari_kprot <- 0.1*(MixGHD::ARI(true_label, kprot_df$cluster))
        
        # PDQ
        pdq_df <- FPDclustering::PDQ(x=df[,-1],
                                     k=cl,
                                     ini='kmd',
                                     dist='gower',
                                     cont=1:con_var,
                                     cat=(con_var + 1):(con_var + nom_var))
        
        
        ari_pdq <- 0.1*(MixGHD::ARI(true_label, pdq_df$label))
        
        # CONVEX K-MEANS
        convk_df <- kamila::gmsClust(conData = as.data.frame(con_df),
                                     catData = kamila::dummyCodeFactorDf(data.frame(apply(cat_df, 2, factor), stringsAsFactors = TRUE)),
                                     nclust = cl,
                                     searchDensity = 20)
        
        ari_convk <- 0.1*(MixGHD::ARI(true_label, convk_df$results$cluster))
        
        # Mixture BN
        bn_df <- clust_bn(df=data.frame(con_df, cat_df, stringsAsFactors = TRUE),
                          K=cl,
                          score_bn = 'bic-cg',
                          multi_start = 10)
        
        ari_bn<- 0.1*(MixGHD::ARI(true_label, bn_df$class))
        
        # LCM
        lcm <- VarSelLCM::VarSelCluster(data.frame(con_df, cat_df, stringsAsFactors=TRUE),
                                        cl,
                                        vbleSelec = FALSE,
                                        crit.varsel = "BIC")
        lcm_label <- VarSelLCM::fitted(lcm, type = "partition")
        
        ari_lcm <- 0.1*(MixGHD::ARI(true_label, lcm_label))
        
        return(c(ari_kam, ari_kprot, ari_pdq, ari_convk, ari_bn, ari_lcm))
      }
    
    ari_mean <- c(
      sum(ari_results[,1]), sum(ari_results[,2]), sum(ari_results[,3]),
      sum(ari_results[,4]), sum(ari_results[,5]), sum(ari_results[,6]),
      sum(ari_results[,7]), sum(ari_results[,8]), sum(ari_results[,9]),
      sum(ari_results[,10]), sum(ari_results[,11]), sum(ari_results[,12]),
      sum(ari_results[,13]), sum(ari_results[,14]), sum(ari_results[,15]),
      sum(ari_results[,16]), sum(ari_results[,17]), sum(ari_results[,18])
      )
    
    print(ari_mean)
    
    ari_sd <- c(
      sd(ari_results[,1])*10, sd(ari_results[,2])*10, sd(ari_results[,3])*10,
      sd(ari_results[,4])*10, sd(ari_results[,5])*10, sd(ari_results[,6])*10,
      sd(ari_results[,7])*10, sd(ari_results[,8])*10, sd(ari_results[,9])*10,
      sd(ari_results[,10])*10, sd(ari_results[,11])*10, sd(ari_results[,12])*10,
      sd(ari_results[,13])*10, sd(ari_results[,14])*10, sd(ari_results[,15])*10,
      sd(ari_results[,16])*10, sd(ari_results[,17])*10, sd(ari_results[,18])*10
      )
    
    return(list(ari_mean, ari_sd))
  }

# ARI for all the methods
tab$ari_mean <- unlist(m1_gaussian_results[c(1,3,5,7,9,11)])
tab$ari_sd <- unlist(m1_gaussian_results[c(2,4,6,8,10,12)])

# End timer
print( Sys.time() - start )

# STOP cluster
parallel::stopCluster(cl = my.cluster)
