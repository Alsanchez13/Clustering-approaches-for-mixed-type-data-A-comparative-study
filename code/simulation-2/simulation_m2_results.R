# TITLE: clustering AH mixed-type data

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
  "mclust",
  "MASS",
  "dplyr",
  'aricode'
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
source("C:\\Users\\a809716\\Desktop\\these_airbus\\article_1\\retour_journal_probability_statistics\\code\\sample_mix_exponential_discrete.R")
source('C:\\Users\\a809716\\Desktop\\these_airbus\\article_1\\codes\\R\\clustering_bayesian_network.R')

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
var <- 12 # 6
cont.ratio <- c(0.33, 0.5, 0.66)
cl <- 2
size <- 300 # 1200
balances <- c(0.5, 0.2)

# Prepare the output table
lab.method <- rep(c("KAMILA", "K Prototyes", "PDQ", "Convex KM", "MBN", "LCM"),6)
lab.var <- rep(rep(c("12"),each=6),6)
lab.contprop <- rep(rep(c("0.33", "0.5", "0.66"),each=6),2)
lab.size <- rep(rep(c("300"),each=6),6)
lab.balance <- rep(rep(c("1/2-1/2", "1/5-4/5"),each=18),1)
K <- rep(cl, 36)
tab<- data.frame(Method=lab.method, Cluster=K, Size=lab.size, Var.num = lab.var, Balance=lab.balance, Cont.Prop = lab.contprop)

View(tab)
# initialize metrics
ari_mean_all <- vector("list", length = length(balances))
ari_sd_all <- vector("list", length = length(balances))
ami_mean_all <- vector("list", length = length(balances))
ami_sd_all <- vector("list", length = length(balances))
#-----------
# START LOOP
#-----------
# Start timer
start <- Sys.time()

ah_clust_balance_results <-
  # foreach(size = sizes, .combine = "c") %:%
  foreach(bal = icount(length(balances)), .combine = "c") %do% {
    print(paste('size:', size))
    print(paste('balance:', balances[bal]))
    
    ari_mean <- vector("list", length = length(cont.ratio))
    ari_sd <- vector("list", length = length(cont.ratio))
    ami_mean <- vector("list", length = length(cont.ratio))
    ami_sd <- vector("list", length = length(cont.ratio))
    
    ari_results <-
      foreach(i = icount(length(cont.ratio)), .combine = "cbind") %do% {
        print(paste('cont.ratio:', cont.ratio[i]))
        
        con_var <- round(var*cont.ratio[i], digits = 0)
        nom_var <- round((var*(1- cont.ratio[i])/2), digits = 0)
        bin_var <- round(var - con_var - nom_var)
        
        ari_cont_prop <-
          foreach(l = icount(10), .combine = "rbind") %dopar% {
            #----------------
            # DATA GENERATION
            #----------------
            df <- generate_mix_exp_disc_k2(
              size=size, con_var=con_var, cat_var=nom_var,
              bin_var=bin_var, mix_weigts = balances[bal]
              )
            #-------------------
            # Data preprocessing
            #-------------------
            # separate and rename continuous variables
            con_df <- df[, 1:(con_var)]
            # separate and rename categorical variables
            cat_df <- df[, (con_var + 1):(con_var + nom_var + bin_var)]
            # true label
            true_label <- df[,ncol(df)]
            
            #---------
            # ANALYSIS
            #---------
            # KAMILA
            kam_df <- kamila::kamila(conVar = as.data.frame(con_df),
                                     catFactor = cat_df,
                                     numClust = cl,
                                     numInit = 20,
                                     maxIter = 25)
            
            ari_kam <- 0.1*(aricode::ARI(true_label, as.factor(kam_df$finalMemb)))
            ami_kam <- 0.1*(aricode::AMI(true_label, as.factor(kam_df$finalMemb)))
            
            # K PROTOTYPES
            kprot_df <- clustMixType::kproto(data.frame(con_df, cat_df),
                                             k = cl,
                                             lambda = clustMixType::lambdaest(
                                               data.frame(con_df,cat_df)),
                                             nstart = 10)
            
            ari_kprot <- 0.1*(aricode::ARI(true_label, as.factor(kprot_df$cluster)))
            ami_kprot <- 0.1*(aricode::AMI(true_label, as.factor(kprot_df$cluster)))
            
            # PDQ
            pdq_df <- FPDclustering::PDQ(x=data.frame(con_df, cat_df),
                                         k=cl,
                                         ini='kmd',
                                         dist='gower',
                                         cont=1:con_var,
                                         cat=(con_var + 1):(con_var + nom_var),
                                         bin = (con_var + nom_var + 1):(con_var + nom_var + bin_var))
            
            ari_pdq <- 0.1*(aricode::ARI(true_label, as.factor(pdq_df$label)))
            ami_pdq <- 0.1*(aricode::AMI(true_label, as.factor(pdq_df$label)))
            
            # CONVEX K-MEANS
            convk_df <- kamila::gmsClust(conData = as.data.frame(con_df),
                                         catData = kamila::dummyCodeFactorDf(cat_df),
                                         nclust = cl,
                                         searchDensity = 20)

            ari_convk <- 0.1*(aricode::ARI(true_label, as.factor(convk_df$results$cluster)))
            ami_convk <- 0.1*(aricode::AMI(true_label, as.factor(convk_df$results$cluster)))
            
            # MBN
            bn_df <- clust_bn(df=data.frame(con_df, cat_df, stringsAsFactors = TRUE),
                              K=cl,
                              score_bn='bic-cg',
                              multi_start = 10)
            
            ari_bn <- 0.1*(aricode::ARI(true_label, as.factor(bn_df$class)))
            ami_bn <- 0.1*(aricode::AMI(true_label, as.factor(bn_df$class)))
            
           # LCM
            lcm <- VarSelLCM::VarSelCluster(data.frame(con_df, cat_df),
                                            cl,
                                            vbleSelec = FALSE,
                                            crit.varsel = "BIC")
            lcm_label <- VarSelLCM::fitted(lcm, type = "partition")

            ari_lcm <- 0.1*(aricode::ARI(true_label, as.factor(lcm_label)))
            ami_lcm <- 0.1*(aricode::AMI(true_label, as.factor(lcm_label)))
            
            return(
              c(ari_kam, ari_kprot, ari_pdq, ari_convk, ari_bn, ari_lcm,
                ami_kam, ami_kprot, ami_pdq, ami_convk, ami_bn, ami_lcm)
              )
          }
        
        ari_mean[[i]] <- c(
          sum(ari_cont_prop[,1]), sum(ari_cont_prop[,2]), sum(ari_cont_prop[,3]),
          sum(ari_cont_prop[,4]), sum(ari_cont_prop[,5]), sum(ari_cont_prop[,6])
          )
        
        ari_sd[[i]] <-c(
          sd(ari_cont_prop[,1])*10, sd(ari_cont_prop[,2])*10, sd(ari_cont_prop[,3])*10,
          sd(ari_cont_prop[,4])*10, sd(ari_cont_prop[,5])*10, sd(ari_cont_prop[,6])*10
          )
        
        ami_mean[[i]] <- c(
          sum(ari_cont_prop[,7]), sum(ari_cont_prop[,8]), sum(ari_cont_prop[,9]),
          sum(ari_cont_prop[,10]), sum(ari_cont_prop[,11]), sum(ari_cont_prop[,12])
          )
        
        ami_sd[[i]] <-c(
          sd(ari_cont_prop[,7])*10, sd(ari_cont_prop[,8])*10, sd(ari_cont_prop[,9])*10,
          sd(ari_cont_prop[,10])*10, sd(ari_cont_prop[,11])*10, sd(ari_cont_prop[,12])*10
        )
      }
    
    ari_mean_all[[bal]] <- unlist(ari_mean)
    ari_sd_all[[bal]] <- unlist(ari_sd)
    ami_mean_all[[bal]] <- unlist(ami_mean)
    ami_sd_all[[bal]] <- unlist(ami_sd)
    
    return(
      list(mean_ari = unlist(ari_mean_all),
           sd_ari = unlist(ari_sd_all),
           mean_ami = unlist(ami_mean_all),
           sd_ami = unlist(ami_sd_all))
    )
  }
# ARI for all the methods
tab$ari_mean <- unlist(unlist(ah_clust_balance_results[[5]]))
tab$ari_sd <- unlist(unlist(ah_clust_balance_results[[6]]))
tab$ami_mean <- unlist(unlist(ah_clust_balance_results[[7]]))
tab$ami_sd <- unlist(unlist(ah_clust_balance_results[[8]]))

# End timer
print( Sys.time() - start )
# STOP cluster
parallel::stopCluster(cl = my.cluster)
#-----------------------------
# VARYING SIZES FOR 3 CLUSTERS
#-----------------------------
# Parameters
var <- 12 # 6
cont.ratio <- c(0.33, 0.5, 0.66)
cl <- 3
size <- 300 # 1200
balances <- list(c(1/3,1/3,1/3), c(1/2, 1/4, 1/4))

# Prepare the output table
lab.method <- rep(c("KAMILA", "K Prototyes", "PDQ", "Convex K Means", "MBN", "LCM"),6)
lab.var <- rep(rep(c("12"),each=6),6)
lab.contprop <- rep(rep(c("0.33", "0.5", "0.66"),each=6),2)
lab.size <- rep(rep(c("300"),each=6),6)
lab.balance <- rep(rep(c("1/2-1/2", "1/5-4/5"),each=18),1)
lab.balance <- rep(rep(c("1/3-1/3-1/3", "1/2-1/4-1/4"),each=18),1)
K <- rep(cl, 36)
tab<- data.frame(Method=lab.method, Cluster=K, Size=lab.size, Var.num = lab.var, Balance=lab.balance, Cont.Prop = lab.contprop)

View(tab)
# initialize metrics
ari_mean_all <- vector("list", length = length(balances))
ari_sd_all <- vector("list", length = length(balances))
ami_mean_all <- vector("list", length = length(balances))
ami_sd_all <- vector("list", length = length(balances))
#-----------
# START LOOP
#-----------
# Start timer
start <- Sys.time()

ah_3_clust_results <-
  # foreach(size = sizes, .combine = "c") %:%
  foreach(bal = icount(length(balances)), .combine = "c") %do% {
    print(paste('size:', size))
    print(paste('balance:', balances[bal]))
    
    ari_mean <- vector("list", length = length(cont.ratio))
    ari_sd <- vector("list", length = length(cont.ratio))
    ami_mean <- vector("list", length = length(cont.ratio))
    ami_sd <- vector("list", length = length(cont.ratio))
    
    ari_results <-
      foreach(i = icount(length(cont.ratio)), .combine = "cbind") %do% {
        print(paste('cont.ratio:', cont.ratio[i]))
        
        con_var <- round(var*cont.ratio[i], digits = 0)
        nom_var <- round((var*(1- cont.ratio[i])/2), digits = 0)
        bin_var <- round(var - con_var - nom_var)
        
        ari_cont_prop <-
          foreach(l = icount(10), .combine = "rbind") %dopar% {
            
            #----------------
            # DATA GENERATION
            #----------------
            df <- generate_mix_exp_disc_k3(
              size=size, con_var=con_var, cat_var=nom_var,
              bin_var=bin_var, mix_weights=balances[[bal]]
              )
            #-------------------
            # Data preprocessing
            #-------------------
            # separate and rename continuous variables
            con_df <- df[, 1:(con_var)]
            # separate and rename categorical variables
            cat_df <- df[, (con_var + 1):(con_var + nom_var + bin_var)]
            # true label
            true_label <- df[,ncol(df)]
            
            #---------
            # ANALYSIS
            #---------
            # KAMILA
            kam_df <- kamila::kamila(conVar = as.data.frame(con_df),
                                     catFactor = cat_df,
                                     numClust = cl,
                                     numInit = 20,
                                     maxIter = 25)
            
            ari_kam <- 0.1*(aricode::ARI(true_label, as.factor(kam_df$finalMemb)))
            ami_kam <- 0.1*(aricode::AMI(true_label, as.factor(kam_df$finalMemb)))
            
            # K PROTOTYPES
            kprot_df <- clustMixType::kproto(data.frame(con_df, cat_df),
                                             k = cl,
                                             lambda = clustMixType::lambdaest(
                                               data.frame(con_df,cat_df)),
                                             nstart = 10)
            
            ari_kprot <- 0.1*(aricode::ARI(true_label, as.factor(kprot_df$cluster)))
            ami_kprot <- 0.1*(aricode::AMI(true_label, as.factor(kprot_df$cluster)))
            
            # PDQ
            pdq_df <- FPDclustering::PDQ(x=data.frame(con_df, cat_df),
                                         k=cl,
                                         ini='kmd',
                                         dist='gower',
                                         cont=1:con_var,
                                         cat=(con_var + 1):(con_var + nom_var),
                                         bin = (con_var + nom_var + 1):(con_var + nom_var + bin_var))
            
            ari_pdq <- 0.1*(aricode::ARI(true_label, as.factor(pdq_df$label)))
            ami_pdq <- 0.1*(aricode::AMI(true_label, as.factor(pdq_df$label)))
            
            # CONVEX K-MEANS
            convk_df <- kamila::gmsClust(conData = as.data.frame(con_df),
                                         catData = kamila::dummyCodeFactorDf(cat_df),
                                         nclust = cl,
                                         searchDensity = 20)
            
            ari_convk <- 0.1*(aricode::ARI(true_label, as.factor(convk_df$results$cluster)))
            ami_convk <- 0.1*(aricode::AMI(true_label, as.factor(convk_df$results$cluster)))
            
            # MBN
            bn_df <- clust_bn(df=data.frame(con_df, cat_df, stringsAsFactors = TRUE),
                              K=cl,
                              score_bn='bic-cg',
                              multi_start = 10)
            
            ari_bn <- 0.1*(aricode::ARI(true_label, as.factor(bn_df$class)))
            ami_bn <- 0.1*(aricode::AMI(true_label, as.factor(bn_df$class)))
            
            # LCM
            lcm <- VarSelLCM::VarSelCluster(data.frame(con_df, cat_df),
                                            cl,
                                            vbleSelec = FALSE,
                                            crit.varsel = "BIC")
            lcm_label <- VarSelLCM::fitted(lcm, type = "partition")
            
            ari_lcm <- 0.1*(aricode::ARI(true_label, as.factor(lcm_label)))
            ami_lcm <- 0.1*(aricode::AMI(true_label, as.factor(lcm_label)))
            
            return(
              c(ari_kam, ari_kprot, ari_pdq, ari_convk, ari_bn, ari_lcm,
                ami_kam, ami_kprot, ami_pdq, ami_convk, ami_bn, ami_lcm)
            )
          }
        ari_mean[[i]] <- c(
          sum(ari_cont_prop[,1]), sum(ari_cont_prop[,2]), sum(ari_cont_prop[,3]),
          sum(ari_cont_prop[,4]), sum(ari_cont_prop[,5]), sum(ari_cont_prop[,6])
          )
        
        ari_sd[[i]] <-c(
          sd(ari_cont_prop[,1])*10, sd(ari_cont_prop[,2])*10, sd(ari_cont_prop[,3])*10,
          sd(ari_cont_prop[,4])*10, sd(ari_cont_prop[,5])*10, sd(ari_cont_prop[,6])*10
        )
        
        ami_mean[[i]] <- c(
          sum(ari_cont_prop[,7]), sum(ari_cont_prop[,8]), sum(ari_cont_prop[,9]),
          sum(ari_cont_prop[,10]), sum(ari_cont_prop[,11]), sum(ari_cont_prop[,12])
        )
        
        ami_sd[[i]] <-c(
          sd(ari_cont_prop[,7])*10, sd(ari_cont_prop[,8])*10, sd(ari_cont_prop[,9])*10,
          sd(ari_cont_prop[,10])*10, sd(ari_cont_prop[,11])*10, sd(ari_cont_prop[,12])*10
        )
      }
    
    ari_mean_all[[bal]] <- unlist(ari_mean)
    ari_sd_all[[bal]] <- unlist(ari_sd)
    ami_mean_all[[bal]] <- unlist(ami_mean)
    ami_sd_all[[bal]] <- unlist(ami_sd)
    
    return(
      list(mean_ari = unlist(ari_mean_all),
           sd_ari = unlist(ari_sd_all),
           mean_ami = unlist(ami_mean_all),
           sd_ami = unlist(ami_sd_all))
    )
  }

# ARI for all the methods
tab$ari_mean <- unlist(unlist(ah_3_clust_results[[5]]))
tab$ari_sd <- unlist(unlist(ah_3_clust_results[[6]]))
tab$ami_mean <- unlist(unlist(ah_3_clust_results[[7]]))
tab$ami_sd <- unlist(unlist(ah_3_clust_results[[8]]))

# End timer
print( Sys.time() - start )
# STOP cluster
parallel::stopCluster(cl = my.cluster)
