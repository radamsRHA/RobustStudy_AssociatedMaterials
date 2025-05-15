################
# load depends #
################
rm(list = ls());library(geiger); library(phybase); library(MASS); library(phangorn); library("caper"); library(nlme); library(piecewiseSEM); library(phyclust); library(distory); library(PRDATR); library(gaussDiff); library(sfsmisc); library(robust); library(robustbase); library(estimatr)

###############################
# set experimental parameters #
###############################
numeric.NumberOfSpecies <- 1000
vector.Lambda <- 1
numeric.NumberOfReps <- 10^5
numeric.Beta <- 0
numeric.TreeDepth <- 10
numeric.p <- 75

vector.Model_01_Theta <- c(1); names(vector.Model_01_Theta) <- c("Sig2")
matrix.RESULTS<- matrix(nrow = numeric.NumberOfReps*length(vector.Lambda), ncol = 36); colnames(matrix.RESULTS) <- c("NumberSpecies", "Lambda", "DeathRate", "TrueBeta", "NumberTraitsP",
                                                                                                                     "RF", "TreeDist", "PRDTR", "MeanBrlen_S", "MeanBrlen_G","MeanAge_S", "MeanAge_G",
                                                                                                                     "Pvalue_SS", "Pvalue_SG","Pvalue_GG", "Pvalue_GS", 
                                                                                                                     "BETA_SS", "BETA_SG","BETA_GG", "BETA_GS", 
                                                                                                                     "RSE_SS", "RSE_SG","RSE_GG", "RSE_GS", 
                                                                                                                     "LogLik_SS", "LogLik_SG","LogLik_GG", "LogLik_GS",
                                                                                                                     "MSE_SS", "MSE_SG","MSE_GG", "MSE_GS", 
                                                                                                                     "RR_P_SS","RR_P_SG","RR_P_GG","RR_P_GS")

#############################
# new trees every replicate #
#############################
numeric.Counter <- 0
for (i in 1:length(vector.Lambda)){
  
  #########################
  # get simulation params #
  #########################
  print(i)
  numeric.Lambda<- vector.Lambda[i]
  numeric.DeathRate <- numeric.Lambda/2

  ###########################
  # loop through replicates #
  ###########################
  for (j in 1:numeric.NumberOfReps){
    
    ################
    # SPECIES TREE #
    ################
    handle.SpeciesTree <- TreeSim::sim.bd.taxa.age(n = numeric.NumberOfSpecies, numbsim = 1, lambda = numeric.Lambda, mu = numeric.DeathRate, age = numeric.TreeDepth, mrca = T)[[1]]
    
    #############
    # GENE TREE #
    #############
    handle.GeneTree <- sim.coaltree.phylo(phy = handle.SpeciesTree, pop.size = 2)
    
    #########################
    # SIM TRAITS PREDICTORS #
    #########################
    handle.TraitData_X_S <- matrix(nrow = numeric.NumberOfSpecies, ncol = numeric.p); rownames(handle.TraitData_X_S) <- handle.SpeciesTree$tip.label; colnames(handle.TraitData_X_S) <- paste0("p_", 1:numeric.p)
    handle.TraitData_X_G <- handle.TraitData_X_S; rownames(handle.TraitData_X_G) <- handle.GeneTree$tip.label;
    for (p in 1:numeric.p){
      handle.TraitData_X_S[,p] <- mvrnorm(n = 1, mu = rep(0, numeric.NumberOfSpecies), Sigma = vcv(handle.SpeciesTree));
      handle.TraitData_X_G[,p] <- mvrnorm(n = 1, mu = rep(0, numeric.NumberOfSpecies), Sigma = vcv(handle.GeneTree));
    }
    
    ######################
    # SIM RESPONSE TRAIT #
    ######################
    vector.Y_S <- handle.TraitData_X_S %*% cbind(rep(numeric.Beta, numeric.p)) + mvrnorm(n = 1, mu = rep(0, numeric.NumberOfSpecies), Sigma = vcv(handle.SpeciesTree));
    vector.Y_G <- handle.TraitData_X_G %*% cbind(rep(numeric.Beta, numeric.p)) + mvrnorm(n = 1, mu = rep(0, numeric.NumberOfSpecies), Sigma = vcv(handle.GeneTree));
    
    handle.RAW_TRAITS_S <- cbind(vector.Y_S, handle.TraitData_X_S); colnames(handle.RAW_TRAITS_S)[1] <- "Y"
    handle.RAW_TRAITS_G <- cbind(vector.Y_G, handle.TraitData_X_G); colnames(handle.RAW_TRAITS_G)[1] <- "Y"
    
    ################
    # COMPUTE PICS #
    ################
    handle.PICS_SS <-  matrix(nrow = (numeric.NumberOfSpecies-1), ncol = (1+numeric.p)); colnames(handle.PICS_SS) <-  c("Y", paste0("p_", 1:numeric.p))
    handle.PICS_SG <- handle.PICS_GS <- handle.PICS_GG <- handle.PICS_SS
    
    for (p in 1:(1+numeric.p)){
      handle.PICS_SS[,p] <- pic(x = handle.RAW_TRAITS_S[,p], phy = handle.SpeciesTree)
      handle.PICS_SG[,p] <- pic(x = handle.RAW_TRAITS_S[,p], phy = handle.GeneTree)
      handle.PICS_GG[,p] <- pic(x = handle.RAW_TRAITS_G[,p], phy = handle.GeneTree)
      handle.PICS_GS[,p] <- pic(x = handle.RAW_TRAITS_G[,p], phy = handle.SpeciesTree)
    }

    ##########################
    # get treespace distance #
    ##########################
    handle.SpeciesTree_GeneTree <- list(); class(handle.SpeciesTree_GeneTree) <- "multiPhylo"; handle.SpeciesTree_GeneTree[[1]] <- handle.SpeciesTree; handle.SpeciesTree_GeneTree[[2]] <- handle.GeneTree
    
    handle.PICS_SS <- data.frame(handle.PICS_SS); handle.PICS_SG <- data.frame(handle.PICS_SG); handle.PICS_GG <- data.frame(handle.PICS_GG); handle.PICS_GS <- data.frame(handle.PICS_GS); 
    
    ##############
    # FIT MODELS #
    ##############
    MODEL_SS <- lm(Y ~ . - 1, data = handle.PICS_SS); MODEL_SS_SUM <- summary(MODEL_SS)
    MODEL_SG <- lm(Y ~ . - 1, data = handle.PICS_SG); MODEL_SG_SUM <- summary(MODEL_SG)
    MODEL_GG <- lm(Y ~ . - 1, data = handle.PICS_GG); MODEL_GG_SUM <- summary(MODEL_GG)
    MODEL_GS <- lm(Y ~ . - 1, data = handle.PICS_GS); MODEL_GS_SUM <- summary(MODEL_GS)
    
    ###############
    # get results #
    ###############
    numeric.P_SS <- 1-pf(q = MODEL_SS_SUM$fstatistic, MODEL_SS_SUM$df[1], MODEL_SS_SUM$df[2])[1]
    numeric.P_SG <- 1-pf(q = MODEL_SG_SUM$fstatistic, MODEL_SG_SUM$df[1], MODEL_SG_SUM$df[2])[1]
    numeric.P_GG <- 1-pf(q = MODEL_GG_SUM$fstatistic, MODEL_GG_SUM$df[1], MODEL_GG_SUM$df[2])[1]
    numeric.P_GS <- 1-pf(q = MODEL_GS_SUM$fstatistic, MODEL_GS_SUM$df[1], MODEL_GS_SUM$df[2])[1]
    
    numeric.MedBeta_SS <- median(MODEL_SS_SUM$coefficients[,1])
    numeric.MedBeta_SG <- median(MODEL_SG_SUM$coefficients[,1])
    numeric.MedBeta_GG <- median(MODEL_GG_SUM$coefficients[,1])
    numeric.MedBeta_GS <- median(MODEL_GS_SUM$coefficients[,1])
    
    ################
    # ROBUST RESID #
    ################
    MODEL_RR_SS <-  lm_robust(Y ~ .-1, data = handle.PICS_SS, se_type = "HC3"); numeric.RR_P_SS <- 1-pf(q = MODEL_RR_SS$fstatistic, MODEL_RR_SS$fstatistic[2], MODEL_RR_SS$fstatistic[3])[1]
    MODEL_RR_SG <-  lm_robust(Y ~ .-1, data = handle.PICS_SG, se_type = "HC3"); numeric.RR_P_SG <- 1-pf(q = MODEL_RR_SG$fstatistic, MODEL_RR_SG$fstatistic[2], MODEL_RR_SG$fstatistic[3])[1]
    MODEL_RR_GG <-  lm_robust(Y ~ .-1, data = handle.PICS_GG, se_type = "HC3"); numeric.RR_P_GG <- 1-pf(q = MODEL_RR_GG$fstatistic, MODEL_RR_GG$fstatistic[2], MODEL_RR_GG$fstatistic[3])[1]
    MODEL_RR_GS <-  lm_robust(Y ~ .-1, data = handle.PICS_GS, se_type = "HC3"); numeric.RR_P_GS <- 1-pf(q = MODEL_RR_GS$fstatistic, MODEL_RR_GS$fstatistic[2], MODEL_RR_GS$fstatistic[3])[1]

    #########
    # PRDTR #
    #########
    list.Model01_BM <- list(handle.Phylogeny = handle.SpeciesTree, string.Model = "BM",vector.Z = rep(0, numeric.NumberOfSpecies), vector.Theta = vector.Model_01_Theta)
    list.Model02_BM <- list(handle.Phylogeny = handle.GeneTree,    string.Model = "BM",vector.Z = rep(0, numeric.NumberOfSpecies), vector.Theta = vector.Model_01_Theta)
    
    ###################
    # COLLECT RESULTS #
    ###################
    vector.RESULTS <- c(numeric.NumberOfSpecies, numeric.Lambda, numeric.DeathRate, numeric.Beta, numeric.p,
                        RF.dist(handle.SpeciesTree, handle.GeneTree, rooted = T, normalize = T), dist.multiPhylo(x = handle.SpeciesTree_GeneTree, method = "geodesic")[1], Function_ComputeHellingerDistances(list.Model_01 = list.Model01_BM, list.Model02_BM, boo.SortNames = T)[1], mean(handle.SpeciesTree$edge.length), mean(handle.GeneTree$edge.length), get.rooted.tree.height(handle.SpeciesTree), get.rooted.tree.height(handle.GeneTree),
                        numeric.P_SS,       numeric.P_SG,        numeric.P_GG,       numeric.P_GS,
                  numeric.MedBeta_SS, numeric.MedBeta_SG,  numeric.MedBeta_GG, numeric.MedBeta_GS,
                        MODEL_SS_SUM$sigma,             MODEL_SG_SUM$sigma,              MODEL_GG_SUM$sigma,             MODEL_GS_SUM$sigma,
                 logLik(MODEL_SS)[1],            logLik(MODEL_SG)[1],             logLik(MODEL_GG)[1],            logLik(MODEL_GS)[1], 
                   mean(MODEL_SS$residuals^2),     mean(MODEL_SG$residuals^2),      mean(MODEL_GG$residuals^2),     mean(MODEL_GS$residuals^2), 
                 numeric.RR_P_SS,                numeric.RR_P_SG,                 numeric.RR_P_GG,                numeric.RR_P_GS)
    
    ######################
    # APPEND SIGNFICANCE #
    ######################
    numeric.Counter <- numeric.Counter + 1
    matrix.RESULTS[numeric.Counter, ] <- vector.RESULTS
 

    
  }
}




