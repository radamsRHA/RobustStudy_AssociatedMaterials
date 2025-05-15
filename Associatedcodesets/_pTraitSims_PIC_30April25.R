################
# load depends #
################
rm(list = ls()); 
library(sfsmisc); 
library(MASS); 
library(estimatr)

###############################
# set experimental parameters #
###############################
vector.Lambda <- 1
numeric.NumberOfReps <- 10^4
numeric.NumberOfSpecies <- 1000
numeric.p <- 100

MATRIX.RESULTS <- matrix(nrow = length(vector.Lambda) * numeric.NumberOfReps, ncol = 10)
colnames(MATRIX.RESULTS) <- c("LAMBDA", "REP", "NumberOfSpecies", "NumberPredictors", "Pval_GS_L2", "Pval_GS_RR", "Pval_Rand_L2", "Pval_Rand_RR","Pval_NoTree_L2", "Pval_NoTree_RR")
string.RESULTS <- "RESULTS_SpeciesXXX_TraitsYYY.txt"; string.RESULTS <- gsub(pattern = "XXX", numeric.NumberOfSpecies, x = string.RESULTS); string.RESULTS <- gsub(pattern = "YYY", numeric.p, x = string.RESULTS)

###############
# SIMULATIONS #
###############
count <- 0
for (lambda in 1:length(vector.Lambda)){
  
  # get lambda #
  print(lambda)
  numeric.Lambda <- vector.Lambda[lambda]
  
  # simulate replicates #
  for (rep in 1:numeric.NumberOfReps){
    
    TREE_SPECIES <- TREE_RAND <- DATA_PIC_GS <- DATA_PIC_RAND <- DATA_PIC_NoTree <- NULL
    TREE_SPECIES <- TreeSim::sim.bd.taxa.age(n = numeric.NumberOfSpecies, numbsim = 1, lambda = numeric.Lambda, mu = numeric.Lambda/2, age = 10, mrca = T)[[1]]
    TREE_RAND <- TreeSim::sim.bd.taxa.age(n = numeric.NumberOfSpecies, numbsim = 1, lambda = numeric.Lambda, mu = numeric.Lambda/2, age = 10, mrca = T)[[1]]
    DATA_PIC_GS <- matrix(nrow = (numeric.NumberOfSpecies-1), ncol = (numeric.p + 1))
    DATA_PIC_RAND <- matrix(nrow = (numeric.NumberOfSpecies-1), ncol = (numeric.p + 1))
    DATA_PIC_NoTree <- matrix(nrow = numeric.NumberOfSpecies, ncol = (numeric.p + 1))
    
    # simulate gene trees #
    for (gene in 1:(numeric.p +1)){
      TREE_GENE <- NA
      TREE_GENE <- phybase::sim.coaltree.phylo(phy = TREE_SPECIES, pop.size = 2, nsamples = 1)
      TRAIT_GENE <- geiger::sim.char(phy = TREE_GENE, par = 1, nsim = 1, model = "BM", root = rnorm(1))[,,1]
      TRAIT_GENE <- TRAIT_GENE[paste0("t", 1:numeric.NumberOfSpecies)]
      TRAIT_PIC <- ape::pic(x = TRAIT_GENE, phy = TREE_SPECIES)
      DATA_PIC_GS[,gene] <- TRAIT_PIC
      
      # PICS under random tree #
      TRAIT_RAND <- ape::pic(x = TRAIT_GENE, phy = TREE_RAND)
      DATA_PIC_RAND[,gene] <- TRAIT_RAND
      
      # No Tree #
      DATA_PIC_NoTree[,gene] <- TRAIT_GENE

    }
    
    # organize data frame #
    colnames(DATA_PIC_GS) <- colnames(DATA_PIC_NoTree) <- colnames(DATA_PIC_RAND) <- c("Y", paste0("X", 1:numeric.p)); DATA_PIC_GS <- data.frame(DATA_PIC_GS); DATA_PIC_RAND <- data.frame(DATA_PIC_RAND); DATA_PIC_NoTree <- data.frame(DATA_PIC_NoTree)
    
    # FIT MODELS #
    MODEL_GS_L2 <- lm(Y ~ 0 + ., DATA_PIC_GS)
    MODEL_GS_RR <- estimatr::lm_robust(Y ~ 0 + ., data = DATA_PIC_GS, se_type = "HC3")
    
    MODEL_RAND_L2 <- lm(Y ~ 0 + ., DATA_PIC_RAND)
    MODEL_RAND_RR <- estimatr::lm_robust(Y ~ 0 + ., data = DATA_PIC_RAND, se_type = "HC3")
    
    MODEL_NoTree_L2 <- lm(Y ~ ., DATA_PIC_NoTree)
    MODEL_NoTree_RR <- estimatr::lm_robust(Y ~ ., data = DATA_PIC_NoTree, se_type = "HC3")
    
    
    # SUMMARY #
    SUMMARY_MODEL_L2 <- summary(MODEL_GS_L2)
    PVAL_MODEL_L2 <- pf(SUMMARY_MODEL_L2$fstatistic[1], SUMMARY_MODEL_L2$fstatistic[2], SUMMARY_MODEL_L2$fstatistic[3], lower.tail = FALSE)
    PVAL_MODEL_RR <- pf(MODEL_GS_RR$fstatistic[1], MODEL_GS_RR$fstatistic[2], MODEL_GS_RR$fstatistic[3], lower.tail = FALSE)
    
    SUMMARY_MODEL_L2_RAND <- summary(MODEL_RAND_L2)
    PVAL_MODEL_L2_RAND <- pf(SUMMARY_MODEL_L2_RAND$fstatistic[1], SUMMARY_MODEL_L2_RAND$fstatistic[2], SUMMARY_MODEL_L2_RAND$fstatistic[3], lower.tail = FALSE)
    PVAL_MODEL_RR_RAND <- pf(MODEL_RAND_RR$fstatistic[1], MODEL_RAND_RR$fstatistic[2], MODEL_RAND_RR$fstatistic[3], lower.tail = FALSE)
    
    SUMMARY_MODEL_L2_NoTree <- summary(MODEL_NoTree_L2)
    PVAL_MODEL_L2_NoTree <- pf(SUMMARY_MODEL_L2_NoTree$fstatistic[1], SUMMARY_MODEL_L2_NoTree$fstatistic[2], SUMMARY_MODEL_L2_NoTree$fstatistic[3], lower.tail = FALSE)
    PVAL_MODEL_RR_NoTree <- pf(MODEL_NoTree_RR$fstatistic[1], MODEL_NoTree_RR$fstatistic[2], MODEL_NoTree_RR$fstatistic[3], lower.tail = FALSE)


  }
}

