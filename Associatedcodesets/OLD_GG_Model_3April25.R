################
# load depends #
################
library(geiger)
library(phybase)
library(ape)
library(phytools)
library(robust)
library(broom)
library(phangorn)
library("distory")

###############################
# Function to simulate traits #
###############################
Function.Simulate_GG_GS <- function(numeric.NumberOfSpecies, 
                                    numeric.NumberOfSpeciesTrees, 
                                    numeric.NumberOfGeneTrees, 
                                    numeric.NumberOfTraitReplicates, 
                                    numeric.Rho,
                                    numeric.NumberTraits,
                                    boo.ScaleDepth){
  
  ######################
  # Define simulations #
  ######################
  vector.RangeRate <- c(10^-4, 10^-3, 10^-2, 10^-1, 10^0, 10^1, 10^2)
  matrix.R <- matrix(numeric.Rho, nrow = numeric.NumberTraits, ncol = numeric.NumberTraits)
  diag(matrix.R) <- 1

  ##################
  # Matrix results #
  ##################
  matrix.RESULTS <- matrix(nrow = length(vector.RangeRate)*numeric.NumberOfSpeciesTrees*numeric.NumberOfGeneTrees*numeric.NumberOfTraitReplicates, ncol = 12)
  vector.ColNames <- c("BirthRate", "SpeciesTree", "GeneTree", "RF", "Geodesic", "TraitRep", "GeneTree_P", "SpeciesTree_P", "GeneTree_Beta", "SpeciesTree_Beta", "GeneTree_RSS", "SpeciesTree_RSS")
  colnames(matrix.RESULTS) <- vector.ColNames
  
  #################################
  # Loop through speciation rates #
  #################################
  counter <- 0
  for (i in 1:length(vector.RangeRate)){
    
    print(i/length(vector.RangeRate))
    
    #######################
    # Get speciation rate #
    #######################
    numeric.SpeciationRate <- NA
    numeric.SpeciationRate <- vector.RangeRate[i]
    
    ########################################
    # Loop through species tree replicates #
    ########################################
    for (j in 1:numeric.NumberOfSpeciesTrees){
      
      #########################
      # Simulate species tree #
      #########################
      handle.SpeciesTree <- NA
      handle.SpeciesTree <- sim.bdtree(b = numeric.SpeciationRate, n = numeric.NumberOfSpecies)
      
      #####################################
      # Loop through gene tree simulation #
      #####################################
      for (g in 1:numeric.NumberOfGeneTrees){
        
        ######################
        # Simulate gene tree #
        ######################
        handle.GeneTree <- NA
        handle.GeneTree <- sim.coaltree.phylo(phy = handle.SpeciesTree, pop.size = 2)
        
        ##################################
        # Loop through trait simulations #
        ##################################
        for (t in 1:numeric.NumberOfTraitReplicates){
          
          ##################################
          # Simulate trait under gene tree #
          ##################################
          vector.TraitData <- NA
          vector.TraitData <- sim.char(handle.GeneTree, par=matrix.R, model="BM", root=0)[,,1]
          
          if (boo.ScaleDepth == T){
            handle.GeneTree_SCALED <- geiger::rescale(x = handle.GeneTree, model = "depth")
            handle.GeneTree_SCALED <- handle.GeneTree_SCALED(depth = 1)
          }
          if (boo.ScaleDepth == F){
            handle.GeneTree_SCALED <- handle.GeneTree
          }
          
      	   if (boo.ScaleDepth == T){
      	     handle.SpeciesTree_SCALED <- geiger::rescale(x = handle.SpeciesTree, model = "depth")
      	     handle.SpeciesTree_SCALED <- handle.SpeciesTree_SCALED(depth = 1)
      	   }
      	   if (boo.ScaleDepth == F){
      	     handle.SpeciesTree_SCALED <- handle.SpeciesTree
      	   }
          
          ###############################
          # Compute PICs with Gene Tree #
          ###############################
          matrix.GeneTree_PICs <- NA
          matrix.GeneTree_PICs <- matrix(nrow = (numeric.NumberOfSpecies-1), ncol = numeric.NumberTraits)
          colnames(matrix.GeneTree_PICs) <- c("GeneTree_PIC_Y", paste0("GeneTree_PIC_X", 1:(numeric.NumberTraits-1)))
          
          for (p in 1:numeric.NumberTraits){
            matrix.GeneTree_PICs[,p] <- pic(x = vector.TraitData[,p], phy = handle.GeneTree_SCALED)
          }
          
          handle.GeneTree_PICs <- data.frame(matrix.GeneTree_PICs) # convert to data frame #
          
          ##################################
          # Compute PICs with Species Tree #
          ##################################
          matrix.SpeciesTree_PICs <- NA
          matrix.SpeciesTree_PICs <- matrix(nrow = (numeric.NumberOfSpecies-1), ncol = numeric.NumberTraits)
          colnames(matrix.SpeciesTree_PICs) <- c("SpeciesTree_PIC_Y", paste0("SpeciesTree_PIC_X", 1:(numeric.NumberTraits-1)))
          
          for (p in 1:numeric.NumberTraits){
            matrix.SpeciesTree_PICs[,p] <- pic(x = vector.TraitData[,p], phy = handle.SpeciesTree_SCALED)
          }
          
          handle.SpeciesTree_PICs <- data.frame(matrix.SpeciesTree_PICs) # convert to data frame #
          
          ############################
          # Run gene tree regression #
          ############################
          handle.Model_GG <- NA
          handle.Model_GG <- lm(GeneTree_PIC_Y ~ . -1, data = handle.GeneTree_PICs)
          
          ############################
          # Run species tree regression #
          ############################
          handle.Model_GS <- NA
          handle.Model_GS <- lm(SpeciesTree_PIC_Y ~ . -1, data = handle.SpeciesTree_PICs)
          
          #####################
          # Summarize results #
          #####################
          counter <- counter + 1
          
          if (numeric.NumberTraits == 2){
            
            ##########
            # get GG #
            ##########
            numeric.p_GG <- summary(handle.Model_GG)$coefficients[4]
            numeric.Beta_GG <- summary(handle.Model_GG)$coefficients[1]
            #numeric.RSS_GG <- deviance(handle.Model_GG)
            numeric.RSS_GG <- mean(handle.Model_GG$residual^2)
            
            ##########
            # get GS #
            ##########
            numeric.p_GS <- summary(handle.Model_GS)$coefficients[4]
            numeric.Beta_GS <- summary(handle.Model_GS)$coefficients[1]
            #numeric.RSS_GS <- deviance(handle.Model_GS)
			numeric.RSS_GS <- mean(handle.Model_GS$residual^2)
            
          }
          
          if (numeric.NumberTraits > 2){
            
            ##########
            # get GG #
            ##########
            numeric.p_GG <- glance(handle.Model_GG)$p.value
            numeric.Beta_GG <- mean(summary(handle.Model_GG)$coefficients[,1])
            #numeric.RSS_GG <- deviance(handle.Model_GG)
            numeric.RSS_GG <- mean(handle.Model_GG$residual^2)
            
            ##########
            # get GS #
            ##########
            numeric.p_GS <- glance(handle.Model_GS)$p.value
            numeric.Beta_GS <- mean(summary(handle.Model_GS)$coefficients[,1])
            #numeric.RSS_GS <- deviance(handle.Model_GS)
            numeric.RSS_GS <- mean(handle.Model_GS$residual^2)
            
          }
          
          ####################
          # Organize results #
          ####################
          list.TREES <- list()
          class(list.TREES) <- "multiPhylo"
          list.TREES[[1]] <- handle.SpeciesTree
          list.TREES[[2]] <- handle.GeneTree
          numeric.GEO <- dist.multiPhylo(x = list.TREES, method = "geodesic")
          numeric.RF <- RF.dist(tree1 = handle.SpeciesTree, handle.GeneTree, normalize = T, rooted = T)
          #vector.RESULTS <- c("BirthRate", "SpeciesTree", "GeneTree", "TraitRep", "GeneTree_P", "SpeciesTree_P", "GeneTree_Beta", "SpeciesTree_Beta", "GeneTree_RR", "SpeciesTree_RR")
          vector.RESULTS <- c(numeric.SpeciationRate, j, g, numeric.RF, numeric.GEO, t, numeric.p_GG, numeric.p_GS, numeric.Beta_GG, numeric.Beta_GS, numeric.RSS_GG, numeric.RSS_GS)
          matrix.RESULTS[counter,] <- vector.RESULTS
        }
      }
    }
  }
  return(matrix.RESULTS)
}
