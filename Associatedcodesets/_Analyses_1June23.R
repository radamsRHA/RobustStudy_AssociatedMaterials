################
# load depends #
################
library(ape)
library(geiger)
library(dplyr)
library(phytools)
library(ggplot2)
library(reshape2)
library(moments)
library(phangorn)
library(robust)

###################
# read input data #
###################
handle.ExpressionData <- as_tibble(read.table(file = 'FILTERED_ProcessedGeneExpressionData_23Sept21.txt', header = T))
handle.ExpressionData[,2:97] <- log(x = (handle.ExpressionData[,2:97] + 1), base = 10)
handle.SpeciesTree_MPEST <- geiger::rescale(x = handle.SpeciesTree_MPEST, model = "depth")
handle.SpeciesTree_MPEST <- handle.SpeciesTree_MPEST(depth = 1)
vector.TipNames <- handle.SpeciesTree_MPEST$tip.label

##############
# read input #
##############
numeric.NumberOfGenes <- length(handle.ExpressionData$human)
string.TISSUE <- "kidney"

###########
# RESULTS #
###########
matrix.RESULTS <- matrix(nrow = numeric.NumberOfGenes, ncol = 19)
colnames(matrix.RESULTS) <- c("ID", "P_SpeciesTree", "P_GeneTree_AA", "P_GeneTree_NT", "P_Delta_AA", "P_Delta_NT", 
                              "B_SpeciesTree", "B_GeneTree_AA", "B_GeneTree_NT", 
                              "Chromosome", "Start", "End", "Mid", "MM_P_SpeciesTree", "MM_P_AA", "MM_P_NT", "LnL_ST", "LnL_NT", "LnL_AA")

matrix.RESULTS_TREE <- matrix(nrow = numeric.NumberOfGenes, ncol = 19)
colnames(matrix.RESULTS_TREE) <- c("ID", "AA_Mean", "AA_SD","AA_Min", "AA_Max","AA_Med", "AA_Skew","AA_Kurt","NT_Mean", "NT_SD","NT_Min", "NT_Max","NT_Med", "NT_Skew", "NT_Kurt","RF_AA", "RF_NT", 
                                   "AA_0", "NT_0")

matrix.PICS_SpeciesTree <- matrix(nrow = numeric.NumberOfGenes, ncol = 17)
matrix.PICS_GeneTree_AA <- matrix(nrow = numeric.NumberOfGenes, ncol = 17)
matrix.PICS_GeneTree_NT <- matrix(nrow = numeric.NumberOfGenes, ncol = 17)

######################
# loop through genes #
######################
for (i in 1:numeric.NumberOfGenes){
  
  #########
  # CLEAN #
  #########
  vector.Data_Gene <- NA
  
  #####################
  # extract gene data #
  #####################
  vector.Data_Gene <- handle.ExpressionData[i,]
  
  ############################
  # filter to correct tissue #
  ############################
  vector.GeneData_Tissue <- vector.Data_Gene[grep(pattern = string.TISSUE, x = names(vector.Data_Gene))]
  vector.GeneData_Tissue_Male <- vector.GeneData_Tissue[grep(pattern = "_male", x = names(vector.GeneData_Tissue))]
  vector.GeneData_Tissue_Female <- vector.GeneData_Tissue[grep(pattern = "_female", x = names(vector.GeneData_Tissue))]
  
  ##############
  # TRIM NAMES #
  ##############
  names(vector.GeneData_Tissue_Male) = gsub(pattern = paste0("_", string.TISSUE, "_male"), replacement = "", x = names(vector.GeneData_Tissue_Male))
  names(vector.GeneData_Tissue_Female) = gsub(pattern = paste0("_", string.TISSUE, "_female"), replacement = "", x = names(vector.GeneData_Tissue_Female))
  
  ##############
  # SORT NAMES #
  ##############
  vector.GeneData_Tissue_Male <- vector.GeneData_Tissue_Male[names(vector.GeneData_Tissue_Male) == vector.TipNames]
  vector.GeneData_Tissue_Female <- vector.GeneData_Tissue_Female[names(vector.GeneData_Tissue_Female) == vector.TipNames]
  
  ###########################
  # CHECK MISSING EXPESSION #
  ###########################
  if (length(vector.GeneData_Tissue_Male[vector.GeneData_Tissue_Male==0]) <= 1){
    
    #############
    # GET TREES #
    #############
    handle.GeneTree_AA <- read.tree(text = vector.Data_Gene$Tree_AA)
    handle.GeneTree_NT <- read.tree(text = vector.Data_Gene$Tree_DNA)
    
    ############################
    # summarize gene tree data #
    ############################
    vector.BS_AminoAcid <- as.numeric(handle.GeneTree_AA$node.label)
    vector.BS_AminoAcid <- vector.BS_AminoAcid[!is.na(vector.BS_AminoAcid)]
    vector.BS_Nucleotide <- as.numeric(handle.GeneTree_NT$node.label)
    vector.BS_Nucleotide <- vector.BS_Nucleotide[!is.na(vector.BS_Nucleotide)]
    
    ###################  
    # summarize trees #
    ###################
    handle.GeneTree_AA <- midpoint.root(handle.GeneTree_AA)
    handle.GeneTree_AA <- geiger::rescale(handle.GeneTree_AA, "depth")
    handle.GeneTree_AA <- handle.GeneTree_AA(depth = 1)
    handle.GeneTree_AA$node.label <- NULL
    
    handle.GeneTree_NT <- midpoint.root(handle.GeneTree_NT)
    handle.GeneTree_NT <- geiger::rescale(handle.GeneTree_NT, "depth")
    handle.GeneTree_NT <- handle.GeneTree_NT(depth = 1)
    handle.GeneTree_NT$node.label <- NULL
    
    ##############################
    # COMPUTE PICS: SPECIES TREE #
    ##############################
    vector.PIC_Male_SpeciesTree <- pic(x = vector.GeneData_Tissue_Male, phy = handle.SpeciesTree_MPEST)
    vector.PIC_FeMale_SpeciesTree <- pic(x = vector.GeneData_Tissue_Female, phy = handle.SpeciesTree_MPEST)
    
    matrix.PICS_SpeciesTree[i,] <- c(vector.Data_Gene$human, vector.PIC_Male_SpeciesTree, vector.PIC_FeMale_SpeciesTree)
    
    #############################
    # COMPUTE PICS: GENETREE_AA #
    #############################
    vector.PIC_Male_GeneTree_AA <- pic(x = vector.GeneData_Tissue_Male, phy = handle.GeneTree_AA)
    vector.PIC_FeMale_GeneTree_AA <- pic(x = vector.GeneData_Tissue_Female, phy = handle.GeneTree_AA)
    
    matrix.PICS_GeneTree_AA[i,] <- c(vector.Data_Gene$human, vector.PIC_Male_GeneTree_AA, vector.PIC_FeMale_GeneTree_AA)
    
    #############################
    # COMPUTE PICS: GENETREE_NT #
    #############################
    vector.PIC_Male_GeneTree_NT <- pic(x = vector.GeneData_Tissue_Male, phy = handle.GeneTree_NT)
    vector.PIC_FeMale_GeneTree_NT <- pic(x = vector.GeneData_Tissue_Female, phy = handle.GeneTree_NT)
    
    matrix.PICS_GeneTree_NT[i,] <- c(vector.Data_Gene$human, vector.PIC_Male_GeneTree_NT, vector.PIC_FeMale_GeneTree_NT)
    
    #########################
    # RUN REGRESSION MODELS #
    #########################
    handle.Model_Regression_SpeciesTree = lm(vector.PIC_FeMale_SpeciesTree ~ vector.PIC_Male_SpeciesTree-1)
    handle.Model_Regression_GeneTree_AA = lm(vector.PIC_FeMale_GeneTree_AA ~ vector.PIC_Male_GeneTree_AA-1)
    handle.Model_Regression_GeneTree_NT = lm(vector.PIC_FeMale_GeneTree_NT ~ vector.PIC_Male_GeneTree_NT-1)
    
    #####################
    # robust regression #
    #####################
    handle.MM_SpeciesTree <- lmRob(vector.PIC_FeMale_SpeciesTree ~ vector.PIC_Male_SpeciesTree-1)
    handle.MM_AA <- lmRob(vector.PIC_FeMale_GeneTree_AA ~ vector.PIC_Male_GeneTree_AA-1)
    handle.MM_NT <- lmRob(vector.PIC_FeMale_GeneTree_NT ~ vector.PIC_Male_GeneTree_NT-1)
    
    ###########
    # RESULTS #
    ###########
    matrix.RESULTS[i,] <- c(vector.Data_Gene$human, 
                            summary(handle.Model_Regression_SpeciesTree)$coefficients[1,4], 
                            summary(handle.Model_Regression_GeneTree_AA)$coefficients[1,4], 
                            summary(handle.Model_Regression_GeneTree_NT)$coefficients[1,4],
                            summary(handle.Model_Regression_GeneTree_AA)$coefficients[1,4] - summary(handle.Model_Regression_SpeciesTree)$coefficients[1,4],
                            summary(handle.Model_Regression_GeneTree_NT)$coefficients[1,4] - summary(handle.Model_Regression_SpeciesTree)$coefficients[1,4],
                            summary(handle.Model_Regression_SpeciesTree)$coefficients[1,1], 
                            summary(handle.Model_Regression_GeneTree_AA)$coefficients[1,1], 
                            summary(handle.Model_Regression_GeneTree_NT)$coefficients[1,1],
                            vector.Data_Gene$Chromosome,vector.Data_Gene$StartPosition,vector.Data_Gene$EndPosition,vector.Data_Gene$StartPosition - vector.Data_Gene$EndPosition, 
                            summary(handle.MM_SpeciesTree)$coefficients[1,4], summary(handle.MM_AA)$coefficients[1,4], summary(handle.MM_NT)$coefficients[1,4], 
                            logLik(handle.Model_Regression_SpeciesTree)[1], logLik(handle.Model_Regression_GeneTree_NT)[1], logLik(handle.Model_Regression_GeneTree_AA)[1])
    
    ###################
    # Append to trees #
    ###################
    matrix.RESULTS_TREE[i,] <- c(vector.Data_Gene$human, 
                             mean(vector.BS_AminoAcid), sd(vector.BS_AminoAcid), min(vector.BS_AminoAcid), max(vector.BS_AminoAcid), median(vector.BS_AminoAcid), skewness(vector.BS_AminoAcid), kurtosis(vector.BS_AminoAcid),
                             mean(vector.BS_Nucleotide), sd(vector.BS_Nucleotide), min(vector.BS_Nucleotide), max(vector.BS_Nucleotide), median(vector.BS_Nucleotide), skewness(vector.BS_Nucleotide), kurtosis(vector.BS_Nucleotide),
                             RF.dist(handle.GeneTree_AA, handle.SpeciesTree_MPEST, normalize = T, rooted = T), RF.dist(tree1 = handle.GeneTree_NT, handle.SpeciesTree_MPEST, normalize = T, rooted = T), 
                             length(vector.BS_AminoAcid[vector.BS_AminoAcid==0]), length(vector.BS_Nucleotide[vector.BS_Nucleotide==0]))

  }
}


########
# PLOT #
########
matrix.PICS_SpeciesTree <- as_tibble(matrix.PICS_SpeciesTree)
matrix.PICS_GeneTree_AA <- as_tibble(matrix.PICS_GeneTree_AA)
matrix.PICS_GeneTree_NT <- as_tibble(matrix.PICS_GeneTree_NT)
colnames(matrix.PICS_SpeciesTree) <- c("ST_ID", paste0("ST_male_", 1:8), paste0("ST_female_", 1:8))
colnames(matrix.PICS_GeneTree_AA) <- c("AA_ID", paste0("AA_male_", 1:8), paste0("AA_female_", 1:8))
colnames(matrix.PICS_GeneTree_NT) <- c("NT_ID", paste0("NT_male_", 1:8), paste0("NT_female_", 1:8))


matrix.RESULTS <- as_tibble(matrix.RESULTS)
matrix.RESULTS2 <- as_tibble(matrix.RESULTS)
matrix.RESULTS2$P_SpeciesTree <- as.numeric(matrix.RESULTS2$P_SpeciesTree)
matrix.RESULTS2$P_GeneTree_AA <- as.numeric(matrix.RESULTS2$P_GeneTree_AA)
matrix.RESULTS2$P_GeneTree_NT <- as.numeric(matrix.RESULTS2$P_GeneTree_NT)
matrix.RESULTS2$P_Delta_AA <- as.numeric(matrix.RESULTS2$P_Delta_AA)
matrix.RESULTS2$P_Delta_NT <- as.numeric(matrix.RESULTS2$P_Delta_NT)
matrix.RESULTS2$MM_P_AA <- as.numeric(matrix.RESULTS2$MM_P_AA)
matrix.RESULTS2$MM_P_NT <- as.numeric(matrix.RESULTS2$MM_P_NT)
matrix.RESULTS2$MM_P_SpeciesTree <- as.numeric(matrix.RESULTS2$MM_P_SpeciesTree)


#matrix.RESULTS3 <- melt(data = matrix.RESULTS2, id.vars = colnames(matrix.RESULTS2)[c(-2:-9)], measure.vars = colnames(matrix.RESULTS2)[c(-1, -8, -9, -10, -11, -12, -13)],)
matrix.RESULTS3 <- melt(data = matrix.RESULTS2, id.vars = colnames(matrix.RESULTS2)[c(-2:-9)], measure.vars = colnames(matrix.RESULTS2)[c(-1,  -10, -11, -12, -13)],)
matrix.RESULTS3$value <- as.numeric(matrix.RESULTS3$value)
matrix.RESULTS3$variable <- factor(matrix.RESULTS3$variable)

matrix.RESULTS4 <- left_join(as_tibble(matrix.RESULTS_TREE), matrix.RESULTS3, "ID")

matrix.RESULTS4 %>% filter(variable %in% c("P_SpeciesTree", "P_GeneTree_AA", "P_GeneTree_NT")) %>% filter(Chromosome =="X") %>% 
  ggplot(aes(x = Mid, y = log(value), color = variable)) + geom_point() + theme_bw()  +
  scale_color_manual(values = c("P_SpeciesTree" = "red",
                                "P_GeneTree_AA" = "steelblue",
                                "P_GeneTree_NT"="darkblue")) 

###################
# plot with trees #
###################
matrix.RESULTS4 %>% filter(variable %in% c("P_Delta_NT", "P_Delta_AA")) %>% filter(Chromosome =="X") %>%
  ggplot(aes(x = AA_Mean, y = value, color = variable)) + geom_point() + theme_bw()  +
  scale_color_manual(values = c("P_Delta_NT" = "red",
                                "P_Delta_AA" = "steelblue"))

############
# examples #
############
matrix.RESULTS3 %>% filter(variable == "P_GeneTree_AA") %>% filter(log(value) < -30 )
matrix.RESULTS3 %>% filter(ID == "ENSG00000124659")

vector.TargetGene <- as.numeric(matrix.PICS_GeneTree_AA %>% filter(AA_ID == "ENSG00000124659"))
plot(vector.TargetGene[10:17] ~ vector.TargetGene[2:9])

summary(lm(vector.TargetGene[10:17] ~ vector.TargetGene[2:9]-1))
summary(lmRob(formula = vector.TargetGene[10:17] ~ vector.TargetGene[2:9]-1))
summary(lmrob(formula = vector.TargetGene[10:17] ~ vector.TargetGene[2:9]-1))
vector.TargetGene[vector.TargetGene < -1000] <- -3

######################
# explore tree stats #
######################
matrix.RESULTS_TREE <- as_tibble(matrix.RESULTS_TREE)
colnames(matrix.RESULTS_TREE)[1] <- "Tree_ID"
matrix.RESULTS_ALL <- cbind(matrix.RESULTS2, matrix.PICS_SpeciesTree, matrix.PICS_GeneTree_AA, matrix.PICS_GeneTree_NT, matrix.RESULTS_TREE)
vector.Y <- log(as.numeric(matrix.RESULTS2$MM_P_AA))
vector.X <- matrix.RESULTS_TREE$AA_0

plot(vector.Y ~ vector.X)
model <- lm(as.numeric(vector.Y) ~ as.numeric(vector.X))
abline(model)


boxplot(log(matrix.RESULTS2$P_SpeciesTree), log(matrix.RESULTS2$P_GeneTree_AA), log(matrix.RESULTS2$P_GeneTree_NT), notch = T)

##############
# FIND GENES #
##############
TARGET_EXAMPLE <- matrix.RESULTS2_HEART %>% filter(B_at >19)
TARGET_EXAMPLE <- TARGET_EXAMPLE[1,]


vector.TargetGene <- as.numeric(matrix.PICS_GeneTree_NT %>% filter(NT_ID == TARGET_EXAMPLE$ID))
plot(vector.TargetGene[10:17] ~ vector.TargetGene[2:9], main = paste0("NT_", TARGET_EXAMPLE$ID))
abline(lm(vector.TargetGene[10:17] ~ vector.TargetGene[2:9]-1))
summary(lm(vector.TargetGene[10:17] ~ vector.TargetGene[2:9]-1))

dim(matrix.RESULTS2_HEART)


#############
# summarize #
#############
dim(matrix.RESULTS2_HEART %>% filter(B_nt >-log10(0.05/length(vector.Chromosomes))))
dim(matrix.RESULTS2_HEART)

vioplot(as.numeric(matrix.RESULTS2_HEART$LnL_ST), as.numeric(matrix.RESULTS2_HEART$LnL_NT), as.numeric(matrix.RESULTS2_HEART$LnL_AA))

#######
# LNL #
#######
matrix.RESULTS2_HEART <- matrix.RESULTS2_HEART %>% mutate(Tissue = "Heart")
matrix.RESULTS2_KIDNEY <- matrix.RESULTS2_KIDNEY %>% mutate(Tissue = "Kidney")
matrix.RESULTS2_BRAIN <- matrix.RESULTS2_BRAIN %>% mutate(Tissue = "Brain")

matrix.RESULTS2_HEART$LnL_AA
vioplot(matrix.RESULTS2_HEART$LnL_AA)
matrix.RESULTS_TREE <- as_tibble(matrix.RESULTS_TREE)

vioplot(as.numeric(matrix.RESULTS2_HEART$LnL_ST), as.numeric(matrix.RESULTS2_HEART$LnL_NT), as.numeric(matrix.RESULTS2_HEART$LnL_AA), rep(0, 100),
        as.numeric(matrix.RESULTS2_KIDNEY$LnL_ST), as.numeric(matrix.RESULTS2_KIDNEY$LnL_NT), as.numeric(matrix.RESULTS2_KIDNEY$LnL_AA),rep(0, 100),
        as.numeric(matrix.RESULTS2_BRAIN$LnL_ST), as.numeric(matrix.RESULTS2_BRAIN$LnL_NT), as.numeric(matrix.RESULTS2_BRAIN$LnL_AA))


prcomp(matrix.RESULTS2_HEART$)
