### Organisation of this file ###
# Functions
  # Extract_MDFs
  # add_Cell_ID_unique <- function(CP_dataframe)
  # get_combined_dataframe
  # load_coloc_data
# Data loading
# Exploratory data analysis 
# Here some things are plot that may be used to tweak thresholding for transfection etc
# Plots of MDF
# Calculations of mito morph measurements

library(Hmisc)
library(dplyr)
library(tidyverse) # ggplot2 etc
library(bayestestR) # to calculate AUC
library(stats)



### FUNCTIONS ###
extract_MDFs <- function(inputData) {
  
  colnames <- colnames(inputData)
  EVF <- c()
  for (i in 1:length(colnames)){
    number <- as.numeric(str_extract_all(colnames[i], "\\d+")[[1]][1])
    EVF[i] <- number
  }
  
  profile <- data.frame('EVF' = c(EVF),
                        'RDF' = c(t(inputData)))
  profile$EVF <- (profile$EVF / max(profile$EVF))*100
  profile <- profile[order(profile$EVF),]
  # the 'profile' dataframe contains a profile 
  # of fractional GFP intensity (RDF) vs the fractional distance from the nucleus (EVF)
  for (i in 1:nrow(profile)){
    cum_curve <- profile[1:i,]
    AUC <- area_under_curve(cum_curve$EVF,cum_curve$RDF, method = c("trapezoid"))
    profile$cum_AUC[i] <- AUC
    #This section estimates the cumulative area under the curve
  }
  MDF <- approx(profile$cum_AUC,profile$EVF, xout = 0.5 , method = "linear")
  MDF <- MDF[[2]]
  # This section estimates (via linear interpolation) the MDF
  # The MDF represents the point at which the cumulative AUC is 50%
  return(MDF)
}

add_Cell_ID_unique <- function(CP_dataframe){
  CP_dataframe$Image_ID <- CP_dataframe$FileName_BlueOrig
  CP_dataframe$Cell_ID_unique <-  paste(CP_dataframe$Image_ID, " Cell", CP_dataframe$ObjectNumber, sep="")
  return(CP_dataframe)
}


get_combined_dataframe <- function(datafolder, condition, mitochannel){
  
  # Load in data from MitoMorphExtractor
  MitoDf <- read.csv( paste0(datafolder,"/AutoCBtiffs/MitoMorphResults_sigma03_new.txt")) 
  MitoDf$uniqueID <- MitoDf$Mito_ID
  MitoDf <- MitoDf %>% separate(Mito_ID, c('Image_ID', 'Cell_ID', 'Mito_ID'),sep="__")
  MitoDf$Mito_NR <- row.names(MitoDf) 
  MitoDf$Cell_ID_unique <- paste(MitoDf$Image_ID, MitoDf$Cell_ID)
  
  # Load in data from CP
  CP_cells <- read.csv(paste0(datafolder,"CP_measurementsAllCells.csv"))
  CP_cyto <- read.csv(paste0(datafolder,"CP_measurementsAllCytoplasm.csv"))
  
  # Append the cytoplasm only measurements to the dataframe
  colnames_cyto <- NULL
  for (i in colnames(CP_cyto)){
    new <- paste0("Cyto_", i)
    colnames_cyto <- c(colnames_cyto, new)
  }
  colnames(CP_cyto) <- colnames_cyto
  CP_cells <- cbind(CP_cells, CP_cyto[, 12:ncol(CP_cyto)])
  CP_cells <- add_Cell_ID_unique(CP_cells)
  df <- CP_cells
  
  # Fix lost cells and subset 
  cells_united <- intersect(unique(MitoDf$Cell_ID_unique), df$Cell_ID_unique)
  cells <- cells_united
  df <- df[which(df$Cell_ID_unique %in% cells_united), ]
  MitoDf <- MitoDf[which(MitoDf$Cell_ID_unique %in% cells_united), ]
  
  # Add conditions are used
  df$condition <- condition
  
  # Add the integrated intensity of the DNA in nucleus and the nuclear GFP signal to df
  CP_nuc <- read.csv(paste0(datafolder,"CP_measurementsfilteredNuclei.csv"))
  CP_nuc <- add_Cell_ID_unique(CP_nuc)
  CP_nuc <- CP_nuc[which(CP_nuc$Cell_ID_unique %in% cells_united), ]
  for (i in 1:nrow(CP_nuc)){
    df$nuclear_Intensity_IntegratedIntensity_BlueOrig[i] <- CP_nuc[which(df$Cell_ID_unique == CP_nuc$Cell_ID_unique[i]), "Intensity_IntegratedIntensity_BlueOrig"]
    df$nuclear_Intensity_MeanIntensity_GreenOrig[i] <- CP_nuc[which(df$Cell_ID_unique == CP_nuc$Cell_ID_unique[i]), "Intensity_MeanIntensity_GreenOrig"]
  }
  
  # Choose how to classify transfected cells
  transf <- df[which(df$Children_TransfectedCells_Count > 0), "Cell_ID_unique"] # based on the CellProfiler classification 
  thresholdTransfected = 1.25*median(df$Intensity_MeanIntensity_GreenOrig) # ARBITRARY
  #transf <- df[which(df$Intensity_MeanIntensity_GreenOrig > thresholdTransfected), "Cell_ID_unique"] 
  untransf <- df$Cell_ID_unique[which(df$Cell_ID_unique %nin% transf)]
  df$transfected <- rep("untransfected")
  df$transfected[which(df$Cell_ID_unique %in% transf)] <- "transfected"
  
  
  # Calculate MDF
  cols <- colnames(df)
  FracDf <- df[ , grepl("RadialDistribution_FracAtD" ,cols) & grepl(mitochannel, cols)]
  RDFs <- FracDf
  for (i in 1:nrow(df)){
    df$MDF_Mito_mean[i] <- as.numeric(extract_MDFs(RDFs[i,]))
    df$MDF_Mito[i] <- as.numeric(extract_MDFs(RDFs[i,]))
  }
  
  # Add transfected status to MitoDf
  for (i in 1:nrow(MitoDf)){
    MitoDf$transfected[i] <- df[which(df$Cell_ID_unique == MitoDf$Cell_ID_unique[i]), "transfected"]
    MitoDf$condition[i] <- df[which(df$Cell_ID_unique == MitoDf$Cell_ID_unique[i]), "condition"]
    MitoDf$Intensity_MeanIntensity_GreenOrig[i] <- df[which(df$Cell_ID_unique == MitoDf$Cell_ID_unique[i]), "Intensity_MeanIntensity_GreenOrig"]
  }
  
  # Total voxels per mito = nr_Branches + Average_Branch_Length + Junction_voxels + End_points
  MitoDf$nr_Total_voxels = MitoDf$nr_Branches*MitoDf$Average_Branch_Length+MitoDf$nr_Junction_voxels+MitoDf$nr_End.point_voxels
  
  # Define some ratios that can be used to further classify the mitochondrial morphology
  MitoDf$ratio_Total_Junctions <- MitoDf$nr_Total_voxels/MitoDf$nr_Junctions
  MitoDf$ratio_Junctions_Endpoints <- (MitoDf$nr_Junctions) / (MitoDf$nr_End.point_voxels)
  MitoDf$adapted_ratio_Junctions_Endpoints <- (MitoDf$nr_End.point_voxels-2)/(MitoDf$nr_Junctions+1)
  
  
  # Classify Mitos complex, according to size and adapted ratio #SUPER ARBTIRARY CAN BE CHANGED
  tubular_rownrs <-  which(MitoDf$nr_Total_voxels >= 100)
  fragmented_rownrs <- which(MitoDf$nr_Total_voxels < 10)

  fused_rownrs <- which(MitoDf$adapted_ratio_Junctions_Endpoints > 0 & MitoDf$adapted_ratio_Junctions_Endpoints <= 0.5 | MitoDf$adapted_ratio_Junctions_Endpoints<0)
  branched_rownrs <- which(MitoDf$adapted_ratio_Junctions_Endpoints > 0.5)
  ringshaped_rownrs <- which( MitoDf$adapted_ratio_Junctions_Endpoints < 0)
  classic_tubular_rownrs <- which(MitoDf$adapted_ratio_Junctions_Endpoints == 0 & MitoDf$nr_Junctions==0)
  

  MitoDf$Mito_classification_simple <- "Intermediate"
  MitoDf$Mito_classification_simple[tubular_rownrs] <- "Tubular (=>100 voxels)"
  MitoDf$Mito_classification_simple[fragmented_rownrs] <- "Fragmented (<10 voxels)"

  MitoDf$Mito_classification <- "Fused ( 0 < AR <= 0.5)"
  MitoDf$Mito_classification[fused_rownrs] <- "Fused ( 0 < AR <= 0.5)"
  MitoDf$Mito_classification[branched_rownrs] <- "Branched (AR > 0.5)"
  MitoDf$Mito_classification[classic_tubular_rownrs] <- "Classic (AR == 0)"
  MitoDf$Mito_classification[ringshaped_rownrs] <- "Hyperfused (AR < 0)"
  
  
  MitoDf$Mito_classification_size <- "1-5 voxels"
  MitoDf$Mito_classification_size[which(MitoDf$nr_Total_voxels >= 5 & MitoDf$nr_Total_voxels <=10) ] <- "5-10 voxels"
  MitoDf$Mito_classification_size[which(MitoDf$nr_Total_voxels >10 & MitoDf$nr_Total_voxels <=30) ] <- "10-30 voxels"
  MitoDf$Mito_classification_size[which(MitoDf$nr_Total_voxels >30 & MitoDf$nr_Total_voxels <=50)] <- "30-50 voxels "
  MitoDf$Mito_classification_size[which(MitoDf$nr_Total_voxels >50 & MitoDf$nr_Total_voxels <=100)] <- "50-100 voxels"
  MitoDf$Mito_classification_size[which(MitoDf$nr_Total_voxels >100)] <- ">100 voxels"
  
  cellDf <- data.frame(unique(df$Cell_ID_unique))
  colnames(cellDf)<- "Cell_ID_unique"
  
  cells=unique(MitoDf$Cell_ID_unique)
  for(cell in cells){
    voxelsizes <- MitoDf$nr_Total_voxels[which(MitoDf$Cell_ID_unique==cell)]
    cellDf$largest_mitochondrion_perCell[which(cellDf$Cell_ID_unique == cell)] <- max(voxelsizes)
    cellDf$Summed_Total_voxels_perCell[which(cellDf$Cell_ID_unique == cell)] <- sum(voxelsizes)
    cellDf$number_of_mitochondria_perCell[which(cellDf$Cell_ID_unique == cell)] <- length(voxelsizes)
    cellDf$variance_mitochondrial_size_perCell[which(cellDf$Cell_ID_unique == cell)] <- var(voxelsizes)
    cellDf$variance_average_branchLength_perCell[which(cellDf$Cell_ID_unique == cell)]<- var(MitoDf$Average_Branch_Length[which(MitoDf$Cell_ID_unique==cell)])
    cellDf$nr_mitos_over50vox[which(cellDf$Cell_ID_unique == cell)] <- length(which(voxelsizes>50))
    
    cellDf$summed_voxels_tubular_perCell[which(cellDf$Cell_ID_unique == cell)]<- sum(MitoDf$nr_Total_voxels[which(MitoDf$Cell_ID_unique ==cell & MitoDf$Mito_classification_simple =="Tubular (=>100 voxels)")])
    cellDf$summed_voxels_intermediate_perCell[which(cellDf$Cell_ID_unique == cell)] <- sum(MitoDf$nr_Total_voxels[which(MitoDf$Cell_ID_unique ==cell & MitoDf$Mito_classification_simple =="Intermediate")])
    cellDf$summed_voxels_fragmented_perCell[which(cellDf$Cell_ID_unique == cell)] <- sum(MitoDf$nr_Total_voxels[which(MitoDf$Cell_ID_unique ==cell & MitoDf$Mito_classification_simple =="Fragmented (<10 voxels)")])

    cellDf$percentage_tubular_voxels_perCell[which(cellDf$Cell_ID_unique == cell)] <-  (cellDf$summed_voxels_tubular_perCell[which(cellDf$Cell_ID_unique == cell)] /cellDf$Summed_Total_voxels_perCell[which(cellDf$Cell_ID_unique == cell)]*100)
    cellDf$percentage_intermediate_voxels_perCell[which(cellDf$Cell_ID_unique == cell)] <-  (cellDf$summed_voxels_intermediate_perCell[which(cellDf$Cell_ID_unique == cell)] /cellDf$Summed_Total_voxels_perCell[which(cellDf$Cell_ID_unique == cell)]*100)
    cellDf$percentage_fragmented_voxels_perCell[which(cellDf$Cell_ID_unique == cell)] <-  (cellDf$summed_voxels_fragmented_perCell[which(cellDf$Cell_ID_unique == cell)] /cellDf$Summed_Total_voxels_perCell[which(cellDf$Cell_ID_unique == cell)] *100)
    
    cellDf$summed_voxels_fused_perCell[which(cellDf$Cell_ID_unique == cell)]<- sum(MitoDf$nr_Total_voxels[which(MitoDf$Cell_ID_unique ==cell & MitoDf$Mito_classification =="Fused ( 0 < AR <= 0.5)")]) 
    cellDf$summed_voxels_branched_perCell[which(cellDf$Cell_ID_unique == cell)] <- sum(MitoDf$nr_Total_voxels[which(MitoDf$Cell_ID_unique ==cell & MitoDf$Mito_classification =="Branched (AR > 0.5)")])
    cellDf$summed_voxels_classic_perCell[which(cellDf$Cell_ID_unique == cell)] <- sum(MitoDf$nr_Total_voxels[which(MitoDf$Cell_ID_unique ==cell & MitoDf$Mito_classification =="Classic (AR == 0)")])
                                                                                  
    cellDf$percentage_fused_voxels_perCell[which(cellDf$Cell_ID_unique == cell)] <-  (cellDf$summed_voxels_fused_perCell[which(cellDf$Cell_ID_unique == cell)] /cellDf$Summed_Total_voxels_perCell[which(cellDf$Cell_ID_unique == cell)]*100)
    cellDf$percentage_branched_voxels_perCell[which(cellDf$Cell_ID_unique == cell)] <-  (cellDf$summed_voxels_branched_perCell[which(cellDf$Cell_ID_unique == cell)] /cellDf$Summed_Total_voxels_perCell[which(cellDf$Cell_ID_unique == cell)]*100)
    cellDf$percentage_classic_voxels_perCell[which(cellDf$Cell_ID_unique == cell)] <-  (cellDf$summed_voxels_classic_perCell[which(cellDf$Cell_ID_unique == cell)] /cellDf$Summed_Total_voxels_perCell[which(cellDf$Cell_ID_unique == cell)] *100)
    
  }
  
  # Merge the new parameters with the other df
  df <- merge(df, cellDf, by="Cell_ID_unique")
  MitoDf <- merge(MitoDf, cellDf, by="Cell_ID_unique")
  
  
  return(list(df, MitoDf))
}


load_coloc_data <- function(datafolder){
  lines <- readLines(paste0(datafolder, "colocFile.txt"))
  IDs <- grep("^DAPI.+__", lines, value=TRUE)
  coloc <- data.frame(IDs)
  M1 <- grep("M1=", lines,  value=TRUE) # M1 = fraction of A(Mito) overlapping B(LD)
  M2 <- grep("M2=", lines,  value=TRUE)
  coloc$M1coef <- M1[lapply(range(1,length(M1)), "%%", 2) == 0]
  coloc$M2coef <- M2[lapply(range(1,length(M2)), "%%", 2) == 0]
  coloc$M1coef <- as.numeric(gsub("[^[:digit:]\\.[:digit:]+]","",coloc$M1))-10
  coloc$M2coef <- as.numeric(gsub("[^[:digit:]\\.[:digit:]+]","",coloc$M2))-20
  coloc <- coloc %>% separate(IDs, c('Image_ID', 'Cell_ID'), sep="__")
  coloc$Cell_ID_unique <- paste(coloc$Image_ID, coloc$Cell_ID)
  
  return(coloc)
  
}


####
# LOAD DATA this example shows how to load (and combine different) outputs
datafolder = "/Volumes/ExtDriveLou/LAB21/"
data <- get_combined_dataframe(datafolder, "null", "Farr")

df <- data[[1]]
MitoDf <- data[[2]]

outfolder <- paste0(datafolder, "Routput")
dir.create(outfolder)





# Load coloc data (if applicable)
coloc <- load_coloc_data(paste0(datafolder, "AutoCBtiffs/"))
unique(coloc$Image_ID)
unique(df$Image_ID)
df$Cell_ID_unique[which(unique(df$Cell_ID_unique) %nin% coloc$Cell_ID_unique)]
df <- df[which(unique(df$Cell_ID_unique) %in% unique(coloc$Cell_ID_unique)), ]
MitoDf <- MitoDf[which(unique(MitoDf$Cell_ID_unique) %in% unique(coloc$Cell_ID_unique)), ]
for (i in 1:nrow(df)){
  df$M1coef[i] <- coloc$M1coef[ which(coloc$Cell_ID_unique == df$Cell_ID_unique[i]) ]
  df$M2coef[i] <- coloc$M2coef[ which(coloc$Cell_ID_unique == df$Cell_ID_unique[i]) ]
}



# Fix conditions (example, change according to your own)
PLIN5_Oleate <- grep("PLIN5-Oleate_LT", df$FileName_BlueOrig, value=FALSE) #grep the pattern that you named your condition
PLIN5_minusOleate <- grep("PLIN5-minusOleate_LT", df$FileName_BlueOrig , value=FALSE)
df$condition[PLIN5_Oleate] <- "PLIN5_Oleate"
df$condition[PLIN5_minusOleate] <- "PLIN5_minusOleate"

df$condition2 <- "nul"
df$condition2[which(df$condition == "PLIN5_Oleate" & df$transfected =="transfected")] <- "PLIN5_Oleate transfected"
df$condition2[which(df$condition == "PLIN5_Oleate" & df$transfected =="untransfected")] <- "PLIN5_Oleate untransfected"
df$condition2[which(df$condition == "PLIN5_minusOleate" & df$transfected =="transfected")] <- "PLIN5_minusOleate transfected"
df$condition2[which(df$condition == "PLIN5_minusOleate" & df$transfected =="untransfected")] <- "PLIN5_minusOleate untransfected"

MitoDf$condition <- "null"
MitoDf$condition2 <- "null"
for (i in 1:nrow(MitoDf)){
  MitoDf$condition[i] <- df$condition[ which(df$Cell_ID_unique == MitoDf$Cell_ID_unique[i]) ]
  MitoDf$condition2[i] <- df$condition2[ which(df$Cell_ID_unique == MitoDf$Cell_ID_unique[i]) ]
}
df %>% count(df$condition2)
MitoDf %>% count(MitoDf$condition2)



my_comparisons <- list( c("PLIN5_Oleate", "PLIN5_minusOleate")) #this is needed for significance
df <- df[which(df$condition %in% c("PLIN5_Oleate", "PLIN5_minusOleate")), ]
MitoDf <- MitoDf[which(MitoDf$condition %in% c("PLIN5_Oleate", "PLIN5_minusOleate")), ]










# Below is some plots code feel free to use or write your own
### Exploratory data analysis ###
# A. Cell cycle

# grouped histogram
ggplot(df, aes(x=Intensity_IntegratedIntensity_BlueOrig, group=condition, color=condition, fill=condition)) +
  geom_histogram(alpha=0.5, position = 'identity', bins=100)

G1 <- which(df$Intensity_IntegratedIntensity_BlueOrig < 150) 
S_early <- which(df$Intensity_IntegratedIntensity_BlueOrig >=150 &df$Intensity_IntegratedIntensity_BlueOrig <=200 )
S_late <- which(df$Intensity_IntegratedIntensity_BlueOrig >200 &df$Intensity_IntegratedIntensity_BlueOrig <=250 )
G2 <- which(df$Intensity_IntegratedIntensity_BlueOrig >250)

df$cell_cycle <- "null"
df$cell_cycle[G1] <- "G1"
df$cell_cycle[S_early] <- "S_early"
df$cell_cycle[S_late] <- "S_late"
df$cell_cycle[G2] <- "G2"


######## Coloc PLOTS ########

df$nuclear_Intensity_MeanIntensity_GreenOrig
ggplot(df, aes(x=Cell_ID_unique, y=nuclear_Intensity_MeanIntensity_GreenOrig)) +
  geom_point() +
  geom_hline(yintercept=0.002)

nuclear <- df[which(df$nuclear_Intensity_MeanIntensity_GreenOrig > 0.002), "Cell_ID_unique"] #ARBITRARY
df$nuclearGFP <- "not nuclear"
df$nuclearGFP[which(df$Cell_ID_unique %in% nuclear)] <- "nuclear GFP"


ggplot(df[which(df$condition=="PLIN5_Oleate"), ], aes(y=M2coef, x=M1coef, color=transfected)) +
  geom_point(size=1.8) +
  theme_grey() +
  xlab("M1 ~ fraction of Mitochondrial signal overlapping LD") +
  ylab("M2 ~ Fraction of LD signal overlapping Mito")


ggplot(df[which(df$condition=="PLIN5_Oleate"), ], aes(y=M2coef, x=M1coef, color=transfected)) +
  geom_point(size=1.8) +
  theme_grey() +
  xlab("M1 ~ fraction of Mitochondrial signal overlapping LD") +
  ylab("M2 ~ Fraction of LD signal overlapping Mito")




ggplot(df[which(df$condition=="PLIN5_minusOleate"), ], aes(y=M2coef, x=M1coef, color=transfected)) +
  geom_point(size=1.8) +
  theme_grey() +
  xlab("M1 ~ fraction of Mitochondrial signal overlapping LD") +
  ylab("M2 ~ Fraction of LD signal overlapping Mito")

ggplot(df, aes(y=M2coef, x=M1coef, color=condition, shape=transfected)) +
  geom_point(size=1.8) +
  theme_grey() +
  xlab("M1 coefficient") +
  ylab("M2 coefficient")

ggplot(df[which(df$condition=="PLIN5_Oleate"), ], aes(y=M1coef, x=log(Cyto_Intensity_IntegratedIntensity_GreenOrig), shape=nuclearGFP)) +
  geom_point(size=1.8, color="#00C0AF") +
  theme_grey() +
  ylab("M1 ~ Proportion of mitochondrial signal overlapping with LD") +
  xlab("log(Integrated cytosolic GFP intensity)") +
  ggtitle("Mito-LD colocalisation in Oleate loaded cos-7 cells")

ggplot(df[which(df$condition=="PLIN5_minusOleate"), ], aes(y=M1coef, x=log(Cyto_Intensity_IntegratedIntensity_GreenOrig), shape=nuclearGFP)) +
  geom_point(size=1.8, color="#FC717F") +
  theme_grey() +
  ylab("M1 ~ Proportion of mitochondrial signal overlapping with LD") +
  xlab("log(Integrated cytosolic GFP intensity)") +
  ggtitle("Mito-LD colocalisation in cos-7 cells without Oleate")





ggplot(df, aes(y=M1coef, x=log(Cyto_Intensity_IntegratedIntensity_GreenOrig), color=condition, shape=nuclearGFP)) +
  geom_point(size=1.8) +
  theme_grey() +
  ylab("M1 coefficient") 

df$nuclear <- rep("not nuclear", nrow(df))
df$nuclear[which(df$Cell_ID_unique %in% nuclearPLIN5)] <- "nuclear"




#############
library(ggpubr)
compare_means(data=df, Summed_Total_voxels_perCell ~ condition2, method="wilcox.test")

my_comparisons2 <- list(c("PLIN5_minusOleate untransfected", "PLIN5_minusOleate transfected"), c("PLIN5_Oleate untransfected", "PLIN5_Oleate transfected"), c("PLIN5_Oleate transfected", "PLIN5_minusOleate transfected"),  c("PLIN5_minusOleate untransfected", "PLIN5_Oleate untransfected"))

# 1a  total amount of mitochondrial voxels
ggplot(df, aes(x=condition2, y = Summed_Total_voxels_perCell, color=transfected, fill=transfected)) +
  geom_boxplot(alpha = .5) +
  labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.7, hjust=0.8)) +
  stat_compare_means(comparisons = my_comparisons2, hide.ns = FALSE,
                     label ="p.signif", label.x = 1.5, label.y = c(2700, 2600, 3000, 3400), size=2, method = "wilcox.test")
  ggtitle("Total mitochondrial voxels per cell between conditions") 
ggsave(filename = paste0(outfolder, "/totalvoxels_perCell_boxplots.png"))

ggplot(df, aes(x=condition2, y = Intensity_IntegratedIntensity_FarrOrig, color=transfected, fill=transfected)) +
  geom_boxplot(alpha = .5) +
  labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.7, hjust=0.8)) +
  stat_compare_means(comparisons = my_comparisons2, hide.ns = FALSE,
                     label ="p.signif", label.x = 1.5, label.y = c(350, 370, 450, 500), size=2, method = "wilcox.test")
ggtitle("Total mitochondrial intensity per cell between conditions") 
ggsave(filename = paste0(outfolder, "/totalvoxels_perCell_boxplots.png"))


# 1b Mito classification size based


ggplot(MitoDf, aes(fill=Mito_classification_size, x=condition2))+ 
  geom_bar( stat="count", position="fill") +
  labs(x = "") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.7, hjust=0.8))
ggsave(filename = paste0(outfolder, "/Mito_classification_distribution_size.png"))


# 1c Mito classification
ggplot(MitoDf, aes(fill=Mito_classification, x=condition2))+ 
  geom_bar( stat="count", position="fill" ) +
  labs(x = "") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.7, hjust=0.8))
ggsave(filename = paste0(outfolder, "/Mito_classification_distribution_conditions.png"))

ggplot(MitoDf, aes(fill=Mito_classification_simple, x=condition))+ 
  geom_bar( stat="count", position="fill" ) +
  labs(x = "") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.7, hjust=0.8))+
  ggtitle("Classification of mitochondria in untransfected cells")
ggsave(filename = paste0(outfolder, "/Mito_classification_distribution_conditions.png"))


ggplot(MitoDf, aes(y=number_of_mitochondria_perCell, x=condition2))+ 
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.7, hjust=0.8))





# 2. does oleate loading change nr of mitos in plin5 cells?

ggplot(MitoDf, aes(x=condition, y = number_of_mitochondria_perCell, color = transfected, fill=transfected)) +
  geom_boxplot(alpha = .5) +
  labs(x = "") +
  stat_compare_means(comparisons = my_comparisons, hide.ns = FALSE,
                     label ="p.signif", label.x = 1.5, label.y = c(100), size=2, method = "wilcox.test")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.7, hjust=0.8)) 
ggsave(filename = paste0(outfolder, "/Branchlength_boxplots.png"))


  
# percentage of voxels is better than absolute numbers because the larger a mito the less it will contribute to absolute nr
# fragmented
ggplot(df, aes(x=condition2 , y =percentage_fragmented_voxels_perCell, fill = condition, color=condition)) +
  geom_boxplot(alpha = .5) +
  labs(x = "") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.7, hjust=0.8)) +
  stat_compare_means(comparisons = my_comparisons2, hide.ns = FALSE,
                     label ="p.signif", label.x = 1.5, label.y = c(90, 80, 70), size=2, method = "wilcox.test")
ggsave(filename = paste0(outfolder, "/Mitocontent_fragmented_boxplots.png"))

# tubular
ggplot(df, aes(x=condition2, y =percentage_tubular_voxels_perCell, color = condition, fill=condition)) +
  geom_boxplot(alpha = .5) +
  labs(x = "") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.7, hjust=0.8)) +
  stat_compare_means(comparisons = my_comparisons2, hide.ns = TRUE,
                     label ="p.signif", label.x = 1.5, label.y = c(90, 80, 70), size=2, method = "wilcox.test")
ggsave(filename = paste0(outfolder, "/Mitocontent_tubular_boxplots.png"))

# intermediate
ggplot(df, aes(x=condition2, y =percentage_intermediate_voxels_perCell, color = condition, fill=condition)) +
  geom_boxplot(alpha = .5) +
  labs(x = "") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.7, hjust=0.8)) +
  stat_compare_means(comparisons = my_comparisons, hide.ns = TRUE,
                     label ="p.signif", label.x = 1.5, label.y = c(90, 80, 70), size=2, method = "wilcox.test")
ggsave(filename = paste0(outfolder, "/Mitocontent_tubular_boxplots.png"))


# Explore if the mitochondria are different 
ggplot(MitoDf, aes(x=nr_Total_voxels, y=adapted_ratio_Junctions_Endpoints, color=condition2, shape=transfected)) +
  geom_point(size=1) +
  theme_grey()

MitoDf$variance_mitochondrial_size_perCell
MitoDf$adapted_ratio_Junctions_Endpoints
MitoDf$Intensity_MeanIntensity_GreenOrig
# Explore if the mitochondria are different   ############
ggplot(MitoDf, aes(x=log(Intensity_MeanIntensity_GreenOrig), y=adapted_ratio_Junctions_Endpoints, color=condition, shape=transfected)) +
  geom_point(size=1) +
  theme_grey()

# absolute summed voxels, percentage above is better i think since the total voxels per cell are different in each group
ggplot(MitoDf, aes(x=condition, y = summed_voxels_tubular_perCell, color = condition, fill=condition, )) +
  geom_boxplot(alpha = .5) +
  labs(x = "") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.7, hjust=0.8)) +
  stat_compare_means(comparisons = my_comparisons, hide.ns = TRUE,
                     label ="p.signif", label.x = 1.5, label.y = c(90, 80, 70), size=2, method = "wilcox.test")
ggsave(filename = paste0(outfolder, "/Mitocontent_tubular_boxplots.png"))


ggplot(MitoDf, aes(x=condition, y = summed_voxels_intermediate_perCell, color = condition, fill=transfected)) +
  geom_boxplot(alpha = .5) +
  labs(x = "") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.7, hjust=0.8)) +
  stat_compare_means(comparisons = my_comparisons, hide.ns = TRUE,
                     label ="p.signif", label.x = 1.5, label.y = c(90, 80, 70), size=2, method = "wilcox.test")
ggsave(filename = paste0(outfolder, "/Mitocontent_tubular_boxplots.png"))


ggplot(MitoDf, aes(x=condition, y = summed_voxels_fragmented_perCell, color = condition, fill=condition)) +
  geom_boxplot(alpha = .5) +
  labs(x = "") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.7, hjust=0.8)) +
  stat_compare_means(comparisons = my_comparisons, hide.ns = TRUE,
                     label ="p.signif", label.x = 1.5, label.y = c(90, 80, 70), size=2, method = "wilcox.test")
ggsave(filename = paste0(outfolder, "/Mitocontent_tubular_boxplots.png"))


ggplot(df, aes(x=condition, y = largest_mitochondrion_perCell, color = condition, fill=condition)) +
  geom_boxplot(alpha = .5) +
  labs(x = "") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.7, hjust=0.8)) +
  ggtitle("largest mitochondrion per cell")+
  stat_compare_means(comparisons = my_comparisons, hide.ns = TRUE,
                     label ="p.signif", label.x = 1.5, label.y = c(90, 80, 70), size=2, method = "wilcox.test")
ggsave(filename = paste0(outfolder, "/largest_boxplots.png"))


# number of mitos in each size category per cell
ggplot(MitoDf[which(MitoDf$transfected=="transfected"), ], aes(fill=Mito_classification_simple, x=Cell_ID_unique) ) + 
         geom_bar( stat="count") +
         facet_wrap(vars(condition), scales = "free_x") +
         theme(axis.text.x=element_blank(),
               axis.ticks.x=element_blank()) +
         labs(x="Cells") 
ggsave(filename = paste0(outfolder, "/Mito_classification_size_distribution_conditions_singlecell.png") , dpi=300)
    

ggplot(MitoDf[which(MitoDf$transfected=="untransfected"), ], aes(fill=Mito_classification_simple, x=Cell_ID_unique) ) + 
  geom_bar( stat="count") +
  facet_wrap(vars(condition), scales = "free_x") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  labs(x="Cells") 
ggsave(filename = paste0(outfolder, "/Mito_classification_size_distribution_conditions_singlecell.png") , dpi=300)


ggplot(MitoDf, aes(fill=Mito_classification, x=Cell_ID_unique) ) + 
  geom_bar( stat="count") +
  facet_wrap(vars(condition), scales = "free_x") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  labs(x="Cells") 
ggsave(filename = paste0(outfolder, "/Mito_classification_size_distribution_conditions_singlecell.png") , dpi=300)




### Exploratory data analysis ###

# Plot the integrated intensity distribution of the DAPI signal within nuclei to check cell cycle (Integrated DAPI ~ DNA content)
hist(df$nuclear_Intensity_IntegratedIntensity_BlueOrig, breaks=100) # it seems as if most cells are in the G1 and S phase
plot(df$nuclear_Intensity_MeanIntensity_GreenOrig, df$nuclear_Intensity_IntegratedIntensity_BlueOrig)  # Nuclear GFP is not correlated with cell cycle

ggplot(df, aes(x=nuclear_Intensity_IntegratedIntensity_BlueOrig, y=nuclear_Intensity_MeanIntensity_GreenOrig, color=transfected))+
  geom_jitter(width = .05, alpha = .8) +
  xlab("nuclear Integrated intensity of DNA") +
  ylab("nuclear Integrated intensity GFP")
ggsave(filename=paste0(outfolder,"NuclearGFP_vs_CellCycle.png"), dpi=300 )

# See the mean intensity of the cells
plot(df$Intensity_MeanIntensity_GreenOrig)  # Can be used to adjust the threshold of classifying cells, see thresholdTransfected above

# Which cells have green in the nucleus
plot(CP_nuc$Intensity_IntegratedIntensity_GreenOrig)
nuclearGFP <- CP_nuc[which(CP_nuc$Intensity_IntegratedIntensity_GreenOrig > 200), "Cell_ID_unique"] #ARBITRARY

df$nuclear <- rep("not nuclear", nrow(df))
df$nuclear[which(df$Cell_ID_unique %in% nuclearGFP)] <- "nuclear"

# Check if I had different acquisition settings which could explain integrated intensity difference but no.
ggplot(df, aes(y=Intensity_MeanIntensity_GreenOrig, x=Cell_ID_unique, color=condition)) +
  geom_point(size=0.8) +
  #theme(legend.position="none") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
ggsave(filename=paste0(outfolder, "/meanGFP_intensity_condition.png"), dpi=300)



# 





### MDF PLOTS ###
datafolderNHS = "/Volumes/ExtDriveLou/LAB21/LAB21clean/NHS/"


# fix conditions

MDf <- read.csv(paste0(datafolderNHS,"CP_measurementsAllCytoplasm.csv"))
tf <-read.csv(paste0(datafolderNHS, "CP_measurements2AllCells.csv"))
tf$Cell_ID_unique <- add_Cell_ID_unique(tf)
tf$transfected <- "untransfected"
tf$transfected[which(tf$Children_TransfectedCells_Count > 0)] <- "transfected"


tf$condition <- "null"
tf$condition[grep("PLIN5-Oleate", tf$FileName_BlueOrig, value=FALSE)] <- "PLIN5_Oleate"
tf$condition[grep("PLIN5-minusOleate", tf$FileName_BlueOrig, value=FALSE)] <- "PLIN5_minusOleate"

#old in cyto
MDf$Cell_ID_unique <- add_Cell_ID_unique(MDf)
MDf$transfected <- "transfected"
MDf$transfected[which(tf$Children_TransfectedCells_Count > 0)] <- "untransfected"
MDf$condition <- "null"
MDf$condition[grep("PLIN5-Oleate", MDf$FileName_BlueOrig, value=FALSE)] <- "PLIN5_Oleate"
MDf$condition[grep("PLIN5-minusOleate", MDf$FileName_BlueOrig, value=FALSE)] <- "PLIN5_minusOleate"
##



tf$Cell_ID_unique <- add_Cell_ID_unique(tf)
tf$transfected <- "transfected"
tf$transfected[which(tf$Children_TransfectedCells_Count > 0)] <- "untransfected"

#FracDf <- MDf[ , grepl("RadialDistribution_FracAtD" ,colnames(MDf)) & grepl("Red", colnames(MDf))]
#FracDf_green <- MDf[ , grepl("RadialDistribution_FracAtD" ,colnames(MDf)) & grepl("GreenOrig", colnames(MDf))]

FracDf <- tf[ , grepl("RadialDistribution_FracAtD" ,colnames(tf)) & grepl("Red", colnames(tf))]
FracDf_green <- tf[ , grepl("RadialDistribution_FracAtD" ,colnames(tf)) & grepl("GreenOrig", colnames(tf))]

RDFs <- FracDf
RDFs_g <-FracDf_green

tf$cell_cycle <-"null"

tf$Cell_ID_unique

for (i in 1:nrow(tf)){
  tf$MDF_mito[i] <- as.numeric(extract_MDFs(RDFs[i,]))
  tf$MDF_GFP[i] <- as.numeric(extract_MDFs(RDFs_g[i,]))
  tf$cell_cycle[i] <- df$cell_cycle[which(df$Cell_ID_unique == tf$Cell_ID_unique[i])]
}



ggplot(tf, aes(x=Intensity_IntegratedIntensity_BlueOrig, group=condition, color=condition, fill=condition)) +
  geom_histogram(alpha=0.5, position = 'identity', bins=100)

G1 <- which(tf$Intensity_IntegratedIntensity_BlueOrig < 150) 
S_early <- which(tf$Intensity_IntegratedIntensity_BlueOrig >=150 &tf$Intensity_IntegratedIntensity_BlueOrig <=200 )
S_late <- which(tf$Intensity_IntegratedIntensity_BlueOrig >200 &tf$Intensity_IntegratedIntensity_BlueOrig <=250 )
G2 <- which(tf$Intensity_IntegratedIntensity_BlueOrig >250)

tf$cell_cycle <- "null"
tf$cell_cycle[G1] <- "G1"
tf$cell_cycle[S_early] <- "S_early"
tf$cell_cycle[S_late] <- "S_late"
tf$cell_cycle[G2] <- "G2"


df$cell_cycle

MDf$Intensity_IntegratedIntensity_GreenOrig
tf$Intensity_IntegratedIntensity_GreenOrig

ggplot(tf, aes(x=log(Intensity_IntegratedIntensity_GreenOrig), y=MDF_mito, color=condition, shape=transfected)) +
  geom_point()
ggsave(filename=paste0(outfolder, "/GFP_intensity_MDF.png"), dpi=300)



ggplot(df, aes(x=MDF_mitochondria, y=Cyto_Intensity_IntegratedIntensity_GreenOrig, color=condition, shape=transfected)) +
  geom_point() +
  labs(y = "mean GFP intensity cytosol",
       x = "MDF mitochondria") + 
  ggtitle("GFP intensity and mitochondrial centricity")
ggsave(filename=paste0(outfolder, "/GFP_intensity_MDF.png"), dpi=300)

tf$Cell_ID_unique
df$md
for (i in 1:nrow(df)){
  df$MDF_mitochondria[i] <- tf$MDF_mito[ which(tf$Cell_ID_unique == df$Cell_ID_unique[i]) ]
  df$MDF_GFP[i] <- tf$MDF_GFP[ which(tf$Cell_ID_unique == df$Cell_ID_unique[i]) ]
}



