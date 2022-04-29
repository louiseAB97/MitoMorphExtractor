library(Hmisc)
library(dplyr)
library(tidyverse) # ggplot2 etc
library(bayestestR) # to calculate AUC

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
  MitoDf <- read.csv( paste0(datafolder,"AutoCBtiffs/MitoMorphResults_sigma03_new.txt")) 
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
  transf <- df[which(df$Children_TransfectedCells_Count > 0), "Cell_ID_unique"] # based on the CP classification
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
  super_fragmented_rownrs <- which(MitoDf$nr_Total_voxels < 4)
  fragmented_rownrs <- which(MitoDf$nr_Total_voxels > 4 & MitoDf$nr_Total_voxels < 15)
  tubular_rownrs <-  which(MitoDf$nr_Total_voxels >= 100)
  fused_rownrs <- which(MitoDf$adapted_ratio_Junctions_Endpoints > 0 & MitoDf$adapted_ratio_Junctions_Endpoints <= 0.5)
  branched_rownrs <- which(MitoDf$adapted_ratio_Junctions_Endpoints > 0.5)
  ringshaped_rownrs <- which( MitoDf$adapted_ratio_Junctions_Endpoints < 0)
  classic_tubular_rownrs <- which(MitoDf$adapted_ratio_Junctions_Endpoints == 0 & MitoDf$nr_Junctions==0)
  
  MitoDf$Mito_classification_simple <- "Intermediate"
  MitoDf$Mito_classification_simple[tubular_rownrs] <- "Tubular"
  MitoDf$Mito_classification_simple[fragmented_rownrs] <- "Fragmented"
  MitoDf$Mito_classification_simple[super_fragmented_rownrs] <- "SuperFragmented / noise"
  
  MitoDf$Mito_classification <- "other"
  MitoDf$Mito_classification[ringshaped_rownrs] <- "Ringshaped"
  MitoDf$Mito_classification[fused_rownrs] <- "Fused"
  MitoDf$Mito_classification[branched_rownrs] <- "Branched"
  MitoDf$Mito_classification[classic_tubular_rownrs] <- "Classic tube"
  
  MitoDf$Mito_classification_size <- "1-5 voxels"
  MitoDf$Mito_classification_size[which(MitoDf$nr_Total_voxels >= 5 & MitoDf$nr_Total_voxels <=10) ] <- "5-10 voxels"
  MitoDf$Mito_classification_size[which(MitoDf$nr_Total_voxels >10 & MitoDf$nr_Total_voxels <=30) ] <- "10-30 voxels"
  MitoDf$Mito_classification_size[which(MitoDf$nr_Total_voxels >30 & MitoDf$nr_Total_voxels <=50)] <- "30-50 voxels "
  MitoDf$Mito_classification_size[which(MitoDf$nr_Total_voxels >50 & MitoDf$nr_Total_voxels <=100)] <- "50-100 voxels"
  MitoDf$Mito_classification_size[which(MitoDf$nr_Total_voxels >100)] <- ">100 voxels"
  
  return(list(df, MitoDf))
}


# LOAD DATA this example shows how to load (and combine different) outputs
datafolder1 = "~/LocalFiles/LAB18(INFECTED)/2hforskolin_infected/"
datafolder2 = "/Users/louise/LocalFiles/LAB18(INFECTED)/PLIN5_infected/"
datafolder3 = "/Users/louise/LocalFiles/LAB17/combined/"

outfolder <- "/Users/louise/LocalFiles/outputLAB1718/"
dir.create(outfolder)

condition1 = "PLIN5_Oleate_2h_forskolin_infected"
condition2 = "PLIN5_Oleate_infected"

data1 <- get_combined_dataframe(datafolder1, condition1, "Red")
data2 <- get_combined_dataframe(datafolder2, condition2, "Red")

df1 <- data1[[1]]
df2 <- data2[[1]]
MitoDf1 <- data1[[2]]
MitoDf2 <- data2[[2]]

df <- rbind(df1, df2)
MitoDf <- rbind(MitoDf1, MitoDf2)

# To load data with different conditions in one ouput file:
df3 <- data3[[1]]
MitoDf3 <- data3[[2]]


PLIN5 <- grep("GFP-PLIN5", df3$FileName_BlueOrig, value=FALSE) #grep the pattern that you named your condition
PLIN5_minOleate <- grep("minusO_DAPI_GFP-PLIN5", df3$FileName_BlueOrig , value=FALSE)
PLIN5_Oleate <- setdiff(PLIN5,PLIN5_minOleate)
PLIN3_Oleate <- grep("GFP-PLIN3", df3$FileName_BlueOrig, value=FALSE)
df3$condition[PLIN3_Oleate] <- "PLIN3_Oleate"
df3$condition[PLIN5_Oleate] <- "PLIN5_Oleate"
df3$condition[PLIN5_minOleate] <- "PLIN5_minOleate"

for (i in 1:nrow(MitoDf3)){
  MitoDf3$condition[i] <- df3$condition[ which(df3$Cell_ID_unique == MitoDf3$Cell_ID_unique[i]) ]
}

df <- rbind(df, df3)
MitoDf <- rbind(MitoDf, MitoDf3)
