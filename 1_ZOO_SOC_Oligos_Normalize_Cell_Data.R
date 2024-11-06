#!/usr/bin/env Rscript --vanilla

#-----------------------------
# Author
#-----------------------------

# The script was written by Tjard Bergmann (2023) and is not licensed.
# Everyone is free to use or edit it under the rules of creative commons (CC BY).

# University of Veterinary Medicine Hannover
# Institute of Zoology, AG Felmy (Neurobiology)
# Website: https://www.tiho-hannover.de/kliniken-institute/institute/institut-fuer-zoologie/arbeitsgruppen/neurobiologie/ag-felmy

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

#-----------------------------
# Reference
#-----------------------------

#Please cite the following manuscript when using this script:

#Developmental profile of oligodendrocyte arrangement, identification and morphology 
#in nuclei of the superior olivary complex  

#-----------------------------
# Description
#-----------------------------

# This R-script normalizes the cell positioning data of 
# animal brain slices focusing on the SOC.

#-----------------------------
# Software dependencies
#-----------------------------

#The script must be executed within the Editor: RStudio
#Source: https://posit.co/

#This script is dependent on the following functions:

#cowplot     
#dplyr        
#EnvStats     
#ggplot2     
#RColorBrewer
#rstudioapi  
#stringr      
#writexl

#-----------------------------
# Guide
#-----------------------------

# To successfully run the program you need to:

# 1) The program is dependent on two input files (Combined_ImageJ-data.csv, ROI-data.csv).
#    The files are zipped in "input_data.zip" on github (https://github.com/M0rph3u2x/Developmental_profile_of_oligodendrocyte_in_nuclei_of_the_superior_olivary_complex/upload/main).
#    Both of these files where created with ImageJ (https://imagej.net/).
#    Combined_ImageJ-data.csv -> Holds all cell coordination data
#    ROI-data.csv             -> Holds the ROI coordinate data for all SOC nuclei
#    Both files must be placed in the same folder as this script!

# 2) Press Ctr+A, then Ctrl+Enter to execute the program

# 3) The script will create a folder with control plot data of each nuclei.
#    Each point in the plot represent a single cell, while the black frame (ROI) represent
#    the outer edge of each nuclei. Cell coordinate data outside of the ROI are discarded.
#    Correct cell data is further processed, normalized and saved in the file (1_Density_Plot_Normalized_data_pixel.xlsx).
#    Discarded cell data is saved in (1_Density_Plot_Discarded_data_pixel.xlsx)

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Load necessary libraries -----------------------------------------------------
#-------------------------------------------------------------------------------

# Define necessary package manager
packages <- c("pacman")

# Install package manager if not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

#Install and/or load all necessary packages via package manager "pacman"
pacman::p_load(cowplot,      # function (plot_grid) creates multiplot figures
               dplyr,        # function(%>%) -> Getting subgroup in dataframe
               EnvStats,     # Sample Size in ggplot2 (function stat_n_text)
               ggplot2,      # Plotting the MNTB data
               RColorBrewer, # Color palette
               rstudioapi,   # Get path from RStudio (function: getSourceEditorContext()$path)
               stringr,      # Extract string (function str_sub)
               writexl)      # function(write_xlsx) [more efficient then write.xlsx!]

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Setup directory where the program is running ---------------------------------
#-------------------------------------------------------------------------------

#Create path to data
library(rstudioapi)    #Get path from RStudio (function: getSourceEditorContext()$path)
R.Directory = sub(pattern = "(.*/).*\\..*$", replacement = "\\1", getSourceEditorContext()$path)

#Create shortcut to main path (Fullpath connects R.directory with follow up folders)
FullPath = function(FileName){ return( paste(R.Directory,FileName,sep="") ) }

# Set the working directory to the folder where the program is running
setwd(R.Directory)

#Submit working directory to variable
work_path <- getwd()

#Create main output folder
out.path <- FullPath("control_plots")
dir.create(out.path, showWarnings = FALSE)

#-------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------
# Calculate cell distribution in Oligodendroytes (color coding cell position density x/y-axis) ---------
#--------------------------------------------------------------------------------------------

#Load csv summary table (contains cummulative data of all cell positions) ------
csv.table.name <- "Combined_ImageJ-data.csv"
#csv.path       <- file.path(out.path, csv.table.name)
cell_data      <- read.csv(file=csv.table.name)

#Load roi summary table (contains cummulative data of all nuclei shape data) ---
roi.table.name <- "ROI-data.csv"
#roi.path       <- file.path(out.path, roi.table.name)
roi_data       <- read.csv(file=roi.table.name)

#Convert pixel to µm [0.586 = Umrechnungsfaktor]
# roi_data[, c(8:13,15,18:19)] <- round(roi_data[, c(8:13,15,18:19)] * 0.586, digits = 2)
# roi_data[,               14] <- round(roi_data[,       14] * 0.586 * 0.586, digits = 2) #Area µm²
# roi_data[,               14] <- roi_data[,14] /1000000 #Area mm²

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Extract and normalize cell position data of the Oligodendrocyte regions:
# "Kontrolle"  "LSO_VGlut-" "LSO_VGlut+" 
#-------------------------------------------------------------------------------

#Extract columns important for distance plot
dist.table <- cell_data[c("Slice_Name","Age","Nucleus","Region","X","Y")]
roi.table  <- roi_data [c("Slice","Age","Nucleus","Region","BX","BY","Width","Height")]

#Create new table that will get new data
DistancePlot.Table <- data.frame()
Discarded_Data     <- data.frame()

#Transform x & y axis data into standardized position data
Region_NR <- 4 #("Kontrolle" ,  "LSO_VGlut-",  "LSO_VGlut+",  "MNTB")

#Get standardized cell position data (X/Y_Norm) and normalized location data (X/Y_Perc)
#for all single region nuclei
for(reg_nr in 1:Region_NR){#reg_nr <-1
  
  #Get cell position data for selected region
  Reg_Name <- unique(dist.table$Region)[reg_nr]
  Reg_Data <- subset(dist.table, Region==Reg_Name)
  
  #Get ROI-Data
  Roi_Data_Reg <- subset(roi.table, Region==Reg_Name)
  
  #Normalize data and calculate percentage values for each slice
  for(slice in unique(Reg_Data$Slice)){#slice <-"P05_15.03.21_T1_OT1_Schnitt6_l"
    
    #Extract slice cell data
    Sli_Cell <- subset(Reg_Data, Slice_Name==slice)
    
    #Extract slice roi data
    Sli_Roi <- subset(Roi_Data_Reg, Slice==slice)
    
    #Calculate ROI min/max values for x/y
    Roi_Min_X <- Sli_Roi$BX
    Roi_Max_X <- Sli_Roi$BX+Sli_Roi$Width
    Roi_Min_Y <- Sli_Roi$BY-Sli_Roi$Height
    Roi_Max_Y <- Sli_Roi$BY
    
    Orientation <- substr(slice,nchar(slice),nchar(slice))
    if(Orientation=="l"){
      X_Norm   <- Roi_Max_X-Sli_Cell$X
    }else{
      X_Norm   <- Sli_Cell$X-Roi_Min_X
    }
    Y_Norm   <-  Sli_Cell$Y-Roi_Min_Y
    
    Sli_Cell$X_Norm <- X_Norm
    Sli_Cell$Y_Norm <- Y_Norm
    Sli_Cell$X_Perc <- round((100/Sli_Roi$Width) *Sli_Cell$X_Norm)
    Sli_Cell$Y_Perc <- round((100/Sli_Roi$Height)*Sli_Cell$Y_Norm)
    
    ############################################################################
    #Control cell position inside/outside ROI-Rectange (Outside=False Data) ----
    
    #Create control plot
    rect_coord  <- data.frame(x1=c(Roi_Min_X,Roi_Min_X,Roi_Max_X,Roi_Max_X,Roi_Min_X),
                              y1=c(Roi_Min_Y,Roi_Max_Y,Roi_Max_Y,Roi_Min_Y,Roi_Min_Y))
    slice_coord <- data.frame(x2=Sli_Cell$X,y2=Sli_Cell$Y)
    
    # Fügen Sie eine Spalte hinzu, um anzugeben, ob der Punkt innerhalb des Rechtecks liegt
    slice_coord$inside <- with(slice_coord, x2 >= Roi_Min_X & x2 <= Roi_Max_X & y2 >= Roi_Min_Y & y2 <= Roi_Max_Y)
    
    # Fügen Sie eine Spalte hinzu, um anzugeben, ob der Punkt innerhalb des Rechtecks liegt
    slice_coord$outside <- with(slice_coord, x2 < Roi_Min_X | x2 > Roi_Max_X | y2 < Roi_Min_Y | y2 > Roi_Max_Y)
    
    # Erstellen Sie den Plot
    label <- paste(Reg_Name,slice,sep="_")
    ggplot() +
      geom_polygon(data = rect_coord, aes(x1, y1), fill = NA, colour = "black") +
      geom_point(data = slice_coord, aes(x2, y2, colour = inside)) +
      scale_color_manual(values = c("TRUE" = "green", "FALSE" = "red")) +
      labs(x = "X-Achse", y = "Y-Achse", title = label)
    if(length(Sli_Cell[slice_coord$inside,1])==length(slice_coord$x2)){
      label <- paste(slice,Reg_Name,"controlplot.tiff",sep="_")
      label <- paste("control_plots",label,sep="/")
    }else{
      label <- paste("Error",slice,Reg_Name,"controlplot.tiff",sep="_")
      label <- paste("control_plots",label,sep="/")
    }
    ggsave(label, dpi=300, width=6, height=6)
    
    #---------------------------------------------------------------------------
    ############################################################################
    
    if(length(Sli_Cell[slice_coord$inside,1])==length(slice_coord$x2)){
      DistancePlot.Table <- rbind(DistancePlot.Table,Sli_Cell)
    }else{
      DistancePlot.Table <- rbind(DistancePlot.Table,Sli_Cell[slice_coord$inside,])
      Discarded_Data     <- rbind(Discarded_Data    ,Sli_Cell[slice_coord$outside,])
    }
  }
}

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Extract, combine and normalize cell position data of "MSO" regions
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Extract, combine and normalize cell position data of the Oligodendrocyte regions:
# ["MSO_lateral"      "MSO_medial"       "MSO_Mitte"] -> "MSO"
#-------------------------------------------------------------------------------

# #Get Nucleus Data (MSO)
Nuc_Name     <- unique(dist.table$Nucleus)[4]
Nuc_Data_MSO <- subset(dist.table, Nucleus==Nuc_Name)

# #Get ROI Data (MSO)
Roi_Data_MSO <- subset(roi.table, Nucleus==Nuc_Name)

# #Get standardized cell position data (X/Y_Norm) and normalized location data (X/Y_Perc)
# #for ("MSO_lateral"      "MSO_medial"       "MSO_Mitte") multi region nuclei

#Normalize data and calculate percentage values for each slice
for(slice in unique(Nuc_Data_MSO$Slice)){#slice <-"P05_15.03.21_T1_OT1_Schnitt6_l"
  
  #Extract slice cell data
  Sli_Cell <- subset(Nuc_Data_MSO, Slice_Name==slice)
  
  #Extract slice roi data
  Sli_Roi <- subset(Roi_Data_MSO, Slice==slice)
  
  #Calculate ROI min/max values for x/y
  Roi_Min_X   <- min(Sli_Roi$BX)
  Roi_Max_X   <- max(Sli_Roi$BX+Sli_Roi$Width)
  Roi_Min_Y   <- min(Sli_Roi$BY-Sli_Roi$Height)
  Roi_Max_Y   <- max(Sli_Roi$BY)
  Roi_Delta_X <- Roi_Max_X-Roi_Min_X
  Roi_Delta_Y <- Roi_Max_Y-Roi_Min_Y
  
  Orientation <- substr(slice,nchar(slice),nchar(slice))
  if(Orientation=="l"){
    X_Norm   <- Roi_Max_X-Sli_Cell$X
  }else{
    X_Norm   <- Sli_Cell$X-Roi_Min_X
  }
  Y_Norm   <-  Sli_Cell$Y-Roi_Min_Y
  
  Sli_Cell$X_Norm <- X_Norm
  Sli_Cell$Y_Norm <- Y_Norm
  Sli_Cell$X_Perc <- round((100/Roi_Delta_X) *Sli_Cell$X_Norm)
  Sli_Cell$Y_Perc <- round((100/Roi_Delta_Y)*Sli_Cell$Y_Norm)
  
  ##############################################################################
  #Control cell position inside/outside ROI-Rectange (Outside=False Data) ------
  
  #Create control plot
  rect_coord  <- data.frame(x1=c(Roi_Min_X,Roi_Min_X,Roi_Max_X,Roi_Max_X,Roi_Min_X),
                            y1=c(Roi_Min_Y,Roi_Max_Y,Roi_Max_Y,Roi_Min_Y,Roi_Min_Y))
  slice_coord <- data.frame(x2=Sli_Cell$X,y2=Sli_Cell$Y)
  
  # Fügen Sie eine Spalte hinzu, um anzugeben, ob der Punkt innerhalb des Rechtecks liegt
  slice_coord$inside <- with(slice_coord, x2 >= Roi_Min_X & x2 <= Roi_Max_X & y2 >= Roi_Min_Y & y2 <= Roi_Max_Y)
  
  # Fügen Sie eine Spalte hinzu, um anzugeben, ob der Punkt innerhalb des Rechtecks liegt
  slice_coord$outside <- with(slice_coord, x2 < Roi_Min_X | x2 > Roi_Max_X | y2 < Roi_Min_Y | y2 > Roi_Max_Y)
  
  # Erstellen Sie den Plot
  label <- paste(Reg_Name,slice,sep="_")
  ggplot() +
    geom_polygon(data = rect_coord, aes(x1, y1), fill = NA, colour = "black") +
    geom_point(data = slice_coord, aes(x2, y2, colour = inside)) +
    scale_color_manual(values = c("TRUE" = "green", "FALSE" = "red")) +
    labs(x = "X-Achse", y = "Y-Achse", title = label)
  if(length(Sli_Cell[slice_coord$inside,1])==length(slice_coord$x2)){
    label <- paste(slice,"MSO","controlplot.tiff",sep="_")
    label <- paste("control_plots",label,sep="/")
  }else{
    label <- paste("Error",slice,"MSO","controlplot.tiff",sep="_")
    label <- paste("control_plots",label,sep="/")
  }
  ggsave(label, dpi=300, width=6, height=6)
  
  #-----------------------------------------------------------------------------
  ##############################################################################
  
  if(length(Sli_Cell[slice_coord$inside,1])==length(slice_coord$x2)){
    DistancePlot.Table <- rbind(DistancePlot.Table,Sli_Cell)
  }else{
    DistancePlot.Table <- rbind(DistancePlot.Table,Sli_Cell[slice_coord$inside,])
    Discarded_Data     <- rbind(Discarded_Data    ,Sli_Cell[slice_coord$outside,])
  }
}  
 
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Format and save DistancePlot Table for cell distribution plots 
#-------------------------------------------------------------------------------

#Bundle data
#DistancePlot.Table_BU <- DistancePlot.Table
#unique(as.numeric(str_sub(DistancePlot.Table$Age, start=2,end=3)))
DistancePlot.Table$Age <- as.numeric(str_sub(DistancePlot.Table$Age, start=2,end=3))
DistancePlot.Table$Age[DistancePlot.Table$Age>=33 & DistancePlot.Table$Age<39] <- 36 
DistancePlot.Table$Age[DistancePlot.Table$Age>49] <- 54
DistancePlot.Table[DistancePlot.Table$Region == "MSO_lateral" |
                   DistancePlot.Table$Region == "MSO_medial" |  
                   DistancePlot.Table$Region == "MSO_mitte","Region"] <- "MSO"


#Extract cell data beyond ROI-Rectangular

# Discarded_Data <- subset(DistancePlot.Table, X_Perc > 100 | X_Perc < 0 | Y_Perc > 100 | Y_Perc < 0)
# sort(unique(DistancePlot.Table$X_Perc))
# sort(unique(DistancePlot.Table$Y_Perc))

#Create csv summary table
#write.table(DistancePlot.Table, file="1_Density_Plot_Normalized_data_pixel.csv", sep=",", row.names=FALSE, col.names=TRUE,  na="NA")

#Create sorted xlsx table of ROI density data

### Function ###

#Create sheets in xlsx file (tuned for cell dataset)
create_xlsx_sheet <- function(sheet.data){
  
  #Create two empty columns with the length of data
  empty.col <- rep("",length(sheet.data[,1]))
  empty.col <- data.frame("Empty1"=empty.col,"Empty2"=empty.col)
  
  #Combine empty columns with data columns
  sheet.data <- cbind(empty.col,sheet.data)
  sheet      <- list(sheet = sheet.data)
  return(sheet)
}

################

#Sort by slice
DistancePlot.Table <- DistancePlot.Table[order(DistancePlot.Table$Slice_Name,DistancePlot.Table$Age),]

#Create sheets
sheets <- list("All_data" = DistancePlot.Table)

#Define sheets
region <- c("MSO","MNTB","LSO_VGlut+","LSO_VGlut-","Kontrolle")

for(reg in region){
  sheet.data <- DistancePlot.Table[DistancePlot.Table$Region==reg,]
  sheet      <- create_xlsx_sheet(sheet.data)
  sheets     <- c(sheets,sheet)
}

#Attach sheet names to xlsx-sheets
names(sheets) <- c("All_data","MSO","MNTB","LSO_VGlut+","LSO_VGlut-","Kontrolle")

#Create xlsx summary table (contains cummulative data of all processed MNTBs)
write_xlsx(sheets, "1_Density_Plot_Normalized_data_pixel.xlsx")

#Create table of discarded data
write_xlsx(Discarded_Data, "1_Density_Plot_Discarded_data_pixel.xlsx")

#-------------------------------------------------------------------------------