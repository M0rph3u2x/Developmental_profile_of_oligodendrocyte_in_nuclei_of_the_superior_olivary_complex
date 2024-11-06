################################################################################
#Note: First execute the R script "1_ZOO_SOC_Oligos_Normalize_Cell_Data.R"
#      This will load all necessary libraries and create the dataframe "DistancePlot.Table"
#      which is needed in order to run this R script!
################################################################################

#-------------------------------------------------------------------------------

#-----------------------------
# Author
#-----------------------------

# The script was written by Tjard Bergmann (2023) and is not licensed.
# Everyone is free to use or edit it under the rules of creative commons (CC BY).

# University of Veterinary Medicine Hannover
# Institute of Zoology, AG Felmy (Neurobiology)
# Website: https://www.tiho-hannover.de/kliniken-institute/institute/institut-fuer-zoologie/arbeitsgruppen/neurobiologie/ag-felmy

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

# 1) Execute the R script "1_ZOO_SOC_Oligos_Normalize_Cell_Data.R" directly before running this script
#    This script is dependent on the dataframe "DistancePlot.Table" which is created when running the
#    first script.

# 2) Press Ctr+A, then Ctrl+Enter to execute this script (it must be located at the same location as the
#    first script).

# 3) The script will create the folder "output_average" which will contain different subfolders
#    containing different steps of processing the data.
#    average_plot: Contains the average plot of all processed nuclea and selected age groups
#    averaged_data: Contains the combined and averaged raw data of all processed nuclea data
#    binned_data: Contains the binned data of the averaged raw data
#    xlsx_grid: Contains 2D grid data of binned averaged ndensity of cells for each nuclea and selected age groups

#-------------------------------------------------------------------------------

#Load function to extract cell distribution positional data and save it in xlsx
create_2d_matrix_density <- function(panel_data,region_type,name){
  sheets <- c()
  
  #Create two dimensional array
  vector1     <- c(NA)
  vector2     <- c(NA)
  panel_table <- array(c(vector1, vector2), dim = c(21, 21, 1))
  
  #Create sorted bin table
  xy_nr <- 1 #Start parameter (first row in panel_data)
  for(x_nr in 1:max(panel_data$xbin)){
    for(y_nr in 1:max(panel_data$ybin)){
      if(xy_nr<=length(panel_data$avg_ndensity)){
        if(x_nr==panel_data$xbin[xy_nr] &&  y_nr==panel_data$ybin[xy_nr]){
          panel_table[y_nr,x_nr,1] <- as.numeric(panel_data$avg_ndensity[xy_nr])
          xy_nr = xy_nr+1 #Skip to next row
        }
      }
    }
  }
  
  #Invert y-axis (inversion is needed to correctly display cell distribution)
  xbin_length <- length(panel_table[1,,1])
  for(xbin_col in 1:xbin_length){
    panel_table[,xbin_col,1] <- rev(panel_table[,xbin_col,1])
  }
  
  #Convert array into dataframe
  panel_dfs <- data.frame("Bin"=length(panel_table[,1,1]):1)
  for(xbin_col in 1:xbin_length){
    panel_df  <- data.frame(x=panel_table[,xbin_col,1]) 
    panel_dfs <- cbind(panel_dfs,panel_df)
  }
  
  #Change column names
  names(panel_dfs) <- c("Bin",1:xbin_length)
  
  #Create sheet
  sheet <- list(panel_dfs)
  sheets<- c(sheets,sheet)  
  
  #Label the sheets with region names
  names(sheets) <- region_type
  
  #Create xlsx summary table
  write_xlsx(sheets,name)
}

create_average_plots <- function(DistancePlot,region,age,slices,name){
  
  #DistancePlot <-DistancePlot.Table
  #region       <-"MSO"
  #age          <-5
  #slices       <-c("P05_15.03.21_T1_OT1_Schnitt6_l","P05_24.03.21_T2_OT1_Schnitt3_l","P05_24.03.21_T2_OT1_Schnitt3_r","P05_24.03.21_T2_OT1_Schnitt4_r","P05_24.04.21_T1_OT1_Schnitt3_r")
  #name         <-"P5_MSO_Average"
  
  
  #-----------------------------------------------------------------------------
  #Bin the region/age specific data by slice
  #-----------------------------------------------------------------------------
  
  #Extract "Region" Data from dataset
  DistancePlot <- subset(DistancePlot, Region==region)

  #Extract "Age" Data from dataset
  DistancePlot <- subset(DistancePlot, Age==age)
  
  all_binned_data <- data.frame(xbin=integer(),
                                ybin=integer(),
                                ncells=integer(),
                                ndensity=integer(),
                                Slice=integer(),
                                Age=integer(),
                                Region=integer())
  
  for(slice in unique(DistancePlot$Slice_Name)){
    
    SliceData <- subset(DistancePlot, Slice_Name==slice)
    
    #Bin data by 21 bins (x/y-Axis)[1 bin = 5% coverage]
    binned_data <- data.frame(xbin=integer(),ybin=integer(),ncells=integer())
    for(xb in 1:21){#xb<-1
      for(yb in 1:21){#yb<-1
        
        #Get limiter
        xmin <- (xb*5)-5
        xmax <- xb*5
        ymin <- (yb*5)-5 
        ymax <- yb*5
        
        find_cell <- subset(SliceData, X_Perc >= xmin & X_Perc < xmax & Y_Perc >= ymin & Y_Perc < ymax)
        ncell     <- length(find_cell$X_Perc)
        binned_data[nrow(binned_data)+1,] <- c(xb,yb,ncell)
      }
    }
    
    #Add ndensity to dataframe [number of cells / highest number of cells per bin intersect]
    binned_data$ndensity <- binned_data$ncells/max(binned_data$ncells)
    binned_data$Slice    <- slice
    binned_data$Region   <- region
    binned_data$Age      <- age
    all_binned_data      <- rbind(all_binned_data,binned_data)
  }
  
  #Label
  label <- paste(name,"binned_data.xlsx",sep="_")
  label <- paste("output_average/binned_data",label,sep="/")
  
  #Create xlsx summary table
  write_xlsx(all_binned_data, label)
  
  #-----------------------------------------------------------------------------
  
  #-----------------------------------------------------------------------------
  # Convert binned slice data into average data per region and age group -------
  #-----------------------------------------------------------------------------
  
  Core <- all_binned_data
  
  #Create new dataframe
  Average <- data.frame(slices=integer(),
                        xbin=integer(),
                        ybin=integer(),
                        count=integer(),
                        ndensity=integer(),
                        sum_count=integer(),
                        sum_ndensity=integer(),
                        avg_count=integer(),
                        avg_ndensity=integer())
  
  #Get number of slices
  slice_nr <- length(unique(Core$Slice))
  
  #Erase empty fields
  Core  <- subset(Core, ndensity!=0)
  
  #Calculate average density
  for(x in 1:21){#x<-1
    for(y in 1:21){#y<-1
      x_set  <- subset(Core, xbin==x)
      xy_set <- subset(x_set, ybin==y)
      if(length(xy_set$xbin)==0){
        slices       =0
        count        =0
        ndensity     =0
        sum_count    =0
        sum_ndensity =0
        avg_count    =0
        avg_ndensity =0
      }else{
        if(length(xy_set$xbin)>1){
          slices       <- paste(unique(xy_set$Slice), collapse = ",")
          count        <- paste(as.character(xy_set$ncells), collapse = ",")
          ndensity     <- paste(as.character(xy_set$ndensity), collapse = ",")
          sum_count    <- sum(xy_set$ncells)
          sum_ndensity <- sum(xy_set$ndensity)
        }else{
          slices       <- unique(xy_set$Slice)
          count        <- xy_set$ncells
          ndensity     <- xy_set$ndensity 
          sum_count    <- xy_set$ncells
          sum_ndensity <- xy_set$ndensity
        }
        avg_count    <- sum(xy_set$ncells)/slice_nr
        avg_ndensity <- sum(xy_set$ndensity)/slice_nr
      }
      Average[nrow(Average)+1,] <- c(slices,x,y,count,ndensity,sum_count,sum_ndensity,avg_count,avg_ndensity)
    }
  }
  
  #Save average dataset
  Average$xbin         <- as.numeric(Average$xbin) 
  Average$ybin         <- as.numeric(Average$ybin)
  Average$avg_count    <- as.numeric(Average$avg_count)
  Average$avg_ndensity <- as.numeric(Average$avg_ndensity)
  
  #Label
  label <- paste("output_average/averaged_data/",name,"_raw_data.xlsx",sep="")
  
  #Create xlsx summary table
  write_xlsx(Average, label)
  
  #-----------------------------------------------------------------------------

  #-----------------------------------------------------------------------------
  #Plot average data -----------------------------------------------------------
  #-----------------------------------------------------------------------------
  
  #Define plot labels
  plot_title = paste(age,region,"average_ndensity",sep="_")
  plot_x     = 'Cell distribution (%) X-Axis (medial(0)->lateral(100))'
  plot_y     = 'Cell distribution (%) Y-Axis (ventral(0)->dorsal(100))'
  
  # Necessary to put Age names into the facet labels
  Age_names <- as_labeller(
    c(`5` = "P5", `7` = "P7",`10` = "P10", 
      `12` = "P12",`14` = "P14", `17` = "P17", `21` = "P21", `36` = "P34-38", `54` = "P50-59"))
  
  #Generate colorpalette (package: colorspace)
  #library(colorspace)
  #wb_test <- sequential_hcl(10, palette = "Grays", rev = TRUE)
  
  color_wb = c("#F9F9F9", "#EBEBEB", "#D6D6D6", "#BEBEBE", "#A4A4A4", 
               "#898989", "#6D6D6D", "#515151", "#363636", "#1B1B1B")
  
  
  #Plot heatmap of normalized cell location
  
  #Filter slice data
  Average_mod <- subset(Average,avg_ndensity!=0.0)
  
  Average_mod$Age    <- age
  Average_mod$Region <- region

  Average_mod %>%
    mutate(Region = factor(Region, levels=region)) %>%
    ggplot(aes(x = xbin, y = ybin)) +
    geom_tile(aes(fill = avg_ndensity)) +
    scale_fill_gradientn(colours = color_wb,
                         values  = seq(0, 1, by = 0.1),
                         breaks  = seq(0, 1, by = 0.1),   # Hier definieren Sie die Ticks
                         labels  = seq(0, 1, by = 0.1),
                         limits = c(0, 1)) + # Hier definieren Sie die Beschriftungen der Ticks
    facet_grid(Age ~ Region, labeller = labeller(Age = Age_names))+
    coord_fixed(ratio = 1) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank())+
    labs(title = plot_title,
         x = plot_x,
         y = plot_y,
         fill = "ndensity") +
    scale_x_continuous(#limits = c(-5, 105),
      breaks = seq(0, 20, by = 5),
      labels = seq(0, 100, by = 25)) +
    scale_y_continuous(#limits = c(-5, 105),
      breaks = seq(0, 20, by = 5),
      labels = seq(0, 100, by = 25)) +
    theme(axis.text.x = element_text(angle = 90)) 
  label <- paste(name,"average_ndensity.tiff",sep="_")
  label <- paste("output_average/average_plot",label,sep="/")
  ggsave(label, dpi=150, width=5, height=5)
  
  #-----------------------------------------------------------------------------
  
  #Create 2d matrix average ndensity grid (xlsx) -------------------------------
  
  #Label
  label <- paste(name,"bin_avg_ndensity.xlsx",sep="_")
  label <- paste("output_average/xlsx_grid/",label,sep="")
  
  create_2d_matrix_density(Average_mod,region,label)
  #-----------------------------------------------------------------------------
}

################################################################################

################################################################################
# Create average density plots and xlsx density grids ##########################
################################################################################

#Create main output folder
out.path <- FullPath("output_average")
dir.create(out.path, showWarnings = FALSE)
out.path <- FullPath("output_average/binned_data")
dir.create(out.path, showWarnings = FALSE)
out.path <- FullPath("output_average/averaged_data")
dir.create(out.path, showWarnings = FALSE)
out.path <- FullPath("output_average/average_plot")
dir.create(out.path, showWarnings = FALSE)
out.path <- FullPath("output_average/xlsx_grid")
dir.create(out.path, showWarnings = FALSE)

################################################################################
################################################################################

##### MSO ######################################################################

#Create average plot for P5 MSO data subset ------------------------------------
region <- "MSO"
age    <- 5
slices <- c("P05_15.03.21_T1_OT1_Schnitt6_l","P05_24.03.21_T2_OT1_Schnitt3_l","P05_24.03.21_T2_OT1_Schnitt3_r",
            "P05_24.03.21_T2_OT1_Schnitt4_r","P05_24.04.21_T1_OT1_Schnitt3_r")
name   <- "P5_MSO_Average"  
create_average_plots(DistancePlot.Table,region,age,slices,name)
      
################################################################################

#Create average plot for P14 MSO data subset -----------------------------------
region <- "MSO"
age    <- 14
slices <- c("P14_11.01.21_T0_OT2_Schnitt1_r","P14_11.01.21_T0_OT2_Schnitt2_r","P14_11.01.21_T0_OT2_Schnitt4_r",
            "P14_11.01.21_T0_OT2_Schnitt5_l","P14_11.01.21_T0_OT2_Schnitt8_l","P14_14.01.21_T2_OT1_Schnitt2_r",
            "P14_14.01.21_T2_OT1_Schnitt3_r","P14_14.01.21_T2_OT1_Schnitt4_l","P14_14.01.21_T2_OT1_Schnitt4_r",
            "P14_14.01.21_T2_OT1_Schnitt5_r","P14_14.01.21_T2_OT1_Schnitt6_l","P14_14.01.21_T2_OT1_Schnitt7_r",
            "P14_14.01.21_T2_OT1_Schnitt8_r","P14_14.01.21_T2_OT1_Schnitt9_l")
name   <- "P14_MSO_Average"  
create_average_plots(DistancePlot.Table,region,age,slices,name)

################################################################################

#Create average plot for P54 MSO data subset -----------------------------------
region <- "MSO"
age    <- 54
slices <- c("P50_04.07.22_T1_OT1_Schnitt1_l","P50_04.07.22_T1_OT1_Schnitt1_r","P50_04.07.22_T1_OT1_Schnitt2_l",
            "P50_04.07.22_T2_OT1_Schnitt1_l","P50_04.07.22_T2_OT1_Schnitt1_r","P50_04.07.22_T2_OT1_Schnitt2_l",
            "P54_08.07.22_T1_OT1_Schnitt1_r","P54_08.07.22_T1_OT1_Schnitt2_l","P54_08.07.22_T1_OT1_Schnitt2_r",
            "P54_08.07.22_T2_OT1_Schnitt1_l","P54_08.07.22_T2_OT1_Schnitt1_r","P54_08.07.22_T2_OT1_Schnitt2_l",
            "P59_30.05.22_T1_OT1_Schnitt1_l","P59_30.05.22_T1_OT1_Schnitt1_r","P59_30.05.22_T1_OT1_Schnitt2_l",
            "P59_30.05.22_T2_OT1_Schnitt1_l","P59_30.05.22_T2_OT1_Schnitt1_r","P59_30.05.22_T2_OT1_Schnitt2_l")
name   <- "P54_MSO_Average"  
create_average_plots(DistancePlot.Table,region,age,slices,name)

################################################################################
################################################################################

##### LSO_VGlut- ###############################################################

#Create average plot for P5 LSO_VGlut- data subset -----------------------------
region <- "LSO_VGlut-"
age    <- 5
slices <- c("P05_15.03.21_T1_OT1_Schnitt6_l","P05_24.03.21_T2_OT1_Schnitt3_l","P05_24.03.21_T2_OT1_Schnitt3_r",
            "P05_24.03.21_T2_OT1_Schnitt4_r","P05_24.04.21_T1_OT1_Schnitt3_r")
name   <- "P5_LSO_VGlut-_Average"  
create_average_plots(DistancePlot.Table,region,age,slices,name)

################################################################################

#Create average plot for P14 LSO_VGlut- data subset ----------------------------
region <- "LSO_VGlut-"
age    <- 14
slices <- c("P14_11.01.21_T0_OT2_Schnitt1_r","P14_11.01.21_T0_OT2_Schnitt2_r","P14_11.01.21_T0_OT2_Schnitt4_r",
            "P14_11.01.21_T0_OT2_Schnitt5_l","P14_11.01.21_T0_OT2_Schnitt8_l","P14_14.01.21_T2_OT1_Schnitt2_r",
            "P14_14.01.21_T2_OT1_Schnitt3_r","P14_14.01.21_T2_OT1_Schnitt4_l","P14_14.01.21_T2_OT1_Schnitt4_r",
            "P14_14.01.21_T2_OT1_Schnitt5_r","P14_14.01.21_T2_OT1_Schnitt6_l","P14_14.01.21_T2_OT1_Schnitt7_r",
            "P14_14.01.21_T2_OT1_Schnitt8_r","P14_14.01.21_T2_OT1_Schnitt9_l")
name   <- "P14_LSO_VGlut-_Average"  
create_average_plots(DistancePlot.Table,region,age,slices,name)

################################################################################

#Create average plot for P54 LSO_VGlut- data subset ----------------------------
region <- "LSO_VGlut-"
age    <- 54
slices <- c("P50_04.07.22_T1_OT1_Schnitt1_l","P50_04.07.22_T1_OT1_Schnitt1_r","P50_04.07.22_T1_OT1_Schnitt2_l",
            "P50_04.07.22_T2_OT1_Schnitt1_l","P50_04.07.22_T2_OT1_Schnitt1_r","P50_04.07.22_T2_OT1_Schnitt2_l",
            "P54_08.07.22_T1_OT1_Schnitt1_r","P54_08.07.22_T1_OT1_Schnitt2_l","P54_08.07.22_T1_OT1_Schnitt2_r",
            "P54_08.07.22_T2_OT1_Schnitt1_l","P54_08.07.22_T2_OT1_Schnitt1_r","P54_08.07.22_T2_OT1_Schnitt2_l",
            "P59_30.05.22_T1_OT1_Schnitt1_l","P59_30.05.22_T1_OT1_Schnitt1_r","P59_30.05.22_T1_OT1_Schnitt2_l",
            "P59_30.05.22_T2_OT1_Schnitt1_l","P59_30.05.22_T2_OT1_Schnitt1_r","P59_30.05.22_T2_OT1_Schnitt2_l")
name   <- "P54_LSO_VGlut-_Average"  
create_average_plots(DistancePlot.Table,region,age,slices,name)

################################################################################
################################################################################

##### LSO_VGlut+ ###############################################################

#Create average plot for P5 LSO_VGlut+ data subset -----------------------------
region <- "LSO_VGlut+"
age    <- 5
slices <- c("P05_15.03.21_T1_OT1_Schnitt6_l","P05_24.03.21_T2_OT1_Schnitt3_l","P05_24.03.21_T2_OT1_Schnitt3_r",
            "P05_24.03.21_T2_OT1_Schnitt4_r","P05_24.04.21_T1_OT1_Schnitt3_r")
name   <- "P5_LSO_VGlut+_Average"  
create_average_plots(DistancePlot.Table,region,age,slices,name)

################################################################################

#Create average plot for P14 LSO_VGlut+ data subset ----------------------------
region <- "LSO_VGlut+"
age    <- 14
slices <- c("P14_11.01.21_T0_OT2_Schnitt1_r","P14_11.01.21_T0_OT2_Schnitt2_r","P14_11.01.21_T0_OT2_Schnitt4_r",
            "P14_11.01.21_T0_OT2_Schnitt5_l","P14_11.01.21_T0_OT2_Schnitt8_l","P14_14.01.21_T2_OT1_Schnitt2_r",
            "P14_14.01.21_T2_OT1_Schnitt3_r","P14_14.01.21_T2_OT1_Schnitt4_l","P14_14.01.21_T2_OT1_Schnitt4_r",
            "P14_14.01.21_T2_OT1_Schnitt5_r","P14_14.01.21_T2_OT1_Schnitt6_l","P14_14.01.21_T2_OT1_Schnitt7_r",
            "P14_14.01.21_T2_OT1_Schnitt8_r","P14_14.01.21_T2_OT1_Schnitt9_l")
name   <- "P14_LSO_VGlut+_Average"  
create_average_plots(DistancePlot.Table,region,age,slices,name)

################################################################################

#Create average plot for P54 LSO_VGlut+ data subset ----------------------------
region <- "LSO_VGlut+"
age    <- 54
slices <- c("P50_04.07.22_T1_OT1_Schnitt1_l","P50_04.07.22_T1_OT1_Schnitt1_r","P50_04.07.22_T1_OT1_Schnitt2_l",
            "P50_04.07.22_T2_OT1_Schnitt1_l","P50_04.07.22_T2_OT1_Schnitt1_r","P50_04.07.22_T2_OT1_Schnitt2_l",
            "P54_08.07.22_T1_OT1_Schnitt1_r","P54_08.07.22_T1_OT1_Schnitt2_l","P54_08.07.22_T1_OT1_Schnitt2_r",
            "P54_08.07.22_T2_OT1_Schnitt1_l","P54_08.07.22_T2_OT1_Schnitt1_r","P54_08.07.22_T2_OT1_Schnitt2_l",
            "P59_30.05.22_T1_OT1_Schnitt1_l","P59_30.05.22_T1_OT1_Schnitt1_r","P59_30.05.22_T1_OT1_Schnitt2_l",
            "P59_30.05.22_T2_OT1_Schnitt1_l","P59_30.05.22_T2_OT1_Schnitt1_r","P59_30.05.22_T2_OT1_Schnitt2_l")
name   <- "P54_LSO_VGlut+_Average"  
create_average_plots(DistancePlot.Table,region,age,slices,name)

################################################################################
################################################################################

##### MNTB #####################################################################

#Create average plot for P5 MNTB data subset ------------------------------------
region <- "MNTB"
age    <- 5
slices <- c("P05_15.03.21_T1_OT1_Schnitt6_l","P05_24.03.21_T2_OT1_Schnitt3_l","P05_24.03.21_T2_OT1_Schnitt3_r",
            "P05_24.03.21_T2_OT1_Schnitt4_r","P05_24.04.21_T1_OT1_Schnitt3_r")
name   <- "P5_MNTB_Average"  
create_average_plots(DistancePlot.Table,region,age,slices,name)

################################################################################

#Create average plot for P14 MNTB data subset -----------------------------------
region <- "MNTB"
age    <- 14
slices <- c("P14_11.01.21_T0_OT2_Schnitt1_r","P14_11.01.21_T0_OT2_Schnitt2_r","P14_11.01.21_T0_OT2_Schnitt4_r",
            "P14_11.01.21_T0_OT2_Schnitt5_l","P14_11.01.21_T0_OT2_Schnitt8_l","P14_14.01.21_T2_OT1_Schnitt2_r",
            "P14_14.01.21_T2_OT1_Schnitt3_r","P14_14.01.21_T2_OT1_Schnitt4_l","P14_14.01.21_T2_OT1_Schnitt4_r",
            "P14_14.01.21_T2_OT1_Schnitt5_r","P14_14.01.21_T2_OT1_Schnitt6_l","P14_14.01.21_T2_OT1_Schnitt7_r",
            "P14_14.01.21_T2_OT1_Schnitt8_r","P14_14.01.21_T2_OT1_Schnitt9_l")
name   <- "P14_MNTB_Average"  
create_average_plots(DistancePlot.Table,region,age,slices,name)

################################################################################

#Create average plot for P54 MNTB data subset -----------------------------------
region <- "MNTB"
age    <- 54
slices <- c("P50_04.07.22_T1_OT1_Schnitt1_l","P50_04.07.22_T1_OT1_Schnitt1_r","P50_04.07.22_T1_OT1_Schnitt2_l",
            "P50_04.07.22_T2_OT1_Schnitt1_l","P50_04.07.22_T2_OT1_Schnitt1_r","P50_04.07.22_T2_OT1_Schnitt2_l",
            "P54_08.07.22_T1_OT1_Schnitt1_r","P54_08.07.22_T1_OT1_Schnitt2_l","P54_08.07.22_T1_OT1_Schnitt2_r",
            "P54_08.07.22_T2_OT1_Schnitt1_l","P54_08.07.22_T2_OT1_Schnitt1_r","P54_08.07.22_T2_OT1_Schnitt2_l",
            "P59_30.05.22_T1_OT1_Schnitt1_l","P59_30.05.22_T1_OT1_Schnitt1_r","P59_30.05.22_T1_OT1_Schnitt2_l",
            "P59_30.05.22_T2_OT1_Schnitt1_l","P59_30.05.22_T2_OT1_Schnitt1_r","P59_30.05.22_T2_OT1_Schnitt2_l")
name   <- "P54_MNTB_Average"  
create_average_plots(DistancePlot.Table,region,age,slices,name)

################################################################################
################################################################################

##### Kontrolle ################################################################

#Create average plot for P5 Kontrolle data subset ------------------------------
region <- "Kontrolle"
age    <- 5
slices <- c("P05_15.03.21_T1_OT1_Schnitt6_l","P05_24.03.21_T2_OT1_Schnitt3_l","P05_24.03.21_T2_OT1_Schnitt3_r",
            "P05_24.03.21_T2_OT1_Schnitt4_r","P05_24.04.21_T1_OT1_Schnitt3_r")
name   <- "P5_Kontrolle_Average"  
create_average_plots(DistancePlot.Table,region,age,slices,name)

################################################################################

#Create average plot for P14 Kontrolle data subset -----------------------------------
region <- "Kontrolle"
age    <- 14
slices <- c("P14_11.01.21_T0_OT2_Schnitt1_r","P14_11.01.21_T0_OT2_Schnitt2_r","P14_11.01.21_T0_OT2_Schnitt4_r",
            "P14_11.01.21_T0_OT2_Schnitt5_l","P14_11.01.21_T0_OT2_Schnitt8_l","P14_14.01.21_T2_OT1_Schnitt2_r",
            "P14_14.01.21_T2_OT1_Schnitt3_r","P14_14.01.21_T2_OT1_Schnitt4_l","P14_14.01.21_T2_OT1_Schnitt4_r",
            "P14_14.01.21_T2_OT1_Schnitt5_r","P14_14.01.21_T2_OT1_Schnitt6_l","P14_14.01.21_T2_OT1_Schnitt7_r",
            "P14_14.01.21_T2_OT1_Schnitt8_r","P14_14.01.21_T2_OT1_Schnitt9_l")
name   <- "P14_Kontrolle_Average"  
create_average_plots(DistancePlot.Table,region,age,slices,name)

################################################################################

#Create average plot for P54 Kontrolle data subset -----------------------------------
region <- "Kontrolle"
age    <- 54
slices <- c("P50_04.07.22_T1_OT1_Schnitt1_l","P50_04.07.22_T1_OT1_Schnitt1_r","P50_04.07.22_T1_OT1_Schnitt2_l",
            "P50_04.07.22_T2_OT1_Schnitt1_l","P50_04.07.22_T2_OT1_Schnitt1_r","P50_04.07.22_T2_OT1_Schnitt2_l",
            "P54_08.07.22_T1_OT1_Schnitt1_r","P54_08.07.22_T1_OT1_Schnitt2_l","P54_08.07.22_T1_OT1_Schnitt2_r",
            "P54_08.07.22_T2_OT1_Schnitt1_l","P54_08.07.22_T2_OT1_Schnitt1_r","P54_08.07.22_T2_OT1_Schnitt2_l",
            "P59_30.05.22_T1_OT1_Schnitt1_l","P59_30.05.22_T1_OT1_Schnitt1_r","P59_30.05.22_T1_OT1_Schnitt2_l",
            "P59_30.05.22_T2_OT1_Schnitt1_l","P59_30.05.22_T2_OT1_Schnitt1_r","P59_30.05.22_T2_OT1_Schnitt2_l")
name   <- "P54_Kontrolle_Average"  
create_average_plots(DistancePlot.Table,region,age,slices,name)

################################################################################

################################################################################