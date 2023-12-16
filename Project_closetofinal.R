

setwd("/Users/allis/OneDrive/Documents/Bioinformatics_Project") 

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("MSstats") 
BiocManager::install("ggplot2") 

#Load Libraries 

library(MSstats) 
library(ggplot2) 

#Import data. Data can be downloaded through the zip file available at https://github.com/allisonv45/Visualizing_Mass_Spectrometry_Data_Final  

raw.data <- read.csv(file="/Users/allis/OneDrive/Documents/Bioinformatics_Project/Choi2017_DDA_Skyline_input.csv")

annot.data <- read.csv(file="/Users/allis/OneDrive/Documents/Bioinformatics_Project/Choi2017_DDA_Skyline_annotation.csv") 

#View the rows and columns in the data set.

head(raw.data)

head(annot.data) 

#Format the Skyline CSV file using the MSstats function which recognizes the Skyline input file. 

input.skyline <- SkylinetoMSstatsFormat(raw.data, 
                                        annotation=annot.data,
                                        removeProtein_with1Feature = TRUE)

#Process the data so that log intensities are created for each protein 
#using the peak area data, the data is also normalized using this function. 


quant.data <- dataProcess(raw = input.skyline,
                          logTrans = 2, 
                          normalization = 'equalizeMedians',
                          summaryMethod = 'TMP',
                          MBimpute = TRUE, 
                          censoredInt = '0',
                          maxQuantileforCensored = 0.999) 

#Extract the protein level data from the object created above.


protein.level <- (quant.data$ProteinLevelData)

head(quant.data$ProteinLevelData)

#Convert the protein name to an integer and add these values to the dataset. 

pnum <- as.numeric(protein.level$Protein) 

protein.level$ProteinNum <- pnum
print(protein.level)
cat('\n\n')



#Create a plot which is modeled after a Manhattan plot and shows all the individual proteins 
#as integers on the x axis and the protein intensity across all conditions and replicates 
#on the y axis. 

png("alldat_visual.png",width=1000,height=1000)
ggplot(data = protein.level,
       aes(x = ProteinNum,
           y = LogIntensities)) +
  geom_jitter(width = 0.2) + 
  labs(x = "Protein", y = "Log_Intensity") +
  geom_line() +
  geom_hline(yintercept = mean(protein.level$LogIntensities), color="blue") 

dev.off() 

png("alldat_visual_withannot.png",width=1000,height=1000)
ggplot(data = protein.level,
       aes(x = ProteinNum,
           y = LogIntensities)) +
  geom_jitter(width = 0.2) + 
  labs(x = "Protein", y = "Log_Intensity") +
  geom_text(
    label= protein.level$ProteinNum,
    nudge_x=0.45, nudge_y=0.1,
    check_overlap=T) + 
  geom_line() +
  geom_hline(yintercept = mean(protein.level$LogIntensities), color="blue") 

dev.off() 

#View the distibution of data. 

mean(protein.level$LogIntensities) 

length(which(protein.level$LogIntensities < 25)) 


length(which(protein.level$LogIntensities > 25)) 

#Create a QC plot and a profile plot using the dataProcessPlots function. 

dataProcessPlots(data = quant.data, 
                 type="QCplot", 
                 which.Protein = 'allonly',
                 address = ("/Users/allis/OneDrive/Documents/Bioinformatics_Project/All_"))

	
dataProcessPlots(data = quant.data, 
                 type="Profileplot", 
                 featureName="NA",
                 width=7, height=7,
                 which.Protein = 'sp|P39521|FHL1_YEAST',
                 address="/Users/allis/OneDrive/Documents/Bioinformatics_Project/FHL1_")


#Create comparison plot for the protein CCW14 to compare expression across all conditions. 

png("protein_comparision.png",width=1000,height=1000) 
ggplot(protein.level[protein.level$Protein == "sp|O13547|CCW14_YEAST",], aes(Condition, LogIntensities)) + 
  geom_point() + 
  stat_summary(fun.y=mean, geom="point", shape=18,
               size=3, color="red") + 
  labs(title="Protein: CCW14_Yeast Condition Comparision",
        x ="Condition Group", y = "Log Intensities")
dev.off()
  



