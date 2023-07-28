##########################################################################
# Code by: grabriella estevam 
# 
# plots RMSD of multiple kinases 
##########################################################################


library(ggplot2)
library("RColorBrewer")
library(readxl)

##########################################################################
# ----------------- Open RMSD excel file and plot  ----------------------
##########################################################################

RMSD <- read.csv("MET_DMS_RMSDs.csv")

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

plot_RMSDs <-ggplot(RMSD)+
  geom_line(aes(x=RMSD$Residue,y=RMSD$AXL),size = 1.5,color="#CC79A7",alpha=0.5) + #pink
  geom_line(aes(x=RMSD$Residue,y=RMSD$IR),size = 1.5,color="#E69F00",alpha=0.5) + # orange
  geom_line(aes(x=RMSD$Residue,y=RMSD$EPHA3),size = 1.5,color="#009E73",alpha=0.5) + #green
  geom_line(aes(x=RMSD$Residue,y=RMSD$KIT),size = 1.5,color="#999999",alpha=0.5) + #grey
  geom_line(aes(x=RMSD$Residue,y=RMSD$RET),size = 1.5,color="#0072B2",alpha=0.5) + #teal
  geom_line(aes(x=RMSD$Residue,y=RMSD$MET),size = 1.5,color="#0072B2") + #orange
  ylab(expression(paste("RMSD (", ring(A),")")))+
  xlab('Residue')+
  theme_bw() +
  theme(text = element_text(size = 15),
        panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank())

plot(plot_RMSDs)
