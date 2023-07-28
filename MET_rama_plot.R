##########################################################################
# MET_rama_plot.R 
# Code written by Gabriella Estevam @ UCSF
# Plots the Beta-5 turn site in Ramachandran space that is equivalent to the MET P1153 proline lock position
# X and Y represent thie Phi and Psi angles of the specific positions - these values are hard coded
# Contour data was obtained from  Lovell et al. 2003
##########################################################################

# packages 
library(ggplot2)
library(bio3d)
library(cowplot)
library("RColorBrewer")

##########################################################################
# ----------------- Open top 500 PDB ramachandran data -------------------
##########################################################################

# datafiles Lovell et al. 2003, converted from .data to .txt for R dataframe formation 
general <- as.data.frame(read.table("rama500-general.txt"))
proline <- as.data.frame(read.table("rama500-pro.txt"))
pre_proline <- as.data.frame(read.table("rama500-prepro.txt"))
glycine <- as.data.frame(read.table("rama500-gly-sym.txt"))


##########################################################################
# -------------------------------Proline----------------------------------
##########################################################################

#### active MET
pdb <- read.pdb("3R7O")
tor <- torsion.pdb(pdb)
#plot(tor$phi, tor$psi)
#head(tor$tbl)
df_tor <- as.data.frame((tor$tbl))
df_tor <- cbind(pos = rownames(df_tor), df_tor)
rownames(df_tor) <- 1:nrow(df_tor)

df <- data.frame(pos=df_tor$pos, 
                 phi = df_tor$phi,
                 psi = df_tor$psi)


plot <- ggplot(df, aes(x=phi, y=psi))+
  geom_point(color="gray",alpha=0.5)+
  #geom_density_2d(geom="polygon", bins=300, alpha=.9)+
  geom_point(aes(x=-63.86839, y=154.2449), colour="black")+
  ggtitle("MET Active")+
  theme_bw()
  #stat_density_2d(geom = "polygon", alpha = 0.1, bins = 1)
plot(plot)

#### inactive MET
pdb2 <- read.pdb("5HTI")
tor2 <- torsion.pdb(pdb2)
#plot(tor$phi, tor$psi)
#head(tor$tbl)
df2_tor <- as.data.frame((tor2$tbl))
df2_tor <- cbind(pos = rownames(df2_tor), df2_tor)
rownames(df2_tor) <- 1:nrow(df2_tor)

df2 <- data.frame(pos=df2_tor$pos, 
                 phi = df2_tor$phi,
                 psi = df2_tor$psi)

plot2 <- ggplot(df2, aes(x=phi, y=psi))+
  geom_point(color="cadetblue3",alpha=0.5)+
  geom_point(aes(x=-64.12199,y=157.30585), colour="black")+
  ggtitle("MET Inactive")
#stat_density_2d(geom = "polygon", alpha = 0.1, bins = 1)

#### ERK2_V101
pdb3 <- read.pdb("4FMQ")
tor3 <- torsion.pdb(pdb3)
#plot(tor$phi, tor$psi)
#head(tor$tbl)
df3_tor <- as.data.frame((tor3$tbl))
df3_tor <- cbind(pos = rownames(df3_tor), df3_tor)
rownames(df3_tor) <- 1:nrow(df3_tor)

df3 <- data.frame(pos=df3_tor$pos, 
                  phi = df3_tor$phi,
                  psi = df3_tor$psi)

plot3 <- ggplot(df3, aes(x=phi, y=psi))+
  geom_point(color="cadetblue3",alpha=0.5)+
  geom_point(aes(x=-127.09523,131.44025), colour="black")+
  ggtitle("ERK2")+
  theme_bw()
#stat_density_2d(geom = "polygon", alpha = 0.1, bins = 1)
#plot_grid(plot,plot2,plot3,ncol=2, nrow=2)

#### ABL_F311
pdb4 <- read.pdb("3CS9")
tor4 <- torsion.pdb(pdb4)
#plot(tor$phi, tor$psi)
#head(tor$tbl)
df4_tor <- as.data.frame((tor4$tbl))
df4_tor <- cbind(pos = rownames(df4_tor), df4_tor)
rownames(df4_tor) <- 1:nrow(df4_tor)

df4 <- data.frame(pos=df4_tor$pos, 
                  phi = df4_tor$phi,
                  psi = df4_tor$psi)

plot4 <- ggplot(df4, aes(x=phi, y=psi))+
  geom_point(color="lightpink",alpha=0.5)+
  geom_point(aes(x=-70.08368, y=149.0016829), colour="black")+
  ggtitle("ABL")+
  theme_bw()
#stat_density_2d(geom = "polygon", alpha = 0.1, bins = 1)

#### CDK6_B.L94
pdb5 <- read.pdb("2EUF")
tor5 <- torsion.pdb(pdb5)
#plot(tor$phi, tor$psi)
#head(tor$tbl)
df5_tor <- as.data.frame((tor5$tbl))
df5_tor <- cbind(pos = rownames(df5_tor), df5_tor)
rownames(df5_tor) <- 1:nrow(df5_tor)

df5 <- data.frame(pos=df5_tor$pos, 
                  phi = df5_tor$phi,
                  psi = df5_tor$psi)

plot5 <- ggplot(df5, aes(x=phi, y=psi))+
  geom_point(color="peru",alpha=0.5)+
  geom_point(aes(x=-118.13874, y=137.2637164), colour="black")+
  ggtitle("CDK6")+
  theme_bw()
plot(plot5)


#### RON_P1157
pdb6 <- read.pdb("3PLS")
tor6 <- torsion.pdb(pdb6)
#plot(tor$phi, tor$psi)
#head(tor$tbl)
df6_tor <- as.data.frame((tor6$tbl))
df6_tor <- cbind(pos = rownames(df6_tor), df6_tor)
rownames(df6_tor) <- 1:nrow(df6_tor)

df6 <- data.frame(pos=df6_tor$pos, 
                  phi = df6_tor$phi,
                  psi = df6_tor$psi)

plot6 <- ggplot(df6, aes(x=phi, y=psi))+
  geom_point(color="peru",alpha=0.5)+
  geom_point(aes(x=-59.09718,y=143.947550), colour="black")+
  ggtitle("RON")+
  theme_bw()
plot(plot6)

#### AXL_P616
pdb7 <- read.pdb("5U6B")
tor7 <- torsion.pdb(pdb7)
#plot(tor$phi, tor$psi)
#head(tor$tbl)
df7_tor <- as.data.frame((tor7$tbl))
df7_tor <- cbind(pos = rownames(df7_tor), df7_tor)
rownames(df7_tor) <- 1:nrow(df7_tor)

df7 <- data.frame(pos=df7_tor$pos, 
                  phi = df7_tor$phi,
                  psi = df7_tor$psi)

plot7 <- ggplot(df7, aes(x=phi, y=psi))+
  geom_point(color="peru",alpha=0.5)+
  geom_point(aes(x=-70.64727,y=146.827600), colour="green")+
  ggtitle("AXL")+
  theme_bw()
plot(plot7)

#### MERTK_P667
pdb8 <- read.pdb("7AB0")
tor8 <- torsion.pdb(pdb8)
#plot(tor$phi, tor$psi)
#head(tor$tbl)
df8_tor <- as.data.frame((tor8$tbl))
df8_tor <- cbind(pos = rownames(df8_tor), df8_tor)
rownames(df8_tor) <- 1:nrow(df8_tor)

df8 <- data.frame(pos=df8_tor$pos, 
                  phi = df8_tor$phi,
                  psi = df8_tor$psi)

plot8 <- ggplot(df8, aes(x=phi, y=psi))+
  geom_point(color="peru",alpha=0.5)+
  geom_point(aes(x=-81.16570, y= 132.74951), colour="black")+
  ggtitle("MERTK")+
  theme_bw()
plot(plot8)

#### Src_I334
pdb9 <- read.pdb("4MXO")
tor9 <- torsion.pdb(pdb9)
#plot(tor$phi, tor$psi)
#head(tor$tbl)
df9_tor <- as.data.frame((tor9$tbl))
df9_tor <- cbind(pos = rownames(df9_tor), df9_tor)
rownames(df9_tor) <- 1:nrow(df9_tor)

df9 <- data.frame(pos=df9_tor$pos, 
                  phi = df9_tor$phi,
                  psi = df9_tor$psi)

plot9 <- ggplot(df9, aes(x=phi, y=psi))+
  geom_point(color="peru",alpha=0.5)+
  geom_point(aes(x=-105.87454, y=144.64973), colour="black")+
  ggtitle("SRC")+
  theme_bw()
plot(plot9)


#### RYK_P407
pdb10 <- read.pdb("6TUA")
tor10 <- torsion.pdb(pdb10)
#plot(tor$phi, tor$psi)
#head(tor$tbl)
df10_tor <- as.data.frame((tor10$tbl))
df10_tor <- cbind(pos = rownames(df10_tor), df10_tor)
rownames(df10_tor) <- 1:nrow(df10_tor)

df10 <- data.frame(pos=df10_tor$pos, 
                  phi = df10_tor$phi,
                  psi = df10_tor$psi)

plot10 <- ggplot(df10, aes(x=phi, y=psi))+
  geom_point(color="peru",alpha=0.5)+
  geom_point(aes(x=-72.26761,y=150.00902), colour="black")+
  ggtitle("RYK")+
  theme_bw()
plot(plot10)

#### RET: 2IVT
pdb11 <- read.pdb("2IVT")
tor11 <- torsion.pdb(pdb11)
#plot(tor$phi, tor$psi)
#head(tor$tbl)
df11_tor <- as.data.frame((tor11$tbl))
df11_tor <- cbind(pos = rownames(df11_tor), df11_tor)
rownames(df11_tor) <- 1:nrow(df11_tor)

df11 <- data.frame(pos=df11_tor$pos, 
                   phi = df11_tor$phi,
                   psi = df11_tor$psi)

plot11 <- ggplot(df11, aes(x=phi, y=psi))+
  geom_point(color="peru",alpha=0.5)+
  geom_point(aes(x=-69.00735, y= 132.072960), colour="black")+
  ggtitle("RET")+
  theme_bw()
plot(plot11)

#### KIT: 1PKG
pdb12 <- read.pdb("1PKG")
tor12 <- torsion.pdb(pdb12)
#plot(tor$phi, tor$psi)
#head(tor$tbl)
df12_tor <- as.data.frame((tor12$tbl))
df12_tor <- cbind(pos = rownames(df12_tor), df12_tor)
rownames(df12_tor) <- 1:nrow(df12_tor)

df12 <- data.frame(pos=df12_tor$pos, 
                   phi = df12_tor$phi,
                   psi = df12_tor$psi)

plot12 <- ggplot(df12, aes(x=phi, y=psi))+
  geom_point(color="peru",alpha=0.5)+
  geom_point(aes(x=-63.71495, y=118.8728049), colour="black")+
  ggtitle("KIT")+
  theme_bw()
plot(plot12)

#### EPHA3: 2QO9 
pdb13 <- read.pdb("2QO9")
tor13 <- torsion.pdb(pdb13)
#plot(tor$phi, tor$psi)
#head(tor$tbl)
df13_tor <- as.data.frame((tor13$tbl))
df13_tor <- cbind(pos = rownames(df13_tor), df13_tor)
rownames(df13_tor) <- 1:nrow(df13_tor)

df13 <- data.frame(pos=df13_tor$pos, 
                   phi = df13_tor$phi,
                   psi = df13_tor$psi)

plot13 <- ggplot(df13, aes(x=phi, y=psi))+
  geom_point(color="peru",alpha=0.5)+
  geom_point(aes(x=-65.01478,y=143.25491), colour="black")+
  ggtitle("EPHA3")+
  theme_bw()
plot(plot13)


#### IR: 4XLV 
pdb14 <- read.pdb("4XLV")
tor14 <- torsion.pdb(pdb14)
#plot(tor$phi, tor$psi)
#head(tor$tbl)
df14_tor <- as.data.frame((tor14$tbl))
df14_tor <- cbind(pos = rownames(df14_tor), df14_tor)
rownames(df14_tor) <- 1:nrow(df14_tor)

df14 <- data.frame(pos=df14_tor$pos, 
                   phi = df14_tor$phi,
                   psi = df14_tor$psi)

plot14 <- ggplot(df14, aes(x=phi, y=psi))+
  geom_point(color="peru",alpha=0.5)+
  geom_point(aes(x=-73.0584,y=126.0359), colour="black")+
  ggtitle("IR")+
  theme_bw()
plot(plot14)


#### EGFR: 2JIT 
pdb15 <- read.pdb("2JIT")
tor15 <- torsion.pdb(pdb15)
#plot(tor$phi, tor$psi)
#head(tor$tbl)
df15_tor <- as.data.frame((tor15$tbl))
df15_tor <- cbind(pos = rownames(df15_tor), df15_tor)
rownames(df15_tor) <- 1:nrow(df15_tor)

df15 <- data.frame(pos=df15_tor$pos, 
                   phi = df15_tor$phi,
                   psi = df15_tor$psi)

plot15 <- ggplot(df15, aes(x=phi, y=psi))+
  geom_point(color="peru",alpha=0.5)+
  geom_point(aes(x=-45.16894, y=131.383636), colour="black")+
  ggtitle("FLT3")+
  theme_bw()
plot(plot15)

#### FLT3 : 4HVS
pdb16 <- read.pdb("4HVS")
tor16 <- torsion.pdb(pdb16)
#plot(tor$phi, tor$psi)
#head(tor$tbl)
df16_tor <- as.data.frame((tor16$tbl))
df16_tor <- cbind(pos = rownames(df16_tor), df16_tor)
rownames(df16_tor) <- 1:nrow(df16_tor)

df16 <- data.frame(pos=df16_tor$pos, 
                   phi = df16_tor$phi,
                   psi = df16_tor$psi)

plot16 <- ggplot(df16, aes(x=phi, y=psi))+
  geom_point(color="peru",alpha=0.5)+
  geom_point(aes(x=-62.77194,y= 124.095419), colour="black")+
  ggtitle("FLT3")+
  theme_bw()
plot(plot16)

#### FGFR2: 6V6Q
pdb17 <- read.pdb("6V6Q")
tor17 <- torsion.pdb(pdb17)
#plot(tor$phi, tor$psi)
#head(tor$tbl)
df17_tor <- as.data.frame((tor17$tbl))
df17_tor <- cbind(pos = rownames(df17_tor), df17_tor)
rownames(df17_tor) <- 1:nrow(df17_tor)

df17 <- data.frame(pos=df17_tor$pos, 
                   phi = df17_tor$phi,
                   psi = df17_tor$psi)

plot17 <- ggplot(df17, aes(x=phi, y=psi))+
  geom_point(color="peru",alpha=0.5)+
  geom_point(aes(x=-69.47166, y= 120.56788), colour="black")+
  ggtitle("FGFR2")+
  theme_bw()
plot(plot17)

#plot_grid(plot,plot3,plot5,plot4,ncol=4, nrow=1)


# General=c(0, 0.0005, 0.02, 1)
##### overlay plots with contour of general

#breaks general=c(0, 0.0005, 0.02, 1)
#breaks glycine=c(0, 0.002,  0.02, 1)
#breaks proline=c(0, 0.002,  0.02, 1),
#breaks pre_proline=c(0, 0.002,  0.02, 1)))

plot_general_and_proline_contour <- ggplot()+
  geom_contour(aes(x=general$V1, y=general$V2, z=general$V3), breaks=c(0, 0.0005, 0.02, 1), color="grey",size=1)+
  geom_contour(aes(x=proline$V1, y=proline$V2, z=proline$V3), breaks=c(0, 0.002, 0.02, 1), color="orange",size=1)+#top 500 general rama data  
  #geom_point(aes(x=-62.77194,y= 124.095419), color="black",size=2.5)+#FLT3
  #geom_point(aes(x=-45.16894, y=131.383636), color="black",size=2.5)+#EGFR
  geom_point(aes(x=-70.08368, y=149.0016829), color="#FFCCD5",size=2.5)+#ABL1
  geom_point(aes(x=-124.27607,y=155.2884189), colour="#FF8FA3",size=2.5)+#ERK2
  geom_point(aes(x=-118.13874, y=137.2637164),color="#FF4D6D",size=2.5)+ #CDK6
  geom_point(aes(x=-105.87454, y=144.64973), colour="#A4133C",size=2.5)+#Src
  geom_point(aes(x=-63.71495, y=118.8728049), colour="#A9D6E5",size=2.5)+#KIT
  geom_point(aes(x=-65.01478, y=143.25491), colour="#89C2D9",size=2.5)+#EPHA3
  geom_point(aes(x=-73.0584, y=126.0359), colour="#61A5C2",size=2.5)+#IR
  geom_point(aes(x=-69.00735, y=132.072960), colour="#468FAF",size=2.5)+#RET
  geom_point(aes(x=-72.26761,y=150.00902), colour="#2C7DA0",size=2.5)+#RYK
  geom_point(aes(x=-81.16570, y= 132.74951), colour="#2A6F97",size=2.5)+#MERTK
  geom_point(aes(x=-70.64727,y=146.827600), colour="#014F86",size=2.5)+#AXL
  geom_point(aes(x=-59.09718,y=143.947550), colour="#01497C",size=2.5)+#RON
  geom_point(aes(x=-64.12199,y=157.30585), colour="#013A63",size=2.5)+#MET_inactive, grey
  geom_point(aes(x=-63.86839, y=154.2449), colour="#012A4A",size=2.5)+#MET_active, black

  
  xlim(-180,180)+
  ylim(-180,180)+
  xlab(expression(paste(Phi)))+
  ylab(expression(paste(Psi)))+
  geom_vline(xintercept = 0, linetype="dashed")+
  geom_hline(yintercept = 0, linetype="dashed")+
  theme_bw() +
  theme(text = element_text(size = 20),
        axis.title.y = element_text(vjust = -3),
        panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank())
plot(plot_general_and_proline_contour)


##########################################################################
# ---------------------------- Pre-Proline -------------------------------
##########################################################################


#### active MET pre proline 1152
#### inactive MET pre proline 1152
#### RON pre proline 1156
#### AXL pre proline 615
#### MERTK pre proline 666
#### RYK pre proline 406
#### ERK2 pre proline 100
#### CDK6 pre proline 93
#### Src pre proline 333
  


plot_pre_proline_contour <- ggplot()+
  geom_contour(aes(x=general$V1, y=general$V2, z=general$V3), breaks=c(0, 0.0005, 0.02, 1), color="grey",size=1)+
  geom_contour(aes(x=pre_proline$V1, y=pre_proline$V2, z=pre_proline$V3), breaks=c(0, 0.002, 0.02, 1), color="#56B4E9",size=1)+
  geom_point(aes(x=-74.91299,y=133.66973), color="#FFB3C1",size=2.5)+#ABL1
  geom_point(aes(x=-124.27607,y=155.2884189), colour="#52B788",size=2.5)+#ERK2
  geom_point(aes(x=-141.23014,y=151.722822),color="#2D6A4F",size=2.5)+ #CDK6
  geom_point(aes(x=-64.74834,y=151.841950), colour="#081C15",size=2.5)+#Src
  geom_point(aes(x=-59.65296, y=156.23086), colour="#A9D6E5",size=2.5)+#KIT
  geom_point(aes(x=-63.56893, y=141.77514), colour="#89C2D9",size=2.5)+#EPHA3
  geom_point(aes(x=-84.03941, y=170.06), colour="#61A5C2",size=2.5)+#IR
  geom_point(aes(x=-60.29415, y= 153.58907), colour="#468FAF",size=2.5)+#RET
  geom_point(aes(x=-76.81846,y=128.139676), colour="#2C7DA0",size=2.5)+#RYK
  geom_point(aes(x=-118.46921,y=92.24596),colour="#2A6F97",size=2.5)+#MERTK
  geom_point(aes(x=-90.97048,y=123.39990), colour="#014F86",size=2.5)+#AXL 
  geom_point(aes(x=-81.42869,y=145.8215), colour="#01497C",size=2.5)+#RON
  geom_point(aes(x=-67.00455,y=145.6883), colour="#013A63",size=2.5)+#MET_inactive
  geom_point(aes(x=-63.86839, y=154.2449), colour="#012A4A",size=2.5)+#MET_active, black
  
  
  xlim(-180,180)+
  ylim(-180,180)+
  xlab(expression(paste(Phi)))+
  ylab(expression(paste(Psi)))+
  geom_vline(xintercept = 0, linetype="dashed")+
  geom_hline(yintercept = 0, linetype="dashed")+
  theme_bw() +
  theme(text = element_text(size = 20),
        axis.title.y = element_text(vjust = -2),
        panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank())
plot(plot_pre_proline_contour)
plot_grid(plot_general_and_proline_contour,plot_pre_proline_contour,ncol=2,nrow=1)


##########################################################################
# ----------------------Supplamental figure -----------------------------
##########################################################################

# VGFR1: 3HNG_905
pdb18 <- read.pdb("3HNG")
tor18 <- torsion.pdb(pdb18)
#plot(tor$phi, tor$psi)
#head(tor$tbl)
df18_tor <- as.data.frame((tor18$tbl))
df18_tor <- cbind(pos = rownames(df18_tor), df18_tor)
rownames(df18_tor) <- 1:nrow(df18_tor)

df18 <- data.frame(pos=df18_tor$pos, 
                   phi = df18_tor$phi,
                   psi = df18_tor$psi)

plot18 <- ggplot(df18, aes(x=phi, y=psi))+
  geom_point(color="peru",alpha=0.5)+
  geom_point(aes(x=-69.91061, y= 135.94265), colour="black")+
  ggtitle("VGFR1")+
  theme_bw()
plot(plot18)

# PTK7 : 6VG3_870
pdb19 <- read.pdb("6VG3")
tor19 <- torsion.pdb(pdb19)
#plot(tor$phi, tor$psi)
#head(tor$tbl)
df19_tor <- as.data.frame((tor19$tbl))
df19_tor <- cbind(pos = rownames(df19_tor), df19_tor)
rownames(df19_tor) <- 1:nrow(df19_tor)

df19 <- data.frame(pos=df19_tor$pos, 
                   phi = df19_tor$phi,
                   psi = df19_tor$psi)

plot19 <- ggplot(df19, aes(x=phi, y=psi))+
  geom_point(color="peru",alpha=0.5)+
  geom_point(aes(x=-76.98176, y=144.887614), colour="black")+
  ggtitle("PTK7")+
  theme_bw()
plot(plot19)

# TRKA_585
pdb20 <- read.pdb("7VKO")
tor20 <- torsion.pdb(pdb20)
#plot(tor$phi, tor$psi)
#head(tor$tbl)
df20_tor <- as.data.frame((tor20$tbl))
df20_tor <- cbind(pos = rownames(df20_tor), df20_tor)
rownames(df20_tor) <- 1:nrow(df20_tor)

df20 <- data.frame(pos=df20_tor$pos, 
                   phi = df20_tor$phi,
                   psi = df20_tor$psi)

plot20 <- ggplot(df20, aes(x=phi, y=psi))+
  geom_point(color="peru",alpha=0.5)+
  geom_point(aes(x=-67.16856, y= 138.99119), colour="black")+
  ggtitle("TRKA")+
  theme_bw()
plot(plot20)


# ROR1_548 
pdb21 <- read.pdb("6TU9")
tor21 <- torsion.pdb(pdb21)
#plot(tor$phi, tor$psi)
#head(tor$tbl)
df21_tor <- as.data.frame((tor21$tbl))
df21_tor <- cbind(pos = rownames(df21_tor), df21_tor)
rownames(df21_tor) <- 1:nrow(df21_tor)

df21 <- data.frame(pos=df21_tor$pos, 
                   phi = df21_tor$phi,
                   psi = df21_tor$psi)

plot21 <- ggplot(df21, aes(x=phi, y=psi))+
  geom_point(color="peru",alpha=0.5)+
  geom_point(aes(x=-70.51445, y= 131.34172), colour="black")+
  ggtitle("ROR1")+
  theme_bw()
plot(plot21)

# DDR1 _6FIO_697
pdb22 <- read.pdb("6FIO")
tor22 <- torsion.pdb(pdb22)
#plot(tor$phi, tor$psi)
#head(tor$tbl)
df22_tor <- as.data.frame((tor22$tbl))
df22_tor <- cbind(pos = rownames(df22_tor), df22_tor)
rownames(df22_tor) <- 1:nrow(df22_tor)

df22 <- data.frame(pos=df22_tor$pos, 
                   phi = df22_tor$phi,
                   psi = df22_tor$psi)

plot22 <- ggplot(df22, aes(x=phi, y=psi))+
  geom_point(color="peru",alpha=0.5)+
  geom_point(aes(x=-68.4778, y= 142.08661), colour="black")+
  ggtitle("DDR1")+
  theme_bw()
plot(plot22)

# ROS_2022
pdb23 <- read.pdb("3ZBF")
tor23 <- torsion.pdb(pdb23)
#plot(tor$phi, tor$psi)
#head(tor$tbl)
df23_tor <- as.data.frame((tor23$tbl))
df23_tor <- cbind(pos = rownames(df23_tor), df23_tor)
rownames(df23_tor) <- 1:nrow(df23_tor)

df23 <- data.frame(pos=df23_tor$pos, 
                   phi = df23_tor$phi,
                   psi = df23_tor$psi)

plot23 <- ggplot(df23, aes(x=phi, y=psi))+
  geom_point(color="peru",alpha=0.5)+
  geom_point(aes(x=-64.38938, y=145.6863), colour="black")+
  ggtitle("ROS")+
  theme_bw()
plot(plot23)

# ALK_1192
pdb24 <- read.pdb("7R7K")
tor24 <- torsion.pdb(pdb24)
#plot(tor$phi, tor$psi)
#head(tor$tbl)
df24_tor <- as.data.frame((tor24$tbl))
df24_tor <- cbind(pos = rownames(df24_tor), df24_tor)
rownames(df24_tor) <- 1:nrow(df24_tor)

df24 <- data.frame(pos=df24_tor$pos, 
                   phi = df24_tor$phi,
                   psi = df24_tor$psi)

plot24 <- ggplot(df24, aes(x=phi, y=psi))+
  geom_point(color="peru",alpha=0.5)+
  geom_point(aes(x=-69.02627, y= 144.91186), colour="black")+
  ggtitle("ROS")+
  theme_bw()
plot(plot24)

# EGFR_761
pdb25 <- read.pdb("3GT8")
tor25 <- torsion.pdb(pdb25)
#plot(tor$phi, tor$psi)
#head(tor$tbl)
df25_tor <- as.data.frame((tor25$tbl))
df25_tor <- cbind(pos = rownames(df25_tor), df25_tor)
rownames(df25_tor) <- 1:nrow(df25_tor)

df25 <- data.frame(pos=df25_tor$pos, 
                   phi = df25_tor$phi,
                   psi = df25_tor$psi)

plot25 <- ggplot(df25, aes(x=phi, y=psi))+
  geom_point(color="peru",alpha=0.5)+
  geom_point(aes(x=-120.06168, y=148.462491), colour="black")+
  ggtitle("EGFR")+
  theme_bw()
plot(plot25)

# TIE2_897
pdb26 <- read.pdb("6MWE")
tor26 <- torsion.pdb(pdb26)
#plot(tor$phi, tor$psi)
#head(tor$tbl)
df26_tor <- as.data.frame((tor26$tbl))
df26_tor <- cbind(pos = rownames(df26_tor), df26_tor)
rownames(df26_tor) <- 1:nrow(df26_tor)

df26 <- data.frame(pos=df26_tor$pos, 
                   phi = df26_tor$phi,
                   psi = df26_tor$psi)

plot26 <- ggplot(df26, aes(x=phi, y=psi))+
  geom_point(color="peru",alpha=0.5)+
  geom_point(aes(x=-96.23733, y=145.43312), colour="black")+
  ggtitle("TIE2")+
  theme_bw()
plot(plot26)

# LMR3_234
pdb27 <- read.pdb("6SEQ")
tor27 <- torsion.pdb(pdb27)
#plot(tor$phi, tor$psi)
#head(tor$tbl)
df27_tor <- as.data.frame((tor27$tbl))
df27_tor <- cbind(pos = rownames(df27_tor), df27_tor)
rownames(df27_tor) <- 1:nrow(df27_tor)

df27 <- data.frame(pos=df27_tor$pos, 
                   phi = df27_tor$phi,
                   psi = df27_tor$psi)

plot27 <- ggplot(df27, aes(x=phi, y=psi))+
  geom_point(color="peru",alpha=0.5)+
  geom_point(aes(x=-75.74247, y=149.3489), colour="black")+
  ggtitle("LMR3")+
  theme_bw()
plot(plot27)



plot_general_and_proline_contour_2 <- ggplot()+
  geom_contour(aes(x=general$V1, y=general$V2, z=general$V3), breaks=c(0, 0.0005, 0.02, 1), color="grey",size=1)+
  geom_contour(aes(x=proline$V1, y=proline$V2, z=proline$V3), breaks=c(0, 0.002, 0.02, 1), color="orange",size=1)+#top 500 general rama data  

  geom_point(aes(x=-64.12199,y=157.30585), colour="#012A4A",size=2.5)+#MET_inactive, grey
  geom_point(aes(x=-70.64727,y=146.827600), colour="#013A63",size=2.5)+#AXL
  geom_point(aes(x=-72.26761,y=150.00902), colour="#01497C",size=2.5)+#RsYK
  geom_point(aes(x=-69.00735, y=132.072960), colour="#014F86",size=2.5)+#RET
  geom_point(aes(x=-73.0584, y=126.0359), colour="#2A6F97",size=2.5)+#IR
  geom_point(aes(x=-65.01478, y=143.26491), colour="#2C7DA0",size=2.5)+#EPHA3
  geom_point(aes(x=-63.71495, y=118.8728049), colour="#468FAF",size=2.5)+#KIT
  geom_point(aes(x=-69.02627, y= 144.91186), colour="#61A5C2",size=2.5)+#ALK
  geom_point(aes(x=-64.38938, y=145.6863), colour="#89C2D9",size=2.5)+#ROS
  geom_point(aes(x=-68.4778, y= 142.08661), colour="#A9D6E5",size=2.5)+#DDR1 
  geom_point(aes(x=-70.51445, y= 131.34172), colour="#E9ECEF",size=2.5)+#ROR1
  geom_point(aes(x=-67.16856, y= 138.99119), colour="#DEE2E6",size=2.5)+#TRKA
  geom_point(aes(x=-76.98176, y=144.887614), colour="#CED4DA",size=2.5)+#PTK7
  geom_point(aes(x=-69.91061, y= 135.94265), colour="#ADB5BD",size=2.5)+#VGFR1
  geom_point(aes(x=-69.47166, y= 120.56788), colour="#6C757D",size=2.5)+#FGFR2
  geom_point(aes(x=-75.74247, y=149.3489), colour="#495057",size=2.5)+#LMR3
  geom_point(aes(x=-96.23733, y=145.43312), colour="#343A40",size=2.5)+#TIE2
  geom_point(aes(x=-120.06168, y=148.462491), colour="#212529",size=2.5)+#EGFR
  geom_point(aes(x=-70.08368, y=149.0016829), color="#FFCCD5",size=2.5)+#ABL1
  geom_point(aes(x=-124.27607,y=155.2884189), colour="#FF8FA3",size=2.5)+#ERK2
  geom_point(aes(x=-118.13874, y=137.2637164),color="#FF4D6D",size=2.5)+ #CDK6
  geom_point(aes(x=-105.87454, y=144.64973), colour="#A4133C",size=2.5)+#Src
  
  xlim(-180,180)+
  ylim(-180,180)+
  xlab(expression(paste(Phi)))+
  ylab(expression(paste(Psi)))+
  geom_vline(xintercept = 0, linetype="dashed")+
  geom_hline(yintercept = 0, linetype="dashed")+
  theme_bw() +
  theme(text = element_text(size = 20),
        axis.title.y = element_text(vjust = -3),
        panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank())
plot(plot_general_and_proline_contour_2)



