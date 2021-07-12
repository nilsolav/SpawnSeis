library(ggplot2)
library(grid)
library(gridExtra)


#library(plyr)
#library(dplyr)
#library(tidyr)

#library(colorspace)
#library(data.table)
#library(ggiraph)

# The purpose of the figures:
# What are the differences between the treatments?
# Are there annual differences?

# Change to csv file location
setwd('D:/DATA/temp')

# Pick data files

# 2020, InnerBay, D2
Data <- data.table::fread("Block1_Treat2_Silent control_D2_Location_InnerBay.csv", sep=";",header=T)
Data$Year = "2020"
Data$Treatment = "Silent"
Data$Location = "InnerBay"

Tmp <- data.table::fread("Block1_Treat3_Boat control_D2_Location_InnerBay.csv", sep=";",header=T)
Tmp$Year = "2020"
Tmp$Treatment = "BoatControl"
Tmp$Location = "InnerBay"
Data <- rbind(Data,Tmp)

Tmp <- data.table::fread("Block2_Treat1_Seismic_D2_Location_InnerBay.csv", sep=";",header=T)
Tmp$Year = "2020"
Tmp$Treatment = "Seismic"
Tmp$Location = "InnerBay"
Data <- rbind(Data,Tmp)

# 2020, OuterBay, D1, Silent and Boat control is missing! Hydrophone replaced with icListen (not worked up)
Tmp <- data.table::fread("Block1_Treat1_Seismic_D1_Location_OuterBay.csv", sep=";",header=T)
Tmp$Year = "2020"
Tmp$Treatment = "Seismic"
Tmp$Location = "OuterBay"
Data <- rbind(Data,Tmp)

# 2021, InnerBay, D15
Tmp <- data.table::fread("Block11_Treat3_Silent control_D15_Location_InnerBay.csv", sep=";",header=T)
Tmp$Year = "2021"
Tmp$Treatment = "Silent"
Tmp$Location = "InnerBay"
Data <- rbind(Data,Tmp)

Tmp <- data.table::fread("Block12_Treat1_Boat control_D15_Location_InnerBay.csv", sep=";",header=T)
Tmp$Year = "2021"
Tmp$Treatment = "BoatControl"
Tmp$Location = "InnerBay"
Data <- rbind(Data,Tmp)

Tmp <- data.table::fread("Block12_Treat3_Seismic_D15_Location_InnerBay.csv", sep=";",header=T)
Tmp$Year = "2021"
Tmp$Treatment = "Seismic"
Tmp$Location = "InnerBay"
Data <- rbind(Data,Tmp)

# 2021, OuterBay, D19
Tmp <- data.table::fread("Block16_Treat1_Silent control_D19_Location_OuterBay.csv", sep=";",header=T)
Tmp$Year = "2021"
Tmp$Treatment = "Silent"
Tmp$Location = "OuterBay"
Data <- rbind(Data,Tmp)

Tmp <- data.table::fread("Block16_Treat3_Boat control_D19_Location_OuterBay.csv", sep=";",header=T)
Tmp$Year = "2021"
Tmp$Treatment = "BoatControl"
Tmp$Location = "OuterBay"
Data <- rbind(Data,Tmp)

Tmp <- data.table::fread("Block16_Treat2_Seismic_D19_Location_OuterBay.csv", sep=";",header=T)
Tmp$Year = "2021"
Tmp$Treatment = "Seismic"
Tmp$Location = "OuterBay"
Data <- rbind(Data,Tmp)

Data$TreatmentYear = paste(Data$Treatment,Data$Year)
n<-cbind(Data$Treatment,Data$Year)
# Overview off all data in one figure (pretty messy :)
ggplot(Data,aes(y=SEL,x=t0, color=TreatmentYear)) + geom_point() + geom_smooth(span=0.1, se=F)
ggplot(Data,aes(y=log(pospeakpressure),x=t0, color=Treatment, shape=Year)) + geom_point()

# SEL Broken into subpanels
Inner2020 <- ggplot(subset(Data, (Location %in% "InnerBay")&(Year %in% "2020")),aes(y=SEL,x=t0, color=TreatmentYear)) + geom_point() + geom_smooth(span=0.1, se=F)
Inner2021 <- ggplot(subset(Data, (Location %in% "InnerBay")&(Year %in% "2021")),aes(y=SEL,x=t0, color=TreatmentYear)) + geom_point() + geom_smooth(span=0.1, se=F)
Outer2020 <- ggplot(subset(Data, (Location %in% "OuterBay")&(Year %in% "2020")),aes(y=SEL,x=t0, color=TreatmentYear)) + geom_point() + geom_smooth(span=0.1, se=F)
Outer2021 <- ggplot(subset(Data, (Location %in% "OuterBay")&(Year %in% "2021")),aes(y=SEL,x=t0, color=TreatmentYear)) + geom_point() + geom_smooth(span=0.1, se=F)

myplot1 <- arrangeGrob(Inner2020, 
                       top = textGrob("(A)", x = unit(0, "npc")
                                      , y   = unit(1, "npc"), just=c("left","top"),
                                      gp=gpar(col="black", fontsize=18)))

myplot2 <- arrangeGrob(Inner2021, top = textGrob("(B)", x = unit(0, "npc")
                                             , y = unit(1, "npc"), just=c("left","top"),
                                             gp=gpar(col="black", fontsize=18)))

myplot3 <- arrangeGrob(Outer2020, top = textGrob("(C)", x = unit(0, "npc")
                                             , y  = unit(1, "npc"), just=c("left","top"),
                                             gp=gpar(col="black", fontsize=18)))

myplot4 <- arrangeGrob(Outer2021, top = textGrob("(D)", x = unit(0, "npc")
                                             , y = unit(1, "npc"), just=c("left","top"),
                                             gp=gpar(col="black",    fontsize=18)))
g <- arrangeGrob(myplot1,myplot2,myplot3,myplot4, ncol=2,nrow=2)
ggsave(file="SEL_1.png", g)

Inner <- ggplot(subset(Data, (Location %in% "InnerBay")),aes(y=SEL,x=t0, color=TreatmentYear)) + geom_point() + geom_smooth(span=0.1, se=F)
Outer <- ggplot(subset(Data, (Location %in% "OuterBay")),aes(y=SEL,x=t0, color=TreatmentYear)) + geom_point() + geom_smooth(span=0.1, se=F)
myplot1 <- arrangeGrob(Inner, top = textGrob("(A) Inner bay", x = unit(0, "npc")
                                                 , y  = unit(1, "npc"), just=c("left","top"),
                                                 gp=gpar(col="black", fontsize=18)))

myplot2 <- arrangeGrob(Outer, top = textGrob("(D) Outer Bay", x = unit(0, "npc")
                                                 , y = unit(1, "npc"), just=c("left","top"),
                                                 gp=gpar(col="black",    fontsize=18)))
g <- arrangeGrob(myplot1,myplot2, ncol=1,nrow=2)
ggsave(file="SEL_2.png", g)


# Pressure
# pospeakpressure broken into subpanels
Inner2020 <- ggplot(subset(Data, (Location %in% "InnerBay")&(Year %in% "2020")),aes(y=pospeakpressure,x=t0, color=Treatment, shape=Year)) + geom_point() + geom_smooth(span=0.1, se=F)
Inner2021 <- ggplot(subset(Data, (Location %in% "InnerBay")&(Year %in% "2021")),aes(y=pospeakpressure,x=t0, color=Treatment, shape=Year)) + geom_point() + geom_smooth(span=0.1, se=F)
Outer2020 <- ggplot(subset(Data, (Location %in% "OuterBay")&(Year %in% "2020")),aes(y=pospeakpressure,x=t0, color=Treatment, shape=Year)) + geom_point() + geom_smooth(span=0.1, se=F)
Outer2021 <- ggplot(subset(Data, (Location %in% "OuterBay")&(Year %in% "2021")),aes(y=pospeakpressure,x=t0, color=Treatment, shape=Year)) + geom_point() + geom_smooth(span=0.1, se=F)

myplot1 <- arrangeGrob(Inner2020, 
                       top = textGrob("(A)", x = unit(0, "npc")
                                      , y   = unit(1, "npc"), just=c("left","top"),
                                      gp=gpar(col="black", fontsize=18)))

myplot2 <- arrangeGrob(Inner2021, top = textGrob("(B)", x = unit(0, "npc")
                                                 , y = unit(1, "npc"), just=c("left","top"),
                                                 gp=gpar(col="black", fontsize=18)))

myplot3 <- arrangeGrob(Outer2020, top = textGrob("(C)", x = unit(0, "npc")
                                                 , y  = unit(1, "npc"), just=c("left","top"),
                                                 gp=gpar(col="black", fontsize=18)))

myplot4 <- arrangeGrob(Outer2021, top = textGrob("(D)", x = unit(0, "npc")
                                                 , y = unit(1, "npc"), just=c("left","top"),
                                                 gp=gpar(col="black",    fontsize=18)))
g <- arrangeGrob(myplot1,myplot2,myplot3,myplot4, ncol=2,nrow=2)

ggsave(file="PeakPressure_1.png", g)


Inner <- ggplot(subset(Data, (Location %in% "InnerBay")),aes(y=pospeakpressure,x=t0, color=TreatmentYear)) + geom_point() + geom_smooth(span=0.1, se=F)
Outer <- ggplot(subset(Data, (Location %in% "OuterBay")),aes(y=pospeakpressure,x=t0, color=TreatmentYear)) + geom_point() + geom_smooth(span=0.1, se=F)
myplot1 <- arrangeGrob(Inner, top = textGrob("(A) Inner bay", x = unit(0, "npc")
                                             , y  = unit(1, "npc"), just=c("left","top"),
                                             gp=gpar(col="black", fontsize=18)))

myplot2 <- arrangeGrob(Outer, top = textGrob("(D) Outer Bay", x = unit(0, "npc")
                                             , y = unit(1, "npc"), just=c("left","top"),
                                             gp=gpar(col="black",    fontsize=18)))
g <- arrangeGrob(myplot1,myplot2, ncol=1,nrow=2)
ggsave(file="PeakPressure_2.png", g)


