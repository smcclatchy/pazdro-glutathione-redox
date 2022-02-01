####################################
# Calico GripData
# March 2018
#  GAC - cleaned a version of the 2017 Grip data at 44wks
#        generate plots to show ratio fallacy
#        grip should not be directly adjusted for BW
####################################

setwd("~/Desktop/GripData_Analysis")

library(tidyverse)
library(lubridate)
library(stringr)
library(forcats)
library(RColorBrewer)

#####################################
#  read grip data
####################################

grip.data <- read.csv("GripData_44wk.csv") %>%
  mutate(Born_Date=mdy(Born_Date),
         Test_Date=mdy(Test_Date)) %>%
  mutate(AgeTest = Test_Date-Born_Date) %>%
  mutate(Group=str_sub(MouseID, 4, 5)) %>%
  separate(Generation, into=c("Generation","Parity"), sep="_") %>%
  mutate(Group=as.factor(Group),
         Generation=as.factor(Generation),
         Parity=as.factor(Parity) )

quartz()
ggplot(grip.data, aes(x=BodyWeight, y=AllPaws_T3, color=Group)) +
  geom_point() +
  geom_smooth(method="lm", se=FALSE) +
  #xlim(0,60) + ylim(0,400) +
  geom_abline(intercept=0, slope=(mean(grip.data$AllPaws_T3)/mean(grip.data$BodyWeight)))
#shows ratio fallacy

quartz()
ggplot(grip.data, aes(x=AllPaws_T2, y=AllPaws_T3, color=Group)) +
  geom_point() +
  geom_smooth(method="lm", se=FALSE) +
  #xlim(0,60) + ylim(0,400) +
  geom_abline(intercept=0, slope=(mean(grip.data$AllPaws_T3)/mean(grip.data$AllPaws_T2)))
# shows regression to the mean effect



grip.data <- mutate(grip.data,
       ForePaws_Avg = rowMeans(select(grip.data, starts_with("ForePaws"))),
       AllPaws_Avg = rowMeans(select(grip.data, starts_with("AllPaws"))) )

quartz()
ggplot(grip.data, aes(x=Group, y=AllPaws_Avg, color=Group)) +
  geom_boxplot() +
  geom_point(position=position_jitter(width=0.25, height=0)) +
  ggtitle("Grip Strength")


quartz()
ggplot(grip.data, aes(x=ForePaws_Avg, y=AllPaws_Avg, color=Group)) +
  geom_point() +
  geom_smooth(method="lm", se=FALSE) +
  geom_abline(intercept=0, slope=1)
# AllPaws = ForePaws + 55g

quartz()
ggplot(grip.data, aes(x=BodyWeight, y=AllPaws_Avg, color=Group)) +
  geom_point() +
  geom_smooth(method="lm", se=FALSE) +
  geom_abline(intercept=0, slope=(mean(grip.data$AllPaws_Avg)/mean(grip.data$BodyWeight)))
#shows ratio fallacy

quartz()
ggplot(grip.data, aes(x=BodyWeight, y=ForePaws_Avg, color=Group)) +
  geom_point() +
  geom_smooth(method="lm", se=FALSE) +
  geom_abline(intercept=0, slope=(mean(grip.data$ForePaws_Avg)/mean(grip.data$BodyWeight)))
#shows ratio fallacy


# calculate adjusted grip using ratio (Adj0)
#         and log-residual (Adj1) methods
grip.data <- mutate(grip.data,
       AllPaws_Adj0 = AllPaws_Avg/BodyWeight,
       AllPaws_Adj1 = residuals(lm(log(AllPaws_Avg) ~ log(BodyWeight), data=grip.data)),
       ForePaws_Adj0 = ForePaws_Avg/BodyWeight,
       ForePaws_Adj1 = residuals(lm(log(ForePaws_Avg) ~ log(BodyWeight), data=grip.data)) )

####
# evaluate effects of treatmemt Group on grip strength
# using different adjustment methods

###
# ratio
anova(lm(AllPaws_Adj0~Group, data=grip.data))
#
quartz()
ggplot(grip.data, aes(x=Group, y=AllPaws_Adj0, color=Group)) +
  geom_boxplot() +
  geom_point(position=position_jitter(width=0.25, height=0)) +
  ggtitle("Ratio Adjusted Grip Strength")

###
# log-residual
anova(lm(AllPaws_Adj1~Group, data=grip.data))
#
quartz()
ggplot(grip.data, aes(x=Group, y=AllPaws_Adj1, color=Group)) +
  geom_boxplot() +
  geom_point(position=position_jitter(width=0.25, height=0))+
  ggtitle("Log-Residual Adjusted Grip Strength")

###
# unadjusted
anova(lm(AllPaws_Avg~Group, data=grip.data))
anova(lm(AllPaws_Adj0~Group, data=grip.data))
anova(lm(AllPaws_Adj1~Group, data=grip.data))
#
quartz()
ggplot(grip.data, aes(x=Group, y=AllPaws_Avg, color=Group)) +
  geom_boxplot() +
  geom_point(position=position_jitter(width=0.25, height=0))

####################################
#####################################

