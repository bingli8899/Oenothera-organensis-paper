# possible package: 
library(modeest)
library(lme4)
library(ggplot2)
library(car)
library(lmerTest)
library(sjPlot)
library(ggpubr)
library(installr)
library(plyr)
library(vegan)
library(TMB)
library(glmmTMB)
library(emmeans)
library(dplyr)
library(boot)
library(lattice)
library(gtable)
library(boot)
library(MASS)
library(performance)
library(Hmsc)
library(tidyverse)

# Based on the comments from reviewer: 
# install.packages("devtools") 
# install.packages("rlang") # Doesn't work 
# install.packages("https://cran.r-project.org/src/contrib/Archive/rlang/rlang_1.0.3.tar.gz", repos = NULL, type="source")
library(devtools)
# install_github("hmsc-r/HMSC")
library(Hmsc)


################## Uploading dataset ##################

# The dataset contains raw  morphological and phenological data: 
df <- read.csv("https://www.dropbox.com/s/six49bn7ncinfpx/2020_Oenothera_Phenology.csv?dl=1")
# The dataset contains the individual plant ID 
Tag <- read.csv("https://www.dropbox.com/s/lpjiaph911jj681/2020_O.organ._Tag.csv?dl=1")
aggregate(Origin ~ Accession.ID, FUN = length, data = Tag) # Appendix table A 
aggregate(Origin ~ Accession.., FUN = length, data = Tag) # Appendix table A

# The dataset contains the summairzed phenological info for each day from df 
day <- read.csv("https://www.dropbox.com/s/loxs9v1nvcdz49l/2020_Phenology.csv?dl=1")

# The dataset contains summarized phenological data from df 
phenf <- read.csv("https://www.dropbox.com/s/wb8uc4nr0dmqfly/Phenology_finaldata.csv?dl=1")

# The dataset contains summarized number of flowers data from df 
fn <- read.csv("https://www.dropbox.com/s/l2eegxpl1v13gtq/2020_Oenothera_number_fl.csv?dl=1")

# Whether each column has the correct format 
is.numeric(df$Flag_number)
is.numeric(df$number_fl)
is.numeric(df$Hypanthium_length)
is.numeric(df$Fl_tube_diamater)
is.numeric(df$style_length)
is.numeric(df$Sugar_content) # No 
is.numeric(df$Nectar_height)
is.numeric(df$Days)
df$Sugar_content <- as.numeric(as.character(df$Sugar_content)) 
is.numeric(df$Sugar_content) 
# Combing two data frame by plant ID, so each plant has its maternal line and origin 
df <- merge(df, Tag, by.x = 'Flag_number', by.y = 'Tag..')
phen <- merge(phenf, Tag, by.x = 'Flag_number', by.y = 'Tag..')





################## The number of flowers ##################

# Change the number of days to the day of the year -- starting from June 21
df$day <- 172 + df$Days # June 21, 2020 is 173 in Julian Date  (leap year)
range(df$day) # Looks good 

# Are data normally distributed? 
pday = ggplot(df, aes(x=Days)) + 
  geom_histogram(colour="black", fill="white")
pday = ggplot(df, aes(x=Days)) +
  geom_histogram(aes(y = ..density..), colour="black", fill="white") +
  stat_function(fun = dnorm, lwd = 2, col = 'red', 
                args = list(mean = mean(df$Days), sd = sd(df$Days)))
pday # normally distributed but slightly skewed 

# Some plot to begin. 
# Number of flowering per day -- using the "day" dataset
p_number_per_day = ggplot(day, aes(x = Julian_day, y = number_.flowering_plants, colour = Origin, group =2)) + 
  geom_line(aes(group = Origin), position = "identity", size  = 1) +
  scale_colour_manual("Origin",   labels = c("Ex situ", "Wild"),
                      values = c("black","dodgerblue3"))+
  theme(
    panel.background = element_rect(fill = "white", colour = "grey",
                                    size = 2, linetype = "solid"),
    panel.grid.major = element_line(), 
    panel.grid.minor = element_line())+
  labs(title = "Number of flowering Plant", size = 13) + 
  xlab("Day of the Year") + ylab("Number of flowering plants")+
  theme(axis.text.x = element_text(size = 14), 
        axis.title.x =element_text(size = 16),
        axis.title.y =element_text(size = 16),
        axis.text.y = element_text(size = 14))+
  theme(plot.title = element_blank())+
  geom_vline(xintercept = 211, color = "red", linetype="dotted", size = 1)+
  geom_vline(xintercept = 248, color = "red", linetype="dotted", size = 1)
p_number_per_day
ggsave("p_number_per_day.tiff", dpi = 700) # The figure is not included in the final manuscript 

# figure of number of flowers per day per plant -- shown on the thesis
p_nperday_perplant = ggplot(day, aes(x = Julian_day, y = AVERAGE.NM, group =2)) + 
  geom_line(aes(group = Origin, colour = Origin), size =1) +
  scale_colour_manual("Origin",   labels = c("Ex situ", "Wild"),
                      values = c("black","dodgerblue3"))+
  theme(
    panel.background = element_rect(fill = "white", colour = "grey",
                                    size = 2, linetype = "solid"),
    panel.grid.major = element_line(), 
    panel.grid.minor = element_line()) + 
  xlab("Day of the year (Julian Date)") + ylab("Average number of flowers per flowering plant")+
  scale_x_continuous(breaks = seq(173,293, by=10))+
  geom_errorbar(aes(y = AVERAGE.NM, ymin = AVERAGE.NM - SE.NM, ymax =  AVERAGE.NM + SE.NM, 
                    colour = Origin), 
                alpha = 0.2, size = 1.5, width = 0) + 
  theme(axis.text.x = element_text(size = 14), 
        axis.title.x =element_text(size = 16),
        axis.title.y =element_text(size = 16),
        axis.text.y = element_text(size = 14))+
  geom_vline(xintercept = 211, color = "red", linetype="dotted", size = 1)+
  geom_vline(xintercept = 248, color = "red", linetype="dotted", size = 1)
p_nperday_perplant 
dev.size("cm")
ggsave("p_nperday_perplant.tiff", dpi = 1000)

# Analysis -- first calculate the three parameters (duration, peak, and first) 
Phen <- ddply(phenf, .(Flag_number)
              , summarise, 
              first = min(Days) + 172,
              duration = max(Days) - min(Days))
print(Phen)

# Calculate peak flowering (the date that 50% of flowers open) 
range(phenf$Flag_number)
is.numeric(phenf$Flag_number)
peak_data = data.frame(matrix(ncol = 2, nrow = 75))
colnames(peak_data) <- c('Flag_number', 'peak')

for (i in phenf$Flag_number) {
    peak = 0 
    data_frame_i <- phenf[phenf$Flag_number == i, ]
    data_frame_i <- data_frame_i %>% mutate(cum = cumsum(number_fl)/sum(number_fl))
    peak = data_frame_i[data_frame_i$cum > 0.5 , 'Days']
    peak_data$Flag_number[i] = i 
    peak_data$peak[i] = peak[1] + 172 # 172 -- Julian Date, and the first number is when it exceeds 50% 
}
print(peak_data) # The dataset has all the peak values. 
# Note that there is no number 66. The label is from 1 to 75, and we have 74 plants. 

Phen <- merge(Phen, peak_data, by.x = 'Flag_number', by.y = 'Flag_number')
print(Phen)

# The three phenological parameters are correlated -- Yes, they are correlated
cor.test(Phen$peak,Phen$first) # p = 0.0007, cor = 0.3833 
cor.test(Phen$peak,Phen$duration) # p = 0.00279, cor = -0.343 
cor.test(Phen$first,Phen$duration) # p<***, cor = -0.712 

# merge the new dataset with accession and add the number of flowers into the dataset 
Phen <- merge(Phen, Tag, by.x = 'Flag_number', by.y = 'Tag..')
dfn <- ddply(df, .(Flag_number), summarise, fn=sum(number_fl))
names(dfn)[1] <- "ID"
Phen <- merge(Phen, dfn, by.x = 'Flag_number', by.y = 'ID')
names(Phen)

# Model - first day of flowering 
first1 <- lmer(first ~ Origin*fn + (1|Accession..), data = Phen)
tab_model(first1) # icc = 0 no random effect 
plot_model(first1, type = "diag")[2] # no random effect
first2 <- lm(first ~ Origin*fn + Accession.., data = Phen)
tab_model(first2) # Accession not significant 
first3 <- lm(first ~ Origin*fn, data = Phen)
first4 <- lm(first ~ Origin + Accession.., data = Phen)
first5 <- lm(first ~ Origin, data = Phen)

# Check models for first day
anova(first3,first2) # P = 0.39 
AIC(first2,first3) # first 3 is better 
AIC(first1,first3) # first 3 is better 
AIC(first3,first4) # first 3 is better
AIC(first3, first5) # first 3 is better 
Anova(first3,type =3) # 
tab_model(first3)

# check model -- okay
par(mfcol = c(2,2))
plot(first3)
m_first <- first3 # This is the model will be used

# Plotting first day
pfirst <- ggplot(aes(x = fn, y = first, color = Origin), data = Phen) + geom_point()+
geom_smooth(aes(linetype = Origin), method = "lm", se = T, alpha = 0.2)+
  theme(
    panel.background = element_rect(fill = "white", colour = "grey",
                                    size = 2, linetype = "solid"),
    panel.grid.major = element_line(), 
    panel.grid.minor = element_line())+
  scale_colour_manual("Origin",   labels = c("Ex situ", "Wild"),
                      values = c("black","dodgerblue3"))+ 
  xlab("Number of flowers") + ylab("First day of flowering") + theme(legend.position = "none")
pfirst

dfirstmean <- Phen
names(dfirstmean)[2] <- "emmean"
mfirst <- as.data.frame(emmeans(first3, specs = "Origin"))
pfirst2 <- ggplot(mfirst, aes(x = Origin, y = emmean), colour = Origin) + 
  geom_jitter(data = dfirstmean, aes(x = Origin, y = emmean, color = Origin), size = 0.7, , alpha = 0.6)+
  geom_point(aes(color = Origin), position = position_dodge(width =0.8), size = 3)+
  geom_errorbar(aes(x = Origin, ymin = emmean - 2*SE, ymax = emmean + 2*SE, color = Origin), position = position_dodge(width =0.8), 
                width = 0.05, size = 1.5)+
  theme(
    panel.background = element_rect(fill = "white", colour = "grey",
                                    size = 2, linetype = "solid"),
    panel.grid.major = element_line(), 
    panel.grid.minor = element_line())+
  scale_colour_manual("Origin",   labels = c("Ex situ", "Wild"),
                      values = c("black","dodgerblue3"))+ 
  xlab("Number of flowers") + ylab("First day of flowering") + theme(legend.position = "none")+
  scale_x_discrete (labels= c("ex situ", "wild"))+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1))+
  scale_colour_manual("Origin",   labels = c("Ex situ", "Wild"),
                      values = c("black","dodgerblue3"))+
  theme(axis.text.x = element_text(size = 12), 
        axis.title.x =element_text(size = 12),
        axis.title.y =element_text(size = 12),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size = 14))
pfirst2

# Model -- duration of flowering 
duration1 <- lmer(duration ~ Origin*fn + (1|Accession..), data = Phen)
tab_model(duration1) # ICC = 0 
plot_model(duration1, type = "diag")[2] # NO random effect
duration2 <- lm(duration ~ Origin*fn + Accession.., data = Phen)
duration3 <- lm(duration ~ Origin*fn , data = Phen)
duration4 <- lm(duration ~ Origin + Accession.., data = Phen)
duration5 <- lm(duration ~ Origin, data = Phen)

# Check models for duration 
anova(duration3,duration2) # P = 0.97 little vairation in maternal line -- duration 3 is chosen
AIC(duration2,duration3) # duration 3 is better
AIC(duration1,duration3) # Duration 3 is better
AIC(duration4,duration3) # model 3 is better 
AIC(duration5,duration3) # model 3 is better 

Anova(duration3,type =3) # p = 0.82
tab_model(duration3)
# check model 
par(mfcol = c(2,2))
plot(duration3)
m_duration <- duration3

# Plotting duration 
pduration <- ggplot(aes(x = fn, y = duration, color = Origin), data = Phen) + geom_point()+
  geom_smooth(aes(linetype = Origin), method = "lm", se = T, alpha = 0.2)+
  theme(
    panel.background = element_rect(fill = "white", colour = "grey",
                                    size = 2, linetype = "solid"),
    panel.grid.major = element_line(), 
    panel.grid.minor = element_line())+
  scale_colour_manual("Origin",   labels = c("Ex situ", "Wild"),
                      values = c("black","dodgerblue3"))+ 
  xlab("Number of flowers") + ylab("Duration of flowering") + theme(legend.position = "none")
pduration

ddurationmean <- Phen
names(ddurationmean)[3] <- "emmean"
mduration <- as.data.frame(emmeans(duration3, specs = "Origin"))
pduration2 <- ggplot(mduration, aes(x = Origin, y = emmean), colour = Origin) + 
  geom_jitter(data = ddurationmean, aes(x = Origin, y = emmean, color = Origin), size = 0.7, , alpha = 0.6)+
  geom_point(aes(color = Origin), position = position_dodge(width =0.8), size = 3)+
  geom_errorbar(aes(x = Origin, ymin = emmean - 2*SE, ymax = emmean + 2*SE, color = Origin), position = position_dodge(width =0.8), 
                width = 0.05, size = 1.5)+
  theme(
    panel.background = element_rect(fill = "white", colour = "grey",
                                    size = 2, linetype = "solid"),
    panel.grid.major = element_line(), 
    panel.grid.minor = element_line())+
  scale_colour_manual("Origin",   labels = c("Ex situ", "Wild"),
                      values = c("black","dodgerblue3"))+ 
  xlab("Number of flowers") + ylab("Duration of flowering") + theme(legend.position = "none")+
  scale_x_discrete (labels= c("ex situ", "wild"))+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1))+
  scale_colour_manual("Origin",   labels = c("Ex situ", "Wild"),
                      values = c("black","dodgerblue3"))+
  theme(axis.text.x = element_text(size = 12), 
        axis.title.x =element_text(size = 12),
        axis.title.y =element_text(size = 12),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size = 14))
pduration2

# Model -- peak day flowering 
peak1 <- lmer(peak ~ Origin*fn + (1|Accession..), data = Phen)
tab_model(peak1) # no random effect 
plot_model(peak1, type = "diag")[2] # no random effect 
peak2 <- lm(peak ~ Origin*fn + Accession.., data = Phen)
tab_model(peak2) # Not significant for accession
peak3 <- lm(peak ~ Origin*fn , data = Phen)
peak4 <- lm(peak ~ Origin + Accession.., data = Phen)
peak5 <- lm(peak ~ Origin, data = Phen)
anova(peak3, peak2) # no difference caused by maternal line 
AIC(peak2,peak3) # peak 3 is better
AIC(peak3,peak1) # peak 3 is better
AIC(peak4,peak3) # peak 3 is better
AIC(peak3,peak5) # peak 3 is better 

Anova(peak3, type =3)
tab_model(peak3)
m_peak <- peak3 # peak 3 is the model we chose based on AIC

# plottig the peak flowering date 
ppeakday <- ggplot(aes(x= fn, y = peak, color = Origin), data = Phen) + geom_point() + 
  geom_smooth(aes(linetype = Origin), method = "lm", se = T, alpha = 0.2) + 
  theme(
    panel.background = element_rect(fill = "white", colour = "grey",
                                    size = 2, linetype = "solid"),
    panel.grid.major = element_line(), 
    panel.grid.minor = element_line())+
  scale_colour_manual("Origin",   labels = c("Ex situ", "Wild"),
                      values = c("black","dodgerblue3"))+ 
  xlab("Number of flowers") + ylab("Peak flowering") + theme(legend.position = "none")
ppeakday

# Plot the difference between Origin
dpeakmean <- Phen

names(dpeakmean)[4] <- "emmean"
mpeak <- as.data.frame(emmeans(peak3, specs = "Origin"))
ppeakday2 <- ggplot(mpeak, aes(x = Origin, y = emmean)) + 
  geom_jitter(data = dpeakmean, aes(x = Origin, y = emmean, color = Origin), 
              size = 0.7, alpha = 0.6)+
  geom_point(aes(color = Origin), position = position_dodge(width =0.8), size = 3)+
  geom_errorbar(aes(x = Origin, ymin = emmean - 2*SE, ymax = emmean + 2*SE, color = Origin), position = position_dodge(width =0.8), 
                width = 0.05, size = 1.5)+
  theme(
    panel.background = element_rect(fill = "white", colour = "grey",
                                    size = 2, linetype = "solid"),
    panel.grid.major = element_line(), 
    panel.grid.minor = element_line())+
  scale_colour_manual("Origin",   labels = c("Ex situ", "Wild"),
                      values = c("black","dodgerblue3"))+ 
  xlab("Number of flowers") + ylab("Peak flowering") + theme(legend.position = "none")+
  scale_x_discrete (labels= c("ex situ", "wild"))+
  theme(axis.text.x = element_text(size = 12), 
        axis.title.x =element_text(size = 12),
        axis.title.y =element_text(size = 12),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size = 14))
  # scale_y_continuous(labels = scales::number_format(accuracy = 0.1))
ppeakday2
# Results above from the modeling of phenological data were reported in Appendix table 7

# arrrange plots 
# This figure below shows the relationship between the number of flowers with phenological time
# This figure is not reported in the main paper, because: 
# 1) This is not of our interest (major reason)
# 2) The plotting does not account for the interaction effect 
pfirst <- pfirst + theme(legend.position = "none") + theme(axis.title.x =element_blank())+
   theme(axis.text.x = element_text(size = 13), 
          axis.title.y =element_text(size = 14),
          axis.text.y = element_text(size = 13))
pduration <- pduration + theme(legend.position = "none") + theme(axis.title.x =element_blank()) + 
  theme(axis.text.x = element_text(size = 13), 
        axis.title.y =element_text(size = 14),
        axis.text.y = element_text(size = 13))
ppeakday <- ppeakday + theme(legend.position = "none") + theme(axis.title.x =element_blank())+
  theme(axis.text.x = element_text(size = 13), 
        axis.title.y =element_text(size = 14),
        axis.text.y = element_text(size = 13))

pphen <- ggarrange(pfirst, pduration, ppeakday,
               labels = c("A", "B", "C"),
               ncol = 3, nrow = 1) 
pphen1 <-annotate_figure (pphen, bottom = text_grob("Number of flowers per plant", size = 14))
pphen1


# arrange the plot -- appendix figure 6
pfirst2 <- pfirst2 +  theme(legend.position = "none") + 
  theme(axis.title.x =element_blank())+
  theme(axis.title.y =element_text(size = 15),
        axis.text.y = element_text(size = 10)) 
pduration2 <- pduration2 + theme(legend.position = "none")+ theme(axis.title.x =element_blank())+
  theme(axis.title.y =element_text(size = 15),
        axis.text.y = element_text(size = 10)) 
ppeakday2 <- ppeakday2 + theme(legend.position = "none")+ theme(axis.title.x =element_blank())+
  theme(axis.title.y =element_text(size = 15),
        axis.text.y = element_text(size = 10)) 

pphen2 <- ggarrange(pfirst2, pduration2, ppeakday2, 
               labels = c("A", "B", "C"), 
               ncol = 3, nrow = 1) 
pphen2
ggsave("pphen.tiff", dpi = 700)

mean(Phen$duration)# On average 94 days 
mean(Phen$first) # On average on the 179th day of the year 
mean(Phen$peak) # On average on the 202th day of the year


################## The number of flowers ##################

# dfn is the dataset with number of flowers per plant
dfn <- merge(dfn, Tag, by.x = 'ID', by.y = 'Tag..') # Add accession and origin information 

# Use posisson distribution
fn_p <- glm(fn ~ Origin, family = "poisson", dfn)
fn_p1 <- glmer(fn ~ Origin + (1|Accession..), family = "poisson", dfn)
fn_p2 <- glm(fn ~ Origin + Accession.., family = "poisson", dfn)
AIC(fn_p,fn_p1)
AIC(fn_p1,fn_p2)
AIC(fn_p2,fn_p)

# fn_p1 is significantly better
tab_model(fn_p1) # ICC is 0.94, meaning that something is wrong. 
# Then, turn to fn_p2, the second best model, which use accession as a 
# fixed factor and thereby can account for the high variation of ICC. 
tab_model(fn_p2)
summary(fn_p2) # the residual and degree of freedom -- overdispersion

# To account for the overdispersion -- change to negative binomial (NB) distribution
fn_nb <- glm.nb(fn ~ Origin,  dfn)
fn_nb1 <- glmer.nb(fn ~ Origin + (1|Accession..), dfn)
fn_nb2 <- glm.nb(fn ~ Origin + Accession.., dfn)

AIC(fn_nb,fn_new)
AIC(fn_nb1,fn_new1)
# fn_nb is much better -- chaning to NB distribution highly improves the likelihood
AIC(fn_nb1,fn_nb) #fn_nb1 is better
AIC(fn_nb1,fn_nb2) # fn_nb1 is better

tab_model(fn_nb1)
summary(fn_nb1) # better -- Thus, we changed "possion" to negative binomial 
getME(fn_nb1,"glmer.nb.theta") 
# Report: theta =  5.41





################## Floral traits ##################
# Function to calculate the coefficient of variance -- This will be used for the below section
co.var <- function(x)
  (
    sd(x)/mean(x)
  )

#### Hypanthium length 

# remove NA
d.hl <- df[!is.na(df$Hypanthium_length),] # create a new dataset d.hl removing NA values from hl columns
hist(d.hl$Hypanthium_length)
names(d.hl)[4] <- "hl"

# three datasets to check whether the shift in sampling method is a problem -- Appendix B 
names(d.hl)
# In column "type", it marks the three different sampling method -- Appendix B 
unique(df$type) 
# "" -- rows that only records the number of flowers without floral trait measurements
# "weekly": the dataset with no ambiguity, we sampled all flowersing plants daily -- appendix B 
# Noted that the weekly dataset is not sampled weekly. We sampled once to twice a week. We used this name for simplicity. 
# "random": the dataset sampled half of the flowering plants randomly after Aug 23 
# "non-random": the dataset sampled half of the flowering plants based on position before Aug 23 

d.hl.weekly = d.hl[d.hl$type == "weekly",] 
d.hl.before23 = d.hl[d.hl$type == "non-random",] 
d.hl.after23 = d.hl[d.hl$type == "random",] 

## Modeling -- Appendix B 

# All data including weekly, randon, and non-random
mhl <- lmer(hl ~ Origin + (1|Days) + (1|Tag) + (1|Accession..), data = d.hl)
Anova(mhl)  # P = 0.0006885 and Chisq = 11.52 and effect estimates = -0.69 for NM 
table_hl<- tab_model(mhl, show.ci = 0.98, show.df = T) 
table_hl
length(d.hl$hl) # 2235 measurements 

# Only weekly data
mhl.weekly <- lmer(hl ~ Origin+ (1|Days) + (1|Tag) + (1|Accession..), data = d.hl.weekly)
Anova(mhl.weekly)  # P = 0.0001628 and Chisq = 14.22 and effect estimates = -0.73 for NM 
table_hl.weekly<- tab_model(mhl.weekly, show.ci = 0.98, show.df = T) 
table_hl.weekly # Table 3
aggregate(hl ~ Origin, FUN = length, d.hl.weekly) # Table 3 
length(d.hl.weekly$hl)
mean(d.hl.weekly$hl)

# Only before 23 data
mhl.before23 <- lmer(hl ~ Origin + (1|Days) + (1|Tag) + (1|Accession..), data = d.hl.before23)
Anova(mhl.before23)  # P = 0.0001461 and Chisq = 14.42 and effect estimates = -0.73 for NM 
table_hl.before23<- tab_model(mhl.before23, show.ci = 0.98, show.df = T) 
table_hl.before23
length(d.hl.before23$hl)

mhl.after23 <- lmer(hl ~ Origin + (1|Days) + (1|Tag) + (1|Accession..), data = d.hl.after23)
Anova(mhl.after23)  # P = 0.3413 and Chisq = 0.9056 and effect estimates = -0.31 for NM 
table_hl.after23<- tab_model(mhl.after23, show.ci = 0.98, show.df = T) 
table_hl.after23 
length(d.hl.after23$hl)
# This is different from our all data analysis 
# Thus, the paper uses weekly dataset only because it was sampled unambiguously 

# Let's check if there is any difference between sampling method 
plot_model(mhl.weekly, type = "diag")[1]  
plot_model(mhl.weekly, type = "diag")[2] 
plot_model(mhl.weekly, type = "diag")[3]  
plot_model(mhl.weekly, type = "diag")[4] 
plot_model(mhl.weekly, type = "re", terms = "Origin", show.data = T)

# Bartlett test -- Variance 
# assign d.hl (all data dataset) to weekly dataset, because we uses the weekly dataset to avoid sampling bias
d.hl <- d.hl.weekly
length(d.hl$hl) # double check if the d.hl dataset is equal to the d.hl.weekly dataset now 

# If data are normally distributed? 
# Across the whole flowering time
p1 = ggplot(d.hl, aes(x=hl)) + 
  geom_histogram(binwidth=.25, colour="black", fill="white")
p1 = ggplot(d.hl, aes(x=hl)) +
  geom_histogram(aes(y = ..density..), binwidth=.25, colour="black", fill="white") +
  stat_function(fun = dnorm, lwd = 2, col = 'red', 
                args = list(mean = mean(d.hl$hl), sd = sd(d.hl$hl)))
p1 # normally distributed

# Compare the variance 
ftest.hl <- bartlett.test(hl ~ Origin, data = d.hl) # Table 4 
ftest.hl
aggregate(hl ~ Origin, FUN = length, d.hl)
aggregate(hl ~ Origin, FUN = var, d.hl) # Calculate variance 
aggregate(hl ~ Origin, FUN = co.var, d.hl) # Calculate coefficient of variance 

# Compare the variance in hl between three time period
# Early season
d.hl1 <- d.hl[d.hl$time == "early",]
phl1 = ggplot(d.hl1, aes(x=hl)) + 
  geom_histogram(binwidth=.25, colour="black", fill="white")
phl1 = ggplot(d.hl1, aes(x=hl)) +
  geom_histogram(aes(y = ..density..), binwidth=.25, colour="black", fill="white") +
  stat_function(fun = dnorm, lwd = 2, col = 'red', 
                args = list(mean = mean(d.hl1$hl), sd = sd(d.hl1$hl)))
phl1 # normally distributed
ftest.hl.1 <- bartlett.test(hl ~ Origin, data = d.hl1)
ftest.hl.1
aggregate(hl ~ Origin, FUN = length, d.hl1)
aggregate(hl ~ Origin, FUN = co.var, d.hl1)

# Middle season 
d.hl2 <- d.hl[d.hl$time == "middle",]
phl2 = ggplot(d.hl2, aes(x=hl)) + 
  geom_histogram(binwidth=.25, colour="black", fill="white")
phl2 = ggplot(d.hl2, aes(x=hl)) +
  geom_histogram(aes(y = ..density..), binwidth=.25, colour="black", fill="white") +
  stat_function(fun = dnorm, lwd = 2, col = 'red', 
                args = list(mean = mean(d.hl2$hl), sd = sd(d.hl2$hl)))
phl2 # normally distributed

ftest.hl.2 <- bartlett.test(hl ~ Origin, data = d.hl2)
ftest.hl.2
aggregate(hl ~ Origin, FUN = length, d.hl2)
aggregate(hl ~ Origin, FUN = co.var, d.hl2)

# Late season 
d.hl3 <- d.hl[d.hl$time == "late",]
phl3 = ggplot(d.hl3, aes(x=hl)) + 
  geom_histogram(binwidth=.25, colour="black", fill="white")
phl3 = ggplot(d.hl3, aes(x=hl)) +
  geom_histogram(aes(y = ..density..), binwidth=.25, colour="black", fill="white") +
  stat_function(fun = dnorm, lwd = 2, col = 'red', 
                args = list(mean = mean(d.hl3$hl), sd = sd(d.hl3$hl)))
phl3#  normally distributed

ftest.hl.3 <- bartlett.test(hl ~ Origin, data = d.hl3)
ftest.hl.3
aggregate(hl ~ Origin, FUN = length, d.hl3)
aggregate(hl ~ Origin, FUN = co.var, d.hl3)

# Plotting with three different time
names(df)[4] <- "hl"
dhlmean <- d.hl.weekly # dhlmean will be used to plot the raw data 
names(dhlmean)[4] <- "emmean" # Change hl to emmean although it is not estimated mean, just easier for coding 
mhl2 <- as.data.frame(emmeans(mhl.weekly, specs = "Origin")) # Get the estimated mean from the model 
print(mhl2)
phl <- ggplot(mhl2, aes(x = Origin, y = emmean), colour = Origin) + 
  geom_boxplot(data = dhlmean, aes(x = Origin, y = emmean, color = Origin), size = 0.7)+ # box plot show the raw data 
  geom_point(aes(colour = Origin), position = position_dodge(width =0.8), size = 3)+ 
  # points show the estimated mean, error bar too small so not shown in the figure
  theme(
    panel.background = element_rect(fill = "white", colour = "grey",
                                    size = 2, linetype = "solid"),
    panel.grid.major = element_line(), 
    panel.grid.minor = element_line()) + xlab("Time") + ylab("Floral tube length (cm)") + 
  theme(legend.position = "top",
        legend.text = element_text(size=10))+
  theme(legend.key = element_rect(fill = alpha("white", 0.0)))+
  theme(axis.text.x = element_text(size = 12), 
        axis.title.x =element_text(size = 12),
        axis.title.y =element_text(size = 12),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size = 14)) +
  theme(legend.text=element_text(size=12), legend.title=element_text(size=12)) +
  scale_x_discrete (labels= c("ex situ", "wild"))+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1))+
  scale_colour_manual("Origin",   labels = c("Ex situ", "Wild"),
                      values = c("black","dodgerblue3"))
phl # figure 5


## Floral tubes 
names(df)[5] <- "ft"
d.ft <- df[!is.na(df$ft),] # remove NA

hist(d.ft$ft) # floral tube typically less than 10, a wrong value is detected
d.ft <- d.ft[d.ft$ft <= 10,]
range(d.ft$ft) # reasonable range 

# To test whether the sampling method has a problem: 
d.ft.weekly <- d.ft[d.ft$type == "weekly", ]
d.ft.before23 <- d.ft[d.ft$type == "non-random", ]
d.ft.after23 <- d.ft[d.ft$type == "random", ]

### Modeling 
# modeling 
m.ft <- lmer(ft ~ Origin + (1|Days) + (1|Tag) + (1|Accession..), data = d.ft) 
Anova(m.ft) # Appendix B 
tab_model(m.ft, show.ci = 0.98) # No significance 
length(d.ft$ft)

# Check whether sampling method has any effect 
m.ft.weekly <- lmer(ft ~ Origin + (1|Days) + (1|Tag) + (1|Accession..), data = d.ft.weekly) 
Anova(m.ft.weekly) # p = 0.76
tab_model(m.ft.weekly, show.ci = 0.98) # effect estimates = -0.02 and No significance
length(d.ft.weekly$ft)

m.ft.before23 <- lmer(ft ~ Origin + (1|Days) + (1|Tag) + (1|Accession..), data = d.ft.before23) 
Anova(m.ft.before23) # p = 0.86 
tab_model(m.ft.before23, show.ci = 0.98) # effect estimates = -0.01 and No significance
length(d.ft.before23$ft)

m.ft.after23 <- lmer(ft ~ Origin + (1|Days) + (1|Tag) + (1|Accession..), data = d.ft.after23) 
Anova(m.ft.after23) # p = 0.87 
tab_model(m.ft.after23, show.ci = 0.98) # effect estimate = -0.01 and no significance 
length(d.ft.after23$ft)
# All those check show that for ft measurements, no effect in sampling method is observed. 
# However, after checking all floral traits, to keep consistency and avoid sampling bias, 
# we used weekly dataset (this needs to be seen after finishing analyzing alla traits, 
# especially check hypanthium length, style length and nectar volumn) 

# Thus, assign d.ft to d.ft.weekly since we used d.ft.weekly as the main dataset (Appendix B)
d.ft <- d.ft.weekly
range(d.ft$ft) # double check if the dataset is correctly assigned 

# Bartlett test 
# Are data normally distributed? 
p2 = ggplot(d.ft, aes(x=ft)) + 
  geom_histogram(binwidth=.25, colour="black", fill="white")
p2 = ggplot(d.ft, aes(x=ft)) +
  geom_histogram(aes(y = ..density..), binwidth=.25, colour="black", fill="white") +
  stat_function(fun = dnorm, lwd = 2, col = 'red', 
                args = list(mean = mean(d.ft$ft), sd = sd(d.ft$ft)))
p2 # slightly skewed
ftest.ft <- bartlett.test(ft ~ Origin, data = d.ft)
ftest.ft # Table 4 
aggregate(ft ~ Origin, FUN = length, d.ft) 
aggregate(ft ~ Origin, FUN = co.var, d.ft) # Table 4 

# EARLY
d.ft1 <- d.ft[d.ft$time == "early",]
pft1 = ggplot(d.ft1, aes(x=ft)) + 
  geom_histogram(binwidth=.25, colour="black", fill="white")
pft1 = ggplot(d.ft1, aes(x=ft)) +
  geom_histogram(aes(y = ..density..), binwidth=.25, colour="black", fill="white") +
  stat_function(fun = dnorm, lwd = 2, col = 'red', 
                args = list(mean = mean(d.ft1$ft), sd = sd(d.ft1$ft)))
pft1 # slightly skewed but nearly normally distributed

ftest.ft.1 <- bartlett.test(ft ~ Origin, data = d.ft1)
ftest.ft.1
aggregate(ft ~ Origin, FUN = length, d.ft1)
aggregate(ft ~ Origin, FUN = co.var, d.ft1)

# MIDDLE
d.ft2 <- d.ft[d.ft$time == "middle",]
pft2 = ggplot(d.ft2, aes(x=ft)) + 
  geom_histogram(binwidth=.25, colour="black", fill="white")
pft2 = ggplot(d.ft2, aes(x=ft)) +
  geom_histogram(aes(y = ..density..), binwidth=.25, colour="black", fill="white") +
  stat_function(fun = dnorm, lwd = 2, col = 'red', 
                args = list(mean = mean(d.ft1$ft), sd = sd(d.ft1$ft)))
pft2 # normally distributed
ftest.ft.2 <- bartlett.test(ft ~ Origin, data = d.ft2)
ftest.ft.2
aggregate(ft ~ Origin, FUN = length, d.ft2)
aggregate(ft ~ Origin, FUN = co.var, d.ft2)

# LATE
d.ft3 <- d.ft[d.ft$time == "late",]
pft3 = ggplot(d.ft3, aes(x=ft)) + 
  geom_histogram(binwidth=.25, colour="black", fill="white")
pft3 = ggplot(d.ft3, aes(x=ft)) +
  geom_histogram(aes(y = ..density..), binwidth=.25, colour="black", fill="white") +
  stat_function(fun = dnorm, lwd = 2, col = 'red', 
                args = list(mean = mean(d.ft3$ft), sd = sd(d.ft3$ft)))
pft3 # skewed but acceptable
ftest.ft.3 <- bartlett.test(ft ~ Origin, data = d.ft3)
ftest.ft.3 
aggregate(ft ~ Origin, FUN = length, d.ft3)
aggregate(ft ~ Origin, FUN = co.var, d.ft3)

# checking model
plot_model(m.ft.weekly, type = "diag")[1]  
plot_model(m.ft.weekly, type = "diag")[2] 
plot_model(m.ft.weekly, type = "diag")[3]  
plot_model(m.ft.weekly, type = "diag")[4] # Not good but acceptable 

plot_model(m.ft.weekly, type = "est")
plot_model(m.ft.weekly, type = "pred", terms = "Origin", show.data = T)
plot_model(m.ft.weekly, type = "pred", terms = "Origin", show.data = T, jitter = T)

# Plot: 
dftmean <- d.ft.weekly # have a new dataset with only non-NA ft measurement 
names(dftmean)[5] <- "emmean" # Raw data but change name for the easiness of coding 
mft2 <- as.data.frame(emmeans(m.ft.weekly, specs = "Origin"))
mft2
pft <- ggplot(mft2, aes(x = Origin, y = emmean, colour = Origin)) + 
  geom_boxplot(data = dftmean, aes(x = Origin, y = emmean, colour = Origin), size = 0.7)+ 
  geom_point(aes(colour = Origin), position = position_dodge(width =0.8), size = 3)+
  scale_colour_manual("Origin",   labels = c("Ex situ", "Wild"),
                      values = c("black","dodgerblue3"))+
  theme(
    panel.background = element_rect(fill = "white", colour = "grey",
                                    size = 2, linetype = "solid"),
    panel.grid.major = element_line(), 
    panel.grid.minor = element_line()) + xlab("Time") + ylab("Floral flare width (mm)") + 
  theme(legend.position = "top",
        legend.text = element_text(size=10))+
  theme(legend.key = element_rect(fill = alpha("white", 0.0)))+
  theme(axis.text.x = element_text(size = 12), 
        axis.title.x =element_text(size = 12),
        axis.title.y =element_text(size = 12),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size = 14)) +
  theme(legend.text=element_text(size=12), legend.title=element_text(size=12)) +
  scale_x_discrete (labels= c("ex situ", "wild"))+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1))
pft


### Style length

names(df)[6] <- "sl"
d.sl <- df[!is.na(df$sl),] # remove NA 
range(d.sl$sl)
d.sl <- d.sl[d.sl$sl >= 8.3,] # longer than the smallest floral tube length 
range(d.sl$sl)

# Build data frame to check the effect of sampling method: 
d.sl.weekly <- d.sl[d.sl$type == "weekly",]
d.sl.before23 <- d.sl[d.sl$type == "non-random",]
d.sl.after23 <- d.sl[d.sl$type == "random",]

## Modeling 

# Check if the sampling method has any bias 
m.sl <- lmer(sl ~ Origin + (1|Days) + (1|Tag) + (1|Accession..), data = d.sl) # All data 
Anova(m.sl) # p = 0.002696 ** 
tab_model(m.sl.2, show.ci = 0.98) # effect estimates = -0.78 for NM 
length(d.sl$sl)

# weekly
m.sl.weekly <- lmer(sl ~ Origin + (1|Days) + (1|Tag) + (1|Accession..), data = d.sl.weekly)
Anova(m.sl.weekly) # p = 0.001953**  # Appendix B 
tab_model(m.sl.weekly, show.ci = 0.98) # Effect estimate = -0.81 
length(d.sl.weekly$sl)
aggregate(sl ~ Origin, FUN = length, d.sl.weekly)

# Before 23
m.sl.before23 <- lmer(sl ~ Origin + (1|Days) + (1|Tag) + (1|Accession..), data = d.sl.before23)
Anova(m.sl.before23) #  0.0003379 *** # Appendix B 
tab_model(m.sl.before23, show.ci = 0.98)  # E.E = -0.81 
length(d.sl.before23$sl)

# after 23
m.sl.after23 <- lmer(sl ~ Origin + (1|Days) + (1|Tag) + (1|Accession..), data = d.sl.after23)
Anova(m.sl.after23) # p = 0.2097 
tab_model(m.sl.after23, show.ci = 0.98) # No significance
length(d.sl.after23$sl)
# It seems that the after-23 dataset -- there is a difference in sampling method 
# Use the weekly dataset 

# Checking model 
plot_model(m.sl.weekly, type = "diag")[1]  
plot_model(m.sl.weekly, type = "diag")[2] 
plot_model(m.sl.weekly, type = "diag")[3]  
plot_model(m.sl.weekly, type = "diag")[4] 

plot_model(m.sl.weekly, type = "est")
plot_model(m.sl.weekly, type = "pred", terms = "Origin", show.data = T)
plot_model(m.sl.weekly, type = "pred", terms = "Origin", show.data = T, jitter = T)
plot_model(m.sl.weekly, type = "re", terms = "Origin", show.data = T)

# Use the weekly dataset, so assign d.sl as the d.sl.weekly 
d.sl <- d.sl.weekly 
length(d.sl$sl)

# Are data normally distributed? 
p3 = ggplot(d.sl, aes(x=sl)) + 
  geom_histogram(binwidth=.25, colour="black", fill="white")
p3 = ggplot(d.sl, aes(x=sl)) +
  geom_histogram(aes(y = ..density..), binwidth=.25, colour="black", fill="white") +
  stat_function(fun = dnorm, lwd = 2, col = 'red', 
                args = list(mean = mean(d.sl$sl), sd = sd(d.sl$sl)))
p3 # slightly skewed but acceptable
ftest.sl <- bartlett.test(sl ~ Origin, data = d.sl)
ftest.sl
aggregate(sl ~ Origin, FUN = length, d.sl)
aggregate(sl ~ Origin, FUN = co.var, d.sl)

# Early
d.sl1 <- d.sl[d.sl$time == "early",]
psl1 = ggplot(d.sl1, aes(x=sl)) + 
  geom_histogram(binwidth=.25, colour="black", fill="white")
psl1 = ggplot(d.sl1, aes(x=sl)) +
  geom_histogram(aes(y = ..density..), binwidth=.25, colour="black", fill="white") +
  stat_function(fun = dnorm, lwd = 2, col = 'red', 
                args = list(mean = mean(d.sl1$sl), sd = sd(d.sl1$sl)))
psl1 # skewed
ftest.sl.1 <-bartlett.test(sl ~ Origin, data = d.sl1)
ftest.sl.1
leveneTest(sl ~ Origin, data = d.sl1) # p = 0.93 no significant difference with the bartlett test 
aggregate(sl ~ Origin, FUN = length, d.sl1)
aggregate(sl ~ Origin, FUN = co.var, d.sl1)

# Middle
d.sl2 <- d.sl[d.sl$time == "middle",]
psl2 = ggplot(d.sl2, aes(x=sl)) + 
  geom_histogram(binwidth=.25, colour="black", fill="white")
psl2 = ggplot(d.sl2, aes(x=sl)) +
  geom_histogram(aes(y = ..density..), binwidth=.25, colour="black", fill="white") +
  stat_function(fun = dnorm, lwd = 2, col = 'red', 
                args = list(mean = mean(d.sl2$sl), sd = sd(d.sl2$sl)))
psl2 # skewed
ftest.sl.2 <-bartlett.test(sl ~ Origin, data = d.sl2)
ftest.sl.2
aggregate(sl ~ Origin, FUN = length, d.sl2)
aggregate(sl ~ Origin, FUN = co.var, d.sl2)

# Late
d.sl3 <- d.sl[d.sl$time == "late",]
psl3 = ggplot(d.sl3, aes(x=sl)) + 
  geom_histogram(binwidth=.25, colour="black", fill="white")
psl3 = ggplot(d.sl2, aes(x=sl)) +
  geom_histogram(aes(y = ..density..), binwidth=.25, colour="black", fill="white") +
  stat_function(fun = dnorm, lwd = 2, col = 'red', 
                args = list(mean = mean(d.sl3$sl), sd = sd(d.sl3$sl)))
psl3 #  skewed but acceptable
ftest.sl.3 <-bartlett.test(sl ~ Origin, data = d.sl3)
ftest.sl.3
aggregate(sl ~ Origin, FUN = length, d.sl3)
aggregate(sl ~ Origin, FUN = co.var, d.sl3)

# Plotting 
dslmean <- d.sl.weekly # used to plot the raw data 
names(dslmean)[6] <- "emmean" # raw data not estimated mean but changed for easiness 
msl2 <- as.data.frame(emmeans(m.sl.weekly, specs = "Origin"))
print(msl2)
psl <- ggplot(msl2, aes(x = Origin, y = emmean), colour = Origin) + 
  geom_boxplot(data = dslmean, aes(x = Origin, y = emmean, colour = Origin))+
  geom_point(aes(color = Origin), position = position_dodge(width =0.8), size = 3)+
  scale_colour_manual("Origin",   labels = c("Ex situ", "Wild"),
                      values = c("black","dodgerblue3"))+
  theme(
    panel.background = element_rect(fill = "white", colour = "grey",
                                    size = 2, linetype = "solid"),
    panel.grid.major = element_line(), 
    panel.grid.minor = element_line()) + xlab("Time") + ylab("Style length (cm)") + 
  theme(legend.position = "top",
        legend.text = element_text(size=10))+
  theme(legend.key = element_rect(fill = alpha("white", 0.0)))+
  theme(axis.text.x = element_text(size = 12), 
        axis.title.x =element_text(size = 12),
        axis.title.y =element_text(size = 12),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size = 14)) +
  theme(legend.text=element_text(size=12), legend.title=element_text(size=12))+
  scale_x_discrete (labels= c("ex situ", "wild"))
psl


## Style extension: style - hypanthium

# Style extension is added because it could depend how pollen deposition to/from the pollinators  
d.sh <- df # create a new dataset 
d.sh <- d.sh[!is.na(d.sh$hl),] # remove hl NA
d.sh <- d.sh[!is.na(d.sh$sl),] # remove sl NA 
d.sh$sh <- d.sh$sl-d.sh$hl # -- calculate the extension that style is longer than hypanthium 
range(d.sh$sh) # Negative values are clearly wrong 
d.sh <- d.sh[d.sh$sh > 0,]
range(d.sh$sh) # double check if there is any negative value 

## Modeling 
# Test whether the shift in sampling method would cause a problem: 
d.sh.weekly <- d.sh[d.sh$type == "weekly",]
d.sh.before23 <- d.sh[d.sh$type == "non-random",]
d.sh.after23 <- d.sh[d.sh$type == "random",]

# All data  
m.sh <- lmer(sh ~Origin + (1|Days) +  + (1|Tag) + (1|Accession..), data = d.sh)
Anova(m.sh) # p = 0.39 
tab_model(m.sh, show.ci = 0.98, show.df = T) 
length(d.sh$sh)

# weekly data 
m.sh.weekly <- lmer(sh ~ Origin + (1|Days) + (1|Tag) + (1|Accession..), data = d.sh.weekly)
Anova(m.sh.weekly) # p = 0.45
tab_model(m.sh.weekly, show.ci = 0.98, show.df = T)  # Non-significant 
length(d.sh.weekly$sh)

# Before 23
m.sh.before23 <- lmer(sh ~ Origin + (1|Days) + (1|Tag) + (1|Accession..), data = d.sh.before23)
Anova(m.sh.before23) # p = 0.48 
tab_model(m.sh.before23, show.ci = 0.98, show.df = T) 
length(d.sh.before23$sh)

# after 23
m.sh.after23 <- lmer(sh ~ Origin + (1|Days) + (1|Tag) + (1|Accession..), data = d.sh.after23)
Anova(m.sh.after23) # p = 0.085 
tab_model(m.sh.after23, show.ci = 0.98, show.df = T) # Marginal significance p = 0.085 
length(d.sh.after23$sh)

# checking 
plot_model(m.sh.weekly, type = "diag")[1]  
plot_model(m.sh.weekly, type = "diag")[2] 
plot_model(m.sh.weekly, type = "diag")[3]  
plot_model(m.sh.weekly, type = "diag")[4] 

plot_model(m.sh, type = "est")
plot_model(m.sh, type = "pred", terms = "Origin", show.data = T, jitter = T, ci.style = "box")
plot_model(m.sh, type = "re", terms = "Origin", show.data = T) 

## Bartlett test 
# Decided to use the weekly data, so assign d.sh as d.sh.weekly 
d.sh <- d.sh.weekly
length(d.sh$sh) # double checked if the dataset is correctly assigned 

# If data are normally distributed? 
psh = ggplot(d.sh, aes(x=sh)) + 
  geom_histogram(binwidth=.25, colour="black", fill="white")
psh = ggplot(d.sh, aes(x=sh)) +
  geom_histogram(aes(y = ..density..), binwidth=.25, colour="black", fill="white") +
  stat_function(fun = dnorm, lwd = 2, col = 'red', 
                args = list(mean = mean(d.sh$sh), sd = sd(d.sh$sh)))
psh # normally distributed 
ftest.sh <- bartlett.test(sh ~ Origin, data = d.sh) # Compare the variance  
ftest.sh
aggregate(sh ~ Origin, FUN = length, d.sh)
aggregate(sh ~ Origin, FUN = co.var, d.sh)

# Compare the variance in hl between three time period
# early 
d.sh1 <- d.sh[d.sh$time == "early",]
psh1 = ggplot(d.sh1, aes(x=sh)) + 
  geom_histogram(binwidth=.25, colour="black", fill="white")
psh1 = ggplot(d.sh1, aes(x=sh)) +
  geom_histogram(aes(y = ..density..), binwidth=.25, colour="black", fill="white") +
  stat_function(fun = dnorm, lwd = 2, col = 'red', 
                args = list(mean = mean(d.sh1$sh), sd = sd(d.sh1$sh)))
psh1 # skwed but acceptable 
ftest.sh.1 <- bartlett.test(sh ~ Origin, data = d.sh1)
ftest.sh.1
leveneTest(sh ~ Origin, data = d.sh1) # match the bartlett test 
aggregate(sh ~ Origin, FUN = length, d.sh1)
aggregate(sh ~ Origin, FUN = co.var, d.sh1)

# middle 
d.sh2 <- d.sh[d.sh$time == "middle",]
psh2 = ggplot(d.sh2, aes(x=sh)) + 
  geom_histogram(binwidth=.25, colour="black", fill="white")
psh2 = ggplot(d.sh2, aes(x=sh)) +
  geom_histogram(aes(y = ..density..), binwidth=.25, colour="black", fill="white") +
  stat_function(fun = dnorm, lwd = 2, col = 'red', 
                args = list(mean = mean(d.sh2$sh), sd = sd(d.sh2$sh)))
psh2 # slightly skwed but still normally distributed
ftest.sh.2 <- bartlett.test(sh ~ Origin, data = d.sh2)
ftest.sh.2
aggregate(sh ~ Origin, FUN = length, d.sh2)
aggregate(sh ~ Origin, FUN = co.var, d.sh2)

# late
d.sh3 <- d.sh[d.sh$time == "late",]
psh3 = ggplot(d.sh3, aes(x=sh)) + 
  geom_histogram(binwidth=.25, colour="black", fill="white")
psh3 = ggplot(d.sh3, aes(x=sh)) +
  geom_histogram(aes(y = ..density..), binwidth=.25, colour="black", fill="white") +
  stat_function(fun = dnorm, lwd = 2, col = 'red', 
                args = list(mean = mean(d.sh3$sh), sd = sd(d.sh3$sh)))
psh3# skewed but acceptable 
ftest.sh.3 <- bartlett.test(sh ~ Origin, data = d.sh3)
ftest.sh.3
leveneTest(sh ~ Origin, data = d.sh3) # match with the Bartlett test 
aggregate(sh ~ Origin, FUN = length, d.sh3)
aggregate(sh ~ Origin, FUN = co.var, d.sh3)

# Plotting: 
dshmean <- d.sh.weekly
names(dshmean)[30] <- "emmean" # plot raw data although it is not emmean
msh2 <- as.data.frame(emmeans(m.sh.weekly, specs = "Origin"))
print(msh2)
psh <- ggplot(msh2, aes(x = Origin, y = emmean)) + 
  geom_boxplot(data = dshmean, aes(x = Origin, y = emmean, colour = Origin))+
  geom_point(aes(colour = Origin), position = position_dodge(width =0.8), size = 3)+
  scale_colour_manual("Origin",   labels = c("Ex situ", "Wild"),
                      values = c("black","dodgerblue3"))+
  theme(
    panel.background = element_rect(fill = "white", colour = "grey",
                                    size = 2, linetype = "solid"),
    panel.grid.major = element_line(), 
    panel.grid.minor = element_line()) + xlab("Time") + ylab("Style extension (cm)") + 
  theme(legend.position = "top",
        legend.text = element_text(size=10))+
  theme(legend.key = element_rect(fill = alpha("white", 0.0)))+
  theme(axis.text.x = element_text(size = 12), 
        axis.title.x =element_text(size = 12),
        axis.title.y =element_text(size = 12),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size = 14)) +
  theme(legend.text=element_text(size=12), legend.title=element_text(size=12)) + 
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1))
psh


## Sugar Content 
names(df)[7] <- "sc"
d.sc <- df
d.sc <- d.sc[!is.na(d.sc$sc),] # remove NAs

## Modeling 
# Test whether the shift in sampling method would make a difference: 
d.sc.weekly <- d.sc[d.sc$type =="weekly",]
d.sc.before23 <- d.sc[d.sc$type =="non-random",]
d.sc.after23 <- d.sc[d.sc$type =="random",]

# All data 
m.sc <- lmer(sc ~ Origin + (1|Days) + (1|Tag) + (1|Accession..), data = d.sc)
Anova(m.sc) # p = 0.15 
tab_model(m.sc, show.ci = 0.98) # EE = 0.52 for NM, non-significant 

# weekly data 
m.sc.weekly <- lmer(sc ~ Origin + (1|Days) + (1|Tag) + (1|Accession..), data = d.sc.weekly)
Anova(m.sc.weekly) # p = 0.16
tab_model(m.sc.weekly, show.ci = 0.98) # E.E. = 0.48 
aggregate(sc ~ Origin, FUN = length, d.sc.weekly)
length(d.sc.weekly$sc)

# before 23
m.sc.before23 <- lmer(sc ~ Origin + (1|Days) + (1|Tag) + (1|Accession..), data = d.sc.before23)
Anova(m.sc.before23) # p = 0.16
tab_model(m.sc.before23, show.ci = 0.98) # E.E. = 0.48
length(d.sc.before23$sc)

# after 23
m.sc.after23 <- lmer(sc ~ Origin + (1|Days) + (1|Tag) + (1|Accession..), data = d.sc.after23)
Anova(m.sc.after23) # p = 0.252
tab_model(m.sc.after23, show.ci = 0.98) # E.E = 0.48 
length(d.sc.after23$sc)
# The shifting in sampling method did not affect the result 
# But to keep the consistency, we use weekly dataset -- Appendix B 

# checking model 
plot_model(m.sc.weekly, type = "diag")[1]  
plot_model(m.sc.weekly, type = "diag")[2] 
plot_model(m.sc.weekly, type = "diag")[3]  
plot_model(m.sc.weekly, type = "diag")[4] 

plot_model(m.sc.weekly, type = "est")
plot_model(m.sc.weekly, type = "pred", terms = "Origin", show.data = T)
plot_model(m.sc.weekly, type = "pred", terms = "Origin", show.data = T, jitter = T)
plot_model(m.sc, type = "re", terms = "Origin", show.data = T)  

# Since we use the weekly dataset, we will change to d.sc as the d.sc.weekly 
d.sc <- d.sc.weekly 
length(d.sc$sc)

# Are data normally distributed? 
p4 = ggplot(d.sc, aes(x=sc)) + 
  geom_histogram(binwidth=.25, colour="black", fill="white")
p4 = ggplot(d.sc, aes(x=sc)) +
  geom_histogram(aes(y = ..density..), binwidth=.25, colour="black", fill="white") +
  stat_function(fun = dnorm, lwd = 2, col = 'red', 
                args = list(mean = mean(d.sc$sc), sd = sd(d.sc$sc)))
p4 #normally distributed 
Ftest.sc <- bartlett.test(sc ~ Origin, data = d.sc) 
Ftest.sc
aggregate(sc ~ Origin, FUN = length, d.sc)
aggregate(sc ~ Origin, FUN = co.var, d.sc)

# Early
d.sc1 <- d.sc[d.sc$time == "early",]
psc1 = ggplot(d.sc1, aes(x=sc)) + 
  geom_histogram(binwidth=.25, colour="black", fill="white")
psc1 = ggplot(d.sc1, aes(x=sc)) +
  geom_histogram(aes(y = ..density..), binwidth=.25, colour="black", fill="white") +
  stat_function(fun = dnorm, lwd = 2, col = 'red', 
                args = list(mean = mean(d.sc1$sc), sd = sd(d.sc1$sc)))
psc1 # acceptable bur skewed 
ftest.sc.1 <-bartlett.test(sc ~ Origin, data = d.sc1)
ftest.sc.1
aggregate(sc ~ Origin, FUN = length, d.sc1)
aggregate(sc ~ Origin, FUN = co.var, d.sc1)

# Middle
d.sc2 <- d.sc[d.sc$time == "middle",]
psc2 = ggplot(d.sc2, aes(x=sc)) + 
  geom_histogram(binwidth=.25, colour="black", fill="white")
psc2 = ggplot(d.sc2, aes(x=sc)) +
  geom_histogram(aes(y = ..density..), binwidth=.25, colour="black", fill="white") +
  stat_function(fun = dnorm, lwd = 2, col = 'red', 
                args = list(mean = mean(d.sc2$sc), sd = sd(d.sc2$sc)))
psc2 # acceptable 
ftest.sc.2 <-bartlett.test(sc ~ Origin, data = d.sc2)
ftest.sc.2
aggregate(sc ~ Origin, FUN = length, d.sc2)
aggregate(sc ~ Origin, FUN = co.var, d.sc2)

# Late
d.sc3 <- d.sc[d.sc$time == "late",]
psc3 = ggplot(d.sc3, aes(x=sc)) + 
  geom_histogram(binwidth=.25, colour="black", fill="white")
psc3 = ggplot(d.sc3, aes(x=sc)) +
  geom_histogram(aes(y = ..density..), binwidth=.25, colour="black", fill="white") +
  stat_function(fun = dnorm, lwd = 2, col = 'red', 
                args = list(mean = mean(d.sc3$sc), sd = sd(d.sc3$sc)))
psc3 # normally distributed 
ftest.sc.3 <-bartlett.test(sc ~ Origin, data = d.sc3)
ftest.sc.3
aggregate(sc ~ Origin, FUN = length, d.sc3)
aggregate(sc ~ Origin, FUN = co.var, d.sc3)

# plotting 
dscmean <- d.sc.weekly # plot raw data 
names(dscmean)[7] <- "emmean" # raw data but change the name for the easiness of coding 
msc2 <- as.data.frame(emmeans(m.sc.weekly, specs = "Origin"))
print(msc2)
psc <- ggplot(msc2, aes(x = Origin, y = emmean), colour = Origin) + 
  geom_boxplot(data = dscmean, aes(x = Origin, y = emmean, colour = Origin))+
  geom_point(aes( colour = Origin), position = position_dodge(width =0.8), size = 3)+
  scale_colour_manual("Origin",   labels = c("Ex situ", "Wild"),
                      values = c("black","dodgerblue3"))+
  theme(
    panel.background = element_rect(fill = "white", colour = "grey",
                                    size = 2, linetype = "solid"),
    panel.grid.major = element_line(), 
    panel.grid.minor = element_line()) + xlab("Time") + ylab("Nectar sugar content (%)") + 
  theme(legend.position = "top",
        legend.text = element_text(size=10))+
  theme(legend.key = element_rect(fill = alpha("white", 0.0)))+
  theme(axis.text.x = element_text(size = 12), 
        axis.title.x =element_text(size = 12),
        axis.title.y =element_text(size = 12),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size = 14)) +
  theme(legend.text=element_text(size=12), legend.title=element_text(size=12))+
  scale_x_discrete (labels= c("ex situ", "wild"))+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1))
psc


# Nectar height/volume 
names(df)[8] <- "nh"
df$nv <- df$nh*3.1415*0.25*0.25 # nectar volume 
# nh is the height, 0.5 is the inner diameter pi*r^2 = pi*0.25*0.25 
d.nh <- df[!is.na(df$nh),]

## Modeling 
# Check if the shift in sampling method would affect the result 
d.nh.weekly <- d.nh[d.nh$type == "weekly",]
d.nh.before23 <- d.nh[d.nh$type == "non-random",]
d.nh.after23 <- d.nh[d.nh$type == "random",]

# Modeling 
# All data 
m.nv <- lmer(nv ~ Origin + (1|Days) + (1|Tag) + (1|Accession..), data = d.nh)
tab_model(m.nv, show.ci = 0.98) # p = 0.015 
Anova(m.nv) # p = 0.0145 
length(d.nh$nv)

# weekly data 
m.nv.weekly <- lmer(nv ~ Origin + (1|Days) + (1|Tag) + (1|Accession..), data = d.nh.weekly)
tab_model(m.nv.weekly, show.ci = 0.98) # p = 0.05
Anova(m.nv.weekly) # p = 0.05 
length(d.nh.weekly$nv)
aggregate(nv ~ Origin, length, data = d.nh.weekly)
bartlett.test(nv ~ Origin, data = d.nh.weekly)

# before 23 
m.nv.before23 <- lmer(nv ~ Origin + (1|Days) + (1|Tag) + (1|Accession..), data = d.nh.before23)
tab_model(m.nv.before23, show.ci = 0.98) # p = 0.013
Anova(m.nv.before23) # p = 0.01226 * significant 
length(d.nh.before23$nv)

# after 23
m.nv.after23 <- lmer(nv ~ Origin + (1|Days) + (1|Tag) + (1|Accession..), data = d.nh.after23)
tab_model(m.nv.after23, show.ci = 0.98) # p = 0.121
Anova(m.nv.after23) # p = 0.119 non significance 
length(d.nh.after23$nv)

# Not consistent -- Appendix B 
# use weekly dataset to avoid sampling bias 

# Checking data 
plot_model(m.nv.weekly, type = "diag")[1]  
plot_model(m.nv.weekly, type = "diag")[2] 
plot_model(m.nv.weekly, type = "diag")[3]  
plot_model(m.nv.weekly, type = "diag")[4] 

plot_model(m.nv.weekly, type = "est")
plot_model(m.nv.weekly, type = "pred", terms = "Origin", show.data = T)
plot_model(m.nv.weekly, type = "pred", terms = "Origin", show.data = T, jitter = T)
plot_model(m.nv.weekly, type = "re", terms = "Origin", show.data = T) 

## Since we will change to use weekly data, so we assigned d.nh as d.nh.weekly 
d.nh <- d.nh.weekly 
length(d.nh.weekly$nv)

## Bartlett test 
# Are data normally distributed? 
p5 = ggplot(d.nh, aes(x=nv)) + 
  geom_histogram(binwidth=.25, colour="black", fill="white")
p5 = ggplot(d.nh, aes(x=nv)) +
  geom_histogram(aes(y = ..density..), binwidth=.25, colour="black", fill="white") +
  stat_function(fun = dnorm, lwd = 2, col = 'red', 
                args = list(mean = mean(d.nh$nv), sd = sd(d.nh$nv)))
p5 # skewed but  normally distributed -- Try Levene test 

# Tried to use squared and squared root but it cannot improve normality 
names(d.nh)
d.nh$sqnv <- (d.nh$nv)*(d.nh$nv)
p5_sq = ggplot(d.nh, aes(x=sqnv)) + 
  geom_histogram(binwidth=.25, colour="black", fill="white")
p5_sq = ggplot(d.nh, aes(x=sqnv)) +
  geom_histogram(aes(y = ..density..), binwidth=.25, colour="black", fill="white") +
  stat_function(fun = dnorm, lwd = 2, col = 'red', 
                args = list(mean = mean(d.nh$sqnv), sd = sd(d.nh$sqnv)))
p5_sq # doesn't look good 

# Tried to use square root but improves normality 
d.nh$sqrtnv <- sqrt(d.nh$nv)
p5_sqrt = ggplot(d.nh, aes(x=sqrtnv)) + 
  geom_histogram(binwidth=.25, colour="black", fill="white")
p5_sqrt = ggplot(d.nh, aes(x=sqrtnv)) +
  geom_histogram(aes(y = ..density..), binwidth=.25, colour="black", fill="white") +
  stat_function(fun = dnorm, lwd = 2, col = 'red', 
                args = list(mean = mean(d.nh$sqrtnv), sd = sd(d.nh$sqrtnv)))
p5_sqrt # looks good

# untransformed data 
Ftest.nh <- bartlett.test(nh ~ Origin, data = d.nh) 
Ftest.nh # p = 0.24 
leveneTest(nh ~ Origin, data = d.nh) # p = 0.52, warning sign 
aggregate(nv ~ Origin, length, data = d.nh)
aggregate(nv ~ Origin, co.var, data = d.nh)
# square rooted data -- better normality 
Ftest.nh.sqrt <- bartlett.test(sqrtnv ~ Origin, data = d.nh) 
Ftest.nh.sqrt # p = 0.17 the p-values match with the untransformed data 

d.nv1 <- d.nh[d.nh$time == "early",]
aggregate(nv ~ Origin, co.var, data = d.nv1)
aggregate(nv ~ Origin, FUN = length, d.nv1)

d.nv2 <- d.nh[d.nh$time == "middle",]
aggregate(nv ~ Origin, co.var, data = d.nv2)
aggregate(nv ~ Origin, FUN = length, d.nv2)

d.nv3 <- d.nh[d.nh$time == "late",]
aggregate(nv ~ Origin, co.var, data = d.nv3)
aggregate(nv ~ Origin, FUN = length, d.nv3)

# Early
pnv1 = ggplot(d.nv1, aes(x=nv)) + 
  geom_histogram(binwidth=.25, colour="black", fill="white")
pnv1 = ggplot(d.nv1, aes(x=nv)) +
  geom_histogram(aes(y = ..density..), binwidth=.25, colour="black", fill="white") +
  stat_function(fun = dnorm, lwd = 2, col = 'red', 
                args = list(mean = mean(d.nv1$nv), sd = sd(d.nv1$nv)))
pnv1 # skwed 
ftest.nv.1 <-bartlett.test(nv ~ Origin, data = d.nv1)
leveneTest(nv ~ Origin, data = d.nv1)
ftest.nv.1 # p = 0.65 
# levene test generates warning sign and has similar p with bartlette test 
ftest.nv.1_sqrt <-bartlett.test(sqrtnv ~ Origin, data = d.nv1)
ftest.nv.1_sqrt # p = 0.39 non-significance match with the untranformed data 

# Middle
pnv2 = ggplot(d.nv2, aes(x=nv)) + 
  geom_histogram(binwidth=.25, colour="black", fill="white")
pnv2 = ggplot(d.nv2, aes(x=nv)) +
  geom_histogram(aes(y = ..density..), binwidth=.25, colour="black", fill="white") +
  stat_function(fun = dnorm, lwd = 2, col = 'red', 
                args = list(mean = mean(d.nv2$nv), sd = sd(d.nv2$nv)))
pnv2 # skwed but acceptable 
ftest.nv.2 <-bartlett.test(nv ~ Origin, data = d.nv2)
leveneTest(nv ~ Origin, data = d.nv2)
ftest.nv.2 #levene test generates warning sign and has similar p with bartlette test 
ftest.nv.2_sqrt <-bartlett.test(sqrtnv ~ Origin, data = d.nv2)
ftest.nv.2_sqrt # p = 0.87 non-significance match with the untranformed data 

# late
pnv3 = ggplot(d.nv3, aes(x=nv)) + 
  geom_histogram(binwidth=.25, colour="black", fill="white")
pnv3 = ggplot(d.nv3, aes(x=nv)) +
  geom_histogram(aes(y = ..density..), binwidth=.25, colour="black", fill="white") +
  stat_function(fun = dnorm, lwd = 2, col = 'red', 
                args = list(mean = mean(d.nv3$nv), sd = sd(d.nv3$nv)))
pnv3 # skwed 
ftest.nv.3 <-bartlett.test(nv ~ Origin, data = d.nv3)
leveneTest(nv ~ Origin, data = d.nv3)
ftest.nv.3 # levene test generates warning sign and has similar p with bartlette test 
ftest.nv.3_sqrt <-bartlett.test(sqrtnv ~ Origin, data = d.nv3)
ftest.nv.3_sqrt # p = 0.28 non-significance match with the untranformed data 

# Although square rooted values of nectar volumes have better normality, 
# the result is matching with the untransformed data 
# thus, the paper reports the untransformed data 

# Plotting 
dnvmean <- d.nh.weekly
names(dnvmean)[30] <- "emmean"
mnv2 <- as.data.frame(emmeans(m.nv.weekly, specs = "Origin"))
print(mnv2)
pnv <- ggplot(mnv2, aes(x = Origin, y = emmean), colour = Origin) + 
  geom_boxplot(data = dnvmean, aes(x = Origin, y = emmean, colour = Origin))+
  geom_point(aes(color = Origin),position = position_dodge(width =0.8), size = 3)+
  scale_colour_manual("Origin",   labels = c("Ex situ", "Wild"),
                      values = c("black","dodgerblue3"))+
  theme(
    panel.background = element_rect(fill = "white", colour = "grey",
                                    size = 2, linetype = "solid"),
    panel.grid.major = element_line(), 
    panel.grid.minor = element_line()) + xlab("Time") + ylab("Nectar volume (uL)") + 

  theme(legend.position = "top",
        legend.text = element_text(size=10))+
  theme(legend.key = element_rect(fill = alpha("white", 0.0)))+
  theme(axis.text.x = element_text(size = 12), 
        axis.title.x =element_text(size = 12),
        axis.title.y =element_text(size = 12),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size = 14)) +
  theme(legend.text=element_text(size=12), legend.title=element_text(size=12))+  
  scale_x_discrete (labels= c("ex situ", "wild"))+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) 
pnv

# Arrange plotting 
phl <- phl +  theme(legend.position = "none") + 
  theme(axis.title.x =element_blank())+
  theme(axis.text.x =element_blank())+
  theme(axis.title.y =element_text(size = 15),
        axis.text.y = element_text(size = 10)) 
pft <- pft + theme(legend.position = "none") + theme(axis.title.x =element_blank())+
  theme(axis.text.x =element_blank())+
  theme(axis.title.y =element_text(size = 15),
        axis.text.y = element_text(size = 10)) 
psc <- psc + theme(legend.position = "none")+ theme(axis.title.x =element_blank())+
  theme(axis.text.x =element_blank())+
  theme(axis.title.y =element_text(size = 15),
        axis.text.y = element_text(size = 10)) 
psl <- psl + theme(legend.position = "none")+ theme(axis.title.x =element_blank())+
  theme(axis.text.x =element_blank())+
  theme(axis.title.y =element_text(size = 15),
        axis.text.y = element_text(size = 10)) 
pnv <- pnv + theme(legend.position = "none")+ theme(axis.title.x =element_blank())+
  theme(axis.text.x =element_blank())+
  theme(axis.title.y =element_text(size = 15),
        axis.text.y = element_text(size = 10)) 
psh <- psh + theme(legend.position = "none")+ theme(axis.title.x =element_blank())+
  theme(axis.text.x =element_blank())+
  theme(axis.title.y =element_text(size = 15),
        axis.text.y = element_text(size = 10)) 

p <- ggarrange(phl, psl, pft, psc, pnv, psh,
               labels = c("A", "B", "C","D","E", "F"), 
               ncol = 2, nrow = 3) 
p
ggsave("p.tiff", dpi = 1200)




#### MDS NON-metric multidimensional scale 

dm = df[df$type == "weekly",] 
# As found by previous analysis, we will only use the weekly data  

# Remove NAs
dm <- dm[!is.na(dm$hl),]
range(dm$hl) # Look good now for hl 

dm <- dm[!is.na(dm$ft),]
range(dm$ft) # should be smaller than 10 at least 
dm <- dm[dm$ft <= 10,]
range(dm$ft) # look good for ft 
hist(dm$ft)

dm <- dm[!is.na(dm$sl),]
range(dm$sl) # looks good

dm <- dm[!is.na(dm$nv),]
range(dm$nv) # looks good

dm <- dm[!is.na(dm$sc),]
range(dm$sc) # looks good

# Data frame dm is ready for the MDS scaling method. 
names(dm) # Include hl, ft, sl, sc, nv
# dm[,c(4,5,6,36)] <- scale(dm[,c(4,5,6,8)]) # doesn't work, include autotransform = T instead
mds <- metaMDS(dm[,c(4,5,6,7,30)],autotransform = T)
str(mds)
dm$nmd1 <- mds$points[,1]
dm$nmd2 <- mds$points[,2]
par(mfcol = c(1,1))

dm[dm$Origin == "1",] <- "NM" # Now the wild is assigned as 1 
dm[dm$Origin == "0",] <- "GH" # Now the ex situ is assigned as 0 

fit <- envfit(mds, dm[,19], perm = 1000) # Here use the origin, origin needs to be changes to 0 and 1 to continue the analysis 
fig <- ordiplot(mds, type = "none")
points(fig, "sites", pch = 1, 
       cex= 1.1, col = "black", bg ="white") 
ordihull(mds, dm$Origin, display = "sites", draw = "polygon")

# ggplot figure
accession <- dm$Accession..
label<- dm$Origin
data.scores <- as.data.frame(scores(mds,display = "sites"))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$site <- rownames(data.scores)  # create a column of site names, from the rownames of data.scores
data.scores$Origin <- label  #  add the grp variable created earlier
data.scores$accession <- as.character(accession)

grp.gh <- data.scores[data.scores$Origin == "GH", ][chull(data.scores[data.scores$Origin== "GH", c("NMDS1", "NMDS2")]), ] 
grp.nm <- data.scores[data.scores$Origin == "NM", ][chull(data.scores[data.scores$Origin == "NM", c("NMDS1", "NMDS2")]), ] 
hull.data <- rbind(grp.gh,grp.nm)

p5<- ggplot() + geom_polygon(data=hull.data, aes(x=NMDS1, y=NMDS2, fill = Origin, group = Origin), 
                             alpha = 0.15) +
  geom_point(data = data.scores, aes(x = NMDS1, y = NMDS2, 
                                     colour = accession, shape = Origin ), size = 2.5)+
  theme(legend.position="right")+
  theme(legend.text=element_text(size=4), legend.title=element_text(size=4))
p5 <- p5 + theme(
  panel.background = element_rect(fill = "white", colour = "grey",
                                  size = 2, linetype = "solid"),
  panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                  colour = "grey"), 
  panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                  colour = "grey"))+ xlab("NMDS1") + ylab("NMDS2")
p5 <- p5 + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_fill_manual("Origin", 
                    label = c("ex situ", "wild"),
                    values = c("dodgerblue1","#66CC33"))+
  scale_colour_manual ("Maternal Line", 
                     labels = c("Wild1", "Wild 2", "Wild 3", "Wild 4", 'Ex situ 1', 'Ex situ 2', 'EX situ 3',
                                'Ex situ 4','Ex situ 5', 'Ex situ 6'), 
                     values = c("#FFCC33", "#EC8B00","#BB9D00","#85AD00","#00B81F","#00BF7D","#00BFC4",
                                "#7997FF","#E76BF3","#FF6C90"))+
  theme(legend.key = element_rect(fill = alpha("white", 0.0)))+
  guides(colour=guide_legend(override.aes=list(shape=15)))+
  theme(legend.position = "top",
        legend.text = element_text(size=12))+ 
  theme(axis.text.x =element_text(size = 15),
        axis.text.y = element_text(size = 15))+
  theme(axis.title.y =element_text(size = 15),
        axis.title.x = element_text(size = 15)) 
p5
ggsave("figure6.tiff", dpi = 1200)


# Correlation coefficient 
# since we only use weekly data, 
# so we assigned the whole data as the weekly data 
df <- df[df$type == "weekly",]
df$sh <- df$sl-df$hl

df.nm <- df[df$Origin == "NM",]
df.gh <- df[df$Origin == "GH",]

# floral tube versus style length
cor.test(df.nm$hl, df.nm$sl)
cor.test(df.gh$hl, df.gh$sl)

# floral tube versus floral flare width
cor.test(df.nm$hl, df.nm$ft)
cor.test(df.gh$hl, df.gh$ft)

# floral tube versus style extension
cor.test(df.nm$hl, df.nm$sh)
cor.test(df.gh$hl, df.gh$sh)

# floral tube versus nectar volume
cor.test(df.nm$hl, df.nm$nv)
cor.test(df.gh$hl, df.gh$nv)

# floral tube versus sugar content
cor.test(df.nm$hl, df.nm$sc)
cor.test(df.gh$hl, df.gh$sc)

# style length versus floral flare
cor.test(df.nm$sl, df.nm$ft)
cor.test(df.gh$sl, df.gh$ft)

# style length versus style extension
cor.test(df.nm$sl, df.nm$sh)
cor.test(df.gh$sl, df.gh$sh)

# style length versus nectar volume
cor.test(df.nm$sl, df.nm$nv)
cor.test(df.gh$sl, df.gh$nv)

# style length versus sugar content
cor.test(df.nm$sl, df.nm$sc)
cor.test(df.gh$sl, df.gh$sc)

# floral flare versus style extension
cor.test(df.nm$ft, df.nm$sh)
cor.test(df.gh$ft, df.gh$sh)

# floral flare versus nectar volume
cor.test(df.nm$ft, df.nm$nv)
cor.test(df.gh$ft, df.gh$nv)

# floral flare versus sugar content
cor.test(df.nm$ft, df.nm$sc)
cor.test(df.gh$ft, df.gh$sc)

# style extension versus nectar volume
cor.test(df.nm$sh, df.nm$nv)
cor.test(df.gh$sh, df.gh$nv)

# style extension versus sugar content
cor.test(df.nm$sh, df.nm$sc)
cor.test(df.gh$sh, df.gh$sc)

# nectar volume versus sugar content
cor.test(df.nm$nv, df.nm$sc)
cor.test(df.gh$nv, df.gh$sc)



# HMSC 
names(df)
# df <- df[df$type == "weekly",] This was done when analyzing correlation coefficient 
# The below code is modified from 
# https://github.com/kate-eisen/oenothera/blob/master/AJB_analyses.Rmd: 

response = df %>% 
  dplyr::select (hl, 
                 ft, 
                 sl, 
                 sc,
                 nv)

Y <- scale(response)
X <- as.data.frame(as.factor(df$Origin))

re1 <- HmscRandomLevel(units = df$Accession..)
re2 <- HmscRandomLevel(units = df$Flag_number)
re3 <- HmscRandomLevel(units = df$Days)

ranlevels <- list(accession = re1, plant = re2, date = re3) 

studyDesign = data.frame(accession = as.factor(df$Accession..),
                         plant= as.factor(df$Flag_number),
                         date = as.factor(df$Days))

m = Hmsc(Y=response, 
         XData       = X, 
         studyDesign = studyDesign,
         ranLevels   = list(accession = re1, 
                            plant = re2,
                            date = re3) )

mod_hmsc <- sampleMcmc(m, samples = 1000, thin = 1, transient = 500, nChains = 4)  

# Plot Beta 
post.Beta <- getPostEstimate(mod_hmsc, parName = "Beta")
plotBeta(mod_hmsc, param = "Mean", post = post.Beta) 
mpost = convertToCodaObject(mod_hmsc)
summary(mpost$Beta) 

# variance decomposition 
vp <- computeVariancePartitioning(mod_hmsc) # the warning is fine
view(vp)
dev.off()
par(mar=c(4,4,4,4))
plotVariancePartitioning(mod_hmsc, VP = vp, las = 2, horiz=F)

#evaluate model explanatory power in terms of R2
preds = computePredictedValues(mod_hmsc)
MF = evaluateModelFit(hM=mod_hmsc, predY=preds)
MF$R2
hist(MF$R2, xlim = c(0,1), main=paste0("Mean = ", round(mean(MF$R2),2)))

#extract values of the variance partitioning
vpvals1 <- as.data.frame(vp$vals)
vpvals  <- as.data.frame(t(vpvals1))
colnames(vpvals) <- rownames(vpvals1)
vpvals$compound  <- colnames(vpvals1)

#extract R2 for each compound (i.e. variance explained by the model)
vpvals$R2<-MF$R2

#create a new variable for the variance not explained by the model
vpvals$Unexplained   <- 1-vpvals$R2

#scale the variance explained by my explanatory variables by the R2
#that way, all the variance explained by each variance component 
#plus the unexplained variance will sum to 1
count <- 0
for (i in 1:4) {      #instead of "4" you should add the number of factors you have 
  res <- vpvals[,i]*vpvals$R2
  count <- count + 1
  vpvals[count] <- res
}


#calculate the mean variance explained by each of my variance components
means <- rep(NA, 4)
count <- 0
for (i in 1:4) {
  mean<-mean(vpvals[,i])
  count <- count + 1
  means[count] <- mean
}

#the means start with the component that is on the left in the VPvals df
means
mean(vpvals$Unexplained)

#check that the scaled variances plus the unexplained variance sum to 1
sum(means)+mean(vpvals$Unexplained) # Yes it is summed up to 1 
