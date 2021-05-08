#Author: Smriti Mehta
#Date: April 2021
#Project: Lalot et al. (2018) Replication

# ############################################################
#                     LOADING LIBRARIES
# ############################################################

library(ggplot2) 
library(dplyr) 
library(tidyr) 
library(car)
library(tidyverse)
library(emmeans)
library(TAM)
library(WrightMap)
library(car)
library(psychometric)
library(paramtest)
library(pwr)
library(plyr)
library(data.table)


# ############################################################
#                    DEFINING NEW FUNCTIONS
# ############################################################

NaMean <- function(x) {
  mean(x,na.rm=T)
}

NaSum <- function(x) {
  sum(x,na.rm=T)
}

rev_code_geb_1 <- function(x) {
  x <- 6 - x
}

rev_code_geb_2 <- function(x) {
  x <- 3 - x
}

subtract_1 <- function(x) {
  x <- x - 1
}

na_to_dot <- function(x) {
  ifelse(x %in% NA, ".", x)
}

### OSC function for converting effect sizes
esComp <-function(
  x,df1,df2,N,esType)
{esComp <-ifelse(esType=="t",
                 sqrt((x^2*(1 / df2)) / (((x^2*1) / df2) + 1)),
                 ifelse(esType=="F",
                        sqrt((x*(df1 / df2)) / 
                               (((x*df1) / df2) + 1))*sqrt(1/df1),
                        ifelse(esType=="r",x,
                               ifelse(esType=="Chi2", 
                                      sqrt(x/N),
                                      ifelse(esType == "z",
                                             tanh(x * sqrt(1/(N-3))),
                                             NA)
                               )
                        )
                 )
)
return(esComp)
}


# ############################################################
#                    LOADING DATA 
# ############################################################

##### If reading in clean data, use this code:
#DF <- read_csv("~data/Lalot et al. (2018) clean data.csv")

DF_full <- read_csv("Lalot et al. (2018).csv") %>% 
  filter(!is.na(workerId)) # removing anyone without an MTurk Id

nrow(DF_full) 
#570 participants

# ############################################################
#       RECODING VARIABLES + RETAINING RELEVANT ONES 
# ############################################################

geb_1_vars <- {c("geb_1", "geb_2", "geb_3", "geb_4" ,
               "geb_5", "geb_6", "geb_7", "geb_8",
               "geb_9", "geb_10", "geb_11", "geb_12",
               "geb_13", "geb_14", "geb_15", "geb_16",
               "geb_17", "geb_18", "geb_19", "geb_20",
               "geb_21", "geb_22", "geb_23", "geb_24",
               "geb_25", "geb_26", "geb_27", "geb_28",
               "geb_29", "geb_30", "geb_31", "geb_32",
               "geb_33", "geb_34", "geb_35"
)}
geb_2_vars <- {c("geb_36", "geb_37", "geb_38",
                 "geb_39", "geb_40", "geb_41",
                 "geb_42", "geb_43", "geb_44",
                 "geb_45", "geb_46", "geb_47",
                 "geb_48", "geb_49")}

geb_1_rev_vars <- { c("geb_3_R", "geb_4_R", "geb_5_R",
                          "geb_8_R", "geb_10_R", "geb_11_R", "geb_12_R",
                          "geb_15_R", "geb_19_R", "geb_21_R", "geb_23_R",
                          "geb_24_R", "geb_25_R", "geb_27_R")}
geb_2_rev_vars <- {c("geb_39_R", "geb_42_R", "geb_43_R", "geb_44_R")
}

xmas_eval_vars <- {
  c("xmas_eval_1", "xmas_eval_2", "xmas_eval_3", "xmas_eval_4")
}

dem_vars <- {
  c("age", "gender", "nationality")}

# recoding and retaining relevant variables,
# filtering out those who did not see a condition
{ DF <-  DF_full %>%
  dplyr::mutate_at(geb_1_rev_vars, rev_code_geb_1) %>%
  dplyr::mutate_at(geb_2_rev_vars, rev_code_geb_2) %>%
  dplyr::rename(`geb_3` = geb_3_R,
         `geb_4` = geb_4_R,
         `geb_5`= geb_5_R,
         `geb_8` = geb_8_R,
         `geb_10` = geb_10_R,
         `geb_11` = geb_11_R,
         `geb_12` = geb_12_R,
         `geb_15` = geb_15_R,
         `geb_19` = geb_19_R,
         `geb_21` = geb_21_R,
         `geb_23` = geb_23_R,
         `geb_24` = geb_24_R,
         `geb_25` = geb_25_R,
         `geb_27` = geb_27_R,
         `geb_39` = geb_39_R,
         `geb_42` = geb_42_R,
         `geb_43` = geb_43_R,
         `geb_44` = geb_44_R,
         `condition` = FL_10_DO) %>%
  dplyr::select(ResponseId, condition, attn_check, all_of(dem_vars),
         all_of(geb_1_vars), all_of(geb_2_vars), all_of(xmas_eval_vars)) %>%
    filter(!is.na(condition))
}

#Third item in the DV scale needs to be recoded
DF$xmas_eval_3 <- 8 - DF$xmas_eval_3

#participants that did not answer any of the items
all.na <- apply(DF[,7:59], 1, function(x){all(is.na(x))})

DF <- DF %>% 
  filter(!all.na) #540 rows

table(DF$condition)
#161 control, 188 majority, 191 minority
  

# ############################################################
#                CREATING GEB DFs FOR CONQUEST
# ############################################################


 GEB <- DF %>%
   dplyr::select(all_of(geb_1_vars), all_of(geb_2_vars))

#Creating file for Conquest
# GEB_C <- GEB %>%
#  mutate_all(subtract_1) %>%
#  mutate_all(na_to_dot)

#write_csv(GEB_C, "geb_data.csv")


 GEB_short <- GEB %>%
   dplyr::select(-c(geb_10, geb_11, geb_18))
 
# GEB_short_C <- GEB_short %>%
#  mutate_all(subtract_1) %>%
#  mutate_all(na_to_dot)

#write_csv(GEB_short_C, "geb_short_data.csv")
  

# ############################################################
#        APPLYING RASCH MODEL TO GEB SCALE
# ############################################################

### GEB Complete Scale
geb_mod <- TAM::tam.mml(GEB, irtmodel="PCM2")

geb_mod_fit <- tam.fit(geb_mod)

range(geb_mod_fit$itemfit$Infit)
#infit values range from .997 - 1.00, so fitting *quite* well!

persons.mod <- tam.wle(geb_mod)
#WLE Reliability= 0.665

#person parameters
wle_estimates <- persons.mod$theta

#adding the rasch model estimate to the DF
DF$geb_wle <- wle_estimates


### GEB Scale minus 3 items

geb_short_mod <- TAM::tam.mml(GEB_short, irtmodel="PCM2")

geb_short_mod_fit <- tam.fit(geb_short_mod)

range(geb_short_mod_fit$itemfit$Infit)
#infit values remain the same

persons.mod.short <- tam.wle(geb_short_mod)
#WLE Reliability= 0.626 

#person parameters
wle_estimates_short <- persons.mod.short$theta

#adding the rasch model estimate to the DF
DF$geb_short_wle <- wle_estimates_short


# ############################################################
#                 WRITING CLEAN DATA FILE
# ############################################################

#First, combining the scores for the DV
 DF$xmas_eval <- rowMeans(DF[, xmas_eval_vars], na.rm = T)

#And creating a standardized version of GEB
DF$geb_stn <- scale(DF$geb_wle, center = T, scale = T)
DF$geb_short_stn <- scale(DF$geb_short_wle, center = T, scale = T)

write_csv(DF, "Lalot et al. (2018) clean data.csv")

# ############################################################
#                  MAIN ANALYSIS
# ############################################################


#Changing condition to a factor and adding the right contrasts
DF$condition <- as.factor(DF$condition)

#Setting contrasts 
contr1 <- c(-1, 2, -1)
contr2 <- c(1, 0, -1)
mat <- cbind(contr1, contr2)

contrasts(DF$condition) <- mat

#Creating the model
mod <- lm(xmas_eval ~  geb_stn + condition + geb_stn*condition, data = DF)

mod.sum <- summary(mod) 

#GEB scale has a significant main effect, 
#along with the interaction of GEB and the contrrast that is NOT of interest 

Anova(mod, type = "III")

# SE and Confidence intervals for the model coefficents
confint(mod)


##### SE FOR MEANS

# Original
o.n.maj <- 69
o.sd.maj <- 0.96
o.se.maj <- round(o.sd.maj/(sqrt(o.n.maj)), 2)

o.n.min <- 71
o.sd.min <- 1.25
o.se.min <- round(o.sd.min/(sqrt(o.n.min)), 2)

o.n.con <- 70
o.sd.con <- 1.29
o.se.con <- round(o.sd.con/(sqrt(o.n.con)), 2)

# Replication
r.n.maj <- nrow(DF[DF$condition %in% "MajorityCondition",])
r.sd.maj <- sd(DF[DF$condition %in% "MajorityCondition",]$xmas_eval, na.rm = T)
r.se.maj <- round(r.sd.maj/(sqrt(r.n.maj)), 2)

r.n.min <- nrow(DF[DF$condition %in% "MinorityCondition",])
r.sd.min <- sd(DF[DF$condition %in% "MinorityCondition",]$xmas_eval, na.rm = T)
r.se.min <- round(r.sd.min/(sqrt(r.n.min)), 2)

r.n.con <- nrow(DF[DF$condition %in% "ControlCondition",])
r.sd.con <- sd(DF[DF$condition %in% "ControlCondition",]$xmas_eval, na.rm = T)
r.se.con <- round(r.sd.con/(sqrt(r.n.con)), 2)

# ############################################################
#                  SIMPLE SLOPES ANALYSIS
# ############################################################


# Simple slopes at GEB = 0 
emtrends(mod, "condition", var="geb_stn")

# Comparing slopes 
emtrends(mod, pairwise ~ condition, var="geb_stn")
#none of the comparisons are significant


(mylist <- list(geb_stn=c(-1, 1), 
                condition=c("ControlCondition",
                            "MinorityCondition","MajorityCondition")))
emmod <- emmeans(mod, ~ geb_stn*condition, at=mylist)
contrast(emmod, "pairwise",by="condition") 

# Plot using emmip
emmip(mod, condition ~ geb_stn, at=mylist,CIs=TRUE)


# save simple slopes as dataset 
moddat <- emmip(mod,geb_stn~condition,at=mylist, CIs=TRUE, plotit=FALSE)

# plot
(p <- ggplot(data=moddat, aes(x=geb_stn,y=yvar, color=condition)) + geom_line())


#Simple simple slopes 

mod.ss <- lm(xmas_eval ~  0 + condition + geb_stn:condition, data = DF)
summary(mod.ss)
confint(mod.ss)



DF$geb.plus <- DF$geb_stn + 1
DF$geb.minus <- DF$geb_stn - 1

m.plus <- lm(xmas_eval ~ geb.plus*condition, data = DF)
summary(m.plus)
confint(m.plus)

m.minus <- lm(xmas_eval ~ geb.minus*condition, data = DF)
summary(m.minus)
confint(m.minus)


summary(lm(geb_wle ~ condition, data = DF))
# ############################################################
#                       EFFECT SIZES
# ############################################################

# **Lalot et al. Study 1**   
#original stats:  
# F(5,204) = 9.09
#t-value for the interaction of interest
# t(204) = −2.96
# η2 = .04 
# N = 210

#replication stats:  
# t(521) = 0.843
# η2 = .01
# N = 540


### original study
Lalot.orig.es <- esComp(x = -2.96, df2 = 204, N = 210, esType = "t")
Lalot.orig.es # 0.2029295

# calculate 95% CI
CIr(r= Lalot.orig.es, n = 210, level = .95) # 0.06944706 - 0.32927320


### replication study
Lalot.rep.es <- esComp(x = 0.843, df2 = 521, N = 540, esType = "t")
Lalot.rep.es #0.03690734

# calculate 95% CI
CIr(r= Lalot.rep.es, n = 540, level = .95) # -0.04761854  0.12090840


Lalot.rep.upper <- round( # store upper bound of CI
  CIr(r = Lalot.rep.es, 
      n = 540, 
      level = .95)[2], # calculate 95% CI, extract upper bound ([2])
  7) # round to 7 digits
#0.1209084



### replication study - Contrast NOT of interest but significant
Lalot.rep.non.es <- esComp(x = 2.133, df2 = 521, N = 540, esType = "t")
Lalot.rep.non.es #0.093

# calculate 95% CI
CIr(r= Lalot.rep.non.es, n = 540, level = .95) # 0.008734039 - 0.176038605

Lalot.rep.non.upper <- round(CIr(r = Lalot.rep.non.es, 
      n = 540, 
      level = .95)[2], # calculate 95% CI, extract upper bound ([2])
  7) # round to 7 digits
#0.1760386

# ############################################################
#                      SMALL TELESCOPES 
# ############################################################


# POINT ESTIMATE 
# original study's power to detect replication effect 
pwr.r.test(n = 210, r = Lalot.rep.es, sig.level = .05) # 0.08296258

# N needed for 80% power to detect effect
pwr.r.test(r = Lalot.rep.es, sig.level = .05, power = .80) # 5758.97


### UPPER BOUND 
# original study's power to detect replication upper bound effect size
pwr.r.test(n = 210, r = Lalot.rep.upper, sig.level = .05) # 0.4175867

# N needed for 80% power to detect effect
pwr.r.test(r = Lalot.rep.upper, sig.level = .05, power = .80) #533.7658


### CONTRAST NOT OF INTEREST

# POINT ESTIMATE 
pwr.r.test(n = 210, r = Lalot.rep.non.es, sig.level = .05) # 0.2697915
# N needed for 80% power to detect effect
pwr.r.test(r = Lalot.rep.non.es, sig.level = .05, power = .80) # 903.5135


### UPPER BOUND 
pwr.r.test(n = 210, r = Lalot.rep.non.upper, sig.level = .05) # 0.7273663
# N needed for 80% power to detect effect
pwr.r.test(r = Lalot.rep.non.upper, sig.level = .05, power = .80) #250.1395


# ############################################################
#                ANALYSIS WITH EXCLUSIONS 
# ############################################################

DF_excl <- DF %>% 
  filter(attn_check %in% 4) #14 people were excluded

#Creating the model
mod.excl <- lm(xmas_eval ~  geb_stn + condition + geb_stn*condition, data = DF_excl)

summary(mod.excl) 

Anova(mod.excl, type = "III")

confint(mod.excl)

# same results, main effect of prior green behavior and a significant interaction (not the one hypothesized).

#effect size
### replication study - contrast of interest
Lalot.rep.excl.es <- esComp(x = 0.862, df2 = 510, N = 526, esType = "t")
Lalot.rep.excl.es #0.03814223

# calculate 95% CI
CIr(r= Lalot.rep.excl.es, n = 526, level = .95) # -0.04750671  0.12323439


### replication study - contrast not of interest
Lalot.rep.excl.es.2 <- esComp(x = 2.156, df2 = 510, N = 526, esType = "t")
Lalot.rep.excl.es.2 #0.09503716

# calculate 95% CI
CIr(r= Lalot.rep.excl.es.2, n = 526, level = .95) #  0.009621314 0.179076175

# ############################################################
#          ANALYSIS AFTER EXCLUDING THREE GEB ITEMS
# ############################################################

#Creating a new model with the short scale
mod.short <- lm(xmas_eval ~ geb_short_stn*condition, data = DF)

summary(mod.short)

Anova(mod.short, type = "III")

confint(mod.short)
#results remain unchanged. 



# ############################################################
#        COMPARING CONQUEST ESTIMATES WITH TAM ESTIMATES
# ############################################################

conquest.est <- read_csv("geb_conquest_estimates.txt", col_names = "con_wle_est")

est.compare <- cbind(conquest.est, wle_estimates)

est.compare$diff <- est.compare$con_wle_est - est.compare$wle_estimates

summary(est.compare)


### Running the model again with EAP estimates
geb_eap <- read_csv("geb_eap.txt", col_names = F)

DF$geb_eap <- geb_eap %>% 
  pull(X1)

mod.eap <- lm(xmas_eval ~  scale(geb_eap) + condition + scale(geb_eap)*condition, data = DF)

summary(mod.eap)

Anova(mod.eap, type = "III")
# Same results

# ############################################################
#                    TABLES AND PLOTS
# ############################################################

NaMean(DF$age)
sd(DF$age)

#gender 1 = Male, 2 = Female, 3 = Non-bi
#4 = prefer to self describe, 5=prefer not to say
DF %>% 
  group_by(gender) %>% 
  dplyr::summarize(N = n(), 
            `%` = round( ((N/nrow(DF))*100 ), 2)) %>% 
  View()

mod.plot <- lm(xmas_eval ~ geb_wle*condition, DF)

Plus_1 <- NaMean(DF$geb_wle) + sd(DF$geb_wle)
Minus_1 <- NaMean(DF$geb_wle) - sd(DF$geb_wle)

newdf <- data.table(condition=rep(c('ControlCondition',
                                    'MajorityCondition', 
                                    'MinorityCondition'), each=2), 
                    geb_wle = c(Minus_1, Plus_1))

newdf$pred <- predict(mod.plot, newdata = newdf,se.fit=TRUE)$fit

newdf <- rbind(newdf[1:2,], newdf[5:6,], newdf[3:4,])

newdf$condition <- ifelse(newdf$condition %in% "ControlCondition", "Control condition", 
       ifelse(newdf$condition %in% "MinorityCondition", "Minority support", 
              ifelse(newdf$condition %in% "MajorityCondition", "Majority support", newdf$condition)))

newdf$condition <- factor(newdf$condition, levels = c("Control condition", "Minority support", "Majority support"))



p <-  ggplot(newdf, aes(x=geb_wle, y=pred, col=condition)) +
  geom_point(size = 2) +
  geom_path(aes(linetype=condition), size = 1) + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"), 
        legend.title = element_blank(), 
        panel.background = element_rect(fill = "white"), 
        panel.grid.minor.y = element_line(color = "gray"), 
        legend.background = element_rect(fill = "white"),
        legend.text = element_text(size = 12), 
        axis.title.x = element_text(size = 15), 
        text = element_text(family = "Times")) + 
  scale_linetype_manual(values=c("dashed", "solid", "solid")) + 
  scale_color_manual(values = c("darkgray", "darkgray", "black")) +
  ylim(c(4, 6.5)) + 
  xlim(c(-.4, .4)) + 
  ggtitle("Attitude towards the collective action (Green Christmas)") + 
  xlab("Green behaviour -1SD      Green behaviour +1SD") + 
  ylab(NULL) 


ggsave( "Interaction Plot.pdf", plot = p, units="in", width=6.3, height=4, dpi=300)
