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
library(TAM)
library(WrightMap)
library(car)
library(psychometric)
library(paramtest)
library(pwr)
library(plyr)


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
DF <- { DF_full %>%
  dplyr::mutate_at(geb_1_rev_vars, rev_code_geb_1) %>%
  dplyr::mutate_at(geb_2_rev_vars, rev_code_geb_2) %>%
  rename(`geb_3` = geb_3_R,
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
  select(ResponseId, condition, attn_check, all_of(dem_vars),
         all_of(geb_1_vars), all_of(geb_2_vars), all_of(xmas_eval_vars)) %>%
    filter(!is.na(condition))
}

#Third items in the DV scale needs to be recoded
DF$xmas_eval_3 <- 8 - DF$xmas_eval_3

#participants that did not answer any of the items
all.na <- apply(DF[,7:59], 1, function(x){all(is.na(x))})

DF <- DF %>% 
  filter(!all.na) #540 rows

table(DF$condition)
#161 control, 188 majority, 191 minority
  

# ############################################################
#                CREATING GEB DF FOR CONQUEST
# ############################################################


GEB <- DF %>%
  select(all_of(geb_1_vars), all_of(geb_2_vars))

#Creating file for Conquest
GEB_C <- GEB %>%
 mutate_all(subtract_1) %>%
 mutate_all(na_to_dot)

#write_csv(GEB_C, "geb_data.csv")


GEB_short <- GEB %>%
  select(-c(geb_10, geb_11, geb_18))

GEB_short_C <- GEB_short %>%
 mutate_all(subtract_1) %>%
 mutate_all(na_to_dot)

#write_csv(GEB_short_C, "geb_short_data.csv")
  

# ############################################################
#        APPLYING RASCH MODEL TO GEB SCALE
# ############################################################

### GEB Complete Scale
geb_mod <- TAM::tam.mml(GEB, irtmodel="PCM2", )

geb_mod_fit <- tam.fit(geb_mod)

range(geb_mod_fit$itemfit$Infit)
#infit values range from .997 - 1.00, so fitting *quite* well!

persons.mod <- tam.wle(geb_mod)
#WLE Reliability= 0.665

#person parameters
wle_estimates <- persons.mod$theta

#adding the rasch model estimate to the DF
DF$geb_wle <-wle_estimates



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
mod <- lm(xmas_eval ~  geb_wle + condition + geb_wle*condition, data = DF)

summary(mod) 

Anova(mod, type = "III")


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
Lalot.orig.es #.20

# calculate 95% CI
CIr(r= Lalot.orig.es, n = 210, level = .95) # 0.06944706 - 0.32927320


### replication study
Lalot.rep.es <- esComp(x = 0.843, df2 = 521, N = 540, esType = "t")
Lalot.rep.es #0.04


# calculate 95% CI
CIr(r= Lalot.rep.es, n = 540, level = .95) # -0.04761854  0.12090840


Lalot.rep.upper <- round( # store upper bound of CI
  CIr(r = Lalot.rep.es, 
      n = 540, 
      level = .95)[2], # calculate 95% CI, extract upper bound ([2])
  7) # round to 7 digits
#0.1209084



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


# ############################################################
#                ANALYSIS WITH EXCLUSIONS 
# ############################################################

DF_excl <- DF %>% 
  filter(attn_check %in% 4) #14 people were excluded

#Creating the model
mod.excl <- lm(xmas_eval ~  geb_wle + condition + geb_wle*condition, data = DF_excl)

summary(mod.excl) 

Anova(mod.excl, type = "III")

#same results, main effect of prior green behavior. 


# ############################################################
#          ANALYSIS AFTER EXCLUDING THREE GEB ITEMS
# ############################################################

#Creating a new model with the short scale
mod.short <- lm(xmas_eval ~  geb_short_wle + condition + geb_short_wle*condition, data = DF)

summary(mod.short)

Anova(mod.short, type = "III")

#results remain unchanged. 




# ############################################################
#        COMPARING CONQUEST ESTIMATES WITH TAM ESTIMATES
# ############################################################

conquest.est <- read_csv("geb_conquest_estimates.txt", col_names = "con_wle_est")

est.compare <- cbind(conquest.est, wle_estimates)

est.compare$diff <- est.compare$con_wle_est - est.compare$wle_estimates

View(est.compare)

