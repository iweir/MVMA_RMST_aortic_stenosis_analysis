# REPRODUCILIBILITY CODE FOR Multivariate meta-analysis model for the difference in restricted mean survival times
# Isabelle Weir, Lu Tian, Ludovic Trinquart


library(metaRMST)

# read in built-in dataset 
data(AorticStenosisTrials)

# test for proportional hazards

## NOTION (trialID==1)
cox1 <- coxph(formula = Surv(Time, Event) ~ Arm, data = AorticStenosisTrials[which(AorticStenosisTrials$trialID==1),]) 
PHtest_p_NOTION<- cox.zph(cox1)$table[3]

## PARTNER (trialID==2)
cox2 <- coxph(formula = Surv(Time, Event) ~ Arm, data = AorticStenosisTrials[which(AorticStenosisTrials$trialID==2),]) 
PHtest_p_PARTNER <- cox.zph(cox2)$table[3]

## SURTAVI (trialID==3)
cox3 <- coxph(formula = Surv(Time, Event) ~ Arm, data = AorticStenosisTrials[which(AorticStenosisTrials$trialID==3),]) 
PHtest_p_SURTAVI<- cox.zph(cox3)$table[3]

## PARTNER2 (trialID==4)
cox4 <- coxph(formula = Surv(Time, Event) ~ Arm, data = AorticStenosisTrials[which(AorticStenosisTrials$trialID==4),]) 
PHtest_p_PARTNER2<- cox.zph(cox4)$table[3]

## USCoreValve (trialID==5)
cox5 <- coxph(formula = Surv(Time, Event) ~ Arm, data = AorticStenosisTrials[which(AorticStenosisTrials$trialID==5),]) 
PHtest_p_USCoreValve<- cox.zph(cox5)$table[3]

## overall test for NPH (fisher method for combining p values

pvalues <- c(PHtest_p_NOTION, PHtest_p_PARTNER, PHtest_p_SURTAVI , PHtest_p_PARTNER2 , PHtest_p_USCoreValve )

C <- (-2)*sum(log(pvalues))
df <- 2*5 # number of trials in meta-anlaysis =5

crit <- qchisq(0.95,df)
Combpval <- 1-pchisq(C, df)
Combpval

# Meta-analysis using each of the 4 methods

# demonstration of meta-analysis to obtain combined effect by multivariate meta-analysis model (method="mvma")
mvma_res <- metaRMSTD(AorticStenosisTrials, time_horizons=c(12,24,36), MA_method="mvma")
mvma_res$REresult

# meta-analysis to obtain combined effect by multivariate meta-analysis model (method="mvma_boot")
boot_res <- metaRMSTD(AorticStenosisTrials, time_horizons=c(12,24,36), MA_method="mvma_boot", nboot=500)
boot_res$REresult

# meta-analysis to obtain combined effect by univariate meta-analysis model (method="uni")
uni_res <- metaRMSTD(AorticStenosisTrials, time_horizons=c(12,24,36), MA_method="uni")
uni_res$result

# meta-analysis to obtain combined effect by univariate meta-analysis model of flexible parametric models (method="uni_flex")
uni_flex_res <- metaRMSTD(AorticStenosisTrials, time_horizons=c(12,24,36), MA_method="uni_flex")
uni_flex_res$result
