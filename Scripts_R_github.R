
#here are the scripts used in the study: 'Maternal trade-off between current and future reproduction over seasonal birth timing in a long-lived primate'

#per main section, and with annotations about the models ran

library(CircStats)
library(circular)
library(blme)
library(ggplot2)
library(mgcv)
library(MuMIn)
library(lubridate)
library(car)
library(DHARMa)

############################################################################################################
######### 1. Tsaobis baboons breed year-round despite living in a seasonal environment
############################################################################################################

#seasonality of conceptions, births and cycle resumptions (taking into account their uncertainties)

#here we only present the code for the case of conceptions, but it the exact same method used for births and cycle resumptions.

datatotal=Conceptions #open the conception excel file. 

i=0
vecteur=c()
for (i in 0:999) {  #we draw 1000 randomised birth dates
  i=i+1  
  j=0
  for (j in 0:nrow(datatotal)-1) {  
    j=j+1
    if (length(datatotal$Conc_Uncertainty_Days[j]) > 0 &&  datatotal$Conc_Uncertainty_Days[j]==0) { 
      vecteur[j]=datatotal$Conception_Radian[j] 
    }
    
    if (length(datatotal$Conc_Uncertainty_Days[j]) > 0 && datatotal$Conc_Uncertainty_Days[j]>0) {  ## ici pour date non precise
      vecteur[j]=runif(1,datatotal$Conception_Radian[j]-datatotal$Uncertainty_Conc_DayRad[j],datatotal$Conception_Radian[j]+datatotal$Uncertainty_Conc_DayRad[j])  #ici, on tire au sort selon le nombre de jour d'incertitude
      if(vecteur[j]>2*pi) {vecteur[j]=vecteur[j]-2*pi } #puis on retombe sur une valeur entre 0 et 2pi, en radian
      if (vecteur[j]<0) {vecteur[j]=vecteur[j]+2*pi}
    }
  }
  
  datatotal=cbind(datatotal,vecteur) #we add the randomised vector in our table. 
  
}


tab_test=data.frame()

for (i in (nrow(datatotal)-999):nrow(datatotal)) {   
  i=i+1
  t=data.frame()  
  
  vecteur=datatotal[,i]
  test=r.test(vecteur)
  r=test$r.bar
  p=test$p.value
  mu=circ.mean(vecteur)
  
  t=data.frame(r,p,mu)
  
  tab_test=rbind(tab_test,t)
}

mean(tab_test$r) #mean value of the resultant length of the rayleigh test

mean(tab_test$p) #mean pvalues of the Rayleigh test
wilcox.test(tab_test$p, conf.int = TRUE, conf.level = 0.95) #95% confidence intervals for the pvalues. 

mean(tab_test$mu) #mean conception date, in radians, which can be converted in days. 




############################################################################################################
######### 2. Distinct birth timings optimize current versus future reproduction 
############################################################################################################

########## Model 1: influence of birth dates on offspring mortality

TAB='Mortality' #open the excel file 'Mortality'

TAB=TAB[which(!is.na(TAB$Death_18months)==TRUE), ] #we only kept the 195 infants for which survival outcome is known
TAB$Death_18months=as.factor(TAB$Death_18months)
TAB$Death_18months=factor(TAB$Death_18months)
TAB$Mother=factor(TAB$Mother)
TAB$Relative_Rank=scale(TAB$Relative_Rank, center=TRUE, scale=TRUE)
TAB$Year_Birth=as.factor(as.character(TAB$Year_Birth))
summary(TAB$Year_Birth)
TAB$DOB=dmy(TAB$DOB)
TAB$First_Not_Seen=dmy(TAB$First_Not_Seen)
TAB$Last_Seen=dmy(TAB$Last_Seen)


#first step: checking which phase of the sine fixed effect (birth date) minimized model's AIC : 

MOD0=glmer(Death_18months ~ sin(Birth_radian+0*pi/12) + Troop + 
             Sex + Parity + Relative_Rank +
             (1|Mother) +  (1|Year_Birth),
           glmerControl(optimizer = 'bobyqa'), data=TAB, family='binomial')
MOD1=glmer(Death_18months ~ sin(Birth_radian+1*pi/12) + Troop + 
             Sex + Parity + Relative_Rank +
             (1|Mother) + (1|Year_Birth),
           glmerControl(optimizer = 'bobyqa'), data=TAB, family='binomial')
MOD2=glmer(Death_18months ~ sin(Birth_radian+2*pi/12) + Troop + 
             Sex + Parity + Relative_Rank +
             (1|Mother) + (1|Year_Birth),
           glmerControl(optimizer = 'bobyqa'), data=TAB, family='binomial')
MOD3=glmer(Death_18months ~ sin(Birth_radian+3*pi/12) + Troop + 
             Sex + Parity + Relative_Rank +
             (1|Mother) + (1|Year_Birth),
           glmerControl(optimizer = 'bobyqa'), data=TAB, family='binomial')
MOD4=glmer(Death_18months ~ sin(Birth_radian+4*pi/12) + Troop + 
             Sex + Parity + Relative_Rank +
             (1|Mother) + (1|Year_Birth),
           glmerControl(optimizer = 'bobyqa'), data=TAB, family='binomial')
MOD5=glmer(Death_18months ~ sin(Birth_radian+5*pi/12) + Troop + 
             Sex + Parity + Relative_Rank +
             (1|Mother) + (1|Year_Birth),
           glmerControl(optimizer = 'bobyqa'), data=TAB, family='binomial')
MOD6=glmer(Death_18months ~ sin(Birth_radian+6*pi/12) + Troop + 
             Sex + Parity + Relative_Rank +
             (1|Mother) + (1|Year_Birth),
           glmerControl(optimizer = 'bobyqa'), data=TAB, family='binomial')
MOD7=glmer(Death_18months ~ sin(Birth_radian+7*pi/12) + Troop + 
             Sex + Parity + Relative_Rank +
             (1|Mother) +  (1|Year_Birth),
           glmerControl(optimizer = 'bobyqa'), data=TAB, family='binomial')
MOD8=glmer(Death_18months ~ sin(Birth_radian+8*pi/12) + Troop + 
             Sex + Parity + Relative_Rank +
             (1|Mother) + (1|Year_Birth),
           glmerControl(optimizer = 'bobyqa'), data=TAB, family='binomial')
MOD9=glmer(Death_18months ~ sin(Birth_radian+9*pi/12) + Troop + 
             Sex + Parity + Relative_Rank +
             (1|Mother) +  (1|Year_Birth),
           glmerControl(optimizer = 'bobyqa'), data=TAB, family='binomial')
MOD10=glmer(Death_18months ~ sin(Birth_radian+10*pi/12) + Troop + 
             Sex + Parity + Relative_Rank +
             (1|Mother) + (1|Year_Birth),
           glmerControl(optimizer = 'bobyqa'), data=TAB, family='binomial')
MOD11=glmer(Death_18months ~ sin(Birth_radian+11*pi/12) + Troop + 
             Sex + Parity + Relative_Rank +
             (1|Mother) +  (1|Year_Birth),
           glmerControl(optimizer = 'bobyqa'), data=TAB, family='binomial')

AIC(MOD0,MOD1,MOD2,MOD3,MOD4,MOD5,MOD6,MOD7,MOD8,MOD9,MOD10,MOD11)

#the model 7 is selected as it minimizes AIC, and we then run this model, with our associated tests: 
summary(MOD7)
Anova(MOD7)
confint(MOD7,method='Wald')

#check model assumptions: 
simulationOutput <- simulateResiduals(fittedModel = MOD7, plot = T)



#here we controlled for uncertainties in the birth dates and for uncertainties in survival outcomes resulting from such uncertainty. 
Pvalue=c()
for (i in 1:1000) {
  
  print(i)
  
  Death=c()
  Rad_Value=c()
  
  for (j in 1:nrow(TAB)) {
    
    DOBj=ymd(as.Date('1970-01-01') +
               days(floor(runif(1,TAB$DOB[j]-TAB$Uncertainty_Birth_Days[j]/2,TAB$DOB[j]+TAB$Uncertainty_Birth_Days[j]/2))))
    Originj=ymd(paste(as.character(year(DOBj)),'01','01',sep='-'))
    Rad_Value=c(Rad_Value,(as.numeric(julian(DOBj, origin = Originj))+1)*2*pi/365.25)
    
    if (TAB$Death_18months[j]=='0') { Death=c(Death,0) }
    if (!TAB$Death_18months[j]=='0') { 
      
      Date_Deathj=ymd(as.Date('1970-01-01') +
                        days(floor(runif(1,TAB$Last_Seen[j],TAB$First_Not_Seen[j]))))
      
      Date_18monthsj= DOBj + 550
      
      if (Date_Deathj < (Date_18monthsj + 1) ) {Death=c(Death,1)}
      if (Date_Deathj > Date_18monthsj ) {Death=c(Death,0)}
      
    }
    
  }
  
  MOD=glmer(Death ~ sin(Rad_Value+7*pi/12) + Troop + 
              Sex + Parity + Relative_Rank +
              (1|Mother) + (1|Year_Birth),
            glmerControl(optimizer = 'bobyqa'), data=TAB, family='binomial')
  Pvalue=c(Pvalue,Anova(MOD)[1,3])
  
}
#we extract the 1000 pvalues of the 1000 models ran controlling for birth date uncertainties
#to confirm the significant effect of birth dates (sin(Rad_Value + 7*pi/12). 
summary(Pvalue)
wilcox.test(Pvalue, conf.int = TRUE, conf.level = 0.95)


########## Model 2: influence of birth dates on maternal IBIs

TAB='IBI' #open the excel file 'IBI

#we first discard some IBIs with first infant who died or second infant who was born dead. 
TAB$DeathInBw=as.factor(TAB$DeathInBw)
TAB=TAB[which(TAB$DeathInBw=='0'),] 
TAB=TAB[which(TAB$Death1=="n"),]
TAB=TAB[which(TAB$Borndead=='n'),] 

TAB$Mother=factor(TAB$Mother)
TAB$Relative_Rank=scale(TAB$Relative_Rank)
TAB$Year_Birth=as.factor(as.character(TAB$Year_Birth))

TAB$DOB1=dmy(TAB$DOB1)
TAB$DOB2=dmy(TAB$DOB2)


#first step: checking which phase of the sine fixed effect (birth date) minimized model's AIC (like for Model 1) : 

MOD0=lmer(IBIDays ~ sin(Birth_Rad) + 
            Sex1 + Parity + Relative_Rank + Troop +
            (1|Mother) + (1|Year_Birth), data=TAB)
MOD1=lmer(IBIDays ~ sin(Birth_Rad + 1*pi/12) + 
            Sex1 + Parity + Relative_Rank + Troop +
            (1|Mother) + (1|Year_Birth), data=TAB)
MOD2=lmer(IBIDays ~ sin(Birth_Rad + 2*pi/12) + 
            Sex1 + Parity + Relative_Rank + Troop +
            (1|Mother) + (1|Year_Birth), data=TAB)
MOD3=lmer(IBIDays ~ sin(Birth_Rad + 3*pi/12) + 
            Sex1 + Parity + Relative_Rank + Troop +
            (1|Mother) + (1|Year_Birth), data=TAB)
MOD4=lmer(IBIDays ~ sin(Birth_Rad + 4*pi/12) + 
            Sex1 + Parity + Relative_Rank + Troop +
            (1|Mother) + (1|Year_Birth), data=TAB)
MOD5=lmer(IBIDays ~ sin(Birth_Rad + 5*pi/12) + 
            Sex1 + Parity + Relative_Rank + Troop +
            (1|Mother) + (1|Year_Birth), data=TAB)
MOD6=lmer(IBIDays ~ sin(Birth_Rad + 6*pi/12) + 
            Sex1 + Parity + Relative_Rank + Troop +
            (1|Mother) + (1|Year_Birth), data=TAB)
MOD7=lmer(IBIDays ~ sin(Birth_Rad + 7*pi/12) + 
            Sex1 + Parity + Relative_Rank + Troop +
            (1|Mother) + (1|Year_Birth), data=TAB)
MOD8=lmer(IBIDays ~ sin(Birth_Rad + 8*pi/12) + 
            Sex1 + Parity + Relative_Rank + Troop +
            (1|Mother) + (1|Year_Birth), data=TAB)
MOD9=lmer(IBIDays ~ sin(Birth_Rad + 9*pi/12) + 
            Sex1 + Parity + Relative_Rank + Troop +
            (1|Mother) + (1|Year_Birth), data=TAB)
MOD10=lmer(IBIDays ~ sin(Birth_Rad + 10*pi/12) + 
             Sex1 + Parity + Relative_Rank + Troop +
             (1|Mother) + (1|Year_Birth), data=TAB)
MOD11=lmer(IBIDays ~ sin(Birth_Rad + 11*pi/12) + 
             Sex1 + Parity + Relative_Rank + Troop +
             (1|Mother) + (1|Year_Birth), data=TAB)

AIC(MOD0,MOD1,MOD2,MOD3,MOD4,MOD5,MOD6,MOD7,MOD8,MOD9,MOD10,MOD11)
#the second model, MOD2, is selected here. 
summary(MOD2)
Anova(MOD2)
confint(MOD2,method='Wald')

#check model assumptions: 
simulationOutput <- simulateResiduals(fittedModel = MOD2, plot = T)

#here we confirm the previous results, controlling for uncertainties in birth dates, resulting in 
#uncertainties in IBI values, & in birth date fixed effect (sine effect)
Pvalue=c()
for (i in 1:1000) {
  
  print(i)
  IBI_Values=c()
  Rad_Value=c()
  
  for (j in 1:nrow(TAB)) {
    
    DOB1j=ymd(as.Date('1970-01-01') +
                days(floor(runif(1,TAB$DOB1[j]-TAB$Uncertainty_Birth_Days[j]/2,TAB$DOB1[j]+TAB$Uncertainty_Birth_Days[j]/2))))
    DOB2j=ymd(as.Date('1970-01-01') +
                days(floor(runif(1,TAB$DOB2[j]-TAB$Uncertainty_Birth_DaysBis[j]/2,TAB$DOB2[j]+TAB$Uncertainty_Birth_DaysBis[j]/2))))
    
    IBI_Values=c(IBI_Values,as.numeric(DOB2j-DOB1j))
    Originj=ymd(paste(as.character(year(DOB1j)),'01','01',sep='-'))
    Rad_Value=c(Rad_Value,(as.numeric(julian(DOB1j, origin = Originj))+1)*2*pi/365.25)
    
  }
  
  MOD=lmer(IBI_Values ~ sin(Rad_Value + 2*pi/12) + 
             Sex1 + Parity + Relative_Rank + Troop +
             (1|Mother) + (1|Year_Birth), data=TAB)
  Pvalue=c(Pvalue,Anova(MOD)[1,3])
  
}
summary(Pvalue)
wilcox.test(Pvalue, conf.int = TRUE, conf.level = 0.95) #confidence interval of the 1000 pvalues of the sine term fixed effect


########## Model 3: individual determinants of giving birth close to the offspring survival optimal timing

TAB='Births' #open the excel file 'Births'

TAB$DOB=dmy(TAB$DOB)
TAB$DOB_Moyen_Offspring=dmy(TAB$DOB_Moyen_Offspring)
summary(TAB$DOB_Moyen_Offspring)
TAB$Relative_Rank=as.numeric(scale(TAB$Relative_Rank))

#the model:
MOD=blmer(abs(Lag_Offspring) ~ Sex + Parity + Relative_Rank +   Troop +
            (1|Mother) +  (1|Year_Birth), data=TAB)
summary(MOD)
Anova(MOD)
confint(MOD, method='Wald')
simulationOutput <- simulateResiduals(fittedModel = MOD, plot = T)

#we control here for uncertainties in birth dates (affecting the response variable) :

#wa had randomized birth dates in our table
for (i in 1:1000) {
  print(i)
  vec=c()
  
  for (j in 1:nrow(TAB)) { 
    
    DATEi=sample(seq(TAB$DOB[j]-0.5*TAB$Uncertainty_Birth_Days[j],TAB$DOB[j]+0.5*TAB$Uncertainty_Birth_Days[j],by='day'),1)
    
    LAGi=as.numeric(DATEi-TAB$DOB_Moyen_Offspring[j])
    
    if (LAGi>182) { LAGi=-183+(LAGi-182) }
    if (LAGi < -182) {LAGi=183+(LAGi+182)}
    vec=c(vec,LAGi)
    
  }
  TAB=cbind(TAB,vec)
}


Estimate_Rank=c()
Estimate_Parite=c()
Estimate_Sex=c()
Pvalue_Rank=c()
Pvalue_Parite=c()
Pvalue_Sex=c()

for (j in (ncol(TAB)-999):ncol(TAB)) {
  print(j)
  
  MOD=blmer(abs(TAB[,j]) ~ Relative_Rank + Parity + Sex + Troop +
              (1|Mother) +  (1|Year_Birth), data=TAB)
  S=summary(MOD)
  Estimate_Rank=c(Estimate_Rank,S$coefficients[2,1])
  Pvalue_Rank=c(Pvalue_Rank,Anova(MOD)[1,3])
  
  Estimate_Parite=c(Estimate_Parite,S$coefficients[3,1])
  Pvalue_Parite=c(Pvalue_Parite,Anova(MOD)[2,3])
  
  Estimate_Sex=c(Estimate_Sex,S$coefficients[4,1])
  Pvalue_Sex=c(Pvalue_Sex,Anova(MOD)[3,3])
  
  
}
#here we can summarize the Pvalues and estimates of our fixed effects. 


########## Model 4: individual determinants of giving birth close to the maternal IBI optimal timing

TAB='Births' #open the excel file 'Births'

TAB$DOB=dmy(TAB$DOB)
TAB$DOB_Moyen_Offspring=dmy(TAB$DOB_Moyen_Mother)
summary(TAB$DOB_Moyen_Mother)
TAB$Relative_Rank=as.numeric(scale(TAB$Relative_Rank))

#the model:
MOD=blmer(abs(Lag_Mother) ~ Sex + Parity + Relative_Rank +   Troop +
            (1|Mother) +  (1|Year_Birth), data=TAB)
summary(MOD)
Anova(MOD)
confint(MOD, method='Wald')
simulationOutput <- simulateResiduals(fittedModel = MOD, plot = T)

#we control here for uncertainties in birth dates (affecting the response variable) :

#wa had randomized birth dates in our table
for (i in 1:1000) {
  print(i)
  vec=c()
  
  for (j in 1:nrow(TAB)) { 
    
    DATEi=sample(seq(TAB$DOB[j]-0.5*TAB$Uncertainty_Birth_Days[j],TAB$DOB[j]+0.5*TAB$Uncertainty_Birth_Days[j],by='day'),1)
    
    LAGi=as.numeric(DATEi-TAB$DOB_Moyen_Mother[j])
    
    if (LAGi>182) { LAGi=-183+(LAGi-182) }
    if (LAGi < -182) {LAGi=183+(LAGi+182)}
    vec=c(vec,LAGi)
    
  }
  TAB=cbind(TAB,vec)
}


Estimate_Rank=c()
Estimate_Parite=c()
Estimate_Sex=c()
Pvalue_Rank=c()
Pvalue_Parite=c()
Pvalue_Sex=c()

for (j in (ncol(TAB)-999):ncol(TAB)) {
  print(j)
  
  MOD=blmer(abs(TAB[,j]) ~ Relative_Rank + Parity + Sex + Troop +
              (1|Mother) +  (1|Year_Birth), data=TAB)
  S=summary(MOD)
  Estimate_Rank=c(Estimate_Rank,S$coefficients[2,1])
  Pvalue_Rank=c(Pvalue_Rank,Anova(MOD)[1,3])
  
  Estimate_Parite=c(Estimate_Parite,S$coefficients[3,1])
  Pvalue_Parite=c(Pvalue_Parite,Anova(MOD)[2,3])
  
  Estimate_Sex=c(Estimate_Sex,S$coefficients[4,1])
  Pvalue_Sex=c(Pvalue_Sex,Anova(MOD)[3,3])
  
  
}
#here we can summarize the Pvalues and estimates of our fixed effects. 


############################################################################################################
######### 3. Birth timings favouring future reproduction intensify mother-offspring conflict
############################################################################################################

########## Model 5: influence of birth dates on suckling probability

TAB='suckling & carrying' #open the excel file. 
TAB=TAB[which(TAB$AgeMonths>1),]
TAB=TAB[which(TAB$AgeMonths<18),]

TAB$Suckling_Scan=as.factor(as.character(TAB$Suckling_Scan)) #
TAB$Focal.ID=factor(TAB$Focal.ID)
TAB$Year_Focal=as.factor(as.character(TAB$Year_Focal))
TAB$ID_Focal=as.factor(as.character(TAB$ID_Focal))
TAB$AgeMonths=as.numeric(scale(TAB$AgeMonths))
TAB$Relative.rank=as.numeric(scale(TAB$Relative.rank))

#First: determine the best age effects on suckling probabilities : 

GAMM_BINOMIAL=gam(Suckling_Scan ~ s(AgeMonths) + s(ID_Focal,bs='re') + s(Focal.ID,bs='re'),
                    family='binomial',data=TAB)
AIC(GAMM_BINOMIAL) #
GAMM_BINOMIAL1=gam(Suckling_Scan ~ AgeMonths + s(ID_Focal,bs='re') + s(Focal.ID,bs='re'),
                     family='binomial',data=TAB)
AIC(GAMM_BINOMIAL1) #
GAMM_BINOMIAL2=gam(Suckling_Scan ~ poly(AgeMonths,2) + s(ID_Focal,bs='re') + s(Focal.ID,bs='re'),
                     family='binomial',data=TAB)
AIC(GAMM_BINOMIAL2) # 
GAMM_BINOMIAL3=gam(Suckling_Scan ~ poly(AgeMonths,3) + s(ID_Focal,bs='re') + s(Focal.ID,bs='re'),
                     family='binomial',data=TAB)
AIC(GAMM_BINOMIAL3) 

#check the best age effect with linear models (once we know it's not a gam)
MOD0=bglmer(Suckling_Scan  ~ AgeMonths + 
              (1|ID_Focal) + (1|Focal.ID), 
            glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 10000)),
            family='binomial',data=TAB)
MOD1=bglmer(Suckling_Scan  ~ poly(AgeMonths,2) + 
              (1|ID_Focal) + (1|Focal.ID), 
            glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 10000)),
            family='binomial',data=TAB)
MOD2=bglmer(Suckling_Scan  ~ poly(AgeMonths,3) + 
              (1|ID_Focal) + (1|Focal.ID), 
            glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 10000)),
            family='binomial',data=TAB)
AIC(MOD0,MOD1,MOD2)



#selection of the best phase for the sine term of the birth date effect 
MOD0=bglmer(Suckling_Scan  ~ AgeMonths + 
              sin(DOB_Radian + 0*pi/12) + 
              Sex + Parity + Relative.rank + Troop + Year_Focal +
              (1|ID_Focal) + (1|Focal.ID), 
            glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 10000)),
            family='binomial',data=TAB)
MOD1=bglmer(Suckling_Scan  ~ AgeMonths + 
              sin(DOB_Radian + 1*pi/12) + 
              Sex + Parity + Relative.rank + Troop + Year_Focal +
              (1|ID_Focal) + (1|Focal.ID), 
            glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 10000)),
            family='binomial',data=TAB)
MOD2=bglmer(Suckling_Scan  ~ AgeMonths + 
              sin(DOB_Radian + 2*pi/12) + 
              Sex + Parity + Relative.rank + Troop + Year_Focal +
              (1|ID_Focal) + (1|Focal.ID), 
            glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 10000)),
            family='binomial',data=TAB)
MOD3=bglmer(Suckling_Scan  ~ AgeMonths + 
              sin(DOB_Radian + 3*pi/12) + 
              Sex + Parity + Relative.rank + Troop + Year_Focal +
              (1|ID_Focal) + (1|Focal.ID), 
            glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 10000)),
            family='binomial',data=TAB)
MOD4=bglmer(Suckling_Scan  ~ AgeMonths + 
              sin(DOB_Radian + 4*pi/12) + 
              Sex + Parity + Relative.rank + Troop + Year_Focal +
              (1|ID_Focal) + (1|Focal.ID), 
            glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 10000)),
            family='binomial',data=TAB)
MOD5=bglmer(Suckling_Scan  ~ AgeMonths + 
              sin(DOB_Radian + 5*pi/12) + 
              Sex + Parity + Relative.rank + Troop + Year_Focal +
              (1|ID_Focal) + (1|Focal.ID), 
            glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 10000)),
            family='binomial',data=TAB)
MOD6=bglmer(Suckling_Scan  ~ AgeMonths + 
              sin(DOB_Radian + 6*pi/12) + 
              Sex + Parity + Relative.rank + Troop + Year_Focal +
              (1|ID_Focal) + (1|Focal.ID), 
            glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 10000)),
            family='binomial',data=TAB)
MOD7=bglmer(Suckling_Scan  ~ AgeMonths + 
              sin(DOB_Radian + 7*pi/12) + 
              Sex + Parity + Relative.rank + Troop + Year_Focal +
              (1|ID_Focal) + (1|Focal.ID), 
            glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 10000)),
            family='binomial',data=TAB)
MOD8=bglmer(Suckling_Scan  ~ AgeMonths + 
              sin(DOB_Radian + 8*pi/12) + 
              Sex + Parity + Relative.rank + Troop + Year_Focal +
              (1|ID_Focal) + (1|Focal.ID), 
            glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 10000)),
            family='binomial',data=TAB)
MOD9=bglmer(Suckling_Scan  ~ AgeMonths + 
              sin(DOB_Radian + 9*pi/12) + 
              Sex + Parity + Relative.rank + Troop + Year_Focal +
              (1|ID_Focal) + (1|Focal.ID), 
            glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 10000)),
            family='binomial',data=TAB)
MOD10=bglmer(Suckling_Scan  ~ AgeMonths + 
               sin(DOB_Radian + 10*pi/12) + 
               Sex + Parity + Relative.rank + Troop + Year_Focal +
               (1|ID_Focal) + (1|Focal.ID), 
             glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 10000)),
             family='binomial',data=TAB)
MOD11=bglmer(Suckling_Scan  ~ AgeMonths + 
               sin(DOB_Radian + 11*pi/12) + 
               Sex + Parity + Relative.rank + Troop + Year_Focal +
               (1|ID_Focal) + (1|Focal.ID), 
             glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 10000)),
             family='binomial',data=TAB)
AIC(MOD0,MOD1,MOD2,MOD3,MOD4,MOD5,MOD6,MOD7,MOD8,MOD9,MOD10,MOD11)
#we select here the MOD9
summary(MOD9)
Anova(MOD9)
confint(MOD9, method='Wald')

simulationOutput <- simulateResiduals(fittedModel = MOD9, plot = T)


#randomization to control for the uncertainty of birth dates on the fixed effect : sine(DOB_Radian)
for (i in 1:1000) {
  print(i)
  vec=c()
  for (j in 1:nrow(TAB)) {
    ff=runif(1,TAB$DOB_Radian[j]-TAB$Uncertainty_Birth_DayRad[j],TAB$DOB_Radian[j]+TAB$Uncertainty_Birth_DayRad[j])
    
    if (ff>2*pi) {ff=ff-2*pi}
    if (ff<0) {ff=ff+2*pi}
    vec[j]=ff
  }
  TAB=cbind(TAB,vec)
}
Pvalue=c()
for (j in (ncol(TAB)-999):ncol(TAB)) {
  print(j)
  
  MOD=bglmer(Suckling_Scan  ~ AgeMonths + 
               sin(TAB[,j] + 9*pi/12) + AgeMonths:sin(TAB[,j] + 9*pi/12) +
               Sex + Parity + Relative.rank + Troop + Year_Focal +
               (1|ID_Focal) + (1|Focal.ID), 
             glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 10000)),
             family='binomial',data=TAB)
  Pvalue=c(Pvalue,Anova(MOD)[8,3]) 
}
summary(Pvalue)





#note here that, for Models 5-7, we similarly ran models with 12 different phase values to test the best effect
#of the sine term of the observation date, as follow : (we won't show these models for mother carrying nor tantrum, but the method used is the same as for suckling)
MOD0=bglmer(Suckling_Scan  ~ AgeMonths + 
              sin(Jour_Obs_Rad + 0*pi/12) + 
              Sex + Parity + Relative.rank + Troop + Year_Focal +
              (1|ID_Focal) + (1|Focal.ID), 
            glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 10000)),
            family='binomial',data=TAB)
MOD1=bglmer(Suckling_Scan  ~ AgeMonths + 
              sin(Jour_Obs_Rad + 1*pi/12) + 
              Sex + Parity + Relative.rank + Troop + Year_Focal +
              (1|ID_Focal) + (1|Focal.ID), 
            glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 10000)),
            family='binomial',data=TAB)
MOD2=bglmer(Suckling_Scan  ~ AgeMonths + 
              sin(Jour_Obs_Rad + 2*pi/12) + 
              Sex + Parity + Relative.rank + Troop + Year_Focal +
              (1|ID_Focal) + (1|Focal.ID), 
            glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 10000)),
            family='binomial',data=TAB)
MOD3=bglmer(Suckling_Scan  ~ AgeMonths + 
              sin(Jour_Obs_Rad + 3*pi/12) + 
              Sex + Parity + Relative.rank + Troop + Year_Focal +
              (1|ID_Focal) + (1|Focal.ID), 
            glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 10000)),
            family='binomial',data=TAB)
MOD4=bglmer(Suckling_Scan  ~ AgeMonths + 
              sin(Jour_Obs_Rad + 4*pi/12) + 
              Sex + Parity + Relative.rank + Troop + Year_Focal +
              (1|ID_Focal) + (1|Focal.ID), 
            glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 10000)),
            family='binomial',data=TAB)
MOD5=bglmer(Suckling_Scan  ~ AgeMonths + 
              sin(Jour_Obs_Rad + 5*pi/12) + 
              Sex + Parity + Relative.rank + Troop + Year_Focal +
              (1|ID_Focal) + (1|Focal.ID), 
            glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 10000)),
            family='binomial',data=TAB)
MOD6=bglmer(Suckling_Scan  ~ AgeMonths + 
              sin(Jour_Obs_Rad + 6*pi/12) + 
              Sex + Parity + Relative.rank + Troop + Year_Focal +
              (1|ID_Focal) + (1|Focal.ID), 
            glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 10000)),
            family='binomial',data=TAB)
MOD7=bglmer(Suckling_Scan  ~ AgeMonths + 
              sin(Jour_Obs_Rad + 7*pi/12) + 
              Sex + Parity + Relative.rank + Troop + Year_Focal +
              (1|ID_Focal) + (1|Focal.ID), 
            glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 10000)),
            family='binomial',data=TAB)
MOD8=bglmer(Suckling_Scan  ~ AgeMonths + 
              sin(Jour_Obs_Rad + 8*pi/12) + 
              Sex + Parity + Relative.rank + Troop + Year_Focal +
              (1|ID_Focal) + (1|Focal.ID), 
            glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 10000)),
            family='binomial',data=TAB)
MOD9=bglmer(Suckling_Scan  ~ AgeMonths + 
              sin(Jour_Obs_Rad + 9*pi/12) + 
              Sex + Parity + Relative.rank + Troop + Year_Focal +
              (1|ID_Focal) + (1|Focal.ID), 
            glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 10000)),
            family='binomial',data=TAB)
MOD10=bglmer(Suckling_Scan  ~ AgeMonths + 
               sin(Jour_Obs_Rad + 10*pi/12) + 
               Sex + Parity + Relative.rank + Troop + Year_Focal +
               (1|ID_Focal) + (1|Focal.ID), 
             glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 10000)),
             family='binomial',data=TAB)
MOD11=bglmer(Suckling_Scan  ~ AgeMonths + 
               sin(Jour_Obs_Rad + 11*pi/12) + 
               Sex + Parity + Relative.rank + Troop + Year_Focal +
               (1|ID_Focal) + (1|Focal.ID), 
             glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 10000)),
             family='binomial',data=TAB)

AIC(MOD0,MOD1,MOD2,MOD3,MOD4,MOD5,MOD6,MOD7,MOD8,MOD9,MOD10,MOD11)

summary(MOD7)
Anova(MOD7)
confint(MOD7,method='Wald')

simulationOutput <- simulateResiduals(fittedModel = MOD7, plot = T)



########## Model 6: influence of birth dates on mother carrying probability

TAB='suckling & carrying' #open the excel file. 
TAB=TAB[which(is.na(TAB$Travel_Mother_Scan)==FALSE),] #only selecting the travelling scan here. 
TAB=TAB[which(TAB$AgeMonths>1),]
TAB=TAB[which(TAB$AgeMonths<18),]

TAB$Travel_Mother_Scan=as.factor(as.character(TAB$Travel_Mother_Scan))
TAB$Focal.ID=factor(TAB$Focal.ID)
TAB$Year_Focal=as.factor(as.character(TAB$Year_Focal))
TAB$ID_Focal=as.factor(as.character(TAB$ID_Focal))
TAB$AgeMonths=as.numeric(scale(TAB$AgeMonths))
TAB$Relative.rank=as.numeric(scale(TAB$Relative.rank))

#First: determine the best age effects on suckling probabilities : 

GAMM_BINOMIAL=gam(Travel_Mother_Scan ~ s(AgeMonths) + s(ID_Focal,bs='re') + s(Focal.ID,bs='re'),
                  family='binomial',data=TAB)
AIC(GAMM_BINOMIAL) #
GAMM_BINOMIAL1=gam(Travel_Mother_Scan ~ AgeMonths + s(ID_Focal,bs='re') + s(Focal.ID,bs='re'),
                     family='binomial',data=TAB)
AIC(GAMM_BINOMIAL1) #
GAMM_BINOMIAL2=gam(Travel_Mother_Scan ~ poly(AgeMonths,2) + s(ID_Focal,bs='re') + s(Focal.ID,bs='re'),
                     family='binomial',data=TAB)
AIC(GAMM_BINOMIAL2) # 
GAMM_BINOMIAL3=gam(Travel_Mother_Scan ~ poly(AgeMonths,3) + s(ID_Focal,bs='re') + s(Focal.ID,bs='re'),
                     family='binomial',data=TAB)
AIC(GAMM_BINOMIAL3) 

#check the best age effect with linear models (once we know it's not a gam)
MOD0=bglmer(Travel_Mother_Scan  ~ AgeMonths + 
              (1|ID_Focal) + (1|Focal.ID), 
            glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 10000)),
            family='binomial',data=TAB)
MOD1=bglmer(Travel_Mother_Scan  ~ poly(AgeMonths,2) + 
              (1|ID_Focal) + (1|Focal.ID), 
            glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 10000)),
            family='binomial',data=TAB)
MOD2=bglmer(Travel_Mother_Scan  ~ poly(AgeMonths,3) + 
              (1|ID_Focal) + (1|Focal.ID), 
            glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 10000)),
            family='binomial',data=TAB)
AIC(MOD0,MOD1,MOD2)


#we found that the linear effect of age is the best. 
#now, we look for the best phase of the sine term of the birth date fixed effect: 
MOD0=bglmer(Travel_Mother_Scan  ~ AgeMonths + sin(DOB_Radian) + Troop + Year_Focal + 
              Sex + Parity + Relative.rank +
              (1|ID_Focal) +  (1|Focal.ID), 
            glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 10000)),
            family='binomial',data=TAB)
MOD1=bglmer(Travel_Mother_Scan  ~ AgeMonths + sin(DOB_Radian + 1*pi/12) + Troop + Year_Focal + 
              Sex + Parity + Relative.rank +
              (1|ID_Focal) +  (1|Focal.ID), 
            glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 10000)),
            family='binomial',data=TAB)
MOD2=bglmer(Travel_Mother_Scan  ~ AgeMonths + sin(DOB_Radian + 2*pi/12) + Troop + Year_Focal + 
              Sex + Parity + Relative.rank +
              (1|ID_Focal) +  (1|Focal.ID), 
            glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 10000)),
            family='binomial',data=TAB)
MOD3=bglmer(Travel_Mother_Scan  ~ AgeMonths + sin(DOB_Radian + 3*pi/12) + Troop + Year_Focal + 
              Sex + Parity + Relative.rank +
              (1|ID_Focal) +  (1|Focal.ID), 
            glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 10000)),
            family='binomial',data=TAB)
MOD4=bglmer(Travel_Mother_Scan  ~ AgeMonths + sin(DOB_Radian + 4*pi/12) + Troop + Year_Focal + 
              Sex + Parity + Relative.rank +
              (1|ID_Focal) +  (1|Focal.ID), 
            glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 10000)),
            family='binomial',data=TAB)
MOD5=bglmer(Travel_Mother_Scan  ~ AgeMonths + sin(DOB_Radian + 5*pi/12) + Troop + Year_Focal + 
              Sex + Parity + Relative.rank +
              (1|ID_Focal) +  (1|Focal.ID), 
            glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 10000)),
            family='binomial',data=TAB)
MOD6=bglmer(Travel_Mother_Scan  ~ AgeMonths + sin(DOB_Radian + 6*pi/12) + Troop + Year_Focal + 
              Sex + Parity + Relative.rank +
              (1|ID_Focal) +  (1|Focal.ID), 
            glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 10000)),
            family='binomial',data=TAB)
MOD7=bglmer(Travel_Mother_Scan  ~ AgeMonths + sin(DOB_Radian + 7*pi/12) + Troop + Year_Focal + 
              Sex + Parity + Relative.rank +
              (1|ID_Focal) +  (1|Focal.ID), 
            glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 10000)),
            family='binomial',data=TAB)
MOD8=bglmer(Travel_Mother_Scan  ~ AgeMonths + sin(DOB_Radian + 8*pi/12) + Troop + Year_Focal + 
              Sex + Parity + Relative.rank +
              (1|ID_Focal) +  (1|Focal.ID), 
            glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 10000)),
            family='binomial',data=TAB)
MOD9=bglmer(Travel_Mother_Scan  ~ AgeMonths + sin(DOB_Radian + 9*pi/12) + Troop + Year_Focal + 
              Sex + Parity + Relative.rank +
              (1|ID_Focal) +  (1|Focal.ID), 
            glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 10000)),
            family='binomial',data=TAB)
MOD10=bglmer(Travel_Mother_Scan  ~ AgeMonths + sin(DOB_Radian + 10*pi/12) + Troop + Year_Focal + 
               Sex + Parity + Relative.rank +
               (1|ID_Focal) +  (1|Focal.ID), 
             glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 10000)),
             family='binomial',data=TAB)
MOD11=bglmer(Travel_Mother_Scan  ~ AgeMonths + sin(DOB_Radian + 11*pi/12) + Troop + Year_Focal + 
               Sex + Parity + Relative.rank +
               (1|ID_Focal) +  (1|Focal.ID), 
             glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 10000)),
             family='binomial',data=TAB)
AIC(MOD0,MOD1,MOD2,MOD3,MOD4,MOD5,MOD6,MOD7,MOD8,MOD9,MOD10,MOD11)
#best model is with the phase=0
summary(MOD0)
Anova(MOD0)
confint(MOD0,method='Wald')
simulationOutput <- simulateResiduals(fittedModel = MOD0, plot = T)


#now we control for birth date uncertaintues affectinf the sine term fixed effect : 


for (i in 1:1000) {
  print(i)
  vec=c()
  for (j in 1:nrow(TAB)) {
    ff=runif(1,TAB$DOB_Radian[j]-TAB$Uncertainty_Birth_DayRad[j],TAB$DOB_Radian[j]+TAB$Uncertainty_Birth_DayRad[j])
    
    if (ff>2*pi) {ff=ff-2*pi}
    if (ff<0) {ff=ff+2*pi}
    vec[j]=ff
  }
  TAB=cbind(TAB,vec)
}

Pvalue=c()
for (j in (ncol(TAB)-999):ncol(TAB)) {
  print(j)
  
  MOD=bglmer(Travel_Mother_Scan  ~ AgeMonths + sin(TAB[,j] + 0*pi/12) +
               Sex + Parity + Relative.rank + AgeMonths:sin(TAB[,j] + 0*pi/12) +
               Troop + Year_Focal +
               (1|ID_Focal) +  (1|Focal.ID), 
             glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 10000)),
             family='binomial',data=TAB)
  Pvalue=c(Pvalue,Anova(MOD)[8,3]) 
  
}
summary(Pvalue)


########## Model 7: influence of birth dates on tantrum probability
TAB='Tantrum' #open excel file

TAB=TAB[which(TAB$AgeMonths<18),]
TAB=TAB[which(TAB$AgeMonths>1),]

TAB$TANTRUM=as.factor(as.character(TAB$TANTRUM)) #
TAB$Focal.ID=factor(TAB$Focal.ID)
TAB$Year_Focal=as.factor(as.character(TAB$Year_Focal))
TAB$T_Obs_Real=as.numeric(scale(TAB$T_Obs_Real))
TAB$Relative_Rank=as.numeric(scale(TAB$Relative_Rank))
TAB$AgeMonths=as.numeric(scale(TAB$AgeMonths))

#first, we assess the best age effect 
MOD_S=gam(TANTRUM ~ s(AgeMonths) +
            T_Obs_Real + s(Focal.ID,bs='re'), 
          family='binomial',data=TAB)
MOD_1=gam(TANTRUM ~ AgeMonths +
            T_Obs_Real + s(Focal.ID,bs='re'), 
          family='binomial',data=TAB)
MOD_2=gam(TANTRUM ~ poly(AgeMonths,2) +
              T_Obs_Real + s(Focal.ID,bs='re'), 
            family='binomial',data=TAB)
MOD_3=gam(TANTRUM ~ poly(AgeMonths,3) +
            T_Obs_Real + s(Focal.ID,bs='re'), 
          family='binomial',data=TAB)
AIC(MOD_S,MOD_1,MOD_2,MOD_3)
#once we know we won't chose the gam model, we ran 3 linear models, and chose the one minimizing AIC: 
MOD1=bglmer(TANTRUM ~ AgeMonths + 
              T_Obs_Real + (1|Focal.ID), 
            glmerControl(optCtrl = list(maxfun = 20000)),
            family='binomial',data=TAB)
MOD2=bglmer(TANTRUM ~ poly(AgeMonths,2) + 
              T_Obs_Real + (1|Focal.ID), 
            glmerControl(optCtrl = list(maxfun = 20000)),
            family='binomial',data=TAB)
MOD3=bglmer(TANTRUM ~ poly(AgeMonths,3) + 
              T_Obs_Real + (1|Focal.ID), 
            glmerControl(optCtrl = list(maxfun = 20000)),
            family='binomial',data=TAB)
AIC(MOD1,MOD2,MOD3)
#we select the polynomial of degree 2

#we chose the best phase of the sine term of the birth date effect:
MOD0=bglmer(TANTRUM ~ poly(AgeMonths,2) + 
              sin(DOB_Radian + 0*pi/12) + Troop + Year_Focal +
              Sex + Parity + Relative_Rank +
              T_Obs_Real + (1|Focal.ID), 
            glmerControl(optCtrl = list(maxfun = 20000)),
            family='binomial',data=TAB)
MOD1=bglmer(TANTRUM ~ poly(AgeMonths,2) + 
              sin(DOB_Radian + 1*pi/12) + Troop + Year_Focal +
              Sex + Parity + Relative_Rank +
              T_Obs_Real + (1|Focal.ID), 
            glmerControl(optCtrl = list(maxfun = 20000)),
            family='binomial',data=TAB)
MOD2=bglmer(TANTRUM ~ poly(AgeMonths,2) + 
                 sin(DOB_Radian + 2*pi/12) + Troop + Year_Focal +
                 Sex + Parity + Relative_Rank +
                 T_Obs_Real + (1|Focal.ID), 
               glmerControl(optCtrl = list(maxfun = 20000)),
               family='binomial',data=TAB)
MOD3=bglmer(TANTRUM ~ poly(AgeMonths,2) + 
              sin(DOB_Radian + 3*pi/12) + Troop + Year_Focal +
              Sex + Parity + Relative_Rank +
              T_Obs_Real + (1|Focal.ID), 
            glmerControl(optCtrl = list(maxfun = 20000)),
            family='binomial',data=TAB)
MOD4=bglmer(TANTRUM ~ poly(AgeMonths,2) + 
              sin(DOB_Radian + 4*pi/12) + Troop + Year_Focal +
              Sex + Parity + Relative_Rank +
              T_Obs_Real + (1|Focal.ID), 
            glmerControl(optCtrl = list(maxfun = 20000)),
            family='binomial',data=TAB)
MOD5=bglmer(TANTRUM ~ poly(AgeMonths,2) + 
              sin(DOB_Radian + 5*pi/12) + Troop + Year_Focal +
              Sex + Parity + Relative_Rank +
              T_Obs_Real + (1|Focal.ID), 
            glmerControl(optCtrl = list(maxfun = 20000)),
            family='binomial',data=TAB)
MOD6=bglmer(TANTRUM ~ poly(AgeMonths,2) + 
              sin(DOB_Radian + 6*pi/12) + Troop + Year_Focal +
              Sex + Parity + Relative_Rank +
              T_Obs_Real + (1|Focal.ID), 
            glmerControl(optCtrl = list(maxfun = 20000)),
            family='binomial',data=TAB)
MOD7=bglmer(TANTRUM ~ poly(AgeMonths,2) + 
              sin(DOB_Radian + 7*pi/12) + Troop + Year_Focal +
              Sex + Parity + Relative_Rank +
              T_Obs_Real + (1|Focal.ID), 
            glmerControl(optCtrl = list(maxfun = 20000)),
            family='binomial',data=TAB)
MOD8=bglmer(TANTRUM ~ poly(AgeMonths,2) + 
              sin(DOB_Radian + 8*pi/12) + Troop + Year_Focal +
              Sex + Parity + Relative_Rank +
              T_Obs_Real + (1|Focal.ID), 
            glmerControl(optCtrl = list(maxfun = 20000)),
            family='binomial',data=TAB)
MOD9=bglmer(TANTRUM ~ poly(AgeMonths,2) + 
              sin(DOB_Radian + 9*pi/12) + Troop + Year_Focal +
              Sex + Parity + Relative_Rank +
              T_Obs_Real + (1|Focal.ID), 
            glmerControl(optCtrl = list(maxfun = 20000)),
            family='binomial',data=TAB)
MOD10=bglmer(TANTRUM ~ poly(AgeMonths,2) + 
              sin(DOB_Radian + 10*pi/12) + Troop + Year_Focal +
              Sex + Parity + Relative_Rank +
              T_Obs_Real + (1|Focal.ID), 
            glmerControl(optCtrl = list(maxfun = 20000)),
            family='binomial',data=TAB)
MOD11=bglmer(TANTRUM ~ poly(AgeMonths,2) + 
              sin(DOB_Radian + 11*pi/12) + Troop + Year_Focal +
              Sex + Parity + Relative_Rank +
              T_Obs_Real + (1|Focal.ID), 
            glmerControl(optCtrl = list(maxfun = 20000)),
            family='binomial',data=TAB)
AIC(MOD0,MOD1,MOD2,MOD3,MOD4,MOD5,MOD6,MOD7,MOD8,MOD9,MOD10,MOD11)
#best model is the one with phase = 2*pi/12

summary(MOD2)
Anova(MOD2)
confint(MOD2,method='Wald')
simulationOutput <- simulateResiduals(fittedModel = MOD2, plot = T)


#Lastly, we control for birth date uncertainties likely to affect our fixed effect: sine term of the birth date
for (i in 1:1000) {
  print(i)
  vec=c()
  for (j in 1:nrow(TAB)) {
    ff=runif(1,TAB$DOB_Radian[j]-TAB$Uncertainty_Birth_DayRad[j],TAB$DOB_Radian[j]+TAB$Uncertainty_Birth_DayRad[j])
    
    if (ff>2*pi) {ff=ff-2*pi}
    if (ff<0) {ff=ff+2*pi}
    vec[j]=ff
  }
  TAB=cbind(TAB,vec)
}
Pvalue=c()
for (j in (ncol(TAB)-999):ncol(TAB)) {
  print(j)
  
  MOD_TEST=bglmer(TANTRUM_conservateur30s_V2 ~ poly(AgeMonths,2) + 
                    sin(TAB[,j] + 2*pi/12) + Troop + Year_Focal +
                    Sex + Parity + Relative_Rank +
                    T_Obs_Real + (1|Focal.ID), 
                  glmerControl(optCtrl = list(maxfun = 20000)),
                  family='binomial',data=TAB)
  Pvalue=c(Pvalue,Anova(MOD_TEST)[2,3]) 
  
} 
summary(Pvalue)
wilcox.test(Pvalue, conf.int = TRUE, conf.level = 0.95) #95% confidence interval of the 1000 randomized pvalues



