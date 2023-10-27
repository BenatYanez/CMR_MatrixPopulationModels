library(dplyr)
library(tidyr)
library(marked)
library(ggplot2)
library(R2ucare)

longdata <- read.table("CMR practical/sparrowrecap.txt",header=T, sep="\t")

length(unique(longdata$id)) #Several individuals appear multiple times becasue of the recapture, there are 1062 unique individuals
table(longdata$sex) #Equal number of males and females
#Marked requires specific dataset arrangement

temp <- longdata[,1:2] # take the first two columns, id and year and put into a temporary dataframe
temp$detect <- 1 # add column for detection (all 1s because these represent captures) 

temp <- temp %>%
  # remove duplicates, which may occur when individuals are caught multiple times in an sampling event
  distinct() %>%
  # spread out data. The fill = 0 adds rows for combinations of id and year where individuals were not observed
  spread(year, detect, fill = 0) %>% 
  # for every individual....
  group_by(id) %>%
  # paste together 0's and 1's using unite()
  # here we are pasting the strings together from the second column (first capture event)
  # to the last capture event ("tail(names(.),1)")
  # use sep="" so there are no characters separating 0's and 1's
  unite("ch", 2:tail(names(.),1), sep = "")

sparrow <- as.data.frame(temp)
#Add bakc information by matching the ids in sparrow to information related to teh id in longdata
sparrow$island <- longdata$island[match(sparrow$id, longdata$id)] 
sparrow$sex <- as.factor(longdata$sex[match(sparrow$id, longdata$id)])
sparrow <- droplevels(subset(sparrow, select = -id)) # remove id column so capture histories appear in first column

#The CJS model estimates survival and detection probability. Using a linear model with a logit scale (bound between 0, 1)
#USes information from capture history where individuals were marked, not observed and the seen at a later date to estimate detection probabilities
#So comparing individuals with history of 101 vs 111 is informative about detection probability in second sampling event
mod1 <- crm(sparrow)
mod1 #Estimate of survival= 0.071, Estimate of detection prob=0.322
mod1 <- cjs.hessian(mod1) # refit model with precision estimates
mod1 #the results are on the logit scale, need to transform them back

mod1$results$reals
plogis(mod1$results$beta$Phi)#0.5177
plogis(mod1$results$beta$p)#0.5799

predict(mod1, newdata=data.frame(sex = c('Female', 'Male')), se=T)

#The model assumes an equal time between each capture event. The assumption can be relaxed by including a vector of time intervals
mod2 <- crm(sparrow, time.intervals = c(1,2,1,1,1,1,1,3,4))
mod2$results$reals
#QUESTIONS: These models assumed constant survival rates and detection probabilities.
#Is this a realistic assumption for this system? 
#If the capture effort is the same every time it is done, we should expect detection probabilities to be the same, or at least to vary very little
#The detection probability might differ between islands
#The survival rate might not be constant and depend on factors that change every year, due to environmental stochasticity
#What term might you include in the model next, and why?
#we might allow detection to change between islands.

sparrow.proc <- process.data(sparrow) # built in function for data processing
str(sparrow.proc)

sparrow.ddl <- make.design.data(sparrow.proc) # built in function for building design matrix 
# specify model formulation: capture probability depends on island
p.island <- list(formula=~island) #Determines which factor p will depend one, could also do =sex to differentiate capture probabilities of males and females
mod3 <- crm(sparrow.proc, 
            sparrow.ddl, 
            model.parameters = list(p = p.island), 
            accumulate=FALSE, hessian = TRUE)
mod3$results$reals
#QUESTION Does it look like detection probability varies among islands? 
#There is variation in detection between Myken and the rest of the islands and Kvaroy and the rest
#Which island has the lowest detection probability? Why might this be?
#Myken has the lowest detection probability
#How might you compare this model with our simpler model above to see if it is a better fit?
#Looking at teh AIC score of each of the models and seeing which is lowest
(mod3$results$AIC) #Better fit than the model with combined detection probabiity
(mod1$results$AIC)
#Look at whether survival changes depending on island
mod4 <- crm(sparrow.proc, 
            sparrow.ddl, 
            model.parameters = list(Phi = p.island, p=p.island), 
            accumulate=FALSE, hessian = TRUE)
mod4$results$reals
(mod4$results$AIC) #Equally good fit as the model with only the detection probability of the islands varying, so the simpler model is favoured (model3)
(mod3$results$AIC)#There is no evidence that survival varies between islands 
#The expectation is that detection is different between islands just due to logistics of the methodology, so we account for this before looking at the biologicallly interesting question of differences in survival
#Does survival vary between sexes?

fit.models <- function() {
  Phi.dot <- list(formula=~1) # constant survival
  Phi.sex <- list(formula=~sex) # survival differs between sexes
  Phi.island <- list(formula=~island) # survival differs between islands
  Phi.sex.island <- list(formula=~sex+island) # survival differs between sexes and islands
  p.dot <- list(formula=~1) # constant detection
  p.sex <- list(formula=~sex) # detection probability differs between sexes
  p.island <- list(formula=~island) # detection probability differs between islands
  p.sex.island <- list(formula=~sex+island) # detection probability differs between sexes and islands
  cml <- create.model.list(c("Phi","p"))
  results <- crm.wrapper(cml, data=sparrow.proc, ddl=sparrow.ddl,
                         external=FALSE, accumulate=FALSE, hessian=TRUE) #Automate running several crm models
  return(results)
}
sparrow.models <- fit.models() # run function 
sparrow.models
#The top model (lowest AIC and simplest model) is that where survival is constant and detection probability differs by island.
mod5 <- sparrow.models[[2]]
ggplot(mod5$results$reals$p, aes(island, estimate, ymin=lcl, ymax=ucl)) + 
  geom_errorbar(width=0.2) + geom_point() + ylim(0,1)
#We can also extract sex and island differences in survival from similarly supported models and plot them to convince ourselves there aren't any big differences
mod6 <- sparrow.models[[10]]
ggplot(mod6$results$reals$Phi, aes(sex, estimate, ymin=lcl, ymax=ucl)) + 
  geom_errorbar(width=0.2) + geom_point() + ylim(0,1)

##We might want to test whether survival changes depending on a factor that differs over time, like weather. 
sparrow.ddl$Phi$cold <- "Cold" # new column
sparrow.ddl$Phi$cold[sparrow.ddl$Phi$time==2 | sparrow.ddl$Phi$time==5 | sparrow.ddl$Phi$time==8] <- "VeryCold" # very cold winters between capture events 2 and 3, 5 and 6, and 8 and 9
Phi.cold <- list(formula=~cold) 
mod8 <- crm(sparrow.proc, 
            sparrow.ddl, 
            model.parameters = list(Phi = Phi.cold, 
                                    p = p.island), 
            accumulate=FALSE, hessian = TRUE)
mod8$results$reals
(mod8$results$AIC)#The model including cold has a higher AIC than the simplest model with the best fit, so there is no evidence that survival varies depending on made-up variable
(mod5$results$AIC)
#Test whether survival and detection prbability varies each year of the study using the built-in time variable
fit.models2 <- function() {
  Phi.dot <- list(formula=~1) # constant survival
  Phi.time <- list(formula=~time) # survival varies over time
  p.island <- list(formula=~island) # detection probability differs between islands
  p.time <- list(formula=~time) # detection probability varies over time
  p.island.time <- list(formula=~island+time) # detection probability varies over time
  cml <- create.model.list(c("Phi","p"))
  results <- crm.wrapper(cml, data=sparrow.proc, ddl=sparrow.ddl,
                         external=FALSE, accumulate=FALSE, hessian=TRUE)
  return(results)
}
sparrow.models2 <- fit.models2() # run function 
sparrow.models2 # display model table
mod9 <- sparrow.models2[[2]] #Seems taht the model with detection probabilities varying over time fits better than one where it does not

(mod9$results$AIC)
ggplot(mod9$results$reals$p, aes(time, estimate, ymin=lcl, ymax=ucl, col=island)) + 
  geom_errorbar(width=0) + geom_point() + ylim(0,1)

#Goodness of fit test
#First reformat dat into a matrix
sparrow.gof <- sparrow$ch %>%
  strsplit('') %>%
  sapply(`[`) %>%
  t() %>%
  unlist() %>%
  as.numeric %>%
  matrix(nrow = nrow(sparrow))
#R2ucare does 3 tests
#Test 1: the overall test. Overall, is there evidence that animals have equal detection probabilities and equal survival?
overall_CJS(sparrow.gof, rep(1,nrow(sparrow)))
#P-value is not significant so we fail to reject the null. No strong evidence for overall lack of fit
#Test 2: Does recapture depend on when an animal was first marked? (Tests the equal detection assumption)
  #Test 2 CT: Is there a difference in p at t+1 between those captured and not captured at t (when animals are known to be alive because are captured later in the study)?
test2ct <- test2ct(sparrow.gof, rep(1,nrow(sparrow))) 
test2ct
#We fail to reject the null so no evidence of a problem with the equal detection assumption
  #Test 2 CL: Is there a difference in the expected time of next recapture between individuals captured and not captured at t when animals are known to be alive?
test2cl <- test2cl(sparrow.gof, rep(1,nrow(sparrow)))
test2cl
#We fail to reject the null so no evidence of a problem with the equal detection assumption

#Test 3: Does marking affect survival? (Tests the equal survival assumption)
  #Test 3 SR: Do individuals with previous marks have different survival rates than first-time captures?
test3sr <- test3sr(sparrow.gof, rep(1,nrow(sparrow)))
test3sr
 #Test 3 SM: For animals seen again, does when they are recaptured depend on whether they were marked on or before t?
test3sm <- test3sm(sparrow.gof, rep(1,nrow(sparrow)))
test3sm
#We fail to reject the null so no evidence of a problem with the unequal survival assumption

