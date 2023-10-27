library(ggplot2)
library(popbio)

survival <- data.frame(stage=factor(c('Juvenile','Yearling','Adult'), levels=c('Juvenile','Yearling','Adult')), estimate=c(0.463, 0.510, 0.559), lcl=c(0.404, 0.445, 0.499), ucl=c(0.524, 0.574, 0.618))

# Plot Survival by stage
ggplot(survival, aes(stage, estimate, ymin=lcl, ymax=ucl)) + 
  geom_errorbar(width=0.2) + geom_point() + ylim(0,1)
#Juveniles have the lowest survival which makes sense as they are still weak

nestdata <- read.table("MPM practical/gjeroynest.txt", header = TRUE, sep = '\t')
head(nestdata)
#To calculate per capita reproduction, the mean number of chicks prior to fledging, since not all will likely fledge the nest successfully.
ClutchNo <- mean(nestdata$clutchno) 
HatchingSuc <- mean(nestdata$hatchingsuc)
FledglingNo <- mean(nestdata$chickno)

#take averge of clutches multiply by probbaility of clutch hatching, and multiply by the expected number of fletching in each clutch
#We only model the female segment so we divide by 2 to get the number of fletching which are female. We assume equal sex ratio.
PerCapitaRepro <-  (ClutchNo * HatchingSuc * FledglingNo) / 2
#NEed the matrix
Phi.juv <- survival$estimate[survival$stage=='Juvenile'] 
Phi.yr <- survival$estimate[survival$stage=='Yearling'] 
Phi.ad <- survival$estimate[survival$stage=='Adult'] 

sparrowMPM <- c(Phi.juv * PerCapitaRepro, Phi.yr * PerCapitaRepro, Phi.ad * PerCapitaRepro, Phi.juv, 0, 0, 0, Phi.yr, Phi.ad)
sparrowMPM<- matrix(sparrowMPM, nrow=3, ncol=3, byrow=T)
#Calculate lamda using the lambda factor in popbio package
lambda(sparrowMPM) #The sparrow population is slightly declining or nearly stable

# project over 15 years
t <- 15
# start with 50 juveniles, 20 yearlings and 30 adults
n0 <- c(50,20,30)

# project dynamics 
projection <- pop.projection(sparrowMPM, n0, iterations = t)
projected <- data.frame(time=1:15, N=projection$pop.sizes)

# plot projected pop size over time
ggplot(projected, aes(time, N)) + 
  geom_line() + ylim(0,150) + ylab('Projected N')

popest <- read.table("MPM practical/popest.txt", header = TRUE, sep = '\t')
head(popest)
# plot N over time
ggplot(popest, aes(year, N)) + 
  geom_line() + ylim(0,200) + ylab('Observed N')
#The prediction does not match the population at all, The initial population of the model is too high even when compared observed population , and the observed population contains both males and females, so we would need to half it.
#That would not change lamda, so the increase in population is probably due to migrants from other islands

#Look at satble stage distribution which is teh long term average relative abundance of different stages
stages  <- c("Juv","Yr","Ad")
colnames(sparrowMPM) <- stages
rownames(sparrowMPM) <- stages
stable.stage(sparrowMPM)
#Look at reproductive values of different satges, the expected distribution of each individual in the stage to future reproduction
reproductive.value(sparrowMPM)

#Sensitivity and Elasticity
# list the vital rates
sparrow.param <- list(Phi.juv = Phi.juv, Phi.yr = Phi.yr, Phi.ad = Phi.ad, R = PerCapitaRepro)
# give the matrix equation 
sparrow.equation <- expression(Phi.juv * R, Phi.yr * R, Phi.ad * R, Phi.juv, 0, 0, 0, Phi.yr, Phi.ad)
# run the sensitivity analysis
sens <- vitalsens(sparrow.equation, sparrow.param)
sens

# plot elasticity of the vital rates 
sens$vitalrate <- factor(c('Phi.juv', 'Phi.yr', 'Phi.ad', 'R'), levels = c('Phi.juv', 'Phi.yr', 'Phi.ad', 'R'))
ggplot(sens, aes(vitalrate, elasticity)) + 
  geom_bar(stat = 'identity') 
#The most important vital rate for growth are reproductive rate and juvenile survival
#This is not similar to the orca example, where the most important aspect was adult survival, this is because unlike orcas the birds reproduce every before they are one year old they are fast-living species
