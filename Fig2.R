set.seed(1234567890)

library(ggplot2)
library(ggridges)
library(patchwork)
library(plyr)
library(deSolve)
library(MASS)
library(fitdistrplus)
library(xlsx)
library(FME)
library(gtools)
library(RColorBrewer)
library(pheatmap)


setwd("~/Desktop/")


Heatmap0 <- matrix(NA, nrow = 8, ncol = 8)
Heatmap1 <- matrix(NA, nrow = 8, ncol = 8)
Heatmap2 <- matrix(NA, nrow = 8, ncol = 8)
Heatmap3 <- matrix(NA, nrow = 8, ncol = 8)
Heatmap4 <- matrix(NA, nrow = 8, ncol = 8)
Heatmap5 <- matrix(NA, nrow = 8, ncol = 8)

Heatmap00 <- matrix(NA, nrow = 8, ncol = 8)
Heatmap11 <- matrix(NA, nrow = 8, ncol = 8)
Heatmap22 <- matrix(NA, nrow = 8, ncol = 8)
Heatmap33 <- matrix(NA, nrow = 8, ncol = 8)
Heatmap44 <- matrix(NA, nrow = 8, ncol = 8)
Heatmap55 <- matrix(NA, nrow = 8, ncol = 8)

Heatmap000 <- matrix(NA, nrow = 8, ncol = 8)
Heatmap111 <- matrix(NA, nrow = 8, ncol = 8)
Heatmap222 <- matrix(NA, nrow = 8, ncol = 8)
Heatmap333 <- matrix(NA, nrow = 8, ncol = 8)
Heatmap444 <- matrix(NA, nrow = 8, ncol = 8)
Heatmap555 <- matrix(NA, nrow = 8, ncol = 8)


for( iii in 1:8 ) {
  
for( jjj in 1:8 ) {
    
    
##########################################################################################
######################################## Epidemic ########################################
##########################################################################################


#################### Time & Parameter setting
Tmin <- 0
Tmax <- 400
step_size <- 1
N <- 10^7
etime <- seq(Tmin,Tmax,step_size)

R0 <- 1.5 + (iii-1)*0.5
ee1 <- 1/2.7    
ee2 <- 1/3.2
ss1 <- 1/6.6
ss2 <- 1/5.3
pp <- (20 + (jjj-1)*10)/100
bb  <- R0/((1-pp)/ss1 + (pp)/ss2) ## Transmission from R0 = 2 in Nature paper


#################### Epidemic model
Covfun<-function(pars){
  
  b  <- as.numeric(pars[1])
  e1 <- as.numeric(pars[2])
  e2 <- as.numeric(pars[3])
  s1 <- as.numeric(pars[4])
  s2 <- as.numeric(pars[5])
  p  <- as.numeric(pars[6])
  
  derivs<-function(time,y,pars){
    
    with(as.list(c(pars,y)),{
      dS  <--b*S*(I1+I2)
      dE1 <-(1-p)*b*S*(I1+I2) - e1*E1
      dE2 <-p*b*S*(I1+I2) - e2*E2
      dI1 <-e1*E1 - s1*I1  ## Symptomatic
      dI2 <-e2*E2 - s2*I2  ## Asymptomatic
      dR  <-s1*I1 + s2*I2
      dC  <-b*S*(I1+I2)
      dC1 <-(1-p)*b*S*(I1+I2)       ## Symptomatic incidence
      dC2 <-p*b*S*(I1+I2)   ## Asymptomatic incidence
      
      return(list(c(dS,dE1,dE2,dI1,dI2,dR,dC,dC1,dC2)))
    })
  }
  
  y<-c(S=1-1/N,E1=0,E2=0,I1=1/N,I2=0,R=0,C=1/N,C1=1/N,C2=0)
  
  times<-etime
  out<-lsoda(y=y,parms=pars,times=times,func=derivs,rtol=0.00004,atol=0.00000000000001)
  out2<-cbind(time=out[,1],S=out[,2],E1=out[,3],E2=out[,4],I1=out[,5],I2=out[,6],R=out[,7],C=out[,8],C1=out[,9],C2=out[,10],
              II=(out[,5]+out[,6]), # Infectious
              LI=(out[,3]+out[,4]+out[,5]+out[,6]) # Latent + Infectious
  )
  as.data.frame(out2)
}

pars <- c(b = bb, e1 = ee1, e2 = ee2, s1 = ss1, s2 = ss2, p = pp)
out <- Covfun(pars)
R1 <- bb/ss1
R2 <- bb/ss2


#################### Incidence_Total
CC <- c()
for( i in 1:(Tmax/step_size) ) {
  CC[i] <- out$C[i+1] - out$C[i]
}
CC <- data.frame(CC)
CC$time <- seq(Tmin,Tmax-step_size,step_size) 


#################### Incidence_Symptomatic
CC1 <- c()
for( i in 1:(Tmax/step_size) ) {
  CC1[i] <- out$C1[i+1] - out$C1[i]
}
CC1 <- data.frame(CC1)
CC1$time <- seq(Tmin,Tmax-step_size,step_size) 


#################### Incidence_Asymptomatic
CC2 <- c()
for( i in 1:(Tmax/step_size) ) {
  CC2[i] <- out$C2[i+1] - out$C2[i]
}
CC2 <- data.frame(CC2)
CC2$time <- seq(0,Tmax-step_size,step_size) 


Incidence <- cbind(round(1000*CC$CC),round(1000*CC2$CC2),round(1000*CC$CC)-round(1000*CC2$CC2))
Incidence <- data.frame(Incidence)
Incidence$time <- seq(Tmin,Tmax-step_size,step_size)
colnames(Incidence) <- c("total","asym","sym","time")


#################### Effective reproduction number
A <- (out$S)
A <- data.frame(A)
Rt <- A*R0
Rt1 <- A*R1
Rt2 <- A*R2
RR <- cbind(Rt,Rt1,Rt2)
RR$time <- etime
colnames(RR) <- c("R0","R1","R2","time")

Peaktime    <- CC$time[which.max(CC$CC)]
CCC <- subset(CC,CC >= max(CC)*0.1)
Growthtime  <- min(CCC$time)
Declinetime <- max(CCC$time)
Unity <- RR$time[min(which(RR$R0<1))]




######################################################################################################
######################################## Viral dynamics model ########################################
######################################################################################################


#################### Setting
Tmin2 <- 0
Tmax2 <- 150
vtime <- seq(Tmin2,Tmax2,step_size)

IT <- 5                          ## Infectiousness threshold
DL <- log10(2*(10^6))            ## Detection limit for low-sensitivity
sample  <- sum(Incidence$total)  ## Total incidence
ssample <- sum(Incidence$sym)    ## Total sym
asample <- sum(Incidence$asym)   ## Total asym

Meanlog <- 1.7624761   ## Incubation period from Epidemics paper
Sdlog   <- 0.4147944


######################################### Viral dynamics model
Covfun<-function(pars){
  
  beta <- as.numeric(pars[1])
  r <- as.numeric(pars[2])
  delta <- as.numeric(pars[3])
  
  derivs<-function(time,y,pars){
    with(as.list(c(pars,y)),{
      dTa<--beta*Ta*V
      dV<-r*Ta*V-delta*V
      
      return(list(c(dTa,dV)))
    })
  }
  y<-c(Ta=1,V=0.01)
  
  times<-c(seq(Tmin,Tmin+Tmax2,step_size))
  out<-lsoda(y=y,parms=pars,times=times,func=derivs,rtol=0.00004,atol=0)
  out2<-cbind(time=out[,1],aV=((log10(out[,3]))))
  as.data.frame(out2)
}


######################################## Symptomatic
#################### Counting incidence
pop1<- read.csv("populationParameters_sym.txt", row.names = 1) ## Symptomatic

Incidence_sym  <- subset(Incidence, Incidence$sym > 0)
Incidence_sym2 <- rep(Incidence_sym$time[1], times = Incidence_sym$sym[1])

for (i in 2:length(Incidence_sym$sym)){
  
  Incidence_sym3  <- rep(Incidence_sym$time[i], times = Incidence_sym$sym[i])
  Incidence_sym2  <- c(Incidence_sym2,Incidence_sym3)
  
}

Incidence_sym2 <- data.frame(Incidence_sym2)
colnames(Incidence_sym2) <- "sym"


#################### Set of randomly sampled parameters for symptomatic 
fitt1<-matrix(NA,nrow=ssample,ncol=3)
j <- 1
while (TRUE) {
  popr1     <- rlnorm(1, log(pop1["r_pop", "value"]), pop1["omega_r", "value"])  ### Log-normal dist
  popbeta1  <- rlnorm(1, log(pop1["beta_pop", "value"]), pop1["omega_beta", "value"])
  popdelta1 <- rlnorm(1, log(pop1["delta_pop", "value"]), pop1["omega_delta", "value"])
  pars <- c(popbeta1,popr1,popdelta1)
  out <- Covfun(pars)
  
  if (max(out$aV) >= IT) {
    fitt1[j,] <- c(popr1,popdelta1,popbeta1)
    j<-j+1
    if (j>ssample)
      break
  }
}
fitt1 <- data.frame(fitt1)
colnames(fitt1) <- c("r","delta","beta")


######################################## Asymptomatic
#################### Counting incidence
pop2<- read.csv("populationParameters_asym.txt", row.names = 1) ## Symptomatic

Incidence_asym <- subset(Incidence, Incidence$asym > 0)
Incidence_asym2 <- rep(Incidence_asym$time[1], times = Incidence_asym$asym[1])

for (i in 2:length(Incidence_asym$asym)){
  
  Incidence_asym3 <- rep(Incidence_asym$time[i], times = Incidence_asym$asym[i])
  Incidence_asym2  <- c(Incidence_asym2,Incidence_asym3)
  
}

Incidence_asym2 <- data.frame(Incidence_asym2)
colnames(Incidence_asym2) <- "asym"


#################### Set of randomly sampled parameters for asymptimatic 
fitt2<-matrix(NA,nrow=asample,ncol=3)
j <- 1
while (TRUE) {
  popr2     <- rlnorm(1, log(pop2["r_pop", "value"]), pop2["omega_r", "value"])  ### Log-normal dist
  popbeta2  <- rlnorm(1, log(pop2["beta_pop", "value"]), pop2["omega_beta", "value"])
  popdelta2 <- rlnorm(1, log(pop2["delta_pop", "value"]), pop2["omega_delta", "value"])
  pars <- c(popbeta2,popr2,popdelta2)
  out <- Covfun(pars)
  
  if (max(out$aV) >= IT) {
    fitt2[j,] <- c(popr2,popdelta2,popbeta2)
    j<-j+1
    if (j>asample)
      break
  }
}
fitt2 <- data.frame(fitt2)
colnames(fitt2) <- c("r","delta","beta")


############################################################################################
######################################## Simulation ########################################
############################################################################################


  

#################### Imaginary symptomatic individuals
df_total1 = data.frame()
Incidence_sym2$sym <- sample(Incidence_sym2$sym,ssample,replace=FALSE)

for( i in 1:ssample ) {
  Tmin  <- Incidence_sym2$sym[i]
  stime <- seq(Tmin,Tmin+Tmax2,step_size)
  pars  <- c(beta=fitt1$beta[i],r=fitt1$r[i],delta=fitt1$delta[i])
  out   <- Covfun(pars)
  
  ppp<-matrix(NA,nrow=length(stime),ncol=1)
  for(jj in 1:length(stime)){
    ii<-jj
    a<-stime[ii]
    dd<-out[out$time==a,]
    ppp[ii]<-dd$aV
  }
  pp<-data.frame(ID=i, time=stime, V=ppp)
  df_total1<-rbind(df_total1,pp)
  
}


#################### Imaginary asymptomatic individuals
df_total2 = data.frame()
Incidence_asym2$asym <- sample(Incidence_asym2$asym,asample,replace=FALSE)

for( i in 1:asample ) {
  Tmin  <- Incidence_asym2$asym[i]
  stime <- seq(Tmin,Tmin+Tmax2,step_size)
  pars  <- c(beta=fitt2$beta[i],r=fitt2$r[i],delta=fitt2$delta[i])
  out   <- Covfun(pars)
  
  ppp<-matrix(NA,nrow=length(stime),ncol=1)
  for(jj in 1:length(stime)){
    ii<-jj
    a<-stime[ii]
    dd<-out[out$time==a,]
    ppp[ii]<-dd$aV
  }
  pp<-data.frame(ID=i, time=stime, V=ppp)
  df_total2<-rbind(df_total2,pp)
  
}

df_total2$g <- rep(Tmax2*100,times=(Tmax2+1)*asample)
df_total2$Time <- rep(vtime, times = asample)




#################### Start Screening
iteration <- 100
screening <- 10+1
testings  <- 4
Timing <- Peaktime


Scenario0 <- matrix(NA, nrow = iteration, ncol = screening)
Scenario1 <- matrix(NA, nrow = iteration, ncol = screening)
Scenario2 <- matrix(NA, nrow = iteration, ncol = screening)
Scenario3 <- matrix(NA, nrow = iteration, ncol = screening)
Scenario4 <- matrix(NA, nrow = iteration, ncol = screening)
Scenario5 <- matrix(NA, nrow = iteration, ncol = screening)


for (k in 1:iteration) {


#################### Incubation period
df_total111    <- round(rlnorm(ssample, Meanlog, Sdlog), digits = 0)
df_total1111   <- df_total111 + Incidence_sym2
df_total11111  <- rep(df_total1111$sym, each = length(vtime))
df_total1$g    <- df_total11111
df_total1$Time <- rep(vtime, times = ssample)


#################### Errors
error1 <- rnorm(ssample*(Tmax2+1), mean=0, sd=pop1$value[9]) ## Constant error 
df_total11 <- df_total1$V+error1
df_total1$obs <- df_total11

error2 <- rnorm(asample*(Tmax2+1), mean=0, sd=pop2$value[9]) ## Constant error 
df_total22 <- df_total2$V+error2
df_total2$obs <- df_total22


Total <- rbind(df_total1,df_total2)
Total$ID <- rep(1:sample, each = length(vtime))
TM <- length(Total$ID)


#################### Denominator
Number <- subset(Total,Total$Time==0)
Number2 <- subset(Number,Number$time <= Timing & Number$g >= Timing)
K <- length(Number2$ID)


#################### COVID-19 patients at t = t_0
Total1  <- subset(Total,Total$ID==Number2$ID[1])

for (i in 2:length(Number2$ID)){
  
  Total1 <- rbind(Total1,subset(Total,Total$ID==Number2$ID[i]))
  
}

Total2 <- subset(Total1,Total1$time >= Timing & Total1$time <= Timing+10)




###########################################################################################
######################################## Scenario0 ########################################
###########################################################################################

for (s in 1:screening){
  
  Total3 <- subset(Total2,Total2$time == Timing + (s-1))
  
  Scenario0[k,s] <- length(Total3$ID[which(Total3$g == Timing + (s-1))])
  
  
}

Scenario0[k,] <- cumsum(Scenario0[k,])/K*100




###########################################################################################
######################################## Scenario1 ########################################
###########################################################################################

SC1 <- c(1,1,1,1,0,0,0,0,0,0,0)

for (s in 1:screening){

    Total3 <- subset(Total2,Total2$time == Timing + (s-1))

    if (SC1[s]==1){

      TT <- assign(paste0('T',s-1), Total3$ID[which(Total3$obs >= DL & Total3$g > Timing+(s-1))])
      SS <- assign(paste0('S',s-1), Total3$ID[which(Total3$g == Timing+(s-1))])

      }else {

      TT <- assign(paste0('T',s-1), Total3$ID[which(Total3$g == Timing+(s-1))])
      SS <- assign(paste0('S',s-1), Total3$ID[which(Total3$g == Timing+(s-1))])

      }

    assign(paste0('U',s-1), union(TT,SS))

}

U0  <- union(U0,U0)
U1  <- union(U0,U1)
U2  <- union(U1,U2)
U3  <- union(U2,U3)
U4  <- union(U3,U4)
U5  <- union(U4,U5)
U6  <- union(U5,U6)
U7  <- union(U6,U7)
U8  <- union(U7,U8)
U9  <- union(U8,U9)
U10 <- union(U9,U10)

Scenario1[k,] <- c(length(U0),
                   length(U1),
                   length(U2),
                   length(U3),
                   length(U4),
                   length(U5),
                   length(U6),
                   length(U7),
                   length(U8),
                   length(U9),
                   length(U10))/K*100

###########################################################################################
######################################## Scenario2 ########################################
###########################################################################################

SC2 <- c(1,0,1,0,1,0,1,0,0,0,0)

for (s in 1:screening){
  
  Total3 <- subset(Total2,Total2$time == Timing + (s-1))
  
  if (SC2[s]==1){
    
    TT <- assign(paste0('T',s-1), Total3$ID[which(Total3$obs >= DL & Total3$g > Timing+(s-1))])
    SS <- assign(paste0('S',s-1), Total3$ID[which(Total3$g == Timing+(s-1))])
    
  }else {
    
    TT <- assign(paste0('T',s-1), Total3$ID[which(Total3$g == Timing+(s-1))])
    SS <- assign(paste0('S',s-1), Total3$ID[which(Total3$g == Timing+(s-1))])
    
  }
  
  assign(paste0('U',s-1), union(TT,SS))
  
}

U0  <- union(U0,U0)
U1  <- union(U0,U1)
U2  <- union(U1,U2)
U3  <- union(U2,U3)
U4  <- union(U3,U4)
U5  <- union(U4,U5)
U6  <- union(U5,U6)
U7  <- union(U6,U7)
U8  <- union(U7,U8)
U9  <- union(U8,U9)
U10 <- union(U9,U10)

Scenario2[k,] <- c(length(U0),
                   length(U1),
                   length(U2),
                   length(U3),
                   length(U4),
                   length(U5),
                   length(U6),
                   length(U7),
                   length(U8),
                   length(U9),
                   length(U10))/K*100




###########################################################################################
######################################## Scenario3 ########################################
###########################################################################################

SC3 <- c(1,0,0,1,0,0,1,0,0,1,0)

for (s in 1:screening){
  
  Total3 <- subset(Total2,Total2$time == Timing + (s-1))
  
  if (SC3[s]==1){
    
    TT <- assign(paste0('T',s-1), Total3$ID[which(Total3$obs >= DL & Total3$g > Timing+(s-1))])
    SS <- assign(paste0('S',s-1), Total3$ID[which(Total3$g == Timing+(s-1))])
    
  }else {
    
    TT <- assign(paste0('T',s-1), Total3$ID[which(Total3$g == Timing+(s-1))])
    SS <- assign(paste0('S',s-1), Total3$ID[which(Total3$g == Timing+(s-1))])
    
  }
  
  assign(paste0('U',s-1), union(TT,SS))
  
}

U0  <- union(U0,U0)
U1  <- union(U0,U1)
U2  <- union(U1,U2)
U3  <- union(U2,U3)
U4  <- union(U3,U4)
U5  <- union(U4,U5)
U6  <- union(U5,U6)
U7  <- union(U6,U7)
U8  <- union(U7,U8)
U9  <- union(U8,U9)
U10 <- union(U9,U10)

Scenario3[k,] <- c(length(U0),
                   length(U1),
                   length(U2),
                   length(U3),
                   length(U4),
                   length(U5),
                   length(U6),
                   length(U7),
                   length(U8),
                   length(U9),
                   length(U10))/K*100




###########################################################################################
######################################## Scenario4 ########################################
###########################################################################################

SC4 <- c(1,0,0,0,0,0,0,0,1,1,1)

for (s in 1:screening){
  
  Total3 <- subset(Total2,Total2$time == Timing + (s-1))
  
  if (SC4[s]==1){
    
    TT <- assign(paste0('T',s-1), Total3$ID[which(Total3$obs >= DL & Total3$g > Timing+(s-1))])
    SS <- assign(paste0('S',s-1), Total3$ID[which(Total3$g == Timing+(s-1))])
    
  }else {
    
    TT <- assign(paste0('T',s-1), Total3$ID[which(Total3$g == Timing+(s-1))])
    SS <- assign(paste0('S',s-1), Total3$ID[which(Total3$g == Timing+(s-1))])
    
  }
  
  assign(paste0('U',s-1), union(TT,SS))
  
}

U0  <- union(U0,U0)
U1  <- union(U0,U1)
U2  <- union(U1,U2)
U3  <- union(U2,U3)
U4  <- union(U3,U4)
U5  <- union(U4,U5)
U6  <- union(U5,U6)
U7  <- union(U6,U7)
U8  <- union(U7,U8)
U9  <- union(U8,U9)
U10 <- union(U9,U10)

Scenario4[k,] <- c(length(U0),
                   length(U1),
                   length(U2),
                   length(U3),
                   length(U4),
                   length(U5),
                   length(U6),
                   length(U7),
                   length(U8),
                   length(U9),
                   length(U10))/K*100




###########################################################################################
######################################## Scenario5 ########################################
###########################################################################################

SC5 <- c(1,0,0,0,0,0,0,0,0,0,0)
random <- sample(2:screening, 3, replace=F)
SC5[random] <- 1 

for (s in 1:screening){
  
  Total3 <- subset(Total2,Total2$time == Timing + (s-1))
  
  if (SC5[s]==1){
    
    TT <- assign(paste0('T',s-1), Total3$ID[which(Total3$obs >= DL & Total3$g > Timing+(s-1))])
    SS <- assign(paste0('S',s-1), Total3$ID[which(Total3$g == Timing+(s-1))])
    
  }else {
    
    TT <- assign(paste0('T',s-1), Total3$ID[which(Total3$g == Timing+(s-1))])
    SS <- assign(paste0('S',s-1), Total3$ID[which(Total3$g == Timing+(s-1))])
    
  }
  
  assign(paste0('U',s-1), union(TT,SS))
  
}

U0  <- union(U0,U0)
U1  <- union(U0,U1)
U2  <- union(U1,U2)
U3  <- union(U2,U3)
U4  <- union(U3,U4)
U5  <- union(U4,U5)
U6  <- union(U5,U6)
U7  <- union(U6,U7)
U8  <- union(U7,U8)
U9  <- union(U8,U9)
U10 <- union(U9,U10)

Scenario5[k,] <- c(length(U0),
                   length(U1),
                   length(U2),
                   length(U3),
                   length(U4),
                   length(U5),
                   length(U6),
                   length(U7),
                   length(U8),
                   length(U9),
                   length(U10))/K*100


}


# Heatmap0[iii,jjj] <- Scenario0[1,screening]
# Heatmap1[iii,jjj] <- Scenario1[1,screening]
# Heatmap2[iii,jjj] <- Scenario2[1,screening]
# Heatmap3[iii,jjj] <- Scenario3[1,screening]
# Heatmap4[iii,jjj] <- Scenario4[1,screening]
# Heatmap5[iii,jjj] <- Scenario5[1,screening]


Heatmap0[iii,jjj] <- mean(Scenario0[,screening])
Heatmap1[iii,jjj] <- mean(Scenario1[,screening])
Heatmap2[iii,jjj] <- mean(Scenario2[,screening])
Heatmap3[iii,jjj] <- mean(Scenario3[,screening])
Heatmap4[iii,jjj] <- mean(Scenario4[,screening])
Heatmap5[iii,jjj] <- mean(Scenario5[,screening])

Heatmap00[iii,jjj] <- quantile(Scenario0[,screening],0.025)
Heatmap11[iii,jjj] <- quantile(Scenario1[,screening],0.025)
Heatmap22[iii,jjj] <- quantile(Scenario2[,screening],0.025)
Heatmap33[iii,jjj] <- quantile(Scenario3[,screening],0.025)
Heatmap44[iii,jjj] <- quantile(Scenario4[,screening],0.025)
Heatmap55[iii,jjj] <- quantile(Scenario5[,screening],0.025)

Heatmap000[iii,jjj] <- quantile(Scenario0[,screening],0.975)
Heatmap111[iii,jjj] <- quantile(Scenario1[,screening],0.975)
Heatmap222[iii,jjj] <- quantile(Scenario2[,screening],0.975)
Heatmap333[iii,jjj] <- quantile(Scenario3[,screening],0.975)
Heatmap444[iii,jjj] <- quantile(Scenario4[,screening],0.975)
Heatmap555[iii,jjj] <- quantile(Scenario5[,screening],0.975)


} #### Proportion
  
} #### R0


Heatmap_0 <- data.frame(R0=rep(seq(1.5,5.0,by=0.5),times=8),Proportion=rep(seq(20,90,by=10),each=8),Value=as.vector(Heatmap0),Min=as.vector(Heatmap00),Max=as.vector(Heatmap000))
Heatmap_1 <- data.frame(R0=rep(seq(1.5,5.0,by=0.5),times=8),Proportion=rep(seq(20,90,by=10),each=8),Value=as.vector(Heatmap1),Min=as.vector(Heatmap11),Max=as.vector(Heatmap111))
Heatmap_2 <- data.frame(R0=rep(seq(1.5,5.0,by=0.5),times=8),Proportion=rep(seq(20,90,by=10),each=8),Value=as.vector(Heatmap2),Min=as.vector(Heatmap22),Max=as.vector(Heatmap222))
Heatmap_3 <- data.frame(R0=rep(seq(1.5,5.0,by=0.5),times=8),Proportion=rep(seq(20,90,by=10),each=8),Value=as.vector(Heatmap3),Min=as.vector(Heatmap33),Max=as.vector(Heatmap333))
Heatmap_4 <- data.frame(R0=rep(seq(1.5,5.0,by=0.5),times=8),Proportion=rep(seq(20,90,by=10),each=8),Value=as.vector(Heatmap4),Min=as.vector(Heatmap44),Max=as.vector(Heatmap444))
Heatmap_5 <- data.frame(R0=rep(seq(1.5,5.0,by=0.5),times=8),Proportion=rep(seq(20,90,by=10),each=8),Value=as.vector(Heatmap5),Min=as.vector(Heatmap55),Max=as.vector(Heatmap555))


write.xlsx(Heatmap_0,"Heatmap0.xlsx")
write.xlsx(Heatmap_1,"Heatmap1.xlsx")
write.xlsx(Heatmap_2,"Heatmap2.xlsx")
write.xlsx(Heatmap_3,"Heatmap3.xlsx")
write.xlsx(Heatmap_4,"Heatmap4.xlsx")
write.xlsx(Heatmap_5,"Heatmap5.xlsx")







