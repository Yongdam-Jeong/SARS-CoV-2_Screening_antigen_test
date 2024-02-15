library(ggplot2)
library(ggridges)
library(patchwork)
library(plyr)
library(deSolve)
library(MASS)
library(fitdistrplus)
library(xlsx)
library(writexl)
library(FME)
library(gtools)
library(RColorBrewer)
library(pheatmap)
library(boot)


setwd("~/Desktop")


############################################################################################
##################################### Epidemic setting #####################################
############################################################################################

##################################### Time
Tmin <- 0
Tmax <- 400
step_size <- 1
N <- 10^7 ## Community size
M <- 10^3 ## Facility size
mtime <- seq(Tmin,Tmax,step_size)


##################################### Parameters
R0 <- 2 ## Reproduction number
e1 <- 1/3.209127601 ## Mean latent period for asym
e2 <- 1/2.745534106 ## Mean latent period for sym
s1 <- 1/6.032712583 ## Mean infectious period for asym
s2 <- 1/7.452044219 ## Mean infectious period for sym
p <- 0.7 ## Proportion of asym
b <- R0/(p/s1 + (1-p)/s2) ## Transmission rate derived from R0


##################################### Model
Epidemic <- function(pars){
  
  b  <- as.numeric(pars[1])
  e1 <- as.numeric(pars[2])
  e2 <- as.numeric(pars[3])
  s1 <- as.numeric(pars[4])
  s2 <- as.numeric(pars[5])
  p  <- as.numeric(pars[6])
  
  derivs<-function(time,y,pars){
    
    with(as.list(c(pars,y)),{
      
      ## Susceptibles
      dS  <- -b*S*(I11+I12+I21+I22)
      
      ## Asymptomatic
      dL11 <- p*b*S*(I11+I12+I21+I22) - 2*e1*L11
      dL12 <- 2*e1*L11 - 2*e1*L12
      dI11 <- 2*e1*L12 - 2*s1*I11
      dI12 <- 2*s1*I11 - 2*s1*I12
      
      ## Symptomatic
      dL21 <- (1-p)*b*S*(I11+I12+I21+I22) - 2*e2*L21
      dL22 <- 2*e2*L21 - 2*e2*L22
      dI21 <- 2*e2*L22 - 2*s2*I21
      dI22 <- 2*s2*I21 - 2*s2*I22
      
      ## Removed 
      dR  <- 2*s1*I12 + 2*s2*I22
      
      dC  <- b*S*(I11+I12+I21+I22)       ## Total
      dC1 <- p*b*S*(I11+I12+I21+I22)     ## Asymptomatic 
      dC2 <- (1-p)*b*S*(I11+I12+I21+I22) ## Symptomatic
      
      
      return(list(c(dS,
                    dL11,dL12,dL21,dL22,
                    dI11,dI12,dI21,dI22,
                    dR,
                    dC,dC1,dC2)))
    })
  }
  
  y <- c(S=1-1/N,
         L11=0,L12=0,L21=0,L22=0,
         I11=1/N,I12=0,I21=0,I22=0,
         R=0,
         C=1/N,C1=1/N,C2=0)
  
  times <- mtime
  out <- lsoda(y=y,parms=pars,times=times,func=derivs,rtol=0.00004,atol=0.00000000000001)
  out2 <- cbind(time=out[,1],S=out[,2],R=out[,11],C=out[,12],C1=out[,13],C2=out[,14])
  as.data.frame(out2)
}


##################################### Incidence over epidemic time #####################################
pars <- c(b,e1,e2,s1,s2,p)
out  <- Epidemic(pars)


#################### Total
IC <- c()
for( i in 1:(Tmax/step_size) ) {
  IC[i] <- out$C[i+1] - out$C[i]
}
IC <- data.frame(IC)
IC$time <- seq(Tmin,Tmax-step_size,step_size) 
colnames(IC) <- c("Incidence","time")


#################### Asymptomatic
IC1 <- c()
for( i in 1:(Tmax/step_size) ) {
  IC1[i] <- out$C1[i+1] - out$C1[i]
}
IC1 <- data.frame(IC1)
IC1$time <- seq(Tmin,Tmax-step_size,step_size) 
colnames(IC1) <- c("Incidence","time")


#################### Symptomatic
IC2 <- c()
for( i in 1:(Tmax/step_size) ) {
  IC2[i] <- out$C2[i+1] - out$C2[i]
}
IC2 <- data.frame(IC2)
IC2$time <- seq(Tmin,Tmax-step_size,step_size) 
colnames(IC2) <- c("Incidence","time")


Incidence <- round(M*cbind(IC1$Incidence,IC2$Incidence))
Incidence <- data.frame(Incidence)
Incidence$time <- seq(Tmin,Tmax-step_size,step_size)
colnames(Incidence) <- c("asym","sym","time")


Peak   <- out$time[which.max(IC$Incidence)]
ICpeak <- IC$Incidence[which(IC$time==Peak)]
IC10   <- subset(IC,IC$Incidence >= ICpeak*0.1)

Growth <- min(IC10$time)
Decline <- max(IC10$time)


##########################################################################################
##################################### Viral dynamics #####################################
##########################################################################################


#################### Setting
Tmin2 <- 0
Tmax2 <- 150
vtime <- seq(Tmin2,Tmax2,step_size)

IT  <- log10(10^5) ## Infectiousness thresholds 
PCR <- log10(10^2) ## Detection limit of PCR testing
NON <- log10(10^0) ## Threshold
DL  <- log10(2*(10^6)) ## Detection limit for low-sensitivity
asample <- sum(Incidence$asym) ## Total asym
ssample <- sum(Incidence$sym) ## Total sym
sample  <- asample + ssample  ## Total incidence


Meanlog <- 1.7624761   ## Distribution of incubation period 
Sdlog   <- 0.4147944


##################################### Model
Covfun <- function(pars){
  
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
  
  times<-c(seq(Tmin2,Tmax2,step_size))
  out<-lsoda(y=y,parms=pars,times=times,func=derivs,rtol=0.00004,atol=0)
  out2<-cbind(time=out[,1],V=((log10(out[,3]))))
  as.data.frame(out2)
}


##################################### Asymptomatic
pop1 <- read.csv("populationParameters_asym.txt", row.names = 1)

Incidence2 <- subset(Incidence, Incidence$asym > 0)
Incidence3 <- rep(Incidence2$time[1], times = Incidence2$asym[1])

for (i in 2:length(Incidence2$asym)){
  
  Incidence4  <- rep(Incidence2$time[i], times = Incidence2$asym[i])
  Incidence3  <- c(Incidence3,Incidence4)
  
}

Incidence3 <- data.frame(Incidence3)
colnames(Incidence3) <- "asym"


fitt <- matrix(NA,nrow=asample,ncol=3)
j <- 1
while (TRUE) {
  
  b <- rlnorm(1, log(pop1["beta_pop", "value"]), pop1["omega_beta", "value"])
  r <- rlnorm(1, log(pop1["r_pop", "value"]), pop1["omega_r", "value"])
  d <- rlnorm(1, log(pop1["delta_pop", "value"]), pop1["omega_delta", "value"])
  
  pars <- c(b,r,d)
  out  <- Covfun(pars)
  
  if (max(out$V) > IT & out$V[which(out$time==Tmax2)] < NON) {
    fitt[j,] <- c(b,r,d)
    j<-j+1
    if (j>asample)
      break
  }
}

fitt <- data.frame(fitt)
colnames(fitt) <- c("b","r","d")

df_total1 = data.frame()
Incidence3$asym <- sample(Incidence3$asym,asample,replace=FALSE)

for ( i in 1:asample ) {
  
  Tmin  <- Incidence3$asym[i]
  stime <- seq(Tmin,Tmin+Tmax2,step_size)
  pars  <- c(fitt$b[i],fitt$r[i],fitt$d[i])
  out   <- Covfun(pars)
  out$time <- stime
  
  out2 <- subset(out,V>=IT)
  latent <- rep(min(out2$time)-step_size, times=length(stime))
  peaktime <- rep(out2$time[which.max(out2$V)], times=length(stime))
  infectious <- rep(max(out2$time), times=length(stime))
  
  detected <- rep(max(out$time[out$V>=NON]), times=length(stime))
  
  ppp   <- matrix(NA,nrow=length(stime),ncol=1)
  
  for (jj in 1:length(stime)) {
    
    ii <- jj
    a  <- stime[ii]
    dd <- out[out$time==a,]
    ppp[ii] <- dd$V
    
  }
  
  pp<-data.frame(ID=i, time=stime, V=ppp, L=latent, P=peaktime, I=infectious, D=detected)
  df_total1<-rbind(df_total1,pp)
  
}

df_total1$g <- rep(Tmax2*100,times=(Tmax2+1)*asample) ## no meaning due to no symptoms
df_total1$Time <- rep(vtime, times = asample) ## within-host time scale


##################################### Symptomatic
pop2 <- read.csv("populationParameters_sym.txt", row.names = 1)

Incidence2 <- subset(Incidence, Incidence$sym > 0)
Incidence3 <- rep(Incidence2$time[1], times = Incidence2$sym[1])

for (i in 2:length(Incidence2$sym)){
  
  Incidence4  <- rep(Incidence2$time[i], times = Incidence2$sym[i])
  Incidence3  <- c(Incidence3,Incidence4)
  
}

Incidence3 <- data.frame(Incidence3)
colnames(Incidence3) <- "sym"


fitt <- matrix(NA,nrow=ssample,ncol=3)
j <- 1
while (TRUE) {
  
  b <- rlnorm(1, log(pop2["beta_pop", "value"]), pop2["omega_beta", "value"])
  r <- rlnorm(1, log(pop2["r_pop", "value"]), pop2["omega_r", "value"])
  d <- rlnorm(1, log(pop2["delta_pop", "value"]), pop2["omega_delta", "value"])
  
  pars <- c(b,r,d)
  out  <- Covfun(pars)
  
  if (max(out$V) >= IT & out$V[which(out$time==Tmax2)] < NON) {
    fitt[j,] <- c(b,r,d)
    j<-j+1
    if (j>ssample)
      break
  }
}
fitt <- data.frame(fitt)
colnames(fitt) <- c("b","r","d")

df_total2 = data.frame()
Incidence3$sym <- sample(Incidence3$sym,ssample,replace=FALSE)

for ( i in 1:ssample ) {
  
  Tmin  <- Incidence3$sym[i]
  stime <- seq(Tmin,Tmin+Tmax2,step_size)
  pars  <- c(fitt$b[i],fitt$r[i],fitt$d[i])
  out   <- Covfun(pars)
  out$time <- stime
 
  out2 <- subset(out,V>=IT)
  latent <- rep(min(out2$time)-step_size, times=length(stime))
  peaktime <- rep(out2$time[which.max(out2$V)], times=length(stime))
  infectious <- rep(max(out2$time), times=length(stime))
  
  detected <- rep(max(out$time[out$V>=NON]), times=length(stime))
  
  ppp   <- matrix(NA,nrow=length(stime),ncol=1)
  
  for (jj in 1:length(stime)) {
    
    ii <- jj
    a  <- stime[ii]
    dd <- out[out$time==a,]
    ppp[ii] <- dd$V
    
  }
  
  pp<-data.frame(ID=i, time=stime, V=ppp, L=latent, P=peaktime, I=infectious, D=detected)
  df_total2<-rbind(df_total2,pp)
  
}



############################################################################################
############################################################################################
######################################## Simulation ########################################
############################################################################################
############################################################################################


#################### Start Screening
iteration <- 100
screening <- 10+1  ## Duration of screening
testings  <- 4     ## # of available tests
Timing <- Peak     ## Initiation timing of screening

######### Identified proportion by symptom and testing 
Schedule0 <- matrix(NA, nrow = iteration, ncol = screening)
Schedule1 <- matrix(NA, nrow = iteration, ncol = screening)
Schedule2 <- matrix(NA, nrow = iteration, ncol = screening)
Schedule3 <- matrix(NA, nrow = iteration, ncol = screening)
Schedule4 <- matrix(NA, nrow = iteration, ncol = screening)
Schedule5 <- matrix(NA, nrow = iteration, ncol = screening)

######### Identified proportion by only testing 
ScheduleT1 <- matrix(NA, nrow = iteration, ncol = screening)
ScheduleT2 <- matrix(NA, nrow = iteration, ncol = screening)
ScheduleT3 <- matrix(NA, nrow = iteration, ncol = screening)
ScheduleT4 <- matrix(NA, nrow = iteration, ncol = screening)
ScheduleT5 <- matrix(NA, nrow = iteration, ncol = screening)

######### Identified number by testing for each number of antigen test
ScheduleN1 <- matrix(NA, nrow = iteration, ncol = testings)
ScheduleN2 <- matrix(NA, nrow = iteration, ncol = testings)
ScheduleN3 <- matrix(NA, nrow = iteration, ncol = testings)
ScheduleN4 <- matrix(NA, nrow = iteration, ncol = testings)
ScheduleN5 <- matrix(NA, nrow = iteration, ncol = testings)

######### Identified number by symptoms, but already used number of antigen test
ScheduleS1 <- matrix(NA, nrow = iteration, ncol = testings)
ScheduleS2 <- matrix(NA, nrow = iteration, ncol = testings)
ScheduleS3 <- matrix(NA, nrow = iteration, ncol = testings)
ScheduleS4 <- matrix(NA, nrow = iteration, ncol = testings)
ScheduleS5 <- matrix(NA, nrow = iteration, ncol = testings)

######### Target & Failed to identify cases
ScheduleF1 <- matrix(NA, nrow = iteration, ncol = 2)
ScheduleF2 <- matrix(NA, nrow = iteration, ncol = 2)
ScheduleF3 <- matrix(NA, nrow = iteration, ncol = 2)
ScheduleF4 <- matrix(NA, nrow = iteration, ncol = 2)
ScheduleF5 <- matrix(NA, nrow = iteration, ncol = 2)


for (k in 1:iteration) {


#################### Incubation period for symptomatic
IP  <- round(rlnorm(ssample, Meanlog, Sdlog), digits = 0)
IP2 <- IP + Incidence3
IP3 <- rep(IP2$sym, each = length(vtime))
df_total2$g <- IP3  ## 
df_total2$Time <- rep(vtime, times = ssample)


#################### Errors
error <- rnorm(asample*(Tmax2+1), mean=0, sd=pop1$value[9]) ## Constant error 
error <- df_total1$V+error
df_total1$obs <- error

error <- rnorm(ssample*(Tmax2+1), mean=0, sd=pop2$value[9]) ## Constant error 
error <- df_total2$V+error
df_total2$obs <- error

Total <- rbind(df_total1,df_total2)
Total$ID <- rep(1:sample, each = length(vtime))
TM <- length(Total$ID)


#################### Denominator (i.e., target for screening at the initiation)
Number <- subset(Total,Total$Time==0)  
Number2 <- subset(Number,Number$time <= Timing & Number$g >= Timing & Number$D >= Timing)
K <- length(Number2$ID)


#################### COVID-19 cases at t = t_0 (i.e., K(t_0))
Total1  <- subset(Total,Total$ID==Number2$ID[1])
for (i in 2:length(Number2$ID)){
  
  Total1 <- rbind(Total1,subset(Total,Total$ID==Number2$ID[i]))
  
}

Total1 <- subset(Total1,Total1$D >= Timing)
Total2 <- subset(Total1,Total1$time >= Timing & Total1$time <= Timing+10)  ## Consider only cases during screening programs for reducing computations 




###########################################################################################
######################################## Schedule0 ########################################
###########################################################################################

for (s in 1:screening){
  
  Total3 <- subset(Total2,Total2$time == Timing + (s-1))
  
  Schedule0[k,s] <- length(Total3$ID[which(Total3$g == Timing + (s-1) & Total3$D >= Timing + (s-1))])
  
  
}

Schedule0[k,] <- cumsum(Schedule0[k,])/K*100




###########################################################################################
######################################## Schedule1 ########################################
###########################################################################################

SC1 <- c(1,1,1,1,0,0,0,0,0,0,0)

for (s in 1:screening){

    Total3 <- subset(Total2,Total2$time == Timing + (s-1))

    if (SC1[s]==1){

      TT <- assign(paste0('T',s-1), Total3$ID[which(Total3$obs >= DL & Total3$g > Timing+(s-1) & Total3$D >= Timing + (s-1))])
      SS <- assign(paste0('S',s-1), Total3$ID[which(Total3$g == Timing+(s-1) & Total3$D >= Timing + (s-1))])

      }else {

      SS <- assign(paste0('S',s-1), Total3$ID[which(Total3$g == Timing+(s-1) & Total3$D >= Timing + (s-1))])

      }

    assign(paste0('U',s-1), union(TT,SS))
    assign(paste0('V',s-1), TT)

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

Schedule1[k,] <- c(length(U0),
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

V0  <- union(V0,V0)
V1  <- union(V0,V1)
V2  <- union(V1,V2)
V3  <- union(V2,V3)
V4  <- union(V3,V4)
V5  <- union(V4,V5)
V6  <- union(V5,V6)
V7  <- union(V6,V7)
V8  <- union(V7,V8)
V9  <- union(V8,V9)
V10 <- union(V9,V10)

ScheduleT1[k,] <- c(length(V0),
                    length(V1),
                    length(V2),
                    length(V3),
                    length(V4),
                    length(V5),
                    length(V6),
                    length(V7),
                    length(V8),
                    length(V9),
                    length(V10))/K*100

Tested <- c(length(V0),
            length(V1)-length(V0),
            length(V2)-length(V1),
            length(V3)-length(V2),
            length(V4)-length(V3),
            length(V5)-length(V4),
            length(V6)-length(V5),
            length(V7)-length(V6),
            length(V8)-length(V7),
            length(V9)-length(V8),
            length(V10)-length(V9))

W0  <- setdiff(U0,V0)
W1  <- setdiff(U1,V1)
W2  <- setdiff(U2,V2)
W3  <- setdiff(U3,V3)
W4  <- setdiff(U4,V4)
W5  <- setdiff(U5,V5)
W6  <- setdiff(U6,V6)
W7  <- setdiff(U7,V7)
W8  <- setdiff(U8,V8)
W9  <- setdiff(U9,V9)
W10 <- setdiff(U10,V10)

Symptoms <- c(length(W0),
              length(W1)-length(W0),
              length(W2)-length(W1),
              length(W3)-length(W2),
              length(W4)-length(W3),
              length(W5)-length(W4),
              length(W6)-length(W5),
              length(W7)-length(W6),
              length(W8)-length(W7),
              length(W9)-length(W8),
              length(W10)-length(W9))

ScheduleN1[k,] <- c(Tested[1],Tested[2],Tested[3],Tested[4])
ScheduleS1[k,] <- c(sum(Symptoms[2]),sum(Symptoms[3]),sum(Symptoms[4]),sum(Symptoms[5:11]))
ScheduleF1[k,] <- c(K,K-length(U10))




###########################################################################################
######################################## Schedule2 ########################################
###########################################################################################

SC2 <- c(1,0,1,0,1,0,1,0,0,0,0)

for (s in 1:screening){
  
  Total3 <- subset(Total2,Total2$time == Timing + (s-1))
  
  if (SC2[s]==1){
    
    TT <- assign(paste0('T',s-1), Total3$ID[which(Total3$obs >= DL & Total3$g > Timing+(s-1) & Total3$D >= Timing + (s-1))])
    SS <- assign(paste0('S',s-1), Total3$ID[which(Total3$g == Timing+(s-1) & Total3$D >= Timing + (s-1))])
    
  }else {
    
    SS <- assign(paste0('S',s-1), Total3$ID[which(Total3$g == Timing+(s-1) & Total3$D >= Timing + (s-1))])
    
  }
  
  assign(paste0('U',s-1), union(TT,SS))
  assign(paste0('V',s-1), TT)
  
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

Schedule2[k,] <- c(length(U0),
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

V0  <- union(V0,V0)
V1  <- union(V0,V1)
V2  <- union(V1,V2)
V3  <- union(V2,V3)
V4  <- union(V3,V4)
V5  <- union(V4,V5)
V6  <- union(V5,V6)
V7  <- union(V6,V7)
V8  <- union(V7,V8)
V9  <- union(V8,V9)
V10 <- union(V9,V10)

ScheduleT2[k,] <- c(length(V0),
                    length(V1),
                    length(V2),
                    length(V3),
                    length(V4),
                    length(V5),
                    length(V6),
                    length(V7),
                    length(V8),
                    length(V9),
                    length(V10))/K*100

Tested <- c(length(V0),
            length(V1)-length(V0),
            length(V2)-length(V1),
            length(V3)-length(V2),
            length(V4)-length(V3),
            length(V5)-length(V4),
            length(V6)-length(V5),
            length(V7)-length(V6),
            length(V8)-length(V7),
            length(V9)-length(V8),
            length(V10)-length(V9))

W0  <- setdiff(U0,V0)
W1  <- setdiff(U1,V1)
W2  <- setdiff(U2,V2)
W3  <- setdiff(U3,V3)
W4  <- setdiff(U4,V4)
W5  <- setdiff(U5,V5)
W6  <- setdiff(U6,V6)
W7  <- setdiff(U7,V7)
W8  <- setdiff(U8,V8)
W9  <- setdiff(U9,V9)
W10 <- setdiff(U10,V10)

Symptoms <- c(length(W0),
              length(W1)-length(W0),
              length(W2)-length(W1),
              length(W3)-length(W2),
              length(W4)-length(W3),
              length(W5)-length(W4),
              length(W6)-length(W5),
              length(W7)-length(W6),
              length(W8)-length(W7),
              length(W9)-length(W8),
              length(W10)-length(W9))

ScheduleN2[k,] <- c(Tested[1],Tested[3],Tested[5],Tested[7])
ScheduleS2[k,] <- c(sum(Symptoms[2:3]),sum(Symptoms[4:5]),sum(Symptoms[6:7]),sum(Symptoms[8:11]))
ScheduleF2[k,] <- c(K,K-length(U10))




###########################################################################################
######################################## Schedule3 ########################################
###########################################################################################

SC3 <- c(1,0,0,1,0,0,1,0,0,1,0)

for (s in 1:screening){
  
  Total3 <- subset(Total2,Total2$time == Timing + (s-1))
  
  if (SC3[s]==1){
    
    TT <- assign(paste0('T',s-1), Total3$ID[which(Total3$obs >= DL & Total3$g > Timing+(s-1) & Total3$D >= Timing + (s-1))])
    SS <- assign(paste0('S',s-1), Total3$ID[which(Total3$g == Timing+(s-1) & Total3$D >= Timing + (s-1))])
    
  }else {
    
    SS <- assign(paste0('S',s-1), Total3$ID[which(Total3$g == Timing+(s-1) & Total3$D >= Timing + (s-1))])
    
  }
  
  assign(paste0('U',s-1), union(TT,SS))
  assign(paste0('V',s-1), TT)
  
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

Schedule3[k,] <- c(length(U0),
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

V0  <- union(V0,V0)
V1  <- union(V0,V1)
V2  <- union(V1,V2)
V3  <- union(V2,V3)
V4  <- union(V3,V4)
V5  <- union(V4,V5)
V6  <- union(V5,V6)
V7  <- union(V6,V7)
V8  <- union(V7,V8)
V9  <- union(V8,V9)
V10 <- union(V9,V10)

ScheduleT3[k,] <- c(length(V0),
                    length(V1),
                    length(V2),
                    length(V3),
                    length(V4),
                    length(V5),
                    length(V6),
                    length(V7),
                    length(V8),
                    length(V9),
                    length(V10))/K*100

Tested <- c(length(V0),
            length(V1)-length(V0),
            length(V2)-length(V1),
            length(V3)-length(V2),
            length(V4)-length(V3),
            length(V5)-length(V4),
            length(V6)-length(V5),
            length(V7)-length(V6),
            length(V8)-length(V7),
            length(V9)-length(V8),
            length(V10)-length(V9))

W0  <- setdiff(U0,V0)
W1  <- setdiff(U1,V1)
W2  <- setdiff(U2,V2)
W3  <- setdiff(U3,V3)
W4  <- setdiff(U4,V4)
W5  <- setdiff(U5,V5)
W6  <- setdiff(U6,V6)
W7  <- setdiff(U7,V7)
W8  <- setdiff(U8,V8)
W9  <- setdiff(U9,V9)
W10 <- setdiff(U10,V10)

Symptoms <- c(length(W0),
              length(W1)-length(W0),
              length(W2)-length(W1),
              length(W3)-length(W2),
              length(W4)-length(W3),
              length(W5)-length(W4),
              length(W6)-length(W5),
              length(W7)-length(W6),
              length(W8)-length(W7),
              length(W9)-length(W8),
              length(W10)-length(W9))

ScheduleN3[k,] <- c(Tested[1],Tested[4],Tested[7],Tested[10])
ScheduleS3[k,] <- c(sum(Symptoms[2:4]),sum(Symptoms[5:7]),sum(Symptoms[8:10]),sum(Symptoms[11]))
ScheduleF3[k,] <- c(K,K-length(U10))


###########################################################################################
######################################## Schedule4 ########################################
###########################################################################################

SC4 <- c(1,0,0,0,0,0,0,0,1,1,1)

for (s in 1:screening){
  
  Total3 <- subset(Total2,Total2$time == Timing + (s-1))
  
  if (SC4[s]==1){
    
    TT <- assign(paste0('T',s-1), Total3$ID[which(Total3$obs >= DL & Total3$g > Timing+(s-1) & Total3$D >= Timing + (s-1))])
    SS <- assign(paste0('S',s-1), Total3$ID[which(Total3$g == Timing+(s-1) & Total3$D >= Timing + (s-1))])
    
  }else {
    
    SS <- assign(paste0('S',s-1), Total3$ID[which(Total3$g == Timing+(s-1) & Total3$D >= Timing + (s-1))])
    
  }
  
  assign(paste0('U',s-1), union(TT,SS))
  assign(paste0('V',s-1), TT)
  
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

Schedule4[k,] <- c(length(U0),
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

V0  <- union(V0,V0)
V1  <- union(V0,V1)
V2  <- union(V1,V2)
V3  <- union(V2,V3)
V4  <- union(V3,V4)
V5  <- union(V4,V5)
V6  <- union(V5,V6)
V7  <- union(V6,V7)
V8  <- union(V7,V8)
V9  <- union(V8,V9)
V10 <- union(V9,V10)

ScheduleT4[k,] <- c(length(V0),
                    length(V1),
                    length(V2),
                    length(V3),
                    length(V4),
                    length(V5),
                    length(V6),
                    length(V7),
                    length(V8),
                    length(V9),
                    length(V10))/K*100

Tested <- c(length(V0),
            length(V1)-length(V0),
            length(V2)-length(V1),
            length(V3)-length(V2),
            length(V4)-length(V3),
            length(V5)-length(V4),
            length(V6)-length(V5),
            length(V7)-length(V6),
            length(V8)-length(V7),
            length(V9)-length(V8),
            length(V10)-length(V9))

W0  <- setdiff(U0,V0)
W1  <- setdiff(U1,V1)
W2  <- setdiff(U2,V2)
W3  <- setdiff(U3,V3)
W4  <- setdiff(U4,V4)
W5  <- setdiff(U5,V5)
W6  <- setdiff(U6,V6)
W7  <- setdiff(U7,V7)
W8  <- setdiff(U8,V8)
W9  <- setdiff(U9,V9)
W10 <- setdiff(U10,V10)

Symptoms <- c(length(W0),
              length(W1)-length(W0),
              length(W2)-length(W1),
              length(W3)-length(W2),
              length(W4)-length(W3),
              length(W5)-length(W4),
              length(W6)-length(W5),
              length(W7)-length(W6),
              length(W8)-length(W7),
              length(W9)-length(W8),
              length(W10)-length(W9))

ScheduleN4[k,] <- c(Tested[1],Tested[9],Tested[10],Tested[11])
ScheduleS4[k,] <- c(sum(Symptoms[2:9]),sum(Symptoms[10]),sum(Symptoms[11]),0)
ScheduleF4[k,] <- c(K,K-length(U10))


###########################################################################################
######################################## Schedule5 ########################################
###########################################################################################

# SC5 <- c(1,0,0,0,0,0,0,0,0,0,0)
# random <- sample(2:screening, 3, replace=F)
# SC5[random] <- 1 

SC5 <- c()
for (i in 1:(length(Total2$ID)/screening)){
  
  RD <- c(1,0,0,0,0,0,0,0,0,0,0)
  random <- sort(sample(2:screening, 3, replace=F))
  RD[random[1]] <- 2
  RD[random[2]] <- 3
  RD[random[3]] <- 4
  
  SC5 <- c(SC5,RD)
  
}
Total2$SC5 <- SC5

# Schedule55 <- matrix(0,nrow=testings,ncol=screening)

for (s in 1:screening){
  
  Total3 <- subset(Total2,Total2$time == Timing + (s-1))
  
  TT1 <- assign(paste0('T',s-1), Total3$ID[which(Total3$obs >= DL & Total3$g > Timing+(s-1) & Total3$D >= Timing + (s-1) & Total3$SC5 == 1)])
  TT2 <- assign(paste0('T',s-1), Total3$ID[which(Total3$obs >= DL & Total3$g > Timing+(s-1) & Total3$D >= Timing + (s-1) & Total3$SC5 == 2)])
  TT3 <- assign(paste0('T',s-1), Total3$ID[which(Total3$obs >= DL & Total3$g > Timing+(s-1) & Total3$D >= Timing + (s-1) & Total3$SC5 == 3)])
  TT4 <- assign(paste0('T',s-1), Total3$ID[which(Total3$obs >= DL & Total3$g > Timing+(s-1) & Total3$D >= Timing + (s-1) & Total3$SC5 == 4)])
  
  TT  <- union(TT4,union(TT3,union(TT2,TT1)))
  
  SS <- assign(paste0('S',s-1), Total3$ID[which(Total3$g == Timing+(s-1) & Total3$D >= Timing + (s-1))])
  
  assign(paste0('U',s-1), union(TT,SS))
  assign(paste0('V',s-1), TT)
  
  assign(paste0('A',s-1), TT1)
  assign(paste0('B',s-1), TT2)
  assign(paste0('C',s-1), TT3)
  assign(paste0('D',s-1), TT4)
  
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

Schedule5[k,] <- c(length(U0),
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

V0  <- union(V0,V0)
V1  <- union(V0,V1)
V2  <- union(V1,V2)
V3  <- union(V2,V3)
V4  <- union(V3,V4)
V5  <- union(V4,V5)
V6  <- union(V5,V6)
V7  <- union(V6,V7)
V8  <- union(V7,V8)
V9  <- union(V8,V9)
V10 <- union(V9,V10)

ScheduleT5[k,] <- c(length(V0),
                    length(V1),
                    length(V2),
                    length(V3),
                    length(V4),
                    length(V5),
                    length(V6),
                    length(V7),
                    length(V8),
                    length(V9),
                    length(V10))/K*100

A0  <- union(A0,A0)
A1  <- union(A0,A1)
A2  <- union(A1,A2)
A3  <- union(A2,A3)
A4  <- union(A3,A4)
A5  <- union(A4,A5)
A6  <- union(A5,A6)
A7  <- union(A6,A7)
A8  <- union(A7,A8)
A9  <- union(A8,A9)
A10 <- union(A9,A10)

B0  <- setdiff(union(B0,B0),A10)
B1  <- setdiff(union(B0,B1),A10)
B2  <- setdiff(union(B1,B2),A10)
B3  <- setdiff(union(B2,B3),A10)
B4  <- setdiff(union(B3,B4),A10)
B5  <- setdiff(union(B4,B5),A10)
B6  <- setdiff(union(B5,B6),A10)
B7  <- setdiff(union(B6,B7),A10)
B8  <- setdiff(union(B7,B8),A10)
B9  <- setdiff(union(B8,B9),A10)
B10 <- setdiff(union(B9,B10),A10)

C0  <- setdiff(union(C0,C0),union(B10,A10))
C1  <- setdiff(union(C0,C1),union(B10,A10))
C2  <- setdiff(union(C1,C2),union(B10,A10))
C3  <- setdiff(union(C2,C3),union(B10,A10))
C4  <- setdiff(union(C3,C4),union(B10,A10))
C5  <- setdiff(union(C4,C5),union(B10,A10))
C6  <- setdiff(union(C5,C6),union(B10,A10))
C7  <- setdiff(union(C6,C7),union(B10,A10))
C8  <- setdiff(union(C7,C8),union(B10,A10))
C9  <- setdiff(union(C8,C9),union(B10,A10))
C10 <- setdiff(union(C9,C10),union(B10,A10))

D0  <- setdiff(union(D0,D0),union(C10,union(B10,A10)))
D1  <- setdiff(union(D0,D1),union(C10,union(B10,A10)))
D2  <- setdiff(union(D1,D2),union(C10,union(B10,A10)))
D3  <- setdiff(union(D2,D3),union(C10,union(B10,A10)))
D4  <- setdiff(union(D3,D4),union(C10,union(B10,A10)))
D5  <- setdiff(union(D4,D5),union(C10,union(B10,A10)))
D6  <- setdiff(union(D5,D6),union(C10,union(B10,A10)))
D7  <- setdiff(union(D6,D7),union(C10,union(B10,A10)))
D8  <- setdiff(union(D7,D8),union(C10,union(B10,A10)))
D9  <- setdiff(union(D8,D9),union(C10,union(B10,A10)))
D10 <- setdiff(union(D9,D10),union(C10,union(B10,A10)))

length(V10) == length(A10) + length(B10) + length(C10) + length(D10)

W0  <- setdiff(U0,V0)
W1  <- setdiff(U1,V1)
W2  <- setdiff(U2,V2)
W3  <- setdiff(U3,V3)
W4  <- setdiff(U4,V4)
W5  <- setdiff(U5,V5)
W6  <- setdiff(U6,V6)
W7  <- setdiff(U7,V7)
W8  <- setdiff(U8,V8)
W9  <- setdiff(U9,V9)
W10 <- setdiff(U10,V10)


Symptoms <- c()
for (l in 1:length(W10)){
  
  TotalS <- subset(Total2,Total2$ID==W10[l])
  TotalS <- subset(TotalS,TotalS$time < TotalS$g)
  
  if (length(TotalS$SC5) == 0){
    Symptoms[l] <- 0
  } else {
    Symptoms[l] <- max(TotalS$SC5)
  }

}

ScheduleN5[k,] <- c(length(A10),length(B10),length(C10),length(D10))
ScheduleS5[k,] <- c(length(which(Symptoms==1)),length(which(Symptoms==2)),length(which(Symptoms==3)),length(which(Symptoms==4)))
ScheduleF5[k,] <- c(K,K-length(U10))

}

