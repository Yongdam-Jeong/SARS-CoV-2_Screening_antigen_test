library(ggplot2)
library(ggridges)
library(ggpubr)
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


setwd("~/Desktop")


HeatmapH <- matrix(NA, nrow = 8, ncol = 8)
HeatmapG <- matrix(NA, nrow = 8, ncol = 8)
HeatmapD <- matrix(NA, nrow = 8, ncol = 8)
Heatmap3 <- matrix(NA, nrow = 8, ncol = 8)
Heatmap5 <- matrix(NA, nrow = 8, ncol = 8)




for( iii in 1:8 ) {
  
  
  for( jjj in 1:8 ) {
    
    
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
  R0 <- 1.5 + (iii-1)*0.5 ## Reproduction number
  e1 <- 1/3.209127601 ## Mean latent period for asym
  e2 <- 1/2.745534106 ## Mean latent period for sym
  s1 <- 1/6.032712583 ## Mean infectious period for asym
  s2 <- 1/7.452044219 ## Mean infectious period for sym
  p <- 0.2 + (jjj-1)*0.1 ## Proportion of asym
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
  # PCR <- log10(10^2) ## Detection limit of PCR testing
  PCR <- 0 ## Detection limit of PCR testing
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
    
    if (max(out$V) > IT & out$V[which(out$time==Tmax2)] < PCR) {
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
    
    detected <- rep(max(out$time[out$V>=PCR]), times=length(stime))
    
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
    
    if (max(out$V) >= IT & out$V[which(out$time==Tmax2)] < PCR) {
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
    
    detected <- rep(max(out$time[out$V>=PCR]), times=length(stime))
    
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
  
  
  #################################################################################################
  ##################################### Screening simulations #####################################
  #################################################################################################
  iteration <- 100
  screening <- 10+1
  
  
  ScheduleH <- matrix(NA, nrow = iteration, ncol = screening)
  ScheduleG <- matrix(NA, nrow = iteration, ncol = screening)
  ScheduleD <- matrix(NA, nrow = iteration, ncol = screening)
  Schedule3 <- matrix(NA, nrow = iteration, ncol = screening)
  Schedule5 <- matrix(NA, nrow = iteration, ncol = screening)
  
  
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
    
    
    
    
    ###########################################################################################
    ######################################## ScheduleH ########################################
    ###########################################################################################
    DL <- log10(6.3*(10^4))
    testings  <- 4
    Timing <- Peak
    
    
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
    
    ScheduleH[k,] <- c(length(U0),
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
    ######################################## ScheduleG ########################################
    ###########################################################################################
    DL <- log10(2*(10^6))
    testings  <- 4
    Timing <- Growth
    
    
    #################### Denominator (i.e., target for screening at the initiation)
    Number <- subset(Total,Total$Time==0)  
    Number2 <- subset(Number,Number$time <= Timing & Number$g >= Timing & Number$D >= Timing)
    K <- length(Number2$ID)
    
    if (K > 0) {
    
    #################### COVID-19 cases at t = t_0 (i.e., K(t_0))
    Total1  <- subset(Total,Total$ID==Number2$ID[1])
    for (i in 2:length(Number2$ID)){
      
      Total1 <- rbind(Total1,subset(Total,Total$ID==Number2$ID[i]))
      
    }
    
    Total1 <- subset(Total1,Total1$D >= Timing)
    Total2 <- subset(Total1,Total1$time >= Timing & Total1$time <= Timing+10)  ## Consider only cases during screening programs for reducing computations 
    
    
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
    
    ScheduleG[k,] <- c(length(U0),
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
    
    
    } else {
      
      ScheduleG[k,]  <- rep(100,times=screening)

    }
    
    
    
    
    ###########################################################################################
    ######################################## ScheduleD ########################################
    ###########################################################################################
    DL <- log10(2*(10^6))
    testings  <- 4
    Timing <- Decline
    
    
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
    
    ScheduleD[k,] <- c(length(U0),
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
    ######################################## Schedule3 ########################################
    ###########################################################################################
    DL <- log10(2*(10^6))
    testings  <- 3
    Timing <- Peak
    
    
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
    
    
    SC1 <- c(1,1,1,0,0,0,0,0,0,0,0)
    
    for (s in 1:screening){
      
      Total3 <- subset(Total2,Total2$time == Timing + (s-1))
      
      if (SC1[s]==1){
        
        TT <- assign(paste0('T',s-1), Total3$ID[which(Total3$obs >= DL & Total3$g > Timing+(s-1) & Total3$D >= Timing + (s-1))])
        SS <- assign(paste0('S',s-1), Total3$ID[which(Total3$g == Timing+(s-1) & Total3$D >= Timing + (s-1))])
        
      }else {
        
        SS <- assign(paste0('S',s-1), Total3$ID[which(Total3$g == Timing+(s-1) & Total3$D >= Timing + (s-1))])
        
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
    
    
    
    
    ###########################################################################################
    ######################################## Schedule5 ########################################
    ###########################################################################################
    DL <- log10(2*(10^6))
    testings  <- 5
    Timing <- Peak
    
    
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
    
    
    SC1 <- c(1,1,1,1,1,0,0,0,0,0,0)
    
    for (s in 1:screening){
      
      Total3 <- subset(Total2,Total2$time == Timing + (s-1))
      
      if (SC1[s]==1){
        
        TT <- assign(paste0('T',s-1), Total3$ID[which(Total3$obs >= DL & Total3$g > Timing+(s-1) & Total3$D >= Timing + (s-1))])
        SS <- assign(paste0('S',s-1), Total3$ID[which(Total3$g == Timing+(s-1) & Total3$D >= Timing + (s-1))])
        
      }else {
        
        SS <- assign(paste0('S',s-1), Total3$ID[which(Total3$g == Timing+(s-1) & Total3$D >= Timing + (s-1))])
        
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
    
    
    
  }
  
  
  HeatmapH[iii,jjj] <- mean(ScheduleH[,screening])
  HeatmapG[iii,jjj] <- mean(ScheduleG[,screening])
  HeatmapD[iii,jjj] <- mean(ScheduleD[,screening])
  Heatmap3[iii,jjj] <- mean(Schedule3[,screening])
  Heatmap5[iii,jjj] <- mean(Schedule5[,screening])
  
  
  } #### Proportion
  
  
} #### R0






