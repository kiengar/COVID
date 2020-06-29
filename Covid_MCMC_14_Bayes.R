install.packages("drc")
install.packages("tidyr")
install.packages("zoo")
install.packages("Hmisc")
install.packages("ggplot2")
install.packages("ggrepel")
install.packages("writexl")


library(drc)
library(tidyr)
library(zoo)
library(Hmisc)
library(ggplot2)
library(ggrepel)
library(writexl)


###############################

#Importdata
df_states<-read.csv("C:\\Users\\laslo\\Documents\\R\\COVID-19-master\\Maj1owid-covid-data.csv")
df_testing<-read.csv("C:\\Users\\laslo\\Documents\\R\\COVID-19-master\\Maj1covid-testing-all-observations.csv")

#Exportdata
explibrary <- "C:\\Users\\laslo\\Documents\\R\\COVID-19-master\\"

#Function:calculate mortality rate with Monte Carlo
##########################################################################
inidelay=8
inistdev=1
inideathrate=0.1

#cases= NA,
#deaths= NA,
Country="France"
run=1000
MA.fit=F
MA.cases=F
MA.deaths=F
case_tracking=T
lognormal=F
calculate_CFR_MCMC <- function (#df_new=df_new,
                                      inidelay=8,
                                      inistdev=1,
                                      inideathrate=0.1,
                            
                                      #cases= NA,
                                      #deaths= NA,
                                      Country=NA,
                                      run=1000,
                                      plot=F,
                                      MA.fit=F,
                                      MA.cases=F,
                                      MA.deaths=F,
                                      case_tracking=T,
                                      lognormal=F) {
  
 require(zoo)
  require(gtools)
  start_time=Sys.time()
   Country_cases <- df_new$new_cases[df_new$location==Country]
  Country_cases
  
  if(MA.cases==T){
  smooth_cases<-zoo::rollapply(Country_cases, 3, FUN = mean, fill = NA, align = "center")
  smooth_cases[1]<-Country_cases[1]
  smooth_cases[length(smooth_cases)]<-Country_cases[length(Country_cases)]
  cases <-smooth_cases
  }else{
    cases<- Country_cases
    }
  
 # plot(Country_cases)
  #lines(smooth_cases, col="green")
  
  
  Country_deaths <- df_new$new_deaths[df_new$location==Country]
  which(Country_deaths>0)[1]
  
  if(MA.deaths==T){
  smooth_deaths <-zoo::rollapply(Country_deaths, 3, FUN = mean,
                                 fill = NA,
                                 align = "center")
  smooth_deaths[1]<-Country_deaths[1]
  smooth_deaths[length(smooth_deaths)]<-Country_deaths[length(Country_deaths)]
  
  deaths <- smooth_deaths
  }else{
    deaths<-Country_deaths
  }
  deaths
  
 
  #__________________________________________________________________________________________________________
  calculate_diff <-function(fdelay,
                            fstdev,
                            fdeathrate){

max_possible_delay=length(Country_deaths)-which(Country_deaths>0)[1]-1
    
if (fdelay <= 0 |
    fdelay>max_possible_delay |
    fstdev <= 0 |
    fdeathrate >= 1 |
    fdeathrate <= (-1)){
  Diff=Inf
}else{

    #set a matrix "m_days" where the produced deaths will be recorded   
    m_days <- matrix(NA, nrow = length(cases), ncol = round(max(cases)))
    
if(case_tracking==T) {
#.........................................................................................  
    if(fdeathrate>=0){ 
      for( i in 1:length(cases)){
        loc_sum=sum(rbinom(round(cases[i]),1,fdeathrate))
        if(loc_sum>0){
          for (z in 1:loc_sum){
            if(lognormal==T){
              m_days[i,z]<- round(i+ fdelay + rlnorm(1,sdlog= fstdev))
            }else{
              m_days[i,z]<- round(i+ fdelay + rnorm(1,sd= fstdev))
            }
          }
        }
        df_proposed_deaths <-as.data.frame(table(m_days))
      }
    } else{
      for( i in 1:length(cases)){
        loc_sum=sum(rbinom(round(cases[i]),1,abs(fdeathrate)))
        if(loc_sum>0){
          for (z in 1:loc_sum){
            if(lognormal==T){
              m_days[i,z]<- round(i+ fdelay + rlnorm(1,sdlog= fstdev))
            }else{
              m_days[i,z]<- round(i+ fdelay + rnorm(1,sd= fstdev))
            }
          }
        }
      }
      df_proposed_deaths <-as.data.frame(table(m_days))
      df_proposed_deaths$Freq<- -df_proposed_deaths$Freq
    }  
  #.........................................................................................  
} else{
  #.........................................................................................  
  if(fdeathrate>=0){ 
    for( i in 1:length(cases)){
      loc_sum=round(cases[i]*fdeathrate)
      if(loc_sum>0){
        for (z in 1:loc_sum){
          if(lognormal==T){
            m_days[i,z]<- round(i+ fdelay + rlnorm(1,sdlog= fstdev))
          }else{
            m_days[i,z]<- round(i+ fdelay + rnorm(1,sd= fstdev))
          }
        }
      }
      df_proposed_deaths <-as.data.frame(table(m_days))
    }
  } else{
    for( i in 1:length(cases)){
      loc_sum=round(cases[i]*abs(fdeathrate))
      if(loc_sum>0){
        for (z in 1:loc_sum){
          if(lognormal==T){
            m_days[i,z]<- round(i+ fdelay + rlnorm(1,sdlog= fstdev))
          }else{
            m_days[i,z]<- round(i+ fdelay + rnorm(1,sd= fstdev))
          }
        }
      }
    }
    df_proposed_deaths <-as.data.frame(table(m_days))
    df_proposed_deaths$Freq<- -df_proposed_deaths$Freq
  }
  #........................................................
}
    
    
    
    
    
    
    #df_proposed_deaths <-as.data.frame(table(m_days))
    
    #This produce a data frame with two columns
    #    1) $m_days - days when deaths occuredu
    #    2) $Freq - how many deaths occured
    
    
    if (nrow(df_proposed_deaths)>0){
      df_proposed_deaths$m_days <- as.numeric(as.character(df_proposed_deaths$m_days))
      df_days <- data.frame(m_days=c(1:max(df_proposed_deaths$m_days)))
      df_days_completed<-merge(df_days, df_proposed_deaths, all=T)
      df_days_completed[is.na(df_days_completed)]<-0
      produced_deaths<-df_days_completed$Freq
    }else{produced_deaths=0}
    
    #produced_deaths is a vector showing the proposed deaths each day
    
    #smoothing with moving average    
    if(MA.fit==T){
      smooth_proposed_deaths <-zoo::rollapply(produced_deaths, 3, FUN = mean, fill = NA, align = "center")
      smooth_proposed_deaths[1]<-produced_deaths[1]
      smooth_proposed_deaths[length(smooth_proposed_deaths)]<-produced_deaths[length(produced_deaths)]
      proposed_deaths<-smooth_proposed_deaths
    }else{
      proposed_deaths <-produced_deaths 
    }
    
    #set the length of the produced_deaths vector equal to that of the reference deaths vector    
    MAX<-length(deaths)
    Proposed_V <-rep(0,MAX)
    Proposed_V <-proposed_deaths[1:MAX]
    Reference_V <-as.numeric(as.character(deaths))
    
    if(plot==T){
      plot(Reference_V, main = Country, xlab="days", ylab="deaths", ylim=c(-5,max(Reference_V)))
      points(Proposed_V, col="red")
      Sys.sleep(0.2)}
   
    #calculate the squered error
    diff<-(Reference_V-Proposed_V)^2
    
    Diff <- sum(diff, na.rm=T)
    Diff
}    
    
    return(Diff)
  }
  #_____________________________________________________________________________________________ 
 
  prior = function(delay, stdev, deathrate){
    
    aprior = dunif(delay, min=0, max=length(Country_deaths)-which(Country_deaths>0)[1], log = F)
    bprior = dunif(stdev, min=0, max=5, log = F)
    cprior = dunif(deathrate, min=0, max=1, log = F)
    return(aprior*bprior*cprior) 
  }
  
  delay <- rep(0, run)
  stdev <- rep(0, run)
  deathrate <- rep(0, run)
  
  delay[1] <- inidelay      # this is just a starting value, which I've set arbitrarily to 3
  stdev[1] <- inistdev
  deathrate[1] <-inideathrate
  
  for (i in 2:run) {
    currentdelay <- delay[i - 1]
    currentstdev <- stdev[i - 1]
    currentdeathrate <-deathrate[i - 1]
    
    proposeddelay <- currentdelay + rnorm(1, mean = 0, sd = 1)
    proposedstdev <- currentstdev + rnorm(1, mean = 0, sd = 0.5)
    proposeddeathrate <- currentdeathrate + rnorm(1, mean = 0, sd = 0.04)
    # while (proposeddeathrate<0){
    # proposeddeathrate <- currentdeathrate + rnorm(1, mean = 0, sd = 0.05)
    # }
    
    # proposeddelay <- ifelse(proposeddelay<0, 0, proposeddelay)
    # max_possible_delay=length(Country_deaths)-which(Country_deaths>0)[1]-1
    # proposeddelay <- ifelse(proposeddelay>max_possible_delay, max_possible_delay, proposeddelay)
    # proposedstdev <- ifelse(proposedstdev<0, 0, proposedstdev)
    # proposeddeathrate <- ifelse(proposeddeathrate< (-0.9), -0.9, proposeddeathrate)
    # proposeddeathrate <- ifelse(proposeddeathrate>0.9, 0.9, proposeddeathrate)

    #print(paste(proposeddelay, proposeddeathrate))
    
        
    proposedres <-calculate_diff(proposeddelay,
                   proposedstdev,
                   proposeddeathrate)
    
    proposedprior <-prior(proposeddelay,
                                 proposedstdev,
                                 proposeddeathrate)
    
    currentres <-calculate_diff(currentdelay,
                   currentstdev,
                   currentdeathrate)
    
    currentprior <-prior(currentdelay,
                                currentstdev,
                                currentdeathrate)
    
    print(proposedres)   
    
    A <- proposedres/currentres
 
 
 
  
    if (runif(1) < 1/A & !is.na(1/A)) {
      delay[i] <- proposeddelay
      stdev[i] <- proposedstdev
      deathrate[i] <- proposeddeathrate
      
    } else {
      delay[i] <- currentdelay
      stdev[i] <- currentstdev
      deathrate[i] <- currentdeathrate
    }
   print(i)
  }
  
  return(list(delay=delay, stdev=stdev, deathrate=deathrate))
}
  


####################################################################




#Select countries with more than 50 deaths

Country<-"France"
df_new=df_states

Country_cases <- df_new$new_cases[df_new$location==Country]
Country_cases
Country_deaths <- df_new$new_deaths[df_new$location==Country]
Country_deaths

plot(Country_cases)
lines(Country_deaths, col="green")



  df_new=df_states
 # test_country <-calculate_testing(Countries[index])
  
  d_selCountry<-NULL
  delay_selCountry<-NULL
  stdev_selCountry <-NULL
  
  MR2_selCountry <-NULL
  
  seldates <- df_new$date[df_new$location==Country & df_new$total_deaths>10]
  length(seldates)
  
  
 # for (gg in length(seldates):1){
     gg=46
      df_new=df_states
    df_selCountry <- df_new[df_new$location==Country,]
    removed.rows<-(nrow(df_selCountry)-gg+1):nrow(df_selCountry)
    df_new <-df_selCountry[-c(removed.rows),]
    
   
    cmr <-calculate_CFR_MCMC(Country=Country,
                                   inidelay = 1,
                                   inistdev = 0.1,
                                   inideathrate =0.1,
                                   run=10000,
                                   MA.cases = T,
                                   MA.deaths = T,
                                   MA.fit = F,
                                   case_tracking=F,
                                   lognormal=F,
                                   plot=F
                             )
    
    dwarnings()
    cmr
    plot(cmr$delay, type="l")
    plot(density(cmr$delay), main="CFR distribution", xlab="CFR")
    
    
    plot(cmr$deathrate, type="l")
    plot(density(cmr$deathrate), main="CFD distribution", xlab="CFD")
    mean(cmr$deathrate [100:1000])
    median((cmr$deathrate [100:1000]))
    
    getmode <- function(v) {
      uniqv <- unique(v)
      uniqv[which.max(tabulate(match(v, uniqv)))]
    }
    
    getmode(cmr$deathrate [500:5000])
    quantile(cmr$deathrate, 0.975)
    quantile(cmr$deathrate, 0.025)
    
    
    plot(cmr$delay)
    plot(density(cmr$delay))
    mean(cmr$delay [500:10000])
    
    plot(cmr$stdev)
    plot(density(cmr$stdev))
    mean(cmr$stdev [500:10000])
    
    
    d_selCountry[gg]=cmr[[3]]
    delay_selCountry[gg]=cmr[[1]]
    stdev_selCountry[gg]=cmr[[2]]
    MR2_selCountry[gg]<-calculate_MR2(Countries[index])
    testing_selCountry[gg]=test_country
    print(paste0(Countries[index],gg))
#  }
  
  
  ggg <- data.frame(Country=Countries[index],
                             Date=rev(seldates),
                             Days=0:-(length(seldates)-1),
                            deathrate=d_selCountry,
                             MR2=MR2_selCountry,
                             delay=delay_selCountry,
                            stdev=stdev_selCountry,
                            testing=testing_selCountry)
  


  ggg2[[index]]<-ggg[[index]][!is.na(ggg[[index]]$deathrate),]
  jpeg(paste0(explibrary, Countries[[index]],"_Mort.jpg"),
       width = 320, height = 240)
  plot(x=ggg2[[index]]$Days, y=ggg2[[index]]$deathrate, ylim=c(0,0.25),
       xlab="Days", ylab = "Mortality rate", main=paste0(Countries[[index]],"Mortality calculations"))
  plx<-predict(loess(ggg2[[index]]$deathrate ~ ggg2[[index]]$Days), se=T)
  
  lines(ggg2[[index]]$Days,plx$fit, col="blue")
  lines(ggg2[[index]]$Days,plx$fit - qt(0.975,plx$df)*plx$se, lty=2, col="blue")
  lines(ggg2[[index]]$Days,plx$fit + qt(0.975,plx$df)*plx$se, lty=2, col="blue")
  points(ggg2[[index]]$Days, ggg2[[index]]$MR2, col="red")
 
  dev.off()
  
  plot(x=ggg2[[index]]$Days, y=ggg2[[index]]$deathrate, ylim=c(0,0.25),
       xlab="Days", ylab = "Mortality rate", main=paste0("Mortality ",Countries[[index]]))
  plx<-predict(loess(ggg2[[index]]$deathrate ~ ggg2[[index]]$Days), se=T)
  
  lines(ggg2[[index]]$Days,plx$fit, col="blue")
  lines(ggg2[[index]]$Days,plx$fit - qt(0.975,plx$df)*plx$se, lty=2, col="blue")
  lines(ggg2[[index]]$Days,plx$fit + qt(0.975,plx$df)*plx$se, lty=2, col="blue")
  points(ggg2[[index]]$Days, ggg2[[index]]$MR2, col="red")
  
  ggg2[[index]]$loess.fit<-plx$fit
  ggg2[[index]]$loess.upperci<-plx$fit + qt(0.975,plx$df)*plx$se
  ggg2[[index]]$loess.lowerci<-plx$fit - qt(0.975,plx$df)*plx$se
  
 



vvv <- Reduce(function(x, y) rbind(x, y), ggg2, accumulate=F)




writexl::write_xlsx(vvv, paste0(explibrary,"latszolagosmotalitastabla5.xlsx"))

################################################Germany


########################################################################################
#Korrelacio
cor(vvv$testing, vvv$deathrate, method="spearman")
rcorr(as.matrix(vvv[,2:3]), type="spearman") 
#Significant correlation


#Korrelacio fitting
library(drc)
model <-drm(vvv$deathrate ~ vvv$testing, fct=EXD.3())

summary(model)
model$coefficients

predict(model)

plot(model,xlab="Testing capacity", ylab = "Apparent death rate",
     main="Correlation")




#f(x) = c + (d-c)(\exp(-x/e))
#death rate if everybody would be tested
#c:(Intercept)   2.1673 

demo.fits <- expand.grid(conc=exp(seq(log(0.01), log(250), length=1000))) 

warnings()
pm <- predict(model, newdata=demo.fits, interval="confidence") 
demo.fits$p <- pm[,1]
demo.fits$pmin <- pm[,2]
demo.fits$pmax <- pm[,3]

library(ggrepel)
ggplot(vvv, aes(x = testing, y = deathrate)) +
  geom_point() +
  #geom_ribbon(data=demo.fits, aes(x=conc, y=p, ymin=pmin, ymax=pmax), alpha=0.2) +
  geom_line(data=demo.fits, aes(x=conc, y=p)) +
  ggtitle("Correlation")+
  geom_text_repel(aes(label=Country),hjust=0, vjust=0)


calculate_diff(16.3, 0.1, -0.004)

dunif(-0.1, min=0, max=1, log = F)

prior(2, 1, -0.1)

0/50
1/0
1/1000
1/0.001

plot=F

calculate_diff <-function(fdelay,
                     fstdev,
                     fdeathrate){
  
  
  #set a matrix "m_days" where the produced deaths will be recorded   
  m_days <- matrix(NA, nrow = length(cases), ncol = round(max(cases)))
  
  if(case_tracking==T) {
    #.........................................................................................  
    if(fdeathrate>=0){ 
      for( i in 1:length(cases)){
        loc_sum=sum(rbinom(round(cases[i]),1,fdeathrate))
        if(loc_sum>0){
          for (z in 1:loc_sum){
            if(lognormal==T){
              m_days[i,z]<- round(i+ fdelay + rlnorm(1,sdlog= fstdev))
            }else{
              m_days[i,z]<- round(i+ fdelay + rnorm(1,sd= fstdev))
            }
          }
        }
        df_proposed_deaths <-as.data.frame(table(m_days))
      }
    } else{
      for( i in 1:length(cases)){
        loc_sum=sum(rbinom(round(cases[i]),1,abs(fdeathrate)))
        if(loc_sum>0){
          for (z in 1:loc_sum){
            if(lognormal==T){
              m_days[i,z]<- round(i+ fdelay + rlnorm(1,sdlog= fstdev))
            }else{
              m_days[i,z]<- round(i+ fdelay + rnorm(1,sd= fstdev))
            }
          }
        }
      }
      df_proposed_deaths <-as.data.frame(table(m_days))
      df_proposed_deaths$Freq<- -df_proposed_deaths$Freq
    }  
    #.........................................................................................  
  } else{
    #.........................................................................................  
    if(fdeathrate>=0){ 
      for( i in 1:length(cases)){
        loc_sum=round(cases[i]*fdeathrate)
        if(loc_sum>0){
          for (z in 1:loc_sum){
            if(lognormal==T){
              m_days[i,z]<- round(i+ fdelay + rlnorm(1,sdlog= fstdev))
            }else{
              m_days[i,z]<- round(i+ fdelay + rnorm(1,sd= fstdev))
            }
          }
        }
        df_proposed_deaths <-as.data.frame(table(m_days))
      }
    } else{
      for( i in 1:length(cases)){
        loc_sum=round(cases[i]*abs(fdeathrate))
        if(loc_sum>0){
          for (z in 1:loc_sum){
            if(lognormal==T){
              m_days[i,z]<- round(i+ fdelay + rlnorm(1,sdlog= fstdev))
            }else{
              m_days[i,z]<- round(i+ fdelay + rnorm(1,sd= fstdev))
            }
          }
        }
      }
      df_proposed_deaths <-as.data.frame(table(m_days))
      df_proposed_deaths$Freq<- -df_proposed_deaths$Freq
    }
    #........................................................
  }
  
  
  
  
  
  
  #df_proposed_deaths <-as.data.frame(table(m_days))
  
  #This produce a data frame with two columns
  #    1) $m_days - days when deaths occuredu
  #    2) $Freq - how many deaths occured
  
  
  if (nrow(df_proposed_deaths)>0){
    df_proposed_deaths$m_days <- as.numeric(as.character(df_proposed_deaths$m_days))
    df_days <- data.frame(m_days=c(1:max(df_proposed_deaths$m_days)))
    df_days_completed<-merge(df_days, df_proposed_deaths, all=T)
    df_days_completed[is.na(df_days_completed)]<-0
    produced_deaths<-df_days_completed$Freq
  }else{produced_deaths=0}
  
  #produced_deaths is a vector showing the proposed deaths each day
  
  #smoothing with moving average    
  if(MA.fit==T){
    smooth_proposed_deaths <-zoo::rollapply(produced_deaths, 3, FUN = mean, fill = NA, align = "center")
    smooth_proposed_deaths[1]<-produced_deaths[1]
    smooth_proposed_deaths[length(smooth_proposed_deaths)]<-produced_deaths[length(produced_deaths)]
    proposed_deaths<-smooth_proposed_deaths
  }else{
    proposed_deaths <-produced_deaths 
  }
  
  #set the length of the produced_deaths vector equal to that of the reference deaths vector    
  MAX<-length(deaths)
  Proposed_V <-rep(0,MAX)
  Proposed_V <-proposed_deaths[1:MAX]
  Reference_V <-as.numeric(as.character(deaths))
  
  #calculate the squered error
  diff<-(Reference_V-Proposed_V)^2
  
  Diff <- sum(diff, na.rm=T)
  Diff
  
  total_proposed_deaths <- sum(Proposed_V, na.rm=T)
  return(Diff)
}

calculate_diff(-1, 9, 0.1)
rbinom (10,1,0.1)
Inf/Inf
Inf/2
0/0
1/0
