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

#USER DEFINIED FUNCTIONS***********************************************************************
#Function:calculate mortality rate with Monte Carlo
##########################################################################
calculate_mortality_rate <- function (#df_new=df_new,
                                      inidelay=8,
                                      inistdev=1,
                                      inideathrate=0.1,
                                      #cases= NA,
                                      #deaths= NA,
                                      lognormal=F,
                                      Country=NA,
                                      MA.cases=T,
                                      MA.deaths=T,
                                      MA.fit=T,
                                      case_tracking=F,
                                      run=100,
                                      plot=T) {
  
 require(zoo)
  Country_cases <- df_new$new_cases[df_new$location==Country]
  Country_cases <- ifelse(Country_cases<0,0,Country_cases)
  Country_cases
  
  if(MA.cases==T){
    smooth_cases<-zoo::rollapply(Country_cases, 3, FUN = mean, fill = NA, align = "center")
    smooth_cases[1]<-Country_cases[1]
    smooth_cases[length(smooth_cases)]<-Country_cases[length(Country_cases)]
    cases <-smooth_cases
  }else{
    cases<-Country_cases
  }
  cases

  
  Country_deaths <- df_new$new_deaths[df_new$location==Country]
  Country_deaths <- ifelse(Country_deaths<0,0,Country_deaths)
  Country_deaths
  
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
  
  if(sum(cases)>10000){
    cases=round(cases/10)
    deaths=round(deaths/10)
  }
  
  if(sum(cases)>100000){
    cases=round(cases/100)
    deaths=round(deaths/100)
  }
  
  if(sum(cases)>1000000){
    cases=round(cases/1000)
    deaths=round(deaths/1000)
  }

  delay=inidelay
  stdev=inistdev
  deathrate=inideathrate
  
  delaynew=inidelay
  stdevnew=inistdev
  deathratenew=inideathrate
  
  Diff <-rep(NA, run)
  Delay<-rep(NA,run)
  Stdev <-rep(NA,run)
  Deathrate <-rep(NA,run)
  Reference_list <-NULL
  Proposed_list <-NULL
  
  #LOOP__________________________________________________________________________________
  for (kk in 1:run){
    delay<-delaynew
    stdev <- stdevnew
    deathrate <- deathratenew

    #set a matrix "m_days" where the produced deaths will be recorded   
    m_days <- matrix(NA, nrow = length(cases), ncol = round(max(cases)))

  if(case_tracking==T) {
    for( i in 1:length(cases)){
      loc_sum=sum(rbinom(round(cases[i]),1,deathrate))
      if(loc_sum>0){
        for (z in 1:loc_sum){
          if(lognormal==T){
            m_days[i,z]<- round(i+ delay + rlnorm(1,sdlog= stdev))
          }else{
            m_days[i,z]<- round(i+ delay + rnorm(1,sd= stdev))
          }
        }
      }
    }
  } else{
    for( i in 1:length(cases)){
      if(round(cases[i]*deathrate)>0){
        for (z in 1:round(cases[i]*deathrate)){
          if(lognormal==T){
            m_days[i,z]<- round(i+ delay + rlnorm(1,sdlog= stdev))
          }else{
            m_days[i,z]<- round(i+ delay + rnorm(1,sd= stdev))
          }
        }
      }
    }
  }  

    df_proposed_deaths <-as.data.frame(table(m_days))

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

    # keep the vectors if plot is required
    if(plot==T){
      Reference_list[[kk]]<- Reference_V
      Proposed_list[[kk]] <- Proposed_V
    }
  
    #calculate the squered error
    diff<-(Reference_V-Proposed_V)^2
    
    #store the sum of squered
    Diff[kk] <- sum(diff, na.rm=T)
    
    Delay[kk]<-delay
    Stdev[kk] <- stdev
    Deathrate[kk]<- deathrate
    
    #go back to the best setting
    if (kk>1 & Diff[kk] > min(Diff, na.rm = T)){
      delay=Delay[which.min(Diff)]
      stdev=Stdev[which.min(Diff)]
      deathrate=Deathrate[which.min(Diff)]
    }
    
    #find new parameters- first check which parametere shouldbe modified
    x1 <- ceiling(runif(1, 0, 3))
    x1
    if (x1==1){
      # delaynew=delay-2+ceiling(runif(1, 0, 3))
      delaynew=runif(1, delay-0.2*delay, delay+0.2*delay )
      delaynew =ifelse(delaynew<0, 0, delaynew)
      stdevnew=stdev
      deathratenew=deathrate
    }
    
    #modify the selected parameter
    if (x1==2){
      stdevnew=runif(1, stdev-0.2*stdev, stdev+0.2*stdev )
      delaynew=delay
      deathratenew=deathrate
    }
    
    if (x1==3){
      stdevnew=stdev
      delaynew=delay
      deathratenew=runif(1, deathrate-0.2*deathrate, deathrate+0.2*deathrate )
    }
    #print(kk)
    
    
  } #LOOP END

  if(plot==T) {
    jpeg(paste0(explibrary, Country,".jpg"),
         width = 320, height = 240)
    plot(Reference_list[[which.min(Diff)]], main = Country, xlab="days", ylab="deaths")
    points(Proposed_list[[which.min(Diff)]], col="red")
    
    dev.off()
    
    plot(Reference_list[[which.min(Diff)]], main = Country, xlab="days", ylab="deaths")
    points(Proposed_list[[which.min(Diff)]], col="red")
    
  } 
 
  #Export data from the function 
  mylist <-list(mCFD=Delay[which.min(Diff)],
                sdCFD=Stdev[which.min(Diff)],
                CFR=Deathrate[which.min(Diff)])
 
  return(mylist)
} 

####################################################################

#Function:Calculate testing
######################################################################

calculate_testing_bydate <-function (Country=Country
                              #,df_new=df_new,
                              #df_testing=df_testing
                              )
  
{
  
  df_testing2=df_testing[max(as.Date(df_new$date,"%d/%m/%Y")) == as.Date(df_testing$Date,"%d/%m/%Y"),]
  
  if(length(grep(Country, df_testing2$Entity))==0){
    postpercases_Country<-NA} else{
    if(length(grep(Country, df_testing2$Entity))>1){
       country_pos<-c(grep(Country, df_testing2$Entity))
       } else{
       country_pos<-grep(Country, df_testing2$Entity)
       }
  
  tpc_Country <-  min(df_testing2$Cumulative.total[country_pos])
  
  casemax=max(df_new$total_cases[df_new$location==Country])
  postpercases_Country=tpc_Country/casemax}

  return(postpercases_Country)
}
########################################################################

#Function:Calculate mortaity rate without delay
###################################################################
calculate_MR2 <-function(Country=Country){
  Country_cases2 <- max(df_new$total_cases[df_new$location==Country])
  Country_deaths2 <- max(df_new$total_deaths[df_new$location==Country])
  crudeCFR=Country_deaths2/Country_cases2 
  return(crudeCFR)
}
#######################################################################





########################################################################################
  #Select countries with more than 50 deaths

  df_selectedcountries <- df_states[df_states$total_deaths>50, ] 
  Countries<- unique(df_selectedcountries$location)
  Countries

  #Calculate parameters from the selected countries
  list_results<-NULL

  #LOOP1_for countries
  for (index in  1:length(Countries)) {
  
  df_new=df_states
  
  CFR_selCountry<-NULL
  mCFD_selCountry<-NULL
  sdCFD_selCountry <-NULL
  test_selCountry <-NULL
  crudeCFR_selCountry <-NULL
  testing_selCountry <- NULL

  #select dates occuring at least 5 deaths in total  
  seldates <- df_new$date[df_new$location==Countries[index] & df_new$total_deaths>4]
  seldates
  length(seldates)
  
    # 2nd LOOP START____for each days(time course)
    for (gg in length(seldates):1){
  
      #set the dataframe   
      df_new=df_states
      df_selCountry <- df_new[df_new$location==Countries[index],]
      removed.rows<-(nrow(df_selCountry)-gg+1):nrow(df_selCountry)
      df_new <-df_selCountry[-c(removed.rows),]
    
      #calculate testing
      test_country <-calculate_testing_bydate(Countries[index])
      test_country

      #calculate CFR, mCFD, sdCFD
      result <-calculate_mortality_rate(Country=Countries[index],
                                   inidelay = 5,
                                   inistdev = 1,
                                   inideathrate =0.1,
                                   run=10000,
                                   lognormal=F,
                                   MA.cases=T,
                                   MA.deaths=T,
                                   MA.fit=T,
                                   case_tracking=T,
                                   plot=T)
    
      CFR_selCountry[gg]=result[[3]]
      mCFD_selCountry[gg]=result[[1]]
      sdCFD_selCountry[gg]=result[[2]]
   
      #add crude CFR 
      crudeCFR_selCountry[gg]<-calculate_MR2(Countries[index])
      testing_selCountry[gg]=test_country
    
      print(paste0(Countries[index],gg))
    }#__2nd LOOP END
  
  list_results[[index]] <- data.frame(Country=Countries[index],
                             Date=rev(seldates),
                             Days=0:-(length(seldates)-1),
                             CFR=CFR_selCountry,
                             crudeCFR=crudeCFR_selCountry,
                             mCFD=mCFD_selCountry,
                             sdCFD=sdCFD_selCountry,
                             testing=testing_selCountry)
  
  #Save a plot
  jpeg(paste0(explibrary, Countries[[index]],"_Mort.jpg"),
       width = 320, height = 240)
  plot(x=list_results[[index]]$Days, y=list_results[[index]]$CFR, ylim=c(0,0.25),
       xlab="Days", ylab = "Mortality rate", main=paste0("Mortality ", Countries[[index]]))
  plx<-predict(loess(list_results[[index]]$CFR ~ list_results[[index]]$Days), se=T)
  
  lines(list_results[[index]]$Days,plx$fit, col="blue")
  lines(list_results[[index]]$Days,plx$fit - qt(0.975,plx$df)*plx$se, lty=2, col="blue")
  lines(list_results[[index]]$Days,plx$fit + qt(0.975,plx$df)*plx$se, lty=2, col="blue")
  points(list_results[[index]]$Days, list_results[[index]]$crudeCFR, col="red")
 
  dev.off()
  
  plot(x=list_results[[index]]$Days, y=list_results[[index]]$CFR, ylim=c(0,0.25),
       xlab="Days", ylab = "Mortality rate", main=paste0("Mortality ",Countries[[index]]))
  plx<-predict(loess(list_results[[index]]$CFR ~ list_results[[index]]$Days), se=T)
  
  lines(list_results[[index]]$Days,plx$fit, col="blue")
  lines(list_results[[index]]$Days,plx$fit - qt(0.975,plx$df)*plx$se, lty=2, col="blue")
  lines(list_results[[index]]$Days,plx$fit + qt(0.975,plx$df)*plx$se, lty=2, col="blue")
  points(list_results[[index]]$Days, list_results[[index]]$crudeCFR, col="red")
  
  list_results[[index]]$loess.fit<-plx$fit
  list_results[[index]]$loess.upperci<-plx$fit + qt(0.975,plx$df)*plx$se
  list_results[[index]]$loess.lowerci<-plx$fit - qt(0.975,plx$df)*plx$se

}

#unite result dataframes
df_full <- Reduce(function(x, y) rbind(x, y), list_results, accumulate=F)

#Save the full dataframe
writexl::write_xlsx(df_full, paste0(explibrary,"MotalitastablaCTr.xlsx"))

#####################################################################################################


########################################################################################
#Import the full table
library("readxl")
df_full <-read_excel("C:\\Users\\laslo\\Documents\\R\\COVID-19-master\\MotalitastablaCTr.xlsx")
df_full


#Fill the missing testing data if it is possible

States<- unique(df_full$Country)
States="Australia"

for (i in 1:length(States)){
  selCountry<-df_full[df_full$Country==States[i],]
  selCountry
  if(  all(is.na(selCountry$testing))==F){
  plot(x=selCountry$Days, y=selCountry$testing, main=States[i])
    }
}

df_full_removed <- na.omit(df_full)
#select a date
df_zeroDay<-df_full[df_full$Days==0,]


#Korrelacio
plot(y=df_full_removed$CFR, x=df_full_removed$testing,
                col=as.factor(df_full_removed$Country), log = "x")

names(df_full_removed)
ggplot(df_full_removed)+
  aes(x=testing, y=CFR, col=Country)+
  geom_point()+
  scale_x_continuous(trans = 'log10') 
  


hcor(df_full_removed$testing, df_full_removed$crude_CFR, method="spearman")
rcorr(as.matrix(df_full_removed[,c(4,8)]), type="spearman") 
#Significant correlation


#Korrelacio fitting
library(drc)
model <-drm(df_full_removed$CFR ~ df_full_removed$testing, fct=EXD.3())

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


