#CAROLINE RICHARDSON
#13358846
#4BCT

library(deSolve)
library(ggplot2)
library(FME)
library(plyr)


runsim <- function(rvec){

  #browser()
  # this is the model function
  print(sprintf("Simulation run number %d",rvec[["RunNo"]]))
  
  model <- function(time, stocks, auxs){
    with(as.list(c(stocks, auxs)),{ 
      
      aEffective.Contact.Rate <- aContactRate * aInfectivity
      aBeta <- aEffective.Contact.Rate / aTotalPopulation
      aLambda <- aBeta * sInfected
      
      #flows
      fIR <- sSusceptible * aLambda
      fRR <- sInfected / aDelay
      fVR <- sSusceptible * aVaccFraction
      
      #stocks
      dS_dt  <- -fIR - fVR
      dI_dt  <- fIR - fRR
      dR_dt  <- fRR + fVR
      
      return (list(c(dS_dt,dI_dt,dR_dt),
                   IR=fIR, RR=fRR, VR=fVR,
                   Beta=aBeta,Lambda=aLambda,DEL=aDelay,
                   Effective.Contact.Rate=aEffective.Contact.Rate,VaccFraction=aVaccFraction,
                   Infectivity=aInfectivity, InitInfected = initInfected))   
    })
  }
  
  # setup the individual simulation run
  START<-0; FINISH<-20; STEP<-0.01; 
  simtime <- seq(START, FINISH, by=STEP)
  
  init<-rvec[["initInfected"]]
  
  a<-c( aContactRate=4, aTotalPopulation=10000, aDelay=2, rvec["aInfectivity"],
       rvec["aVaccFraction"], rvec["initInfected"])
  
  stocks  <- c(sSusceptible=10000-init, sInfected=init, sRecovered=0)
  
  o<-data.frame(ode(y=stocks, simtime, func = model, 
                    parms=a, method="euler"))
  o$RunNumber<-rvec["RunNo"]
  o
}


aInfectivity.MIN <- 0;     aInfectivity.MAX <- 1.0;
aVaccFraction.MIN <- 0.0;  aVaccFraction.MAX <- 1.0;
initInfected.MIN <- 1.0;   initInfected.MAX <- 25.0;

parRange<-data.frame(
  min=c(aInfectivity.MIN, aVaccFraction.MIN, initInfected.MIN),
  max=c(aInfectivity.MAX, aVaccFraction.MAX, initInfected.MAX)
)


rownames(parRange)<-c("aInfectivity","aVaccFraction","initInfected")

set.seed(1234)

#200 runs
NRUNS <- 200

p<-data.frame(RunNo=1:NRUNS,Latinhyper(parRange,NRUNS))


out<-apply(p,1,function(x){
  #browser()
  df <- runsim(x)
  df
})

out2<-lapply(out, function(x){
  x[which.max(x$sInfected), c("sInfected", "VaccFraction", "Infectivity")]
})


df<-rbind.fill(out)

df2<-rbind.fill(out2) 


#PLOTTING

ggplot(df,aes(x=time,y=sInfected,color=RunNumber)) + 
  geom_path() + scale_colour_gradientn(colours=rainbow(10))+
  ylab("Infected") +
  xlab("Time (Days)") + guides(color=FALSE)

ggplot(df2,aes(x=VaccFraction, y=Infectivity)) +
  ggtitle("Sensitivity Analysis") +
  xlab("VF")+
  ylab("INF")+
  geom_point(aes(size=sInfected, color=sInfected )) +
  scale_color_gradientn(name="sInfected", colours=rainbow(10)) +
  scale_size(range=c(3,10))
