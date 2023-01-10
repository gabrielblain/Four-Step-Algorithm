#Table S1 R-code for calculating the four-Step algorithm
rm(list = ls())
if (! require (gamlss)) {install.packages ("gamlss")}
if (! require (gamlss.dist)) {install.packages ("gamlss.dist")}
library(gamlss)
library(gamlss.dist)
#setwd("C:/mydocuments/") #setting the work directory
data.matrix=read.csv("mydata.csv", header=TRUE,sep=",")
#maydata.csv is 3-column a csv file (the input file).
#1st column is the years;
#2nd is the months (1 to 12)
#3rd is the monthly rainfall amounts
#For example
#Years Month Rain
#1900 1 150
#[...] [...] [...]
#2020 12 200
mumber.time.scales=12 # Must range between 1 and 12
n.locals=length(data.matrix[1,])-2
for (ts in 1:mumber.time.scales){
  time.scale=ts
  suitable.dist.deltaAICc=matrix(NA,12,n.locals)
  N=length(data.matrix[,1])
  moderate=matrix(NA,12,n.locals)
  severe=matrix(NA,12,n.locals)
  extreme=matrix(NA,12,n.locals)
  for (local in 1:n.locals){
    rain=as.matrix(data.matrix[,(local+2)])
    colnames(rain)=c("rain")
    year=data.matrix[,1]
    month=data.matrix[,2]
    data.orig=as.data.frame(cbind(year,month,rain))
    nt=length(data.orig$month)
    data=matrix(NA,(nt-12),3)
    rain=matrix(NA,nt,1)
    for (end in time.scale:(nt)){
      rain[end,1]=(sum(data.orig$rain[(end-(time.scale-1)):end]))}
    data=na.omit(as.data.frame(cbind(data.orig$year, data.orig$month,
                                     rain)))
    colnames(data)=c("year","month","rain")
    data=data[order(data$month),]
    for (month in 1:12){
      i.month=which(data$month==month)
      data.by.month=as.matrix(na.omit(data[min(i.month):max(i.month),]))
      n=length(data.by.month[,1])
      time=as.matrix(seq(1:n))
      x.axis=as.matrix(qnorm(ppoints(n)))
      SPI.gamma=matrix(NA,n,1)
      probacum.gamma=matrix(NA,n,1)
      non.rain.data.boolean=data.by.month[,3]>0
      id=which(data.by.month[,3]>0)
      time.nonzero=as.matrix(time[id])
      non.rain.data=as.matrix(data.by.month[,3][non.rain.data.boolean])
      np=length(non.rain.data); nz=n-np; data.by.month[,3][data.by.month[,
                                                                         3]<0]=0
      probzero=(nz+1)/(2*(n+1))
      t.gam=gamlss(non.rain.data~1,family=GA, mu.link = "identity",
                   sigma.link ="log")
      t.gam.ns10=gamlss(non.rain.data~poly(time.nonzero,1),family=GA,
                        mu.link = "identity", sigma.link ="log")
      t.gam.ns11=gamlss(non.rain.data~poly(time.nonzero,1), sigma.formula
                        =~poly(time.nonzero,1), family=GA, mu.link = "identity", sigma.link ="log")
      probacum.gamma=as.matrix(probzero+(1-probzero)*pGA(data.by.month[,3],
                                                         mu = t.gam$mu.fv[1], sigma = t.gam$sigma.fv[1], lower.tail = TRUE, log.p =
                                                           FALSE))
      for (t in 1:n){
        if (probacum.gamma[t,1]<=0.5){sinal=sqrt(log(1/probacum.gamma[t,
                                                                      1]^2))
        SPI.gamma[t,1]=-1*(sinal-
                             (2.515517+0.802853*sinal+0.010328*sinal^2)/
                             (1+1.432788*sinal+0.189269*sinal^2+0.001308*sinal^3))}
        else{sinal=sqrt(log(1/(1-probacum.gamma[t,1])^2))
        SPI.gamma[t,1]=sinal-(2.515517+0.802853*sinal+0.010328*sinal^2)/
          (1+1.432788*sinal+0.189269*sinal^2+0.001308*sinal^3)}
      }
      AICc.gamma=AIC(t.gam, k=2, c=TRUE)
      AICc.nsgamma10=AIC(t.gam.ns10, k=2, c=TRUE)
      AICc.nsgamma11=AIC(t.gam.ns11, k=2, c=TRUE)
      akaikec=cbind(AICc.gamma,AICc.nsgamma10,AICc.nsgamma11)
      delta.akaikec=as.data.frame(akaikec-akaikec[which.min(akaikec)])
      best.dist.AICc=delta.akaikec[which(delta.akaikec==0)]
      suitable.dist.AICc=colnames(best.dist.AICc)
      LRT=(1-pchisq(-(t.gam.ns10$P.deviance-t.gam$P.deviance),1))
      if (LRT>0.05){suitable.dist.LRT="model1"}else{
        LRT=(1-pchisq(-(t.gam.ns11$P.deviance-t.gam.ns10$P.deviance),1))
        if(LRT>0.05){suitable.dist.LRT="model2"}
        else{suitable.dist.LRT="model3"}}
      delta.test.AICc=which(delta.akaikec<=2)
      if (length(delta.test.AICc)==3)
      {suitable.dist.deltaAICc[month,local]=suitable.dist.LRT}
      if (length(delta.test.AICc)==1){
        if (suitable.dist.AICc=="AICc.gamma")
        {suitable.dist.deltaAICc[month,local]="model1"}
        if (suitable.dist.AICc=="AICc.nsgamma10")
        {suitable.dist.deltaAICc[month,local]="model2"}
        if (suitable.dist.AICc=="AICc.nsgamma11")
        {suitable.dist.deltaAICc[month,local]="model3"}}
      if (length(delta.test.AICc)==2){
        if (delta.test.AICc[1]==1 && delta.test.AICc[2]==2){
          if (t.gam.ns10$G.deviance<t.gam$G.deviance){
            LRT=LR.test(t.gam,t.gam.ns11, print="FALSE")
            if (LRT$p.val>=0.05){suitable.dist.deltaAICc[month,local]="model1"}
            else{suitable.dist.deltaAICc[month,local]="model2"}}
          else{suitable.dist.deltaAICc[month,local]="model1"}}
        if (delta.test.AICc[1]==1 && delta.test.AICc[2]==3){
          if (t.gam.ns11$G.deviance<t.gam$G.deviance){
            LRT=LR.test(t.gam,t.gam.ns11, print="FALSE")
            if (LRT$p.val>=0.05)
            {suitable.dist.deltaAICc[month,local]="model1"}
            else{suitable.dist.deltaAICc[month,local]="model3"}}
          else{suitable.dist.deltaAICc[month,local]="model1"}}
        if (delta.test.AICc[1]==2 && delta.test.AICc[2]==3){
          if (t.gam.ns11$G.deviance<t.gam.ns10$G.deviance){
            LRT=LR.test(t.gam.ns10,t.gam.ns11, print="FALSE")
            if (LRT$p.val>=0.05)
            {suitable.dist.deltaAICc[month,local]="model2"}
            else{suitable.dist.deltaAICc[month,local]="model3"}}
          else{suitable.dist.deltaAICc[month,local]="model2"}}}
      p.moderate.drought=(pnorm(-1, mean=0, sd=1))
      p.severe.drought=(pnorm(-1.5, mean=0, sd=1))
      p.extreme.drought=(pnorm(-2, mean=0, sd=1))
      p.chuva.moderate=(p.moderate.drought-probzero)/(1-probzero)
      p.chuva.severe=(p.severe.drought-probzero)/(1-probzero)
      p.chuva.extreme=(p.extreme.drought-probzero)/(1-probzero)
      if(p.chuva.moderate>0){
        moderate.rain=qGA(p.chuva.moderate,mu = t.gam$mu.fv[1], sigma =
                            (t.gam$sigma.fv[1]),lower.tail = TRUE, log.p = FALSE)}
      else{moderate.rain=0}
      if(p.chuva.severe>0){
        severe.rain=qGA(p.chuva.severe,mu = t.gam$mu.fv[1], sigma =
                          t.gam$sigma.fv[1],lower.tail = TRUE, log.p = FALSE)}
      else{severe.rain=0}
      if(p.chuva.extreme>0){
        extreme.rain=qGA(p.chuva.extreme,mu = t.gam$mu.fv[1], sigma =
                           t.gam$sigma.fv[1],lower.tail = TRUE, log.p = FALSE)}
      else{extreme.rain=0}
      if (suitable.dist.deltaAICc[month,local]=="model1"){
        actual.prob.moderate=100*pnorm(-1, mean=0, sd=1, lower.tail=TRUE)
        actual.prob.severe=100*pnorm(-1.5, mean=0, sd=1, lower.tail=TRUE)
        actual.prob.extreme=100*pnorm(-2, mean=0, sd=1, lower.tail=TRUE)
        mu.trend=0; sigma.trend=0
      }
      if (suitable.dist.deltaAICc[month,local]=="model2"){
        f=length(t.gam.ns10$mu.fv)
        timevarying.moderate=as.matrix(pGA(moderate.rain, mu =
                                             t.gam.ns10$mu.fv, sigma = t.gam.ns10$sigma.fv, lower.tail = TRUE, log.p =
                                             FALSE))
        timevarying.severe=as.matrix(pGA(severe.rain, mu =
                                           t.gam.ns10$mu.fv, sigma = t.gam.ns10$sigma.fv, lower.tail = TRUE, log.p =
                                           FALSE))
        timevarying.extreme=as.matrix(pGA(extreme.rain, mu =
                                            t.gam.ns10$mu.fv, sigma = t.gam.ns10$sigma.fv, lower.tail = TRUE, log.p =
                                            FALSE))
        mu.trend=t.gam.ns10$mu.coefficients[2]; sigma.trend=0
        actual.prob.moderate=100*(probzero+(1-
                                              probzero)*(timevarying.moderate[f]))
        actual.prob.severe=100*(probzero+(1-
                                            probzero)*(timevarying.severe[f]))
        actual.prob.extreme=100*(probzero+(1-
                                             probzero)*(timevarying.extreme[f]))
      }
      if (suitable.dist.deltaAICc[month,local]=="model3"){
        f=length(t.gam.ns10$mu.fv)
        timevarying.moderate=as.matrix(pGA(moderate.rain, mu =
                                             t.gam.ns11$mu.fv, sigma = t.gam.ns11$sigma.fv, lower.tail = TRUE, log.p =
                                             FALSE))
        timevarying.severe=as.matrix(pGA(severe.rain, mu =
                                           t.gam.ns11$mu.fv, sigma = t.gam.ns11$sigma.fv, lower.tail = TRUE, log.p =
                                           FALSE))
        timevarying.extreme=as.matrix(pGA(extreme.rain, mu =
                                            t.gam.ns11$mu.fv, sigma = t.gam.ns11$sigma.fv, lower.tail = TRUE, log.p =
                                            FALSE))
        mu.trend=t.gam.ns11$mu.coefficients[2];
        sigma.trend=t.gam.ns11$sigma.coefficients[2]
        actual.prob.moderate=100*(probzero+(1-
                                              probzero)*(timevarying.moderate[f]))
        actual.prob.severe=100*(probzero+(1-
                                            probzero)*(timevarying.severe[f]))
        actual.prob.extreme=100*(probzero+(1-
                                             probzero)*(timevarying.extreme[f]))
      }
      if(p.moderate.drought>probzero){
        moderate[month,local]=actual.prob.moderate}
      else{moderate[month,local]="too many zeros"}
      if(p.severe.drought>probzero){
        severe[month,local]=actual.prob.severe}
      else{severe[month,local]="too many zeros"}
      if(p.extreme.drought>probzero){
        extreme[month,local]=actual.prob.extreme}
      else{extreme[month,local]="too many zeros"}
    }}
  colnames(moderate)=c("Moderate")
  colnames(severe)=c("Severe")
  colnames(extreme)=c("Extreme")
  if (ts==1){finalts1=format(cbind(moderate,severe,extreme),digits = 2)}
  if (ts==2){finalts2=format(cbind(moderate,severe,extreme),digits = 2)}
  if (ts==3){finalts3=format(cbind(moderate,severe,extreme),digits = 2)}
  if (ts==4){finalts4=format(cbind(moderate,severe,extreme),digits = 2)}
  if (ts==5){finalts5=format(cbind(moderate,severe,extreme),digits = 2)}
  if (ts==6){finalts6=format(cbind(moderate,severe,extreme),digits = 2)}
  if (ts==7){finalts7=format(cbind(moderate,severe,extreme),digits = 2)}
  if (ts==8){finalts8=format(cbind(moderate,severe,extreme),digits = 2)}
  if (ts==9){finalts9=format(cbind(moderate,severe,extreme),digits = 2)}
  if (ts==10){finalts10=format(cbind(moderate,severe,extreme),digits = 2)}
  if (ts==11){finalts11=format(cbind(moderate,severe,extreme),digits = 2)}
  if (ts==12){finalts12=format(cbind(moderate,severe,extreme),digits = 2)}
}
final=as.matrix(cbind(finalts1,finalts2,
                      finalts3,finalts4,finalts5,finalts6,
                      finalts7,finalts8,finalts9,
                      finalts10,finalts11,finalts12))
row.names(final)=c("Jan","Feb","Mar","Apr","May","Jun",
                   "Jul","Ago","Sep","Oct","Nov","Dec")
name.ts=c("","Time Scale 1","",
          "","Time Scale 2","",
          "","Time Scale 3","",
          "","Time Scale 4","",
          "","Time Scale 5","",
          "","Time Scale 6","",
          "","Time Scale 7","",
          "","Time Scale 8","",
          "","Time Scale 9","",
          "","Time Scale 10","",
          "","Time Scale 11","",
          "","Time Scale 12","")
final=rbind(name.ts,final)
write.csv(final,"Currently_Expected_Frequency.csv")