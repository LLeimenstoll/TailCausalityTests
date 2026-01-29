#Simulation Study for different configurations

setwd()
simulation_path<-"Simulation_Study/"
source("functions.R")
source("testing_strategy.R")

rep<-10000
n<-10000
k1<-floor(n^0.4)

cases<-0
deltas<-0
df1s<-0
df2s<-0
dfhs<-0

results<-data.frame()

for(r in 1:rep){
  df1<-sample(c(2, 3, 4), 1)
  df2<-sample(c(2, 3, 4), 1)
  dfh<-sample(c(2, 3, 4), 1)
  beta_h1<-runif(1,0.1,0.9)
  beta_h2<-runif(1,0.1,0.9)
  beta_1<-runif(1,0.1,0.9)
  e1<-rt(n,df1)
  e2<-rt(n,df2)
  conf<-rbinom(1,1,0.5)
  connect<-rbinom(1,1,0.5)
  case <- ifelse(conf == 0,
                 ifelse(connect == 0, 1, 3),
                 ifelse(connect == 0, 2, 4))
  H<-conf*rt(n,dfh)
  
  X1<-e1+beta_h1*H
  X2<-e2+connect*beta_1*X1+beta_h2*H
  
  cases[r]<-case
  df1s[r]<-df1
  df2s[r]<-df2
  dfhs[r]<-dfh
  # run testing strategy
  results[r,1:12]<-c(testing_strategy(X1,X2,k=k1),beta_1,beta_h1,beta_h2)
  # run testing strategy including confounder
  if(case %in% c(2,4)){
    res_conf<-testing_strategy(X1,X2,H=H,k=k1)
    results[r,13:17]<-c(res_conf)
  }
  print(r)
}

df<-cbind.data.frame(cases,results, df1s,df2s,dfhs)

colnames(df)<-c("case","dir","pretest","conftest","indtest","c12","c21",
                "hill1","hill2","p_val","beta1","beta_h1","beta_h2","dirH",
                "pretestH","conftestH","c12H","c21H","df1","df2","dfh")
head(df)

write.csv(df, file = paste0(simulation_path,"output/results_cases_k39.csv"),row.names =FALSE)

