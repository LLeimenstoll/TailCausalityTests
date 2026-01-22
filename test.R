#Test code

X<-rnorm(100)
Y<-X+rnorm(100)
l1<-lm(Y~X)

summary(l1)
