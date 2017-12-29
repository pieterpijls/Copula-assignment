## ----global_options, include=FALSE---------------------------------------
knitr::opts_chunk$set(fig.pos = 'H')

## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----  include=FALSE-----------------------------------------------------
#Set seed and clear memory
rm(list=ls())
set.seed(0387948)

## ---- include=FALSE------------------------------------------------------
#Load packages
packages <- c("copula","mgcv", "fCopulae", "Ecdat", "fGarch", "MASS", "tidyverse", "ggthemes","plotly","ggplot2","xtable","knitr","statmod", "actuar", "fitdistrplus", "PerformanceAnalytics", "gridExtra", "grid", "scales", "quadprog","quantmod","plyr","reshape","kableExtra","xts","float")

suppressMessages(packages <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x)
    library(x, character.only = TRUE)
  }
}))

## ----  include=FALSE-----------------------------------------------------
# Load data
load("~/Documents/KU LEUVEN/MFAE/STATISTICS FOR FINANCE&INSURANCE/data_copula_387948.RData")

## ----  include=FALSE-----------------------------------------------------
# create boolean parameter
data$rc <- !(data$cens==0)

## ----  include=FALSE-----------------------------------------------------
#head(data)
summary(data)

## ----fig1, echo=FALSE, fig.cap="\\label{fig:fig1}Histogram of the loss and expense variable"----
# histograms
p1 <- ggplot(data, aes(x = loss)) +
  geom_histogram(binwidth = 50000) +
  theme_hc()

p2 <- ggplot(data, aes(x = expense)) +
  geom_histogram(binwidth = 10000 ) +
  theme_hc()

grid.arrange(p1,p2,ncol=2)

## ----fig2, echo=FALSE, fig.cap="\\label{fig:fig2}Scatterplot of loss and expense variable"----
# scatter
par(mfrow=c(1,1))
z <- ggplot(data, aes(x = log(loss), y=log(expense),colour=rc)) +
  geom_jitter() + theme_hc() + scale_colour_discrete(name="Censored?", breaks=c(FALSE,TRUE),labels=c("Not censored","Censored"))
z

## ---- include=F----------------------------------------------------------
#z <- ggplot(data, aes(x = log(loss), y=log(expense),group=)) +
  #geom_jitter(x=log(data$loss[data$rc]),color="red") + geom_jitter(x=log(data$loss[!data$rc]))+ theme_hc()

## ---- include=F----------------------------------------------------------
#Pearson correlation
cor(data$loss,data$expense)

## ---- include=FALSE------------------------------------------------------
#Summary stats
mean(data$loss)
mean(data$expense)

sd(data$loss)
sd(data$expense)

skewness(data$loss)
skewness(data$expense)

kurtosis(data$loss)
kurtosis(data$expense)

dt <- data.frame(matrix(ncol = 2, nrow = 4))
colnames(dt) <- c("Loss","Expense")
rownames(dt) <- c("Mean","Standard Deviation","Skewness","Kurtosis")
dt[1,] <- c(mean(data$loss),mean(data$expense))
dt[2,] <- c(sd(data$loss),sd(data$expense))
dt[3,] <- c(skewness(data$loss),skewness(data$expense))
dt[4,] <- c(kurtosis(data$loss),kurtosis(data$expense))

## ---- echo=FALSE---------------------------------------------------------
# Table
kable(dt, format = "latex", caption = "Summary Statistics", booktabs = T) %>% kable_styling(latex_options = c("striped", "hold_position"))

## ---- include=FALSE------------------------------------------------------
# Log likelihood functions
x  <- data$loss 
rc <- data$rc
loglik.pareto.loss <- function(par){ 
  sum(dpareto(x[!rc],par[1],par[2],log=T)) + sum(ppareto(x[rc],par[1],par[2],log.p=T,lower.tail=F)) 
} #Loss

loglik.pareto.expense <- fitdist(data$expense, "pareto") #Expenses

## ---- include=FALSE------------------------------------------------------
# Optimize
oo <- optim(c(1,13161),loglik.pareto.loss,control=list(fnscale=-1))
par.pareto.loss <- oo$par

## ----  include=FALSE-----------------------------------------------------
# Parameters
par.pareto.loss
loglik.pareto.expense$estimate

## ---- echo=FALSE---------------------------------------------------------
# Table with parameters
dt <- data.frame(matrix(ncol = 2, nrow = 2))
colnames(dt) <- c("Loss","Expense")
rownames(dt) <- c("alpha","beta")
dt[1,] <- c(par.pareto.loss[1],loglik.pareto.expense$estimate[1])
dt[2,] <- c(par.pareto.loss[2],loglik.pareto.expense$estimate[2])

kable(dt, format = "latex", caption = "Parameters Pareto Model", booktabs = T) %>% kable_styling(latex_options = c("striped", "hold_position"))

## ----fig3, echo=FALSE, fig.cap="\\label{fig:fig3}Comparison of Empirical CDF and Pareto CDF"----
par(mfrow=c(1,2))

#Loss
plot(ecdf(data$loss), do.points = FALSE, xlab = 'Loss', ylab = 'CDF', main = '',  xlim = c(0, max(data$loss)), lwd = 2)
# Fitted CDF
x<-data$loss 
curve(ppareto(x, shape=par.pareto.loss[1], scale=par.pareto.loss[2]), ylab="Probability", main="Pareto CDF",col=2, lwd=2, add=T)
legend('right', c('Empirical CDF', 'Fitted CDF'), col = c(1, 2), lwd = 2)

#Expense
plot(ecdf(data$expense), do.points = FALSE, xlab = 'Expense', ylab = 'CDF', main = '', xlim = c(0, max(data$ expense)), lwd = 2)
# Fitted CDF
curve(ppareto(x, shape=loglik.pareto.expense$estimate[1], scale=loglik.pareto.expense$estimate[2]), ylab="Probability", main="Pareto CDF",col=2, lwd=2, add=T)
legend('right', c('Empirical CDF', 'Fitted CDF'), col = c(1, 2), lwd = 2)

## ---- echo=F-------------------------------------------------------------
# AIC and BIC
dt <- data.frame(matrix(ncol = 2, nrow = 2))
colnames(dt) <- c("Loss","Expense")
rownames(dt) <- c("AIC","BIC")

AIC.pareto.loss <- 2*2-2*oo$value
BIC.pareto.loss <- log(nrow(data))*2-2*oo$value

dt[1,] <- c(AIC.pareto.loss,loglik.pareto.expense$aic)
dt[2,] <- c(BIC.pareto.loss,loglik.pareto.expense$bic)

kable(dt, format = "latex", caption = "Information Criteria", booktabs = T) %>% kable_styling(latex_options = c("striped", "hold_position"))

## ----fig4, echo=FALSE, fig.cap="\\label{fig:fig4}Histogram of the the pseudo-observations"----
# Pseudo-observations
u1 <- ppareto(data$loss, shape=par.pareto.loss[1], scale = par.pareto.loss[2])
u2 <- ppareto(data$expense, shape=loglik.pareto.expense$estimate[1], scale = loglik.pareto.expense$estimate[2])

#uniform-transformed 
U.hat = cbind(u1,u2)

#histograms
par(mfrow=c(1,2), cex.axis=1.2, cex.lab=1.2, cex.main=1.2)
hist(u1, main="Loss", xlab=expression(hat(U)[1]), freq = FALSE)
hist(u2, main="Expense", xlab=expression(hat(U)[2]), freq = FALSE)

## ----fig6, echo=FALSE, fig.cap="\\label{fig:fig6}Bivariate scatterplot of the uniform-transformed data"----
# Scatter 
U.hatgraph = as.data.frame(U.hat)
par(mfrow=c(1,1))
z <- ggplot(U.hatgraph, aes(x = U.hatgraph$u1 ,y=U.hatgraph$u2)) +
  geom_jitter() + theme_hc()
z

## ---- include=F----------------------------------------------------------
loss    <- data[,1] 
expense <- data[,2]

# create unfirom from zero to one
z <- seq(0,1,length.out =100)

# Kendall Distribution Function for Archimedean Copulas
Kcop <- function(x) {
  Kn(x,cbind(loss,expense),method="GR")
}

# calculate Kendells Tau
Kendall   <- cor(loss, expense, method="kendall") # Kendell tau estimate
Kendall

# Calculate theta of Gumbel copula
# tau = 1-(1/theta) <=> theta = 1/(1-tau)
Gumbel.theta  <- (-1)/(Kendall-1)
Gumbel.theta

# Calculate theta of Clayton copula
# tau = theta/(theta+2) <=> theta = 2*tau/(1-tau) 
Clayton.theta <- (2*Kendall)/(1-Kendall)
Clayton.theta

## ---- include=F----------------------------------------------------------
# Enpirical cdf Clayton
Clayton.K <- function(x){
  x-(x^(-Clayton.theta)-1)/(-Clayton.theta*x^(-Clayton.theta-1))
    }

# Empirical cdf Gumbel
Gumbel.K <- function(x){
   x - ((-log(x))^Gumbel.theta)/(-Gumbel.theta * (-log(x))^(Gumbel.theta-1)/x)
    }


## ----fig7, echo=FALSE, fig.height=4,fig.width=10 ,fig.cap="\\label{fig:fig7}Archimidean Copula Identification: Non-parametric estimation (black), Non-parametric Clayton (red), Non-parametric Gumbel (green)"----
plot(z, Kcop(z),col=1 , "l", ylim=c(0,1))
curve(Clayton.K(x),add=T,col=2)
curve(Gumbel.K(x),add=T,col=3)

## ---- include=F----------------------------------------------------------
# empirical cumulative distribution function
U.np=cbind((rank(data$loss)-0.5)/(nrow(data)),(rank(data$expense)-0.5)/(nrow(data)))
colnames(U.np)=c("u1","u2")
head(U.np)

## ---- include=F----------------------------------------------------------
# Clayton
Ccl = fitCopula(copula=claytonCopula(1, dim=2), data=U.np, method="ml")
Ccl@estimate

## ---- include=F----------------------------------------------------------
# Gumbel 
#Parametric copulas were fit to uniform-transformed flows (Gumbel gives error since data shows negative dependence)
Cgu = fitCopula(copula=gumbelCopula(), data=U.np, method="ml")
Cgu@estimate

## ----fig8, echo=FALSE,fig.height=4,fig.width=10 , fig.cap="\\label{fig:fig8}Empricial Copula contours (black) and Non-parametric Clayton and Gumbel (red)"----
#empirical copula
n = nrow(data)
Udex = (1:n)/(n+1)
Cn = C.n(u = cbind(rep(Udex, n), rep(Udex, each=n)) , U.np, offset=0, method="C")
EmpCop = expression(contour(Udex,Udex,matrix(Cn,n,n), col=2, add=T))

par(mfrow=c(1,2))

# Non parametric Clayton
contour(claytonCopula(param=Ccl@estimate[1], dim = 2), 
        pCopula, main = "Non parametric Clayton",
        xlab = expression(hat(U)[1]), ylab = expression(hat(U)[2]))
eval(EmpCop)

# Non parametric Gumbel
contour(gumbelCopula(param=Cgu@estimate[1], dim = 2), 
        pCopula, main = "Non-Parametric Gumbel",
        xlab = expression(hat(U)[1]), ylab = expression(hat(U)[2]))
eval(EmpCop)

## ---- include=F----------------------------------------------------------
#Parametric copula were U.hat is modeled with Pareto parameters taking censoring into account
Cgu.parametric = fitCopula(copula=gumbelCopula(), data=U.hat, method="ml")
Cgu.parametric@estimate

## ----fig9, echo=FALSE,fig.height=4,fig.width=10 , fig.cap="\\label{fig:fig9}Empricial Copula contours (black) and non-parametric Gumbel versus Parametric Gumbel (red)"----
#empirical copula
par(mfrow=c(1,2))

# Non parametric
contour(gumbelCopula(param=Cgu@estimate[1], dim = 2), 
        pCopula, main = "Non-Parametric Gumbel",
        xlab = expression(hat(U)[1]), ylab = expression(hat(U)[2]))
eval(EmpCop)

# Parametric
contour(gumbelCopula(param=Cgu.parametric@estimate[1], dim = 2), 
        pCopula, main = "Parametric Gumbel",
        xlab = expression(hat(U)[1]), ylab = expression(hat(U)[2]))
eval(EmpCop)

## ---- echo=FALSE---------------------------------------------------------
# AIC
dt <- data.frame(matrix(ncol = 1, nrow = 3))
colnames(dt) <- c("AIC")
rownames(dt) <- c("Non-Parametric Clayton","Non-Parametric Gumbel","Parametric Gumbel")

AIC.pareto.loss <- 2*2-2*oo$value
BIC.pareto.loss <- log(nrow(data))*2-2*oo$value

dt[1,] <- AIC(Ccl) # Non-parametric Clayton
dt[2,] <- AIC(Cgu) # Non-parametric Gumbel
dt[3,] <- AIC(Cgu.parametric) # Parametric Clayton

kable(dt, format = "latex", caption = "AIC Models", booktabs = T) %>% kable_styling(latex_options = c("striped", "hold_position"))

## ---- include=F----------------------------------------------------------
#simulating observations with the given copula
U <- rCopula(1000, copula = gumbelCopula(param=Cgu.parametric@estimate , dim = 2))
plot(U[,1:2], xlab = expression(U[1]), ylab = expression(U[2]))

## ---- include=F----------------------------------------------------------
#In the following step take the inverse of the Pareto cdf. 
loss <- pinvpareto(U[,1], shape=par.pareto.loss[1],scale=par.pareto.loss[2])
expense <- pinvpareto(U[,2], shape=loglik.pareto.expense$estimate[1],scale=loglik.pareto.expense$estimate[2])

randomsample <- cbind(loss,expense)

plot(randomsample[,1:2], xlab = "Loss", ylab ="Expense")

## ---- include=F----------------------------------------------------------
# scatter (no log transformation)
par(mfrow=c(1,1))
z <- ggplot(data, aes(x = loss, y=expense,colour=rc)) +
  geom_jitter() + theme_hc() + scale_colour_discrete(name="Censored?", breaks=c(FALSE,TRUE),labels=c("Not censored","Censored"))
z

