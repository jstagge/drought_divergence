# *------------------------------------------------------------------
# | PROGRAM NAME: continental_diverge
# | FILE NAME: 01_continental_diverge.R
# | DATE: 
# | CREATED BY:  Jim Stagge         
# *----------------------------------------------------------------
# | PURPOSE:  This script calculates the trend of drought occurence
# |				using splines at the continental scale for Europe.
# |				Reproduces Figure 1 from "Observed Drought Indices ..."
# |
# *------------------------------------------------------------------
# | COMMENTS:               
# |
# |  1:  
# |  2:  
# |  3: 
# |*------------------------------------------------------------------
# | DATA USED:               
# | 
# |
# |*------------------------------------------------------------------
# | CONTENTS:               
# |
# |  PART 1:  
# |  PART 2: 
# |  PART 3: 
# *-----------------------------------------------------------------
# | UPDATES:               
# |
# |
# *------------------------------------------------------------------
 
### Clear any existing data or functions.
rm(list=ls())

###########################################################################
## Set the Paths
###########################################################################
### Path for Data and Output	
data_path <- "./data"
output_path <- "./output"
function_path <- "./functions"

### Set output location
output_name <- "continental_diverge"

write_output_path <- file.path(output_path,output_name)
write_fig_path <- file.path(output_path, "figures")

###########################################################################
###  Load functions
###########################################################################
### Load these functions for this unique project
require(mgcv)
require(ggplot2)

### Load project specific functions
file.sources = list.files(function_path, pattern="*.R", recursive=TRUE)
sapply(file.path(function_path, file.sources),source)

###########################################################################
## Set Initial Values
###########################################################################


################################################
###  Read in Data
#################################################

#### Read in Monthly SPI
WFD_SPI.df <- read.csv(file.path(data_path, "WFD/WFD_AreaBelow20_SPI_monthly_clusters6.csv"))
WFD_SPI.df$Data <- "WFD"
WFDEI_SPI.df <- read.csv(file.path(data_path, "WFDEI/WFDEI_AreaBelow20_SPI_monthly_clusters6.csv"))
WFDEI_SPI.df$Data <- "WFDEI"
SPI.df <- rbind(WFD_SPI.df, WFDEI_SPI.df)
SPI.df$Data <- factor(SPI.df$Data)
SPI.df$Date <- as.Date(SPI.df$Date)
SPI.df$YearFrac <- as.numeric(SPI.df$Date)
SPI.df$Year <- as.numeric(substr(SPI.df$Date,1,4))
SPI.df$Europe <- round(SPI.df$Europe, digits = 5)

#### Read in SPEI
WFD_SPEI.df <- read.csv(file.path(data_path, "WFD/WFD_AreaBelow20_SPEI_monthly_clusters6.csv"))
WFD_SPEI.df$Data <- "WFD"
WFDEI_SPEI.df <- read.csv(file.path(data_path, "WFDEI/WFDEI_AreaBelow20_SPEI_monthly_clusters6.csv"))
WFDEI_SPEI.df$Data <- "WFDEI"
SPEI.df <- rbind(WFD_SPEI.df, WFDEI_SPEI.df)
SPEI.df$Data <- factor(SPEI.df$Data)
SPEI.df$Date <- as.Date(SPEI.df$Date)
SPEI.df$YearFrac <- as.numeric(SPEI.df$Date)
SPEI.df$Year <- as.numeric(substr(SPEI.df$Date,1,4))
SPEI.df$Europe <- round(SPEI.df$Europe, digits = 5)

### Calculate Diff
Diff.df <- SPI.df
Diff.df[,2:8] <- SPEI.df[,2:8] - SPI.df[,2:8]

n.obs <- dim(SPI.df)[1]
################################################
###  Fit SPI Model
#################################################
## Fit a smoother for Year to the data
# Optimizer with upper/lower bounds
#ctrl <- list(niterEM = 0, msVerbose = TRUE, optimMethod="L-BFGS-B")
ctrl <- list(niterEM = 0, msVerbose = TRUE, optimMethod="BFGS")

### Fit plain model
m1.SPI <- gamm(Europe ~ s(YearFrac, k = 10) + Data, data = SPI.df,   control = ctrl, method="REML", family=binomial, weights=rep(1e5, n.obs))
summary(m1.SPI$gam)
plot(m1.SPI$gam, residuals=TRUE)

## Account for potential annual cycle
m2.SPI <- gamm(Europe ~ s(YearFrac, k = 10) + s(Month, bs = "cc", k = 12) + Data, data = SPI.df,   control = ctrl, method="REML", family=binomial, weights=rep(1e5, n.obs))
summary(m2.SPI$gam)
plot(m2.SPI$gam, residuals=TRUE)
 
## Compare models
#anova(m1.SPI$lme, m2.SPI$lme)	
AIC(m1.SPI$lme, m2.SPI$lme)	 
BIC(m1.SPI$lme, m2.SPI$lme)	 
  
  
## look at autocorrelation in residuals:
acf(resid(m2.SPI$lme, type = "normalized"))
## ...wait... look at the plot, only then do...
pacf(resid(m2.SPI$lme, type = "normalized"))
## seems like like sharp cut-off in ACF and PACF - AR terms probably best
 
## Fit the AR1
m3.SPI.monthly <- try(gamm(Europe ~ s(YearFrac, k = 10) + s(Month, bs = "cc", k = 12) + Data, data = SPI.df, correlation = corARMA(form = ~ 1|Year, p = 1),   control = ctrl, method="REML", family=binomial(link="logit"), weights=rep(1e5, n.obs)))
## ...and fit the AR2
m4.SPI.monthly <- try(gamm(Europe ~ s(YearFrac, k = 10) + s(Month, bs = "cc", k = 12) + Data, data = SPI.df, correlation = corARMA(form = ~ 1|Year, p = 2),   control = ctrl, method="REML", family=binomial(link="logit"), weights=rep(1e5, n.obs)))
## ...and fit the AR3
m5.SPI.monthly <- try(gamm(Europe ~ s(YearFrac, k = 10) + s(Month, bs = "cc", k = 12) + Data, data = SPI.df, correlation = corARMA(form = ~ 1|Year, p = 3),   control = ctrl, method="REML", family=binomial(link="logit"), weights=rep(1e5, n.obs)))

#anova(m1.SPI$lme, m2.SPI$lme)	
#### Adding monthly doesn't do anything   
#anova(m1.SPI$lme, m3.SPI.monthly$lme, m4.SPI.monthly$lme, m5.SPI.monthly$lme)	
AIC(m1.SPI$lme, m3.SPI.monthly$lme, m4.SPI.monthly$lme, m5.SPI.monthly$lme)	
BIC(m1.SPI$lme, m3.SPI.monthly$lme, m4.SPI.monthly$lme, m5.SPI.monthly$lme)	


#######################
###  It looks like monthly pattern is not useful
############################
## ...so fit the AR1
m3.SPI <- try(gamm(Europe ~ s(YearFrac, k = 10) + Data, data = SPI.df, correlation = corARMA(form = ~ 1|Year, p = 1),   control = ctrl, method="REML", family=binomial, weights=rep(1e5, n.obs)))
summary(m3.SPI$gam)
#### No noticeable difference between WFD and WFDEI, drop the Data term
## ...and fit the AR2
m4.SPI <- try(gamm(Europe ~ s(YearFrac, k = 3), data = SPI.df, correlation = corARMA(form = ~ 1|Year, p = 2),   control = ctrl, method="REML", family=binomial, weights=rep(1e5, n.obs)))
## ...and fit the AR3
m5.SPI <- try(gamm(Europe ~ s(YearFrac, k = 10), data = SPI.df, correlation = corARMA(form = ~ 1|Year, p = 3),   control = ctrl, method="REML", family=binomial, weights=rep(1e5, n.obs)))	

### Compare models with AR terms   
#anova(m1$lme, m3.SPI$lme, m3.SPI.monthly$lme)
#anova(m1$lme, m3.SPI$lme, m4.SPI$lme, m5.SPI$lme)	
AIC(m1.SPI$lme, m3.SPI$lme, m3.SPI.monthly$lme)
AIC(m1.SPI$lme, m3.SPI$lme, m4.SPI$lme, m5.SPI$lme)

### AR2 without monthly pattern is the best
m.final <- m4.SPI

summary(m.final$gam)
plot(m.final$gam, residuals=TRUE)

## look at autocorrelation in residuals:
acf(resid(m.final$lme, type = "normalized"))
pacf(resid(m.final$lme, type = "normalized"))
## No more temporal autocorrelation

#################################
###  Calculate significance
#########################################
#### No significance test
### Year has a p-vale of 0.0205  < 0.025
### Significant trend
summary(m.final$gam)

just.WFD <- subset(SPI.df, Data=="WFD")
mfinal.d <- Deriv(m.final$gam, newdata = just.WFD)
plot(mfinal.d, sizer = TRUE, alpha = 0.05)

plot(Europe ~ Date, data = just.WFD, type = "p")
p1 <- predict(m.final$gam, newdata = just.WFD, type="response")
lines(p1 ~ just.WFD$Date, col="black")

CI <- confint(mfinal.d, alpha = 0.05)
S <- signifD(p1, mfinal.d$Year$deriv, CI$Year$upper, CI$Year$lower, eval = 0)
lines(S$incr ~ Date, data = just.WFD, lwd = 3, col = "red")
lines(S$decr ~ Date, data = just.WFD, lwd = 3, col = "red")
Pred.WFD <- data.frame(Date=just.WFD$Date, Pred=p1, Inc=S$incr, Dec=S$decr, Data="WFD")

just.WFDEI <- subset(SPI.df, Data=="WFDEI")
mfinal.d <- Deriv(m.final$gam, newdata = just.WFDEI)
plot(mfinal.d, sizer = TRUE, alpha = 0.05)

plot(Europe ~ Date, data = just.WFDEI, type = "p")
p1 <- predict(m.final$gam, newdata = just.WFDEI, type="response")
lines(p1 ~ just.WFDEI$Date, col="black")

CI <- confint(mfinal.d, alpha = 0.05)
S <- signifD(p1, mfinal.d$Year$deriv, CI$Year$upper, CI$Year$lower, eval = 0)
lines(S$incr ~ Date, data = just.WFDEI, lwd = 3, col = "red")
lines(S$decr ~ Date, data = just.WFDEI, lwd = 3, col = "red")
Pred.WFDEI <- data.frame(Date=just.WFDEI$Date, Pred=p1, Inc=S$incr, Dec=S$decr, Data="WFDEI")
Pred.df <- rbind(Pred.WFD, Pred.WFDEI)



Pred.comb <- SPI.df[!duplicated(SPI.df$Date),]
Pred.comb$Data <- "WFD"
mfinal.d <- Deriv(m.final$gam, newdata = Pred.comb)
plot(mfinal.d, sizer = TRUE, alpha = 0.05)

plot(Europe ~ Date, data = Pred.comb, type = "p")
p1 <- predict(m.final$gam, newdata = Pred.comb, type="response")
lines(p1 ~ Pred.comb$Date, col="black")

CI <- confint(mfinal.d, alpha = 0.05)
S <- signifD(p1, mfinal.d$Year$deriv, CI$Year$upper, CI$Year$lower, eval = 0)
lines(S$incr ~ Date, data = Pred.comb, lwd = 3, col = "red")
lines(S$decr ~ Date, data = Pred.comb, lwd = 3, col = "red")
Pred.comb <- data.frame(Date=Pred.comb$Date, Pred=p1, Inc=S$incr, Dec=S$decr, Data="WFD")



######

#p <- ggplot(SPI.df, aes(x=Date, group=Data))
p <- ggplot(SPI.df, aes(x=Date, group=Data))
p <- p + geom_hline(yintercept=0.2, linetype="longdash", colour="grey30")
p <- p + geom_line(aes(y=Europe), colour="grey70", size=0.3)
#p <- p + geom_line(aes(y=Europe, colour=Data), size=0.3,alpha=0.7)
p <- p + geom_line(data=Pred.df, aes(y=Pred), colour="black", size = 0.75, colour="grey20")
p <- p + geom_line(data=Pred.df, aes(y=Inc), colour="blue", size=1)
p <- p + geom_line(data=Pred.df, aes(y=Dec), colour="red", size=1)
#p <- p + geom_segment(aes(x = as.Date("1971-01-01"), y = 0.025, xend = as.Date("2000-12-31"), yend = 0.025), arrow=arrow(length = unit(0.1, "inches"), angle = 15, ends = "both", type = "closed"))
p <- p + geom_segment(aes(x = as.Date("1971-01-01"), y = 0.015, xend = as.Date("2000-12-31"), yend = 0.015), size=2, alpha=0.8, colour="grey80")
p <- p + scale_y_continuous(name=expression(paste("SPI6 Drought Area (", A[SPI6], ", %)", sep="")), labels = percent)
p <- p + theme_classic(9) + theme(
axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
p <- p +  theme(axis.title.y=element_text(margin=margin(0,12,0,0)))

p <- p + theme(legend.position=c(0.25,0.85), legend.direction="horizontal")
p <- p + scale_colour_brewer(type="qual", palette="Set1")
p <- p + coord_cartesian(ylim=c(0,0.5))
p

ggsave(file.path(write_fig_path, "SPI6_Version1.png"), p, width=7, height=4, dpi=400)
ggsave(file.path(write_fig_path, "SPI6_Version1.pdf"), p, width=7, height=4)




################################################
###  Fit SPEI Model
#################################################
m1.SPEI <- gamm(Europe ~ s(YearFrac, k = 10) + Data, data = SPEI.df,   control = ctrl, method="REML", family=binomial, weights=rep(1e5, n.obs))
summary(m1.SPEI$gam)
plot(m1.SPEI$gam, residuals=TRUE)

## Account for potential annual cycle
m2.SPEI <- gamm(Europe ~ s(YearFrac, k = 10) + s(Month, bs = "cc", k = 12) + Data, data = SPEI.df,   control = ctrl, method="REML", family=binomial, weights=rep(1e5, n.obs))
summary(m2.SPEI$gam)
plot(m2.SPEI$gam, residuals=TRUE)
### Weirdly month term is slightly significant

 ## Compare models
anova(m1.SPEI$lme, m2.SPEI$lme)	
AIC(m1.SPEI$lme, m2.SPEI$lme)	
BIC(m1.SPEI$lme, m2.SPEI$lme)	
 
## look at autocorrelation in residuals:
acf(resid(m2.SPEI$lme, type = "normalized"))
## ...wait... look at the plot, only then do...
pacf(resid(m2.SPEI$lme, type = "normalized"))
## seems like like sharp cut-off in ACF and PACF - AR terms probably best
 
 
## ...so fit the AR1
m3.SPEI.monthly <- try(gamm(Europe ~ s(YearFrac, k = 10) + s(Month, bs = "cc", k = 12) + Data, data = SPEI.df, correlation = corARMA(form = ~ 1|Year, p = 1),   control = ctrl, method="REML", family=binomial, weights=rep(1e5, n.obs)))
## ...and fit the AR2
m4.SPEI.monthly <- try(gamm(Europe ~ s(YearFrac, k = 10) + s(Month, bs = "cc", k = 12) + Data, data = SPEI.df, correlation = corARMA(form = ~ 1|Year, p = 2),   control = ctrl, method="REML", family=binomial, weights=rep(1e5, n.obs)))
## ...and fit the AR3
m5.SPEI.monthly <- try(gamm(Europe ~ s(YearFrac, k = 10) + s(Month, bs = "cc", k = 12) + Data, data = SPEI.df, correlation = corARMA(form = ~ 1|Year, p = 3),   control = ctrl, method="REML", family=binomial, weights=rep(1e5, n.obs)))	

anova(m1.SPEI.monthly$lme, m2.SPEI.monthly$lme)	
#### Adding monthly doesn't do anything   
anova(m1.SPEI.monthly$lme, m3.SPEI.monthly$lme, m4.SPEI.monthly$lme, m5.SPEI.monthly$lme)	
AIC(m1.SPEI.monthly$lme, m2.SPEI.monthly$lme)	
AIC(m1.SPEI.monthly$lme, m2.SPEI.monthly$lme,m3.SPEI.monthly$lme, m4.SPEI.monthly$lme, m5.SPEI.monthly$lme)	
BIC(m1.SPEI.monthly$lme, m2.SPEI.monthly$lme)	
BIC(m1.SPEI.monthly$lme, m2.SPEI.monthly$lme,m3.SPEI.monthly$lme, m4.SPEI.monthly$lme, m5.SPEI.monthly$lme)	

1/(1+exp(-(m4.SPEI.monthly$gam$coefficients[[1]] + m4.SPEI.monthly$gam$coefficients[[2]])))
1/(1+exp(-(m4.SPEI.monthly$gam$coefficients[[1]])))

#### M4 seems best again, try without monthly term
# ...so fit the AR1
m3.SPEI <- try(gamm(Europe ~ s(YearFrac, k = 10) + Data, data = SPEI.df, correlation = corARMA(form = ~ 1|Year, p = 1),   control = ctrl, method="REML", family=binomial, weights=rep(1e5, n.obs)))
## ...and fit the AR2
m4.SPEI <- try(gamm(Europe ~ s(YearFrac, k = 10) + Data, data = SPEI.df, correlation = corARMA(form = ~ 1|Year, p = 2),   control = ctrl, method="REML", family=binomial, weights=rep(1e5, n.obs)))
## ...and fit the AR3
m5.SPEI <- try(gamm(Europe ~ s(YearFrac, k = 10) + Data, data = SPEI.df, correlation = corARMA(form = ~ 1|Year, p = 3),   control = ctrl, method="REML", family=binomial, weights=rep(1e5, n.obs)))	

#anova(m1$lme, m3$lme, m4$lme, m5$lme)	
 AIC(m1.SPEI$lme, m2.SPEI$lme)	
AIC(m1.SPEI$lme, m2.SPEI$lme,m3.SPEI.monthly$lme, m4.SPEI.monthly$lme, m5.SPEI.monthly$lme,m3.SPEI$lme, m4.SPEI$lme, m5.SPEI$lme)	
BIC(m1.SPEI$lme, m2.SPEI$lme,m3.SPEI.monthly$lme, m4.SPEI.monthly$lme, m5.SPEI.monthly$lme,m3.SPEI$lme, m4.SPEI$lme, m5.SPEI$lme)	

 
### Looks like AR2 is the best
m.final <- m4.SPEI

summary(m.final$gam)
plot(m.final$gam, residuals=TRUE)
vis.gam(m.final$gam)
 
## look at autocorrelation in residuals:
acf(resid(m.final$lme, type = "normalized"))
pacf(resid(m.final$lme, type = "normalized"))
## No more temporal autocorrelation

layout(matrix(1:2, ncol = 2))
plot(m.final$gam, scale = 0)
layout(1)

dev.new()
layout(matrix(1:2, ncol = 2))
res <- resid(m.final$lme, type = "normalized")
acf(res, lag.max = 36, main = "ACF - AR(2) errors")
pacf(res, lag.max = 36, main = "pACF- AR(2) errors")
layout(1)


#################################
###  Calculate significance
#########################################
#### No significance test
### Year has a p-vale of 0.0205  < 0.025
### Significant trend


just.WFD <- subset(SPEI.df, Data=="WFD")
mfinal.d <- Deriv(m.final$gam, newdata = just.WFD)
plot(mfinal.d, sizer = TRUE, alpha = 0.05)

plot(Europe ~ Date, data = just.WFD, type = "p")
p1 <- predict(m.final$gam, newdata = just.WFD, type="response")
lines(p1 ~ just.WFD$Date, col="black")

CI <- confint(mfinal.d, alpha = 0.05)
S <- signifD(p1, mfinal.d$Year$deriv, CI$Year$upper, CI$Year$lower, eval = 0)
lines(S$incr ~ Date, data = just.WFD, lwd = 3, col = "red")
lines(S$decr ~ Date, data = just.WFD, lwd = 3, col = "red")
Pred.WFD <- data.frame(Date=just.WFD$Date, Pred=p1, Inc=S$incr, Dec=S$decr, Data="WFD")

just.WFDEI <- subset(SPEI.df, Data=="WFDEI")
mfinal.d <- Deriv(m.final$gam, newdata = just.WFDEI)
plot(mfinal.d, sizer = TRUE, alpha = 0.05)

plot(Europe ~ Date, data = just.WFDEI, type = "p")
p1 <- predict(m.final$gam, newdata = just.WFDEI, type="response")
lines(p1 ~ just.WFDEI$Date, col="black")

CI <- confint(mfinal.d, alpha = 0.05)
S <- signifD(p1, mfinal.d$Year$deriv, CI$Year$upper, CI$Year$lower, eval = 0)
lines(S$incr ~ Date, data = just.WFDEI, lwd = 3, col = "red")
lines(S$decr ~ Date, data = just.WFDEI, lwd = 3, col = "blue")
Pred.WFDEI <- data.frame(Date=just.WFDEI$Date, Pred=p1, Inc=S$incr, Dec=S$decr, Data="WFDEI")
Pred.df <- rbind(Pred.WFD, Pred.WFDEI)


Pred.comb <- SPEI.df[!duplicated(SPEI.df$Date),]
Pred.comb$Data <- "WFD"
mfinal.d <- Deriv(m.final$gam, newdata = Pred.comb)
plot(mfinal.d, sizer = TRUE, alpha = 0.05)

plot(Europe ~ Date, data = Pred.comb, type = "p")
p1 <- predict(m.final$gam, newdata = Pred.comb, type="response")
lines(p1 ~ Pred.comb$Date, col="black")

CI <- confint(mfinal.d, alpha = 0.05)
S <- signifD(p1, mfinal.d$Year$deriv, CI$Year$upper, CI$Year$lower, eval = 0)
lines(S$incr ~ Date, data = Pred.comb, lwd = 3, col = "red")
lines(S$decr ~ Date, data = Pred.comb, lwd = 3, col = "blue")
Pred.comb <- data.frame(Date=Pred.comb$Date, Pred=p1, Inc=S$incr, Dec=S$decr, Data="WFD")





#############

p <- ggplot(SPEI.df, aes(x=Date, group=Data))
p <- p + geom_hline(yintercept=0.2, linetype="longdash", colour="grey30")
p <- p + geom_line(aes(y=Europe), colour="grey70", size=0.3)
#p <- p + geom_line(aes(y=Europe, colour=Data), size=0.3,alpha=0.7)
p <- p + geom_line(data=Pred.df, aes(y=Pred), colour="black", size = 0.75, colour="grey20")
p <- p + geom_line(data=Pred.df, aes(y=Inc), colour="red", size=1)
p <- p + geom_line(data=Pred.df, aes(y=Dec), colour="blue", size=1)
#p <- p + geom_segment(aes(x = as.Date("1971-01-01"), y = 0.025, xend = as.Date("2000-12-31"), yend = 0.025), arrow=arrow(length = unit(0.1, "inches"), angle = 15, ends = "both", type = "closed"))
p <- p + geom_segment(aes(x = as.Date("1971-01-01"), y = 0.015, xend = as.Date("2000-12-31"), yend = 0.015), size=2, alpha=0.8, colour="grey80")
p <- p + scale_y_continuous(name=expression(paste("SPEI6 Drought Area (", A[SPEI6], ", %)", sep="")), labels = percent)
p <- p + theme_classic(9) + theme(
axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
p <- p +  theme(axis.title.y=element_text(margin=margin(0,12,0,0)))

p <- p + theme(legend.position=c(0.25,0.85), legend.direction="horizontal")
p <- p + scale_colour_brewer(type="qual", palette="Set1")
p <- p + coord_cartesian(ylim=c(0,0.5))
p

ggsave(file.path(write_fig_path, "SPEI6_Version1.png"), p, width=7, height=4, dpi=400)
ggsave(file.path(write_fig_path, "SPEI6_Version1.pdf"), p, width=7, height=4)




#################################
## Redo WFDEI with adjustment for difference in data 
#################################
just.WFDEI <- subset(SPEI.df, Data=="WFDEI")
just.WFDEI$Data <- "WFD"
just.WFDEI$Data <- factor(just.WFDEI$Data)
p1 <- predict(m.final$gam, newdata = just.WFDEI, type="response")

WFDEI.adj <- just.WFDEI$Europe
## Put into logit space
WFDEI.adj <- log(WFDEI.adj/(1-WFDEI.adj))
WFDEI.adj <- WFDEI.adj - m.final$gam$coefficients[[2]]
WFDEI.adj <- 1/(1+exp(-WFDEI.adj))
just.WFDEI$Europe <- WFDEI.adj
just.WFDEI$Data <- "WFDEI"
SPEI.df.adj <- rbind(just.WFD, just.WFDEI)



p <- ggplot(SPEI.df.adj, aes(x=Date, group=Data))
p <- p + geom_hline(yintercept=0.2, linetype="longdash", colour="grey30")
p <- p + geom_line(aes(y=Europe), colour="grey70", size=0.3)
#p <- p + geom_line(aes(y=Europe, colour=Data), size=0.3,alpha=0.7)
p <- p + geom_line(data=Pred.comb, aes(y=Pred), colour="black", size = 0.75, colour="grey20")
p <- p + geom_line(data=Pred.comb, aes(y=Inc), colour="red", size=1)
p <- p + geom_line(data=Pred.comb, aes(y=Dec), colour="blue", size=1)
#p <- p + geom_segment(aes(x = as.Date("1971-01-01"), y = 0.025, xend = as.Date("2000-12-31"), yend = 0.025), arrow=arrow(length = unit(0.1, "inches"), angle = 15, ends = "both", type = "closed"))
p <- p + geom_segment(aes(x = as.Date("1971-01-01"), y = 0.015, xend = as.Date("2000-12-31"), yend = 0.015), size=2, alpha=0.8, colour="grey80")
p <- p + scale_y_continuous(name=expression(paste("SPEI6 Drought Area (", A[SPEI6], ", %)", sep="")), labels = percent)
p <- p + theme_classic(9) + theme(
axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
p <- p +  theme(axis.title.y=element_text(margin=margin(0,12,0,0)))

p <- p + theme(legend.position=c(0.25,0.85), legend.direction="horizontal")
p <- p + scale_colour_brewer(type="qual", palette="Set1")
p <- p + coord_cartesian(ylim=c(0,0.5))
p

ggsave(file.path(write_fig_path, "SPEI6_Version1_adj.png"), p, width=7, height=4, dpi=400)
ggsave(file.path(write_fig_path, "SPEI6_Version1_adj.pdf"), p, width=7, height=4)




################################################
###  Fit Diff Model
#################################################
m1.Diff <- gamm(Europe ~ s(YearFrac, k = 10) + Data, data = Diff.df,   control = ctrl, method="REML")
summary(m1.Diff$gam)
plot(m1.Diff$gam, residuals=TRUE)

## Account for potential annual cycle
m2.Diff <- gamm(Europe ~ s(YearFrac, k = 10) + s(Month, bs = "cc", k = 12) + Data, data = Diff.df,   control = ctrl, method="REML")
summary(m2.Diff$gam)
plot(m2.Diff$gam, residuals=TRUE)
### Weirdly month term is slightly significant

 ## Compare models
anova(m1.Diff$lme, m2.Diff$lme)	
AIC(m1.Diff$lme, m2.Diff$lme)	
 
## look at autocorrelation in residuals:
acf(resid(m2.Diff$lme, type = "normalized"))
## ...wait... look at the plot, only then do...
pacf(resid(m2.Diff$lme, type = "normalized"))
## seems like like sharp cut-off in ACF and PACF - AR terms probably best
 
 
## ...so fit the AR1
m3.Diff.month <- try(gamm(Europe ~ s(YearFrac, k = 10) + s(Month, bs = "cc", k = 12) + Data, data = Diff.df, correlation = corARMA(form = ~ 1|Year, p = 1),   control = ctrl, method="REML"))
## ...and fit the AR2
m4.Diff.month <- try(gamm(Europe ~ s(YearFrac, k = 10) + s(Month, bs = "cc", k = 12) + Data, data = Diff.df, correlation = corARMA(form = ~ 1|Year, p = 2),   control = ctrl, method="REML"))
## ...and fit the AR3
m5.Diff.month <- try(gamm(Europe ~ s(YearFrac, k = 10) + s(Month, bs = "cc", k = 12) + Data, data = Diff.df, correlation = corARMA(form = ~ 1|Year, p = 3),   control = ctrl, method="REML"))	

anova(m1.Diff$lme, m2.Diff$lme)	
#### Adding monthly doesn't do anything   
anova(m1.Diff$lme, m3.Diff.month$lme, m4.Diff.month$lme, m5.Diff.month$lme)	
 
#### M4 seems best again, try without monthly term
# ...so fit the AR1
m3.Diff <- try(gamm(Europe ~ s(YearFrac, k = 10) + Data, data = Diff.df, correlation = corARMA(form = ~ 1|Year, p = 1),   control = ctrl, method="REML"))
## ...and fit the AR2
m4.Diff <- try(gamm(Europe ~ s(YearFrac, k = 10) + Data, data = Diff.df, correlation = corARMA(form = ~ 1|Year, p = 2),   control = ctrl, method="REML"))
## ...and fit the AR3
m5.Diff <- try(gamm(Europe ~ s(YearFrac, k = 10) + Data, data = Diff.df, correlation = corARMA(form = ~ 1|Year, p = 3),   control = ctrl, method="REML"))	

anova(m1.Diff$lme, m3.Diff$lme, m4.Diff$lme, m5.Diff$lme)	
AIC(m1.Diff$lme, m3.Diff$lme, m4.Diff$lme, m5.Diff$lme, m3.Diff.month$lme, m4.Diff.month$lme, m5.Diff.month$lme)  
BIC(m1.Diff$lme, m3.Diff$lme, m4.Diff$lme, m5.Diff$lme, m3.Diff.month$lme, m4.Diff.month$lme, m5.Diff.month$lme)  
  
### Looks like AR2 is the best
m.final <- m4.Diff.month

summary(m.final$gam)
plot(m.final$gam, residuals=TRUE)

## look at autocorrelation in residuals:
acf(resid(m.final$lme, type = "normalized"))
pacf(resid(m.final$lme, type = "normalized"))
## No more temporal autocorrelation

#################################
###  Calculate significance
#########################################
#### No significance test
### Year has a p-vale of 0.0205  < 0.025
### Significant trend
month.list <- seq(1,12,0.01)
test2.data <- data.frame(Month=month.list, Data="WFD", YearFrac=0)
p1 <- predict(m.final$gam, newdata = test2.data)
p1 <- p1 - median(p1)
plot(month.list, p1, type="l")

month.list[which.min(abs(p1))]
month.list[which.min(p1)]
month.list[which.max(p1)]


just.WFD <- subset(Diff.df, Data=="WFD")
just.WFD$Month <- 6.68
mfinal.d <- Deriv(m.final$gam, newdata = just.WFD)
plot(mfinal.d, sizer = TRUE, alpha = 0.05)

plot(Europe ~ Date, data = just.WFD, type = "p")
p1 <- predict(m.final$gam, newdata = just.WFD)
lines(p1 ~ just.WFD$Date, col="black")

CI <- confint(mfinal.d, alpha = 0.05)
S <- signifD(p1, mfinal.d$Year$deriv, CI$Year$upper, CI$Year$lower, eval = 0)
lines(S$incr ~ Date, data = just.WFD, lwd = 3, col = "red")
lines(S$decr ~ Date, data = just.WFD, lwd = 3, col = "red")
Pred.WFD <- data.frame(Date=just.WFD$Date, Pred=p1, Inc=S$incr, Dec=S$decr, Data="WFD")

just.WFDEI <- subset(Diff.df, Data=="WFDEI")
just.WFDEI$Month <- 6.68
mfinal.d <- Deriv(m.final$gam, newdata = just.WFDEI)
plot(mfinal.d, sizer = TRUE, alpha = 0.05)

plot(Europe ~ Date, data = just.WFDEI, type = "p")
p1 <- predict(m.final$gam, newdata = just.WFDEI)
lines(p1 ~ just.WFDEI$Date, col="black")

CI <- confint(mfinal.d, alpha = 0.05)
S <- signifD(p1, mfinal.d$Year$deriv, CI$Year$upper, CI$Year$lower, eval = 0)
lines(S$incr ~ Date, data = just.WFDEI, lwd = 3, col = "red")
lines(S$decr ~ Date, data = just.WFDEI, lwd = 3, col = "red")
Pred.WFDEI <- data.frame(Date=just.WFDEI$Date, Pred=p1, Inc=S$incr, Dec=S$decr, Data="WFDEI")
Pred.df <- rbind(Pred.WFD, Pred.WFDEI)



Pred.comb <- Diff.df[!duplicated(Diff.df$Date),]
Pred.comb$Data <- "WFD"
Pred.comb$Month <- 6.68
mfinal.d <- Deriv(m.final$gam, newdata = Pred.comb)
plot(mfinal.d, sizer = TRUE, alpha = 0.05)

plot(Europe ~ Date, data = Pred.comb, type = "p")
p1 <- predict(m.final$gam, newdata = Pred.comb, type="response")
lines(p1 ~ Pred.comb$Date, col="black")

CI <- confint(mfinal.d, alpha = 0.05)
S <- signifD(p1, mfinal.d$Year$deriv, CI$Year$upper, CI$Year$lower, eval = 0)
lines(S$incr ~ Date, data = Pred.comb, lwd = 3, col = "red")
lines(S$decr ~ Date, data = Pred.comb, lwd = 3, col = "blue")
Pred.comb <- data.frame(Date=Pred.comb$Date, Pred=p1, Inc=S$incr, Dec=S$decr, Data="WFD")



#######################################

p <- ggplot(Diff.df, aes(x=Date, group=Data))
p <- p + geom_hline(yintercept=0, linetype="longdash", colour="grey30")
p <- p + geom_line(aes(y=Europe), colour="grey70", size=0.3)
#p <- p + geom_line(aes(y=Europe, colour=Data), size=0.3,alpha=0.7)
p <- p + geom_line(data=Pred.df, aes(y=Pred), colour="black", size = 0.75, colour="grey20")
p <- p + geom_line(data=Pred.df, aes(y=Inc), colour="red", size=1)
p <- p + geom_line(data=Pred.df, aes(y=Dec), colour="red", size=1)
#p <- p + geom_segment(aes(x = as.Date("1971-01-01"), y = 0.025, xend = as.Date("2000-12-31"), yend = 0.025), arrow=arrow(length = unit(0.1, "inches"), angle = 15, ends = "both", type = "closed"))
p <- p + geom_segment(aes(x = as.Date("1971-01-01"), y = -0.09, xend = as.Date("2000-12-31"), yend = -0.09), size=2, alpha=0.8, colour="grey80")
p <- p + scale_y_continuous(name="SPEI6 Drought Area", labels = percent)
p <- p + theme_classic(9)
p <- p + theme(legend.position=c(0.25,0.85), legend.direction="horizontal")
p <- p + scale_colour_brewer(type="qual", palette="Set1")
p <- p + coord_cartesian(ylim=c(-0.1,0.3))
p

ggsave(file.path(write_fig_path, "Diff6_Version1.png"), p, width=7, height=4, dpi=400)
ggsave(file.path(write_fig_path, "Diff6_Version1.pdf"), p, width=7, height=4)

#################################
## Redo WFDEI with adjustment for difference in data 
#################################
just.WFDEI <- subset(Diff.df, Data=="WFDEI")
just.WFDEI$Data <- "WFD"
just.WFDEI$Month <- 6.68

p1 <- predict(m.final$gam, newdata = just.WFDEI, type="response")

mfinal.d <- Deriv(m.final$gam, newdata = just.WFDEI)
CI <- confint(mfinal.d, alpha = 0.05)
S <- signifD(p1, mfinal.d$Year$deriv, CI$Year$upper, CI$Year$lower, eval = 0)

WFDEI.adj <- just.WFDEI$Europe
WFDEI.adj <- WFDEI.adj - m.final$gam$coefficients[[2]]

just.WFDEI$Europe <- WFDEI.adj
just.WFDEI$Data <- "WFDEI"
Diff.df.adj <- rbind(just.WFD, just.WFDEI)

Pred.WFDEI <- data.frame(Date=just.WFDEI$Date, Pred=p1, Inc=S$incr, Dec=S$decr, Data="WFDEI")
Pred.df <- rbind(Pred.WFD, Pred.WFDEI)


p <- ggplot(Diff.df.adj, aes(x=Date, group=Data))
p <- p + geom_hline(yintercept=0, linetype="longdash", colour="grey30")
p <- p + geom_line(aes(y=Europe), colour="grey70", size=0.3)
#p <- p + geom_line(aes(y=Europe, colour=Data), size=0.3,alpha=0.7)
p <- p + geom_line(data=Pred.comb, aes(y=Pred), colour="black", size = 0.75, colour="grey20")
p <- p + geom_line(data=Pred.comb, aes(y=Inc), colour="red", size=1)
p <- p + geom_line(data=Pred.comb, aes(y=Dec), colour="red", size=1)
#p <- p + geom_segment(aes(x = as.Date("1971-01-01"), y = 0.025, xend = as.Date("2000-12-31"), yend = 0.025), arrow=arrow(length = unit(0.1, "inches"), angle = 15, ends = "both", type = "closed"))
p <- p + geom_segment(aes(x = as.Date("1971-01-01"), y = -0.11, xend = as.Date("2000-12-31"), yend = -0.11), size=2, alpha=0.8, colour="grey80")
p <- p + scale_y_continuous(name="Diff6 Drought Area", labels = percent, breaks=seq(-1,1,0.05))
p <- p + theme_classic(9)
p <- p + theme(legend.position=c(0.25,0.85), legend.direction="horizontal")
p <- p + scale_colour_brewer(type="qual", palette="Set1")
p <- p + coord_cartesian(ylim=c(-0.12,0.21))
p

ggsave(file.path(write_fig_path, "Diff6_Version1_adj.png"), p, width=7, height=4, dpi=400)
ggsave(file.path(write_fig_path, "Diff6_Version1_adj.pdf"), p, width=7, height=4)



