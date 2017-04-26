### This function clears the memory in R, in case there are old variables from previous work.
rm(list=ls())



require(ggplot2)
require(grid)
require(zoo)
require(scales)
require(gridExtra)
require(zyp)
require(Kendall)

require(rkt)
require(lubridate)

## load the packages and code we need
require(mgcv)
#require(nlme)

#source("gam_method_month.R")
source("derivFun.R")


###########################################################################
##  Prepare time sequence without leap days
###########################################################################
WFD.dates <- seq(as.Date("1958/1/1"), as.Date("2001/12/31"), "days")
WFDEI.dates <- seq(as.Date("1979/1/1"), as.Date("2014/12/31"), "days")
Comb.dates <- seq(as.Date("1958/1/1"), as.Date("2014/12/31"), "days")

WFD.month.day <- substring(WFD.dates,6,10)
WFDEI.month.day <- substring(WFDEI.dates,6,10)
Comb.month.day <- substring(Comb.dates,6,10)

### Create time index with no leap years
WFD.noleap <- WFD.dates[WFD.month.day != "02-29"]
WFDEI.noleap <- WFDEI.dates[WFDEI.month.day != "02-29"]
Comb.noleap <- Comb.dates[Comb.month.day != "02-29"]

Comb.months <- as.numeric(substr(Comb.noleap, 6,7))

monthly.dates <- c(Comb.noleap[which(diff(Comb.months,1)!=0)], as.Date("2014-12-31"))
 
################################################
###  Read in Data
#################################################
#### Read in SPI
Tair.df <- read.csv("Tair_Option2_Europe.csv")
Tair.df$Date <- as.Date(Tair.df$Date)

Tair.df$Data <- factor(Tair.df$Data)
Tair.df$YearFrac <- as.numeric(Tair.df$Date)
Tair.df$Year <- as.numeric(substr(Tair.df$Date,1,4))
Tair.df$Month <- as.numeric(substring(Tair.df$Date,6,7))

monthly.test <- Tair.df$Date %in% monthly.dates
Tair.monthly.df <- Tair.df[monthly.test,]


################################################
###  Fit Tair Model
#################################################
## Fit a smoother for Year to the data
# Optimizer with upper/lower bounds
#ctrl <- list(niterEM = 0, msVerbose = TRUE, optimMethod="L-BFGS-B")
ctrl <- list(niterEM = 0, msVerbose = TRUE, optimMethod="BFGS")

### Fit plain model
m1 <- gamm(Value ~ s(YearFrac, k = 10) + Data, data = Tair.monthly.df,   control = ctrl, method="REML")
summary(m1$gam)
plot(m1$gam, residuals=TRUE)

## Account for potential annual cycle
m2 <- gamm(Value ~ s(YearFrac, k = 10) + s(Month, bs = "cc", k = 12) + Data, data = Tair.monthly.df,   control = ctrl, method="REML")
summary(m2$gam)
plot(m2$gam, residuals=TRUE)
 
## Compare models
anova(m1$lme, m2$lme)	
AIC(m1$lme, m2$lme)	 
BIC(m1$lme, m2$lme)	 
  
#### Looks like no significnat monthly pattern  

## look at autocorrelation in residuals:
acf(resid(m2$lme, type = "normalized"))
## ...wait... look at the plot, only then do...
pacf(resid(m2$lme, type = "normalized"))
## seems like like sharp cut-off in ACF and PACF - AR terms probably best
 
 
## Fit the AR1
m3.monthly <- try(gamm(Value ~ s(YearFrac, k = 10) + s(Month, bs = "cc", k = 12) + Data, data = Tair.monthly.df, correlation = corARMA(form = ~ 1|Year, p = 1),   control = ctrl, method="REML"))
## ...and fit the AR2
m4.monthly <- try(gamm(Value ~ s(YearFrac, k = 10) + s(Month, bs = "cc", k = 12) + Data, data = Tair.monthly.df, correlation = corARMA(form = ~ 1|Year, p = 2),   control = ctrl, method="REML"))
## ...and fit the AR3
m5.monthly <- try(gamm(Value ~ s(YearFrac, k = 10) + s(Month, bs = "cc", k = 12) + Data, data = Tair.monthly.df, correlation = corARMA(form = ~ 1|Year, p = 3),   control = ctrl, method="REML"))

anova(m1$lme, m2$lme)	
#### Adding monthly doesn't do anything   
anova(m1$lme, m3.monthly$lme, m4.monthly$lme, m5.monthly$lme)	
AIC(m1$lme, m3.monthly$lme, m4.monthly$lme, m5.monthly$lme)	
BIC(m1$lme, m3.monthly$lme, m4.monthly$lme, m5.monthly$lme)	




##################################
####  Try without monthly term
####################################
## Fit the AR1
m3 <- try(gamm(Value ~ s(YearFrac, k = 10) + Data, data = Tair.monthly.df, correlation = corARMA(form = ~ 1|Year, p = 1),   control = ctrl, method="REML"))
## ...and fit the AR2
m4 <- try(gamm(Value ~ s(YearFrac, k = 10) + Data, data = Tair.monthly.df, correlation = corARMA(form = ~ 1|Year, p = 2),   control = ctrl, method="REML"))
## ...and fit the AR3
m5 <- try(gamm(Value ~ s(YearFrac, k = 10) + Data, data = Tair.monthly.df, correlation = corARMA(form = ~ 1|Year, p = 3),   control = ctrl, method="REML"))

anova(m1$lme, m2$lme)	
#### Adding monthly doesn't do anything  
anova(m1$lme, m3$lme,  m4$lme, m5$lme)	
AIC(m1$lme, m3$lme,  m4$lme, m5$lme)	
BIC(m1$lme, m3$lme,  m4$lme, m5$lme)	

anova(m1$lme, m4$lme,  m4.monthly$lme)	
AIC(m1$lme, m3$lme,  m3.monthly$lme, m4$lme, m4.monthly$lme,m5$lme, m5.monthly$lme)	
BIC(m1$lme, m3$lme,  m3.monthly$lme, m4$lme, m4.monthly$lme,m5$lme, m5.monthly$lme)	


m4.yup <- try(gamm(Value ~ s(YearFrac, k = 10), data = Tair.monthly.df, correlation = corARMA(form = ~ 1|Year, p = 2),   control = ctrl, method="REML"))
anova(m1$lme, m4.yup$lme, m4$lme,  m4.monthly$lme)	

##################################
####  Calculate final model
####################################
m.final <- m4.yup
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

just.WFD <- subset(Tair.df, Data=="WFD")
mfinal.d <- Deriv(m.final$gam, newdata = just.WFD)
plot(mfinal.d, sizer = TRUE, alpha = 0.05)

plot(Value ~ Date, data = just.WFD, type = "p")
p1 <- predict(m.final$gam, newdata = just.WFD)
lines(p1 ~ just.WFD$Date, col="black")

CI <- confint(mfinal.d, alpha = 0.05)
S <- signifD(p1, mfinal.d$Year$deriv, CI$Year$upper, CI$Year$lower, eval = 0)
lines(S$incr ~ Date, data = just.WFD, lwd = 3, col = "blue")
lines(S$decr ~ Date, data = just.WFD, lwd = 3, col = "red")
Pred.WFD <- data.frame(Date=just.WFD$Date, Pred=p1, Inc=S$incr, Dec=S$decr, Data="WFD")

just.WFDEI <- subset(Tair.df, Data=="WFDEI")
mfinal.d <- Deriv(m.final$gam, newdata = just.WFDEI)
plot(mfinal.d, sizer = TRUE, alpha = 0.05)

plot(Value ~ Date, data = just.WFDEI, type = "p")
p1 <- predict(m.final$gam, newdata = just.WFDEI)
lines(p1 ~ just.WFDEI$Date, col="black")

CI <- confint(mfinal.d, alpha = 0.05)
S <- signifD(p1, mfinal.d$Year$deriv, CI$Year$upper, CI$Year$lower, eval = 0)
lines(S$incr ~ Date, data = just.WFDEI, lwd = 3, col = "blue")
lines(S$decr ~ Date, data = just.WFDEI, lwd = 3, col = "red")
Pred.WFDEI <- data.frame(Date=just.WFDEI$Date, Pred=p1, Inc=S$incr, Dec=S$decr, Data="WFDEI")
Pred.df <- rbind(Pred.WFD, Pred.WFDEI)



which.min(Pred.df$Inc)
Pred.df[7721,]

Pred.df$Pred[29200] - Pred.df$Pred[7721]
Pred.df$Date[29200] - Pred.df$Date[7721]
(1.401728/(13093/365.25))*10

Tair.df$Date[23726:29200][Tair.df$Value[23726:29200] < 0]



Pred.comb <- Tair.monthly.df[!duplicated(Tair.monthly.df$Date),]
Pred.comb$Data <- "WFD"
Pred.comb$Month <- zero.point
mfinal.d <- Deriv(m.final$gam, newdata = Pred.comb)
plot(mfinal.d, sizer = TRUE, alpha = 0.05)

plot(Value ~ Date, data = Pred.comb, type = "p")
p1 <- predict(m.final$gam, newdata = Pred.comb, type="response")
lines(p1 ~ Pred.comb$Date, col="black")

CI <- confint(mfinal.d, alpha = 0.05)
S <- signifD(p1, mfinal.d$Year$deriv, CI$Year$upper, CI$Year$lower, eval = 0)
lines(S$incr ~ Date, data = Pred.comb, lwd = 3, col = "red")
lines(S$decr ~ Date, data = Pred.comb, lwd = 3, col = "red")
Pred.comb <- data.frame(Date=Pred.comb$Date, Pred=p1, Inc=S$incr, Dec=S$decr, Data="WFD")



#######################################

p <- ggplot(Tair.monthly.df, aes(x=Date, group=Data))
p <- p + geom_hline(yintercept=0, linetype="longdash", colour="grey30")
p <- p + geom_line(aes(y=Value), colour="grey70", size=0.3)
#p <- p + geom_line(aes(y=Europe, colour=Data), size=0.3,alpha=0.7)
p <- p + geom_line(data=Pred.df, aes(y=Pred), colour="black", size = 0.75, colour="grey20")
p <- p + geom_line(data=Pred.df, aes(y=Inc), colour="red", size=1)
p <- p + geom_line(data=Pred.df, aes(y=Dec), colour="red", size=1)
p <- p + geom_segment(aes(x = as.Date("1971-01-01"), y = -2, xend = as.Date("2000-12-31"), yend = -2), size=2, alpha=0.8, colour="grey80")
p <- p + scale_y_continuous(name="European 6 month Tavg anomaly (C)")
p <- p + theme_classic(9)
p <- p + theme(legend.position=c(0.25,0.85), legend.direction="horizontal")
p <- p + scale_colour_brewer(type="qual", palette="Set1")
p <- p + coord_cartesian(ylim=c(-2.2, 2.3))
p

ggsave("Tair6_Version1.png", p, width=7, height=4, dpi=400)
ggsave("Tair6_Version1.pdf", p, width=7, height=4)

p <- ggplot(Tair.monthly.df, aes(x=Date, group=Data))
p <- p + geom_hline(yintercept=0, linetype="longdash", colour="grey30")
#p <- p + geom_line(aes(y=Europe), colour="grey70", size=0.3)
p <- p + geom_line(aes(y=Value, colour=Data), size=0.3,alpha=0.7)
p <- p + geom_line(data=Pred.df, aes(y=Pred), size = 0.75, colour="grey20", linetype="dashed")
p <- p + geom_line(data=Pred.df, aes(y=Inc), colour="black", size=1)
p <- p + geom_line(data=Pred.df, aes(y=Dec), colour="black", size=1)
p <- p + geom_segment(aes(x = as.Date("1971-01-01"), y = -2, xend = as.Date("2000-12-31"), yend = -2), size=2, alpha=0.8, colour="grey80")
p <- p + scale_y_continuous(name="European 6 month Tavg anomaly (C)")
p <- p + theme_classic(9)
p <- p + theme(legend.position=c(0.25,0.85), legend.direction="horizontal")
p <- p + scale_colour_brewer(type="qual", palette="Set1")
p <- p + coord_cartesian(ylim=c(-2.2, 2.3))
p

ggsave("Tair6_Version2.png", p, width=7, height=4, dpi=400)
ggsave("Tair6_Version2.pdf", p, width=7, height=4)





p <- ggplot(Tair.monthly.df, aes(x=Date, group=Data))
p <- p + geom_hline(yintercept=0, linetype="longdash", colour="grey30")
#p <- p + geom_line(aes(y=Europe), colour="grey70", size=0.3)
p <- p + geom_line(aes(y=Value, colour=Data), size=0.5,alpha=0.7)
p <- p + geom_line(data=Pred.comb, aes(y=Pred), size = 0.9, colour="grey20", linetype="dashed")
p <- p + geom_line(data=Pred.comb, aes(y=Inc), colour="black", size=1.2)
p <- p + geom_line(data=Pred.comb, aes(y=Dec), colour="black", size=1.2)
p <- p + geom_segment(aes(x = as.Date("1971-01-01"), y = -2, xend = as.Date("2000-12-31"), yend = -2), size=2, alpha=0.8, colour="grey80")
p <- p + scale_y_continuous(name=expression(paste("6 month ",T[mean], "  anomaly (C)", sep="")))
p <- p + theme_classic(11) + theme(
axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
p <- p +  theme(axis.title.y=element_text(margin=margin(0,12,0,0)))
p <- p + theme(legend.position="none")
p <- p + scale_colour_brewer(type="qual", palette="Set1")
p <- p + coord_cartesian(ylim=c(-2.2, 2.3))
p

ggsave("Tair6_Version2_poster2.png", p, width=7, height=4, dpi=400)
ggsave("Tair6_Version2_poster.pdf", p, width=7, height=4)




p <- ggplot(Tair.monthly.df, aes(x=Date, group=Data))
p <- p + geom_hline(yintercept=0, linetype="longdash", colour="grey30")
p <- p + geom_line(aes(y=Value), colour="grey70", size=0.3)
#p <- p + geom_line(aes(y=Europe, colour=Data), size=0.3,alpha=0.7)
p <- p + geom_line(data=Pred.df, aes(y=Pred), size = 0.75, colour="grey20", linetype="dashed")
p <- p + geom_line(data=Pred.df, aes(y=Inc), colour="black", size=1)
p <- p + geom_line(data=Pred.df, aes(y=Dec), colour="black", size=1)
p <- p + geom_segment(aes(x = as.Date("1971-01-01"), y = -2, xend = as.Date("2000-12-31"), yend = -2), size=2, alpha=0.8, colour="grey80")
p <- p + scale_y_continuous(name="European 6 month Tavg anomaly (C)")
p <- p + theme_classic(9)
p <- p + theme(legend.position=c(0.25,0.85), legend.direction="horizontal")
p <- p + scale_colour_brewer(type="qual", palette="Set1")
p <- p + coord_cartesian(ylim=c(-2.2, 2.3))
p


ggsave("Tair6_Version3.png", p, width=7, height=4, dpi=400)
ggsave("Tair6_Version3.pdf", p, width=7, height=4)

