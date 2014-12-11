library(reshape2)
library(plyr)
library(ggplot2)

t <- data.frame(h = seq(0,5,.01))
t$ls <- ifelse(t$h < 1, 1-t$h, 0)
t$s <- ifelse(t$h < 1, 1-(3/2)*t$h+(1/2)*t$h^3, 0)
t$m <- ifelse(t$h == 0, 1, log(t$h +1)/t$h)
t$e <- exp(-t$h)

tm <- melt(t, id.vars = c('h'))
tm$variable <- mapvalues(tm$variable, from = c('ls','s','m','e'), to = c('Linear with sill', 'Spherical', 'MARIAH', 'Exponential'))

#Plot with base graphics
png('AutoCorFunPlots.png', width = 8, height = 6, units = "in", res = 300)
par(mfrow=c(2,2))
plot(t$h,t$ls,type="l",main="Linear with sill",ylim=c(0,1),xlab="Hydrologic distance",ylab='Autocorrelation')
plot(t$h, t$s,type = "l",main="Spherical",ylim=c(0,1),xlab="Hydrologic distance",ylab='Autocorrelation')
plot(t$h, t$m,type = "l",main="MARIAH",ylim=c(0,1),xlab="Hydrologic distance",ylab='Autocorrelation')
plot(t$h, t$e,type="l",main="Exponential",ylim=c(0,1),xlab="Hydrologic distance",ylab='Autocorrelation')
dev.off()

#Same plot with ggplot
p <- ggplot(tm) 
p + geom_line(aes(h, value), size = 1) + 
  facet_wrap(~ variable) + 
  xlab("Hydrologic distance") +
  ylab("Autocorrelation") + 
  theme_bw()
ggsave("AutoCorFunPlots_ggplot.png", width = 8, height = 6, units = "in")


