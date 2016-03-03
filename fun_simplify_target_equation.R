simplify_target_equation <- function(betahat, ss, station) {
  options(warn = -1)
  bsub <- betahat[betahat$STATION_KEY == station,]
  bsub <- plyr::rename(bsub, c('HDWTR100' = 'HDWTR'))
  bsubm <- melt(bsub, id.vars = 'STATION_KEY')
  vals <- ss[ss$STATION_KEY == station, names(bsub)[3:11]]
  vals <- melt(vals, measure.vars = 1:9)
  inter <- merge(bsubm, vals, by = 'variable', suffixes = c('.betahat',''))
  inter$inter <- inter$value.betahat*as.numeric(inter$value)
  BSTI <- ss[ss$STATION_KEY == station, 'log10_BSTI']
  Z <- BSTI - bsub$`(Intercept)`[1] - sum(inter$inter)
  
  inter2 <- inter[inter$variable != 'sum_1095_days',]
  b <- bsub$`(Intercept)`[1] + sum(inter2$inter) + Z
  b_u <- (b*max_log10_bsti)/100
  m_u <- (bsub$sum_1095_days*max_log10_bsti)/100
  
  return(list("b_u" = b_u, "m_u" = m_u))
  options(warn = 0)
}

