simplify_target_equation <- function(betahat, ss, station) {
  options(warn = -1)
  betahat <- plyr::rename(betahat, c('HDWTR100' = 'HDWTR'))
  bsubm <- melt(betahat)
  vals <- ss[ss$STATION_KEY == station, names(betahat)[2:10]]
  vals <- melt(vals, measure.vars = 1:9)
  vals$value <- as.numeric(vals$value)
  inter <- merge(bsubm, vals, by = 'variable', suffixes = c('.betahat',''))
  inter$inter <- inter$value.betahat * inter$value
  BSTI <- ss[ss$STATION_KEY == station, 'log10_BSTI']
  Z <- BSTI - betahat$`(Intercept)`[1] - sum(inter$inter)
  
  inter2 <- inter[inter$variable != 'sum_1095_days',]
  b <- betahat$`(Intercept)`[1] + sum(inter2$inter) + Z
  b_u <- (b*max_log10_bsti)/100
  m_u <- (betahat$sum_1095_days*max_log10_bsti)/100
  
  return(list("b_u" = b_u, "m_u" = m_u))
  options(warn = 0)
}

