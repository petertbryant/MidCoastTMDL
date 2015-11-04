source('funCorrelationPlots.R')

#### Correlation plots ####
#names(fss2.s2.col) <- fss2.s2.col

#Method: Look for r<=0.2. If lower importance variable has r>=0.2 with any higher importance variables select the highest 
#importance variable unless that variable is already conceptually represented.

# Precip  
#png('precip_cor.png')
pairs(fss2.s2[, fss2.s2.col[fss2.s2.col$Category == 'Precipitation', 'var']],
      lower.panel = panel.smooth, upper.panel = panel.cor, 
      diag.panel = panel.hist)
#dev.off()
# keep "sum_1095_days"

#Disturbance
pairs(fss2.s2[ ,fss2.s2.col[fss2.s2.col$Category == 'Disturbance', 'var']],
      lower.panel = panel.smooth, upper.panel = panel.cor, 
      diag.panel = panel.hist)
# everything is coorelated. checked the median sd of importance scores for
# disturbance variables and no overlap exists between DIS_1YR_PARSA (vi = 2.84,
# visd = 0.17) and DIS_10YR_PRSA (vi = 1.21, visd = 0.16)
# keep "DIS_1YR_PARSA"

# Lithology/soils
pairs(fss2.s2[, fss2.s2.col[fss2.s2.col$Category == 'Lithology and soils', 
                            'var']], 
      lower.panel = panel.smooth, upper.panel = panel.cor, 
      diag.panel = panel.hist)
# almost everything is coorelated
# Keep "EROD_PARCA", "SILT_CLAY_PARCA"

# Ownership 
pairs(fss2.s2[, fss2.s2.col[fss2.s2.col$Category == 'Land use', 'var']], 
      lower.panel = panel.smooth, upper.panel = panel.cor, 
      diag.panel = panel.hist)
# Keep"ROADLEN_DRSA","OWN_FED_PRCA","POP_DARCA","OWN_PRI_PRCA","OWN_AGR_PARCA"

#Susceptibility
pairs(fss2.s2[, fss2.s2.col[fss2.s2.col$Category == 'Landslide susceptibility', 
                            'var']], 
      lower.panel = panel.smooth, upper.panel = panel.cor,
      diag.panel = panel.hist)
#All correlated. 
#Keep "SUSCEP5_PARCA"

#Location
pairs(fss2.s2[, fss2.s2.col[fss2.s2.col$Category == 'Location', 'var']], 
      lower.panel = panel.smooth, upper.panel = panel.cor,
      diag.panel = panel.hist)
#Keep "MIN_Z"

#Stream attributes
pairs(fss2.s2[, fss2.s2.col[fss2.s2.col$Category == 'Stream attributes', 
                            'var']],
      lower.panel = panel.smooth, upper.panel = panel.cor,
      diag.panel = panel.hist)
#Keep "STRMPWR","XSLOPE_MAP"