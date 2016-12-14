library(randomForest)
library(reshape)
library(plyr)

options(stringsAsFactors = FALSE)

#RUN1
vi_median_name <- "bsti_vi_median_20161213_1603.Rdata"
bsti_name <- "bsti_20161213_1603.Rdata"

#RUN2
# vi_median_name <- "bsti_vi_median_20160205_1147.Rdata"
# bsti_name <- "bsti_20160205_1147.Rdata"

load(vi_median_name)
load(bsti_name)

#### Variable selection ####
# Values drop off and then level out. Arbitrarily going with 50% of the variables.
# grab all variable names with median values > 1.004880e-04 = 50% of the data
# This 50% of the data reflects 50% of the original list of variables prior to scaling
# Scaling had the effect of dropping variables that were all 0s anyway.
bsti.s2.col <- bsti.vi.median[1:ceiling(nrow(bsti.vi.median)) / 2, ][, 1]
bsti.vi.median <- bsti.vi.median[1:ceiling(nrow(bsti.vi.median)) / 2, ]
#bsti.s2.col <- c("FSS_26Aug14",(bsti.vi.median[,'var_name']))
bsti.s2 <- bsti[, colnames(bsti) %in% c(bsti.s2.col,'HDWTR')]

# bsti.s2.col <- vars[vars$var %in% names(bsti.s2),]
# #bsti.s2.col <- bsti.s2.col[bsti.s2.col$var != 'FSS_26Aug14',]
# bsti.s2.col <- merge(bsti.s2.col, bsti.vi.median[, c('var_name','median')], 
#                      by.x = 'var', by.y = 'var_name', all.x = TRUE)
# bsti.s2.col <- arrange(bsti.s2.col, desc(median))


# # #All together
correlation_threshold <- 0.4
pcor <- cor(bsti[,setdiff(bsti.vi.median$var_name,"DATE")])
pnames <- attr(pcor, "dimnames")[[1]]
pkeep <- pnames
tracking <- vector("list", length(pnames))
names(tracking) <- pnames
for (i in length(pnames):1) {
  tracking[[i]] <- pcor[order(abs(pcor[,i])),i]#[apply(pcor, MARGIN = 1, FUN = function(x) abs(x) >= correlation_threshold)[,i]]
  if (any(round(abs(pcor[i,][-i]),2) >= correlation_threshold)) {
    if (i != 1) {
      pkeep <- pkeep[-i]
      pcor <- pcor[-i,-i,drop=FALSE]
    }
  } 
}

#Further remove variables to reduce the influence of correlation on raising variable importance
#bsti.s2 <- bsti.s2[,colnames(bsti.s2) %in% keeps.s2]
bsti.s2 <- bsti.s2[, c(pkeep,'HDWTR')]
colnames(bsti.s2)


# remove any NAs
bsti.s2 <- data.frame(na.omit(bsti.s2))

