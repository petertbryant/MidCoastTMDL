library(party)
library(reshape)

options(stringsAsFactors = FALSE)
vars <- read.csv("VarNames_RF.csv")
bugs <- read.csv("ssn_RF_data.csv")

# -----------------------------------------------------------
# Correlation Matrix Functions
# source
# http://stackoverflow.com/questions/15271103/how-to-modify-this-correlation-matrix-plot
# -----------------------------------------------------------
panel.cor <- function(x, y, digits=2, cex.cor)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  test <- cor.test(x,y)
  Signif <- ifelse(round(test$p.value,3)<0.001,"p<0.001",paste("p=",round(test$p.value,3)))  
  text(0.5, 0.25, paste("r=",txt))
  text(.5, .75, Signif)
}

panel.smooth<-function (x, y, col = "black", bg = NA, pch = 18, 
                        cex = 0.8, col.smooth = "red", span = 2/3, iter = 3, ...) 
{
  points(x, y, pch = pch, col = col, bg = bg, cex = cex)
  ok <- is.finite(x) & is.finite(y)
  if (any(ok)) 
    lines(stats::lowess(x[ok], y[ok], f = span, iter = iter), 
          col = col.smooth, ...)
}

panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col="gray", ...)
}

# ----------------------------------------------------------- #
# FSS1 - Random forests including the physical habitat data
# ----------------------------------------------------------- #
# Step 1, per Genuer et al (2000), use brute force to identify important variables. 
# Since we have a large amount of variables random forest is run 50 times with a 
# very high number of trees per forest (ntree). This yields a distribution of importance scores.
# Removal is based on these distributions. We keep the variables with the highest scores.

vars.fss1.s1 <- vars[vars$fss1.rf_keep == 1,]
fss1.s1 <- bugs[,colnames(bugs) %in% vars.fss1.s1$var]

# remove NAs in response variable
fss1.s1 <- fss1.s1[(!is.na(fss1.s1$FSS_26Aug14)),]

# remove any NAs
#fss1.s1 <- data.frame(na.omit(fss1.s1))

colnames(fss1.s1)

# mtry and ntree values 
mtry.fss1.s1 <- as.integer(((ncol(fss1.s1)-1) / 3),0)

# initialize the variable importance df
fss1.s1.vi <- data.frame(matrix(, nrow = ncol(fss1.s1)-1, ncol = 50))
fss1.s1.col <- colnames(fss1.s1)
fss1.s1.col <- fss1.s1.col[!(fss1.s1.col == "FSS_26Aug14")]

# WARNING - Takes a long time to run.
print(Sys.time())
set.seed(42)
for (i in 1:50) {
  print(i)
  fss1.s1.cf <- cforest(FSS_26Aug14 ~ ., data = fss1.s1, controls = cforest_unbiased(ntree = 2000, mtry = mtry.fss1.s1))
  fss1.s1.vi[,i]<- varimp(fss1.s1.cf, conditional=FALSE)
}
print(Sys.time())

# Add var names and index
fss1.s1.vi[,51]<- fss.col
fss1.s1.vi[,52]<-c(1:length(fss.col))
colnames(fss1.s1.vi)[51] <- "var_name"
colnames(fss1.s1.vi)[52] <- "var_index"



# ----------
# Save the df with a timestamp so we don't accidently overwrite it.
timestamp <- format(Sys.time(), "%Y%m%d_%H%M")
save(fss1.s1.vi, file=paste0("fss1_vi_s1_",timestamp,".RData"))

#load("fss1_vi_s1_20141019_1451.RData")
# ----------

fss1.s1.vi.l <- melt(fss1.s1.vi, id=c("var_name","var_index"))

bymedian <- with(fss1.s1.vi.l, reorder(var_index, -value, median))
boxplot(value ~ bymedian, data = fss1.s1.vi.l,
        xlab = "Variable index", ylab = "Importance", 
        varwidth = TRUE,
        col = "lightgray")

fss1.s1.vi.median <- cast(fss1.s1.vi.l,var_name + var_index ~ ., value ='value', median)
colnames(fss1.s1.vi.median )[3] <- "median"

# sort the data so largest at the top
fss1.s1.vi.median <- fss1.s1.vi.median[with(fss1.s1.vi.median, order(-median)), ]



# ----------------------------------------------------------- #
# FSS2 - Random forests excluding the physical habitat data
# ----------------------------------------------------------- #
# Step 1, per Genuer et al (2000), use brute force to identify important variables. 
# Since we have a large amount of variables random forest is run 50 times with a 
# very high number of trees per forest (ntree). This yields a distribution of importance scores.
# Removal is based on these distributions. We keep the variables with the highest scores.

vars.fss2.s1 <- vars[vars$fss2.rf_keep == 1,]
fss2.s1 <- bugs[,colnames(bugs) %in% vars.fss2.s1$var]

# remove NAs in response variable
fss2.s1 <- fss2.s1[(!is.na(fss2.s1$FSS_26Aug14)),]

colnames(fss2.s1)

# mtry and ntree values 
mtry.fss2.s1 <- as.integer(((ncol(fss2.s1)-1) / 3),0)

# initialize the variable importance df
fss2.s1.vi <- data.frame(matrix(, nrow = ncol(fss2.s1)-1, ncol = 50))
fss2.s1.col <- colnames(fss2.s1)
fss2.s1.col <- fss2.s1.col[!(fss2.s1.col == "FSS_26Aug14")]

# WARNING - Takes a long time to run.
set.seed(42)
for (i in 1:50) {
  print(i)
  fss2.s1.cf <- cforest(FSS_26Aug14 ~ ., data = fss2.s1, controls = cforest_unbiased(ntree = 2000, mtry = mtry.fss2.s1))
  fss2.s1.vi[,i]<- varimp(fss2.s1.cf, conditional=FALSE)
}

test <- ctree(FSS_26Aug14 ~ ., data = fss2.s1)
plot(test)

# Add var names and index
fss2.s1.vi[,51]<- fss.col
fss2.s1.vi[,52]<-c(1:length(fss.col))
colnames(fss2.s1.vi)[51] <- "var_name"
colnames(fss2.s1.vi)[52] <- "var_index"

# ----------
# Save the df with a timestamp so we don't accidently overwrite it.
timestamp <- format(Sys.time(), "%Y%m%d_%H%M")
save(fss2.s1.vi, file=paste0("fss2_vi_s1_",timestamp,".RData"))

load("fss2_vi_s1_20141019_1451.RData")
# ----------

fss2.s1.vi.l <- melt(fss2.s1.vi, id=c("var_name","var_index"))

bymedian <- with(fss2.s1.vi.l, reorder(var_index, -value, median))
boxplot(value ~ bymedian, data = fss2.s1.vi.l,
        xlab = "Variable index", ylab = "Importance", 
        varwidth = TRUE,
        col = "lightgray")

fss2.s1.vi.median <- cast(fss2.s1.vi.l,var_name + var_index ~ ., value ='value', median)
colnames(fss2.s1.vi.median )[3] <- "median"

# sort the data so largest at the top
fss2.s1.vi.median <- fss2.s1.vi.median[with(fss2.s1.vi.median, order(-median)), ]

# Variable removal
# After the first 15 largest median values starts to flatten out so 
# we will take the top 30 to step 2.

# grab all variable names with median values > 0.5
fss2.s2.col <- fss2.s1.vi.median[fss2.s1.vi.median$median > 0.5,][,1]
fss2.s2.col <- c("FSS_26Aug14",fss2.s2.col)
fss2.s2 <- fss2.s1[,colnames(fss2.s1) %in% fss2.s2.col]

# remove any NAs
fss2.s2 <- data.frame(na.omit(fss2.s2))

#write.csv(fss2.s2, 'fss2_s2_data.csv')

# -----------------------------------------------------------
# Correlation plots
colnames(fss2.s2)

# Precip
pairs(fss2.s2[,c(2:6)],
      lower.panel=panel.smooth, upper.panel=panel.cor,diag.panel=panel.hist)

# Disturb
pairs(fss2.s2[,c(14:16,27)],
      lower.panel=panel.smooth, upper.panel=panel.cor,diag.panel=panel.hist)

# Lithology/soils
pairs(fss2.s2[,c(11:13,20:26)],
      lower.panel=panel.smooth, upper.panel=panel.cor,diag.panel=panel.hist)

# Ownership
pairs(fss2.s2[,c(17:19,28:31)],
      lower.panel=panel.smooth, upper.panel=panel.cor,diag.panel=panel.hist)

# Disturb and ownership
pairs(fss2.s2[,c(14:16,27,17:19,28:31)],
      lower.panel=panel.smooth, upper.panel=panel.cor,diag.panel=panel.hist)


# STEP 2. Here we follow reccomendations by Strobl et al (2008) and Strobl et al (2009) 
# and calculate conditional variable importance. This reduces importance scores on 
# variables that get hight socres from step 1 becuase they highly correlated to 
# to other ones. Probably applies to some of the climate variables.
# We couldn't do this in step 1 becuase we had too many variables and NAs.
# Not enough computer memory. A smaller dataset is ideal.



# initialize the variable importance df
fss2.s2.vi <- data.frame(matrix(, nrow = ncol(fss2.s2)-1, ncol = 1))
fss2.s2.col <- fss2.s2.col[!(fss2.s2.col == "FSS_26Aug14")]

# mtry and ntree values 
mtry.fss2.s2 <- as.integer(((ncol(fss2.s2)-1) / 3),0)

set.seed(55)
fss2.s2.cf <- cforest(FSS_26Aug14 ~ ., data = fss2.s2, controls = cforest_unbiased(mtry = mtry.fss2.s2, savesplitstats = TRUE))
#fss2.s2.vi[,1] <- varimp(fss2.s2.cf, conditional=TRUE, threshold = 0.8)
test2 <- varimp(fss2.s2.cf, conditional=FALSE)
