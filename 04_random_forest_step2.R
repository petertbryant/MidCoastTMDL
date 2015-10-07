# --- #
##### Random Forest Step 2 ####
colnames(fss2.s2)

# mtry and ntree values 
mtry.fss2.s2 <- as.integer(((ncol(fss2.s2)-1) / 3),0)

# initialize the variable importance df
fss2.s2.vi <- data.frame(matrix(nrow = ncol(fss2.s2) - 1, ncol = 50))
fss2.s2.visd <- data.frame(matrix(nrow = ncol(fss2.s2) - 1, ncol = 50))

fss2.s2.col <- colnames(fss2.s2)
fss2.s2.col <- fss2.s2.col[!(fss2.s2.col == "FSS_26Aug14")]

beg <- Sys.time()
set.seed(100)
for (i in 1:50) {
  fss2.s2.rf <- randomForest(FSS_26Aug14 ~ ., 
                             data = fss2.s2, 
                             ntree = 1000, 
                             keep.forest = TRUE, 
                             importance = TRUE)
  fss2.s2.vi[, i] <- importance(fss2.s2.rf, type = 1, conditional = TRUE)
  fss2.s2.visd[, i] <- fss2.s2.rf$importanceSD
}
print(Sys.time() - beg)

# Add var names and index
fss2.s2.vi[, 51]<- fss2.s2.col
fss2.s2.vi[, 52]<-c(1:length(fss2.s2.col))
colnames(fss2.s2.vi)[51] <- "var_name"
colnames(fss2.s2.vi)[52] <- "var_index"
fss2.s2.visd[, 51]<- fss2.s2.col
fss2.s2.visd[, 52]<-c(1:length(fss2.s2.col))
colnames(fss2.s2.visd)[51] <- "var_name"
colnames(fss2.s2.visd)[52] <- "var_index"

# --- #
# Save the df with a timestamp so we don't accidently overwrite it.
timestamp <- format(Sys.time(), "%Y%m%d_%H%M")
save(fss2.s2.vi, file = paste0("fss2_s2_vi_", timestamp, ".RData"))
save(fss2.s2.visd, file = paste0("fss2_s2_visd_", timestamp, ".RData"))
timestamp
load("fss2_s2_vi_20141027_2010.RData")
load("fss2_s2_visd_20141027_2010.RData")

#### s2 boxplot ####

fss2.s2.vi.l <- melt(fss2.s2.vi, id = c("var_name", "var_index"))

# bymedian <- sort(sapply(fss2.s2.rm.rf.vi, median))
# index.merge <- data.frame('variable' = names(bymedian), 'index' = 1:11)
# fm <- melt(fss2.s2.rm.rf.vi)
# fm <- merge(fm, index.merge, by = 'variable', all.x = TRUE)

fss2.vars <- fss2.s1.vi[fss2.s1.vi$var_name %in% setdiff(names(fss2.s2),
                                                         "FSS_26Aug14"),]
fss2.s2.vi.l <- melt(fss2.vars, id = c("var_name", "var_index"))
png('varImp_s2.png', width = 960, height = 960)
bymedian <- with(fss2.s2.vi.l, reorder(var_name, value, median))
par(yaxt = "n",mar = c(5, 8, 4, 5), cex = 2)
boxplot(value ~ bymedian, data = fss2.s2.vi.l,
        xlab = "% Increase MSE", 
        varwidth = TRUE,
        col = "lightgray", 
        horizontal = TRUE)
lablist.y<-levels(bymedian)
axis(2, labels = FALSE)
text(y = 1:100, par("usr")[1], labels = lablist.y, pos = 2, xpd = TRUE)
dev.off()

png('varImpALL_s2.png', width = 960, height = 960)
bymedian <- with(fss2.s2.vi.l, reorder(var_name, value, median))
par(yaxt="n", mar=c(5, 8, 4, 5), cex=2)
boxplot(value ~ bymedian, data = fss2.s2.vi.l,
        xlab = "% Increase MSE", 
        varwidth = TRUE,
        col = "lightgray", 
        horizontal = TRUE)
lablist.y<-levels(bymedian)
axis(2, labels = FALSE)
text(y = 1:14, par("usr")[1], labels = lablist.y, pos = 2, xpd = TRUE)
dev.off()

#R2
1 - sum((fss2.s2$FSS_26Aug14 - predict(fss2.s2.rf))^2) / 
  sum((fss2.s2$FSS_26Aug14 - mean(fss2.s2$FSS_26Aug14))^2)
#0.5658064
#This corroborates the use of these variables by showing that with a smaller subest of variables
#we achieve essentially the same R2. When we look at the bottom 2/3 of the variables we get a much smaller R2.

#### Median df creation ####

fss2.s2.vi.median <- cast(fss2.s2.vi.l,var_name + var_index ~ ., 
                          value ='value', median)
colnames(fss2.s2.vi.median )[3] <- "median"

# sort the data so largest at the top
fss2.s2.vi.median <- fss2.s2.vi.median[with(fss2.s2.vi.median, 
                                            order(-median)), ]

#### Partial Dependence Plots ####
for (i in 1:length(fss2.s2.vi$var_name)) {
  filename <- paste("partialPlot_", fss2.s2.vi$var_name[i], ".png", sep = "")
  png(filename, width = 960, height = 960)
  partialPlot(fss2.s2.rf, 
              fss2.s2, 
              x.var = fss2.s2.vi$var_name[i], 
              ylab = 'Mean FSS', 
              ylim = c(9, 18))
  dev.off()  
}
