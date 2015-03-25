####trying rf on bottom 2/3 of data####
fss2.s3.col <- fss2.s1.vi.median[fss2.s1.vi.median$median <= 5.477664e-01,][,1]
fss2.s3.col <- c("FSS_26Aug14",fss2.s3.col)
fss2.s3 <- fss2.s1[,colnames(fss2.s1) %in% fss2.s3.col]

# mtry and ntree values 
mtry.fss2.s3 <- as.integer(((ncol(fss2.s3)-1) / 3),0)

# initialize the variable importance df
fss2.s3.vi <- data.frame(matrix(, nrow = ncol(fss2.s3)-1, ncol = 50))
fss2.s3.visd <- data.frame(matrix(, nrow = ncol(fss2.s3)-1, ncol = 50))

fss2.s3.col <- colnames(fss2.s3)
fss2.s3.col <- fss2.s3.col[!(fss2.s3.col == "FSS_26Aug14")]

beg <- Sys.time()
set.seed(100)
for (i in 1:50) {
  fss2.s3.rf <- randomForest(FSS_26Aug14 ~ ., 
                             data = fss2.s3, 
                             ntree = 1000, 
                             keep.forest = TRUE, 
                             importance = TRUE)
  fss2.s3.vi[,i] <- importance(fss2.s3.rf, type = 1, conditional = TRUE)
  fss2.s3.visd[,i] <- fss2.s3.rf$importanceSD
}
print(Sys.time() - beg)

# Add var names and index
fss2.s3.vi[,51]<- fss2.s3.col
fss2.s3.vi[,52]<-c(1:length(fss2.s3.col))
colnames(fss2.s3.vi)[51] <- "var_name"
colnames(fss2.s3.vi)[52] <- "var_index"
fss2.s3.visd[,51]<- fss2.s3.col
fss2.s3.visd[,52]<-c(1:length(fss2.s3.col))
colnames(fss2.s3.visd)[51] <- "var_name"
colnames(fss2.s3.visd)[52] <- "var_index"

fss2.s3.vi.l <- melt(fss2.s3.vi, id=c("var_name","var_index"))

# bymedian <- sort(sapply(fss2.s2.rm.rf.vi, median))
# index.merge <- data.frame('variable' = names(bymedian), 'index' = 1:11)
# fm <- melt(fss2.s2.rm.rf.vi)
# fm <- merge(fm, index.merge, by = 'variable', all.x = TRUE)

#png('varImpALL_s2.png', width = 960, height = 960)
bymedian <- with(fss2.s3.vi.l, reorder(var_name, value, median))
par(yaxt="n",mar=c(5, 8, 4, 5))
boxplot(value ~ bymedian, data = fss2.s3.vi.l,
        xlab = "% Increase MSE", 
        varwidth = TRUE,
        col = "lightgray", horizontal = TRUE)
lablist.y<-levels(bymedian)
axis(2, labels = FALSE)
text(y = 1:63, par("usr")[1], labels = lablist.y, pos = 2, xpd = TRUE)

#R2
1-sum((fss2.s3$FSS_26Aug14-predict(fss2.s3.rf))^2)/sum((fss2.s3$FSS_26Aug14-mean(fss2.s3$FSS_26Aug14))^2)
#0.4067842
