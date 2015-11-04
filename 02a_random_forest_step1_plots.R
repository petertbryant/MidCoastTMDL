#png('varImpALL.png', width = 960, height = 960)
bymedian <- with(fss2.s1.vi.l, reorder(var_index, value, median))
boxplot(value ~ bymedian, data = fss2.s1.vi.l,
        ylab = "Variable index", xlab = "% Increase MSE", 
        varwidth = TRUE,
        col = "lightgray", horizontal = TRUE)
#dev.off()

#R2
1 - sum((fss2.s1$FSS_26Aug14 - predict(fss2.s1.rf))^2) / 
  sum((fss2.s1$FSS_26Aug14 - mean(fss2.s1$FSS_26Aug14))^2)
#0.5491