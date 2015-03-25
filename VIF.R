#stepwise VIF function used below
vif_func<-function(in_frame,ob.ssn){
  
  vif_init<-NULL
  for(val in names(in_frame)){
    form_in<-formula(paste(val,' ~ ',paste(setdiff(names(in_frame),val),collapse = ' + ')))
    tmp.ssn <- glmssn(form_in, 
                      ob.ssn,
                      EstMeth = "ML",
                      CorModels = c("locID","Exponential.tailup", "Exponential.taildown", "Exponential.Euclid"),
                      addfunccol = "afvArea",
                      family = "Gaussian")
    ssn.summary <- summary(tmp.ssn)
    VIF.val <- 1/(1 - ssn.summary$Rsquared)
    new.row <- c(val, VIF.val)
    vif_init <- rbind(vif_init,new.row)
  }
  return(vif_init)
}
#   
#   for(val in names(in_frame)){
#     form_in<-formula(paste(val,' ~ .'))
#     vif_init<-rbind(vif_init,c(val,VIF(lm(form_in,data=in_frame))))
#   }
#   vif_max<-max(as.numeric(vif_init[,2]))
#   
#   if(vif_max < thresh){
#     if(trace==T){ #print output of each iteration
#       prmatrix(vif_init,collab=c('var','vif'),rowlab=rep('',nrow(vif_init)),quote=F)
#       cat('\n')
#       cat(paste('All variables have VIF < ', thresh,', max VIF ',round(vif_max,2), sep=''),'\n\n')
#     }
#     return(names(in_frame))
#   }
#   else{
#     
#     in_dat<-in_frame
#     
#     #backwards selection of explanatory variables, stops when all VIF values are below 'thresh'
#     while(vif_max >= thresh){
#       
#       vif_vals<-NULL
#       
#       for(val in names(in_dat)){
#         form_in<-formula(paste(val,' ~ .'))
#         vif_add<-VIF(lm(form_in,data=in_dat))
#         vif_vals<-rbind(vif_vals,c(val,vif_add))
#       }
#       max_row<-which(vif_vals[,2] == max(as.numeric(vif_vals[,2])))[1]
#       
#       vif_max<-as.numeric(vif_vals[max_row,2])
#       
#       if(vif_max<thresh) break
#       
#       if(trace==T){ #print output of each iteration
#         prmatrix(vif_vals,collab=c('var','vif'),rowlab=rep('',nrow(vif_vals)),quote=F)
#         cat('\n')
#         cat('removed: ',vif_vals[max_row,1],vif_max,'\n\n')
#         flush.console()
#       }
#       
#       in_dat<-in_dat[,!names(in_dat) %in% vif_vals[max_row,1]]
#       
#     }
#     
#     return(names(in_dat))
#     
#   }
#   
# }