library(SSN)

options(scipen = 100)

#Get the model object
load('ssn1_glmssn_SSE.Rdata')

#get the SSN
ssn1 <- importSSN("//deqhq1/TMDL/TMDL_WR/MidCoast/Models/Sediment/SSN/LSN05/lsn.ssn", o.write = TRUE)

#Run a scenario with zero disturbance
ssn1.glmssn.SSE.0 <- ssn1.glmssn.SSE

preds.0 <- getSSNdata.frame(ssn1, Name = "preds")
