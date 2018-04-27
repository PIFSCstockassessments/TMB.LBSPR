#' @import parallel
#' @import data.table

get_F30 <- function(D,NumCores){

  #D <- Final  #THIS IS FOR TESTING PURPOSES, COMMENT-OUT WHEN DONE.

  D <- data.table(D)
  L <- split(D, seq(nrow(D)))

  #D <- L[[69]]  # FOR TESTING ONLY, COMMENT OUT

  FindF30 <- function(L){ # Function returns a F30 for a single line of data

    #L <- D # TESTING ONLY, COMMENT OUT.

    anF    <- (L$SPR/0.3)*L$Fmort  # Starting guess for F30
    SPR <- 0
    counter <- 0
    while(!(SPR>=0.29&SPR<=0.31)){ # This loop searches for F30

      SPR <- GetSPR(L,anF,anLc=NA)

      if(D$Fmort<0){anF=-9999;break}

      #print(c(round(counter,0),round(SPR,2),round(anF,5)))
      if(anF>=6 | counter>20 | anF<0){anF=-9999; break} # This break is to insure no infinite loop

      if(SPR<0.24)            {anF=anF*0.90}
      if(SPR>=0.24&SPR<0.29)  {anF=anF*0.98}
      if(SPR>=0.29&SPR<=0.31) {break;}
      if(SPR>0.31&SPR<0.36)   {anF=anF*1.02}
      if(SPR>=0.36&SPR<0.6)   {anF=anF*1.10}
      if(SPR>=0.6&SPR<0.7)    {anF=anF*1.50}
      if(SPR>=0.7)            {anF=anF*2.00}

      #if(abs(SPR-0.3)>0.15){ F_increment=anF*0.35}else if(abs(SPR-0.3)>0.15){F_increment=anF*0.2}else{F_increment=anF*0.1}
      #if(SPR>0.31){ anF = anF+F_increment}else if(SPR<0.29){anF=anF-F_increment}
      counter<-counter+1

    }

    return(anF)
  }

  # Execute parallel processing
  no_cores <- detectCores()-1
  #cl <- makeCluster(no_cores)
  cl <- makeCluster(NumCores)
  start<-proc.time()[3]
  clusterEvalQ(cl,require(TMB.LBSPR))
  Out <- parLapply(cl,L,FindF30)
  print((proc.time()[3]-start)/60)
  stopCluster(cl)
  #beep(sound=3);

  # Process data
  F30    <- vector(length=nrow(D))
  for(i in 1:length(F30)){
    F30[i] <- Out[[i]][[1]]
  }

  F30 <- data.table(F30)
  return(F30)

}
