#' @import parallel
#' @import data.table

get_Lc30 <- function(D,NumCores){

  #D <- Final  #THIS IS FOR TESTING PURPOSES, COMMENT-OUT WHEN DONE.

  D <- data.table(D)
  L <- split(D, seq(nrow(D)))

  FindLc30 <- function(D){

    #D <- L[[564]]  # FOR TESTING ONLY, COMMENT OUT

    anLc <- 5
    SPR     <- 0
    counter <- 0
    Llambda <- D$Linf*(1-exp(-D$K*D$Amax)) # Calculate Llambda
    while(!(SPR>=0.28&SPR<=0.32)){ # This loop searches for F30

      SPR <- get_SPR(D,aFmort=NA,anLc=anLc)

      if(D$Fmort<0){anLc=-9999;break}

      #print(c(round(counter,0),round(SPR,2),round(anLc,5)))
      if(counter>20){anLc=-9999; break} # This break is to insure no infinite loop
      if(anLc<=0){anLc=0;break;}

      if(SPR<0.24)            {anLc=anLc+Llambda*0.1}
      if(SPR>=0.24&SPR<0.29)  {anLc=anLc+Llambda*0.02}
      if(SPR>=0.29&SPR<=0.31) {break;}
      if(SPR>0.31&SPR<0.36)   {anLc=anLc-Llambda*0.02}
      if(SPR>=0.36&SPR<0.6)   {anLc=anLc-Llambda*0.1}
      if(SPR>=0.6&SPR<0.7)    {anLc=anLc-Llambda*0.2}
      if(SPR>=0.7)            {anLc=anLc-Llambda*0.3}

      counter<-counter+1
    }

    return(anLc)
  }

  #FindLc30(L) # REMOVE, FOR TESTING ONLY

  # Execute parallel processing
  #no_cores <- detectCores()-1
  #cl <- makeCluster(no_cores)
  cl <- makeCluster(NumCores)
  start<-proc.time()[3]
  clusterEvalQ(cl,require(TMB.LBSPR))
  Out <- parLapply(cl,L,FindLc30)
  print((proc.time()[3]-start)/60)
  stopCluster(cl)
  #beep(sound=3);

  # Process data
  Lc30    <- vector(length=nrow(D))
  for(i in 1:length(Lc30)){
    Lc30[i] <- Out[[i]][[1]]
  }

  Lc30 <- data.table(Lc30)
  return(Lc30)

}
