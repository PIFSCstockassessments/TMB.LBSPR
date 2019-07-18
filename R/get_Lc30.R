#' @import parallel
#' @importFrom data.table data.table

get_Lc30 <- function(D,NumCores){

  #D <- Final  #THIS IS FOR TESTING PURPOSES, COMMENT-OUT WHEN DONE.

  D <- data.table(D)
  L <- split(D, seq(nrow(D)))

  FindLc30 <- function(L){

    anLc <- 5
    SPR     <- 0
    counter <- 0
    Llambda <- L$Linf*(1-exp(-L$K*L$Amax)) # Calculate Llambda
    while(!(SPR>=0.28&SPR<=0.32)){ # This loop searches for F30

      SPR <- get_SPR(L,aFmort=NA,anLc=anLc)

      if(L$Fmort<0){anLc=-9999;break}

      #print(c(round(counter,0),round(SPR,2),round(anLc,5))) # COMMENT OUT
      #if(counter>25){anLc=-9999; break} # This break is to insure no infinite loop
      if(counter>25){anLc=0; break} # This break is to insure no infinite loop
      if(anLc<=0){anLc=0;break;}

      if(SPR<0.24)            {anLc=min(L$Llambda,anLc+Llambda*0.1)}
      if(SPR>=0.24&SPR<0.29)  {anLc=min(L$Llambda,anLc+Llambda*0.02)}
      if(SPR>=0.29&SPR<=0.31) {break;}
      if(SPR>0.31&SPR<0.36)   {anLc=max(1,anLc-Llambda*0.02)}
      if(SPR>=0.36&SPR<0.6)   {anLc=max(1,anLc-Llambda*0.1)}
      if(SPR>=0.6&SPR<0.7)    {anLc=max(1,anLc-Llambda*0.2)}
      if(SPR>=0.7)            {anLc=max(1,anLc-Llambda*0.3)}

      counter<-counter+1
    }

    return(anLc)
  }

  # Execute parallel processing
  #no_cores <- detectCores()-1
  #cl <- makeCluster(no_cores)
  cl <- makeCluster(NumCores)
  start<-proc.time()[3]
  clusterEvalQ(cl,require(TMB.LBSPR))
  Out <- parLapply(cl,L,FindLc30)
  #Out <- sapply(L,FindLc30)  # COMMENT OUT - TESTING
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
