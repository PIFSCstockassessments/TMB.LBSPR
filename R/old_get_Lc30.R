old_get_Lc30 <- function(D){ # This function loops through F values to find F30 for each iterations of the LB-SPR output

  L <- data.table(D)
  L <- split(L, seq(nrow(L)))

  Lc30Loop <- function(L){

    timesteps <- 12 #Run model at different time steps per year

    Linf   <- L$Linf
    CVLinf <- L$CVLinf
    K      <- L$K/timesteps
    A0     <- L$A0*timesteps
    M      <- L$M/timesteps
    LS50   <- L$LS50
    LS95   <- L$LS95
    Lmat50 <- L$Lmat50
    Lmat95 <- L$Lmat95
    Amax   <- L$Amax*timesteps*1.5 #Run model longer than specified longevity
    beta   <- L$beta
    alpha  <- 1 # No need for the alpha in the L-W relationship
    Fvalue <- L$Fmort

    # Some guess starting values
    anLc <- (L$LS95+L$LS50)/2
    if(L$SPR>0.5){ anLc<- (L$LS95+L$LS50)/2*0.5 } else if (L$SPR<0.1) { anLc<- L$Linf - L$Linf*0.8   }
    if(L$SPR>0.2&L$SPR<0.4){ anLc <- (L$LS95+L$LS50)/2 }

    #print(paste("Start: ",anLc))

    Lc_increment <- anLc*0.01

    SPR     <- 0
    counter <- 0
    Linf.vect <-  rnorm(n=20,mean=Linf,sd=Linf*CVLinf)
    while(!(SPR>=0.28&SPR<=0.32)){ # This loop searches for F30

      if(L$SPR > 1 | L$Fmort<0) {anLC <- -9999; break}

      GTG_SPR <- vector(length=20) # Stores SPR for each GTG
      for(GTG in 1:20){  # This loop runs through 20 growth-type groups

        aLinf <- Linf.vect[GTG]
        # Calculate maturity vector
        Maturity <- vector(length=Amax)
        for(i in 1:Amax){
          aLength      <- Linf*(1-exp(-K*(i-A0)))
          Maturity[i]  <- 1/(1 + exp(-log(19)*(aLength-Lmat50)/(Lmat95-Lmat50) )  )
        }

        # Calculate pristine spawner biomass
        SBP=0
        N    <- vector(length=Amax)
        B    <- vector(length=Amax)
        N[1] <- 1000 # Recruitment set at 1000 recruits
        for (i in 2:Amax){
          N[i] <- N[i-1]*exp(-M)
          B[i] <- N[i]*alpha*(Linf*(1-exp(-K*(i-A0))))^beta # Updates N to B
          SBP  <- SBP+B[i]*Maturity[i]
        }

        # Calculate selectivity vector
        Selectivity <- vector(length=Amax)
        for(i in 1:Amax){
          aLength         <- aLinf*(1-exp(-K*(i-A0)))
          Selectivity[i]  <- 1/(1 + exp(-log(19)*(aLength-anLc+1)/(anLc+1-anLc) )  )
        }

        # Calculate exploited spawner biomass
        SBE=0
        N    <- vector(length=Amax)
        B    <- vector(length=Amax)
        N[1] <- 1000 # Recruitment set at 1000 recruits
        for (i in 2:Amax){
          N[i] <- N[i-1]*exp(-(M+Fvalue*Selectivity[i]))
          B[i] <- N[i]*alpha*(aLinf*(1-exp(-K*(i-A0))))^beta # Updates N to B
          SBE  <- SBE+B[i]*Maturity[i]
        }
        GTG_SPR[GTG] <- SBE/SBP

      } # End of GTG loop

      SPR <- mean(GTG_SPR)
      #print(paste("SPR: ",SPR))
      if(is.na(SPR)){ break }

      if(anLc>=Linf | counter>12){anLc=-9999; break} # This break is to insure no infinite loop
      if(abs(SPR-0.3)>0.20){ Lc_increment=anLc*0.15}else if(abs(SPR-0.3)>0.13){Lc_increment=anLc*0.1}else{Lc_increment=anLc*0.01}
      if(SPR>0.32){ anLc = anLc-Lc_increment}else if(SPR<0.28){anLc=anLc+Lc_increment}
      counter<-counter+1
    } # End of while loop

    return(anLc)

  } # End of iteration loop function

  # Execute parallel processing
  no_cores <- detectCores()-1
  cl <- makeCluster(no_cores)
  start<-proc.time()[3]
  Out <- parLapply(cl,L,Lc30Loop)
  #Out <- lapply(cl,L,Lc30Loop)
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

} # End of GetLc30 function

