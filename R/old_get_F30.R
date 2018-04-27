old_get_F30 <- function(D){ # This function loops through F values to find F30 for each iterations of the LB-SPR output

  D <- data.table(D)
  L <- split(D, seq(nrow(D)))

  F30Loop <- function(L){

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

    #anF         <- M #Starting point for estimating F30

    anF    <- (L$SPR/0.3)*L$Fmort/timesteps


    F_increment <- M*0.01

    SPR     <- 0
    counter <- 0
    Linf.vect <-  rnorm(n=20,mean=Linf,sd=Linf*CVLinf)
    while(!(SPR>=0.29&SPR<=0.31)){ # This loop searches for F30

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
          Selectivity[i]  <- 1/(1 + exp(-log(19)*(aLength-LS50)/(LS95-LS50) )  )
        }

        # Calculate exploited spawner biomass
        SBE=0
        N    <- vector(length=Amax)
        B    <- vector(length=Amax)
        N[1] <- 1000 # Recruitment set at 1000 recruits
        for (i in 2:Amax){
          N[i] <- N[i-1]*exp(-(M+anF*Selectivity[i]))
          B[i] <- N[i]*alpha*(aLinf*(1-exp(-K*(i-A0))))^beta # Updates N to B
          SBE  <- SBE+B[i]*Maturity[i]
        }
        GTG_SPR[GTG] <- SBE/SBP

      } # End of GTG loop

      SPR <- mean(GTG_SPR)

      if(anF*timesteps>=2.5 | counter>20 | anF<0){anF=-9999/timesteps; break} # This break is to insure no infinite loop
      if(abs(SPR-0.3)>0.15){ F_increment=anF*0.35}else if(abs(SPR-0.3)>0.15){F_increment=anF*0.2}else{F_increment=anF*0.1}
      if(SPR>0.31){ anF = anF+F_increment}else if(SPR<0.29){anF=anF-F_increment}
      counter<-counter+1
    } # End of while loop

    return(anF*timesteps)

  } # End of iteration loop function

  # Execute parallel processing
  no_cores <- detectCores()-1
  cl <- makeCluster(no_cores)
  start<-proc.time()[3]
  Out <- parLapply(cl,L,F30Loop)
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

} # End of GetF30 function
