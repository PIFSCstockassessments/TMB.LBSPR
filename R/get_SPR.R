get_SPR <- function(D,aFmort,anLc){ # Function to search F30 or Lc30

  # Note: Divide by 10 to convert mm to cm to make this run faster.

  if(is.na(aFmort)){  # Search for Lc30 knowing Fmort
    LS50   <- (anLc/10)-1
    LS95   <- anLc/10
    Fmort  <- D$Fmort
  }

  if(is.na(anLc)){   # Search for F30 knowing selectivity
    LS50   <- D$LS50/10
    LS95   <- D$LS95/10
    Fmort  <- aFmort
  }

  n_GTG <- 20

  Linf   <- D$Linf/10
  CVLinf <- D$CVLinf
  K      <- D$K
  M      <- D$M

  Lmat50 <- D$Lmat50/10
  Lmat95 <- D$Lmat95/10
  beta   <- D$beta

  GTG_vec <- round(rnorm(n=20,Linf,Linf*CVLinf),0)
  n_l     <- max(GTG_vec)

  L_vec    <- vector(length=n_l)
  Sel_vec  <- vector(length=n_l)
  M_vec    <- vector(length=n_l)
  F_vec    <- vector(length=n_l)
  Z_vec    <- vector(length=n_l)
  NL_vec   <- vector(length=n_l)
  NL_prist_vec <- vector(length=n_l)
  Fec      <- vector(length=n_l)
  SSB_vec  <- vector(length=n_GTG)
  SSB0_vec <- vector(length=n_GTG)
  DL_vec   <- vector(length=n_l)
  DLS_mat  <- matrix(nrow=n_l,ncol=n_GTG)
  DL_final <- vector(length=n_l)
  DLS_final<- vector(length=n_l)

  # Populate the length vector
  L_vec[1]=1
  for(i in 2:n_l){
    L_vec[i]=L_vec[i-1]+1
  }

  # Calculate selectivity at length vector
  for(i in 1:n_l){
    Sel_vec[i]=1/(1+exp(-log(19)*(L_vec[i]-LS50)/(LS95-LS50)))
  }

  # Populate M vector
  for(i in 1:n_l){
    M_vec[i]=M
  }

  # Populate F vector
  for(i in 1:n_l){
    F_vec[i]=Fmort*Sel_vec[i]
  }

  # Populate Z vector
  for(i in 1:n_l){
    Z_vec[i]=M_vec[i]+F_vec[i];
  }


  # Start of GTG loop
  for(GTG in 1:n_GTG){

    #GTG <- 1  # REMOVE THIS WHEN DONE

    aLinf=GTG_vec[GTG]

    # Populate Number at length vector
    NL_vec[1]=1*((aLinf-L_vec[1]-1)/(aLinf-L_vec[1]))^(Z_vec[1]/K)
    for(i in 2:n_l){
      NL_vec[i]<-NL_vec[i-1]*((aLinf-L_vec[i]-1)/(aLinf-L_vec[i]))^(Z_vec[i]/K)
      if(L_vec[i]>=(aLinf-1)){NL_vec[i]<-0}
    }

    # Populate pristine Number at length vector
    NL_prist_vec[1]=1*((aLinf-L_vec[1]-1)/(aLinf-L_vec[1]))^(M_vec[1]/K)
    for(i in 2:n_l){
      NL_prist_vec[i]<-NL_prist_vec[i-1]*((aLinf-L_vec[i]-1)/(aLinf-L_vec[i]))^(M_vec[i]/K)
      if(L_vec[i]>=(aLinf-1)){NL_prist_vec[i]<-0;}
    }

    # Calculate fecundity at length vector
    for(i in 1:n_l-1){
      Fec[i] <- 1/(1+exp(-log(19)*(L_vec[i]-Lmat50)/(Lmat95-Lmat50)))*L_vec[i]^beta
    }

    # Calculate spawning per recruit
    aSSB <- 0
    for(i in 1:(n_l-1)){
      aSSB <- aSSB+1/Z_vec[i]*(NL_vec[i]-NL_vec[i+1])*Fec[i]
    }
    SSB_vec[GTG]=aSSB

    # Calculate pristine spawning per recruit
    aSSB0 <- 0
    for(i in 1:(n_l-1)){
      aSSB0<-aSSB0+1/M_vec[i]*(NL_prist_vec[i]-NL_prist_vec[i+1])*Fec[i];
    }
    SSB0_vec[GTG]=aSSB0

    # Populate Density at length vector
    for(i in 1:(n_l-1)){
      DL_vec[i]<-1/Z_vec[i]*(NL_vec[i]-NL_vec[i+1]);
    }

    # Standardize DL_vec column
    ColSum<-0;
    for(i in 1:(n_l-1)){
      ColSum<-ColSum+DL_vec[i];
    }

    for(i in 1:(n_l-1)){
      DLS_mat[i,GTG]<-DL_vec[i]/ColSum;
    }

  } # End of GTG loop

  # Sum density by length for all GTGs
  for(i in 1:n_l){
    DL_final[i]=0;
    for(GTG in 1:n_GTG){
      DL_final[i] <- DL_final[i]+DLS_mat[i,GTG];
    }
    DL_final[i] <- DL_final[i]*Sel_vec[i];
  }

  # Standardize DL_final column
  ColSum <- 0
  for(i in 1:n_l){
    ColSum=ColSum+DL_final[i];
  }

  for(i in 1:n_l){
    DLS_final[i]<-DL_final[i]/ColSum
  }

  # Calculate SPR
  SSB<-0
  SSB0<-0
  for(GTG in 1:n_GTG){
    SSB=SSB+SSB_vec[GTG];
    SSB0=SSB0+SSB0_vec[GTG];
  }
  SPR <- SSB/SSB0

  return(SPR)

} # End of GetSPR function
