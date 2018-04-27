RunModel <- function(D, Species, n_iteration,n_GTG,starting,ManageF30, ManageLc30, NumCores,Seed){
 
  require(TMB); require(data.table); require(parallel); require(truncnorm); require(MASS); require(ggplot2)
  
  set.seed(Seed)
  
  Results <- RunLBSPR(INP, Species, n_iteration, n_GTG, starting, NumCores)
  
  Final <- Results[[1]]
  
  #Final <- test[[1]]
  
  # Estimate biomass from catch
  Final$Bio.catch  <- Final$Catch/Final$Fmort/(Final$Fmort+Final$M)*(1-exp(-(Final$Fmort+Final$M)))
  Final$Bio.catch  <- ifelse(Final$Fmort<0,-9999,Final$Bio.catch)
  
  # Managment option calculations
  
  if(ManageF30==T)  { 
    F30      <- NewGetF30(Final,NumCores)
    Final    <- cbind(Final,F30)
    
    Final$C30.survey <- Final$Bio.survey*Final$F30/(Final$F30+Final$M)*(1-exp(-(Final$F30+Final$M)))
    Final$C30.catch  <- Final$Bio.catch*Final$F30/(Final$F30+Final$M)*(1-exp(-(Final$F30+Final$M)))  
    Final$F_Fmsy     <- Final$Fmort/Final$F30
    
    Final$C30.survey <- ifelse(Final$F30<0,-9999,Final$C30.survey)
    Final$C30.catch  <- ifelse(Final$F30<0|Final$Fmort<0,-9999,Final$C30.catch)
    Final$F_Fmsy     <- ifelse(Final$Fmort<0|Final$F30<0,-9999,Final$F_Fmsy)
    
  }
  
  if(ManageLc30==T) { 
    Lc30  <- NewGetLc30(Final,NumCores)
    Final <- cbind(Final,Lc30)
  }
  
  Results[[1]] <- Final
  GetModelFit(Results)
  
  return(Results)
  
}

GetModelFit <- function(Results){  # Residual graphs an preliminary results
  
  
  
  Final <- Results[[1]]
  
  # Create fitting and residual graphs
  Length_bins <- Results[[2]]
  Prop_obs    <- Results[[3]]
  Prop_exp    <- Results[[4]]
  
  #Exp_prop <- Exp_prop[rowSums(Exp_prop)<=1&rowSums(Exp_prop)>0.95,] # Removes anomalous outputs
  #Exp_prop[Exp_prop>0.5] <- NA # Removes anomalous outputs
  
  List     <- Final$Fmort>0&(rowSums(Prop_exp)<=1&rowSums(Prop_exp)>0.95) #Make a list of runs resulting in negative F (for filter) and remove anomalous outputs
  
  TotalCount <- INP[,11:length(INP)]
  TotalCount <- t(TotalCount)
  TotalCount <- apply(TotalCount,2,FUN=median,na.rm=TRUE)
  TotalCount <- sum(TotalCount)
  
  Prop_obs  <- Prop_obs[List,]
  Prop_obs  <- Prop_obs/rowSums(Prop_obs)
  Count_obs <- Prop_obs*TotalCount
  Count_obs <- Count_obs[]
  Final_count_obs <- apply(Count_obs,2,FUN=median,na.rm=TRUE)
  
  
  Prop_exp <- Prop_exp[List,]
  #Exp_prop <- Exp_prop[rowSums(Exp_prop)<=1&rowSums(Exp_prop)>0.95,] # Removes anomalous outputs
  #Exp_prop[Exp_prop>0.5] <- NA # Removes anomalous outputs
  Prop_exp        <- Prop_exp/rowSums(Prop_exp)
  Count_exp       <- Prop_exp*TotalCount
  Count_exp       <- Count_exp[]
  Final_count_exp <- apply(Count_exp,2,FUN=median,na.rm=TRUE)
  
  #Residual calculations
  resid        <- Count_obs-Count_exp
  resid.median <- apply(resid,2,FUN=median,na.rm=TRUE)
  resid.SD     <- apply(resid,2,FUN=sd,na.rm=TRUE)
  
  Data <- as.data.frame(cbind(Length_bins,Final_count_obs, Final_count_exp,resid.median,resid.SD))
  setnames(Data,c(1:5),c("Length_obs","Count_obs","Expected","Resid","Resid.SD"))
  
  # Graphs
  theme <- theme(axis.title=element_text(size=10),
                 axis.text=element_text(color="black"),
                 legend.text=element_text(size=7),
                 legend.title=element_text(size=7),
                 legend.margin=margin(unit(0,"cm")),
                 legend.position="none" )
  
  fit.plot   <- ggplot(data=Data)+geom_col(aes(x=Length_obs,y=Count_obs),fill="cadetblue3",col="black")+
                geom_line(aes(x=Length_obs,y=Expected),col="red",size=1.3)+
                scale_y_continuous(expand=c(0,0))+     
                labs(x="Observed fork length (mm)",y=bquote("Count"))+
                theme_classic()+theme
  
  
  res.plot <- ggplot(data=Data,aes(x=Length_obs))+geom_point(aes(y=resid.median))+
              geom_errorbar(aes(ymin=(resid.median-resid.SD),ymax=(resid.median+resid.SD)))+
              labs(x="Observed fork length (mm)",y=bquote("Residual (mm)"))+
              geom_hline(yintercept=0,col="red")+
              theme_classic()+theme
  
  Fgraph <- ggplot(data=Final,aes(x=Fmort))+geom_histogram()+theme_classic()+scale_x_continuous(limits=c(-1,2))+theme
  
  
  Final  <- subset(Final,Fmort>0)
  
  print(median(Final$Fmort))
  print(median(Final$SPR,na.rm=TRUE))
  print(median(Final$LS50))
  print(median(Final$LS95))
  print(median(Final$Bio.catch))
  print(median(Final$Bio.survey))
 
  mypath <- paste("1_OUTPUTS/","FIT",".tiff",sep="")
  ggsave(mypath,plot=fit.plot,units="cm",height=5.5,width=8.5,pointsize=5, dpi=300, compression="lzw")
  mypath <- paste("1_OUTPUTS/","RES",".tiff",sep="")
  ggsave(mypath,plot=res.plot,units="cm",height=5.5,width=8.5,pointsize=5, dpi=300, compression="lzw")
  
  print(fit.plot)
  print(Fgraph)
  
}

RunLBSPR <- function(D, Species, n_iteration,n_GTG,starting,NumCores){

# Specific parameters
D <- data.table(D)

#These 3 lines are for testing purposes. Gray out when not testing code.
#n_iteration <- 100
#n_GTG <- 20
#D <- data.table(  read.xlsx("C:\\Users\\Marc.Nadon\\Documents\\Work docs\\01_Projects\\01_Guam OFL\\0_R_Guam OFL\\Data\\LH PAR.xlsx",sheet="CAME")  )

LH.source         <- D$LH.source[1]
Family            <- D$Family[1]
Lmax              <- D[1,]$Val1
Lmax.sd           <- D[1,]$Val2
beta              <- D[10,]$Val1

n_iter_extra      <- n_iteration*3 

# Obtain life history parameters
if(LH.source=="Stepwise"){
  ParDist <- GetStepDist(Family, Lmax,Lmax.sd, n_iter_extra)
  ParDist$Lmat50 <- ParDist$Lmat-1
  ParDist$Lmat95 <- ParDist$Lmat
  ParDist$Amat   <- ParDist$A0-1/ParDist$K*log(1-ParDist$Lmat50/ParDist$Linf)
  
  # Calculate M obtained using S=0.05 to S=0.04
  M_004     <- -log(0.04)/ParDist$Amax
  ParDist$M <- M_004
}

if(LH.source=="Study"){
  
  ParDist        <- data.table(Linf=rep(1,n_iter_extra)) # Initialize
  ParDist$K      <- -9999
  ParDist$Lmax   <- -9999
  
  ParDist$Linf   <- RandomSample(n_iter_extra,D[2,]$Dist,D[2,]$Val1,D[2,]$Val2)
  ParDist$K      <- RandomSample(n_iter_extra,D[4,]$Dist,D[4,]$Val1,D[4,]$Val2)
  
  #cov        <- D[2,]$Val2*D[4,]$Val2*-0.66 # Inserting a -0.66 correlation coefficient between Linf and K
  
  #for(i in 1:n_iter_extra){
  #  repeat{
  #  GrowthDist <- mvrnorm(n=1,mu=rbind(D[2,]$Val1,D[4,]$Val1),Sigma=rbind(c(D[2,]$Val2,cov),c(cov,D[4,]$Val2)))
  #  if(GrowthDist[1]>0&GrowthDist[2]>0){break}
  #  }
    
  #  ParDist[i,1]    <- GrowthDist[1]
  #  ParDist[i,2]    <- GrowthDist[2]
  #}
  
  ParDist$A0     <- D[5,]$Val1
  
  if(!is.na(D[8,]$Val1)){
    ParDist$M    <- RandomSample(n_iter_extra,D[8,]$Dist,D[8,]$Val1,D[8,]$Val2)
    ParDist$Amax <- -log(0.04)/ParDist$M
  }
  
  if(is.na(D[8,]$Val1)){
    ParDist$Amax <- RandomSample(n_iter_extra,D[9,]$Dist,D[9,]$Val1,D[9,]$Val2)
    ParDist$M    <- -log(0.04)/ParDist$Amax
  }
  
  
  LmatDistance   <- D[7,]$Val1-D[6,]$Val1
  ParDist$Lmat50 <- RandomSample(n_iter_extra,D[6,]$Dist,D[6,]$Val1,D[6,]$Val2)
  ParDist$Lmat95 <- ParDist$Lmat50+LmatDistance
  ParDist$Amat   <- ParDist$A0-1/ParDist$K*log(1-ParDist$Lmat50/ParDist$Linf)
}

ParDist$CVLinf     <- RandomSample(n_iter_extra,D[3,]$Dist,D[3,]$Val1,D[3,]$Val2)
ParDist$beta       <- beta
ParDist$Bio.survey <- RandomSample(n_iter_extra,D[11,]$Dist,D[11,]$Val1,D[11,]$Val2)
ParDist$Catch      <- RandomSample(n_iter_extra,D[12,]$Dist,D[12,]$Val1,D[12,]$Val2)

# Remove problematic iterations and re-sample to get iteration count = n_iteration
FilteredParDist <- subset(ParDist,
                        Linf   >= D[2,]$Min & Linf   <= D[2,]$Max &
                        K      >= D[4,]$Min & K      <= D[4,]$Max &
                        Lmat95 >= D[7,]$Min & Lmat95 <= D[7,]$Max &
                        Amax   >= D[9,]$Min & Amax   <= D[9,]$Max &
                        Amat   < Amax & Amat > 0)

ParDist <- NULL
if(nrow(FilteredParDist)>=n_iteration){ ParDist <- FilteredParDist[1:n_iteration,] }else{
  n_missing    <- n_iteration-nrow(FilteredParDist)
  ParDist      <- FilteredParDist
  ExtraParDist <- FilteredParDist[sample(n_missing,replace=T),] 
  ParDist      <- rbind(FilteredParDist,ExtraParDist)
}

ParDist        <- ParDist[,c("Lmax","Linf","CVLinf","K","A0","beta","Lmat50","Lmat95","Amat","Amax","M","Bio.survey","Catch")]

# Calculate the size of the length bins used for the observed size structure
bin_diff <- vector(length=length(D$Length_obs)-1)
for(i in 2:length(D$Length_obs)){  bin_diff[i] <- D$Length_obs[i]-D$Length_obs[i-1]   }
counts <- table(bin_diff)
bin_size <- as.numeric(names(counts)[which.max(counts)])

# Fill missing zeroe counts in observed length data
Length_obs_full <-   D[,10:length(D)]
currentbin <- min(D$Length_obs)
numbin <- (max(D$Length_obs)-min(D$Length_obs))/bin_size
for(i in 1:(numbin)){
  
  if(Length_obs_full[i,1]==currentbin){currentbin <- currentbin+bin_size; next} # Skip
  
  newrow <- numeric(length(Length_obs_full))
  newrow[1] <- currentbin
  newrow <- data.frame(t(newrow))
  colnames(newrow) <- colnames(Length_obs_full)
  Length_obs_full <- rbind(Length_obs_full[1:i-1,],newrow,Length_obs_full[-(1:i-1),])
  
  currentbin <- currentbin+bin_size
  i <- i+1
}


# Create a list with each Monte Carlo inputs for the LBSPR model
Input    <- list(n_iteration)
for(i in 1:n_iteration){

  SD    <- round(ParDist$Linf[i]*ParDist$CVLinf[i],0)
  Min   <- round(ParDist$Linf[i]-1.645*SD,0)
  Max   <- round(ParDist$Linf[i]+1.645*SD,0)
  Range <- round(Max-Min,0)
  Incr  <- round(Range/(n_GTG-1),0)
  
  GTG_vec <- vector(length=n_GTG)
  GTG_vec[1] <- Min
  for(G in 2:n_GTG){ GTG_vec[G] <- round(GTG_vec[G-1]+Incr,0)       }
  n_l    <- max(GTG_vec)
  R0_vec <- dnorm(GTG_vec,ParDist$Linf[i],SD)

  # Draw a random size structure
  if(length(D)>10){
  aRandomSizeStruct <- sample(2:length(Length_obs_full),size=1, replace=T)
  aCount_obs        <- Length_obs_full[[aRandomSizeStruct]]
  } else { aCount_obs <- D[[10]] }
  
  aLength_obs <- Length_obs_full[[1]]
  
  Input[[i]] <- list(length_obs=aLength_obs,count_obs=aCount_obs,GTG_vec=GTG_vec,R0_vec=R0_vec,
                     M=ParDist$M[i],K=ParDist$K[i],Linf=ParDist$Linf[i],Mat50=ParDist$Lmat50[i],Mat95=ParDist$Lmat95[i],beta=beta,
                     CVLinf=ParDist$CVLinf[i],n_l=n_l,starting=starting)
}

RunGTG <- function(info){
  
  data       <- info
  
  #data <- Input[[1]] # FOR TESTING PURPOSES, DELETE AFTER
  #starting=list(Fmort=0.2, LS50=200, LS95=250) # SAME AS ABOVE
  #dyn.load(dynlib("GTG")) #SAME
  
  parameters <- starting
  
  # Run TMB model
  tryCatch({
  model <- MakeADFun(data, parameters, DLL="GTG")
  fit   <- nlminb(model$par, model$fn, model$gr)
  }, error = function(e) return ("Error!"))
    
  Results <- list(model$report()$LS50,model$report()$LS95,model$report()$Fmort,model$report()$SPR,model$report()$prop_obs,model$report()$prop_expect)
  return(Results)
}

# Execute parallel processing
no_cores <- detectCores()-1
cl <- makeCluster(NumCores)
clusterEvalQ(cl,require(TMB))
clusterEvalQ(cl,dyn.load(dynlib("GTG")))

start<-proc.time()[3]
Out <- parLapply(cl,Input,RunGTG)
print((proc.time()[3]-start)/60)
stopCluster(cl)
#beep(sound=3);

# Process data
Final    <- matrix(nrow=n_iteration, ncol=4); colnames(Final) <- c("LS50","LS95","Fmort","SPR")
Prop_exp <- matrix(ncol=length(aCount_obs),nrow=n_iteration)
Prop_obs <- matrix(ncol=length(aCount_obs),nrow=n_iteration)
for(i in 1:n_iteration){
  
  Final[i,1] <- Out[[i]][[1]] 
  Final[i,2] <- Out[[i]][[2]] 
  Final[i,3] <- Out[[i]][[3]] 
  Final[i,4] <- Out[[i]][[4]] 
  
  Prop_obs[i,] <- Out[[i]][[5]]
  Prop_exp[i,] <- Out[[i]][[6]]
}

Final    <- cbind(ParDist,Final)

aList <- list(Final,aLength_obs,Prop_obs,Prop_exp)

return(aList)

} # End of RunLBSPR function

NewGetF30 <- function(D,NumCores){
  
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
  
NewGetLc30 <- function(D,NumCores){
  
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
    
    SPR <- GetSPR(D,aFmort=NA,anLc=anLc)
    
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
  
GetSPR <- function(D,aFmort,anLc){ # Function to search F30 or Lc30
  
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

GetF30 <- function(D){ # This function loops through F values to find F30 for each iterations of the LB-SPR output
  
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

GetLc30 <- function(D){ # This function loops through F values to find F30 for each iterations of the LB-SPR output
  
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

