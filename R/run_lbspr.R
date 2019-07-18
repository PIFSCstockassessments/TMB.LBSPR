#' @import parallel
#' @import TMB
#' @import StepwiseLH
#' @importFrom data.table data.table
#' @useDynLib TMB.LBSPR
#'
run_lbspr <- function(D, Species, n_iteration, n_GTG, starting, NumCores){

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
    ParDist <- Get_distributions(Family_Input=Family, Lmax.mean=Lmax,Lmax.SD=Lmax.sd,M_method=0.04, n_iter=n_iter_extra)
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

    ParDist$Linf   <- StepwiseLH::GenerateRandom(n_iter_extra,D[2,]$Dist,D[2,]$Val1,D[2,]$Val2)
    ParDist$K      <- StepwiseLH::GenerateRandom(n_iter_extra,D[4,]$Dist,D[4,]$Val1,D[4,]$Val2)

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
      ParDist$M    <- StepwiseLH::GenerateRandom(n_iter_extra,D[8,]$Dist,D[8,]$Val1,D[8,]$Val2)
      ParDist$Amax <- -log(0.04)/ParDist$M
    }

    if(is.na(D[8,]$Val1)){
      ParDist$Amax <- StepwiseLH::GenerateRandom(n_iter_extra,D[9,]$Dist,D[9,]$Val1,D[9,]$Val2)
      ParDist$M    <- -log(0.04)/ParDist$Amax
    }


    LmatDistance   <- D[7,]$Val1-D[6,]$Val1
    ParDist$Lmat50 <- StepwiseLH::GenerateRandom(n_iter_extra,D[6,]$Dist,D[6,]$Val1,D[6,]$Val2)
    ParDist$Lmat95 <- ParDist$Lmat50+LmatDistance
    ParDist$Amat   <- ParDist$A0-1/ParDist$K*log(1-ParDist$Lmat50/ParDist$Linf)
  }

  ParDist$CVLinf     <- StepwiseLH::GenerateRandom(n_iter_extra,D[3,]$Dist,D[3,]$Val1,D[3,]$Val2)
  ParDist$beta       <- beta
  ParDist$Bio.survey <- StepwiseLH::GenerateRandom(n_iter_extra,D[11,]$Dist,D[11,]$Val1,D[11,]$Val2)
  ParDist$Catch      <- StepwiseLH::GenerateRandom(n_iter_extra,D[12,]$Dist,D[12,]$Val1,D[12,]$Val2)

  # Initalize variables for used to in the data.table::subset function to deal with the
  # "undefined global functions or variables" issue when checking the R Package.
  Linf <- NULL
  K <- NULL
  Lmat95 <- NULL
  Amax <- NULL
  Amat <- NULL
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
      model <- MakeADFun(data, parameters, DLL="TMB.LBSPR")
      fit   <- nlminb(model$par, model$fn, model$gr)
    }, error = function(e) return ("Error!"))

    Results <- list(model$report()$LS50,model$report()$LS95,model$report()$Fmort,model$report()$SPR,model$report()$prop_obs,model$report()$prop_expect)
    return(Results)
  }

  # Execute parallel processing
  no_cores <- detectCores()-1
  if(NumCores<=0)    {cl <- makeCluster(no_cores)} # If the number of core is not manually specified, use total cores-1
  else if(NumCores>0){cl <- makeCluster(NumCores)}
  clusterEvalQ(cl,require(TMB))
  #clusterEvalQ(cl,dyn.load(dynlib("GTG")))
  #clusterEvalQ(cl,dyn.load(dynlib("TMB.LBSPR")))


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

