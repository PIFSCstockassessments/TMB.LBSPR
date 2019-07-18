
#' Model Fit results
#'
#' Takes the model residuals from run_analyses and generates diagnostic plots
#'
#' @param Results TMB.LBSPR Results
#' @param INP Input Data
#' @param outdir Location to store outputs to file. Defaults to TMB.LBSPR directory in the system's Home.
#'
#' @import ggplot2
#' @importFrom grDevices dev.off tiff
#' @importFrom graphics hist
#' @importFrom stats IQR dnorm median nlminb quantile rnorm sd
#' @importFrom utils write.csv
#' @importFrom data.table data.table setnames
#'
model_fit <- function(Results, INP, outdir){  # Residual graphs an preliminary results

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

  fit.plot   <- ggplot(data=Data)+geom_col(aes_(x=~Length_obs,y=~Count_obs),fill="cadetblue3",col="black")+
    geom_line(aes_(x=~Length_obs,y=~Expected),col="red",size=1.3)+
    scale_y_continuous(expand=c(0,0))+
    labs(x="Observed fork length (mm)",y=bquote("Count"))+
    theme_classic()+theme


  res.plot <- ggplot(data=Data,aes_(x=~Length_obs))+geom_point(aes(y=resid.median))+
    geom_errorbar(aes(ymin=(resid.median-resid.SD),ymax=(resid.median+resid.SD)))+
    labs(x="Observed fork length (mm)",y=bquote("Residual (mm)"))+
    geom_hline(yintercept=0,col="red")+
    theme_classic()+theme

  Fgraph <- ggplot(data=Final,aes_(x=~Fmort))+geom_histogram()+theme_classic()+scale_x_continuous(limits=c(-1,2))+theme


  Final <- Final[Final$Fmort > 0,]

  print(median(Final$Fmort))
  print(median(Final$SPR,na.rm=TRUE))
  print(median(Final$LS50))
  print(median(Final$LS95))
  print(median(Final$Bio.catch))
  print(median(Final$Bio.survey))


  filename <- paste0("FIT",".tiff")
  mypath <- file.path(outdir,filename)
  ggsave(mypath,plot=fit.plot,units="cm",height=5.5,width=8.5,pointsize=5, dpi=300, compression="lzw")

  filename <- paste0("RES",".tiff")
  mypath <- file.path(outdir,filename)
  ggsave(mypath,plot=res.plot,units="cm",height=5.5,width=8.5,pointsize=5, dpi=300, compression="lzw")

  print(fit.plot)
  print(Fgraph)

}
