#' Run TMB-LBSPR Analysis
#'
#' Runs TMB LBSPR model and calls process_results function.
#'
#' @param D Input Data
#' @param Species Species name
#' @param n_iteration Number of Iterations
#' @param n_GTG Number of growth-type groups used in simulation (default at 20)
#' @param starting List of starting values for the 3 estimated parameters
#' @param NumCores Number of CPU cores to use (Default?)
#' @param Seed Psuedo-Random Number Seed
#' @param ManageF30 Option to return max fishing mortality (F30) management recommendation
#' @param ManageLc30 Option to return minimum size (Lc30) management recommendation.
#'
#' @export

run_analyses <- function(D, Species, n_iteration,n_GTG,starting,ManageF30=TRUE, ManageLc30=TRUE, NumCores=-999,Seed=1){

  set.seed(Seed)

  Results <- run_lbspr(D, Species, n_iteration, n_GTG, starting, NumCores)

  Final <- Results[[1]]

  # Estimate biomass from catch
  Final$Bio.catch  <- Final$Catch/Final$Fmort/(Final$Fmort+Final$M)*(1-exp(-(Final$Fmort+Final$M)))
  Final$Bio.catch  <- ifelse(Final$Fmort<0,-9999,Final$Bio.catch)

  # Managment option calculations

  Final$F30 <- -9999
  if(ManageF30==T)  {
    Final$F30      <- get_F30(Final,NumCores)
    #Final    <- cbind(Final,F30)

    Final$C30.survey <- Final$Bio.survey*Final$F30/(Final$F30+Final$M)*(1-exp(-(Final$F30+Final$M)))
    Final$C30.catch  <- Final$Bio.catch*Final$F30/(Final$F30+Final$M)*(1-exp(-(Final$F30+Final$M)))
    Final$F_Fmsy     <- Final$Fmort/Final$F30

    Final$C30.survey <- ifelse(Final$F30<0,-9999,Final$C30.survey)
    Final$C30.catch  <- ifelse(Final$F30<0|Final$Fmort<0,-9999,Final$C30.catch)
    Final$F_Fmsy     <- ifelse(Final$Fmort<0|Final$F30<0,-9999,Final$F_Fmsy)

  }

  Final$Lc30 <- -9999
  if(ManageLc30==T) { Final$Lc30  <- get_Lc30(Final,NumCores)  }
  #Final <- cbind(Final,Lc30)

  # Figure out proper save path
  project_home <- file.path(home=Sys.getenv("HOME"), "TMB.LBSPR")
  if(!dir.exists(project_home)){
    message("Creating directory for TMB.LBSPR runs: ", project_home ," ...")
    dir.create(project_home)
  }

  outdir <- file.path(project_home,paste0("LBSPR_",Species))
  dir.create(outdir)

  Results[[1]] <- Final
  model_fit(Results, D, outdir)

if(ManageLc30==T&ManageF30==T){  # Skip graphics processing if management is turned off to prevent crash
  if(D$Val1[11]!=9999){
    process_outputs(Results[[1]],"Both",SHOW.LC=TRUE,outdir)
  }else{process_outputs(Results[[1]],"Catch only",SHOW.LC=TRUE,outdir)}
}
  return(Results)

}







