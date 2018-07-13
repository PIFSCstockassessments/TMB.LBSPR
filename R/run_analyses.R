#' @import data.table
#' @export
run_analyses <- function(D, Species, n_iteration,n_GTG,starting,ManageF30, ManageLc30, NumCores,Seed){

  set.seed(Seed)

  Results <- run_lbspr(INP, Species, n_iteration, n_GTG, starting, NumCores)

  Final <- Results[[1]]

  # Estimate biomass from catch
  Final$Bio.catch  <- Final$Catch/Final$Fmort/(Final$Fmort+Final$M)*(1-exp(-(Final$Fmort+Final$M)))
  Final$Bio.catch  <- ifelse(Final$Fmort<0,-9999,Final$Bio.catch)

  # Managment option calculations

  if(ManageF30==T)  {
    F30      <- get_F30(Final,NumCores)
    Final    <- cbind(Final,F30)

    Final$C30.survey <- Final$Bio.survey*Final$F30/(Final$F30+Final$M)*(1-exp(-(Final$F30+Final$M)))
    Final$C30.catch  <- Final$Bio.catch*Final$F30/(Final$F30+Final$M)*(1-exp(-(Final$F30+Final$M)))
    Final$F_Fmsy     <- Final$Fmort/Final$F30

    Final$C30.survey <- ifelse(Final$F30<0,-9999,Final$C30.survey)
    Final$C30.catch  <- ifelse(Final$F30<0|Final$Fmort<0,-9999,Final$C30.catch)
    Final$F_Fmsy     <- ifelse(Final$Fmort<0|Final$F30<0,-9999,Final$F_Fmsy)

  }

  if(ManageLc30==T) {
    Lc30  <- get_Lc30(Final,NumCores)
    Final <- cbind(Final,Lc30)
  }


  # Figure out proper save path
  project_home <- file.path(home=Sys.getenv("HOME"), "TMB.LBSPR")
  if(!dir.exists(project_home)){
    message("Creating directory for TMB.LBSPR runs: ", project_home ," ...")
    dir.create(project_home)
  }

  outdir <- file.path(project_home,paste0("LBSPR_",Species))
  dir.create(outdir)

  Results[[1]] <- Final
  model_fit(Results, outdir)

  if(D$Val1[11]!=9999){
    process_outputs(Results[[1]],"Both",SHOW.LC=TRUE,outdir)
  }else{process_outputs(Results[[1]],"Catch only",SHOW.LC=TRUE,outdir)}

  return(Results)

}







