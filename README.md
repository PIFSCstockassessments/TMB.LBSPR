### This is how to run the TMB.LBSPR package using the "Example" dataset included in the package.

require(data.table); require(openxlsx); require(TMB.LBSPR)

INP      <- Example
aSpecies <- "TestSpecies"

# Run model
Results   <- run_analyses(INP,aSpecies,
                      n_iteration=1000,
                      n_GTG=20,
                      starting=list(Fmort=0.2, LS50=200, LS95=250),
                      ManageF30 =T,
                      ManageLc30=T,
                      NumCores=7,
                      Seed=1)
