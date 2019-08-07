# TMB.LBSPR

Summary Description

## Pacakge Installation 



## Example 

This is how to run the TMB.LBSPR package using the "Example" dataset included in the package.

```
require(data.table) 
require(openxlsx) 
require(TMB.LBSPR)

INP      <- Example
aSpecies <- "TestSpecies"
```

To Run the model: 

```
# Run model
Results   <- run_analyses(INP,aSpecies,
                      n_iteration=1000,
                      n_GTG=20,
                      starting=list(Fmort=0.2, LS50=200, LS95=250),
                      ManageF30 =T,
                      ManageLc30=T,
                      NumCores=7,
                      Seed=1)
```

# Output

By default, ouptut from `run_analyses` will be written to `TMB.LBSPR` directory, located in the users HOME directory. 

For Windows users, the HOME directory is typically located at `C:/users/[USERNAME]`. Where `[USERNAME]` is the system's user account name.


## Github Disclaimer

This repository is a scientific product and is not official communication of the National Oceanic and Atmospheric Administration, or the United States Department of Commerce. All NOAA GitHub project code is provided on an ‘as is’ basis and the user assumes responsibility for its use. Any claims against the Department of Commerce or Department of Commerce bureaus stemming from the use of this GitHub project will be governed by all applicable Federal law. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or favoring by the Department of Commerce. The Department of Commerce seal and logo, or the seal and logo of a DOC bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by DOC or the United States Government.
