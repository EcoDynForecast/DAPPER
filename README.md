Welcome to the Git Repository for the DAPPER (Data Assimilation for Predicting Productivity in Ecosystem and Regions) approach.

The approach uses the language R to initiate and analyze a simulation and uses Fortran to execute a simulation.  We are assuming that you have R and a fortran compiler on your computer or cluster.

**Step 1: Compile the fortran code**

Obtain a fortran compiler for your computer. We have tested the code for gfortran on Macs and the Linux cluster at Virginia Tech.  


At the command line use the appropriate commands below.  We encourage the use of openmp to parallelize the data assimilation.  The parallelization works by assigning different plots to different cores during an iteration of the MCMC.

```
### For MAC or LINUX: gfortran without openmp #### 

gfortran -c -fPIC prob_functions_mod.f90
gfortran -c -fPIC R3PG_MODEL_MOD.f90
gfortran -c -fPIC  DAPPER_plot_mod.f90
gfortran -shared  -fPIC -o DAPPER_MCMC.so DAPPER_MCMC_mod.f90 DAPPER_plot_mod.o R3PG_MODEL_MOD.o prob_functions_mod.o
gfortran -shared  -fPIC -o r3pg_interface.so R3PG_R_INTERFACE.f90 R3PG_MODEL_MOD.o

### For MAC or LINUX: gfortran with openmp #### 

gfortran -c -fPIC prob_functions_mod.f90 -fopenmp
gfortran -c -fPIC R3PG_MODEL_MOD.f90 -fopenmp
gfortran -c -fPIC  DAPPER_plot_mod.f90 -fopenmp
gfortran -shared  -fPIC -o DAPPER_MCMC.so DAPPER_MCMC_mod.f90 DAPPER_plot_mod.o R3PG_MODEL_MOD.o prob_functions_mod.o -fopenmp
gfortran -shared  -fPIC -o r3pg_interface.so R3PG_R_INTERFACE.f90 R3PG_MODEL_MOD.o -fopenmp

### For MAC or LINUX:  intel without openmp #### 

ifort -c -fPIC prob_functions_mod.f90
ifort -c -fPIC R3PG_MODEL_MOD.f90
ifort -c -fPIC  DAPPER_plot_mod.f90
ifort -shared  -fPIC -o DAPPER_MCMC.so DAPPER_MCMC_mod.f90 DAPPER_plot_mod.o R3PG_MODEL_MOD.o prob_functions_mod.o
ifort -shared  -fPIC -o r3pg_interface.so R3PG_R_INTERFACE.f90 R3PG_MODEL_MOD.o

### For MAC or LINUX:  intel with openmp #### 

ifort -c -fPIC prob_functions_mod.f90 -openmp
ifort -c -fPIC R3PG_MODEL_MOD.f90 -openmp
ifort -c -fPIC  DAPPER_plot_mod.f90 -openmp
ifort -shared  -fPIC -o DAPPER_MCMC.so DAPPER_MCMC_mod.f90 DAPPER_plot_mod.o R3PG_MODEL_MOD.o prob_functions_mod.o  -openmp
ifort -shared  -fPIC -o r3pg_interface.so R3PG_R_INTERFACE.f90 R3PG_MODEL_MOD.o -openmp

### For PC:  gfortran without openmp #### 

gfortran -c -fPIC prob_functions_mod.f90
gfortran -c -fPIC R3PG_MODEL_MOD.f90
gfortran -c -fPIC  DAPPER_plot_mod.f90
gfortran -shared  -fPIC -o DAPPER_MCMC.so DAPPER_MCMC_mod.f90 DAPPER_plot_mod.o R3PG_MODEL_MOD.o prob_functions_mod.o
gfortran -shared  -fPIC -o r3pg_interface.dll R3PG_R_INTERFACE.f90 R3PG_MODEL_MOD.o
```

**Step 2: Change the paths in the `run_DAPPER.R` script**

The 'run_DAPPER.R` script is runs the full analysis.  Change the following paths for the DAPPER code and DAPPER input data to match the paths on your computer

```{r}
working_directory =  '/Users/quinn/Dropbox/Research/DAPPER_papers/regional_forecasting/DAPPER'
input_directory = '/Users/quinn/Dropbox/Research/DAPPER_papers/regional_forecasting/DAPPER_inputdata'
```
**Step 3: Change your MCMC chain options**

Set the number of iterations that you want to run.  The total length of the MCMC chain will be the number of iterations x the number of fit parameters.

```{r}
niter = 20000
```

Set the number of iterations that you want to throw away as the burn in

```{r}
burn = 10000
```

Set the interval between save iterations in the chain.  This will thin the chain and reduce the size of your saved chain

```{r}
thin_interval = 2
```

**Step 4: Set your priors file**

Your prior file is a csv file that is located in the working_directory/priors folder

```{r}
priors_file = 'default_priors.csv'
```

**Step 5:  Name your run**

```{r}
run_name = 'SS_analysis'
```

**Step 6: Define starting point for MCMC chain**

If you set `restart_from_chain = FALSE` then the chain will start at the initial value set in the priors csv file.

If you set `restart_from_chain = TRUE` then the chain will start using the final iteration of a previously run chain: `restart_chain =  'SS_val6.1.2017-09-03.08.22.08.Rdata'`

Note:  the previously run chain has to have the same parameters and plots as the current chain for the restart to work without error

**Step 7: Define set of plots used in the assimilation**

Change the `obs_set` variable to the number corresponding to the plot set that you want to run.  See lines 163 through 253 in the `working_directory/scripts/prepare_obs.R` script for a list of the supported sets of plots.  As an example:

```{r}
obs_set = 14 #This is for plots at the Duke Site
obs_set = 21 #This is the full set of 294 plots
```

Note that `focal_plotID` allow you to assimilate a single plot but you have to know its ID in the input data.

Make sure that your plots are in the input_directory and the directory with in the input_directory is included in the `all_studies` vector (and not commented out).  For example the following vector would run the the 294 plots in Thomas et al. 2017:

```{r}
all_studies = c(
  '/SETRES/TIER4_SETRES',
  '/PINEMAP/TIER3_PINEMAP',
  '/NC2/TIER4_NC2',
  '/Duke/TIER4_Duke',
 '/Waycross/TIER4_Waycross',
  '/FMC_Thinning/TIER1_FMC_Thinning',
  '/FPC_RW18/TIER2_RW18'
)
```

**Step 8: Set other options**

`create_plot`: (TRUE/FALSE) If TRUE then a PDF analyzing the assimilation will be produced

`only_create_plot`: (TRUE/FALSE) If TRUE, then use the restart_chain to make PDF analyzing the assimilation with running a new assimilation

`val_set`: (integer) If you want to hold a subset of plots out from the assimilation as validation, define the validation set here.  The validation set is an integer that corresponds to a set of plots defined in the `/working_directory/validation_sets/validation_sets.csv` file.

`fr_model`: (1 or 2) 1 = estimate the plot-specific soil fertility parameter (FR or lambda) for each plot; 2 = use the empirical model to estimate the parameter for non-fertilized plots (currently a function of site index and mean annual temperature)

`FR_fert_assumption`: (0 or 1) 0 = assume that fertilization plots have an FR of 1; 1 = do not assume that fertilization plots have an FR of 1.

`FR_separate_npar_groups`: (TRUE/FALSE) This relates to the efficiency of the fitting by grouping plots where the FR values are simultaneously picked from the proposal distribution.

`use_fol`: (TRUE/FALSE) You may not think the foliage 'observations' should be used in the assimilation because they are derived from allometric relationships.  If so, set to FALSE so that they are not used in the assimilation.

`use_dk_pars`: (0 or 1) If 0 then do not use there separate parameters for the Duke Site (see Thomas et al. 2017 Biogeosciences), If 1 then use three separate parameters

`use_age_edc`, `use_sm_edc`, `use_fr_edc`: Don't worry about and don't use

`nstreams`: Do not change unless you increase the number of data streams

`state_space`: Do not change.  Should be equal to 1.

`tracked_plotnum`: Define the plotID of the plot that you want to the latent states

`windows_machine`: (TRUE/FALSE) set to TRUE if you are using a windows machine because the source code will have a .dll rather than a .so suffix.

`plot_WSx1000`,`plot_thinpower`,`plot_mort_rate`: Development stage, always set to FALSE

**Step 9: Run run_DAPPER.R script**

Run the entire `run_DAPPER.R` script.  This will run all the other scripts and code.

**Step 10: Analyze assimilation**

Your chain will be located in the `working_directory/chains/` directory

A PDF with the chains, marginal histograms of the parameters, plot predictions with observations, and latent state for the focal plot is found in `working_directory/figures/`
