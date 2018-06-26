**Welcome to the Git Repository for the DAPPER (Data Assimilation for Predicting Productivity in Ecosystem and Regions) approach.**

The approach uses the language R to initiate and analyze a simulation and uses Fortran to execute a simulation.  We are assuming that you have R and a fortran compiler on your computer or cluster.

The repository has the latest updates.  If you are interested in the code used in Thomas et al. 2018 Ecological Applications, please download Tag 'Regional_Forecasting_Paper_2017'

Please use the following citations to reference the DAPPER approach

Thomas, R. Q., E. B. Brooks, , E. J. Ward, R. H. Wynne, T. J. Albaugh, H. Dinon-Aldridge, H. E. Burkhart, J.-C. Domec, T. R. Fox, C. A. Gonzalez-Benecke, T. A. Martin, A. Noormets, D. A. Sampson, and R. O. Teskey. 2017. Leveraging 35 years of Pinus taeda research in the southeastern US to constrain forest carbon cycle predictions: regional data assimilation using ecosystem experiments. Biogeosciences 14:3525-3547.

Thomas, R. Q., A. L. Jersild, E. B. Brooks, V. A. Thomas, R. H. Wynne. A mid-century ecological forecast with partitioned uncertainty predicts increases in loblolly pine forest productivity.  Accepted at Ecological Applications.

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

**Step 2: Prepare your input directory**

You need a separate directory to hold the input files.  The directory should contain a folder for each study and a folder for CO2.  The folder for each study must have three files
- Study_name_met.csv
- Study_name_obs.csv
- Study_name_plotlist.csv

The first is the meterology drivers, the second is the observations that are compared to the model predictions, and the third is characteristics of each plot including the initial conditions and soil information.

The *DAPPER_inputdata_public* repository (https://github.com/EcoDynForecast/DAPPER_inputdata_public) provides an example of the input directory using the plots at the Duke site (McCarthy et al. 2010 New Phytologist).

**Step 3: Create your working directory**

You need a directory that exists outside the the DAPPER directory that includes files that are unique to a particular run of DAPPER.  Currently the directory requires a priors file, a valdiation set file, and the run_DAPPER.R script. Examples of all required files can be found in the `DAPPER/example_run/` directory.

**Step 4: Change the paths in the `run_DAPPER.R` script**

The 'run_DAPPER.R` script runs the full analysis.  Change the following paths for the DAPPER code and DAPPER input data to match the paths on your computer

```{r}
working_directory = "/Users/quinn/Dropbox/Research/test_DAPPER_run/"
DAPPER_directory =  "/Users/quinn/Dropbox/Research/DAPPER/"
input_directory = "/Users/quinn/Dropbox/Research/DAPPER_inputdata_public"
```
**Step 5: Change your MCMC chain options**

Set the number of iterations that you want to run.  The total length of the MCMC chain will be the number of iterations x the number of fit parameters.

```{r}
niter = 20000
```

Set the number of iterations that you want to throw away as the burn in.  Value must be 1 or greater.

```{r}
burn = 10000
```

Set the interval between iterations in the chain taht you want to save. This will thin the chain and reduce the size of your saved chain

```{r}
thin_interval = 2
```

**Step 6: Set your priors file**

Your prior file is a csv file that is located in the working_directory folder

```{r}
priors_file = 'default_priors.csv'
```

**Step 7:  Name your run**

```{r}
run_name = 'SS_analysis'
```

**Step 8: Define starting point for MCMC chain**

If you set `restart_from_chain = FALSE` then the chain will start at the initial value set in the priors csv file.

If you set `restart_from_chain = TRUE` then the chain will start using the final iteration of a previously run chain:

```{r}
restart_chain =  'SS_val6.1.2017-09-03.08.22.08.Rdata'
```

Note:  the previously run chain has to have the same parameters and plots as the current chain for the restart to work without error.

**Step 9: Define set of plots used in the assimilation**

Change the `obs_set` variable to the integer corresponding to the set of plots that you want to use in the assimilation.  See lines 163 through 253 in the `DAPPER_directory/scripts/prepare_obs.R` script for a list of the supported sets of plots.  As an example:

```{r}
obs_set = 14 #This is for plots at the Duke Site
```

Note that `focal_plotID` allows you to assimilate a single plot but you have to know its ID in the input data.

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

**Step 10: Set other options**

`validation_set_file`: (NA or file name) If you want to hold a subset of plots out from the assimilation as validation, define the validation set file name here.

`create_plot`: (TRUE/FALSE) If TRUE, then a PDF analyzing the assimilation will be produced

`only_create_plot`: (TRUE/FALSE) If TRUE, then use the restart_chain to make PDF analyzing the chain without running a new assimilation

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

`use_gep` (0 or 1) use gross ecosystem productivity in likelihood

`use_et` (0 or 1) use evapotranspiration in likelihood

`use_ctrans` (0 or 1) use canopy transpiration in likelihood

`use_gep_uncert` (0 or 1) include observational uncertainity in GEP

`use_et_uncert` (0 or 1) include observational uncertainity in ET

`use_ctrans_uncert` (0 or 1) include observational uncertainity in Ctrans

**Step 11: Run run_DAPPER.R script**

Run the entire `run_DAPPER.R` script.  This will run all the other scripts and code.

**Step 12: Analyze assimilation**

Your chain will be located in the `working_directory` directory

A PDF with the chains, marginal histograms of the parameters, plot predictions with observations, and latent state for the focal plot is found in `working_directory`

**EXAMPLE**

As a example we provide the `run_DAPPER.R` script for assimilating observations from the Duke site.  This assimilation uses observations of GEP and ET from the Ameriflux database, biomass observations from McCarthy et al. 2010 (New Phytologist), and LAI produced by Eric Ward. 

Create a working directory (e.g., DAPPER_example_run) and copy the files out of the `/DAPPER/example_run/` directory to your working directory. Create the run_DAPPER.R script in your working directory and add the code below to the script.

```{r}
rm(list = ls())
#---CONTROL INFORMATION----------------------------
working_directory = "/Users/quinn/Dropbox/Research/DAPPER_example_run"
DAPPER_directory =  "/Users/quinn/Dropbox/Research/DAPPER/"
input_directory = "/Users/quinn/Dropbox/Research/DAPPER_inputdata_public"
niter = 100
chain_number = 1
burn =  1
thin_interval = 1
run_name = 'Your_duke_assimilation'
restart_from_chain = FALSE
restart_chain =  NA
priors_file = 'example_priors.csv'
validation_set_file = NA #file name of .csv that defines plot numbers for plots that are not fit but are compared to the obs.
create_plot = TRUE
only_create_plot = FALSE
obs_set = 14 #Select which plots are used in analysis.  See prepare_obs.R for number guide 
focal_plotID = NA #30001 #Setting a value here causes only a single plot to be simulated and fit
fr_model = 1  # 1 = estimate FR for each plot, 2 = empirical FR model
FR_fert_assumption = 0 #0 = assume fertilization plots have FR = 1, 1 = do not assume fertilization plots have FR = 1
FR_separate_npar_groups = 2  #Assigns a different parameter group to groups of FR values: 0 = all one group, 1 = separate groups, 2 = all plots separate
use_fol = TRUE  #TRUE= use allometric estimates of foliage biomass in fitting
use_dk_pars = 1  #0 = do not use 3 specific parameters for the Duke site, 1 = use the 3 specific parameters
use_age_edc = 0  #0 = do not use an ecological constraint on the age function (see code); 1 = use the constraint
use_sm_edc = 0  #0 = do not use an ecological constraint on the soil moisture function (see code); 1 = use the constraint
use_fr_edc = 0   #0 = do note use an ecological constraint on the SI - FR function (see code); 1 = use the constraint
nstreams = 19
use_gep = 1
use_et = 1  
use_ctrans =1 
use_gep_uncert = 1
use_et_uncert = 1
use_ctrans_uncert = 1
state_space = 1
tracked_plotnum = 1
windows_machine = FALSE
#----------------------------------------------------
all_studies = c(
  '/Duke/TIER4_Duke'
)

#---SELECT COMPONENTS THAT ARE ALLOWED TO HAVE UNCERTAINITY--
plot_WSx1000 = FALSE  #include plot specific WSx1000 parameter
plot_thinpower = FALSE #include plot specific thinpower parameter
plot_mort_rate = FALSE #include plot specific mortality rate parameter

setwd(paste(DAPPER_directory,'/scripts/',sep=''))
source('prepare_DAPPER_MCMC.R')
```
