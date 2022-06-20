# pyCHARMM Workshop
## June 20-24, 2022, 4:00 - 5:30 PM*, Chem 1706

# pyCHARMM: A python instantiation of CHARMM

## Objectives: In this workshop you will be introduced to basic setup and use of pyCHARMM for biomolecule modeling and simulation, including topics of ligand docking, ligand and protein preparation for MSλD simulations. The integration of other python functions, including machine learning derived energy functions will be discussed. Advanced topics including the use of MPI through python for parallel simulation control and integration of complex workflows involving integration of other python tools with CHARMM functionality. The workshop attendees will utilize jupyter-lab to explore pyCHARMM basic tutorials. Bring your laptops to connect to the gollum cluster if you have an account or if you are remote use your local resources. Tutorial examples will not, in general, be compute intensive, but will utilize pymol, CHARMM/OpenMM and CHARMM/BLaDE.

## Preliminaries: We will utilize pyCHARMM from the current developers stream. If you are working on gollum you can link to the pyCHARMM version I compiled using (tsch):
> _setenv CHARMM_LIB_DIR /users/brookscl/charmm/c47-dev-release/install-pycharmm-nompi/lib_
>  _pip install /users/brookscl/charmm/c47-dev-release/tool/pycharmm_

## Building pyCHARMM: If you are building pyCHARMM in your local environment (assuming Linux, gcc (I use 7.3), NVIDIA drivers (I use 10)/gpu access with OpenMM already installed), in the directory you want to build charmm (I use build_charmm):
> _../configure —with-blade —with-fftdock —without-mpi —as-library -p ../install-pycharmm-nompi_

## Following this you will need to set the environment variable __CHARMM_LIB_DIR__ as noted above to point to your install library directory, and install the pyCHARMM modules using pip as noted above.

## Required modules/codes: 
> 1. __MMTSB Toolset__ is utilized and should be installed (https://github.com/mmtsb/toolset).
> 2. __propKa__ (https://github.com/jensengroup/propka) will be utilized.
> 3. __pymol__ will be utilzied locally to view structures.

## Using jupyter-lab: You can use jupyter lab to run the tutorial examples (see pyCHARMMWorkshop.pdf in Notes).

# Workshop Agenda

## Monday 4:30 - 5:30 PM 
> __Introduction to CHARMM, how to build pyCHARMM and first steps in setting up molecular simulations with pyCHARMM, selecting atoms and basic manipulations of the PSF using pyCHARMM.__
## Tuesday 4:00 PM - 5:30 PM 
> __Setting up solvated molecular dynamics simulations, using OpenMM, BLaDE GPU accelerated simulation engines for molecular dynamics.__
## Wednesday 4:30 - 5:30 PM 
> __Using pyCHARMM CDOCKER for rigid-receptor, flexible ligand and flexible-receptor, flexible ligand docking.__
## Thursday 4:00 - 5:30 PM 
> __Free energy simulations i) absolute solvation free energy calculations; ii) setting-up MSλD with msld-py-prep and pyCHARMM.__
## Friday 4:00 - 5:30 PM 
> __Advanced topics including integration of python tools with pyCHARMM, utilizing machine-learned potentials in pyCHARMM, parallelism in python with MPI.__
