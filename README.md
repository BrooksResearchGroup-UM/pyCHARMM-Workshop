# pyCHARMM Workshop

Date: June 20-24, 2022

Time: 4:00 - 5:30 PM EDT

Location: Chem 1706

----------------------------------------------
## pyCHARMM: A python instantiation of CHARMM

## Objectives
In this workshop, you will be introduced to basic setup and use of [pyCHARMM](https://charmm-dev.org/wiki/index.php/Pycharmm/Installing) for biomolecule modeling and simulation, including topics of ligand docking, ligand and protein preparation for MSλD simulations. The integration of other python functions, including machine learning derived energy functions will be discussed. Advanced topics including the use of MPI through python for parallel simulation control and integration of complex workflows involving integration of other python tools with CHARMM functionality. 

The workshop attendees will utilize jupyter-lab to explore pyCHARMM basic tutorials. Bring your laptops to connect to the *gollum* cluster if you have an account or if you are remote use your local resources. Tutorial examples will not, in general, be compute intensive, but will utilize pymol, CHARMM/OpenMM and CHARMM/BLaDE.

-------------------------------------------------

## Preliminaries
Clone this repo and create conda environment for pycharmm

```shell
cd pyCHARMM-Workshop/
conda env create -f environment.yml -y
conda activate pyCHARMM
```

We will utilize pyCHARMM from the current developers stream. If you are working on *gollum* you can link to the pyCHARMM version I compiled using (tsch):
```shell
export CHARMM_HOME=/home/brookscl/charmm/c47-dev-release/
export CHARMM_LIB_DIR=$CHARMM_HOME/install-pycharmm-nompi/lib
```
Install pyCHARMM in the conda environment using pip
```
pip install $CHARMM_HOME/tool/pycharmm
```

## Building pyCHARMM 
If you are building pyCHARMM in your local environment (assuming Linux, gcc (I use v7.3), NVIDIA drivers (I use cuda v10.0)/gpu access with OpenMM already installed), in the **CHARMM directory**:

```shell
mkdir -p build_charmm; cd build_charmm
../configure --with-blade --with-fftdock --without-mpi \
--as-library -p ../install-pycharmm-nompi
make -j 2 -C install
```

Following this you will need to set the environment variable `CHARMM_LIB_DIR` as noted above to point to your install library directory, and install the pyCHARMM modules using pip as noted above.
```
pip install ../tool/pycharmm
export CHARMM_LIB_DIR=/path/to/charmm_dir/install-pycharmm-nompi/lib
```

## Required modules/codes: 
1. [**MMTSB Toolset**](https://github.com/mmtsb/toolset) is utilized and should be installed .
2. [**propKa**](https://github.com/jensengroup/propka) will be utilized.
3. [**pymol**](https://pymol.org/2/) will be utilzied locally to view structures.

## Using jupyter-lab
You can use jupyter lab to run the tutorial examples to start a jupyter-server on *gollum/satyr*. 

**NOTE:** Since each jupyter server needs to bind to a dedicated port on gollum, you will need to choose your own port number. It should be a Dynamic port in the range 49152 to 65535. (Some port numbers are already taken e.g. 65432)

**On Gollum**
```shell
ssh gollum
module load pycharmm/0.3
jupyter lab --no-browser --port xxxxx 
```
**On Local Machine**
```shell
ssh -N -f -L localhost:xxxxx:localhost:xxxxx gollum
```
Where `xxxxx` is the port number of your choice. You should be able to see a url like `http://localhost:xxxxx/lab` from the gollum stdout. Copy and paste the url in a web browser, and you should be connected to jupyter lab on gollum with pycharmm kernel.


------------------------------------------------------
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
