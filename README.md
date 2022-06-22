# pyCHARMM Workshop

Date: June 20-24, 2022

Time: 4:00 - 5:30 PM EDT

Location: Chem 1706

----------------------------------------------

## pyCHARMM: A python instantiation of CHARMM

## Objectives

In this workshop, you will be introduced to basic setup and use of [pyCHARMM](https://charmm-dev.org/wiki/index.php/Pycharmm/Installing) for biomolecule modeling and simulation, including topics of ligand docking, ligand and protein preparation for MSλD simulations. The integration of other python functions, including machine learning derived energy functions will be discussed. Advanced topics including the use of MPI through python for parallel simulation control and integration of complex workflows involving integration of other python tools with CHARMM functionality.

The workshop attendees will utilize jupyter-lab to explore pyCHARMM basic tutorials. Bring your laptops to connect to the *gollum* cluster if you have an account or if you are remote use your local resources. Tutorial examples will not, in general, be compute intensive, but will utilize pymol, CHARMM/OpenMM and CHARMM/BLaDE.

----------------------------------------------

## Using Pre-built Environment with jupyter-lab on Gollum

If you have access to ***gollum***, you can use jupyter lab to run the tutorial examples to start a jupyter-server on ***gollum/satyr***.

**NOTE:** Since each jupyter server needs to bind to a dedicated port on gollum, you will need to choose a unique port number for your own jupyter server. It should be a dynamic port in the range from 49152 to 65535. (Some port numbers are already taken e.g. 65432)

To use a compute node for GPU computation with jupyter lab, you will need to login to your favorite compute node on gollum. It is always better to check whether the node is busy by checking the status with

```bash
# to check partition
sinfo -N -o '%10N %10T %18E %.8e %.8O %.15C' -p gpu    
# To check a single node
sinfo -N -o '%10N %10T %18E %.8e %.8O %.15C' -n gollum152 
```

For this example, I use ***gollum152***.

### On Gollum

```shell
# login to gollum 152 from head node
ssh gollum152 
intra_ip=`ifconfig | grep "inet 192.*" | awk '{print $2;}'`
module load pycharmm/0.4
jupyter lab --no-browser --ip $intra_ip --port xxxxx 
```

### On Local Machine

```shell
ssh -N -f -L localhost:xxxxx:gollum152:xxxxx gollum
```

Where `xxxxx` is the port number of your choice, and you should also substitute ***gollum152*** with the compute node of your choice (But of course you can also just use gollum152, and we can share the same compute node). You should be able to see an url like `http://127.0.0.1:xxxxx/lab` from the jupyter lab stdout. Copy and paste the url in a web browser, and you should be connected to jupyter lab on gollum with a pycharmm kernel (Do not use the url starting with 192). To check if you are on the correct node as well as the current GPU usage, you can type in a jupyter cell

```jupyter
%%bash
hostname
nvidia-smi
```

----------------------------------------------

## Preliminaries

For those who do not have access to *gollum* or would like to set up their own conda environment for pycharmm, clone this repo and create the conda environment by

```shell
cd pyCHARMM-Workshop/
conda create -n pyCHARMM -c conda-forge -c schrodinger \
pymol rdkit jupyterlab openbabel openmm
conda activate pyCHARMM
```

We will utilize pyCHARMM from the current developers stream. If you are working on *gollum* you can link to the pyCHARMM version I compiled using (tsch):

```shell
export CHARMM_HOME=/home/brookscl/charmm/c47-dev-release/
export CHARMM_LIB_DIR=$CHARMM_HOME/install-pycharmm-nompi/lib
```

Install pyCHARMM in the conda environment using pip

```shell
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

Following this you will need to set the environment variable `CHARMM_LIB_DIR` as noted above to point to your install library directory, and install the pyCHARMM modules using pip.

```shell
pip install ../tool/pycharmm
export CHARMM_LIB_DIR=/path/to/charmm_dir/install-pycharmm-nompi/lib
```

## Required modules/codes:

1. [**MMTSB Toolset**](https://github.com/mmtsb/toolset) is utilized and should be installed .
2. [**propKa**](https://github.com/jensengroup/propka) will be utilized.
3. [**pymol**](https://pymol.org/2/) will be utilzied locally to view structures.

## Using jupyter-lab locally

If you have pycharmm installed locally in the `pyCHARMM` virtual environment, you can start a jupyter-lab by typing `jupyter lab` in the same conda environment.

----------------------------------------------

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
