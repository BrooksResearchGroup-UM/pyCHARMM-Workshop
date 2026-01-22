# Preface
## This repository contains a working set of examples, current documentation and details of compiling and running pyCHARMM. The materials here are part of a growing repository for examples and tutorials on using pyCHARMM. As we evolve the program, expand its example base and run new workshops based on pyCHARMM, we anticiapte this repository will grow. 

## We note that pyCHARMM is now part of CHARMM distribution - _free of charge_ to academic and government/non-profit labs. The CHARMM/pyCHARMM source is available through academic licensing at [academiccharmm.org](https://academiccharmm.org/program).

# pyCHARMM Workshop

Date: June 20-24, 2022

Time: 4:00 - 5:30 PM EDT

Location: Chem 1706

## Preliminaries

For those who do not have access to *gollum* or would like to set up their own conda environment for pycharmm see the pages under **0Install_tools** for instructions on setting up the environment needed on your local machine.

We will utilize pyCHARMM from the current developers stream. If you are working on *gollum* you can link to the CHARMM/pyCHARMM version in the *gollum* module:

## Required modules/codes:
### In the workshop examples you will find below the following tools are used.

1. [**MMTSB Toolset**](https://github.com/mmtsb/toolset) is utilized and should be installed .

## Using jupyter-lab locally

If you have pycharmm installed locally in a virtual environment, you can start a jupyter-lab by typing `jupyter lab` in the same conda environment.

## If working through examples use the following order
1. **FirstExample** - build an alanine dipeptide and solvate it
2. **ProteinDynamics** - set-up solvated protein system and run dynamics on it
3. **SimpleMPIExample** - use mpi to compute alanine dipeptide energy surface and plot it
4. **AbsoluteSolvation** - Distributed (mpi) calculation of the solvent half of the absolute solvation free energy using MBAR and BAR
5. **Aladipeptide_HFBString_MPI** - Construct minimum energy path for $\phi,\psi$ conformational transition in alanine dipeptide using HFB-String approach and mpi task distribution
6. **ReplicaExchangeExample** - Compute the free energy change for charging a residue in vacuum with replica exchange $\lambda$-windowing distributed simulations using mpi
7. **CDOCKER_Tutorial** - Tutorial on using CDOCKER pyCHARMM implementation for flexible ligand - rigid receptor and flexible ligand - flexible receptor docking
8. **DrudeExample** - Simple example of setting-up and using Drude oscillator polarizable model
9. **NeuralNetExample** - Use of TorchANI machine learned quantum chemistry potentials in CHARMM calculations through pyCHARMM
10. **msld_py_prep_Tutorial** - Tutorial on setting-up and using MS $\lambda$ D free energy calculations with *py_prep*


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
