# pyCHARMM Workshop
# June 20-24, 2022
# 4:00 - 5:30 PM*
# Chem 1706

# pyCHARMM: A python instantiation of CHARMM

Objectives: In this workshop you will be introduced to basic setup and use of pyCHARMM for biomolecule modeling and simulation, including topics of ligand docking, ligand and protein preparation for MSλD simulations. The integration of other python functions, including machine learning derived energy functions will be discussed. Advanced topics including the use of MPI through python for parallel simulation control and integration of complex workflows involving integration of other python tools with CHARMM functionality. The workshop attendees will utilize jupyter-lab to explore pyCHARMM basic tutorials. Bring your laptops to connect to the gollum cluster if you have an account or if you are remote use your local resources. Tutorial examples will not, in general, be compute intensive, but will utilize pymol, CHARMM/OpenMM and CHARMM/BLaDE.

Preliminaries: We will utilize pyCHARMM from the current developers stream. If you are working on gollum you can link to the pyCHARMM version I compiled using (tsch):
setenv CHARMM_LIB_DIR /users/brookscl/charmm/c47-dev-release/install-pycharmm-nompi/lib
pip install /users/brookscl/charmm/c47-dev-release/tool/pycharmm

Building pyCHARMM: If you are building pyCHARMM in your local environment (assuming Linux, gcc (I use 7.3), NVIDIA drivers (I use 10)/gpu access with OpenMM already installed), in the directory you want to build charmm (I use build_charmm):
../configure —with-blade —with-fftdock —without-mpi —as-library -p ../install-pycharmm-nompi

Following this you will need to set the environment variable CHARMM_LIB_DIR as noted above to point to your install library directory, and install the pyCHARMM modules using pip as noted above.

Using jupyter-lab: Truman from the group has provided information on running jupyter-lab through an interactive queue on gollum. If you are using gollum during the workshop, ask Truman. If you are using local resources you’ll need to configure jupyter-lab to use in your environment. Also, David Braun investigated how to use jupyter-lab remotely from your computer via the following procedure:
To run jupyter lab on the cluster and view it in a browser on your computer

1) ssh satyr
2) module load anaconda/2021 ...
3) conda activate <your virtual env>
(Note this environment should have pyCHARMM installed as well as any other python modules/applications you want to use.)
4) jupyter lab —no-browser
You should see the following somewhere in the output created:
[I 2022-06-20 09:31:36.674 ServerApp] Jupyter Server 1.17.1 is running at:
[I 2022-06-20 09:31:36.674 ServerApp] http://localhost:39004/lab
[I 2022-06-20 09:31:36.674 ServerApp]  or http://127.0.0.1:39004/lab
[I 2022-06-20 09:31:36.674 ServerApp] Use Control-C to stop this server and shut down all kernels (twice to skip confirmation).

5) on your home computer  run "ssh -N -f -L localhost:39004:localhost:39004 satyr" 
you will need to substitute the port your notebook is running on where I have 39004.
This step will only work if you use an ssh config file and 'ssh satyr' normally logs you in.

# Workshop Agenda

# Monday 4:30 - 5:30 PM 
Introduction to CHARMM, how to build pyCHARMM and first steps in setting up molecular simulations with pyCHARMM, selecting atoms and basic manipulations of the PSF using pyCHARMM.
# Tuesday 4:00 PM - 5:30 PM 
Setting up solvated molecular dynamics simulations, using OpenMM, BLaDE GPU accelerated simulation engines for molecular dynamics.
# Wednesday 4:30 - 5:30 PM 
Using pyCHARMM CDOCKER for rigid-receptor, flexible ligand and flexible-receptor, flexible ligand docking
# Thursday 4:00 - 5:30 PM 
Free energy simulations i) absolute solvation free energy calculations; ii) setting-up MSλD with msld-py-prep and pyCHARMM
# Friday 4:00 - 5:30 PM 
Advanced topics including integration of python tools with pyCHARMM, utilizing machine-learned potentials in pyCHARMM, parallelism in python with MPI.
