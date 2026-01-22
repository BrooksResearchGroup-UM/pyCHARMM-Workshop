<span style="color:red">
<h3><b>Note:</b> Apple hardware no longer supports CUDA, but does support OpenCL. Thus, we cannot build CHARMM/pyCHARMM versions that include BLaDe or DOMDEC GPU kernels. However, OpenMM supports OpenCL and thus can be installed and used with any GPU support on your Apple computer. At present we are working to make a fully gnu compiler supported version through a conda environment, but for present we still utilize <i>clang</i> and <i>clang++</i> from the MacOS Xcode toolchain.</h3></span>    

## 0. Installing needed tools for CHARMM/pyCHARMM
### In order to use CHARMM/pyCHARMM you will need to:
- **Create a conda environment capable of building CHARMM, pyCHARMM (Part 1)**
- **Install gfortran and OpenMPI with MacPorts or Homebrew (Part 2) (No longer necessary on MacOSX with Apple hardware)**
- **Install the MMTSB ToolSet from [MMTSB](https://feig.bch.msu.edu/mmtsb/Main_Page). Follow the instructions to install the package.**
- **Obtain the CHARMM software (free to academics and government labs) from [Academic CHARMM](https://academiccharmm.org/program). Follow the directions below to build a conda environment capable of installing CHARMM/pyCHARMM.**
- **Install CHARMM and pyCHARMM (Part 3)**

## 1. Creating conda environment to install and use CHARMM/pyCHARMM
- **You will need a base anaconda/miniconda/conda-forge (conda-forge recommended) installation: see [conda-forge installation](https://conda-forge.org/download/#:~:text=Miniforge%20is%20the%20preferred%20conda,conda%20create%20or%20mamba%20create%20.).**
- **Follow steps in 1 to create the conda environment from the included YAML file.**
> **Note: We have added new support for [OpenMM PyTorch plugin](https://github.com/openmm/openmm-torch).**

### 1. Building the CHARMM/pyCHARMM compatable environment with a YAML file
 
`charmm_cpu_env.yml`
 
```YAML
name: charmm_cpu_env # This represents the name you want to use for your conda environment
channels:        # This YAML file is for cpu-only implementation hardware
  - conda-forge
  - pytorch
  - defaults
dependencies:
  # =====================================================
  # Compilers (GCC for CHARMM compatibility)
  # =====================================================
  #- gcc  # gcc only for Linux at present
  #- gxx  # gcc only for Linux at present
  - gfortran
  - binutils
  - make
  - cmake
  # =====================================================
  # MPI (OpenMPI with CUDA awareness)
  # =====================================================
  - openmpi>=5.0.8
  - mpi4py
  # =====================================================
  # Python
  # =====================================================
  - python=3.12
  # =====================================================
  # OpenMM (includes CUDA platform when cuda-toolkit present)
  # =====================================================
  - openmm>=8.2
  - openmm-torch
  # =====================================================
  # Math libraries needed for CHARMM
  # =====================================================
  - fftw
  - clfft
  # =====================================================
  # TorchANI for support of ML  QM/ML potentials
  # =====================================================
  - torchani
  # =====================================================
  # Utilities usefull for CHARMM/pyCHARMM
  # =====================================================
  - readline
  - scipy
  - jsonpickle
  - pip
  - pdoc  # needed for documentation creation
  # =====================================================
  # Jupyter & Development
  # =====================================================
  - jupyter
  - ipywidgets
  - ipympl
  - ipyparallel
##########################################################
##########################################################
# The above is required or highly desired for basic CHARMM
# pyCHARMM build on a Linux cluster supporting CUDA
##########################################################
# File chm_addon.txt
##########################################################
# What's below this line is recommended for building
# modeling/computational workflows and their analyses
# with CHARMM/pyCHARMM
##########################################################
# =====================================================
# MD Analysis & Simulation Tools
# =====================================================
  - mdanalysis
  - mdtraj
  - biopython
  - parmed
  - pdbfixer
# =====================================================
# Free Energy & Enhanced Sampling
# =====================================================
  - pymbar
  - alchemlyb
# =====================================================
# Visualization
# =====================================================
  - matplotlib
  - nglview
  - py3dmol
  - pymol-open-source
# =====================================================
# Chemistry Tools
# =====================================================
  - rdkit
  - openbabel
  - mendeleev
# =====================================================
# pip Tools
# =====================================================
  - pip:
    - crimm # note latest crimm offers many structure preperation features
    - fastmbar
```

- **You can edit this YAML file to add/change the installed packages. You can install this new conda environment with the command:**


`conda env create -n <your_charmm_environment_name> -f charmm_cpu_env.yml`

- **If you do not specify `-n <your_charmm_environment_name>` it will default to the name at the top of the YAML file.**
- 
## 2. Install gfortran and OpenMPI using MacPorts or Homebrew.
### You have two choices, if you are already using either MacPorts or Homebrew, skip down to installing gfortran and OpenMPI  below. (Note thos is unncessary any longer since gfortran and openmpi are now available in conda-forge. 
### MacPorts install:
> - **[Install](https://www.macports.org/install.php) MacPorts for your operating system Follow the directions at the link above.**
### `gfortran` and `OpenMPI` install
> `sudo port install gcc13.4`

> `sudo port install openmpi-gcc13.4`

> `sudo port select --set gcc mp-gcc13.4`

> `sudo port select --set mpi openmpi-gcc13.4-fortran`

<div class="alert alert-block alert-warning">
## Homebrew install:
> - **[Install](https://brew.sh/) Homebrew for your operating system Follow the directions at the link above.**
### `gfortran` and `OpenMPI` install

> `brew install gcc`

> `brew install open-mpi`

<div class="alert alert-block alert-info">
<h4><b>Note:</b> If you want to use a particular gcc (current default seems to be gcc 13.1), then it appears you need the following instead.</h4>
</div>  

#### In bash

> `brew install gcc@12`

> `gfortran --version` # prints the version of gfortran

> `export HOMEBREW_CC=gcc-12.x` # use the same version as just queried from above

> `export HOMEBREW_CXX=g++-12.x`

> `brew install open-mpi --build-from-source`

## Once you've installed gcc and OpenMPI you can move onto building CHARMM and pyCHARMM (Part 3)


## 3. CHARMM and pyCHARMM installation once conda environment is installed and active and gfortran and OpenMPI are installed with MacPorts/Homebrew.
### Go to CHARMM source root and build CHARMM/pyCHARMM with configure

```csh
conda activate charmm_cpu_env
cd <charmm_root>
mkdir build_charmm
cd build_charmm
# Build CHARMM with FFTDOCK, DOMDEC (default) and OpenMM  and OpenMM-Torch (default)
../configure --without-cuda --with-fftdock -p <charmm_install_path>
make -j install
```

- **_charmm_cpu_env_ should be replaced with the name of your conda virtual environment**
- **\<charmm_root\> is the path to the charmm top level tree**
- **\<charmm_install_path\> is the path where you want the CHARMM installation to reside**
- **Note: CHARMM html documentation and pyCHARMM html files are created in the directory \<charmm_install_path\>/html_doc and \<charmm_install_path\>/html_doc\/pycharmm respectively.***

### CHARMM/pyCHARMM are built from the same source if you want to use pyCHARMM with mpi4py (the python api for mpi control at the python level use the following build

```csh
conda activate charmm_cpu_env
cd <charmm_root>
cd build_charmm
rm -r *   # Clean the build directory
# Build CHARMM with FFTDOCK and OpenMM (default)
../configure --without-cuda --with-fftdock --without-mpi -p <pycharmm_install_path>
make -j install
cd <charmm_root>
export CHARMM_LIB_DIR=<pycharmm_install_path>/lib # bash syntax
setenv CHARMM_LIB_DIR <pycharmm_install_path>/lib # csh syntax
conda env config vars set CHARMM_LIB_DIR=<pycharmm_install_path>/lib  # every time when this conda environment (charmm_env) is activated, the environmental variable CHARMM_LIB_DIR is there automatically.
```

## Note: *<pycharmm_install_path>* is the path where you want the pyCHARMM installation to reside. <br>
## Note that *<pycharmm_install_path>* can be the same as *<charmm_install_path>*, i.e., you can install both CHARMM and pyCHARMM in the same install folder *<charmm_install_path>*.<br>
## Note: *CHARMM_LIB_DIR* should automatically be derived from installation in an environment.

