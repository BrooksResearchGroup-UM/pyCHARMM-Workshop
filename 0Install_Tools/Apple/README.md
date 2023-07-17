<div class="alert alert-block alert-info">
<h3><b>Note:</b> Apple hardware no longer supports CUDA, but does support OpenCL. Thus, we cannot build CHARMM/pyCHARMM versions that include BLaDe or DOMDEC GPU kernels. However, OpenMM supports OpenCL and thus can be installed and used with any GPU support on your Apple computer.</h3></div>    

## 0. Installing needed tools for CHARMM/pyCHARMM
### In order to use CHARMM/pyCHARMM you will need to:
- **Create a conda environment capable of building CHARMM, pyCHARMM (Part 1)**
- **Install gfortran and OpenMPI with MacPorts or Homebrew (Part 2)**
- **Install the MMTSB ToolSet from [MMTSB](https://feig.bch.msu.edu/mmtsb/Main_Page). Follow the instructions to install the package.**
- **Obtain the CHARMM software (free to academics and government labs) from [AcademicCHARMM](https://academiccharmm.org/program). Follow the directions below to build a conda environment capable of installing CHARMM/pyCHARMM.**
- **Install CHARMM and pyCHARMM (Part 3)**

## 1. Creating conda environment to install and use CHARMM/pyCHARMM
- **You will need a base anaconda/miniconda installation: see [anaconda installation](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html).**
- **Follow steps in 1a to create the conda environment manually OR follow steps in 1b to create the conda environment from a YAML file.**
### 1a. Create a conda environment manually
- **Make a conda environment (See below for a shortcut using `conda env create -f <name_of_environment>.yml`). In this example we will call the enviroment _charmm_env_**


`conda create -y -n charmm_env python=3.11` # note python can be > 3.9


- **Activate this environment**


`conda activate charmm_env`


- **Install mamba as a faster conda**


`conda install -y -c conda-forge mamba`
 


- **Install needed packages to build CHARMM and pyCHARMM**


`mamba install -y -c conda-forge make cmake fftw clfft openmpi openmm mpi4py readline==8.2 rdkit openbabel pandas pytorch jupyter_core jupyter_client jupyterlab jupyterlab_widgets jupyter_server jupyterlab_server jupyter_console jupyter jupytext biopython py3dmol mdtraj nglview jsonpickle pymol-open-source clfft pip`
    


### 1b. Building the CHARMM/pyCHARMM compatable environment with a YAML file
 
`charmm_env.yml`
 
```YAML
name: charmm_env # This represents the name you want to use for your conda environment
channels:        # This YAML file is for Apple hardware
  - conda-forge
  - defaults
dependencies:
  - python==3.11
  - mamba
  - make
  - cmake
  - fftw
  - openmpi
  - openmm
  - mpi4py
  - readline==8.2
  - rdkit
  - openbabel
  - pandas
  - pytorch
  - jupyter_core
  - jupyter_client
  - jupyterlab
  - jupyterlab_widgets
  - jupyter_server
  - jupyterlab_server
  - jupyter_console
  - jupyter
  - jupytext
  - biopython
  - py3dmol
  - mdtraj
  - nglview
  - jsonpickle
  - pymol-open-source
  - pip
  - ca-certificates
  - certifi
  - openssl
  - clfft
variables:  # Note one would need to change this to point to their path
  CHARMM_LIB_DIR: Users/brookscl/charmm/c47-dev-release/install_charmm_env/lib
prefix: /Users/brookscl/opt/anaconda3/envs/charmm_env
```

- **You can edit this YAML file to add/change the installed packages. You can install this new conda environment with the command:**


`conda env create -f charmm_env.yml`


## 2. Install gfortran and OpenMPI using MacPorts or Home.
### Uou have two choices, if you are already using either MacPorts or HomeBrew, skip down to installing gfortran and OpenMPI  below.
### MacPorts install:
> - **[Install](https://www.macports.org/install.php) MacPorts for your operating system Follow the directions at the link above.**
### `gfortran` and `OpenMPI` install
> `sudo port install gcc12`

> `sudo port install openmpi-gcc12`


### Homebrew install:
> - **[Install](https://brew.sh/) Homebrew for your operating system Follow the directions at the link above.**
### `gfortran` and `OpenMPI` install

> `brew install gcc`

> `brew install open-mpi`

<div class="alert alert-block alert-info">
<h3><b>Note:</b> If you want to use a particular gcc (current default seems to be gcc 13.1), then it appears you need the following instead.</h3>
</div>  

#### In bash

> `brew install gcc@12`

> `gfortran --version` # prints the version of gfortran

> `export HOMEBREW_CC=gcc-12.x` # use the same version as just queried from above

> `export HOMEBREW_CXX=g++-12.x`

> `brew install open-mpi --build-from-source`

## Once you've installed gcc and OpenMPI you can move onto building CHARMM and pyCHARMM (Part 3)


## 3. CHARMM and pyCHARMM installation once conda environment is installed and active and gfortran and OpenMPI are installed with MacPorts/Homebrew.
### Go to CHARMM source root and build CHARMM with configure

<blockquote>

```csh
conda activate charmm_env
cd <charmm_root>
mkdir build_charmm
cd build_charmm
# Build CHARMM with FFTDOCK, DOMDEC (default) and OpenMM (default)
../configure --without-cuda --with-fftdock -p <charmm_install_path>
make -j <n> install
```

</blockquote>

- **_charmm_env_ should be replaced with the name of your conda virtual environment**
- **\<charmm_root\> is the path to the charmm top level tree**
- **\<charmm_install_path\> is the path where you want the CHARMM installation to reside**
- **\<n\> is the number of cores to use in compiling the code**

### pyCHARMM is built from the same source and can be built in the same build directory

<blockquote>

```csh
conda activate charmm_wcuda12
cd <charmm_root>
cd build_charmm
rm -r *   # Clean the build directory
# Build CHARMM with FFTDOCK and OpenMM (default)
../configure --without-cuda --with-fftdock --without-mpi --as-library -p <pycharmm_install_path>
make -j <n> install
cd <charmm_root>
pip install `pwd`/tool/pycharmm  # Installs the pyCHARMM modules in your current environment
export CHARMM_LIB_DIR=<pycharmm_install_path>/lib # bash syntax
setenv CHARMM_LIB_DIR <pycharmm_install_path>/lib # csh syntax
conda env config vars set CHARMM_LIB_DIR=<pycharmm_install_path>/lib  # every time when this conda environment (charmm_env) is activated, the environmental variable CHARMM_LIB_DIR is there automatically.
```

</blockquote>

<div class="alert alert-block alert-info">
    <b>Note:</b> <i>$<$pycharmm_install_path$>$</i> is the path where you want the pyCHARMM installation to reside. <br><b>Note</b> that <i>$<$pycharmm_install_path$>$</i> can be the same as <i>$<$charmm_install_path$>$</i>, i.e., you can install both charmm and pyCHARMM in the same install folder.
</div>

```python

```
