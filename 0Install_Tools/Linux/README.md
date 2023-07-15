## 0. Installing needed tools for CHARMM/pyCHARMM
### In order to use CHARMM, pyCHARMM, and pyALF you will need to:
- **Create a conda environment capable of building CHARMM, pyCHARMM, and pyALF (Part 1)**
- **Install the MMTSB ToolSet from [MMTSB](https://feig.bch.msu.edu/mmtsb/Main_Page). Follow the instructions to install the package.**
- **Obtain the CHARMM software (free to academics and government labs) from [AcademicCHARMM](https://academiccharmm.org/program). Follow the directions below to build a conda environment capable of installing CHARMM/pyCHARMM.**
- **Install CHARMM and pyCHARMM (Part 2)**
- **Obtain [ALF version 3.2](https://github.com/ryanleehayes/alf) from github**
- **Install pyALF (Part 3)**

## 1. Creating conda environment to install and use CHARMM/pyCHARMM
- **You will need a base anaconda/miniconda installation: see [anaconda installation](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html).**
- **Follow steps in 1a to create the conda environment manually OR follow steps in 1b to create the conda environment from a YAML file.**
### 1a. Create a conda environment manually
- **Make a conda environment (See below for a shortcut using `conda env create -f <name_of_environment>.yml`)**<p>
`conda create -y -n <name_of_environment> python=3.9` # note python can be > 3.9
- **Activate this environment**<p>
`conda activate <name_of_environment>`
- **Install mamba as a faster conda**<p>
`conda install -y -c conda-forge mamba`
- **Install CUDA from NVIDIA. Pick one version compatible with your drivers as described below. You can see available [CUDA Toolkit packages (Table 3)](https://docs.nvidia.com/cuda/cuda-toolkit-release-notes/index.html#title-new-cuda-tools)**<p>
`mamba install -y -c nvidia cuda` # note this should install CUDA 12.1.1<p>
`mamba install -y -c "nvidia/label/cuda-12.0.0" cuda` # note this should install CUDA 12.0<p>
 - **Install needed packages to build CHARMM and pyCHARMM**<p>
`mamba install -y -c conda-forge gcc gxx gfortran make cmake binutils fftw openmpi openmm mpi4py sysroot_linux-64==2.17 readline==8.2 rdkit openbabel pandas pytorch jupyter_core jupyter_client jupyterlab jupyterlab_widgets jupyter_server jupyterlab_server jupyter_console jupyter jupytext biopython py3dmol mdtraj nglview jsonpickle pymol-open-source`
    
<div class="alert alert-block alert-warning">
<b>Note on CUDA Toolkit/Driver and Compiler Compatabilities:</b> In choosing the CUDA Toolkit you need to coordinate with the compatable CUDA Driver and compilers. The table below outlines these requirements. You should check with your systems manager regarding the installed CUDA Driver on the computer cluster/machine on which you plan to install CHARMM/pyCHARMM. However, you can also glean this information by running the command <i>nvidia-smi</i> on one of the nodes of your GPU-equipped computers. In this case the CUDA Driver will be displayed at the top of the output created from this command:</div><p>

`nvidia-smi`

```
+-----------------------------------------------------------------------------+
| NVIDIA-SMI 525.85.05    Driver Version: 525.85.05    CUDA Version: 12.0     |
|-------------------------------+----------------------+----------------------+
```

<div class="alert alert-block alert-warning">
Specifying that the Driver Version is 525.85.05. Thus, as seen from the table below, this Driver is compatable with CUDA 12.0, GCC >= 12.1 or Intel Compilers 2021.6.
</div>

| Toolkit Version | Minimum Required Driver | Recommended GCC | Recommended Intel Compilers
| :-: | :-: | :-: | :-:     
|CUDA 12.1.x|>=530.30.02|12.2|2021.7
|CUDA 12.0.0|>=525.85.05|12.1|2021.6
|CUDA 11.8.x|>=520.61.05|11|2021
|CUDA 11.7.x|>=515.48.07|11|2021
|CUDA 11.6.x|>=510.47.03|11|2021
|CUDA 11.5.x|>=495.29.05|11|2021
|CUDA 11.4.x|>=470.82.01|10.3|19.1
|CUDA 11.3.x|>=465.19.01|10.3|19.1
|CUDA 11.2.x|>=460.32.03|10.3|19.1
|CUDA 11.1 (11.1.0)|>=455.32|10.3|19.1
|CUDA 11.0 (11.0.3)|>=450.51.06|10.3|19.1
|CUDA 10|>= 440.33|10.2|18.0
|CUDA 9|>= 396.37|4.8.5|17.0
|CUDA 8|>= 375.26|4.8.2|15, 16

<div class="alert alert-block alert-warning">
<b>Newer drivers work with older CUDA versions, but older drivers do not work with newer CUDA version. Thus, if your driver is older than 525.85.05, the oldest available CUDA version, 11.3.1, should be compatible with the largest number of drivers. You can install it using
</div><p>

`mamba install -y -c "nvidia/label/cuda-11.3.1" cuda`

<div class="alert alert-block alert-warning">
<b>This CUDA version is incompatible with current versions of gcc, but version 10.4 works well so replace "gcc gxx gfortran" with "gcc==10.4 gxx==10.4 gfortran==10.4"
</div><p>

`mamba install -y -c conda-forge gcc==10.4 gxx==10.4 gfortran==10.4 make cmake binutils fftw openmpi openmm mpi4py sysroot_linux-64==2.17 readline==8.2 rdkit openbabel pandas pytorch jupyter_core jupyter_client jupyterlab jupyterlab_widgets jupyter_server jupyterlab_server jupyter_console jupyter jupytext biopython py3dmol mdtraj nglview jsonpickle pymol-open-source`

### 1b. Building the CHARMM/pyCHARMM compatable environment with a YAML file
 
`charmm_wcuda12.yml`
 
```YAML
name: charmm_wcuda12  # This represents the name you want to use for your conda environment
channels:
  - defaults
  - conda-forge
  - nvidia/label/cuda-12.0.0  # This tells conda to explicitly load cuda 12.0.0
  #- nvidia/label/cuda-12.1.1
  #- nvidia/label/cuda-11.8.0
  #- nvidia/label/cuda-11.7.1
dependencies:
  - python=3.9
  - mamba
  - cuda
  - ca-certificates
  - certifi
  - openssl
  - gcc
  - gxx
  - gfortran
  - make
  - cmake
  - binutils
  - fftw
  - openmpi
  - openmm
  - mpi4py
  - sysroot_linux-64==2.17
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
prefix: /home/brookscl/.conda/envs/charmm_wcuda12  # This corresponds to the path to this environment
```

- **You can edit this YAML file to add/change the cuda version as noted above. You can install this new conda environment with the command:**


`conda env create -f charmm_wcuda12.yml`


## 2. CHARMM and pyCHARMM installation once conda environment is installed and active.
- **Go to CHARMM source root and build CHARMM with configure**

<blockquote>

```csh
conda activate charmm_wcuda12
cd <charmm_root>
mkdir build_charmm
cd build_charmm
# Build CHARMM with BLaDe, FFTDOCK,DOMDEC (-u) and OpenMM (default)
export FFTW_HOME=$CONDA_PREFIX # bash syntax
setenv FFTW_HOME $CONDA_PREFIX # csh syntax
../configure --with-blade --with-fftdock -u  -D nvcc_ptx_target=52 -p <charmm_install_path>
make -j <n> install
```

</blockquote>

- **charmm_wcuda12 should be replaced with the name of your conda virtual environemnt**
- **\<charmm_root\> is the path to the charmm top level tree**
- **\<charmm_install_path\> is the path where you want the CHARMM installation to reside**
- **\<n\> is the number of cores to use in compiling the code**
- **-D nvcc_ptx_target=52 is required for older GPUs like GTX980s**

### pyCHARMM is built from the same source and can be built in the same build directory

<blockquote>

```csh
conda activate charmm_wcuda12
cd <charmm_root>
cd build_charmm
rm -r *   # Clean the build directory
# Build CHARMM with BLaDe, FFTDOCK,DOMDEC (-u) and OpenMM (default)
export FFTW_HOME=$CONDA_PREFIX # bash syntax
setenv FFTW_HOME $CONDA_PREFIX # csh syntax
../configure --with-blade --with-fftdock --without-mpi --as-library  -D nvcc_ptx_target=52 -p <pycharmm_install_path>
make -j <n> install
cd <charmm_root>
pip install `pwd`/tool/pycharmm  # Installs the pyCHARMM modules in your current environment
export CHARMM_LIB_DIR=<pycharmm_install_path>/lib # bash syntax
setenv CHARMM_LIB_DIR <pycharmm_install_path>/lib # csh syntax
conda env config vars set CHARMM_LIB_DIR=<pycharmm_install_path>/lib  # every time when this conda environment (charmm_wcuda12) is activated, the environmental variable CHARMM_LIB_DIR is there automatically.
```

</blockquote>
    
- **\<pycharmm_install_path\> is the path where you want the pyCHARMM installation to reside**

## 3. pyALF installation
- **Download [ALF version 3.2](https://github.com/ryanleehayes/alf) from github.**
- **Go to pyALF source root and build as follows:**

<blockquote>

```bash
conda activate charmm_wcuda12
cd <alf_root>
export ALF_SOURCE_DIR=`pwd` # bash syntax
setenv ALF_SOURCE_DIR `pwd` # csh syntax
# Compile ALF executables in wham and dca
cd $ALF_SOURCE_DIR/alf/wham
bash Clean.sh
cmake ./
make wham
cd $ALF_SOURCE_DIR/alf/dca
bash Clean.sh
cmake ./
make all
# Install pyALF in current conda virtual environment
cd $ALF_SOURCE_DIR
pip install -e .
python -c "import alf"
```

</blockquote>

- **charmm_wcuda12 should be replaced with the name of your conda virtual environemnt**
- **\<alf_root\> is the path to the alf top level tree**
