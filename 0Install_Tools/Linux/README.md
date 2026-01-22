## 0. Installing needed tools for CHARMM/pyCHARMM
### In order to use CHARMM, pyCHARMM, and pyALF you will need to:
- **Create a conda environment capable of building CHARMM, pyCHARMM, and pyALF (Part 1)**
- **Install the MMTSB ToolSet from [MMTSB](https://feig.bch.msu.edu/mmtsb/Main_Page). Follow the instructions to install the package.**
- **Obtain the CHARMM software (free to academics and government labs) from [AcademicCHARMM](https://academiccharmm.org/program). Follow the directions below to build a conda environment capable of installing CHARMM/pyCHARMM.**
- **Install CHARMM and pyCHARMM (Part 2)**
- **Obtain [ALF version 3.2](https://github.com/ryanleehayes/alf) from github**
- **Install pyALF (Part 3)**

## 1. Creating conda environment to install and use CHARMM/pyCHARMM
- **You will need a base anaconda/miniconda/conda-forge (conda-forge recommended) installation: see [conda-forge installation](https://conda-forge.org/download/#:~:text=Miniforge%20is%20the%20preferred%20conda,conda%20create%20or%20mamba%20create%20.).**
- **Follow steps in 1 to create the conda environment from the included YAML file.**
> **Note: We have added new support for [OpenMM PyTorch plugin](https://github.com/openmm/openmm-torch).**

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
|CUDA 13.x|>=580.xx|12.x,13.x,14.x|unknown
|CUDA 12.x|>=526.60.13|12.x,13.x|unknown
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

### 1. Building the CHARMM/pyCHARMM compatable environment with a YAML file
 
`charmm_gpu_env.yml`
 
```YAML
# File charmm_gpu_env.yml
name: charmm_gpu_env
channels:
  - conda-forge
  - nvidia
  - pytorch
dependencies:
  # =====================================================
  # CUDA 12.6 Toolkit (PINNED - do not update)
  # =====================================================
  - cuda-version=12.6.*
  - cuda-toolkit=12.6.*
  - cuda-cudart
  - cuda-nvcc
  - cuda-libraries
  - cuda-nvrtc
  - cudnn
  # =====================================================
  # Compilers (GCC 12.x for CHARMM compatibility)
  # =====================================================
  - gcc=12.2.0
  - gxx=12.2.0
  - gfortran=12.2.0
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
  - openmm-torch=*=*cuda126*
  # =====================================================
  # Math libraries needed for CHARMM
  # =====================================================
  - fftw
  - clfft
  # =====================================================
  # TorchANI for support of ML  QM/ML potentials
  # =====================================================
  - torchani=*=*cuda126*
  # =====================================================
  # Utilities usefull for CHARMM/pyCHARMM
  # =====================================================
  - readline
  - scipy
  - jsonpickle
  - pip
  - pdoc
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
    - crimm
    - fastmbar
```

- **You can edit this YAML file to add/change the cuda version as noted above. You can install this new conda environment with the command:**

<blockquote>
 
```csh 
export CONDA_OVERRIDE_CUDA="12.9" # enables build on head nodes w/o GPU
setenv CONDA_OVERRIDE_CODA 12.9   # enables build on head nodes w/o GPU
conda env create -n <your_charmm_environment_name> -f charmm_gpu_env.yml
```

</blockquote>

## 2. CHARMM and pyCHARMM installation once conda environment is installed and active.
- **Go to CHARMM source root and build CHARMM with configure**

<blockquote>

```csh
conda activate charmm_gpu_env
cd <charmm_root>
mkdir build_charmm
cd build_charmm
# Build CHARMM with BLaDe, FFTDOCK,DOMDEC (-u) and OpenMM (default)
export FFTW_HOME=$CONDA_PREFIX # bash syntax
setenv FFTW_HOME $CONDA_PREFIX # csh syntax
../configure --with-blade --with-fftdock -u  -D nvcc_ptx_target="52:75" -p <charmm_install_path>
make -j install
```

</blockquote>

- **charmm__gpu_env should be replaced with the name of your conda virtual environemnt**
- **\<charmm_root\> is the path to the charmm top level tree**
- **\<charmm_install_path\> is the path where you want the CHARMM installation to reside**
- **`-D nvcc_ptx_target=52` is required for older GPUs like GTX980s default is 75 so `-D nvcc_ptx_target` is not necessary if building for more modern Nvidia GPUs.**

### pyCHARMM is built from the same source and can be built in the same build directory. If you want to ues pyCHARMM with MPI in python use this build

<blockquote>

```csh
conda activate charmm_gpu_env
cd <charmm_root>
cd build_charmm
rm -r *   # Clean the build directory
# Build CHARMM with BLaDe, FFTDOCK,DOMDEC (-u) and OpenMM (default)
export FFTW_HOME=$CONDA_PREFIX # bash syntax
setenv FFTW_HOME $CONDA_PREFIX # csh syntax
../configure --with-blade --with-fftdock --without-mpi --as-library  -D nvcc_ptx_target=52 -p <pycharmm_install_path>
make -j install
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

- **charmm_gpu_env should be replaced with the name of your conda virtual environemnt**
- **\<alf_root\> is the path to the alf top level tree**
