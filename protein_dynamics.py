#!/usr/bin/env python
# coding: utf-8

# # Run protein dynamics with PBCs using PME
# ## This notebook builds on the previous protein_setup.ipynb and sets up a molecular dynamics simulation using Particle-Mesh Ewald running on CHARMM/OpenMM or CHARMM/BLaDE.
# ## Note you need top have run the previous notebook and have created the solvated psf and pdb to use this notebook.

# In[ ]:


# This script provides a simple example of running molecular dynamics of 
# a solvated protein structure created in an earlier notebook tutorial
# further illustrating the functionality exposed in pyCHARMM.
#  copyright C.L. Brooks III, June 20, 2022

# These are general python modules needed for this  tutorial
import os
import sys
import numpy as np

# These are a subset of the pycharmm modules that were installed when
# pycharmm was installed in your python environment
import pycharmm
import pycharmm.generate as gen
import pycharmm.ic as ic
import pycharmm.coor as coor
import pycharmm.energy as energy
import pycharmm.dynamics as dyn
import pycharmm.nbonds as nbonds
import pycharmm.minimize as minimize
import pycharmm.crystal as crystal
import pycharmm.image as image
import pycharmm.psf as psf
import pycharmm.read as read
import pycharmm.write as write
import pycharmm.settings as settings
import pycharmm.cons_harm as cons_harm
import pycharmm.cons_fix as cons_fix
import pycharmm.select as select
import pycharmm.shake as shake

from pycharmm.lib import charmm as libcharmm


# ## The following are some helper functions for setting-up and running the dynamics

# In[ ]:


################################
# Ensure that FFT grid is product of small primes 2, 3, 5
def is_factor(n):
    if (n % 2 != 0): return False  # favors even number
    while n:
        flag = False
        for x in (2,3,5):
            if n % x == 0:
               n = n / x
               flag = True
               break

        if flag: continue
        break

    if n == 1: return True
    return False

def checkfft(n, margin = 5):
    n = int(n) + margin
    while 1:
        if is_factor(n): break
        else: n += 1
    return n


# In[ ]:


#################################
def setup_PBC(boxhalf=0.0, protein_segments=[],solvent_resname='TIP3',ions=[],blade=False):
    """input: boxhalf [0.0]
              solute  []
              solvent_resname ['']
              ions []
              blade [False]
    defines the periodic boundary conditions for a cubic volume of boxsize. 
    Uses: crystal_define_cubic(), crystal.build(), image.setup_residue,
    image.setup_segment to construct symmetry operations. 

    If global variable openmm is true
    the image centering is at [boxhalf,boxhalf,boxhalf] otherwise at [0,0,0].
    """
    crystal.define_cubic(boxhalf*2)
    crystal.build(boxhalf)

    if blade: boxhalf = 0.0 # center at origin for blade
    for segment in protein_segments:
        image.setup_segment(boxhalf,boxhalf, boxhalf, segment)
    if len(solvent_resname)>0: image.setup_residue(boxhalf,boxhalf, boxhalf, solvent_resname)
    for ion in ions:
        image.setup_residue(boxhalf, boxhalf, boxhalf, ion)
    # for systems using centering not at origin, translate coordinates by halfbox
    xyz = coor.get_positions()
    xyz += boxhalf
    coor.set_positions(xyz)
    print('Coordinates translated by {} A in each dimension to be consistent with image centering'\
          .format(boxhalf))

    return


# ## Set-up system - topology and parameter files

# In[ ]:


# Read in the topology (rtf) and parameter file (prm) for proteins
# equivalent to the CHARMM scripting command: read rtf card name toppar/top_all36_prot.rtf
read.rtf('toppar/top_all36_prot.rtf')
# equivalent to the CHARMM scripting command: read param card flexible name toppar/par_all36m_prot.prm
read.prm('toppar/par_all36m_prot.prm', flex=True)

# stream in the water/ions parameter using the pycharmm.lingo module
# equivalent to the CHARMM scripting command: stream toppar/toppar_water_ions.str
pycharmm.lingo.charmm_script('stream toppar/toppar_water_ions.str')


# ## Read in psf and coordinates of solvated systems
# ## Find the segment names for the solute (first n-2 segments assumed), solvent (penultimate segment assumed) and ions (last segment)

# In[ ]:


# Some choices for pdbids
pdbid = '4pti'  # bovine pancreatic trypsin inhibitor - an x-ray structuture
#pdbid = '6pti'  # bovine pancreatic trypsin inhibitor - an x-ray structuture
#pdbid = '5wyo'  # HDEA a dimeric (two chain) protein - an NMR structure

# Read the psf and coordinates for the solvated peptide
# Read psf card name pdb/adp+wat.psf
read.psf_card('pdb/{}+wat.psf'.format(pdbid))
# read coor pdb name pdb/adp+wat_min.pdb
read.pdb('pdb/{}+wat_min.pdb'.format(pdbid),resid=True)


# ## Set-up box size and periodic boundary conditions using the python function setup_PBC

# In[ ]:


# Now setup periodic boundaries
# boxsize
stats = coor.stat()
xsize = stats['xmax'] - stats['xmin']
ysize = stats['ymax'] - stats['ymin']
zsize = stats['zmax'] - stats['zmin']
boxsize = max(xsize, ysize, zsize)

# half box size
boxhalf = boxsize / 2.0
# Note we could probably do something to extract the information passed to setup_PPC using pyCHARMM functions
# but I didn't have the time so I just made some generic thing
setup_PBC(boxhalf=boxhalf, protein_segments=['PROA', 'PROB', 'PROC'],
          solvent_resname='TIP3',ions=['CLA', 'SOD', 'POT'],blade=False)


# ## In this section we set-up our nonbonded parameters with cutoff schemes, etc.
# ## We illustrate two ways to input non-bonded parameters as well as the use of the pyCHARMM _settings.set_warn_level(wrnlev)_ command to alter the warning levels and reduce the output. 
# ## Note we use two different cutoff schemes (vswitch and vfswitch)

# In[ ]:


# Set-up non-bonded parameters
# Now specify nonbonded cutoffs for solvated box
cutnb = min(boxhalf,12)
cutim = cutnb
ctofnb = cutnb - 1.0
ctonnb = cutnb - 3.0
# Determine the appropriate cubic fft grid for this boxsize
fft = checkfft(n=np.ceil(boxhalf)*2,margin=0)
# Set-up the parameters
nb_wPME_vsw = pycharmm.NonBondedScript(cutnb=cutnb, cutim=cutim,
                                       ctonnb=ctonnb, ctofnb=ctofnb,
                                       cdie=True, eps=1,
                                       atom=True, vatom=True,
                                       switch=True, vfswitch=False, vswitch=True,
                                       inbfrq=-1, imgfrq=-1,
                                       ewald=True,pmewald=True,kappa=0.32,
                                       fftx=fft,ffty=fft,fftz=fft,order=4)
# Let's set the wrnlev to 0 to avoid the large output
old_wrnlev = settings.set_warn_level(0)
nb_wPME_vsw.run()
settings.set_warn_level(old_wrnlev)
energy.show()
# Let's also set-up a set of nonbonded parameters using vfswitch instead of vswitch
nb_wPME_vfsw_dict = {'cutnb':cutnb, 
                     'cutim':cutim,
                     'ctonnb':ctonnb, 
                     'ctofnb':ctofnb,
                     'cdie':True,
                     'eps':1,
                     'atom':True, 'vatom':True,
                     'switch':True, 'vfswitch':True, 'vswitch':False,
                     'inbfrq':-1, 'imgfrq':-1,
                     'ewald':True,'pmewald':True,'kappa':0.32,
                     'fftx':fft,'ffty':fft,'fftz':fft,'order':4}
# Let's set the wrnlev to 0 to avoid the large output
old_wrnlev = settings.set_warn_level(0)
pycharmm.NonBondedScript(**nb_wPME_vfsw_dict).run()
settings.set_warn_level(old_wrnlev)
energy.show()
# Now go back to the original nonbonded parameters
# Let's set the wrnlev to 0 to avoid the large output
old_wrnlev = settings.set_warn_level(0)
nb_wPME_vsw.run()
settings.set_warn_level(old_wrnlev)
energy.show()


# ## This function will run a short md start and restart using either CHARMM/OpenMM interface or the CHARMM/BLaDE interface in pyCHARMM

# In[ ]:


def run_md(useomm=False,useblade=False,nequil=1000,nsteps=5000,nsavc=100,leap=True,lang=True):
    if useomm: append='omm'
    elif useblade: append='blade'
    dyn.set_fbetas(np.full((psf.get_natom()),1.0,dtype=float))
   
    res_file = pycharmm.CharmmFile(file_name='res/{}.res'.format(pdbid), file_unit=2,
                                   formatted=True,read_only=False)
    lam_file = pycharmm.CharmmFile(file_name='res/{}.lam'.format(pdbid), 
                                   file_unit=3,
                                   formatted=False,read_only=False)
    my_dyn = pycharmm.DynamicsScript(leap=leap, lang=lang, start=True,
                                     nstep=nequil, timest=0.002,
                                     firstt=298.0, finalt=298.0, tbath=298.0,
                                     tstruc=298.0,
                                     teminc=0.0, twindh=0.0, twindl=0.0,
                                     iunwri=res_file.file_unit,
                                     iunlam=lam_file.file_unit,
                                     inbfrq=-1, imgfrq=-1,
                                     iasors=0, iasvel=1, ichecw=0, iscale=0,
                                     iscvel=0,echeck=-1, nsavc=0, nsavv=0, nsavl=0, ntrfrq=0,
                                     isvfrq=nsavc,
                                     iprfrq=2*nsavc, nprint=nsavc, ihtfrq=0, ieqfrq=0,
                                     ilbfrq=0,ihbfrq=0,
                                     omm=useomm, blade=useblade)
    my_dyn.run()

    res_file.close()
    lam_file.close()
    # open unit 2 write form name res/{}.res
    res_file = pycharmm.CharmmFile(file_name='res/{}.res'.format(pdbid), file_unit=2,
                                   formatted=True,read_only=False)
    lam_file = pycharmm.CharmmFile(file_name='res/{}.lam'.format(pdbid), 
                                   file_unit=3,
                                   formatted=False,read_only=False)
    # open unit 1 write file name dcd/{}.dcd
    dcd_file = pycharmm.CharmmFile(file_name='dcd/{}_{}.dcd'.format(pdbid,append), file_unit=1,
                                   formatted=False,read_only=False)

    my_dyn = pycharmm.DynamicsScript(leap=leap, lang=lang, start=False, restart = True,
                                     nstep=nsteps, timest=0.002,
                                     firstt=298.0, finalt=298.0, tbath=298.0,
                                     tstruc=298.0,
                                     teminc=0.0, twindh=0.0, twindl=0.0,
                                     iunwri=res_file.file_unit,
                                     iunrea=res_file.file_unit,
                                     iuncrd=dcd_file.file_unit,
                                     iunlam=lam_file.file_unit,
                                     inbfrq=-1, imgfrq=-1,
                                     iasors=0, iasvel=1, ichecw=0, iscale=0,
                                     iscvel=0,echeck=-1, nsavc=nsavc, nsavv=0, nsavl=0, ntrfrq=0,
                                     isvfrq=nsavc,
                                     iprfrq=2*nsavc, nprint=nsavc, ihtfrq=0, ieqfrq=0,
                                     ilbfrq=0,ihbfrq=0,
                                     omm=useomm, blade=useblade)
    my_dyn.run()

    res_file.close()
    lam_file.close()
    dcd_file.close()
    return


# # Using CHARMM/BLaDE interface through pyCHARMM
# ## In this section we set-up a short dynamics run equilibration and production using the CHARMM/BLaDE gpu accelerated engine in CHARMM (if a gpu isn't present this should be skipped - we use torch functions to check whether gpu is present).

# In[ ]:


# Set-up short dynamics
if not os.path.isdir('res'): os.system('mkdir res')
if not os.path.isdir('dcd'): os.system('mkdir dcd')

# Check to see if cuda is available to run BLaDE
import torch
cuda  = torch.cuda.is_available()
num_devices = torch.cuda.device_count()
if cuda:
    print('Running CHARMM/BLaDE MD example on computer with {} CUDA devices'.format(num_devices))
    print('Running on device {} which is a {}'.format(torch.cuda.current_device(),torch.cuda.get_device_name(torch.cuda.current_device())))
    run_md(useblade = 'prmc pref 1 iprs 100 prdv 100')
else: print('Example not run, no CUDA devices available')


# # View your trajectory
# ## Using vmd in your terminal window issue the following command:
# > *vmd pdb/4pti+wat.psf dcd/4pti_blade.dcd*

# # Using CHARMM/OpenMM interface through pyCHARMM
# ## In this section we set-up a short dynamics run equilibration using the CHARMM/OpenMM gpu accelerated engine in CHARMM (if a gpu isn't present this should choose the platform _'cpu'_ and run there.

# In[ ]:


# Set-up short dynamics
if not os.path.isdir('res'): os.system('mkdir res')
if not os.path.isdir('dcd'): os.system('mkdir dcd')
# Check to see if cuda is available to run BLaDE
import torch
cuda  = torch.cuda.is_available()
num_devices = torch.cuda.device_count()
if cuda:
    print('Running CHARMM/OpenMM MD example on computer with {} CUDA devices'.format(num_devices))
    print('Running on device {} which is a {}'.format(torch.cuda.current_device(),torch.cuda.get_device_name(torch.cuda.current_device())))
else: print('No CUDA devices available, using either CPU or OpenCL')
run_md(useomm = 'gamma 2 prmc pref 1 iprsfrq 100')


# # View your trajectory
# ## Using vmd in your terminal window issue the following command:
# > *vmd pdb/4pti+wat.psf dcd/4pti_omm.dcd*
