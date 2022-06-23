#! /usr/bin/env python

####
#### Executable script to build MSLD ready ligand files
#### JV,LC 01/2022
####

import msld_chk
import msld_mcs
import msld_crn
import msld_prm
import msld_wrt
import glob

###
### This script is executed in 2 steps to (1st) build the MCSS
### and (2nd) perform charge renormalization. This allows the user
### to check that the identified MCSS is correct with vis_check.py.
### Thus, the user should manually call msld_py_prep.py twice.
###
### All ligand structure files (mol2) and toppar files must be 
### available prior to running this script
###
### The inFrag list of lists allows you to move core atoms into 
### alchemical fragments at specific sites. Each nested list 
### corresponds to a single site attached to the ligand core.
###
### The AnCore list of lists allows you to move (non-DUM) atoms
### listed as "anchor atoms" (atoms connecting core and fragment
### components) into the core upon charge renormalization.
###


#####################################################################
## (1) Define System and File Variables

sysname = "jnk1"                          # name of future output files
molfile = "mol_list.txt"                  # list of mol2 file names

mcsout = 'MCS_for_MSLD.txt'               # MCS output filename
outdir = 'build.'+sysname                 # MSLD output directory

cgenff=True                               # Are CGenFF/ParamChem parameters being used?

inFrag=[[]]  # reflig core atoms to include in each fragment at each site (list of nsub lists)
AnCore=[[]]  # anchor atoms at each site to include in the core (list of nsub lists)

#####################################################################
if len(glob.glob(mcsout)) == 0:
    ## (2) Check molfile and toppar files before getting started
    msld_chk.MsldCHK(molfile)
    print("chk finished")
    
    #####################################################################
    ## (3) Maximum Common SubStruct Search with bonded-environments
    ## "mcsout" = results of the search; edit this file to manual edit ligand splicing
    ## cutoff = RMSD & distance cutoff to differentiate different atoms
    ## change debug to True to get more stdout printed 
    
    reflig = msld_mcs.MsldMCS(molfile,mcsout,cutoff=0.8,debug=False)
    print("MCS results printed to "+mcsout)
    print("Reference Ligand is "+reflig)
    quit()


#####################################################################
## (4) Perform Charge-Renormalization 
## To manually move atoms from the core into alchemical fragments, 
## update the "inFrag" variable above. "Anchor atoms" (and connected Hs)
## are automatically included in each alchemical fragment unless 
## specifically stated to be "AnCore"

msld_crn.MsldCRN(mcsout,outdir,inFrag,AnCore,ChkQChange=True,verbose=True,debug=False)


#####################################################################
## (5) Write Ligand Parameters & the Charmm ALF input scripts
msld_prm.MsldPRM(outdir,cgenff,verbose=False,debug=False)
msld_wrt.writeALF_Files(sysname,outdir,cgenff)


## Final Notes to the user
print("default TOPPAR parameters copied into build."+sysname+". Check to make sure these work for your system!")


## FINISHED

