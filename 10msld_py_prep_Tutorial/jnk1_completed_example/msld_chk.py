#! /usr/bin/env python

##
## Check molfile and toppar files before getting started
##

import numpy as np
import glob,os

class CHK_Error(Exception):
    import sys
    sys.exit

def MsldCHK(molfile):
    """
    Check the following items:
     - each atom is uniquely named (within a single mol2 file)
     - each mol2 file has an associated rtf file
       (parse ParamChem stream files into rtf/prm files)
    """

    # Load in the file names
    fp=open(molfile,'r')
    mols=[]
    for line in fp:
        mols.append(line.split()[0]) 
    fp.close()

    # Some quick file checks:
    for mol in range(len(mols)):
        # Check to make sure the mol2 files exist
        if len(glob.glob(mols[mol]+'.mol2')) == 0:
            raise CHK_Error(mols[mol]+'.mol2 does not exist - check and resubmit')
        # Check to see if every atom is uniquely named - exit if not
        chk=0
        fp=open(mols[mol]+'.mol2','r')
        line=fp.readline()
        aname=[]
        while line:
            if line[0:13] == '@<TRIPOS>ATOM':
                line=fp.readline()
                while line[0:13] != '@<TRIPOS>BOND':
                    tmp=line.split()[1]
                    if not (tmp in aname):
                        aname.append(tmp)
                    else:
                        chk=1
                    line=fp.readline()
                break
            else:
                line=fp.readline()
        fp.close()
        if chk:
            raise CHK_Error(mols[mol]+'.mol2 atom names are NOT unique - fix this in both structure and toppar files - then resubmit')

        ## If atom names are not unique - you must rename them before

        # Check to make sure rtf/prm files exist (parse str files)
        if len(glob.glob(mols[mol]+'.rtf')) == 0:
            # check to see if a ParamChem stream file exists:
            if len(glob.glob(mols[mol]+'.str')) == 0:
                raise CHK_Error(mols[mol]+'.rtf does not exist - check and resubmit')
            else:
                # assume the stream file is from paramchem
                fp=open(mols[mol]+'.str','r')
                rp=open(mols[mol]+'.rtf','w')
                pp=open(mols[mol]+'.prm','w')

                line=fp.readline()
                while line:
                    # rtf first
                    if line[0:8] == 'read rtf':
                        line=fp.readline()
                        while line[0:3] != 'END':
                            rp.write(line)
                            line=fp.readline()
                        rp.write(line)
                        line=fp.readline()
                    # prm next
                    elif line[0:10] == 'read param':
                        line=fp.readline()
                        while line[0:3] != 'END':
                            pp.write(line)
                            line=fp.readline()
                        pp.write(line)
                        line=fp.readline()
                    else:
                        line=fp.readline()

                fp.close()
                rp.close()
                pp.close()

        # Check that atom names in the mol2 file match the atom names in the rtf
        a2name=[]
        fp=open(mols[mol]+'.rtf','r')
        for line in fp:
            if line[0:4] == 'ATOM':
                a2name.append(line.split()[1])
        fp.close()

        cnt=0
        for at1 in aname:
            for at2 in a2name:
                if at1 == at2:
                    cnt+=1
        if cnt != len(aname):
            raise CHK_Error('Atom names in '+mols[mol]+'.mol2 do not match those found in '+mols[mol]+'.rtf! Check and resubmit.')

        # mols[mol] ready to go

    return

