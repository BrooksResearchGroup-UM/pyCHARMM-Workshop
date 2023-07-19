#! /usr/bin/env python

##
## Maximum Common Substructure search with bonded-environments
##

import numpy as np
from copy import deepcopy


### Requirments for successful use:
###     ** Molecules MUST be spacially aligned
###     ** Atom names (within each mol2 file) must be unique for each atom
###        (atom names can be similar between different mol2 files)
###     ** Atom names cannot have the "+" symbol in their name



def MsldMCS(molfile,mcsout,cutoff=0.8,debug=False):
    """
    Use bond-connectivity, atomtype definitions, and distance metrics to identify the
    core and fragment atoms for a set of supplied molecules.

    ** Molecules should be in a mol2 format; Atomtype parameters should be in a CHARMM RTF format
    ** Molecules must be spacially aligned (such that common core atoms are within 1.0A of each other)
    ** Atom names must be unique for each atom per molecule; cross molecule atom name similarities are ok
    ** Atom names cannot have the "+" character in their name

    """

    ms=[]    # ms is a list of mol2 files
    atoms=[] # atoms is a list of lists of all atom names
    heavy=[] # heavy is a list of lists of HEAVY atom names only
    xyzs=[]  # xyzs is a list of dictionaries for coordinates
    bonds=[] # bonds is a list of np.arrays of bonds (0/1 = no/yes bond)
    types=[] # types is a list of dictionaries of atom names -> atom types
    
    ## (1) Read "mol_list"
    fp=open(molfile,'r')
    for line in fp:
        ms.append(line.rstrip())
    fp.close()
    
    ## (2) Read each mol file to load in the coordinates for each atom and its bonds
    # loop over each file
    for m in range(len(ms)):
        atoms.append([])
        heavy.append([])
        xyzs.append({})
        fp=open(ms[m]+'.mol2','r')
        line=fp.readline()
        while line:
            if line.rstrip() == '@<TRIPOS>ATOM':
                # extract atom names and coordinates
                line=fp.readline()
                while line.rstrip() != '@<TRIPOS>BOND':
                    tmp=line.split()
                    atoms[m].append(tmp[1])
                    xyzs[m][tmp[1]]={'X':float(tmp[2]),'Y':float(tmp[3]),'Z':float(tmp[4]),'Elem':tmp[5].split('.')[0]}
                    if tmp[5].split('.')[0] != 'H':
                        heavy[m].append(tmp[1])
                    line=fp.readline()
            if line.rstrip() == '@<TRIPOS>BOND':
                bonds.append(np.zeros((len(atoms[m]),len(atoms[m])),dtype=int))
                line=fp.readline()
                # create a bond matrix
                while line:
                    if line == '\n':
                        # empty line
                        break
                    if line[0] == '@':
                        # signifies end of "BOND" section of mol2
                        break
                    tmp=line.split()
                    bonds[m][int(tmp[1])-1][int(tmp[2])-1]=1
                    bonds[m][int(tmp[2])-1][int(tmp[1])-1]=1
                    line=fp.readline()
            line=fp.readline()
        
        fp.close()
    
    # #debug#
    # if debug:
    #     for m in range(len(ms)):
    #         print("Molecule",ms[m])
    #         print(atoms[m])
    #         print(xyzs[m])
    #         #print(bonds[m])
    #         for l in range(bonds[m].shape[0]):
    #             print(bonds[m][l])
    #         print("")
    #         print("")
    # #debug#

    ## (3) Read the str file to get the atom types
    for m in range(len(ms)):
        types.append({})
        numLPs=0
        #fp=open(ms[m]+'.str','r')  # msld_chk splits str files into rtf/prm files
        fp=open(ms[m]+'.rtf','r')
        line=fp.readline()
        while line:
            if line[0:4] == 'ATOM':
                tmp=line.split()
                # check for LonePair site
                if tmp[1][0:2] == 'LP':
                    numLPs+=1
                else:
                    types[m][tmp[1]] = tmp[2]
            line=fp.readline()
        # chk for same number of atoms
        if len(types[m]) != len(atoms[m]):
            print("ERROR: Inconsistent # of atoms between types[m] and atoms[m] after reading RTF file",ms[m]+'.str')
            print("len(types[m])=",len(types[m]),"; len(atoms[m])=",len(atoms[m]))
            quit()
    
    # #debug#
    # if debug:
    #     for m in range(len(types)):
    #         print(types[m])
    #         print("")
    # #debug#
    
    
    ## (4) Start finding the core (heavy atom only search)
    # function to create an atom's bond pattern 
    def getBonded(atomnum,molnum):
        # first pull out the atom indices that "atomnum" is bonded to
        b1=[]
        for at in range(bonds[molnum].shape[0]):
            if bonds[molnum][atomnum][at] == 1:
                b1.append(at)
        # now convert those indices to atom types
        b2=[]
        for at in range(len(atoms[molnum])):
            if at in b1:
                b2.append(types[molnum][atoms[molnum][at]])
        # make the dictionary current atomtype -> bonding neighbors
        bonded=['1','2']
        bonded[0]=types[molnum][atoms[molnum][atomnum]]  # bonded[0] = current atom's type
        bonded[1]=b2                                     # bonded[1] = neighbor types
    
        return bonded
    
    
    ## (i) Find all atom type matches
    matches=[]
    for m1 in range(len(ms)):
        matches.append({})
        for at1 in range(len(atoms[m1])):
            # compare only heavy atoms
            if xyzs[m1][atoms[m1][at1]]['Elem'] == 'H':
                continue
            matches[-1][atoms[m1][at1]]=[]
            patt1=getBonded(at1,m1)
            for m2 in range(len(ms)):
                matches[-1][atoms[m1][at1]].append([])
                if m2 == m1:
                    continue
                for at2 in range(len(atoms[m2])):
                    if xyzs[m2][atoms[m2][at2]]['Elem'] == 'H':
                        continue
                    patt2=getBonded(at2,m2)
                    # compare patterns 1 and 2 - store if they match
                    # if at1 and at2 types match and they both have the same number of neighbors
                    if (patt1[0] == patt2[0]) and (len(patt1[1]) == len(patt2[1])):
                        atsum=0
                        for at in patt1[1]:
                            if at in patt2[1]:
                                atsum+=1
                                # remove at from patt2[1] to prevent duplicate matches
                                for i in range(len(patt2[1])):
                                    if at == patt2[1][i]:
                                        patt2[1].pop(i)
                                        break
                        if atsum == len(patt1[1]):
                            # store possible match
                            matches[m1][atoms[m1][at1]][m2].append(atoms[m2][at2])
    
    # debug functions
    def printListDict(struct):
        for m in range(len(ms)):
            print("Molecule",m)
            keys = struct[m].keys()
            for k in keys:
                print(k,struct[m][k])
            print("")
    def printListList(struct):
        for m in range(len(ms)):
            print("Molecule",m)
            for l in struct[m]:
                print(l)
            print("")
    
    if debug:
        print("VERY BEGINNING")
        printListDict(matches)
    
                
    ## (ii) start sorting through the matches
    ## (iiA) remove entries with more than one empty list 
    for m1 in range(len(ms)):
        droplist=[]
        ibuff=-1
        for at1 in range(len(atoms[m1])):
            empty=0
            if xyzs[m1][atoms[m1][at1]]['Elem'] == 'H':
                continue
            ibuff+=1
            for m2 in range(len(ms)):
                if matches[m1][atoms[m1][at1]][m2] == []:
                    empty+=1
            if empty > 1:
                droplist.append(atoms[m1][at1])
        # pop off items from matches dict
        for i in droplist:
            matches[m1].pop(i)
    
    # debug
    # if debug:
    #     printListDict(matches)
                
    ## (iiB) remove entries with exactly one match and add them to cores
    cores=[] # core atoms and their matches
    for m1 in range(len(ms)):
        droplist=[]
        cores.append({})
        for at1 in range(len(atoms[m1])):
            if atoms[m1][at1] in matches[m1].keys():
                chk=0
                for m2 in range(len(ms)):
                    if m1 != m2:
                        if len(matches[m1][atoms[m1][at1]][m2]) == 1:
                            chk+=1
                if chk == len(ms) - 1: # then we have a core atom match
                    droplist.append(atoms[m1][at1])
        for i in droplist:
            at=matches[m1].pop(i)
            cores[m1][i]=at
    
    # # debug
    if debug:
        print("AFTER QUICK EMPTY or 1 MATCH CHECKS")
        print("")
        printListDict(cores)
        print("")
        print("")
        printListDict(matches)
        print("")
        print("")
    
    ## (iiC) Make sure all core lists have the same number of atoms: Correct inconsistencies
    chkcores=[]
    for m1 in range(len(ms)):
        chkcores.append([])
        m1cores=list(cores[m1].keys())
        for m2 in range(len(ms)):
            chkcores[m1].append([])
            if m1 != m2:
                m2cores=list(cores[m2].keys())
                # generate list of m1 core matches from cores[m2] key-values
                m2matches=[]
                for at2 in m2cores:
                    m2matches.append(cores[m2][at2][m1][0])
                for at1 in m1cores:
                    if not (at1 in m2matches):
                        chkcores[m1][m2].append(at1)
                        # at1 = key in m1
                        line=deepcopy(cores[m1][at1])
                        line[m1]=[at1]
                        cores[m2][line[m2][0]]=line
                        cores[m2][line[m2][0]][m2]=[]
    
    ## (iiD) for remaining atoms, try to match by common bonds to core atoms
    for m1 in range(len(ms)):
        for at1 in range(len(atoms[m1])):
            m1list=[]
            if atoms[m1][at1] in matches[m1].keys():
                # identify bonded core atoms
                for at2 in range(len(atoms[m1])):
                    if atoms[m1][at1] == atoms[m1][at2]:
                        continue
                    if atoms[m1][at2] in cores[m1].keys():
                        if bonds[m1][at1][at2] == 1: # meaning there's a bond
                            m1list.append(atoms[m1][at2])
                m2list=[]
                for m2 in range(len(ms)):
                    m2list.append([])
                    if m1 != m2:
                        for at2 in range(len(atoms[m2])):
                            if atoms[m2][at2] in matches[m1][atoms[m1][at1]][m2]:
                                m2list[-1].append([])
                                # identify bonded core atoms
                                for at3 in range(len(atoms[m2])):
                                    if atoms[m2][at2] == atoms[m2][at3]:
                                        continue
                                    if atoms[m2][at3] in cores[m2].keys():
                                        if bonds[m2][at2][at3] == 1: # there's a bond
                                            m2list[-1][-1].append(atoms[m2][at3])
                # Now start comparing core-bondedness similarities
                droplist=[]
                for m2 in range(len(ms)):
                    droplist.append([])
                    if m1 != m2:
                        ibuff=-1
                        for at2 in range(len(atoms[m2])):
                            if atoms[m2][at2] in matches[m1][atoms[m1][at1]][m2]:
                                ibuff+=1
                                if len(m2list[m2][ibuff]) != len(m1list):
                                    droplist[-1].append(atoms[m2][at2])
                                    continue
                                chk=0
                                for c1 in range(len(m1list)):
                                    mcomp=cores[m1][m1list[c1]][m2][0]  # this is the core atom in m2 that matches m1list[c1]
                                    if mcomp in m2list[m2][ibuff]:
                                        chk+=1
                                if chk != len(m1list):
                                    droplist[-1].append(atoms[m2][at2])
                                    continue
    
                #if debug:
                #    print(atoms[m1][at1],m1list)
                #    print(m2list)
                #    print('droplist',droplist,'\n')
    
                # Now rebuild matches, while removing the droplist items (matches[molecule][key] = list of list)
                tmplist = []
                for m2 in range(len(ms)):
                    tmplist.append([])
                    if m1 != m2:
                        for at2 in range(len(atoms[m2])):
                            if atoms[m2][at2] in matches[m1][atoms[m1][at1]][m2]:
                                if not (atoms[m2][at2] in droplist[m2]):
                                    tmplist[m2].append(atoms[m2][at2])
                # interpret tmplist for potential removal from matches
                empty=0
                chk=0
                for m2 in range(len(ms)):
                    if tmplist[m2] == []:
                        empty+=1
                    elif len(tmplist[m2]) == 1:
                        chk+=1
                    else:
                        pass
                if empty > 1:
                    matches[m1].pop(atoms[m1][at1])
                elif chk == len(ms)-1:
                    matches[m1].pop(atoms[m1][at1])
                    cores[m1][atoms[m1][at1]]=deepcopy(tmplist)
                else:
                    matches[m1][atoms[m1][at1]]=deepcopy(tmplist)
    
    # # debug
    if debug:
        print("BEFORE CORE CHK")
        print("")
        printListDict(matches)
        print("")
        print("")
        printListDict(cores)
        print("")
        print("")
    
    ## (iiE) Double check that new two (or more) core atoms do not have the same match. If they do,
    ## then move them out of cores and back into matches. Let next step (distances) try to fix things
    for m1 in range(len(ms)):
        m1keys=list(cores[m1].keys())
        for m2 in range(len(ms)):
            if m1 != m2:
                droplist=[]
                for k1 in range(len(m1keys)):
                    at=cores[m1][m1keys[k1]][m2][0]
                    for k2 in range(k1+1,len(m1keys)):
                        if at == cores[m1][m1keys[k2]][m2][0]:
                            # then we have a match that we need to fix
                            if not (m1keys[k1] in droplist):
                                droplist.append(m1keys[k1])
                            if not (m1keys[k2] in droplist):
                                droplist.append(m1keys[k2])
                if len(droplist) > 0:
                    for k in droplist:
                        matches[m1][k]=cores[m1][k]
                        cores[m1].pop(k)
                    m1keys=list(cores[m1].keys())
               
              
    
    # # debug
    if debug:
        print("AFTER CORE CHK, BUT BEFORE DIST CHECKS")
        print("")
        printListDict(matches)
        print("")
        print("")
        printListDict(cores)
        print("")
        print("")
    
    
    ## (iiF) RMSD analysis between remaining match atoms (last resort)
    ## *** !!! This requires the molecules to be ALIGNED to work correctly !!! ***
    ## To check that all molecules are suitably aligned, do RMSD of corrent CORE atoms
    ## (exit with error message if RMSD is > cutoff)
    
    
    #cutoff = 0.8    # disregard atom pairs with RMSDs above this value ## uncomment to hardcode this
    
       
    # For each core atom, compute an "average" xyz position, and then calc the RMSD for all
    # real atomic coordinates to this average "reference" position
    rmsd={}
    for m1 in range(1):   # treat first molecule as the "reference molecule"
        for k in cores[m1].keys():
            avg={"X":0.00,"Y":0.00,"Z":0.00}
            rmsd[k]=0.00
            # calc avg position first
            avg['X']+=xyzs[m1][k]["X"]
            avg['Y']+=xyzs[m1][k]["Y"]
            avg['Z']+=xyzs[m1][k]["Z"]
            for m2 in range(1,len(ms)):
                avg['X']+=xyzs[m2][cores[m1][k][m2][0]]["X"]
                avg['Y']+=xyzs[m2][cores[m1][k][m2][0]]["Y"]
                avg['Z']+=xyzs[m2][cores[m1][k][m2][0]]["Z"]
            avg['X']=avg['X']/float(len(ms))
            avg['Y']=avg['Y']/float(len(ms))
            avg['Z']=avg['Z']/float(len(ms))
    
            # then calc rmsd for all atoms pairs to this average position
            tmp=xyzs[m1][k]["X"]-avg['X']
            rmsd[k]+=tmp*tmp
            tmp=xyzs[m1][k]["Y"]-avg['Y']
            rmsd[k]+=tmp*tmp
            tmp=xyzs[m1][k]["Z"]-avg['Z']
            rmsd[k]+=tmp*tmp
            for m2 in range(1,len(ms)):
                tmp=xyzs[m2][cores[m1][k][m2][0]]["X"]-avg['X']
                rmsd[k]+=tmp*tmp
                tmp=xyzs[m2][cores[m1][k][m2][0]]["Y"]-avg['Y']
                rmsd[k]+=tmp*tmp
                tmp=xyzs[m2][cores[m1][k][m2][0]]["Z"]-avg['Z']
                rmsd[k]+=tmp*tmp
            rmsd[k]=np.sqrt(rmsd[k]/float(len(ms)))
    
    chk=0
    for k in rmsd.keys():
        if rmsd[k] > cutoff:
            chk+=1
            print("Large RMSD for CORE atom",k,"- RMSD >",cutoff,"("+str(round(rmsd[k],3))+").")
    if chk == len(cores[0].keys()) and len(cores[0].keys()) != 0:
        print("ERROR: Large RMSDs detected for ALL core atoms.\n     ",
              "This suggests that the molecules are not properly aligned!. Check and resubmit")
        quit()
    
    ## (iiG) Assuming the rmsd calculation has successfully completed and there were no alignment problems
    ## Move on to using distance to try to deduce if the remaining atoms pairs can be matched and moved into the core
    ## Calc distance between atom pairs and eliminate options above the cutoff value (cutoff+0.3)
    
    for m1 in range(len(ms)):
        for at1 in range(len(atoms[m1])):
            if atoms[m1][at1] in matches[m1].keys():
                droplist=[]
                for m2 in range(len(ms)):
                    droplist.append([])
                    if m1 != m2:
                        for at2 in range(len(atoms[m2])):
                            if atoms[m2][at2] in matches[m1][atoms[m1][at1]][m2]:
                                dist=0.00
                                tmp=xyzs[m1][atoms[m1][at1]]['X']-xyzs[m2][atoms[m2][at2]]['X']
                                dist+=tmp*tmp
                                tmp=xyzs[m1][atoms[m1][at1]]['Y']-xyzs[m2][atoms[m2][at2]]['Y']
                                dist+=tmp*tmp
                                tmp=xyzs[m1][atoms[m1][at1]]['Z']-xyzs[m2][atoms[m2][at2]]['Z']
                                dist+=tmp*tmp
                                dist=np.sqrt(dist)
                                if dist > cutoff+0.3:
                                    droplist[-1].append(atoms[m2][at2])
    
                # Now rebuild matches, while removing the droplist items (matches[molecule][key] = list of list)
                tmplist = []
                for m2 in range(len(ms)):
                    tmplist.append([])
                    if m1 != m2:
                        for at2 in range(len(atoms[m2])):
                            if atoms[m2][at2] in matches[m1][atoms[m1][at1]][m2]:
                                if not (atoms[m2][at2] in droplist[m2]):
                                    tmplist[m2].append(atoms[m2][at2])
                # interpret tmplist for potential removal from matches
                empty=0
                chk=0
                for m2 in range(len(ms)):
                    if tmplist[m2] == []:
                        empty+=1
                    elif len(tmplist[m2]) == 1:
                        chk+=1
                    else:
                        pass
                if empty > 1:
                    matches[m1].pop(atoms[m1][at1])
                elif chk == len(ms)-1:
                    matches[m1].pop(atoms[m1][at1])
                    cores[m1][atoms[m1][at1]]=deepcopy(tmplist)
                else:
                    matches[m1][atoms[m1][at1]]=deepcopy(tmplist)
    
    
    # debug
    if debug:
        print("AFTER DIST CHKS")
        print("MATCHES")
        printListDict(matches)
        print("")
        print("CORE ATOMS")
        printListDict(cores)
        print("")
        print("")
    
    
    ## (iiH) Try a variety of checks to ensure that all atoms in matches are either (i) classified
    ## as a core atom or discarded as a future fragment atom
    def emptyMatches(mollist):
        """ check if the matches structure is empty for all molecules in mollist """
        echk=0
        for mol in range(len(mollist)):
            if len(matches[m]) != 0:
                echk+=1
        return echk
    
    # skip checks if matches is already empty
    chk=emptyMatches(ms)
    if chk == 0:
        pass
    else:
        # for remaining atoms, search their matches for atoms that might have been identified as a core atom already;
        # if found, remove those atoms from matches[m][at]; classify matches[m][at] as core or discarded atom if possible
        for m1 in range(len(ms)):
            for at1 in range(len(atoms[m1])):
                if atoms[m1][at1] in matches[m1].keys():
                    droplist=[]
                    for m2 in range(len(ms)):
                        droplist.append([])
                        if m1 != m2:
                            # append the index of at1 in matches[m1][atoms[m1][at1]][m2]
                            for at2 in range(len(matches[m1][atoms[m1][at1]][m2])):
                                if matches[m1][atoms[m1][at1]][m2][at2] in cores[m2].keys():
                                    droplist[m2].append(at2)
                            # remove droplist
                            droplist[m2].sort(reverse=True)
                            for d in droplist[m2]:
                                matches[m1][atoms[m1][at1]][m2].pop(d)
                    # check for easy to remove match lines
                    empty=0
                    chk=0
                    for m2 in range(len(ms)):
                        if matches[m1][atoms[m1][at1]][m2] == []:
                            empty+=1
                        elif len(tmplist[m2]) == 1:
                            chk+=1
                        else:
                            pass
                    if empty > 1:
                        matches[m1].pop(atoms[m1][at1])
                    elif chk == len(ms)-1:
                        cores[m1][atoms[m1][at1]]=matches[m1].pop(atoms[m1][at1])
                    else:
                        pass
        
    
    ## (iiI) Make sure all atoms are out of matches - exit with error if not
    chk=emptyMatches(ms)
    if chk > 0:
        print("ERROR: Atoms have not been classified as core or fragment atoms!\n    ",
              "User input is required to continue running. Please specify the matching atom(s)",
              "for each atom in each molecule\n")
        for m1 in range(len(ms)):
            if len(matches[m1]) > 0:
                print("### Molecule",m1,"("+ms[m1]+")")
                print(matches[m1])
            for at1 in range(len(atoms[m1])):
                if atoms[m1][at1] in matches[m1].keys():
                    print("Molecule",m1,"("+ms[m1]+") - ATOM",atoms[m1][at1],"matches with:")
                    for m2 in range(len(ms)):
                        if m1 != m2:
                            if len(matches[m1][atoms[m1][at1]][m2]) > 1:
                                useratom=input("Molecule "+str(m2)+" ("+ms[m2]+"),("+str(matches[m1][atoms[m1][at1]][m2])+"): ")
                                if useratom == 'exit' or useratom == 'EXIT' or useratom == 'Exit':
                                    print("Exiting...")
                                    quit()
                                matches[m1][atoms[m1][at1]][m2]=[useratom]
                    print("")
            print("\n")
        printListDict(matches)
        quit()
    
    ## (iiJ) Make sure all core lists have the same number of atoms: Correct inconsistencies
    chkcores=[]
    for m1 in range(len(ms)):
        chkcores.append([])
        m1cores=list(cores[m1].keys())
        for m2 in range(len(ms)):
            chkcores[m1].append([])
            if m1 != m2:
                m2cores=list(cores[m2].keys())
                # generate list of m1 core matches from cores[m2] key-values
                m2matches=[]
                for at2 in m2cores:
                    m2matches.append(cores[m2][at2][m1][0])
                for at1 in m1cores:
                    if not (at1 in m2matches):
                        chkcores[m1][m2].append(at1)
                        # at1 = key in m1
                        line=deepcopy(cores[m1][at1])
                        line[m1]=[at1]
                        cores[m2][line[m2][0]]=line
                        cores[m2][line[m2][0]][m2]=[]
    
    # debug
    if debug:
        print("AFTER CORE CHKS")
        print("MATCHES")
        printListDict(matches)
        print("")
        print("CORE ATOMS")
        printListDict(cores)
        print("")
        print("")
    
    
    ## (5) Figure out the number of fragments in each molecule, pair them together between molecules,
    ## and remove redundancies
    fragbonds=deepcopy(bonds)         ## no H bonds, no core-core bonds
    for m1 in range(len(ms)):
        for at1 in range(len(atoms[m1])):
            for at2 in range(at1+1,len(atoms[m1])):
                if (atoms[m1][at1] in cores[m1].keys()) and (atoms[m1][at2] in cores[m1].keys()):
                        fragbonds[m1][at1][at2]=0
                        fragbonds[m1][at2][at1]=0
                if (xyzs[m1][atoms[m1][at1]]['Elem'] == 'H') or (xyzs[m1][atoms[m1][at2]]['Elem'] == 'H'):
                        fragbonds[m1][at1][at2]=0
                        fragbonds[m1][at2][at1]=0
        ##print("Molecule",m1)
        ##for at1 in range(len(atoms[m1])):
        ##    for at2 in range(len(atoms[m1])):
        ##        if at2 < at1:
        ##            print("  ",end='')
        ##        else:
        ##            print(" "+str(fragbonds[m1][at1][at2]),end='')
        ##    print("")
    
    # figure out the number of core-to-(not core) attachments
    nsites=[]
    Aatoms=[]
    for m in range(len(ms)):
        Aatoms.append({})
        chk=0
        for at1 in range(len(atoms[m])):
            if atoms[m][at1] in cores[m].keys():
                for at2 in range(len(atoms[m])):
                    if xyzs[m][atoms[m][at2]]['Elem'] != 'H':
                        if fragbonds[m][at1][at2] == 1:
                            chk+=1
                            if atoms[m][at1] in Aatoms[m].keys():
                                Aatoms[m][atoms[m][at1]].append(atoms[m][at2])
                            else:
                                Aatoms[m][atoms[m][at1]]=[atoms[m][at2]]
        nsites.append(chk)
    
    # initial nsites and Aatoms
    if debug:
        print("NSITES and AATOMS")
        print(nsites)
        print(Aatoms)
    
    
    # now take these sites and start expounding on them to list out all fragment atoms attached to them
    # (remember to look for rings (where one fragment branch joins another nsite/Aatom point)
    
    def GetIndex(atomname,molnum):
        """ given an atomname for a molecule number, return the atom index"""
        for atomindex in range(len(atoms[molnum])):
            if atoms[molnum][atomindex] == atomname:
                break
        return atomindex
    
    Fatoms=[]
    for m in range(len(ms)):
        Fatoms.append([])
        for k in Aatoms[m].keys():
            for fragatom in Aatoms[m][k]:
                Fatoms[m].append([])
                Fatoms[m][-1].append(fragatom)
                # find the atom index for fragatom
                at1 = GetIndex(fragatom,m)
                # enumerate all bonds to this atom (excluding core and H atoms)
                chk=1
                atlist=[]
                while chk > 0:
                    i=fragbonds[m][at1,:]
                    for j in range(len(i)):
                        if i[j] == 1:
                            if not (atoms[m][j] in cores[m].keys()) and not (atoms[m][j] in Fatoms[m][-1]):
                                atlist.append(j)
                    if atlist==[]:
                        chk=0
                    else:
                        at1=atlist.pop()
                        if not (atoms[m][at1] in Fatoms[m][-1]):
                            Fatoms[m][-1].append(atoms[m][at1])
        # check for redundant/identical fragments & merge if found
        droplist=[]   # list the indices to drop from Fatoms[m][-1]
        tsites=nsites[m]
        for k1 in range(tsites-1):
            for k2 in range(k1+1,tsites):
                chk=0
                tmpfrag=deepcopy(Fatoms[m][k2])
                for at1 in Fatoms[m][k1]:
                    for at2 in range(len(tmpfrag)):
                        if at1 == tmpfrag[at2]:
                            chk+=1
                            tmpfrag.pop(at2)
                            break
                if chk == len(Fatoms[m][k1]):
                    # then we have a match
                    if not (k2 in droplist):
                        nsites[m]=nsites[m]-1
                        droplist.append(k2)
        # merge found matches; modify: (i) Fatoms, (ii) Aatoms (as part of next routine)
        if len(droplist) > 0:
            tmp=[]
            droplist.sort(reverse=True) # to work highest index to lowest
            for i in droplist:
                tmp=Fatoms[m].pop(i)  # pop in reverse direction so we don't remove things we actually want
        # replace Aatom values with frag atom lists
        for frag in Fatoms[m]:
            droplist=[]
            for k in Aatoms[m].keys():
                for f in range(len(Aatoms[m][k])):
                    if Aatoms[m][k][f] in frag:
                        droplist.append(k)
            if len(droplist) > 1:
                # merge key_names and replace with new frag list
                newname=""
                for at in range(len(droplist)):
                    if at == 0:
                        newname+=droplist[at]
                    else:
                        newname+="+"+droplist[at]
    
                    if len(Aatoms[m][droplist[at]]) == 1:
                        Aatoms[m].pop(droplist[at]) # doesn't work if Aatoms[m][k] has more than one atom in it's nested list
                    else:
                        # try making a new Aatoms[m][at] list without the matching frag atom
                        for f in range(len(Aatoms[m][droplist[at]])):
                            if Aatoms[m][droplist[at]][f] in frag:
                                break
                        newlist=[]
                        for f2 in range(len(Aatoms[m][droplist[at]])):
                            if f != f2:
                                newlist.append(Aatoms[m][droplist[at]][f2])
                        Aatoms[m][droplist[at]]=newlist
                Aatoms[m][newname] = deepcopy(frag)
            elif len(droplist) == 1:
                Aatoms[m][droplist[0]] = deepcopy(frag)
            else: #len(droplist) < 1:
                print("ERROR: frag mismatch between Fatoms list and Aatoms match")
                quit()
    
    # Do some quick checks:
    for m in range(len(ms)):
        if nsites[m] != len(Aatoms[m]) or nsites[m] != len(Fatoms[m]):
            print("ERROR: mismatch in the number of sites for molecule",m)
            print("MOLECULE",m)
            print("NSITES[M]",nsites[m])
            print("AATOMS[M]",Aatoms[m])
            print("FATOMS[M]",Fatoms[m])
            quit()
    
    
    if debug:
        # nsites, Aatoms, and Fatoms after sorting is complete
        print("")
        print("NSITES:",nsites,'\n')
        print("AATOMS:_______________________________")
        printListDict(Aatoms)
        print("FATOMS:_______________________________")
        printListList(Fatoms)
        print("")
    
    ## (6) Identify the reference ligand (either hardcode it as the first molecule, or pick the molecules with the smallest number of fragment atoms
    
    #refnum=0  # uncomment if you want to hardcode this value
    minfragatoms=-1
    for m in range(len(ms)):
        sumchk=0
        for k in Aatoms[m].keys():
            sumchk+=len(Aatoms[m][k])
        if minfragatoms == -1:
            minfragatoms=sumchk
            refnum=m
        else:
            if sumchk < minfragatoms:
                minfragatoms=sumchk
                refnum=m
    
    if debug:
        print("REFNUM =",refnum,"\n\n")
    
    
    ## Make sure all ligands have the same number of sites
    chk=0
    for m in range(len(ms)):
        if nsites[m] != nsites[refnum]:
            chk+=1
    if chk != 0:
        print("ERROR: Not all molecules have the same number of sites defined! Check and resubmit")
        quit()
    
    ## Order the fragments so that everything is consistent when we do redundant checks
    # use the refnum molecule to decide what fragment is first
    Ftemplate=[]
    for k in Aatoms[refnum].keys():
        Ftemplate.append(k)
    
    # generate the expected keys based off what's in cores, then check to make sure everything is correct
    def chkMergedKey(keytocheck):
        """ Check if a key is a merged key (name+name). Split into a list if yes, return value if not """
        tmp=keytocheck.split('+')
        if len(tmp) > 1:
            keyname = tmp
        else:
            keyname=[keytocheck]
        return keyname
    
    Forder=[]
    for m1 in range(len(ms)):
        Forder.append([])
        for temp in Ftemplate:
            tt=chkMergedKey(temp)
            if len(tt) > 1: # then it is a merged key
                tmp=""
                for t3 in range(len(tt)):
                    if m1 == refnum:
                        i=tt[t3]
                    else:
                        i=cores[refnum][tt[t3]][m1][0]
                    if t3 == 0:
                        tmp+=i
                    else:
                        tmp+="+"+i
                Forder[m1].append(tmp)
    
            else:           # it is a single value key
                if m1 == refnum:
                    Forder[m1].append(tt[0])
                else:
                    Forder[m1].append(cores[refnum][tt[0]][m1][0])
        # make sure no two keys are the same
        chk = 0
        for f1 in range(len(Forder[m1])):
            for f2 in range(len(Forder[m1])):
                if f1 != f2:
                    if Forder[m1][f1] == Forder[m1][f2]:
                        chk+=1
        if chk > 0:
            print("ERROR: Non-unique keys found linking core to fragment atoms!")
            quit()
        # make sure the generated key actually matches what's in Aatoms[m1]
        # and account for order differences
        chk = 0
        tmp=list(Aatoms[m1].keys())
        for f in Forder[m1]:
            f1=chkMergedKey(f)
            if len(f1) == 1:
                for f2 in range(len(tmp)):
                    if tmp[f2] == f:
                        chk+=1
                        tmp.pop(f2)
                        break
            else:
                # consider different orders; if found, update Aatoms keyname
                for f2 in range(len(tmp)):
                    f3=chkMergedKey(tmp[f2])
                    if len(f3) == len(f1):
                        chk2=0
                        for f4 in f1:
                            for f5 in range(len(f3)):
                                if f4 == f3[f5]:
                                    chk2+=1
                                    f3.pop(f5)
                                    break
                        if chk2 == len(f1):
                            chk+=1
                            kk=tmp.pop(f2)
                            break
                # update Aatoms key
                if kk != f:
                    Aatoms[m1][f]=Aatoms[m1][kk]
                    Aatoms[m1].pop(kk)
        if chk != len(Forder[m1]):
            print("ERROR: Mismatch between generated and actual Aatom keys!")
            print("Molecule",m1)
            print("AATOMS",Aatoms[m1].keys())
            print("FORDER",Forder[m1])
            quit()
    
    if debug:
        print("FORDER:")
        printListList(Forder)
    
    ## (7) Now look for redundancies in the fragments in each molecule so that we only print out a single fragment per site
    # refnum fragments are added by default
    
    # shortest distance function
    def findShortestDistance(atname1,atname2,molnum):
        """ Explore bonds to find the minimum distance between two atoms """
        at1=GetIndex(atname1,molnum)
        at2=GetIndex(atname2,molnum)
    
        ## at1 = "anchor atom"; start here and tree down until we find the atom we're looking for (at2)
        chk = 0
        atlist=[at1]
        while True:
            chk+=1
            bdlist=[]
            for i in range(len(atlist)):
                ii=bonds[molnum][atlist[i],:]
                for j in range(len(ii)):
                    if ii[j] == 1 and xyzs[molnum][atoms[molnum][j]]['Elem'] != 'H':
                        if not (j in bdlist):
                            bdlist.append(j)
            if at2 in bdlist:
                mindist = chk
                break
            else:
                atlist=deepcopy(bdlist)
    
        return mindist
    
    # use a compare function:
    def fragCompare(ufidx,ufkey,qidx,qkey):
        """ Compare two fragments & try to figure out if they are the same or not.
            Return "TRUE" if they match, otherwise, return "FALSE". """
        uflist=Aatoms[ufidx][ufkey]
        qlist=Aatoms[qidx][qkey]
    
        # reject if list lenghts are different
        if len(uflist) != len(qlist):
            return False
    
        # reject if lists have atoms with different bonding patterns or unequal #s of atoms with same bonding patterns
        #   i. calc bonding patterns for each atom in uflist and qlist
        ufpatt=[]
        for at in uflist:
            ufpatt.append(getBonded(GetIndex(at,ufidx),ufidx))
        qpatt=[]
        for at in qlist:
            qpatt.append(getBonded(GetIndex(at,qidx),qidx))
    
        #   ii. find shortest distance to a single core atom (first core atom in key if its a merged key)
        atname1=ufkey.split('+')[0]
        for atname2 in range(len(uflist)):
            mindist=findShortestDistance(atname1,uflist[atname2],ufidx)
            ufpatt[atname2].append(mindist)
        atname1=qkey.split('+')[0]
        for atname2 in range(len(qlist)):
            mindist=findShortestDistance(atname1,qlist[atname2],qidx)
            qpatt[atname2].append(mindist)
            
        qtemp=deepcopy(qpatt)
        
        #   iii. now compare types
        amatch=0
        for at1 in range(len(uflist)):
            for at2 in range(len(qlist)):
                if ufpatt[at1][0] == qtemp[at2][0] and ufpatt[at1][2] == qtemp[at2][2]: # atom type matches
                    bdsum=0
                    for bd1 in ufpatt[at1][1]:
                        if bd1 in qtemp[at2][1]:
                            bdsum+=1
                            #remove bd1 atomtype from qtemp[at2][1] to prevent duplicities
                            for i in range(len(qtemp[at2][1])):
                                if qtemp[at2][1][i] == bd1:
                                    qtemp[at2][1].pop(i)
                                    break
                    if bdsum == len(ufpatt[at1][1]): # then they match
                        amatch+=1
                        qtemp[at2][0]='MATCHED' # prevents duplicities
                        break
        
        if amatch != len(uflist):
            return False
    
        # else: they match
        return True
    
    
    UFrag = []  # list of Unique Fragment ms indices
    for site in range(nsites[refnum]):
        UFrag.append([refnum])
    
    for site in range(nsites[refnum]):
        skip=[refnum] # refnum is included in UFrag by default, so we want to skip it in our comparisons
        ifrag=-1
        while len(skip) != len(ms):
            for frag in range(ifrag+1,len(UFrag[site])):
                newfrag=-1
                for m1 in range(len(ms)):
                    if not (m1 in skip):
                        #if COMPARISON says the groups are the same, then add m1 to skip (b/c it's frag is already represented)
                        #otherwise (else), add first non-match to UFrag and do the comparison over again
                        #repeat until all ms indices are in skip (b/c all unique fragments are in UFrag[site]
                        matched=fragCompare(UFrag[site][frag],Forder[UFrag[site][frag]][site],m1,Forder[m1][site])
                        if matched: # matched == True
                            skip.append(m1)
                        else:
                            if newfrag==-1:
                                newfrag = m1
                                UFrag[site].append(newfrag)
                                skip.append(newfrag)
            ifrag=frag
        
    if debug:
        print("UFRAG")
        print(UFrag,"\n")
    
    ## (8) Print A formatted file of information to pass onto CRN (separate script)
    ## This also allows you to modify things by hand 
    #fp=open('MCS_for_MSLD.txt','w')
    fp=open(mcsout,'w')
    fp.write('# Maximum Common Substructure Search for Multisite Lambda Dynamics (JV 2022)\n')
    fp.write('# %d molecules processed\n\n' % (len(ms)))
    # (1) Print nsubs info
    fp.write("NSUBS")
    for site in range(nsites[refnum]):
        fp.write(" %d" % (len(UFrag[site])))
    fp.write("\n\n")
    fp.write("REFLIG %s\n\n" %(ms[refnum]))
    # (2) Print Core Atoms
    fp.write('CORE \n')
    for mol in range(len(ms)):
        fp.write("%s" % (ms[mol]))
        # print heavy atoms in the core
        for k in cores[refnum].keys():
            if mol == refnum:
                at=k
            else:
                at=cores[refnum][k][mol][0]
            fp.write(" %s" % (at))
            # print attached hydrogen atoms too
            at1=GetIndex(at,mol)
            for at2 in range(len(atoms[mol])):
                if bonds[mol][at1][at2] == 1:
                    if xyzs[mol][atoms[mol][at2]]['Elem'] == 'H':
                        fp.write(" %s" % (atoms[mol][at2]))
        fp.write("\n")
    fp.write("\n")
    # (3) Print Anchor Atoms
    fp.write('ANCHOR ATOMS\n')
    for mol in range(len(ms)):
        fp.write("%s" % (ms[mol]))
        for site in range(nsites[refnum]):
            aatom=Forder[mol][site]
            chk=0
            for c in range(len(aatom)):
                if aatom[c] == '+':
                    # we have a "merged" Aatom
                    chk=1
            if chk == 1:
                fp.write(" %s" % ("DUM"))
            else:
                fp.write(" %s" % (aatom))
        fp.write("\n")
    fp.write("\n")
    # (4) Print Fragments
    for site in range(nsites[refnum]):
        fp.write('SITE '+str(site+1)+' FRAGMENTS\n')
        for uf in UFrag[site]:
            fp.write("%s" % (ms[uf]))
            for at in Aatoms[uf][Forder[uf][site]]:
                fp.write(" %s" % (at))
                # also print attached hydrogen atoms
                at1=GetIndex(at,uf)
                for at2 in range(len(atoms[uf])):
                    if bonds[uf][at1][at2] == 1:
                        if xyzs[uf][atoms[uf][at2]]['Elem'] == 'H':
                            fp.write(" %s" % (atoms[uf][at2]))
            fp.write("\n")
        fp.write("\n")
    
    fp.write("END  \n")
    fp.close()
    
    
    # finished
    return ms[refnum]








