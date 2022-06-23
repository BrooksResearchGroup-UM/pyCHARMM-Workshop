#! /usr/bin/env python

##
## Write out ligand parameters (non-charge FF parameters)
##

import numpy as np
import pandas as pd
from copy import deepcopy
import os

class PRM_Error(Exception):
    import sys
    sys.exit

def MsldPRM(outdir,cgenff,verbose=False,debug=False):
    """
    Using the information within the outdir rtf files, build outdir/full_ligand.prm 
    so that all ligand parameters are within a single file.
  
    If cgenff=True, the we expect missing parameters, otherwise
    write out warnings if explicit parameters are not available!

    Use verbose=True to get extra output
    Use debug=True to get LOTS of extra output
    """

    # load base information
    fp=open(outdir+'/nsubs','r')
    nsubs=[]
    for sub in fp.readline().split():
        nsubs.append(int(sub))
    nsites=len(nsubs)
    fp.close()

    if verbose:
        print("\nBuild full_ligand.prm for a series of ligands with")
        print(nsites,'sites (',nsubs,'substituents per site)\n')

    
    #################################################################
    ## Read the information from each *.rtf

    def readRTF(filename):
        """ Read and process a given rtf file """
        tmp={'NAME':'','MASS':[],'ATTYPE':{},'BOND':[],'IMPR':[]}
        masses=[]
        atoms=[]
        types={}
        bonds=[]
        imprs=[]

        fp=open(outdir+'/'+filename,'r')
        line=fp.readline()
        # get the name of the original file from the 2nd line
        line=fp.readline()
        tmp['NAME']=line.split()[-1][:-1]
        # get the rest of the rtf info
        while line:
            if line[0:4] == 'MASS':
                tmp['MASS'].append(line)
                masses.append(line)
            if line[0:4] == 'ATOM':
                lns=line.split()
                tmp['ATTYPE'][lns[1]]=lns[2]
                atoms.append(lns[1])
                types[lns[1]]=lns[2]
            if line[0:4] == 'BOND':
                tmp['BOND'].append(line.split()[1:])
                bonds.append(line.split()[1:])
            if line[0:4] == 'IMPR' or line[0:4] == '!IMP':
                tmp['IMPR'].append(line.split()[1:])
                imprs.append(line.split()[1:])
            if line[0:3] == 'END':
                break
            line=fp.readline()
        fp.close()
        return tmp,masses,atoms,types,bonds,imprs

    # initialize empty structures
    rtfinfo=[]
    allmasses=[]
    allatoms=[]
    alltypes={}
    allbonds=[]
    allimprs=[]

    # read core.rtf first
    (tmp,Cmass,Catoms,Ctypes,Cbonds,Cimprs)=readRTF('core.rtf')
    rtfinfo.append(tmp)
    allmasses=allmasses+Cmass
    allatoms=allatoms+Catoms
    alltypes.update(Ctypes)
    allbonds=allbonds+Cbonds
    allimprs=allimprs+Cimprs

    # read frags next
    for site in range(nsites):
        for frag in range(nsubs[site]):
            sname='site'+str(site+1)+'_sub'+str(frag+1)+'_pres.rtf'
            (tmp,Fmass,Fatoms,Ftypes,Fbonds,Fimprs)=readRTF(sname)
            rtfinfo.append(tmp)
            allmasses=allmasses+Fmass
            allatoms=allatoms+Fatoms
            alltypes.update(Ftypes)
            allbonds=allbonds+Fbonds
            allimprs=allimprs+Fimprs

    #if debug:
    #    print("MsldPRM structures after reading rtf files:")
    #    print("allmasses\n",allmasses)
    #    print("allatoms\n",allatoms)
    #    print("alltypes\n",alltypes)
    #    print("allbonds\n",allbonds)
    #    print("allimprs\n",allimprs)
    #    print()
    ##debug#


    #################################################################
    ## Read the information from each *.prm
    ## Open all the .prm files and read in the parameters
    ## Store the parameters in a dictionary structure for rapid checking
    duplprm={'BOND':[],'ANGL':[],'DIHE':[],'IMPR':[],'NONB':[]}
    prmbonds={}
    prmangs={}
    prmphis={}
    prmimpr={}
    prmnb={}
    for rtf in range(len(rtfinfo)):
        fp=open(rtfinfo[rtf]['NAME']+'.prm','r')
        line=fp.readline()
        while line:
            # look for BOND terms
            if line[0:4] == 'BOND':
                line=fp.readline()
                while (line != '\n') and (line[0:4] != 'ANGL'):
                    lns=line.split()
                    namef=lns[0]+' '+lns[1]
                    namer=lns[1]+' '+lns[0]
                    # check to make sure the bond isn't already in prmbonds
                    if not ((namef in prmbonds) or (namer in prmbonds)):
                        prmbonds[namef]=lns[2]+' '+lns[3]
                    else: # chk for same types but different parameters
                        if (namef in prmbonds) and (prmbonds[namef] != lns[2]+' '+lns[3]):
                            print(' ** Duplicate BOND with unique parameters ('+namef+') ** UNcomment as necessary') 
                            duplprm['BOND'].append('!'+namef+' '+lns[2]+' '+lns[3])
                        elif (namer in prmbonds) and (prmbonds[namer] != lns[2]+' '+lns[3]):
                            print(' ** Duplicate BOND with unique parameters ('+namer+') ** UNcomment as necessary') 
                            duplprm['BOND'].append('!'+namer+' '+lns[2]+' '+lns[3])
                        else:
                            pass
                    line=fp.readline()
            duplprm['BOND']=list(set(duplprm['BOND']))
            # look for ANGLE terms
            if line[0:4] == 'ANGL':
                line=fp.readline()
                while (line != '\n') and (line[0:4] != 'DIHE'):
                    lns=line.split()
                    # chk for UB terms
                    if len(lns) > 5:
                        if lns[5] == '!':
                            ub=False
                        else:
                            ub=True
                    else:
                        ub=False
                    namef=lns[0]+' '+lns[1]+' '+lns[2]
                    namer=lns[2]+' '+lns[1]+' '+lns[0]
                    # check to make sure the angle isn't already in prmangs
                    if not ((namef in prmangs) or (namer in prmangs)):
                        if ub:
                            prmangs[namef]=lns[3]+' '+lns[4]+' '+lns[5]+' '+lns[6]
                        else:
                            prmangs[namef]=lns[3]+' '+lns[4]
                    else: # chk for same types but different parameters
                        if ub:
                            if (namef in prmangs) and (prmangs[namef] != lns[3]+' '+lns[4]+' '+lns[5]+' '+lns[6]):
                                print(' ** Duplicate ANGLE with unique parameters ('+namef+') ** UNcomment as necessary') 
                                duplprm['ANGL'].append('!'+namef+' '+lns[3]+' '+lns[4]+' '+lns[5]+' '+lns[6])
                            elif (namer in prmangs) and (prmangs[namer] != lns[3]+' '+lns[4]+' '+lns[5]+' '+lns[6]):
                                print(' ** Duplicate ANGLE with unique parameters ('+namer+') ** UNcomment as necessary') 
                                duplprm['ANGL'].append('!'+namer+' '+lns[3]+' '+lns[4]+' '+lns[5]+' '+lns[6])
                            else:
                                pass
                        else:
                            if (namef in prmangs) and (prmangs[namef] != lns[3]+' '+lns[4]):
                                print(' ** Duplicate ANGLE with unique parameters ('+namef+') ** UNcomment as necessary') 
                                duplprm['ANGL'].append('!'+namef+' '+lns[3]+' '+lns[4])
                            elif (namer in prmangs) and (prmangs[namer] != lns[3]+' '+lns[4]):
                                print(' ** Duplicate ANGLE with unique parameters ('+namer+') ** UNcomment as necessary') 
                                duplprm['ANGL'].append('!'+namer+' '+lns[3]+' '+lns[4])
                            else:
                                pass
                    line=fp.readline()
            duplprm['ANGL']=list(set(duplprm['ANGL']))
            # look for DIHEDRAL terms
            if line[0:4] == 'DIHE':
                line=fp.readline()
                while (line != '\n') and (line[0:4] != 'IMPR'):
                    lns=line.split()
                    # add the multiplicity into the name
                    namef=lns[0]+' '+lns[1]+' '+lns[2]+' '+lns[3]+' '+lns[5]
                    namer=lns[3]+' '+lns[2]+' '+lns[1]+' '+lns[0]+' '+lns[5]
                    # check to make sure the phi isn't already in prmphis
                    if not ((namef in prmphis) or (namer in prmphis)):
                        prmphis[namef]=lns[4]+' '+lns[5]+' '+lns[6]
                    else: # chk for same types but different parameters
                        if (namef in prmphis) and (prmphis[namef] != lns[4]+' '+lns[5]+' '+lns[6]):
                            print(' ** Duplicate PHI with unique parameters ('+namef+') ** UNcomment as necessary') 
                            namef=lns[0]+' '+lns[1]+' '+lns[2]+' '+lns[3]+' '
                            duplprm['DIHE'].append('!'+namef+' '+lns[4]+' '+lns[5]+' '+lns[6])
                        elif (namer in prmphis) and (prmphis[namer] != lns[4]+' '+lns[5]+' '+lns[6]):
                            print(' ** Duplicate PHI with unique parameters ('+namer+') ** UNcomment as necessary') 
                            namer=lns[3]+' '+lns[2]+' '+lns[1]+' '+lns[0]+' '
                            duplprm['DIHE'].append('!'+namer+' '+lns[4]+' '+lns[5]+' '+lns[6])
                        else:
                            pass
                    line=fp.readline()
            duplprm['DIHE']=list(set(duplprm['DIHE']))
            # look for IMPROPER terms
            if line[0:4] == 'IMPR':
                line=fp.readline()
                while (line != '\n') and ((line[0:3] != 'NON') or (line[0:3] != 'END')):
                    lns=line.split()
                    # add the multiplicity into the name
                    namef=lns[0]+' '+lns[1]+' '+lns[2]+' '+lns[3]
                    namer=lns[3]+' '+lns[2]+' '+lns[1]+' '+lns[0]
                    # check to make sure the impr isn't already in prmimpr
                    if not ((namef in prmimpr) or (namer in prmimpr)):
                        prmimpr[namef]=lns[4]+' '+lns[5]+' '+lns[6]
                    else: # chk for same types but different parameters
                        if (namef in prmimpr) and (prmimpr[namef] != lns[4]+' '+lns[5]+' '+lns[6]):
                            print(' ** Duplicate IMPR with unique parameters ('+namef+') ** UNcomment as necessary') 
                            duplprm['IMPR'].append('!'+namef+' '+lns[4]+' '+lns[5]+' '+lns[6])
                        elif (namer in prmimpr) and (prmimpr[namer] != lns[4]+' '+lns[5]+' '+lns[6]):
                            print(' ** Duplicate IMPR PHI with unique parameters ('+namer+') ** UNcomment as necessary') 
                            duplprm['IMPR'].append('!'+namer+' '+lns[4]+' '+lns[5]+' '+lns[6])
                        else:
                            pass
                    line=fp.readline()
            duplprm['IMPR']=list(set(duplprm['IMPR']))
            # look for NONBONDED terms
            if not cgenff:
                if line[0:3] == 'NON':
                    line=fp.readline()
                    if line[0:3] == 'cut':
                        line=fp.readline() # skip the second nbond line
                    while (line != '\n') and (line[0:3] != 'END'):
                        lns=line.split()
                        if not (lns[0] in prmnb):
                            prmnb[lns[0]]=line
                        else:
                            if prmnb[lns[0]].split()[2:3] != lns[2:3]:
                                print(' ** Duplicate NONBONDED with unique parameters ('+namef+') ** UNcomment as necessary')
                                duplprm['NONB'].append('!'+line)
                        line=fp.readline() # skip the second nbond line
            duplprm['NONB']=list(set(duplprm['NONB']))
            # look for the END
            if line[0:3] == 'END':
                break
            line=fp.readline()
        fp.close()

    #if debug:
    #    print(len(prmbonds),"Total bond parameters")
    #    print(prmbonds,'\n')
    #    print(len(prmangs),"Total angle parameters")
    #    print(prmangs,'\n')
    #    print(len(prmphis),"Total dihedral parameters")
    #    print(prmphis,'\n')
    #    print(len(prmimpr),"Total improper dihedral parameters")
    #    print(prmimpr,'\n')
    #    print(len(prmnb),"Total nonbonded parameters")
    #    print(prmnb,'\n')
    #    print()
    #    print('Duplicate parameters:')
    #    print(duplprm)
    #    print()
    ##debug#


    #################################################################
    ## Build a Bond Matrix and Use it to Construct Connection Trees

    # Make a bond matrix
    bondmatx=pd.DataFrame(np.zeros((len(allatoms),len(allatoms))),columns=allatoms,index=allatoms,dtype=int)
    for bd in allbonds:
        bondmatx.loc[bd[0]][bd[1]]=1
        bondmatx.loc[bd[1]][bd[0]]=1

    if debug:
        print(bondmatx)

    # initialize new structures to be written out (lists keep the order, dicts for redundancy chks)
    newbondL=[]
    newbondD={}
    newangL=[]
    newangD={}
    newphiL=[]
    newphiD={}
    newimprL=[]
    newimprD={}
    newnbL=[]
    newnbD={}

    # build a tree for each atom and build parameter lists
    for at1 in allatoms:
        tree={}
        tree[at1]=list(bondmatx[bondmatx.loc[at1][:]==1].index)
        for at2 in tree[at1]:
            tree[at2]=list(bondmatx[(bondmatx.loc[at2][:]==1) & (bondmatx.loc[at2][:].index!=at1)].index)
            for at3 in tree[at2]:
                tree[at3]=list(bondmatx[(bondmatx.loc[at3][:]==1) & (bondmatx.loc[at3][:].index!=at2)].index)
        ##if debug:
        #if debug and at1 == allatoms[0]:
        #    print('Bond Tree for',at1)
        #    print(at1,tree[at1])
        #    for at2 in tree[at1]:
        #        print(at2,tree[at2])
        #        for at3 in tree[at2]:
        #            print(at3,tree[at3])
        ##debug#

        # cycle through the tree to identify bond, angle, and dihedral parameters
        for at2 in tree[at1]:
            # bond parameters
            namef=alltypes[at1]+' '+alltypes[at2]
            namer=alltypes[at2]+' '+alltypes[at1]
            if (namef in prmbonds) and not (namef in newbondD):
                newbondL.append(namef)
                newbondD[namef]=prmbonds[namef]
            elif (namer in prmbonds) and not (namer in newbondD):
                newbondL.append(namer)
                newbondD[namer]=prmbonds[namer]
            else:
                if not ((namef in prmbonds) or (namer in prmbonds)):
                    if not cgenff:
                        print("No BOND parameter found for "+at1+'-'+at2+' ('+namef+')')
                    else:
                        if verbose:
                            print("No BOND parameter found for "+at1+'-'+at2+' ('+namef+')')

            # nonbond parameters
            if (alltypes[at2] in prmnb) and not (alltypes[at2] in newnbD):
                newnbL.append(alltypes[at2])
                newnbD[alltypes[at2]]=prmnb[alltypes[at2]]
            else:
                pass # - no warning for missing nonbonded parameters

            # angle parameters
            for at3 in tree[at2]:
                namef=alltypes[at1]+' '+alltypes[at2]+' '+alltypes[at3]
                namer=alltypes[at3]+' '+alltypes[at2]+' '+alltypes[at1]
                if (namef in prmangs) and not (namef in newangD):
                    newangL.append(namef)
                    newangD[namef]=prmangs[namef]
                elif (namer in prmangs) and not (namer in newangD):
                    newangL.append(namer)
                    newangD[namer]=prmangs[namer]
                else:
                    if not ((namef in prmangs) or (namer in prmangs)):
                        if not cgenff:
                            print("No ANGLE parameter found for "+at1+'-'+at2+'-'+at3+' ('+namef+')')
                        else:
                            if verbose:
                                print("No ANGLE parameter found for "+at1+'-'+at2+'-'+at3+' ('+namef+')')

                # dihedral parameters
                for at4 in tree[at3]:
                    maxmultiplicity=6
                    chk=0
                    for m in range(1,maxmultiplicity+1): # check different multiplicities (cgenff goes up to 6)
                        namef=alltypes[at1]+' '+alltypes[at2]+' '+alltypes[at3]+' '+alltypes[at4]+' '+str(m)
                        namer=alltypes[at4]+' '+alltypes[at3]+' '+alltypes[at2]+' '+alltypes[at1]+' '+str(m)
                        if (namef in prmphis) and not (namef in newphiD):
                            newphiL.append(namef)
                            newphiD[namef]=prmphis[namef]
                        elif (namer in prmphis) and not (namer in newphiD):
                            newphiL.append(namer)
                            newphiD[namer]=prmphis[namer]
                        else:
                            if not ((namef in prmphis) or (namer in prmphis)):
                                chk+=1
                    if chk == maxmultiplicity:
                        if not cgenff:
                            print("No DIHEDRAL parameter found for "+at1+'-'+at2+'-'+at3+'-'+at4+' ('+namef[:-2]+')')
                        else:
                            if verbose:
                                print("No DIHEDRAL parameter found for "+at1+'-'+at2+'-'+at3+'-'+at4+' ('+namef[:-2]+')')

    #if debug:
    #    print("Bond Prms in Use:")
    #    #print(newbondL)
    #    print(newbondD)
    #    print("Angle Prms in Use:")
    #    #print(newangL)
    #    print(newangD)
    #    print("Dihedral Prms in Use:")
    #    #print(newphiL)
    #    print(newphiD)
    ##debug#

    # cycle through the impropers separately (since they can't be deduced from the tree)
    for impr in allimprs:
        namef=alltypes[impr[0]]+' '+alltypes[impr[1]]+' '+alltypes[impr[2]]+' '+alltypes[impr[3]]
        namer=alltypes[impr[3]]+' '+alltypes[impr[2]]+' '+alltypes[impr[1]]+' '+alltypes[impr[0]]
        if (namef in prmimpr) and not (namef in newimprD):
            newimprL.append(namef)
            newimprD[namef]=prmimpr[namef]
        elif (namer in prmimpr) and not (namer in newimprD):
            newimprL.append(namer)
            newimprD[namer]=prmimpr[namer]
        else:
            if not ((namef in prmimpr) or (namer in prmimpr)):
                if not cgenff:
                    print("No IMPROPER parameter found for "+impr[0]+'-'+impr[1]+'-'+impr[2]+'-'+impr[3]+' ('+namef+')')
                else:
                    if verbose:
                        print("No IMPROPER parameter found for "+impr[0]+'-'+impr[1]+'-'+impr[2]+'-'+impr[3]+' ('+namef+')')
    #if debug:
    #    print("Improper Prms in Use:")
    #    #print(newimprL)
    #    print(newimprD)


    #################################################################
    ## Write full_ligand.prm

    # Add in bonds and make a bond matrix
    fp=open(outdir+'/full_ligand.prm','w')
    fp.write("* MSLD ligand prm file generated with py_prep (JV,LC)\n* \n\n")

    fp.write("BONDS\n")
    for bond in newbondL:
        fp.write("%s %s\n" % (bond,newbondD[bond]))
    fp.write("\nANGLES\n")
    for angl in newangL:
        fp.write("%s %s\n" % (angl,newangD[angl]))
    fp.write("\nDIHEDRALS\n")
    for phi in newphiL:
        fp.write("%s %s\n" % (phi[:-2],newphiD[phi]))
    fp.write("\nIMPROPERS\n")
    for impr in newimprL:
        fp.write("%s %s\n" % (impr,newimprD[impr]))
    if len(newnbL) != 0:
        fp.write("\n")
        fp.write("NONBONDED nbxmod 5 atom cdiel switch vatom vdistance vswitch -\n")
        fp.write("cutnb 14.0 ctofnb 12.0 ctonnb 11.5 eps 1.0 e14fac 0.5  geom\n")
        for nb in newnbL:
            fp.write("%s" %(newnbD[nb])) # no new line character
            #fp.write("%s\n" %(newnbD[nb]))
    fp.write("\nEND")
    fp.close()

    return

