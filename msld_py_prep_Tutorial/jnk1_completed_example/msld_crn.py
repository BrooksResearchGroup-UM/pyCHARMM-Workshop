#! /usr/bin/env python

##
## Perform Charge Renormalization
##

import numpy as np
import pandas as pd
from copy import deepcopy
import os

class CRN_Error(Exception):
    import sys
    sys.exit

def MsldCRN(mcsout,outdir,inFrag,AnCore,ChkQChange=True,verbose=False,debug=False,ll=True):
    """
    Using the information within mcsout, as well as previously supplied mol2 and toppar
    files, perform charge renormalization to generate MSLD suitable force field parameters

    Use ChkQChange=True to check for charge perturbations between different substituents
    Use verbose=True to get extra output
    Use debug=True to get LOTS of extra output
    Use ll=True to build the "large_lig.pdb" file needed for easy solvation with Lg_Solvate.sh
    """

    #################################################################
    ## Read the information from mcsout
    mols=[]      # list of file names
    cores=[]     # list of lists of core atom names (same indexing as mols)
    Aatoms=[]    # list of lists of anchor atom names (same indexing as mols)
    frags=[]     # list of lists of file names for the fragments
    Fatoms=[]    # list of lists of the fragment atom names
    fp=open(mcsout,'r')
    line=fp.readline()
    nsites=0
    while line:
        if line[0:4] == 'NSUB':
            nsubs=[int(x) for x in line.split()[1:]]
            line=fp.readline()
        if line[0:4] == 'REFL':
            reflig=line.split()[1]
            line=fp.readline()
        if line[0:4] == 'CORE':
            line=fp.readline()
            while (line != '\n') and (line[0:4] != 'ANCH'): 
                lns=line.split()
                mols.append(lns[0])
                cores.append(lns[1:])
                line=fp.readline()
        if line[0:4] == 'ANCH':
            line=fp.readline()
            while (line != '\n') and (line[0:4] != 'SITE'): 
                lns=line.split()
                if lns[0] != mols[len(Aatoms)]:
                    raise CRN_Error('File names do not match between CORE and ANCHOR ATOMS - check and resubmit')
                Aatoms.append(lns[1:])
                line=fp.readline()
        if line[0:4] == 'SITE':
            nsites+=1
            frags.append([])
            Fatoms.append([])
            line=fp.readline()
            while (line != '\n') and (line[0:4] != 'SITE') and (line[0:3] != 'END'):
                lns=line.split()
                frags[nsites-1].append(lns[0])
                Fatoms[nsites-1].append(lns[1:])
                line=fp.readline()
        if line[0:3] == 'END':
            break
        line=fp.readline()
    fp.close()
    if nsites != len(nsubs):
        raise CRN_Error('Mismatch on the # of sites - check and resubmit')

    # transform cores into a DataFrame (DF) with reflig column headers and mols indices
    refnum=0
    for mol in range(len(mols)):
        if mols[mol] == reflig:
            refnum=mol
            break
    coreheader=cores[refnum]
    cores=pd.DataFrame(cores,columns=coreheader,index=mols,dtype=str)

    if debug:
        print("nsubs = ",nsubs,'\n')
        print("reflig = ",reflig,'\n')
        print("cores =  ("+str(cores.shape[0])+" total)\n",cores,'\n')
        print("anchor atoms =  ("+str(len(Aatoms))+" total)\n",Aatoms,'\n')
        for site in range(nsites):
            print("Site "+str(site+1)+" Fragments =  ("+str(len(frags[site]))+" total)")
            for frag in range(len(frags[site])):
                print(frags[site][frag],Fatoms[site][frag])
            print("\n")
    #debug#

    #################################################################
    ## Read the information from each *.rtf
    rtfinfo=[]   # list of dictionary infomation
    rtfvers1=36  # rtf version info
    rtfvers2=1

    for mol in mols:
        #initialize with blank lists
        #rtfinfo.append({'NAME':mol,'QNET':0,'MASS':[],'ATOM':[],'BOND':[],'IMPR':[],'ATTYPE':{},'ATQ':{}})
        rtfinfo.append({'NAME':mol,'QNET':0,'MASS':[],'ATTYPE':{},'ATQ':{},'BOND':[],'IMPR':[],'LP':[]})
    
        fp=open(mol+'.rtf','r')
        line=fp.readline()
        if line[0] == '*' and mol == refnum:
            while line:  # look for the charmm ff version # (assuming it directly follows the title lines)
                line=fp.readline()
                if line[0] != '*':
                    lns=line.split()
                    rtfvers1=lns[0]
                    rtfvers2=lns[1]
                    break
        while line:
            if line[0:4] == 'MASS':
                while line[0:4] == 'MASS':
                    lns=line.split()
                    rtfinfo[-1]['MASS'].append(lns[2:])
                    line=fp.readline()
            if line[0:4] == 'RESI':
                rtfinfo[-1]['QNET']=line.split()[2]
                line=fp.readline()
            if line[0:4] == 'ATOM':
                while line[0:4] == 'ATOM':
                    lns=line.split()
                    #rtfinfo[-1]['ATOM'].append(lns[1:4])
                    rtfinfo[-1]['ATTYPE'][lns[1]]=lns[2]
                    rtfinfo[-1]['ATQ'][lns[1]]=lns[3]
                    line=fp.readline()
            if line[0:4] == 'BOND':
                while line[0:4] == 'BOND':
                    lns=line.split()
                    # account for multiple bond entries on a single line
                    bdnum = len(lns) - 1
                    if bdnum == 2:
                        rtfinfo[-1]['BOND'].append(lns[1:3])
                    else:
                        for i in range(int(bdnum/2)):
                            j=2*i+1
                            rtfinfo[-1]['BOND'].append(lns[j:(j+2)])
                    line=fp.readline()
            if line[0:4] == 'IMPR':
                while line[0:4] == 'IMPR':
                    lns=line.split()
                    rtfinfo[-1]['IMPR'].append(lns[1:5])
                    line=fp.readline()
            if line[0:4] == 'LONE':
                while line[0:4] == 'LONE':
                    lns=line.split()
                    rtfinfo[-1]['LP'].append(lns[2:])
                    line=fp.readline()
            if line[0:3] == 'END':
                break
            line=fp.readline()

        fp.close()

    # Deduce if LPs go into the core or into a fragment
    coreLPs=[]
    for mol in range(len(mols)):
        coreLPs.append([])
        for lp in rtfinfo[mol]['LP']:
            #if lp[1] in coreheader:
            if lp[1] in list(cores.loc[mols[mol]][:].values):
                #add to the corelist DataFrame
                coreLPs[mol].append(lp[0])
            else: # it is in a fragment
                 # figure out which fragment lp[1] is in and then add lp[0] into the list of fragment atoms
                for site in range(nsites):
                    for frag in range(len(frags[site])):
                        if frags[site][frag] == rtfinfo[mol]['NAME']:
                            if lp[1] in Fatoms[site][frag]:
                                Fatoms[site][frag].append(lp[0])

    # Add Core LPs into cores
    coreLPs=pd.DataFrame(coreLPs,columns=coreLPs[refnum],index=mols,dtype=str)
    if not coreLPs.empty:
        for lp in range(coreLPs.shape[1]):
            cores[coreLPs.iloc[refnum][lp]]=coreLPs.iloc[:,lp]

    if debug:
        # print LP info
        if not coreLPs.empty:
            print("New cores with LP atoms:\n",cores,'\n')
        for mol in range(len(mols)):
            print(rtfinfo[mol]['NAME'],rtfinfo[mol]['QNET'])
            #for field in ['MASS','ATOM','BOND','IMPR','ATTYPE','ATQ']:
            for field in ['MASS','ATTYPE','ATQ','BOND','IMPR','LP']:
                print(" -- "+field+' -- ')
                print(rtfinfo[mol][field])
                print()
            print('----------')
    #debug#

    # for the reference ligand, identify H atoms bonded to heavy atoms (in the core)
    Hcore={}
    #for at in cores.loc[reflig][:]:
    for at in coreheader:
        Hcore[at]=[]
        for bd in rtfinfo[refnum]['BOND']:
            if at in bd:
                if at[0] != 'H' and bd[0][0] == 'H':
                    Hcore[at].append(bd[0])
                elif at[0] != 'H' and bd[1][0] == 'H':
                    Hcore[at].append(bd[1])
                # readable both ways (heavy atom <-> H atom)
                elif at[0] == 'H' and bd[0][0] != 'H':
                    Hcore[at].append(bd[0])
                elif at[0] == 'H' and bd[1][0] != 'H':
                    Hcore[at].append(bd[1])
                else:
                    pass
    if debug:
        print("H atom bonds in the core:")
        print(Hcore)
    #debug#


    #################################################################
    ## Account for inFrag and AnCore specifications (modifies cores)
    ## "AnCore" is really only useful to specify Aatoms to be left in
    ## the core - otherwise, they are (by default) moved onto the frags

    # check for correct formats for inFrag and AnCore
    crtfrm='['
    for site in range(nsites):
        crtfrm+='[]'
        if site != (nsites-1):
            crtfrm+=','
    crtfrm+=']'
    if len(inFrag) != nsites:
        # if it isn't right, print an error message and setup a default inFrag variable to continue with
        print("\ninFrag is not specified correctly!! It should look like this: "+crtfrm)
        print("A default (empty) value will be used, but this may not yield the results you want")
        inFrag=[]
        for site in range(nsites):
            inFrag.append([])
    if len(AnCore) != nsites:
        # if it isn't right, print an error message and setup a default inFrag variable to continue with
        print("\nAnCore is not specified correctly!! It should look like this: "+crtfrm)
        print("A default (empty) value will be used, but this may not yield the results you want")
        AnCore=[]
        for site in range(nsites):
            AnCore.append([])
    # create list of atoms to move out of the core (add in Aatoms by default)
    droplist=[]
    for at in range(len(Aatoms[refnum])):
        droplist.append([])
        if Aatoms[refnum][at] != 'DUM':
            droplist[-1].append(Aatoms[refnum][at])
    # add in any inFrag atoms 
    for site in range(nsites):
        if len(inFrag[site]) > 0:
            for at in inFrag[site]:
                if not (at in droplist[site]):
                    droplist[site].append(at)
    # remove any AnCore atoms
    for site in range(nsites):
        if len(AnCore[site]) > 0:
            for at in AnCore[site]:
                if at in droplist[site]:
                    droplist[site].remove(at)

    # add in any H atoms bonded to droplist atoms
    tmp=deepcopy(droplist)
    for site in range(nsites):
        for at in tmp[site]:
            if at[0:1] == 'H':
                pass # don't check for atoms bonded to a H
            else:
                if len(Hcore[at]) > 0:
                    for hat in Hcore[at]:
                        if hat in coreheader and not (hat in droplist[site]):
                            droplist[site].append(hat)

    # modify existing variable structures
    flatlist=[at for site in droplist for at in site]
    drops=cores[flatlist]               # DF of droplist atoms extracted from cores
    cores=cores.drop(columns=flatlist)  # cores is modified to remove droplist atoms
    coreheader=list(cores.columns.values)
    for site in range(nsites):
        for frag in range(len(frags[site])):
            for at in droplist[site]:
                Fatoms[site][frag].append(drops.loc[frags[site][frag]][at])

    if debug:
        print('\nlist of atoms to drop from the core:\n',droplist,'\n')
        print('atoms cut from the cores DF:\n',drops,'\n')
        print('modified fragments are:')
        for site in range(nsites):
            print("Site "+str(site+1)+" Fragments =  ("+str(len(frags[site]))+" total)")
            for frag in range(len(frags[site])):
                print(frags[site][frag],Fatoms[site][frag])
            print("\n")
        print("\n")
    #debug#


    #################################################################
    ## Perform Charge Renormalization (CRN)

    offset=0.000001  # the amount by which charges are renormalized
    dec=6            # the precision of final (written) charges
    Qcut=5.0         # provide a warning for large charge diffs in the core

    # check for charge perturbations
    qnet=float(rtfinfo[refnum]['QNET'])
    delQ=[]    # list of charge change boolean
    for mol in range(len(mols)):
        delQ.append(False) # Charge Change boolean
        if float(rtfinfo[mol]['QNET']) != float(qnet):
            print("We have a charge change between the reflig ("+mols[refnum]+' and '+mols[mol]+').',\
                  "ref_qnet = "+str(qnet),"; mol_qnet = "+rtfinfo[mol]['QNET'])
            if not ChkQChange:
                print("    ***  Recommended to turn ON the ChkQChange option !!  ***")
            delQ[mol]=True

    #(1) Gather and Average Core Charges 
    Qcore=pd.DataFrame(np.zeros((cores.shape[0]+2,cores.shape[1])),columns=coreheader,index=mols+['mean','stdev'],dtype=float)
    for mol in range(len(mols)):
        for at in range(len(cores.iloc[mol][:])):
            Qcore.iloc[mol][at]=rtfinfo[mol]['ATQ'][cores.iloc[mol][at]]
    for col in Qcore:
        Qcore.loc['mean',col]=float(str(np.around(np.average(Qcore.loc[mols,col]),decimals=dec)))
        Qcore.loc['stdev',col]=float(str(np.around(np.std(Qcore.loc[mols,col]),decimals=dec)))
    if debug:
        print(Qcore)
        print()
    #debug#

    #(2) Gather Fragment totals & CRN
    Qfrag=[]
    Qint=[]
    siteavg=[]
    if verbose:
        print("\nFragment Charge Differences Following Charge ReNormalization:")
    for site in range(nsites):
        sitesum=[]
        Qfrag.append([])
        for frag in range(len(frags[site])):
            # collect all the charges into Qfrag (list of lists of Series of frag atom charges)
            Qfrag[site].append(pd.Series(np.zeros(len(Fatoms[site][frag])),index=Fatoms[site][frag]))
            for mol in range(len(mols)):
                if mols[mol] == frags[site][frag]:
                    break
            for at in Fatoms[site][frag]:
                Qfrag[site][frag][at]=rtfinfo[mol]['ATQ'][at]
            # get the total charge of each fragment at each site, the mean, and start to renormalize charge
            sitesum.append(float(str(np.around(Qfrag[site][frag].sum(),decimals=dec))))

        if ChkQChange and verbose:
            print("Charge Change Chk:","old sitesum",sitesum)
        # check for charge changes; if found - remove the nearest int charge and crn
        QQ=[0 for q in sitesum]
        if ChkQChange:
            orig_sitesum=deepcopy(sitesum)
            molQs = [rtfinfo[mol]['QNET'] == rtfinfo[refnum]['QNET'] for mol in range(len(mols))]
            if not all(molQs):
                QQ=np.rint(sitesum) # round to nearest integer with numpy
                QQ=[int(q) for q in QQ]
                sitesum=[sitesum[q]-float(QQ[q]) for q in range(len(sitesum))]
                # get majority Qint
                Qint.append(float(str(np.around(np.rint(np.asarray(QQ).mean()),decimals=2))))
        if ChkQChange and verbose:
            print("Charge Change Chk:","Qs round to these integers",QQ)
            print("Charge Change Chk:","new sitesum",sitesum)
            print(" ")

        # avg the sitesums, and progressively CRN each fragment atom
        siteavg.append(float(str(np.around(np.asarray(sitesum).mean(),decimals=dec))))
        if verbose:
            print("Site "+str(site+1)+" Q(avg) = "+str(siteavg[site]))
        for frag in range(len(frags[site])):
            qdiff=float(str(np.around(siteavg[site]-sitesum[frag],decimals=dec)))
            foffset=offset
            if qdiff < 0.0:
                foffset=foffset*-1.0
            nsteps=int(np.around(qdiff/foffset,decimals=0)) # should be a positive int
            atom=0
            for step in range(nsteps):
                # cycle back to the beginning of the molecule
                if atom == len(Fatoms[site][frag]):
                    atom=0
                # if atom name starts with LP = skip it
                if Fatoms[site][frag][atom][0:2] == 'LP':
                    while Fatoms[site][frag][atom][0:2] == 'LP': # account for many LPs next to each other
                        atom+=1
                        if atom == len(Fatoms[site][frag]):
                            atom=0
                # add the offset to a non-LP atom
                Qfrag[site][frag][atom]+=foffset
                atom+=1
            # check that total frag charge matches siteavg[site]
            if ChkQChange:
                qchk=float(str(np.around(Qfrag[site][frag].sum()-float(QQ[frag]),decimals=dec)))
            else:
                qchk=float(str(np.around(Qfrag[site][frag].sum(),decimals=dec)))
            #print(frags[site][frag],qchk,str(qchk),str(qchk)=='-0.0')
            if str(qchk) == '-0.0':    
                qchk = float(0.0)   
            #print(frags[site][frag],qchk,str(qchk))
            if str(qchk) != str(siteavg[site]):
                raise CRN_Error('Error for fragment CRN for site '+str(site+1)+', frag = '+frags[site][frag]+
                ". Total charge ("+str(qchk)+") doesn't match the site average ("+str(siteavg[site])+") ")
            else:
                if verbose:
                    if sitesum[frag] == 0.0:
                        pdiff=float(str(np.around(qdiff*100.0,decimals=2)))
                    else:
                        pdiff=float(str(np.around(qdiff/sitesum[frag]*100.0,decimals=2)))
                    print("  "+frags[site][frag]+" Q(diff from avg) = "+str(qdiff)+
                          " Q(orig) = ",end='')
                    if ChkQChange:
                        print(str(orig_sitesum[frag]),end='')
                    else:
                        print(str(sitesum[frag]),end='')
                    print(" (a "+str(pdiff)+"% diff from orig charges)",end='')
                    if abs(pdiff) > Qcut:
                        print(" ** CHECK")
                    else:
                        print("")
        if verbose:
            print("")


    #if debug:
    #    for site in range(nsites):
    #    #for site in range(1):
    #        print('\nCharges for fragments at site '+str(site+1)+':')
    #        for frag in range(len(frags[site])):
    #        #for frag in range(1):
    #            print(frags[site][frag])
    #            print(Qfrag[site][frag])
    #            print()
    #        print()
    ##debug#


    #(3) CRN Core charges to neutralize ligands
    # extract current core (mean) charges
    QQ=pd.Series(Qcore.loc['mean'][:])
    # figure out the difference between the reflig's charge and current site charges
    sitesum=float(str(np.around(np.sum(siteavg),decimals=dec)))  # sum of all site charges
    intsum=float(str(np.around(np.sum(Qint),decimals=dec)))
    QQsum=float(str(np.around(QQ.sum(),decimals=dec)))      # sum of core (mean) charges
    tsum=float(str(np.around(QQsum+sitesum+intsum,decimals=dec)))  # total charge of msld ligand
    qdiff=float(str(np.around(float(rtfinfo[refnum]['QNET'])-tsum,decimals=dec)))

    coffset=offset
    if qdiff < 0.0:
        coffset=coffset*-1.0
    nsteps=int(np.around(qdiff/coffset,decimals=0)) # should be a positive int

    if debug: 
        print("ideal total charge = ",rtfinfo[refnum]['QNET'])
        print("site averages = ",siteavg)
        print("sitesum = ",sitesum)
        print("intsum = ",intsum)
        print("QQsum = ",QQsum)
        print("tsum = ",tsum)
        print("qdiff = ",qdiff)
        print("coffset = ",coffset)
        print("nsteps = ",nsteps)
        print()
    #debug#

    if debug:
        QQold=QQ.copy()
    ##debug#

    # do charge renormalization
    atom=0
    for step in range(nsteps):
        # cycle back to the beginning of the molecule
        if atom == QQ.shape[0]:
            atom=0
        # if atom name starts with LP = skip it
        if QQ.index[atom] == 'LP':
            while QQ.index[atom] == 'LP':
                atom+=1
                if atom == QQ.shape[0]:
                    atom=0
        # add the offset to a non-LP atom
        QQ[atom]+=coffset
        atom+=1
    # check that total core charges match the neg of sitesum (+qnet accounts for charged molecules)
    qchk=float(str(np.around(-1*QQ.sum(),decimals=dec)))
    sitechk=float(str(np.around(sitesum+intsum-qnet,decimals=dec)))
    if debug: 
        print("QCHK",qchk) 
        print("SiteCHK",sitechk) 
        print("Qnet",qnet) 
    if str(qchk) == '-0.0':
        qchk = 0.0
    if str(qchk) != str(sitechk):
        raise CRN_Error("Error for core CRN. Total charge ("+str(qchk*-1)+
                        ") doesn't neutralize site charge sums ("+str(sitesum)+") ")
    else:
        if verbose:
            print("Core Q(sum) = "+str(qchk*-1))
            if QQsum == 0.0:
                pdiff=0.0
            else:
                pdiff=float(str(np.around(qdiff/QQsum*100,decimals=2)))
            print("  Q(orig) = "+str(QQsum)+" (a "+str(pdiff)+"% diff from orig charges)",end='')
            if abs(pdiff) > Qcut:
                print(" ** CHECK")
            else:
                print("")
    if verbose:
        print("")

    if debug:
        print("CRN Core Charges:\nName,Old-Q,CRN-Q")
        for at in range(QQ.shape[0]):
            print(QQ.index[at],QQold[at],QQ[at])
    ##debug##
    
    #################################################################
    ## Check that the output directory exists - mkdir if not
    os.system('if [ ! -d '+outdir+' ]; then mkdir '+outdir+'; fi')


    #################################################################
    ## Read the information from each *.mol2 & write the *.pdb files
    ## and rename all atom names to a standardized format
    ## (only read the information you need:
    ##    - core from reflig/refnum
    ##    - frags from frags[site] files

    segid='LIG'     # segid for the ligand
    resname=segid   # residue name for the ligand

    def getElementSymbol(atomname):
        """ Read a Sybil atom type out of a mol2 file and return the atomic symbol """
        # assume that if there isn't a '.' that atomname is good as is (H, Cl, Br, etc)
        if atomname.find('.') == -1:
            newname=atomname
        # otherwise, figure out where the '.' is and use the letters before it
        else:
            newname=atomname.split('.')[0]
        return newname

    # read in the XYZ coords for the core
    coreXYZ={}
    coreEmt={} # store sybil atom types for better element deduction
    fp=open(reflig+'.mol2','r')
    line=fp.readline()
    while line:
        if line[0:13] == '@<TRIPOS>ATOM':
            line=fp.readline()
            while line[0:13] != '@<TRIPOS>BOND':
                lns=line.split()
                if lns[1] in coreheader:
                    coreXYZ[lns[1]]=lns[2:5]
                    coreEmt[lns[1]]=lns[5]
                line=fp.readline()
            break
        else:
            line=fp.readline()
    fp.close()
    # translate atom names into standardized format
    Ctrans={}  # translation dictionary between old (key) and new (value)
    if len(coreheader) > 200:
        raise CRN_Error("Can't handle more than 200 core atoms")
    pnum=[1,65,48]  # chr(65) == 'A' (ord('A') == 65); chr(48) = '0'
    for at in coreheader:
        # use function to get atomic symbol
        if at[0:2] == 'LP':
            newname='LP'
        else:
            newname=getElementSymbol(coreEmt[at])
        if len(newname) > 3:
            raise CRN_Error("Error in determining atomic symbol for atom "+at+" in the core")
        char=len(newname)
        if char == 1: # single char element
            if pnum[0] < 10:
                newat=newname+'00'+str(pnum[0])
            elif pnum[0] > 9 and pnum[0] < 100:
                newat=newname+'0'+str(pnum[0])
            else:
                newat=newname+str(pnum[0])
            pnum[0]+=1
        elif char == 2: # double char element
            if pnum[1] > 90:
                raise CRN_Error("Too many double char elements in core (max=26)")
            else:
                newat=newname+chr(pnum[2])+chr(pnum[1])
            pnum[1]+=1
        else:
            raise CRN_Error("Error in Core atom name translation")
        Ctrans[at]=newat
    # write core.pdb (atom order does not matter - but put H's after heavy atoms anyways)
    added=[]
    fp=open(outdir+'/core.pdb','w')
    row=0
    for at in coreheader:
        if at[0] == 'H':
            pass
        elif at[0] == 'L': # elif at[0:2] == 'LP':
            # no LP atoms added!
            pass
        else:
            if not (at in added):
                added.append(at)
                row+=1
                fp.write("ATOM  %5d %-4s %4s%5d    %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s\n" % (\
                         row,Ctrans[at],segid,1,float(coreXYZ[at][0]),float(coreXYZ[at][1]),float(coreXYZ[at][2]),1,0,segid))
                # check for bonded H's
                if len(Hcore[at]) > 0:
                    for at2 in Hcore[at]:
                        added.append(at2)
                        row+=1
                        fp.write("ATOM  %5d %-4s %4s%5d    %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s\n" % (\
                                 row,Ctrans[at2],segid,1,float(coreXYZ[at2][0]),float(coreXYZ[at2][1]),float(coreXYZ[at2][2]),1,0,segid))
    fp.write("TER\nEND")
    fp.close()

    # work on the fragments next
    Hfrag=[]
    Ftrans=[]
    for site in range(nsites):
        Hfrag.append([])
        Ftrans.append([])
        for frag in range(len(frags[site])):
            # read in the fragment XYZ coords
            fragXYZ={}
            fragEmt={} # store sybil atom types for better element deduction
            fp=open(frags[site][frag]+'.mol2','r')
            line=fp.readline()
            while line:
                if line[0:13] == '@<TRIPOS>ATOM':
                    line=fp.readline()
                    while line[0:13] != '@<TRIPOS>BOND':
                        lns=line.split()
                        if lns[1] in Fatoms[site][frag]:
                            fragXYZ[lns[1]]=lns[2:5]
                            fragEmt[lns[1]]=lns[5]
                        line=fp.readline()
                    break
                else:
                    line=fp.readline()
            fp.close()
            # identify mol # that equals frag #
            for mol in range(len(mols)):
                if mols[mol] == frags[site][frag]:
                    break
            # identify H atoms bonded to heavy atoms (in each fragment)
            Hfrag[site].append({})
            for at in Fatoms[site][frag]:
                Hfrag[site][frag][at]=[]
                for bd in rtfinfo[mol]['BOND']:
                    if at in bd:
                        if at[0] != 'H' and bd[0][0] == 'H':
                            Hfrag[site][frag][at].append(bd[0])
                        elif at[0] != 'H' and bd[1][0] == 'H':
                            Hfrag[site][frag][at].append(bd[1])
                        # readable both ways (heavy atom <-> H atom)
                        elif at[0] == 'H' and bd[0][0] != 'H':
                            Hfrag[site][frag][at].append(bd[0])
                        elif at[0] == 'H' and bd[1][0] != 'H':
                            Hfrag[site][frag][at].append(bd[1])
                        else:
                            pass
            #if debug:
            #    if frag == 0 and site == 0:
            #        print(frags[site][frag]+" XYZ coordinates")
            #        print(fragXYZ)
            #        print("H atom bonds in "+frags[site][frag])
            #        print(Hfrag[site][frag])
            #        print()
            ##debug#

            # translate atom names into standardized format
            Ftrans[site].append({})
            for at in Fatoms[site][frag]:
                # use function to get atomic symbol
                if at[0:2] == 'LP':
                    newname='LP'
                else:
                    newname=getElementSymbol(fragEmt[at])
                if len(newname) > 3:
                    raise CRN_Error("Error in determining atomic symbol for atom "+at+" in "+frags[site][frag])
                char=len(newname)
                if char == 1: # single char element
                    if pnum[0] > 999:
                        raise CRN_Error("Too many single char elements in ligand (max=1000)")
                    if pnum[0] < 10:
                        newat=newname+'00'+str(pnum[0])
                    elif pnum[0] > 9 and pnum[0] < 100:
                        newat=newname+'0'+str(pnum[0])
                    else:
                        newat=newname+str(pnum[0])
                    pnum[0]+=1
                elif char == 2: # double char element
                    if pnum[1] > 90:
                        #reset pnum[1] back to 'A' and increment pnum[2] by 1
                        pnum[1]=65
                        pnum[2]+=1
                    if pnum[2] > 57:
                        raise CRN_Error("Too many double char elements in ligand (max=260)")
                    else:
                        newat=newname+chr(pnum[2])+chr(pnum[1])
                    pnum[1]+=1
                else:
                    raise CRN_Error("Error in Frag atom name translation")
                Ftrans[site][frag][at]=newat

            # write the site[site]_sub[frag].pdb file
            added=[]
            fp=open(outdir+'/site'+str(site+1)+'_sub'+str(frag+1)+'_frag.pdb','w')
            row=0
            for at in Fatoms[site][frag]:
                if at[0] == 'H':
                    pass # H's follow the heavy atom
                elif at[0] == 'L': # elif at[0:2] == 'LP':
                    # no LP atoms added!
                    pass
                else:
                    if not (at in added):
                        added.append(at)
                        row+=1
                        fp.write("ATOM  %5d %-4s %4s%5d    %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s\n" % (\
                                 row,Ftrans[site][frag][at],segid,1,float(fragXYZ[at][0]),float(fragXYZ[at][1]),float(fragXYZ[at][2]),1,0,segid))
                        # check for bonded H's
                        if len(Hfrag[site][frag][at]) > 0:
                            for at2 in Hfrag[site][frag][at]:
                                added.append(at2)
                                row+=1
                                fp.write("ATOM  %5d %-4s %4s%5d    %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s\n" % (\
                                         row,Ftrans[site][frag][at2],segid,1,float(fragXYZ[at2][0]),float(fragXYZ[at2][1]),float(fragXYZ[at2][2]),1,0,segid))
            fp.write("TER\nEND")
            fp.close()
 
    # write translation file for old to new atom names
    filename='translation.txt'
    fp=open(filename,'w')
    fp.write("Original Atom Name -> New Atom Name\n")
    fp.write("CORE\n")
    for at in coreheader:
        fp.write("%s %s\n" % (at,Ctrans[at]))
    fp.write("\n")
    # also do it for the fragments
    for site in range(nsites):
        for frag in range(len(frags[site])):
            fp.write("SITE %d %s (from %s)\n" % (site+1,'site'+str(site+1)+'_sub'+str(frag+1),frags[site][frag]))
            for at in Fatoms[site][frag]:
                fp.write("%s %s\n" % (at,Ftrans[site][frag][at]))
            fp.write("\n")
        fp.write("\n")
    fp.close()
    if verbose:
        print("\nOld to New Atom Name Translations found in: "+filename)
        print("MSLD Files Written Into: "+outdir)
        print()

    # make a "large_lig.pdb" file for easier solv_prep preparation
    if ll:
        fp=open(outdir+'/large_lig.pdb','w')
        # add the core
        ip=open(outdir+'/core.pdb','r')
        for line in ip:
            if line[0:4] == 'ATOM':
                fp.write(line)
        ip.close()
        # add the largest fragment at each site
        for site in range(nsites):
            maxnum=0
            maxfrag=1
            for frag in range(len(frags[site])):
                if len(Fatoms[site][frag]) > maxnum:
                    print("maxnum from",'site'+str(site+1)+'_sub'+str(frag+1)+'_frag.pdb')
                    maxnum=len(Fatoms[site][frag])
                    maxfrag=frag
            ip=open(outdir+'/site'+str(site+1)+'_sub'+str(maxfrag+1)+'_frag.pdb','r')
            for line in ip:
                if line[0:4] == 'ATOM':
                    fp.write(line)
            ip.close()
        fp.write("TER\nEND")
        fp.close()


    #################################################################
    ## Write *.rtf 

    # write core.rtf
    fp=open(outdir+'/core.rtf','w')
    fp.write('* ligand core rtf file generated with msld_py_prep for MSLD (JV,LC)\n')
    fp.write('* (core from %s)\n* \n' % (reflig))
    fp.write('  %d %d\n' % (rtfvers1,rtfvers2))

    lp=open(outdir+'/lpsites.inp','w')
    lp.write("* Load LP Site Definitions (if applicable)\n*\n\n")

    # add in mass statements if there are any and remove redundancies
    # (assumes consistant atom typing across the different molecules)

    # core mass statements first
    addlist=[]
    tmplist=[x[0] for x in rtfinfo[refnum]['MASS']]  # generate tmp list of all atom types with mass statements
    for at in coreheader:
        for type in range(len(tmplist)):
            if rtfinfo[refnum]['ATTYPE'][at] == tmplist[type]:
                if not (rtfinfo[refnum]['ATTYPE'][at] in addlist):
                    addlist.append(rtfinfo[refnum]['ATTYPE'][at])
                    fp.write("MASS -1 %s %s %s\n" % (rtfinfo[refnum]['MASS'][type][0],
                             rtfinfo[refnum]['MASS'][type][1],rtfinfo[refnum]['MASS'][type][2]))
    # frag mass statements second
    for site in range(nsites):
        for frag in range(len(frags[site])):
            # figure mol equiv for frag
            for mol in range(len(mols)):
                if mols[mol] == frags[site][frag]:
                    break
            tmplist=[x[0] for x in rtfinfo[mol]['MASS']]
            # loop over atoms in each frag
            for at in Fatoms[site][frag]:
                for type in range(len(tmplist)):
                    if rtfinfo[mol]['ATTYPE'][at] == tmplist[type]:
                        if not (rtfinfo[mol]['ATTYPE'][at] in addlist):
                            addlist.append(rtfinfo[mol]['ATTYPE'][at])
                            fp.write("MASS -1 %s %s %s\n" % (rtfinfo[mol]['MASS'][type][0],
                                     rtfinfo[mol]['MASS'][type][1],rtfinfo[mol]['MASS'][type][2]))
    # continue with core.rtf atom block
    fp.write("\n")
    fp.write("RESI  %s    %5.3f\n" % (segid,qchk*-1))  # Writes the core net charge only 
    fp.write("GROUP \n")  # Not subdivided (an exercise left to the user if needed)
    for at in coreheader:
        if at[0] == 'H':
            pass
        elif at[0] == 'L': # if at[0:2] == 'LP':
            # LP atoms added, but no Hcore checks
            #fp.write("ATOM %-4s %-6s %9.5f \n" % (Ctrans[at],rtfinfo[refnum]['ATTYPE'][at],QQ.loc[at]))
            fp.write("ATOM %-4s %-6s %10.6f \n" % (Ctrans[at],rtfinfo[refnum]['ATTYPE'][at],QQ.loc[at]))
            # Figure out the "lonepair coli" data and print it to lpsites.inp
            for ln in rtfinfo[refnum]['LP']:
                if ln[0] == at:
                    lp.write("LONEPAIR COLI sele atom @ligseg @resnum %s end -\n" % (Ctrans[ln[0]]))
                    lp.write("              sele atom @ligseg @resnum %s end -\n" % (Ctrans[ln[1]]))
                    lp.write("              sele atom @ligseg @resnum %s end -\n" % (Ctrans[ln[2]]))
                    #lp.write("         DIST %s SCAL %s\ncoor shake\n\n" % (ln[4],ln[6])) # doesn't work with CGenFF
                    lp.write("         DIST %s SCAL %s\ncoor shake\n\n" % (ln[4],'0.00'))
        else:
            #fp.write("ATOM %-4s %-6s %9.5f \n" % (Ctrans[at],rtfinfo[refnum]['ATTYPE'][at],QQ.loc[at]))
            fp.write("ATOM %-4s %-6s %10.6f \n" % (Ctrans[at],rtfinfo[refnum]['ATTYPE'][at],QQ.loc[at]))
            # check for bonded H's
            if len(Hcore[at]) > 0:
                for at2 in Hcore[at]:
                    #fp.write("ATOM %-4s %-6s %9.5f \n" % (Ctrans[at2],rtfinfo[refnum]['ATTYPE'][at2],QQ.loc[at2]))
                    fp.write("ATOM %-4s %-6s %10.6f \n" % (Ctrans[at2],rtfinfo[refnum]['ATTYPE'][at2],QQ.loc[at2]))
    # core.rtf bond & impr lines
    for bd in rtfinfo[refnum]['BOND']:
        if (bd[0] in coreheader) and (bd[1] in coreheader):
            fp.write("BOND %-4s %-4s\n" % (Ctrans[bd[0]],Ctrans[bd[1]]))
    for impr in rtfinfo[refnum]['IMPR']:
        if (impr[0] in coreheader) and (impr[1] in coreheader) and \
           (impr[2] in coreheader) and (impr[3] in coreheader):
            fp.write("IMPR %-4s %-4s %-4s %-4s\n" % (Ctrans[impr[0]],Ctrans[impr[1]],\
                     Ctrans[impr[2]],Ctrans[impr[3]]))
    fp.write('PATCH FIRST NONE LAST NONE\n\nEND')
    fp.close()

    # write site[site]_sub[frag].rtf files
    for site in range(nsites):
        for frag in range(len(frags[site])):
            # figure out mol #
            for mol in range(len(mols)):
                if mols[mol] == frags[site][frag]:
                    break
            # write atom block
            fp=open(outdir+'/site'+str(site+1)+'_sub'+str(frag+1)+'_pres.rtf','w')
            fp.write('* fragment patch rtf file generated with py_prep for MSLD (JV,LC)\n')
            fp.write('* (fragment from %s)\n* \n' % (frags[site][frag]))
            fp.write('  %d %d\n\n' % (rtfvers1,rtfvers2))
            if ChkQChange:             # Calc frag net charge
                qsum=0.0
                for at in Fatoms[site][frag]:
                    qsum+=Qfrag[site][frag][at]
                qsum=float(str(np.around(qsum,decimals=dec)))
                fp.write('PRES p%d_%d    %5.3f\n' % (site+1,frag+1,qsum))
            else:
                fp.write('PRES p%d_%d    %5.3f\n' % (site+1,frag+1,siteavg[site]))
            fp.write('GROUP \n')
            for at in Fatoms[site][frag]:
                if at[0] == 'H':
                    pass
                elif at[0] == 'L': # if at[0:2] == 'LP':
                    # LP atoms added, but no Hfrag checks
                    #fp.write("ATOM %-4s %-6s %9.5f \n" % (Ftrans[site][frag][at],rtfinfo[mol]['ATTYPE'][at],Qfrag[site][frag][at]))
                    fp.write("ATOM %-4s %-6s %10.6f \n" % (Ftrans[site][frag][at],rtfinfo[mol]['ATTYPE'][at],Qfrag[site][frag][at]))
                    # Figure out the "lonepair coli" data and print it to lpsites.inp
                    for ln in rtfinfo[mol]['LP']:
                        if ln[0] == at:
                            lp.write("LONEPAIR COLI sele atom @ligseg @resnum %s end -\n" % (Ftrans[site][frag][ln[0]]))
                            lp.write("              sele atom @ligseg @resnum %s end -\n" % (Ftrans[site][frag][ln[1]]))
                            lp.write("              sele atom @ligseg @resnum %s end -\n" % (Ftrans[site][frag][ln[2]]))
                            #lp.write("         DIST %s SCAL %s\ncoor shake\n\n" % (ln[4],ln[6])) # doesn't work with CGenFF
                            lp.write("         DIST %s SCAL %s\ncoor shake\n\n" % (ln[4],'0.00'))
                else:
                    #fp.write("ATOM %-4s %-6s %9.5f \n" % (Ftrans[site][frag][at],rtfinfo[mol]['ATTYPE'][at],Qfrag[site][frag][at]))
                    fp.write("ATOM %-4s %-6s %10.6f \n" % (Ftrans[site][frag][at],rtfinfo[mol]['ATTYPE'][at],Qfrag[site][frag][at]))
                    # check for bonded H's
                    if len(Hfrag[site][frag][at]) > 0:
                        for at2 in Hfrag[site][frag][at]:
                            #fp.write("ATOM %-4s %-6s %9.5f \n" % (Ftrans[site][frag][at2],rtfinfo[mol]['ATTYPE'][at2],Qfrag[site][frag][at2]))
                            fp.write("ATOM %-4s %-6s %10.6f \n" % (Ftrans[site][frag][at2],rtfinfo[mol]['ATTYPE'][at2],Qfrag[site][frag][at2]))
            # write bond and impr lines (does NOT account for neighboring sites directly bonded to this site!)
            for bd in rtfinfo[mol]['BOND']:
                coreseries=cores.loc[mols[mol]][:]
                corelist=list(coreseries.values)
                if ((bd[0] in corelist) and (bd[1] in Fatoms[site][frag])) or \
                   ((bd[1] in corelist) and (bd[0] in Fatoms[site][frag])) or \
                   ((bd[0] in Fatoms[site][frag]) and (bd[1] in Fatoms[site][frag])):
                    # translate atom names then print
                    if bd[0] in corelist:
                        bdtmp=coreseries[coreseries.isin([bd[0]])].index[0]
                        at1 = Ctrans[bdtmp]
                    else:
                        at1 = Ftrans[site][frag][bd[0]]
                    if bd[1] in corelist:
                        bdtmp=coreseries[coreseries.isin([bd[1]])].index[0]
                        at2 = Ctrans[bdtmp]
                    else:
                        at2 = Ftrans[site][frag][bd[1]]
                    fp.write("BOND %-4s %-4s\n" % (at1,at2))
            for impr in rtfinfo[mol]['IMPR']:
                chk=0 # boolean for if we write it or not
                for at in impr:
                    # if at least one atom is a part of the fragment - then we want to print it
                    if at in Fatoms[site][frag]:
                        chk=1
                if chk:
                    # core atoms for this fragment specifically
                    coreseries=cores.loc[mols[mol]][:]
                    corelist=list(coreseries.values)
                    # chk to see if the impr spans two sites
                    crossite=[]
                    for at in impr:
                        if at == 'SKIP': 
                            break
                        if (not(at in corelist)) and (not(at in Fatoms[site][frag])):
                            chk=0   # at is in another site
                            if not (at in Aatoms[mol]):
                                impr.append('SKIP') # can't do crosssite IMPRs if the atoms aren't Aatoms
                            else:
                                # find 2nd site #
                                for site2 in range(nsites):
                                    if at == Aatoms[mol][site2]:
                                        crossite.append(site2)
                    if len(crossite) > 0:
                        crossite=list(set(crossite))
                        if len(crossite) > 1:
                            print(" !! Can't print IMPR line between three sites! (skipped)")
                            break
                    # write the impr line
                    if chk:
                        fp.write("IMPR")
                        for at in impr:
                            if at in corelist:
                                attmp=coreseries[coreseries.isin([at])].index[0]
                                fp.write(" %-4s" % (Ctrans[attmp]))
                            else:
                                fp.write(" %-4s" % (Ftrans[site][frag][at]))
                        fp.write("\n")
                    else: # cross site impropers written as comments
                        if impr[-1] == 'SKIP':
                            print(' !! Cannot print Cross Site IMPR between sites involving non-anchor fragment atoms')
                        else:
                            print(' ** Cross Site IMPR Found ('+frags[site][frag]+') ** Manually UNcomment as needed')
                            for frag2 in range(len(frags[crossite[0]])):
                                fp.write("!IMPR")
                                for at in impr:
                                    if at in corelist:
                                        attmp=coreseries[coreseries.isin([at])].index[0]
                                        fp.write(" %-4s" % (Ctrans[attmp]))
                                    elif at in Ftrans[site][frag]:
                                        fp.write(" %-4s" % (Ftrans[site][frag][at]))
                                    else: # it is in other site's fragment
                                        # figure out mol2 #
                                        for mol2 in range(len(mols)):
                                            if mols[mol2] == frags[crossite[0]][frag2]:
                                                break
                                        fp.write(" %-4s" % (Ftrans[crossite[0]][frag2][Aatoms[mol2][crossite[0]]]))
                                fp.write("\n")

            fp.write('\nEND')
            fp.close()
    lp.close()
    print()



    #################################################################
    ## Write CHARMM (prep) input script - full_ligand.prm not needed to write this file

    fp=open(outdir+'/nsubs','w')
    nblocks=0
    for site in range(nsites):
        fp.write("%d " % (len(frags[site])))
        nblocks+=len(frags[site])
    fp.close()
    fp=open(outdir+'/nblocks','w')
    fp.write("%d" % (nblocks))
    fp.close()
    fp=open(outdir+'/nreps','w')
    fp.write("1")
    fp.close()


    return

