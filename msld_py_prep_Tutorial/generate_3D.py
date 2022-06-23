#! /usr/bin/env python

from numpy import genfromtxt
import collections
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit import ML
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import TemplateAlign
#from rdkit.Chem import MCS
from rdkit.Chem import rdFMCS
from rdkit.ML.Cluster import ClusterVis
from rdkit.ML.Cluster import ClusterUtils
from copy import deepcopy
import numpy as np
import pandas as pd
import sys,os
import pickle

####
#### Run with rdenv environment activated!
#### (source activate rdenv)
####

#############################################
# (1) Load in mol2 files 
#############################################

def load_sdf(fn):
    """
    Load sdf specified file fn into rdkit and returns list of 
    RDKit mol objects for each molecule in sdf. OpenBabel
    kekulizes molecules in the sdf file when converting from a
    mol2 file. This facilitates input into RDKit.
    """
    suppl = Chem.SDMolSupplier(fn,removeHs=False)
    if None in suppl:
        print("COULD NOT LOAD MOLECULE IN %s" % fn)
        quit()
    return [x for x in suppl]

def load_mol2(fn):
    """
    Load mol2 specified file fn into rdkit and returns mol object.
    If molecule is kekulized in mol2, this is the better way of
    loading into RDKit.
    """
    mol = Chem.MolFromMol2File(fn,removeHs=False)
    if mol == None:
        print("Could not load %s into RDKit successfully." % fn)
        quit()
    return mol

def load_smiles(mols):
    """
    Load smiles specified in list into RDKit and return list of
    RDKit objects
    """ 
    molss = []
    for idx, mol in enumerate(mols):
        # print(idx,mol)
        molobj = Chem.MolFromSmiles(mol)
        if not molobj:
            print("Could not load SMILES INTO RDKit.\nCheck\
                    input SMILES: "+str(mol))
            quit()
        molss.append(molobj)
    return molss

class RDKit_Tools:
    def __init__(self,f,smiles_col, names_col, delim=",", rng="all"):
        """
        Input:
        - `f`        : csv file name
        - `coln`     : SMILES column name in `f`
        - `rng`      : row range for which molecules to draw
        """
        self.rng = rng
        self.df = pd.read_csv(f,sep=delim)
        self.smiles_col_name = smiles_col
        self.names_col_name = names_col 
        self.smiles_list = list(self.df[smiles_col])
        self.mol_obj_list = load_smiles(self.smiles_list)
        self.names_col = self.df[names_col]
        if rng != "all":
            self.smiles_list = self.smiles_list[rng]
            self.mol_obj_list = self.mol_obj_list[rng]

    def get_dataframe(self,save=False, filter_idx=[],delim=","):
        """
        Returns dataframe from Class Object.
        In:
        - `save`       : Name for dataframe csv file to
          (str)          be generated. If None, then just
                         return df

        - `filter_idx` : Filter Class object df by
          (int list)     row index list
        """
        if filter_idx:
            # Make sure all indices in list are integers
            filter_idx = list(map(int,filter_idx))
            dfout =  self.df.iloc[filter_idx]
            dfout =  dfout[[self.names_col_name,self.smiles_col_name]] 
        else: 
            dfout = self.df[[self.names_col_name,self.smiles_col_name]]

        if save:
            dfout.to_csv(save,sep=delim,index=False)
        
        return dfout

    def draw_smiles(self,img_fname, mol_list="all"):
        if mol_list == "all":
            mol_list = self.mol_obj_list

        if not mol_list:
            print("Could not plot empty molecule object list")
            quit()

        img=Draw.MolsToGridImage(mol_list,molsPerRow=4,subImgSize=(200,200),returnPNG=False)
        return img.save(img_fname)

    def filter_molecules(self, pat, hyd=True, mol_list="all"):
        """
        Filter molecules based on pattern.
        In:
        - `pat` : smiles pattern to match

        Out:
        - `matches` : list of RDK mol objects containing `pat`
        - `indices` : list of mol indices containing `pat`
        """
        print("Matching molecules to pattern...")
        if mol_list == "all":
            mol_list = self.mol_obj_list
            if hyd:
                mol_list = [Chem.AddHs(x) for x in mol_list]

        params = Chem.SmilesParserParams()
        params.removeHs= not hyd # draw and work with explicit Hs
        pat = Chem.MolFromSmiles(pat, params)
        if hyd:
            pat = Chem.AddHs(pat,explicitOnly=hyd)

        matches = []
        indices = []
        for i, mol in enumerate(mol_list):
            if mol.HasSubstructMatch(pat):
                matches.append(mol)
                indices.append(i)

        matches = [Chem.RemoveHs(match) for match in matches]        
        print("Done.\n")
        return matches, indices 


    def write_mol_block(self, molblock, fname):
        """
        Write molblock to fname with extension .mol
        """
        with open(fname+".mol", "w") as f:
            f.write(molblock) 

    def find_bond_groups(self,mol):
        """
        Helper function for get_num_rotatable()
        Taken from https://www.rdkit.org/docs/Cookbook.html
        'Find groups of contiguous rotatable bonds and return them sorted by decreasing size'
        """
        rot_atom_pairs = mol.GetSubstructMatches(RotatableBondSmarts)
        rot_bond_set = set([mol.GetBondBetweenAtoms(*ap).GetIdx() for ap in rot_atom_pairs])
        rot_bond_groups = []
        while (rot_bond_set):
            i = rot_bond_set.pop()
            connected_bond_set = set([i])
            stack = [i]
            while (stack):
                i = stack.pop()
                b = mol.GetBondWithIdx(i)
                bonds = []
                for a in (b.GetBeginAtom(), b.GetEndAtom()):
                    bonds.extend([b.GetIdx() for b in a.GetBonds() if (
                        (b.GetIdx() in rot_bond_set) and (not (b.GetIdx() in connected_bond_set)))])
                connected_bond_set.update(bonds)
                stack.extend(bonds)
            rot_bond_set.difference_update(connected_bond_set)
            rot_bond_groups.append(tuple(connected_bond_set))
        return tuple(sorted(rot_bond_groups, reverse = True, key = lambda x: len(x)))

    def get_num_rotatable(self,mol_list="all", mol_names="all", save=True, typerot="RDK",strict=True):
        """
        Get number of rotatable bonds for a list of 

        In:
        -`typerot` : definition of rotatable bond. Can be "Lip" (Lipinski), "RDK" (RDKit), or
                     "contig" (Contiguous. Yields largest number of contiguous rotatable bonds
                     for every molecule in mol_list.
        """
        # Make sure molecule object and name lists are same length  
        if len(mol_names) != len(mol_list):
            print("List of names and list of molecules are not the same length.")
            quit()
        
        # Get molecule object and name lists 
        if mol_list == "all":
            mol_list = self.mol_obj_list
            mol_names = self.names_col

        
        # if typerot == "Lip":
        #     rotable = [Chem.Lipinski.NumRotatableBonds(x) for x in mol_list]

        if typerot == "RDK":
            rotable = [Chem.rdMolDescriptors.CalcNumRotatableBonds(x,strict=strict) for x in mol_list]
        
        elif typerot == "contig":
            # Find groups of contiguous rotatable bonds in mol
            bond_groups = [self.find_bond_groups(mol) for mol in mol_list]
            # As bond groups are sorted by decreasing size, the size of the first group (if any)
            # is the largest number of contiguous rotatable bonds in mol
            rotable = [len(bond_group[0]) if bond_group else 0 for bond_group in bond_groups]

        else:
            print("Could not identify type of rotatable bond definition.")
            quit()

        return rotable



    def generate_3D(self, mol_list="all", mol_names="all", mini=True, ref=None):
        """
        Generate 3D coordinates in mol format.
        """
        print("Generating 3D coordinates...")
        generated_mols = []

        # Make sure molecule object and name lists are same length  
        if len(mol_names) != len(mol_list):
            print("List of names and list of molecules are not the same length.")
            quit()

        # Get Reference 2D and 3D Structures (if specified)
        if ref:
            if isinstance(ref,int):
                generated_mols.append(mol_names[ref])
                ref=mol_list[ref]
            
            ref_2D = ref
            ref_3D = Chem.AddHs(ref_2D)
            AllChem.EmbedMolecule(ref_3D)
            if mini:
                AllChem.MMFFOptimizeMolecule(ref_3D)

        # Get molecule object and name lists 
        if mol_list == "all":
            mol_list = self.mol_obj_list
            mol_names = self.names_col

        # Check if molecule name list has correct element type.
        # Coerce to strings if not.
        if not any([isinstance(name, str) for name in mol_names]):
            mol_names = list(map(str,mol_names))

        # Generate 3D Coords, minimize (if True), align to ref (if given)
        for i,(name,mol) in enumerate(zip(mol_names, mol_list)):
            m2 = Chem.AddHs(mol)
            AllChem.EmbedMolecule(m2)
            if mini:
                AllChem.MMFFOptimizeMolecule(m2)
            if ref:
                if mol == ref_2D:
                    self.write_mol_block(Chem.MolToMolBlock(m2),name)
  
                else:
                    at_matches = mol.GetSubstructMatch(ref_2D)
                    if not at_matches:
                        try:
                            Chem.rdMolAlign.AlignMol(m2,ref_3D)
                            self.write_mol_block(Chem.MolToMolBlock(m2),name)
                            generated_mols.append(name)
                        except:
                            print("Molecule %s in list does not contain matching atoms to ref." % name)

                    else:
                        atmap = [(j,k) for j,k in zip(at_matches, list(range(len(at_matches))))]
                        # print(m2)
                        # print(ref_3D)
                        # self.draw_smiles("testtt.png",[m2,ref_3D])
                        # quit()
                        Chem.rdMolAlign.AlignMol(m2,ref_3D,atomMap=atmap)
                        self.write_mol_block(Chem.MolToMolBlock(m2),name)
                        generated_mols.append(name)

            else:
                self.write_mol_block(Chem.MolToMolBlock(m2),name)
                generated_mols.append(name)


        return generated_mols


    def get_group_bins(self, fnames):
        """ 
        `fnames` : list of file names
        """
        # Coerce fnames to str elements
        fnames = list(map(str, fnames)) 
        lists = []
        for i in fnames:
            lists.append(genfromtxt(i+".chkmol", delimiter=';',dtype=str))

        lists = list(map(list, lists))
        lists = [list(filter(lambda x: x != '', el)) for el in lists]

        unique_groups = []
        bins = []
        for mol_groups in lists:
            for group in mol_groups:
                if group not in unique_groups:
                    unique_groups.append(group)
                    bins.append(0)
                bins[unique_groups.index(group)] += 1
        return bins, unique_groups


    def histogram(self, fout, bins, labels):
        """
        Plots histogram based on bins and labels for those bins
        """
        plt.bar(labels, bins)
        plt.suptitle("Chemical Group Distribution of Filtered Molecules")
        plt.xticks(rotation='82.5')
        plt.subplots_adjust(bottom=0.35)
        plt.savefig(fout, dpi=400)
        plt.close()     


# dataset = "CHEMBL_Dataset.csv"
# dfout= "filtered_df.csv"
# dfout2= "filtered_df1.csv"
# histname = "chemgroupfreq.pdf"
# pattern2 = "C1(=O)NC(=O)CC1" 
# pattern1 = "C1C=CC=C2C(=O)C=COC2=1" 
# names_col = "ChEMBL ID"
# smiles_col = "Smiles"
# delim=";"
# 
# 
# mols = RDKit_Tools(dataset, smiles_col=smiles_col, names_col=names_col,delim=delim) # Load molecules into RDKit
# matches,indices = mols.filter_molecules(pattern1) # Filter based on smiles pattern
# pattern1 = load_smiles([pattern1])[0]
# mols.draw_smiles("filtered_mols.png",matches[:50]) # Visualize filtered molecules
# filtered_names = mols.generate_3D(matches, indices, ref=pattern1) # Generate mol files after minimization
# print(filtered_names)
# mols.get_dataframe(save=dfout, filter_idx=filtered_names,delim=delim)


# mols1 = RDKit_Tools(dfout, smiles_col=smiles_col, names_col=names_col,delim=delim) # Load molecules into RDKit
# matches,indices = mols1.filter_molecules(pattern2) # Filter based on smiles pattern
# mols1.draw_smiles("filtered_mols1.png",matches[:50]) # Visualize filtered molecules
# pattern2 = load_smiles([pattern2])[0]
# filtered_names = mols1.generate_3D(matches, indices, ref=pattern2) # Generate mol files after minimization
# mols1.get_dataframe(save=dfout2, filter_idx=filtered_names,delim=delim)
# print(filtered_names)

# print(mols.get_num_rotatable(matches, indices))
