# msld-py-prep

*msld-py-prep* scripts for creating a multiple topology model for Multisite λ-Dynamics (MSλD) in CHARMM. 

These scripts identify common atoms (similar partial atomic charge and identical atom types) across different compounds of interest with a maximum common substructure search. A charge renormalization algorithm is then implemented to generate a set of partial atomic charges suitable for multisite sampling of many chemical functional groups with MSλD. 

The output is a directory of CHARMM compatible files to use for MSλD simulations, including an example CHARMM input script. 


# Installation 
You can start using the scripts by cloning this repository and following the tutorial in the `examples` directory:

`git clone https://github.com/Vilseck-Lab/msld-py-prep.git`<br>
`cd msld-py-prep`

# Dependencies
- msld_py_prep requires Python 3 (or later versions), numpy, and pandas (or alternatively Anaconda) <br>
- Lg_Solvate.sh uses convpdb.pl from the [MMTSB toolset](https://github.com/mmtsb/toolset) <br>
- vis_check.py is written for use with [PyMOL](https://github.com/schrodinger/pymol-open-source) <br>

# Usage
Please see the examples directory for a detailed tutorial.

# Citation
Please cite the following reference:
"Optimizing Multisite λ-Dynamics Throughput with Charge Renormalization"
Jonah Z. Vilseck, Luis F. Cervantes, Ryan L. Hayes, and Charles L. Brooks
*Journal of Chemical Information and Modeling* **2022** *62 (6)*, 1479-1488;
DOI: 10.1021/acs.jcim.2c00047 

# Disclaimer
These scripts are provided as is and are subject to future modification.
