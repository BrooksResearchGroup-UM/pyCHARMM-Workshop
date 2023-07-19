# Preparing Proteins for Simulation with crimm
## `crimm` is a package under development in the Brooks research group by Ziqiao (Truman) Xu. It is available for download on [GitHub](https://github.com/BrooksResearchGroup-UM/crimm). It is quite useful in preparing files from the RCSB for simulations with [CHARMM, pyCHARMM](https://academiccharmm.org) or [OpenMM](https://openmm.org). It inegrates a number of existing tools, including [Biopython](https://biopython.org), [propka](https://github.com/jensengroup/propka), [nglview](https://nglviewer.org), and others.
## crimm enables a number of tasks to be integrated into simulation workflows that include:
<h3>
<ol start="1">
    <li><p>directly parsing mmCIF files from the RCSB</p><br>
    <li><p>completing/filling in missing loops and gaps in protein strucutres by:</p>
    <ul>
    <li><p>scanning existing pdb strucures for homologues sequences with matching intact regions needed for completion</p></li>
    <li><p>importing the complete model from AlphaFold2 and utilizing the needed gap and loops</p></li>
    </ul><br>      
    <li><p>visualization of structural models within Jupyter Notebook/Jupyter Lab via [NGLView](http://nglviewer.org/nglview/latest/)</p><br>
    <li><p>calculate the pKa values of residues and prepare patches based on the CHARMM topology files for those residues to be titrated given the input pH and propKa computed pKa</p><br>
    <li><p>create CHARMM patches for disulfide bonds and patch them into the structure</p><br>
    <li><p>evolving adaptors that fully integrate pyCHARMM and OpenMM (in progress)</p><br>
    </li>
    </ol>
</h3>

## To explore this Jupyter Notebook Tutorial on preparing proteins for simiulation you will need to install the `crimm` package and other needed python modules in the environment you created for this project, as detailed in the <i>0Install</i> folder.

### To added the needed package (propka):
```csh
conda activate <name_of_environment>
mamba install -c conda-forge propka
pip install crimm
```

### Now you are all set to run the tutorial and crimm is part of your working packages in this conda environment.
