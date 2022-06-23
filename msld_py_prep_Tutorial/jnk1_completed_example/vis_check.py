##
## visual check in PyMOL
##
## Usage:
## 	open PyMOL, then run this script in the command line ("run vis_check.py")
##

mols=[]   # list of file names
cores=[]  # list of core atoms in those files
mcs_filename = 'MCS_for_MSLD.txt'  # check that this matches your filename

# read in info
fp=open(mcs_filename,'r')
line=fp.readline()
while line:
    if line[0:4] == 'CORE':
        line=fp.readline()
        while (line != '\n') and (line[0:4] != 'ANCH'): 
            lns=line.split()
            mols.append(lns[0])
            cores.append(lns[1:])
            line=fp.readline()
    line=fp.readline()
    if line[0:4] == 'ANCH':
        break

# make pretty pictures
cmd.set('sphere_scale',value=0.25)

for mol in range(len(mols)):
    # load in the mol2 file as an independent object
    cmd.load(mols[mol]+'.mol2',mols[mol])
    cmd.hide("everything",mols[mol])
    cmd.show("lines",mols[mol])

    # read *.core file and highlight the "core" atoms 
    sele_ats=cores[mol]
    for at in sele_ats:
        cmd.show(representation="spheres",selection="name "+at+' and '+mols[mol])

cmd.label(expression="name")
cmd.set('label_size',value=20)
cmd.bg_color(color="grey50")

# finished
