#!/bin/bash
# ##############################################################
# # This is used to test flexible cdocker
# # Use T4 L99A 181L as the test receptor
# ##############################################################

mkdir -p fcdocker
for pdb in `cat ligands/list`; do
	if ! [ -s fcdocker/${pdb}/clusterResult.csv ]; then
		rm -rf fcdocker/${pdb}
		cp ligands/${pdb}_ligand.pdb ./ligand.pdb
		cp ligands/${pdb}_ligandrtf ./ligandrtf
		python standard.py > out.file
		mv out.file dockresult
		mv dockresult fcdocker/${pdb}
		rm ligand.pdb ligandrtf
		fi
	done
