#! /bin/bash

##
## Solvate the largest ligand in a cubic box of water with convpdb.pl (from MMTSB toolset)
## the optional cutoff_distance defines the minimum distance from molecule to edge of solvation box
##
## Usage:
##	./Lg_Solvate.sh [cutoff_distance]


cutoff=$1
if [[ $cutoff == '' ]] ; then
   cutoff=12
fi
echo "Using a solvate cutoff of $cutoff"

name=`more build*/name`

InDir=./build.${name}
OutDir=solv_prep
Dir=`pwd`
LgLig=large_lig.pdb

cp -r $InDir ./$OutDir

convpdb.pl -solvate -cubic -cutoff $cutoff $InDir/$LgLig | grep TIP3 | convpdb.pl -segnames > $OutDir/solvent.pdb

convpdb.pl -solvate -cubic -cutoff $cutoff $InDir/$LgLig | grep LIG | convpdb.pl -segnames | sed "s/PRO0/LIG/g" > $OutDir/LgLig.pdb

XYZ=`awk 'BEGIN {CC=0} {if ($1=="ATOM" && $2==1) {if (CC==0) {CC=1; X=$6; Y=$7; Z=$8} else {print X-$6,Y-$7,Z-$8}}}' $OutDir/LgLig.pdb $InDir/$LgLig`

echo $XYZ


cd $OutDir
FILES=`ls core.pdb site*.pdb`
cd $Dir

for FILE in $FILES
do

  convpdb.pl -translate $XYZ -segnames $InDir/$FILE | sed "s/PRO0/LIG/g" > $OutDir/$FILE

done
