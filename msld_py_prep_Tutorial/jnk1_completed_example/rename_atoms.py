#! /usr/bin/env python

#
# Take a mol2 file and rename all the atom names to be unique!
#

import sys,os

if len(sys.argv) != 2:
    print("\n This script needs a command line argument to work:")
    print("    1st = name of the mol2 file to process")
    quit()

filename = sys.argv[1]

fp=open(filename,'r')
gp=open('tmp.mol2','w')
line=fp.readline()
while line:
    if line[0:13] == '@<TRIPOS>ATOM':
        gp.write(line)
        line=fp.readline()
        while line[0:13] != '@<TRIPOS>BOND':
            tmp=line.split()
            #print(len(tmp),tmp)
            if int(tmp[0]) < 10:
                tmp[1] = tmp[1]+'00'+tmp[0]
            elif int(tmp[0]) > 9 and int(tmp[0]) < 100:
                tmp[1] = tmp[1]+'0'+tmp[0]
            elif int(tmp[0]) > 99 and int(tmp[0]) < 1000:
                tmp[1] = tmp[1]+tmp[0]
            else:
                raise ValueError("More than 1000 atoms not handled")
            gp.write("%7s %-6s  %10s%10s%10s %-6s%3s  %-4s     %9s\n" % (tmp[0],tmp[1],
                     tmp[2],tmp[3],tmp[4],tmp[5],tmp[6],tmp[7],tmp[8]))
            line=fp.readline()

        #gp.write(line)
    else:
        gp.write(line)
        line=fp.readline()

fp.close()
gp.close()

# overwrite the original file
os.system('mv tmp.mol2 '+filename)

# finished
