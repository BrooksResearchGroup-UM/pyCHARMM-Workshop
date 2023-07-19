# pycharmm flexible cdocker test case 

## Import module
import numpy as np
import pandas as pd
import pycharmm
import pycharmm.lib as lib
import pycharmm.read as read
import pycharmm.lingo as lingo
import pycharmm.settings as settings

from pycharmm.cdocker import Flexible_CDOCKER 


################################################################
# #
# #		Begin of pyCHARMM Flexible CDOCKER
# #
# ###############################################################

## Topology and parameter files 
settings.set_bomb_level(-1)
read.rtf('"../Toppar/top_all36_prot.rtf"')
read.rtf('"../Toppar/top_all36_cgenff.rtf"', append = True)
read.prm('"../Toppar/par_all36m_prot.prm"', flex = True)
read.prm('"../Toppar/par_all36_cgenff.prm"', append = True, flex = True)
settings.set_bomb_level(0)
read.stream('ligandrtf')

## Grid box information 
xcen = np.loadtxt('xcen', dtype = float)
ycen = np.loadtxt('ycen', dtype = float)
zcen = np.loadtxt('zcen', dtype = float)
maxlen = np.loadtxt('maxlen', dtype = float)

## Read in the receptor flexible side chain selection
flexchain = pd.read_csv('flexchain.csv', sep = '\t', index_col = 0)

## Flexible CDOCKER standard docking protocol 
sortedResult, dockResult = Flexible_CDOCKER(xcen = xcen, ycen = ycen, zcen = zcen, 
                                           maxlen = maxlen, flexchain = flexchain) 
print(sortedResult)
print(dockResult)

