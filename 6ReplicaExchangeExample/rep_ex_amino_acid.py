'''
Hamiltonian replica exchange simulation of 2-state amino acid in vacuum.
Only perform exchange in lambda_chem space.
Run NVT simulations using CPU (w/o PME, since "BLOCK not implemented for COLFFT" in CPU).
'''

from mpi4py import MPI
import FastMBAR

import os
import sys
import subprocess
import numpy as np
import pandas as pd

import pycharmm
import pycharmm.read as read
import pycharmm.lingo as lingo
import pycharmm.generate as gen
import pycharmm.settings as settings
import pycharmm.write as write
import pycharmm.nbonds as nbonds
import pycharmm.ic as ic
import pycharmm.coor as coor
import pycharmm.energy as energy
import pycharmm.dynamics as dyn
import pycharmm.minimize as minimize
import pycharmm.crystal as crystal
import pycharmm.select as select
import pycharmm.image as image
import pycharmm.psf as psf
import pycharmm.param as param
import pycharmm.cons_harm as cons_harm
import pycharmm.cons_fix as cons_fix
import pycharmm.shake as shake
import pycharmm.scalar as scalar
import pycharmm.charmm_file as charmm_file

# user defined variables
nl_chem=int(sys.argv[1]) # number of windows in lambda_chem
tsim=298.15              # simulation temperature
ncycles=1000             # number of cycles
nsteps_per_cyc=1000      # number of steps per cycle 

# residue name and segment id specific to this system
comp='lys'
segid='PROT'

# if we find that additional simulation time is needed after the first run, 
# we can restart/continue our simulation by chaning variable itt to 2
itt=1                    # run index
ittm1=itt-1              # index of previous run

# for free energy calculation, we need to recompute energies for dcd trajectories.
# if we saved dcd too frequently and would like to reduce the frequency in energy
# calculations, we could increase variable skip_dcd.
skip_dcd=1

# boltzmann constant and inverse temperature 
kB=0.0019872042586408316 # kcal/mol/K
beta=1/kB/tsim

###################################
# mpi
comm=MPI.COMM_WORLD
nproc=comm.Get_size()
rank=comm.Get_rank()
repl_id=rank

###################################
# directories
param_dir='toppar'
wrk_dir='.'
out_dir='aa'+str(rank)

###################################
# input and output files
inp_psf_fn=wrk_dir+'/'+comp+'.psf'
inp_pdb_fn=wrk_dir+'/'+comp+'.pdb'
inp_box_fn=wrk_dir+'/box_size.dat'
dcd_fn=out_dir+'/dcd/prod'+str(itt)+'.dcd'
rst_fn=out_dir+'/res/prod'+str(itt)+'.res'
rpr_fn=out_dir+'/res/prod'+str(ittm1)+'.res'
pdb_fn=out_dir+'/pdb/prod'+str(itt)+'.pdb'
log_fn=out_dir+'/out/prod'+str(itt)+'.out'
exc_fn=out_dir+'/exc/exch'+str(itt)+'.dat'
his_fn=out_dir+'/his/cond_his'+str(itt)+'.dat'
cnd_fn=out_dir+'/../cond.dat'
exch_all_fn=out_dir+'/../exchange-all-'+str(itt)+'.dat'
exch_all_pre_fn=out_dir+'/../exchange-all-'+str(ittm1)+'.dat'
ratio_fn=out_dir+'/../ratio-'+str(itt)+'.dat'
fe_fn=out_dir+'/../mbar.out'


# lambda_chem
l_g='{0:.4f}'.format(1/(nl_chem-1)*repl_id)

# file units
log_unit=30
dcd_unit=40
rst_unit=50
rpr_unit=70  # unit previous restart file

# create directories
os.system('mkdir -p '+out_dir+'/dcd')
os.system('mkdir -p '+out_dir+'/res')
os.system('mkdir -p '+out_dir+'/pdb')
os.system('mkdir -p '+out_dir+'/out')
os.system('mkdir -p '+out_dir+'/exc')
os.system('mkdir -p '+out_dir+'/his')

def open_clog():
    '''
    specify charmm output file for a given system
    '''
    clog=charmm_file.CharmmFile(file_name=log_fn,file_unit=log_unit,read_only=False,formatted=True)
    lingo.charmm_script('outu '+str(log_unit))
    return clog

def close_clog(clog):
    lingo.charmm_script('outu 6')
    clog.close()

def read_param():
    read.rtf(param_dir+'/top_all36_prot_hedi_xrliu.rtf')
    read.prm(param_dir+'/par_all36m_prot.prm',flex=True)

def read_init():
    read.psf_card(inp_psf_fn)
    read.pdb(inp_pdb_fn,resid=True)

def titr_grp(resn):
    '''
    atom names in the titratable group for a given amino acid
    '''
    if resn == 'ASP': 
       type_list=['CB','HB1','HB2','CG','OD1','OD2','HD1','HD2']
    elif resn == 'GLU':
       type_list=['CG','HG1','HG2','CD','OE1','OE2','HE1','HE2']
    elif resn == 'HSP' or resn == 'HSD' or resn == 'HSE':
       type_list=['CB','HB1','HB2','CG','ND1','HD1','CE1','HE1','CD2','HD2','NE2','HE2']
    elif resn == 'LYS':
       type_list=['CE','HE1','HE2','NZ','HZ1','HZ2','HZ3']
    return type_list

def define_sub(aa,ires):
    '''
    define substituents
    '''
    resn=aa.upper()
    resid=str(ires)
    atom_list=titr_grp(resn)
    sele_w=~pycharmm.SelectAtoms(select_all=True)
    sele_m=~pycharmm.SelectAtoms(select_all=True)
    sele_p=~pycharmm.SelectAtoms(select_all=True)
    sele_u=~pycharmm.SelectAtoms(select_all=True)
    for name in atom_list:
        sele_w = sele_w | pycharmm.SelectAtoms().by_res_and_type(seg_id=segid,res_id=resid,atom_type=name+'W')
        sele_m = sele_m | pycharmm.SelectAtoms().by_res_and_type(seg_id=segid,res_id=resid,atom_type=name+'M')
        sele_p = sele_p | pycharmm.SelectAtoms().by_res_and_type(seg_id=segid,res_id=resid,atom_type=name+'P')
        sele_u = sele_u | pycharmm.SelectAtoms().by_res_and_type(seg_id=segid,res_id=resid,atom_type=name+'U')
    sele_w.store('site{}{}subW'.format(resn,resid))
    sele_m.store('site{}{}subM'.format(resn,resid))
    sele_p.store('site{}{}subP'.format(resn,resid))
    sele_u.store('site{}{}subU'.format(resn,resid))

class set_block:
    '''
    set up block for simulation
    '''
    def __init__(self,aa='',ires='1',nstates='2',lmd=''):
        self.resn=aa.upper()
        self.resid=str(ires)
        self.nstat=int(nstates)
        self.lmd=float(lmd)
    def run(self):
        self.set_name()
        self.reset()
        self.call()
    def set_name(self):
        sele_name_w='site%s%ssubW'%(self.resn,self.resid)
        sele_name_m='site%s%ssubM'%(self.resn,self.resid)
        sele_name_p='site%s%ssubP'%(self.resn,self.resid)
        sele_name_u='site%s%ssubU'%(self.resn,self.resid)
        if self.resn == 'ASP' or self.resn == 'GLU':
           self.sele_name=[sele_name_w,sele_name_m,sele_name_p]
        elif self.resn == 'HSD':
           self.sele_name=[sele_name_w,sele_name_u]
        elif self.resn == 'HSE':
           self.sele_name=[sele_name_w,sele_name_m]
        elif self.resn == 'LYS':
           self.sele_name=[sele_name_w,sele_name_m]
        elif self.resn == 'HSP':
           self.sele_name=[sele_name_w,sele_name_u,sele_name_m]
    def reset(self):
        lingo.charmm_script('''
              BLOCK %s
                    clear 
              END
        '''%(self.nstat+1))
    def call(self):
        excl_list=''
        for i in range(0,self.nstat):
            for j in range(i+1,self.nstat):
                excl_list=excl_list+str(i+2)+' '+str(j+2)+' '
        lingo.charmm_script('''
              BLOCK {nblk}
                    call 2 sele {sele1} end
                    call 3 sele {sele2} end

                    coef 1 1 1.0
                    coef 1 2 {lmbd1}
                    coef 1 3 {lmbd2}

                    coef 2 2 {lmbd1}
                    coef 2 3 0.0

                    coef 3 3 {lmbd2}

                    excl {excl}

                    rmla bond thet impr
                    msma
              END
        '''.format(nblk=self.nstat+1,
                   sele1=self.sele_name[0],
                   sele2=self.sele_name[1],
                   lmbd1=self.lmd,
                   lmbd2=1.0-self.lmd,
                   excl=excl_list))

nbond={'elec': True,
       'atom': True,
       'cdie': True,
       'fswitch': True,
       'eps': 1,
       'vdw': True,
       'vatom': True,
       'vfswitch': True,
       'cutnb': 14,
       'cutim': 14,
       'ctofnb': 12,
       'ctonnb': 10,
       'ewald': False
      }

def setup_nb(nb):
    nb_nopme=pycharmm.NonBondedScript(**nb)
    nb_nopme.run()
    energy.show()

def set_pbc():
    size=np.loadtxt(inp_box_fn,usecols=(0))
    crystal.define_cubic(size)
    crystal.build(12)
    # shift box center to origin
    pos=coor.get_positions()
    coor.set_positions(pos-size*0.5)
    offset=0
    image.setup_segment(offset,offset,offset,segid)
    # write coordinates
    write.coor_pdb(out_dir+'/init.pdb')

class resd_restrain:
    '''
    restrain all or heavy atoms between analogous groups of the 2 or 3 protonation states using resd restraint.
    dynamics doesn't work if analogous atoms are on top of each other.
    '''
    def __init__(self,aa='',ires='',heavy=True):
        self.heavy=heavy
        self.resn=aa.upper()
        self.resid=str(ires)
    def run(self):
        self.reset()
        self.restrain()
        self.show()
    def reset(self):
        lingo.charmm_script('RESD reset')
    def resd(self,name):
        sel_w=pycharmm.SelectAtoms().by_res_and_type(seg_id=segid,res_id=self.resid,atom_type=name+'W')
        sel_m=pycharmm.SelectAtoms().by_res_and_type(seg_id=segid,res_id=self.resid,atom_type=name+'M')
        sel_p=pycharmm.SelectAtoms().by_res_and_type(seg_id=segid,res_id=self.resid,atom_type=name+'P')
        sel_u=pycharmm.SelectAtoms().by_res_and_type(seg_id=segid,res_id=self.resid,atom_type=name+'U')
        w=sel_w.get_n_selected()
        m=sel_m.get_n_selected()
        p=sel_p.get_n_selected()
        u=sel_u.get_n_selected()
        if w==1 and m==1:
           lingo.charmm_script('resd kval 100 rval 0.0 eval 2 ival 1 1.0 %s %s %s %s %s %s'%
                               (segid,self.resid,name+'W',segid,self.resid,name+'M'))
        if w==1 and p==1:
           lingo.charmm_script('resd kval 100 rval 0.0 eval 2 ival 1 1.0 %s %s %s %s %s %s'%
                               (segid,self.resid,name+'W',segid,self.resid,name+'P'))
        if w==1 and u==1:
           lingo.charmm_script('resd kval 100 rval 0.0 eval 2 ival 1 1.0 %s %s %s %s %s %s'%
                               (segid,self.resid,name+'W',segid,self.resid,name+'U'))
    def restrain(self):
        atom_list=titr_grp(self.resn)
        for atom_name in atom_list:
            if self.heavy==True:  # restrain heavy atoms
               if atom_name[0]!='H':
                  self.resd(atom_name)
            else:                 # restrain all atoms
               self.resd(atom_name)
    def show(self):
        lingo.charmm_script('print resdistances')
        #energy.show()
        lingo.charmm_script('energy omm')

def mini():
    minimize.run_sd(nstep=50,nprint=10,tolenr=1e-3,tolgrd=1e-3)
    energy.show()
    # write coordinates
    write.coor_pdb(out_dir+'/mini.pdb')

def dyn_init():
    lingo.charmm_script('''
          faster on
    ''')
    shake.on(bonh=True, fast=True, tol=1e-7)
    n = psf.get_natom()
    scalar.set_fbetas([1.0] * n)
    

# CPT simulation using omm
dyn_dict={
    'leap': True,
    'verlet': False,
    'cpt': False,
    'new': False,
    'langevin': True,
    'omm': False,
    'start': True,
    'timestep': 0.002,
    'nstep': nsteps_per_cyc,
    'nsavc': 0,
    'nsavv': 0,
    'nsavl': 0,      # frequency for saving lambda values in lamda-dynamics
    'nprint': nsteps_per_cyc, # Frequency to write to output
    'iprfrq': nsteps_per_cyc, # Frequency to calculate averages
    'isvfrq': nsteps_per_cyc, # Frequency to save restart file
    'ntrfrq': 5000,
    'inbfrq':-1,
    'ihbfrq':0,
    'ilbfrq':0,
    'imgfrq':-1,
    'iunrea':-1,
    'iunwri':rst_unit,
    'iuncrd':-1,
    'iunldm':-1,
    'firstt': tsim,
    'finalt': tsim,
    'tstruct': tsim,
    'tbath': tsim,
    'iasors': 1, # assign velocities
    'iasvel': 1, # method for assignment of velocities during heating & equil when IASORS is nonzero.
                 # This option also controls the initial assignment of velocities 
    'iscale': 0, # not scale velocities on a restart
    'scale': 1,  # scaling factor for velocity scaling
    'ichecw': 0, # not check temperature
    'echeck': -1 # not check energy
}


def nvt(run):
    rst=charmm_file.CharmmFile(file_name=rst_fn,file_unit=rst_unit,read_only=False,formatted='formatted')

    if itt>=2:
       if run==0:
          dyn_dict['start']=False
          dyn_dict['restart']=True
          dyn_dict['iunrea']=rpr_unit
          rpr=charmm_file.CharmmFile(file_name=rpr_fn,file_unit=rpr_unit,read_only=True,formatted='formatted')
       else:
          dyn_dict['start']=True
          dyn_dict['restart']=False
          dyn_dict['iunrea']=-1

    npt_prod=pycharmm.DynamicsScript(**dyn_dict)
    npt_prod.run()
    rst.close()
    if itt>=2 and run==0:
       rpr.close()

    # last snapshot
    write.coor_pdb(pdb_fn)

def write_cond_file():
    sendbuf=np.array([float(l_g)],dtype='f')
    recvbuf=None
    if rank==0:
       recvbuf = np.empty([nproc, 1], dtype='f')
    comm.Gather(sendbuf,recvbuf,root=0)
    if rank==0:
       np.savetxt(cnd_fn,recvbuf,fmt='%.4f')
    comm.barrier()
    # broadcast condition array
    if rank==0:
        cond_array=recvbuf
    else:
        cond_array=np.empty([nproc, 1], dtype='f')
    comm.Bcast(cond_array, root=0)
    return cond_array

def get_neighbor(current_cond_id,run_id):
    if run_id%2==0:
       if current_cond_id%2==0:
          neighbor_cond_id=current_cond_id+1
       else:
          neighbor_cond_id=current_cond_id-1
    else:
       if current_cond_id%2==0:
          neighbor_cond_id=current_cond_id-1
       else:
          neighbor_cond_id=current_cond_id+1
    if neighbor_cond_id<0:
       neighbor_cond_id=0
    if neighbor_cond_id>=nproc:
       neighbor_cond_id=nproc-1
    return neighbor_cond_id

class rep_ex:
      '''
      perform Hamiltonian replica exchange (lambda_chem only)
      '''
      def __init__(self,ncycle='',condid='',lmd='',conditions=''):
          self.cond_id=int(condid)
          self.cond_array=conditions
          self.lmd=float(lmd)
          self.nrun=int(ncycle)
      def run(self):
          self.open_files()
          for i in range(0,self.nrun):
              #print("run {}, rank {} lambdas are intra {} and inter {}".format(i,rank,self.intra,self.inter))
              nvt(i)
              self.write_trj(i)
              self.cond_data=self.gather_condid(i)
              self.lmd_new=self.swap_neighbor(i)
              #print("after exchange, rank {} lambdas are l {} ".format(rank,self.lmd_new))
              self.get_ener_new()
              self.metropolis()
              #print('run ',i, 'cond ',self.cond_id,'Accept',self.Accept)
              self.write_exch_all(i)
              self.update_cond()
          self.close_files()
      def open_files(self):
          # all exchange information
          if rank==0:
             self.exch_all=open(exch_all_fn,'w',buffering=1)
             self.exch_all.write('# replica temp. ener. neighbor ntemp nene prob p success? newrep\n')
          # exchange information for each replica
          self.exch=open(exc_fn,'w',buffering=1)
          self.exch.write('#run repl_id ener_im neighbor_repl_id ener_in prob rand accept replica_next\n')
          # condition history for each replica
          self.hist=open(his_fn,'w',buffering=1)
          # dcd for each replica
          self.dcd=charmm_file.CharmmFile(file_name=dcd_fn,file_unit=dcd_unit,read_only=False,formatted=False)
          lingo.charmm_script('''
          traj IWRITE {} NWRITE 1 NFILE {}
          * title
          *
          '''.format(dcd_unit,self.nrun))
      def close_files(self):
          if rank==0:
             self.exch_all.close()
          self.exch.close()
          self.hist.close()
          self.dcd.close()
          comm.barrier()
      def write_trj(self,run):
          lingo.charmm_script('traj write')
          self.ener_im=lingo.get_energy_value('ENER')
          self.exch.write('%5d %2d %15f '%(run,repl_id,self.ener_im))
          self.hist.write('%i %s \n'%(run,self.lmd))
          comm.barrier()
      def gather_condid(self,run):
          data=self.cond_id   
          data=comm.gather(data,root=0)
          data=comm.bcast(data,root=0)
          #print('run',run,'rank',rank,data)
          return data
      def swap_neighbor(self,run):
          self.new_cond_id=get_neighbor(self.cond_id,run)
          self.neighbor_repl_id=self.cond_data.index(self.new_cond_id)
          lmd=self.cond_array[self.new_cond_id][0]
          return lmd
      def get_ener_new(self):
          blk=set_block(comp,ires='1',nstates='2',lmd=self.lmd_new)
          blk.run()
          energy.show()
          self.ener_in=lingo.get_energy_value('ENER')
          self.exch.write('%2d %15f '%(self.neighbor_repl_id,self.ener_in))
          comm.barrier()
      def metropolis(self):
          self.Accept=True
          self.prob=1
          self.rand=0
          if self.cond_id > self.new_cond_id:
             comm.send(self.ener_im, dest=self.neighbor_repl_id, tag=11)
             comm.send(self.ener_in, dest=self.neighbor_repl_id, tag=12)
          elif self.cond_id < self.new_cond_id:
             ener_jn=comm.recv(source=self.neighbor_repl_id,tag=11)
             ener_jm=comm.recv(source=self.neighbor_repl_id,tag=12)
             delta=-1*beta*(ener_jm+self.ener_in-self.ener_im-ener_jn)
             #print("beta:",beta,"delta:",delta,"ener_jm:",ener_jm,"ener_in:",self.ener_in,"ener_im:",self.ener_im,"ener_jn:",ener_jn)
             if delta<0:
                self.prob=np.exp(delta)
                #print('prob ',self.prob)
                self.rand=np.random.uniform(low=0.0, high=1.0)
                if self.prob < self.rand:
                   self.Accept=False
                   #print('reject')
          comm.barrier()
          if self.cond_id < self.new_cond_id:
             comm.send(self.Accept,dest=self.neighbor_repl_id,tag=13)
          elif self.cond_id > self.new_cond_id:
             self.Accept=comm.recv(source=self.neighbor_repl_id,tag=13)
          comm.barrier()
          self.exch.write('%5f %5f %5s '%(self.prob,self.rand,str(self.Accept)[0]))
      def write_exch_all(self,run):
          # exchange file for each replica
          if self.Accept == True:
             replica_next=self.neighbor_repl_id
          else:
             replica_next=repl_id
          self.exch.write('%2d \n'%(replica_next))
          exch_str='{:2d} {:12.6f} {:15.6f} {:2d} {:12.6f} {:15.6f} {:5.3f} {:5.3f} {:1s} {:2d}\n'.format(repl_id,self.cond_id,self.ener_im,self.neighbor_repl_id,self.new_cond_id,self.ener_in,self.prob,self.rand,str(self.Accept)[0],replica_next)
          # gather all exchange info (a string) to an array
          comm.barrier()
          exch_str_all=comm.gather(exch_str,root=0)
          comm.barrier()
          if rank==0:
             exch_list_all=[]
             for i in range(0,nproc):
                 exch_list_all.append(exch_str_all[i].split())
             exch_all_array=np.array(exch_list_all)
             exch_array_sort=pd.DataFrame(exch_all_array,columns=['repl_id','condid_m','ener_im','neighbor_repl_id','condid_n','ener_in','prob','rand','accept','neighbor_repl_id']).astype({'condid_m':float}).sort_values(by='condid_m',ascending=True).to_numpy()
             self.exch_all.write('# Exchange %15d: STEP %12d: REPEAT     1\n'%(run+1,run+1))
             for i in range(0,nproc):
                 self.exch_all.write('%2d %12f %15f %2d %12f %15f %5.3f %5.3f %1s %2d\n'%
                 (int(exch_array_sort[i,0])+1,float(exch_array_sort[i,1]),float(exch_array_sort[i,2]),
                  int(exch_array_sort[i,3])+1,float(exch_array_sort[i,4]),float(exch_array_sort[i,5]),
                  float(exch_array_sort[i,6]),float(exch_array_sort[i,7]),exch_array_sort[i,8],int(exch_array_sort[i,9])+1))
      def update_cond(self):
          # update condition 
          if self.Accept == False:
             self.new_cond_id=self.cond_id
             self.lmd=self.cond_array[self.new_cond_id][0]
             blk=set_block(comp,ires='1',nstates='2',lmd=self.lmd)
             blk.run()
          else:
             self.lmd=self.lmd_new
          self.cond_id=self.new_cond_id
      def trj_unmixing(self):
          if rank==0:
             #lingo.charmm_script('prnlev 10')
             exch_all=charmm_file.CharmmFile(file_name=exch_all_fn,read_only=True,formatted=True)
             for i in range(0,nproc):
                 unit_inp=7+i
                 unit_out=48+i
                 fn_inp='aa'+str(i)+'/dcd/prod'+str(itt)+'.dcd'
                 os.system('mkdir -p cond'+str(i)+'/dcd')
                 fn_out='cond'+str(i)+'/dcd/prod'+str(itt)+'.dcd'
                 lingo.charmm_script('open read file unit {} name {}'.format(unit_inp,fn_inp))
                 lingo.charmm_script('open writ file unit {} name {}'.format(unit_out,fn_out))
             lingo.charmm_script('merge firstu 7 nunit {} outp 48 RTOTemp excu {} NEXChange {} nrpl {} nrep 1'.format(nproc,exch_all.file_unit,ncycles,nproc))
             exch_all.close
      def accept_ratio(self):
          if rank==0:
             dat=np.loadtxt(exch_all_fn,usecols=(1,4,8),dtype={'names':('i','j','accept'),'formats':(float,float,'|S1')})
             df=pd.DataFrame(dat)
             #print(df)
             df_ex=df.groupby(['i','j']).size().reset_index(name='allex')
             df_ac=df.groupby(['i','j','accept']).size().reset_index(name='success')
             df_ac_t=df_ac[df_ac['accept']==b'T'].reset_index()
             print(df_ex)
             print(df_ac_t)
             npairs=df_ex.shape[0]
             ratio_list=[]
             #print(npairs)
             for i in range(0,npairs):
                 cond_i=df_ex.loc[i,['i']][0]
                 cond_j=df_ex.loc[i,['j']][0]
                 if cond_i < cond_j:
                    if np.sum((df_ac_t['i']==cond_i) & (df_ac_t['j']==cond_j))==0:
                       #print(cond_i,cond_j,'no successful exchange')
                       accepted=0
                    else:
                       accepted=df_ac_t[(df_ac_t['i']==cond_i) & (df_ac_t['j']==cond_j)]['success'].to_numpy()[0]
                    exchanged=df_ex[(df_ex['i']==cond_i) & (df_ex['j']==cond_j)]['allex'].to_numpy()[0]
                    ratio=accepted/exchanged
                    #print(cond_i,cond_j,accepted,exchanged,ratio)
                    ratio_list.append([cond_i,cond_j,accepted,exchanged,ratio])
             ratio_array=np.array(ratio_list)
             #print(ratio_array)
             np.savetxt(ratio_fn,ratio_array,fmt='%.6f')

def get_latest_lmd(cond_arr):
    exch_all_pre=np.loadtxt(exch_all_pre_fn,usecols=(9),dtype=int)
    repid_new=exch_all_pre[-1*nproc:]-1
    l='{0:.4f}'.format(cond_arr[repid_new==rank][0][0])
    condid=np.where(repid_new==rank)[0][0]
    #print(rank+1,intra,inter,condid)
    return l,condid

def calc_ener(cond_array):
    dcd_cond_fn='cond'+str(rank)+'/dcd/prod'+str(itt)+'.dcd'
    os.system('mkdir -p cond'+str(rank)+'/ener')
    for i in range(0,nproc):
        l_new=cond_array[i][0]
        blk=set_block(comp,ires=1,nstates='2',lmd=l_new)
        blk.run()
        #print(l_new)
        ener_fn='cond'+str(rank)+'/ener/ener'+str(itt)+'-'+str(l_new)+'.dat'
        dcd_cond=charmm_file.CharmmFile(file_name=dcd_cond_fn,file_unit=7,read_only=True,formatted=False)
        lingo.charmm_script('traj query unit 7')
        nf=lingo.get_energy_value('NFILE')
        skip=lingo.get_energy_value('SKIP')
        lingo.charmm_script('traj firstu 7 nunit 1 skip %s'%(skip*skip_dcd))
        ener=open(ener_fn,'w',buffering=1)
        settings.set_verbosity(0)
        for iframe in range(0,int(nf/skip_dcd)):
            lingo.charmm_script('''
                       traj read
                       energy
            ''')
            e=lingo.get_energy_value('ENER')
            ener.write('%7i %13.4f\n'%(iframe,e))
        ener.close()
        dcd_cond.close()

def mbar(combine,cond_array):
    if rank==0:
       if combine==True:
          excl=10
          configs=int(itt*ncycles/skip_dcd)-excl
          beg=1
       else:
          excl=10
          configs=int(ncycles/skip_dcd)-excl
          beg=itt
       ener=np.zeros((nproc,configs,nproc), dtype = float)
       for i in range(0,nproc):
           # protein calc energies loaded in steps
           for step in range(0,nproc):
               l_new=cond_array[step][0]
               #fileName = "cond"+str(i)+"/ener/ener-"+str(l_new)+".dat"
               #tmp = np.loadtxt(fileName,dtype=float,usecols=(1),skiprows=0,max_rows=configs)
               for ii in range(beg,itt+1):
                   fileName = "cond"+str(i)+"/ener/ener"+str(ii)+"-"+str(l_new)+".dat"
                   tmpi = np.loadtxt(fileName,dtype=float,usecols=(1),skiprows=0,max_rows=ncycles)
                   if ii>beg:
                      tmp=np.concatenate((tmp,tmpi))
                   elif ii==beg:
                      tmp=tmpi
               if step > 0:
                   tmp2 = np.column_stack((tmp2,tmp))
               else:
                   tmp2=np.copy(tmp)
           ener[i,:,:] = tmp2[excl:,:]

       # calculate free energy difference using FastMBAR
       N_k = [configs for i in range(nproc)]
       N_k = np.asarray(N_k)
       tmp = ener.reshape((-1,nproc))
       u_kn = tmp.T * beta
       print(u_kn)
       print(N_k)
       mb = FastMBAR.FastMBAR(energy = u_kn, num_conf = N_k, cuda = False, bootstrap=True, verbose=True)
       F = mb.F
       F_std = mb.F_std
       print("F is:",F)
       print("F_std is:",F_std)

       # write to disk
       fp=open(fe_fn,'w')
       for i in range(len(F)):
           l=cond_array[i][0]
           # convert reduced energies back to regular energies
           fp.write("%.6f %.6f %.6f\n" % (l,F[i]/beta,F_std[i]/beta))
       fp.close()



def main():
    cond_arr=write_cond_file()
    if itt >= 2:
       l,condid=get_latest_lmd(cond_arr)
    else:
       condid=rank
       l=l_g
    charmm_log=open_clog()
    read_param()
    read_init()
    set_pbc()
    define_sub(comp,ires=1)
    blk=set_block(comp,ires=1,nstates='2',lmd=l)
    blk.run()
    setup_nb(nbond)

    # NOE and IC restraints don't work in omm dynamics, although it works in omm energy.
    # b/c in energy calculation, if a term is not present in omm, charmm will be used to compute this term.
    restrain=resd_restrain(comp,ires=1)
    restrain.run()

    if itt == 1:
       mini()
    dyn_init()
    rex=rep_ex(ncycle=ncycles,condid=condid,lmd=l,conditions=cond_arr)
    rex.run()
    close_clog(charmm_log)
    rex.trj_unmixing()
    comm.barrier()
    rex.accept_ratio()
    calc_ener(cond_arr)
    comm.barrier()
    mbar(combine=True,cond_array=cond_arr)
    comm.barrier()

if __name__ == '__main__':
   main()
