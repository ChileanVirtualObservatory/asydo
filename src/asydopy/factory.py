# Create several cubes in parallel, using random parameters
from multiprocessing import Pool
from vu import *
from db import *
import math
import random
import numpy as np
import sys
import os
import time
import matplotlib.pyplot as plt
import copy

#SPEED_OF_LIGHT = 299792458.0
#KILO = 1000

class IMCConf:
   def  __init__(self,number,dbpath="ASYDO",mol_list="all",mol_prob=0.3,x_pos=0.0,y_pos=0.0,f_pos=300000,spa_pix=5,spe_pix=100,fov=500,bw=2000,rvel=(150,1000),temp=(50,500),semiaxis=(10,300),fwhm=(10,50),angle=(0,math.pi),rot=(50,500),curtosis=(-5,5)):
      self.rvel=rvel
      self.number=number
      self.mol_prob=mol_prob
      self.mol_list=mol_list
      self.dbpath=dbpath
      self.x_pos=x_pos
      self.y_pos=y_pos
      self.f_pos=f_pos
      self.spa_pix=spa_pix
      self.spe_pix=spe_pix
      self.bw=bw
      self.fov=fov
      self.temp=temp
      self.semiaxis=semiaxis
      self.fwhm=fwhm
      self.rot=rot
      self.angle=angle
      self.curtosis=curtosis
      self.force_list=list()
      self.ban_list=list()
   
   def set_params(self,template):
      self.rvel    =template.rvel
      self.mol_prob=template.mol_prob
      self.mol_list=template.mol_list
      self.dbpath  =template.dbpath
      self.x_pos   =template.x_pos
      self.y_pos   =template.y_pos
      self.f_pos   =template.f_pos
      self.spa_pix =template.spa_pix
      self.spe_pix =template.spe_pix
      self.bw      =template.bw
      self.fov     =template.fov
      self.temp    =template.temp
      self.semiaxis=template.semiaxis
      self.fwhm    =template.fwhm
      self.rot     =template.rot
      self.angle   =template.angle
      self.curtosis=template.curtosis
      self.force_list=copy.deepcopy(template.force_list)
      self.ban_list  =copy.deepcopy(template.ban_list)

def unitary_IMC(conf):
    def rget(val):
         if isinstance(val,tuple):
             return random.uniform(val[0],val[1])
         else:
             return val
    
    print "Generating cube", conf.number
    db=DataBase(conf.dbpath)
    db.connect()
    try:
      os.mkdir("logs")
    except OSError:
      pass
    log=open('logs/exp-c'+str(conf.number)+'.log', 'w')
    univ=Universe(log)
    xpos=rget(conf.x_pos)
    ypos=rget(conf.y_pos)
    univ.createSource('AutoGenCube-'+str(conf.number),xpos,ypos)
    fpos=rget(conf.f_pos)
    bw=rget(conf.bw)
    rv=rget(conf.rvel)
    lf=(fpos - bw/2.0)*math.sqrt((1 + rv*KILO/SPEED_OF_LIGHT)/(1 - rv*KILO/SPEED_OF_LIGHT))
    uf=(fpos + bw/2.0)
    if (conf.molecules=="all"):
        mlist=db.getMoleculeList(lf,uf)
        chList=list()
        for mol in mlist:
           chList.append(mol[0])
    else:
        chList=conf.molecules
    temp=rget(conf.temp)
    #print chList
    # HERE Random selection of molecules
    for chName in chList:
       if chName in conf.ban_list:
          log.write("Mol: "+chName[0]+" banned!")
          continue
       if not (chName in conf.force_list):
          if random.random() > conf.mol_prob:
             continue
       molist=db.getSpeciesList(chName,lf,uf)
       s_x=rget(conf.semiaxis)
       s_y=rget(conf.semiaxis)
       angle=rget(conf.angle)
       rot=rget(conf.rot)
       fw=rget(conf.fwhm)
       curt=rget(conf.curtosis)
       mms=""
       for mol in molist:
          if mms!="":
             mms+=","
          mms+=str(mol)
       model=IMCM(log,conf.dbpath,mol,temp,('normal',s_x,s_y,angle),('skew',fw,curt),('linear',angle,rot))
       model.set_radial_velocity(rv)
       univ.add_component('AutoGenCube-'+str(n),model)
    fov=rget(conf.fov)
    bw=rget(conf.bw)
    cube=univ.gen_cube('AutoGenCube-'+str(n),xpos,ypos,fpos,fov/conf.spa_pix,fov,bw/conf.spe_pix,bw)
    db.disconnect()
    log.close()
    return cube

def unitary_IMC_Wrap(n):
   try:
      return(unitaryGen(n))
   except:
      print('%d: %s' % (n, traceback.format_exc()))

def gen_IMC_cubes(confs,nproc):
   p = Pool(len(confs))
   result=p.map(unitary_IMC_Wrap,confs)
   return result

#sample_size=30000
#CFParams.ban_list.append('Phosphapropynylidyne')
#negative=np.asarray(parallelGen(sample_size/2,4))
#CFParams.ban_list=list()
#CFParams.force_list.append('Phosphapropynylidyne')
#positive=np.asarray(parallelGen(sample_size/2,4))
#result=np.vstack((negative,positive))
#np.save('exp1.npy', result)

#import cProfile
#import re
#cProfile.run('unitaryGen(0)')
