# Create several cubes in parallel, using random parameters
from multiprocessing import Pool
from Universe import *
from IMCM import *
from DataBase import *
import math
import random
import numpy as np
import sys
import os
import time
import matplotlib.pyplot as plt

SPEED_OF_LIGHT = 299792458.0
KILO = 1000

global CFParams


def rget(val):
   if isinstance(val,tuple):
      return random.uniform(val[0],val[1])
   else:
      return val

class CFP:
   def  __init__(self,dbpath,mol_prob=0.3,x_pos=0.0,y_pos=0.0,f_pos=300000,spa_pix=5,spe_pix=100,fov=500,bw=2000,rvel=(150,1000),temp=(50,500),semiaxis=(10,300),fwhm=(10,50),angle=(0,math.pi),rot=(50,500),curtosis=(-5,5)):
      self.rvel=rvel
      self.mol_prob=mol_prob
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

   def forceMolecule(self,name):
      self.force_list.append(name)

   def banMolecule(self,name):
      self.ban_list.append(name)

CFParams=CFP('ASYDO')

def unitaryGen(n):
      global CFParams
      print "Generating cube", n
      db=DataBase(CFParams.dbpath)
      db.connect()
      try:
        os.mkdir("logs")
      except OSError:
        pass
      log=open('logs/cube-'+str(n)+'.log', 'w')
      univ=Universe(log)
      xpos=rget(CFParams.x_pos)
      ypos=rget(CFParams.y_pos)
      univ.createSource('AutoGenCube-'+str(n),xpos,ypos)
      fpos=rget(CFParams.f_pos)
      bw=rget(CFParams.bw)
      rv=rget(CFParams.rvel)
      lf=(fpos - bw/2.0)*math.sqrt((1 + rv*KILO/SPEED_OF_LIGHT)/(1 - rv*KILO/SPEED_OF_LIGHT))
      uf=(fpos + bw/2.0)
      chList=db.getMoleculeList(lf,uf)
      temp=rget(CFParams.temp)
      #print chList
      # HERE Random selection of molecules
      for chName in chList:
         if chName in CFParams.ban_list:
            continue
         if CFParams.mol_prob < random.random() and chName not in CFParams.force_list:
            continue
         molist=db.getSpeciesList(chName[0],lf,uf)
         s_x=rget(CFParams.semiaxis)
         s_y=rget(CFParams.semiaxis)
         angle=rget(CFParams.angle)
         rot=rget(CFParams.rot)
         fw=rget(CFParams.fwhm)
         curt=rget(CFParams.curtosis)
         mms=""
         for mol in molist:
            if mms!="":
               mms+=","
            mms+=str(mol[0])
         model=IMCM(log,mol[0],temp,('normal',s_x,s_y,angle),('skew',fw,curt),('linear',angle,rot))
         model.setRadialVelocity(rv)
         univ.addComponent('AutoGenCube-'+str(n),model)
      fov=rget(CFParams.fov)
      bw=rget(CFParams.bw)
      cspec=CubeSpec(xpos,ypos,fpos,fov/CFParams.spa_pix,fov,bw/CFParams.spe_pix,bw)
      cube=univ.genCube('AutoGenCube-'+str(n),cspec)
      db.disconnect()
      log.close()
      return cube.data.flatten()

def parallelGen(samples,nproc):
   p = Pool(nproc)
   result=p.map(unitaryGen,range(samples))
   return result

      #info = np.asarray(result)
      #np.save('exp1.npy', info)
#cubelst=parallelGen(100,8)
import cProfile
import re
cProfile.run('unitaryGen(0)')
