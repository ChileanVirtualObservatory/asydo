import math


SPEED_OF_LIGHT = 299792458.0
KILO = 1000

class Component:
    """Abstract Component Model"""
    def __init__(self,log,z_base=0.0):
        self.log=log
        self.z=z_base 
        self.rv=SPEED_OF_LIGHT/KILO*((self.z**2 + 2*self.z)/(self.z**2 + 2*self.z + 2))

    def setRadialVelocity(self,rvel):
        """Set radial velocity in km/s"""
        self.rv=rvel
        self.z=math.sqrt((1 + self.rv*KILO/SPEED_OF_LIGHT)/(1 - self.rv*KILO/SPEED_OF_LIGHT)) -1

    def info(self):
        return "(none)"
    
    def register(self,comp_name,alpha,delta):
        self.comp_name=comp_name
        self.alpha=alpha
        self.delta=delta

    def project(self,cube):
        pass

