import numpy as np
from scipy import signal
import math
from matplotlib import pyplot as plt
from asydo import factory, vu
import astropy

al_axis= np.arange(-2,2,0.01)
de_axis= np.arange(-2,2,0.01)
freq_axis= np.arange(50,150,0.1)
S1=np.zeros((400,400))
S2=np.zeros((400,400))
S,bord=vu.gen_surface(('normal',3000.0,500.0,0.1), 0.0, 0.0,al_axis,de_axis)
for x in range(bord[0][0],bord[0][1]):
   for y in range(bord[1][0],bord[1][1]):
      S1[x,y]=S[x-bord[0][0],y-bord[1][0]]
S,bord=vu.gen_surface(('normal',2990.0,499.0,0.1), 0.0, 0.0,al_axis,de_axis)
for x in range(bord[0][0],bord[0][1]):
   for y in range(bord[1][0],bord[1][1]):
      S2[x,y]=0.9952*S[x-bord[0][0],y-bord[1][0]]
plt.imshow(S1-S2)
plt.show()
print freq_axis
[distro,win]=vu.gen_line(('skew',50000.0,6), 100, freq_axis)
print distro
print win
print len(distro)
print len(freq_axis[win[0]:win[1]])
plt.plot(freq_axis[win[0]:win[1]+1],distro)
plt.show()
