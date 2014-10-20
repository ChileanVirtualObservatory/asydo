import pickle
import numpy as np
from scipy import signal
import math
from matplotlib import pyplot as plt
from asydopy import factory, vu
import astropy
from mpl_toolkits.mplot3d import Axes3D
from sklearn.cluster import AffinityPropagation
from sklearn.cluster import KMeans
from sklearn import metrics

fwhm_bords=(0.8, 10)
f_pos=290000
bw=1000.0
spa_pix=8
spe_pix=500
wlvl=6
wstp=10

dbpath = "../../ASYDO"
template = factory.IMCConf(0, dbpath,
                           mol_list="all",
                           mol_prob=0.3,
                           x_pos=0.0,
                           y_pos=0.0,
                           f_pos=f_pos,
                           spa_pix=spa_pix,
                           spe_pix=spe_pix,
                           fov=500,
                           bw=bw,
                           rvel=(100, 1000),
                           temp=(50, 500),
                           semiaxis=(10, 300),
                           fwhm=fwhm_bords,
                           angle=(0, math.pi),
                           rot=(1,10),
                           curtosis=(-3, 3))
def unique(seq):
    seen = list()
    seen_add = seen.append
    return [ x for x in seq if not (x in seen or seen_add(x))]

def extract_data(cube):
    real_lines=[]
    real_fw=[]
    real_temps=[]
    real_imgs=[]
    real_rv=[]
    real_mol=[]
    for hdu in cube.hdulist:
        if isinstance(hdu,astropy.io.fits.hdu.table.BinTableHDU):
            print hdu.data
            real_mol.append(hdu.data.field(2).tolist())
            real_lines.append(hdu.data.field(3).tolist())
            real_rv.append(hdu.data.field(4).tolist())
            real_fw.append(hdu.data.field(5).tolist())
            real_temps.append(hdu.data.field(6).tolist())
        if isinstance(hdu,astropy.io.fits.hdu.image.ImageHDU):
            real_imgs.append(hdu.data)
    real_lines=[item for sublist in real_lines for item in sublist]
    real_fw=[item for sublist in real_fw for item in sublist]
    real_rv=[item for sublist in real_rv for item in sublist]
    real_temps=[item for sublist in real_temps for item in sublist]
    real_mol=[item for sublist in real_mol for item in sublist]
    return real_lines,real_fw,real_temps,real_imgs,real_rv,real_mol


cube = factory.unitary_IMC_cube(template)
real_lines,real_fw,real_temps,real_imgs,real_rv,real_mol=extract_data(cube)
mols=dict()
mcon=0
#print real_mol
for i in range(len(real_lines)):
   if real_mol[i] not in mols.keys():
      st=real_mol[i]
      mols[st]=mcon
      mcon+=1
   real_lines[i]=vu.freq_correct(real_lines[i],real_rv[i])
   real_fw[i]=vu.S_FACTOR*vu.fwhm2sigma(real_lines[i],real_fw[i])

for i in range(len(real_mol)):
   real_mol[i]=mols[real_mol[i]]


cs=len(real_imgs)/6+1
col=3
if cs > 3:
  col=cs
  cs=3

#aplt.imshow(real_imgs[0].transpose())
#plt.axis('off')
#X,Y=np.meshgrid(np.arange(len(real_imgs[0])),np.arange(len(real_imgs[1][0])),sparse=False, indexing='xy')
#plt.scatter(X.flatten(),Y.flatten(),marker='x',s= 100*real_imgs[1].flatten(),c='w')
#plt.scatter(X.flatten(),Y.flatten(),marker='.',s=-100*real_imgs[1].flatten(),c='k')
#plt.tight_layout()
#plt.show()

for i in range(0,len(real_imgs),2):
   plt.subplot(col,cs,i/2+1)
   plt.imshow(real_imgs[i].transpose())
   plt.axis('off')
   X,Y=np.meshgrid(np.arange(len(real_imgs[i+1])),np.arange(len(real_imgs[i+1][0])),sparse=False, indexing='xy')
   plt.scatter(X.flatten(),Y.flatten(),marker='x',s= 5*real_imgs[i+1].flatten(),c='w')
   plt.scatter(X.flatten(),Y.flatten(),marker='.',s=-10*real_imgs[i+1].flatten(),c='k')
plt.tight_layout()
plt.show()


# Analysis


sigma_low=vu.fwhm2sigma(f_pos-bw/2.0,fwhm_bords[0])
sigma_up=vu.fwhm2sigma(f_pos+bw/2.0,fwhm_bords[1])
ini=1*sigma_low*float(spe_pix)/bw
end=10*sigma_up*float(spe_pix)/bw
sigma_base=sigma_low
twopisig=np.sqrt(2*np.pi)*sigma_base
#print sigma_low,sigma_up
#print ini,end
candidates=[]
ll=np.linspace(ini,end,wlvl)
X=[[] for i in range(wlvl)]
Y=[[] for i in range(wlvl)]
F=[[] for i in range(wlvl)]
A=[[] for i in range(wlvl)]
ALL=[[] for i in range(wlvl)]
ALL=[]
#[[] for i in range(wlvl)]

for y in range(spa_pix):
   candidates_y=[]
   for x in range(spa_pix):
      spectra = cube.data[:,x,y]
      peakind=[]
      for i in range(wlvl-1):
          dt=signal.find_peaks_cwt(spectra, np.linspace(ll[i],2*ll[i+1],wstp),min_snr=2)
          for peak in dt:
             amp=spectra[peak]
             if amp < 0.1:
                continue
             fr=cube.freq_axis[peak]
             X[i].append(x)
             Y[i].append(y)
             F[i].append(fr)
             A[i].append(amp)
             ALL.append([x,y,fr,amp,i])
          peakind.append(dt)
      #print x,y
      #print peakind
      candidates_y.append(peakind)
   candidates.append(candidates_y)

n_clusters=10
#for i in range(wlvl-1):
#   print i
DAT=np.array(ALL)
#print ALL
labels=KMeans(init='k-means++', n_clusters=n_clusters).fit_predict(DAT) 
#af = AffinityPropagation().fit(DAT)
#cluster_centers_indices = af.cluster_centers_indices_
#labels = af.labels_
#n_clusters = len(cluster_centers_indices)
print('Estimated number of clusters: %d' % n_clusters)
#print("Homogeneity: %0.3f" % metrics.homogeneity_score(labels_true, labels))
#print("Completeness: %0.3f" % metrics.completeness_score(labels_true, labels))
#print("V-measure: %0.3f" % metrics.v_measure_score(labels_true, labels))
#print("Adjusted Rand Index: %0.3f"
#      % metrics.adjusted_rand_score(labels_true, labels))
#print("Adjusted Mutual Information: %0.3f"
#      % metrics.adjusted_mutual_info_score(labels_true, labels))
#print("Silhouette Coefficient: %0.3f"
#      % metrics.silhouette_score(X, labels, metric='sqeuclidean'))



   #ax = fig.add_subplot(3,n_clusters/3,i+1)#, projection='3d')
for i in range(n_clusters):
   fig = plt.figure()
   img=np.zeros((spa_pix,spa_pix))
   count=0
   for lab in labels:
      if lab==i:
        img[int(DAT[count,1]),int(DAT[count,0])]+=DAT[count,3]
      count+=1
   plt.imshow(img)
   plt.show()
   
#ax = fig.add_subplot(1,1,1, projection='3d')
#ax.scatter(X[0], Y[0], F[0], c=labels)


   #ax.set_xlabel('X Label')
   #ax.set_ylabel('Y Label')
   #ax.set_zlabel('Z Label')


#print candidates




#samp=len(y_real)
#pks=len(peakind)
#means=x[peakind]
#sigmas=sigma_base*np.ones(pks)
#amps=y_real[peakind]/twopisig
#amps=amps.clip(min=0)
#
#print "peaks GS"
#print pks 
#
#print "real GS"
#print len(real_lines)
#
#
#modlist=[]
#
#for i in range(pks):
#   gauss = GaussianModel(prefix='g'+str(i)+'_')
#   if i==0:
#      pars=gauss.make_params()
#   else:
#      pars.update(gauss.make_params())
#   pars['g'+str(i)+'_center'].set(means[i],min=f_pos-bw/2.0,max=f_pos+bw/2.0)
#   #pars['g'+str(i)+'_center'].set(means[i])
#   pars['g'+str(i)+'_sigma'].set(sigmas[i], min=0.0001, max=10*sigma_up)
#   pars['g'+str(i)+'_amplitude'].set(amps[i], min=0)
#   if i==0:
#      mod=gauss
#   else:
#      mod = mod  + gauss
#
##print pars
#
#out = mod.fit(y_real,pars, x=x)
#fq=[]
#tp=[]
#fw=[]
#pdic=out.params.valuesdict()
#for i in range(pks):
#   prefix='g'+str(i)+'_'
#   fq.append(pdic['g'+str(i)+'_center'])
#   fww=pdic['g'+str(i)+'_fwhm']
#   fw.append(fww)
#   #tp.append(pdic['g'+str(i)+'_amplitude']/(np.sqrt(2*np.pi)*fww/vu.S_FACTOR))
#   tp.append(pdic['g'+str(i)+'_amplitude'])
#
#fig, axarr = plt.subplots(4, sharex=True)
#axarr[0].plot(x, y_real, label='Real Data')
#areaor = np.pi * 30.0 *np.array(real_fw)
#areafit = np.pi * 30.0 *np.array(fw)
#axarr[1].scatter(real_lines,real_temps,c=real_mol,s=areaor,alpha=0.7)
#axarr[2].plot(x, out.best_fit, 'g', label='Fitted')
#axarr[3].scatter(fq,tp,s=areafit,alpha=0.7)
##axarr[3].imshow(matr,interpolation='none',origin='lower',extent=[template.f_pos - template.bw / 2,template.f_pos + template.bw / 2,ini,5*end])
##axarr[3].autoscale(False)
##axarr[3].set_adjustable('box-forced')
##plt.colorbar(im2)
#plt.legend()
#plt.show()
#
#
#
#
#
#
