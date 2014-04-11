import sys
import urllib
import urllib2
import requests
<<<<<<< HEAD
import atpy


slap_serv='https://find.nrao.edu/splata-slap/slap'
#alma_band_freq={'3':[88,116],'4':[125,163],'6':[211,275],'7':[275,373],'8':[385,500],'9':[602,720]}
alma_band_freq={'7':[275,373]}
segment_freq={'0':[1,50],'1':[50,100],'2':[100,150],'3':[150,200],'4':[200,250],'5':[250,300],'6':[300,350],'7':[350,400],'8':[400,450],\
				  '9':[450,500],'10':[500,550],'11':[550,600],'12':[600,650],'13':[650,700],'14':[700,750],'15':[750,800],'16':[800,850],\
				  '17':[850,900],'18':[900,950],'19':[950,1000]}

def getLines(log,prefix,freqs):
   c=299792458.0
   for band in freqs:
      log.write('band='+band+'\n')
      w_init=c/(freqs[band][0]*1000000000.0);
      w_end= c/(freqs[band][1]*1000000000.0);
      data = '?REQUEST=queryData&WAVELENGTH='+str(w_init)+'/'+str(w_end)+'&VERB=3'
      curl=slap_serv + data.encode('utf-8')
      log.write('  -> Downloading lines via SLAP:\n')
      log.write('  -> '+curl+'\n')
      req = urllib2.Request(curl)
      response = urllib2.urlopen(req)
      votable = response.read()
      stri = prefix+band
      location = './votables/'+stri+'.xml'
      f = open(location, 'w')
      f.write(votable)
      f.close()
=======

slap_serv = 'https://find.nrao.edu/splata-slap/slap'
alma_band_freq = {'3': [88, 116], '4': [125, 163], '6': [211, 275], '7': [275, 373], '8': [385, 500], '9': [602, 720]}
segment_freq = {'0': [1, 50], '1': [50, 100], '2': [100, 150], '3': [150, 200], '4': [200, 250], '5': [250, 300],
                '6': [300, 350], '7': [350, 400], '8': [400, 450], \
                '9': [450, 500], '10': [500, 550], '11': [550, 600], '12': [600, 650], '13': [650, 700],
                '14': [700, 750], '15': [750, 800], '16': [800, 850], \
                '17': [850, 900], '18': [900, 950], '19': [950, 1000]}


def getLines(log, prefix, freqs):
    c = 299792458.0
    for band in freqs:
        log.write('band=' + band + '\n')
        w_init = c / (freqs[band][0] * 1000000000.0);
        w_end = c / (freqs[band][1] * 1000000000.0);
        data = '?REQUEST=queryData&WAVELENGTH=' + str(w_init) + '/' + str(w_end) + '&VERB=3'
        curl = slap_serv + data.encode('utf-8')
        log.write('  -> Downloading lines via SLAP:\n')
        log.write('  -> ' + curl + '\n')
        req = urllib2.Request(curl)
        response = urllib2.urlopen(req)
        votable = response.read()
        stri = prefix + band
        location = './votables/' + stri + '.xml'
        f = open(location, 'w')
        f.write(votable)
        f.close()

>>>>>>> 4b41c0edde51a499fa9065ca36401a7d7a5001b0


#def getAllLines(log):
#   data = '?REQUEST=queryData&WAVELENGTH=0.0/1.0&VERB=3'
#   curl=slap_serv + data.encode('utf-8')
#   log.write('  -> Downloading lines via SLAP:\n')
#   log.write('  -> '+curl+'\n')
#   req = urllib2.Request(curl)
#   response = urllib2.urlopen(req)
#   votable = response.read()
#   location = './votables/all.xml'
#   f = open(location, 'w')
#   f.write(votable)
#   f.close()



getLines(sys.stdout, 'band', alma_band_freq)
#getLines(sys.stdout,'segment',segment_freq)

