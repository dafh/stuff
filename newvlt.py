#newvlt.py
import sys
import os
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

#inputs
class image:
      def __init__(self, name):
          self.name = name
          hdulist = fits.open(self.name)
          self.header = hdulist[0].header
          self.data = hdulist[0].data
          if name.find('sens')==-1:
             wat0 = hdulist[0].header['WAT2_001'].split('spec1 = "')[1].split()
             wat1 = hdulist[0].header['WAT2_002'].split('spec2 = "')[1].split()
             self.len0 = float(wat0[5])
             self.len1 = float(wat1[5])
             self.beam0bin = float(wat0[4])
             self.beam1bin = float(wat1[4])
             self.first0 = float(wat0[3])
             self.first1 = float(wat0[3])
          else:
              self.crval = float(hdulist[0].header['CRVAL1'])
              self.crpix = float(hdulist[0].header['CRPIX1'])
              self.cdelt = float(hdulist[0].header['CDELT1'])
          hdulist.close()


      def get_key(self,key):
          keys = {}
          for keyword in self.header:
              aux = keyword.split('=')[0]
              if aux.find(key) != -1:
                 keys[aux] = self.header[aux]
                 new_key = aux
          return keys, new_key

def inputs ():
    SN = []
    sens = []
    dictionary = {}
    if (len(sys.argv)!=8):
       print 'USAGE: key=arg; data, qispMW, uqispMW, wispMW, redshiftkps, AvMW, RvMW as keywords'
       sys.exit()
    else:
        for i in range(len(sys.argv)):
            if i>0:
               aux =  sys.argv[i].split('=')
               key = aux[0]
               arg = aux[1]
               dictionary[key] = arg
        SNfile = dictionary['data']
        if os.path.isfile(SNfile):
           inputs = open(SNfile,'r')
           for line in inputs:
                image = line.split('\n')[0]
                if os.path.isfile(image):
                   if line.find('sens')== -1:
                      SN.append(image)
                   else:
                      sens.append(image)
                else:
                    print 'image %s not found' % image
           inputs.close()
        else:
            print '%s not found' % SNfile
            sys.exit()
        return SN, sens, dictionary


if __name__ == '__main__':
   SN, sens, dictionary = inputs()
   wispMW = float(dictionary['wispMW'])/10000.
   #We check if POSANG angle is the same
   SN_im = []
   sens_im = []
   exptimes = []
   for sn in SN:
       aux = image(sn)
       pos_key, kpos = aux.get_key('WOLL POSANG')
       exp_key, kexp = aux.get_key('EXPTIME')
       exptimes.append(exp_key[kexp])
       if pos_key[kpos] >= 0.05:
          print 'bad POSANG'
          sys.exit()
       SN_im.append(aux)
   for se in sens:
       aux = image(se)
       sens_im.append(aux)

   file0 = SN_im[0].data
   file45 = SN_im[1].data
   file22 = SN_im[2].data
   file67 = SN_im[3].data

   #Para imitar el codigo de IDL no mas
   beam0bin = np.array([SN_im[0].beam0bin, SN_im[1].beam0bin, SN_im[2].beam0bin, SN_im[3].beam0bin])
   beam1bin = np.array([SN_im[0].beam1bin, SN_im[1].beam1bin, SN_im[2].beam1bin, SN_im[3].beam1bin])
   len0 = np.array([SN_im[0].len0, SN_im[1].len0, SN_im[2].len0, SN_im[3].len0])
   len1 = np.array([SN_im[0].len1, SN_im[1].len1, SN_im[2].len1, SN_im[3].len1])
   first0 = np.array([SN_im[0].first0, SN_im[1].first0, SN_im[2].first0, SN_im[3].first0])
   first1 = np.array([SN_im[0].first1, SN_im[1].first1, SN_im[2].first1, SN_im[3].first1])


   lam00 = first0[0] + beam0bin[0]*(np.arange(int(len0[0])) + 1)
   lam01 = first1[0] + beam1bin[0]*(np.arange(int(len1[0])) + 1)
   lam450 = first0[1] + beam0bin[1]*(np.arange(int(len0[1])) + 1)
   lam451 = first1[1] + beam1bin[1]*(np.arange(int(len1[1])) + 1)
   lam220 = first0[2] + beam0bin[2]*(np.arange(int(len0[2])) + 1)
   lam221 = first1[2] + beam1bin[2]*(np.arange(int(len1[2])) + 1)
   lam670 = first0[3] + beam0bin[3]*(np.arange(int(len0[3])) + 1)
   lam671 = first1[3] + beam1bin[3]*(np.arange(int(len1[3])) + 1)

   key_gain, gaux = SN_im[0].get_key('OUT1 GAIN')
   key_ron, ronaux = SN_im[0].get_key('OUT1 RON')

   GAIN = float(key_gain[gaux])
   RON = float(key_ron[ronaux])

   n_pix = len(lam00)
   fluxindex = 0
   skyindex = 1

   # matrix with elements arrays of len = 1833, following  wang code: beam index (object), beam index (sky), waveplate index
   f = np.empty([n_pix,2,2,4])
   # angle 0 ap 1
   f[:,0,0,0] = file0[0,fluxindex]
   f[:,0,1,0] = file0[0,skyindex]
   #angle 0 ap 2
   if np.array_equal(lam00,lam01):
      f[:,1,0,0] = file0[1,fluxindex]
      f[:,1,1,0] = file0[1,skyindex]
   else:
       f[:,1,0,0] = interpolate.interp1d(lam01,file0[1,fluxindex],kind='linear',bounds_error=False,fill_value='extrapolate')(lam00)
       f[:,1,1,0] = interpolate.interp1d(lam01,file0[1,skyindex],kind='linear',bounds_error=False,fill_value='extrapolate')(lam00)
   # angle22 ap 1 & 2
   if np.array_equal(lam00,lam220):
       f[:,0,0,1] = file22[0,fluxindex]
       f[:,0,1,1] = file22[0,skyindex]
   else:
       f[:,0,0,1] = interpolate.interp1d(lam220,file22[0,fluxindex],kind='linear',bounds_error=False,fill_value='extrapolate')(lam00)
       f[:,0,1,1] = interpolate.interp1d(lam220,file22[0,skyindex],kind='linear',bounds_error=False,fill_value='extrapolate')(lam00)

   if np.array_equal(lam00,lam221):
        f[:,1,0,1] = f22[1,fluxindex]
        f[:,1,1,1] = f22[1,skyindex]
   else:
        f[:,1,0,1] = interpolate.interp1d(lam221,file22[1,fluxindex],kind='linear',bounds_error=False,fill_value='extrapolate')(lam00)
        f[:,1,1,1] = interpolate.interp1d(lam221,file22[1,skyindex],kind='linear',bounds_error=False,fill_value='extrapolate')(lam00)
    # angle45 ap 1 & 2
   if np.array_equal(lam00,lam450):
       f[:,0,0,2] = f45[0,fluxindex]
       f[:,0,1,2] = f45[0,skyindex]
   else:
        f[:,0,0,2] =  interpolate.interp1d(lam450,file45[0,fluxindex],kind='linear',bounds_error=False,fill_value='extrapolate')(lam00)
        f[:,0,1,2] = interpolate.interp1d(lam450,file45[0,skyindex],kind='linear',bounds_error=False,fill_value='extrapolate')(lam00)
   if np.array_equal(lam00,lam451):
       f[:,1,0,2] = f45[1,fluxindex]
       f[:,1,1,2] = f45[1,skyindex]
   else:
        f[:,1,0,2] = interpolate.interp1d(lam451,file45[1,fluxindex],kind='linear',bounds_error=False,fill_value='extrapolate')(lam00)
        f[:,1,1,2] = interpolate.interp1d(lam451,file45[1,skyindex],kind='linear',bounds_error=False,fill_value='extrapolate')(lam00)
   #angle67 ap 1 & 2
   if np.array_equal(lam00,lam670):
       f[:,0,0,3] = f67[0,fluxindex]
       f[:,0,1,3] = f67[0,skyindex]
   else:
      f[:,0,0,3] = interpolate.interp1d(lam670,file67[0,fluxindex],kind='linear',bounds_error=False,fill_value='extrapolate')(lam00)
      f[:,0,1,3] = interpolate.interp1d(lam670,file67[0,skyindex],kind='linear',bounds_error=False,fill_value='extrapolate')(lam00)
   if np.array_equal(lam00,lam671):
	  f[:,1,0,3] = f67[1,fluxindex]
	  f[:,1,1,3] = f67[1,skyindex]
   else:
       f[:,1,0,3] = interpolate.interp1d(lam671,file67[1,fluxindex],kind='linear',bounds_error=False,fill_value='extrapolate')(lam00)
       f[:,1,1,3] = interpolate.interp1d(lam671,file67[1,skyindex],kind='linear',bounds_error=False,fill_value='extrapolate')(lam00)

   #Reminder: valuesspecfile[5] == len(lam00) = n_pix
   #Calculando derivadas
   deriv = np.empty([n_pix,2,4])
   for wave in range(4):
       for beam in range(2):
           aux = ((np.diff(f[:,beam,0,wave])[1:] + np.diff(f[:,beam,0,wave])[:n_pix-2]) + (np.diff(f[:,beam,1,wave])[1:] + np.diff(f[:,beam,1,wave])[:n_pix-2]))/(np.diff(lam00)[1:] + np.diff(lam00)[:n_pix-2])
           aux = np.insert(aux,0,aux[0])
           aux = np.append(aux,aux[n_pix-2])
           deriv[:,beam,wave] = aux

   #Calculando sigma beam
   errlam = 0.1 #default en wang ??
   errbeam = np.sqrt(f[:,:,0,:] + f[:,:,1,:] + deriv[:,:,:]*errlam + RON**2)/np.sqrt(GAIN)

   #sensitivity and flux calibrations
   sens0 = sens_im[0].data
   sens1 = sens_im[1].data
   lamsens0 = sens_im[0].crval + np.arange(len(sens0)) + 1 - sens_im[0].crpix*sens_im[0].cdelt
   lamsens1 = sens_im[1].crval + np.arange(len(sens0)) + 1 - sens_im[1].crpix*sens_im[1].cdelt

   #calibdir? --> CALIB
   if os.path.isdir('/CALIB'):
      print 'wahaha'

   else:
       print 'NO CALIB'
       sys.exit()
   #plt.plot(lam01,file0[1,fluxindex],lam00,f[:,1,0,0])
   #plt.show()
