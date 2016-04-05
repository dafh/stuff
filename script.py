import subprocess
from SCRIPT import dir_path
import glob
from astropy.io import fits
import os
import numpy as np



abspath = os.path.abspath(__file__)


#Remember last inputs!
if os.path.isfile('.dirname'):
    with open('.dirname','r') as out:
         dirnames = out.readlines()
    dirnames = ' '.join(dirnames).replace('\n','').split()

else:
    dirnames = ['']*6

if os.path.isfile('.inputs'):
    with open('.inputs','r') as out:
        inputs = out.readlines()
    inputs = ' '.join(inputs).replace('\n','').split()
else:
    inputs = ['']*6

bias_folder = raw_input('BIAS_FOLDER (/%s): ' % dirnames[0])
if bias_folder != '':
   dirnames[0] = bias_folder
bias_path = os.getcwd() + '/' + dirnames[0] + '/'
dir_path(bias_path)
subprocess.call('rm BIAS*',shell=True)
image_bias = glob.glob('FORS*.fits')
for image in image_bias:
    hdulist = fits.open(image)
    header = hdulist[0].header
    axis2 = header['NAXIS2']
    archive = 'BIAS_%sp' % axis2
    if os.path.isfile(archive):
       outfile = open(archive,'a')
    else:
        outfile = open(archive,'w')
    outfile.write(image + '\n')
    outfile.close()
    hdulist.close()

BIAS = glob.glob('BIAS*')
print BIAS
#if len(BIAS) == 1:
#    inputs[0]=BIAS[0]
ask = raw_input('which to use?(%s): ' % inputs[0])
if ask != '':
   inputs[0]= ask
bias = inputs[0]
os.chdir(os.path.dirname(abspath))

flat_folder = raw_input('FLAT FOLDER (/%s): ' % dirnames[1])
if flat_folder != '':
   dirnames[1] = flat_folder
flat_path = os.getcwd() + '/' + dirnames[1] + '/'
dir_path(flat_path)
subprocess.call('rm FLAT*',shell=True)
image_flat = glob.glob('FORS*.fits')
for image in image_flat:
    hdulist = fits.open(image)
    header = hdulist[0].header
    if header['OBJECT'].find('FLAT') != -1:
       axis2 = header['NAXIS2']
       time = int(header['EXPTIME'])
       archive = 'FLAT_%ss_%sp' % (time,axis2)
       if os.path.isfile(archive):
          outfile = open(archive,'a')
       else:
          outfile = open(archive,'w')
       outfile.write(image + '\n')
       outfile.close()
    hdulist.close()

FLAT = glob.glob('FLAT*')
print FLAT
#if len(FLAT) == 1:
#    inputs[1]=FLAT[0]
ask = raw_input('which to use?(%s): ' % inputs[1])
if ask != '':
   inputs[1]= ask
flat = inputs[1]
os.chdir(os.path.dirname(abspath))


std_folder = raw_input('STD FOLDER (/%s): ' % dirnames[2])
if std_folder != '':
   dirnames[2] = std_folder
std_path = os.getcwd() + '/' + dirnames[2] + '/'
dir_path(std_path)
subprocess.call('rm STD*',shell=True)
image_std = glob.glob('FORS*.fits')
for image in image_std:
    hdulist = fits.open(image)
    header = hdulist[0].header
    if header['OBJECT'].find('STD') != -1:
       axis2 = header['NAXIS2']
       time = int(header['EXPTIME'])
       archive = 'STD_%ss_%sp' % (time,axis2)
       if os.path.isfile(archive):
          outfile = open(archive,'a')
       else:
          outfile = open(archive,'w')
       outfile.write(image + '\n')
       outfile.close()
    hdulist.close()

STD = glob.glob('STD*')
print STD
#if len(STD) == 1:
#    inputs[2]=STD[0]
ask = raw_input('which to use?(%s): ' % inputs[2])
if ask != '':
   inputs[2]= ask
std = inputs[2]
os.chdir(os.path.dirname(abspath))


wave_folder = raw_input('WAVE FOLDER (/%s): ' % dirnames[3])
if wave_folder != '':
   dirnames[3] = wave_folder
wave_path = os.getcwd() + '/' + dirnames[3] + '/'
dir_path(wave_path)
subprocess.call('rm WAVE*',shell=True)
image_wave = glob.glob('FORS*.fits')
for image in image_wave:
    hdulist = fits.open(image)
    header = hdulist[0].header
    if header['OBJECT'].find('WAVE') != -1:
       axis2 = header['NAXIS2']
       time = int(header['EXPTIME'])
       archive = 'WAVE_%ss_%sp' % (time,axis2)
       if os.path.isfile(archive):
          outfile = open(archive,'a')
       else:
          outfile = open(archive,'w')
       outfile.write(image + '\n')
       outfile.close()
    hdulist.close()

WAVE = glob.glob('WAVE*')
print WAVE
#if len(WAVE) == 1:
#    inputs[3]=WAVE[0]
ask = raw_input('which to use?(%s): ' % inputs[3])
if ask != '':
   inputs[3]= ask
wave = inputs[3]
os.chdir(os.path.dirname(abspath))

flux_folder = raw_input('FLUX FOLDER (/%s): ' % dirnames[4])
if flux_folder != '':
   dirnames[4] = flux_folder
flux_path = os.getcwd() + '/' + dirnames[4] + '/'
dir_path(flux_path)
subprocess.call('rm flux*',shell=True)
image_flux = glob.glob('FORS*.fits')
for image in image_flux:
    hdulist = fits.open(image)
    header = hdulist[0].header
    if header['OBJECT'].find('TBD') != -1 and header['OBJECT'].find('SN') == -1:
       axis2 = header['NAXIS2']
       time = int(header['EXPTIME'])
       archive = 'flux_%ss_%sp' % (time,axis2)
       if os.path.isfile(archive):
          outfile = open(archive,'a')
       else:
          outfile = open(archive,'w')
       outfile.write(image + '\n')
       outfile.close()
    hdulist.close()

FLUX = glob.glob('flux*')
print FLUX
#if len(FLUX) == 1:
#    inputs[4]=FLUX[0]
ask = raw_input('which to use?(%s): ' % inputs[4])
if ask != '':
   inputs[4]= ask
flux = inputs[4]
os.chdir(os.path.dirname(abspath))

science_folder = raw_input('SCIENCE FOLDER (/%s): ' % dirnames[5])
if science_folder != '':
   dirnames[5] = science_folder
science_path = os.getcwd() + '/' + dirnames[5] + '/'
dir_path(science_path)
subprocess.call('rm SN*',shell=True)
image_science = glob.glob('FORS*.fits')
for image in image_science:
    hdulist = fits.open(image)
    header = hdulist[0].header
    if header['OBJECT'].find('SN') != -1:
       axis2 = header['NAXIS2']
       time = int(header['EXPTIME'])
       archive = 'SN_%ss_%sp' % (time,axis2)
       if os.path.isfile(archive):
          outfile = open(archive,'a')
       else:
          outfile = open(archive,'w')
       outfile.write(image + '\n')
       outfile.close()
    hdulist.close()

SCIENCE = glob.glob('SN*')
print SCIENCE
#if len(SCIENCE) == 1:
#    inputs[5]=SCIENCE[0]
ask = raw_input('which to use?(%s): ' % inputs[5])
if ask != '':
   inputs[5]= ask
science = inputs[5]
os.chdir(os.path.dirname(abspath))

dirnames = np.asarray(dirnames)
np.savetxt('.dirname',dirnames,delimiter=' ',fmt='%10s')
np.savetxt('.inputs',inputs,delimiter=' ',fmt='%10s')
subprocess.call('python SCRIPT.py bias=%s flat=%s std=%s wave=%s flux=%s science=%s' % (bias_path + bias,flat_path + flat,std_path + std,wave_path + wave,flux_path + flux,science_path + science),shell=True)
