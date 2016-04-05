import sys
import os
from PIL import Image
import numpy as np
from astropy.io import fits
import collections
import time
import subprocess
import pyfits
import glob.glob

from pyraf import iraf
from pyraf.iraf import imred, ccdred, noao, images, tv, stsdas, twodspec, imutil, onedspec, longslit, kpnoslit, plot


def ds9_calling():
    subprocess.call('ds9 -geometry 2048x2048 -tile &',shell=True)
    time.sleep(5)

def header(image):
    hdulist = fits.open(image)
    hdulist.info()
    hdulist.close()

def spectra(image):
    leer = open(image,'r')
    escribir = open('ms.' + image,'w')
    for line in leer:
        line = line.split('\n')
        linea = line[0].split('fits')
        escribir.write(linea[0] + 'ms.fits'+'\n')
        iraf.kpnoslit.dispaxi = 1.0
        iraf.kpnoslit.apall.format='multispec'
        iraf.kpnoslit.apall.interact='yes'
        iraf.kpnoslit.apall.find='yes'
        iraf.kpnoslit.apall.recente='yes'
        iraf.kpnoslit.apall.resize='no'
        iraf.kpnoslit.apall.edit='yes'
        iraf.kpnoslit.apall.trace='yes'
        iraf.kpnoslit.apall.fittrac='yes'
        iraf.kpnoslit.apall.extrac='yes'
        iraf.kpnoslit.apall.extras='yes'
        iraf.kpnoslit.apall.reference= ''
        iraf.kpnoslit.apall.review='yes'
        iraf.kpnoslit.apall.nfind= 2.0
        iraf.kpnoslit.apall.nsum = 10
        iraf.kpnoslit.apall.t_func='spline3'
        iraf.kpnoslit.apall.t_order= 2.0
        iraf.kpnoslit.apall.t_niter = 1.0
        iraf.kpnoslit.apall.backgro='fit'
        iraf.kpnoslit.apall.b_funct = 'chebyshev'
        iraf.kpnoslit.apall.b_order = 2.0
        iraf.kpnoslit.apall.t_low_r  = 3.0
        iraf.kpnoslit.apall.t_high = 3.0
        iraf.kpnoslit.apall.readnoi='HIERARCH ESO DET OUT1 RON'
        iraf.kpnoslit.apall.gain = 'HIERARCH ESO DET OUT1 GAIN'
        iraf.kpnoslit.apall.input = line[0]
        iraf.kpnoslit.apall.output = ''
        #iraf.lpar(iraf.kpnoslit.apall)
        iraf.kpnoslit.apall(mode='h')
    leer.close()
    escribir.close()

def wave_spectra(wave,science,wave_folder):


    WAVE = open(wave,'r')
    SCIENCE = open(science,'r')
    w = collections.OrderedDict()
    s = collections.OrderedDict()
    Output = open('w.ms.'+ science,'w')

    REF = 'w.ms.' + SCIENCE.readline().split('\n')[0]

    #REWIND THE FILE!!!!!!!
    SCIENCE.seek(0)


    for line in WAVE:
        line = line.split('\n')
        hdulist = fits.open(wave_folder + line[0])
        w[line[0]] = hdulist[0].header['HIERARCH ESO INS RETA2 ROT']
        hdulist.close()


    for lines in SCIENCE:
        lines = lines.split('\n')
        hdulist = fits.open(lines[0])
        s[lines[0]] = hdulist[0].header['HIERARCH ESO INS RETA2 ROT']
        hdulist.close()



    cont = 1
    for sciences in s:
        #print sciences
        for waves in w:
            if s[sciences] == w[waves]:
               if cont == 1:
                  ref = ''
               else:
                  ref=REF
               print waves, sciences, ref
               output = wave_extraction(wave_folder + waves,sciences)
               Output.write(output + '\n')
               identify(output,ref)
               cont+=1
               break


    WAVE.close()
    SCIENCE.close()
    Output.close()
    return ref

def wave_extraction(wave,ref_science):
    iraf.kpnoslit.dispaxi = 1.0
    iraf.kpnoslit.apall.format='multispec'
    iraf.kpnoslit.apall.interact='no'
    iraf.kpnoslit.apall.find='yes'
    iraf.kpnoslit.apall.recente='no'
    iraf.kpnoslit.apall.resize='no'
    iraf.kpnoslit.apall.edit='no'
    iraf.kpnoslit.apall.trace='no'
    iraf.kpnoslit.apall.fittrac='no'
    iraf.kpnoslit.apall.extrac='yes'
    iraf.kpnoslit.apall.extras='yes'
    #iraf.kpnoslit.apall.review='yes'
    #iraf.kpnoslit.apall.nfind= 2.0
    #iraf.kpnoslit.apall.t_func='legendre'
    #iraf.kpnoslit.apall.t_order= 4.0
    iraf.kpnoslit.apall.backgro='none'
    iraf.kpnoslit.apall.readnoi='HIERARCH ESO DET OUT1 RON'
    iraf.kpnoslit.apall.gain = 'HIERARCH ESO DET OUT1 GAIN'
    iraf.kpnoslit.apall.input = wave
    iraf.kpnoslit.apall.reference = ref_science
    iraf.kpnoslit.apall.output = 'w.ms.' + ref_science
    if os.path.isfile('w.ms.' + ref_science):
       subprocess.call('rm w.ms.' + ref_science,shell=True)
    iraf.kpnoslit.apall(mode='h')
    return 'w.ms.' + ref_science

def identify (output,REF):

    if REF=='':
       #print 'estas en el' + output

       iraf.kpnoslit.identify.images = output
       iraf.kpnoslit.identify.section= 'middle line'
       iraf.kpnoslit.identify.database = 'database'
       iraf.kpnoslit.identify.coordlist = 'linelists$fors.dat'
       iraf.kpnoslit.identify.function = 'spline3'
       iraf.kpnoslit.identify.order = '5'
       #iraf.kpnoslit.identify.function = 'spline3'
       #iraf.kpnoslit.identify.order = '5'
       #iraf.kpnoslit.identify.fwidth = 5
       #iraf.lpar(iraf.kpnoslit.identify)
       iraf.kpnoslit.identify(mode='h')

       #change this because header of the first one  is stupid


    else:
       #print 'estas en el ' + output

       iraf.kpnoslit.reidentify.images = output
       iraf.kpnoslit.reidentify.reference = REF
       iraf.kpnoslit.reidentify.section = 'middle line'
       iraf.kpnoslit.reidentify.database = 'database'
       iraf.kpnoslit.reidentify.answer = 'no'
       iraf.kpnoslit.reidentify.refit = 'no'
       iraf.kpnoslit.reidentify.verbose = 'yes'
       #iraf.kpnoslit.reidentify.shift = 'INDEF'
       #iraf.kpnoslit.reidentify.search = 'INDEF'
       iraf.kpnoslit.reidentify(mode='h')

def edit_header (input1,input2):
    abrir1 = open(input1,'r')
    abrir2 = open(input2,'r')
    for line in abrir1:
        linea = abrir2.readline().split('\n')[0]
        line = line.split('\n')[0]
        iraf.hedit.images = line
        iraf.hedit.fields = 'REFSPEC1'
        iraf.hedit.add  = 'yes'
        iraf.hedit.value = linea
        iraf.hedit.show = 'yes'
        iraf.hedit.ver = 'no'
        iraf.hedit(mode='h')

    abrir1.close()
    abrir2.close()

def wave_cal (images):
    iraf.dispcor.glo = 'no'
    iraf.dispcor.flux = 'yes'
    iraf.dispcor.input = '@'+ images
    iraf.dispcor.output = 'd//@' + images
    iraf.lpar(iraf.dispcor)
    iraf.dispcor(mode = 'h')

def sens(std):
    iraf.sensfunc.standard = std
    iraf.sensfunc.sensitiv = 'sens'
    iraf.sensfunc.interac= 'yes'
    iraf.sensfunc.ignorea= 'no'
    iraf.sensfunc.funct = 'spline3'
    iraf.sensfunc.order = 6
    iraf.sensfunc.observatory = 'observatory'
    iraf.sensfunc.extinc = 'onedstds$lasilla.dat'
    #iraf.lpar(iraf.sensfunc)
    iraf.sensfunc(mode='h')

def standards(images):
    iraf.observatory.command = 'set'
    iraf.observatory.obsid = 'esovlt'
    iraf.observatory(mode='h')
    outfile = open(images,'r')

    for line in outfile:
         dimage = line.split('\n')[0]

         hdulist = fits.open(dimage)
         exptime = hdulist[0].header['EXPTIME']
         airmass = hdulist[0].header['HIERARCH ESO TEL AIRM END']
         hdulist.close()
         iraf.standard.input = dimage
         iraf.standard.output = 'std'
         iraf.standard.extinction =  'onedstds$lasilla.dat'
         iraf.standard.caldir = 'onedstds$ctionewcal/'
         iraf.standard.observa = 'observatory'
         #iraf.lpar(iraf.standard)
         iraf.standard.airmass = airmass
         iraf.standard.exptime = exptime
         #iraf.standard.star_name = star
         iraf.standard(mode='h')

    outfile.close()

def inputs(sysargv):
    dictionary = {}
    path = {}
    if len(sysargv)!= 7:
        print 'USAGE: bias=bias std=std wave=wave flux=flux science=science flat=flat'
        sys.exit()
    for i in range(len(sysargv)-1):
        line = sysargv[i+1].split('=')
        key = line[0]
        aux = line[1].split('/')
        name = aux[-1]
        paths = '/'.join(aux[0:-1]) + '/'
        path[key] = paths
        dictionary[key] = name
    return dictionary, path

def BIAS(bias):
    #iraf.unlearn(iraf.zerocombine)
    path_bias = os.path.dirname(bias)
    bias = bias.replace(path_bias + '/','')

    iraf.chdir(path_bias)
    os.chdir(path_bias)
    iraf.zerocombine.input = '@' + bias
    if os.path.isfile('MasterBias.fits'):
       subprocess.call(['rm','MasterBias.fits'])
    iraf.zerocombine.output = 'MasterBias.fits'
    iraf.zerocombine.rdnoise = 'HIERARCH ESO DET OUT1 RON'
    iraf.zerocombine.gain = 'HIERARCH ESO DET OUT1 GAIN'
    iraf.zerocombine.ccdtype = ''
    #iraf.lpar(iraf.zerocombine)
    iraf.zerocombine(mode='h')
    #iraf.unlearn(iraf.zerocombine)

def FLAT(masterbias,flat): # FALTA NORMALIZAR EL MASTERFLAT

    complete_reduction(flat,masterbias,'')

    path_flat = os.path.dirname(flat)
    flat = flat.replace(path_flat + '/','')

    iraf.chdir(path_flat)
    os.chdir(path_flat)

    iraf.flatcombine.input = '@' + flat
    if os.path.isfile('MasterFlat.fits'):
       subprocess.call(['rm', 'MasterFlat.fits'])
    iraf.flatcombine.output = 'MasterFlat.fits'
    iraf.flatcombine.process = 'no'
    iraf.flatcombine.subsets = 'no'
    iraf.flatcombine.rdnoise = 'HIERARCH ESO DET OUT1 RON'
    iraf.flatcombine.gain = 'HIERARCH ESO DET OUT1 GAIN'
    iraf.flatcombine.ccdtype = ''
    iraf.flatcombine.statsec = '[800:900,5:105]'
    iraf.flatcombine.reject = 'minmax'
    #iraf.lpar(iraf.flatcombine)
    iraf.flatcombine(mode='h')
    #data = pyfits.getdata('MasterFlat.fits')
    #maxim = np.max(data)
    #iraf.imarith.operand1 = 'MasterFlat.fits'
    #iraf.imarith.operand2 = maxim
    #iraf.imarith.op = '/'
    #iraf.imarith.result = 'MasterFlat_norm.fits'
    #iraf.imarith(mode='h')

def complete_reduction (lista,masterbias,masterflat,cosmic='no'):
    full_bias = masterbias
    full_flat = masterflat
    full_lista = lista
    path_lista = os.path.dirname(full_lista)
    lista = full_lista.replace(path_lista + '/','')

    if masterbias != '':
       masterbias = full_bias.replace(os.path.dirname(full_bias) + '/','')
       if os.path.isfile(path_lista + '/' + masterbias) and os.path.dirname(full_bias) != path_lista:
           subprocess.call('rm ' + path_lista + '/' + masterbias,shell=True)
           subprocess.call('cp ' + full_bias  + ' ' + path_lista, shell=True)
       if not os.path.isfile(path_lista + '/' + masterbias):
           subprocess.call('cp ' + full_bias  + ' ' + path_lista, shell=True)
    if masterflat != '':
       masterflat = full_flat.replace(os.path.dirname(full_flat) + '/','')
       if os.path.isfile(path_lista + '/'+ masterflat) and os.path.dirname(full_flat) != path_lista:
          subprocess.call('rm ' + path_lista + '/' + masterflat,shell=True)
          subprocess.call('cp ' + full_flat + ' ' + path_lista,shell=True)
       if not os.path.isfile(path_lista + '/' + masterflat):
          subprocess.call('cp ' + full_flat  + ' ' + path_lista, shell=True)
    os.chdir(path_lista)
    iraf.chdir(path_lista)

    print 'Reduction ' + lista
    lis = open(lista,'r')
    for line in lis:
        line = line.split('\n')
        hdulist = fits.open(line[0])
        X = hdulist[0].header['HIERARCH ESO DET CHIP1 X'] #Start of X = 1
        Y = hdulist[0].header['HIERARCH ESO DET CHIP1 Y'] #Start of Y = 1
        NX = hdulist[0].header['HIERARCH ESO DET OUT1 NX'] #Valid in X
        NY = hdulist[0].header['HIERARCH ESO DET OUT1 NY'] #Valid in Y
        bias_sec = ''
        trim_sec= ''
        Y_MIN = 910  # You need to set this!
        Y_MAX = 1131 # You need to set this!
        X_TRIM = 200
        PRESCANX = hdulist[0].header['HIERARCH ESO DET OUT1 PRSCX'] #PRESCAN X
        PRESCANY = hdulist[0].header['HIERARCH ESO DET OUT1 PRSCY'] #PRESCAN Y
        OVERSCANX = hdulist[0].header['HIERARCH ESO DET OUT1 OVSCX'] #OVERSCAN X
        OVERSCANY = hdulist[0].header['HIERARCH ESO DET OUT1 OVSCY'] #OVERSCAN Y

        iraf.ccdproc.images = line[0]
        iraf.ccdproc.output = ''
        iraf.ccdproc.ccdtype = ''
        iraf.ccdproc.fixpix = 'no'
        iraf.ccdproc.flatcor = 'no'
        iraf.ccdproc.zerocor = 'no'
        iraf.ccdproc.darkcor = 'no'

        if str(PRESCANX) == '0' and str(PRESCANY) == '0' and str(OVERSCANX) == '0' and str(OVERSCANY) == '0':
           iraf.ccdproc.overscan = 'no'
           print "NO OVERSCAN"
           if NY > Y_MAX:
              iraf.ccdproc.trim = 'yes'
              print 'TRIMMIN'
              trim_sec = '[' + str(X_TRIM) + ':' + str(NX - OVERSCANX) + ',' + str(Y) + ':' + str(NY) + ']'
              iraf.ccdproc(mode='h')
           else:
              iraf.ccdproc.trim = 'no'
              print 'NO TRIMMIN'
        else:
            iraf.ccdproc.overscan = 'yes'
            iraf.ccdproc.trim = 'yes'
            if PRESCANX > 0 and OVERSCANX >= 0:
               if PRESCANY > 0 and OVERSCANY >= 0:
                  bias_sec = '[' + str(X) + ':' + str(PRESCANX) + ',' +  str(Y) + ':' + str(PRESCANY) + ']'
               elif PRESCANY >= 0 and OVERSCANY >0:
                    bias_sec = '[' + str(X) + ':' + str(PRESCANX) + ',' +  str(NY- OVERSCANY) + ':' + str(NY) + ']'
               else:
                    bias_sec = '[' + str(X) + ':' + str(PRESCANX) + ',*]'
            elif PRESCANX >= 0 and OVERSCANX > 0:
                 if PRESCANY > 0 and OVERSCANY >= 0:
                    bias_sec = '[' + str(NX-OVERSCANX) + ':' + str(NX) + ',' +  str(Y) + ':' + str(PRESCANY) + ']'
                 elif PRESCANY >= 0 and OVERSCANY >0:
                      bias_sec = '['+ str( NX- OVERSCANX) + ':' + str(NX) +  ',' +  str(NY- OVERSCANY) + ':' + str(NY) + ']'
                 else:
                      bias_sec = '[' + str(NX-OVERSCANX) + ':' + str(NX) + ',*]'
            else:
                 if PRESCANY > 0 and OVERSCANY >= 0:
                    bias_sec = '[*,' +  str(Y) + ':' + str(PRESCANY) + ']'
                 elif PRESCANY >= 0 and OVERSCANY >0:
                      bias_sec = '[*,' +  str(NY-OVERSCANY) + ':' + str(NY) +']'
            if NY > Y_MAX:
               trim_sec = '[' + str(X_TRIM) + ':' + str(NX - OVERSCANX) + ',' + str(Y_MIN) + ':' + str(Y_MAX) + ']'
            else:
               trim_sec = '[' + str(X_TRIM) + ':' + str(NX - OVERSCANX) + ',' + str(Y) + ':' + str(NY) + ']'

            print 'OVERSCAN & TRIMMIN \n %s %s' % (bias_sec, trim_sec)
            iraf.ccdproc.biassec = bias_sec
            iraf.ccdproc.trimsec = trim_sec
            iraf.ccdproc(mode='h')

        if masterbias != '' and masterflat != '':
           if cosmic == 'yes':
             c_out = single_cosmic(line[0])
             with open('c'+lista,'a') as outfile:
                  outfile.write(c_out + '\n')
           if not os.path.isfile('response.fits'):
              response(masterflat)
              apply_imarith(masterflat,'/','response.fits')
           apply_imarith(c_out,'/','response.fits')

        iraf.ccdproc.flatcor = 'yes'
        if masterflat == '':
           iraf.ccdproc.flatcor = 'no'
        if masterbias == '':
           iraf.ccdproc.zerocor = 'no'
        iraf.ccdproc.darkcor= 'no'
        iraf.ccdproc.zero = masterbias
        iraf.ccdproc.flat = masterflat
        #iraf.lpar(iraf.ccdproc)
        iraf.ccdproc(mode='h')
        hdulist.close()
    lis.close()
    if os.path.isfile(masterbias) and os.path.dirname(full_bias) != path_lista:
       subprocess.call(['rm',masterbias])
    if os.path.isfile(masterflat) and os.path.dirname(full_flat) != path_lista:
       subprocess.call(['rm',masterflat])

def cosmic_rays (List_image):
    full_path = List_image
    path = os.path.dirname(full_path)
    List_image = full_path.replace(path + '/','')
    os.chdir(path)
    iraf.chdir(path)

    image_list = open(List_image,'r')
    image_out = open('c'+List_image,'w')
    for line in image_list:
        line = line.split('\n')
        print 'Cosmic rays removal: ' + line[0]
        hdulist = fits.open(line[0])
        #print hdulist[0].header['HIERARCH ESO DET OUT1 RON']
        iraf.lacos_spec.input =  line[0]
        iraf.lacos_spec.output = 'c' + line[0]
        if os.path.isfile('c' + line[0]):
           os.system('rm c'+line[0])
        if os.path.isfile('mask_'+line[0]):
           os.system('rm mask_'+line[0])
        iraf.lacos_spec.outmask= 'mask_'+line[0]
        iraf.lacos_spec.readn= hdulist[0].header['HIERARCH ESO DET OUT1 RON']
        iraf.lacos_spec.gain = hdulist[0].header['HIERARCH ESO DET OUT1 GAIN']
        iraf.lacos_spec.verbose = 'no'
        iraf.lacos_spec(mode='h')
        os.system('rm mask_'+line[0])
        image_out.write('c'+line[0]+'\n')
        hdulist.close()
    image_list.close()
    image_out.close()

def dir_path (path):
    os.chdir(path)
    iraf.chdir(path)

def response(FLAT):
    iraf.longslit.response.calib = FLAT
    iraf.longslit.response.norma = FLAT
    iraf.longslit.response.response = 'response.fits'
    iraf.longslit.response.interac = 'yes'
    iraf.longslit.response.funct = 'spline3'
    iraf.longslit.response.order = '20'
    iraf.longslit.response(mode='h')

def apply_imarith (op1,op,op2):
    if op1.find('fits')!= -1:
       operand = op1
       result = ''
    else:
        operand = '@' + op1
        result = operand
    iraf.imarith.operand1 = operand
    iraf.imarith.op = op
    iraf.imarith.operand2 = op2
    iraf.imarith.result = result
    #print operand,op,op2
    iraf.imarith(mode='h')
    print 'Applying %s to %s' % (op2,op1)

def single_cosmic (image):
    hdulist = fits.open(image)
    iraf.lacos_spec.input =  image
    iraf.lacos_spec.output = 'c' + image
    if os.path.isfile('c' + image):
       os.system('rm c'+image)
    if os.path.isfile('mask_'+image):
       os.system('rm mask_'+image)
    iraf.lacos_spec.outmask= 'mask_'+image
    iraf.lacos_spec.readn= hdulist[0].header['HIERARCH ESO DET OUT1 RON']
    iraf.lacos_spec.gain = hdulist[0].header['HIERARCH ESO DET OUT1 GAIN']
    iraf.lacos_spec.verbose = 'no'
    iraf.lacos_spec(mode='h')
    os.system('rm mask_'+ image)
    hdulist.close()
    return 'c' + image

def wimage_ms_ref(waves,list_image,wave_path,ref_ID):
    image_ms = open(list_image,'r')
    w_image_ms = open('w.ms.' + list_image,'w')
    wave_list = open(waves,'r')
    for line in image_ms:
        im = line.split('\n')[0]
        hdulist = fits.open(im)
        angle_rot = hdulist[0].header['HIERARCH ESO INS RETA2 ROT']
        for linea in wave_list:
             wave = linea.split('\n')[0]
             hdu_aux = fits.open(wave)
             angle_rot_2 = hdu_aux[0].header['HIERARCH ESO INS RETA2 ROT']
             if angle_rot == angle_rot_2:
                out = wave_extraction(wave,im)
                w_image_ms.write(out + '\n')
                identify(out,ref_ID)
                w_image_ms.seek(0)
                break
    image_ms.close()
    w_image_ms.close()
    wave_list.close()
    return 'w.ms.' + list_image


if __name__ == '__main__':

  iraf.task(lacos_spec= "lacos_spec.cl")  #LACOSMIC
  abspath = os.path.abspath(__file__) #Definir main_directory
  dictionary, path = inputs(sys.argv)       #definir inputs
  #Setting variables and lacos_spec in directory
  wave = dictionary['wave']
  flat = dictionary['flat']
  std = dictionary['std']
  science = dictionary['science']
  flux = dictionary['flux']
  bias = dictionary['bias']
  #COSMIC REMOVAL
  c_wave = 'c' + wave
  c_flat = 'c' + flat
  c_science = 'c' + science
  c_flux = 'c' + flux
  #c_std = 'c' + std
  #AFTER APALL
  c_ms_science = 'ms.'+ c_science
  c_ms_flux = 'ms.' + c_flux
  #c_ms_std = 'ms.' + c_std
  #WITH WAVE CAL
  w_ms_science = 'w.' + c_ms_science
  #w_ms_std = 'w.' + c_ms_std
  #w_ms_flux = 'w.' + c_ms_flux
  #Copying lacos_spec

  subprocess.call('echo %s %s %s %s %s | xargs -n 1 cp lacos_spec.cl' % (path['flat'],path['flux'],path['science'],path['wave'],path['std']),shell=True)

  reduce_data = raw_input('Reduce data?(y|N): ' )
  if reduce_data == 'y':
       #LACOSMIC
     cosmic = 'no'
     remove_cosmics = raw_input('Remove_cosmics now?(y|N): ')
     if remove_cosmics =='y':
        cosmic = 'yes'

     BIAS(path['bias'] + bias)
     FLAT(path['bias'] +"MasterBias.fits",path['flat'] + flat)

     complete_reduction(path['science'] + science, path['bias'] +"MasterBias.fits", path['flat'] + "MasterFlat.fits",cosmic=cosmic)
     complete_reduction(path['wave'] + wave,  path['bias'] + "MasterBias.fits", path['flat'] + "MasterFlat.fits",cosmic=cosmic)
     complete_reduction(path['flux'] + flux, path['bias'] + "MasterBias.fits",  path['flat'] + "MasterFlat.fits",cosmic=cosmic)

     if remove_cosmics != 'y':
        cosmic_rays(path['science'] + science) #c_science
        cosmic_rays(path['wave'] + wave) #c_wave
        cosmic_rays(path['flux'] + flux) #c_flux

 #SPECTRA for SCIENCE
  dir_path(path['science'])
  spectra(c_science) # output = c_ms_science

  dir_path(os.path.dirname(abspath))
  image= Image.open('fors.jpeg').show()
  dir_path(path['science'])

  ref_sci_ID = wave_spectra( path['wave'] + c_wave, c_science, path['wave']) # output = w_ms_science, ref = image with ID
  edit_header(c_ms_science,w_ms_science)
  wave_cal(c_ms_science) # output = dc_ms_science
  subprocess.call('mv d*.fits ' + os.path.dirname(abspath),shell=True)

  #subprocess.call('rm -R database',shell=True)

#SPECTRA FOR FLUX
  dir_path(path['flux'])
  spectra(c_flux)
#####wave_spectra(path['wave'] + c_wave, c_flux, path['wave'])
  #edit_header(c_ms_flux, w_ms_science)
  w_ms_flux = wimage_ms_ref(path['wave'] + c_wave, c_flux, path['wave'],ref_sci_ID)
  edit_header(c_ms_flux,w_ms_flux)
  wave_cal(c_ms_flux)
  subprocess.call('ls d*.fits > dflux',shell=True)
  #LOOK THIS!!!!!!!!!!!!!!!!
  standards('dflux')
  sens('std')
  subprocess.call('mv sens*.fits ' + os.path.dirname(abspath), shell=True)
  #subprocess.call('cp -R database ' + os.path.dirname(abspath), shell=True)
  #subprocess.call('rm -R database c*.fits w*.fits %s %s %s %s %s %s %s dflux d*.fits std' % (c_flux,c_wave,c_science,c_ms_science,w_ms_flux,w_ms_science,c_ms_flux),shell=True)
  dir_path(os.path.dirname(abspath))

  print '\nThis is the End: Output'
  subprocess.call('ls *.fits*',shell=True)

## generating input for .pro

  with open('vltpol','w') as outfile:
       dfits = glob.glob('d*.fits')
       sen = glob.glob('sens*')
       new_dfits = np.zeros(4)
       for j in range(len(dfits)):
           hdu_list = fits.open(dfits[j])
           ang = hdu_list[0].header['HIERARCH ESO INS RETA2 ROT']
           if ang == '0.0':
              new_dfits[0] = dfits[j]
           elif ang == '45.0':
              new_dfits[1] = dfits[j]
           elif ang == '22.5':
               new_dfits[2] = dfits[j]
           else:
               new_dfits[3] = dfits[j]

       for k in range(len(new_dfits)):
           outfile.write('%s\n' % new_dfits[k])

       for i in range(len(sen)):
           outfile.write('%s\n' % sen[i])
