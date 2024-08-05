import numpy as np
from astropy.io import fits
from astropy.wcs import WCS


def gaussian(x, mu, sigma):
    return np.exp(-np.power(x - mu, 2.) / (2. * np.power(sigma, 2.))) # / (sigma * np.sqrt(2. * np.pi))



def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def read_config(config, availgrism='both'): 
    configfile = open(config, 'r')  
    config_pars = {} 
    for line in configfile:
        if ( (line[0] != '#') & (len(line)>1)): 
            name = line.split()[0] 
            
            if name == 'node_wave': 
                tmpnode = [float(s) for s in line.split()[1::]]
                if availgrism.lower() == 'g102':
                    nodelam = [x for x in tmpnode if x <= config_pars['transition_wave']]
                elif availgrism.lower() == 'g141':
                    nodelam = [x for x in tmpnode if x > config_pars['transition_wave']]
                else:
                    nodelam = tmpnode
                config_pars.setdefault('node_wave', []) 
                for l in nodelam:  
                    config_pars['node_wave'].append(l)
            elif name == 'mask_region1':
                masklam = [float(s) for s in line.split()[1::]]
                config_pars.setdefault('mask_region1', []) 
                config_pars['mask_region1'].append(masklam[0]) 
                config_pars['mask_region1'].append(masklam[1]) 
            elif name == 'mask_region2': 
                masklam = [float(s) for s in line.split()[1::]]
                config_pars.setdefault('mask_region2', []) 
                config_pars['mask_region2'].append(masklam[0]) 
                config_pars['mask_region2'].append(masklam[1])
            elif name == 'mask_region3': 
                masklam = [float(s) for s in line.split()[1::]]
                config_pars.setdefault('mask_region3', []) 
                config_pars['mask_region3'].append(masklam[0]) 
                config_pars['mask_region3'].append(masklam[1])  
            else:  
                val  = line.split()[1] 
                if is_number(val) == True:  val = float(val)
                if val == 'True': val = True
                if val == 'False': val = False 
                config_pars[line.split()[0]] = val 
    configfile.close()
    return config_pars    



# Creates .reg files for direct and dispersed images
# Written by Mason Huberty and adapted by KVN
def create_regions(parno, path_to_wisp_data):
    
    #Direct image
    hdul = fits.open('Par28_speccat.fits')
    cat=hdul[1].data
    f = open("Par28v2direct.reg",'a')
    for i in range(len(cat['ra'])):
        f.write("WCS;circle("+str(cat['ra'][i])+','+str(cat['dec'][i])+',0.5") # color=green text={'+str(cat['id'][i])+' z='+str(round(cat['redshift'][i],3))+'} font="times 10 bold italic" textangle=30\n')
    f.close()

    #This and subsequent are for the first order beams. Offsets taken from the config files
    f = open("Par28F115r_grism.reg",'a')
    w = WCS('Par28_228_f115w-gr150r_drz_sci.fits')
    for i in range(len(cat['ra'])):
        ra, dec = (cat['ra'][i], cat['dec'][i])
        x, y = w.all_world2pix(ra, dec, 1)
        f.write("box("+str(x-0.6548681566074263)+','+str(y-33.73739138173772)+',10.0,93.54,0.0) # color=red text={'+str(cat['id'][i])+'} font="times 10 bold italic" textangle=30\n')
    f.close()
    
    
    f = open("Par28F115c_grism.reg",'a')
    w = WCS('Par28_228_f115w-gr150c_drz_sci.fits')
    for i in range(len(cat['ra'])):
        ra, dec = (cat['ra'][i], cat['dec'][i])
        x, y = w.all_world2pix(ra, dec, 1)
        f.write("box("+str(x-31.91107156101387)+','+str(y-1.3922939626209256)+',97.28751330105166,10.0,0.0) # color=red text={'+str(cat['id'][i])+'} font="times 10 bold italic" textangle=30\n')
    f.close()


    f = open("Par28F150r_grism.reg",'a')
    w = WCS('Par28_228_f150w-gr150r_drz_sci.fits')
    for i in range(len(cat['ra'])):
        ra, dec = (cat['ra'][i], cat['dec'][i])
        x, y = w.all_world2pix(ra, dec, 1)   
        f.write("box("+str(x-0.6548681566074263)+','+str(y-106.79254657227568)+',10.0,93.54,0.0) # color=red text={'+str(cat['id'][i])+'} font="times 10 bold italic" textangle=30\n')
    f.close()
    
    
    f = open("Par28F150c_grism.reg",'a')
    w = WCS('Par28_228_f150w-gr150c_drz_sci.fits')
    for i in range(len(cat['ra'])): 
        ra, dec = (cat['ra'][i], cat['dec'][i])
        x, y = w.all_world2pix(ra, dec, 1)
        f.write("box("+str(x-96.44444)+','+str(y-0.6548681566074263)+',93.54,10.0,0.0) # color=red text={'+str(cat['id'][i])+'} font="times 10 bold italic" textangle=30\n')
    f.close()
    
    
    f = open("Par28F200r_grism.reg",'a')
    w = WCS('Par28_228_f200w-gr150r_drz_sci.fits')
    for i in range(len(cat['ra'])):
        ra, dec = (cat['ra'][i], cat['dec'][i])
        x, y = w.all_world2pix(ra, dec, 1)
        f.write("box("+str(x-0.6548681566074263)+','+str(y-204.8370874255101)+',10.0,131.78,0.0) # color=red text={'+str(cat['id'][i])+'} font="times 10 bold italic" textangle=30\n')
    f.close()
    
    
    f = open("Par28F200c_grism.reg",'a')
    w = WCS('Par28_228_f200w-gr150c_drz_sci.fits')
    for i in range(len(cat['ra'])):
        ra, dec = (cat['ra'][i], cat['dec'][i])
        x, y = w.all_world2pix(ra, dec, 1)
        f.write("box("+str(x-200.9228)+','+str(y-0.6548681566074263)+',127.806,10.0,0.0) # color=red text={'+str(cat['id'][i])+'} font="times 10 bold italic" textangle=30\n')
    f.close()
    
    
