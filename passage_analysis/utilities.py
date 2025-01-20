import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import glob
import os
from astropy.table import Table
import pandas as pd

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
            elif name == 'mask_region4': 
                masklam = [float(s) for s in line.split()[1::]]
                config_pars.setdefault('mask_region4', []) 
                config_pars['mask_region4'].append(masklam[0]) 
                config_pars['mask_region4'].append(masklam[1])  
            elif name == 'mask_region5': 
                masklam = [float(s) for s in line.split()[1::]]
                config_pars.setdefault('mask_region5', []) 
                config_pars['mask_region5'].append(masklam[0]) 
                config_pars['mask_region5'].append(masklam[1])  
            elif name == 'mask_region6': 
                masklam = [float(s) for s in line.split()[1::]]
                config_pars.setdefault('mask_region6', []) 
                config_pars['mask_region6'].append(masklam[0]) 
                config_pars['mask_region6'].append(masklam[1])  
            elif name == 'mask_region7': 
                masklam = [float(s) for s in line.split()[1::]]
                config_pars.setdefault('mask_region7', []) 
                config_pars['mask_region7'].append(masklam[0]) 
                config_pars['mask_region7'].append(masklam[1])  
            elif name == 'mask_region8': 
                masklam = [float(s) for s in line.split()[1::]]
                config_pars.setdefault('mask_region8', []) 
                config_pars['mask_region8'].append(masklam[0]) 
                config_pars['mask_region8'].append(masklam[1])  
            else:  
                val  = line.split()[1] 
                if is_number(val) == True:  val = float(val)
                if val == 'True': val = True
                if val == 'False': val = False 
                config_pars[line.split()[0]] = val 
    configfile.close()
    return config_pars    



# Creates .reg files for direct and dispersed images
# Written by Mason Huberty and adapted by KVN to work as part of the line finding software

def write_obj_region(parno, path_to_data, catalog, regfile_name, xoffset, yoffset, w, b_width, b_length): 
    file = open(path_to_data + "/Par" + str(parno) + "/DATA/" + "Par" + str(parno) + regfile_name, 'a')
    for i in range(len(catalog)):
        ra, dec = catalog['ra'][i], catalog['dec'][i]
        x, y = w.all_world2pix(ra, dec, 1)
        file.write("box("+str(x-xoffset)+','+str(y-yoffset)+','+str(b_width)+','+str(b_length)+',0.0) # color=red text={'+str(catalog['id'][i])+'} font="times 10 bold italic" textangle=30\n')
    file.close



# NOTE - The x & y offset values are specific for NIRISS. 
#        These will need to be updated for NIRCam data!
def create_regions(parno, path_to_data):

    spec_cat = glob.glob(path_to_data + "/Par" + str(parno) + "/DATA/DIRECT_GRISM/Par*spec*.fits")
    hdul = fits.open(spec_cat[0])
    cat=hdul[1].data

    # Direct image region files
    f = open(path_to_data + "/Par" + str(parno) + "/DATA/" + "Par" + str(parno) + 'regions_phot.reg','a')
    for i in range(len(cat)):
        f.write("WCS;circle("+str(cat['ra'][i])+','+str(cat['dec'][i])+',0.5") # color=green text={'+str(cat['id'][i])+' z='+str(round(cat['redshift'][i],3))+'} font="times 10 bold italic" textangle=30\n')
    f.close()

    #This and subsequent are for the first order beams. Offsets taken from the config files
    # KVN :  Adding code to first check if the paths exist. Only create region file if filter/orientation is available 
    f115grism_R = glob.glob(path_to_data + "/Par" + str(parno) + "/DATA/*f115w*gr150r_drz_sci.fits")
    if len(f115grism_R) != 0:
        write_obj_region(parno, path_to_data, cat, "F115r_grism.reg", 0.6548681566074263, 33.73739138173772, w = WCS(f115grism_R[0]), 
                        b_width = 10.0, b_length = 93.54)
    
    f115grism_C = glob.glob(path_to_data + "/Par" + str(parno) + "/DATA/*f115w*gr150c_drz_sci.fits") # Added KVN 19-Aug-2024
    if len(f115grism_C) != 0:
        write_obj_region(parno, path_to_data, cat, "F115c_grism.reg", 31.91107156101387, 1.3922939626209256, w = WCS(f115grism_C[0]), 
                        b_width = 97.28751330105166, b_length = 10.0)

    f150grism_R = glob.glob(path_to_data + "/Par" + str(parno) + "/DATA/*f150w*gr150r_drz_sci.fits") # Added KVN 19-Aug-2024
    if len(f150grism_R) != 0:
        write_obj_region(parno, path_to_data, cat, "F150r_grism.reg", 0.6548681566074263, 106.79254657227568, w = WCS(f150grism_R[0]), 
                        b_width = 10.0, b_length = 93.54)
    
    f150grism_C = glob.glob(path_to_data + "/Par" + str(parno) + "/DATA/*f150w*gr150c_drz_sci.fits") # Added KVN 19-Aug-2024
    if len(f150grism_C) != 0:
        write_obj_region(parno, path_to_data, cat, "F150c_grism.reg", 96.44444, 0.6548681566074263, w = WCS(f150grism_C[0]), 
                        b_width = 93.54, b_length = 10.0)
    
    f200grism_R = glob.glob(path_to_data + "/Par" + str(parno) + "/DATA/*f200w*gr150r_drz_sci.fits") # Added KVN 19-Aug-2024
    if len(f200grism_R) != 0:    
        write_obj_region(parno, path_to_data, cat, "F200r_grism.reg", 0.6548681566074263, 204.8370874255101, w = WCS(f200grism_R[0]), 
                        b_width = 10.0, b_length = 131.78)        
    
    f200grism_C = glob.glob(path_to_data + "/Par" + str(parno) + "/DATA/*f200w*gr150c_drz_sci.fits") # Added KVN 19-Aug-2024
    if len(f200grism_C) != 0:
        write_obj_region(parno, path_to_data, cat, "F200c_grism.reg", 200.9228, 0.6548681566074263, w = WCS(f200grism_C[0]), 
                        b_width = 127.806, b_length = 10.0)        
    

def add_header_keyword(parno, path_to_data):

    main_directory = os.path.join(path_to_data, f"Par{parno:s}")
    files = glob.glob(os.path.join(main_directory, 'spec2D/*.fits'))
    test_file = fits.open(files[0])
    header = test_file[2].header

    if header['RADESYS'] != 'ICRS':
        print('Headers are not updated. Updating all now.')
        for j in files:
            # Open the file header for viewing and load the header
            hdulist = fits.open(j)
            for i in range(len(hdulist)):
                header = hdulist[i].header
                try: header['RADESYS'] = 'ICRS'
                except: print('no RADESYS in header')

            hdulist.writeto(j, overwrite='True')

def make_spectra_dat_files(parno, path_to_data):

    main_directory = os.path.join(path_to_data, f"Par{parno:s}")
    spec1D_directory = os.path.join(main_directory, "spec1D")

    os.system(f"mkdir -p {os.path.join(main_directory, 'Spectra'):s}")

    files = sorted(glob.glob(os.path.join(spec1D_directory, '*1D.fits')))

    # Check if converted files already exist. If they do not, carry on with the conversion.
    # Otherwise, this step can be skipped.
    # Here, I just check if there are more converted files than objects (there can be 1-3 per object)
    # depending on how many filters are available for each object

    print('\nThere are ' + str(len(files)) + ' PASSAGE files to convert.\n')
    for f in files:

        fff = fits.open(f)
        print(f)

        for ext in range(1, len(fff)):

            tb = Table(fff[ext].data).to_pandas()
            t_out = pd.DataFrame(fff[ext].data)

            ### !!! IMPORTANT !!!
            ### VM: This is temporary while we fix this in the pipeline
            ### The R/C spectra have already been treated for the following operations
            if "EXTVER" not in fff[ext].header:

                t_out['wave'] = tb['wave']
                t_out['flux'] = tb['flux']/tb['flat']
                t_out['error'] = tb['err']/tb['flat']
                t_out['contam'] = tb['contam']/tb['flat']
                t_out['zeroth'] = np.zeros(len(tb['wave'])).astype('int')

            t_out = Table.from_pandas(t_out)
            ### VM: Explicitly fix nans to 0s
            for col in t_out.columns:
                t_out[col][np.isnan(t_out[col])] = 0
            t_out = t_out.filled(0.0) # Replace nans with zeros

            # Spectra dispersed beyond the chip have zero fluxes that must be replaced to prevent crashes in fitting.
            t_out['flux'][np.where(t_out['flux'] == 0.0)] = np.median(t_out['flux'][np.where(t_out['flux'] != 0.0)])
            t_out['error'][np.where(t_out['error'] == 0.0)]=np.median(t_out['error'][np.where(t_out['error'] != 0.0)])

            for filt in ["115", "150", "200"]:

                if fff[ext].header['EXTNAME'] == f'F{filt:s}W' and "EXTVER" not in fff[ext].header:
                    t_out.write(os.path.join(main_directory, "Spectra",
                                os.path.basename(f).replace('1D.fits', f'G{filt:s}_1D.dat')),
                                    format='ascii.fixed_width_two_line', overwrite=True)

                elif fff[ext].header['EXTNAME'] == f'F{filt:s}W' and "EXTVER" in fff[ext].header:
                    orient = fff[ext].header["FILTER"][-1]
                    t_out.write(os.path.join(main_directory, "Spectra",
                                os.path.basename(f).replace('1D.fits', f'G{filt:s}_1D_{orient:s}.dat')),
                            format='ascii.fixed_width_two_line', overwrite=True)

