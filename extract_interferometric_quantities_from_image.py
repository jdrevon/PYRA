# -*- coding: utf-8 -*-
"""
Created on Thu Nov  7 00:16:05 2024

@author: jdrevon
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Oct 23 23:07:15 2024

@author: jdrevon
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 08:19:27 2024

@author: jdrevon
"""

from DATA_reading import OIFITS_READING_concatenate
import matplotlib.pyplot as plt
import numpy as np
from numpy.fft import ifftshift, fftshift, fftfreq, fft2
from astropy.io import fits
from scipy.interpolate import RegularGridInterpolator, interpn
from astropy import units as units
import math
from scipy.ndimage import shift

def orientation_image(file_path):

    # Ouvrir le fichier FITS
    with fits.open(file_path) as hdul:
        image = hdul[0].data  # Charger les données de l'image
        header = hdul[0].header
        
        # Lire les valeurs de CDELT1 et CDELT2
        cdelt1 = header.get('CDELT1', 1.0)  # Valeur par défaut 1.0 si absente
        cdelt2 = header.get('CDELT2', 1.0)
        
        # Appliquer les flips sans modifier le header
        if cdelt1 < 0:  # Flip horizontal si CDELT1 est négatif
            image = np.flip(image, axis=1)
        if cdelt2 < 0:  # Flip vertical si CDELT2 est négatif
            image = np.flip(image, axis=0)
    
    return image


def reading_info_OIFITS(DATA_OBS, path_image):

    # Read data
    
    wavel            = DATA_OBS['VIS2']['WAVEL'][np.invert(DATA_OBS['VIS2']['FLAG'])]
    U                = DATA_OBS['VIS2']['U'][np.invert(DATA_OBS['VIS2']['FLAG'])]
    V                = DATA_OBS['VIS2']['V'][np.invert(DATA_OBS['VIS2']['FLAG'])]
    V2_MATISSE       = DATA_OBS['VIS2']['VIS2'][np.invert(DATA_OBS['VIS2']['FLAG'])]
    V2_MATISSE_ERR   = DATA_OBS['VIS2']['VIS2_ERR'][np.invert(DATA_OBS['VIS2']['FLAG'])]
    q_u              = U / wavel
    q_v              = V / wavel
    
    U1, U2           = DATA_OBS['T3']['U1'][np.invert(DATA_OBS['T3']['FLAG'])],DATA_OBS['T3']['U2'][np.invert(DATA_OBS['T3']['FLAG'])]
    U3               = U1+U2
    
    V1,V2            = DATA_OBS['T3']['V1'][np.invert(DATA_OBS['T3']['FLAG'])], DATA_OBS['T3']['V2'][np.invert(DATA_OBS['T3']['FLAG'])]
    V3               = V1+V2
    
    T3               = DATA_OBS['T3']['T3'][np.invert(DATA_OBS['T3']['FLAG'])]
    T3_ERR           = DATA_OBS['T3']['T3_ERR'][np.invert(DATA_OBS['T3']['FLAG'])]
    WL               = DATA_OBS['T3']['WAVEL'][np.invert(DATA_OBS['T3']['FLAG'])]
    
    q_u1, q_u2, q_u3 = U1/WL, U2/WL, U3/WL 
    q_v1, q_v2, q_v3 = V1/WL, V2/WL, V3/WL

    image, header = read_fits_image(path_image) #Read the image
    image = orientation_image(path_image)
    return q_u, q_v, q_u1, q_u2, q_u3, q_v1, q_v2, q_v3, V2_MATISSE, V2_MATISSE_ERR, T3, image, header, T3_ERR

def mas_to_rad(values):
    """Convert milliarcseconds to radians."""
    return values * (np.pi / (180 * 3600 * 1000))

def rad_to_mas(values):
    """Convert radians to milliarcseconds."""
    return values / (np.pi / (180 * 3600 * 1000))

def create_coordinate_arrays(image_shape, x_center, y_center, x_scale, y_scale):
    x_size, y_size = image_shape
    # x_image = (np.arange(x_size) - x_center) * x_scale
    # y_image = (np.arange(y_size) - y_center) * y_scale
    x = np.linspace(-0.5, 0.5, x_size)
    x_image = x*(x_size*x_scale)
    y_image = x_image

    # return np.linspace(-x_size*x_scale/2,x_size*x_scale/2,x_size), np.linspace(-y_size*y_scale/2,y_size*y_scale/2,y_size)
    return x_image, y_image



def pad_image(image, padding):
    """Pads an image with additional zeros for Fourier transform."""
    dimy, dimx = image.shape
    padx = (dimx * padding - dimx) // 2
    pady = (dimy * padding - dimy) // 2
    return np.pad(image, ((pady, pady), (padx, padx)), 'constant', constant_values=0)

def spatial_to_frequency_rad(image, pixelsize):
    """Convert spatial lengths in rad to spatial frequencies in rad^-1 for a 2D grid."""
    # delta_x = x_image_rad[1] - x_image_rad[0]
    # delta_y = y_image_rad[1] - y_image_rad[0]

    # freq_x = fftshift(fftfreq(len(x_image_rad), d=delta_x)) 
    # freq_y = fftshift(fftfreq(len(y_image_rad), d=delta_y))

    freq  = fftshift(fftfreq(image.shape[0], pixelsize))

    return freq

def center_on_photocenter(image):
    """
    Centers an image based on its photocenter (center of mass of intensity).
    :param image: 2D numpy array of the image.
    :return: Centered image.
    """
    y, x = np.indices(image.shape)
    total_intensity = np.sum(image)
    center_y = np.sum(y * image) / total_intensity
    center_x = np.sum(x * image) / total_intensity
    
    shift_y = image.shape[0] / 2 - center_y
    shift_x = image.shape[1] / 2 - center_x
    
    return shift(image, shift=(shift_y, shift_x), mode='constant', cval=0)



def image_to_FFT(image):
    """Compute the 2D Fourier Transform of the image."""
    # return   np.fft.ifftshift(np.fft.fft2(np.fft.fftshift(image, axes=[-2, -1]), axes=[-2, -1]), axes=[-2, -1])
    return   ifftshift(fft2(fftshift(image)))
    # return np.abs(np.fft.fftshift(np.fft.ifft2(np.fft.fftshift(image, axes=[-2, -1]), axes=[-2, -1]), axes=[-2, -1]))

def read_fits_image(filename):
    """Read the FITS file."""
    
    hdul = fits.open(filename)
    header = hdul[0].header
    unitx= header['CDELT1']
    unity= header['CDELT2']
    
    if np.abs(unitx) == np.abs(unity):
        
        if np.logical_and(unitx > 0, unity>0):
            
            hdul = fits.open(filename)
            image = hdul[0].data
            header = hdul[0].header
            hdul.close()

        elif np.logical_and(unitx < 0, unity>0):

            hdul = fits.open(filename)
            image = np.flip(hdul[0].data,axis=1)
            # image = hdul[0].data
            header = hdul[0].header
            hdul.close()            

        elif np.logical_and(unitx > 0, unity<0):

            hdul = fits.open(filename)
            image = np.flip(hdul[0].data,axis=0)
            header = hdul[0].header
            hdul.close()            

    else:
        print('Not a regular grid')
        exit() 
    
    hdul.close()
    return image, header

def extract_header_info(header):
    """Extract necessary information from the FITS header and convert them in radians."""
    
    unit= header['CUNIT1']
    
    if unit == 'mas':
        
        if header['CDELT1'] < 0:
        
            x_scale  = -header['CDELT1']*units.mas.to(units.rad)
            y_scale  = header['CDELT2']*units.mas.to(units.rad)

        elif header['CDELT2'] < 0:
        
            x_scale  = header['CDELT1']*units.mas.to(units.rad)
            y_scale  = -header['CDELT2']*units.mas.to(units.rad)

        else:
            
            x_scale  = header['CDELT1']*units.mas.to(units.rad)
            y_scale  = header['CDELT2']*units.mas.to(units.rad)


    elif unit == 'rad':

        if header['CDELT1'] < 0:
            x_scale  = -header['CDELT1']
            y_scale  = header['CDELT2']

        elif header['CDELT2'] < 0:
            x_scale  = header['CDELT1']
            y_scale  = -header['CDELT2']

        else:            
            x_scale  = header['CDELT1']
            y_scale  = header['CDELT2']
        
        
    x_center = header['CRPIX1'] #-1  # FITS headers are 1-indexed
    y_center = header['CRPIX2'] #- 1  # FITS headers are 1-indexed
    
    return x_center, y_center, x_scale, y_scale


def IMAGE_to_V2(DATA_OBS, path_image, padding):
    """Convert an image to visibility squared (V2) data and extract the CP values"""

    # #Center the original image on the brightest pixel
    # image_centered = center_image_on_brightest_pixel(image)
    
    q_u_interp, q_v_interp, q_u1, q_u2, q_u3, q_v1, q_v2, q_v3, V2_MATISSE, V2_MATISSE_ERR, CP_MATISSE, image, header, CP_MATISSE_ERR = reading_info_OIFITS(DATA_OBS, path_image)
    
    def image_is_even(image):
        rows, cols = image.shape
        return rows % 2 == 0 and cols % 2 == 0
    
    def pad_image_to_even(image):
        rows, cols = image.shape
        pad_y = 0 if rows % 2 == 0 else 1
        pad_x = 0 if cols % 2 == 0 else 1
        
        padded_image = np.pad(image, ((0, pad_y), (0, pad_x)), mode='constant', constant_values=0)
        
        return padded_image

    if image_is_even(image):  # Exemple de vérification
        print("Image paire")
    else:
        print("Image impaire")
        image = pad_image_to_even(image)

    x_center, y_center, x_scale, y_scale = extract_header_info(header) #read header info
    x_image, y_image = create_coordinate_arrays(image.shape, x_center, y_center, x_scale, y_scale) #create the x and y-axis list

    
    # plt.figure()
    # plt.imshow(image/np.max(image), extent=(min(rad_to_mas(x_image)), max(rad_to_mas(x_image)), min(rad_to_mas(y_image)), max(rad_to_mas(y_image))), origin='lower')
    # ax=plt.gca()
    # ax.set_xlim(max(rad_to_mas(x_image)), min(rad_to_mas(x_image)))
    # plt.colorbar(label="Normalized Intensity")

    if padding > 0:
        padded_image = pad_image(image, padding)
        FFT_image = image_to_FFT(padded_image)
        q_image = spatial_to_frequency_rad(padded_image, x_scale)
        
    else:
        FFT_image = image_to_FFT(image)
        q_image = spatial_to_frequency_rad(image, x_scale)
 

    # FFT_final = np.abs(FFT_image/ np.amax(FFT_image))**2
    # FFT_final = np.abs(FFT_image/ np.amax(FFT_image))**2
    # FFT_final = np.abs(dft_result**2 / np.amax(dft_result)**2)
    
    
    # plt.figure()
    # plt.imshow(FFT_final, extent=(min(q_x), max(q_x), min(q_y), max(q_y)), origin='lower', norm=colors.PowerNorm(gamma=0.5))
    
    # plt.xlabel("Spatial frequency (cycles/rad)")
    # plt.ylabel("Spatial frequency (cycles/rad)")
    # plt.title("Fourier Transform of the Image")
    # plt.colorbar(label="Normalized Intensity")

    FFT_real = np.real(FFT_image)
    FFT_imag = np.imag(FFT_image)
    
    interpolator_real = RegularGridInterpolator((q_image, q_image), FFT_real, method='linear', bounds_error=False, fill_value=None)
    interpolator_imag = RegularGridInterpolator((q_image, q_image), FFT_imag, method='linear', bounds_error=False, fill_value=None)    

    def MIRA_complex(y, x):
          return interpolator_real((y, x)) + 1j*interpolator_imag((y, x))

    # interpolator = RegularGridInterpolator((q_image, q_image), FFT_final, method='linear', bounds_error=False, fill_value=None)
    V2_MiRA     = np.abs([MIRA_complex(q_v_interp[b], q_u_interp[b])/MIRA_complex(0, 0) for b in range(len(q_u_interp))])**2

    # FFT_real = np.real(FFT_image)
    # FFT_imag = np.imag(FFT_image)

    # interpolator_real = RegularGridInterpolator((q_image, q_image), FFT_real, method='linear', bounds_error=False, fill_value=None)
    # interpolator_imag = RegularGridInterpolator((q_image, q_image), FFT_imag, method='linear', bounds_error=False, fill_value=None)

    # def MIRA_complex(y, x):
    #       return interpolator_real((y, x)) + 1j*interpolator_imag((y, x))

    q1 = np.sqrt(q_u1**2+q_v1**2)
    q2 = np.sqrt(q_u2**2+q_v2**2)
    q3 = np.sqrt(q_u3**2+q_v3**2)

    Vis_B1 = np.array([MIRA_complex(q_v1[b], q_u1[b]) for b in range(len(q_u1))])
    Vis_B2 = np.array([MIRA_complex(q_v2[b], q_u2[b]) for b in range(len(q_u2))])
    Vis_B3 = np.array([MIRA_complex(q_v3[b], q_u3[b]) for b in range(len(q_u3))])
        
    W = Vis_B1*Vis_B2*np.conj(Vis_B3)/MIRA_complex(0, 0)**3
            
    CP_MiRA = [np.angle(W[b], deg=True) for b in range(len(W))]
    # CP_MiRA = np.array(CP)/np.pi*180
    
    q_CP = np.amax ([q1,q2,q3], axis=0)

    q = np.sqrt(q_u_interp**2+q_v_interp**2)
        
    return q, V2_MATISSE, V2_MATISSE_ERR, V2_MiRA, q_CP, CP_MiRA, CP_MATISSE

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
