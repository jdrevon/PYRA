# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 21:45:17 2024

@author: jdrevon
"""

import numpy as np
from astropy.io import fits 


def trouver_position_sous_array(liste_arrays, sous_array):
    for i, array in enumerate(liste_arrays):
        # Convertir les arrays en listes pour pouvoir utiliser la méthode sorted
        array_liste = list(array)
        sous_array_liste = list(sous_array)

        # Vérifier si les listes triées sont équivalentes
        if sorted(array_liste) == sorted(sous_array_liste):
            # Retourner la position du sous-array trouvé
            return i

    # Retourner -1 si le sous-array n'est pas trouvé
    return -1

def remove_spaces_in_string_array(string_array):
    # Use list comprehension to remove spaces from each string in the array
    result_array = np.array([string.replace(' ', '') for string in string_array])
    
    return result_array

def OIFITS_READING_concatenate(filename):
    
    OIFITS_TOT = np.zeros(1, dtype=object)
        
    with fits.open(filename, memmap=False) as fichier:

        name_HDU = np.array([fichier[t].name for t in range(len(fichier))])
        
        index_vis   = np.where(name_HDU=='OI_VIS2')[0]
        index_t3    = np.where(name_HDU=='OI_T3')[0]
        index_wvl   = np.where(name_HDU=='OI_WAVELENGTH')[0]
        index_flux = np.where(name_HDU=='OI_FLUX')[0]
        index_array = np.where(name_HDU=='OI_ARRAY')[0]

        insame_wvl  = np.array([fichier[t].header['INSNAME'] for t in index_wvl])
        
        dic = {}
        dic['NAME']                    = filename         

               
        vis2     = []
        vis2_err = []
        u_coord = []
        v_coord = []
        vis2_flag = []
        wavel     = []
        flag = []

        t3     = []
        t3_err = []
        t3_flag = []
        u1_coord = []
        v1_coord = []
        u2_coord = []
        v2_coord = []
        u3_coord = []
        v3_coord = []
        flag = []
        
        flux = np.empty((0, 0))  # Tableau vide 2D
        flux_err = np.empty((0, 0))
        flux_flag = np.empty((0, 0))
        flux_TEL = np.empty((0, 0))
        wavel_flux = np.empty((0, 0))
        
        STA_INDEX_tmp = []
        TEL_NAME_tmp  = []        
        
        for s in range(len(index_array)):

            STA_INDEX_tmp = np.append(STA_INDEX_tmp,np.array(fichier[index_array[s]].data['STA_INDEX']))               
            TEL_NAME_tmp = np.append(TEL_NAME_tmp,[np.array(fichier[index_array[s]].data['TEL_NAME'])])                

        STA_INDEX_tmp = np.reshape(STA_INDEX_tmp, (-1,4))
        TEL_NAME_tmp = np.reshape(TEL_NAME_tmp, (-1,4))

        
        for i in range(len(index_vis)):
            
            wavel_index =  np.where(insame_wvl == fichier[index_vis[i]].header['INSNAME'])[0][0] 
            len_wavel_vis = fichier[index_wvl[wavel_index]].header['NAXIS2']
            
            vis2 = np.append(vis2,np.array(fichier[index_vis[i]].data['VIS2DATA']))
            vis2_err = np.append(vis2_err,np.array(fichier[index_vis[i]].data['VIS2ERR']))
            vis2_flag = np.append(vis2_flag,np.array(fichier[index_vis[i]].data['FLAG']))
            u_coord = np.append(u_coord,np.array([fichier[index_vis[i]].data['UCOORD']]*len_wavel_vis).T)
            v_coord = np.append(v_coord,np.array([fichier[index_vis[i]].data['VCOORD']]*len_wavel_vis).T)
            
            wavel = np.append(wavel,[fichier[index_wvl[wavel_index]].data['EFF_WAVE']]*len(fichier[index_vis[i]].data['VIS2ERR']))
                
            vis2     = np.array(vis2, dtype='object').astype(float)
            vis2_err = np.array(vis2_err, dtype='object').astype(float)
            vis2_flag = np.array(vis2_flag, dtype='bool') # ==> MAGIC TRICK! Like this all the NaN value that has been added to fit the wavelengths will automatically be flagged! But we will conserve the exact shape for all data, MOUHAHAHA (yes it's 2am and I'm proud of it)
            u_coord  = np.array(u_coord, dtype='object').astype(float)
            v_coord  = np.array(v_coord, dtype='object').astype(float)
            wavel    = np.array(wavel, dtype='object').astype(float)
            flag     = np.array(flag, dtype='object').astype(float)
        
        dic['VIS2'] = {}
        dic['VIS2']['WAVEL']      = wavel               
        dic['VIS2']['BASELINE']   = (u_coord**2+v_coord**2)**(1/2)

        dic['VIS2']['U']          = u_coord    
        dic['VIS2']['V']          = v_coord  


        dic['VIS2']['VIS2']       = vis2  
        dic['VIS2']['VIS2_ERR']   = vis2_err
        dic['VIS2']['FLAG']       = vis2_flag

        wavel     = []

        for i in range(len(index_t3)):

            wavel_index =  np.where(insame_wvl == fichier[index_t3[i]].header['INSNAME'])[0][0] 
            len_wavel_t3 = fichier[index_wvl[wavel_index]].header['NAXIS2']
            
            t3       = np.append(t3,np.array(fichier[index_t3[i]].data['T3PHI']))
            t3_err   = np.append(t3_err,np.array(fichier[index_t3[i]].data['T3PHIERR']))
            t3_flag  = np.append(t3_flag,np.array(fichier[index_t3[i]].data['FLAG']))
            u1_coord = np.append(u1_coord,np.array([fichier[index_t3[i]].data['U1COORD']]*len_wavel_t3).T)
            v1_coord = np.append(v1_coord,np.array([fichier[index_t3[i]].data['V1COORD']]*len_wavel_t3).T)
            u2_coord = np.append(u2_coord,np.array([fichier[index_t3[i]].data['U2COORD']]*len_wavel_t3).T)
            v2_coord = np.append(v2_coord,np.array([fichier[index_t3[i]].data['V2COORD']]*len_wavel_t3).T)
            u3_coord = np.append(u3_coord,np.array([fichier[index_t3[i]].data['U1COORD']+fichier[index_t3[i]].data['U2COORD']]*len_wavel_t3).T)
            v3_coord = np.append(v3_coord,np.array([fichier[index_t3[i]].data['V1COORD']+fichier[index_t3[i]].data['V2COORD']]*len_wavel_t3).T)
            
            wavel = np.append(wavel,[fichier[index_wvl[wavel_index]].data['EFF_WAVE']]*len(fichier[index_t3[i]].data['T3PHI']))
            
            
            
            t3     = np.array(t3, dtype='object').astype(float)
            t3_err = np.array(t3_err, dtype='object').astype(float)
            t3_flag = np.array(t3_flag, dtype='bool') # ==> MAGIC TRICK! Like this all the NaN value that has been added to fit the wavelengths will automatically be flagged! But we will conserve the exact shape for all data, MOUHAHAHA (yes it's 2am and I'm proud of it)
            u1_coord  = np.array(u1_coord, dtype='object').astype(float)
            v1_coord  = np.array(v1_coord, dtype='object').astype(float)
            u2_coord  = np.array(u2_coord, dtype='object').astype(float)
            v2_coord  = np.array(v2_coord, dtype='object').astype(float)
            u3_coord  = np.array(u3_coord, dtype='object').astype(float)
            v3_coord  = np.array(v3_coord, dtype='object').astype(float)
            wavel    = np.array(wavel, dtype='object').astype(float)
            flag     = np.array(flag, dtype='object').astype(float)
    
        dic['T3'] = {}
        dic['T3']['WAVEL']      = wavel               
        dic['T3']['B1']   = (u1_coord**2+v1_coord**2)**(1/2)
        dic['T3']['B2']   = (u2_coord**2+v2_coord**2)**(1/2)
        dic['T3']['B3']   = (u3_coord**2+v3_coord**2)**(1/2)

        dic['T3']['U1']          = u1_coord    
        dic['T3']['V1']          = v1_coord  
        dic['T3']['U2']          = u2_coord    
        dic['T3']['V2']          = v2_coord  
        dic['T3']['U3']          = u3_coord    
        dic['T3']['V3']          = v3_coord  


        dic['T3']['T3']       = t3  
        dic['T3']['T3_ERR']   = t3_err
        dic['T3']['FLAG']     = t3_flag



        for i in range(len(index_flux)):
            # Récupérer l'index et la longueur d'onde correspondante
            wavel_index = np.where(insame_wvl == fichier[index_t3[i]].header['INSNAME'])[0][0]
            len_wavel_flux = fichier[index_wvl[wavel_index]].header['NAXIS2']
        
            # Extraire les données à ajouter
            new_flux = np.array(fichier[index_flux[i]].data['FLUXDATA'])
            new_flux_err = np.array(fichier[index_flux[i]].data['FLUXERR'])
            new_flux_flag = np.array(fichier[index_flux[i]].data['FLAG'])
            new_wavel_flux = np.array([fichier[index_wvl[wavel_index]].data['EFF_WAVE']] * 4)
        
            # Si flux est vide, affecter directement
            if flux.size == 0:
                flux = new_flux
                flux_err = new_flux_err
                flux_flag = new_flux_flag
                wavel_flux = new_wavel_flux
            else:
                # Vérifier les dimensions et ajuster si nécessaire pour flux, flux_err, flux_flag
                if flux.shape[1] != new_flux.shape[1]:
                    max_len = max(flux.shape[1], new_flux.shape[1])
                    flux = np.pad(flux, ((0, 0), (0, max_len - flux.shape[1])), constant_values=np.nan)
                    flux_err = np.pad(flux_err, ((0, 0), (0, max_len - flux_err.shape[1])), constant_values=np.nan)
                    flux_flag = np.pad(flux_flag, ((0, 0), (0, max_len - flux_flag.shape[1])), constant_values=np.nan)
                    new_flux = np.pad(new_flux, ((0, 0), (0, max_len - new_flux.shape[1])), constant_values=np.nan)
                    new_flux_err = np.pad(new_flux_err, ((0, 0), (0, max_len - new_flux_err.shape[1])), constant_values=np.nan)
                    new_flux_flag = np.pad(new_flux_flag, ((0, 0), (0, max_len - new_flux_flag.shape[1])), constant_values=np.nan)
        
                # Ajouter les nouvelles données
                flux = np.append(flux, new_flux, axis=0)
                flux_err = np.append(flux_err, new_flux_err, axis=0)
                flux_flag = np.append(flux_flag, new_flux_flag, axis=0)
        
                # Ajuster wavel_flux si nécessaire
                if wavel_flux.shape[1] != new_wavel_flux.shape[1]:
                    max_len_wvl = max(wavel_flux.shape[1], new_wavel_flux.shape[1])
                    wavel_flux = np.pad(wavel_flux, ((0, 0), (0, max_len_wvl - wavel_flux.shape[1])), constant_values=np.nan)
                    new_wavel_flux = np.pad(new_wavel_flux, ((0, 0), (0, max_len_wvl - new_wavel_flux.shape[1])), constant_values=np.nan)
        
                # Concaténer wavel_flux
                wavel_flux = np.append(wavel_flux, new_wavel_flux, axis=0)
        
            # Gestion des STA_INDEX et TEL_NAME
            index_good_conf = trouver_position_sous_array(STA_INDEX_tmp, np.array(fichier[index_flux[i]].data['STA_INDEX']))
            new_flux_TEL = np.array(
                [np.array(remove_spaces_in_string_array(TEL_NAME_tmp[index_good_conf])[
                    np.nonzero(np.array(fichier[index_flux[i]].data['STA_INDEX'])[:, None] == STA_INDEX_tmp[index_good_conf])[1]
                ])] * len_wavel_flux
            ).T
        
            # Ajuster flux_TEL si nécessaire
            if flux_TEL.size == 0:
                flux_TEL = new_flux_TEL
            else:
                if flux_TEL.shape[1] != new_flux_TEL.shape[1]:
                    max_len_tel = max(flux_TEL.shape[1], new_flux_TEL.shape[1])
                    flux_TEL = np.pad(flux_TEL, ((0, 0), (0, max_len_tel - flux_TEL.shape[1])), constant_values=np.nan)
                    new_flux_TEL = np.pad(new_flux_TEL, ((0, 0), (0, max_len_tel - new_flux_TEL.shape[1])), constant_values=np.nan)
        
                # Concaténer flux_TEL
                flux_TEL = np.append(flux_TEL, new_flux_TEL, axis=0)
        
        # Mise à jour du dictionnaire final
        dic['FLUX'] = {}
        dic['FLUX']['WAVEL'] = wavel_flux
        dic['FLUX']['FLUX'] = flux
        dic['FLUX']['FLUX_ERR'] = flux_err
        dic['FLUX']['FLAG'] = flux_flag
        dic['FLUX']['AT_NUMBER'] = flux_TEL

            

    OIFITS_TOT   = dic

    fichier.close()


    print('OIFITS READING OVER')
    
    return OIFITS_TOT

# filename = 'C:/Users/jdrevon/Desktop/Margaux/ALL_LM_TOT.fits'

# # filename = 'C:/Users/jdrevon/Desktop/SUPERVISION/2024/VLAD/data_issue/xtra_n_band2.fits'
# OIFITS_TOT = OIFITS_READING_concatenate(filename)

# import matplotlib.pyplot as plt

# plt.figure()
# plt.scatter(OIFITS_TOT['VIS2']['BASELINE']/OIFITS_TOT['VIS2']['WAVEL'],OIFITS_TOT['VIS2']['VIS2'], s=2)
# ax=plt.gca()
# ax.set_yscale('log')