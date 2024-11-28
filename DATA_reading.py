# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 21:45:17 2024

@author: jdrevon
"""

import numpy as np
from astropy.io import fits 

def OIFITS_READING_concatenate(filename):
    
    OIFITS_TOT = np.zeros(1, dtype=object)
        
    with fits.open(filename, memmap=False) as fichier:

        name_HDU = np.array([fichier[t].name for t in range(len(fichier))])
        
        index_vis   = np.where(name_HDU=='OI_VIS2')[0]
        index_t3    = np.where(name_HDU=='OI_T3')[0]
        index_wvl   = np.where(name_HDU=='OI_WAVELENGTH')[0]

        lambda_fichier   =  fichier[index_wvl[0]].data['EFF_WAVE']
        lambda_bandwidth =  fichier[index_wvl[0]].data['EFF_BAND']
        dic                            = {'WAVEL': lambda_fichier}
        dic['NAME']                    = filename         
        dic['BANDWIDTH']               = lambda_bandwidth        

               
        wavel_shape = np.array([len(fichier[t].data['EFF_WAVE']) for t in index_wvl])

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
        
        
        
        for i in range(len(index_vis)):
            
            if type(fichier[index_vis[i]].data['VIS2DATA'][0]) == np.float64:
                len_wavel_vis = 1
            else:
                len_wavel_vis = len(fichier[index_vis[i]].data['VIS2DATA'][0])
            
            if len_wavel_vis == max(wavel_shape):
                vis2= np.append(vis2,np.array(fichier[index_vis[i]].data['VIS2DATA']))
                vis2_err = np.append(vis2_err,np.array(fichier[index_vis[i]].data['VIS2ERR']))
                vis2_flag = np.append(vis2_flag,np.array(fichier[index_vis[i]].data['FLAG']))
                u_coord = np.append(u_coord,np.array([fichier[index_vis[i]].data['UCOORD']]*len_wavel_vis).T)
                v_coord = np.append(v_coord,np.array([fichier[index_vis[i]].data['VCOORD']]*len_wavel_vis).T)
                
                index=index_wvl[wavel_shape==len_wavel_vis][0] 
                
                wavel = np.append(wavel,[fichier[index].data['EFF_WAVE']]*len(fichier[index_vis[i]].data['VIS2ERR']))
            
            else:
                diff = max(wavel_shape)-len_wavel_vis
        
                shape = np.shape(np.array(fichier[index_vis[i]].data['VIS2DATA']))[0]
                add = np.empty((shape,diff))*np.nan
        
                vis2 = np.append(vis2,np.hstack((np.array(fichier[index_vis[i]].data['VIS2DATA']),add)), axis=0)
                vis2_err = np.append(vis2_err,np.hstack((np.array(fichier[index_vis[i]].data['VIS2ERR']),add)), axis=0)
                vis2_flag = np.append(vis2_flag,np.hstack((np.array(fichier[index_vis[i]].data['FLAG']),add)), axis=0)
                u_coord = np.append(u_coord,np.hstack((np.array([fichier[index_vis[i]].data['UCOORD']]*len_wavel_vis).T,add)), axis=0)
                v_coord = np.append(v_coord,np.hstack((np.array([fichier[index_vis[i]].data['VCOORD']]*len_wavel_vis).T,add)), axis=0)

                index=index_wvl[wavel_shape==len_wavel_vis][0] 
                
                wavel = np.append(wavel,np.hstack((np.array([fichier[index].data['EFF_WAVE']]*len(fichier[index_vis[i]].data['VIS2ERR'])),add)), axis=0)
            
            
            
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

            if type(fichier[index_t3[i]].data['T3PHI'][0]) == np.float64:
                len_wavel_t3 = 1
            else:
                len_wavel_t3 = len(fichier[index_t3[i]].data['T3PHI'][0])
            
            if len_wavel_t3 == max(wavel_shape):
                t3       = np.append(t3,np.array(fichier[index_t3[i]].data['T3PHI']))
                t3_err   = np.append(t3_err,np.array(fichier[index_t3[i]].data['T3PHIERR']))
                t3_flag  = np.append(t3_flag,np.array(fichier[index_t3[i]].data['FLAG']))
                u1_coord = np.append(u1_coord,np.array([fichier[index_t3[i]].data['U1COORD']]*len_wavel_t3).T)
                v1_coord = np.append(v1_coord,np.array([fichier[index_t3[i]].data['V1COORD']]*len_wavel_t3).T)
                u2_coord = np.append(u2_coord,np.array([fichier[index_t3[i]].data['U2COORD']]*len_wavel_t3).T)
                v2_coord = np.append(v2_coord,np.array([fichier[index_t3[i]].data['V2COORD']]*len_wavel_t3).T)
                u3_coord = np.append(u3_coord,np.array([fichier[index_t3[i]].data['U1COORD']+fichier[index_t3[i]].data['U2COORD']]*len_wavel_t3).T)
                v3_coord = np.append(v3_coord,np.array([fichier[index_t3[i]].data['V1COORD']+fichier[index_t3[i]].data['V2COORD']]*len_wavel_t3).T)
                
                index = index_wvl[wavel_shape==len_wavel_t3][0] 
                                
                wavel = np.append(wavel,[fichier[index].data['EFF_WAVE']]*len(fichier[index_t3[i]].data['T3PHI']))
            else:
                diff = max(wavel_shape)-len((np.array(fichier[index_t3[i]].data['T3PHI'])[0]))
        
                shape = np.shape(np.array(fichier[index_vis[i]].data['T3PHI']))[0]
                add = np.empty((shape,diff))*np.nan
        
                t3 = np.append(t3,np.hstack((np.array(fichier[index_t3[i]].data['T3PHI']),add)), axis=0)
                t3_err = np.append(t3_err,np.hstack((np.array(fichier[index_t3[i]].data['T3PHIERR']),add)), axis=0)
                t3_flag = np.append(t3_flag,np.hstack((np.array(fichier[index_t3[i]].data['FLAG']),add)), axis=0)
                u1_coord = np.append(u1_coord,np.hstack((np.array([fichier[index_t3[i]].data['U1COORD']]*len_wavel_t3).T,add)), axis=0)
                v1_coord = np.append(v1_coord,np.hstack((np.array([fichier[index_t3[i]].data['V1COORD']]*len_wavel_t3).T,add)), axis=0)
                u2_coord = np.append(u2_coord,np.hstack((np.array([fichier[index_t3[i]].data['U2COORD']]*len_wavel_t3).T,add)), axis=0)
                v2_coord = np.append(v2_coord,np.hstack((np.array([fichier[index_t3[i]].data['V2COORD']]*len_wavel_t3).T,add)), axis=0)
                u3_coord = np.append(u3_coord,np.hstack((np.array([fichier[index_t3[i]].data['U1COORD']+fichier[index_t3[i]].data['U2COORD']]*len_wavel_t3).T,add)), axis=0)
                v3_coord = np.append(v3_coord,np.hstack((np.array([fichier[index_t3[i]].data['V1COORD']+fichier[index_t3[i]].data['V2COORD']]*len_wavel_t3).T,add)), axis=0)

                index = index_wvl[wavel_shape==len_wavel_t3][0] 
                
                wavel = np.append(wavel,[fichier[index].data['EFF_WAVE']]*len(fichier[index_t3[i]].data['T3PHI']))
            
            
            
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


            

    OIFITS_TOT   = dic

    fichier.close()


    print('OIFITS READING OVER')
    
    return OIFITS_TOT

# filename = 'C:/partage/LM_ALL_4.096_4.11.fits'

# # filename = 'C:/Users/jdrevon/Desktop/SUPERVISION/2024/VLAD/data_issue/xtra_n_band2.fits'
# OIFITS_TOT = OIFITS_READING_concatenate(filename)

# import matplotlib.pyplot as plt

# plt.figure()
# plt.scatter(OIFITS_TOT['VIS2']['BASELINE']/OIFITS_TOT['VIS2']['WAVEL'],OIFITS_TOT['VIS2']['VIS2'], s=2)
# ax=plt.gca()
# ax.set_yscale('log')