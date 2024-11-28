import numpy as np
import os 
from astropy.io import fits 
import matplotlib
matplotlib.use('Agg') 
import random
import math
from extract_interferometric_quantities_from_image import IMAGE_to_V2
import matplotlib.pyplot as plt
import concurrent.futures
from DATA_reading import OIFITS_READING_concatenate
import sys
from scipy.ndimage import gaussian_filter

def rad_to_mas(rad):
    y = rad*1E3/(np.pi/(180*3600))
    return y    


path_data = "/home2/jdrevon/DATA_rec/FULL_3.05_3.12.fits"
path_results_tmp = "/home2/jdrevon/R_For_rec/"
target_name = 'RFor310'


num_cores= 30
hyperparameter_numb = 10 
h_min, h_max = 3, 9
iter_rec = 100

prior_image = "Dirac" #Dirac, random, path of the image
regularization = ['compactness'] # 'hyperbolic' or 'compactness' 

DATA_OBS = OIFITS_READING_concatenate(path_data)

mira_max_pix_size = rad_to_mas(0.5/(max(max(abs(DATA_OBS['VIS2']['U'])/DATA_OBS['VIS2']['WAVEL']),max(abs(DATA_OBS['VIS2']['V'])/DATA_OBS['VIS2']['WAVEL']))))

obs_res       =  rad_to_mas(0.5/(max(max(abs(DATA_OBS['VIS2']['U'])/DATA_OBS['VIS2']['WAVEL']),max(abs(DATA_OBS['VIS2']['V'])/DATA_OBS['VIS2']['WAVEL']))))
obs_hyper_res =  rad_to_mas(0.25/(max(np.sqrt(DATA_OBS['VIS2']['U']**2+DATA_OBS['VIS2']['V']**2)/DATA_OBS['VIS2']['WAVEL'])))
obs_FoV       =  rad_to_mas(1.0/(min(np.sqrt(DATA_OBS['VIS2']['U']**2+DATA_OBS['VIS2']['V']**2)/DATA_OBS['VIS2']['WAVEL'])))

print(obs_res, obs_hyper_res, obs_FoV)

# angular_size = np.floor(np.linspace(0.3*obs_FoV, 0.8*obs_FoV,5)).astype(int) # expected angular size of your object  
# pixelsize    = np.floor(np.linspace(obs_res/3, obs_res,5)*10**1)/10**1 #[0.3,0.5,0.7,1.0,1.3,1.5] #mas
# FoV          = np.floor(np.linspace(0.7*obs_FoV, 1.0*obs_FoV,4)).astype(int) #mas

# if all(element <= mira_max_pix_size for element in pixelsize):
    
#     print("All the pixelsize elements are inferior to the maximum pixelsize: %.5fmas computed by MiRA to avoid aliasing, we can continue."%mira_max_pix_size)

# else:

#     print("All the pixelsize elements ARE NOT inferior to the maximum pixelsize: %.5fmas computed by MiRA to avoid aliasing, please modify the values, run ABORTED."%mira_max_pix_size)
#     sys.exit(1)



# angular_size = [5] # expected angular size of your object  
# pixelsize = [0.3,0.5] #mas
# FoV = [40,50] #mas
# hyperparameter = [1E4,1E5]


ftol = 0
gtol = 0
xtol = 0
ftot = 1 #sum of total flux
maxeval = 20000
verb = 500


import shutil

def random_tau(tau, n_points=100):
    # Générer une liste de valeurs balayées log10 entre tau et tau * 10^4
    log10_sweep = np.logspace(np.log10(tau), np.log10(tau * 10**4), num=n_points)
    
    # Choisir une valeur aléatoire parmi ces valeurs balayées
    tau_final = np.random.choice(log10_sweep)
    
    return tau_final

def lister_dossiers(repertoire):
    dossiers = []
    for element in os.listdir(repertoire):
        if os.path.isdir(os.path.join(repertoire, element)):
            dossiers.append(element)
    return dossiers

def regroup_image(general_path):
    # general_path = path_results_tmp
    # general_path = "C:/Users/jdrevon/Desktop/MiRA_images/"
    
    noms_des_objets = lister_dossiers(general_path)
    
    for i in range(len(noms_des_objets)):
    
        reg = lister_dossiers(general_path + noms_des_objets[i])
        
        for j in range(len(reg)):
            
            path = general_path + noms_des_objets[i] + '/' + reg[j] + '/'
    
            path_param = lister_dossiers(path)
            
            if not os.path.exists(path+'ALL_images/'):
                os.mkdir(path+'ALL_images/')
            
            for k in range(len(path_param)):
                try:
                    shutil.copy(path+path_param[k]+'/reg_'+reg[j]+'_'+path_param[k]+"_image.png", path+'ALL_images/')
                except: 
                    None

    return

def mscatter(x,y,ax=None, m=None, **kw):
    import matplotlib.markers as mmarkers
    if not ax: ax=plt.gca()
    sc = ax.scatter(x,y,**kw)
    if (m is not None) and (len(m)==len(x)):
        paths = []
        for marker in m:
            if isinstance(marker, mmarkers.MarkerStyle):
                marker_obj = marker
            else:
                marker_obj = mmarkers.MarkerStyle(marker)
            path = marker_obj.get_path().transformed(
                        marker_obj.get_transform())
            paths.append(path)
        sc.set_paths(paths)
    return sc

# Définir une fonction pour exécuter le calcul pour une valeur d'hyperparamètre donnée
def run_calculation_compactness(reg, gam, pix, fov, h):

    path_param = path_reg + 'FoV_%.1f_pixsize_%.1f_gamma_%.1f/' % (fov, pix, gam)
    path_hyper = path_param + 'reg_%s_FoV_%.1f_pixsize_%.1f_gamma_%.1f_mu_%.1E.fits' % (reg, fov, pix, gam, h)

    os.system("ymira -pixelsize=%.1fmas -fov=%.1fmas -regul=%s -mu=%.1E -gamma=%imas" % (
        pix, fov, reg, h, gam) + \
              " -ftol=%.1f -gtol=%.1f -xtol=%.1f -maxeval=%i -verb=%i -overwrite -save_visibilities -save_initial " % (
                  ftol, gtol, xtol, maxeval, verb) + \
              " -flux=%.1f --use_vis=none -save_dirty_beam -save_dirty_map -save_residual_map --use_vis2=all --use_t3=all -initial=%s -recenter %s %s" % (ftot, prior_image, path_data, path_hyper))
    return path_hyper



def run_calculation_hyperbolic(reg, tau, pix, fov, h):

    path_param = path_reg + 'FoV_%.1f_pixsize_%.1f_tau_%.1E/' % (fov, pix, tau)
    path_hyper = path_param + 'reg_%s_FoV_%.1f_pixsize_%.1f_tau_%.1E_mu_%.1E.fits' % (reg, fov, pix, tau, h)

    os.system("ymira -pixelsize=%.1fmas -fov=%.1fmas -regul=%s -mu=%.1E -tau=%.1E" % (
        pix, fov, reg, h, tau) + \
              " -ftol=%.1f -gtol=%.1f -xtol=%.1f -maxeval=%i -verb=%i -overwrite -save_visibilities -save_initial" % (
                  ftol, gtol, xtol, maxeval, verb) + \
              " -flux=%.1f --use_vis=none -save_dirty_beam -save_dirty_map -save_residual_map --use_vis2=all --use_t3=all -initial=%s -recenter %s %s" % (ftot, prior_image, path_data, path_hyper))
    return path_hyper


# Creating the first folder in the folder architecture #regularization
for r in range(len(regularization)):
    path_reg = path_results_tmp + target_name + "/" + regularization[r] + "/"
    if not os.path.exists(path_reg):
        os.makedirs(path_reg)

    # Creating the second folder in the folder architecture #angular_size_pixel_size_FoV
    if regularization[r] == 'compactness':
        
        angular_size_p   = [] 
        pixsize_p        = []
        FoV_p            = []
        hyperparameter_p = []
        chi2_p           = []

        path_save = path_results_tmp + '/' + target_name + '/' + regularization[r] + '/'
        
        if not os.path.exists(path_save+'ALL_images/'):
            os.mkdir(path_save+'ALL_images/')
        

        for i in range(iter_rec):
            angular_size,pixelsize,FoV = round(random.uniform(0.3*obs_FoV, 0.8*obs_FoV), 3), math.floor(round(random.uniform(obs_res/3, obs_res), 3)*10**1)/10**1, round(random.uniform(0.7*obs_FoV, 1.0*obs_FoV), 3)
            path_param = path_reg + 'FoV_%.1f_pixsize_%.1f_gamma_%.1f/' % (FoV, pixelsize, angular_size)
                           
            if not os.path.exists(path_param):
                os.makedirs(path_param)

            hyperparameter=[]
            for _ in range(hyperparameter_numb):
                # Générer une puissance aléatoire entre les bornes
                log_aleatoire = random.uniform(h_min, h_max)
                
                # Calculer le nombre en utilisant 10 à la puissance du logarithme aléatoire
                nombre_aleatoire = 10 ** log_aleatoire
                
                # Ajouter le nombre à la liste
                hyperparameter.append(nombre_aleatoire)

            hyperparameter.sort()
            # Exécuter le calcul pour chaque valeur d'hyperparamètre en parallèle
            with concurrent.futures.ThreadPoolExecutor(max_workers=num_cores) as executor:
                futures = [executor.submit(run_calculation_compactness, regularization[r], angular_size, pixelsize,
                                           FoV, h) for h in hyperparameter]

                # Attendre que tous les calculs soient terminés
                for future in concurrent.futures.as_completed(futures):
                    try:
                        result = future.result()
                    except Exception as exc:
                        print(exc)
                        # Gérer toute exception survenue pendant le calcul

            # Tracer la figure
            
            path = []
            chi2 = []   
            image_tot = []
            converge = []
            neval = []
            conv_test = []
            color = []
            plot_style = []
            style = []
            gpnorm = []

            num_images = len(hyperparameter)
            num_cols = 6  # Nombre d'images par lignes
            num_rows = len(hyperparameter) # Calcul du nombre de lignes nécessaires                

            fig1, axes = plt.subplots(num_rows, num_cols, figsize=(45, 5 * num_rows))
#                    fig2, axes2 = plt.subplots(num_rows, num_cols, figsize=(20, 5 * num_rows))
#                    fig3, axes3 = plt.subplots(num_rows, num_cols, figsize=(20, 5 * num_rows))
#                    fig4, axes4 = plt.subplots(num_rows, num_cols, figsize=(20, 5 * num_rows))
            #fig5, axes5 = plt.subplots((num_rows+1)//2, 2, figsize=(10, 5 * (num_rows+1)//2))
#                    fig6, axes6 = plt.subplots(num_rows, num_cols, figsize=(10, 5 * num_rows))
            
            
            k=-1                

            for h in hyperparameter:
                k+=1
                
                
                path_hyper = path_param + 'reg_%s_FoV_%.1f_pixsize_%.1f_gamma_%.1f_mu_%.1E.fits' % (
                    regularization[r], FoV, pixelsize, angular_size, h)


                
                path.append(path_hyper)
                path_results= path_hyper
                hdul = fits.open(path_results)
                
                chi2.append(hdul[2].header['CHISQ'])
                converge.append(hdul[2].header['CONVERGE'])
                neval.append(hdul[2].header['NEVAL'])
                gpnorm.append(hdul[2].header['GPNORM'])
                
                angular_size_p.append(angular_size)
                pixsize_p.append(pixelsize)
                FoV_p.append(FoV)
                hyperparameter_p.append(h)
                chi2_p.append(hdul[2].header['CHISQ'])

                
                gpnorm_cond     = hdul[2].header['GPNORM']
                converge_cond   = hdul[2].header['CONVERGE'] 
                neval_cond      = hdul[2].header['NEVAL'] 
                
                #GREEN if converge, if num eval < max eval -1 , and if GPNORM < 1 
                
                note = ""
                
                if np.logical_and(np.logical_and(converge_cond == True, gpnorm_cond<1), neval_cond<maxeval-1) == True:
                    
                    conv_test.append(True)
                    color.append('chartreuse')
                    style.append('o')                                                   
                    note = 'C'
                    
                #ORANGE, WARNING if converge and if num eval >= max eval -1 or if GPNORM > 1 
                
                elif np.logical_and(converge_cond == True, np.logical_or(neval_cond >= maxeval-1, gpnorm_cond>1))== True:
                    
                    conv_test.append(False)
                    color.append('orange')
                    style.append('o')
                    note = 'W'
                #RED, non-converged if converge and if num eval >= max eval -1 or if GPNORM > 1 
 
                elif converge_cond == False :
                  
                    conv_test.append(False)
                    color.append('red')
                    style.append('x')
                    note = 'NC'


                q_MiRA, V2_MATISSE, V2_MATISSE_ERR, V2_MiRA, q_CP, CP_MiRA, CP_MATISSE = IMAGE_to_V2(path_data, path_results, 8) #compute the visibilities and closure phase

                image=fits.getdata(path_results)    
                shape = np.shape(image)
                size_image = shape[0]
                x_image = np.linspace(-shape[0]/2,shape[0]/2,shape[0],endpoint=False)*pixelsize
                
                row = k // 2
                col = k % 2
#                        axs  = axes[row, col] if num_rows > 1 else axes[col]
#                        axs2 = axes2[row, col] if num_rows > 1 else axes2[col]
#                        axs3 = axes3[row, col] if num_rows > 1 else axes3[col]
#                        axs4 = axes4[row, col] if num_rows > 1 else axes4[col]
                #axs5 = axes5[row, col] if num_rows > 1 else axes5[col]
#                        axs6 = axes6[row, col] if num_rows > 1 else axes6[col]

                axes[k,0].imshow(image,extent=[min(x_image),max(x_image),min(x_image),max(x_image)], cmap = 'hot')#,norm=colors.PowerNorm(gamma=0.5))#, norm=colors.PowerNorm(gamma=1.5))
                axes[k,0].set_xlabel(r'$\alpha$ [mas]',fontsize = 12)
                axes[k,0].set_ylabel(r'$\delta$ [mas]',fontsize = 12)
                axes[k,0].text(max(x_image)*40/100, max(x_image)*80/100, r'$\mu = %.1E $'%h, color='white')
                axes[k,0].text(min(x_image)+max(x_image)*10/100, max(x_image)*80/100, r'$\chi^2 = %.2E, %s $'%(hdul[2].header['CHISQ'], note), color='white')                        
                axes[k,0].minorticks_on()
                axes[k,0].tick_params(axis='x', labelsize=13)
                axes[k,0].tick_params(axis='y', labelsize=13)
                axes[k,0].set_xlim([min(x_image),max(x_image)])
                axes[k,0].set_ylim([min(x_image),max(x_image)])

                axes[k,2].errorbar(q_MiRA, V2_MATISSE, V2_MATISSE_ERR, label='Observations', alpha=0.05, c='tab:blue',  fmt='o', ecolor='red', ms=2)
                axes[k,2].scatter(q_MiRA, V2_MiRA, s=3, label = 'MiRA', c='tab:orange')
                axes[k,2].set_xlabel(r'B/$\lambda$ [rad$^-1$]',fontsize = 12)
                axes[k,2].set_ylabel(r'V$^2$',fontsize = 12)
                axes[k,2].legend()
                plt.tight_layout()

                axes[k,3].errorbar(q_MiRA, V2_MATISSE, V2_MATISSE_ERR, label='Observations', alpha=0.05, c='tab:blue',  fmt='o', ecolor='red', ms=2)
                axes[k,3].scatter(q_MiRA, V2_MiRA, s=3, label = 'MiRA', c='tab:orange')
                axes[k,3].set_xlabel(r'B/$\lambda$ [rad$^-1$]',fontsize = 12)
                axes[k,3].set_ylabel(r'V$^2$',fontsize = 12)
                axes[k,3].set_yscale('log')
                axes[k,3].legend()
                plt.tight_layout()

                axes[k,4].scatter(q_CP,np.array(CP_MATISSE), s=2, label='Observation', alpha=0.1)
                axes[k,4].scatter(q_CP,np.array(CP_MiRA), s=2, label='MiRA')
                axes[k,4].set_xlabel(r'B/$\lambda$ [rad$^-1$]',fontsize = 12)
                axes[k,4].set_ylabel('CP',fontsize = 12)
                axes[k,4].legend()
                plt.tight_layout()

                axes[k,5].scatter(q_MiRA, (V2_MATISSE-V2_MiRA)/V2_MATISSE_ERR, s=3)
                axes[k,5].set_xlabel(r'B/$\lambda$ [rad$^-1$]',fontsize = 12)
                axes[k,5].set_ylabel(r'V2 Pearson residuals',fontsize = 12)
                plt.tight_layout()

                p_mas = pixelsize  # pixel in mas
                FWHM_mas = obs_res  # mas
                FWHM_pixels = FWHM_mas / p_mas
                
                # Calcul de sigma en pixels
                sigma_pixels = FWHM_pixels / 2.35482
                
                
                image_conv = gaussian_filter(image, sigma=sigma_pixels)
        
                axes[k,1].imshow(image_conv,extent=[min(x_image),max(x_image),min(x_image),max(x_image)], cmap = 'hot')#,norm=colors.PowerNorm(gamma=0.5))#, norm=colors.PowerNorm(gamma=1.5))
                
                axes[k,1].set_xlabel(r'$\alpha$ [mas]',fontsize = 12)
                axes[k,1].set_ylabel(r'$\delta$ [mas]',fontsize = 12)
                axes[k,1].text(max(x_image)*40/100, max(x_image)*80/100, r'$\mu = %.1E $'%h, color='white')
                axes[k,1].text(min(x_image)+max(x_image)*10/100, max(x_image)*80/100, r'$\chi^2 = %.2E, %s $'%(hdul[2].header['CHISQ'], note), color='white')
                
            
                axes[k,1].minorticks_on()
                axes[k,1].tick_params(axis='x', labelsize=13)
                axes[k,1].tick_params(axis='y', labelsize=13)
                axes[k,1].set_xlim([min(x_image),max(x_image)])
                axes[k,1].set_ylim([min(x_image),max(x_image)])
                plt.tight_layout()
                
                #For Astropy Version < 6.0 fixed a bug
                
                with fits.open(path_results, mode='update', ignore_missing_end=True) as hdul:

                    header = hdul[0].header.copy()
                    
                    hdul[3].header.insert(5, ('PCOUNT', 0))
                    hdul[3].header.insert(6, ('GCOUNT', 1))                     
              
                    hdul[4].header.insert(5, ('PCOUNT', 0))
                    hdul[4].header.insert(6, ('GCOUNT', 1))
              
                    hdul[5].header.insert(5, ('PCOUNT', 0))
                    hdul[5].header.insert(6, ('GCOUNT', 1))
                    

                    data = image_conv  
                    new_image_hdu = fits.ImageHDU(data, header=header, name='CONVOLVED_IMAGE')

                    hdul.append(new_image_hdu)


            column_titles = ["Image", "Image Convolved", 'V2', 'V2 log', 'CP', 'V2 residuals']

            for c in range(len(column_titles)):
              axes[0,c].set_title(column_titles[c], fontsize=20)

            
            plt.subplots_adjust(left=0.4, right=0.99, top=0.92, bottom=0.4, wspace=0.1, hspace=0.2)
            plt.tight_layout(pad=1.0, rect=[0, 0, 1, 0.96])
            output_dir = path_param + 'reg_%s_FoV_%.1f_pixsize_%.1f_gamma_%.1f_image.png'%(regularization[r],FoV, pixelsize, angular_size)
            fig1.savefig(output_dir)

#                    output_dir2 = path_param + 'reg_%s_FoV_%.1f_pixsize_%.1f_gamma_%.1f_V2linear.png'%(regularization[r],FoV[f], pixelsize[p], gamma)
#                    fig2.savefig(output_dir2)

#                    output_dir3 = path_param + 'reg_%s_FoV_%.1f_pixsize_%.1f_gamma_%.1f_V2log.png'%(regularization[r],FoV[f], pixelsize[p], gamma)
#                    fig3.savefig(output_dir3)

#                    output_dir4 = path_param + 'reg_%s_FoV_%.1f_pixsize_%.1f_gamma_%.1f_residuals.png'%(regularization[r],FoV[f], pixelsize[p], gamma)
#                    fig4.savefig(output_dir4)

#                    output_dir5 = path_param + 'reg_%s_FoV_%.1f_pixsize_%.1f_gamma_%.1f_image_conv.png'%(regularization[r],FoV[f], pixelsize[p], gamma)
#                    fig5.savefig(output_dir5)

#                    output_dir6 = path_param + 'reg_%s_FoV_%.1f_pixsize_%.1f_gamma_%.1f_CP.png'%(regularization[r],FoV[f], pixelsize[p], gamma)
#                    fig6.savefig(output_dir6)


            plt.close(fig1)
#                    plt.close(fig2)
#                    plt.close(fig3)
#                    plt.close(fig4)
#                    plt.close(fig5)

#                    plt.close(fig6)
                
            # fig=plt.figure()
            # mscatter(hyperparameter,chi2, c=color, m=style)
            # ax=plt.gca()
            # ax.set_xlabel('Hyperparameter')
            # ax.set_ylabel(r'$\chi^2$')
            # ax.set_xscale('log')
            # ax.set_yscale('log')
            # fig.savefig(path_param + 'reg_%s_FoV_%.1f_pixsize_%.1f_gamma_%.1f_l_curve.png'%(regularization[r],FoV, pixelsize, angular_size), dpi=300)
            # plt.close(fig)
            
            
            
                
        fig31=plt.figure()
        plt.scatter(hyperparameter_p,chi2_p, c=np.array(angular_size_p), s=2, cmap='rainbow')
        plt.colorbar()
        ax=plt.gca()
        ax.set_xlabel('Hyperparameter')
        ax.set_ylabel(r'$\chi^2$')
        ax.set_xscale('log')
        ax.set_yscale('log')
        fig31.savefig(path_save + 'l_curve_gamma.png', dpi=300)
        plt.close(fig31)
        
        fig32=plt.figure()
        plt.scatter(hyperparameter_p,chi2_p, c=np.array(pixsize_p), s=2, cmap='rainbow')
        plt.colorbar()
        ax=plt.gca()
        ax.set_xlabel('Hyperparameter')
        ax.set_ylabel(r'$\chi^2$')
        ax.set_xscale('log')
        ax.set_yscale('log')
        fig32.savefig(path_save + 'l_curve_pixsize.png', dpi=300)
        plt.close(fig32)
                
        fig33=plt.figure()
        plt.scatter(hyperparameter_p,chi2_p, c=np.array(FoV_p), s=2, cmap='rainbow')
        plt.colorbar()
        ax=plt.gca()
        ax.set_xlabel('Hyperparameter')
        ax.set_ylabel(r'$\chi^2$')
        ax.set_xscale('log')
        ax.set_yscale('log')
        fig33.savefig(path_save + 'l_curve_FoV.png', dpi=300)
        plt.close(fig33)
                    

                                    
    if regularization[r]=='hyperbolic':

        path_reg = path_results_tmp + target_name + "/" + regularization[r] + "/"
        if not os.path.exists(path_reg):
            os.makedirs(path_reg)
        
        path_save = path_results_tmp + '/' + target_name + '/' + regularization[r] + '/'
        
        if not os.path.exists(path_save+'ALL_images/'):
            os.mkdir(path_save+'ALL_images/')

        tau_p            = [] 
        pixsize_p        = []
        FoV_p            = []
        hyperparameter_p = []
        chi2_p           = []
        
        for i in range(iter_rec):

          pixelsize,FoV = math.floor(round(random.uniform(obs_res/3, obs_res), 3)*10**1)/10**1, round(random.uniform(0.7*obs_FoV, 1.0*obs_FoV), 3)
                     
          tau_tmp = (pixelsize / FoV)**2 * 5E1
          
          tau = random_tau(tau_tmp)
          
          if tau == 0 : tau=1
  
          path_param = path_reg + 'FoV_%.1f_pixsize_%.1f_tau_%.1E/' % (FoV, pixelsize, tau)
          

              
          if not os.path.exists(path_param):
              os.makedirs(path_param)
  
          hyperparameter=[]
          for _ in range(hyperparameter_numb):
              # Générer une puissance aléatoire entre les bornes
              log_aleatoire = random.uniform(h_min, h_max)
              
              # Calculer le nombre en utilisant 10 à la puissance du logarithme aléatoire
              nombre_aleatoire = 10 ** log_aleatoire
              
              # Ajouter le nombre à la liste
              hyperparameter.append(nombre_aleatoire)
              
          hyperparameter.sort()
          
          # Exécuter le calcul pour chaque valeur d'hyperparamètre en parallèle
          with concurrent.futures.ThreadPoolExecutor(max_workers=num_cores) as executor:
              futures = [executor.submit(run_calculation_hyperbolic, regularization[r], tau, pixelsize,
                                         FoV, h) for h in hyperparameter]
  
              # Attendre que tous les calculs soient terminés
              for future in concurrent.futures.as_completed(futures):
                  try:
                      result = future.result()
                  except Exception as exc:
                      print(exc)
                      # Gérer toute exception survenue pendant le calcul
  
          # Tracer la figure
  
          path = []
          chi2 = []   
          image_tot = []
          converge = []
          neval = []
          conv_test = []
          color = []
          plot_style = []
          style = []
          gpnorm = []
  
  
          num_images = len(hyperparameter)
          num_cols = 6  # Nombre d'images par lignes
          num_rows = len(hyperparameter) # Calcul du nombre de lignes nécessaires                
  
          fig1, axes = plt.subplots(num_rows, num_cols, figsize=(45, 5 * num_rows))
  #                    fig2, axes2 = plt.subplots(num_rows, num_cols, figsize=(20, 5 * num_rows))
  #                    fig3, axes3 = plt.subplots(num_rows, num_cols, figsize=(20, 5 * num_rows))
  #                    fig4, axes4 = plt.subplots(num_rows, num_cols, figsize=(20, 5 * num_rows))
          #fig5, axes5 = plt.subplots((num_rows+1)//2, 2, figsize=(10, 5 * (num_rows+1)//2))
  #                    fig6, axes6 = plt.subplots(num_rows, num_cols, figsize=(10, 5 * num_rows))
              
          
          
          k=-1                
  
          for h in hyperparameter:
              k+=1     
  
              path_hyper = path_param + 'reg_%s_FoV_%.1f_pixsize_%.1f_tau_%.1E_mu_%.1E.fits' % (
                  regularization[r], FoV, pixelsize, tau, h)
                            
    
              path.append(path_hyper)
              path_results = path_hyper
              hdul = fits.open(path_results)
              
              chi2.append(hdul[2].header['CHISQ'])
              converge.append(hdul[2].header['CONVERGE'])
              neval.append(hdul[2].header['NEVAL'])
              gpnorm.append(hdul[2].header['GPNORM'])
              
              
              tau_p.append(tau)
              pixsize_p.append(pixelsize)
              FoV_p.append(FoV)
              hyperparameter_p.append(h)
              chi2_p.append(hdul[2].header['CHISQ'])
              
              
              gpnorm_cond     = hdul[2].header['GPNORM']
              converge_cond   = hdul[2].header['CONVERGE'] 
              neval_cond      = hdul[2].header['NEVAL'] 
              
              #GREEN if converge, if num eval < max eval -1 , and if GPNORM < 1 
              
              note = ""
              
              if np.logical_and(np.logical_and(converge_cond == True, gpnorm_cond<1), neval_cond<maxeval-1) == True:
                  
                  conv_test.append(True)
                  color.append('chartreuse')
                  style.append('o')                                                   
                  note = 'C'
                  
              #ORANGE, WARNING if converge and if num eval >= max eval -1 or if GPNORM > 1 
              
              elif np.logical_and(converge_cond == True, np.logical_or(neval_cond >= maxeval-1, gpnorm_cond>1))== True:
                  
                  conv_test.append(False)
                  color.append('orange')
                  style.append('o')
                  note = 'W'
              #RED, non-converged if converge and if num eval >= max eval -1 or if GPNORM > 1 
  
              elif converge_cond == False :
                
                  conv_test.append(False)
                  color.append('red')
                  style.append('x')
                  note = 'NC'
  
  
              q_MiRA, V2_MATISSE, V2_MATISSE_ERR, V2_MiRA, q_CP, CP_MiRA, CP_MATISSE = IMAGE_to_V2(path_data, path_results, 8) #compute the visibilities and closure phase
  
              image=fits.getdata(path_results)    
              shape = np.shape(image)
              size_image = shape[0]
              x_image = np.linspace(-shape[0]/2,shape[0]/2,shape[0],endpoint=False)*pixelsize
              
              row = k // 2
              col = k % 2
  #                        axs  = axes[row, col] if num_rows > 1 else axes[col]
  #                        axs2 = axes2[row, col] if num_rows > 1 else axes2[col]
  #                        axs3 = axes3[row, col] if num_rows > 1 else axes3[col]
  #                        axs4 = axes4[row, col] if num_rows > 1 else axes4[col]
              #axs5 = axes5[row, col] if num_rows > 1 else axes5[col]
  #                        axs6 = axes6[row, col] if num_rows > 1 else axes6[col]
  
              axes[k,0].imshow(image,extent=[min(x_image),max(x_image),min(x_image),max(x_image)], cmap = 'hot')#,norm=colors.PowerNorm(gamma=0.5))#, norm=colors.PowerNorm(gamma=1.5))
              axes[k,0].set_xlabel(r'$\alpha$ [mas]',fontsize = 12)
              axes[k,0].set_ylabel(r'$\delta$ [mas]',fontsize = 12)
              axes[k,0].text(max(x_image)*40/100, max(x_image)*80/100, r'$\mu = %.1E $'%h, color='white')
              axes[k,0].text(min(x_image)+max(x_image)*10/100, max(x_image)*80/100, r'$\chi^2 = %.2E, %s $'%(hdul[2].header['CHISQ'], note), color='white')                        
              axes[k,0].minorticks_on()
              axes[k,0].tick_params(axis='x', labelsize=13)
              axes[k,0].tick_params(axis='y', labelsize=13)
              axes[k,0].set_xlim([min(x_image),max(x_image)])
              axes[k,0].set_ylim([min(x_image),max(x_image)])
  
              axes[k,2].errorbar(q_MiRA, V2_MATISSE, V2_MATISSE_ERR, label='Observations', alpha=0.05, c='tab:blue',  fmt='o', ecolor='red', ms=2)
              axes[k,2].scatter(q_MiRA, V2_MiRA, s=3, label = 'MiRA', c='tab:orange')
              axes[k,2].set_xlabel(r'B/$\lambda$ [rad$^-1$]',fontsize = 12)
              axes[k,2].set_ylabel(r'V$^2$',fontsize = 12)
              axes[k,2].legend()
              plt.tight_layout()
  
              axes[k,3].errorbar(q_MiRA, V2_MATISSE, V2_MATISSE_ERR, label='Observations', alpha=0.05, c='tab:blue',  fmt='o', ecolor='red', ms=2)
              axes[k,3].scatter(q_MiRA, V2_MiRA, s=3, label = 'MiRA', c='tab:orange')
              axes[k,3].set_xlabel(r'B/$\lambda$ [rad$^-1$]',fontsize = 12)
              axes[k,3].set_ylabel(r'V$^2$',fontsize = 12)
              axes[k,3].set_yscale('log')
              axes[k,3].legend()
              plt.tight_layout()
  
              axes[k,4].scatter(q_CP,np.array(CP_MATISSE), s=2, label='Observation', alpha=0.1)
              axes[k,4].scatter(q_CP,np.array(CP_MiRA), s=2, label='MiRA')
              axes[k,4].set_xlabel(r'B/$\lambda$ [rad$^-1$]',fontsize = 12)
              axes[k,4].set_ylabel('CP',fontsize = 12)
              axes[k,4].legend()
              plt.tight_layout()
  
              axes[k,5].scatter(q_MiRA, (V2_MATISSE-V2_MiRA)/V2_MATISSE_ERR, s=3)
              axes[k,5].set_xlabel(r'B/$\lambda$ [rad$^-1$]',fontsize = 12)
              axes[k,5].set_ylabel(r'V2 Pearson residuals',fontsize = 12)
              plt.tight_layout()
  
              p_mas = pixelsize  # pixel in mas
              FWHM_mas = obs_res  # mas
              FWHM_pixels = FWHM_mas / p_mas
              
              # Calcul de sigma en pixels
              sigma_pixels = FWHM_pixels / 2.35482
                            
              image_conv = gaussian_filter(image, sigma=sigma_pixels)
      
              axes[k,1].imshow(image_conv,extent=[min(x_image),max(x_image),min(x_image),max(x_image)], cmap = 'hot')#,norm=colors.PowerNorm(gamma=0.5))#, norm=colors.PowerNorm(gamma=1.5))
              
              axes[k,1].set_xlabel(r'$\alpha$ [mas]',fontsize = 12)
              axes[k,1].set_ylabel(r'$\delta$ [mas]',fontsize = 12)
              axes[k,1].text(max(x_image)*40/100, max(x_image)*80/100, r'$\mu = %.1E $'%h, color='white')
              axes[k,1].text(min(x_image)+max(x_image)*10/100, max(x_image)*80/100, r'$\chi^2 = %.2E, %s $'%(hdul[2].header['CHISQ'], note), color='white')
              
          
              axes[k,1].minorticks_on()
              axes[k,1].tick_params(axis='x', labelsize=13)
              axes[k,1].tick_params(axis='y', labelsize=13)
              axes[k,1].set_xlim([min(x_image),max(x_image)])
              axes[k,1].set_ylim([min(x_image),max(x_image)])
              plt.tight_layout()
  
  
              with fits.open(path_results, mode='update', verify='fix') as hdul:
                # Copier le header du premier HDU
                header = hdul[0].header.copy()
                      
                hdul[3].header.insert(5, ('PCOUNT', 0))
                hdul[3].header.insert(6, ('GCOUNT', 1))
                              
                hdul[4].header.insert(5, ('PCOUNT', 0))
                hdul[4].header.insert(6, ('GCOUNT', 1))
                
                hdul[5].header.insert(5, ('PCOUNT', 0))
                hdul[5].header.insert(6, ('GCOUNT', 1))
                      
                # Créer une nouvelle image avec le header copié et spécifier le nom de l'extension
                data = image_conv  # Image 2D de 100x100 pixels
                new_image_hdu = fits.ImageHDU(data, header=header, name='CONVOLVED_IMAGE')
                
                      
                # Ajouter la nouvelle image à la liste des HDUs
                hdul.append(new_image_hdu)
  
  
  
  
          column_titles = ["Image", "Image Convolved", 'V2', 'V2 log', 'CP', 'V2 residuals']
  
          for c in range(len(column_titles)):
            axes[0,c].set_title(column_titles[c], fontsize=20)
  
          
          plt.subplots_adjust(left=0.4, right=0.99, top=0.92, bottom=0.4, wspace=0.1, hspace=0.2)
          plt.tight_layout(pad=1.0, rect=[0, 0, 1, 0.96])
          output_dir = path_param + 'reg_%s_FoV_%.1f_pixsize_%.1f_tau_%.1E_image.png'%(regularization[r], FoV, pixelsize, tau)
          fig1.savefig(output_dir)
  
  #                output_dir2 = path_param + 'reg_%s_FoV_%.1f_pixsize_%.1f_tau_%.1E_V2linear.png'%(regularization[r], FoV[f], pixelsize[p], tau)
  #                fig2.savefig(output_dir2)
  
  #                output_dir3 = path_param + 'reg_%s_FoV_%.1f_pixsize_%.1f_tau_%.1E_V2log.png'%(regularization[r], FoV[f], pixelsize[p], tau)
  #                fig3.savefig(output_dir3)
  
  #                output_dir4 = path_param + 'reg_%s_FoV_%.1f_pixsize_%.1f_tau_%.1E_residuals.png'%(regularization[r], FoV[f], pixelsize[p], tau)
  #                fig4.savefig(output_dir4)
  
  #                output_dir5 = path_param + 'reg_%s_FoV_%.1f_pixsize_%.1f_tau_%.1E_image_conv.png'%(regularization[r], FoV[f], pixelsize[p], tau)
  #                fig5.savefig(output_dir5)
  
  #                output_dir6 = path_param + 'reg_%s_FoV_%.1f_pixsize_%.1f_tau_%.1E_CP.png'%(regularization[r], FoV[f], pixelsize[p], tau)
  #                fig6.savefig(output_dir6)
  
  
          plt.close(fig1)
  #                plt.close(fig2)
  #                plt.close(fig3)
  #                plt.close(fig4)
  #                plt.close(fig5)
  #                plt.close(fig6)
              
      
          # fig=plt.figure()
          # mscatter(hyperparameter,chi2, c=color, m=style)
          # ax=plt.gca()
          # ax.set_xlabel('Hyperparameter')
          # ax.set_ylabel(r'$\chi^2$')
          # ax.set_xscale('log')
          # ax.set_yscale('log')
          # fig.savefig(path_param + 'reg_%s_FoV_%.1f_pixsize_%.1f_tau_%.1E_l_curve.png'%(regularization[r], FoV, pixelsize, tau), dpi=300)
          # plt.close(fig)
        
        
        fig31=plt.figure()
        plt.scatter(hyperparameter_p,chi2_p, c=np.array(tau_p), s=2, cmap='rainbow')
        plt.colorbar()
        ax=plt.gca()
        ax.set_xlabel('Hyperparameter')
        ax.set_ylabel(r'$\chi^2$')
        ax.set_xscale('log')
        ax.set_yscale('log')
        fig31.savefig(path_save + 'l_curve_tau.png', dpi=300)
        plt.close(fig31)
        
        fig32=plt.figure()
        plt.scatter(hyperparameter_p,chi2_p, c=np.array(pixsize_p), s=2, cmap='rainbow')
        plt.colorbar()
        ax=plt.gca()
        ax.set_xlabel('Hyperparameter')
        ax.set_ylabel(r'$\chi^2$')
        ax.set_xscale('log')
        ax.set_yscale('log')
        fig32.savefig(path_save + 'l_curve_pixsize.png', dpi=300)
        plt.close(fig32)
                
        fig33=plt.figure()
        plt.scatter(hyperparameter_p,chi2_p, c=np.array(FoV_p), s=2, cmap='rainbow')
        plt.colorbar()
        ax=plt.gca()
        ax.set_xlabel('Hyperparameter')
        ax.set_ylabel(r'$\chi^2$')
        ax.set_xscale('log')
        ax.set_yscale('log')
        fig33.savefig(path_save + 'l_curve_FoV.png', dpi=300)
        plt.close(fig33)
    
        
        
regroup_image(path_results_tmp)                              