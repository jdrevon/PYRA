# PYRA
Python for MiRA

This code is a Python layer to make running MiRA. 
It needs Python 3.8+ installed and a Linux/Ubuntu environment to work.
The code has been created to scan randomly and in the more efficient way the parameter space to make image reconstruction using the MiRA code from the JMMC. 

The parameters that will be scan by the software are: 

the pixelsize, the field of view (FoV), the hyperparameter, and depending on the regularization if it is either the compactness or hyperbolic the code will vary respectively the gamma (full width at half maximum (FWHM) of the prior distribution of light) paramater and the tau (edge preserving threshold) parameter (for more information about those two last parameters, please look at the MiRA github webpage).

The pixelsize search interval is defined based on the interferometric resolution of the user observations. Assuming the hyper-resolution to be lambda_min/4B_max, and knowing that to avoid aliasing issue, MiRA does not accept pixelsize larger than lambda_min/2B_max, the parameter is randomly chosen in between [lambda/6B, lambda/2B]. By experience, it is useful also to probe smaller pixelsize than the interferometric pixel-size at lambda/4B that is why lambda/6B is used,

The FoV search interval is defined based on the interferometric field of view of the user observations defined as FoV=lambda_max/B_min. The parameter is then randomly chosen in between [0.7*FoV,1*FoV],

The tau (edge preserving threshold), if we assume an uniform image, the mean pixel value should be equal to flux*(pixelsize/FoV)^2, assuming a nromalized flux, we have then (pixelsize/FoV)^2. So this value is used as the lower boundary of the interval. The higher boundary has been arbitrarly defined as (pixelsize/FoV)^2 * 10^4. So the tau value is randomly searched in between [(pixelsize/FoV)^2, (pixelsize/FoV)^2*10^4],

The gamma (full width at half maximum (FWHM) of the prior distribution of light), the boundary as been arbitraly defined as: [0.3*FoV,0.8*FoV].

The hyperparameter (mu) is the weight given to the regularization function in the given formula for the chi2 minimization: chi2_TOT = chi2_DATA + mu*f_regularization. The mu interval is defined by the user but the ideal value is usually between 10^2 and 10^6. We usually first want to scan in between [10^3,10^9] to have a broad overview of the impact of our hyperparameter on the output of the image reconstruction. For a given s
For a given set of parameter (pixelsize,FoV,gamma) or (pixelsize,FoV,tau), the code will scan randomly the hyperparameter space as many times as the user wants using the variable hyperparameter_numb described below.  

How to use:

The main file is OBJECT_INSTRU_MiRA_v3.py:

The variable that the user can changed are :

path_data : full path of the .fits file containing all the observations (if the data are separated in several .fits files please merge them using your favorite python code or using OIFits Explorer software from the JMMC),
path_results_tmp : full path where you want to store the result from the image reconstructions,
target_name : this will be the name of the subfolder where all the different output will be stored,
num_cores: number of cores that will be dedicated to the routine,
hyperparameters_numb : this is the number of hyperparameters that will be generated for a single set of parameter. Up to now it is useless to put associated to the num_cores variables a number which is higher than the hyperprameters_numb because the paralalization occurs during the hyperparameter iterations within a given set of parameters.
h_min, h_max: this are the power of ten number that will be used by the code to define the boundaries of the hyperprameter space,
iter_rec: this is the number of set of parameter that you want the code to generate. Ex: hyperparameters_numb=10, iter_rec = 100, the code will generate in total 1000 different models. 
prior_image : the user can defined either he wants to use a "Dirac" or either he wants to use a model image under a .fits file.
regularization: here is the name of the regularization the user wants to probe, you can chose to put one ['compactness'] or ['hyperbolic'] or both ['compactness','hyperbolic'].
maxeval: this is the maximum of evaluation function to be done if the code has not converged yet before MiRA decides to stop (default value in this code: 20000).

Once you are done with the variable, just run the script.

All the output are stored in path_result_tmp/regularization/ and are stored in a subfolder with a name containing the variable names and their associate values.
You have also in ALL_IMAGES folder all the output in images with all the image reconstruction, the image reconstruction convolved by the interferometric beam, comparison between the reconstructed and the observed squared visibilities and closure phase as well as the residuals of the squared visibilities for each set of parameters and for each hyperparameter values. 
