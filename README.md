# PYRA: Python for MiRA

**PYRA** is a Python wrapper designed to facilitate and optimize the use of **MiRA** (Multi-aperture Image Reconstruction Algorithm), developed by the [JMMC](https://github.com/emmt/MiRA?tab=readme-ov-file). This tool enables efficient and random scanning of the parameter space for image reconstruction. 

**Requirements:**
- Python 3.8+
- Linux/Ubuntu environment

---

## Overview
PYRA helps automate and enhance MiRA's functionality by exploring key parameters related to image reconstruction. The parameters scanned include:
- **Pixel size**
- **Field of View (FoV)**
- **Hyperparameter (\( \mu \))**
- **Regularization-specific parameters:**
  - For *Compactness*: **gamma** (FWHM of the prior distribution of light)
  - For *Hyperbolic*: **tau** (edge-preserving threshold)

These parameters are selected within specific intervals, ensuring compatibility with interferometric observations.

---

## Parameter Search Intervals

# PYRA: Python for MiRA

**Parameter Search Intervals**

1. **Pixel Size (\(\text{pixelsize}\))**
   Defined by the interferometric resolution:
   `[λ_min / 6B_max, λ_min / 2B_max]`

The pixelsize search interval is defined based on the interferometric resolution of the user observations. Assuming the hyper-resolution to be lambda_min/4B_max, and knowing that to avoid aliasing issue, MiRA does not accept pixelsize larger than lambda_min/2B_max, the parameter is randomly chosen in between [lambda/6B, lambda/2B]. By experience, it is useful also to probe smaller pixelsize than the interferometric pixel-size at lambda/4B that is why lambda/6B is used,


2. **Field of View (FoV)**
   Defined based on the interferometric field of view:
   `[0.7 × FoV, 1 × FoV]`

The FoV search interval is defined based on the interferometric field of view of the user observations defined as FoV=lambda_max/B_min. The parameter is then randomly chosen in between [0.7*FoV,1*FoV],

3. **Tau (\(\tau\))**
   Edge-preserving threshold for hyperbolic regularization:
   `[(pixelsize/FoV)^2, (pixelsize/FoV)^2 × 10^4]`

The tau (edge preserving threshold), if we assume an uniform image, the mean pixel value should be equal to flux*(pixelsize/FoV)^2, assuming a nromalized flux, we have then (pixelsize/FoV)^2. So this value is used as the lower boundary of the interval. The higher boundary has been arbitrarly defined as (pixelsize/FoV)^2 * 10^4. So the tau value is randomly searched in between [(pixelsize/FoV)^2, (pixelsize/FoV)^2*10^4],

4. **Gamma (\(\gamma\))**
   FWHM of prior light distribution for compactness regularization:
   `[0.3 × FoV, 0.8 × FoV]`
   
The gamma (full width at half maximum (FWHM) of the prior distribution of light), the boundary as been arbitraly defined as: [0.3*FoV,0.8*FoV].

5. **Hyperparameter (\(\mu\))**
   Weight for the regularization function in the total chi-square minimization:
   Typically `[10^3, 10^9]`

The hyperparameter (mu) is the weight given to the regularization function in the given formula for the chi2 minimization: chi2_TOT = chi2_DATA + mu*f_regularization. The mu interval is defined by the user but the ideal value is usually between 10^2 and 10^6. We usually first want to scan in between [10^3,10^9] to have a broad overview of the impact of our hyperparameter on the output of the image reconstruction. For a given s
For a given set of parameter (pixelsize,FoV,gamma) or (pixelsize,FoV,tau), the code will scan randomly the hyperparameter space as many times as the user wants using the variable hyperparameter_numb described below.  

---

## How to Use

### Main Script: `OBJECT_INSTRU_MiRA_v3.py`
The primary script contains variables you need to configure before running. Below is a list of key variables and their purposes:

| **Variable**             | **Description**                                                                                                                                                                                                                              |
|--------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `path_data`              | Full path to the `.fits` file containing observations. If multiple `.fits` files exist, merge them before running the script (e.g., using Python or [OIFits Explorer](https://www.jmmc.fr/)).                                               |
| `path_results_tmp`       | Path to the directory for storing results.                                                                                                                                                                                                  |
| `target_name`            | Subfolder name where outputs will be stored.                                                                                                                                                                                                |
| `num_cores`              | Number of CPU cores dedicated to the routine.                                                                                                                                                                                               |
| `hyperparameters_numb`   | Number of hyperparameter iterations for each parameter set. Ensure this value is greater than `num_cores`.                                                                                                                                  |
| `h_min`, `h_max`         | Define the bounds of the hyperparameter space as powers of 10.                                                                                                                                                                              |
| `iter_rec`               | Number of parameter sets to generate. For example, `hyperparameters_numb=10` and `iter_rec=100` will generate 1,000 models.                                                                                                                 |
| `prior_image`            | Specify the prior image: `"Dirac"` or a `.fits` file containing the model image.                                                                                                                                                            |
| `regularization`         | Choose one or more regularization methods: `['compactness']`, `['hyperbolic']`, or both.                                                                                                                                                    |
| `maxeval`                | Maximum number of evaluations for MiRA to stop if convergence is not achieved (default: 20,000).                                                                                                                                           |

---

### Running the Script
1. Configure the above variables in the script.
2. Run the script:
 ```bash
 python3 OBJECT_INSTRU_MiRA_v3.py
```
## Outputs
Results are stored in subfolders within path_results_tmp, organized by:

- Regularization methods
- Parameter values

Visualizations in the ALL_IMAGES folder include:

- Reconstructed images
- Beam-convolved images
- Comparisons of reconstructed vs. observed visibilities
- Residuals of squared visibilities
