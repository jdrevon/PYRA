# PYRA: Python for MiRA

**PYRA** is a Python wrapper designed to facilitate and optimize the use of **MiRA** (Model-independent Image Reconstruction Algorithm), developed by the [JMMC](https://github.com/emmt/MiRA?tab=readme-ov-file). This tool enables efficient and random scanning of the parameter space for image reconstruction. 

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
   $$
   \frac{λ_{\text{min}}}{6B_{\text{max}}} \leq \text{pixelsize} \leq \frac{λ_{\text{min}}}{2B_{\text{max}}}
   $$

3. **Field of View (FoV)**
   Defined based on the interferometric field of view:
   $$
   \text{FoV} = \frac{λ_{\text{max}}}{B_{\text{min}}}, \quad \text{search range: } [0.7 \times \text{FoV}, 1.0 \times \text{FoV}]
   $$

4. **Tau (\(\tau\))**
   Edge-preserving threshold for hyperbolic regularization:
   $$
   \left(\frac{\text{pixelsize}}{\text{FoV}}\right)^2 \leq \tau \leq \left(\frac{\text{pixelsize}}{\text{FoV}}\right)^2 \times 10^4
   $$

       if we assume an uniform image, the mean pixel value should be equal to flux*(pixelsize/FoV)^2, assuming a nromalized flux, we have then (pixelsize/FoV)^2. So this value is used as the lower boundary of the interval. The higher boundary has been arbitrarly defined as (pixelsize/FoV)^2 * 10^4. 

6. **Gamma (\(\gamma\))**
   FWHM of prior light distribution for compactness regularization:
   $$
   [0.3 \times \text{FoV}, 0.8 \times \text{FoV}]
   $$

7. **Hyperparameter (\(\mu\))**
   Weight for the regularization function in the total chi-square minimization:
   $$
   \chi^2_{\text{TOT}} = \chi^2_{\text{DATA}} + \mu \times f_{\text{regularization}}
   $$


- Broad range for exploration: \( 10^3 \) to \( 10^9 \).

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
