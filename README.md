# Spherical-Spline-Interpolation

## Overview

This project implements **spherical spline interpolation** to address approximation problems on the sphere, specifically for **gravity potential data**. The interpolation uses three distinct kernel functions—**Abel-Poisson**, **Singularity**, and **Logarithmic**—to interpolate data from a coarser 6-degree grid to a finer 3-degree grid. The work is particularly relevant to fields such as **geophysics**, **physical geodesy**, and **environmental sciences**, where high-resolution, spatially distributed data are critical.

The study also explores the effects of **ill-conditioning** in the design matrix and applies **regularization methods** (TSVD, Tikhonov, Cholesky decomposition) to achieve stable solutions. Finally, the implementation evaluates the interpolation performance under **noise conditions**.

---

### Kernel Functions
Three kernel functions are used:
1. **Abel-Poisson Kernel**:

2. **Singularity Kernel**:

3. **Logarithmic Kernel**:

### Regularization Methods
- **Cholesky Decomposition**: For numerically stable solutions in well-conditioned cases.
- **TSVD (Truncated Singular Value Decomposition)**: For handling ill-conditioned matrices by truncating small singular values.
- **Tikhonov Regularization**: Uses a regularization parameter derived from **Variance Component Estimation (VCE)**.

### Dataset
The gravity potential data used in this project was obtained from the **International Center for Global Earth Models (ICGEM)**. The original dataset is sampled on a 6-degree global grid and interpolated to a finer 3-degree grid.

## Results

### Without Noise
The interpolation was performed using all three kernel functions. Below is a summary of the results:

#### Abel-Poisson Kernel
| Method              | Mean Error (\(m^2/s^2\)) | Norm of Errors (\(m^2/s^2\)) |
|---------------------|---------------------------|-------------------------------|
| Cholesky            | 0.1612                   | 2014.1234                     |
| TSVD                | 0.1766                   | 4175.8890                     |
| Tikhonov (VCE)      | -2.4909                  | 3931.1049                     |

#### Singularity Kernel
| Method              | Mean Error (\(m^2/s^2\)) | Norm of Errors (\(m^2/s^2\)) |
|---------------------|---------------------------|-------------------------------|
| Cholesky            | 0.1814                   | 2014.3146                     |
| TSVD                | 0.1696                   | 2016.4105                     |
| Tikhonov (VCE)      | 2.0778                   | 3931.1049                     |

#### Logarithmic Kernel
| Method              | Mean Error (\(m^2/s^2\)) | Norm of Errors (\(m^2/s^2\)) |
|---------------------|---------------------------|-------------------------------|
| Cholesky            | 0.1503                   | 2005.1889                     |
| TSVD                | 0.1746                   | 2013.2417                     |
| Tikhonov (VCE)      | 0.4764                   | 5436.0035                     |

### With Added Noise
White noise with a standard deviation of \(200 \, \text{m}^2/\text{s}^2\) was added to the input data. Below are the results for the **Abel-Poisson Kernel**:

| Method              | Mean Error (\(m^2/s^2\)) | Norm of Errors (\(m^2/s^2\)) |
|---------------------|---------------------------|-------------------------------|
| Cholesky            | -3.1558                  | 5548.1711                     |
| TSVD                | -3.1328                  | 12137.6909                    |
| Tikhonov (VCE)      | 71.1711                  | 8200.7749                     |

Tikhonov regularization provided smoother and more realistic results in noisy conditions compared to Cholesky and TSVD.


