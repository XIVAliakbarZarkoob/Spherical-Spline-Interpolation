# Spherical-Spline-Interpolation

## Overview

This project implements **spherical spline interpolation** to address approximation problems on the sphere, specifically for **gravity potential data**. The interpolation uses three distinct kernel functions—**Abel-Poisson**, **Singularity**, and **Logarithmic**—to interpolate data from a coarser 6-degree grid to a finer 3-degree grid. The work is particularly relevant to fields such as **geophysics**, **physical geodesy**, and **environmental sciences**, where high-resolution, spatially distributed data are critical.

The study also explores the effects of **ill-conditioning** in the design matrix and applies **regularization methods** (TSVD, Tikhonov, Cholesky decomposition) to achieve stable solutions. Finally, the implementation evaluates the interpolation performance under **noise conditions**.

---

## Key Features

### Spherical Spline Interpolation
The interpolation follows the equation:

\[ S(y) = \sum_{i=1}^{N} a_i K(x_i, y) \]

Where:
- \( y \): Desired point.
- \( x_i \): Input points.
- \( a_i \): Coefficients for spherical spline interpolation.
- \( K(x_i, y) \): Kernel function.

The coefficients \( a_i \) are computed by solving a linear system:

\[ A a = F \]

Here, \( A \) is the design matrix formed by the kernel functions, and \( F \) is the vector of function values at the input points. Given the **ill-conditioning** of \( A \), regularization methods are used to compute \( a_i \).

### Kernel Functions
Three kernel functions are used:
1. **Abel-Poisson Kernel**:
   \[
   K_{Abel-Poisson}(x, y) = \frac{1}{4\pi} \frac{1-h^2}{(L_h(x, y))^{3/2}}
   \]

2. **Singularity Kernel**:
   \[
   K_{Singularity}(x, y) = \frac{1}{2\pi} \frac{1}{(L_h(x, y))^{1/2}}
   \]

3. **Logarithmic Kernel**:
   \[
   K_{Logarithmic}(x, y) = \frac{1}{2\pi h} \ln \left( 1 + \frac{2h}{(L_h(x, y))^{1/2} + 1 - h} \right)
   \]

Where:
- \( L_h(x, y) = 1 + h^2 - 2h \langle x, y \rangle \).
- \( h \): Smoothing parameter (\( 0 < h < 1 \)).

### Regularization Methods
- **Cholesky Decomposition**: For numerically stable solutions in well-conditioned cases.
- **TSVD (Truncated Singular Value Decomposition)**: For handling ill-conditioned matrices by truncating small singular values.
- **Tikhonov Regularization**: Uses a regularization parameter derived from **Variance Component Estimation (VCE)**.

---

## Implementation Details

### Dataset
The gravity potential data used in this project was obtained from the **International Center for Global Earth Models (ICGEM)**. The original dataset is sampled on a 6-degree global grid and interpolated to a finer 3-degree grid.

### Error Definition
The interpolation error is defined as:

\[ e = l - \hat{l} \]

Where:
- \( l \): True value (obtained from the XGM2019 model).
- \( \hat{l} \): Interpolated value.

### Handling Ill-Conditioned Design Matrix
The design matrix \( A \) is ill-conditioned with a condition number of \( 7.15 \times 10^{49} \). Without regularization, solutions are unstable, as demonstrated by the **Picard plot** and interpolation results.

---

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


