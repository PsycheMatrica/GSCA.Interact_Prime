# InteractGSCA_Prime

### Author:
Heungsun Hwang & Gyeongcheol Cho

## Description:
- The **InteractGSCA_Prime** package enables users to estimate and evaluate GSCA models with component interactions.

## Features:
- Estimate model parameters and calculate their standard errors (SE) along with 95% confidence intervals (CI).
- Enable parallel computing for bootstrap sampling.
- Allow users to determine sign-fixing indicators for components.

## Installation:
To use this package in MATLAB:
1. Clone or download the repository:
   ```bash
   git clone https://github.com/GyeongcheolCho/InteractGSCA_Prime.git
   ```
2. Add the package to your MATLAB path:
   ```matlab
    addpath(genpath('InteractGSCA_Prime'))
   ```

## Usage:
- For examples on how to use the package, refer to the `Run_Example_InteractGSCA.m` file. This file demonstrates the implementation of `InteractGSCA()` using the ACSI dataset.

## Compatibility:
- Tested on MATLAB R2023b.
- Likely compatible with earlier MATLAB versions.

### Citation (APA):
- If you use **InteractGSCA_Prime** in your research or publications, please cite it in APA format as follows:

```plaintext
Hwang, H. & Cho, G. (2024). InteractGSCA_Prime: A package for generalized structured component analysis with component interactions [Computer software]. GitHub. https://github.com/GyeongcheolCho/InteractGSCA_Prime
```
