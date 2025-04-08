# Structural Dynamics - MDED Time Integration Method

This repository contains MATLAB code and resources related to the paper:

**"A new model-dependent time integration scheme with effective numerical damping for dynamic analysis"**  
*Amir Hossein Namadchi, Emadodin Jandaghi, Javad Alamatian*  
Published in *Engineering with Computers (2021)*  
[DOI: 10.1007/s00366-020-00960-w](https://doi.org/10.1007/s00366-020-00960-w)

## Summary

This project introduces a new semi-explicit time integration method called **MDED** that:
- Inherits properties from the Bathe method  
- Has effective numerical damping  
- Requires no time-step subdivision  
- Reduces computation time in nonlinear analysis  
- Matches the Bathe method in accuracy for linear systems  

## Repository Contents

- `main.m` – Example simulations using the MDED method  
- `functions/` – MATLAB functions for MDED formulation  
- `models/` – Example structural models used in the paper (e.g., circular arch, shallow dome, square plate)  
- `figures/` – Output plots comparing MDED with other methods  

## How to Use

1. Open `main.m` in MATLAB.
2. Choose the example to run (linear or nonlinear).
3. Run the script to simulate and compare results.

## Citing

If you use this code, please cite the original paper:

```
Namadchi AH, Jandaghi E, Alamatian J. (2021). A new model-dependent time integration scheme with effective numerical damping for dynamic analysis. *Engineering with Computers*, 37, 2543–2558.
```
