# -*- coding: utf-8 -*-
"""
Created on Thu Dec  7 08:01:51 2023

@author: urosz
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import minimize
import CD_model

#Given thermodynamic parameters (dG,dH,dCp) and T-range, this function returns molar residual ellipticity (MRE) as a function of T.
def calculate_cd_signal(dg,dh,dcp,t,n,cd_params_v,matrike_spectro,coef_matrix_polys_total,coef_matrix_polys_double,):
    """
    Calculate the circular dichroism (CD) model using thermodynamic parameters.

    Parameters:
    -----------
    dg : float
        Gibbs free energy at reference temperature (273.15 K) in units kcal/(mol res).
    dh : float
        Enthalpy at reference temperature (273.15 K)  in units kcal/(mol res).
    dcp : float
        Heat capacity at reference temperature (273.15 K) in units kcal/(mol K res).
    t : array
        Array of temperature values in °C.
    n : int
        Number of residues in the peptide.
    cd_params_v : list
        Fixed CD parameters (h2, dh2, k, c, dc) plus nucelation contant (v).
    matrike_spectro : list
        List of spectroscopic contribution matrices [M_H2, M_H1, M_double_H2, M_double_H1].
    coef_matrix_polys_total : array
        polynomial coefficient matrix for total.
    coef_matrix_polys_double : array
        Coefficient matrix for double helices polynomial.

    Returns:
    --------
    cd_model : array
        Modeled CD values as a function of temperature in units 10^3 deg cm^2/(dmol res).
    """

    h2, dh2, k, c, dc, v = cd_params_v  # Unpack fixed parameters
    M_H2, M_H1, M_double_H2, M_double_H1 = matrike_spectro  # Unpack spectroscopic matrices

    # Convert temperatures to Kelvin
    t=np.array(t)
    t_K = t + 273.15
    
    # Calculate temperature-dependent enthalpy
    H_t_w = dh + dcp * (t_K - T0)

    # Temperature dependence of the w parameter
    w = np.exp((-1/R) * (dg / T0 + H_t_w * ((1/t_K) - (1/T0)) + dcp * (1 - (T0 / t_K) - np.log(t_K / T0))))
    
    # Construct the V x W matrix
    VW = CD_model.matrix_vw(v, w, v_pow_max, w_pow_max)
    
    # Apply polynomials to the VW matrix
    matrix_polys_tot = coef_matrix_polys_total * VW
    matrix_polys_2 = coef_matrix_polys_double * VW
    
    # Calculate partition function (q) and normalize for probabilities -> P_matrixs_tot,P_matrixs_2
    q_sums_t=np.sum(matrix_polys_tot,axis=1)
    
    P_matrixs_tot = matrix_polys_tot/q_sums_t[:,None] # each term divided by total sum of polynomial- normalization -> PROBABILITIES
    P_matrixs_2 = matrix_polys_2/q_sums_t[:,None] # each term divided by total sum of polynomial- normalization -> PROBABILITIES
        
    # Compute the CD signal
    cd_model = np.sum((P_matrixs_tot-P_matrixs_2)*( M_H2[None,:]*(h2 + dh2*t)[:,None] +M_H1[None,:]*((h2 + dh2*t)*(1-k/6))[:,None] +(np.full(M_H2.shape, n+1, dtype=int)-(M_H2+M_H1))[None,:] *(c  + dc *t)[:,None])+
               
                (P_matrixs_2)*( M_double_H2[None,:]*(h2 + dh2*t)[:,None] +M_double_H1[None,:]*((h2 + dh2*t)*(1-k/6))[:,None] +(np.full(M_double_H2.shape, n+1, dtype=int)-(M_double_H2+M_double_H1))[None,:] *(c  + dc *t)[:,None]), axis=1)/(n+1)

    
    fh = ((P_matrixs_tot-P_matrixs_2)*(M_H2[None,:] + M_H1[None,:]) +(P_matrixs_2)*(M_double_H2[None,:] + M_double_H1[None,:])).sum(axis=1)/(n+1)
    
    return (cd_model,fh)

def load_cd_data(file_path, n):
    """Load CD data."""
    D_cd = pd.read_csv(file_path, sep='\t').dropna()
    temp_cd = D_cd[f't_{n}']
    mre_cd = D_cd[f'mre_{n}']
    
    return temp_cd, mre_cd

#RETURNS a matrix with all the coefficient of the terms in partiton function and seperatly for only dobule-helix terms.
def load_polynomial_matrices(n, v_pow_max, w_pow_max):
    """Load polynomial matrices."""
    poly_total = open(f'./Q_total/Q_total_{n}.txt').read()
    poly_double = open(f'./Q_double_H/Q_doubleH_{n}.txt').read()
    
    coef_matrix_polys_total = CD_model.polynomial_to_matrix(poly_total, v_pow_max, w_pow_max).transpose().flatten()
    coef_matrix_polys_double = CD_model.polynomial_to_matrix(poly_double, v_pow_max, w_pow_max).transpose().flatten()
    
    return coef_matrix_polys_total, coef_matrix_polys_double

"""

          MAIN PROGRAM
          
"""

# DEFINE CONSTANTS AND PARAMETERS
v_pow_max, w_pow_max = 5, 38  # Max powers of v and w for the VxW matrix
T0 = 273.15  # Reference temperature in Kelvin
R = 1.987 * 10**(-3)  # Gas constant [kcal K-1 mol-1]

# FIXED CD PARAMETERS AND NUCLEATION CONSTANT
# These are baseline CD parameters and the nucleation constant from Zavrtanik et al. (2024)
cd_params_v_fix = [-41.0, 0.1, 3.4, 2.1, -0.045, 0.07]

# INPUT DATA
n = 32  # Number of peptide residues

# Load experimental CD data: temperature (in °C) and MRE values (in 10^3 deg cm^2/(dmol res))
temp_cd, mre = load_cd_data(f'test_{n}.txt', n)

# LOAD REQUIRED MATRICES FOR CD SIGNAL CALCULATION
# Load polynomial matrices for partition functions (Q_total and Q_doubleH)
coef_matrix_polys_total, coef_matrix_polys_double = load_polynomial_matrices(n, v_pow_max, w_pow_max)

# Load spectroscopic contribution matrices
M_H2, M_H1, M_double_H2, M_double_H1 = CD_model.spectro_matrices(v_pow_max, w_pow_max)
matrike_spectro = [M_H2, M_H1, M_double_H2, M_double_H1]

# OBJECTIVE FUNCTION FOR OPTIMIZATION
def fit_objective(params, params_fix, tcd, n, data_cd):
    """
    Objective function to minimize during fitting.

    Parameters:
    -----------
    params : list
        List of fitting parameters [dg, dh, dcp].
    params_fix : list
        List of fixed CD baseline parameters and nucleation constant.
    tcd : numpy array
        Array of temperature data (°C).
    n : int
        Number of peptide residues.
    data_cd : numpy array
        Experimental CD data (MRE) in units 10^3 deg cm^2/(dmol res).

    Returns:
    --------
    residuals : numpy array
        The residuals between the model and the experimental CD data.
    """
    dg, dh, dcp = params
    
    model_cd = calculate_cd_signal(dg, dh, dcp, tcd, n, params_fix, matrike_spectro, coef_matrix_polys_total, coef_matrix_polys_double)[0]
    
    residuals = model_cd - data_cd
    return np.sum(residuals ** 2)

# PERFORM FITTING TO OPTIMIZE PARAMETERS
initial_guess = [-0.22, -1.3, 0.002]  # Initial guess for [dg, dh, dcp]

result = minimize(fit_objective, initial_guess, args=(cd_params_v_fix, temp_cd, n, mre), method='Nelder-Mead')

# Extract optimized parameters from the result
dg_opt, dh_opt, dcp_opt = result.x

# PRINT RESULTS - optimized TD parameters
print(f'Optimized parameters:\n dG = {dg_opt:.3f} kcal/(mol res)\n dH = {dh_opt:.2f} kcal/(mol res)\n dCp = {dcp_opt:.3f} kcal/(mol K res)')

# CALCULATE THE FITTED MODEL
fit_model_cd = calculate_cd_signal(dg_opt, dh_opt, dcp_opt, temp_cd, n, cd_params_v_fix, matrike_spectro, coef_matrix_polys_total, coef_matrix_polys_double)


# PLOT RESULTS
fig, ax1 = plt.subplots()

# Plot MRE vs Temperature (left y-axis)
color = 'tab:blue'
ax1.set_xlabel('Temperature (°C)')
ax1.set_ylabel('MRE (10^3 deg cm^2/dmol res)', color=color)
ax1.plot(temp_cd, fit_model_cd[0], label='Fitted Model', color=color)
ax1.scatter(temp_cd, mre, label='Experimental Data', color='red')
ax1.tick_params(axis='y', labelcolor=color)
ax1.legend(loc="upper left")

# Create a secondary y-axis to plot fH vs Temperature
ax2 = ax1.twinx()
color = 'black'  # Change this to 'black' instead of 'tab:black'
ax2.set_ylabel('fH', color=color)
ax2.plot(temp_cd, fit_model_cd[1], color=color)
ax2.tick_params(axis='y', labelcolor=color)

# Show plot
plt.title('MRE and fH vs Temperature')
plt.show()



