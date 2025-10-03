"""
Python code to perform original CFM-Phyto model (Inomura et al, 2020) with 
parameters from 2020 study
"""

# Import packages
import numpy as np

# GLOBAL VARIABLES  
# Model parameters (fixed constants)
E = 7.74E-1
Q_C_DNA = 9.41E-4
Q_RNA_P_min = 2.23E-4
Q_max_P = 9.79E-3
# Model parameters (MCMC fit for Healey, 1985)
m = 3.93E-1
VI_max = 2.77E2
AI = 8.63E-3
A_Pho = 1.60E1
A_Bio = 2.71E-1
Q_C_Pro_Other = 2.40E-1
A_P_RNA = 4.23E-3
A_Pho_PChl = 2.82E-2
Q_P_Other0 = 6.53E-4
Q_C_Other0 = 1.82E-2
# Elemental ratios of macro-molecules: Inomura et al (2020) Table 1
Y_Chl_NC = 4.0/55.0
Y_Pro_NC = 1.0/4.49
Y_RNA_CP = 10.7/1.0
Y_RNA_NP = 3.8/1.0
Y_DNA_NC = 3.8/11.1
Y_DNA_PC = 1.0/11.1
Y_Plip_CP = 40.0/1.0
Y_Nstore_CN = 2.0/1.0
# Storage (0 because only nutrient limited cases)
Q_N_Sto = 0
Q_P_Sto = 0


def pho(df):
    """Calculate Chl:C ratio"""
    tI = VI_max * (1-np.exp(-AI * df['Irrad.']))
    ChlC = (1+E)/tI * df['Specific.growth'] + m/tI
    return ChlC


def mu_max(df):
    """Calculate maximum growth rate for given environmental conditions"""
    # Photosynthesis
    tI = VI_max * (1-np.exp(-AI * df['Irrad.']))
    A_Chl = (1+E)/tI
    B_Chl = m/tI
	# Coefficients of quadratic
    Q_C_Other = Q_C_Other0 + Q_C_DNA + Q_C_Pro_Other
    a_M = (A_Pho*A_Chl +  A_Bio)*Y_RNA_CP*A_P_RNA
    b_M = (1.0 + A_Pho + Y_Plip_CP*A_Pho_PChl)*A_Chl + A_Bio + (A_Pho*B_Chl + Q_C_Pro_Other)*Y_RNA_CP*A_P_RNA
    c_M = (1.0 + A_Pho + Y_Plip_CP*A_Pho_PChl)*B_Chl + Q_C_Other + Y_RNA_CP*Q_RNA_P_min + Q_N_Sto*Y_Nstore_CN - 1.0
	# Evaluate root of quadratic
    return (- b_M + (b_M*b_M - 4.0*a_M*c_M)**0.5)/(2.0*a_M)


def NC_eval(df):
    """Calculate N:C ratio under N-limitation (Q_N_Sto=0)"""
    # Photosynthesis
    tI = VI_max * (1-np.exp(-AI * df['Irrad.']))
    A_Chl = (1+E)/tI
    B_Chl = m/tI
    # Coefficients of polynomial
    a_N = (A_Pho*A_Chl + A_Bio)*Y_RNA_NP*A_P_RNA
    b_N = Y_Pro_NC*(A_Bio + A_Pho*A_Chl) + (A_Pho*B_Chl + Q_C_Pro_Other)*Y_RNA_NP*A_P_RNA + Y_Chl_NC*A_Chl
    c_N = Y_Pro_NC*(A_Pho*B_Chl + Q_C_Pro_Other) + Y_Chl_NC*B_Chl + Y_RNA_NP*Q_RNA_P_min + Y_DNA_NC*Q_C_DNA + Q_N_Sto
	# Evaluate polynomial
    return a_N*df['Specific.growth']*df['Specific.growth'] + b_N*df['Specific.growth'] + c_N


def PC_eval(df):
    """Calculate P:C ratio under P-limitation (Q_P_Sto=0)"""
    # Photosynthesis
    tI = VI_max * (1-np.exp(-AI * df['Irrad.']))
    A_Chl = (1+E)/tI
    B_Chl = m/tI
    # Coefficients of polynomial
    a_P = (A_Pho*A_Chl + A_Bio)*A_P_RNA
    b_P = (A_Pho*B_Chl + Q_C_Pro_Other)*A_P_RNA + A_Pho_PChl*A_Chl
    c_P = A_Pho_PChl*B_Chl + Q_RNA_P_min + Y_DNA_PC*Q_C_DNA + Q_P_Other0 + Q_P_Sto
	# Evaluate polynomial
    return a_P*df['Specific.growth']*df['Specific.growth'] + b_P*df['Specific.growth'] + c_P