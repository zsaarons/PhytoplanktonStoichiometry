"""
Python code to read in stoichiometry database and run adapted 
CFM-Phyto model (modified from Inomura et al, 2020)
"""

# Import packages
import numpy as np
import pandas as pd
import pathlib

# GLOBAL VARIABLES  
# File Locations
RootFolder = pathlib.Path(__file__).parent.resolve()
DataFile = 'StoichiometryDatabase.xlsx'
DataSheet = 'Data'
OutputFolder = 'MCMC_Output'
OutputFolder_Trace = 'Trace_NetCDFs'
OutputFolder_Figs = 'MCMC_Figures'
# Assigned model parameters
A_RNA_P = 7.15E-3
Q10_sigma = 2.0
Q_RNA_P_min = 2.23E-4
Q10_kcat_bio = 2.52
Q_max_P = 9.79E-3
Q_C_DNA = 9.41E-4
E = 7.74E-1
Q10_m = 2.5
T_st = 25  #Reference temperature for Q10 calculations
# Elemental ratios of macro-molecules: Inomura et al (2020) Table 1
Y_Chl_NC = 4.0/55.0
Y_Pro_NC = 1.0/4.49
Y_RNA_CP = 10.7/1.0
Y_RNA_NP = 3.8/1.0
Y_DNA_NC = 3.8/11.1
Y_DNA_PC = 1.0/11.1
Y_Plip_CP = 40.0/1.0
Y_Nstore_CN = 2.0/1.0
# Storage (Assume 0 storage for nutrient limited cases)
Q_N_Sto = 0
Q_P_Sto = 0
# Allometric coefficients for model pools with fixed sizes and m
V_ProOther = 0.8 - (172*0.811 + 155*0.939)/327
V_m = -0.189
V_DNA = 0.75 - 0.939
V_OtherC = 0.8 - (103*0.811 + 124*0.939)/227
V_OtherP = 0.76 - (80*0.811 + 104*0.939)/184


def readData():
    """Read Data from excel database and filter for necessary information"""
    df = pd.read_excel(RootFolder/DataFile, sheet_name=DataSheet, header=0)

    # If no cell volume (or volume from reference literature) in database, 
    # calculate volume from cell carbon using Menden-Deuer and Lessard, 2000. 
    # Apply different calculations for diatoms and for all other phyla
    df['V'] = df['Vol.']
    df['V'] = df['V'].fillna(df['Vol.literature'])
    df.loc[(df['Phylum']=='Bacillariophyceae'),'V'] = \
        df.loc[(df['Phylum']=='Bacillariophyceae'),'V'].fillna(4.641*(df.loc[(df['Phylum']=='Bacillariophyceae'),'Qc']*12.011*1E12)**1.233)
    df['V'] = df['V'].fillna(5.114*(df['Qc']*12.011*1E12)**1.065)

    # Remove rows without necessary data
    df = df.dropna(subset=['Temp.','Irrad.','Specific.growth','Day.length','V'])

    # For N,P media concentrations that are not floats (e.g. ~x,<x) set as nan
    df['N'] = pd.to_numeric(df['N'], errors='coerce')
    df['P'] = pd.to_numeric(df['P'], errors='coerce')

    # Calculate observed C:N and C:P ratios
    df['CN'] = 1/df['NC']
    df['CP'] = 1/df['PC']
    return df


def phoForm(df,theta):
    """Calculate rate of photosynthesis"""
    AI,Pm,m = theta #unpacks parameters

    #PI equation
    expo = (AI * df['Irrad.'])/Pm
    tI = Pm * (1-np.exp(-expo)) * df['Day.length']/24
    return tI


def pho(df, theta):
    """Calculate Chl:C ratio"""
    AI,Pm,m = theta #unpacks parameters

    #Maintenance/respiration
    m_temp =  Q10_m**((df['Temp.']-T_st)/10)
    m_vol = df['V']**V_m

    tI = phoForm(df,theta)
    chl = (1+E)/tI * df['Specific.growth'] + (m*m_temp*m_vol)/tI
    return chl


def mu_max(df,theta):
    """Calculate maximum growth rate for given environmental conditions"""
    AI,Pm,m,A_Bio,A_Pho,A_Pho_PChl,Q_C_Other0,Q_C_Pro_Other = theta #unpacks parameters
    #Temp dependent
    Q10_exp = ((df['Temp.']-T_st)/10)
    RNA_T = A_RNA_P/(Q10_sigma**Q10_exp)
    Bio_T = A_Bio/(Q10_kcat_bio**Q10_exp)
    Mtemp = m * (df['V']**V_m) * (Q10_m**Q10_exp)
    #Photosynthesis
    tI = phoForm(df,theta[0:3])
    A_Chl = (1+E)/tI
    B_Chl = Mtemp/tI
    # Volume dependent
    Vdep_ProOther = df['V']**V_ProOther
    Vdep_DNA = df['V']**V_DNA
    Vdep_OtherC = df['V']**V_OtherC
	# coefficients of quadratic
    Q_C_Other = Q_C_Other0*Vdep_OtherC + Q_C_DNA*Vdep_DNA + Q_C_Pro_Other*Vdep_ProOther
    a_M = (A_Pho*A_Chl + Bio_T)*Y_RNA_CP*RNA_T
    b_M = (1.0 + A_Pho + Y_Plip_CP*A_Pho_PChl)*A_Chl + Bio_T + (A_Pho*B_Chl + Q_C_Pro_Other*Vdep_ProOther)*Y_RNA_CP*RNA_T
    c_M = Q_C_Other + (1.0 + A_Pho + Y_Plip_CP*A_Pho_PChl)*B_Chl + Y_RNA_CP*Q_RNA_P_min + Q_N_Sto*Y_Nstore_CN - 1.0
	# evaluate root of quadratic
    return ((- c_M*a_M*4.0 + b_M*b_M)**0.5 - b_M)/(2.0*a_M)


def NC_eval(df,theta,return_coeffs=False):
    """Calculate N:C ratio under N-limitation (Q_N_Sto=0)"""
    AI,Pm,m,A_Bio,A_Pho,Q_C_Pro_Other = theta #unpacks parameters
    #Temp dependent
    Q10_exp = ((df['Temp.']-T_st)/10)
    RNA_T = A_RNA_P/(Q10_sigma**Q10_exp)
    Bio_T = A_Bio/(Q10_kcat_bio**Q10_exp)
    Mtemp = m * (df['V']**V_m) * (Q10_m**Q10_exp)
    #Photosynthesis
    tI = phoForm(df,theta[0:3])
    A_Chl = (1+E)/tI
    B_Chl = Mtemp/tI
    # Volume dependent
    Vdep_ProOther = df['V']**V_ProOther
    Vdep_DNA = df['V']**V_DNA
    # coefficients of polynomial
    a_N = (A_Pho*A_Chl + Bio_T)*Y_RNA_NP*RNA_T
    b_N = Y_Pro_NC*(Bio_T + A_Pho*A_Chl) + (A_Pho*B_Chl + Q_C_Pro_Other*Vdep_ProOther)*Y_RNA_NP*RNA_T + Y_Chl_NC*A_Chl
    c_N = Y_Pro_NC*(A_Pho*B_Chl + Q_C_Pro_Other*Vdep_ProOther) + Y_Chl_NC*B_Chl + Y_RNA_NP*Q_RNA_P_min + Y_DNA_NC*Q_C_DNA*Vdep_DNA + Q_N_Sto
    # Option to return coefficients a_N, b_N, c_N
    if return_coeffs==True: return a_N, b_N, c_N
	# evaluate polynomial
    return a_N*df['Specific.growth']*df['Specific.growth'] + b_N*df['Specific.growth'] + c_N


def PC_eval(df,theta):
    """Calculate P:C ratio under P-limitation (Q_P_Sto=0)"""
    AI,Pm,m,A_Bio,A_Pho,A_Pho_PChl,Q_C_Pro_Other,Q_P_Other0 = theta
    #Temp dependent
    Q10_exp = ((df['Temp.']-T_st)/10)
    RNA_T = A_RNA_P/(Q10_sigma**Q10_exp)
    Bio_T = A_Bio/(Q10_kcat_bio**Q10_exp)
    Mtemp = m * (df['V']**V_m) * (Q10_m**Q10_exp)
    #Photosynthesis
    tI = phoForm(df,theta[0:3])
    A_Chl = (1+E)/tI
    B_Chl = Mtemp/tI
    # Volume dependent
    Vdep_ProOther = df['V']**V_ProOther
    Vdep_DNA = df['V']**V_DNA
    Vdep_OtherP = df['V']**V_OtherP
    # coefficients of polynomial
    a_P = (A_Pho*A_Chl + Bio_T)*RNA_T
    b_P = (A_Pho*B_Chl + Q_C_Pro_Other*Vdep_ProOther)*RNA_T + A_Pho_PChl*A_Chl
    c_P = A_Pho_PChl*B_Chl + Q_RNA_P_min + Y_DNA_PC*Q_C_DNA*Vdep_DNA + Q_P_Other0*Vdep_OtherP + Q_P_Sto
	# evaluate polynomial
    return a_P*df['Specific.growth']*df['Specific.growth'] + b_P*df['Specific.growth'] + c_P


def C_N(df,theta):
    """Calculate carbon biomass in N-limited chemostat"""
    Q_N_NonSto = NC_eval(df,theta)
    N_in_Nlim = np.asarray(df['N']).astype('float32') * 1E-3 #convert N input into chemostat from umol L-1 to mol m-3
    return N_in_Nlim/Q_N_NonSto


def C_P(df,theta):
    """Calculate carbon biomass in P-limited chemostat"""
    Q_P_NonSto = PC_eval(df,theta)
    P_in_Plim = np.asarray(df['P']).astype('float32') * 1E-3 #convert P input into chemostat from umol L-1 to mol m-3
    return P_in_Plim/Q_P_NonSto