"""
Apply cell model to compiled stoichiometry database
"""

# Import packages
import pandas as pd
import numpy as np
import scipy as sp
import itertools
from pylr2 import regress2
import CellModel
import Inomura2020_Model


def readParam(FitData):
    """Read in parameters fit using MCMC method"""
    MCMC = pd.read_csv('%s/ParamFit_%s.csv' % (CellModel.RootFolder/CellModel.OutputFolder, FitData), header=0)
    return MCMC['mean']


def LinReg_RMA(df, varx, vary, xseq):
    """Calculate RMA linear regression and r statistic"""
    rma = regress2(df[varx], df[vary], _method_type_2="reduced major axis")
    return rma['intercept']+rma['slope']*xseq, rma['r']


def curve_func(x, a, b, c):
    return a*(x**2) + b*x + c


def QuadReg(df, varx, vary, xseq):
    """Calculate quadratic regression"""
    mask = ~np.isnan(df[varx]) & ~np.isnan(df[vary])
    popt = sp.optimize.curve_fit(curve_func,df[varx][mask],df[vary][mask],
                     bounds=(0,np.inf))[0] # Curve parameters >= 0
    return curve_func(xseq,*popt)


def ModelStoich(df):
    """Calculate model output (N:C and P:C) for several model formulations"""   
    # Modified CFM-Phyto model fit to Healey, 1985 (primary fit for analysis), 
    # nutrient limited data, and nutrient limited data excluding cyanobacteria
    for fit in ['Healey','NutLim','NutLim_NoCyano']:
        param = readParam(fit)
        label = '' if fit == 'Healey' else '_%s'%(fit)
        df['Mod%s_NC'%(label)] = CellModel.NC_eval(df,[*param[0:5],param[7]])
        df['Mod%s_PC'%(label)] = CellModel.PC_eval(df,[*param[0:6],*param[7:9]])

    # CFM-Phyto from Inomura et al, 2020 study with Healey, 1985 parameters
    df['Mod_Inomura_NC'] = Inomura2020_Model.NC_eval(df)
    df['Mod_Inomura_PC'] = Inomura2020_Model.PC_eval(df)

    # Linear regression model for growth rate vs N:C under N limitation
    df_Nlim = df[(df['Manipulation']==2)].dropna(subset=['NC']).reset_index()
    df['Mod_LinRegRMA_NC'] = LinReg_RMA(df_Nlim, 'Specific.growth',
                                        'NC', df['Specific.growth'])[0]
    # Linear regression model using only Healey, 1985 data
    df_Nlim_Healey = df_Nlim[(df_Nlim['Ref.']=='(Healey, 1985)')].reset_index()
    df['Mod_LinRegRMAHealey_NC'] = LinReg_RMA(df_Nlim_Healey, 'Specific.growth',
                                              'NC', df['Specific.growth'])[0]

    # Quadratic regression model for growth rate vs P:C under P limitation
    df_Plim = df[(df['Manipulation']==3)].dropna(subset=['PC']).reset_index()
    df['Mod_QuadReg_PC'] = QuadReg(df_Plim, 'Specific.growth',
                                   'PC', df['Specific.growth'])
    return df


def Model_df():
    """Return dataframe with literature observations and model calculations"""
    df = CellModel.readData() # Read in filtered data
    df = ModelStoich(df) # Calculate model N:C and P:C ratios
    return df


def GrowthRate_Redfield(df,theta):
    """Calculate growth rate when cell N:C equals Redfield ratio (16:106)"""
    # Get coefficients for growth rate to N:C relationship 
    a_N, b_N, c_N = CellModel.NC_eval(df,theta,return_coeffs=True)
    # Set left side of equation to 0, given Redfield N:C = a_N*µ^2 + b_N*µ + c_N
    c_Rf = c_N - 16/106
    # Solve for growth rate
    return (- b_N + (b_N*b_N - 4.0*a_N*c_Rf)**0.5)/(2.0*a_N)


def OptimalGrowth_CN():
    """
    Calculate cellular C:N ratio under optimal growth, (i.e., no nutrient
    or light limitation)
    """
    # Build dataframe with range of conditions for calculations
    T = [10,20,30]
    V = [1,100]
    D = [12,24]
    I = [np.inf] # No photoinhibition, so maximize photosynthesis (equals Pmax)
    #lists = [T, V, D, I]
    #df = pd.DataFrame(list(itertools.product(*lists)), 
    #                  columns=['Temp.', 'V', 'Day.length', 'Irrad.'])
    lists = [V, T, D, I]
    df = pd.DataFrame(list(itertools.product(*lists)), 
                      columns=['V', 'Temp.', 'Day.length', 'Irrad.'])

    # Calculate maximum growth rate and corresponding N:C
    params = readParam('Healey')
    df['mu_max'] = CellModel.mu_max(df,params[0:8])
    df['Specific.growth'] = df['mu_max']
    df['CN'] = 1/CellModel.NC_eval(df,[*params[0:5],params[7]])

    # Calculate growth rate assuming cell N:C equals Redfield
    #df['mu_Redfield'] = GrowthRate_Redfield(df,[*params[0:5],params[7]])
    
    df = df.drop(['Irrad.','Specific.growth'],axis=1)
    return df