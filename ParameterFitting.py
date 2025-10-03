"""
Fit model parameters for temperature and volume dependent adaptation of
CFM-Phyto using Bayesian Markov Chain Monte Carlo (MCMC) approach
"""

# Import packages
import pymc as pm
import arviz as az
import numpy as np
import pandas as pd
from scipy.optimize import least_squares
import CellModel
from pathlib import Path

#Random number generator for least squares initial values
rng = np.random.default_rng(1234)


def Output_Folders():
    """Check if necessary output folders exist and build directories if not"""
    for dir in ['',CellModel.OutputFolder_Figs,CellModel.OutputFolder_Trace]:
        Path(CellModel.RootFolder / CellModel.OutputFolder / 
             dir).mkdir(parents=True, exist_ok=True)
    return None


def Filtered_data(df):    
    """Get subsets of data needed for fitting parameters"""
    df_Replete = df.loc[(df['Manipulation'] == 1)]
    df_NLim = df.loc[(df['Manipulation']==2)].dropna(subset=['NC']) #for N limited
    df_PLim = df.loc[(df['Manipulation']==3)].dropna(subset=['PC']) #for P limited

    # Chemostat required for calculating stoichiometry from culture
    # carbon and nutrient input
    df_NLim_Ccell = df.loc[((df['Manipulation']==2)&(df['Growth.conditions']==3))].dropna(subset=['NC','N','Ccell']) #N limited
    df_PLim_Ccell = df.loc[((df['Manipulation']==3)&(df['Growth.conditions']==3))].dropna(subset=['PC','P','Ccell']) #P limited
    return df_Replete, df_NLim, df_PLim, df_NLim_Ccell, df_PLim_Ccell


def Residual_Chl(theta,df):
    """Residual for least squares for Chl:C calculation"""
    residual = df['ChlC'] - CellModel.pho(df,theta)
    return residual


def Residual_CNP(theta,Chl_params,df_Replete,df_NLim,df_PLim):
    """Residual for least squares for C:N:P calculations"""
    NC = df_NLim['NC'] - CellModel.NC_eval(df_NLim,[*Chl_params,*theta[0:2],theta[4]])
    PC = df_PLim['PC'] - CellModel.PC_eval(df_PLim,[*Chl_params,*theta[0:3],*theta[4:6]])
    umax = df_Replete['Specific.growth'] - CellModel.mu_max(df_Replete,[*Chl_params,*theta[0:5]])
    return np.concatenate((NC,PC,umax))


def LeastSquares_Fit(MCMC_step,fun_input):
    """Run least squares fit to get initial estimate for parameters"""
    # Provide initial values/bounds/input and call relevant function to calculate residuals
    if MCMC_step == 'Chl':   # 1st step of MCMC: photosynthesis/respiration parameters
        results = least_squares(Residual_Chl, x0=[rng.random()]*3, bounds=(1E-6,1E3), args=(fun_input))

    elif MCMC_step == 'CNP':     # 2nd step of MCMC: C:N:P stoichiometry parameters
        results = least_squares(Residual_CNP, x0=[rng.random()]*6, bounds=(1E-6,1E3), args=(fun_input))

    # Print and return solution and cost function
    theta = results.x; print("Least Squares: ", theta)
    cost = results.cost; print("Error: ", cost)
    return theta, cost


def MCMC_Output(MCMC_step,FitData,trace,rows,cols):
    """Save/return output of MCMC parameter fitting"""
    # Save to file
    trace.to_netcdf('%s/%s_Trace_%s.nc' % (CellModel.RootFolder/
                                           CellModel.OutputFolder/
                                           CellModel.OutputFolder_Trace,
                                           MCMC_step,FitData))

    # Plot posterior and save
    axes = az.plot_posterior(trace, grid=(rows,cols), figsize=(cols*5,rows*5))
    fig = axes.ravel()[0].figure
    fig.savefig('%s/%s_PosteriorPlot_%s' % (CellModel.RootFolder/
                                            CellModel.OutputFolder/
                                            CellModel.OutputFolder_Figs,
                                            MCMC_step,FitData))

    # Print and return output of MCMC fit
    params = az.summary(trace, round_to=10); print(params) 
    return params


def Run_Chl_Model(df,FitData):
    """Run MCMC model to fit photosynthesis and maintenance parameters"""
    # Remove rows that do not contain any of N:C, P:C, or N:P
    df = df.dropna(subset=['NC','PC','NP'],how='all')

    df_Chl = df.dropna(subset=['ChlC']) # Filter data (must include observed Chl:C values)

    # Least squares solution for Chl:C which will be used to inform priors for PyMC model
    theta, cost = LeastSquares_Fit('Chl',[df_Chl]) #calculate lsq values and error
    
    #Set up model and run MCMC fitting using NUTS sampler
    with pm.Model() as Model_Chl:
        AI = pm.TruncatedNormal("AI", mu=theta[0], sigma=10*theta[0], lower=0)#, upper=10) #upper = 100 (experiment for C:N:P output)
        Pm = pm.TruncatedNormal("Pm", mu=theta[1], sigma=10*theta[1], lower=0)
        m = pm.TruncatedNormal("m", mu=theta[2], sigma=10*theta[2], lower=0, upper=0.9)
        pm.Normal("Chl", mu=CellModel.pho(df_Chl,[AI,Pm,m]), sigma=cost, observed=list(df_Chl['ChlC']))
        trace_Chl = pm.sample_prior_predictive()

    with Model_Chl:
        step = pm.NUTS(target_accept=0.9) 
        # For tunes and draws, choose values that lead to rhat=1
        trace_Chl.extend(pm.sample(step=step,chains=4,tune=10000,draws=5000,cores=1)) 
    
    # Process and save output
    Chl_Params = MCMC_Output('ChlC',FitData,trace_Chl,1,3)
    return Chl_Params  


def Run_CNP_Model(df, Chl_params, FitData):
    """Run MCMC model to fit remaining parameters"""
    df_Replete, df_NLim, df_PLim, df_NLim_Ccell, df_PLim_Ccell = Filtered_data(df)
    
    # PyMC Model for remaining params given output from chl:C and using least squares to inform priors
    theta, cost = LeastSquares_Fit('CNP',[Chl_params,df_Replete,df_NLim,df_PLim]) 

    with pm.Model() as Model_CNP:
        A_Bio = pm.TruncatedNormal("A_Bio", mu=theta[0], sigma=1E6*theta[0], lower=0)
        A_Pho = pm.TruncatedNormal("A_Pho", mu=theta[1], sigma=1E6*theta[1], lower=0)
        A_Pho_PChl = pm.TruncatedNormal("A_Pho_PChl", mu=theta[2], sigma=1E6*theta[2], lower=0)
        Q_C_Other0 = pm.TruncatedNormal("Q_C_Other0", mu=theta[3], sigma=1E6*theta[3], lower=0)
        Q_C_Pro_Other = pm.TruncatedNormal("Q_C_Pro_Other", mu=theta[4], sigma=1E6*theta[4], lower=0)
        Q_P_Other0 = pm.TruncatedNormal("Q_P_Other0", mu=theta[5], sigma=1E6*theta[5], lower=0)

        theta = [A_Bio,A_Pho,A_Pho_PChl,Q_C_Other0,Q_C_Pro_Other,Q_P_Other0]

        µmax = CellModel.mu_max(df_Replete,[*Chl_params,*theta[0:5]])
        pm.Normal("µ_max", mu=µmax, sigma=cost, observed=df_Replete['Specific.growth'])

        NC_N = CellModel.NC_eval(df_NLim,[*Chl_params,*theta[0:2],theta[4]])
        pm.Normal("NC_Nlim", mu=NC_N, sigma=cost, observed=df_NLim['NC'])
        PC_P = CellModel.PC_eval(df_PLim,[*Chl_params,*theta[0:3],*theta[4:6]])
        pm.Normal("PC_Plim", mu=PC_P, sigma=cost, observed=df_PLim['PC'])

        Ccell_N = CellModel.C_N(df_NLim_Ccell,[*Chl_params,*theta[0:2],theta[4]])
        pm.Normal("C_Cell_Nlim", mu=Ccell_N, sigma=cost, observed=df_NLim_Ccell['Ccell'])
        Ccell_P = CellModel.C_P(df_PLim_Ccell,[*Chl_params,*theta[0:3],*theta[4:6]])
        pm.Normal("C_Cell_Plim", mu=Ccell_P, sigma=cost, observed=df_PLim_Ccell['Ccell'])
    
        trace_CNP = pm.sample_prior_predictive()
    
    with Model_CNP:
        step = pm.NUTS(target_accept=0.9) 
        # For tunes and draws, choose values that lead to rhat=1
        trace_CNP.extend(pm.sample(step=step,chains=4,tune=10000,draws=5000,cores=1)) 

    # Process and save output
    CNP_Params = MCMC_Output('CNP',FitData,trace_CNP,2,4)
    return CNP_Params


def main():
    """Fit parameters and save"""
    # Ensure folders exist to save output
    Output_Folders()
    
    # Select which data will be used to fit model
    # Note "Q_C_Other0" param is not calculated for nutrient limited options
    FittingDatasets = ['Healey', 'NutLim', 'NutLim_NoCyano']
    for FittingData in FittingDatasets:
        df = CellModel.readData()
        # Filter to only include data used to fit model with options:
        if FittingData == 'Healey':
                df = df[(df['Ref.'] == '(Healey, 1985)')]
        elif FittingData == 'NutLim':
            df = df[(df['Manipulation'] > 1)]
            df = df.drop(df.loc[df['Ref.'].isin(['(Marañón et al, 2018)'])].index)
        elif FittingData == 'NutLim_NoCyano':
            df = df[(df['Manipulation'] > 1) & (df['Class'] != 'Cyanophyceae')]
            df = df.drop(df.loc[df['Ref.'].isin(['(Marañón et al, 2018)'])].index)

        # Fit model parameters in two steps using MCMC method
        # Fit photosynthesis and maintenance parameters and return output
        MCMCFit_Chl = Run_Chl_Model(df, FittingData)
        # Fit rest of parameters and return output
        MCMCFit_CNP = Run_CNP_Model(df, MCMCFit_Chl['mean'], FittingData)
        # Concatenate into single dataframe
        MCMCFit = pd.concat([MCMCFit_Chl, MCMCFit_CNP])

        # Save parameter fits
        MCMCFit.to_csv('%s/ParamFit_%s.csv' % (CellModel.RootFolder/
                                               CellModel.OutputFolder,
                                               FittingData))
    return None

main()