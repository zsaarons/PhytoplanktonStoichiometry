"""
Produce figures for main text and supporting information
"""

# Import packages
import numpy as np
import matplotlib.pyplot as plt
import Plotting
import ApplyModel


####################################
# Main Text
####################################

def Fig1_ObsStoich(df,pdf):
    """
    Produce Figure 1: Histograms of observed stoichiometry (N:C, P:C, N:P)
    across nutrient manipulations (N-limited, P-limited, nutrient replete)
    """
    # Set up figure
    Plotting.figset()
    fig, axes = plt.subplots(3,3,layout='constrained',figsize=(5.5,4.5))
    Plotting.LabelSubplot(axes, xcoord=-.13) # label axes

    # Plot histograms
    Plotting.Histogram_ObservedCNP(axes, df, ['NC','PC','NP'])
    for i in range(len(axes)): # Set y-axes to be uniform within each row
        for j in range(len(axes[0])):
            xmax = 2 if j<2 else 1
            ymax = 2-j
            axes[i,j].set_ylim(bottom=0,top=axes[ymax,j].get_ylim()[1])
            axes[i,j].set_xlim(left=0,right=axes[xmax,j].get_xlim()[1])

    # Create legend
    Plotting.LegendContents(fig, ['Redfield line', 'Mean ratio'], (0.9,0.98))
    if pdf==True: Plotting.saveFig(fig,'Fig1_HistCNP')
    return None


def Figs_2_5_ObsCNP_Irrad(df,ratio,pdf):
    """
    Produce Figures 2 and 5: Growth rate vs N:C and P:C under N and P
    limitation respectively, colored by irradiance
    """
    # Set up figure
    Plotting.figset()
    fig = plt.figure(layout='constrained',
                     figsize=(3.2 if ratio=='NC' else 3.3, 2.5))
    ax1=plt.gca()

    # Plot ratio under relevant limitation
    lim = 2 if ratio=='NC' else 3
    cbar_lim = 9.9 if ratio=='PC' else None
    df_lim = df[(df['Manipulation']==lim)].dropna(subset=[ratio]).reset_index()
    Plotting.Scatterplot_GrowthVsStoich(ax1,df_lim,ratio,'I_norm',
                                        'Linear' if ratio=='NC' else 'Quadratic',
                                        cbar_min=cbar_lim)

    # Create legend
    colorbar_labels = [10,100,1000] if ratio=='NC' else [10,100,600]
    Plotting.Add_Colorbar(fig,ax1,df_lim,colorbar_labels,cbar_min=cbar_lim) 
    regress_leg = 'Linear regression' if ratio=='NC' else 'Quadratic regression'
    Plotting.LegendContents(fig, ['Healey, 1985 data',regress_leg,'Redfield line'],
                            (0.97 if ratio=='NC' else 0.99,0.33))
    if pdf==True: Plotting.saveFig(fig,f'Fig{'2' if ratio=='NC' else '5'}_{ratio}_Irrad')
    return None


def Figs_3_6_RelevantLim(df,ratio,pdf):
    """
    Produce Figures 3 and 8: Observed vs model N:C and P:C under N and P
    limitation respectively,
    """
    # Set up figure
    Plotting.figset()
    fig = plt.figure(layout='constrained', figsize=(3.5, 2.25)) 
    ax1=plt.gca()

    # Plot ratio under relevant limitation
    lim = 2 if ratio=='NC' else 3
    df_lim = df[(df['Manipulation']==lim)].dropna(subset=[ratio]).reset_index()
    Plotting.Scatterplot_ObsVsModelStoich(ax1,df_lim,ratio,'Phylum_Plotting',
                                          (0.1,0.88 if ratio=='NC' else 0.68))

    # Create legend
    legend_items = ['Phylum:', *sorted(df_lim['Phylum_Plotting'].unique()), '',
                    'Healey, 1985 data', 'Linear regression', '1:1 line',
                    'Redfield star']
    Plotting.LegendContents(ax1, legend_items, (1.5,1.07))
    if pdf==True: Plotting.saveFig(fig,'Fig3_NC_Nlim' if ratio=='NC' else 'Fig6_PC_Plim')
    return None



def Fig4_NC_Replete(df,pdf):
    """
    Produce Figure 6: N:C under nutrient replete conditions, where panel a is
    growth rate vs N:C and panels b and c show macromolecular composition under
    light and non-light limited growth respectively
    """
    # Set up figure
    Plotting.figset()
    fig, axes = plt.subplot_mosaic([['obs','light'], ['obs','nonlight']],
                                   layout='constrained',width_ratios=[1.8,1],
                                   height_ratios=[1.5,1],figsize=(5,3))
    fig.set_constrained_layout_pads(h_pad=0.02)
    Plotting.LabelSubplot(np.array(list(axes.values())),
                          xcoord=[-0.18,-0.12,-0.12],
                          ycoord = [0.97, 0.96, 0.96]) # label axes

    # Get replete data (Fernandez-Reiriz et al, 1989 excluded due to outliers)
    df_Rep = df[(df['Manipulation']==1)].drop(df.loc[df['Ref.'].isin(['(Fernández-Reiriz et al, 1989)'])].index).dropna(subset=['NC']).reset_index()

    # Plot growth rate vs N:C and annotate
    Plotting.Scatterplot_GrowthVsStoich(axes['obs'],df_Rep,'NC','Phylum_Plotting')
    Plotting.Annotate_NC_Replete(axes['obs'])

    # Plot N allocation schematics for light and non-light limited growth
    Plotting.Allocation_Schematic(axes['light'], None, 0,
                                  [[0.03,0.03], [0,0.055], [0.2011,0.0287]],
                                  [-0.0587,0.2311], [0.5,0.12], [1.5,0.05])
    Plotting.Allocation_Schematic(axes['nonlight'], 'Growth Rate', 2,
                                  [[0.03,0.03], [0,0.055], [0.0010,0.071]],
                                  [0.0625,0.031], [1.2,0.09], [1.5,0.05])

    # Create legend
    legend_items = ['Phylum:', *sorted(df_Rep['Phylum_Plotting'].unique()), '',
                    'Redfield line']
    Plotting.LegendContents(fig, legend_items, (0.53,0), cols=4)
    if pdf==True: Plotting.saveFig(fig,'Fig4_NC_Replete')
    return None


####################################
# Supporting Information
####################################

def FigS1_Healey(df,pdf):
    """
    Produce Figure S1: Growth rate vs N:C and P:C under N and P limitation
    respectively for Healey, 1985 dataset, with CFM-Phyto fit for given ratio
    and maximum growth rate
    """
    # Set up figure
    Plotting.figset()
    fig, axes = plt.subplots(1,2,layout='constrained',figsize=(6.5,2.5))
    Plotting.LabelSubplot(axes) # label axes

    # Get Healey data
    dfHealey = df[(df['Ref.']=='(Healey, 1985)')] #.assign(Marker=1) #REMOVE THIS
    dfHealey_Nlim = dfHealey[(dfHealey['Manipulation']==2)].dropna(subset=['NC']).reset_index()
    dfHealey_Plim = dfHealey[(dfHealey['Manipulation']==3)].dropna(subset=['PC']).reset_index()

    # Plot ratio under relevant limitation
    Plotting.Scatterplot_GrowthVsStoich(axes[0],dfHealey_Nlim,'NC','Irrad.',Rf_ratio=False)
    Plotting.Scatterplot_GrowthVsStoich(axes[1],dfHealey_Plim,'PC','Irrad.',Rf_ratio=False)

    # Add CFM-Phyto predictions
    Plotting.Healey_ModelFit(axes[0],dfHealey_Nlim,'NC')
    Plotting.Healey_ModelFit(axes[1],dfHealey_Plim,'PC')

    # Create legend
    h,l = axes[1].get_legend_handles_labels()
    [ax.legend().set_visible(False) for ax in axes]
    Plotting.LegendContents(axes[1], [r'Irradiance ($\mathrm{µE \; m^{-2} \; s^{-1}}$):',
                                      *l,'','Healey, 1985 data','Model stoichiometry',
                                      'Model $\\mu_{max}$','Redfield star'],
                                      (1.45,1.05), new_l=l, new_h=h)
    if pdf==True: Plotting.saveFig(fig,'FigS1_HealeyCNP')
    return None


def FigS2_ObsStoich_Reverse(df,pdf):
    """
    Produce Figure S2: Histograms of observed stoichiometry (C:N, C:P, P:N)
    across nutrient manipulations (N limited, P limited, nutrient replete)
    """
    # Calculate P:N ratio
    df['PN'] = 1/df['NP']
    # Set up figure
    Plotting.figset()
    fig, axes = plt.subplots(3,3,layout='constrained',figsize=(5.5,4.5))
    Plotting.LabelSubplot(axes, xcoord=-.13) # label axes

    # Plot histograms
    Plotting.Histogram_ObservedCNP(axes, df, ['CN','CP','PN'])
    for i in range(len(axes)): # Set y-axes to be uniform within each row
        for j in range(len(axes[0])):
            xmax = j if j<2 else 0
            ymax = 2 if j<2 else 1
            axes[i,j].set_ylim(bottom=0,top=axes[ymax,j].get_ylim()[1])
            axes[i,j].set_xlim(left=0,right=axes[xmax,j].get_xlim()[1])

    # Create legend
    Plotting.LegendContents(fig, ['Redfield line', 'Mean ratio'], (0.9,0.98))
    if pdf==True: Plotting.saveFig(fig,'FigS2_HistCNP_Reverse')
    return None


def FigS3_S8_CultureType(df,ratio,pdf):
    """
    Produce Figures S3 and S8: a) Growth rate vs N:C (P:C) under N (P) limitation
    and b) observed vs model N:C (P:C) under N (P) limitation, colored by culture
    type (batch, semi-continuous, chemostat)
    """
    # Set up figure
    Plotting.figset()
    fig, axes = plt.subplots(1,2,layout='constrained',figsize=(6.5,2.5))
    Plotting.LabelSubplot(axes) # label axes

    # Plot ratio under relevant limitation
    lim = 2 if ratio=='NC' else 3
    cbar_lim = 9.9 if ratio=='PC' else None
    df_lim = df[(df['Manipulation']==lim)].dropna(subset=[ratio]).reset_index()

    Plotting.Scatterplot_GrowthVsStoich(axes[0],df_lim,ratio,'Growth.conditions',
                                        'Linear' if ratio=='NC' else 'Quadratic',
                                        cbar_min=cbar_lim)
    Plotting.Scatterplot_ObsVsModelStoich(axes[1],df_lim,ratio,'Growth.conditions',
                                          (0.1,0.88 if ratio=='NC' else 0.68))
    
    # Create legend
    legend_items = ['Culture system:','Batch','Semi-continuous','Chemostat', '',
                    'Healey, 1985 data', 'Linear regression' if ratio=='NC' else 'Regression',
                    '1:1 line', 'Redfield star']
    Plotting.LegendContents(axes[1], legend_items, (1.45,1.05))
    if pdf==True: Plotting.saveFig(fig,'FigS3_NC_Nlim' if ratio=='NC' else 'FigS8_PC_Plim')
    return None


def FigS4_NC_Nlim_Phylum(df,pdf):
    """
    Produce Figure S4: Observed vs model N:C under N limitation, seperated by
    phylum
    """
    # Set up figure
    Plotting.figset()
    fig, axes = plt.subplots(2,4,layout='constrained',figsize=(7.3,4))
    fig.set_constrained_layout_pads(h_pad=0.01)
    axes[-1,-1].axis('off')
    
    #Filter info for subpanels
    df_Nlim = df[(df['Manipulation']==2)].dropna(subset=['NC']).reset_index()
    phyla = df_Nlim['Phylum_Plotting'].sort_values().unique()
    Plotting.LabelSubplot(axes,titles=[*phyla,''],ycoord=1.05)
    axes[-1,-1].texts[-1].set_text('')

    # Plot N:C across phyla/subpanels
    maxAx = np.max([df_Nlim['NC'].max(), df_Nlim['Mod_NC'].max()])
    for i, ax in enumerate(axes.flat):
        if i<7: Plotting.Scatterplot_ObsVsModelStoich(ax,df_Nlim[(df_Nlim['Phylum_Plotting']==phyla[i])],
                                                      'NC','Phylum_Plotting',(0.1,0.85),r2_loc=(0.65,0.25),
                                                      axes_max=maxAx)
        if i%4>0: ax.set(ylabel=None)
        ax.set(xlabel=None if i<3 else 'Obs. '+r'$(\mathrm{mol \; N}) \: (\mathrm{mol \; C})^{-1}$')

    legend_items = ['Healey, 1985 data', 'Linear regression', '1:1 line', 'Redfield star']
    Plotting.LegendContents(fig, legend_items, (0.9,0.4))
    if pdf==True: Plotting.saveFig(fig,'FigS4_NC_Nlim_Phylum')
    return None


def Figs_S5_S9_ModelOpts(df,ratio,pdf):
    """
    Produce Figures S5 and S9: Observed vs model N:C and P:C under N and P
    limitation respectively for 2 model formulations
    """
    # Dictionary containing model options
    ModOpts = {'NC': {'titles': ['Inomura et al., 2020','Nutrient-limited fit'],
                      'models': ['Inomura','NutLim']},
               'PC': {'titles': ['Nutrient-limited fit','No Cyanobacteria fit'],
                      'models': ['NutLim','NutLim_NoCyano']}}

    # Set up figure
    Plotting.figset()
    fig, axes = plt.subplots(1,2,layout='constrained',figsize=(6,2.5))
    Plotting.LabelSubplot(axes,titles=ModOpts[ratio]['titles'], ycoord=1.07)

    rmse = [(0.55,0.25),(0.1,0.8)] if ratio=='NC' else [(0.15,0.75),(0.5,0.29)]
    r2 = [(0.42,0.91),(0.7,0.25)] if ratio=='NC' else [(0.7,0.29),(0.3,0.75)]

    # Plot fit for each model formulation
    for i in range(len(ModOpts[ratio]['models'])):
        lim = 2 if ratio=='NC' else 3
        df_lim = df[(df['Manipulation']==lim)].dropna(subset=[ratio]).reset_index()
        df_lim[f'Mod_{ratio}'] = df_lim[f'Mod_{ModOpts[ratio]['models'][i]}_{ratio}']
        Plotting.Scatterplot_ObsVsModelStoich(axes[i],df_lim,ratio,'Phylum_Plotting',
                                              rmse[i],r2_loc=r2[i],
                                              axes_max=0.06 if (ratio=='PC' and i==1) else None)
        if i>0: axes[i].set(ylabel=None)

    # Create legend
    legend_items = ['Phylum:', *sorted(df_lim['Phylum_Plotting'].unique()), '',
                    'Healey, 1985 data', 'Linear regression', '1:1 line', 'Redfield star']
    Plotting.LegendContents(axes[1], legend_items, (1.4,1.1))
    if pdf==True: Plotting.saveFig(fig,f'FigS{'5' if ratio=='NC' else '9'}_{ratio}_ModelOpts')
    return None


def FigS6_NCvsGrowth_Plim(df,pdf):
    """
    Produce Figure S6: Growth rate vs N:C under P limitation, colored
    by irradiance
    """
    # Set up figure
    Plotting.figset()
    fig = plt.figure(layout='constrained', figsize=(3.2, 2.5))
    ax1=plt.gca()

    # Plot ratio under relevant limitation
    df_lim = df[(df['Manipulation']==3)].dropna(subset=['NC']).reset_index()
    Plotting.Scatterplot_GrowthVsStoich(ax1,df_lim,'NC','I_norm','Linear',cbar_min=9.9)
    
    # Create legend
    Plotting.Add_Colorbar(fig,ax1,df_lim,[10,100,600],cbar_min=9.9)
    Plotting.LegendContents(fig,['Healey, 1985 data','Linear regression',
                                 'Redfield line'],(0.97,0.33))
    if pdf==True: Plotting.saveFig(fig,'FigS6_NC_Plim_Irrad')
    return None


def FigS7_NC_Plim(df,pdf):
    """
    Produce Figure S7: Observed vs model N:C under P limitation, with annotation
    to show magnitude of potential N storage
    """
    # Set up figure
    Plotting.figset()
    fig = plt.figure(layout='constrained',figsize=(3.5,2.25))
    ax1=plt.gca()

    # Plot N:C under P limitation
    df_Plim = df[(df['Manipulation']==3)].dropna(subset=['NC']).reset_index()
    Plotting.Scatterplot_ObsVsModelStoich(ax1,df_Plim,'NC','Phylum_Plotting',(0.1,0.88))
    Plotting.Annotate_NSto(ax1, df_Plim) # Annotate potential N storage

    # Create legend
    legend_items = ['Phylum:', *sorted(df_Plim['Phylum_Plotting'].unique()), '',
                    'Healey, 1985 data', 'Linear regression', '1:1 line', 'Redfield star']
    Plotting.LegendContents(ax1, legend_items, (1.5,1.05))
    if pdf==True: Plotting.saveFig(fig,'FigS7_NC_Plim')
    return None


####################################
# Run Figures
####################################

def runFigures_MainText(df,pdf):
    """Produce and save figures in main text of paper"""
    Fig1_ObsStoich(df,pdf)
    Figs_2_5_ObsCNP_Irrad(df,'NC',pdf)
    Figs_3_6_RelevantLim(df,'NC',pdf)
    Fig4_NC_Replete(df,pdf)
    Figs_2_5_ObsCNP_Irrad(df,'PC',pdf)
    Figs_3_6_RelevantLim(df,'PC',pdf)
    return None


def runFigures_SI(df,pdf):
    """Produce and save figures in supporting information"""
    FigS1_Healey(df,pdf)
    FigS2_ObsStoich_Reverse(df,pdf)
    FigS3_S8_CultureType(df,'NC',pdf)
    FigS4_NC_Nlim_Phylum(df,pdf)
    Figs_S5_S9_ModelOpts(df,'NC',pdf)
    FigS6_NCvsGrowth_Plim(df,pdf)
    FigS7_NC_Plim(df,pdf)
    FigS3_S8_CultureType(df,'PC',pdf)
    Figs_S5_S9_ModelOpts(df,'PC',pdf)
    return None


def runFigures(main_text=True,SI=True,save=True):
    """Get model output and produce figures in paper"""
    # Get model output/analysis and prep for plotting
    df = ApplyModel.Model_df()
    df = Plotting.PrepData(df)

    # Produce and save figures in main text and SI
    if main_text==True: runFigures_MainText(df,save)
    if SI==True: runFigures_SI(df,save)
    return None