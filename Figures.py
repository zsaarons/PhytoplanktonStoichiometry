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

def Fig1_ObsStoich(df,png):
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
    if png==True: Plotting.saveFig(fig,'Fig1_HistCNP')
    return None


def Figs_2_7_ObsCNP_Irrad(df,ratio,png):
    """
    Produce Figures 2 and 7: Growth rate vs N:C and P:C under N and P
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
    if png==True: Plotting.saveFig(fig,f'Fig{'2' if ratio=='NC' else '7'}_{ratio}_Irrad')
    return None


def Figs_3_8_RelevantLim(df,ratio,png):
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
    Plotting.Scatterplot_ObsVsModelStoich(ax1,df_lim,ratio,
                                          (0.1,0.88 if ratio=='NC' else 0.68))

    # Create legend
    legend_items = ['Phylum:', *sorted(df_lim['Phylum_Plotting'].unique()), '',
                    'Healey, 1985 data', 'Linear regression', '1:1 line',
                    'Redfield star']
    Plotting.LegendContents(ax1, legend_items, (1.5,1.07))
    if png==True: Plotting.saveFig(fig,'Fig3_NC_Nlim' if ratio=='NC' else 'Fig8_PC_Plim')
    return None


def Fig4_NC_Nlim_Phylum(df,png):
    """
    Produce Figure 4: Observed vs model N:C under N limitation, seperated by
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
        if i<7: Plotting.Scatterplot_ObsVsModelStoich(ax,df_Nlim[(df_Nlim['Phylum_Plotting']==phyla[i])],'NC',(0.1,0.85),r2_loc=(0.65,0.25),axes_max=maxAx)
        if i%4>0: ax.set(ylabel=None)
        ax.set(xlabel=None if i<3 else 'Obs. '+r'$(\mathrm{mol \; N}) \: (\mathrm{mol \; C})^{-1}$')

    legend_items = ['Healey, 1985 data', 'Linear regression', '1:1 line', 'Redfield star']
    Plotting.LegendContents(fig, legend_items, (0.9,0.4))
    if png==True: Plotting.saveFig(fig,'Fig4_NC_Nlim_Phylum')
    return None


def Fig5_NC_Plim(df,png):
    """
    Produce Figure 5: Observed vs model N:C under P limitation, with annotation
    to show magnitude of potential N storage
    """
    # Set up figure
    Plotting.figset()
    fig = plt.figure(layout='constrained',figsize=(3.5,2.25))
    ax1=plt.gca()

    # Plot N:C under P limitation
    df_Plim = df[(df['Manipulation']==3)].dropna(subset=['NC']).reset_index()
    Plotting.Scatterplot_ObsVsModelStoich(ax1,df_Plim,'NC',(0.1,0.88))
    Plotting.Annotate_NSto(ax1, df_Plim) # Annotate potential N storage

    # Create legend
    legend_items = ['Phylum:', *sorted(df_Plim['Phylum_Plotting'].unique()), '',
                    'Healey, 1985 data', 'Linear regression', '1:1 line', 'Redfield star']
    Plotting.LegendContents(ax1, legend_items, (1.5,1.05))
    if png==True: Plotting.saveFig(fig,'Fig5_NC_Plim')
    return None


def Fig6_NC_Replete(df,png):
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
    df_Rep = df[(df['Manipulation']==1)].drop(df.loc[df['Ref.'].isin(['(Fernandez-Reiriz et al, 1989)'])].index).dropna(subset=['NC']).reset_index()

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
    if png==True: Plotting.saveFig(fig,'Fig6_NC_Replete')
    return None


####################################
# Supporting Information
####################################

def FigS1_Healey(df,png):
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
    if png==True: Plotting.saveFig(fig,'FigS1_HealeyCNP')
    return None


def FigS2_ObsStoich_Reverse(df,png):
    """
    Produce Figure S2: Histograms of observed stoichiometry (C:N, C:P)
    across nutrient manipulations (N limited, P limited, nutrient replete)
    """
    # Set up figure
    Plotting.figset()
    fig, axes = plt.subplots(3,2,layout='constrained',figsize=(4,4.5))
    Plotting.LabelSubplot(axes, xcoord=-.13) # label axes

    # Plot histograms
    Plotting.Histogram_ObservedCNP(axes, df, ['CN','CP'])
    for i in range(len(axes)): # Set y-axes to be uniform within each row
        for j in range(len(axes[0])):
            xmax = j; ymax = 2
            axes[i,j].set_ylim(bottom=0,top=axes[ymax,j].get_ylim()[1])
            axes[i,j].set_xlim(left=0,right=axes[xmax,j].get_xlim()[1])

    # Create legend
    Plotting.LegendContents(fig, ['Redfield line', 'Mean ratio'], (0.9,0.98))
    if png==True: Plotting.saveFig(fig,'FigS2_HistCNP_Reverse')
    return None


def Figs_S3_S5_ModelOpts(df,ratio,png):
    """
    Produce Figures S3 and S5: Observed vs model N:C and P:C under N and P
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
        Plotting.Scatterplot_ObsVsModelStoich(axes[i],df_lim,ratio,rmse[i],r2_loc=r2[i])
        if i>0: axes[i].set(ylabel=None)

    # Create legend
    legend_items = ['Phylum:', *sorted(df_lim['Phylum_Plotting'].unique()), '',
                    'Healey, 1985 data', 'Linear regression', '1:1 line', 'Redfield star']
    Plotting.LegendContents(axes[1], legend_items, (1.4,1.1))
    if png==True: Plotting.saveFig(fig,f'FigS{'3' if ratio=='NC' else '5'}_{ratio}_ModelOpts')
    return None


def FigS4_NCvsGrowth_Plim(df,png):
    """
    Produce Figure S4: Growth rate vs N:C under P limitation, colored
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
    if png==True: Plotting.saveFig(fig,'FigS4_NC_Plim')
    return None


def FigS6_PC_NlimReplete(df,png):
    """
    Produce Figure S6: Growth rate vs P:C under N limitation (panel a) and
    nutrient replete conditions (panel b)
    """
    # Set up figure
    Plotting.figset()
    fig, axes = plt.subplots(1,2,layout='constrained',figsize=(6,2.5))
    Plotting.LabelSubplot(axes,titles=['N-limited','Nutrient replete'], ycoord=1.07)

    # Plot P:C under N limitation and replete conditions
    df_Nlim = df[(df['Manipulation']==2)].dropna(subset=['PC']).reset_index()
    Plotting.Scatterplot_GrowthVsStoich(axes[0],df_Nlim,'PC','Phylum_Plotting')
    # Exclude outliers (P:C>0.07 in nutrient replete; 2 data pts from 
    # 'Finkel et al, 2006' and 'McKew et al, 2015'
    df_Rep = df[(df['Manipulation']==1)&(df['PC']<0.07)].dropna(subset=['PC']).reset_index()
    Plotting.Scatterplot_GrowthVsStoich(axes[1],df_Rep,'PC','Phylum_Plotting')
    axes[1].set(ylabel=None)

    # Create legend
    legend_items = ['Phylum:', *sorted(df_Rep['Phylum_Plotting'].unique()), '',
                    'Healey, 1985 data', 'Redfield line']
    Plotting.LegendContents(axes[1], legend_items, (1.4,1.05))
    if png==True: Plotting.saveFig(fig,'FigS6_PC_NlimReplete')
    return None


####################################
# Run Figures
####################################

def runFigures_MainText(df,png):
    """Produce and save figures in main text of paper"""
    Fig1_ObsStoich(df,png)
    Figs_2_7_ObsCNP_Irrad(df,'NC',png)
    Figs_3_8_RelevantLim(df,'NC',png)
    Fig4_NC_Nlim_Phylum(df,png)
    Fig5_NC_Plim(df,png)
    Fig6_NC_Replete(df,png)
    Figs_2_7_ObsCNP_Irrad(df,'PC',png)
    Figs_3_8_RelevantLim(df,'PC',png)
    return None


def runFigures_SI(df,png):
    """Produce and save figures in supporting information"""
    FigS1_Healey(df,png)
    FigS2_ObsStoich_Reverse(df,png)
    Figs_S3_S5_ModelOpts(df,'NC',png)
    FigS4_NCvsGrowth_Plim(df,png)
    Figs_S3_S5_ModelOpts(df,'PC',png)
    FigS6_PC_NlimReplete(df,png)
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