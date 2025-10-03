"""
Set of functions used for plotting figures
"""

# Import packages
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.cm import ScalarMappable
import seaborn as sns
from sklearn.metrics import root_mean_squared_error
from matplotlib.legend_handler import HandlerTuple
import ApplyModel
import CellModel
from pathlib import Path

# Global variables
Redfield = {'C':106, 'N':16, 'P':1}
FigFolder = 'Figures'


def figset():
    """Matplotlib figure settings"""
    plt.rcParams.update({
                 'font.family': 'Times New Roman',
                 'axes.linewidth' : 0.5,
                 'lines.linewidth' : 1,
                 'lines.markersize': 4,
                 'lines.markeredgewidth':0.5,
                 'font.size': 11,
                 'axes.titlesize': 11,
                 'axes.labelsize': 11,
                 'axes.labelpad': 3,
                 'axes.spines.top': False,
                 'axes.spines.right': False,
                 'figure.subplot.wspace': 0.25,
                 'figure.subplot.hspace': 0.35,
                 'figure.constrained_layout.w_pad': 0.05,
                 'figure.constrained_layout.h_pad': 0.05,
                 'xtick.labelsize': 10,
                 'ytick.labelsize': 10,
                 'legend.fontsize': 9,
                 'legend.title_fontsize': 9,
                 'legend.markerscale': 1,
                 'legend.borderpad': 0.3,
                 'legend.labelspacing': 0.3,
                 'legend.handlelength': 1.5,
                 'legend.handletextpad': 0.3,
                 'legend.columnspacing': 1.4,
                 'savefig.bbox': 'tight',
                 'savefig.pad_inches': 0.05,
                 'figure.dpi': 600,
                 'text.latex.preamble': r'\usepackage[greek,english]{babel}',
                 'text.latex.preamble': r'\usepackage{amsmath}',
                 "mathtext.fontset" : "stix"})
    return None


def saveFig(fig, fileName):
    """Save figure as png file"""
    # Check that output folder exists and build if not
    Path(CellModel.RootFolder/FigFolder).mkdir(parents=True,
                                               exist_ok=True)
    
    # Save figure and clear
    plt.savefig('%s/%s'%(CellModel.RootFolder/FigFolder,fileName))
    fig.clf()
    return None


def ColorPalette_Phylum():
    """Palette that assigns color to phylum,
    following Okabe and Ito, 2008 (https://jfly.uni-koeln.de/color/#pallet)"""
    OkabeIto_palette = {'Bacillariophyta': "#0072B2",
                        'Chlorophyta': "#009E73",
                        'Cryptophyta': "#F0E442",
                        'Cyanobacteria': "#D55E00",
                        'Dinoflagellata': "#CC79A7",
                        'Haptophyta': "#56B4E9",
                        'Ochrophyta': "#E69F00",
                        'Other': "#000000"}
    return OkabeIto_palette


def PrepData(df):
    """Add additional columns to dataframe used for plotting"""
    # Condense phyla for plotting
    df['Phylum_Plotting'] = df['Phylum']
    df.loc[df['Phylum'].isin(['Euglenozoa','Rhodophyta','Chlorarachniophyta',
                              'Glaucophyta']), 'Phylum_Plotting'] = 'Other'

    # Different marker for Healey data points
    df['Marker'] = 1
    df.loc[df['Ref.'] == '(Healey, 1985)', 'Marker'] = 2

    # Log of irradiance normalized to day length
    df['I_norm'] = df['Irrad.'] * df['Day.length']/24
    return df


def LabelSubplot(axes, xcoord=-0.1, ycoord=1.02, titles=None):
    """Add alphabetical labels to subplots"""
    axs = axes.flat
    for n, ax in enumerate(axs):
        if type(xcoord) == float:
            xcoord = [xcoord]*len(axs); ycoord = [ycoord]*len(axs)
        ax.text(xcoord[n], ycoord[n], '%s)'%(chr(ord('a')+n)) + '' if titles==None
                                    else '%s) '%(chr(ord('a')+n)) + titles[n],
                                    fontsize=11, transform=ax.transAxes)
    return axes


def LegendContents(axes, labels, loc, cols=1, new_l=[], new_h=[]):
    """Create legend from dictionary of label and handle options"""
    # Dictionary with options for legends and handles
    LegDict = {'Phylum:': plt.plot([],marker='', ls='')[0],
               '': plt.plot([],marker='', ls='')[0], # Empty line
               'Healey, 1985 data': plt.Line2D([],[],marker='D',linestyle='None',
                                               markeredgewidth=0.32,c='k'),
               'Regression': plt.Line2D([],[],linestyle='dashed',c='k'),
               'Linear regression': plt.Line2D([],[],linestyle='dashed',c='k'),
               'Quadratic regression': plt.Line2D([],[],linestyle='dashed',c='k'),
               '1:1 line': plt.Line2D([],[],c='k'),
               'Mean ratio': plt.Line2D([],[],linestyle='dashdot',c='k'),
               'Model $\\mu_{max}$': plt.Line2D([],[],linestyle='dashed',c='k'),
               'Model stoichiometry': plt.Line2D([],[],c='k'),
               r'Irradiance ($\mathrm{µE \; m^{-2} \; s^{-1}}$):': plt.plot([],marker='',
                                                                            ls='')[0],
               'Redfield line': plt.Line2D([],[],linestyle='dashed',c='r'),
               'Redfield star': (plt.Line2D([],[],marker='$\u2217$',ls='',c='r'),
                                 plt.Line2D([],[],linestyle='dashed',c='r'))}

    for phylum, hex in ColorPalette_Phylum().items():
        LegDict[phylum] = plt.Line2D([],[],marker='o',markeredgewidth=0.32,
                                     markeredgecolor='k',markerfacecolor=hex,
                                     linestyle='None')
    for lab, hand in zip(new_l,new_h):
        LegDict[lab] = hand
    
    # Build legend labels and handles from function input
    h = [LegDict[item] for item in labels]
    l = ['Redfield ratio' if 'Redfield' in item else item for item in labels]

    # Add legend to figure
    leg = axes.legend(h,l,loc='upper center',bbox_to_anchor=loc,frameon=False,
                      ncol=cols,handler_map={tuple: HandlerTuple(ndivide=None)})
    if 'Phylum:' in labels: leg.texts[0].set_position((-100, 0)) 
    if r'Irradiance ($\mathrm{µE \; m^{-2} \; s^{-1}}$):' in labels: leg.texts[0].set_position((-150, 0)) 
    return axes


def Add_Colorbar(fig, axes, df, labs, cbar_min=None):
    """Add colorbar for normalized irradiance"""
    cbar_range = LogNorm(df['I_norm'].min() if cbar_min==None else cbar_min,
                         df['I_norm'].max())
    sm =  ScalarMappable(norm=cbar_range,
                         cmap=sns.color_palette("magma", as_cmap=True))
    cax = axes.inset_axes([1.1, 0.22, 0.05, 0.68]) #([1.1, 0.2, 0.05, 0.7])
    cbar = fig.colorbar(sm, cax=cax)
    cbar.set_ticklabels([])
    [cbar.ax.text(2,l*0.97,l,ha='left',va='center',fontsize=10) for l in labs]
    cbar.set_label(r'$\mathrm{µE \; m^{-2} \; s^{-1}}$', size=10, labelpad=20)
    fig.text(0.8,0.95,'Normalized irradiance:',ha='left',va='center',fontsize=9)
    return axes


def CalcOutliers(dfs, ratio):
    """Calculate outliers 2.5x interquartile range beyond 25th/75th quartile"""
    for i in range(len(dfs)):
        for j in range(len(ratio)):
            q1,q3=dfs[i][ratio[j]].quantile([0.25,0.75])
            IQR=q3-q1
            dfs[i].loc[((dfs[i][ratio[j]]<(q1-2.5*IQR))|(dfs[i][ratio[j]]>(q3+2.5*IQR))),
                       ratio[j]] = np.nan
    return dfs


def Histogram_ObservedCNP(axes, df, ratio):
    """Plot histogram of observed stoichiometry for each culture condition"""
    # Figure variables
    pal = sns.cubehelix_palette(3, rot=-.25, light=.65)
    titles = ['N-limited','P-limited','Nutrient replete']

    # Split data by culture condition and calculate outliers
    dfs = [df[(df['Manipulation']==2)], df[(df['Manipulation']==3)],
           df[(df['Manipulation']==1)]]
    dfs_withOutliers = dfs.copy()
    dfs = CalcOutliers(dfs_withOutliers, ratio)

    # Loop through combinations of culture conditions and ratios
    for i in range(len(dfs)):
        for j in range(len(ratio)):
            sns.histplot(data=dfs[i].dropna(subset=[ratio[j]]), x=ratio[j],
                         stat='density', element="step", kde=True,
                         color=pal[i], ax=axes[i,j]) # Plot histogram
            # Add Redfield ratio line
            Rf = Redfield[ratio[j][0]]/Redfield[ratio[j][1]]
            if i==2: axes[i,j].axvline(Rf, ymin=0, ymax=3.56, color='r',
                                       linestyle='dashed', clip_on=False)
            # Add mean line (include outliers in calculating mean)
            kdeline = axes[i,j].lines[0]
            mean = dfs_withOutliers[i].dropna(subset=[ratio[j]])[ratio[j]].mean()
            xs = kdeline.get_xdata(); ys = kdeline.get_ydata()
            height = np.interp(mean, xs, ys)
            axes[i,j].vlines(mean, 0, height, color='k', linestyle='dashdot',
                             label='Mean ratio')
            # Add median line (include outliers in calculating mean)
            #median = dfs_withOutliers[i].dropna(subset=[ratio[j]])[ratio[j]].median()
            #height = np.interp(median, xs, ys)
            #axes[i,j].vlines(median, 0, height, color='k', linestyle='dotted',
            #                 label='Median ratio')

            # Format axes  
            axes[i,j].set(xlabel=r'$(\mathrm{mol \: %s}) \, (\mathrm{mol \: %s})^{-1}$'%(ratio[j][0],ratio[j][1]) if i==2 else None)       
            if j==0: axes[i,j].text(-0.5, 0.5, titles[i], color=pal[i],
                                    fontweight='bold', ha='center', va='center',
                                    rotation='vertical', transform=axes[i,j].transAxes)
            if j>0: axes[i,j].set(ylabel="")
            if i<2: axes[i,j].axes.get_xaxis().set_visible(False)
    return axes


def Histogram_ObservedCNP_V2(axes, df, ratio):
    """Plot histogram of observed stoichiometry for each culture condition"""
    # Figure variables
    pal = sns.cubehelix_palette(3, rot=-.25, light=.65)
    titles = ['Nutrient replete','N-limited','P-limited']

    # Split data by culture condition and calculate outliers
    dfs = [df[(df['Manipulation']==1)], df[(df['Manipulation']==2)],
           df[(df['Manipulation']==3)]]
    dfs_withOutliers = dfs.copy()
    dfs = CalcOutliers(dfs_withOutliers, ratio)

    # Loop through combinations of culture conditions and ratios
    for i in range(len(dfs)):
        for j in range(len(ratio)):
            sns.histplot(data=dfs[i].dropna(subset=[ratio[j]]), x=ratio[j],
                         stat='density', element="step", kde=True,
                         color=pal[i], ax=axes[i,j]) # Plot histogram
            # Add Redfield ratio line
            Rf = Redfield[ratio[j][0]]/Redfield[ratio[j][1]]
            if i==2: axes[i,j].axvline(Rf, ymin=0, ymax=3.56, color='r',
                                       linestyle='dashed', clip_on=False)
            # Add mean line (include outliers in calculating mean)
            kdeline = axes[i,j].lines[0]
            mean = dfs_withOutliers[i].dropna(subset=[ratio[j]])[ratio[j]].mean()
            xs = kdeline.get_xdata(); ys = kdeline.get_ydata()
            height = np.interp(mean, xs, ys)
            axes[i,j].vlines(mean, 0, height, color='k', linestyle='dashdot',
                             label='Mean ratio')

            # Format axes  
            axes[i,j].set(xlabel=r'$(\mathrm{mol \: %s}) \, (\mathrm{mol \: %s})^{-1}$'%(ratio[j][0],ratio[j][1]) if i==2 else None)       
            if j==0: axes[i,j].text(-0.5, 0.5, titles[i], color=pal[i],
                                    fontweight='bold', ha='center', va='center',
                                    rotation='vertical', transform=axes[i,j].transAxes)
            if j>0: axes[i,j].set(ylabel="")
            if i<2: axes[i,j].axes.get_xaxis().set_visible(False)
    return axes


def PlotRegression(axes, type, df, varx, vary, xlim, r_loc):
    """Plot regression line (and statistic) for given variables"""
    # Range to calculate regression fit over
    xseq = np.linspace(xlim[0], xlim[1], num=100)

    # Get linear or quadratic regression
    if type == 'Linear':
        regress, r = ApplyModel.LinReg_RMA(df, varx, vary, xseq)
    elif type == 'Quadratic':
        regress = ApplyModel.QuadReg(df, varx, vary, xseq)
    
    # Plot output
    axes.plot(xseq, regress, color='k', alpha=0.8, linestyle='dashed',
              label='%s regression'%(type))
    
    # Annotate with relevant statistic
    r_stat = 'r' if varx == 'Specific.growth' else '$R^2$'
    n = 1 if varx == 'Specific.growth' else 2
    if r_loc != None:
        axes.text(r_loc[0], r_loc[1], f"{r_stat}: {r**n:.2f}", fontsize=10,
                  color='k', alpha=0.8, ha='left', va='bottom',
                  transform=axes.transAxes)
    return axes


def Scatterplot_GrowthVsStoich(axes, df, ratio, hue_col, regression=None,
                               Rf_ratio=True, cbar_min=None):
    """Scatterplot of growth rate vs given stoichiometric ratio"""
    # Palette for plotting (either phylum or irradiance)
    if hue_col == 'Phylum_Plotting':
        pal = ColorPalette_Phylum()
    elif hue_col == 'Irrad.':
        pal = sns.color_palette("magma", len(df['Irrad.'].unique()))
    elif hue_col == 'I_norm':
        pal = sns.color_palette("magma", as_cmap=True)
        hnorm = LogNorm(df['I_norm'].min() if cbar_min==None else cbar_min,
                        df['I_norm'].max())

    # Plot growth rate vs ratio
    sns.scatterplot(data=df, x='Specific.growth', y=ratio,hue=hue_col,
                    palette=pal, hue_norm=hnorm if hue_col=='I_norm' else None,
                    edgecolor='k', style='Marker', markers={1:'o',2:'D'},
                    legend=False, ax=axes)
    
    # Axes limits
    xlim = [df['Specific.growth'].min()-0.06,
            round(1.05 * df['Specific.growth'].max())] #ADJUST AS NEEDED FOR DIFF FIGURES
    axes.set_xlim(left=xlim[0], right=xlim[1]) #ADJUST AS NEEDED FOR DIFF FIGURES
    axes.set_ylim(bottom=0, top=df[ratio].max()*1.05) #ADJUST AS NEEDED FOR DIFF FIGURES

    # Plot regression as specified by inputs
    r_loc = (0.7, 0.2) if regression=='Linear' else None #ADJUST X,Y AS NEEDED FOR DIFF FIGURES
    if regression != None: PlotRegression(axes, regression, df, 'Specific.growth',
                                           ratio, xlim, r_loc)

    # Add marker for Redfield Ratio
    Rf = Redfield[ratio[0]]/Redfield[ratio[1]]
    if Rf_ratio==True: axes.axline((0, Rf), color='r', linestyle='dashed',
                                   slope=0, label='Redfield ratio')

    # Axes labels
    xlab = r'Growth Rate $\left( \mathrm{d^{-1}} \right)$'
    ylab = r'Observed $(\mathrm{mol \; %s}) \: (\mathrm{mol \; %s})^{-1}$'%(ratio[0],ratio[1])
    axes.set(xlabel=xlab, ylabel=ylab)
    return axes


def Scatterplot_ObsVsModelStoich(axes, df, ratio, rmse_loc, r2_loc=(0.7, 0.3), axes_max=None):
    """Scatterplot of observed vs model values for given stoichiometric ratio"""
    # Plot stoichiometry
    sns.scatterplot(data=df, x=ratio, y='Mod_%s'%(ratio), hue='Phylum_Plotting',
                    palette=ColorPalette_Phylum(), edgecolor='k', style='Marker',
                    markers={1:'o',2:'D'}, legend=False, ax=axes)
    
    # Plot 1:1 line between observed and model stoichiometry
    axes.axline((0, 0), color='k', slope=1)

    # Calculate and print RMSE value
    rmse = root_mean_squared_error(df[ratio], df['Mod_%s'%ratio])
    axes.text(rmse_loc[0], rmse_loc[1], f"$RMSE$: {rmse:.3f}", fontsize=10,
                 ha="left", va="bottom", transform=axes.transAxes)
    
    # Axes limits
    max = np.max([df[ratio].max(),df['Mod_%s'%ratio].max()]) if axes_max==None else axes_max
    axes.set_xlim(left=0, right=max*1.05)
    axes.set_ylim(bottom=0, top=max*1.05)

    # Plot linear regression and R2 value
    PlotRegression(axes, 'Linear', df, ratio, 'Mod_%s'%ratio, [0, max*1.05],
                   r2_loc) #ADJUST X,Y AS NEEDED FOR DIFF FIGURES)
    
    # Add markers for Redfield Ratio
    Rf = Redfield[ratio[0]]/Redfield[ratio[1]]
    axes.text(0, Rf, '$\u2217$-', color='r', ha='right', va='center')
    axes.text(Rf, 0, '-$\u2217$', color='r', rotation=270, ha='center', va='top')

    # Axes labels
    lab = r'$(\mathrm{mol \; %s}) \: (\mathrm{mol \; %s})^{-1}$'%(ratio[0],ratio[1])
    axes.set(xlabel='Observed %s'%(lab), ylabel='Model %s'%(lab))
    return axes


def Annotate_NSto(axes, df):
    """Add annotation showing maximum potential N storage"""
    # Plot sample points at determined coordinates
    pts_x = [0.205,0.205]
    linreg = ApplyModel.LinReg_RMA(df, 'NC', 'Mod_NC', pts_x[0])[0]
    pts_y = [linreg, linreg + 0.035]
    sns.scatterplot(x=pts_x,y=pts_y,color='#9da09e',edgecolor='k',linewidth=0.6,s=24,zorder=3,legend=False,ax=axes)

    # Add arrow and text to denote magnitude of potential N storage
    axes.annotate('', xy=(pts_x[1],pts_y[1]), xytext=(pts_x[0],pts_y[0]),
                  arrowprops = dict(color='k', width=0.4, headwidth=2.5,
                                    headlength=2.5, shrink=0.215))
    axes.annotate(r'$+$ max N', xy=(pts_x[1],pts_y[1]), xytext=(pts_x[0]+.01,pts_y[0]+0.03),
                  verticalalignment='center', horizontalalignment='left', fontsize=10,
                  bbox=dict(boxstyle='square,pad=0', fc='white', alpha=1, edgecolor='none'))
    axes.annotate('storage', xy=(pts_x[1],pts_y[1]), xytext=(pts_x[0]+.02,pts_y[0]+0.014),
                  verticalalignment='center', horizontalalignment='left', fontsize=10,
                  bbox=dict(boxstyle='square,pad=0.05', fc='white', alpha=1, edgecolor='none'))
    return axes


def Annotate_NC_Replete(axes):
    """Add annotations to plot of growth vs N:C in nutrient replete cultures"""
    # Add arrows
    arrow = dict(color='k', width=1.5, headwidth=6, headlength=5, shrink=0.03)
    axes.annotate('', xy=(2.3,0.14), xytext=(0,0.02), arrowprops=arrow)
    axes.annotate('', xy=(2.3,0.16), xytext=(0,0.31), arrowprops=arrow) #xytext=(0,0.305)

    #Add Labels for limitation
    axes.text(1, 0.27, 'Light limited growth', verticalalignment='center',
              horizontalalignment='left', fontsize=10) 
    axes.text(1, 0.04, 'Non-light limited growth', verticalalignment='center',
              horizontalalignment='left', fontsize=10) 
    return axes


def Allocation_Schematic(axes, xlab, ylab_loc, y_NC, mu, photo, bio):
    """Schematic of growth rate vs cell N quota in response to limitation"""
    # Axes Settings
    axes.tick_params(left=False, labelleft=False, bottom=False, labelbottom=False) 
    #axes.set(xlabel=xlab, ylabel='N Quota')
    axes.set(xlabel=xlab, ylabel='N allocation'+' '*ylab_loc) #'N investment','N usage'
    axes.set_xlim(left=0,right=2)

    # Plot N quota
    pal = ['#b8b8b8','#a00000','#298c8c']  # color palette
    axes.stackplot([0,2],y_NC,ec='w',lw=0.5,colors=(pal[0],pal[1],pal[2]))

    # Annotate photosynthetic and biosynthetic protein pools
    axes.text(photo[0],photo[1],'Photo',verticalalignment='center',
              horizontalalignment='left',fontsize=10) 
    axes.text(bio[0],bio[1],'Bio',verticalalignment='center',
              horizontalalignment='left',fontsize=10) 

    # Add arrow for growth rate
    axes.plot([0,1.95],[mu[1],(1.95*mu[0]+mu[1])],color='k',lw=2.3)
    axes.annotate('', xy=(2,(2*mu[0]+mu[1])), xytext=(1.8,(1.8*mu[0]+mu[1])),
                  arrowprops=dict(color='k',width=1,headwidth=6,headlength=5))

    # Add broken axes in 'other' allocation
    axes.vlines(-0.004,0.0115,0.016,color='w',lw=0.6,zorder=3,clip_on=False)
    axes.plot([0,0],[0.0115,0.016],marker=[(-1, -0.5),(1, 0.5)],markersize=6,
              linestyle="none",color='k',mew=0.7,mec='k',zorder=4,
              clip_on=False)  # remove white space between axis break marks
    return axes


def Extrapolate_Healey(df):
    """Extrapolate model output values for Healey, 1985"""
    # Model parameters for calculations
    param = ApplyModel.readParam('Healey')

    # Calculate model N:C and P:C at each irradiance from growth of 0 to µ_max
    df['Mod_umax'] = CellModel.mu_max(df,param[0:8])
    mu = [np.linspace(0,x,round(x/0.05)) for x in df['Mod_umax'].unique()]
    Irr = np.repeat(df['Irrad.'].unique(),[len(x) for x in mu])
    df_Mod = pd.DataFrame({'Specific.growth': np.concatenate(mu),
                           'Irrad.': Irr,
                           'Temp.': [df['Temp.'].unique()[0]] * len(Irr),
                           'V': [df['V'].unique()[0]] * len(Irr),
                           'Day.length': [df['Day.length'].unique()[0]] * len(Irr)})
    df_Mod['Mod_NC'] = CellModel.NC_eval(df_Mod,[*param[0:5],param[7]])
    df_Mod['Mod_PC'] = CellModel.PC_eval(df_Mod,[*param[0:6],*param[7:9]])

    # Calculate model N:C and P:C at each µ_max for irradiances from 2-180
    I = np.arange(2,180,2)
    df_umax = pd.DataFrame({'Irrad.': I,
                            'Temp.': [df['Temp.'].unique()[0]] * len(I),
                            'V': [df['V'].unique()[0]] * len(I),
                            'Day.length': [df['Day.length'].unique()[0]] * len(I)})
    df_umax['Specific.growth'] = CellModel.mu_max(df_umax,param[0:8])
    df_umax['Mod_NC'] = CellModel.NC_eval(df_umax,[*param[0:5],param[7]])
    df_umax['Mod_PC'] = CellModel.PC_eval(df_umax,[*param[0:6],*param[7:9]])
    return df_Mod, df_umax


def Healey_ModelFit(axes, df, ratio):
    """Add model output for Healey, 1985 data"""
    # Get model predictions
    df_Mod, df_umax = Extrapolate_Healey(df)

    # Plot model prediction for N:C or P:C
    pal = sns.color_palette("magma", len(df['Irrad.'].unique()))
    sns.lineplot(data=df_Mod, x='Specific.growth', y=f'Mod_{ratio}',
                 hue='Irrad.', palette=pal, ax=axes)

    # Plot model maximum growth rate
    sns.lineplot(data=df_umax,x='Specific.growth', y=f'Mod_{ratio}', color='k',
                 linestyle='--', ax=axes)
    
    # Add marker for Redfield Ratio
    Rf = Redfield[ratio[0]]/Redfield[ratio[1]]
    axes.text(0, Rf, '$\u2217$-', color='r', ha='right', va='center')
    
    # Axes limits
    axes.autoscale()
    axes.set_xlim(xmin=0); axes.set_ylim(ymin=0)
    return axes