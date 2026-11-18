"""
Produce tables used in main text and supporting information as .tex files
"""

# Import packages
import numpy as np
import pandas as pd
import re
from IPython.display import display
import ApplyModel
from pathlib import Path
import CellModel

TableFolder = 'Tables'


def saveTable(table, fileName):
    """Save table as latex file"""
    # Check that output folder exists and build if not
    Path(CellModel.RootFolder/TableFolder).mkdir(parents=True,
                                                 exist_ok=True)
    # Save table
    with open('%s/%s.tex'%(CellModel.RootFolder/TableFolder,fileName),
              'w') as table_file:
        table_file.write(table)
    return None


####################################
#  Main Text
####################################

def Table1_StoichDatabase(df,save):
    """
    Produce Table 1: Summary of database contents
    """
    # Remove rows that do not contain any of N:C, P:C, or N:P
    df = df.copy().dropna(subset=['NC','PC','NP'],how='all')
    data_filter = ['Phylum','NC','PC','NP']
    col = ['Nutrient Manipulation','All Data','Contains N:C','Contains P:C','Contains N:P']
    manip = ['N-Limited','P-Limited','Nutrient Replete','Total']; manip_num = [2,3,1]

    # Extract data for table
    rows = []
    for i in range(len(manip)): # For each nutrient manipulation
        df_manip = df[(df['Manipulation'] == manip_num[i])] if i<3 else df
        data_pts = ['{:d} ({:d})'.format(len(df_manip.dropna(subset=[data_filter[j]])),len(df_manip.dropna(subset=[data_filter[j]])['Ref.'].unique())) for j in range(4)]
        rows.append([manip[i],*data_pts])

    # Compile data
    table = pd.DataFrame(rows,columns=col).style.hide(level=0, axis=0)
    tex = table.to_latex(caption='Summary of data points (studies) in stoichiometric data compilation.',
                         column_format='lcccc',
                         position='H',
                         hrules=True)
    tex = tex.replace("\\begin{tabular}", "\\small \\begin{tabular}")
    if save==True: saveTable(tex,'Table1_StoichiometryDatabase')
    #display(table)
    return None


def Table2_MeanStoich(df,save):
    """
    Produce Table 2: Mean values (and interquartile ranges) of C:N:P in data
    compilations for all nutrient regimes as compared to Redfield ratios
    """
    ratio = ['NC','PC','NP'] # Stoichiometries to plot in order
    Rf = {'C':106, 'N':16, 'P':1} # To calculate Redfield ratio
    nut_manipulation = ['N-Limited','P-Limited','Nutrient Replete']
    manip = [2,3,1]
    col = ['Nutrient Manipulation',
           r'$(\mathrm{mol \: N}) (\mathrm{mol \: C})^{-1}$', #N:C
           r'$(\mathrm{mol \: P}) (\mathrm{mol \: C})^{-1}$', #P:C
           r'$(\mathrm{mol \: N}) (\mathrm{mol \: P})^{-1}$'] #N:P

    # For each nutrient manipulation
    rows = []
    num_length = {'NC':3,'PC':4,'NP':1}    
    for i in range(len(nut_manipulation)):
        df_manip = df[(df['Manipulation'] == manip[i])]
        stoich = ['{:.{prec}f} ({:.{prec}f})'.format(np.mean(df_manip[j]),df_manip[j].quantile(0.75)-df_manip[j].quantile(0.25),prec=num_length[j]) for j in ratio]
        rows.append([nut_manipulation[i],*stoich])
            
    # Add row for Redfield Ratio
    rows.append(['Redfield Ratio','{:.3f}'.format(Rf['N']/Rf['C']),'{:.4f}'.format(Rf['P']/Rf['C']),'{:.1f}'.format(Rf['N']/Rf['P'])])

    # Compile data
    table = pd.DataFrame(rows,columns=col).style.hide(level=0, axis=0)
    tex = table.to_latex(caption='Mean cellular ratios (interquartile range) in literature compilation as'\
                         ' compared to Redfield ratio.',
                         column_format='lccc',
                         position='H',
                         hrules=True)
    tex = tex.replace("\\begin{tabular}", "\\small \\begin{tabular}")
    if save==True: saveTable(tex,'Table2_MeanObservedStoich')
    #display(table)
    return None


####################################
#  Supporting Information
####################################

def TableS1_StoichDatabasePhylum(df,save):
    """
    Produce Table S1: Summary of database contents by phylum
    """
    # Remove rows that do not contain any of N:C, P:C, or N:P
    df = df.copy().dropna(subset=['NC','PC','NP'],how='all')

    df['Phylum_Condensed'] = df['Phylum']
    df.loc[df['Phylum'].isin(['Cryptophyta','Euglenozoa','Ochrophyta',
                              'Rhodophyta','Chlorarachniophyta','Glaucophyta']),
                              'Phylum_Condensed'] = 'Other'
    PFT = np.sort(df['Phylum_Condensed'].unique())

    data_filter = ['','NC','PC','NP']
    data_subset = ['All Data','Contains N:C','Contains P:C','Contains N:P']
    manip = ['Nutrient Replete','N-Limited','P-Limited']

    col = pd.MultiIndex.from_arrays([[ *['\\it{Phylum}']*len(PFT)],[*PFT]])
    idx = ['\\rule{0pt}{2ex} \\it{Nutrient Replete}',*['\\quad %s'%j for j in data_subset],
           '\\noalign{\\vspace{0.25ex}} \\cline{1-7} \\rule{0pt}{2.6ex} \\it{N-Limited}',*['\\quad %s'%j for j in data_subset],
           '\\noalign{\\vspace{0.25ex}} \\cline{1-7} \\rule{0pt}{2.6ex} \\it{P-Limited}',*['\\quad %s'%j for j in data_subset],
           '\\noalign{\\vspace{0.25ex}} \\cline{1-7} \\rule{0pt}{2.6ex} \\it{Total}']
               
    # Extract data for table
    rows = []
    for i in range(len(manip)): # For each nutrient manipulation
        rows.append([*['']*(len(PFT))])
        df_manip = df[(df['Manipulation'] == i+1)]
        for j in range(len(data_filter)): # Containing given stoichiometric ratio
            df_filter = df_manip.dropna(subset=[data_filter[j]]) if j > 0 else df_manip
            taxa = ['{:d} ({:d})'.format(len(df_filter[(df_filter['Phylum_Condensed'] == k)]),len(df_filter[(df_filter['Phylum_Condensed'] == k)]['Ref.'].unique())) for k in PFT]
            rows.append(taxa)
    #rows.append([len(df[(df['Phylum_Condensed'] == k)]) for k in PFT])
    rows.append(['{:d} ({:d})'.format(len(df[(df['Phylum_Condensed'] == k)]),len(df[(df['Phylum_Condensed'] == k)]['Ref.'].unique())) for k in PFT])
    
    # Compile data
    table = pd.DataFrame(rows,columns=col,index=idx).style
    tex = table.to_latex(caption='Number of data points (studies) by phylum in stoichiometric data compilation.',
                          column_format='@{}lcccccc',
                          multicol_align="c",
                          position='H',
                          hrules=True).replace('& Studies','\\cmidrule{2-7} & Studies')
    tex = tex.replace("\\begin{tabular}", "\\small \\begin{tabular}").replace("\\begin{table}[H]", "\\begin{sidewaystable}").replace("\\end{table}", "\\end{sidewaystable}")
    if save==True: saveTable(tex,'TableS1_StoichiometryDatabase_PhylumSummary')
    #display(table)
    return None


def TableS2_ModelConstants(save):
    """
    Produce Table S1: Values for fixed model parameters
    """
    # Table info
    param = {r'$E$': [7.74E-1,'--'],
              r'$A^P_{RNA}$': [7.15E-3,r'mol P mol C\textsuperscript{-1} d'],
              r'$Q_{P,min}^{RNA}$': [2.23E-4,r'mol P mol C\textsuperscript{-1}'],
              r'$Q_C^{DNA}$': [9.41E-4,r'mol C mol C\textsuperscript{-1}'],
              r'$T_{ref}$': [25,r'{\degree}C'],
              r'$Q_{10(m)}$': [2.5,'--'],
              r'$Q_{10(\sigma_a))}$': [2.0,'--'],
              r'$Q_{10(Bio)}$': [2.52,'--'],
              r'$V_{ref}$': [1,r'{\unit{\um}}\textsuperscript{3}'],
              r'$b_m$': [-0.189,'--'],
              r'$b_{DNA}$': [-0.189,'--'],
              r'$b_{Pro-Other}$': [-0.717,'--'],
              r'$b_{OtherP}$': [-0.123,'--'],
              r'$b_{OtherC}$': [-0.809,'--']}

    #make table
    params = pd.DataFrame({'Parameter': param.keys(),
                           r'\; Value': [value[0] for value in param.values()],
                           'Units': [value[1] for value in param.values()]})

    # Style and convert to latex
    table = params.style.format({r'\; Value':'{:.2E}'}).hide(axis=0)
    tex = table.to_latex(caption=r'Values used for fixed parameters in CFM-Phyto.',
                         column_format='llc',
                         position='H',
                         hrules=True)
    tex = re.sub(r'(\d.\d{2})E(-?\+?\d{2})', r'\\num{\1d\2}', tex)
    tex = re.sub(r'(& \\)', r'& \\;\\', tex)
    if save==True: saveTable(tex,'TableS2_ModelParamTable_Constants')
    #display(table)
    return None


def TableS3_ModelFitOptions(save):
    """
    Produce Table S2: Values for model parameters fit to 3 datasets (Healey
    1985 study; all nutrient replete data in compilation; all nutrient replete
    data excluding cultures of Cyanobacteria)
    """
    # Get data and compile table
    param_unit = {r'$A_I$': r'$\frac{\text{mol C (mol C in Chl)\textsuperscript{-1} d\textsuperscript{-1}}}{\text{{\unit{\umol}} m\textsuperscript{-2} s\textsuperscript{-1}}}$',
                  r'$P_{max}$': r'mol C (mol C in Chl)\textsuperscript{-1} d\textsuperscript{-1}',
                  r'$m$': r'd\textsuperscript{-1}',
                  r'$A_{Bio}$': r'mol C mol C\textsuperscript{-1}',
                  r'$A_{Pho}$': r'mol C (mol C in Chl)\textsuperscript{-1}',
                  r'$A_{Pho}^{P:Chl}$': r'mol P (mol C in Chl)\textsuperscript{-1}',
                  r'$Q_{C}^{Other0}$': r'mol C mol C\textsuperscript{-1}',
                  r'$Q_{C}^{Pro-Other}$': r'mol C mol C\textsuperscript{-1}',
                  r'$Q_{P}^{Other0}$': r'mol P mol C\textsuperscript{-1}'}
    
    params = pd.DataFrame({'Parameter': param_unit.keys(),
                           'Healey, 1985': ApplyModel.readParam('Healey')[0:9],
                           'Nutrient-limited': ApplyModel.readParam('NutLim')[0:9],
                           'No Cyano': ApplyModel.readParam('NutLim_NoCyano')[0:9],
                           'Unit': param_unit.values()})
    
    # Remove because cant get actual value (derived from calculations of in u_max,
    # which requires nutrient replete data)
    params.loc[6,'Nutrient-limited'] = float(-1)
    params.loc[6,'No Cyano'] = float(-1)
    # Swap order of rows
    params = params.reindex([1,0,2,4,5,3,7,8,6])

    # Style and convert to latex
    table = params.style.format({'Healey, 1985':'{:.2E}',
                                 'Nutrient-limited':'{:.2E}',
                                 'No Cyano':'{:.2E}'}).hide(axis=0)
    tex = table.to_latex(caption=r'Values of CFM-Phyto parameters as fit to several datasets.',
                         column_format='llllc',
                         position='H',
                         hrules=True).replace('-1.00E+00','--')
    tex = re.sub(r'(\d.\d{2})E(-?\+?\d{2})', r'\\num{\1d\2}', tex)
    tex = re.sub(r'(\n\$.*?&\s)(.*?&\s)(.*?&\s)', r'\1\, \2\\quad\:\3\;', tex)
    tex = re.sub(r'(&\s)([^\d]*?)(--.*?&\s)([^\d]*?)(--)',
                 r'\1\\qquad\\quad \3\\qquad \5', tex)
    if save==True: saveTable(tex,'TableS3_ModelParamTable_FitOpt')
    #display(table)
    return None


def TableS4_OptimalGrowthCN(df,save):
    """
    Produce Table S3: Model prediction of maximum growth rate and C:N ratio in
    nutrient replete cultures with given temperature, day length, and cell size
    """
    # Rename columns
    table = df.rename(columns={'V': r'Cell Volume ($\mathrm{\mu m}^{3}$)',
                               'Temp.': r'Temperature (\degree C)',
                               'Day.length': r'Day length (hrs)',
                               'mu_max': r'$\mu_{max}$',
                               'CN': r'$(\mathrm{mol \: C}) (\mathrm{mol \: N})^{-1}$'}).style.hide(level=0, axis=0)
    table = table.format({r'$\mu_{max}$':'{:.2f}',
                          r'$(\mathrm{mol \: C}) (\mathrm{mol \: N})^{-1}$':'{:.2f}'}).hide(axis=0)
               
    # Write latex table
    table_caption = 'Model output for phytoplankton maximum growth rate and associated C:N \\\\ratio at given cell volume and environmental conditions (temperature and day length).'
    tex = table.to_latex(caption=table_caption,
                         column_format='cccccc',
                         position='H',
                         hrules=True)
    tex = tex.replace("\\begin{tabular}", "\\small \\begin{tabular}")
    if save==True: saveTable(tex,'TableS4_OptimalGrowthCN')
    #display(table)
    return None


def TableS5_StoichDatabase_CultureType(df,save):
    """
    Produce Table: Summary of culture types in database
    """
    # Remove rows that do not contain any of N:C, P:C, or N:P
    df = df.copy().dropna(subset=['NC','PC','NP'],how='all')

    data_filter = ['','NC','PC','NP']
    data_subset = ['All Data','Contains N:C','Contains P:C','Contains N:P']
    manip = ['Nutrient Replete','N-Limited','P-Limited']

    col = pd.MultiIndex.from_arrays([[ *['\\it{Culture System}']*4],
                                     [r'\;Batch\;',r'\;Semi-Continuous\;',r'\;Chemostat\;',r'\;Turbidostat\;']])
    idx = ['\\rule{0pt}{2ex} \\it{Nutrient Replete \\quad\\;}',*['\\quad %s'%j for j in data_subset],
           '\\noalign{\\vspace{0.25ex}} \\cline{1-5} \\rule{0pt}{2.6ex} \\it{N-Limited}',*['\\quad %s'%j for j in data_subset],
           '\\noalign{\\vspace{0.25ex}} \\cline{1-5} \\rule{0pt}{2.6ex} \\it{P-Limited}',*['\\quad %s'%j for j in data_subset],
           '\\noalign{\\vspace{0.25ex}} \\cline{1-5} \\rule{0pt}{2.6ex} \\it{Total}']
               
    # Extract data for table
    rows = []
    for i in range(len(manip)): # For each nutrient manipulation
        rows.append([*['']*4])
        df_manip = df[(df['Manipulation'] == i+1)]
        for j in range(len(data_filter)): # Containing given stoichiometric ratio
            df_filter = df_manip.dropna(subset=[data_filter[j]]) if j > 0 else df_manip
            culture = ['{:d} ({:d})'.format(len(df_filter[(df_filter['Growth.conditions'] == k+1)]),len(df_filter[(df_filter['Growth.conditions'] == k+1)]['Ref.'].unique())) for k in range(4)]
            rows.append(culture)
    #rows.append([len(df[(df['Growth.conditions'] == k+1)]) for k in range(4)])
    rows.append(['{:d} ({:d})'.format(len(df[(df['Growth.conditions'] == k+1)]),len(df[(df['Growth.conditions'] == k+1)]['Ref.'].unique())) for k in range(4)])

    # Compile data
    table = pd.DataFrame(rows,columns=col,index=idx).style
    tex = table.to_latex(caption='Number of data points (studies) by culture system in stoichiometric \\\\data compilation.',
                          column_format='@{}lcccc',
                          multicol_align="c",
                          position='H',
                          hrules=True).replace('& Batch','\\cmidrule{2-5} & Batch')
    tex = tex.replace("\\begin{tabular}", "\\small \\begin{tabular}")
    if save==True: saveTable(tex,'TableS5_StoichiometryDatabase_CultureType')
    #display(table)
    return None


####################################
# Run Tables
####################################

def runTables(save=True):
    """Produce and save tables in paper"""
    # Get model output/analysis
    df = ApplyModel.Model_df()
    df_optimalCN = ApplyModel.OptimalGrowth_CN()

    # Main text
    Table1_StoichDatabase(df,save)
    Table2_MeanStoich(df,save)

    # Supporting Information
    TableS1_StoichDatabasePhylum(df,save)
    TableS2_ModelConstants(save)
    TableS3_ModelFitOptions(save)
    TableS4_OptimalGrowthCN(df_optimalCN,save)
    TableS5_StoichDatabase_CultureType(df,save)
    return None
