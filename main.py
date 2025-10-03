"""
Run data compilation/model analysis and produce all figures and tables in study.
Adjust inputs to main() to determine what output to save. Default is to save all
figures and tables to project folder. Alternatively, can provide main() with
optional input parameters:
"saveFigures=False" can be used to display output instead of saving
"main_figs=False" or "SI_figs=False" will suppress figures from main or supplement
"saveTables=False" will stop tables from being saved to project folder
"""

import ApplyModel
import Figures
import Tables


def main(saveFigures=True,main_figs=True,SI_figs=True,saveTables=True):
    """Run model using calculated parameters and produce tables/figures"""
    # Get compiled data and apply macromolecular model (multiple formulations)
    df = ApplyModel.Model_df()

    # Produce tables and figures
    Figures.runFigures(main_text=main_figs,SI=SI_figs,save=saveFigures)
    Tables.runTables(save=saveTables)

    return None

# Default: run/save all figures and tables
# Optional: set input(s) to False to suppress output (see above)
main()