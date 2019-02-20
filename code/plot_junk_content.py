#----------------------------------------------------------------------------------------
# This script creates plots for interface mapping description and fraction of junk PPIs
# for all structural interactomes together.
# Run script 'predict_edgotypes_and_junkPPIs.py' before running this script.
#----------------------------------------------------------------------------------------

import os
import pickle
import numpy as np
from pathlib import Path
from plot_tools import bar_plot

def main():
    
    # prediction methods for which results can be plotted
    # options: geometry, physics
    #pred_method = 'geometry'
    pred_method = 'physics'
    
    # method of calculating mutation ∆∆G for which results will be used in physics-based approach
    # options: 'bindprofx' or 'foldx'
    ddg_method = 'bindprofx'
    #ddg_method = 'foldx'
    
    # reference interactome names
    interactome_names = ['HI-II-14', 'IntAct', 'experiment']
    
    # structural interactome names for plot labels
    struc_interactome_names = ['Y2H-SI', 'IntAct-SI', 'experiment']
    
    # interactome colors for plot
    #interactome_colors = ['blue', 'green', 'red']
    interactome_colors = ['black', 'black', 'black']
    
    # plot confidence interval for the fraction of junk PPIs
    plotConfidenceIntervals = True
    
    # show figures
    showFigs = False
    
    dataDir = Path('/Volumes/MG_Samsung/junk_ppi_content/data')
    
    # directory to save processed data shared by all interactomes
    inDir = dataDir / 'processed'
    
    # figure directory
    figDir = Path('../figures') / 'combined'
    
    # create output directories if not existing
    if not figDir.exists():
        os.makedirs(figDir)
    
    if pred_method is 'physics':
        pred_method = '_'.join((pred_method, ddg_method))
    
    # input data files
    junkPPIFile_names = ('fraction_junk_PPIs_%s.pkl' % pred_method,
                         'fraction_junk_PPIs_%s.pkl' % pred_method,
                         'fraction_junk_PPIs_experiment.pkl')
    
    pN_E_results, pN_E_bounds = [], []
    for interactome_name, junkPPIFile_name in zip(interactome_names, junkPPIFile_names):
        junkPPIFile = inDir / interactome_name / junkPPIFile_name
        with open(junkPPIFile, 'rb') as f:
            all_results = pickle.load(f)
        pN_E_results.append( 100 * all_results['P(N|E)'] )
        if 'P(N|E)_CI' in all_results:
            lower, upper = all_results['P(N|E)_CI']
            pN_E_bounds.append( (100 * lower, 100 * upper) )
        else:
            pN_E_bounds.append( (0, 0) )
    
    if plotConfidenceIntervals:
        upper = [p + upper for p, (lower, upper) in zip(pN_E_results, pN_E_bounds)]
        maxY = max(upper)
    else:
        maxY = max(pN_E_results)
    maxY = 5 * np.ceil(maxY / 5)
    bar_plot(pN_E_results,
             pN_E_bounds if plotConfidenceIntervals else [0],
             xlabels = struc_interactome_names,
             ylabels = np.arange(0, maxY + 5, 5),
             ylabel = 'Fraction of junk PPIs (%)',
             colors = 'white',
             fmt = 'k.',
             capsize = 10 if plotConfidenceIntervals else 0,
             msize = 11,
             ewidth = 1.25,
             ecolors = interactome_colors,
             #ecolors = 'black',
             fontsize = 18,
             xlim = [0.8, 3.1],
             ylim = [0, maxY],
             yMinorTicks = True,
             adjustBottom = 0.2,
             shiftBottomAxis = -0.1,
             xbounds = (1,3),
             show = showFigs,
             figdir = figDir,
             figname = 'Fraction_junk_PPIs_%s' % pred_method)

if __name__ == "__main__":
    main()
