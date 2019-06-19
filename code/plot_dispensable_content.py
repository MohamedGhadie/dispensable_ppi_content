#----------------------------------------------------------------------------------------
# Plot the fraction of dispensable PPIs calculated from predictions and experiments.
#
# Run the following scripts before running this script:
# - calculate_dispensable_content_geometry.py
# - calculate_dispensable_content_physics.py
# - calculate_dispensable_content_experiment.py
#----------------------------------------------------------------------------------------

import os
import pickle
import numpy as np
from pathlib import Path
from plot_tools import curve_plot

def main():
    
    # prediction methods for which results can be plotted
    # options: geometry, physics
    pred_method = 'geometry'
    
    # method of calculating mutation ∆∆G for which results will be used in physics-based approach
    # options: 'bindprofx' or 'foldx'
    ddg_method = 'bindprofx'
    
    # set to True to plot dispensable PPI content using fraction of mono-edgetic mutations 
    # instead of edgetic mutations
    mono_edgetic = False
    
    # reference interactome names
    interactome_names = ['HI-II-14', 'IntAct', 'experiment']
    
    # structural interactome names for plot labels
    struc_interactome_names = ['Y2H-SI', 'IntAct-SI', 'Experiment']
    
    # plot confidence interval for the fraction of dispensable PPIs
    plotConfidenceIntervals = True
    
    # show figures
    showFigs = False
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # parent directory of processed data files
    procDir = dataDir / 'processed'
        
    # figure directory
    figDir = Path('../figures') / 'combined'
    
    if pred_method is 'physics':
        pred_method = '_'.join((pred_method, ddg_method))
    if mono_edgetic:
        pred_method = pred_method + '_monoedgetic'
    
    # input data files
    dispensablePPIFile_names = ('fraction_disp_PPIs_%s.pkl' % pred_method,
                                'fraction_disp_PPIs_%s.pkl' % pred_method,
                                'fraction_disp_PPIs_experiment.pkl')
    
    # create output directories if not existing
    if not figDir.exists():
        os.makedirs(figDir)
    
    pN_E_results, pN_E_bounds = [], []
    for interactome_name, dispensablePPIFile_name in zip(interactome_names, dispensablePPIFile_names):
        dispensablePPIFile = procDir / interactome_name / dispensablePPIFile_name
        with open(dispensablePPIFile, 'rb') as f:
            all_results = pickle.load(f)
        pN_E_results.append( 100 * all_results['P(N|E)'] )
        if 'P(N|E)_CI' in all_results:
            lower, upper = all_results['P(N|E)_CI']
            pN_E_bounds.append( (100 * lower, 100 * upper) )
        else:
            pN_E_bounds.append( (0, 0) )
    
    numInteractomes = len(interactome_names)
    if plotConfidenceIntervals:
        upper = [p + upper for p, (lower, upper) in zip(pN_E_results, pN_E_bounds)]
        maxY = max(upper)
    else:
        maxY = max(pN_E_results)
    maxY = 5 * np.ceil(maxY / 5)
    curve_plot (pN_E_results,
                error = pN_E_bounds if plotConfidenceIntervals else None,
                xlim = [0.8, numInteractomes + 0.1],
                ylim = [0, maxY],
                styles = '.k',
                capsize = 10 if plotConfidenceIntervals else 0,
                msize = 16,
                ewidth = 1.25,
                ecolors = 'k',
                ylabel = 'Fraction of dispensable PPIs (%)',
                yMinorTicks = 4,
                xticks = list(np.arange(1, numInteractomes + 1)),
                xticklabels = struc_interactome_names,
                yticklabels = list(np.arange(0, maxY + 5, 5)),
                fontsize = 20,
                adjustBottom = 0.2,
                shiftBottomAxis = -0.1,
                xbounds = (1, numInteractomes - 1) if mono_edgetic else (1, numInteractomes),
                show = showFigs,
                figdir = figDir,
                figname = 'Fraction_disp_PPIs_%s' % pred_method)

if __name__ == "__main__":
    main()
