#----------------------------------------------------------------------------------------
# Plot the fraction of junk PPIs calculated from predictions and experiments
#----------------------------------------------------------------------------------------

import os
import pickle
import numpy as np
from pathlib import Path
from stat_tools import sderror
from plot_tools import multi_bar_plot

def main():
    
    # tissue expression databases
    expr_db = ['Illumina', 'GTEx', 'HPA', 'Fantom5']
    
    # reference interactome name
    ref_interactome_name = 'IntAct'
    
    # interactome names
    interactome_names = ['Random interactions', 'Reference interactome', 'Structural interactome']
    
    # interactome colors
    interactome_colors = ['green', 'blue', 'red']
    
    # show figures
    showFigs = False
    
    # directory of processed data files shared by all interactomes
    dataDir = Path('../data') / 'processed'
    
    # directory of processed data files specific to interactome
    interactomeDir = dataDir / ref_interactome_name
    
    # figure directory
    figDir = Path('../figures') / ref_interactome_name
    
    # create output directories if not existing
    if not figDir.exists():
        os.makedirs(figDir)
    
    # input data files
    coexprFiles = [interactomeDir / ('interactome_coexpr_%s.pkl' % db) for db in expr_db]
                         
    allcoexpr = {}
    for db, coexprFile in zip(expr_db, coexprFiles):
        with open(coexprFile, 'rb') as f:
            allcoexpr[db] = pickle.load(f)
    
    means, errors = [], []
    for interactome_name in interactome_names:
        coexpr = [allcoexpr[db][interactome_name]["coexpr"] for db in expr_db]
        means.append([np.mean(c) for c in coexpr])
        errors.append([sderror(c) for c in coexpr])
    
    multi_bar_plot(means,
                   errors,
                   xlabels = expr_db,
                   ylabel = 'Tissue co-expression of\ninteraction partners',
                   ylabels = [round(x, 1) for x in np.arange(0, 0.9, 0.2)],
                   colors = interactome_colors,
                   edgecolor = 'k',
                   barwidth = 0.2,
                   fontsize = 20,
                   leg = interactome_names,
                   show = showFigs,
                   figdir = figDir,
                   figname = 'interactome_coexpr')
    
if __name__ == '__main__':
    main()
