#----------------------------------------------------------------------------------------
# Plot tissue coexpression of interaction partners
#
# Run the following script before running this script:
# - compare_interactome_tissue_coexpr.py
#----------------------------------------------------------------------------------------

import os
import pickle
import numpy as np
from pathlib import Path
from stat_tools import sderror
from plot_tools import multi_bar_plot

def main():
    
    # reference interactome name
    # options: HI-II-14, IntAct
    ref_interactome_name = 'HI-II-14'
    
    # tissue expression databases
    expr_db = ['Illumina', 'GTEx', 'HPA', 'Fantom5']
    
    # structural interactome names
    struc_name = {'HI-II-14':'Y2H-SI', 'IntAct':'IntAct-SI'}
    
    # interactome names
    interactome_names = ['Random interactions', 'Reference interactome', 'Structural interactome']
    
    # interactome colors
    interactome_colors = ['limegreen', 'steelblue', 'orangered']
    
    # show figures
    showFigs = False
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # parent directory of processed data files
    procDir = dataDir / 'processed'
    
    # directory of processed data files specific to interactome
    interactomeDir = procDir / ref_interactome_name
    
    # figure directory
    figDir = Path('../figures') / ref_interactome_name / 'coexpr'
    
    # input data files
    coexprFiles = [interactomeDir / 'coexpr' / ('interactome_coexpr_%s.pkl' % db) for db in expr_db]
    
    # create output directories if not existing
    if not figDir.exists():
        os.makedirs(figDir)
    
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
                   errors = errors,
                   xlabels = expr_db,
                   ylabel = 'Tissue co-expression of\ninteraction partners',
                   ylabels = [round(x, 1) for x in np.arange(0, 0.9, 0.2)],
                   colors = interactome_colors,
                   edgecolor = 'k',
                   ewidth = 1.5,
                   barwidth = 0.2,
                   bargap = 0.03,
                   fontsize = 18,
                   capsize = 5,
                   msize = 8,
                   leg = ['Random pairs',
                          '%s reference interactome' % ref_interactome_name,
                          struc_name[ref_interactome_name]],
                   show = showFigs,
                   figdir = figDir,
                   figname = 'interactome_coexpr')
    
if __name__ == '__main__':
    main()
