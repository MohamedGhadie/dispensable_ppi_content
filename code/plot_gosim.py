#----------------------------------------------------------------------------------------
# Plot gene ontology (GO) similarity of interaction partners
#
# Run the following script before running this script:
# - compare_interactome_gosim.py
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
    
    # similarity measure used to calculate GO similarity
    sim_measure = 'SimGIC'
    
    # root ontologies on which GO similarity was calculated
    ont_root = ['biological_process', 'molecular_function', 'cellular_component']
    
    # root ontology labels used to label output files and figures
    ont_abv = {'biological_process':'bp', 'molecular_function':'mf', 'cellular_component':'cc'}
    
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
    figDir = Path('../figures') / ref_interactome_name / 'gosim'
    
    # input data files
    gosimFiles = [interactomeDir / 'gosim' / ('allPPI_gosim_%s_%s.pkl' % (ont_abv[root], sim_measure)) for root in ont_root]
    
    # create output directories if not existing
    if not figDir.exists():
        os.makedirs(figDir)
    
    allgosim = {}
    for root, gosimFile in zip(ont_root, gosimFiles):
        with open(gosimFile, 'rb') as f:
            allgosim[root] = pickle.load(f)
    
    means, errors = [], []
    for interactome_name in interactome_names:
        gosim = [allgosim[root][interactome_name]["gosim"] for root in ont_root]
        means.append([np.mean(s) for s in gosim])
        errors.append([sderror(s) for s in gosim])
    
    multi_bar_plot(means,
                   errors = errors,
                   xlabels = [r.title().replace('_','\n') for r in ont_root],
                   ylabel = 'Gene ontology similarity\nof interaction partners',
                   ylabels = [round(x, 1) for x in np.arange(0, 0.6, 0.1)],
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
                   figname = 'interactome_gosim_%s' % sim_measure)
    
if __name__ == '__main__':
    main()
