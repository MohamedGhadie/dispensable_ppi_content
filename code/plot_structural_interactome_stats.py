#----------------------------------------------------------------------------------------
# Plot properties of the structural interactomes combined.
#
# Run the following scripts for each interactome before running this script:
# - produce_data_mappings.py
# - produce_structural_interactome.py
# - structural_interactome_stats.py
#----------------------------------------------------------------------------------------

import os
import pickle
import numpy as np
from pathlib import Path
from plot_tools import multi_bar_plot, multi_histogram_plot

def main():
    
    # reference interactome names
    interactome_names = ['HI-II-14', 'IntAct']
    
    # structural interactome names
    struc_interactome_names = ['Y2H-SI', 'IntAct-SI']
    
    # interactome colors
    interactome_colors = ['darkviolet', 'gray']
    
    # show figures
    showFigs = False
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # parent directory of processed data files
    procDir = dataDir / 'processed'
    
    # figure directory
    figDir = Path('../figures') / 'combined'
    
    # input data files
    proteinSeqFile = procDir / 'human_reference_sequences.pkl'
    numPPIsPerPDBFiles = [procDir / name / 'numPPIsPerPDB.pkl' for name in interactome_names]
    interfaceFiles = [procDir / name / 'proteinInterfaces.pkl' for name in interactome_names]
    
    # create output directories if not existing
    if not figDir.exists():
        os.makedirs(figDir)
    
    #------------------------------------------------------------------------------------
    # plot distribution of the number of PPIs modeled per PDB structure
    #------------------------------------------------------------------------------------
    
    maxVal = 10
    distribution = []
    for inPath in numPPIsPerPDBFiles:
        with open(inPath, 'rb') as f:
            numPPIsPerPDB = pickle.load(f)
        numPPIsPerPDB_capped = [min(n, maxVal) for n in numPPIsPerPDB]
        distribution.append([numPPIsPerPDB_capped.count(i) for i in np.arange(1, maxVal + 1)])
    
    multi_bar_plot(distribution,
                   xlabels = list(map(str, np.arange(1, maxVal))) + [ str(maxVal) + '+'],
                   xlabel = 'Number of PPIs modeled by PDB structure',
                   ylabel = 'Number of PDB structures',
                   ylabels = [0, 500, 1000, 1500, 2000],
                   ylim = [0, 2000],
                   colors = interactome_colors,
                   edgecolor = 'k',
                   barwidth = 0.4,
                   fontsize = 22,
                   leg = struc_interactome_names,
                   show = showFigs,
                   figdir = figDir,
                   figname = 'numPPIsPerPDB')
    
    #------------------------------------------------------------------------------------
    # plot distribution of protein interface length and relative length
    #------------------------------------------------------------------------------------
    
    with open(proteinSeqFile, 'rb') as f:
        proteinSeq = pickle.load(f)
    
    interfaceLen, interfaceRatio = [], []
    for inPath in interfaceFiles:
        with open(inPath, 'rb') as f:
            interfaces = pickle.load(f)
        interfaceLen.append([len(i) for i in interfaces.values()])
        interfaceRatio.append([len(i)/len(proteinSeq[p]) for (p, _), i in interfaces.items()])
    
    multi_histogram_plot (list(reversed(interfaceLen)),
                          colors = list(reversed(interactome_colors)),
                          xlabel = 'Number of interfacial residues',
                          ylabel = 'Number of proteins',
                          ylabels = [0, 250, 500, 750, 1000, 1250],
                          ylim = [0, 1250],
                          edgecolor = 'k',
                          fontsize = 22,
                          bins = 25,
                          alpha = 1,
                          leg = struc_interactome_names,
                          show = showFigs,
                          figdir = figDir,
                          figname = 'interface_length')
    
    multi_histogram_plot (list(reversed(interfaceRatio)),
                          colors = list(reversed(interactome_colors)),
                          xlabel = 'Fraction of interfacial residues',
                          ylabel = 'Number of proteins',
                          edgecolor = 'k',
                          fontsize = 24,
                          bins = 25,
                          alpha = 1,
                          leg = struc_interactome_names,
                          show = showFigs,
                          figdir = figDir,
                          figname = 'interface_frac')

if __name__ == '__main__':
    main()
