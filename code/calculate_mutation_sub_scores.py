#------------------------------------------------------------------------------------
# Calculate substitution score for all mutations
#------------------------------------------------------------------------------------
    
import os
import pickle
import numpy as np
from pathlib import Path
from mutation_processing_tools import remove_mutation_overlaps
from id_mapping import produce_substitution_matrix
from stat_tools import sderror, bootstrap_test
from plot_tools import box_plot

def main():
    
    # substitution matrix name
    matrixName = 'PAM30'
    
    # show figures
    showFigs = False
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # figure directory
    figDir = Path('../figures') / 'combined'
    
    # input data files
    naturalMutationsFile = procDir / 'dbsnp_mutations4.txt'
    diseaseMutationsFile = procDir / 'clinvar_mutations6.txt'
    
    # output data files
    subsMatrixFile = procDir / (matrixName + '.pkl')
    
    # create output directories if not existing
    if not procDir.exists():
        os.makedirs(procDir)
    if not figDir.exists():
        os.makedirs(figDir)
    
    #------------------------------------------------------------------------------------
    # produce substitution matrix
    #------------------------------------------------------------------------------------
    
    if not subsMatrixFile.is_file():
        produce_substitution_matrix (matrixName, subsMatrixFile)
    with open(subsMatrixFile, 'rb') as f:
        subsTable = pickle.load(f)
    
    #------------------------------------------------------------------------------------
    # load  mutations
    #------------------------------------------------------------------------------------
    
    naturalMutations, diseaseMutations = remove_mutation_overlaps (naturalMutationsFile, diseaseMutationsFile)
    
    #------------------------------------------------------------------------------------
    # calculate mutation substitution scores
    #------------------------------------------------------------------------------------
    
    natMutScore = [subsTable[x] for x in zip(naturalMutations["wt_res"], naturalMutations["mut_res"])]
    disMutScore = [subsTable[x] for x in zip(diseaseMutations["wt_res"], diseaseMutations["mut_res"])]
    
    print()
    print( 'Avg. score for natural mutations: %.1f (SE = %g, n = %d)' % (np.mean(natMutScore),
                                                                         sderror(natMutScore),
                                                                         len(natMutScore)) )
    print( 'Avg. score for disease mutations: %.1f (SE = %g, n = %d)' % (np.mean(disMutScore),
                                                                         sderror(disMutScore),
                                                                         len(disMutScore)) )
    
    bootstrap_test (natMutScore, disMutScore, 10000)
    
    box_plot([natMutScore, disMutScore],
             xlabels = ('Non-disease\nmutations','Disease\nmutations'),
             ylabels = [-15, -10, -5, 0, 5],
             ylabel = 'PAM30 substitution score',
             ylim = [-16, 5],
             colors = ['turquoise', 'magenta'],
             fontsize = 26,
             show = showFigs,
             figdir = figDir,
             figname = 'substitution_score')

if __name__ == '__main__':
    main()
