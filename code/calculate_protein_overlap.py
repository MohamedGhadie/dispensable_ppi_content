#----------------------------------------------------------------------------------------
# Calculate overlap in proteins covered by predictions and experiments
#----------------------------------------------------------------------------------------

import os
import pickle
import pandas as pd
from pathlib import Path
from plot_tools import venn2_plot

def main():
    
    # reference interactome name. Options: 'HI-II-14' or 'IntAct'
    interactome_name = 'HI-II-14'
    
    # show figures
    showFigs = False
    
    # interactome colors
    interactome_color = {'HI-II-14':'royalblue', 'IntAct':'limegreen'}
    
    # directory of processed data files shared by all interactomes
    dataDir = Path('../data') / 'processed'
    
    # directory of processed data files specific to interactome
    interactomeDir = dataDir / interactome_name
    
    # directory of processed data files specific to experiment
    experimentDir = dataDir / 'experiment'
    
    # figure directory
    figDir = Path('../figures') / interactome_name
    
    # create output directories if not existing
    if not figDir.exists():
        os.makedirs(figDir)
    
    # input data files
    natMutPredEdgotypeFile = interactomeDir / 'nondisease_mutation_edgotype_geometry.txt'
    disMutPredEdgotypeFile = interactomeDir / 'disease_mutation_edgotype_geometry.txt'
    natMutExpEdgotypeFile = experimentDir / 'nondisease_mutation_edgotype_experiment.txt'
    disMutExpEdgotypeFile = experimentDir / 'disease_mutation_edgotype_experiment.txt'
    UniprotIDmapFile = dataDir / 'to_human_uniprotID_map.pkl'
    
    #----------------------------------------------------------------------------------------
    # process proteins covered by predictions
    #----------------------------------------------------------------------------------------

    naturalMutations = pd.read_table (natMutPredEdgotypeFile, sep='\t')
    diseaseMutations = pd.read_table (disMutPredEdgotypeFile, sep='\t')
    
    natMutProteins = set(naturalMutations["protein"].values)
    disMutProteins = set(diseaseMutations["protein"].values)
    allMutProteins = natMutProteins | disMutProteins
    
    print( '\n' + 'Prediction mutations:' )
    print('Non-disease: %d mutations on %d proteins' % (len(naturalMutations), len(natMutProteins)))
    print('Disease: %d mutations on %d proteins' % (len(diseaseMutations), len(disMutProteins)))
    print('Total mutated proteins: %d' % len(allMutProteins))
    
    #----------------------------------------------------------------------------------------
    # process proteins covered by experiments
    #----------------------------------------------------------------------------------------
    
    expNaturalMutations = pd.read_table (natMutExpEdgotypeFile, sep='\t')
    expDiseaseMutations = pd.read_table (disMutExpEdgotypeFile, sep='\t')
    
    with open(UniprotIDmapFile, 'rb') as f:
        uniprotID = pickle.load(f)
    expNaturalMutations["protein"] = expNaturalMutations["Symbol"].apply(lambda x:
                                                                         uniprotID[x] if x in uniprotID 
                                                                         else x + 'unknown')
    expDiseaseMutations["protein"] = expDiseaseMutations["Symbol"].apply(lambda x:
                                                                         uniprotID[x] if x in uniprotID 
                                                                         else x + 'unknown')
    
    natMutProteins_exp = set(expNaturalMutations["protein"].values)
    disMutProteins_exp = set(expDiseaseMutations["protein"].values)
    allMutProteins_exp = natMutProteins_exp | disMutProteins_exp
    
    print( '\n' + 'Experiment mutations:' )
    print('Non-disease: %d mutations on %d proteins' % (len(expNaturalMutations), len(natMutProteins_exp)))
    print('Disease: %d mutations on %d proteins' % (len(expDiseaseMutations), len(disMutProteins_exp)))
    print('Total mutated proteins: %d' % len(allMutProteins_exp))
    
    #----------------------------------------------------------------------------------------
    # calculate overlap in protein space
    #----------------------------------------------------------------------------------------
    
    natMutProteins_overlap = natMutProteins & natMutProteins_exp
    disMutProteins_overlap = disMutProteins & disMutProteins_exp
    allMutProteins_overlap = allMutProteins & allMutProteins_exp
    
    print( '\n' + 'Overlap in proteins covered by all mutations:' )
    print( 'prediction proteins: %.1f %% (%d out of %d)' % ( 100 * len(allMutProteins_overlap) 
                                                                 / len(allMutProteins),
                                                             len(allMutProteins_overlap),
                                                             len(allMutProteins) ) )
    print( 'experiment proteins: %.1f %% (%d out of %d)' % ( 100 * len(allMutProteins_overlap) 
                                                                 / len(allMutProteins_exp),
                                                             len(allMutProteins_overlap),
                                                             len(allMutProteins_exp) ) )
    
    print( '\n' + 'Overlap in proteins covered by non-disease mutations:' )
    print( 'prediction proteins: %.1f %% (%d out of %d)' % ( 100 * len(natMutProteins_overlap) 
                                                                 / len(natMutProteins),
                                                             len(natMutProteins_overlap),
                                                             len(natMutProteins) ) )
    print( 'experiment proteins: %.1f %% (%d out of %d)' % ( 100 * len(natMutProteins_overlap) 
                                                                 / len(natMutProteins_exp),
                                                             len(natMutProteins_overlap),
                                                             len(natMutProteins_exp) ) )
    
    print( '\n' + 'Overlap in proteins covered by disease mutations:' )
    print( 'prediction proteins: %.1f %% (%d out of %d)' % ( 100 * len(disMutProteins_overlap) 
                                                                 / len(disMutProteins),
                                                             len(disMutProteins_overlap),
                                                             len(disMutProteins) ) )
    print( 'experiment proteins: %.1f %% (%d out of %d)' % ( 100 * len(disMutProteins_overlap) 
                                                                 / len(disMutProteins_exp),
                                                             len(disMutProteins_overlap),
                                                             len(disMutProteins_exp) ) )
    
    # plot overlap in proteins covered by nondisease mutations
    venn2_plot([natMutProteins, natMutProteins_exp],
               #labels = ['Prediction', 'Experiment'],
               colors = [interactome_color[interactome_name], 'orangered'],
               hideNum = True,
               show = showFigs,
               figdir = figDir,
               figname = 'non_disease_protein_overlap')
    # plot overlap in proteins covered by disease mutations
    venn2_plot([disMutProteins, disMutProteins_exp],
               #labels = ['Prediction', 'Experiment'],
               colors = [interactome_color[interactome_name], 'orangered'],
               hideNum = True,
               show = showFigs,
               figdir = figDir,
               figname = 'disease_protein_overlap')
    # plot legend circles
    venn2_plot([{1}, {2}],
               #labels = ['Prediction', 'Experiment'],
               colors = [interactome_color[interactome_name], 'orangered'],
               hideNum = True,
               show = showFigs,
               figdir = figDir,
               figname = 'legend_overlap')

if __name__ == '__main__':
    main()
