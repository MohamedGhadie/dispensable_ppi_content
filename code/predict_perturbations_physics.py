#----------------------------------------------------------------------------------------
# This script constructs a structural interactome from a reference interactome by mapping 
# PPI binding interfaces, at atomic resolution, from experimentally determined 
# three-dimensional structural models in PDB onto PPIs in the reference interactome. 
# Next, Mendelian disease-causing mutations and common neutral mutations not associated 
# with disease are mapped onto the structural interactome, and PPI perturbations by 
# mutations are predicted. Mutation edgotypes ("edgetic" and "non-edgetic") per mutation 
# and per PPI are predicted using predicted PPI perturbations. Then, the fraction of junk 
# PPIs (PPIs neutral upon perturbation) in the structural interactome is estimated using 
# predicted mutation edgotypes, and compared to the fraction of junk PPIs estimated using 
# mutation edgotypes determined from experiments.
#
# Run script 'data_initial_processing.py' for preprocessing select data files before 
# running this script.
#----------------------------------------------------------------------------------------

import os
import pickle
import pandas as pd
import numpy as np
from pathlib import Path
from stat_tools import t_test, fisher_test, sderror, sderror_on_fraction
from plot_tools import bar_plot, multi_histogram_plot
from ddg_tools import read_protein_mutation_ddg
from mutation_interface_edgotype import energy_based_perturbation

def main():
    
    # reference interactome name. Options: 'HI-II-14' or 'IntAct'
    interactome_name = 'IntAct'
    
    # method of calculating mutation ∆∆G for which results will be used
    # options: 'bindprofx' or 'foldx'
    ddg_method = 'foldx'
    
    # Minimum reduction in binding free energy DDG required for interaction perturbation
    ddgCutoff = 0.5
    
    # calculate confidence interval for the fraction of junk PPIs
    computeConfidenceIntervals = True
    
    # % confidence interval
    CI = 95
    
    # show figures
    showFigs = False
    
    # directory of processed data files shared by all interactomes
    dataDir = Path('../data') / 'processed'
    
    # directory of processed data files specific to interactome
    interactomeDir = dataDir / interactome_name
    
    # directory to save output data files
    outDir = interactomeDir
    
    # figure directory
    figDir = Path('../figures') / interactome_name
    
    # create output directories if not existing
    if not interactomeDir.exists():
        os.makedirs(interactomeDir)
    if not outDir.exists():
        os.makedirs(outDir)
    if not figDir.exists():
        os.makedirs(figDir)
    
    # input data files
    geometryPerturbsFile = interactomeDir / 'unique_mutation_perturbs_geometry.pkl'
    natMutDDGinFile = interactomeDir / ('nondisease_mutations_%s_ddg.txt' % ddg_method)
    disMutDDGinFile = interactomeDir / ('disease_mutations_%s_ddg.txt' % ddg_method)
    
    # output data files
    natMutDDGoutFile = interactomeDir / ('nondisMut_%s_∆∆G_used.txt' % ddg_method)
    disMutDDGoutFile = interactomeDir / ('disMut_%s_∆∆G_used.txt' % ddg_method)
    physicsPerturbsFile = outDir / ('mutation_perturbs_physics_%s.pkl' % ddg_method)
    
    #------------------------------------------------------------------------------------
    # Fraction of mutation-targeted PPIs with ∆∆G exceeding a specified cutoff
    #------------------------------------------------------------------------------------
    
    # read change in binding free energy for interfacial mutations mapped on PDB chains
    naturalMutationsDDG = read_protein_mutation_ddg(natMutDDGinFile, 'binding')
    diseaseMutationsDDG = read_protein_mutation_ddg(disMutDDGinFile, 'binding')
    
    naturalMutations = pd.DataFrame(columns=["protein", "partner", "protein_pos", "mut_res", 
                                             "pdb_id", "chain_id", "chain_partner", "chain_mut", "ddg"])
    for i, item in enumerate(naturalMutationsDDG.items()):
        naturalMutations.loc[i] = item[0] + item[1]
    
    diseaseMutations = pd.DataFrame(columns=["protein", "partner", "protein_pos", "mut_res", 
                                             "pdb_id", "chain_id", "chain_partner", "chain_mut", "ddg"])
    for i, item in enumerate(diseaseMutationsDDG.items()):
        diseaseMutations.loc[i] = item[0] + item[1]
    
    naturalMutations.to_csv(natMutDDGoutFile, index=False, sep='\t')
    diseaseMutations.to_csv(disMutDDGoutFile, index=False, sep='\t')
    
    natMutDDGs = naturalMutations["ddg"].values
    disMutDDGs = diseaseMutations["ddg"].values
    
    numNatural_ddg_considered = len(natMutDDGs)
    numDisease_ddg_considered = len(disMutDDGs)
    
    print( '\n' + 'Avg. change in binding energy (∆∆G) for mutation-targeted PPIs:' )
    print( 'Non-disease: %.1f (SE = %g, n = %d)' % (np.mean(natMutDDGs),
                                                    sderror(natMutDDGs),
                                                    numNatural_ddg_considered) )
    
    print( 'Disease: %.1f (SE = %g, n = %d)' % (np.mean(disMutDDGs),
                                                sderror(disMutDDGs),
                                                numDisease_ddg_considered) )
    # Statistical significance of difference in means
    t_test(natMutDDGs, disMutDDGs)
    
    multi_histogram_plot ([disMutDDGs, natMutDDGs],
                          ['red', 'green'],
                          xlabel = 'Change in binding free energy (∆∆G)',
                          ylabel = 'Number of mutations',
                          leg = ['Disease interfacial mutations', 'Non-disease interfacial mutations'],
                          bins = 25,
                          alpha = 0.7,
                          fontsize = 24,
                          show = showFigs,
                          figdir = figDir,
                          figname = 'mut_ddg_histogram_%s' % ddg_method)
    
    numNatural_ddg = sum(natMutDDGs > ddgCutoff)
    numDisease_ddg = sum(disMutDDGs > ddgCutoff)
    fracNatural_ddg = numNatural_ddg / numNatural_ddg_considered
    fracDisease_ddg = numDisease_ddg / numDisease_ddg_considered
    fracNatural_ddg_error = sderror_on_fraction( numNatural_ddg, numNatural_ddg_considered )
    fracDisease_ddg_error = sderror_on_fraction( numDisease_ddg, numDisease_ddg_considered )
    
    print( '\n' + 'Fraction of mutation-targeted PPIs with ∆∆G > %.1f kcal/mol:' % ddgCutoff )
    print( 'Non-disease: %.3f (SE = %g, ntot = %d)' % (fracNatural_ddg,
                                                       fracNatural_ddg_error,
                                                       numNatural_ddg_considered) )
    
    print( 'Disease: %.3f (SE = %g, ntot = %d)' % (fracDisease_ddg,
                                                   fracDisease_ddg_error,
                                                   numDisease_ddg_considered) )
    
    # Statistical significance of difference in fractions
    fisher_test([numNatural_ddg, numNatural_ddg_considered - numNatural_ddg],
                [numDisease_ddg, numDisease_ddg_considered - numDisease_ddg])
    
    bar_plot([fracNatural_ddg, fracDisease_ddg],
             error = [fracNatural_ddg_error, fracDisease_ddg_error],
             xlabels = ('PPIs with\nnon-disease\nmutations\nat interface',
                        'PPIs with\ndisease\nmutations\nat interface'),
             ylabel = ('Fraction with ∆∆G > %.1f kcal/mol' % ddgCutoff),
             ylabels = [0, 0.2, 0.4, 0.6, 0.8],
             ylim = [0, 0.8],
             colors = ['turquoise', 'orangered'],
             barwidth = 0.6,
             fontsize = 24,
             show = showFigs,
             figdir = figDir,
             figname = 'mut_ddg_frac_>%.1f_%s' % (ddgCutoff, ddg_method))
    
    #------------------------------------------------------------------------------------
    # predict PPI perturbations
    #------------------------------------------------------------------------------------

    if geometryPerturbsFile.is_file():
        print( '\n' + 'Loading geometry-based PPI perturbation predictions' )
        with open(geometryPerturbsFile, 'rb') as f:
            naturalPerturbs, diseasePerturbs = pickle.load(f)
    else:
        print( '\n' + 'Geometry-based PPI perturbation prediction file not found' )
        return

    print( '\n' + 'Performing physics-based edgotype prediction for non-disease mutations' )
    naturalPerturbs["perturbations"], knownDDG, unknownDDG = energy_based_perturbation (naturalPerturbs,
                                                                                        naturalMutationsDDG,
                                                                                        ddgCutoff)
    print( '\n' + 'Performing physics-based edgotype prediction for disease mutations' )
    diseasePerturbs["perturbations"], knownDDG, unknownDDG = energy_based_perturbation (diseasePerturbs,
                                                                                        diseaseMutationsDDG,
                                                                                        ddgCutoff)
    with open(physicsPerturbsFile, 'wb') as fOut:
        pickle.dump([naturalPerturbs, diseasePerturbs], fOut)

if __name__ == "__main__":
    main()
