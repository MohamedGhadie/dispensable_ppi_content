#----------------------------------------------------------------------------------------
# Calculate the fraction of junk protein-protein interactions (PPIs) in the human 
# interactome, i.e., the fraction of PPIs that are effectively neutral upon perturbation.
# The fraction of junk PPIs is calculated from experimentally determined PPI perturbations.
#----------------------------------------------------------------------------------------

import os
import pickle
import pandas as pd
import numpy as np
from pathlib import Path
from text_tools import write_list_table
from stat_tools import fisher_test, sderror_on_fraction, proportion_ratio_CI
from plot_tools import pie_plot

def main():
    
    # reference interactome name
    interactome_name = 'experiment'
    
    # treat experiment quasi-null mutation as edgetic
    consider_exp_QuasiNull_perturbs = False
    
    # calculate confidence interval for the fraction of junk PPIs
    computeConfidenceIntervals = True
    
    # % confidence interval
    CI = 95
    
    # show figures
    showFigs = False
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # directory of data files from external sources
    extDir = Path('/Volumes/MG_Samsung/junk_ppi_content/data/external')
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # directory of processed data files specific to interactome
    interactomeDir = procDir / interactome_name
    
    # figure directory
    figDir = Path('../figures') / interactome_name
    
    if consider_exp_QuasiNull_perturbs:
        interactomeDir = interactomeDir / 'quasiNull_perturbs_considered'
        figDir = figDir / 'quasiNull_perturbs_considered'
    
    # input data files
    edgotypeFile = extDir / "Sahni_2015_Table_S3.xlsx"
    
    # output data files
    natMutEdgotypeFile = interactomeDir / 'nondisease_mutation_edgotype_experiment.txt'
    disMutEdgotypeFile = interactomeDir / 'disease_mutation_edgotype_experiment.txt'
    junkPPIFile = interactomeDir / 'fraction_junk_PPIs_experiment.pkl'
    
    # create output directories if not existing
    if not procDir.exists():
        os.makedirs(procDir)
    if not interactomeDir.exists():
        os.makedirs(interactomeDir)
    if not figDir.exists():
        os.makedirs(figDir)
    
    #------------------------------------------------------------------------------------
    # process experimental edgotype data
    #------------------------------------------------------------------------------------
    
    expMutations = pd.read_excel(edgotypeFile, sheet_name = 'Table S3C')
    expMutations_ppi = pd.read_excel(edgotypeFile, sheet_name = 'Table S3A')
    
    expMutations = expMutations[ expMutations["Edgotype_class"]
                                 .apply(lambda x: x in {'Edgetic',
                                                        'Quasi-null',
                                                        'Quasi-wild-type'})
                               ].reset_index( drop = True )
    
    expMutations["partners"] = expMutations["Allele_ID"].apply(
        lambda x: expMutations_ppi.loc[ expMutations_ppi["Allele_ID"] == x,
                                        "Interactor_Gene_ID" ].values )
    expMutations["Mut_interaction"] = expMutations["Allele_ID"].apply(
        lambda x: expMutations_ppi.loc[ expMutations_ppi["Allele_ID"] == x,
                                        "Y2H_score" ].values )
    
    expMutations["WT_interaction"] = expMutations.apply(
        lambda x: np.array( [ expMutations_ppi.loc[ (expMutations_ppi["Category"] == 'Wild-type') & 
                                                    (expMutations_ppi["Entrez_Gene_ID"] == x["Entrez_Gene_ID"]) & 
                                                    (expMutations_ppi["Interactor_Gene_ID"] == p),
                                                    "Y2H_score"].item()
                              for p in x["partners"] ] ),
        axis=1 )
    
    expMutations["perturbations"] = expMutations.apply(
        lambda x:  0 + ( ( x["WT_interaction"] == 1 ) & ( x["Mut_interaction"] == 0 ) ), axis=1 )
    
    expNaturalMutations = expMutations[ expMutations["Category"]
                                        == 'Non-disease variant' ].reset_index( drop = True )
    expDiseaseMutations = expMutations[ expMutations["Category"]
                                        == 'Disease mutation' ].reset_index( drop = True )
    
    # write mutation edgotypes to tab-delimited file
    write_list_table (expNaturalMutations, ["partners", "Mut_interaction", "WT_interaction", "perturbations"], natMutEdgotypeFile)
    write_list_table (expDiseaseMutations, ["partners", "Mut_interaction", "WT_interaction", "perturbations"], disMutEdgotypeFile)
    
    #------------------------------------------------------------------------------------
    # Fraction of edgetic mutations in experiments
    #------------------------------------------------------------------------------------
    
    if consider_exp_QuasiNull_perturbs:
        # Number of edgetic (E) and quasi-null (Q) mutations in experiment
        numNaturalMut_edgetic = sum( (expNaturalMutations["Edgotype_class"] == 'Edgetic') |
                                         (expNaturalMutations["Edgotype_class"] == 'Quasi-null') )
        numDiseaseMut_edgetic = sum( (expDiseaseMutations["Edgotype_class"] == 'Edgetic') |
                                         (expDiseaseMutations["Edgotype_class"] == 'Quasi-null') )
        
        # Number of non-edgetic mutations (quasi-wild-type) in experiment
        numNaturalMut_nonedgetic = sum( expNaturalMutations["Edgotype_class"] == 'Quasi-wild-type' )
        numDiseaseMut_nonedgetic = sum( expDiseaseMutations["Edgotype_class"] == 'Quasi-wild-type' )
    else:
        # Number of edgetic (E) mutations in experiment
        numNaturalMut_edgetic = sum( expNaturalMutations["Edgotype_class"] == 'Edgetic' )
        numDiseaseMut_edgetic = sum( expDiseaseMutations["Edgotype_class"] == 'Edgetic' )
        
        # Number of non-edgetic mutations (quasi-wild-type or quasi-null) in experiment
        numNaturalMut_nonedgetic = sum( (expNaturalMutations["Edgotype_class"] == 'Quasi-null') |
                                            (expNaturalMutations["Edgotype_class"] == 'Quasi-wild-type') )
        numDiseaseMut_nonedgetic = sum( (expDiseaseMutations["Edgotype_class"] == 'Quasi-null') |
                                            (expDiseaseMutations["Edgotype_class"] == 'Quasi-wild-type') )    
    # Total number of mutations tested in experiment
    numNaturalMut_considered = numNaturalMut_edgetic + numNaturalMut_nonedgetic
    numDiseaseMut_considered = numDiseaseMut_edgetic + numDiseaseMut_nonedgetic
    
    # Fraction of mutations that are edgetic (E) in experiment
    fracNaturalMut_edgetic = numNaturalMut_edgetic / numNaturalMut_considered
    fracDiseaseMut_edgetic = numDiseaseMut_edgetic / numDiseaseMut_considered
    
    # Standard error of the fraction of edgetic mutations in experiment
    fracNaturalMut_error = sderror_on_fraction(numNaturalMut_edgetic, numNaturalMut_considered)
    fracDiseaseMut_error = sderror_on_fraction(numDiseaseMut_edgetic, numDiseaseMut_considered)
    
    print( '\n' + 'Fraction of experimental edgetic mutations:' )
    print('Non-disease mutations: %.3f (SE = %g, %d out of %d)' % (fracNaturalMut_edgetic,
                                                                   fracNaturalMut_error,
                                                                   numNaturalMut_edgetic,
                                                                   numNaturalMut_considered))
    
    print('Disease mutations: %.3f (SE = %g, %d out of %d)' % (fracDiseaseMut_edgetic,
                                                               fracDiseaseMut_error,
                                                               numDiseaseMut_edgetic,
                                                               numDiseaseMut_considered))
    
    fisher_test([numNaturalMut_edgetic, numNaturalMut_nonedgetic],
                [numDiseaseMut_edgetic, numDiseaseMut_nonedgetic])
    
    pie_plot([numNaturalMut_nonedgetic, numNaturalMut_edgetic],
             angle = 90,
             colors = ['mediumslateblue', 'red'],
             edgewidth = 2,
             show = showFigs,
             figdir = figDir,
             figname = 'non_disease_edgetic_mutations_exp')
    pie_plot([numDiseaseMut_nonedgetic, numDiseaseMut_edgetic],
             angle = 90,
             colors = ['mediumslateblue', 'red'],
             edgewidth = 2,
             show = showFigs,
             figdir = figDir,
             figname = 'disease_edgetic_mutations_exp')
    
    #------------------------------------------------------------------------------------
    # apply Bayes' theorem to calculate the fraction of PPIs that are junk, i.e., 
    # effectively neutral under perturbation
    #------------------------------------------------------------------------------------
    
    # Probability for new missense mutations to be neutral (N)
    pN = 0.27

    # Probability for new missense mutations to be mildly deleterious (M)
    pM = 0.53

    # Probability for new missense mutations to be strongly detrimental (S)
    pS = 0.20
    
    # Probability for neutral mutations (N) to be edgetic (E)
    pE_N = fracNaturalMut_edgetic
    
    # Probability for mildly deleterious mutations (M) to be edgetic (E)
    pE_M = fracDiseaseMut_edgetic
    
    # Probability for strongly detrimental mutations (S) to be edgetic (E)
    pE_S = 0
        
    # Probability for a new missense mutation to be edgetic
    pE = (pE_N * pN) + (pE_M * pM) + (pE_S * pS)
    
    # Probability for edgetic mutations to be effectively neutral
    pN_E = pE_N * pN / pE
    
    allresults = {}
    allresults['P(N|E)'] = pN_E
    
    print('')
    print( 'P(N) = %.1f %%' % (100 * pN) )
    print( 'P(M) = %.1f %%' % (100 * pM) )
    print( 'P(S) = %.1f %%' % (100 * pS) )
    print( 'P(E|N) = %.1f %%' % (100 * pE_N) )
    print( 'P(E|M) = %.1f %%' % (100 * pE_M) )
    print( 'P(E|S) = %.1f %%' % (100 * pE_S) )
    print( 'P(E) = P(E|N)*P(N) + P(E|M)*P(M) + P(E|S)*P(S) = %.1f %%' % (100 * pE) )
    print( 'Fraction of junk PPIs P(N|E) = P(E|N)*P(N)/P(E) = %.1f %%' % (100 * pN_E) )
    
    # calculate 95% confidence interval
    if computeConfidenceIntervals:
        n_N, n_M = numNaturalMut_considered, numDiseaseMut_considered
        k_obs_N, k_obs_M = numNaturalMut_edgetic, numDiseaseMut_edgetic
        pE_M_pE_N_lower, pE_M_pE_N_upper = proportion_ratio_CI (k_obs_M,
                                                                n_M,
                                                                k_obs_N,
                                                                n_N,
                                                                conf = CI)
        pN_E_lower = 1 / ( pE_M_pE_N_upper * (pM / pN) + 1 )
        pN_E_upper = 1 / ( pE_M_pE_N_lower * (pM / pN) + 1 )
        print( '%.1f%% confidence interval for P(N|E) = (%f, %f)' 
                % (CI, 100 * pN_E_lower, 100 * pN_E_upper) )
        if not (np.isnan(pN_E_lower) or np.isnan(pN_E_upper)):
            allresults['P(N|E)_CI'] = [pN_E - pN_E_lower, pN_E_upper - pN_E]
    with open(junkPPIFile, 'wb') as fOut:
        pickle.dump(allresults, fOut)

if __name__ == "__main__":
    main()
