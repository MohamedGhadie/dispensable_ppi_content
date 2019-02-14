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
from stat_tools import fisher_test, sderror_on_fraction, proportion_ratio_CI
from plot_tools import pie_plot

def main():
    
    run_locally = True
    
    interactome_name = 'experiment'
    
    # calculate confidence interval for the fraction of junk PPIs
    computeConfidenceIntervals = True
    
    # treat experiment quasi-null mutation as edgetic
    consider_exp_QuasiNull_perturbs = False
    
    # % confidence interval
    CI = 95
    
    # show figures
    showFigs = False
    
    dataDir = Path('/Volumes/MG_Samsung/junk_ppi_content/data')
    
    # directory for data files from external sources
    inDir = dataDir / 'external'
    
    # directory to save processed data shared by all interactomes
    outDir = dataDir / 'processed' / interactome_name
    
    # figure directory
    if consider_exp_QuasiNull_perturbs:
        figDir = Path('../figures') / interactome_name / 'exp_QuasiNull_perturbs_considered'
    else:
        figDir = Path('../figures') / interactome_name
    
    # directory for PDB structure files
    if run_locally:
        # directory for PDB structure files on local computer
        pdbDir = Path('/Volumes/MG_Samsung/pdb_files')
    else:
        # directory for PDB structure files on server
        pdbDir = Path('../../pdb_files')
    
    # create output directories if not existing
    if not outDir.exists():
        os.makedirs(outDir)
    if not figDir.exists():
        os.makedirs(figDir)
    
    # input data files
    edgotypeFile = inDir / "Sahni_2015_Table_S3.xlsx"
    
    # output data files
    junkPPIFile = outDir / 'fraction_junk_PPIs_experiment.pkl'
                   
    #------------------------------------------------------------------------------------
    # process experimental edgotype data
    #------------------------------------------------------------------------------------
    
    expMutations = pd.read_excel( edgotypeFile,
                                  sheet_name = 'Table S3C' )
    expMutations_ppi = pd.read_excel( inDir / "Sahni_2015_Table_S3.xlsx",
                                      sheet_name = 'Table S3A' )
    
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
             #labels = ['Nonedgetic', 'Edgetic'],
             colors = ['thistle', 'mediumslateblue'],
             show = showFigs,
             figdir = figDir,
             figname = 'non_disease_edgetic_mutations_exp')
    pie_plot([numDiseaseMut_nonedgetic, numDiseaseMut_edgetic],
             angle = 90,
             #labels = ['Nonedgetic', 'Edgetic'],
             colors = ['thistle', 'mediumslateblue'],
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
