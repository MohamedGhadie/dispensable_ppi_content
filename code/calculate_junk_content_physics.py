#----------------------------------------------------------------------------------------
# Calculate the fraction of junk protein-protein interactions (PPIs) in the structural
# interactome, i.e., the fraction of PPIs that are effectively neutral upon perturbation.
# The fraction of junk PPIs is calculated from physics-based predictions of PPI 
# perturbations.
#
# Run the following script before running this script:
# - predict_perturbations_physics.py
#----------------------------------------------------------------------------------------

import os
import pickle
import numpy as np
from pathlib import Path
from text_tools import write_list_table
from stat_tools import fisher_test, sderror_on_fraction, proportion_ratio_CI
from plot_tools import pie_plot
from mutation_interface_edgotype import assign_edgotypes

def main():
    
    # reference interactome name
    # options: HI-II-14, IntAct
    interactome_name = 'IntAct'
    
    # method of calculating mutation ∆∆G for which results will be used
    # options: bindprofx, foldx
    ddg_method = 'bindprofx'
    
    # set to True to calculate junk PPI content using fraction of mono-edgetic mutations 
    # instead of fraction of edgetic mutations
    mono_edgetic = False
    
    # Minimum reduction in binding free energy DDG required for interaction perturbation
    ddgCutoff = 0.5
    
    # calculate confidence interval for the fraction of junk PPIs
    computeConfidenceIntervals = True
    
    # % confidence interval
    CI = 95
    
    # show figures
    showFigs = False
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # directory of processed data files specific to interactome
    interactomeDir = procDir / interactome_name
    
    # figure directory
    figDir = Path('../figures') / interactome_name
    
    # create output directories if not existing
    if not procDir.exists():
        os.makedirs(procDir)
    if not interactomeDir.exists():
        os.makedirs(interactomeDir)
    if not figDir.exists():
        os.makedirs(figDir)
        
    # input data files
    physicsPerturbsFile = interactomeDir / ('mutation_perturbs_physics_%s.pkl' % ddg_method)
    
    # output data files
    natMutEdgotypeFile = interactomeDir / ('nondisease_mutation_edgotype_physics_%s.txt' % ddg_method)
    disMutEdgotypeFile = interactomeDir / ('disease_mutation_edgotype_physics_%s.txt' % ddg_method)
    junkPPIFile = interactomeDir / ('fraction_junk_PPIs_physics_%s_%s.pkl' % (ddg_method, 'monoedgetic' if mono_edgetic else ''))
    
    #------------------------------------------------------------------------------------
    # Load interactome perturbations
    #------------------------------------------------------------------------------------
    
    with open(physicsPerturbsFile, 'rb') as f:
        naturalPerturbs, diseasePerturbs = pickle.load(f)
        
    #------------------------------------------------------------------------------------
    # Assign mutation edgotypes
    #------------------------------------------------------------------------------------
    
    print( '\n' + 'Labeling mutation edgotypes' )
    naturalPerturbs["edgotype"] = assign_edgotypes (naturalPerturbs["perturbations"].tolist(),
                                                    mono_edgetic = False)
    diseasePerturbs["edgotype"] = assign_edgotypes (diseasePerturbs["perturbations"].tolist(),
                                                    mono_edgetic = False)
    
    if mono_edgetic:
        print( '\n' + 'Labeling mono-edgetic mutations' )
        naturalPerturbs["mono-edgotype"] = assign_edgotypes (naturalPerturbs["perturbations"].tolist(),
                                                             mono_edgetic = True)
        diseasePerturbs["mono-edgotype"] = assign_edgotypes (diseasePerturbs["perturbations"].tolist(),
                                                             mono_edgetic = True)
    
    # write predicted mutation edgotypes to tab-delimited file
    write_list_table (naturalPerturbs, ["partners", "perturbations"], natMutEdgotypeFile)
    write_list_table (diseasePerturbs, ["partners", "perturbations"], disMutEdgotypeFile)
    
    #------------------------------------------------------------------------------------
    # Fraction of predicted edgetic mutations
    #------------------------------------------------------------------------------------
    
    if mono_edgetic:
        numNaturalMut_edgetic = sum(naturalPerturbs["mono-edgotype"] == 'mono-edgetic')
        numNaturalMut_nonedgetic = sum(naturalPerturbs["mono-edgotype"].apply(lambda x: 
                                                        x in ('non-edgetic', 'edgetic')))
        numDiseaseMut_edgetic = sum(diseasePerturbs["mono-edgotype"] == 'mono-edgetic')
        numDiseaseMut_nonedgetic = sum(diseasePerturbs["mono-edgotype"].apply(lambda x: 
                                                        x in ('non-edgetic', 'edgetic')))
    else:
        numNaturalMut_edgetic = sum(naturalPerturbs["edgotype"] == 'edgetic')
        numNaturalMut_nonedgetic = sum(naturalPerturbs["edgotype"] == 'non-edgetic')
        numDiseaseMut_edgetic = sum(diseasePerturbs["edgotype"] == 'edgetic')
        numDiseaseMut_nonedgetic = sum(diseasePerturbs["edgotype"] == 'non-edgetic')
    
    numNaturalMut_considered = numNaturalMut_edgetic + numNaturalMut_nonedgetic
    numDiseaseMut_considered = numDiseaseMut_edgetic + numDiseaseMut_nonedgetic
    
    fracNaturalMut_edgetic = numNaturalMut_edgetic / numNaturalMut_considered
    fracDiseaseMut_edgetic = numDiseaseMut_edgetic / numDiseaseMut_considered
    fracNaturalMut_error = sderror_on_fraction(numNaturalMut_edgetic, numNaturalMut_considered)
    fracDiseaseMut_error = sderror_on_fraction(numDiseaseMut_edgetic, numDiseaseMut_considered)
    
    label = 'monoedgetic' if mono_edgetic else 'edgetic'
    print( '\n' + 'Fraction of predicted %s mutations:' % label )
    print( 'Non-disease mutations: %.3f (SE = %g, %d out of %d)' % (fracNaturalMut_edgetic,
                                                                    fracNaturalMut_error,
                                                                    numNaturalMut_edgetic,
                                                                    numNaturalMut_considered) )
    
    print( 'Disease mutations: %.3f (SE = %g, %d out of %d)' % (fracDiseaseMut_edgetic,
                                                                fracDiseaseMut_error,
                                                                numDiseaseMut_edgetic,
                                                                numDiseaseMut_considered) )
    
    fisher_test([numNaturalMut_edgetic, numNaturalMut_nonedgetic],
                [numDiseaseMut_edgetic, numDiseaseMut_nonedgetic])
    
    pie_plot([numNaturalMut_nonedgetic, numNaturalMut_edgetic],
             angle = 90,
             colors = ['thistle', 'mediumslateblue'],
             show = showFigs,
             figdir = figDir,
             figname = 'non_disease_%s_mutations_physics' % label)
    pie_plot([numDiseaseMut_nonedgetic, numDiseaseMut_edgetic],
             angle=90,
             colors = ['thistle', 'mediumslateblue'],
             show = showFigs,
             figdir = figDir,
             figname = 'disease_%s_mutations_physics' % label)
    
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
