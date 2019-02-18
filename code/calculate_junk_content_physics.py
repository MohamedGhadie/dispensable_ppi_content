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
from ddg_tools import read_protein_mutation_ddg
from stat_tools import t_test, fisher_test, sderror, sderror_on_fraction, proportion_ratio_CI
from plot_tools import bar_plot, multi_histogram_plot
from mutation_interface_edgotype import energy_based_perturbation

def main():
    
    # valid reference interactome names
    interactome_names = ['HI-II-14', 'IntAct']
    
    # choose interactome (index in interactome_names)
    interactome_choise = 0
    
    # Minimum reduction in binding free energy DDG required for interaction perturbation
    ddgCutoff = 0.5
    
    # calculate confidence interval for the fraction of junk PPIs
    computeConfidenceIntervals = True
    
    # % confidence interval
    CI = 95
    
    # show figures
    showFigs = False
    
    # selected reference interactome
    interactome_name = interactome_names [interactome_choise]
    
    dataDir = Path('/Volumes/MG_Samsung/junk_ppi_content/data')
    
    # directory to save processed data shared by all interactomes
    inDir = dataDir / 'processed' / interactome_name
    
    # directory to save processed data specific to interactome
    outDir = dataDir / 'processed' / interactome_name
    
    # figure directory
    figDir = Path('../figures') / interactome_name
    
    # create output directories if not existing
    if not outDir.exists():
        os.makedirs(outDir)
    if not figDir.exists():
        os.makedirs(figDir)
    
    # input data files
    naturalMutationsDDGFile = inDir / 'nondisease_mutations_onchains_ddg.txt'
    diseaseMutationsDDGFile = inDir / 'disease_mutations_onchains_ddg.txt'
    geometryPerturbFile = inDir / 'filtered_mutation_perturbs_geometry_hgmd.pkl'
    
    # output data files
    physicsPerturbFile = outDir / 'mutation_perturbs_physics.pkl'
    ddgOutputFile = outDir / 'mutation_∆∆G.xlsx'
    edgotypeFile = outDir / 'mutation_edgotypes_physics.xlsx'
    junkPPIFile = outDir / 'fraction_junk_PPIs_physics.pkl'
        
    #------------------------------------------------------------------------------------
    # Fraction of mutation-targeted PPIs with ∆∆G exceeding a specified cutoff
    #------------------------------------------------------------------------------------
    
    # read change in binding free energy for interfacial mutations mapped on PDB chains
    naturalMutationsDDG = read_protein_mutation_ddg(naturalMutationsDDGFile, 'binding')
    diseaseMutationsDDG = read_protein_mutation_ddg(diseaseMutationsDDGFile, 'binding')
    
    naturalMutations = pd.DataFrame(columns=["protein", "partner", "protein_pos", "pdb_id", 
                                             "chain_mut", "partner_chain", "ddg"])
    for i, item in enumerate(naturalMutationsDDG.items()):
        naturalMutations.loc[i] = item[0] + item[1]
    
    diseaseMutations = pd.DataFrame(columns=["protein", "partner", "protein_pos", "pdb_id", 
                                             "chain_mut", "partner_chain", "ddg"])
    for i, item in enumerate(diseaseMutationsDDG.items()):
        diseaseMutations.loc[i] = item[0] + item[1]
    
    writer = pd.ExcelWriter(str(ddgOutputFile))
    naturalMutations.to_excel(writer, sheet_name = 'nondisease_mutations_∆∆G')
    diseaseMutations.to_excel(writer, sheet_name = 'disease_mutations_∆∆G')
    writer.save()
    
    natMutDDGs = naturalMutations["ddg"].apply(float).values
    disMutDDGs = diseaseMutations["ddg"].apply(float).values
    
    numNatural_ddg_considered = len(natMutDDGs)
    numDisease_ddg_considered = len(disMutDDGs)
    
    print( '\n' + 'Avg. change in binding energy (∆∆G) for mutation-targeted PPIs:' )
    print( 'Non-disease: %.1f (SE = %g, n = %d)' % (np.mean(natMutDDGs),
                                                    sderror(natMutDDGs),
                                                    numNatural_ddg_considered) )
    
    print( 'Disease: %.1f (SE = %g, n = %d)' % (np.mean(disMutDDGs),
                                                sderror(disMutDDGs),
                                                numDisease_ddg_considered) )
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
                          figname = 'mut_ddg_histogram')
    
    numNatural_ddg = sum( natMutDDGs > ddgCutoff )
    numDisease_ddg = sum( disMutDDGs > ddgCutoff )
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
             figname = 'Sample_fraction_ddg_>%.1f' % ddgCutoff)
    
    #------------------------------------------------------------------------------------
    # Fractions of predicted edgetic and mono-edgetic mutations using physics-based 
    # prediction of PPI perturbation (i.e., using binding energy change)
    #------------------------------------------------------------------------------------
    
    if physicsPerturbFile.is_file():
        print( '\n' + 'Loading physics-based mutation edgotype predictions' )
        with open(physicsPerturbFile, 'rb') as f:
            naturalPerturbs, diseasePerturbs = pickle.load(f)
    else:
        if geometryPerturbFile.is_file():
            print( '\n' + 'Loading geometry-based PPI perturbation predictions' )
            with open(geometryPerturbFile, 'rb') as f:
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
        with open(physicsPerturbFile, 'wb') as fOut:
            pickle.dump([naturalPerturbs, diseasePerturbs], fOut)
        
    perturbations = [[p for p in pert if not np.isnan(p)] for pert in naturalPerturbs["perturbations"].values]
    perturbations = [p for p in perturbations if p]
    edgetic = [sum(np.array(p) > 0) > 0 for p in perturbations]
    numNaturalMut_edgetic = sum(edgetic)
    numNaturalMut_considered = len(edgetic)
    fracNaturalMut_edgetic = numNaturalMut_edgetic / numNaturalMut_considered
    
    print( 'Fraction of predicted non-disease edgetic mutations: %f (SEF = %g, ntot = %d)' 
            % (fracNaturalMut_edgetic,
               sderror_on_fraction(numNaturalMut_edgetic, numNaturalMut_considered),
               numNaturalMut_considered) )
    
    perturbations = [[p for p in pert if not np.isnan(p)] for pert in diseasePerturbs["perturbations"].values]
    perturbations = [p for p in perturbations if p]
    edgetic = [sum(np.array(p) > 0) > 0 for p in perturbations]
    numDiseaseMut_edgetic = sum(edgetic)
    numDiseaseMut_considered = len(edgetic)
    fracDiseaseMut_edgetic = numDiseaseMut_edgetic / numDiseaseMut_considered
    
    print( 'Fraction of predicted disease edgetic mutations: %f (SEF = %g, ntot = %d)' 
            % (fracDiseaseMut_edgetic,
               sderror_on_fraction(numDiseaseMut_edgetic, numDiseaseMut_considered),
               numDiseaseMut_considered) )
    
    naturalPerturbsWrite = naturalPerturbs.copy()
    naturalPerturbsWrite["partners"] = naturalPerturbsWrite["partners"].apply(lambda x: ','.join(x))
    naturalPerturbsWrite["perturbations"] = naturalPerturbsWrite["perturbations"].apply(lambda x: ','.join(map(str, x)))
    
    diseasePerturbsWrite = diseasePerturbs.copy()
    diseasePerturbsWrite["partners"] = diseasePerturbsWrite["partners"].apply(lambda x: ','.join(x))
    diseasePerturbsWrite["perturbations"] = diseasePerturbsWrite["perturbations"].apply(lambda x: ','.join(map(str, x)))
    
    writer = pd.ExcelWriter(str(edgotypeFile))
    naturalPerturbsWrite.to_excel(writer, sheet_name = 'nondisease_mutations')
    diseasePerturbsWrite.to_excel(writer, sheet_name = 'disease_mutations')
    writer.save()
        
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
