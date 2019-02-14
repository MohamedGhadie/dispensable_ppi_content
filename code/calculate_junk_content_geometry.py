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
from interactome_tools import read_interface_annotated_interactome
from stat_tools import fisher_test, sderror_on_fraction, proportion_ratio_CI
from plot_tools import pie_plot
from mutation_interface_edgotype import mutation_PPI_interface_perturbations

def main():
    
    # valid reference interactome names
    interactome_names = ['HI-II-14', 'IntAct']
    
    # choose interactome (index in interactome_names)
    interactome_choise = 1
    
    # database options for processed disease mutations
    disMut_dbs = ['HGMD', 'ClinVar']
    
    # disease mutation database choice (index in disease_databases)
    disMut_db_choice = 1
    
    # selected disease mutation database
    disMut_db = disMut_dbs [disMut_db_choice]
    
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
    inDir = dataDir / 'processed'
    
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
    ProteinSeqFile = inDir / 'human_reference_sequences.pkl'    
    chainSeqFile = inDir / 'chain_sequences.pkl'
    chainListFile = inDir / 'pdb_seqres_chains.list'
    pdbChainsFile = inDir / 'pdb_seqres_chains.pkl'
    chainResOrderFile = inDir / 'chain_seqres_order.pkl'
    chainMapFile = inDir / 'human_pdb_chain_map_filtered.txt'
    dbsnpMutationsFile = inDir / 'dbsnp_mutations4.txt'
    hgmdMutationsFile = inDir / 'HGMD_missense_DM_mutations_matched.txt'
    clinvarMutationsFile = inDir / 'clinvar_mutations6.txt'
    interfaceAnnotatedInteractomeFile = outDir / 'human_interface_annotated_interactome.txt'
    
    # output data files
    mutationPerturbFile = outDir / 'mutation_perturbs_geometry.pkl'
    filteredMutationPerturbFile = outDir / 'filtered_mutation_perturbs_geometry.pkl'
    natMutEdgotypeFile = outDir / 'nondisease_mutation_edgotype_geometry.txt'
    disMutEdgotypeFile = outDir / 'disease_mutation_edgotype_geometry.txt'
    edgotypeFile = outDir / 'mutation_edgotypes_geometry.xlsx'
    junkPPIFile = outDir / 'fraction_junk_PPIs_geometry.pkl'
    
    #------------------------------------------------------------------------------------
    # further process mutations
    #------------------------------------------------------------------------------------
    
    # load processed mutations
    print( '\n' + 'Reading processed mutations' )
    naturalMutations = pd.read_table(dbsnpMutationsFile, sep='\t')
    if disMut_db is 'HGMD':
        diseaseMutations = pd.read_table(hgmdMutationsFile, sep='\t')
    elif disMut_db is 'ClinVar':
        diseaseMutations = pd.read_table(clinvarMutationsFile, sep='\t')
    
    naturalSequenceMatch = sum(naturalMutations["seq_match"])
    diseaseSequenceMatch = sum(diseaseMutations["Seq_match"])
    print( '\n' + '%d non-disease mutations matching to UniProt sequence (%.1f %%)'
            % (naturalSequenceMatch, 100 * naturalSequenceMatch / len(naturalMutations)) )
    print( '%d disease mutations matching to UniProt sequence (%.1f %%)'
            % (diseaseSequenceMatch, 100 * diseaseSequenceMatch / len(diseaseMutations)) )
    
    # keep only mutations matching to UniProt sequence
    naturalMutations = naturalMutations[naturalMutations["seq_match"]].reset_index(drop=True)
    diseaseMutations = diseaseMutations[diseaseMutations["Seq_match"]].reset_index(drop=True)
    
    # rename columns
    naturalMutations.rename(columns={"protein":"Protein",
                                     "residue":"Mut_res",
                                     "aa_position":"Mutation_Position",
                                     "context":"WT_context",
                                     "mut_context_pos":"Mut_context_pos",
                                     "allele":"Frequency_reporting_allele",
                                     "allele.1": "Variation_allele"}, inplace=True)
    
    naturalMutations["wt_res"] = naturalMutations.apply(lambda x:
                                                        x["WT_context"][x["Mut_context_pos"]-1],
                                                        axis=1)
    
    if disMut_db is 'HGMD':
        diseaseMutations.rename(columns={"Codon_number":"Mutation_Position"}, inplace=True)
        diseaseMutations["Mut_res"] = diseaseMutations.apply(lambda x:
                                                             x["Mut_context"][x["Mut_context_pos"] - 1],
                                                             axis=1)
        diseaseMutations["wt_res"] = diseaseMutations.apply(lambda x: 
                                                            x["WT_context"][x["Mut_context_pos"]-1],
                                                            axis=1)
    elif disMut_db is 'ClinVar':
        diseaseMutations.rename(columns={"mut_res":"Mut_res",
                                         "aa_position":"Mutation_Position",
                                         "mut_context_pos":"Mut_context_pos"}, inplace=True)
    
    # identify common mutations overlapping in position with disease mutations
    df = diseaseMutations[ ["Protein", "Mutation_Position"] ].copy()
    df2 = naturalMutations[ ["Protein", "Mutation_Position"] ].copy()
    df2 = df2.drop_duplicates().reset_index(drop=True)
    df = df.append(df2, ignore_index=True)
    duplicates = df.duplicated(keep='first')[ len(diseaseMutations) : ]
    duplicates = duplicates.reset_index(drop=True)
    print( '\n' + '%d unique non-disease variants overlap in location with disease mutations' % sum(duplicates) )
    
    # remove common mutations overlapping in position with disease mutations
    for i, row in df2[ duplicates ].iterrows():
        todrop = ( ( naturalMutations["Protein"] == row.Protein ) & 
                   ( naturalMutations["Mutation_Position"] == row.Mutation_Position ) )
        naturalMutations = naturalMutations[todrop == False]
    print( 'Overlaping non-disease mutations removed' )
    
    numNaturalMutations = len(naturalMutations)
    numDiseaseMutations = len(diseaseMutations)
    
    # remove invalid mutations, i.e., those with mutation residue similar to wild type
    naturalMutations = naturalMutations [naturalMutations["wt_res"] != naturalMutations["Mut_res"]]
    diseaseMutations = diseaseMutations [diseaseMutations["wt_res"] != diseaseMutations["Mut_res"]]
    
    print( '\n' + 'Number of invalid mutations removed (WT = Mut)' )
    print( 'non-disease invalid mutations: %d' % (numNaturalMutations - len(naturalMutations)) )
    print( 'disease invalid mutations: %d' % (numDiseaseMutations - len(diseaseMutations)) )
    
    numNaturalMutations = len(naturalMutations)
    numDiseaseMutations = len(diseaseMutations)
    
    print( '\n' + 'Number of mutations after removing invalid mutations:' )
    print( 'non-disease mutations: %d' % numNaturalMutations )
    print( 'disease mutations: %d' % numDiseaseMutations )
    
    # remove duplicate mutations, by position and mutant residue
    naturalMutations = naturalMutations.drop_duplicates(subset=["Protein",
                                                                "Mutation_Position",
                                                                "Mut_res"]).reset_index(drop=True)
    diseaseMutations = diseaseMutations.drop_duplicates(subset=["Protein",
                                                                "Mutation_Position",
                                                                "Mut_res"]).reset_index(drop=True)
    
    numNaturalMutations = len(naturalMutations)
    numDiseaseMutations = len(diseaseMutations)
    
    print( '\n' + 'Number of mutations after removing duplicates by position and mutant residue:' )
    print( 'non-disease mutations: %d' % numNaturalMutations )
    print( 'disease mutations: %d' % numDiseaseMutations )
    
    # plot fractions of disease and non-disease mutations
    pie_plot([numNaturalMutations, numDiseaseMutations],
             labels=['Non-disease\nmutations', 'Disease\nmutations'],
             angle=90,
             colors = ['turquoise', 'orangered'],
             show = showFigs,
             figdir = figDir,
             figname = 'mutations')
    
    writer = pd.ExcelWriter(str(outDir / 'all_disease_mutations.xlsx'))
    diseaseMutations.to_excel (writer, sheet_name = 'disease_mutations')
    writer.save()
    diseaseMutations.to_csv(outDir / 'all_disease_mutations.txt', index=False, sep='\t')
    
    #------------------------------------------------------------------------------------
    # predict PPI perturbations
    #------------------------------------------------------------------------------------
    
    # Consider PPI perturbations only for PPIs with this maximum number of interfaces 
    # mapped from distinct PDB binding chain pairs.
    # This parameter is irrelevent if the flag "merge_interfaces" is True.
    # Set to inf for unlimited number of interfaces.
    maxInterfaces = np.inf
    
    # Predict PPI perturbation if mutation is this number of residues away in sequence 
    # from an interface residue. Set to 0 if mutation must be exactly at interface residue.
    numResFromInterface = 0
    
    # Consider PPI perturbations only for PPIs with this minimum number of partners
    minPartners = 1
    
    if mutationPerturbFile.is_file():
        print( '\n' + 'Loading geometry-based PPI perturbation predictions' )
        with open(mutationPerturbFile, 'rb') as f:
            naturalMutations_perturb, diseaseMutations_perturb = pickle.load(f)
    else:
        print( '\n' + 'Predicting PPI perturbations based on geometry' )
        annotatedInteractome = read_interface_annotated_interactome(interfaceAnnotatedInteractomeFile)
        naturalMutations_perturb = mutation_PPI_interface_perturbations(naturalMutations,
                                                                        annotatedInteractome,
                                                                        maxInterfaces,
                                                                        numResFromInterface)
        diseaseMutations_perturb = mutation_PPI_interface_perturbations(diseaseMutations,
                                                                        annotatedInteractome,
                                                                        maxInterfaces,
                                                                        numResFromInterface)
        with open(mutationPerturbFile, 'wb') as fOut:
            pickle.dump([naturalMutations_perturb, diseaseMutations_perturb], fOut)
    
    naturalPerturbs = naturalMutations.copy()
    naturalPerturbs["partners"] = naturalMutations_perturb.apply(lambda x: x[0])
    naturalPerturbs["perturbations"] = naturalMutations_perturb.apply(lambda x: x[1])
    naturalPerturbs = naturalPerturbs[naturalPerturbs["partners"].apply(len) 
                                      >= minPartners].reset_index(drop=True)
        
    diseasePerturbs = diseaseMutations.copy()
    diseasePerturbs["partners"] = diseaseMutations_perturb.apply(lambda x: x[0])
    diseasePerturbs["perturbations"] = diseaseMutations_perturb.apply(lambda x: x[1])
    diseasePerturbs = diseasePerturbs[diseasePerturbs["partners"].apply(len) 
                                      >= minPartners].reset_index(drop=True)
        
    # drop duplicate mutations based on location, regardless of residue type 
    naturalPerturbs = naturalPerturbs.drop_duplicates(subset=["Protein",
                                                              "Mutation_Position"]).reset_index(drop=True)
    diseasePerturbs = diseasePerturbs.drop_duplicates(subset=["Protein",
                                                              "Mutation_Position"]).reset_index(drop=True)
    
    print( '\n' + 'Number of mutations with PPI perturbation predictions after removing duplicates by position' )
    print( 'non-disease: %d' % len(naturalPerturbs) )
    print( 'disease: %d' % len(diseasePerturbs) )
    
    with open(filteredMutationPerturbFile, 'wb') as fOut:
        pickle.dump([naturalPerturbs, diseasePerturbs], fOut)
    
    #------------------------------------------------------------------------------------
    # Assign mutation edgotypes
    #------------------------------------------------------------------------------------
    
    print( '\n' + 'Labeling mutation edgotypes' )
    naturalPerturbs["Edgotype"] = naturalPerturbs["perturbations"].apply(lambda x:
                                                                         'Edgetic' if sum(x > 0) > 0
                                                                         else 'Non-edgetic')
    diseasePerturbs["Edgotype"] = diseasePerturbs["perturbations"].apply(lambda x:
                                                                         'Edgetic' if sum(x > 0) > 0
                                                                         else 'Non-edgetic')
    
    # write predicted interactome perturbations by non-disease mutations to excel file
    naturalPerturbs_output = naturalPerturbs.copy()
    naturalPerturbs_output["partners"] = naturalPerturbs_output["partners"].apply(lambda x: ','.join(x))
    naturalPerturbs_output["perturbations"] = ( naturalPerturbs_output["perturbations"]
                                                .apply(lambda x: ','.join(map(str, [int(np.ceil(k)) for k in x]))) )
    
    # write predicted interactome perturbations by disease mutations to excel file
    diseasePerturbs_output = diseasePerturbs.copy()
    diseasePerturbs_output["partners"] = diseasePerturbs_output["partners"].apply(lambda x: ','.join(x))
    diseasePerturbs_output["perturbations"] = ( diseasePerturbs_output["perturbations"]
                                                .apply(lambda x: ','.join(map(str, [int(np.ceil(k)) for k in x]))) )
    
    naturalPerturbs_output.to_csv(natMutEdgotypeFile, index=False, sep='\t')
    diseasePerturbs_output.to_csv(disMutEdgotypeFile, index=False, sep='\t')
    
    writer = pd.ExcelWriter(str(edgotypeFile))
    naturalPerturbs_output.to_excel (writer, sheet_name = 'Nondisease_mutations')
    #output_columns = diseasePerturbs_output.columns.tolist()
    diseasePerturbs_output.to_excel (writer, sheet_name = 'Disease_mutations')
    writer.save()
        
    #------------------------------------------------------------------------------------
    # Fraction of predicted edgetic mutations
    #------------------------------------------------------------------------------------
        
    numNaturalMut_edgetic = sum(naturalPerturbs["Edgotype"] == 'Edgetic')
    numNaturalMut_nonedgetic = sum(naturalPerturbs["Edgotype"] == 'Non-edgetic')
    numDiseaseMut_edgetic = sum(diseasePerturbs["Edgotype"] == 'Edgetic')
    numDiseaseMut_nonedgetic = sum(diseasePerturbs["Edgotype"] == 'Non-edgetic')
    
    numNaturalMut_considered = numNaturalMut_edgetic + numNaturalMut_nonedgetic
    numDiseaseMut_considered = numDiseaseMut_edgetic + numDiseaseMut_nonedgetic
    
    fracNaturalMut_edgetic = numNaturalMut_edgetic / numNaturalMut_considered
    fracDiseaseMut_edgetic = numDiseaseMut_edgetic / numDiseaseMut_considered
    fracNaturalMut_error = sderror_on_fraction(numNaturalMut_edgetic, numNaturalMut_considered)
    fracDiseaseMut_error = sderror_on_fraction(numDiseaseMut_edgetic, numDiseaseMut_considered)
    
    print( '\n' + 'Fraction of predicted edgetic mutations:' )
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
             figname = 'non_disease_edgetic_mutations_geometrybased')
    pie_plot([numDiseaseMut_nonedgetic, numDiseaseMut_edgetic],
             angle=90,
             colors = ['thistle', 'mediumslateblue'],
             show = showFigs,
             figdir = figDir,
             figname = 'disease_edgetic_mutations_geometrybased')
            
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
