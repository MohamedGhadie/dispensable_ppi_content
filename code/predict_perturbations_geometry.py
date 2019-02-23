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
import numpy as np
from pathlib import Path
from interactome_tools import read_interface_annotated_interactome
from mutation_processing_tools import remove_mutation_overlaps
from mutation_interface_edgotype import mutation_PPI_interface_perturbations

def main():
    
    # reference interactome name. Options: 'HI-II-14' or 'IntAct'
    interactome_name = 'IntAct'
    
    # directory of processed data files shared by all interactomes
    dataDir = Path('../data') / 'processed'
    
    # directory of processed data files specific to interactome
    interactomeDir = dataDir / interactome_name
    
    # directory to save output data files
    outDir = dataDir / interactome_name
    
    # create output directories if not existing
    if not interactomeDir.exists():
        os.makedirs(interactomeDir)
    if not outDir.exists():
        os.makedirs(outDir)
    
    # input data files
    naturalMutationsFile = dataDir / 'dbsnp_mutations4.txt'
    diseaseMutationsFile = dataDir / 'clinvar_mutations6.txt'
    interfaceAnnotatedInteractomeFile = interactomeDir / 'human_interface_annotated_interactome.txt'
    
    # output data files
    uniqueMutationPerturbsFile = outDir / 'unique_mutation_perturbs_geometry.pkl'
    
    #------------------------------------------------------------------------------------
    # further process mutations
    #------------------------------------------------------------------------------------
    
    naturalMutations, diseaseMutations = remove_mutation_overlaps (naturalMutationsFile,
                                                                   diseaseMutationsFile)
    
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
    
    annotatedInteractome = read_interface_annotated_interactome(interfaceAnnotatedInteractomeFile)
    
    print( '\n' + 'Predicting PPI perturbations based on geometry' )
    naturalMutation_perturbs = mutation_PPI_interface_perturbations(naturalMutations,
                                                                    annotatedInteractome,
                                                                    maxInterfaces = maxInterfaces,
                                                                    dist = numResFromInterface)
    diseaseMutation_perturbs = mutation_PPI_interface_perturbations(diseaseMutations,
                                                                    annotatedInteractome,
                                                                    maxInterfaces,
                                                                    dist = numResFromInterface)
    
    naturalPerturbs = naturalMutations.copy()
    naturalPerturbs["partners"], naturalPerturbs["perturbations"] = zip(* naturalMutation_perturbs)
    naturalPerturbs = naturalPerturbs[naturalPerturbs["partners"].apply(len) >= minPartners]
        
    diseasePerturbs = diseaseMutations.copy()
    diseasePerturbs["partners"], diseasePerturbs["perturbations"] = zip(* diseaseMutation_perturbs)
    diseasePerturbs = diseasePerturbs[diseasePerturbs["partners"].apply(len) >= minPartners]
        
    # drop duplicate mutations based on location, regardless of residue type 
    naturalPerturbs = naturalPerturbs.drop_duplicates(subset=["protein", "mut_position"]).reset_index(drop=True)
    diseasePerturbs = diseasePerturbs.drop_duplicates(subset=["protein", "mut_position"]).reset_index(drop=True)
    
    print( '\n' + 'Number of mutations with PPI perturbation predictions after removing duplicates by position' )
    print( 'non-disease: %d' % len(naturalPerturbs) )
    print( 'disease: %d' % len(diseasePerturbs) )
    
    with open(uniqueMutationPerturbsFile, 'wb') as fOut:
        pickle.dump([naturalPerturbs, diseasePerturbs], fOut)

if __name__ == "__main__":
    main()
