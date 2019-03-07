#----------------------------------------------------------------------------------------
# Predict interactome perturbations based on geometry. Given a structural interactome 
# with PPI binding interface annotations, common neutral mutations not associated 
# with disease as well as Mendelian disease-causing mutations are mapped onto the 
# structural interactome. Then, PPI perturbations caused by mutations are predicted based 
# on mutation location relative to interaction interface.
#
# Run the following scripts before running this script:
# - produce_structural_interactome.py
# - process_dbsnp_mutations.py
# - process_clinvar_mutations.py
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
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # directory of processed data files specific to interactome
    interactomeDir = procDir / interactome_name
    
    # create output directories if not existing
    if not procDir.exists():
        os.makedirs(procDir)
    if not interactomeDir.exists():
        os.makedirs(interactomeDir)
    
    # input data files
    naturalMutationsFile = procDir / 'dbsnp_mutations4.txt'
    diseaseMutationsFile = procDir / 'clinvar_mutations6.txt'
    interfaceAnnotatedInteractomeFile = interactomeDir / 'human_interface_annotated_interactome.txt'
    
    # output data files
    uniqueMutationPerturbsFile = interactomeDir / 'unique_mutation_perturbs_geometry.pkl'
    
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
    
    annotatedInteractome = read_interface_annotated_interactome (interfaceAnnotatedInteractomeFile)
    
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
