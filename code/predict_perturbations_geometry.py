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
from mutation_interface_edgotype import mutation_PPI_interface_perturbations, create_perturbed_network
from plot_tools import network_plot

def main():
    
    # reference interactome name
    # options: HI-II-14, IntAct
    interactome_name = 'HI-II-14'
    
    # plot perturbed interactome and produce files for use by Cytoscape
    plot_perturbations = False
    
    # show figures
    showFigs = False
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # directory of processed data files specific to interactome
    interactomeDir = procDir / interactome_name
    
    # directory of network perturbation output data files for use by Cytoscape
    cytoscapeDir = interactomeDir / 'cytoscape'
        
    # figure directory
    figDir = Path('../figures') / interactome_name
    
    # input data files
    naturalMutationsFile = procDir / 'dbsnp_mutations4.txt'
    diseaseMutationsFile = procDir / 'clinvar_mutations6.txt'
    interfaceAnnotatedInteractomeFile = interactomeDir / 'human_interface_annotated_interactome.txt'
    
    # output data files
    uniqueMutationPerturbsFile = interactomeDir / 'unique_mutation_perturbs_geometry.pkl'
    naturalMutEdgeFile = cytoscapeDir / 'nondiseaseMut_perturbed_edges_geometry'
    naturalMutNodeFile = cytoscapeDir / 'nondiseaseMut_node_colors_geometry'
    diseaseMutEdgeFile = cytoscapeDir / 'diseaseMut_perturbed_edges_geometry'
    diseaseMutNodeFile = cytoscapeDir / 'diseaseMut_node_colors_geometry'
    
    # create output directories if not existing
    if not interactomeDir.exists():
        os.makedirs(interactomeDir)
    if not cytoscapeDir.exists():
        os.makedirs(cytoscapeDir)
    if not figDir.exists():
        os.makedirs(figDir)
    
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
    
    structuralInteractome = read_interface_annotated_interactome (interfaceAnnotatedInteractomeFile)
    
    print( '\n' + 'Predicting PPI perturbations based on geometry' )
    naturalMutation_perturbs = mutation_PPI_interface_perturbations(naturalMutations,
                                                                    structuralInteractome,
                                                                    maxInterfaces = maxInterfaces,
                                                                    dist = numResFromInterface)
    diseaseMutation_perturbs = mutation_PPI_interface_perturbations(diseaseMutations,
                                                                    structuralInteractome,
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
    
    #------------------------------------------------------------------------------------
    # plot network perturbations
    #------------------------------------------------------------------------------------
    
    if plot_perturbations:
        print( '\n' + 'Creating network perturbed by non-disease mutations' )
        nodes, edges, nodeColors, edgeColors = create_perturbed_network (structuralInteractome,
                                                                         naturalPerturbs,
                                                                         naturalMutEdgeFile,
                                                                         naturalMutNodeFile)
        network_plot (edges,
                      nodes = nodes,
                      nodeSizes = [20] * len(nodes),
                      edgeWidth = 1,
                      nodeColors = nodeColors,
                      edgeColors = edgeColors,
                      show = showFigs,
                      figdir = figDir,
                      figname = 'nondisease_mut_perturbed_interactome_geometry')
    
        print( '\n' + 'Creating network perturbed by disease mutations' )
        nodes, edges, nodeColors, edgeColors = create_perturbed_network (structuralInteractome,
                                                                         diseasePerturbs,
                                                                         diseaseMutEdgeFile,
                                                                         diseaseMutNodeFile)
        network_plot (edges,
                      nodes = nodes,
                      nodeSizes = [20] * len(nodes),
                      edgeWidth = 1,
                      nodeColors = nodeColors,
                      edgeColors = edgeColors,
                      show = showFigs,
                      figdir = figDir,
                      figname = 'disease_mut_perturbed_interactome_geometry')

if __name__ == "__main__":
    main()
