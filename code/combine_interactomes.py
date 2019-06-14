
import os
import pandas as pd
from pathlib import Path
from interactome_tools import remove_duplicate_PPIs, read_single_interface_annotated_interactome
from structural_annotation import filter_chain_annotations_by_protein

def main():
    
    # reference interactome name
    interactome_name = 'combined_interactome'
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # directory of processed data files specific to interactome
    interactomeDir = procDir / interactome_name
    
    # input data files
    chainMapFile = procDir / 'human_pdb_chain_map_filtered.txt'
    
    # output data files
    interfaceAnnotatedInteractomeFile = interactomeDir / 'human_interface_annotated_interactome.txt'
    strucInteractomeChainMapFile = interactomeDir / 'struc_interactome_pdb_chain_map.txt'
    
    if not interactomeDir.exists():
        os.makedirs(interactomeDir)
    
    Y2H_SI = pd.read_table(procDir / 'HI-II-14' / 'human_interface_annotated_interactome.txt')
    IntAct_SI = pd.read_table(procDir / 'IntAct' / 'human_interface_annotated_interactome.txt')
    
    Y2H_SI = Y2H_SI[["Protein_1", "Protein_2", "Interfaces", "Chain_pairs"]]
    IntAct_SI = IntAct_SI[["Protein_1", "Protein_2", "Interfaces", "Chain_pairs"]]
    
    combined = pd.concat([Y2H_SI, IntAct_SI], ignore_index=True)
    print('length with duplicates: %d' % len(combined))
    
    combined = remove_duplicate_PPIs (combined)
    print('length without duplicates: %d' % len(combined))
    
    combined.to_csv (interfaceAnnotatedInteractomeFile, index=False, sep='\t')
    
    structuralInteractome = read_single_interface_annotated_interactome (interfaceAnnotatedInteractomeFile)
    interactomeProteins = list(set(structuralInteractome[["Protein_1", "Protein_2"]].values.flatten()))
    print( '\n' + 'Structural interactome:' )
    print( '%d PPIs' % len(structuralInteractome) )
    print( '%d proteins' % len(interactomeProteins) )
    print()
    
    if not strucInteractomeChainMapFile.is_file():
        print('filtering chain annotations by structural interactome proteins')
        filter_chain_annotations_by_protein (chainMapFile,
                                             interactomeProteins,
                                             strucInteractomeChainMapFile)

if __name__ == "__main__":
    main()
