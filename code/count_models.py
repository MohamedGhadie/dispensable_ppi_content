#----------------------------------------------------------------------------------------
# Count the number of PPIs with co-crystal structures.
#
# Run the following scripts before running this script:
# - 
# - 
# - 
#----------------------------------------------------------------------------------------

import pandas as pd
from pathlib import Path
from interactome_tools import read_single_interface_annotated_interactome
from structural_annotation import is_cocrystal

def main():
    
    # reference interactome name
    # options: HI-II-14, IntAct
    interactome_name = 'IntAct'
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # directory of processed data files specific to interactome
    interactomeDir = procDir / interactome_name
    
    # input data files
    structuralInteractomeFile = interactomeDir / 'human_interface_annotated_interactome.txt'
    chainMapFile = interactomeDir / 'struc_interactome_pdb_chain_map.txt'
    
    #------------------------------------------------------------------------------------
    # load reference and structural interactomes
    #------------------------------------------------------------------------------------
    
    interactome = read_single_interface_annotated_interactome (structuralInteractomeFile)
    chainMap = pd.read_table (chainMapFile, sep='\t')
    
    crystals = 0
    for _, row in interactome.iterrows():
        for chain_1, chain_2 in row.Chain_pairs:
            if is_cocrystal (row.Protein_1, chain_1, chainMap):
                if is_cocrystal (row.Protein_2, chain_2, chainMap):
                    crystals += 1
                    break
    
    models = len(interactome) - crystals
    print('Number of PPIs with co-crystal structures: %d' % crystals)
    print('Number of PPIs with homology models: %d' % models)

if __name__ == '__main__':
    main()
