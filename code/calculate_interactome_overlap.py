#----------------------------------------------------------------------------------------
# Calculate overlap in PPIs and proteins between two interactomes.
#
# Run the following scripts for each interactome before running this script:
# - produce_data_mappings.py
# - process_interactome.py
#
# If comparing structural interactomes, also run the following script:
# - produce_structural_interactome.py
#----------------------------------------------------------------------------------------

import pandas as pd
from pathlib import Path

def main():
    
    # names of reference interactomes
    interactome_names = ['HI-II-14', 'IntAct']
    
    # type of interactomes to compare
    # options: reference, structural
    interactome_type = 'structural'
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # input data files
    if interactome_type is 'reference':
        interactome_File1 = procDir / interactome_names[0] / 'human_interactome.txt'
        interactome_File2 = procDir / interactome_names[1] / 'human_interactome.txt'    
    elif interactome_type is 'structural':
        interactome_File1 = procDir / interactome_names[0] / 'human_interface_annotated_interactome.txt'
        interactome_File2 = procDir / interactome_names[1] / 'human_interface_annotated_interactome.txt'
    
    interactome1 = pd.read_table (interactome_File1, sep='\t')
    interactome2 = pd.read_table (interactome_File2, sep='\t')
    
    ppis_1 = {tuple(sorted(ppi)) for ppi in interactome1[["Protein_1","Protein_2"]].values}
    ppis_2 = {tuple(sorted(ppi)) for ppi in interactome2[["Protein_1","Protein_2"]].values}
    shared = len(ppis_1 & ppis_2)
    
    print('\n' + 'Shared PPIs:')
    print('%s: %d out of %d PPIs (%f %%)' % (interactome_names[0],
                                             shared,
                                             len(ppis_1),
                                             100 * shared / len(ppis_1)))
    print('%s: %d out of %d PPIs (%f %%)' % (interactome_names[1],
                                             shared,
                                             len(ppis_2),
                                             100 * shared / len(ppis_2)))
    
    proteins_1 = set(interactome1[["Protein_1","Protein_2"]].values.flatten())
    proteins_2 = set(interactome2[["Protein_1","Protein_2"]].values.flatten())
    shared = len(proteins_1 & proteins_2)
    
    print('\n' + 'Shared proteins')
    print('%s = %d out of %d (%f %%)' % (interactome_names[0],
                                         shared,
                                         len(proteins_1),
                                         100 * shared / len(proteins_1)))
    print('%s = %d out of %d (%f %%)' % (interactome_names[1],
                                         shared,
                                         len(proteins_2),
                                         100 * shared / len(proteins_2)))
    
if __name__ == '__main__':
    main()
