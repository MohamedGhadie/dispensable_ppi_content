from pathlib import Path
from text_tools import read_list_table
from interactome_tools import read_single_interface_annotated_interactome

def main():
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # directory of processed data files specific to interactome
    interactomeDir = procDir / 'HI-II-14'
    
    oldInterfaceFile = procDir / 'pdb_interfaces_reserved.txt'
    newInterfaceFile = procDir / 'pdb_interfaces.txt'
    
    old_interfaces = {}
    interface_df = read_list_table (oldInterfaceFile, "Chain1_interface", int, '\t')
    for _, row in interface_df.iterrows():
        chainKey = row.Chain_1 + '-' + row.Chain_2
        if -1 in row.Chain1_interface:
            old_interfaces[chainKey] = []
        else:
            old_interfaces[chainKey] = row.Chain1_interface
    
    new_interfaces = {}
    interface_df = read_list_table (newInterfaceFile, "Chain1_interface", int, '\t')
    for _, row in interface_df.iterrows():
        chainKey = row.Chain_1 + '-' + row.Chain_2
        if -1 in row.Chain1_interface:
            new_interfaces[chainKey] = []
        else:
            new_interfaces[chainKey] = row.Chain1_interface
    
    match = 0
    mismatch = 0
    notfound = 0
    for k, v in new_interfaces.items():
        if k in old_interfaces:
            if v == old_interfaces[k]:
                match += 1
            else:
                mismatch += 1
        else:
            notfound += 1
    
    print('\n' + 'Chain interfaces')
    print('Not found: %d' % notfound)
    print('Matches: %d' % match)
    print('Mismatches: %d' % mismatch)
    
    oldInteractomeFile = interactomeDir / 'human_interface_annotated_interactome_reserved.txt'
    newInteractomeFile = interactomeDir / 'human_interface_annotated_interactome.txt'
    
    oldInteractome = read_single_interface_annotated_interactome (oldInteractomeFile)
    newInteractome = read_single_interface_annotated_interactome (newInteractomeFile)
    
    match = 0
    mismatch = 0
    notfound = 0
    for _, row in newInteractome.iterrows():
        ind = (oldInteractome["Protein_1"] == row.Protein_1) & (oldInteractome["Protein_2"] == row.Protein_2)
        if len(ind) > 0:
            if row.Interfaces == oldInteractome.loc[ind, "Interfaces"].item():
                match += 1
            else:
                mismatch += 1
        else:
            notfound += 1
    
    print('\n' + 'Structural interactome interfaces')
    print('Not found: %d' % notfound)
    print('Matches: %d' % match)
    print('Mismatches: %d' % mismatch)
    
if __name__ == '__main__':
    main()
