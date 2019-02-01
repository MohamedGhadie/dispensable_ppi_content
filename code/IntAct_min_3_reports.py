import pandas as pd
from pathlib import Path
from interactome_tools import (read_single_interface_annotated_interactome,
                               write_single_interface_annotated_interactome)

def main():
    
    # directory where interactome processed data is saved
    interactomeDir = Path('../data/processed/IntAct')
    
    interactomeDupFile = interactomeDir / 'human_interactome_all.txt'
    interactomeDup = pd.read_table( interactomeDupFile )
    print( '\n' + 'Size of reference interactome with duplicates: %d' % len( interactomeDup ) )
    
    structuralInteractomeFile = interactomeDir / 'human_interface_annotated_interactome.txt'
    structuralInteractome = read_single_interface_annotated_interactome( structuralInteractomeFile )
    print( '\n' + 'Size of structural interactome: %d' % len( structuralInteractome ) )
    
    dupPPIstrings = interactomeDup.apply(lambda x:
                            '-'.join( sorted( [ x["Protein_1"], x["Protein_2"] ] ) ), axis=1)
    strucPPIstrings = structuralInteractome.apply(lambda x:
                            '-'.join( sorted( [ x["Protein_1"], x["Protein_2"] ] ) ), axis=1)
    
    counts = dupPPIstrings.value_counts()
    dup = strucPPIstrings.apply(lambda x: counts[ x ] > 2 if x in counts else 0)
    
    structuralInteractome = structuralInteractome[ dup ].reset_index( drop = True )
    
    outPath = interactomeDir / 'human_interface_annotated_interactome_min3reports.txt'
    write_single_interface_annotated_interactome ( structuralInteractome, outPath )

if __name__ == "__main__":
    main()
