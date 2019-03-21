# Process reference interactome from protein-protein interaction dataset.
#
# Run the following scripts before running this script:
# - produce_data_mappings.py
#----------------------------------------------------------------------------------------

import os
import pandas as pd
from pathlib import Path
from text_tools import parse_HI_II_14_interactome, parse_IntAct_interactions
from id_mapping import produce_protein_interaction_dict
from interactome_tools import (write_interactome_sequences,
                               remove_interactions_reported_once,
                               duplicated_PPIs,
                               remove_duplicate_PPIs)

def main():
    
    # reference interactome name
    # options: HI-II-14, IntAct
    interactome_name = 'IntAct'
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # directory of data files from external sources
    extDir = Path('/Volumes/MG_Samsung/junk_ppi_content/data/external')
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # directory of processed data files specific to interactome
    interactomeDir = procDir / interactome_name
    
    # input data files
    HI_II_14_datafile = extDir / 'HI-II-14.tsv'
    IntAct_datafile = extDir / 'intact.txt'
    UniprotIDmapFile = procDir / 'to_human_uniprotID_map.pkl'
    UniqueGeneSwissProtIDFile = procDir / 'uniprot_unique_gene_reviewed_human_proteome.list'
    GeneMapFile = procDir / 'to_human_geneName_map.pkl'
    UniqueGeneSequenceFile = procDir / 'human_unique_gene_reference_sequences.txt'
    
    # output data files
    InteractomeFile1 = interactomeDir / 'human_interactome_all.txt'
    InteractomeFile = interactomeDir / 'human_interactome.txt'
    InteractomeSequenceFile = interactomeDir / 'interactome_sequences.fasta'
    ProteinPartnersFile = interactomeDir / 'protein_interaction_partners.pkl'
    
    # create output directories if not existing
    if not procDir.exists():
        os.makedirs(procDir)
    if not interactomeDir.exists():
        os.makedirs(interactomeDir)
    
    if not InteractomeFile1.is_file():
        print('parsing interactome dataset')
        if interactome_name is 'HI-II-14':
            parse_HI_II_14_interactome(HI_II_14_datafile,
                                       UniprotIDmapFile,
                                       InteractomeFile1,
                                       selfPPIs = False)
        elif interactome_name is 'IntAct':
            parse_IntAct_interactions(IntAct_datafile, 
                                      UniqueGeneSwissProtIDFile,
                                      InteractomeFile1,
                                      geneMapFile = GeneMapFile,
                                      selfPPIs = False)
        else:
            print('Error: interactome name %s not recognized. Quiting...' % interactome_name)
            return
    
    if not InteractomeFile.is_file():
        interactome = pd.read_table(InteractomeFile1)
        print('Initial interactome size: %d PPIs' % len(interactome))
        if interactome_name is 'IntAct':
            interactome = remove_interactions_reported_once(interactome)
            #interactome = duplicated_PPIs(interactome)
            print('Interactome size after removing PPIs reported only once: %d PPIs' % len(interactome))
        interactome = remove_duplicate_PPIs(interactome)
        print('Interactome size after removing duplicate PPIs: %d PPIs' % len(interactome))
        interactome.to_csv(InteractomeFile, index=False, sep='\t')
    
    if not InteractomeSequenceFile.is_file():
        print('writing interactome protein sequences')
        write_interactome_sequences(InteractomeFile,
                                    UniqueGeneSequenceFile,
                                    InteractomeSequenceFile)
    
    if not ProteinPartnersFile.is_file():
        print('producing protein interaction partners dictionary')
        produce_protein_interaction_dict(InteractomeFile,
                                         ProteinPartnersFile)

if __name__ == "__main__":
    main()
