import os
import pandas as pd
from pathlib import Path
from text_tools import (parse_HI_II_14_interactome,
                        parse_IntAct_interactions)
from id_mapping import produce_protein_interaction_dict
from interactome_tools import (write_interactome_sequences,
                               remove_interactions_reported_once,
                               duplicated_PPIs,
                               remove_duplicate_PPIs)

def main():
    
    # valid reference interactome names
    interactome_names = ['HI-II-14', 'IntAct']
    
    # choose interactome (index in interactome_names)
    interactome_choise = 0
        
    # select reference interactome
    interactome_name = interactome_names[ interactome_choise ]
    
    dataDir = Path('/Volumes/MG_Samsung/junk_ppi_content/data')
    
    # directory for external input data files
    extInDir = dataDir / 'external'
    
    # directory for processed input data files
    inDir = dataDir / 'processed'
    
    # directory to save processed data specific to interactome
    outDir = dataDir / 'processed' / interactome_name
         
    # create output directories if not existing
    if not outDir.exists():
        os.makedirs(outDir)
    
    # input files from external sources
    HI_II_14_rawfile = extInDir / 'HI-II-14.tsv'
    IntAct_rawfile = extInDir / 'intact.txt'
    
    # processed input data files
    UniprotIDmapFile = inDir / 'to_human_uniprotID_map.pkl'
    UniqueGeneSwissProtIDFile = inDir / 'uniprot_unique_gene_reviewed_human_proteome.list'
    GeneMapFile = inDir / 'to_human_geneName_map.pkl'
    UniqueGeneSequenceFile = inDir / 'human_unique_gene_reference_sequences.txt'
    
    # output data files
    InteractomeFile1 = outDir / 'human_interactome_all.txt'
    InteractomeFile = outDir / 'human_interactome.txt'
    InteractomeSequenceFile = outDir / 'interactome_sequences.fasta'
    ProteinPartnersFile = outDir / 'protein_interaction_partners.pkl'
    
    if not InteractomeFile1.is_file():
        print('parsing interactome dataset')
        if interactome_name is 'HI-II-14':
            parse_HI_II_14_interactome(HI_II_14_rawfile,
                                       UniprotIDmapFile,
                                       InteractomeFile1,
                                       selfPPIs = False)
        elif interactome_name is 'IntAct':
            parse_IntAct_interactions(IntAct_rawfile, 
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
        if interactome_name == 'IntAct':
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
