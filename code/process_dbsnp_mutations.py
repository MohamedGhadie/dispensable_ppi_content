#----------------------------------------------------------------------------------------
# This script processes mutations from the dbSNP database.
#
# Run the following scripts before running this script:
# - produce_data_mappings.py
#----------------------------------------------------------------------------------------

import os
from pathlib import Path
from mutation_processing_tools import (parse_dbsnp_flatfile_keys,
                                       parse_dbsnp_flatfile,
                                       filter_and_merge_dbsnp_mutations,
                                       get_flanking_sequences,
                                       match_flanking_sequences,
                                       remove_synon_nonsense_mutations)

def main():
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # directory of dbSNP intermediate processed data files
    dbsnpDir = procDir / 'dbsnp_intermediate'
    
    # number of residues included in each side of mutation flanking sequence
    flankingSeqSideLen = 10
        
    # input data files
    UniprotIDmapFile = procDir / 'to_human_uniprotID_map.pkl'
    refseqFile = procDir / 'refseq_human_protein_sequences.txt'
    UniqueGeneSequenceFile = procDir / 'human_unique_gene_reference_sequences.txt'
    
    # output data files
    dbsnpMutationsFile1 = procDir / 'dbsnp_mutations1.txt'
    dbsnpMutationsFile2 = procDir / 'dbsnp_mutations2.txt'
    dbsnpMutationsFile3 = procDir / 'dbsnp_mutations3.txt'
    dbsnpMutationsFile4 = procDir / 'dbsnp_mutations4.txt'
    
    # create output directories if not existing
    if not procDir.exists():
        os.makedirs(procDir)
    
    #------------------------------------------------------------------------------------
    # process dbSNP mutations
    #------------------------------------------------------------------------------------
    
    print('merging dbSNP mutation files')
    filter_and_merge_dbsnp_mutations (dbsnpDir,
                                      UniprotIDmapFile,
                                      dbsnpMutationsFile1,
                                      pausetime = 30)
    
    print('reading mutation flanking sequences from RefSeq transcripts')
    get_flanking_sequences (dbsnpMutationsFile1,
                            refseqFile,
                            flankingSeqSideLen,
                            dbsnpMutationsFile2)

    print('mapping flanking sequences onto UniProt sequences')
    match_flanking_sequences (dbsnpMutationsFile2,
                              UniqueGeneSequenceFile,
                              dbsnpMutationsFile3,
                              mask = True)
    
    print('removing synonymous and nonsense mutations')
    remove_synon_nonsense_mutations (dbsnpMutationsFile3, dbsnpMutationsFile4)

if __name__ == "__main__":
    main()

