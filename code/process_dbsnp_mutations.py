#----------------------------------------------------------------------------------------
# This script processes mutations from the dbSNP database.
#----------------------------------------------------------------------------------------

from pathlib import Path
from mutation_processing_tools import (filter_and_merge_dbsnp_mutations,
                                       get_flanking_sequences,
                                       match_flanking_sequences,
                                       remove_synon_nonsense_mutations)

def main():
    
    dataDir = Path('../data/processed')
    
    # directory of pre-processed dbSNP mutation files
    inDir = dataDir / 'dbsnp_intermediate'
    
    # directory to save output data files
    outDir = dataDir
    
    # number of residues included in each side of mutation flanking sequence
    flankingSeqSideLen = 10
    
    # pausing time in seconds after filtering each dbSNP chromosome file
    pausetime = 60
        
    # input data files
    UniprotIDmapFile = dataDir / 'to_human_uniprotID_map.pkl'
    rnaToProteinRefseqIDMapFile = dataDir / 'human_rnaToProtein_refseqID_map.pkl'
    refseqFile = dataDir / 'refseq_human_protein_sequences.txt'
    UniqueGeneSequenceFile = dataDir / 'human_unique_gene_reference_sequences.txt'
    
    # output data files
    dbsnpMutationsFile1 = outDir / 'dbsnp_mutations1.txt'
    dbsnpMutationsFile2 = outDir / 'dbsnp_mutations2.txt'
    dbsnpMutationsFile3 = outDir / 'dbsnp_mutations3.txt'
    dbsnpMutationsFile4 = outDir / 'dbsnp_mutations4.txt'
    
    #------------------------------------------------------------------------------------
    # process dbSNP mutations
    #------------------------------------------------------------------------------------
    
    if not dbsnpMutationsFile1.is_file():
        print( 'Merging dbSNP mutation files' )
        filter_and_merge_dbsnp_mutations (inDir,
                                          UniprotIDmapFile,
                                          pausetime,
                                          dbsnpMutationsFile1)
    
    print( 'Reading mutation flanking sequences from RefSeq transcripts' )
    get_flanking_sequences (dbsnpMutationsFile1,
                            refseqFile,
                            flankingSeqSideLen,
                            dbsnpMutationsFile2)

    print( 'Mapping flanking sequences onto UniProt sequences' )
    match_flanking_sequences (dbsnpMutationsFile2,
                              UniqueGeneSequenceFile,
                              dbsnpMutationsFile3,
                              mask = True)
    
    print( 'Removing synonymous and nonsense mutations' )
    remove_synon_nonsense_mutations (dbsnpMutationsFile3,
                                     dbsnpMutationsFile4)

if __name__ == "__main__":
    main()

