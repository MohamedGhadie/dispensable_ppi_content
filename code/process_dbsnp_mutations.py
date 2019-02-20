#----------------------------------------------------------------------------------------
# This script processes mutations from the ClinVar database.
#----------------------------------------------------------------------------------------

from pathlib import Path
from text_tools import (filter_and_merge_dbsnp_mutations,
                        get_flanking_sequences,
                        remove_dbsnp_synon_nonsense_mutations,
                        match_masked_flanking_sequences)

def main():
    
    dataDir = Path('/Volumes/MG_Samsung/junk_ppi_content/data/processed')
    
    # directory for data files from external sources
    inDir = dataDir / 'dbsnp_intermediate'
    
    # directory to save processed data
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
    
    print( 'Reading dbSNP mutation flanking sequences from RefSeq transcripts' )
    get_flanking_sequences (dbsnpMutationsFile1,
                            refseqFile,
                            flankingSeqSideLen,
                            dbsnpMutationsFile2)
    
    print( 'Removing dbSNP synonymous and nonsense mutations' )
    remove_dbsnp_synon_nonsense_mutations (dbsnpMutationsFile2,
                                           dbsnpMutationsFile3)

    print( 'Mapping dbSNP flanking sequences onto UniProt sequences' )
    match_masked_flanking_sequences (dbsnpMutationsFile3,
                                     UniqueGeneSequenceFile,
                                     dbsnpMutationsFile4)

if __name__ == "__main__":
    main()

