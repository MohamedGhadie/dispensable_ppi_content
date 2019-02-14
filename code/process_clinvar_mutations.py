#----------------------------------------------------------------------------------------
# This script processes mutations from the ClinVar database.
#----------------------------------------------------------------------------------------

import os
from pathlib import Path
from text_tools import (filter_clinvar_mutations,
                        decompose_clinvar_snp_mutations,
                        map_clinvar_protein_refseq_IDs,
                        remove_synon_nonsense_mutations,
                        get_flanking_sequences,
                        match_clinvar_flanking_sequences)

def main():
    
    dataDir = Path('/Volumes/MG_Samsung/junk_ppi_content/data')
    
    # directory for data files from external sources
    inDir = dataDir / 'external'
    
    # directory to save processed data
    outDir = dataDir / 'processed'
        
    # create output directories if not existing
    if not outDir.exists():
        os.makedirs( outDir )
    
    # number of residues included in each side of mutation flanking sequence
    flankingSeqSideLen = 10
    
    # input data files
    allClinvarMutationsFile = inDir / 'variant_summary.txt'
    UniprotIDmapFile = outDir / 'to_human_uniprotID_map.pkl'
    rnaToProteinRefseqIDMapFile = outDir / 'human_rnaToProtein_refseqID_map.pkl'
    refseqFile = outDir / 'refseq_human_protein_sequences.txt'
    UniqueGeneSequenceFile = outDir / 'human_unique_gene_reference_sequences.txt'
    
    # output data files
    clinvarMutationsFile1 = outDir / 'clinvar_mutations1.txt'
    clinvarMutationsFile2 = outDir / 'clinvar_mutations2.txt'
    clinvarMutationsFile3 = outDir / 'clinvar_mutations3.txt'
    clinvarMutationsFile4 = outDir / 'clinvar_mutations4.txt'
    clinvarMutationsFile5 = outDir / 'clinvar_mutations5.txt'
    clinvarMutationsFile6 = outDir / 'clinvar_mutations6.txt'
    
    #------------------------------------------------------------------------------------
    # process disease mutations
    #------------------------------------------------------------------------------------
    
    if not clinvarMutationsFile1.is_file():
        print( 'Filtering ClinVar mutation file' )
        filter_clinvar_mutations (allClinvarMutationsFile,
                                  clinvarMutationsFile1,
                                  assembly = 'GRCh38',
                                  origin = 'germline',
                                  type = 'single nucleotide variant',
                                  incClinSig = 'Pathogenic',
                                  excClinSig = ['Pathogenic/Likely pathogenic', 'Likely pathogenic',
                                                'drug response', 'Affects', 'risk factor', 'association',
                                                'protective', 'Benign/Likely benign', 'Likely benign',
                                                'Benign', 'Uncertain significance', 'not provided',
                                                'other', 'no interpretation for the single variant',
                                                'conflicting data from submitters',
                                                'Conflicting interpretations of pathogenicity'],
                                  status = ['practice guideline', 'reviewed by expert panel',
                                            'criteria provided, multiple submitters, no conflicts',
                                            'criteria provided, single submitter'],
                                  uniprotIDmapFile = UniprotIDmapFile)
    
    if not clinvarMutationsFile2.is_file():
        print( 'Decomposing mutation names' )
        decompose_clinvar_snp_mutations (clinvarMutationsFile1, clinvarMutationsFile2)
    
    if not clinvarMutationsFile3.is_file():
        print( 'Mapping protein refseq IDs' )
        map_clinvar_protein_refseq_IDs (clinvarMutationsFile2,
                                        rnaToProteinRefseqIDMapFile,
                                        clinvarMutationsFile3)
    
    if not clinvarMutationsFile4.is_file():
        print( 'Removing synonymous and nonsense mutations' )
        remove_synon_nonsense_mutations (clinvarMutationsFile3,
                                         clinvarMutationsFile4)
    
    if not clinvarMutationsFile5.is_file():
        print( 'Reading mutation flanking sequences from protein RefSeq transcripts' )
        get_flanking_sequences (clinvarMutationsFile4,
                                refseqFile,
                                flankingSeqSideLen,
                                clinvarMutationsFile5)
    
    if not clinvarMutationsFile6.is_file():
        print( 'Matching RefSeq flanking sequences to UniProt sequences' )
        match_clinvar_flanking_sequences (clinvarMutationsFile5,
                                          UniqueGeneSequenceFile,
                                          clinvarMutationsFile6)

if __name__ == "__main__":
    main()

