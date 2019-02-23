#----------------------------------------------------------------------------------------
# This script processes mutations from the ClinVar database.
#----------------------------------------------------------------------------------------

import os
from pathlib import Path
from mutation_processing_tools import (filter_clinvar_mutations,
                                       decompose_clinvar_snp_mutations,
                                       map_clinvar_protein_refseq_IDs,
                                       get_flanking_sequences,
                                       match_flanking_sequences,
                                       remove_synon_nonsense_mutations)

def main():
    
    dataDir = Path('../data')
    
    # directory of data files from external sources
    extDir = dataDir / 'external'
    
    # directory of processed data files shared by all interactomes
    inDir = dataDir / 'processed'
    
    # directory to save output data files
    outDir = dataDir / 'processed'
        
    # create output directories if not existing
    if not outDir.exists():
        os.makedirs( outDir )
    
    # number of residues included in each side of mutation flanking sequence
    flankingSeqSideLen = 10
    
    # input data files
    allClinvarMutationsFile = extDir / 'variant_summary.txt'
    UniprotIDmapFile = inDir / 'to_human_uniprotID_map.pkl'
    rnaToProteinRefseqIDMapFile = inDir / 'human_rnaToProtein_refseqID_map.pkl'
    refseqFile = inDir / 'refseq_human_protein_sequences.txt'
    UniqueGeneSequenceFile = inDir / 'human_unique_gene_reference_sequences.txt'
    
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
    
    print( 'Decomposing mutation names' )
    decompose_clinvar_snp_mutations (clinvarMutationsFile1, clinvarMutationsFile2)

    print( 'Mapping protein refseq IDs' )
    map_clinvar_protein_refseq_IDs (clinvarMutationsFile2,
                                    rnaToProteinRefseqIDMapFile,
                                    clinvarMutationsFile3)

    print( 'Reading mutation flanking sequences from protein RefSeq transcripts' )
    get_flanking_sequences (clinvarMutationsFile3,
                            refseqFile,
                            flankingSeqSideLen,
                            clinvarMutationsFile4)

    print( 'Matching RefSeq flanking sequences to UniProt sequences' )
    match_flanking_sequences (clinvarMutationsFile4,
                              UniqueGeneSequenceFile,
                              clinvarMutationsFile5)

    print( 'Removing synonymous and nonsense mutations' )
    remove_synon_nonsense_mutations (clinvarMutationsFile5,
                                     clinvarMutationsFile6)

if __name__ == "__main__":
    main()

