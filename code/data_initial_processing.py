#----------------------------------------------------------------------------------------
# This script is for initial processing of some data files prior to being processed by 
# other scripts.
#----------------------------------------------------------------------------------------

from pathlib import Path
from text_tools import (parse_refSeq_fasta,
                        parse_dbsnp_flatfile_keys,
                        parse_dbsnp_flatfile,
                        reduce_fasta_headers)

def main():
    
    # reduce headers in PDB fasta chain sequence file
    reduce_pdb_fasta_file = False
    
    # reduce headers in UniProt fasta sequence file
    reduce_uniprot_fasta_file = False
    
    # parse RefSeq fasta sequence files and save to tables
    read_refseq_fasta_files = False
    
    # parse dbSNP files for all chromosomes and save to tables
    read_dbsnp_flatfiles = False
    
    # directory for data files from external sources
    inDir = Path('../data/external')
    
    # directory to save processed data
    outDir = Path('../data/processed')
    
    # directory for RefSeq raw data files
    refseqInDir = inDir / 'refSeq'
    
    # directory to save RefSeq processed data files
    refseqOutDir = outDir / 'refSeq_intermediate'
    
    # directory for dbSNP raw data files
    dbsnpInDir = inDir / 'dbsnp'
    
    # directory to save dbSNP processed data files
    dbsnpOutDir = outDir / 'dbsnp_intermediate'
    
    # create output directories if not existing
    if not outDir.exists():
        os.makedirs(outDir)
    if not refseqOutDir.exists():
        os.makedirs(refseqOutDir)
    if not dbsnpOutDir.exists():
        os.makedirs(dbsnpOutDir)
    
    if reduce_pdb_fasta_file:
        reduce_fasta_headers (inDir / 'pdbaa',
                              ' ',
                              1,
                              1,
                              inDir / 'pdb_reduced.fasta')
    
    if reduce_uniprot_fasta_file:
        reduce_fasta_headers (inDir / 'human_reference_sequences.fasta',
                              '|',
                              2,
                              2,
                              inDir / 'human_reference_sequences_reduced.fasta')
    
    if read_refseq_fasta_files:
        fileNum = list( np.arange(1, 23) )
        for num in fileNum:
            parse_refSeq_fasta (refseqInDir / ('human.' + str(num) + '.protein.faa'),
                                refseqOutDir / ('refseq_human_protein_' + str(num) + '.txt'))
    
    if read_dbsnp_flatfiles:
        chromosomeList = list( np.arange(1, 23) ) + ['X', 'Y']
        parse_dbsnp_flatfile_keys (dbsnpInDir / ('ds_flat_ch1.flat'),
                                   'us-ascii',
                                   150,
                                   dbsnpOutDir)
        for chromosome in chromosomeList:
            parse_dbsnp_flatfile (dbsnpInDir / ('ds_flat_ch' + str(chromosome) + '.flat'),
                                  dbsnpOutDir,
                                  'us-ascii',
                                  150,
                                  dbsnpOutDir / ('dbsnp_chr' + str(chromosome) + '.txt'))

if __name__ == "__main__":
    main()
