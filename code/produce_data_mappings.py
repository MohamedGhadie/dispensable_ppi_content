#----------------------------------------------------------------------------------------
# Process several ID and sequence files.
#----------------------------------------------------------------------------------------

import os
from pathlib import Path
from text_tools import parse_fasta_file, reduce_fasta_headers
from id_mapping import (produce_geneName_dict,
                        produce_uniqueGene_swissProtIDs,
                        produce_uniqueGene_sequences,
                        produce_proteinSeq_dict,
                        produce_uniprotID_dict,
                        produce_rnaToProtein_refseqID_dict,
                        produce_chainSeq_dict,
                        produce_chain_dict,
                        produce_chain_strucRes_dict)
from pdb_tools import produce_chain_list, produce_pdb_ids

def main():
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # directory of data files from external sources
    extDir = Path('/Volumes/MG_Samsung/junk_ppi_content/data/external')
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # input data files
    uniprotRefSeqFile = extDir / 'UP000005640_9606.fasta'
    idmapFile = extDir / 'HUMAN_9606_idmapping.dat'
    proteomeListFile = extDir / 'uniprot_reviewed_human_proteome.list'
    refseqIDFile = extDir / 'LRG_RefSeqGene'
    pdbSeqresFile = extDir / 'pdb_seqres.txt'
    chainResAnnotFile = extDir / 'ss_dis.txt'
    
    # output data files
    refSeqFile = extDir / 'human_reference_sequences.fasta'
    SequenceFile = procDir / 'human_reference_sequences.txt'
    GeneMapFile = procDir / 'to_human_geneName_map.pkl'
    UniqueGeneSwissProtIDFile = procDir / 'uniprot_unique_gene_reviewed_human_proteome.list'
    UniqueGeneSequenceFile = procDir / 'human_unique_gene_reference_sequences.txt'
    ProteinSeqFile = procDir / 'human_reference_sequences.pkl'
    UniprotIDmapFile = procDir / 'to_human_uniprotID_map.pkl'
    rnaToProteinRefseqIDMapFile = procDir / 'human_rnaToProtein_refseqID_map.pkl'
    seqresFile = extDir / 'pdb_seqres_reduced.fasta'
    chainSeqFile = procDir / 'chain_sequences.pkl'
    chainListFile = procDir / 'pdb_seqres_chains.list'
    pdbidFile = procDir / 'seqres_pdb_IDs.list'
    pdbChainsFile = procDir / 'pdb_seqres_chains.pkl'
    chainStrucResFile = procDir / 'chain_strucRes.pkl'
    
    # create output directories if not existing
    if not procDir.exists():
        os.makedirs(procDir)
    
    #------------------------------------------------------------------------------------
    # process data files
    #------------------------------------------------------------------------------------
    
    if not refSeqFile.is_file():
        print('Reducing headers in protein sequence fasta file')
        reduce_fasta_headers (uniprotRefSeqFile, '|', 2, 2, refSeqFile)
    
    if not SequenceFile.is_file():
        print('reading protein sequence fasta file')
        parse_fasta_file (refSeqFile, SequenceFile)
    
    if not GeneMapFile.is_file():
        print('producing UniProtID-to-geneName dictionary')
        produce_geneName_dict (idmapFile, proteomeListFile, GeneMapFile)
    
    if not UniqueGeneSwissProtIDFile.is_file():
        print('producing list of unique-gene UniProt IDs')
        produce_uniqueGene_swissProtIDs (proteomeListFile, GeneMapFile, UniqueGeneSwissProtIDFile)
    
    if not UniqueGeneSequenceFile.is_file():
        print('producing sequence file for unique-gene UniProt IDs')
        produce_uniqueGene_sequences (SequenceFile, UniqueGeneSwissProtIDFile, GeneMapFile, UniqueGeneSequenceFile)
    
    if not ProteinSeqFile.is_file():
        print('producing protein sequence dictionary')
        produce_proteinSeq_dict (refSeqFile, ProteinSeqFile)
    
    if not UniprotIDmapFile.is_file():
        print('producing to-UniProt-ID dictionary')
        produce_uniprotID_dict (idmapFile, UniqueGeneSwissProtIDFile, UniprotIDmapFile)
    
    if not rnaToProteinRefseqIDMapFile.is_file():
        print('producing rna to protein RefSeq ID dictionary')
        produce_rnaToProtein_refseqID_dict (refseqIDFile, rnaToProteinRefseqIDMapFile)
    
    if not seqresFile.is_file():
        print('reducing headers in PDB chain sequence file')
        reduce_fasta_headers (pdbSeqresFile, ' ', 1, 1, seqresFile)
    
    if not chainSeqFile.is_file():
        print('producing PDB chain sequence dictionary from fasta records')
        produce_chainSeq_dict (seqresFile, chainSeqFile)
    
    if not chainListFile.is_file():
        print('producing PDB chain ID file from fasta records')
        produce_chain_list (seqresFile, chainListFile)
    
    if not pdbidFile.is_file():
        print('producing unique PDB ID file from chain list file')
        produce_pdb_ids (chainListFile, pdbidFile)
    
    if not pdbChainsFile.is_file():    
        print('producing PDB chain dictionary from chain list file')
        produce_chain_dict (chainListFile, pdbChainsFile)
    
    if not chainStrucResFile.is_file():
        print('parsing PDB chain structured-residue order file')
        produce_chain_strucRes_dict (chainResAnnotFile, chainStrucResFile)

if __name__ == "__main__":
    main()
