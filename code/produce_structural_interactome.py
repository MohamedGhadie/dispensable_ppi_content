#----------------------------------------------------------------------------------------
# This script constructs a structural interactome from a reference interactome by mapping 
# interaction binding interfaces at amino acid resolution from experimentally determined 
# three-dimensional structural models in PDB onto interactions in the reference interactome.
#
# Run the following scripts before running this script:
# - produce_data_mappings.py
# - process_interactome.py
#----------------------------------------------------------------------------------------

import os
import pickle
import pandas as pd
from random import sample
from pathlib import Path
from pdb_tools import download_structures
from interactome_tools import (write_chain_annotated_interactome_to_excel,
                               read_single_interface_annotated_interactome,
                               write_single_interface_annotated_interactome_to_excel)
from structural_annotation import (parse_blast_file,
                                   locate_alignments,
                                   filter_chain_annotations,
                                   produce_protein_chain_dict,
                                   produce_alignment_evalue_dict,
                                   produce_chain_annotated_interactome,
                                   produce_interface_annotated_interactome,
                                   merge_interactome_interface_annotations,
                                   remove_duplicate_interface_annotations,
                                   read_chain_annotated_interactome,
                                   filter_chain_annotations_by_protein)
from plot_tools import network_plot

def main():
    
    # reference interactome name
    # options: HI-II-14, IntAct
    interactome_name = 'HI-II-14'
        
    # consider only interfaces with this minimum coverage fraction successfully mapped onto PPI
    mapCutoff = 0.5
    
    # max binding distance for interface residues in PDB structure
    bindingDist = 5
    
    # If True, merge all mapped interfaces of a PPI into one interface
    mergeInterfaces = True
    
    # max number of interfaces mapped from distinct chain-pair annotations for each PPI
    maxInterfaces = 5
    
    # max number of chain-pair interface calculations per PPI, including known interfaces
    maxAttempts = 100
    
    # If True, randomly sample from a PPI's chain pair annotations, otherwise start
    # from first annotation
    randChainPairs = False
    
    # download missing PDB structures whose chain pairs map onto interactome
    download_missing_structures = False
    
    # show figures
    showFigs = False
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # directory of data files from external sources
    extDir = dataDir / 'external'
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # directory of processed data files specific to interactome
    interactomeDir = procDir / interactome_name
    
    # figure directory
    figDir = Path('../figures') / interactome_name
        
    # directory for PDB structure files
    pdbDir = Path('/Volumes/MG_Samsung/pdb_files')
    
    if not procDir.exists():
        os.makedirs(procDir)
    if not interactomeDir.exists():
        os.makedirs(interactomeDir)
    if not figDir.exists():
        os.makedirs(figDir)
    
    # input data files
    pdbBlastFile = extDir / 'human_pdb_e-5'
    chainStrucResFile = procDir / 'chain_strucRes.pkl'
    proteinSeqFile = procDir / 'human_reference_sequences.pkl'
    pdbChainsFile = procDir / 'pdb_seqres_chains.pkl'
    chainListFile = procDir / 'pdb_seqres_chains.list'
    interactomeFile = interactomeDir / 'human_interactome.txt'
    
    # output data files
    chainMapFile1 = procDir / 'human_pdb_alignment.txt'
    chainMapFile2 = procDir / 'human_pdb_chain_map.txt'
    chainMapFile = procDir / 'human_pdb_chain_map_filtered.txt'
    proteinChainsFile = procDir / 'protein_chains.pkl'
    alignmentEvalueFile = procDir / 'human_protein_chain_min_alignment_evalues.pkl'
    chainInterfaceFile = procDir / 'pdb_interfaces.txt'
    chainAnnotatedInteractomeFile = interactomeDir / 'human_chain_annotated_interactome.txt'
    interactomeChainIDFile = interactomeDir / 'interactome_chainIDs.txt'
    interactomePDBIDFile = interactomeDir / 'interactome_pdbIDs.txt'
    interfaceAnnotatedInteractomeFile1 = interactomeDir / 'human_interface_annotated_interactome_withDuplicates.txt'
    interfaceAnnotatedInteractomeFile = interactomeDir / 'human_interface_annotated_interactome.txt'
    interactomeChainMapFile = interactomeDir / 'struc_interactome_pdb_chain_map.txt'
    chainAnnotationsFile = interactomeDir / ( interactome_name + '_PPI_chain_annotations.xlsx' )
    interfaceAnnotationsFile = interactomeDir / ( interactome_name + '_structural_interactome.xlsx' )
    
    # pausing time in seconds after processing each 50 million lines in BLAST file
    pausetime = 0
    if not chainMapFile1.is_file():
        print( 'Parsing BLAST protein-chain alignment file' )
        parse_blast_file(pdbBlastFile,
                         'us-ascii',
                         pausetime,
                         chainMapFile1)
    
    # If True, chain residue is required to match aligned protein residue for their 
    # positions to be considered aligned
    resMatch = False
    # pausing time in seconds between locating alignments on queries and subjects
    pausetime = 120
    if not chainMapFile2.is_file():
        print( 'Locating aligned residues on protein and chain sequences' )
        locate_alignments(chainMapFile1,
                          chainMapFile2,
                          resMatch = resMatch,
                          pausetime = pausetime)
    
    # Maximum e-value cutoff to filter out protein chain annotations
    evalue = 1e-10
    # Minimum chain coverage fraction required for annotation to be kept
    chainCoverage = 0
    if not chainMapFile.is_file():
        print( 'Filtering chain annotations' )
        filter_chain_annotations(chainMapFile2,
                                 evalue,
                                 chainCoverage,
                                 chainMapFile)
    
        print( 'Producing protein chains dictionary' )
        produce_protein_chain_dict(chainMapFile,
                                   proteinChainsFile)
        
        print( 'Producing protein-chain alignment evalue dictionary' )
        produce_alignment_evalue_dict (chainMapFile,
                                       alignmentEvalueFile,
                                       method = 'min')
    
    if not chainAnnotatedInteractomeFile.is_file():
        print('producing chain-annotated interactome')
        produce_chain_annotated_interactome(interactomeFile,
                                            proteinChainsFile,
                                            chainAnnotatedInteractomeFile,
                                            alignmentEvalueFile = alignmentEvalueFile)
    
    chainAnnotatedInteractome = read_chain_annotated_interactome ( chainAnnotatedInteractomeFile )
    uniqueChains = set()
    for ls in chainAnnotatedInteractome["Mapping_chains"].values:
        for pair in ls:
            uniqueChains.update( pair )
    uniqueChains = pd.Series( list( uniqueChains ) )
    uniquePDBs = uniqueChains.apply(lambda x: x.split('_')[ 0 ]).drop_duplicates()
    print( '\n' + 'Interactome chain-pair annotations:' )
    print('%d unique chains in %d unique PDB structures' % (len(uniqueChains), len(uniquePDBs)))
    uniqueChains.to_csv(interactomeChainIDFile, index=False)
    uniquePDBs.to_csv(interactomePDBIDFile, index=False)
    print('Unique chain IDs and unique PDB IDs written to file')
    
    if download_missing_structures:
        print('downloading missing structures for PDB IDs mapping onto interactome')
        download_structures( interactomeDir / 'interactome_pdbIDs.txt',
                             pdbDir )
    
    if not interfaceAnnotatedInteractomeFile1.is_file():
        print('mapping chain interfaces onto chain-annotated interactome')
        produce_interface_annotated_interactome(chainAnnotatedInteractomeFile,
                                                pdbDir,
                                                chainMapFile,
                                                chainInterfaceFile,
                                                chainStrucResFile,
                                                maxInterfaces,
                                                maxAttempts,
                                                randChainPairs,
                                                mapCutoff,
                                                bindingDist,
                                                interfaceAnnotatedInteractomeFile1)
        if mergeInterfaces:
            print('merging interface annotations for each PPI')
            merge_interactome_interface_annotations(interfaceAnnotatedInteractomeFile1,
                                                    interfaceAnnotatedInteractomeFile)
        else:
            print('removing duplicate interface annotations for each PPI without merging')
            remove_duplicate_interface_annotations(interfaceAnnotatedInteractomeFile1,
                                                   interfaceAnnotatedInteractomeFile)
    
    interactome = pd.read_table( interactomeFile )
    print( '\n' + 'Reference interactome:' )
    print( '%d PPIs' % len( interactome ) )
    print( '%d proteins' % len( set( interactome[ ["Protein_1", "Protein_2"] ].values.flatten() ) ) )
    
    with open(proteinSeqFile, 'rb') as f:
        proteinSeq = pickle.load(f)
    print( '\n' + 'Protein sequences:' )
    print( '%d sequences' % len( proteinSeq.keys() ) )
    
    with open(pdbChainsFile, 'rb') as f:
        pdbChains = pickle.load(f)
    with open(chainListFile, 'r') as f:
        chainIDs = set( f.read().split() )
    print( '\n' + 'PDB structures available:' )
    print( '%d structures' % len( pdbChains.keys() ) )
    print( '%d chains' % len( chainIDs ) )
    
    chainAnnotatedInteractome = read_chain_annotated_interactome( chainAnnotatedInteractomeFile )
    print( '\n' + 'Chain-annotated interactome:' )
    print( '%d PPIs' % len( chainAnnotatedInteractome ) )
    print( '%d proteins' % len( set( chainAnnotatedInteractome[ ["Protein_1", "Protein_2"] ].values.flatten() ) ) )
    
    write_chain_annotated_interactome_to_excel (
                    chainAnnotatedInteractome,
                    chainAnnotationsFile,
                    sheet_name = 'PPI_chain_annotations')
    
    structuralInteractome = read_single_interface_annotated_interactome( interfaceAnnotatedInteractomeFile )
    interactomeProteins = list(set(structuralInteractome[["Protein_1", "Protein_2"]].values.flatten()))
    print( '\n' + 'Structural interactome:' )
    print( '%d PPIs' % len(structuralInteractome) )
    print( '%d proteins' % len(interactomeProteins) )
    
    write_single_interface_annotated_interactome_to_excel (
                    structuralInteractome,
                    interfaceAnnotationsFile,
                    sheet_name = interactome_name + '_structural_interactome')
    
    print( '\n' + 'Plotting structural interactome' )
    edges = list(structuralInteractome[["Protein_1", "Protein_2"]].values)
    network_plot(edges,
                 show = showFigs,
                 figdir = figDir,
                 figname = 'structural_interactome')
    
    if not interactomeChainMapFile.is_file():
        print('\n' + 'filtering chain annotations by structural interactome proteins')
        filter_chain_annotations_by_protein (chainMapFile,
                                             interactomeProteins,
                                             interactomeChainMapFile)

if __name__ == "__main__":
    main()
