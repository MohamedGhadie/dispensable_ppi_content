import os
import pickle
import pandas as pd
from random import sample
from pathlib import Path
from pdb_tools import download_structures
from interactome_tools import (write_chain_annotated_interactome_to_excel,
                               read_single_interface_annotated_interactome,
                               write_single_interface_annotated_interactome_to_excel)
from structural_annotation import (produce_chain_annotated_interactome,
                                   produce_interface_annotated_interactome,
                                   merge_interactome_interface_annotations,
                                   remove_duplicate_interface_annotations,
                                   read_chain_annotated_interactome,
                                   filter_chain_annotations_by_protein)
from plot_tools import network_plot

def main():
    
    # valid reference interactome names
    interactome_names = ['HI-II-14', 'IntAct']
    
    # choose interactome (index in interactome_names)
    interactome_choise = 1
    
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
    
    # show figures
    showFigs = False
    
    # download missing PDB structures whose chain pairs map onto interactome
    download_missing_structures = False
        
    # select reference interactome
    interactome_name = interactome_names[ interactome_choise ]
    
    dataDir = Path('/Volumes/MG_Samsung/junk_ppi_content/data')
    
    # input data directory
    inDir = dataDir / 'processed'
    
    # directory to save processed data specific to interactome
    outDir = dataDir / 'processed' / interactome_name
    
    # figure directory
    figDir = Path('../figures') / interactome_name
        
    # directory for PDB structure files
    pdbDir = Path('/Volumes/MG_Samsung/pdb_files')
    
    if not figDir.exists():
        os.makedirs(figDir)
    if not outDir.exists():
        os.makedirs(outDir)
    
    # input data files
    chainMapFile = inDir / 'human_pdb_chain_map_filtered.txt'
    proteinChainsFile = inDir / 'protein_chains.pkl'
    chainStrucResFile = inDir / 'chain_strucRes.pkl'
    ProteinSeqFile = inDir / 'human_reference_sequences.pkl'
    pdbChainsFile = inDir / 'pdb_seqres_chains.pkl'
    chainListFile = inDir / 'pdb_seqres_chains.list'
    alignmentEvalueFile = inDir / 'human_protein_chain_min_alignment_evalues.pkl'
    InteractomeFile = outDir / 'human_interactome.txt'
    
    # output data files
    chainAnnotatedInteractomeFile = outDir / 'human_chain_annotated_interactome.txt'
    interfaceAnnotatedInteractomeFile1 = outDir / 'human_interface_annotated_interactome_withDuplicates.txt'
    interfaceAnnotatedInteractomeFile = outDir / 'human_interface_annotated_interactome.txt'
    interactomeChainMapFile = outDir / 'struc_interactome_pdb_chain_map.txt'
    chainAnnotationsFile = outDir / ( interactome_name + '_PPI_chain_annotations.xlsx' )
    interfaceAnnotationsFile = outDir / ( interactome_name + '_structural_interactome.xlsx' )
    chainInterfaceFile = inDir / 'pdb_interfaces.txt'
    
    if not chainAnnotatedInteractomeFile.is_file():
        print('producing chain-annotated interactome')
        produce_chain_annotated_interactome(InteractomeFile,
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
    uniqueChains.to_csv(outDir / 'interactome_chainIDs.txt', index=False)
    uniquePDBs.to_csv(outDir / 'interactome_pdbIDs.txt', index=False)
    print('Unique chain IDs and unique PDB IDs written to file')
    
    if download_missing_structures:
        print('downloading missing structures for PDB IDs mapping onto interactome')
        download_structures( interactomeOutDir / 'interactome_pdbIDs.txt',
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
    
    interactome = pd.read_table( InteractomeFile )
    print( '\n' + 'Reference interactome:' )
    print( '%d PPIs' % len( interactome ) )
    print( '%d proteins' % len( set( interactome[ ["Protein_1", "Protein_2"] ].values.flatten() ) ) )
    
    with open(ProteinSeqFile, 'rb') as f:
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
