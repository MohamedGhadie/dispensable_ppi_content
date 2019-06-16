#----------------------------------------------------------------------------------------
# Produce structural models for PPIs in SKEMPI database.
#
# Run the following scripts before running this script:
# - 
# - 
# - 
#----------------------------------------------------------------------------------------

import os
from pathlib import Path
from text_tools import parse_blast_file
from pdb_tools import download_structures
from id_mapping import produce_protein_chain_dict
from interactome_tools import read_chain_annotated_interactome, read_single_interface_annotated_interactome
from structural_annotation import (locate_alignments,
                                   filter_skempi_chain_annotations,
                                   produce_alignment_evalue_dict,
                                   produce_chain_annotated_interactome,
                                   produce_interface_annotated_interactome,
                                   merge_interactome_interface_annotations,
                                   remove_duplicate_interface_annotations)

def main():
    
    # consider only interfaces with this minimum fraction mapped onto PPI
    mapCutoff = 0.5
    
    # max binding distance for interface residues in PDB structure
    bindingDist = 5
    
    # consider only models with this minimum coverage for proteins and chains
    minCov = 0.5
    
    # consider only models with this maximun coverage for proteins and chains
    maxCov = 0.9
    
    # download missing PDB structures whose chain pairs map onto interactome
    download_missing_structures = True
    
    # allow downloading of PDB structures while constructing the structural interactome
    allow_pdb_downloads = True
    
    # suppress PDB warnings when constructing the structural interactome
    suppress_pdb_warnings = True
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # directory of data files from external sources
    extDir = dataDir / 'external'
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    skempiDir = procDir / 'skempi'
    
    # directory for PDB structure files
    pdbDir = Path('/Volumes/MG_Samsung/pdb_files')
    
    # input data files
    pdbBlastFile = extDir / 'skempi_pdb_e1'
    chainSeqFile = procDir / 'chain_sequences.pkl'
    chainStrucResFile = procDir / 'chain_strucRes.pkl'
    interactomeFile = skempiDir / 'skempi_interactome.txt'
    
    # output data files
    chainMapFile1 = skempiDir / 'skempi_pdb_alignment.txt'
    chainMapFile2 = skempiDir / 'skempi_pdb_chain_map.txt'
    chainMapFile3 = skempiDir / 'skempi_pdb_chain_map_filtered.txt'
    proteinChainsFile = skempiDir / 'protein_chains.pkl'
    alignmentEvalueFile = skempiDir / 'skempi_protein_chain_min_alignment_evalues.pkl'
    chainInterfaceFile = skempiDir / 'pdb_interfaces.txt'
    chainAnnotatedInteractomeFile = skempiDir / 'skempi_chain_annotated_interactome.txt'
    chainIDFile = skempiDir / 'interactome_chainIDs.txt'
    pdbIDFile = skempiDir / 'interactome_pdbIDs.txt'
    interfaceAnnotatedInteractomeFile1 = skempiDir / 'skempi_interface_annotated_interactome_withDuplicates.txt'
    interfaceAnnotatedInteractomeFile = skempiDir / 'skempi_interface_annotated_interactome.txt'
    
    # create directories if not existing
    if not skempiDir.exists():
        os.makedirs(skempiDir)
    if not pdbDir.exists():
        os.makedirs(pdbDir)
    
    if not chainMapFile1.is_file():
        print('parsing BLAST protein-chain alignment file')
        parse_blast_file (pdbBlastFile, chainMapFile1)
    
    if not chainMapFile2.is_file():
        print('locating aligned residues on protein and chain sequences')
        locate_alignments (chainMapFile1,
                           chainMapFile2,
                           resMatch = False,
                           pausetime = 0)
    
    if not chainMapFile3.is_file():
        filter_skempi_chain_annotations (chainMapFile2,
                                         chainMapFile3,
                                         minCov = minCov,
                                         maxCov = maxCov)

    if not proteinChainsFile.is_file():
        print('producing protein chains dictionary')
        produce_protein_chain_dict (chainMapFile3, proteinChainsFile)
    
    if not alignmentEvalueFile.is_file():
        print('producing protein-chain alignment evalue dictionary')
        produce_alignment_evalue_dict (chainMapFile3, alignmentEvalueFile, method = 'min')
    
    if not chainAnnotatedInteractomeFile.is_file():
        print('producing chain-annotated interactome')
        produce_chain_annotated_interactome (interactomeFile,
                                             proteinChainsFile,
                                             chainAnnotatedInteractomeFile,
                                             alignmentEvalueFile = alignmentEvalueFile)
    
    chainAnnotatedInteractome = read_chain_annotated_interactome (chainAnnotatedInteractomeFile)
    interactomeProteins = list(set(chainAnnotatedInteractome[["Protein_1", "Protein_2"]].values.flatten()))
    print('\n' + 'Chain-annotated interactome:')
    print('%d PPIs' % len(chainAnnotatedInteractome))
    print('%d proteins' % len(interactomeProteins))
    print()
    
    uniqueChains = set()
    for ls in chainAnnotatedInteractome["Mapping_chains"].values:
        for pair in ls:
            uniqueChains.update(pair)
    uniquePDBs = {id.split('_')[0] for id in uniqueChains}
    print('\n' + 'Interactome chain-pair annotations:')
    print('%d unique chains in %d unique PDB structures' % (len(uniqueChains), len(uniquePDBs)))
    
    with open(chainIDFile, 'w') as f:
        for i in sorted(uniqueChains):
            f.write("%s\n" % i)
    with open(pdbIDFile, 'w') as f:
        for i in sorted(uniquePDBs):
            f.write("%s\n" % i)
    
    if download_missing_structures:
        print('downloading missing structures for PDB IDs mapping onto interactome')
        download_structures (pdbIDFile, pdbDir)
    
    if not interfaceAnnotatedInteractomeFile1.is_file():
        print('mapping chain interfaces onto chain-annotated interactome')
        # merge all mapped interfaces of a PPI into one interface
        mergeInterfaces = True
        # max number of interfaces mapped from distinct chain-pair annotations for each PPI
        maxInterfaces = 1
        # max number of chain-pair interface calculations per PPI, including known interfaces
        maxAttempts = 100
        # if True, randomly sample from a PPI's chain pair annotations, otherwise start
        # from first annotation
        randChainPairs = False
        # produce structural interactome
        produce_interface_annotated_interactome (chainAnnotatedInteractomeFile,
                                                 pdbDir,
                                                 chainSeqFile,
                                                 chainMapFile3,
                                                 chainInterfaceFile,
                                                 chainStrucResFile,
                                                 maxInterfaces,
                                                 maxAttempts,
                                                 randChainPairs,
                                                 mapCutoff,
                                                 bindingDist,
                                                 interfaceAnnotatedInteractomeFile1,
                                                 downloadPDB = allow_pdb_downloads,
                                                 suppressWarnings = suppress_pdb_warnings)
        if mergeInterfaces:
            print('merging interface annotations for each PPI')
            merge_interactome_interface_annotations (interfaceAnnotatedInteractomeFile1,
                                                     interfaceAnnotatedInteractomeFile)
        else:
            print('removing duplicate interface annotations for each PPI without merging')
            remove_duplicate_interface_annotations (interfaceAnnotatedInteractomeFile1,
                                                    interfaceAnnotatedInteractomeFile)
    
    structuralInteractome = read_single_interface_annotated_interactome (interfaceAnnotatedInteractomeFile)
    interactomeProteins = list(set(structuralInteractome[["Protein_1", "Protein_2"]].values.flatten()))
    print('\n' + 'Structural interactome:')
    print('%d PPIs' % len(structuralInteractome))
    print('%d proteins' % len(interactomeProteins))
    print()
    
if __name__ == '__main__':
    main()
