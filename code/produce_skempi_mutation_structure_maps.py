
import os
import pandas as pd
from pathlib import Path
from structural_annotation import (process_skempi_mutations,
                                   write_skempi_mutation_structure_maps,
                                   write_skempi_mutation_crystal_maps)

def main():
    
    # use crystal or homology model structures 
    # options: model, crystal
    structure = 'model'
    
    filter_by_model_ddg = True
    
    # method of calculating mutation ∆∆G
    # options: bindprofx, foldx
    ddg_method = 'bindprofx'
    
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
    skempiFile = extDir / 'skempi_v2.csv'
    interactomeFile = skempiDir / 'skempi_interface_annotated_interactome.txt'
    chainMapFile = skempiDir / 'skempi_pdb_chain_map_filtered.txt'
    chainSeqFile = procDir / 'chain_sequences.pkl'
    chainStrucResFile = procDir / 'chain_strucRes.pkl'
    chainInterfaceFile = skempiDir / 'pdb_interfaces.txt'
#     modelddgFile = skempiDir / ('skempi_model_mutations_%s_ddg.txt' % ddg_method)
    modelddgFile = skempiDir / 'skempi_model_mutations_foldx_ddg.txt'
    
    # output data files
    processedSkempiFile = skempiDir / 'processed_skempi_mutations.txt'
    mutationModelMapFile = skempiDir / ('skempi_%s_mutations_%s_ddg.txt' % (structure, ddg_method))
    
    # create directories if not existing
    if not skempiDir.exists():
        os.makedirs(skempiDir)
    
    # load mutations from SKEMPI file
    mutations = pd.read_table (skempiFile, sep=';')
    
    if not processedSkempiFile.is_file():
        process_skempi_mutations (skempiFile,
                                  chainSeqFile,
                                  chainStrucResFile,
                                  pdbDir,
                                  processedSkempiFile,
                                  downloadPDB = allow_pdb_downloads,
                                  suppressWarnings = suppress_pdb_warnings)
    if structure is 'model':
        if not mutationModelMapFile.is_file():
            write_skempi_mutation_structure_maps (processedSkempiFile,
                                                  interactomeFile,
                                                  chainMapFile,
                                                  chainSeqFile,
                                                  chainStrucResFile,
                                                  pdbDir,
                                                  mutationModelMapFile,
                                                  chainInterfaceFile = chainInterfaceFile,
                                                  downloadPDB = allow_pdb_downloads,
                                                  suppressWarnings = suppress_pdb_warnings)
    elif structure is 'crystal':
        if not mutationModelMapFile.is_file():
            if filter_by_model_ddg:
                write_skempi_mutation_crystal_maps (processedSkempiFile,
                                                    mutationModelMapFile,
                                                    modelddgFile = modelddgFile)
            else:
                write_skempi_mutation_crystal_maps (processedSkempiFile, mutationModelMapFile)
        

if __name__ == '__main__':
    main()
