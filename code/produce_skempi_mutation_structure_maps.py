#----------------------------------------------------------------------------------------
# Map single-point mutations in SKEMPI onto structures to be submitted for 
# ∆∆G calculations.
#
# Run the following scripts before running this script:
# - produce_data_mappings.py
# - process_skempi_file.py
#----------------------------------------------------------------------------------------

import os
from pathlib import Path
from structural_annotation import write_skempi_mutation_structure_maps, write_skempi_mutation_crystal_maps

def main():
    
    # use crystal or homology model structures 
    # options: model, crystal
    structure = 'model'
    
    filter_by_model_ddg = False
    
    # method of calculating mutation ∆∆G
    # options: bindprofx, foldx
    ddg_method = 'foldx'
    
    # allow downloading of PDB structures while constructing the structural interactome
    allow_pdb_downloads = True
    
    # suppress PDB warnings when constructing the structural interactome
    suppress_pdb_warnings = True
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    skempiDir = procDir / 'skempi'
    
    # directory for PDB structure files
    pdbDir = Path('../pdb_files')
    
    # input data files
    processedSkempiFile = skempiDir / 'processed_skempi_mutations.txt'
    interactomeFile = skempiDir / 'skempi_interface_annotated_interactome.txt'
    chainMapFile = skempiDir / 'skempi_pdb_chain_map_filtered.txt'
    chainSeqFile = procDir / 'chain_sequences.pkl'
    chainStrucResFile = procDir / 'chain_strucRes.pkl'
    chainInterfaceFile = procDir / 'pdb_interfaces.txt'
#     modelddgFile = skempiDir / ('skempi_model_mutations_%s_ddg.txt' % ddg_method)
    modelddgFile = skempiDir / 'skempi_model_mutations_foldx_ddg.txt'
    
    # output data files
    mutationMapFile = skempiDir / ('skempi_%s_mutations_%s_ddg.txt' % (structure, ddg_method))
    
    # create directories if not existing
    if not skempiDir.exists():
        os.makedirs(skempiDir)
    if not pdbDir.exists():
        os.makedirs(pdbDir)
    
    if structure is 'model':
        if not mutationMapFile.is_file():
            if not chainInterfaceFile.is_file():
                chainInterfaceFile = None
            write_skempi_mutation_structure_maps (processedSkempiFile,
                                                  interactomeFile,
                                                  chainMapFile,
                                                  chainSeqFile,
                                                  chainStrucResFile,
                                                  pdbDir,
                                                  mutationMapFile,
                                                  chainInterfaceFile = chainInterfaceFile,
                                                  downloadPDB = allow_pdb_downloads,
                                                  suppressWarnings = suppress_pdb_warnings)
    elif structure is 'crystal':
        if not mutationMapFile.is_file():
            if filter_by_model_ddg and modelddgFile.is_file():
                write_skempi_mutation_crystal_maps (processedSkempiFile,
                                                    mutationMapFile,
                                                    modelddgFile = modelddgFile)
            else:
                write_skempi_mutation_crystal_maps (processedSkempiFile, mutationMapFile)

if __name__ == '__main__':
    main()
