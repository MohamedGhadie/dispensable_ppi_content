#----------------------------------------------------------------------------------------
# Map mutations in the structural interactome onto structural models to be submitted for 
# ∆∆G calculations.
#
# Run the following scripts before running this script:
# - produce_data_mappings.py
# - produce_structural_interactome.py
# - calculate_junk_content_geometry.py
#----------------------------------------------------------------------------------------

import os
from pathlib import Path
from text_tools import read_list_table
from structural_annotation import write_mutation_structure_maps

def main():
    
    # reference interactome name
    # options: HI-II-14, IntAct
    interactome_name = 'HI-II-14'
    
    # allow PDB structure files to be downloaded
    download_pdbs = True
    
    # suppress PDB warnings
    suppress_pdb_warnings = True
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # directory of processed data files specific to interactome
    interactomeDir = procDir / interactome_name
    
    # directory for PDB structure files
    pdbDir = Path('/Volumes/MG_Samsung/pdb_files')
    
    # input data files
    ProteinSeqFile = procDir / 'human_reference_sequences.pkl'    
    chainSeqFile = procDir / 'chain_sequences.pkl'
    chainStrucResFile = procDir / 'chain_strucRes.pkl'
    chainInterfaceFile = procDir / 'pdb_interfaces.txt'
    chainMapFile = interactomeDir / 'struc_interactome_pdb_chain_map.txt'
    natMutEdgotypeFile = interactomeDir / 'nondisease_mutation_edgotype_geometry.txt'
    disMutEdgotypeFile = interactomeDir / 'disease_mutation_edgotype_geometry.txt'
    interfaceAnnotatedInteractomeFile = interactomeDir / 'human_interface_annotated_interactome.txt'
    
    # output data files
    natural_mutations_onchains_file = interactomeDir / 'nondisease_mutations_onchains.txt'
    disease_mutations_onchains_file = interactomeDir / 'disease_mutations_onchains.txt'
    
    # create output directories if not existing
    if not procDir.exists():
        os.makedirs(procDir)
    if not interactomeDir.exists():
        os.makedirs(interactomeDir)
    
    #------------------------------------------------------------------------------------
    # write edgetic mutations mapped onto PPI structural models to file
    #------------------------------------------------------------------------------------
    
    naturalMutations = read_list_table (natMutEdgotypeFile,
                                        cols = ["partners", "perturbations"],
                                        dtyp = [str, int])
    diseaseMutations = read_list_table (disMutEdgotypeFile,
                                        cols = ["partners", "perturbations"],
                                        dtyp = [str, int])
    
    naturalMutations = naturalMutations[naturalMutations["Edgotype"] == 'Edgetic'].reset_index(drop=True)
    diseaseMutations = diseaseMutations[diseaseMutations["Edgotype"] == 'Edgetic'].reset_index(drop=True)
    
    # write edgetic non-disease mutations
    if not natural_mutations_onchains_file.is_file():
        print( '\n' + 'Writing edgetic non-disease mutations with structure mapping to file ' 
                + str( natural_mutations_onchains_file ) )
        write_mutation_structure_maps (naturalMutations,
                                       interfaceAnnotatedInteractomeFile,
                                       chainMapFile,
                                       chainSeqFile,
                                       ProteinSeqFile,
                                       chainStrucResFile,
                                       pdbDir,
                                       natural_mutations_onchains_file,
                                       chainInterfaceFile = chainInterfaceFile,
                                       downloadPDB = download_pdbs,
                                       suppressWarnings = suppress_pdb_warnings)
    
    # write edgetic disease mutations
    if not disease_mutations_onchains_file.is_file():
        print( '\n' + 'Writing edgetic disease mutations with structure mapping to file ' 
                + str( disease_mutations_onchains_file ) )
        write_mutation_structure_maps (diseaseMutations,
                                       interfaceAnnotatedInteractomeFile,
                                       chainMapFile,
                                       chainSeqFile,
                                       ProteinSeqFile,
                                       chainStrucResFile,
                                       pdbDir,
                                       disease_mutations_onchains_file,
                                       chainInterfaceFile = chainInterfaceFile,
                                       downloadPDB = download_pdbs,
                                       suppressWarnings = suppress_pdb_warnings)

if __name__ == "__main__":
    main()
