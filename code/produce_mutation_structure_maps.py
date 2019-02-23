import os
from pathlib import Path
from text_tools import read_list_table
from structural_annotation import write_mutation_structure_maps

def main():
    
    # reference interactome name. Options: 'HI-II-14' or 'IntAct'
    interactome_name = 'HI-II-14'
    
    # allow PDB structure files to be downloaded
    download_pdbs = True
    
    # suppress PDB warnings
    suppress_pdb_warnings = True
    
    # directory of processed data files shared by all interactomes
    dataDir = Path('../data') / 'processed'
    
    # directory of processed data files specific to interactome
    interactomeDir = dataDir / interactome_name
    
    # directory to save output data files
    outDir = dataDir / interactome_name
    
    # directory for PDB structure files
    pdbDir = Path('/Volumes/MG_Samsung/pdb_files')
    
    # create output directories if not existing
    if not outDir.exists():
        os.makedirs(outDir)
    if not interactomeDir.exists():
        os.makedirs(interactomeDir)
    
    # input data files
    ProteinSeqFile = dataDir / 'human_reference_sequences.pkl'    
    chainSeqFile = dataDir / 'chain_sequences.pkl'
    chainStrucResFile = dataDir / 'chain_strucRes.pkl'
    chainInterfaceFile = dataDir / 'pdb_interfaces.txt'
    chainMapFile = interactomeDir / 'struc_interactome_pdb_chain_map.txt'
    natMutEdgotypeFile = interactomeDir / 'nondisease_mutation_edgotype_geometry.txt'
    disMutEdgotypeFile = interactomeDir / 'disease_mutation_edgotype_geometry.txt'
    interfaceAnnotatedInteractomeFile = interactomeDir / 'human_interface_annotated_interactome.txt'
    
    # output data files
    natural_mutations_onchains_file = outDir / 'nondisease_mutations_onchains.txt'
    disease_mutations_onchains_file = outDir / 'disease_mutations_onchains.txt'
    
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
