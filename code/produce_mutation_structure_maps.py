import os
from pathlib import Path
from text_tools import read_list_table
from structural_annotation import write_pdb_mapped_mutations

def main():
    
    # valid reference interactome names
    interactome_names = ['HI-II-14', 'IntAct']
    
    # choose interactome (index in interactome_names)
    interactome_choise = 1
    
    download_pdbs = True
    
    suppress_pdb_warnings = False
    
    # selected reference interactome
    interactome_name = interactome_names [interactome_choise]
    
    dataDir = Path('/Volumes/MG_Samsung/junk_ppi_content/data/processed')
    
    # directory to save processed data specific to interactome
    outDir = dataDir / interactome_name
    
    # directory for PDB structure files
    pdbDir = Path('/Volumes/MG_Samsung/pdb_files')
    
    # create output directories if not existing
    if not outDir.exists():
        os.makedirs(outDir)
    
    # input data files
    natMutEdgotypeFile = outDir / 'nondisease_mutation_edgotype_geometry.txt'
    disMutEdgotypeFile = outDir / 'disease_mutation_edgotype_geometry.txt'
    ProteinSeqFile = dataDir / 'human_reference_sequences.pkl' 
    chainMapFile = outDir / 'struc_interactome_pdb_chain_map.txt'   
    chainSeqFile = dataDir / 'chain_sequences.pkl'
    chainStrucResFile = dataDir / 'chain_strucRes.pkl'
    interfaceAnnotatedInteractomeFile = outDir / 'human_interface_annotated_interactome.txt'
    chainInterfaceFile = dataDir / 'pdb_interfaces.txt'
    
    # output data files
    natural_mutations_onchains_file = outDir / 'nondisease_mutations_onchains.txt'
    disease_mutations_onchains_file = outDir / 'disease_mutations_onchains.txt'
    
    #------------------------------------------------------------------------------------
    # write edgetic mutations with structure mapping to file
    #------------------------------------------------------------------------------------
    
    naturalMutations = read_list_table (natMutEdgotypeFile,
                                        cols = ["partners", "perturbations"],
                                        dtyp = [str, int])
    diseaseMutations = read_list_table (disMutEdgotypeFile,
                                        cols = ["partners", "perturbations"],
                                        dtyp = [str, int])
    
    naturalMutations = naturalMutations[naturalMutations["Edgotype"] == 'Edgetic'].reset_index(drop=True)
    diseaseMutations = diseaseMutations[diseaseMutations["Edgotype"] == 'Edgetic'].reset_index(drop=True)
    
    # write to file edgetic non-disease mutations and perturbed protein pair after mapping 
    # mutation onto PDB chains
    if not natural_mutations_onchains_file.is_file():
        print( '\n' + 'Writing edgetic non-disease mutations with structure mapping to file ' 
                + str( natural_mutations_onchains_file ) )
        write_pdb_mapped_mutations (naturalMutations,
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
    
    # write to file edgetic disease mutations and perturbed protein pairs after mapping 
    # mutation onto PDB chains
    if not disease_mutations_onchains_file.is_file():
        print( '\n' + 'Writing edgetic disease mutations with structure mapping to file ' 
                + str( disease_mutations_onchains_file ) )
        write_pdb_mapped_mutations (diseaseMutations,
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
