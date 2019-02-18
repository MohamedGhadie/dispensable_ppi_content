import os
from pathlib import Path
from ddg_tools import read_foldx_results, produce_foldx_jobs, write_mutation_ddg_tofile

def main():
    
    # valid reference interactome names
    interactome_names = ['HI-II-14', 'IntAct']
    
    # choose interactome (index in interactome_names)
    interactome_choise = 0
    
    # select reference interactome
    interactome_name = interactome_names[ interactome_choise ]
    
    # directory for input data files
    dataDir = Path('/Volumes/MG_Samsung/junk_ppi_content/data/processed') / interactome_name
    
    # directory for input data files
    inDir = dataDir / 'foldx_disMut' / 'results_1'
    
    # directory for output data files
    outDir = dataDir / 'foldx_disMut'
    
    # directory for PDB structure files
    pdbDir = Path('/Volumes/MG_Samsung/pdb_files')
    
    # create output directories if not existing
    if not outDir.exists():
        os.makedirs(outDir)
    
    processed, unprocessed = read_foldx_results (inDir)
    
    write_mutation_ddg_tofile (processed,
                               dataDir / 'nondisease_mutations_onchains_foldx_ddg.txt',
                               dataDir / 'nondisease_mutations_onchains_foldx_ddg_2.txt',
                               'binding')
    write_mutation_ddg_tofile (processed,
                               dataDir / 'disease_mutations_onchains.txt',
                               dataDir / 'disease_mutations_onchains_foldx_ddg.txt',
                               'binding')
    
    produce_foldx_jobs (unprocessed,
                        pdbDir,
                        outDir,
                        write_hpc_jobfiles = True,
                        nodes = 1,
                        ppn = 1,
                        pmem = 7700,
                        walltime = '1:00:00:00',
                        rapid = 'evf-115-aa',
                        username = 'ghadie84',
                        serverDataDir = '/sf1/project/evf-115-aa/ghadie84/foldx/data')

if __name__ == "__main__":
    main()