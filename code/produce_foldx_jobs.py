import os
from pathlib import Path
from ddg_tools import read_unprocessed_ddg_mutations, produce_foldx_jobs

def main():
    
    # valid reference interactome names
    interactome_names = ['HI-II-14', 'IntAct']
    
    # choose interactome (index in interactome_names)
    interactome_choise = 0
    
    # select reference interactome
    interactome_name = interactome_names[ interactome_choise ]
    
    # directory for input data files
    dataDir = Path('/Volumes/MG_Samsung/junk_ppi_content/data/processed') / interactome_name
    
    # directory for output data files
    outDir = dataDir / 'foldx_disMut'
    
    # directory for PDB structure files
    pdbDir = Path('/Volumes/MG_Samsung/pdb_files')
    
    # input file containing mutations to submit to bindprofx
    mutationsFile = dataDir / 'disease_mutations_onchains_foldx_ddg.txt'
    
    # create output directories if not existing
    if not outDir.exists():
        os.makedirs( outDir )
    
    mutations = read_unprocessed_ddg_mutations (mutationsFile, 'binding')
    
    produce_foldx_jobs (mutations,
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