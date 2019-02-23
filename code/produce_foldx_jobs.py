import os
from pathlib import Path
from ddg_tools import read_unprocessed_ddg_mutations, produce_foldx_jobs

def main():
    
    # reference interactome name. Options: 'HI-II-14' or 'IntAct'
    interactome_name = 'IntAct'
    
    # directory for input data files
    dataDir = Path('../data') / 'processed'
    
    # directory of processed data files specific to interactome
    interactomeDir = dataDir / interactome_name
    
    # directory for output data files
    outDir = interactomeDir / 'foldx'
    
    # directory for PDB structure files
    pdbDir = Path('/Volumes/MG_Samsung/pdb_files')
    
    # input file containing mutations to submit to bindprofx
    mutationsFile = interactomeDir / 'disease_mutations_foldx_ddg.txt'
    
    # create output directories if not existing
    if not outDir.exists():
        os.makedirs(outDir)
    
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
