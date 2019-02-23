import os
from pathlib import Path
from ddg_tools import read_foldx_results, produce_foldx_jobs, write_mutation_ddg_tofile

def main():
    
    # reference interactome name. Options: 'HI-II-14' or 'IntAct'
    interactome_name = 'IntAct'
    
    # directory for input data files
    dataDir = Path('../data') / 'processed'
    
    # directory of processed data files specific to interactome
    interactomeDir = dataDir / interactome_name
    
    # directory for input data files
    inDir = interactomeDir / 'foldx' / 'results_all_2'
    
    # directory for output data files
    outDir = interactomeDir / 'foldx'
    
    # directory for PDB structure files
    pdbDir = Path('/Volumes/MG_Samsung/pdb_files')
    
    # create output directories if not existing
    if not outDir.exists():
        os.makedirs(outDir)
    
    processed, unprocessed = read_foldx_results (inDir)
    
    write_mutation_ddg_tofile (processed,
                               interactomeDir / 'nondisease_mutations_foldx_ddg.txt',
                               interactomeDir / 'nondisease_mutations_foldx_ddg_2.txt',
                               'binding')
    write_mutation_ddg_tofile (processed,
                               interactomeDir / 'disease_mutations_foldx_ddg.txt',
                               interactomeDir / 'disease_mutations_foldx_ddg_2.txt',
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
