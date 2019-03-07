#----------------------------------------------------------------------------------------
# This script processes results from foldx calculations. Since each foldx result file 
# contains only one mutation, no second-round jobs will be produced.
#
# Requirements:
# Files must be in format produced by script produce_mutation_structure_maps.py
#----------------------------------------------------------------------------------------

import os
from pathlib import Path
from ddg_tools import read_foldx_results, produce_foldx_jobs, write_mutation_ddg_tofile

def main():
    
    # reference interactome name
    # options: HI-II-14, IntAct
    interactome_name = 'IntAct'
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # directory of processed data files specific to interactome
    interactomeDir = procDir / interactome_name
    
    # directory of foldx results
    inDir = interactomeDir / 'foldx' / 'results_all_2'
    
    # directory of foldx output jobs
    outDir = interactomeDir / 'foldx'
    
    # directory of PDB structure files
    pdbDir = Path('/Volumes/MG_Samsung/pdb_files')
    
    # create output directories if not existing
    if not procDir.exists():
        os.makedirs(procDir)
    if not interactomeDir.exists():
        os.makedirs(interactomeDir)
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
