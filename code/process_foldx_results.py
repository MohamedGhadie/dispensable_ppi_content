#----------------------------------------------------------------------------------------
# This script processes results from foldx calculations. Since each foldx result file 
# contains only one mutation, no second-round jobs will be produced.
#
# Requirements:
# Files must be in format produced by script produce_mutation_structure_maps.py
#----------------------------------------------------------------------------------------

import os
from pathlib import Path
from ddg_tools import (read_foldx_results,
                       write_mutation_ddg_tofile,
                       produce_foldx_and_beluga_jobs)

def main():
    
    # reference interactome name
    # options: HI-II-14, IntAct
    interactome_name = 'HI-II-14'
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # directory of processed data files specific to interactome
    interactomeDir = procDir / interactome_name
    
    # directory of foldx results
    inDir = interactomeDir / 'foldx' / 'results'
    
    # directory of foldx output jobs
    outDir = interactomeDir / 'foldx'
    
    # directory of PDB structure files
    pdbDir = Path('../pdb_files')
    
    # create output directories if not existing
    if not interactomeDir.exists():
        os.makedirs(interactomeDir)
    if not outDir.exists():
        os.makedirs(outDir)
    
    processed, unprocessed = read_foldx_results (inDir)
    
    write_mutation_ddg_tofile (processed,
                               interactomeDir / 'nondisease_mutations_foldx_ddg.txt',
                               interactomeDir / 'nondisease_mutations_foldx_ddg_2.txt',
                               type = 'binding')
    write_mutation_ddg_tofile (processed,
                               interactomeDir / 'disease_mutations_foldx_ddg.txt',
                               interactomeDir / 'disease_mutations_foldx_ddg_2.txt',
                               type = 'binding')
    
    os.remove (interactomeDir / 'nondisease_mutations_foldx_ddg.txt')
    os.remove (interactomeDir / 'disease_mutations_foldx_ddg.txt')
    os.rename (interactomeDir / 'nondisease_mutations_foldx_ddg_2.txt',
               interactomeDir / 'nondisease_mutations_foldx_ddg.txt')
    os.rename (interactomeDir / 'disease_mutations_foldx_ddg_2.txt',
               interactomeDir / 'disease_mutations_foldx_ddg.txt')
    
    produce_foldx_and_beluga_jobs (unprocessed,
                                   pdbDir,
                                   outDir,
                                   account = 'ctb-yxia',
                                   walltime = '1-00',
                                   ntasks = 1,
                                   nodes = 1,
                                   ntasks_per_node = 1,
                                   cpus_per_task = 1,
                                   mem_per_cpu = '8000M',
                                   outputfile = '/project/def-yxia/ghadie84/foldx/data/%x-%j.out',
                                   username = 'ghadie84',
                                   serverDataDir = '/project/def-yxia/ghadie84/foldx/data')

if __name__ == "__main__":
    main()
