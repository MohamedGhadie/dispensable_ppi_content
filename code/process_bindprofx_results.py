#----------------------------------------------------------------------------------------
# This script processes results from bindprofx calculations, and produces second-round jobs
# from jobs that failed in the first-round of bindprofx calculations. Each produced 
# second-round job contains only one mutation.
#
# This script also processes results from bindprofx second-round calculations.
#
# Requirements:
# Files must be in format produced by script produce_mutation_structure_maps.py
#----------------------------------------------------------------------------------------

import os
from pathlib import Path
from ddg_tools import (read_bindprofx_results,
                       write_mutation_ddg_tofile,
                       produce_bindprofx_and_guillimin_jobs)

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
    
    # directory of bindprofx results
    inDir = interactomeDir / 'bindprofx' / 'results'
    
    # directory of bindprofx output jobs
    outDir = interactomeDir / 'bindprofx'
    
    # directory of PDB structure files
    pdbDir = Path('../pdb_files')
    
    # create output directories if not existing
    if not outDir.exists():
        os.makedirs(outDir)
    
    processed, unprocessed = read_bindprofx_results (inDir)
    
    write_mutation_ddg_tofile (processed,
                               interactomeDir / 'nondisease_mutations_bindprofx_ddg.txt',
                               interactomeDir / 'nondisease_mutations_bindprofx_ddg_2.txt',
                               'binding')
    write_mutation_ddg_tofile (processed,
                               interactomeDir / 'disease_mutations_bindprofx_ddg.txt',
                               interactomeDir / 'disease_mutations_bindprofx_ddg_2.txt',
                               'binding')
    
    os.remove (interactomeDir / 'nondisease_mutations_bindprofx_ddg.txt')
    os.remove (interactomeDir / 'disease_mutations_bindprofx_ddg.txt')
    os.rename (interactomeDir / 'nondisease_mutations_bindprofx_ddg_2.txt',
               interactomeDir / 'nondisease_mutations_bindprofx_ddg.txt')
    os.rename (interactomeDir / 'disease_mutations_bindprofx_ddg_2.txt',
               interactomeDir / 'disease_mutations_bindprofx_ddg.txt')
    
    produce_bindprofx_and_guillimin_jobs (unprocessed,
                                          pdbDir,
                                          outDir,
                                          nodes = 1,
                                          ppn = 1,
                                          pmem = 7700,
                                          walltime = '1:00:00:00',
                                          rapid = 'evf-115-aa',
                                          username = 'ghadie84',
                                          extraCommands = ['source /home/ghadie84/venv/bin/activate'],
                                          serverDataDir = '/sf1/project/evf-115-aa/ghadie84/bindprofx/data')

if __name__ == "__main__":
    main()
