#----------------------------------------------------------------------------------------
# This script produces jobs for mutations mapped onto structures to be submitted to 
# bindprofx for ∆∆G calculations. Mutations in the same structure are combined in one job. 
# Mutations with existing ∆∆G values in the input file are skipped.
#
# Requirements:
# Files must be in format produced by script produce_mutation_structure_maps.py
#----------------------------------------------------------------------------------------

import os
from pathlib import Path
from ddg_tools import (append_mutation_ddg_files,
                       read_unprocessed_ddg_mutations,
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
    
    # directory of bindprofx output jobs
    outDir = interactomeDir / 'bindprofx'
    
    # directory of PDB structure files
    pdbDir = Path('../pdb_files')
    
    # input data files
    nondiseaseMutFile = interactomeDir / 'nondisease_mutations_bindprofx_ddg.txt'
    diseaseMutFile = interactomeDir / 'disease_mutations_bindprofx_ddg.txt'
    
    # temporary output files
    allMutFile = interactomeDir / 'all_mutations_bindprofx_ddg.txt'
    
    # create output directories if not existing
    if not outDir.exists():
        os.makedirs(outDir)
    
    append_mutation_ddg_files (nondiseaseMutFile, diseaseMutFile, allMutFile)
    mutations = read_unprocessed_ddg_mutations (allMutFile, 'binding')
    
    produce_bindprofx_and_guillimin_jobs (mutations,
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
    os.remove(allMutFile)

if __name__ == "__main__":
    main()
