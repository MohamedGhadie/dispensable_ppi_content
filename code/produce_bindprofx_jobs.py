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
from ddg_tools import read_unprocessed_ddg_mutations, produce_bindprofx_jobs

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
    
    # directory of bindprofx output jobs
    outDir = interactomeDir / 'bindprofx'
    
    # directory of PDB structure files
    pdbDir = Path('/Volumes/MG_Samsung/pdb_files')
    
    # input file containing mutations to submit to bindprofx
    mutationsFile = interactomeDir / 'nondisease_mutations_bindprofx_ddg.txt'
    
    # create output directories if not existing
    if not outDir.exists():
        os.makedirs(outDir)
    
    mutations = read_unprocessed_ddg_mutations (mutationsFile, 'binding')
    
    produce_bindprofx_jobs (mutations,
                            pdbDir,
                            outDir,
                            write_hpc_jobfiles = True,
                            nodes = 1,
                            ppn = 1,
                            pmem = 7700,
                            walltime = '1:00:00:00',
                            rapid = 'evf-115-aa',
                            username = 'ghadie84',
                            hpcCommands = ['source /home/ghadie84/venv/bin/activate'],
                            serverDataDir = '/sf1/project/evf-115-aa/ghadie84/bindprofx/data')

if __name__ == "__main__":
    main()
